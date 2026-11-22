"""Run ORCA jobs locally, with resume, timeouts and MPI fallback."""
from __future__ import annotations

import json
import os
import signal
import subprocess
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from rich.progress import (BarColumn, MofNCompleteColumn, Progress, SpinnerColumn,
                           TextColumn, TimeElapsedColumn)

from . import parsers
from .theme import PALETTE, console, note, warn

#: Live ORCA child processes, so an interrupt can stop them instead of orphaning them.
_LIVE: dict[int, subprocess.Popen] = {}
_LIVE_LOCK = threading.Lock()
_CANCELLED = threading.Event()


def _register(proc: subprocess.Popen) -> None:
    with _LIVE_LOCK:
        _LIVE[proc.pid] = proc


def _unregister(proc: subprocess.Popen) -> None:
    with _LIVE_LOCK:
        _LIVE.pop(proc.pid, None)


def stop_all(grace: float = 5.0) -> int:
    """Terminate every running ORCA process. Returns how many were stopped.

    ORCA spawns helper processes, so each job is started in its own process group and the
    whole group is signalled - otherwise the children outlive the parent and keep burning
    cores after the pipeline has exited.
    """
    _CANCELLED.set()
    with _LIVE_LOCK:
        procs = list(_LIVE.values())
    for p in procs:
        try:
            _kill_tree(p, signal.SIGTERM)
        except (ProcessLookupError, PermissionError, OSError):
            try:
                p.terminate()
            except Exception:
                pass
    deadline = time.time() + grace
    for p in procs:
        try:
            p.wait(timeout=max(0.1, deadline - time.time()))
        except Exception:
            try:
                _kill_tree(p, signal.SIGKILL)
            except Exception:
                pass
    return len(procs)


class JobResult(dict):
    @property
    def succeeded(self) -> bool:
        return bool(self.get("ok"))


def _kill_tree(proc, sig) -> None:
    """Stop a job and everything it spawned.

    ORCA launches helper processes (mpirun, orca_gtoint, ...), so killing only the
    parent leaves them running and holding the scratch files. On POSIX each job is
    started in its own process group and the group is signalled. Windows has no
    process groups in this sense; terminate()/kill() plus taskkill /T is the
    equivalent, so the behaviour is the same even though the mechanism is not.
    """
    import os
    import signal as _signal
    import subprocess as _sp
    if proc is None or proc.poll() is not None:
        return
    if os.name == "posix":
        try:
            os.killpg(os.getpgid(proc.pid), sig)
            return
        except (ProcessLookupError, PermissionError, OSError):
            pass
    else:
        try:
            _sp.run(["taskkill", "/F", "/T", "/PID", str(proc.pid)],
                    capture_output=True, timeout=15)
            return
        except Exception:
            pass
    try:
        proc.kill() if sig == getattr(_signal, "SIGKILL", 9) else proc.terminate()
    except Exception:
        pass


def _orca_tools(orca_path: str) -> dict:
    p = Path(orca_path)
    return {"orca": str(p), "orca_plot": str(p.parent / "orca_plot"),
            "orca_vpot": str(p.parent / "orca_vpot"),
            "orca_vpot_mpi": str(p.parent / "orca_vpot_mpi")}


def _strip_pal(inp: Path) -> None:
    """Remove the %pal block so ORCA runs serially (MPI fallback)."""
    txt = inp.read_text()
    out = "\n".join(l for l in txt.splitlines() if not l.strip().lower().startswith("%pal"))
    inp.write_text(out + "\n")


def _set_pal(inp: Path, nprocs: int) -> None:
    """Rewrite the %pal line to a different process count."""
    if nprocs <= 1:
        _strip_pal(inp)
        return
    lines = []
    for l in inp.read_text().splitlines():
        if l.strip().lower().startswith("%pal"):
            lines.append(f"%pal nprocs {nprocs} end")
        else:
            lines.append(l)
    inp.write_text("\n".join(lines) + "\n")


def _keyword_line_of(path: Path) -> str | None:
    """The ORCA '!' line, read from an input file or from its echo in an output."""
    try:
        text = path.read_text(errors="replace")
    except OSError:
        return None
    for raw in text.splitlines():
        line = raw.strip()
        # ORCA echoes the input into the output as "|  1> ! B3LYP ...".
        if line.startswith("|") and ">" in line:
            line = line.split(">", 1)[1].strip()
        if line.startswith("!"):
            return " ".join(line.split())
        if line.startswith("*") or line.lower().startswith("%coords"):
            break            # past the keyword block
    return None


def _keywords_changed(inp: Path, out: Path) -> tuple[str, str] | None:
    """(old, new) when a finished output was produced at a different level of theory.

    `resume` reuses any job that terminated normally, which is right for an interrupted
    run and wrong after the functional or basis is changed: the numbers would stay at
    the old level while the figure footnote claims the new one. Only the keyword line is
    compared, so a change in process count still resumes.
    """
    want, have = _keyword_line_of(inp), _keyword_line_of(out)
    if want is None or have is None or want == have:
        return None
    return have, want


def run_job(cfg: dict, inp: Path, label: str = "") -> JobResult:
    """Execute one ORCA input. Returns a JobResult; never raises on calculation failure."""
    inp = Path(inp)
    out = inp.with_suffix(".out")
    tools = _orca_tools(cfg["execution"]["orca_path"])

    if cfg["execution"]["resume"] and out.exists() and parsers.terminated_normally(out):
        stale = _keywords_changed(inp, out)
        if stale is None:
            return JobResult(ok=True, input=str(inp), output=str(out), reused=True,
                             wall_s=parsers.wall_time_s(out) or 0.0, label=label)
        note(f"{label or inp.stem}: recomputing, the level of theory changed\n"
             f"      was: {stale[0]}\n      now: {stale[1]}")

    if cfg["execution"]["dry_run"]:
        return JobResult(ok=True, input=str(inp), output=None, dry_run=True, label=label)

    if not Path(tools["orca"]).exists():
        return JobResult(ok=False, input=str(inp), error=f"ORCA not found at {tools['orca']}",
                         label=label)

    env = dict(os.environ)
    scratch = Path(cfg["execution"]["scratch_dir"]).resolve()
    scratch.mkdir(parents=True, exist_ok=True)
    env.setdefault("ORCA_SCRDIR", str(scratch))
    env["PATH"] = str(Path(tools["orca"]).parent) + os.pathsep + env.get("PATH", "")

    timeout_s = int(cfg["execution"]["timeout_min"]) * 60
    t0 = time.time()
    for attempt in ("parallel", "serial"):
        if _CANCELLED.is_set():
            return JobResult(ok=False, input=str(inp), output=str(out), label=label,
                             error="cancelled before it started", wall_s=0.0, cancelled=True)
        proc = None
        try:
            with open(out, "w") as fh:
                # start_new_session: the job gets its own process group, so an interrupt
                # can stop ORCA and every helper it spawned.
                proc = subprocess.Popen([tools["orca"], inp.name], cwd=str(inp.parent),
                                        stdout=fh, stderr=subprocess.STDOUT, env=env,
                                        start_new_session=(os.name == 'posix'))
                _register(proc)
                proc.wait(timeout=timeout_s)
        except subprocess.TimeoutExpired:
            try:
                _kill_tree(proc, signal.SIGKILL)
            except Exception:
                pass
            return JobResult(ok=False, input=str(inp), output=str(out), label=label,
                             error=f"timed out after {cfg['execution']['timeout_min']} min",
                             wall_s=time.time() - t0)
        finally:
            if proc is not None:
                _unregister(proc)

        if _CANCELLED.is_set():
            return JobResult(ok=False, input=str(inp), output=str(out), label=label,
                             error="cancelled", wall_s=time.time() - t0, cancelled=True)

        if parsers.terminated_normally(out):
            break
        # One automatic retry serially if MPI startup is what failed. The output file may
        # not exist at all if ORCA never started, so check before reading it.
        produced = out.read_text(errors="replace") if out.exists() else ""
        if attempt == "parallel" and "topology discovery failed" in produced:
            warn(f"{label or inp.name}: Open MPI could not start (no CPU topology visible). "
                 "Retrying this job serially.")
            _strip_pal(inp)
            continue
        # Asking for more MPI slots than there are physical cores is a startup failure,
        # not a calculation failure: nothing has been computed, so retrying at the real
        # core count costs only the seconds already spent printing the banner.
        if attempt == "parallel" and "not enough slots available" in produced:
            from .config import physical_cores
            phys = physical_cores()
            asked = int(cfg["execution"]["nprocs"])
            warn(f"{label or inp.name}: Open MPI refused {asked} slots - a slot is a "
                 f"physical core and this machine has {phys}. Retrying this job on "
                 f"{phys} core(s). Set execution.nprocs to {phys} to avoid the false "
                 f"start on every job.")
            _set_pal(inp, phys)
            continue
        break

    wall = time.time() - t0
    if not parsers.terminated_normally(out):
        return JobResult(ok=False, input=str(inp), output=str(out), label=label,
                         error=parsers.error_message(out), wall_s=wall)

    res = JobResult(ok=True, input=str(inp), output=str(out), wall_s=wall, label=label,
                    energy=parsers.final_energy(out))
    cutoff = float(cfg["theory"]["imaginary_cutoff"])
    imag = parsers.imaginary_modes(out, cutoff)
    if imag:
        res["imaginary_modes"] = imag
    return res


def run_many(cfg: dict, jobs: list[tuple[Path, str]]) -> list[JobResult]:
    """Run jobs with execution.max_parallel_jobs concurrency and a live progress display.

    The display matters: an optimization+frequency job on a drug-sized molecule runs for
    an hour or more, and a console that prints nothing until the first one finishes is
    indistinguishable from a hang.
    """
    nmax = max(1, int(cfg["execution"]["max_parallel_jobs"]))
    results: list[JobResult] = []
    if not jobs:
        return results

    _CANCELLED.clear()
    console.print(f"  [muted]{len(jobs)} job(s), {nmax} at a time, "
                  f"{cfg['execution']['nprocs']} core(s) each[/muted]")
    console.print("  [muted]ORCA prints nothing until a job ends; the timer below is live. "
                  "Ctrl+C stops cleanly and `dft run` resumes where it left off.[/muted]")

    progress = Progress(
        SpinnerColumn(style=PALETTE["stage"]),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(bar_width=18, complete_style=PALETTE["ok"]),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
        console=console, transient=False, refresh_per_second=4,
    )

    def _execute(inp: Path, label: str) -> JobResult:
        tid = progress.add_task(f"[info]{label}[/info]", total=None)
        try:
            r = run_job(cfg, inp, label)
        finally:
            progress.remove_task(tid)
        return r

    with progress:
        overall = progress.add_task("[stage]all jobs[/stage]", total=len(jobs))
        try:
            if nmax == 1:
                for inp, label in jobs:
                    r = _execute(inp, label)
                    _report(r, progress.console)
                    results.append(r)
                    progress.advance(overall)
            else:
                with ThreadPoolExecutor(max_workers=nmax) as ex:
                    futs = {ex.submit(_execute, inp, label): label for inp, label in jobs}
                    try:
                        for fut in _as_completed(futs):
                            r = fut.result()
                            _report(r, progress.console)
                            results.append(r)
                            progress.advance(overall)
                    except KeyboardInterrupt:
                        _interrupt(ex, futs)
                        raise
        except KeyboardInterrupt:
            progress.stop()
            _cancel_notice(results, len(jobs))
            raise

    return results


def _as_completed(futs):
    from concurrent.futures import as_completed
    return as_completed(futs)


def _interrupt(ex, futs) -> None:
    """Stop ORCA, then drop anything still queued."""
    n = stop_all()
    for f in futs:
        f.cancel()
    try:
        ex.shutdown(wait=False, cancel_futures=True)
    except TypeError:                      # Python < 3.9
        ex.shutdown(wait=False)
    if n:
        warn(f"terminated {n} running ORCA process(es)")


def _cancel_notice(results, total: int) -> None:
    done = sum(1 for r in results if r.get("ok"))
    console.print()
    warn(f"interrupted - {done} of {total} job(s) finished and are saved")
    console.print("  [muted]nothing is lost: completed jobs are skipped on the next run "
                  "(execution.resume), and the interrupted one restarts from its input.[/muted]")
    console.print("  [muted]resume with:[/muted] [value]dft run[/value]   "
                  "[muted]check progress with:[/muted] [value]dft status[/value]\n")


def _report(r: JobResult, out=None) -> None:
    out = out or console
    if r.get("cancelled"):
        return
    if r.get("dry_run"):
        out.print(f"  [ok]OK[/ok]  input written (dry run): {Path(r['input']).name}")
        return
    if r.get("reused"):
        out.print(f"  [ok]OK[/ok]  reusing completed output "
                  f"({r['wall_s']:.0f} s originally): {r.get('label','')}")
        return
    if r.succeeded:
        msg = f"done in {r.get('wall_s', 0):.0f} s"
        if r.get("energy") is not None:
            msg += f", E = {r['energy']:.8f} Eh"
        out.print(f"  [ok]OK[/ok]  {r.get('label','')} - {msg}")
        if r.get("imaginary_modes"):
            out.print(f"  [bad]X[/bad]  {r.get('label','')}: IMAGINARY FREQUENCIES "
                      f"{', '.join(f'{f:.1f}' for f in r['imaginary_modes'])} cm^-1 - "
                      "a saddle point, not a minimum")
    else:
        out.print(f"  [bad]X[/bad]  {r.get('label','')} failed: {r.get('error', 'unknown')}")


def write_manifest(cfg: dict, results: list[JobResult], stage: str) -> Path:
    path = Path(cfg["execution"]["work_dir"]) / f"manifest_{stage}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps([dict(r) for r in results], indent=2, default=str))
    return path
