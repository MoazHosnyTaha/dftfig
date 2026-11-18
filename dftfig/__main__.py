"""Command-line entry point.

    python -m dftfig wizard                 interactive configuration
    python -m dftfig params [--key K]       parameter reference
    python -m dftfig all -c config.yaml     run everything
    python -m dftfig <stage> -c config.yaml run one stage

Stages, in order: build, opt, sp, fukui, cubes, descriptors, figures.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import config as C
from . import cubes as cubemod
from . import descriptors as descmod
from . import geom, orca_io, properties, runner
from .theme import banner, console, ok, step, warn, bad

STAGES = ["build", "opt", "sp", "fukui", "cubes", "descriptors", "figures"]


def _mols(cfg, announce: bool = True):
    """Resolve the molecule set, with each compound's own charge and multiplicity."""
    from . import molecules as M
    specs = M.load_molecules(cfg)
    if announce:
        M.report(specs)
        path = M.write_manifest(cfg, specs)
        console.print(f"  [muted]written:[/muted] [value]{path}[/value]")
    return specs


def stage_build(cfg, specs):
    banner("build", "input structures -> starting geometry")
    return geom.build_all(cfg, specs)


def stage_opt(cfg, specs):
    banner("opt", orca_io.stage_label("opt") + f"  [{_level(cfg, 'opt')}]")
    jobs = []
    for name, spec in specs.items():
        for solvent in cfg["solvation"]["solvents"]:
            inp = orca_io.write_input(cfg, spec, solvent, "opt")
            # Parentheses, not brackets: rich would parse "[q+1 m1]" as a markup tag
            # and silently delete it from the label.
            q = f"(q{spec.charge:+d} m{spec.multiplicity})"
            jobs.append((inp, f"{name} {q} / {solvent} / optimization"))
    res = runner.run_many(cfg, jobs)
    runner.write_manifest(cfg, res, "opt")
    _summarize(res, "optimization")

    outcomes = _escape_saddles(cfg, specs)
    if outcomes:
        _record_saddle_outcomes(cfg, outcomes)
    return res


def _escape_saddles(cfg, specs):
    """Re-optimize anything that came back a saddle point, until it is a minimum.

    A frequency job that reports an imaginary mode has found a structure that is not a
    stable species. Everything downstream - the gap, the thermochemistry, the Fukui
    indices, the MEP - would then describe a geometry that does not exist, so the fix
    belongs here rather than in a warning the reader has to act on.
    """
    from . import parsers, refine

    if not cfg["theory"].get("run_freq", True):
        return []
    max_attempts = int(cfg["theory"].get("saddle_max_attempts", 0) or 0)
    cutoff = float(cfg["theory"].get("imaginary_cutoff", -10.0))
    outcomes = []
    for name, spec in specs.items():
        for solvent in cfg["solvation"]["solvents"]:
            work = Path(cfg["execution"]["work_dir"]) / name / solvent
            out = work / f"{name}_{solvent}_opt.out"
            if not out.exists():
                continue
            imag = parsers.imaginary_modes(out, cutoff)
            outcome = refine.Outcome(name=name, solvent=solvent)
            outcome.attempts.append(refine.Attempt(
                0, "as optimized", len(imag),
                min(parsers.frequencies(out), default=None)))
            if not imag:
                outcome.resolved = True
                outcomes.append(outcome)
                continue
            if max_attempts <= 0:
                bad(f"{name}/{solvent}: {len(imag)} imaginary mode(s), lowest "
                    f"{imag[0]:.2f} cm-1 - NOT a minimum. Set theory.saddle_max_attempts "
                    f"above 0 to have the pipeline fix this automatically.")
                outcomes.append(outcome)
                continue
            outcomes.append(_retry_until_minimum(cfg, spec, solvent, outcome,
                                                max_attempts, cutoff))
    return outcomes


def _retry_until_minimum(cfg, spec, solvent, outcome, max_attempts, cutoff):
    from . import parsers, refine

    name = spec.name
    work = Path(cfg["execution"]["work_dir"]) / name / solvent
    canonical_out = work / f"{name}_{solvent}_opt.out"
    canonical_xyz = work / f"{name}_{solvent}_opt.xyz"
    step_a = float(cfg["theory"].get("saddle_step_angstrom", 0.30))
    small = float(cfg["theory"].get("saddle_small_mode_cm", 50.0))

    def _score(out_path):
        """Lower is better: fewer imaginary modes first, then a smaller |frequency|."""
        im = parsers.imaginary_modes(out_path, cutoff)
        return (len(im), abs(im[0]) if im else 0.0)

    source_out, source_xyz = canonical_out, canonical_xyz
    best = (_score(canonical_out), None)          # (score, tag) - None means the original
    for attempt in range(max_attempts):
        imag = parsers.imaginary_modes(source_out, cutoff)
        if not imag:
            break
        recipe = refine.plan(imag[0], attempt, small)
        if recipe is None:
            warn(f"{name}/{solvent}: no strategy left after {attempt} attempt(s)")
            break
        label, needs_disp, sign, extra = recipe

        if not source_xyz.exists():
            warn(f"{name}/{solvent}: {source_xyz.name} missing - cannot restart")
            break
        symbols, coords = refine.read_xyz(source_xyz)
        if needs_disp:
            vecs = refine.imaginary_mode_vectors(source_out, symbols)
            if not vecs:
                warn(f"{name}/{solvent}: the output has no NORMAL MODES block, so the "
                     f"imaginary mode cannot be followed - retrying with a finer grid "
                     f"instead")
                label, needs_disp, extra = ("finer grid, tighter optimization", False,
                                            ["DefGrid3", "TightOpt", "TightSCF"])
            else:
                coords = refine.displaced_geometry(symbols, coords, vecs[0][1],
                                                   step_a, sign)

        tag = f"r{attempt + 1}"
        start = work / f"{name}_{solvent}_opt_{tag}_start.xyz"
        refine.write_xyz(start, symbols, coords,
                         f"dftfig saddle-point restart {tag}: {label}")
        step(f"{name}/{solvent}: {len(imag)} imaginary mode(s), lowest {imag[0]:.2f} cm-1 "
             f"-> restart {attempt + 1}: {label}")

        inp = orca_io.write_input(cfg, spec, solvent, "opt",
                                  geom_xyz=start.name, extra_keywords=extra,
                                  basename=f"{name}_{solvent}_opt_{tag}")
        runner.run_many(cfg, [(inp, f"{name} / {solvent} / saddle restart {attempt + 1}")])

        new_out = work / f"{name}_{solvent}_opt_{tag}.out"
        new_xyz = work / f"{name}_{solvent}_opt_{tag}.xyz"
        if not new_out.exists():
            warn(f"{name}/{solvent}: restart {attempt + 1} produced no output")
            break
        new_imag = parsers.imaginary_modes(new_out, cutoff)
        outcome.attempts.append(refine.Attempt(
            attempt + 1, label, len(new_imag),
            min(parsers.frequencies(new_out), default=None)))

        # Keep the best structure seen, not the most recent one. A restart can land on a
        # worse stationary point than it started from - toluene's methyl rotor went
        # -31 -> -28 -> -53 cm-1 - and promoting the last attempt would ship that.
        if _score(new_out) < best[0]:
            best = (_score(new_out), tag)
        source_out, source_xyz = new_out, new_xyz
        if not new_imag:
            outcome.resolved = True
            _promote(work, name, solvent, tag)
            ok(f"{name}/{solvent}: {outcome.summary()}")
            return outcome

    # Nothing reached a clean minimum. Promote whichever attempt came closest, and be
    # specific about what is left.
    if best[1] is not None:
        _promote(work, name, solvent, best[1])
    resid = best[0][1]
    outcome.best_residual_cm = resid
    if resid and resid < small:
        outcome.flat_mode = True
        warn(f"{name}/{solvent}: a residual imaginary mode of -{resid:.1f} cm-1 survived "
             f"{len(outcome.attempts) - 1} restart(s). A mode this small is a nearly free "
             f"internal rotation - a methyl or hydroxyl turning - not a reaction saddle. "
             f"The geometry is a minimum for every practical purpose; the harmonic "
             f"frequencies are simply the wrong model for that one coordinate, and the "
             f"thermochemistry needs a hindered-rotor treatment rather than more "
             f"optimization.")
    else:
        bad(f"{name}/{solvent}: {outcome.summary()} - the reported numbers describe a "
            f"saddle point, not a stable species")
    return outcome


def _promote(work: Path, name: str, solvent: str, tag: str) -> None:
    """Make the successful restart the geometry every later stage reads.

    The original saddle-point files are kept under *_saddle.* so the history stays on
    disk and nothing is silently overwritten.
    """
    import shutil
    for ext in (".out", ".xyz", ".gbw", ".hess", ".inp"):
        canon = work / f"{name}_{solvent}_opt{ext}"
        new = work / f"{name}_{solvent}_opt_{tag}{ext}"
        if not new.exists():
            continue
        if canon.exists():
            keep = work / f"{name}_{solvent}_opt_saddle{ext}"
            if not keep.exists():
                shutil.copy2(canon, keep)
        shutil.copy2(new, canon)


def _record_saddle_outcomes(cfg, outcomes) -> None:
    """Write the restart history next to the other tables."""
    import csv
    work = Path(cfg["execution"]["work_dir"])
    work.mkdir(parents=True, exist_ok=True)
    path = work / "saddle_history.csv"
    with path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["molecule", "medium", "attempt", "strategy", "n_imaginary",
                    "lowest_freq_cm-1", "resolved", "flat_mode", "summary"])
        for o in outcomes:
            for a in o.attempts:
                w.writerow([o.name, o.solvent, a.index, a.strategy, a.n_imaginary,
                            "" if a.lowest_freq is None else f"{a.lowest_freq:.2f}",
                            o.resolved, o.flat_mode, o.summary()])
    fixed = [o for o in outcomes if o.resolved and o.n_restarts]
    if fixed:
        ok(f"saddle points escaped for {', '.join(o.name for o in fixed)} "
           f"-> {path.name}")


def stage_sp(cfg, specs):
    if not cfg["property"]["enabled"]:
        warn("property.enabled is false - skipping the property single point")
        return []
    banner("sp", orca_io.stage_label("sp") + f"  [{_level(cfg, 'sp')}]")
    jobs = []
    for name, spec in specs.items():
        for solvent in cfg["solvation"]["solvents"]:
            d = Path(cfg["execution"]["work_dir"]) / name / solvent
            optxyz = d / f"{name}_{solvent}_opt.xyz"
            if not optxyz.exists():
                warn(f"{name}/{solvent}: optimized geometry {optxyz.name} missing - "
                     "run the opt stage first")
                continue
            inp = orca_io.write_input(cfg, spec, solvent, "sp", geom_xyz=optxyz.name)
            jobs.append((inp, f"{name} / {solvent} / property single point"))
    res = runner.run_many(cfg, jobs)
    runner.write_manifest(cfg, res, "sp")
    _summarize(res, "property single point")
    return res


def stage_fukui(cfg, specs):
    if not cfg["fukui"]["enabled"]:
        warn("fukui.enabled is false - skipping the ion single points")
        return []
    banner("fukui", "N-1 and N+1 single points at the reference geometry")
    jobs = []
    for name, spec in specs.items():
        for solvent in cfg["solvation"]["solvents"]:
            d = Path(cfg["execution"]["work_dir"]) / name / solvent
            optxyz = d / f"{name}_{solvent}_opt.xyz"
            if not optxyz.exists():
                continue
            for stg, delta in (("cat", -1), ("ani", +1)):
                q, mult = spec.ion(delta)
                inp = orca_io.write_input(cfg, spec, solvent, stg, geom_xyz=optxyz.name)
                what = "cation (one electron removed)" if delta < 0 else "anion (one electron added)"
                jobs.append((inp, f"{name} / {solvent} / {what}, q={q:+d} m={mult}"))
    res = runner.run_many(cfg, jobs)
    runner.write_manifest(cfg, res, "fukui")
    _summarize(res, "ion single points")
    return res


def stage_cubes(cfg):
    if not cfg["cubes"]["enabled"]:
        warn("cubes.enabled is false - skipping volumetric data")
        return {}
    banner("cubes", "electrostatic potential on the electron-density grid")
    return cubemod.generate_all(cfg)


def stage_descriptors(cfg):
    banner("descriptors", "global reactivity indices and condensed Fukui functions")
    df = descmod.collect(cfg)
    descmod.fukui(cfg)
    if not df.empty:
        _print_descriptor_table(cfg, df)
    return df


def stage_figures(cfg, specs, df):
    banner("figures", ", ".join(cfg["figures"]["which"]))
    from . import mlp as mlpmod
    mols = {n: s.mol for n, s in specs.items()}
    rd = {n: mlpmod.rdkit_properties(m) for n, m in mols.items()}
    prop_df = properties.load_properties(cfg, rd, mols=mols)
    out = Path(cfg["execution"]["work_dir"]) / "figures"
    out.mkdir(parents=True, exist_ok=True)
    prop_df.to_csv(out / "properties_resolved.csv", index=False, float_format="%.4f")

    made = []
    which = cfg["figures"]["which"]
    if "fig2" in which:
        from .figures import fig2
        made += fig2.build(cfg, mols, df)
    if "fig12" in which:
        from .figures import fig12
        made += fig12.build(cfg, mols, prop_df,
                            charges={n: s.charge for n, s in specs.items()})
    if "fig13" in which:
        from .figures import fig13
        made += fig13.build(cfg, mols, prop_df)
    if "fig14" in which:
        from .figures import fig14
        made += fig14.build(cfg, mols, prop_df)
    return made


def _level(cfg, stage):
    return orca_io.keyword_line(cfg, stage, cfg["figures"]["reference_solvent"])[2:]


def _summarize(res, what):
    if not res:
        return
    good = [r for r in res if r.get("ok")]
    failed = [r for r in res if not r.get("ok")]
    imag = [r for r in good if r.get("imaginary_modes")]
    total = sum(r.get("wall_s", 0) or 0 for r in res)
    console.print(f"  [muted]{what}: {len(good)}/{len(res)} succeeded, "
                  f"{total / 60:.1f} min total[/muted]")
    for r in failed:
        bad(f"{r.get('label')}: {r.get('error')}")
    for r in imag:
        bad(f"{r.get('label')}: NOT a minimum - imaginary modes at "
            f"{', '.join(f'{f:.1f}' for f in r['imaginary_modes'])} cm^-1")


def _print_descriptor_table(cfg, df):
    from rich.table import Table
    ref = cfg["figures"]["reference_solvent"]
    sub = df[df["medium"] == ref]
    if sub.empty:
        return
    suffix = "_koopmans" if "IE_eV_koopmans" in df.columns else "_deltascf"
    t = Table(title=f"global reactivity descriptors in {ref}  ({suffix[1:]})",
              header_style="brand", box=None, padding=(0, 2), title_style="brand")
    cols = ["molecule", "e_homo_eV", "e_lumo_eV", "gap_eV",
            f"IE_eV{suffix}", f"EA_eV{suffix}", f"chi_eV{suffix}",
            f"eta_eV{suffix}", f"omega_eV{suffix}"]
    heads = ["molecule", "E(HOMO)", "E(LUMO)", "gap", "IE", "EA", "chi", "eta", "omega"]
    for h in heads:
        t.add_column(h, justify="right" if h != "molecule" else "left")
    for _, r in sub.iterrows():
        cells = [str(r["molecule"])]
        for c in cols[1:]:
            v = r.get(c)
            cells.append("-" if v is None or (isinstance(v, float) and v != v) else f"{v:.3f}")
        t.add_row(*cells)
    console.print(t)
    console.print("  [muted]all values in eV; see descriptors.csv for the full table "
                  "including both IE/EA derivations[/muted]")


def stage_status(cfg, specs=None):
    from rich.table import Table
    from . import parsers
    if specs is None:
        specs = _mols(cfg, announce=False)
    work = Path(cfg["execution"]["work_dir"])
    t = Table(header_style="brand", box=None, padding=(0, 2))
    for h in ["molecule", "medium", "opt", "freq", "sp", "cat", "ani", "cubes", "esp"]:
        t.add_column(h)
    for name, spec in specs.items():
        for solvent in cfg["solvation"]["solvents"]:
            d = work / name / solvent
            cells = [f"{name} [muted]q{spec.charge:+d}[/muted]", solvent]
            for stg in ("opt", "sp", "cat", "ani"):
                f = d / f"{name}_{solvent}_{stg}.out"
                if stg == "opt":
                    if f.exists() and parsers.terminated_normally(f):
                        cells.append("[ok]done[/ok]")
                        imag = parsers.imaginary_modes(
                            f, float(cfg["theory"]["imaginary_cutoff"]))
                        freqs = parsers.frequencies(f)
                        cells.append("[bad]imaginary[/bad]" if imag else
                                     ("[ok]minimum[/ok]" if freqs else "[muted]-[/muted]"))
                    else:
                        cells += ["[bad]no[/bad]", "[muted]-[/muted]"]
                else:
                    cells.append("[ok]done[/ok]" if f.exists() and parsers.terminated_normally(f)
                                 else "[bad]no[/bad]")
            cells.append("[ok]yes[/ok]" if (d / "eldens.cube").exists() else "[bad]no[/bad]")
            cells.append("[ok]yes[/ok]" if (d / "esp.cube").exists() else "[bad]no[/bad]")
            t.add_row(*cells)
    console.print(t)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="dftfig", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stage", choices=STAGES + ["all", "wizard", "params", "status", "validate"])
    ap.add_argument("-c", "--config", default="config.yaml")
    ap.add_argument("--key", default=None, help="with 'params': show one parameter")
    ap.add_argument("--no-advanced", action="store_true",
                    help="with 'params': hide advanced parameters")
    args = ap.parse_args(argv)

    if args.stage == "wizard":
        from .wizard import run_wizard
        run_wizard(args.config)
        return 0
    if args.stage == "params":
        from .wizard import show_params
        show_params(args.key, advanced=not args.no_advanced)
        return 0

    if not Path(args.config).exists():
        bad(f"no config file at {args.config}. Create one with: python -m dftfig wizard")
        return 2
    try:
        cfg = C.load(args.config)
    except C.ConfigError as e:
        bad(str(e))
        return 2

    if args.stage == "validate":
        ok("configuration is valid")
        est = C.estimate(cfg)
        for k, v in est.items():
            console.print(f"  [muted]{k}:[/muted] [value]{v}[/value]")
        mc = C.machine()
        console.print(f"  [muted]this machine:[/muted] [value]{mc['cores']} cores"
                      + (f", {mc['ram_gb']:.1f} GB RAM" if mc["ram_gb"] else "") + "[/value]")
        for a in C.advisories(cfg):
            warn(a)
        return 0
    if args.stage == "status":
        try:
            stage_status(cfg)
        except (ValueError, FileNotFoundError) as e:
            bad(str(e))
            return 2
        return 0

    try:
        specs = _mols(cfg)
    except (ValueError, FileNotFoundError) as e:
        bad(str(e))
        return 2
    df = None
    todo = STAGES if args.stage == "all" else [args.stage]
    for s in todo:
        if s == "build":
            stage_build(cfg, specs)
        elif s == "opt":
            stage_opt(cfg, specs)
        elif s == "sp":
            stage_sp(cfg, specs)
        elif s == "fukui":
            stage_fukui(cfg, specs)
        elif s == "cubes":
            stage_cubes(cfg)
        elif s == "descriptors":
            df = stage_descriptors(cfg)
        elif s == "figures":
            if df is None:
                df = descmod.collect(cfg)
            made = stage_figures(cfg, specs, df)
            for p in made:
                ok(str(p))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        from . import runner as _r
        _r.stop_all()
        sys.exit(130)
