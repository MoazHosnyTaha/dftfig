"""Provenance: tying every number and every figure back to the calculation behind it.

A figure in a manuscript is a claim, and a reader has to be able to check it. Three
records carry that, each answering a different question:

  descriptors.csv columns   what produced THIS ROW - the dftfig and ORCA versions, the
                            literal keyword lines of the two jobs behind it, and a hash
                            of the settings that could have changed the numbers
  runs/manifest.json        what the run as a whole produced - every input, output, cube
                            and figure file with its size and SHA-256
  the Figure 2 stamp        the two identifiers (version, config hash) that let a printed
                            figure be matched back to the run directory it came from

`dft verify` re-checks the manifest against what is on disk and exits non-zero if
anything is missing, altered or unaccounted for.
"""
from __future__ import annotations

import hashlib
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import yaml

from . import __version__
from . import parsers
from .theme import bad, ok, warn

MANIFEST_NAME = "manifest.json"

#: Sections that cannot change a number, excluded from `config_sha` wholesale.
#:
#:   render   how a cube is drawn - camera, isovalue, colours, image size. Changes the
#:            picture, never the underlying potential or orbital. Also where pymol_path
#:            lives.
#:   figures  layout, DPI, output formats, panel labels, and which figures to draw.
#:
#: One honest caveat, stated here rather than buried: `figures.reference_solvent`
#: selects WHICH medium the Fukui table and the figures are built from. It does not
#: change any number in descriptors.csv - that table has a row per molecule x medium
#: regardless - but two runs with the same config_sha can have produced Fukui CSVs for
#: different media. The medium is in the file name and in the `medium` column, so it is
#: never ambiguous in practice.
HASH_EXCLUDED_SECTIONS = ("render", "figures")

#: Individual keys excluded for the same reason: they decide where the work happens and
#: how fast, not what comes out.
#:
#:   work_dir, scratch_dir   where files are written
#:   orca_path               which copy of ORCA on this machine; the VERSION that ran is
#:                           recorded separately, parsed from the output itself
#:   nprocs, maxcore_mb,     parallelism and memory. This is the case the hash exists
#:   max_parallel_jobs,      for: the same chemistry run on a laptop and on a cluster
#:   cubes.esp_nprocs        must hash identically, or nobody can use the hash to decide
#:                           whether two results are comparable.
#:   resume                  whether finished jobs are reused or redone
#:
#: Everything else is in, including theory, basis sets, solvation, population scheme,
#: grid, convergence thresholds, conformer search settings and the molecule set.
HASH_EXCLUDED_KEYS = (
    "execution.work_dir",
    "execution.scratch_dir",
    "execution.orca_path",
    "execution.nprocs",
    "execution.maxcore_mb",
    "execution.max_parallel_jobs",
    "execution.resume",
    "cubes.esp_nprocs",
)

#: File kinds the manifest tracks - everything a number or a figure rests on.
#:
#: .gbw and .hess are deliberately absent. They are binary ORCA restart files, often
#: hundreds of megabytes, and ORCA does not write them byte-identically for an identical
#: calculation, so including them would make `dft verify` report differences that mean
#: nothing. The .out file is the record of what the calculation found.
TRACKED_SUFFIXES = (".inp", ".out", ".xyz", ".cube", ".csv", ".pml", ".pdb",
                    ".png", ".pdf", ".svg")


# --------------------------------------------------------------------- configuration
def hashable_config(cfg: dict) -> dict:
    """The part of the configuration that can change a number."""
    out: dict = {}
    for section, params in cfg.items():
        if section in HASH_EXCLUDED_SECTIONS:
            continue
        if not isinstance(params, dict):
            out[section] = params
            continue
        out[section] = {k: v for k, v in params.items()
                        if f"{section}.{k}" not in HASH_EXCLUDED_KEYS}
    return out


def config_sha(cfg: dict, length: int = 12) -> str:
    """Short hash of the number-affecting configuration.

    Serialized through YAML with sorted keys so the hash depends on the values and not
    on the order a wizard or a hand-edited file happened to write them in.
    """
    text = yaml.safe_dump(hashable_config(cfg), sort_keys=True, default_flow_style=False)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()[:length]


# ------------------------------------------------------------------------ ORCA version
def orca_version_of(path) -> str:
    """ORCA's version from one output file, read from its banner.

    Only the head of the file is read: the banner is in the first page and an ORCA
    output can be tens of megabytes.
    """
    p = Path(path)
    try:
        with p.open(errors="replace") as fh:
            head = fh.read(16384)
    except OSError:
        return ""
    return parsers.program_version(head)


def orca_version(cfg: dict) -> str:
    """The ORCA version(s) that produced the outputs in this run directory.

    Joined with '+' when a resumed run contains more than one - that is a fact about the
    data worth seeing, not something to silently pick a winner from.
    """
    work = Path(cfg["execution"]["work_dir"])
    seen: list[str] = []
    for out in sorted(work.rglob("*.out")):
        v = orca_version_of(out)
        if v and v not in seen:
            seen.append(v)
    return "+".join(seen)


def keyword_lines(cfg: dict, name: str, solvent: str) -> dict[str, str]:
    """The literal '!' lines of the opt and property jobs for one molecule and medium.

    Preference is the output file, because that is the input ORCA actually consumed; the
    .inp is the fallback for a dry run, where no output exists yet.
    """
    d = Path(cfg["execution"]["work_dir"]) / name / solvent
    out: dict[str, str] = {}
    for stage, column in (("opt", "orca_keywords_opt"), ("sp", "orca_keywords_sp")):
        line = ""
        o = d / f"{name}_{solvent}_{stage}.out"
        if o.exists():
            line = parsers.input_keyword_line(o)
        if not line:
            i = d / f"{name}_{solvent}_{stage}.inp"
            if i.exists():
                for raw in i.read_text(errors="replace").splitlines():
                    if raw.strip().startswith("!"):
                        line = raw.strip()
                        break
        out[column] = line
    return out


def stamp(cfg: dict) -> dict:
    """The run-level identifiers that go on every row of descriptors.csv and on Figure 2."""
    return {"dftfig_version": __version__, "config_sha": config_sha(cfg)}


# ---------------------------------------------------------------------- run manifest
def _tracked_files(cfg: dict) -> list[Path]:
    work = Path(cfg["execution"]["work_dir"])
    scratch = Path(cfg["execution"]["scratch_dir"]).resolve()
    files = []
    for p in sorted(work.rglob("*")):
        if not p.is_file() or p.is_symlink():
            continue
        # A scratch directory placed inside the run directory holds ORCA's temporary
        # files, which are deleted on completion - listing them guarantees a failing
        # verify.
        if scratch in p.resolve().parents:
            continue
        if p.name == MANIFEST_NAME or p.suffix.lower() not in TRACKED_SUFFIXES:
            continue
        files.append(p)
    return files


def sha256_file(path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def build(cfg: dict) -> dict:
    """Assemble the manifest for the run in cfg['execution']['work_dir']."""
    work = Path(cfg["execution"]["work_dir"])
    files = []
    for p in _tracked_files(cfg):
        files.append({
            "path": p.relative_to(work).as_posix(),
            "bytes": p.stat().st_size,
            "sha256": sha256_file(p),
        })
    return {
        "dftfig_version": __version__,
        "orca_version": orca_version(cfg),
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds").replace(
            "+00:00", "Z"),
        "platform": {
            "system": platform.platform(),
            "machine": platform.machine(),
            "python": platform.python_version(),
        },
        "config_sha": config_sha(cfg),
        "config_sha_excludes": {
            "sections": list(HASH_EXCLUDED_SECTIONS),
            "keys": list(HASH_EXCLUDED_KEYS),
        },
        "config": cfg,
        "tracked_suffixes": list(TRACKED_SUFFIXES),
        "files": files,
    }


def write(cfg: dict) -> Path:
    """Write runs/manifest.json and return its path."""
    work = Path(cfg["execution"]["work_dir"])
    work.mkdir(parents=True, exist_ok=True)
    path = work / MANIFEST_NAME
    man = build(cfg)
    path.write_text(json.dumps(man, indent=2, default=str))
    ok(f"provenance manifest: {len(man['files'])} file(s) -> {path}")
    return path


def load(cfg: dict) -> dict | None:
    path = Path(cfg["execution"]["work_dir"]) / MANIFEST_NAME
    if not path.exists():
        return None
    return json.loads(path.read_text())


def verify(cfg: dict) -> dict:
    """Compare the manifest with the files on disk.

    Returns {'missing', 'changed', 'unlisted', 'ok'} - lists of relative paths plus the
    count that matched. `unlisted` is restricted to the file kinds the manifest tracks,
    so a user's own notes.txt in the run directory is not reported as a discrepancy.
    """
    man = load(cfg)
    if man is None:
        raise FileNotFoundError(
            f"no {MANIFEST_NAME} in {cfg['execution']['work_dir']} - run a stage first, "
            f"or the run predates provenance tracking")

    work = Path(cfg["execution"]["work_dir"])
    listed = {e["path"]: e for e in man.get("files", [])}
    # Verify against the suffixes the manifest itself declares: a manifest written by an
    # older version that tracked fewer kinds must not have every new kind reported as
    # unlisted.
    suffixes = tuple(man.get("tracked_suffixes") or TRACKED_SUFFIXES)
    on_disk = {p.relative_to(work).as_posix() for p in _tracked_files(cfg)
               if p.suffix.lower() in suffixes}

    missing, changed, matched = [], [], 0
    for rel, entry in sorted(listed.items()):
        p = work / rel
        if not p.exists():
            missing.append(rel)
            continue
        if p.stat().st_size != entry.get("bytes") or sha256_file(p) != entry.get("sha256"):
            changed.append(rel)
            continue
        matched += 1
    unlisted = sorted(on_disk - set(listed))
    return {"missing": missing, "changed": changed, "unlisted": unlisted,
            "ok": matched, "manifest": man}


def verify_report(cfg: dict) -> int:
    """Print the result of `verify` and return a shell exit code (0 = everything agrees)."""
    from .theme import console

    try:
        r = verify(cfg)
    except FileNotFoundError as e:
        bad(str(e))
        return 2

    man = r["manifest"]
    console.print(f"  [muted]manifest written[/muted] [value]{man.get('generated_utc')}[/value] "
                  f"[muted]by dftfig[/muted] [value]{man.get('dftfig_version')}[/value] "
                  f"[muted]with ORCA[/muted] [value]{man.get('orca_version') or 'unknown'}[/value]")
    console.print(f"  [muted]config_sha[/muted] [value]{man.get('config_sha')}[/value]")

    for rel in r["missing"]:
        bad(f"missing: {rel}")
    for rel in r["changed"]:
        bad(f"changed since the run: {rel}")
    for rel in r["unlisted"]:
        warn(f"present but not in the manifest: {rel}")

    n_bad = len(r["missing"]) + len(r["changed"])
    if not n_bad and not r["unlisted"]:
        ok(f"{r['ok']} file(s) verified - every one matches the manifest")
        return 0
    if not n_bad:
        # Extra files are not corruption: re-rendering a figure or adding a medium is
        # normal. They still mean the manifest is stale, so this is a failure for a
        # script that is checking a result is exactly what was published.
        warn(f"{r['ok']} file(s) match, {len(r['unlisted'])} not listed - the manifest is "
             f"out of date; re-run the stage that wrote them")
        return 1
    bad(f"{n_bad} file(s) missing or altered since the manifest was written - the results "
        f"in this directory are NOT the ones it describes")
    return 1


def main(argv=None) -> int:  # pragma: no cover - thin wrapper
    """`python -m dftfig.provenance [config.yaml]`, for scripting without the CLI."""
    from . import config as C
    path = (argv or sys.argv[1:] or ["config.yaml"])[0]
    return verify_report(C.load(path))


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
