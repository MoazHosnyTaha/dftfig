"""A finished calculation, bundled, so dftfig can be run without ORCA.

Every other way into this package needs a quantum chemistry program and hours of compute
before there is anything to look at. That is a poor way to evaluate a tool, and an
impossible one for anyone deciding whether it is worth installing ORCA in the first
place.

What ships here is a real run, not a mock: aniline and phenol at B3LYP/def2-SVP with a
def2-TZVP property single point, D3BJ, CPCM water, cubes on a 50-point grid. The ORCA
outputs are the genuine files that calculation produced, and the cubes are its real
volumetric data, gzipped to keep the package small.

The pair is chosen, not arbitrary. Aniline is a base and phenol is an acid, so the
speciation figure shows a curve rising and a curve falling on the same pH axis; and two
molecules is the smallest set that demonstrates what the figures are for, which is
comparing compounds on a shared scale. `dft demo` unpacks them into a working directory and
runs the two stages that need no quantum chemistry - descriptors and figures - against
them.

So the demo exercises the same code a real run does, from the same kind of inputs. The
only thing it skips is the part that needs ORCA, which has already happened.
"""
from __future__ import annotations

import gzip
import shutil
from pathlib import Path

from . import config as C
from . import params as P
from .theme import note, ok, warn

DATA = Path(__file__).parent / "demo_data"

# The level of theory the bundled outputs were actually run at. These are not defaults
# or suggestions: the figures stamp their footnotes from this config, so a mismatch here
# would print a method line that does not describe the data underneath it.
LEVEL = {
    "theory": {"functional": "B3LYP", "basis_opt": "def2-SVP", "dispersion": "D3BJ",
               "ri_approx": "RIJCOSX", "aux_basis": "def2/J", "scf_conv": "TightSCF",
               "grid": "DefGrid2", "run_freq": True},
    "property": {"enabled": True, "basis_sp": "def2-TZVP", "population": "Hirshfeld",
                 "koopmans_mode": "both"},
    "solvation": {"model": "CPCM", "solvents": ["water"]},
    "cubes": {"enabled": True, "grid_points": 50, "orbitals": "homo_lumo", "esp": True,
              "keep_cubes": True},
    "fukui": {"enabled": True},
}


def unpack(dest: Path) -> Path:
    """Lay the bundled run out under `dest` and write a config.yaml describing it."""
    dest = Path(dest).expanduser().resolve()
    (dest / "molecules").mkdir(parents=True, exist_ok=True)
    (dest / "data").mkdir(parents=True, exist_ok=True)

    for sdf in sorted(DATA.glob("*.sdf")):
        shutil.copy(sdf, dest / "molecules" / sdf.name)
    shutil.copy(DATA / "properties.csv", dest / "data" / "properties.csv")

    for src in sorted((DATA / "runs").rglob("*")):
        if not src.is_file():
            continue
        rel = src.relative_to(DATA / "runs")
        out = dest / "runs" / rel
        out.parent.mkdir(parents=True, exist_ok=True)
        if src.suffix == ".gz":
            out = out.with_suffix("")          # eldens.cube.gz -> eldens.cube
            if not out.exists():
                with gzip.open(src, "rb") as f, open(out, "wb") as g:
                    shutil.copyfileobj(f, g)
        else:
            shutil.copy(src, out)

    cfg = P.defaults()
    for section, values in LEVEL.items():
        cfg[section].update(values)
    cfg["molecules"].update({"structure_source": "sdf",
                             "sdf_dir": str(dest / "molecules"),
                             "sdf_pattern": "*.sdf", "name_from": "filename"})
    cfg["execution"].update({"work_dir": str(dest / "runs"),
                             "scratch_dir": str(dest / "scratch")})
    cfg["figures"].update({"reference_solvent": "water",
                           "properties_csv": str(dest / "data" / "properties.csv")})
    C.save(cfg, dest / "config.yaml")
    return dest


def run(dest: Path | str = "dftfig-demo") -> int:
    """Unpack the bundled run and build its descriptors and figures.

    Returns a process exit code: zero when at least the descriptor tables were written.
    """
    from . import __main__ as M
    from . import molecules as _M

    dest = unpack(Path(dest))
    note(f"bundled run unpacked into {dest}")
    note("aniline and phenol, B3LYP/def2-SVP with a def2-TZVP property single point, "
         "CPCM water - the ORCA outputs and cubes in there are the real ones those "
         "calculations produced.")

    cfg = C.load(dest / "config.yaml")
    specs = _M.load_molecules(cfg)
    df = M.stage_descriptors(cfg)
    M.stage_figures(cfg, specs, df)

    figures = sorted((dest / "runs" / "figures").glob("*.png"))
    if figures:
        ok(f"{len(figures)} figure(s) in {dest / 'runs' / 'figures'}")
        for f in figures:
            note(f"  {f.name}")
    else:
        warn("no figures were produced - see the reasons above. The descriptor tables "
             "below do not need PyMOL and should be there regardless.")

    for table in ("descriptors.csv", "atom_descriptors.csv"):
        path = dest / "runs" / table
        if path.exists():
            ok(f"{table} -> {path}")

    note(f"the config that produced all of this is {dest / 'config.yaml'} - edit it and "
         f"re-run any stage with `dft <stage> -c {dest / 'config.yaml'}`.")
    return 0 if df is not None and not df.empty else 1
