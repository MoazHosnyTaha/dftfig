"""The `dft` command: an interactive front door to the pipeline.

Running `dft` with no arguments opens a menu. Every menu entry maps to a subcommand,
so anything you can click you can also script:

    dft                  interactive menu
    dft wizard           build or edit config.yaml
    dft run              run every stage
    dft run opt sp       run selected stages
    dft status           what is done and what is missing
    dft params [KEY]     the parameter reference
    dft figures          rebuild the figures only
    dft init             scaffold config.yaml, molecules/ and data/ in this folder
    dft where            show the installation and config paths
    dft doctor           check ORCA, PyMOL, cores and packages on this machine
"""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

from rich.align import Align
from rich.panel import Panel
from rich.prompt import Confirm, Prompt
from rich.table import Table
from rich.text import Text

from . import __version__
from . import config as C
from .theme import PALETTE, bad, console, ok, warn

STAGES = ["build", "opt", "sp", "fukui", "cubes", "descriptors", "figures"]

STAGE_BLURB = {
    "build": "read the structures, assign charge and multiplicity, make 3D geometries",
    "opt": "geometry optimization and frequencies (the slow one)",
    "sp": "single point at the property basis, plus the cation and anion",
    "fukui": "condensed Fukui indices from the three charge states",
    "cubes": "orbital, density and electrostatic-potential grids",
    "descriptors": "reactivity indices, the warnings column, descriptors.csv",
    "figures": "render figures 2, 12, 13 and 14 - no DFT, safe to repeat",
}

MENU = [
    ("1", "Run the pipeline", "guided",
     "Set up, pick the compounds, check the cost, run - all in one pass"),
    ("2", "Run one stage", "stage",
     "A single stage - useful for re-rendering figures without redoing the DFT"),
    ("3", "Settings", "wizard",
     "Walk every parameter interactively and write config.yaml"),
    ("4", "Status", "status",
     "Which compounds are done, which are missing, which are not minima"),
    ("5", "Check settings", "validate",
     "Validate config.yaml and print the job, disk and memory estimate"),
    ("6", "Parameter reference", "params",
     "What every parameter does and how to choose"),
    ("7", "Open results folder", "open",
     "Show where the figures, cubes and tables were written"),
    ("8", "Check this machine", "doctor",
     "ORCA, PyMOL, cores, memory and Python packages - run this first if something fails"),
    ("q", "Quit", "quit", ""),
]


def _config_path(explicit: str | None = None) -> Path:
    """config.yaml in the current directory, unless one was named explicitly."""
    if explicit:
        return Path(explicit).expanduser()
    env = os.environ.get("DFTFIG_CONFIG")
    return Path(env).expanduser() if env else Path.cwd() / "config.yaml"


def _banner(cfg_path: Path) -> None:
    exists = cfg_path.exists()
    state = (f"[ok]config:[/ok] [value]{cfg_path}[/value]" if exists
             else f"[warn]no config yet[/warn] [muted]({cfg_path})[/muted]")
    console.print()
    console.print(Panel(
        Align.center(Text.from_markup(
            f"[brand]dft[/brand]  [muted]v{__version__}[/muted]  "
            f"[muted]|[/muted]  automated DFT + figure pipeline for ORCA\n"
            f"{state}\n"
            f"[muted]working directory: {Path.cwd()}[/muted]"
        )),
        border_style=PALETTE["brand"], padding=(1, 3),
    ))


def _menu(cfg_path: Path) -> str:
    from . import menu as M
    picked = M.select("What would you like to do?",
                      [M.Option(key, label, blurb) for key, label, _, blurb in MENU],
                      default=0, allow_quit=True)
    return picked or "q"


def _run_stage(cfg_path: Path, stages: list[str]) -> int:
    from .__main__ import main as pipeline_main
    try:
        for s in stages:
            rc = pipeline_main([s, "-c", str(cfg_path)])
            if rc:
                return rc
    except KeyboardInterrupt:
        # runner.run_many has already stopped ORCA and explained how to resume;
        # no traceback is useful here.
        return 130
    return 0


def _ensure_config(cfg_path: Path) -> bool:
    """Offer to create a config when none exists. Returns True if one is present."""
    if cfg_path.exists():
        return True
    warn(f"no configuration at {cfg_path}")
    if not Confirm.ask("  [brand]create one now?[/brand]", default=True, console=console):
        return False
    from .wizard import run_wizard
    run_wizard(str(cfg_path))
    return cfg_path.exists()


def _open_results(cfg_path: Path) -> None:
    if not cfg_path.exists():
        warn("no config, so no results directory is defined yet")
        return
    cfg = C.load(str(cfg_path))
    work = Path(cfg["execution"]["work_dir"]).resolve()
    figs = work / "figures"
    t = Table(show_header=False, box=None, padding=(0, 2))
    t.add_column(style=PALETTE["muted"]); t.add_column(style=PALETTE["value"])
    t.add_row("results", str(work))
    t.add_row("figures", str(figs))
    t.add_row("descriptors", str(work / "descriptors.csv"))
    console.print(t)
    if figs.exists():
        made = sorted(p.name for p in figs.iterdir() if p.suffix in (".png", ".pdf", ".svg"))
        console.print(f"  [muted]figures present:[/muted] "
                      f"[value]{', '.join(made) if made else 'none yet'}[/value]")
    opener = "xdg-open" if sys.platform.startswith("linux") else "open"
    if figs.exists() and Confirm.ask("  [brand]open the figures folder?[/brand]",
                                     default=False, console=console):
        subprocess.Popen([opener, str(figs)], stdout=subprocess.DEVNULL,
                         stderr=subprocess.DEVNULL)


def _where() -> None:
    from . import __file__ as pkg_file
    from . import params as P
    t = Table(show_header=False, box=None, padding=(0, 2))
    t.add_column(style=PALETTE["muted"]); t.add_column(style=PALETTE["value"])
    t.add_row("dft version", __version__)
    t.add_row("package", str(Path(pkg_file).parent))
    t.add_row("python", sys.executable)
    t.add_row("working directory", str(Path.cwd()))
    t.add_row("config it would use", str(_config_path()))
    console.print(Panel(t, title="[brand]installation[/brand]",
                        border_style=PALETTE["muted"], expand=False))
    try:
        from .render import find_pymol
        exe = find_pymol(P.defaults())
        (ok if exe else warn)(f"PyMOL: {exe or 'not found'}")
    except Exception as e:  # pragma: no cover - diagnostic only
        warn(f"PyMOL check failed: {e}")
    from . import discover
    found = discover.find_orca()
    if found:
        ok(f"ORCA: {found}")
    else:
        bad("ORCA: not found. Install it (free for academic use), then put it on your "
            "PATH, set $ORCA_PATH, or set execution.orca_path.")


def _doctor() -> int:
    """Report what this machine provides, so a new user can diagnose their own setup."""
    import platform
    import sys as _sys
    from . import discover
    from . import config as C

    console.print()
    t = Table(show_header=False, box=None, padding=(0, 2))
    t.add_column(style=PALETTE["muted"]); t.add_column(style=PALETTE["value"])
    t.add_row("dftfig", __version__)
    t.add_row("python", f"{platform.python_version()}  ({_sys.executable})")
    t.add_row("platform", f"{platform.system()} {platform.release()}  {platform.machine()}")
    m = C.machine()
    t.add_row("cores", f"{m['cores_physical']} physical / {m['cores']} logical"
                       "   [MPI slots are allocated per physical core]")
    if m["ram_gb"]:
        t.add_row("memory", f"{m['ram_gb']:.1f} GB")
    console.print(Panel(t, title="[brand]this machine[/brand]",
                        border_style=PALETTE["muted"], expand=False))

    rep = discover.report()
    if rep["orca"]:
        ok(f"ORCA: {rep['orca']}")
        ver = discover.orca_version(rep["orca"])
        console.print(f"      [muted]version[/muted] [value]{ver or 'could not be determined'}[/value]")
    else:
        bad("ORCA: not found")
        console.print("      [muted]dftfig checked $ORCA_PATH, your PATH, and ~/orca*, "
                      "/opt/orca*, /usr/local/orca*, /share/apps/orca*.\n"
                      "      A candidate counts only if orca_plot and orca_vpot sit beside "
                      "it - on most Linux\n      desktops /usr/bin/orca is the GNOME screen "
                      "reader, not the chemistry program.[/muted]")
    if len(rep["orca_candidates"]) > 1:
        console.print("      [muted]other candidates:[/muted] "
                      + ", ".join(rep["orca_candidates"][1:]))

    (ok if rep["pymol"] else warn)(
        f"PyMOL: {rep['pymol']}" if rep["pymol"] else
        "PyMOL: not found - 3D panels will be skipped. "
        "conda install -c conda-forge pymol-open-source")

    missing = []
    for mod in ("rdkit", "numpy", "scipy", "pandas", "matplotlib", "yaml", "PIL", "rich"):
        try:
            __import__(mod)
        except ImportError:
            missing.append(mod)
    (bad if missing else ok)(
        f"python packages: missing {', '.join(missing)}" if missing
        else "python packages: all present")

    console.print(f"\n  [muted]remembered paths:[/muted] [value]{rep['user_config']}[/value]\n")
    return 0 if (rep["orca"] and not missing) else 1


def _init(cfg_path: Path) -> int:
    """Create an editable working directory: config.yaml, molecules/ and data templates."""
    tpl = Path(__file__).parent / "templates"
    data = Path.cwd() / "data"
    (data / "som").mkdir(parents=True, exist_ok=True)
    moldir = Path.cwd() / "molecules"
    moldir.mkdir(parents=True, exist_ok=True)

    made = []
    hint = moldir / "PUT_YOUR_SDF_FILES_HERE.txt"
    if not hint.exists() and not any(moldir.glob("*.sdf")):
        hint.write_text(
            "Put one .sdf file per compound in this folder.\n"
            "\n"
            "  * The file name becomes the molecule name, which becomes the output\n"
            "    directory and the figure panel label. Keep it short.\n"
            "  * Each compound's charge is read from the formal charges in the file,\n"
            "    so a neutral and a cation can run side by side in the same batch.\n"
            "  * If a file already contains 3D coordinates they are used as the starting\n"
            "    geometry; otherwise a conformer is generated.\n"
            "  * A multi-record SDF also works - each record becomes its own molecule.\n"
            "\n"
            "Then:  dft validate   (check what was detected)\n"
            "       dft run        (run the pipeline)\n"
            "\n"
            "Delete this file once you have added your compounds.\n"
        )
        made.append(str(hint))
    target = data / "properties.csv"
    if target.exists():
        warn(f"{target} already exists - left untouched")
    else:
        # Headers only. The pipeline derives pKa, logD, TPSA, MW, HBD and clogP from
        # the structures; this file exists only so a measured value can override one.
        header = (tpl / "properties.csv").read_text().splitlines()[0]
        target.write_text(
                "# Optional. Everything here overrides what the pipeline works out from\n"
                "# the structures. Leave it empty and pKa, logD, TPSA, MW, HBD and clogP\n"
                "# are all derived; fill a cell in and that value is used instead.\n"
                "# Permeability, clearance, PPB, logBB, HIA and P-gp have no structural\n"
                "# estimate - supply them from an ADMET model or leave the panels blank.\n"
                "# One row per compound; the molecule column must match the SDF file name.\n"
                + header + "\n")
        made.append(str(target))
    readme = data / "som" / "README.txt"
    if not readme.exists():
        readme.write_bytes((tpl / "som_README.txt").read_bytes())
        made.append(str(readme))

    if cfg_path.exists():
        warn(f"{cfg_path} already exists - left untouched")
    else:
        from . import params as P
        cfg = P.defaults()
        cfg["molecules"]["structure_source"] = "sdf"
        cfg["molecules"]["sdf_dir"] = "molecules"
        cfg["theory"]["functional"] = "B3LYP"
        cfg["theory"]["basis_opt"] = "6-31+G(d,p)"
        cfg["property"]["basis_sp"] = "def2-TZVP"
        # One medium is what the figures are built from; five is five times the compute
        # for supplementary panels. Add them back with:
        #   solvation.solvents: [water, methanol, chloroform, dmso, gas]
        cfg["solvation"]["solvents"] = ["water"]
        cfg["property"]["koopmans_mode"] = "both"
        C.save(cfg, cfg_path)
        made.append(str(cfg_path))

    # Advisories are best-effort here. A freshly scaffolded directory has an empty
    # molecules/ folder, so full validation necessarily fails - that is the user's next
    # step, not an error in what init just created.
    try:
        for n in C.advisories(C.load(str(cfg_path))):
            warn(n)
    except C.ConfigError:
        pass

    for m in made:
        ok(f"created {m}")
    mdir = Path(cfg_path).parent / "molecules"
    n_sdf = len(list(mdir.glob("*.sdf"))) if mdir.exists() else 0
    if not n_sdf:
        console.print(f"  [muted]next:[/muted] put one [value].sdf[/value] per compound in "
                      f"[value]{mdir}[/value][muted], then run[/muted] [value]dft run[/value]\n")
        console.print("  [muted]structures are the only required input: charge, pKa, logD, "
                      "TPSA and the\n  CNS MPO score are all worked out from them.[/muted]\n")
    else:
        console.print(f"  [muted]next:[/muted] [value]dft run[/value] "
                      f"[muted]- it walks you through every setting, then runs[/muted]\n")
    return 0


def interactive(cfg_path: Path) -> int:
    while True:
        _banner(cfg_path)
        choice = _menu(cfg_path)
        action = dict((m[0], m[2]) for m in MENU)[choice]

        if action == "quit":
            console.print("  [muted]bye[/muted]\n")
            return 0
        if action == "guided":
            guided(cfg_path)
        elif action == "wizard":
            from .wizard import run_wizard
            run_wizard(str(cfg_path))
        elif action == "run":
            if _ensure_config(cfg_path):
                _run_stage(cfg_path, STAGES)
        elif action == "stage":
            if _ensure_config(cfg_path):
                from . import menu as M
                s = M.select("Which stage?",
                             [M.Option(x, x, STAGE_BLURB.get(x, "")) for x in STAGES],
                             default=STAGES.index("figures"))
                if s:
                    _run_stage(cfg_path, [s])
        elif action == "status":
            if _ensure_config(cfg_path):
                _run_stage(cfg_path, ["status"])
        elif action == "validate":
            if _ensure_config(cfg_path):
                _run_stage(cfg_path, ["validate"])
        elif action == "params":
            from .wizard import show_params
            key = Prompt.ask("  [brand]parameter key[/brand] "
                             "[muted](blank = all)[/muted]", default="", console=console)
            show_params(key or None)
        elif action == "doctor":
            _doctor()
        elif action == "open":
            _open_results(cfg_path)
        elif action == "init":
            _init(cfg_path)

        console.print()
        from . import menu as M
        if not M.confirm("Back to the menu?", default=True):
            return 0


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)

    cfg_path = _config_path()
    if "-c" in argv:
        i = argv.index("-c")
        if i + 1 < len(argv):
            cfg_path = _config_path(argv[i + 1])
            del argv[i:i + 2]
    elif "--config" in argv:
        i = argv.index("--config")
        if i + 1 < len(argv):
            cfg_path = _config_path(argv[i + 1])
            del argv[i:i + 2]

    if not argv:
        try:
            return guided(cfg_path)
        except (KeyboardInterrupt, EOFError):
            console.print("\n  [muted]interrupted[/muted]\n")
            return 130

    cmd = argv[0]
    rest = argv[1:]

    if cmd in ("-h", "--help", "help"):
        console.print(Panel(Text(__doc__ or "", style=PALETTE["info"]),
                            title="[brand]dft[/brand]",
                            border_style=PALETTE["brand"], expand=False))
        return 0
    if cmd in ("-V", "--version", "version"):
        console.print(f"dft {__version__}")
        return 0
    if cmd == "doctor":
        return _doctor()
    if cmd == "where":
        _where()
        return 0
    if cmd == "init":
        return _init(cfg_path)
    if cmd == "menu":
        return interactive(cfg_path)
    if cmd == "wizard":
        from .wizard import run_wizard
        run_wizard(str(cfg_path))
        return 0
    if cmd == "params":
        from .wizard import show_params
        show_params(rest[0] if rest else None)
        return 0
    if cmd == "run":
        stages = [s for s in rest if s in STAGES]
        if stages:                               # explicit stages: no questions asked
            if not _ensure_config(cfg_path):
                return 2
            return _run_stage(cfg_path, stages)
        try:
            return guided(cfg_path)
        except (KeyboardInterrupt, EOFError):
            console.print("\n  [muted]interrupted[/muted]\n")
            return 130
    if cmd in STAGES or cmd in ("status", "validate"):
        if cmd not in ("status", "validate") and not _ensure_config(cfg_path):
            return 2
        return _run_stage(cfg_path, [cmd])

    bad(f"unknown command: {cmd}")
    console.print("  [muted]try:[/muted] [value]dft --help[/value]  "
                  "[muted]or just[/muted] [value]dft[/value] [muted]for the menu[/muted]")
    return 2


# ----------------------------------------------------------------------- guided run
def guided(cfg_path: Path) -> int:
    """`dft run` with no arguments: set up, choose, check, and run - in one pass.

    Every question has a sensible default on the arrow keys, so the whole pipeline is
    enter-enter-enter for someone who just wants the figures.
    """
    from . import config as C
    from . import menu as M

    _banner(cfg_path)

    # 1. a place to work -------------------------------------------------------------
    if not cfg_path.exists():
        console.print("  [muted]no config here yet - setting the folder up first[/muted]\n")
        _init(cfg_path)

    moldir = Path(cfg_path).parent / "molecules"
    try:
        cfg = C.load(str(cfg_path))
    except C.ConfigError:
        cfg = None
    if cfg is not None:
        moldir = Path(cfg["molecules"]["sdf_dir"]).expanduser()
        if not moldir.is_absolute():
            moldir = Path(cfg_path).parent / moldir

    # 2. structures ------------------------------------------------------------------
    sdfs = sorted(moldir.glob(cfg["molecules"]["sdf_pattern"])) if cfg else []
    if not sdfs and (cfg is None or cfg["molecules"]["structure_source"] == "sdf"):
        bad(f"no structures found in {moldir}")
        console.print(
            "\n  Put one [value].sdf[/value] file per compound in that folder and run "
            "[value]dft run[/value] again.\n"
            "  [muted]The file name becomes the compound name. Charge is read from the "
            "formal charges in\n  the file, so neutrals, cations and anions can run "
            "together in one batch. Nothing else\n  is required: pKa, logD, TPSA, MW and "
            "the CNS MPO score are all derived from the structure.[/muted]\n")
        return 2

    # 3. what was detected -------------------------------------------------------------
    from . import molecules as Mol
    try:
        specs = Mol.load_molecules(cfg)
    except Exception as exc:                       # noqa: BLE001 - shown, not swallowed
        bad(f"could not read the structures: {exc}")
        return 2
    Mol.report(specs)

    which = M.select_many(
        "Which compounds should run?",
        [M.Option(n, n, f"{s.formula}   charge {s.charge:+d}")
         for n, s in specs.items()],
        subtitle="space toggles, enter accepts - all are selected",
        preselected=list(specs))
    if not which:
        console.print("  [muted]nothing selected[/muted]")
        return 0
    if set(which) != set(specs):
        cfg["molecules"]["only"] = ", ".join(which)

    # 4. how to run it -----------------------------------------------------------------
    # No presets. Every setting is put to the user, step by step, in the order the
    # pipeline will use it. The structure source and the compound selection were just
    # resolved above, so those keys are not re-asked - changing them here would
    # invalidate the detection the user has already seen and approved.
    from . import wizard
    console.print()
    console.print(Panel(Text(
        "Now the settings, one step at a time, in the order the pipeline runs them.\n"
        "Each one shows what it does, how to choose it and what it costs. The value "
        "already in use is pre-selected, so Enter keeps it and arrow keys (or the "
        "number) change it.",
        style=PALETTE["info"]),
        title="[brand]settings[/brand]", border_style=PALETTE["brand"], expand=False))
    cfg = wizard.interview(cfg, skip=(
        "molecules.structure_source", "molecules.sdf_dir", "molecules.sdf_pattern",
        "molecules.xyz_dir", "molecules.name_from", "molecules.only",
        "molecules.limit", "molecules.definitions"))
    for n in C.repair(cfg):
        warn(n)
    C.save(cfg, cfg_path)

    for n in C.advisories(cfg):
        warn(n)

    # 5. what it will cost -------------------------------------------------------------
    _print_estimate(cfg)
    if not M.confirm("Start the run?", default=True):
        console.print("  [muted]nothing run; the settings are saved[/muted]")
        return 0

    # 6. go ----------------------------------------------------------------------------
    rc = _run_stage(cfg_path, list(STAGES))
    if rc == 0:
        _show_outputs(cfg)
    return rc


def _print_estimate(cfg) -> None:
    from . import config as C
    est = C.estimate(cfg)
    t = Table(show_header=False, box=None, padding=(0, 2))
    t.add_column(style=PALETTE["muted"]); t.add_column(style=PALETTE["value"])
    for k in ("molecules", "solvents", "total_orca_jobs", "cube_files",
              "cube_disk_mb", "peak_ram_gb"):
        if k in est:
            t.add_row(k.replace("_", " "), str(est[k]))
    console.print()
    console.print(Panel(t, title="[brand]what this will cost[/brand]",
                        border_style=PALETTE["muted"], expand=False))


def _show_outputs(cfg) -> None:
    work = Path(cfg["execution"]["work_dir"])
    figs = work / "figures"
    console.print()
    ok(f"results in {work.resolve()}")
    if figs.exists():
        for f in sorted(figs.glob("*.png")):
            console.print(f"    [value]{f.name}[/value]")
    console.print(f"    [value]descriptors.csv[/value] "
                  f"[muted]- every number behind the figures, with a warnings column"
                  f"[/muted]\n")


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        from . import runner as _r
        _r.stop_all()
        console.print("\n  [muted]stopped[/muted]\n")
        sys.exit(130)
