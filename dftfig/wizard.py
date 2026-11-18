"""Interactive, colour-coded parameter wizard.

Renders the registry in `params.py`: for every parameter it shows what the parameter does,
how to decide, what it costs, and - for categorical parameters - a table of every allowed
value with its own explanation and preset markers.
"""
from __future__ import annotations

import re
from pathlib import Path

from rich.align import Align
from rich.markup import escape
from rich.columns import Columns
from rich.panel import Panel
from rich.prompt import Confirm, FloatPrompt, IntPrompt, Prompt

from . import menu as M
from rich.rule import Rule
from rich.table import Table
from rich.text import Text

from . import config as C
from . import params as P
from .theme import PALETTE, console

TAG_MARK = {
    P.TAG_RECOMMENDED: ("suggested", "ok"),
    P.TAG_FAST: ("cheap", "warn"),
}


def _marks(tags) -> tuple:
    """Only the markers that still mean something: the paper tag is not shown anywhere.

    The registry still carries it on individual choices, but it is internal bookkeeping -
    the pipeline no longer offers to reproduce the manuscript, so advertising which value
    it used would only push a choice that has nothing to recommend it on its own.
    """
    return tuple(TAG_MARK[t][0] for t in tags if t in TAG_MARK)


def _tagstr(tags) -> Text:
    t = Text()
    for tag in tags:
        if tag in TAG_MARK:
            label, style = TAG_MARK[tag]
            if len(t):
                t.append(" ")
            t.append(label, style=style)
    return t


def _header() -> None:
    console.print()
    console.print(Panel(
        Align.center(Text.from_markup(
            "[brand]dftfig[/brand]  [muted]|[/muted]  automated DFT + figure pipeline for ORCA\n"
            "[muted]optimized geometries -> frontier orbitals -> MEP -> lipophilicity -> MPO[/muted]"
        )),
        border_style=PALETTE["brand"], padding=(1, 4),
    ))


# --------------------------------------------------------------------------- the steps
# The interview walks the pipeline in the order it actually executes, so the questions
# arrive in the order the calculation will use the answers.
STEP_PURPOSE = {
    "Molecules":
        "Reads your .sdf files and decides what each compound IS: its name, its charge, "
        "its spin state, and whether the coordinates in the file are usable as a starting "
        "geometry.",
    "Conformer search":
        "Only used when a structure has no usable 3D coordinates. RDKit generates "
        "candidate conformers, force-field minimizes them, and the lowest one is handed "
        "to ORCA as the starting point. It does not decide the final geometry - DFT does.",
    "Level of theory (geometry optimization)":
        "The slowest step and the one every later number inherits. ORCA relaxes each "
        "structure to a stationary point and then runs frequencies to prove it is a "
        "minimum. The functional and basis chosen here set the geometry, the "
        "thermochemistry, and the starting point for everything downstream.",
    "Implicit solvation":
        "The continuum model wrapped around the molecule. It shifts geometries slightly "
        "and energies a lot, and it is the single thing reviewers most often find "
        "mis-named in a methods section.",
    "Property single point":
        "A second, larger-basis calculation on the geometry that was just optimized. "
        "Orbital energies, the population analysis and the ionization/electron-affinity "
        "numbers come from here, not from the optimization.",
    "Fukui functions":
        "Two extra single points, at charge +1 and -1, on the same geometry. The "
        "difference in atomic populations gives the condensed Fukui indices - where the "
        "molecule is attacked by an electrophile, a nucleophile or a radical.",
    "Volumetric data (cube files)":
        "The grids the 3D pictures are drawn from: HOMO and LUMO, the total density, and "
        "the electrostatic potential mapped onto it. This is where Figure 2 comes from.",
    "3D rendering":
        "How those grids are turned into images - isovalues, colours, transparency, "
        "camera.",
    "Figure assembly":
        "Which figures are built and how the panels are laid out, labelled and scaled.",
    "Execution":
        "Cores, memory, where files are written, and whether a part-finished batch "
        "resumes or restarts.",
}

# Groups the user has already tuned to their liking: offered as one keep-or-change
# question rather than parameter by parameter.
FIGURE_GROUPS = ("3D rendering", "Figure assembly")

# A parameter that cannot take effect given the other answers is not asked about.
def _relevant(key: str, cfg: dict) -> bool:
    g = lambda k: C.get(cfg, k)
    rules = {
        "molecules.charge":            lambda: g("molecules.charge_mode") != "auto",
        "molecules.multiplicity":      lambda: g("molecules.multiplicity_mode") != "auto",
        "theory.imaginary_cutoff":     lambda: bool(g("theory.run_freq")),
        "theory.saddle_max_attempts":  lambda: bool(g("theory.run_freq")),
        "theory.saddle_step_angstrom": lambda: bool(g("theory.run_freq"))
                                               and g("theory.saddle_max_attempts") > 0,
        "theory.saddle_small_mode_cm": lambda: bool(g("theory.run_freq"))
                                               and g("theory.saddle_max_attempts") > 0,
        "theory.aux_basis":            lambda: g("theory.ri_approx") != "none",
        "property.basis_sp":           lambda: bool(g("property.enabled")),
        "property.population":         lambda: bool(g("property.enabled")),
        "property.koopmans_mode":      lambda: bool(g("property.enabled")),
        "fukui.geometry":              lambda: bool(g("fukui.enabled")),
        "fukui.smear_scheme":          lambda: bool(g("fukui.enabled")),
        "cubes.grid_points":           lambda: bool(g("cubes.enabled")),
        "cubes.orbitals":              lambda: bool(g("cubes.enabled")),
        "cubes.esp":                   lambda: bool(g("cubes.enabled")),
        "cubes.esp_nprocs":            lambda: bool(g("cubes.enabled")) and bool(g("cubes.esp")),
        "cubes.keep_cubes":            lambda: bool(g("cubes.enabled")),
        "molecules.xyz_dir":           lambda: g("molecules.structure_source") == "xyz",
        "molecules.sdf_dir":           lambda: g("molecules.structure_source") == "sdf",
        "molecules.sdf_pattern":       lambda: g("molecules.structure_source") == "sdf",
        "molecules.use_sdf_geometry":  lambda: g("molecules.structure_source") == "sdf",
    }
    try:
        return rules.get(key, lambda: True)()
    except Exception:                               # noqa: BLE001 - never block a question
        return True


def _step_map(titles: list[str], here: int) -> None:
    """Where this step sits in the pipeline, so the walk never feels open-ended."""
    t = Text()
    for i, name in enumerate(titles):
        mark = "> " if i == here else ("  " if i > here else "  ")
        style = (PALETTE["brand"] if i == here
                 else PALETTE["ok"] if i < here else PALETTE["muted"])
        t.append(f"  {mark}{i + 1}. {name}\n", style=style)
    console.print(Panel(t, title="[brand]the run, step by step[/brand]",
                        border_style=PALETTE["muted"], expand=False))


def _param_panel(p: P.Param, current) -> None:
    body = Text()
    body.append(p.help.strip() + "\n", style=PALETTE["info"])
    if p.guidance:
        body.append("\n")
        body.append("How to choose: ", style="bold")
        body.append(p.guidance.strip() + "\n", style="default")
    if p.cost:
        body.append("\n")
        body.append("Cost: ", style=PALETTE["warn"])
        body.append(p.cost.strip() + "\n", style=PALETTE["muted"])

    title = f"[value]{p.label}[/value] [muted]({p.key}"
    title += f", {p.unit}" if p.unit else ""
    title += ")[/muted]"
    console.print(Panel(body, title=title, title_align="left",
                        border_style=PALETTE["muted"], padding=(0, 2)))

    if p.choices:
        t = Table(show_header=True, header_style="brand", box=None, padding=(0, 2))
        t.add_column("value", style=PALETTE["value"], no_wrap=True)
        t.add_column("name")
        t.add_column("what it means")
        t.add_column("note")
        for c in p.choices:
            marker = "*" if (c.value == current or
                             (isinstance(current, list) and c.value in current)) else " "
            t.add_row(f"{marker} {escape(str(c.value))}", c.label, c.detail or "",
                      _tagstr(c.tags))
        console.print(t)
        console.print(f"  [muted]* = current value[/muted]")


def _ask(p: P.Param, current):
    _param_panel(p, current)
    if p.kind == "bool":
        return M.confirm(p.label, default=bool(current))
    if p.kind == "choice":
        opts = [M.Option(c.value, f"{c.value}   -   {c.label}", c.detail, _marks(c.tags))
                for c in p.choices]
        default = next((i for i, c in enumerate(p.choices) if c.value == current), 0)
        picked = M.select(p.label, opts, default=default, allow_quit=False)
        return current if picked is None else picked
    if p.kind == "multichoice":
        opts = [M.Option(c.value, f"{c.value}   -   {c.label}", c.detail, _marks(c.tags))
                for c in p.choices]
        while True:
            picked = M.select_many(p.label, opts, subtitle="space toggles, enter accepts",
                                   preselected=list(current))
            if picked:
                return picked
            console.print("  [bad]pick at least one[/bad]")
    if p.kind == "int":
        while True:
            v = IntPrompt.ask(f"  [brand]{p.label}[/brand]", default=int(current), console=console)
            if p.minimum is not None and v < p.minimum:
                console.print(f"  [bad]minimum is {p.minimum}[/bad]"); continue
            if p.maximum is not None and v > p.maximum:
                console.print(f"  [bad]maximum is {p.maximum}[/bad]"); continue
            return v
    if p.kind == "float":
        while True:
            v = FloatPrompt.ask(f"  [brand]{p.label}[/brand]", default=float(current), console=console)
            if p.minimum is not None and v < p.minimum:
                console.print(f"  [bad]minimum is {p.minimum}[/bad]"); continue
            if p.maximum is not None and v > p.maximum:
                console.print(f"  [bad]maximum is {p.maximum}[/bad]"); continue
            return v
    if p.kind == "kvlist":
        console.print("  [muted]current entries:[/muted]")
        for k, v in current.items():
            # escape(): a SMILES like [nH] or [NH+] is valid rich markup and would
            # otherwise be silently swallowed from the display.
            console.print(f"    [value]{escape(str(k))}[/value] = {escape(str(v))}")
        if not M.confirm("replace this molecule set?", default=False):
            return current
        out = {}
        console.print("  [muted]enter NAME=SMILES one per line; blank line to finish[/muted]")
        while True:
            raw = Prompt.ask("    entry", default="", console=console)
            if not raw.strip():
                break
            if "=" not in raw:
                console.print("    [bad]expected NAME=SMILES[/bad]"); continue
            name, smi = raw.split("=", 1)
            out[name.strip()] = smi.strip()
        return out or current
    # str / path
    return Prompt.ask(f"  [brand]{p.label}[/brand]", default=str(current), console=console)


def _summary(cfg: dict) -> None:
    est = C.estimate(cfg)
    t = Table(show_header=False, box=None, padding=(0, 2))
    t.add_column(style=PALETTE["muted"]); t.add_column(style=PALETTE["value"])
    t.add_row("molecules", str(est["molecules"]))
    t.add_row("media", str(est["solvents"]))
    t.add_row("OPT+FREQ jobs", str(est["opt_freq_jobs"]))
    t.add_row("property single points", str(est["property_jobs"]))
    t.add_row("Fukui single points", str(est["fukui_jobs"]))
    t.add_row("total ORCA jobs", f"[brand]{est['total_orca_jobs']}[/brand]")
    t.add_row("cube files", str(est["cube_files"]))
    t.add_row("cube disk", f"{est['cube_disk_mb']} MB")
    t.add_row("peak RAM request", f"{est['peak_ram_gb']} GB "
                                  f"({cfg['execution']['maxcore_mb']} MB x "
                                  f"{cfg['execution']['nprocs']} cores x "
                                  f"{cfg['execution']['max_parallel_jobs']} jobs)")
    console.print(Panel(t, title="[brand]what this will run[/brand]",
                        border_style=PALETTE["stage"], expand=False))

    kw = _keyword_preview(cfg)
    console.print(Panel(Text(kw, style=PALETTE["value"]),
                        title="[brand]ORCA keyword line (optimization step)[/brand]",
                        border_style=PALETTE["muted"], expand=False))

    mc = C.machine()
    console.print(f"  [muted]this machine:[/muted] [value]{mc['cores']} cores"
                  + (f", {mc['ram_gb']:.1f} GB RAM" if mc["ram_gb"] else "") + "[/value]")
    notes = C.advisories(cfg)
    for n in notes:
        console.print(f"  [warn]![/warn] {n}")
    if not notes:
        console.print("  [ok]these settings fit this machine[/ok]")


def _keyword_preview(cfg: dict) -> str:
    from .orca_io import keyword_line
    return keyword_line(cfg, stage="opt", solvent=cfg["figures"]["reference_solvent"])


def _sections_in(message: str) -> list[str]:
    """Which config sections an error message refers to, in order of appearance."""
    found: list[str] = []
    for m in re.finditer(r"\b([a-z_]+)\.([a-z_0-9]+)\b", message):
        key = f"{m.group(1)}.{m.group(2)}"
        if key in P.ALL_PARAMS and key not in found:
            found.append(key)
    return found


def _settle(cfg: dict) -> None:
    """Auto-repair what can be repaired, then re-ask only the parameters still at fault.

    The point is that nothing already answered is thrown away: a single bad combination
    should cost one prompt, not the whole wizard.
    """
    fixed = C.repair(cfg)
    if fixed:
        console.print(Panel("\n".join(fixed),
                            title="[ok]adjusted automatically[/ok]",
                            border_style=PALETTE["ok"], expand=False))

    for _ in range(6):
        try:
            C.validate(cfg)
            console.print("  [ok]configuration is valid[/ok]")
            return
        except C.ConfigError as e:
            msg = str(e)
            console.print(Panel(msg, title="[bad]still needs a decision[/bad]",
                                border_style=PALETTE["bad"], expand=False))

            # A missing folder the config itself names is worth offering to create.
            missing = re.search(r"(molecules\.sdf_dir|molecules\.xyz_dir) = (\S+) does not exist",
                                msg)
            if missing and M.confirm(f"create {missing.group(2)} now?",
                                       default=True, console=console):
                Path(missing.group(2)).expanduser().mkdir(parents=True, exist_ok=True)
                console.print(f"  [ok]created {missing.group(2)}[/ok] "
                              f"[muted]- put one .sdf per compound in it[/muted]")
                continue

            keys = _sections_in(msg)
            if not keys:
                break
            console.print(f"  [muted]re-asking only: {', '.join(keys)}[/muted]")
            for key in keys:
                p = P.ALL_PARAMS[key]
                section, name = key.split(".", 1)
                console.print()
                cfg[section][name] = _ask(p, cfg[section][name])
            C.repair(cfg)

    console.print("  [warn]saving anyway - run `dft validate` once you have fixed "
                  "the remaining items[/warn]")


def interview(cfg: dict, skip: tuple = ()) -> dict:
    """Walk the pipeline step by step, asking about every parameter that applies.

    There is no preset to pick from. The walk starts at whatever `cfg` already holds -
    the registry defaults for a new project - and every setting is put to the user with
    what it does, how to choose it, and what it costs. Parameters that cannot take
    effect given the other answers are not asked about (`_relevant`), and the two
    figure-styling groups are offered as one keep-or-change question because their
    current values are a deliberate look rather than a scientific choice.

    `skip` names keys the caller has already settled - `dft run` resolves the structure
    source and the compound selection before the interview starts, and re-asking them
    would invalidate the detection it just showed.
    """
    show_advanced = M.confirm(
        "Include the advanced settings?   (scratch paths, convergence tolerances, "
        "iteration limits - sensible defaults, rarely touched)",
        default=False)

    titles = [g.title for g in P.GROUPS]
    _step_map(titles, -1)

    for gi, g in enumerate(P.GROUPS):
        visible = [p for p in g.params
                   if (show_advanced or not p.advanced)
                   and p.key not in skip
                   and _relevant(p.key, cfg)]
        if not visible:
            continue

        console.print()
        console.print(Rule(f"[brand]step {gi + 1} of {len(titles)}  -  {g.title}[/brand]",
                           style=PALETTE["brand"]))
        purpose = STEP_PURPOSE.get(g.title, g.blurb)
        console.print(Panel(Text(purpose, style=PALETTE["info"]),
                            border_style=PALETTE["stage"], padding=(0, 2)))

        # The figure look is something the user has already settled; do not make them
        # walk 17 styling questions to get to the next calculation step.
        if g.title in FIGURE_GROUPS:
            cur = Table(show_header=False, box=None, padding=(0, 2))
            cur.add_column(style=PALETTE["muted"]); cur.add_column(style=PALETTE["value"])
            for p in visible:
                v = C.get(cfg, p.key)
                cur.add_row(p.label, escape(str(v)))
            console.print(cur)
            choice = M.select(
                g.title,
                [M.Option("keep", "keep these as they are",
                          "the figure style you already have - nothing changes",
                          ("suggested",)),
                 M.Option("edit", f"go through all {len(visible)} settings",
                          "isovalues, colours, transparency, DPI, panel labels")],
                default=0, allow_quit=False)
            if choice == "keep":
                continue

        console.print(f"  [muted]{len(visible)} setting(s) in this step. "
                      f"The value already in use is pre-selected - press Enter to keep "
                      f"it.[/muted]")

        for pi, p in enumerate(visible):
            if not _relevant(p.key, cfg):     # an earlier answer in this step may have
                continue                      # just made this one moot
            console.print()
            console.print(f"  [stage]step {gi + 1}.{pi + 1}[/stage]  "
                          f"[muted]({pi + 1} of {len(visible)} in {g.title})[/muted]")
            section, name = p.key.split(".", 1)
            cfg[section][name] = _ask(p, cfg[section][name])

    return cfg


def run_wizard(out_path: str = "config.yaml", start: str | None = None,
               base: dict | None = None) -> dict:
    """Walk the whole registry interactively and write a config file.

    `base` starts the walk from a config that is already in hand - `dft run` passes the
    one it has just assembled, so the review step edits that rather than resetting to a
    preset the user did not pick.
    """
    _header()
    if base is not None:
        cfg = base
    else:
        # Editing settings means editing the settings that are there, not resetting to
        # the registry defaults and making the user re-enter everything they chose.
        try:
            cfg = C.load(out_path)
            console.print(f"  [ok]editing:[/ok] [value]{Path(out_path).resolve()}[/value]")
        except Exception:                           # noqa: BLE001 - absent or unreadable
            cfg = P.defaults()
            console.print("  [muted]no usable config here - starting from the "
                          "defaults[/muted]")
    interview(cfg, skip=())
    console.print()
    _settle(cfg)

    _summary(cfg)
    C.save(cfg, out_path)
    console.print(f"\n  [ok]written:[/ok] [value]{Path(out_path).resolve()}[/value]")
    console.print("  [muted]run it with:[/muted] [value]python -m dftfig all -c "
                  f"{out_path}[/value]")
    console.print("  [muted]or one stage at a time:[/muted] "
                  "[value]build | opt | sp | fukui | cubes | descriptors | render | figures[/value]\n")
    return cfg


def show_params(key: str | None = None, advanced: bool = True) -> None:
    """Print the parameter reference without prompting for anything."""
    _header()
    if key:
        p = P.ALL_PARAMS.get(key)
        if not p:
            console.print(f"  [bad]no such parameter:[/bad] {key}")
            matches = [k for k in P.ALL_PARAMS if key in k]
            if matches:
                console.print("  [muted]did you mean:[/muted] " + ", ".join(matches))
            return
        _param_panel(p, p.default)
        console.print(f"  [muted]default:[/muted] [value]{p.default}[/value]")
        return
    for g in P.GROUPS:
        console.print()
        console.print(Rule(f"[brand]{g.title}[/brand]", style=PALETTE["brand"]))
        console.print(f"  [muted]{g.blurb}[/muted]\n")
        for p in g.params:
            if p.advanced and not advanced:
                continue
            _param_panel(p, p.default)
            console.print(f"  [muted]default:[/muted] [value]{p.default}[/value]\n")
