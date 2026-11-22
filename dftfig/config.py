"""Config loading, validation and merging."""
from __future__ import annotations

import copy
from pathlib import Path
from typing import Any

import yaml

from . import params as P


def deep_merge(base: dict, over: dict) -> dict:
    """Merge a user config over the defaults.

    The config is exactly two levels deep - section, then parameter - so the merge stops
    at the parameter. A parameter whose VALUE is itself a mapping (molecules.definitions)
    is replaced wholesale, never merged: a user who lists two molecules means two, not
    two plus whatever the defaults contained.
    """
    out = copy.deepcopy(base)
    for section, params in (over or {}).items():
        if isinstance(params, dict) and isinstance(out.get(section), dict):
            for key, value in params.items():
                out[section][key] = copy.deepcopy(value)
        else:
            out[section] = copy.deepcopy(params)
    return out


def load(path: str | Path) -> dict:
    """Load a YAML config, filling any missing key from the registry defaults."""
    user = yaml.safe_load(Path(path).read_text()) or {}
    cfg = deep_merge(P.defaults(), user)
    from .theme import note as _note, warn
    for found in resolve_executables(cfg):
        _note(found)
    for fixed in repair(cfg):
        warn(f"auto-corrected {fixed}")
    validate(cfg)
    return cfg


def save(cfg: dict, path: str | Path) -> None:
    Path(path).write_text(
        "# dftfig configuration\n"
        "# Every key is documented in the parameter registry: run `dftfig params`\n"
        "# or `python -m dftfig params --key theory.functional` for one entry.\n\n"
        + yaml.safe_dump(cfg, sort_keys=False, default_flow_style=False, width=100)
    )


def get(cfg: dict, key: str) -> Any:
    section, name = key.split(".", 1)
    return cfg[section][name]


class ConfigError(ValueError):
    pass


_ANNOUNCED: set = set()


def resolve_executables(cfg: dict) -> list[str]:
    """Fill in empty external-program paths by searching this machine.

    Kept out of `repair` because it touches the filesystem and should not run inside
    unit tests of pure config logic. Returns notes describing what was located.
    """
    from . import discover
    notes: list[str] = []
    # Several stages each load the config, so without this the same discovery line
    # is printed three or four times in one run.
    global _ANNOUNCED

    if not str(cfg["execution"].get("orca_path") or "").strip():
        found = discover.find_orca()
        if found:
            cfg["execution"]["orca_path"] = str(found)
            discover.remember("orca_path", str(found))
            if ("orca", str(found)) not in _ANNOUNCED:
                _ANNOUNCED.add(("orca", str(found)))
                notes.append(f"found ORCA at {found}")
    if not str(cfg["render"].get("pymol_path") or "").strip():
        found = discover.find_pymol()
        if found:
            cfg["render"]["pymol_path"] = str(found)
            discover.remember("pymol_path", str(found))
            if ("pymol", str(found)) not in _ANNOUNCED:
                _ANNOUNCED.add(("pymol", str(found)))
                notes.append(f"found PyMOL at {found}")
    return notes


def repair(cfg: dict) -> list[str]:
    """Fix mechanically-resolvable inconsistencies in place; return what was changed.

    These are combinations with exactly one sensible resolution, so asking the user to
    retype an answer would be busywork. Anything needing a real decision (an empty
    molecule set, a missing SDF folder) is left for `validate` to report.
    """
    changed: list[str] = []

    def set_(section, key, value, why):
        if cfg[section][key] != value:
            changed.append(f"{section}.{key}: {cfg[section][key]!r} -> {value!r}  ({why})")
            cfg[section][key] = value

    # Open MPI allocates one slot per PHYSICAL core. Asking for more does not run
    # slowly - it does not run at all, failing in Startup before any chemistry happens.
    # There is exactly one sensible resolution, so apply it rather than asking again.
    phys = physical_cores()
    nproc = int(cfg["execution"]["nprocs"])
    if nproc > phys:
        import os as _os
        set_("execution", "nprocs", phys,
             f"Open MPI offers one slot per physical core; this machine has {phys} "
             f"({_os.cpu_count()} hardware threads), and ORCA aborts at startup if asked "
             f"for more")

    func = str(cfg["theory"]["functional"])
    if cfg["theory"]["dispersion"] != "none" and (
            "D3BJ" in func or "D4" in func or "3c" in func):
        set_("theory", "dispersion", "none",
             f"{func} already includes dispersion; adding it twice is an error")

    if cfg["theory"]["ri_approx"] == "RIJK" and "/JK" not in str(cfg["theory"]["aux_basis"]):
        set_("theory", "aux_basis", "def2/JK", "RIJK needs a /JK auxiliary basis")

    if cfg["theory"]["run_freq"] and cfg["theory"]["opt_level"] == "LooseOpt":
        set_("theory", "opt_level", "Opt",
             "LooseOpt geometries give spurious imaginary frequencies")

    if cfg["solvation"]["model"] == "none" and cfg["solvation"]["solvents"] != ["gas"]:
        set_("solvation", "solvents", ["gas"], "no continuum model, so only gas phase applies")

    if cfg["figures"]["reference_solvent"] not in cfg["solvation"]["solvents"]:
        first = cfg["solvation"]["solvents"][0] if cfg["solvation"]["solvents"] else "water"
        set_("figures", "reference_solvent", first,
             "the figure solvent must be one that was actually calculated")

    if cfg["property"]["population"] == "Mulliken" and "+" in str(cfg["property"]["basis_sp"]):
        set_("property", "population", "Hirshfeld",
             "Mulliken charges are unreliable with diffuse functions")

    if (cfg["fukui"]["enabled"] or cfg["cubes"]["enabled"]) and not cfg["property"]["enabled"]:
        set_("property", "enabled", True,
             "Fukui functions and cube files are produced by the property single point")

    if "fig2" in cfg["figures"]["which"]:
        if not cfg["cubes"]["enabled"]:
            set_("cubes", "enabled", True, "Figure 2 panels B and C are built from cube files")
        if cfg["render"]["engine"] == "none":
            set_("render", "engine", "pymol", "Figure 2 needs a rendering engine")

    return changed


def discovered_figures(cfg: dict | None = None) -> list[str]:
    """Every figure name this installation can draw. Empty if discovery itself fails.

    Config validation must not be the thing that breaks because somebody's plugin has a
    syntax error - figureapi already warns about each one it could not load.
    """
    try:
        from . import figureapi
        return sorted(figureapi.available(cfg))
    except Exception:
        return []


def validate(cfg: dict) -> None:
    """Check values against the registry and catch combinations that cannot work."""
    errs: list[str] = []

    for key, p in P.ALL_PARAMS.items():
        section, name = key.split(".", 1)
        if section not in cfg or name not in cfg[section]:
            errs.append(f"{key}: missing")
            continue
        v = cfg[section][name]
        if p.kind == "choice":
            allowed = [c.value for c in p.choices]
            if allowed and v not in allowed:
                errs.append(f"{key}: {v!r} is not one of {allowed}")
        elif p.kind == "multichoice":
            allowed = [c.value for c in p.choices]
            if key == "figures.which":
                # The figure set is open: a plugin name is as legal as a built-in, and
                # validating against the registry rather than the static choice list is
                # what makes a user's own figure usable without editing dftfig.
                allowed = sorted(set(allowed) | set(discovered_figures(cfg)))
            if not isinstance(v, list):
                errs.append(f"{key}: expected a list, got {type(v).__name__}")
            elif allowed:
                bad = [x for x in v if x not in allowed]
                if bad:
                    errs.append(f"{key}: {bad} not in {allowed}")
                if not v:
                    errs.append(f"{key}: at least one value required")
        elif p.kind in ("int", "float"):
            if not isinstance(v, (int, float)) or isinstance(v, bool):
                errs.append(f"{key}: expected a number, got {v!r}")
            else:
                if p.minimum is not None and v < p.minimum:
                    errs.append(f"{key}: {v} is below the minimum {p.minimum}")
                if p.maximum is not None and v > p.maximum:
                    errs.append(f"{key}: {v} is above the maximum {p.maximum}")
        elif p.kind == "bool" and not isinstance(v, bool):
            errs.append(f"{key}: expected true/false, got {v!r}")

    # --- cross-parameter checks: the combinations that waste an afternoon -------------
    if cfg["figures"]["reference_solvent"] not in cfg["solvation"]["solvents"]:
        errs.append(
            f"figures.reference_solvent = {cfg['figures']['reference_solvent']!r} is not in "
            f"solvation.solvents = {cfg['solvation']['solvents']}. The figures would have no data."
        )
    if cfg["solvation"]["model"] == "none" and cfg["solvation"]["solvents"] != ["gas"]:
        errs.append(
            "solvation.model = 'none' but solvents other than 'gas' are requested. "
            "Either pick CPCM/SMD or set solvents to ['gas']."
        )
    func = str(cfg["theory"]["functional"])
    disp = cfg["theory"]["dispersion"]
    if disp != "none" and ("D3BJ" in func or "D4" in func or "3c" in func):
        errs.append(
            f"theory.functional = {func!r} already contains a dispersion correction and "
            f"theory.dispersion = {disp!r} would add a second one. Set dispersion to 'none'."
        )
    if cfg["theory"]["ri_approx"] == "RIJK" and "/JK" not in str(cfg["theory"]["aux_basis"]):
        errs.append("theory.ri_approx = 'RIJK' needs a /JK auxiliary basis, e.g. 'def2/JK'.")
    if cfg["theory"]["run_freq"] and cfg["theory"]["opt_level"] == "LooseOpt":
        errs.append(
            "theory.run_freq with opt_level = 'LooseOpt' will produce spurious imaginary "
            "frequencies. Use 'Opt' or 'TightOpt'."
        )
    if cfg["fukui"]["enabled"] and not cfg["property"]["enabled"]:
        errs.append("fukui.enabled requires property.enabled (the neutral reference single point).")
    if cfg["cubes"]["enabled"] and not cfg["property"]["enabled"]:
        errs.append("cubes.enabled requires property.enabled (cubes are written by the property job).")
    if "fig2" in cfg["figures"]["which"] and not cfg["cubes"]["enabled"]:
        errs.append("figures 'fig2' needs cubes.enabled (panels B and C are built from cube files).")
    if "fig2" in cfg["figures"]["which"] and cfg["render"]["engine"] == "none":
        errs.append("figures 'fig2' needs a render.engine other than 'none'.")
    if cfg["property"]["population"] == "Mulliken" and "+" in str(cfg["property"]["basis_sp"]):
        errs.append(
            "property.population = 'Mulliken' with a diffuse ('+') basis gives unreliable, "
            "sometimes negative atomic populations and will corrupt the Fukui indices. "
            "Use Hirshfeld or Loewdin."
        )
    if cfg["molecules"]["multiplicity"] < 1:
        errs.append("molecules.multiplicity must be >= 1.")
    source = cfg["molecules"]["structure_source"]
    if source == "sdf":
        d = Path(cfg["molecules"]["sdf_dir"]).expanduser()
        if not d.exists():
            errs.append(
                f"molecules.sdf_dir = {d} does not exist. Create it and put one .sdf per "
                f"compound inside (or run `dft init`, which scaffolds it)."
            )
        elif not sorted(d.glob(cfg["molecules"]["sdf_pattern"])):
            errs.append(
                f"no files matching {cfg['molecules']['sdf_pattern']!r} in {d} - "
                f"nothing to calculate."
            )
    elif not cfg["molecules"]["definitions"]:
        errs.append(
            f"molecules.structure_source = {source!r} but molecules.definitions is empty - "
            f"nothing to calculate. Add NAME: SMILES entries, or switch to 'sdf'."
        )

    if errs:
        raise ConfigError("Configuration problems:\n  - " + "\n  - ".join(errs))


def physical_cores() -> int:
    """Physical cores, not hardware threads.

    This is the number that matters for ORCA: Open MPI allocates one slot per physical
    core by default, so `%pal nprocs` above this count fails at startup with "not enough
    slots available" before a single integral is computed. os.cpu_count() reports threads
    (12 on a 6-core laptop with SMT), which is exactly the wrong number to offer.
    """
    import os
    import platform
    import subprocess

    system = platform.system()

    if system == "Darwin":
        # macOS has no /proc. hw.physicalcpu is exactly the number we want.
        try:
            out = subprocess.run(["sysctl", "-n", "hw.physicalcpu"],
                                 capture_output=True, text=True, timeout=5)
            n = int(out.stdout.strip())
            if n > 0:
                return n
        except Exception:
            pass

    if system == "Windows":
        # NUMBER_OF_PROCESSORS counts logical CPUs; wmic gives physical cores.
        try:
            out = subprocess.run(
                ["wmic", "cpu", "get", "NumberOfCores", "/value"],
                capture_output=True, text=True, timeout=10)
            vals = [int(l.split("=")[1]) for l in out.stdout.splitlines()
                    if l.strip().startswith("NumberOfCores=")]
            if vals:
                return sum(vals)
        except Exception:
            pass

    try:
        seen = set()
        pkg = core = None
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.startswith("physical id"):
                pkg = line.split(":")[1].strip()
            elif line.startswith("core id"):
                core = line.split(":")[1].strip()
            elif not line.strip():
                if pkg is not None and core is not None:
                    seen.add((pkg, core))
                pkg = core = None
        if pkg is not None and core is not None:
            seen.add((pkg, core))
        if seen:
            return len(seen)
    except Exception:
        pass
    return os.cpu_count() or 1


def machine() -> dict:
    """Cores and total RAM of the machine actually running the pipeline."""
    import os
    info = {"cores": os.cpu_count() or 1, "cores_physical": physical_cores(),
            "ram_gb": None}
    try:
        for line in Path("/proc/meminfo").read_text().splitlines():
            if line.startswith("MemTotal:"):
                info["ram_gb"] = int(line.split()[1]) / 1024 / 1024
                break
    except Exception:
        pass
    return info


def advisories(cfg: dict) -> list[str]:
    """Non-fatal warnings about settings that do not fit this machine.

    These are separate from `validate`: an oversubscribed machine still runs, it just
    runs badly or gets a job killed by the OOM killer halfway through an overnight batch.
    """
    out: list[str] = []
    m = machine()
    nproc = int(cfg["execution"]["nprocs"])
    njobs = int(cfg["execution"]["max_parallel_jobs"])
    maxcore = int(cfg["execution"]["maxcore_mb"])

    phys = m["cores_physical"]
    if nproc > phys:
        out.append(
            f"nprocs ({nproc}) exceeds the {phys} physical cores on this machine "
            f"({m['cores']} hardware threads). Open MPI allocates one slot per physical "
            f"core, so ORCA will abort in Startup with 'not enough slots available' "
            f"before it computes anything. Lower nprocs to {phys} or less."
        )

    cores_asked = nproc * njobs
    if cores_asked > phys:
        out.append(
            f"nprocs ({nproc}) x max_parallel_jobs ({njobs}) = {cores_asked} cores, but this "
            f"machine has {phys} physical cores. The jobs will contend and finish later "
            f"than if you asked for fewer."
        )

    if m["ram_gb"]:
        asked = maxcore * nproc * njobs / 1024
        budget = m["ram_gb"] * 0.75
        if asked > budget:
            fit = int(budget * 1024 / max(nproc * njobs, 1))
            out.append(
                f"maxcore_mb ({maxcore}) x nprocs ({nproc}) x max_parallel_jobs ({njobs}) asks "
                f"for {asked:.1f} GB, but this machine has {m['ram_gb']:.1f} GB total. ORCA "
                f"overshoots this target, so a job will likely be killed. Lower maxcore_mb to "
                f"about {fit} MB, or reduce nprocs / max_parallel_jobs."
            )

    if nproc > 16:
        out.append(
            f"nprocs = {nproc}: hybrid DFT with RIJCOSX stops scaling usefully past roughly "
            f"8-16 cores for a drug-sized molecule. Running more molecules concurrently "
            f"(max_parallel_jobs) converts the extra cores into throughput instead."
        )

    sp = str(cfg["property"]["basis_sp"])
    diffuse = ("++" in sp or sp.lower().startswith("aug-") or "TZVPD" in sp
               or "QZVPD" in sp)
    # Every cubes.orbitals setting in the registry plots the LUMO; the names checked here
    # were "frontier_pm1" and "all_valence", which no longer exist as choices, so the
    # advisory was silently skipped for homo_lumo_pm1 and homo_lumo_pm2 - the two
    # settings that plot the MOST unbound virtual orbitals.
    if (diffuse and cfg["cubes"]["enabled"]
            and cfg["cubes"]["orbitals"] in ("homo_lumo", "homo_lumo_pm1", "homo_lumo_pm2")):
        out.append(
            f"property.basis_sp = {sp} carries diffuse functions, and the LUMO of a neutral "
            f"organic molecule is usually unbound (positive orbital energy). With a diffuse "
            f"set that unbound state is described as a discretised continuum function, so "
            f"the plotted 'LUMO' can be a diffuse cloud detached from the molecule rather "
            f"than the pi* orbital you want to show. Use def2-TZVP for the orbital figure, "
            f"or keep the diffuse set and read the LUMO panel with that in mind."
        )

    from . import discover
    orca = str(cfg["execution"].get("orca_path") or "")
    if not orca:
        out.append(
            "no ORCA binary found on this machine. dftfig searched $ORCA_PATH, your PATH "
            "and the usual install locations. Install ORCA (free for academic use) and "
            "either put it on your PATH, set $ORCA_PATH, or set execution.orca_path to "
            "the absolute path of the binary. Everything that does not need ORCA - "
            "`dft params`, `dft validate`, and the figures stage on an existing run - "
            "works without it."
        )
    elif not discover.is_real_orca(orca):
        out.append(
            f"execution.orca_path = {orca} does not look like an ORCA installation: "
            f"orca_plot and orca_vpot are not beside it. On most Linux desktops "
            f"/usr/bin/orca is the GNOME Orca screen reader, not the quantum chemistry "
            f"program. Check the path, or run `dft doctor` to see the candidates found."
        )

    if cfg["molecules"]["charge_mode"] == "fixed":
        q = int(cfg["molecules"]["charge"])
        out.append(
            f"molecules.charge_mode = 'fixed' forces charge {q:+d} on every compound. If the "
            f"set is not uniformly charged, switch to 'auto' - a wrong charge changes the "
            f"electron count and makes every orbital and descriptor meaningless."
        )

    # An extra electron on an already-anionic species is often simply unbound without
    # diffuse functions, and the SCF then describes a continuum state rather than the
    # anion. This shows up as a negative electron affinity.
    if cfg["fukui"]["enabled"]:
        basis = str(cfg["property"]["basis_sp"])
        diffuse = ("+" in basis) or basis.lower().startswith(("aug-", "ma-"))
        if not diffuse:
            try:
                from . import molecules as _M
                specs = _M.load_molecules(cfg)
                anionic = [s.name for s in specs.values() if s.charge < 0]
            except Exception:
                anionic = []
            if anionic:
                verb = "carries" if len(anionic) == 1 else "carry"
                out.append(
                    f"{', '.join(anionic)} {verb} a negative charge, and the Fukui step adds "
                    f"one more electron. The property basis {basis} has no diffuse "
                    f"functions, so that extra electron is likely unbound - you will get a "
                    f"negative electron affinity describing a continuum state rather than a "
                    f"real anion. Use a basis with diffuse functions (6-311++G(d,p), "
                    f"ma-def2-TZVP) or read f+ for these species with caution."
                )

    est = estimate(cfg)
    if est["total_orca_jobs"] > 24:
        out.append(
            f"{est['total_orca_jobs']} ORCA jobs queued across "
            f"{est['molecules']} molecule(s) x {est['solvents']} medium(s). At roughly an hour "
            f"per optimization+frequency job for a 26-atom molecule, budget accordingly, or "
            f"trim solvation.solvents to the one medium your figures actually use."
        )
    return out


def estimate(cfg: dict) -> dict:
    """Rough job accounting so the cost of a settings choice is visible before launching."""
    from . import molecules as _M
    nmol = _M.count(cfg)
    nsolv = len(cfg["solvation"]["solvents"])
    opt = nmol * nsolv
    sp = nmol * nsolv if cfg["property"]["enabled"] else 0
    fuk = 2 * nmol * nsolv if cfg["fukui"]["enabled"] else 0
    n_orb = {"homo_lumo": 2, "homo_lumo_pm1": 4, "homo_lumo_pm2": 6}.get(cfg["cubes"]["orbitals"], 2)
    g = cfg["cubes"]["grid_points"]
    cubes_mb = (n_orb + 1) * nmol * nsolv * (g ** 3) * 13 / 1e6 if cfg["cubes"]["enabled"] else 0
    return {
        "molecules": nmol,
        "solvents": nsolv,
        "opt_freq_jobs": opt,
        "property_jobs": sp,
        "fukui_jobs": fuk,
        "total_orca_jobs": opt + sp + fuk,
        "cube_files": (n_orb + 1 + (1 if cfg["cubes"]["esp"] else 0)) * nmol * nsolv if cfg["cubes"]["enabled"] else 0,
        "cube_disk_mb": round(cubes_mb, 1),
        "peak_ram_gb": round(cfg["execution"]["maxcore_mb"] * cfg["execution"]["nprocs"]
                             * cfg["execution"]["max_parallel_jobs"] / 1024, 1),
    }
