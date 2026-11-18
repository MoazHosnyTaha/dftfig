"""ORCA input generation.

Three kinds of job are produced per molecule per medium:

  opt     OPT (+FREQ)  at theory.basis_opt         -> optimized geometry, minimum check
  sp      SP           at property.basis_sp        -> orbital energies, populations, cubes
  cat/ani SP on the sp geometry at N-1 / N+1       -> condensed Fukui functions, delta-SCF IE/EA
"""
from __future__ import annotations

from pathlib import Path


#: ORCA's internal names for the continuum solvents we expose.
SOLVENT_NAME = {
    "water": "Water", "methanol": "Methanol", "chloroform": "Chloroform",
    "dmso": "DMSO", "acetonitrile": "Acetonitrile", "toluene": "Toluene", "thf": "THF",
}
#: Dielectric constants, quoted in figure captions and used to order the solvent axis.
DIELECTRIC = {
    "gas": 1.0, "water": 78.36, "methanol": 32.61, "chloroform": 4.71,
    "dmso": 46.70, "acetonitrile": 36.60, "toluene": 2.37, "thf": 7.43,
}

POPULATION_KEY = {
    "Hirshfeld": "P_Hirshfeld", "Mulliken": "P_Mulliken", "Loewdin": "P_Loewdin",
    "Mayer": "P_Mayer",
}


def keyword_line(cfg: dict, stage: str, solvent: str, multiplicity: int = 1) -> str:
    """Assemble the ORCA '!' keyword line for one stage."""
    t = cfg["theory"]
    parts = [str(t["functional"])]

    if t["dispersion"] != "none":
        parts.append(str(t["dispersion"]))

    # r2SCAN-3c is a composite method that carries its own basis set.
    composite = "3c" in str(t["functional"])
    if not composite:
        parts.append(str(cfg["property"]["basis_sp"] if stage in ("sp", "cat", "ani")
                         else t["basis_opt"]))

    if t["ri_approx"] == "RIJCOSX":
        parts += ["RIJCOSX", str(t["aux_basis"])]
    elif t["ri_approx"] == "RIJK":
        parts += ["RIJK", str(t["aux_basis"])]

    if stage == "opt":
        parts.append(str(t["opt_level"]))
        if t["run_freq"]:
            parts.append("FREQ")
    else:
        parts.append("SP")

    parts.append(str(t["scf_conv"]))
    parts.append(str(t["grid"]))

    if cfg["solvation"]["model"] == "CPCM" and solvent != "gas":
        parts.append(f"CPCM({SOLVENT_NAME[solvent]})")

    # Open-shell species need unrestricted Kohn-Sham: the N+-1 ions always, and the
    # neutral too whenever its own multiplicity is not a singlet.
    if stage in ("cat", "ani") or multiplicity != 1:
        parts.append("UKS")

    return "! " + " ".join(parts)


def _blocks(cfg: dict, stage: str, solvent: str, spec=None) -> list[str]:
    out: list[str] = []
    ex = cfg["execution"]
    out.append(f"%pal nprocs {int(ex['nprocs'])} end")
    out.append(f"%maxcore {int(ex['maxcore_mb'])}")

    if cfg["solvation"]["model"] == "SMD" and solvent != "gas":
        out.append(f'%cpcm\n  smd true\n  SMDsolvent "{SOLVENT_NAME[solvent]}"\nend')

    if stage == "opt" and cfg["theory"]["run_freq"]:
        out.append(f"%freq  Temp {float(cfg['theory']['temperature'])}  end")

    if stage in ("sp", "cat", "ani"):
        pop = POPULATION_KEY[cfg["property"]["population"]]
        out.append(f"%output\n  Print[{pop}] 1\n  Print[P_Mulliken] 1\n  Print[P_Loewdin] 1\nend")

    if stage == "sp" and cfg["cubes"]["enabled"] and spec is not None:
        block = _plots_block(cfg, spec)
        if block:
            out.append(block)

    return out


def _plots_block(cfg: dict, spec) -> str:
    """%plots block requesting frontier-orbital and electron-density cubes.

    ORCA 6 has no ElPot entry in %plots - the electrostatic potential is evaluated
    separately by orca_vpot on the grid of the density cube (see cubes.py).
    """
    try:
        h = spec.homo_index()
    except ValueError:
        # Open shell: alpha and beta frontier orbitals differ, so a single HOMO/LUMO
        # pair is not well defined. The density cube is still produced below.
        from .theme import warn
        warn(f"{spec.name}: open-shell species (multiplicity {spec.multiplicity}) - "
             "frontier-orbital cubes skipped; the density and MEP are still generated")
        n = int(cfg["cubes"]["grid_points"])
        return ("%plots\n"
                f"  dim1 {n}\n  dim2 {n}\n  dim3 {n}\n  Format Gaussian_Cube\n"
                '  ElDens("eldens.cube");\nend')
    span = {"homo_lumo": 0, "homo_lumo_pm1": 1, "homo_lumo_pm2": 2}[cfg["cubes"]["orbitals"]]
    n = int(cfg["cubes"]["grid_points"])
    lines = ["%plots", f"  dim1 {n}", f"  dim2 {n}", f"  dim3 {n}", "  Format Gaussian_Cube"]
    for d in range(span, -1, -1):
        idx = h - d
        tag = "homo" if d == 0 else f"homo-{d}"
        lines.append(f'  MO("{tag}.cube", {idx}, 0);')
    for d in range(0, span + 1):
        idx = h + 1 + d
        tag = "lumo" if d == 0 else f"lumo+{d}"
        lines.append(f'  MO("{tag}.cube", {idx}, 0);')
    lines.append('  ElDens("eldens.cube");')
    lines.append("end")
    return "\n".join(lines)


def write_input(cfg: dict, spec, solvent: str, stage: str,
                geom_xyz: str | None = None, extra_keywords=(),
                basename: str | None = None) -> Path:
    """Write one ORCA input file for one molecule and return its path.

    Charge and multiplicity come from the molecule, not from the run, so a neutral and a
    cation in the same batch each get their own correct electron count.

    `extra_keywords` are appended to the '!' line and `basename` overrides the output file
    stem; both exist so a saddle-point restart can run a tighter optimization alongside the
    original without overwriting it.
    """
    name = spec.name
    work = Path(cfg["execution"]["work_dir"]) / name / solvent
    work.mkdir(parents=True, exist_ok=True)
    base = basename or f"{name}_{solvent}_{stage}"

    if stage == "cat":
        charge, mult = spec.ion(-1)      # one electron removed
    elif stage == "ani":
        charge, mult = spec.ion(+1)      # one electron added
    else:
        charge, mult = spec.charge, spec.multiplicity

    kw = keyword_line(cfg, stage, solvent, mult)
    for extra in extra_keywords:
        # Replace rather than duplicate: "! ... Opt ... TightOpt" is ambiguous, and so is
        # giving two grid keywords.
        family = {"Opt": ("Opt", "TightOpt", "LooseOpt", "VeryTightOpt"),
                  "SCF": ("NormalSCF", "TightSCF", "VeryTightSCF", "SloppySCF"),
                  "Grid": ("DefGrid1", "DefGrid2", "DefGrid3")}
        for members in family.values():
            if extra in members:
                kw = " ".join(w for w in kw.split()
                              if w not in members or w == extra)
                break
        if extra not in kw.split():
            kw += f" {extra}"
    body = [kw, ""]
    body += _blocks(cfg, stage, solvent, spec)
    body.append("")

    if geom_xyz is None:
        geom_xyz = f"{name}_start.xyz"
        src = Path(cfg["execution"]["work_dir"]) / name / f"{name}_start.xyz"
        if src.exists():
            (work / geom_xyz).write_text(src.read_text())
    body.append(f"* xyzfile {charge} {mult} {geom_xyz}")
    body.append("")

    path = work / f"{base}.inp"
    path.write_text("\n".join(body))
    return path


def stage_label(stage: str) -> str:
    return {"opt": "geometry optimization + frequencies",
            "sp": "property single point",
            "cat": "cation single point (N-1)",
            "ani": "anion single point (N+1)"}[stage]
