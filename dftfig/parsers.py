"""Parse ORCA output files.

Every function takes the text of a .out file (or a path) and is tolerant of missing
sections - a job that crashed half way still yields whatever it managed to print.
"""
from __future__ import annotations

import re
from pathlib import Path

HARTREE_TO_EV = 27.211386245988
HARTREE_TO_KCAL = 627.509474


def _text(src) -> str:
    """Accept either a path to an output file or the text of one.

    `Path(src).exists()` is not a safe test on arbitrary text: a filename component may
    not exceed NAME_MAX (255 bytes on Linux), so asking whether a 170 kB ORCA output
    "exists" raises OSError 36 rather than returning False. That is what used to take a
    whole run down when a job failed - error_message() resolved the path to text and then
    handed that text back to terminated_normally(), which tried to open it as a file.
    A real path is short and has no newline, so check those before touching the
    filesystem, and treat any OSError as "this was text, not a path".
    """
    if isinstance(src, Path):
        try:
            return src.read_text(errors="replace") if src.is_file() else str(src)
        except OSError:
            return str(src)
    if isinstance(src, str) and len(src) <= 4096 and "\n" not in src and "\x00" not in src:
        try:
            q = Path(src)
            if q.is_file():
                return q.read_text(errors="replace")
        except OSError:
            pass
    return str(src)


def terminated_normally(src) -> bool:
    return "****ORCA TERMINATED NORMALLY****" in _text(src)


def final_energy(src) -> float | None:
    hits = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", _text(src))
    return float(hits[-1]) if hits else None


def orbital_energies(src) -> list[dict]:
    """Return the LAST 'ORBITAL ENERGIES' table as [{no, occ, eh, ev}, ...]."""
    txt = _text(src)
    blocks = txt.split("ORBITAL ENERGIES")
    if len(blocks) < 2:
        return []
    rows: list[dict] = []
    for line in blocks[-1].splitlines():
        m = re.match(r"\s*(\d+)\s+(\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$", line)
        if m:
            rows.append({"no": int(m.group(1)), "occ": float(m.group(2)),
                         "eh": float(m.group(3)), "ev": float(m.group(4))})
        elif rows and line.strip() and not line.strip()[0].isdigit():
            break
    return rows


def homo_lumo(src) -> dict:
    """Frontier-orbital energies in eV, identified from the occupation column."""
    rows = orbital_energies(src)
    if not rows:
        return {}
    occ = [r for r in rows if r["occ"] > 0.5]
    vir = [r for r in rows if r["occ"] <= 0.5]
    if not occ or not vir:
        return {}
    homo, lumo = occ[-1], vir[0]
    out = {
        "homo_index": homo["no"], "lumo_index": lumo["no"],
        "e_homo_eV": homo["ev"], "e_lumo_eV": lumo["ev"],
        "e_homo_Eh": homo["eh"], "e_lumo_Eh": lumo["eh"],
        "gap_eV": lumo["ev"] - homo["ev"],
    }
    if len(occ) >= 2:
        out["e_homo_minus1_eV"] = occ[-2]["ev"]
        out["homo_degeneracy_eV"] = abs(occ[-1]["ev"] - occ[-2]["ev"])
    if len(vir) >= 2:
        out["e_lumo_plus1_eV"] = vir[1]["ev"]
        out["lumo_degeneracy_eV"] = abs(vir[1]["ev"] - vir[0]["ev"])
    return out


def frequencies(src) -> list[float]:
    """All vibrational frequencies in cm^-1 from the last FREQ block (zeros included)."""
    txt = _text(src)
    if "VIBRATIONAL FREQUENCIES" not in txt:
        return []
    block = txt.split("VIBRATIONAL FREQUENCIES")[-1]
    block = block.split("NORMAL MODES")[0]
    return [float(m) for m in re.findall(r"^\s*\d+:\s+(-?\d+\.\d+)\s*cm\*\*-1", block, re.M)]


def imaginary_modes(src, cutoff: float = -10.0) -> list[float]:
    """Frequencies more negative than `cutoff` - i.e. genuine imaginary modes."""
    return [f for f in frequencies(src) if f < cutoff]


_POP_HEADERS = {
    "Hirshfeld": "HIRSHFELD ANALYSIS",
    "Mulliken": "MULLIKEN ATOMIC CHARGES",
    "Loewdin": "LOEWDIN ATOMIC CHARGES",
}


def populations(src, scheme: str = "Hirshfeld") -> dict[int, float]:
    """Atomic charges as {0-based atom index: charge}."""
    txt = _text(src)
    header = _POP_HEADERS.get(scheme)
    if header is None or header not in txt:
        return {}
    block = txt.split(header)[-1]
    charges: dict[int, float] = {}
    if scheme == "Hirshfeld":
        # "   0 O   -0.319030    0.000000"
        for line in block.splitlines():
            m = re.match(r"\s*(\d+)\s+([A-Za-z]{1,2})\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$", line)
            if m:
                charges[int(m.group(1))] = float(m.group(3))
            elif charges and line.strip().startswith("-"):
                break
    else:
        # "   0 O :   -0.319030"
        for line in block.splitlines():
            m = re.match(r"\s*(\d+)\s+([A-Za-z]{1,2})\s*:\s*(-?\d+\.\d+)", line)
            if m:
                charges[int(m.group(1))] = float(m.group(3))
            elif charges and ("Sum of atomic charges" in line or line.strip().startswith("--")):
                break
    return charges


def atom_symbols(src) -> list[str]:
    """Element symbols in input order, from the last Cartesian coordinate block."""
    txt = _text(src)
    key = "CARTESIAN COORDINATES (ANGSTROEM)"
    if key not in txt:
        return []
    block = txt.split(key)[-1]
    syms = []
    for line in block.splitlines()[2:]:
        m = re.match(r"\s*([A-Za-z]{1,2})\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", line)
        if m:
            syms.append(m.group(1))
        elif syms:
            break
    return syms


def geometry(src) -> tuple[list[str], list[tuple[float, float, float]]]:
    """Last Cartesian geometry in Angstrom."""
    txt = _text(src)
    key = "CARTESIAN COORDINATES (ANGSTROEM)"
    if key not in txt:
        return [], []
    block = txt.split(key)[-1]
    syms, xyz = [], []
    for line in block.splitlines()[2:]:
        m = re.match(r"\s*([A-Za-z]{1,2})\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", line)
        if m:
            syms.append(m.group(1))
            xyz.append((float(m.group(2)), float(m.group(3)), float(m.group(4))))
        elif syms:
            break
    return syms, xyz


def dipole(src) -> float | None:
    hits = re.findall(r"Magnitude \(Debye\)\s*:\s*(-?\d+\.\d+)", _text(src))
    return float(hits[-1]) if hits else None


def thermochemistry(src) -> dict:
    txt = _text(src)
    out: dict[str, float] = {}
    pat = {
        "zpe_Eh": r"Zero point energy\s*\.\.\.\s*(-?\d+\.\d+)",
        "enthalpy_Eh": r"Total enthalpy\s*\.\.\.\s*(-?\d+\.\d+)",
        "gibbs_Eh": r"Final Gibbs free energy\s*\.\.\.\s*(-?\d+\.\d+)",
        "entropy_term_Eh": r"Total entropy correction\s*\.\.\.\s*(-?\d+\.\d+)",
        "temperature_K": r"Temperature\s*\.\.\.\s*(-?\d+\.\d+)\s*K",
    }
    for k, rx in pat.items():
        m = re.findall(rx, txt)
        if m:
            out[k] = float(m[-1])
    return out


def wall_time_s(src) -> float | None:
    m = re.findall(r"TOTAL RUN TIME:\s*(\d+) days (\d+) hours (\d+) minutes (\d+) seconds (\d+) msec",
                   _text(src))
    if not m:
        return None
    d, h, mi, s, ms = (int(x) for x in m[-1])
    return d * 86400 + h * 3600 + mi * 60 + s + ms / 1000.0


def convergence_failed(src) -> bool:
    txt = _text(src)
    return ("The optimization did not converge" in txt
            or "SCF NOT CONVERGED" in txt
            or "SERIOUS PROBLEM IN SOSCF" in txt)


def error_message(src) -> str:
    """A short diagnosis for a failed job."""
    txt = _text(src)
    if terminated_normally(txt):
        return ""
    for pat, msg in [
        (r"topology discovery failed", "Open MPI could not read the CPU topology - run serially "
                                      "(execution.nprocs = 1) or fix the MPI installation"),
        (r"not enough slots available",
         "execution.nprocs asks for more MPI slots than Open MPI will allocate. A slot is a "
         "PHYSICAL core, not a hardware thread, so a 6-core/12-thread laptop offers 6. "
         "Lower execution.nprocs to the physical core count"),
        (r"There are not enough slots",
         "execution.nprocs exceeds the physical core count - lower it"),
        (r"Error \(ORCA_SCF\): SERIOUS PROBLEM", "SCF failed to converge"),
        (r"The optimization did not converge", "geometry optimization did not converge"),
        (r"Unknown identifier", "an unrecognized keyword in the input file"),
        (r"not enough memory|insufficient memory|bad_alloc",
         "ran out of memory - lower execution.maxcore_mb"),
        (r"ORCA finished by error termination in (\w+)", "error termination"),
    ]:
        m = re.search(pat, txt)
        if m:
            return msg
    tail = [l for l in txt.splitlines() if l.strip()][-3:]
    return " / ".join(t.strip() for t in tail) if tail else "no output produced"
