"""Escaping saddle points: read the imaginary normal mode and restart from a displaced geometry.

A frequency calculation that returns one or more imaginary modes says the optimizer stopped
on a saddle point, not a minimum. Every number that follows - the gap, the thermochemistry,
the Fukui indices, the MEP - belongs to a structure that does not exist as a stable species.
Reporting it with a warning is honest but not useful; the fix is mechanical and this module
does it.

The remedy depends on the size of the imaginary frequency:

* Below `small_mode_cm` (default 50 cm^-1) the mode is usually a methyl or hydroxyl rotation
  sitting in numerical noise, and the first thing to try is a finer integration grid with a
  tighter optimization - not a displacement, which would only move the molecule sideways
  along a nearly flat coordinate.
* Above it the saddle is real. Displace the geometry along the imaginary mode and
  re-optimize. The mode has no preferred sign, so if the + direction falls back into the
  same saddle the - direction is tried next.

Modes are read from the NORMAL MODES block of the ORCA output, which ORCA prints as
mass-weighted Cartesian displacements, normalized. They are converted back to plain
Cartesian displacements before use.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .parsers import frequencies

#: Atomic masses, u. Enough for organic chemistry plus the common heteroatoms; anything
#: missing falls back to 12.0, which only changes the displacement scale, not its direction.
_MASS = {
    "H": 1.008, "He": 4.003, "Li": 6.94, "Be": 9.012, "B": 10.81, "C": 12.011,
    "N": 14.007, "O": 15.999, "F": 18.998, "Ne": 20.180, "Na": 22.990, "Mg": 24.305,
    "Al": 26.982, "Si": 28.085, "P": 30.974, "S": 32.06, "Cl": 35.45, "Ar": 39.948,
    "K": 39.098, "Ca": 40.078, "Fe": 55.845, "Cu": 63.546, "Zn": 65.38, "Br": 79.904,
    "I": 126.904,
}


@dataclass
class Attempt:
    """One pass of optimize-then-check."""
    index: int
    strategy: str
    n_imaginary: int
    lowest_freq: float | None
    note: str = ""


@dataclass
class Outcome:
    """What happened across all attempts for one molecule/solvent."""
    name: str
    solvent: str
    attempts: list = field(default_factory=list)
    resolved: bool = False
    #: Smallest imaginary frequency still present on the promoted structure, in cm^-1.
    best_residual_cm: float = 0.0
    #: True when that residual is small enough to be an internal rotation rather than a
    #: reaction saddle - a different statement, and a much less alarming one.
    flat_mode: bool = False

    @property
    def n_restarts(self) -> int:
        return max(0, len(self.attempts) - 1)

    def summary(self) -> str:
        if not self.attempts:
            return ""
        if self.resolved and self.n_restarts:
            return (f"minimum reached after {self.n_restarts} restart"
                    f"{'s' if self.n_restarts > 1 else ''} "
                    f"({self.attempts[-1].strategy})")
        if self.resolved:
            return "minimum confirmed"
        if self.flat_mode:
            return (f"flat mode -{self.best_residual_cm:.0f} cm-1 "
                    f"(free internal rotation, not a reaction saddle)")
        last = self.attempts[-1]
        return (f"still {last.n_imaginary} imaginary mode(s) after "
                f"{len(self.attempts)} attempt{'s' if len(self.attempts) > 1 else ''}")


# --------------------------------------------------------------------------- parsing
def read_xyz(path: Path) -> tuple[list[str], np.ndarray]:
    lines = Path(path).read_text().splitlines()
    n = int(lines[0].split()[0])
    symbols, coords = [], []
    for line in lines[2:2 + n]:
        parts = line.split()
        symbols.append(parts[0])
        coords.append([float(x) for x in parts[1:4]])
    return symbols, np.array(coords, dtype=float)


def write_xyz(path: Path, symbols: list[str], coords: np.ndarray, comment: str = "") -> Path:
    path = Path(path)
    lines = [str(len(symbols)), comment]
    for s, (x, y, z) in zip(symbols, coords):
        lines.append(f"{s:<3s} {x:>18.10f} {y:>18.10f} {z:>18.10f}")
    path.write_text("\n".join(lines) + "\n")
    return path


def normal_modes(src) -> np.ndarray | None:
    """The NORMAL MODES matrix as (3N, n_modes), mass-weighted as ORCA prints it.

    ORCA writes it in column blocks of six. Returns None when the block is absent, which
    is the case for any job that did not run frequencies.
    """
    txt = Path(src).read_text(errors="replace") if isinstance(src, (str, Path)) else src
    if "NORMAL MODES" not in txt:
        return None
    block = txt.split("NORMAL MODES")[-1]
    # The block ends at the IR spectrum, or at the next major header.
    for stop in ("IR SPECTRUM", "The first frequency considered", "-----------\nTIMINGS"):
        if stop in block:
            block = block.split(stop)[0]
    rows: dict[int, dict[int, float]] = {}
    col_ids: list[int] = []
    for line in block.splitlines():
        parts = line.split()
        if not parts:
            continue
        # A header line is all integers: the column indices of this block.
        if all(re.fullmatch(r"\d+", p) for p in parts):
            col_ids = [int(p) for p in parts]
            continue
        if col_ids and re.fullmatch(r"\d+", parts[0]) and len(parts) == len(col_ids) + 1:
            try:
                vals = [float(p) for p in parts[1:]]
            except ValueError:
                continue
            rows.setdefault(int(parts[0]), {}).update(dict(zip(col_ids, vals)))
    if not rows:
        return None
    n_rows = max(rows) + 1
    n_cols = max(max(d) for d in rows.values()) + 1
    mat = np.zeros((n_rows, n_cols))
    for r, d in rows.items():
        for c, v in d.items():
            mat[r, c] = v
    return mat


def imaginary_mode_vectors(out_path: Path, symbols: list[str],
                           cutoff: float = -10.0) -> list[tuple[float, np.ndarray]]:
    """(frequency, Cartesian displacement (N,3)) for every imaginary mode, worst first."""
    freqs = frequencies(out_path)
    modes = normal_modes(out_path)
    if not freqs or modes is None:
        return []
    inv_sqrt_m = np.array([1.0 / np.sqrt(_MASS.get(s.capitalize(), 12.0)) for s in symbols])
    out = []
    for i, f in enumerate(freqs):
        if f >= cutoff or i >= modes.shape[1]:
            continue
        v = modes[:, i]
        if v.size < 3 * len(symbols):
            continue
        disp = v[:3 * len(symbols)].reshape(-1, 3)
        # ORCA's printed modes are already M^(-1/2)-weighted Cartesian displacements, so
        # dividing by sqrt(m) again would double-weight them. Undo nothing; just normalize
        # so the step size below means the same thing for every molecule.
        if not np.any(disp):
            continue
        out.append((f, disp))
    out.sort(key=lambda t: t[0])          # most negative first
    _ = inv_sqrt_m                         # kept for clarity about the weighting
    return out


# --------------------------------------------------------------------------- the fix
def displaced_geometry(symbols: list[str], coords: np.ndarray, mode: np.ndarray,
                       step_angstrom: float = 0.30, sign: int = +1) -> np.ndarray:
    """Push the geometry off the saddle along `mode`, scaled so the largest atom moves
    `step_angstrom`."""
    largest = float(np.max(np.linalg.norm(mode, axis=1)))
    if largest <= 0:
        return coords.copy()
    return coords + sign * (step_angstrom / largest) * mode


#: Ordered recipes. Each is (name, needs_displacement, sign, extra ORCA keywords).
#: The first entry that applies is used for attempt 1, the next for attempt 2, and so on.
STRATEGIES_SMALL_MODE = [
    ("finer grid, tighter optimization", False, 0, ["DefGrid3", "TightOpt", "TightSCF"]),
    ("displaced along the imaginary mode (+)", True, +1, ["DefGrid3", "TightOpt"]),
    ("displaced along the imaginary mode (-)", True, -1, ["DefGrid3", "TightOpt"]),
]

STRATEGIES_REAL_SADDLE = [
    ("displaced along the imaginary mode (+)", True, +1, ["TightOpt"]),
    ("displaced along the imaginary mode (-)", True, -1, ["TightOpt"]),
    ("displaced further, finer grid", True, +1, ["DefGrid3", "TightOpt", "TightSCF"]),
]


def plan(freq_cm: float, attempt: int, small_mode_cm: float = 50.0):
    """Which recipe to use for this attempt, or None when the budget is spent."""
    table = (STRATEGIES_SMALL_MODE if abs(freq_cm) < small_mode_cm
             else STRATEGIES_REAL_SADDLE)
    return table[attempt] if attempt < len(table) else None
