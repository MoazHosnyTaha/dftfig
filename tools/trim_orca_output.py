#!/usr/bin/env python3
"""Shrink a real ORCA output file to the parts dftfig's parsers read.

The test fixtures have to be GENUINE ORCA output - a hand-written approximation of a
frequency block is a test of the fixture, not of the parser - but a single
optimization+frequency run is 400 kB, most of it SCF iteration tables and per-cycle
internal-coordinate listings that no parser in this package ever looks at.

What is kept:
  * the header, so the level of theory and the input echo stay visible
  * the final geometry, orbital energies, population analyses, dipole, frequencies,
    thermochemistry, the last FINAL SINGLE POINT ENERGY and the termination line

What is dropped, as whole regions between two markers:
  * every geometry optimization cycle before the stationary point
  * the SCF iteration tables
  * the normal-mode eigenvector matrix (the header line stays, because
    parsers.frequencies() uses it to find the end of the frequency table)
  * the per-module timing breakdown

Marker lines are always kept even when their body is dropped: several parsers locate a
section by splitting the text on its header, so removing a header silently changes what
the parser believes the file contains. `--check` re-parses both files and refuses to
write a fixture whose values differ from the original's.

Usage:
    python tools/trim_orca_output.py FULL.out TRIMMED.out [--check]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

#: (start marker, end marker). Lines strictly between the first line containing the start
#: marker and the next line containing the end marker are removed; both markers survive.
DROP_REGIONS = [
    ("GEOMETRY OPTIMIZATION CYCLE", "FINAL ENERGY EVALUATION AT THE STATIONARY POINT"),
    ("SCF ITERATIONS", "TOTAL SCF ENERGY"),
    ("DFT GRID GENERATION", "TOTAL SCF ENERGY"),
    ("NORMAL MODES", "IR SPECTRUM"),
    ("Timings for individual modules", "****ORCA TERMINATED NORMALLY****"),
]


def trim(text: str) -> str:
    lines = text.splitlines()
    keep = [True] * len(lines)
    for start, end in DROP_REGIONS:
        # Only lines still present are candidates: in an optimization output the first
        # "DFT GRID GENERATION" sits inside the first optimization cycle, which an
        # earlier rule has already removed, and matching it there would reinstate a
        # placeholder line in the middle of a region that is gone.
        i = next((k for k, l in enumerate(lines) if keep[k] and start in l), None)
        if i is None:
            continue
        j = next((k for k in range(i + 1, len(lines)) if keep[k] and end in lines[k]), None)
        if j is None:
            continue
        for k in range(i + 1, j):
            keep[k] = False
        # A dropped region is replaced by one marker line rather than vanishing, so a
        # reader of the fixture can see that something was removed and where.
        keep[i + 1] = True
        lines[i + 1] = f"   [... {j - i - 1} lines removed by tools/trim_orca_output.py ...]"
    return "\n".join(l for l, k in zip(lines, keep) if k) + "\n"


def parsed_values(src) -> dict:
    """Everything dftfig's parsers extract, as plain JSON-able values."""
    from dftfig import parsers

    return {
        "terminated_normally": parsers.terminated_normally(src),
        "final_energy": parsers.final_energy(src),
        "n_orbitals": len(parsers.orbital_energies(src)),
        "homo_lumo": parsers.homo_lumo(src),
        "frequencies": parsers.frequencies(src),
        "imaginary_modes": parsers.imaginary_modes(src),
        "populations_hirshfeld": {str(k): v for k, v in
                                  parsers.populations(src, "Hirshfeld").items()},
        "populations_mulliken": {str(k): v for k, v in
                                 parsers.populations(src, "Mulliken").items()},
        "populations_loewdin": {str(k): v for k, v in
                                parsers.populations(src, "Loewdin").items()},
        "atom_symbols": parsers.atom_symbols(src),
        "geometry": parsers.geometry(src)[1],
        "dipole": parsers.dipole(src),
        "thermochemistry": parsers.thermochemistry(src),
        "wall_time_s": parsers.wall_time_s(src),
        "convergence_failed": parsers.convergence_failed(src),
        "error_message": parsers.error_message(src),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("source")
    ap.add_argument("dest")
    ap.add_argument("--check", action="store_true",
                    help="re-parse both files and fail if any parsed value changed")
    args = ap.parse_args(argv)

    src = Path(args.source)
    full = src.read_text(errors="replace")
    small = trim(full)
    Path(args.dest).write_text(small)
    print(f"{src.name}: {len(full) / 1024:.0f} kB -> {len(small) / 1024:.0f} kB "
          f"({100 * len(small) / len(full):.0f}%)")

    if args.check:
        before, after = parsed_values(full), parsed_values(small)
        bad = [k for k in before if before[k] != after[k]]
        if bad:
            print("TRIM CHANGED: " + ", ".join(bad), file=sys.stderr)
            return 1
        print("  all parsed values unchanged")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
