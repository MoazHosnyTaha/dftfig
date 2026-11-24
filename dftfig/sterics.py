"""Steric and shape descriptors from an optimized geometry.

Everything here is computed from element symbols and Cartesian coordinates alone, so it
runs on the structure ORCA already optimized: no extra quantum chemistry, no external
program, no conformer ensemble. That is a deliberate boundary. The electronic half of a
descriptor set is produced by `descriptors.py` from the ORCA single points at the DFT
level the run was configured with; this module adds the steric half beside it.

Two sources, kept separate because they answer different questions:

  MORFEUS   solvent-accessible surface, Hirshfeld-partitioned dispersion (P_int) and
            buried volume - how much space an atom occupies and how exposed it is.
  RDKit     principal moments of inertia and the shape indices derived from them -
            how the mass of the whole molecule is distributed.

Atom indices follow MORFEUS's own 1-based convention internally and are reported in both
conventions in the output table, matching `descriptors.fukui`.

A note on pyramidalization: MORFEUS picks neighbours by distance when it is not told
which atoms are bonded, which gives every hydrogen three "neighbours" and a meaningless
P near 2. It is computed here only for atoms with exactly three bonded neighbours, from
the covalent-radius connectivity matrix, and left empty everywhere else.
"""
from __future__ import annotations

from collections.abc import Sequence

import numpy as np

# Column order for the molecular block, so the CSV is stable across runs.
MOLECULAR_COLUMNS = [
    "sasa_area_A2", "sasa_volume_A3",
    "disp_area_A2", "disp_volume_A3", "disp_p_int", "disp_p_min", "disp_p_max",
    "max_diameter_A",
    "pmi1", "pmi2", "pmi3", "npr1", "npr2",
    "asphericity", "eccentricity", "inertial_shape_factor",
    "radius_of_gyration_A", "spherocity_index", "plane_best_fit_A",
]

ATOMIC_COLUMNS = [
    "sasa_area_A2", "disp_area_A2", "disp_p_int", "disp_p_min", "disp_p_max",
    "buried_volume_A3", "free_volume_A3", "buried_volume_fraction",
    "coordination_number", "pyramidalization_P", "pyramidalization_angle_deg",
]


def available() -> tuple[bool, str]:
    """(usable, reason). Steric descriptors are optional: a run without MORFEUS still
    produces the electronic table, with the steric columns absent rather than empty."""
    try:
        import morfeus  # noqa: F401
    except ImportError as exc:
        return False, f"morfeus-ml is not installed ({exc})"
    return True, ""


def _xyz_block(symbols: Sequence[str], coords: np.ndarray) -> str:
    lines = [f"{len(symbols)}", ""]
    for s, (x, y, z) in zip(symbols, coords):
        lines.append(f"{s} {x:.8f} {y:.8f} {z:.8f}")
    return "\n".join(lines)


def _shape(symbols: Sequence[str], coords: np.ndarray) -> dict:
    """RDKit's inertial shape indices, evaluated on the DFT geometry.

    Built from an XYZ block rather than the input molecule: these descriptors need
    coordinates and atomic masses, not bonds, and going through XYZ keeps the numbers
    tied to the geometry that was actually optimized.
    """
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors as rd

    mol = Chem.MolFromXYZBlock(_xyz_block(symbols, coords))
    if mol is None or mol.GetNumConformers() == 0:
        return {}
    fns = {
        "pmi1": rd.CalcPMI1, "pmi2": rd.CalcPMI2, "pmi3": rd.CalcPMI3,
        "npr1": rd.CalcNPR1, "npr2": rd.CalcNPR2,
        "asphericity": rd.CalcAsphericity, "eccentricity": rd.CalcEccentricity,
        "inertial_shape_factor": rd.CalcInertialShapeFactor,
        "radius_of_gyration_A": rd.CalcRadiusOfGyration,
        "spherocity_index": rd.CalcSpherocityIndex,
        "plane_best_fit_A": rd.CalcPBF,
    }
    out = {}
    for key, fn in fns.items():
        try:
            out[key] = float(fn(mol))
        except Exception:        # a linear or single-atom case trips some of these
            out[key] = float("nan")
    return out


def molecular(symbols: Sequence[str], coords) -> dict:
    """Whole-molecule steric and shape descriptors. Empty dict if MORFEUS is missing."""
    usable, _ = available()
    if not usable or not symbols:
        return {}
    from morfeus import SASA, Dispersion

    c = np.asarray(coords, dtype=float)
    sasa = SASA(symbols, c)
    disp = Dispersion(symbols, c)
    row = {
        "sasa_area_A2": float(sasa.area),
        "sasa_volume_A3": float(sasa.volume),
        "disp_area_A2": float(disp.area),
        "disp_volume_A3": float(disp.volume),
        "disp_p_int": float(disp.p_int),
        "disp_p_min": float(disp.p_min),
        "disp_p_max": float(disp.p_max),
        "max_diameter_A": _max_diameter(c),
    }
    row.update(_shape(symbols, c))
    return {k: row.get(k, float("nan")) for k in MOLECULAR_COLUMNS}


def _max_diameter(coords: np.ndarray) -> float:
    """Largest interatomic distance - the molecule's longest dimension."""
    if len(coords) < 2:
        return 0.0
    d = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)
    return float(d.max())


def atomic(symbols: Sequence[str], coords) -> list[dict]:
    """Atom-centred steric descriptors, one dict per atom in input order."""
    usable, _ = available()
    if not usable or not symbols:
        return []
    from morfeus import SASA, BuriedVolume, Dispersion, Pyramidalization
    from morfeus.utils import get_connectivity_matrix

    c = np.asarray(coords, dtype=float)
    sasa = SASA(symbols, c)
    disp = Dispersion(symbols, c)
    conn = get_connectivity_matrix(c, elements=list(symbols))

    rows = []
    for i in range(1, len(symbols) + 1):        # MORFEUS indexes atoms from one
        neighbours = [int(j) + 1 for j in np.flatnonzero(conn[i - 1])]
        row = {
            "atom_index_0based": i - 1,
            "atom_index_1based": i,
            "element": symbols[i - 1],
            "sasa_area_A2": float(sasa.atom_areas.get(i, float("nan"))),
            "disp_area_A2": float(disp.atom_areas.get(i, float("nan"))),
            "disp_p_int": float(disp.atom_p_int.get(i, float("nan"))),
            "disp_p_min": float(disp.atom_p_min.get(i, float("nan"))),
            "disp_p_max": float(disp.atom_p_max.get(i, float("nan"))),
            "coordination_number": len(neighbours),
            "pyramidalization_P": float("nan"),
            "pyramidalization_angle_deg": float("nan"),
        }
        try:
            # include_hs is off by default in MORFEUS, which is right for a metal centre
            # in a phosphine cone but wrong here: on an organic molecule the hydrogens
            # are most of what crowds an atom.
            bv = BuriedVolume(symbols, c, i, include_hs=True)
            row["buried_volume_A3"] = float(bv.buried_volume)
            row["free_volume_A3"] = float(bv.free_volume)
            row["buried_volume_fraction"] = float(bv.fraction_buried_volume)
        except Exception:
            row["buried_volume_A3"] = float("nan")
            row["free_volume_A3"] = float("nan")
            row["buried_volume_fraction"] = float("nan")
        if len(neighbours) == 3:
            try:
                # A planar centre sums to 360 degrees exactly, and MORFEUS takes the
                # square root of (360 - sum) for its angle variant - which is a tiny
                # negative number in floating point. The value reported here is the
                # distance-based P, so the warning that produces is noise.
                with np.errstate(invalid="ignore"):
                    p = Pyramidalization(c, i, neighbor_indices=neighbours)
                row["pyramidalization_P"] = float(p.P)
                row["pyramidalization_angle_deg"] = float(p.alpha)
            except Exception:
                pass
        rows.append(row)
    return rows
