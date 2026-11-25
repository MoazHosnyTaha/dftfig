"""Putting a congeneric series into one shared viewpoint.

`render.orientation: principal_plane` orients each molecule on its own inertia tensor,
which is right for flat fused rings but drifts as soon as the substituents differ: the
tensor is a property of the whole molecule, so adding a chain tilts the scaffold.
`first_molecule` instead fits every molecule onto the first one through the substructure
they have in common, so the scaffold sits in the same place in every panel and the
reader's eye compares the chemistry rather than the viewpoint.

The fit has to be done with the camera rather than by moving the atoms. The orbital and
electrostatic-potential panels are isosurfaces through a cube, and a cube's values live
on a grid fixed in space - rotating the molecule would leave its own density behind.
Rotating the camera moves everything in the scene together.

So each molecule gets a 3x3 rotation, written into the PyMOL view matrix:

  R_ref  the reference frame, from the principal axes of the common core of the first
         molecule, with a fixed sign convention so the result does not depend on which
         eigenvector LAPACK happened to return
  U_k    the rotation that superimposes molecule k's core onto the first molecule's,
         by least squares over the atom correspondence the substructure match gives
  R_k    R_ref . U_k, which is what the camera uses

Panel one therefore shows the core face-on, and every other panel shows its own core in
that same frame.
"""
from __future__ import annotations

import numpy as np

from .theme import note, warn

# Atoms whose correspondence is meaningful for a fit. Hydrogens are matched by the MCS
# too, but they are the atoms most likely to be placed differently by the optimizer, and
# they carry little of the scaffold's geometry.
MIN_CORE_ATOMS = 3


def common_core(mols: dict) -> dict[str, tuple[int, ...]]:
    """Atom indices of the maximum common substructure, per molecule.

    The tuples are returned in the *pattern's* atom order, so element i of one molecule's
    tuple corresponds to element i of every other's - which is the correspondence the
    least-squares fit needs. Empty dict if the molecules share nothing usable.
    """
    from rdkit import Chem
    from rdkit.Chem import rdFMCS

    names = list(mols)
    if len(names) < 2:
        return {}
    try:
        res = rdFMCS.FindMCS([mols[n] for n in names], timeout=20,
                             ringMatchesRingOnly=True, completeRingsOnly=True)
    except Exception as exc:                     # pragma: no cover - RDKit internals
        warn(f"maximum common substructure search failed ({exc})")
        return {}
    if res.canceled or not res.smartsString:
        warn("maximum common substructure search timed out - falling back to "
             "per-molecule principal axes")
        return {}

    patt = Chem.MolFromSmarts(res.smartsString)
    cores = {}
    for name in names:
        match = mols[name].GetSubstructMatch(patt)
        if not match:
            warn(f"{name}: does not contain the common core - falling back to its own "
                 "principal axes")
            return {}
        cores[name] = tuple(match)
    if res.numAtoms < MIN_CORE_ATOMS:
        warn(f"the common core is only {res.numAtoms} atom(s), too few to define an "
             "orientation - falling back to per-molecule principal axes")
        return {}
    note(f"common core: {res.numAtoms} atoms, {res.numBonds} bonds")
    return cores


def principal_frame(coords: np.ndarray) -> np.ndarray:
    """Rotation whose rows are the principal axes of `coords`, largest spread first.

    The sign of an eigenvector is arbitrary, and so is the handedness of the frame built
    from three of them. Both are pinned here - the largest component of each axis is made
    positive and the third axis is replaced by the cross product of the first two - so
    that two runs on the same coordinates, or on two molecules with the same core, agree.
    """
    x = np.asarray(coords, dtype=float)
    x = x - x.mean(axis=0)
    # Eigenvectors of the covariance matrix, descending by spread. The covariance shares
    # its eigenvectors with the inertia tensor, which is what PyMOL's orient uses.
    _, vecs = np.linalg.eigh(x.T @ x)
    axes = vecs[:, ::-1].T                       # rows, largest variance first

    for i, axis in enumerate(axes):
        if axis[np.argmax(np.abs(axis))] < 0:
            axes[i] = -axis
    axes[2] = np.cross(axes[0], axes[1])         # right-handed by construction
    return axes


def kabsch(mobile: np.ndarray, target: np.ndarray) -> np.ndarray:
    """Rotation U minimizing |U.mobile - target|, both centred first (Kabsch 1976)."""
    p = np.asarray(mobile, dtype=float)
    q = np.asarray(target, dtype=float)
    p = p - p.mean(axis=0)
    q = q - q.mean(axis=0)
    v, _, wt = np.linalg.svd(p.T @ q)
    # Guard against the reflection that SVD will happily return: a mirror image is a
    # better numerical fit than a rotation whenever the points are near-coplanar, and a
    # mirrored molecule in a figure is a drawing of a different stereoisomer.
    d = np.sign(np.linalg.det(v @ wt))
    return (wt.T @ np.diag([1.0, 1.0, d]) @ v.T)


def camera_rotations(names: list[str], cores: dict[str, tuple[int, ...]],
                     coords: dict[str, np.ndarray]) -> dict[str, list[float]]:
    """{name: nine floats} for the PyMOL view matrix, aligned onto the first molecule.

    Returns an empty dict when there is nothing to align to, which leaves the caller on
    its existing per-molecule orientation.
    """
    usable = [n for n in names if n in cores and n in coords]
    if len(usable) < 2:
        return {}

    ref = usable[0]
    ref_core = np.asarray(coords[ref])[list(cores[ref])]
    r_ref = principal_frame(ref_core)

    out = {ref: [float(v) for v in r_ref.ravel()]}
    for name in usable[1:]:
        core = np.asarray(coords[name])[list(cores[name])]
        if len(core) != len(ref_core):           # pragma: no cover - MCS guarantees this
            continue
        out[name] = [float(v) for v in (r_ref @ kabsch(core, ref_core)).ravel()]
    return out
