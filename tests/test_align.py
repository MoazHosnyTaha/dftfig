"""Putting a congeneric series in one viewpoint.

A wrong rotation here does not look wrong. Every panel still shows a plausible molecule
from a plausible angle, which is precisely the failure `render.orientation:
first_molecule` exists to prevent, so there is nothing in the picture to notice it by.
These tests therefore check the geometry numerically, and the PyMOL convention against
PyMOL itself.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from dftfig import align, render


def _rot(axis: str, degrees: float) -> np.ndarray:
    t = math.radians(degrees)
    c, s = math.cos(t), math.sin(t)
    if axis == "x":
        return np.array([[1, 0, 0], [0, c, -s], [0, s, c]], dtype=float)
    if axis == "y":
        return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]], dtype=float)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]], dtype=float)


# A flat hexagon in the xy-plane, which is what a shared aromatic scaffold amounts to.
RING = np.array([[math.cos(k * math.pi / 3), math.sin(k * math.pi / 3), 0.0]
                 for k in range(6)])


# ------------------------------------------------------------------------ the fit
def test_kabsch_recovers_a_known_rotation():
    """Rotate the ring by a known amount; the fit back onto it must be the inverse."""
    truth = _rot("y", 37.0) @ _rot("z", 12.0)
    found = align.kabsch((truth @ RING.T).T, RING)
    assert np.allclose(found, truth.T, atol=1e-9)


def test_kabsch_superimposes_to_machine_precision():
    rotated = (_rot("x", 61.0) @ RING.T).T + np.array([4.0, -2.0, 7.0])
    u = align.kabsch(rotated, RING)
    moved = (u @ (rotated - rotated.mean(axis=0)).T).T
    assert np.allclose(moved, RING - RING.mean(axis=0), atol=1e-9)


def test_kabsch_never_returns_a_reflection():
    """A mirror image fits coplanar points at least as well as a rotation does, and SVD
    will return one. A mirrored molecule in a figure is a drawing of the other
    enantiomer, so the reflection has to be ruled out rather than accepted."""
    mirrored = RING * np.array([1.0, 1.0, -1.0])
    for target in (RING, mirrored, (_rot("z", 25.0) @ RING.T).T):
        assert np.linalg.det(align.kabsch(RING, target)) == pytest.approx(1.0, abs=1e-9)


# ------------------------------------------------------------------ finding the core
def _mols(**smiles):
    from rdkit import Chem
    return {name: Chem.AddHs(Chem.MolFromSmiles(s)) for name, s in smiles.items()}


def test_common_core_finds_the_shared_scaffold():
    """Aniline and phenol share the ring and its hydrogens; the heteroatom is not part
    of the core, and the two matches are in the same correspondence order."""
    cores = align.common_core(_mols(ANIL="Nc1ccccc1", PHEN="Oc1ccccc1"))
    assert set(cores) == {"ANIL", "PHEN"}
    assert len(cores["ANIL"]) == len(cores["PHEN"]) >= 6


def test_common_core_spans_a_congeneric_series():
    cores = align.common_core(_mols(A="Nc1ccccc1", B="Nc1ccc(C)cc1", C="Nc1ccc(Cl)cc1"))
    assert len({len(v) for v in cores.values()}) == 1, "one correspondence per molecule"


def test_common_core_declines_when_there_is_nothing_shared():
    """Two molecules with no common substructure have no shared frame to align to, and
    saying so is what puts every panel back on its own principal axes."""
    assert align.common_core(_mols(A="O=C=O", B="N#N")) == {}


def test_common_core_declines_a_core_too_small_to_orient_by():
    """Water and ammonia have one atom in common - a hydrogen. A point defines no
    orientation, so there is nothing to align to."""
    assert align.common_core(_mols(A="O", B="N")) == {}


def test_a_methyl_group_is_a_large_enough_core():
    """The guard is about having an orientation at all, not about the core being a ring:
    four atoms of a methyl are enough to fit one molecule onto another."""
    cores = align.common_core(_mols(A="CO", B="CN"))
    assert len(cores) == 2 and len(cores["A"]) == len(cores["B"]) == 4


def test_common_core_needs_more_than_one_molecule():
    assert align.common_core(_mols(A="Nc1ccccc1")) == {}
    assert align.common_core({}) == {}


# ------------------------------------------------------------------- the reference
def test_principal_frame_is_orthonormal_and_right_handed():
    frame = align.principal_frame(RING + np.array([3.0, 3.0, 3.0]))
    assert np.allclose(frame @ frame.T, np.eye(3), atol=1e-9)
    assert np.linalg.det(frame) == pytest.approx(1.0, abs=1e-9)


def test_principal_frame_puts_a_flat_ring_face_on():
    """The axis of least spread becomes the view axis, so a planar core is seen flat."""
    tilted = (_rot("x", 40.0) @ RING.T).T
    seen = (align.principal_frame(tilted) @ tilted.T).T
    assert np.ptp(seen[:, 2]) == pytest.approx(0.0, abs=1e-9)


def test_principal_frame_is_reproducible_under_an_eigenvector_sign_flip():
    """Two molecules with the same core must get the same frame, and the sign an
    eigensolver returns is arbitrary - so the convention has to be pinned."""
    a = align.principal_frame(RING)
    b = align.principal_frame(RING[::-1] + 1e-12)      # same points, different order
    assert np.allclose(a, b, atol=1e-6)


# ------------------------------------------------------------- the composed rotations
def test_camera_rotations_bring_every_core_into_one_frame():
    cores = {"a": tuple(range(6)), "b": tuple(range(6))}
    coords = {"a": RING,
              "b": (_rot("y", 73.0) @ RING.T).T + np.array([10.0, 0.0, -4.0])}
    rots = align.camera_rotations(["a", "b"], cores, coords)
    assert set(rots) == {"a", "b"}

    seen = {}
    for name, flat in rots.items():
        core = coords[name][list(cores[name])]
        seen[name] = (np.array(flat).reshape(3, 3) @ (core - core.mean(axis=0)).T).T
    assert np.allclose(seen["a"], seen["b"], atol=1e-9)


def test_camera_rotations_need_at_least_two_molecules():
    assert align.camera_rotations(["a"], {"a": (0, 1, 2)}, {"a": RING}) == {}
    assert align.camera_rotations([], {}, {}) == {}


def test_the_first_molecule_defines_the_frame():
    """Panel one is the reference, so its rotation is the frame itself - swapping the
    order must change what everything is aligned to."""
    cores = {"a": tuple(range(6)), "b": tuple(range(6))}
    coords = {"a": RING, "b": (_rot("z", 50.0) @ RING.T).T}
    first = align.camera_rotations(["a", "b"], cores, coords)
    second = align.camera_rotations(["b", "a"], cores, coords)
    assert not np.allclose(first["a"], second["a"], atol=1e-6)


# -------------------------------------------------------------- the PyMOL convention
def test_the_rotation_is_transposed_for_pymol():
    """align returns the world-to-camera rotation with the view axes as rows; PyMOL's
    view matrix holds them as columns. Emitting it untransposed gives every panel a
    different, plausible-looking orientation."""
    cfg = {"render": {"orientation": "first_molecule"}}
    rotation = list(range(9))                     # 0..8, so the order is readable
    lines = render._orient(cfg, "mol", view_rotation=rotation)

    body = next(line for line in lines if line.startswith("_view[0:9]"))
    emitted = [float(v) for v in body.split("= [", 1)[1].rstrip("]").split(",")]
    assert emitted == [0, 3, 6, 1, 4, 7, 2, 5, 8]


def test_the_rotation_is_set_before_the_zoom():
    """zoom fits the molecule as it currently sits, so it has to come after the view is
    replaced or the framing is computed for the wrong orientation."""
    cfg = {"render": {"orientation": "first_molecule"}}
    lines = render._orient(cfg, "mol", view_rotation=[1, 0, 0, 0, 1, 0, 0, 0, 1])
    assert [i for i, line in enumerate(lines) if "set_view" in line][0] < \
           [i for i, line in enumerate(lines) if "zoom" in line][0]


def test_without_a_rotation_it_falls_back_to_per_molecule_axes():
    cfg = {"render": {"orientation": "first_molecule"}}
    assert any("cmd.orient" in line for line in render._orient(cfg, "mol", None))


@pytest.mark.pymol
def test_pymol_stores_the_view_axes_as_columns(tmp_path):
    """Pin the convention against PyMOL rather than against our own belief about it.

    The unit test above only checks that we transpose; this checks that transposing is
    the right thing to do. `cmd.orient` is PyMOL's own "show this face-on", so after it
    runs, the matrix read the correct way round must be the one that flattens a planar
    molecule along the view axis. That is the measurement the original bug failed.
    """
    import json
    import subprocess

    from dftfig import demo, discover

    xyz = demo.DATA / "runs" / "ANIL" / "water" / "ANIL_water_sp.xyz"
    script = tmp_path / "probe.pml"
    out = tmp_path / "probe.json"
    script.write_text(
        "import json\n"
        "import numpy as np\n"
        f"cmd.load({str(xyz)!r}, 'mol')\n"
        "cmd.orient('mol')\n"
        "R = np.array(cmd.get_view()[0:9]).reshape(3, 3)\n"
        "x = np.array(cmd.get_coords('mol'), dtype=float)\n"
        "x = x - x.mean(axis=0)\n"
        # One statement per line: a .pml is read line by line, so a call split across
        # lines is silently dropped rather than reported.
        "depth = {'as_stored': float(np.ptp((R @ x.T).T[:, 2])), "
        "'transposed': float(np.ptp((R.T @ x.T).T[:, 2]))}\n"
        f"json.dump(depth, open({str(out)!r}, 'w'))\n"
    )
    subprocess.run([discover.find_pymol(), "-cq", str(script)],
                   check=True, capture_output=True, timeout=300)

    depth = json.loads(out.read_text())
    # Aniline is planar apart from the amine hydrogens, so face-on is well under an
    # angstrom of depth while any other axis spans most of the ring.
    assert depth["transposed"] < 0.8
    assert depth["as_stored"] > 2.0


def test_orientation_none_does_not_rotate_at_all():
    lines = render._orient({"render": {"orientation": "none"}}, "mol", [1.0] * 9)
    assert not any("orient" in line or "set_view" in line for line in lines)
