"""Steric and shape descriptors.

The geometries here are exact by construction - a regular tetrahedron, a planar
trigonal centre, a linear triatomic - so the expected values are known from the
geometry rather than copied from a previous run of this code.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from dftfig import sterics

morfeus = pytest.importorskip("morfeus", reason="steric descriptors need morfeus-ml")


# Methane: C at the origin, four H on alternating cube corners at 1.09 A.
_d = 1.09 / math.sqrt(3.0)
METHANE_SYMBOLS = ["C", "H", "H", "H", "H"]
METHANE_COORDS = np.array([
    [0.0, 0.0, 0.0],
    [_d, _d, _d], [_d, -_d, -_d], [-_d, _d, -_d], [-_d, -_d, _d],
])

# Formaldehyde-like planar trigonal carbon, all in z = 0.
PLANAR_SYMBOLS = ["C", "O", "H", "H"]
PLANAR_COORDS = np.array([
    [0.000, 0.000, 0.0],
    [0.000, 1.205, 0.0],
    [0.943, -0.588, 0.0],
    [-0.943, -0.588, 0.0],
])


def test_available_reports_a_reason_when_it_says_no():
    usable, why = sterics.available()
    assert usable is True and why == ""


def test_molecular_returns_every_declared_column():
    row = sterics.molecular(METHANE_SYMBOLS, METHANE_COORDS)
    assert list(row) == sterics.MOLECULAR_COLUMNS


def test_molecular_is_empty_without_atoms():
    assert sterics.molecular([], np.zeros((0, 3))) == {}
    assert sterics.atomic([], np.zeros((0, 3))) == []


def test_max_diameter_is_the_longest_interatomic_distance():
    # The four methane hydrogens sit on alternating cube corners, so the H...H
    # distance is the face diagonal of that cube: 2 * d * sqrt(2).
    expected = 2.0 * _d * math.sqrt(2.0)
    row = sterics.molecular(METHANE_SYMBOLS, METHANE_COORDS)
    assert row["max_diameter_A"] == pytest.approx(expected, abs=1e-6)


def test_methane_is_spherical_and_formaldehyde_is_not():
    """A tetrahedral centre has three equal principal moments: NPR1 = NPR2 = 1."""
    m = sterics.molecular(METHANE_SYMBOLS, METHANE_COORDS)
    assert m["npr1"] == pytest.approx(1.0, abs=1e-3)
    assert m["npr2"] == pytest.approx(1.0, abs=1e-3)
    assert m["asphericity"] == pytest.approx(0.0, abs=1e-3)


def test_a_flat_molecule_satisfies_the_planarity_relation():
    """For a planar molecule the perpendicular moment is the sum of the in-plane ones,
    I3 = I1 + I2, which is NPR1 + NPR2 = 1 once both are divided by I3."""
    p = sterics.molecular(PLANAR_SYMBOLS, PLANAR_COORDS)
    assert p["npr1"] + p["npr2"] == pytest.approx(1.0, abs=1e-6)
    assert p["npr1"] < p["npr2"]
    assert p["plane_best_fit_A"] == pytest.approx(0.0, abs=1e-6)


def test_atomic_covers_every_atom_in_input_order():
    rows = sterics.atomic(METHANE_SYMBOLS, METHANE_COORDS)
    assert len(rows) == len(METHANE_SYMBOLS)
    assert [r["element"] for r in rows] == METHANE_SYMBOLS
    assert [r["atom_index_0based"] for r in rows] == [0, 1, 2, 3, 4]
    assert [r["atom_index_1based"] for r in rows] == [1, 2, 3, 4, 5]
    for key in sterics.ATOMIC_COLUMNS:
        assert key in rows[0], key


def test_coordination_number_comes_from_bonds_not_proximity():
    rows = sterics.atomic(METHANE_SYMBOLS, METHANE_COORDS)
    assert rows[0]["coordination_number"] == 4           # carbon
    assert [r["coordination_number"] for r in rows[1:]] == [1, 1, 1, 1]


def test_pyramidalization_is_computed_only_for_three_coordinate_atoms():
    """MORFEUS picks neighbours by distance unless told otherwise, which hands every
    hydrogen three of them and a meaningless P. Only real three-coordinate atoms
    should carry a value."""
    rows = sterics.atomic(METHANE_SYMBOLS, METHANE_COORDS)
    # Four-coordinate carbon and one-coordinate hydrogens: nothing is pyramidalizable.
    assert all(math.isnan(r["pyramidalization_P"]) for r in rows)

    planar = sterics.atomic(PLANAR_SYMBOLS, PLANAR_COORDS)
    carbon = planar[0]
    assert carbon["coordination_number"] == 3
    # A flat trigonal centre is the definition of unpyramidalized.
    assert carbon["pyramidalization_P"] == pytest.approx(0.0, abs=1e-6)
    assert all(math.isnan(r["pyramidalization_P"]) for r in planar[1:])


def test_buried_and_free_volume_partition_the_probe_sphere():
    """Buried plus free is the whole 3.5 A sphere MORFEUS measures in, and the fraction
    is the first divided by that total."""
    sphere = 4.0 / 3.0 * math.pi * 3.5 ** 3
    for row in sterics.atomic(METHANE_SYMBOLS, METHANE_COORDS):
        total = row["buried_volume_A3"] + row["free_volume_A3"]
        assert total == pytest.approx(sphere, rel=1e-6)
        assert row["buried_volume_fraction"] == pytest.approx(
            row["buried_volume_A3"] / sphere, rel=1e-6)
        assert 0.0 < row["buried_volume_fraction"] < 1.0


def test_buried_volume_counts_hydrogens():
    """MORFEUS excludes hydrogens by default, which is right for a metal centre and
    wrong here: on an organic molecule they are most of what crowds an atom. With them
    excluded, methane's carbon would see nothing at all inside its sphere."""
    from morfeus import BuriedVolume

    without_hs = BuriedVolume(METHANE_SYMBOLS, METHANE_COORDS, 1).buried_volume
    with_hs = sterics.atomic(METHANE_SYMBOLS, METHANE_COORDS)[0]["buried_volume_A3"]
    assert without_hs == pytest.approx(0.0, abs=1e-9)
    assert with_hs > 1.0


def _run_dir(cfg, fixtures, solvents=("water",)):
    """Lay out a finished ANIL run from the stored ORCA outputs."""
    import shutil
    from pathlib import Path

    work = Path(cfg["execution"]["work_dir"])
    for solvent in solvents:
        d = work / "ANIL" / solvent
        d.mkdir(parents=True, exist_ok=True)
        for stage in ("sp", "cat", "ani", "opt"):
            src = fixtures / f"ANIL_water_{stage}.out"
            if src.exists():
                shutil.copy(src, d / f"ANIL_{solvent}_{stage}.out")
    return work


def test_atom_table_writes_one_row_per_atom(make_cfg, fixtures):
    from dftfig import descriptors

    cfg = make_cfg()
    work = _run_dir(cfg, fixtures)
    table = descriptors.atom_table(cfg)

    assert len(table) == 14                      # aniline, C6H5NH2
    assert (work / "atom_descriptors.csv").exists()
    assert list(table["element"][:2]) == ["N", "C"]
    for column in sterics.ATOMIC_COLUMNS:
        assert column in table.columns, column


def test_atom_table_carries_the_fukui_block_for_the_reference_medium(make_cfg, fixtures):
    """Fukui functions need the ion single points, which are run in the reference
    medium only - so those columns are filled there and absent elsewhere."""
    from dftfig import descriptors

    cfg = make_cfg()
    cfg["figures"]["reference_solvent"] = "water"
    _run_dir(cfg, fixtures)
    tables = descriptors.fukui(cfg)
    table = descriptors.atom_table(cfg, tables)

    assert "f_plus" in table.columns
    water = table[table["medium"] == "water"]
    assert water["f_plus"].notna().all()
    assert water["dual_descriptor"].notna().all()


def test_atom_table_reads_a_cached_fukui_file_when_it_is_not_handed_one(make_cfg, fixtures):
    """`dft descriptors` on a finished run should not need the Fukui tables passed in."""
    from dftfig import descriptors

    cfg = make_cfg()
    _run_dir(cfg, fixtures)
    descriptors.fukui(cfg)                        # leaves fukui_ANIL_water.csv behind
    assert "f_plus" in descriptors.atom_table(cfg).columns


def test_sterics_can_be_switched_off(make_cfg, fixtures):
    from pathlib import Path

    from dftfig import descriptors

    cfg = make_cfg()
    cfg["property"]["sterics"] = False
    work = _run_dir(cfg, fixtures)

    assert descriptors.atom_table(cfg).empty
    assert not (Path(work) / "atom_descriptors.csv").exists()
    row = descriptors.collect(cfg)
    assert not row.empty
    assert "sasa_area_A2" not in row.columns     # electronic half still produced


def test_collect_appends_the_molecular_steric_block(make_cfg, fixtures):
    from dftfig import descriptors

    cfg = make_cfg()
    _run_dir(cfg, fixtures)
    row = descriptors.collect(cfg)

    assert len(row) == 1
    for column in sterics.MOLECULAR_COLUMNS:
        assert column in row.columns, column
    assert row["sasa_area_A2"].iloc[0] > 0.0
    # The steric block must describe the same geometry as the electronic one, so the
    # molecular surface has to match what the per-atom table sums to.
    atoms = descriptors.atom_table(cfg)
    assert atoms["sasa_area_A2"].sum() == pytest.approx(
        row["sasa_area_A2"].iloc[0], rel=1e-6)


def test_atom_areas_sum_to_the_molecular_surface():
    rows = sterics.atomic(METHANE_SYMBOLS, METHANE_COORDS)
    molecule = sterics.molecular(METHANE_SYMBOLS, METHANE_COORDS)
    assert sum(r["sasa_area_A2"] for r in rows) == pytest.approx(
        molecule["sasa_area_A2"], rel=1e-6)
    assert sum(r["disp_area_A2"] for r in rows) == pytest.approx(
        molecule["disp_area_A2"], rel=1e-6)
