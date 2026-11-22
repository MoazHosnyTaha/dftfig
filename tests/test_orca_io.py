from pathlib import Path



def _cfg(**over):
    from dftfig import params
    c = params.defaults()
    for dotted, v in over.items():
        sec, key = dotted.split("__")
        c[sec][key] = v
    return c


def test_property_step_defaults_to_the_optimization_functional():
    """An empty property.functional_sp must not change any existing run."""
    from dftfig.orca_io import keyword_line, stage_functional
    c = _cfg(theory__functional="B3LYP")
    assert stage_functional(c, "sp") == "B3LYP"
    assert keyword_line(c, "sp", "water").startswith("! B3LYP ")


def test_property_step_can_use_its_own_functional():
    from dftfig.orca_io import keyword_line
    c = _cfg(theory__functional="B3LYP", property__functional_sp="wB97X-D4")
    assert keyword_line(c, "opt", "water").startswith("! B3LYP ")
    assert keyword_line(c, "sp", "water").startswith("! wB97X-D4 ")


def test_the_ion_single_points_follow_the_property_functional():
    """Koopmans descriptors compare the neutral SP with the ions; they must match."""
    from dftfig.orca_io import stage_functional
    c = _cfg(theory__functional="B3LYP", property__functional_sp="PBE0 D3BJ")
    assert all(stage_functional(c, s) == "PBE0 D3BJ" for s in ("sp", "cat", "ani"))
    assert stage_functional(c, "opt") == "B3LYP"


def test_a_composite_property_functional_brings_its_own_basis():
    """r2SCAN-3c names no basis, whatever property.basis_sp says."""
    from dftfig.orca_io import keyword_line
    c = _cfg(theory__functional="B3LYP", property__functional_sp="r2SCAN-3c",
             property__basis_sp="def2-TZVP")
    sp = keyword_line(c, "sp", "water")
    assert sp.startswith("! r2SCAN-3c ")
    assert "def2-TZVP" not in sp
    assert "6-31+G(d,p)" in keyword_line(c, "opt", "water")


def test_dispersion_is_not_added_twice_on_the_property_step():
    """theory.dispersion must not land on a functional that already carries one."""
    from dftfig.orca_io import keyword_line
    c = _cfg(theory__functional="B3LYP", theory__dispersion="D3BJ",
             property__functional_sp="wB97X-D4")
    assert "D3BJ" not in keyword_line(c, "sp", "water")
    assert "D3BJ" in keyword_line(c, "opt", "water")


def test_the_figure_stamp_reports_both_halves():
    import pandas as pd

    from dftfig.figureapi import FigureContext
    c = _cfg(theory__functional="r2SCAN-3c", theory__dispersion="none",
             property__functional_sp="wB97X-D4")
    ctx = FigureContext(c, {}, pd.DataFrame(), pd.DataFrame(), Path("x"), {})
    stamp = ctx.level_of_theory()
    assert "r2SCAN-3c opt" in stamp          # composite: no basis printed
    assert "wB97X-D4/def2-TZVP property" in stamp
