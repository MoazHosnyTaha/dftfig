"""The figure plugin contract: the context, discovery, and skipping on unmet needs.

These are the tests a third-party figure author depends on. If one of them fails, a
plugin written against the documented contract no longer works - which is what
API_VERSION exists to make visible.
"""
from __future__ import annotations

import importlib.metadata
import importlib.util
import shutil
import sys
from pathlib import Path

import matplotlib
import pandas as pd
import pytest

matplotlib.use("Agg")

from dftfig import figureapi                                          # noqa: E402
from dftfig import mlp as mlpmod                                      # noqa: E402
from dftfig import molecules as M                                     # noqa: E402
from dftfig import properties                                         # noqa: E402

EXAMPLES = Path(__file__).resolve().parent.parent / "examples"

#: A whole figure plugin, as a user would drop it into ./dft_figures.
TOY_PLUGIN = """
import matplotlib.pyplot as plt
from dftfig.figureapi import FigureContext

NAME = "toy"
TITLE = "A toy panel"
REQUIRES = ()

def build(ctx: FigureContext) -> list:
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, len(ctx.names)])
    return ctx.save(fig, "toy")
"""


@pytest.fixture
def project(make_cfg, fixtures):
    """A config plus the molecule set, descriptor and property tables of a real run."""
    cfg = make_cfg(**{
        "molecules.definitions": {"ANIL": "Nc1ccccc1", "PHEN": "Oc1ccccc1"},
        "figures.formats": ["png"],
        "figures.dpi": 110,
        "conformers.n_confs": 5,
    })
    specs = M.load_molecules(cfg)
    mols = {n: s.mol for n, s in specs.items()}
    rd = {n: mlpmod.rdkit_properties(m) for n, m in mols.items()}
    prop_df = properties.load_properties(cfg, rd, mols=mols)
    desc = pd.read_csv(fixtures / "descriptors.csv")
    return cfg, specs, desc, prop_df


@pytest.fixture
def ctx(project):
    cfg, specs, desc, prop_df = project
    return figureapi.build_context(cfg, specs, desc, prop_df)


@pytest.fixture
def plugin_dir(project, tmp_path):
    """An empty ./dft_figures for this project, already wired into the config."""
    cfg = project[0]
    d = tmp_path / "dft_figures"
    d.mkdir()
    cfg["figures"]["plugin_dir"] = str(d)
    return d


# ------------------------------------------------------------------ the data contract
def test_the_context_carries_everything_the_contract_promises(ctx):
    assert ctx.api_version == figureapi.API_VERSION
    assert ctx.names == ["ANIL", "PHEN"]
    assert set(ctx.mols) == set(ctx.names)

    c = ctx.compound("ANIL")
    assert c.name == "ANIL" and c.formula == "C6H7N"
    assert (c.charge, c.multiplicity) == (0, 1)
    assert c.mol.GetNumAtoms() == 14                      # hydrogens explicit
    assert c.run_dir.name == "ANIL"

    assert ctx.reference_solvent == "water"
    assert ctx.dielectric() > 70                          # water
    assert "opt //" in ctx.level_of_theory()
    assert ctx.output_dir.name == "figures"
    assert set(ctx.provenance) == {"dftfig_version", "config_sha"}
    assert ctx.param("figures.dpi") == 110
    assert ctx.param("figures.not_a_key", "fallback") == "fallback"


def test_descriptor_and_property_access(ctx):
    assert ctx.descriptor("ANIL", "gap_eV") == pytest.approx(5.335, abs=1e-3)
    assert ctx.descriptor("ANIL", "gap_eV", medium="dmso") is None    # not calculated
    assert ctx.descriptor("ANIL", "no_such_column", default=0) == 0
    assert ctx.descriptor_row("NOPE") is None

    assert ctx.property("ANIL", "TPSA") == pytest.approx(26.02, abs=0.01)
    assert ctx.property("ANIL", "no_such_column", default="-") == "-"
    assert ctx.property_row("ANIL")["molecule"] == "ANIL"


def test_rows_hold_python_scalars_not_numpy_ones(ctx):
    """Regression: `df.iloc[0]["is_minimum"]` is numpy.bool_, and `value is True` fails.

    the frontier-orbitals figure decides whether to print "minimum confirmed" with exactly that test, and the
    caption silently disappeared when the row came from .iloc instead of .iterrows.
    """
    row = ctx.descriptor_row("ANIL")
    assert row["is_minimum"] is True
    assert isinstance(row["molecule"], str)


def test_the_reference_solvent_is_the_default_medium(ctx):
    assert ctx.descriptor_row("ANIL").equals(ctx.descriptor_row("ANIL", "water"))


# ------------------------------------------------------------------------ file lookup
def test_cubes_are_located_by_quantity_not_by_file_name(ctx):
    assert ctx.cube("ANIL", "homo") is None               # nothing computed in this run
    work = ctx.run_dir("ANIL", "water")
    work.mkdir(parents=True, exist_ok=True)
    (work / "homo.cube").write_text("not really a cube")
    assert ctx.cube("ANIL", "homo").name == "homo.cube"


def test_an_unknown_cube_kind_is_an_error_not_a_silent_none(ctx):
    with pytest.raises(figureapi.FigureError, match="unknown cube kind"):
        ctx.cube("ANIL", "spin_density")


def test_an_unknown_image_kind_is_an_error(ctx):
    with pytest.raises(figureapi.FigureError, match="unknown image kind"):
        ctx.image("ANIL", "ribbon")


def test_fukui_and_som_tables_are_found_through_the_context(ctx):
    assert ctx.fukui("ANIL") is None
    assert ctx.som_table("ANIL") is None
    d = ctx.run_dir("ANIL")
    d.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"atom_index_0based": [0], "element": ["N"], "f_zero": [0.19]}).to_csv(
        d / "fukui_ANIL_water.csv", index=False)
    assert float(ctx.fukui("ANIL")["f_zero"].iloc[0]) == pytest.approx(0.19)


# ----------------------------------------------------------------------- requirements
def test_requirements_are_reported_as_sentences(ctx):
    ctx.config["render"]["engine"] = "none"
    unmet = ctx.unmet((figureapi.REQ_RENDERING,))
    assert len(unmet) == 1 and "rendering engine" in unmet[0]
    assert ctx.unmet((figureapi.REQ_PROPERTIES,)) == []


def test_an_empty_descriptor_table_is_an_unmet_requirement(project):
    cfg, specs, _, prop_df = project
    bare = figureapi.build_context(cfg, specs, None, prop_df)
    assert bare.unmet((figureapi.REQ_DESCRIPTORS,))


def test_a_figure_with_an_unmet_requirement_is_skipped_not_raised(project, plugin_dir,
                                                                  capsys):
    """The whole point of REQUIRES: a message, not a traceback from inside the renderer."""
    cfg, specs, desc, _ = project
    (plugin_dir / "needy.py").write_text(
        TOY_PLUGIN.replace("REQUIRES = ()", 'REQUIRES = ("rendering",)')
                  .replace('"toy"', '"needy"'))
    cfg["render"]["engine"] = "none"
    cfg["figures"]["which"] = ["needy"]

    from dftfig import __main__ as cli
    assert cli.stage_figures(cfg, specs, desc) == []
    assert "skipped" in capsys.readouterr().out


def test_a_figure_named_in_the_config_but_not_installed_is_reported(project, plugin_dir,
                                                                    capsys):
    cfg, specs, desc, _ = project
    cfg["figures"]["which"] = ["does_not_exist"]
    from dftfig import __main__ as cli
    assert cli.stage_figures(cfg, specs, desc) == []
    out = capsys.readouterr().out
    assert "does_not_exist" in out and "frontier_orbitals" in out          # lists what IS available


# --------------------------------------------------------------------------- discovery
def test_the_four_built_ins_are_discovered_with_their_declarations():
    found = figureapi.available()
    assert {"frontier_orbitals", "speciation", "lipophilicity", "mpo_metabolism"} <= set(found)
    assert found["frontier_orbitals"].requires == ("rendering", "cubes", "descriptors")
    assert found["lipophilicity"].requires == ("rendering",)
    assert all(found[n].builtin and found[n].title for n in ("frontier_orbitals", "speciation"))


def test_a_plain_py_file_in_the_project_directory_is_discovered(project, plugin_dir):
    cfg = project[0]
    (plugin_dir / "toy.py").write_text(TOY_PLUGIN)
    found = figureapi.available(cfg)
    assert "toy" in found
    assert found["toy"].title == "A toy panel"
    assert "toy.py" in found["toy"].origin
    assert not found["toy"].builtin


def test_a_project_file_overrides_a_built_in_of_the_same_name(project, plugin_dir):
    cfg = project[0]
    (plugin_dir / "mine.py").write_text(TOY_PLUGIN.replace('"toy"', '"frontier_orbitals"'))
    assert not figureapi.available(cfg)["frontier_orbitals"].builtin


def test_an_installed_package_is_discovered_through_the_entry_point_group(monkeypatch,
                                                                          plugin_dir):
    """What `pip install somebody-elses-figures` has to achieve, without installing one."""
    path = plugin_dir / "packaged.py"
    path.write_text(TOY_PLUGIN.replace('"toy"', '"packaged"'))
    loader = importlib.util.spec_from_file_location("thirdparty.figures.packaged", path)
    module = importlib.util.module_from_spec(loader)
    sys.modules["thirdparty.figures.packaged"] = module
    loader.loader.exec_module(module)

    class FakeEntryPoint:
        name = "packaged"
        value = "thirdparty.figures.packaged"

        def load(self):
            return module

    monkeypatch.setattr(
        importlib.metadata, "entry_points",
        lambda group=None: ([FakeEntryPoint()]
                            if group == figureapi.ENTRY_POINT_GROUP else []))
    found = figureapi.available()
    assert "packaged" in found
    assert found["packaged"].origin == "package 'thirdparty'"
    assert "frontier_orbitals" in found                        # the built-ins are still there


def test_a_broken_plugin_is_reported_and_skipped(project, plugin_dir, capsys):
    cfg = project[0]
    (plugin_dir / "broken.py").write_text("this is not python(")
    (plugin_dir / "toy.py").write_text(TOY_PLUGIN)
    found = figureapi.available(cfg)
    assert "toy" in found                       # one bad file does not stop the rest
    assert "broken.py" in capsys.readouterr().out


def test_a_module_without_build_is_rejected(project, plugin_dir):
    cfg = project[0]
    (plugin_dir / "nobuild.py").write_text("NAME = 'nobuild'\n")
    assert "nobuild" not in figureapi.available(cfg)


def test_a_plugin_asking_for_a_future_contract_is_refused(project, plugin_dir, capsys):
    cfg = project[0]
    (plugin_dir / "future.py").write_text(
        TOY_PLUGIN + f"\nAPI_VERSION = {figureapi.API_VERSION + 1}\n")
    assert "toy" not in figureapi.available(cfg)
    assert "figure API" in capsys.readouterr().out


def test_an_unknown_requirement_name_is_refused(project, plugin_dir, capsys):
    cfg = project[0]
    (plugin_dir / "weird.py").write_text(
        TOY_PLUGIN.replace("REQUIRES = ()", 'REQUIRES = ("a_pony",)'))
    assert "toy" not in figureapi.available(cfg)
    assert "unknown requirement" in capsys.readouterr().out


def test_load_names_what_is_available_when_asked_for_something_that_is_not():
    with pytest.raises(figureapi.FigureError, match="frontier_orbitals"):
        figureapi.load("nonexistent")


def test_a_discovered_plugin_is_a_legal_value_for_figures_which(project, plugin_dir):
    """Config validation has to accept plugin names, or they cannot be used at all."""
    from dftfig import config as C
    cfg = project[0]
    (plugin_dir / "toy.py").write_text(TOY_PLUGIN)
    cfg["figures"]["which"] = ["speciation", "toy"]
    C.validate(cfg)                                     # raises ConfigError on failure

    cfg["figures"]["which"] = ["speciation", "not_a_figure"]
    with pytest.raises(C.ConfigError):
        C.validate(cfg)


# ------------------------------------------------------------------- the example plugin
def test_the_shipped_example_plugin_is_discovered_and_renders(project, plugin_dir):
    """The tutorial's whole loop: drop the file in, see it listed, run it, get a PNG."""
    cfg, specs, desc, prop_df = project
    shutil.copy2(EXAMPLES / "fukui_bars" / "fukui_bars.py", plugin_dir / "fukui_bars.py")

    spec = figureapi.available(cfg)["fukui_bars"]
    assert spec.requires == ("fukui",)
    assert "fukui_bars.py" in spec.origin

    # The condensed Fukui table the figure reads, as the descriptors stage writes it.
    for name in ("ANIL", "PHEN"):
        d = Path(cfg["execution"]["work_dir"]) / name
        d.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({
            "atom_index_0based": range(6),
            "atom_index_1based": range(1, 7),
            "element": ["C"] * 5 + ["N" if name == "ANIL" else "O"],
            "f_plus": [0.05, 0.11, 0.07, 0.13, 0.04, 0.21],
            "f_minus": [0.09, 0.06, 0.12, 0.05, 0.10, 0.17],
        }).to_csv(d / f"fukui_{name}_water.csv", index=False)

    context = figureapi.build_context(cfg, specs, desc, prop_df)
    assert context.unmet(spec.requires) == []
    made = spec.build(context)
    assert [Path(p).name for p in made] == ["fukui_bars.png"]
    assert Path(made[0]).stat().st_size > 10_000


def test_the_example_plugin_draws_nothing_when_no_fukui_table_exists(project, plugin_dir):
    cfg, specs, desc, prop_df = project
    shutil.copy2(EXAMPLES / "fukui_bars" / "fukui_bars.py", plugin_dir / "fukui_bars.py")
    spec = figureapi.available(cfg)["fukui_bars"]
    assert spec.build(figureapi.build_context(cfg, specs, desc, prop_df)) == []
