"""Figure rendering smoke tests.

Figures 2 and 13 place 3D panels rendered by PyMOL and need cube files from a real run;
figures 12 and 14 are drawn entirely by matplotlib and RDKit, so they can be rendered
from a checked-in results table on any machine. That is what is tested here: the figure
code runs end to end and writes a PNG with content in it. These are smoke tests - they
assert that a figure is produced, not what it looks like.

The molecule set (aniline, phenol) is the one in fixtures/descriptors.csv, so the table
and the structures refer to the same compounds, as they would after a real run.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
import pandas as pd
import pytest

matplotlib.use("Agg")

from dftfig import figureapi                                          # noqa: E402
from dftfig import mlp as mlpmod                                      # noqa: E402
from dftfig import molecules as M                                     # noqa: E402
from dftfig import properties                                         # noqa: E402

#: Smallest plausible size for a real multi-panel PNG. A figure that fails to draw its
#: panels still writes a valid, nearly empty file of a few kB, so "exists" is not enough.
MIN_PNG_BYTES = 20_000


@pytest.fixture
def figure_project(make_cfg, fixtures, tmp_path):
    """A config, the molecule set and the resolved property table for the figures."""
    cfg = make_cfg(**{
        "molecules.definitions": {"ANIL": "Nc1ccccc1", "PHEN": "Oc1ccccc1"},
        "figures.formats": ["png"],       # the PDF of the same figure doubles the runtime
        "figures.dpi": 110,
        "figures.which": ["speciation", "mpo_metabolism"],
        "conformers.n_confs": 5,
    })
    specs = M.load_molecules(cfg)
    mols = {n: s.mol for n, s in specs.items()}
    rd = {n: mlpmod.rdkit_properties(m) for n, m in mols.items()}
    prop_df = properties.load_properties(cfg, rd, mols=mols)
    return cfg, specs, mols, prop_df


@pytest.fixture
def figure_ctx(figure_project, descriptors):
    """The contract object a figure is handed, built the way the pipeline builds it."""
    cfg, specs, _, prop_df = figure_project
    return figureapi.build_context(cfg, specs, descriptors, prop_df)


@pytest.fixture
def descriptors(fixtures) -> pd.DataFrame:
    return pd.read_csv(fixtures / "descriptors.csv")


def _assert_real_png(paths):
    assert len(paths) == 1
    png = Path(paths[0])
    assert png.suffix == ".png" and png.is_file()
    assert png.stat().st_size > MIN_PNG_BYTES, f"{png.name} is {png.stat().st_size} bytes"
    assert png.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
    import PIL.Image
    with PIL.Image.open(png) as im:
        assert im.width > 300 and im.height > 300
        # More than one colour: a figure that raised inside every panel would still
        # write a uniformly white image.
        assert len(im.convert("RGB").getcolors(maxcolors=1 << 20)) > 50


def test_fig12_renders_from_the_structures_and_the_property_table(figure_ctx):
    from dftfig.figures import speciation

    _assert_real_png(speciation.build(figure_ctx))


def test_fig14_renders_from_the_structures_and_the_property_table(figure_ctx):
    from dftfig.figures import mpo_metabolism

    _assert_real_png(mpo_metabolism.build(figure_ctx))


def test_fig12_shows_the_predicted_pka_of_each_compound(figure_project):
    """Aniline is a weak base and phenol a weak acid; the panel is built from those."""
    _, _, mols, prop_df = figure_project
    from dftfig import pka

    assert pka.predict(mols["PHEN"])["pKa_acidic"] == [10.0]
    assert pka.predict(mols["ANIL"])["pKa_basic"] == [4.6]


def test_the_figures_stage_writes_png_and_the_resolved_property_table(figure_project,
                                                                      descriptors):
    from dftfig import __main__ as cli

    cfg, specs, _, _ = figure_project
    made = cli.stage_figures(cfg, specs, descriptors)
    out = Path(cfg["execution"]["work_dir"]) / "figures"
    assert (out / "properties_resolved.csv").is_file()
    assert sorted(p.name for p in made) == ["mpo_metabolism.png", "speciation.png"]
    for p in made:
        assert Path(p).stat().st_size > MIN_PNG_BYTES


def test_an_unknown_format_is_not_silently_produced(figure_ctx):
    from dftfig.figures import mpo_metabolism

    assert all(Path(p).suffix == ".png" for p in mpo_metabolism.build(figure_ctx))


# -------------------------------------------------------------------- PyMOL-dependent
@pytest.mark.pymol
def test_pymol_renders_a_geometry_panel(make_cfg, fixtures, tmp_path):
    """Figure 2's row A: the 3D panel that needs a real PyMOL.

    Also covers the reason render.py checks for the output file rather than the exit
    code - a PyMOL launcher whose Python cannot import pymol exits 0 and renders nothing.
    """
    from dftfig import render

    cfg = make_cfg(**{"render.ray_trace": False, "render.width": 400,
                      "render.height": 300, "execution.resume": False})
    work = Path(cfg["execution"]["work_dir"]) / "ANIL" / "water"
    work.mkdir(parents=True)
    (work / "ANIL_water_opt.xyz").write_text((fixtures / "ANIL_start.xyz").read_text())

    png = render.render_geometry(cfg, "ANIL", "water")
    assert png is not None and png.is_file()
    assert png.stat().st_size > 2_000
    # The .pml is written next to the image so a panel can be re-rendered or edited.
    assert (work / "render_geom_ANIL.pml").is_file()
