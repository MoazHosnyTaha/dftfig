"""The flagship figure's assembly, without PyMOL.

This module was the least-tested code in the package at 12% while being the figure most
people will look at, and both bugs found in it - a matplotlib backend that blocked
forever on import, and a canvas a third wider than its panels - were on lines no test
executed. The renders themselves need PyMOL and a cube; the layout around them does not,
and that is where the mistakes were. So the context here hands back no images, which
`show_image` draws as a placeholder, and everything else in `build` runs for real.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import pandas as pd
import pytest

from dftfig import params
from dftfig.figureapi import Compound, FigureContext
from dftfig.figures import frontier_orbitals


class _NoRenderContext(FigureContext):
    """A context whose 3D panels are unavailable, as on a machine with no PyMOL.

    `save` is overridden as well, for a duller reason: the real one closes the figure, so
    a test that wants to look at what was drawn has to hold on to it first.
    """

    def image(self, name, kind, medium=None, vrange=None):
        return None

    def surface_esp_stats(self, name, medium=None):
        return self.config.get("_fake_esp", {}).get(name)

    def save(self, fig, stem):
        self.captured.append(fig)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        path = self.output_dir / f"{stem}.png"
        fig.savefig(path, dpi=int(self.config["figures"]["dpi"]))
        return [path]


def _descriptors(rows: list[dict]) -> pd.DataFrame:
    base = {"medium": "water", "e_homo_eV": -5.6, "e_lumo_eV": -0.3, "gap_eV": 5.3,
            "n_imaginary": 0, "is_minimum": True, "saddle_restarts": 0,
            "flat_mode": False, "lowest_freq_cm-1": 410.0}
    return pd.DataFrame([{**base, **r} for r in rows])


def _context(tmp_path, names=("ANIL", "PHEN"), descriptors=None, **overrides):
    from rdkit import Chem

    cfg = params.defaults()
    cfg["execution"]["work_dir"] = str(tmp_path / "runs")
    cfg["figures"]["reference_solvent"] = "water"
    cfg["figures"]["formats"] = ["png"]
    cfg["solvation"]["solvents"] = ["water"]
    for key, value in overrides.items():
        section, _, name = key.partition("__")
        cfg[section][name] = value

    smiles = {"ANIL": "Nc1ccccc1", "PHEN": "Oc1ccccc1", "TOLU": "Cc1ccccc1",
              "XYLE": "Cc1ccccc1C"}
    compounds = {}
    for n in names:
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles.get(n, "c1ccccc1")))
        compounds[n] = Compound(name=n, mol=mol, charge=0, multiplicity=1,
                                formula="", run_dir=tmp_path / "runs" / n)
    if descriptors is None:
        descriptors = _descriptors([{"molecule": n} for n in names])
    out = tmp_path / "runs" / "figures"
    out.mkdir(parents=True, exist_ok=True)
    ctx = _NoRenderContext(cfg, compounds, descriptors, pd.DataFrame(), out,
                           {"dftfig_version": "0.0.0", "config_sha": "deadbeef"})
    ctx.captured = []
    return ctx


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


def _build(ctx):
    """Run the figure and hand back the matplotlib figure it drew."""
    paths = frontier_orbitals.build(ctx)
    assert paths and paths[0].exists(), "the figure was not written"
    assert len(ctx.captured) == 1
    return ctx.captured[0]


# ------------------------------------------------------------------- shape of the figure
def test_build_returns_one_figure_with_a_column_per_molecule(tmp_path):
    fig = _build(_context(tmp_path))
    # Four panel rows - geometry, LUMO, HOMO, MEP - plus the shared colour bar. Band B
    # of the printed figure is two rows, the orbital pair with the gap between them.
    assert len(fig.axes) == 4 * 2 + 1


@pytest.mark.parametrize("n", [1, 2, 3, 4])
def test_panel_count_follows_the_molecule_set(tmp_path, n):
    names = ("ANIL", "PHEN", "TOLU", "XYLE")[:n]
    fig = _build(_context(tmp_path, names=names))
    assert len(fig.axes) == 4 * n + 1


def test_the_canvas_is_sized_by_its_panels_not_its_footnote(tmp_path, monkeypatch):
    """The style saves with bbox='tight', which grows the output until it contains every
    artist - so a fixed-length footnote wider than the panels silently sets the width and
    pads the figure with empty space. Two molecules made it a third empty."""
    for n in (1, 2, 3, 4):
        names = ("ANIL", "PHEN", "TOLU", "XYLE")[:n]
        ctx = _context(tmp_path / f"w{n}", names=names)
        fig = _build(ctx)

        from PIL import Image
        # The file `_build` already wrote, at the configured save dpi - which is not
        # fig.get_dpi(), the on-screen one.
        with Image.open(ctx.output_dir / "frontier_orbitals.png") as im:
            saved_in = im.width / int(ctx.config["figures"]["dpi"])
        # Within a tenth of the canvas: the save padding plus whatever the right-anchored
        # "electron-poor" label puts past the last tick. Not exact, because chasing the
        # last tenth of an inch is not what this guards - the failure it guards against
        # was a canvas 73% wider than its panels at one molecule, and 42% empty.
        assert saved_in < fig.get_figwidth() * 1.10, \
            f"{n} molecule(s): saved {saved_in:.2f} in for a {fig.get_figwidth():.2f} in canvas"
        plt.close(fig)


def test_the_footnote_never_overlaps_the_colour_bar_label(tmp_path):
    """The footnote wraps, and the bottom margin has to be sized for the number of lines
    it actually took rather than for one."""
    from dftfig.figstyle import wrap_fields

    for n in (1, 4):
        names = ("ANIL", "PHEN", "TOLU", "XYLE")[:n]
        fig = _build(_context(tmp_path, names=names))
        texts = [t.get_text() for t in fig.texts]
        footnote_lines = [t for t in texts if "isosurface" in t]
        assert footnote_lines, "the method footnote is missing"
        # However many lines it wrapped to, the top of the block must sit below the axes.
        lines = wrap_fields(fig, footnote_lines[0].split("\n"), size=6.5)
        assert len(lines) >= 1
        plt.close(fig)


# -------------------------------------------------------------------- what it reports
def test_a_confirmed_minimum_is_stated(tmp_path):
    fig = _build(_context(tmp_path, names=("ANIL",)))
    assert any("minimum confirmed" in t.get_text() for t in fig.texts)


def test_an_imaginary_mode_is_flagged_as_not_a_minimum(tmp_path):
    desc = _descriptors([{"molecule": "ANIL", "n_imaginary": 2, "is_minimum": False,
                          "lowest_freq_cm-1": -31.0}])
    fig = _build(_context(tmp_path, names=("ANIL",), descriptors=desc))
    flags = [t.get_text() for t in fig.texts if "NOT A MINIMUM" in t.get_text()]
    assert flags and "2 imaginary mode(s)" in flags[0]


def test_a_flat_mode_is_not_called_a_saddle(tmp_path):
    """A residual low imaginary frequency from a free internal rotation is not a reaction
    saddle, and labelling it NOT A MINIMUM in red overstates what was found."""
    desc = _descriptors([{"molecule": "ANIL", "n_imaginary": 1, "is_minimum": False,
                          "flat_mode": True, "lowest_freq_cm-1": -27.8}])
    fig = _build(_context(tmp_path, names=("ANIL",), descriptors=desc))
    texts = [t.get_text() for t in fig.texts]
    assert any("flat mode" in t for t in texts)
    assert not any("NOT A MINIMUM" in t for t in texts)


def test_a_restarted_optimization_says_so(tmp_path):
    """A structure that needed a saddle escape is a valid minimum, but the published
    geometry is not the one the first optimization returned."""
    desc = _descriptors([{"molecule": "ANIL", "saddle_restarts": 2}])
    fig = _build(_context(tmp_path, names=("ANIL",), descriptors=desc))
    assert any("minimum after 2 restarts" in t.get_text() for t in fig.texts)


def test_orbital_energies_and_the_gap_are_printed(tmp_path):
    desc = _descriptors([{"molecule": "ANIL", "e_homo_eV": -5.644,
                          "e_lumo_eV": -0.309, "gap_eV": 5.335}])
    fig = _build(_context(tmp_path, names=("ANIL",), descriptors=desc))
    texts = " ".join(t.get_text() for t in fig.texts) + " ".join(
        t.get_text() for ax in fig.axes for t in ax.texts) + " ".join(
        ax.get_title() for ax in fig.axes)
    assert "-5.644" in texts and "-0.309" in texts and "5.335" in texts


def test_an_unbound_lumo_is_called_out(tmp_path):
    """A positive LUMO lies above the vacuum level, so the gap is not an excitation
    energy and the Koopmans electron affinity is meaningless."""
    desc = _descriptors([{"molecule": "ANIL", "e_lumo_eV": +0.84, "gap_eV": 6.5}])
    fig = _build(_context(tmp_path, names=("ANIL",), descriptors=desc))
    assert any("LUMO unbound" in t.get_text() for t in fig.texts) or \
        any("LUMO unbound" in t.get_text() for ax in fig.axes for t in ax.texts)


def test_a_missing_descriptor_row_does_not_stop_the_figure(tmp_path):
    """A compound whose property single point failed still gets its column, with the
    panels drawn and the energies left blank, rather than taking the figure down."""
    desc = _descriptors([{"molecule": "ANIL"}])          # PHEN absent
    fig = _build(_context(tmp_path, names=("ANIL", "PHEN"), descriptors=desc))
    assert len(fig.axes) == 4 * 2 + 1


def test_the_panel_titles_use_the_display_labels(tmp_path):
    """figures.labels renames a compound in the figure only - the run directory keeps the
    name the calculation was filed under."""
    ctx = _context(tmp_path, names=("ANIL", "PHEN"))
    ctx.config["figures"]["labels"] = {"ANIL": "Aniline"}
    fig = _build(ctx)
    titles = [ax.get_title() for ax in fig.axes]
    assert "Aniline" in titles
    assert "PHEN" in titles                               # unlabelled, so its own name


# ------------------------------------------------------------------ the colour scale
def test_the_shared_scale_spans_the_measured_range(tmp_path):
    ctx = _context(tmp_path)
    ctx.config["_fake_esp"] = {"ANIL": {"v_lo": -0.0371, "v_hi": 0.0402},
                               "PHEN": {"v_lo": -0.0288, "v_hi": 0.0617}}
    fig = _build(ctx)
    bar = fig.axes[-1]
    lo, hi = bar.get_xlim()
    assert lo <= -0.0371 and hi >= 0.0617
    # The label wraps to the canvas, so where the line break lands depends on how many
    # molecules there are; compare the words, not the layout.
    assert "shared by all panels" in " ".join(bar.get_xlabel().split())


def test_per_panel_scaling_says_the_panels_are_not_comparable(tmp_path):
    ctx = _context(tmp_path, render__esp_scale="per_panel")
    ctx.config["_fake_esp"] = {"ANIL": {"v_lo": -0.04, "v_hi": 0.04},
                               "PHEN": {"v_lo": -0.02, "v_hi": 0.09}}
    fig = _build(ctx)
    assert "NOT comparable" in " ".join(fig.axes[-1].get_xlabel().split())


def test_without_measured_surfaces_it_falls_back_to_the_configured_range(tmp_path):
    ctx = _context(tmp_path)                               # no _fake_esp
    fig = _build(ctx)
    r = float(ctx.param("render.esp_range"))
    lo, hi = fig.axes[-1].get_xlim()
    assert (lo, hi) == pytest.approx((-r, r), abs=1e-9)


# ------------------------------------------------------------------------ provenance
def test_the_version_and_settings_hash_are_stamped_on_the_figure(tmp_path):
    """What lets a printed figure be matched back to a row of descriptors.csv."""
    fig = _build(_context(tmp_path))
    stamp = " ".join(t.get_text() for t in fig.texts)
    assert "0.0.0" in stamp and "deadbeef" in stamp
