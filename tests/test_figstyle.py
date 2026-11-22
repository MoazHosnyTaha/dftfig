

def test_interactive_backends_are_recognised():
    """A Qt or Tk backend blocks on QApplication creation in a headless batch run."""
    from dftfig.figstyle import needs_file_backend
    assert needs_file_backend("QtAgg")
    assert needs_file_backend("TkAgg")
    assert not needs_file_backend("Agg")
    assert not needs_file_backend("pdf")


def test_unresolved_backend_sentinel_is_left_alone():
    """Before matplotlib resolves it, rcParams['backend'] is a sentinel, not a name."""
    from matplotlib import rcsetup

    from dftfig.figstyle import needs_file_backend
    assert not needs_file_backend(rcsetup._auto_backend_sentinel)


def test_package_import_pins_a_file_backend():
    """Importing dftfig must pin the backend before anything resolves it.

    Asserted on the environment variable rather than on get_backend(), because
    get_backend() reports whatever the test machine can actually import - it says "agg"
    on a box with no Qt whether or not the fix is present, which would pass for the
    wrong reason.
    """
    import os
    import subprocess
    import sys
    env = {k: v for k, v in os.environ.items() if k != "MPLBACKEND"}
    out = subprocess.run(
        [sys.executable, "-c", "import os, dftfig; print(os.environ['MPLBACKEND'])"],
        capture_output=True, text=True, check=True, env=env)
    assert out.stdout.strip() == "Agg"


def test_an_empty_backend_variable_counts_as_unset():
    """MPLBACKEND='' makes matplotlib auto-select, so it must not be left alone."""
    import os
    import subprocess
    import sys
    out = subprocess.run(
        [sys.executable, "-c", "import os, dftfig; print(os.environ['MPLBACKEND'])"],
        capture_output=True, text=True, check=True,
        env={**os.environ, "MPLBACKEND": ""})
    assert out.stdout.strip() == "Agg"


def test_a_deliberate_backend_choice_is_respected():
    import os
    import subprocess
    import sys
    out = subprocess.run(
        [sys.executable, "-c", "import os, dftfig; print(os.environ['MPLBACKEND'])"],
        capture_output=True, text=True, check=True,
        env={**os.environ, "MPLBACKEND": "pdf"})
    assert out.stdout.strip() == "pdf"


def _fig(width, height=8.2):
    import matplotlib.pyplot as plt
    from dftfig.figstyle import apply_style
    apply_style(300)
    return plt.figure(figsize=(width, height))


FIELDS = ["B3LYP-D3BJ/def2-SVP opt // B3LYP-D3BJ/def2-TZVP property, CPCM",
          "orbital isosurface 0.02 a.u.",
          "density isosurface 0.0004 a.u.",
          "MEP range -0.037 to +0.062 a.u., shared"]


def test_a_narrow_figure_wraps_the_footnote():
    """Two molecules give a 4.9 in canvas; the method line does not fit on one line."""
    import matplotlib.pyplot as plt

    from dftfig.figstyle import wrap_fields
    fig = _fig(2.0 * 2 + 0.9)
    assert len(wrap_fields(fig, FIELDS, size=6.5)) > 1
    plt.close(fig)


def test_a_wide_figure_keeps_the_footnote_on_one_line():
    import matplotlib.pyplot as plt

    from dftfig.figstyle import wrap_fields
    fig = _fig(2.0 * 6 + 0.9)
    assert len(wrap_fields(fig, FIELDS, size=6.5)) == 1
    plt.close(fig)


def test_wrapping_never_drops_a_field():
    import matplotlib.pyplot as plt

    from dftfig.figstyle import wrap_fields
    fig = _fig(2.0 * 2 + 0.9)
    joined = "".join(wrap_fields(fig, FIELDS, size=6.5))
    for field in FIELDS:
        assert field in joined
    plt.close(fig)


def test_an_oversized_title_is_shrunk_to_fit():
    """With bbox_inches='tight' an over-wide title widens the file instead of clipping."""
    import matplotlib.pyplot as plt

    from dftfig.figstyle import shrink_to_width
    fig = _fig(4.9)
    t = fig.suptitle("Optimized geometry, frontier orbitals and electrostatic "
                     "potential in water (eps = 78.36)", fontsize=9.5, x=0.055,
                     ha="left")
    assert shrink_to_width(fig, t) < 9.5
    plt.close(fig)


def test_a_title_that_already_fits_is_left_alone():
    import matplotlib.pyplot as plt

    from dftfig.figstyle import shrink_to_width
    fig = _fig(8.9)
    t = fig.suptitle("Short title", fontsize=9.5, x=0.055, ha="left")
    assert shrink_to_width(fig, t) == 9.5
    plt.close(fig)
