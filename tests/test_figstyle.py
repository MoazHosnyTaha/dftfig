

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
