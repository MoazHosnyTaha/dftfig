"""Shared fixtures.

The ORCA outputs under `fixtures/` are genuine ORCA 6.1 output for aniline in water
(B3LYP-D3BJ/def2-SVP opt+freq // def2-TZVP property, CPCM), trimmed with
`tools/trim_orca_output.py`; see fixtures/README.md. Nothing in the suite needs ORCA,
PyMOL or a network connection unless it carries the matching marker.
"""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

import pytest

# The package is imported from the checkout rather than from site-packages, so the suite
# tests the working tree even when an older dftfig is installed in the same environment.
# An explicit insert is required because pytest's rootdir is not added to sys.path when
# the interpreter runs with PYTHONSAFEPATH, as the conda launcher here does.
ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

FIXTURES = Path(__file__).parent / "fixtures"

#: The four real outputs, by the stage that produced them.
REAL_OUTPUTS = {
    "opt": "ANIL_water_opt.out",
    "sp": "ANIL_water_sp.out",
    "cat": "ANIL_water_cat.out",
    "ani": "ANIL_water_ani.out",
}


@pytest.fixture
def child_env():
    """Build the environment for a subprocess that has to `import dftfig`.

    A few things can only be tested honestly by starting a fresh interpreter: the
    matplotlib backend is pinned at import time, and the demo's whole claim is that it
    runs with no ORCA anywhere. Those children do not inherit the import path pytest
    arranges for the parent, so without this they pass or fail according to whether the
    package happens to be pip-installed in the current environment rather than according
    to the behaviour under test.

    Pass `KEY=None` to remove a variable from the child's environment.
    """
    import os

    def _env(**overrides):
        env = dict(os.environ)
        root = str(Path(__file__).resolve().parents[1])
        existing = env.get("PYTHONPATH")
        env["PYTHONPATH"] = os.pathsep.join([root, existing]) if existing else root
        for key, value in overrides.items():
            if value is None:
                env.pop(key, None)
            else:
                env[key] = value
        return env
    return _env


def pytest_runtest_setup(item):
    """Skip marked tests when the external program they drive is not installed.

    Doing this here rather than with a skipif decorator on each test keeps the
    availability check in one place, and means a contributor adding a PyMOL test only has
    to remember the marker.
    """
    from dftfig import discover

    if list(item.iter_markers(name="pymol")) and discover.find_pymol() is None:
        pytest.skip("PyMOL not installed")
    if list(item.iter_markers(name="orca")) and discover.find_orca() is None:
        pytest.skip("no verified ORCA installation on this machine")


@pytest.fixture(scope="session")
def fixtures() -> Path:
    return FIXTURES


@pytest.fixture(scope="session")
def full_values() -> dict:
    """What the parsers returned for the UNTRIMMED outputs, recorded before trimming."""
    return json.loads((FIXTURES / "full_file_values.json").read_text())


@pytest.fixture(params=sorted(REAL_OUTPUTS), scope="session")
def real_output(request) -> Path:
    return FIXTURES / REAL_OUTPUTS[request.param]


@pytest.fixture(scope="session")
def opt_out() -> Path:
    return FIXTURES / REAL_OUTPUTS["opt"]


@pytest.fixture(scope="session")
def sp_out() -> Path:
    return FIXTURES / REAL_OUTPUTS["sp"]


@pytest.fixture(scope="session")
def cat_out() -> Path:
    return FIXTURES / REAL_OUTPUTS["cat"]


@pytest.fixture
def make_cfg(tmp_path):
    """A registry-default config with dotted-key overrides, writing inside tmp_path.

    `config.load()` is deliberately NOT used: it searches the filesystem for ORCA and
    PyMOL and writes to the user's ~/.config, which would make these tests depend on the
    machine they run on.
    """
    from dftfig import params as P

    def _make(**overrides):
        cfg = copy.deepcopy(P.defaults())
        cfg["execution"]["work_dir"] = str(tmp_path / "runs")
        cfg["execution"]["scratch_dir"] = str(tmp_path / "scratch")
        cfg["execution"]["dry_run"] = True
        cfg["molecules"]["structure_source"] = "smiles"
        cfg["molecules"]["definitions"] = {"ANIL": "Nc1ccccc1"}
        cfg["figures"]["properties_csv"] = str(tmp_path / "properties.csv")
        for key, value in overrides.items():
            section, name = key.split(".", 1)
            cfg[section][name] = value
        return cfg

    return _make
