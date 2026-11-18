"""Executable discovery.

The test that matters: on most Linux desktops /usr/bin/orca is the GNOME Orca SCREEN
READER. `shutil.which("orca")` finds it, ORCA-the-quantum-chemistry-program is never
reached, and the resulting failure blames the user's input file. A candidate is accepted
only when the helper binaries that ship with a real ORCA installation sit beside it.
"""
from __future__ import annotations

import stat

import pytest

from dftfig import discover


def _executable(path, body="#!/bin/sh\nexit 0\n"):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(body)
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return path


@pytest.fixture
def screen_reader(tmp_path):
    """What /usr/bin/orca actually is on a Linux desktop: one binary, no siblings."""
    return _executable(tmp_path / "usr_bin" / "orca")


@pytest.fixture
def real_orca(tmp_path):
    d = tmp_path / "orca_6_1_0"
    exe = _executable(d / "orca")
    for sibling in discover.ORCA_SIBLINGS:
        _executable(d / sibling)
    return exe


def test_a_bare_executable_named_orca_is_rejected(screen_reader):
    assert discover.is_real_orca(screen_reader) is False


def test_an_installation_with_orca_plot_and_orca_vpot_is_accepted(real_orca):
    assert discover.is_real_orca(real_orca) is True


def test_one_missing_sibling_is_enough_to_reject(real_orca):
    (real_orca.parent / "orca_vpot").unlink()
    assert discover.is_real_orca(real_orca) is False


def test_a_non_executable_file_is_rejected(tmp_path):
    d = tmp_path / "orca_6_1_0"
    d.mkdir()
    (d / "orca").write_text("not executable")
    for sibling in discover.ORCA_SIBLINGS:
        _executable(d / sibling)
    assert discover.is_real_orca(d / "orca") is False


@pytest.mark.parametrize("candidate", [None, "", "/nowhere/at/all/orca"])
def test_nothing_and_missing_paths_are_rejected(candidate):
    assert discover.is_real_orca(candidate) is False


def test_find_orca_prefers_a_verified_installation_over_the_screen_reader(
        monkeypatch, real_orca, screen_reader, tmp_path):
    monkeypatch.setenv("ORCA_PATH", str(real_orca))
    monkeypatch.setattr(discover.shutil, "which", lambda name: str(screen_reader))
    monkeypatch.setattr(discover, "_remembered", dict)
    monkeypatch.setattr(discover, "ORCA_HINTS", ())
    assert discover.find_orca() == real_orca


def test_find_orca_returns_none_when_only_the_screen_reader_exists(
        monkeypatch, screen_reader):
    monkeypatch.delenv("ORCA_PATH", raising=False)
    monkeypatch.delenv("ORCA_BIN", raising=False)
    monkeypatch.setattr(discover.shutil, "which", lambda name: str(screen_reader))
    monkeypatch.setattr(discover, "_remembered", dict)
    monkeypatch.setattr(discover, "ORCA_HINTS", ())
    assert discover.find_orca() is None


def test_a_directory_is_accepted_as_well_as_the_binary(real_orca):
    assert discover.find_orca(str(real_orca.parent)) == real_orca


def test_an_explicit_path_is_trusted_even_unverified(screen_reader):
    """The user who sets execution.orca_path has said what they mean.

    It is still reported by `advisories`, but overriding an explicit setting silently
    would be worse than running what was asked for.
    """
    assert discover.find_orca(str(screen_reader)) == screen_reader


def test_remember_and_read_back_the_user_level_path(tmp_path, monkeypatch, real_orca):
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "config"))
    discover.remember("orca_path", str(real_orca))
    assert discover.user_config_path().is_file()
    assert discover._remembered()["orca_path"] == str(real_orca)


def test_remember_never_raises_when_the_cache_cannot_be_written(monkeypatch, tmp_path):
    # A read-only HOME is not a reason to fail a 10-hour batch.
    blocked = tmp_path / "blocked"
    blocked.write_text("not a directory")
    monkeypatch.setenv("XDG_CONFIG_HOME", str(blocked))
    discover.remember("orca_path", "/somewhere")


def test_report_lists_what_this_machine_has():
    r = discover.report()
    assert set(r) == {"orca", "orca_candidates", "pymol", "user_config"}
    assert all(discover.is_real_orca(c) for c in r["orca_candidates"])


def test_find_pymol_prefers_the_interpreter_s_own_environment(monkeypatch, tmp_path):
    """A system-wide /usr/bin/pymol often starts a Python that cannot import pymol."""
    import sys

    envbin = tmp_path / "env" / "bin"
    _executable(envbin / "pymol")
    monkeypatch.setattr(sys, "executable", str(envbin / "python"))
    monkeypatch.setattr(discover.shutil, "which", lambda name: "/usr/bin/pymol")
    assert discover.find_pymol() == envbin / "pymol"


@pytest.mark.orca
def test_the_installed_orca_reports_its_version():
    exe = discover.find_orca()
    version = discover.orca_version(exe)
    assert version.startswith("6.") or version.startswith("5.")
