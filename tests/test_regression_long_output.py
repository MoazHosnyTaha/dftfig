"""Regression: a failed job with a long output must not take the run down.

`parsers._text()` used to decide "is this a path or is this output text?" with
`Path(src).exists()`. A filename component may not exceed NAME_MAX (255 bytes on Linux),
so that call raises OSError errno 36 - "File name too long" - instead of returning False
when handed the text of an ORCA output. `error_message()` resolved its argument to text
and then passed that text to `terminated_normally()`, which called `_text()` on it again,
so EVERY failed job whose output exceeded 255 bytes - which is every failed job - ended
the whole batch with a traceback instead of a diagnosis. Fixed in 1.4.0.

Both halves of the trigger are covered: a very long single-component string (no path
separator, so the whole thing is one filename component) and a realistic multi-line
output of more than 100 kB.
"""
from __future__ import annotations

import errno

import pytest

from dftfig import parsers


@pytest.fixture
def long_output(fixtures) -> str:
    """A >100 kB string that looks exactly like what error_message() used to be given."""
    text = (fixtures / "fail_scf_not_converged.out").read_text()
    while len(text) < 100_000:
        text += text
    assert len(text) > 100_000
    return text


def test_terminated_normally_accepts_a_long_output_string(long_output):
    assert parsers.terminated_normally(long_output) is False


def test_error_message_accepts_a_long_output_string(long_output):
    assert parsers.error_message(long_output) == "SCF failed to converge"


def test_error_message_result_can_be_fed_back_through_the_parsers(long_output):
    """The exact call sequence that crashed: error_message -> _text -> terminated_normally."""
    assert parsers.terminated_normally(parsers.error_message(long_output)) is False


def test_a_single_component_longer_than_name_max_is_treated_as_text():
    """No separators and 4000 bytes: os.stat() on this raises ENAMETOOLONG, not False."""
    blob = "x" * 4000
    assert parsers.terminated_normally(blob) is False
    assert parsers.final_energy(blob) is None
    assert parsers.frequencies(blob) == []
    assert parsers.populations(blob) == {}


def test_the_underlying_os_call_really_does_raise(tmp_path):
    """Guards the assumption the fix rests on, so the test is not merely tautological."""
    import os

    with pytest.raises(OSError) as exc:
        os.stat(str(tmp_path / ("y" * 4000)))
    assert exc.value.errno in (errno.ENAMETOOLONG, errno.ENOENT)


def test_a_real_path_is_still_read_as_a_file(opt_out):
    assert parsers.terminated_normally(str(opt_out)) is True


def test_a_short_nonexistent_path_is_treated_as_text(tmp_path):
    missing = str(tmp_path / "not_here.out")
    assert parsers.terminated_normally(missing) is False
    assert parsers.error_message(missing).endswith("not_here.out")


@pytest.mark.parametrize("func", [parsers.terminated_normally, parsers.error_message,
                                  parsers.convergence_failed, parsers.final_energy,
                                  parsers.wall_time_s])
def test_every_entry_point_survives_a_long_string(func, long_output):
    func(long_output)
