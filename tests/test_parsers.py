"""Parser coverage against real ORCA 6.1 output.

Every public function in `dftfig.parsers` is exercised here. The reference numbers are
the values the parsers returned for the UNTRIMMED outputs, recorded in
fixtures/full_file_values.json before the fixtures were cut down, so a trim that
changed what a parser sees fails the suite rather than quietly moving the baseline.
"""
from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path

import pytest

from dftfig import parsers

ANILINE_SYMBOLS = ["N", "C", "C", "C", "C", "C", "C", "H", "H", "H", "H", "H", "H", "H"]


def _trim_tool():
    """Load tools/trim_orca_output.py, which is a script rather than part of the package."""
    path = Path(__file__).resolve().parent.parent / "tools" / "trim_orca_output.py"
    spec = importlib.util.spec_from_file_location("trim_orca_output", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# --------------------------------------------------------------------------- trimming
def test_trimmed_fixture_parses_identically_to_the_full_output(real_output, full_values):
    """The fixtures are cut down; the values the parsers read from them are not.

    This is the test that licenses the trimming: fixtures/full_file_values.json holds what
    each parser returned for the original 70-400 kB file, so any future edit to a fixture
    that drops a section a parser depends on shows up here.
    """
    recorded = full_values[real_output.name]["values"]
    # Through JSON: the recorded values lost the tuple/list distinction on the way to
    # disk, and a coordinate triple is the same geometry either way.
    got = json.loads(json.dumps(_trim_tool().parsed_values(real_output.read_text())))
    differing = {k: (recorded[k], got[k]) for k in recorded if recorded[k] != got[k]}
    assert not differing, f"trimming changed {sorted(differing)}"


def test_fixtures_are_substantially_smaller_than_the_originals(real_output, full_values):
    full = full_values[real_output.name]["source_bytes"]
    assert real_output.stat().st_size <= full


# ------------------------------------------------------------------------ termination
def test_terminated_normally_on_every_successful_job(real_output):
    assert parsers.terminated_normally(real_output) is True


@pytest.mark.parametrize("name", ["fail_scf_not_converged.out", "fail_mpi_slots.out",
                                  "fail_opt_not_converged.out"])
def test_terminated_normally_is_false_for_failures(fixtures, name):
    assert parsers.terminated_normally(fixtures / name) is False


def test_text_accepts_a_path_a_path_string_and_the_output_itself(opt_out):
    assert parsers.final_energy(opt_out) == parsers.final_energy(str(opt_out))
    assert parsers.final_energy(opt_out) == parsers.final_energy(opt_out.read_text())


# ----------------------------------------------------------------------------- energy
def test_final_energy_is_the_last_one_in_the_file(opt_out):
    # The optimization prints one FINAL SINGLE POINT ENERGY per cycle; only the value at
    # the stationary point is the result.
    assert parsers.final_energy(opt_out) == pytest.approx(-287.250216764344, abs=1e-9)


def test_final_energy_of_each_ion(cat_out, sp_out):
    assert parsers.final_energy(sp_out) == pytest.approx(-287.566886455514, abs=1e-9)
    assert parsers.final_energy(cat_out) == pytest.approx(-287.366123957532, abs=1e-9)


def test_final_energy_is_none_when_nothing_was_computed(fixtures):
    assert parsers.final_energy(fixtures / "fail_mpi_slots.out") is None


# ------------------------------------------------------------------- orbital energies
def test_orbital_energies_table(sp_out):
    rows = parsers.orbital_energies(sp_out)
    assert len(rows) == 36
    assert [r["no"] for r in rows] == list(range(36))
    assert rows[0]["eh"] == pytest.approx(-14.314052)
    assert rows[0]["ev"] == pytest.approx(-389.5052)
    # The eV column is the Eh column converted, which is what makes it safe for
    # descriptors.py to use either.
    assert rows[0]["ev"] == pytest.approx(rows[0]["eh"] * parsers.HARTREE_TO_EV, abs=1e-3)


def test_homo_lumo_indices_follow_the_occupation_column(sp_out):
    fr = parsers.homo_lumo(sp_out)
    # Aniline is C6H7N: 50 electrons, 25 doubly occupied orbitals, so ORCA's 0-based
    # HOMO is number 24.
    assert fr["homo_index"] == 24
    assert fr["lumo_index"] == 25
    assert fr["e_homo_eV"] == pytest.approx(-5.644)
    assert fr["e_lumo_eV"] == pytest.approx(-0.309)
    assert fr["gap_eV"] == pytest.approx(fr["e_lumo_eV"] - fr["e_homo_eV"])
    assert fr["gap_eV"] == pytest.approx(5.335)
    assert fr["e_homo_minus1_eV"] < fr["e_homo_eV"] < fr["e_lumo_eV"] < fr["e_lumo_plus1_eV"]


def test_homo_lumo_of_an_open_shell_output_is_the_spin_up_set(cat_out):
    """A UKS job prints SPIN UP then SPIN DOWN under one ORBITAL ENERGIES header.

    The row scan stops at the "SPIN DOWN ORBITALS" line because it is neither a data row
    nor blank, so the frontier pair reported for an ion is the ALPHA one. That is the set
    the Koopmans descriptors want; this test pins the behaviour so it cannot change
    silently into a mix of the two spins.
    """
    rows = parsers.orbital_energies(cat_out)
    assert all(r["occ"] in (0.0, 1.0) for r in rows)
    fr = parsers.homo_lumo(cat_out)
    assert fr["homo_index"] == 24 and fr["lumo_index"] == 25
    assert fr["e_homo_eV"] == pytest.approx(-7.6457)


def test_homo_lumo_is_empty_without_an_orbital_table(fixtures):
    assert parsers.homo_lumo(fixtures / "fail_mpi_slots.out") == {}


# ------------------------------------------------------------------------ frequencies
def test_frequencies_include_the_six_zero_modes(opt_out):
    freqs = parsers.frequencies(opt_out)
    # 14 atoms -> 3N = 42 printed modes, the first six being translations and rotations.
    assert len(freqs) == 42
    assert freqs[:6] == [0.0] * 6
    assert freqs[6] == pytest.approx(221.95)
    assert freqs == sorted(freqs)


def test_a_converged_minimum_has_no_imaginary_modes(opt_out):
    assert parsers.imaginary_modes(opt_out) == []


def test_a_single_point_has_no_frequency_block(sp_out):
    assert parsers.frequencies(sp_out) == []


def test_imaginary_modes_respect_the_cutoff(fixtures):
    """The cutoff is what separates a reaction saddle from numerical noise.

    The fixture carries a genuine -312 cm^-1 mode and a -4.21 cm^-1 artefact of the
    integration grid. At the default -10 cm^-1 only the first is imaginary; a user who
    lowers the cutoff sees both.
    """
    out = fixtures / "imaginary_mode.out"
    assert parsers.imaginary_modes(out) == [-312.45]
    assert parsers.imaginary_modes(out, cutoff=-1.0) == [-312.45, -4.21]
    assert parsers.imaginary_modes(out, cutoff=-500.0) == []


# --------------------------------------------------------------------------- geometry
def test_atom_symbols_and_geometry_come_from_the_last_coordinate_block(opt_out):
    syms = parsers.atom_symbols(opt_out)
    assert syms == ANILINE_SYMBOLS
    syms2, xyz = parsers.geometry(opt_out)
    assert syms2 == syms and len(xyz) == 14
    # The optimized nitrogen, not the starting one: the input geometry had it at
    # 2.272019 only after eight cycles of refinement.
    assert xyz[0] == pytest.approx((2.272019, -0.351253, 0.304614), abs=1e-6)


def test_geometry_is_empty_when_the_job_never_printed_coordinates(fixtures):
    assert parsers.geometry(fixtures / "fail_mpi_slots.out") == ([], [])


# ------------------------------------------------------------------------ populations
@pytest.mark.parametrize("scheme", ["Hirshfeld", "Mulliken", "Loewdin"])
def test_populations_sum_to_the_molecular_charge(sp_out, scheme):
    q = parsers.populations(sp_out, scheme)
    assert len(q) == 14
    assert sum(q.values()) == pytest.approx(0.0, abs=5e-4)
    # Nitrogen is atom 0 and carries negative charge in every scheme, though the
    # magnitude is scheme-dependent (Loewdin puts only -0.05 on it where Hirshfeld
    # puts -0.12): a sign test is the strongest claim common to all three.
    assert q[0] < 0.0


@pytest.mark.parametrize("out,charge", [("cat", 1.0), ("ani", -1.0)])
def test_ion_populations_sum_to_the_ion_charge(fixtures, out, charge):
    q = parsers.populations(fixtures / f"ANIL_water_{out}.out", "Hirshfeld")
    assert sum(q.values()) == pytest.approx(charge, abs=5e-4)


def test_populations_are_empty_for_a_scheme_the_job_did_not_print(opt_out):
    # The optimization input asks for Mulliken and Loewdin only; Hirshfeld is requested
    # in the property single point. An absent scheme is empty, never a partial table.
    assert parsers.populations(opt_out, "Hirshfeld") == {}
    assert len(parsers.populations(opt_out, "Mulliken")) == 14


def test_unknown_population_scheme_returns_empty(sp_out):
    assert parsers.populations(sp_out, "NotAScheme") == {}


# ---------------------------------------------------------------------------- dipole
def test_dipole_magnitude(sp_out):
    assert parsers.dipole(sp_out) == pytest.approx(2.452952842)


# -------------------------------------------------------------------- thermochemistry
def test_thermochemistry(opt_out):
    th = parsers.thermochemistry(opt_out)
    assert set(th) == {"zpe_Eh", "enthalpy_Eh", "gibbs_Eh", "entropy_term_Eh",
                       "temperature_K"}
    assert th["temperature_K"] == pytest.approx(298.15)
    assert th["zpe_Eh"] == pytest.approx(0.11677716)
    assert th["gibbs_Eh"] == pytest.approx(-287.16254792)
    # G = H - TS, with ORCA printing the TS term as "Total entropy correction".
    assert th["gibbs_Eh"] == pytest.approx(th["enthalpy_Eh"] + th["entropy_term_Eh"],
                                           abs=1e-6)


def test_thermochemistry_is_empty_without_a_frequency_job(sp_out):
    assert parsers.thermochemistry(sp_out) == {}


# --------------------------------------------------------------------------- walltime
def test_wall_time(opt_out, sp_out):
    assert parsers.wall_time_s(opt_out) == pytest.approx(6 * 60 + 2.139)
    assert parsers.wall_time_s(sp_out) == pytest.approx(45.442)


def test_wall_time_is_none_when_the_job_did_not_finish(fixtures):
    assert parsers.wall_time_s(fixtures / "fail_scf_not_converged.out") is None


# ------------------------------------------------------------- convergence and errors
def test_convergence_failed_detects_both_kinds(fixtures, opt_out):
    assert parsers.convergence_failed(fixtures / "fail_scf_not_converged.out") is True
    assert parsers.convergence_failed(fixtures / "fail_opt_not_converged.out") is True
    assert parsers.convergence_failed(opt_out) is False


def test_mpi_startup_failure_is_not_a_convergence_failure(fixtures):
    # Nothing was computed at all, so reporting it as a convergence problem would send
    # the user to the wrong setting.
    assert parsers.convergence_failed(fixtures / "fail_mpi_slots.out") is False


@pytest.mark.parametrize("name,expected", [
    ("fail_scf_not_converged.out", "SCF failed to converge"),
    ("fail_opt_not_converged.out", "geometry optimization did not converge"),
])
def test_error_message_diagnoses_each_failure(fixtures, name, expected):
    assert parsers.error_message(fixtures / name) == expected


def test_error_message_for_the_mpi_slot_failure_names_the_setting(fixtures):
    msg = parsers.error_message(fixtures / "fail_mpi_slots.out")
    assert "execution.nprocs" in msg
    assert "PHYSICAL core" in msg


def test_error_message_is_empty_for_a_successful_job(real_output):
    assert parsers.error_message(real_output) == ""


def test_error_message_falls_back_to_the_last_lines(tmp_path):
    out = tmp_path / "odd.out"
    out.write_text("something\nunrecognised happened\nand then it stopped\n")
    msg = parsers.error_message(out)
    assert "and then it stopped" in msg


def test_error_message_for_an_empty_output(tmp_path):
    out = tmp_path / "empty.out"
    out.write_text("")
    assert parsers.error_message(out) == "no output produced"


def test_unit_constants():
    assert parsers.HARTREE_TO_EV == pytest.approx(27.211386245988)
    assert math.isclose(parsers.HARTREE_TO_KCAL, 627.509474)
