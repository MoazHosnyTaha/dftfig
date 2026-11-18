"""Config repair, validation and advisories.

The combinations covered here are the ones that have actually cost a run: each either
aborts ORCA at startup, makes it compute the wrong thing, or produces a figure that
misleads. `repair` resolves the ones with a single sensible answer, `validate` refuses
the rest, and `advisories` reports what merely fits this machine badly.
"""
from __future__ import annotations

import pytest

from dftfig import config as C


@pytest.fixture
def four_cores(monkeypatch):
    """Pin the core count: repair's clamp must be tested, not the runner's CPU."""
    monkeypatch.setattr(C, "physical_cores", lambda: 4)
    monkeypatch.setattr(C, "machine", lambda: {"cores": 8, "cores_physical": 4,
                                               "ram_gb": 16.0})
    return 4


# ------------------------------------------------------------------------- the basics
def test_defaults_validate(make_cfg):
    C.validate(make_cfg())


def test_defaults_need_no_repair(make_cfg, four_cores):
    assert C.repair(make_cfg(**{"execution.nprocs": 4})) == []


def test_deep_merge_replaces_a_mapping_valued_parameter_wholesale():
    base = {"molecules": {"definitions": {"A": "C", "B": "N"}, "limit": 0}}
    over = {"molecules": {"definitions": {"C": "O"}}}
    merged = C.deep_merge(base, over)
    # A user who lists one molecule means one, not one plus whatever was in the defaults.
    assert merged["molecules"]["definitions"] == {"C": "O"}
    assert merged["molecules"]["limit"] == 0
    assert base["molecules"]["definitions"] == {"A": "C", "B": "N"}   # not mutated


# ------------------------------------------------------------------------------ nprocs
def test_nprocs_above_the_physical_core_count_is_clamped(make_cfg, four_cores):
    """Open MPI allocates one slot per physical core; asking for more aborts in Startup."""
    cfg = make_cfg(**{"execution.nprocs": 64})
    changed = C.repair(cfg)
    assert cfg["execution"]["nprocs"] == 4
    assert any("nprocs" in c for c in changed)
    C.validate(cfg)


def test_nprocs_at_or_below_the_core_count_is_left_alone(make_cfg, four_cores):
    cfg = make_cfg(**{"execution.nprocs": 2})
    C.repair(cfg)
    assert cfg["execution"]["nprocs"] == 2


# -------------------------------------------------------------------------- dispersion
@pytest.mark.parametrize("functional", ["B3LYP D3BJ", "wB97X-D4", "r2SCAN-3c"])
def test_dispersion_on_a_functional_that_already_has_it_is_removed(make_cfg, four_cores,
                                                                   functional):
    cfg = make_cfg(**{"theory.functional": functional, "theory.dispersion": "D3BJ",
                      "execution.nprocs": 4})
    changed = C.repair(cfg)
    assert cfg["theory"]["dispersion"] == "none"
    assert any("dispersion" in c for c in changed)


@pytest.mark.parametrize("functional", ["B3LYP D3BJ", "wB97X-D4", "r2SCAN-3c"])
def test_double_dispersion_is_a_validation_error_if_not_repaired(make_cfg, functional):
    cfg = make_cfg(**{"theory.functional": functional, "theory.dispersion": "D4"})
    with pytest.raises(C.ConfigError, match="dispersion"):
        C.validate(cfg)


def test_dispersion_on_a_bare_functional_is_kept(make_cfg, four_cores):
    cfg = make_cfg(**{"theory.functional": "B3LYP", "theory.dispersion": "D3BJ",
                      "execution.nprocs": 4})
    C.repair(cfg)
    assert cfg["theory"]["dispersion"] == "D3BJ"


# --------------------------------------------------------------------------- RI / RIJK
def test_rijk_without_a_jk_auxiliary_basis_is_repaired(make_cfg, four_cores):
    cfg = make_cfg(**{"theory.ri_approx": "RIJK", "theory.aux_basis": "def2/J",
                      "execution.nprocs": 4})
    C.repair(cfg)
    assert cfg["theory"]["aux_basis"] == "def2/JK"
    C.validate(cfg)


def test_rijk_with_a_jk_basis_is_untouched(make_cfg, four_cores):
    cfg = make_cfg(**{"theory.ri_approx": "RIJK", "theory.aux_basis": "def2/JK",
                      "execution.nprocs": 4})
    C.repair(cfg)
    assert cfg["theory"]["aux_basis"] == "def2/JK"


def test_rijk_without_a_jk_basis_fails_validation(make_cfg):
    cfg = make_cfg(**{"theory.ri_approx": "RIJK", "theory.aux_basis": "def2/J"})
    with pytest.raises(C.ConfigError, match="/JK"):
        C.validate(cfg)


# ------------------------------------------------------------------- frequencies / opt
def test_frequencies_on_a_loose_optimization_are_repaired(make_cfg, four_cores):
    """A LooseOpt geometry is not converged enough for the harmonic approximation.

    The symptom is spurious imaginary modes, which then send the saddle-point escape
    chasing a minimum that was already there.
    """
    cfg = make_cfg(**{"theory.run_freq": True, "theory.opt_level": "LooseOpt",
                      "execution.nprocs": 4})
    C.repair(cfg)
    assert cfg["theory"]["opt_level"] == "Opt"


def test_loose_opt_without_frequencies_is_allowed(make_cfg, four_cores):
    cfg = make_cfg(**{"theory.run_freq": False, "theory.opt_level": "LooseOpt",
                      "execution.nprocs": 4})
    C.repair(cfg)
    assert cfg["theory"]["opt_level"] == "LooseOpt"
    C.validate(cfg)


def test_frequencies_on_loose_opt_fail_validation(make_cfg):
    cfg = make_cfg(**{"theory.run_freq": True, "theory.opt_level": "LooseOpt"})
    with pytest.raises(C.ConfigError, match="LooseOpt"):
        C.validate(cfg)


# -------------------------------------------------------------------------- population
def test_mulliken_with_a_diffuse_basis_is_switched_to_hirshfeld(make_cfg, four_cores):
    """Mulliken populations with diffuse functions can go negative and corrupt Fukui."""
    cfg = make_cfg(**{"property.population": "Mulliken",
                      "property.basis_sp": "6-311++G(d,p)", "execution.nprocs": 4})
    C.repair(cfg)
    assert cfg["property"]["population"] == "Hirshfeld"


def test_mulliken_without_diffuse_functions_is_kept(make_cfg, four_cores):
    cfg = make_cfg(**{"property.population": "Mulliken",
                      "property.basis_sp": "def2-TZVP", "execution.nprocs": 4})
    C.repair(cfg)
    assert cfg["property"]["population"] == "Mulliken"


def test_mulliken_with_a_diffuse_basis_fails_validation(make_cfg):
    cfg = make_cfg(**{"property.population": "Mulliken",
                      "property.basis_sp": "6-311++G(d,p)"})
    with pytest.raises(C.ConfigError, match="Mulliken"):
        C.validate(cfg)


# ------------------------------------------------------------------- stage dependencies
@pytest.mark.parametrize("stage", ["cubes", "fukui"])
def test_a_stage_that_needs_the_property_job_re_enables_it(make_cfg, four_cores, stage):
    cfg = make_cfg(**{"property.enabled": False, "execution.nprocs": 4})
    cfg["cubes"]["enabled"] = stage == "cubes"
    cfg["fukui"]["enabled"] = stage == "fukui"
    C.repair(cfg)
    assert cfg["property"]["enabled"] is True
    C.validate(cfg)


@pytest.mark.parametrize("stage", ["cubes", "fukui"])
def test_that_dependency_is_a_validation_error_if_not_repaired(make_cfg, stage):
    cfg = make_cfg(**{"property.enabled": False})
    cfg["cubes"]["enabled"] = stage == "cubes"
    cfg["fukui"]["enabled"] = stage == "fukui"
    with pytest.raises(C.ConfigError, match="property.enabled"):
        C.validate(cfg)


def test_figure_2_turns_on_the_cubes_and_a_render_engine(make_cfg, four_cores):
    cfg = make_cfg(**{"figures.which": ["fig2"], "render.engine": "none",
                      "execution.nprocs": 4})
    cfg["cubes"]["enabled"] = False
    C.repair(cfg)
    assert cfg["cubes"]["enabled"] is True
    assert cfg["render"]["engine"] == "pymol"


# ---------------------------------------------------------------------------- solvation
def test_no_continuum_model_means_gas_phase_only(make_cfg, four_cores):
    cfg = make_cfg(**{"solvation.model": "none", "solvation.solvents": ["water"],
                      "execution.nprocs": 4})
    C.repair(cfg)
    assert cfg["solvation"]["solvents"] == ["gas"]
    assert cfg["figures"]["reference_solvent"] == "gas"
    C.validate(cfg)


def test_the_figure_solvent_must_be_one_that_was_calculated(make_cfg, four_cores):
    cfg = make_cfg(**{"solvation.solvents": ["dmso"], "execution.nprocs": 4})
    cfg["figures"]["reference_solvent"] = "water"
    C.repair(cfg)
    assert cfg["figures"]["reference_solvent"] == "dmso"


def test_a_figure_solvent_outside_the_solvent_list_fails_validation(make_cfg):
    cfg = make_cfg(**{"solvation.solvents": ["dmso"]})
    cfg["figures"]["reference_solvent"] = "water"
    with pytest.raises(C.ConfigError, match="reference_solvent"):
        C.validate(cfg)


# ------------------------------------------------------------------- registry checking
def test_validate_rejects_a_value_outside_the_registry_choices(make_cfg):
    cfg = make_cfg(**{"theory.functional": "PM3"})
    with pytest.raises(C.ConfigError, match="not one of"):
        C.validate(cfg)


def test_validate_rejects_a_wrong_type(make_cfg):
    cfg = make_cfg(**{"execution.nprocs": "four"})
    with pytest.raises(C.ConfigError, match="expected a number"):
        C.validate(cfg)


def test_validate_rejects_a_value_below_the_minimum(make_cfg):
    cfg = make_cfg(**{"execution.nprocs": 0})
    with pytest.raises(C.ConfigError, match="minimum"):
        C.validate(cfg)


def test_validate_rejects_an_empty_molecule_set(make_cfg):
    cfg = make_cfg(**{"molecules.definitions": {}})
    with pytest.raises(C.ConfigError, match="nothing to calculate"):
        C.validate(cfg)


def test_validate_reports_a_missing_sdf_directory(make_cfg, tmp_path):
    cfg = make_cfg(**{"molecules.structure_source": "sdf",
                      "molecules.sdf_dir": str(tmp_path / "absent")})
    with pytest.raises(C.ConfigError, match="does not exist"):
        C.validate(cfg)


def test_validate_reports_several_problems_at_once(make_cfg):
    cfg = make_cfg(**{"theory.ri_approx": "RIJK", "theory.aux_basis": "def2/J",
                      "theory.run_freq": True, "theory.opt_level": "LooseOpt"})
    with pytest.raises(C.ConfigError) as exc:
        C.validate(cfg)
    assert "/JK" in str(exc.value) and "LooseOpt" in str(exc.value)


# --------------------------------------------------------------------------- advisories
def test_advisory_when_nprocs_exceeds_the_physical_cores(make_cfg, four_cores):
    # repair() normally clamps this; advisories() has to say so too, because a config
    # written by hand and passed straight to `dft validate` never goes through repair.
    msgs = C.advisories(make_cfg(**{"execution.nprocs": 12}))
    assert any("not enough slots" in m for m in msgs)


def test_advisory_when_the_jobs_together_oversubscribe_the_machine(make_cfg, four_cores):
    msgs = C.advisories(make_cfg(**{"execution.nprocs": 4,
                                    "execution.max_parallel_jobs": 4}))
    assert any("contend" in m for m in msgs)


def test_advisory_when_the_memory_request_exceeds_the_machine(make_cfg, four_cores):
    msgs = C.advisories(make_cfg(**{"execution.nprocs": 4, "execution.maxcore_mb": 8000,
                                    "execution.max_parallel_jobs": 2}))
    assert any("maxcore_mb" in m for m in msgs)


@pytest.mark.parametrize("basis", ["6-311++G(d,p)", "def2-TZVPD"])
@pytest.mark.parametrize("orbitals", ["homo_lumo", "homo_lumo_pm1", "homo_lumo_pm2"])
def test_advisory_for_a_diffuse_property_basis_with_a_lumo_cube(make_cfg, four_cores,
                                                                basis, orbitals):
    """Every cubes.orbitals setting plots the LUMO, so every one earns the advisory.

    The LUMO of a neutral organic molecule is usually unbound, and a diffuse set
    describes that unbound state as a discretised continuum function - a cloud detached
    from the molecule rather than the pi* orbital the figure is supposed to show.
    """
    cfg = make_cfg(**{"property.basis_sp": basis, "execution.nprocs": 4})
    cfg["cubes"]["enabled"] = True
    cfg["cubes"]["orbitals"] = orbitals
    assert any("diffuse functions" in m for m in C.advisories(cfg))


def test_no_lumo_advisory_without_diffuse_functions(make_cfg, four_cores):
    cfg = make_cfg(**{"property.basis_sp": "def2-TZVP", "execution.nprocs": 4})
    cfg["cubes"]["enabled"] = True
    assert not any("diffuse functions" in m for m in C.advisories(cfg))


def test_advisory_when_no_orca_was_found(make_cfg, four_cores):
    msgs = C.advisories(make_cfg(**{"execution.orca_path": "", "execution.nprocs": 4}))
    assert any("no ORCA binary found" in m for m in msgs)


def test_advisory_when_orca_path_is_not_really_orca(make_cfg, four_cores, tmp_path):
    fake = tmp_path / "orca"
    fake.write_text("#!/bin/sh\necho screen reader\n")
    fake.chmod(0o755)
    msgs = C.advisories(make_cfg(**{"execution.orca_path": str(fake),
                                    "execution.nprocs": 4}))
    assert any("screen reader" in m for m in msgs)


def test_advisory_when_a_fixed_charge_is_forced_on_every_compound(make_cfg, four_cores):
    msgs = C.advisories(make_cfg(**{"molecules.charge_mode": "fixed",
                                    "molecules.charge": 1, "execution.nprocs": 4}))
    assert any("charge_mode" in m for m in msgs)


def test_advisory_for_an_anion_without_diffuse_functions_in_the_fukui_step(make_cfg,
                                                                          four_cores):
    cfg = make_cfg(**{"molecules.definitions": {"BENZOATE": "[O-]C(=O)c1ccccc1"},
                      "property.basis_sp": "def2-TZVP", "execution.nprocs": 4})
    cfg["fukui"]["enabled"] = True
    assert any("unbound" in m for m in C.advisories(cfg))


# ----------------------------------------------------------------------------- machine
def test_physical_cores_is_at_least_one_and_no_more_than_the_thread_count():
    import os
    phys = C.physical_cores()
    assert 1 <= phys <= (os.cpu_count() or 1)


def test_machine_reports_both_core_counts():
    m = C.machine()
    assert m["cores_physical"] <= m["cores"]


def test_estimate_counts_the_jobs_a_run_will_submit(make_cfg):
    cfg = make_cfg(**{"molecules.definitions": {"A": "c1ccccc1", "B": "Nc1ccccc1"},
                      "solvation.solvents": ["water", "dmso"]})
    est = C.estimate(cfg)
    # 2 molecules x 2 media: one opt, one property single point, two ion single points.
    assert est["opt_freq_jobs"] == 4
    assert est["property_jobs"] == 4
    assert est["fukui_jobs"] == 8
    assert est["total_orca_jobs"] == 16


# ---------------------------------------------------------------------- save and reload
def test_save_then_load_round_trips(make_cfg, tmp_path, monkeypatch):
    monkeypatch.setattr(C, "resolve_executables", lambda cfg: [])
    # nprocs = 1 so that load()'s repair step has nothing to clamp: a CI runner with two
    # physical cores would otherwise rewrite the default 4 and the round trip would
    # "fail" for a reason that has nothing to do with saving and loading.
    cfg = make_cfg(**{"execution.nprocs": 1})
    path = tmp_path / "config.yaml"
    C.save(cfg, path)
    assert C.load(path) == cfg


def test_load_fills_missing_keys_from_the_registry(tmp_path, monkeypatch):
    monkeypatch.setattr(C, "resolve_executables", lambda cfg: [])
    path = tmp_path / "config.yaml"
    path.write_text("molecules:\n  structure_source: smiles\n"
                    "  definitions:\n    ANIL: Nc1ccccc1\n")
    cfg = C.load(path)
    assert cfg["theory"]["functional"] == "B3LYP"
    assert cfg["molecules"]["definitions"] == {"ANIL": "Nc1ccccc1"}


def test_get_reads_a_dotted_key(make_cfg):
    assert C.get(make_cfg(), "theory.functional") == "B3LYP"
