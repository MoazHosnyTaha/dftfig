"""Full build + opt stages with execution.dry_run, on a machine with no ORCA.

dry_run generates every input file and skips only the subprocess call, so the part of
the pipeline that decides WHAT ORCA will compute - keyword line, parallel and memory
blocks, solvation, and above all the charge and multiplicity of each compound - is
testable without a quantum chemistry program installed. A wrong charge is the failure
that costs the most: it changes the electron count, so every orbital, descriptor and
Fukui index that follows describes a different species, and nothing about the output
looks wrong.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from dftfig import __main__ as cli
from dftfig import molecules as M
from dftfig import orca_io, runner

#: A neutral molecule and a cation in one run, which is what per-compound charge
#: detection exists for.
MOLECULES = {"ANIL": "Nc1ccccc1", "BENZAMM": "[NH3+]Cc1ccccc1"}


@pytest.fixture
def project(make_cfg):
    cfg = make_cfg(**{
        "molecules.definitions": MOLECULES,
        "conformers.n_confs": 5,          # the geometry is not what is under test here
        "execution.nprocs": 4,
        "execution.maxcore_mb": 900,
        "execution.max_parallel_jobs": 1,
        "cubes.enabled": False,           # no %plots block in the way
    })
    specs = M.load_molecules(cfg)
    return cfg, specs


@pytest.fixture
def built(project):
    cfg, specs = project
    cli.stage_build(cfg, specs)
    return cfg, specs


def inp_text(cfg, name, solvent, stage) -> str:
    path = (Path(cfg["execution"]["work_dir"]) / name / solvent
            / f"{name}_{solvent}_{stage}.inp")
    assert path.is_file(), f"{path} was not written"
    return path.read_text()


# ------------------------------------------------------------------------------ build
def test_build_writes_one_starting_geometry_per_molecule(built):
    cfg, specs = built
    work = Path(cfg["execution"]["work_dir"])
    for name in specs:
        xyz = work / name / f"{name}_start.xyz"
        assert xyz.is_file()
        assert int(xyz.read_text().splitlines()[0]) == specs[name].n_atoms


def test_charge_and_multiplicity_are_detected_per_compound(project):
    _, specs = project
    assert specs["ANIL"].charge == 0 and specs["ANIL"].multiplicity == 1
    assert specs["BENZAMM"].charge == 1 and specs["BENZAMM"].multiplicity == 1
    assert specs["BENZAMM"].charge_origin == "detected from formal charges"


# -------------------------------------------------------------------------------- opt
def test_opt_stage_writes_inputs_and_runs_nothing(built):
    cfg, specs = built
    results = cli.stage_opt(cfg, specs)
    assert len(results) == 2
    assert all(r["ok"] and r.get("dry_run") for r in results)
    work = Path(cfg["execution"]["work_dir"])
    assert sorted(p.name for p in work.rglob("*.inp")) == [
        "ANIL_water_opt.inp", "BENZAMM_water_opt.inp"]
    # Nothing ran: no ORCA output anywhere under the project.
    assert list(work.rglob("*.out")) == []


def test_opt_keyword_line(built):
    cfg, specs = built
    cli.stage_opt(cfg, specs)
    first = inp_text(cfg, "ANIL", "water", "opt").splitlines()[0]
    assert first == ("! B3LYP 6-31+G(d,p) RIJCOSX def2/J Opt FREQ TightSCF DefGrid2 "
                     "CPCM(Water)")


def test_parallel_and_memory_blocks(built):
    cfg, specs = built
    cli.stage_opt(cfg, specs)
    txt = inp_text(cfg, "ANIL", "water", "opt")
    assert "%pal nprocs 4 end" in txt
    assert "%maxcore 900" in txt
    # ORCA's %maxcore is per process, so the total is nprocs x maxcore - that is the
    # number advisories() compares against the machine's RAM.
    assert "%maxcore 3600" not in txt


def test_frequency_block_carries_the_temperature(built):
    cfg, specs = built
    cli.stage_opt(cfg, specs)
    assert "%freq  Temp 298.15  end" in inp_text(cfg, "ANIL", "water", "opt")


def test_charge_and_multiplicity_reach_the_coordinate_line(built):
    cfg, specs = built
    cli.stage_opt(cfg, specs)
    assert "* xyzfile 0 1 ANIL_start.xyz" in inp_text(cfg, "ANIL", "water", "opt")
    assert "* xyzfile 1 1 BENZAMM_start.xyz" in inp_text(cfg, "BENZAMM", "water", "opt")


def test_the_starting_geometry_is_copied_beside_the_input(built):
    cfg, specs = built
    cli.stage_opt(cfg, specs)
    work = Path(cfg["execution"]["work_dir"])
    assert (work / "ANIL" / "water" / "ANIL_start.xyz").is_file()


def test_a_second_dry_run_is_idempotent(built):
    cfg, specs = built
    cli.stage_opt(cfg, specs)
    before = inp_text(cfg, "ANIL", "water", "opt")
    cli.stage_opt(cfg, specs)
    assert inp_text(cfg, "ANIL", "water", "opt") == before


# ------------------------------------------------------------------------- solvation
def test_smd_is_a_block_not_a_keyword(built):
    cfg, specs = built
    cfg["solvation"]["model"] = "SMD"
    cli.stage_opt(cfg, specs)
    txt = inp_text(cfg, "ANIL", "water", "opt")
    assert "CPCM(Water)" not in txt.splitlines()[0]
    assert "smd true" in txt and 'SMDsolvent "Water"' in txt


def test_gas_phase_has_no_solvation_block(built):
    cfg, specs = built
    cfg["solvation"]["solvents"] = ["gas"]
    cfg["figures"]["reference_solvent"] = "gas"
    cli.stage_opt(cfg, specs)
    txt = inp_text(cfg, "ANIL", "gas", "opt")
    assert "CPCM" not in txt and "smd" not in txt


# --------------------------------------------------------- single points and the ions
def test_the_property_single_point_uses_the_property_basis(built):
    cfg, specs = built
    path = orca_io.write_input(cfg, specs["ANIL"], "water", "sp",
                               geom_xyz="ANIL_water_opt.xyz")
    first = path.read_text().splitlines()[0]
    assert "def2-TZVP" in first and "6-31+G(d,p)" not in first
    assert first.endswith("CPCM(Water)")
    assert " SP " in first
    # Hirshfeld is what descriptors.py reads for the condensed Fukui functions.
    assert "Print[P_Hirshfeld] 1" in path.read_text()


@pytest.mark.parametrize("stage,charge,mult", [("cat", 1, 2), ("ani", -1, 2)])
def test_the_ion_single_points_are_unrestricted_and_correctly_charged(
        built, stage, charge, mult):
    cfg, specs = built
    path = orca_io.write_input(cfg, specs["ANIL"], "water", stage,
                               geom_xyz="ANIL_water_opt.xyz")
    txt = path.read_text()
    # N-1 and N+1 electrons: open shell, so UKS, and a doublet.
    assert "UKS" in txt.splitlines()[0]
    assert f"* xyzfile {charge} {mult} ANIL_water_opt.xyz" in txt


def test_the_cation_of_a_cation_is_doubly_charged(built):
    cfg, specs = built
    path = orca_io.write_input(cfg, specs["BENZAMM"], "water", "cat",
                               geom_xyz="BENZAMM_water_opt.xyz")
    assert "* xyzfile 2 2 BENZAMM_water_opt.xyz" in path.read_text()


# ------------------------------------------------------------------- keyword assembly
def test_a_composite_method_brings_its_own_basis_set(make_cfg):
    cfg = make_cfg(**{"theory.functional": "r2SCAN-3c"})
    line = orca_io.keyword_line(cfg, "opt", "water")
    assert "r2SCAN-3c" in line
    assert "6-31+G(d,p)" not in line and "def2-SVP" not in line


def test_an_open_shell_neutral_also_gets_uks(make_cfg):
    cfg = make_cfg()
    assert "UKS" in orca_io.keyword_line(cfg, "opt", "water", multiplicity=3)
    assert "UKS" not in orca_io.keyword_line(cfg, "opt", "water", multiplicity=1)


def test_an_extra_keyword_replaces_its_family_rather_than_duplicating_it(built):
    """A saddle-point restart asks for TightOpt; "Opt ... TightOpt" is ambiguous input."""
    cfg, specs = built
    path = orca_io.write_input(cfg, specs["ANIL"], "water", "opt",
                               extra_keywords=["TightOpt", "DefGrid3"],
                               basename="ANIL_water_opt_r1")
    words = path.read_text().splitlines()[0].split()
    assert "TightOpt" in words and "Opt" not in words
    assert "DefGrid3" in words and "DefGrid2" not in words
    assert path.name == "ANIL_water_opt_r1.inp"


def test_cubes_add_a_plots_block_with_the_right_orbital_indices(built):
    cfg, specs = built
    cfg["cubes"]["enabled"] = True
    cfg["cubes"]["orbitals"] = "homo_lumo"
    cfg["cubes"]["grid_points"] = 60
    txt = orca_io.write_input(cfg, specs["ANIL"], "water", "sp",
                              geom_xyz="ANIL_water_opt.xyz").read_text()
    # Aniline has 50 electrons, so the HOMO is MO 24 and the LUMO MO 25.
    assert '  MO("homo.cube", 24, 0);' in txt
    assert '  MO("lumo.cube", 25, 0);' in txt
    assert "  dim1 60" in txt


# ------------------------------------------------------------------------- the runner
def test_run_job_in_dry_run_mode_never_looks_for_orca(make_cfg, tmp_path):
    cfg = make_cfg(**{"execution.orca_path": "/definitely/not/here/orca"})
    inp = tmp_path / "x.inp"
    inp.write_text("! HF STO-3G\n* xyz 0 1\nH 0 0 0\nH 0 0 0.74\n*\n")
    res = runner.run_job(cfg, inp, "dry")
    assert res["ok"] and res["dry_run"]
    assert not inp.with_suffix(".out").exists()


def test_run_job_reports_a_missing_orca_when_not_in_dry_run(make_cfg, tmp_path):
    cfg = make_cfg(**{"execution.dry_run": False,
                      "execution.orca_path": "/definitely/not/here/orca"})
    inp = tmp_path / "x.inp"
    inp.write_text("! HF STO-3G\n")
    res = runner.run_job(cfg, inp, "no orca")
    assert res["ok"] is False
    assert "ORCA not found" in res["error"]


def test_a_finished_job_is_reused_rather_than_rerun(make_cfg, fixtures, tmp_path):
    """execution.resume is what makes an interrupted overnight batch cheap to restart."""
    cfg = make_cfg(**{"execution.dry_run": False, "execution.resume": True})
    inp = tmp_path / "ANIL_water_opt.inp"
    inp.write_text("! B3LYP\n")
    inp.with_suffix(".out").write_text((fixtures / "ANIL_water_opt.out").read_text())
    res = runner.run_job(cfg, inp, "reuse")
    assert res["ok"] and res["reused"]
    assert res["wall_s"] == pytest.approx(362.139)


def test_strip_pal_is_what_the_mpi_fallback_does(tmp_path):
    inp = tmp_path / "x.inp"
    inp.write_text("! B3LYP SP\n%pal nprocs 12 end\n%maxcore 900\n* xyz 0 1\n*\n")
    runner._set_pal(inp, 4)
    assert "%pal nprocs 4 end" in inp.read_text()
    runner._set_pal(inp, 1)
    txt = inp.read_text()
    assert "%pal" not in txt and "%maxcore 900" in txt


@pytest.mark.orca
def test_a_real_orca_job_runs(make_cfg, tmp_path):
    """Only with -m orca: this is the one test that needs the program itself."""
    from dftfig import discover, parsers

    cfg = make_cfg(**{"execution.dry_run": False, "execution.nprocs": 1,
                      "execution.orca_path": str(discover.find_orca())})
    inp = tmp_path / "h2.inp"
    inp.write_text("! HF STO-3G SP\n* xyz 0 1\nH 0 0 0\nH 0 0 0.74\n*\n")
    res = runner.run_job(cfg, inp, "H2")
    assert res["ok"], res.get("error")
    assert parsers.final_energy(inp.with_suffix(".out")) == pytest.approx(-1.117, abs=0.02)
