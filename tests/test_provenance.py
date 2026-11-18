"""Provenance: the config hash, the run manifest, and reading ORCA's own version.

None of these tests start ORCA. The fixture in tests/data is a real ORCA 6.1.0 output
(a B3LYP-D3BJ/def2-TZVP CPCM(Water) single point on aniline) precisely so the parsers
are exercised against what ORCA actually prints rather than a hand-written sample.
"""
from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from dftfig import __version__
from dftfig import params as P
from dftfig import parsers, provenance

FIXTURE = Path(__file__).parent / "data" / "ANIL_water_sp.out"


@pytest.fixture
def cfg(tmp_path):
    """Registry defaults pointed at a temporary run directory.

    `params.defaults()` rather than `config.load()`: the hash is pure dictionary logic,
    and load() searches the machine for ORCA and PyMOL, which would make the result
    depend on whoever is running the suite.
    """
    c = P.defaults()
    c["execution"]["work_dir"] = str(tmp_path / "runs")
    c["execution"]["scratch_dir"] = str(tmp_path / "scratch")
    c["molecules"]["structure_source"] = "smiles"
    c["molecules"]["definitions"] = {"ANIL": "Nc1ccccc1"}
    return c


# ------------------------------------------------------------------ the ORCA version
def test_orca_version_parsed_from_a_real_output():
    assert parsers.program_version(FIXTURE) == "6.1.0"
    assert provenance.orca_version_of(FIXTURE) == "6.1.0"


def test_orca_version_is_empty_when_absent(tmp_path):
    p = tmp_path / "truncated.out"
    p.write_text("ORCA died before printing its banner\n")
    assert parsers.program_version(p) == ""


def test_keyword_line_comes_back_verbatim():
    kw = parsers.input_keyword_line(FIXTURE)
    assert kw == ("! B3LYP D3BJ def2-TZVP RIJCOSX def2/J SP TightSCF DefGrid2 "
                  "CPCM(Water)")


def test_keyword_lines_prefer_the_output_over_the_input(cfg):
    """The line that ran, not the line that was written.

    The runner rewrites an input in place when MPI refuses the requested slot count, and
    a saddle-point restart appends keywords, so when both files exist the output's echo
    is the authority.
    """
    d = Path(cfg["execution"]["work_dir"]) / "ANIL" / "water"
    d.mkdir(parents=True)
    (d / "ANIL_water_sp.out").write_text(FIXTURE.read_text())
    (d / "ANIL_water_sp.inp").write_text("! SOMETHING ELSE ENTIRELY\n")
    (d / "ANIL_water_opt.inp").write_text("! B3LYP def2-SVP Opt FREQ TightSCF\n")

    kw = provenance.keyword_lines(cfg, "ANIL", "water")
    assert "def2-TZVP" in kw["orca_keywords_sp"]
    assert "SOMETHING ELSE" not in kw["orca_keywords_sp"]
    # No opt output yet - the input is the only record, and is used.
    assert kw["orca_keywords_opt"] == "! B3LYP def2-SVP Opt FREQ TightSCF"


# -------------------------------------------------------------------- the config hash
@pytest.mark.parametrize("section,key,value", [
    ("execution", "nprocs", 1),
    ("execution", "maxcore_mb", 16000),
    ("execution", "max_parallel_jobs", 8),
    ("execution", "work_dir", "/somewhere/else"),
    ("execution", "scratch_dir", "/fast/scratch"),
    ("execution", "orca_path", "/opt/orca_6_1_0/orca"),
    ("execution", "resume", False),
    ("cubes", "esp_nprocs", 4),
    ("render", "pymol_path", "/usr/bin/pymol"),
    ("render", "width", 2400),
    ("render", "orbital_isovalue", 0.05),
    ("figures", "dpi", 600),
    ("figures", "which", ["fig2"]),
])
def test_config_sha_ignores_what_cannot_change_a_number(cfg, section, key, value):
    """The whole point of the hash: same chemistry on a laptop and a cluster hashes alike."""
    before = provenance.config_sha(cfg)
    other = copy.deepcopy(cfg)
    other[section][key] = value
    assert other[section][key] != cfg[section][key], "the test changed nothing"
    assert provenance.config_sha(other) == before


@pytest.mark.parametrize("section,key,value", [
    ("theory", "functional", "PBE0"),
    ("theory", "basis_opt", "def2-SVP"),
    ("theory", "dispersion", "D4"),
    ("theory", "grid", "DefGrid3"),
    ("theory", "scf_conv", "VeryTightSCF"),
    ("property", "basis_sp", "ma-def2-TZVP"),
    ("property", "population", "Loewdin"),
    ("solvation", "model", "SMD"),
    ("solvation", "solvents", ["dmso"]),
    ("cubes", "grid_points", 120),
    ("molecules", "charge", -1),
    ("conformers", "random_seed", 1234),
])
def test_config_sha_changes_with_anything_that_can(cfg, section, key, value):
    before = provenance.config_sha(cfg)
    other = copy.deepcopy(cfg)
    other[section][key] = value
    assert other[section][key] != cfg[section][key], "the test changed nothing"
    assert provenance.config_sha(other) != before


def test_config_sha_is_insensitive_to_key_order(cfg):
    """Serialized with sorted keys, so a re-ordered config.yaml is the same config."""
    shuffled = {s: dict(reversed(list(v.items()))) if isinstance(v, dict) else v
                for s, v in reversed(list(cfg.items()))}
    assert provenance.config_sha(shuffled) == provenance.config_sha(cfg)


def test_config_sha_is_short_and_hexadecimal(cfg):
    sha = provenance.config_sha(cfg)
    assert len(sha) == 12
    assert all(c in "0123456789abcdef" for c in sha)


# ----------------------------------------------------------------- the run manifest
def _populate(cfg) -> Path:
    """A minimal but realistic run directory: one input, one output, a cube and a figure."""
    work = Path(cfg["execution"]["work_dir"])
    d = work / "ANIL" / "water"
    d.mkdir(parents=True)
    (d / "ANIL_water_sp.out").write_text(FIXTURE.read_text())
    (d / "ANIL_water_sp.inp").write_text("! B3LYP def2-TZVP SP\n")
    (d / "ANIL_water_opt.xyz").write_text("1\n\nH 0.0 0.0 0.0\n")
    (d / "homo.cube").write_text("cube\n" * 10)
    (work / "descriptors.csv").write_text("molecule,medium\nANIL,water\n")
    figs = work / "figures"
    figs.mkdir()
    (figs / "figure2.png").write_bytes(b"\x89PNG not really")
    return work


def test_manifest_round_trips_through_verify(cfg):
    _populate(cfg)
    path = provenance.write(cfg)
    assert path.exists()

    man = json.loads(path.read_text())
    assert man["dftfig_version"] == __version__
    assert man["orca_version"] == "6.1.0"          # read from the output, not from ORCA
    assert man["config_sha"] == provenance.config_sha(cfg)
    assert man["generated_utc"].endswith("Z")
    assert man["config"]["theory"]["functional"] == cfg["theory"]["functional"]

    listed = {e["path"] for e in man["files"]}
    assert "ANIL/water/ANIL_water_sp.out" in listed
    assert "ANIL/water/homo.cube" in listed
    assert "figures/figure2.png" in listed
    assert "descriptors.csv" in listed
    assert provenance.MANIFEST_NAME not in listed   # it cannot contain its own hash
    for e in man["files"]:
        assert len(e["sha256"]) == 64
        assert e["bytes"] == (Path(cfg["execution"]["work_dir"]) / e["path"]).stat().st_size

    r = provenance.verify(cfg)
    assert (r["missing"], r["changed"], r["unlisted"]) == ([], [], [])
    assert r["ok"] == len(man["files"])
    assert provenance.verify_report(cfg) == 0


def test_verify_detects_an_altered_file(cfg):
    work = _populate(cfg)
    provenance.write(cfg)
    (work / "descriptors.csv").write_text("molecule,medium\nANIL,dmso\n")

    r = provenance.verify(cfg)
    assert r["changed"] == ["descriptors.csv"]
    assert provenance.verify_report(cfg) != 0


def test_verify_detects_a_missing_file(cfg):
    work = _populate(cfg)
    provenance.write(cfg)
    (work / "figures" / "figure2.png").unlink()

    r = provenance.verify(cfg)
    assert r["missing"] == ["figures/figure2.png"]
    assert provenance.verify_report(cfg) != 0


def test_verify_detects_a_file_nobody_recorded(cfg):
    work = _populate(cfg)
    provenance.write(cfg)
    (work / "figures" / "figure13.png").write_bytes(b"added afterwards")

    r = provenance.verify(cfg)
    assert r["unlisted"] == ["figures/figure13.png"]
    assert (r["missing"], r["changed"]) == ([], [])
    assert provenance.verify_report(cfg) != 0


def test_verify_ignores_files_it_never_claimed_to_track(cfg):
    """A reader's own notes next to the results are not a provenance failure."""
    work = _populate(cfg)
    provenance.write(cfg)
    (work / "notes.txt").write_text("remember to redo the dmso run")
    (work / "ANIL" / "water" / "ANIL_water_sp.gbw").write_bytes(b"\x00" * 64)

    r = provenance.verify(cfg)
    assert (r["missing"], r["changed"], r["unlisted"]) == ([], [], [])


def test_verify_exits_non_zero_without_a_manifest(cfg):
    _populate(cfg)
    with pytest.raises(FileNotFoundError):
        provenance.verify(cfg)
    assert provenance.verify_report(cfg) == 2


def test_verify_through_the_dft_command(cfg, tmp_path, monkeypatch):
    """`dft verify` is meant for scripts, so the exit code is the contract."""
    from dftfig import cli
    from dftfig import config as C

    _populate(cfg)
    # An explicit (if fictitious) ORCA path stops config.load() searching this machine
    # and caching whatever it finds in the user's ~/.config.
    cfg["execution"]["orca_path"] = str(tmp_path / "orca")
    cfg["render"]["pymol_path"] = str(tmp_path / "pymol")
    cfg_path = tmp_path / "config.yaml"
    C.save(cfg, cfg_path)
    provenance.write(cfg)

    assert cli.main(["verify", "-c", str(cfg_path)]) == 0

    (Path(cfg["execution"]["work_dir"]) / "descriptors.csv").write_text("tampered\n")
    assert cli.main(["verify", "-c", str(cfg_path)]) == 1
