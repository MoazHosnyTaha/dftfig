

def _write(tmp_path, kw_in, kw_out):
    inp = tmp_path / "j.inp"
    out = tmp_path / "j.out"
    inp.write_text(f"{kw_in}\n*xyz 0 1\nH 0 0 0\n*\n")
    out.write_text("INPUT FILE\n|  1> " + kw_out + "\n|  2> *xyz 0 1\n\n"
                   "****ORCA TERMINATED NORMALLY****\nTOTAL RUN TIME: 0 sec\n")
    return inp, out


def test_an_unchanged_keyword_line_still_resumes(tmp_path):
    from dftfig.runner import _keywords_changed
    inp, out = _write(tmp_path, "! B3LYP def2-SVP Opt", "! B3LYP def2-SVP Opt")
    assert _keywords_changed(inp, out) is None


def test_whitespace_alone_does_not_force_a_recompute(tmp_path):
    from dftfig.runner import _keywords_changed
    inp, out = _write(tmp_path, "!  B3LYP   def2-SVP  Opt", "! B3LYP def2-SVP Opt")
    assert _keywords_changed(inp, out) is None


def test_a_changed_functional_is_detected(tmp_path):
    """Changing property.functional_sp must recompute, not reuse the old numbers."""
    from dftfig.runner import _keywords_changed
    inp, out = _write(tmp_path, "! wB97X-D4 def2-TZVP SP", "! B3LYP def2-TZVP SP")
    changed = _keywords_changed(inp, out)
    assert changed is not None
    assert changed[0].startswith("! B3LYP") and changed[1].startswith("! wB97X-D4")


def test_an_output_with_no_echo_is_left_alone(tmp_path):
    """Never force a recompute on evidence we do not have."""
    from dftfig.runner import _keyword_line_of, _keywords_changed
    inp = tmp_path / "j.inp"
    inp.write_text("! B3LYP def2-SVP Opt\n*xyz 0 1\nH 0 0 0\n*\n")
    out = tmp_path / "j.out"
    out.write_text("****ORCA TERMINATED NORMALLY****\n")
    assert _keyword_line_of(out) is None
    assert _keywords_changed(inp, out) is None
