# Contributing to dftfig

Contributions are welcome, including from people who have never contributed to a
scientific Python package before. If something here is unclear, that is a bug in this
document — please open an issue saying so.

## Getting help, reporting a problem, asking for a feature

| what | where |
|---|---|
| Something does not work | [open an issue](https://github.com/MoazHosnyTaha/dftfig/issues) |
| A question about how to use it | open an issue with the `question` label |
| A figure looks wrong | open an issue and attach the figure plus `runs/descriptors.csv` |
| A chemistry question about the results | open an issue — these are often the most useful reports |

**Before opening an issue about a failure, run `dft doctor` and paste its output.** It
reports your ORCA and PyMOL paths, core counts and installed packages, which is what
almost every diagnosis starts from.

If ORCA itself failed, attach the relevant `runs/<molecule>/<medium>/*.out`. The
failure is usually in the last 40 lines.

## Setting up a development environment

```bash
git clone https://github.com/MoazHosnyTaha/dftfig && cd dftfig
conda env create -f environment.yml
conda activate dftfig
pip install -e ".[dev]"
pytest
```

You do **not** need ORCA installed to work on most of the package. The test suite runs
against checked-in ORCA output files, and `execution.dry_run: true` generates every
input file without launching anything. Tests that genuinely need ORCA are marked:

```bash
pytest -m "not orca"      # what CI runs; the default on a machine without ORCA
pytest -m orca            # only the ones that need it
```

Tests that render a 3D panel are marked `pymol` and skip themselves when PyMOL is not
importable, so a pip-only environment still runs everything else.

## What a good pull request looks like

1. **One change per pull request.** A packaging fix and a new figure are two.
2. **A test that fails before and passes after.** For a parser fix, that usually means
   adding a trimmed ORCA output to `tests/fixtures/` that reproduces the problem.
3. **Explain the chemistry, not just the code.** If you change a default, say what
   physical argument supports the new value. Defaults in this package are documented
   with reasons, and reviewers will ask for yours.
4. `ruff check .` passes.

## Adding a figure

Figures are plugins, and the interface is documented in `dftfig/figureapi.py`. A figure
module declares `NAME`, `TITLE` and `REQUIRES`, and exposes `build(ctx)` taking a
`FigureContext` and returning the list of files it wrote. It reads nothing but that
context: no ORCA output, no path arithmetic, no `cfg["execution"]["work_dir"]`.

Start from [`docs/writing-a-figure.md`](docs/writing-a-figure.md), which builds one from
scratch, and from `examples/fukui_bars/` which is the same figure ready to run. The four
built-ins under `dftfig/figures/` are the reference implementations: fig12 for a figure
that needs only the structures and the property table, fig13 for one colour scale shared
across rendered panels, fig2 for the full case.

You do not need a pull request to use your own figure - drop it in `./dft_figures` or
ship it in a package declaring the `dftfig.figures` entry-point group. Send a pull
request when the figure is one other people in your field would want.

If your figure needs data the pipeline does not currently expose, say so in the issue
rather than reaching around the interface — that is a gap in the contract, and other
figures probably need the same thing.

## Things that would genuinely help

- **Windows support for the arrow-key menus.** The key reader is termios-based, so
  Windows currently falls back to numbered prompts. An `msvcrt` backend in
  `dftfig/menu.py` would fix it. It is deliberately not written yet because it cannot
  be tested on the maintainer's machines.
- **More pKa reference data.** The SMARTS model in `dftfig/pka.py` is validated against
  17 compounds in `tests/test_pka.py` (15 within 1.5 log units, carbazole correctly
  non-ionizable, warfarin a documented xfail). Adding well-measured values, especially for
  classes it currently gets wrong, directly improves Figure 12.
- **Basis sets and elements beyond organics.** ECPs and transition metals are untested.
- **Other quantum chemistry backends.** The ORCA-specific code is confined to
  `orca_io.py`, `parsers.py`, `runner.py` and `cubes.py`.

## Code of conduct

By participating you agree to abide by the [Code of Conduct](CODE_OF_CONDUCT.md).
