# Changelog

All notable changes to this project are documented here.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project uses [semantic versioning](https://semver.org/).

## [Unreleased]

### Added
- A documented, versioned figure plugin API (`dftfig/figureapi.py`, `API_VERSION = 1`).
  A figure is now a module that declares `NAME`, `TITLE` and `REQUIRES` and exposes
  `build(ctx)`, where `ctx` is a `FigureContext` carrying the resolved config, the
  molecules (RDKit mol in a defined atom order, plus charge, multiplicity, formula and
  run directory), the DFT descriptor table, the resolved property table, the output
  directory and the provenance block - with accessors for the things figures used to
  work out for themselves: locating a cube file, obtaining a rendered 3D panel, reading
  a descriptor for one molecule and medium, the reference solvent. A figure never builds
  a path. The docstring states what a change to the contract would mean for existing
  plugins.
- Figure discovery from two sources: the `dftfig.figures` entry-point group, which the
  four built-ins now use themselves and which any pip-installed package can declare, and
  plain `.py` files in a project-local directory (`figures.plugin_dir`, `./dft_figures`
  by default) for the common case where a user does not want to package anything. Both
  appear in `dft params figures.which` and in the interview, with where they came from;
  a project-local file overrides a built-in of the same name.
- `REQUIRES` is checked before a figure is drawn, so a figure needing PyMOL on a machine
  without PyMOL is skipped with that sentence instead of raising from inside the
  renderer. A broken or malformed plugin is reported and skipped rather than killing the
  run, and one declaring a higher `API_VERSION` than this dftfig provides is refused.
- `docs/writing-a-figure.md`: a complete worked example built from scratch - a
  condensed-Fukui bar panel, under 60 lines, no PyMOL - covering the whole loop from
  writing the file to the PNG. The same figure ships runnable in
  `examples/fukui_bars/`, and a test loads it through directory discovery and renders it.
- `figures.plugin_dir` (advanced): where project-local figure modules are looked for.

### Changed
- fig2, fig12, fig13 and fig14 consume only the documented contract. They no longer
  import `_save` out of `fig2`, construct run-directory paths, call `cubes`/`render`/`mlp`
  directly, or read `cfg["execution"]["work_dir"]`. `stage_figures` dispatches through
  the registry instead of four hard-coded `if` branches with four different argument
  lists. The rendered output is unchanged - verified by regenerating all four figures
  from a finished run and comparing them pixel by pixel against the previous
  implementation.
- `PYMOL_RGB` moved from `figures/fig2.py` to `render.py`, next to the palettes it
  belongs to, and the matplotlib colour-bar ramps are now `figstyle.esp_colormap` /
  `figstyle.mlp_colormap`, reachable from any plugin as `ctx.esp_colormap()` /
  `ctx.mlp_colormap()`.

### Fixed
- `FigureContext.descriptor_row` / `property_row` return Python scalars. `df.iloc[0]`
  hands back `numpy.bool_`, for which `value is True` is False - the test Figure 2 uses
  to decide whether to print "minimum confirmed". The caption disappeared the moment the
  row stopped coming from `DataFrame.iterrows`, which is covered by a test now.
- fig12 no longer evaluates `int(permanent_charge or 0)` on a missing value. `NaN` is
  truthy, so that expression reached `int(nan)` and raised `ValueError` for any compound
  whose properties CSV left the column blank.
- fig2's "minimum after N restarts" branch rebound the loop variable `n` - the compound
  name - to the restart count. Harmless only because nothing read it afterwards.

- A pytest suite: 207 tests over the parsers, the config layer, executable discovery,
  the pKa model, input generation under `execution.dry_run` and figure rendering. It
  needs neither ORCA nor a network connection; tests that drive ORCA or PyMOL carry the
  `orca` / `pymol` markers and skip themselves when the program is absent.
- Test fixtures of genuine ORCA 6.1 output (aniline in water, opt+freq / property /
  cation / anion), trimmed to 19-82% of their original size by
  `tools/trim_orca_output.py`. `tests/fixtures/full_file_values.json` records what every
  parser returned for the untrimmed files, and a test compares the trimmed fixtures
  against it, so the trimming is verified rather than trusted.
- The accuracy claim for the SMARTS pKa model is now a parametrised test over 17
  reference compounds with their literature sources, not a sentence in the README.
  Warfarin is an `xfail`: its vinylogous acid is read as a phenol and comes out 5.0 log
  units high.
- GitHub Actions: the suite on Ubuntu and macOS across Python 3.10-3.13, a conda job
  that installs PyMOL so the rendering tests actually run, a ruff lint job, and a build
  job that checks the wheel with twine and runs the installed `dft` command.
- Provenance (`dftfig/provenance.py`), so every number and figure can be traced to the
  calculation that produced it:
  - `descriptors.csv` gains `dftfig_version`, `orca_version`, `orca_keywords_opt`,
    `orca_keywords_sp` and `config_sha`. The ORCA version and both keyword lines are
    read from the output files themselves, never by starting ORCA again.
  - `config_sha` is a 12-character hash of the resolved configuration, excluding what
    cannot change a number: the whole `render` and `figures` sections, plus
    `execution.work_dir`, `execution.scratch_dir`, `execution.orca_path`,
    `execution.nprocs`, `execution.maxcore_mb`, `execution.max_parallel_jobs`,
    `execution.resume` and `cubes.esp_nprocs`. The same chemistry run on a laptop and on
    a cluster hashes identically; a different functional, basis, solvation model,
    population scheme, grid, threshold or molecule set does not.
  - `runs/manifest.json`, written at the end of every run: versions, an ISO-8601 UTC
    timestamp, the platform, the resolved config, the `config_sha`, and the path, size
    and SHA-256 of every input, output, geometry, cube, table, PyMOL script and figure.
  - `dft verify` (`python -m dftfig verify`) re-checks that manifest against the disk and
    reports files that are missing, changed or unlisted. It exits 0 only when everything
    agrees, 1 on any disagreement and 2 when there is no manifest, so it can gate a
    publication script.
  - Figure 2's method stamp carries a second line with the dftfig version and the
    `config_sha`, which is what links a printed figure to a row of `descriptors.csv`.
- `parsers.program_version()` and `parsers.input_keyword_line()`: the ORCA version and the
  literal `!` keyword line, both read out of an existing output file.
- `tests/test_provenance.py` with a real ORCA 6.1.0 output as a fixture
  (`tests/fixtures/ANIL_water_sp.out`): hash stability under core-count changes, hash
  sensitivity to functional and basis, the manifest round-trip through `dft verify`, and
  the version parser.
- `dft doctor`: reports ORCA, PyMOL, physical and logical core counts, memory and
  Python packages, so a new user can diagnose their own machine before asking for help.
- Automatic discovery of ORCA and PyMOL (`dftfig/discover.py`). Checks `$ORCA_PATH`,
  `PATH` and the usual install locations, and remembers what it finds in
  `~/.config/dftfig/paths.yaml`. A candidate is accepted only if `orca_plot` and
  `orca_vpot` sit beside it, which rules out the GNOME Orca screen reader that
  occupies `/usr/bin/orca` on most Linux desktops.
- `LICENSE` (MIT) and `THIRD_PARTY.md` recording every dependency's licence and
  stating explicitly that ORCA is not redistributed.
- macOS and Windows branches for physical-core detection.

### Changed
- `execution.orca_path` now defaults to empty, meaning "search this machine", instead
  of an absolute path inside the original author's home directory.
- Packaging: the version has one source (`dftfig.__version__`); RDKit is declared as the
  dependency it has always been; classifiers, keywords, URLs and test extras added.
  `pip install dftfig` now produces a working installation.
- `install.sh` and `uninstall.sh` no longer contain hardcoded paths. The `dft` launcher
  resolves its own location through symlinks and finds a usable interpreter.
- Process cancellation goes through one `_kill_tree()` helper: process groups on POSIX,
  `taskkill /T` on Windows.
- The arrow-key menu is used only where termios exists; elsewhere every question falls
  back automatically to a numbered prompt.

### Fixed
- The advisory about plotting a LUMO with a diffuse basis checked `cubes.orbitals`
  against `frontier_pm1` and `all_valence`, which are not values the registry offers, so
  it was never issued for `homo_lumo_pm1` or `homo_lumo_pm2` - the two settings that plot
  the most weakly bound virtual orbitals.
- `properties.load_properties` wrote the `formula` column with `.at` before the column
  existed, which pandas 4 turns from a silent insertion into a KeyError.

### Removed
- `run_all.sh`, superseded by `dft run`.

## [1.4.0]

### Fixed
- A failed ORCA job no longer ends the whole run with `OSError: File name too long`.
  `parsers._text()` used `Path(src).exists()` to tell a path from output text; a
  filename component may not exceed 255 bytes, so that call raised rather than
  returning False on a long output, and `error_message()` fed it exactly such a string.
- `execution.nprocs` above the physical core count is clamped on load, and a job that
  still hits Open MPI's "not enough slots available" is retried at the real count.
  Open MPI allocates one slot per physical core, not per hardware thread.

### Changed
- `property.basis_sp` defaults to `def2-TZVP` instead of `6-311++G(d,p)`. A neutral
  organic molecule has a negative electron affinity, so a diffuse basis describes its
  LUMO as a discretised continuum state: measured on benzene, phenol, aniline and
  toluene, `6-311++G(d,p)` gave LUMO energies of +1.11 to +1.21 eV and 9.4-10.8 eV gaps
  with detached lobes, against bound pi* LUMOs and 5.5-6.7 eV without diffuse functions.
  Added `ma-def2-TZVP` and `def2-TZVPD`, and an advisory when a diffuse basis is
  combined with a LUMO cube.
- The molecule table prints once per run rather than once per stage.

## [1.3.0]

### Changed
- Presets removed. `dft run` walks all ten pipeline steps and asks about every setting
  that applies, in execution order. The two figure-styling steps are offered as a single
  keep-or-change question. Parameters that cannot take effect are not asked about.

### Fixed
- `dft` menu "Settings" edited the registry defaults instead of the saved `config.yaml`.

## [1.2.0]

### Fixed
- The MEP is computed with serial `orca_vpot` by default. `orca_vpot_mpi` re-reads
  `<base>.SHARKINP.tmp`, which ORCA deletes on normal termination, so the parallel path
  worked only after a job that had crashed - producing an empty MEP row while the single
  point reported success.
- The saddle-point escape promotes the best attempt, not the last one.
- A mode smaller than `theory.saddle_small_mode_cm` is reported as a free internal
  rotation rather than a reaction saddle; no amount of further optimization removes one.
- Figure 12 uses one depiction scale across the row; Figure 13 captions are anchored to
  their grid cells; Figure 14 de-overlaps labels of compounds with identical TPSA.

## [1.1.1]

### Fixed
- The interactive menu used `tty.setraw`, which also clears `OPOST`, so newlines stopped
  carriage-returning and the display walked diagonally across the terminal. Now
  `tty.setcbreak`.
