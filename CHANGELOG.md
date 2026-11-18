# Changelog

All notable changes to this project are documented here.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project uses [semantic versioning](https://semver.org/).

## [Unreleased]

### Added
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
