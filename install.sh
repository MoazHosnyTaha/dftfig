#!/usr/bin/env bash
# Convenience installer for dftfig.
#
#     bash install.sh
#
# It does exactly what the README describes, nothing hidden:
#   1. creates a conda environment from environment.yml (this is how PyMOL and
#      RDKit are obtained; the open-source PyMOL build is not on PyPI)
#   2. pip-installs this checkout into it
#   3. runs `dft doctor` so you can see what was found
#
# ORCA is NOT installed. It is a separate program, free for academic use, which
# you download and install yourself; dftfig locates it at runtime.
#
# If you would rather do it by hand:
#     conda env create -f environment.yml && conda activate dftfig && pip install -e .
set -euo pipefail

HERE="$(cd -P "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_NAME="${DFTFIG_ENV_NAME:-dftfig}"

say()  { printf '  %s\n' "$*"; }
ok()   { printf '  \033[1;32mOK\033[0m  %s\n' "$*"; }
warn() { printf '  \033[1;33m!\033[0m   %s\n' "$*"; }
die()  { printf '  \033[1;31mX\033[0m   %s\n' "$*" >&2; exit 1; }

printf '\n\033[1;35mdftfig\033[0m - installer\n\n'

CONDA=""
for c in micromamba mamba conda; do
  command -v "$c" >/dev/null 2>&1 && { CONDA="$c"; break; }
done

if [ -z "$CONDA" ]; then
  warn "no conda/mamba/micromamba found."
  say  "dftfig itself installs fine with pip:"
  say  "    pip install -e '$HERE'"
  say  "but the 3D panels need PyMOL, whose open-source build is distributed"
  say  "through conda-forge rather than PyPI. Without it the pipeline still runs"
  say  "and simply skips those panels."
  printf '\n  install with pip only? [y/N] '
  read -r reply </dev/tty || reply=n
  case "$reply" in
    [yY]*) python3 -m pip install -e "$HERE" && ok "installed" ;;
    *) die "nothing done" ;;
  esac
else
  ok "using $CONDA"
  if "$CONDA" env list 2>/dev/null | grep -qE "^\s*${ENV_NAME}\s"; then
    warn "environment '$ENV_NAME' already exists - installing into it"
  else
    say "creating environment '$ENV_NAME' from environment.yml ..."
    "$CONDA" env create -n "$ENV_NAME" -f "$HERE/environment.yml" \
      || die "environment creation failed"
    ok "created '$ENV_NAME'"
  fi
  ENVPY="$("$CONDA" run -n "$ENV_NAME" python -c 'import sys; print(sys.executable)')" \
    || die "cannot run python inside '$ENV_NAME'"
  "$ENVPY" -m pip install -e "$HERE" || die "pip install failed"
  ok "dftfig installed into '$ENV_NAME'"
fi

printf '\n'
if command -v dft >/dev/null 2>&1; then
  ok "dft resolves to: $(command -v dft)"
else
  warn "'dft' is not on your PATH in this shell."
  say  "activate the environment first:   conda activate $ENV_NAME"
  say  "or run it directly:               $HERE/dft"
fi

printf '\n  checking what this machine provides:\n\n'
(command -v dft >/dev/null 2>&1 && dft doctor) || "$HERE/dft" doctor || true

cat <<'MSG'

  Getting started
  ---------------
    mkdir -p ~/my_project && cd ~/my_project
    dft init        set up this folder
    dft run         walk every setting, then run the pipeline
    dft doctor      check ORCA, PyMOL, cores and packages
    dft --help      all commands

  ORCA is not included. Download it (free for academic use) from the ORCA forum,
  then either put it on your PATH or set $ORCA_PATH. `dft doctor` will confirm.

MSG
