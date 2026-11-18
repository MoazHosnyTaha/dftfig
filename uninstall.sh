#!/usr/bin/env bash
# Remove dftfig. Nothing is deleted without telling you what and where first.
set -euo pipefail
ENV_NAME="${DFTFIG_ENV_NAME:-dftfig}"
CFG="${XDG_CONFIG_HOME:-$HOME/.config}/dftfig"

cat <<MSG

  To remove dftfig:

    pip uninstall dftfig

  If you installed with install.sh, the conda environment can go too:

    conda env remove -n $ENV_NAME

  Remembered executable paths (safe to delete):

    $CFG

  Your project folders - config.yaml, molecules/, runs/ - are NOT touched by
  any of the above. Delete those yourself if you want them gone.

MSG
