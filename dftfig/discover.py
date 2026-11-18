"""Locate the external programs dftfig drives: ORCA and PyMOL.

Neither is a Python dependency. ORCA is proprietary (free for academic use) and the
user installs it themselves; PyMOL usually comes from conda-forge. So the paths cannot
be baked into the package, and a default pointing at one developer's home directory is
worse than no default at all - it fails on every other machine with a message that
blames the user's installation.

Discovery order is deliberate, and verification is not optional:

* On most Linux desktops ``/usr/bin/orca`` is the GNOME Orca SCREEN READER, not the
  quantum chemistry program. A bare ``shutil.which("orca")`` finds it, and the failure
  that follows is baffling. Every candidate is therefore checked for the sibling helper
  binaries (``orca_plot``, ``orca_vpot``) that only a real ORCA installation has.
* A path the user confirmed once is remembered in a user-level file, so the question is
  asked at most once per machine rather than once per project.
"""
from __future__ import annotations

import os
import shutil
from pathlib import Path

# A genuine ORCA installation ships these next to the main binary. The screen reader
# does not, and neither does anything else that happens to be called "orca".
ORCA_SIBLINGS = ("orca_plot", "orca_vpot")

# Where people actually unpack ORCA. Globs, newest-looking last.
ORCA_HINTS = (
    "~/orca*", "~/ORCA*", "~/opt/orca*", "~/software/orca*", "~/bin/orca*",
    "/opt/orca*", "/opt/ORCA*", "/usr/local/orca*", "/usr/local/ORCA*",
    "/share/apps/orca*", "/opt/software/orca*", "/Applications/orca*",
)


def user_config_path() -> Path:
    """Machine-level settings, shared by every project on this computer."""
    base = os.environ.get("XDG_CONFIG_HOME") or "~/.config"
    return Path(base).expanduser() / "dftfig" / "paths.yaml"


def _remembered() -> dict:
    p = user_config_path()
    if not p.is_file():
        return {}
    try:
        import yaml
        return yaml.safe_load(p.read_text()) or {}
    except Exception:
        return {}


def remember(key: str, value: str) -> None:
    """Record a confirmed path so it is never asked for again on this machine."""
    p = user_config_path()
    data = _remembered()
    data[key] = str(value)
    try:
        import yaml
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(yaml.safe_dump(data, sort_keys=True))
    except Exception:
        pass            # a cache that cannot be written is not worth failing a run over


def is_real_orca(candidate) -> bool:
    """True when this binary is the quantum chemistry ORCA.

    Checked by the presence of the helper executables that ship with it, which is both
    cheap and decisive: it rules out the GNOME screen reader without executing anything.
    """
    if not candidate:
        return False
    p = Path(candidate).expanduser()
    if not (p.is_file() and os.access(p, os.X_OK)):
        return False
    return all((p.parent / s).is_file() for s in ORCA_SIBLINGS)


def orca_candidates() -> list[Path]:
    """Every plausible ORCA binary on this machine, verified, best first."""
    seen: list[Path] = []

    def add(c) -> None:
        if not c:
            return
        p = Path(c).expanduser()
        if p.is_dir():
            p = p / "orca"
        if is_real_orca(p) and p.resolve() not in [q.resolve() for q in seen]:
            seen.append(p)

    add(_remembered().get("orca_path"))
    add(os.environ.get("ORCA_PATH"))          # the variable ORCA's own docs use
    add(os.environ.get("ORCA_BIN"))
    add(shutil.which("orca"))                 # may be the screen reader; add() filters it
    for pattern in ORCA_HINTS:
        root = Path(pattern).expanduser()
        try:
            for hit in sorted(root.parent.glob(root.name)):
                add(hit / "orca")
                add(hit)
        except (OSError, IndexError):
            continue
    return seen


def find_orca(explicit: str | None = None) -> Path | None:
    """The ORCA binary to use, or None if this machine has none we can verify."""
    if explicit:
        p = Path(explicit).expanduser()
        if p.is_dir():
            p = p / "orca"
        if is_real_orca(p):
            return p
        if p.is_file():
            return p      # trust an explicit path even unverified; warn elsewhere
    cands = orca_candidates()
    return cands[0] if cands else None


def orca_version(exe) -> str:
    """ORCA's own version string, for provenance. Empty when it cannot be determined."""
    import re
    import subprocess
    import tempfile
    exe = Path(exe)
    try:
        with tempfile.TemporaryDirectory() as d:
            inp = Path(d) / "v.inp"
            inp.write_text("! HF STO-3G\n* xyz 0 1\nH 0 0 0\nH 0 0 0.74\n*\n")
            r = subprocess.run([str(exe), str(inp)], capture_output=True, text=True,
                               timeout=120, cwd=d)
            m = re.search(r"Program Version\s+(\S+)", r.stdout or "")
            return m.group(1) if m else ""
    except Exception:
        return ""


def find_pymol(explicit: str | None = None) -> Path | None:
    """A PyMOL launcher, preferring the one in this interpreter's own environment.

    A system-wide /usr/bin/pymol often starts a different Python that cannot import the
    pymol module, and - unhelpfully - still exits 0 while rendering nothing.
    """
    import sys
    if explicit:
        p = Path(explicit).expanduser()
        if p.is_file():
            return p
    local = Path(sys.executable).parent / "pymol"
    if local.is_file():
        return local
    for c in (_remembered().get("pymol_path"), shutil.which("pymol"),
              shutil.which("pymol.exe")):
        if c and Path(c).is_file():
            return Path(c)
    return None


def report() -> dict:
    """What this machine has, for `dft where` and for the diagnostics in `dft doctor`."""
    orca = find_orca()
    pymol = find_pymol()
    return {
        "orca": str(orca) if orca else "",
        "orca_candidates": [str(c) for c in orca_candidates()],
        "pymol": str(pymol) if pymol else "",
        "user_config": str(user_config_path()),
    }
