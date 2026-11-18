# Third-party software

dftfig is MIT licensed. It depends on the following, none of which imposes
copyleft obligations on this package.

## Imported into the dftfig process

| package | licence |
|---|---|
| RDKit | BSD-3-Clause |
| NumPy | BSD-3-Clause (with 0BSD, MIT, Zlib, CC0 components) |
| SciPy | BSD-3-Clause |
| pandas | BSD-3-Clause |
| matplotlib | matplotlib licence (BSD-compatible, PSF-derived) |
| PyYAML | MIT |
| Pillow | MIT-CMU (HPND) |
| rich | MIT |

## Invoked as external programs

These are run as separate processes. dftfig never imports them, links against
them, or distributes any part of them.

| program | licence | how it is obtained |
|---|---|---|
| ORCA | proprietary; free for academic use under Max Planck Institute terms | the user downloads and installs it themselves from the ORCA forum |
| PyMOL (open-source) | permissive Schrodinger, LLC open-source licence | `conda install -c conda-forge pymol-open-source` |

**ORCA is not redistributed with dftfig and is not a Python dependency.** You
must obtain your own copy and accept its licence. dftfig locates the binary at
runtime and calls it; it contains no ORCA code.
