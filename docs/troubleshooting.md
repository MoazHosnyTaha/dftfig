# Troubleshooting

Symptoms, causes and fixes. Run `dft doctor` first — it reports your ORCA and PyMOL paths, core counts and installed packages, which is where most diagnoses start.

<!-- Moved out of README.md to keep it short. Everything here was in the README up to v0.6.0. -->

---

## Troubleshooting

| symptom | cause and fix |
|---|---|
| `ORCA not found` | run `dft doctor`; install ORCA, put it on your `PATH` or set `$ORCA_PATH`, or set `execution.orca_path` |
| `topology discovery failed` | Open MPI cannot read the CPU topology; dftfig retries the job serially by itself. Set `execution.nprocs: 1` to skip the false start |
| `There are not enough slots available in the system` | `execution.nprocs` is above the **physical** core count. Open MPI allocates one slot per physical core, not per hardware thread, so a 6-core/12-thread laptop offers 6. The pipeline now clamps `nprocs` to the physical count when the config is loaded, and retries any job that still hits this at the real count |
| `Unknown identifier in PLOTS block` | an ORCA version older than 6; `%plots` syntax differs |
| job dies with `bad_alloc` / out of memory | lower `execution.maxcore_mb`; remember it is per process |
| `PyMOL executable not found` | `conda install -c conda-forge pymol-open-source`, or set `render.pymol_path` |
| MEP panel blank | `cubes.esp` was false, or `orca_vpot` failed — check for `esp.cube`; if you had set `cubes.esp_nprocs` above 1, set it back to 1 (see [why the MEP is computed on one core](methods.md#why-the-mep-is-computed-on-one-core)) |
| LUMO positive, lobes floating off the molecule | a diffuse basis (`++`, `aug-`, `-TZVPD`) has no bound LUMO for a neutral molecule — see [which basis for the single point](methods.md#which-basis-for-the-single-point). `property.basis_sp` now defaults to `def2-TZVP` for this reason |
| a failed job ends the whole run with `OSError: File name too long` | fixed in 1.4.0 — `parsers._text` treated a long output string as a filename. Update the package |
| `lipophilicity` panel missing, "atom ORDER differs" | the ORCA output's atom order does not match the RDKit molecule, so the Crippen logP values cannot be mapped; the message names the first differing position |
| optimization will not converge | raise `conformers.n_confs`; try `TightSCF`; inspect `*_trj.xyz` in Avogadro |
| small imaginary modes (−5 to −15 cm⁻¹) | usually grid noise in a methoxy rotation — try `grid: DefGrid3` and `opt_level: TightOpt` |
| genuine imaginary mode (< −50 cm⁻¹) | a real saddle point: distort along the mode and re-optimize |
| figure panels misaligned | `render.orientation: principal_plane` for planar systems; check each `render_*.png` individually |
| Fukui values all zero | the N±1 jobs did not finish — `python -m dftfig status` |

---

---

[← back to the README](../README.md)
