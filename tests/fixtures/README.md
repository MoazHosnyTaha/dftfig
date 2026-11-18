# Test fixtures

## Real ORCA output

`ANIL_water_opt.out`, `ANIL_water_sp.out`, `ANIL_water_cat.out` and
`ANIL_water_ani.out` are genuine **ORCA 6.1.0** output for aniline (C6H7N) in water,
produced by dftfig itself at

    B3LYP-D3BJ/def2-SVP  Opt + Freq  //  B3LYP-D3BJ/def2-TZVP  SP,  CPCM(Water)

with the cation (N-1 electrons) and anion (N+1) single points at the optimized
geometry. `ANIL.sdf` is the input structure and `ANIL_start.xyz` the ETKDG starting
geometry for the same run; `descriptors.csv` is a real descriptors table from a two
compound run (aniline and phenol).

The four outputs are **trimmed**, not synthesised. A full optimization+frequency output
is 392 kB, almost all of it SCF iteration tables and per-cycle internal coordinates that
no parser reads. `tools/trim_orca_output.py` removes whole regions between marker lines
and leaves every section the parsers look at:

    python tools/trim_orca_output.py FULL.out tests/fixtures/TRIMMED.out --check

`--check` re-parses both files and refuses to write a fixture whose parsed values
differ. `full_file_values.json` records what every parser returned for the ORIGINAL
files, with their byte sizes, and
`test_parsers.py::test_trimmed_fixture_parses_identically_to_the_full_output` compares
the trimmed fixtures against it on every run - so the trimming is verified by the suite
rather than trusted.

| fixture | original | trimmed |
| --- | --- | --- |
| `ANIL_water_opt.out` | 392 kB | 74 kB |
| `ANIL_water_sp.out`  | 69 kB  | 55 kB |
| `ANIL_water_cat.out` | 80 kB  | 65 kB |
| `ANIL_water_ani.out` | 80 kB  | 65 kB |

## Synthesised failure modes

These four are written by hand, because a failing run is not something to produce on
purpose and the failure text is short enough to reproduce exactly. Each reproduces the
message ORCA or Open MPI actually prints:

| fixture | what it reproduces |
| --- | --- |
| `fail_scf_not_converged.out` | `SCF NOT CONVERGED AFTER 125 CYCLES` + `SERIOUS PROBLEM IN SOSCF` |
| `fail_mpi_slots.out` | Open MPI's "not enough slots available" startup failure - the one caused by `execution.nprocs` above the physical core count |
| `fail_opt_not_converged.out` | an optimization that hit the cycle limit |
| `imaginary_mode.out` | a genuine -312.45 cm^-1 mode next to a -4.21 cm^-1 grid artefact, which is what `theory.imaginary_cutoff` exists to separate |

Total fixture size is kept under 400 kB; run `du -sh tests/fixtures` after adding
anything.
