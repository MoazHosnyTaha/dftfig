# What dftfig computes, and how

The science behind the numbers: which choices matter, what is derived rather than calculated, what the pipeline checks, and what it cannot do.

<!-- Moved out of README.md to keep it short. Everything here was in the README up to v0.6.0. -->

---

## The parameters that matter most

A short orientation; `python -m dftfig params --key …` has the full text for each.

**`theory.functional` + `theory.dispersion`.** The manuscript used bare B3LYP, which has no London
dispersion — while arguing throughout that the compounds work by π-stacking. If you are
reproducing, keep `dispersion: none`; if you are doing new work, use `D3BJ`. It costs nothing.
Do not set `dispersion` when the functional already contains it (ωB97X-D4, r2SCAN-3c) — validation
will stop you.

**`solvation.model`.** ORCA has **no IEF-PCM**. CPCM is the closest equivalent and agrees with it
to well under 0.1 kcal/mol in water and DMSO; the two diverge more in chloroform (ε = 4.7). SMD
adds the non-electrostatic term and is the better model, at the price of no longer matching the
published numbers exactly.

**`property.koopmans_mode`.** Koopmans' theorem is a Hartree–Fock result. In Kohn–Sham DFT the
HOMO equals the exact ionization energy only for the exact functional, and B3LYP undershoots by
1–3 eV — the manuscript's IE of ~5.5 eV for a carbazole should be read against ~7.6 eV measured.
`both` reports Koopmans and ΔSCF side by side; the ΔSCF values are free if `fukui.enabled`.

**`render.esp_range` and `render.esp_autoscale`.** The most abused setting in published MEP
figures. Auto-scaling each panel to its own range — what most GUIs do by default — means the same
colour encodes a different potential in each panel, and the figure is not a comparison. The
default fixes one range across all molecules and prints it on the colour bar.

**`render.orbital_isovalue` / `render.density_isovalue`.** 0.02 a.u. for orbitals and 0.0004 a.u.
for the density surface are the literature conventions. Whatever you choose, the pipeline stamps
them on the figure, because the manuscript states neither.

**`execution.maxcore_mb`.** Memory **per MPI process**, not in total. Real usage is
`maxcore_mb × nprocs × max_parallel_jobs`. This is the most common cause of ORCA crashes.

**`figures.som_source`.** `fukui_radical` builds `mpo_metabolism` panels B–E from the radical Fukui index this
pipeline already computes — self-contained, but a reactivity proxy that knows nothing about how
the substrate sits in the CYP active site. `csv` reads a trained prediction you export from
XenoSite or ADMETlab. Say which one you used in the caption.

---

### Which basis for the single point

The property single point is where the orbital energies, the populations, the cubes and the
MEP all come from, and one basis has to serve all four. The tension is this:

* **Diffuse functions** (`++` in 6-311++G(d,p), `aug-` in aug-cc-pVTZ, the `D` in def2-TZVPD)
  describe density far from the nuclei. That is what the MEP tail is made of, and what an
  anion needs.
* **A neutral organic molecule has a negative electron affinity.** Its LUMO is not a bound
  state. Give the calculation diffuse functions and that unbound state is described as a
  discretised piece of the continuum: the orbital energy comes out positive, the plotted
  "LUMO" is a vague cloud sitting beside the molecule rather than the π* orbital, and the
  HOMO–LUMO gap is a property of the basis set rather than of the compound.

Measured on benzene, phenol, aniline and toluene:

| basis | LUMO energy | HOMO–LUMO gap | what the LUMO panel shows |
|---|---|---|---|
| 6-311++G(d,p) | +1.11 to +1.21 eV | 9.4–10.8 eV | detached diffuse lobes |
| def2-SVP | −0.09 to −0.24 eV | 5.47–6.73 eV | π* with the correct nodal pattern |

So `property.basis_sp` defaults to **def2-TZVP**: bound π*, stable orbital energies,
well-behaved populations. Use **ma-def2-TZVP** when one single point has to give both a
converged MEP and a sensible LUMO, and a fully diffuse set (**def2-TZVPD**) only when you
are deliberately studying an anion. The pipeline warns when a diffuse basis is combined
with a LUMO cube.

Note that this is about the *property* step only. The geometry optimization basis
(`theory.basis_opt`) has no such constraint — diffuse functions there mostly cost time.

---


## How pKa and logD are obtained

`dftfig/pka.py` holds a curated SMARTS table: about forty ionizable classes, each with a
class-average pKa from the standard compilations (Perrin; Serjeant & Dempsey; the ranges in
Manallack, *Perspect. Medicin. Chem.* 2007, **1**, 25). Matching runs on a neutralized copy
of the structure, so a compound drawn as its hydrochloride still reports the amine.

Checked against measured values, it reproduces the dominant pKa within about one log unit
across a reference set spanning carboxylic acids, phenols, anilines, azines, azoles,
amidines and amino acids — aspirin 4.2 (3.5), ibuprofen 4.5 (4.4), phenol 10.0 (10.0),
propranolol 9.5 (9.5), pyridine 5.2 (5.2), imidazole 7.0 (7.0), benzimidazole 5.5 (5.5),
adenine 5.0 (4.2), metformin 12.4 (12.4), lidocaine 7.8 (7.9), glycine 2.2/10.6 (2.35/9.6).

It gets four things deliberately right that a naive SMARTS list gets wrong:

1. a carbazole, pyrrole or indole N–H is **not** ionizable between pH 0 and 14, so the
   chlocarbazomycins come back with no basic centre at all — the manuscript's basic pKa of
   4.87–5.51 for those compounds has no chemical basis
2. an amino group on a pyridine or pyrimidine does not protonate on the exocyclic nitrogen;
   the ring nitrogen does, and only once per fused ring system
3. a tetrazole's three remaining nitrogens are not basic azole nitrogens
4. a biguanide matches the guanidine pattern twice but has one strongly basic centre

That accuracy claim is a test, not a paragraph: `tests/test_pka.py` carries the reference
table with its literature sources and asserts the dominant pKa to within 1.5 log units for
every compound in it. Run `pytest tests/test_pka.py` to reproduce it.

Known weaknesses are listed in `pka.KNOWN_WEAKNESSES`: vinylogous acids (warfarin) are
reported as ordinary phenols and come out several units too high, strong ring substituent
effects are not modelled, only the first ionization of a polyprotic acid is taken, and
tautomers are used as drawn. Warfarin is in the test file as an `xfail` with that
explanation attached, so the limitation is visible rather than discovered.

logD(7.4) follows from Crippen clogP and the neutral fraction under the independent-site
Henderson–Hasselbalch model. Only the neutral microspecies is taken to partition, which is
standard and which makes the value a **lower bound for a strong base** — ion pairing adds
roughly a log unit in practice. A permanently charged compound has no neutral form at all,
so its logD is reported at a floor and flagged.

Set `figures.pka_source: csv` to use only measured values, or `none` to switch the model
off and draw each compound at the formal charge of the structure you supplied.

---


## Putting a series in one viewpoint

`render.orientation` decides how each molecule is turned before it is rendered, and for a
congeneric series the choice changes what the figure is for.

**`principal_plane`** (the default) orients each molecule on its own inertia tensor. For
flat fused rings that puts the molecular plane in the screen plane, which is usually what
you want. But the tensor is a property of the *whole* molecule, so it drifts as soon as
the substituents differ — add a chain to one analogue and its scaffold tilts relative to
the others.

**`first_molecule`** fits the series into one frame instead. The maximum common
substructure is found across the whole set (RDKit, rings matched whole), the first
molecule's copy of it defines the reference frame from that core's principal axes, and
every other molecule is least-squares fitted onto it (Kabsch). The scaffold then sits in
the same place in every panel and the reader's eye compares the chemistry rather than the
viewpoint. If the set shares nothing large enough to orient by, dftfig says so and each
panel falls back to its own principal axes.

Three details decided here rather than left to a library default:

- **The camera moves, not the molecule.** The orbital and electrostatic-potential panels
  are isosurfaces through a cube, and a cube's values live on a grid fixed in space —
  rotating the molecule would leave its own density behind. The fit is written into
  PyMOL's view matrix.
- **Reflections are rejected.** A mirror image fits near-coplanar points at least as well
  as a rotation does, and a singular value decomposition will happily return one. A
  mirrored molecule in a figure is a drawing of the other enantiomer, so the determinant
  is forced positive.
- **The reference frame has a fixed sign convention.** The sign of an eigenvector is
  arbitrary and so is the handedness of a frame built from three of them, so the largest
  component of each axis is made positive and the third axis is replaced by the cross
  product of the first two. Without that, two runs on the same molecules could disagree.

A wrong rotation here does not look wrong: every panel still shows a plausible molecule
from a plausible angle, which is exactly the failure this setting exists to prevent. The
tests therefore check the geometry numerically — that the cores coincide once the camera
rotation is applied — and check PyMOL's view-matrix convention against PyMOL itself.

---


## Steric and shape descriptors

Everything in this section is computed from element symbols and Cartesian coordinates
alone — the geometry ORCA already optimized. No extra quantum chemistry, no external
program, no second conformer. The electronic half of the descriptor set comes from the
single points at the level of theory the run was configured with; this is the steric half
beside it.

**Whole-molecule, appended to `descriptors.csv`:**

| column | meaning |
|---|---|
| `sasa_area_A2`, `sasa_volume_A3` | solvent-accessible surface and the volume it encloses |
| `disp_area_A2`, `disp_volume_A3` | the same on the dispersion surface |
| `disp_p_int`, `disp_p_min`, `disp_p_max` | dispersion descriptor P<sub>int</sub> and its range over the surface — how dispersion-attractive the molecule is |
| `max_diameter_A` | largest interatomic distance: the longest dimension |
| `pmi1`–`pmi3`, `npr1`, `npr2` | principal moments of inertia and the normalised ratios that place a molecule on the rod–disc–sphere triangle |
| `asphericity`, `eccentricity`, `spherocity_index`, `inertial_shape_factor` | shape indices derived from those moments |
| `radius_of_gyration_A` | mass-weighted spread about the centre of mass |
| `plane_best_fit_A` | mean deviation from the best-fitting plane: how flat the molecule is |

**Per atom, written to `atom_descriptors.csv`** beside the Fukui functions and charges:
`sasa_area_A2`, `disp_area_A2`, `disp_p_int`, `disp_p_min`, `disp_p_max`,
`buried_volume_A3`, `free_volume_A3`, `buried_volume_fraction`, `coordination_number`,
`pyramidalization_P` and `pyramidalization_angle_deg`.

This is the table to hand a QSAR or machine-learning model: one row per atom, the local
electronic response and the local steric environment side by side, with the molecule and
medium named in the first two columns.

Three choices in there are not the library defaults, and are worth knowing about:

- **Buried volume counts hydrogens.** MORFEUS excludes them by default, which is right
  for a phosphine cone around a metal centre and wrong here — on an organic molecule the
  hydrogens are most of what crowds an atom. Excluding them gives methane's carbon a
  buried volume of exactly zero.
- **Pyramidalization is computed from bonds, not proximity.** Asked for the neighbours of
  an atom, MORFEUS takes the three nearest ones, which hands every hydrogen a
  pyramidalization near 2. It is reported here only for atoms with exactly three bonded
  neighbours, from a covalent-radius connectivity matrix, and left empty elsewhere.
- **Fukui columns are filled for the reference medium only.** The condensed Fukui
  functions need the N−1 and N+1 single points, which are run in the reference solvent;
  the steric columns are geometric and exist for every medium that optimized. An empty
  Fukui cell in another medium means "not run there", not "failed".

**These describe one conformer.** The geometry is the lowest-energy conformer from the
force-field search, optimized at the DFT level — not a Boltzmann-weighted average over an
ensemble. Workflows that report Boltzmann-weighted descriptors run a cheap semi-empirical
method over every conformer to get the weights; dftfig does not, because the same
molecules are then carried into DFT single points and figures where one geometry is what
is wanted. If you need ensemble averages, generate the conformers separately and average
the per-conformer tables yourself.

Set `property.sterics: false` to skip the whole block. It needs
[MORFEUS](https://github.com/digital-chemistry-laboratory/morfeus); without it the
electronic descriptors are still produced and the steric columns are simply absent.

---


## Properties the pipeline cannot calculate

pKa, logD, permeability and clearance are not quantum-chemical observables. They come from
`data/properties.csv`, one row per molecule:

| column | used by | source |
|---|---|---|
| `pKa_acidic`, `pKa_basic` | `speciation` (both panels) | ADMETlab, MarvinSketch, your own model — semicolon-separated for several sites |
| `pKa_basic_main` | CNS MPO | the most basic pKa |
| `logP`, `logD74` | `mpo_metabolism` panel A, CNS MPO | ADMETlab 3.0 |
| `TPSA`, `MW`, `HBD` | `lipophilicity`, 14A, CNS MPO | filled from RDKit if left blank |
| `Papp_MDCK`, `Papp_Caco2` | `mpo_metabolism` panel A bubble area | ADMETlab / Deep-PK |
| `CLint_u`, `PPB`, `logBB`, `HIA`, `Pgp_substrate` | carried into `properties_resolved.csv` | ADMETlab / Deep-PK |

The shipped file is pre-filled with the manuscript's published values so the figures reproduce out
of the box. Two cautions about those numbers, both documented in `computational_review.md`:

* The **basic pKa values (2.99–5.51)** are not chemically credible for a carbazole. The
  pyrrole-type N9 is an extremely weak base (conjugate-acid pKa around −6); its N–H is in fact
  weakly *acidic*. A basic pKa near 5 would imply appreciable cationic character at physiological
  pH, which drives both the `speciation` charge curve and one of the six CNS MPO terms. Treat them as
  placeholders and replace them.
* The paper's **`logBB` and `BBB` rows contradict each other** by two orders of magnitude
  (logBB = −1.56 ⇒ C_brain/C_blood = 0.028, against a tabulated "BBB" of 6.91). Only `logBB` is
  carried through here.

Missing values are handled, not faked: a panel that needs a number it does not have says so on the
figure rather than quietly dropping a molecule.

---


## What the pipeline checks for you

* **Every optimized structure is tested for imaginary frequencies, and a saddle point is
  fixed rather than merely reported.** The manuscript states frequencies were computed but
  never reports the result — which is the one thing that calculation was needed for. See
  section 9b.
* **One atom-numbering scheme throughout.** Fukui indices, metabolism maps and optional atom
  labels all key off the RDKit atom order of the SMILES you supplied, which is also the ORCA input
  order. The manuscript uses three incompatible schemes across sections and quotes dihedrals over
  atoms (C10, C12, N11) that do not exist in carbazole numbering.
* **Isovalues and colour ranges are stamped on every figure** and held fixed across panels.
* **The level of theory is stamped on `frontier_orbitals`**, including the two-step `opt // property` protocol
  that the manuscript implies but never describes.
* **Both IE/EA derivations are reported** so the Koopmans approximation is visible rather than
  silent.
* **A positive LUMO energy is flagged.** With a diffuse basis (`6-311++G(d,p)`, `aug-cc-pVTZ`) the
  lowest virtual orbitals of a molecule whose true electron affinity is negative are discretised
  continuum states, not valence antibonding orbitals. The HOMO–LUMO gap is then not an excitation
  energy and the Koopmans EA is meaningless. `descriptors.csv` records the warning, the terminal
  prints it, and the gap annotation on the frontier-orbitals figure turns amber and says so.
* **Koopmans and ΔSCF electron affinities are compared in sign.** When they disagree, the orbital
  picture and the energy-difference picture are telling opposite stories and only ΔSCF is usable.
* **Formal charge is carried into `speciation`.** A permanently charged compound is plotted at its own
  charge, not at zero, and the panel says which charge it is drawn at.
* **MEP colour ranges are measured, not assumed** (`render.esp_scale: auto_shared`). A fixed
  ±0.05 a.u. range saturates for an ion — the potential around a cation is positive everywhere —
  and the panel becomes a featureless blob.

### Escaping a saddle point

An imaginary frequency means the optimizer stopped at a saddle point. The structure is not
a stable species, so the gap, the thermochemistry, the Fukui indices and the MEP computed
from it all describe something that does not exist. `theory.saddle_max_attempts` (default
**2**) makes the pipeline fix it:

1. read the imaginary normal mode out of the ORCA output
2. displace the geometry along it — scaled so the atom that moves most moves
   `theory.saddle_step_angstrom` (0.3 Å)
3. optimize again with `TightOpt`, and check the frequencies
4. if it fell back into the same saddle, try the opposite direction; then a finer grid

An imaginary mode smaller than `theory.saddle_small_mode_cm` (50 cm⁻¹) is usually a methyl
or hydroxyl rotation sitting in integration-grid noise, so that case tries `DefGrid3` with
a tighter optimization *first* and only displaces if the mode survives.

The successful geometry becomes `<name>_<solvent>_opt.xyz`, which is what every later
stage reads; the original saddle-point files are kept as `<name>_<solvent>_opt_saddle.*`
and each restart as `..._r1.*`, so nothing is overwritten silently. The whole history goes
to `runs/saddle_history.csv`, `descriptors.csv` gains a `saddle_restarts` column, and
the frontier-orbitals figure says **"minimum after 1 restart"** instead of a plain "minimum confirmed" — a
geometry that needed rescuing should not pass as one that did not.

Verified on planar ammonia, the textbook case:

| | planar start | after 1 restart |
|---|---|---|
| mean H–N–H | 120.00° | 105.91° (experiment 106.7°) |
| out of plane | 0.00° | 22.83° |
| lowest frequency | **−887.6 cm⁻¹** | none imaginary |
| energy | −56.50023506 Eh | −56.51028566 Eh |

The 6.31 kcal/mol the pipeline climbed down is the NH₃ inversion barrier, measured at
5.8 kcal/mol.

Set `theory.saddle_max_attempts: 0` to get the old behaviour — a loud warning and nothing
else. Each restart is a full OPT+FREQ, so budget for it on a large batch.

**The best attempt is kept, not the last one.** Each restart is scored on the number of
imaginary modes first and on the magnitude of the largest one second, and the winner is the
geometry that gets promoted into the figures and `descriptors.csv`. Without this, a restart
that landed *further* from a minimum than the structure it started from would silently
replace it — in one test run toluene went −31.0 → −27.8 → −52.9 cm⁻¹ and the −52.9 cm⁻¹
geometry was the one published.

**A flat mode is not a saddle point.** A methyl, methoxy or ethyl rotor turns almost freely;
the harmonic approximation assigns it a frequency near zero, and numerical noise in the
integration grid pushes that value to either sign. When the surviving mode is smaller in
magnitude than `theory.saddle_small_mode_cm` (default 50 cm⁻¹), the pipeline stops restarting
and says so explicitly:

```
flat mode -27.8 cm-1 (free internal rotation, not a reaction saddle)
```

the frontier-orbitals figure shows this as an amber badge rather than the red `NOT A MINIMUM`, and
`descriptors.csv` carries a `flat_mode` column. The distinction matters because no amount of
further optimization will remove such a mode — the correct treatment is a hindered-rotor
partition function, which is outside what this pipeline does. The thermochemistry it prints
for such a molecule carries a harmonic-rotor error of roughly 0.1–0.5 kcal/mol in the entropy
term.

---

### Why the MEP is computed on one core

The molecular electrostatic potential comes from `orca_vpot`, which is the only ORCA utility
that produces one — `%plots` has no `ElPot` type and `orca_plot`'s 22 plot types do not
include it either (both verified on ORCA 6.1.0).

`orca_vpot` has an MPI build, `orca_vpot_mpi`, and it is a trap. The MPI binary re-reads
`<basename>.SHARKINP.tmp`, and **ORCA deletes that file when a job terminates normally**. So
the parallel path works only after a crashed job and fails after a successful one — the
failure mode is an empty MEP row in the frontier-orbitals figure with the single-point job reporting success.

`cubes.esp_nprocs` therefore defaults to `1`. Serial evaluation is not the bottleneck people
expect: at wB97X-D4/6-311++G(d,p) with SMD water, a 512 000-point grid on benzene took 36.9 s.
If you raise `esp_nprocs` above 1 the pipeline warns, tries the MPI binary anyway, and falls
back to serial automatically when it fails.

---

### The warnings column

`descriptors.csv` carries a `warnings` column, empty for a clean result. Everything the pipeline
flags on a figure also lands there, so a batch of 40 compounds can be triaged with one `grep`:

```bash
awk -F, 'NR==1 || $NF != ""' runs/descriptors.csv
```

---


## Provenance — tracing a number back to its calculation

A figure in a manuscript is a claim. Six months later, a reviewer asks which functional produced
panel B, or a co-author sends back `descriptors.csv` with one value that looks wrong — and nobody
can tell whether it came from the run that made the figure. Every run therefore records what
produced it, in three places.

**In `descriptors.csv`.** Five columns at the end of every row:

| column | what it is |
|---|---|
| `dftfig_version` | the version of this package that wrote the row |
| `orca_version` | parsed from that molecule's own output banner — not from the ORCA installed today |
| `orca_keywords_opt` | the literal `!` line of the optimization job |
| `orca_keywords_sp` | the literal `!` line of the property single point |
| `config_sha` | short hash of the settings that can change a number |

The two keyword lines are read back out of the ORCA output's own echo of its input, so they are
what ran, not what the config intended: the runner rewrites `%pal` when Open MPI refuses the
requested slot count, and a saddle-point restart appends keywords ([escaping a saddle point](#escaping-a-saddle-point)).

**`config_sha`** answers one question: *are these two results comparable?* It hashes the resolved
configuration — serialized to YAML with sorted keys — after dropping everything that cannot change
a number:

- the whole `render` and `figures` sections (camera, isovalues, colours, DPI, formats, layout);
- `execution.nprocs`, `execution.maxcore_mb`, `execution.max_parallel_jobs`, `cubes.esp_nprocs`
  (parallelism and memory), `execution.work_dir`, `execution.scratch_dir`, `execution.orca_path`
  and `execution.resume`.

So the same chemistry run on a six-core laptop and on a 64-core node gives the **same**
`config_sha`, and changing the functional, a basis set, the solvation model, the population
scheme, the integration grid, an SCF threshold or the molecule set gives a **different** one. If
the hash matches, the two numbers were produced by the same recipe. One caveat worth knowing:
`figures.reference_solvent` is excluded with the rest of its section, so it selects which medium
the Fukui table was built for without changing the hash — but the medium is in the file name and
in the `medium` column.

**In `runs/manifest.json`**, written at the end of every run (including a single re-run stage,
because a stale manifest is worse than none):

```json
{
  "dftfig_version": "1.4.0",
  "orca_version": "6.1.0",
  "generated_utc": "2026-03-04T11:42:07Z",
  "platform": {"system": "Linux-6.8.0-x86_64", "machine": "x86_64", "python": "3.11.8"},
  "config_sha": "20cc2d7027f8",
  "config": { "...the whole resolved configuration..." },
  "files": [
    {"path": "CCB1/water/CCB1_water_sp.out", "bytes": 1840221, "sha256": "9f86d0…"},
    {"path": "figures/frontier_orbitals.png",          "bytes":  612004, "sha256": "e3b0c4…"}
  ]
}
```

Every input, output, geometry, cube, table, PyMOL script and figure is listed with its size and
SHA-256. `.gbw` and `.hess` are not: they are large binary restart files that ORCA does not write
byte-identically twice, so hashing them would report differences that mean nothing.

**On `frontier_orbitals`.** The method stamp along the bottom carries a second line, `dftfig 0.6.0 | config
20cc2d7027f8`. That is enough to match a printed figure to a row in `descriptors.csv` and to the
run directory it came from.

### Checking a run with `dft verify`

```bash
dft verify                     # uses ./config.yaml
dft verify -c path/to/config.yaml
```

It re-hashes everything the manifest lists and reports three kinds of disagreement:

- **missing** — a listed file is gone;
- **changed** — a file's bytes differ from what was recorded (a figure re-saved by hand, a CSV
  edited in a spreadsheet, a corrupted copy);
- **not in the manifest** — a tracked file exists that the manifest does not list, which means the
  manifest is out of date.

It exits `0` only when everything agrees, so it can be used as a gate:

```bash
dft verify && cp runs/figures/frontier_orbitals.pdf ~/manuscript/
```

Any disagreement exits `1`; a missing manifest exits `2`. The usual cause of "changed" is honest —
somebody re-rendered a figure without re-running the stage — and the fix is to re-run the stage so
the manifest describes what is actually on disk.

---


## Citing the methods

Whatever you run, the caption needs: functional, both basis sets, the solvation model **named as
CPCM or SMD rather than IEF-PCM**, the population scheme, the isovalues for every surface, the MEP
colour range, and whether IE/EA are Koopmans or ΔSCF. `descriptors.csv` and the stamp on `frontier_orbitals`
carry all of it.

---

---

[← back to the README](../README.md)
