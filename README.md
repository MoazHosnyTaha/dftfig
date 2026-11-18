# dftfig — automated DFT and figure pipeline for ORCA

[![tests](https://github.com/moazhassanein/dftfig/actions/workflows/tests.yml/badge.svg)](https://github.com/moazhassanein/dftfig/actions/workflows/tests.yml)
[![build](https://github.com/moazhassanein/dftfig/actions/workflows/build.yml/badge.svg)](https://github.com/moazhassanein/dftfig/actions/workflows/build.yml)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%E2%80%933.13-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

Give it structures; it gives you back the figures. dftfig drives ORCA through geometry
optimization, frequencies, a property single point and the N±1 single points, then renders
the result as publication-ready figures — frontier orbitals, electrostatic potential,
lipophilicity surfaces and physicochemical space — with every parameter explained as you
choose it.

It exists because the gap between "the calculation finished" and "the figure is in the
manuscript" is usually a week of manual Avogadro, GaussView, PyMOL and Illustrator work
that nobody records, and that nobody can reproduce afterwards.

```
 .sdf  ->  charge + multiplicity  ->  OPT+FREQ  ->  property single point  ->  N±1 single points
               (auto, per file)          |                   |                        |
                                   minimum check      orbitals, MEP,            condensed Fukui
                                                      populations               functions, ΔSCF
                                                           |
          Figure 2        Figure 12        Figure 13        Figure 14
      geometry/HOMO/      pKa + net       lipophilicity    CNS MPO space +
        LUMO/MEP         charge vs pH        surfaces      metabolism maps
```

## Structures are the only required input

Drop one `.sdf` per compound into `molecules/` and type `dft run`. Everything else is worked
out for you:

| quantity | where it comes from |
|---|---|
| charge, multiplicity | the formal charges in each file — neutrals, cations and anions run together in one batch |
| starting geometry | the file's 3D coordinates if it has them, otherwise a conformer search |
| MW, formula, TPSA, HBD, clogP | the connection table (RDKit) |
| **pKa (acidic and basic)** | the SMARTS ionization model in `dftfig/pka.py` — about ±1 log unit |
| **logD at pH 7.4** | clogP and that ionization, neutral species only |
| CNS MPO score | Wager's six desirability functions, from the above |
| HOMO/LUMO, gap, IE, EA, Fukui | the ORCA calculations |

The only things it will *not* invent are permeability, clearance, plasma protein binding,
logBB and P-gp: those need a trained ADMET model. Leave them out and the figures say so;
supply them in `data/properties.csv` and they are used. Anything you put in that file
overrides what was derived.

---

## 0. Installing

### Requirements

| | |
|---|---|
| Python | 3.10 or newer |
| [ORCA](https://orcaforum.kofo.mpg.de) | 6.x. Free for academic use; **you install it yourself** — it is not redistributed here and is not a Python dependency |
| PyMOL | optional, for the 3D panels. The open-source build comes from conda-forge |

Linux and macOS are fully supported. Windows is expected to work but is not part of
the tested matrix: the arrow-key menus need termios, so there every question falls back
automatically to a numbered prompt that takes the same answers. Reports and fixes
welcome.

### With conda (recommended — this is how you get PyMOL)

```bash
git clone https://github.com/moazhassanein/dftfig && cd dftfig
conda env create -f environment.yml
conda activate dftfig
pip install -e .
```

### With pip only

```bash
pip install dftfig
```

RDKit and everything else comes with it. PyMOL does not — the open-source build is not on
PyPI — so the 3D panels of Figure 2 and Figure 13 are skipped with a message rather than
failing. Add it with `conda install -c conda-forge pymol-open-source` when you want them.

### Then check what your machine provides

```bash
dft doctor
```

```
  OK  ORCA: /opt/orca_6_1_0/orca
      version 6.1.0
  OK  PyMOL: /home/you/miniconda3/envs/dftfig/bin/pymol
  OK  python packages: all present
```

**dftfig finds ORCA by itself.** It checks `$ORCA_PATH`, then your `PATH`, then the usual
install locations (`~/orca*`, `/opt/orca*`, `/usr/local/orca*`, `/share/apps/orca*`), and
remembers what it finds in `~/.config/dftfig/paths.yaml`. A candidate only counts if
`orca_plot` and `orca_vpot` sit beside it — on most Linux desktops `/usr/bin/orca` is the
GNOME **screen reader**, and running that instead produces a failure that makes no sense.
Set `execution.orca_path` explicitly if you keep several ORCA versions.

### Everyday use

```bash
mkdir -p ~/my_project && cd ~/my_project
dft init          # config.yaml + editable data/ templates in this folder
dft               # interactive menu
dft run           # walk every setting, then run the whole pipeline
dft run opt sp    # selected stages only, no questions
dft figures       # rebuild the figures without redoing DFT
dft status        # what is done, what is missing, what is not a minimum
dft verify        # re-check every file against runs/manifest.json (§3b)
dft validate      # check settings against this machine's cores and RAM
dft doctor        # ORCA, PyMOL, cores, memory, packages
dft params        # the full parameter reference
dft params theory.functional     # one parameter
dft --help
```

`dft` with no arguments opens a menu; every menu entry maps to a subcommand, so anything
you can pick you can also script.

**Config resolution**: `dft` uses `config.yaml` in the current directory. Override with
`dft -c other.yaml ...`, or set `DFTFIG_CONFIG`. One folder per project keeps runs separate.

**Running from a checkout without installing**: the `dft` script at the repository root
works out where it lives and which interpreter can import the dependencies. Symlink it
anywhere, or override with `DFTFIG_HOME` and `DFTFIG_PYTHON`.

---

## 1. One command

```bash
dft run
```

That is the whole interface. It sets the folder up if it is empty, reads your structures,
and then walks you through the run:

1. which compounds to run (all of them, by default — space toggles)
2. whether to include the advanced settings (no, by default)
3. **every setting that applies, one step at a time**, in the order the pipeline uses them
4. start the run? — shown with the job count, disk and memory it will need

There are no presets. Nothing is chosen for you behind a label like *recommended* or
*fast*, and nothing reproduces anybody else's published settings. The walk starts from
the registry defaults and puts each decision to you with what it does, how to choose it
and what it costs. The value already in use is pre-selected in every menu, so **Enter
keeps it** and you can go as fast or as slowly as you like.

Every choice is numbered. Move with the **arrow keys** and pick with **enter**, or just
**type the number** — on a single-choice menu the number picks it immediately, on the
compound list it toggles that compound. `j`/`k` work as up/down, **q** or **Esc** cancels.

If the menus misbehave in your terminal, set `DFTFIG_NO_TUI=1` and every question becomes
a plain numbered prompt with no cursor control at all:

```bash
DFTFIG_NO_TUI=1 dft run
```

Piped input falls back to the same plain prompts automatically, so `dft run < answers.txt`
works in a script.

For scripted use, naming stages skips every question:

```bash
dft run figures          # re-render only; no DFT, no prompts
dft run opt sp           # just those two stages
```

### The ten steps

The walk follows the pipeline in execution order, so each answer arrives when the
calculation is about to use it. Every step opens with what it does and how many settings
it holds:

| step | what it decides | settings |
|---|---|---|
| 1. Molecules | what each compound *is*: charge, spin, whether the file's coordinates are usable | 3 |
| 2. Conformer search | the starting geometry, when the file has no usable 3D coordinates | 4 |
| 3. Level of theory (geometry optimization) | functional, basis, dispersion, RI, SCF and geometry convergence, grid, frequencies, saddle escape | 9 |
| 4. Implicit solvation | CPCM or SMD, and which media | 2 |
| 5. Property single point | the larger basis, the population scheme, Koopmans vs ΔSCF | 4 |
| 6. Fukui functions | whether to run the N±1 points and on which geometry | 2 |
| 7. Volumetric data | the grids Figure 2 is drawn from | 4 |
| 8. 3D rendering | isovalues, colours, transparency, camera | keep, or 18 |
| 9. Figure assembly | which figures, DPI, labels, palettes, MPO variant | keep, or 17 |
| 10. Execution | cores, memory, scratch, resume | 8 |

Steps 8 and 9 are the figure *look* rather than the science, so they are offered as one
keep-or-change question with the current values listed. Say *keep* and the figures come
out exactly as they do now; say *change* and you get all 35 settings individually.

A setting that cannot take effect is never asked about: a fixed charge when charge
detection is automatic, the saddle-escape step size when frequencies are off, anything
under `cubes` when cubes are disabled.

For each parameter it prints:

* **what it does**, in plain language;
* **how to choose** — and when the default is the wrong answer;
* **cost** — the effect on wall time, memory or disk;
* for categorical parameters, a **table of every allowed value** with its own explanation,
  and a `suggested` or `cheap` badge where one value is the obvious starting point.

Advanced parameters (scratch paths, tolerances, iteration limits) are hidden unless you
ask for them in step 0.

At the end it prints the job count, disk and peak-RAM estimate, and **the literal ORCA keyword
line** it will write — so you see the level of theory before committing to it.

### Reading the reference without being prompted

```bash
python -m dftfig params                               # the whole registry
python -m dftfig params --key theory.functional       # one parameter
python -m dftfig params --key render.esp_range
python -m dftfig params --no-advanced                 # skip the fiddly ones
```

### Editing `config.yaml` by hand

Perfectly fine — it is plain YAML and every key maps to a registry entry. Validate before running:

```bash
python -m dftfig validate -c config.yaml
```

Validation catches the combinations that otherwise waste an afternoon: a functional that already
contains dispersion with `dispersion` set as well; a reference solvent you never calculated;
`RIJK` without a `/JK` auxiliary basis; frequencies on top of a `LooseOpt`; Mulliken populations
with a diffuse basis; Figure 2 requested with cubes disabled.

---

## 2. Stages

`dft run` executes all of them in order. You can also run one at a time — useful because the figures are quick to iterate on once the DFT is done:

```bash
python -m dftfig build       -c config.yaml   # SMILES -> starting geometries
python -m dftfig opt         -c config.yaml   # OPT (+FREQ), minimum check
python -m dftfig sp          -c config.yaml   # property single point, cubes
python -m dftfig fukui       -c config.yaml   # N-1 / N+1 single points
python -m dftfig cubes       -c config.yaml   # MEP via orca_vpot
python -m dftfig descriptors -c config.yaml   # descriptors.csv, fukui_*.csv
python -m dftfig figures     -c config.yaml   # the four figures
python -m dftfig status      -c config.yaml   # what is done, what is missing
```

`execution.resume: true` (the default) means a stage re-run skips anything that already
terminated normally, so you can restyle a figure without repeating hours of DFT.

**Before a long campaign**, set `execution.dry_run: true` and run `all`. Every input file is
written and nothing is launched; read the `.inp` files, confirm the keyword line, then turn it off.

---

## 3. What comes out

```
runs/
  CCB1/
    CCB1_start.xyz                  conformer-search output
    water/
      CCB1_water_opt.inp/.out/.xyz  optimization + frequencies
      CCB1_water_sp.inp/.out        property single point
      CCB1_water_cat.out            N-1
      CCB1_water_ani.out            N+1
      homo.cube  lumo.cube          frontier orbitals
      eldens.cube  esp.cube         density and electrostatic potential
      render_*.png  render_*.pml    panel images and the PyMOL script that drew each one
      CCB1_mlp.pdb                  per-atom lipophilicity in the B-factor column
    fukui_CCB1_water.csv            f+, f-, f0, Δf, Δω per atom
  descriptors.csv                   one row per molecule × medium
  manifest_{opt,sp,fukui}.json      per-job status, wall time, energies
  manifest.json                     provenance: versions, settings hash, SHA-256 per file
  figures/
    figure2.{png,pdf}  figure12.*  figure13.*  figure14.*
    properties_resolved.csv         what the figures actually used
```

Every `.cube` opens directly in **Avogadro** or **GaussView** if you want to look around
interactively, and every `.pml` can be opened in PyMOL (`pymol render_mep_CCB1.pml`) to adjust a
camera angle by hand and re-save.

---

## 3b. Provenance — tracing a number back to its calculation

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
requested slot count, and a saddle-point restart appends keywords (§8b).

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
    {"path": "figures/figure2.png",          "bytes":  612004, "sha256": "e3b0c4…"}
  ]
}
```

Every input, output, geometry, cube, table, PyMOL script and figure is listed with its size and
SHA-256. `.gbw` and `.hess` are not: they are large binary restart files that ORCA does not write
byte-identically twice, so hashing them would report differences that mean nothing.

**On Figure 2.** The method stamp along the bottom carries a second line, `dftfig 1.4.0 | config
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
dft verify && cp runs/figures/figure2.pdf ~/manuscript/
```

Any disagreement exits `1`; a missing manifest exits `2`. The usual cause of "changed" is honest —
somebody re-rendered a figure without re-running the stage — and the fix is to re-run the stage so
the manifest describes what is actually on disk.

---

## 4. The parameters that matter most

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

**`figures.som_source`.** `fukui_radical` builds Figure 14 B–E from the radical Fukui index this
pipeline already computes — self-contained, but a reactivity proxy that knows nothing about how
the substrate sits in the CYP active site. `csv` reads a trained prediction you export from
XenoSite or ADMETlab. Say which one you used in the caption.

---

### 4b. Which basis for the single point

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

## 5. Using your own molecules

### The normal route: one SDF per compound

```bash
mkdir -p ~/my_project/molecules && cd ~/my_project
dft init                      # scaffolds config.yaml, molecules/ and data/
cp /path/to/*.sdf molecules/  # one file per compound
dft validate                  # see exactly what was detected
dft run
```

The file name becomes the molecule name, which becomes the output directory and the figure
panel label — so keep it short. A multi-record SDF also works: each record becomes its own
molecule.

**There are no built-in compounds.** `molecules.definitions` is empty by default. To reproduce
the manuscript's chlocarbazomycins, put their SDF files in `molecules/` like any
other compounds — there is no built-in copy of the published settings.

### Charge is detected per compound

This matters more than it sounds. A batch containing a neutral nucleoside and a quaternary
ammonium cation has two different electron counts; running both at one global charge computes a
different chemical species for one of them — wrong SCF, wrong orbitals, meaningless descriptors.

With `molecules.charge_mode: auto` (the default) each compound's charge is read from its own
formal charges, and they are printed before any job starts:

```
  molecule     formula       atoms   charge   mult   electrons   geometry    source
  ACET_ANION   C2H3O2-           7       -1      1          32   from file   ACET_ANION.sdf
  FTU          C11H13FN4O4      33        0      1         148   from file   FTU.sdf
  R016         C19H25N4+        48       +1      1         166   from file   R016.sdf
  charge: detected from formal charges
  mixed charges in this run: -1, 0, +1  - each molecule keeps its own
```

The same table is written to `runs/molecules.csv`. **Check it before a long run** — an unexpected
charge there is the cheapest error you will ever catch.

Compounds with different charges run concurrently; `max_parallel_jobs` schedules across the whole
set regardless of charge state.

Multiplicity follows the same logic: `auto` picks the lowest spin state consistent with the
electron count (singlet if even, doublet if odd). A multiplicity that is impossible for the
electron count is rejected up front rather than failing inside ORCA. Open-shell species get UKS
automatically; their frontier-orbital cubes are skipped with a warning, because alpha and beta
HOMOs differ and a single pair is not well defined.

The N±1 Fukui calculations are taken relative to each molecule's own charge — a cation's "anion"
is its neutral form, not an anion:

| molecule | charge | N−1 | N+1 |
|---|---|---|---|
| acetate | −1 | 0 | −2 |
| water | 0 | +1 | −1 |
| ammonium | +1 | +2 | 0 |

`dft validate` warns if an already-anionic compound is about to gain another electron on a basis
with no diffuse functions — that electron is usually unbound, and you get a negative electron
affinity describing a continuum state rather than a real anion.

### How many compounds to run

```yaml
molecules:
  limit: 3        # 0 = all of them
```

Files are taken in sorted order, so the selection is reproducible, and the names that were
selected are printed. Use it to time one compound before committing to forty.

### 3D geometries

If an SDF already contains 3D coordinates they are used directly as the DFT starting geometry —
a docking pose or a crystal structure encodes information a conformer generator cannot recover.
Records with no 3D coordinates get an ETKDG conformer. Set `molecules.use_sdf_geometry: false` to
force a fresh conformer search for everything.

### SMILES instead of SDF

```yaml
molecules:
  structure_source: smiles
  definitions:
    MyCompound1: "COc1ccc2[nH]c3ccccc3c2c1Cl"
    MyCompound2: "[NH+]1(CCCC1)CC=CCNc2ccnc3c2c4c(CCC4)cn3"
```

Charge detection works identically — the `[NH+]` above is read as +1.

### Existing .xyz geometries

```yaml
molecules:
  structure_source: xyz
  xyz_dir: input_xyz          # expects input_xyz/MyCompound1.xyz
```

A SMILES is still required per molecule: it supplies the bond orders, hydrogen-bond perception and
Crippen lipophilicity increments that an XYZ cannot. The pipeline warns if the atom count
disagrees, and refuses to colour a lipophilicity surface it cannot map atom-for-atom.

---

## 5b. How pKa and logD are obtained

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

## 6. Properties the pipeline cannot calculate

pKa, logD, permeability and clearance are not quantum-chemical observables. They come from
`data/properties.csv`, one row per molecule:

| column | used by | source |
|---|---|---|
| `pKa_acidic`, `pKa_basic` | Figure 12 (both panels) | ADMETlab, MarvinSketch, your own model — semicolon-separated for several sites |
| `pKa_basic_main` | CNS MPO | the most basic pKa |
| `logP`, `logD74` | Figure 14A, CNS MPO | ADMETlab 3.0 |
| `TPSA`, `MW`, `HBD` | Figure 13, 14A, CNS MPO | filled from RDKit if left blank |
| `Papp_MDCK`, `Papp_Caco2` | Figure 14A bubble area | ADMETlab / Deep-PK |
| `CLint_u`, `PPB`, `logBB`, `HIA`, `Pgp_substrate` | carried into `properties_resolved.csv` | ADMETlab / Deep-PK |

The shipped file is pre-filled with the manuscript's published values so the figures reproduce out
of the box. Two cautions about those numbers, both documented in `computational_review.md`:

* The **basic pKa values (2.99–5.51)** are not chemically credible for a carbazole. The
  pyrrole-type N9 is an extremely weak base (conjugate-acid pKa around −6); its N–H is in fact
  weakly *acidic*. A basic pKa near 5 would imply appreciable cationic character at physiological
  pH, which drives both the Figure 12B curve and one of the six CNS MPO terms. Treat them as
  placeholders and replace them.
* The paper's **`logBB` and `BBB` rows contradict each other** by two orders of magnitude
  (logBB = −1.56 ⇒ C_brain/C_blood = 0.028, against a tabulated "BBB" of 6.91). Only `logBB` is
  carried through here.

Missing values are handled, not faked: a panel that needs a number it does not have says so on the
figure rather than quietly dropping a molecule.

---

## 7. Wall-clock expectations

Measured on a 25-atom carbazole, **serial** (one core):

| stage | def2-SVP | 6-31+G(d,p) |
|---|---|---|
| OPT (LooseOpt, no FREQ) | ~6 min | ~20–35 min |
| OPT + analytic FREQ | ~15 min | ~60–100 min |
| property single point | <1 min | ~3–6 min |
| each N±1 single point | <1 min | ~3–6 min |
| MEP (`orca_vpot`, 80³ grid) | ~1–3 min | ~5–10 min |
| rendering + figure assembly | ~1 min | ~1 min |

Multiply by molecules × media. Parallel scaling notes:

* `nprocs` is per job. Hybrid DFT with RIJCOSX scales usefully to about 8–16 cores for a molecule
  this size; past that, communication eats the gain.
* **Running several molecules concurrently usually beats one molecule on all cores** — set
  `max_parallel_jobs: 2` with `nprocs: 4` rather than `max_parallel_jobs: 1` with `nprocs: 12`.
* Keep `max_parallel_jobs × nprocs` at or below your core count, and watch that
  `maxcore_mb × nprocs × max_parallel_jobs` fits in RAM. **`dft validate` now checks both against
  the cores and RAM this machine actually reports** and tells you what to lower — the defaults
  (4 × 2 × 900 MB = 7.0 GB) are sized for a 12-core, 11.5 GB laptop.

### If Open MPI will not start

A message like `topology discovery failed --> Returned value Not supported (-8)` means Open MPI
cannot read the CPU topology — common inside containers and sandboxes where `/sys/devices/system/cpu`
is not mounted. The pipeline **detects this and automatically retries the job serially**, with a
warning. On a normal desktop install this does not occur; if it does, run with `nprocs: 1` and use
`max_parallel_jobs` for throughput instead.

---

## 8. What the pipeline checks for you

* **Every optimized structure is tested for imaginary frequencies, and a saddle point is
  fixed rather than merely reported.** The manuscript states frequencies were computed but
  never reports the result — which is the one thing that calculation was needed for. See
  section 9b.
* **One atom-numbering scheme throughout.** Fukui indices, metabolism maps and optional atom
  labels all key off the RDKit atom order of the SMILES you supplied, which is also the ORCA input
  order. The manuscript uses three incompatible schemes across sections and quotes dihedrals over
  atoms (C10, C12, N11) that do not exist in carbazole numbering.
* **Isovalues and colour ranges are stamped on every figure** and held fixed across panels.
* **The level of theory is stamped on Figure 2**, including the two-step `opt // property` protocol
  that the manuscript implies but never describes.
* **Both IE/EA derivations are reported** so the Koopmans approximation is visible rather than
  silent.
* **A positive LUMO energy is flagged.** With a diffuse basis (`6-311++G(d,p)`, `aug-cc-pVTZ`) the
  lowest virtual orbitals of a molecule whose true electron affinity is negative are discretised
  continuum states, not valence antibonding orbitals. The HOMO–LUMO gap is then not an excitation
  energy and the Koopmans EA is meaningless. `descriptors.csv` records the warning, the terminal
  prints it, and the gap annotation on Figure 2 turns amber and says so.
* **Koopmans and ΔSCF electron affinities are compared in sign.** When they disagree, the orbital
  picture and the energy-difference picture are telling opposite stories and only ΔSCF is usable.
* **Formal charge is carried into Figure 12.** A permanently charged compound is plotted at its own
  charge, not at zero, and the panel says which charge it is drawn at.
* **MEP colour ranges are measured, not assumed** (`render.esp_scale: auto_shared`). A fixed
  ±0.05 a.u. range saturates for an ion — the potential around a cation is positive everywhere —
  and the panel becomes a featureless blob.

### 8b. Escaping a saddle point

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
Figure 2 says **"minimum after 1 restart"** instead of a plain "minimum confirmed" — a
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

Figure 2 shows this as an amber badge rather than the red `NOT A MINIMUM`, and
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
failure mode is an empty MEP row in Figure 2 with the single-point job reporting success.

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

## 9. Interrupting a run

`Ctrl+C` is safe and is the intended way to stop. The pipeline terminates every running
ORCA process (each job runs in its own process group, so helper processes go too), reports
how many jobs had finished, and exits with status 130 — no traceback.

```
  !   terminated 2 running ORCA process(es)
  !   interrupted - 3 of 16 job(s) finished and are saved
  nothing is lost: completed jobs are skipped on the next run (execution.resume),
  and the interrupted one restarts from its input.
  resume with: dft run   check progress with: dft status
```

Completed jobs are detected by ORCA's own `TERMINATED NORMALLY` marker, so a half-written
output from an interrupted job is correctly re-run rather than mistaken for a result.

**A long job is not a hung job.** ORCA writes nothing to the console until a job ends, and an
optimization plus frequencies on a drug-sized molecule takes an hour or more. The live display
shows a spinner and elapsed time per running job plus an overall counter, so you can always
tell the difference:

```
  2 job(s), 2 at a time, 4 core(s) each
  ORCA prints nothing until a job ends; the timer below is live.
  ⠏ FTU (q+0 m1) / water / optimization     0:12:47
  ⠏ R016 (q+1 m1) / water / optimization    0:11:52
    all jobs ━━━━━━━━━━         0/2        0:12:47
```

To watch the calculation itself, tail ORCA's output in another terminal:

```bash
tail -f runs/<molecule>/water/<molecule>_water_opt.out
```

---

## 10. Troubleshooting

| symptom | cause and fix |
|---|---|
| `ORCA not found` | run `dft doctor`; install ORCA, put it on your `PATH` or set `$ORCA_PATH`, or set `execution.orca_path` |
| `topology discovery failed` | Open MPI cannot read the CPU topology; dftfig retries the job serially by itself. Set `execution.nprocs: 1` to skip the false start |
| `There are not enough slots available in the system` | `execution.nprocs` is above the **physical** core count. Open MPI allocates one slot per physical core, not per hardware thread, so a 6-core/12-thread laptop offers 6. The pipeline now clamps `nprocs` to the physical count when the config is loaded, and retries any job that still hits this at the real count |
| `Unknown identifier in PLOTS block` | an ORCA version older than 6; `%plots` syntax differs |
| job dies with `bad_alloc` / out of memory | lower `execution.maxcore_mb`; remember it is per process |
| `PyMOL executable not found` | `conda install -c conda-forge pymol-open-source`, or set `render.pymol_path` |
| MEP panel blank | `cubes.esp` was false, or `orca_vpot` failed — check for `esp.cube`; if you had set `cubes.esp_nprocs` above 1, set it back to 1 (see §8b) |
| LUMO positive, lobes floating off the molecule | a diffuse basis (`++`, `aug-`, `-TZVPD`) has no bound LUMO for a neutral molecule — see §4b. `property.basis_sp` now defaults to `def2-TZVP` for this reason |
| a failed job ends the whole run with `OSError: File name too long` | fixed in 1.4.0 — `parsers._text` treated a long output string as a filename. Update the package |
| Figure 13 panel missing, "atom ORDER differs" | the ORCA output's atom order does not match the RDKit molecule, so the Crippen logP values cannot be mapped; the message names the first differing position |
| optimization will not converge | raise `conformers.n_confs`; try `TightSCF`; inspect `*_trj.xyz` in Avogadro |
| small imaginary modes (−5 to −15 cm⁻¹) | usually grid noise in a methoxy rotation — try `grid: DefGrid3` and `opt_level: TightOpt` |
| genuine imaginary mode (< −50 cm⁻¹) | a real saddle point: distort along the mode and re-optimize |
| figure panels misaligned | `render.orientation: principal_plane` for planar systems; check each `render_*.png` individually |
| Fukui values all zero | the N±1 jobs did not finish — `python -m dftfig status` |

---

## 11. Citing the methods

Whatever you run, the caption needs: functional, both basis sets, the solvation model **named as
CPCM or SMD rather than IEF-PCM**, the population scheme, the isovalues for every surface, the MEP
colour range, and whether IE/EA are Koopmans or ΔSCF. `descriptors.csv` and the stamp on Figure 2
carry all of it.

---

## 12. Tests

```bash
pip install -e ".[dev]"
pytest -m "not orca"
```

The suite needs neither ORCA nor a network connection. It runs against genuine (trimmed)
ORCA 6.1 output checked into `tests/fixtures/`, and against `execution.dry_run: true`,
which generates every input file without launching anything.

| file | what it covers |
| --- | --- |
| `test_parsers.py` | every public function in `parsers.py`, against the real aniline outputs; also checks the trimmed fixtures still parse to the values the full files gave |
| `test_regression_long_output.py` | the `OSError: File name too long` crash of 1.4.0 — a failed job with a long output |
| `test_config.py` | `repair`, `validate` and `advisories` for the combinations that have actually cost a run |
| `test_discover.py` | ORCA verification, including the GNOME screen reader at `/usr/bin/orca` |
| `test_pka.py` | the pKa reference table, with literature values and sources |
| `test_integration_dryrun.py` | build + opt end to end: keyword line, `%pal`/`%maxcore`, solvation, per-compound charge |
| `test_figures.py` | Figures 12 and 14 render a real PNG from a checked-in results table |

Markers: `orca` needs an ORCA installation, `pymol` needs PyMOL. Both skip themselves
when the program is absent, and CI runs `pytest -m "not orca"`.
