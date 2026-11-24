# Using dftfig

Installing, running, choosing settings, and getting your molecules in.

<!-- Moved out of README.md to keep it short. Everything here was in the README up to v0.6.0. -->

---

## Installing

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
git clone https://github.com/MoazHosnyTaha/dftfig && cd dftfig
conda env create -f environment.yml
conda activate dftfig
pip install -e .
```

### With pip only

```bash
pip install dftfig
```

RDKit and everything else comes with it. PyMOL does not — the open-source build is not on
PyPI — so the 3D panels of `frontier_orbitals` and `lipophilicity` are skipped with a message rather than
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
dft verify        # re-check every file against runs/manifest.json
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


## One command

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
| 7. Volumetric data | the grids the frontier-orbitals figure is drawn from | 4 |
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
with a diffuse basis; the frontier-orbitals figure requested with cubes disabled.

---


## Stages

`dft run` executes all of them in order. You can also run one at a time — useful because the figures are quick to iterate on once the DFT is done:

```bash
python -m dftfig build       -c config.yaml   # SMILES -> starting geometries
python -m dftfig opt         -c config.yaml   # OPT (+FREQ), minimum check
python -m dftfig sp          -c config.yaml   # property single point, cubes
python -m dftfig fukui       -c config.yaml   # N-1 / N+1 single points
python -m dftfig cubes       -c config.yaml   # MEP via orca_vpot
python -m dftfig descriptors -c config.yaml   # descriptors.csv, fukui_*.csv, atom_descriptors.csv
python -m dftfig figures     -c config.yaml   # the four figures
python -m dftfig status      -c config.yaml   # what is done, what is missing
```

`execution.resume: true` (the default) means a stage re-run skips anything that already
terminated normally, so you can restyle a figure without repeating hours of DFT.

**Before a long campaign**, set `execution.dry_run: true` and run `all`. Every input file is
written and nothing is launched; read the `.inp` files, confirm the keyword line, then turn it off.

---


## Using your own molecules

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


## Wall-clock expectations

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


## Interrupting a run

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


## What comes out

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
  descriptors.csv                   one row per molecule × medium: electronic + steric
  atom_descriptors.csv              one row per atom: Fukui + steric descriptors
  manifest_{opt,sp,fukui}.json      per-job status, wall time, energies
  manifest.json                     provenance: versions, settings hash, SHA-256 per file
  figures/
    frontier_orbitals.{png,pdf}  speciation.*  lipophilicity.*  mpo_metabolism.*
    properties_resolved.csv         what the figures actually used
```

Every `.cube` opens directly in **Avogadro** or **GaussView** if you want to look around
interactively, and every `.pml` can be opened in PyMOL (`pymol render_mep_CCB1.pml`) to adjust a
camera angle by hand and re-save.

---

---

[← back to the README](../README.md)
