# dftfig — complete parameter reference
Generated from the registry in `dftfig/params.py` by `tools/gen_param_docs.py`. `dft params` shows exactly this information, one parameter at a time, in colour.
`paper` = reproduces the chlocarbazomycin manuscript · `recommended` = what I would run for new work · `fast` = plumbing test

---

## Molecules  (`molecules`)
*What to calculate. One SDF per compound is the normal route; SMILES also works.*

### `molecules.structure_source`
**Where structures come from** · type `choice` · default `sdf`

The input format for your compounds.

**How to choose.** 'sdf' is the normal route: drop one .sdf per compound into a folder and the pipeline reads them all, taking each compound's charge and 3D geometry from the file. 'smiles' is convenient for a handful of compounds typed by hand. 'xyz' is for geometries you already trust and do not want a conformer search to disturb - it still needs a SMILES per molecule for bond orders.

| value | name | what it means | preset |
|---|---|---|---|
| `sdf` | SDF files (one per compound) | Reads every .sdf in molecules.sdf_dir. Charge is taken from the formal charges in the file, and 3D coordinates are used directly when present. | `recommended` `fast` |
| `smiles` | SMILES typed in the config | Name-to-SMILES entries under molecules.definitions. A conformer is generated for each. | `paper` |
| `xyz` | Existing .xyz geometries | Reads <name>.xyz from molecules.xyz_dir and skips the conformer search. Still requires molecules.definitions for the chemistry perception. |  |

### `molecules.sdf_dir`
**Folder of SDF files** · type `path` · default `molecules`

Directory searched for compound files when structure_source = 'sdf'. Relative paths are resolved from the folder you run `dft` in.

**How to choose.** One .sdf per compound is cleanest - the file name becomes the molecule name, which becomes the output directory and the figure panel label. A multi-record SDF also works: each record becomes its own molecule. Export from Avogadro, ChemDraw, MOE or RDKit; keep 3D coordinates if you have them and the pipeline will start the optimization from them.

### `molecules.sdf_pattern`
**File pattern** · type `str` · default `*.sdf`

Glob matched inside sdf_dir.

**How to choose.** Widen to '*.sd*' to pick up .sd and .sdf together, or narrow to a prefix like 'series1_*.sdf' to run one subset without moving files around.

### `molecules.name_from`
**Molecule naming** · type `choice` · default `filename`

Where each molecule's name comes from.

**How to choose.** The name becomes a directory name, a filename prefix and a figure panel label, so short and filesystem-safe is best. 'filename' is predictable; 'sdf_title' uses the title line inside the file, which is handy when the files came from a database export with meaningless names. Characters that are awkward in paths are replaced with underscores either way.

| value | name | what it means | preset |
|---|---|---|---|
| `filename` | The .sdf file name | compound_17.sdf -> compound_17. | `recommended` `fast` |
| `sdf_title` | The title line inside the file | The first line of the SDF record, when it is not empty. |  |

### `molecules.only`
**Run only these compounds** · type `str` · default ``

Comma-separated compound names. Empty means every structure found.

**How to choose.** Names are the SDF file names without the extension - the same names that appear in runs/molecules.csv and on the figure panels. Use this to re-run one compound that failed without touching the rest, or to produce a figure from a subset of a large batch. `dft run` sets it for you when you deselect compounds in the list it shows.

### `molecules.limit`
**How many compounds to run** — compounds · type `int` · default `0` · range `0–`

Run only the first N compounds found. 0 means all of them.

**How to choose.** Set this before a long campaign so you know exactly what is queued: with 40 SDF files in a folder, limit = 3 runs three of them and tells you which. Use it to time one compound before committing to the batch - `dft validate` reports the resulting job count. Files are taken in sorted order, so the selection is reproducible.

### `molecules.charge_mode`
**Charge handling** · type `choice` · default `auto`

How each compound's net charge is decided.

**How to choose.** 'auto' reads the formal charges in the structure, per compound. This is what you want whenever a set is not uniformly neutral - a quaternary ammonium in the same batch as a neutral nucleoside has a different electron count, and running it at charge 0 computes a different species entirely (wrong electron count, wrong SCF, meaningless orbitals). 'fixed' forces one charge on every compound and is only correct for a uniform set. The detected charge of every molecule is printed before any job starts and written to molecules.csv, so an unexpected value is visible immediately.

| value | name | what it means | preset |
|---|---|---|---|
| `auto` | Detect per compound | Sum of formal charges in the SMILES or SDF. Different compounds in one run may carry different charges. | `recommended` `fast` |
| `fixed` | Force molecules.charge on all of them | Every compound runs at the same charge. Correct only for a uniform set. | `paper` |

### `molecules.charge`
**Fixed charge** · type `int` · default `0`

The charge applied to every compound when charge_mode = 'fixed'. Ignored in 'auto' mode.

**How to choose.** Use this to force a charge state the structure does not encode - for example to study the anion of a neutral acid without editing the input file. Remember it then applies to the whole set.

### `molecules.multiplicity_mode`
**Spin multiplicity handling** · type `choice` · default `auto`

How each compound's spin multiplicity is decided.

**How to choose.** 'auto' picks the lowest spin state consistent with the electron count: a singlet for an even number of electrons, a doublet for an odd number. That is right for closed-shell organics and for radical ions. Use 'fixed' only when you specifically want a higher spin state, such as a triplet; the pipeline rejects a multiplicity that is impossible for the electron count rather than letting ORCA fail halfway through.

| value | name | what it means | preset |
|---|---|---|---|
| `auto` | Lowest spin state for the electron count | Singlet if even, doublet if odd. | `recommended` `fast` |
| `fixed` | Force molecules.multiplicity | For triplets and other deliberately high-spin states. | `paper` |

### `molecules.multiplicity`
**Fixed multiplicity** · type `int` · default `1` · range `1–`

2S+1 applied to every compound when multiplicity_mode = 'fixed'. 1 = closed-shell singlet, 2 = doublet, 3 = triplet.

**How to choose.** Anything other than 1 switches ORCA to unrestricted DFT (UKS), which doubles the orbital and population output and makes the frontier-orbital panels ambiguous (alpha and beta sets differ).

### `molecules.definitions`
**Molecule definitions (SMILES)** · type `kvlist` · default `{}`

Name-to-SMILES map, used when structure_source is 'smiles' or 'xyz'. Ignored for 'sdf'.

**How to choose.** Empty by default - this pipeline has no built-in compounds. Add your own as NAME: SMILES, or switch to SDF input. Panel order in every figure follows the order given here. Keep names short: they become directory names and figure labels. For the chlocarbazomycin reference set from the manuscript, use the `paper` preset or run `dft init --demo`.

### `molecules.use_sdf_geometry`
**Use 3D coordinates from the file** · type `bool` · default `True`

When an SDF record already contains a 3D conformer, start the DFT optimization from it instead of generating a new one.

**How to choose.** Keep this on when your SDFs came from a docking pose, a crystal structure or a careful manual build - those geometries encode information a conformer generator cannot recover. Turn it off to force a fresh conformer search for every compound, which is the right choice if the SDFs were produced from 2D depictions or you do not trust their provenance. Records with no 3D coordinates always get a generated conformer regardless.

### `molecules.xyz_dir`  *(advanced)*
**Directory of .xyz files** · type `path` · default `input_xyz`

Only read when structure_source = 'xyz'. Expects one <molecule-name>.xyz per entry in molecules.definitions.

**How to choose.** Plain XYZ format (atom count, comment line, then element plus Cartesian coordinates in Angstrom). Avogadro and GaussView both export this.

## Conformer search  (`conformers`)
*Finding a sensible starting structure before DFT.*

### `conformers.n_confs`
**Conformers to generate** — structures · type `int` · default `50` · range `1–`

How many trial 3D structures RDKit embeds before force-field ranking. The lowest-energy survivor is what DFT then optimizes.

**How to choose.** For rigid fused-ring systems like carbazoles only the methoxy rotation and the hydroxyl orientation matter, so 20-50 is plenty. Push to 300+ for molecules with more than about 5 rotatable bonds, otherwise you risk starting DFT from a geometry that is not the global minimum - a far more expensive mistake than the conformer search itself.

**Cost.** Seconds. Never the bottleneck.

### `conformers.forcefield`
**Force field for pre-optimization** · type `choice` · default `MMFF94s`

The classical force field used to relax and rank the embedded conformers.

**How to choose.** MMFF94s is the right default for drug-like organics. Switch to UFF only if MMFF refuses to parameterize your molecule (unusual elements, exotic valences) - UFF geometries are noticeably worse but DFT will fix them.

| value | name | what it means | preset |
|---|---|---|---|
| `MMFF94s` | MMFF94s | Static variant of MMFF94 with planarized sp2 nitrogens. Best general choice for aromatic heterocycles - it will keep the carbazole NH planar instead of pyramidalizing it. | `paper` `recommended` `fast` |
| `MMFF94` | MMFF94 | Original Merck force field. Slightly different nitrogen out-of-plane treatment. |  |
| `UFF` | UFF | Universal Force Field. Covers the whole periodic table but with poor accuracy for organics. Fallback only. |  |

### `conformers.prune_rms`
**Duplicate pruning RMSD** — A · type `float` · default `0.5` · range `0.0–`

Conformers closer than this heavy-atom RMSD are treated as the same structure and collapsed to one.

**How to choose.** 0.5 A removes trivial duplicates without discarding genuinely different rotamers. Lower it toward 0.1 if you want to inspect fine differences; raise it to 1.0 to keep only clearly distinct conformers.

### `conformers.random_seed`
**Random seed** · type `int` · default `42`

Seed for the ETKDG embedding. Fixing it makes the whole pipeline bit-reproducible.

**How to choose.** Keep it fixed for anything you intend to publish, so a reviewer rerunning the pipeline gets your exact geometries. Change it only to test that your result is not an artefact of one particular embedding.

### `conformers.ff_max_iters`  *(advanced)*
**Force-field iterations** · type `int` · default `2000`

Maximum optimization steps per conformer at the force-field level.

**How to choose.** 2000 converges essentially any drug-like molecule. Raise only if the log reports unconverged conformers.

## Level of theory (geometry optimization)  (`theory`)
*The functional, basis set and convergence settings for the OPT+FREQ step.*

### `theory.functional`
**Exchange-correlation functional** · type `choice` · default `B3LYP`

The DFT approximation. Controls accuracy of geometries, orbital energies and every descriptor derived from them.

**How to choose.** B3LYP reproduces the manuscript. Be aware that bare B3LYP has NO London dispersion, which matters because the paper's structural argument is about pi-stacking - see theory.dispersion below. If you are starting fresh rather than reproducing, wB97X-D4 or r2SCAN-3c are better value for the same money.

**Cost.** Hybrids (B3LYP, PBE0, M06-2X) need exact exchange and are ~3-5x a pure GGA. Range-separated hybrids (wB97X) add another ~30%.

| value | name | what it means | preset |
|---|---|---|---|
| `B3LYP` | B3LYP | The functional used in the manuscript. Reliable ground-state geometries for organics; underestimates barriers; no dispersion of its own. | `paper` |
| `B3LYP D3BJ` | B3LYP-D3BJ | B3LYP plus Grimme's D3 dispersion with Becke-Johnson damping. Essentially free, and strictly better than bare B3LYP for anything involving stacking or intramolecular contacts. |  |
| `PBE0 D3BJ` | PBE0-D3BJ | Parameter-free hybrid, 25% exact exchange. Often better orbital energies and gaps than B3LYP. |  |
| `wB97X-D4` | wB97X-D4 | Range-separated hybrid with D4 dispersion. The most reliable general-purpose choice today for geometries, gaps and non-covalent interaction energies. | `recommended` |
| `M062X` | M06-2X | Heavily parameterized meta-hybrid, good for main-group thermochemistry and non-covalent complexes. Needs a dense integration grid (DefGrid3). |  |
| `r2SCAN-3c` | r2SCAN-3c | Composite method: meta-GGA with its own fitted triple-zeta basis and dispersion baked in. Very fast and surprisingly accurate - ignores the basis_opt setting. | `fast` |

### `theory.basis_opt`
**Basis set (optimization)** · type `choice` · default `6-31+G(d,p)`

The atomic-orbital basis used for the geometry optimization and frequencies.

**How to choose.** 6-31+G(d,p) is the manuscript's choice: split-valence double-zeta, polarization on all atoms, diffuse functions on heavy atoms. The diffuse set is what makes it slow, and it is what you need for anions and for LUMO energies to mean anything. def2-SVP is faster and usually gives indistinguishable geometries.

**Cost.** Scales steeply: def2-SVP ~ 1x, 6-31+G(d,p) ~ 2-3x, def2-TZVP ~ 6-10x for the same molecule.

| value | name | what it means | preset |
|---|---|---|---|
| `6-31+G(d,p)` | 6-31+G(d,p) | The manuscript's basis. Diffuse functions on heavy atoms (+), polarization on heavy atoms and hydrogens (d,p). | `paper` |
| `6-31G(d)` | 6-31G(d) | No diffuse functions, no hydrogen polarization. Noticeably cheaper; fine for geometries, poor for LUMO energies and anions. |  |
| `def2-SVP` | def2-SVP | Karlsruhe double-zeta. Modern, well-balanced, pairs with the def2/J auxiliary basis for RIJCOSX. Good geometry-only choice. | `fast` |
| `def2-TZVP` | def2-TZVP | Triple-zeta. Use when orbital energies or energetics - not just geometry - are the deliverable. | `recommended` |
| `ma-def2-SVP` | ma-def2-SVP | Minimally augmented def2-SVP: diffuse functions at a fraction of the cost of a full aug- set. A good substitute for 6-31+G(d,p). |  |

### `theory.dispersion`
**Dispersion correction** · type `choice` · default `none`

An additive London-dispersion term. Density functionals without it systematically miss attraction between non-bonded fragments.

**How to choose.** Set to 'none' to reproduce the manuscript exactly. Set to D3BJ or D4 for any new work: it costs essentially nothing and the paper's own pi-stacking narrative needs it. Leave at 'none' if the functional already includes dispersion by construction (wB97X-D4, r2SCAN-3c) - adding it twice is a real error.

**Cost.** Negligible - a pairwise atom-atom sum, microseconds.

| value | name | what it means | preset |
|---|---|---|---|
| `none` | None | No dispersion. What the manuscript did. Acceptable for covalent geometries of a single rigid molecule; not acceptable for stacking or conformer energetics. | `paper` |
| `D3BJ` | D3(BJ) | Grimme's D3 with Becke-Johnson damping. The standard pairing for B3LYP and PBE0. | `recommended` `fast` |
| `D4` | D4 | Newer, charge-dependent dispersion. Slightly better, slightly slower, needs parameters for your functional. |  |

### `theory.ri_approx`
**Integral approximation** · type `choice` · default `RIJCOSX`

How the Coulomb (J) and exchange (K) integrals are evaluated. RI = resolution of the identity, an auxiliary-basis expansion; COSX = chain-of-spheres numerical exchange.

**How to choose.** RIJCOSX makes hybrid DFT 5-20x faster with errors of ~0.1 kcal/mol on relative energies and under 0.001 eV on orbital energies - invisible at the precision anyone reports. Set to 'none' only if you want integrals handled exactly as Gaussian 09 would have, for a strict numerical comparison with the manuscript.

**Cost.** RIJCOSX typically cuts hybrid-DFT SCF time by 5-20x. 'none' is the slow, exact path.

| value | name | what it means | preset |
|---|---|---|---|
| `RIJCOSX` | RIJCOSX | RI for Coulomb, chain-of-spheres for exchange. The standard production setting for hybrid DFT in ORCA. Requires an auxiliary basis (see theory.aux_basis). | `recommended` `fast` |
| `none` | None (exact, conventional) | Four-index integrals computed exactly, as Gaussian 09 does by default. Slow, but removes any approximation when comparing numbers against the published values. | `paper` |
| `RIJK` | RIJK | RI for both J and K. More accurate than COSX for exchange but needs a /JK auxiliary basis and more memory. |  |

### `theory.aux_basis`  *(advanced)*
**Auxiliary basis** · type `str` · default `def2/J`

The fitting basis used by RIJCOSX/RIJK. Ignored when ri_approx = none.

**How to choose.** def2/J works with every Karlsruhe and Pople basis for main-group elements and is the right answer in almost all cases. Use def2/JK if you selected RIJK.

### `theory.scf_conv`
**SCF convergence** · type `choice` · default `TightSCF`

How tightly the self-consistent field is converged at each geometry step.

**How to choose.** TightSCF is the sane default and is required for reliable frequencies. NormalSCF can leave enough noise in the gradient to stall a geometry optimization. VeryTightSCF only for numerical frequencies or very flat potential surfaces.

**Cost.** Each tier is roughly 10-25% more SCF iterations than the one below.

| value | name | what it means | preset |
|---|---|---|---|
| `TightSCF` | TightSCF | Energy change < 1e-8 Eh. Required for clean analytic frequencies. | `paper` `recommended` |
| `NormalSCF` | NormalSCF | Energy change < 1e-6 Eh. Fine for a quick geometry. | `fast` |
| `VeryTightSCF` | VeryTightSCF | Energy change < 1e-9 Eh. For numerical gradients/Hessians and flat surfaces. |  |

### `theory.grid`
**DFT integration grid** · type `choice` · default `DefGrid2`

The numerical grid on which the exchange-correlation functional is integrated.

**How to choose.** DefGrid2 is ORCA 6's default and is adequate for B3LYP. Move to DefGrid3 for Minnesota functionals (M06-2X and relatives), which are notoriously grid-sensitive, or if frequencies come out with small spurious imaginary modes.

**Cost.** DefGrid3 is roughly 2-3x DefGrid2 in the XC build step.

| value | name | what it means | preset |
|---|---|---|---|
| `DefGrid2` | DefGrid2 | ORCA 6 default. Fine for B3LYP, PBE0, wB97X. | `paper` `recommended` `fast` |
| `DefGrid1` | DefGrid1 | Coarse. Pre-screening only - do not report frequencies from it. |  |
| `DefGrid3` | DefGrid3 | Dense. Required for M06-2X; cures spurious tiny imaginary modes. |  |

### `theory.opt_level`
**Geometry convergence** · type `choice` · default `Opt`

The thresholds on energy change, gradient and displacement that end the optimization.

**How to choose.** 'Opt' matches Gaussian's default and is right for almost everything. TightOpt is worth it when you are about to compute frequencies on a very flat system and want to avoid residual imaginary modes.

| value | name | what it means | preset |
|---|---|---|---|
| `Opt` | Opt (normal) | Gradient < 3e-4 Eh/bohr. Standard. | `paper` `recommended` |
| `LooseOpt` | LooseOpt | Quick and dirty; do not follow with a frequency job. | `fast` |
| `TightOpt` | TightOpt | Gradient < 1e-4 Eh/bohr. Flat surfaces and pre-frequency polish. |  |

### `theory.run_freq`
**Run frequency calculation** · type `bool` · default `True`

Computes the Hessian after optimization. This is the ONLY way to prove your structure is a genuine minimum rather than a saddle point.

**How to choose.** Keep this on. The pipeline parses the frequencies and flags any structure with imaginary modes instead of silently carrying it forward - which is exactly the check the manuscript says it performed but never reports. Turn off only for a plumbing test.

**Cost.** Analytic Hessian for a hybrid functional roughly doubles to triples the total job time.

### `theory.imaginary_cutoff`  *(advanced)*
**Imaginary-mode tolerance** — cm^-1 · type `float` · default `-10.0`

Modes more negative than this are reported as genuine imaginary frequencies; small negative values above it are treated as grid noise.

**How to choose.** -10 cm^-1 is the usual convention. A mode at -5 cm^-1 is almost always numerical noise in a floppy methoxy rotation; a mode at -150 cm^-1 means you are on a saddle point and must re-optimize along that mode.

### `theory.saddle_max_attempts`
**Automatic restarts off a saddle point** — restarts · type `int` · default `2` · range `0–5`

How many times to re-optimize a structure that came back with imaginary frequencies.

**How to choose.** A frequency job reporting an imaginary mode has found a saddle point, not a minimum - the gap, thermochemistry, Fukui indices and MEP all then describe a structure that does not exist as a stable species. With this above 0 the pipeline fixes it instead of just reporting it: it reads the imaginary normal mode, pushes the geometry along it and optimizes again, trying the opposite direction if the first falls back into the same saddle. 2 clears most cases - an amino group that needs to pyramidalise, a methyl stuck eclipsed, a ring at its flat inversion point. Set 0 to get the old behaviour of a loud warning and nothing else. Each restart is a full OPT+FREQ, so budget for it.

**Cost.** One extra OPT+FREQ per restart, for affected molecules only.

### `theory.saddle_step_angstrom`  *(advanced)*
**Restart displacement** — A · type `float` · default `0.3` · range `0.05–1.0`

How far the atom that moves most is pushed along the imaginary mode before re-optimizing.

**How to choose.** 0.3 A is enough to leave the saddle without landing somewhere unrelated. Too small and the optimizer walks straight back; too large and you can cross into a different conformer, which is a different molecule for every purpose except the formula. Raise it only if a restart keeps returning the same imaginary mode.

### `theory.saddle_small_mode_cm`  *(advanced)*
**Small-mode threshold** — cm^-1 · type `float` · default `50.0` · range `0.0–300.0`

Below this magnitude an imaginary mode is treated as a numerical artifact first and a real saddle second.

**How to choose.** An imaginary frequency of a few tens of wavenumbers is usually a methyl or hydroxyl rotation sitting in integration-grid noise, and the first thing to try is a finer grid with a tighter optimization rather than a displacement along a nearly flat coordinate. Above the threshold the saddle is taken as real and the geometry is displaced straight away.

### `theory.temperature`  *(advanced)*
**Thermochemistry temperature** — K · type `float` · default `298.15`

Temperature for the ideal-gas/rigid-rotor/harmonic-oscillator free-energy corrections.

**How to choose.** 298.15 K for standard-state reporting. Use 310.15 K if you want the thermochemistry to match the 310 K the manuscript used for its MD simulations.

## Implicit solvation  (`solvation`)
*The continuum-solvent treatment, applied to every stage.*

### `solvation.model`
**Solvation model** · type `choice` · default `CPCM`

The implicit-solvent model. The solute sits in a cavity inside a polarizable continuum with the dielectric constant of the chosen solvent.

**How to choose.** IMPORTANT: ORCA does not implement Gaussian's IEF-PCM. CPCM (conductor-like PCM) is the closest equivalent and the two agree to well under 0.1 kcal/mol in high-dielectric solvents such as water (eps 78) and DMSO (eps 47). They diverge more in low-dielectric media - chloroform (eps 4.7) is where you would see the difference. SMD adds a non-electrostatic cavitation/dispersion/repulsion term and is the better model for solvation free energies, though it changes orbital energies slightly relative to the published CPCM/IEF-PCM numbers.

| value | name | what it means | preset |
|---|---|---|---|
| `CPCM` | CPCM | Conductor-like PCM. ORCA's closest analogue to the manuscript's IEF-PCM; electrostatics only. | `paper` `fast` |
| `SMD` | SMD | Truhlar's solvation model: CPCM electrostatics plus an empirical non-electrostatic term. Best choice for solvation free energies and for logP-type reasoning. | `recommended` |
| `none` | Gas phase only | No continuum. Use when you want bare vertical ionization energies and electron affinities, which is the only context where Koopmans-type descriptors are physically interpretable. |  |

### `solvation.solvents`
**Solvent list** · type `multichoice` · default `['water']`

One full set of calculations is run per entry. 'gas' adds a vacuum reference.

**How to choose.** The manuscript reports vacuum plus water, methanol, chloroform and DMSO. Each solvent multiplies the total cost by one, so start with water alone to validate the workflow, then add the rest. Figures 2, 12, 13 and 14 are all built from the solvent named in figures.reference_solvent.

**Cost.** Linear: 4 solvents = 4x the DFT wall time.

| value | name | what it means | preset |
|---|---|---|---|
| `water` | Water (eps 78.4) | The manuscript's reference medium for Figure 2. | `paper` `recommended` `fast` |
| `methanol` | Methanol (eps 32.6) | Polar protic. | `paper` |
| `chloroform` | Chloroform (eps 4.7) | Low dielectric; where CPCM and IEF-PCM differ most. | `paper` |
| `dmso` | DMSO (eps 46.7) | Polar aprotic. | `paper` |
| `gas` | Gas phase | No continuum. The only rigorous reference for ionization energies. | `paper` |
| `acetonitrile` | Acetonitrile (eps 36.6) | Polar aprotic alternative. |  |
| `toluene` | Toluene (eps 2.4) | Apolar aromatic. |  |
| `thf` | THF (eps 7.3) | Ethereal, low dielectric. |  |

## Property single point  (`property`)
*A second, larger-basis calculation on the optimized geometry - the source of all descriptors.*

### `property.enabled`
**Run property single point** · type `bool` · default `True`

Performs a single-point calculation at a larger basis on the optimized geometry, and takes the orbital energies, populations and cubes from it.

**How to choose.** Keep on. This is the step the manuscript implies but never describes: its Methods give 6-31+G(d,p) while the MEP section quotes 6-311++G(d,p), with no single point in between. Making it explicit is what lets you label every figure with the level of theory that actually produced it.

**Cost.** One SCF, no gradient - typically 5-15% of the optimization cost.

### `property.basis_sp`
**Basis set (property step)** · type `choice` · default `def2-TZVP`

Basis for the single point that generates orbital energies, populations, cubes and MEP.

**How to choose.** This choice is a trade-off between two things you want from the same calculation, and the right answer depends on which one the figure is for. Fully diffuse sets (the '++' in 6-311++G(d,p), the 'aug-' in aug-cc-pVTZ) describe density far from the nuclei, which is where the MEP tail and any genuine anion live. But a neutral organic molecule has a NEGATIVE electron affinity: its LUMO is not a bound state at all, and a diffuse set describes that unbound state as a discretised piece of the continuum. The plotted 'LUMO' then comes out as a vague cloud sitting off the molecule instead of the pi* orbital, and the HOMO-LUMO gap it reports is an artifact of the basis rather than a property of the compound - benzene, phenol, aniline and toluene all gave positive LUMO energies and 9-11 eV gaps at 6-311++G(d,p), against bound pi* LUMOs and 5.5-6.7 eV at def2-TZVP. So: def2-TZVP for an orbital figure, ma-def2-TZVP if you also want a well-converged MEP, and a fully diffuse set only when you are deliberately studying an anion.

**Cost.** Triple-zeta single point on an optimized geometry: minutes for a drug-sized molecule. The diffuse sets are roughly 1.5-2x the cost of def2-TZVP and make the SCF harder to converge.

| value | name | what it means | preset |
|---|---|---|---|
| `def2-TZVP` | def2-TZVP | Karlsruhe triple-zeta, no diffuse functions. Bound pi* LUMO with the right nodal pattern, stable orbital energies, well-behaved populations. The right default for an orbital or gap figure. | `recommended` |
| `ma-def2-TZVP` | ma-def2-TZVP | Minimally augmented: one diffuse shell per atom. Recovers the MEP tail and anion energetics without the continuum-LUMO problem. The compromise choice when one single point has to serve both the MEP and the orbitals. |  |
| `def2-TZVPD` | def2-TZVPD | Fully diffuse Karlsruhe set. For anions and electron affinities, where the extra diffuse functions describe a state that really is bound. |  |
| `6-311++G(d,p)` | 6-311++G(d,p) | Pople triple-zeta, diffuse on heavy atoms and hydrogens. Common in the literature. Expect an unbound LUMO for neutral molecules. |  |
| `def2-SVP` | def2-SVP | Double-zeta. Too small for published orbital energies; useful for checking that the pipeline runs. | `fast` |

### `property.population`
**Population analysis** · type `choice` · default `Hirshfeld`

The scheme that partitions the total electron density into per-atom charges. These charges are what the condensed Fukui functions are built from.

**How to choose.** Hirshfeld is the manuscript's choice and the right one: it is basis-set stable, always yields positive atomic populations, and is the standard for condensed Fukui indices. Note that Gaussian 09 - the code the paper cites - has no native Hirshfeld analysis, so those numbers must have come from an external tool that is never named. ORCA computes it directly.

| value | name | what it means | preset |
|---|---|---|---|
| `Hirshfeld` | Hirshfeld | Stockholder partitioning against promolecular atomic densities. Basis-set robust; the standard for condensed Fukui functions. | `paper` `recommended` |
| `Loewdin` | Loewdin | Symmetrically orthogonalized Mulliken. More stable than Mulliken, still basis-dependent. |  |
| `Mulliken` | Mulliken | The classic. Badly behaved with diffuse functions - can give negative populations. Avoid with any '+' or 'aug-' basis. | `fast` |
| `Mayer` | Mayer | Gives bond orders and valences in addition to charges. |  |

### `property.koopmans_mode`
**Descriptor derivation** · type `choice` · default `both`

How ionization energy and electron affinity are obtained for the global reactivity descriptors.

**How to choose.** 'koopmans' takes IE = -e(HOMO) and EA = -e(LUMO), which is what the manuscript did. Be aware this is a Hartree-Fock theorem: in Kohn-Sham DFT the HOMO equals the exact IE only for the exact functional, and B3LYP typically underestimates IE by 1-3 eV (the paper's 5.5 eV for a carbazole should be compared with ~7.6 eV experimental). 'delta_scf' instead computes real energy differences from the N+1 and N-1 calculations - physically meaningful, and free if you are already running Fukui. 'both' reports the two side by side, which is the honest option.

**Cost.** delta_scf needs two extra single points per molecule - already paid for if fukui.enabled.

| value | name | what it means | preset |
|---|---|---|---|
| `both` | Both (recommended) | Report Koopmans and delta-SCF values together so the approximation is visible. | `recommended` |
| `koopmans` | Koopmans only | IE = -e(HOMO), EA = -e(LUMO). Reproduces the manuscript's Table S1 approach. | `paper` `fast` |
| `delta_scf` | Delta-SCF only | IE = E(cation) - E(neutral), EA = E(neutral) - E(anion), all at the neutral geometry. Vertical, and physically defensible. |  |

## Fukui functions  (`fukui`)
*Local reactivity: which atoms are electrophilic or nucleophilic.*

### `fukui.enabled`
**Compute condensed Fukui functions** · type `bool` · default `True`

Runs the N+1 and N-1 single points needed for f+, f-, f0, the dual descriptor and the multiphilic descriptor.

**How to choose.** Required for the Fukui panels and the supplementary tables. Also supplies the delta-SCF ionization energies and electron affinities at no extra cost.

**Cost.** Two additional single points per molecule per solvent - roughly 2x the property step.

### `fukui.geometry`
**Geometry for the ion calculations** · type `choice` · default `frozen`

Whether the cation and anion are computed at the neutral geometry or re-optimized.

**How to choose.** 'frozen' gives VERTICAL quantities, which is what Fukui functions are defined from and what you want for reactivity indices. 'relaxed' gives adiabatic values, which are the right thing for redox potentials but NOT for Fukui analysis. Keep frozen unless you specifically want adiabatic ionization energies.

| value | name | what it means | preset |
|---|---|---|---|
| `frozen` | Frozen (vertical) | Single points at the optimized neutral geometry. Correct for Fukui functions and vertical IE/EA. | `paper` `recommended` `fast` |
| `relaxed` | Relaxed (adiabatic) | Re-optimize each ion. Gives adiabatic IE/EA; much more expensive and not what Fukui indices are defined from. |  |

### `fukui.smear_scheme`  *(advanced)*
**Degenerate-orbital handling** · type `choice` · default `none`

Whether to average Fukui indices over near-degenerate frontier orbitals.

**How to choose.** Leave at 'none' for closed-shell organics. Symmetric molecules with degenerate frontier orbitals can give misleading single-orbital pictures; in that case the finite-difference route used here is already the robust answer.

| value | name | what it means | preset |
|---|---|---|---|
| `none` | None | Straight finite differences of the atomic populations. | `paper` `recommended` |

## Volumetric data (cube files)  (`cubes`)
*The 3D grids behind the orbital and MEP panels of Figure 2.*

### `cubes.enabled`
**Generate cube files** · type `bool` · default `True`

Writes Gaussian-format .cube files for the frontier orbitals, the electron density and the electrostatic potential.

**How to choose.** Needed for Figure 2 panels B and C and for Figure 13. The cubes are also directly loadable in Avogadro and GaussView if you want to look around interactively.

**Cost.** Cube generation is cheap; the MEP evaluation is the slow part (see cubes.esp).

### `cubes.grid_points`
**Grid points per axis** — points · type `int` · default `80` · range `20–200`

Resolution of the volumetric grid. The cube contains this value cubed points in total.

**How to choose.** 80 gives smooth publication-quality isosurfaces for a 30-atom molecule. 40 is visibly faceted and should only be used for testing. Going past 120 produces files of hundreds of MB for no visible gain. Remember the cost is cubic: 80 -> 512,000 points, 120 -> 1,728,000.

**Cost.** Cubic in this number. Dominates MEP evaluation time and disk usage.

### `cubes.orbitals`
**Which orbitals to plot** · type `choice` · default `homo_lumo`

Selects the molecular orbitals written as cubes.

**How to choose.** homo_lumo covers Figure 2B. Extend to homo_lumo_pm1 when a frontier orbital is nearly degenerate with the one below/above it and a single-orbital picture would mislead - check the orbital energies in descriptors.csv first.

| value | name | what it means | preset |
|---|---|---|---|
| `homo_lumo` | HOMO and LUMO | The two frontier orbitals. What Figure 2B shows. | `paper` `recommended` `fast` |
| `homo_lumo_pm1` | HOMO-1 ... LUMO+1 | Four orbitals - use when frontier orbitals are near-degenerate. |  |
| `homo_lumo_pm2` | HOMO-2 ... LUMO+2 | Six orbitals. Rarely needed for a figure. |  |

### `cubes.esp`
**Compute electrostatic potential** · type `bool` · default `True`

Evaluates V(r) at every grid point with orca_vpot and writes a matching ESP cube, which is then painted onto the density isosurface to make the MEP map.

**How to choose.** Required for Figure 2C. Note this is a genuine quantum-mechanical MEP from the SCF density, not a point-charge approximation.

**Cost.** Scales as (basis functions)^2 x (grid points). Roughly 5-10 min for a 30-atom molecule on an 80^3 grid at triple-zeta; use cubes.esp_nprocs to parallelize.

### `cubes.esp_nprocs`  *(advanced)*
**Parallel processes for MEP** — processes · type `int` · default `1` · range `1–64`

Number of MPI processes for orca_vpot. 1 runs the serial binary, which is the only form that works reliably after a completed job.

**How to choose.** Leave this at 1. The parallel build of orca_vpot re-reads the SHARK basis scratch file <job>.SHARKINP.tmp, and ORCA deletes that file whenever a calculation terminates normally - so the parallel path fails for exactly the jobs you want a MEP from, with a 'Failed to read NAtoms, NBasis, NShells' error. The serial binary rebuilds everything from the .gbw and is fast enough: 512,000 grid points on a 12-atom molecule takes about 37 s. If you do raise it, a failure falls back to serial automatically.

### `cubes.keep_cubes`  *(advanced)*
**Keep cube files after rendering** · type `bool` · default `True`

Whether to retain the .cube files once the figures have been produced.

**How to choose.** Keep them if you want to re-render with different isovalues or inspect the orbitals in Avogadro later - regenerating means rerunning the SCF. Set to false only if disk space is tight: an 80^3 cube is about 14 MB and you get four or more per molecule per solvent.

## 3D rendering  (`render`)
*How the orbital, MEP and lipophilicity surfaces are drawn.*

### `render.engine`
**Rendering engine** · type `choice` · default `pymol`

The program that turns cube files into images.

**How to choose.** PyMOL is scripted headlessly, needs no clicking, and produces consistent camera angles across panels - which is what makes a multi-molecule figure readable. VMD is an alternative if you prefer its isosurface shading. 'none' skips rendering and leaves you the cubes to open in Avogadro or GaussView yourself.

| value | name | what it means | preset |
|---|---|---|---|
| `pymol` | PyMOL | Scripted via generated .pml files and ray-traced offscreen. Best control over colour ramps and camera. | `paper` `recommended` `fast` |
| `vmd` | VMD | Scripted via generated .tcl files with Tachyon rendering. Nicer default isosurface shading; fiddlier colour ramps. |  |
| `none` | None - cubes only | Stop after cube generation. Open the .cube files in Avogadro or GaussView and render by hand. |  |

### `render.width`
**Image width** — px · type `int` · default `1600` · range `200–`

Pixel width of each rendered 3D panel before assembly into the figure.

**How to choose.** 1600 px per panel gives a 4-panel row at roughly 300 dpi across a journal page. Raise to 2400 for a full-page figure; ray-trace time scales with the pixel count.

### `render.height`
**Image height** — px · type `int` · default `1200` · range `200–`

Pixel height of each rendered 3D panel.

**How to choose.** Keep the aspect ratio close to the shape of your molecules. Carbazoles are wide and flat, so 4:3 or 16:9 wastes less space than a square.

### `render.ray_trace`
**Ray-trace images** · type `bool` · default `True`

Use PyMOL's ray tracer rather than the OpenGL preview for the final image.

**How to choose.** On for anything going into a manuscript - it gives proper shadows, smooth surfaces and anti-aliased edges. Off while you are iterating on camera angles, where it just costs you time.

**Cost.** Seconds to a minute per panel depending on pixel count and surface complexity.

### `render.orbital_isovalue`
**Orbital isosurface value** — e^(1/2) bohr^(-3/2) · type `float` · default `0.02` · range `0.001–`

The contour level at which molecular-orbital lobes are drawn. Both +value and -value surfaces are rendered, in two different colours.

**How to choose.** 0.02 is the near-universal convention in the literature and what most published HOMO/LUMO figures use; 0.05 gives tighter, smaller lobes and 0.01 gives diffuse, sprawling ones that can obscure the molecular skeleton. Whatever you pick, use the SAME value for every molecule in the figure or the comparison is meaningless - and state it in the caption, which the manuscript does not.

### `render.density_isovalue`
**Density isosurface value** — e/bohr^3 · type `float` · default `0.0004` · range `1e-05–`

The electron-density contour that defines the molecular surface onto which the electrostatic potential is painted for the MEP map.

**How to choose.** 0.0004 a.u. is Bader's definition of the molecular surface and encloses roughly 99% of the electron density - the standard choice for MEP maps and the one that corresponds to the van der Waals surface. 0.002 a.u. gives a tighter surface closer to the nuclei, with systematically larger |V| values. Must be identical across panels, and must be stated in the caption.

### `render.esp_range`
**MEP colour-scale limits** — a.u. · type `float` · default `0.05` · range `0.001–`

The electrostatic potential is coloured from -value (red, electron-rich) to +value (blue, electron-poor). 0.05 a.u. is about 31 kcal/mol.

**How to choose.** This is the single most abused setting in published MEP figures. If each panel is auto-scaled to its own min/max - which is what most GUIs do by default - then panels CANNOT be compared, because the same colour means a different potential in each one. Fixing one range across all molecules is what makes Figure 2C a real comparison. 0.05 a.u. suits neutral organics with hydroxyl and methoxy groups; widen it if your series includes charged species.

### `render.esp_scale`
**MEP scaling mode** · type `choice` · default `auto_shared`

How the electrostatic-potential colour range is chosen across the panels.

**How to choose.** A fixed range is right for a set of neutral molecules and is what makes panels comparable. It breaks for ions: the potential around a +1 cation is positive essentially everywhere and far outside +-0.05 a.u., so the panel saturates into a featureless blue blob and all the chemical detail is lost. 'auto_shared' measures the potential on each molecular surface and picks ONE range that covers them all - still a fair comparison, and it prints the range on the colour bar. 'per_panel' scales each molecule separately: use it only to show one molecule's own detail, never to compare, and the figure then labels each panel with its own range. 'fixed' uses render.esp_range exactly.

| value | name | what it means | preset |
|---|---|---|---|
| `auto_shared` | One range covering every panel, measured from the data | Comparable across panels and safe for mixed charges. | `recommended` `fast` |
| `fixed` | Exactly render.esp_range | Reproducible and comparable, but saturates for charged species. | `paper` |
| `per_panel` | Each molecule scaled to its own range | Maximum detail per molecule; panels are NOT comparable. |  |

### `render.esp_autoscale`
**Auto-scale MEP per molecule (legacy)** · type `bool` · default `False`

If true, each MEP panel is scaled to its own potential range instead of the fixed esp_range.

**How to choose.** Leave this OFF for any comparative figure. It exists only because reviewers sometimes ask to see a single molecule's full dynamic range. If you turn it on, the pipeline writes the per-panel range into the figure so the reader is not misled.

### `render.esp_palette`
**MEP colour palette** · type `choice` · default `paper`

The colour ramp used for the electrostatic potential surface.

**How to choose.** 'paper' reproduces the manuscript's blue-green-yellow-orange-red sequence. 'rwb' (red-white-blue) is the cleaner scientific choice because white marks the zero of the potential unambiguously, and it survives greyscale printing. Rainbow ramps are not perceptually uniform and exaggerate mid-range differences.

| value | name | what it means | preset |
|---|---|---|---|
| `paper` | Blue-green-yellow-orange-red | The manuscript's ramp. Red = negative potential (electron-rich), blue = positive. | `paper` |
| `rwb` | Red-white-blue | Diverging, perceptually sane, zero is unmistakable. The better default. | `recommended` |
| `rainbow` | Rainbow (PyMOL default) | Included for continuity with older figures. Not recommended. |  |

### `render.orbital_colors`
**Orbital lobe colours** · type `str` · default `firebrick,skyblue`

Two comma-separated colour names for the positive and negative phase of each orbital.

**How to choose.** Any PyMOL colour name works. Keep the two phases strongly contrasting; avoid red/green for colour-blind readers. The classic choices are red/blue or orange/cyan.

### `render.background`
**Background** · type `choice` · default `white`

Background of the rendered panels.

**How to choose.** 'transparent' is the most flexible for figure assembly - it composites onto any page colour. 'white' is safest if the journal flattens alpha channels.

| value | name | what it means | preset |
|---|---|---|---|
| `white` | White | Opaque white. Safe everywhere. | `paper` `fast` |
| `transparent` | Transparent | Alpha channel preserved; best for assembly. | `recommended` |
| `black` | Black | For slides and posters. |  |

### `render.orientation`
**Molecular orientation** · type `choice` · default `principal_plane`

How each molecule is rotated before rendering.

**How to choose.** 'principal_plane' puts the molecular plane in the screen plane using the inertia tensor, which is what you want for flat fused-ring systems - every panel ends up oriented the same way, so the reader compares chemistry rather than viewpoint. 'first_molecule' aligns everything onto the first molecule's frame, which is better for a congeneric series that shares a scaffold. 'none' keeps the ORCA orientation, which is essentially arbitrary.

| value | name | what it means | preset |
|---|---|---|---|
| `principal_plane` | Principal plane (inertia tensor) | Molecular plane parallel to the screen. Best for planar aromatics. | `paper` `recommended` `fast` |
| `first_molecule` | Align to first molecule | Least-squares fit of every molecule onto the first by maximum common substructure. Best for a congeneric series. |  |
| `none` | As calculated | Whatever orientation ORCA left. Not recommended for figures. |  |

### `render.stick_radius`  *(advanced)*
**Stick radius** — A · type `float` · default `0.12`

Bond cylinder radius in the ball-and-stick representation.

**How to choose.** 0.12-0.15 A reads well at figure scale. Thinner sticks disappear when the figure is reduced to a single column.

### `render.sphere_scale`  *(advanced)*
**Atom sphere scale** · type `float` · default `0.22`

Atomic sphere radius as a fraction of the van der Waals radius.

**How to choose.** 0.22 gives the conventional ball-and-stick look. 1.0 would give a space-filling model, which hides the orbital lobes.

### `render.surface_quality`
**Surface mesh quality** · type `int` · default `2` · range `0–4`

Tessellation density of the molecular and isodensity surfaces. Higher values subdivide the mesh further, so curved regions read as smooth rather than faceted.

**How to choose.** 2 is the right setting for a figure; 0-1 leave visible facets on small molecules, which is what makes a surface look blocky. 3 and above rarely changes anything you can see and makes ray tracing noticeably slower. If a surface still looks angular at 2, the problem is usually colour banding rather than geometry - raise render.surface_color_smoothing instead.

**Cost.** Each step roughly doubles the vertex count and the ray-trace time for that surface.

### `render.surface_color_smoothing`
**Surface colour smoothing** · type `float` · default `0.8` · range `0.0–1.0`

How far per-atom colours are blended across the molecular surface (PyMOL's surface_color_smoothing_threshold).

**How to choose.** This controls colour banding, not geometry. Each surface vertex takes its colour from the nearest atom, so with little smoothing you get hard-edged patches wherever two atoms of different value meet - the blocky look on a lipophilicity map. 0.30 blends them into a continuous gradient, which is the honest rendering of a fragment-based potential anyway, since the underlying quantity is a distance-weighted sum rather than a per-atom step function. Compared side by side on a carbazole: PyMOL's default of 0.05 and a value of 0.30 both leave visible angular wedges; 0.80 gives a continuous gradient while each feature is still clearly attributable to its atom; 1.00 begins to wash out small localized features. Lower it toward 0.05 if you want to see exactly which atom owns which patch. Affects only per-atom colouring (Figure 13), not the ramp-coloured MEP surface, which is already continuous.

### `render.surface_transparency`
**Surface transparency** · type `float` · default `0.25` · range `0.0–1.0`

Transparency of the MEP and lipophilicity surfaces, 0 = opaque, 1 = invisible.

**How to choose.** 0.2-0.3 lets the molecular skeleton show through the surface, which is what makes a MEP map interpretable - otherwise the reader cannot tell which atom a red patch belongs to. Set to 0 if you want the classic opaque surface.

### `render.show_hydrogens`
**Show hydrogens** · type `bool` · default `True`

Whether hydrogen atoms appear in the 3D panels.

**How to choose.** Keep them: the MEP panels of this manuscript specifically discuss the positive potential over particular hydrogens, and hiding them would make that unreadable.

### `render.pymol_path`  *(advanced)*
**PyMOL executable** · type `path` · default ``

Full path to the PyMOL binary. Empty means search PATH, then the directory of the Python interpreter running the pipeline.

**How to choose.** Only set this if you have several PyMOL installations and need a specific one - for example a licensed Incentive PyMOL rather than the open-source build.

### `render.vmd_path`  *(advanced)*
**VMD executable** · type `path` · default ``

Full path to the VMD binary. Only used when render.engine = 'vmd'.

**How to choose.** Empty means search PATH.

### `render.label_atoms`
**Atom-number labels** · type `str` · default ``

Comma-separated 1-based atom indices to label in the geometry panel, e.g. '3,4,9'.

**How to choose.** Use this to fix the manuscript's worst reproducibility problem - it quotes dihedrals over atoms C10, C12 and N11, which do not exist in carbazole numbering, and uses three mutually incompatible schemes across sections. Labelling the atoms you actually discuss removes all ambiguity. Leave empty for no labels.

## Figure assembly  (`figures`)
*Which figures to build and how they are typeset.*

### `figures.which`
**Figures to produce** · type `multichoice` · default `['fig2', 'fig12', 'fig13', 'fig14']`

Selects which of the manuscript-style figures are assembled.

**How to choose.** fig2 needs the full DFT chain. fig12, fig13 and fig14 need only the optimized geometries plus the properties CSV, so they are quick to iterate on once the DFT has run once.

| value | name | what it means | preset |
|---|---|---|---|
| `fig2` | Figure 2 - geometry / HOMO-LUMO / MEP | Three rows: optimized structures, frontier orbitals with the gap, MEP surfaces. Requires DFT, cubes and rendering. | `paper` `recommended` `fast` |
| `fig12` | Figure 12 - pKa and net charge vs pH | Acid-base analysis: annotated 2D structures plus a speciation curve. Needs pKa values in the properties CSV. | `paper` `recommended` |
| `fig13` | Figure 13 - lipophilicity surfaces | Molecular lipophilicity potential painted on the molecular surface, with TPSA. Requires optimized geometries and rendering. | `paper` `recommended` |
| `fig14` | Figure 14 - CNS MPO and metabolism maps | logD/TPSA bubble plot with MPO scoring, plus site-of-metabolism sensitivity maps. Needs the properties and SOM CSVs. | `paper` `recommended` |

### `figures.reference_solvent`
**Solvent shown in the figures** · type `choice` · default `water`

Which solvent's results are used for the figure panels when several were calculated.

**How to choose.** The manuscript shows water in Figure 2 and relegates the other solvents to the supplement. Must be one of the solvents you actually calculated.

| value | name | what it means | preset |
|---|---|---|---|
| `water` | Water | The manuscript's choice for Figure 2. | `paper` `recommended` `fast` |
| `gas` | Gas phase | The only medium in which Koopmans descriptors are physically interpretable. |  |
| `methanol` | Methanol |  |  |
| `chloroform` | Chloroform |  |  |
| `dmso` | DMSO |  |  |

### `figures.dpi`
**Output resolution** — dpi · type `int` · default `300` · range `72–`

Raster resolution of the assembled figures.

**How to choose.** 300 dpi is the floor for most journals; 600 for line art if the publisher asks. Vector formats ignore this for the plot elements but still use it for the embedded 3D renders.

### `figures.formats`
**Output formats** · type `multichoice` · default `['png', 'pdf']`

File formats written for each figure.

**How to choose.** PNG for looking at and for submission systems that demand raster; PDF for vector text and lines that stay sharp at any zoom. TIFF only if a journal insists.

| value | name | what it means | preset |
|---|---|---|---|
| `png` | PNG | Raster, universally accepted. | `paper` `recommended` `fast` |
| `pdf` | PDF | Vector plot elements, embedded raster renders. | `recommended` |
| `svg` | SVG | Editable vector - open in Inkscape to tweak labels by hand. |  |
| `tiff` | TIFF | Large files; some journals require it. |  |

### `figures.panel_label_style`
**Panel letter style** · type `choice` · default `upper`

How panel letters are drawn in the assembled figures.

**How to choose.** Match your target journal. Molecular Diversity uses bare capitals (A, B, C), which is the default here.

| value | name | what it means | preset |
|---|---|---|---|
| `upper` | A B C | Bare capitals, as in the manuscript. | `paper` `recommended` `fast` |
| `lower` | a b c | Lower case. |  |
| `paren_upper` | (A) (B) (C) | Parenthesized capitals. |  |
| `none` | No letters | For a single-panel figure. |  |

### `figures.ph_range`
**pH range for Figure 12B** · type `str` · default `0,14`

Minimum and maximum pH of the speciation curve, comma separated.

**How to choose.** 0-14 shows the full picture. Narrow to 4,10 if you only care about the physiologically relevant window and want the transitions to fill the axes.

### `figures.physiological_ph`
**Physiological pH marker** · type `float` · default `7.4`

The pH at which the predominant species is annotated in Figure 12B.

**How to choose.** 7.4 is blood plasma. Use 6.5 for the small-intestine lumen or 1.5 for gastric fluid if you are arguing about oral absorption.

### `figures.mlp_palette`
**Lipophilicity palette** · type `choice` · default `paper`

Colour ramp for the molecular lipophilicity potential surface in Figure 13.

**How to choose.** 'paper' matches the manuscript's caption exactly - red for the least lipophilic fragments through to blue for the most lipophilic. Note this is the OPPOSITE of the convention most medicinal chemists expect (blue usually means polar), so whichever you choose, the caption has to say so. 'conventional' uses blue = hydrophilic, orange = lipophilic.

| value | name | what it means | preset |
|---|---|---|---|
| `paper` | Red (polar) to blue (lipophilic) | The manuscript's stated ramp. | `paper` |
| `conventional` | Blue (polar) to orange (lipophilic) | The more widely used convention in medicinal chemistry. | `recommended` |

### `figures.mlp_transparency`
**Lipophilicity surface transparency** · type `float` · default `0.35` · range `0.0–1.0`

Transparency of the Figure 13 surfaces: 0 = opaque, 1 = invisible.

**How to choose.** A fully opaque surface hides the molecule, so the reader cannot tell which atom a coloured patch belongs to - the same objection that applies to an unlabelled MEP map. 0.3-0.4 lets the skeleton read through while keeping the colours saturated. Push past 0.5 and the colour washes out; set 0 only if you specifically want the classic solid surface.

### `figures.mlp_show_structure`
**Show the molecule inside the surface** · type `bool` · default `True`

Draws the optimized geometry as thin sticks underneath the translucent lipophilicity surface.

**How to choose.** Keep this on. It is what makes the panel interpretable: a blue patch means nothing until you can see it sits over the chlorine rather than the methoxy. Turn it off only for a purely qualitative shape comparison.

### `figures.mpo_variant`
**CNS MPO implementation** · type `choice` · default `wager2010`

Which published set of desirability functions is used to score CNS drug-likeness.

**How to choose.** wager2010 is the original six-parameter CNS MPO (clogP, clogD, MW, TPSA, HBD, most-basic pKa), each mapped to a 0-1 desirability and summed to a 0-6 score. Note that the manuscript's Methods state the HBD threshold as <=1, whereas Wager's published desirability function for HBD is centred on <=0.5 - the pipeline uses the published function and will therefore differ slightly from the paper's Table 3.

| value | name | what it means | preset |
|---|---|---|---|
| `wager2010` | Wager et al. 2010 (6 parameters) | The original and still standard CNS MPO. | `paper` `recommended` `fast` |

### `figures.mpo_binary_verdict`
**Show accept/reject verdict** · type `bool` · default `False`

Whether to print a pass/fail label alongside the CNS MPO score.

**How to choose.** Leave OFF. CNS MPO is a continuous desirability score and its authors explicitly propose MPO >= 4 as an enrichment criterion, not a filter. The manuscript's Table 3 binarizes it into 'Accepted'/'Rejected', which puts a score of 3.76 on the wrong side of a boundary the method does not have. Turn it on only if a co-author insists.

### `figures.properties_csv`
**Physicochemical properties file** · type `path` · default `data/properties.csv`

Optional CSV of measured or externally predicted values. Anything in it overrides what the pipeline works out for itself.

**How to choose.** You do not need this file. Give the pipeline SDF structures and it computes MW, TPSA, HBD and Crippen logP from the connection table, predicts pKa from the SMARTS ionization model, and derives logD at pH 7.4 from the two. Supply the CSV when you have measured pKa or logD, or when you want permeability, clearance, plasma protein binding or P-gp on the figures - those need a trained ADMET model (ADMETlab 3.0, Deep-PK, your own) and the pipeline will not invent them. Columns are matched by name; blanks are filled, not overwritten.

### `figures.pka_source`
**Where pKa comes from** · type `choice` · default `predict`

How the acidic and basic pKa values behind Figure 12 and logD are obtained.

**How to choose.** 'predict' needs nothing from you: a curated SMARTS table assigns a class average pKa to every ionizable group it recognises - carboxylic acid, phenol, aliphatic amine, pyridine, imidazole, tetrazole and about thirty more - and reports the groups it found alongside the number. Expect +-1 to 1.5 log units, which settles whether a compound is charged at pH 7.4 but is not a value to quote in a paper. 'csv' uses only what you supply and leaves the rest blank. 'none' switches ionization off entirely, which draws every compound at the formal charge of the structure you gave.

| value | name | what it means | preset |
|---|---|---|---|
| `predict` | Predict from structure | SMARTS class averages; fills any pKa the CSV does not give. | `recommended` `fast` |
| `csv` | Only what the properties CSV gives | Nothing is inferred. Panels without a value say so. | `paper` |
| `none` | No ionization model | Every compound is drawn at its formal charge. |  |

### `figures.logd_source`
**Where logD(7.4) comes from** · type `choice` · default `derive`

How the buffer lipophilicity on Figure 14 is obtained.

**How to choose.** 'derive' computes logD from Crippen logP and the ionization model: only the neutral microspecies is taken to partition, which is the standard approximation. It is a lower bound for strong bases, where ion pairing adds roughly one unit, and it is a floor rather than a prediction for a permanently charged compound, which has no neutral form at all. 'csv' uses only supplied values. 'none' leaves the axis empty unless the CSV fills it.

| value | name | what it means | preset |
|---|---|---|---|
| `derive` | Derive from clogP and the predicted pKa | Self-contained; labelled on the figure as derived. | `recommended` `fast` |
| `csv` | Only what the properties CSV gives | For measured logD or an external ADMET prediction. | `paper` |
| `none` | No logD | Figure 14 panel A falls back to clogP. |  |

### `figures.som_source`
**Site-of-metabolism source** · type `choice` · default `fukui_radical`

Where the per-atom metabolic-lability weights for Figure 14 B-E come from.

**How to choose.** 'fukui_radical' uses the condensed radical Fukui index f0 that this pipeline already computes - susceptibility to radical attack is a recognized quantum-chemical proxy for the hydrogen abstraction and aromatic oxidation steps that CYP450 performs, so the panel is fully reproducible from your own DFT with no external service. It is a reactivity proxy, NOT a trained site-of-metabolism model: it knows nothing about how the substrate binds in the enzyme active site, which is often what decides the real site. 'csv' reads measured or predicted probabilities you export from XenoSite or ADMETlab; use that when you have them and say in the caption which source you used. 'both' prefers the CSV and silently falls back to Fukui where a file is missing.

| value | name | what it means | preset |
|---|---|---|---|
| `fukui_radical` | Radical Fukui index f0 (from this pipeline's DFT) | Self-contained and reproducible; a reactivity proxy, blind to binding mode. | `recommended` `fast` |
| `csv` | CSV export from XenoSite / ADMETlab | A trained site-of-metabolism model. What the manuscript used. | `paper` |
| `both` | CSV where available, Fukui otherwise | Prefer the trained model, fall back to the DFT proxy. |  |

### `figures.som_dir`
**Site-of-metabolism directory** · type `path` · default `data/som`

Directory holding <molecule>.csv files mapping atom index to metabolism probability, used for Figure 14 panels B-E.

**How to choose.** Export these from XenoSite or ADMETlab's SOM module. Format: two columns, atom_index (1-based, matching the SMILES atom order) and probability (0-1). If a file is missing, that panel is drawn greyed out with a note rather than silently omitted.

## Execution  (`execution`)
*Where ORCA lives, how much hardware to use, and restart behaviour.*

### `execution.orca_path`
**ORCA executable** · type `path` · default `/home/moazhosnytaha/orca_6_1_0/orca`

Full path to the ORCA binary. The helper tools orca_plot and orca_vpot are located next to it automatically.

**How to choose.** Must be the absolute path, not just 'orca' - ORCA requires this for its own parallel startup. Verified present at the default on this machine.

### `execution.nprocs`
**Parallel processes** — cores · type `int` · default `4` · range `1–`

Number of MPI processes ORCA uses per job (the %pal block).

**How to choose.** Leave 2-4 cores free for the operating system and for the figure rendering that runs alongside. Hybrid DFT with RIJCOSX scales reasonably to about 8-16 cores for a 30-atom molecule; beyond that the communication overhead eats the gain. Running two molecules on 6 cores each usually beats one molecule on 12.

### `execution.maxcore_mb`
**Memory per process** — MB · type `int` · default `900` · range `200–`

The %maxcore value - memory ORCA may use PER MPI PROCESS, not in total.

**How to choose.** This is the single most common cause of ORCA crashes, because it is memory PER MPI PROCESS, not per job. Real usage is maxcore_mb x nprocs x max_parallel_jobs - so 3000 MB on 8 cores with 2 concurrent jobs asks for 48 GB, not 3 GB. Set it to roughly 75% of (usable RAM / (nprocs x max_parallel_jobs)); ORCA treats the value as a target rather than a hard cap and will overshoot. `dft validate` checks the product against the RAM this machine actually has and warns if it does not fit.

### `execution.scratch_dir`
**Scratch directory** · type `path` · default `scratch`

Working directory for ORCA's large temporary integral files.

**How to choose.** Put it on the fastest local disk you have, ideally an SSD, and never on a network mount - ORCA's integral files are written and reread constantly and NFS will dominate your runtime. Needs several GB free per concurrent job.

### `execution.resume`
**Skip completed jobs** · type `bool` · default `True`

If a job's output already exists and terminated normally, reuse it instead of recalculating.

**How to choose.** Keep on. It makes the pipeline restartable after an interruption and lets you iterate on the figures without repeating hours of DFT. The check is for ORCA's own 'TERMINATED NORMALLY' marker, so a crashed or truncated job is correctly rerun.

### `execution.max_parallel_jobs`
**Concurrent ORCA jobs** — jobs · type `int` · default `2` · range `1–`

How many ORCA calculations run at the same time, each using nprocs cores.

**How to choose.** Total core usage is max_parallel_jobs x nprocs, which must not exceed your core count. For a set of small molecules, 3 jobs x 4 cores finishes sooner than 1 job x 12 cores because ORCA's parallel efficiency falls off. Watch memory: each job separately claims maxcore_mb x nprocs.

### `execution.timeout_min`  *(advanced)*
**Per-job timeout** — minutes · type `int` · default `480`

Wall-clock limit after which a single ORCA job is killed and marked failed.

**How to choose.** A generous ceiling that catches genuinely stuck jobs - an optimization oscillating between two geometries will otherwise run forever. 480 min is ample for a 30-atom molecule; raise it for anything over 60 atoms.

### `execution.work_dir`
**Output directory** · type `path` · default `runs`

Root directory for all generated inputs, outputs, cubes, figures and logs.

**How to choose.** One subdirectory per molecule per solvent is created underneath. Keep separate work directories for separate projects so the resume logic cannot reuse the wrong output.

### `execution.dry_run`
**Dry run** · type `bool` · default `False`

Generate every input file and script but do not launch ORCA.

**How to choose.** The fastest way to check your settings before committing hours of compute: inspect the generated .inp files, confirm the keyword line is what you intended, then rerun with this off. Highly recommended before any long campaign.
