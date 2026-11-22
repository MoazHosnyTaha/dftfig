"""
Parameter registry — the single source of truth for every knob in the pipeline.

Each `Param` carries, besides its default:

  * ``help``     - what the parameter actually does, in plain language.
  * ``guidance`` - how to decide, and when the default is the wrong answer.
  * ``cost``     - the effect on wall-clock time or memory, so the trade-off is visible.
  * ``choices``  - for categorical parameters, every allowed value with its own explanation.

The interactive wizard renders exactly this information, so documentation and behaviour
cannot drift apart. Add a parameter here and it appears in the wizard, the YAML config,
and ``dftfig params`` automatically.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

# Tags mark choices that belong to a named preset.
TAG_PAPER = "paper"            # reproduces the chlocarbazomycin manuscript
TAG_RECOMMENDED = "best"       # what I would actually run
TAG_FAST = "fast"              # for testing that the plumbing works


@dataclass
class Choice:
    value: Any
    label: str
    detail: str
    tags: tuple = ()


@dataclass
class Param:
    key: str                       # dotted path into the config dict, e.g. "theory.functional"
    label: str
    kind: str                      # choice | multichoice | str | int | float | bool | path | kvlist
    default: Any
    help: str
    guidance: str = ""
    unit: str = ""
    cost: str = ""
    choices: list = field(default_factory=list)
    advanced: bool = False
    minimum: float | None = None
    maximum: float | None = None


@dataclass
class Group:
    name: str
    title: str
    blurb: str
    params: list


# ----------------------------------------------------------------------------- molecules
#: The chlocarbazomycins, reconstructed from Figure 1 of Mol. Divers. 2026, 30:5593.
#: Loaded only by the `paper` preset or `dft init --demo` - not a default.
REFERENCE_SET = {
    "CCB1": "COc1ccc2[nH]c3ccccc3c2c1Cl",
    "CCB2": "COc1ccc2[nH]c3ccc(O)cc3c2c1Cl",
    "CCB3": "COc1ccc2[nH]c3c(O)cccc3c2c1Cl",
    "CCB4": "COn1c2ccccc2c2c1ccc(OC)c2Cl",
}

G_MOLECULES = Group(
    "molecules", "Molecules",
    "What to calculate. One SDF per compound is the normal route; SMILES also works.",
    [
        Param(
            "molecules.structure_source", "Where structures come from", "choice", "sdf",
            help="The input format for your compounds.",
            guidance="'sdf' is the normal route: drop one .sdf per compound into a folder and "
                     "the pipeline reads them all, taking each compound's charge and 3D geometry "
                     "from the file. 'smiles' is convenient for a handful of compounds typed by "
                     "hand. 'xyz' is for geometries you already trust and do not want a conformer "
                     "search to disturb - it still needs a SMILES per molecule for bond orders.",
            choices=[
                Choice("sdf", "SDF files (one per compound)",
                       "Reads every .sdf in molecules.sdf_dir. Charge is taken from the formal "
                       "charges in the file, and 3D coordinates are used directly when present.",
                       (TAG_RECOMMENDED, TAG_FAST)),
                Choice("smiles", "SMILES typed in the config",
                       "Name-to-SMILES entries under molecules.definitions. A conformer is "
                       "generated for each.", (TAG_PAPER,)),
                Choice("xyz", "Existing .xyz geometries",
                       "Reads <name>.xyz from molecules.xyz_dir and skips the conformer search. "
                       "Still requires molecules.definitions for the chemistry perception."),
            ],
        ),
        Param(
            "molecules.sdf_dir", "Folder of SDF files", "path", "molecules",
            help="Directory searched for compound files when structure_source = 'sdf'. "
                 "Relative paths are resolved from the folder you run `dft` in.",
            guidance="One .sdf per compound is cleanest - the file name becomes the molecule "
                     "name, which becomes the output directory and the figure panel label. A "
                     "multi-record SDF also works: each record becomes its own molecule. Export "
                     "from Avogadro, ChemDraw, MOE or RDKit; keep 3D coordinates if you have "
                     "them and the pipeline will start the optimization from them.",
        ),
        Param(
            "molecules.sdf_pattern", "File pattern", "str", "*.sdf",
            help="Glob matched inside sdf_dir.",
            guidance="Widen to '*.sd*' to pick up .sd and .sdf together, or narrow to a prefix "
                     "like 'series1_*.sdf' to run one subset without moving files around.",
        ),
        Param(
            "molecules.name_from", "Molecule naming", "choice", "filename",
            help="Where each molecule's name comes from.",
            guidance="The name becomes a directory name, a filename prefix and a figure panel "
                     "label, so short and filesystem-safe is best. 'filename' is predictable; "
                     "'sdf_title' uses the title line inside the file, which is handy when the "
                     "files came from a database export with meaningless names. Characters that "
                     "are awkward in paths are replaced with underscores either way.",
            choices=[
                Choice("filename", "The .sdf file name",
                       "compound_17.sdf -> compound_17.", (TAG_RECOMMENDED, TAG_FAST)),
                Choice("sdf_title", "The title line inside the file",
                       "The first line of the SDF record, when it is not empty."),
            ],
        ),
        Param(
            "molecules.only", "Run only these compounds", "str", "",
            help="Comma-separated compound names. Empty means every structure found.",
            guidance="Names are the SDF file names without the extension - the same names "
                     "that appear in runs/molecules.csv and on the figure panels. Use this "
                     "to re-run one compound that failed without touching the rest, or to "
                     "produce a figure from a subset of a large batch. `dft run` sets it "
                     "for you when you deselect compounds in the list it shows.",
        ),
        Param(
            "molecules.limit", "How many compounds to run", "int", 0, unit="compounds",
            help="Run only the first N compounds found. 0 means all of them.",
            guidance="Set this before a long campaign so you know exactly what is queued: with "
                     "40 SDF files in a folder, limit = 3 runs three of them and tells you which. "
                     "Use it to time one compound before committing to the batch - `dft validate` "
                     "reports the resulting job count. Files are taken in sorted order, so the "
                     "selection is reproducible.",
            minimum=0,
        ),
        Param(
            "molecules.charge_mode", "Charge handling", "choice", "auto",
            help="How each compound's net charge is decided.",
            guidance="'auto' reads the formal charges in the structure, per compound. This is "
                     "what you want whenever a set is not uniformly neutral - a quaternary "
                     "ammonium in the same batch as a neutral nucleoside has a different electron "
                     "count, and running it at charge 0 computes a different species entirely "
                     "(wrong electron count, wrong SCF, meaningless orbitals). 'fixed' forces one "
                     "charge on every compound and is only correct for a uniform set. The "
                     "detected charge of every molecule is printed before any job starts and "
                     "written to molecules.csv, so an unexpected value is visible immediately.",
            choices=[
                Choice("auto", "Detect per compound",
                       "Sum of formal charges in the SMILES or SDF. Different compounds in one "
                       "run may carry different charges.", (TAG_RECOMMENDED, TAG_FAST)),
                Choice("fixed", "Force molecules.charge on all of them",
                       "Every compound runs at the same charge. Correct only for a uniform set.",
                       (TAG_PAPER,)),
            ],
        ),
        Param(
            "molecules.charge", "Fixed charge", "int", 0,
            help="The charge applied to every compound when charge_mode = 'fixed'. Ignored in "
                 "'auto' mode.",
            guidance="Use this to force a charge state the structure does not encode - for "
                     "example to study the anion of a neutral acid without editing the input "
                     "file. Remember it then applies to the whole set.",
        ),
        Param(
            "molecules.multiplicity_mode", "Spin multiplicity handling", "choice", "auto",
            help="How each compound's spin multiplicity is decided.",
            guidance="'auto' picks the lowest spin state consistent with the electron count: a "
                     "singlet for an even number of electrons, a doublet for an odd number. That "
                     "is right for closed-shell organics and for radical ions. Use 'fixed' only "
                     "when you specifically want a higher spin state, such as a triplet; the "
                     "pipeline rejects a multiplicity that is impossible for the electron count "
                     "rather than letting ORCA fail halfway through.",
            choices=[
                Choice("auto", "Lowest spin state for the electron count",
                       "Singlet if even, doublet if odd.", (TAG_RECOMMENDED, TAG_FAST)),
                Choice("fixed", "Force molecules.multiplicity",
                       "For triplets and other deliberately high-spin states.", (TAG_PAPER,)),
            ],
        ),
        Param(
            "molecules.multiplicity", "Fixed multiplicity", "int", 1,
            help="2S+1 applied to every compound when multiplicity_mode = 'fixed'. "
                 "1 = closed-shell singlet, 2 = doublet, 3 = triplet.",
            guidance="Anything other than 1 switches ORCA to unrestricted DFT (UKS), which "
                     "doubles the orbital and population output and makes the frontier-orbital "
                     "panels ambiguous (alpha and beta sets differ).",
            minimum=1,
        ),
        Param(
            "molecules.definitions", "Molecule definitions (SMILES)", "kvlist", {},
            help="Name-to-SMILES map, used when structure_source is 'smiles' or 'xyz'. Ignored "
                 "for 'sdf'.",
            guidance="Empty by default - this pipeline has no built-in compounds. Add your own "
                     "as NAME: SMILES, or switch to SDF input. Panel order in every figure "
                     "follows the order given here. Keep names short: they become directory "
                     "names and figure labels. For the chlocarbazomycin reference set from the "
                     "manuscript, use the `paper` preset or run `dft init --demo`.",
        ),
        Param(
            "molecules.use_sdf_geometry", "Use 3D coordinates from the file", "bool", True,
            help="When an SDF record already contains a 3D conformer, start the DFT optimization "
                 "from it instead of generating a new one.",
            guidance="Keep this on when your SDFs came from a docking pose, a crystal structure "
                     "or a careful manual build - those geometries encode information a "
                     "conformer generator cannot recover. Turn it off to force a fresh conformer "
                     "search for every compound, which is the right choice if the SDFs were "
                     "produced from 2D depictions or you do not trust their provenance. Records "
                     "with no 3D coordinates always get a generated conformer regardless.",
        ),
        Param(
            "molecules.xyz_dir", "Directory of .xyz files", "path", "input_xyz",
            help="Only read when structure_source = 'xyz'. Expects one <molecule-name>.xyz "
                 "per entry in molecules.definitions.",
            guidance="Plain XYZ format (atom count, comment line, then element plus Cartesian "
                     "coordinates in Angstrom). Avogadro and GaussView both export this.",
            advanced=True,
        ),
    ],
)

# ---------------------------------------------------------------------------- conformers
G_CONFORMERS = Group(
    "conformers", "Conformer search", "Finding a sensible starting structure before DFT.",
    [
        Param(
            "conformers.n_confs", "Conformers to generate", "int", 50, unit="structures",
            help="How many trial 3D structures RDKit embeds before force-field ranking. The lowest-energy "
                 "survivor is what DFT then optimizes.",
            guidance="For rigid fused-ring systems like carbazoles only the methoxy rotation and the "
                     "hydroxyl orientation matter, so 20-50 is plenty. Push to 300+ for molecules with "
                     "more than about 5 rotatable bonds, otherwise you risk starting DFT from a geometry "
                     "that is not the global minimum - a far more expensive mistake than the conformer "
                     "search itself.",
            cost="Seconds. Never the bottleneck.",
            minimum=1,
        ),
        Param(
            "conformers.forcefield", "Force field for pre-optimization", "choice", "MMFF94s",
            help="The classical force field used to relax and rank the embedded conformers.",
            guidance="MMFF94s is the right default for drug-like organics. Switch to UFF only if MMFF "
                     "refuses to parameterize your molecule (unusual elements, exotic valences) - UFF "
                     "geometries are noticeably worse but DFT will fix them.",
            choices=[
                Choice("MMFF94s", "MMFF94s",
                       "Static variant of MMFF94 with planarized sp2 nitrogens. Best general choice for "
                       "aromatic heterocycles - it will keep the carbazole NH planar instead of pyramidalizing it.",
                       (TAG_PAPER, TAG_RECOMMENDED, TAG_FAST)),
                Choice("MMFF94", "MMFF94",
                       "Original Merck force field. Slightly different nitrogen out-of-plane treatment."),
                Choice("UFF", "UFF",
                       "Universal Force Field. Covers the whole periodic table but with poor accuracy for "
                       "organics. Fallback only."),
            ],
        ),
        Param(
            "conformers.prune_rms", "Duplicate pruning RMSD", "float", 0.5, unit="A",
            help="Conformers closer than this heavy-atom RMSD are treated as the same structure and "
                 "collapsed to one.",
            guidance="0.5 A removes trivial duplicates without discarding genuinely different rotamers. "
                     "Lower it toward 0.1 if you want to inspect fine differences; raise it to 1.0 to "
                     "keep only clearly distinct conformers.",
            minimum=0.0,
        ),
        Param(
            "conformers.random_seed", "Random seed", "int", 42,
            help="Seed for the ETKDG embedding. Fixing it makes the whole pipeline bit-reproducible.",
            guidance="Keep it fixed for anything you intend to publish, so a reviewer rerunning the "
                     "pipeline gets your exact geometries. Change it only to test that your result is "
                     "not an artefact of one particular embedding.",
        ),
        Param(
            "conformers.ff_max_iters", "Force-field iterations", "int", 2000,
            help="Maximum optimization steps per conformer at the force-field level.",
            guidance="2000 converges essentially any drug-like molecule. Raise only if the log reports "
                     "unconverged conformers.",
            advanced=True,
        ),
    ],
)

# -------------------------------------------------------------------------------- theory
G_THEORY = Group(
    "theory", "Level of theory (geometry optimization)",
    "The functional, basis set and convergence settings for the OPT+FREQ step.",
    [
        Param(
            "theory.functional", "Exchange-correlation functional", "choice", "B3LYP",
            help="The DFT approximation. Controls accuracy of geometries, orbital energies and every "
                 "descriptor derived from them.",
            guidance="B3LYP reproduces the manuscript. Be aware that bare B3LYP has NO London dispersion, "
                     "which matters because the paper's structural argument is about pi-stacking - see "
                     "theory.dispersion below. If you are starting fresh rather than reproducing, "
                     "wB97X-D4 or r2SCAN-3c are better value for the same money.",
            cost="Hybrids (B3LYP, PBE0, M06-2X) need exact exchange and are ~3-5x a pure GGA. "
                 "Range-separated hybrids (wB97X) add another ~30%.",
            choices=[
                Choice("B3LYP", "B3LYP",
                       "The functional used in the manuscript. Reliable ground-state geometries for "
                       "organics; underestimates barriers; no dispersion of its own.", (TAG_PAPER,)),
                Choice("B3LYP D3BJ", "B3LYP-D3BJ",
                       "B3LYP plus Grimme's D3 dispersion with Becke-Johnson damping. Essentially free, "
                       "and strictly better than bare B3LYP for anything involving stacking or "
                       "intramolecular contacts.", ()),
                Choice("PBE0 D3BJ", "PBE0-D3BJ",
                       "Parameter-free hybrid, 25% exact exchange. Often better orbital energies and gaps "
                       "than B3LYP."),
                Choice("wB97X-D4", "wB97X-D4",
                       "Range-separated hybrid with D4 dispersion. The most reliable general-purpose "
                       "choice today for geometries, gaps and non-covalent interaction energies.",
                       (TAG_RECOMMENDED,)),
                Choice("M062X", "M06-2X",
                       "Heavily parameterized meta-hybrid, good for main-group thermochemistry and "
                       "non-covalent complexes. Needs a dense integration grid (DefGrid3)."),
                Choice("r2SCAN-3c", "r2SCAN-3c",
                       "Composite method: meta-GGA with its own fitted triple-zeta basis and dispersion "
                       "baked in. Very fast and surprisingly accurate - ignores the basis_opt setting.",
                       (TAG_FAST,)),
            ],
        ),
        Param(
            "theory.basis_opt", "Basis set (optimization)", "choice", "6-31+G(d,p)",
            help="The atomic-orbital basis used for the geometry optimization and frequencies.",
            guidance="6-31+G(d,p) is the manuscript's choice: split-valence double-zeta, polarization on "
                     "all atoms, diffuse functions on heavy atoms. The diffuse set is what makes it slow, "
                     "and it is what you need for anions and for LUMO energies to mean anything. "
                     "def2-SVP is faster and usually gives indistinguishable geometries.",
            cost="Scales steeply: def2-SVP ~ 1x, 6-31+G(d,p) ~ 2-3x, def2-TZVP ~ 6-10x for the same molecule.",
            choices=[
                Choice("6-31+G(d,p)", "6-31+G(d,p)",
                       "The manuscript's basis. Diffuse functions on heavy atoms (+), polarization on "
                       "heavy atoms and hydrogens (d,p).", (TAG_PAPER,)),
                Choice("6-31G(d)", "6-31G(d)",
                       "No diffuse functions, no hydrogen polarization. Noticeably cheaper; fine for "
                       "geometries, poor for LUMO energies and anions."),
                Choice("def2-SVP", "def2-SVP",
                       "Karlsruhe double-zeta. Modern, well-balanced, pairs with the def2/J auxiliary "
                       "basis for RIJCOSX. Good geometry-only choice.", (TAG_FAST,)),
                Choice("def2-TZVP", "def2-TZVP",
                       "Triple-zeta. Use when orbital energies or energetics - not just geometry - are "
                       "the deliverable.", (TAG_RECOMMENDED,)),
                Choice("ma-def2-SVP", "ma-def2-SVP",
                       "Minimally augmented def2-SVP: diffuse functions at a fraction of the cost of a "
                       "full aug- set. A good substitute for 6-31+G(d,p)."),
            ],
        ),
        Param(
            "theory.dispersion", "Dispersion correction", "choice", "none",
            help="An additive London-dispersion term. Density functionals without it systematically miss "
                 "attraction between non-bonded fragments.",
            guidance="Set to 'none' to reproduce the manuscript exactly. Set to D3BJ or D4 for any new "
                     "work: it costs essentially nothing and the paper's own pi-stacking narrative needs "
                     "it. Leave at 'none' if the functional already includes dispersion by construction "
                     "(wB97X-D4, r2SCAN-3c) - adding it twice is a real error.",
            cost="Negligible - a pairwise atom-atom sum, microseconds.",
            choices=[
                Choice("none", "None",
                       "No dispersion. What the manuscript did. Acceptable for covalent geometries of a "
                       "single rigid molecule; not acceptable for stacking or conformer energetics.", (TAG_PAPER,)),
                Choice("D3BJ", "D3(BJ)",
                       "Grimme's D3 with Becke-Johnson damping. The standard pairing for B3LYP and PBE0.",
                       (TAG_RECOMMENDED, TAG_FAST)),
                Choice("D4", "D4",
                       "Newer, charge-dependent dispersion. Slightly better, slightly slower, needs "
                       "parameters for your functional."),
            ],
        ),
        Param(
            "theory.ri_approx", "Integral approximation", "choice", "RIJCOSX",
            help="How the Coulomb (J) and exchange (K) integrals are evaluated. RI = resolution of the "
                 "identity, an auxiliary-basis expansion; COSX = chain-of-spheres numerical exchange.",
            guidance="RIJCOSX makes hybrid DFT 5-20x faster with errors of ~0.1 kcal/mol on relative "
                     "energies and under 0.001 eV on orbital energies - invisible at the precision anyone "
                     "reports. Set to 'none' only if you want integrals handled exactly as Gaussian 09 "
                     "would have, for a strict numerical comparison with the manuscript.",
            cost="RIJCOSX typically cuts hybrid-DFT SCF time by 5-20x. 'none' is the slow, exact path.",
            choices=[
                Choice("RIJCOSX", "RIJCOSX",
                       "RI for Coulomb, chain-of-spheres for exchange. The standard production setting "
                       "for hybrid DFT in ORCA. Requires an auxiliary basis (see theory.aux_basis).",
                       (TAG_RECOMMENDED, TAG_FAST)),
                Choice("none", "None (exact, conventional)",
                       "Four-index integrals computed exactly, as Gaussian 09 does by default. Slow, but "
                       "removes any approximation when comparing numbers against the published values.", (TAG_PAPER,)),
                Choice("RIJK", "RIJK",
                       "RI for both J and K. More accurate than COSX for exchange but needs a /JK "
                       "auxiliary basis and more memory."),
            ],
        ),
        Param(
            "theory.aux_basis", "Auxiliary basis", "str", "def2/J",
            help="The fitting basis used by RIJCOSX/RIJK. Ignored when ri_approx = none.",
            guidance="def2/J works with every Karlsruhe and Pople basis for main-group elements and is "
                     "the right answer in almost all cases. Use def2/JK if you selected RIJK.",
            advanced=True,
        ),
        Param(
            "theory.scf_conv", "SCF convergence", "choice", "TightSCF",
            help="How tightly the self-consistent field is converged at each geometry step.",
            guidance="TightSCF is the sane default and is required for reliable frequencies. NormalSCF "
                     "can leave enough noise in the gradient to stall a geometry optimization. "
                     "VeryTightSCF only for numerical frequencies or very flat potential surfaces.",
            cost="Each tier is roughly 10-25% more SCF iterations than the one below.",
            choices=[
                Choice("TightSCF", "TightSCF",
                       "Energy change < 1e-8 Eh. Required for clean analytic frequencies.",
                       (TAG_PAPER, TAG_RECOMMENDED)),
                Choice("NormalSCF", "NormalSCF", "Energy change < 1e-6 Eh. Fine for a quick geometry.", (TAG_FAST,)),
                Choice("VeryTightSCF", "VeryTightSCF",
                       "Energy change < 1e-9 Eh. For numerical gradients/Hessians and flat surfaces."),
            ],
        ),
        Param(
            "theory.grid", "DFT integration grid", "choice", "DefGrid2",
            help="The numerical grid on which the exchange-correlation functional is integrated.",
            guidance="DefGrid2 is ORCA 6's default and is adequate for B3LYP. Move to DefGrid3 for "
                     "Minnesota functionals (M06-2X and relatives), which are notoriously grid-sensitive, "
                     "or if frequencies come out with small spurious imaginary modes.",
            cost="DefGrid3 is roughly 2-3x DefGrid2 in the XC build step.",
            choices=[
                Choice("DefGrid2", "DefGrid2", "ORCA 6 default. Fine for B3LYP, PBE0, wB97X.",
                       (TAG_PAPER, TAG_RECOMMENDED, TAG_FAST)),
                Choice("DefGrid1", "DefGrid1", "Coarse. Pre-screening only - do not report frequencies from it."),
                Choice("DefGrid3", "DefGrid3", "Dense. Required for M06-2X; cures spurious tiny imaginary modes."),
            ],
        ),
        Param(
            "theory.opt_level", "Geometry convergence", "choice", "Opt",
            help="The thresholds on energy change, gradient and displacement that end the optimization.",
            guidance="'Opt' matches Gaussian's default and is right for almost everything. TightOpt is "
                     "worth it when you are about to compute frequencies on a very flat system and want "
                     "to avoid residual imaginary modes.",
            choices=[
                Choice("Opt", "Opt (normal)", "Gradient < 3e-4 Eh/bohr. Standard.", (TAG_PAPER, TAG_RECOMMENDED)),
                Choice("LooseOpt", "LooseOpt", "Quick and dirty; do not follow with a frequency job.", (TAG_FAST,)),
                Choice("TightOpt", "TightOpt", "Gradient < 1e-4 Eh/bohr. Flat surfaces and pre-frequency polish."),
            ],
        ),
        Param(
            "theory.run_freq", "Run frequency calculation", "bool", True,
            help="Computes the Hessian after optimization. This is the ONLY way to prove your structure "
                 "is a genuine minimum rather than a saddle point.",
            guidance="Keep this on. The pipeline parses the frequencies and flags any structure with "
                     "imaginary modes instead of silently carrying it forward - which is exactly the check "
                     "the manuscript says it performed but never reports. Turn off only for a plumbing test.",
            cost="Analytic Hessian for a hybrid functional roughly doubles to triples the total job time.",
        ),
        Param(
            "theory.imaginary_cutoff", "Imaginary-mode tolerance", "float", -10.0, unit="cm^-1",
            help="Modes more negative than this are reported as genuine imaginary frequencies; small "
                 "negative values above it are treated as grid noise.",
            guidance="-10 cm^-1 is the usual convention. A mode at -5 cm^-1 is almost always numerical "
                     "noise in a floppy methoxy rotation; a mode at -150 cm^-1 means you are on a saddle "
                     "point and must re-optimize along that mode.",
            advanced=True,
        ),
        Param(
            "theory.saddle_max_attempts", "Automatic restarts off a saddle point", "int", 2,
            unit="restarts", minimum=0, maximum=5,
            help="How many times to re-optimize a structure that came back with imaginary "
                 "frequencies.",
            guidance="A frequency job reporting an imaginary mode has found a saddle point, "
                     "not a minimum - the gap, thermochemistry, Fukui indices and MEP all "
                     "then describe a structure that does not exist as a stable species. "
                     "With this above 0 the pipeline fixes it instead of just reporting it: "
                     "it reads the imaginary normal mode, pushes the geometry along it and "
                     "optimizes again, trying the opposite direction if the first falls back "
                     "into the same saddle. 2 clears most cases - an amino group that needs "
                     "to pyramidalise, a methyl stuck eclipsed, a ring at its flat inversion "
                     "point. Set 0 to get the old behaviour of a loud warning and nothing "
                     "else. Each restart is a full OPT+FREQ, so budget for it.",
            cost="One extra OPT+FREQ per restart, for affected molecules only.",
        ),
        Param(
            "theory.saddle_step_angstrom", "Restart displacement", "float", 0.30,
            unit="A", minimum=0.05, maximum=1.0, advanced=True,
            help="How far the atom that moves most is pushed along the imaginary mode "
                 "before re-optimizing.",
            guidance="0.3 A is enough to leave the saddle without landing somewhere "
                     "unrelated. Too small and the optimizer walks straight back; too large "
                     "and you can cross into a different conformer, which is a different "
                     "molecule for every purpose except the formula. Raise it only if a "
                     "restart keeps returning the same imaginary mode.",
        ),
        Param(
            "theory.saddle_small_mode_cm", "Small-mode threshold", "float", 50.0,
            unit="cm^-1", minimum=0.0, maximum=300.0, advanced=True,
            help="Below this magnitude an imaginary mode is treated as a numerical artifact "
                 "first and a real saddle second.",
            guidance="An imaginary frequency of a few tens of wavenumbers is usually a "
                     "methyl or hydroxyl rotation sitting in integration-grid noise, and the "
                     "first thing to try is a finer grid with a tighter optimization rather "
                     "than a displacement along a nearly flat coordinate. Above the "
                     "threshold the saddle is taken as real and the geometry is displaced "
                     "straight away.",
        ),
        Param(
            "theory.temperature", "Thermochemistry temperature", "float", 298.15, unit="K",
            help="Temperature for the ideal-gas/rigid-rotor/harmonic-oscillator free-energy corrections.",
            guidance="298.15 K for standard-state reporting. Use 310.15 K if you want the thermochemistry "
                     "to match the 310 K the manuscript used for its MD simulations.",
            advanced=True,
        ),
    ],
)

# ----------------------------------------------------------------------------- solvation
G_SOLVATION = Group(
    "solvation", "Implicit solvation", "The continuum-solvent treatment, applied to every stage.",
    [
        Param(
            "solvation.model", "Solvation model", "choice", "CPCM",
            help="The implicit-solvent model. The solute sits in a cavity inside a polarizable continuum "
                 "with the dielectric constant of the chosen solvent.",
            guidance="IMPORTANT: ORCA does not implement Gaussian's IEF-PCM. CPCM (conductor-like PCM) is "
                     "the closest equivalent and the two agree to well under 0.1 kcal/mol in high-dielectric "
                     "solvents such as water (eps 78) and DMSO (eps 47). They diverge more in low-dielectric "
                     "media - chloroform (eps 4.7) is where you would see the difference. SMD adds a "
                     "non-electrostatic cavitation/dispersion/repulsion term and is the better model for "
                     "solvation free energies, though it changes orbital energies slightly relative to the "
                     "published CPCM/IEF-PCM numbers.",
            choices=[
                Choice("CPCM", "CPCM",
                       "Conductor-like PCM. ORCA's closest analogue to the manuscript's IEF-PCM; "
                       "electrostatics only.", (TAG_PAPER, TAG_FAST)),
                Choice("SMD", "SMD",
                       "Truhlar's solvation model: CPCM electrostatics plus an empirical non-electrostatic "
                       "term. Best choice for solvation free energies and for logP-type reasoning.",
                       (TAG_RECOMMENDED,)),
                Choice("none", "Gas phase only",
                       "No continuum. Use when you want bare vertical ionization energies and electron "
                       "affinities, which is the only context where Koopmans-type descriptors are "
                       "physically interpretable."),
            ],
        ),
        Param(
            "solvation.solvents", "Solvent list", "multichoice", ["water"],
            help="One full set of calculations is run per entry. 'gas' adds a vacuum reference.",
            guidance="The manuscript reports vacuum plus water, methanol, chloroform and DMSO. Each "
                     "solvent multiplies the total cost by one, so start with water alone to validate "
                     "the workflow, then add the rest. Figures 2, 12, 13 and 14 are all built from the "
                     "solvent named in figures.reference_solvent.",
            cost="Linear: 4 solvents = 4x the DFT wall time.",
            choices=[
                Choice("water", "Water (eps 78.4)", "The manuscript's reference medium for Figure 2.",
                       (TAG_PAPER, TAG_RECOMMENDED, TAG_FAST)),
                Choice("methanol", "Methanol (eps 32.6)", "Polar protic.", (TAG_PAPER,)),
                Choice("chloroform", "Chloroform (eps 4.7)", "Low dielectric; where CPCM and IEF-PCM differ most.", (TAG_PAPER,)),
                Choice("dmso", "DMSO (eps 46.7)", "Polar aprotic.", (TAG_PAPER,)),
                Choice("gas", "Gas phase", "No continuum. The only rigorous reference for ionization energies.", (TAG_PAPER,)),
                Choice("acetonitrile", "Acetonitrile (eps 36.6)", "Polar aprotic alternative."),
                Choice("toluene", "Toluene (eps 2.4)", "Apolar aromatic."),
                Choice("thf", "THF (eps 7.3)", "Ethereal, low dielectric."),
            ],
        ),
    ],
)

# ------------------------------------------------------------------------------ property
G_PROPERTY = Group(
    "property", "Property single point",
    "A second, larger-basis calculation on the optimized geometry - the source of all descriptors.",
    [
        Param(
            "property.enabled", "Run property single point", "bool", True,
            help="Performs a single-point calculation at a larger basis on the optimized geometry, and "
                 "takes the orbital energies, populations and cubes from it.",
            guidance="Keep on. This is the step the manuscript implies but never describes: its Methods "
                     "give 6-31+G(d,p) while the MEP section quotes 6-311++G(d,p), with no single point "
                     "in between. Making it explicit is what lets you label every figure with the level "
                     "of theory that actually produced it.",
            cost="One SCF, no gradient - typically 5-15% of the optimization cost.",
        ),
        Param(
            "property.basis_sp", "Basis set (property step)", "choice", "def2-TZVP",
            help="Basis for the single point that generates orbital energies, populations, cubes and MEP.",
            guidance="This choice is a trade-off between two things you want from the same calculation, "
                     "and the right answer depends on which one the figure is for. Fully diffuse sets "
                     "(the '++' in 6-311++G(d,p), the 'aug-' in aug-cc-pVTZ) describe density far from "
                     "the nuclei, which is where the MEP tail and any genuine anion live. But a neutral "
                     "organic molecule has a NEGATIVE electron affinity: its LUMO is not a bound state at "
                     "all, and a diffuse set describes that unbound state as a discretised piece of the "
                     "continuum. The plotted 'LUMO' then comes out as a vague cloud sitting off the "
                     "molecule instead of the pi* orbital, and the HOMO-LUMO gap it reports is an artifact "
                     "of the basis rather than a property of the compound - benzene, phenol, aniline and "
                     "toluene all gave positive LUMO energies and 9-11 eV gaps at 6-311++G(d,p), against "
                     "bound pi* LUMOs and 5.5-6.7 eV at def2-TZVP. So: def2-TZVP for an orbital figure, "
                     "ma-def2-TZVP if you also want a well-converged MEP, and a fully diffuse set only "
                     "when you are deliberately studying an anion.",
            cost="Triple-zeta single point on an optimized geometry: minutes for a drug-sized "
                 "molecule. The diffuse sets are roughly 1.5-2x the cost of def2-TZVP and make the "
                 "SCF harder to converge.",
            choices=[
                Choice("def2-TZVP", "def2-TZVP",
                       "Karlsruhe triple-zeta, no diffuse functions. Bound pi* LUMO with the right "
                       "nodal pattern, stable orbital energies, well-behaved populations. The right "
                       "default for an orbital or gap figure.",
                       (TAG_RECOMMENDED,)),
                Choice("ma-def2-TZVP", "ma-def2-TZVP",
                       "Minimally augmented: one diffuse shell per atom. Recovers the MEP tail and "
                       "anion energetics without the continuum-LUMO problem. The compromise choice "
                       "when one single point has to serve both the MEP and the orbitals."),
                Choice("def2-TZVPD", "def2-TZVPD",
                       "Fully diffuse Karlsruhe set. For anions and electron affinities, where the "
                       "extra diffuse functions describe a state that really is bound."),
                Choice("6-311++G(d,p)", "6-311++G(d,p)",
                       "Pople triple-zeta, diffuse on heavy atoms and hydrogens. Common in the "
                       "literature. Expect an unbound LUMO for neutral molecules."),
                Choice("def2-SVP", "def2-SVP",
                       "Double-zeta. Too small for published orbital energies; useful for checking "
                       "that the pipeline runs.", (TAG_FAST,)),
            ],
        ),
        Param(
            "property.population", "Population analysis", "choice", "Hirshfeld",
            help="The scheme that partitions the total electron density into per-atom charges. These "
                 "charges are what the condensed Fukui functions are built from.",
            guidance="Hirshfeld is the manuscript's choice and the right one: it is basis-set stable, "
                     "always yields positive atomic populations, and is the standard for condensed Fukui "
                     "indices. Note that Gaussian 09 - the code the paper cites - has no native Hirshfeld "
                     "analysis, so those numbers must have come from an external tool that is never named. "
                     "ORCA computes it directly.",
            choices=[
                Choice("Hirshfeld", "Hirshfeld",
                       "Stockholder partitioning against promolecular atomic densities. Basis-set robust; "
                       "the standard for condensed Fukui functions.", (TAG_PAPER, TAG_RECOMMENDED)),
                Choice("Loewdin", "Loewdin",
                       "Symmetrically orthogonalized Mulliken. More stable than Mulliken, still "
                       "basis-dependent."),
                Choice("Mulliken", "Mulliken",
                       "The classic. Badly behaved with diffuse functions - can give negative populations. "
                       "Avoid with any '+' or 'aug-' basis.", (TAG_FAST,)),
                Choice("Mayer", "Mayer",
                       "Gives bond orders and valences in addition to charges."),
            ],
        ),
        Param(
            "property.koopmans_mode", "Descriptor derivation", "choice", "both",
            help="How ionization energy and electron affinity are obtained for the global reactivity "
                 "descriptors.",
            guidance="'koopmans' takes IE = -e(HOMO) and EA = -e(LUMO), which is what the manuscript did. "
                     "Be aware this is a Hartree-Fock theorem: in Kohn-Sham DFT the HOMO equals the exact "
                     "IE only for the exact functional, and B3LYP typically underestimates IE by 1-3 eV "
                     "(the paper's 5.5 eV for a carbazole should be compared with ~7.6 eV experimental). "
                     "'delta_scf' instead computes real energy differences from the N+1 and N-1 "
                     "calculations - physically meaningful, and free if you are already running Fukui. "
                     "'both' reports the two side by side, which is the honest option.",
            cost="delta_scf needs two extra single points per molecule - already paid for if fukui.enabled.",
            choices=[
                Choice("both", "Both (recommended)",
                       "Report Koopmans and delta-SCF values together so the approximation is visible.",
                       (TAG_RECOMMENDED,)),
                Choice("koopmans", "Koopmans only",
                       "IE = -e(HOMO), EA = -e(LUMO). Reproduces the manuscript's Table S1 approach.", (TAG_PAPER, TAG_FAST)),
                Choice("delta_scf", "Delta-SCF only",
                       "IE = E(cation) - E(neutral), EA = E(neutral) - E(anion), all at the neutral "
                       "geometry. Vertical, and physically defensible."),
            ],
        ),
    ],
)

# --------------------------------------------------------------------------------- fukui
G_FUKUI = Group(
    "fukui", "Fukui functions", "Local reactivity: which atoms are electrophilic or nucleophilic.",
    [
        Param(
            "fukui.enabled", "Compute condensed Fukui functions", "bool", True,
            help="Runs the N+1 and N-1 single points needed for f+, f-, f0, the dual descriptor and the "
                 "multiphilic descriptor.",
            guidance="Required for the Fukui panels and the supplementary tables. Also supplies the "
                     "delta-SCF ionization energies and electron affinities at no extra cost.",
            cost="Two additional single points per molecule per solvent - roughly 2x the property step.",
        ),
        Param(
            "fukui.geometry", "Geometry for the ion calculations", "choice", "frozen",
            help="Whether the cation and anion are computed at the neutral geometry or re-optimized.",
            guidance="'frozen' gives VERTICAL quantities, which is what Fukui functions are defined from "
                     "and what you want for reactivity indices. 'relaxed' gives adiabatic values, which "
                     "are the right thing for redox potentials but NOT for Fukui analysis. Keep frozen "
                     "unless you specifically want adiabatic ionization energies.",
            choices=[
                Choice("frozen", "Frozen (vertical)",
                       "Single points at the optimized neutral geometry. Correct for Fukui functions and "
                       "vertical IE/EA.", (TAG_PAPER, TAG_RECOMMENDED, TAG_FAST)),
                Choice("relaxed", "Relaxed (adiabatic)",
                       "Re-optimize each ion. Gives adiabatic IE/EA; much more expensive and not what "
                       "Fukui indices are defined from."),
            ],
        ),
        Param(
            "fukui.smear_scheme", "Degenerate-orbital handling", "choice", "none",
            help="Whether to average Fukui indices over near-degenerate frontier orbitals.",
            guidance="Leave at 'none' for closed-shell organics. Symmetric molecules with degenerate "
                     "frontier orbitals can give misleading single-orbital pictures; in that case the "
                     "finite-difference route used here is already the robust answer.",
            advanced=True,
            choices=[
                Choice("none", "None", "Straight finite differences of the atomic populations.", (TAG_PAPER, TAG_RECOMMENDED)),
            ],
        ),
    ],
)

# --------------------------------------------------------------------------------- cubes
G_CUBES = Group(
    "cubes", "Volumetric data (cube files)",
    "The 3D grids behind the orbital and MEP panels of Figure 2.",
    [
        Param(
            "cubes.enabled", "Generate cube files", "bool", True,
            help="Writes Gaussian-format .cube files for the frontier orbitals, the electron density and "
                 "the electrostatic potential.",
            guidance="Needed for Figure 2 panels B and C and for Figure 13. The cubes are also directly "
                     "loadable in Avogadro and GaussView if you want to look around interactively.",
            cost="Cube generation is cheap; the MEP evaluation is the slow part (see cubes.esp).",
        ),
        Param(
            "cubes.grid_points", "Grid points per axis", "int", 80, unit="points",
            help="Resolution of the volumetric grid. The cube contains this value cubed points in total.",
            guidance="80 gives smooth publication-quality isosurfaces for a 30-atom molecule. 40 is "
                     "visibly faceted and should only be used for testing. Going past 120 produces files "
                     "of hundreds of MB for no visible gain. Remember the cost is cubic: 80 -> 512,000 "
                     "points, 120 -> 1,728,000.",
            cost="Cubic in this number. Dominates MEP evaluation time and disk usage.",
            minimum=20, maximum=200,
        ),
        Param(
            "cubes.orbitals", "Which orbitals to plot", "choice", "homo_lumo",
            help="Selects the molecular orbitals written as cubes.",
            guidance="homo_lumo covers Figure 2B. Extend to homo_lumo_pm1 when a frontier orbital is "
                     "nearly degenerate with the one below/above it and a single-orbital picture would "
                     "mislead - check the orbital energies in descriptors.csv first.",
            choices=[
                Choice("homo_lumo", "HOMO and LUMO", "The two frontier orbitals. What Figure 2B shows.",
                       (TAG_PAPER, TAG_RECOMMENDED, TAG_FAST)),
                Choice("homo_lumo_pm1", "HOMO-1 ... LUMO+1", "Four orbitals - use when frontier orbitals are near-degenerate."),
                Choice("homo_lumo_pm2", "HOMO-2 ... LUMO+2", "Six orbitals. Rarely needed for a figure."),
            ],
        ),
        Param(
            "cubes.esp", "Compute electrostatic potential", "bool", True,
            help="Evaluates V(r) at every grid point with orca_vpot and writes a matching ESP cube, which "
                 "is then painted onto the density isosurface to make the MEP map.",
            guidance="Required for Figure 2C. Note this is a genuine quantum-mechanical MEP from the SCF "
                     "density, not a point-charge approximation.",
            cost="Scales as (basis functions)^2 x (grid points). Roughly 5-10 min for a 30-atom molecule "
                 "on an 80^3 grid at triple-zeta; use cubes.esp_nprocs to parallelize.",
        ),
        Param(
            "cubes.esp_nprocs", "Parallel processes for MEP", "int", 1, unit="processes",
            minimum=1, maximum=64,
            help="Number of MPI processes for orca_vpot. 1 runs the serial binary, which "
                 "is the only form that works reliably after a completed job.",
            guidance="Leave this at 1. The parallel build of orca_vpot re-reads the SHARK "
                     "basis scratch file <job>.SHARKINP.tmp, and ORCA deletes that file "
                     "whenever a calculation terminates normally - so the parallel path "
                     "fails for exactly the jobs you want a MEP from, with a 'Failed to "
                     "read NAtoms, NBasis, NShells' error. The serial binary rebuilds "
                     "everything from the .gbw and is fast enough: 512,000 grid points on "
                     "a 12-atom molecule takes about 37 s. If you do raise it, a failure "
                     "falls back to serial automatically.",
            advanced=True,
        ),
        Param(
            "cubes.keep_cubes", "Keep cube files after rendering", "bool", True,
            help="Whether to retain the .cube files once the figures have been produced.",
            guidance="Keep them if you want to re-render with different isovalues or inspect the orbitals "
                     "in Avogadro later - regenerating means rerunning the SCF. Set to false only if disk "
                     "space is tight: an 80^3 cube is about 14 MB and you get four or more per molecule "
                     "per solvent.",
            advanced=True,
        ),
    ],
)

# -------------------------------------------------------------------------------- render
G_RENDER = Group(
    "render", "3D rendering", "How the orbital, MEP and lipophilicity surfaces are drawn.",
    [
        Param(
            "render.engine", "Rendering engine", "choice", "pymol",
            help="The program that turns cube files into images.",
            guidance="PyMOL is scripted headlessly, needs no clicking, and produces consistent camera "
                     "angles across panels - which is what makes a multi-molecule figure readable. VMD is "
                     "an alternative if you prefer its isosurface shading. 'none' skips rendering and "
                     "leaves you the cubes to open in Avogadro or GaussView yourself.",
            choices=[
                Choice("pymol", "PyMOL",
                       "Scripted via generated .pml files and ray-traced offscreen. Best control over "
                       "colour ramps and camera.", (TAG_PAPER, TAG_RECOMMENDED, TAG_FAST)),
                Choice("vmd", "VMD",
                       "Scripted via generated .tcl files with Tachyon rendering. Nicer default "
                       "isosurface shading; fiddlier colour ramps."),
                Choice("none", "None - cubes only",
                       "Stop after cube generation. Open the .cube files in Avogadro or GaussView and "
                       "render by hand."),
            ],
        ),
        Param(
            "render.width", "Image width", "int", 1600, unit="px",
            help="Pixel width of each rendered 3D panel before assembly into the figure.",
            guidance="1600 px per panel gives a 4-panel row at roughly 300 dpi across a journal page. "
                     "Raise to 2400 for a full-page figure; ray-trace time scales with the pixel count.",
            minimum=200,
        ),
        Param(
            "render.height", "Image height", "int", 1200, unit="px",
            help="Pixel height of each rendered 3D panel.",
            guidance="Keep the aspect ratio close to the shape of your molecules. Carbazoles are wide and "
                     "flat, so 4:3 or 16:9 wastes less space than a square.",
            minimum=200,
        ),
        Param(
            "render.ray_trace", "Ray-trace images", "bool", True,
            help="Use PyMOL's ray tracer rather than the OpenGL preview for the final image.",
            guidance="On for anything going into a manuscript - it gives proper shadows, smooth surfaces "
                     "and anti-aliased edges. Off while you are iterating on camera angles, where it just "
                     "costs you time.",
            cost="Seconds to a minute per panel depending on pixel count and surface complexity.",
        ),
        Param(
            "render.orbital_isovalue", "Orbital isosurface value", "float", 0.02, unit="e^(1/2) bohr^(-3/2)",
            help="The contour level at which molecular-orbital lobes are drawn. Both +value and -value "
                 "surfaces are rendered, in two different colours.",
            guidance="0.02 is the near-universal convention in the literature and what most published "
                     "HOMO/LUMO figures use; 0.05 gives tighter, smaller lobes and 0.01 gives diffuse, "
                     "sprawling ones that can obscure the molecular skeleton. Whatever you pick, use the "
                     "SAME value for every molecule in the figure or the comparison is meaningless - and "
                     "state it in the caption, which the manuscript does not.",
            minimum=0.001,
        ),
        Param(
            "render.density_isovalue", "Density isosurface value", "float", 0.0004, unit="e/bohr^3",
            help="The electron-density contour that defines the molecular surface onto which the "
                 "electrostatic potential is painted for the MEP map.",
            guidance="0.0004 a.u. is Bader's definition of the molecular surface and encloses roughly 99% "
                     "of the electron density - the standard choice for MEP maps and the one that "
                     "corresponds to the van der Waals surface. 0.002 a.u. gives a tighter surface closer "
                     "to the nuclei, with systematically larger |V| values. Must be identical across "
                     "panels, and must be stated in the caption.",
            minimum=0.00001,
        ),
        Param(
            "render.esp_range", "MEP colour-scale limits", "float", 0.05, unit="a.u.",
            help="The electrostatic potential is coloured from -value (red, electron-rich) to +value "
                 "(blue, electron-poor). 0.05 a.u. is about 31 kcal/mol.",
            guidance="This is the single most abused setting in published MEP figures. If each panel is "
                     "auto-scaled to its own min/max - which is what most GUIs do by default - then "
                     "panels CANNOT be compared, because the same colour means a different potential in "
                     "each one. Fixing one range across all molecules is what makes Figure 2C a real "
                     "comparison. 0.05 a.u. suits neutral organics with hydroxyl and methoxy groups; "
                     "widen it if your series includes charged species.",
            minimum=0.001,
        ),
        Param(
            "render.esp_scale", "MEP scaling mode", "choice", "auto_shared",
            help="How the electrostatic-potential colour range is chosen across the panels.",
            guidance="A fixed range is right for a set of neutral molecules and is what makes "
                     "panels comparable. It breaks for ions: the potential around a +1 cation is "
                     "positive essentially everywhere and far outside +-0.05 a.u., so the panel "
                     "saturates into a featureless blue blob and all the chemical detail is lost. "
                     "'auto_shared' measures the potential on each molecular surface and picks ONE "
                     "range that covers them all - still a fair comparison, and it prints the range "
                     "on the colour bar. 'per_panel' scales each molecule separately: use it only "
                     "to show one molecule's own detail, never to compare, and the figure then "
                     "labels each panel with its own range. 'fixed' uses render.esp_range exactly.",
            choices=[
                Choice("auto_shared", "One range covering every panel, measured from the data",
                       "Comparable across panels and safe for mixed charges.",
                       (TAG_RECOMMENDED, TAG_FAST)),
                Choice("fixed", "Exactly render.esp_range",
                       "Reproducible and comparable, but saturates for charged species.",
                       (TAG_PAPER,)),
                Choice("per_panel", "Each molecule scaled to its own range",
                       "Maximum detail per molecule; panels are NOT comparable."),
            ],
        ),
        Param(
            "render.esp_autoscale", "Auto-scale MEP per molecule (legacy)", "bool", False,
            help="If true, each MEP panel is scaled to its own potential range instead of the fixed "
                 "esp_range.",
            guidance="Leave this OFF for any comparative figure. It exists only because reviewers "
                     "sometimes ask to see a single molecule's full dynamic range. If you turn it on, the "
                     "pipeline writes the per-panel range into the figure so the reader is not misled.",
        ),
        Param(
            "render.esp_palette", "MEP colour palette", "choice", "paper",
            help="The colour ramp used for the electrostatic potential surface.",
            guidance="'paper' reproduces the manuscript's blue-green-yellow-orange-red sequence. "
                     "'rwb' (red-white-blue) is the cleaner scientific choice because white marks the "
                     "zero of the potential unambiguously, and it survives greyscale printing. Rainbow "
                     "ramps are not perceptually uniform and exaggerate mid-range differences.",
            choices=[
                Choice("paper", "Blue-green-yellow-orange-red",
                       "The manuscript's ramp. Red = negative potential (electron-rich), blue = positive.",
                       (TAG_PAPER,)),
                Choice("rwb", "Red-white-blue",
                       "Diverging, perceptually sane, zero is unmistakable. The better default.", (TAG_RECOMMENDED,)),
                Choice("rainbow", "Rainbow (PyMOL default)", "Included for continuity with older figures. Not recommended."),
            ],
        ),
        Param(
            "render.orbital_colors", "Orbital lobe colours", "str", "firebrick,skyblue",
            help="Two comma-separated colour names for the positive and negative phase of each orbital.",
            guidance="Any PyMOL colour name works. Keep the two phases strongly contrasting; avoid "
                     "red/green for colour-blind readers. The classic choices are red/blue or orange/cyan.",
        ),
        Param(
            "render.background", "Background", "choice", "white",
            help="Background of the rendered panels.",
            guidance="'transparent' is the most flexible for figure assembly - it composites onto any "
                     "page colour. 'white' is safest if the journal flattens alpha channels.",
            choices=[
                Choice("white", "White", "Opaque white. Safe everywhere.", (TAG_PAPER, TAG_FAST)),
                Choice("transparent", "Transparent", "Alpha channel preserved; best for assembly.", (TAG_RECOMMENDED,)),
                Choice("black", "Black", "For slides and posters."),
            ],
        ),
        Param(
            "render.orientation", "Molecular orientation", "choice", "principal_plane",
            help="How each molecule is rotated before rendering.",
            guidance="'principal_plane' puts the molecular plane in the screen plane using the inertia "
                     "tensor, which is what you want for flat fused-ring systems - every panel ends up "
                     "oriented the same way, so the reader compares chemistry rather than viewpoint. "
                     "'first_molecule' aligns everything onto the first molecule's frame, which is better "
                     "for a congeneric series that shares a scaffold. 'none' keeps the ORCA orientation, "
                     "which is essentially arbitrary.",
            choices=[
                Choice("principal_plane", "Principal plane (inertia tensor)",
                       "Molecular plane parallel to the screen. Best for planar aromatics.",
                       (TAG_PAPER, TAG_RECOMMENDED, TAG_FAST)),
                Choice("first_molecule", "Align to first molecule",
                       "Least-squares fit of every molecule onto the first by maximum common substructure. "
                       "Best for a congeneric series."),
                Choice("none", "As calculated", "Whatever orientation ORCA left. Not recommended for figures."),
            ],
        ),
        Param(
            "render.stick_radius", "Stick radius", "float", 0.12, unit="A",
            help="Bond cylinder radius in the ball-and-stick representation.",
            guidance="0.12-0.15 A reads well at figure scale. Thinner sticks disappear when the figure is "
                     "reduced to a single column.",
            advanced=True,
        ),
        Param(
            "render.sphere_scale", "Atom sphere scale", "float", 0.22,
            help="Atomic sphere radius as a fraction of the van der Waals radius.",
            guidance="0.22 gives the conventional ball-and-stick look. 1.0 would give a space-filling "
                     "model, which hides the orbital lobes.",
            advanced=True,
        ),
        Param(
            "render.surface_quality", "Surface mesh quality", "int", 2,
            help="Tessellation density of the molecular and isodensity surfaces. Higher values "
                 "subdivide the mesh further, so curved regions read as smooth rather than faceted.",
            guidance="2 is the right setting for a figure; 0-1 leave visible facets on small "
                     "molecules, which is what makes a surface look blocky. 3 and above rarely "
                     "changes anything you can see and makes ray tracing noticeably slower. If a "
                     "surface still looks angular at 2, the problem is usually colour banding "
                     "rather than geometry - raise render.surface_color_smoothing instead.",
            cost="Each step roughly doubles the vertex count and the ray-trace time for that surface.",
            minimum=0, maximum=4,
        ),
        Param(
            "render.surface_color_smoothing", "Surface colour smoothing", "float", 0.80,
            help="How far per-atom colours are blended across the molecular surface "
                 "(PyMOL's surface_color_smoothing_threshold).",
            guidance="This controls colour banding, not geometry. Each surface vertex takes its "
                     "colour from the nearest atom, so with little smoothing you get hard-edged "
                     "patches wherever two atoms of different value meet - the blocky look on a "
                     "lipophilicity map. 0.30 blends them into a continuous gradient, which is the "
                     "honest rendering of a fragment-based potential anyway, since the underlying "
                     "quantity is a distance-weighted sum rather than a per-atom step function. "
                     "Compared side by side on a carbazole: PyMOL's default of 0.05 and a value of "
                     "0.30 both leave visible angular wedges; 0.80 gives a continuous gradient "
                     "while each feature is still clearly attributable to its atom; 1.00 begins to "
                     "wash out small localized features. Lower it toward 0.05 if you want to see "
                     "exactly which atom owns which patch. Affects only per-atom colouring "
                     "(Figure 13), not the ramp-coloured MEP surface, which is already continuous.",
            minimum=0.0, maximum=1.0,
        ),
        Param(
            "render.surface_transparency", "Surface transparency", "float", 0.25,
            help="Transparency of the MEP and lipophilicity surfaces, 0 = opaque, 1 = invisible.",
            guidance="0.2-0.3 lets the molecular skeleton show through the surface, which is what makes a "
                     "MEP map interpretable - otherwise the reader cannot tell which atom a red patch "
                     "belongs to. Set to 0 if you want the classic opaque surface.",
            minimum=0.0, maximum=1.0,
        ),
        Param(
            "render.show_hydrogens", "Show hydrogens", "bool", True,
            help="Whether hydrogen atoms appear in the 3D panels.",
            guidance="Keep them: the MEP panels of this manuscript specifically discuss the positive "
                     "potential over particular hydrogens, and hiding them would make that unreadable.",
        ),
        Param(
            "render.pymol_path", "PyMOL executable", "path", "",
            help="Full path to the PyMOL binary. Empty means search PATH, then the directory of "
                 "the Python interpreter running the pipeline.",
            guidance="Only set this if you have several PyMOL installations and need a specific one - "
                     "for example a licensed Incentive PyMOL rather than the open-source build.",
            advanced=True,
        ),
        Param(
            "render.vmd_path", "VMD executable", "path", "",
            help="Full path to the VMD binary. Only used when render.engine = 'vmd'.",
            guidance="Empty means search PATH.",
            advanced=True,
        ),
        Param(
            "render.label_atoms", "Atom-number labels", "str", "",
            help="Comma-separated 1-based atom indices to label in the geometry panel, e.g. '3,4,9'.",
            guidance="Use this to fix the manuscript's worst reproducibility problem - it quotes "
                     "dihedrals over atoms C10, C12 and N11, which do not exist in carbazole numbering, "
                     "and uses three mutually incompatible schemes across sections. Labelling the atoms "
                     "you actually discuss removes all ambiguity. Leave empty for no labels.",
        ),
    ],
)

# ------------------------------------------------------------------------------- figures
G_FIGURES = Group(
    "figures", "Figure assembly", "Which figures to build and how they are typeset.",
    [
        Param(
            "figures.which", "Figures to produce", "multichoice",
            ["frontier_orbitals", "speciation", "lipophilicity", "mpo_metabolism"],
            help="Selects which figures are assembled. The four below are built in; any "
                 "plugin discovered in figures.plugin_dir or registered by an installed "
                 "package can be named here too.",
            guidance="frontier_orbitals needs the full DFT chain: cubes and 3D rendering. The "
                     "other three need only the optimized geometries plus the derived "
                     "properties, so they are quick to iterate on once the DFT has run once. "
                     "The older names fig2 / fig12 / fig13 / fig14 still work and are migrated "
                     "automatically.",
            choices=[
                Choice("frontier_orbitals", "Geometry, frontier orbitals and MEP",
                       "Three rows: optimized structures with the minimum check, HOMO and LUMO "
                       "with the gap, and the electrostatic potential mapped on the electron "
                       "density. Requires DFT, cubes and rendering.",
                       (TAG_RECOMMENDED, TAG_FAST)),
                Choice("speciation", "Acid-base speciation vs pH",
                       "Annotated 2D structures showing every ionizable centre, plus net "
                       "molecular charge across pH 0-14. pKa is predicted from structure unless "
                       "you supply measured values.", (TAG_RECOMMENDED,)),
                Choice("lipophilicity", "Lipophilicity surfaces",
                       "Per-atom Crippen lipophilicity painted on the molecular surface, with "
                       "TPSA and clogP. Requires optimized geometries and rendering.",
                       (TAG_RECOMMENDED,)),
                Choice("mpo_metabolism", "CNS MPO space and metabolism maps",
                       "logD against TPSA with the CNS MPO score, plus site-of-metabolism maps "
                       "weighted by the radical Fukui index.", (TAG_RECOMMENDED,)),
            ],
        ),
        Param(
            "figures.plugin_dir", "Where your own figures live", "path", "dft_figures",
            help="Directory scanned for plain .py figure modules, relative to the project "
                 "directory. Anything found there appears in figures.which alongside the "
                 "built-ins.",
            guidance="This is the no-packaging route: drop a file defining NAME, REQUIRES "
                     "and build(ctx) into ./dft_figures and it is discovered on the next "
                     "run. A figure you want to share instead goes in a pip-installable "
                     "package that registers the 'dftfig.figures' entry-point group - both "
                     "routes land in the same registry. See docs/writing-a-figure.md.",
            advanced=True,
        ),
        Param(
            "figures.reference_solvent", "Solvent shown in the figures", "choice", "water",
            help="Which solvent's results are used for the figure panels when several were calculated.",
            guidance="The manuscript shows water in Figure 2 and relegates the other solvents to the "
                     "supplement. Must be one of the solvents you actually calculated.",
            choices=[
                Choice("water", "Water", "The manuscript's choice for Figure 2.", (TAG_PAPER, TAG_RECOMMENDED, TAG_FAST)),
                Choice("gas", "Gas phase", "The only medium in which Koopmans descriptors are physically interpretable."),
                Choice("methanol", "Methanol", ""),
                Choice("chloroform", "Chloroform", ""),
                Choice("dmso", "DMSO", ""),
            ],
        ),
        Param(
            "figures.dpi", "Output resolution", "int", 300, unit="dpi",
            help="Raster resolution of the assembled figures.",
            guidance="300 dpi is the floor for most journals; 600 for line art if the publisher asks. "
                     "Vector formats ignore this for the plot elements but still use it for the embedded "
                     "3D renders.",
            minimum=72,
        ),
        Param(
            "figures.formats", "Output formats", "multichoice", ["png", "pdf"],
            help="File formats written for each figure.",
            guidance="PNG for looking at and for submission systems that demand raster; PDF for vector "
                     "text and lines that stay sharp at any zoom. TIFF only if a journal insists.",
            choices=[
                Choice("png", "PNG", "Raster, universally accepted.", (TAG_PAPER, TAG_RECOMMENDED, TAG_FAST)),
                Choice("pdf", "PDF", "Vector plot elements, embedded raster renders.", (TAG_RECOMMENDED,)),
                Choice("svg", "SVG", "Editable vector - open in Inkscape to tweak labels by hand."),
                Choice("tiff", "TIFF", "Large files; some journals require it."),
            ],
        ),
        Param(
            "figures.panel_label_style", "Panel letter style", "choice", "upper",
            help="How panel letters are drawn in the assembled figures.",
            guidance="Match your target journal. Molecular Diversity uses bare capitals (A, B, C), which "
                     "is the default here.",
            choices=[
                Choice("upper", "A B C", "Bare capitals, as in the manuscript.", (TAG_PAPER, TAG_RECOMMENDED, TAG_FAST)),
                Choice("lower", "a b c", "Lower case."),
                Choice("paren_upper", "(A) (B) (C)", "Parenthesized capitals."),
                Choice("none", "No letters", "For a single-panel figure."),
            ],
        ),
        Param(
            "figures.ph_range", "pH range for Figure 12B", "str", "0,14",
            help="Minimum and maximum pH of the speciation curve, comma separated.",
            guidance="0-14 shows the full picture. Narrow to 4,10 if you only care about the "
                     "physiologically relevant window and want the transitions to fill the axes.",
        ),
        Param(
            "figures.physiological_ph", "Physiological pH marker", "float", 7.4,
            help="The pH at which the predominant species is annotated in Figure 12B.",
            guidance="7.4 is blood plasma. Use 6.5 for the small-intestine lumen or 1.5 for gastric fluid "
                     "if you are arguing about oral absorption.",
        ),
        Param(
            "figures.mlp_palette", "Lipophilicity palette", "choice", "paper",
            help="Colour ramp for the molecular lipophilicity potential surface in Figure 13.",
            guidance="'paper' matches the manuscript's caption exactly - red for the least lipophilic "
                     "fragments through to blue for the most lipophilic. Note this is the OPPOSITE of the "
                     "convention most medicinal chemists expect (blue usually means polar), so whichever "
                     "you choose, the caption has to say so. 'conventional' uses blue = hydrophilic, "
                     "orange = lipophilic.",
            choices=[
                Choice("paper", "Red (polar) to blue (lipophilic)",
                       "The manuscript's stated ramp.", (TAG_PAPER,)),
                Choice("conventional", "Blue (polar) to orange (lipophilic)",
                       "The more widely used convention in medicinal chemistry.", (TAG_RECOMMENDED,)),
            ],
        ),
        Param(
            "figures.mlp_transparency", "Lipophilicity surface transparency", "float", 0.35,
            help="Transparency of the Figure 13 surfaces: 0 = opaque, 1 = invisible.",
            guidance="A fully opaque surface hides the molecule, so the reader cannot tell which "
                     "atom a coloured patch belongs to - the same objection that applies to an "
                     "unlabelled MEP map. 0.3-0.4 lets the skeleton read through while keeping the "
                     "colours saturated. Push past 0.5 and the colour washes out; set 0 only if you "
                     "specifically want the classic solid surface.",
            minimum=0.0, maximum=1.0,
        ),
        Param(
            "figures.mlp_show_structure", "Show the molecule inside the surface", "bool", True,
            help="Draws the optimized geometry as thin sticks underneath the translucent "
                 "lipophilicity surface.",
            guidance="Keep this on. It is what makes the panel interpretable: a blue patch means "
                     "nothing until you can see it sits over the chlorine rather than the methoxy. "
                     "Turn it off only for a purely qualitative shape comparison.",
        ),
        Param(
            "figures.mpo_variant", "CNS MPO implementation", "choice", "wager2010",
            help="Which published set of desirability functions is used to score CNS drug-likeness.",
            guidance="wager2010 is the original six-parameter CNS MPO (clogP, clogD, MW, TPSA, HBD, "
                     "most-basic pKa), each mapped to a 0-1 desirability and summed to a 0-6 score. Note "
                     "that the manuscript's Methods state the HBD threshold as <=1, whereas Wager's "
                     "published desirability function for HBD is centred on <=0.5 - the pipeline uses the "
                     "published function and will therefore differ slightly from the paper's Table 3.",
            choices=[
                Choice("wager2010", "Wager et al. 2010 (6 parameters)",
                       "The original and still standard CNS MPO.", (TAG_PAPER, TAG_RECOMMENDED, TAG_FAST)),
            ],
        ),
        Param(
            "figures.mpo_binary_verdict", "Show accept/reject verdict", "bool", False,
            help="Whether to print a pass/fail label alongside the CNS MPO score.",
            guidance="Leave OFF. CNS MPO is a continuous desirability score and its authors explicitly "
                     "propose MPO >= 4 as an enrichment criterion, not a filter. The manuscript's Table 3 "
                     "binarizes it into 'Accepted'/'Rejected', which puts a score of 3.76 on the wrong "
                     "side of a boundary the method does not have. Turn it on only if a co-author insists.",
        ),
        Param(
            "figures.properties_csv", "Physicochemical properties file", "path",
            "data/properties.csv",
            help="Optional CSV of measured or externally predicted values. Anything in it "
                 "overrides what the pipeline works out for itself.",
            guidance="You do not need this file. Give the pipeline SDF structures and it "
                     "computes MW, TPSA, HBD and Crippen logP from the connection table, "
                     "predicts pKa from the SMARTS ionization model, and derives logD at "
                     "pH 7.4 from the two. Supply the CSV when you have measured pKa or "
                     "logD, or when you want permeability, clearance, plasma protein "
                     "binding or P-gp on the figures - those need a trained ADMET model "
                     "(ADMETlab 3.0, Deep-PK, your own) and the pipeline will not invent "
                     "them. Columns are matched by name; blanks are filled, not overwritten.",
        ),
        Param(
            "figures.pka_source", "Where pKa comes from", "choice", "predict",
            help="How the acidic and basic pKa values behind Figure 12 and logD are obtained.",
            guidance="'predict' needs nothing from you: a curated SMARTS table assigns a class "
                     "average pKa to every ionizable group it recognises - carboxylic acid, "
                     "phenol, aliphatic amine, pyridine, imidazole, tetrazole and about thirty "
                     "more - and reports the groups it found alongside the number. Expect +-1 "
                     "to 1.5 log units, which settles whether a compound is charged at pH 7.4 "
                     "but is not a value to quote in a paper. 'csv' uses only what you supply "
                     "and leaves the rest blank. 'none' switches ionization off entirely, which "
                     "draws every compound at the formal charge of the structure you gave.",
            choices=[
                Choice("predict", "Predict from structure",
                       "SMARTS class averages; fills any pKa the CSV does not give.",
                       (TAG_RECOMMENDED, TAG_FAST)),
                Choice("csv", "Only what the properties CSV gives",
                       "Nothing is inferred. Panels without a value say so.", (TAG_PAPER,)),
                Choice("none", "No ionization model",
                       "Every compound is drawn at its formal charge."),
            ],
        ),
        Param(
            "figures.logd_source", "Where logD(7.4) comes from", "choice", "derive",
            help="How the buffer lipophilicity on Figure 14 is obtained.",
            guidance="'derive' computes logD from Crippen logP and the ionization model: only "
                     "the neutral microspecies is taken to partition, which is the standard "
                     "approximation. It is a lower bound for strong bases, where ion pairing "
                     "adds roughly one unit, and it is a floor rather than a prediction for a "
                     "permanently charged compound, which has no neutral form at all. 'csv' "
                     "uses only supplied values. 'none' leaves the axis empty unless the CSV "
                     "fills it.",
            choices=[
                Choice("derive", "Derive from clogP and the predicted pKa",
                       "Self-contained; labelled on the figure as derived.",
                       (TAG_RECOMMENDED, TAG_FAST)),
                Choice("csv", "Only what the properties CSV gives",
                       "For measured logD or an external ADMET prediction.", (TAG_PAPER,)),
                Choice("none", "No logD", "Figure 14 panel A falls back to clogP."),
            ],
        ),
        Param(
            "figures.som_source", "Site-of-metabolism source", "choice", "fukui_radical",
            help="Where the per-atom metabolic-lability weights for Figure 14 B-E come from.",
            guidance="'fukui_radical' uses the condensed radical Fukui index f0 that this pipeline "
                     "already computes - susceptibility to radical attack is a recognized "
                     "quantum-chemical proxy for the hydrogen abstraction and aromatic oxidation "
                     "steps that CYP450 performs, so the panel is fully reproducible from your own "
                     "DFT with no external service. It is a reactivity proxy, NOT a trained "
                     "site-of-metabolism model: it knows nothing about how the substrate binds in the "
                     "enzyme active site, which is often what decides the real site. 'csv' reads "
                     "measured or predicted probabilities you export from XenoSite or ADMETlab; use "
                     "that when you have them and say in the caption which source you used. 'both' "
                     "prefers the CSV and silently falls back to Fukui where a file is missing.",
            choices=[
                Choice("fukui_radical", "Radical Fukui index f0 (from this pipeline's DFT)",
                       "Self-contained and reproducible; a reactivity proxy, blind to binding mode.",
                       (TAG_RECOMMENDED, TAG_FAST)),
                Choice("csv", "CSV export from XenoSite / ADMETlab",
                       "A trained site-of-metabolism model. What the manuscript used.", (TAG_PAPER,)),
                Choice("both", "CSV where available, Fukui otherwise",
                       "Prefer the trained model, fall back to the DFT proxy."),
            ],
        ),
        Param(
            "figures.som_dir", "Site-of-metabolism directory", "path", "data/som",
            help="Directory holding <molecule>.csv files mapping atom index to metabolism probability, "
                 "used for Figure 14 panels B-E.",
            guidance="Export these from XenoSite or ADMETlab's SOM module. Format: two columns, "
                     "atom_index (1-based, matching the SMILES atom order) and probability (0-1). If a "
                     "file is missing, that panel is drawn greyed out with a note rather than silently "
                     "omitted.",
        ),
    ],
)

# ----------------------------------------------------------------------------- execution
G_EXECUTION = Group(
    "execution", "Execution", "Where ORCA lives, how much hardware to use, and restart behaviour.",
    [
        Param(
            "execution.orca_path", "ORCA executable", "path", "",
            help="Full path to the ORCA binary. Leave empty to search this machine automatically. "
                 "The helper tools orca_plot and orca_vpot are located next to it.",
            guidance="Empty means auto-discovery, which is almost always what you want: dftfig "
                     "checks $ORCA_PATH, then your PATH, then the usual install locations "
                     "(~/orca*, /opt/orca*, /usr/local/orca*, /share/apps/orca*), and remembers "
                     "what it finds in ~/.config/dftfig/paths.yaml so it only has to look once. "
                     "A candidate is accepted only if orca_plot and orca_vpot sit beside it - on "
                     "most Linux desktops /usr/bin/orca is the GNOME Orca SCREEN READER, and "
                     "running that instead produces a failure that makes no sense. Set this "
                     "explicitly if you have several ORCA versions and want a particular one. It "
                     "must be the absolute path, never bare 'orca': ORCA needs it for its own "
                     "parallel startup. Run `dft doctor` to see what was found.",
        ),
        Param(
            "execution.nprocs", "Parallel processes", "int", 4, unit="cores",
            help="Number of MPI processes ORCA uses per job (the %pal block).",
            guidance="Leave 2-4 cores free for the operating system and for the figure rendering that "
                     "runs alongside. Hybrid DFT with RIJCOSX scales reasonably to about 8-16 cores for a "
                     "30-atom molecule; beyond that the communication overhead eats the gain. Running two "
                     "molecules on 6 cores each usually beats one molecule on 12.",
            minimum=1,
        ),
        Param(
            "execution.maxcore_mb", "Memory per process", "int", 900, unit="MB",
            help="The %maxcore value - memory ORCA may use PER MPI PROCESS, not in total.",
            guidance="This is the single most common cause of ORCA crashes, because it is memory PER "
                     "MPI PROCESS, not per job. Real usage is maxcore_mb x nprocs x max_parallel_jobs - "
                     "so 3000 MB on 8 cores with 2 concurrent jobs asks for 48 GB, not 3 GB. Set it to "
                     "roughly 75% of (usable RAM / (nprocs x max_parallel_jobs)); ORCA treats the value "
                     "as a target rather than a hard cap and will overshoot. `dft validate` checks the "
                     "product against the RAM this machine actually has and warns if it does not fit.",
            minimum=200,
        ),
        Param(
            "execution.scratch_dir", "Scratch directory", "path", "scratch",
            help="Working directory for ORCA's large temporary integral files.",
            guidance="Put it on the fastest local disk you have, ideally an SSD, and never on a network "
                     "mount - ORCA's integral files are written and reread constantly and NFS will "
                     "dominate your runtime. Needs several GB free per concurrent job.",
        ),
        Param(
            "execution.resume", "Skip completed jobs", "bool", True,
            help="If a job's output already exists and terminated normally, reuse it instead of "
                 "recalculating.",
            guidance="Keep on. It makes the pipeline restartable after an interruption and lets you "
                     "iterate on the figures without repeating hours of DFT. The check is for ORCA's own "
                     "'TERMINATED NORMALLY' marker, so a crashed or truncated job is correctly rerun.",
        ),
        Param(
            "execution.max_parallel_jobs", "Concurrent ORCA jobs", "int", 2, unit="jobs",
            help="How many ORCA calculations run at the same time, each using nprocs cores.",
            guidance="Total core usage is max_parallel_jobs x nprocs, which must not exceed your core "
                     "count. For a set of small molecules, 3 jobs x 4 cores finishes sooner than 1 job x "
                     "12 cores because ORCA's parallel efficiency falls off. Watch memory: each job "
                     "separately claims maxcore_mb x nprocs.",
            minimum=1,
        ),
        Param(
            "execution.timeout_min", "Per-job timeout", "int", 480, unit="minutes",
            help="Wall-clock limit after which a single ORCA job is killed and marked failed.",
            guidance="A generous ceiling that catches genuinely stuck jobs - an optimization oscillating "
                     "between two geometries will otherwise run forever. 480 min is ample for a 30-atom "
                     "molecule; raise it for anything over 60 atoms.",
            advanced=True,
        ),
        Param(
            "execution.work_dir", "Output directory", "path", "runs",
            help="Root directory for all generated inputs, outputs, cubes, figures and logs.",
            guidance="One subdirectory per molecule per solvent is created underneath. Keep separate work "
                     "directories for separate projects so the resume logic cannot reuse the wrong output.",
        ),
        Param(
            "execution.dry_run", "Dry run", "bool", False,
            help="Generate every input file and script but do not launch ORCA.",
            guidance="The fastest way to check your settings before committing hours of compute: inspect "
                     "the generated .inp files, confirm the keyword line is what you intended, then rerun "
                     "with this off. Highly recommended before any long campaign.",
        ),
    ],
)

GROUPS = [G_MOLECULES, G_CONFORMERS, G_THEORY, G_SOLVATION, G_PROPERTY,
          G_FUKUI, G_CUBES, G_RENDER, G_FIGURES, G_EXECUTION]

ALL_PARAMS = {p.key: p for g in GROUPS for p in g.params}


def defaults() -> dict:
    """Nested dict of every parameter's default value."""
    cfg: dict = {}
    for key, p in ALL_PARAMS.items():
        section, name = key.split(".", 1)
        cfg.setdefault(section, {})[name] = p.default
    return cfg


def preset(tag: str) -> dict:
    """Config in which every categorical parameter takes the choice carrying ``tag``.

    Parameters with no tagged choice keep their default.
    """
    cfg = defaults()
    for key, p in ALL_PARAMS.items():
        section, name = key.split(".", 1)
        if p.kind == "choice":
            for c in p.choices:
                if tag in c.tags:
                    cfg[section][name] = c.value
                    break
        elif p.kind == "multichoice":
            tagged = [c.value for c in p.choices if tag in c.tags]
            if tagged:
                cfg[section][name] = tagged
    if tag == TAG_PAPER:
        # The manuscript reports vacuum plus four solvents and no dispersion, and its four
        # compounds are neutral - so a single fixed charge is correct there.
        cfg["theory"]["run_freq"] = True
        cfg["render"]["esp_autoscale"] = False
        cfg["render"]["esp_scale"] = "fixed"
        cfg["figures"]["pka_source"] = "csv"
        cfg["figures"]["logd_source"] = "csv"
        cfg["molecules"]["definitions"] = dict(REFERENCE_SET)
        cfg["molecules"]["structure_source"] = "smiles"
        cfg["molecules"]["charge_mode"] = "fixed"
        cfg["molecules"]["charge"] = 0
        cfg["molecules"]["multiplicity_mode"] = "fixed"
        cfg["molecules"]["multiplicity"] = 1
    # A preset is assembled choice by choice, so two independently "recommended" values
    # can contradict each other - wB97X-D4 carries its own dispersion, yet D3BJ is the
    # recommended dispersion for functionals that do not. Resolve that here rather than
    # shipping a preset its own validator rejects.
    func = str(cfg["theory"]["functional"])
    if "D3BJ" in func or "D4" in func or "3c" in func:
        cfg["theory"]["dispersion"] = "none"

    if tag == TAG_FAST:
        cfg["conformers"]["n_confs"] = 10
        cfg["cubes"]["grid_points"] = 40
        cfg["theory"]["run_freq"] = False
        cfg["render"]["ray_trace"] = False
        cfg["render"]["width"] = 600
        cfg["render"]["height"] = 450
        cfg["figures"]["dpi"] = 150
    return cfg
