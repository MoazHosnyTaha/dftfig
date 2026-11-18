"""Structure-derived ionization: which groups ionize, at what pKa, and the logD that follows.

The values are class averages from standard medicinal-chemistry pKa compilations (Perrin,
"Dissociation Constants of Organic Bases in Aqueous Solution"; Serjeant & Dempsey,
"Ionisation Constants of Organic Acids in Aqueous Solution"; and the ranges tabulated in
Manallack, Perspect. Medicin. Chem. 2007, 1, 25-38). They are NOT measurements and they are
NOT a trained predictor: expect +-1 to 1.5 log units, more for a heavily substituted ring.

That accuracy is enough for the question the figures ask - is this compound charged at pH 7.4,
and by how much does that cut its lipophilicity - and not enough to quote in a paper. Anything
in the properties CSV overrides what is predicted here.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

#: Floor on the neutral fraction when converting logP to logD. A permanently charged
#: compound has no neutral microspecies at all; without a floor logD is -infinity.
_MIN_NEUTRAL_FRACTION = 1e-4


@dataclass(frozen=True)
class Site:
    name: str
    smarts: str
    kind: str          # "acid" (loses H+, goes negative) or "base" (gains H+, goes positive)
    pka: float
    anchor: int        # index into the SMARTS match of the atom that ionizes
    note: str = ""
    #: Which atoms this rule takes off the table for later rules. "anchor" claims only the
    #: ionizing atom; "heteroatoms" claims every N/O/S/P in the match, which is what stops a
    #: biguanide matching the guanidine pattern twice and a tetrazole's three remaining
    #: nitrogens being counted as basic azole nitrogens.
    claim: str = "anchor"
    #: True for the arylamine rules: an amino group on a ring that already has a basic
    #: nitrogen does not protonate on the exocyclic N, so the rule must stand down.
    veto_if_azine: bool = False


# Ordered most specific first: the first rule to claim an atom wins, so a sulfonamide
# nitrogen is never also counted as a plain amine.
SITES: tuple[Site, ...] = (
    # ---------------------------------------------------------------- acids
    Site("sulfonic acid", "[SX4](=O)(=O)[OX2H1]", "acid", -1.5, 3,
         "Fully ionized at every physiological pH."),
    Site("phosphonic/phosphate OH", "[PX4](=O)[OX2H1]", "acid", 2.0, 2,
         "First ionization only; the second is around pKa 7 and is not modelled."),
    Site("acyl sulfonamide", "[CX3](=O)[NX3H1][SX4](=O)(=O)", "acid", 4.5, 2,
         "A carboxylic-acid bioisostere and about as acidic."),
    Site("alpha-amino acid carboxyl", "[NX3][CX4][CX3](=O)[OX2H1]", "acid", 2.2, 4,
         "The neighbouring ammonium pulls this down from 4.5; glycine is 2.35."),
    Site("benzoic acid", "c[CX3](=O)[OX2H1]", "acid", 4.2, 3, ""),
    Site("carboxylic acid", "[CX3](=O)[OX2H1]", "acid", 4.5, 2, ""),
    Site("tetrazole (1H)", "[nX3H1]1[nX2][nX2][nX2][cX3]1", "acid", 4.9, 0,
         "Carboxylic-acid bioisostere.", claim="heteroatoms"),
    Site("tetrazole (2H)", "[nX2]1[nX3H1][nX2][nX2][cX3]1", "acid", 4.9, 1, "",
         claim="heteroatoms"),
    Site("thiophenol", "c[SX2H1]", "acid", 6.6, 1, ""),
    Site("hydroxamic acid", "[CX3](=O)[NX3H1][OX2H1]", "acid", 9.0, 3, ""),
    Site("imide", "[NX3H1](C=O)C=O", "acid", 9.5, 0,
         "Two flanking carbonyls make the N-H genuinely acidic; a single amide N-H is not."),
    Site("phenol", "c[OX2H1]", "acid", 10.0, 1,
         "Class average. Electron-withdrawing substituents lower this a long way - "
         "4-nitrophenol is 7.2, pentafluorophenol 5.5."),
    Site("aryl sulfonamide", "c[SX4](=O)(=O)[NX3H1]", "acid", 10.1, 4, ""),
    Site("sulfonamide", "[SX4](=O)(=O)[NX3H1]", "acid", 10.5, 3, ""),
    Site("thiol", "[CX4][SX2H1]", "acid", 10.5, 1, ""),

    # ---------------------------------------------------------------- bases
    Site("biguanide", "[NX3][CX3](=[NX2])[NX3][CX3](=[NX2])[NX3]", "base", 12.4, 2,
         "Metformin and relatives: one strongly basic centre, not two.",
         claim="heteroatoms"),
    Site("guanidine", "[NX3][CX3](=[NX2])[NX3]", "base", 13.6, 2,
         "Protonated at every physiological pH.", claim="heteroatoms"),
    Site("amidine", "[NX3][CX3]=[NX2]", "base", 11.5, 2, "", claim="heteroatoms"),
    Site("morpholine N", "[NX3;H0,H1]1[CX4][CX4][OX2][CX4][CX4]1", "base", 8.4, 0,
         "The ring oxygen pulls this down about 1.5 units from an ordinary tertiary amine.",
         claim="heteroatoms"),
    Site("piperazine (distal N)", "[NX3;H0,H1]1[CX4][CX4][NX3][CX4][CX4]1", "base", 5.6, 3,
         "The second nitrogen of a piperazine; the first stays near 9."),
    Site("amine alpha to a carbonyl",
         "[NX3;H0,H1;!$(Nc);!$(N*=[O,N,S])][CX4][CX3](=[OX1])[!#8]", "base", 7.8, 0,
         "The glycinamide motif of lidocaine and the local anaesthetics; the carbonyl costs "
         "about two units. Excluded when the carbonyl is a carboxylic acid, where the "
         "neighbouring carboxylate raises the amine instead.", claim="anchor"),
    Site("beta-amino alcohol", "[NX3;H1,H2;!$(Nc)][CX4][CX4][OX2H1]", "base", 9.5, 0,
         "The beta-blocker motif; the hydroxyl costs roughly one unit of basicity."),
    Site("primary aliphatic amine",
         "[NX3;H2;!$(N[!#6]);!$(N*=[O,N,S]);!$(Nc)][CX4]", "base", 10.6, 0, ""),
    Site("secondary aliphatic amine",
         "[NX3;H1;!$(N[!#6]);!$(N*=[O,N,S]);!$(Nc)]([CX4])[CX4]", "base", 10.8, 0, ""),
    Site("tertiary aliphatic amine",
         "[NX3;H0;!$(N[!#6]);!$(N*=[O,N,S]);!$(Nc)]([CX4])([CX4])[CX4]", "base", 9.8, 0, ""),
    Site("imidazole (isolated ring)", "[nX2;!R2]1[cX3;!R2][nX3H1;!R2][cX3;!R2][cX3;!R2]1",
         "base", 7.0, 0,
         "The histidine side chain; half protonated near physiological pH. Restricted to an "
         "isolated ring - fusing it to a benzene drops the pKa to 5.5 (benzimidazole) and to "
         "a pyrimidine to 2.4 (purine).", claim="heteroatoms"),
    Site("benzimidazole N", "[nX2]1[cX3][nX3H1][cX3]2[cX3][cX3][cX3][cX3][cX3]12",
         "base", 5.5, 0, "Fusion to benzene costs the isolated imidazole about 1.5 units.",
         claim="heteroatoms"),
    Site("amino-substituted azine", "[nX2;$(n:c[NX3;H1,H2])]", "base", 5.0, 0,
         "An amino group feeding into a pyridine or pyrimidine ring raises the ring "
         "nitrogen's basicity: adenine 4.2, 2-aminopyridine 6.7."),
    Site("N,N-dialkylaniline", "[NX3;H0]([CX4])([CX4])c", "base", 5.1, 0, "",
         veto_if_azine=True),
    Site("N-alkylaniline", "[NX3;H1]([CX4])c", "base", 5.0, 0, "", veto_if_azine=True),
    Site("aniline", "[NX3;H2]c", "base", 4.6, 0, "", veto_if_azine=True),
    Site("pyridine", "[nX2]1[cX3][cX3][cX3][cX3][cX3]1", "base", 5.2, 0,
         "Quinoline 4.9, isoquinoline 5.4, pyridine 5.2 - all within the spread of this value."),
    Site("pyrimidine/pyrazine N", "[nX2]1[cX3][nX2][cX3][cX3][cX3]1", "base", 1.3, 0,
         "A second ring nitrogen removes roughly four units of basicity."),
    Site("azole ring N", "[nX2;r5]", "base", 2.5, 0,
         "Pyrazole 2.5, 1,2,4-triazole 2.2, thiazole 2.5. Not basic at physiological pH."),

)

#: Groups deliberately NOT treated as ionizable in water between pH 0 and 14, with the
#: reason. Printed by `dft params properties.pka_source` so the omissions are visible.
NOT_IONIZABLE = {
    "aliphatic alcohol": "pKa about 16",
    "amide N-H": "pKa about 17",
    "pyrrole / indole / carbazole N-H": "pKa 17 to 20 as an acid, about -6 as a base",
    "diarylamine": "conjugate acid well below pH 0",
    "aromatic C-H, ether, ester, nitrile, halide": "no ionizable proton",
}

#: Classes the model is known to get wrong, so they can be overridden deliberately.
KNOWN_WEAKNESSES = (
    "Vinylogous acids - 4-hydroxycoumarin (warfarin), tetronic and tetramic acids - are "
    "reported as ordinary phenols or alcohols and come out several units too high.",
    "Ring substituent effects are not modelled: a nitro or polyfluoro substituent can move "
    "a phenol or an aniline by 3 units and the class average will not follow it.",
    "Only the first ionization of a polyprotic acid is modelled.",
    "Tautomers are taken as drawn; the model does not search for the most acidic tautomer.",
)

_COMPILED = tuple((s, Chem.MolFromSmarts(s.smarts)) for s in SITES)
_UNCHARGER = rdMolStandardize.Uncharger()


def _prepare(mol: Chem.Mol) -> Chem.Mol:
    """Heavy-atom view of the molecule.

    Explicit hydrogens have to go first. With them in the graph a protonated amine has
    degree 4 and zero implicit hydrogens, which makes it indistinguishable from a
    quaternary ammonium - the molecule is then read as permanently charged and its real
    basic pKa is never found.
    """
    try:
        m = Chem.RemoveHs(Chem.Mol(mol))
        Chem.SanitizeMol(m)
        return m
    except Exception:
        return mol


def _neutral(mol: Chem.Mol) -> Chem.Mol:
    """Strip formal charges so an already-protonated amine still matches the amine rules."""
    try:
        m = _UNCHARGER.uncharge(_prepare(mol))
        Chem.SanitizeMol(m)
        return m
    except Exception:
        return _prepare(mol)


def permanent_charge(mol: Chem.Mol) -> int:
    """Charge that no pH can remove: quaternary nitrogen, and similar.

    A pyrrolidinium N-H is NOT permanent - it is a protonated amine and comes back as a
    basic site. A quaternary ammonium with four carbon substituents is.
    """
    q = 0
    for a in _prepare(mol).GetAtoms():
        if a.GetSymbol() == "N" and a.GetFormalCharge() == 1 and a.GetTotalNumHs() == 0 \
                and a.GetDegree() == 4:
            q += 1
        elif a.GetSymbol() == "B" and a.GetFormalCharge() == -1 and a.GetDegree() == 4:
            q -= 1
    return q


def _ring_systems(mol: Chem.Mol) -> list[set[int]]:
    """Fused aromatic ring systems, as sets of atom indices."""
    systems: list[set[int]] = []
    for ring in mol.GetRingInfo().AtomRings():
        s = set(ring)
        merged = [t for t in systems if t & s]
        for t in merged:
            s |= t
            systems.remove(t)
        systems.append(s)
    return systems


def _attached_ring_is_azine(mol: Chem.Mol, n_idx: int, systems) -> bool:
    """True if the arylamine's ring system carries its own aromatic nitrogen.

    An amino group on a pyridine or pyrimidine does not protonate on the exocyclic
    nitrogen - the ring nitrogen takes the proton, and does so more readily because the
    amino group feeds electron density into the ring. Counting both would give the
    molecule two basic centres where it has one.
    """
    atom = mol.GetAtomWithIdx(n_idx)
    for nb in atom.GetNeighbors():
        if not nb.GetIsAromatic():
            continue
        for sysset in systems:
            if nb.GetIdx() in sysset and any(
                    mol.GetAtomWithIdx(i).GetSymbol() == "N"
                    and mol.GetAtomWithIdx(i).GetIsAromatic() for i in sysset):
                return True
    return False


def predict_sites(mol: Chem.Mol) -> list[dict]:
    """Every ionizable centre found, as dicts with atom index, kind, pKa and the rule name."""
    if mol is None:
        return []
    neutral = _neutral(mol)
    systems = _ring_systems(neutral)
    claimed: set[int] = set()
    found: list[dict] = []
    for site, patt in _COMPILED:
        if patt is None:
            continue
        for match in neutral.GetSubstructMatches(patt, uniquify=True):
            if site.anchor >= len(match):
                continue
            idx = match[site.anchor]
            if idx in claimed:
                continue
            if site.veto_if_azine and _attached_ring_is_azine(neutral, idx, systems):
                continue
            if site.claim == "heteroatoms":
                hetero = {i for i in match
                          if neutral.GetAtomWithIdx(i).GetSymbol() in ("N", "O", "S", "P")}
                if hetero & claimed:      # overlaps a group already assigned - e.g. the
                    continue              # second guanidine of a biguanide
            if site.claim == "heteroatoms":
                claimed.update(i for i in match
                               if neutral.GetAtomWithIdx(i).GetSymbol() in ("N", "O", "S", "P"))
            else:
                claimed.add(idx)
            found.append({"atom": idx, "kind": site.kind, "pKa": site.pka,
                          "rule": site.name, "note": site.note})

    # One basic centre per fused aromatic ring system. A purine has three aromatic
    # nitrogens that each match a rule, but protonating one deactivates the rest; keeping
    # them all would put the molecule at charge +3 in acid.
    best: dict[int, dict] = {}
    keep: list[dict] = []
    for s in found:
        atom = neutral.GetAtomWithIdx(s["atom"])
        if s["kind"] != "base" or not atom.GetIsAromatic():
            keep.append(s)
            continue
        sys_id = next((k for k, sysset in enumerate(systems) if s["atom"] in sysset), None)
        if sys_id is None:
            keep.append(s)
        elif sys_id not in best or s["pKa"] > best[sys_id]["pKa"]:
            best[sys_id] = s
    keep.extend(best.values())
    keep.sort(key=lambda d: (d["kind"], -d["pKa"]))
    return keep


def predict(mol: Chem.Mol) -> dict:
    """Acidic and basic pKa lists plus the permanent charge, for one molecule."""
    sites = predict_sites(mol)
    return {
        "pKa_acidic": [s["pKa"] for s in sites if s["kind"] == "acid"],
        "pKa_basic": [s["pKa"] for s in sites if s["kind"] == "base"],
        "sites": sites,
        "permanent_charge": permanent_charge(mol) if mol is not None else 0,
    }


def neutral_fraction(ph: float, acidic, basic) -> float:
    """Fraction of the molecule in its fully neutral microspecies (independent sites)."""
    f = 1.0
    for pka in acidic:
        f /= 1.0 + 10.0 ** (ph - pka)
    for pka in basic:
        f /= 1.0 + 10.0 ** (pka - ph)
    return float(f)


def logd_from_logp(logp: float, ph: float, acidic, basic, permanent_charge: int = 0):
    """logD at a given pH from logP and the ionization model.

    Only the neutral microspecies is taken to partition into octanol, which is the standard
    approximation. Ion pairing makes the real logD of a strong base roughly 1 unit higher
    than this, so treat the value as a lower bound for amines.
    """
    if logp is None or (isinstance(logp, float) and np.isnan(logp)):
        return np.nan, False
    f = neutral_fraction(ph, acidic, basic)
    if permanent_charge:
        f = 0.0
    floored = f < _MIN_NEUTRAL_FRACTION
    f = max(f, _MIN_NEUTRAL_FRACTION)
    return float(logp + np.log10(f)), floored
