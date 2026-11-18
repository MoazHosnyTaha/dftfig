"""Validation of the SMARTS pKa model against measured values.

dftfig does not compute pKa - it assigns a class average to each ionizable group it
recognises (see dftfig/pka.py). The claim that goes with it is "+-1 to 1.5 log units for
the dominant ionization", and that claim is only worth making if it is checked, so the
reference set lives here as a test rather than in a paragraph of the README.

Experimental values are the standard compilations:

  * D. D. Perrin, "Dissociation Constants of Organic Bases in Aqueous Solution",
    Butterworths, London, 1965 (and the 1972 supplement).
  * E. P. Serjeant, B. Dempsey, "Ionisation Constants of Organic Acids in Aqueous
    Solution", IUPAC Chemical Data Series No. 23, Pergamon, Oxford, 1979.
  * D. T. Manallack, "The pKa Distribution of Drugs: Application to Drug Discovery",
    Perspect. Medicin. Chem. 2007, 1, 25-38 - the tabulated ranges per functional class.

Each compound is checked on its DOMINANT ionization only: the lowest acidic pKa (the
first proton to leave) or the highest basic pKa (the centre that protonates first).
Minor, far-off centres - losartan's azole nitrogen at 2.5, for instance - are not part
of the claim and are not asserted.

TOLERANCE: 1.5 log units, the stated accuracy of the model. It is NOT to be widened to
make a compound pass. A compound the class averages cannot reach is marked xfail with
the reason, which is the honest way to ship a known limitation.
"""
from __future__ import annotations

import pytest
from rdkit import Chem

from dftfig import pka

#: 1.5 log units is the accuracy claimed in the pka module docstring and in the README.
TOLERANCE = 1.5

#: (name, SMILES, experimental pKa, kind)
REFERENCE = [
    ("aspirin",          "CC(=O)Oc1ccccc1C(=O)O",                      3.5,  "acid"),
    ("ibuprofen",        "CC(C)Cc1ccc(cc1)C(C)C(=O)O",                 4.4,  "acid"),
    ("phenol",           "Oc1ccccc1",                                 10.0,  "acid"),
    ("paracetamol",      "CC(=O)Nc1ccc(O)cc1",                         9.5,  "acid"),
    ("losartan",         "CCCCc1nc(Cl)c(CO)n1Cc1ccc(cc1)-c1ccccc1-"
                         "c1nnn[nH]1",                                 4.9,  "acid"),
    ("diphenhydramine",  "CN(C)CCOC(c1ccccc1)c1ccccc1",                9.0,  "base"),
    ("propranolol",      "CC(C)NCC(O)COc1cccc2ccccc12",                9.5,  "base"),
    ("atenolol",         "CC(C)NCC(O)COc1ccc(CC(N)=O)cc1",             9.6,  "base"),
    ("pyridine",         "c1ccncc1",                                   5.2,  "base"),
    ("quinoline",        "c1ccc2ncccc2c1",                             4.9,  "base"),
    ("imidazole",        "c1c[nH]cn1",                                 7.0,  "base"),
    ("benzimidazole",    "c1ccc2[nH]cnc2c1",                           5.5,  "base"),
    ("adenine",          "Nc1ncnc2[nH]cnc12",                          4.2,  "base"),
    ("metformin",        "CN(C)C(=N)NC(N)=N",                         12.4,  "base"),
    ("lidocaine",        "CCN(CC)CC(=O)Nc1c(C)cccc1C",                 7.9,  "base"),
]

#: Compounds the class-average model is known to get wrong, kept in the suite so the
#: failure is documented and measured rather than discovered by a user.
KNOWN_FAILURES = {
    "warfarin": (
        "CC(=O)CC(c1ccccc1)C1=C(O)c2ccccc2OC1=O", 5.0, "acid",
        "The 4-hydroxycoumarin of warfarin is a VINYLOGOUS acid: the enol sits in "
        "conjugation with the lactone carbonyl, so the anion is delocalised onto that "
        "carbonyl oxygen and the measured pKa is 5.0. The SMARTS `c[OX2H1]` sees only "
        "an aromatic C-OH and assigns the phenol average of 10.0 - about 5 units too "
        "high, the largest error in the reference set. Fixing it needs a rule for the "
        "vinylogous-acid motif (also tetronic and tetramic acids), not a wider "
        "tolerance; it is listed in pka.KNOWN_WEAKNESSES."
    ),
}


def dominant_pka(mol, kind: str):
    """The pKa that dominates ionization: lowest acidic, or highest basic."""
    pred = pka.predict(mol)
    values = pred["pKa_acidic"] if kind == "acid" else pred["pKa_basic"]
    if not values:
        return None
    return min(values) if kind == "acid" else max(values)


@pytest.mark.parametrize("name,smiles,experimental,kind",
                         REFERENCE, ids=[r[0] for r in REFERENCE])
def test_dominant_pka_is_within_tolerance(name, smiles, experimental, kind):
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, f"{name}: unparseable reference SMILES"
    got = dominant_pka(mol, kind)
    assert got is not None, f"{name}: no {kind}ic site found at all"
    assert abs(got - experimental) <= TOLERANCE, (
        f"{name}: predicted {got:.1f}, experimental {experimental:.1f} "
        f"({got - experimental:+.1f} log units)"
    )


@pytest.mark.parametrize("name,smiles,experimental,kind",
                         REFERENCE, ids=[r[0] for r in REFERENCE])
def test_the_dominant_site_is_of_the_expected_kind(name, smiles, experimental, kind):
    """A base predicted as an acid would flip the sign of the charge at pH 7.4."""
    pred = pka.predict(Chem.MolFromSmiles(smiles))
    assert pred["pKa_acidic"] if kind == "acid" else pred["pKa_basic"]


@pytest.mark.xfail(strict=True, reason=KNOWN_FAILURES["warfarin"][3])
def test_warfarin_vinylogous_acid():
    smiles, experimental, kind, _ = KNOWN_FAILURES["warfarin"]
    got = dominant_pka(Chem.MolFromSmiles(smiles), kind)
    assert abs(got - experimental) <= TOLERANCE


def test_warfarin_failure_is_exactly_the_documented_one():
    """Pin the size and cause of the known failure, so a change to it is noticed."""
    mol = Chem.MolFromSmiles(KNOWN_FAILURES["warfarin"][0])
    sites = pka.predict_sites(mol)
    assert [s["rule"] for s in sites] == ["phenol"]
    assert dominant_pka(mol, "acid") == pytest.approx(10.0)
    assert any("Vinylogous" in w for w in pka.KNOWN_WEAKNESSES)


def test_carbazole_has_no_ionizable_group_between_ph_0_and_14():
    """The carbazole N-H is pKa ~20 as an acid and far below 0 as a base.

    Reporting it as ionizable would put a charge on every carbazole-containing compound
    in the logD panel, which is where this check comes from.
    """
    mol = Chem.MolFromSmiles("c1ccc2c(c1)[nH]c1ccccc12")
    pred = pka.predict(mol)
    assert pred["pKa_acidic"] == []
    assert pred["pKa_basic"] == []
    assert pred["permanent_charge"] == 0
    assert pka.neutral_fraction(0.0, [], []) == 1.0
    assert pka.neutral_fraction(14.0, [], []) == 1.0


# ------------------------------------------------------------------ the model's edges
def test_a_protonated_amine_is_still_read_as_a_basic_site():
    """Drawn protonated or drawn neutral, it is the same compound with the same pKa."""
    neutral = pka.predict(Chem.MolFromSmiles("CCN(CC)CC"))
    cation = pka.predict(Chem.MolFromSmiles("CC[NH+](CC)CC"))
    assert neutral["pKa_basic"] == cation["pKa_basic"]
    assert cation["permanent_charge"] == 0


def test_a_quaternary_ammonium_is_permanently_charged():
    pred = pka.predict(Chem.MolFromSmiles("C[N+](C)(C)CCO"))
    assert pred["permanent_charge"] == 1
    assert pred["pKa_basic"] == []


def test_one_basic_centre_per_fused_aromatic_system():
    """Protonating one ring nitrogen deactivates the rest; adenine is +1 in acid, not +3."""
    pred = pka.predict(Chem.MolFromSmiles("Nc1ncnc2[nH]cnc12"))
    assert len(pred["pKa_basic"]) == 1


def test_a_biguanide_is_one_basic_centre_not_two():
    pred = pka.predict(Chem.MolFromSmiles("CN(C)C(=N)NC(N)=N"))
    assert pred["pKa_basic"] == [12.4]


def test_neutral_fraction_is_a_half_at_the_pka():
    assert pka.neutral_fraction(4.5, [4.5], []) == pytest.approx(0.5)
    assert pka.neutral_fraction(9.8, [], [9.8]) == pytest.approx(0.5)


def test_logd_is_below_logp_for_an_ionized_acid():
    logd, floored = pka.logd_from_logp(3.0, 7.4, [4.4], [])
    assert logd < 3.0
    assert not floored


def test_logd_of_a_permanently_charged_compound_is_floored_not_minus_infinity():
    logd, floored = pka.logd_from_logp(1.0, 7.4, [], [], permanent_charge=1)
    assert floored is True
    assert logd == pytest.approx(1.0 - 4.0)      # the 1e-4 neutral-fraction floor


def test_logd_equals_logp_for_a_neutral_compound():
    logd, floored = pka.logd_from_logp(2.5, 7.4, [], [])
    assert logd == pytest.approx(2.5)
    assert not floored
