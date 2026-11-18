"""Physicochemical properties, CNS MPO scoring and acid-base speciation.

Structures are the only required input. MW, TPSA, HBD and Crippen logP come straight from
the connection table; pKa comes from the SMARTS ionization model in `pka.py`; logD follows
from logP and that ionization. Each is labelled with where it came from in
`properties_resolved.csv`, and anything supplied in the properties CSV overrides it.

Quantities that genuinely require a trained ADMET model - permeability, clearance, plasma
protein binding, P-gp - are NOT invented here. They stay blank unless you supply them, and
the figures say so rather than plotting a guess.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from . import pka
from .theme import note, warn


# ----------------------------------------------------------------- CNS MPO (Wager 2010)
def _ramp(x: float, one_below: float, zero_above: float) -> float:
    """Desirability 1 below `one_below`, falling linearly to 0 at `zero_above`."""
    if x <= one_below:
        return 1.0
    if x >= zero_above:
        return 0.0
    return (zero_above - x) / (zero_above - one_below)


def _hump(x: float, rise0: float, rise1: float, fall0: float, fall1: float) -> float:
    """Desirability 0 below rise0, 1 on [rise1, fall0], 0 above fall1."""
    if x <= rise0 or x >= fall1:
        return 0.0
    if x < rise1:
        return (x - rise0) / (rise1 - rise0)
    if x <= fall0:
        return 1.0
    return (fall1 - x) / (fall1 - fall0)


#: The six Wager (2010) desirability functions, ACS Chem. Neurosci. 1, 435-449.
MPO_FUNCTIONS = {
    "clogP":  lambda v: _ramp(v, 3.0, 5.0),
    "clogD":  lambda v: _ramp(v, 2.0, 4.0),
    "MW":     lambda v: _ramp(v, 360.0, 500.0),
    "TPSA":   lambda v: _hump(v, 20.0, 40.0, 90.0, 120.0),
    "HBD":    lambda v: _ramp(v, 0.5, 3.5),
    "pKa":    lambda v: _ramp(v, 8.0, 10.0),
}

MPO_LABELS = {
    "clogP": "intrinsic lipophilicity (clogP)",
    "clogD": "buffer lipophilicity (clogD pH 7.4)",
    "MW": "molecular weight",
    "TPSA": "topological polar surface area",
    "HBD": "hydrogen-bond donors",
    "pKa": "most basic pKa",
}


def cns_mpo(clogP: float, clogD: float, MW: float, TPSA: float, HBD: float,
            pKa_basic: float) -> dict:
    """CNS MPO: six desirability terms summed to a 0-6 score, plus the breakdown."""
    vals = {"clogP": clogP, "clogD": clogD, "MW": MW, "TPSA": TPSA,
            "HBD": HBD, "pKa": pKa_basic}
    terms = {}
    for k, f in MPO_FUNCTIONS.items():
        v = vals[k]
        terms[k] = float("nan") if v is None or (isinstance(v, float) and np.isnan(v)) else f(float(v))
    score = float(np.nansum(list(terms.values())))
    return {"cns_mpo": score, **{f"d_{k}": v for k, v in terms.items()}}


# -------------------------------------------------------------------- acid-base speciation
def parse_pka_list(raw) -> list[float]:
    if raw is None or (isinstance(raw, float) and np.isnan(raw)):
        return []
    s = str(raw).strip()
    if not s or s.lower() in ("none", "nan", "-"):
        return []
    out = []
    for tok in s.replace(",", ";").split(";"):
        tok = tok.strip()
        if tok:
            try:
                out.append(float(tok))
            except ValueError:
                warn(f"ignoring unparsable pKa value {tok!r}")
    return out


def net_charge(ph, acidic: list[float], basic: list[float],
               baseline_charge: int = 0) -> np.ndarray:
    """Net charge vs pH under the independent-site Henderson-Hasselbalch approximation.

    Each acidic site (HA -> A- + H+) contributes -1 x 1/(1+10^(pKa-pH)).
    Each basic site (B + H+ -> BH+) contributes +1 x 1/(1+10^(pH-pKa)).

    `baseline_charge` is the charge that no pH can change - a quaternary ammonium, a
    sulfonium, a borate. It is NOT the formal charge of the structure you drew: a
    pyrrolidine drawn as its hydrochloride salt has formal charge +1, but that +1 is
    exactly what the basic pKa term already describes, and adding both counts it twice.
    When no pKa is known at all the caller passes the formal charge here and the curve
    stays flat at it, which is the honest statement that nothing was modelled.

    Sites are treated as independent, which is the usual assumption behind a predicted
    net-charge curve; strongly coupled ionizations need explicit microspeciation.
    """
    ph = np.atleast_1d(np.asarray(ph, dtype=float))
    q = np.full_like(ph, float(baseline_charge))
    # Not `pka`: that name is the pka MODULE at the top of this file, and shadowing it
    # inside a loop is how a later edit here would reach for pka.predict and get a float.
    for pka_acid in acidic:
        q -= 1.0 / (1.0 + 10.0 ** (pka_acid - ph))
    for pka_base in basic:
        q += 1.0 / (1.0 + 10.0 ** (ph - pka_base))
    return q


def dominant_species(ph: float, acidic: list[float], basic: list[float],
                     baseline_charge: int = 0) -> str:
    q = float(net_charge(ph, acidic, basic, baseline_charge)[0])
    if q < -0.5:
        return "anion"
    if q > 0.5:
        return "cation"
    if abs(q) < 0.05:
        return "neutral"
    return "mostly neutral"


# ---------------------------------------------------------------------------- CSV loading
REQUIRED = ["molecule"]
OPTIONAL_NUMERIC = ["logP", "logD74", "TPSA", "MW", "HBD", "pKa_basic_main",
                    "Papp_Caco2", "Papp_MDCK", "CLint_u", "PPB", "logBB", "HIA",
                    "Pgp_substrate"]


def load_properties(cfg: dict, rdkit_props: dict[str, dict] | None = None,
                    mols: dict | None = None) -> pd.DataFrame:
    """Read the properties CSV and fill structure-derived columns from RDKit."""
    path = Path(cfg["figures"]["properties_csv"])
    names = list(rdkit_props) if rdkit_props else list(cfg["molecules"]["definitions"])
    if not path.exists() and cfg["figures"].get("pka_source", "predict") == "csv":
        # The packaged template holds the manuscript's published values for CCB1-4. Fall
        # back to it ONLY when the config has asked for CSV-sourced values and the file is
        # missing - never when the pipeline is predicting, where silently substituting
        # someone else's numbers for a name match would be indefensible.
        packaged = Path(__file__).parent / "templates" / "properties.csv"
        if packaged.exists():
            template = pd.read_csv(packaged, comment="#")
            if set(names) & set(template["molecule"]):
                warn(f"{path} not found - using the properties template shipped with dftfig "
                     f"(the manuscript's published values for the reference compounds). "
                     f"Run `dft init` to copy an editable copy into this directory.")
                path = packaged
    if path.exists():
        # comment="#": the skeleton `dft init` writes carries an explanatory
        # header, and users annotate their own files the same way.
        df = pd.read_csv(path, comment="#")
        df.columns = [c.strip() for c in df.columns]
    else:
        warn(f"{path} not found - using structure-derived properties only; panels that need "
             "pKa, logD or permeability will say so on the figure")
        df = pd.DataFrame({"molecule": names})

    df = df[df["molecule"].isin(names)].copy()
    missing = [n for n in names if n not in set(df["molecule"])]
    for n in missing:
        df = pd.concat([df, pd.DataFrame([{"molecule": n}])], ignore_index=True)
    df["_order"] = df["molecule"].map({n: i for i, n in enumerate(names)})
    df = df.sort_values("_order").drop(columns="_order").reset_index(drop=True)

    for col in OPTIONAL_NUMERIC + ["pKa_acidic", "pKa_basic"]:
        if col not in df.columns:
            df[col] = np.nan

    if rdkit_props:
        for i, row in df.iterrows():
            rp = rdkit_props.get(row["molecule"])
            if not rp:
                continue
            for src, dst in [("MW", "MW"), ("TPSA", "TPSA"), ("HBD", "HBD"),
                             ("clogP_crippen", "clogP_rdkit")]:
                if dst not in df.columns:
                    df[dst] = np.nan
                if pd.isna(df.at[i, dst]) if dst in df.columns else True:
                    df.at[i, dst] = rp[src]
            # The column has to exist before .at writes into it: pandas 4 raises KeyError
            # instead of creating it, and a user CSV without a `formula` column is the
            # normal case, not the exception.
            if "formula" not in df.columns:
                df["formula"] = ""
            df.at[i, "formula"] = rp["formula"]
            df.at[i, "clogP_rdkit"] = rp["clogP_crippen"]

    # ---------------------------------------------------------------- ionization
    # Nothing beyond the structures is required: any pKa left blank is predicted from the
    # SMARTS model in pka.py, and logD follows from clogP and that ionization. Values given
    # in the CSV always win - the prediction only fills gaps.
    pka_source = cfg["figures"].get("pka_source", "predict")
    logd_source = cfg["figures"].get("logd_source", "derive")
    ph74 = float(cfg["figures"].get("physiological_ph", 7.4))
    for col in ("pKa_source", "logD74_source", "ionizable_groups"):
        if col not in df.columns:
            df[col] = ""
    if "permanent_charge" not in df.columns:
        df["permanent_charge"] = 0
    # pKa lists are written back as strings ("4.9; 10.0"); an all-NaN column arrives as
    # float64 and pandas refuses a string into it.
    for col in ("pKa_acidic", "pKa_basic", "pKa_source", "logD74_source",
                "ionizable_groups"):
        df[col] = df[col].astype(object)
    df["permanent_charge"] = pd.to_numeric(df["permanent_charge"],
                                           errors="coerce").fillna(0).astype(int)
    df["logD74"] = pd.to_numeric(df["logD74"], errors="coerce")

    predicted_any, floored_any = [], []
    for i, row in df.iterrows():
        name = row["molecule"]
        mol = (mols or {}).get(name)
        acidic = parse_pka_list(row.get("pKa_acidic"))
        basic = parse_pka_list(row.get("pKa_basic"))
        have_csv = bool(acidic or basic)
        df.at[i, "pKa_source"] = "csv" if have_csv else "none"

        if pka_source == "predict" and not have_csv and mol is not None:
            pred = pka.predict(mol)
            acidic, basic = pred["pKa_acidic"], pred["pKa_basic"]
            df.at[i, "pKa_acidic"] = "; ".join(f"{v:.1f}" for v in acidic)
            df.at[i, "pKa_basic"] = "; ".join(f"{v:.1f}" for v in basic)
            df.at[i, "permanent_charge"] = pred["permanent_charge"]
            df.at[i, "ionizable_groups"] = "; ".join(
                f"{s['rule']} {s['pKa']:.1f}" for s in pred["sites"])
            df.at[i, "pKa_source"] = "predicted"
            if acidic or basic or pred["permanent_charge"]:
                predicted_any.append(name)

        if pd.isna(row.get("logD74")):
            df.at[i, "logD74_source"] = "none"
            if logd_source == "derive":
                logp = row.get("logP")
                if pd.isna(logp):
                    logp = df.at[i, "clogP_rdkit"] if "clogP_rdkit" in df.columns else np.nan
                val, floored = pka.logd_from_logp(
                    logp, ph74, acidic, basic, int(df.at[i, "permanent_charge"] or 0))
                if not pd.isna(val):
                    df.at[i, "logD74"] = val
                    df.at[i, "logD74_source"] = "derived"
                    if floored:
                        floored_any.append(name)
        else:
            df.at[i, "logD74_source"] = "csv"

    if predicted_any:
        note(f"pKa predicted from structure for {', '.join(predicted_any)} "
             f"(SMARTS class averages, about +-1 log unit - put measured values in "
             f"{cfg['figures']['properties_csv']} to override)")
    if floored_any:
        warn(f"{', '.join(floored_any)} carry a permanent charge or are fully ionized at "
             f"pH {ph74:g}: the neutral fraction is below 1e-4, so logD is reported at that "
             f"floor and is a lower bound, not a prediction")

    # CNS MPO from whichever logP is available.
    rows = []
    for _, row in df.iterrows():
        clogP = row.get("logP")
        if pd.isna(clogP):
            clogP = row.get("clogP_rdkit", np.nan)
        clogD = row.get("logD74", np.nan)
        pka_b = row.get("pKa_basic_main", np.nan)
        if pd.isna(pka_b):
            bl = parse_pka_list(row.get("pKa_basic"))
            pka_b = max(bl) if bl else np.nan
        rows.append(cns_mpo(clogP, clogD, row.get("MW", np.nan), row.get("TPSA", np.nan),
                            row.get("HBD", np.nan), pka_b))
    df = pd.concat([df.reset_index(drop=True), pd.DataFrame(rows)], axis=1)
    return df


def write_template(path: Path, names: list[str]) -> Path:
    """Write a properties CSV skeleton with every column the figures can use."""
    cols = (["molecule", "pKa_acidic", "pKa_basic", "pKa_basic_main", "logP", "logD74",
             "TPSA", "MW", "HBD"] + OPTIONAL_NUMERIC[6:])
    df = pd.DataFrame({c: ["" for _ in names] for c in cols})
    df["molecule"] = names
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    return path
