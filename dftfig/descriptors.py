"""Global reactivity descriptors and condensed Fukui functions.

Conventions, stated explicitly because the literature is not consistent:

  IE  ionization energy          Koopmans: -e(HOMO)    delta-SCF: E(N-1) - E(N)
  EA  electron affinity          Koopmans: -e(LUMO)    delta-SCF: E(N) - E(N+1)
  chi electronegativity          (IE + EA) / 2
  mu  chemical potential         -chi
  eta chemical hardness          (IE - EA) / 2   [Parr-Pearson]
      `eta_full` = IE - EA is also reported, since some papers use that definition and
      the resulting electrophilicity index differs by a factor of two.
  S   chemical softness          1 / eta
  w   electrophilicity index     mu^2 / (2 eta)
  eps nucleophilicity index      1 / w

Condensed Fukui functions from atomic charges q (finite difference, frozen geometry):

  f+_k = q_k(N)   - q_k(N+1)      susceptibility to nucleophilic attack
  f-_k = q_k(N-1) - q_k(N)        susceptibility to electrophilic attack
  f0_k = (f+ + f-) / 2            radical attack
  df_k = f+_k - f-_k              dual descriptor   (>0 electrophilic site)
  dw_k = w * df_k                 multiphilic descriptor
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from . import molecules as _M
from . import parsers
from .orca_io import DIELECTRIC
from .parsers import HARTREE_TO_EV
from .theme import bad, ok, warn


def global_descriptors(ie: float, ea: float) -> dict:
    """All global reactivity indices from an ionization energy and electron affinity in eV."""
    chi = (ie + ea) / 2.0
    mu = -chi
    eta = (ie - ea) / 2.0
    eta_full = ie - ea
    out = {"IE_eV": ie, "EA_eV": ea, "chi_eV": chi, "mu_eV": mu,
           "eta_eV": eta, "eta_full_eV": eta_full}
    if abs(eta) > 1e-9:
        out["S_eV-1"] = 1.0 / eta
        out["omega_eV"] = mu ** 2 / (2.0 * eta)
        out["epsilon_eV-1"] = 1.0 / out["omega_eV"] if abs(out["omega_eV"]) > 1e-12 else float("nan")
    else:
        out["S_eV-1"] = float("nan")
        out["omega_eV"] = float("nan")
        out["epsilon_eV-1"] = float("nan")
    return out


def collect(cfg: dict) -> pd.DataFrame:
    """Build the descriptor table over every molecule x medium that completed."""
    work = Path(cfg["execution"]["work_dir"])
    mode = cfg["property"]["koopmans_mode"]
    _restarts = _saddle_restarts(work)
    _flat = _flat_modes(work)
    rows = []

    for name in _M.names(cfg):
        for solvent in cfg["solvation"]["solvents"]:
            d = work / name / solvent
            sp = d / f"{name}_{solvent}_sp.out"
            if not sp.exists() or not parsers.terminated_normally(sp):
                warn(f"{name}/{solvent}: no completed property single point - skipped")
                continue

            hl = parsers.homo_lumo(sp)
            row = {
                "molecule": name, "medium": solvent,
                "dielectric": DIELECTRIC.get(solvent, float("nan")),
                "E_sp_Eh": parsers.final_energy(sp),
                "dipole_D": parsers.dipole(sp),
                **hl,
            }

            opt = d / f"{name}_{solvent}_opt.out"
            if opt.exists() and parsers.terminated_normally(opt):
                freqs = parsers.frequencies(opt)
                imag = parsers.imaginary_modes(opt, float(cfg["theory"]["imaginary_cutoff"]))
                row["n_imaginary"] = len(imag)
                row["lowest_freq_cm-1"] = min((f for f in freqs if abs(f) > 1e-6), default=float("nan"))
                row["is_minimum"] = (len(imag) == 0) if freqs else None
                row.update(parsers.thermochemistry(opt))

            # --- Koopmans ----------------------------------------------------------------
            if mode in ("koopmans", "both") and hl:
                k = global_descriptors(-hl["e_homo_eV"], -hl["e_lumo_eV"])
                row.update({f"{key}_koopmans": v for key, v in k.items()})

            # --- delta-SCF ---------------------------------------------------------------
            e_n = parsers.final_energy(sp)
            cat, ani = d / f"{name}_{solvent}_cat.out", d / f"{name}_{solvent}_ani.out"
            have_ions = all(p.exists() and parsers.terminated_normally(p) for p in (cat, ani))
            if mode in ("delta_scf", "both") and have_ions and e_n is not None:
                e_cat, e_ani = parsers.final_energy(cat), parsers.final_energy(ani)
                ie = (e_cat - e_n) * HARTREE_TO_EV
                ea = (e_n - e_ani) * HARTREE_TO_EV
                k = global_descriptors(ie, ea)
                row.update({f"{key}_deltascf": v for key, v in k.items()})

            # --- sanity checks that decide whether the numbers mean anything ---------
            notes = []
            if hl and hl.get("e_lumo_eV", -1) > 0:
                notes.append("LUMO above vacuum level")
            k_ea = row.get("EA_eV_koopmans")
            d_ea = row.get("EA_eV_deltascf")
            if k_ea is not None and d_ea is not None and k_ea * d_ea < 0:
                notes.append("Koopmans and delta-SCF electron affinities disagree in sign")
            if row.get("n_imaginary"):
                lo = abs(float(row.get("lowest_freq_cm-1") or 0.0))
                if _flat.get((row.get("molecule"), row.get("medium"))):
                    notes.append(f"residual flat mode at -{lo:.0f} cm-1 - a free internal "
                                 f"rotation, not a reaction saddle; use a hindered-rotor "
                                 f"treatment for the thermochemistry")
                else:
                    notes.append(f"{row['n_imaginary']} imaginary mode(s) - not a minimum")
            # A structure that needed restarts is a valid minimum, but the reader should
            # know the published geometry is not the one the first optimization returned.
            r_n = _restarts.get((row.get("molecule"), row.get("medium")))
            if r_n:
                row["saddle_restarts"] = r_n
                notes.append(f"minimum reached after {r_n} saddle-point restart"
                             f"{'s' if r_n > 1 else ''}")
            else:
                row["saddle_restarts"] = 0
            row["flat_mode"] = bool(_flat.get((row.get("molecule"), row.get("medium"))))
            row["warnings"] = "; ".join(notes)
            rows.append(row)

    df = pd.DataFrame(rows)
    if not df.empty:
        out = work / "descriptors.csv"
        df.to_csv(out, index=False, float_format="%.6f")
        ok(f"descriptors for {len(df)} molecule/medium combinations -> {out}")

        for _, r in df.iterrows():
            if not r.get("warnings"):
                continue
            bad(f"{r['molecule']}/{r['medium']}: {r['warnings']}")

        if (df.get("e_lumo_eV", pd.Series(dtype=float)) > 0).any():
            warn("A positive LUMO energy means the orbital lies above the vacuum level: the "
                 "electron is not bound. With diffuse ('++'/'aug-') functions the lowest "
                 "virtuals are discretised continuum states, not valence antibonding "
                 "orbitals - so the HOMO-LUMO gap is not an excitation energy and the "
                 "Koopmans electron affinity is meaningless. Use the delta-SCF columns, "
                 "drop the diffuse functions for the orbital picture, or run TD-DFT for a "
                 "real excitation energy.")
    return df


def fukui(cfg: dict) -> dict[str, pd.DataFrame]:
    """Condensed Fukui functions per molecule, one CSV each. Keyed by molecule name."""
    work = Path(cfg["execution"]["work_dir"])
    scheme = cfg["property"]["population"]
    ref = cfg["figures"]["reference_solvent"]
    tables: dict[str, pd.DataFrame] = {}

    desc = work / "descriptors.csv"
    omega_lookup: dict[str, float] = {}
    if desc.exists():
        dd = pd.read_csv(desc)
        col = "omega_eV_koopmans" if "omega_eV_koopmans" in dd.columns else (
            "omega_eV_deltascf" if "omega_eV_deltascf" in dd.columns else None)
        if col:
            sub = dd[dd["medium"] == ref]
            omega_lookup = dict(zip(sub["molecule"], sub[col]))

    for name in _M.names(cfg):
        d = work / name / ref
        files = {s: d / f"{name}_{ref}_{s}.out" for s in ("sp", "cat", "ani")}
        if not all(p.exists() and parsers.terminated_normally(p) for p in files.values()):
            warn(f"{name}: Fukui needs completed sp/cat/ani single points in {ref} - skipped")
            continue

        q_n = parsers.populations(files["sp"], scheme)
        q_cat = parsers.populations(files["cat"], scheme)
        q_ani = parsers.populations(files["ani"], scheme)
        if not (q_n and q_cat and q_ani):
            warn(f"{name}: {scheme} charges missing from one of the ion calculations - skipped")
            continue

        syms = parsers.atom_symbols(files["sp"])
        omega = float(omega_lookup.get(name, float("nan")))
        rows = []
        for i in sorted(q_n):
            if i not in q_cat or i not in q_ani:
                continue
            fp = q_n[i] - q_ani[i]
            fm = q_cat[i] - q_n[i]
            df_ = fp - fm
            rows.append({
                "atom_index_0based": i,
                "atom_index_1based": i + 1,
                "element": syms[i] if i < len(syms) else "?",
                f"q_neutral_{scheme}": q_n[i],
                f"q_cation_{scheme}": q_cat[i],
                f"q_anion_{scheme}": q_ani[i],
                "f_plus": fp, "f_minus": fm, "f_zero": (fp + fm) / 2.0,
                "dual_descriptor": df_, "multiphilic_descriptor": omega * df_,
            })
        t = pd.DataFrame(rows)
        out = work / name / f"fukui_{name}_{ref}.csv"
        t.to_csv(out, index=False, float_format="%.6f")
        tables[name] = t
        top_e = t.nlargest(3, "f_minus")["atom_index_1based"].tolist()
        top_n = t.nlargest(3, "f_plus")["atom_index_1based"].tolist()
        ok(f"{name}: Fukui -> {out.name}  (most electrophilic sites: atoms {top_n}; "
           f"most nucleophilic: atoms {top_e})")
    return tables


def _saddle_restarts(work) -> dict:
    """{(molecule, medium): restarts} from the history the opt stage wrote."""
    path = Path(work) / "saddle_history.csv"
    if not path.exists():
        return {}
    try:
        h = pd.read_csv(path)
    except Exception:
        return {}
    out = {}
    for (name, medium), g in h.groupby(["molecule", "medium"]):
        if bool(g["resolved"].iloc[-1]):
            out[(name, medium)] = int(g["attempt"].max())
    return out


def _flat_modes(work) -> dict:
    """{(molecule, medium): True} where the leftover imaginary mode is an internal rotation."""
    path = Path(work) / "saddle_history.csv"
    if not path.exists():
        return {}
    try:
        h = pd.read_csv(path)
    except Exception:
        return {}
    if "flat_mode" not in h.columns:
        return {}
    out = {}
    for (name, medium), g in h.groupby(["molecule", "medium"]):
        out[(name, medium)] = bool(g["flat_mode"].iloc[-1])
    return out
