"""Figure 12 - acid-base analysis: pKa annotation and net charge versus pH."""
from __future__ import annotations

import io
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D

from .. import properties as props
from ..figstyle import apply_style, caption, panel_letter
from ..theme import ok, warn
from .fig2 import _save

HUES = ["#7C4DFF", "#00A5A5", "#FF8A3D", "#C2185B", "#2E7D32", "#5D4037"]


def _annotated_2d(mol: Chem.Mol, acidic: list[float], basic: list[float],
                  size=(420, 330)) -> np.ndarray:
    """2D depiction with hydrogen-bond donors/acceptors marked A/D and pKa values shown."""
    m = Chem.RemoveHs(Chem.Mol(mol))
    Chem.rdDepictor.Compute2DCoords(m)

    # Atomic-number primitives (#7/#8/#16), not N/O/S: the uppercase element symbols in
    # SMARTS mean ALIPHATIC, so [N;!H0] silently misses an aromatic carbazole N-H.
    donor = Chem.MolFromSmarts("[$([#7;!H0]),$([#8,#16;H1])]")
    acceptor = Chem.MolFromSmarts("[$([#8,#16;H0;v2]),$([#8,#16;-]),$([#7;v3;H0])]")
    d_idx = {i for (i,) in m.GetSubstructMatches(donor)}
    a_idx = {i for (i,) in m.GetSubstructMatches(acceptor)}

    for atom in m.GetAtoms():
        i = atom.GetIdx()
        tag = ""
        if i in d_idx:
            tag = "D"
        elif i in a_idx:
            tag = "A"
        if tag:
            atom.SetProp("atomNote", tag)

    # Attach each acidic pKa to the most acidic donor (phenol OH before the ring NH).
    oh = [i for i in d_idx if m.GetAtomWithIdx(i).GetSymbol() == "O"]
    nh = [i for i in d_idx if m.GetAtomWithIdx(i).GetSymbol() == "N"]
    for pka, i in zip(sorted(acidic), oh + nh):
        a = m.GetAtomWithIdx(i)
        a.SetProp("atomNote", f"D p$K_a$ {pka:.1f}" if False else f"D {pka:.1f}")

    drawer = rdMolDraw2D.MolDraw2DCairo(*size)
    o = drawer.drawOptions()
    o.addAtomIndices = False
    o.annotationFontScale = 0.70
    o.baseFontSize = 0.55
    # Without this RDKit scales every molecule to fill its canvas, so benzene comes out
    # three times the size of a substituted ring beside it and the row reads as though
    # the compounds differ in size. A fixed bond length gives one scale across the row;
    # RDKit still shrinks anything that would not otherwise fit.
    o.fixedBondLength = 46.0
    highlight = sorted(d_idx | a_idx)
    cols = {}
    for i in highlight:
        cols[i] = (0.85, 0.75, 1.0) if i in d_idx else (0.75, 0.93, 0.93)
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, m, highlightAtoms=highlight,
                                       highlightAtomColors=cols)
    drawer.FinishDrawing()
    import PIL.Image
    return np.array(PIL.Image.open(io.BytesIO(drawer.GetDrawingText())))


def build(cfg: dict, mols: dict, prop_df: pd.DataFrame,
          charges: dict | None = None) -> list[Path]:
    apply_style(int(cfg["figures"]["dpi"]))
    names = list(mols)
    ph_lo, ph_hi = [float(x) for x in str(cfg["figures"]["ph_range"]).split(",")]
    ph_phys = float(cfg["figures"]["physiological_ph"])
    ph = np.linspace(ph_lo, ph_hi, 600)

    charges = charges or {}
    # Does any compound actually have a pKa? That decides the layout, because the
    # explanatory note needs its own strip of figure under the axes - drawn inside
    # them it would sit on top of a curve and hide it.
    has_pka = False
    for n in names:
        row = prop_df[prop_df["molecule"] == n]
        if len(row) and (props.parse_pka_list(row["pKa_acidic"].iloc[0])
                         or props.parse_pka_list(row["pKa_basic"].iloc[0])):
            has_pka = True
            break

    ncol = len(names)
    fig = plt.figure(figsize=(2.15 * ncol + 0.6, 6.4 if has_pka else 6.9))
    gs = GridSpec(2, ncol, figure=fig, height_ratios=[1.0, 1.25], hspace=0.30,
                  wspace=0.04, left=0.08, right=0.985, top=0.93,
                  bottom=0.11 if has_pka else 0.20)

    any_pka = False
    predicted_names: set = set()
    for j, n in enumerate(names):
        ax = fig.add_subplot(gs[0, j])
        ax.set_axis_off()
        row = prop_df[prop_df["molecule"] == n]
        acidic = props.parse_pka_list(row["pKa_acidic"].iloc[0]) if len(row) else []
        basic = props.parse_pka_list(row["pKa_basic"].iloc[0]) if len(row) else []
        any_pka = any_pka or bool(acidic or basic)
        mol = mols.get(n)
        if mol is not None:
            ax.imshow(_annotated_2d(mol, acidic, basic))
        q = int(charges.get(n, 0))
        src = str(row["pKa_source"].iloc[0]) if len(row) and "pKa_source" in row else ""
        perm = int(row["permanent_charge"].iloc[0] or 0) if (
            len(row) and "permanent_charge" in row) else 0
        txt = []
        if acidic:
            txt.append("acidic p$K_a$ " + ", ".join(f"{v:.1f}" for v in acidic))
        if basic:
            txt.append("basic p$K_a$ " + ", ".join(f"{v:.1f}" for v in basic))
        if perm:
            txt.append(f"permanent charge {perm:+d}")
        # The input's formal charge is a drawing convention, not a result: a pyrrolidine
        # drawn as its hydrochloride is the same compound as one drawn neutral, and panel B
        # works out the charge from pH. Say which form was supplied and leave it there.
        # Second line, not a suffix: a long compound name plus a parenthetical runs
        # straight into the neighbouring panel's title.
        if q and src in ("predicted", "csv"):
            sub = f"supplied as the {q:+d} ion"
        elif q:
            sub = f"formal charge {q:+d}"
        else:
            sub = ""
        ax.set_title(n + (f"\n{sub}" if sub else ""), fontsize=9,
                     linespacing=1.35)
        if txt:
            colour = "#4A4A4A" if src == "csv" else "#2F6FA8"
            line = "   ".join(txt) + ("   (predicted)" if src == "predicted" else "")
            caption(ax, line, size=7, color=colour)
            predicted_names.add(n) if src == "predicted" else None
        elif src == "predicted":
            # A real, checkable statement now: the model looked and found nothing that
            # ionizes between pH 0 and 14. Carbazole and its ethers are genuinely in
            # this class - the ring N-H is around pKa 20.
            caption(ax, "no ionizable group between pH 0 and 14", size=7, color="#2E7D32")
            predicted_names.add(n)
        else:
            caption(ax, "no p$K_a$ supplied - shown at its formal charge", size=7,
                    color="#C77800")

    axb = fig.add_subplot(gs[1, :])
    ends = []
    qmin, qmax = 0.0, 0.0
    dashed_any = False
    for i, n in enumerate(names):
        row = prop_df[prop_df["molecule"] == n]
        acidic = props.parse_pka_list(row["pKa_acidic"].iloc[0]) if len(row) else []
        basic = props.parse_pka_list(row["pKa_basic"].iloc[0]) if len(row) else []
        # Baseline: the charge no pH can remove. For a compound drawn in its protonated
        # form that is zero - the +1 is already described by the basic pKa, and adding the
        # input's formal charge on top would double-count it.
        rowp = prop_df[prop_df["molecule"] == n]
        src_n = str(rowp["pKa_source"].iloc[0]) if (
            len(rowp) and "pKa_source" in rowp) else ""
        if src_n in ("predicted", "csv") and "permanent_charge" in rowp:
            base_q = int(rowp["permanent_charge"].iloc[0] or 0)
        else:
            base_q = int(charges.get(n, 0))
        qcurve = props.net_charge(ph, acidic, basic, base_q)
        style = "-" if (acidic or basic or base_q) else "--"
        dashed_any = dashed_any or style == "--"
        axb.plot(ph, qcurve, color=HUES[i % len(HUES)], label=n, lw=1.6, ls=style)
        ends.append([float(qcurve[-1]), n, HUES[i % len(HUES)]])
        qmin, qmax = min(qmin, float(qcurve.min())), max(qmax, float(qcurve.max()))

    axb.axvline(ph_phys, color="#B0B0B0", lw=0.9, ls="--", zorder=0)
    axb.annotate(f"pH {ph_phys:g}", xy=(ph_phys, 0.015), xycoords=("data", "axes fraction"),
                 xytext=(3, 0), textcoords="offset points", fontsize=7,
                 color="#707070", va="bottom")
    axb.axhline(0, color="#D8D8D8", lw=0.8, zorder=0)
    axb.set_xlabel("pH")
    axb.set_ylabel("net molecular charge")
    axb.set_xlim(ph_lo, ph_hi + (ph_hi - ph_lo) * 0.07)
    # Explicit limits: a permanently charged compound sits at +-1 or beyond, and
    # autoscaling alone has repeatedly left such a curve under an annotation.
    qlo = min(qmin, min(e[0] for e in ends))
    qhi = max(qmax, max(e[0] for e in ends))
    pad = max(0.18, 0.12 * (qhi - qlo or 1.0))
    axb.set_ylim(qlo - pad, qhi + pad)

    # Stagger end-of-line labels that would otherwise print on top of each other where
    # several curves converge on the same final charge. Separation is measured against the
    # axis, not against the spread of the end values, which is often zero.
    span = (qhi + pad) - (qlo - pad)
    ends.sort(key=lambda e: e[0])
    for k in range(1, len(ends)):
        if ends[k][0] - ends[k - 1][0] < 0.07 * span:
            ends[k][0] = ends[k - 1][0] + 0.07 * span
    for y, n, colour in ends:
        axb.annotate(n, xy=(ph[-1], y), xytext=(4, 0), textcoords="offset points",
                     fontsize=8, color=colour, va="center")
    axb.set_title("Net charge as a function of pH (independent-site Henderson\u2013Hasselbalch)",
                  fontsize=9)
    if not any_pka:
        # Below the axes, never over them.
        fig.text(0.5, 0.025,
                 "No p$K_a$ values supplied, so no ionization is modelled: each compound is "
                 "drawn flat (dashed)\nat the formal charge of the structure you provided. "
                 "Add pKa_acidic / pKa_basic to\n" + str(cfg["figures"]["properties_csv"])
                 + " to get real speciation curves.",
                 ha="center", va="bottom", fontsize=8, color="#C77800",
                 bbox=dict(boxstyle="round,pad=0.4", fc="#FFF6E5", ec="#C77800", lw=0.7))
    elif dashed_any:
        axb.plot([], [], color="#9E9E9E", ls="--", lw=1.4,
                 label="dashed = nothing ionizes in this range; flat at the formal charge")
        axb.legend(loc="best", fontsize=7, framealpha=0.85)

    if predicted_names:
        fig.text(0.5, 0.012,
                 "p$K_a$ predicted from structure (SMARTS class averages, about "
                 "$\\pm$1 log unit) - supply measured values in the properties CSV to "
                 "override", ha="center", va="bottom", fontsize=7, color="#2F6FA8")

    panel_letter(fig.axes[0], "A", cfg["figures"]["panel_label_style"], dx=-0.10, dy=1.02)
    panel_letter(axb, "B", cfg["figures"]["panel_label_style"], dx=-0.055, dy=1.02)
    fig.suptitle("Acid\u2013base behaviour of the derivatives", fontsize=9.5,
                 x=0.08, y=0.995, ha="left")
    return _save(cfg, fig, "figure12")
