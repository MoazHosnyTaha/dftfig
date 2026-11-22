"""CNS MPO physicochemical space and site-of-metabolism maps.

Reference implementation of a figure that reads a per-atom table - the condensed Fukui
indices - and paints it onto a 2D depiction. The context locates the table; folding the
hydrogens onto their heavy atoms is this figure's own chemistry and stays here.
"""
from __future__ import annotations

import io

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from rdkit import Chem
from rdkit.Chem.Draw import SimilarityMaps, rdMolDraw2D

from ..figstyle import caption, panel_letter
from ..figureapi import REQ_PROPERTIES, FigureContext
from ..theme import warn

NAME = "mpo_metabolism"
TITLE = "CNS MPO space and predicted sites of metabolism"
REQUIRES = (REQ_PROPERTIES,)

HUES = ["#7C4DFF", "#00A5A5", "#FF8A3D", "#C2185B", "#2E7D32", "#5D4037"]


def _som_image(mol: Chem.Mol, weights: np.ndarray, size=(380, 320)):
    m = Chem.RemoveHs(Chem.Mol(mol))
    Chem.rdDepictor.Compute2DCoords(m)
    d = rdMolDraw2D.MolDraw2DCairo(*size)
    d.drawOptions().baseFontSize = 0.55
    SimilarityMaps.GetSimilarityMapFromWeights(m, list(weights), draw2d=d)
    d.FinishDrawing()
    import PIL.Image
    return np.array(PIL.Image.open(io.BytesIO(d.GetDrawingText())))


def _fukui_radical(ctx: FigureContext, name: str, mol: Chem.Mol):
    """Per-heavy-atom radical Fukui index f0, mapped onto the heavy-atom depiction order."""
    t = ctx.fukui(name)
    if t is None:
        return None, "run the descriptors stage first"
    if "f_zero" not in t.columns:
        return None, "f_zero column missing"
    # The contract guarantees the Fukui table is indexed by the same atom order as
    # Compound.mol, hydrogens included.
    full = np.zeros(mol.GetNumAtoms())
    for _, r in t.iterrows():
        i = int(r["atom_index_0based"])
        if 0 <= i < len(full):
            full[i] = float(r["f_zero"])
    # Fold each hydrogen's value onto the heavy atom it is bonded to, then drop hydrogens:
    # CYP abstraction at a C-H is reported at the carbon.
    heavy_ids = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    folded = {hid: full[hid] for hid in heavy_ids}
    for a in mol.GetAtoms():
        if a.GetAtomicNum() == 1:
            for nb in a.GetNeighbors():
                if nb.GetIdx() in folded:
                    folded[nb.GetIdx()] += full[a.GetIdx()]
    return np.array([folded[i] for i in heavy_ids]), None


def _external_som(ctx: FigureContext, name: str, n_heavy: int):
    df = ctx.som_table(name)
    if df is None:
        return None, f"{name}.csv not found"
    cols = {c.lower().strip(): c for c in df.columns}
    icol = cols.get("atom_index") or cols.get("atom") or list(df.columns)[0]
    pcol = cols.get("probability") or cols.get("prob") or cols.get("score") or list(df.columns)[1]
    w = np.zeros(n_heavy)
    bad = 0
    for _, r in df.iterrows():
        try:
            i = int(r[icol]) - 1
        except (ValueError, TypeError):
            bad += 1
            continue
        if 0 <= i < n_heavy:
            w[i] = float(r[pcol])
        else:
            bad += 1
    if bad:
        warn(f"{name}: {bad} SOM row(s) fall outside the molecule's "
             f"{n_heavy} heavy atoms and were ignored")
    return w, None


def build(ctx: FigureContext) -> list:
    ctx.apply_style()
    names = ctx.names
    nmol = len(names)

    fig = plt.figure(figsize=(max(7.2, 1.85 * nmol + 0.4), 7.0))
    gs = GridSpec(2, nmol, figure=fig, height_ratios=[1.35, 1.0], hspace=0.30,
                  wspace=0.04, left=0.095, right=0.985, top=0.935, bottom=0.075)

    # ---------------- Panel A: logD vs TPSA, CNS MPO desirability space ---------------
    axa = fig.add_subplot(gs[0, :])
    # Wager's desirable region: TPSA 40-90 A^2 and clogD <= 2 fully desirable, <= 4 tolerated.
    axa.axhspan(40, 90, xmin=0, xmax=1, color="#7C4DFF", alpha=0.07, zorder=0)
    axa.axvline(2.0, color="#00A5A5", lw=0.9, ls="--", zorder=1)
    axa.axvline(4.0, color="#C62828", lw=0.9, ls="--", zorder=1)
    axa.annotate("clogD $\\leq$ 2\nfully desirable", xy=(2.0, 0.015),
                 xycoords=("data", "axes fraction"), xytext=(-4, 4),
                 textcoords="offset points", fontsize=6.5, color="#00A5A5",
                 ha="right", va="bottom")
    axa.annotate("clogD > 4\nmetabolically unstable", xy=(4.0, 0.015),
                 xycoords=("data", "axes fraction"), xytext=(4, 4),
                 textcoords="offset points", fontsize=6.5, color="#C62828",
                 ha="left", va="bottom")

    xs, ys, ss, cs, labels, mpos = [], [], [], [], [], []
    missing_logd, substituted, derived = [], [], []
    for n in names:
        if ctx.property_row(n) is None:
            continue
        logd = ctx.property(n, "logD74", np.nan)
        tpsa = ctx.property(n, "TPSA", np.nan)
        if pd.isna(logd):
            # No measured or predicted logD: fall back to the calculated neutral logP
            # so the panel still carries information. These are NOT the same quantity -
            # for an ionizable or permanently charged compound logD(7.4) is lower, often
            # by several log units - so the substitution is labelled on the figure.
            clogp = ctx.property(n, "clogP_rdkit", np.nan)
            if not pd.isna(clogp):
                logd = clogp
                substituted.append(n)
        if pd.isna(logd) or pd.isna(tpsa):
            missing_logd.append(n)
            continue
        if str(ctx.property(n, "logD74_source", "")) == "derived":
            derived.append(n)
        papp = ctx.property(n, "Papp_MDCK", np.nan)
        size = 90.0 if pd.isna(papp) else 60.0 + 900.0 * float(papp) / 1e-4
        xs.append(float(logd)); ys.append(float(tpsa)); ss.append(size)
        mpo = float(ctx.property(n, "cns_mpo", np.nan))
        cs.append(mpo); labels.append(ctx.label(n)); mpos.append(mpo)

    if xs:
        sc = axa.scatter(xs, ys, s=ss, c=cs, cmap="viridis", vmin=0, vmax=6,
                         edgecolor="white", linewidth=0.9, zorder=3)
        cb = fig.colorbar(sc, ax=axa, pad=0.012, fraction=0.035)
        cb.set_label("CNS MPO score (0\u20136)", fontsize=7.5)
        cb.ax.tick_params(labelsize=7)
        # Positional isomers can share a point exactly (CCB2 and CCB3 have identical TPSA
        # and near-identical logD). Merge coincident points into one label rather than
        # printing two strings on top of each other.
        xr = (max(xs) - min(xs)) or 1.0
        yr = (max(ys) - min(ys)) or 1.0
        clusters: list[dict] = []
        for x, y, n, m in zip(xs, ys, labels, mpos):
            for c in clusters:
                if (abs(x - c["x"]) / xr < 0.04) and (abs(y - c["y"]) / yr < 0.04):
                    c["members"].append((n, m))
                    break
            else:
                clusters.append({"x": x, "y": y, "members": [(n, m)]})

        # Two compounds at the same TPSA but different logD are not merged into one
        # label, and their texts still overlap if both are offset the same way.
        #
        # The vertical side is decided by where the point sits, not by its index: a
        # label pushed DOWN from a point near the bottom of the axes lands on the
        # "clogD <= 2 / fully desirable" threshold annotation, which is pinned at
        # y = 0.015 in axes fraction. So low points are labelled above and high points
        # below, which also keeps every label inside the axes. Neighbours are then
        # separated horizontally instead, alternating along the x order.
        clusters.sort(key=lambda c: c["x"])
        axa.margins(0.22)                      # before the transform is queried
        for k, c in enumerate(clusters):
            lines = [f"{n}  MPO {m:.2f}" if not np.isnan(m) else n for n, m in c["members"]]
            y_frac = axa.transAxes.inverted().transform(
                axa.transData.transform((c["x"], c["y"])))[1]
            up = y_frac < 0.55
            dyk = 11 if up else -11 - 7 * (len(lines) - 1)
            dxk, ha = (9, "left") if k % 2 == 0 else (-9, "right")
            axa.annotate("\n".join(lines), xy=(c["x"], c["y"]), xytext=(dxk, dyk),
                         textcoords="offset points", fontsize=7.5, color="#2A2A2A",
                         ha=ha, va="bottom" if up else "top",
                         arrowprops={"arrowstyle": "-", "lw": 0.6, "color": "#B0B0B0"})
        if derived and not substituted:
            axa.annotate("logD derived from Crippen clogP and the predicted ionization "
                         "(neutral species only); a lower bound for strong bases",
                         xy=(0.0, -0.155), xycoords="axes fraction", fontsize=6.3,
                         color="#2F6FA8", ha="left", va="top")
        if substituted:
            axa.annotate("clogP substituted for logD (no logD supplied): "
                         + ", ".join(substituted)
                         + " - the real logD is lower for anything that ionizes",
                         xy=(0.0, -0.155), xycoords="axes fraction", fontsize=6.3,
                         color="#C77800", ha="left", va="top")
    else:
        axa.text(0.5, 0.5, "logD and TPSA are required for this panel\n"
                           "(supply logD74 in the properties CSV)",
                 transform=axa.transAxes, ha="center", va="center",
                 fontsize=9, color="#C62828")

    if missing_logd:
        axa.annotate("not plotted (no logD): " + ", ".join(missing_logd),
                     xy=(0.0, -0.13), xycoords="axes fraction", fontsize=6.5,
                     color="#9E9E9E", ha="left")

    axa.set_xlabel("buffer lipophilicity, clogD at pH 7.4"
                   + ("   (clogP substituted where logD was not supplied)" if substituted else ""))
    axa.set_ylabel("topological polar surface area / $\\AA^2$")
    axa.set_title("Alignment of lipophilicity and polarity against the CNS MPO desirable space",
                  fontsize=9)
    # Inside the axes, top right: the strip under the axes carries the provenance note.
    axa.annotate("bubble area $\\propto$ $P_{app,\\,A\\rightarrow B}$ (MDCK)",
                 xy=(0.988, 0.975), xycoords="axes fraction", fontsize=6.5,
                 ha="right", va="top", color="#707070")
    axa.annotate("TPSA 40\u201390 $\\AA^2$", xy=(0.012, 65), xycoords=("axes fraction", "data"),
                 fontsize=6.5, color="#7C4DFF", va="center")

    # ---------------- Panels B-E: site-of-metabolism sensitivity maps -----------------
    letters = "BCDEFGHI"
    source = ctx.param("figures.som_source")
    for j, n in enumerate(names):
        ax = fig.add_subplot(gs[1, j])
        ax.set_axis_off()
        compound = ctx.compound(n)
        if compound is None:
            continue
        mol = compound.mol
        heavy = Chem.RemoveHs(Chem.Mol(mol))
        w, err, used = None, None, ""
        if source in ("csv", "both"):
            w, err = _external_som(ctx, n, heavy.GetNumAtoms())
            used = "XenoSite/ADMETlab export"
        if w is None and source in ("fukui_radical", "both"):
            w, err = _fukui_radical(ctx, n, mol)
            used = "radical Fukui index $f^{0}$"
        ax.set_title(f"{letters[j]}  {ctx.label(n)}", fontsize=9)
        if w is None:
            ax.text(0.5, 0.5, f"no metabolism weights\n\n{err}",
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=7, color="#9E9E9E")
            continue
        ax.imshow(_som_image(mol, w))
        k = int(np.argmax(w))
        unit = "p" if used.startswith("XenoSite") else "$f^{0}$"
        caption(ax, f"most labile: {heavy.GetAtomWithIdx(k).GetSymbol()}{k + 1}"
                    f"  ({unit} = {w.max():.3f})", size=7)
        if j == 0:
            ax.annotate(f"weights from {used}", xy=(0.0, -0.22), xycoords="axes fraction",
                        fontsize=6.5, color="#707070", ha="left")

    panel_letter(axa, "A", ctx.param("figures.panel_label_style"), dx=-0.075, dy=1.02)
    fig.suptitle("Multiparameter optimization and predicted sites of metabolism",
                 fontsize=9.5, x=0.095, y=0.995, ha="left")
    return ctx.save(fig, "mpo_metabolism")
