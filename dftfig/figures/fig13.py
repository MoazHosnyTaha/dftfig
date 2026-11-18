"""Figure 13 - molecular lipophilicity potential surfaces with TPSA."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.gridspec import GridSpec

from .. import mlp as mlpmod
from .. import render
from ..figstyle import (apply_style, cell_caption,
                        pad_to_common_aspect, show_image, trim_white)
from ..render import MLP_PALETTES
from ..theme import warn
from .fig2 import PYMOL_RGB, _save


def _mlp_colormap(palette: str) -> LinearSegmentedColormap:
    cols = MLP_PALETTES[palette]
    return LinearSegmentedColormap.from_list(f"mlp_{palette}", [PYMOL_RGB[c] for c in cols])


def build(cfg: dict, mols: dict, prop_df: pd.DataFrame) -> list[Path]:
    apply_style(int(cfg["figures"]["dpi"]))
    names = list(mols)

    prepared = {}
    lo_all, hi_all = [], []
    for n in names:
        mol = mols.get(n)
        if mol is None:
            continue
        pdb, lo, hi = mlpmod.prepare_mlp(cfg, n, mol)
        if pdb is None:
            continue
        prepared[n] = pdb
        lo_all.append(lo); hi_all.append(hi)

    if not prepared:
        warn("no lipophilicity surfaces could be prepared - skipping Figure 13")
        return []

    # One shared scale across every panel, so colours mean the same thing throughout.
    lo = float(np.min(lo_all)); hi = float(np.max(hi_all))
    bound = max(abs(lo), abs(hi))
    lo, hi = -bound, bound

    images = {}
    for n, pdb in prepared.items():
        p = render.render_mlp(cfg, n, pdb, lo, hi)
        if p and Path(p).exists() and cfg["render"]["background"] == "transparent":
            trim_white(p)
        images[n] = p

    pad_to_common_aspect([images.get(n) for n in names])

    ncol = len(names)
    fig = plt.figure(figsize=(2.1 * ncol + 0.5, 3.5))
    gs = GridSpec(2, ncol, figure=fig, height_ratios=[1.0, 0.075], hspace=0.30,
                  wspace=0.03, left=0.03, right=0.985, top=0.86, bottom=0.14)

    letters = "ABCDEFGH"
    for j, n in enumerate(names):
        ax = fig.add_subplot(gs[0, j])
        row = prop_df[prop_df["molecule"] == n]
        tpsa = float(row["TPSA"].iloc[0]) if len(row) and not pd.isna(row["TPSA"].iloc[0]) else np.nan
        clogp = row["clogP_rdkit"].iloc[0] if len(row) and "clogP_rdkit" in row else np.nan
        sub = []
        if not pd.isna(tpsa):
            sub.append(f"TPSA {tpsa:.2f} $\\AA^2$")
        if not pd.isna(clogp):
            sub.append(f"clogP {float(clogp):.2f}")
        # Title and caption go on the grid cell, not the axes. A panel whose surface
        # failed to render has no image, so its axes keeps the full cell height while its
        # neighbours shrink to their image aspect - and anything anchored to the axes then
        # sits at a different height from the rest of the row.
        show_image(ax, images.get(n))
        cell_caption(fig, ax, f"{letters[j]}  {n}", size=9, color="#2A2A2A", dy=0.965,
                     va="bottom")
        cell_caption(fig, ax, "    ".join(sub), size=7.5)

    cax = fig.add_subplot(gs[1, :])
    cmap = _mlp_colormap(cfg["figures"]["mlp_palette"])
    cax.imshow(np.linspace(0, 1, 512)[None, :], aspect="auto", cmap=cmap,
               extent=(lo, hi, 0, 1))
    cax.set_yticks([])
    cax.set_xticks(np.linspace(lo, hi, 5))
    cax.set_xticklabels([f"{t:+.2f}" for t in np.linspace(lo, hi, 5)])
    cax.set_xlabel("per-atom Crippen lipophilicity contribution (log$P$ units)",
                   fontsize=7, labelpad=1)
    for s in ("top", "right", "left"):
        cax.spines[s].set_visible(False)
    left_lab, right_lab = ("least lipophilic", "most lipophilic")
    cax.text(lo, 1.35, left_lab, fontsize=6.5, ha="left", va="bottom")
    cax.text(hi, 1.35, right_lab, fontsize=6.5, ha="right", va="bottom")

    fig.suptitle("Molecular lipophilicity potential mapped on the molecular surface",
                 fontsize=9.5, x=0.03, y=0.995, ha="left")
    return _save(cfg, fig, "figure13")
