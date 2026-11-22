"""Molecular lipophilicity potential surfaces with TPSA.

Reference implementation of a figure that needs one colour scale shared across panels:
the range is collected from every compound first, then handed back when each image is
requested.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec

from ..figstyle import cell_caption, pad_to_common_aspect, show_image
from ..figureapi import REQ_RENDERING, FigureContext
from ..theme import warn

NAME = "lipophilicity"
TITLE = "Molecular lipophilicity potential on the molecular surface"
REQUIRES = (REQ_RENDERING,)


def build(ctx: FigureContext) -> list:
    ctx.apply_style()
    names = ctx.names

    ranges = {}
    for n in names:
        r = ctx.lipophilicity_range(n)
        if r is not None:
            ranges[n] = r

    if not ranges:
        warn("no lipophilicity surfaces could be prepared - skipping Figure 13")
        return []

    # One shared scale across every panel, so colours mean the same thing throughout.
    bound = max(abs(min(r[0] for r in ranges.values())),
                abs(max(r[1] for r in ranges.values())))
    lo, hi = -bound, bound

    images = {n: ctx.image(n, "lipophilicity", vrange=(lo, hi)) for n in ranges}
    pad_to_common_aspect([images.get(n) for n in names])

    ncol = len(names)
    fig = plt.figure(figsize=(2.1 * ncol + 0.5, 3.5))
    gs = GridSpec(2, ncol, figure=fig, height_ratios=[1.0, 0.075], hspace=0.30,
                  wspace=0.03, left=0.03, right=0.985, top=0.86, bottom=0.14)

    letters = "ABCDEFGH"
    for j, n in enumerate(names):
        ax = fig.add_subplot(gs[0, j])
        tpsa = ctx.property(n, "TPSA", np.nan)
        clogp = ctx.property(n, "clogP_rdkit", np.nan)
        sub = []
        if not pd.isna(tpsa):
            sub.append(f"TPSA {float(tpsa):.2f} $\\AA^2$")
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
    cax.imshow(np.linspace(0, 1, 512)[None, :], aspect="auto", cmap=ctx.mlp_colormap(),
               extent=(lo, hi, 0, 1))
    cax.set_yticks([])
    cax.set_xticks(np.linspace(lo, hi, 5))
    cax.set_xticklabels([f"{t:+.2f}" for t in np.linspace(lo, hi, 5)])
    cax.set_xlabel("per-atom Crippen lipophilicity contribution (log$P$ units)",
                   fontsize=7, labelpad=1)
    for s in ("top", "right", "left"):
        cax.spines[s].set_visible(False)
    cax.text(lo, 1.35, "least lipophilic", fontsize=6.5, ha="left", va="bottom")
    cax.text(hi, 1.35, "most lipophilic", fontsize=6.5, ha="right", va="bottom")

    fig.suptitle("Molecular lipophilicity potential mapped on the molecular surface",
                 fontsize=9.5, x=0.03, y=0.995, ha="left")
    return ctx.save(fig, "lipophilicity")
