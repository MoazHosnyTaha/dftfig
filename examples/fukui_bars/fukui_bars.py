"""Condensed Fukui indices as a bar panel, one column per compound.

Drop this file into ./dft_figures in a dftfig project and add `fukui_bars` to
figures.which. Nothing else is needed - no packaging, no edit to dftfig.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from dftfig.figureapi import REQ_FUKUI, FigureContext
from dftfig.figstyle import panel_letter

NAME = "fukui_bars"
TITLE = "Condensed Fukui indices per heavy atom"
REQUIRES = (REQ_FUKUI,)


def build(ctx: FigureContext) -> list:
    ctx.apply_style()
    tables = {n: ctx.fukui(n) for n in ctx.names}
    tables = {n: t for n, t in tables.items() if t is not None}
    if not tables:
        return []                      # nothing computed yet; the stage says so already

    # One y scale for every panel: the point of the figure is that one compound is more
    # reactive than another, and per-panel scaling would hide exactly that.
    top = max(float(t[["f_plus", "f_minus"]].to_numpy().max()) for t in tables.values())

    names = list(tables)
    fig, axes = plt.subplots(1, len(names), figsize=(2.6 * len(names) + 0.8, 3.2),
                             squeeze=False, sharey=True)
    for k, (ax, name) in enumerate(zip(axes[0], names)):
        t = tables[name][tables[name]["element"] != "H"]
        x = np.arange(len(t))
        ax.bar(x - 0.2, t["f_plus"], width=0.4, color="#7C4DFF",
               label="$f^{+}$  (nucleophilic attack here)")
        ax.bar(x + 0.2, t["f_minus"], width=0.4, color="#FF8A3D",
               label="$f^{-}$  (electrophilic attack here)")
        ax.set_xticks(x)
        ax.set_xticklabels([f"{e}{i}" for e, i in zip(t["element"], t["atom_index_1based"])],
                           rotation=90, fontsize=6)
        ax.set_ylim(0, top * 1.15)
        # The panel letter is drawn just above the axes, so the title needs to
        # sit below it rather than share the line.
        ax.set_title(f"{name}   {ctx.compound(name).formula}", fontsize=9, pad=3)
        ax.spines[["top", "right"]].set_visible(False)
        panel_letter(ax, "ABCDEFGH"[k], ctx.param("figures.panel_label_style"),
                     dy=1.13)
        if k == 0:
            ax.set_ylabel("condensed Fukui index")
            ax.legend(loc="upper left", fontsize=6.5, framealpha=0.9)

    fig.suptitle(f"Site reactivity in {ctx.reference_solvent} "
                 f"({ctx.param('property.population')} populations)",
                 fontsize=9.5, x=0.02, ha="left")
    fig.text(0.02, 0.005, f"{ctx.level_of_theory()}   |   dftfig "
                          f"{ctx.provenance['dftfig_version']}   |   "
                          f"config {ctx.provenance['config_sha']}",
             fontsize=6.5, color="#707070", ha="left", va="bottom")
    fig.tight_layout(rect=(0, 0.03, 1, 0.94))
    return ctx.save(fig, "fukui_bars")
