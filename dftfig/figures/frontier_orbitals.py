"""Optimized geometries, frontier orbitals with the gap, and MEP surfaces.

Reference implementation of the plugin contract for the hardest case: three rows of 3D
panels, a shared colour scale measured from the data, and a provenance stamp. Everything
it reads comes from the :class:`~dftfig.figureapi.FigureContext`.
"""
from __future__ import annotations

import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec

from ..figstyle import band_letter, cell_caption, pad_to_common_aspect, show_image
from ..figureapi import REQ_CUBES, REQ_DESCRIPTORS, REQ_RENDERING, FigureContext
from ..theme import warn

NAME = "frontier_orbitals"
TITLE = "Optimized geometry, frontier orbitals and electrostatic potential"
REQUIRES = (REQ_RENDERING, REQ_CUBES, REQ_DESCRIPTORS)


def build(ctx: FigureContext) -> list:
    """Render every panel, then assemble the figure."""
    ctx.apply_style()
    names = ctx.names
    solvent = ctx.reference_solvent

    # --- MEP colour range -------------------------------------------------------------
    # The potential around an ion is dominated by its monopole: a +1 cation is positive
    # essentially everywhere and far outside the +-0.05 a.u. that suits neutral organics,
    # so a fixed range renders it as a saturated, featureless blob. Measure the surface
    # potential and scale accordingly.
    mode = ctx.param("render.esp_scale", "auto_shared")
    if ctx.param("render.esp_autoscale"):          # legacy flag still honoured
        mode = "per_panel"
    r0 = float(ctx.param("render.esp_range"))
    stats = {n: ctx.surface_esp_stats(n) for n in names}
    measured = {n: (s["v_lo"], s["v_hi"]) for n, s in stats.items() if s}

    def _round(lo, hi):
        span = hi - lo
        q = 10 ** (math.floor(math.log10(span)) - 1) if span > 0 else 0.001
        return math.floor(lo / q) * q, math.ceil(hi / q) * q

    vrange = (-r0, r0)
    per_panel_range: dict[str, tuple[float, float]] = {}
    if mode == "per_panel" and measured:
        per_panel_range = {n: _round(*v) for n, v in measured.items()}
    elif mode == "auto_shared" and measured:
        vrange = _round(min(v[0] for v in measured.values()),
                        max(v[1] for v in measured.values()))
    elif mode == "fixed" and measured:
        over = [n for n, v in measured.items()
                if v[0] < -r0 * 1.3 or v[1] > r0 * 1.3]
        if over:
            warn(f"the MEP surface of {', '.join(over)} runs well beyond the fixed range "
                 f"+-{r0} a.u., so those panels saturate and lose all detail. This is what "
                 f"happens to charged species - a +1 cation is positive everywhere. "
                 f"Set render.esp_scale: auto_shared.")

    # The context writes the geometry the renderer reads and trims the white border when
    # the background is transparent - both used to be repeated here by hand.
    panels: dict[str, dict] = {}
    for n in names:
        panels[n] = {
            "geom": ctx.image(n, "geometry"),
            "homo": ctx.image(n, "homo"),
            "lumo": ctx.image(n, "lumo"),
            "mep": ctx.image(n, "mep", vrange=per_panel_range.get(n, vrange)),
        }

    # Equalize panel aspect within each row so the labels align across columns.
    for role in ("geom", "homo", "lumo", "mep"):
        pad_to_common_aspect([panels[n][role] for n in names])

    # --- energies ---------------------------------------------------------------------
    en: dict[str, dict] = {}
    for n in names:
        r = ctx.descriptor_row(n, solvent)
        if r is None:
            continue
        en[n] = {
            "homo": r.get("e_homo_eV"), "lumo": r.get("e_lumo_eV"),
            "gap": r.get("gap_eV"),
            "n_imag": r.get("n_imaginary"),
            "is_min": r.get("is_minimum"),
            "restarts": int(r.get("saddle_restarts") or 0),
            "flat_mode": bool(r.get("flat_mode")),
            "lowest_freq": r.get("lowest_freq_cm-1"),
            "lumo_unbound": (r.get("e_lumo_eV") is not None
                             and not pd.isna(r.get("e_lumo_eV"))
                             and r.get("e_lumo_eV") > 0),
        }

    ncol = len(names)
    fig = plt.figure(figsize=(2.0 * ncol + 0.9, 8.2))
    gs = GridSpec(5, ncol, figure=fig,
                  height_ratios=[1.05, 1.0, 1.0, 1.15, 0.085],
                  hspace=0.16, wspace=0.03,
                  left=0.055, right=0.985, top=0.945, bottom=0.075)
    label_style = ctx.param("figures.panel_label_style")

    # Row A - optimized geometries, flagged if they are not true minima
    for j, n in enumerate(names):
        ax = fig.add_subplot(gs[0, j])
        e = en.get(n, {})
        nim = e.get("n_imag")
        suspect = bool(nim) and not pd.isna(nim)
        show_image(ax, panels[n]["geom"], title=ctx.label(n))
        if suspect and e.get("flat_mode"):
            # A few tens of wavenumbers on a methyl or hydroxyl rotor is a flat coordinate,
            # not a reaction saddle. Calling it "NOT A MINIMUM" in red overstates it.
            lo = e.get("lowest_freq")
            txt = ("flat mode " + (f"{float(lo):.0f} cm$^{{-1}}$" if lo else "")
                   + " - free internal rotation")
            cell_caption(fig, ax, txt, size=7, color="#C77800")
        elif suspect:
            badge = cell_caption(fig, ax, f"NOT A MINIMUM - {int(nim)} imaginary mode(s)",
                                 size=7.5, color="white",
                                 bbox={"boxstyle": "round,pad=0.25", "fc": "#C62828",
                                           "ec": "none"})
            if badge is not None:
                badge.set_fontweight("bold")
        elif e.get("is_min") is True and e.get("restarts"):
            # The structure is a genuine minimum, but it is not the one the first
            # optimization returned - say so rather than let the restart disappear.
            # This counter used to be called `n`, shadowing the compound name for the rest
            # of the iteration; it survived only because nothing read `n` after it.
            k = e["restarts"]
            cell_caption(fig, ax, f"minimum after {k} restart{'s' if k > 1 else ''}",
                         size=7, color="#2F6FA8")
        elif e.get("is_min") is True:
            cell_caption(fig, ax, "minimum confirmed", size=7, color="#2E7D32")

    # Row B - LUMO above HOMO, with the gap between them
    for j, n in enumerate(names):
        axl = fig.add_subplot(gs[1, j])
        e = en.get(n, {})
        lab = "LUMO" if pd.isna(e.get("lumo", np.nan)) else f"LUMO  {e['lumo']:.3f} eV"
        show_image(axl, panels[n]["lumo"], title=lab, title_size=8)

        hlab = "HOMO" if pd.isna(e.get("homo", np.nan)) else f"HOMO  {e['homo']:.3f} eV"
        axh = fig.add_subplot(gs[2, j])
        show_image(axh, panels[n]["homo"], below=hlab, below_size=8)

        if not pd.isna(e.get("gap", np.nan)):
            unbound = bool(e.get("lumo_unbound"))
            axl.annotate(
                f"$\\Delta E_{{gap}}$ = {e['gap']:.3f} eV"
                + ("\nLUMO unbound - not an excitation energy" if unbound else ""),
                xy=(0.5, -0.06), xycoords="axes fraction", ha="center", va="top",
                fontsize=8 if not unbound else 7,
                color="#8A4B00" if unbound else "black",
                bbox={"boxstyle": "round,pad=0.22",
                          "fc": "#FFF6E5" if unbound else "#F2EEFF",
                          "ec": "#C77800" if unbound else "#7C4DFF", "lw": 0.7})

    # Row C - MEP
    for j, n in enumerate(names):
        ax = fig.add_subplot(gs[3, j])
        t = ""
        if n in per_panel_range:
            t = f"{per_panel_range[n][0]:+.3f} to {per_panel_range[n][1]:+.3f} a.u."
        show_image(ax, panels[n]["mep"], title=t, title_size=7)

    # Shared MEP colour bar
    cax = fig.add_subplot(gs[4, :])
    grad = np.linspace(-1, 1, 512)[None, :]
    cax.imshow(grad, aspect="auto", cmap=ctx.esp_colormap(),
               extent=(vrange[0], vrange[1], 0, 1))
    cax.set_yticks([])
    ticks = np.linspace(vrange[0], vrange[1], 5)
    cax.set_xticks(ticks)
    cax.set_xticklabels([f"{t:+.3f}" for t in ticks])
    note = {"per_panel": "   (each panel scaled separately - NOT comparable)",
            "auto_shared": "   (range measured from the surfaces, shared by all panels)",
            "fixed": ""}.get(mode, "")
    cax.set_xlabel("electrostatic potential V(r) / a.u." + note, fontsize=7, labelpad=1)
    for s in ("top", "right", "left"):
        cax.spines[s].set_visible(False)
    cax.text(vrange[0], 1.35, "electron-rich", fontsize=6.5, ha="left", va="bottom")
    cax.text(vrange[1], 1.35, "electron-poor", fontsize=6.5, ha="right", va="bottom")

    band_letter(fig, "A", 0.952, label_style)
    band_letter(fig, "B", 0.755, label_style)
    band_letter(fig, "C", 0.330, label_style)

    fig.suptitle(
        f"Optimized geometry, frontier orbitals and electrostatic potential in "
        f"{solvent} ($\\epsilon$ = {ctx.dielectric():g})",
        fontsize=9.5, y=0.985, x=0.055, ha="left")
    mep_note = ("scaled per panel" if mode == "per_panel"
                else f"range {vrange[0]:+.3f} to {vrange[1]:+.3f} a.u., shared")
    fig.text(0.055, 0.021,
             f"{ctx.level_of_theory()}   |   orbital isosurface "
             f"$\\pm${ctx.param('render.orbital_isovalue')} a.u."
             f"   |   density isosurface {ctx.param('render.density_isovalue')} a.u."
             f"   |   MEP {mep_note}",
             fontsize=6.5, ha="left", va="bottom", color="#707070")
    # Provenance on its own line rather than appended to the method line: that line is
    # already near the figure width at print size, and the two identifiers below are what
    # a reader needs to find the run - the config hash matches the config_sha column in
    # descriptors.csv and the one in runs/manifest.json. Two short fields, no date and no
    # file paths: anything longer stops being read.
    fig.text(0.055, 0.006,
             f"dftfig {ctx.provenance['dftfig_version']}   |   "
             f"config {ctx.provenance['config_sha']}",
             fontsize=6.5, ha="left", va="bottom", color="#707070")

    return ctx.save(fig, "frontier_orbitals")
