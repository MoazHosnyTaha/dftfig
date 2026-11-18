"""Figure 2 - optimized geometries, frontier orbitals with the gap, and MEP surfaces."""
from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.gridspec import GridSpec

from .. import cubes, render
from ..figstyle import (apply_style, band_letter, caption, cell_caption,
                        pad_to_common_aspect,
                        show_image, trim_white)
from ..orca_io import DIELECTRIC
from ..render import ESP_PALETTES
from ..theme import ok, step, warn

PYMOL_RGB = {
    "red": (1.0, 0.0, 0.0), "orange": (1.0, 0.5, 0.0), "yellow": (1.0, 1.0, 0.0),
    "green": (0.0, 1.0, 0.0), "cyan": (0.0, 1.0, 1.0), "blue": (0.0, 0.0, 1.0),
    "white": (1.0, 1.0, 1.0), "purple": (0.5, 0.0, 0.5),
}


def _esp_colormap(palette: str) -> LinearSegmentedColormap:
    stops, colors = ESP_PALETTES[palette]
    pos = [(s + 1.0) / 2.0 for s in stops]
    return LinearSegmentedColormap.from_list(
        f"esp_{palette}", list(zip(pos, [PYMOL_RGB[c] for c in colors])))


def build(cfg: dict, mols: dict, descriptors: pd.DataFrame | None = None) -> list[Path]:
    """Render every panel, then assemble Figure 2."""
    apply_style(int(cfg["figures"]["dpi"]))
    names = list(mols)
    solvent = cfg["figures"]["reference_solvent"]
    work = Path(cfg["execution"]["work_dir"])

    # --- MEP colour range -------------------------------------------------------------
    # The potential around an ion is dominated by its monopole: a +1 cation is positive
    # essentially everywhere and far outside the +-0.05 a.u. that suits neutral organics,
    # so a fixed range renders it as a saturated, featureless blob. Measure the surface
    # potential and scale accordingly.
    mode = cfg["render"].get("esp_scale", "auto_shared")
    if cfg["render"].get("esp_autoscale"):          # legacy flag still honoured
        mode = "per_panel"
    r0 = float(cfg["render"]["esp_range"])
    stats = {n: cubes.surface_esp_stats(cfg, n, solvent) for n in names}
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

    panels: dict[str, dict] = {}
    for n in names:
        render.write_sp_geometry(cfg, n, solvent)
        panels[n] = {
            "geom": render.render_geometry(cfg, n, solvent),
            "homo": render.render_orbital(cfg, n, solvent, "homo"),
            "lumo": render.render_orbital(cfg, n, solvent, "lumo"),
            "mep": render.render_mep(cfg, n, solvent, per_panel_range.get(n, vrange)),
        }
        for p in panels[n].values():
            if p and Path(p).exists() and cfg["render"]["background"] == "transparent":
                trim_white(p)

    # Equalize panel aspect within each row so the labels align across columns.
    for role in ("geom", "homo", "lumo", "mep"):
        pad_to_common_aspect([panels[n][role] for n in names])

    # --- energies ---------------------------------------------------------------------
    en: dict[str, dict] = {}
    if descriptors is not None and not descriptors.empty:
        sub = descriptors[descriptors["medium"] == solvent]
        for _, r in sub.iterrows():
            en[r["molecule"]] = {
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

    # Row A - optimized geometries, flagged if they are not true minima
    for j, n in enumerate(names):
        ax = fig.add_subplot(gs[0, j])
        e = en.get(n, {})
        nim = e.get("n_imag")
        suspect = bool(nim) and not pd.isna(nim)
        show_image(ax, panels[n]["geom"], title=n)
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
                                 bbox=dict(boxstyle="round,pad=0.25", fc="#C62828",
                                           ec="none"))
            if badge is not None:
                badge.set_fontweight("bold")
        elif e.get("is_min") is True and e.get("restarts"):
            # The structure is a genuine minimum, but it is not the one the first
            # optimization returned - say so rather than let the restart disappear.
            n = e["restarts"]
            cell_caption(fig, ax, f"minimum after {n} restart{'s' if n > 1 else ''}",
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
                bbox=dict(boxstyle="round,pad=0.22",
                          fc="#FFF6E5" if unbound else "#F2EEFF",
                          ec="#C77800" if unbound else "#7C4DFF", lw=0.7))

    # Row C - MEP
    for j, n in enumerate(names):
        ax = fig.add_subplot(gs[3, j])
        t = ""
        if n in per_panel_range:
            t = f"{per_panel_range[n][0]:+.3f} to {per_panel_range[n][1]:+.3f} a.u."
        show_image(ax, panels[n]["mep"], title=t, title_size=7)

    # Shared MEP colour bar
    cax = fig.add_subplot(gs[4, :])
    cmap = _esp_colormap(cfg["render"]["esp_palette"])
    grad = np.linspace(-1, 1, 512)[None, :]
    cax.imshow(grad, aspect="auto", cmap=cmap, extent=(vrange[0], vrange[1], 0, 1))
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

    band_letter(fig, "A", 0.952, cfg["figures"]["panel_label_style"])
    band_letter(fig, "B", 0.755, cfg["figures"]["panel_label_style"])
    band_letter(fig, "C", 0.330, cfg["figures"]["panel_label_style"])

    eps = DIELECTRIC.get(solvent, float("nan"))
    level = _level_string(cfg)
    fig.suptitle(
        f"Optimized geometry, frontier orbitals and electrostatic potential in "
        f"{solvent} ($\\epsilon$ = {eps:g})",
        fontsize=9.5, y=0.985, x=0.055, ha="left")
    fig.text(0.055, 0.008,
             f"{level}   |   orbital isosurface $\\pm${cfg['render']['orbital_isovalue']} a.u."
             f"   |   density isosurface {cfg['render']['density_isovalue']} a.u."
             f"   |   MEP {'scaled per panel' if mode == 'per_panel' else 'range %+.3f to %+.3f a.u., shared' % vrange}",
             fontsize=6.5, ha="left", va="bottom", color="#707070")

    return _save(cfg, fig, "figure2")


def _level_string(cfg: dict) -> str:
    t = cfg["theory"]
    disp = "" if t["dispersion"] == "none" else f"-{t['dispersion']}"
    solv = cfg["solvation"]["model"]
    opt = f"{t['functional']}{disp}/{t['basis_opt']}"
    if cfg["property"]["enabled"]:
        return (f"{opt} opt // {t['functional']}{disp}/{cfg['property']['basis_sp']} "
                f"property, {solv}")
    return f"{opt}, {solv}"


def _save(cfg: dict, fig, stem: str) -> list[Path]:
    outdir = Path(cfg["execution"]["work_dir"]) / "figures"
    outdir.mkdir(parents=True, exist_ok=True)
    paths = []
    for fmt in cfg["figures"]["formats"]:
        p = outdir / f"{stem}.{fmt}"
        fig.savefig(p, dpi=int(cfg["figures"]["dpi"]))
        paths.append(p)
    plt.close(fig)
    ok(f"{stem} -> {', '.join(str(p) for p in paths)}")
    return paths
