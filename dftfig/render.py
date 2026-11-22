"""Headless 3D rendering of orbitals, MEP and lipophilicity surfaces.

Scripts are written to disk next to the images they produce, so every panel in a figure
can be traced back to - and re-run from - the exact commands that drew it.
"""
from __future__ import annotations

import subprocess
import time
from pathlib import Path

from .theme import bad, console, note, warn

#: The named PyMOL colours the ramps below are written in, as RGB. A matplotlib colour
#: bar has to be built from the same numbers the 3D panel was painted with, or the bar
#: under the figure means something slightly different from the surfaces above it.
PYMOL_RGB = {
    "red": (1.0, 0.0, 0.0), "orange": (1.0, 0.5, 0.0), "yellow": (1.0, 1.0, 0.0),
    "green": (0.0, 1.0, 0.0), "cyan": (0.0, 1.0, 1.0), "blue": (0.0, 0.0, 1.0),
    "white": (1.0, 1.0, 1.0), "purple": (0.5, 0.0, 0.5),
}

#: ESP colour ramps. Values are fractions of render.esp_range.
ESP_PALETTES = {
    # The manuscript's sequence: red = electron-rich (negative V), blue = electron-poor.
    "paper":   ([-1.0, -0.5, 0.0, 0.5, 1.0], ["red", "yellow", "green", "cyan", "blue"]),
    # Diverging, perceptually sane, unambiguous zero.
    "rwb":     ([-1.0, 0.0, 1.0], ["red", "white", "blue"]),
    "rainbow": ([-1.0, -0.5, 0.0, 0.5, 1.0], ["red", "orange", "green", "blue", "purple"]),
}

MLP_PALETTES = {
    "paper":        ["red", "white", "blue"],        # red = polar, blue = lipophilic
    "conventional": ["blue", "white", "orange"],     # blue = polar, orange = lipophilic
}


def find_pymol(cfg: dict) -> str | None:
    """Locate a PyMOL executable (one implementation, in `discover`)."""
    from .discover import find_pymol as _find
    p = _find(cfg["render"].get("pymol_path") or "")
    return str(p) if p else None


def _run_pymol(cfg: dict, script: Path, expect: Path | None = None,
               label: str = "") -> bool:
    """Run a .pml headlessly and confirm it actually produced its image.

    PyMOL launchers can exit 0 without rendering anything, so the exit code is not
    sufficient evidence - the output file is.

    PyMOL is silent until it finishes, and a dense cube on a drug-sized molecule can
    take a while, so a spinner runs for the duration. Without it the figures stage
    looks identical to a hang - which is exactly how it was first reported.
    """
    exe = find_pymol(cfg)
    if exe is None:
        bad("PyMOL executable not found. Set render.pymol_path, or install it with "
            "`conda install -c conda-forge pymol-open-source`.")
        return False
    what = label or script.stem
    t0 = time.time()
    with console.status(f"[muted]rendering {what}[/muted]", spinner="dots"):
        proc = subprocess.run([exe, "-cq", script.name], capture_output=True, text=True,
                              cwd=str(script.parent))
    dt = time.time() - t0
    # Only the slow ones are worth a line; a fast render would just add noise.
    if dt >= 5.0:
        note(f"rendered {what} in {dt:.0f} s")
    tail = (proc.stderr or "") + (proc.stdout or "")
    if proc.returncode != 0:
        bad(f"PyMOL exited {proc.returncode} on {script.name}: {tail[-400:]}")
        return False
    if expect is not None and not Path(expect).exists():
        bad(f"PyMOL ran but produced no image for {script.name} (using {exe}). "
            f"Last output: {tail[-400:]}")
        return False
    return True


def _preamble(cfg: dict) -> list[str]:
    r = cfg["render"]
    bg = {"white": "white", "black": "black", "transparent": "white"}[r["background"]]
    lines = [
        "from pymol import cmd",
        "cmd.reinitialize()",
        f'cmd.bg_color("{bg}")',
        f'cmd.set("ray_opaque_background", {0 if r["background"] == "transparent" else 1})',
        'cmd.set("antialias", 2)',
        'cmd.set("orthoscopic", 1)',
        'cmd.set("ray_trace_mode", 0)',
        'cmd.set("specular", 0.2)',
        'cmd.set("ambient", 0.25)',
        'cmd.set("direct", 0.5)',
        'cmd.set("reflect", 0.4)',
        'cmd.set("depth_cue", 0)',
        'cmd.set("ray_shadow", 0)',
        f'cmd.set("surface_quality", {int(r["surface_quality"])})',
        'cmd.set("mesh_quality", 2)',
        'cmd.set("transparency_mode", 2)',
        'cmd.set("surface_color_smoothing", 1)',
        f'cmd.set("surface_color_smoothing_threshold", {float(r["surface_color_smoothing"])})',
        # A translucent surface is seen from both sides; without this its interior renders
        # unlit and the molecule inside looks like it is sitting in a dark cavity.
        'cmd.set("two_sided_lighting", 1)',
        'cmd.set("ray_interior_color", "grey70")',
        'cmd.set("ray_transparency_contrast", 0.6)',
        f'cmd.set("stick_radius", {float(r["stick_radius"])})',
        f'cmd.set("sphere_scale", {float(r["sphere_scale"])})',
        'cmd.set("label_size", 22)',
        'cmd.set("label_color", "black")',
        'cmd.set("label_font_id", 7)',
    ]
    return lines


def _load_molecule(cfg: dict, xyz: Path, obj: str = "mol") -> list[str]:
    r = cfg["render"]
    lines = [
        f'cmd.load("{xyz.name}", "{obj}")',
        f'cmd.hide("everything", "{obj}")',
        f'cmd.show("sticks", "{obj}")',
        f'cmd.show("spheres", "{obj}")',
        f'cmd.color("grey30", "{obj} and elem C")',
        f'cmd.color("white", "{obj} and elem H")',
        f'cmd.color("tv_blue", "{obj} and elem N")',
        f'cmd.color("firebrick", "{obj} and elem O")',
        f'cmd.color("green", "{obj} and elem Cl")',
    ]
    if not r["show_hydrogens"]:
        lines.append(f'cmd.hide("everything", "{obj} and elem H")')
    return lines


def _orient(cfg: dict, obj: str = "mol", mcs_indices: list[int] | None = None) -> list[str]:
    mode = cfg["render"]["orientation"]
    if mode == "none":
        return [f'cmd.zoom("{obj}", 2.0)']
    if mode == "first_molecule" and mcs_indices:
        sel = "+".join(str(i) for i in mcs_indices)
        return [f'cmd.orient("{obj} and index {sel}")', f'cmd.zoom("{obj}", 2.0)']
    return [f'cmd.orient("{obj}")', f'cmd.zoom("{obj}", 2.0)']


def _png(cfg: dict, out: Path) -> list[str]:
    r = cfg["render"]
    w, h = int(r["width"]), int(r["height"])
    if r["ray_trace"]:
        return [f'cmd.png("{out.name}", width={w}, height={h}, dpi=300, ray=1)']
    return [f'cmd.png("{out.name}", width={w}, height={h}, dpi=300, ray=0)']


def _labels(cfg: dict, obj: str = "mol") -> list[str]:
    raw = str(cfg["render"]["label_atoms"]).strip()
    if not raw:
        return []
    idx = [s.strip() for s in raw.split(",") if s.strip().isdigit()]
    if not idx:
        return []
    sel = "+".join(idx)
    return [f'cmd.label("{obj} and index {sel}", "\'%s%s\' % (elem, rank+1)")']


def render_geometry(cfg: dict, name: str, solvent: str, mcs_indices=None) -> Path | None:
    """Ball-and-stick render of the DFT-optimized geometry (the frontier_orbitals figure, row A)."""
    work = Path(cfg["execution"]["work_dir"]) / name / solvent
    xyz = work / f"{name}_{solvent}_sp.xyz"
    if not xyz.exists():
        xyz = work / f"{name}_{solvent}_opt.xyz"
    out = work / f"render_geom_{name}.png"
    if not xyz.exists():
        warn(f"{name}/{solvent}: no optimized geometry to render")
        return None
    if out.exists() and cfg["execution"]["resume"]:
        return out

    s = _preamble(cfg) + _load_molecule(cfg, xyz) + _labels(cfg)
    s += _orient(cfg, "mol", mcs_indices) + _png(cfg, out)
    script = work / f"render_geom_{name}.pml"
    script.write_text("\n".join(s) + "\n")
    return out if _run_pymol(cfg, script, out, f"{name} geometry") else None


def render_orbital(cfg: dict, name: str, solvent: str, which: str,
                   mcs_indices=None) -> Path | None:
    """Two-phase isosurface of one molecular orbital (the frontier_orbitals figure, row B)."""
    work = Path(cfg["execution"]["work_dir"]) / name / solvent
    cube = work / f"{which}.cube"
    xyz = work / f"{name}_{solvent}_sp.xyz"
    out = work / f"render_{which}_{name}.png"
    if not cube.exists():
        warn(f"{name}/{solvent}: {cube.name} missing")
        return None
    if out.exists() and cfg["execution"]["resume"]:
        return out

    pos, neg = [c.strip() for c in str(cfg["render"]["orbital_colors"]).split(",")][:2]
    iso = float(cfg["render"]["orbital_isovalue"])
    s = _preamble(cfg)
    if xyz.exists():
        s += _load_molecule(cfg, xyz)
    s += [
        f'cmd.load("{cube.name}", "orb")',
        f'cmd.isosurface("orb_p", "orb", {iso})',
        f'cmd.isosurface("orb_n", "orb", {-iso})',
        f'cmd.color("{pos}", "orb_p")',
        f'cmd.color("{neg}", "orb_n")',
        'cmd.set("transparency", 0.15, "orb_p")',
        'cmd.set("transparency", 0.15, "orb_n")',
    ]
    s += _orient(cfg, "mol" if xyz.exists() else "orb_p", mcs_indices) + _png(cfg, out)
    script = work / f"render_{which}_{name}.pml"
    script.write_text("\n".join(s) + "\n")
    return out if _run_pymol(cfg, script, out, f"{name} {which.upper()}") else None


def render_mep(cfg: dict, name: str, solvent: str,
               vrange: float | tuple[float, float] | None = None,
               mcs_indices=None) -> Path | None:
    """Electrostatic potential painted on the electron-density isosurface (the frontier_orbitals figure, row C)."""
    work = Path(cfg["execution"]["work_dir"]) / name / solvent
    dens, esp = work / "eldens.cube", work / "esp.cube"
    xyz = work / f"{name}_{solvent}_sp.xyz"
    out = work / f"render_mep_{name}.png"
    if not (dens.exists() and esp.exists()):
        warn(f"{name}/{solvent}: need both eldens.cube and esp.cube for the MEP")
        return None
    if out.exists() and cfg["execution"]["resume"]:
        return out

    if vrange is None:
        vrange = float(cfg["render"]["esp_range"])
    # A symmetric range is a special case of (lo, hi). Ions need the asymmetric form:
    # the potential on a cation's surface never reaches zero, so forcing symmetry about
    # zero throws away the whole dynamic range.
    if isinstance(vrange, (tuple, list)):
        lo_v, hi_v = float(vrange[0]), float(vrange[1])
    else:
        lo_v, hi_v = -float(vrange), float(vrange)
    stops, colors = ESP_PALETTES[cfg["render"]["esp_palette"]]
    vals = [f"{lo_v + (s + 1.0) * 0.5 * (hi_v - lo_v):.6f}" for s in stops]
    cols = ", ".join(f'"{c}"' for c in colors)
    iso = float(cfg["render"]["density_isovalue"])
    transp = float(cfg["render"]["surface_transparency"])

    s = _preamble(cfg)
    if xyz.exists():
        s += _load_molecule(cfg, xyz)
        s += ['cmd.set("stick_radius", 0.07)', 'cmd.set("sphere_scale", 0.12)']
    s += [
        f'cmd.load("{dens.name}", "dens")',
        f'cmd.load("{esp.name}", "esppot")',
        f'cmd.isosurface("surf", "dens", {iso})',
        f'cmd.ramp_new("espramp", "esppot", [{", ".join(vals)}], [{cols}])',
        'cmd.set("surface_color", "espramp", "surf")',
        f'cmd.set("transparency", {transp}, "surf")',
        'cmd.disable("espramp")',
    ]
    s += _orient(cfg, "mol" if xyz.exists() else "surf", mcs_indices) + _png(cfg, out)
    script = work / f"render_mep_{name}.pml"
    script.write_text("\n".join(s) + "\n")
    return out if _run_pymol(cfg, script, out, f"{name} electrostatic potential") else None


def render_mlp(cfg: dict, name: str, pdb: Path, lo: float, hi: float,
               mcs_indices=None) -> Path | None:
    """Molecular surface coloured by per-atom lipophilicity held in the B-factor column."""
    work = Path(pdb).parent
    out = work / f"render_mlp_{name}.png"
    if out.exists() and cfg["execution"]["resume"]:
        return out
    colors = MLP_PALETTES[cfg["figures"]["mlp_palette"]]
    ramp = "_".join(colors)
    transp = float(cfg["figures"]["mlp_transparency"])
    show_struct = bool(cfg["figures"]["mlp_show_structure"])
    name_pdb = Path(pdb).name

    # The surface and the skeleton are separate objects: cmd.spectrum colours every
    # representation of an object, so a single object would ramp the sticks as well and
    # lose the element colours that make the panel readable.
    s = _preamble(cfg) + [
        f'cmd.load("{name_pdb}", "surf")',
        'cmd.hide("everything", "surf")',
        'cmd.show("surface", "surf")',
        f'cmd.spectrum("b", "{ramp}", "surf", minimum={lo:.4f}, maximum={hi:.4f})',
        f'cmd.set("transparency", {transp}, "surf")',
    ]
    if show_struct:
        s += [
            f'cmd.load("{name_pdb}", "mol")',
            'cmd.hide("everything", "mol")',
            'cmd.show("sticks", "mol")',
            'cmd.set("stick_radius", 0.07, "mol")',
            'cmd.color("grey25", "mol and elem C")',
            'cmd.color("grey60", "mol and elem H")',
            'cmd.color("tv_blue", "mol and elem N")',
            'cmd.color("firebrick", "mol and elem O")',
            'cmd.color("green", "mol and elem Cl")',
        ]
    s += _orient(cfg, "surf", mcs_indices) + _png(cfg, out)
    script = work / f"render_mlp_{name}.pml"
    script.write_text("\n".join(s) + "\n")
    return out if _run_pymol(cfg, script, out, f"{name} lipophilicity surface") else None


def write_sp_geometry(cfg: dict, name: str, solvent: str) -> Path | None:
    """Write the optimized geometry used by the property step as an .xyz for rendering."""
    from . import parsers
    work = Path(cfg["execution"]["work_dir"]) / name / solvent
    src = work / f"{name}_{solvent}_sp.out"
    if not src.exists():
        src = work / f"{name}_{solvent}_opt.out"
    if not src.exists():
        return None
    syms, xyz = parsers.geometry(src)
    if not syms:
        return None
    target = work / f"{name}_{solvent}_sp.xyz"
    lines = [str(len(syms)), f"{name} optimized geometry ({solvent})"]
    for sym, (x, y, z) in zip(syms, xyz):
        lines.append(f"{sym:<3s} {x:>14.8f} {y:>14.8f} {z:>14.8f}")
    target.write_text("\n".join(lines) + "\n")
    return target
