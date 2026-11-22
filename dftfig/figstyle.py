"""Self-contained figure styling, so the package has no external style dependency."""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl

#: Role-mapped font ladder: (base, annotation, tick).
SIZES = (9, 8, 7)
GREY = "#4A4A4A"


def esp_colormap(palette: str):
    """Matplotlib colormap for the MEP colour bar, built from the PyMOL ramp itself."""
    from matplotlib.colors import LinearSegmentedColormap

    from .render import ESP_PALETTES, PYMOL_RGB
    stops, colors = ESP_PALETTES[palette]
    pos = [(s + 1.0) / 2.0 for s in stops]
    return LinearSegmentedColormap.from_list(
        f"esp_{palette}", list(zip(pos, [PYMOL_RGB[c] for c in colors])))


def mlp_colormap(palette: str):
    """Matplotlib colormap for the lipophilicity colour bar."""
    from matplotlib.colors import LinearSegmentedColormap

    from .render import MLP_PALETTES, PYMOL_RGB
    return LinearSegmentedColormap.from_list(
        f"mlp_{palette}", [PYMOL_RGB[c] for c in MLP_PALETTES[palette]])


def apply_style(dpi: int = 300, sizes: tuple = SIZES) -> None:
    base, anno, tick = sizes
    mpl.rcParams.update({
        "figure.dpi": 110,
        "savefig.dpi": dpi,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
        "font.size": base,
        "axes.titlesize": base,
        "axes.labelsize": base,
        "axes.titleweight": "normal",
        "axes.titlelocation": "left",
        "axes.labelcolor": GREY,
        "axes.edgecolor": GREY,
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": False,
        "legend.fontsize": anno,
        "legend.frameon": False,
        "legend.handlelength": 1.4,
        "xtick.labelsize": tick,
        "ytick.labelsize": tick,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.color": GREY,
        "ytick.color": GREY,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "lines.linewidth": 1.4,
        "lines.markersize": 4.5,
        "text.color": GREY,
    })


def panel_letter(ax, letter: str, style: str = "upper", dx: float = -0.04, dy: float = 1.04):
    if style == "none":
        return
    txt = {"upper": letter.upper(), "lower": letter.lower(),
           "paren_upper": f"({letter.upper()})"}.get(style, letter.upper())
    ax.text(dx, dy, txt, transform=ax.transAxes, fontsize=SIZES[0] + 3,
            fontweight="bold", va="bottom", ha="left", color="black")


def band_letter(fig, letter: str, y: float, style: str = "upper", x: float = 0.012):
    if style == "none":
        return
    txt = {"upper": letter.upper(), "lower": letter.lower(),
           "paren_upper": f"({letter.upper()})"}.get(style, letter.upper())
    fig.text(x, y, txt, fontsize=SIZES[0] + 3, fontweight="bold",
             va="top", ha="left", color="black")


def caption(ax, text: str, size: int | None = None, color: str = GREY, dy: float = -0.035):
    """Text under an axes. Needed because set_xlabel is invisible once the axis is off."""
    if not text:
        return
    ax.text(0.5, dy, text, transform=ax.transAxes, ha="center", va="top",
            fontsize=size or SIZES[1], color=color)


def show_image(ax, path, title: str = "", title_size: int | None = None,
               below: str = "", below_size: int | None = None,
               below_color: str = GREY):
    """Place a rendered PNG in an axes with no frame."""
    import matplotlib.image as mpimg
    ax.set_axis_off()
    if path is None:
        ax.text(0.5, 0.5, "not rendered", transform=ax.transAxes, ha="center", va="center",
                fontsize=SIZES[1], color="#B0B0B0")
    else:
        ax.imshow(mpimg.imread(str(path)))
    if title:
        ax.set_title(title, fontsize=title_size or SIZES[0], loc="center", pad=3)
    caption(ax, below, below_size, below_color)


def cell_caption(fig, ax, text: str, size: int | None = None, color: str = GREY,
                 bbox: dict | None = None, dy: float = 0.004, va: str = "bottom"):
    """Caption pinned to the bottom of the axes' GRID CELL, not to the axes box.

    An imshow axes shrinks to its image's aspect ratio, so for a small or wide molecule
    the axes box sits well above the bottom of its grid cell. Text placed just below the
    axes then floats up into the row beneath and lands on its labels. Anchoring to the
    cell keeps every badge on the same baseline whatever shape the molecule is.
    """
    if not text:
        return None
    spec = ax.get_subplotspec()
    if spec is None:
        return ax.text(0.5, -0.02, text, transform=ax.transAxes, ha="center", va="top",
                       fontsize=size or SIZES[1], color=color, bbox=bbox)
    pos = spec.get_position(fig)
    # dy below 0.5 is a small offset above the cell floor; at or above it, a fraction of
    # the cell height - which is how a title is placed at the top of the same cell.
    y = pos.y0 + dy if dy < 0.5 else pos.y0 + dy * (pos.y1 - pos.y0)
    return fig.text(0.5 * (pos.x0 + pos.x1), y, text, ha="center", va=va,
                    fontsize=size or SIZES[1], color=color, bbox=bbox)


def trim_white(path, pad: int = 8):
    """Crop uniform white/transparent margins from a rendered panel."""
    import numpy as np
    from PIL import Image
    im = Image.open(path)
    arr = np.array(im.convert("RGBA"))
    alpha = arr[..., 3]
    rgb = arr[..., :3]
    content = (alpha > 8) & ~np.all(rgb > 247, axis=-1)
    if not content.any():
        return path
    ys, xs = np.where(content)
    y0, y1 = max(ys.min() - pad, 0), min(ys.max() + pad + 1, arr.shape[0])
    x0, x1 = max(xs.min() - pad, 0), min(xs.max() + pad + 1, arr.shape[1])
    im.crop((x0, y0, x1, y1)).save(path)
    return path


def pad_to_common_aspect(paths, pad_frac: float = 0.02):
    """Pad a row of panel images to one common aspect ratio, in place.

    Matplotlib sets an imshow axes to the image's aspect and shrinks the axes box to
    match. Panels of differing aspect therefore get differently sized boxes, and every
    title, caption and annotation anchored in axes coordinates drifts out of alignment
    across the row. Padding the images to a common aspect first keeps all the boxes
    identical, so the labels line up.
    """
    from PIL import Image

    real = [p for p in paths if p is not None and Path(p).exists()]
    if len(real) < 2:
        return
    sizes = {}
    for p in real:
        with Image.open(p) as im:
            sizes[p] = im.size                      # (w, h)
    target = max(w / h for w, h in sizes.values())
    target = max(target, 1e-6)

    for p in real:
        w, h = sizes[p]
        tw, th = w, h
        if w / h < target:
            tw = int(round(h * target))
        else:
            th = int(round(w / target))
        tw = int(round(tw * (1 + pad_frac)))
        th = int(round(th * (1 + pad_frac)))
        im = Image.open(p).convert("RGBA")
        canvas = Image.new("RGBA", (tw, th), (255, 255, 255, 0))
        # Paste without an alpha mask: using the image as its own mask blends its RGB
        # against the empty canvas and washes out semi-transparent surfaces.
        canvas.paste(im, ((tw - w) // 2, (th - h) // 2))
        canvas.save(p)
