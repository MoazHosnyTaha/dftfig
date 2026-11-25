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


def needs_file_backend(current: object) -> bool:
    """Is matplotlib currently pointed at a backend that wants a display?

    ``current`` is ``rcParams["backend"]``, which is a sentinel object rather than a
    string until matplotlib resolves it - in that case there is nothing to correct yet
    and the environment variable set in ``dftfig/__init__`` decides.
    """
    return isinstance(current, str) and current.lower() in _interactive_backends()


def _interactive_backends() -> frozenset[str]:
    """Names of the backends that want a window, across matplotlib versions.

    matplotlib 3.9 moved this list from ``rcsetup.interactive_bk`` into the backend
    registry, and dftfig has to run on whatever matplotlib the user's environment has.
    """
    from matplotlib import rcsetup
    names = getattr(rcsetup, "interactive_bk", None)
    if names is None:
        try:
            from matplotlib.backends.registry import BackendFilter, backend_registry
            names = backend_registry.list_builtin(BackendFilter.INTERACTIVE)
        except ImportError:          # neither available: fall back to the usual names
            names = ("qtagg", "qt5agg", "tkagg", "wxagg", "macosx", "gtk3agg", "gtk4agg")
    return frozenset(n.lower() for n in names)


def force_file_backend() -> None:
    """Switch away from an interactive backend, late, if one slipped through.

    The environment variable set at package import handles the normal case. This is the
    backstop for when something imported matplotlib first and fixed the backend before
    dftfig was imported at all.
    """
    if needs_file_backend(mpl.rcParams.get("backend")):
        from .theme import note
        was = mpl.rcParams["backend"]
        mpl.use("Agg", force=True)
        note(f"matplotlib was set to the interactive {was} backend; switched to Agg "
             "because dftfig writes figures to files and never opens a window")


def apply_style(dpi: int = 300, sizes: tuple = SIZES) -> None:
    force_file_backend()
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


def footnote_step(fig, size: float) -> float:
    """Line spacing for :func:`footnote`, in figure fraction."""
    return size * 1.6 / (fig.get_figheight() * 72.0)


def wrap_fields(fig, fields, *, size: float, left: float = 0.055,
                right: float = 0.985, sep: str = "   |   ") -> list[str]:
    """Pack ``fields`` into the fewest lines that fit the figure width.

    Separate from :func:`footnote` because a caller has to know how many lines the
    footnote will take *before* it fixes the layout: the lines occupy the bottom
    margin, and a margin sized for one line puts the second through the colour-bar
    label.
    """
    renderer = fig.canvas.get_renderer()
    avail = right - left

    def too_wide(text: str) -> bool:
        probe = fig.text(0, 0, text, fontsize=size)
        w = probe.get_window_extent(renderer=renderer).width / fig.bbox.width
        probe.remove()
        return w > avail

    lines: list[str] = []
    current = ""
    for field in fields:
        candidate = f"{current}{sep}{field}" if current else str(field)
        if current and too_wide(candidate):
            lines.append(current)
            current = str(field)
        else:
            current = candidate
    if current:
        lines.append(current)
    return lines


def wrap_text(fig, text: str, *, size: float, left: float = 0.055,
              right: float = 0.985) -> str:
    """Break `text` on word boundaries into lines that fit the figure width.

    :func:`wrap_fields` packs whole fields and cannot help when one field is itself
    wider than the figure, which is what happens to a sentence-length axis label on a
    narrow canvas. With ``bbox_inches="tight"`` that label does not overflow - it widens
    the saved file, so the panels end up in a canvas sized by a caption.
    """
    renderer = fig.canvas.get_renderer()
    avail = right - left

    def width_of(s: str) -> float:
        probe = fig.text(0, 0, s, fontsize=size)
        w = probe.get_window_extent(renderer=renderer).width / fig.bbox.width
        probe.remove()
        return w

    lines: list[str] = []
    current = ""
    for word in str(text).split():
        candidate = f"{current} {word}" if current else word
        if current and width_of(candidate) > avail:
            lines.append(current)
            current = word
        else:
            current = candidate
    if current:
        lines.append(current)
    return "\n".join(lines)


def shrink_to_width(fig, artist, *, left: float = 0.055, right: float = 0.985,
                    min_size: float = 7.0) -> float:
    """Reduce a text artist's size until it fits the figure width. Returns that size.

    Same reason as :func:`footnote`: with ``bbox_inches="tight"`` an over-wide title
    does not overflow, it silently widens the saved figure, and the panels then sit in
    a canvas sized by a sentence rather than by the data.
    """
    renderer = fig.canvas.get_renderer()
    avail = (right - left) * fig.bbox.width
    have = artist.get_window_extent(renderer=renderer).width
    size = artist.get_fontsize()
    if have > avail:
        size = max(min_size, size * avail / have)
        artist.set_fontsize(size)
    return size


def footnote(fig, fields, *, x: float = 0.055, y: float = 0.021, size: float = 6.5,
             color: str = GREY, right: float = 0.985, sep: str = "   |   ") -> float:
    """Draw a field-separated footnote, wrapped to the width of the figure.

    A long single-line footnote silently decides how wide the saved file is: the style
    saves with ``bbox_inches="tight"``, which grows the output until it contains every
    artist, so a method line wider than the panels pads the figure with empty space on
    the right. How much depends on how many molecules there are, which is why the same
    footnote looks right on four and leaves a third of the canvas blank on two.

    Fields are packed into as many lines as fit and drawn upward from ``y``, so the
    anchor keeps its meaning whatever the wrap produces. Returns the topmost line's y.
    """
    lines = wrap_fields(fig, fields, size=size, left=x, right=right, sep=sep)
    step = footnote_step(fig, size)
    for i, line in enumerate(reversed(lines)):
        fig.text(x, y + i * step, line, fontsize=size, ha="left", va="bottom",
                 color=color)
    return y + (len(lines) - 1) * step


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
