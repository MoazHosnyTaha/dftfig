"""Colour theme and console helpers shared by the whole pipeline."""
from __future__ import annotations

from rich.console import Console
from rich.theme import Theme

#: Palette used for every panel, table and prompt in the wizard.
#: Change these six entries to restyle the whole interface at once.
PALETTE = {
    "brand":    "bold #7C4DFF",   # headings / pipeline name
    "ok":       "bold #00C853",   # success, converged, "recommended"
    "warn":     "bold #FFAB00",   # caution, cost warnings
    "bad":      "bold #FF5252",   # errors, imaginary frequencies
    "info":     "#4FC3F7",        # parameter help text
    "muted":    "#90A4AE",        # units, defaults, secondary notes
    "value":    "bold #FFD740",   # current/chosen values
    "paper":    "bold #FF80AB",   # reserved
    "stage":    "bold #26C6DA",   # pipeline stage banners
    "suggested": "bold #00C853",  # the choice to take if you have no reason not to
    "cheap":    "bold #FFAB00",   # the choice that costs least, not the one to publish
}

THEME = Theme({
    "brand": PALETTE["brand"], "ok": PALETTE["ok"], "warn": PALETTE["warn"],
    "bad": PALETTE["bad"], "info": PALETTE["info"], "muted": PALETTE["muted"],
    "value": PALETTE["value"], "paper": PALETTE["paper"], "stage": PALETTE["stage"],
    "suggested": PALETTE["suggested"], "cheap": PALETTE["cheap"],
})

console = Console(theme=THEME, highlight=False)


def banner(text: str, subtitle: str = "") -> None:
    """Print a stage banner, e.g. at the start of the geometry-optimization stage."""
    from rich.panel import Panel
    body = f"[stage]{text}[/stage]"
    if subtitle:
        body += f"\n[muted]{subtitle}[/muted]"
    console.print(Panel(body, border_style=PALETTE["stage"], expand=False))


def ok(msg: str) -> None:
    console.print(f"  [ok]OK[/ok]  {msg}")


def warn(msg: str) -> None:
    console.print(f"  [warn]! [/warn]  {msg}")


def bad(msg: str) -> None:
    console.print(f"  [bad]X [/bad]  {msg}")


def step(msg: str) -> None:
    console.print(f"  [info]->[/info]  {msg}")


def note(msg: str) -> None:
    """A value the pipeline supplied for you. Not a problem, but you should know."""
    console.print(f"  [info]i [/info]  {msg}")
