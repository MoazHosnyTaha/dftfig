"""Arrow-key menus, drawn with rich and read straight from the terminal.

No third-party dependency: the environment this ships into is read-only, so the key
handling is stdlib termios/tty and the drawing is the `rich` that the package already
uses. Everything degrades to a numbered prompt when there is no terminal to read from -
a pipe, a CI job, `dft run < answers.txt`.
"""
from __future__ import annotations

import os
import sys
from dataclasses import dataclass

from rich.console import Group
from rich.live import Live
from rich.panel import Panel
from rich.text import Text

from .theme import PALETTE, console

UP, DOWN, ENTER, QUIT, TOGGLE, HOME, END = "up", "down", "enter", "quit", "toggle", "home", "end"


@dataclass
class Option:
    value: object
    label: str
    detail: str = ""
    tags: tuple = ()
    disabled: bool = False


def interactive_terminal() -> bool:
    """True when we can both read keys and draw over them.

    The key reader below is built on termios/tty, which exist only on POSIX. On
    Windows the module imports fine but `tty.setcbreak` does not, so the check
    includes the platform rather than discovering it halfway through a menu: there,
    and anywhere stdin is a pipe, every question degrades to a numbered prompt that
    answers to the same keystrokes.

    Set DFTFIG_NO_TUI=1 to force that fallback on any platform - useful when a
    terminal emulator mangles the redraw, and required for `dft run < answers.txt`.
    """
    if os.environ.get("DFTFIG_NO_TUI"):
        return False
    if os.name != "posix":
        return False
    try:
        import termios  # noqa: F401 - presence is the test
        import tty      # noqa: F401
    except ImportError:
        return False
    try:
        return sys.stdin.isatty() and sys.stdout.isatty()
    except Exception:
        return False


# --------------------------------------------------------------------------- key reading
class _Keys:
    """Raw-mode keyboard for the lifetime of one menu.

    Raw mode is entered once, not once per keystroke. Toggling it per key leaves the
    terminal in canonical mode in between, where the line discipline holds everything up
    until a newline and then delivers it in one burst - a burst that a single read would
    consume and discard all but the first key of. Anything read beyond the current key is
    kept in `_buf` instead.
    """

    def __init__(self) -> None:
        self.fd = sys.stdin.fileno()
        self._buf = b""
        self._old = None

    def __enter__(self) -> _Keys:
        import termios
        import tty
        self._old = termios.tcgetattr(self.fd)
        # cbreak, NOT raw. Both give unbuffered single-key reads, but setraw also clears
        # OPOST, which is what translates "\n" into carriage-return + line-feed on the way
        # out. Without it every line printed while the menu is open starts at the column
        # the previous line ended on, and the display walks diagonally across the screen.
        # cbreak clears ICANON and ECHO only, so output behaves normally.
        tty.setcbreak(self.fd)
        return self

    def __exit__(self, *exc) -> None:
        import termios
        if self._old is not None:
            termios.tcsetattr(self.fd, termios.TCSADRAIN, self._old)

    def _fill(self, blocking: bool = True) -> None:
        import select
        if not blocking and not select.select([self.fd], [], [], 0)[0]:
            return
        chunk = os.read(self.fd, 64)
        if chunk:
            self._buf += chunk

    def key(self) -> str:
        while True:
            if not self._buf:
                self._fill()
                continue
            if self._buf[:1] == b"\x1b":
                if len(self._buf) == 1:
                    self._fill(blocking=False)
                if len(self._buf) == 1:          # a lone Esc, nothing following
                    self._buf = b""
                    return QUIT
                if self._buf[1:2] in (b"[", b"O"):
                    code, self._buf = self._buf[2:3], self._buf[3:]
                    return {b"A": UP, b"B": DOWN, b"H": HOME, b"F": END}.get(code, "")
                self._buf = self._buf[1:]
                return ""
            ch, self._buf = self._buf[:1], self._buf[1:]
            return _classify(ch)


def _classify(data: bytes) -> str:
    if data in (b"\r", b"\n"):
        return ENTER
    if data == b" ":
        return TOGGLE
    if data in (b"\x03", b"\x04"):           # Ctrl-C, Ctrl-D
        raise KeyboardInterrupt
    ch = data.decode("utf-8", "ignore")
    if ch in ("k", "K"):
        return UP
    if ch in ("j", "J"):
        return DOWN
    if ch in ("q", "Q"):
        return QUIT
    return ch


# --------------------------------------------------------------------------- rendering
def _tag_text(tags) -> Text:
    t = Text()
    for tag in tags:
        t.append(f" {tag} ", style=f"reverse {PALETTE.get(tag, PALETTE['brand'])}")
        t.append(" ")
    return t


def _render(title: str, subtitle: str, options: list[Option], cursor: int,
            chosen: set[int] | None, footer: str) -> Panel:
    rows = []
    for i, o in enumerate(options):
        here = i == cursor
        line = Text()
        line.append("  ")
        line.append(f"{i + 1:>2}. " if i < 9 else "    ",
                    style=PALETTE["value"] if here else PALETTE["muted"])
        if chosen is not None:
            line.append("[x] " if i in chosen else "[ ] ",
                        style=PALETTE["ok"] if i in chosen else PALETTE["muted"])
        line.append("> " if here else "  ", style=PALETTE["brand"] if here else "")
        style = PALETTE["muted"] if o.disabled else (PALETTE["brand"] if here else "")
        line.append(o.label, style=style)
        if o.tags:
            line.append("  ")
            line.append_text(_tag_text(o.tags))
        rows.append(line)
        if here and o.detail:
            d = Text("        " + o.detail, style=PALETTE["muted"])
            rows.append(d)
    body = Group(*rows)
    head = Text(title, style=PALETTE["stage"])
    if subtitle:
        head.append("\n" + subtitle, style=PALETTE["muted"])
    return Panel(Group(head, Text(""), body, Text(""), Text("  " + footer,
                                                           style=PALETTE["muted"])),
                 border_style=PALETTE["brand"], expand=False)


# --------------------------------------------------------------------------- fallbacks
def _numbered(title: str, options: list[Option], default: int) -> int:
    from rich.prompt import Prompt
    console.print(f"\n  [stage]{title}[/stage]")
    for i, o in enumerate(options, 1):
        console.print(f"   {i:>2}. {o.label}" + (f"   [muted]{o.detail}[/muted]"
                                                 if o.detail else ""))
    raw = Prompt.ask("  [brand]choose[/brand]",
                     choices=[str(i) for i in range(1, len(options) + 1)],
                     default=str(default + 1), console=console)
    return int(raw) - 1


# --------------------------------------------------------------------------- public API
def select(title: str, options, subtitle: str = "", default: int = 0,
           allow_quit: bool = True):
    """One choice from a list. Returns the chosen Option.value, or None if quit."""
    options = [o if isinstance(o, Option) else Option(o, str(o)) for o in options]
    if not options:
        return None
    if not interactive_terminal():
        return options[_numbered(title, options, default)].value

    cursor = max(0, min(default, len(options) - 1))
    foot = ("up/down move   enter select   or just type the number"
            + ("   q cancel" if allow_quit else ""))
    with _Keys() as keys, Live(_render(title, subtitle, options, cursor, None, foot),
                               console=console, auto_refresh=False,
                               transient=True) as live:
        while True:
            key = keys.key()
            if key == UP:
                cursor = (cursor - 1) % len(options)
            elif key == DOWN:
                cursor = (cursor + 1) % len(options)
            elif key == HOME:
                cursor = 0
            elif key == END:
                cursor = len(options) - 1
            elif key.isdigit() and key != "0" and int(key) <= len(options):
                cursor = int(key) - 1
                live.update(_render(title, subtitle, options, cursor, None, foot), refresh=True)
                return options[cursor].value
            elif key == ENTER:
                if not options[cursor].disabled:
                    return options[cursor].value
            elif key == QUIT and allow_quit:
                return None
            live.update(_render(title, subtitle, options, cursor, None, foot), refresh=True)


def select_many(title: str, options, subtitle: str = "", preselected=(),
                allow_empty: bool = False):
    """Several choices from a list. Space toggles, enter accepts."""
    options = [o if isinstance(o, Option) else Option(o, str(o)) for o in options]
    chosen = {i for i, o in enumerate(options) if o.value in preselected}
    if not interactive_terminal():
        from rich.prompt import Prompt
        console.print(f"\n  [stage]{title}[/stage]")
        for i, o in enumerate(options, 1):
            console.print(f"   {i:>2}. {o.label}")
        raw = Prompt.ask("  [brand]choose (comma-separated numbers)[/brand]",
                         default=",".join(str(i + 1) for i in sorted(chosen)),
                         console=console)
        idx = [int(x) - 1 for x in raw.replace(" ", "").split(",") if x.isdigit()]
        return [options[i].value for i in idx if 0 <= i < len(options)]

    cursor = 0
    foot = ("up/down move   space toggle   enter accept"
            "   or type a number to toggle it")
    with _Keys() as keys, Live(_render(title, subtitle, options, cursor, chosen, foot),
                               console=console, auto_refresh=False,
                               transient=True) as live:
        while True:
            key = keys.key()
            if key == UP:
                cursor = (cursor - 1) % len(options)
            elif key == DOWN:
                cursor = (cursor + 1) % len(options)
            elif key == TOGGLE:
                chosen.symmetric_difference_update({cursor})
            elif key.isdigit() and key != "0" and int(key) <= len(options):
                chosen.symmetric_difference_update({int(key) - 1})
            elif key == ENTER and (chosen or allow_empty):
                return [options[i].value for i in sorted(chosen)]
            live.update(_render(title, subtitle, options, cursor, chosen, foot), refresh=True)


def confirm(question: str, default: bool = True) -> bool:
    """Yes/no, numbered the same way in both modes.

    The fallback is a numbered prompt rather than rich's y/n so that one set of answers
    drives the pipeline whether or not there is a terminal - `printf '1,3\n1\n2\n2\n' |
    dft run` means the same thing either way.
    """
    opts = [Option(True, "yes"), Option(False, "no")]
    idx = 0 if default else 1
    if not interactive_terminal():
        return bool(opts[_numbered(question, opts, idx)].value)
    return bool(select(question, opts, default=idx, allow_quit=False))
