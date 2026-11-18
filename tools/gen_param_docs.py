"""Regenerate docs/parameter_reference.md from the registry in dftfig/params.py.

Run from the package root:  python tools/gen_param_docs.py
The registry is the single source of truth; this file only formats it.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dftfig.params import GROUPS, TAG_FAST, TAG_PAPER, TAG_RECOMMENDED  # noqa: E402

TAG_NAME = {TAG_PAPER: "paper", TAG_RECOMMENDED: "recommended", TAG_FAST: "fast"}


def fmt_param(p) -> list[str]:
    head = f"### `{p.key}`" + ("  *(advanced)*" if p.advanced else "")
    meta = [f"**{p.label}**"]
    if p.unit:
        meta.append(f"\u2014 {p.unit}")
    meta.append(f"\u00b7 type `{p.kind}`")
    meta.append(f"\u00b7 default `{p.default}`")
    if p.minimum is not None or p.maximum is not None:
        rng = f"{p.minimum if p.minimum is not None else ''}\u2013{p.maximum if p.maximum is not None else ''}"
        meta.append(f"\u00b7 range `{rng}`")
    out = [head, " ".join(meta), ""]
    if p.help:
        out += [p.help, ""]
    if p.guidance:
        out += [f"**How to choose.** {p.guidance}", ""]
    if p.cost:
        out += [f"**Cost.** {p.cost}", ""]
    if p.choices:
        out += ["| value | name | what it means | preset |", "|---|---|---|---|"]
        for c in p.choices:
            tags = " ".join(f"`{TAG_NAME.get(t, t)}`" for t in c.tags)
            detail = c.detail.replace("|", "\\|")
            out.append(f"| `{c.value}` | {c.label} | {detail} | {tags} |")
        out.append("")
    return out


def main() -> None:
    lines = [
        "# dftfig \u2014 complete parameter reference",
        "Generated from the registry in `dftfig/params.py` by `tools/gen_param_docs.py`. "
        "`dft params` shows exactly this information, one parameter at a time, in colour.",
        "`paper` = reproduces the chlocarbazomycin manuscript \u00b7 "
        "`recommended` = what I would run for new work \u00b7 `fast` = plumbing test",
        "", "---", "",
    ]
    n = 0
    for g in GROUPS:
        lines += [f"## {g.title}  (`{g.name}`)", f"*{g.blurb}*", ""]
        for p in g.params:
            lines += fmt_param(p)
            n += 1
    out = Path(__file__).resolve().parent.parent / "docs" / "parameter_reference.md"
    out.parent.mkdir(exist_ok=True)
    out.write_text("\n".join(lines).rstrip() + "\n")
    print(f"{n} parameters -> {out}")


if __name__ == "__main__":
    main()
