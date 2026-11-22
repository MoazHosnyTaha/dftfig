# Writing a figure

dftfig's four built-in figures are a medicinal-chemistry story: frontier orbitals, pKa,
lipophilicity surfaces, CNS MPO. Your field almost certainly wants different panels. This
page builds one from scratch - a condensed-Fukui bar chart - and runs it. The whole
figure is under 60 lines, and dftfig is not modified at all.

The contract itself is documented in `dftfig/figureapi.py`; `dft params figures.which`
lists everything your installation can currently draw.

---

## 1. What a figure is

A figure is a Python module that declares four things:

```python
NAME     = "fukui_bars"                      # the string you put in figures.which
TITLE    = "Condensed Fukui indices per heavy atom"
REQUIRES = (REQ_FUKUI,)                      # what the run must have produced
def build(ctx): ...                          # draws it, returns the files written
```

`build` takes exactly one argument - a `FigureContext` - and returns the list of paths it
wrote. Everything it is allowed to read comes off that context. It never opens a file by
name and never reads `cfg["execution"]["work_dir"]`: where a run keeps its files is
dftfig's business, and the context is how that stays true.

`REQUIRES` is checked **before** `build` is called. A figure that needs PyMOL on a
machine without PyMOL is skipped with that sentence printed, rather than raising from
somewhere inside the renderer twenty minutes into a run. The names are:

| requirement | satisfied when |
| --- | --- |
| `REQ_RENDERING` | a 3D engine is configured and found (PyMOL) |
| `REQ_CUBES` | `cubes.enabled` - volumetric data was produced |
| `REQ_DESCRIPTORS` | the DFT descriptor table is not empty |
| `REQ_PROPERTIES` | the physicochemical property table is not empty |
| `REQ_FUKUI` | `fukui.enabled` - condensed Fukui indices were computed |

Declare only what you genuinely cannot draw without. A panel that degrades gracefully
should check at draw time instead and say so on the figure.

## 2. What the context gives you

```python
ctx.names                      # compounds in run order - your column order
ctx.compound(name)             # .mol (RDKit, explicit H), .charge, .multiplicity,
                               #   .formula, .run_dir
ctx.reference_solvent          # the medium the figures are built from
ctx.dielectric()               # its permittivity
ctx.level_of_theory()          # "B3LYP-D3BJ/def2-SVP opt // .../def2-TZVP property, CPCM"
ctx.provenance                 # {"dftfig_version": ..., "config_sha": ...}
ctx.param("figures.dpi")       # any config value, by dotted key

ctx.descriptor_row(name)       # a row of descriptors.csv (Python scalars, not numpy)
ctx.descriptor(name, "gap_eV") # one value
ctx.property_row(name)         # the resolved physicochemical row
ctx.property(name, "TPSA")     # one value
ctx.fukui(name)                # per-atom Fukui table, or None
ctx.som_table(name)            # an external site-of-metabolism export, or None

ctx.cube(name, "homo")         # locate a cube file, or None - never build the path
ctx.surface_esp_stats(name)    # V(r) on the molecular surface
ctx.image(name, "mep", vrange=(lo, hi))     # a rendered 3D panel as a PNG, or None
ctx.lipophilicity_range(name)  # so every panel can share one colour scale

ctx.apply_style()              # dftfig's matplotlib style at figures.dpi
ctx.esp_colormap()             # colour bar matching the 3D MEP panels
ctx.save(fig, "fukui_bars")    # writes every figures.formats entry; returns the paths
```

Atom order matters and the contract pins it down: `compound.mol` is in the same order as
the Fukui table and the geometry ORCA reported, hydrogens included. A per-atom value can
therefore be painted onto a depiction without a mapping step.

## 3. The figure

Condensed Fukui indices are a good first figure: the `descriptors` stage already computed
them, so there is nothing to run, and the panel needs no 3D rendering.

Create `dft_figures/fukui_bars.py` in your project directory:

```python
"""Condensed Fukui indices as a bar panel, one column per compound."""
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
```

A copy ready to run is in `examples/fukui_bars/`.

## 4. See it discovered

```console
$ dft params figures.which
  ...
  value         name                                        what it means
* fig2          Figure 2 - geometry / HOMO-LUMO / MEP       ...
  ...
  fukui_bars    Condensed Fukui indices per heavy atom      plugin, from
                                                            dft_figures/fukui_bars.py;
                                                            needs fukui
```

Nothing was installed and nothing in dftfig was edited. The interview
(`dft run`, step "Figure assembly") lists it in the same place.

If it does not appear, dftfig will have printed why - a syntax error, a missing `build`,
or an unknown requirement name. A broken plugin is reported and skipped; it never stops
a run.

## 5. Run it

```yaml
# config.yaml
figures:
  which: [fig2, fig12, fig13, fig14, fukui_bars]
```

```console
$ dft figures
 figures
 fig2, fig12, fig13, fig14, fukui_bars
  ...
  OK  fukui_bars -> runs/figures/fukui_bars.png
```

That is the whole loop: write the file, see it listed, run it, get a PNG at
`figures.dpi` in every format in `figures.formats`.

## 6. Sharing it as a package

A local file is right for a figure that belongs to one project. To ship one, register the
`dftfig.figures` entry-point group - the same mechanism dftfig uses for its own four:

```toml
# pyproject.toml of your package
[project.entry-points."dftfig.figures"]
fukui_bars = "mypkg.figures.fukui_bars"
```

`pip install mypkg` and it appears in `dft params figures.which` with
`plugin, from package 'mypkg'`. A project-local file of the same `NAME` still wins, so a
user can override your figure without uninstalling it.

## 7. Versioning

`dftfig.figureapi.API_VERSION` is the contract's version, currently **1**. It is bumped
only when an existing plugin could stop working - a method removed or renamed, an
argument's meaning changed, a field dropped from `Compound`, `build`'s signature altered.
New accessors, new requirement names and new optional keyword arguments do not bump it.

If your figure needs a guarantee, declare it:

```python
API_VERSION = 1      # refuse to load on a dftfig that provides an older contract
```

dftfig refuses to load a plugin that asks for a version higher than it provides, and says
so, rather than failing later with an AttributeError.

## 8. Rules of thumb

1. If you need something the context does not expose, that is a gap in the contract -
   open an issue. Reaching around it into `cfg["execution"]["work_dir"]` works today and
   breaks on the next release.
2. Declare `REQUIRES` honestly, and handle *partial* data inside `build`: a run where one
   compound out of five has no Fukui table is normal, and the panel should say so rather
   than vanish.
3. Return every file you wrote. `dft verify` hashes them, and the provenance manifest is
   built from that list.
4. Put the method and the provenance stamp on the figure. A figure that leaves the
   repository without them cannot be checked by a reader.
