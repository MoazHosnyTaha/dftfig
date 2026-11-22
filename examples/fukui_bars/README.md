# Example plugin: condensed Fukui bars

A complete dftfig figure in under 60 lines of code. It reads the condensed Fukui indices
the `descriptors` stage already wrote and draws one bar panel per compound. No PyMOL, no
new calculation.

## Run it

From a project directory that already has a finished run:

```bash
mkdir -p dft_figures
cp /path/to/dftfig/examples/fukui_bars/fukui_bars.py dft_figures/
dft params figures.which          # fukui_bars is now listed, "plugin, from dft_figures/..."
```

Add it to `config.yaml`:

```yaml
figures:
  which: [frontier_orbitals, speciation, lipophilicity, mpo_metabolism, fukui_bars]
```

and draw it:

```bash
dft figures
# -> runs/figures/fukui_bars.png
```

## What it demonstrates

| contract feature | used for |
| --- | --- |
| `NAME`, `TITLE`, `REQUIRES` | the pipeline can list and skip the figure without importing its data |
| `ctx.fukui(name)` | finds the per-atom table; the figure never builds a path |
| `ctx.compound(name).formula` | per-compound metadata straight off the contract |
| `ctx.level_of_theory()`, `ctx.provenance` | the same footnote the built-in figures carry |
| `ctx.save(fig, stem)` | writes every format in `figures.formats` at `figures.dpi` |

The full walkthrough, written from scratch, is in
[docs/writing-a-figure.md](../../docs/writing-a-figure.md).
