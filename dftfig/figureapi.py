"""The figure plugin contract: what a figure module is handed, and what it must expose.

A figure in dftfig is an ordinary Python module. The pipeline finds it, asks it what it
needs, and - if the run can supply that - calls its ``build`` function with one argument:
a :class:`FigureContext`. Everything the figure is allowed to read comes from that
object. A figure never constructs a path, never opens a cube file by name, never reaches
into ``cfg["execution"]["work_dir"]``, and never imports a private helper out of another
figure. Those were all real habits of the four built-ins before this module existed, and
every one of them is a reason a third party could not write a figure against dftfig.

WHAT A FIGURE MODULE MUST EXPOSE
--------------------------------
``NAME``      short identifier, the string a user puts in ``figures.which``.
``TITLE``     one human-readable line, shown by ``dft params figures.which``.
``REQUIRES``  tuple of requirement names (see :data:`REQUIREMENTS`). The pipeline checks
              these BEFORE calling ``build`` and skips the figure with a specific message
              when one is unmet, instead of letting it crash halfway through a run.
``build(ctx)`` the function that draws it. Returns the list of files it wrote.

The minimum is ``NAME`` and ``build``; ``TITLE`` falls back to the module docstring's
first line and ``REQUIRES`` to the empty tuple.

VERSIONING
----------
:data:`API_VERSION` is the contract's version. It is bumped when an existing plugin could
stop working: a method removed or renamed, an argument's meaning changed, a field
dropped from :class:`Compound`, or ``build``'s signature altered. Adding a new accessor
method, a new requirement name or a new optional keyword argument does NOT bump it -
those cannot break code that does not use them. A plugin that needs a guarantee can read
``ctx.api_version``; the pipeline refuses to load a plugin that declares an
``API_VERSION`` higher than its own, because that plugin is asking for something this
dftfig does not have.

DISCOVERY
---------
Two sources, in this order (the later wins on a name collision, which is what lets a
user override a built-in with their own file):

1. the ``dftfig.figures`` entry-point group - dftfig's own four are registered there, and
   so is any separately pip-installed package that declares the same group;
2. plain ``.py`` files in the project's plugin directory (``figures.plugin_dir``,
   ``./dft_figures`` by default) - the case that needs no packaging at all.
"""
from __future__ import annotations

import importlib
import importlib.util
import sys
from dataclasses import dataclass, field
from pathlib import Path
from collections.abc import Callable
from typing import Any

import pandas as pd

#: Bumped only on a change that can break an existing plugin. See the module docstring.
API_VERSION = 1

#: Entry-point group third-party packages register their figures under.
ENTRY_POINT_GROUP = "dftfig.figures"

#: dftfig's own four, imported directly as well as declared as entry points.
#:
#: The direct import is not redundant. `importlib.metadata.entry_points` reads installed
#: distribution metadata, and a source checkout run through `python -m dftfig` - which is
#: how the bundled launcher and the test suite run it - has none. Without this the
#: built-ins would vanish from a checkout while continuing to work from a wheel, which is
#: exactly the kind of difference that only shows up in somebody else's CI.
BUILTIN_MODULES = ("frontier_orbitals", "speciation", "lipophilicity",
                   "mpo_metabolism")

# The built-ins were once named after the figure numbers of the paper this
# pipeline was first written to reproduce. Those numbers mean nothing to anyone
# else, so the names now describe the content. A config written against the old
# names still works: the alias is resolved here and migrated by config.repair,
# which tells the user what it changed.
LEGACY_NAMES = {"fig2": "frontier_orbitals", "fig12": "speciation",
                "fig13": "lipophilicity", "fig14": "mpo_metabolism"}

REQ_RENDERING = "rendering"
REQ_CUBES = "cubes"
REQ_DESCRIPTORS = "descriptors"
REQ_PROPERTIES = "properties"
REQ_FUKUI = "fukui"

#: Every requirement a figure may declare, with the sentence shown when it is unmet.
#: The text is the whole point: "skipped" on its own sends the user looking in the wrong
#: place, which is what happened when a missing PyMOL surfaced as a traceback from deep
#: inside the renderer.
REQUIREMENTS: dict[str, str] = {
    REQ_RENDERING: "no 3D rendering engine is available (render.engine is 'none', or "
                   "PyMOL could not be found - install it with "
                   "`conda install -c conda-forge pymol-open-source`, or set "
                   "render.pymol_path)",
    REQ_CUBES: "cubes.enabled is false, so no volumetric data was produced",
    REQ_DESCRIPTORS: "the DFT descriptor table is empty - run the descriptors stage first",
    REQ_PROPERTIES: "the physicochemical property table is empty",
    REQ_FUKUI: "fukui.enabled is false, so no condensed Fukui indices were computed",
}

#: The 3D panels the context can produce. Anything else raises, rather than silently
#: returning None and leaving an empty panel the author then has to explain.
IMAGE_KINDS = ("geometry", "homo", "lumo", "mep", "lipophilicity")

#: Volumetric fields written by the property single point, by the name a figure asks for.
CUBE_KINDS = {"homo": "homo.cube", "lumo": "lumo.cube",
              "density": "eldens.cube", "potential": "esp.cube"}


class FigureError(RuntimeError):
    """A figure module is malformed, or asks for something this dftfig cannot give."""


# ------------------------------------------------------------------- the data contract
@dataclass(frozen=True)
class Compound:
    """One molecule of the run, as a figure sees it.

    ``mol`` is the in-memory RDKit molecule with explicit hydrogens, in the atom order
    that every per-atom table in the run is indexed by: the Fukui CSV, the Crippen
    contributions and the geometry ORCA reports. A figure that paints per-atom values
    onto a depiction may rely on that order.
    """

    name: str
    mol: Any                     # rdkit.Chem.Mol, typed loosely to keep this module light
    charge: int
    multiplicity: int
    formula: str
    run_dir: Path                # where this compound's calculations live


@dataclass
class FigureContext:
    """Everything a figure is given, and the only thing it may use.

    The accessor methods exist because the built-ins used to compute each of them by hand
    from ``cfg["execution"]["work_dir"]``. Keeping that arithmetic here means the layout
    of a run directory is dftfig's business, not a plugin author's, and it can change
    without breaking a single figure.
    """

    config: dict
    compounds: dict[str, Compound]
    descriptors: pd.DataFrame
    properties: pd.DataFrame
    output_dir: Path
    provenance: dict
    api_version: int = API_VERSION

    #: Prepared lipophilicity surfaces, cached so a figure that asks for the scale and
    #: then for the image does not build the same PDB twice.
    _mlp_cache: dict = field(default_factory=dict, repr=False)
    _geometry_written: set = field(default_factory=set, repr=False)

    # ---------------------------------------------------------------- the molecule set
    @property
    def names(self) -> list[str]:
        """Compound names in run order. This is the column order of a multi-panel row."""
        return list(self.compounds)

    @property
    def mols(self) -> dict:
        """``{name: RDKit Mol}``, in the same order as :attr:`names`."""
        return {n: c.mol for n, c in self.compounds.items()}

    def compound(self, name: str) -> Compound | None:
        return self.compounds.get(name)

    # --------------------------------------------------------------- run-wide settings
    @property
    def reference_solvent(self) -> str:
        """The medium the figures are built from - ``figures.reference_solvent``."""
        return self.config["figures"]["reference_solvent"]

    def dielectric(self, medium: str | None = None) -> float:
        """Relative permittivity of a medium; NaN when dftfig has no value for it."""
        from .orca_io import DIELECTRIC
        return DIELECTRIC.get(medium or self.reference_solvent, float("nan"))

    def level_of_theory(self) -> str:
        """The method line a figure prints as its footnote."""
        t = self.config["theory"]
        disp = "" if t["dispersion"] == "none" else f"-{t['dispersion']}"
        opt = f"{t['functional']}{disp}/{t['basis_opt']}"
        solv = self.config["solvation"]["model"]
        if self.config["property"]["enabled"]:
            return (f"{opt} opt // {t['functional']}{disp}/"
                    f"{self.config['property']['basis_sp']} property, {solv}")
        return f"{opt}, {solv}"

    def param(self, key: str, default: Any = None) -> Any:
        """One configuration value by dotted key, e.g. ``ctx.param("figures.dpi")``."""
        section, _, name = key.partition(".")
        return self.config.get(section, {}).get(name, default)

    # ---------------------------------------------------------------------- the numbers
    def descriptor_row(self, name: str, medium: str | None = None):
        """The descriptors.csv row for one compound in one medium, or None.

        Values are Python scalars, not numpy ones. That is a promise, not an accident:
        ``df.iloc[0]["is_minimum"]`` hands back ``numpy.bool_``, for which the natural
        test ``value is True`` is False, and a plugin author has no reason to expect it.
        """
        if self.descriptors is None or self.descriptors.empty:
            return None
        medium = medium or self.reference_solvent
        sub = self.descriptors[(self.descriptors["molecule"] == name)
                               & (self.descriptors["medium"] == medium)]
        return None if sub.empty else _native_row(sub)

    def descriptor(self, name: str, column: str, medium: str | None = None,
                   default: Any = None) -> Any:
        """One DFT descriptor for one compound in one medium."""
        row = self.descriptor_row(name, medium)
        if row is None or column not in row.index:
            return default
        value = row[column]
        return default if value is None or _isna(value) else value

    def property_row(self, name: str):
        """The resolved physicochemical row for one compound, or None.

        Python scalars, for the same reason as :meth:`descriptor_row`.
        """
        if self.properties is None or self.properties.empty:
            return None
        sub = self.properties[self.properties["molecule"] == name]
        return None if sub.empty else _native_row(sub)

    def property(self, name: str, column: str, default: Any = None) -> Any:
        """One physicochemical value - measured where supplied, predicted otherwise."""
        row = self.property_row(name)
        if row is None or column not in row.index:
            return default
        value = row[column]
        return default if value is None or _isna(value) else value

    def fukui(self, name: str, medium: str | None = None) -> pd.DataFrame | None:
        """Condensed Fukui indices for one compound, one row per atom, or None.

        Columns include ``atom_index_0based``, ``element``, ``f_plus``, ``f_minus``,
        ``f_zero`` and ``dual_descriptor``. Atom indices are into :attr:`Compound.mol`,
        hydrogens included.
        """
        medium = medium or self.reference_solvent
        path = self.run_dir(name) / f"fukui_{name}_{medium}.csv"
        if not path.exists():
            return None
        return pd.read_csv(path)

    def som_table(self, name: str) -> pd.DataFrame | None:
        """An externally predicted site-of-metabolism export for one compound, or None."""
        path = Path(self.config["figures"]["som_dir"]) / f"{name}.csv"
        if not path.exists():
            return None
        return pd.read_csv(path)

    # ------------------------------------------------------------------------ the files
    def run_dir(self, name: str, medium: str | None = None) -> Path:
        """Where one compound's calculations live. ``medium=None`` gives the parent."""
        c = self.compounds.get(name)
        base = c.run_dir if c else Path(self.config["execution"]["work_dir"]) / name
        return base if medium is None else base / medium

    def cube(self, name: str, kind: str, medium: str | None = None) -> Path | None:
        """Locate one volumetric file; None when it was never produced.

        ``kind`` is one of ``homo``, ``lumo``, ``density``, ``potential`` - the physical
        quantity, not the file name, so the on-disk naming stays dftfig's business.
        """
        if kind not in CUBE_KINDS:
            raise FigureError(f"unknown cube kind {kind!r}; expected one of "
                              f"{', '.join(sorted(CUBE_KINDS))}")
        path = self.run_dir(name, medium or self.reference_solvent) / CUBE_KINDS[kind]
        return path if path.exists() else None

    def surface_esp_stats(self, name: str, medium: str | None = None) -> dict:
        """Electrostatic potential statistics on the molecular surface; ``{}`` if absent.

        Keys: ``v_min``, ``v_max``, ``v_lo``, ``v_hi`` (2nd/98th percentiles - the honest
        colour limits), ``v_mean``, ``n_surface_points``.
        """
        from . import cubes
        return cubes.surface_esp_stats(self.config, name, medium or self.reference_solvent)

    def lipophilicity_range(self, name: str) -> tuple[float, float] | None:
        """Min/max per-atom Crippen contribution, or None if the surface cannot be built.

        Returned separately from the image so a figure can fix ONE colour scale across
        every panel before any of them is rendered - without that each panel would be
        scaled to itself and the colours would not mean the same thing twice.
        """
        pdb, lo, hi = self._lipophilicity(name)
        return None if pdb is None else (lo, hi)

    def image(self, name: str, kind: str, medium: str | None = None,
              vrange: tuple[float, float] | None = None) -> Path | None:
        """A rendered 3D panel as a PNG, or None when it could not be produced.

        ``kind`` is one of :data:`IMAGE_KINDS`. ``vrange`` sets the colour limits for
        ``mep`` and ``lipophilicity`` and is ignored by the others. The context also does
        the two chores every caller used to repeat: writing out the geometry the renderer
        reads, and trimming the white border when the background is transparent.
        """
        if kind not in IMAGE_KINDS:
            raise FigureError(f"unknown image kind {kind!r}; expected one of "
                              f"{', '.join(IMAGE_KINDS)}")
        from . import render
        medium = medium or self.reference_solvent

        if kind == "lipophilicity":
            pdb, lo, hi = self._lipophilicity(name)
            if pdb is None:
                return None
            lo, hi = vrange if vrange is not None else (lo, hi)
            path = render.render_mlp(self.config, name, pdb, lo, hi)
        else:
            self._ensure_geometry(name, medium)
            if kind == "geometry":
                path = render.render_geometry(self.config, name, medium)
            elif kind in ("homo", "lumo"):
                path = render.render_orbital(self.config, name, medium, kind)
            else:
                if vrange is None:
                    r = float(self.config["render"]["esp_range"])
                    vrange = (-r, r)
                path = render.render_mep(self.config, name, medium, vrange)

        if path and Path(path).exists() and self.config["render"]["background"] == "transparent":
            from .figstyle import trim_white
            trim_white(path)
        return path

    def _ensure_geometry(self, name: str, medium: str) -> None:
        if (name, medium) not in self._geometry_written:
            from . import render
            render.write_sp_geometry(self.config, name, medium)
            self._geometry_written.add((name, medium))

    def _lipophilicity(self, name: str):
        if name not in self._mlp_cache:
            c = self.compounds.get(name)
            if c is None:
                self._mlp_cache[name] = (None, 0.0, 0.0)
            else:
                from . import mlp as mlpmod
                self._mlp_cache[name] = mlpmod.prepare_mlp(self.config, name, c.mol)
        return self._mlp_cache[name]

    # ------------------------------------------------------------------------- drawing
    def apply_style(self) -> None:
        """Install dftfig's matplotlib style at the configured resolution."""
        from .figstyle import apply_style
        apply_style(int(self.config["figures"]["dpi"]))

    def esp_colormap(self):
        """Matplotlib colormap matching the 3D MEP panels, for the shared colour bar."""
        from .figstyle import esp_colormap
        return esp_colormap(self.config["render"]["esp_palette"])

    def mlp_colormap(self):
        """Matplotlib colormap matching the 3D lipophilicity panels."""
        from .figstyle import mlp_colormap
        return mlp_colormap(self.config["figures"]["mlp_palette"])

    def save(self, fig, stem: str) -> list[Path]:
        """Write a finished figure in every configured format. Returns the paths."""
        import matplotlib.pyplot as plt

        from .theme import ok
        self.output_dir.mkdir(parents=True, exist_ok=True)
        paths = []
        for fmt in self.config["figures"]["formats"]:
            p = self.output_dir / f"{stem}.{fmt}"
            fig.savefig(p, dpi=int(self.config["figures"]["dpi"]))
            paths.append(p)
        plt.close(fig)
        ok(f"{stem} -> {', '.join(str(p) for p in paths)}")
        return paths

    # -------------------------------------------------------------------- requirements
    def unmet(self, requires) -> list[str]:
        """Which of ``requires`` this run cannot satisfy, as sentences for the user."""
        out = []
        for req in requires:
            if req not in REQUIREMENTS:
                out.append(f"declares an unknown requirement {req!r}")
            elif not self.provides(req):
                out.append(REQUIREMENTS[req])
        return out

    def provides(self, requirement: str) -> bool:
        """Whether this run can satisfy one requirement."""
        if requirement == REQ_RENDERING:
            if self.config["render"]["engine"] == "none":
                return False
            from . import render
            return render.find_pymol(self.config) is not None
        if requirement == REQ_CUBES:
            return bool(self.config["cubes"]["enabled"])
        if requirement == REQ_DESCRIPTORS:
            return self.descriptors is not None and not self.descriptors.empty
        if requirement == REQ_PROPERTIES:
            return self.properties is not None and not self.properties.empty
        if requirement == REQ_FUKUI:
            return bool(self.config["fukui"]["enabled"])
        return False


def _native_row(sub: pd.DataFrame):
    """First row of ``sub`` with Python scalars rather than numpy ones."""
    return sub.astype(object).iloc[0]


def _isna(value) -> bool:
    try:
        return bool(pd.isna(value))
    except (TypeError, ValueError):
        return False


def build_context(cfg: dict, specs: dict, descriptors, properties) -> FigureContext:
    """Assemble the context the pipeline hands to every figure of one run."""
    from . import provenance
    work = Path(cfg["execution"]["work_dir"])
    compounds = {
        name: Compound(name=name, mol=s.mol, charge=s.charge,
                       multiplicity=s.multiplicity, formula=s.formula,
                       run_dir=work / name)
        for name, s in specs.items()
    }
    return FigureContext(
        config=cfg, compounds=compounds,
        descriptors=pd.DataFrame() if descriptors is None else descriptors,
        properties=pd.DataFrame() if properties is None else properties,
        output_dir=work / "figures", provenance=provenance.stamp(cfg))


# --------------------------------------------------------------------------- discovery
@dataclass(frozen=True)
class FigureSpec:
    """What the pipeline knows about one figure before it draws anything."""

    name: str
    title: str
    requires: tuple
    build: Callable[[FigureContext], list]
    origin: str                  # where it came from, for `dft params figures.which`
    module: str = ""

    @property
    def builtin(self) -> bool:
        return self.module.startswith("dftfig.figures.")


def describe(module, origin: str = "") -> FigureSpec:
    """Read a module's declaration, checking it before the pipeline relies on it."""
    build = getattr(module, "build", None)
    if not callable(build):
        raise FigureError(f"{module.__name__} has no build() function")
    declared = int(getattr(module, "API_VERSION", API_VERSION))
    if declared > API_VERSION:
        raise FigureError(f"{module.__name__} needs figure API v{declared}, but this "
                          f"dftfig provides v{API_VERSION}")
    name = getattr(module, "NAME", None) or module.__name__.rsplit(".", 1)[-1]
    doc = (module.__doc__ or "").strip().splitlines()
    title = getattr(module, "TITLE", None) or (doc[0] if doc else name)
    requires = tuple(getattr(module, "REQUIRES", ()))
    unknown = [r for r in requires if r not in REQUIREMENTS]
    if unknown:
        raise FigureError(f"{module.__name__} declares unknown requirement(s) {unknown}; "
                          f"known names are {sorted(REQUIREMENTS)}")
    return FigureSpec(name=name, title=title, requires=requires, build=build,
                      origin=origin, module=module.__name__)


def plugin_dir(cfg: dict | None = None) -> Path:
    """The project-local directory scanned for plain ``.py`` figures."""
    raw = (cfg or {}).get("figures", {}).get("plugin_dir") or "dft_figures"
    return Path(raw).expanduser()


def _from_builtins() -> dict[str, FigureSpec]:
    from .theme import warn
    found: dict[str, FigureSpec] = {}
    for stem in BUILTIN_MODULES:
        try:
            module = importlib.import_module(f".figures.{stem}", package=__package__)
            spec = describe(module, "built in")
        except Exception as exc:
            warn(f"built-in figure {stem!r} could not be loaded: {exc}")
            continue
        found[spec.name] = spec
    return found


def _from_entry_points() -> dict[str, FigureSpec]:
    from importlib.metadata import entry_points

    from .theme import warn
    found: dict[str, FigureSpec] = {}
    for ep in entry_points(group=ENTRY_POINT_GROUP):
        try:
            module = ep.load()
            origin = ("built in" if str(getattr(module, "__name__", "")).startswith(
                "dftfig.figures.") else f"package '{ep.value.split('.')[0]}'")
            spec = describe(module, origin)
        except Exception as exc:                 # a broken plugin must not kill a run
            warn(f"figure plugin {ep.name!r} could not be loaded: {exc}")
            continue
        found[ep.name or spec.name] = spec
    return found


def _from_directory(directory: Path) -> dict[str, FigureSpec]:
    """Import every ``.py`` in a directory as a figure module.

    Loaded by file location rather than by putting the directory on ``sys.path``: two
    projects may each hold a ``dft_figures/panel.py``, and a path insertion would make
    whichever ran first shadow the other for the rest of the process.
    """
    found: dict[str, FigureSpec] = {}
    if not directory.is_dir():
        return found
    from .theme import warn
    for path in sorted(directory.glob("*.py")):
        if path.name.startswith("_"):
            continue
        mod_name = f"dftfig_plugin_{path.stem}"
        try:
            loader_spec = importlib.util.spec_from_file_location(mod_name, path)
            module = importlib.util.module_from_spec(loader_spec)
            sys.modules[mod_name] = module
            loader_spec.loader.exec_module(module)
            spec = describe(module, f"{directory}/{path.name}")
        except Exception as exc:
            warn(f"figure plugin {path} could not be loaded: {exc}")
            continue
        found[spec.name] = spec
    return found


def available(cfg: dict | None = None) -> dict[str, FigureSpec]:
    """Every figure this installation can draw, keyed by the name used in ``figures.which``.

    A project-local file wins over an installed package of the same name, which is what
    makes overriding a built-in possible without uninstalling anything.
    """
    found = _from_builtins()
    found.update(_from_entry_points())
    found.update(_from_directory(plugin_dir(cfg)))
    return found


def load(name: str, cfg: dict | None = None) -> FigureSpec:
    """Resolve a figure name to its spec, accepting the superseded names."""
    name = LEGACY_NAMES.get(name, name)
    found = available(cfg)
    if name not in found:
        raise FigureError(
            f"no figure named {name!r}. Available: " + ", ".join(sorted(found))
        )
    return found[name]
