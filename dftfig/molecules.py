"""Molecule input: SMILES, SDF or XYZ, with per-molecule charge and multiplicity.

Charge is a property of each compound, not of the run. A set containing a neutral
nucleoside and a quaternary ammonium cation has two different electron counts, and
handing ORCA a single global charge would silently compute the wrong species for one of
them. Every molecule therefore carries its own charge and multiplicity, detected from the
structure unless you override it.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

from .theme import bad, ok, step, warn


@dataclass
class MoleculeSpec:
    name: str
    mol: Chem.Mol                 # hydrogens explicit
    charge: int
    multiplicity: int
    source: str                   # smiles | sdf | xyz
    source_path: str | None = None
    smiles: str = ""
    charge_origin: str = ""
    multiplicity_origin: str = ""
    has_3d: bool = False
    notes: list = field(default_factory=list)

    @property
    def n_electrons(self) -> int:
        return sum(a.GetAtomicNum() for a in self.mol.GetAtoms()) - self.charge

    @property
    def formula(self) -> str:
        return rdMolDescriptors.CalcMolFormula(self.mol)

    @property
    def n_atoms(self) -> int:
        return self.mol.GetNumAtoms()

    def homo_index(self) -> int:
        """0-based index of the highest occupied MO (closed shell only)."""
        ne = self.n_electrons
        if self.multiplicity != 1 or ne % 2:
            raise ValueError(f"{self.name}: frontier-orbital index is defined here for "
                             f"closed-shell species only (multiplicity {self.multiplicity}, "
                             f"{ne} electrons)")
        return ne // 2 - 1

    def ion(self, delta: int) -> tuple[int, int]:
        """(charge, multiplicity) for the species with `delta` electrons added."""
        charge = self.charge - delta
        ne = sum(a.GetAtomicNum() for a in self.mol.GetAtoms()) - charge
        return charge, lowest_multiplicity(ne)


def lowest_multiplicity(n_electrons: int) -> int:
    """Lowest spin state consistent with the electron count: singlet or doublet."""
    return 1 if n_electrons % 2 == 0 else 2


def _is_3d(mol: Chem.Mol) -> bool:
    if mol.GetNumConformers() == 0:
        return False
    conf = mol.GetConformer()
    if conf.Is3D():
        # A "3D" flag with every z at zero is a flattened depiction, not a geometry.
        zs = [abs(conf.GetAtomPosition(i).z) for i in range(mol.GetNumAtoms())]
        return max(zs) > 1e-3
    return False


def _charge_for(mol: Chem.Mol, cfg: dict, name: str) -> tuple[int, str]:
    mode = cfg["molecules"]["charge_mode"]
    if mode == "fixed":
        return int(cfg["molecules"]["charge"]), "fixed in config"
    q = Chem.GetFormalCharge(mol)
    return q, "detected from formal charges"


def _multiplicity_for(mol: Chem.Mol, charge: int, cfg: dict, name: str) -> tuple[int, str]:
    ne = sum(a.GetAtomicNum() for a in mol.GetAtoms()) - charge
    auto = lowest_multiplicity(ne)
    if cfg["molecules"]["multiplicity_mode"] == "fixed":
        m = int(cfg["molecules"]["multiplicity"])
        if (ne % 2 == 0) != (m % 2 == 1):
            raise ValueError(
                f"{name}: multiplicity {m} is impossible for {ne} electrons. "
                f"An even electron count needs an odd multiplicity (1, 3, 5...) and an odd "
                f"count needs an even one (2, 4...). The structure implies {auto}."
            )
        return m, "fixed in config"
    return auto, f"lowest spin state for {ne} electrons"


def _from_smiles(name: str, smiles: str, cfg: dict) -> MoleculeSpec:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"{name}: RDKit cannot parse the SMILES {smiles!r}")
    mol = Chem.AddHs(mol)
    q, qorig = _charge_for(mol, cfg, name)
    m, morig = _multiplicity_for(mol, q, cfg, name)
    return MoleculeSpec(name=name, mol=mol, charge=q, multiplicity=m, source="smiles",
                        smiles=Chem.MolToSmiles(Chem.RemoveHs(mol)),
                        charge_origin=qorig, multiplicity_origin=morig, has_3d=False)


def _sdf_records(path: Path, cfg: dict):
    """Yield (name, mol) for every record in one SDF file."""
    supplier = Chem.SDMolSupplier(str(path), removeHs=False, sanitize=True)
    stem = path.stem
    records = []
    for i, mol in enumerate(supplier):
        if mol is None:
            warn(f"{path.name}: record {i + 1} could not be parsed - skipped")
            continue
        title = (mol.GetProp("_Name").strip() if mol.HasProp("_Name") else "")
        if cfg["molecules"]["name_from"] == "sdf_title" and title:
            name = title
        else:
            name = stem
        records.append((name, mol))
    if len(records) > 1:
        # Disambiguate several records sharing one name.
        seen: dict[str, int] = {}
        out = []
        for name, mol in records:
            seen[name] = seen.get(name, 0) + 1
            out.append((name if seen[name] == 1 else f"{name}_{seen[name]}", mol))
        return out
    return records


def _from_sdf(path: Path, cfg: dict) -> list[MoleculeSpec]:
    specs = []
    for name, raw in _sdf_records(path, cfg):
        mol = raw
        had3d = _is_3d(mol)
        mol = Chem.AddHs(mol, addCoords=had3d)
        q, qorig = _charge_for(mol, cfg, name)
        m, morig = _multiplicity_for(mol, q, cfg, name)
        try:
            smi = Chem.MolToSmiles(Chem.RemoveHs(Chem.Mol(mol)))
        except Exception:
            smi = ""
        spec = MoleculeSpec(name=name, mol=mol, charge=q, multiplicity=m, source="sdf",
                            source_path=str(path), smiles=smi,
                            charge_origin=qorig, multiplicity_origin=morig, has_3d=had3d)
        if not had3d:
            spec.notes.append("no 3D coordinates in the SDF - a conformer will be generated")
        specs.append(spec)
    return specs


def _sanitize_name(name: str) -> str:
    keep = "".join(ch if (ch.isalnum() or ch in "-_") else "_" for ch in name.strip())
    return keep.strip("_") or "mol"


def load_molecules(cfg: dict) -> dict[str, MoleculeSpec]:
    """Build the molecule set from whichever source the config names."""
    source = cfg["molecules"]["structure_source"]
    specs: list[MoleculeSpec] = []

    if source == "sdf":
        d = Path(cfg["molecules"]["sdf_dir"]).expanduser()
        pattern = cfg["molecules"]["sdf_pattern"]
        if not d.exists():
            raise FileNotFoundError(
                f"molecules.sdf_dir = {d} does not exist. Create it and drop one .sdf per "
                f"compound inside, or run `dft init` to scaffold it."
            )
        files = sorted(d.glob(pattern))
        if not files:
            raise FileNotFoundError(
                f"no files matching {pattern!r} in {d}. Put one .sdf per compound there."
            )
        for f in files:
            specs.extend(_from_sdf(f, cfg))

    elif source == "smiles":
        defs = cfg["molecules"]["definitions"]
        if not defs:
            raise ValueError(
                "molecules.structure_source = 'smiles' but molecules.definitions is empty. "
                "Add NAME: SMILES entries, or switch to 'sdf' and point molecules.sdf_dir "
                "at a folder of .sdf files."
            )
        for name, smi in defs.items():
            specs.append(_from_smiles(name, smi, cfg))

    elif source == "xyz":
        d = Path(cfg["molecules"]["xyz_dir"]).expanduser()
        defs = cfg["molecules"]["definitions"]
        if not defs:
            raise ValueError(
                "structure_source = 'xyz' still needs molecules.definitions: the SMILES "
                "supplies bond orders, hydrogen-bond perception and the Crippen "
                "lipophilicity increments that the XYZ cannot."
            )
        for name, smi in defs.items():
            spec = _from_smiles(name, smi, cfg)
            spec.source = "xyz"
            spec.source_path = str(d / f"{name}.xyz")
            if not Path(spec.source_path).exists():
                raise FileNotFoundError(f"{name}: {spec.source_path} is missing")
            specs.append(spec)
    else:
        raise ValueError(f"unknown molecules.structure_source: {source!r}")

    # --- de-duplicate names, then apply the run limit --------------------------------
    seen: dict[str, int] = {}
    for s in specs:
        base = _sanitize_name(s.name)
        seen[base] = seen.get(base, 0) + 1
        s.name = base if seen[base] == 1 else f"{base}_{seen[base]}"

    # --- an explicit subset, if one was named ----------------------------------------
    only = cfg["molecules"].get("only") or ""
    wanted = selected_names(cfg)
    if wanted:
        have = {s.name for s in specs}
        missing = [w for w in wanted if w not in have]
        if missing:
            warn(f"molecules.only names {', '.join(missing)}, which "
                 f"{'is' if len(missing) == 1 else 'are'} not among the structures found "
                 f"({', '.join(sorted(have))})")
        specs = [s for s in specs if s.name in set(wanted)]
        if not specs:
            raise ValueError(
                f"molecules.only = {only!r} selected nothing. Clear it to run every "
                f"structure, or use one of: {', '.join(sorted(have))}")

    limit = int(cfg["molecules"]["limit"])
    if limit and limit < len(specs):
        warn(f"{len(specs)} molecule(s) found; molecules.limit = {limit}, so only the first "
             f"{limit} will run: {', '.join(s.name for s in specs[:limit])}")
        specs = specs[:limit]
    elif limit > len(specs):
        warn(f"molecules.limit = {limit} but only {len(specs)} molecule(s) were found")

    return {s.name: s for s in specs}


_REPORTED: set = set()


def report(specs: dict[str, MoleculeSpec], force: bool = False) -> None:
    """Print what was detected, so a wrong charge is visible before hours of compute.

    Each stage re-reads the structures, so without this guard the same table is printed
    three or four times in one run and the thing it exists to make visible - an
    unexpected charge - is lost in the repetition.
    """
    from rich.table import Table
    from .theme import PALETTE, console

    key = tuple(sorted((n, s.charge, s.multiplicity) for n, s in specs.items()))
    if not force and key in _REPORTED:
        return
    _REPORTED.add(key)

    t = Table(header_style="brand", box=None, padding=(0, 2))
    for h in ["molecule", "formula", "atoms", "charge", "mult", "electrons", "geometry", "source"]:
        t.add_column(h, justify="right" if h in ("atoms", "charge", "mult", "electrons") else "left")
    charges = set()
    for s in specs.values():
        charges.add(s.charge)
        qstyle = "value" if s.charge == 0 else "warn"
        t.add_row(s.name, s.formula, str(s.n_atoms),
                  f"[{qstyle}]{s.charge:+d}[/{qstyle}]" if s.charge else "0",
                  str(s.multiplicity), str(s.n_electrons),
                  "from file" if s.has_3d else "generated",
                  Path(s.source_path).name if s.source_path else s.source)
    console.print(t)

    origins = {s.charge_origin for s in specs.values()}
    console.print(f"  [muted]charge: {', '.join(sorted(origins))}[/muted]")
    if len(charges) > 1:
        shown = ", ".join(f"{c:+d}" if c else "0" for c in sorted(charges))
        console.print(f"  [ok]mixed charges in this run:[/ok] [value]{shown}[/value]"
                      f"  [muted]- each molecule keeps its own[/muted]")
    for s in specs.values():
        for n in s.notes:
            console.print(f"  [muted]{s.name}: {n}[/muted]")


def write_manifest(cfg: dict, specs: dict[str, MoleculeSpec]) -> Path:
    import pandas as pd
    work = Path(cfg["execution"]["work_dir"])
    work.mkdir(parents=True, exist_ok=True)
    rows = [{
        "molecule": s.name, "formula": s.formula, "n_atoms": s.n_atoms,
        "charge": s.charge, "charge_origin": s.charge_origin,
        "multiplicity": s.multiplicity, "multiplicity_origin": s.multiplicity_origin,
        "n_electrons": s.n_electrons, "geometry_from_file": s.has_3d,
        "source": s.source, "source_path": s.source_path or "", "smiles": s.smiles,
    } for s in specs.values()]
    path = work / "molecules.csv"
    pd.DataFrame(rows).to_csv(path, index=False)
    return path


def names(cfg: dict) -> list:
    """Molecule names for this run, in order.

    Prefers the manifest written by the current run, so downstream stages see exactly the
    set that was calculated (including any molecules.limit), and falls back to resolving
    the sources again.
    """
    manifest = Path(cfg["execution"]["work_dir"]) / "molecules.csv"
    if manifest.exists():
        import pandas as pd
        try:
            return [str(n) for n in pd.read_csv(manifest)["molecule"].tolist()]
        except Exception:
            pass
    try:
        return list(load_molecules(cfg))
    except Exception:
        return list(cfg["molecules"]["definitions"])


def selected_names(cfg: dict) -> list:
    """molecules.only as a list of names, empty when every structure runs."""
    only = cfg["molecules"].get("only") or ""
    if isinstance(only, str):
        return [x.strip() for x in only.split(",") if x.strip()]
    return [str(x).strip() for x in only]


def count(cfg: dict) -> int:
    """How many compounds this configuration will actually run."""
    source = cfg["molecules"]["structure_source"]
    limit = int(cfg["molecules"]["limit"])
    wanted = selected_names(cfg)
    if wanted:
        # An explicit subset decides the count on its own; the cost estimate has to match
        # what the run will do, not how many files are in the folder.
        n = len(wanted)
    elif source == "sdf":
        d = Path(cfg["molecules"]["sdf_dir"]).expanduser()
        n = len(sorted(d.glob(cfg["molecules"]["sdf_pattern"]))) if d.exists() else 0
    else:
        n = len(cfg["molecules"]["definitions"])
    return min(n, limit) if limit else n
