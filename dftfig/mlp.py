"""Molecular lipophilicity potential and physicochemical properties.

The MLP of Eq. 5 in the manuscript, MLP = sum_i F_i f(d_ik), is a distance-weighted sum
of per-fragment lipophilicity contributions. Here the contributions F_i are RDKit's
Crippen atomic logP increments and the distance weighting is performed by the renderer:
each atom's contribution is written into the B-factor column of a PDB and PyMOL
interpolates it across the solvent-excluded surface. This is the standard practical
construction of an MLP surface and keeps the whole figure reproducible from the
optimized geometry alone.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import Crippen, rdMolDescriptors

from . import parsers
from .theme import warn


def crippen_atom_contributions(mol: Chem.Mol) -> np.ndarray:
    """Per-atom Crippen logP increments, in RDKit atom order (hydrogens included)."""
    contribs = rdMolDescriptors._CalcCrippenContribs(mol)
    return np.array([c[0] for c in contribs])


def tpsa_atom_contributions(mol: Chem.Mol) -> np.ndarray:
    """Per-atom TPSA contributions (zero for carbon and hydrogen)."""
    return np.array(rdMolDescriptors._CalcTPSAContribs(mol))


def rdkit_properties(mol: Chem.Mol) -> dict:
    """Everything the structure alone determines - no external predictor needed."""
    heavy = Chem.RemoveHs(mol)
    return {
        "formula": rdMolDescriptors.CalcMolFormula(heavy),
        "MW": rdMolDescriptors.CalcExactMolWt(heavy),
        "TPSA": rdMolDescriptors.CalcTPSA(heavy),
        "clogP_crippen": Crippen.MolLogP(heavy),
        "MR_crippen": Crippen.MolMR(heavy),
        "HBD": rdMolDescriptors.CalcNumHBD(heavy),
        "HBA": rdMolDescriptors.CalcNumHBA(heavy),
        "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(heavy),
        "aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(heavy),
        "heavy_atoms": heavy.GetNumHeavyAtoms(),
    }


ELEMENT_ORDER: dict[str, int] = {}


def write_pdb_with_bfactors(name: str, symbols: list[str], coords, values,
                            path: Path) -> Path:
    """Write a PDB whose B-factor column carries one scalar per atom."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"REMARK   dftfig  {name}  B-factor column = per-atom lipophilicity contribution"]
    counts: dict[str, int] = {}
    for i, (sym, (x, y, z), v) in enumerate(zip(symbols, coords, values), 1):
        counts[sym] = counts.get(sym, 0) + 1
        aname = f"{sym}{counts[sym]}"[:4]
        lines.append(
            f"HETATM{i:5d} {aname:<4s}LIG A   1    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00{v:6.2f}          {sym:>2s}"
        )
    lines.append("END")
    path.write_text("\n".join(lines) + "\n")
    return path


def prepare_mlp(cfg: dict, name: str, mol: Chem.Mol) -> tuple[Path | None, float, float]:
    """Build the B-factor PDB for one molecule from its optimized geometry.

    Returns (pdb_path, min_contribution, max_contribution).
    """
    solvent = cfg["figures"]["reference_solvent"]
    work = Path(cfg["execution"]["work_dir"]) / name / solvent
    src = work / f"{name}_{solvent}_sp.out"
    if not src.exists():
        src = work / f"{name}_{solvent}_opt.out"
    if not src.exists():
        warn(f"{name}: no optimized geometry for the lipophilicity surface")
        return None, 0.0, 0.0

    syms, coords = parsers.geometry(src)
    if not syms:
        warn(f"{name}: could not read the optimized geometry")
        return None, 0.0, 0.0

    contribs = crippen_atom_contributions(mol)
    if len(contribs) != len(syms):
        warn(f"{name}: {len(contribs)} RDKit atoms vs {len(syms)} ORCA atoms - "
             "the lipophilicity surface would be misassigned, skipping")
        return None, 0.0, 0.0

    mol_syms = [a.GetSymbol() for a in mol.GetAtoms()]
    if mol_syms != syms:
        # The per-atom Crippen contributions are indexed by the in-memory molecule. If
        # that molecule's atom order is not the order ORCA reports, there is no mapping
        # to hand and painting the surface anyway would put a hydroxyl's lipophilicity on
        # a ring carbon. Skip, and say precisely where the orders diverge so the cause can
        # be found - including whether the geometry actually sent to ORCA already differed
        # from the molecule, which localises the drift to before or after the input was
        # written.
        start = Path(cfg["execution"]["work_dir"]) / name / f"{name}_start.xyz"
        start_syms = []
        if start.exists():
            lines = start.read_text().splitlines()
            try:
                start_syms = [ln.split()[0] for ln in lines[2:2 + int(lines[0].split()[0])]]
            except (ValueError, IndexError):
                start_syms = []
        first = next((i for i, (a, b) in enumerate(zip(mol_syms, syms)) if a != b),
                     min(len(mol_syms), len(syms)))
        where = ""
        if start_syms:
            where = (" (the geometry sent to ORCA matches its output, so the in-memory "
                     "molecule is the one that changed)" if start_syms == syms else
                     " (the geometry sent to ORCA already differed from its output)")
        warn(f"{name}: atom ORDER differs between RDKit and ORCA at position {first} - "
             f"RDKit has {''.join(mol_syms[first:first + 4])}..., ORCA has "
             f"{''.join(syms[first:first + 4])}...{where}. Skipping the lipophilicity "
             f"surface rather than colouring the wrong atoms.")
        return None, 0.0, 0.0

    pdb = write_pdb_with_bfactors(name, syms, coords, contribs, work / f"{name}_mlp.pdb")
    return pdb, float(contribs.min()), float(contribs.max())
