"""Build 3D starting geometries from SMILES."""
from __future__ import annotations

from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

from .theme import ok, step, warn

PERIODIC = Chem.GetPeriodicTable()


def build_conformer(name: str, smiles: str, cfg: dict) -> Chem.Mol:
    """Embed conformers, force-field optimize, and return the mol with the best one first."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"{name}: RDKit cannot parse the SMILES {smiles!r}")
    mol = Chem.AddHs(mol)

    ps = AllChem.ETKDGv3()
    ps.randomSeed = int(cfg["conformers"]["random_seed"])
    ps.pruneRmsThresh = float(cfg["conformers"]["prune_rms"])
    ps.useSmallRingTorsions = True
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=int(cfg["conformers"]["n_confs"]), params=ps)
    if not cids:
        raise RuntimeError(f"{name}: conformer embedding produced nothing")

    ff = cfg["conformers"]["forcefield"]
    maxit = int(cfg["conformers"]["ff_max_iters"])
    if ff.startswith("MMFF"):
        variant = "MMFF94s" if ff == "MMFF94s" else "MMFF94"
        if not AllChem.MMFFHasAllMoleculeParams(mol):
            warn(f"{name}: MMFF cannot parameterize this molecule; falling back to UFF")
            res = AllChem.UFFOptimizeMoleculeConfs(mol, maxIters=maxit)
        else:
            res = AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=maxit, mmffVariant=variant)
    else:
        res = AllChem.UFFOptimizeMoleculeConfs(mol, maxIters=maxit)

    energies = [(e, cid) for cid, (conv, e) in zip(cids, res)]
    n_unconv = sum(1 for conv, _ in res if conv != 0)
    if n_unconv:
        warn(f"{name}: {n_unconv}/{len(res)} conformers did not converge at the force-field level")
    energies.sort()
    best_e, best_cid = energies[0]
    step(f"{name}: {len(cids)} conformers, best force-field energy {best_e:.3f} kcal/mol")

    out = Chem.Mol(mol)
    out.RemoveAllConformers()
    out.AddConformer(mol.GetConformer(best_cid), assignId=True)
    out.SetProp("_Name", name)
    out.SetProp("ff_energy", f"{best_e:.6f}")
    out.SetProp("n_conformers", str(len(cids)))
    return out


def write_xyz(mol: Chem.Mol, path: Path, comment: str = "") -> Path:
    conf = mol.GetConformer()
    lines = [str(mol.GetNumAtoms()), comment]
    for i, atom in enumerate(mol.GetAtoms()):
        p = conf.GetAtomPosition(i)
        lines.append(f"{atom.GetSymbol():<3s} {p.x:>14.8f} {p.y:>14.8f} {p.z:>14.8f}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")
    return path


def read_xyz(path: Path) -> tuple[list[str], list[tuple[float, float, float]]]:
    lines = Path(path).read_text().splitlines()
    n = int(lines[0].split()[0])
    syms, xyz = [], []
    for ln in lines[2:2 + n]:
        f = ln.split()
        syms.append(f[0])
        xyz.append((float(f[1]), float(f[2]), float(f[3])))
    return syms, xyz


def n_electrons(mol: Chem.Mol, charge: int = 0) -> int:
    return sum(a.GetAtomicNum() for a in mol.GetAtoms()) - charge


def homo_index(mol: Chem.Mol, charge: int = 0, multiplicity: int = 1) -> int:
    """0-based index of the highest occupied MO for a closed-shell system."""
    ne = n_electrons(mol, charge)
    if multiplicity != 1 or ne % 2:
        raise ValueError("homo_index is defined here for closed-shell systems only")
    return ne // 2 - 1


def build_all(cfg: dict, specs: dict) -> dict:
    """Write one starting geometry per molecule to <work>/<name>/<name>_start.xyz.

    Three routes, in order of preference:
      * an .xyz the user supplied            - copied verbatim
      * 3D coordinates already in the SDF    - used as-is
      * anything else                        - ETKDG conformer search
    """
    work = Path(cfg["execution"]["work_dir"])
    use_sdf_geom = bool(cfg["molecules"]["use_sdf_geometry"])

    for name, spec in specs.items():
        outdir = work / name
        outdir.mkdir(parents=True, exist_ok=True)

        if spec.source == "xyz":
            src = Path(spec.source_path)
            syms, _ = read_xyz(src)
            (outdir / f"{name}_start.xyz").write_text(src.read_text())
            step(f"{name}: using supplied geometry {src.name} ({len(syms)} atoms)")
            if spec.mol.GetNumAtoms() != len(syms):
                warn(f"{name}: the SMILES has {spec.mol.GetNumAtoms()} atoms but the xyz has "
                     f"{len(syms)}; atom-indexed outputs (Fukui, metabolism maps) would be "
                     "misaligned")
            continue

        if spec.has_3d and use_sdf_geom:
            write_xyz(spec.mol, outdir / f"{name}_start.xyz",
                      comment=f"{name} charge {spec.charge:+d} mult {spec.multiplicity} "
                              f"(3D geometry from {Path(spec.source_path).name})")
            Chem.MolToMolFile(spec.mol, str(outdir / f"{name}.mol"))
            step(f"{name}: 3D geometry taken from the input file ({spec.n_atoms} atoms)")
            continue

        mol = build_conformer(name, spec.smiles or Chem.MolToSmiles(Chem.RemoveHs(spec.mol)), cfg)
        if mol.GetNumAtoms() != spec.mol.GetNumAtoms():
            warn(f"{name}: conformer has {mol.GetNumAtoms()} atoms vs {spec.mol.GetNumAtoms()} "
                 "in the source - using the conformer, atom indices follow it")
        spec.mol = mol
        write_xyz(mol, outdir / f"{name}_start.xyz",
                  comment=f"{name} charge {spec.charge:+d} mult {spec.multiplicity} "
                          f"ETKDGv3/{cfg['conformers']['forcefield']}")
        Chem.MolToMolFile(mol, str(outdir / f"{name}.mol"))

    ok(f"built {len(specs)} starting geometries in {work}")
    return specs
