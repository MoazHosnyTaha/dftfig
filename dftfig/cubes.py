"""Cube-file handling and quantum-mechanical electrostatic potential evaluation.

ORCA's %plots block writes molecular-orbital and electron-density cubes directly during
the property single point. It has no entry for the electrostatic potential, so the MEP is
built here in three steps:

  1. read the grid definition out of the electron-density cube,
  2. hand every grid point to orca_vpot, which evaluates V(r) from the SCF density,
  3. write the potentials back out as a cube sharing the density cube's grid,

so that a rendering program can paint the second onto an isosurface of the first.
"""
from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import molecules as _M
from .theme import ok, step, warn

BOHR_PER_ANGSTROM = 1.8897261254578281


@dataclass
class Cube:
    comment1: str
    comment2: str
    origin: np.ndarray           # (3,) in Bohr
    shape: tuple                 # (nx, ny, nz)
    axes: np.ndarray             # (3, 3) voxel vectors in Bohr
    atom_z: np.ndarray           # (natoms,)
    atom_q: np.ndarray           # (natoms,)
    atom_xyz: np.ndarray         # (natoms, 3) in Bohr
    data: np.ndarray             # (nx, ny, nz)

    @property
    def natoms(self) -> int:
        return len(self.atom_z)

    def grid_points(self) -> np.ndarray:
        """(N, 3) array of every grid point in Bohr, in cube storage order."""
        nx, ny, nz = self.shape
        i = np.arange(nx)[:, None, None]
        j = np.arange(ny)[None, :, None]
        k = np.arange(nz)[None, None, :]
        pts = (self.origin[None, None, None, :]
               + i[..., None] * self.axes[0]
               + j[..., None] * self.axes[1]
               + k[..., None] * self.axes[2])
        return pts.reshape(-1, 3)


def read_cube(path: str | Path) -> Cube:
    lines = Path(path).read_text().splitlines()
    c1, c2 = lines[0], lines[1]
    f = lines[2].split()
    natoms = int(f[0])
    origin = np.array([float(x) for x in f[1:4]])
    shape, axes = [], []
    for i in range(3):
        g = lines[3 + i].split()
        shape.append(int(g[0]))
        axes.append([float(x) for x in g[1:4]])
    az, aq, axyz = [], [], []
    for i in range(abs(natoms)):
        g = lines[6 + i].split()
        az.append(int(g[0])); aq.append(float(g[1]))
        axyz.append([float(x) for x in g[2:5]])
    body = " ".join(lines[6 + abs(natoms):]).split()
    data = np.array(body, dtype=float).reshape(shape)
    return Cube(c1, c2, origin, tuple(shape), np.array(axes),
                np.array(az), np.array(aq), np.array(axyz), data)


def write_cube(cube: Cube, path: str | Path, comment2: str | None = None) -> Path:
    path = Path(path)
    out = [cube.comment1, comment2 if comment2 is not None else cube.comment2]
    out.append(f"{cube.natoms:5d} {cube.origin[0]:12.6f} {cube.origin[1]:12.6f} {cube.origin[2]:12.6f}")
    for i in range(3):
        out.append(f"{cube.shape[i]:5d} {cube.axes[i][0]:12.6f} {cube.axes[i][1]:12.6f} "
                   f"{cube.axes[i][2]:12.6f}")
    for z, q, p in zip(cube.atom_z, cube.atom_q, cube.atom_xyz):
        out.append(f"{int(z):5d} {q:12.6f} {p[0]:12.6f} {p[1]:12.6f} {p[2]:12.6f}")
    flat = cube.data.reshape(-1)
    chunk = []
    for i, v in enumerate(flat, 1):
        chunk.append(f"{v:13.5E}")
        if i % 6 == 0:
            out.append("".join(chunk)); chunk = []
    if chunk:
        out.append("".join(chunk))
    path.write_text("\n".join(out) + "\n")
    return path


def esp_cube(cfg: dict, name: str, solvent: str) -> Path | None:
    """Evaluate V(r) on the density-cube grid with orca_vpot and write esp.cube."""
    work = Path(cfg["execution"]["work_dir"]) / name / solvent
    dens = work / "eldens.cube"
    base = f"{name}_{solvent}_sp"
    gbw = work / f"{base}.gbw"
    target = work / "esp.cube"

    if target.exists() and cfg["execution"]["resume"]:
        ok(f"{name}/{solvent}: reusing {target.name}")
        return target
    if not dens.exists():
        warn(f"{name}/{solvent}: {dens.name} missing - cannot build the MEP")
        return None
    if not gbw.exists():
        warn(f"{name}/{solvent}: {gbw.name} missing - cannot build the MEP")
        return None

    cube = read_cube(dens)
    pts = cube.grid_points()
    step(f"{name}/{solvent}: evaluating the electrostatic potential at {len(pts):,} grid points")

    gridfile = work / "vpot_grid.xyz"
    with open(gridfile, "w") as fh:
        fh.write(f"{len(pts)}\n")
        np.savetxt(fh, pts, fmt="%.8f")

    potfile = work / "vpot.out"
    orca = Path(cfg["execution"]["orca_path"])
    nproc = int(cfg["cubes"]["esp_nprocs"]) or 1
    vpot = orca.parent / "orca_vpot"
    vpot_mpi = orca.parent / "orca_vpot_mpi"

    def _run_serial() -> subprocess.CompletedProcess:
        return subprocess.run(
            [str(vpot), gbw.name, f"{base}.scfp", gridfile.name, potfile.name],
            cwd=str(work), capture_output=True, text=True)

    proc = None
    if nproc > 1 and vpot_mpi.exists():
        # The MPI build re-reads the SHARK basis scratch (<base>.SHARKINP.tmp), and ORCA
        # deletes that file when a job terminates normally - so this path fails for any
        # calculation that actually finished. It is kept because it works when the scratch
        # survives, but a failure here is expected rather than exceptional.
        drv = work / "vpot_driver.inp"
        drv.write_text(f"{nproc}\n{gbw.name}\n{base}.scfp\n{gridfile.name}\n"
                       f"{potfile.name}\n")
        proc = subprocess.run([str(vpot), drv.name], cwd=str(work),
                              capture_output=True, text=True)
        if not potfile.exists():
            warn(f"{name}/{solvent}: the parallel orca_vpot could not run - it needs "
                 f"{base}.SHARKINP.tmp, which ORCA removes when a job ends normally. "
                 f"Evaluating the potential serially instead.")
            proc = _run_serial()
    else:
        proc = _run_serial()

    if not potfile.exists():
        detail = (proc.stdout or proc.stderr or "").strip()[-300:]
        warn(f"{name}/{solvent}: orca_vpot failed - {detail}")
        return None

    vals = np.loadtxt(potfile, skiprows=1, usecols=3)
    if vals.size != pts.shape[0]:
        warn(f"{name}/{solvent}: orca_vpot returned {vals.size} values for "
             f"{pts.shape[0]} points - MEP skipped")
        return None

    esp = Cube(cube.comment1, "Electrostatic potential (a.u.) on the electron-density grid",
               cube.origin, cube.shape, cube.axes, cube.atom_z, cube.atom_q, cube.atom_xyz,
               vals.reshape(cube.shape))
    write_cube(esp, target)
    gridfile.unlink(missing_ok=True)
    potfile.unlink(missing_ok=True)
    ok(f"{name}/{solvent}: {target.name}  "
       f"(V range {vals.min():+.4f} to {vals.max():+.4f} a.u.)")
    return target


def surface_esp_stats(cfg: dict, name: str, solvent: str) -> dict:
    """Min/max/mean of V(r) restricted to the molecular surface, for honest colour scaling."""
    work = Path(cfg["execution"]["work_dir"]) / name / solvent
    d, e = work / "eldens.cube", work / "esp.cube"
    if not (d.exists() and e.exists()):
        return {}
    dens, esp = read_cube(d), read_cube(e)
    iso = float(cfg["render"]["density_isovalue"])
    band = (dens.data > iso * 0.8) & (dens.data < iso * 1.25)
    if not band.any():
        return {}
    v = esp.data[band]
    # The raw extremes can sit on a single grid point; the 2nd/98th percentiles are the
    # honest colour limits and are what the figure scales to.
    lo, hi = np.percentile(v, [2.0, 98.0])
    return {"v_min": float(v.min()), "v_max": float(v.max()),
            "v_lo": float(lo), "v_hi": float(hi),
            "v_mean": float(v.mean()), "n_surface_points": int(band.sum())}


def generate_all(cfg: dict) -> dict:
    """Produce the ESP cube for every molecule in the figure solvent."""
    out = {}
    if not cfg["cubes"]["enabled"] or not cfg["cubes"]["esp"]:
        return out
    solvent = cfg["figures"]["reference_solvent"]
    for name in _M.names(cfg):
        p = esp_cube(cfg, name, solvent)
        if p:
            out[name] = p
    return out
