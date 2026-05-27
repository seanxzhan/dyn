"""Procedural beam geometry for the rig-space physics lessons.

Everything here is pure NumPy: we build a rectangular *beam* (an elongated
cuboid), split it into tetrahedra analytically, tag which vertices are on the
surface / interior / clamped, and compute a lumped (diagonal) mass.  No external
mesher is needed, and because the geometry is analytic we can later check every
formula against closed form.

This is exactly the object the papers start from (their §"spatial
discretization"): a tetrahedral mesh with a lumped mass matrix, whose vertices
split into surface nodes ``s`` (driven by the rig) and interior nodes ``n``.
See ``tutorials/01_polyscope_and_mesh.md`` for the narrative.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# The 8 corners of one voxel are addressed by (di, dj, dk) in {0,1}^3.  We split
# every voxel into 6 tetrahedra that all share the main diagonal (0,0,0)->(1,1,1).
# This "Kuhn / Freudenthal" subdivision is *consistent* between neighbouring
# voxels, so the triangular faces on a shared wall match up and the global tet
# mesh is conforming (no cracks).
_VOXEL_TETS = (
    ((0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1)),
    ((0, 0, 0), (1, 0, 0), (1, 0, 1), (1, 1, 1)),
    ((0, 0, 0), (0, 1, 0), (1, 1, 0), (1, 1, 1)),
    ((0, 0, 0), (0, 1, 0), (0, 1, 1), (1, 1, 1)),
    ((0, 0, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1)),
    ((0, 0, 0), (0, 0, 1), (0, 1, 1), (1, 1, 1)),
)


@dataclass
class Beam:
    """A tetrahedralized cuboid beam plus the bookkeeping the lessons need."""

    verts: np.ndarray         # (N, 3) rest positions
    tets: np.ndarray          # (T, 4) int, positively oriented (signed vol > 0)
    surf_faces: np.ndarray    # (F, 3) int, boundary triangles (surface rendering)
    mass: np.ndarray          # (N,) lumped nodal masses
    surface_idx: np.ndarray   # vertices on the boundary  (the rig drives these)
    interior_idx: np.ndarray  # vertices strictly inside   (free FEM nodes later)
    clamp_idx: np.ndarray     # clamped Dirichlet vertices (default: the x=0 face)
    res: tuple                # (nx, ny, nz) voxel resolution
    lengths: tuple            # (Lx, Ly, Lz)

    @property
    def n(self) -> int:
        return self.verts.shape[0]


def tet_volumes(verts: np.ndarray, tets: np.ndarray) -> np.ndarray:
    """Signed volume of each tet, V = det[e1 e2 e3] / 6 with ei = vi - v0.

    Positive when the four nodes form a right-handed (correctly oriented) tet.
    """
    v0 = verts[tets[:, 0]]
    e1 = verts[tets[:, 1]] - v0
    e2 = verts[tets[:, 2]] - v0
    e3 = verts[tets[:, 3]] - v0
    return np.einsum("ij,ij->i", np.cross(e1, e2), e3) / 6.0


def _orient_positive(verts: np.ndarray, tets: np.ndarray) -> np.ndarray:
    """Swap two nodes of any negatively-oriented tet so every signed vol > 0."""
    tets = tets.copy()
    flip = tet_volumes(verts, tets) < 0
    tets[flip] = tets[flip][:, [0, 1, 3, 2]]  # swapping nodes 2,3 negates det
    return tets


def _boundary_faces(tets: np.ndarray) -> np.ndarray:
    """Triangles that belong to exactly one tet are on the boundary surface."""
    count: dict[tuple, int] = {}
    for t in tets:
        for tri in ((t[0], t[1], t[2]), (t[0], t[1], t[3]),
                    (t[0], t[2], t[3]), (t[1], t[2], t[3])):
            key = tuple(sorted(int(x) for x in tri))
            count[key] = count.get(key, 0) + 1
    return np.array([k for k, v in count.items() if v == 1], dtype=np.int64)


def lumped_mass(verts: np.ndarray, tets: np.ndarray, density: float = 1.0) -> np.ndarray:
    """Lumped (diagonal) mass: split each tet's mass equally to its 4 nodes.

    M_ii = density * (1/4) * sum of volumes of tets touching node i.  This is the
    diagonal "lumped mass matrix" the papers use (cheap, and good enough here).
    """
    vol = np.abs(tet_volumes(verts, tets))
    m = np.zeros(verts.shape[0])
    np.add.at(m, tets.reshape(-1), np.repeat(density * vol / 4.0, 4))
    return m


def make_beam(res=(10, 3, 3), lengths=(4.0, 1.0, 1.0), density: float = 1.0) -> Beam:
    """Build a cantilever beam: a cuboid of ``lengths`` split into ``res`` voxels.

    The beam runs along +x; the x=0 face is the clamped (fixed) end.
    """
    nx, ny, nz = res
    Lx, Ly, Lz = lengths
    xs, ys, zs = np.linspace(0, Lx, nx + 1), np.linspace(0, Ly, ny + 1), np.linspace(0, Lz, nz + 1)

    # Vertices on a regular grid; vid maps grid coords (i,j,k) -> a flat index.
    def vid(i, j, k):
        return (i * (ny + 1) + j) * (nz + 1) + k

    ijk = np.array([(i, j, k) for i in range(nx + 1)
                    for j in range(ny + 1) for k in range(nz + 1)])
    verts = np.column_stack([xs[ijk[:, 0]], ys[ijk[:, 1]], zs[ijk[:, 2]]]).astype(float)

    # Tetrahedra: 6 per voxel, using the consistent Kuhn subdivision above.
    tets = np.array(
        [[vid(i + di, j + dj, k + dk) for (di, dj, dk) in tet]
         for i in range(nx) for j in range(ny) for k in range(nz)
         for tet in _VOXEL_TETS],
        dtype=np.int64,
    )
    tets = _orient_positive(verts, tets)

    # Node tags from grid coordinates: on the surface if any index is at a face.
    on_surf = ((ijk[:, 0] == 0) | (ijk[:, 0] == nx)
               | (ijk[:, 1] == 0) | (ijk[:, 1] == ny)
               | (ijk[:, 2] == 0) | (ijk[:, 2] == nz))

    return Beam(
        verts=verts,
        tets=tets,
        surf_faces=_boundary_faces(tets),
        mass=lumped_mass(verts, tets, density),
        surface_idx=np.where(on_surf)[0],
        interior_idx=np.where(~on_surf)[0],
        clamp_idx=np.where(ijk[:, 0] == 0)[0],
        res=res,
        lengths=lengths,
    )
