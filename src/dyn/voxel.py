"""Turn a surface mesh into a volumetric tet mesh — no external mesher.

Rig-space physics needs a *volumetric* body (tets with a strain energy), but an
FBX gives only a surface.  We build the volume the dependency-free way: test a
regular grid of points for being inside the surface (via the **generalized
winding number**, robust to non-watertight meshes), keep the inside voxels, and
split each into 6 tetrahedra with the same Kuhn subdivision as ``mesh.py``.  The
result is a blocky but honest volumetric proxy of the character.

We then transfer the surface's skinning weights onto the proxy vertices, so the
same skeleton drives the proxy — and the whole Lessons 4–8 pipeline runs on it.
"""

from __future__ import annotations

import numpy as np

from dyn.mesh import Beam, _VOXEL_TETS, _boundary_faces, _orient_positive, lumped_mass


def winding_number(query: np.ndarray, tri: np.ndarray, chunk: int = 512) -> np.ndarray:
    """Generalized winding number of points ``query`` (G,3) w.r.t. triangles
    ``tri`` (T,3,3).  ~+1 strictly inside a closed surface, ~0 outside.
    Uses the Van Oosterom–Strackee solid-angle formula, summed over triangles.
    """
    a, b, c = tri[:, 0], tri[:, 1], tri[:, 2]
    out = np.zeros(len(query))
    for s in range(0, len(query), chunk):
        q = query[s:s + chunk][:, None, :]                 # (g,1,3)
        va, vb, vc = a[None] - q, b[None] - q, c[None] - q  # (g,T,3)
        la = np.linalg.norm(va, axis=2)
        lb = np.linalg.norm(vb, axis=2)
        lc = np.linalg.norm(vc, axis=2)
        det = np.einsum("gti,gti->gt", va, np.cross(vb, vc))
        denom = (la * lb * lc
                 + np.einsum("gti,gti->gt", va, vb) * lc
                 + np.einsum("gti,gti->gt", vb, vc) * la
                 + np.einsum("gti,gti->gt", vc, va) * lb)
        out[s:s + chunk] = np.arctan2(det, denom).sum(1) / (2.0 * np.pi)
    return out


def voxelize(V: np.ndarray, faces: np.ndarray, res: int = 16,
             density: float = 1.0, iso: float = 0.5) -> Beam:
    """Voxelize the interior of a surface into a tet ``Beam`` proxy.

    ``res`` = number of voxels along the longest bounding-box axis.
    """
    tri = V[faces]
    lo, hi = V.min(0), V.max(0)
    h = (hi - lo).max() / res
    lo = lo - 0.5 * h                                       # pad by half a voxel
    n = np.maximum(np.ceil((hi + 0.5 * h - lo) / h).astype(int), 1)

    ijk = np.array([(i, j, k) for i in range(n[0]) for j in range(n[1]) for k in range(n[2])])
    centers = lo + (ijk + 0.5) * h
    inside = np.abs(winding_number(centers, tri)) > iso
    cells = ijk[inside]
    if len(cells) == 0:
        raise RuntimeError("no interior voxels found — increase `res`")

    # collect the corner nodes of kept voxels, de-duplicated by integer coords
    node_id: dict[tuple, int] = {}
    verts: list = []

    def nid(coord):
        key = tuple(int(t) for t in coord)
        if key not in node_id:
            node_id[key] = len(verts)
            verts.append(lo + np.array(key) * h)
        return node_id[key]

    tets = []
    for cx, cy, cz in cells:
        corner = {(di, dj, dk): nid((cx + di, cy + dj, cz + dk))
                  for di in (0, 1) for dj in (0, 1) for dk in (0, 1)}
        for tet in _VOXEL_TETS:
            tets.append([corner[c] for c in tet])

    verts = np.array(verts, dtype=float)
    tets = _orient_positive(verts, np.array(tets, dtype=np.int64))
    surf = _boundary_faces(tets)
    surf_set = set(int(i) for i in surf.reshape(-1))
    surface_idx = np.array(sorted(surf_set), dtype=int)
    interior_idx = np.array([i for i in range(len(verts)) if i not in surf_set], dtype=int)

    return Beam(verts=verts, tets=tets, surf_faces=surf,
                mass=lumped_mass(verts, tets, density),
                surface_idx=surface_idx, interior_idx=interior_idx,
                clamp_idx=np.array([], dtype=int), res=tuple(int(x) for x in n),
                lengths=tuple((hi - lo)))


def transfer_weights(src_verts: np.ndarray, src_weights: np.ndarray,
                     dst_verts: np.ndarray, k: int = 3) -> np.ndarray:
    """Skin the proxy: each proxy vertex inherits an inverse-distance blend of
    the ``k`` nearest surface vertices' weights.  Rows are renormalized to 1.
    """
    out = np.zeros((len(dst_verts), src_weights.shape[1]))
    for i, q in enumerate(dst_verts):
        d2 = ((src_verts - q) ** 2).sum(1)
        nn = np.argpartition(d2, min(k, len(d2) - 1))[:k]
        wts = 1.0 / (d2[nn] + 1e-12)
        out[i] = (wts[:, None] * src_weights[nn]).sum(0)
    rs = out.sum(1)
    return out / np.where(rs > 0, rs, 1.0)[:, None]
