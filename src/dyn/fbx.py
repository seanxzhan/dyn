"""A minimal, dependency-free reader for binary FBX files.

FBX is a tree of typed nodes.  We parse it with only the standard library
(``struct`` + ``zlib``) into a lightweight ``FBXNode`` tree, then pull out the
pieces rig-space physics needs: the surface mesh, the skeleton (bones + bind
poses), and linear-blend-skinning weights.  No Autodesk SDK, no assimp.

Reference for the binary layout: the (reverse-engineered) Blender FBX format.
Files of version >= 7500 use 64-bit node headers; ours are 7700.
"""

from __future__ import annotations

import struct
import zlib
from dataclasses import dataclass

import numpy as np

from dyn.rig import finite_difference_jacobian

_ARRAY_DTYPE = {"f": "<f4", "d": "<f8", "l": "<i8", "i": "<i4", "b": "<i1"}
_SCALAR = {"Y": ("<h", 2), "C": ("<?", 1), "I": ("<i", 4),
           "F": ("<f", 4), "D": ("<d", 8), "L": ("<q", 8)}


class FBXNode:
    __slots__ = ("name", "props", "children")

    def __init__(self, name: str):
        self.name = name
        self.props: list = []
        self.children: list[FBXNode] = []

    def first(self, name: str):
        for c in self.children:
            if c.name == name:
                return c
        return None

    def find(self, name: str):
        return [c for c in self.children if c.name == name]

    def __repr__(self):
        return f"FBXNode({self.name!r}, {len(self.props)} props, {len(self.children)} kids)"


def _read_property(buf, off):
    t = chr(buf[off])
    off += 1
    if t in _SCALAR:
        fmt, size = _SCALAR[t]
        v = struct.unpack_from(fmt, buf, off)[0]
        return v, off + size
    if t in ("S", "R"):
        (ln,) = struct.unpack_from("<I", buf, off)
        off += 4
        raw = bytes(buf[off:off + ln])
        return (raw.decode("utf-8", "replace") if t == "S" else raw), off + ln
    if t in _ARRAY_DTYPE:
        arrlen, enc, comp = struct.unpack_from("<III", buf, off)
        off += 12
        payload = bytes(buf[off:off + comp])
        off += comp
        if enc == 1:
            payload = zlib.decompress(payload)
        return np.frombuffer(payload, dtype=np.dtype(_ARRAY_DTYPE[t]), count=arrlen), off
    raise ValueError(f"unknown FBX property type {t!r} at offset {off - 1}")


def _read_node(buf, off, is64):
    hdr = 24 if is64 else 12
    if is64:
        end_off, nprops, _plen = struct.unpack_from("<QQQ", buf, off)
    else:
        end_off, nprops, _plen = struct.unpack_from("<III", buf, off)
    name_len = buf[off + hdr]
    if end_off == 0 and nprops == 0 and name_len == 0:      # null/sentinel record
        return None, off + hdr + 1
    off += hdr + 1
    node = FBXNode(bytes(buf[off:off + name_len]).decode("utf-8", "replace"))
    off += name_len
    for _ in range(nprops):
        v, off = _read_property(buf, off)
        node.props.append(v)
    if off < end_off:                                       # nested child records
        while True:
            child, off = _read_node(buf, off, is64)
            if child is None:
                break
            node.children.append(child)
    return node, end_off


def parse(path: str) -> tuple[FBXNode, int]:
    """Parse an FBX file into a root FBXNode. Returns (root, version)."""
    with open(path, "rb") as f:
        buf = f.read()
    if buf[:20] != b"Kaydara FBX Binary  ":
        raise ValueError("not a binary FBX file")
    (version,) = struct.unpack_from("<I", buf, 23)
    is64 = version >= 7500
    off = 27
    root = FBXNode("__root__")
    null = 25 if is64 else 13
    while off < len(buf) - null:
        node, off = _read_node(buf, off, is64)
        if node is None:
            break
        root.children.append(node)
    return root, version


# ----------------------------------------------------------------------------
# Extracting a rigged character: surface mesh + skeleton + skin weights
# ----------------------------------------------------------------------------

def _mat4(node) -> np.ndarray:
    """FBX 4x4: stored column-major with row-vector convention (translation in
    the last row), so the usable column-vector matrix is the transpose."""
    return np.asarray(node.props[0], dtype=float).reshape(4, 4).T


def _triangulate(pvi: np.ndarray) -> np.ndarray:
    """PolygonVertexIndex -> triangle fan list.  Last index of each polygon is
    bit-flipped (negative): real index = ~idx = -idx-1."""
    faces, poly = [], []
    for raw in pvi:
        idx = int(raw)
        last = idx < 0
        poly.append(-idx - 1 if last else idx)
        if last:
            for k in range(1, len(poly) - 1):
                faces.append((poly[0], poly[k], poly[k + 1]))
            poly = []
    return np.array(faces, dtype=np.int64)


@dataclass
class Character:
    V: np.ndarray            # (Nv,3) rest world vertex positions (bind pose)
    faces: np.ndarray        # (Ntri,3) triangulated surface
    weights: np.ndarray      # (Nv, nb) linear-blend-skinning weights, rows sum to 1
    bone_names: list         # length nb
    bone_parent: np.ndarray  # (nb,) parent bone index, or -1 for a root
    B: np.ndarray            # (nb,4,4) bind-pose world transforms (TransformLink)
    order: list              # bone indices in topological order (parents first)


def load_character(path: str, mesh_index: int = 0) -> Character:
    """Parse a skinned FBX into a Character (mesh + skeleton + LBS weights)."""
    root, _ = parse(path)
    objs = root.first("Objects")
    conns = root.first("Connections")

    def oid(n):
        return n.props[0]

    def subtype(n):
        return n.props[2] if len(n.props) > 2 and isinstance(n.props[2], str) else ""

    def name_of(n):
        return n.props[1].split("\x00")[0] if len(n.props) > 1 and isinstance(n.props[1], str) else ""

    geoms = [c for c in objs.children if c.name == "Geometry" and subtype(c) == "Mesh"]
    bones = [c for c in objs.children if c.name == "Model" and subtype(c) == "LimbNode"]
    clusters = [c for c in objs.children if c.name == "Deformer" and subtype(c) == "Cluster"]

    geom = geoms[mesh_index]
    v_geom = np.asarray(geom.first("Vertices").props[0], dtype=float).reshape(-1, 3)
    faces = _triangulate(np.asarray(geom.first("PolygonVertexIndex").props[0]))

    # bone bookkeeping
    bone_idx = {oid(b): i for i, b in enumerate(bones)}
    nb = len(bones)
    B = np.tile(np.eye(4), (nb, 1, 1))
    bone_parent = np.full(nb, -1, dtype=int)
    for cc in conns.children:                       # bone -> bone parent edges
        if cc.props[0] == "OO" and cc.props[1] in bone_idx and cc.props[2] in bone_idx:
            bone_parent[bone_idx[cc.props[1]]] = bone_idx[cc.props[2]]

    # cluster -> bone:  connection is (child=bone Model, parent=cluster Deformer)
    clu_id = {oid(c) for c in clusters}
    clu2bone = {}
    for cc in conns.children:
        if cc.props[0] == "OO" and cc.props[1] in bone_idx and cc.props[2] in clu_id:
            clu2bone[cc.props[2]] = bone_idx[cc.props[1]]

    # mesh bind transform (same on every cluster) maps geometry -> world bind pose
    mesh_xf = np.eye(4)
    for c in clusters:
        if c.first("Transform") is not None:
            mesh_xf = _mat4(c.first("Transform"))
            break
    V = (np.hstack([v_geom, np.ones((len(v_geom), 1))]) @ mesh_xf.T)[:, :3]

    # skin weights + bind world transforms from each cluster
    weights = np.zeros((len(v_geom), nb))
    for c in clusters:
        b = clu2bone.get(oid(c))
        if b is None:
            continue
        if c.first("TransformLink") is not None:
            B[b] = _mat4(c.first("TransformLink"))
        idx_n, w_n = c.first("Indexes"), c.first("Weights")
        if idx_n is not None and w_n is not None and len(idx_n.props) and len(w_n.props):
            weights[np.asarray(idx_n.props[0], dtype=int), b] = np.asarray(w_n.props[0], dtype=float)

    # normalize; vertices with no influence ride the first root bone rigidly
    rowsum = weights.sum(1)
    roots = np.where(bone_parent < 0)[0]
    if len(roots):
        weights[rowsum == 0, roots[0]] = 1.0
        rowsum = weights.sum(1)
    weights /= np.where(rowsum > 0, rowsum, 1.0)[:, None]

    # topological order (parents before children)
    order, seen = [], set()

    def visit(b):
        if b in seen:
            return
        if bone_parent[b] >= 0:
            visit(bone_parent[b])
        seen.add(b)
        order.append(b)

    for b in range(nb):
        visit(b)

    return Character(V=V, faces=faces, weights=weights,
                     bone_names=[name_of(b) for b in bones],
                     bone_parent=bone_parent, B=B, order=order)


def _rodrigues(w: np.ndarray) -> np.ndarray:
    """Axis-angle vector -> 3x3 rotation (exponential map)."""
    theta = float(np.linalg.norm(w))
    if theta < 1e-12:
        return np.eye(3)
    k = w / theta
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) + np.sin(theta) * K + (1 - np.cos(theta)) * (K @ K)


class SkeletonRig:
    """Linear-blend-skinning rig driven by per-joint rotations.

    Implements the artist's skeleton as the black-box rig map s(p): each
    controllable bone rotates about one or more axes (propagated to descendants
    through the hierarchy), and every vertex follows

        x_v(p) = Σ_b w_{vb} · M_b(p) · B_b⁻¹ · p0_v,

    with B_b the bind-pose world transform and M_b(p) the posed world transform.
    At p=0, M_b = B_b so x = p0 (the bind pose) — verified in the lesson.  The
    Jacobian is taken by finite differences, exactly as the papers treat a
    black-box rig.

    ``axes`` gives, per control bone, the rotation axes it may use (each a unit
    3-vector); a bone with one axis is a 1-DOF hinge.  Default is full 3-DOF
    (the world axes) per bone.
    """

    def __init__(self, verts, weights, char: Character, control_bones, axes=None, amp: float = 1.0):
        self.W = np.asarray(weights, dtype=float)               # (M, nb)
        self.p0h = np.hstack([verts, np.ones((len(verts), 1))])  # (M,4)
        self.B = char.B
        self.Binv = np.linalg.inv(char.B)
        self.parent = char.bone_parent
        self.order = char.order
        self.amp = amp
        # bind-pose local transforms  L_b = B_parent⁻¹ B_b  (= B_b at a root)
        self.L = np.empty_like(char.B)
        for b in range(len(char.B)):
            par = self.parent[b]
            self.L[b] = char.B[b] if par < 0 else self.Binv[par] @ char.B[b]
        self.control = list(control_bones)
        if axes is None:
            axes = [np.eye(3)] * len(self.control)               # 3-DOF per bone
        self.axes = [np.atleast_2d(np.asarray(a, dtype=float)) for a in axes]
        self.axes = [a / np.linalg.norm(a, axis=1, keepdims=True) for a in self.axes]
        self.offset = np.cumsum([0] + [a.shape[0] for a in self.axes])
        self.dim = int(self.offset[-1])
        self.names = [f"{char.bone_names[b]}.{j}"
                      for bi, b in enumerate(self.control) for j in range(self.axes[bi].shape[0])]
        self.rest_params = np.zeros(self.dim)

    def _world_transforms(self, p):
        p = np.asarray(p, dtype=float)
        dR = np.tile(np.eye(4), (len(self.B), 1, 1))
        for bi, b in enumerate(self.control):
            ang = self.amp * p[self.offset[bi]:self.offset[bi + 1]]
            R = np.eye(3)
            for j, theta in enumerate(ang):
                R = R @ _rodrigues(theta * self.axes[bi][j])
            dR[b, :3, :3] = R
        M = np.empty_like(self.B)
        for b in self.order:
            par = self.parent[b]
            base = self.L[b] if par < 0 else M[par] @ self.L[b]
            M[b] = base @ dR[b]
        return M

    def s(self, p):
        S = self._world_transforms(p) @ self.Binv                  # (nb,4,4)
        skinned = np.einsum("bij,mj->bmi", S, self.p0h)            # (nb,M,4)
        return np.einsum("mb,bmi->mi", self.W, skinned[..., :3])   # (M,3)

    def jacobian(self, p):
        return finite_difference_jacobian(self, p)

