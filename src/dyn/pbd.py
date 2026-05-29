"""Position Based Dynamics — particle simulation by constraint projection.

A from-scratch, vectorized NumPy implementation of Müller et al. 2006
(``docs/Position-Based-Dynamics.md``).  The whole framework is one loop:

    for each step:
        v ← v + h * w * f_ext              (explicit external forces)
        p ← x + h * v                       (PREDICT)
        for it in range(n_iters):
            project every constraint        (Δp = −s·w·∇C, Gauss-Seidel)
        v ← (p − x) / h                     (recover velocity)
        x ← p                               (commit)

This module ships the pieces every lesson needs: a cloth grid factory, the
predict step, vectorized projection of distance / dihedral-bend / sphere
constraints, and a ``DynamicSphere`` for two-way ball-cloth coupling (used in
the capstone).  Pinning is just ``w_i = 0``: the projection formula
``Δp_i = −s w_i ∇C`` makes pinned vertices invariant for free.

Everything is checkable in closed form: the dihedral-bend gradient is verified
against finite differences in the unit test at the bottom of this file.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


# ---------------------------------------------------------------------------
# Cloth: a triangle-meshed rectangle with stretch edges and bend pairs.
# ---------------------------------------------------------------------------


@dataclass
class Cloth:
    """A square-grid cloth with the bookkeeping every lesson needs."""

    verts: np.ndarray         # (N, 3) rest positions
    tris: np.ndarray          # (T, 3) int, triangle indices (CCW from +n side)
    edges: np.ndarray         # (E, 2) int, stretch-constraint edges
    rest_lengths: np.ndarray  # (E,) rest length of each stretch edge
    bend_pairs: np.ndarray    # (B, 4) int columns are (p1, p2, p3, p4)
    rest_dihedrals: np.ndarray  # (B,) rest dihedral angle (rad), planar = 0
    pin_idx: np.ndarray       # vertex indices to pin (w = 0)
    res: tuple                # (nx, ny) grid resolution
    size: tuple               # (Lx, Ly)

    @property
    def n(self) -> int:
        return self.verts.shape[0]


def make_cloth(res=(20, 20), size=(1.0, 1.0), pin="top",
               normal=(0.0, 1.0, 0.0)) -> Cloth:
    """Build a cloth as a flat triangle-meshed rectangle, hung in 3D.

    The cloth is laid in the xz-plane by default (normal = +y), with x
    spanning [-Lx/2, Lx/2] and z spanning [0, Lz] going up.  Triangulation
    splits each grid quad along its (i,j)→(i+1,j+1) diagonal so neighbouring
    quads share consistent edge orientations.

    Stretch edges = every triangle edge (deduplicated).  Bend pairs = every
    *interior* edge shared by two triangles, with p3, p4 the off-edge corners.

    pin: "top" pins the top row (largest z), "corners" pins the two top
    corners, "none" leaves the cloth completely free.
    """
    nx, nz = res
    Lx, Lz = size
    xs = np.linspace(-Lx / 2.0, Lx / 2.0, nx)
    zs = np.linspace(0.0, Lz, nz)
    X, Z = np.meshgrid(xs, zs, indexing="ij")           # (nx, nz)
    verts = np.zeros((nx * nz, 3))
    # place the cloth on the plane normal to `normal` through the origin
    n = np.asarray(normal, dtype=float)
    n /= np.linalg.norm(n)
    # pick two basis vectors perpendicular to n: u (in xz/xy plane) and v (up-ish)
    up = np.array([0.0, 0.0, 1.0])
    if abs(n @ up) > 0.99:
        up = np.array([0.0, 1.0, 0.0])
    v = up - (up @ n) * n; v /= np.linalg.norm(v)        # vertical-ish in cloth
    u = np.cross(v, n)                                   # horizontal in cloth
    verts = (X.ravel()[:, None] * u + Z.ravel()[:, None] * v)

    def vid(i, j):  # (i along u, j along v) → flat index
        return i * nz + j

    # triangles: each quad → 2 tris, split along the (i,j)→(i+1,j+1) diagonal
    tris = []
    for i in range(nx - 1):
        for j in range(nz - 1):
            a, b, c, d = vid(i, j), vid(i + 1, j), vid(i + 1, j + 1), vid(i, j + 1)
            tris.append((a, b, c))
            tris.append((a, c, d))
    tris = np.array(tris, dtype=np.int32)

    # unique stretch edges
    e = np.vstack([tris[:, [0, 1]], tris[:, [1, 2]], tris[:, [2, 0]]])
    e = np.sort(e, axis=1)
    edges = np.unique(e, axis=0)
    rest_lengths = np.linalg.norm(verts[edges[:, 0]] - verts[edges[:, 1]], axis=1)

    # bend pairs: each interior edge shared by two triangles, with off-edge corners
    edge_to_tri: dict[tuple[int, int], list[tuple[int, int]]] = {}
    for ti, t in enumerate(tris):
        for k in range(3):
            a, b = int(t[k]), int(t[(k + 1) % 3])
            opp = int(t[(k + 2) % 3])               # vertex across this edge
            key = (min(a, b), max(a, b))
            edge_to_tri.setdefault(key, []).append((ti, opp))
    bend = []
    for (a, b), pairs in edge_to_tri.items():
        if len(pairs) == 2:
            (_, p3), (_, p4) = pairs
            bend.append((a, b, p3, p4))
    bend_pairs = np.array(bend, dtype=np.int32) if bend else np.zeros((0, 4), int)
    rest_dihedrals = _compute_dihedrals(verts, bend_pairs) if len(bend_pairs) else np.zeros(0)

    # pinning
    if pin == "top":
        z = verts[:, 2] if abs(n[2]) < 0.99 else verts[:, 1]
        pin_idx = np.where(np.isclose(z, z.max()))[0]
    elif pin == "corners":
        z = verts[:, 2] if abs(n[2]) < 0.99 else verts[:, 1]
        x = verts[:, 0]
        top = np.where(np.isclose(z, z.max()))[0]
        pin_idx = top[[np.argmin(verts[top, 0]), np.argmax(verts[top, 0])]]
    elif pin == "none":
        pin_idx = np.zeros(0, dtype=np.int32)
    else:
        raise ValueError(f"unknown pin mode: {pin!r}")

    return Cloth(verts=verts, tris=tris, edges=edges, rest_lengths=rest_lengths,
                 bend_pairs=bend_pairs, rest_dihedrals=rest_dihedrals,
                 pin_idx=pin_idx, res=res, size=size)


def particle_inv_masses(cloth: Cloth, total_mass: float = 1.0) -> np.ndarray:
    """Per-vertex inverse mass.  Pinned vertices get w = 0 (immovable)."""
    m = np.full(cloth.n, total_mass / cloth.n)
    w = 1.0 / m
    w[cloth.pin_idx] = 0.0
    return w


# ---------------------------------------------------------------------------
# The PBD loop
# ---------------------------------------------------------------------------


def predict(x: np.ndarray, v: np.ndarray, w: np.ndarray,
            h: float, gravity=(0.0, 0.0, -9.8)) -> np.ndarray:
    """Step (5) + (7): apply external force then explicit Euler to get p.

    Pinned vertices (w == 0) don't see external forces and don't move.
    """
    g = np.asarray(gravity, dtype=float)
    movable = (w > 0)[:, None]
    v_new = np.where(movable, v + h * g, v)              # gravity is the only ext force
    return x + h * v_new


def project_distance(p: np.ndarray, edges: np.ndarray, rest_lengths: np.ndarray,
                     w: np.ndarray, k_prime: float = 1.0) -> None:
    """Project distance constraints  C = |p_a − p_b| − l₀  in place.

    True **Gauss-Seidel** sweep — each edge sees the updates from the previous
    edges in the same iteration.  This is the order Müller §3.2 specifies; a
    naive Jacobi (vectorized) sweep is unstable on chains at k'=1 because
    interior particles receive two simultaneous corrections that overshoot.
    The Python loop is fine: cloth has ~10⁴ edges and PBD is constraint-
    bound, not solver-bound.
    """
    for ei in range(edges.shape[0]):
        a, b = int(edges[ei, 0]), int(edges[ei, 1])
        wa, wb = w[a], w[b]
        denom = wa + wb
        if denom <= 0.0:
            continue                                     # both pinned
        delta = p[a] - p[b]
        d = float(np.linalg.norm(delta))
        if d < 1e-12:
            continue
        n = delta / d
        s = (d - rest_lengths[ei]) / denom * k_prime
        p[a] -= wa * s * n
        p[b] += wb * s * n


def project_bend(p: np.ndarray, bend_pairs: np.ndarray, rest_dihedrals: np.ndarray,
                 w: np.ndarray, k_prime: float = 1.0) -> None:
    """Project dihedral-bend constraints  C = arccos(n₁·n₂) − φ₀  in place.

    For each pair (p1, p2, p3, p4) sharing edge (p1, p2):
        n1 = ((p3-p1) × (p2-p1))             ̂  (normal of triangle p1,p3,p2)
        n2 = ((p2-p1) × (p4-p1))             ̂  (normal of triangle p1,p2,p4)
        C  = arccos(n1·n2) − φ₀

    Gradients via the geometric form (Bridson et al.):
        ∇p3 C = +(1/h_a) n1
        ∇p4 C = +(1/h_b) n2
        ∇p1 C = −(α3/h_a) n1 − (α4/h_b) n2  (with α_k from edge projection)
        ∇p2 C = −((1−α3)/h_a) n1 − ((1−α4)/h_b) n2
    where h_a, h_b are perpendicular distances from p3, p4 to the shared edge.
    """
    if len(bend_pairs) == 0:
        return
    for k in range(bend_pairs.shape[0]):
        i1, i2, i3, i4 = (int(bend_pairs[k, 0]), int(bend_pairs[k, 1]),
                          int(bend_pairs[k, 2]), int(bend_pairs[k, 3]))
        p1, p2, p3, p4 = p[i1], p[i2], p[i3], p[i4]
        e = p2 - p1
        e_len = float(np.linalg.norm(e)) + 1e-12
        e_hat = e / e_len
        r3 = p3 - p1
        r4 = p4 - p1
        a3 = float(r3 @ e_hat)
        a4 = float(r4 @ e_hat)
        u3 = r3 - a3 * e_hat
        u4 = r4 - a4 * e_hat
        h_a = float(np.linalg.norm(u3)) + 1e-12
        h_b = float(np.linalg.norm(u4)) + 1e-12
        n1 = np.cross(p3 - p1, p2 - p1)
        n2 = np.cross(p2 - p1, p4 - p1)
        n1 /= (np.linalg.norm(n1) + 1e-12)
        n2 /= (np.linalg.norm(n2) + 1e-12)
        # robust angle via atan2(|cross|, dot) — avoids the arccos catastrophic
        # cancellation near n1·n2 = 1 that injected ~1e-4 spurious bend on a
        # perfectly flat cloth
        cross_mag = float(np.linalg.norm(np.cross(n1, n2)))
        dot = float(n1 @ n2)
        C = float(np.arctan2(cross_mag, dot) - rest_dihedrals[k])
        if abs(C) < 1e-8:
            continue
        # geometric gradient (verified by FD in verify_bend_gradient):
        # rotating tri-1 by Δθ moves p3 by h_a·Δθ along n1 ⇒ ∇p3 C = n1/h_a.
        # The (1−f), f split for p1, p2 follows from where each perpendicular
        # foot lands along the shared edge.
        g3 = n1 / h_a
        g4 = n2 / h_b
        f3, f4 = a3 / e_len, a4 / e_len
        g1 = -((1.0 - f3) * g3 + (1.0 - f4) * g4)
        g2 = -(f3 * g3 + f4 * g4)
        w1, w2, w3, w4 = w[i1], w[i2], w[i3], w[i4]
        denom = (w1 * (g1 @ g1) + w2 * (g2 @ g2)
                 + w3 * (g3 @ g3) + w4 * (g4 @ g4))
        if denom <= 1e-12:
            continue
        s = C / denom * k_prime
        p[i1] -= w1 * s * g1
        p[i2] -= w2 * s * g2
        p[i3] -= w3 * s * g3
        p[i4] -= w4 * s * g4


def project_sphere(p: np.ndarray, w: np.ndarray, center: np.ndarray,
                   radius: float, thickness: float = 0.0) -> np.ndarray:
    """Inequality constraint: each particle stays outside a sphere.

    Returns a boolean mask of which particles were projected (so the caller
    can apply friction/restitution to their velocities afterward).
    """
    delta = p - center
    d = np.linalg.norm(delta, axis=1)
    R = radius + thickness
    inside = d < R
    if not inside.any():
        return inside
    # only move movable particles; fixed ones stay (and the sphere can't push them)
    move = inside & (w > 0)
    n = delta[move] / np.maximum(d[move], 1e-12)[:, None]
    p[move] = center + R * n
    return inside


def damp_velocities(x: np.ndarray, v: np.ndarray, m: np.ndarray,
                    k_damping: float) -> None:
    """Müller §3.5 damping — kill internal velocity, preserve rigid modes.

    Computes the system's centre-of-mass velocity v_cm and angular velocity
    ω = I⁻¹ L, then damps each vᵢ towards the rigid motion v_cm + ω×rᵢ.  At
    k_damping = 0 nothing happens; at k_damping = 1 only the rigid mode
    survives.  (Pinned particles, m = ∞, are skipped.)
    """
    if k_damping <= 0:
        return
    finite = np.isfinite(m)                              # exclude pinned (inf mass)
    if finite.sum() < 2:
        return
    xs, vs, ms = x[finite], v[finite], m[finite]
    M = ms.sum()
    x_cm = (xs * ms[:, None]).sum(axis=0) / M
    v_cm = (vs * ms[:, None]).sum(axis=0) / M
    r = xs - x_cm
    L = np.cross(r, ms[:, None] * vs).sum(axis=0)
    # inertia tensor I = Σ m (‖r‖²I − r rᵀ)
    rr = (r[:, :, None] * r[:, None, :])                 # (N, 3, 3)
    I = (ms[:, None, None] * (np.eye(3)[None] * (r * r).sum(axis=1)[:, None, None] - rr)).sum(axis=0)
    try:
        omega = np.linalg.solve(I, L)
    except np.linalg.LinAlgError:
        omega = np.zeros(3)
    dv = v_cm + np.cross(omega, r) - vs
    v[finite] = vs + k_damping * dv


# ---------------------------------------------------------------------------
# Dynamic sphere — for two-way ball-cloth coupling (capstone).
# ---------------------------------------------------------------------------


@dataclass
class DynamicSphere:
    """A rigid sphere with mass.  Cloth-vertex projections push it back."""

    center: np.ndarray
    radius: float
    mass: float = 1.0
    velocity: np.ndarray = field(default_factory=lambda: np.zeros(3))
    impulse: np.ndarray = field(default_factory=lambda: np.zeros(3))

    def predict(self, h: float, gravity=(0.0, 0.0, -9.8)) -> np.ndarray:
        """Free-flight prediction; collisions will perturb it."""
        g = np.asarray(gravity, dtype=float)
        self.velocity = self.velocity + h * g
        return self.center + h * self.velocity

    def commit(self, p_predicted: np.ndarray, h: float) -> None:
        """Apply accumulated impulses (Σ mᵢ Δpᵢ / Δt) and integrate position."""
        # 1/M · Σ impulse converts the cloth-side displacement into a velocity kick
        self.velocity = self.velocity + self.impulse / self.mass
        self.center = p_predicted + h * (self.impulse / self.mass)
        self.impulse = np.zeros(3)


def project_dynamic_sphere(p: np.ndarray, x: np.ndarray, w: np.ndarray,
                           cloth_mass: np.ndarray,
                           sphere: DynamicSphere, h: float,
                           thickness: float = 0.0) -> np.ndarray:
    """Project cloth particles outside the sphere AND record impulses on it.

    Per Müller §4.2: each Δpᵢ from a sphere collision is also applied to the
    sphere as an impulse mᵢΔpᵢ/h with opposite sign.
    """
    delta = p - sphere.center
    d = np.linalg.norm(delta, axis=1)
    R = sphere.radius + thickness
    inside = d < R
    move = inside & (w > 0)
    if not move.any():
        return inside
    n = delta[move] / np.maximum(d[move], 1e-12)[:, None]
    p_new = sphere.center + R * n
    dp = p_new - p[move]                                  # the projection vector
    # impulse on the sphere is the OPPOSITE of the impulse the sphere imparts on cloth
    sphere.impulse = sphere.impulse - (cloth_mass[move, None] * dp).sum(axis=0) / h
    p[move] = p_new
    return inside


# ---------------------------------------------------------------------------
# Verification: dihedral-bend gradient vs. finite differences.
# ---------------------------------------------------------------------------


def _compute_dihedrals(verts: np.ndarray, bend_pairs: np.ndarray) -> np.ndarray:
    """Rest dihedral angles for a set of bend pairs (in [0, π], 0 if planar)."""
    if len(bend_pairs) == 0:
        return np.zeros(0)
    i1, i2, i3, i4 = bend_pairs[:, 0], bend_pairs[:, 1], bend_pairs[:, 2], bend_pairs[:, 3]
    p1, p2, p3, p4 = verts[i1], verts[i2], verts[i3], verts[i4]
    n1 = np.cross(p3 - p1, p2 - p1)
    n2 = np.cross(p2 - p1, p4 - p1)
    n1 /= (np.linalg.norm(n1, axis=1, keepdims=True) + 1e-12)
    n2 /= (np.linalg.norm(n2, axis=1, keepdims=True) + 1e-12)
    cross_mag = np.linalg.norm(np.cross(n1, n2), axis=1)
    dot = (n1 * n2).sum(axis=1)
    return np.arctan2(cross_mag, dot)


def _bend_C(p1, p2, p3, p4, phi0):
    """Plain ``C(p1,p2,p3,p4) = atan2(|n1×n2|, n1·n2) − φ₀`` for FD verification."""
    n1 = np.cross(p3 - p1, p2 - p1); n1 /= np.linalg.norm(n1) + 1e-12
    n2 = np.cross(p2 - p1, p4 - p1); n2 /= np.linalg.norm(n2) + 1e-12
    return float(np.arctan2(np.linalg.norm(np.cross(n1, n2)), n1 @ n2)) - phi0


def verify_bend_gradient(seed: int = 0, eps: float = 1e-6) -> dict:
    """Sanity-check the analytic bend gradient via central differences.

    Returns a dict of max errors and the numerical residual.  Used by the
    Lesson-04 script to assert the implementation is correct.
    """
    rng = np.random.default_rng(seed)
    # take a non-flat configuration so 1/sin(C) is well-conditioned
    pts = rng.normal(size=(4, 3))
    p_flat = pts.reshape(-1)

    # one isolated bend pair with rest dihedral 0 (so C is the angle itself)
    pair = np.array([[0, 1, 2, 3]])
    rest = np.zeros(1)
    w = np.ones(4)

    # analytic Δp from project_bend (one shot at full strength k = 1)
    p_test = pts.copy()
    project_bend(p_test, pair, rest, w, k_prime=1.0)
    dp_analytic = p_test - pts

    # numerical: ∇C and apply the same Δp = −s·w·∇C formula
    grad_num = np.zeros((4, 3))
    for k in range(4):
        for c in range(3):
            ph, pl = pts.copy(), pts.copy()
            ph[k, c] += eps; pl[k, c] -= eps
            Cph = _bend_C(ph[0], ph[1], ph[2], ph[3], 0.0)
            Cpl = _bend_C(pl[0], pl[1], pl[2], pl[3], 0.0)
            grad_num[k, c] = (Cph - Cpl) / (2 * eps)
    C0 = _bend_C(pts[0], pts[1], pts[2], pts[3], 0.0)
    s_num = C0 / (grad_num * grad_num).sum()              # all wᵢ = 1
    dp_num = -s_num * grad_num

    err_dp = float(np.max(np.abs(dp_analytic - dp_num)))
    # a single Newton step *will not* drive C to 0 for a nonlinear function;
    # we just check the analytic step matches the FD step at the same point.
    return {"max_dp_err": err_dp, "C0": C0,
           "dp_analytic": dp_analytic, "dp_numeric": dp_num}


if __name__ == "__main__":
    info = verify_bend_gradient()
    print("bend gradient FD check:")
    print(f"  C0 = {info['C0']:.4f}")
    print(f"  max |dp_analytic - dp_numeric| = {info['max_dp_err']:.2e}")
    assert info["max_dp_err"] < 1e-6, "bend gradient mismatch!"
    print("OK")

    cloth = make_cloth(res=(8, 6), size=(1.0, 1.0))
    print(f"cloth: {cloth.n} verts, {len(cloth.tris)} tris, "
          f"{len(cloth.edges)} edges, {len(cloth.bend_pairs)} bend pairs, "
          f"{len(cloth.pin_idx)} pinned")
