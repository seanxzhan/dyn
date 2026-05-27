"""Rigs: maps from a few parameters ``p`` to all vertex positions ``s(p)``.

A *rig* is the low-dimensional handle an animator controls.  Rig-space physics
runs the simulation in the coordinates ``p`` (a handful of numbers) instead of
the ~thousands of vertices.  The object that makes this work is the **rig
Jacobian** ``J = ds/dp`` — it pushes per-vertex forces down into rig space as
``f_p = Jᵀ f_s`` (the mathematical heart of both papers).

Two rigs, chosen to expose the central distinction:

    LinearRig : s(p) = s0 + B p     ->  J = B is CONSTANT, curvature d²s/dp² = 0
    BendRig   : one rotating "bone"  ->  J(θ) VARIES, curvature d²s/dθ² ≠ 0

Constant ``J`` is exactly the regime the 2013 paper *engineers* ("linearise the
rig once per step", killing the O(p²) term).  Varying ``J`` is where the 2012
curvature term ``∂_p J`` lives.  See ``tutorials/02`` and ``tutorials/03``.

Convention: a vertex field ``(N,3)`` is flattened **row-major** to a vector of
length ``3N`` = [x0,y0,z0, x1,y1,z1, ...].  A Jacobian therefore has shape
``(3N, dim)``; its i-th column, reshaped to ``(N,3)``, is the displacement field
"which way every vertex moves when you nudge p[i]".
"""

from __future__ import annotations

import numpy as np


class Rig:
    """Interface: ``dim`` parameters, ``s(p) -> (N,3)``, ``jacobian(p) -> (3N,dim)``."""

    dim: int
    names: list[str]
    rest_params: np.ndarray

    def s(self, p):  # pragma: no cover - interface
        raise NotImplementedError

    def jacobian(self, p):  # pragma: no cover - interface
        raise NotImplementedError


class LinearRig(Rig):
    """A *linear* (blendshape-style) rig: ``s(p) = s0 + B p``.

    Each parameter adds a fixed displacement field (a column of ``B``).  Because
    the map is affine, the Jacobian ``J = B`` does not depend on ``p`` and the
    curvature ``d²s/dp²`` is exactly zero.  This is the cleanest possible rig and
    the regime the 2013 solver reduces every rig to, locally, each timestep.

    The built-in modes use a quadratic profile ``(x/Lx)²`` along the beam, which
    vanishes *with zero slope* at the clamped end x=0 — so the clamp stays fixed.
    """

    _PROFILES = {
        "bend_z": (2, lambda xn: xn ** 2),          # tip lifts in +z (clamp fixed)
        "sway_y": (1, lambda xn: xn ** 2),          # tip sways in +y (clamp fixed)
        "stretch_x": (0, lambda xn: xn),            # uniform stretch along the beam
        "lift_z": (2, lambda xn: np.ones_like(xn)),  # rigid translation in z (moves clamp too)
        "shift_y": (1, lambda xn: np.ones_like(xn)),  # rigid translation in y (moves clamp too)
    }

    def __init__(self, beam, modes=("bend_z", "sway_y"), amp: float = 1.0):
        self.beam = beam
        self.names = list(modes)
        self.dim = len(modes)
        self.rest_params = np.zeros(self.dim)
        self._s0 = beam.verts.reshape(-1).copy()

        xn = beam.verts[:, 0] / beam.lengths[0]   # normalized position along length, 0..1
        cols = []
        for m in modes:
            axis, profile = self._PROFILES[m]
            d = np.zeros_like(beam.verts)
            d[:, axis] = amp * profile(xn)
            cols.append(d.reshape(-1))
        self.B = np.column_stack(cols)            # (3N, dim), constant

    def s(self, p):
        p = np.asarray(p, dtype=float)
        return (self._s0 + self.B @ p).reshape(-1, 3)

    def jacobian(self, p=None):
        return self.B                              # constant — ignores p


# --- rotation about the y-axis and its derivatives w.r.t. the angle ----------
# R_y(a) rotates the (x,z) plane; the lessons bend the beam in this plane.
def _rot_y(a):       # (M,) angles -> (M,3,3) rotation matrices
    c, s, z, o = np.cos(a), np.sin(a), np.zeros_like(a), np.ones_like(a)
    return np.stack([np.stack([c, z, s], -1),
                     np.stack([z, o, z], -1),
                     np.stack([-s, z, c], -1)], axis=1)


def _rot_y_d(a):     # dR_y/da
    c, s, z = np.cos(a), np.sin(a), np.zeros_like(a)
    return np.stack([np.stack([-s, z, c], -1),
                     np.stack([z, z, z], -1),
                     np.stack([-c, z, -s], -1)], axis=1)


def _rot_y_dd(a):    # d²R_y/da²
    c, s, z = np.cos(a), np.sin(a), np.zeros_like(a)
    return np.stack([np.stack([-c, z, -s], -1),
                     np.stack([z, z, z], -1),
                     np.stack([s, z, -c], -1)], axis=1)


class BendRig(Rig):
    """A *nonlinear* one-bone bend (linear-blend skinning to a single joint).

    A single angle ``θ`` rotates each vertex about the y-axis through a pivot at
    the clamped end.  Each vertex turns by its *own* angle ``α(x) = θ·w(x)``,
    where ``w`` is a smoothstep weight ramping 0→1 from the clamp to the tip
    (``w(0)=0, w'(0)=0`` so the clamp stays put).  This is literally LBS with one
    bone and a smooth weight — the simplest cousin of a "jiggle bone".

        s_v(θ) = pivot + R_y(α_v) · (v_rest − pivot),     α_v = θ·w_v

    Because ``R_y`` is built from sin/cos, the map is nonlinear in ``θ``:

        ds_v/dθ   = w_v · R_y'(α_v) · offset_v
        d²s_v/dθ² = w_v² · R_y''(α_v) · offset_v        (nonzero -> rig curvature)

    The second derivative is the term that costs the 2012 method O(p²) rig
    evaluations and that the 2013 method discards by linearising the rig.
    """

    def __init__(self, beam, amp: float = 1.0):
        self.beam = beam
        self.dim = 1
        self.names = ["bend θ"]
        self.rest_params = np.zeros(1)
        self.amp = amp

        Lx = beam.lengths[0]
        t = np.clip(beam.verts[:, 0] / Lx, 0.0, 1.0)
        self.w = 3 * t ** 2 - 2 * t ** 3                     # smoothstep weight, (N,)
        self.pivot = np.array([0.0, beam.verts[:, 1].mean(), beam.verts[:, 2].mean()])
        self.offset = beam.verts - self.pivot                # (N,3)

    def _alpha(self, p):
        theta = float(np.asarray(p).reshape(-1)[0])
        return self.amp * theta * self.w                     # per-vertex angle (N,)

    def s(self, p):
        R = _rot_y(self._alpha(p))                           # (N,3,3)
        return self.pivot + np.einsum("nij,nj->ni", R, self.offset)

    def jacobian(self, p):
        dR = _rot_y_d(self._alpha(p))
        # chain rule: dα/dθ = amp·w, applied per vertex
        Jv = np.einsum("nij,nj->ni", dR, self.offset) * (self.amp * self.w)[:, None]
        return Jv.reshape(-1, 1)                              # (3N, 1)

    def second_derivative(self, p):
        """d²s/dθ² as an (N,3) field — the rig *curvature* you can visualize."""
        ddR = _rot_y_dd(self._alpha(p))
        return np.einsum("nij,nj->ni", ddR, self.offset) * ((self.amp * self.w) ** 2)[:, None]


def finite_difference_jacobian(rig: Rig, p, eps: float = 1e-6) -> np.ndarray:
    """Central-difference estimate of ``J = ds/dp`` — the black-box recipe.

    The papers compute ``J`` exactly this way when the rig has no analytic
    derivative (Maya, deformers, ...).  We use it in the lessons to *verify* our
    closed-form Jacobians: ``J_fd ≈ J_analytic`` to machine precision.
    """
    p = np.asarray(p, dtype=float)
    cols = []
    for i in range(rig.dim):
        pp, pm = p.copy(), p.copy()
        pp[i] += eps
        pm[i] -= eps
        cols.append(((rig.s(pp) - rig.s(pm)) / (2 * eps)).reshape(-1))
    return np.column_stack(cols)


def as_vertex_vectors(column: np.ndarray) -> np.ndarray:
    """Reshape a length-3N vector (e.g. one Jacobian column) into (N,3)."""
    return np.asarray(column).reshape(-1, 3)
