"""Solving the physics *in rig space*: reduced Newton, statics, and dynamics.

Everything here works in the reduced coordinate ``p`` (a handful of numbers).
The bridge is the rig Jacobian: a vertex-space gradient ``g_x`` becomes a
rig-space gradient ``g_p = Jᵀ g_x`` (Lesson 3).  Because ``p`` is tiny, we form
the reduced Hessian by simply finite-differencing the reduced gradient — and
that automatically includes the rig-curvature term ``∂_p J`` (the 2012 O(p²)
cost).  Freeze the rig's linear approximation (``LinearizedRig``) and the same
code drops the curvature — that is exactly the 2013 trick (Lesson 7).
"""

from __future__ import annotations

import numpy as np


def numerical_jacobian(f, x, eps: float = 1e-6) -> np.ndarray:
    """Central-difference Jacobian of f: R^d -> R^d, returned as (d,d)."""
    x = np.asarray(x, dtype=float)
    d = x.shape[0]
    J = np.zeros((d, d))
    for i in range(d):
        xp, xm = x.copy(), x.copy()
        xp[i] += eps
        xm[i] -= eps
        J[:, i] = (f(xp) - f(xm)) / (2.0 * eps)
    return J


def newton_minimize(energy, grad, p0, max_iter: int = 30, tol: float = 1e-8,
                    hess=None, line_search: bool = True, record=False):
    """Minimize ``energy(p)`` by Newton's method with a backtracking line search.

    ``grad`` returns the reduced gradient; ``hess`` (optional) the reduced
    Hessian — if omitted we finite-difference ``grad`` (the full 2012 Hessian,
    curvature included).  Returns ``(p, history)`` where history is a list of
    ``(p, energy, |grad|)`` per iteration (empty unless ``record``).
    """
    p = np.array(p0, dtype=float)
    history = []
    for _ in range(max_iter):
        g = grad(p)
        gn = float(np.linalg.norm(g))
        if record:
            history.append((p.copy(), float(energy(p)), gn))
        if gn < tol:
            break
        H = hess(p) if hess is not None else numerical_jacobian(grad, p)
        H = 0.5 * (H + H.T)
        # Levenberg-style fallback: nudge toward PD until the solve gives descent.
        dp = None
        for ridge in (0.0, 1e-6, 1e-3, 1e-1, 1e1):
            try:
                cand = np.linalg.solve(H + ridge * np.eye(H.shape[0]), -g)
            except np.linalg.LinAlgError:
                continue
            if cand @ g < 0:               # a descent direction
                dp = cand
                break
        if dp is None:
            dp = -g                        # plain gradient descent
        # backtracking (Armijo)
        a, E0, slope = 1.0, float(energy(p)), float(g @ dp)
        if line_search:
            for _ in range(40):
                if float(energy(p + a * dp)) <= E0 + 1e-4 * a * slope:
                    break
                a *= 0.5
        p = p + a * dp
    if record:
        history.append((p.copy(), float(energy(p)), float(np.linalg.norm(grad(p)))))
    return p, history


class LinearizedRig:
    """Frozen first-order rig: ``s(p) = s(p0) + J0 (p − p0)`` with ``J0`` constant.

    Wrapping a nonlinear rig in this each timestep is the 2013 "linearise the rig"
    trick: the curvature ``∂²s/∂p²`` becomes exactly zero, so the reduced Hessian
    is ``J0ᵀ K J0`` with no O(p²) rig term.
    """

    def __init__(self, rig, p0):
        self.dim = rig.dim
        self.names = rig.names
        self.p0 = np.array(p0, dtype=float)
        self.s0 = rig.s(self.p0)
        self.J0 = rig.jacobian(self.p0)
        self.rest_params = self.p0.copy()

    def s(self, p):
        p = np.asarray(p, dtype=float)
        return self.s0 + (self.J0 @ (p - self.p0)).reshape(-1, 3)

    def jacobian(self, p=None):
        return self.J0


def make_static_problem(rig, elastic, mass=None, gravity=None):
    """Build ``(energy, grad)`` for the rig-space *statics* problem.

    Minimizes elastic energy plus a gravity potential over the rig parameters:
    ``E(p) = W(s(p)) − Σ_i m_i g·s_i(p)``.  The minimizer is the physically
    relaxed pose the rig can reach under gravity.
    """
    g_vec = np.zeros(3) if gravity is None else np.asarray(gravity, dtype=float)
    fg = None if mass is None else mass[:, None] * g_vec        # gravity force (N,3)

    def energy(p):
        x = rig.s(p)
        E = elastic.energy(x)
        if fg is not None:
            E -= float(np.sum(fg * x))
        return E

    def grad(p):
        x = rig.s(p)
        gx = elastic.gradient(x)
        if fg is not None:
            gx = gx - fg
        return rig.jacobian(p).T @ gx.reshape(-1)

    return energy, grad


class ImplicitEuler:
    """Backward-Euler dynamics in rig space, solved as an energy minimization.

    Each step minimizes the *incremental potential*

        Φ(p) = 1/(2h²) (s(p) − x*)ᵀ M (s(p) − x*) + W(s(p)) − Σ m_i g·s_i(p),

    where ``x* = 2 x_n − x_{n-1}`` is the inertial (free-flight) prediction.  The
    stationary point is implicit Euler; minimizing over ``p`` keeps the motion in
    rig space.  Ad-hoc velocity damping models dissipation.
    """

    def __init__(self, rig, elastic, mass, h: float = 0.02,
                 gravity=(0.0, 0.0, -9.8), damping: float = 0.04, p0=None,
                 free_idx=None, relinearize: bool = False):
        self._base_rig = rig
        self.rig, self.elastic = rig, elastic
        self.M = mass[:, None]
        self.h = float(h)
        self.fg = self.M * np.asarray(gravity, dtype=float)
        self.damping = float(damping)
        self.relinearize = relinearize              # 2013: freeze the rig each step
        self.p = rig.rest_params.copy() if p0 is None else np.array(p0, dtype=float)
        # Which parameters the physics solves for; the rest are *driven* (primary
        # animation) and prescribed each step via step(driven_values=...).
        self.free_idx = np.arange(rig.dim) if free_idx is None else np.asarray(free_idx, dtype=int)
        self.driven_idx = np.array([i for i in range(rig.dim)
                                    if i not in set(self.free_idx.tolist())], dtype=int)
        self.x = rig.s(self.p)
        self.x_prev = self.x.copy()                 # start at rest (v = 0)

    def _prediction(self):
        return 2.0 * self.x - self.x_prev           # x* = x_n + h v_n

    def incremental_energy(self, p):
        x = self.rig.s(p)
        d = x - self._prediction()
        inertia = 0.5 / self.h ** 2 * float(np.sum(self.M * d * d))
        return inertia + self.elastic.energy(x) - float(np.sum(self.fg * x))

    def incremental_gradient(self, p):
        x = self.rig.s(p)
        gx = self.M * (x - self._prediction()) / self.h ** 2 \
            + self.elastic.gradient(x) - self.fg
        return self.rig.jacobian(p).T @ gx.reshape(-1)

    def step(self, driven_values=None, newton_iter: int = 8):
        if driven_values is not None and self.driven_idx.size:
            self.p[self.driven_idx] = np.asarray(driven_values, dtype=float)
        if self.relinearize:                        # 2013: linearize rig about p_n
            self.rig = LinearizedRig(self._base_rig, self.p)

        def embed(free):                            # free sub-vector -> full p
            q = self.p.copy()
            q[self.free_idx] = free
            return q

        energy = lambda free: self.incremental_energy(embed(free))         # noqa: E731
        grad = lambda free: self.incremental_gradient(embed(free))[self.free_idx]  # noqa: E731
        free_new, _ = newton_minimize(energy, grad, self.p[self.free_idx], max_iter=newton_iter)

        p_new = embed(free_new)
        x_new = self.rig.s(p_new)
        v = (x_new - self.x) / self.h
        v *= (1.0 - self.damping)                   # ad-hoc damping
        self.x_prev = x_new - self.h * v            # bake damped velocity into history
        self.x = x_new
        self.p = p_new
        return p_new

    def kinetic_energy(self) -> float:
        v = (self.x - self.x_prev) / self.h
        return 0.5 * float(np.sum(self.M * v * v))

    def pluck(self, vel_field):
        """Inject a velocity (N,3): rewind history so x* carries this velocity."""
        self.x_prev = self.x - self.h * np.asarray(vel_field, dtype=float)
