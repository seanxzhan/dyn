"""Lesson 05 — Statics in rig space: minimize energy over the rig parameters.

    Run :  python src/dyn/tutorials/05_statics_in_rig_space.py
    Read:  src/dyn/tutorials/05_statics_in_rig_space.md

The beam sags under gravity, but only within rig space.  We Newton-minimize
E(p) = W(s(p)) − gravity over the 2 rig parameters.  The energy landscape is
drawn as a bowl beside the beam, with the Newton iterates walking down to the
minimum.  Scrub the "Newton iteration" slider to watch the pose converge.
"""

import os

import numpy as np

from dyn import viz
from dyn.energy import StVK
from dyn.mesh import make_beam
from dyn.rig import LinearRig
from dyn.solver import make_static_problem, newton_minimize

P0R, P1R = (-1.3, 0.5), (-0.9, 0.9)          # plotted ranges for the landscape
GRID = 25
OFFSET = np.array([0.0, 3.0, 0.0])           # place the bowl beside the beam
SX, SY = 4.0, 3.0                            # landscape footprint


def _landscape(energy):
    p0s, p1s = np.linspace(*P0R, GRID), np.linspace(*P1R, GRID)
    Eg = np.array([[energy(np.array([a, b])) for b in p1s] for a in p0s])
    sz = 2.0 / (Eg.max() + 1e-9)
    u = (p0s - p0s[0]) / (p0s[-1] - p0s[0])
    w = (p1s - p1s[0]) / (p1s[-1] - p1s[0])
    V = np.zeros((GRID * GRID, 3))
    S = np.zeros(GRID * GRID)
    idx = lambda i, j: i * GRID + j           # noqa: E731
    for i in range(GRID):
        for j in range(GRID):
            V[idx(i, j)] = OFFSET + [u[i] * SX, w[j] * SY, Eg[i, j] * sz]
            S[idx(i, j)] = Eg[i, j]
    F = []
    for i in range(GRID - 1):
        for j in range(GRID - 1):
            a, b, c, d = idx(i, j), idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)
            F += [[a, b, c], [a, c, d]]
    return V, np.array(F), S, sz


def _to_landscape(p, E, sz):
    u = (p[0] - P0R[0]) / (P0R[1] - P0R[0])
    w = (p[1] - P1R[0]) / (P1R[1] - P1R[0])
    return OFFSET + np.array([u * SX, w * SY, E * sz])


def main(show=True):
    beam = make_beam()
    rig = LinearRig(beam, modes=("bend_z", "sway_y"))
    st = {"young": 200.0, "grav": 9.8}

    def solve():
        mat = StVK(beam, young=st["young"], poisson=0.40)
        energy, grad = make_static_problem(rig, mat, mass=beam.mass, gravity=(0, 0, -st["grav"]))
        pstar, hist = newton_minimize(energy, grad, np.zeros(2), record=True)
        return mat, energy, pstar, hist

    mat, energy, pstar, hist = solve()
    print("=== Lesson 5: statics = bottom of the energy landscape ===")
    print(f"converged in {len(hist) - 1} Newton iterations to p* = "
          f"{np.array2string(pstar, precision=4)}")
    print("  iter |        p            energy E(p)   |grad|")
    for k, (p, e, gn) in enumerate(hist):
        print(f"  {k:4d} | {np.array2string(p, precision=3):>16}  {e:11.4f}   {gn:.2e}")

    if not show:
        return pstar, hist

    import polyscope as ps
    import polyscope.imgui as psim

    # live beam morphs to each Newton iterate; rest pose shown as a translucent ghost
    v = viz.register_beam(beam)                  # this also initializes Polyscope
    ps.register_volume_mesh("rest (ghost)", beam.verts, tets=beam.tets,
                            color=(0.6, 0.6, 0.6), transparency=0.25, enabled=True)
    state = {"V": None, "F": None, "S": None, "sz": None, "hist": hist, "k": 0, "dirty": True}

    def rebuild():
        nonlocal energy, pstar
        mat, energy, pstar, hist = solve()
        V, F, S, sz = _landscape(energy)
        ps.register_surface_mesh("energy landscape", V, F, transparency=0.85,
                                 enabled=True).add_scalar_quantity(
            "E(p)", S, enabled=True, cmap="viridis")
        path = np.array([_to_landscape(p, e, sz) for (p, e, _) in hist])
        edges = np.array([[i, i + 1] for i in range(len(path) - 1)])
        ps.register_curve_network("Newton path", path, edges if len(edges) else "line",
                                  color=(1.0, 0.4, 0.0), radius=0.004)
        ps.register_point_cloud("iterates", path, color=(1.0, 0.4, 0.0), radius=0.01)
        state.update(V=V, F=F, S=S, sz=sz, hist=hist, dirty=False)
        state["k"] = min(state["k"], len(hist) - 1)

    def callback():
        psim.Text("Lesson 5 — statics: minimize  E(p) = W(s(p)) − Σ m_i g·s_i(p)")
        psim.Separator()
        ch1, st["young"] = psim.SliderFloat("Young's modulus E", float(st["young"]), 40.0, 800.0)
        ch2, st["grav"] = psim.SliderFloat("gravity |g|", float(st["grav"]), 0.0, 20.0)
        if ch1 or ch2 or state["dirty"]:
            rebuild()
        hist = state["hist"]
        kmax = len(hist) - 1
        _, state["k"] = psim.SliderInt("Newton iteration", int(state["k"]), 0, kmax)
        if psim.Button("re-solve from p = 0"):
            state["dirty"] = True

        p_k, e_k, gn_k = hist[state["k"]]
        v.update_vertex_positions(rig.s(p_k))
        # highlight the current iterate on the bowl
        cur = _to_landscape(p_k, e_k, state["sz"])[None, :]
        ps.register_point_cloud("current iterate", cur, color=(0.1, 0.9, 0.2), radius=0.02)
        psim.Separator()
        psim.Text(f"iterate {state['k']}/{kmax}:  p = {np.array2string(p_k, precision=3)}")
        psim.Text(f"energy E = {e_k:.4f}    |grad| = {gn_k:.2e}")
        psim.Text("The green dot rolls to the bottom of the bowl = the relaxed pose p*.")

    viz.show(callback)
    return pstar, hist


if __name__ == "__main__":
    main(show=not os.environ.get("DYN_NO_SHOW"))
