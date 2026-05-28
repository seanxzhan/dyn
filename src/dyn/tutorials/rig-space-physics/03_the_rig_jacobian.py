"""Lesson 03 — The rig Jacobian J = ds/dp and the projection f_p = Jᵀ f_s.

    Run :  python src/dyn/tutorials/03_the_rig_jacobian.py
    Read:  src/dyn/tutorials/03_the_rig_jacobian.md

Each column of J is "which way every vertex moves when you nudge p[i]" — drawn as
a per-vertex arrow field.  For the LinearRig the arrows are CONSTANT; switch to
the nonlinear BendRig and they visibly change as you bend (rig curvature, the
2012 O(p²) term).  We also verify J against finite differences, and show how a
tip force f_s projects to a rig-space force f_p = Jᵀ f_s.
"""

import os

import numpy as np

from dyn import viz
from dyn.mesh import make_beam
from dyn.rig import BendRig, LinearRig, finite_difference_jacobian

_COLORS = [(0.95, 0.55, 0.15), (0.20, 0.55, 0.95), (0.45, 0.80, 0.30)]


def tip_vertex(beam):
    """Index of a vertex at the free tip (max x, centered in y,z)."""
    cand = np.where(np.isclose(beam.verts[:, 0], beam.verts[:, 0].max()))[0]
    mid = beam.verts[cand][:, 1:].mean(0)
    return cand[np.argmin(((beam.verts[cand][:, 1:] - mid) ** 2).sum(1))]


def _check_jacobians(lin, bend):
    print("=== Lesson 3: J = ds/dp, verified against finite differences ===")
    for p in (np.array([0.4, -0.3]),):
        err = np.abs(lin.jacobian(p) - finite_difference_jacobian(lin, p)).max()
        print(f"LinearRig p={p}: max|J_analytic - J_fd| = {err:.2e}  (J is constant)")
        assert err < 1e-6
    for th in (0.0, 0.7):
        p = np.array([th])
        err = np.abs(bend.jacobian(p) - finite_difference_jacobian(bend, p)).max()
        curv = np.linalg.norm(bend.second_derivative(p))
        print(f"BendRig θ={th:+.2f}: max|J_an - J_fd| = {err:.2e}   "
              f"||d²s/dθ²|| = {curv:.2f}")
        assert err < 1e-5
    print("LinearRig curvature d²s/dp² = 0  ->  J constant  (the 2013 regime).")
    print("BendRig curvature ≠ 0      ->  J varies with θ   (the 2012 O(p²) term).\n")


def main(show=True):
    beam = make_beam()
    lin = LinearRig(beam, modes=("bend_z", "sway_y"))
    bend = BendRig(beam)
    _check_jacobians(lin, bend)

    if not show:
        return lin, bend

    import polyscope.imgui as psim

    v = viz.register_beam(beam)
    tip = tip_vertex(beam)
    state = {"use_bend": False, "p_lin": np.zeros(lin.dim), "p_bend": np.zeros(1),
             "show_J": True, "show_force": False, "fz": -2.0}
    drawn: set[str] = set()

    def clear():
        for name in list(drawn):
            try:
                v.remove_quantity(name)
            except Exception:
                pass
        drawn.clear()

    def draw_vec(name, vecs, color, length=0.5):
        viz.add_vectors(v, name, vecs, color=color, length=length)
        drawn.add(name)

    def callback():
        psim.Text("Lesson 3 — the rig Jacobian J = ds/dp")
        psim.Separator()
        _, state["use_bend"] = psim.Checkbox("use BendRig (nonlinear: J varies)", state["use_bend"])
        rig = bend if state["use_bend"] else lin
        key = "p_bend" if state["use_bend"] else "p_lin"
        p = state[key]

        for i in range(rig.dim):
            _, p[i] = psim.SliderFloat(f"p[{i}]  ({rig.names[i]})", float(p[i]), -1.5, 1.5)
        _, state["show_J"] = psim.Checkbox("show Jacobian columns (arrows)", state["show_J"])
        _, state["show_force"] = psim.Checkbox("show force projection  f_p = Jᵀ f_s", state["show_force"])
        _, state["fz"] = psim.SliderFloat("tip force f_s (z component)", float(state["fz"]), -4.0, 4.0)

        # Re-evaluate the rig and its Jacobian at the current parameters.
        v.update_vertex_positions(rig.s(p))
        J = rig.jacobian(p)                        # (3N, dim)
        clear()

        if state["show_J"]:
            for i in range(rig.dim):
                col = J[:, i].reshape(-1, 3)        # i-th column as a vertex field
                draw_vec(f"J column {i} ({rig.names[i]})", col, _COLORS[i % len(_COLORS)])

        if state["show_force"]:
            f_s = np.zeros((beam.n, 3))
            f_s[tip, 2] = state["fz"]               # a point force at the tip, in z
            f_p = J.T @ f_s.reshape(-1)             # generalized force in rig space (dim,)
            response = (J @ f_p).reshape(-1, 3)     # the motion the rig takes under it
            draw_vec("applied force f_s (tip)", f_s, (0.90, 0.10, 0.10), length=0.4)
            draw_vec("rig response  J Jᵀ f_s", response, (0.10, 0.70, 0.90), length=0.5)
            psim.Text(f"f_p = Jᵀ f_s = {np.array2string(f_p, precision=3, suppress_small=True)}")

        psim.Separator()
        psim.Text("LinearRig: arrows are fixed.  BendRig: bend the beam -> arrows turn.")

    viz.show(callback)
    return lin, bend


if __name__ == "__main__":
    main(show=not os.environ.get("DYN_NO_SHOW"))
