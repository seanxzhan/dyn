"""Lesson 07 — Nonlinear rig: the 2012 curvature term vs. the 2013 linearization.

    Run :  python src/dyn/tutorials/07_nonlinear_rig_2012_vs_2013.py
    Read:  src/dyn/tutorials/07_nonlinear_rig_2012_vs_2013.md

Two dynamics solves of the SAME nonlinear BendRig run together:
  * "2012 full"  uses the true rig each Newton step (curvature ∂²s/∂θ² included),
  * "2013 linear" freezes s ≈ s(θ_n) + J·Δθ each step (curvature dropped).
They stay within ~1% of beam height — the 2013 paper's "imperceptible" claim —
while the linear version needs no O(d²) rig-curvature evaluations.  The magnified
difference field shows where (little) disagreement lives.
"""

import os

import numpy as np

from dyn import viz
from dyn.energy import StVK
from dyn.mesh import make_beam
from dyn.rig import BendRig
from dyn.solver import ImplicitEuler

MAG = 25.0  # magnification for the difference field


def main(show=True):
    beam = make_beam()
    mat = StVK(beam, young=200.0, poisson=0.40)
    H = beam.lengths[2]

    def make_pair(grav, damp):
        full = ImplicitEuler(BendRig(beam), mat, beam.mass, h=0.02,
                             gravity=(0, 0, -grav), damping=damp, relinearize=False)
        lin = ImplicitEuler(BendRig(beam), mat, beam.mass, h=0.02,
                            gravity=(0, 0, -grav), damping=damp, relinearize=True)
        return full, lin

    print("=== Lesson 7: 2012 (full curvature) vs 2013 (linearized rig) ===")
    full, lin = make_pair(9.8, 0.04)
    worst = 0.0
    for _ in range(150):
        full.step()
        lin.step()
        worst = max(worst, np.abs(full.x - lin.x).max())
    print(f"over 150 steps: max vertex difference = {worst:.3e} "
          f"({100 * worst / H:.3f}% of beam height)")
    print(f"final bend angle: full θ = {full.p[0]:+.4f}   linearized θ = {lin.p[0]:+.4f}")
    print("=> linearizing the rig barely changes the motion, but removes the")
    print(f"   curvature term ||d²s/dθ²|| = {np.linalg.norm(BendRig(beam).second_derivative(full.p)):.2f}")

    if not show:
        return full, lin

    import polyscope as ps
    import polyscope.imgui as psim

    full, lin = make_pair(9.8, 0.04)
    v = viz.register_beam(beam, name="2012 full (solid)")
    ghost = ps.register_volume_mesh("2013 linearized (ghost)", beam.verts, tets=beam.tets,
                                    color=(0.3, 0.85, 0.4), transparency=0.45)
    st = {"play": True, "grav": 9.8, "damp": 0.04, "worst": 0.0, "n": 0}

    def callback():
        psim.Text("Lesson 7 — does Newton need the rig curvature ∂²s/∂θ²?")
        psim.Separator()
        _, st["play"] = psim.Checkbox("Play", st["play"])
        if psim.Button("reset"):
            f2, l2 = make_pair(st["grav"], st["damp"])
            full.__dict__.update(f2.__dict__)
            lin.__dict__.update(l2.__dict__)
            st["worst"] = 0.0
            st["n"] = 0
        if psim.Button("pluck"):
            kick = np.zeros((beam.n, 3))
            kick[:, 2] = 8.0 * (beam.verts[:, 0] / beam.lengths[0]) ** 2
            full.pluck((full.x - full.x_prev) / full.h + kick)
            lin.pluck((lin.x - lin.x_prev) / lin.h + kick)
        _, st["grav"] = psim.SliderFloat("gravity |g|", float(st["grav"]), 0.0, 25.0)
        _, st["damp"] = psim.SliderFloat("damping", float(st["damp"]), 0.0, 0.3)
        for s in (full, lin):
            s.fg = s.M * np.array([0.0, 0.0, -st["grav"]])
            s.damping = st["damp"]

        if st["play"]:
            full.step()
            lin.step()
            st["n"] += 1
            st["worst"] = max(st["worst"], float(np.abs(full.x - lin.x).max()))

        v.update_vertex_positions(full.x)
        ghost.update_vertex_positions(lin.x)
        viz.add_vectors(v, f"difference x{MAG:.0f} (full − linearized)",
                        (full.x - lin.x) * MAG, color=(0.9, 0.1, 0.5), length=0.5)

        curv = np.linalg.norm(BendRig(beam).second_derivative(full.p))
        psim.Separator()
        psim.Text(f"step {st['n']}   θ_full = {full.p[0]:+.4f}   θ_lin = {lin.p[0]:+.4f}")
        psim.Text(f"max vertex diff so far = {st['worst']:.3e}  "
                  f"({100 * st['worst'] / beam.lengths[2]:.3f}% of beam height)")
        psim.Text(f"rig curvature ||∂²s/∂θ²|| = {curv:.2f}  "
                  f"(the 2012 O(d²) term the 2013 method drops)")

    viz.show(callback)
    return full, lin


if __name__ == "__main__":
    main(show=not os.environ.get("DYN_NO_SHOW"))
