"""Lesson 06 — Dynamics: implicit Euler as energy minimization (secondary motion).

    Run :  python src/dyn/tutorials/06_dynamics_implicit_euler.py
    Read:  src/dyn/tutorials/06_dynamics_implicit_euler.md

Press Play: released from rest, the beam falls under gravity, overshoots, and
jiggles — that damped transient around the rig pose IS secondary motion.  It
settles onto the static equilibrium from Lesson 5 (shown as a ghost).  Pluck it,
change damping/stiffness, or turn physics off to compare.
"""

import os

import numpy as np

from dyn import viz
from dyn.energy import StVK
from dyn.mesh import make_beam
from dyn.rig import LinearRig
from dyn.solver import ImplicitEuler, make_static_problem, newton_minimize

TRAIL = 240


def main(show=True):
    beam = make_beam()
    rig = LinearRig(beam, modes=("bend_z", "sway_y"))
    mat = StVK(beam, young=200.0, poisson=0.40)
    Lx = beam.lengths[0]
    _cand = np.where(np.isclose(beam.verts[:, 0], beam.verts[:, 0].max()))[0]
    tip = int(_cand[np.argmin(((beam.verts[_cand][:, 1:] - beam.verts[_cand][:, 1:].mean(0)) ** 2).sum(1))])

    def equilibrium(young, grav):
        m = StVK(beam, young=young, poisson=0.40)
        en, gr = make_static_problem(rig, m, mass=beam.mass, gravity=(0, 0, -grav))
        return newton_minimize(en, gr, np.zeros(2))[0]

    print("=== Lesson 6: implicit-Euler dynamics ===")
    sim = ImplicitEuler(rig, mat, beam.mass, h=0.02, gravity=(0, 0, -9.8), damping=0.04)
    pstar = equilibrium(200.0, 9.8)
    for _ in range(200):
        sim.step()
    print(f"released from rest -> settled at p = {np.array2string(sim.p, precision=3)}")
    print(f"static equilibrium  p* = {np.array2string(pstar, precision=3)}  (they match)")

    if not show:
        return sim

    import polyscope as ps
    import polyscope.imgui as psim

    sim = ImplicitEuler(rig, mat, beam.mass, h=0.02, gravity=(0, 0, -9.8), damping=0.04)
    v = viz.register_beam(beam)
    ghost = ps.register_volume_mesh("static equilibrium (ghost)", rig.s(pstar), tets=beam.tets,
                                    color=(0.3, 0.8, 0.4), transparency=0.25)
    st = {"play": True, "physics": True, "young": 200.0, "grav": 9.8,
          "damp": 0.04, "n": 0, "trail": []}

    def reset():
        nonlocal sim
        sim = ImplicitEuler(rig, mat, beam.mass, h=0.02,
                            gravity=(0, 0, -st["grav"]), damping=st["damp"])
        st["n"] = 0
        st["trail"] = []

    def callback():
        nonlocal pstar
        psim.Text("Lesson 6 — minimize Φ(p) = inertia + W − gravity, once per step")
        psim.Separator()
        _, st["play"] = psim.Checkbox("Play", st["play"])
        psim.SameLine()
        _, st["physics"] = psim.Checkbox("physics on (off = rig stays at rest)", st["physics"])
        if psim.Button("reset"):
            reset()
        psim.SameLine()
        if psim.Button("pluck (flick the tip up)"):
            kick = np.zeros((beam.n, 3))
            kick[:, 2] = 7.0 * (beam.verts[:, 0] / Lx) ** 2
            sim.pluck((sim.x - sim.x_prev) / sim.h + kick)

        ch_g, st["grav"] = psim.SliderFloat("gravity |g|", float(st["grav"]), 0.0, 20.0)
        ch_y, st["young"] = psim.SliderFloat("Young's modulus E", float(st["young"]), 40.0, 800.0)
        _, st["damp"] = psim.SliderFloat("damping", float(st["damp"]), 0.0, 0.3)

        # live parameter updates (no reset): material, gravity, damping
        sim.elastic = StVK(beam, young=st["young"], poisson=0.40)
        sim.fg = sim.M * np.array([0.0, 0.0, -st["grav"]])
        sim.damping = st["damp"]
        if ch_g or ch_y:
            pstar = equilibrium(st["young"], st["grav"])
            ghost.update_vertex_positions(rig.s(pstar))

        if st["play"] and st["physics"]:
            sim.step()
            st["n"] += 1
            st["trail"].append(sim.x[tip].copy())
            st["trail"] = st["trail"][-TRAIL:]
        if not st["physics"]:
            sim.p[:] = rig.rest_params
            sim.x = rig.s(sim.p)
            sim.x_prev = sim.x.copy()

        v.update_vertex_positions(sim.x)
        if len(st["trail"]) > 2:
            pts = np.array(st["trail"])
            ps.register_curve_network("tip trajectory", pts, "line",
                                      color=(0.95, 0.55, 0.1), radius=0.003)
        psim.Separator()
        psim.Text(f"step {st['n']}   KE = {sim.kinetic_energy():.4f}   "
                  f"p = {np.array2string(sim.p, precision=3)}")
        psim.Text("Overshoot then settle onto the green ghost = damped secondary motion.")

    viz.show(callback)
    return sim


if __name__ == "__main__":
    main(show=not os.environ.get("DYN_NO_SHOW"))
