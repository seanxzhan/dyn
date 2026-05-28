"""Lesson 08 (capstone) — Secondary motion on a primary animation.

    Run :  python src/dyn/tutorials/08_capstone_primary_secondary.py
    Read:  src/dyn/tutorials/08_capstone_primary_secondary.md

The production setting: an animator keyframes some rig parameters (PRIMARY) and
physics adds the rest (SECONDARY).  Here a "lift" parameter is driven up and down
(the keyframe), while a free "bend" parameter is simulated.  Inertia makes the
beam lag and overshoot the handle — secondary motion — for free, in rig space.
Sweep the drive frequency through the beam's resonance to see it amplify.
"""

import os

import numpy as np

from dyn import viz
from dyn.energy import StVK
from dyn.mesh import make_beam
from dyn.rig import LinearRig
from dyn.solver import ImplicitEuler

TRAIL = 240


def main(show=True):
    beam = make_beam()
    # p[0] = lift_z (DRIVEN primary, moves whole beam), p[1] = bend_z (FREE secondary)
    rig = LinearRig(beam, modes=("lift_z", "bend_z"))
    mat = StVK(beam, young=200.0, poisson=0.40)
    _c = np.where(np.isclose(beam.verts[:, 0], beam.verts[:, 0].max()))[0]
    tip = int(_c[np.argmin(((beam.verts[_c][:, 1:] - beam.verts[_c][:, 1:].mean(0)) ** 2).sum(1))])

    print("=== Lesson 8: primary (driven lift) + secondary (simulated bend) ===")
    sim = ImplicitEuler(rig, mat, beam.mass, h=0.02, gravity=(0, 0, 0.0),
                        damping=0.06, free_idx=[1])
    A, w, h = 0.4, 5.0, 0.02
    peak = 0.0
    for n in range(400):
        sim.step(driven_values=[A * np.sin(w * n * h)])
        peak = max(peak, abs(sim.p[1]))
    print(f"drive A={A}, ω={w}: peak secondary bend |p[1]| = {peak:.3f}  "
          f"(0 would mean no secondary motion)")

    if not show:
        return sim

    import polyscope as ps
    import polyscope.imgui as psim

    sim = ImplicitEuler(rig, mat, beam.mass, h=0.02, gravity=(0, 0, 0.0),
                        damping=0.06, free_idx=[1])
    v = viz.register_beam(beam)
    ghost = ps.register_volume_mesh("primary only (no physics)", beam.verts, tets=beam.tets,
                                    color=(0.6, 0.6, 0.6), transparency=0.3)
    st = {"play": True, "physics": True, "A": 0.4, "w": 5.0, "damp": 0.06,
          "young": 200.0, "grav": 0.0, "n": 0, "trail": []}

    def reset():
        nonlocal sim
        sim = ImplicitEuler(rig, mat, beam.mass, h=0.02, gravity=(0, 0, -st["grav"]),
                            damping=st["damp"], free_idx=[1])
        st["n"], st["trail"] = 0, []

    def callback():
        psim.Text("Lesson 8 — PRIMARY (driven lift) + SECONDARY (simulated bend)")
        psim.Separator()
        _, st["play"] = psim.Checkbox("Play", st["play"])
        psim.SameLine()
        _, st["physics"] = psim.Checkbox("secondary physics on", st["physics"])
        if psim.Button("reset"):
            reset()
        _, st["A"] = psim.SliderFloat("primary amplitude A", float(st["A"]), 0.0, 1.0)
        _, st["w"] = psim.SliderFloat("primary frequency ω (sweep for resonance)", float(st["w"]), 0.5, 18.0)
        _, st["damp"] = psim.SliderFloat("damping", float(st["damp"]), 0.0, 0.3)
        ch_y, st["young"] = psim.SliderFloat("Young's modulus E", float(st["young"]), 40.0, 800.0)
        _, st["grav"] = psim.SliderFloat("gravity |g|", float(st["grav"]), 0.0, 15.0)

        sim.elastic = StVK(beam, young=st["young"], poisson=0.40) if ch_y else sim.elastic
        sim.fg = sim.M * np.array([0.0, 0.0, -st["grav"]])
        sim.damping = st["damp"]

        t = st["n"] * sim.h
        lift = st["A"] * np.sin(st["w"] * t)
        if st["play"]:
            if st["physics"]:
                sim.step(driven_values=[lift])
            else:
                sim.p[:] = [lift, 0.0]            # primary only, bend frozen at 0
                sim.x = rig.s(sim.p)
                sim.x_prev = sim.x.copy()
            st["n"] += 1
            st["trail"].append(sim.x[tip].copy())
            st["trail"] = st["trail"][-TRAIL:]

        v.update_vertex_positions(sim.x)
        ghost.update_vertex_positions(rig.s([sim.p[0], 0.0]))   # where the rig alone would be
        if len(st["trail"]) > 2:
            ps.register_curve_network("tip trajectory", np.array(st["trail"]), "line",
                                      color=(0.95, 0.55, 0.1), radius=0.003)
        psim.Separator()
        psim.Text(f"primary lift p[0] = {sim.p[0]:+.3f}   secondary bend p[1] = {sim.p[1]:+.3f}")
        psim.Text("Grey ghost = rig with bend=0 (no physics). The gap is the secondary motion.")
        psim.Text("Sweep ω toward the beam's natural frequency to amplify the lag/overshoot.")

    viz.show(callback)
    return sim


if __name__ == "__main__":
    main(show=not os.environ.get("DYN_NO_SHOW"))
