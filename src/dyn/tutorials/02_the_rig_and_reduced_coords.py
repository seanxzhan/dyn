"""Lesson 02 — The rig and reduced coordinates.

    Run :  python src/dyn/tutorials/02_the_rig_and_reduced_coords.py
    Read:  src/dyn/tutorials/02_the_rig_and_reduced_coords.md

A LinearRig maps a few parameters p to every vertex via s(p) = s0 + B p.  Drag
the sliders: a handful of numbers drive all 3N vertex DOFs.  That low-dimensional
set {s(p)} is "rig space" — the subspace the physics will live in.
"""

import os

import numpy as np

from dyn import viz
from dyn.mesh import make_beam
from dyn.rig import LinearRig


def main(show=True):
    beam = make_beam()
    rig = LinearRig(beam, modes=("bend_z", "sway_y"))

    print("=== Lesson 2: the rig as a low-dimensional map ===")
    print(f"rig parameters p : {rig.names}  (dim = {rig.dim})")
    print(f"full vertex DOFs : {beam.n * 3}")
    print(f"compression      : {beam.n * 3} DOFs driven by {rig.dim} numbers "
          f"(~{beam.n * 3 // rig.dim}x)")
    # The rig must reproduce the rest pose at p = 0.
    assert np.allclose(rig.s(np.zeros(rig.dim)), beam.verts)

    if not show:
        return rig

    import polyscope.imgui as psim

    v = viz.register_beam(beam)
    state = {"p": rig.rest_params.copy()}
    rng = 1.5  # slider range for each parameter

    def callback():
        psim.Text("Lesson 2 — a rig is a map  p  ->  s(p) = s0 + B p")
        psim.Text(f"{rig.dim} parameters drive all {beam.n * 3} vertex coordinates.")
        psim.Separator()
        changed = False
        for i in range(rig.dim):
            ch, state["p"][i] = psim.SliderFloat(
                f"p[{i}]  ({rig.names[i]})", float(state["p"][i]), -rng, rng)
            changed = changed or ch
        if psim.Button("reset to rest pose (p = 0)"):
            state["p"][:] = rig.rest_params
            changed = True
        if changed:
            v.update_vertex_positions(rig.s(state["p"]))   # re-evaluate the rig
        psim.Separator()
        psim.Text(f"p = {np.array2string(state['p'], precision=2, suppress_small=True)}")
        psim.Text("Every reachable shape lies in 'rig space' = span(B) around s0.")

    viz.show(callback)
    return rig


if __name__ == "__main__":
    main(show=not os.environ.get("DYN_NO_SHOW"))
