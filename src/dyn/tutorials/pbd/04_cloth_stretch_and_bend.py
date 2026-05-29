"""Lesson 04 — Cloth: stretch + dihedral bend, the same PBD loop on a 2D mesh.

    Run :  python src/dyn/tutorials/pbd/04_cloth_stretch_and_bend.py
    Read:  src/dyn/tutorials/pbd/04_cloth_stretch_and_bend.md

A grid of particles, top row pinned, hangs under gravity.  Two constraint
types — distance per edge, dihedral angle per shared-edge triangle pair —
projected sequentially each step.  Slide bend stiffness from 0 (limp paper
towel) to 1 (cardstock) to feel why a length-independent bending term is the
paper's key cloth contribution.
"""

from __future__ import annotations

import os

import numpy as np

from dyn import viz
from dyn.pbd import (make_cloth, particle_inv_masses, predict, project_distance,
                     project_bend, damp_velocities, verify_bend_gradient)


def main(show=True):
    cloth = make_cloth(res=(20, 15), size=(1.5, 1.0), pin="top",
                       normal=(0.0, 1.0, 0.0))
    w = particle_inv_masses(cloth, total_mass=1.0)
    print("=== Lesson 4: cloth (stretch + dihedral bend) ===")
    print(f"verts={cloth.n}  tris={len(cloth.tris)}  "
          f"edges={len(cloth.edges)}  bend_pairs={len(cloth.bend_pairs)}  "
          f"pinned={len(cloth.pin_idx)}")
    info = verify_bend_gradient()
    print(f"bend gradient FD check: max |dp_analytic - dp_numeric| "
          f"= {info['max_dp_err']:.2e}  (must be ~1e-10)")
    assert info["max_dp_err"] < 1e-6

    if not show:
        return cloth

    import polyscope as ps
    import polyscope.imgui as psim

    viz.init_once()
    sm = ps.register_surface_mesh("cloth", cloth.verts, cloth.tris,
                                  color=(0.85, 0.55, 0.20), edge_width=0.5,
                                  smooth_shade=True)
    sm.set_back_face_policy("identical")
    pin_pc = ps.register_point_cloud("pinned", cloth.verts[cloth.pin_idx],
                                     color=(0.90, 0.20, 0.20), radius=0.012)
    h = 0.016                                            # ~60 fps
    state = {"play": True, "iters": 12, "k_stretch": 1.0, "k_bend": 0.4,
             "damp": 0.04, "g": 9.8, "wind": 0.0, "n": 0,
             "x": cloth.verts.copy(), "v": np.zeros_like(cloth.verts)}

    def reset():
        state["x"] = cloth.verts.copy()
        state["v"][:] = 0.0
        state["n"] = 0

    def callback():
        psim.Text("Lesson 4 — cloth = particles + stretch edges + dihedral bend pairs")
        psim.Separator()
        _, state["play"] = psim.Checkbox("Play", state["play"])
        if psim.Button("reset"):
            reset()
        psim.SameLine()
        if psim.Button("wind kick"):
            state["wind"] = 1.0
        _, state["iters"] = psim.SliderInt("iterations", int(state["iters"]), 1, 30)
        _, state["k_stretch"] = psim.SliderFloat("stretch stiffness k",
                                                  float(state["k_stretch"]), 0.1, 1.0)
        _, state["k_bend"] = psim.SliderFloat("bend stiffness k",
                                               float(state["k_bend"]), 0.0, 1.0)
        _, state["damp"] = psim.SliderFloat("damping", float(state["damp"]), 0.0, 0.5)
        _, state["g"] = psim.SliderFloat("|gravity|", float(state["g"]), 0.0, 20.0)

        if state["play"]:
            iters = max(int(state["iters"]), 1)
            kp_s = 1.0 - (1.0 - state["k_stretch"]) ** (1.0 / iters)
            kp_b = 1.0 - (1.0 - state["k_bend"]) ** (1.0 / iters) if state["k_bend"] > 0 else 0.0
            # one-step external "wind" impulse on the bottom-half particles
            f_ext = np.zeros_like(state["x"])
            f_ext[:, 2] = -state["g"]                    # gravity
            if state["wind"] > 0:
                lower = state["x"][:, 2] < state["x"][:, 2].mean()
                f_ext[lower, 1] += 80.0 * state["wind"]
                state["wind"] = 0.0
            # predict by hand because we have a per-particle f_ext (with wind)
            v_new = state["v"] + h * (f_ext * (w > 0)[:, None])
            p = state["x"] + h * v_new
            for _ in range(iters):
                project_distance(p, cloth.edges, cloth.rest_lengths, w, k_prime=kp_s)
                if kp_b > 0:
                    project_bend(p, cloth.bend_pairs, cloth.rest_dihedrals, w, k_prime=kp_b)
            state["v"] = (p - state["x"]) / h
            # rigid-mode-preserving damping
            m = np.where(w > 0, 1.0, np.inf)
            damp_velocities(p, state["v"], m, state["damp"])
            state["x"] = p
            state["n"] += 1

        sm.update_vertex_positions(state["x"])
        pin_pc.update_point_positions(state["x"][cloth.pin_idx])
        # show stretch ratio so you can see softness/decoupling visually
        seg = np.linalg.norm(state["x"][cloth.edges[:, 0]]
                              - state["x"][cloth.edges[:, 1]], axis=1)
        psim.Separator()
        psim.Text(f"step {state['n']}   mean stretch ratio = "
                  f"{(seg / cloth.rest_lengths).mean():.3f}   "
                  f"(1.0 = inextensible)")
        psim.Text("k_bend → 0 collapses the cloth into folds; k_bend → 1 makes it cardboard.")

    viz.show(callback)
    return state


if __name__ == "__main__":
    main(show=not os.environ.get("DYN_NO_SHOW"))
