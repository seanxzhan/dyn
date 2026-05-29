"""Lesson 05 — Cloth meets a static sphere: inequality constraints + response.

    Run :  python src/dyn/tutorials/pbd/05_collisions_and_response.py
    Read:  src/dyn/tutorials/pbd/05_collisions_and_response.md

A pinned cloth swings under gravity into a static sphere positioned in its
swing arc.  Same predict/project loop as Lesson 4, plus one inequality
constraint (`project_sphere`) inside the iteration.  After projection,
friction/restitution shape the recovered velocity at collided particles.
"""

from __future__ import annotations

import os

import numpy as np

from dyn import viz
from dyn.pbd import (make_cloth, particle_inv_masses, project_distance,
                     project_bend, project_sphere, damp_velocities)


def main(show=True):
    cloth = make_cloth(res=(20, 20), size=(1.2, 1.2), pin="top",
                       normal=(0.0, 1.0, 0.0))
    w = particle_inv_masses(cloth, total_mass=1.0)
    # sphere sits in the cloth's swing path, slightly forward of the rest plane
    sphere_center = np.array([0.0, 0.5, 0.5])
    sphere_radius = 0.30

    print("=== Lesson 5: cloth + static sphere ===")
    print(f"verts={cloth.n}  tris={len(cloth.tris)}  edges={len(cloth.edges)}  "
          f"bend_pairs={len(cloth.bend_pairs)}  pinned={len(cloth.pin_idx)}")

    if not show:
        # numerical sanity: simulate, check no committed state ever penetrates
        x = cloth.verts.copy()
        v = np.zeros_like(x)
        h = 0.016
        thickness = 0.01
        worst_penetration = 0.0
        for step in range(200):
            v[:, 2] -= h * 9.8
            v[w == 0] = 0.0
            p = x + h * v
            for _ in range(10):
                project_distance(p, cloth.edges, cloth.rest_lengths, w, k_prime=1.0)
                project_bend(p, cloth.bend_pairs, cloth.rest_dihedrals, w, k_prime=0.4)
                project_sphere(p, w, sphere_center, sphere_radius, thickness=thickness)
            v = (p - x) / h
            x = p
            d = np.linalg.norm(x[w > 0] - sphere_center, axis=1).min()
            worst_penetration = min(worst_penetration, d - (sphere_radius + thickness))
        print(f"  after 200 steps: worst penetration = {worst_penetration:.2e}  "
              f"(must be ≥ -1e-6 — projection clamps to surface)")
        assert worst_penetration > -1e-4, f"penetration {worst_penetration} too large"
        return cloth

    import polyscope as ps
    import polyscope.imgui as psim

    viz.init_once()
    sm = ps.register_surface_mesh("cloth", cloth.verts, cloth.tris,
                                  color=(0.85, 0.55, 0.20), edge_width=0.5,
                                  smooth_shade=True)
    sm.set_back_face_policy("identical")
    sphere_v, sphere_t = _icosphere(sphere_center, sphere_radius, subdiv=2)
    sphere_mesh = ps.register_surface_mesh("sphere", sphere_v, sphere_t,
                                            color=(0.30, 0.55, 0.85),
                                            smooth_shade=True)
    pin_pc = ps.register_point_cloud("pinned", cloth.verts[cloth.pin_idx],
                                      color=(0.90, 0.20, 0.20), radius=0.012)
    h = 0.016
    state = {"play": True, "iters": 12, "k_stretch": 1.0, "k_bend": 0.4,
             "damp": 0.04, "g": 9.8, "mu": 0.4, "rest": 0.0,
             "thickness": 0.01, "n": 0,
             "x": cloth.verts.copy(), "v": np.zeros_like(cloth.verts)}

    def reset():
        state["x"] = cloth.verts.copy()
        state["v"][:] = 0.0
        state["n"] = 0

    def callback():
        psim.Text("Lesson 5 — cloth swings into a static sphere")
        psim.Separator()
        _, state["play"] = psim.Checkbox("Play", state["play"])
        if psim.Button("reset"):
            reset()
        _, state["iters"] = psim.SliderInt("iterations", int(state["iters"]), 1, 30)
        _, state["k_stretch"] = psim.SliderFloat("stretch stiffness k",
                                                  float(state["k_stretch"]), 0.1, 1.0)
        _, state["k_bend"] = psim.SliderFloat("bend stiffness k",
                                               float(state["k_bend"]), 0.0, 1.0)
        _, state["mu"] = psim.SliderFloat("friction μ (tangent)", float(state["mu"]), 0.0, 1.0)
        _, state["rest"] = psim.SliderFloat("restitution e (normal bounce)",
                                             float(state["rest"]), 0.0, 1.0)
        _, state["damp"] = psim.SliderFloat("internal damping", float(state["damp"]), 0.0, 0.5)

        if state["play"]:
            iters = max(int(state["iters"]), 1)
            kp_s = 1.0 - (1.0 - state["k_stretch"]) ** (1.0 / iters)
            kp_b = 1.0 - (1.0 - state["k_bend"]) ** (1.0 / iters) if state["k_bend"] > 0 else 0.0
            x, v = state["x"], state["v"]
            v_new = v.copy()
            v_new[:, 2] -= h * state["g"]
            v_new[w == 0] = 0.0
            p = x + h * v_new
            collided_any = np.zeros(cloth.n, dtype=bool)
            for _ in range(iters):
                project_distance(p, cloth.edges, cloth.rest_lengths, w, k_prime=kp_s)
                if kp_b > 0:
                    project_bend(p, cloth.bend_pairs, cloth.rest_dihedrals, w, k_prime=kp_b)
                hit = project_sphere(p, w, sphere_center, sphere_radius,
                                     thickness=state["thickness"])
                collided_any |= hit
            v_post = (p - x) / h
            if collided_any.any():
                idx = np.where(collided_any & (w > 0))[0]
                if len(idx):
                    n_hat = (p[idx] - sphere_center)
                    n_hat /= np.maximum(np.linalg.norm(n_hat, axis=1, keepdims=True), 1e-12)
                    vn = (v_post[idx] * n_hat).sum(axis=1, keepdims=True) * n_hat
                    vt = v_post[idx] - vn
                    v_post[idx] = (1.0 - state["mu"]) * vt - state["rest"] * vn
            m = np.where(w > 0, 1.0, np.inf)
            damp_velocities(p, v_post, m, state["damp"])
            state["v"] = v_post
            state["x"] = p
            state["n"] += 1

        sm.update_vertex_positions(state["x"])
        pin_pc.update_point_positions(state["x"][cloth.pin_idx])
        d = np.linalg.norm(state["x"] - sphere_center, axis=1)
        in_contact = (d < sphere_radius + state["thickness"] + 5e-3).sum()
        psim.Separator()
        psim.Text(f"step {state['n']}   particles in contact: {in_contact}   "
                  f"min |p-c|={d.min():.3f}   R+τ={sphere_radius + state['thickness']:.3f}")
        psim.Text("μ=0 slips off, μ=1 sticks; e>0 bounces (set e=0 for cloth).")

    viz.show(callback)
    return state


def _icosphere(center, radius, subdiv=2):
    """A small icosphere builder so we don't need a mesh file for the collider."""
    t = (1.0 + 5.0 ** 0.5) / 2.0
    v = np.array([
        [-1, t, 0], [1, t, 0], [-1, -t, 0], [1, -t, 0],
        [0, -1, t], [0, 1, t], [0, -1, -t], [0, 1, -t],
        [t, 0, -1], [t, 0, 1], [-t, 0, -1], [-t, 0, 1],
    ], dtype=float)
    f = np.array([
        [0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11],
        [1, 5, 9], [5, 11, 4], [11, 10, 2], [10, 7, 6], [7, 1, 8],
        [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9],
        [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1],
    ], dtype=np.int32)
    for _ in range(subdiv):
        edge_mid = {}
        new_f = []
        new_v = list(v)

        def mid(a, b):
            key = (min(a, b), max(a, b))
            if key in edge_mid:
                return edge_mid[key]
            new_v.append((new_v[a] + new_v[b]) * 0.5)
            edge_mid[key] = len(new_v) - 1
            return len(new_v) - 1

        for tri in f:
            a, b, c = tri
            ab, bc, ca = mid(a, b), mid(b, c), mid(c, a)
            new_f += [[a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]]
        v = np.array(new_v)
        f = np.array(new_f, dtype=np.int32)
    v = v / np.linalg.norm(v, axis=1, keepdims=True) * radius + np.asarray(center)
    return v, f


if __name__ == "__main__":
    main(show=not os.environ.get("DYN_NO_SHOW"))
