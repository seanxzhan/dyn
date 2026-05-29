"""Lesson 06 — Capstone: a ball thrown at a hanging cloth (two-way coupling).

    Run :  python src/dyn/tutorials/pbd/06_capstone_ball_through_cloth.py
    Read:  src/dyn/tutorials/pbd/06_capstone_ball_through_cloth.md

A pinned cloth hangs as a curtain.  A dynamic sphere — mass, gravity,
velocity — falls onto it or is launched into it.  The ball pushes the cloth;
the cloth pushes back via accumulated impulses (`project_dynamic_sphere`).
Sub-stepping is exposed as a slider so you can observe the tunneling /
catching transition without changing the visible frame rate.
"""

from __future__ import annotations

import os

import numpy as np

from dyn import viz
from dyn.pbd import (make_cloth, particle_inv_masses, project_distance,
                     project_bend, project_dynamic_sphere, damp_velocities,
                     DynamicSphere)


def main(show=True):
    cloth = make_cloth(res=(24, 24), size=(1.6, 1.6), pin="top",
                       normal=(0.0, 1.0, 0.0))
    w = particle_inv_masses(cloth, total_mass=1.0)
    cloth_mass = np.full(cloth.n, 1.0 / cloth.n)
    sphere = DynamicSphere(center=np.array([0.0, -1.5, 0.5]), radius=0.18,
                            mass=0.5,
                            velocity=np.zeros(3))

    print("=== Lesson 6: capstone — ball through cloth ===")
    print(f"verts={cloth.n}  tris={len(cloth.tris)}  edges={len(cloth.edges)}  "
          f"bend_pairs={len(cloth.bend_pairs)}  pinned={len(cloth.pin_idx)}")
    print(f"sphere mass={sphere.mass}  R={sphere.radius}")

    if not show:
        # numerical sanity: throw the ball at the cloth's middle, check it
        # transfers momentum (cloth slows it, sphere keeps moving but slower).
        x = cloth.verts.copy()
        v = np.zeros_like(x)
        # let cloth settle for a few steps under gravity first (with damping
        # so its in-plane wobble dies before we throw the ball)
        h = 0.016
        n_iters = 10
        thickness = 0.01
        m_arr = np.where(w > 0, 1.0, np.inf)
        for _ in range(40):
            v[:, 2] -= h * 9.8
            v[w == 0] = 0.0
            p = x + h * v
            for _ in range(n_iters):
                project_distance(p, cloth.edges, cloth.rest_lengths, w, k_prime=1.0)
                project_bend(p, cloth.bend_pairs, cloth.rest_dihedrals, w, k_prime=0.4)
            v = (p - x) / h
            damp_velocities(p, v, m_arr, 0.3)
            x = p
        # aim at the cloth's middle (~z=0.8) from y=-0.6
        sphere.center = np.array([0.0, -0.6, 0.8])
        sphere.velocity = np.array([0.0, 6.0, 0.0])
        sphere.impulse = np.zeros(3)
        v0 = sphere.velocity.copy()
        max_pen = 0.0
        # 3 sub-steps per frame: at v=6 m/s and grid spacing ~7 cm, h_sub=0.005
        # keeps the ball moving < half a cell between collision checks.
        n_sub = 3
        h_sub = h / n_sub
        for step in range(60):
            for _ in range(n_sub):
                v[:, 2] -= h_sub * 9.8
                v[w == 0] = 0.0
                p = x + h_sub * v
                P = sphere.predict(h_sub)
                sphere.center = P
                for _ in range(n_iters):
                    project_distance(p, cloth.edges, cloth.rest_lengths, w, k_prime=1.0)
                    project_bend(p, cloth.bend_pairs, cloth.rest_dihedrals, w, k_prime=0.4)
                    project_dynamic_sphere(p, x, w, cloth_mass, sphere, h_sub,
                                           thickness=thickness)
                v_post = (p - x) / h_sub
                v = v_post
                x = p
                sphere.commit(P, h_sub)
            d = np.linalg.norm(x[w > 0] - sphere.center, axis=1).min()
            max_pen = min(max_pen, d - (sphere.radius + thickness))
        speed_y_drop = v0[1] - sphere.velocity[1]
        print(f"  ball v_y: started {v0[1]:.2f}, ended {sphere.velocity[1]:.2f}  "
              f"(slowed by {speed_y_drop:.2f}) — cloth absorbed momentum")
        print(f"  worst penetration during throw: {max_pen:.2e}  (must be ≥ -1e-4)")
        assert speed_y_drop > 0.5, "ball didn't lose meaningful y-momentum to the cloth"
        assert max_pen > -1e-3, f"penetration {max_pen} too large even with sub-stepping"
        return cloth, sphere

    import polyscope as ps
    import polyscope.imgui as psim

    viz.init_once()
    sm = ps.register_surface_mesh("cloth", cloth.verts, cloth.tris,
                                  color=(0.85, 0.55, 0.20), edge_width=0.5,
                                  smooth_shade=True)
    sm.set_back_face_policy("identical")
    sphere_v_local, sphere_t = _icosphere(np.zeros(3), sphere.radius, subdiv=2)
    sphere_mesh = ps.register_surface_mesh("ball", sphere_v_local + sphere.center,
                                            sphere_t, color=(0.30, 0.55, 0.85),
                                            smooth_shade=True)
    pin_pc = ps.register_point_cloud("pinned", cloth.verts[cloth.pin_idx],
                                      color=(0.90, 0.20, 0.20), radius=0.012)

    h = 0.016
    state = {"play": True, "iters": 12, "substeps": 2,
             "k_stretch": 1.0, "k_bend": 0.4, "damp": 0.04,
             "g": 9.8, "thickness": 0.01, "throw_speed": 5.0,
             "ball_mass": 0.5, "ball_locked": True, "n": 0,
             "x": cloth.verts.copy(), "v": np.zeros_like(cloth.verts)}

    def reset(also_ball=True):
        state["x"] = cloth.verts.copy()
        state["v"][:] = 0.0
        state["n"] = 0
        if also_ball:
            sphere.center = np.array([0.0, -1.5, 0.5])
            sphere.velocity = np.zeros(3)
            sphere.impulse = np.zeros(3)
            state["ball_locked"] = True

    def callback():
        psim.Text("Lesson 6 (capstone) — ball through cloth, two-way coupled")
        psim.Separator()
        _, state["play"] = psim.Checkbox("Play", state["play"])
        if psim.Button("reset (cloth + ball)"):
            reset(also_ball=True)
        psim.SameLine()
        if psim.Button("release ball (drop)"):
            sphere.center = np.array([0.0, 0.0, 1.4])
            sphere.velocity = np.zeros(3)
            sphere.impulse = np.zeros(3)
            state["ball_locked"] = False
        psim.SameLine()
        if psim.Button("throw ball"):
            sphere.center = np.array([0.0, -1.5, 0.5])
            sphere.velocity = np.array([0.0, state["throw_speed"], 0.0])
            sphere.impulse = np.zeros(3)
            state["ball_locked"] = False

        _, state["iters"] = psim.SliderInt("iterations", int(state["iters"]), 1, 30)
        _, state["substeps"] = psim.SliderInt("sub-steps per frame",
                                                int(state["substeps"]), 1, 8)
        _, state["k_stretch"] = psim.SliderFloat("stretch stiffness k",
                                                  float(state["k_stretch"]), 0.1, 1.0)
        _, state["k_bend"] = psim.SliderFloat("bend stiffness k",
                                               float(state["k_bend"]), 0.0, 1.0)
        _, state["damp"] = psim.SliderFloat("internal damping",
                                              float(state["damp"]), 0.0, 0.5)
        _, state["throw_speed"] = psim.SliderFloat("throw speed (m/s)",
                                                     float(state["throw_speed"]), 0.5, 12.0)
        _, state["ball_mass"] = psim.SliderFloat("ball mass (kg)",
                                                   float(state["ball_mass"]), 0.05, 5.0)
        sphere.mass = state["ball_mass"]

        if state["play"]:
            n_sub = max(int(state["substeps"]), 1)
            h_sub = h / n_sub
            iters = max(int(state["iters"]), 1)
            kp_s = 1.0 - (1.0 - state["k_stretch"]) ** (1.0 / iters)
            kp_b = 1.0 - (1.0 - state["k_bend"]) ** (1.0 / iters) if state["k_bend"] > 0 else 0.0
            for _ in range(n_sub):
                _step(state, sphere, cloth, w, cloth_mass, h_sub, iters, kp_s, kp_b)
            state["n"] += 1

        sm.update_vertex_positions(state["x"])
        pin_pc.update_point_positions(state["x"][cloth.pin_idx])
        sphere_mesh.update_vertex_positions(sphere_v_local + sphere.center)
        psim.Separator()
        psim.Text(f"step {state['n']}   ball pos = {sphere.center.round(2)}   "
                  f"v = {sphere.velocity.round(2)}")
        psim.Text("throw, then watch: too few iters/sub-steps and the ball passes through.")

    viz.show(callback)
    return state


def _step(state, sphere, cloth, w, cloth_mass, h_sub, iters, kp_s, kp_b):
    x, v = state["x"], state["v"]
    v_new = v.copy()
    v_new[:, 2] -= h_sub * state["g"]
    v_new[w == 0] = 0.0
    p = x + h_sub * v_new
    if state["ball_locked"]:
        # ball hasn't been released — keep it static, no integration
        P = sphere.center.copy()
        old_center = sphere.center.copy()
        sphere.impulse = np.zeros(3)
        for _ in range(iters):
            project_distance(p, cloth.edges, cloth.rest_lengths, w, k_prime=kp_s)
            if kp_b > 0:
                project_bend(p, cloth.bend_pairs, cloth.rest_dihedrals, w, k_prime=kp_b)
            project_dynamic_sphere(p, x, w, cloth_mass, sphere, h_sub,
                                   thickness=state["thickness"])
        sphere.impulse = np.zeros(3)                     # discard impulse on locked ball
    else:
        P = sphere.predict(h_sub, gravity=(0.0, 0.0, -state["g"]))
        old_center = sphere.center.copy()
        sphere.center = P                                # use predicted center for collisions
        for _ in range(iters):
            project_distance(p, cloth.edges, cloth.rest_lengths, w, k_prime=kp_s)
            if kp_b > 0:
                project_bend(p, cloth.bend_pairs, cloth.rest_dihedrals, w, k_prime=kp_b)
            project_dynamic_sphere(p, x, w, cloth_mass, sphere, h_sub,
                                   thickness=state["thickness"])
        sphere.center = old_center                       # restore for commit()
        sphere.commit(P, h_sub)
    v_post = (p - x) / h_sub
    m = np.where(w > 0, 1.0, np.inf)
    damp_velocities(p, v_post, m, state["damp"])
    state["v"] = v_post
    state["x"] = p


def _icosphere(center, radius, subdiv=2):
    """A small icosphere builder so we don't need a mesh file for the ball."""
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
