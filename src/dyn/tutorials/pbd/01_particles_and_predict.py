"""Lesson 01 — Particles and the predict step (the loop without projection).

    Run :  python src/dyn/tutorials/pbd/01_particles_and_predict.py
    Read:  src/dyn/tutorials/pbd/01_particles_and_predict.md

A 6x6 grid of point masses falls under gravity.  With no constraints to
project, the loop reduces to plain explicit Euler — particles are independent
parabolas.  This sets the data layout (x, v, w) and the predict-then-commit
shell that Lessons 2-6 keep, only swapping in different projection steps.
"""

from __future__ import annotations

import os

import numpy as np

from dyn import viz
from dyn.pbd import predict


TRAIL = 240


def main(show=True):
    # 6x6 grid of particles in the xz-plane (y = 0).  No constraints yet.
    nx, nz = 6, 6
    xs = np.linspace(-0.5, 0.5, nx)
    zs = np.linspace(0.0, 1.0, nz)
    X, Z = np.meshgrid(xs, zs, indexing="ij")
    x = np.stack([X.ravel(), np.zeros(nx * nz), Z.ravel()], axis=1)
    v = np.zeros_like(x)
    w = np.ones(x.shape[0])                       # all unit mass, all free

    print("=== Lesson 1: particles + predict (no constraints) ===")
    print(f"particles: {x.shape[0]}  (= {nx} * {nz})")
    print(f"shapes   : x {x.shape}, v {v.shape}, w {w.shape}")
    # one step of pure explicit Euler under gravity, by hand
    h = 0.02
    p = predict(x.copy(), v, w, h, gravity=(0, 0, -9.8))
    v_after = (p - x) / h
    print(f"after one step (h={h}): max |v_z| = {np.max(np.abs(v_after[:, 2])):.4f}  "
          f"(= h * |g| = {h * 9.8:.4f})")
    assert np.allclose(v_after[:, 2], -h * 9.8)            # gravity-only freefall: v = h·g

    if not show:
        return x, v, w

    import polyscope as ps
    import polyscope.imgui as psim

    viz.init_once()
    pc = ps.register_point_cloud("particles", x, radius=0.012, color=(0.85, 0.55, 0.20))
    state = {"play": True, "predict_on": True, "g": 9.8, "n": 0,
             "x": x.copy(), "v": v.copy(), "w": w.copy(), "trail": []}
    track_idx = (nx // 2) * nz + (nz - 1)                 # one corner particle

    def reset():
        state["x"] = x.copy()
        state["v"] = v.copy()
        state["n"] = 0
        state["trail"] = []

    def callback():
        psim.Text("Lesson 1 — predict-then-commit, no projection step yet")
        psim.Separator()
        _, state["play"] = psim.Checkbox("Play", state["play"])
        psim.SameLine()
        _, state["predict_on"] = psim.Checkbox("predict", state["predict_on"])
        if psim.Button("reset"):
            reset()
        _, state["g"] = psim.SliderFloat("|gravity|", float(state["g"]), 0.0, 20.0)

        if state["play"]:
            if state["predict_on"]:
                p = predict(state["x"], state["v"], state["w"], h,
                            gravity=(0, 0, -state["g"]))
                state["v"] = (p - state["x"]) / h          # line 13
                state["x"] = p                              # line 14
            state["n"] += 1
            state["trail"].append(state["x"][track_idx].copy())
            state["trail"] = state["trail"][-TRAIL:]

        pc.update_point_positions(state["x"])
        if len(state["trail"]) > 2:
            ps.register_curve_network("trail (one particle)", np.array(state["trail"]),
                                      "line", color=(0.95, 0.55, 0.1), radius=0.003)
        psim.Separator()
        psim.Text(f"step {state['n']}   z range "
                  f"[{state['x'][:, 2].min():+.2f}, {state['x'][:, 2].max():+.2f}]")
        psim.Text("Without constraints these particles are independent — no cloth.")
        psim.Text("Lesson 2 adds distance constraints; the same loop becomes a rope.")

    viz.show(callback)
    return state


if __name__ == "__main__":
    main(show=not os.environ.get("DYN_NO_SHOW"))
