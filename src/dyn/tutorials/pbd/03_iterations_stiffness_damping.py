"""Lesson 03 — Iterations, stiffness k', and rigid-mode-preserving damping.

    Run :  python src/dyn/tutorials/pbd/03_iterations_stiffness_damping.py
    Read:  src/dyn/tutorials/pbd/03_iterations_stiffness_damping.md

Two ropes side by side share gravity, mass, and rest length but diverge on
iteration count and damping mode.  Slide the iters knob to watch one rope go
from mushy to inextensible while the other holds steady; flip damping mode to
see Müller's rigid-mode trick keep the global swing alive when the naive
damper would kill it.
"""

from __future__ import annotations

import os

import numpy as np

from dyn import viz
from dyn.pbd import predict, project_distance, damp_velocities


def build_rope(n=24, length=1.5, x_offset=0.0):
    x = np.zeros((n, 3))
    x[:, 0] = np.linspace(0.0, length, n) + x_offset
    x[:, 2] = 0.5
    edges = np.stack([np.arange(n - 1), np.arange(1, n)], axis=1)
    rest = np.full(n - 1, length / (n - 1))
    w = np.ones(n)
    w[0] = 0.0                                          # pin leftmost
    return x, edges, rest, w


def main(show=True):
    h = 0.02
    nseg = 24
    L = 1.5
    rope_a = build_rope(nseg, L, x_offset=0.0)
    rope_b = build_rope(nseg, L, x_offset=2.5)

    # numerical demonstration: with k' correction the stretch ratio at convergence
    # is independent of iteration count, while naive k * delta is not
    print("=== Lesson 3: iterations, stiffness, damping ===")
    for iters in (1, 5, 30):
        x, edges, rest, w = build_rope(nseg, L)
        v = np.zeros_like(x)
        for _ in range(80):
            p = predict(x, v, w, h, gravity=(0, 0, -9.8))
            for _ in range(iters):
                project_distance(p, edges, rest, w, k_prime=1.0)
            v = (p - x) / h
            x = p
        seg = np.linalg.norm(x[1:] - x[:-1], axis=1)
        print(f"  iters={iters:>2}: stretch ratio at t=1.6s = {seg.sum() / rest.sum():.3f}")

    if not show:
        return

    import polyscope as ps
    import polyscope.imgui as psim

    viz.init_once()
    xa, ea, ra, wa = rope_a
    xb, eb, rb, wb = rope_b
    cn_a = ps.register_curve_network("rope A (low iters)", xa, ea,
                                     color=(0.85, 0.55, 0.20), radius=0.012)
    cn_b = ps.register_curve_network("rope B (high iters)", xb, eb,
                                     color=(0.20, 0.55, 0.85), radius=0.012)

    state = {"play": True, "iters_a": 1, "iters_b": 30, "k": 1.0, "use_kp": True,
             "damp": 0.04, "rigid_damp": True, "g": 9.8, "n": 0,
             "x_a": xa.copy(), "v_a": np.zeros_like(xa),
             "x_b": xb.copy(), "v_b": np.zeros_like(xb)}

    def reset():
        state["x_a"] = xa.copy(); state["v_a"][:] = 0.0
        state["x_b"] = xb.copy(); state["v_b"][:] = 0.0
        state["n"] = 0

    def step_one(x, v, w, edges, rest, iters):
        p = predict(x, v, w, h, gravity=(0, 0, -state["g"]))
        kp = (1.0 - (1.0 - state["k"]) ** (1.0 / max(iters, 1))) if state["use_kp"] else state["k"]
        for _ in range(iters):
            project_distance(p, edges, rest, w, k_prime=kp)
        v_new = (p - x) / h
        # damping (rigid-mode-preserving vs naive)
        if state["rigid_damp"]:
            m = np.where(w > 0, 1.0 / np.maximum(w, 1e-12), np.inf)
            damp_velocities(p, v_new, m, state["damp"])
        else:
            v_new *= (1.0 - state["damp"])
        return p, v_new

    def callback():
        psim.Text("Lesson 3 — iterations vs stiffness vs damping mode")
        psim.Separator()
        _, state["play"] = psim.Checkbox("Play", state["play"])
        if psim.Button("reset"):
            reset()
        _, state["iters_a"] = psim.SliderInt("iters (rope A, orange)", state["iters_a"], 1, 40)
        _, state["iters_b"] = psim.SliderInt("iters (rope B, blue)", state["iters_b"], 1, 40)
        _, state["k"] = psim.SliderFloat("stiffness k", float(state["k"]), 0.05, 1.0)
        _, state["use_kp"] = psim.Checkbox("use k' (iteration-independent)", state["use_kp"])
        _, state["damp"] = psim.SliderFloat("damping coefficient", float(state["damp"]), 0.0, 1.0)
        _, state["rigid_damp"] = psim.Checkbox("rigid-mode-preserving damping (off = naive v ← (1−c)v)",
                                               state["rigid_damp"])
        _, state["g"] = psim.SliderFloat("|gravity|", float(state["g"]), 0.0, 20.0)

        if state["play"]:
            state["x_a"], state["v_a"] = step_one(state["x_a"], state["v_a"], wa, ea, ra,
                                                   state["iters_a"])
            state["x_b"], state["v_b"] = step_one(state["x_b"], state["v_b"], wb, eb, rb,
                                                   state["iters_b"])
            state["n"] += 1

        cn_a.update_node_positions(state["x_a"])
        cn_b.update_node_positions(state["x_b"])
        seg_a = np.linalg.norm(state["x_a"][1:] - state["x_a"][:-1], axis=1).sum()
        seg_b = np.linalg.norm(state["x_b"][1:] - state["x_b"][:-1], axis=1).sum()
        psim.Separator()
        psim.Text(f"step {state['n']}   stretch A: {seg_a / ra.sum():.3f}   "
                  f"stretch B: {seg_b / rb.sum():.3f}  (1.0 = inextensible)")
        psim.Text("Toggle k' to feel the timestep-trap; toggle rigid-mode damping at k_d=1.")

    viz.show(callback)


if __name__ == "__main__":
    main(show=not os.environ.get("DYN_NO_SHOW"))
