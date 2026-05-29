"""Lesson 02 — The distance constraint, projection, and a hanging rope.

    Run :  python src/dyn/tutorials/pbd/02_distance_constraint.py
    Read:  src/dyn/tutorials/pbd/02_distance_constraint.md

A chain of N particles, top pinned (w=0), connected by distance constraints.
The same predict-then-project loop as Lesson 1, but now `project_distance`
keeps each segment at its rest length.  Crank the iterations slider: Jacobi
corrections take many sweeps to propagate down a long chain, so the rope is
stretchy at low iters and rigid at high.
"""

from __future__ import annotations

import os

import numpy as np

from dyn import viz
from dyn.pbd import predict, project_distance


TRAIL = 240


def build_rope(n=20, length=1.5):
    """N particles in a horizontal line, top index 0 will be pinned."""
    x = np.zeros((n, 3))
    x[:, 0] = np.linspace(0.0, length, n)               # along +x
    x[:, 2] = 0.5                                       # at z = 0.5
    edges = np.stack([np.arange(n - 1), np.arange(1, n)], axis=1)
    rest = np.full(n - 1, length / (n - 1))
    w = np.ones(n)
    w[0] = 0.0                                          # pin the leftmost
    return x, edges, rest, w


def main(show=True):
    x0, edges, rest, w = build_rope(n=20, length=1.5)
    print("=== Lesson 2: distance constraint, hanging rope ===")
    print(f"particles N = {len(x0)}, edges E = {len(edges)}, pinned = {(w == 0).sum()}")
    # one Jacobi sweep on a stretched edge fully closes the gap (single edge)
    p = np.array([[0.0, 0, 0], [3.0, 0, 0]])
    edges1 = np.array([[0, 1]])
    rest1 = np.array([1.0])
    w1 = np.ones(2)
    project_distance(p, edges1, rest1, w1, k_prime=1.0)
    print(f"sanity: stretched 3.0 → rest 1.0 in one sweep gives |p0-p1| = "
          f"{np.linalg.norm(p[0] - p[1]):.4f}")
    assert np.isclose(np.linalg.norm(p[0] - p[1]), 1.0)

    if not show:
        return x0, edges, rest, w

    import polyscope as ps
    import polyscope.imgui as psim

    viz.init_once()
    cn = ps.register_curve_network("rope", x0, edges, color=(0.85, 0.55, 0.20),
                                   radius=0.012)
    pin_pc = ps.register_point_cloud("pin", x0[w == 0], color=(0.90, 0.20, 0.20),
                                     radius=0.018)
    state = {"play": True, "project_on": True, "pin": True, "iters": 8,
             "g": 9.8, "x": x0.copy(), "v": np.zeros_like(x0), "n": 0,
             "trail": []}
    h = 0.02
    tip = len(x0) - 1                                   # the tip vertex (free end)

    def reset():
        state["x"] = x0.copy()
        state["v"][:] = 0.0
        state["n"] = 0
        state["trail"] = []

    def callback():
        psim.Text("Lesson 2 — distance constraints make the chain a rope")
        psim.Separator()
        _, state["play"] = psim.Checkbox("Play", state["play"])
        psim.SameLine()
        _, state["project_on"] = psim.Checkbox("projection", state["project_on"])
        psim.SameLine()
        ch_pin, state["pin"] = psim.Checkbox("pin top", state["pin"])
        if psim.Button("reset"):
            reset()
        _, state["iters"] = psim.SliderInt("iterations", int(state["iters"]), 1, 40)
        _, state["g"] = psim.SliderFloat("|gravity|", float(state["g"]), 0.0, 20.0)
        if ch_pin:
            w[0] = 0.0 if state["pin"] else 1.0

        if state["play"]:
            p = predict(state["x"], state["v"], w, h, gravity=(0, 0, -state["g"]))
            if state["project_on"]:
                for _ in range(state["iters"]):
                    project_distance(p, edges, rest, w, k_prime=1.0)
            state["v"] = (p - state["x"]) / h
            state["x"] = p
            state["n"] += 1
            state["trail"].append(state["x"][tip].copy())
            state["trail"] = state["trail"][-TRAIL:]

        cn.update_node_positions(state["x"])
        if w[0] == 0.0:
            pin_pc.update_point_positions(state["x"][w == 0])
        if len(state["trail"]) > 2:
            ps.register_curve_network("tip trail", np.array(state["trail"]),
                                      "line", color=(0.95, 0.55, 0.1), radius=0.003)
        # measure rope length so we can see stretching at low iter counts
        seg = np.linalg.norm(state["x"][1:] - state["x"][:-1], axis=1)
        psim.Separator()
        psim.Text(f"step {state['n']}   rope length / rest = "
                  f"{seg.sum():.3f} / {rest.sum():.3f}   "
                  f"(ratio {seg.sum() / rest.sum():.3f})")
        psim.Text("Low iterations → stretchy.  High iterations → near-inextensible.")

    viz.show(callback)
    return state


if __name__ == "__main__":
    main(show=not os.environ.get("DYN_NO_SHOW"))
