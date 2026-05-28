"""Lesson 09 — Rig-space physics on a real FBX character.

    Run :  python src/dyn/tutorials/09_fbx_character.py
    Read:  src/dyn/tutorials/09_fbx_character.md

Everything from Lessons 1-8, now driven by an artist's skeleton instead of a toy
rig.  We load data/2072.fbx (a fish) with our pure-Python FBX reader, voxelize it
into a tet proxy (no external mesher), skin the proxy with the same bones, then
run the 2013 (linearized-rig) solver.  By default the fish is held by its
mid-body and the rear half + tail droop and jiggle under gravity — the cantilever
beam of Lesson 6, now as a character.  Optionally drive the body to make it swim.

The rig is a black box (linear-blend skinning); its Jacobian comes from finite
differences, exactly as in the papers.

Scale note: this fish model is ~1 unit long and very light (proxy mass ~0.02),
so a large gravity constant (default 150) is needed to load it visibly — this is
equivalent to a softer material.  Treat "gravity" here as a load knob.
"""

import os

import numpy as np

from dyn import viz
from dyn.energy import StVK
from dyn.fbx import SkeletonRig, load_character
from dyn.solver import ImplicitEuler
from dyn.voxel import transfer_weights, voxelize

# Spine-driven swim / droop: c_body_1 is the PRIMARY (driven) bone; the rest of
# the spine + fins are SECONDARY (simulated).  Only bones that meaningfully
# deform the coarse voxel proxy can be free DOFs — a thin fin barely strains any
# tets, giving a near-zero-stiffness mode that blows up — so we filter candidates
# by how much proxy skin weight they carry (MIN_PROXY_WEIGHT).
PRIMARY = "c_body_1"
SECONDARY = ["c_body_2", "c_body_3", "c_caudal_fin", "c_dorsal_fin"]
TAIL_BONES = ["c_caudal_fin", "c_body_3", "c_body_2"]   # for the trajectory trail
MIN_PROXY_WEIGHT = 2.0
# Each free joint is a 1-DOF HINGE about world x (the horizontal axis ⟂ to the
# fish's length z and to gravity y), so the rear droops cleanly in the vertical
# plane instead of buckling sideways.  The driven bone yaws about y (the swim).
SWIM_AXIS = np.array([[0.0, 1.0, 0.0]])
BEND_AXIS = np.array([[1.0, 0.0, 0.0]])
YOUNG, GRAVITY, DAMPING = 15.0, 90.0, 0.02
TRAIL = 220


def _pick(names, wanted):
    idx = {n: i for i, n in enumerate(names)}
    return [idx[w] for w in wanted if w in idx]


def main(show=True, path="data/2072.fbx", res=22):
    here = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
    path = os.path.join(here, path) if not os.path.isabs(path) else path

    ch = load_character(path)
    print(f"=== Lesson 9: rig-space physics on {os.path.basename(path)} ===")
    print(f"surface: {ch.V.shape[0]} verts, {ch.faces.shape[0]} tris, {len(ch.bone_names)} bones")
    proxy = voxelize(ch.V, ch.faces, res=res)
    Wp = transfer_weights(ch.V, ch.weights, proxy.verts, k=3)
    print(f"voxel proxy: {proxy.verts.shape[0]} verts, {len(proxy.tets)} tets")

    driven = _pick(ch.bone_names, [PRIMARY])
    free_b = [b for b in _pick(ch.bone_names, SECONDARY) if Wp[:, b].sum() >= MIN_PROXY_WEIGHT]
    control = driven + free_b
    axes = [SWIM_AXIS] * len(driven) + [BEND_AXIS] * len(free_b)   # 1 DOF per bone
    free_idx = list(range(len(driven), len(control)))             # one param per bone
    print(f"PRIMARY (driven): {[ch.bone_names[b] for b in driven]}")
    print(f"SECONDARY (free): {[ch.bone_names[b] for b in free_b]}  -> {len(free_idx)} simulated DOFs")

    rig_proxy = SkeletonRig(proxy.verts, Wp, ch, control_bones=control, axes=axes)   # physics rig
    rig_surf = SkeletonRig(ch.V, ch.weights, ch, control_bones=control, axes=axes)   # rendering rig
    assert np.allclose(rig_proxy.s(np.zeros(rig_proxy.dim)), proxy.verts, atol=1e-9)

    # tail vertex (most-influenced by the tail-most free bone) for the trail
    tail_bone = next((b for b in _pick(ch.bone_names, TAIL_BONES) if b in control), control[-1])
    tip = int(np.argmax(ch.weights[:, tail_bone]))
    body_len = float(np.ptp(ch.V[:, 2]))

    mat = StVK(proxy, young=YOUNG, poisson=0.40)

    def new_sim(grav, damp):
        return ImplicitEuler(rig_proxy, mat, proxy.mass, h=0.02, gravity=(0, -grav, 0),
                             damping=damp, free_idx=free_idx, relinearize=True)

    sim = new_sim(GRAVITY, DAMPING)
    for _ in range(200):                                   # release-under-gravity check
        sim.step(newton_iter=5)
    disp = np.linalg.norm(rig_surf.s(sim.p)[tip] - ch.V[tip]) / body_len
    print(f"held mid-body, rear droops: tail settled at {100 * disp:.1f}% of body length")

    if not show:
        return sim

    import polyscope as ps
    import polyscope.imgui as psim

    viz.init_once()
    if hasattr(ps, "set_up_dir"):
        ps.set_up_dir("y_up")                              # fish lies flat; gravity is -y
    sim = new_sim(GRAVITY, DAMPING)
    surf = ps.register_surface_mesh("character (physics)", ch.V, ch.faces,
                                    color=(0.55, 0.65, 0.85), smooth_shade=True)
    ghost = ps.register_surface_mesh("rig pose (no physics)", ch.V, ch.faces,
                                     color=(0.6, 0.6, 0.6), transparency=0.3, enabled=True)
    proxy_vm = ps.register_volume_mesh("physics proxy (tets)", proxy.verts, tets=proxy.tets,
                                       color=(0.95, 0.6, 0.2), transparency=0.4, enabled=False)
    st = {"play": True, "physics": True, "A": 0.0, "w": 5.0, "damp": DAMPING,
          "grav": GRAVITY, "young": YOUNG, "show_proxy": False, "show_ghost": True,
          "n": 0, "trail": []}

    def pluck():
        v = (sim.x - sim.x_prev) / sim.h
        kick = np.zeros_like(sim.x)
        rear = proxy.verts[:, 2] < np.median(proxy.verts[:, 2])    # tail half
        kick[rear, 1] = 3.0                                        # flick up
        sim.pluck(v + kick)

    def callback():
        psim.Text("Lesson 9 — secondary motion on a real skeleton (linear-blend-skinning rig)")
        psim.Text("Held by the mid-body; the rear + tail droop and jiggle (a cantilever fish).")
        psim.Separator()
        _, st["play"] = psim.Checkbox("Play", st["play"]); psim.SameLine()
        _, st["physics"] = psim.Checkbox("secondary physics on", st["physics"])
        _, st["show_proxy"] = psim.Checkbox("show physics proxy (tets)", st["show_proxy"])
        psim.SameLine()
        _, st["show_ghost"] = psim.Checkbox("show rig-pose ghost", st["show_ghost"])
        if psim.Button("reset"):
            sim.__dict__.update(new_sim(st["grav"], st["damp"]).__dict__)
            st["n"], st["trail"] = 0, []
        psim.SameLine()
        if psim.Button("pluck tail"):
            pluck()
        ch_g, st["grav"] = psim.SliderFloat("gravity load", float(st["grav"]), 0.0, 400.0)
        ch_y, st["young"] = psim.SliderFloat("stiffness E", float(st["young"]), 5.0, 300.0)
        _, st["damp"] = psim.SliderFloat("damping", float(st["damp"]), 0.0, 0.2)
        _, st["A"] = psim.SliderFloat("swim: primary amplitude A", float(st["A"]), 0.0, 1.0)
        _, st["w"] = psim.SliderFloat("swim: frequency ω", float(st["w"]), 0.5, 18.0)

        if ch_y:
            sim.elastic = StVK(proxy, young=st["young"], poisson=0.40)
        sim.fg = sim.M * np.array([0.0, -st["grav"], 0.0])
        sim.damping = st["damp"]

        yaw = st["A"] * np.sin(st["w"] * st["n"] * sim.h)      # primary swim (1 DOF)
        if st["play"]:
            if st["physics"]:
                sim.step(driven_values=[yaw], newton_iter=5)
            else:
                sim.p[:] = 0.0
                sim.p[sim.driven_idx] = yaw                # primary only, secondary frozen
                sim.x = rig_proxy.s(sim.p)
                sim.x_prev = sim.x.copy()
            st["n"] += 1

        # render the smooth surface skinned by the *simulated* joint angles
        surf.update_vertex_positions(rig_surf.s(sim.p))
        prim = sim.p.copy(); prim[free_idx] = 0.0
        ghost.update_vertex_positions(rig_surf.s(prim))    # rig pose with NO secondary
        ghost.set_enabled(st["show_ghost"])
        proxy_vm.set_enabled(st["show_proxy"])
        if st["show_proxy"]:
            proxy_vm.update_vertex_positions(sim.x)

        st["trail"].append(rig_surf.s(sim.p)[tip].copy())
        st["trail"] = st["trail"][-TRAIL:]
        if len(st["trail"]) > 2:
            ps.register_curve_network("tail trajectory", np.array(st["trail"]), "line",
                                      color=(0.95, 0.55, 0.1), radius=0.0025)
        psim.Separator()
        tail_disp = np.linalg.norm(rig_surf.s(sim.p)[tip] - ch.V[tip]) / body_len
        psim.Text(f"step {st['n']}   tail displaced {100 * tail_disp:.1f}% of body length   "
                  f"peak free-rot {np.abs(sim.p[free_idx]).max():+.2f} rad")
        psim.Text("Grey ghost = rig with secondary=0 (no physics). The gap + tail trail is the physics.")

    viz.show(callback)
    return sim


if __name__ == "__main__":
    main(show=not os.environ.get("DYN_NO_SHOW"))
