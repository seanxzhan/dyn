"""Lesson 04 — Elastic energy (FEM / StVK): the potential W and its forces.

    Run :  python src/dyn/tutorials/04_elastic_energy.py
    Read:  src/dyn/tutorials/04_elastic_energy.md

Pose the beam with the rig sliders and watch the StVK energy density light up
where the beam is strained (highest near the clamp).  Toggle the elastic forces
(-dW/dx) to see them pull the mesh back toward rest.  Energy is rotation
invariant: a rigid rotation stores none.
"""

import os

import numpy as np

from dyn import viz
from dyn.energy import StVK
from dyn.mesh import make_beam
from dyn.rig import LinearRig


def _rigid_rotation(verts, deg=35.0):
    a = np.radians(deg)
    R = np.array([[np.cos(a), 0, np.sin(a)], [0, 1, 0], [-np.sin(a), 0, np.cos(a)]])
    return verts @ R.T


def main(show=True):
    beam = make_beam()
    mat = StVK(beam, young=200.0, poisson=0.40)
    rig = LinearRig(beam, modes=("bend_z", "sway_y"))

    print("=== Lesson 4: elastic energy W(x) ===")
    print(f"Lamé params: mu = {mat.mu:.2f},  lambda = {mat.lam:.2f}  "
          f"(from E={mat.young}, nu={mat.poisson})")
    print(f"rest energy  W(rest) = {mat.energy(beam.verts):.3e}   (F=I -> E=0)")
    print(f"rigid-rotation energy = {mat.energy(_rigid_rotation(beam.verts)):.3e}   "
          f"(StVK is rotation invariant)")
    # directional check that the analytic gradient is really dW/dx
    x = beam.verts + 0.02 * np.random.default_rng(1).standard_normal(beam.verts.shape)
    d = np.random.default_rng(2).standard_normal(beam.verts.shape)
    fd = (mat.energy(x + 1e-6 * d) - mat.energy(x - 1e-6 * d)) / 2e-6
    an = float(np.sum(mat.gradient(x) * d))
    print(f"gradient directional check: analytic {an:.5f} vs finite-diff {fd:.5f}")
    assert abs(an - fd) < 1e-3 * (abs(fd) + 1)

    if not show:
        return mat

    import polyscope.imgui as psim

    v = viz.register_beam(beam)
    state = {"p": np.array([-0.6, 0.3]), "show_forces": True, "force_scale": 1.0}

    def callback():
        psim.Text("Lesson 4 — elastic potential W(x) = Σ_e Ψ(F_e)·V0_e")
        psim.Separator()
        changed = False
        for i in range(rig.dim):
            ch, state["p"][i] = psim.SliderFloat(f"p[{i}] ({rig.names[i]})",
                                                 float(state["p"][i]), -1.5, 1.5)
            changed |= ch
        _, state["show_forces"] = psim.Checkbox("show elastic forces -dW/dx", state["show_forces"])

        x = rig.s(state["p"])
        v.update_vertex_positions(x)
        # per-tet energy density (heatmap) and total energy
        Fd = mat.deformation_gradient(x)
        psi = mat.energy_density(Fd)
        v.add_scalar_quantity("energy density Ψ (per tet)", psi, defined_on="cells",
                              datatype="standard", enabled=True, cmap="reds")
        if state["show_forces"]:
            viz.add_vectors(v, "elastic force  -dW/dx", mat.forces(x),
                            color=(0.1, 0.3, 0.9), length=0.4)
        else:
            try:
                v.remove_quantity("elastic force  -dW/dx")
            except Exception:
                pass
        psim.Separator()
        psim.Text(f"total elastic energy  W = {mat.energy(x):.3f}")
        psim.Text("Strain (and energy) concentrates where the beam bends most — near the clamp.")

    viz.show(callback)
    return mat


if __name__ == "__main__":
    main(show=not os.environ.get("DYN_NO_SHOW"))
