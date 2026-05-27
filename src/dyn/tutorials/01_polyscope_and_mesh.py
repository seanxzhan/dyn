"""Lesson 01 — The discretized beam and the DOF split.

    Run :  python src/dyn/tutorials/01_polyscope_and_mesh.py
    Read:  src/dyn/tutorials/01_polyscope_and_mesh.md

What you should see: a tetrahedral beam, colored by "node type".  The clamped
end (red points) is the Dirichlet boundary.  This sets up the vocabulary —
vertices ``x`` split into surface nodes ``s`` (later driven by the rig) and
interior nodes ``n`` (free FEM nodes) — used by every later lesson.
"""

import os

import numpy as np

from dyn import viz
from dyn.mesh import make_beam, tet_volumes


def main(show=True):
    beam = make_beam(res=(10, 3, 3), lengths=(4.0, 1.0, 1.0), density=1.0)
    vol = tet_volumes(beam.verts, beam.tets)

    print("=== Lesson 1: the discretized beam ===")
    print(f"vertices      N = {beam.n}")
    print(f"tets          T = {len(beam.tets)}  (all positively oriented: {bool(np.all(vol > 0))})")
    print(f"surface  |s| = {len(beam.surface_idx):4d}   <- the rig will drive these")
    print(f"interior |n| = {len(beam.interior_idx):4d}   <- free FEM nodes (later lessons)")
    print(f"clamped      = {len(beam.clamp_idx):4d}   <- fixed Dirichlet boundary (x=0 face)")
    print(f"total mass     = {beam.mass.sum():.4f}  (= density * volume = {vol.sum():.4f})")
    # The lumped mass must conserve total mass; tets must not be inverted.
    assert np.isclose(beam.mass.sum(), vol.sum())
    assert np.all(vol > 0)

    if not show:
        return beam

    import polyscope as ps
    import polyscope.imgui as psim

    v = viz.register_beam(beam)
    viz.add_node_groups(v, beam)                       # categorical color: interior/surface/clamp
    ps.register_point_cloud("clamped nodes", beam.verts[beam.clamp_idx],
                            color=(0.90, 0.20, 0.20), radius=0.012)

    def callback():
        psim.Text("Lesson 1 — the discretized continuum")
        psim.Separator()
        psim.Text("The beam is a tetrahedral mesh (a lumped-mass FEM body).")
        psim.Text("Its vertices split into  x = {n} interior  U  {s} surface.")
        psim.Text("Enable the 'node type' quantity to see the split.")
        psim.Text("Red points = clamped boundary; they never move.")

    viz.show(callback)
    return beam


if __name__ == "__main__":
    main(show=not os.environ.get("DYN_NO_SHOW"))
