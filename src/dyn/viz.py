"""Thin Polyscope helpers shared by the lessons.

The goal is to keep visualization boilerplate out of the lesson scripts so each
one reads as the *idea* plus a small callback.  Nothing here is conceptually
deep — it just wraps registration, node tagging, and the show loop.

Headless verification: set ``DYN_SMOKE=1`` to use Polyscope's mock GL backend
and tick a few frames instead of opening a blocking window.  This exercises the
full registration + callback path in CI without a display.
"""

from __future__ import annotations

import os

import numpy as np
import polyscope as ps

_INITED = False


def smoke() -> bool:
    return bool(os.environ.get("DYN_SMOKE"))


def init_once() -> None:
    """Initialize Polyscope exactly once (mock backend under DYN_SMOKE)."""
    global _INITED
    if _INITED:
        return
    ps.init(backend="openGL_mock") if smoke() else ps.init()
    # The beam bends in z; make z up and keep a soft ground shadow.
    if hasattr(ps, "set_up_dir"):
        ps.set_up_dir("z_up")
    if hasattr(ps, "set_ground_plane_mode"):
        ps.set_ground_plane_mode("shadow_only")
    _INITED = True


def register_beam(beam, positions=None, name="beam", color=(0.85, 0.80, 0.70), **kw):
    """Register the beam as a Polyscope volume mesh (renders its boundary)."""
    init_once()
    pos = beam.verts if positions is None else positions
    return ps.register_volume_mesh(name, pos, tets=beam.tets, color=color, **kw)


def add_node_groups(vol, beam, enabled=True):
    """Categorical scalar tagging vertices: 0 interior, 1 surface, 2 clamped."""
    label = np.zeros(beam.n)
    label[beam.surface_idx] = 1
    label[beam.clamp_idx] = 2
    vol.add_scalar_quantity(
        "node type (0=interior, 1=surface, 2=clamp)",
        label, defined_on="vertices", datatype="categorical", enabled=enabled,
    )
    return label


def add_vectors(struct, name, vecs, defined_on="vertices", enabled=True, **kw):
    """Add a per-vertex vector quantity (PointCloud has no ``defined_on``)."""
    if defined_on is None:
        struct.add_vector_quantity(name, vecs, enabled=enabled, **kw)
    else:
        struct.add_vector_quantity(name, vecs, defined_on=defined_on, enabled=enabled, **kw)


def show(callback=None, frames: int = 3) -> None:
    """Open the interactive window — or tick a few mock frames under DYN_SMOKE."""
    init_once()
    if callback is not None:
        ps.set_user_callback(callback)
    if smoke():
        for _ in range(frames):
            ps.frame_tick()
    else:
        ps.show()
