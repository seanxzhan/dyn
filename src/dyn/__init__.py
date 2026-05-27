"""dyn — a from-scratch, visual tour of rig-space physics.

The package is deliberately small and readable: each module maps to one idea
in the Hahn et al. papers (see ``docs/``).  The runnable lessons live in
``src/dyn/tutorials`` and open an interactive Polyscope window.

    mesh : the discretized continuum (a tetrahedral beam + lumped mass)
    rig  : maps from a few parameters p to all vertices, s(p), with Jacobians
    viz  : thin Polyscope helpers shared by the lessons

Submodules are imported lazily (``from dyn.mesh import make_beam``) so that the
pure-NumPy parts never drag in Polyscope.
"""

__version__ = "0.1.0"
