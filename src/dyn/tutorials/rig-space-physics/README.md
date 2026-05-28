# Rig-Space Physics — a step-by-step, visual course

A from-scratch tour of the two Hahn et al. papers (`docs/`), built on a
procedural beam so every formula is checkable in closed form, and visualized
interactively with [Polyscope](https://polyscope.run). Each lesson is a runnable
script paired with a math note; shared code lives one level up in
[`../mesh.py`](../mesh.py), [`../rig.py`](../rig.py), [`../viz.py`](../viz.py).

## Setup & running

```bash
source activate.sh          # conda activate dyn   (numpy + polyscope already installed)
python src/dyn/tutorials/01_polyscope_and_mesh.py
```

Each script opens a Polyscope window with sliders/checkboxes; read the matching
`.md` alongside it. Close the window to exit.

- `DYN_SMOKE=1` — run headless (mock GL backend, a few frame ticks) to verify a
  lesson without opening a window. Used for the checks below.
- `DYN_NO_SHOW=1` — run only the NumPy part (math + asserts), skip Polyscope.

## The arc

Each lesson adds exactly one idea. The thread mirrors the literature: start in a
**linear** rig (constant Jacobian — the regime the 2013 method engineers), then a
**nonlinear** rig where the rig *curvature* appears (the 2012 bottleneck).

| # | Lesson | One idea | Status |
|---|--------|----------|--------|
| 01 | [Discretization & DOF split](01_polyscope_and_mesh.md) | continuum → tet mesh, lumped mass, $x=\{s\}\cup\{n\}$ | ✅ built |
| 02 | [The rig & reduced coords](02_the_rig_and_reduced_coords.md) | $s(p)=s_0+Bp$; rig space is a subspace | ✅ built |
| 03 | [The rig Jacobian](03_the_rig_jacobian.md) | $J=\partial s/\partial p$; forces project as $f_p=J^{\mathsf T}f_s$ | ✅ built |
| 04 | [Elastic energy (FEM/StVK)](04_elastic_energy.md) | $F,\;E=\tfrac12(F^{\mathsf T}F-I),\;W$; forces & stiffness | ✅ built |
| 05 | [Statics in rig space](05_statics_in_rig_space.md) | minimize $W(s(p))$+gravity via Newton with $J^{\mathsf T}KJ$ | ✅ built |
| 06 | [Dynamics (implicit Euler)](06_dynamics_implicit_euler.md) | $\Phi=\tfrac1{2h^2}\text{inertia}+W$; secondary motion | ✅ built |
| 07 | [Nonlinear rig: 2012 vs 2013](07_nonlinear_rig_2012_vs_2013.md) | the $\partial^2 s/\partial p^2$ cost; linearize the rig | ✅ built |
| 08 | [Capstone: primary + secondary](08_capstone_primary_secondary.md) | driven keyframes + simulated jiggle | ✅ built |
| 09 | [Real FBX character](09_fbx_character.md) | the artist's skeleton as the rig: LBS + voxel proxy + 2013 solver | ✅ built |

The library adds `energy.py` (StVK) and `solver.py` (reduced Newton, statics,
implicit-Euler stepper, linearized-rig) for Lessons 4–8, then `fbx.py` (a pure
stdlib FBX reader + an LBS `SkeletonRig`) and `voxel.py` (winding-number
voxelizer + skin-weight transfer) for Lesson 9. Further extensions (interior DOFs
+ volumetric skinning $q=Ws$ on the real surface, TetGen/fTetWild meshing,
PCA-reduced parameters) are sketched at the end of Lessons 8–9.

## Verify the built lessons

```bash
# math only (fast, no GUI): e.g. tets oriented, J matches finite diff, Newton converges
DYN_NO_SHOW=1 python src/dyn/tutorials/05_statics_in_rig_space.py
# full GUI code path, headless, for every lesson:
for f in src/dyn/tutorials/0*.py; do DYN_SMOKE=1 python "$f"; done
```

## Why a beam and not the FBX in `data/`?

For *learning* the math, a procedural beam wins: its tets are generated
analytically and its rig is closed-form, so the exact Jacobian is known and
finite differences can be checked against it. The FBX characters need an FBX
parser **and** a tet mesher (neither installed) — plumbing unrelated to the
physics. They return as the optional capstone (Lesson 08) once the mechanics are
clear.
