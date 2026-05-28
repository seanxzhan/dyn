# Lesson 9 — Rig-space physics on a real FBX character

> Run `python src/dyn/tutorials/09_fbx_character.py`.
> Code: [`09_fbx_character.py`](09_fbx_character.py) · Library: [`../fbx.py`](../fbx.py), [`../voxel.py`](../voxel.py)

Everything in Lessons 1–8 used a rig we *invented*. The whole pitch of the papers
is that the rig is whatever the **animator already built** — a skeleton, a stack
of deformers, a black box. This lesson closes that gap: we load one of the
`data/*.fbx` characters, use its real skeleton as the rig, and run the exact same
secondary-motion solver. No code in `energy.py` or `solver.py` changes.

It's built from three dependency-free pieces (numpy only): an FBX reader, a
voxelizer, and a linear-blend-skinning rig.

## 1. What an FBX rig actually contains

[`fbx.py`](../fbx.py) parses the binary FBX (a tree of typed nodes; stdlib
`struct`+`zlib`) and pulls out:

- a **surface mesh** (`Geometry`: vertices + polygons), and
- a **skeleton**: `Model`/`LimbNode` bones, their parent hierarchy (from
  `Connections`), and each bone's **bind-pose world transform** $B_b$
  (the cluster's `TransformLink`), and
- **skin weights**: each `Deformer`/`Cluster` lists the vertices a bone
  influences and by how much.

For `2072.fbx` that's an 865-vertex fish with a 13-bone skeleton:
`origin → skl_root → c_body_1 → c_body_2 → c_body_3 → c_caudal_fin`, plus
pectoral/dorsal fins and a head/mouth chain.

## 2. The rig map: linear-blend skinning

The skeleton is the rig. Its parameters $p$ are **joint rotations**. Here each
free joint is a **1-DOF hinge** about the horizontal axis (world $x$) — the axis
perpendicular to the fish's length ($z$) and to gravity ($y$) — so the rear
droops cleanly in the vertical plane instead of buckling sideways. (Letting each
joint rotate in full 3-DOF, under a heavy load, lets it twist out of plane; a
hinge matches how a tail actually bends.) The map $s(p)$ is standard
**linear-blend skinning**:

$$ x_v(p) \;=\; \sum_b w_{vb}\;\, M_b(p)\, B_b^{-1}\, x^0_v, $$

where $x^0_v$ is the bind-pose vertex, $B_b$ the bone's bind world transform, and
$M_b(p)$ its *posed* world transform, obtained by walking the hierarchy
$M_b = M_{\text{parent}}\, L_b\, \Delta R_b(\omega_b)$ with bind-local
$L_b = B_{\text{parent}}^{-1}B_b$. At $p=0$ every $\Delta R=I$, so $M_b=B_b$ and
$x_v=x^0_v$ — the bind pose. The script asserts this `s(0) == mesh` to machine
precision; it is the check that the whole transform/bind bookkeeping is correct.

This `SkeletonRig` is a genuine black box: there is no tidy analytic Jacobian, so
`jacobian(p)` is taken by **finite differences** — exactly what Hahn et al. do
when their rig lives in Maya. Our toy `LinearRig`/`BendRig` were stand-ins for
*this*.

## 3. Getting a volume: voxelization (no mesher)

The elastic energy (Lesson 4) needs a *volumetric* tet mesh, but the FBX is only
a surface. [`voxel.py`](../voxel.py) builds the volume without any external
mesher: it samples a regular grid, tests each point for being inside the surface
with the **generalized winding number**

$$ w(q) = \frac{1}{4\pi}\sum_{\text{tris}} \Omega_{\text{tri}}(q)\;\;(\approx 1\ \text{inside},\ 0\ \text{outside}), $$

keeps the inside voxels, and splits each into 6 tets (the same Kuhn subdivision as
Lesson 1). The character's surface skin weights are then transferred to the proxy
vertices by an inverse-distance blend of the nearest surface vertices. The result
(≈1400 tets for the fish) is a blocky but honest volumetric stand-in that the
same skeleton drives.

## 4. The solve: 2013 secondary motion on the skeleton

We reuse `ImplicitEuler(relinearize=True)` — the **2013 linearized-rig** method
from Lesson 7, which is what makes a black-box LBS rig affordable (one
finite-difference Jacobian per step, no $\mathcal O(d^2)$ curvature term). The rig
parameters are split as in Lesson 8:

- **primary (driven):** the body bone `c_body_1` is held (or yaw-animated
  $A\sin\omega t$ if you turn on "swim");
- **secondary (free):** the rest of the spine and the captured fins are
  simulated.

By default `c_body_1` is fixed, so the fish is effectively **held by its
mid-body** and the rear half + tail **droop and jiggle under gravity** — exactly
the cantilever beam of Lesson 6, now a character. Released from the bind pose,
the tail falls ~14% of a body length, **overshoots, and rings down** (purely in
the vertical plane — the tail $x$ stays put). Turn on "swim" to add the keyframed
primary motion and watch the tail lag and whip instead.

**A note on scale.** This model is ~1 unit long and very light (proxy mass
≈0.02), so the elastic stiffness dwarfs ordinary gravity. We therefore use a
large "gravity load" (default 150) to bend it visibly — equivalently, a soft
material. Joint *rotations* are scale-invariant, so the physics is identical;
only the constant looks unusual. Treat the gravity slider as a load knob.

## A necessary safeguard

A free joint only has restoring stiffness if rotating it actually strains the
proxy's tets. A *thin* feature (a fin) may carry almost no proxy volume, giving a
near-zero-stiffness mode that blows up. So the lesson keeps a secondary bone only
if it carries enough proxy skin weight (`MIN_PROXY_WEIGHT`). At the default
`res=22` the body spine, caudal fin and dorsal fin all qualify; the wispy
pectoral fins are filtered out. (Raise `res` to include more, at the cost of
speed.)

## What to look at

- **Play.** Released from the bind pose, the rear half and tail **droop and
  jiggle**, falling ~14% of a body length, overshooting, then ringing down — the
  same behavior you saw on the beam. The orange curve traces the tail vertex: a
  clean vertical arc (no sideways swing).
- **rig-pose ghost.** The grey fish is the rig with all secondary parameters
  forced to $0$ (no physics). With the body held, that is just the straight bind
  pose; the gap between it and the solid fish is exactly what rig-space physics
  adds. (Earlier, with no gravity, the only effect was the tail *lagging*, so the
  physics fish looked *less* bent than the keyframe — correct, but undramatic.)
- **pluck tail.** Flick the tail with a velocity impulse and watch it oscillate
  back — the character version of plucking the beam.
- **show physics proxy.** Reveal the orange tet mesh the simulation actually runs
  on — the smooth fish is just that proxy's joint angles applied to the render
  mesh.
- **swim (A, ω).** Animate the body and the tail lags and whips; push ω toward
  the natural frequency for resonance. **gravity / stiffness / damping** trade off
  droop depth, frequency, and how long it rings.

## How faithful is this?

The structure is exactly the papers': a black-box rig, forces projected by
$J^{\mathsf T}$, implicit-Euler energy minimization in rig space, the
linearized-rig speedup, output as editable rig-parameter curves. The honest
simplifications here are the **coarse voxel proxy** (a real pipeline would tet-mesh
the actual surface with TetGen/fTetWild and keep the original surface as the
deforming domain), no collision handling, and our $q=Ws$ interior-skinning step
(Lesson 8's note) being skipped because the rig drives all proxy vertices. Those
are the natural next investments — but the physics you are watching is the real
method, on a real rig.

References: 2012 §3 (black-box rig, finite-difference Jacobian), 2013 §4.1
(linear rig) — see
[`../../../docs/Rig-Space-Physics.md`](../../../docs/Rig-Space-Physics.md),
[`../../../docs/Efficient-Secondary-Motion-Rig-Space.md`](../../../docs/Efficient-Secondary-Motion-Rig-Space.md).
