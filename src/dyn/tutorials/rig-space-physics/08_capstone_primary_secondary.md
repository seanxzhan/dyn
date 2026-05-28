# Lesson 8 (capstone) — Secondary motion on a primary animation

> Run `python src/dyn/tutorials/08_capstone_primary_secondary.py`.
> Code: [`08_capstone_primary_secondary.py`](08_capstone_primary_secondary.py) · Library: [`../solver.py`](../solver.py)

This is the scenario the papers are actually built for, and it ties together
everything from Lessons 1–7.

## 1. Primary vs. secondary

An animator keyframes a few rig parameters — the **primary** animation. Physics
should add **secondary** motion (jiggle, lag, follow-through) to the *rest* of the
rig, automatically, and return it as more rig-parameter curves the animator can
still edit. So we split the rig parameters:

$$ p = (\,p_{\text{driven}},\; p_{\text{free}}\,), $$

prescribe $p_{\text{driven}}(t)$ from the keyframes, and solve the rig-space
dynamics (Lesson 6) over $p_{\text{free}}$ **only**:

$$ p_{\text{free}}^{\,n+1} = \arg\min_{p_{\text{free}}}\;
   \Phi\big(p_{\text{driven}}^{\,n+1},\, p_{\text{free}}\big). $$

The coupling is inertial: $\Phi$'s prediction $x^\*=2x_n-x_{n-1}$ carries the
momentum of the driven motion, so accelerating the driven part exerts an inertial
force on the free part. (`ImplicitEuler(free_idx=...)` minimizes over the free
indices each step while `step(driven_values=...)` injects the primary.)

## 2. The demo

- **Primary:** `lift_z` — a rigid up/down translation of the whole beam, driven by
  $p_0(t)=A\sin(\omega t)$. Think of grabbing the clamp and shaking it. Being a
  rigid mode it stores no elastic energy on its own.
- **Secondary:** `bend_z` — the free parameter. As the handle accelerates, the
  beam's mass makes the tip lag and overshoot; the elastic energy provides the
  restoring force. The result is a bend that oscillates around the primary motion
  — secondary motion, expressed as a rig-parameter curve $p_1(t)$.

The **grey ghost** is the rig with $p_1=0$: where the beam would be with the
primary animation alone and no physics. The gap between the solid beam and the
ghost *is* the secondary motion.

## What to look at

- **physics on/off.** Off, the solid beam coincides with the ghost — it just
  rides the handle rigidly. On, the tip lags and whips; the orange trail shows the
  loopy path inertia produces.
- **Sweep ω (frequency).** Low $\omega$: the beam quasi-statically follows the
  handle, little lag. As $\omega$ approaches the beam's **natural frequency** the
  secondary bend grows dramatically — resonance — then falls off above it. The
  console's "peak secondary bend" quantifies this; the slider lets you hunt the
  resonance live.
- **damping / stiffness / amplitude.** Damping tames resonance; a stiffer beam
  (larger $E$) raises the resonant frequency and shrinks the lag; bigger $A$
  scales the whole effect. Add **gravity** to also sag while it shakes.

## What you've built

Across the course you implemented, from scratch and verified against finite
differences, the complete rig-space physics method:

1. a tet-discretized body with a lumped mass and a surface/interior split (L1);
2. a rig $s(p)$ and its Jacobian $J=\partial s/\partial p$, with $J^{\mathsf T}$
   projecting forces into rig space (L2–L3);
3. StVK elastic energy $W$, its forces and stiffness (L4);
4. statics as $\min_p W(s(p))+$gravity, by Newton with $J^{\mathsf T}KJ$ (L5);
5. dynamics as per-step minimization of the incremental potential — secondary
   motion (L6);
6. the linear-rig shortcut that drops the $\mathcal O(d^2)$ curvature term, at
   imperceptible quality cost (L7);
7. driven primary + simulated secondary, the production loop (here).

## Honest next steps (beyond this course)

- **Free interior DOFs + volumetric skinning $q=Ws$** — the 2013 paper's other
  big idea: keep interior nodes as unknowns (richer dynamics) but eliminate them
  with a learned linear map fit from "shaking" examples. Here the rig drove *all*
  vertices, sidestepping this.
- **A real character** — load one of `data/*.fbx` (skeleton + skin), tetrahedralize
  the volume, and use the skeleton as the rig. This needs an FBX loader and a tet
  mesher (e.g. TetGen / fTetWild), which we deliberately avoided to keep the math
  in closed form. This is the natural capstone-of-the-capstone.
- **Reduced rig parameters (PCA)** and **per-parameter material control** — the
  2012 paper's §5–6 extensions.

Reference: 2012 §6 (interior secondary-motion toggle), 2013 §4.2 (skinning) —
[`../../../docs/Rig-Space-Physics.md`](../../../docs/Rig-Space-Physics.md),
[`../../../docs/Efficient-Secondary-Motion-Rig-Space.md`](../../../docs/Efficient-Secondary-Motion-Rig-Space.md).
