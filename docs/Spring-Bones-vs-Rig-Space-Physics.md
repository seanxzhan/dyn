# Spring Bones / Jiggle Bones vs. Rig-Space Physics

> How the ubiquitous real-time "jiggle bone" technique relates to the two rig-space physics papers in this
> folder: [Rig-Space Physics (2012)](./Rig-Space-Physics.md) and
> [Efficient Simulation of Secondary Motion in Rig-Space (2013)](./Efficient-Secondary-Motion-Rig-Space.md).

## TL;DR

They solve the **same problem** — add physically-plausible secondary motion to keyframed animation, kept
in an editable rig representation — at **very different points on the fidelity/generality/cost spectrum**.
Jiggle bones are the cheap, decoupled, hand-tuned, real-time approximation; rig-space physics is the
coupled, continuum-elastic, automatic, rig-agnostic (but heavier) version. The papers' energy
`H = inertia + elastic W` is the continuum, coupled generalization of the jiggle bone's
`½mv² + ½kx²`.

## What spring / jiggle bones are

A jiggle (a.k.a. spring / wiggle) bone is an **extra bone the rigger adds by hand** that is *not*
keyframed. It is driven by a **lumped mass-spring-damper** that lags behind its parent's motion —
essentially a damped harmonic oscillator:

```
m·a = −k·(x − x_target) − c·v + (drag from parent motion) + gravity
```

- `x_target` = where the bone "should" be per the primary animation
- `k` = stiffness, `c` = damping, `m` = mass

Properties:

- Integrated per-frame with **explicit / semi-implicit Euler or Verlet**.
- Each bone solved **independently** (or in simple parent→child chains) — no global coupling.
- Bone transforms drive the mesh via ordinary **linear blend skinning**.
- Hand-tuned parameters: bounciness, speed, gravity, angle limits.
- Universal engine support: Unity `ConfigurableJoint`, Unreal PhAT, Godot `SpringJoint3D`,
  Source `$jigglebone`.

## The shared core idea

Both jiggle bones and the two papers **simulate secondary motion inside a low-dimensional subspace tied to
the rig**, and output animation that lives in the rig representation so artists can still edit it. Neither
dumps raw per-vertex simulation onto the animator. Jiggle bones are the decades-old, real-time
approximation of exactly the rig-space physics thesis.

The papers' incremental potential

```
H = (h²/2)·(inertia term)  +  W(elastic energy)
```

is the **continuum, coupled generalization** of a jiggle bone's `½mv² + ½kx²`. A single jiggle bone is a
1-DOF, linearized, decoupled special case of the same inertia-plus-restoring-force minimization the papers
solve.

## How they differ

| Aspect | Jiggle / spring bones | Rig-Space Physics (2012) | Efficient Secondary Motion (2013) |
|---|---|---|---|
| **Subspace** | extra bones added **by hand** | the **existing arbitrary rig** (black box: skeletons, blendshapes, cages, FFD) | same |
| **Dynamics model** | per-bone mass-spring-damper (ad hoc `k`, `c`) | continuum FEM elasticity (StVK / Neo-Hookean) on a volumetric mesh | same energy `W` |
| **Coupling** | **none** — each bone independent | **full** coupling via elastic energy + volume preservation | full, via skinned interior |
| **Material meaning** | spring constants, no real material | true Young's-modulus-like material; per-parameter stiffness control (§5) | inherits FEM, but **drops** material control |
| **Volume / flesh bulging** | no | yes | yes |
| **Integration** | explicit / Verlet, one-way driven | implicit Euler as optimization (Newton + BFGS) | implicit Euler + **linearized rig** + skinning |
| **Cost** | trivially real-time | seconds/frame (offline) | ~<1 s/frame (near-interactive) |
| **Generality** | only the bone abstraction | **any** rig | any rig |

## The most useful framing

The two papers are the **principled, general, automatic version of jiggle bones**:

- **Jiggle bones** = a *manually authored* subspace (you place the bones) + *decoupled, linear* lumped
  springs + explicit integration. Great for hair, cloth, tails, and appendages where coupling and volume
  don't matter.
- **Rig-space physics** = the *rig you already have* as the subspace + a *coupled continuum elastic* model
  + implicit optimization. Captures flesh jiggle, bulging, volume preservation, and cross-parameter
  coupling that independent springs fundamentally cannot — and works with deformers that aren't bones.

Two sharp connections:

1. **The 2013 sumo with 174 "secondary-motion rig parameters"** is conceptually *"strap on lots of jiggle
   handles"* — except they are driven by a single **coupled FEM elastic solve** with volume preservation,
   not 174 independent springs. That example competes directly on jiggle-bones' home turf (many secondary
   DOFs, near-interactive) while keeping physical coupling.
2. **The 2013 "physics-based volumetric skinning" `q = W·s`** uses the *same machinery* as bone skinning —
   sparse linear blend weights with partition-of-unity — but points it the opposite way: surface→interior,
   purely to eliminate interior DOFs. Same tool, different job.

The lineage is still active: [*Real-Time Secondary Animation with Spring Decomposed Skinning* (CGF
2025)](https://onlinelibrary.wiley.com/doi/10.1111/cgf.70209) decomposes skinning into springs for
real-time secondary motion — sitting between the jiggle-bone and rig-space-physics worlds.

## When each wins

- **Jiggle bones:** real-time, trivial to implement, universal engine support, ideal for loosely-coupled
  appendages. Cost: hand-tuned, decoupled, no volume preservation, can't model coupled flesh deformation,
  tied to the bone abstraction.
- **Rig-space physics:** physically correct coupled/volumetric deformation, automatic, rig-agnostic,
  editable output. Cost: heavier compute, needs a tet mesh and material tuning.

## Sources

- [Wayline — Jiggle Physics implementation guide (mass-spring-damper)](https://www.wayline.io/blog/jiggle-physics-implementation-guide)
- [naelstrof/JigglePhysics (Unity addon)](https://github.com/naelstrof/JigglePhysics)
- [Valve Developer Community — `$jigglebone`](https://developer.valvesoftware.com/wiki/$jigglebone)
- [Real-Time Secondary Animation with Spring Decomposed Skinning (CGF 2025)](https://onlinelibrary.wiley.com/doi/10.1111/cgf.70209)
