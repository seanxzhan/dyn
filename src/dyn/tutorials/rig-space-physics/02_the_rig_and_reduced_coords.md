# Lesson 2 — The rig and reduced coordinates

> Run `python src/dyn/tutorials/02_the_rig_and_reduced_coords.py`.
> Code: [`02_the_rig_and_reduced_coords.py`](02_the_rig_and_reduced_coords.py) · Library: [`../rig.py`](../rig.py)

## 1. What a rig *is*, mathematically

An animator never moves 528 vertex coordinates by hand. They move a few
**rig parameters** $p \in \mathbb{R}^d$ — joint angles, blendshape weights, cage
handles — and the rig expands those into a full deformation. Formally a rig is a
map

$$ s : \mathbb{R}^d \longrightarrow \mathbb{R}^{3N_s}, \qquad p \mapsto s(p), $$

from a *few* parameters to the *many* surface vertices. In the papers this map is
a **black box**: it could be a skeleton with linear-blend skinning, a stack of
blendshapes, a cage with harmonic coordinates, a lattice FFD. Rig-space physics
treats it abstractly — all it ever needs from the rig is to *evaluate* $s(p)$ and
(Lesson 3) to differentiate it.

The set of all shapes the rig can produce,

$$ \mathcal{M} = \{\, s(p) : p \in \mathbb{R}^d \,\} \subset \mathbb{R}^{3N_s}, $$

is **rig space** — a $d$-dimensional surface ("manifold") sitting inside the huge
$3N_s$-dimensional space of all possible vertex configurations. The whole point
of the method: *run the physics inside $\mathcal{M}$*, solving for $p$ (a handful
of numbers) instead of for every vertex. $p$ are the **reduced coordinates**.

## 2. The simplest rig: linear

The cleanest possible rig is **linear** (a blendshape rig):

$$ s(p) = s_0 + B\,p, \qquad B \in \mathbb{R}^{3N\times d}. $$

$s_0$ is the rest shape and each column of $B$ is one **mode** — a fixed
displacement field added in proportion to its parameter. Our `LinearRig` ships
two modes, each with a profile $(x/L_x)^2$ along the beam:

- `bend_z` — the tip lifts in $+z$,
- `sway_y` — the tip swings sideways in $+y$.

The quadratic profile vanishes *with zero slope* at the clamp ($x=0$), so the
clamped end stays put and the deformation grows smoothly toward the tip — exactly
the shape a clamped beam wants. Here rig space $\mathcal{M}$ is just the affine
plane $s_0 + \mathrm{span}(B)$: a flat 2-D slice through configuration space.

Why start linear? Because it makes rig space a literal flat subspace and the
Jacobian (next lesson) a *constant*. That is precisely the local picture the 2013
paper engineers on purpose — it linearizes any rig to $s(p_n) + J\,(p - p_n)$ once
per timestep. Master the linear case and you understand the regime the fast
method runs in.

## 3. The connection to jiggle bones

A classic "jiggle bone" is the same idea stripped to one decoupled spring: a
hand-placed extra bone whose 1-D offset is a reduced coordinate, driving the mesh
through linear-blend skinning. Rig-space physics generalizes this — it reuses the
*whole rig* as the subspace and couples everything through real elasticity
instead of one independent spring per bone. See
[`../../../docs/Spring-Bones-vs-Rig-Space-Physics.md`](../../../docs/Spring-Bones-vs-Rig-Space-Physics.md).

## What to look at

- Drag `p[0]` and `p[1]`. **Two numbers** smoothly reshape all 176 vertices.
  Every shape you can make is a point of the 2-D rig space; you *cannot* leave it
  no matter how you set the sliders — there is no slider that, say, dimples one
  face, because that motion is not in $\mathrm{span}(B)$.
- The clamped end stays glued while the tip moves most — the $(x/L_x)^2$ profile.
- `reset` returns to $p=0$, i.e. $s(0)=s_0$ (the script asserts this).

## The takeaway for the physics

Optimizing over $p \in \mathbb{R}^2$ instead of $x \in \mathbb{R}^{528}$ is what
makes rig-space physics tractable and keeps the output editable (it *is* rig
parameters). The price is that motion is confined to $\mathcal{M}$ — fidelity is
only as rich as the rig. Next we differentiate the map: the rig Jacobian
$J=\partial s/\partial p$ is the bridge that lets forces in vertex space talk to
the reduced coordinate $p$.
