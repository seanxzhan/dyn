# Lesson 1 — Particles and the predict step

> Run `python src/dyn/tutorials/pbd/01_particles_and_predict.py`.
> Code: [`01_particles_and_predict.py`](01_particles_and_predict.py) · Library: [`../../pbd.py`](../../pbd.py)

This lesson has almost no physics. Its only job is to fix the data layout and
the **inner loop** every later lesson runs.

## 1. The PBD data: positions, velocities, inverse masses

A PBD body is just **N point masses**. We store

$$
x \in \mathbb{R}^{N\times 3}, \qquad
v \in \mathbb{R}^{N\times 3}, \qquad
w \in \mathbb{R}^{N},\quad w_i = 1/m_i.
$$

That's it. No tetrahedra, no FEM basis functions, no rig. Triangles and bend
pairs (Lesson 4) just record *which particles are constrained together* — the
particles themselves are still naked points.

**Inverse mass** $w_i$ instead of mass $m_i$ for one practical reason: pinning
a vertex is just $w_i = 0$. Every formula in PBD multiplies corrections by
$w_i$, so a pinned particle is automatically immovable, with no special-case
code path. We use this everywhere.

## 2. The simulation loop, half-built

The full algorithm (Müller §3.1) is the famous 8-line loop:

```
v ← v + h · w · f_ext            # gravity
p ← x + h · v                     # PREDICT next position
project constraints on p          # ← the heart of PBD (Lessons 2–5)
v ← (p − x) / h                   # recover velocity
x ← p                             # commit
```

**Today we run only lines 1, 2, 4, 5** — *no constraints*. With no constraints
to project, "predict" is just plain explicit Euler with $f = mg$, and the
velocity recovery $(p - x)/h$ trivially gives back what we already had.
Particles fall like points in vacuum.

The point of running this gutted version: every later lesson keeps the same
shell and only edits the projection step. Once you can read the code in this
lesson, you have read the outer loop of *every* PBD demo on the planet.

## 3. Why predict-then-project beats compute-forces-then-integrate

A classical force-based step looks like

$$
\text{accumulate } f_{\text{ext}}+f_{\text{int}} \;\to\; a=f/m \;\to\; v += a\,h \;\to\; x += v\,h .
$$

When something *constrains* the result — a wall, a rope, a pinned vertex —
forces have to predict the *future* state correctly so that velocities, when
integrated, end up where we want. Constraint forces are derivatives of an
implicit equation; getting them right is exactly the implicit Euler /
backward Euler problem (the topic of the rig-space-physics tutorial in
[`../rig-space-physics`](../rig-space-physics)).

PBD inverts the order. It first commits to a **predicted position** assuming no
internal constraints, then **moves the prediction** until constraints are
satisfied. Velocities are *recovered* from how far the projection moved each
particle. Two consequences:

- **Pinning**, **collision response**, **inextensibility** all become "snap the
  point to a valid place" — no force balance, no time-step shrinking.
- The integrator is **unconditionally stable**: line 5 commits to a
  configuration the projector already validated, never extrapolating into a
  diverging future.

The trade is that `solverIterations` is now a stiffness knob (Lesson 3) and
material stiffness is timestep-dependent — the canonical PBD flaw, fixed
later by XPBD.

## What to look at

- Hit **Play.** A 6×6 grid of particles falls under gravity. With no
  constraints, there is no cloth — just 36 independent points in freefall. The
  trail on one particle traces a parabola.
- Toggle **predict** off. The loop runs `v ← v + h·w·g`, then `x ← p` with
  `p = x + h·v_old`, but skipping the prediction means `p = x` and the
  particles never move (you've removed the only line that integrates).
- **gravity = 0.** The grid hangs in space. Without internal constraints,
  PBD has nothing to *do*.

## Where this is going

Lesson 2 puts the constraint projection in. Two particles connected by a
distance constraint will, under gravity, swing as a rigid pendulum — exactly
because the projector keeps the distance equal to its rest length. From there
every later constraint type is the same recipe with a different $C$.

Reference: Müller et al. 2006 §3.1 — see
[`../../../../docs/Position-Based-Dynamics.md`](../../../../docs/Position-Based-Dynamics.md).
