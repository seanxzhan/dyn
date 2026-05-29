# Lesson 5 — Collisions: inequality constraints, friction, restitution

> Run `python src/dyn/tutorials/pbd/05_collisions_and_response.py`.
> Code: [`05_collisions_and_response.py`](05_collisions_and_response.py) · Library: [`../../pbd.py`](../../pbd.py)

A cloth falls onto a static sphere. Same projection loop as Lesson 4, with
**one new constraint**: each particle must stay outside the sphere. The
twist — collisions are *inequalities*, only active when violated.

## 1. Equality vs inequality constraints

Every constraint we've seen so far is an **equality**: distance edges *want*
$|p_a - p_b| = \ell_0$, dihedrals *want* the rest angle. The projection
formula is symmetric — overshoot and undershoot are corrected the same way.

Collisions are different. A particle near the sphere is fine; a particle
*inside* the sphere is not. So we project only when violated:

$$
C_{\text{coll}}(p) \;=\; \|p - c\| - r \;\geq\; 0
\quad\Longrightarrow\quad
\text{project iff } \|p - c\| < r.
$$

When violated, the gradient is just the outward normal $\hat n = (p-c)/\|p-c\|$
and the cleanest projection is to push $p$ exactly to the surface:

$$
p \;\leftarrow\; c + r\,\hat n.
$$

This is `project_sphere` in [`../../pbd.py`](../../pbd.py). It returns the
boolean mask of particles that hit so the velocity step can apply
material-specific response.

## 2. Where in the loop?

Collisions go **inside the projection loop** alongside the cloth constraints,
so each iteration sees the latest positions:

```
predict  →  p
for it in range(n_iters):
    project_distance(p, ...)    # equality
    project_bend(p, ...)        # equality
    project_sphere(p, ...)      # INEQUALITY — only when violated
v ← (p − x) / h
```

Putting collisions **last** in the iteration is a deliberate choice: stretch
and bend may push a particle into the sphere; the sphere projection then
pushes it back out. On the next iteration, stretch/bend may again pull it
inward; the sphere kicks again. The system converges to "particle on
surface, fabric tensioned tangentially around it" — exactly what cloth
draped on a ball should look like.

## 3. Friction & restitution: a velocity post-process

PBD recovers velocity at the end of the step from positions: $v = (p - x)/h$.
Naively, this gives a **perfectly inelastic, frictionless** collision: the
particle ends up on the sphere surface with whatever tangential velocity it
had, no bounce. To get realistic behaviour, decompose $v$ at each colliding
particle into normal and tangential components and modify them:

$$
v_n = (v\cdot\hat n)\,\hat n,\qquad v_t = v - v_n,
$$

then

$$
\boxed{\;
v \;\leftarrow\;
\underbrace{(1 - \mu)\,v_t}_{\text{Coulomb-like friction}}
\;-\;\underbrace{e\,v_n}_{\text{restitution}}
\;}
$$

- $\mu \in [0, 1]$ — friction coefficient. $\mu = 0$ is ice (cloth slides
  off), $\mu = 1$ is grip (cloth sticks tangentially).
- $e \in [0, 1]$ — restitution. $e = 0$ is "stick to it", $e = 1$ is
  "perfectly bouncy". For cloth, you almost always want $e = 0$ (cloth
  doesn't bounce off things; it drapes).

This *isn't* the rigorous Coulomb cone, just the cheap version that looks
right. PBD's contributions are positional; the velocity post-process is just
"don't violate basic intuition".

## 4. Collision thickness

A common practical trick: project particles to $r + \tau$ rather than $r$,
where $\tau$ is a small "thickness". This prevents the cloth from
**z-fighting** the sphere visually (you can see the surface mesh punching
through the sphere mesh) and absorbs small numerical jitter. In real
production cloth, $\tau$ is also where you encode garment offset (a sweater
sits *off* the body, not *on* it).

## What to look at

- **drop the cloth.** Hit `release` to unpin the top row. The cloth falls,
  hits the sphere, drapes around it — pure equality + inequality
  constraints, no special-case code.
- **friction slider.** With $\mu = 0$ the cloth slides right off the sphere
  into a heap. With $\mu = 1$ it sticks where it lands. Most "real fabric"
  feels are in the $0.2$ – $0.6$ range.
- **bend stiffness.** Now you can see why bend matters for collisions: with
  bend $= 0$, cloth wraps the sphere so tightly it forms cusps; turn bend up
  and it stays smoother, more "canvas-like" against the curvature.
- **stretch slider.** Drop it low, and the cloth visibly stretches around
  the sphere — gravity pulls the unsupported part down while the supported
  part is fixed by the contact.

## Where this is going

Lesson 6 (capstone) makes the sphere **dynamic**: it has mass, falls under
gravity, and the cloth can push it back. That requires *one* more change —
when projecting cloth particles outside the sphere, also accumulate the
opposite impulse on the sphere itself. Then "throw a ball through cloth"
is just a button that sets the sphere's initial velocity.

Reference: Müller et al. 2006 §4.3 (cloth-rigid collisions) — see
[`../../../../docs/Position-Based-Dynamics.md`](../../../../docs/Position-Based-Dynamics.md).
