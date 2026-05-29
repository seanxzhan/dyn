# Lesson 6 — Capstone: a ball thrown at a hanging cloth

> Run `python src/dyn/tutorials/pbd/06_capstone_ball_through_cloth.py`.
> Code: [`06_capstone_ball_through_cloth.py`](06_capstone_ball_through_cloth.py) · Library: [`../../pbd.py`](../../pbd.py)

The full demo. A pinned cloth hangs like a curtain; a dynamic sphere with
mass falls (or is launched) into it. Stretch + bend + collisions + two-way
momentum coupling — every piece from Lessons 1–5, in the same loop.

## 1. Two-way coupling — what changes

Lesson 5's sphere was static: cloth particles got pushed out, the sphere
didn't notice. Now the sphere has mass $M$ and velocity $V$, so when the
cloth pushes a particle out by $\Delta p_i$, the sphere must receive an
**equal and opposite impulse**:

$$
J_{\text{sphere}} \;=\; -\sum_{i\,\in\,\text{contacts}} m_i \,\frac{\Delta p_i}{h}.
$$

(Mass × Δposition / step = impulse. The sign is negative because the cloth's
projection vector $\Delta p_i$ points *outward* from the sphere; the
reaction on the sphere points *inward*.)

This is straight Newton's third law expressed in the position-displacement
language of PBD. After all cloth-vs-sphere projections, we update the
sphere's velocity:

$$
V \;\leftarrow\; V + \frac{J_{\text{sphere}}}{M}.
$$

`project_dynamic_sphere` in [`../../pbd.py`](../../pbd.py) does exactly this:
it projects each colliding particle to the sphere surface (Lesson 5) **and**
accumulates the corresponding impulse on the sphere.

## 2. The full step

```
                                    cloth                   ball
                                    ─────                   ────
predict      v ← v + h·g            (per particle)          V ← V + h·g
             p ← x + h·v            (predicted)             P ← C + h·V
                                                            (predicted center)

iterate                                                     ┌── these stay in
  project_distance(p, ...)                                  │   sync because
  project_bend(p, ...)                                      │   they share P:
  project_dynamic_sphere(p, ...) ────► impulse on ball ◄────┤   each iter sees
                                                            │   the ball's new
                                                            │   position too
                                                            └──

commit       v ← (p − x) / h        (cloth velocity)        V ← V + impulse/M
             x ← p                  (cloth position)        C ← P + h·(impulse/M)
```

The sphere's centre is **fixed during the iteration** (we use its predicted
$P = C + h V$ for collisions). Apply the impulse and integrate the centre
*after* the iteration — Müller §4.2's "shock-propagation" approach. This
keeps the iteration cheap (no per-iter rebuild of sphere geometry) and
preserves linear momentum to the precision of the impulse sum.

(Angular momentum on the sphere is *not* preserved here — we treat it as a
point mass. For a billiard-ball or rolling-ball look you'd accumulate
torque from off-centre contacts and integrate angular velocity. That's a
useful add-on, not a foundational change.)

## 3. What the iteration count buys you, again

This is where iteration count earns its keep. Per Lesson 3, information
travels one constraint per sweep. With a fast-moving ball:

- **few iterations**: contact "wave" hasn't propagated outward through the
  cloth yet, so the ball's nose pushes a few neighbours forward while the
  rest of the cloth doesn't react. Stretch ratios near contact spike, the
  cloth visibly distends, and the ball squeezes through.
- **many iterations**: each step, the impulse propagates several segments
  outward; the cloth tensions tangentially around the impact, distributes
  the load, and stops the ball.

The same trade-off cloth solvers have argued about for 30 years, distilled
into a slider.

## 4. Sub-stepping — the trick that saves the demo

For $V \approx 5$ m/s and $h = 0.016$ s, the ball moves $\sim 8$ cm per step.
If the cloth's grid spacing is $\sim 5$ cm, the ball can cross a whole
quad in one step — and PBD's projection only sees the *committed* state, so
mid-step tunneling goes unprojected.

The standard PBD fix is **sub-stepping**: take $k$ inner steps of size
$h/k$. Each inner step still does a full predict-project-commit, so contact
is rechecked $k$ times per visible frame. $k = 2$ – $4$ is usually enough.
The visible difference between $k = 1$ and $k = 4$ on a fast throw is
"ball goes through cloth" vs "ball gets caught" — even with the same
iteration count.

(XPBD goes further: it explicitly couples the sub-stepping to compliance,
so the *physical* stiffness is preserved as $h$ changes. PBD has no such
guarantee, which is why these demos always need tuning.)

## 5. What to look at

- **drop test.** Hit `reset` then `release ball`. The ball falls onto the
  draped portion, and the cloth catches it like a hammock. With high bend
  stiffness it stays nearly flat; with low bend it pockets around the ball.
- **throw test.** Hit `throw` to give the ball forward velocity into the
  cloth's face. Vary speed (slider) and watch the spectrum:
  - low speed: the cloth catches it, swings forward, swings back.
  - medium speed: the cloth deforms a lot, the ball squeezes through if
    iterations / sub-steps are too few.
  - high speed: tunneling — the ball pops out the back. Crank sub-steps to
    fix it without changing $h$.
- **iteration vs sub-step.** Try (iterations = 4, sub-steps = 1) vs
  (iterations = 1, sub-steps = 4). Same total work, but the second tunnels
  *less* on a fast throw — the bottleneck is collision detection frequency,
  not constraint convergence.
- **mass ratio.** A 0.1 kg ball into 1 kg cloth bounces; a 10 kg ball just
  punches through and drags the cloth with it. The impulse formula
  $\Delta V = J/M$ does the work — heavier ball, smaller velocity change
  per contact.

## 6. Where this stops being PBD

This demo is the natural endpoint of "small NumPy PBD": it does the headline
thing — interactive cloth — with about 350 lines of physics code total.
What it doesn't do, and where the field has gone since 2006:

- **XPBD** (Macklin et al. 2016) — adds compliance $\alpha = 1/(k\,h^2)$ to
  every constraint; the $k'$ trick goes away and stiffness becomes
  timestep-independent in the limit. Drop-in replacement.
- **Self-collision and CCD** — cloth-vs-itself, broadphase via spatial hash,
  continuous collision tests. The current demo will let cloth pass through
  itself.
- **Strain limiting** — clamp stretch ratio per edge as a hard cap. PBD's
  $k'$ scheme leaves residual stretch; for hero garments you want a
  guaranteed bound. (Note: this is *limiting*, not *breaking* — the edge
  stays connected, just refuses to stretch past a threshold.)
- **Reduced cloth** (Hahn et al. 2014, "Subspace Clothing Simulation") —
  same sub-space trick rig-space-physics uses, applied to cloth. The cloth
  mesh becomes a thin shell driven by a body; PBD updates a tiny set of
  basis coefficients instead of every vertex.

The throughline: positions over forces, projection over differentiation,
and "constraint" as the lingua franca. Once you've seen it, every
follow-up paper is an extension or a refinement of the same loop.

Reference: Müller et al. 2006 §4.2 (rigid-body interaction) — see
[`../../../../docs/Position-Based-Dynamics.md`](../../../../docs/Position-Based-Dynamics.md).
