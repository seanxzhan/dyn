# Lesson 6 — Dynamics: implicit Euler as energy minimization

> Run `python src/dyn/tutorials/06_dynamics_implicit_euler.py`.
> Code: [`06_dynamics_implicit_euler.py`](06_dynamics_implicit_euler.py) · Library: [`../solver.py`](../solver.py)

This is the heart of the course. Statics found where the beam rests; now we add
**inertia** and watch it get there over time — overshooting and jiggling. That
transient is the secondary motion the papers exist to produce.

## 1. From Newton's law to a per-step minimization

Time-discretize $M\ddot{x} = f(x) = -\partial_x W$ with **backward (implicit)
Euler**. Writing the prediction $x^\* = 2x_n - x_{n-1}$ (where the beam would
drift by inertia alone, $=x_n + h\,v_n$), implicit Euler is

$$ M\,\frac{x_{n+1} - x^\*}{h^2} \;=\; -\,\partial_x W(x_{n+1}). $$

The magic: this is exactly the stationarity condition of an **incremental
potential**

$$ \Phi(x) = \underbrace{\frac{1}{2h^2}\,(x - x^\*)^{\mathsf T} M\,(x - x^\*)}_{\text{inertia}}
          \;+\; \underbrace{W(x)}_{\text{elastic}} \;-\; \sum_i m_i\,g\cdot x_i. $$

So **each timestep is a minimization** — the same kind we did for statics, plus
the inertia term that pulls $x$ toward the free-flight prediction $x^\*$. Implicit
Euler is unconditionally stable, which is why we can take big steps and still not
explode.

## 2. Doing it in rig space

Restrict to the rig, $x=s(p)$, and minimize $\Phi$ over $p$:

$$ \frac{\partial \Phi}{\partial p}
   = J^{\mathsf T}\!\Big[\tfrac{1}{h^2} M\,(s(p) - x^\*) + \partial_x W - f_{\text{grav}}\Big] = 0. $$

This is the 2012/2013 reduced equation of motion. It is literally the static
gradient (Lesson 5) **plus one inertial force term**, all projected by
$J^{\mathsf T}$. `ImplicitEuler` builds this $\Phi$ and its reduced gradient and
calls the same `newton_minimize` each step; afterwards it reads off the new
velocity and applies a little damping. Each frame of the viewer is one such
step.

## 3. Why this is "secondary motion"

A jiggle bone is a 1-DOF damped oscillator hung off the rig. What you are running
is the *same idea done properly*: a coupled, volumetric, implicitly-integrated
oscillation living in the rig's own coordinates. The beam doesn't snap to the
gravity pose — it falls past it and rings down, exactly like flesh or a tail
follows through after the primary motion. Compare:

- **jiggle bone:** $m\ddot{x} = -k(x-x_{\text{tgt}}) - c\dot x$, per bone, decoupled;
- **here:** $M\ddot x = -\partial_x W - $ damping, with $W$ a *coupled* continuum
  energy, solved in rig space.

## What to look at

- **Play.** Released flat, the beam sags, **overshoots past** the green ghost
  (the static equilibrium from Lesson 5), swings back, and rings down onto it.
  The orange curve is the tip's path — a decaying oscillation.
- **physics on/off.** Off, the beam just stays at its rest rig pose (no sag at
  all). The *difference* between the two is precisely what the dynamics adds.
- **pluck.** Flick the tip: a velocity impulse re-excites the oscillation. Watch
  it ring down again to the same equilibrium.
- **damping.** Near $0$ the beam oscillates for a long time (underdamped); larger
  values settle quickly. **Stiffness $E$:** stiffer → higher frequency, smaller
  sag. **Gravity:** deeper sag, the ghost moves with it.

## A note on "primary + secondary"

Here gravity is the only forcing, so the whole motion is "secondary". In
production the rig is split: an animator keyframes some parameters (the
**primary** animation) and physics adds **secondary** motion on the rest. That is
exactly the capstone (Lesson 8): we drive one rig parameter and let the inertial
coupling make the others jiggle.

## Where this is going

Everything so far used a **linear** rig, where $J$ is constant. Lesson 7 switches
to the **nonlinear** `BendRig` and asks: do we really need the expensive rig
curvature term in the Hessian? Answer (the 2013 contribution): linearize the rig
once per step and the motion is essentially unchanged — for a fraction of the
cost.

Reference: 2012 eq. 2–3 (incremental potential), 2013 eq. 3 — see
[`../../../docs/Rig-Space-Physics.md`](../../../docs/Rig-Space-Physics.md).
