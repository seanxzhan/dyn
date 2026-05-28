# Lesson 5 — Statics in rig space

> Run `python src/dyn/tutorials/05_statics_in_rig_space.py`.
> Code: [`05_statics_in_rig_space.py`](05_statics_in_rig_space.py) · Library: [`../solver.py`](../solver.py)

We now do our first real solve: where does the beam *settle* under gravity, given
that it can only move within rig space? This is statics — and it is just "find
the bottom of an energy landscape", restricted to the rig parameters.

## 1. The energy to minimize

Add gravity (an external potential $W_{\text{ext}}=-\sum_i m_i\,g\cdot x_i$) to
the elastic energy and substitute the rig $x=s(p)$:

$$ E(p) \;=\; W\!\big(s(p)\big) \;-\; \sum_i m_i\,g\cdot s_i(p). $$

The physical equilibrium is the **minimizer** $p^\star=\arg\min_p E(p)$: the pose
where elastic restoring forces exactly balance gravity, *within what the rig can
express*. (Full-space physics would minimize over all vertices; rig-space
physics minimizes over the handful of $p$, and the answer stays an editable rig
pose.)

## 2. First-order condition — forces projected by $J^{\mathsf T}$

Setting $\partial E/\partial p = 0$ and using the chain rule $\partial_p =
J^{\mathsf T}\partial_x$ (Lesson 3):

$$ \frac{\partial E}{\partial p}
   \;=\; J^{\mathsf T}\!\Big(\underbrace{\partial_x W}_{\text{elastic}}
   \;-\; \underbrace{f_{\text{grav}}}_{m_i g}\Big) \;=\; 0. $$

This is exactly the papers' static equilibrium condition: take the ordinary
vertex-space forces and **project them into rig space with $J^{\mathsf T}$**; at
equilibrium the projected force vanishes. (`make_static_problem` builds this
$E$ and its reduced gradient.)

## 3. Solving it — Newton in rig space

Newton's method needs the reduced gradient $g_p$ and the reduced **Hessian**

$$ H_{pp} = \frac{\partial g_p}{\partial p}
          = J^{\mathsf T} K\, J \;+\; (\partial_p J)\,(\partial_x W - f_{\text{grav}}), $$

then steps $\;p \leftarrow p + \alpha\,\Delta p,\;\; H_{pp}\Delta p = -g_p\;$ with
a line search. Two notes:

- For a **linear** rig $\partial_p J = 0$, so $H_{pp}=J^{\mathsf T}KJ$ exactly —
  the clean projected stiffness.
- We form $H_{pp}$ by **finite-differencing the reduced gradient** over $p$
  (`numerical_jacobian` in [`../solver.py`](../solver.py)). Since $p$ has only 2
  entries this is 4 evaluations, and — crucially — it captures the curvature term
  $\partial_p J$ automatically. (That convenience is exactly the $\mathcal{O}(d^2)$
  cost the 2012 paper pays for a black-box rig, and the subject of Lesson 7.)

The console prints the iterates: energy drops and $\lvert g_p\rvert\to 0$ in
about 5 steps — Newton's quadratic convergence near the minimum.

## What to look at

- **The bowl.** Beside the beam is the energy landscape $E(p)$ over the two rig
  parameters, drawn as a height field colored by energy. Statics = roll to the
  bottom. The **orange path** is the Newton trajectory; the **green dot** is the
  current iterate.
- **Scrub "Newton iteration".** The beam morphs from rest ($p=0$, the grey ghost)
  down to the sagged equilibrium as the green dot descends into the bowl. You are
  watching the optimizer and the pose at the same time.
- **Gravity / Young's modulus sliders.** More gravity (or a softer material =
  smaller $E$) deepens the bowl and moves its minimum — the beam sags further.
  The bowl is steeper for stiffer materials (bigger curvature $K$).
- Note the minimum sits at $p^\star\approx(-0.5,\,0)$: it sags in $-z$ and, by
  symmetry, not at all sideways.

## Where this is going

Statics found *where* the beam ends up. Lesson 6 adds **inertia** and asks how it
gets there over time — the same energy plus a kinetic term, minimized once per
timestep. The beam will overshoot this equilibrium and jiggle around it: that
transient is the secondary motion we are after, and it settles onto exactly the
$p^\star$ you found here.

Reference: 2012 §3 (first-order conditions, Newton, Schur complement) — see
[`../../../docs/Rig-Space-Physics.md`](../../../docs/Rig-Space-Physics.md).
