# Lesson 7 — Nonlinear rig: the 2012 curvature term vs. the 2013 linearization

> Run `python src/dyn/tutorials/07_nonlinear_rig_2012_vs_2013.py`.
> Code: [`07_nonlinear_rig_2012_vs_2013.py`](07_nonlinear_rig_2012_vs_2013.py) · Library: [`../solver.py`](../solver.py)

Until now the rig was linear, so its Jacobian $J$ was constant and the second
derivative $\partial^2 s/\partial p^2$ was exactly zero. Real rigs (skeletons,
deformers) are **nonlinear** — and that one fact is the difference between the
2012 paper and its much faster 2013 successor. This lesson makes it concrete with
the `BendRig`.

## 1. What the curvature costs (2012)

Newton on the incremental potential needs the reduced **Hessian**
$\partial^2\Phi/\partial p^2$. Differentiating the reduced gradient
$g_p = J^{\mathsf T}\!\big[\tfrac1{h^2}M(s-x^\*) + \partial_x W - f_{\text{grav}}\big]$
by the product rule gives two kinds of terms:

$$ H_{pp} = \underbrace{J^{\mathsf T}\Big(\tfrac1{h^2}M + K\Big)J}_{\text{projected mass + stiffness}}
   \;+\; \underbrace{(\partial_p J)\;\big[\tfrac1{h^2}M(s-x^\*) + \partial_x W - f_{\text{grav}}\big]}_{\text{rig curvature} \;\times\; \text{residual force}}. $$

That second term contains $\partial_p J = \partial^2 s/\partial p^2$. For a
black-box rig with $d$ parameters you get it by finite-differencing $J$ — which
itself costs $\mathcal O(d)$ rig evaluations — so the curvature term costs
$\mathcal O(d^2)$ rig evaluations **per Newton iteration**. For the 2012 paper's
36–174-parameter rigs that is the dominant cost and the reason it is far from
interactive.

(In our solver this term is hidden inside `numerical_jacobian(grad)`: because the
true `BendRig.s` is nonlinear, finite-differencing the reduced gradient over $p$
automatically picks up the curvature. Convenient — and exactly the expensive
thing.)

## 2. The 2013 fix: linearize the rig once per step

Approximate the rig by its tangent at the start of the step (`LinearizedRig`):

$$ s(p) \;\approx\; s(p_n) + J(p_n)\,(p - p_n). $$

Now $J$ is **constant within the step**, so $\partial_p J \equiv 0$: the entire
curvature term vanishes and the Hessian collapses to the clean projected form
$J^{\mathsf T}(\tfrac1{h^2}M + K)J$. You evaluate the real rig's Jacobian *once*
per step instead of $\mathcal O(d^2)$ times per iteration. Crucially they
linearize the **rig**, not the elastic forces — so $W$ stays nonlinear in $p$ and
the integrator keeps its stability.

`ImplicitEuler(..., relinearize=True)` does exactly this: at each step it rebuilds
`LinearizedRig` about the current $p$, solves, advances, and re-linearizes next
step.

## 3. Does it matter? Watch them overlap

The lesson runs both solvers on the same `BendRig` under gravity:

- **2012 full** (solid) — true rig every iteration, curvature included;
- **2013 linearized** (green ghost) — frozen $J$ per step.

The console reports the worst vertex disagreement over 150 steps: about **1% of
the beam's height** — right in line with the paper's reported avg/max error of
0.08% / 1.4%, i.e. *imperceptible*. The two beams track each other through the
overshoot and ring-down. The magenta arrows are the difference field
**magnified ×25** so you can even see where the (tiny) disagreement concentrates
— at the most-bent region, where the rig is most nonlinear.

## What to look at

- **Overlap.** Solid and ghost stay locked together as the beam bends and
  jiggles. Pluck it; they re-excite together.
- **Difference field (×25).** The only visible disagreement is during fast,
  high-curvature motion, and it is still a fraction of the geometry.
- **Curvature readout.** $\lVert\partial^2 s/\partial\theta^2\rVert$ is clearly
  nonzero — this is the term the linearized solver simply never computes.
- Because our `BendRig` has $d=1$, the cost gap ($\mathcal O(d^2)$ vs $\mathcal
  O(1)$) is conceptual here; at $d=174$ (the 2013 sumo) it is the difference
  between minutes-per-frame and under a second.

## The third 2013 trick (briefly)

The paper goes further with **deferred Jacobian evaluation**: even one $J$ per
step is reused across several steps, guarded by a cheap kinetic-energy error
check, and only recomputed when the rig has moved enough to matter. We don't
implement it, but it is the same spirit — exploit temporal coherence to touch the
black-box rig as little as possible.

## Where this is going

You now have the full method: physics in rig space, statics and dynamics by
minimizing an energy, forces moved by $J^{\mathsf T}$, and the linear-rig
shortcut that makes it fast. Lesson 8 puts it together in the production
framing — an animator drives a **primary** parameter while physics adds
**secondary** motion to the rest.

Reference: 2013 §4.1 (linear rig), §4.3 (deferred Jacobian) — see
[`../../../docs/Efficient-Secondary-Motion-Rig-Space.md`](../../../docs/Efficient-Secondary-Motion-Rig-Space.md).
