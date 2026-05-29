# Lesson 2 — The distance constraint, and the projection formula

> Run `python src/dyn/tutorials/pbd/02_distance_constraint.py`.
> Code: [`02_distance_constraint.py`](02_distance_constraint.py) · Library: [`../../pbd.py`](../../pbd.py)

This is the only piece of math you have to learn. Once it lands, every other
PBD constraint type is the same recipe with a different $C$.

## 1. A constraint as a scalar function

A constraint is a scalar function $C(p_1,\dots,p_n)=0$ over the positions of
the few particles it touches. The simplest one is the **distance constraint**:

$$ C_{\text{dist}}(p_1,p_2) \;=\; \lVert p_1 - p_2 \rVert - d , $$

zero when the two particles are exactly $d$ apart. PBD's only job each
iteration is to nudge $p_1, p_2$ until $C = 0$.

## 2. The projection formula

We want $\Delta p$ such that $C(p+\Delta p)=0$. Linearize:

$$ C(p+\Delta p) \;\approx\; C(p) + \nabla C(p) \cdot \Delta p \;=\; 0. $$

There are infinitely many $\Delta p$ that satisfy this — pick the one that
**preserves linear and angular momentum** of the system. For an internal
constraint $C$ is invariant under rigid motion, so $\nabla C$ is automatically
**perpendicular to the rigid-body modes** (translation and rotation). If we
restrict the correction to lie along $\nabla C$, momentum is conserved for
free. So write

$$ \Delta p \;=\; \lambda\,\nabla C(p), $$

substitute, solve for $\lambda$, weight each particle by its inverse mass
$w_i$ (so a heavy particle moves less than a light one for the same impulse,
and a pinned $w_i = 0$ doesn't move at all):

$$
\boxed{\;\Delta p_i \;=\; -\,s\,w_i\,\nabla_{p_i} C(p), \qquad
       s \;=\; \frac{C(p)}{\sum_j w_j\,\lVert\nabla_{p_j} C(p)\rVert^2}\;}
$$

That formula, with the right $C$ and $\nabla C$, is the projection step for
*every* constraint in the paper. Memorise it. Lessons 4–5 just plug in
different $C$'s.

## 3. Worked example — the distance constraint

For $C(p_1,p_2)=\lVert p_1-p_2\rVert - d$ with unit normal
$n=(p_1-p_2)/\lVert p_1-p_2\rVert$,

$$ \nabla_{p_1} C = n, \qquad \nabla_{p_2} C = -n,
\qquad \lVert\nabla_{p_1}C\rVert^2 = \lVert\nabla_{p_2}C\rVert^2 = 1. $$

Plug into the boxed formula:

$$
\Delta p_1 = -\frac{w_1}{w_1+w_2}\,(\lVert p_1-p_2\rVert - d)\,n,
\qquad
\Delta p_2 = +\frac{w_2}{w_1+w_2}\,(\lVert p_1-p_2\rVert - d)\,n .
$$

Two readings of the same equations:
- **Symmetric**: both endpoints share the displacement, weighted by mass.
- **Pinning is automatic**: set $w_1 = 0$ and $\Delta p_1 = 0$ and $\Delta p_2$
  alone closes the gap (because $w_1+w_2 = w_2$). No special branch in code.

This is the formula Jakobsen [Jak01] used for his Verlet rope. The PBD paper
recovers it as a special case; everything else follows the same recipe.

## 4. The full one-step PBD loop, now alive

Lesson 1's gutted loop now has the projection step we glossed over:

```
v ← v + h · w · g                            # external forces
p ← x + h · v                                 # PREDICT
for it in range(n_iters):
    project_distance(p, edges, rest, w)       # ← Lesson 2 ships this
v ← (p − x) / h                               # recover
x ← p
```

`project_distance` in [`../../pbd.py`](../../pbd.py) is one Jacobi sweep over
all edges using the formula above. Lesson 3 will discuss what changes when we
run it Gauss-Seidel style or repeat it many times.

## 5. The demo: a hanging rope

We chain $N$ particles, pin the top one ($w_0 = 0$), connect each consecutive
pair with a distance constraint, and let gravity do the rest. Released from
horizontal, the rope swings, each segment fighting to stay at its rest
length. Without constraints (Lesson 1) the particles would just fall
parallel; *with* the projection step they swing as a coupled chain.

## What to look at

- Hit **Play.** Released horizontal under gravity, the rope swings down,
  loses some energy to numerical damping, and settles toward straight-down.
  Watch the trail: chaotic but bounded.
- **iterations slider.** With 1 iteration (Jacobi), corrections only
  propagate one edge per step, so the rope feels stretchy — you can see it
  elongate during fast swings. Crank to 20+ and the chain becomes nearly
  inextensible. (This is exactly the iteration → stiffness correspondence we
  unpack in Lesson 3.)
- **pin top off / on.** Off: $w_0 \ne 0$ and the whole rope drops as a
  group, internally cohering as a falling chain. On: classical pendulum.
- Toggle **projection** off: the chain immediately disintegrates as in
  Lesson 1 — the only thing holding it together *is* the projection step.

## Where this is going

We have one constraint type and a working pendulum. Lesson 3 looks at the
iteration / stiffness story (and the linearized stiffness fix $k'$ that makes
the user-facing knob behave). Lesson 4 swaps in two new $C$ functions
(stretch and dihedral bend) and gets cloth.

Reference: Müller et al. 2006 §3.3 — see
[`../../../../docs/Position-Based-Dynamics.md`](../../../../docs/Position-Based-Dynamics.md).
