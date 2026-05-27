# Lesson 4 — Elastic energy: the potential $W$, its forces and stiffness

> Run `python src/dyn/tutorials/04_elastic_energy.py`.
> Code: [`04_elastic_energy.py`](04_elastic_energy.py) · Library: [`../energy.py`](../energy.py)

Lessons 1–3 were kinematics — *where* things are. Now we add the physics: an
energy that says how much a deformation *costs*. Its gradient is force, its
Hessian is stiffness, and the rest of the course is minimizing it.

## 1. Measuring deformation: the deformation gradient

Inside one tet the map from rest material coordinates to deformed world
coordinates is affine, so its gradient is a single constant matrix per tet, the
**deformation gradient** $F \in \mathbb{R}^{3\times3}$. With edge matrices
$D_m=[\,V_1{-}V_0,\,V_2{-}V_0,\,V_3{-}V_0\,]$ (rest) and $D_s$ (deformed),

$$ F = D_s\,D_m^{-1}. $$

At rest $D_s=D_m$ so $F=I$. $F$ encodes everything: stretch, shear, rotation.

## 2. Strain must ignore rotation

A material doesn't store energy when you merely *rotate* it. But $F$ changes
under rotation ($F\to RF$), so energy can't be a function of $F$ directly. The
fix is the **Green (Lagrangian) strain**

$$ E = \tfrac12\big(F^{\mathsf T}F - I\big). $$

Because $F^{\mathsf T}F$ kills any rotation ($(RF)^{\mathsf T}(RF)=F^{\mathsf
T}F$), $E=0$ for *every* rigid motion. The script verifies this: rotating the
whole beam $35°$ stores ~$0$ energy. (A linear strain
$\tfrac12(F+F^{\mathsf T})-I$ would wrongly charge for rotation — that is why we
need the nonlinear $E$.)

## 3. The St. Venant–Kirchhoff energy

The simplest rotation-invariant material is **StVK**: energy density quadratic
in $E$,

$$ \Psi(F) = \mu\,\lVert E\rVert_F^2 + \tfrac{\lambda}{2}\,\mathrm{tr}(E)^2,
\qquad W = \sum_e \Psi(F_e)\,V^0_e. $$

$\mu$ (shear) and $\lambda$ (volume) are the **Lamé parameters**, derived from a
Young's modulus $E_{\text{young}}$ and Poisson ratio $\nu$
(`lame_params`). The two terms penalize *shape* change and *volume* change; this
is the continuum version of the springs in a jiggle bone, and it is the $W$ that
appears in both papers' incremental potential.

## 4. Force = $-\nabla W$, stiffness = $\nabla^2 W$

Differentiating $W$ gives the elastic **force** $f=-\partial W/\partial x$. The
clean route is the **first Piola–Kirchhoff stress**

$$ P(F) = F\,\big(2\mu E + \lambda\,\mathrm{tr}(E)\,I\big), $$

after which the gradient on a tet's nodes $1,2,3$ is the matrix
$V^0\,P\,D_m^{-\mathsf T}$ (its columns), and node $0$ gets minus their sum. The
script confirms this analytic gradient against a finite-difference directional
derivative of $W$ (they agree to ~5 digits). The **stiffness** $K=\partial^2
W/\partial x^2$ — how forces change as nodes move — is what Newton will need;
we never assemble the full $K$, because in rig space we only ever need the tiny
$J^{\mathsf T}KJ$ (Lessons 5–6).

## What to look at

- **Energy-density heatmap.** Bend the beam (`p[0]`) and the per-tet $\Psi$ (red
  cells) lights up **near the clamp**, where the bend curvature — and hence the
  strain — is largest. The tip, which moves the most, is barely strained because
  it moves almost rigidly. Strain is about *differential* motion, not
  displacement.
- **Forces.** Toggle the blue arrows: $-\partial W/\partial x$ points everywhere
  back toward the rest shape — the restoring force the physics will balance
  against gravity and inertia.
- **Total $W$.** Watch the readout grow as you push the sliders away from rest;
  it is exactly $0$ at $p=0$.

## Where this is going

We now have $W$, its gradient, and (implicitly) its Hessian. Two ingredients
remain trivial wrappers around them:

- **Statics (Lesson 5):** minimize $W(s(p))$ + gravity over $p$. Forces project
  to rig space as $J^{\mathsf T}\partial_x W$; stiffness as $J^{\mathsf T}KJ$.
- **Dynamics (Lesson 6):** add the inertia term and minimize the incremental
  potential each timestep.

Reference: 2013 eq. 2 (the modified StVK they use) — see
[`../../../docs/Efficient-Secondary-Motion-Rig-Space.md`](../../../docs/Efficient-Secondary-Motion-Rig-Space.md).
