# Lesson 4 — Cloth: stretch + dihedral bend

> Run `python src/dyn/tutorials/pbd/04_cloth_stretch_and_bend.py`.
> Code: [`04_cloth_stretch_and_bend.py`](04_cloth_stretch_and_bend.py) · Library: [`../../pbd.py`](../../pbd.py)

The hanging-rope demo of Lesson 2, generalised to a 2D triangle mesh, is
already cloth. The only new thing this lesson introduces is a **second
constraint type** for bending — without it the cloth is "limp paper towel";
with it, "canvas".

## 1. Cloth = particles + stretch edges + bend pairs

`make_cloth(res=(M,N))` in [`../../pbd.py`](../../pbd.py) generates:

- $M{\times}N$ particles arranged on a regular grid in the cloth plane;
- a triangulation: each grid quad split along one diagonal into two
  triangles (consistent across the whole mesh, so the diagonals point the
  same way);
- one **stretch edge** per unique triangle edge, each with rest length =
  initial length;
- one **bend pair** per *interior* edge — i.e. each edge shared by two
  triangles — with the four indices $(p_1, p_2, p_3, p_4)$ where $(p_1,p_2)$
  is the shared edge and $p_3, p_4$ are the off-edge corners of the two
  adjacent triangles. Rest dihedral $\varphi_0 = 0$ for a flat cloth.

Pin the top row by setting $w_i = 0$ there. Hang under gravity. Done.

## 2. The dihedral-bend constraint

For a pair of triangles sharing edge $(p_1, p_2)$ with off-edge vertices
$p_3, p_4$:

$$
n_1 \;=\; \widehat{(p_3-p_1)\times(p_2-p_1)},\qquad
n_2 \;=\; \widehat{(p_2-p_1)\times(p_4-p_1)},
$$

$$
C_{\text{bend}}(p_1,p_2,p_3,p_4) \;=\; \arccos(n_1\cdot n_2) - \varphi_0.
$$

Two key properties:

1. **Length-independent.** $C$ depends only on the *angle* between triangle
   normals, not on edge lengths. So bending stiffness can be set
   independently of stretching stiffness — exactly the paper's headline
   advantage over distance-based bending (Provot 1995) or
   discrete-shells-style bend (Grinspun et al. 2003) which couple the two.
2. **Rigid-motion invariant.** $C$ is unchanged under translation +
   rotation, so $\nabla C$ is perpendicular to rigid-body modes ⇒ the
   projection conserves momentum.

The gradients have a clean geometric form: rotating triangle 1 around the
shared edge by $\Delta\theta$ moves $p_3$ by $h_a\,\Delta\theta$
perpendicular to the triangle, in direction $\pm n_1$, where $h_a$ is the
perpendicular distance from $p_3$ to the edge. So

$$
\nabla_{p_3} C = \frac{n_1}{h_a},\qquad
\nabla_{p_4} C = \frac{n_2}{h_b},
$$

and the gradients on $p_1, p_2$ follow from translation invariance — the
four gradients must sum to zero, with the split between $p_1, p_2$ given by
the along-edge position of each perpendicular foot. The full formula is in
`project_bend` in [`../../pbd.py`](../../pbd.py); the script asserts it
matches finite differences to ~$10^{-10}$ via `verify_bend_gradient`.

(Aside: the Müller paper's Appendix A derives essentially the same gradient
via the chain rule on the normalized cross product; the geometric form above
is equivalent and a bit cleaner to vectorize.)

## 3. The PBD step is unchanged

Two constraint types, same loop:

```
predict
for it in range(n_iters):
    project_distance(p, edges, rest, w, k_prime=k_stretch')
    project_bend(p, bend_pairs, rest_dihedrals, w, k_prime=k_bend')
recover velocity, commit
```

Order matters slightly (Gauss-Seidel) — projecting stretch first and bend
second tends to look better because stretch is "stiffer" geometrically and
should win, but either order is valid.

## What to look at

- **bend stiffness slider.** At $k_{\text{bend}}=0$ the cloth is limp — it
  collapses into a thin column under gravity, with sharp folds. Crank
  upwards: the cloth resists folding and starts behaving like canvas, then
  cardstock. **Stretch stiffness is held at 1** the whole time so you can
  see the two are decoupled.
- **stretch stiffness slider.** Drop it well below 1 and the cloth
  *elongates* under its own weight while still resisting folds — the
  decoupling lets you fake stretchy fabric.
- **iterations** above ~10 give a reasonable feel for a $20{\times}15$ grid.
  More iterations = stiffer for both stretch and bend (the $k'$ correction
  decouples user-facing stiffness from iteration count).
- **wind kick.** A short impulse pushes the cloth horizontally, then it
  swings back. The bend stiffness controls how "papery" the wave looks.

## Where this is going

The cloth has internal physics. Lesson 5 adds **collisions**: a static
sphere collider with friction and restitution, which is just one more
inequality constraint inside the same projection loop. From there the
capstone (Lesson 6) makes the sphere dynamic and throws it.

Reference: Müller et al. 2006 §4.1 (cloth representation) and Appendix A
(bend gradient) — see
[`../../../../docs/Position-Based-Dynamics.md`](../../../../docs/Position-Based-Dynamics.md).
