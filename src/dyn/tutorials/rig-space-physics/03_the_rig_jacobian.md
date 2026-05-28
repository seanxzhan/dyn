# Lesson 3 — The rig Jacobian and the force projection $J^{\mathsf T}$

> Run `python src/dyn/tutorials/03_the_rig_jacobian.py`.
> Code: [`03_the_rig_jacobian.py`](03_the_rig_jacobian.py) · Library: [`../rig.py`](../rig.py)

Lesson 2 gave us the rig map $s(p)$. To do physics in rig space we need its
**derivative**. This single object — the rig Jacobian — is the mathematical heart
of both papers.

## 1. The Jacobian: a basis for allowed motion

The rig Jacobian is

$$ J(p) \;=\; \frac{\partial s}{\partial p} \;\in\; \mathbb{R}^{3N \times d}. $$

Column $i$ is $\partial s/\partial p_i$, a vector in $\mathbb{R}^{3N}$ — i.e. a
**displacement field over all vertices**: "if I nudge $p_i$ by a hair, this is
the direction (and relative speed) every vertex moves." The lesson draws each
column as a per-vertex arrow field.

The columns of $J$ span the **tangent space** to rig space $\mathcal{M}$ at $p$:
the only velocities the rig can produce are $\dot s = J\dot p$. Forces that try to
push the mesh in any direction *outside* $\mathrm{col}(J)$ simply cannot be
followed — the rig has no coordinate for them.

A first-order (Taylor) view, used everywhere later:

$$ s(p + \Delta p) \;\approx\; s(p) + J(p)\,\Delta p. $$

## 2. Linear vs. nonlinear: is $J$ constant?

**LinearRig.** $s(p)=s_0+Bp$, so

$$ J = B \quad\text{(constant)}, \qquad \frac{\partial^2 s}{\partial p^2} = 0. $$

Move the sliders: the arrow field never changes. Rig space is flat, so its
tangent basis is the same everywhere.

**BendRig.** One angle $\theta$ rotates each vertex about the clamp by its own
$\alpha_v=\theta\,w_v$ (a smoothstep weight $w$ from clamp to tip — literally
one-bone linear-blend skinning). With $R_y$ the rotation about $y$,

$$ s_v(\theta) = \text{pivot} + R_y(\alpha_v)\,\text{offset}_v, $$
$$ \frac{\partial s_v}{\partial \theta} = w_v\,R_y'(\alpha_v)\,\text{offset}_v,
\qquad
\frac{\partial^2 s_v}{\partial \theta^2} = w_v^2\,R_y''(\alpha_v)\,\text{offset}_v. $$

Now $J(\theta)$ **depends on $\theta$** (it contains $\sin\alpha,\cos\alpha$), and
the second derivative is **nonzero** — the curvature of rig space. Bend the beam
in the viewer and watch the arrows rotate with it.

This is the whole 2012-vs-2013 story in one slider:

- The curvature $\partial^2 s/\partial p^2$ ("$\partial_p J$") is what the **2012**
  Newton solve needs, and evaluating it for a black-box rig costs
  $\mathcal{O}(d^2)$ rig calls *per iteration* — the bottleneck.
- The **2013** method just *linearizes the rig* once per step
  ($s\approx s(p_n)+J\Delta p$), making this term vanish. Constant-$J$ is the
  `LinearRig` regime — you are looking at why the fast method chose it.

## 3. Verifying $J$ — finite differences

For a real (black-box) rig you rarely have a formula for $J$. You get it by
**finite differences**, exactly as the papers do:

$$ J_{:,i} \approx \frac{s(p+\varepsilon e_i) - s(p-\varepsilon e_i)}{2\varepsilon}. $$

`finite_difference_jacobian` does this; the script asserts it matches our
closed-form $J$ to ~$10^{-10}$ for both rigs. (For a black box, $d$ parameters
cost $\mathcal{O}(d)$ rig evaluations for $J$, and $\mathcal{O}(d^2)$ for the
curvature — hence the cost story above.)

## 4. The punchline: $J^{\mathsf T}$ projects forces into rig space

How does a force on the *mesh* become a force on the *parameters*? Through the
**transpose** of the Jacobian. If $f_s$ is a force in vertex space, the work it
does under a rig motion $\delta p$ is

$$ \delta W = f_s^{\mathsf T}\,\delta s = f_s^{\mathsf T} J\,\delta p
           = \underbrace{(J^{\mathsf T} f_s)}_{f_p}{}^{\mathsf T}\,\delta p . $$

So the **generalized force** conjugate to $p$ is

$$ \boxed{\,f_p = J^{\mathsf T} f_s\,}. $$

This is the operator that appears in both papers' equilibrium condition
$\partial_p H = J^{\mathsf T}\!\big[M_s a_s + \partial_s W\big] = 0$: take the
ordinary surface forces and **project them with $J^{\mathsf T}$** into the few
rig coordinates. Likewise the rig-space stiffness is $J^{\mathsf T} K\, J$. Once
you have $J$, *all* of the physics gets pulled into rig space by sandwiching with
$J$ and $J^{\mathsf T}$.

Enable **force projection** in the viewer: a downward point force $f_s$ at the
tip (red) yields the rig-space force $f_p=J^{\mathsf T}f_s$ (printed), and the
blue field $J f_p = J J^{\mathsf T} f_s$ shows the motion the rig actually takes
— the input force *projected onto what the rig can do*. With the LinearRig the
`bend_z` component dominates $f_p$; with the BendRig, $f_p$ changes as you bend,
because $J$ itself rotates.

## What to look at

- **LinearRig** + Jacobian columns: two fixed fields (`bend_z` lifts in $z$,
  `sway_y` pushes in $y$); sliders don't change the arrows.
- **BendRig** + Jacobian columns: bend the beam and the single arrow field turns
  — you are seeing rig curvature directly.
- **Force projection**: note the blue response never points exactly along $f_s$;
  it is $f_s$ filtered through $\mathrm{col}(J)$.

## Where this is going

We can now move forces between vertex space and rig space. Lesson 4 supplies the
forces themselves — the elastic energy $W$, its gradient (forces), and Hessian
(stiffness) — and then Lessons 5–6 assemble statics and dynamics in rig space by
wrapping them in $J^{\mathsf T}(\cdot)J$.

Reference: Hahn et al. 2012 §3 and 2013 §4 — see
[`../../../docs/Rig-Space-Physics.md`](../../../docs/Rig-Space-Physics.md) and
[`../../../docs/Efficient-Secondary-Motion-Rig-Space.md`](../../../docs/Efficient-Secondary-Motion-Rig-Space.md).
