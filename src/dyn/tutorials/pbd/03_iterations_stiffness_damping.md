# Lesson 3 — Iterations, stiffness, and damping

> Run `python src/dyn/tutorials/pbd/03_iterations_stiffness_damping.py`.
> Code: [`03_iterations_stiffness_damping.py`](03_iterations_stiffness_damping.py) · Library: [`../../pbd.py`](../../pbd.py)

Lesson 2 left you with two knobs that act suspiciously like material
properties: **iteration count** and **stiffness $k$**. This lesson explains
what each really does, why naive $k$ is timestep-dependent, and how to add
**damping** without freezing the whole rope solid.

## 1. Iteration count = signal speed

PBD's projection is Gauss-Seidel-style: each constraint sees the previous
ones' updates. But corrections only propagate **one constraint per
iteration**. On a long chain (Lesson 2's 20-particle rope) information
travels at one segment per sweep. So:

- **1 iteration / step**: the rope is mushy. A force at one end takes
  $\sim N$ frames to be "felt" at the other.
- **N iterations / step**: information traverses the whole chain within a
  single timestep — the chain behaves rigidly.

This is the fundamental PBD scaling: stiffness is bounded by your iteration
budget. It is exactly why dense, stiff structures (long ropes, tight cloth)
need many iterations or fall back to compliance-aware variants like XPBD.

## 2. Stiffness $k$ and the timestep trap

The cheap way to make a constraint *softer* than fully-projected is to
multiply each correction by $k\in[0,1]$:

$$ \Delta p \;\leftarrow\; k \cdot \Delta p . $$

Run for $n_s$ iterations and the residual error scales as
$\Delta p\,(1-k)^{n_s}$. So $k$ does **not** linearly control how stiff the
material feels — doubling iterations changes the effective stiffness even at
fixed $k$. Müller's fix:

$$ \boxed{\;k' \;=\; 1 - (1-k)^{1/n_s}\;} $$

Now the residual is $\Delta p\,(1-k')^{n_s} = \Delta p\,(1-k)$ — independent
of $n_s$ at fixed $k$. The user-facing slider behaves like a material
parameter.

(Even with $k'$, the *physical* stiffness still depends on $\Delta t$ —
this is PBD's well-known flaw, fixed by XPBD which gives every constraint a
**compliance** $\alpha = 1/(\text{stiffness} \cdot \Delta t^2)$ that converges
to a true implicit-Euler solution as iterations $\to \infty$.)

## 3. Damping that does not freeze rigid motion

The naive damper, $v \leftarrow (1-c)\,v$, kills *all* motion — including
useful global swing. Müller §3.5 proposes a smart alternative: compute the
system's centre-of-mass linear velocity $v_{cm}$ and angular velocity
$\omega = I^{-1}L$, and damp only the **deviation** from the rigid mode:

$$
\Delta v_i \;=\; v_{cm} + \omega\times r_i - v_i,
\qquad
v_i \;\leftarrow\; v_i + k_{\text{damping}}\,\Delta v_i .
$$

Set $k_{\text{damping}}=1$ and only the rigid translation+rotation survives —
the rope locks into a swinging stick. Set it small (~0.05) and internal
jitter dies quickly while the global swing keeps energy. The script lets you
toggle between this and naive damping; the difference is dramatic on a long
rope.

This is one of those small details that separates PBD demos that *feel*
right from ones that don't.

## What to look at

- **Two ropes side by side** with the same physics but different
  `iterations`. Crank one to 1, the other to 30: the low-iter rope is
  visibly stretchy during fast swings; the high-iter one stays at its rest
  length.
- **stiffness $k$** without the $k'$ correction. Drag the iteration slider
  with $k=0.5$: the rope's perceived stiffness *changes* — it gets stiffer
  with more iterations even though "the constraint" hasn't moved. Toggle
  $k'$ on and the dependency vanishes.
- **damping mode**: switch between rigid-mode-preserving and naive. With
  rigid-mode damping, even at $k_{\text{damping}}=1$ the rope still swings
  as a unit — just rigidly. With naive damping at the same value, *everything*
  including the swing is killed.

## Where this is going

We now have a working 1D PBD demo with the right knobs. Lesson 4 takes the
exact same loop, swaps in a **2D triangle mesh**, and adds **dihedral bend**
constraints — instant cloth.

Reference: Müller et al. 2006 §3.5 (damping) and the $k'$ derivation in §3.3
— see [`../../../../docs/Position-Based-Dynamics.md`](../../../../docs/Position-Based-Dynamics.md).
