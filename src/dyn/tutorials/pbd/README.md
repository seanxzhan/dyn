# Position Based Dynamics — a step-by-step, visual course

A from-scratch tour of Müller et al. 2006 (`docs/Position-Based-Dynamics.md`),
built on a procedural cloth so every formula is checkable in closed form, and
visualized interactively with [Polyscope](https://polyscope.run). Each lesson
is a runnable script paired with a math note; shared code lives in
[`../../pbd.py`](../../pbd.py) (the library) and [`../viz.py`](../viz.py).

The arc ends at the **capstone** (Lesson 6): a ball thrown at a hanging cloth,
two-way coupled, every piece from 1–5 in the same loop.

## Setup & running

```bash
source activate.sh          # conda activate dyn   (numpy + polyscope already installed)
python src/dyn/tutorials/pbd/01_particles_and_predict.py
```

Each script opens a Polyscope window with sliders/checkboxes; read the matching
`.md` alongside it. Close the window to exit.

- `DYN_SMOKE=1` — run headless (mock GL backend, a few frame ticks) to verify a
  lesson without opening a window. Used for the checks below.
- `DYN_NO_SHOW=1` — run only the NumPy part (math + asserts), skip Polyscope.

## The arc

Each lesson adds exactly one idea. The thread mirrors the paper: start with
particles + Euler, layer in distance constraints, then iteration tuning, then
the second constraint type (bend), then collisions, then dynamic-sphere
coupling.

| # | Lesson | One idea | Status |
|---|--------|----------|--------|
| 01 | [Particles + the predict step](01_particles_and_predict.md) | $x, v, w$; $v\!\leftarrow\! v + h g$; $p\!\leftarrow\! x + h v$ — no projection yet | ✅ built |
| 02 | [The distance constraint](02_distance_constraint.md) | $C(p_a,p_b)=\|p_a-p_b\|-\ell_0$, $\Delta p_i = -s w_i \nabla C$, hung rope | ✅ built |
| 03 | [Iterations, stiffness $k'$, damping](03_iterations_stiffness_damping.md) | iters = signal speed; $k'=1-(1-k)^{1/n_s}$; rigid-mode-preserving damping | ✅ built |
| 04 | [Cloth: stretch + dihedral bend](04_cloth_stretch_and_bend.md) | a 2D triangle mesh, dihedral angle as constraint, length-independent bend | ✅ built |
| 05 | [Collisions and response](05_collisions_and_response.md) | inequality constraints, friction $\mu$ and restitution $e$ as velocity post-process | ✅ built |
| 06 | [**Capstone**: ball through cloth](06_capstone_ball_through_cloth.md) | dynamic sphere, two-way impulse coupling, sub-stepping for tunneling | ✅ built |

The library [`../../pbd.py`](../../pbd.py) ships:

- a `Cloth` factory (grid + stretch edges + bend pairs in one call)
- `predict`, `project_distance`, `project_bend`, `project_sphere` — the four
  building blocks of the loop
- `damp_velocities` — Müller §3.5's rigid-mode-preserving damper
- `DynamicSphere` + `project_dynamic_sphere` — the capstone's two-way coupler
- `verify_bend_gradient` — finite-difference check the analytic bend gradient
  matches FD to ~$10^{-11}$, asserted by the cloth lesson

## Verify all lessons

```bash
# math only (fast, no GUI): builds, asserts FD bend gradient, asserts no collision penetration
for f in src/dyn/tutorials/pbd/0*.py; do DYN_NO_SHOW=1 python "$f"; done
# full GUI code path, headless, for every lesson:
for f in src/dyn/tutorials/pbd/0*.py; do DYN_SMOKE=1 python "$f"; done
```

## Why a procedural cloth and not a real garment?

The cloth grid factory is closed-form: vertex positions are $X u + Z v$ on a
chosen plane, triangulation is a regular split, rest lengths are exact, rest
dihedrals are 0. This means every constraint can be **verified analytically**
— `verify_bend_gradient` checks the dihedral gradient against finite
differences to machine precision, the distance constraints reduce to a
solvable hanging-rope geometry, and the inequality of Lesson 5 is a sphere
test that evaluates in closed form. A real garment buys none of this and
costs an FBX parser + UV-aware mesher.

The same cloth + ball setup, scaled up to a 100x100 grid and 10 iterations,
is what production cloth was in 2006. The capstone runs interactively at the
20-30 grid scale on a laptop CPU — about 4 orders of magnitude removed from a
movie-quality drape, but every piece of physics is the same.

## What this course does not cover

- **Self-collision.** Cloth-vs-itself, broadphase via spatial hash,
  continuous collision detection. Listed as future work in the paper; the
  field has since converged on signed-distance / repulsion + impulse-based
  responses.
- **XPBD.** The compliance trick that decouples stiffness from $\Delta t$
  (Macklin et al. 2016). It's a one-line drop-in over the loop here; mostly
  you replace $k'$ with $\alpha = 1/(k h^2)$ and add a Lagrange multiplier
  per constraint.
- **Strain limiting.** Hard cap on per-edge stretch. Mentioned in the
  capstone notes as the next thing you'd add for hero-quality cloth.

The throughline of this course is *positions over forces*, *projection over
differentiation*, and *constraint as the lingua franca*. Once you have it,
every paper since 2006 is an extension or refinement of this loop.
