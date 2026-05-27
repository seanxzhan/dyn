# Rig-Space Physics

> Fabian Hahn, Sebastian Martin, Bernhard Thomaszewski, Robert Sumner, Stelian Coros, Markus Gross.
> *ACM Transactions on Graphics 31(4), SIGGRAPH 2012.* ETH Zurich / Disney Research Zurich.
> PDF: [`Hahn et al. - Rig-Space Physics.pdf`](./Hahn%20et%20al.%20-%20Rig-Space%20Physics.pdf)

## TL;DR

Run physics-based simulation **inside the deformation subspace defined by an animator's rig** instead
of in the full space of mesh vertices. Implicit Euler is cast as an energy minimization, minimized over
the interior FEM nodes plus the rig parameters, with the rig Jacobian projecting nodal forces/stiffness
into rig space. The output is ordinary **animation curves on the rig parameters** — editable like
keyframes — letting artists dial continuously between hand animation and free simulation.

## The problem

- A **rig** maps a small set of artist controls `p` to a high-resolution mesh `s`.
  Expressive and editable, but tedious for physical effects (jiggle, sag, follow-through).
- **Full physics** gives those effects for free, but produces raw per-vertex motion that does not live
  on the rig — animators can't edit it and it doesn't fit the pipeline.
- Rig-space physics bridges the two: physical motion expressed *as rig-parameter curves*.

## Method

### 1. Dynamics as a minimization

Continuum EOM `ρẍ = f(x)` with `f = −∂ₓW` (potential `W = W_int + W_ext`).
Discretize in time with implicit (backward) Euler:

```
ρ( (x_{n+1} − xₙ)/h² − vₙ/h ) = f(x_{n+1})
```

Because `f = −∂ₓW`, this is the stationarity condition of an **incremental potential** (eq. 2):

```
H[x_{n+1}] = (ρh²/2)·‖ (x_{n+1} − xₙ)/h² − vₙ/h ‖²  +  W(x_{n+1})
```

- First term: inertial penalty pulling `x` toward the free-flight prediction `xₙ + h·vₙ`.
- Second term: elastic + external (gravity, penalty collisions) energy.

Each timestep is therefore a minimization — which makes the rig subspace trivial to add: just minimize the
same `H` over fewer variables.

### 2. Spatial discretization with a split DOF set

Tetrahedralize the character. DOFs split into:

```
x = {n} ∪ {s}
```

- `s` = surface nodes, **tied to the rig**: `s = s(p)`. The map `p ↦ s(p)` is a **black box**
  (skeletons, blendshapes, cages, FFD — anything); derivatives `∂s/∂p`, `∂²s/∂p²` via finite differences
  if not supplied.
- `n` = interior nodes, kept as **free, independent** FEM DOFs (the rig only deforms the surface).

Linear tet basis, StVK or Neo-Hookean `W_int`, lumped masses `M_n, M_s`. The incremental potential (eq. 3):

```
H(n,p) = (h²/2)·[(n−nₙ)/h² − vₙ/h]ᵀ M_n [(n−nₙ)/h² − vₙ/h]
       + (h²/2)·[(s(p)−sₙ)/h² − wₙ/h]ᵀ M_s [(s(p)−sₙ)/h² − wₙ/h]
       + W(n, s(p))
```

Each timestep minimizes `H` over `(n, p)`.

### 3. First-order conditions — the rig Jacobian

`[∂_nH, ∂_pH]ᵀ = 0`. Interior equation = ordinary FEM EOM:

```
∂_nH = M_n[(n−nₙ)/h² − vₙ/h] + ∂_nW = 0
```

Rig equation = surface EOM **projected into rig space by `J_s = ∂s/∂p`**:

```
∂_pH = J_sᵀ M_s[(s(p)−sₙ)/h² − wₙ/h] + J_sᵀ ∂_sW = 0
```

The `J_sᵀ(·)` projection (generalized forces in reduced coordinates) is the mathematical heart of the paper.

### 4. Newton system and Hessian

Solve `[Δn, Δp]` from (eq. 4) with blocks (eq. 5):

```
H_nn = (1/h²)M_n + ∂_nnW                              (sparse: FEM tangent + mass)
H_pn = J_sᵀ ∂_snW                                      (interior↔rig coupling)
H_pp = J_sᵀ[ (1/h²)M_s + ∂_ssW ]J_s                    ← projected surface tangent + mass
     + ∂_pJ_sᵀ · M_s[(s(p)−sₙ)/h² − wₙ/h]               ← rig-curvature × inertia residual
     + ∂_pJ_sᵀ · ∂_sW                                   ← rig-curvature × surface force
```

The last two `H_pp` terms involve `∂_pJ_s = ∂²s/∂p²` (rig curvature). With a black-box rig these cost
**O(p²) rig evaluations per Newton iteration** — the bottleneck. Line search uses cubic interpolation
satisfying the Wolfe conditions.

### 5. Two tricks for tractability

- **BFGS on the `H_pp` block only.** Avoids the O(p²) curvature evals via a rank-2 quasi-Newton update;
  `H_nn` (large, sparse, cheap) stays analytic. Trades convergence rate for fewer black-box rig calls.
- **Schur complement (block Gauss elimination).** `H_nn` is sparse but `H_pn, H_pp` are dense (the `J_s`
  projection). Condense the interior out:

  ```
  (H_pp − H_pn H_nn⁻¹ H_np) Δp = ∂_pH − H_pn H_nn⁻¹ ∂_nH
  H_nn Δn = ∂_nH − H_np Δp
  ```

  Prefactor `H_nn` once; the reduced dense system in `Δp` is tiny (few rig params).

## Extensions

### Rig-space material control (§5)

Artists set one stiffness scalar `S_i` per rig parameter; the system infers a per-element stiffness scale
`µ_e`. Uses **static condensation** `q(p) = (n(p), s(p))` (interior in static equilibrium given the
boundary), `J_q = ∂q/∂p`. At rest the rig-space stiffness is

```
H_pp = J_qᵀ ∂_qqW_int J_q = Σ_e ( J_qᵀ ∂_qqW_e J_q )
```

Demand it match the desired scaling `S` (eq. 6):

```
S·H_pp = Σ_e ( J_qᵀ ∂_qqW_e J_q )·µ_e
```

`p²` equations in `#elements` unknowns; solved as a constrained least-squares QP with `µ_e ≥ 0` and an
`L²` regularizer `Σ_e(µ_e−1)²` (prefer default stiffness). Handles overlapping rig-parameter influence
regions correctly (Fig. 3).

### Physics-based high-level rig parameters (§6)

Collect simulated rig states `P = [p₁ … p_m]`, take right singular vectors `U`, substitute `p = U·r`
(span, **not** mean-centered, so rest `p=0` stays representable). Lower-dimensional `r` → big speedup
(30 → 6 DOFs, near-imperceptible difference). "Kernel-PCA in spirit": linear analysis in the already
nonlinearly-warped rig space.

### Physics-based inverse kinematics (§6)

Physical energy resolves redundancy:

```
min_{p,n}  H(n,s(p)) + α·Σ_{k∈H} ‖s_k(p) − h_k‖²
```

Works static (drop momentum) or dynamic.

### Interior secondary motion toggle (§6)

`n` and `s` are decoupled, so drop `M_n` (interior → static equilibrium each step, stiffer objects) or
keep it (waves propagate, softer objects). Same switch on `M_s` toggles dynamic vs. static rig params.

## Results (Table 1, Core i7-930, seconds per frame — not real time)

| Model | dim p | dim n | dim s | t/frame (s) |
|---|---|---|---|---|
| elephant belly | 36 | 3990 | 2253 | 46.47 |
| elephant belly (PCA) | 6 | 3990 | 2253 | 3.77 |
| elephant trunk | 13 | 3990 | 2253 | 7.24 |
| flower | 24 | 1089 | 987 | 10.54 |
| sphere scaling | 3 | 5430 | 2580 | 1.41 |

- Newton converged in <5 iterations; BFGS in 10–20.
- Bottleneck = black-box Maya interface + O(p²) finite-difference rig evals. Analytic Jacobians >3×
  faster; in-process rig vs. round-tripping through Maya ~18× faster.

## Limitations

- Not interactive for many rig parameters (Maya data transfer + O(p²) rig evals).
- Fast input motion → instability (boundary conditions jump between steps); they slowed input, simulated,
  then resampled.
- Material-stiffness analysis is **local at rest**; nonlinear rigs may not behave intuitively elsewhere.

## Relationship to the follow-up

The 2013 paper [Efficient Simulation of Secondary Motion in Rig-Space](./Efficient-Secondary-Motion-Rig-Space.md)
removes the interior DOFs (learned linear skinning), linearizes the rig per step (kills `∂²s/∂p²` and the
O(p²) evals), and defers Jacobian updates — 1–2 orders of magnitude faster at the same quality.
