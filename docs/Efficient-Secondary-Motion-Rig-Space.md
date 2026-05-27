# Efficient Simulation of Secondary Motion in Rig-Space

> Fabian Hahn, Bernhard Thomaszewski, Stelian Coros, Robert W. Sumner, Markus Gross.
> *Symposium on Computer Animation (SCA) 2013.* ETH Zurich / Disney Research Zurich.
> PDF: [`Hahn et al. - 2013 - Efficient simulation of secondary motion in rig-space.pdf`](./Hahn%20et%20al.%20-%202013%20-%20Efficient%20simulation%20of%20secondary%20motion%20in%20rig-space.pdf)

## TL;DR

A 1–2 orders-of-magnitude faster version of [Rig-Space Physics (2012)](./Rig-Space-Physics.md) at the
**same quality**. It keeps the 2012 implicit-Euler energy but (1) **linearizes the rig** once per step
(killing the expensive `∂²s/∂p²` term), (2) **eliminates the interior vertices** via a learned linear
"physics-based volumetric skinning" map `q = W·s`, and (3) **defers Jacobian re-evaluation** under a
kinetic-energy error guard. The simulation then runs purely in the free rig parameters.

## Starting point (recap of 2012, new notation)

Same black-box rig `p ↦ s(p)`. Interior nodes are now `q` (2012 called them `n`). Modified St. Venant–
Kirchhoff per-element energy (eq. 2):

```
W_e(xᵉ) = ½μ‖E‖_F + λ(1 − Vᵉ/V₀ᵉ),     E = ½(FᵀF − I),   F = d·D⁻¹
```

`F` = deformation gradient (`d`, `D` = 3×3 deformed/undeformed tet-edge matrices), `E` = Green strain.
Total `W = Σ_e W_e·V₀ᵉ = W(s, q)`. The 2012 implicit-Euler functional, restated (eq. 3):

```
H(p,q) = (h²/2)(aₛᵀMₛaₛ + a_qᵀM_q a_q) + W(s(p), q)
aₛ = (1/h²)(s(p) − 2sₙ + sₙ₋₁),   a_q = (1/h²)(q − 2qₙ + qₙ₋₁)
```

Two stated bottlenecks: Newton needs **first and second** finite-difference rig derivatives every
iteration, and the system carries `q` (~thousands) as unknowns.

## Contribution 1 — Linear rig approximation (§4.1)

Linearize the rig **once at the start of each timestep**:

```
s(p) ≈ s(pₙ) + J(pₙ)·(p − pₙ),     J = ∂s/∂p          (eq. 4)
```

Key choice: linearize the **rig**, not the elastic forces (unlike Baraff–Witkin semi-implicit). The EOM
stay nonlinear in `p` (better stability); the only nonlinear re-evaluations in Newton are the cheap
elastic gradient/Hessian. Consequences:

- `J` evaluated **once per timestep**, not per Newton iteration.
- The rig's second derivative **`∂²s/∂p²` vanishes entirely** — removing the 2012 `H_pp` curvature terms
  that cost O(p²) rig evals.

## Contribution 2 — Physics-based volumetric skinning (§4.2)

Goal: drop `q` as DOFs. Ideal map (interior in static equilibrium given the surface):

```
q(p) = argmin_q̃ W(s(p), q̃)                              (eq. 5)
```

Implicit/nonlinear → too slow. **Approximate by an explicit linear map `q = W·s`**, learned from
examples. (Cage coords like harmonic/Green are rejected: they need a hand-rigged cage and place interior
vertices far from elastic equilibrium, overestimating `W` and corrupting dynamics.)

### Example generation ("shaking")

Take ~5–10 artist calisthenics poses; apply impulse vectors to the surface as initial velocities; run a
few steps of the *full* dynamic sim (eq. 3). Yields `m` matched (surface, interior) configs
`xᵉ = {sᵉ, qᵉ}` that reflect how simulated params drive the interior.

### Per-vertex weight fit (eq. 6)

For each interior vertex `q_j`, fit surface weights `wʲ` across the `m` examples:

```
            ‖ ⎡ s₁¹ … sₙ¹ ⎤        ⎡ q_j¹ ⎤ ‖²
wʲ = argmin ‖ ⎢  ⋮   ⋱  ⋮ ⎥·w̃ʲ  −  ⎢  ⋮   ⎥ ‖
       w̃ʲ   ‖ ⎢ s₁ᵐ … sₙᵐ ⎥        ⎢ q_jᵐ ⎥ ‖
            ‖ ⎣  1  …  1  ⎦        ⎣  1   ⎦ ‖
```

Each pose-row = 3 equations (x/y/z); the `1…1` row is a **partition-of-unity** constraint (weights sum to
1) for correct behavior under rigid surface motion. Stack into skinning matrix
`W = [w¹,w¹,w¹,…,w^{nq},w^{nq},w^{nq}]ᵀ`, giving `q(s) = W·s` ("physics-based volumetric skinning").

### Sparse correspondences (Algorithm 1)

Raw eq. 6 overfits (negative weights; remote surface vertices). Fix with an **L1 regularizer** on weights
(Schmidt et al. 2007) to force a sparse, local set with bounded residual:

```
for each interior vertex j:
  (wʲ, r₀) = solveNNLS(j, S₀)              # non-neg least squares on full/conservative set
  r = r₀, S = S₀
  while  r < δ  or  r/r₀ < 1.5:            # while fit still acceptable
     (S̃, w̃ʲ) = reduceCorrespondenceSet(S, wʲ)   # drop smallest-weight surface verts
     S = S̃
     (wʲ, r) = solveL1(j, S, w̃ʲ)          # L1-regularized re-solve (Newton)
  (wʲ, r_f) = solveNNLS(j, S)              # final weights WITHOUT L1 (debias)
```

Stops before residual exceeds `δ`. Elephant: from 20 candidates → **~5.5 correspondences/vertex avg**,
tracking ground-truth elastic energy far better than NNLS at 5/20/500 correspondences (Fig. 5). Computed
once in preprocessing.

## Contribution 3 — Deferred Jacobian evaluation (§4.3)

Even one `J`/step dominates cost; temporal coherence means `J` barely changes. Reuse it, guarded by an
indicator: compare linear prediction `s̃ = s(pₙ) + J·Δp` to the true `s(pₙ+Δp)` via the kinetic energy of
the discrepancy (eq. 7):

```
ΔE_kin = (1/2h)·(s̃ − s(pₙ+Δp))ᵀ Mₛ (s̃ − s(pₙ+Δp))      (eq. 7)
```

Costs only **one** rig eval. Optimistic loop: step reusing old `J`, check `ΔE_kin`; if too nonlinear, roll
back, recompute `J(pₙ)`, re-simulate. To avoid excessive rollbacks on rapid motion, add a cheap
**predictive** indicator using an extrapolation of past states:

```
p̃ = pₙ + (1/h)(pₙ − pₙ₋₁)
```

On the elephant **trunk** (nonlinear wire deformer) the predictor flags 280/360 frames for update; without
it, 281 frames became costly rollbacks. On the **belly** (nearly linear rig) deferral pays off fully.

## Resulting gradient & Hessian (§4.4)

With skinning, `H = H(p)` alone. ⚠️ Notation collision: bold `W` = **skinning matrix**; `∂W/∂s` etc. =
**elastic energy**.

```
∂H/∂p = Jᵀ[ Mₛaₛ + ∂W/∂s + Wᵀ( M_q a_q + ∂W/∂q ) ]
```

```
∂²H/∂p² = Jᵀ[ (1/h²)Mₛ + ∂²W/∂s² + Wᵀ( (1/h²)M_q + ∂²W/∂q² )W ] J
        + Jᵀ[ Wᵀ ∂²W/∂q∂s + ∂²W/∂s∂q W ] J
```

Chain rule outward: interior terms `(M_q a_q + ∂W/∂q)` pulled back to the surface by `Wᵀ` (since
`∂q/∂s = W`), combined with surface terms, then projected to rig space by `Jᵀ`. `J` is **constant** within
the step → **no rig second-derivative term** (the 2012 `∂_pJ_s` terms are gone), and there are **no `q`
unknowns**.

## Results (Table 1, Core i7-930, seconds per frame)

`IV1`=static interior solve, `IV2`=skinning; `JE1`=per-frame Jacobian, `JE2`=deferred.

| Solver | Trunk (13 p) | Belly (36 p) | Sumone (174 p) |
|---|---|---|---|
| RSP (2012) | 7.24 ×1 | 46.5 ×1 | — (out of reach) |
| IV1, JE1 | 0.37 ×19 | 0.53 ×86 | 2.82 ×1.0 |
| IV1, JE2 | 0.39 ×18 | 0.39 ×118 | 2.48 ×1.1 |
| IV2, JE1 | 0.12 ×56 | 0.27 ×168 | 1.04 ×2.7 |
| IV2, JE2 | 0.13 ×53 | **0.13 ×348** | **0.58 ×4.9** |

- Quality vs. 2012: **avg/max vertex error < 0.08% / 1.4%** of character height, all frames (imperceptible).
- The 174-parameter sumo (belly, chest, cheeks, hair) is infeasible in 2012; here <1 s/frame.
- Skinning validation (Fig. 5): sparse method tracks reference elastic energy; NNLS at 5/20/500
  correspondences overfits and spikes to ~5× reference.
- Skinning is precompute-only (elephant: 1328 interior + 761 surface, 74 poses, ~30 s/vertex; sumo: 967 +
  1302, 77 poses, ~16.5 s/vertex).

## Limitations

- **No rig-space material control** (the 2012 per-parameter stiffness QP is dropped).
- Materials that are soft at rest *and* stable under fast motion are hard to tune.
- Fixed `h = 0.01 s`; rapid motion needs smaller steps (adaptive stepping = future work).
- Single rest mesh → extreme poses can invert tets (remeshing/meshless = future work).

## Relationship to the original

See [Rig-Space Physics (2012)](./Rig-Space-Physics.md). In short: 2012 minimizes `H` over
`(interior nodes, rig params)` with a nonlinear black-box rig (free interior DOFs + Schur complement +
BFGS to dodge O(p²) rig evals). 2013 minimizes `H(p)` over **rig params only** — interior eliminated by
`q = Ws`, rig replaced per step by `s ≈ s(pₙ) + JΔp`, and `J` recomputed lazily.
