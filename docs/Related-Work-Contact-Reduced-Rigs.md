# Related Work — Contact in Reduced & Rig-Space Models

> Literature landscape for **contact-reactive secondary motion on a rig / reduced subspace**, gathered while
> scoping the [research direction](./Research-Direction.md). Reusable as a related-work section. Organized by
> how each work relates to the target idea (contact-reactive secondary rig motion with no volumetric solve).

## Verdict

**The mechanism exists in pieces, and the *bare* combination ("secondary rig motion + contact + no FEM")
already ships in game engines.** The clean, principled version — contact mapped onto the *existing* rig with
an explicit answer to where the rig can't represent the contact, all without a volumetric elastic solve — is
not occupied, but it is hemmed in tightly on three sides. Novelty cannot be "spring bones that collide"; it
must be the **surface→bone contact mapping** and the **representability-gap handling**.

## 1. Rig-space / subspace physics (the lineage)

- **Hahn et al. 2012, Rig-Space Physics** — physics inside the animator's rig subspace; contact enters only
  as a penalty in `W_ext`, projected to rig space by `Jᵀ`. The `Jᵀ`-projection of contact force is the
  mechanism "contact drives the bones." Digest: [Rig-Space-Physics.md](./Rig-Space-Physics.md).
- **Hahn et al. 2013, Efficient Simulation of Secondary Motion in Rig-Space** — 1–2 orders faster; linearized
  rig + learned interior skinning. Digest: [Efficient-Secondary-Motion-Rig-Space.md](./Efficient-Secondary-Motion-Rig-Space.md).

## 2. Spring/handle secondary motion, no FEM — the substrate (but no contact)

- **Akyürek & Sahillioğlu 2025, Spring Decomposed Skinning (CGF/SGP)** — spring bones on the *existing* rig +
  helper bones, PBD, no volumetric mesh, real-time. **Contact is a documented-but-empty slot** in its PBD
  force term. Closest substrate to the target idea. Digest:
  [Spring-Decomposed-Skinning.md](./Spring-Decomposed-Skinning.md).
  <https://onlinelibrary.wiley.com/doi/10.1111/cgf.70209> · code <https://github.com/bartuakyurek/Spring-Decomposed-Skinning>
- **Rohmer et al. 2021, Velocity Skinning** — velocity-driven per-vertex secondary deformation, no simulation,
  no temporal integration. Explicitly **does not handle collisions**.
  <https://arxiv.org/pdf/2104.04934>

## 3. Contact-reactive spring bones with no FEM — already shipping (but crude)

Game/DCC spring-bone systems push bones out of sphere/capsule/box colliders; the skinned mesh follows. This
*is* "contact-reactive secondary rig motion, no elastic solve" — but **one-way** (collider repels bone, no
back-coupling), on **extra hand-placed bones**, **decoupled per-chain**, **hand-tuned**, and it is bone-position
*avoidance*, not surface-contact-driven deformation. It never confronts the representability gap.

- Magica Cloth / BoneCloth — <https://magicasoft.jp/en/paramater-collider-collision-2/>
- BoneDynamics Pro (Blender) — <https://superhivemarket.com/products/bonedynamics>
- Valve `$jigglebone` — <https://developer.valvesoftware.com/wiki/$jigglebone>
- Unity Dynamic Bone, Unreal physics-asset spring bones (e.g. God of War Ragnarök muscle/tissue heuristics).

## 4. Two-way contact on a rig with secondary DOFs — but with a reduced elastic solve

Closest *principled* matches to "secondary DOFs absorb contact + skeleton reacts," but all **solve a reduced
elastic energy** (not "no solve"):

- **Tapia, Romero, Pérez, Otaduy 2021, Parametric Skeletons with Reduced Soft-Tissue Deformations (CGF)** —
  SMPL/MANO + subspace soft tissue, **two-way coupling**, frictional contact, interactive. Uses a
  physics-based skinning-subspace elastic model. <https://mslab.es/projects/PSfSTD/>
- **Teng et al. 2015, Fully momentum-conserving reduced deformable bodies with collision, contact,
  articulation, and skinning (SCA)** — <https://dl.acm.org/doi/10.1145/2786784.2786787>
- **Capell et al. 2002, Interactive skeleton-driven dynamic deformations** —
  <https://dl.acm.org/doi/10.1145/566654.566622>

## 5. Subspace contact + the representability gap (why it's hard)

These define the central tension: a global subspace can't express *local* contact deformation, so they add
local full-space DOFs (a solve) — the opposite of "no solve."

- **Harmon & Zorin 2013, Subspace Integration with Local Deformations** — *the canonical representability-gap
  paper.* Global modal subspace augmented on-the-fly with **local basis vectors derived from analytic
  point-load (Boussinesq) solutions** at contact points; output-sensitive (#local DOFs ∝ #contact regions),
  runtime cubature for the changing basis, basis-change continuity via `q̄ = ŪᵀUq`, ~2 orders over full FEM.
  **Still a reduced-FEM solve** (needs tet mesh + precomputed modes; global part is physics modes, not an
  artist rig) → squarely the "solve" side, but its **analytic Boussinesq local function is the mechanism the
  chosen direction borrows** (see §11 and [Research-Direction](./Research-Direction.md)). PDF in repo:
  [`Harmon and Zorin - 2013 - Subspace integration with local deformations.pdf`](./Harmon%20and%20Zorin%20-%202013%20-%20Subspace%20integration%20with%20local%20deformations.pdf).
  Limitation (their §8): great for *very local + global*, degrades in the intermediate regime as the contact
  patch approaches object size.
- **Teng, Otaduy, Kim 2014, Simulating Articulated Subspace Self-Contact** — articulation makes self-contact
  predictable/low-rank; interactive on 100K–400K elements.
  <http://www.tkim.graphics/SASS/TengOtaduyKim2014.pdf> · <https://dl.acm.org/doi/10.1145/2601097.2601181>
- **Teng, Meyer, DeRose, Kim 2015, Subspace Condensation** — activate full-space DOFs locally near novel
  collisions, two-way coupled to the subspace. <http://www.tkim.graphics/CONDENSE/> ·
  <https://dl.acm.org/doi/10.1145/2766904>
- **Xu & Barbič 2016, Pose-Space Subspace Dynamics** — secondary FEM dynamics under rigged motion at ms/frame,
  self-contact via modes/derivatives under contact constraints.
  <https://viterbi-web.usc.edu/~jbarbic/multiModal/XuBarbic-SIGGRAPH2016.pdf>
- **Trading Spaces 2024, Adaptive Subspace Time Integration for Contacting Elastodynamics** —
  <https://dl.acm.org/doi/10.1145/3687946>
- **Embedded IPC 2024** — reduced elasticity + embedded high-res collision, intersection-free.
  <https://arxiv.org/pdf/2409.16385>

## 6. Learning-based contact corrections (the "learn the gap" answer)

- **Romero, Casas, Pérez, Otaduy 2021, Learning Contact Corrections for Handle-Based Subspace Dynamics
  (SIGGRAPH)** — handle/rig-like subspace dynamics + learned nonlinear corrections, decoupling internal vs.
  external contact-driven corrections. **Most directly comparable to the target idea.**
  <http://crisrom002.gitlab.io/files/papers/RCPO21.pdf> · <https://dl.acm.org/doi/10.1145/3450626.3459875> ·
  code <https://github.com/crisrom002/hsd-learning-contact-corrections>
- **Romero et al. 2022, Contact-Centric Deformation Learning** —
  <https://mslab.es/projects/ContactCentricLearning/contents/Romero_SIG2022_final.pdf>

## 7. PBD / shape-matching / oriented particles (contact via projection, but particles *are* the model)

- **Müller & Chentanez 2011, Adding Physics to Animated Characters with Oriented Particles** — skeleton +
  oriented particles, collisions, self-collision, volume; shape matching *is* the (cheap) elastic solve, on
  particles not bones. <https://matthias-research.github.io/pages/publications/animParticles.pdf>
- **Müller et al. 2007, Position Based Dynamics** — the integration framework SDS and many others build on.
  <https://matthias-research.github.io/pages/publications/posBasedDyn.pdf>

## 8. Geometric / implicit contact (no dynamics)

- **Vaillant et al. 2013, Implicit Skinning** — real-time skin-on-skin contact and bulging via implicit
  surfaces, no solve at all; the "simulate without simulating" geometric extreme.
  <https://inria.hal.science/hal-00819270/>

## 9. Complementary dynamics (came up in SDS comparisons)

- **Zhang et al. 2020, Complementary Dynamics** — secondary motion orthogonal to the rig's deformation.
- **Benchekroun et al. 2023, Fast Complementary Dynamics via Skinning Eigenmodes** — real-time variant.
- **Wu & Umetani 2023** — PBD on a tet mesh in a complementary subspace; SDS reports ~10× faster and argues it
  is more stable (no mesh simulation).

## 10. Adaptive DOF enrichment (candidate gap-fix; prior art for "spawn DOFs at contact")

- **Grinspun et al. 2002, CHARMS** — refine *basis functions*, not elements; a localized spring bone is a
  refinable basis function.
- **Debunne et al. 2001** — adaptive real-time multiresolution soft body.
- **Harmon & Zorin 2013** — contact-triggered local enrichment of a subspace (detailed in §5; closest to "add
  DOFs at the contact patch").
- **Narain et al. 2012** — adaptive anisotropic remeshing for cloth.

## 11. Analytic / precomputed-response contact — the chosen route's ancestors

The "supply the local contact deformation in closed form, no solve" family. **This is the lineage the chosen
direction must beat** (see [Research-Direction](./Research-Direction.md)).

- **James & Pai 1999, ArtDefo: Accurate Real-Time Deformable Objects** — *the real ancestor.* Precomputed
  boundary-element **Green's functions** give exact-feeling, real-time contact deformation of an elastic body.
  What it lacks (our wedge): an **editable animation rig**, **secondary dynamics**, and **mesh-free /
  on-demand analytic** contact (it precomputes per-body Green's functions).
  <https://graphics.stanford.edu/papers/artdefo/>
- **James & Pai 2003, Multiresolution Green's Function Methods for Interactive Simulation of Large-Scale
  Elastostatic Objects** — multiresolution scaling of the same Green's-function idea.
- **Pauly, Pai, Guibas 2004, Quasi-Rigid Objects in Contact** — uses the **Boussinesq** half-space solution
  for contact on point clouds; direct precedent for "analytic contact field as the local response."
- **Boussinesq / contact mechanics (Johnson 1987)** — the closed-form half-space point/area-load solution
  itself; the material-aware (`G`, `ν`) displacement field the chosen route evaluates instead of simulating.

## Where it's open

| | existing rig (not extra bones) | two-way / contact-driven | handles rig-can't-represent case | no elastic solve |
|---|---|---|---|---|
| Spring Decomposed Skinning '25 | ✅ | ❌ (no contact) | ❌ | ✅ |
| Game spring bones + colliders | ➖ (extra bones) | ❌ (one-way) | ❌ | ✅ |
| Tapia '21 | ✅ | ✅ | partial | ❌ (reduced FEM) |
| Oriented particles '11 | ➖ | ✅ | ➖ | ❌ (shape-match) |
| Subspace Condensation / Trading Spaces | ✅ | ✅ | ✅ (adds full DOFs) | ❌ |
| Harmon & Zorin '13 | ❌ (FEM modes) | ✅ | ✅ (analytic local basis) | ❌ (reduced FEM) |
| ArtDefo '99 | ❌ (no rig/anim) | ✅ | ✅ (Green's functions) | ✅ (closed-form) |
| Romero '21 (learned corrections) | ✅ | ✅ | ✅ (learned) | ➖ (subspace + net) |
| **chosen direction** | **✅** | **✅** | **✅ (closed-form field)** | **✅** |

**No prior row hits all four; the chosen direction aims to.** The unclaimed slot: contact resolved as a
surface→bone projection on the **existing editable rig** (the `Jᵀ` / skinning-weight-transpose bridge),
letting **secondary DOFs absorb the gross reaction**, while the local deformation the rig can't represent is
supplied by a **closed-form contact-mechanics field** (Boussinesq) — no volumetric elastic energy. ArtDefo
gets the last two columns but has no rig/animation; Harmon-Zorin gets the response but pays for a reduced-FEM
solve on precomputed modes.
