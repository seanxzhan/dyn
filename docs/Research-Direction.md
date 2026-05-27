# Research Direction — Contact-Reactive Rig-Space Dynamics

> Working notes for the `dyn` project: a SIGGRAPH-aimed paper on real-time, physically-grounded
> secondary motion that lives on a character's rig. This is the opinionated "where we are / where we're
> going" doc. Neutral literature survey lives in
> [Related-Work — Contact in Reduced & Rig-Space Models](./Related-Work-Contact-Reduced-Rigs.md).
> Paper digests: [Rig-Space Physics](./Rig-Space-Physics.md) ·
> [Efficient Secondary Motion](./Efficient-Secondary-Motion-Rig-Space.md) ·
> [Spring Decomposed Skinning](./Spring-Decomposed-Skinning.md) ·
> [Spring Bones vs. Rig-Space Physics](./Spring-Bones-vs-Rig-Space-Physics.md).

## The goal

Secondary motion (jiggle, sway, squash, follow-through) that is **all three** of:

1. **Real-time** — engine-deployable, cost tied to #handles not #vertices.
2. **Physically grounded** — real elastic behavior, not ad-hoc per-bone springs.
3. **Editable / automatic** — output lives on the rig; artists keep control; works on arbitrary
   production rigs.

No existing method is all three at once ("pick-any-two" gap): FEM rig-space physics is physical+editable but
offline; jiggle/spring bones are real-time+editable but not physical and hand-tuned; neural emulators are
real-time+physical but not editable and shape-specific.

## The taste filter (apply to every idea)

> **Does the idea do something hand bones fundamentally *cannot* (a new capability or phenomenon), or does it
> just do the same thing more *correctly*?**

"Same effect as hand bones but physically correct" is a **weak pitch** — correctness is not the unmet need.
Artists keep hand bones for *control and stylization*, not physics, and [DyRT (James & Pai 2002)] already
won the "physically-correct jiggle from precomputed modes" argument two decades ago. Lead with a single
clean insight and an obvious "who cares," not a stack of qualifiers (see the [[research-taste-clean-ideas]]
principle). Resist re-stacking correction modules onto a clean kernel just to manufacture novelty.

## The contact angle (current focus)

The capability hand bones genuinely lack — and the cleanest reason secondary motion would *need* real
physics — is **contact responsiveness**: when the mesh touches something, the rig reacts and that drives a
plausible deformation, *without* a full simulation. The pitch: exploit how fast rigged deformation is to
make it **responsive to contact**.

### Why this is defensible (and where it isn't)

- **Contact in reduced/subspace models is already well studied** — subspace self-contact, reduced bodies
  with contact, complementary dynamics, neural subspace + collision. Self-collision *avoidance* is **not** a
  novel headline. (See the survey doc.)
- **Contact-reactive spring bones with no FEM already *ship*** — game engines push spring bones out of
  sphere/capsule colliders (Magica Cloth, Dynamic Bone, Unreal). So "spring bones that collide" is not novel.
- **The closest academic substrate is [Spring Decomposed Skinning](./Spring-Decomposed-Skinning.md)
  (CGF 2025):** spring bones on the *existing* rig + helper bones, PBD, no FEM, real-time. Crucially, its
  PBD force term is *documented to admit contact forces but they omit them* → **"just add contact to the
  force" is engineering they already sanctioned, not a contribution.**

### Where the actual contribution must live

SDS gives the substrate for free, so the novelty has to be in the two things it can't:

1. **The surface→bone contact bridge.** Contact happens on *mesh vertices*; the DOFs are *bones*. Mapping a
   contact force back to bones is the **transpose of the skinning map** (`Σⱼ wᵢⱼ`) — the LBS analogue of
   Hahn's `Jᵀ` generalized-force projection in [Rig-Space Physics](./Rig-Space-Physics.md). SDS never builds
   this because it never has surface forces. Concrete and unclaimed.
2. **The representability-gap handling.** A 3-DOF translational spring bone cannot form a local dent where a
   finger presses if no bone lives there. A global/low-rank rig fundamentally *cannot express local,
   high-rank contact deformation*. This is THE central tension. Everyone who takes it seriously today either
   adds a solve (Subspace Condensation, Trading Spaces) or learns a correction (Romero 2021). Our answer:
   don't represent it on the rig at all — supply the local deformation as a **closed-form contact-mechanics
   field** (Boussinesq), so it costs no DOFs and no solve (see the pitch below).

### One-line pitch

> **Stop trying to make the rig represent contact — split it: let the fast, editable rig absorb the gross
> global reaction (via a surface→bone force projection) and let a closed-form contact-mechanics field supply
> the local, material-aware squish, giving physically-grounded contact deformation in real time with no
> volume ever simulated.**

**The single insight (one idea, not a stack):** contact deformation has a known structure — a smooth global
reaction plus a fast-decaying local boundary layer — and a rig is the *wrong* basis for the boundary layer
but already perfect for the global part; so don't *simulate* the boundary layer, *evaluate* it in closed form
(Boussinesq / contact mechanics) and hand the rest to the rig. Both halves fall out of that one reframe.

**Why it survives the lit review:**

| closest work | what it has | what this adds |
|---|---|---|
| [Spring Decomposed Skinning](./Spring-Decomposed-Skinning.md) (2025) | mesh-free spring-bone secondary motion on the rig | the contact it explicitly lacks |
| game spring bones + colliders | contact-reactive, no FEM | a *material-aware deformation shape*, not bone-position avoidance |
| [Harmon & Zorin 2013](./Related-Work-Contact-Reduced-Rigs.md) | analytic local contact enrichment (Boussinesq) | done on an **editable rig with no FEM backbone**, not a reduced-FEM modal solve |
| **ArtDefo (James & Pai 1999) — the real ancestor** | closed-form (Green's-function) contact deformation in real time | an **animatable rig + secondary dynamics**, not a fixed precomputed elastostatic body |

**Passes the taste filter:** it does something hand bones *fundamentally cannot* — a dimple whose shape
depends on the material's stiffness / Poisson ratio, plus a two-way reaction that absorbs energy and changes
the jiggle — rather than doing jiggle "more correctly" (the DyRT trap).

**Honest caveats (the fine print):**

- **Nearest ancestor is ArtDefo** (James & Pai 1999, *ArtDefo: Accurate Real-Time Deformable Objects*):
  precomputed boundary-element **Green's functions** give real-time contact deformation of an elastic body.
  We must cite it and beat it on three axes it lacks — an **editable animation rig**, **secondary dynamics**,
  and **mesh-free / on-demand analytic** contact (no per-body Green's-function precompute). Analytic-response
  lineage to position within: James & Pai 2003 (multiresolution Green's functions), Pauly et al. 2004
  (Boussinesq for point-cloud contact), Harmon & Zorin 2013 (Boussinesq basis in a reduced-FEM subspace).
- **Sweet-spot regime only** (per Harmon–Zorin §8): crisp for a fingertip poke; degrades as the contact patch
  approaches body size (e.g. sitting flat on a bench), where the local-boundary-layer assumption breaks.

## The internal tension to confront honestly

"No elastic solve" and "handle where the rig can't represent contact" pull against each other: the reason the
rig *can't* represent the contact is exactly the high-rank local deformation that normally needs extra DOFs
(a solve) to express. Three ways out, in order of preference:

- **Closed-form local field (the chosen route — this is the pitch).** Don't put the local deformation on the
  rig *or* solve for it — *evaluate* it from the analytic contact-mechanics solution (Boussinesq point/area
  load), aligned to the contact normal and scaled by the material's shear modulus / Poisson ratio, layered on
  the rig-driven surface. Physically grounded, material-aware, no DOFs, no solve. This is
  [Harmon & Zorin 2013]'s analytic-basis idea, but used as a *displacement layer on a mesh-free editable rig*
  instead of an added vector in a reduced-FEM modal subspace. Risk to own: H&Z can fall back to the exact
  `Ku=f`; we have no `K`, so we lean entirely on the analytic field's validity (good in the local sweet spot,
  not beyond).
- **Concede gracefully** — where even the analytic field is wrong (large/conforming contact), project the
  contact onto the nearest rig-expressible motion and *characterize* when pure rig-space contact looks right.
- **Adaptive enrichment** — spawn extra spring-bone DOFs *at the contact patch* on the fly, retire on release;
  the principled fix when the contact must alter the **global dynamics** (energy absorption) — which a
  cosmetic layer can't capture (Harmon–Zorin's bunny and bumpy-box examples show this matters). A localized
  spring bone *is* a refinable basis function (à la CHARMS). Prior art for adaptive DOFs (CHARMS 2002, Debunne
  2001, Harmon–Zorin 2013, adaptive cloth 2012) means "adaptive DOFs" alone is **not** novel; the open white
  space is contact-triggered enrichment of a real-time + editable + arbitrary-rig model with
  **guaranteed-continuous spawn/retire** (no popping; momentum/energy transfer — H&Z's `q̄ = ŪᵀUq`
  reprojection is a starting point) and a **cheap/learned trigger** (error indicators usually need the full
  solve; contact location is cheap). Reviewer-proofing: show adaptive overhead beats a fixed richer basis,
  especially at many simultaneous contacts.

## Other live clean-kernel candidates (not contact-specific)

- **K1 — Spring bones as a compilation target.** Auto-compile an FEM character into an engine-native
  (Unity/Unreal) spring-bone rig within an error bound. Best *adoption* story.
- **K3 — Two-way spring bones.** Passive back-coupling so the secondary mass pushes back on primary motion
  (everything today is strictly one-way driven). Most novel *problem*; scope risk toward character control.
- Leading bets before the contact focus were **K1 (adoption)** and **K3 (real gap)**. Contact responsiveness
  is the differentiator that makes the *forward model* worth more than jiggle bones; K3 is closely related
  (contact is one source of back-coupling).

## The physics throughline (how to stay principled)

**Modal analysis** is the physically-correct bridge between cheap spring bones and FEM: each low-frequency
eigenmode (`Kφ = ω²Mφ`) is a decoupled spring-damper, so a *physically-calibrated* spring bone fits each
spring's stiffness (`k ← ω²`) and damping (`c ← Rayleigh`) to a real elastic mode and places bones at mode
antinodes. Pose generalization is the open problem (linear modes are valid only near rest); fixes = modal
derivatives (Barbič–James), corotational/modal warping, or pose-conditioned neural modes.

## Downstream applications to keep in view (not required for paper #1)

A fast, pose-generalizable, *differentiable* forward model of contact-reactive rig dynamics also powers
inverse design / fabrication: soft-robot design+control, compliant mechanisms/metamaterials, animatronics
(inverse-design internal stiffness so passive jiggle matches a CG target), prosthetics, garments. Contact is
*essential* for the soft-robotics control/design tie (gripper conforming to a target grasp; animatronic
flesh squishing under load). Differentiable contact in reduced rig-space (non-smooth → hard) is the natural
"new problem" extension.

## Open questions / next steps

- Novelty check the front-runner against EG 2025/26, SCA, SIGGRAPH Asia, and arXiv before building.
- Decide the representability-gap stance: **concede-and-characterize** vs. **adaptive enrichment**.
- Prototype the surface→bone contact projection on top of the SDS substrate (its code is public).
- Pressure-test against the taste filter: name the capability hand bones can't match.
