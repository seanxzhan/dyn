# Real-Time Secondary Animation with Spring Decomposed Skinning

> Bartu Akyürek, Yusuf Sahillioğlu.
> *Computer Graphics Forum 44(5), Eurographics Symposium on Geometry Processing (SGP) 2025.*
> Dept. of Computer Engineering, METU, Turkey.
> PDF: [`Akyürek and Sahillioğlu - 2025 - Real-Time Secondary Animation with Spring Decomposed Skinning.pdf`](./Aky%C3%BCrek%20and%20Sahillio%C4%9Flu%20-%202025%20-%20Real-Time%20Secondary%20Animation%20with%20Spring%20Decomposed%20Skinning.pdf)
> Code: <https://github.com/bartuakyurek/Spring-Decomposed-Skinning>

## TL;DR

Add secondary motion to a skinned character by **simulating springs on the rig bones themselves** — never on a
volumetric mesh. The artist marks some bones as **spring bones**; each carries a mass-spring whose head
("fixed mass") is driven by the keyframed animation and whose tail ("free mass") is integrated with
**Position-Based Dynamics**. Simulated bone tips are converted back to bone transforms via a
Rotate-Scale-Translate fit and fed to ordinary **linear blend skinning**. Cost scales with **#bones, not
#vertices** → real-time (4–14 ms/frame). Explicitly the *academic formalization of the industry "jiggle
bone."* **No FEM, no tetrahedralization, and — important for us — no contact handling.**

## The problem & positioning

- Skinning (LBS) is fast and intuitive but **static**: no jiggle, sway, squash, follow-through.
- FEM secondary-motion methods (incl. [Rig-Space Physics](./Rig-Space-Physics.md) and its
  [follow-up](./Efficient-Secondary-Motion-Rig-Space.md)) add real dynamics but pay for a **volumetric tet
  mesh**; cost grows with mesh resolution and tetrahedralizing an existing rigged asset is cumbersome.
- SDS keeps physics **in the space of rig handles** (like the rig-space line) but replaces the FEM continuum
  with **decoupled Hookean springs on bones** — trading physical fidelity (no volume preservation, no coupled
  flesh) for resolution-independent real-time speed and plug-and-play pipeline compatibility.

This is the same fidelity/speed trade as [Spring Bones vs. Rig-Space Physics](./Spring-Bones-vs-Rig-Space-Physics.md),
but SDS is the *principled, rig-integrated* version of jiggle bones rather than the ad-hoc engine feature.

## Method

### Spring bones: fixed + free mass (§3, Fig. 2)

Each spring bone is a mass-spring along the bone segment:

- **Fixed mass** at the head — **not simulated**; its position is set by forward kinematics / keyframes +
  kinematic constraints. Carries the *primary* animation. (Mass `m_i = 0` ⇒ skipped in the solver.)
- **Free mass** at the tail — **simulated** each frame; carries the *secondary* motion.
- **Point spring bone**: both masses at the tail via a zero-rest-length spring → enables squash/stretch.

Two artist roles for spring bones:

- **Primary bones** converted to spring bones → **global** secondary motion (whole-body floppy/squash).
- **Helper bones** added by the artist → **local** secondary motion (e.g. belly-fat jiggle on SMPL). Helper
  bones are placed nearly *orthogonal* to the primary bones, echoing Complementary Dynamics' orthogonality idea.

### Forces & PBD integration (§3.1, Alg. 1)

Hookean spring + damping on particle `p₁` (spring vector `p_{1,2}=p₂−p₁`, `n = p_{1,2}/‖p_{1,2}‖`):

```
f_s = k_s(‖p_{1,2}‖ − l₀)·n          (eq. 1)
f_d = −k_d[ n·(v₂+v₁) ]·n            (eq. 2)
f   = f_s + f_d        (on p₂: −f_s)
```

Integrated with **PBD** [Müller et al. 2007] for unconditional stability under large steps; extra velocity
damping `dₛ ∈ (0,1)`. An optional **stretch constraint** (eq. 3) re-projects spring length toward `l₀`.

> ⚠️ **The hook for our work:** the paper states the external force `f` *"can include gravity and other
> **contact forces**"* — but they **omit both**, "focusing solely on the effect of spring forces." Contact is
> a named-but-unused slot in their force term.

### Kinematic constraints — applied *after* the solve (§3.2, Alg. 2, Fig. 4)

Simulating only free masses **breaks bone connectivity** (a child's head no longer sits at its parent's tail).
They restore the kinematic tree (and optionally rest bone lengths via `fixed_scale`) with constraints
projected **outside** the PBD loop. Key empirical finding: projecting these *inside* PBD causes **abrupt
vibrations** (Fig. 4b); projecting them *after* yields smooth jiggling (Fig. 4c). Bending of a chain emerges
from per-bone 3-DOF translation + this connectivity restoration, not from any per-spring bending DOF.

### Back to bone transforms: RST fit (§3.3, Appendix A)

Skinning needs matrices `Mᵢ`, but the solver outputs *positions*. They fit a **Rotate-Scale-Translate**
affine map taking each rest bone segment to its simulated segment (translation = tail displacement; rotation
via **Rodrigues' formula** aligning rest→posed bone vectors; scale = length ratio). Then standard LBS:

```
v_j = Σ_i w_{ij} Mᵢ v_j⁰            (eq. 4)
```

`Mᵢ = MᵢᴰMᵢᴿ` keeps the rigid FK part `Mᴿ` for non-spring bones and overrides with the dynamic part `Mᴰ`
for spring bones (Alg. 3).

## Results

- **Real-time, resolution-independent** (Table 1, 2.5 GHz i7): Cloth 12.2 ms, Monstera 14.1 ms, Duck 9.5 ms,
  SMPL 4.7 ms, Elephant ~8–12 ms. Cost tracks #bones, not #vertices (Duck has 12.9k verts but only 5–6 bones).
- Range of effects (Fig. 1): plant-leaf sway (composite rigid), cow/duck stretch-squash (elastic), SMPL
  belly-fat jiggle (soft tissue), whole-elephant wobble (global), cloth.
- **~10× faster than Wu & Umetani 2023** (PBD on a tet mesh in a complementary subspace), which it argues is
  prone to instability and detail-dependent stiffness because it simulates the mesh.
- Compared against **Fast Complementary Dynamics** [Benchekroun et al. 2023] and **Complementary Dynamics**
  [Zhang et al. 2020]: SDS gives more direct artistic control by simulating bones rather than orthogonal mesh
  motion (Fig. 18).

## Limitations (their §5.3)

- **No collision detection / no external (incl. contact) forces — left as future work.** Long helper-bone
  chains self-intersect with nothing to stop them.
- **Only 3 DOF per spring** (x/y/z translation); no per-spring twist or bend.
- Hookean springs can still **explode** under bad parameters even within PBD (position-based, not force-based).
- **Hand-tuned** stiffness/damping; automating them is open.
- Single rest mesh; decoupled springs ⇒ **no volume preservation, no coupled flesh** deformation.

## Relationship to our direction (contact-reactive secondary rig motion)

SDS is the **closest existing substrate to the idea we're chasing**, and reading it sharpens the gap:

1. **It already does "secondary motion on the existing rig, no elastic solve, real-time."** Primary bones →
   global, helper bones → local; free masses absorb the dynamics while fixed masses carry the artist's
   animation. That is exactly "let secondary/auxiliary DOFs absorb the motion," minus contact.
2. **Contact is a one-line slot they left empty.** Their PBD force `f` admits contact forces; they just never
   use them. So naively "add contact penalties to `f`" is *engineering, not a contribution* — the paper
   already says it's possible.
3. **The real open problem is the surface↔bone gap under contact.** Contact lives on *mesh vertices*; the DOFs
   are *bones*. Mapping a contact force from the surface to the bones is the **transpose of the skinning map**
   (the `Jᵀ` projection from [Rig-Space Physics](./Rig-Space-Physics.md), here `Σ_j w_{ij}` instead of an FEM
   Jacobian). SDS never builds this bridge because it never has surface forces.
4. **The representational-mismatch question survives intact.** A spring bone is 3-DOF translation; it cannot
   form a local dent where a finger presses if no bone lives there. SDS sidesteps this (no contact); we can't.
   The clean contribution is *what to do at that mismatch without a volumetric solve* — the thing
   [Subspace Condensation / Trading Spaces] answer by adding full-space DOFs and [Romero et al. 2021] answer by
   learning a correction.
5. **Kinematic-constraint-after-solve is a useful trick to inherit:** contact response would likewise want to
   avoid being fought by the connectivity/IK projection — order of operations matters.

**Net:** SDS is the natural baseline and code to build on, and the framing "SDS + principled contact" is
concrete — but the novelty must be in the surface→bone contact mapping and the mismatch handling, not in the
mere act of adding contact forces.
