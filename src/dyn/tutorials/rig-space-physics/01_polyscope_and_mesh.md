# Lesson 1 — The discretized continuum and the DOF split

> Run `python src/dyn/tutorials/01_polyscope_and_mesh.py`.
> Code: [`01_polyscope_and_mesh.py`](01_polyscope_and_mesh.py) · Library: [`../mesh.py`](../mesh.py)

This lesson has almost no physics. Its only job is to build the object every
later lesson stands on, and to fix the vocabulary the papers use.

## 1. From a continuum to a tet mesh

A soft body is a chunk of material occupying a volume $\Omega$. Its motion obeys
Newton's law as a *continuum*,

$$ \rho\,\ddot{x} = f(x), \qquad f = -\frac{\partial W}{\partial x}, $$

where $\rho$ is density, $x$ is the (infinite-dimensional) deformation field, and
$W$ is potential energy (elastic + external). We cannot store an infinite field,
so we **discretize space**: tile $\Omega$ with tetrahedra and track only the
vertex positions. Inside each tet the deformation is interpolated linearly from
its 4 corners (linear tet basis functions). The continuous unknown becomes a
finite vector $x \in \mathbb{R}^{3N}$ of $N$ vertex positions.

We use a **beam** (an elongated cuboid) because it is the cleanest soft body to
reason about: clamp one end and it sags/jiggles like a diving board. `make_beam`
builds it analytically — a regular grid of vertices, each voxel split into 6
tetrahedra by the **Kuhn subdivision** (the same diagonal in every voxel, so the
shared faces match and the mesh has no cracks). No external mesher; the geometry
is known in closed form, which we exploit later to check our math.

## 2. Lumped mass

The mass matrix $M$ comes from integrating $\rho$ against the basis functions.
The exact ("consistent") mass matrix couples neighbouring vertices. Like the
papers, we use the cheaper **lumped mass**: a diagonal matrix where each vertex
gets a quarter of the mass of every tet that touches it,

$$ M_{ii} = \rho \cdot \tfrac14 \!\!\sum_{e \ni i}\! V_e . $$

The script checks $\sum_i M_{ii} = \rho\,\mathrm{vol}(\Omega)$ — mass is
conserved. We also orient every tet so its signed volume
$V_e = \tfrac16\det[\,v_1-v_0,\;v_2-v_0,\;v_3-v_0\,] > 0$; consistent orientation
matters once we compute strain energy.

## 3. The DOF split: surface vs. interior

This is the one idea to take away. Partition the vertices into

$$ x = \{s\} \;\cup\; \{n\}. $$

- **Surface nodes $s$** lie on the boundary $\partial\Omega$. In rig-space
  physics *these are the ones the rig moves*: the animator's rig is a map
  $p \mapsto s(p)$ to the surface. (Lesson 2.)
- **Interior nodes $n$** are strictly inside. The rig says nothing about them;
  they are free FEM degrees of freedom that the physics resolves. (They become
  central in the dynamics lessons; the 2013 paper eliminates them with a learned
  skinning map $q = Ws$.)

A handful of vertices are additionally **clamped** (Dirichlet): the $x=0$ face is
pinned, turning the beam into a cantilever. Clamped nodes never move, so they
drop out of the unknowns.

For the default `res=(10,3,3)` beam: $N=176$ vertices, $T=540$ tets, with 140
surface / 36 interior / 16 clamped nodes.

## What to look at

- Turn on the **node type** quantity (it colors interior=0, surface=1, clamp=2).
  Notice the interior nodes form an inner core; the rig will only ever touch the
  outer shell.
- The **red point cloud** is the clamped end. Hold this picture: in later lessons
  the free end will lift, bend, and jiggle while these points stay frozen.

## Where this is going

| symbol | meaning | introduced |
|---|---|---|
| $x \in \mathbb{R}^{3N}$ | all vertex positions | here |
| $s,\;n$ | surface / interior nodes | here |
| $p$ | rig parameters (reduced coords) | Lesson 2 |
| $J = \partial s/\partial p$ | rig Jacobian | Lesson 3 |
| $W,\;\partial W/\partial x$ | elastic energy, forces | Lesson 4 |
| $H = \tfrac{h^2}{2}\text{(inertia)} + W$ | incremental potential | Lesson 6 |

Reference: Hahn et al. 2012, *Rig-Space Physics* — see
[`../../../docs/Rig-Space-Physics.md`](../../../docs/Rig-Space-Physics.md).
