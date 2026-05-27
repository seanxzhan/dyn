"""Elastic energy on the tet mesh: the StVK material model.

This is the potential ``W(x)`` whose negative gradient is the elastic force and
whose Hessian is the stiffness.  Everything is vectorized over tets and checked
against finite differences in ``tutorials/04``.

Per tetrahedron, with rest/deformed edge matrices ``Dm, Ds`` (columns
``v_i - v_0``):

    F = Ds Dm⁻¹                      deformation gradient (3x3)
    E = ½(FᵀF − I)                   Green strain
    Ψ = μ‖E‖_F² + ½λ tr(E)²          St. Venant–Kirchhoff energy density
    W = Σ_e Ψ_e · V0_e               total elastic energy

Forces come from the first Piola–Kirchhoff stress ``P = F·(2μE + λ tr(E) I)``;
the standard result is that the energy gradient w.r.t. nodes 1,2,3 is the matrix
``V0 · P Dm⁻ᵀ`` (columns), and the node-0 gradient is minus their sum.
"""

from __future__ import annotations

import numpy as np

_I3 = np.eye(3)


def lame_params(young: float, poisson: float) -> tuple[float, float]:
    """Convert (Young's modulus E, Poisson ratio ν) to Lamé (μ, λ)."""
    mu = young / (2.0 * (1.0 + poisson))
    lam = young * poisson / ((1.0 + poisson) * (1.0 - 2.0 * poisson))
    return mu, lam


class StVK:
    """St. Venant–Kirchhoff elastic energy for a fixed tet mesh."""

    def __init__(self, beam, young: float = 200.0, poisson: float = 0.40):
        self.tets = beam.tets
        self.mu, self.lam = lame_params(young, poisson)
        self.young, self.poisson = young, poisson

        v = beam.verts
        # Rest edge matrix Dm (T,3,3): columns (V1-V0, V2-V0, V3-V0).
        self.Dm = np.stack([v[self.tets[:, 1]] - v[self.tets[:, 0]],
                            v[self.tets[:, 2]] - v[self.tets[:, 0]],
                            v[self.tets[:, 3]] - v[self.tets[:, 0]]], axis=-1)
        self.DmInv = np.linalg.inv(self.Dm)
        self.V0 = np.abs(np.linalg.det(self.Dm)) / 6.0          # rest volumes (T,)

    # --- kinematics --------------------------------------------------------
    def deformation_gradient(self, x: np.ndarray) -> np.ndarray:
        """F = Ds Dm⁻¹ for every tet.  x:(N,3) -> F:(T,3,3).  Rest pose -> F=I."""
        v0 = x[self.tets[:, 0]]
        Ds = np.stack([x[self.tets[:, 1]] - v0,
                       x[self.tets[:, 2]] - v0,
                       x[self.tets[:, 3]] - v0], axis=-1)
        return Ds @ self.DmInv

    def green_strain(self, F: np.ndarray) -> np.ndarray:
        """E = ½(FᵀF − I), the nonlinear (rotation-invariant) strain."""
        FtF = np.einsum("tji,tjk->tik", F, F)
        return 0.5 * (FtF - _I3)

    # --- energy, gradient (force), and per-tet density ---------------------
    def energy_density(self, F: np.ndarray) -> np.ndarray:
        """Ψ per tet (energy per unit rest volume)."""
        E = self.green_strain(F)
        trE = np.trace(E, axis1=1, axis2=2)
        normE2 = np.einsum("tij,tij->t", E, E)
        return self.mu * normE2 + 0.5 * self.lam * trE ** 2

    def energy(self, x: np.ndarray) -> float:
        return float(np.sum(self.energy_density(self.deformation_gradient(x)) * self.V0))

    def gradient(self, x: np.ndarray) -> np.ndarray:
        """∂W/∂x as an (N,3) field (elastic force is the negative of this)."""
        F = self.deformation_gradient(x)
        E = self.green_strain(F)
        trE = np.trace(E, axis1=1, axis2=2)
        S = 2.0 * self.mu * E + self.lam * trE[:, None, None] * _I3        # 2nd PK stress
        P = np.einsum("tij,tjk->tik", F, S)                                # 1st PK stress
        # G (T,3,3): columns are ∂W/∂x1, ∂W/∂x2, ∂W/∂x3  =  V0 · P Dm⁻ᵀ
        G = self.V0[:, None, None] * np.einsum("tij,tkj->tik", P, self.DmInv)
        grad = np.zeros_like(x)
        for a in range(3):
            np.add.at(grad, self.tets[:, a + 1], G[:, :, a])
        np.add.at(grad, self.tets[:, 0], -(G[:, :, 0] + G[:, :, 1] + G[:, :, 2]))
        return grad

    def forces(self, x: np.ndarray) -> np.ndarray:
        return -self.gradient(x)

    def per_vertex_energy(self, x: np.ndarray) -> np.ndarray:
        """Spread each tet's energy to its 4 nodes (for visualization)."""
        we = self.energy_density(self.deformation_gradient(x)) * self.V0
        out = np.zeros(x.shape[0])
        np.add.at(out, self.tets.reshape(-1), np.repeat(we / 4.0, 4))
        return out
