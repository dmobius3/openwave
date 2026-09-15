"""
Nonlinear restoring forces. Each processor adds a term to the acceleration
accumulator psi_new. Runs in the UPDATE stage after Laplacian, before Leapfrog.
"""

from ..pipeline import BaseProcessor, Stage

import taichi as ti

from .features import PsiField, WaveGrid


class NonlinearCubic(BaseProcessor):
    """
    Cubic self-interaction: adds -gamma * |psi|^2 * psi to psi_new.

    EWT coupling: gamma = 1/eps_M = N_geom * pi^3 (order 2.4e4 in EWT units).
    Start with small gamma (0.01 .. 1.0) to see the qualitative effect,
    scale up once the solver is proven stable.

    Must run AFTER a processor that initializes psi_new (i.e. Laplacian),
    and BEFORE Leapfrog. The order=15 slot is reserved for that.
    """

    name = "NonlinearCubic"
    stage = Stage.UPDATE
    order = 15
    requires = (WaveGrid, PsiField)

    def __init__(self, gamma: float):
        self.gamma = gamma

    def process(self, ctx) -> None:
        grid = ctx.data.require(WaveGrid)
        field = ctx.data.require(PsiField)
        _apply_cubic(field.psi, field.psi_new, self.gamma, grid.nx, grid.ny, grid.nz)


@ti.kernel
def _apply_cubic(
    psi: ti.template(),
    accel: ti.template(),
    gamma: ti.f32,
    nx: ti.i32,
    ny: ti.i32,
    nz: ti.i32,
):
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        u = psi[i, j, k].norm_sqr()
        accel[i, j, k] += -gamma * u * psi[i, j, k]
