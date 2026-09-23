"""Outer boundary condition. Dirichlet psi = 0 for now."""

from ..pipeline import BaseProcessor, Stage

import taichi as ti

from .features import PsiLongField, WaveGrid


class DirichletBoundaryProcessor(BaseProcessor):
    """Zero the outer shell of the grid every step."""

    name = "DirichletBoundary"
    stage = Stage.POST_UPDATE
    order = 10
    requires = (WaveGrid, PsiLongField)

    def process(self, ctx) -> None:
        grid = ctx.data.require(WaveGrid)
        field = ctx.data.require(PsiLongField)
        _dirichlet(field.psi, field.psi_prev, grid.nx, grid.ny, grid.nz)


@ti.kernel
def _dirichlet(
    psi: ti.template(),
    prev: ti.template(),
    nx: ti.i32,
    ny: ti.i32,
    nz: ti.i32,
):
    z = ti.Vector([0.0, 0.0, 0.0])
    for i, j, k in ti.ndrange(nx, ny, nz):
        on_edge = i == 0 or i == nx - 1 or j == 0 or j == ny - 1 or k == 0 or k == nz - 1
        if on_edge:
            psi[i, j, k] = z
            prev[i, j, k] = z
