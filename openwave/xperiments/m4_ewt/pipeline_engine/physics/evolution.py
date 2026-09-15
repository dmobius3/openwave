"""
Wave evolution: Laplacian then leapfrog.

Two separate processors, because the leapfrog needs the Laplacian already
written into psi_new as scratch. This keeps each processor single-purpose
and lets us swap the laplacian stencil independently of the time integrator.
"""

from ..pipeline import BaseProcessor, Stage

import taichi as ti

from .features import PsiField, WaveGrid


class LaplacianProcessor(BaseProcessor):
    """Writes the vector Laplacian of psi into psi_new (scratch)."""
    name = "Laplacian"
    stage = Stage.UPDATE
    order = 10
    requires = (WaveGrid, PsiField)

    def process(self, ctx) -> None:
        grid = ctx.data.require(WaveGrid)
        field = ctx.data.require(PsiField)
        _laplacian(field.psi, field.psi_new,
                   grid.nx, grid.ny, grid.nz, grid.dx)


@ti.kernel
def _laplacian(
    psi: ti.template(),
    out: ti.template(),
    nx: ti.i32, ny: ti.i32, nz: ti.i32,
    dx: ti.f32,
):
    inv_dx2 = 1.0 / (dx * dx)
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        face_sum = (
            psi[i + 1, j, k] + psi[i - 1, j, k]
            + psi[i, j + 1, k] + psi[i, j - 1, k]
            + psi[i, j, k + 1] + psi[i, j, k - 1]
        )
        out[i, j, k] = (face_sum - 6.0 * psi[i, j, k]) * inv_dx2


class LeapfrogProcessor(BaseProcessor):
    """
    Leapfrog update: psi_new = 2*psi - psi_prev + (c*dt)^2 * laplacian.
    Then swaps: prev <- psi, psi <- new.
    Requires the Laplacian already written into psi_new.
    """
    name = "Leapfrog"
    stage = Stage.UPDATE
    order = 20
    requires = (WaveGrid, PsiField)

    def process(self, ctx) -> None:
        grid = ctx.data.require(WaveGrid)
        field = ctx.data.require(PsiField)
        c2dt2 = (grid.c * ctx.sim.dt) ** 2
        _leapfrog(field.psi, field.psi_prev, field.psi_new,
                  grid.nx, grid.ny, grid.nz, c2dt2)


@ti.kernel
def _leapfrog(
    psi: ti.template(),
    prev: ti.template(),
    new: ti.template(),
    nx: ti.i32, ny: ti.i32, nz: ti.i32,
    c2dt2: ti.f32,
):
    # Update: read laplacian from new, write psi(t+dt) back into new.
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        lap = new[i, j, k]
        new[i, j, k] = 2.0 * psi[i, j, k] - prev[i, j, k] + c2dt2 * lap
    # Swap time levels.
    for i, j, k in ti.ndrange(nx, ny, nz):
        prev[i, j, k] = psi[i, j, k]
        psi[i, j, k] = new[i, j, k]