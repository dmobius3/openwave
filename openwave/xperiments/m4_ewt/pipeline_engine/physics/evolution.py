"""
Wave evolution: spatial operator then leapfrog.

Contract: psi_new is the acceleration accumulator between UPDATE processors.
Each UPDATE processor adds its contribution to psi_new:
    Laplacian   -> +c^2 * Laplacian(psi)
    Nonlinearity-> -gamma * |psi|^2 * psi   (additive, runs after Laplacian)
    Leapfrog    -> integrates: psi_new = 2*psi - psi_prev + dt^2 * accel
"""

from ..pipeline import BaseProcessor, Stage

import taichi as ti

from .features import PsiField, WaveGrid


class LaplacianProcessor(BaseProcessor):
    """
    Writes c^2 * Laplacian(psi) into psi_new (acceleration contribution
    from the spatial operator). Overwrites; must run before any additive
    processor in the UPDATE stage.
    """
    name = "Laplacian"
    stage = Stage.UPDATE
    order = 10
    requires = (WaveGrid, PsiField)

    def process(self, ctx) -> None:
        grid = ctx.data.require(WaveGrid)
        field = ctx.data.require(PsiField)
        _laplacian(field.psi, field.psi_new,
                   grid.nx, grid.ny, grid.nz, grid.dx, grid.c)


@ti.kernel
def _laplacian(
    psi: ti.template(),
    out: ti.template(),
    nx: ti.i32, ny: ti.i32, nz: ti.i32,
    dx: ti.f32,
    c: ti.f32,
):
    inv_dx2 = 1.0 / (dx * dx)
    c2 = c * c
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        face_sum = (
            psi[i + 1, j, k] + psi[i - 1, j, k]
            + psi[i, j + 1, k] + psi[i, j - 1, k]
            + psi[i, j, k + 1] + psi[i, j, k - 1]
        )
        out[i, j, k] = c2 * (face_sum - 6.0 * psi[i, j, k]) * inv_dx2


class LeapfrogProcessor(BaseProcessor):
    """
    Integrates: psi_new = 2*psi - psi_prev + dt^2 * psi_new
    where psi_new holds the accumulated acceleration (spatial + forces).
    Then swaps time levels: prev <- psi, psi <- new.
    """
    name = "Leapfrog"
    stage = Stage.UPDATE
    order = 20
    requires = (WaveGrid, PsiField)

    def process(self, ctx) -> None:
        grid = ctx.data.require(WaveGrid)
        field = ctx.data.require(PsiField)
        dt2 = ctx.sim.dt * ctx.sim.dt
        _leapfrog(field.psi, field.psi_prev, field.psi_new,
                  grid.nx, grid.ny, grid.nz, dt2)


@ti.kernel
def _leapfrog(
    psi: ti.template(),
    prev: ti.template(),
    new: ti.template(),
    nx: ti.i32, ny: ti.i32, nz: ti.i32,
    dt2: ti.f32,
):
    # Read accumulated acceleration, produce psi(t+dt) into the same buffer.
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        accel = new[i, j, k]
        new[i, j, k] = 2.0 * psi[i, j, k] - prev[i, j, k] + dt2 * accel
    # Swap time levels.
    for i, j, k in ti.ndrange(nx, ny, nz):
        prev[i, j, k] = psi[i, j, k]
        psi[i, j, k] = new[i, j, k]