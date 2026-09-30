"""
Wave evolution: spatial operator then leapfrog.

Contract: psi_new is the acceleration accumulator between UPDATE processors.
Each UPDATE processor adds its contribution to psi_new:
    Laplacian   -> +c^2 * Laplacian(psi)
    Nonlinearity-> -gamma * |psi|^2 * psi   (additive, runs after Laplacian)
    Leapfrog    -> integrates: psi_new = 2*psi - psi_prev + dt^2 * accel

Note: no `from __future__ import annotations` here. Taichi kernel
definitions need live type objects, not stringified annotations.
"""

from ..pipeline import BaseProcessor, Stage

import taichi as ti

from .features import PsiLongField, WaveGrid


class LaplacianProcessor(BaseProcessor):
    """
    c^2 * Laplacian(psi) into psi_new. Overwrites; must run before any
    additive processor in the UPDATE stage.

    field_type selects which triple-buffer field to operate on. Default
    is PsiLongField, preserving the behaviour of the pre-parameterised
    version. Pass PsiBaseField (or any PsiTripleBuffer subclass) to run
    the same operator on a different field.
    """

    name = "Laplacian"
    stage = Stage.UPDATE
    order = 10
    provides = ()

    def __init__(self, field_type=PsiLongField):
        self.field_type = field_type
        self.requires = (WaveGrid, field_type)

    def process(self, ctx) -> None:
        grid = ctx.data.require(WaveGrid)
        field = ctx.data.require(self.field_type)
        _laplacian(field.psi, field.psi_new, grid.nx, grid.ny, grid.nz, grid.dx, grid.c)


@ti.kernel
def _laplacian(
    psi: ti.template(),
    out: ti.template(),
    nx: ti.i32,
    ny: ti.i32,
    nz: ti.i32,
    dx: ti.f32,
    c: ti.f32,
):
    inv_dx2 = 1.0 / (dx * dx)
    c2 = c * c
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        face_sum = (
            psi[i + 1, j, k]
            + psi[i - 1, j, k]
            + psi[i, j + 1, k]
            + psi[i, j - 1, k]
            + psi[i, j, k + 1]
            + psi[i, j, k - 1]
        )
        out[i, j, k] = c2 * (face_sum - 6.0 * psi[i, j, k]) * inv_dx2


class LeapfrogProcessor(BaseProcessor):
    """
    psi_new = 2*psi - psi_prev + dt^2 * psi_new, then swaps time levels.

    field_type selects which triple-buffer field to integrate. Default
    is PsiLongField. Two fields evolved independently means two
    LeapfrogProcessor instances with different field_type, not one
    instance handling both.
    """

    name = "Leapfrog"
    stage = Stage.UPDATE
    order = 20
    provides = ()

    def __init__(self, field_type=PsiLongField):
        self.field_type = field_type
        self.requires = (WaveGrid, field_type)

    def process(self, ctx) -> None:
        grid = ctx.data.require(WaveGrid)
        field = ctx.data.require(self.field_type)
        dt2 = ctx.sim.dt * ctx.sim.dt
        _leapfrog(field.psi, field.psi_prev, field.psi_new, grid.nx, grid.ny, grid.nz, dt2)


@ti.kernel
def _leapfrog(
    psi: ti.template(),
    prev: ti.template(),
    new: ti.template(),
    nx: ti.i32,
    ny: ti.i32,
    nz: ti.i32,
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