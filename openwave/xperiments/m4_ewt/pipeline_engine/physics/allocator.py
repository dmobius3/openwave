"""Wave field allocator. Lifecycle-only processor."""

from __future__ import annotations

from ..pipeline import BaseProcessor

import taichi as ti

from .features import PsiField, WaveGrid, WaveStats


class AllocateWaveField(BaseProcessor):
    """
    Allocates the wave field, the grid feature, and the stats holder.
    Runs once at setup. Provides: WaveGrid, PsiField, WaveStats.
    """

    name = "AllocateWaveField"
    stage = None
    order = 10
    provides = (WaveGrid, PsiField, WaveStats)

    def __init__(self, nx: int, ny: int, nz: int, dx: float, c: float):
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.dx = dx
        self.c = c

    def setup(self, ctx) -> None:
        grid = WaveGrid(nx=self.nx, ny=self.ny, nz=self.nz, dx=self.dx, c=self.c)
        shape = (self.nx, self.ny, self.nz)
        ctx.data.set(grid)
        ctx.data.set(
            PsiField(
                psi=ti.Vector.field(3, dtype=ti.f32, shape=shape),
                psi_prev=ti.Vector.field(3, dtype=ti.f32, shape=shape),
                psi_new=ti.Vector.field(3, dtype=ti.f32, shape=shape),
            )
        )
        ctx.data.set(WaveStats())
