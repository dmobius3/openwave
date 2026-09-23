"""Amplitude and mass trackers. Write into WaveStats."""

from ..pipeline import BaseProcessor, Stage

import taichi as ti

from .features import PsiLongField, WaveGrid, WaveStats


class AmplitudeTracker(BaseProcessor):
    """Max |psi|^2 and total mass (sum of |psi|^2) into WaveStats."""

    name = "AmplitudeTracker"
    stage = Stage.MEASURE
    order = 10
    requires = (WaveGrid, PsiLongField, WaveStats)

    def __init__(self, every: int = 10):
        self.every = max(1, every)

    def setup(self, ctx) -> None:
        self._max = ti.field(dtype=ti.f32, shape=())
        self._mass = ti.field(dtype=ti.f32, shape=())

    def process(self, ctx) -> None:
        if ctx.sim.step % self.every != 0:
            return
        grid = ctx.data.require(WaveGrid)
        field = ctx.data.require(PsiLongField)
        stats = ctx.data.require(WaveStats)
        _reduce(field.psi, self._max, self._mass, grid.nx, grid.ny, grid.nz)
        stats.amp_max = float(self._max[None])
        stats.mass = float(self._mass[None])


@ti.kernel
def _reduce(
    psi: ti.template(),
    out_max: ti.template(),
    out_mass: ti.template(),
    nx: ti.i32,
    ny: ti.i32,
    nz: ti.i32,
):
    out_max[None] = 0.0
    out_mass[None] = 0.0
    for i, j, k in ti.ndrange(nx, ny, nz):
        s = psi[i, j, k].norm_sqr()
        ti.atomic_max(out_max[None], s)
        ti.atomic_add(out_mass[None], s)
