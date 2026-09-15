"""Initial condition: a Gaussian pulse of radial displacement."""

from ..pipeline import BaseProcessor, Stage

import taichi as ti

from .features import PsiField, WaveGrid


class SeedPulse(BaseProcessor):
    """
    Seeds psi(0) and psi_prev(-dt) with a Gaussian radial pulse.
    Runs once, at the first PRE_UPDATE step.
    """
    name = "SeedPulse"
    stage = Stage.PRE_UPDATE
    order = 10
    requires = (WaveGrid, PsiField)

    def __init__(self, amplitude: float = 1.0, radius: float = 4.0):
        self.amplitude = amplitude
        self.radius = radius

    def process(self, ctx) -> None:
        if ctx.sim.step > 0:
            return
        grid = ctx.data.require(WaveGrid)
        field = ctx.data.require(PsiField)
        _seed_pulse(
            field.psi, field.psi_prev,
            grid.nx, grid.ny, grid.nz,
            self.amplitude, self.radius,
        )


@ti.kernel
def _seed_pulse(
    psi: ti.template(),
    prev: ti.template(),
    nx: ti.i32, ny: ti.i32, nz: ti.i32,
    amplitude: ti.f32, radius: ti.f32,
):
    cx = nx * 0.5
    cy = ny * 0.5
    cz = nz * 0.5
    inv_r2 = 1.0 / (radius * radius)
    for i, j, k in ti.ndrange(nx, ny, nz):
        d = ti.Vector([i - cx, j - cy, k - cz])
        r = d.norm()
        v = ti.Vector([0.0, 0.0, 0.0])
        if r > 1e-6 and r < 3.0 * radius:
            env = amplitude * ti.exp(-r * r * inv_r2)
            v = env * (d / r)
        psi[i, j, k] = v
        prev[i, j, k] = v