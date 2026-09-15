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
            field.psi,
            field.psi_prev,
            grid.nx,
            grid.ny,
            grid.nz,
            self.amplitude,
            self.radius,
        )


@ti.kernel
def _seed_pulse(
    psi: ti.template(),
    prev: ti.template(),
    nx: ti.i32,
    ny: ti.i32,
    nz: ti.i32,
    amplitude: ti.f32,
    radius: ti.f32,
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


from .wc_types import WCState


class SeedMultiCenter(BaseProcessor):
    """
    Seeds the wave field with a superposition of localized packets,
    one per active wave center in WCState.

    Each packet is a Gaussian radial displacement:
        psi += A * cos(phase) * exp(-r^2 / radius^2) * r_hat
    summed over all active centers. Released from rest.

    Runs once, at the first PRE_UPDATE step.
    """

    name = "SeedMultiCenter"
    stage = Stage.PRE_UPDATE
    order = 10
    requires = (WaveGrid, PsiField, WCState)

    def __init__(self, radius: float = 4.0):
        self.radius = radius

    def process(self, ctx) -> None:
        if ctx.sim.step > 0:
            return
        grid = ctx.data.require(WaveGrid)
        field = ctx.data.require(PsiField)
        wc_state = ctx.data.require(WCState)

        _zero_field(field.psi, field.psi_prev, grid.nx, grid.ny, grid.nz)
        for wc in wc_state.centers:
            if not wc.active:
                continue
            _add_pulse(
                field.psi,
                field.psi_prev,
                grid.nx,
                grid.ny,
                grid.nz,
                wc.x,
                wc.y,
                wc.z,
                wc.amplitude,
                self.radius,
                wc.phase,
            )


@ti.kernel
def _zero_field(
    psi: ti.template(),
    prev: ti.template(),
    nx: ti.i32,
    ny: ti.i32,
    nz: ti.i32,
):
    z = ti.Vector([0.0, 0.0, 0.0])
    for i, j, k in ti.ndrange(nx, ny, nz):
        psi[i, j, k] = z
        prev[i, j, k] = z


@ti.kernel
def _add_pulse(
    psi: ti.template(),
    prev: ti.template(),
    nx: ti.i32,
    ny: ti.i32,
    nz: ti.i32,
    cx: ti.f32,
    cy: ti.f32,
    cz: ti.f32,
    amplitude: ti.f32,
    radius: ti.f32,
    phase: ti.f32,
):
    inv_r2 = 1.0 / (radius * radius)
    phase_cos = ti.cos(phase)
    for i, j, k in ti.ndrange(nx, ny, nz):
        d = ti.Vector([i - cx, j - cy, k - cz])
        r = d.norm()
        if r > 1e-6 and r < 3.0 * radius:
            env = amplitude * phase_cos * ti.exp(-r * r * inv_r2)
            v = env * (d / r)
            psi[i, j, k] += v
            prev[i, j, k] += v
