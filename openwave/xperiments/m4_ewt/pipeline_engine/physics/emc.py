"""
EMC density dynamics and wave speed modulation.

B4a instantaneous density:
    rho_norm = 1.0 - beta_rho * |psi|^2
where rho_norm is normalised to the statutory background
(1.0 = vacuum, 0 = fully depleted). The absolute density rho_0 is
not stored in f32 (it is ~3.3e52, beyond f32 range); the normalised
form is used throughout the EMC chain.

Wave speed modulation B5a:
    c^2(rho) = c0^2 * rho_norm
so depleted regions have a lower local wave speed, which is what
makes the soliton self-trapping and gives the c^2(rho) nonlinearity
the plan's Section 1.5 describes.

Both processors run in Stage.UPDATE, before LaplacianVariableCoeff:
    ClearAcceleration   (0)
    UpdateEMCDensity    (5)
    UpdateWaveSpeed     (6)
    LaplacianVarCoeff   (10)
    NonlinearCubic      (15, optional)
    Leapfrog            (20)

Note: no `from __future__ import annotations` here. Taichi kernel
definitions need live type objects.
"""

from ..pipeline import BaseProcessor, Stage
from .features import (
    EMCDensityField,
    PsiLongField,
    WaveGrid,
    WaveSpeedField,
)
from .units import UnitSystem

import taichi as ti


class UpdateEMCDensityProcessor(BaseProcessor):
    """
    B4a instantaneous EMC density:
        rho_norm = 1.0 - beta_rho * |psi|^2

    Writes normalised rho into EMCDensityField.rho. The result is
    immediately consistent with the current psi; the density does not
    evolve on its own (that is B4b/B4c, deferred).

    beta_rho sets the push-out strength, normalised to the statutory
    background: the processor here uses the plan's beta_rho / rho_0 as a
    single dimensionless coefficient, because rho is stored normalised.
    beta_rho = 0 gives uniform statutory density (no push-out);
    beta_rho = 1 gives full depletion at |psi|^2 = 1.

    No floor is applied to rho. The regime beta_rho |psi|^2 > 1 gives
    rho < 0, hence c2_local < 0, which is outside the model's domain of
    validity. A pipeline that reaches this regime is misconfigured;
    the negative c2 is the signal. Clamping would hide the parameter
    error.
    """

    name = "UpdateEMCDensity"
    stage = Stage.UPDATE
    order = 5
    provides = ()

    def __init__(self, field_type=PsiLongField, beta_rho: float = 1.0):
        self.field_type = field_type
        self.beta_rho = float(beta_rho)
        self.requires = (WaveGrid, field_type, EMCDensityField)

    def process(self, ctx) -> None:
        grid = ctx.data.require(WaveGrid)
        psi = ctx.data.require(self.field_type).psi
        rho = ctx.data.require(EMCDensityField).rho
        _update_rho(psi, rho, self.beta_rho, grid.nx, grid.ny, grid.nz)


@ti.kernel
def _update_rho(
    psi: ti.template(),
    rho: ti.template(),
    beta_rho: ti.f32,
    nx: ti.i32,
    ny: ti.i32,
    nz: ti.i32,
):
    for i, j, k in ti.ndrange(nx, ny, nz):
        rho[i, j, k] = 1.0 - beta_rho * psi[i, j, k].norm_sqr()


class UpdateWaveSpeedProcessor(BaseProcessor):
    """
    B5a wave speed modulation:
        c2_local = c0^2 * rho_norm

    Reads EMCDensityField.rho (normalised), writes WaveSpeedField.c2_local.
    c0 comes from the UnitSystem feature. A depleted region (rho_norm < 1)
    has a lower local c^2, which is what makes the nonlinearity
    self-trapping.
    """

    name = "UpdateWaveSpeed"
    stage = Stage.UPDATE
    order = 6
    provides = ()

    def __init__(self):
        self.requires = (WaveGrid, UnitSystem, EMCDensityField, WaveSpeedField)

    def process(self, ctx) -> None:
        grid = ctx.data.require(WaveGrid)
        units = ctx.data.require(UnitSystem)
        rho = ctx.data.require(EMCDensityField).rho
        c2 = ctx.data.require(WaveSpeedField).c2_local
        _update_c2(rho, c2, float(units.c) ** 2, grid.nx, grid.ny, grid.nz)


@ti.kernel
def _update_c2(
    rho: ti.template(),
    c2: ti.template(),
    c0_sq: ti.f32,
    nx: ti.i32,
    ny: ti.i32,
    nz: ti.i32,
):
    for i, j, k in ti.ndrange(nx, ny, nz):
        c2[i, j, k] = c0_sq * rho[i, j, k]