"""
Wave evolution: accumulator clear, spatial operator, leapfrog.

Contract:
    Stage.UPDATE processors write to psi_new (the acceleration
    accumulator) and must run in this order:
        ClearAcceleration  (order 0)   psi_new = 0 on the full grid
        Laplacian          (order 10)  psi_new += c^2 * Laplacian(psi)
        NonlinearCubic     (order 15)  psi_new += -gamma |psi|^2 psi
        Leapfrog           (order 20)  psi_new = 2 psi - psi_prev + dt^2 psi_new
    Every processor after ClearAcceleration uses +=. The order among the
    additive processors does not matter for correctness (addition is
    commutative); only the position of Clear first and Leapfrog last.

Note: no `from __future__ import annotations` here. Taichi kernel
definitions need live type objects, not stringified annotations.
"""

from ..pipeline import BaseProcessor, Stage

import taichi as ti

from .features import PsiLongField, WaveGrid, WaveSpeedField
from .units import UnitSystem

# ======================================================================
# Accumulator clear
# ======================================================================


class ClearAccelerationProcessor(BaseProcessor):
    """
    Zeroes psi_new on the full grid at the start of the UPDATE stage.

    Without this, psi_new would retain the acceleration from the previous
    step, and any pipeline that omits the Laplacian (a pure nonlinear
    test, a test with only an external force) would silently integrate
    stale acceleration. With it, every UPDATE processor that follows may
    use += safely.

    field_type selects which triple-buffer field's accumulator to clear.
    """

    name = "ClearAcceleration"
    stage = Stage.UPDATE
    order = 0
    provides = ()

    def __init__(self, field_type=PsiLongField):
        self.field_type = field_type
        self.requires = (WaveGrid, field_type)

    def process(self, ctx) -> None:
        grid = ctx.data.require(WaveGrid)
        field = ctx.data.require(self.field_type)
        _clear_accel(field.psi_new, grid.nx, grid.ny, grid.nz)


@ti.kernel
def _clear_accel(
    new: ti.template(),
    nx: ti.i32,
    ny: ti.i32,
    nz: ti.i32,
):
    z = ti.Vector([0.0, 0.0, 0.0])
    for i, j, k in ti.ndrange(nx, ny, nz):
        new[i, j, k] = z


# ======================================================================
# Laplacian
# ======================================================================


class LaplacianProcessor(BaseProcessor):
    """
    Adds c^2 * Laplacian(psi) to psi_new. Uses += so it composes with
    other additive processors; requires ClearAccelerationProcessor to
    have zeroed psi_new first (order 0 before order 10).

    The wave speed c comes from the UnitSystem feature, not from the
    grid. field_type selects which triple-buffer field to operate on.
    Default is PsiLongField.
    """

    name = "Laplacian"
    stage = Stage.UPDATE
    order = 10
    provides = ()

    def __init__(self, field_type=PsiLongField):
        self.field_type = field_type
        self.requires = (WaveGrid, field_type, UnitSystem)

    def process(self, ctx) -> None:
        grid = ctx.data.require(WaveGrid)
        units = ctx.data.require(UnitSystem)
        field = ctx.data.require(self.field_type)
        _laplacian(
            field.psi, field.psi_new,
            grid.nx, grid.ny, grid.nz, grid.dx,
            float(units.c),
        )


# ======================================================================
# Leapfrog
# ======================================================================

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
        out[i, j, k] += c2 * (face_sum - 6.0 * psi[i, j, k]) * inv_dx2

class LeapfrogProcessor(BaseProcessor):
    """
    psi_new = 2*psi - psi_prev + dt^2 * psi_new, then swaps time levels
    on the interior only. Boundary voxels of psi and psi_prev are left
    to the boundary processor in POST_UPDATE.

    field_type selects which triple-buffer field to integrate. Default
    is PsiLongField. Two fields evolved independently means two
    LeapfrogProcessor instances with different field_type.
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
    # First pass: integrate the acceleration into new, interior only.
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        accel = new[i, j, k]
        new[i, j, k] = 2.0 * psi[i, j, k] - prev[i, j, k] + dt2 * accel
    # Second pass: swap time levels, interior only. Boundary voxels are
    # owned by the boundary processor in POST_UPDATE; the swap must not
    # touch them, or it would overwrite what that processor wrote.
    #
    # Note on performance: this is a physical copy, not a pointer
    # rotation. A pointer swap at the Python level would invalidate
    # references held through the FeatureBag and is deferred to a
    # separate contract change.
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        prev[i, j, k] = psi[i, j, k]
        psi[i, j, k] = new[i, j, k]

class LaplacianVariableCoeffProcessor(BaseProcessor):
    """
    Adds div(c^2(rho) grad psi) to psi_new, using the conservative flux
    form on half-grids:
        c^2_{i+1/2} = (c^2_i + c^2_{i+1}) / 2
        flux_x = c^2_{i+1/2} (psi_{i+1} - psi_i) / dx
        div(c^2 grad psi) = (flux_x[i+1/2] - flux_x[i-1/2]) / dx + ...
    Uses += so it composes with other additive processors; requires
    ClearAccelerationProcessor to have zeroed psi_new first (order 0
    before order 10).

    The naive form c^2_i * laplacian(psi) is not self-adjoint when c^2
    varies, so it breaks the staggered leapfrog invariant. The flux
    form is the discrete variation of
        E_grad = sum (1/2) c^2_{i+1/2} ((psi_{i+1} - psi_i)/dx)^2 dx
    and conserves the invariant. Plan item 1.23 requires the flux form;
    the naive form is a bug.

    Requires WaveSpeedField to have been populated with c^2 values.
    """

    name = "LaplacianVariableCoeff"
    stage = Stage.UPDATE
    order = 10
    provides = ()

    def __init__(self, field_type=PsiLongField):
        self.field_type = field_type
        self.requires = (WaveGrid, field_type, WaveSpeedField)

    def process(self, ctx) -> None:
        grid = ctx.data.require(WaveGrid)
        field = ctx.data.require(self.field_type)
        c2 = ctx.data.require(WaveSpeedField).c2_local
        _laplacian_var_coeff(
            field.psi, c2, field.psi_new,
            grid.nx, grid.ny, grid.nz, grid.dx,
        )


@ti.kernel
def _laplacian_var_coeff(
    psi: ti.template(),
    c2: ti.template(),
    out: ti.template(),
    nx: ti.i32,
    ny: ti.i32,
    nz: ti.i32,
    dx: ti.f32,
):
    inv_dx2 = 1.0 / (dx * dx)
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        c_ii = c2[i, j, k]
        c2_xp = 0.5 * (c_ii + c2[i + 1, j, k])
        c2_xm = 0.5 * (c_ii + c2[i - 1, j, k])
        c2_yp = 0.5 * (c_ii + c2[i, j + 1, k])
        c2_ym = 0.5 * (c_ii + c2[i, j - 1, k])
        c2_zp = 0.5 * (c_ii + c2[i, j, k + 1])
        c2_zm = 0.5 * (c_ii + c2[i, j, k - 1])

        flux_x = (
            c2_xp * (psi[i + 1, j, k] - psi[i, j, k])
            - c2_xm * (psi[i, j, k] - psi[i - 1, j, k])
        )
        flux_y = (
            c2_yp * (psi[i, j + 1, k] - psi[i, j, k])
            - c2_ym * (psi[i, j, k] - psi[i, j - 1, k])
        )
        flux_z = (
            c2_zp * (psi[i, j, k + 1] - psi[i, j, k])
            - c2_zm * (psi[i, j, k] - psi[i, j, k - 1])
        )
        out[i, j, k] += (flux_x + flux_y + flux_z) * inv_dx2