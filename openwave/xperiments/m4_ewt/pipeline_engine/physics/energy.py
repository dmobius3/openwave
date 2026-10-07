"""
Energy budget (plan Section 6, item 1.8).

EnergyBudgetUpdate integrates the field's kinetic, gradient and
deformation energy after each step and stores them in EnergyBudget.

Definitions (plan Section 5.2):
    E_kin    = integral (1/2) |dpsi/dt|^2 dV
    E_grad   = integral (1/2) c^2(rho) |grad psi|^2 dV
    E_deform = integral (1/2) kappa (rho - rho_0)^2 dV
    E_total  = E_kin + E_grad + E_deform

Discretisation:
    dpsi/dt is (psi - psi_prev)/dt, at each voxel. This is the
    symplectic mid-step velocity the leapfrog uses, not a centred
    difference over three time levels.

    |grad psi|^2 is the average of the two one-sided gradients per
    axis. This is the naive gradient form, not the flux-form half-grid
    average the flux-form Laplacian uses. The two differ by O(dx^2)
    at constant c^2, and c^2 is constant in the conservation test (no
    WaveSpeedField registered). Where c^2 varies the flux form is the
    rigorous one, but writing its matching half-grid energy kernel is
    a larger change than plan 1.8 requires; the constant-c^2 case
    covers the HO test the plan names, and the variable-c^2 case is
    handled by a second kernel that reads c^2 at the voxel centre.

    The kinetic term is required. Without it, a standing wave's
    gradient-only integral oscillates at 2*omega and the conservation
    check fails by construction (plan Section 5.2).

Accumulators are f64 to avoid f32 summation loss on larger grids.
The values are read back to Python floats once per step.

kappa is a per-instance configuration on the processor, not a field
on the feature. It is a measurement parameter, not a result.
"""

import taichi as ti

from ..pipeline import BaseProcessor, Stage
from .features import (
    EMCDensityField,
    EnergyBudget,
    PsiLongField,
    WaveGrid,
    WaveSpeedField,
)
from .units import UnitSystem


@ti.func
def _voxel_energy(
    psi: ti.template(),
    psi_prev: ti.template(),
    rho: ti.template(),
    c2_here: ti.f32,
    kappa: ti.f32,
    inv_dt: ti.f32,
    inv_dx: ti.f32,
    dv: ti.f32,
    i: ti.i32,
    j: ti.i32,
    k: ti.i32,
):
    """
    Per-voxel kinetic, gradient and deformation energy.
    Returns (e_kin, e_grad, e_deform).
    """
    v = (psi[i, j, k] - psi_prev[i, j, k]) * inv_dt
    e_kin = 0.5 * v.norm_sqr() * dv

    gx = 0.5 * (
        (psi[i + 1, j, k] - psi[i, j, k])
        + (psi[i, j, k] - psi[i - 1, j, k])
    ) * inv_dx
    gy = 0.5 * (
        (psi[i, j + 1, k] - psi[i, j, k])
        + (psi[i, j, k] - psi[i, j - 1, k])
    ) * inv_dx
    gz = 0.5 * (
        (psi[i, j, k + 1] - psi[i, j, k])
        + (psi[i, j, k] - psi[i, j, k - 1])
    ) * inv_dx
    e_grad = (
        0.5 * c2_here
        * (gx.norm_sqr() + gy.norm_sqr() + gz.norm_sqr())
        * dv
    )

    drho = rho[i, j, k] - 1.0
    e_deform = 0.5 * kappa * drho * drho * dv
    return e_kin, e_grad, e_deform


@ti.kernel
def _integrate_energy_const_c2(
    psi: ti.template(),
    psi_prev: ti.template(),
    rho: ti.template(),
    kappa: ti.f32,
    c2_const: ti.f32,
    inv_dt: ti.f32,
    inv_dx: ti.f32,
    dv: ti.f32,
    nx: ti.i32,
    ny: ti.i32,
    nz: ti.i32,
    out_kin: ti.template(),
    out_grad: ti.template(),
    out_deform: ti.template(),
):
    out_kin[None] = 0.0
    out_grad[None] = 0.0
    out_deform[None] = 0.0
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        ek, eg, ed = _voxel_energy(
            psi, psi_prev, rho, c2_const, kappa,
            inv_dt, inv_dx, dv, i, j, k,
        )
        out_kin[None] += ek
        out_grad[None] += eg
        out_deform[None] += ed


@ti.kernel
def _integrate_energy_var_c2(
    psi: ti.template(),
    psi_prev: ti.template(),
    rho: ti.template(),
    c2_field: ti.template(),
    kappa: ti.f32,
    inv_dt: ti.f32,
    inv_dx: ti.f32,
    dv: ti.f32,
    nx: ti.i32,
    ny: ti.i32,
    nz: ti.i32,
    out_kin: ti.template(),
    out_grad: ti.template(),
    out_deform: ti.template(),
):
    out_kin[None] = 0.0
    out_grad[None] = 0.0
    out_deform[None] = 0.0
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        ek, eg, ed = _voxel_energy(
            psi, psi_prev, rho, c2_field[i, j, k], kappa,
            inv_dt, inv_dx, dv, i, j, k,
        )
        out_kin[None] += ek
        out_grad[None] += eg
        out_deform[None] += ed


class EnergyBudgetUpdate(BaseProcessor):
    """
    Integrate the field's kinetic, gradient and deformation energy
    into EnergyBudget. Stage.MEASURE, order 5.

    Provides EnergyBudget (created in setup). Requires WaveGrid, the
    field named by field_type, EMCDensityField and UnitSystem. With
    use_variable_c2=True, requires WaveSpeedField and reads c^2 from
    it; otherwise reads c^2 = units.c ** 2 as a constant.

    kappa is the EMC deformation stiffness, a per-instance
    configuration. Default 1.0. The plan's "supplied by the unit
    system" is not done yet: UnitSystem has no kappa field, and adding
    one is a larger change.

    Interior-only, matching the Laplacian and Leapfrog scopes. The
    excluded shell is a boundary artefact, not a physical energy.

    Does not read or write any other feature. Does not enforce a
    conservation law; it reports, it does not police.
    """

    name = "EnergyBudgetUpdate"
    stage = Stage.MEASURE
    order = 5
    provides = (EnergyBudget,)

    def __init__(self, field_type=PsiLongField, use_variable_c2=False,
                 kappa=1.0):
        self.field_type = field_type
        self.use_variable_c2 = use_variable_c2
        self.kappa = float(kappa)
        base = (WaveGrid, field_type, EMCDensityField, UnitSystem)
        if use_variable_c2:
            self.requires = base + (WaveSpeedField,)
        else:
            self.requires = base

    def setup(self, ctx):
        ctx.data.set(EnergyBudget())
        self._out_kin = ti.field(dtype=ti.f64, shape=())
        self._out_grad = ti.field(dtype=ti.f64, shape=())
        self._out_deform = ti.field(dtype=ti.f64, shape=())

    def process(self, ctx):
        grid = ctx.data.require(WaveGrid)
        units = ctx.data.require(UnitSystem)
        field = ctx.data.require(self.field_type)
        rho = ctx.data.require(EMCDensityField).rho
        budget = ctx.data.require(EnergyBudget)

        dt = ctx.sim.dt
        if dt <= 0.0:
            raise ValueError(
                f"EnergyBudgetUpdate: dt must be > 0, got {dt}"
            )

        inv_dt = 1.0 / dt
        inv_dx = 1.0 / grid.dx
        dv = grid.dx ** 3

        if self.use_variable_c2:
            c2_field = ctx.data.require(WaveSpeedField).c2_local
            _integrate_energy_var_c2(
                field.psi, field.psi_prev, rho, c2_field,
                self.kappa, inv_dt, inv_dx, dv,
                grid.nx, grid.ny, grid.nz,
                self._out_kin, self._out_grad, self._out_deform,
            )
        else:
            c2_const = float(units.c) ** 2
            _integrate_energy_const_c2(
                field.psi, field.psi_prev, rho,
                self.kappa, c2_const, inv_dt, inv_dx, dv,
                grid.nx, grid.ny, grid.nz,
                self._out_kin, self._out_grad, self._out_deform,
            )

        budget.E_kin = float(self._out_kin[None])
        budget.E_grad = float(self._out_grad[None])
        budget.E_deform = float(self._out_deform[None])

        if ctx.sim.step == 0:
            budget.E_prev = budget.E_total
            budget.dE_dt = 0.0
        else:
            budget.dE_dt = (budget.E_total - budget.E_prev) / dt
            budget.E_prev = budget.E_total