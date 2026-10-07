"""
Tests for EnergyBudget and EnergyBudgetUpdate (plan item 1.8).

A standing wave in a 3D box, seeded at rest, is integrated with the
production Clear / Laplacian / Leapfrog chain. The kinetic + gradient
+ deformation total is conserved to O(dt^2), per plan Section 5.2:
the formal E_soliton is not the leapfrog's machine-precision
invariant (that is the staggered expression in _test_variable_coeff).

The negative control (gradient-only) shows that E_kin is required:
without it, the integral oscillates at 2*omega and the drift is far
larger than the tolerance the positive test uses.

Run from the project root:
    python -m openwave.xperiments.m4_ewt.pipeline_engine.physics._test_energy_budget

Note: no `from __future__ import annotations` here. Taichi kernel
definitions need live type objects.
"""

import sys
import traceback

import numpy as np
import taichi as ti

from .features import (
    EMCDensityField,
    EnergyBudget,
    PsiLongField,
    WaveGrid,
)
from ..pipeline import (
    BaseProcessor,
    ErrorPolicy,
    Pipeline,
    PipelineError,
    Stage,
)
from .allocator import AllocateWaveField
from .emc import UpdateEMCDensityProcessor, UpdateWaveSpeedProcessor
from .energy import EnergyBudgetUpdate
from .evolution import (
    ClearAccelerationProcessor,
    LaplacianProcessor,
    LeapfrogProcessor,
)
from .units import NaturalUnitSystem, UnitSystem


_TI_INITIALIZED = False


def _ti_init():
    global _TI_INITIALIZED
    if not _TI_INITIALIZED:
        ti.init(arch=ti.cpu, log_level=ti.ERROR)
        _TI_INITIALIZED = True


# ======================================================================
# Fixtures
# ======================================================================


class _SeedStandingWave(BaseProcessor):
    """
    Seed PsiLong with a 3D standing wave, rest start (psi_prev = psi).
    Modes are mode * pi / (n - 1) per axis, zero on the outer shell.
    """

    name = "_SeedStandingWave"
    stage = Stage.PRE_UPDATE
    order = 5
    requires = (WaveGrid, PsiLongField)

    def __init__(self, amp=1.0, mode=1):
        self.amp = float(amp)
        self.mode = int(mode)

    def process(self, ctx):
        if ctx.sim.step > 0:
            return
        grid = ctx.data.require(WaveGrid)
        field = ctx.data.require(PsiLongField)
        _seed_standing_wave(
            field.psi, field.psi_prev,
            grid.nx, grid.ny, grid.nz,
            self.amp, self.mode,
        )


@ti.kernel
def _seed_standing_wave(
    psi: ti.template(),
    prev: ti.template(),
    nx: ti.i32,
    ny: ti.i32,
    nz: ti.i32,
    amp: ti.f32,
    mode: ti.i32,
):
    kx = ti.cast(mode, ti.f32) * ti.math.pi / ti.cast(nx - 1, ti.f32)
    ky = ti.cast(mode, ti.f32) * ti.math.pi / ti.cast(ny - 1, ti.f32)
    kz = ti.cast(mode, ti.f32) * ti.math.pi / ti.cast(nz - 1, ti.f32)
    for i, j, k in ti.ndrange(nx, ny, nz):
        s = (
            ti.sin(kx * ti.cast(i, ti.f32))
            * ti.sin(ky * ti.cast(j, ti.f32))
            * ti.sin(kz * ti.cast(k, ti.f32))
        )
        v = ti.Vector([amp * s, 0.0, 0.0])
        psi[i, j, k] = v
        prev[i, j, k] = v


class _SeedRhoUniform(BaseProcessor):
    """Set rho = 1.0 everywhere at step 0."""

    name = "_SeedRhoUniform"
    stage = Stage.PRE_UPDATE
    order = 6
    requires = (EMCDensityField,)

    def process(self, ctx):
        if ctx.sim.step > 0:
            return
        ctx.data.require(EMCDensityField).rho.fill(1.0)


class _SeedRhoDeficit(BaseProcessor):
    """
    Set rho = 1 - 0.5 * exp(-r^2/8) at step 0. Sharp deficit at centre.
    Used only by test_deficit_contributes_to_deform.
    """

    name = "_SeedRhoDeficit"
    stage = Stage.PRE_UPDATE
    order = 6
    requires = (EMCDensityField, WaveGrid)

    def process(self, ctx):
        if ctx.sim.step > 0:
            return
        grid = ctx.data.require(WaveGrid)
        ctx.data.require(EMCDensityField).rho.from_numpy(
            _rho_deficit_array(grid.nx, grid.ny, grid.nz)
        )


def _rho_deficit_array(nx, ny, nz):
    cx, cy, cz = (nx - 1) / 2.0, (ny - 1) / 2.0, (nz - 1) / 2.0
    i, j, k = np.meshgrid(
        np.arange(nx), np.arange(ny), np.arange(nz), indexing="ij")
    r2 = (i - cx) ** 2 + (j - cy) ** 2 + (k - cz) ** 2
    return (1.0 - 0.5 * np.exp(-r2 / 8.0)).astype(np.float32)


def _run(grid_n=16, dt=0.05, max_steps=100, amp=1.0, mode=1,
         rho_def=False, kappa=1.0):
    """
    Minimal energy-conservation pipeline: allocate, seed standing wave,
    seed rho (uniform or deficit), Clear, Laplacian, Leapfrog,
    EnergyBudgetUpdate.
    """
    class P(Pipeline):
        def __init__(self):
            super().__init__(
                error_policy=ErrorPolicy.FAIL_FAST,
                external_provides=(UnitSystem,),
            )
            self.add(AllocateWaveField(
                nx=grid_n, ny=grid_n, nz=grid_n, dx=1.0))
            self.add(_SeedStandingWave(amp=amp, mode=mode))
            if rho_def:
                self.add(_SeedRhoDeficit())
            else:
                self.add(_SeedRhoUniform())
            self.add(ClearAccelerationProcessor(field_type=PsiLongField))
            self.add(LaplacianProcessor(field_type=PsiLongField))
            self.add(LeapfrogProcessor(field_type=PsiLongField))
            self.add(EnergyBudgetUpdate(
                field_type=PsiLongField, kappa=kappa))

    from ..runner import Runner
    from ..sinks import InMemorySink

    return Runner({"session": InMemorySink()}).run(
        P(), name="energy_budget_test", params={}, dt=dt,
        max_steps=max_steps,
        initial_features=[NaturalUnitSystem(grid_voxels_per_lambda=20,
                                            cfl_safety=0.5)],
    )


# ======================================================================
# Positive tests
# ======================================================================


def test_uniform_rho_gives_zero_deform():
    """
    rho = 1.0 everywhere -> (rho - 1.0) = 0 everywhere -> E_deform = 0.

    Mutation caught: drho computed from |rho| rather than (rho - 1),
    or the deformation term dropped entirely.
    """
    _ti_init()
    ctx = _run(grid_n=12, dt=0.05, max_steps=3)
    assert ctx.diag.errors == [], ctx.diag.errors

    b = ctx.data.require(EnergyBudget)
    assert b.E_deform == 0.0, b.E_deform
    assert b.E_kin > 0.0, b.E_kin
    assert b.E_grad > 0.0, b.E_grad


def test_deficit_contributes_to_deform():
    """
    rho < 1 in the centre -> E_deform > 0, and its value matches an
    independent numpy integral over the interior.

    Mutation caught: E_deform reading a constant instead of the field,
    or the accumulation loop touching the outer shell.
    """
    _ti_init()
    n = 12
    ctx = _run(grid_n=n, dt=0.05, max_steps=2, rho_def=True, kappa=1.0)
    assert ctx.diag.errors == [], ctx.diag.errors

    b = ctx.data.require(EnergyBudget)
    rho_arr = ctx.data.require(EMCDensityField).rho.to_numpy()
    interior = (slice(1, -1),) * 3
    drho = rho_arr[interior] - 1.0
    expected = 0.5 * (drho * drho).sum()  # dx = 1, dv = 1, kappa = 1
    assert abs(b.E_deform - expected) < 1e-3, (b.E_deform, expected)


def test_standing_wave_conserves_total_energy():
    """
    A 3D standing wave at rest in a Dirichlet box, integrated with the
    production Clear / Laplacian / Leapfrog, conserves E_total to
    O(dt^2). Plan Section 5.2 states the bar for the formal E_soliton
    is O(dt^2) convergence, not machine precision.

    The measure is the fraction of E_total moved in one step:
        |dE_dt * dt| / |E_total|

    Mutation caught: the kinetic term dropped from E_total (drift
    rises above the tolerance; see test_gradient_only_energy_oscillates
    for the same failure mode made explicit).
    """
    _ti_init()
    ctx = _run(grid_n=16, dt=0.05, max_steps=50)
    assert ctx.diag.errors == [], ctx.diag.errors

    b = ctx.data.require(EnergyBudget)
    assert b.E_total > 0.0, b.E_total
    step_fraction = abs(b.dE_dt * 0.05) / abs(b.E_total)
    assert step_fraction < 1e-2, (
        f"step fraction drift {step_fraction} vs 1e-2"
    )


def test_dE_dt_zero_on_static_state():
    """
    A field with psi == psi_prev and no gradient has E_total = 0 on
    both steps, so dE_dt = 0 after the first step.

    Mutation caught: dE_dt computed with a sign flip (would not matter
    here, since 0 - 0 = 0), or E_prev updated before the difference
    instead of after (also would not matter here). This is a smoke
    check that the reporting path runs without error; the meaningful
    sign and order checks are in test_standing_wave_conserves_total_energy.
    """
    _ti_init()

    class _SeedZero(BaseProcessor):
        name = "_SeedZero"
        stage = Stage.PRE_UPDATE
        order = 5
        requires = (PsiLongField,)

        def process(self, ctx):
            if ctx.sim.step > 0:
                return
            f = ctx.data.require(PsiLongField)
            f.psi.fill(0.0)
            f.psi_prev.fill(0.0)

    class P(Pipeline):
        def __init__(self):
            super().__init__(
                error_policy=ErrorPolicy.FAIL_FAST,
                external_provides=(UnitSystem,),
            )
            self.add(AllocateWaveField(nx=8, ny=8, nz=8, dx=1.0))
            self.add(_SeedZero())
            self.add(_SeedRhoUniform())
            self.add(EnergyBudgetUpdate(field_type=PsiLongField))

    from ..runner import Runner
    from ..sinks import InMemorySink
    ctx = Runner({"session": InMemorySink()}).run(
        P(), name="static_test", params={}, dt=0.05, max_steps=4,
        initial_features=[NaturalUnitSystem()],
    )
    assert ctx.diag.errors == [], ctx.diag.errors

    b = ctx.data.require(EnergyBudget)
    assert b.E_total == 0.0, b.E_total
    assert b.dE_dt == 0.0, b.dE_dt


def test_variable_c2_path_differs_from_const_path():
    """
    With c^2 read from WaveSpeedField, E_grad differs from the same
    field scored at constant c^2 = units.c ** 2. UpdateEMCDensity and
    UpdateWaveSpeed run first, so c^2 varies with |psi|^2.

    Mutation caught: use_variable_c2=True silently falls back to the
    constant-c^2 kernel; the two E_grad values agree.
    """
    _ti_init()
    from .allocator import AllocateWaveSpeed

    class P(Pipeline):
        def __init__(self, use_var):
            super().__init__(
                error_policy=ErrorPolicy.FAIL_FAST,
                external_provides=(UnitSystem,),
            )
            self.add(AllocateWaveField(nx=12, ny=12, nz=12, dx=1.0))
            self.add(AllocateWaveSpeed())
            self.add(_SeedStandingWave(amp=1.0, mode=1))
            self.add(UpdateEMCDensityProcessor(
                field_type=PsiLongField, beta_rho=0.1))
            self.add(UpdateWaveSpeedProcessor())
            self.add(ClearAccelerationProcessor(field_type=PsiLongField))
            self.add(LaplacianProcessor(field_type=PsiLongField))
            self.add(LeapfrogProcessor(field_type=PsiLongField))
            self.add(EnergyBudgetUpdate(
                field_type=PsiLongField, use_variable_c2=use_var))

    from ..runner import Runner
    from ..sinks import InMemorySink

    def run(use_var):
        return Runner({"session": InMemorySink()}).run(
            P(use_var), name=f"var_c2_{use_var}", params={},
            dt=0.05, max_steps=1,
            initial_features=[NaturalUnitSystem(grid_voxels_per_lambda=20,
                                                cfl_safety=0.5)],
        )

    ctx_const = run(False)
    ctx_var = run(True)
    assert ctx_const.diag.errors == [], ctx_const.diag.errors
    assert ctx_var.diag.errors == [], ctx_var.diag.errors

    e_grad_const = ctx_const.data.require(EnergyBudget).E_grad
    e_grad_var = ctx_var.data.require(EnergyBudget).E_grad

    assert e_grad_const > 0.0, e_grad_const
    assert e_grad_var > 0.0, e_grad_var
    rel = abs(e_grad_var - e_grad_const) / abs(e_grad_const)
    assert rel > 1e-3, (
        f"variable-c2 E_grad {e_grad_var} vs const {e_grad_const}, "
        f"relative difference {rel} < 1e-3"
    )


# ======================================================================
# Negative control: kinetic term is required
# ======================================================================


def test_gradient_only_energy_oscillates():
    """
    Negative control. The same standing wave, scored by E_grad alone,
    shows a step-fraction drift far above the tolerance the positive
    test uses. The plan says a gradient-only integral oscillates at
    2*omega; the observed drift is the discrete signature of that
    oscillation.

    The check compares two drift figures on the same run:
      frac_grad_only = |E_grad(end) - E_grad(start)| / |E_grad(start)|
      frac_full      = |dE_dt * dt| / |E_total|

    The gradient-only figure must be at least an order of magnitude
    larger, or the positive test's tolerance is not doing any work.

    Mutation caught: E_kin silently dropped from E_total, which makes
    frac_full match frac_grad_only and the ratio collapses to ~1.
    """
    _ti_init()
    n = 16
    dt = 0.05
    steps = 50

    # Initial E_grad on the seeded field, computed in numpy float64
    # with the same per-axis average the kernel uses.
    x = np.arange(n)
    kx = np.pi / (n - 1)
    i, j, k = np.meshgrid(x, x, x, indexing="ij")
    s = np.sin(kx * i) * np.sin(kx * j) * np.sin(kx * k)
    psi0 = s.astype(np.float64)

    def grad_energy(psi):
        gx = 0.5 * (
            (psi[2:, 1:-1, 1:-1] - psi[1:-1, 1:-1, 1:-1])
            + (psi[1:-1, 1:-1, 1:-1] - psi[:-2, 1:-1, 1:-1])
        )
        gy = 0.5 * (
            (psi[1:-1, 2:, 1:-1] - psi[1:-1, 1:-1, 1:-1])
            + (psi[1:-1, 1:-1, 1:-1] - psi[1:-1, :-2, 1:-1])
        )
        gz = 0.5 * (
            (psi[1:-1, 1:-1, 2:] - psi[1:-1, 1:-1, 1:-1])
            + (psi[1:-1, 1:-1, 1:-1] - psi[1:-1, 1:-1, :-2])
        )
        return 0.5 * (gx * gx + gy * gy + gz * gz).sum()

    ctx = _run(grid_n=n, dt=dt, max_steps=steps)
    assert ctx.diag.errors == [], ctx.diag.errors
    b = ctx.data.require(EnergyBudget)

    e_grad_initial = grad_energy(psi0)
    e_grad_final = b.E_grad

    frac_grad_only = abs(e_grad_final - e_grad_initial) / abs(e_grad_initial)
    frac_full = abs(b.dE_dt * dt) / abs(b.E_total)

    assert frac_grad_only > 10.0 * frac_full, (
        f"gradient-only drift {frac_grad_only} vs full {frac_full}; "
        f"the negative control is not distinguishable from the positive"
    )


# ======================================================================
# Build-time
# ======================================================================


def test_budget_update_requires_all_features():
    """
    EnergyBudgetUpdate.requires names its features, and the two
    flavours (constant c^2 and variable c^2) declare different tuples.

    Mutation caught: any of the required features dropped from
    requires, or WaveSpeedField added to the constant-c^2 flavour.
    """
    _ti_init()
    from .features import WaveSpeedField

    proc = EnergyBudgetUpdate(field_type=PsiLongField)
    assert WaveGrid in proc.requires, proc.requires
    assert PsiLongField in proc.requires, proc.requires
    assert EMCDensityField in proc.requires, proc.requires
    assert UnitSystem in proc.requires, proc.requires
    assert WaveSpeedField not in proc.requires, proc.requires

    proc_var = EnergyBudgetUpdate(
        field_type=PsiLongField, use_variable_c2=True)
    assert WaveSpeedField in proc_var.requires, proc_var.requires


def main() -> int:
    tests = [
        test_uniform_rho_gives_zero_deform,
        test_deficit_contributes_to_deform,
        test_standing_wave_conserves_total_energy,
        test_dE_dt_zero_on_static_state,
        test_variable_c2_path_differs_from_const_path,
        test_gradient_only_energy_oscillates,
        test_budget_update_requires_all_features,
    ]
    passed = 0
    for t in tests:
        try:
            t()
        except AssertionError as e:
            print(f"FAIL: {t.__name__}: {e}")
            traceback.print_exc()
        except Exception as e:
            print(f"ERROR: {t.__name__}: {type(e).__name__}: {e}")
            traceback.print_exc()
        else:
            print(f"PASS: {t.__name__}")
            passed += 1
    print(f"\n{passed}/{len(tests)} tests passed")
    return 0 if passed == len(tests) else 1


if __name__ == "__main__":
    sys.exit(main())