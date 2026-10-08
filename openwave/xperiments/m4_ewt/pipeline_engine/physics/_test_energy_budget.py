"""
Tests for EnergyBudget and EnergyBudgetUpdate (plan item 1.8).

The energy kernel uses the forward-difference edge sum of plan item
1.23, with the mean of the two leapfrog levels in time. Tests here
pin the shape (edge sum, shell edges included), the conservation
(spread over a full period and O(dt^2) convergence), and the dE_dt
reporting against a finite difference of the recorded E_total.

The negative control (gradient-only) shows why E_kin is required.

Run from the project root:
    python -m openwave.xperiments.m4_ewt.pipeline_engine.physics._test_energy_budget

Note: no `from __future__ import annotations` here. Taichi kernel
definitions need live type objects.
"""

import sys
import traceback
from dataclasses import dataclass

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
from .allocator import AllocateWaveField, AllocateWaveSpeed
from .emc import UpdateEMCDensityProcessor, UpdateWaveSpeedProcessor
from .energy import EnergyBudgetUpdate
from .evolution import (
    ClearAccelerationProcessor,
    LaplacianProcessor,
    LeapfrogProcessor,
)
from .units import UnitSystem


_TI_INITIALIZED = False


def _ti_init():
    global _TI_INITIALIZED
    if not _TI_INITIALIZED:
        ti.init(arch=ti.cpu, log_level=ti.ERROR)
        _TI_INITIALIZED = True


# ======================================================================
# Fixtures
# ======================================================================


@dataclass(frozen=True)
class _UnitsWithC(UnitSystem):
    """
    Local unit system with a chosen c. The rest is fixed so a test can
    set c != 1 without touching the shipped implementations.
    """
    c_value: float = 1.0

    @property
    def c(self): return self.c_value
    @property
    def wavelength(self): return 1.0
    @property
    def dx(self): return 0.05
    @property
    def dt(self): return 0.01
    @property
    def rho_0(self): return 1.0
    def to_physical_length(self, x): return x
    def to_physical_time(self, t): return t
    def to_physical_energy(self, E): return E
    def to_physical_density(self, r): return r


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
    """Set rho = 1 - 0.5 * exp(-r^2/8) at step 0. Sharp deficit at centre."""

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


class _RunResult:
    """
    Lightweight container for a run and its recorded histories.

    Written as a plain class rather than a namedtuple so the field
    names appear in autocomplete and read like the other fixtures in
    this file.
    """

    __slots__ = ("ctx", "E_total", "E_grad", "dE_dt")

    def __init__(self, ctx, E_total, E_grad, dE_dt):
        self.ctx = ctx
        self.E_total = E_total
        self.E_grad = E_grad
        self.dE_dt = dE_dt


def _run(grid_n=16, dx=1.0, dt=0.05, max_steps=100, amp=1.0, mode=1,
         rho_def=False, kappa=1.0, use_variable_c2=False, c_value=1.0):
    """
    Minimal energy pipeline: allocate, seed standing wave, seed rho,
    Clear, Laplacian, Leapfrog, EnergyBudgetUpdate, plus a recorder
    that captures E_total, E_grad and dE_dt each step.

    Returns a _RunResult with the ctx and the three histories.
    """
    E_total_hist = []
    E_grad_hist = []
    dE_dt_hist = []

    class _Record(BaseProcessor):
        name = "_Record"
        stage = Stage.MEASURE
        order = 6
        requires = (EnergyBudget,)

        def process(self, ctx):
            b = ctx.data.require(EnergyBudget)
            E_total_hist.append(b.E_total)
            E_grad_hist.append(b.E_grad)
            dE_dt_hist.append(b.dE_dt)

    class P(Pipeline):
        def __init__(self):
            super().__init__(
                error_policy=ErrorPolicy.FAIL_FAST,
                external_provides=(UnitSystem,),
            )
            self.add(AllocateWaveField(
                nx=grid_n, ny=grid_n, nz=grid_n, dx=dx))
            if use_variable_c2:
                self.add(AllocateWaveSpeed())
            self.add(_SeedStandingWave(amp=amp, mode=mode))
            if use_variable_c2:
                self.add(UpdateEMCDensityProcessor(
                    field_type=PsiLongField, beta_rho=0.1))
                self.add(UpdateWaveSpeedProcessor())
            elif rho_def:
                self.add(_SeedRhoDeficit())
            else:
                self.add(_SeedRhoUniform())
            self.add(ClearAccelerationProcessor(field_type=PsiLongField))
            self.add(LaplacianProcessor(field_type=PsiLongField))
            self.add(LeapfrogProcessor(field_type=PsiLongField))
            self.add(EnergyBudgetUpdate(
                field_type=PsiLongField,
                use_variable_c2=use_variable_c2,
                kappa=kappa))
            self.add(_Record())

    from ..runner import Runner
    from ..sinks import InMemorySink

    ctx = Runner({"session": InMemorySink()}).run(
        P(), name="energy_budget_test", params={}, dt=dt,
        max_steps=max_steps,
        initial_features=[_UnitsWithC(c_value)],
    )
    return _RunResult(ctx, E_total_hist, E_grad_hist, dE_dt_hist)


# ======================================================================
# Structural tests
# ======================================================================


def test_uniform_rho_gives_zero_deform():
    """
    rho = 1.0 everywhere -> (rho - 1.0) = 0 everywhere -> E_deform = 0.

    Mutation caught by this test: the deformation term reading |rho|
    instead of (rho - 1). The "deformation term dropped entirely" arm
    is caught by test_deficit_contributes_to_deform.
    """
    _ti_init()
    r = _run(grid_n=12, dt=0.05, max_steps=3)
    assert r.ctx.diag.errors == [], r.ctx.diag.errors

    b = r.ctx.data.require(EnergyBudget)
    assert b.E_deform == 0.0, b.E_deform
    assert b.E_kin > 0.0, b.E_kin
    assert b.E_grad > 0.0, b.E_grad


def test_deficit_contributes_to_deform():
    """
    rho < 1 in the centre -> E_deform > 0, and its value matches an
    independent numpy integral over the full grid.

    Mutation caught: E_deform reading a constant instead of the field.
    """
    _ti_init()
    n = 12
    r = _run(grid_n=n, dt=0.05, max_steps=2, rho_def=True, kappa=1.0)
    assert r.ctx.diag.errors == [], r.ctx.diag.errors

    b = r.ctx.data.require(EnergyBudget)
    rho_arr = r.ctx.data.require(EMCDensityField).rho.to_numpy()
    drho = rho_arr - 1.0
    expected = 0.5 * (drho * drho).sum()  # dx = 1, dv = 1, kappa = 1
    assert abs(b.E_deform - expected) < 1e-3, (b.E_deform, expected)


def test_energy_total_includes_all_three_components():
    """
    Structural check: E_total is the sum of its three components, and
    it exceeds E_grad + E_deform by exactly E_kin on an arena where
    E_kin > 0.

    Mutation caught: E_total defined as E_grad + E_deform (kinetic
    term dropped), or any component substituted for another.
    """
    _ti_init()
    r = _run(grid_n=12, dt=0.05, max_steps=5)
    assert r.ctx.diag.errors == [], r.ctx.diag.errors

    b = r.ctx.data.require(EnergyBudget)
    assert b.E_kin > 0.0, b.E_kin
    assert b.E_grad > 0.0, b.E_grad

    assert b.E_total > b.E_grad + b.E_deform, (
        f"E_total {b.E_total} not greater than E_grad + E_deform "
        f"{b.E_grad + b.E_deform}"
    )
    assert abs(b.E_total - (b.E_kin + b.E_grad + b.E_deform)) < 1e-9, (
        f"E_total {b.E_total} vs sum {b.E_kin + b.E_grad + b.E_deform}"
    )


def test_dE_dt_zero_on_static_state():
    """
    A field with psi == psi_prev and no gradient has E_total = 0 on
    both steps, so dE_dt = 0 after the first step. Smoke check that
    the reporting path runs without error.
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
        initial_features=[_UnitsWithC(1.0)],
    )
    assert ctx.diag.errors == [], ctx.diag.errors

    b = ctx.data.require(EnergyBudget)
    assert b.E_total == 0.0, b.E_total
    assert b.dE_dt == 0.0, b.dE_dt


def test_variable_c2_path_differs_from_const_path():
    """
    With c^2 read from WaveSpeedField, E_grad differs from the same
    field scored at constant c^2. The arena has a non-uniform rho via
    UpdateEMCDensity from the standing wave.

    Mutation caught: use_variable_c2=True silently falls back to the
    constant-c^2 kernel; the two E_grad values agree.
    """
    _ti_init()

    r_const = _run(grid_n=12, dt=0.05, max_steps=1,
                   use_variable_c2=False)
    r_var = _run(grid_n=12, dt=0.05, max_steps=1,
                 use_variable_c2=True)
    assert r_const.ctx.diag.errors == [], r_const.ctx.diag.errors
    assert r_var.ctx.diag.errors == [], r_var.ctx.diag.errors

    e_grad_const = r_const.ctx.data.require(EnergyBudget).E_grad
    e_grad_var = r_var.ctx.data.require(EnergyBudget).E_grad

    assert e_grad_const > 0.0, e_grad_const
    assert e_grad_var > 0.0, e_grad_var
    rel = abs(e_grad_var - e_grad_const) / abs(e_grad_const)
    assert rel > 1e-3, (
        f"variable-c2 E_grad {e_grad_var} vs const {e_grad_const}, "
        f"relative difference {rel} < 1e-3"
    )


# ======================================================================
# Conservation: per-step recording, spread, order
# ======================================================================


def test_standing_wave_conserves_total_energy():
    """
    A 3D standing wave at rest, integrated with the production chain,
    conserves E_total over one period to the O(dt^2) bar of plan item
    1.23. The measure is the spread of the per-step E_total over the
    run, (max - min) / mean.

    Mutation caught: any of the kernel bugs a single-step read could
    not see (kinetic dropped, edge sum reverted to centered, shell
    edges dropped, dE_dt hard-wired), because the spread is measured
    over the whole run.

    Reference (reviewer, N=16, dt=0.05, mean form, edges, shell
    included): spread 1.64e-4. Threshold 5e-4 leaves a 3x margin.
    """
    _ti_init()
    n = 16
    dt = 0.05
    steps = 350  # about one period for mode-1 in a 16^3 box

    r = _run(grid_n=n, dt=dt, max_steps=steps, mode=1)
    assert r.ctx.diag.errors == [], r.ctx.diag.errors
    assert len(r.E_total) == steps, len(r.E_total)

    arr = np.array(r.E_total)
    spread = (arr.max() - arr.min()) / abs(arr.mean())
    assert spread < 5e-4, (
        f"E_total spread over {steps} steps = {spread}, "
        f"expected < 5e-4 (mean form, edge sum)"
    )


def test_energy_conservation_is_second_order_in_dt():
    """
    The same standing wave at dt and dt/2. Under the mean form the
    spread is O(dt^2), so the ratio must be about 4.

    Mutation caught: the kernel using the level-n+1-only gradient
    (first order in dt, ratio about 2), or a mixer that is not the
    mean. A one-sided gradient gives spread about 1.8e-2 at dt=0.05,
    a hundred times the mean's 1.6e-4.

    Threshold: ratio >= 3, so a first-order path fails.
    """
    _ti_init()
    n = 16
    dt = 0.05
    steps_coarse = 350
    steps_fine = 700

    r1 = _run(grid_n=n, dt=dt, max_steps=steps_coarse, mode=1)
    r2 = _run(grid_n=n, dt=dt / 2, max_steps=steps_fine, mode=1)
    assert r1.ctx.diag.errors == [], r1.ctx.diag.errors
    assert r2.ctx.diag.errors == [], r2.ctx.diag.errors

    a1 = np.array(r1.E_total)
    a2 = np.array(r2.E_total)
    spread1 = (a1.max() - a1.min()) / abs(a1.mean())
    spread2 = (a2.max() - a2.min()) / abs(a2.mean())

    ratio = spread1 / spread2 if spread2 > 0.0 else float("inf")
    assert ratio >= 3.0, (
        f"spread(dt)={spread1:.3e}, spread(dt/2)={spread2:.3e}, "
        f"ratio {ratio:.2f} < 3 (expected ~4 for O(dt^2))"
    )


def test_dE_dt_matches_finite_difference():
    """
    dE_dt reported by the processor equals the finite difference of
    the recorded E_total, sign included, for every step after the
    first.

    Mutation caught: dE_dt sign flipped, dE_dt divided by dt twice,
    E_prev updated before the difference, or the whole dE_dt path
    short-circuited.
    """
    _ti_init()
    n = 12
    dt = 0.05
    steps = 20

    r = _run(grid_n=n, dt=dt, max_steps=steps, mode=1)
    assert r.ctx.diag.errors == [], r.ctx.diag.errors
    assert len(r.E_total) == steps
    assert len(r.dE_dt) == steps

    for step in range(1, steps):
        expected = (r.E_total[step] - r.E_total[step - 1]) / dt
        got = r.dE_dt[step]
        assert abs(got - expected) < 1e-9, (
            f"step {step}: dE_dt {got} vs finite difference {expected}"
        )


def test_conservation_at_nonunit_dx():
    """
    Same standing wave, run at dx = 0.5 (four times smaller cells, so
    four times more edges per axis). The mean edge sum must still
    conserve to the same O(dt^2) bar.

    Mutation caught: any kernel that drops the dx factor, hard-codes
    dx = 1, or reads the wrong exponent for dv.

    Threshold: 1e-3 here, looser than the 5e-4 the dx=1 test uses.
    The discrete dispersion contributes a (k dx)^2 term whose relative
    share grows as dx shrinks, so the spread at dx=0.5 is measured
    near 6.6e-4 against the same arena's 1.6e-4 at dx=1. A kernel
    that dropped dx, hard-coded dx=1, or read the wrong dv exponent
    gives a spread near 1e-1, two orders of magnitude above this bar.
    """
    _ti_init()
    n = 16
    dx = 0.5
    dt = 0.05
    steps = 350

    r = _run(grid_n=n, dx=dx, dt=dt, max_steps=steps, mode=1)
    assert r.ctx.diag.errors == [], r.ctx.diag.errors

    arr = np.array(r.E_total)
    spread = (arr.max() - arr.min()) / abs(arr.mean())
    assert spread < 1e-3, (
        f"E_total spread at dx=0.5 over {steps} steps = {spread}, "
        f"expected < 1e-3"
    )


# ======================================================================
# Negative control: kinetic term is required
# ======================================================================


def test_gradient_only_energy_oscillates():
    """
    Negative control. The same standing wave, scored by E_grad alone,
    has a per-step spread far above the tolerance the positive test
    uses. Plan item 1.23 says a gradient-only integral oscillates at
    2*omega; the measured spread is the discrete signature of that
    oscillation.

    Both figures are spread over the same run, so they are directly
    comparable. The gradient-only spread must be at least ten times
    the full E_total spread.

    Mutation caught: E_kin silently dropped from E_total. That makes
    the two spreads agree and the ratio collapse to 1.
    """
    _ti_init()
    n = 16
    dt = 0.05
    steps = 350

    r = _run(grid_n=n, dt=dt, max_steps=steps, mode=1)
    assert r.ctx.diag.errors == [], r.ctx.diag.errors

    g = np.array(r.E_grad)
    e = np.array(r.E_total)
    spread_grad_only = (g.max() - g.min()) / abs(g.mean())
    spread_full = (e.max() - e.min()) / abs(e.mean())

    assert spread_grad_only > 10.0 * spread_full, (
        f"gradient-only spread {spread_grad_only:.3e} vs full "
        f"{spread_full:.3e}; the negative control is not "
        f"distinguishable from the positive"
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
        test_energy_total_includes_all_three_components,
        test_dE_dt_zero_on_static_state,
        test_variable_c2_path_differs_from_const_path,
        test_standing_wave_conserves_total_energy,
        test_energy_conservation_is_second_order_in_dt,
        test_dE_dt_matches_finite_difference,
        test_conservation_at_nonunit_dx,
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