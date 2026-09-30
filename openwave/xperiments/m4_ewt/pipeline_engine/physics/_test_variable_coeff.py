"""
Tests for LaplacianVariableCoeffProcessor and the _laplacian_var_coeff
kernel: flux form correctness, accumulation contract, and the
discriminating energy-conservation test (flux form vs naive form).

Run from the project root:
    python -m openwave.xperiments.m4_ewt.pipeline_engine.physics._test_variable_coeff

Note: no `from __future__ import annotations` here. Taichi kernel
definitions need live type objects.
"""

import sys
import traceback

import numpy as np
import taichi as ti

from .features import PsiLongField, WaveGrid, WaveSpeedField
from ..pipeline import BaseProcessor, Pipeline, Stage
from .allocator import AllocateWaveField, AllocateWaveSpeed
from .emc import UpdateEMCDensityProcessor, UpdateWaveSpeedProcessor
from .evolution import (
    ClearAccelerationProcessor,
    LaplacianVariableCoeffProcessor,
    LeapfrogProcessor,
    _clear_accel,
    _laplacian_var_coeff,
)
from .units import NaturalUnitSystem, UnitSystem


_TI_INITIALIZED = False


def _ti_init():
    global _TI_INITIALIZED
    if not _TI_INITIALIZED:
        ti.init(arch=ti.cpu, log_level=ti.ERROR)
        _TI_INITIALIZED = True


# ======================================================================
# Kernel-level tests
# ======================================================================


def test_var_coeff_const_matches_const_laplacian():
    """
    c^2 = c0^2 everywhere -> _laplacian_var_coeff must match the
    constant-coefficient _laplacian at every interior voxel.

    Mutation: half-grid average wrong (e.g. c2[i+1] used directly)
    -> outputs differ, check fails.
    """
    _ti_init()
    from .evolution import _laplacian

    N = 16
    dx = 1.0
    c0 = 2.0
    c2_val = c0 * c0

    psi = ti.Vector.field(3, dtype=ti.f32, shape=(N, N, N))
    c2 = ti.field(dtype=ti.f32, shape=(N, N, N))
    acc_var = ti.Vector.field(3, dtype=ti.f32, shape=(N, N, N))
    acc_const = ti.Vector.field(3, dtype=ti.f32, shape=(N, N, N))

    @ti.kernel
    def _fill(p: ti.template(), c: ti.template()):
        kx = 2.0 * ti.math.pi / 8.0
        for i, j, k in p:
            p[i, j, k] = ti.Vector([
                ti.sin(kx * ti.cast(i, ti.f32)),
                ti.cos(kx * ti.cast(j, ti.f32)),
                0.5 * ti.sin(kx * ti.cast(k, ti.f32)),
            ])
            c[i, j, k] = c2_val

    _fill(psi, c2)
    acc_var.fill(0.0)
    acc_const.fill(0.0)

    _laplacian_var_coeff(psi, c2, acc_var, N, N, N, dx)
    _laplacian(psi, acc_const, N, N, N, dx, c0)

    a = acc_var.to_numpy()
    b = acc_const.to_numpy()
    diff = np.abs(a - b).max()
    assert diff < 1e-5, f"max diff {diff}"


def test_var_coeff_accumulates_not_overwrites():
    """
    Applying twice must double the result at each interior voxel.

    Mutation: += changed to = -> second application matches first.
    """
    _ti_init()

    N = 16
    dx = 1.0

    psi = ti.Vector.field(3, dtype=ti.f32, shape=(N, N, N))
    c2 = ti.field(dtype=ti.f32, shape=(N, N, N))
    acc = ti.Vector.field(3, dtype=ti.f32, shape=(N, N, N))

    @ti.kernel
    def _fill(p: ti.template(), c: ti.template()):
        kx = 2.0 * ti.math.pi / 8.0
        for i, j, k in p:
            p[i, j, k] = ti.Vector([
                ti.sin(kx * ti.cast(i, ti.f32)), 0.0, 0.0,
            ])
            c[i, j, k] = 1.0 + 0.5 * ti.sin(2.0 * kx * ti.cast(i, ti.f32))

    _fill(psi, c2)
    acc.fill(0.0)

    _laplacian_var_coeff(psi, c2, acc, N, N, N, dx)
    one = acc.to_numpy()[8, 8, 8].copy()
    _laplacian_var_coeff(psi, c2, acc, N, N, N, dx)
    two = acc.to_numpy()[8, 8, 8].copy()

    assert np.allclose(two, 2.0 * one, atol=1e-5), (one, two)


def test_var_coeff_missing_half_average_changes_result():
    """
    A kernel that drops the half-grid average (uses c2[i+1] directly
    instead of 0.5*(c2[i]+c2[i+1])) must give a different result on a
    variable c2 field.

    Mutation: production kernel drops the average -> production and the
    reference (correctly averaged) kernel diverge, check fails.
    """
    _ti_init()

    N = 16
    dx = 1.0

    psi = ti.Vector.field(3, dtype=ti.f32, shape=(N, N, N))
    c2 = ti.field(dtype=ti.f32, shape=(N, N, N))
    acc_prod = ti.Vector.field(3, dtype=ti.f32, shape=(N, N, N))
    acc_ref = ti.Vector.field(3, dtype=ti.f32, shape=(N, N, N))

    @ti.kernel
    def _fill(p: ti.template(), c: ti.template()):
        kx = 2.0 * ti.math.pi / 8.0
        for i, j, k in p:
            p[i, j, k] = ti.Vector([
                ti.sin(kx * ti.cast(i, ti.f32)),
                ti.cos(kx * ti.cast(j, ti.f32)),
                0.0,
            ])
            # c2 varies strongly with x
            c[i, j, k] = 1.0 + 0.5 * ti.sin(kx * ti.cast(i, ti.f32))

    _fill(psi, c2)
    acc_prod.fill(0.0)
    acc_ref.fill(0.0)

    _laplacian_var_coeff(psi, c2, acc_prod, N, N, N, dx)
    _laplacian_var_coeff_no_average(psi, c2, acc_ref, N, N, N, dx)

    a = acc_prod.to_numpy()
    b = acc_ref.to_numpy()
    diff = np.abs(a - b).max()
    assert diff > 1e-3, f"missing-average not detected, max diff {diff}"


@ti.kernel
def _laplacian_var_coeff_no_average(
    psi: ti.template(),
    c2: ti.template(),
    out: ti.template(),
    nx: ti.i32,
    ny: ti.i32,
    nz: ti.i32,
    dx: ti.f32,
):
    # Deliberately omits the half-grid average: uses c2 at the neighbor
    # directly. Different from production when c2 varies.
    inv_dx2 = 1.0 / (dx * dx)
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        c2_xp = c2[i + 1, j, k]
        c2_xm = c2[i - 1, j, k]
        c2_yp = c2[i, j + 1, k]
        c2_ym = c2[i, j - 1, k]
        c2_zp = c2[i, j, k + 1]
        c2_zm = c2[i, j, k - 1]

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


# ======================================================================
# Energy conservation
# ======================================================================


@ti.kernel
def _seed_zero_boundary(
    psi: ti.template(),
    psi_prev: ti.template(),
    c2: ti.template(),
    n: ti.i32,
):
    kx = ti.math.pi / ti.cast(n - 1, ti.f32)
    ky = ti.math.pi / ti.cast(n - 1, ti.f32)
    kz = ti.math.pi / ti.cast(n - 1, ti.f32)
    for i, j, k in psi:
        val = (
            0.1
            * ti.sin(kx * ti.cast(i, ti.f32))
            * ti.sin(ky * ti.cast(j, ti.f32))
            * ti.sin(kz * ti.cast(k, ti.f32))
        )
        v = ti.Vector([val, 0.0, 0.0])
        psi[i, j, k] = v
        psi_prev[i, j, k] = v
        c2[i, j, k] = 1.0 + 0.3 * ti.cos(2.0 * kx * ti.cast(i, ti.f32))

@ti.kernel
def _leapfrog_var_coeff_step(
    psi: ti.template(),
    psi_prev: ti.template(),
    psi_new: ti.template(),
    c2: ti.template(),
    nx: ti.i32, ny: ti.i32, nz: ti.i32,
    dx: ti.f32, dt2: ti.f32,
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
        lap = inv_dx2 * (
            c2_xp * (psi[i + 1, j, k] - psi[i, j, k])
            - c2_xm * (psi[i, j, k] - psi[i - 1, j, k])
            + c2_yp * (psi[i, j + 1, k] - psi[i, j, k])
            - c2_ym * (psi[i, j, k] - psi[i, j - 1, k])
            + c2_zp * (psi[i, j, k + 1] - psi[i, j, k])
            - c2_zm * (psi[i, j, k] - psi[i, j, k - 1])
        )
        psi_new[i, j, k] = (
            2.0 * psi[i, j, k] - psi_prev[i, j, k] + dt2 * lap
        )
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        psi_prev[i, j, k] = psi[i, j, k]
        psi[i, j, k] = psi_new[i, j, k]


@ti.kernel
def _leapfrog_naive_step(
    psi: ti.template(),
    psi_prev: ti.template(),
    psi_new: ti.template(),
    c2: ti.template(),
    nx: ti.i32, ny: ti.i32, nz: ti.i32,
    dx: ti.f32, dt2: ti.f32,
):
    inv_dx2 = 1.0 / (dx * dx)
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        lap = (
            psi[i + 1, j, k] + psi[i - 1, j, k]
            + psi[i, j + 1, k] + psi[i, j - 1, k]
            + psi[i, j, k + 1] + psi[i, j, k - 1]
            - 6.0 * psi[i, j, k]
        ) * inv_dx2
        psi_new[i, j, k] = (
            2.0 * psi[i, j, k] - psi_prev[i, j, k]
            + dt2 * c2[i, j, k] * lap
        )
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        psi_prev[i, j, k] = psi[i, j, k]
        psi[i, j, k] = psi_new[i, j, k]


def _staggered_energy(psi, psi_prev, c2, scratch, N, dx, dt):
    """
    E_stag = 0.5 * ||(psi^{n+1} - psi^n)/dt||^2
           + 0.5 * <psi^n, M psi^{n+1}>
    M = -L_var (flux form).
    """
    _clear_accel(scratch, N, N, N)
    _laplacian_var_coeff(psi, c2, scratch, N, N, N, dx)
    M_psi_np1 = -scratch.to_numpy()
    p_n = psi_prev.to_numpy()
    p_np1 = psi.to_numpy()
    v = (p_np1 - p_n) / dt
    return 0.5 * (v * v).sum() + 0.5 * (p_n * M_psi_np1).sum()


def test_flux_form_conserves_staggered_invariant():
    """
    Leapfrog + flux-form LaplacianVarCoeff, variable c^2, zero boundary
    seed. The staggered invariant is conserved over 50 steps.

    Mutation: half-grid average wrong, or += semantics broken.
    """
    _ti_init()

    N = 32
    dx = 1.0
    dt = 0.01
    n_steps = 50

    psi = ti.Vector.field(3, dtype=ti.f32, shape=(N, N, N))
    psi_prev = ti.Vector.field(3, dtype=ti.f32, shape=(N, N, N))
    psi_new = ti.Vector.field(3, dtype=ti.f32, shape=(N, N, N))
    c2 = ti.field(dtype=ti.f32, shape=(N, N, N))

    _seed_zero_boundary(psi, psi_prev, c2, N)

    _leapfrog_var_coeff_step(psi, psi_prev, psi_new, c2, N, N, N, dx, dt * dt)
    E0 = _staggered_energy(psi, psi_prev, c2, psi_new, N, dx, dt)

    Es = [E0]
    for _ in range(n_steps - 1):
        _leapfrog_var_coeff_step(psi, psi_prev, psi_new, c2, N, N, N, dx, dt * dt)
        Es.append(_staggered_energy(psi, psi_prev, c2, psi_new, N, dx, dt))

    drift = max(abs(e - E0) for e in Es) / abs(E0)
    print(f"  flux form relative drift: {drift:.3e}")
    assert drift < 1e-4, f"flux form drift too large: {drift}"


def test_naive_form_breaks_staggered_invariant():
    """
    Leapfrog + naive form c^2_i * laplacian(psi), same seed. Drift
    must be larger than the flux-form drift by a clear margin.

    Mutation: this is the discriminating test, unchanged.

    Measured on N=32, dt=0.05, 50 steps, zero-boundary seed:
        flux form drift  5.6e-06   (f32 rounding)
        naive form drift 8.2e-04   (operator asymmetry)
    Ratio ~146x. The naive threshold 5e-4 leaves a 1.6x margin over
    the measured value; if the flux form were accidentally substituted,
    naive would drop below 1e-5 and the check would fail.
    """
    _ti_init()

    N = 32
    dx = 1.0
    dt = 0.01
    n_steps = 50

    psi = ti.Vector.field(3, dtype=ti.f32, shape=(N, N, N))
    psi_prev = ti.Vector.field(3, dtype=ti.f32, shape=(N, N, N))
    psi_new = ti.Vector.field(3, dtype=ti.f32, shape=(N, N, N))
    c2 = ti.field(dtype=ti.f32, shape=(N, N, N))

    _seed_zero_boundary(psi, psi_prev, c2, N)

    _leapfrog_naive_step(psi, psi_prev, psi_new, c2, N, N, N, dx, dt * dt)
    E0 = _staggered_energy(psi, psi_prev, c2, psi_new, N, dx, dt)

    Es = [E0]
    for _ in range(n_steps - 1):
        _leapfrog_naive_step(psi, psi_prev, psi_new, c2, N, N, N, dx, dt * dt)
        Es.append(_staggered_energy(psi, psi_prev, c2, psi_new, N, dx, dt))

    drift = max(abs(e - E0) for e in Es) / abs(E0)
    print(f"  naive form relative drift: {drift:.3e}")
    assert drift > 5e-4, f"naive form drift unexpectedly small: {drift}"


# ======================================================================
# Pipeline integration
# ======================================================================


class _SeedHarmonicLong(BaseProcessor):
    name = "_SeedHarmonicLong"
    stage = Stage.PRE_UPDATE
    order = 10
    requires = (WaveGrid, PsiLongField)
    provides = ()

    def __init__(self, amp: float = 0.3):
        self.amp = amp

    def process(self, ctx) -> None:
        if ctx.sim.step > 0:
            return
        grid = ctx.data.require(WaveGrid)
        field = ctx.data.require(PsiLongField)
        _seed_harmonic_long(field.psi, field.psi_prev,
                            grid.nx, grid.ny, grid.nz, self.amp)


@ti.kernel
def _seed_harmonic_long(psi: ti.template(), prev: ti.template(),
                        nx: ti.i32, ny: ti.i32, nz: ti.i32, amp: ti.f32):
    kx = 2.0 * ti.math.pi / 8.0
    for i, j, k in ti.ndrange(nx, ny, nz):
        v = ti.Vector([amp * ti.sin(kx * ti.cast(i, ti.f32)), 0.0, 0.0])
        psi[i, j, k] = v
        prev[i, j, k] = v


def test_pipeline_runs_with_variable_coeff_chain():
    """
    Full chain: ClearAcceleration -> UpdateEMCDensity -> UpdateWaveSpeed
    -> LaplacianVariableCoeff -> Leapfrog. Runs 20 steps without error.

    Mutation: any processor's requires wrong -> PipelineError at build.
    """
    _ti_init()

    class P(Pipeline):
        def __init__(self):
            super().__init__(external_provides=(UnitSystem,))
            self.add(AllocateWaveField(nx=16, ny=16, nz=16, dx=1.0, c=1.0))
            self.add(AllocateWaveSpeed())
            self.add(_SeedHarmonicLong(amp=0.3))
            self.add(ClearAccelerationProcessor())
            self.add(UpdateEMCDensityProcessor(beta_rho=0.1))
            self.add(UpdateWaveSpeedProcessor())
            self.add(LaplacianVariableCoeffProcessor())
            self.add(LeapfrogProcessor())

    from ..runner import Runner
    from ..sinks import InMemorySink
    runner = Runner({"session": InMemorySink()})
    ctx = runner.run(
        P(), name="var_coeff_chain", params={}, dt=0.05, max_steps=20,
        initial_features=[NaturalUnitSystem()],
    )
    assert ctx.diag.errors == [], ctx.diag.errors
    assert ctx.sim.step == 20


# ======================================================================
# Runner
# ======================================================================


def main() -> int:
    tests = [
        test_var_coeff_const_matches_const_laplacian,
        test_var_coeff_accumulates_not_overwrites,
        test_var_coeff_missing_half_average_changes_result,
        test_flux_form_conserves_staggered_invariant,
        test_naive_form_breaks_staggered_invariant,
        test_pipeline_runs_with_variable_coeff_chain,
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