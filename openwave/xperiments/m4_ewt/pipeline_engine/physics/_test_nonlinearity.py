"""
Tests for NonlinearCubic: additive contract, interior-only scope,
composition with Laplacian.

Every test names the mutation it catches.

Run from the project root:
    python -m openwave.xperiments.m4_ewt.pipeline_engine.physics._test_nonlinearity

Note: no `from __future__ import annotations` here. Taichi kernel
definitions need live type objects.
"""

import sys
import traceback

import numpy as np
import taichi as ti

from .features import PsiLongField, WaveGrid
from ..pipeline import BaseProcessor, Pipeline, Stage
from .allocator import AllocateWaveField
from .evolution import (
    ClearAccelerationProcessor,
    LaplacianProcessor,
    LeapfrogProcessor,
)
from .nonlinearity import NonlinearCubic
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


def test_apply_cubic_adds_to_existing_accel():
    """
    The kernel must use += so it composes with Laplacian.

    Mutation: += changed to = -> accel at interior becomes only the
    nonlinear term, not (pre-filled sentinel + nonlinear term), check
    fails.
    """
    _ti_init()
    from .nonlinearity import _apply_cubic

    N = 16
    psi = ti.Vector.field(3, dtype=ti.f32, shape=(N, N, N))
    accel = ti.Vector.field(3, dtype=ti.f32, shape=(N, N, N))

    @ti.kernel
    def _fill_psi(f: ti.template()):
        for i, j, k in f:
            f[i, j, k] = ti.Vector([0.1, 0.2, 0.3])

    _fill_psi(psi)
    accel.fill(7.0)

    gamma = 2.0
    _apply_cubic(psi, accel, gamma, N, N, N)

    c = 8
    out = accel.to_numpy()[c, c, c]
    # Nonlinear contribution: -gamma * |psi|^2 * psi = -2 * 0.14 * [.1,.2,.3]
    # = [-0.028, -0.056, -0.084]
    expected = np.array([7.0 - 0.028, 7.0 - 0.056, 7.0 - 0.084], dtype=np.float32)
    assert np.allclose(out, expected, atol=1e-5), (out, expected)


def test_apply_cubic_accumulates_on_repeated_application():
    """
    Applying twice must double the result.

    Mutation: += changed to = -> second application gives the same
    value as the first, check fails.
    """
    _ti_init()
    from .nonlinearity import _apply_cubic

    N = 16
    psi = ti.Vector.field(3, dtype=ti.f32, shape=(N, N, N))
    accel = ti.Vector.field(3, dtype=ti.f32, shape=(N, N, N))

    @ti.kernel
    def _fill_psi(f: ti.template()):
        for i, j, k in f:
            f[i, j, k] = ti.Vector([0.1, 0.2, 0.3])

    _fill_psi(psi)
    accel.fill(0.0)

    gamma = 2.0
    _apply_cubic(psi, accel, gamma, N, N, N)
    one = accel.to_numpy()[8, 8, 8].copy()
    _apply_cubic(psi, accel, gamma, N, N, N)
    two = accel.to_numpy()[8, 8, 8].copy()

    assert np.allclose(two, 2.0 * one, atol=1e-5), (one, two)


def test_apply_cubic_touches_interior_only():
    """
    Boundary voxels of accel must be untouched.

    Mutation: ti.ndrange((1, nx-1), ...) changed to (nx, ny, nz) ->
    boundary voxels change from the sentinel, check fails.
    """
    _ti_init()
    from .nonlinearity import _apply_cubic

    N = 16
    psi = ti.Vector.field(3, dtype=ti.f32, shape=(N, N, N))
    accel = ti.Vector.field(3, dtype=ti.f32, shape=(N, N, N))

    @ti.kernel
    def _fill_psi(f: ti.template()):
        for i, j, k in f:
            f[i, j, k] = ti.Vector([0.1, 0.2, 0.3])

    _fill_psi(psi)
    accel.fill(999.0)

    _apply_cubic(psi, accel, 2.0, N, N, N)
    arr = accel.to_numpy()

    for i, j, k in [(0, 8, 8), (15, 8, 8), (8, 0, 8), (8, 15, 8), (8, 8, 0), (8, 8, 15)]:
        assert np.allclose(arr[i, j, k], 999.0, atol=1e-6), (i, j, k, arr[i, j, k])


# ======================================================================
# Pipeline-level tests
# ======================================================================


def test_nonlinear_zero_gamma_in_pipeline_is_noop():
    """
    With gamma=0 the nonlinearity contributes nothing; the pipeline
    runs and the field is unchanged compared to a pipeline without it.

    Mutation: kernel writes something regardless of gamma -> the two
    runs diverge, check fails.
    """
    _ti_init()

    def _build(with_nonlinear):
        class P(Pipeline):
            def __init__(self):
                super().__init__(external_provides=(UnitSystem,))
                self.add(AllocateWaveField(nx=8, ny=8, nz=8, dx=1.0))
                self.add(_SeedLinear())
                self.add(ClearAccelerationProcessor())
                self.add(LaplacianProcessor())
                if with_nonlinear:
                    self.add(NonlinearCubic(gamma=0.0))
                self.add(LeapfrogProcessor())

        return P()

    ctx_a = _run(_build(False), max_steps=5)
    ctx_b = _run(_build(True), max_steps=5)
    a = ctx_a.data.require(PsiLongField).psi.to_numpy()
    b = ctx_b.data.require(PsiLongField).psi.to_numpy()
    assert np.abs(a - b).max() < 1e-6, np.abs(a - b).max()


def test_nonlinear_composes_with_laplacian_in_pipeline():
    """
    After Clear + Laplacian + NonlinearCubic, accel at an interior voxel
    equals (Laplacian contribution) + (nonlinear contribution), computed
    independently from a harmonic seed with non-zero Laplacian.

    Mutation caught: NonlinearCubic uses = instead of +=. With the
    harmonic seed the Laplacian contribution is non-zero, so a dropped
    term changes the result. A ramp seed would have zero Laplacian and
    this test could not see the difference.
    """
    _ti_init()

    class P(Pipeline):
        def __init__(self):
            super().__init__(external_provides=(UnitSystem,))
            self.add(AllocateWaveField(nx=16, ny=16, nz=16, dx=1.0))
            self.add(_SeedLinear())
            self.add(ClearAccelerationProcessor())
            self.add(LaplacianProcessor())
            self.add(NonlinearCubic(gamma=0.5))

    ctx = _run(P(), max_steps=1)
    assert ctx.diag.errors == [], ctx.diag.errors

    long_field = ctx.data.require(PsiLongField)
    psi = long_field.psi.to_numpy()
    accel = long_field.psi_new.to_numpy()

    # Independent computation at an interior voxel.
    c = 8
    dx2 = 1.0
    lap = (
        psi[c + 1, c, c]
        + psi[c - 1, c, c]
        + psi[c, c + 1, c]
        + psi[c, c - 1, c]
        + psi[c, c, c + 1]
        + psi[c, c, c - 1]
        - 6.0 * psi[c, c, c]
    ) / dx2
    u = float(np.dot(psi[c, c, c], psi[c, c, c]))
    nl = -0.5 * u * psi[c, c, c]
    expected = lap + nl
    assert np.allclose(accel[c, c, c], expected, atol=1e-4), (accel[c, c, c], expected)


# ======================================================================
# Helpers
# ======================================================================


class _SeedLinear(BaseProcessor):
    """Linear ramp seed, deterministic and asymmetric per axis."""

    name = "_SeedLinear"
    stage = Stage.PRE_UPDATE
    order = 10
    requires = (WaveGrid, PsiLongField)
    provides = ()

    def process(self, ctx) -> None:
        if ctx.sim.step > 0:
            return
        grid = ctx.data.require(WaveGrid)
        field = ctx.data.require(PsiLongField)
        _seed_linear(field.psi, field.psi_prev, grid.nx, grid.ny, grid.nz)


@ti.kernel
def _seed_linear(psi: ti.template(), prev: ti.template(), nx: ti.i32, ny: ti.i32, nz: ti.i32):
    """
    Harmonic seed, asymmetric per axis.

    A ramp seed (i, 2j, 3k) has zero Laplacian, so the composition
    test cannot tell += from = after a clear. This seed has a
    non-zero discrete Laplacian on the interior, which is what makes
    the composition visible.
    """
    kx = 2.0 * ti.math.pi / 8.0
    ky = 2.0 * ti.math.pi / 6.0
    kz = 2.0 * ti.math.pi / 5.0
    for i, j, k in ti.ndrange(nx, ny, nz):
        v = ti.Vector(
            [
                0.3 * ti.sin(kx * ti.cast(i, ti.f32)),
                0.2 * ti.cos(ky * ti.cast(j, ti.f32)),
                0.1 * ti.sin(kz * ti.cast(k, ti.f32)),
            ]
        )
        psi[i, j, k] = v
        prev[i, j, k] = v


def _run(pipeline, max_steps=1):
    from ..runner import Runner
    from ..sinks import InMemorySink

    return Runner({"session": InMemorySink()}).run(
        pipeline,
        name="nonlinear_test",
        params={},
        dt=0.1,
        max_steps=max_steps,
        initial_features=[NaturalUnitSystem()],
    )


# ======================================================================
# Runner
# ======================================================================


def main() -> int:
    tests = [
        test_apply_cubic_adds_to_existing_accel,
        test_apply_cubic_accumulates_on_repeated_application,
        test_apply_cubic_touches_interior_only,
        test_nonlinear_zero_gamma_in_pipeline_is_noop,
        test_nonlinear_composes_with_laplacian_in_pipeline,
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
