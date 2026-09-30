"""
Tests for 1.3a: parameterised LaplacianProcessor and LeapfrogProcessor.

Every test names the mutation it catches, or states explicitly that it
is a smoke check and cannot catch a specific mutation.

Run from the project root:
    python -m openwave.xperiments.m4_ewt.pipeline_engine.physics._test_evolution_multifield

Note: no `from __future__ import annotations` here. Taichi kernel
definitions need live type objects.
"""

import math
import sys
import traceback

import numpy as np
import taichi as ti

from .features import (
    PsiBaseField,
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
from .evolution import LaplacianProcessor, LeapfrogProcessor


# ======================================================================
# Taichi init
# ======================================================================

_TI_INITIALIZED = False


def _ti_init():
    global _TI_INITIALIZED
    if not _TI_INITIALIZED:
        ti.init(arch=ti.cpu, log_level=ti.ERROR)
        _TI_INITIALIZED = True


# ======================================================================
# Test fixture: deterministic seeds
# ======================================================================


class _SeedBothFields(BaseProcessor):
    """
    Seeds PsiBase and PsiLong with DIFFERENT patterns on the first step.
    Base is non-linear (sin/cos), long is linear (asymmetric ramp).
    Different patterns are deliberate: if a processor touches the wrong
    field, tests can see it.
    """
    name = "_SeedBothFields"
    stage = Stage.PRE_UPDATE
    order = 10
    requires = (WaveGrid, PsiBaseField, PsiLongField)
    provides = ()

    def process(self, ctx) -> None:
        if ctx.sim.step > 0:
            return
        grid = ctx.data.require(WaveGrid)
        base = ctx.data.require(PsiBaseField)
        long = ctx.data.require(PsiLongField)
        _seed_both(base.psi, base.psi_prev, long.psi, long.psi_prev,
                   grid.nx, grid.ny, grid.nz)


@ti.kernel
def _seed_both(base: ti.template(), base_prev: ti.template(),
               long: ti.template(), long_prev: ti.template(),
               nx: ti.i32, ny: ti.i32, nz: ti.i32):
    for i, j, k in ti.ndrange(nx, ny, nz):
        kx = 2.0 * ti.math.pi / ti.cast(nx, ti.f32)
        ky = 2.0 * ti.math.pi / ti.cast(ny, ti.f32)
        kz = 3.0 * ti.math.pi / ti.cast(nz, ti.f32)
        b = ti.Vector([
            ti.sin(kx * i) + ti.cos(ky * j),
            ti.sin(kx * i) + ti.sin(kz * k),
            ti.cos(ky * j) + ti.sin(kz * k),
        ])
        l = ti.Vector([
            ti.cast(i, ti.f32),
            ti.cast(2 * j, ti.f32),
            ti.cast(3 * k, ti.f32),
        ])
        base[i, j, k] = b
        base_prev[i, j, k] = b
        long[i, j, k] = l
        long_prev[i, j, k] = l


# ======================================================================
# Pipeline builders
# ======================================================================


def _build_default_pipeline():
    class P(Pipeline):
        def __init__(self):
            super().__init__()
            self.add(AllocateWaveField(nx=16, ny=16, nz=16, dx=1.0, c=1.0))
            self.add(_SeedBothFields())
            self.add(LaplacianProcessor())
            self.add(LeapfrogProcessor())
    return P()


def _build_parameterised_pipeline():
    class P(Pipeline):
        def __init__(self):
            super().__init__()
            self.add(AllocateWaveField(nx=16, ny=16, nz=16, dx=1.0, c=1.0))
            self.add(_SeedBothFields())
            self.add(LaplacianProcessor(field_type=PsiBaseField))
            self.add(LeapfrogProcessor(field_type=PsiBaseField))
    return P()


def _run(pipeline, max_steps=1):
    from ..runner import Runner
    from ..sinks import InMemorySink
    return Runner({"session": InMemorySink()}).run(
        pipeline, name="multifield_test", params={}, dt=0.1, max_steps=max_steps,
    )


# ======================================================================
# Pure-Python tests (no Taichi needed)
# ======================================================================


def test_default_field_type_is_long():
    """
    Mutation: default changed to PsiBaseField or PsiTransField -> fails.
    """
    assert LaplacianProcessor().field_type is PsiLongField
    assert LeapfrogProcessor().field_type is PsiLongField


def test_requires_reflects_field_type():
    """
    Mutation: requires hardcoded to (WaveGrid, PsiLongField) on the
    class -> fails for a parameterised instance.
    """
    lp = LaplacianProcessor(field_type=PsiBaseField)
    lf = LeapfrogProcessor(field_type=PsiBaseField)
    assert WaveGrid in lp.requires
    assert PsiBaseField in lp.requires
    assert PsiLongField not in lp.requires
    assert WaveGrid in lf.requires
    assert PsiBaseField in lf.requires
    assert PsiLongField not in lf.requires


def test_pipeline_error_when_parameterised_field_missing():
    """
    Mutation: Pipeline._validate ignores instance requires -> no
    PipelineError, then KeyError later. Check fails on missing error.
    """
    class NeverProvided:
        pass

    class Bad(Pipeline):
        def __init__(self):
            super().__init__(error_policy=ErrorPolicy.FAIL_FAST)
            self.add(AllocateWaveField(nx=8, ny=8, nz=8, dx=1.0, c=1.0))
            self.add(LaplacianProcessor(field_type=NeverProvided))

    try:
        _run(Bad(), max_steps=1)
    except PipelineError as e:
        assert "NeverProvided" in str(e), str(e)
        return
    raise AssertionError("expected PipelineError, none raised")


# ======================================================================
# Taichi tests
# ======================================================================


def test_laplacian_default_writes_long():
    """
    Smoke check: default pipeline runs one step without error, and
    PsiLong.psi_new has the vector-field shape (16, 16, 16, 3).

    This test cannot catch a no-op LaplacianProcessor.process. The
    default field_type is PsiLongField and the PsiLong seed is linear,
    so a correct Laplacian and a no-op both leave psi_new at zero.
    Kept as a shape and no-crash check.
    """
    _ti_init()
    ctx = _run(_build_default_pipeline(), max_steps=1)
    assert ctx.diag.errors == [], ctx.diag.errors

    long_field = ctx.data.require(PsiLongField)
    arr = long_field.psi_new.to_numpy()
    assert arr.shape == (16, 16, 16, 3), arr.shape


def test_laplacian_on_base_does_not_touch_long():
    """
    Mutation: Laplacian writes to PsiLong regardless of field_type ->
    PsiLong.psi changes from its seed, check fails.
    """
    _ti_init()
    ctx = _run(_build_parameterised_pipeline(), max_steps=1)
    assert ctx.diag.errors == [], ctx.diag.errors

    long_field = ctx.data.require(PsiLongField)
    arr = long_field.psi.to_numpy()
    for i in (1, 8, 14):
        for j in (1, 8, 14):
            for k in (1, 8, 14):
                expected = np.array([i, 2 * j, 3 * k], dtype=np.float32)
                assert np.allclose(arr[i, j, k], expected, atol=1e-6), \
                    (i, j, k, arr[i, j, k], expected)


def test_leapfrog_on_base_leaves_long_bit_identical():
    """
    Mutation: Leapfrog integrates both fields -> PsiLong.psi drifts
    from the seed after several steps, check fails.
    """
    _ti_init()
    ctx = _run(_build_parameterised_pipeline(), max_steps=5)
    assert ctx.diag.errors == [], ctx.diag.errors

    long_field = ctx.data.require(PsiLongField)
    arr = long_field.psi.to_numpy()
    for i in (1, 8, 14):
        for j in (1, 8, 14):
            for k in (1, 8, 14):
                expected = np.array([i, 2 * j, 3 * k], dtype=np.float32)
                assert np.allclose(arr[i, j, k], expected, atol=1e-6), \
                    (i, j, k, arr[i, j, k], expected)


def test_laplacian_matches_analytic_on_asymmetric_seed():
    """
    Seed psi = sin(kx i) + cos(ky j) + sin(kz k) with three distinct
    wavenumbers. Laplacian is computed analytically using the discrete
    second-derivative eigenvalue -4 sin^2(k/2).

    Mutation: any axis swap (i<->j, i<->k, j<->k) in the kernel reads
    a different profile than the analytic expectation, check fails.
    """
    _ti_init()

    N = 32
    dx = 1.0
    c = 1.0

    field = ti.field(dtype=ti.f32, shape=(N, N, N))
    scratch = ti.field(dtype=ti.f32, shape=(N, N, N))

    @ti.kernel
    def _fill(f: ti.template()):
        kx = 2.0 * ti.math.pi / ti.cast(N, ti.f32)
        ky = 2.0 * ti.math.pi / ti.cast(N, ti.f32)
        kz = 3.0 * ti.math.pi / ti.cast(N, ti.f32)
        for i, j, k in f:
            f[i, j, k] = ti.sin(kx * i) + ti.cos(ky * j) + ti.sin(kz * k)

    _fill(field)

    from .evolution import _laplacian
    _laplacian(field, scratch, N, N, N, dx, c)

    out = scratch.to_numpy()
    kx = 2.0 * math.pi / N
    ky = 2.0 * math.pi / N
    kz = 3.0 * math.pi / N
    f_i = -4.0 * math.sin(kx / 2.0) ** 2
    f_j = -4.0 * math.sin(ky / 2.0) ** 2
    f_k = -4.0 * math.sin(kz / 2.0) ** 2
    worst = 0.0
    for i in (5, 16, 25):
        for j in (5, 16, 25):
            for k in (5, 16, 25):
                expected = (f_i * math.sin(kx * i)
                            + f_j * math.cos(ky * j)
                            + f_k * math.sin(kz * k))
                actual = out[i, j, k]
                worst = max(worst, abs(actual - expected))
    assert worst < 1e-3, f"worst analytic error {worst}"


def test_two_fields_evolve_independently():
    """
    Run A evolves PsiBase only. Run B evolves both PsiBase and PsiLong.
    PsiBase.psi after Run A must equal PsiBase.psi after Run B.

    Mutation: any cross-talk between fields makes Run B's PsiBase
    differ from Run A's, check fails.
    """
    _ti_init()

    ctx_a = _run(_build_parameterised_pipeline(), max_steps=10)
    assert ctx_a.diag.errors == [], ctx_a.diag.errors
    base_a = ctx_a.data.require(PsiBaseField).psi.to_numpy()

    class P(Pipeline):
        def __init__(self):
            super().__init__()
            self.add(AllocateWaveField(nx=16, ny=16, nz=16, dx=1.0, c=1.0))
            self.add(_SeedBothFields())
            self.add(LaplacianProcessor(field_type=PsiBaseField))
            self.add(LaplacianProcessor(field_type=PsiLongField))
            self.add(LeapfrogProcessor(field_type=PsiBaseField))
            self.add(LeapfrogProcessor(field_type=PsiLongField))

    ctx_b = _run(P(), max_steps=10)
    assert ctx_b.diag.errors == [], ctx_b.diag.errors
    base_b = ctx_b.data.require(PsiBaseField).psi.to_numpy()

    max_diff = float(np.abs(base_a - base_b).max())
    assert max_diff < 1e-6, f"cross-talk: max diff {max_diff}"


# ======================================================================
# Runner
# ======================================================================


def main() -> int:
    tests = [
        test_default_field_type_is_long,
        test_requires_reflects_field_type,
        test_pipeline_error_when_parameterised_field_missing,
        test_laplacian_default_writes_long,
        test_laplacian_on_base_does_not_touch_long,
        test_leapfrog_on_base_leaves_long_bit_identical,
        test_laplacian_matches_analytic_on_asymmetric_seed,
        test_two_fields_evolve_independently,
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