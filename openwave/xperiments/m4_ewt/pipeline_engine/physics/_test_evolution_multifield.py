"""
Tests for 1.3a: parameterised LaplacianProcessor, LeapfrogProcessor,
and ClearAccelerationProcessor. Plus two structural tests for the Q3
migration of the wave speed c from WaveGrid to UnitSystem.

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
from .evolution import (
    ClearAccelerationProcessor,
    LaplacianProcessor,
    LeapfrogProcessor,
)
from .units import NaturalUnitSystem, UnitSystem

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
        _seed_both(base.psi, base.psi_prev, long.psi, long.psi_prev, grid.nx, grid.ny, grid.nz)


@ti.kernel
def _seed_both(
    base: ti.template(),
    base_prev: ti.template(),
    long: ti.template(),
    long_prev: ti.template(),
    nx: ti.i32,
    ny: ti.i32,
    nz: ti.i32,
):
    for i, j, k in ti.ndrange(nx, ny, nz):
        kx = 2.0 * ti.math.pi / ti.cast(nx, ti.f32)
        ky = 2.0 * ti.math.pi / ti.cast(ny, ti.f32)
        kz = 3.0 * ti.math.pi / ti.cast(nz, ti.f32)
        b = ti.Vector(
            [
                ti.sin(kx * i) + ti.cos(ky * j),
                ti.sin(kx * i) + ti.sin(kz * k),
                ti.cos(ky * j) + ti.sin(kz * k),
            ]
        )
        l = ti.Vector(
            [
                ti.cast(i, ti.f32),
                ti.cast(2 * j, ti.f32),
                ti.cast(3 * k, ti.f32),
            ]
        )
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
            super().__init__(external_provides=(UnitSystem,))
            self.add(AllocateWaveField(nx=16, ny=16, nz=16, dx=1.0))
            self.add(_SeedBothFields())
            self.add(ClearAccelerationProcessor())
            self.add(LaplacianProcessor())
            self.add(LeapfrogProcessor())

    return P()


def _build_parameterised_pipeline():
    class P(Pipeline):
        def __init__(self):
            super().__init__(external_provides=(UnitSystem,))
            self.add(AllocateWaveField(nx=16, ny=16, nz=16, dx=1.0))
            self.add(_SeedBothFields())
            self.add(ClearAccelerationProcessor(field_type=PsiBaseField))
            self.add(LaplacianProcessor(field_type=PsiBaseField))
            self.add(LeapfrogProcessor(field_type=PsiBaseField))

    return P()


def _run(pipeline, max_steps=1):
    from ..runner import Runner
    from ..sinks import InMemorySink

    return Runner({"session": InMemorySink()}).run(
        pipeline,
        name="multifield_test",
        params={},
        dt=0.1,
        max_steps=max_steps,
        initial_features=[NaturalUnitSystem()],
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
    assert ClearAccelerationProcessor().field_type is PsiLongField


def test_requires_reflects_field_type():
    """
    Mutation: requires hardcoded to (WaveGrid, PsiLongField) on the
    class -> fails for a parameterised instance.
    """
    for proc_class in (ClearAccelerationProcessor, LaplacianProcessor, LeapfrogProcessor):
        p = proc_class(field_type=PsiBaseField)
        assert WaveGrid in p.requires, proc_class.__name__
        assert PsiBaseField in p.requires, proc_class.__name__
        assert PsiLongField not in p.requires, proc_class.__name__


def test_wavegrid_has_no_c_attribute():
    """
    Q3 structural test. The wave speed c lives in UnitSystem, not in
    WaveGrid.

    Mutation: c added back to WaveGrid -> the two sources of c return,
    and a future LaplacianProcessor could silently read the wrong one.
    """
    grid = WaveGrid(nx=4, ny=4, nz=4, dx=0.1)
    assert not hasattr(grid, "c")


def test_laplacian_requires_unitsystem():
    """
    Q3 structural test. LaplacianProcessor declares UnitSystem in its
    requires, so Pipeline._validate forces the caller to provide one.

    Mutation: UnitSystem dropped from requires -> the processor can no
    longer read c from the canonical source, and a future
    reintroduction of grid.c would go unnoticed.
    """
    lp = LaplacianProcessor()
    assert UnitSystem in lp.requires


def test_pipeline_error_when_parameterised_field_missing():
    """
    Mutation: Pipeline._validate ignores instance requires -> no
    PipelineError, then KeyError later. Check fails on missing error.
    """

    class NeverProvided:
        pass

    class Bad(Pipeline):
        def __init__(self):
            super().__init__(error_policy=ErrorPolicy.FAIL_FAST, external_provides=(UnitSystem,))
            self.add(AllocateWaveField(nx=8, ny=8, nz=8, dx=1.0))
            self.add(LaplacianProcessor(field_type=NeverProvided))

    try:
        _run(Bad(), max_steps=1)
    except PipelineError as e:
        assert "NeverProvided" in str(e), str(e)
        return
    raise AssertionError("expected PipelineError, none raised")


# ======================================================================
# Taichi tests: accumulator contract
# ======================================================================


def test_clear_acceleration_kernel_zeros_field():
    """
    Direct kernel test. Full grid, interior and boundary, must end at zero.

    Mutation: _clear_accel writes nothing -> the field retains 999.0.
    """
    _ti_init()
    from .evolution import _clear_accel

    N = 8
    field = ti.Vector.field(3, dtype=ti.f32, shape=(N, N, N))
    field.fill(999.0)
    _clear_accel(field, N, N, N)

    arr = field.to_numpy()
    assert np.abs(arr).max() < 1e-6, np.abs(arr).max()


def test_laplacian_accumulates_on_repeated_application():
    """
    Contract test for += semantics: applying _laplacian twice to a
    non-zero accumulator doubles the result at each interior voxel.

    Mutation: _laplacian changed back to = (overwrite) -> the second
    application gives the same value as the first, check fails.
    """
    _ti_init()
    from .evolution import _laplacian

    N = 16
    field = ti.field(dtype=ti.f32, shape=(N, N, N))
    acc = ti.field(dtype=ti.f32, shape=(N, N, N))

    @ti.kernel
    def _fill(f: ti.template()):
        for i, j, k in f:
            f[i, j, k] = ti.sin(0.3 * ti.cast(i, ti.f32)) + 0.5 * ti.cos(0.2 * ti.cast(j, ti.f32))

    _fill(field)
    # acc starts at zero (Taichi default). Apply once, snapshot.
    _laplacian(field, acc, N, N, N, 1.0, 1.0)
    one = acc.to_numpy()[8, 8, 8]
    # Apply again without clearing.
    _laplacian(field, acc, N, N, N, 1.0, 1.0)
    two = acc.to_numpy()[8, 8, 8]

    assert abs(two - 2.0 * one) < 1e-4, (one, two)


# ======================================================================
# Taichi tests: parameterisation and isolation
# ======================================================================


def test_laplacian_default_writes_long():
    """
    Smoke check: the default pipeline runs one step without error, and
    PsiLong.psi_new has the vector-field shape (16, 16, 16, 3).

    This test cannot catch a no-op LaplacianProcessor.process. Kept as
    a shape and no-crash check.
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
    PsiLong.psi_new moves off zero, check fails.

    PsiLong.psi and PsiLong.psi_prev are seeded by _SeedBothFields to
    the asymmetric ramp [i, 2j, 3k]; PsiLong.psi_new is left at zero by
    the seed and stays zero under a correct field_type=PsiBaseField
    Laplacian. A wiring error that writes to PsiLong.psi_new would
    move it off zero.
    """
    _ti_init()
    ctx = _run(_build_parameterised_pipeline(), max_steps=1)
    assert ctx.diag.errors == [], ctx.diag.errors

    long_field = ctx.data.require(PsiLongField)
    psi_arr = long_field.psi.to_numpy()
    prev_arr = long_field.psi_prev.to_numpy()
    new_arr = long_field.psi_new.to_numpy()
    for i in (1, 8, 14):
        for j in (1, 8, 14):
            for k in (1, 8, 14):
                expected = np.array([i, 2 * j, 3 * k], dtype=np.float32)
                assert np.allclose(psi_arr[i, j, k], expected, atol=1e-6), (
                    i,
                    j,
                    k,
                    psi_arr[i, j, k],
                    expected,
                )
                assert np.allclose(prev_arr[i, j, k], expected, atol=1e-6), (
                    i,
                    j,
                    k,
                    prev_arr[i, j, k],
                    expected,
                )
                assert np.allclose(new_arr[i, j, k], 0.0, atol=1e-6), (
                    i,
                    j,
                    k,
                    new_arr[i, j, k],
                    "must be 0",
                )


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
                assert np.allclose(arr[i, j, k], expected, atol=1e-6), (
                    i,
                    j,
                    k,
                    arr[i, j, k],
                    expected,
                )


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
    # Taichi initialises to zero, but make it explicit for the += kernel.
    scratch.fill(0.0)

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
                expected = f_i * math.sin(kx * i) + f_j * math.cos(ky * j) + f_k * math.sin(kz * k)
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
            super().__init__(external_provides=(UnitSystem,))
            self.add(AllocateWaveField(nx=16, ny=16, nz=16, dx=1.0))
            self.add(_SeedBothFields())
            self.add(ClearAccelerationProcessor(field_type=PsiBaseField))
            self.add(ClearAccelerationProcessor(field_type=PsiLongField))
            self.add(LaplacianProcessor(field_type=PsiBaseField))
            self.add(LaplacianProcessor(field_type=PsiLongField))
            self.add(LeapfrogProcessor(field_type=PsiBaseField))
            self.add(LeapfrogProcessor(field_type=PsiLongField))

    ctx_b = _run(P(), max_steps=10)
    assert ctx_b.diag.errors == [], ctx_b.diag.errors
    base_b = ctx_b.data.require(PsiBaseField).psi.to_numpy()

    max_diff = float(np.abs(base_a - base_b).max())
    assert max_diff < 1e-6, f"cross-talk: max diff {max_diff}"


def test_boundary_voxels_unchanged_by_leapfrog_swap():
    """
    Boundary voxels of PsiLong.psi must stay at their seed value after
    several Leapfrog steps, because the swap loop is interior-only.

    Mutation: the swap loop reverted to full grid -> boundary voxels of
    psi receive values from psi_new, which is never written there by
    Laplacian (interior-only), so they would become 0 (initial new) and
    the check fails.
    """
    _ti_init()

    class P(Pipeline):
        def __init__(self):
            super().__init__(external_provides=(UnitSystem,))
            self.add(AllocateWaveField(nx=16, ny=16, nz=16, dx=1.0))
            self.add(_SeedBothFields())
            self.add(ClearAccelerationProcessor())
            self.add(LaplacianProcessor())
            self.add(LeapfrogProcessor())

    ctx = _run(P(), max_steps=5)
    assert ctx.diag.errors == [], ctx.diag.errors

    long_field = ctx.data.require(PsiLongField)
    arr = long_field.psi.to_numpy()
    # Boundary voxels: i = 0 and i = 15, etc. Seed value there is
    # [0, 0, 0] for i=0 and [15, 30, 45] for i=15 etc.
    for i, j, k in [(0, 8, 8), (15, 8, 8), (8, 0, 8), (8, 15, 8), (8, 8, 0), (8, 8, 15)]:
        expected = np.array([i, 2 * j, 3 * k], dtype=np.float32)
        assert np.allclose(arr[i, j, k], expected, atol=1e-6), (i, j, k, arr[i, j, k], expected)


def test_leapfrog_kernel_matches_numpy():
    """
    Direct kernel test of _leapfrog against a numpy reference.

    Mutation: each of the following leaves the pipeline-level tests
    green but changes the output on this arena:
      - sign error in the dt2 * accel term
      - full-grid scope where the correct scope is interior-only
      - time-level swap order inverted
    """
    _ti_init()
    from .evolution import _leapfrog

    nx, ny, nz = 6, 8, 10
    dt2 = 0.37
    rng = np.random.default_rng(seed=607_001)

    psi_np = rng.standard_normal((nx, ny, nz, 3)).astype(np.float32)
    prev_np = rng.standard_normal((nx, ny, nz, 3)).astype(np.float32)
    new_np = rng.standard_normal((nx, ny, nz, 3)).astype(np.float32)

    # Arena invariant: prev != psi on the boundary. Random data
    # guarantees it; assert so a future change to the seed does not
    # silently remove the arm that catches a full-grid swap.
    boundary = np.zeros((nx, ny, nz), dtype=bool)
    boundary[0, :, :] = True
    boundary[-1, :, :] = True
    boundary[:, 0, :] = True
    boundary[:, -1, :] = True
    boundary[:, :, 0] = True
    boundary[:, :, -1] = True
    assert not np.allclose(psi_np[boundary], prev_np[boundary], atol=1e-6), (
        "arena error: prev == psi on boundary; a full-grid swap would "
        "not change any boundary voxel and this test could not see it"
    )

    psi = ti.Vector.field(3, dtype=ti.f32, shape=(nx, ny, nz))
    prev = ti.Vector.field(3, dtype=ti.f32, shape=(nx, ny, nz))
    new = ti.Vector.field(3, dtype=ti.f32, shape=(nx, ny, nz))
    psi.from_numpy(psi_np)
    prev.from_numpy(prev_np)
    new.from_numpy(new_np)

    _leapfrog(psi, prev, new, nx, ny, nz, dt2)

    psi_out = psi.to_numpy()
    prev_out = prev.to_numpy()
    new_out = new.to_numpy()

    # Numpy reference: interior-only integration, interior-only swap.
    # new_ref carries the integrated values, not the input new.
    interior = (slice(1, -1),) * 3
    psi_ref = psi_np.copy()
    prev_ref = prev_np.copy()
    new_ref = new_np.copy()
    new_ref[interior] = 2.0 * psi_np[interior] - prev_np[interior] + dt2 * new_np[interior]
    prev_ref[interior] = psi_np[interior]
    psi_ref[interior] = new_ref[interior]

    # Interior: all three buffers must match the reference.
    assert np.allclose(
        psi_out[interior], psi_ref[interior], atol=1e-5
    ), f"psi interior max diff {np.abs(psi_out[interior] - psi_ref[interior]).max()}"
    assert np.allclose(
        prev_out[interior], prev_ref[interior], atol=1e-5
    ), f"prev interior max diff {np.abs(prev_out[interior] - prev_ref[interior]).max()}"
    assert np.allclose(
        new_out[interior], new_ref[interior], atol=1e-5
    ), f"new interior max diff {np.abs(new_out[interior] - new_ref[interior]).max()}"

    # Boundary: untouched. This is where a full-grid scope fires.
    assert np.allclose(
        psi_out[boundary], psi_np[boundary], atol=1e-6
    ), "psi boundary changed; the swap is not interior-only"
    assert np.allclose(
        prev_out[boundary], prev_np[boundary], atol=1e-6
    ), "prev boundary changed; the swap is not interior-only"
    assert np.allclose(
        new_out[boundary], new_np[boundary], atol=1e-6
    ), "new boundary changed; the integration is not interior-only"


# ======================================================================
# Runner
# ======================================================================


def main() -> int:
    tests = [
        test_default_field_type_is_long,
        test_requires_reflects_field_type,
        test_wavegrid_has_no_c_attribute,
        test_laplacian_requires_unitsystem,
        test_pipeline_error_when_parameterised_field_missing,
        test_clear_acceleration_kernel_zeros_field,
        test_laplacian_accumulates_on_repeated_application,
        test_laplacian_default_writes_long,
        test_laplacian_on_base_does_not_touch_long,
        test_leapfrog_on_base_leaves_long_bit_identical,
        test_laplacian_matches_analytic_on_asymmetric_seed,
        test_two_fields_evolve_independently,
        test_boundary_voxels_unchanged_by_leapfrog_swap,
        test_leapfrog_kernel_matches_numpy,
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
