"""
Tests for 1.2: TrackerFields allocation, shared-buffer contract with
EMCDensityField.rho, and ThreePlaneSampler.

Includes a contract fixture (_MockTrackersUpdate) that exercises the
full tracker surface with trivial formulas. Real formulas belong to
the physics layer.

Run from the project root:
    python -m openwave.xperiments.m4_ewt.pipeline_engine.physics._test_sampling

Note: no `from __future__ import annotations` here. Taichi kernel
definitions (_mock_update) need live type objects.
"""

import sys
import traceback

import taichi as ti

from .features import (
    EMCDensityField,
    PsiLongField,
    TrackerFields,
    WaveGrid,
)
from ..pipeline import BaseProcessor, Stage
from ..utils.sampling import ThreePlaneSampler

# ======================================================================
# Taichi init
# ======================================================================

_TI_INITIALIZED = False


def _ti_init():
    """Initialise Taichi once per process (Taichi 1.7.x re-init warnings)."""
    global _TI_INITIALIZED
    if not _TI_INITIALIZED:
        ti.init(arch=ti.cpu, log_level=ti.ERROR)
        _TI_INITIALIZED = True


# ======================================================================
# Contract fixture
# ======================================================================


class _MockTrackersUpdate(BaseProcessor):
    """
    Contract fixture. Exercises the full 1.2 surface with trivial
    formulas:

      - reads WaveGrid, PsiLongField.psi, EMCDensityField.rho
      - full grid loop writing every per-voxel field
      - writes the three scalar globals via ThreePlaneSampler
        (amp_global as the RMS, sqrt(<amp_local^2>))

    Real formulas (EMA, zero-crossing, rho*V*(fA)^2) belong to the
    physics layer (2.5 / 2.6) and are deliberately not implemented here.
    """

    name = "_MockTrackersUpdate"
    stage = Stage.MEASURE
    order = 100
    requires = (WaveGrid, PsiLongField, EMCDensityField, TrackerFields)
    provides = ()
    stateless = True

    def setup(self, ctx) -> None:
        grid = ctx.data.require(WaveGrid)
        self._sampler = ThreePlaneSampler(grid.nx, grid.ny, grid.nz)
        self._amp_sq = ti.field(dtype=ti.f32, shape=(grid.nx, grid.ny, grid.nz))

    def process(self, ctx) -> None:
        grid = ctx.data.require(WaveGrid)
        psi = ctx.data.require(PsiLongField).psi
        rho = ctx.data.require(EMCDensityField).rho
        trk = ctx.data.require(TrackerFields)

        _mock_update(
            psi,
            rho,
            trk.amp_local,
            trk.freq_local,
            trk.energy_long_local,
            trk.energy_trans_local,
            trk.rho_local,
            trk.last_crossing,
            grid.nx,
            grid.ny,
            grid.nz,
            float(ctx.sim.step),
        )

        _square(trk.amp_local, self._amp_sq)
        trk.amp_global[None] = self._sampler.average(self._amp_sq) ** 0.5
        trk.freq_global[None] = self._sampler.average(trk.freq_local)
        trk.energy_global[None] = self._sampler.average(trk.energy_long_local)


@ti.kernel
def _mock_update(
    psi: ti.template(),
    rho: ti.template(),
    amp: ti.template(),
    freq: ti.template(),
    energy_long: ti.template(),
    energy_trans: ti.template(),
    rho_local: ti.template(),
    last_crossing: ti.template(),
    nx: ti.i32,
    ny: ti.i32,
    nz: ti.i32,
    step: ti.f32,
):
    for i, j, k in ti.ndrange(nx, ny, nz):
        s = psi[i, j, k].norm_sqr()
        # Nonzero sentinels: every field is allocated at zero, so a zero
        # expectation cannot tell a write from no write.
        # Odd i adds 2.0, so amp_local varies over the planes and its RMS
        # differs from its mean.
        amp[i, j, k] = s + step + 2.0 * ti.cast(i % 2, ti.f32)
        freq[i, j, k] = 1.0
        energy_long[i, j, k] = s + 2.0
        energy_trans[i, j, k] = 3.0
        rho_local[i, j, k] = 1.0 + step
        last_crossing[i, j, k] = step


@ti.kernel
def _square(src: ti.template(), dst: ti.template()):
    for I in ti.grouped(src):
        dst[I] = src[I] * src[I]


# ======================================================================
# Pipeline helpers
# ======================================================================


def _make_pipeline():
    from ..pipeline import Pipeline
    from .allocator import AllocateTrackers, AllocateWaveField

    class P(Pipeline):
        def __init__(self) -> None:
            super().__init__()
            self.add(AllocateWaveField(nx=16, ny=16, nz=16, dx=1.0, c=1.0))
            self.add(AllocateTrackers())
            self.add(_MockTrackersUpdate())

    return P()


def _run_mock(max_steps: int = 10):
    from ..runner import Runner
    from ..sinks import InMemorySink

    runner = Runner({"session": InMemorySink()}, check_stateless=True)
    return runner.run(
        _make_pipeline(),
        name="mock_tracker_test",
        params={},
        dt=0.1,
        max_steps=max_steps,
    )


# ======================================================================
# Tests
# ======================================================================


def test_tracker_fields_keys_independently():
    from ..context import FeatureBag

    class _S:
        pass

    bag = FeatureBag()
    tf = TrackerFields(
        amp_local=_S(),
        freq_local=_S(),
        energy_long_local=_S(),
        energy_trans_local=_S(),
        rho_local=_S(),
        last_crossing=_S(),
        amp_global=_S(),
        freq_global=_S(),
        energy_global=_S(),
    )
    bag.set(tf)
    assert bag.require(TrackerFields) is tf


def test_rho_local_shares_buffer_with_emc():
    _ti_init()

    from ..pipeline import Pipeline
    from ..runner import Runner
    from ..sinks import InMemorySink
    from .allocator import AllocateTrackers, AllocateWaveField

    class P(Pipeline):
        def __init__(self) -> None:
            super().__init__()
            self.add(AllocateWaveField(nx=8, ny=8, nz=8, dx=1.0, c=1.0))
            self.add(AllocateTrackers())

    runner = Runner({"session": InMemorySink()})
    ctx = runner.run(P(), name="rho_share_test", params={}, max_steps=0)

    emc = ctx.data.require(EMCDensityField)
    trk = ctx.data.require(TrackerFields)
    assert trk.rho_local is emc.rho


def test_mock_runs_in_pipeline():
    _ti_init()
    ctx = _run_mock(max_steps=10)
    assert ctx.sim.step == 10, ctx.sim.step
    assert ctx.diag.errors == [], ctx.diag.errors


def test_mock_wrote_all_pervoxel_fields():
    _ti_init()
    ctx = _run_mock(max_steps=10)

    trk = ctx.data.require(TrackerFields)
    c = 8
    last_step_value = 9.0  # Pipeline.step increments AFTER process()
    assert abs(trk.amp_local[c, c, c] - last_step_value) < 1e-5, trk.amp_local[c, c, c]
    odd = trk.amp_local[c + 1, c, c]
    assert abs(odd - (last_step_value + 2.0)) < 1e-5, odd
    assert abs(trk.freq_local[c, c, c] - 1.0) < 1e-6, trk.freq_local[c, c, c]
    assert abs(trk.energy_long_local[c, c, c] - 2.0) < 1e-6, trk.energy_long_local[c, c, c]
    assert abs(trk.energy_trans_local[c, c, c] - 3.0) < 1e-6, trk.energy_trans_local[c, c, c]
    assert abs(trk.last_crossing[c, c, c] - last_step_value) < 1e-5, trk.last_crossing[c, c, c]


def test_mock_wrote_rho_via_shared_buffer():
    _ti_init()
    ctx = _run_mock(max_steps=10)

    emc = ctx.data.require(EMCDensityField)
    c = 8
    expected = 1.0 + 9.0
    assert abs(emc.rho[c, c, c] - expected) < 1e-5, emc.rho[c, c, c]


def test_mock_wrote_scalar_globals():
    _ti_init()
    ctx = _run_mock(max_steps=10)

    trk = ctx.data.require(TrackerFields)
    for scalar in (trk.amp_global, trk.freq_global, trk.energy_global):
        v = float(scalar[None])
        assert v == v, "NaN in scalar global"
    # psi is never seeded, so at the last step (step 9) freq_local and
    # energy_long_local are uniform and their globals equal their values.
    # amp_local is 9 at even i and 11 at odd i: on the 16^3 grid one third
    # of the sampled voxels read 11 (half of the xy and xz planes, none of
    # the yz plane at i = 8). The RMS and the plain mean of those values
    # must differ, or this check cannot tell the two aggregates apart.
    amp_rms = ((2 * 9.0**2 + 11.0**2) / 3) ** 0.5
    amp_mean = (2 * 9.0 + 11.0) / 3
    assert abs(amp_rms - amp_mean) > 1e-2
    assert abs(float(trk.amp_global[None]) - amp_rms) < 1e-4, float(trk.amp_global[None])
    assert abs(float(trk.freq_global[None]) - 1.0) < 1e-6, float(trk.freq_global[None])
    assert abs(float(trk.energy_global[None]) - 2.0) < 1e-6, float(trk.energy_global[None])


def test_mock_requires_declared():
    _ti_init()

    from ..pipeline import ErrorPolicy, Pipeline, PipelineError
    from ..runner import Runner
    from ..sinks import InMemorySink
    from .allocator import AllocateWaveField

    class P(Pipeline):
        def __init__(self) -> None:
            super().__init__(error_policy=ErrorPolicy.FAIL_FAST)
            self.add(AllocateWaveField(nx=8, ny=8, nz=8, dx=1.0, c=1.0))
            self.add(_MockTrackersUpdate())

    runner = Runner({"session": InMemorySink()})
    try:
        runner.run(P(), name="missing_trackers", params={}, max_steps=1)
    except PipelineError as e:
        assert "TrackerFields" in str(e), str(e)
        return
    raise AssertionError("expected PipelineError, none raised")


def test_mock_stateless_guard_passes():
    _ti_init()
    ctx = _run_mock(max_steps=10)
    assert not any("_MockTrackersUpdate" in e for e in ctx.diag.errors), ctx.diag.errors


# ----------------------------------------------------------------------
# Sampler
# ----------------------------------------------------------------------


def test_three_plane_sampler_uniform_field():
    _ti_init()

    field = ti.field(dtype=ti.f32, shape=(32, 32, 32))
    field.fill(2.5)
    sampler = ThreePlaneSampler(32, 32, 32)
    assert abs(sampler.average(field) - 2.5) < 1e-6


def test_three_plane_sampler_anisotropic_grid():
    """
    nx != ny != nz with f(i,j,k) = i + 2j + 3k: catches axis swaps in
    the copy kernels. A uniform field cannot, since every voxel a
    swapped kernel reads holds the same value.
    """
    _ti_init()

    nx, ny, nz = 16, 32, 64
    field = ti.field(dtype=ti.f32, shape=(nx, ny, nz))

    @ti.kernel
    def _fill(f: ti.template()):
        for i, j, k in f:
            f[i, j, k] = ti.cast(i + 2 * j + 3 * k, ti.f32)

    _fill(field)
    sampler = ThreePlaneSampler(nx, ny, nz)
    avg = sampler.average(field)

    # Analytic, computed outside the sampler. Each plane fixes one axis
    # at n // 2 and averages the other two over the [1:-1] interior.
    def m(n):
        return (n - 1) / 2  # mean of the interior indices 1..n-2

    planes = [
        ((nx - 2) * (ny - 2), m(nx) + 2 * m(ny) + 3 * (nz // 2)),  # xy at mid_z
        ((nx - 2) * (nz - 2), m(nx) + 2 * (ny // 2) + 3 * m(nz)),  # xz at mid_y
        ((ny - 2) * (nz - 2), (nx // 2) + 2 * m(ny) + 3 * m(nz)),  # yz at mid_x
    ]
    expected = sum(n * v for n, v in planes) / sum(n for n, _ in planes)
    assert abs(avg - expected) < 1e-3, (avg, expected)


def test_three_plane_sampler_odd_dimensions():
    """Odd cubic grid: midpoint = nx // 2 stays inside the [1:-1] interior."""
    _ti_init()

    field = ti.field(dtype=ti.f32, shape=(15, 15, 15))
    field.fill(-1.25)
    sampler = ThreePlaneSampler(15, 15, 15)
    assert abs(sampler.average(field) - (-1.25)) < 1e-6


def test_three_plane_sampler_odd_anisotropic():
    """Both odd dimensions and unequal axes."""
    _ti_init()

    nx, ny, nz = 15, 21, 33
    field = ti.field(dtype=ti.f32, shape=(nx, ny, nz))
    field.fill(0.5)
    sampler = ThreePlaneSampler(nx, ny, nz)
    assert abs(sampler.average(field) - 0.5) < 1e-6


def test_three_plane_sampler_linear_field_strict():
    """
    Strict analytic check for f(i,j,k) = i+j+k on a 32^3 grid.

    Interior [1:-1] indices run 1..30, mean 15.5. Each centre plane
    fixes one axis at nx//2 = 16, so the plane mean is
    15.5 + 15.5 + 16 = 47.0. Off-by-one shifts this by >= 0.5.
    """
    _ti_init()

    nx = ny = nz = 32
    field = ti.field(dtype=ti.f32, shape=(nx, ny, nz))

    @ti.kernel
    def _fill(f: ti.template()):
        for i, j, k in f:
            f[i, j, k] = ti.cast(i + j + k, ti.f32)

    _fill(field)
    sampler = ThreePlaneSampler(nx, ny, nz)
    avg = sampler.average(field)
    # Analytic, computed outside the sampler.
    interior_mean = (1 + 30) / 2  # = 15.5
    expected = interior_mean + interior_mean + (nx // 2)
    assert abs(avg - expected) < 1e-4, (avg, expected)


def test_three_plane_sampler_excludes_edges():
    """
    Outer shell at 1000, interior at 1. Any edge voxel in the average
    moves it off 1. The linear field above cannot show this: its edge
    and interior means are both 15.5 per axis.
    """
    _ti_init()

    n = 16
    field = ti.field(dtype=ti.f32, shape=(n, n, n))
    field.fill(1000.0)

    @ti.kernel
    def _fill_interior(f: ti.template()):
        for i, j, k in ti.ndrange((1, n - 1), (1, n - 1), (1, n - 1)):
            f[i, j, k] = 1.0

    _fill_interior(field)
    sampler = ThreePlaneSampler(n, n, n)
    avg = sampler.average(field)
    assert abs(avg - 1.0) < 1e-6, avg


def test_three_plane_sampler_rejects_tiny_grid():
    for dims in ((1, 8, 8), (8, 1, 8), (8, 8, 1)):
        try:
            ThreePlaneSampler(*dims)
        except ValueError:
            continue
        raise AssertionError(f"expected ValueError for dims={dims}")


# ======================================================================
# Runner
# ======================================================================


def main() -> int:
    tests = [
        test_tracker_fields_keys_independently,
        test_rho_local_shares_buffer_with_emc,
        test_mock_runs_in_pipeline,
        test_mock_wrote_all_pervoxel_fields,
        test_mock_wrote_rho_via_shared_buffer,
        test_mock_wrote_scalar_globals,
        test_mock_requires_declared,
        test_mock_stateless_guard_passes,
        test_three_plane_sampler_uniform_field,
        test_three_plane_sampler_anisotropic_grid,
        test_three_plane_sampler_odd_dimensions,
        test_three_plane_sampler_odd_anisotropic,
        test_three_plane_sampler_linear_field_strict,
        test_three_plane_sampler_excludes_edges,
        test_three_plane_sampler_rejects_tiny_grid,
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
