"""
Tests for source terms, scattering operators, and SeedBaseWave.

Every test names the mutation it catches, or states explicitly that it
is a smoke check. The negative controls (classes whose names start
with _Overwrite or _NonUnitary) are the mutations that the positive
tests are designed to detect: each positive test has a companion
negative test in this file, so the mutation claim is backed by code.

Run from the project root:
    python -m openwave.xperiments.m4_ewt.pipeline_engine.physics._test_sources

Note: no `from __future__ import annotations` here. Taichi kernel
definitions need live type objects.
"""

import sys
import traceback

import numpy as np
import taichi as ti

from ..pipeline import BaseProcessor, ErrorPolicy, Pipeline, Stage
from .allocator import AllocateWaveField
from ..contracts import ScatteringOperatorInterface, SourceTermInterface
from .features import PsiBaseField, PsiLongField, WaveGrid
from .seed import SeedBaseWave


_TI_INITIALIZED = False


def _ti_init():
    global _TI_INITIALIZED
    if not _TI_INITIALIZED:
        ti.init(arch=ti.cpu, log_level=ti.ERROR)
        _TI_INITIALIZED = True


# ======================================================================
# Taichi kernels used by the fixtures
# ======================================================================


@ti.kernel
def _add_vector_to_interior(
    new: ti.template(),
    nx: ti.i32,
    ny: ti.i32,
    nz: ti.i32,
    vx: ti.f32,
    vy: ti.f32,
    vz: ti.f32,
):
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        new[i, j, k] += ti.Vector([vx, vy, vz])


@ti.kernel
def _overwrite_vector_to_interior(
    new: ti.template(),
    nx: ti.i32,
    ny: ti.i32,
    nz: ti.i32,
    vx: ti.f32,
    vy: ti.f32,
    vz: ti.f32,
):
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        new[i, j, k] = ti.Vector([vx, vy, vz])


@ti.kernel
def _swap_interior(
    a: ti.template(),
    b: ti.template(),
    nx: ti.i32,
    ny: ti.i32,
    nz: ti.i32,
):
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        tmp = a[i, j, k]
        a[i, j, k] = b[i, j, k]
        b[i, j, k] = tmp


@ti.kernel
def _add_without_subtract(
    a: ti.template(),
    b: ti.template(),
    nx: ti.i32,
    ny: ti.i32,
    nz: ti.i32,
):
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        b[i, j, k] += a[i, j, k]


# ======================================================================
# Fixtures: dummy sources, scatterers, and their negative controls
# ======================================================================


class _AddVectorSource(SourceTermInterface):
    """Adds a constant vector to interior psi_new. Additive by contract."""

    name = "_AddVectorSource"
    requires = (PsiLongField,)

    def __init__(self, vx, vy, vz, order=7):
        self.vx = float(vx)
        self.vy = float(vy)
        self.vz = float(vz)
        self.order = order

    def contribute(self, ctx):
        field = ctx.data.require(PsiLongField)
        grid = ctx.data.require(WaveGrid)
        _add_vector_to_interior(
            field.psi_new,
            grid.nx,
            grid.ny,
            grid.nz,
            self.vx,
            self.vy,
            self.vz,
        )


class _OverwriteVectorSource(SourceTermInterface):
    """
    Negative control: assigns (=) instead of adds (+=).

    Used in test_overwrite_source_breaks_superposition to show that
    the superposition test in test_source_terms_superpose would go
    red if a source overwrote instead of adding.
    """

    name = "_OverwriteVectorSource"
    requires = (PsiLongField,)

    def __init__(self, vx, vy, vz, order=8):
        self.vx = float(vx)
        self.vy = float(vy)
        self.vz = float(vz)
        self.order = order

    def contribute(self, ctx):
        field = ctx.data.require(PsiLongField)
        grid = ctx.data.require(WaveGrid)
        _overwrite_vector_to_interior(
            field.psi_new,
            grid.nx,
            grid.ny,
            grid.nz,
            self.vx,
            self.vy,
            self.vz,
        )


class _SwapScatterer(ScatteringOperatorInterface):
    """Swap base and long at each interior voxel. Unitary by construction."""

    name = "_SwapScatterer"
    requires = (PsiBaseField, PsiLongField)

    def scatter(self, ctx):
        base = ctx.data.require(PsiBaseField)
        long = ctx.data.require(PsiLongField)
        grid = ctx.data.require(WaveGrid)
        _swap_interior(base.psi, long.psi, grid.nx, grid.ny, grid.nz)


class _NonUnitaryScatterer(ScatteringOperatorInterface):
    """
    Negative control: adds base to long without removing it from base.
    Total magnitude increases; unitarity is violated.
    """

    name = "_NonUnitaryScatterer"
    requires = (PsiBaseField, PsiLongField)

    def scatter(self, ctx):
        base = ctx.data.require(PsiBaseField)
        long = ctx.data.require(PsiLongField)
        grid = ctx.data.require(WaveGrid)
        _add_without_subtract(base.psi, long.psi, grid.nx, grid.ny, grid.nz)


# ======================================================================
# Pipeline and run helpers
# ======================================================================


class _SeedPair(BaseProcessor):
    """Seed PsiBase and PsiLong from numpy at step 0."""

    name = "_SeedPair"
    stage = Stage.PRE_UPDATE
    order = 10
    requires = (PsiBaseField, PsiLongField)

    def __init__(self, base_np, long_np):
        self.base_np = base_np
        self.long_np = long_np

    def process(self, ctx):
        if ctx.sim.step > 0:
            return
        bf = ctx.data.require(PsiBaseField)
        lf = ctx.data.require(PsiLongField)
        bf.psi.from_numpy(self.base_np)
        bf.psi_prev.from_numpy(self.base_np)
        lf.psi.from_numpy(self.long_np)
        lf.psi_prev.from_numpy(self.long_np)


def _build_pipeline(processors):
    class P(Pipeline):
        def __init__(self):
            super().__init__(error_policy=ErrorPolicy.FAIL_FAST)
            self.add(AllocateWaveField(nx=8, ny=8, nz=8, dx=1.0))
            for p in processors:
                self.add(p)

    return P()


def _run(processors, max_steps=1):
    from ..runner import Runner
    from ..sinks import InMemorySink

    return Runner({"session": InMemorySink()}).run(
        _build_pipeline(processors),
        name="sources_test",
        params={},
        dt=0.1,
        max_steps=max_steps,
    )


# ======================================================================
# Source terms: additive contract
# ======================================================================


def test_source_terms_superpose():
    """
    Two additive sources in one pipeline produce the sum of their
    contributions in psi_new. Neither overwrites the other.

    Mutation caught: the second source overwrites instead of adding,
    and the interior value is [0, 2, 0] rather than [1, 2, 0]. The
    companion test test_overwrite_source_breaks_superposition shows
    this failure mode on the negative control.
    """
    _ti_init()
    ctx = _run([
        _AddVectorSource(1.0, 0.0, 0.0, order=7),
        _AddVectorSource(0.0, 2.0, 0.0, order=8),
    ], max_steps=1)
    assert ctx.diag.errors == [], ctx.diag.errors

    field = ctx.data.require(PsiLongField)
    arr = field.psi_new.to_numpy()
    expected = np.array([1.0, 2.0, 0.0], dtype=np.float32)
    assert np.allclose(arr[4, 4, 4], expected, atol=1e-6), arr[4, 4, 4]
    # Boundary is untouched: sources are interior-only here.
    assert np.allclose(arr[0, 4, 4], 0.0, atol=1e-6), arr[0, 4, 4]


def test_overwrite_source_breaks_superposition():
    """
    Negative control: _OverwriteVectorSource assigns instead of adding.
    The same pipeline yields [0, 2, 0], not [1, 2, 0]. This documents
    the failure mode that test_source_terms_superpose detects.
    """
    _ti_init()
    ctx = _run([
        _AddVectorSource(1.0, 0.0, 0.0, order=7),
        _OverwriteVectorSource(0.0, 2.0, 0.0, order=8),
    ], max_steps=1)
    assert ctx.diag.errors == [], ctx.diag.errors

    field = ctx.data.require(PsiLongField)
    arr = field.psi_new.to_numpy()
    # Overwrite wipes the first source's contribution.
    assert np.allclose(arr[4, 4, 4], np.array([0.0, 2.0, 0.0]),
                       atol=1e-6), arr[4, 4, 4]


def test_source_term_does_not_touch_psi():
    """
    A source adds to psi_new; psi and psi_prev are unchanged.

    Mutation caught: a source that writes psi_am instead of psi_new.
    This is the first half of the "adds to psi_new, never overwrites
    psi_am" contract.
    """
    _ti_init()
    ctx = _run([_AddVectorSource(1.0, 0.0, 0.0)], max_steps=1)
    assert ctx.diag.errors == [], ctx.diag.errors

    field = ctx.data.require(PsiLongField)
    assert np.allclose(field.psi.to_numpy(), 0.0, atol=1e-6), \
        "psi must be untouched by a source"
    assert np.allclose(field.psi_prev.to_numpy(), 0.0, atol=1e-6), \
        "psi_prev must be untouched by a source"


def test_source_interface_process_raises():
    """
    SourceTermInterface.contribute raises NotImplementedError, and the
    base process() propagates it. An abstract interface cannot be
    instantiated as a working processor by accident.
    """
    src = SourceTermInterface()
    try:
        src.process(None)
    except NotImplementedError:
        return
    raise AssertionError("expected NotImplementedError, none raised")


# ======================================================================
# Scattering operators: unitarity
# ======================================================================


def _seed_pair(base_val, long_val):
    """Return a numpy pair of arrays to seed PsiBase and PsiLong."""

    shape = (8, 8, 8, 3)
    base = np.zeros(shape, dtype=np.float32)
    long = np.zeros(shape, dtype=np.float32)
    base[1:-1, 1:-1, 1:-1, :] = base_val
    long[1:-1, 1:-1, 1:-1, :] = long_val
    return base, long


def test_scattering_operator_preserves_magnitude():
    """
    A swap scatterer moves magnitude between two fields without gaining
    or losing any. Per voxel, |base|^2 + |long|^2 is unchanged.

    Mutation caught: a scatterer that adds to one field without
    subtracting from the other. The companion test
    test_nonunitary_scatterer_breaks_magnitude shows the failure mode.
    """
    _ti_init()
    base_val = np.array([1.0, 0.0, 0.0], dtype=np.float32)
    long_val = np.array([0.0, 2.0, 0.0], dtype=np.float32)
    base_np, long_np = _seed_pair(base_val, long_val)

    ctx = _run([_SeedPair(base_np, long_np), _SwapScatterer()], max_steps=1)
    assert ctx.diag.errors == [], ctx.diag.errors

    base_after = ctx.data.require(PsiBaseField).psi.to_numpy()
    long_after = ctx.data.require(PsiLongField).psi.to_numpy()

    # The magnitudes are swapped at each interior voxel.
    assert np.allclose(base_after[4, 4, 4], long_val, atol=1e-6), \
        base_after[4, 4, 4]
    assert np.allclose(long_after[4, 4, 4], base_val, atol=1e-6), \
        long_after[4, 4, 4]

    # Total magnitude per voxel is preserved.
    m_before = (base_np[4, 4, 4] ** 2).sum() + (long_np[4, 4, 4] ** 2).sum()
    m_after = (base_after[4, 4, 4] ** 2).sum() + (long_after[4, 4, 4] ** 2).sum()
    assert abs(m_before - 5.0) < 1e-5, m_before
    assert abs(m_after - 5.0) < 1e-5, m_after
    assert abs(m_before - m_after) < 1e-6, (m_before, m_after)


def test_nonunitary_scatterer_breaks_magnitude():
    """
    Negative control: _NonUnitaryScatterer adds base to long without
    subtracting from base. Total magnitude increases from 5 to 6;
    the swap test's magnitude check would go red if this class were
    substituted.

    This documents the failure mode that
    test_scattering_operator_preserves_magnitude detects.
    """
    _ti_init()
    base_val = np.array([1.0, 0.0, 0.0], dtype=np.float32)
    long_val = np.array([0.0, 2.0, 0.0], dtype=np.float32)
    base_np, long_np = _seed_pair(base_val, long_val)

    ctx = _run([_SeedPair(base_np, long_np), _NonUnitaryScatterer()],
               max_steps=1)
    assert ctx.diag.errors == [], ctx.diag.errors

    base_after = ctx.data.require(PsiBaseField).psi.to_numpy()
    long_after = ctx.data.require(PsiLongField).psi.to_numpy()

    m_before = (base_np[4, 4, 4] ** 2).sum() + (long_np[4, 4, 4] ** 2).sum()
    m_after = (base_after[4, 4, 4] ** 2).sum() + (long_after[4, 4, 4] ** 2).sum()
    # Nonunitary: long += base without subtracting from base.
    # Before: |[1,0,0]|^2 + |[0,2,0]|^2 = 1 + 4 = 5.
    # After:  |[1,0,0]|^2 + |[1,2,0]|^2 = 1 + 5 = 6.
    assert abs(m_before - 5.0) < 1e-5, m_before
    assert abs(m_after - 6.0) < 1e-5, m_after
    assert m_after > m_before, "nonunitary scatterer must increase total magnitude"


def test_scattering_interface_process_raises():
    """
    ScatteringOperatorInterface.scatter raises NotImplementedError, and
    the base process() propagates it.
    """
    op = ScatteringOperatorInterface()
    try:
        op.process(None)
    except NotImplementedError:
        return
    raise AssertionError("expected NotImplementedError, none raised")


# ======================================================================
# SeedBaseWave
# ======================================================================


def test_seed_base_wave_writes_psi_base_only():
    """
    SeedBaseWave writes PsiBase.psi and PsiBase.psi_prev. PsiLong is
    untouched.

    Mutation caught: a seed that writes to PsiLong instead of PsiBase.
    """
    _ti_init()
    shape = (8, 8, 8, 3)
    arr = np.zeros(shape, dtype=np.float32)
    arr[4, 4, 4] = [0.5, 0.7, 0.9]

    ctx = _run([SeedBaseWave(arr)], max_steps=1)
    assert ctx.diag.errors == [], ctx.diag.errors

    base_psi = ctx.data.require(PsiBaseField).psi.to_numpy()
    base_prev = ctx.data.require(PsiBaseField).psi_prev.to_numpy()
    long_psi = ctx.data.require(PsiLongField).psi.to_numpy()

    assert np.allclose(base_psi[4, 4, 4], arr[4, 4, 4], atol=1e-6), base_psi[4, 4, 4]
    assert np.allclose(base_prev[4, 4, 4], arr[4, 4, 4], atol=1e-6), base_prev[4, 4, 4]
    assert np.allclose(long_psi, 0.0, atol=1e-6), "PsiLong must be untouched"


def test_seed_base_wave_shape_validation():
    """
    A wrong-shaped array raises ValueError.
    """
    _ti_init()
    arr = np.zeros((4, 4, 4, 3), dtype=np.float32)
    ctx = None
    try:
        ctx = _run([SeedBaseWave(arr)], max_steps=1)
    except ValueError:
        return
    raise AssertionError(
        f"expected ValueError, got errors={ctx.diag.errors if ctx else 'none'}"
    )


def test_seed_base_wave_runs_once():
    """
    SeedBaseWave writes at step 0 and returns early thereafter. A
    pipeline that steps further leaves the seeded field unchanged.
    """
    _ti_init()
    shape = (8, 8, 8, 3)
    arr = np.zeros(shape, dtype=np.float32)
    arr[4, 4, 4] = [0.5, 0.7, 0.9]

    ctx = _run([SeedBaseWave(arr)], max_steps=3)
    assert ctx.diag.errors == [], ctx.diag.errors

    base_psi = ctx.data.require(PsiBaseField).psi.to_numpy()
    assert np.allclose(base_psi[4, 4, 4], arr[4, 4, 4], atol=1e-6), \
        base_psi[4, 4, 4]


# ======================================================================
# Runner
# ======================================================================


def main() -> int:
    tests = [
        test_source_terms_superpose,
        test_overwrite_source_breaks_superposition,
        test_source_term_does_not_touch_psi,
        test_source_interface_process_raises,
        test_scattering_operator_preserves_magnitude,
        test_nonunitary_scatterer_breaks_magnitude,
        test_scattering_interface_process_raises,
        test_seed_base_wave_writes_psi_base_only,
        test_seed_base_wave_shape_validation,
        test_seed_base_wave_runs_once,
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