"""
Tests for BoundaryCondition and BoundaryProcessor (plan item 1.7).

Three kinds: dirichlet, periodic, reflecting. Each has a test on a
prepared arena that makes the kind's behaviour visible.

No energy budget is checked. Flux accounting through the boundary is
plan item 1.8 (Phase C) and does not exist yet.

Run from the project root:
    python -m openwave.xperiments.m4_ewt.pipeline_engine.physics._test_boundary

Note: no `from __future__ import annotations` here. Taichi kernel
definitions need live type objects.
"""

import sys
import traceback

import numpy as np
import taichi as ti

from .features import BoundaryCondition, PsiLongField, WaveGrid
from ..pipeline import (
    BaseProcessor,
    ErrorPolicy,
    Pipeline,
    PipelineError,
    Stage,
)
from .allocator import AllocateWaveField
from .boundary import BoundaryProcessor


_TI_INITIALIZED = False


def _ti_init():
    global _TI_INITIALIZED
    if not _TI_INITIALIZED:
        ti.init(arch=ti.cpu, log_level=ti.ERROR)
        _TI_INITIALIZED = True


class _SeedField(BaseProcessor):
    """Seed PsiLong.psi and psi_prev from a numpy array at step 0."""

    name = "_SeedField"
    stage = Stage.PRE_UPDATE
    order = 10
    requires = (PsiLongField,)

    def __init__(self, arr):
        self.arr = arr

    def process(self, ctx):
        if ctx.sim.step > 0:
            return
        f = ctx.data.require(PsiLongField)
        f.psi.from_numpy(self.arr)
        f.psi_prev.from_numpy(self.arr)


def _build_pipeline(seed_arr, bc, max_steps=1):
    class P(Pipeline):
        def __init__(self):
            super().__init__(
                error_policy=ErrorPolicy.FAIL_FAST,
                external_provides=(BoundaryCondition,),
            )
            self.add(AllocateWaveField(nx=6, ny=6, nz=6, dx=1.0))
            self.add(_SeedField(seed_arr))
            self.add(BoundaryProcessor(field_type=PsiLongField))

    from ..runner import Runner
    from ..sinks import InMemorySink

    return Runner({"session": InMemorySink()}).run(
        P(), name="boundary_test", params={}, dt=0.1, max_steps=max_steps,
        initial_features=[bc],
    )


# ======================================================================
# BoundaryCondition feature
# ======================================================================


def test_boundary_condition_accepts_three_kinds():
    for kind in ("dirichlet", "periodic", "reflecting"):
        bc = BoundaryCondition(kind=kind)
        assert bc.kind == kind


def test_boundary_condition_rejects_unknown_kind():
    for kind in ("banana", "absorbing"):
        try:
            BoundaryCondition(kind=kind)
        except ValueError:
            continue
        raise AssertionError(f"expected ValueError for kind={kind!r}")


def test_boundary_condition_rejects_nonpositive_r_domain():
    for bad in (0.0, -1.0):
        try:
            BoundaryCondition(kind="dirichlet", r_domain=bad)
        except ValueError:
            continue
        raise AssertionError(f"expected ValueError for r_domain={bad}")


# ======================================================================
# BoundaryProcessor: build-time and behaviour
# ======================================================================


def test_boundary_processor_requires_boundary_condition():
    """
    A pipeline that registers BoundaryProcessor but not the feature
    fails at build time.

    Mutation caught: requires missing the BoundaryCondition feature.
    """
    _ti_init()

    class Bad(Pipeline):
        def __init__(self):
            super().__init__(error_policy=ErrorPolicy.FAIL_FAST)
            self.add(AllocateWaveField(nx=6, ny=6, nz=6, dx=1.0))
            self.add(BoundaryProcessor(field_type=PsiLongField))

    from ..runner import Runner
    from ..sinks import InMemorySink

    try:
        Runner({"session": InMemorySink()}).run(
            Bad(), name="missing_bc", params={}, dt=0.1, max_steps=1,
        )
    except PipelineError as e:
        assert "BoundaryCondition" in str(e), str(e)
        return
    raise AssertionError("expected PipelineError, none raised")


def test_dirichlet_zeroes_outer_shell():
    """
    Dirichlet kind: outer shell of psi and psi_prev is zero after one
    step, interior unchanged.

    Mutation caught: dirichlet zeros the interior, or does not zero
    the shell.
    """
    _ti_init()
    arr = np.ones((6, 6, 6, 3), dtype=np.float32)
    ctx = _build_pipeline(arr, BoundaryCondition("dirichlet"), max_steps=1)
    assert ctx.diag.errors == [], ctx.diag.errors

    psi = ctx.data.require(PsiLongField).psi.to_numpy()
    # Boundary shell is zero.
    assert np.allclose(psi[0, :, :], 0.0, atol=1e-6)
    assert np.allclose(psi[-1, :, :], 0.0, atol=1e-6)
    assert np.allclose(psi[:, 0, :], 0.0, atol=1e-6)
    assert np.allclose(psi[:, -1, :], 0.0, atol=1e-6)
    assert np.allclose(psi[:, :, 0], 0.0, atol=1e-6)
    assert np.allclose(psi[:, :, -1], 0.0, atol=1e-6)
    # Interior is one (the seed value).
    assert np.allclose(psi[3, 3, 3], 1.0, atol=1e-6)


def test_periodic_wraps_opposite_face():
    """
    Periodic kind: shell value is copied from the opposite interior
    face. Prepare a field with distinct values at x=1 and x=nx-2; the
    shell at x=0 should equal x=nx-2, and shell at x=nx-1 should equal
    x=1.

    Mutation caught: periodic copies from the wrong face.
    """
    _ti_init()
    arr = np.zeros((6, 6, 6, 3), dtype=np.float32)
    arr[1, 3, 3] = [5.0, 0.0, 0.0]
    arr[4, 3, 3] = [7.0, 0.0, 0.0]
    ctx = _build_pipeline(arr, BoundaryCondition("periodic"), max_steps=1)
    assert ctx.diag.errors == [], ctx.diag.errors

    psi = ctx.data.require(PsiLongField).psi.to_numpy()
    # x=0 gets x=4 value; x=5 gets x=1 value.
    assert np.allclose(psi[0, 3, 3], [7.0, 0.0, 0.0], atol=1e-6), psi[0, 3, 3]
    assert np.allclose(psi[5, 3, 3], [5.0, 0.0, 0.0], atol=1e-6), psi[5, 3, 3]


def test_reflecting_mirrors_interior():
    """
    Reflecting kind: shell value is copied from the immediate interior
    neighbour. Prepare arr[1,3,3] = 5; shell at x=0 should be 5.

    Mutation caught: reflecting copies from the wrong distance.
    """
    _ti_init()
    arr = np.zeros((6, 6, 6, 3), dtype=np.float32)
    arr[1, 3, 3] = [5.0, 0.0, 0.0]
    arr[4, 3, 3] = [7.0, 0.0, 0.0]
    ctx = _build_pipeline(arr, BoundaryCondition("reflecting"), max_steps=1)
    assert ctx.diag.errors == [], ctx.diag.errors

    psi = ctx.data.require(PsiLongField).psi.to_numpy()
    assert np.allclose(psi[0, 3, 3], [5.0, 0.0, 0.0], atol=1e-6), psi[0, 3, 3]
    assert np.allclose(psi[5, 3, 3], [7.0, 0.0, 0.0], atol=1e-6), psi[5, 3, 3]


def test_kind_override_beats_feature():
    """
    kind_override on the processor wins over the feature's kind. This
    is the test-only path documented in BoundaryProcessor.

    Mutation caught: kind_override ignored.
    """
    _ti_init()
    arr = np.ones((6, 6, 6, 3), dtype=np.float32)

    # Feature says periodic; processor override says dirichlet.
    # With dirichlet override, the shell must be zeroed.
    class P(Pipeline):
        def __init__(self):
            super().__init__(
                error_policy=ErrorPolicy.FAIL_FAST,
                external_provides=(BoundaryCondition,),
            )
            self.add(AllocateWaveField(nx=6, ny=6, nz=6, dx=1.0))
            self.add(_SeedField(arr))
            self.add(BoundaryProcessor(
                field_type=PsiLongField, kind_override="dirichlet"))

    from ..runner import Runner
    from ..sinks import InMemorySink
    ctx = Runner({"session": InMemorySink()}).run(
        P(), name="override_test", params={}, dt=0.1, max_steps=1,
        initial_features=[BoundaryCondition("periodic")],
    )
    assert ctx.diag.errors == [], ctx.diag.errors
    psi = ctx.data.require(PsiLongField).psi.to_numpy()
    # Dirichlet won: shell is zero.
    assert np.allclose(psi[0, :, :], 0.0, atol=1e-6)


def _numpy_boundary(kind, p):
    """Reference for the three kernels: every face, any grid shape."""
    p = p.copy()
    if kind == "dirichlet":
        p[0] = p[-1] = 0.0
        p[:, 0] = p[:, -1] = 0.0
        p[:, :, 0] = p[:, :, -1] = 0.0
    elif kind == "periodic":
        p[0], p[-1] = p[-2], p[1]
        p[:, 0], p[:, -1] = p[:, -2], p[:, 1]
        p[:, :, 0], p[:, :, -1] = p[:, :, -2], p[:, :, 1]
    else:
        p[0], p[-1] = p[1], p[-2]
        p[:, 0], p[:, -1] = p[:, 1], p[:, -2]
        p[:, :, 0], p[:, :, -1] = p[:, :, 1], p[:, :, -2]
    return p


def test_kinds_match_numpy_reference_all_faces():
    """
    Each kind against a numpy reference on a non-cubic grid (5, 7, 9),
    with random psi and a different random psi_prev, applied to
    PsiBaseField while PsiLongField carries a random sentinel near 1e6.
    The sentinel is not uniform, so a copy kernel run on PsiLong moves it.

    The tests above seed values on the x axis of a cubic grid, so they
    see neither the y and z faces nor an nx/ny/nz mix-up.

    Mutation caught: a y or z face dropped or copied from the wrong
    cell; ny and nz swapped; psi_prev left alone; the processor
    ignoring field_type (the sentinel moves).
    """
    _ti_init()
    from .features import PsiBaseField

    nx, ny, nz = 5, 7, 9
    rng = np.random.default_rng(seed=615)
    psi0 = rng.standard_normal((nx, ny, nz, 3)).astype(np.float32)
    prev0 = rng.standard_normal((nx, ny, nz, 3)).astype(np.float32)
    sentinel = (1e6 + rng.standard_normal((nx, ny, nz, 3))).astype(np.float32)

    class _SeedBaseAndSentinel(BaseProcessor):
        name = "_SeedBaseAndSentinel"
        stage = Stage.PRE_UPDATE
        order = 10
        requires = (PsiBaseField, PsiLongField)

        def process(self, ctx):
            if ctx.sim.step > 0:
                return
            base = ctx.data.require(PsiBaseField)
            long = ctx.data.require(PsiLongField)
            base.psi.from_numpy(psi0)
            base.psi_prev.from_numpy(prev0)
            long.psi.from_numpy(sentinel)
            long.psi_prev.from_numpy(sentinel)

    class P(Pipeline):
        def __init__(self):
            super().__init__(
                error_policy=ErrorPolicy.FAIL_FAST,
                external_provides=(BoundaryCondition,),
            )
            self.add(AllocateWaveField(nx=nx, ny=ny, nz=nz, dx=1.0))
            self.add(_SeedBaseAndSentinel())
            self.add(BoundaryProcessor(field_type=PsiBaseField))

    from ..runner import Runner
    from ..sinks import InMemorySink

    for kind in ("dirichlet", "periodic", "reflecting"):
        ctx = Runner({"session": InMemorySink()}).run(
            P(), name=f"reference_{kind}", params={}, dt=0.1, max_steps=1,
            initial_features=[BoundaryCondition(kind)],
        )
        assert ctx.diag.errors == [], ctx.diag.errors
        base = ctx.data.require(PsiBaseField)
        long = ctx.data.require(PsiLongField)
        d_psi = np.abs(base.psi.to_numpy() - _numpy_boundary(kind, psi0)).max()
        d_prev = np.abs(base.psi_prev.to_numpy() - _numpy_boundary(kind, prev0)).max()
        d_long = np.abs(long.psi.to_numpy() - sentinel).max()
        assert d_psi == 0.0, f"{kind}: psi differs from reference by {d_psi}"
        assert d_prev == 0.0, f"{kind}: psi_prev differs from reference by {d_prev}"
        assert d_long == 0.0, f"{kind}: PsiLong sentinel moved by {d_long}"


# ======================================================================
# Runner
# ======================================================================


def main() -> int:
    tests = [
        test_boundary_condition_accepts_three_kinds,
        test_boundary_condition_rejects_unknown_kind,
        test_boundary_condition_rejects_nonpositive_r_domain,
        test_boundary_processor_requires_boundary_condition,
        test_dirichlet_zeroes_outer_shell,
        test_periodic_wraps_opposite_face,
        test_reflecting_mirrors_interior,
        test_kind_override_beats_feature,
        test_kinds_match_numpy_reference_all_faces,
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