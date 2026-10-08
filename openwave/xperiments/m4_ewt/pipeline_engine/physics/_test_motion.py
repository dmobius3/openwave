"""
Tests for WCMotionProcessor and the two shipped drift rules (plan 1.6).

Two rules: NoOpRule (control) and GradRhoRule (F = -grad rho). Each
has a positive test on a prepared arena and a negative-control sibling
that documents the failure mode the positive test detects.

Run from the project root:
    python -m openwave.xperiments.m4_ewt.pipeline_engine.physics._test_motion

Note: no `from __future__ import annotations` here. Taichi kernel
definitions need live type objects.
"""

import sys
import traceback

import numpy as np
import taichi as ti

from .features import EMCDensityField, WaveGrid
from ..pipeline import (
    BaseProcessor,
    ErrorPolicy,
    Pipeline,
    PipelineError,
    Stage,
)
from .allocator import AllocateWaveField
from .motion import (
    GradRhoRule,
    NoOpRule,
    WCDriftRule,
    WCMotionProcessor,
)
from .wc_types import WC, WCState

_TI_INITIALIZED = False


def _ti_init():
    global _TI_INITIALIZED
    if not _TI_INITIALIZED:
        ti.init(arch=ti.cpu, log_level=ti.ERROR)
        _TI_INITIALIZED = True


# ======================================================================
# Test-local fixtures
# ======================================================================


class _SeedRho(BaseProcessor):
    """Seed EMCDensityField.rho from numpy at step 0. Test fixture."""

    name = "_SeedRho"
    stage = Stage.PRE_UPDATE
    order = 5
    requires = (EMCDensityField,)

    def __init__(self, arr):
        self.arr = arr

    def process(self, ctx):
        if ctx.sim.step > 0:
            return
        ctx.data.require(EMCDensityField).rho.from_numpy(self.arr)


class _ConstantRule(WCDriftRule):
    """Return a fixed acceleration. Test fixture."""

    def __init__(self, ax, ay, az):
        self.a = (float(ax), float(ay), float(az))

    def __call__(self, ctx, wc_index):
        return self.a


class _IndexRule(WCDriftRule):
    """
    Return an acceleration that depends on the WC's index. Test fixture
    for "the processor passes the right index".
    """

    def __call__(self, ctx, wc_index):
        return (float(wc_index), 0.0, 0.0)


class _ReverseGradRhoRule(WCDriftRule):
    """
    Negative control: F = +grad rho. Moves WCs away from the deficit.
    Used in test_reverse_grad_rho_moves_away to show that the sign in
    GradRhoRule matters and the positive test is not trivially green.
    """

    def __call__(self, ctx, wc_index):
        grid = ctx.data.require(WaveGrid)
        wc_state = ctx.data.require(WCState)
        rho = ctx.data.require(EMCDensityField).rho.to_numpy()
        wc = wc_state.centers[wc_index]
        i = max(1, min(grid.nx - 2, int(round(wc.x))))
        j = max(1, min(grid.ny - 2, int(round(wc.y))))
        k = max(1, min(grid.nz - 2, int(round(wc.z))))
        gx = 0.5 * (rho[i + 1, j, k] - rho[i - 1, j, k])
        gy = 0.5 * (rho[i, j + 1, k] - rho[i, j - 1, k])
        gz = 0.5 * (rho[i, j, k + 1] - rho[i, j, k - 1])
        return (float(gx), float(gy), float(gz))


def _build_pipeline(wcs, rule, rho_arr=None, dt=0.1, max_steps=1):
    """
    Build a minimal pipeline: allocate, seed rho if given, run the
    motion processor. WCState is passed via initial_features.
    """

    class P(Pipeline):
        def __init__(self):
            super().__init__(
                error_policy=ErrorPolicy.FAIL_FAST,
                external_provides=(WCState,),
            )
            self.add(AllocateWaveField(nx=8, ny=8, nz=8, dx=1.0))
            if rho_arr is not None:
                self.add(_SeedRho(rho_arr))
            self.add(WCMotionProcessor(rule=rule))

    from ..runner import Runner
    from ..sinks import InMemorySink

    return Runner({"session": InMemorySink()}).run(
        P(),
        name="motion_test",
        params={},
        dt=dt,
        max_steps=max_steps,
        initial_features=[WCState(centers=wcs)],
    )


# ======================================================================
# No-op rule
# ======================================================================


def test_noop_rule_leaves_position_unchanged():
    """
    Control: NoOpRule passed explicitly keeps the position at its
    starting value for any number of steps.

    The processor's default rule is NoOpRule by construction, but this
    test passes NoOpRule explicitly, so it does not exercise the
    default path. Mutation caught: a NoOpRule that returns a nonzero
    acceleration.
    """
    _ti_init()
    wc = WC(x=4.0, y=4.0, z=4.0, velocity=(0.0, 0.0, 0.0))
    ctx = _build_pipeline([wc], NoOpRule(), max_steps=5)
    assert ctx.diag.errors == [], ctx.diag.errors

    state = ctx.data.require(WCState)
    c = state.centers[0]
    assert (c.x, c.y, c.z) == (4.0, 4.0, 4.0), (c.x, c.y, c.z)
    assert c.velocity == (0.0, 0.0, 0.0), c.velocity
    assert c.force == (0.0, 0.0, 0.0), c.force


def test_inactive_wc_is_skipped():
    """
    An inactive WC is not updated: position, velocity, force all stay
    at their starting values.

    Mutation caught: the active check dropped, so an inactive WC moves.
    """
    _ti_init()
    wc = WC(x=4.0, y=4.0, z=4.0, active=False)
    ctx = _build_pipeline([wc], _ConstantRule(1.0, 0.0, 0.0), max_steps=3)
    assert ctx.diag.errors == [], ctx.diag.errors

    c = ctx.data.require(WCState).centers[0]
    assert (c.x, c.y, c.z) == (4.0, 4.0, 4.0), (c.x, c.y, c.z)
    assert c.velocity == (0.0, 0.0, 0.0), c.velocity
    assert c.force == (0.0, 0.0, 0.0), c.force


# ======================================================================
# Constant rule: symplectic Euler integration
# ======================================================================


def test_constant_rule_integration_matches_analytic():
    """
    Constant acceleration (1, 0, 0), dt = 0.1, three steps. The
    integration is symplectic Euler:
        v_n = n * a * dt          = 0.1, 0.2, 0.3
        x_n = x0 + v_n * dt       = 4.01, 4.03, 4.06

    Mutation caught: a wrong dt power, or a velocity update that
    writes the new velocity to the position instead of the old.
    """
    _ti_init()
    wc = WC(x=4.0, y=4.0, z=4.0)
    ctx = _build_pipeline([wc], _ConstantRule(1.0, 0.0, 0.0), dt=0.1, max_steps=3)
    assert ctx.diag.errors == [], ctx.diag.errors

    c = ctx.data.require(WCState).centers[0]
    assert abs(c.x - (4.0 + 0.06)) < 1e-6, c.x
    assert abs(c.velocity[0] - 0.3) < 1e-6, c.velocity
    assert c.force == (1.0, 0.0, 0.0), c.force


def test_index_rule_passes_correct_index_per_center():
    """
    The processor passes each center its own index. A rule that
    returns (i, 0, 0) gives WC[0] zero drift and WC[1] a nonzero one.

    Mutation caught: the processor passes a fixed index, or passes
    the same index to every center.
    """
    _ti_init()
    wcs = [
        WC(x=3.0, y=4.0, z=4.0),
        WC(x=5.0, y=4.0, z=4.0),
    ]
    ctx = _build_pipeline(wcs, _IndexRule(), dt=0.1, max_steps=2)
    assert ctx.diag.errors == [], ctx.diag.errors

    c0, c1 = ctx.data.require(WCState).centers
    # Index 0 -> a = 0, no motion.
    assert c0.x == 3.0, c0.x
    assert c0.velocity == (0.0, 0.0, 0.0), c0.velocity
    # Index 1 -> a = 1, v = 0.2 after two steps, x = 5 + 0.03.
    assert abs(c1.x - (5.0 + 0.03)) < 1e-6, c1.x
    assert abs(c1.velocity[0] - 0.2) < 1e-6, c1.velocity


# ======================================================================
# GradRhoRule: direction and negative control
# ======================================================================


def _rho_with_minimum_at(i0, j0, k0, nx=8, ny=8, nz=8):
    """
    rho field with a sharp Gaussian deficit at (i0, j0, k0), 1.0
    everywhere else. Simple, smooth, monotonic toward the centre.
    """
    i, j, k = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz), indexing="ij")
    r2 = (i - i0) ** 2 + (j - j0) ** 2 + (k - k0) ** 2
    rho = 1.0 - 0.5 * np.exp(-r2 / 2.0)
    return rho.astype(np.float32)


def test_grad_rho_moves_wc_toward_deficit():
    """
    rho has its minimum at (3, 4, 4). A WC at (5, 4, 4) under
    GradRhoRule drifts in the -x direction: its x coordinate after
    several steps is less than its starting x.

    Mutation caught: the sign in GradRhoRule flipped. The companion
    test test_reverse_grad_rho_moves_away uses the flipped rule as a
    positive control for the direction claim.
    """
    _ti_init()
    rho = _rho_with_minimum_at(i0=3, j0=4, k0=4)
    wc = WC(x=5.0, y=4.0, z=4.0)
    ctx = _build_pipeline([wc], GradRhoRule(), rho_arr=rho, max_steps=5)
    assert ctx.diag.errors == [], ctx.diag.errors

    c = ctx.data.require(WCState).centers[0]
    assert c.x < 5.0, c.x
    assert c.velocity[0] < 0.0, c.velocity


def test_reverse_grad_rho_moves_away():
    """
    Negative control: +grad rho on the same arena moves the WC away
    from the deficit. This documents that the direction claim in the
    positive test is not a property of the arena; it is a property of
    the rule's sign.
    """
    _ti_init()
    rho = _rho_with_minimum_at(i0=3, j0=4, k0=4)
    wc = WC(x=5.0, y=4.0, z=4.0)
    ctx = _build_pipeline([wc], _ReverseGradRhoRule(), rho_arr=rho, max_steps=5)
    assert ctx.diag.errors == [], ctx.diag.errors

    c = ctx.data.require(WCState).centers[0]
    assert c.x > 5.0, c.x
    assert c.velocity[0] > 0.0, c.velocity


# ======================================================================
# Build-time: requires inspection and the WCState arm
# ======================================================================


def test_motion_processor_requires_all_three_features():
    """
    WCMotionProcessor.requires names all three of its features, and a
    pipeline missing the one the build-time graph can test in isolation
    (WCState, provided via initial_features rather than by a processor)
    fails at build time.

    The other two, WaveGrid and EMCDensityField, are pinned by direct
    inspection: AllocateWaveField provides both together, so a
    build-time arm that removes one removes the other.

    Mutation caught: any of the three dropped from requires.
    """
    _ti_init()

    proc = WCMotionProcessor(rule=NoOpRule())
    assert WaveGrid in proc.requires, proc.requires
    assert WCState in proc.requires, proc.requires
    assert EMCDensityField in proc.requires, proc.requires

    # Build-time arm for the feature the pipeline can be made to miss.
    class NoWCState(Pipeline):
        def __init__(self):
            super().__init__(error_policy=ErrorPolicy.FAIL_FAST)
            self.add(AllocateWaveField(nx=8, ny=8, nz=8, dx=1.0))
            self.add(WCMotionProcessor(rule=NoOpRule()))

    from ..runner import Runner
    from ..sinks import InMemorySink

    try:
        Runner({"session": InMemorySink()}).run(
            NoWCState(),
            name="missing_wcstate",
            params={},
            dt=0.1,
            max_steps=1,
        )
    except PipelineError as e:
        assert "WCState" in str(e), str(e)
        return
    raise AssertionError("expected PipelineError for missing WCState")


# ======================================================================
# Runner
# ======================================================================


def main() -> int:
    tests = [
        test_noop_rule_leaves_position_unchanged,
        test_inactive_wc_is_skipped,
        test_constant_rule_integration_matches_analytic,
        test_index_rule_passes_correct_index_per_center,
        test_grad_rho_moves_wc_toward_deficit,
        test_reverse_grad_rho_moves_away,
        test_motion_processor_requires_all_three_features,
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
