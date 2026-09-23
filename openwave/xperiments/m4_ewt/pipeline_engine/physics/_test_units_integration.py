"""
Integration test: UnitSystem passes through the pipeline engine.

Builds a minimal pipeline with a processor that requires UnitSystem,
runs it with a Runner, and verifies the processor reads the right
constants from the feature.
"""

from __future__ import annotations

import sys
from pathlib import Path

_THIS = Path(__file__).resolve()
_PROJECT_ROOT = _THIS.parents[5]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from openwave.xperiments.m4_ewt.pipeline_engine import (
    BaseProcessor,
    Pipeline,
    Runner,
    Stage,
)
from openwave.xperiments.m4_ewt.pipeline_engine.sinks import InMemorySink
from openwave.xperiments.m4_ewt.pipeline_engine.physics.units import (
    NaturalUnitSystem,
    OpenWaveUnitSystem,
    SIUnitSystem,
    UnitSystem,
)


class CapturedConstants:
    """Holder for the constants the dummy processor saw."""

    def __init__(self) -> None:
        self.c: float = 0.0
        self.wavelength: float = 0.0
        self.dx: float = 0.0
        self.dt: float = 0.0
        self.gamma: float = 0.0


class CaptureUnitSystem(BaseProcessor):
    name = "CaptureUnitSystem"
    stage = Stage.MEASURE
    order = 10
    requires = (UnitSystem,)
    provides = (CapturedConstants,)

    def setup(self, ctx) -> None:
        ctx.data.set(CapturedConstants())

    def process(self, ctx) -> None:
        units = ctx.data.require(UnitSystem)
        out = ctx.data.require(CapturedConstants)
        out.c = units.c
        out.wavelength = units.wavelength
        out.dx = units.dx
        out.dt = units.dt
        out.gamma = units.gamma


class DummyPipeline(Pipeline):
    def __init__(self) -> None:
        super().__init__(external_provides=(UnitSystem,))
        self.add(CaptureUnitSystem())


def _run_with(units: UnitSystem):
    sinks = {"session": InMemorySink()}
    runner = Runner(sinks)
    return runner.run(
        DummyPipeline(),
        name=f"units_integration_{type(units).__name__}",
        params={},
        dt=units.dt,
        max_steps=3,
        initial_features=[units],
    )


def test_natural_runs_through_pipeline():
    units = NaturalUnitSystem()
    ctx = _run_with(units)
    captured = ctx.data.require(CapturedConstants)
    assert captured.c == units.c
    assert captured.wavelength == units.wavelength
    assert captured.dx == units.dx
    assert captured.dt == units.dt
    assert captured.gamma == units.gamma
    assert not ctx.diag.errors


def test_openwave_runs_through_pipeline():
    units = OpenWaveUnitSystem()
    ctx = _run_with(units)
    captured = ctx.data.require(CapturedConstants)
    assert captured.c == units.c
    assert captured.gamma == units.gamma
    assert not ctx.diag.errors


def test_si_runs_through_pipeline():
    units = SIUnitSystem()
    ctx = _run_with(units)
    captured = ctx.data.require(CapturedConstants)
    assert captured.c == units.c
    assert not ctx.diag.errors


def main() -> int:
    tests = [
        test_natural_runs_through_pipeline,
        test_openwave_runs_through_pipeline,
        test_si_runs_through_pipeline,
    ]
    passed = 0
    for t in tests:
        try:
            t()
        except AssertionError as e:
            print(f"FAIL: {t.__name__}: {e}")
        except Exception as e:  # noqa: BLE001
            print(f"ERROR: {t.__name__}: {type(e).__name__}: {e}")
        else:
            print(f"PASS: {t.__name__}")
            passed += 1
    print(f"\n{passed}/{len(tests)} tests passed")
    return 0 if passed == len(tests) else 1


if __name__ == "__main__":
    sys.exit(main())
