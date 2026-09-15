"""
Smoke test for the pipeline engine. No Taichi, no physics.

Run:
    python -m openwave.xperiments.m4_ewt.pipeline_engine._smoke_test

or, from inside the folder:
    python _smoke_test.py
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping

# Local imports. If run as a script from inside the folder, use relative paths.
try:
    from . import (
        BaseProcessor,
        Context,
        ErrorPolicy,
        InMemorySink,
        JsonSessionSink,
        LiveJsonSink,
        LogProcessor,
        Pipeline,
        Runner,
        Stage,
    )
except ImportError:
    # Running as a plain script from inside the folder.
    from context import (  # type: ignore[no-redef]
        Context,
    )
    from pipeline import BaseProcessor, ErrorPolicy, Pipeline, Stage  # type: ignore[no-redef]
    from sinks import InMemorySink, JsonSessionSink, LiveJsonSink  # type: ignore[no-redef]
    from loggers import LogProcessor  # type: ignore[no-redef]
    from runner import Runner  # type: ignore[no-redef]


# ================================================================
# Data features (runtime state lives in ctx.data)
# ================================================================

@dataclass
class Counter:
    value: int = 0


@dataclass
class History:
    values: list[float] = field(default_factory=list)


# ================================================================
# Processors
# ================================================================

class AllocateState(BaseProcessor):
    name = "AllocateState"
    stage = None
    order = 10
    provides = (Counter, History)

    def setup(self, ctx: Context) -> None:
        ctx.data.set(Counter())
        ctx.data.set(History())


class Increment(BaseProcessor):
    name = "Increment"
    stage = Stage.UPDATE
    order = 10
    requires = (Counter,)

    def process(self, ctx: Context) -> None:
        ctx.data.require(Counter).value += 1


class RecordValue(BaseProcessor):
    name = "RecordValue"
    stage = Stage.MEASURE
    order = 10
    requires = (Counter, History)

    def process(self, ctx: Context) -> None:
        c = ctx.data.require(Counter)
        h = ctx.data.require(History)
        h.values.append(float(c.value % 7))


class BadProcessor(BaseProcessor):
    """Forces a failure to demonstrate ErrorPolicy.SOFT_STOP."""
    name = "BadProcessor"
    stage = Stage.UPDATE
    order = 90

    def process(self, ctx: Context) -> None:
        raise RuntimeError("intentional failure")


# ================================================================
# Pipeline + payload
# ================================================================

def _payload(ctx: Context) -> Mapping[str, Any]:
    c = ctx.data.require(Counter)
    h = ctx.data.require(History)
    return {
        "counter": c.value,
        "last": h.values[-1] if h.values else 0.0,
        "history_len": len(h.values),
    }


class DemoPipeline(Pipeline):
    def __init__(self, *, inject_failure: bool = False) -> None:
        super().__init__(error_policy=ErrorPolicy.SOFT_STOP)
        self.add(AllocateState())
        self.add(Increment())
        if inject_failure:
            self.add(BadProcessor())
        self.add(RecordValue())
        self.add(LogProcessor("SessionLog", "session", _payload, every=10, order=50))
        self.add(LogProcessor("LiveLog",    "live",    _payload, every=1,  order=60))


# ================================================================
# Runner helper
# ================================================================

def _print_summary(label: str, ctx: Context) -> None:
    print(f"--- {label} ---")
    print(f"steps simulated : {ctx.sim.step}")
    print(f"counter         : {ctx.data.require(Counter).value}")
    print(f"errors          : {ctx.diag.errors}")
    print(f"warnings        : {ctx.diag.warnings}")
    slowest = sorted(
        ((k, sum(v) / len(v) * 1e6) for k, v in ctx.diag.timings.items()),
        key=lambda x: -x[1],
    )[:5]
    print("slowest operations (avg us):")
    for name, us in slowest:
        print(f"  {name:<28} {us:>10.2f}")
    print()


def _run_demo(
    *,
    inject_failure: bool,
    use_memory: bool,
    out_dir: Path,
) -> Context:
    sinks: dict[str, Any] = {}
    if use_memory:
        sinks["session"] = InMemorySink()
        sinks["live"] = InMemorySink()
    else:
        sinks["session"] = JsonSessionSink(out_dir / "session.json", flush_every=100)
        sinks["live"] = LiveJsonSink(out_dir / "live.json", history=200)

    pipeline = DemoPipeline(inject_failure=inject_failure)
    runner = Runner(sinks, check_stateless=True)

    return runner.run(
        pipeline,
        name="demo",
        params={"note": "engine-only demo", "amplitude": 1.0},
        output_dir=out_dir,
        dt=1.0,
        max_steps=200,
    )


# ================================================================
# Entry point
# ================================================================

def main() -> None:
    out_dir = Path("out_engine_smoke")

    print("=" * 64)
    print("Pipeline engine smoke test")
    print("=" * 64)
    print()

    # 1) Normal run, in-memory sinks
    ctx = _run_demo(inject_failure=False, use_memory=True, out_dir=out_dir)
    _print_summary("normal run (in-memory)", ctx)
    assert ctx.sim.step == 200, f"expected 200 steps, got {ctx.sim.step}"
    assert ctx.data.require(Counter).value == 200
    assert ctx.diag.errors == []

    # 2) Run with injected failure, in-memory sinks
    ctx = _run_demo(inject_failure=True, use_memory=True, out_dir=out_dir)
    _print_summary("failure run (in-memory)", ctx)
    assert ctx.sim.step == 1, f"expected 1 step, got {ctx.sim.step}"
    assert any("BadProcessor" in e for e in ctx.diag.errors)

    # 3) Normal run, file sinks
    ctx = _run_demo(inject_failure=False, use_memory=False, out_dir=out_dir)
    _print_summary("normal run (file sinks)", ctx)
    session_path = out_dir / "session.json"
    assert session_path.exists(), f"expected session log at {session_path}"
    print(f"session log written to: {session_path}")
    print()

    print("=" * 64)
    print("SMOKE TEST PASSED")
    print("=" * 64)


if __name__ == "__main__":
    main()