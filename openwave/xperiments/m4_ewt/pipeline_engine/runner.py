"""
The default runner: owns sinks, builds the Context, runs the pipeline.

Contains:
    IRunner  -- the runner protocol
    Runner   -- the default implementation

Domain-agnostic: no imports of physics, no imports of Taichi.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping, Protocol, runtime_checkable

from .context import (
    Context,
    Diagnostics,
    FeatureBag,
    LogContext,
    Params,
    RunContext,
    SimContext,
)
from .pipeline import Pipeline
from .sinks import ILogSink


@runtime_checkable
class IRunner(Protocol):
    """
    Runner interface. Implementations decide how to run a pipeline.
    The default Runner runs synchronously, in-process.
    """

    def run(
        self,
        pipeline: Pipeline,
        *,
        name: str,
        params: Mapping[str, Any] | None = None,
        output_dir: Path | str = "out",
        dt: float = 1.0,
        max_steps: int,
    ) -> Context: ...


class Runner:
    """Default runner. Owns sinks and context; runs the pipeline synchronously."""

    def __init__(
        self,
        sinks: Mapping[str, ILogSink],
        *,
        check_stateless: bool = False,
    ) -> None:
        self._sinks: dict[str, ILogSink] = dict(sinks)
        self._check_stateless = check_stateless

    def run(
        self,
        pipeline: Pipeline,
        *,
        name: str,
        params: Mapping[str, Any] | None = None,
        output_dir: Path | str = "out",
        dt: float = 1.0,
        max_steps: int,
        initial_features: list[object] | None = None,
    ) -> Context:
        run = RunContext(name=name, output_dir=Path(output_dir))
        ctx = Context(
            run=run,
            params=Params(params),
            sim=SimContext(dt=dt),
            data=FeatureBag(),
            log=LogContext(sinks=dict(self._sinks)),
            diag=Diagnostics(),
        )

        # Seed any features provided by the caller (e.g. WCState).
        if initial_features:
            for f in initial_features:
                ctx.data.set(f)

        meta: Mapping[str, Any] = {
            "run": run.name,
            "params": dict(params or {}),
            "output_dir": str(run.output_dir),
        }

        opened: list[ILogSink] = []
        try:
            for sink in self._sinks.values():
                sink.open(meta)
                opened.append(sink)

            pipeline.setup(ctx)
            for _ in range(max_steps):
                if ctx.sim.should_stop:
                    break
                pipeline.step(ctx)
        finally:
            try:
                pipeline.teardown(ctx)
            finally:
                for sink in reversed(opened):
                    try:
                        sink.flush()
                    except Exception as e:
                        ctx.diag.error(f"sink.flush: {type(e).__name__}: {e}")
                    try:
                        sink.close()
                    except Exception as e:
                        ctx.diag.error(f"sink.close: {type(e).__name__}: {e}")

        return ctx
