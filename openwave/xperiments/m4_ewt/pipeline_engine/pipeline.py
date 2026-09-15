"""
Pipeline, stages, and the processor contract.

Contains:
    Stage             -- PRE_UPDATE, UPDATE, POST_UPDATE, MEASURE
    ErrorPolicy       -- FAIL_FAST, SOFT_STOP, CONTINUE
    PipelineStopSignal -- internal, used by SOFT_STOP
    PipelineError     -- raised on graph validation failure
    IProcessor        -- the processor protocol
    BaseProcessor     -- convenience base
    Pipeline          -- staged processor list with validation and error policy
"""

from __future__ import annotations

import time
from enum import Enum, auto
from typing import Callable, Protocol, runtime_checkable

from .context import Context

# ================================================================
# Stage and error policy
# ================================================================


class Stage(Enum):
    PRE_UPDATE = auto()
    UPDATE = auto()
    POST_UPDATE = auto()
    MEASURE = auto()


class ErrorPolicy(Enum):
    FAIL_FAST = auto()  # re-raise, stop the run immediately
    SOFT_STOP = auto()  # log, set should_stop, teardown, do not re-raise
    CONTINUE = auto()  # log, keep going with the next processor


class PipelineStopSignal(Exception):
    """Internal signal used by SOFT_STOP. Never escapes Pipeline.step()."""


class PipelineError(RuntimeError):
    """Raised on pipeline construction or validation failures."""


# ================================================================
# Processor
# ================================================================


@runtime_checkable
class IProcessor(Protocol):
    """
    Stateless processor. Configuration in __init__; runtime state in ctx.data.

    stage=None means lifecycle-only: only setup/teardown are called.
    stateless=False opts out of the stateless enforcement in Pipeline._safe.
    """

    name: str
    stage: Stage | None
    order: int
    requires: tuple[type, ...]
    provides: tuple[type, ...]
    stateless: bool

    def setup(self, ctx: Context) -> None: ...
    def process(self, ctx: Context) -> None: ...
    def teardown(self, ctx: Context) -> None: ...


class BaseProcessor:
    """Convenience base. Override what you need."""

    name: str = "BaseProcessor"
    stage: Stage | None = None
    order: int = 0
    requires: tuple[type, ...] = ()
    provides: tuple[type, ...] = ()
    stateless: bool = True

    def setup(self, ctx: Context) -> None:
        pass

    def process(self, ctx: Context) -> None:
        pass

    def teardown(self, ctx: Context) -> None:
        pass


# ================================================================
# Pipeline
# ================================================================


class Pipeline:
    """
    Staged processor list.

    build()  validates the requires/provides graph.
    setup()  runs all processors' setup(), with rollback on failure.
    step()   runs the four stages in order.
    teardown() runs all processors' teardown() in reverse, never masking errors.
    """

    def __init__(
        self,
        *,
        error_policy: ErrorPolicy = ErrorPolicy.SOFT_STOP,
        external_provides: tuple[type, ...] = (),
    ) -> None:
        self._all: list[IProcessor] = []
        self._by_stage: dict[Stage, list[IProcessor]] = {s: [] for s in Stage}
        self._built = False
        self._error_policy = error_policy
        self._external_provides = external_provides

    # --- Builder ---

    def add(self, p: IProcessor) -> "Pipeline":
        self._all.append(p)
        if p.stage is not None:
            self._by_stage[p.stage].append(p)
        self._built = False
        return self

    # --- Validation ---

    def build(self) -> None:
        self._all.sort(key=lambda p: p.order)
        for s in Stage:
            self._by_stage[s].sort(key=lambda p: p.order)
        self._validate()
        self._built = True

    def _validate(self) -> None:
        available: set[type] = set(self._external_provides)
        # Lifecycle-only processors first, in global order.
        for p in self._all:
            if p.stage is None:
                self._check_requires(p, available)
                available.update(p.provides)
        # Then per-stage, in stage order.
        for stage in Stage:
            for p in self._by_stage[stage]:
                self._check_requires(p, available)
                available.update(p.provides)

    @staticmethod
    def _check_requires(p: IProcessor, available: set[type]) -> None:
        for r in p.requires:
            if r not in available:
                raise PipelineError(
                    f"{p.name}: requires '{r.__name__}', " f"not provided by any earlier processor"
                )

    # --- Lifecycle ---

    def setup(self, ctx: Context) -> None:
        if not self._built:
            self.build()
        done: list[IProcessor] = []
        try:
            for p in self._all:
                self._safe(lambda p=p: p.setup(ctx), p, "setup", ctx)
                done.append(p)
        except Exception:
            # Roll back what has been set up.
            for p in reversed(done):
                try:
                    p.teardown(ctx)
                except Exception as e:
                    ctx.diag.error(f"{p.name}.teardown(after setup failure): {e}")
            raise

    def step(self, ctx: Context) -> None:
        stopped = False
        try:
            for stage in Stage:
                for p in self._by_stage[stage]:
                    self._safe(lambda p=p: p.process(ctx), p, "process", ctx)
                    if ctx.sim.should_stop:
                        stopped = True
                        break
                if stopped:
                    break
        except PipelineStopSignal:
            ctx.sim.should_stop = True
        # Advance the clock exactly once, regardless of how the step ended.
        ctx.sim.step += 1
        ctx.sim.t += ctx.sim.dt

    def teardown(self, ctx: Context) -> None:
        for p in reversed(self._all):
            try:
                p.teardown(ctx)
            except Exception as e:
                # Never let a teardown failure mask the others.
                ctx.diag.error(f"{p.name}.teardown: {type(e).__name__}: {e}")

    # --- Error handling and stateless enforcement ---

    def _safe(
        self,
        fn: Callable[[], None],
        p: IProcessor,
        op: str,
        ctx: Context,
    ) -> None:
        t0 = time.perf_counter()
        try:
            if op == "process" and getattr(p, "stateless", True):
                self._stateless_guard(fn, p)
            else:
                fn()
        except PipelineStopSignal:
            raise
        except Exception as e:
            msg = f"{p.name}.{op}: {type(e).__name__}: {e}"
            ctx.diag.error(msg)
            if self._error_policy is ErrorPolicy.FAIL_FAST:
                raise
            if self._error_policy is ErrorPolicy.CONTINUE and op == "process":
                return
            raise PipelineStopSignal(msg) from e
        finally:
            ctx.diag.time(f"{p.name}.{op}", time.perf_counter() - t0)

    @staticmethod
    def _stateless_guard(fn: Callable[[], None], p: IProcessor) -> None:
        before = {k: id(v) for k, v in p.__dict__.items()}
        fn()
        after = {k: id(v) for k, v in p.__dict__.items()}
        if before != after:
            changed = sorted(k for k in before | after if before.get(k) != after.get(k))
            raise PipelineError(
                f"{p.name}: mutated instance fields in process(): {changed}. "
                f"Runtime state must live in ctx.data. "
                f"Set `stateless = False` if memoization is required."
            )
