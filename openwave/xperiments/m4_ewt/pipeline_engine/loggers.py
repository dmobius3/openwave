"""
Generic logging processor.

LogProcessor is domain-agnostic: it reads a payload from the Context and
writes it to a named sink. What to log and how often is configured per instance.
"""
from __future__ import annotations

from typing import Any, Callable, Mapping

from .context import Context
from .pipeline import BaseProcessor, Stage


class LogProcessor(BaseProcessor):
    """
    Reads a payload from ctx and writes to a named sink.

    payload_fn(ctx) -> Mapping[str, Any]
    every: write only every N steps (>= 1)
    """

    stage = Stage.MEASURE
    provides: tuple[type, ...] = ()

    def __init__(
        self,
        name: str,
        sink_name: str,
        payload_fn: Callable[[Context], Mapping[str, Any]],
        every: int = 1,
        order: int = 0,
    ) -> None:
        self.name = name
        self.sink_name = sink_name
        self.payload_fn = payload_fn
        self.every = max(1, every)
        self.order = order

    def setup(self, ctx: Context) -> None:
        # Fail fast if the sink is missing.
        ctx.log.get(self.sink_name)

    def process(self, ctx: Context) -> None:
        if self.every > 1 and ctx.sim.step % self.every != 0:
            return
        sink = ctx.log.get(self.sink_name)
        sink.write(ctx.sim.step, self.payload_fn(ctx))