"""
Log sinks. All implement the ILogSink protocol.

Contains:
    ILogSink         -- the protocol
    JsonSessionSink  -- full session log, atomic write on flush/close
    LiveJsonSink     -- small rolling file, overwritten on each write
    InMemorySink     -- for tests and post-mortem inspection
"""
from __future__ import annotations

import json
from collections import deque
from pathlib import Path
from typing import Any, Mapping, Protocol, runtime_checkable


@runtime_checkable
class ILogSink(Protocol):
    def open(self, meta: Mapping[str, Any]) -> None: ...
    def write(self, step: int, payload: Mapping[str, Any]) -> None: ...
    def flush(self) -> None: ...
    def close(self) -> None: ...


class JsonSessionSink:
    """Full session log. Atomic write on flush/close. Read after the run."""

    def __init__(self, path: Path | str, flush_every: int = 0) -> None:
        self.path = Path(path)
        self.flush_every = flush_every
        self._meta: dict[str, Any] = {}
        self._records: list[dict[str, Any]] = []
        self._n = 0

    def open(self, meta: Mapping[str, Any]) -> None:
        self._meta = dict(meta)
        self.path.parent.mkdir(parents=True, exist_ok=True)

    def write(self, step: int, payload: Mapping[str, Any]) -> None:
        self._records.append({"step": step, **payload})
        self._n += 1
        if self.flush_every and self._n % self.flush_every == 0:
            self.flush()

    def flush(self) -> None:
        self._dump()

    def close(self) -> None:
        self._dump()

    def _dump(self) -> None:
        tmp = self.path.with_suffix(self.path.suffix + ".tmp")
        tmp.write_text(json.dumps({"meta": self._meta, "records": self._records}))
        tmp.replace(self.path)


class LiveJsonSink:
    """Small rolling file, overwritten on each write. Read by a monitor process."""

    def __init__(self, path: Path | str, history: int = 2000) -> None:
        self.path = Path(path)
        self.history = history
        self._meta: dict[str, Any] = {}
        self._records: deque[dict[str, Any]] = deque(maxlen=history)

    def open(self, meta: Mapping[str, Any]) -> None:
        self._meta = dict(meta)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        # Create the file immediately so a waiting monitor does not time out.
        self._dump()

    def write(self, step: int, payload: Mapping[str, Any]) -> None:
        self._records.append({"step": step, **payload})
        self._dump()

    def flush(self) -> None:
        pass

    def close(self) -> None:
        try:
            self.path.unlink()
        except FileNotFoundError:
            pass

    def _dump(self) -> None:
        tmp = self.path.with_suffix(self.path.suffix + ".tmp")
        tmp.write_text(json.dumps({"meta": self._meta, "records": list(self._records)}))
        tmp.replace(self.path)


class InMemorySink:
    """For tests and post-mortem inspection."""

    def __init__(self) -> None:
        self.meta: dict[str, Any] = {}
        self.records: list[dict[str, Any]] = []

    def open(self, meta: Mapping[str, Any]) -> None:
        self.meta = dict(meta)

    def write(self, step: int, payload: Mapping[str, Any]) -> None:
        self.records.append({"step": step, **payload})

    def flush(self) -> None:
        pass

    def close(self) -> None:
        pass