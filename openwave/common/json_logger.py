"""
COMMON JSON LOGGER - model-agnostic instrumentation sink.

Every model (M3, M4, ...) calls init_session() once with a flat
metadata dict, then log_timestep() for each data row.  The logger
writes one JSON file per session, named:

    <model>_<xperiment>[_K<k>].json

where the optional _K<k> suffix is added when the metadata contains
a key "K" (wave-centre count).
"""

import json
from pathlib import Path
from datetime import datetime
from typing import Any, Dict, List, Optional

# ---------------------------------------------------------------------------
# Internal state
# ---------------------------------------------------------------------------
_buffer: List[Dict[str, Any]] = []
_meta: Dict[str, Any] = {}
_filename: Optional[str] = None
_data_dir: Path = Path("data")
_flush_interval: int = 100

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def init_session(
    metadata: Dict[str, Any],
    data_dir: Path = Path("data"),
    flush_interval: int = 100,
) -> None:
    """
    Start a new logging session.

    Parameters
    ----------
    metadata : dict
        Flat dictionary with session metadata. Required keys:
            "model"       - e.g. "M3", "M4"
            "xperiment"   - xparameter file name without .py
        Optional keys:
            "K"           - wave-centre count (used for filename)
            "engine"      - sub-dict with engine parameters
            ... anything else you want recorded.
    data_dir : Path
        Directory where JSON files will be stored.
    flush_interval : int
        Auto-flush after this many rows.
    """
    global _buffer, _meta, _filename, _data_dir, _flush_interval

    _buffer = []
    _meta = metadata.copy()
    _meta.setdefault("timestamp", datetime.now().isoformat())
    _data_dir = Path(data_dir)
    _flush_interval = flush_interval

    # Build filename with timestamp to avoid overwriting previous runs
    model = _meta.get("model", "unknown")
    xper = _meta.get("xperiment", "unknown")
    k = _meta.get("K")
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    if k is not None:
        _filename = f"{model}_{xper}_K{k}_{ts}.json"
    else:
        _filename = f"{model}_{xper}_{ts}.json"


def log_timestep(data: Dict[str, Any]) -> None:
    """
    Append one timestep record.  Must contain at least "timestep".
    """
    global _buffer
    _buffer.append(data)
    if len(_buffer) >= _flush_interval:
        _flush()


def finalize() -> None:
    """Flush any remaining records.  Call at end of simulation."""
    _flush()
    _buffer = []
    _meta = {}


# ---------------------------------------------------------------------------
# Internal
# ---------------------------------------------------------------------------

def _flush() -> None:
    global _buffer, _meta, _filename, _data_dir
    if not _buffer or _filename is None:
        return

    _data_dir.mkdir(parents=True, exist_ok=True)
    log_path = _data_dir / _filename

    # Merge with existing file if it exists (resume after crash)
    if log_path.exists():
        try:
            with open(log_path, "r") as fh:
                doc = json.load(fh)
            doc["data"].extend(_buffer)
        except (json.JSONDecodeError, OSError):
            # File is empty or corrupted – overwrite it
            doc = {"metadata": _meta, "data": list(_buffer)}
    else:
        doc = {"metadata": _meta, "data": list(_buffer)}

    with open(log_path, "w") as fh:
        json.dump(doc, fh, indent=2)

    print(f"[json_logger] Flushed {len(_buffer)} records -> {log_path}")
    _buffer.clear()