"""Append-only output ledger per protocol §13."""
import json
import time
import hashlib
import numpy as np


class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, (np.bool_,)):
            return bool(obj)
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)

LEDGER_PATH = "../ledger/OUTPUT_LEDGER.jsonl"

def _timestamp():
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())

def append_record(record):
    record["timestamp"] = _timestamp()
    line = json.dumps(record, separators=(",", ":"), cls=NumpyEncoder)
    with open(LEDGER_PATH, "a") as f:
        f.write(line + "\n")
    return line

def gate_record(gate_id, arena_id, rung, parent_status, mutation_status,
                measured_values, law_exception=False, extra=None):
    rec = {
        "type": "GATE",
        "gate_id": gate_id,
        "arena_id": arena_id,
        "rung": rung,
        "parent_status": parent_status,
        "mutation_status": mutation_status,
        "measured_values": measured_values,
        "law_exception": law_exception,
    }
    if extra:
        rec.update(extra)
    return append_record(rec)

def resource_record(gate_id, wall_clock_seconds, cumulative_seconds, ceiling_seconds=48*3600):
    return append_record({
        "type": "RESOURCE",
        "gate_id": gate_id,
        "wall_clock_seconds": wall_clock_seconds,
        "cumulative_seconds": cumulative_seconds,
        "ceiling_seconds": ceiling_seconds,
        "within_ceiling": cumulative_seconds < ceiling_seconds,
    })

def arm_allocation_record(gate_id, allocations):
    return append_record({
        "type": "ARM-ALLOCATION",
        "gate_id": gate_id,
        "allocations": allocations,
    })

def disposition_record(disposition_type, data):
    return append_record({
        "type": "DISPOSITION",
        "disposition": disposition_type,
        "data": data,
    })

def coverage_record(all_exercised, prohibited_calls, details):
    return append_record({
        "type": "COVERAGE",
        "all_exercised": all_exercised,
        "prohibited_calls": prohibited_calls,
        "details": details,
    })

def hash_ledger():
    h = hashlib.sha256()
    try:
        with open(LEDGER_PATH, "rb") as f:
            for chunk in iter(lambda: f.read(8192), b""):
                h.update(chunk)
        return h.hexdigest()
    except FileNotFoundError:
        return None
