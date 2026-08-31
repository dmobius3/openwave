"""Output ledger utilities: append-only JSON records per § 13."""
import json
import time
import os

LEDGER_PATH = 'ledger/OUTPUT_LEDGER.jsonl'

_start_time = None
_cumulative_seconds = 0.0


def init_clock():
    global _start_time
    _start_time = time.time()


def _wall_since_init():
    if _start_time is None:
        return 0.0
    return time.time() - _start_time


class _NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        import numpy as np
        if isinstance(obj, (np.bool_,)):
            return bool(obj)
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)


def append_record(record):
    """Append one JSON record to the output ledger."""
    with open(LEDGER_PATH, 'a') as f:
        f.write(json.dumps(record, ensure_ascii=False, cls=_NumpyEncoder) + '\n')


def gate_record(gate_id, arena_id, rung, parent_status, mutation_status,
                measured_values, wall_clock_seconds, law_exception=False,
                **extra):
    """Write a GATE record per § 13."""
    global _cumulative_seconds
    _cumulative_seconds += wall_clock_seconds
    rec = {
        'type': 'GATE',
        'gate_id': gate_id,
        'arena_id': arena_id,
        'rung': rung,
        'parent_status': parent_status,
        'mutation_status': mutation_status,
        'measured_values': measured_values,
        'law_exception': law_exception,
        'wall_clock_seconds': round(wall_clock_seconds, 3),
        'cumulative_seconds': round(_cumulative_seconds, 3),
    }
    rec.update(extra)
    append_record(rec)
    return rec


def resource_record(gate_id, wall_clock_seconds):
    """Write a RESOURCE record per § 13."""
    global _cumulative_seconds
    _cumulative_seconds += wall_clock_seconds
    rec = {
        'type': 'RESOURCE',
        'gate_id': gate_id,
        'wall_clock_seconds': round(wall_clock_seconds, 3),
        'cumulative_seconds': round(_cumulative_seconds, 3),
    }
    append_record(rec)
    return rec


def arm_allocation_record(allocations):
    """Write an ARM-ALLOCATION record per § 13."""
    rec = {
        'type': 'ARM-ALLOCATION',
        'allocations': allocations,
    }
    append_record(rec)
    return rec


def disposition_record(disposition_type, **fields):
    """Write a DISPOSITION record per § 13."""
    rec = {'type': 'DISPOSITION', 'disposition': disposition_type}
    rec.update(fields)
    append_record(rec)
    return rec


def coverage_record(registry_coverage, prohibited_check):
    """Write a COVERAGE record for gate 10."""
    rec = {
        'type': 'COVERAGE',
        'registry_coverage': registry_coverage,
        'prohibited_check': prohibited_check,
    }
    append_record(rec)
    return rec
