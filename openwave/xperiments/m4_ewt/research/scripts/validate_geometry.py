#!/usr/bin/env python3
"""
Functional validation of the locked 1-3-6 geometry and the
wave-centre stability metric.

Usage:
    python research/scripts/validate_geometry.py
"""

import math
import sys
import types
from pathlib import Path

# Ensure the repo root is on sys.path so that `openwave` is importable
_REPO_ROOT = Path(__file__).resolve().parents[5]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from openwave.common import constants
from openwave.xperiments.m4_ewt.xparameters.utils import geometry as geom
from openwave.xperiments.m4_ewt.utils import instrumentation as inst

PASS = 0
FAIL = 1


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def dist(p, q):
    """Euclidean distance between two 3-tuples."""
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(p, q)))


def is_on_well(d, n, tol=1e-9):
    """Return True if *d* is *n* times lambda within *tol*."""
    return abs(d - n) <= tol


# ---------------------------------------------------------------------------
# Tests for tetrahedron_10_locked
# ---------------------------------------------------------------------------


def test_locked_inner_on_well():
    """Inner 3 WCs are exactly 1-lambda from the centre."""
    for edge in (1e-12, 2e-15, 4e-15):
        lam = constants.EWAVE_LENGTH / edge
        pts = geom.tetrahedron_10_locked(edge)
        centre = pts[0]
        for i in (1, 2, 3):
            d = dist(centre, pts[i]) / lam
            if not is_on_well(d, 1.0):
                print(f"  FAIL inner WC {i} at univ_edge={edge}: {d:.6f} lambda")
                return FAIL
    return PASS


def test_locked_outer_on_well():
    """Outer 6 WCs are exactly 2-lambda from the centre."""
    for edge in (1e-12, 2e-15, 4e-15):
        lam = constants.EWAVE_LENGTH / edge
        pts = geom.tetrahedron_10_locked(edge)
        centre = pts[0]
        for i in range(4, 10):
            d = dist(centre, pts[i]) / lam
            if not is_on_well(d, 2.0):
                print(f"  FAIL outer WC {i} at univ_edge={edge}: {d:.6f} lambda")
                return FAIL
    return PASS


def test_locked_scale_invariant():
    """Separations in units of lambda do not depend on univ_edge."""
    ref = None
    for edge in (1e-12, 2e-15, 4e-15):
        lam = constants.EWAVE_LENGTH / edge
        pts = geom.tetrahedron_10_locked(edge)
        d = sorted(
            [dist(pts[i], pts[j]) / lam for i in range(len(pts)) for j in range(i + 1, len(pts))]
        )
        if ref is None:
            ref = d
        else:
            if any(abs(a - b) > 1e-9 for a, b in zip(ref, d)):
                print(f"  FAIL scale invariance at univ_edge={edge}")
                return FAIL
    return PASS


# ---------------------------------------------------------------------------
# Tests for generate_positions_by_EWT_geometry_locked
# ---------------------------------------------------------------------------


def test_locked_dispatcher_k10():
    """K=10 returns the same as tetrahedron_10_locked."""
    for edge in (2e-15, 4e-15):
        pts1 = geom.tetrahedron_10_locked(edge)
        pts2 = geom.generate_positions_by_EWT_geometry_locked(edge, 10)
        for i, (p, q) in enumerate(zip(pts1, pts2)):
            if any(abs(a - b) > 1e-12 for a, b in zip(p, q)):
                print(f"  FAIL K=10 mismatch at index {i}, univ_edge={edge}")
                return FAIL
    return PASS


def test_locked_dispatcher_delegates():
    """For K != 10 the locked dispatcher returns bit-identical points to the legacy one."""
    for K in (2, 3, 4, 5, 6, 7, 8, 9):
        for edge in (2e-15, 4e-15):
            pts1 = geom.generate_positions_by_EWT_geometry(edge, K)
            pts2 = geom.generate_positions_by_EWT_geometry_locked(edge, K)
            for i, (p, q) in enumerate(zip(pts1, pts2)):
                if any(abs(a - b) > 1e-12 for a, b in zip(p, q)):
                    print(f"  FAIL K={K} mismatch at index {i}, univ_edge={edge}")
                    return FAIL
    return PASS


# ---------------------------------------------------------------------------
# Tests for log_stability_metrics (B2 fix)
# ---------------------------------------------------------------------------


class _MockWC:
    """Minimal wave-centre stub for testing the drift metric."""

    def __init__(self, positions, active=None):
        self.num_sources = len(positions)
        self.position_float = positions
        self.active = active or [1] * len(positions)


def test_stability_metric_survives_deactivation():
    """log_stability_metrics must not raise when a WC is deactivated."""
    # silence the JSON logger
    inst.json_logger = types.SimpleNamespace(log_timestep=lambda *a, **k: None)

    positions = [[i * 0.01, 0.0, 0.0] for i in range(10)]
    wc_all = _MockWC(positions)
    _, n = inst.log_stability_metrics(1, wc_all)
    if n != 10:
        print(f"  FAIL expected 10 active, got {n}")
        return FAIL

    # deactivate the last WC
    active9 = [1] * 10
    active9[9] = 0
    wc9 = _MockWC(positions, active9)
    try:
        _, n = inst.log_stability_metrics(2, wc9)
    except Exception as exc:
        print(f"  FAIL raised {type(exc).__name__}: {exc}")
        return FAIL
    if n != 9:
        print(f"  FAIL expected 9 active, got {n}")
        return FAIL
    return PASS


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------


def main():
    tests = [
        ("tetrahedron_10_locked: inner 3 on 1*lambda well", test_locked_inner_on_well),
        ("tetrahedron_10_locked: outer 6 on 2*lambda well", test_locked_outer_on_well),
        ("tetrahedron_10_locked: scale invariant", test_locked_scale_invariant),
        ("locked dispatcher: K=10 -> tetrahedron_10_locked", test_locked_dispatcher_k10),
        ("locked dispatcher: K!=10 delegates to legacy", test_locked_dispatcher_delegates),
        (
            "log_stability_metrics: survives WC deactivation",
            test_stability_metric_survives_deactivation,
        ),
    ]

    passed = 0
    for name, func in tests:
        print(f"  {name} ... ", end="", flush=True)
        res = func()
        if res == PASS:
            print("PASS")
            passed += 1
        else:
            print("FAIL")
            # func already printed details

    print(f"\n{passed}/{len(tests)} tests passed")
    sys.exit(0 if passed == len(tests) else 1)


if __name__ == "__main__":
    main()
