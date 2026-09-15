"""
Test build_wc_state: pure Python, no simulation.

Run from anywhere:
    python openwave/xperiments/m4_ewt/pipeline_engine/physics/_test_wc_factory.py
"""

import math
import sys
from pathlib import Path

# Bootstrap: add m4_ewt/ to sys.path so 'pipeline_engine' is importable
# as a proper package, and 'physics' is its subpackage.
_THIS = Path(__file__).resolve()
_ENGINE_DIR = _THIS.parents[1]  # .../pipeline_engine
_M4EWT_DIR = _ENGINE_DIR.parent  # .../m4_ewt
sys.path.insert(0, str(_M4EWT_DIR))

from pipeline_engine.physics.wc_factory import build_wc_state
from pipeline_engine.physics.wc_types import WCState


def report(name, wc_state):
    print(f"--- {name} ---")
    print(f"  K       = {wc_state.K}")
    print(f"  active  = {wc_state.active_count}")
    for i, wc in enumerate(wc_state.centers):
        print(
            f"  WC[{i:2d}]: pos=({wc.x:7.3f}, {wc.y:7.3f}, {wc.z:7.3f}) "
            f"phase={wc.phase:6.3f} active={wc.active}"
        )
    n = wc_state.K
    if n >= 2:
        dmin = float("inf")
        dmax = 0.0
        for i in range(n):
            for j in range(i + 1, n):
                a = wc_state.centers[i]
                b = wc_state.centers[j]
                d = math.sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2)
                dmin = min(dmin, d)
                dmax = max(dmax, d)
        print(f"  pair dist: min={dmin:.3f} max={dmax:.3f}")
    print()


def main():
    grid = 64
    wavelength = 12.0

    cases = [
        (1, "single"),
        (2, "pair"),
        (3, "line"),
        (5, "line"),
        (5, "golden"),
        (10, "golden"),
        (10, "tetrahedron_10_locked"),
    ]

    for K, geometry in cases:
        state = build_wc_state(
            K=K,
            geometry=geometry,
            nx=grid,
            ny=grid,
            nz=grid,
            wavelength=wavelength,
        )
        report(f"K={K}, geometry='{geometry}'", state)

    print("--- error cases ---")
    unexpected = 0
    for K, geometry in [(1, "pair"), (2, "single"), (1, "line"), (11, "tetrahedron_10_locked")]:
        try:
            build_wc_state(K=K, geometry=geometry, nx=grid, ny=grid, nz=grid)
            print(f"  K={K}, geometry='{geometry}': NO ERROR (unexpected)")
            unexpected += 1
        except ValueError as e:
            print(f"  K={K}, geometry='{geometry}': ValueError -> {e}")

    if unexpected:
        print(f"FAILED: {unexpected} error case(s) did not raise ValueError")
        sys.exit(1)
    print("all error cases raised ValueError")


if __name__ == "__main__":
    main()
