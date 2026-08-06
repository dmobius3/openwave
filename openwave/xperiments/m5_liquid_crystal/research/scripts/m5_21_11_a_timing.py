"""M5.21.11 D2 timing probes (affordability measurement ONLY).

Runs 1-2 SHORT relaxations of the census instrument of record
(m5_21_2b_a_instrument.py, T2 / sym / pinned) purely to measure
wall-clock cost per FIRE iteration at the grid sizes the pre-registered
ladder will use, so the rung budget in m5_21_11_framework.md § 2 is
MEASURED, not assumed.

The endpoints are NOT physics: maxit is far below convergence depth,
no energy from these runs enters any fit, and the tags are
t11timing_* so they can never be confused with ladder rungs.

MODES: run | report
Out: ../data/m5_21_11_timing.json
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT = os.path.join(DATA, "m5_21_11_timing.json")

# (label, args, maxit) -- exactly the instrument-of-record settings the
# ladder will freeze (T2, sym, pinned, w2 = the recorded census match).
PROBES = [
    ("n48_d0.1", ["seed=A", "term=T2", "stencil=sym", "eps=0", "n=48",
                  "delta=0.1", "maxit=400", "w2=0.002758100",
                  "tag=t11timing_n48"], 400),
    ("n64_d0.3", ["seed=A", "term=T2", "stencil=sym", "eps=0", "n=64",
                  "delta=0.3", "maxit=200", "w2=0.002758100",
                  "tag=t11timing_n64"], 200),
]

RUNG_MAXIT = 12000          # the census production depth (m5_21_2b c48 runs)


def run() -> None:
    rows = {}
    for label, args, maxit in PROBES:
        t0 = time.time()
        subprocess.run(
            [sys.executable, os.path.join(HERE, "m5_21_2b_a_instrument.py"),
             "relax", *args],
            check=True,
        )
        wall = time.time() - t0
        rows[label] = {
            "args": " ".join(args),
            "maxit": maxit,
            "wall_s": round(wall, 2),
            "s_per_iter": round(wall / maxit, 4),
            "proj_rung_12k_min": round(wall / maxit * RUNG_MAXIT / 60.0, 1),
        }
        print(f"[timing] {label}: {wall:.1f}s for {maxit} it "
              f"-> {rows[label]['proj_rung_12k_min']} min per 12k rung")
    with open(OUT, "w") as f:
        json.dump(rows, f, indent=2)
    print(f"[timing] wrote {OUT}")


def report() -> None:
    with open(OUT) as f:
        rows = json.load(f)
    n48 = rows["n48_d0.1"]["s_per_iter"]
    n64 = rows["n64_d0.3"]["s_per_iter"]
    print(json.dumps(rows, indent=2))
    print(f"refine cost ratio n64/n48 (per iter) = {n64 / n48:.2f} "
          f"(volume scaling predicts (64/48)^3 = {(64 / 48) ** 3:.2f})")


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "run"
    {"run": run, "report": report}[mode]()
