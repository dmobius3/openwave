"""Stage 2, informational: sensitivity of the T4 protocol to the reading of "uniform refinement".

Both rooms refine by midpoint insertion (the level-0 graded nodes are kept).  The other reading
regrades every level: y-nodes at pi/2 -+ (pi/2)(k/N)^2 with N = 16 * 2^level.  This script runs
Room A's assembly (copy in lib_a/, loaded by s2_audit.load_rooms) with that node set, applies the
same estimates, and compares with my stage-1 bottoms.  Writes s2_regrade.json.  Not a check:
nothing here prints PASS.
"""
import json
import math
import os
import sys

import numpy as np
import scipy.sparse.linalg as spla

HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, "lib_a", "t4_fem.py")).read()
src = src[:src.index("results = {}")]
A = {"__name__": "room_a_t4", "__file__": os.path.join(HERE, "lib_a", "t4_fem.py")}
sys.argv = ["t4_fem.py"]
exec(compile(src, "lib_a/t4_fem.py", "exec"), A)


def ynodes_regraded(level):
    N = 16 * 2**level
    k = np.arange(N + 1)
    side = (math.pi / 2) * (k / N) ** 2
    return np.concatenate([math.pi / 2 - side[::-1], math.pi / 2 + side[1:]])


A["ynodes"] = ynodes_regraded
stage1 = json.load(open(os.path.join(HERE, "stage1.json")))["v2"]
WIDTHS = {"1/4": 0.25, "1/2": 0.5, "1": 1.0, "7/5": 1.4, "3/2": 1.5, "pi/2": math.pi / 2, "17/10": 1.7}
out = {}
for lab, W in WIDTHS.items():
    lv = []
    for lev in range(5):
        K, M, *_ = A["assemble"](W, lev)
        lv.append(float(np.min(spla.eigsh(K, k=1, M=M, sigma=0.0, which="LM", return_eigenvectors=False))))
    a, b, c = lv[2], lv[3], lv[4]
    p = math.log2((a - b) / (b - c))
    ex = c + (c - b) / (2**p - 1)
    mine = stage1[lab]["bottom"]
    out[lab] = {"bottom_levels": lv, "p": p, "extrap": ex, "err": abs(ex - c), "rel_vs_mine": (ex - mine) / mine}
    print(f"W={lab:6s} regraded: levels " + " ".join(f"{v:.10f}" for v in lv)
          + f"  p={p:.6f} extrap={ex:.12f} err={abs(ex-c):.3e} rel vs mine {(ex-mine)/mine:+.2e}", flush=True)
with open(os.path.join(HERE, "s2_regrade.json"), "w") as fh:
    json.dump(out, fh, indent=1)
