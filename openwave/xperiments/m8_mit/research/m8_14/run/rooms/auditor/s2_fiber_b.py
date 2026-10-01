"""Stage 2: Room B's claim that on triangles with an edge on the fiber the free node's w-gradient
is "exactly zero", so the large 1/|cos y| integral is multiplied by 0.

Room B computes gradients from np.linalg.inv of the 3x3 barycentric matrix, so "exactly" is a
floating-point question.  For every fiber-edge triangle at levels 0..4 (W = 17/10 and 1/4) this
records |gw_free| and the product |gw_free|^2 * Iic, the stiffness entry that would leak in, and
compares it with the K diagonal.  Writes s2_fiber_b.json.

Check F1 PASS iff |gw_free| * h_w < 1e-12 (zero to roundoff, h_w the w-cell; a genuine gradient
is O(1/h_w)) and max leak / min K-diagonal < 1e-12.  Planted: gw_free + 1e-6 (a gradient that is
not zero) must FAIL.  (A first version tested only the leak ratio; its planted twin did not fire,
because a 1e-6 gradient leaks only 5e-13 of a diagonal entry; the criterion was restated as above.)
"""
import json
import math
import os
import sys
import importlib.util

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("room_b_t4", os.path.join(HERE, "lib_b", "t4_fem.py"))
B = importlib.util.module_from_spec(spec)
sys.argv = ["t4_fem.py"]
spec.loader.exec_module(B)
out, lines = {}, []


def run(W, planted=False):
    ys, nw = B.base_mesh()
    ws = np.linspace(-W, W, nw + 1)
    res = []
    for lev in range(5):
        ifib = int(np.argmin(np.abs(ys - math.pi / 2)))
        worst_g, worst_leak = 0.0, 0.0
        for i in (ifib - 1, ifib):
            for j in range(len(ws) - 1):
                nodes = [(i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)]
                for tri in ([0, 1, 2], [0, 2, 3]):
                    nd = [nodes[t] for t in tri]
                    on = [a for a, _ in nd].count(ifib)
                    if on != 2:
                        continue
                    P = np.array([[ys[a], ws[b]] for a, b in nd])
                    Ic, Iic, M, C = B.tri_integrals(P, B.weight)
                    q = [t for t, (a, _) in enumerate(nd) if a != ifib][0]
                    g = abs(C[2][q]) + (1e-6 if planted else 0.0)
                    worst_g = max(worst_g, g * (ws[1] - ws[0]))
                    worst_leak = max(worst_leak, g * g * Iic)
        K, M, idx, sgn = B.assemble(ys, ws) if lev <= 2 else (None, None, None, None)
        kd = float(K.diagonal().min()) if K is not None else float("nan")
        res.append({"level": lev, "max_abs_gw_free_times_hw": worst_g, "max_leak": worst_leak, "min_Kdiag": kd})
        ys, ws = B.refine(ys), B.refine(ws)
        if lev >= 2 and planted:
            break
    kmin = min(r["min_Kdiag"] for r in res if r["min_Kdiag"] == r["min_Kdiag"])
    ratio = max(r["max_leak"] for r in res) / kmin
    gmax = max(r["max_abs_gw_free_times_hw"] for r in res)
    tag = "PASS" if ratio < 1e-12 and gmax < 1e-12 else "FAIL"
    s = f"{'[planted] ' if planted else ''}{tag}  F1 room b W={W}: max |gw_free|*h_w {gmax:.1e}, " \
        f"max leak {max(r['max_leak'] for r in res):.1e}, ratio to min K-diagonal {ratio:.1e}"
    print(s, flush=True)
    lines.append(s)
    return res, (ratio, gmax)


for W in (1.7, 0.25):
    out[str(W)] = run(W)
r, (ratio, gmax) = run(1.7, planted=True)
out["planted_fired"] = bool(ratio >= 1e-12 or gmax >= 1e-12)
out["lines"] = lines
with open(os.path.join(HERE, "s2_fiber_b.json"), "w") as fh:
    json.dump(out, fh, indent=1)
