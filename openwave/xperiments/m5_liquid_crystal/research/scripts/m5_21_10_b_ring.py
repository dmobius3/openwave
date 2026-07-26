"""M5.21.10 Read 3: the point-vs-charged-ring tie-breaker at next
rigor.

The M5.21.2 tie was instrument-limited (ring -3.7% under fwd vs +23%
under the 2h re-read, pre-2b instrument). Next rigor here: all three
seeds (A hedgehog; R aring=4; R aring=6) relaxed on the CERTIFIED
symmetrized instrument (T2, sym, n=32 pinned, maxit=16000), then

    read 1  native (sym) converged energies + stop reasons
    read 2  cross-stencil re-read: the SAME converged endpoint's
            energy re-evaluated under fwd and 2h stencils (no
            re-relax): if the ordering is stencil-stable, the tie
            verdict is instrument-robust
    read 3  state distance: rel Frobenius distance between endpoints
            + core spectra: same-basin (the 2b merge) or distinct?

Consumes m5_21_6_end_t32_{A,R4,R6}.npz; writes m5_21_10_ring.json.
"""
from __future__ import annotations

import importlib.util
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")

_spec = importlib.util.spec_from_file_location(
    "ins", os.path.join(HERE, "m5_21_2b_a_instrument.py"))
INS = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(INS)


def spec_core(M, cfg, r_c=6.0):
    X, Y, Z = INS.coords(cfg["n"], cfg["h"])
    r = np.sqrt(X * X + Y * Y + Z * Z)
    lam = np.linalg.eigvalsh(M)
    return np.mean(lam[r < r_c], axis=0)


def main():
    tags = ("t32_A", "t32_R4", "t32_R6")
    out = {"tags": list(tags), "rows": {}}
    Ms = {}
    for tag in tags:
        z = np.load(os.path.join(DATA, f"m5_21_6_end_{tag}.npz"))
        M = z["M"].astype(float)
        Ms[tag] = M
        n, h = int(z["n"]), float(z["h"])
        row_meta = json.load(open(os.path.join(
            DATA, f"m5_21_6_row_{tag}.json")))
        r = {"stop": row_meta["stop"], "maxit": row_meta["maxit"],
             "E_end_recorded": row_meta["E_end"]}
        for st in ("sym", "fwd", "2h"):
            cfg = INS.base_cfg(term="T2", n=n, L=n * h, bc="pinned",
                               stencil=st)
            e_u, e_d, e_v = INS.e_parts(M, cfg)
            r[f"E_{st}"] = float(e_u + e_d + e_v)
        cfg = INS.base_cfg(term="T2", n=n, L=n * h, bc="pinned")
        r["spec_core"] = [float(x) for x in spec_core(M, cfg)]
        out["rows"][tag] = r
        print(tag, json.dumps(r))
    # pairwise state distances (raw, no alignment: pinned frame)
    out["dist"] = {}
    for a in tags:
        for b in tags:
            if a < b:
                d = float(np.linalg.norm(Ms[a] - Ms[b]) /
                          max(np.linalg.norm(Ms[a]), 1e-300))
                out["dist"][f"{a}__{b}"] = d
    # verdict scaffolding: orderings per stencil
    for st in ("sym", "fwd", "2h"):
        es = {t: out["rows"][t][f"E_{st}"] for t in tags}
        order = sorted(tags, key=lambda t: es[t])
        rel = {t: (es[t] - es[order[0]]) /
               max(abs(es[order[0]]), 1e-300) for t in tags}
        out[f"order_{st}"] = {"order": order,
                              "rel_above_min": rel}
        print(st, "order", order,
              {t: f"{rel[t]:.2e}" for t in tags})
    with open(os.path.join(DATA, "m5_21_10_ring.json"), "w") as f:
        json.dump(out, f, indent=1)


if __name__ == "__main__":
    main()
