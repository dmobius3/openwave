"""M5.21.5 G2: the matched T2 descent family for the mu ladder.

Fresh A-hedgehog descents at n = 24 (L = 36) and n = 48 (L = 72),
h = 1.5, EXACTLY the t32_A recipe otherwise (T2, sym, pinned,
maxit = 16000, delta = 0.3; t32_A row of record: E_end = 4.767567,
stop = f_tol, wall 724 s). The n = 32 rung is m5_21_6_end_t32_A.npz
(on disk, M5.21.10 Read 3). Writer is m5_21_5-named so nothing of
the m5_21_6 family is touched.

Usage: python3 m5_21_5_b_ladder.py n=24
       python3 m5_21_5_b_ladder.py n=48
Out:   data/m5_21_5_end_t<n>_A.npz + data/m5_21_5_row_t<n>_A.json
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")

_spec = importlib.util.spec_from_file_location(
    "ins", os.path.join(HERE, "m5_21_2b_a_instrument.py"))
INS = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(INS)


def run(n):
    cfg = INS.base_cfg(term="T2", n=n, L=1.5 * n, bc="pinned",
                       seed="A", stencil="sym", maxit=16000)
    tag = f"t{n}_A"
    M0 = INS.make_seed(cfg)
    free = ~INS.pin_shell(n, cfg["h"])
    e0 = INS.e_parts(M0, cfg)
    M, states, info = INS.fire(M0, cfg, free, max_iter=cfg["maxit"],
                               log_every=500, snaps=(), tag=tag)
    e_u, e_d, e_v = INS.e_parts(M, cfg)
    row = {k: cfg[k] for k in ("seed", "term", "stencil", "eps", "n",
                               "L", "h", "delta", "bc", "maxit")}
    row.update({
        "tag": tag, "E_end": float(e_u + e_d + e_v),
        "E_u": float(e_u), "E_d": float(e_d), "E_v": float(e_v),
        "E_seed": float(sum(e0)),
        "r_half": INS.r_half(M, cfg),
        "min_gap_end": INS.min_gap(M),
        "stop": info["stop"], "trace": info["trace"][-4:],
        "wall_s": info["wall_s"]})
    os.makedirs(DATA, exist_ok=True)
    np.savez_compressed(
        os.path.join(DATA, f"m5_21_5_end_{tag}.npz"),
        M=M.astype(np.float32), delta=cfg["delta"], h=cfg["h"],
        n=cfg["n"], bc=cfg["bc"], term=cfg["term"],
        maxit=cfg["maxit"])
    with open(os.path.join(DATA, f"m5_21_5_row_{tag}.json"),
              "w") as f:
        json.dump(row, f, indent=1)
    print(json.dumps({k: row[k] for k in
                      ("tag", "E_end", "E_u", "E_v", "r_half",
                       "stop", "wall_s")}), flush=True)
    return row


if __name__ == "__main__":
    kw = dict(a.split("=", 1) for a in sys.argv[1:])
    run(int(kw.get("n", 24)))
