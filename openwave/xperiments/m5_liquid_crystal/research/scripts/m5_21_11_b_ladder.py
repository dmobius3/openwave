"""M5.21.11 ladder runner: the frozen route-(b) rung program (COMPUTE ONLY).

Executes the FROZEN protocol of findings/m5_21_11_framework.md §§ 2-3
against the instrument of record m5_21_2b_a_instrument.py, IMPORTED
UNCHANGED (this script introduces no energy code: seeds, FIRE, e_parts,
consistency, virial all come from the instrument module; the frozen
§ 1.1 "no new energy code" clause is honored by construction).

Frozen specification executed here (framework § 2-3):
    production rungs  delta in {0.30, 0.25, 0.20, 0.15, 0.12, 0.09,
                      0.07, 0.05}, branches A/C/B, N = 48, L = 48,
                      T2/sym/eps=0, w2 = 0.002758100 (the 2b § 1b
                      calibration pin), pinned, FIRE maxit 12000
    continuation      rung k+1 seeds from the rung-k endpoint of the
                      SAME branch and SAME N (first rung: the analytic
                      instrument seed); identity is NEVER assumed from
                      continuation (the § 5 readers decide, separately)
    refinement        delta in {0.30, 0.12, 0.05} x N in {32, 64}
                      (N = 48 = the production rungs), same scheme
                      within each N (pre-run interpretation pin,
                      recorded in the task doc: the framework's
                      continuation row applied uniformly per N)

Modes:
    chain branch=A n=48 deltas=0.30,0.25,...   sequential continuation
        chain, one background process per (branch, n); idempotent:
        completed rungs (row json present) are skipped on resume
    status                                     table of finished rungs

Out (all local, dataset policy): ../data/m5_21_11_row_<tag>.json,
../data/m5_21_11_end_<tag>.npz, tag = t11lad_<br>_n<n>_d<delta>.
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")

_spec = importlib.util.spec_from_file_location(
    "ins2b", os.path.join(HERE, "m5_21_2b_a_instrument.py"))
INS = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(INS)

W2_PIN = 0.002758100          # the 2b § 1b calibration pin
MAXIT = 12000                 # the census production depth (frozen § 2)


def tag_of(branch, n, delta):
    return f"t11lad_{branch}_n{n}_d{delta:g}"


def run_one(branch, n, delta, M_seed=None):
    """one frozen rung: instrument-of-record relaxation + endpoint row.

    Mirrors INS.relax exactly (same calls, same row fields) with one
    difference: M0 may be a continuation endpoint instead of the
    analytic seed. Energies, gates, FIRE all come from INS unchanged.
    """
    cfg = INS.base_cfg(seed=branch, term="T2", stencil="sym", eps=0.0,
                       n=n, L=48.0, delta=delta, bc="pinned",
                       maxit=MAXIT, w2=W2_PIN)
    tag = tag_of(branch, n, delta)
    if M_seed is None:
        M0 = INS.make_seed(cfg)
        seed_kind = "analytic"
    else:
        # continuation carries the INTERIOR only: the pinned shell is
        # part of the instrument (INS.relax pins the CURRENT-delta
        # analytic far field), so the shell cells are re-seeded at this
        # rung's own vacuum. Inheriting the old-delta shell adds a
        # spurious (delta_prev - delta)^2 penalty over the shell volume
        # (caught 2026-08-07 on the first contaminated rung: virial
        # residual 1.02 vs the census +0.03 grade).
        M0 = np.ascontiguousarray(M_seed.astype(np.float64))
        shell = INS.pin_shell(n, cfg["h"])
        M0[shell] = INS.make_seed(cfg)[shell]
        seed_kind = "continuation"
    free = ~INS.pin_shell(n, cfg["h"])
    e0 = INS.e_parts(M0, cfg)
    t0 = time.time()
    M, states, info = INS.fire(M0, cfg, free, max_iter=cfg["maxit"],
                               log_every=500, snaps=(), tag=tag)
    e_u, e_d, e_v = INS.e_parts(M, cfg)
    row = {k: cfg[k] for k in ("seed", "term", "stencil", "eps", "n",
                               "L", "h", "delta", "bc", "maxit", "w2")}
    row.update({
        "tag": tag, "branch": branch, "seed_kind": seed_kind,
        "E_end": float(e_u + e_d + e_v),
        "E_u": float(e_u), "E_d": float(e_d), "E_v": float(e_v),
        "E_seed": float(sum(e0)),
        "virial_resid": float((-e_u + e_d + 3 * e_v)
                              / max(e_u + e_d + e_v, 1e-300)),
        "u_over_3v": float(e_u / max(3 * e_v, 1e-300)),
        "r_half": INS.r_half(M, cfg),
        "ring": INS.ring_read(M, cfg),
        "consistency": INS.consistency(M, cfg),
        "min_gap_end": INS.min_gap(M),
        "stop": info["stop"], "trace": info["trace"][-6:],
        "wall_s": time.time() - t0})
    os.makedirs(DATA, exist_ok=True)
    np.savez_compressed(
        os.path.join(DATA, f"m5_21_11_end_{tag}.npz"),
        M=M, n=n, delta=delta, h=cfg["h"], branch=branch,
        w2=W2_PIN, maxit=MAXIT, seed_kind=seed_kind)
    with open(os.path.join(DATA, f"m5_21_11_row_{tag}.json"),
              "w") as f:
        json.dump(row, f, indent=1)
    print(json.dumps({k: row[k] for k in
                      ("tag", "E_end", "E_u", "E_v", "virial_resid",
                       "min_gap_end", "stop", "wall_s")}
                     | {"xratio":
                        row["consistency"]["xstencil_ratio"]}),
          flush=True)
    return M, row


def chain(branch, n, deltas):
    M_prev = None
    for i, d in enumerate(deltas):
        tag = tag_of(branch, n, d)
        rowp = os.path.join(DATA, f"m5_21_11_row_{tag}.json")
        endp = os.path.join(DATA, f"m5_21_11_end_{tag}.npz")
        if os.path.exists(rowp) and os.path.exists(endp):
            print(f"[chain {branch} n{n}] rung d={d:g} already done, "
                  f"loading endpoint", flush=True)
            M_prev = np.load(endp)["M"].astype(np.float64)
            continue
        print(f"[chain {branch} n{n}] rung {i + 1}/{len(deltas)} "
              f"d={d:g} seed="
              f"{'analytic' if M_prev is None else 'continuation'}",
              flush=True)
        M_prev, _ = run_one(branch, n, d, M_seed=M_prev)
    print(f"[chain {branch} n{n}] COMPLETE", flush=True)


def status():
    import glob
    for p in sorted(glob.glob(os.path.join(DATA,
                                           "m5_21_11_row_t11lad_*"))):
        r = json.load(open(p))
        print(f"{r['tag']:28s} E {r['E_end']:10.4f} "
              f"vir {r['virial_resid']:+.3f} "
              f"xr {r['consistency']['xstencil_ratio']:5.2f} "
              f"gap {r['min_gap_end']:.4f} stop {r['stop']} "
              f"[{r['wall_s']:.0f}s]")


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "status"
    kw = dict(a.split("=", 1) for a in sys.argv[2:])
    if mode == "chain":
        chain(kw["branch"], int(kw["n"]),
              [float(x) for x in kw["deltas"].split(",")])
    else:
        status()
