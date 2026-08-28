"""M5.32 R3 arm (ii) INDEPENDENT ADVERSARIAL AUDIT of the relaxed
boost-dressed pair verdict (light, data-only: no relaxation is run).

Inputs (read only): ../data/m5_32_r3_pair.json (the 44 relaxation rows)
and the end states ../data/m5_32_r3_ii/*.npz. The producer's script is
NOT read; its tables are recomputed from the row energies and its end
states are re-evaluated with the certified oracles (m5_32_lagrangian.py,
m5_32_r2_b_bounded.py energy_grad / reads).

EQUATIONS FIRST
---------------
Candidate (arm (ii) of R3):   L_lambda = -4 [(1 - lambda) I1 + lambda I1_h] - V4
Static lattice energy (the R2.b instrument, omega = 0):
    E_lambda(M) = 4 h^3 sum_br wt sum_cells sum_{i<j} q_lambda(F_ij) + V4(M)
    F_ij = A_i eta A_j - A_j eta A_i,  A_i = d_i M (sym stencil branches)
    q_lambda(F) = (1 - lambda) <F, F>_eta + lambda <F, F>_h
    <F, F>_eta = eta_ac eta_bd F^ab F^cd,   <F, F>_h = h_ac h_bd F^ab F^cd
    h_ab = eta_ab + 2 (eta u)_a (eta u)_b,  u = timelike unit eigenvector
           of M eta (u^T eta u = -1)

(A1) Interaction energies from the row energies E(row):
    E_int(kind, d; lambda, dressed) = E(pair) - 2 E(single)          [same lambda, dressing, box]
    DP(kind, d; lambda) = E_int(dressed) - E_int(undressed)          [the dressing part]
    undressed fields have M_0i = 0 so u = e0, h = 1 and q_lambda = q_eta
    for every lambda: the undressed E_int is lambda-independent and the
    lambda = 0 undressed rows serve every lambda (checked on the lambda = 1
    undressed rows that exist, 0 difference expected).
    outer slope   s = [X(24) - X(18)] / 6
    fit_inv       X(d) ~ A + B/d                    (least squares, 4 points)
    fit_log       X(d) ~ A + B/d + C ln d / d
    log-log       ln|X| ~ a + b ln d                (only if X has one sign)
    Sign convention: like charges REPEL when E_int FALLS with d (s < 0, B > 0).

(A2) The S / T / V4 split (the auditor's own decomposition, u-frame;
    F_ij is antisymmetric, so F'_00 = 0 and T_h = -T_eta identically:
    E_lambda = E_S + (1 - 2 lambda) E_Teta + V4, checked numerically):
    P^a_b = delta^a_b + u^a (eta u)_b       (P u = 0; projector onto the
                                             eta-orthogonal complement of u)
    F_S = P F P^T                            (the purely spatial part in the
                                             frame where u = e0)
    S(F)     = <F_S, F_S>_eta = <F_S, F_S>_h      (both metrics agree on F_S)
    T_eta(F) = <F, F>_eta - S(F)                  (time-row part, eta metric)
    T_h(F)   = <F, F>_h   - S(F)                  (time-row part, flip metric)
    E_lambda = E_S + (1 - lambda) E_Teta + lambda E_Th + V4,
    each E_X = 4 h^3 sum_br wt sum_cells sum_{i<j} X(F_ij).
    Check: E_lambda(re-evaluated) == stored row E to <= 1e-8 relative and
    == the oracle energy_grad(M, cfg, lambda)[0].
    In the u-frame T_h - T_eta = 2 sum_i (F'_0i^2 + F'_i0^2) >= 0: the flip
    metric raises exactly the mixed time-space rows.

(A3) Uptick vs residual heal drift: every row stops at the accepted-step
    budget still falling. The stored trace E(acc) is fit linearly over the
    whole run (R2_lin) and over its first and second halves (slopes s1, s2;
    deceleration ratio s2 / s1): a ratio near 1 means a constant-rate
    slide with no minimum in sight, and then no extrapolation to E_inf is
    possible from the data. Per row the one-budget drift D = -s2 x 1500.
    For an 18 -> 24 uptick U of the dressing part:
        drift-stability: |U| vs the RATE DIFFERENCE of the participating
            rows over one more budget, dD = [D(p24) - D(u24)] - [D(p18) - D(u18)]
        convergence: |U| vs D of each pair row (the unresolved remainder)
    U is called SIGNIFICANT only if |U| > 2 max(D(p18), D(p24)) AND the
    pair rows decelerate (s2 / s1 < 0.5); DRIFT-STABLE if |U| > 3 |dD|.
    The seed-level uptick is recomputed from the stored seed energies.

(A4) Dressing amplitude: max |M_0i| over the grid, in a ball of radius
    r_ball = 5 around each core (cores = the two largest peaks of the
    curvature density q_eta(F) separated by >= 6 cells) and in the far
    field (distance > 8 from both cores), at the end state (npz) vs the
    seed (the row's seed_amp, and the oracle-reconstructed single seed
    dressed_electron(cfg, scale)).

(A5) The static like-charge control: E_int(same, d; undressed) from the
    rows; a RISE with d means the static baseline is ATTRACTIVE on this
    instrument.

Runtime: ~35 s, one process (n = 32 and n = 48 end states re-evaluated).
Out: ../data/m5_32_r3_audit_pair.json
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
NPZ = os.path.join(DATA, "m5_32_r3_ii")
OUT = os.path.join(DATA, "m5_32_r3_audit_pair.json")
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


R2B = _load("m5_32_r2_b_bounded", "m5_32_r2_b_bounded.py")
LAG = R2B.LAG
B3 = R2B.B3
ETA = B3.ETA

LAMS = (0.0, 0.75, 1.0)
DS = (10.0, 14.0, 18.0, 24.0)


# ================= A1: row-level recompute =================
def row_key(r):
    return (float(r["lam"]), r["kind"], float(r["d"]), bool(r["dressed"]),
            int(r["n"]))


def fit_inv(ds, xs):
    X = np.column_stack([np.ones(len(ds)), 1.0 / np.asarray(ds)])
    c, *_ = np.linalg.lstsq(X, np.asarray(xs), rcond=None)
    pred = X @ c
    ss = np.sum((np.asarray(xs) - np.mean(xs)) ** 2)
    r2 = 1.0 - np.sum((np.asarray(xs) - pred) ** 2) / ss if ss > 0 else 1.0
    return {"A": float(c[0]), "B": float(c[1]), "R2": float(r2)}


def fit_log(ds, xs):
    d = np.asarray(ds)
    X = np.column_stack([np.ones(len(d)), 1.0 / d, np.log(d) / d])
    c, *_ = np.linalg.lstsq(X, np.asarray(xs), rcond=None)
    pred = X @ c
    ss = np.sum((np.asarray(xs) - np.mean(xs)) ** 2)
    r2 = 1.0 - np.sum((np.asarray(xs) - pred) ** 2) / ss if ss > 0 else 1.0
    return {"A": float(c[0]), "B": float(c[1]), "C": float(c[2]),
            "R2": float(r2)}


def fit_loglog(ds, xs):
    x = np.asarray(xs)
    if not (np.all(x > 0) or np.all(x < 0)):
        return {"note": "mixed sign, no log-log fit"}
    X = np.column_stack([np.ones(len(ds)), np.log(ds)])
    c, *_ = np.linalg.lstsq(X, np.log(np.abs(x)), rcond=None)
    return {"a": float(c[0]), "b_exponent": float(c[1])}


def sequence_report(ds, xs):
    xs = [float(v) for v in xs]
    diffs = np.diff(xs)
    mono = ("decreasing" if np.all(diffs < 0) else
            "increasing" if np.all(diffs > 0) else "non-monotone")
    slope = (xs[-1] - xs[-2]) / (ds[-1] - ds[-2])
    return {"points": [[float(d), x] for d, x in zip(ds, xs)],
            "outer_slope_18_24": float(slope),
            "outer_sign": "REPULSIVE" if slope < 0 else "ATTRACTIVE",
            "monotone": mono,
            "fit_inv": fit_inv(ds, xs), "fit_log": fit_log(ds, xs),
            "fit_loglog": fit_loglog(ds, xs)}


def stage_a1(rows):
    E = {row_key(r): float(r["E"]) for r in rows}
    out = {"n_rows": len(rows), "lambda_independence_undressed": {},
           "Eint": {}, "sequences": {}, "n48_ladder": {}}
    # undressed lambda-independence (rows that exist at lambda = 1)
    for k, v in E.items():
        lam, kind, d, dr, n = k
        if not dr and lam != 0.0:
            k0 = (0.0, kind, d, False, n)
            if k0 in E:
                out["lambda_independence_undressed"][
                    f"{kind}_d{d:g}_n{n}"] = {
                    "E_lam0": E[k0], f"E_lam{lam:g}": v,
                    "abs_diff": abs(E[k0] - v)}

    def eint(lam, kind, d, dr, n):
        ks = (lam, kind, d, dr, n)
        if ks not in E:
            return None
        # undressed singles are lambda-independent: fall back to lambda 0
        for sl in (lam, 0.0):
            if (sl, "single", 0.0, dr, n) in E:
                return E[ks] - 2.0 * E[(sl, "single", 0.0, dr, n)]
        return None

    for lam in LAMS:
        for kind in ("same", "anti"):
            seq = {"static_undressed": [], "dressed_total": [],
                   "dressing_part": []}
            for d in DS:
                st = eint(0.0, kind, d, False, 32)   # lambda-independent
                dt = eint(lam, kind, d, True, 32)
                out["Eint"][f"lam{lam:g}_{kind}_d{d:g}_n32"] = {
                    "Eint_static_undressed": st, "Eint_dressed_total": dt,
                    "dressing_part": None if (st is None or dt is None)
                    else dt - st}
                seq["static_undressed"].append(st)
                seq["dressed_total"].append(dt)
                seq["dressing_part"].append(
                    None if (st is None or dt is None) else dt - st)
            for name, xs in seq.items():
                if all(x is not None for x in xs):
                    out["sequences"][f"lam{lam:g}_{kind}_{name}"] = \
                        sequence_report(DS, xs)
    # the n = 48, L = 72 ladder (lambda 1, same, d 14 / 24)
    lad = {}
    for d in (14.0, 24.0):
        st = eint(1.0, "same", d, False, 48)
        dt = eint(1.0, "same", d, True, 48)
        lad[f"d{d:g}"] = {"Eint_static_undressed": st,
                          "Eint_dressed_total": dt,
                          "dressing_part": dt - st}
    dp14, dp24 = lad["d14"]["dressing_part"], lad["d24"]["dressing_part"]
    lad["slope_14_24_n48"] = (dp24 - dp14) / 10.0
    e32 = out["Eint"]
    dp32_14 = e32["lam1_same_d14_n32"]["dressing_part"]
    dp32_24 = e32["lam1_same_d24_n32"]["dressing_part"]
    lad["slope_14_24_n32"] = (dp32_24 - dp32_14) / 10.0
    lad["slope_ratio_n48_over_n32"] = lad["slope_14_24_n48"] / \
        lad["slope_14_24_n32"]
    lad["dressing_part_n48_over_n32"] = {
        "d14": dp14 / dp32_14, "d24": dp24 / dp32_24}
    out["n48_ladder"] = lad
    return out


# ================= A2: re-evaluation + the S/T/V4 split =================
def st_split(M, cfg):
    """(E_S, E_Teta, E_Th, V4) with the u-frame projector."""
    h3 = cfg["h"] ** 3
    u0, _, _, _, _, ok, gap = R2B.tl_eig(M)
    assert np.all(ok), "timelike eigenvector missing somewhere"
    hh = R2B.h_of(u0)
    eu = u0 @ ETA                                        # (eta u)_b
    P = np.eye(4) + u0[..., :, None] * eu[..., None, :]  # P^a_b
    assert np.max(np.abs(np.einsum("...ab,...b->...a", P, u0))) < 1e-9
    accS = accE = accH = 0.0
    chk = 0.0
    for br, wt in B3.branches(cfg["stencil"]):
        A = [B3.d1(M, ax, cfg["h"], br) for ax in range(3)]
        for i in range(3):
            for j in range(i + 1, 3):
                F = B3.comm_eta(A[i], A[j])
                FS = np.einsum("...ac,...cd,...bd->...ab", P, F, P)
                s_eta = R2B.q_eta(FS)
                s_h = R2B.q_h(FS, hh)
                chk = max(chk, float(np.max(np.abs(s_eta - s_h))))
                accS += wt * np.sum(s_eta)
                accE += wt * np.sum(R2B.q_eta(F) - s_eta)
                accH += wt * np.sum(R2B.q_h(F, hh) - s_eta)
    _, ev = B3.e_parts(M, cfg)
    return {"E_S": float(4 * h3 * accS), "E_Teta": float(4 * h3 * accE),
            "E_Th": float(4 * h3 * accH), "V4": float(ev),
            "max_abs_S_eta_minus_S_h": chk, "min_gap": float(np.min(gap))}


def e_of_split(sp, lam):
    return sp["E_S"] + (1 - lam) * sp["E_Teta"] + lam * sp["E_Th"] + sp["V4"]


def stage_a2(rows):
    E = {r["tag"]: float(r["E"]) for r in rows}
    tags = [f"lam{lam:g}_dr1_same_d{d:g}_n32" for lam in LAMS for d in DS] + \
           [f"lam{lam:g}_dr1_single_d0_n32" for lam in LAMS] + \
           [f"lam0_un0_same_d{d:g}_n32" for d in DS] + \
           ["lam0_un0_single_d0_n32",
            "lam1_dr1_same_d14_n48", "lam1_dr1_same_d24_n48",
            "lam1_dr1_single_d0_n48", "lam1_un0_same_d14_n48",
            "lam1_un0_same_d24_n48", "lam1_un0_single_d0_n48"]
    splits, checks = {}, {}
    for tag in tags:
        n = 48 if tag.endswith("n48") else 32
        L = 72.0 if n == 48 else 48.0
        cfg = R2B.cfg_of(n, L)
        M = np.load(os.path.join(NPZ, tag + ".npz"))["M"]
        lam = float(tag.split("_")[0][3:])
        sp = st_split(M, cfg)
        splits[tag] = sp
        e_or, _, _ = R2B.energy_grad(M, cfg, lam)
        e_me = e_of_split(sp, lam)
        checks[tag] = {"lam": lam, "E_stored": E[tag], "E_oracle": e_or,
                       "E_split": e_me,
                       "rel_err_oracle_vs_stored": abs(e_or - E[tag]) / E[tag],
                       "rel_err_split_vs_stored": abs(e_me - E[tag]) / E[tag],
                       "E_split_lam0": e_of_split(sp, 0.0),
                       "E_split_lam1": e_of_split(sp, 1.0)}
        log(f"A2 {tag}: stored {E[tag]:.6f} oracle {e_or:.6f} "
            f"split {e_me:.6f} rel {checks[tag]['rel_err_split_vs_stored']:.1e}")

    def dp_split(lam, d, n):
        tg = f"lam{lam:g}_dr1_same_d{d:g}_n{n}"
        sg = f"lam{lam:g}_dr1_single_d0_n{n}"
        ulam = 0 if n == 32 else 1
        ug = f"lam{ulam}_un0_same_d{d:g}_n{n}"
        us = f"lam{ulam}_un0_single_d0_n{n}"
        out = {}
        for comp in ("E_S", "E_Teta", "E_Th", "V4"):
            out[comp] = (splits[tg][comp] - 2 * splits[sg][comp]) - \
                (splits[ug][comp] - 2 * splits[us][comp])
        out["T_lambda"] = (1 - lam) * out["E_Teta"] + lam * out["E_Th"]
        out["dressing_part_total"] = out["E_S"] + out["T_lambda"] + out["V4"]
        return out

    dp = {}
    for lam in LAMS:
        for d in DS:
            dp[f"lam{lam:g}_same_d{d:g}_n32"] = dp_split(lam, d, 32)
        k18, k24 = f"lam{lam:g}_same_d18_n32", f"lam{lam:g}_same_d24_n32"
        dp[f"lam{lam:g}_same_delta_24_minus_18_n32"] = {
            c: dp[k24][c] - dp[k18][c] for c in dp[k18]}
    for d in (14.0, 24.0):
        dp[f"lam1_same_d{d:g}_n48"] = dp_split(1.0, d, 48)
    dp["lam1_same_delta_24_minus_14_n48"] = {
        c: dp["lam1_same_d24_n48"][c] - dp["lam1_same_d14_n48"][c]
        for c in dp["lam1_same_d14_n48"]}
    return {"splits": splits, "checks": checks, "dressing_part_split": dp}


# ================= A3: uptick vs residual drift =================
def trace_drift(trace, budget=1500):
    acc = np.array([t["acc"] for t in trace], float)
    e = np.array([t["E"] for t in trace], float)
    X = np.column_stack([np.ones(len(acc)), acc])
    c, *_ = np.linalg.lstsq(X, e, rcond=None)
    pred = X @ c
    ss = np.sum((e - e.mean()) ** 2)
    r2 = 1.0 - np.sum((e - pred) ** 2) / ss if ss > 0 else 1.0
    mid = acc[-1] / 2
    s1 = np.polyfit(acc[acc <= mid], e[acc <= mid], 1)[0]
    s2 = np.polyfit(acc[acc > mid], e[acc > mid], 1)[0]
    return {"E_end": float(e[-1]), "slope_full": float(c[1]),
            "R2_linear": float(r2), "slope_first_half": float(s1),
            "slope_second_half": float(s2),
            "decel_ratio_s2_over_s1": float(s2 / s1) if s1 != 0 else None,
            "one_budget_drift": float(-s2 * budget)}


def stage_a3(rows, rec, a1):
    R = {r["tag"]: r for r in rows}
    ext = {}
    for tag, r in R.items():
        ext[tag] = trace_drift(r["descent"]["trace"])
        ext[tag]["E_drop_total"] = float(r["descent"]["E_drop"])
        ext[tag]["last_quarter_dE"] = float(r["descent"]["last_quarter_dE"])
    # seed-level dressing parts from the stored seed energies
    SE = {k: float(v["E_total"]) for k, v in rec["seed"]["seed_energies"].items()}
    seed_dp = {}
    for lam in LAMS:
        for kind in ("same", "anti"):
            pts = []
            for d in DS:
                ei_d = SE[f"lam{lam:g}_dr1_{kind}_d{d:g}_n32"] - \
                    2 * SE[f"lam{lam:g}_dr1_single_d0_n32"]
                ei_u = SE[f"lam{lam:g}_un0_{kind}_d{d:g}_n32"] - \
                    2 * SE[f"lam{lam:g}_un0_single_d0_n32"]
                pts.append([d, ei_d - ei_u])
            seed_dp[f"lam{lam:g}_{kind}"] = {
                "points": pts, "uptick_18_24": pts[3][1] - pts[2][1]}
    out = {"per_row": ext, "seed_dressing_part": seed_dp, "upticks": {}}
    for lam in (0.75, 1.0):
        for kind in ("same", "anti"):
            seq = a1["sequences"][f"lam{lam:g}_{kind}_dressing_part"]["points"]
            dp18 = [p[1] for p in seq if p[0] == 18.0][0]
            dp24 = [p[1] for p in seq if p[0] == 24.0][0]
            U = dp24 - dp18
            tags = {"p18": f"lam{lam:g}_dr1_{kind}_d18_n32",
                    "p24": f"lam{lam:g}_dr1_{kind}_d24_n32",
                    "u18": f"lam0_un0_{kind}_d18_n32",
                    "u24": f"lam0_un0_{kind}_d24_n32"}
            D = {k: ext[t]["one_budget_drift"] for k, t in tags.items()}
            dD = (D["p24"] - D["u24"]) - (D["p18"] - D["u18"])
            dec = max(ext[tags["p18"]]["decel_ratio_s2_over_s1"],
                      ext[tags["p24"]]["decel_ratio_s2_over_s1"])
            lq = [ext[tags["p18"]]["last_quarter_dE"],
                  ext[tags["p24"]]["last_quarter_dE"]]
            out["upticks"][f"lam{lam:g}_{kind}"] = {
                "dressing_part_d18": dp18, "dressing_part_d24": dp24,
                "uptick_end": U,
                "uptick_seed": seed_dp[f"lam{lam:g}_{kind}"]["uptick_18_24"],
                "one_budget_drift_rows": D,
                "rate_difference_one_budget": dD,
                "pair_rows_R2_linear": [ext[tags["p18"]]["R2_linear"],
                                        ext[tags["p24"]]["R2_linear"]],
                "pair_rows_decel_ratio_max": dec,
                "pair_rows_last_quarter_dE": lq,
                "pair_rows_E_drop_total": [ext[tags["p18"]]["E_drop_total"],
                                           ext[tags["p24"]]["E_drop_total"]],
                "drift_stable": bool(abs(U) > 3 * abs(dD)),
                "significant": bool(abs(U) > 2 * max(D["p18"], D["p24"])
                                    and dec < 0.5)}
    return out


# ================= A4: amplitudes =================
def find_cores(M, cfg, min_sep_cells=6):
    """two largest peaks of the eta curvature density."""
    dens = 0.0
    for br, wt in B3.branches(cfg["stencil"]):
        A = [B3.d1(M, ax, cfg["h"], br) for ax in range(3)]
        for i in range(3):
            for j in range(i + 1, 3):
                F = B3.comm_eta(A[i], A[j])
                dens = dens + wt * np.abs(R2B.q_eta(F))
    n = cfg["n"]
    idx = np.indices((n, n, n))
    c1 = np.unravel_index(np.argmax(dens), dens.shape)
    dist1 = np.sqrt(sum((idx[k] - c1[k]) ** 2 for k in range(3)))
    masked = np.where(dist1 >= min_sep_cells, dens, -np.inf)
    c2 = np.unravel_index(np.argmax(masked), dens.shape)
    return c1, c2, float(dens[c1]), float(dens[c2])


def amp_report(M, cfg, cores, r_ball=5.0, r_far=8.0):
    n, h = cfg["n"], cfg["h"]
    idx = np.indices((n, n, n))
    m0i = np.max(np.abs(M[..., 0, 1:]), axis=-1)
    out = {"grid_max_abs_M0i": float(np.max(m0i))}
    dmin = np.inf
    for k, c in enumerate(cores):
        dist = h * np.sqrt(sum((idx[a] - c[a]) ** 2 for a in range(3)))
        ball = dist <= r_ball
        out[f"core{k}_ijk"] = [int(v) for v in c]
        out[f"core{k}_ball_max_abs_M0i"] = float(np.max(m0i[ball]))
        out[f"core{k}_ball_mean_abs_M0i"] = float(np.mean(m0i[ball]))
        dmin = np.minimum(dmin, dist)
    far = dmin > r_far
    out["far_max_abs_M0i"] = float(np.max(m0i[far]))
    out["far_mean_abs_M0i"] = float(np.mean(m0i[far]))
    return out


def stage_a4(rows):
    R = {r["tag"]: r for r in rows}
    out = {}
    for tag in ["lam1_dr1_same_d18_n32", "lam1_dr1_same_d24_n32",
                "lam1_dr1_anti_d24_n32", "lam0.75_dr1_same_d24_n32",
                "lam0_dr1_same_d18_n32", "lam1_dr1_single_d0_n32",
                "lam0_dr1_single_d0_n32"]:
        r = R[tag]
        cfg = R2B.cfg_of(int(r["n"]), float(r["L"]))
        M = np.load(os.path.join(NPZ, tag + ".npz"))["M"]
        c1, c2, p1, p2 = find_cores(M, cfg)
        cores = (c1, c2) if r["kind"] != "single" else (c1,)
        rep = amp_report(M, cfg, cores)
        if len(cores) == 2:
            rep["core_separation_units"] = float(cfg["h"] * np.sqrt(
                sum((c1[k] - c2[k]) ** 2 for k in range(3))))
            rep["d_nominal"] = float(r["d"])
        rep["seed_grid_max_abs_M0i"] = r["seed_amp"]["grid_max_abs_M0i"]
        rep["producer_end_grid_max_abs_M0i"] = r["end_amp"]["grid_max_abs_M0i"]
        rep["rel_change_grid_max_seed_to_end"] = (
            rep["grid_max_abs_M0i"] - rep["seed_grid_max_abs_M0i"]) / \
            rep["seed_grid_max_abs_M0i"]
        rep["seed_top_ball_max_norm_M0i"] = r["seed_amp"].get(
            "top_ball_max_norm_M0i")
        if r["kind"] == "single":
            Ms, _ = R2B.dressed_electron(cfg, scale=float(r["scale"]))
            reps = amp_report(Ms, cfg, cores)
            rep["oracle_seed_grid_max_abs_M0i"] = reps["grid_max_abs_M0i"]
            rep["oracle_seed_core0_ball_max_abs_M0i"] = \
                reps["core0_ball_max_abs_M0i"]
            rep["oracle_seed_far_max_abs_M0i"] = reps["far_max_abs_M0i"]
            rep["end_minus_oracle_seed_max_abs_field"] = float(
                np.max(np.abs(M - Ms)))
        out[tag] = rep
        log(f"A4 {tag}: grid max {rep['grid_max_abs_M0i']:.4f} "
            f"(seed {rep['seed_grid_max_abs_M0i']:.4f}) far "
            f"{rep['far_max_abs_M0i']:.4f}")
    return out


# ================= A5: static control =================
def stage_a5(a1):
    seq = a1["sequences"]["lam0_same_static_undressed"]
    seqa = a1["sequences"]["lam0_anti_static_undressed"]
    return {"same": seq, "anti": seqa,
            "same_rises_with_d": seq["monotone"] == "increasing",
            "same_10_to_24": [seq["points"][0][1], seq["points"][-1][1]]}


# ================= main =================
def main():
    with open(os.path.join(DATA, "m5_32_r3_pair.json")) as f:
        rec = json.load(f)
    rows = rec["rows"]
    log(f"loaded {len(rows)} rows")
    a1 = stage_a1(rows)
    for k, v in a1["sequences"].items():
        log(f"A1 {k}: {[round(p[1], 2) for p in v['points']]} slope "
            f"{v['outer_slope_18_24']:+.2f} {v['outer_sign']} "
            f"B_inv {v['fit_inv']['B']:+.1f} B_log {v['fit_log']['B']:+.1f}")
    log(f"A1 n48 ladder {json.dumps(a1['n48_ladder'])}")
    a5 = stage_a5(a1)
    a3 = stage_a3(rows, rec, a1)
    for k, v in a3["upticks"].items():
        log(f"A3 {k}: uptick end {v['uptick_end']:+.1f} seed "
            f"{v['uptick_seed']:+.1f} drift/budget "
            f"{[round(x, 1) for x in v['one_budget_drift_rows'].values()]} "
            f"rate-diff {v['rate_difference_one_budget']:+.1f} decel "
            f"{v['pair_rows_decel_ratio_max']:.2f} R2lin "
            f"{[round(x, 4) for x in v['pair_rows_R2_linear']]} "
            f"drift_stable {v['drift_stable']} significant {v['significant']}")
    a2 = stage_a2(rows)
    for k, v in a2["dressing_part_split"].items():
        log(f"A2 DP split {k}: S {v['E_S']:+.1f} Teta {v['E_Teta']:+.1f} "
            f"Th {v['E_Th']:+.1f} V4 {v['V4']:+.3f} total "
            f"{v['dressing_part_total']:+.1f}")
    a4 = stage_a4(rows)
    out = {"task": "M5.32 R3 arm (ii) independent audit (data-only)",
           "inputs": ["m5_32_r3_pair.json", "m5_32_r3_ii/*.npz"],
           "oracles": ["m5_32_r2_b_bounded.energy_grad / tl_eig / h_of / "
                       "dressed_electron", "m5_21_3_a_4d (B3) stencil + V4"],
           "A1": a1, "A2": a2, "A3": a3, "A4": a4, "A5": a5,
           "wall_s": round(time.time() - T0, 1)}
    with open(OUT, "w") as f:
        json.dump(out, f, indent=1, default=float)
    log(f"wrote {OUT}")


if __name__ == "__main__":
    main()
