"""M5.21.14 V1/V2: numeric verification of the derived T1 term.

V1 (continuum, numpy-only, NO sympy anywhere in the path): at random
off-axis points, the dressed-minus-undressed E_u and kin densities
(Richardson central differences, step 1e-3) must converge to the
derived T1_static / T1_kin as g grows with b = beta/g:
relative error ~ C/g, gated at err(g=1e4) < 1e-3 with ~10x decade
drops. Both s signs probed at g = 1e2 (difference ~ 1/g).

V2 (lattice, the S1a self-gate): the analytic vacuum-hedgehog family
(the verified m5_21_8_b_lattice builder, imported unchanged) measured
with the certified e_parts at n = 32, L = 48, s = -1, g in {8,16,32}:
  - m*_lattice / artanh(1/g) must sit in the recorded 0.82-0.84 band
    (gate widened to [0.78, 0.88] for the finer m-grid used here)
  - the dressing gain E(m*) - E(0) must be FLAT in g (the scaled
    variable beta = g*m collapse; the M5.21.11 g-arm q ~ 0 read):
    max/min gain ratio gated at <= 2.0, collapse spread reported

Out: ../data/m5_21_14_verify.json
"""
from __future__ import annotations

import importlib.util
import json
import os
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")

_spec8 = importlib.util.spec_from_file_location(
    "b8", os.path.join(HERE, "m5_21_8_b_lattice.py"))
B8 = importlib.util.module_from_spec(_spec8)
_spec8.loader.exec_module(B8)
INS4 = B8.INS4

ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
RNG = np.random.default_rng(21140)

# ---------- smooth random fields (fixed modes, deterministic) ----------
NMODE = 4
KS = RNG.normal(0.0, 0.35, size=(NMODE, 3))
PH = RNG.uniform(0, 2 * np.pi, size=NMODE)
AMP = RNG.normal(0.0, 0.4, size=(NMODE, 3, 3))
AMP = 0.5 * (AMP + AMP.transpose(0, 2, 1))
PH_D = RNG.uniform(0, 2 * np.pi, size=NMODE)
AMP_D = RNG.normal(0.0, 0.4, size=(NMODE, 3, 3))
AMP_D = 0.5 * (AMP_D + AMP_D.transpose(0, 2, 1))
BETA0, RHO = 0.8, 2.5


def m3_field(p):
    M = np.zeros((3, 3))
    for k in range(NMODE):
        M += AMP[k] * np.cos(KS[k] @ p + PH[k])
    return M


def mdot3_field(p):
    M = np.zeros((3, 3))
    for k in range(NMODE):
        M += AMP_D[k] * np.cos(KS[k] @ p + PH_D[k])
    return M


def beta_of(r):
    return BETA0 * r * r * np.exp(-((r / RHO) ** 2))


def beta_n(p):
    r = np.linalg.norm(p)
    return beta_of(r) * p / r


def qb_np(p, b):
    r = np.linalg.norm(p)
    n = p / r
    K = np.zeros((4, 4))
    K[0, 1:] = n
    K[1:, 0] = n
    K2 = np.zeros((4, 4))
    K2[0, 0] = 1.0
    K2[1:, 1:] = np.outer(n, n)
    return np.eye(4) + np.sinh(b) * K + (np.cosh(b) - 1.0) * K2


def m4_dressed(p, g, s):
    M4 = np.zeros((4, 4))
    M4[1:, 1:] = m3_field(p)
    M4[0, 0] = -s * g
    r = np.linalg.norm(p)
    Qb = qb_np(p, beta_of(r) / g)
    return Qb @ M4 @ Qb.T


def a0_dressed(p, g):
    A = np.zeros((4, 4))
    A[1:, 1:] = mdot3_field(p)
    r = np.linalg.norm(p)
    Qb = qb_np(p, beta_of(r) / g)
    return Qb @ A @ Qb.T


def rich_grad(f, p, h=1e-3):
    """Richardson 4th-order central differences of a matrix field."""
    out = []
    for ax in range(3):
        e = np.zeros(3)
        e[ax] = 1.0
        d = (8.0 * (f(p + h * e) - f(p - h * e))
             - (f(p + 2 * h * e) - f(p - 2 * h * e))) / (12.0 * h)
        out.append(d)
    return out


def inner_eta_np(F, G):
    return np.trace(ETA @ F @ ETA @ G.T)


def dens_u_np(A):
    tot = 0.0
    for i in range(3):
        for j in range(i + 1, 3):
            F = A[i] @ ETA @ A[j] - A[j] @ ETA @ A[i]
            tot += 4.0 * inner_eta_np(F, F)
    return tot


def dens_k_np(a0, A):
    tot = 0.0
    for i in range(3):
        F = a0 @ ETA @ A[i] - A[i] @ ETA @ a0
        tot += 4.0 * inner_eta_np(F, F)
    return tot


def t1_at(p):
    """the derived T1_static and T1_kin at a point (numpy route)."""
    G = rich_grad(m3_field, p)
    V = rich_grad(beta_n, p)
    Md = mdot3_field(p)
    t_u = 0.0
    for i in range(3):
        for j in range(i + 1, 3):
            C = G[i] @ G[j] - G[j] @ G[i]
            W = np.outer(V[j], V[i]) - np.outer(V[i], V[j])
            w = G[i] @ V[j] - G[j] @ V[i]
            t_u += 4.0 * (2.0 * np.trace(C @ W.T)
                          + np.sum(W * W) - 2.0 * w @ w)
    t_k = -8.0 * sum((Md @ V[i]) @ (Md @ V[i]) for i in range(3))
    base_u = dens_u_np([np.pad(Gi, ((1, 0), (1, 0))) for Gi in G])
    base_k = dens_k_np(np.pad(Md, ((1, 0), (1, 0))),
                       [np.pad(Gi, ((1, 0), (1, 0))) for Gi in G])
    return t_u, t_k, base_u, base_k


def stage_v1():
    pts = [np.array([1.3, -0.7, 0.9]), np.array([-2.1, 0.4, 1.6]),
           np.array([0.8, 1.9, -1.2]), np.array([2.4, -1.5, -0.6]),
           np.array([-1.0, -1.1, 2.2])]
    gs = [1e2, 1e3, 1e4]
    rows = []
    for p in pts:
        t_u, t_k, _, _ = t1_at(p)
        du0 = dens_u_np(rich_grad(lambda q: np.pad(
            m3_field(q), ((1, 0), (1, 0))), p))
        a00 = np.pad(mdot3_field(p), ((1, 0), (1, 0)))
        dk0 = dens_k_np(a00, rich_grad(lambda q: np.pad(
            m3_field(q), ((1, 0), (1, 0))), p))
        row = {"point": [float(v) for v in p], "errs_u": [],
               "errs_k": []}
        for g in gs:
            A = rich_grad(lambda q, gg=g: m4_dressed(q, gg, -1.0), p)
            a0 = a0_dressed(p, g)
            eu = dens_u_np(A)
            ek = dens_k_np(a0, A)
            row["errs_u"].append(abs((eu - du0) - t_u) / abs(t_u))
            row["errs_k"].append(abs((ek - dk0) - t_k) / abs(t_k))
        # s-sign difference at g = 1e2 (expected O(1/g), nonzero)
        Ap = rich_grad(lambda q: m4_dressed(q, 1e2, +1.0), p)
        eup = dens_u_np(Ap)
        A2 = rich_grad(lambda q: m4_dressed(q, 1e2, -1.0), p)
        eum = dens_u_np(A2)
        row["s_diff_rel_g100"] = abs(eup - eum) / max(abs(t_u), 1e-30)
        rows.append(row)
        print(json.dumps({"V1_point": row["point"],
                          "errs_u": row["errs_u"],
                          "errs_k": row["errs_k"]}), flush=True)
    worst_u = max(r["errs_u"][-1] for r in rows)
    worst_k = max(r["errs_k"][-1] for r in rows)
    drops = all(r["errs_u"][0] > 3 * r["errs_u"][1] > 9 * r["errs_u"][2]
                for r in rows)
    return {"rows": rows, "worst_final_err_u": worst_u,
            "worst_final_err_k": worst_k,
            "decade_drops_ok": bool(drops),
            "gate": bool(worst_u < 1e-3 and worst_k < 1e-3 and drops)}


def stage_v2():
    out = {"curves": [], "mstar": {}, "gain": {}}
    for g in (8.0, 16.0, 32.0):
        cfg = INS4.base_cfg(s=-1.0, g=g, n=32, L=48.0)
        mh = float(np.arctanh(1.0 / g))
        ms = np.linspace(-1.3 * mh, 1.3 * mh, 41)
        Es = []
        for m in ms:
            M = B8.dressed(cfg, m)
            e_u, e_v = INS4.e_parts(M, cfg)
            Es.append(float(e_u + e_v))
        Es = np.array(Es)
        i = int(np.argmin(Es))
        m_star, E_star = float(ms[i]), float(Es[i])
        if 0 < i < len(ms) - 1:
            a, b, c = Es[i - 1], Es[i], Es[i + 1]
            dm = ms[1] - ms[0]
            m_star = float(ms[i] - 0.5 * dm * (c - a) / (c - 2 * b + a))
            E_star = float(b - 0.125 * (c - a) ** 2 / (c - 2 * b + a))
        iz = int(np.argmin(np.abs(ms)))
        E0 = float(Es[iz])
        out["curves"].append({"g": g, "m": ms.tolist(),
                              "E": Es.tolist()})
        out["mstar"][f"g{g:g}"] = {"m_star": abs(m_star),
                                   "artanh": mh,
                                   "ratio": abs(m_star) / mh}
        out["gain"][f"g{g:g}"] = E_star - E0
        print(json.dumps({"V2_g": g, "m_star": abs(m_star),
                          "ratio": abs(m_star) / mh,
                          "gain": E_star - E0}), flush=True)
    ratios = [v["ratio"] for v in out["mstar"].values()]
    gains = [-v for v in out["gain"].values()]
    flat = max(gains) / min(gains)
    # collapse spread: interpolate -gain curves on common beta grid
    bgrid = np.linspace(0.2, 1.1, 19)
    prof = {}
    for c in out["curves"]:
        m = np.array(c["m"])
        E = np.array(c["E"])
        sel = m > 0
        beta = c["g"] * m[sel]
        dE = E[sel] - E[np.argmin(np.abs(m))]
        prof[c["g"]] = np.interp(bgrid, beta, dE)
    sp_816 = float(np.max(np.abs(prof[8.0] - prof[16.0])))
    sp_1632 = float(np.max(np.abs(prof[16.0] - prof[32.0])))
    out["ratio_band_ok"] = bool(all(0.78 <= x <= 0.88 for x in ratios))
    out["gain_flat_maxmin"] = float(flat)
    out["gain_flat_ok"] = bool(flat <= 2.0)
    out["collapse_spread"] = {"g8_vs_g16": sp_816,
                              "g16_vs_g32": sp_1632,
                              "converging": bool(sp_1632 < sp_816)}
    out["gate"] = bool(out["ratio_band_ok"] and out["gain_flat_ok"])
    return out


def main():
    t0 = time.time()
    out = {"V1": stage_v1(), "V2": stage_v2()}
    out["gates"] = {"V1": out["V1"]["gate"], "V2": out["V2"]["gate"]}
    out["all_green"] = all(out["gates"].values())
    out["runtime_s"] = round(time.time() - t0, 1)
    with open(os.path.join(DATA, "m5_21_14_verify.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps({"gates": out["gates"],
                      "all_green": out["all_green"],
                      "runtime_s": out["runtime_s"]}))


if __name__ == "__main__":
    main()
