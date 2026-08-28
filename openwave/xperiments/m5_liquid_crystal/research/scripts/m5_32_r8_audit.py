"""M5.32 R8 ADVERSARIAL AUDIT: classes C5 and C6 under the R7 IR obstruction.

An INDEPENDENT rebuild of claims C1..C10. The producer's scripts
(m5_32_r8_a_quartics.py, m5_32_r8_b_ir_theorem.py, m5_32_r8_c_collect.py)
were NOT read and are NOT imported. Every quartic density below is
implemented from the written definition in the brief. Oracles: the
certified stack (m5_21_3_a_4d.py: coords / base_cfg / vac4 / e_parts /
kin_of) and the m5_21_8_b_lattice.py hedgehog + clock builders, plus my
own predecessor auditor m5_32_r7_audit.py (differencing and shell
harness re-derived here, not imported).

EQUATIONS FIRST
---------------
eta = diag(-1, 1, 1, 1), index 0 = time. M real symmetric 4x4, RAW
CONTRAVARIANT internal entries, M -> L M L^T. Jets A_mu = d_mu M with
A_0 = omega a0, a0 = B8.a0_unit(cfg, 0) (the REALIZED clock channel).

    F_munu       = A_mu eta A_nu - A_nu eta A_mu
    <F,G>_eta    = tr(eta F eta G^T) = sum_ab w_a w_b F_ab G_ab,  w = diag(eta)
    T_munu       = tr(A_mu eta A_nu eta) = sum_ab w_a w_b (A_mu)_ab (A_nu)_ba
    I1           = sum_{mu<nu} w_mu w_nu <F_munu, F_munu>_eta
    R[nu,a]      = sum_mu F[mu,nu][a,mu]
    I4           = sum_{nu,a} w_nu w_a R[nu,a]^2

    Q_I1sq       = I1^2
    Q_I4sq       = I4^2
    Q_Fpair      = sum_{mu<nu} sum_{rho<sig} w_mu w_nu w_rho w_sig
                       <F_munu, F_rhosig>_eta^2
    Q_C6a        = [sum_mu w_mu T_mumu]^2
    Q_C6b        = sum_{mu,nu} w_mu w_nu T_munu^2
    Q_BI         = b^2 (sqrt(1 + 2 I1 / b^2) - 1),  b^2 = 1e4

Stencil: certified sym = per-branch density, then the 1/2 + 1/2 weighted
average of the DENSITIES (a quartic does not commute with that average;
the alternative reading is audited explicitly in mode `robust`).

Legendre: I(omega) = A + C2 omega^2 + C4 omega^4  ->  H = C2 omega^2 +
3 C4 omega^4 - A, and for the certified control L = -4 I1 - V4 the
extensive inertia is kin = -4 C2(I1).

Modes (each writes a checkpoint under ../checkpoints/m5_32_r8/):
    c1      degree-6 sampling of every term, odd + degree-6 content (C1)
    ladder  the three-box ladder, A / C2 / C4 exponents (C2, C3, C4, C5)
    shell   shell integrands, the z-axis string, h-refinement (C4, C6)
    robust  branch read, boundary overlap, fit window, mutations
    theorem generators, degenerate vacua, the annihilation cost (C7)
    tail    jets / deviation / spectrum / a0 in the far field (C8)
    coeff   the coefficient ladder, negative c5, delta ladder (C9, C10)
    merge   assembles ../data/m5_32_r8_audit.json
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
CKPT = os.path.join(HERE, "..", "checkpoints", "m5_32_r8")
os.makedirs(CKPT, exist_ok=True)
T0 = time.time()


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


B3 = _load("m5_21_3_a_4d", "m5_21_3_a_4d.py")            # certified stack
B8 = _load("m5_21_8_b_lattice", "m5_21_8_b_lattice.py")  # hedgehog + clock

ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
W = np.diag(ETA)
WW = W[:, None] * W[None, :]
G, DELTA, S = 32.0, 0.3, -1.0
B2_BI = 1.0e4
TERMS = ("Q_I1sq", "Q_I4sq", "Q_Fpair", "Q_C6a", "Q_C6b", "Q_BI", "I1")
# MY omega grid: 9 samples (the producer used 5), symmetric, degree-6 capable
OMEGAS = np.array([-1.5, -1.0, -0.5, -0.25, 0.0, 0.25, 0.5, 1.0, 1.5])
BRANCHES = (("bwd", 0.5), ("fwd", 0.5))          # MY branch order (reversed)


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def dump(tag, obj):
    with open(os.path.join(CKPT, f"{tag}.json"), "w") as f:
        json.dump(obj, f, indent=1, default=float)
    log(f"checkpoint {tag}.json written")


def rel(a, b):
    d = max(abs(a), abs(b), 1e-300)
    return abs(a - b) / d


def loglog(xs, ys):
    xs = np.asarray(xs, float)
    ys = np.asarray(ys, float)
    m = ys != 0
    if m.sum() < 2:
        return float("nan")
    return float(np.polyfit(np.log(xs[m]), np.log(np.abs(ys[m])), 1)[0])


def slope_err(xs, ys):
    """least-squares slope of log y vs log x with its standard error."""
    x = np.log(np.asarray(xs, float))
    y = np.log(np.abs(np.asarray(ys, float)))
    n = len(x)
    p, res, *_ = np.polyfit(x, y, 1, full=True)
    if n <= 2 or len(res) == 0:
        return float(p[0]), float("nan")
    s2 = res[0] / (n - 2)
    se = float(np.sqrt(s2 / np.sum((x - x.mean()) ** 2)))
    return float(p[0]), se


# ===================== MY differencing / algebra =====================
def mydiff(f, ax, h, br):
    out = np.zeros_like(f)
    sl = [slice(None)] * f.ndim

    def at(i):
        s = list(sl); s[ax] = i; return tuple(s)
    d = (f[at(slice(1, None))] - f[at(slice(0, -1))]) / h
    if br == "fwd":
        out[at(slice(0, -1))] = d
    elif br == "bwd":
        out[at(slice(1, None))] = d
    else:
        raise ValueError(br)
    return out


def aeta_b(A, Bm):
    """A eta B, batched."""
    return (A * W[None, :]) @ Bm


def comm_eta(A, Bm):
    return aeta_b(A, Bm) - aeta_b(Bm, A)


def inner_eta(F, Gm):
    """<F,G>_eta = tr(eta F eta G^T) = sum_ab w_a w_b F_ab G_ab."""
    return np.einsum("ab,...ab,...ab->...", WW, F, Gm, optimize=True)


def tr_AetaBeta(A, Bm):
    """tr(A eta B eta) = sum_ab w_a w_b A_ab B_ba."""
    return np.einsum("ab,...ab,...ba->...", WW, A, Bm, optimize=True)


# ===================== the term densities (mine) =====================
PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]


def densities(A, terms=TERMS, b2=B2_BI):
    """per-cell densities of every term from A = [A_0, A_1, A_2, A_3]."""
    Fs = {}
    for (mu, nu) in PAIRS:
        Fs[(mu, nu)] = comm_eta(A[mu], A[nu])

    def getF(mu, nu):
        if mu == nu:
            return None
        return Fs[(mu, nu)] if (mu, nu) in Fs else -Fs[(nu, mu)]

    out = {}
    need_I1 = ("I1" in terms) or ("Q_I1sq" in terms) or ("Q_BI" in terms)
    if need_I1:
        i1 = 0.0
        for (mu, nu) in PAIRS:
            i1 = i1 + W[mu] * W[nu] * inner_eta(Fs[(mu, nu)], Fs[(mu, nu)])
        if "I1" in terms:
            out["I1"] = i1
        if "Q_I1sq" in terms:
            out["Q_I1sq"] = i1 * i1
        if "Q_BI" in terms:
            out["Q_BI"] = b2 * (np.sqrt(np.maximum(1.0 + 2.0 * i1 / b2, 0.0)) - 1.0)
    if "Q_I4sq" in terms:
        sh = A[0].shape[:-2]
        R = np.zeros(sh + (4, 4))
        for nu in range(4):
            for mu in range(4):
                if mu == nu:
                    continue
                R[..., nu, :] += getF(mu, nu)[..., :, mu]
        i4 = np.einsum("na,...na,...na->...", WW, R, R, optimize=True)
        out["Q_I4sq"] = i4 * i4
    if "Q_Fpair" in terms:
        tot = 0.0
        for (mu, nu) in PAIRS:
            wp = W[mu] * W[nu]
            for (rh, sg) in PAIRS:
                wq = W[rh] * W[sg]
                v = inner_eta(Fs[(mu, nu)], Fs[(rh, sg)])
                tot = tot + wp * wq * v * v
        out["Q_Fpair"] = tot
    if ("Q_C6a" in terms) or ("Q_C6b" in terms):
        Tm = {}
        for mu in range(4):
            for nu in range(mu, 4):
                Tm[(mu, nu)] = tr_AetaBeta(A[mu], A[nu])
        if "Q_C6a" in terms:
            s = 0.0
            for mu in range(4):
                s = s + W[mu] * Tm[(mu, mu)]
            out["Q_C6a"] = s * s
        if "Q_C6b" in terms:
            tot = 0.0
            for mu in range(4):
                for nu in range(4):
                    t = Tm[(mu, nu)] if (mu, nu) in Tm else Tm[(nu, mu)]
                    tot = tot + W[mu] * W[nu] * t * t
            out["Q_C6b"] = tot
    return out


# ===================== fitting machinery =====================
def _fitmat(omegas, deg):
    V = np.vander(np.asarray(omegas, float), deg + 1, increasing=True)
    return V, np.linalg.pinv(V)


def sample_box(n, L, delta=DELTA, g=G, omegas=OMEGAS, terms=TERMS,
               branch_read="per_branch", a0_scale=None, a0_override=None,
               eta_mutant=False):
    """per-cell density samples I_term(omega, x), h^3 already folded in.

    branch_read: 'per_branch' = certified (density per branch, then the
    weighted average); 'avg_jets' = average the SPATIAL JETS first, then
    one density (the alternative reading audited in C4-a)."""
    cfg = B3.base_cfg(s=S, g=g, n=n, L=L, delta=delta)
    M = B8.dressed(cfg, 0.0)
    a0 = B8.a0_unit(cfg, 0.0) if a0_override is None else a0_override
    if a0_scale is not None:
        a0 = a0 * a0_scale[..., None, None]
    h3 = cfg["h"] ** 3
    nom = len(omegas)
    vals = {t: np.zeros((nom,) + M.shape[:-2]) for t in terms}
    if branch_read == "avg_jets":
        Asp = [sum(wt * mydiff(M, ax, cfg["h"], br) for br, wt in BRANCHES)
               for ax in range(3)]
        brs = [(Asp, 1.0)]
    else:
        brs = [([mydiff(M, ax, cfg["h"], br) for ax in range(3)], wt)
               for br, wt in BRANCHES]
    global WW, W
    W_save, WW_save = W, WW
    if eta_mutant:                       # MUTATION: all-plus internal metric
        W = np.ones(4)
        WW = W[:, None] * W[None, :]
    try:
        for Asp, wt in brs:
            for k, om in enumerate(omegas):
                d = densities([om * a0] + list(Asp), terms=terms)
                for t in terms:
                    vals[t][k] += wt * d[t]
    finally:
        W, WW = W_save, WW_save
    for t in terms:
        vals[t] *= h3
    return cfg, M, a0, vals


def coeffs_of(vals, omegas=OMEGAS, deg=6):
    """per-cell polynomial coefficients c_0..c_deg of I(omega, x)."""
    V, P = _fitmat(omegas, deg)
    flat = vals.reshape(len(omegas), -1)
    c = P @ flat
    resid = np.abs(V @ c - flat)
    return c.reshape((deg + 1,) + vals.shape[1:]), float(resid.max())


def integrals(vals):
    return vals.reshape(len(vals), -1).sum(axis=1)


# =============================== C1 ===============================
def audit_c1(n=32, L=48.0):
    """MORE omegas than the producer (9 vs 5) and a degree-6 fit."""
    out = {"n": n, "L": L, "omegas": OMEGAS.tolist(), "deg": 6, "terms": {}}
    cfg, M, a0, vals = sample_box(n, L)
    for t in TERMS:
        I = integrals(vals[t])
        c6, res6 = coeffs_of(vals[t], deg=6)
        c4, res4 = coeffs_of(vals[t], deg=4)
        C = [float(c6[p].sum()) for p in range(7)]
        Ceven = [float(c4[p].sum()) for p in range(5)]
        # odd content read off the INTEGRAL, the way the producer would
        odd = 0.0
        for k, om in enumerate(OMEGAS):
            if om <= 0:
                continue
            j = int(np.argmin(np.abs(OMEGAS + om)))
            odd = max(odd, abs(I[k] - I[j]) / max(abs(I[k]), abs(I[j]), 1e-300))
        scale = max(abs(v) for v in I)
        # the 5-sample degree-4 extraction the producer used
        sub = np.array([0, 2, 4, 6, 8])          # -1.5,-0.5,0,0.5,1.5
        V5, P5 = _fitmat(OMEGAS[sub], 4)
        c5 = P5 @ vals[t].reshape(len(OMEGAS), -1)[sub]
        C4_5samp = float(c5[4].sum())
        out["terms"][t] = {
            "I_of_omega": I.tolist(),
            "A": C[0], "C1_odd": C[1], "C2": C[2], "C3_odd": C[3],
            "C4": C[4], "C5_odd": C[5], "C6_deg6": C[6],
            "C4_deg4fit": Ceven[4], "C4_5sample_deg4": C4_5samp,
            "odd_rel_integral": float(odd),
            "deg6_rel_to_C4": abs(C[6]) / max(abs(C[4]), 1e-300),
            "deg6_rel_to_scale": abs(C[6]) / max(scale, 1e-300),
            "resid_deg6_max_percell": res6,
            "resid_deg4_max_percell": res4,
            "resid_deg4_rel": res4 / max(np.abs(vals[t]).max(), 1e-300),
            "C4_5samp_vs_deg6_rel": rel(C4_5samp, C[4])}
    # MUTATION: inject a known odd part and a known degree-6 part into the
    # SAMPLES; both detectors must fire at order 1
    v = vals["Q_I1sq"].copy()
    scale = np.abs(v).max()
    v_odd = v + 0.05 * scale * OMEGAS[:, None, None, None] ** 3
    v_d6 = v + 0.05 * scale * OMEGAS[:, None, None, None] ** 6
    Io = integrals(v_odd)
    odd_mut = max(abs(Io[k] - Io[len(OMEGAS) - 1 - k])
                  / max(abs(Io[k]), 1e-300) for k in range(4))
    c_od, _ = coeffs_of(v_odd, deg=6)
    c_d6, _ = coeffs_of(v_d6, deg=6)
    out["mutation_detectors"] = {
        "odd_injected_rel": float(odd_mut),
        "odd_coeff_C3_recovered": float(c_od[3].sum()),
        "deg6_injected_C6_recovered": float(c_d6[6].sum()),
        "deg6_expected": float(0.05 * scale * v[0].size),
        "clean_C3": out["terms"]["Q_I1sq"]["C3_odd"],
        "clean_C6": out["terms"]["Q_I1sq"]["C6_deg6"],
        "detectors_fire": bool(odd_mut > 1e-3
                               and abs(float(c_d6[6].sum())) > 1e-3)}
    out["worst_odd_rel"] = max(v["odd_rel_integral"] for v in out["terms"].values())
    out["worst_deg6_rel_to_scale"] = max(
        v["deg6_rel_to_scale"] for v in out["terms"].values())
    out["worst_deg6_rel_to_C4_polynomial_terms"] = max(
        v["deg6_rel_to_C4"] for k, v in out["terms"].items() if k != "I1")
    dump("c1_degree", out)
    return out


# ========================= the box ladder =========================
BOXES = ((32, 48.0), (48, 72.0), (64, 96.0))


def audit_ladder(boxes=BOXES, branch_read="per_branch", tag="ladder"):
    out = {"branch_read": branch_read, "boxes": [], "terms": {}}
    store = {t: {"A": [], "C2": [], "C4": [], "C6": []} for t in TERMS}
    Ls = []
    for (n, L) in boxes:
        t0 = time.time()
        cfg, M, a0, vals = sample_box(n, L, branch_read=branch_read)
        Ls.append(L)
        row = {"n": n, "L": L, "h": cfg["h"]}
        for t in TERMS:
            c, _ = coeffs_of(vals[t], deg=6)
            for k, p in (("A", 0), ("C2", 2), ("C4", 4), ("C6", 6)):
                store[t][k].append(float(c[p].sum()))
        row["wall_s"] = round(time.time() - t0, 2)
        out["boxes"].append(row)
        log(f"{tag}: box n={n} L={L} done ({row['wall_s']}s) "
            f"I1 C2 {store['I1']['C2'][-1]:.5f} Q_I1sq C4 {store['Q_I1sq']['C4'][-1]:.5f}")
        dump(f"{tag}_partial", {"out": out, "store": store, "L": Ls})
    for t in TERMS:
        s = store[t]
        out["terms"][t] = {
            "L": Ls, "A": s["A"], "C2": s["C2"], "C4": s["C4"], "C6": s["C6"],
            "A_exponent_in_L": loglog(Ls, s["A"]),
            "C2_exponent_in_L": loglog(Ls, s["C2"]),
            "C4_exponent_in_L": loglog(Ls, s["C4"]),
            "C4_pair_exponents": [
                float(np.log(abs(s["C4"][i + 1] / s["C4"][i]))
                      / np.log(Ls[i + 1] / Ls[i])) if s["C4"][i] else None
                for i in range(len(Ls) - 1)],
            "kin_from_C2": [-4.0 * v for v in s["C2"]]}
    out["I1_control"] = {
        "kin_L48": -4.0 * store["I1"]["C2"][0],
        "certified_kin_L48": 426.5070121483972,
        "rel_dev": rel(-4.0 * store["I1"]["C2"][0], 426.5070121483972),
        "C4_of_I1_L48": store["I1"]["C4"][0]}
    dump(tag, out)
    return out


def audit_ladder_mutation():
    """MUTATION: the all-plus internal metric must break the I1 control."""
    cfg, M, a0, vals = sample_box(32, 48.0, terms=("I1",), eta_mutant=True)
    c, _ = coeffs_of(vals["I1"], deg=6)
    out = {"mutant": "internal eta -> identity in <,>_eta and comm_eta",
           "kin_mutant": -4.0 * float(c[2].sum()),
           "certified_kin": 426.5070121483972}
    out["rel_dev"] = rel(out["kin_mutant"], out["certified_kin"])
    dump("ladder_mutation", out)
    return out


# ===================== the crux: where the C4 lives =====================
def shell_profile(dens, r, h, rmax, dr=1.5, mask=None):
    d = dens if mask is None else np.where(mask, dens, 0.0)
    edges = np.arange(0.0, rmax + dr, dr)
    idx = np.digitize(r.ravel(), edges) - 1
    tot = np.bincount(idx, weights=d.ravel(), minlength=len(edges))
    cnt = np.bincount(idx, minlength=len(edges))
    mid = 0.5 * (edges[:-1] + edges[1:])
    return mid, tot[:len(mid)] / dr, cnt[:len(mid)]


def audit_shell(n=64, L=96.0):
    out = {"n": n, "L": L}
    cfg, M, a0, vals = sample_box(n, L)
    h = cfg["h"]
    X, Y, Z = B3.coords(n, h)
    r = np.sqrt(X * X + Y * Y + Z * Z)
    rho = np.sqrt(X * X + Y * Y)                 # distance to the z axis
    tube = rho < 1.01 * h                        # the 4 near-axis columns
    coef = {}
    for t in ("I1", "Q_I1sq", "Q_I4sq", "Q_Fpair", "Q_C6a", "Q_C6b"):
        c, _ = coeffs_of(vals[t], deg=6)
        coef[t] = {"A": c[0].copy(), "C2": c[2].copy(), "C4": c[4].copy()}
    # the net omega^2 density of I1 and its square (the C6 mechanism)
    c2 = coef["I1"]["C2"]
    out["sign_definite"] = {
        "c2_max": float(c2.max()), "c2_min": float(c2.min()),
        "frac_cells_positive": float((c2 > 0).mean()),
        "all_le_zero": bool(c2.max() <= 0.0)}
    # shells, restricted to the fully-inside region r < L/2
    rmax = L / 2.0
    fields = {"net_c2": c2, "c2_squared": c2 * c2,
              "c2_Q_I1sq": coef["Q_I1sq"]["C2"],
              "c2_Q_Fpair": coef["Q_Fpair"]["C2"],
              "static_Q_I1sq": coef["Q_I1sq"]["A"],
              "c4_Q_I1sq": coef["Q_I1sq"]["C4"],
              "c4_Q_I4sq": coef["Q_I4sq"]["C4"],
              "c4_Q_Fpair": coef["Q_Fpair"]["C4"],
              "c4_Q_C6a": coef["Q_C6a"]["C4"],
              "static_I1": coef["I1"]["A"]}
    out["shells"] = {}
    for nm, f in fields.items():
        mid, per_dr, cnt = shell_profile(f, r, h, rmax)
        midT, per_drT, _ = shell_profile(f, r, h, rmax, mask=tube)
        midB, per_drB, _ = shell_profile(f, r, h, rmax, mask=~tube)
        win = (mid > 6.0) & (mid < 0.85 * rmax)
        sl, se = slope_err(mid[win], per_dr[win])
        slB, seB = slope_err(midB[win], per_drB[win])
        slT, seT = slope_err(midT[win], per_drT[win])
        # a wider and a narrower window (the fit-window attack)
        w2 = (mid > 3.0) & (mid < rmax)
        w3 = (mid > 12.0) & (mid < 0.85 * rmax)
        out["shells"][nm] = {
            "r_mid": mid.tolist(), "per_unit_r": per_dr.tolist(),
            "per_unit_r_off_axis": per_drB.tolist(),
            "per_unit_r_axis_tube": per_drT.tolist(),
            "cells_per_shell": cnt.tolist(),
            "slope": sl, "slope_stderr": se,
            "slope_off_axis": slB, "slope_off_axis_stderr": seB,
            "slope_axis_tube": slT, "slope_axis_tube_stderr": seT,
            "slope_window_3_to_wall": slope_err(mid[w2], per_dr[w2])[0],
            "slope_window_12_to_085": slope_err(mid[w3], per_dr[w3])[0],
            "total": float(f.sum()),
            "total_off_axis": float(f[~tube].sum()),
            "axis_tube_fraction": float(f[tube].sum() / f.sum()) if f.sum() else None}
    out["axis_string"] = {
        "note": ("rho < h picks the 4 lattice columns adjacent to the z axis, "
                 "where the Euler-angle frame of the hedgehog is discontinuous"),
        "n_tube_cells": int(tube.sum()), "n_cells": int(tube.size),
        "tube_cell_fraction": float(tube.mean())}
    # how wide is the string? exclude a growing tube and watch the slope
    sweep = {}
    win = None
    for wid in (0.0, 1.01, 2.01, 3.01, 4.01, 6.01):
        msk = rho >= wid * h
        row = {}
        for nm in ("c2_squared", "c4_Q_I1sq", "c4_Q_Fpair", "c4_Q_I4sq"):
            f = fields[nm]
            mid, per_dr, _ = shell_profile(f, r, h, rmax, mask=msk)
            if win is None:
                win = (mid > 6.0) & (mid < 0.85 * rmax)
            sl, se = slope_err(mid[win], per_dr[win])
            row[nm] = {"slope": sl, "stderr": se,
                       "total": float(f[msk].sum()),
                       "frac_of_full": float(f[msk].sum() / f.sum())}
        row["cells_kept_frac"] = float(msk.mean())
        sweep[f"exclude_rho_lt_{wid:g}h"] = row
    out["tube_sweep"] = sweep
    out["naive_prediction"] = (
        "|c2| falls as r^-2 (jets ~ 1/r, a0 flat), so c2^2 ~ r^-4 and its "
        "shell integrand ~ r^-2, which CONVERGES; the measured slow decay is "
        "carried by the z-axis coordinate string, not by the bulk tail")
    dump("shell", out)
    return out


def audit_string(zs=(5.0, 10.0, 20.0, 40.0), rhos=(1.0, 0.1, 0.01, 1e-4)):
    """is the hedgehog frame really discontinuous on the z axis?
    Sample M on a phi ring of shrinking radius at fixed z: if the spread
    over the ring does not go to zero, M has a line discontinuity."""
    d4 = np.diag([G, 1.0, DELTA, 0.0])
    rows = []
    for z in zs:
        for rr in rhos:
            Ms = []
            for phi in np.linspace(0.0, 2 * np.pi, 17)[:-1]:
                x, y = rr * np.cos(phi), rr * np.sin(phi)
                ph = np.arctan2(y, x)
                th = -np.arctan2(z, np.sqrt(x * x + y * y))
                Q = np.einsum(
                    "...ab,...bc->...ac",
                    B8.rot_field(B8.G3, np.full((1, 1, 1), ph)),
                    B8.rot_field(B8.G2, np.full((1, 1, 1), th)))
                Ms.append(np.einsum("...ab,bc,...dc->...ad", Q, d4, Q)[0, 0, 0])
            Ms = np.array(Ms)
            rows.append({"z": z, "rho": rr,
                         "ring_spread_max": float(np.abs(Ms - Ms.mean(0)).max()),
                         "ring_spread_rms": float(np.sqrt(
                             ((Ms - Ms.mean(0)) ** 2).sum(axis=(1, 2)).mean()))})
    out = {"rows": rows,
           "spread_at_smallest_rho": [r for r in rows if r["rho"] == rhos[-1]],
           "note": ("a spread that does NOT shrink with rho proves the ansatz "
                    "M = Qh d4 Qh^T is discontinuous along the whole z axis, "
                    "so the jets there are set by the lattice spacing, not by r")}
    dump("string", out)
    return out


def audit_h_refine(L=48.0, ns=(32, 64, 96)):
    """halve h at FIXED L: an IR quantity must be h-stable; the C4 of a
    line singularity diverges as h -> 0."""
    out = {"L": L, "rows": []}
    for n in ns:
        cfg, M, a0, vals = sample_box(n, L, terms=("I1", "Q_I1sq", "Q_I4sq",
                                                   "Q_Fpair", "Q_C6a"))
        h = cfg["h"]
        X, Y, Z = B3.coords(n, h)
        rho = np.sqrt(X * X + Y * Y)
        tube = rho < 1.01 * h
        row = {"n": n, "h": h}
        for t in ("I1", "Q_I1sq", "Q_I4sq", "Q_Fpair", "Q_C6a"):
            c, _ = coeffs_of(vals[t], deg=6)
            row[f"{t}_A"] = float(c[0].sum())
            row[f"{t}_C2"] = float(c[2].sum())
            row[f"{t}_C4"] = float(c[4].sum())
            row[f"{t}_C4_off_axis"] = float(c[4][~tube].sum())
        row["kin"] = -4.0 * row["I1_C2"]
        out["rows"].append(row)
        log(f"h-refine n={n} h={h}: kin {row['kin']:.4f} "
            f"Q_I1sq C4 {row['Q_I1sq_C4']:.5f}")
        dump("h_refine_partial", out)
    hs = [r["h"] for r in out["rows"]]
    out["exponents_in_h"] = {}
    for k in ("kin", "Q_I1sq_C4", "Q_I4sq_C4", "Q_Fpair_C4", "Q_C6a_C4",
              "Q_I1sq_C4_off_axis", "Q_I4sq_C4_off_axis"):
        out["exponents_in_h"][k] = loglog(hs, [r[k] for r in out["rows"]])
    dump("h_refine", out)
    return out


def audit_overlap(nsmall=32, nbig=64, h=1.5):
    """R7-N2: the same cells in two boxes. If the shared interior cells
    agree, the L growth is added volume, not a boundary artifact."""
    out = {}
    cfg1, _, _, v1 = sample_box(nsmall, nsmall * h,
                                terms=("I1", "Q_I1sq", "Q_C6a"))
    cfg2, _, _, v2 = sample_box(nbig, nbig * h,
                                terms=("I1", "Q_I1sq", "Q_C6a"))
    o = (nbig - nsmall) // 2
    X1, Y1, Z1 = B3.coords(nsmall, h)
    X2, Y2, Z2 = B3.coords(nbig, h)
    out["grid_aligned_max_dev"] = float(
        np.abs(X2[o:o + nsmall, o:o + nsmall, o:o + nsmall] - X1).max())
    interior = np.zeros((nsmall,) * 3, dtype=bool)
    interior[1:-1, 1:-1, 1:-1] = True
    for t in ("I1", "Q_I1sq", "Q_C6a"):
        c1, _ = coeffs_of(v1[t], deg=6)
        c2, _ = coeffs_of(v2[t], deg=6)
        row = {}
        for k, p in (("C2", 2), ("C4", 4)):
            a = c1[p]
            b = c2[p][o:o + nsmall, o:o + nsmall, o:o + nsmall]
            row[f"{k}_max_abs_diff_interior"] = float(np.abs(a - b)[interior].max())
            row[f"{k}_rel_dev_sum_interior"] = float(
                abs(b[interior].sum() - a[interior].sum())
                / max(abs(a[interior].sum()), 1e-300))
            row[f"{k}_small_interior_sum"] = float(a[interior].sum())
            row[f"{k}_small_total"] = float(a.sum())
            row[f"{k}_big_total"] = float(c2[p].sum())
            row[f"{k}_big_outside_small"] = float(c2[p].sum() - b.sum())
        out[t] = row
    dump("overlap", out)
    return out


def audit_robust():
    out = {}
    out["branch_read_avg_jets"] = audit_ladder(
        branch_read="avg_jets", tag="ladder_avgjets")["terms"]
    out["overlap"] = audit_overlap()
    out["h_refine"] = audit_h_refine()
    dump("robust", out)
    return out


# =============================== C7 ===============================
def gen_boost(k):
    K = np.zeros((4, 4)); K[0, k] = K[k, 0] = 1.0
    return K


def gen_rot(i, j):
    J = np.zeros((4, 4)); J[i, j] = -1.0; J[j, i] = 1.0
    return J


GENS = {"rot_12": gen_rot(1, 2), "rot_13": gen_rot(1, 3), "rot_23": gen_rot(2, 3),
        "boost_1": gen_boost(1), "boost_2": gen_boost(2), "boost_3": gen_boost(3)}


def tangent(Xg, M):
    """the ONLY tangent of M -> L M L^T: X M + M X^T (symmetric)."""
    return Xg @ M + M @ Xg.T


def tangent_wrong(Xg, M):
    """the map X M - M X^T (what reproduces the producer's numbers)."""
    return Xg @ M - M @ Xg.T


def audit_theorem():
    out = {"vacuum": [G, 1.0, DELTA, 0.0], "generators": {}}
    Mv = np.diag([G, 1.0, DELTA, 0.0])
    for nm, Xg in GENS.items():
        t = tangent(Xg, Mv)
        tw = tangent_wrong(Xg, Mv)
        out["generators"][nm] = {
            "tangent_max_abs": float(np.abs(t).max()),
            "tangent_is_symmetric": bool(np.abs(t - t.T).max() < 1e-14),
            "annihilates": bool(np.abs(t).max() < 1e-14),
            "wrong_map_max_abs": float(np.abs(tw).max()),
            "wrong_map_is_symmetric": bool(np.abs(tw - tw.T).max() < 1e-14),
            "tangent": t.tolist()}
    out["algebra"] = (
        "for diagonal M and a rotation J_ij (J_ij = -1, J_ji = +1) the tangent "
        "J M + M J^T = J M - M J has both nonzero entries equal to m_i - m_j, "
        "so it vanishes iff m_i = m_j; for a boost K_k the tangent K M + M K^T "
        "= {K, M} has entries m_0 + m_k, which vanishes iff m_k = -m_0, "
        "impossible for m_0 = 32 > 0 and m_k >= 0")
    out["min_over_generators"] = min(
        v["tangent_max_abs"] for v in out["generators"].values())
    out["producer_numbers_reproduced_by"] = (
        "X M - M X^T, which is ANTISYMMETRIC for diagonal M and therefore is "
        "not in the tangent space of the symmetric field M at all")
    # the CONVERSE: build degenerate vacua and check the annihilation
    conv = []
    for (dl, plane) in ((1.0, "rot_12"), (0.0, "rot_23")):
        Mv2 = np.diag([G, 1.0, dl, 0.0])
        t = tangent(GENS[plane], Mv2)
        conv.append({"delta": dl, "plane": plane,
                     "vacuum": np.diag(Mv2).tolist(),
                     "tangent_max_abs": float(np.abs(t).max()),
                     "annihilates": bool(np.abs(t).max() < 1e-14)})
    out["converse"] = conv
    # what the escape COSTS on the actual hedgehog
    cost = []
    for dl in (0.0, 0.05, 0.3, 1.0):
        cfg = B3.base_cfg(s=S, g=G, n=32, L=48.0, delta=dl)
        M = B8.dressed(cfg, 0.0)
        a0 = B8.a0_unit(cfg, 0.0)                    # the REALIZED clock
        eu, ev = B3.e_parts(M, cfg)
        X, Y, Z = B3.coords(32, cfg["h"])
        r = np.sqrt(X * X + Y * Y + Z * Z)
        Mvac = B3.vac4(cfg)
        dev = np.sqrt(((M - Mvac) ** 2).sum(axis=(-1, -2)))
        row = {"delta": dl, "E_u": float(eu), "E_V4": float(ev),
               "kin_certified": float(B3.kin_of(M, a0, cfg)),
               "a0_norm_max": float(np.linalg.norm(a0, axis=(-2, -1)).max()),
               "dev_from_vac_rms_r12": float(dev[(r > 12) & (r < 15)].mean()),
               "dev_from_vac_rms_r21": float(dev[(r > 21) & (r < 24)].mean())}
        # the CONJUGATED annihilating generator on the hedgehog:
        # a0_conj = Q (X d4 + d4 X^T) Q^T, exactly zero when X annihilates d4
        for plane in ("rot_12", "rot_23"):
            Xg = GENS[plane]
            d4 = B3.vac4(cfg)
            core = tangent(Xg, d4)
            row[f"conj_{plane}_core_max_abs"] = float(np.abs(core).max())
            # the FIXED (unconjugated) generator on the hedgehog
            a0f = Xg @ M + M @ np.swapaxes(np.broadcast_to(Xg, M.shape), -1, -2)
            nf = np.linalg.norm(a0f, axis=(-2, -1))
            row[f"fixed_{plane}_a0_r12"] = float(nf[(r > 12) & (r < 15)].mean())
            row[f"fixed_{plane}_a0_r21"] = float(nf[(r > 21) & (r < 24)].mean())
        cost.append(row)
    out["degenerate_cost"] = cost
    out["escape_theorem"] = (
        "M(x) = Q(x) d4 Q(x)^T with Q in SO(1,3) satisfies Q^T (Q^-1)^T = I, so "
        "for the CO-MOVING generator X_loc = Q X Q^-1 the clock flow is "
        "Q (X d4 + d4 X^T) Q^T: if X annihilates the vacuum it annihilates the "
        "whole hedgehog and kin = 0 (no clock); if X does not annihilate the "
        "vacuum the flow is a nonzero constant at infinity. Degeneracy removes "
        "the obstruction only by removing the clock.")
    dump("theorem", out)
    return out


# =============================== C8 ===============================
def audit_tail(n=64, L=96.0):
    cfg = B3.base_cfg(s=S, g=G, n=n, L=L, delta=DELTA)
    M = B8.dressed(cfg, 0.0)
    a0 = B8.a0_unit(cfg, 0.0)
    h = cfg["h"]
    X, Y, Z = B3.coords(n, h)
    r = np.sqrt(X * X + Y * Y + Z * Z)
    Mvac = B3.vac4(cfg)
    A = [mydiff(M, ax, h, "fwd") for ax in range(3)]
    jet = np.sqrt(sum((a * a).sum(axis=(-1, -2)) for a in A))
    dev = np.sqrt(((M - Mvac) ** 2).sum(axis=(-1, -2)))
    an = np.linalg.norm(a0, axis=(-2, -1))
    spec = np.linalg.eigvals(M @ ETA).real
    spec = np.sort(spec, axis=-1)
    vspec = np.sort(np.linalg.eigvals(Mvac @ ETA).real)
    rows = []
    for rc in np.arange(6.0, L / 2.0, 3.0):
        m = (r >= rc) & (r < rc + 3.0)
        if m.sum() < 20:
            continue
        rows.append({"r": float(rc + 1.5),
                     "jet_rms": float(np.sqrt((jet[m] ** 2).mean())),
                     "dev_rms": float(np.sqrt((dev[m] ** 2).mean())),
                     "a0_rms": float(np.sqrt((an[m] ** 2).mean())),
                     "spec_max_dev_from_vacuum":
                         float(np.abs(spec[m] - vspec).max())})
    rr = [x["r"] for x in rows]
    out = {"n": n, "L": L, "rows": rows,
           "vacuum_spectrum_of_M_eta": vspec.tolist(),
           "jet_exponent": loglog(rr, [x["jet_rms"] for x in rows]),
           "dev_exponent": loglog(rr, [x["dev_rms"] for x in rows]),
           "a0_exponent": loglog(rr, [x["a0_rms"] for x in rows]),
           "dev_rms_level": rows[0]["dev_rms"],
           "a0_level": rows[0]["a0_rms"],
           "delta_sqrt2": float(DELTA * np.sqrt(2.0)),
           "spec_max_dev_over_all_shells":
               max(x["spec_max_dev_from_vacuum"] for x in rows)}
    # MUTATION: a field that DOES relax to the constant vacuum must show a
    # strongly negative dev exponent (the detector has to be able to see decay)
    env = np.exp(-(r / 12.0))
    Mm = Mvac + env[..., None, None] * (M - Mvac)
    devm = np.sqrt(((Mm - Mvac) ** 2).sum(axis=(-1, -2)))
    dm = []
    for rc in np.arange(6.0, L / 2.0, 3.0):
        m = (r >= rc) & (r < rc + 3.0)
        if m.sum() < 20:
            continue
        dm.append(float(np.sqrt((devm[m] ** 2).mean())))
    out["mutation_decaying_field_dev_exponent"] = loglog(rr[:len(dm)], dm)
    dump("tail", out)
    return out


# =========================== C9 and C10 ===========================
def audit_coeff(lad):
    """the coefficient ladder from MY OWN C2 / C4."""
    out = {"per_term": {}, "note":
           "H2 = -4 C2(I1) + c5 C2(q); flip needs c5 >= 4|C2_I1|/|C2_q| "
           "when C2_q < 0"}
    Ls = lad["terms"]["I1"]["L"]
    c2_i1 = lad["terms"]["I1"]["C2"]
    kin = [-4.0 * v for v in c2_i1]
    for t in ("Q_I1sq", "Q_Fpair", "Q_I4sq"):
        c2q = lad["terms"][t]["C2"]
        c4q = lad["terms"][t]["C4"]
        req = [4.0 * abs(a) / abs(b) for a, b in zip(c2_i1, c2q)]
        om = [float(np.sqrt(k / (2.0 * c))) for k, c in zip(kin, c4q)]
        out["per_term"][t] = {
            "L": Ls, "c5_required_to_flip_H2": req,
            "exponent_in_L": loglog(Ls, req),
            "certified_inertia_H2": kin,
            "quartic_C4": c4q,
            "omega_where_quartic_matches_at_c5_1": om,
            "radiation_window_omega_max": 0.786,
            "omega_match_over_window": [v / 0.786 for v in om],
            "C2q_sign": [float(np.sign(v)) for v in c2q]}
    # the reading the producer did not test: c5 < 0
    out["negative_c5"] = {
        "H2_effect": ("C2_q < 0 for every C5 term, so c5 < 0 makes H2 = "
                      "-4 C2_I1 + c5 C2_q MORE positive: no flip, ever"),
        "H4_effect": ("H4 = 3 c5 C4_q with C4_q > 0, so c5 < 0 makes the "
                      "energy's quartic NEGATIVE and H unbounded below at "
                      "large omega: the free-omega problem has no minimum"),
        "C4_signs": {t: float(np.sign(lad["terms"][t]["C4"][0]))
                     for t in ("Q_I1sq", "Q_I4sq", "Q_Fpair", "Q_C6a", "Q_BI")}}
    bi = lad["terms"]["Q_BI"]
    out["born_infeld"] = {
        "C4": bi["C4"], "C2": bi["C2"], "A": bi["A"],
        "C4_sign": float(np.sign(bi["C4"][0])),
        "C2_matches_I1_rel": [rel(a, b) for a, b in zip(bi["C2"], c2_i1)],
        "reading": ("Born-Infeld at b^2 = 1e4 is I1 plus a tiny NEGATIVE "
                    "quartic; it saturates the wrong way for a clock, and "
                    "3 C4 < 0 makes H turn over at large omega")}
    # the test the producer stated but did not solve: the ACTUAL fixed-J
    # omega*, from J = 2 C2_tot omega + 4 C4_tot omega^3 with
    # C2_tot = kin + c5 C2_q and C4_tot = c5 C4_q
    fj = {}
    for t in ("Q_I1sq", "Q_I4sq", "Q_Fpair"):
        c2q = lad["terms"][t]["C2"]
        c4q = lad["terms"][t]["C4"]
        rows = []
        for c5 in (0.0, 1.0, 10.0, 100.0, 300.0):
            om, omL = [], []
            for k, L in enumerate(Ls):
                c2t = kin[k] + c5 * c2q[k]
                c4t = c5 * c4q[k]
                J = 200.0
                rt = np.roots([4.0 * c4t, 0.0, 2.0 * c2t, -J]) if c4t else \
                    np.array([J / (2.0 * c2t)])
                rr = [float(z.real) for z in np.atleast_1d(rt)
                      if abs(np.imag(z)) < 1e-9 and z.real > 0]
                w = min(rr) if rr else float("nan")
                om.append(w); omL.append(w * L)
            rows.append({"c5": c5, "omega_star": om, "omega_star_times_L": omL,
                         "omega_star_L_spread":
                             (max(omL) / min(omL) if min(omL) > 0 else None),
                         "exponent_in_L": loglog(Ls, om)})
        fj[t] = rows
    out["fixed_J_omega_star"] = {
        "J": 200.0, "L": Ls, "rows": fj,
        "reading": ("omega* L constant across the ladder means omega* ~ 1/L; "
                    "the cubic term is negligible at every c5 that keeps "
                    "C2_tot > 0, and C2_tot <= 0 removes the fixed-J minimum")}
    # unboundedness under a negative quartic
    unb = []
    for c5 in (-1.0, 1.0):
        c2t = kin[0] + c5 * lad["terms"]["Q_I1sq"]["C2"][0]
        c4t = c5 * lad["terms"]["Q_I1sq"]["C4"][0]
        H = [float(c2t * w ** 2 + 3.0 * c4t * w ** 4) for w in (1.0, 10.0, 30.0)]
        unb.append({"c5": c5, "C2_tot": c2t, "H4_coeff": 3.0 * c4t,
                    "H_at_omega_1_10_30": H,
                    "unbounded_below": bool(c4t < 0)})
    out["negative_c5"]["measured"] = unb
    # does the conclusion survive at other delta?
    dl_rows = []
    for dl in (0.1, 0.3, 0.6, 1.0):
        cfg, M, a0, vals = sample_box(32, 48.0, delta=dl,
                                      terms=("I1", "Q_I1sq"))
        ci, _ = coeffs_of(vals["I1"], deg=6)
        cq, _ = coeffs_of(vals["Q_I1sq"], deg=6)
        c2i, c2q = float(ci[2].sum()), float(cq[2].sum())
        c4q = float(cq[4].sum())
        dl_rows.append({"delta": dl, "kin": -4.0 * c2i, "C2_Q_I1sq": c2q,
                        "C4_Q_I1sq": c4q,
                        "c5_required": 4.0 * abs(c2i) / abs(c2q),
                        "omega_match": float(np.sqrt(-4.0 * c2i / (2.0 * c4q)))})
        log(f"delta {dl}: kin {-4.0 * c2i:.3f} c5_req {dl_rows[-1]['c5_required']:.2f}")
    # does the L^1.03 growth of the required c5 survive at delta = 1?
    dl_lad = []
    for dl in (0.1, 1.0):
        req, om4 = [], []
        for (n, L) in BOXES:
            _, _, _, v = sample_box(n, L, delta=dl, terms=("I1", "Q_I1sq"))
            ci, _ = coeffs_of(v["I1"], deg=6)
            cq, _ = coeffs_of(v["Q_I1sq"], deg=6)
            req.append(4.0 * abs(float(ci[2].sum())) / abs(float(cq[2].sum())))
            om4.append(float(np.sqrt(-4.0 * float(ci[2].sum())
                                     / (2.0 * float(cq[4].sum())))))
        dl_lad.append({"delta": dl, "L": [L for _, L in BOXES],
                       "c5_required": req,
                       "exponent_in_L": loglog([L for _, L in BOXES], req),
                       "omega_match": om4,
                       "omega_match_over_window": [v / 0.786 for v in om4]})
        log(f"delta {dl} L-ladder: c5_req {req} exp {dl_lad[-1]['exponent_in_L']:.4f}")
    out["delta_L_ladder"] = dl_lad
    out["delta_ladder"] = {
        "rows": dl_rows,
        "c5_required_spread":
            max(r["c5_required"] for r in dl_rows)
            / min(r["c5_required"] for r in dl_rows),
        "kin_exponent_in_delta": loglog([r["delta"] for r in dl_rows],
                                        [r["kin"] for r in dl_rows]),
        "c5_required_exponent_in_delta": loglog(
            [r["delta"] for r in dl_rows], [r["c5_required"] for r in dl_rows])}
    dump("coeff", out)
    return out


# =============================== main ===============================
def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    if mode in ("c1", "all"):
        audit_c1()
    if mode in ("ladder", "all"):
        audit_ladder()
        audit_ladder_mutation()
    if mode in ("shell", "all"):
        audit_string()
        audit_shell()
    if mode in ("robust", "all"):
        audit_robust()
    if mode in ("theorem", "all"):
        audit_theorem()
    if mode in ("tail", "all"):
        audit_tail()
    if mode in ("coeff", "all"):
        with open(os.path.join(CKPT, "ladder.json")) as f:
            audit_coeff(json.load(f))
    if mode in ("merge", "all"):
        _load("merge", "m5_32_r8_audit_merge.py").main()
    log(f"mode {mode} done")


if __name__ == "__main__":
    main()
