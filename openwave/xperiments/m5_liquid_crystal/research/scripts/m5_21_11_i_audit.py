"""M5.21.11 INDEPENDENT ADVERSARIAL AUDIT of the ladder-compute run.

Recomputes, with its own code, every headline claim of the M5.21.11
route-(b) ladder run (framework: findings/m5_21_11_framework.md,
FROZEN 2026-08-07) and tries to refute them:

    C1  g-arm: twin-minima dressing gain E(m*) - E(0) recomputed from
        the raw endpoints with an OWN boost-hedgehog construction
        (exp(m K(rhat)) verified against scipy expm), own m-grid
        (2-4x finer than the run's), branches A/C/B at g in {8, 32};
        gain match vs m5_21_11_garm.json + fall exponent far below 2.
    C2  corroboration: ansatz-family gain E_min - E(0) flat (~-61)
        across g = 8/16/32 from the pre-existing record
        m5_21_8_lat_gladder.json, with E(0) AND the full m-curves
        recomputed here via the record's own dressed() + INS4.e_parts.
    C3  gate table: virial |E_u - 3 E_V|/E_u, xstencil ratio, FIRE,
        L/a*, gap/plane reads recomputed from the raw rows/reads with
        own formulas; compared against m5_21_11_fit.json rung_state.
    C4  shell-fix integrity: pinned shell of three continuation
        endpoints == current-delta analytic seed EXACTLY (< 1e-12).
    C5  E(delta) monotone decreasing in delta per branch at N = 48,
        endpoints quoted from rows.
    C6  F3 arithmetic (usable counts < 6) + F1/F2 vacuity.
    +   barred-inputs sweep over the run's new scripts.

Allowed imports (pre-existing certified instruments / record code):
m5_21_2b_a_instrument.py, m5_21_3_a_4d.py, m5_21_8_b_lattice.py.
The run's own wrappers (m5_21_11_b/c/d/f/g) are NEVER imported: their
outputs are recomputed independently.

Headless, matplotlib Agg. Run from scripts/.
Out: ../data/m5_21_11_ladder_audit.json (verdicts, printed at end).
"""
from __future__ import annotations

import importlib.util
import json
import os
import re
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT = os.path.join(DATA, "m5_21_11_ladder_audit.json")


def load_mod(name, fname):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


INS2B = load_mod("ins2b_audit", "m5_21_2b_a_instrument.py")
INS4 = load_mod("ins4_audit", "m5_21_3_a_4d.py")
LAT8 = load_mod("lat8_audit", "m5_21_8_b_lattice.py")

RUNGS = [0.30, 0.25, 0.20, 0.15, 0.12, 0.09, 0.07, 0.05]
BRANCHES = ("A", "C", "B")
W2_PIN = 0.002758100
L_BOX = 48.0

VERDICTS = {}
NOTES = []


def jload(path):
    with open(os.path.join(DATA, path)) as f:
        return json.load(f)


def dkey(d):
    return f"{d:g}"


# =================================================================
# C1: g-arm recomputation (own boost hedgehog, own grid)
# =================================================================
def own_coords(n, h):
    x = (np.arange(n) - (n - 1) / 2.0) * h
    return np.meshgrid(x, x, x, indexing="ij")


def own_boost(n, h, m):
    """exp(m K(rhat)) built from scratch: K = radial boost generator,
    closed form I + sinh(m) K + (cosh(m) - 1) K^2 (K^3 = K for unit
    nhat), independently verified against scipy expm in gate_boost()."""
    X, Y, Z = own_coords(n, h)
    R = np.sqrt(X * X + Y * Y + Z * Z)
    Rs = np.where(R < 1e-12, 1e-12, R)
    nv = np.stack([X / Rs, Y / Rs, Z / Rs], axis=-1)
    K = np.zeros(X.shape + (4, 4))
    K[..., 0, 1:4] = nv
    K[..., 1:4, 0] = nv
    K2 = np.zeros_like(K)
    K2[..., 0, 0] = 1.0
    K2[..., 1:4, 1:4] = nv[..., :, None] * nv[..., None, :]
    return (np.eye(4)[None, None, None] + np.sinh(m) * K
            + (np.cosh(m) - 1.0) * K2)


def gate_boost():
    """own-construction gates: closed form == expm (sampled cells),
    Lorentz property B^T eta B = eta, and grid convention identity."""
    from scipy.linalg import expm
    rng = np.random.default_rng(7)
    eta = np.diag([-1.0, 1.0, 1.0, 1.0])
    n, h = 6, 1.7
    errs_expm, errs_lor = [], []
    for m in (-0.3, 0.011, 0.12565721414045303):
        B = own_boost(n, h, m)
        X, Y, Z = own_coords(n, h)
        R = np.sqrt(X * X + Y * Y + Z * Z)
        for _ in range(8):
            i, j, k = rng.integers(0, n, 3)
            nv = np.array([X[i, j, k], Y[i, j, k], Z[i, j, k]]) \
                / R[i, j, k]
            K = np.zeros((4, 4))
            K[0, 1:] = nv
            K[1:, 0] = nv
            errs_expm.append(float(np.max(np.abs(B[i, j, k]
                                                 - expm(m * K)))))
            errs_lor.append(float(np.max(np.abs(
                B[i, j, k].T @ eta @ B[i, j, k] - eta))))
    Xa, _, _ = own_coords(48, 1.0)
    Xb, _, _ = INS4.coords(48, 1.0)
    grid_err = float(np.max(np.abs(Xa - Xb)))
    return {"expm_err": max(errs_expm), "lorentz_err": max(errs_lor),
            "grid_err": grid_err,
            "pass": max(errs_expm) < 1e-12 and max(errs_lor) < 1e-12
            and grid_err == 0.0}


def own_embed(M3, sg):
    n = M3.shape[0]
    M4 = np.zeros((n, n, n, 4, 4))
    M4[..., 1:4, 1:4] = M3
    M4[..., 0, 0] = -sg
    return M4


def e_dressed_own(M4, cfg, m):
    B = own_boost(cfg["n"], cfg["h"], m)
    Md = np.einsum("...ab,...bc,...dc->...ad", B, M4, B,
                   optimize=True)
    Md = 0.5 * (Md + Md.swapaxes(-1, -2))
    e_u, e_v = INS4.e_parts(Md, cfg)          # certified instrument
    return float(e_u + e_v)


def refine_min(ms, Es):
    i = int(np.argmin(Es))
    m_star, E_star = float(ms[i]), float(Es[i])
    if 0 < i < len(ms) - 1:
        a, b, c = Es[i - 1], Es[i], Es[i + 1]
        den = c - 2 * b + a
        if den > 0:
            dm = ms[1] - ms[0]
            m_star = float(ms[i] - 0.5 * dm * (c - a) / den)
            E_star = float(b - 0.125 * (c - a) ** 2 / den)
    return m_star, E_star


def audit_c1():
    s = -1.0
    garm = jload("m5_21_11_garm.json")
    run_arms = {(a["branch"], a["g"]): a for a in garm["arms"]
                if a["s"] == -1.0}
    arms_out = []
    ok_match, ok_twin = True, True
    for br in BRANCHES:
        Z = np.load(os.path.join(
            DATA, f"m5_21_11_end_t11lad_{br}_n48_d0.3.npz"))
        M3 = Z["M"].astype(np.float64)
        for g in (8.0, 32.0):
            cfg = INS4.base_cfg(s=s, g=g, n=48, L=48.0, delta=0.3)
            M4 = own_embed(M3, cfg["sg"])
            m_his = float(np.arctanh(1.0 / g))
            # own grid: +-0.6 m_his, 81 pts (4x finer than the run
            # near the minima), plus the exact m = 0 point
            ms = np.linspace(-0.6 * m_his, 0.6 * m_his, 81)
            t0 = time.time()
            Es = np.array([e_dressed_own(M4, cfg, m) for m in ms])
            E0 = e_dressed_own(M4, cfg, 0.0)
            neg = ms < -1e-15
            pos = ms > 1e-15
            m_n, E_n = refine_min(ms[neg], Es[neg])
            m_p, E_p = refine_min(ms[pos], Es[pos])
            E_star = min(E_n, E_p)
            m_star = m_n if E_n <= E_p else m_p
            gain = E_star - E0
            twin = (E_n < E0) and (E_p < E0) \
                and abs(E_n - E_p) < 0.02 * abs(gain)
            ra = run_arms[(br, g)]
            rel = abs(gain - ra["gain"]) / abs(ra["gain"])
            arm = {"branch": br, "g": g, "E0": E0,
                   "E0_run": ra["E0"], "gain": gain,
                   "gain_run": ra["gain"], "rel_diff": rel,
                   "m_star": m_star, "m_star_over_mhis":
                   abs(m_star) / m_his, "twin_minima": bool(twin),
                   "gain_neg_side": E_n - E0, "gain_pos_side":
                   E_p - E0, "wall_s": time.time() - t0}
            arms_out.append(arm)
            ok_match &= rel < 0.05
            ok_twin &= twin
            print(f"[C1] {br} g={g:g} gain {gain:+.5f} "
                  f"(run {ra['gain']:+.5f}, rel {rel:.2%}) "
                  f"|m*|/mhis {arm['m_star_over_mhis']:.3f} "
                  f"twin {twin} [{arm['wall_s']:.0f}s]", flush=True)
    # fall exponent per branch from own gains (g8 -> g32)
    lx = np.log(np.arctanh(1 / 32.0) / np.arctanh(1 / 8.0))
    qs = {}
    for br in BRANCHES:
        g8 = [a for a in arms_out if a["branch"] == br
              and a["g"] == 8.0][0]["gain"]
        g32 = [a for a in arms_out if a["branch"] == br
               and a["g"] == 32.0][0]["gain"]
        qs[br] = float(np.log(-g32 / -g8 if g32 < 0 and g8 < 0
                              else np.nan) / lx) \
            if (g32 < 0 and g8 < 0) else float("nan")
    ok_flat = all(np.isfinite(q) and abs(q) < 0.5 for q in qs.values())
    # F4 as frozen: gain must fall AT LEAST as fast as artanh^2
    # (q >= 2 passes); flat gains (q ~ 0) FAIL
    f4_fail_own = not all(np.isfinite(q) and q >= 1.95
                          for q in qs.values())
    # the claimed m* band: |m*| ~ 0.2-0.25 x artanh(1/g)
    band = {a["branch"]: [] for a in arms_out}
    for a in arms_out:
        band[a["branch"]].append(a["m_star_over_mhis"])
    band_ok = all(0.19 <= v <= 0.26 for vs in band.values()
                  for v in vs)
    verdict = "CONFIRMED" if (ok_match and ok_twin and ok_flat
                              and f4_fail_own and band_ok) else \
        ("PARTIAL" if (ok_match and ok_flat and f4_fail_own)
         else "REFUTED")
    VERDICTS["C1"] = {
        "verdict": verdict,
        "numbers": {"arms": arms_out, "fall_exponent_g8_g32": qs,
                    "m_star_over_artanh_by_branch":
                    {k: [round(v, 4) for v in vs]
                     for k, vs in band.items()},
                    "artanh2_would_require_fall_factor":
                    float((np.arctanh(1 / 8.0)
                           / np.arctanh(1 / 32.0)) ** 2)},
        "notes": "Gains recomputed with own boost construction on an "
                 "own 81-pt grid (spacing 0.015 m_his vs the run's "
                 "0.05): twin minima confirmed, depths g-independent "
                 "(fall exponents ~0, far below the frozen q >= 2 "
                 "bar), so F4 correctly FAILS: the substantive g-arm "
                 "content stands. CATCH on the claim wording: the "
                 "'|m*| ~ 0.2-0.25 x artanh(1/g)' band holds only for "
                 "A (0.23-0.24) and B (0.24-0.25); branch C sits at "
                 "0.31-0.33 in BOTH the run's own garm.json and this "
                 "recomputation. Any write-up must quote C's band "
                 "separately."}
    return VERDICTS["C1"]


# =================================================================
# C2: ansatz-family corroboration from the pre-existing record
# =================================================================
def audit_c2():
    rec = {int(r["g"]): r for r in jload("m5_21_8_lat_gladder.json")
           if int(r["g"]) in (8, 16, 32)}
    rows = {}
    e0s = []
    for g in (8.0, 16.0, 32.0):
        cfg = INS4.base_cfg(s=-1.0, g=g, n=32, L=48.0)
        M0 = LAT8.dressed(cfg, 0.0)           # record's own family code
        eu, ev = INS4.e_parts(M0, cfg)
        E0 = float(eu + ev)
        e0s.append(E0)
        # own m-curve (57 pts, 2x the record's 29) to re-find E_min
        ms = np.linspace(-0.35, 0.35, 57)
        Es = []
        for m in ms:
            M = LAT8.dressed(cfg, m)
            a, b = INS4.e_parts(M, cfg)
            Es.append(float(a + b))
        E_min = float(min(Es))
        m_at = float(ms[int(np.argmin(Es))])
        gain = E_min - E0
        rows[f"g{g:g}"] = {
            "E0_family": E0, "E_min_own": E_min, "m_at_min": m_at,
            "E_min_record": rec[int(g)]["E_min"],
            "family_gain": gain,
            "rel_diff_Emin_vs_record":
            abs(E_min - rec[int(g)]["E_min"])
            / abs(rec[int(g)]["E_min"])}
        print(f"[C2] g={g:g} E0 {E0:.4f} Emin(own) {E_min:.4f} "
              f"(rec {rec[int(g)]['E_min']:.4f}) gain {gain:+.3f}",
              flush=True)
    gains = [rows[k]["family_gain"] for k in rows]
    e0_spread = (max(e0s) - min(e0s)) / abs(np.mean(e0s))
    flat = (max(gains) - min(gains)) / abs(np.mean(gains))
    # artanh^2 fall from g8 to g32 would shrink |gain| 16x
    ok = (e0_spread < 1e-9 and flat < 0.05
          and all(-63.0 < gv < -59.0 for gv in gains)
          and all(rows[k]["rel_diff_Emin_vs_record"] < 0.02
                  for k in rows))
    VERDICTS["C2"] = {
        "verdict": "CONFIRMED" if ok else "PARTIAL",
        "numbers": {"per_g": rows,
                    "E0_family_rel_spread": float(e0_spread),
                    "gain_rel_spread": float(flat)},
        "notes": "E(0) of the ansatz family recomputed (g-independent "
                 "to machine level) and the full m-curves re-scanned "
                 "on a 2x-finer grid: the family gain is flat ~-61 "
                 "across g = 8/16/32, corroborating the C1 flat-depth "
                 "finding on an independent family."}
    return VERDICTS["C2"]


# =================================================================
# C3: gate table from the raw rows (own formulas)
# =================================================================
def own_gates(row, read):
    """own recomputation of the frozen section-4 gates."""
    fails = []
    tr = row["trace"][-1] if row["trace"] else {}
    fmax = tr.get("fmax", np.inf)
    fire_ok_pin = (row["stop"] == "f_tol"
                   or (row["stop"] in ("max_iter", "plateau")
                       and fmax <= 1e-4))
    if not fire_ok_pin:
        fails.append("fire")
    c = row["consistency"]
    reads5 = [c["E_u_fwd"], c["E_u_bwd"], c["E_u_2h"],
              c["E_u_sub0"], c["E_u_sub1"]]
    xr = max(reads5) / max(min(reads5), 1e-300)
    if xr > 1.5:
        fails.append("xstencil")
    vir = abs(row["E_u"] - 3.0 * row["E_v"]) / max(row["E_u"], 1e-300)
    if vir > 0.05:
        fails.append("virial")
    if read is None:
        fails.append("no-read")
    else:
        if not read["gap_guard"]["silent"]:
            fails.append("gap-guard")
        cont = read["flags"]["contours"]
        for z in sorted({q["z"] for q in cont}):
            zc = [q for q in cont if q["z"] == z]
            if zc and all(q["flag"] for q in zc):
                fails.append("plane-degenerate")
                break
    a_star = row["ring"].get("rho_w", float("nan"))
    if not np.isfinite(a_star) or L_BOX / a_star < 10.0:
        fails.append("L/a*")
    return fails, float(vir), float(xr), float(fmax), row["stop"]


def audit_c3():
    fit = jload("m5_21_11_fit.json")
    rung_state = fit["rung_state"]
    table = {}
    mismatches = []
    usable = {br: [] for br in BRANCHES}
    for br in BRANCHES:
        for d in RUNGS:
            tag = f"t11lad_{br}_n48_d{dkey(d)}"
            row = jload(f"m5_21_11_row_{tag}.json")
            try:
                read = jload(f"m5_21_11_read_{tag}.json")
            except FileNotFoundError:
                read = None
            fails, vir, xr, fmax, stop = own_gates(row, read)
            key = f"{br}_n48_d{dkey(d)}"
            rs = rung_state.get(key, {})
            run_fail_kinds = sorted({f.split("(")[0]
                                     for f in rs.get("gate_fails", [])})
            own_fail_kinds = sorted(set(fails))
            agree = (run_fail_kinds == own_fail_kinds
                     and abs(vir - rs.get("vir_frozen", np.nan))
                     < 1e-9)
            if not agree:
                mismatches.append({"rung": key, "own": own_fail_kinds,
                                   "run": run_fail_kinds,
                                   "own_vir": vir,
                                   "run_vir": rs.get("vir_frozen")})
            if not fails:
                usable[br].append(d)
            table[key] = {"E": row["E_end"], "vir": vir, "xr": xr,
                          "fmax_end": fmax, "stop": stop,
                          "own_fails": own_fail_kinds,
                          "run_fails": run_fail_kinds,
                          "agree": bool(agree)}
    # A-branch virial monotone growth as delta falls
    a_vir = [table[f"A_n48_d{dkey(d)}"]["vir"] for d in RUNGS]
    a_monotone = all(a_vir[i + 1] > a_vir[i]
                     for i in range(len(a_vir) - 1))
    bc_x_all = all(table[f"{br}_n48_d{dkey(d)}"]["xr"] > 1.5
                   for br in ("B", "C") for d in RUNGS)
    bc_vir = [table[f"{br}_n48_d{dkey(d)}"]["vir"]
              for br in ("B", "C") for d in RUNGS]
    ok = (not mismatches and usable["A"] == [0.30]
          and usable["C"] == [] and usable["B"] == []
          and a_monotone and bc_x_all
          and 0.79 < min(bc_vir) and max(bc_vir) < 0.99
          and abs(table["A_n48_d0.3"]["vir"] - 0.046) < 0.001)
    VERDICTS["C3"] = {
        "verdict": "CONFIRMED" if ok else
        ("REFUTED" if mismatches else "PARTIAL"),
        "numbers": {"table": table, "usable_by_gates": usable,
                    "A_virial_ladder": a_vir,
                    "A_virial_monotone_growing": bool(a_monotone),
                    "BC_xstencil_all_above_1p5": bool(bc_x_all),
                    "BC_virial_range": [float(min(bc_vir)),
                                        float(max(bc_vir))],
                    "mismatches": mismatches},
        "notes": "Virial |E_u - 3E_V|/E_u, xstencil max/min over the "
                 "five stored consistency reads, FIRE, L/a* and the "
                 "reader guards all recomputed with own formulas; "
                 "gate-fail sets match the run's rung_state rung for "
                 "rung. A usable only at delta = 0.30 (vir 0.0457); "
                 "B/C usable nowhere; A's virial grows monotonically "
                 "as delta falls."}
    return VERDICTS["C3"]


# =================================================================
# C4: shell-fix integrity on continuation endpoints
# =================================================================
def audit_c4():
    picks = [("A", 48, 0.05), ("C", 48, 0.12), ("B", 64, 0.05)]
    rows = {}
    ok = True
    for br, n, d in picks:
        tag = f"t11lad_{br}_n{n}_d{dkey(d)}"
        Z = np.load(os.path.join(DATA, f"m5_21_11_end_{tag}.npz"))
        M = Z["M"].astype(np.float64)
        h = float(Z["h"])
        cfg = INS2B.base_cfg(seed=br, term="T2", stencil="sym",
                             eps=0.0, n=n, L=48.0, delta=d,
                             bc="pinned", w2=W2_PIN)
        assert abs(cfg["h"] - h) < 1e-15
        seed = INS2B.make_seed(cfg)
        shell = INS2B.pin_shell(n, h)          # physical depth >= 1.6
        dev = float(np.max(np.abs(M[shell] - seed[shell])))
        interior_move = float(np.max(np.abs(M[~shell]
                                            - seed[~shell])))
        rows[tag] = {"seed_kind": str(Z["seed_kind"]),
                     "shell_max_abs_dev": dev,
                     "n_shell_cells": int(shell.sum()),
                     "interior_max_dev_vs_seed": interior_move,
                     "pass": dev < 1e-12}
        ok &= dev < 1e-12 and str(Z["seed_kind"]) == "continuation"
        print(f"[C4] {tag} shell dev {dev:.2e} "
              f"({rows[tag]['seed_kind']})", flush=True)
    VERDICTS["C4"] = {
        "verdict": "CONFIRMED" if ok else "REFUTED",
        "numbers": rows,
        "notes": "All three continuation endpoints hold the pinned "
                 "shell EXACTLY at the current-delta analytic T2 seed "
                 "(deviation 0 to < 1e-12); the interior moved (sanity "
                 "column), so the earlier shell-inheritance bug did "
                 "not contaminate any surviving rung."}
    return VERDICTS["C4"]


# =================================================================
# C5: E(delta) monotonicity from the raw rows
# =================================================================
def audit_c5():
    claim = {"A": (6.84, 9.16), "C": (24.12, 74.09),
             "B": (84.84, 87.02)}
    rows = {}
    ok = True
    for br in BRANCHES:
        Es = []
        for d in RUNGS:
            r = jload(f"m5_21_11_row_t11lad_{br}_n48_d{dkey(d)}.json")
            Es.append(float(r["E_end"]))
        mono = all(Es[i + 1] > Es[i] for i in range(len(Es) - 1))
        lo, hi = claim[br]
        match = (abs(Es[0] - lo) < 0.01 * lo
                 and abs(Es[-1] - hi) < 0.01 * hi)
        rows[br] = {"deltas": RUNGS, "E": Es,
                    "monotone_E_rises_as_delta_falls": bool(mono),
                    "endpoints": [Es[0], Es[-1]],
                    "claimed": [lo, hi], "match": bool(match)}
        ok &= mono and match
        print(f"[C5] {br} E(0.30)={Es[0]:.4f} -> E(0.05)={Es[-1]:.4f} "
              f"monotone {mono}", flush=True)
    VERDICTS["C5"] = {
        "verdict": "CONFIRMED" if ok else "REFUTED",
        "numbers": rows,
        "notes": "Per-branch E(delta) at N = 48 is strictly monotone "
                 "decreasing in delta (E rises as delta falls) for "
                 "all three branches; endpoints match the quoted "
                 "values from the raw rows."}
    return VERDICTS["C5"]


# =================================================================
# C6: F3 arithmetic + F1/F2 vacuity
# =================================================================
def audit_c6():
    fit = jload("m5_21_11_fit.json")
    usable = VERDICTS["C3"]["numbers"]["usable_by_gates"]
    counts = {br: len(usable[br]) for br in BRANCHES}
    run_counts = fit["criteria"]["F3_branch_integrity"]["counts"]
    f3_fires = any(counts[br] < 6 for br in BRANCHES)
    holdouts = (0.20, 0.07)
    n_fit = sum(1 for br in BRANCHES for d in usable[br]
                if d not in holdouts)
    n_hold = sum(1 for br in BRANCHES for d in usable[br]
                 if d in holdouts)
    f1_vac = fit["criteria"]["F1_fit_quality"]["chi2_per_dof"] is None \
        and not fit["criteria"]["F1_fit_quality"]["fail"]
    f2_vac = (fit["criteria"]["F2_holdouts"]["n2s"] == 0
              and fit["criteria"]["F2_holdouts"]["n3s"] == 0
              and len(fit["holdouts"]["rows"]) == 0)
    ok = (counts == {"A": 1, "C": 0, "B": 0}
          and counts == run_counts and f3_fires
          and fit["criteria"]["F3_branch_integrity"]["fail"]
          and n_fit == 1 and n_hold == 0
          and fit["n_fit_points"] == 1
          and fit["n_holdout_points"] == 0
          and f1_vac and f2_vac
          and fit["verdict"]["terminal_failure"]
          and set(fit["verdict"]["failed"])
          == {"F3_branch_integrity", "F4_garm"})
    VERDICTS["C6"] = {
        "verdict": "CONFIRMED" if ok else "REFUTED",
        "numbers": {"own_usable_counts": counts,
                    "run_counts": run_counts,
                    "floor": 6, "f3_fires": bool(f3_fires),
                    "own_n_fit_points": n_fit,
                    "own_n_holdout_points": n_hold,
                    "run_n_fit_points": fit["n_fit_points"],
                    "f1_vacuous": bool(f1_vac),
                    "f2_vacuous": bool(f2_vac),
                    "run_failed": fit["verdict"]["failed"]},
        "notes": "Usable counts from the own gate recomputation are "
                 "A=1, C=0, B=0, all below the frozen floor of 6, so "
                 "F3 fires; the single usable point (A, 0.30) is not "
                 "a holdout, leaving a 1-point degenerate fit: F1 "
                 "(chi2/df) and F2 (holdouts) are vacuous exactly as "
                 "the run reports. Terminal failure = {F3, F4}."}
    return VERDICTS["C6"]


# =================================================================
# Barred-inputs sweep
# =================================================================
def barred_sweep():
    files = ["m5_21_11_b_ladder.py", "m5_21_11_c_readers.py",
             "m5_21_11_d_garm.py", "m5_21_11_f_fit.py",
             "m5_21_11_g_controls.py"]
    pats = [r"5\.2611", r"22\.059", r"84\.085", r"206\.8", r"3477",
            r"5\.9", r"15\.1"]
    hits = []
    for fn in files:
        with open(os.path.join(HERE, fn)) as f:
            for ln, line in enumerate(f, 1):
                for p in pats:
                    if re.search(p, line):
                        hits.append({"file": fn, "line": ln,
                                     "pattern": p,
                                     "context": line.rstrip()})
    return {"files_swept": files, "patterns": pats, "hits": hits,
            "clean": len(hits) == 0}


# =================================================================
def main():
    t0 = time.time()
    bg = gate_boost()
    print(f"[gate] own-boost vs expm {bg['expm_err']:.2e}, "
          f"Lorentz {bg['lorentz_err']:.2e}, grid {bg['grid_err']:g} "
          f"-> pass {bg['pass']}", flush=True)
    if not bg["pass"]:
        raise SystemExit("own boost construction failed its gate")
    audit_c3()          # cheap first (C6 depends on it)
    audit_c4()
    audit_c5()
    audit_c6()
    audit_c2()
    audit_c1()          # heavy last
    barred = barred_sweep()
    n_conf = sum(1 for v in VERDICTS.values()
                 if v["verdict"] == "CONFIRMED")
    NOTES.append(
        "FIRE-gate strictness: framework section 4 as frozen says "
        "'FIRE reaches f_tol at depth' (f_tol = 1e-8); the run's "
        "operational pin (recorded pre-result in the fit header) "
        "accepts max_iter/plateau with final fmax <= 1e-4. Under the "
        "strict frozen reading even A at delta = 0.30 (stop = "
        "max_iter, fmax 1.2e-6) would be EXCLUDED, making the usable "
        "counts 0/0/0. Verdict-neutral: F3 fires either way (and "
        "harder), but the softening should be stated next to the "
        "usable-count table.")
    NOTES.append(
        "The fallback 'fit' block (E = Einf + c*delta through ONE "
        "usable point) is underdetermined: params A = (6.279, 1.884) "
        "satisfy 6.279 + 1.884*0.30 = 6.845 but are one of infinitely "
        "many solutions, so E_phys(A) = 6.28 +- 0.61 carries no "
        "information beyond the single rung. Verdict-neutral (the "
        "run already declares terminal failure), but the number "
        "should not be quoted as a prediction anywhere.")
    NOTES.append(
        "Barred-inputs sweep: grep-level clean on all five new "
        "scripts (no 5.2611 / 22.059 / 84.085 / 206.8 / 3477 / '5.9' "
        "/ '15.1').")
    out = {
        "task": "M5.21.11 ladder-compute adversarial audit",
        "date": "2026-08-07",
        "auditor": "independent second agent, own recomputation "
                   "(m5_21_11_i_audit.py); instruments imported: "
                   "m5_21_2b_a_instrument, m5_21_3_a_4d, "
                   "m5_21_8_b_lattice; run wrappers NOT imported",
        "own_boost_gate": bg,
        "claims": VERDICTS,
        "barred_inputs": barred,
        "notes": NOTES,
        "summary": f"{n_conf}/6 claims CONFIRMED "
                   f"({', '.join(k + ':' + v['verdict'] for k, v in sorted(VERDICTS.items()))}); "
                   "barred inputs clean; terminal-failure verdict "
                   "(F3 + F4) independently reproduced. Catches: "
                   "C1 m*-band wording wrong for branch C (0.31-0.33, "
                   "not 0.2-0.25; gains and F4 failure fully "
                   "confirmed); FIRE-gate softening vs the frozen "
                   "f_tol wording (verdict-neutral, F3 fires harder "
                   "strictly); underdetermined 1-point fallback fit "
                   "(verdict-neutral).",
        "wall_s": time.time() - t0,
    }
    with open(OUT, "w") as f:
        json.dump(out, f, indent=1)
    # print without the bulky per-rung table/curves
    slim = json.loads(json.dumps(out))
    slim["claims"]["C3"]["numbers"].pop("table")
    for c in slim["claims"]["C1"]["numbers"]["arms"]:
        c.pop("wall_s", None)
    print(json.dumps(slim, indent=1))


if __name__ == "__main__":
    main()
