"""M5.21.9 pre-gate verify: the author's negative-delta suggestion
(reply 2026-07-19 14:53, tasks/m5_21_convo.md), checked before anything
runs on it (the M5.21.8 verify-first pattern).

The author's suggestion: the omega^2-coefficient region that Minimize
routes to -inf ("0 < delta <= 2 && g > 1") "could be for
(g > 1 && delta <= 0), so maybe delta should be negative? (not an
issue, especially when shifting spectrum)".

CHECKS (on the ALREADY-AUDITED M5.21.8 instruments, imported):
    N1 exact factorization of the transcribed omega^2 coefficient:
        coeff(g, d) = 8 d^3 [(g-1) d - 2 g] / (g-1)
        => for g > 1: coeff > 0 iff d < 0 or d > 2g/(g-1).
        Verified numerically against the raw transcription.
    N2 sign map extended to delta < 0 (the author's new region):
        expect POSITIVE everywhere at g > 1, delta < 0.
    N3 independent-pipeline confirmation (ansatz-level, no formula):
        dE/d(omega^2) at probe points with delta < 0, m = -|m*|.
    N4 boundedness is NOT a clock: Hm (Out[54]) is LINEAR in omega^2,
        so a positive coefficient parks the minimum at omega* = 0
        exactly (= the Out[51] static value). Checked on the formula
        (fine omega grid) and on the pipeline E(omega) curve.
    N5 the physical corner: sign + value at (g, d) = (1e10, -1e-10)
        vs (1e10, +1e-10), transcribed AND factored forms.
    N6 the m-sector at delta < 0: does the m* minimum survive?
        m-scan at (8, -1e-4) (the delta -> 0 limit) and (8, -0.3).
    N7 the "shifting spectrum" remark is EXACT: M -> M + c I under
        d -> d + c I, so every derivative and F is invariant;
        negative delta is a genuine reordering, not a shift artifact.

Out: data/m5_21_9_negdelta.json + printed tables.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
sys.path.insert(0, HERE)

from m5_21_8_a_verify import (  # noqa: E402
    E_cone, Hm_formula, boost_hedge, fin_formula, mstar, om2_coeff,
    rot, G1, G2, G3, XI,
)


def coeff_factored(g, dl):
    """N1: the closed factorization of om2_coeff."""
    return 8.0 * dl ** 3 * ((g - 1.0) * dl - 2.0 * g) / (g - 1.0)


def n1():
    rng = np.random.default_rng(219)
    worst = 0.0
    for _ in range(4000):
        g = float(np.exp(rng.uniform(np.log(1.01), np.log(1e10))))
        dl = float(rng.uniform(-3.0, 3.0))
        if abs(dl) < 1e-6:
            continue
        a, b = om2_coeff(g, dl), coeff_factored(g, dl)
        rel = abs(a - b) / max(abs(a), abs(b), 1e-300)
        worst = max(worst, rel)
    row = {"samples": 4000, "worst_rel": float(worst),
           "factored_ok": bool(worst < 1e-10),
           "sign_law": "g>1: coeff>0 iff d<0 or d>2g/(g-1)"}
    print("N1", row, flush=True)
    return row


def n2():
    gs = np.array([1.5, 3.0, 8.0, 30.0, 1e3, 1e10])
    dls = np.array([-3.0, -1.9, -1.0, -0.3, -0.05, -1e-3, -1e-10,
                    1e-10, 0.3])
    grid = [[float(np.sign(om2_coeff(g, d))) for d in dls]
            for g in gs]
    print("N2 om2-coeff sign map incl. delta<0 (rows g, cols delta):")
    for g, row in zip(gs, grid):
        print(f"  g={g:8g}: {row}")
    neg_cols = [j for j, d in enumerate(dls) if d < 0]
    all_pos = all(grid[i][j] == 1.0 for i in range(len(gs))
                  for j in neg_cols)
    print("N2 all delta<0 cells positive:", all_pos, flush=True)
    return {"gs": gs.tolist(), "dls": dls.tolist(), "sign": grid,
            "all_neg_delta_positive": bool(all_pos)}


def n3():
    rows = []
    for (g, dl) in ((8.0, -0.3), (8.0, -0.05), (3.0, -0.3),
                    (100.0, -1e-3)):
        m = -abs(mstar(g))
        r = 1.3
        om1, om2 = 0.2, 0.4
        dE = (E_cone(g, dl, m, om2, r=r) - E_cone(g, dl, m, om1, r=r)) \
            / (om2 ** 2 - om1 ** 2)
        rows.append({"g": g, "dl": dl, "m_used": float(m),
                     "dE_dom2_pipeline": float(dE),
                     "om2_coeff_formula": float(om2_coeff(g, dl)),
                     "same_sign": bool(np.sign(dE)
                                       == np.sign(om2_coeff(g, dl))),
                     "positive": bool(dE > 0)})
        print("N3", rows[-1], flush=True)
    return rows


def n4():
    g, dl = 8.0, -0.3
    # formula level: linear in omega^2, min at 0, = Out[51] static value
    oms = np.linspace(0.0, 2.0, 161)
    vals = np.array([Hm_formula(g, dl, 1.0, om) for om in oms])
    co = np.polyfit(oms ** 2, vals, 1)
    resid = float(np.max(np.abs(np.polyval(co, oms ** 2) - vals))
                  / np.max(np.abs(vals)))
    row_f = {"g": g, "dl": dl,
             "linear_in_om2_resid": resid,
             "argmin_om": float(oms[int(np.argmin(vals))]),
             "hm_at_0": float(Hm_formula(g, dl, 1.0, 0.0)),
             "his_fin_formula": float(fin_formula(g, dl, 1.0)),
             "slope_positive": bool(co[0] > 0)}
    print("N4 formula", row_f, flush=True)
    # pipeline level: E(omega) curve at m = -|m*|. The fit uses the
    # omega > 0 points only: the instrument's omega = 0 branch is a
    # t = 0 SNAPSHOT (Ha_avg special case) while omega > 0 points are
    # clock-phase-averaged, so E(0) carries a known evaluator offset.
    m = -abs(mstar(g))
    oms_p = np.linspace(0.0, 0.6, 13)
    Es = np.array([E_cone(g, dl, m, om, r=1.3) for om in oms_p])
    pos = oms_p > 0
    co_p = np.polyfit(oms_p[pos] ** 2, Es[pos], 1)
    resid_p = float(np.max(np.abs(np.polyval(co_p, oms_p[pos] ** 2)
                                  - Es[pos]))
                    / max(np.max(np.abs(Es[pos])), 1e-300))
    row_p = {"oms": oms_p.tolist(), "E": Es.tolist(),
             "linear_in_om2_resid_pos": resid_p,
             "argmin_om_pos": float(
                 oms_p[pos][int(np.argmin(Es[pos]))]),
             "slope": float(co_p[0]),
             "slope_positive": bool(co_p[0] > 0),
             "intercept_fit": float(co_p[1]),
             "E_at_0_snapshot": float(Es[0]),
             "snapshot_offset": float(Es[0] - co_p[1])}
    print("N4 pipeline", {k: v for k, v in row_p.items()
                          if k not in ("oms", "E")}, flush=True)
    return {"formula": row_f, "pipeline": row_p}


def n5():
    rows = []
    for dl in (-1e-10, 1e-10):
        g = 1e10
        rows.append({"g": g, "dl": dl,
                     "om2_coeff_transcribed": float(om2_coeff(g, dl)),
                     "om2_coeff_factored": float(coeff_factored(g, dl)),
                     "sign": float(np.sign(om2_coeff(g, dl)))})
        print("N5", rows[-1], flush=True)
    return rows


def n6():
    rows = []
    for dl in (-1e-4, -0.3):
        ms = np.linspace(-0.5, 0.5, 41)
        Es = [E_cone(8.0, dl, m, 0.3) for m in ms]
        i = int(np.argmin(Es))
        m_fit = ms[i]
        if 0 < i < len(ms) - 1:
            a, b, c = Es[i - 1], Es[i], Es[i + 1]
            m_fit = ms[i] - 0.5 * (ms[1] - ms[0]) * (c - a) / \
                (c - 2 * b + a)
        rows.append({"g": 8.0, "dl": dl,
                     "m_min_num": float(m_fit),
                     "m_star_abs": float(abs(mstar(8.0))),
                     "interior_min": bool(0 < i < len(ms) - 1)})
        print("N6", rows[-1], flush=True)
    return rows


def _M_spec(t, x, y, z, spec, m, om):
    """the ansatz with a GENERAL diagonal spectrum (for N7 only)."""
    phi = np.arctan2(y, x)
    th = -np.arctan2(z, np.sqrt(x * x + y * y))
    Qh = rot(G3, phi) @ rot(G2, th) @ rot(G1, om * t)
    Qb = boost_hedge(m, x, y, z)
    Q = Qb @ Qh
    return Q @ np.diag(spec) @ Q.T


def _H_spec(t, x, y, z, spec, m, om, h=1e-5):
    dM = []
    for k, v in enumerate((t, x, y, z)):
        args = [t, x, y, z]
        a_p = list(args); a_p[k] = v + h
        a_m = list(args); a_m[k] = v - h
        dM.append((_M_spec(*a_p, spec, m, om)
                   - _M_spec(*a_m, spec, m, om)) / (2 * h))
    H = 0.0
    for i in range(4):
        for j in range(4):
            F = dM[i] @ XI @ dM[j] - dM[j] @ XI @ dM[i]
            H += (F[1, 2] ** 2 + F[1, 3] ** 2 + F[2, 3] ** 2
                  - F[0, 1] ** 2 - F[0, 2] ** 2 - F[0, 3] ** 2)
    return H


def n7():
    """shift invariance: EXACT on the undressed (m = 0, orthogonal Q)
    family, and expected BROKEN on the boosted one: Qb is symmetric,
    not orthogonal, so d -> d + c I adds c Qb(x)^2, position-dependent.
    """
    rng = np.random.default_rng(220)
    rows = []
    for m_mode in (0.0, None):
        for _ in range(4):
            g, dl = 8.0, -0.3
            m = 0.0 if m_mode == 0.0 else rng.uniform(-0.3, 0.3)
            om = rng.uniform(0.05, 0.6)
            t = rng.uniform(0, 5.0)
            x, z = rng.uniform(0.8, 3.0), rng.uniform(0.2, 2.0)
            c = rng.uniform(-2.0, 2.0)
            base = [g, 1.0, dl, 0.0]
            H0 = _H_spec(t, x, 1e-12, z, base, m, om)
            Hc = _H_spec(t, x, 1e-12, z, [s + c for s in base], m, om)
            rel = abs(H0 - Hc) / max(abs(H0), 1e-300)
            rows.append({"m": float(m), "c": float(c),
                         "H0": float(H0),
                         "rel_shift_change": float(rel)})
            print("N7", rows[-1], flush=True)
    w_rot = max(r["rel_shift_change"] for r in rows if r["m"] == 0.0)
    w_boo = max(r["rel_shift_change"] for r in rows if r["m"] != 0.0)
    rows.append({"worst_rel_undressed_m0": float(w_rot),
                 "worst_rel_boosted": float(w_boo),
                 "invariant_undressed": bool(w_rot < 1e-8),
                 "broken_by_boost": bool(w_boo > 1e-3)})
    print("N7", rows[-1], flush=True)
    return rows


def n8():
    """does negative delta also heal the dropped cone channel? the
    M5.21.8 C4 read at (8, +0.3): E(Delta) FALLS as the cone shrinks
    (negatively divergent). Same ladder at (8, -0.3)."""
    g, dl = 8.0, -0.3
    m = -abs(mstar(g))
    Ds = (0.4, 0.2, 0.1, 0.05)
    Es = [E_cone(g, dl, m, 0.3, Delta=D) for D in Ds]
    d43 = (Es[3] - Es[2]) / (np.log(Ds[2]) - np.log(Ds[3]))
    row = {"g": g, "dl": dl,
           "E_at_Delta": dict(zip(map(str, Ds), map(float, Es))),
           "log_slope_last": float(d43),
           "still_negatively_divergent": bool(d43 < 0)}
    print("N8", row, flush=True)
    return row


if __name__ == "__main__":
    out = {"N1_factorization": n1(),
           "N2_sign_map_negdelta": n2(),
           "N3_pipeline_probes": n3(),
           "N4_linearity_omstar0": n4(),
           "N5_physical_corner": n5(),
           "N6_mstar_negdelta": n6(),
           "N7_shift_invariance": n7(),
           "N8_cone_channel": n8()}
    os.makedirs(DATA, exist_ok=True)
    with open(os.path.join(DATA, "m5_21_9_negdelta.json"), "w") as f:
        json.dump(out, f, indent=1)
    print("saved m5_21_9_negdelta.json")
