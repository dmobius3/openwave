"""Independent adversarial audit of the m5_21_9 round-2 Larmor claims.

Refits saved JSON data only. No simulations. Writes one JSON verdict file.
All fits and unwraps here are the auditor's own, not the primary script's.
"""

import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
OUT = DATA / "m5_21_9_audit_r2.json"


def load(name):
    with open(DATA / name) as f:
        return json.load(f)


def finite_rows(rows):
    keep = []
    for r in rows:
        vals = [r["phi_xy"], r["E"], r["KE"],
                r["J"]["x"], r["J"]["y"], r["J"]["z"]]
        if all(isinstance(v, (int, float)) and math.isfinite(v)
               for v in vals):
            keep.append(r)
        else:
            break
    return keep


def arrs(rows):
    t = np.array([r["t"] for r in rows])
    phi = np.unwrap(np.array([r["phi_xy"] for r in rows]))
    e = np.array([r["E"] for r in rows])
    jx = np.array([r["J"]["x"] for r in rows])
    jy = np.array([r["J"]["y"] for r in rows])
    jz = np.array([r["J"]["z"] for r in rows])
    return t, phi, e, jx, jy, jz


def lin_slope(t, y):
    if len(t) < 2:
        return float("nan")
    return float(np.polyfit(t, y, 1)[0])


def quad_fit_origin(t, y):
    """Fit y = a*t + b*t^2 through the origin. Return a, b, se_a, se_b."""
    x = np.column_stack([t, t * t])
    coef, res, _, _ = np.linalg.lstsq(x, y, rcond=None)
    n, p = len(t), 2
    yhat = x @ coef
    rss = float(np.sum((y - yhat) ** 2))
    dof = max(n - p, 1)
    s2 = rss / dof
    cov = s2 * np.linalg.inv(x.T @ x)
    se = np.sqrt(np.diag(cov))
    ss_tot = float(np.sum(y ** 2))
    r2 = 1.0 - rss / ss_tot if ss_tot > 0 else float("nan")
    return coef[0], coef[1], se[0], se[1], rss, r2


def lin_fit_origin(t, y):
    """Fit y = a*t through the origin. Return a, rss."""
    a = float(np.sum(t * y) / np.sum(t * t))
    rss = float(np.sum((y - a * t) ** 2))
    return a, rss


def pair_dphi(fwd, rev):
    """Common finite prefix, offset-removed pair phase difference."""
    rf = finite_rows(fwd["rows"])
    rr = finite_rows(rev["rows"])
    n = min(len(rf), len(rr))
    tf, pf, _, _, _, _ = arrs(rf[:n])
    tr, pr, _, _, _, _ = arrs(rr[:n])
    assert np.allclose(tf, tr)
    dphi = pf - pr
    dphi = dphi - dphi[0]
    return tf, dphi


def split_half(t, y):
    n = len(t)
    h = n // 2
    s1 = lin_slope(t[:h], y[:h])
    s2 = lin_slope(t[h:], y[h:])
    ratio = s2 / s1 if s1 != 0 else float("inf")
    return s1, s2, ratio


report = {}

# ---------------- load everything ----------------
ctrl = load("m5_21_9_larmor_r2ctrl.json")
p02f = load("m5_21_9_larmor_r2p02f.json")
p02r = load("m5_21_9_larmor_r2p02r.json")
p03f = load("m5_21_9_larmor_r2p03f.json")
p03r = load("m5_21_9_larmor_r2p03r.json")
env03 = load("m5_21_9_larmor_env0.03.json")
env03rev = load("m5_21_9_larmor_env0.03rev.json")
env03jm = load("m5_21_9_larmor_env0.03_jm.json")

# ---------------- R1: lifetimes + control health ----------------
r1 = {"runs": {}}
for name, d in [("r2p02f", p02f), ("r2p02r", p02r),
                ("r2p03f", p03f), ("r2p03r", p03r),
                ("r2ctrl", ctrl)]:
    fr = finite_rows(d["rows"])
    r1["runs"][name] = {
        "n_rows": len(d["rows"]),
        "n_finite": len(fr),
        "last_finite_t": fr[-1]["t"] if fr else None,
        "E_first": fr[0]["E"] if fr else None,
        "E_last": fr[-1]["E"] if fr else None,
    }

fc = finite_rows(ctrl["rows"])
tc, pc, ec, jxc, jyc, jzc = arrs(fc)
de = np.diff(ec)
r1["ctrl_E_monotone_decline"] = bool(np.all(de < 0))
r1["ctrl_n_positive_dE_snaps"] = int(np.sum(de > 0))
r1["ctrl_max_positive_dE"] = float(np.max(de))
r1["ctrl_max_abs_dE_per_snap"] = float(np.max(np.abs(de)))
r1["ctrl_E_range"] = [float(np.min(ec)), float(np.max(ec))]
r1["ctrl_max_abs_dJ_per_snap"] = float(max(
    np.max(np.abs(np.diff(jxc))),
    np.max(np.abs(np.diff(jyc))),
    np.max(np.abs(np.diff(jzc)))))
r1["ctrl_phi_drift_total"] = float(pc[-1] - pc[0])
r1["ctrl_omega_fit"] = lin_slope(tc, pc)

# claim-faithful checks (the claim says trend decline, not monotone)
ok_field = all(6 <= r1["runs"][k]["n_finite"] <= 7
               and 10.0 <= r1["runs"][k]["last_finite_t"] <= 12.0
               for k in ("r2p02f", "r2p02r", "r2p03f", "r2p03r"))
ok_ctrl = (r1["runs"]["r2ctrl"]["n_finite"] == 41
           and abs(r1["runs"]["r2ctrl"]["E_first"] - 6.74) < 0.01
           and abs(r1["runs"]["r2ctrl"]["E_last"] - 5.54) < 0.01
           and r1["ctrl_E_range"][0] > 5.4
           and r1["ctrl_E_range"][1] < 6.8
           and r1["ctrl_max_abs_dJ_per_snap"] < 1e-3)
r1["nuance"] = (
    "Confirmed as stated, two sharpenings. (1) The control E decline "
    "is a trend with small oscillation, not monotone: 16 of 40 snap "
    "steps tick up, largest +0.116; E stays inside [5.54, 6.74] and "
    "J steps stay below 8.5e-4, so no blowup. (2) Two field runs are "
    "already corrupted before the NaN: r2p02f last finite row has "
    "E = -246.4 (KE 258.5) and r2p03r has E = -26.9 (KE 54.2), so "
    "the usable field-on prefix is t <= 8..10, one snap shorter than "
    "the finite prefix. Field-run E swings wildly from t = 2 on "
    "(e.g. 53.4 -> 24.7 -> 9.0 -> 44.0), so the arena is disturbed "
    "from the first field-on snap.")
r1["verdict"] = "CONFIRMED" if (ok_field and ok_ctrl) else "REFUTED"
report["R1"] = r1

# ---------------- R2: pair dphi accelerates ----------------
r2 = {}
for tag, fwd, rev in [("p02", p02f, p02r), ("p03", p03f, p03r)]:
    t, dphi = pair_dphi(fwd, rev)
    s1, s2, ratio = split_half(t, dphi)
    r2[tag] = {
        "n_common": len(t),
        "t_last": float(t[-1]),
        "dphi_last": float(dphi[-1]),
        "half1_slope": s1,
        "half2_slope": s2,
        "split_half_ratio": ratio,
    }
accel_02 = r2["p02"]["split_half_ratio"]
accel_03 = r2["p03"]["split_half_ratio"]
r2["claimed_ratios"] = {"p02": 34.0, "p03": 40.0}
match_02 = 20.0 < accel_02 < 60.0
match_03 = 20.0 < accel_03 < 60.0
r2["verdict"] = ("CONFIRMED" if (match_02 and match_03)
                 else ("NUANCE" if accel_02 > 5 and accel_03 > 5
                       else "REFUTED"))
report["R2"] = r2

# ---------------- R3: J-flip run is an exact symmetry image ----------
r3 = {}
tp, pp, ep, jxp, jyp, jzp = arrs(env03["rows"])
tm, pm, em, jxm, jym, jzm = arrs(env03jm["rows"])
r3["omega_plusJ"] = lin_slope(tp, pp)
r3["omega_minusJ"] = lin_slope(tm, pm)
r3["omega_rel_diff"] = float(
    abs(r3["omega_plusJ"] - r3["omega_minusJ"])
    / abs(r3["omega_plusJ"]))
r3["max_abs_dE_rowwise"] = float(np.max(np.abs(ep - em)))
r3["max_abs_phi_image_err"] = float(
    np.max(np.abs((pp - np.pi) - pm)))
r3["max_abs_J_image_err"] = float(max(
    np.max(np.abs(jxp + jxm)),
    np.max(np.abs(jyp + jym)),
    np.max(np.abs(jzp + jzm))))
same_omega = r3["omega_rel_diff"] < 1e-6
same_e = r3["max_abs_dE_rowwise"] < 5e-4
r3["nuance"] = (
    "Stronger than claimed: not merely 3-decimal E agreement but a "
    "machine-precision image. Row by row, E identical to 1.8e-15, "
    "phi_jm = phi_plus - pi to 3.3e-16, J_jm = -J_plus to 3.5e-18. "
    "The J-flip run is the same trajectory under an exact symmetry "
    "map and carries zero independent information about a J-driven "
    "precession sign.")
r3["verdict"] = "CONFIRMED" if (same_omega and same_e) else "REFUTED"
report["R3"] = r3

# ---------------- R4: attack the final verdict ----------------
r4 = {"pairs": {}}
b02 = p02f["B_meas"]
b03 = p03f["B_meas"]
r4["B_ratio_p03_over_p02"] = float(b03 / b02)
for tag, fwd, rev in [("p02", p02f, p02r), ("p03", p03f, p03r)]:
    t, dphi = pair_dphi(fwd, rev)
    a, b, se_a, se_b, rss_q, r2q = quad_fit_origin(t, dphi)
    a_lin, rss_l = lin_fit_origin(t, dphi)
    n = len(t)
    f_stat = ((rss_l - rss_q) / 1.0) / (rss_q / (n - 2))
    r4["pairs"][tag] = {
        "quad_a": float(a),
        "quad_a_se": float(se_a),
        "quad_a_tstat": float(a / se_a) if se_a > 0 else None,
        "quad_b": float(b),
        "quad_b_se": float(se_b),
        "quad_b_tstat": float(b / se_b) if se_b > 0 else None,
        "quad_rss": rss_q,
        "quad_r2": r2q,
        "lin_only_a": a_lin,
        "lin_only_rss": rss_l,
        "rss_ratio_lin_over_quad": (rss_l / rss_q
                                    if rss_q > 0 else float("inf")),
        "f_stat_quad_vs_lin": float(f_stat),
        "f_crit_0p05_1_dof": 7.71 if n == 6 else None,
    }
a02 = r4["pairs"]["p02"]["quad_a"]
a03 = r4["pairs"]["p03"]["quad_a"]
b02q = r4["pairs"]["p02"]["quad_b"]
b03q = r4["pairs"]["p03"]["quad_b"]
r4["linear_term_ratio_p03_over_p02"] = float(a03 / a02) if a02 else None
r4["quad_term_ratio_p03_over_p02"] = float(b03q / b02q) if b02q else None
r4["dphi_last_ratio_p03_over_p02"] = float(
    r2["p03"]["dphi_last"] / r2["p02"]["dphi_last"])

# round-1 pair (dt = 0.005, 21 rows): same test at finer sampling
t1, dphi1 = pair_dphi(env03, env03rev)
s1a, s1b, sr1 = split_half(t1, dphi1)
a1, b1, se_a1, se_b1, rss_q1, r2q1 = quad_fit_origin(t1, dphi1)
r4["round1_pair_env0.03"] = {
    "n_common": len(t1),
    "dphi_last": float(dphi1[-1]),
    "split_half_ratio": sr1,
    "quad_a": float(a1),
    "quad_a_se": float(se_a1),
    "quad_a_tstat": float(a1 / se_a1) if se_a1 > 0 else None,
    "quad_b": float(b1),
    "quad_b_tstat": float(b1 / se_b1) if se_b1 > 0 else None,
}

# fieldless control drift vs the round-1 single-run Omega claim
r4["ctrl_omega_fieldless"] = r1["ctrl_omega_fit"]
r4["round1_omega_claimed"] = 8.249e-4
r4["ctrl_drift_over_round1_omega"] = float(
    r1["ctrl_omega_fit"] / 8.249e-4)

# a resolved linear rate needs: positive a, significant, B-scaled 2.25x
a_pos = a02 > 0 and a03 > 0
lin_sig_02 = abs(r4["pairs"]["p02"]["quad_a_tstat"] or 0) > 2.0
lin_sig_03 = abs(r4["pairs"]["p03"]["quad_a_tstat"] or 0) > 2.0
a_ratio = r4["linear_term_ratio_p03_over_p02"] or 0.0
b_scaled = abs(a_ratio - 2.25) < 0.5
resolved_linear = a_pos and lin_sig_02 and lin_sig_03 and b_scaled
r4["resolved_linear_rate_supported"] = bool(resolved_linear)
r4["attack_summary"] = (
    "No reading of the data supports a resolved linear rate. "
    "(1) Constant-rate-plus-noise fails: the data is deterministic, "
    "and the quadratic term beats the linear-only model by F = 24.9 "
    "(p02) and 21.0 (p03) against F_crit(1,4) = 7.71. "
    "(2) The fitted linear terms are NEGATIVE (a = -7.4e-5 and "
    "-2.0e-4, |t| ~ 2.6-2.8), the wrong sign for the pre-registered "
    "pair read, and their ratio 2.77 does not match the B ratio "
    "2.25; the quadratic term ratio 2.68 and endpoint ratio 2.63 do "
    "not either. Everything scales super-linearly in eps, as a "
    "preparation-texture divergence would. "
    "(3) The finer round-1 pair (21 points, dt 0.005) shows the same "
    "acceleration: split-half ratio 25.2, a = -1.6e-4 (negative, "
    "t = -5.0), b t-stat +8.8. "
    "(4) The fieldless control azimuth drifts at 9.3e-4, i.e. 1.13x "
    "the round-1 claimed Omega of 8.249e-4, so a single-run linear "
    "read is background-dominated even without a field. The note's "
    "final verdict (unmeasurable in this configuration; dominant "
    "systematic is preparation-texture divergence; adiabatic ramp-on "
    "from a common fieldless state is the designed fix) survives "
    "the attack.")
r4["verdict"] = "REFUTED" if resolved_linear else "CONFIRMED"
report["R4"] = r4

with open(OUT, "w") as f:
    json.dump(report, f, indent=2)
print("wrote", OUT)
print(json.dumps(report, indent=2))
