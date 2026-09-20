"""M5.32 R23-1b adversarial audit: the c ladder on three boxes (claims D1 to D9).

EQUATIONS FIRST
---------------
Field: real symmetric 4 x 4 M(x), block-diagonal, spatial block S (3 x 3), sorted eigenvalues
(small, mid, top), vacuum S = (delta, delta, 1), delta = 0.3, M_00 = g = 8 in the vacuum.
    E      = E_curv + V4 + E_lin                       (definitions: m5_32_r23_1_audit.py)
    eps(x) = (mid - small) / 2,   eps(r) = shell mean of eps over |r - R| < half_width
    r_eps  = first R > r_ref where eps(R) R^p falls to level * eps(r_ref) r_ref^p,
             p = (sqrt(5) - 1) / 2, r_ref = 4.5, level = 1/2 (linear interpolation)
    R_rms  = sqrt( sum eps^2 r^2 / sum eps^2 ) over free cells with 5 <= r <= L/2 - 2h
    linear one-mode halo:  eps'' = eps / r^2 + beta r^2 eps,  beta = 5.81 c / (2 K),
             K = 8 (1 - delta)^2,  decaying solution sqrt(r) K_nu(sqrt(beta) r^2 / 2),
             nu = sqrt(5) / 4, cutoff radius beta^(-1/4)
    potential on the split (delta + eps, delta - eps):  V4 = k4 eps^4 (M_00 re-solved),
             -c L = 5.81 c eps^2 - c eps^4 / 2; the quartic restoring force exceeds the
             quadratic one where eps > eps_x = sqrt(5.81 c / (2 k4))
    uniaxial hedgehog S = delta I + (1 - delta) rhat rhat^T: curvature density
             8 (1 - delta)^4 / r^4 in the continuum (checked on the lattice here)

INDEPENDENCE
------------
Nothing is imported from m5_32_r23_1_cscan.py or any other script of the record. Reused from
the earlier independent auditor (m5_32_r23_1_audit.py, its own code, valid for any n): Par,
energy_parts (curv_density, pot_density), gradient_S, solve_u, make_seed, r_half_of,
degree_read, volume_disclinations, fit_candidates. Written here: the row matcher, all shell
profiles, the r_eps reader and its variants, R_rms, the linear-theory prediction, the split
coefficient k4, the hedgehog box integral, the perturbed L-BFGS stability probe, the label
replica, every threshold. All thresholds were fixed before the verdicts were read.

Run: python3 m5_32_r23_1b_audit.py     (about 2 min, at most 4 threads, no pools)
Writes: ../data/m5_32_r23_1b_audit.json
"""

import os

for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_k] = "4"

import json  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
import warnings  # noqa: E402

import numpy as np  # noqa: E402
from scipy.optimize import minimize  # noqa: E402
from scipy.special import kve  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import m5_32_r23_1_audit as EA  # noqa: E402  (the earlier independent auditor's own code)

DATA = os.path.join(HERE, "..", "data")
FIELDS = os.path.join(DATA, "m5_32_r23_1")
ROW_FILES = ("m5_32_r23_1_cscan.json", "m5_32_r23_1_cscan_extra.json",
             "m5_32_r23_1_cscan_extra2.json")  # fmt: skip
COLLECT_JSON = os.path.join(DATA, "m5_32_r23_1_cscan_collect.json")
OUT_JSON = os.path.join(DATA, "m5_32_r23_1b_audit.json")

DELTA, W1S, G8 = 0.3, 25.0, 8.0
P_TAIL = 0.5 * (np.sqrt(5.0) - 1.0)
NU = np.sqrt(5.0) / 4.0
K_LIN = 8.0 * (1.0 - DELTA) ** 2
A_LIN = (G8 + DELTA) * (1.0 - DELTA)
R_CORE = 4.0
BOXES = {
    "n32_L48": (32, 48.0, (0.0, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2)),
    "n48_L72": (48, 72.0, (3e-5, 1e-4, 3e-4, 1e-3, 3e-3)),
    "n48_L48": (48, 48.0, (3e-4, 3e-3, 1e-2, 3e-2)),
}
WINDOW = (3e-4, 1e-3, 3e-3)
CLAIM_R_EPS = {
    "n48_L72": {3e-5: 32.98, 1e-4: 28.11, 3e-4: 22.17, 1e-3: 16.11, 3e-3: 9.39},
    "n32_L48": {3e-5: None, 1e-4: None, 3e-4: None, 1e-3: 15.81, 3e-3: 9.36, 1e-2: 6.67,
                3e-2: 6.33},
    "n48_L48": {3e-4: 20.92, 3e-3: 8.82, 1e-2: 7.60, 3e-2: 7.23},
}  # fmt: skip
T0 = time.time()
CHECKS = []


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def add(cid, claim, method, value, expectation, fails_if, ok):
    CHECKS.append(
        {
            "id": cid,
            "claim": claim,
            "method": method,
            "value": value,
            "expectation": expectation,
            "fails_if": fails_if,
            "verdict": "PASS" if ok else "FAIL",
        }
    )
    log(f"{cid}: {'PASS' if ok else 'FAIL'}")


def tag_of(c, n, L, src=""):
    return f"rad_pin_d{DELTA:.4g}_w{W1S:.4g}_c{c:g}_n{n}_L{L:g}{src}"


def ck(c):
    return f"{c:g}"


# ======================= my own reads =======================
def shell_profile(eps, r, h, L, half_width, stat="mean"):
    """[R, eps shell statistic] on the 0.75 grid from 3.0 to L/2 - 2h."""
    out = []
    for R in np.arange(3.0, 0.5 * L, 0.75):
        if R > 0.5 * L - 2.0 * h + 1e-9:
            break
        sh = np.abs(r - R) < half_width
        if sh.any():
            v = eps[sh]
            if stat == "mean":
                s = v.mean()
            elif stat == "median":
                s = np.median(v)
            else:
                s = np.sqrt(np.mean(v * v))
            out.append([float(R), float(s)])
    return np.array(out)


def r_eps_read(prof, p_tail=P_TAIL, r_ref=4.5, level=0.5):
    R = prof[:, 0]
    y = prof[:, 1] * R**p_tail
    y0 = float(np.interp(r_ref, R, y))
    if y0 <= 0:
        return None
    k = np.where((R > r_ref) & (y < level * y0))[0]
    if len(k) == 0:
        return None
    i = k[0]
    return float(R[i - 1] + (level * y0 - y[i - 1]) * (R[i] - R[i - 1]) / (y[i] - y[i - 1]))


def seg_slopes(cs, vals):
    """local log slopes between neighbors, the least-squares slope, None if a value is missing."""
    if any(v is None or not np.isfinite(v) or v <= 0 for v in vals):
        return None
    lc, lv = np.log(np.array(cs)), np.log(np.array(vals, dtype=float))
    seg = [float(x) for x in np.diff(lv) / np.diff(lc)]
    lsq = float(np.polyfit(lc, lv, 1)[0]) if len(cs) >= 2 else None
    return {"values": [float(v) for v in vals], "segments": seg, "lsq": lsq}


def linear_theory_r_eps(c, grid):
    beta = A_LIN * c / (2.0 * K_LIN)
    z = 0.5 * np.sqrt(beta) * grid**2
    y = np.sqrt(grid) * kve(NU, z) * np.exp(-z)
    return r_eps_read(np.stack([grid, y], axis=1)), float(beta**-0.25)


def window_fit(prof, L, eps_cap=None, r_hi=None):
    """the earlier auditor's two-parameter candidates on 5 <= r <= L/2 - 6, eps above 1e-3 of
    its r = 5 value; eps_cap restricts to the linear tail (eps below the cap)."""
    e5 = float(np.interp(5.0, prof[:, 0], prof[:, 1]))
    hi = 0.5 * L - 6.0 if r_hi is None else r_hi
    k = (prof[:, 0] >= 5.0) & (prof[:, 0] <= hi + 1e-9) & (prof[:, 1] > 1e-3 * e5)
    if eps_cap is not None:
        k &= prof[:, 1] < eps_cap
    if k.sum() < 5:
        return None
    with warnings.catch_warnings(), np.errstate(all="ignore"):
        warnings.simplefilter("ignore")
        f = EA.fit_candidates(prof[k, 0], prof[k, 1])
    return {
        "points": int(k.sum()),
        "r_first_last": [float(prof[k, 0][0]), float(prof[k, 0][-1])],
        "knu_R_log": f["knu_cutoff"]["R_cut_log_fit"],
        "knu_R_lin": f["knu_cutoff"]["R_cut_lin_fit"],
        "knu_rms_log": f["knu_cutoff"]["rms_log"],
        "yukawa_R_log": f["yukawa"]["R_log_fit"],
        "yukawa_R_lin": f["yukawa"]["R_lin_fit"],
        "yukawa_rms_log": f["yukawa"]["rms_log"],
        "power_exponent_log": f["power"]["exponent_log_fit"],
        "power_rms_log": f["power"]["rms_log"],
    }


def load_row(tag, n, L, c, all_rows):
    P = EA.Par(DELTA, W1S, c, n, L)
    M = np.load(os.path.join(FIELDS, tag + ".npz"))["M"]
    S = np.ascontiguousarray(M[..., 1:, 1:])
    m00 = np.ascontiguousarray(M[..., 0, 0])
    ec, v4, lin = EA.energy_parts(S, m00, P)
    E = float((ec + v4 + lin).sum())
    best = None
    for fname, row in all_rows.get(tag, []):
        rel = abs(row["chunks"][-1]["E"] - E) / abs(E)
        if best is None or rel < best[0]:
            best = (rel, fname, row)
    lam, vec = np.linalg.eigh(S)
    return {
        "tag": tag, "P": P, "M": M, "S": S, "m00": m00, "ec": ec, "v4": v4, "lin": lin,
        "e": ec + v4 + lin, "E": E, "row": best[2], "row_file": best[1], "lam": lam, "vec": vec,
        "eps": 0.5 * (lam[..., 1] - lam[..., 0]), "last": best[2]["chunks"][-1],
        "n_candidates": len(all_rows.get(tag, [])),
    }  # fmt: skip


def main():
    all_rows = {}
    for fname in ROW_FILES:
        with open(os.path.join(DATA, fname)) as f:
            for t, r in json.load(f)["rows"].items():
                if r.get("status") == "OK" and r.get("chunks"):
                    all_rows.setdefault(t, []).append((fname, r))
    with open(COLLECT_JSON) as f:
        collect = json.load(f)

    D, missing = {}, []
    for box, (n, L, cs) in BOXES.items():
        for c in cs:
            t = tag_of(c, n, L)
            if not os.path.exists(os.path.join(FIELDS, t + ".npz")) or t not in all_rows:
                missing.append(t)
                continue
            D[(box, c)] = load_row(t, n, L, c, all_rows)
    t22 = tag_of(3e-4, 32, 48.0, "_from_r22s")
    D[("n32_L48_from_r22s", 3e-4)] = load_row(t22, 32, 48.0, 3e-4, all_rows)
    log(f"{len(D)} end fields loaded")

    # reads on every row
    for key, d in D.items():
        P = d["P"]
        h, L = P.h, P.L
        st = np.array(d["last"]["eps_profile_r_eps_top"])[:, :2]
        d["prof_stored"] = st
        d["prof"] = {
            "record_h/2": shell_profile(d["eps"], P.r, h, L, 0.5 * h),
            "fixed_0.75": shell_profile(d["eps"], P.r, h, L, 0.75),
            "wide_h": shell_profile(d["eps"], P.r, h, L, h),
            "narrow_0.375": shell_profile(d["eps"], P.r, h, L, 0.375),
            "median_h/2": shell_profile(d["eps"], P.r, h, L, 0.5 * h, "median"),
            "rms_h/2": shell_profile(d["eps"], P.r, h, L, 0.5 * h, "rms"),
        }
        pm = d["prof"]["record_h/2"]
        d["r_eps"] = {k: r_eps_read(p) for k, p in d["prof"].items()}
        d["r_eps"]["from_stored_profile"] = r_eps_read(st)
        d["r_eps_variants"] = {
            "ref4.5_half": r_eps_read(pm),
            "ref6_half": r_eps_read(pm, r_ref=6.0),
            "ref4.5_1/e": r_eps_read(pm, level=float(np.exp(-1.0))),
            "ref6_1/e": r_eps_read(pm, r_ref=6.0, level=float(np.exp(-1.0))),
            "ref4.5_half_no_power_weight": r_eps_read(pm, p_tail=0.0),
        }
        m = P.free & (P.r >= 5.0) & (P.r <= 0.5 * L - 2.0 * h)
        e2 = d["eps"][m] ** 2
        d["R_rms"] = float(np.sqrt(np.sum(e2 * P.r[m] ** 2) / np.sum(e2)))
        rh = EA.r_half_of(d["e"], P.r)
        d["r_half"], d["r_half_interp"] = rh[0], rh[1]
        d["discl"] = EA.volume_disclinations(d["vec"][..., :, 2], P)
        pf = d["prof"]["fixed_0.75"]
        d["eps_at"] = {
            f"{R:g}": {
                "shell_h/2": float(np.interp(R, pm[:, 0], pm[:, 1])),
                "shell_0.75": float(np.interp(R, pf[:, 0], pf[:, 1])),
            }
            for R in (4.5, 6.0, 9.0)
        }
        d["fit"] = window_fit(pm, L) if key[1] > 0 else None
    log("reads done")

    # ---------------- D1 ----------------
    claimed = [(b, c) for b in ("n48_L72", "n48_L48") for c in BOXES[b][2]]
    claimed += [("n32_L48", 1e-2), ("n32_L48", 3e-2)]
    blk = max(
        float(max(np.abs(d["M"][..., 0, 1:]).max(), np.abs(d["M"][..., 1:, 0]).max()))
        for d in D.values()
    )
    asym = max(float(np.abs(d["S"] - np.swapaxes(d["S"], -1, -2)).max()) for d in D.values())
    add(
        "D1.rows_exist",
        "every end field of the ladder exists, is block-diagonal and symmetric",
        "tags rebuilt from the claim text; max |M_0i| and max |S - S^T| over all cells",
        {"missing": missing, "max_M0i": blk, "max_asym": asym, "rows": len(D)},
        "no missing tag, M_0i exactly zero, asymmetry under 1e-15",
        "a missing end field, a nonzero M_0i or an asymmetric spatial block",
        not missing and blk == 0.0 and asym < 1e-15,
    )
    rep, worst = {}, 0.0
    for key in claimed:
        d = D[key]
        mine = {"E": d["E"], "E_curv": float(d["ec"].sum()), "V4": float(d["v4"].sum()),
                "E_lin": float(d["lin"].sum())}  # fmt: skip
        rel = {k: abs(v - d["last"][k]) / abs(d["last"][k]) for k, v in mine.items()}
        rep[d["tag"]] = {"mine": mine, "rel_diff": rel, "matched_row_file": d["row_file"]}
        worst = max(worst, max(rel.values()))
    add(
        "D1.energies",
        "stored E, E_curv, V4, E_lin of every n48 row and of the two new n32 rows reproduce "
        "from the end fields",
        "the earlier auditor's energy code (one-sided stencils, per-cell potential) on the end "
        "field against the last chunk of the row whose E is nearest",
        {"worst_rel_diff": worst, "rows": rep},
        "every relative difference under 1e-8",
        "any energy part of any of the 11 rows off by 1e-8 relative or more",
        worst < 1e-8,
    )
    fm, ok_f = {}, True
    for key in claimed:
        d = D[key]
        G = EA.gradient_S(d["S"], d["m00"], d["P"])
        mine = float((np.abs(G) * d["P"].free[..., None, None]).max())
        st = d["last"]["fmax_spatial"]
        fm[d["tag"]] = {"mine": mine, "stored": st, "rel": abs(mine - st) / st,
                        "gate": d["P"].gate, "label": d["row"]["label"]}  # fmt: skip
        ok_f &= abs(mine - st) / st < 1e-6
    add(
        "D1.fmax",
        "the stored fmax_spatial of the same 11 rows reproduces",
        "the earlier auditor's analytic gradient, max |G_ab| over free cells (2 pinned layers)",
        fm,
        "relative difference under 1e-6 on every row",
        "a stored fmax that my gradient does not reproduce to 1e-6",
        ok_f,
    )
    pins, ok_p = {}, True
    for key, d in D.items():
        P = d["P"]
        seed, _ = EA.make_seed(P, "rad")
        pin = ~P.free
        k = int(np.ceil(EA.PIN_DEPTH / P.h))
        inner = np.zeros_like(P.free)
        inner[k:-k, k:-k, k:-k] = True
        inner[k + 1 : -k - 1, k + 1 : -k - 1, k + 1 : -k - 1] = False
        dp = float(np.abs(d["S"][pin] - seed[pin]).max())
        df = float(np.abs(d["S"][inner] - seed[inner]).max())
        pins[d["tag"]] = {
            "pinned_layers": k,
            "pinned_depth": k * P.h,
            "max_diff_on_pinned": dp,
            "max_diff_on_first_free_layer": df,
        }
        ok_p &= dp < 1e-13 and df > 1e-6
    add(
        "D1.pinned_shell",
        "the outer shell is pinned at the analytic seed; its depth is 2 cells on every lattice "
        "(3.0 at h 1.5, 2.0 at h 1.0; the record's '1.6' is the nominal depth before rounding up)",
        "own seed against the end field on the 2 outer layers and on the first free layer",
        pins,
        "pinned cells equal the seed to 1e-13 and the first free layer has moved",
        "a pinned cell off the seed, or a first free layer still at the seed (a deeper pin)",
        ok_p,
    )
    cm, ok_c = {}, True
    for key, d in D.items():
        s = collect["rows"].get(d["tag"])
        if s is None:
            continue
        rel = abs(s["E"] - d["E"]) / abs(d["E"])
        cm[d["tag"]] = {
            "collect_E": s["E"],
            "end_field_E": d["E"],
            "rel": rel,
            "json_candidates": d["n_candidates"],
            "matched": d["row_file"],
        }
        ok_c &= rel < 1e-8
    add(
        "D1.collect_quotes_end_fields",
        "the summary the record quotes (collect JSON) carries, for every row, the E of the end "
        "field on disk (one row exists twice in the row files)",
        "collect E against my E of the end field",
        cm,
        "relative difference under 1e-8 on every row",
        "a collect row quoting a chunk that is not the end field on disk",
        ok_c,
    )

    # ---------------- D2 ----------------
    pd_, ok_pd = {}, True
    for key, d in D.items():
        a, b = d["prof"]["record_h/2"], d["prof_stored"]
        same = len(a) == len(b)
        diff = float(np.abs(a[:, 1] - b[:, 1]).max()) if same else None
        pd_[d["tag"]] = {"points": len(a), "stored_points": len(b), "max_abs_diff": diff}
        ok_pd &= same and diff < 1e-10
    add(
        "D2.profile_reproduces",
        "the stored eps profile of every row is the shell mean of (mid - small) / 2 on "
        "|r - R| < h / 2, R on the 0.75 grid up to L/2 - 2h",
        "own eigenvalues and own shells on the end field against the last chunk",
        pd_,
        "max abs difference under 1e-10 and the same number of points",
        "a stored profile I cannot rebuild from the end field",
        ok_pd,
    )
    rv, ok_rv = {}, True
    for box, table in CLAIM_R_EPS.items():
        for c, want in table.items():
            got = D[(box, c)]["r_eps"]["record_h/2"]
            rv[f"{box} c{ck(c)}"] = {"claimed": want, "mine": got}
            if want is None:
                ok_rv &= got is None
            else:
                ok_rv &= got is not None and abs(got - want) < 0.006
    add(
        "D2.r_eps_values",
        "the 16 quoted r_eps values (and the three None reads on n32 L48)",
        "own reader on own shell profile (half-width h / 2, the record's definition)",
        rv,
        "each within 0.006 of the quote, None where the record says None",
        "any value off by more than 0.006 or a None / number mismatch",
        ok_rv,
    )
    mean_defs = ("fixed_0.75", "wide_h", "narrow_0.375")
    sens = {}
    for box, table in CLAIM_R_EPS.items():
        for c, want in table.items():
            if want is None:
                continue
            d = D[(box, c)]
            base = d["r_eps"]["record_h/2"]
            alts = {k: d["r_eps"][k] for k in mean_defs + ("median_h/2", "rms_h/2")}
            dev = [None if alts[k] is None else abs(alts[k] - base) / base for k in mean_defs]
            sens[(box, c)] = {
                "record_h/2": base,
                "alternatives": alts,
                "max_rel_dev_mean_shells": None if any(v is None for v in dev) else max(dev),
                "a_mean_shell_reads_None": any(v is None for v in dev),
            }
    inner = [(b, c) for (b, c) in sens if c in WINDOW and not (b == "n48_L48" and c == 3e-4)]
    edge = [k for k in sens if k not in inner]

    def sens_ok(keys):
        return all(
            not sens[k]["a_mean_shell_reads_None"] and sens[k]["max_rel_dev_mean_shells"] < 0.05
            for k in keys
        )

    add(
        "D2.shell_sensitivity_window",
        "the r_eps values inside the window c 3e-4 to 3e-3 are properties of the field, not of "
        "the shell definition (they are quoted to 4 figures)",
        "the same reader on shell means with half-width 0.75, h and 0.375 (the record uses "
        "h / 2); the shell median and the shell rms are listed too but do not enter the verdict",
        {f"{b} c{ck(c)}": sens[(b, c)] for (b, c) in inner},
        "every mean-shell variant within 5 percent of the quoted value",
        "a shell width that moves a window value by 5 percent or more (that is 0.05 on a local "
        "slope over a factor 3 in c)",
        sens_ok(inner),
    )
    add(
        "D2.shell_sensitivity_edges",
        "the r_eps values at the ceiling (c <= 1e-4, and n48 L48 c 3e-4) and at the floor "
        "(c >= 1e-2) are stable against the shell definition",
        "same as above",
        {f"{b} c{ck(c)}": sens[(b, c)] for (b, c) in edge},
        "every mean-shell variant within 5 percent and never None",
        "a shell width that moves a value by 5 percent or more, or turns it into None (the "
        "read then sits on the end of the stored window)",
        sens_ok(edge),
    )
    nn = {}
    for c in (3e-5, 1e-4, 3e-4):
        a = D[("n32_L48", c)]
        p = a["prof"]["record_h/2"]
        y = p[:, 1] * p[:, 0] ** P_TAIL
        nn[ck(c)] = {
            "n32_L48_window_end": float(p[-1, 0]),
            "weighted_eps_at_window_end_over_value_at_4.5": float(
                y[-1] / np.interp(4.5, p[:, 0], y)
            ),
            "n48_L72_r_eps": D[("n48_L72", c)]["r_eps"]["record_h/2"],
        }
    add(
        "D2.none_means_beyond_window",
        "the three None reads on n32 L48 mean 'beyond the stored window (r = 21)', not 'no halo'",
        "the weighted profile at the last stored radius, and the L72 read at the same c",
        nn,
        "weighted profile still above one half at r = 21 and the L72 read above 21",
        "a None row whose L72 read lies inside the L48 window",
        all(v["n48_L72_r_eps"] > v["n32_L48_window_end"] for v in nn.values())
        and all(v["weighted_eps_at_window_end_over_value_at_4.5"] > 0.5 for v in nn.values()),
    )

    # ---------------- D3 ----------------
    fl = {}
    for box in ("n32_L48", "n48_L48"):
        for c in (3e-3, 1e-2, 3e-2):
            d = D[(box, c)]
            fl[f"{box} c{ck(c)}"] = {
                "r_eps": d["r_eps"]["record_h/2"],
                "r_half_interp": d["r_half_interp"],
                "R_rms_r>=5": d["R_rms"],
                "knu_R_log": d["fit"]["knu_R_log"],
                "yukawa_R_log": d["fit"]["yukawa_R_log"],
                "eps_at": d["eps_at"],
                "interior_spectrum_r1.5": d["last"]["shells_r_small_mid_top"][0],
            }
    top = [fl[f"{b} c{ck(c)}"] for b in ("n32_L48", "n48_L48") for c in (1e-2, 3e-2)]
    add(
        "D3.saturation_numbers",
        "for c >= 1e-2 r_eps reads 6.3 to 7.6 and the interpolated half-energy radius 5.6 to 5.8",
        "own reads on the four end fields",
        fl,
        "r_eps in [6.30, 7.65] and interpolated r_half in [5.55, 5.85] on all four rows",
        "a read outside the quoted range",
        all(6.30 <= v["r_eps"] <= 7.65 and 5.55 <= v["r_half_interp"] <= 5.85 for v in top),
    )
    sat = {}
    for box in ("n32_L48", "n48_L48"):
        for read in ("r_eps", "R_rms_r>=5", "knu_R_log", "yukawa_R_log", "r_half_interp"):
            a, b = fl[f"{box} c0.01"][read], fl[f"{box} c0.03"][read]
            sat[f"{box} {read}"] = float(np.log(b / a) / np.log(3.0))
    add(
        "D3.no_read_responds",
        "the top two ladder points carry no halo-size information: no size read still responds "
        "to c between 1e-2 and 3e-2",
        "local log slope between c 1e-2 and 3e-2 of five size reads on both spacings",
        sat,
        "every |slope| under 0.1 (a -1/4 law would give -0.25)",
        "a size read that still falls with a slope of magnitude 0.1 or more (the rows would "
        "then carry information and the floor would be the reader's, not the field's)",
        all(abs(v) < 0.1 for v in sat.values()),
    )
    col = {}
    for box in ("n32_L48", "n48_L48"):
        a, b = D[(box, 3e-3)]["eps_at"], D[(box, 1e-2)]["eps_at"]
        col[box] = {
            f"eps({R})_c3e-3_over_c1e-2": a[R]["shell_h/2"] / b[R]["shell_h/2"]
            for R in ("4.5", "6", "9")
        }
    lin_max = float(np.exp(0.5 * 4.5**2 * (10.0 / 3.0 - 1.0) / 5.12**2))
    col["largest_drop_at_r4.5_a_pure_retreat_allows"] = {
        "assumption": "K_nu envelope exp(-r^2 / 2 R^2), R = 5.12 at c 3e-3 shrinking as "
        "c^(-1/2) (the steepest law on the table)",
        "factor": lin_max,
    }
    add(
        "D3.retreat_or_collapse",
        "WORDING: 'the halo has retreated into the seed's core scale' (a halo of the same "
        "amplitude with a smaller radius)",
        "shell-mean eps at r = 4.5, 6, 9 at c = 3e-3 against c = 1e-2, both spacings; a "
        "retreating cutoff envelope bounds the drop at r = 4.5",
        col,
        "eps(4.5) falls by less than 5 (twice the bound of a pure retreat)",
        "eps at the core edge falls by 5 or more between c 3e-3 and 1e-2: the split does not "
        "shrink, its amplitude collapses (the halo is gone, not retreated)",
        all(col[b]["eps(4.5)_c3e-3_over_c1e-2"] < 5.0 for b in ("n32_L48", "n48_L48")),
    )
    lat = {}
    for c in (3e-3, 1e-2, 3e-2):
        a, b = D[("n32_L48", c)]["eps_at"], D[("n48_L48", c)]["eps_at"]
        lat[ck(c)] = {
            f"eps({R})_h1.5_over_h1.0": a[R]["shell_0.75"] / b[R]["shell_0.75"]
            for R in ("4.5", "6", "9")
        }
    lat["h_squared_ratio"] = 2.25
    add(
        "D3.floor_core_or_lattice",
        "the residual split at c >= 1e-2 is set by the core (a physical, spacing-independent "
        "profile), not by the lattice",
        "shell-mean eps (equal physical half-width 0.75) at r = 4.5, 6, 9 on L48 at h 1.5 "
        "against h 1.0; c = 3e-3 is the control where a real halo is present",
        lat,
        "ratio h 1.5 / h 1.0 inside [0.8, 1.25] at c 1e-2 and 3e-2, as it is at the control",
        "a ratio outside [0.8, 1.25] at c >= 1e-2: the residual eps scales with the spacing "
        "(about h^2 = 2.25), so it is discretization biaxiality of a uniaxial hedgehog",
        all(0.8 <= v <= 1.25 for c in ("0.01", "0.03") for v in lat[c].values()),
    )
    # stability probe: is the c = 1e-2 state an artifact of the symmetric seed?
    d = D[("n32_L48", 1e-2)]
    P, S0 = d["P"], d["S"]
    mcur = [d["m00"].copy()]

    def fg(x):
        Sx = x.reshape(S0.shape)
        mm = EA.solve_u(Sx, P, mcur[0], iters=12)
        mcur[0] = mm
        e = sum(float(q.sum()) for q in EA.energy_parts(Sx, mm, P))
        return e, (EA.gradient_S(Sx, mm, P) * P.free[..., None, None]).ravel()

    env = np.exp(-((P.r / 6.0) ** 2)) * P.free
    pert = np.zeros_like(S0)
    pert[..., 0, 0], pert[..., 1, 1] = env, -env
    lam_p = np.linalg.eigvalsh(S0 + 0.05 * pert)
    e_p = 0.5 * (lam_p[..., 1] - lam_p[..., 0])
    opts = {"maxiter": 300, "maxcor": 20, "ftol": 1e-15, "gtol": 1e-9}
    res = minimize(fg, (S0 + 0.05 * pert).ravel(), jac=True, method="L-BFGS-B", options=opts)
    lam_e = np.linalg.eigvalsh(res.x.reshape(S0.shape))
    e_e = 0.5 * (lam_e[..., 1] - lam_e[..., 0])
    sh6 = np.abs(P.r - 6.0) < 0.75
    stab = {
        "E_stored_field": d["E"],
        "E_after_perturbed_descent": float(res.fun),
        "rel_change": float((res.fun - d["E"]) / d["E"]),
        "iterations": int(res.nit),
        "max_abs_gradient_end": float(np.abs(res.jac).max()),
        "eps(6)_stored_field": float(d["eps"][sh6].mean()),
        "eps(6)_perturbed_start": float(e_p[sh6].mean()),
        "eps(6)_after_descent": float(e_e[sh6].mean()),
    }
    add(
        "D3.floor_state_is_stable",
        "the almost uniaxial state at c = 1e-2 is a minimum, not a symmetric stationary point "
        "that the uniaxial seed never left",
        "own L-BFGS descent (own reduced energy and gradient, M_00 re-solved) from the n32 end "
        "field plus a biaxial core perturbation 0.05 exp(-(r/6)^2) diag(1, -1, 0)",
        stab,
        "the descent returns: E not lower than the stored E by more than 1e-5 relative, and "
        "eps(6) within 10 percent of the stored field",
        "a lower energy state with a larger split is found (the floor would be a seed artifact)",
        stab["rel_change"] > -1e-5
        and abs(stab["eps(6)_after_descent"] / stab["eps(6)_stored_field"] - 1.0) < 0.10,
    )

    # ---------------- D4 ----------------
    cl = {}
    for c in (3e-5, 1e-4, 3e-4, 1e-3, 3e-3):
        d = D[("n48_L72", c)]
        free_half = 0.5 * d["P"].L - np.ceil(EA.PIN_DEPTH / d["P"].h) * d["P"].h
        cl[ck(c)] = {
            "r_eps": d["r_eps"]["record_h/2"],
            "free_half_box": float(free_half),
            "r_eps_over_free_half_box": d["r_eps"]["record_h/2"] / free_half,
            "last_stored_profile_radius": float(d["prof_stored"][-1, 0]),
            "label": d["row"]["label"],
            "fmax": d["last"]["fmax_spatial"],
            "gate": d["P"].gate,
        }
    ref = D[("n32_L48", 1e-3)]
    ref_ratio = ref["r_eps"]["record_h/2"] / (24.0 - 3.0)
    ref_agree = abs(ref["r_eps"]["record_h/2"] / D[("n48_L72", 1e-3)]["r_eps"]["record_h/2"] - 1)
    cl["demonstrated_safe_ratio"] = {
        "n32_L48_c1e-3_r_eps_over_free_half_box": ref_ratio,
        "its_disagreement_with_L72": float(ref_agree),
    }
    add(
        "D4.ceiling_rows",
        "on n48 L72 the c = 3e-5 and 1e-4 reads are wall-limited and those rows are not "
        "converged (FALLING, fmax 5.0e-3 and 3.5e-3 against the gate 1e-3)",
        "r_eps against the free half-box (36 minus the 2 pinned layers = 33); the largest "
        "ratio demonstrated harmless is the L48 c = 1e-3 row (agrees with L72 to 2 percent)",
        cl,
        "both reads above the demonstrated-safe ratio, both rows FALLING with fmax over gate",
        "a read at or below the demonstrated-safe ratio (it would not be wall-limited by "
        "this evidence) or a row at the gate",
        all(
            cl[k]["r_eps_over_free_half_box"] > ref_ratio
            and cl[k]["label"] == "FALLING"
            and cl[k]["fmax"] > cl[k]["gate"]
            for k in ("3e-05", "0.0001")
        )
        and abs(cl["3e-05"]["fmax"] - 5.0e-3) < 1e-4
        and abs(cl["0.0001"]["fmax"] - 3.5e-3) < 1e-4,
    )
    ag = {}
    for c in (3e-5, 1e-4, 3e-4, 1e-3, 3e-3):
        out = {}
        for nm in ("record_h/2", "wide_h"):
            a, b = D[("n32_L48", c)]["prof"][nm], D[("n48_L72", c)]["prof"][nm]
            for lo, hi in ((5.0, 12.0), (5.0, 15.0), (5.0, 21.0)):
                m = (a[:, 0] >= lo) & (a[:, 0] <= hi)
                rel = a[m, 1] / np.interp(a[m, 0], b[:, 0], b[:, 1]) - 1.0
                i = int(np.abs(rel).argmax())
                out[f"{nm} r{lo:g}-{hi:g}"] = {
                    "max_abs_rel_diff": float(np.abs(rel).max()),
                    "at_r": float(a[m, 0][i]),
                    "mean_rel_diff": float(rel.mean()),
                }
        ag[ck(c)] = out
    key15 = "record_h/2 r5-15"
    add(
        "D4.box_agreement_profile",
        "inside the readable window (c >= 3e-4) the two boxes agree on eps(r) within 10 percent "
        "over 5 <= r <= 15 (same spacing 1.5, same shells)",
        "max over shells of |eps_L48 / eps_L72 - 1| per c, also on 5-12 and 5-21 and with wide "
        "shells",
        ag,
        "under 10 percent at c = 3e-4, 1e-3 and 3e-3",
        "a window c where the boxes differ by 10 percent or more on 5 <= r <= 15",
        all(ag[ck(c)][key15]["max_abs_rel_diff"] < 0.10 for c in WINDOW),
    )

    # ---------------- D5 ----------------
    def col_of(box, fn, cs=WINDOW):
        return [fn(D[(box, c)]) for c in cs]

    Pk = EA.Par(DELTA, W1S, 0.0, 4, 6.0)
    k4s = {}
    for e in (0.02, 0.05, 0.10, 0.15):
        Ssp = np.diag([DELTA + e, DELTA - e, 1.0])[None]
        mm = EA.solve_u(Ssp, Pk, np.array([G8]))
        k4s[f"{e:g}"] = float(EA.pot_density(Ssp, mm, Pk)[0][0] / Pk.h**3 / e**4)
    k4 = k4s["0.05"]
    eps_x = {c: float(np.sqrt(A_LIN * c / (2.0 * k4))) for c in (3e-5, 1e-4) + WINDOW}
    TAIL = {
        c: window_fit(D[("n48_L72", c)]["prof"]["record_h/2"], 72.0, eps_cap=0.5 * eps_x[c])
        for c in eps_x
    }
    L72_18 = {}
    for c in WINDOW:
        L72_18[c] = window_fit(D[("n48_L72", c)]["prof"]["record_h/2"], 72.0, r_hi=18.0)
    reads = {
        "r_eps_ref4.5_half (the record's)": lambda d: d["r_eps_variants"]["ref4.5_half"],
        "r_eps_ref6_half": lambda d: d["r_eps_variants"]["ref6_half"],
        "r_eps_ref4.5_1/e": lambda d: d["r_eps_variants"]["ref4.5_1/e"],
        "r_eps_ref6_1/e": lambda d: d["r_eps_variants"]["ref6_1/e"],
        "r_eps_no_power_weight": lambda d: d["r_eps_variants"]["ref4.5_half_no_power_weight"],
        "r_eps_wide_shell_h": lambda d: d["r_eps"]["wide_h"],
        "r_eps_narrow_shell_0.375": lambda d: d["r_eps"]["narrow_0.375"],
        "R_rms_eps2_weighted_r>=5": lambda d: d["R_rms"],
        "r_half_interpolated": lambda d: d["r_half_interp"],
        "disclination_network_r_max": lambda d: d["discl"]["r_max"],
        "knu_R_log_fit": lambda d: d["fit"]["knu_R_log"],
        "knu_R_lin_fit": lambda d: d["fit"]["knu_R_lin"],
        "yukawa_R_log_fit": lambda d: d["fit"]["yukawa_R_log"],
    }
    table = {}
    for nm, fn in reads.items():
        table[nm] = {
            "n48_L72_window_3e-4_1e-3_3e-3": seg_slopes(WINDOW, col_of("n48_L72", fn)),
            "n48_L72_four_points_1e-4_to_3e-3": seg_slopes(
                (1e-4,) + WINDOW, col_of("n48_L72", fn, (1e-4,) + WINDOW)
            ),
            "n32_L48_two_points_1e-3_3e-3": seg_slopes(
                (1e-3, 3e-3), col_of("n32_L48", fn, (1e-3, 3e-3))
            ),
        }

    def bounded(f):
        """a K_nu cutoff, None when the fit did not run or ran away (above 400, no cutoff)."""
        return None if f is None or f["knu_R_log"] > 400.0 else f["knu_R_log"]

    for nm, src in (
        ("knu_R_log_fit_window_5_to_18", L72_18),
        ("knu_R_log_fit_linear_tail_only", TAIL),
    ):
        table[nm] = {
            "n48_L72_window_3e-4_1e-3_3e-3": seg_slopes(WINDOW, [bounded(src[c]) for c in WINDOW]),
            "n48_L72_four_points_1e-4_to_3e-3": None,
            "n32_L48_two_points_1e-3_3e-3": None,
        }
    rec = table["r_eps_ref4.5_half (the record's)"]["n48_L72_window_3e-4_1e-3_3e-3"]
    knu = table["knu_R_log_fit"]["n48_L72_window_3e-4_1e-3_3e-3"]
    add(
        "D5.quoted_slopes",
        "on n48 L72 the local log slopes of r_eps are -0.27 and -0.49, those of the fitted "
        "K_nu cutoff (10.83, 7.67, 5.12) are -0.29 and -0.37",
        "own r_eps; the earlier auditor's log-space K_nu fit on my profile (5 <= r <= 30, eps "
        "above 1e-3 of its r = 5 value)",
        {"r_eps": rec, "knu_R_log_fit": knu},
        "slopes within 0.01 of the quotes, K_nu radii within 1 percent",
        "a slope off by more than 0.01 or a K_nu radius off by more than 1 percent",
        abs(rec["segments"][0] + 0.27) < 0.01
        and abs(rec["segments"][1] + 0.49) < 0.01
        and abs(knu["segments"][0] + 0.29) < 0.01
        and abs(knu["segments"][1] + 0.37) < 0.01
        and all(abs(a / b - 1) < 0.01 for a, b in zip(knu["values"], (10.83, 7.67, 5.12))),
    )
    fq = {ck(c): D[("n48_L72", c)]["fit"] for c in WINDOW}
    met = {ck(c): fq[ck(c)]["knu_R_lin"] / fq[ck(c)]["knu_R_log"] for c in WINDOW}
    add(
        "D5.knu_fit_metric",
        "the fitted K_nu cutoff is a read of the profile (it does not depend on whether the "
        "least squares run in log space or in linear space)",
        "the same two-parameter K_nu form fitted both ways on the same points",
        {
            "fits": fq,
            "R_lin_over_R_log": met,
            "slopes_lin_fit": table["knu_R_lin_fit"]["n48_L72_window_3e-4_1e-3_3e-3"],
        },
        "the two radii within 15 percent on all three window rows, RMS log residual under 0.15",
        "a row where the radius moves by 15 percent or more with the metric, or where the "
        "K_nu form misses the profile by more than 15 percent RMS (the 'cutoff' is then a "
        "property of the fit, and so are its slopes)",
        all(abs(v - 1.0) < 0.15 for v in met.values())
        and all(fq[ck(c)]["knu_rms_log"] < 0.15 for c in WINDOW),
    )
    kb = {}
    for c in (1e-3, 3e-3):
        kb[ck(c)] = {
            "n32_L48_window_5_to_18": D[("n32_L48", c)]["fit"]["knu_R_log"],
            "n48_L72_window_5_to_30": D[("n48_L72", c)]["fit"]["knu_R_log"],
            "n48_L72_window_5_to_18": L72_18[c]["knu_R_log"],
        }
    kb["n48_L72_window_5_to_18_at_c3e-4"] = {
        "knu_R_log": L72_18[3e-4]["knu_R_log"],
        "note": "above 400 means no cutoff inside r <= 18: the fit runs away",
    }
    kb["slope_1e-3_to_3e-3"] = {
        k: float(np.log(kb["0.003"][k] / kb["0.001"][k]) / np.log(3.0)) for k in kb["0.001"]
    }
    add(
        "D5.knu_box_agreement",
        "the K_nu cutoff agrees between the boxes where r_eps does (c = 1e-3 and 3e-3)",
        "K_nu log fit on both boxes; the L72 fit repeated on the L48 window (5 to 18) to "
        "separate the window from the box",
        kb,
        "L48 against L72 (each on its own window, as the record quotes them) within 10 percent",
        "a K_nu cutoff that differs by 10 percent or more between boxes whose r_eps agree to "
        "2 percent",
        all(
            abs(kb[k]["n32_L48_window_5_to_18"] / kb[k]["n48_L72_window_5_to_30"] - 1) < 0.10
            for k in ("0.001", "0.003")
        ),
    )
    single = {}
    for nm, t in table.items():
        s = t["n48_L72_window_3e-4_1e-3_3e-3"]
        if s is not None:
            single[nm] = {
                "segments": s["segments"],
                "difference": abs(s["segments"][1] - s["segments"][0]),
                "lsq": s["lsq"],
                "compatible_with_one_exponent": abs(s["segments"][1] - s["segments"][0]) <= 0.10,
            }
    allseg = [x for v in single.values() for x in v["segments"]]
    add(
        "D5.not_a_power_law",
        "inside the window the radius is not a power law of c: no size read gives one exponent",
        "two local slopes per read on n48 L72 (15 reads: r_eps with other reference radii, "
        "levels, weights and shells; the eps^2-weighted rms radius; the interpolated "
        "half-energy radius; the disclination network's outer radius; K_nu and Yukawa scales, "
        "the K_nu fit also on the linear tail alone and on the 5 to 18 window, which has no "
        "cutoff at c = 3e-4 and drops out)",
        {
            "reads": single,
            "slope_range_over_all_reads": [min(allseg), max(allseg)],
            "full_table": table,
        },
        "every read has local slopes differing by more than 0.10 (the 2 percent convergence "
        "uncertainty of a radius is 0.02 on a slope)",
        "a read whose two local slopes agree within 0.10: that read cannot exclude a single "
        "exponent",
        not any(v["compatible_with_one_exponent"] for v in single.values()),
    )
    lo_, hi_ = min(rec["segments"]), max(rec["segments"])
    lo_k, hi_k = min(knu["segments"]), max(knu["segments"])
    add(
        "D5.inside_the_spread",
        "WORDING: 'neither -1/4 nor -1/2 is certified and both lie inside the spread'",
        "the literal intervals spanned by the quoted slopes",
        {
            "r_eps_interval": [lo_, hi_],
            "knu_interval": [lo_k, hi_k],
            "distance_of_-1/4_outside_r_eps_interval": -0.25 - hi_,
            "distance_of_-1/2_outside_r_eps_interval": lo_ + 0.5,
            "all_reads_interval": [min(allseg), max(allseg)],
        },
        "-0.25 and -0.5 inside [min, max] of the two quoted r_eps slopes",
        "either exponent outside the interval of the quoted r_eps slopes (true only of the "
        "union over other reads, or within the 0.02 convergence uncertainty)",
        lo_ <= -0.5 <= hi_ and lo_ <= -0.25 <= hi_,
    )
    floor_r = D[("n32_L48", 3e-2)]["r_eps"]["record_h/2"]
    amp = {ck(c): D[("n48_L72", c)]["eps_at"]["4.5"]["shell_h/2"] for c in (1e-4,) + WINDOW}
    net = {ck(c): D[("n48_L72", c)]["discl"] for c in (1e-4,) + WINDOW}
    ends = {
        "upper_end_c3e-3": {
            "r_eps": rec["values"][2],
            "floor_r_eps_same_spacing": floor_r,
            "ratio_to_floor": rec["values"][2] / floor_r,
            "eps(4.5)_by_c": amp,
            "eps(4.5)_c3e-3_over_c3e-4": amp["0.003"] / amp["0.0003"],
        },
        "lower_end_c3e-4": {
            "r_eps_over_free_half_box": rec["values"][0] / 33.0,
            "disclination_network_by_c": net,
        },
    }
    add(
        "D5.window_ends_clean",
        "the readable window is c in [3e-4, 3e-3]: both end points are free of the floor and "
        "of structure in the read region",
        "upper end: r_eps against the floor value at the same spacing, and the split amplitude "
        "at r = 4.5 against the plateau; lower end: the outer radius of the disclination "
        "network against the start of the read region (r_ref 4.5, fit window from 5)",
        ends,
        "r_eps(3e-3) at least twice the floor, eps(4.5) within 20 percent of its c = 3e-4 "
        "value, the network inside r = 5 on all three rows",
        "an end point within a factor 2 of the floor, an amplitude already down by more than "
        "20 percent, or a network reaching into the read region",
        ends["upper_end_c3e-3"]["ratio_to_floor"] >= 2.0
        and ends["upper_end_c3e-3"]["eps(4.5)_c3e-3_over_c3e-4"] > 0.8
        and all((net[ck(c)]["r_max"] or 0.0) < 5.0 for c in WINDOW),
    )
    dr = {}
    for c in WINDOW:
        ch = D[("n48_L72", c)]["row"]["chunks"]
        last_it = ch[-1]["iters"]
        win = [q for q in ch if q["iters"] >= last_it - 1000 and q["iters"] > 0]
        vals = [r_eps_read(np.array(q["eps_profile_r_eps_top"])[:, :2]) for q in win]
        dr[ck(c)] = {
            "iters": [q["iters"] for q in win],
            "r_eps": vals,
            "E": [q["E"] for q in win],
            "fmax": [q.get("fmax_spatial") for q in win],
            "change_over_window_rel": (vals[-1] - vals[0]) / vals[-1],
            "last_chunk_step_rel": (vals[-1] - vals[-2]) / vals[-1],
            "monotone": bool(np.all(np.diff(vals) >= 0) or np.all(np.diff(vals) <= 0)),
            "label": D[("n48_L72", c)]["row"]["label"],
        }
    x48 = D[("n32_L48", 1e-3)]["r_eps"]["record_h/2"]
    dr["cross_check_c1e-3"] = {
        "n32_L48_AT_GATE_r_eps": x48,
        "n48_L72_FALLING_r_eps": rec["values"][1],
        "rel_diff": abs(x48 / rec["values"][1] - 1.0),
    }
    add(
        "D5.falling_rows_radii",
        "the FALLING L72 rows at c = 3e-4 and 1e-3 (fmax 1.1e-3 and 6.7e-3) are converged "
        "enough for their radii",
        "r_eps recomputed from every stored chunk of the last 1000 iterations; the converged "
        "L48 row at c = 1e-3 as a cross-check",
        dr,
        "change over the last 1000 iterations under 2 percent, last chunk step under 0.5 "
        "percent, L48 and L72 within 3 percent at c = 1e-3",
        "a radius still moving by 2 percent per 1000 iterations or more",
        all(
            abs(dr[ck(c)]["change_over_window_rel"]) < 0.02
            and abs(dr[ck(c)]["last_chunk_step_rel"]) < 0.005
            for c in (3e-4, 1e-3)
        )
        and dr["cross_check_c1e-3"]["rel_diff"] < 0.03,
    )
    bx = {}
    for c in (3e-5, 1e-4):
        a, b = D[("n32_L48", c)], D[("n48_L72", c)]
        sa = collect["rows"][a["tag"]]
        sb = collect["rows"][b["tag"]]
        bx[ck(c)] = {
            "r_half_cell_reader": [
                a["r_half"],
                b["r_half"],
                abs(a["r_half"] - b["r_half"]) / b["r_half"],
            ],
            "r_half_interpolated": [
                a["r_half_interp"],
                b["r_half_interp"],
                abs(a["r_half_interp"] - b["r_half_interp"]) / b["r_half_interp"],
            ],
            "labels": [a["row"]["label"], b["row"]["label"]],
            "fmax_over_gate": [a["last"]["fmax_spatial"] / 1e-3, b["last"]["fmax_spatial"] / 1e-3],
            "collect_r_half_drift": [sa["r_half_drift"], sb["r_half_drift"]],
            "collect_both_converged": collect["box_check_L48_vs_L72"][str(c)]["both_converged"],
        }
    add(
        "D5.box_check_numbers",
        "the L48 and L72 half-energy radii differ by 11.7 percent at c = 3e-5 and by 7.1 "
        "percent at c = 1e-4 (L72 row unconverged); the pre-registered label is "
        "HALO_BOX_LIMITED",
        "own half-energy radius (cell reader and interpolated) on the four end fields",
        {"rows": bx, "collect_label": collect["label_R23_1"]},
        "cell-reader differences within 0.002 of 0.117 and 0.071, label HALO_BOX_LIMITED",
        "a difference I cannot reproduce or another label in the collect JSON",
        abs(bx["3e-05"]["r_half_cell_reader"][2] - 0.117) < 0.002
        and abs(bx["0.0001"]["r_half_cell_reader"][2] - 0.071) < 0.002
        and collect["label_R23_1"] == "HALO_BOX_LIMITED",
    )
    # replica of the pre-registered label rule with the c = 1e-4 drift gate relaxed
    cf = {}
    rows32 = {c: collect["rows"][tag_of(c, 32, 48.0)] for c in BOXES["n32_L48"][2] if c > 0}
    for read in ("r_half", "r_eps"):
        pts = [
            (c, s[read])
            for c, s in rows32.items()
            if c >= 1e-4 and s["CONVERGED_FOR_FIT"] and s[read] is not None
        ]
        cf[read] = {
            "points": pts,
            "slope": float(
                np.polyfit(np.log([p[0] for p in pts]), np.log([p[1] for p in pts]), 1)[0]
            ),
        }
    sl = [cf[k]["slope"] for k in ("r_half", "r_eps")]
    if all(-0.6 <= v <= -0.4 for v in sl):
        cf_label = "HALO_HALF"
    elif all(-0.3 <= v <= -0.2 for v in sl):
        cf_label = "HALO_QUARTER"
    else:
        cf_label = "HALO_OTHER"
    cf["label_if_the_c1e-4_box_check_had_counted"] = cf_label
    cf["margin"] = {
        "L72_c1e-4_r_half_drift": collect["rows"][tag_of(1e-4, 48, 72.0)]["r_half_drift"],
        "drift_gate": 0.02,
        "box_difference_at_c1e-4": bx["0.0001"]["r_half_cell_reader"][2],
        "box_gate": 0.10,
    }
    add(
        "D5.label_margin",
        "WORDING: HALO_BOX_LIMITED is a robust outcome of the pre-registered rule",
        "replica of the rule's last step from the collect summary: the c = 1e-4 box check "
        "passes its 10 percent gate (7.1) and is excluded only by the drift flag of the L72 "
        "row (r_half drift 2.35 percent against 2.0, read on a cell-quantized radius); the "
        "label the rule gives if that one flag flips",
        cf,
        "the same label HALO_BOX_LIMITED",
        "another label: the outcome then hangs on a 0.35 point margin of one drift read, and "
        "the alternative label would rest on a fit through the floor rows c >= 1e-2",
        cf_label == "HALO_BOX_LIMITED",
    )

    # ---------------- D6 ----------------
    grid_r = np.arange(3.0, 60.1, 0.75)
    grid_f = np.arange(3.0, 60.0, 0.001)
    quoted = {3e-5: 18.5, 1e-4: 14.0, 3e-4: 11.0, 1e-3: 8.6, 3e-3: 7.2}
    lt = {}
    for c, q in quoted.items():
        a, rc = linear_theory_r_eps(c, grid_r)
        b, _ = linear_theory_r_eps(c, grid_f)
        meas = D[("n48_L72", c)]["r_eps"]["record_h/2"]
        lt[ck(c)] = {"quoted": q, "on_reader_grid": a, "on_fine_grid": b, "beta_-1/4": rc,
                     "measured_L72": meas, "measured_over_predicted": meas / a}  # fmt: skip
    ratios = [v["measured_over_predicted"] for v in lt.values()]
    add(
        "D6.arithmetic",
        "the linearized one-mode equation predicts r_eps = 18.5, 14.0, 11.0, 8.6, 7.2; the "
        "L72 reads are 1.3 to 2.0 times larger; a factor 2 in R is a factor 16 in beta",
        "the record's reader on sqrt(r) K_nu(sqrt(beta) r^2 / 2), beta = 5.81 c / (2 * 3.92), "
        "on the reader's 0.75 grid and on a 0.001 grid",
        {"rows": lt, "ratio_range": [min(ratios), max(ratios)], "two_to_the_fourth": 2.0**4},
        "predictions within 0.05 of the quotes; ratio range within 0.05 of [1.3, 2.0]",
        "a prediction or a ratio I cannot reproduce",
        all(abs(v["on_reader_grid"] - v["quoted"]) < 0.05 for v in lt.values())
        and abs(min(ratios) - 1.3) < 0.05
        and abs(max(ratios) - 2.0) < 0.05,
    )
    nl = {}
    for c in quoted:
        d = D[("n48_L72", c)]
        ex = eps_x[c]
        p = d["prof"]["record_h/2"]
        below = np.where((p[:, 0] >= 4.5) & (p[:, 1] < ex))[0]
        e6 = d["eps_at"]["6"]["shell_h/2"]
        nl[ck(c)] = {
            "eps_x": ex,
            "eps(4.5)": d["eps_at"]["4.5"]["shell_h/2"],
            "eps(6)": e6,
            "quartic_over_quadratic_force_at_r6": (e6 / ex) ** 2,
            "eps(6)_over_gap_0.7": e6 / (1.0 - DELTA),
            "minus_cL_eps4_over_eps2_term_at_r6": e6**2 / (2.0 * A_LIN),
            "first_shell_radius_with_eps_below_eps_x": (
                float(p[below[0], 0]) if len(below) else None
            ),
            "r_eps": d["r_eps"]["record_h/2"],
        }
    add(
        "D6.reason_nonlinearity",
        "possible reason 1 is absent: in the read region the split is small enough for the "
        "linearized potential term (5.81 c eps^2 dominates V4 = k4 eps^4)",
        "k4 from my own per-cell potential on the split with M_00 re-solved; eps_x = sqrt(5.81 "
        "c / (2 k4)) is where the two restoring forces are equal; shell-mean eps at r = 6",
        {"k4_by_eps": k4s, "rows": nl},
        "force ratio (eps / eps_x)^2 at r = 6 under 0.25 at every c of the five",
        "a c where the quartic V4 force at r = 6 is a quarter of the linear one or more: the "
        "linear equation cannot hold from r_ref on, and nonlinearity is a live reason for the "
        "prefactor",
        all(v["quartic_over_quadratic_force_at_r6"] < 0.25 for v in nl.values()),
    )
    ov = {ck(c): D[("n48_L72", c)]["discl"] for c in quoted}
    add(
        "D6.reason_disclinations",
        "possible reason 2 is absent: the disclination network of the top eigenvector stays "
        "inside the read region's inner edge (r_ref 4.5, fit window from 5)",
        "frustrated elementary squares over the free region (the earlier auditor's counter), "
        "their largest radius",
        ov,
        "r_max under 5 at every c of the five",
        "a c where the network reaches r = 5 or beyond (the one-mode, defect-free ansatz is "
        "then not the state that was measured)",
        all((v["r_max"] or 0.0) < 5.0 for v in ov.values()),
    )
    tf = {}
    for c in quoted:
        d = D[("n48_L72", c)]
        f = TAIL[c]
        rc = lt[ck(c)]["beta_-1/4"]
        tf[ck(c)] = {
            "tail_fit": f,
            "beta_-1/4": rc,
            "R_tail_over_theory": None if f is None else f["knu_R_log"] / rc,
            "R_full_window_over_theory": d["fit"]["knu_R_log"] / rc,
            "implied_beta_over_theory_beta": None if f is None else (f["knu_R_log"] / rc) ** -4,
        }
    rt = [v["R_tail_over_theory"] for v in tf.values() if v["R_tail_over_theory"] is not None]
    add(
        "D6.verdict_linear_tail",
        "'the linearized one-mode equation does not describe the measured halo quantitatively' "
        "survives a fairer test: the K_nu cutoff fitted ONLY on the linear tail (eps below half "
        "of eps_x) against beta^(-1/4)",
        "log-space K_nu fit restricted to shells with eps < eps_x / 2, r >= 5; ratio of the "
        "fitted cutoff to the theoretical one",
        tf,
        "the record's statement stands: at least one ratio outside [0.9, 1.1]",
        "every tail ratio inside [0.9, 1.1]: the equation would then describe the tail and the "
        "factor 1.3 to 2.0 would be the r_eps reader starting inside the nonlinear plateau",
        any(abs(v - 1.0) > 0.10 for v in rt),
    )

    # ---------------- D7 ----------------
    a, b, s = D[("n32_L48", 3e-3)], D[("n48_L72", 3e-3)], D[("n48_L48", 3e-3)]
    vals = {
        "r_eps": [x["r_eps"]["record_h/2"] for x in (a, b, s)],
        "E": [x["E"] for x in (a, b, s)],
        "E_L72_minus_E_L48": b["E"] - a["E"],
        "rel": (b["E"] - a["E"]) / a["E"],
    }
    add(
        "D7.values",
        "at c = 3e-3: r_eps 9.36, 9.39, 8.82 and E 6.2314, 6.5104, 6.3330; the L48 to L72 "
        "energy difference is 0.279 (4.5 percent)",
        "own reads and own energies",
        vals,
        "r_eps within 0.006, E within 1e-4, difference within 0.001 of 0.279",
        "any quoted number off",
        all(abs(x - y) < 0.006 for x, y in zip(vals["r_eps"], (9.36, 9.39, 8.82)))
        and all(abs(x - y) < 1e-4 for x, y in zip(vals["E"], (6.2314, 6.5104, 6.3330)))
        and abs(vals["E_L72_minus_E_L48"] - 0.279) < 0.001,
    )
    hh = {}
    for nm, x in (("n32_L48", a), ("n48_L72", b)):
        _, ext = EA.make_seed(x["P"], "rad")
        hh[nm] = EA.curv_density(ext, x["P"].h)
    r72 = b["P"].r
    mid = (r72 > 10.0) & (r72 < 18.0)
    dens = float(np.mean(hh["n48_L72"][mid] * r72[mid] ** 4 / b["P"].h ** 3) / (1.0 - DELTA) ** 4)
    add(
        "D7.hedgehog_density",
        "the curvature energy density of the uniaxial hedgehog is 8 (1 - delta)^4 / r^4",
        "the earlier auditor's lattice curvature density on S = delta I + (1 - delta) rhat "
        "rhat^T, mean of density r^4 / (1 - delta)^4 over cells with 10 < r < 18 at h = 1.5",
        {"coefficient": dens},
        "coefficient within 3 percent of 8",
        "a lattice coefficient off 8 by more than 3 percent",
        abs(dens / 8.0 - 1.0) < 0.03,
    )
    N = 240
    xg = (np.arange(N) + 0.5) / N
    Xg = np.stack(np.meshgrid(xg, xg, xg, indexing="ij"), axis=-1)
    rg = np.sqrt(np.sum(Xg**2, axis=-1))
    k_cube = float(4.0 * np.pi - 8.0 * np.sum((rg > 1.0) * rg**-4.0) / N**3)
    cont = 8.0 * (1.0 - DELTA) ** 4 * k_cube * (1.0 / 24.0 - 1.0 / 36.0)
    dec = {"continuum_cube_factor_k": k_cube, "continuum_cube_minus_cube": float(cont)}
    for R in (9.0, 12.0, 15.0):
        ia, ib = a["P"].r < R, b["P"].r < R
        dec[f"split_at_r{R:g}"] = {
            "rows_inside_diff": float(b["e"][ib].sum() - a["e"][ia].sum()),
            "rows_outside_diff": float(b["e"][~ib].sum() - a["e"][~ia].sum()),
            "hedgehog_outside_diff_lattice": float(
                hh["n48_L72"][~ib].sum() - hh["n32_L48"][~ia].sum()
            ),
        }
    outc = np.abs(b["P"].X).max(axis=-1) > 24.0
    dec["L72_cells_outside_the_L48_cube"] = {
        "rows": float(b["e"][outc].sum()),
        "hedgehog": float(hh["n48_L72"][outc].sum()),
    }
    hd = dec["split_at_r12"]["hedgehog_outside_diff_lattice"]
    add(
        "D7.coulomb_tail_accounts",
        "the 0.279 energy difference between L48 and L72 at c = 3e-3 is the Coulomb tail of "
        "the hedgehog between the two walls",
        "lattice hedgehog curvature energy outside r = 12 on each lattice (boundary stencils "
        "and pinned shells as the rows have them), L72 minus L48; the continuum cube-minus-"
        "cube integral; the rows' own difference split inside / outside r = 12",
        dec,
        "the hedgehog difference within 5 percent of 0.279 and the rows' difference inside "
        "r = 12 under 5 percent of 0.279",
        "a hedgehog tail that misses 0.279 by more than 5 percent, or an interior that "
        "differs between the boxes by more than 0.014",
        abs(hd / vals["E_L72_minus_E_L48"] - 1.0) < 0.05
        and abs(dec["split_at_r12"]["rows_inside_diff"]) < 0.05 * 0.279,
    )

    # ---------------- D8 ----------------
    deg = {}
    ok_clean, ok_low, ok_class, flips = True, True, True, {}
    for key, d in D.items():
        P = d["P"]
        topv = d["vec"][..., :, 2]
        gap = d["lam"][..., 2] - d["lam"][..., 1]
        rec_ = {"c": key[1], "box": key[0], "network": d["discl"], "cubes": {}}
        any_frus = False
        for hv in (6.0, 9.0, 12.0):
            k = int(round(hv / P.h))
            c0 = P.n // 2
            m1 = EA.degree_read(topv, gap, P, c0 - k, c0 + k)
            m2 = EA.degree_read(topv, gap, P, c0 - k - 1, c0 + k)
            st = d["last"]["degree_top"][f"{hv:g}"]
            rec_["cubes"][f"{hv:g}"] = {
                "stored_degree_conflicts": st,
                "record_cube": {q: m1[q] for q in ("frustrated_plaquettes", "conflicts_bfs",
                                                   "conflicts_rhat", "degree_bfs")},
                "box_centered_cube": {q: m2[q] for q in ("frustrated_plaquettes",
                                                         "conflicts_bfs", "conflicts_rhat",
                                                         "degree_bfs")},
            }  # fmt: skip
            any_frus |= m1["frustrated_plaquettes"] > 0 or m2["frustrated_plaquettes"] > 0
            ok_class &= (st[1] == 0) == (m1["frustrated_plaquettes"] == 0)
            if (m1["frustrated_plaquettes"] == 0) != (m2["frustrated_plaquettes"] == 0):
                flips[f"{d['tag']} r{hv:g}"] = [
                    m1["frustrated_plaquettes"],
                    m2["frustrated_plaquettes"],
                    m2["degree_bfs"],
                ]
            if key[1] >= 1e-3:
                for m in (m1, m2):
                    ok_clean &= (
                        m["frustrated_plaquettes"] == 0
                        and m["conflicts_bfs"] == 0
                        and m["degree_bfs"] is not None
                        and abs(abs(m["degree_bfs"]) - 1.0) < 1e-6
                    )
                ok_clean &= abs(abs(st[0]) - 1.0) < 1e-6 and st[1] == 0
        if key[1] <= 3e-4:
            ok_low &= any_frus
        deg[d["tag"]] = rec_
    add(
        "D8.clean_degree_high_c",
        "clean unit |degree| with zero orientation conflicts on the r 6, 9, 12 cubes for every "
        "rad-seed row with c >= 1e-3 on every box and spacing",
        "the earlier auditor's reader (gauge-invariant frustrated squares, flood-fill "
        "orientation, solid-angle degree) on the record's cube and on the box-centered cube",
        deg,
        "zero frustrated squares, zero conflicts, |degree| = 1 to 1e-6 on both cubes, 9 rows",
        "a frustrated square, a conflict or a non-unit degree on any cube of a c >= 1e-3 row",
        ok_clean,
    )
    add(
        "D8.conflicts_only_low_c",
        "conflicts appear only at c <= 3e-4, and every row there has them; the stored zero / "
        "nonzero pattern is the gauge-invariant one",
        "disclination piercings per cube on all 18 rows (values in D8.clean_degree_high_c)",
        {"rows_c<=3e-4_all_pierced": ok_low, "stored_pattern_matches_piercings": ok_class},
        "every c <= 3e-4 row pierced on at least one cube; stored conflicts zero exactly where "
        "my piercing count is zero",
        "a clean c <= 3e-4 row, or a stored clean cube with piercings (or the reverse)",
        ok_low and ok_class,
    )
    add(
        "D8.cube_centering",
        "WORDING: the clean / pierced class of a cube does not depend on the record's cube "
        "being centered half a cell off the box center",
        "the same read on the box-centered cube (one cell wider on the low side)",
        {"cubes_that_change_class [record, centered, centered degree]": flips},
        "no cube changes class",
        "a cube pierced on the record's cube and clean on the box-centered one, or the reverse",
        not flips,
    )
    netr = {d["tag"]: d["discl"]["r_max"] for d in D.values()}
    jump = {box: {ck(c): D[(box, c)]["discl"]["r_max"] for c in BOXES[box][2]} for box in BOXES}
    a_, b_ = jump["n48_L72"]["0.0003"], jump["n48_L72"]["0.001"]
    add(
        "D8.network_radius_is_smooth",
        "WORDING: the c ladder is one family of states (the disclination network shrinks "
        "smoothly with c, like the halo radius)",
        "outer radius of the network per row; the largest local log slope between neighbors",
        {
            "r_max_by_row": netr,
            "by_box": jump,
            "L72_slope_3e-4_to_1e-3": float(np.log(b_ / a_) / np.log(1e-3 / 3e-4)),
        },
        "local slope of the network radius between c 3e-4 and 1e-3 no steeper than -0.6",
        "a slope steeper than -0.6 between two window points: the network collapses into the "
        "core between them, so the window straddles a structural change",
        float(np.log(b_ / a_) / np.log(1e-3 / 3e-4)) > -0.6,
    )

    # ---------------- D9 ----------------
    Siso = ((1.0 + 2.0 * DELTA) / 3.0 * np.eye(3))[None]
    V0 = float(EA.pot_density(Siso, np.array([G8]), Pk)[0][0] / Pk.h**3)
    c_phys = 3.1e-3 * V0 / (2.0 * A_LIN)
    ex = {}
    anchors = {
        "r_eps": dict(zip(WINDOW, rec["values"])),
        "knu_R_log": dict(zip(WINDOW, knu["values"])),
    }
    for read, tab in anchors.items():
        for c, r0 in tab.items():
            for s_ in (-0.25, -0.5):
                v = r0 * (c_phys / c) ** s_
                ex[f"{read} anchor c{ck(c)} slope {s_}"] = {
                    "lattice_units": v,
                    "core_radii": v / R_CORE,
                }
    q = ex["r_eps anchor c0.001 slope -0.25"], ex["r_eps anchor c0.001 slope -0.5"]
    add(
        "D9.arithmetic",
        "V_0 = 0.019253, c_phys = 5.1e-6; anchored at r_eps(1e-3) = 16.1 the slopes -1/4 and "
        "-1/2 give 60 and 225 lattice units, 15 and 56 core radii",
        "own V4 of the isotropic block at M_00 = g; the conversion as stated",
        {
            "V_0": V0,
            "c_phys": c_phys,
            "anchored_at_c1e-3": q,
            "decades_below_the_smallest_window_c": float(np.log10(3e-4 / c_phys)),
        },
        "V_0 within 1e-6, c_phys within 0.05e-6 of 5.1e-6, 60 and 225 within 1, 15 and 56 "
        "within 0.5",
        "any of the quoted numbers off",
        abs(V0 - 0.019253) < 1e-6
        and abs(c_phys - 5.1e-6) < 0.05e-6
        and abs(q[0]["lattice_units"] - 60) < 1
        and abs(q[1]["lattice_units"] - 225) < 1
        and abs(q[0]["core_radii"] - 15) < 0.5
        and abs(q[1]["core_radii"] - 56) < 0.5,
    )
    lo_all = min(v["core_radii"] for v in ex.values())
    hi_all = max(v["core_radii"] for v in ex.values())
    re_ = [v["core_radii"] for k, v in ex.items() if k.startswith("r_eps")]
    add(
        "D9.anchor_robustness",
        "WORDING: '15 to 56 core radii' is a fair statement of the extrapolation range",
        "the same two exponents anchored at each of the three window points, with r_eps and "
        "with the K_nu cutoff",
        {
            "all": ex,
            "r_eps_range_core_radii": [min(re_), max(re_)],
            "all_reads_range_core_radii": [lo_all, hi_all],
        },
        "both ends of the range move by less than 20 percent under the anchor and the reader",
        "an end of the range that moves by 20 percent or more (the range then states the "
        "anchor, not the extrapolation)",
        abs(lo_all / 15.0 - 1.0) < 0.20
        and abs(hi_all / 56.0 - 1.0) < 0.20
        and abs(min(re_) / 15.0 - 1.0) < 0.20
        and abs(max(re_) / 56.0 - 1.0) < 0.20,
    )
    add(
        "D9.nothing_measured",
        "all extrapolated radii lie far outside both boxes, so nothing at the physical point "
        "was measured",
        "the smallest extrapolated radius over all anchors and readers against the larger "
        "free half-box (33)",
        {
            "smallest_lattice_units": lo_all * R_CORE,
            "larger_free_half_box": 33.0,
            "note": "the K_nu anchors give 25 to 30 lattice units, inside the L72 box; a cutoff "
            "radius is not an r_eps (r_eps is about 2 cutoff radii), so this does not put the "
            "physical point in reach, but 'all far outside both boxes' holds for r_eps only",
        },
        "every r_eps extrapolation above 33",
        "an r_eps extrapolation that fits inside the L72 free region",
        min(re_) * R_CORE > 33.0,
    )

    n_pass = sum(1 for c in CHECKS if c["verdict"] == "PASS")
    out = {
        "task": "M5.32 R23-1b adversarial audit (claims D1 to D9)",
        "independence": "no import from the record; energy, gradient, M_00 solve, seed, "
        "half-energy reader, degree reader, disclination counter and fit forms reused from the "
        "earlier independent auditor (m5_32_r23_1_audit.py); shells, r_eps reader and variants, "
        "rms radius, linear theory, k4, box integral, stability probe, label replica are new",
        "checks": CHECKS,
        "counts": {"PASS": n_pass, "FAIL": len(CHECKS) - n_pass, "total": len(CHECKS)},
        "wall_s": round(time.time() - T0, 1),
    }
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=1, default=float)
    log(f"PASS {n_pass} FAIL {len(CHECKS) - n_pass} -> {os.path.relpath(OUT_JSON, HERE)}")
    for c in CHECKS:
        if c["verdict"] == "FAIL":
            log(f"FAIL {c['id']}")


if __name__ == "__main__":
    main()
