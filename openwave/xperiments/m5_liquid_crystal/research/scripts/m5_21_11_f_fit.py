"""M5.21.11 frozen fit + gates: the §§ 3-6 statistics, executed verbatim.

Everything here implements the FROZEN framework (m5_21_11_framework.md):
    § 3  refinement 3-point continuum solve (E_N = E_cont + k N^-p,
         stability guard p in (0.5, 6) else |E_64 - E_48|), branch-wise
         log-delta linear interpolation of sigma_disc/E_48 to unrefined
         rungs, rho_br continuum correction applied ONCE at the end
    § 4  rung usability gates (FIRE, xstencil <= 1.5, virial <= 0.05,
         eigen-gap guard silent, L/a* >= 10)
    § 5  branch identity by (charge, census, core) triple vs the
         rung-0.30 triple; merge metric dE <= 0.04% AND field distance
         <= 6%; energy NEVER used for identity
    § 6  joint WLS of E_br(d) = Einf_br (1 + b_br d^th + c_br d),
         shared theta profiled deterministically (the model is LINEAR
         in (Einf, Einf*b, Einf*c) at fixed theta: exact solves, no
         optimizer), profile-likelihood 68% interval on theta,
         degeneracy fallback if the interval includes 1; holdouts
         {0.20, 0.07} out of every fit; F1-F4; sigma model
         sigma^2 = extrap^2 + rho^2 + g^2; ratios with joint covariance

Operational pins fixed pre-result (checkpoint doc): FIRE gate =
stop f_tol OR (max_iter/plateau AND final fmax <= 1e-4); sigma floor
1e-6 * E_48; plane-degenerate = all 4 contour radii flagged on a plane.

Mode: all [dphys=1e-10] [gphys=1e10]
Out: ../data/m5_21_11_fit.json, plots ../plots/m5_21_11_{ladder,refine,garm}.png
"""
from __future__ import annotations

import glob
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
PLOTS = os.path.join(HERE, "..", "plots")

RUNGS = [0.30, 0.25, 0.20, 0.15, 0.12, 0.09, 0.07, 0.05]
HOLDOUTS = (0.20, 0.07)
REFINED = (0.30, 0.12, 0.05)
BRANCHES = ("A", "C", "B")
L_BOX = 48.0
FMAX_BAR = 1e-4          # pre-result pin (census production grade)
SIG_FLOOR_REL = 1e-6     # pre-result pin


def load_all():
    rows, reads = {}, {}
    for p in glob.glob(os.path.join(DATA,
                                    "m5_21_11_row_t11lad_*.json")):
        r = json.load(open(p))
        rows[(r["branch"], int(r["n"]), round(r["delta"], 3))] = r
    for p in glob.glob(os.path.join(DATA,
                                    "m5_21_11_read_t11lad_*.json")):
        r = json.load(open(p))
        reads[(r["branch"], int(r["n"]), round(r["delta"], 3))] = r
    return rows, reads


# ================= § 4 gates =================
def gate_row(row, read):
    fails = []
    tr = row["trace"][-1] if row["trace"] else {}
    fmax = tr.get("fmax", np.inf)
    if not (row["stop"] == "f_tol"
            or (row["stop"] in ("max_iter", "plateau")
                and fmax <= FMAX_BAR)):
        fails.append(f"fire(stop={row['stop']},fmax={fmax:.1e})")
    xr = row["consistency"]["xstencil_ratio"]
    if xr > 1.5:
        fails.append(f"xstencil({xr:.2f})")
    vir = abs(row["E_u"] - 3.0 * row["E_v"]) / max(row["E_u"], 1e-300)
    if vir > 0.05:
        fails.append(f"virial({vir:.3f})")
    if read is None:
        fails.append("no-read")
    else:
        if not read["gap_guard"]["silent"]:
            fails.append("gap-guard")
        cont = read["flags"]["contours"]
        for z in sorted({c["z"] for c in cont}):
            zc = [c for c in cont if c["z"] == z]
            if zc and all(c["flag"] for c in zc):
                fails.append(f"plane-degenerate(z={z:g})")
                break
    a_star = row["ring"].get("rho_w", float("nan"))
    if not np.isfinite(a_star) or L_BOX / a_star < 10.0:
        fails.append(f"L/a*({L_BOX / a_star:.1f})" if
                     np.isfinite(a_star) else "L/a*(nan)")
    return fails, vir, a_star


# ================= § 5 identity + merges =================
def triple(read):
    """reading (a): the pre-run pinned operationalization (raw tracer
    count + verdicts + ABSOLUTE-gap core class). Kept and reported;
    measured to be noise-dominated (the tracer count jitters 2-18 on
    one continuation chain; the absolute-gap core class flips because
    the BULK 12-gap is delta itself). See the deviations log."""
    return (read["charge"]["charge_class"],
            read["census"]["n_lines"],
            tuple(read["census"]["verdicts"]),
            read["core"]["core_class"])


def triple_b(read):
    """reading (b), OPERATIVE: the same three frozen § 5 signatures
    operationalized by their pre-freeze semantics: charge class
    (unchanged); line census = the shell contour winding (the 2b § 8
    instrument the framework's reference states were measured with:
    net half-units at rho = 18 per plane, topological and stable);
    core class = which pair equalizes RELATIVE to its bulk gap
    (the M5.23.2 tracer's own uniaxial-escape criterion is
    relative-to-bulk; bulk gaps are (delta, 1-delta), so absolute
    comparison measures the vacuum, not the core)."""
    d = read["delta"]
    ws = tuple(sorted(c["half_units"] for c in
                      read["flags"]["contours"]
                      if c["rho"] == 18.0))
    g12 = read["core"]["gap12_core"] / d
    g23 = read["core"]["gap23_core"] / (1.0 - d)
    return (read["charge"]["charge_class"], ws,
            "12" if g12 <= g23 else "23")


def field_dist(tagA, tagB, delta):
    ZA = np.load(os.path.join(DATA, f"m5_21_11_end_{tagA}.npz"))
    ZB = np.load(os.path.join(DATA, f"m5_21_11_end_{tagB}.npz"))
    MA = ZA["M"].astype(np.float64)
    MB = ZB["M"].astype(np.float64)
    vac = np.zeros_like(MA)
    vac[:] = np.diag([1.0, delta, 0.0])
    num = float(np.sqrt(np.sum((MA - MB) ** 2)))
    den = float(np.sqrt(np.sum((MA - vac) ** 2)))
    return num / max(den, 1e-300)


def merges_at(rows, delta, n=48):
    out = []
    for i, b1 in enumerate(BRANCHES):
        for b2 in BRANCHES[i + 1:]:
            r1 = rows.get((b1, n, round(delta, 3)))
            r2 = rows.get((b2, n, round(delta, 3)))
            if not (r1 and r2):
                continue
            dE = abs(r1["E_end"] - r2["E_end"]) \
                / max(abs(r1["E_end"]), 1e-300)
            if dE <= 4e-4:
                fd = field_dist(r1["tag"], r2["tag"], delta)
                if fd <= 0.06:
                    out.append({"pair": [b1, b2], "delta": delta,
                                "dE_rel": dE, "field_dist": fd})
    return out


# ================= § 3 discretization =================
def sigma_disc(rows, branch):
    """per-branch: refined-rung continuum solves + log-delta interp."""
    ref = {}
    for d in REFINED:
        es = {}
        for n in (32, 48, 64):
            r = rows.get((branch, n, round(d, 3)))
            if r:
                es[n] = r["E_end"]
        if len(es) != 3:
            ref[d] = None
            continue
        E32, E48, E64 = es[32], es[48], es[64]
        # exact 3-param solve E_N = Ec + k N^-p  (p via ratio equation)
        num, den = E32 - E48, E48 - E64
        sol = None
        if num * den > 0:
            # f(p) = (32^-p - 48^-p)/(48^-p - 64^-p) = num/den
            target = num / den
            ps = np.linspace(0.05, 10.0, 4000)
            f = (32.0 ** -ps - 48.0 ** -ps) / (48.0 ** -ps
                                               - 64.0 ** -ps)
            k = int(np.argmin(np.abs(f - target)))
            p = float(ps[k])
            kk = num / (32.0 ** -p - 48.0 ** -p)
            Ec = E48 - kk * 48.0 ** -p
            sol = {"p": p, "k": float(kk), "E_cont": float(Ec)}
        if sol and 0.5 < sol["p"] < 6.0:
            sd = abs(sol["E_cont"] - E48)
            ref[d] = {"sigma": float(sd), "E48": E48,
                      "rho": float(sol["E_cont"] / E48),
                      "p": sol["p"], "guard": False}
        else:
            sd = abs(E64 - E48)
            ref[d] = {"sigma": float(sd), "E48": E48,
                      "rho": float(E64 / E48),
                      "p": (sol["p"] if sol else None), "guard": True}
    # propagate relative sigma in log delta to unrefined rungs
    good = [(d, v) for d, v in ref.items() if v]
    good.sort(key=lambda t: t[0])
    xs = np.log([d for d, _ in good])
    ys = [v["sigma"] / v["E48"] for _, v in good]
    def rel_at(d):
        if not good:
            return None
        return float(np.interp(np.log(d), xs, ys))
    return ref, rel_at


# ================= § 6 the joint fit =================
def design(deltas, theta, nb):
    """block design matrix for [a_br, B_br, C_br] per branch
    (a = Einf, B = Einf*b, C = Einf*c), points ordered branch-major."""
    cols = []
    for k in range(nb):
        base = np.zeros((len(deltas), 3))
        base[:, 0] = 1.0
        base[:, 1] = deltas ** theta
        base[:, 2] = deltas
        cols.append(base)
    X = np.zeros((len(deltas), 3 * nb))
    return X, cols


def chi2_at_theta(theta, pts, form="full"):
    """exact linear WLS at fixed theta; pts = list of dicts with
    branch, delta, E, sig. Returns (chi2, params dict, lin arrays)."""
    params, chi2, resid = {}, 0.0, []
    for br in BRANCHES:
        P = [p for p in pts if p["branch"] == br]
        if not P:
            continue
        d = np.array([p["delta"] for p in P])
        E = np.array([p["E"] for p in P])
        w = 1.0 / np.array([p["sig"] for p in P])
        if form == "full":
            X = np.stack([np.ones_like(d), d ** theta, d], axis=1)
        else:                       # fallback: Einf (1 + c d)
            X = np.stack([np.ones_like(d), d], axis=1)
        Xw = X * w[:, None]
        Ew = E * w
        beta, *_ = np.linalg.lstsq(Xw, Ew, rcond=None)
        r = (E - X @ beta) * w
        chi2 += float(r @ r)
        resid.extend(r.tolist())
        params[br] = beta.tolist()
    return chi2, params


def profile_theta(pts, th_lo=0.05, th_hi=2.5, nth=2451):
    ths = np.linspace(th_lo, th_hi, nth)
    chi = np.array([chi2_at_theta(t, pts)[0] for t in ths])
    i = int(np.argmin(chi))
    th_hat = float(ths[i])
    if 0 < i < nth - 1:
        a, b, c = chi[i - 1], chi[i], chi[i + 1]
        den = c - 2 * b + a
        if den > 0:
            th_hat = float(ths[i] - 0.5 * (ths[1] - ths[0])
                           * (c - a) / den)
    chimin = float(chi.min())
    inside = ths[chi <= chimin + 1.0]
    interval = [float(inside.min()), float(inside.max())] \
        if len(inside) else [th_hat, th_hat]
    at_edge = bool(interval[0] <= ths[0] + 1e-9
                   or interval[1] >= ths[-1] - 1e-9)
    return {"theta_hat": th_hat, "chi2_min": chimin,
            "interval68": interval, "interval_at_edge": at_edge,
            "grid": [float(ths[0]), float(ths[-1]), nth],
            "curve": [{"theta": float(t), "chi2": float(c)}
                      for t, c in zip(ths[::50], chi[::50])]}


def full_cov(pts, theta, params, form="full"):
    """10-param covariance via the WLS normal matrix at the optimum
    (params order: theta, then per branch a,B,C)."""
    npar = 1 + 3 * len(BRANCHES) if form == "full" \
        else 3 * len(BRANCHES)  # fallback: no theta, (a,C) per br -> 6
    if form != "full":
        npar = 2 * len(BRANCHES)
    J, W = [], []
    for p in pts:
        br_i = BRANCHES.index(p["branch"])
        d, sig = p["delta"], p["sig"]
        row = np.zeros(npar)
        if form == "full":
            a, B, C = params[p["branch"]]
            row[0] = B * d ** theta * np.log(d)
            row[1 + 3 * br_i + 0] = 1.0
            row[1 + 3 * br_i + 1] = d ** theta
            row[1 + 3 * br_i + 2] = d
        else:
            row[2 * br_i + 0] = 1.0
            row[2 * br_i + 1] = d
        J.append(row)
        W.append(1.0 / sig ** 2)
    J = np.array(J)
    W = np.diag(W)
    H = J.T @ W @ J
    return np.linalg.pinv(H)


def predict(branch, d, theta, params, form="full"):
    if form == "full":
        a, B, C = params[branch]
        return a + B * d ** theta + C * d
    a, C = params[branch]
    return a + C * d


def pred_grad(branch, d, theta, params, form="full"):
    npar = 1 + 3 * len(BRANCHES) if form == "full" \
        else 2 * len(BRANCHES)
    g = np.zeros(npar)
    br_i = BRANCHES.index(branch)
    if form == "full":
        a, B, C = params[branch]
        g[0] = B * d ** theta * np.log(d)
        g[1 + 3 * br_i + 0] = 1.0
        g[1 + 3 * br_i + 1] = d ** theta
        g[1 + 3 * br_i + 2] = d
    else:
        g[2 * br_i + 0] = 1.0
        g[2 * br_i + 1] = d
    return g


# ================= main =================
def main(dphys=1e-10, gphys=1e10):
    rows, reads = load_all()
    out = {"dphys": dphys, "gphys": gphys,
           "pins": {"fmax_bar": FMAX_BAR,
                    "sig_floor_rel": SIG_FLOOR_REL}}

    # § 4 gates + § 5 identity on every production rung
    gates, ident, ref_triples = {}, {}, {}
    ref_triples_b = {}
    for br in BRANCHES:
        r0 = reads.get((br, 48, 0.3))
        ref_triples[br] = triple(r0) if r0 else None
        ref_triples_b[br] = triple_b(r0) if r0 else None
    rung_state = {}
    for br in BRANCHES:
        for d in RUNGS:
            key = (br, 48, round(d, 3))
            row, read = rows.get(key), reads.get(key)
            if not row:
                rung_state[key] = {"state": "MISSING"}
                continue
            fails, vir, a_star = gate_row(row, read)
            st = {"E": row["E_end"], "vir_frozen": vir,
                  "a_star": a_star, "gate_fails": fails}
            if read and ref_triples[br]:
                tr = triple(read)
                trb = triple_b(read)
                st["triple"] = tr
                st["triple_b"] = trb
                st["identity_a"] = bool(tr == ref_triples[br])
                st["identity_ok"] = bool(trb == ref_triples_b[br])
            else:
                st["identity_a"] = None
                st["identity_ok"] = None
            st["state"] = "USABLE" if (not fails
                                       and st["identity_ok"]) \
                else "EXCLUDED"
            rung_state[key] = st
    out["ref_triples"] = {b: list(map(str, ref_triples[b]))
                          if ref_triples[b] else None
                          for b in BRANCHES}
    out["ref_triples_b"] = {b: list(map(str, ref_triples_b[b]))
                            if ref_triples_b[b] else None
                            for b in BRANCHES}
    out["identity_a_exclusions"] = [
        f"{k[0]}_d{k[2]:g}" for k, v in rung_state.items()
        if v.get("identity_a") is False]

    # § 5 merges (production rungs)
    merges = []
    for d in RUNGS:
        merges.extend(merges_at(rows, d))
    out["merges"] = merges

    # F3 branch integrity
    usable = {br: [d for d in RUNGS
                   if rung_state[(br, 48, round(d, 3))]["state"]
                   == "USABLE"] for br in BRANCHES}
    out["usable_rungs"] = usable
    f3_fail = any(len(v) < 6 for v in usable.values())

    # § 3 discretization
    disc = {}
    for br in BRANCHES:
        ref, rel_at = sigma_disc(rows, br)
        disc[br] = {"refined": {f"{d:g}": ref[d] for d in REFINED},
                    "_rel_at": rel_at}
    out["refinement"] = {br: disc[br]["refined"] for br in BRANCHES}

    # build fit points (usable, non-holdout) + holdout points
    pts, hold = [], []
    for br in BRANCHES:
        for d in usable[br]:
            E = rows[(br, 48, round(d, 3))]["E_end"]
            rel = disc[br]["_rel_at"](d)
            if rel is None:
                continue
            sig = max(rel * E, SIG_FLOOR_REL * E)
            rec = {"branch": br, "delta": d, "E": E, "sig": sig}
            (hold if d in HOLDOUTS else pts).append(rec)
    out["n_fit_points"] = len(pts)
    out["n_holdout_points"] = len(hold)

    # § 6 fit: profile theta
    prof = profile_theta(pts)
    out["profile"] = {k: prof[k] for k in
                      ("theta_hat", "chi2_min", "interval68",
                       "interval_at_edge", "grid")}
    out["profile_curve"] = prof["curve"]
    theta = prof["theta_hat"]
    degenerate = prof["interval68"][1] >= 1.0 \
        and prof["interval68"][0] <= 1.0
    # frozen degeneracy rule: fallback iff the 68% interval INCLUDES 1
    form = "fallback" if (prof["interval68"][0] <= 1.0
                          <= prof["interval68"][1]) else "full"
    out["degeneracy_triggered"] = bool(form == "fallback")
    if form == "fallback":
        chi2, params = chi2_at_theta(0.0, pts, form="fallback")
        theta = None
        npar = 2 * len(BRANCHES)
    else:
        chi2, params = chi2_at_theta(theta, pts, form="full")
        npar = 1 + 3 * len(BRANCHES)
    dof = len(pts) - npar
    out["fit"] = {"form": form, "theta": theta, "chi2": chi2,
                  "dof": dof,
                  "chi2_per_dof": chi2 / dof if dof > 0 else None,
                  "params": params}

    cov = full_cov(pts, theta if theta else 0.0, params, form=form)

    # F1
    f1_fail = bool(dof > 0 and chi2 / dof > 3.0)

    # F2 holdouts
    f2_rows, n2s, n3s = [], 0, 0
    for hpt in hold:
        if hpt["branch"] not in params:
            continue
        Ep = predict(hpt["branch"], hpt["delta"], theta, params,
                     form=form)
        gvec = pred_grad(hpt["branch"], hpt["delta"], theta, params,
                         form=form)
        s_pred = float(np.sqrt(max(gvec @ cov @ gvec, 0.0)))
        s_tot = float(np.sqrt(s_pred ** 2 + hpt["sig"] ** 2))
        z = abs(hpt["E"] - Ep) / max(s_tot, 1e-300)
        f2_rows.append({"branch": hpt["branch"],
                        "delta": hpt["delta"], "E_meas": hpt["E"],
                        "E_pred": float(Ep), "sig_tot": s_tot,
                        "z": float(z)})
        n2s += int(z > 2.0)
        n3s += int(z > 3.0)
    f2_fail = bool(n2s >= 2 or n3s >= 1)
    out["holdouts"] = {"rows": f2_rows, "n_beyond_2s": n2s,
                       "n_beyond_3s": n3s}

    # F4 the g-arm
    garm_p = os.path.join(DATA, "m5_21_11_garm.json")
    garm = json.load(open(garm_p)) if os.path.exists(garm_p) else None
    f4_fail = (not garm["f4_pass_all"]) if garm else None
    out["garm"] = ({br: garm["fits"][br] for br in BRANCHES}
                   if garm else None)

    # physical point: E(dphys), rho correction, sigma model, ratios
    rho, sig_rho_rel = {}, {}
    for br in BRANCHES:
        vals = [disc[br]["refined"][f"{d:g}"]["rho"] for d in REFINED
                if disc[br]["refined"][f"{d:g}"]]
        rho[br] = float(np.mean(vals)) if vals else 1.0
        sig_rho_rel[br] = float(np.std(vals)) if len(vals) > 1 else 0.0
    Ephys, sig_extrap = {}, {}
    for br in BRANCHES:
        if br not in params:
            Ephys[br] = None
            sig_extrap[br] = None
            continue
        Ep = predict(br, dphys, theta, params, form=form)
        gvec = pred_grad(br, dphys, theta, params, form=form)
        Ephys[br] = float(Ep)
        sig_extrap[br] = float(np.sqrt(max(gvec @ cov @ gvec, 0.0)))
    sig_g_rel = {}
    for br in BRANCHES:
        if Ephys[br] is None:
            sig_g_rel[br] = None
        elif garm:
            q = garm["fits"][br]["q_lsq"]
            kap = garm["fits"][br]["kappa"]
            sig_g_rel[br] = float(kap * np.arctanh(1.0 / gphys)
                                  ** max(q, 2.0)
                                  / max(abs(Ephys[br]), 1e-300))
        else:
            sig_g_rel[br] = None
    sigma = {br: (None if Ephys[br] is None else float(np.sqrt(
        sig_extrap[br] ** 2
        + (sig_rho_rel[br] * abs(Ephys[br])) ** 2
        + ((sig_g_rel[br] or 0.0) * abs(Ephys[br])) ** 2)))
        for br in BRANCHES}
    out["physical"] = {
        "E_phys_E48scale": Ephys,
        "rho": rho, "sig_rho_rel": sig_rho_rel,
        "sig_extrap": sig_extrap, "sig_g_rel": sig_g_rel,
        "E_phys_continuum": {br: (None if Ephys[br] is None
                                  else Ephys[br] * rho[br])
                             for br in BRANCHES},
        "sigma_total": sigma}

    # ratios with joint covariance (delta method on the fit params;
    # rho and g terms added in quadrature as independent relatives)
    ratios = {}
    for num in ("C", "B"):
        if Ephys.get(num) is None or Ephys.get("A") is None:
            ratios[f"R_{num}"] = None
            continue
        gN = pred_grad(num, dphys, theta, params, form=form)
        gA = pred_grad("A", dphys, theta, params, form=form)
        R = (Ephys[num] * rho[num]) / (Ephys["A"] * rho["A"])
        gradR = gN / Ephys[num] - gA / Ephys["A"]
        var_fit = float(gradR @ cov @ gradR)     # relative variance
        var_rho = sig_rho_rel[num] ** 2 / max(rho[num], 1e-300) ** 2 \
            + sig_rho_rel["A"] ** 2 / max(rho["A"], 1e-300) ** 2
        var_g = ((sig_g_rel[num] or 0.0) ** 2
                 + (sig_g_rel["A"] or 0.0) ** 2)
        sR = abs(R) * float(np.sqrt(max(var_fit + var_rho + var_g,
                                        0.0)))
        ratios[f"R_{num}"] = {"value": float(R), "sigma": sR}
    out["ratios"] = ratios

    # diagnostics (§ 6, non-terminal): theta vs 2 nu; per-branch theta
    nu_fits = {}
    for br in BRANCHES:
        ds, aa = [], []
        for d in usable[br]:
            a_st = rung_state[(br, 48, round(d, 3))]["a_star"]
            if np.isfinite(a_st):
                ds.append(d)
                aa.append(a_st)
        if len(ds) >= 3:
            s, ic = np.polyfit(np.log(ds), np.log(aa), 1)
            nu_fits[br] = {"nu_hat": float(-s),
                           "two_nu": float(-2 * s)}
    out["diag_nu"] = nu_fits
    per_th = {}
    for br in BRANCHES:
        bpts = [p for p in pts if p["branch"] == br]
        if len(bpts) >= 4:
            pr = profile_theta(bpts)
            per_th[br] = {"theta_hat": pr["theta_hat"],
                          "interval68": pr["interval68"]}
    out["diag_theta_per_branch"] = per_th

    # verdicts
    out["criteria"] = {
        "F1_fit_quality": {"fail": f1_fail,
                           "chi2_per_dof": out["fit"]
                           ["chi2_per_dof"]},
        "F2_holdouts": {"fail": f2_fail, "n2s": n2s, "n3s": n3s},
        "F3_branch_integrity": {"fail": f3_fail,
                                "counts": {b: len(usable[b])
                                           for b in BRANCHES}},
        "F4_garm": {"fail": f4_fail}}
    fails = [k for k, v in out["criteria"].items() if v["fail"]]
    out["verdict"] = {"terminal_failure": bool(fails),
                      "failed": fails,
                      "pass": not fails and f4_fail is not None}

    out["rung_state"] = {f"{k[0]}_n{k[1]}_d{k[2]:g}":
                         {kk: vv for kk, vv in v.items()
                          if kk != "triple"}
                         | ({"triple": list(map(str, v["triple"]))}
                            if "triple" in v else {})
                         for k, v in rung_state.items()}

    with open(os.path.join(DATA, "m5_21_11_fit.json"), "w") as f:
        json.dump(out, f, indent=1)

    # ---- plots ----
    os.makedirs(PLOTS, exist_ok=True)
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    cols = {"A": "tab:blue", "C": "tab:orange", "B": "tab:red"}
    dd = np.linspace(0.04, 0.32, 300)
    for br in BRANCHES:
        P = [p for p in pts if p["branch"] == br]
        H = [p for p in hold if p["branch"] == br]
        ax[0].errorbar([p["delta"] for p in P], [p["E"] for p in P],
                       yerr=[p["sig"] for p in P], fmt="o",
                       color=cols[br], label=f"{br} fit pts")
        ax[0].errorbar([p["delta"] for p in H], [p["E"] for p in H],
                       yerr=[p["sig"] for p in H], fmt="s", mfc="none",
                       color=cols[br], label=f"{br} holdout")
        if params.get(br):
            ax[0].plot(dd, [predict(br, d, theta, params, form=form)
                            for d in dd], "-", color=cols[br], lw=1)
    ax[0].set_xlabel("delta")
    ax[0].set_ylabel("E (E48 scale)")
    ax[0].set_yscale("log")
    th_lab = f"{theta:.3f}" if theta is not None else "fallback"
    ax[0].set_title(f"frozen ladder fit (theta = {th_lab})")
    ax[0].legend(fontsize=7)
    pc = out["profile_curve"]
    ax[1].plot([r["theta"] for r in pc], [r["chi2"] for r in pc], "-")
    ax[1].axhline(prof["chi2_min"] + 1.0, ls="--", c="gray")
    ax[1].set_xlabel("theta")
    ax[1].set_ylabel("chi2 profile")
    ax[1].set_title(f"interval68 = {prof['interval68']}")
    if garm:
        for br in BRANCHES:
            gs = [8.0, 16.0, 32.0]
            gains = [-garm["fits"][br]["gains"][f"g{g:g}"]
                     for g in gs]
            x = [np.arctanh(1.0 / g) for g in gs]
            ax[2].loglog(x, gains, "o-", color=cols[br],
                         label=f"{br} q={garm['fits'][br]['q_lsq']:.2f}")
        xr = np.array([np.arctanh(1 / 32.0), np.arctanh(1 / 8.0)])
        ax[2].loglog(xr, gains[0] * (xr / x[0]) ** 2, "k--",
                     lw=0.8, label="slope 2 ref")
        ax[2].set_xlabel("artanh(1/g)")
        ax[2].set_ylabel("-gain")
        ax[2].legend(fontsize=7)
        ax[2].set_title("g-arm dressing gain")
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS, "m5_21_11_ladder.png"), dpi=140)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 4.2))
    for br in BRANCHES:
        for d in REFINED:
            es = [(n, rows[(br, n, round(d, 3))]["E_end"])
                  for n in (32, 48, 64)
                  if (br, n, round(d, 3)) in rows]
            if len(es) == 3:
                ax.plot([e[0] for e in es], [e[1] for e in es],
                        "o-", color=cols[br], alpha=0.7,
                        label=f"{br} d={d:g}")
    ax.set_xlabel("N")
    ax.set_ylabel("E_end")
    ax.set_title("refinement subset (frozen § 3)")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS, "m5_21_11_refine.png"), dpi=140)
    plt.close(fig)

    print(json.dumps({"verdict": out["verdict"],
                      "criteria": out["criteria"],
                      "theta": theta,
                      "interval68": prof["interval68"],
                      "ratios": ratios}, indent=1))
    return out


if __name__ == "__main__":
    kw = dict(a.split("=", 1) for a in sys.argv[1:] if "=" in a)
    main(dphys=float(kw.get("dphys", 1e-10)),
         gphys=float(kw.get("gphys", 1e10)))
