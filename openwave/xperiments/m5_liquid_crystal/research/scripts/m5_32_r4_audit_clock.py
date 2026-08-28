"""M5.32 R4 (the clock) ADVERSARIAL AUDIT: an independent rebuild of the
fixed-J electron claims of the candidate

    L_lambda = -4 [(1 - lambda) I1 + lambda I1_h] - V4,
    I1_h = I1 with eta -> h_cov = eta + 2 (eta u)(eta u)^T,

at the toy point g = 32, delta = 0.3, s = -1. The producer's script
(m5_32_r4_clock.py) and its JSON were NOT read. Oracles: the registries
(m5_32_lagrangian.py: I1, V4; m5_32_terms_ext.py: h_cov_np), the
certified stack m5_21_3_a_4d.py (sym stencil, e_parts, kin_of), the
family builders m5_21_8_b_lattice.py (dressed hedgehog, a0_unit) and
the audit's own R2 instrument m5_32_r2_audit_lattice.py (qb_field,
lattice_parts, h_from_Q, the rigid Family), reused as-is.

EQUATIONS FIRST
---------------
Jets A_mu = d_mu M, curvature F_{mu nu} = A_mu eta A_nu - A_nu eta A_mu.
    q_eta(F) = tr(eta F eta F^T),  q_h(F) = tr(h F h F^T)
    q_lam    = (1 - lam) q_eta + lam q_h                (LINEAR in lam)
Legendre energy of L_lambda with A_0 = omega a0:
    E(omega) = E_stat + omega^2 kin
    E_stat   = h^3 sum_x [ 4 sum_{i<j} q_lam(F_ij) ] + E_V
    kin      = h^3 sum_x   4 sum_i   q_lam(comm_eta(a0, A_i))
Fixed-J electron (Legendre in omega at fixed J = dE/domega = 2 omega kin):
    E_J(amp, R) = E_stat(amp, R) + J^2 / (4 kin(amp, R)),
    omega* = J / (2 kin*),   R* = argmin_R min_amp E_J.
Boost-dressing families about the time axis on the m5_21_8 hedgehog Mb
(m = 0, rotation-only; eigenvalues pinned so E_V = 0 on every member):
    M(amp, R) = Qb Mb Qb^T,  a0 = Qb a0_unit Qb^T,  h = Qb^-T Qb^-1,
    Qb = radial boost with rapidity b(r) = amp tanh(r/2) w(r; R)
    w = exp(-(r/R)^p) (p = 1, 2)                 the producer's family
    w = 1 (r < R), cos^2(pi (r - R) / 4) ramp to 0 at R + 2   hard support
    w = (1 + r/R)^-q, q = 1, 2, 3               power-law tail
    w = 1                                        rigid (the R2 instrument)
K2 decomposition at fixed amp = amp*(R) (envelope theorem: the partial
derivative at amp* is the total derivative of the minimized E_J):
    dE_J/dR = dE_stat/dR + d(J^2 / 4 kin)/dR    (central differences)
Dilation: b_mu(r) = amp tanh(r / (2 mu)) w(r / mu; R), dE_J/dmu at mu = 1.
K5 plane-wave probe. Single-phase profile M(phase), phase = k.x - w t:
    A_mu = k_mu M'(phase)  =>  F_{mu nu} = k_mu k_nu [M', M']_eta = 0
identically: I1 and I1_h vanish for ANY single plane wave, any amplitude
(checked numerically with exact Frechet jets). Crossed waves M = e^W vac
e^W^T, W = eps [G_a cos(k1.x - w1 t) + G_b cos(k2.x - w2 t)], k1 != k2:
F = O(eps^2), density O(eps^4), e(2 eps) / e(eps) -> 16 (quartic
kernel). Non-Lorentz control M = vac + eps D cos(phase): V4 quadratic
(ratio 4, k-independent mass-like), the I1 part still zero / quartic.

Modes (each writes ../data/m5_32_r4_audit_clock_<mode>.json; `merge`
assembles ../data/m5_32_r4_audit_clock.json):
    box n=32 L=48 kinds=exp2,exp1   (K1/K3/K4 grid; any n, L, kinds)
    alt                              (K2 alternative families, n32 L48)
    decomp                           (K2 derivative decomposition + dilation)
    pw                               (K5 plane-wave probe)
    merge
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys
import time

import numpy as np
from scipy.linalg import expm, expm_frechet

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT_STEM = "m5_32_r4_audit_clock"


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


R2 = _load("m5_32_r2_audit_lattice", "m5_32_r2_audit_lattice.py")
L0 = R2.L0
TX = R2.TX
B3 = R2.B3
B8 = R2.B8
ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
G, DELTA, S = 32.0, 0.3, -1.0
T0 = time.time()

AMPS = np.round(np.concatenate([np.arange(0.0, 0.0601, 0.0025),
                                [0.07, 0.08, 0.1, 0.12, 0.15, 0.2, 0.3]]), 6)
LAMS = (0.0, 0.75, 1.0)
JS = (5.0, 10.0, 50.0, 200.0, 800.0)


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


# ================= the localized dressing families =================
def window(r, R, kind):
    if kind == "rigid":
        return np.ones_like(r)
    if kind == "exp2":
        return np.exp(-(r / R) ** 2)
    if kind == "exp1":
        return np.exp(-(r / R))
    if kind == "hard":
        x = np.clip((r - R) / 2.0, 0.0, 1.0)
        return np.cos(0.5 * np.pi * x) ** 2
    if kind.startswith("pow"):
        q = float(kind[3:])
        return (1.0 + r / R) ** (-q)
    raise ValueError(kind)


class LocFamily:
    def __init__(self, n, L):
        self.n, self.L = n, L
        self.cfg = B3.base_cfg(s=S, g=G, n=n, L=L, delta=DELTA)
        self.Mb = B8.dressed(self.cfg, 0.0)
        self.a0 = B8.a0_unit(self.cfg, 0.0)
        X, Y, Z = B3.coords(n, self.cfg["h"])
        self.R = np.sqrt(X * X + Y * Y + Z * Z)
        self.cache = {}

    def b_of(self, amp, R, kind, mu=1.0):
        r = self.R / mu
        return amp * np.tanh(r / 2.0) * window(r, R, kind)

    def parts(self, amp, R, kind, mu=1.0, registry_u=False):
        """(E_u0, E_u1, E_V, kin0, kin1), h^3-weighted."""
        key = (round(amp, 12), round(R, 9), kind, round(mu, 9), registry_u)
        if key in self.cache:
            return self.cache[key]
        Qb = R2.qb_field(self.cfg, self.b_of(amp, R, kind, mu))
        Md = B3.sym4(R2.conj(Qb, self.Mb))
        a0d = B3.sym4(R2.conj(Qb, self.a0))
        h = R2.h_of(Md) if registry_u else R2.h_from_Q(Qb)
        out = tuple(float(v) for v in R2.lattice_parts(Md, self.cfg, a0d, h=h))
        self.cache[key] = out
        return out


def es_kin(p, lam):
    eu0, eu1, ev, k0, k1 = p
    return R2.mix(eu0, eu1, lam) + ev, R2.mix(k0, k1, lam)


def min_over_amp(rows, J, lam, kin_guard=0.0):
    """rows: list of (amp, parts). E_J on the amp grid, parabolic refinement
    around an interior argmin; kin <= kin_guard -> E = inf (excluded)."""
    amps = np.array([a for a, _ in rows])
    es = np.array([es_kin(p, lam)[0] for _, p in rows])
    kk = np.array([es_kin(p, lam)[1] for _, p in rows])
    ok = kk > kin_guard
    E = np.where(ok, es + J * J / (4.0 * np.where(ok, kk, 1.0)), np.inf)
    i = int(np.argmin(E))
    if not np.isfinite(E[i]):
        return None
    a_star, E_star = amps[i], E[i]
    interior = 0 < i < len(rows) - 1 and np.isfinite(E[i - 1]) and np.isfinite(E[i + 1])
    refined = False
    if interior:
        a, b, c = E[i - 1], E[i], E[i + 1]
        da = amps[i + 1] - amps[i]
        if c - 2 * b + a > 1e-14 and abs(amps[i] - amps[i - 1] - da) < 1e-12:
            a_star = amps[i] - 0.5 * da * (c - a) / (c - 2 * b + a)
            E_star = b - 0.125 * (c - a) ** 2 / (c - 2 * b + a)
            refined = True
    k_star = float(np.interp(a_star, amps, kk))
    es_star = float(np.interp(a_star, amps, es))
    edge = "interior" if interior else ("amp_min" if i == 0 else "amp_max")
    if ok.sum() and i == int(np.where(ok)[0][-1]) and not ok.all():
        edge = "kin_guard_edge"
    return {"amp_star": float(a_star), "E_J": float(E_star), "E_stat": es_star,
            "kin": k_star, "J2_over_4kin": float(J * J / (4.0 * k_star)),
            "omega_star": float(J / (2.0 * k_star)), "grid_index": i,
            "amp_edge": edge, "parabolic_refined": refined,
            "n_kin_positive": int(ok.sum()), "n_amps": len(rows)}


def interior_test(Rs, EJ):
    """argmin over the R grid; interior if strictly inside."""
    EJ = np.array([np.inf if e is None else e for e in EJ])
    i = int(np.argmin(EJ))
    mono = bool(np.all(np.diff(EJ[np.isfinite(EJ)]) < 0))
    return {"R_star": float(Rs[i]), "argmin_index": i, "interior": bool(0 < i < len(Rs) - 1),
            "at_wall": bool(i == len(Rs) - 1), "monotone_decreasing_in_R": mono,
            "E_J_of_R": [None if not np.isfinite(e) else float(e) for e in EJ]}


# ================= K1 / K3 / K4: the (R x amp) grid per box =================
def mode_box(n, L, kinds, tag):
    fam = LocFamily(n, L)
    h = fam.cfg["h"]
    Rs = [6.0, 9.0, 12.0, 18.0, L / 2.0]
    # certified cross-check at amp = 0 (lam = 0 parts == e_parts / kin_of)
    eu, ev = B3.e_parts(fam.Mb, fam.cfg)
    kk = B3.kin_of(fam.Mb, fam.a0, fam.cfg)
    p0 = fam.parts(0.0, 6.0, "exp2")
    cert = {"E_u_cert": float(eu), "E_V_cert": float(ev), "kin_cert": float(kk),
            "E_u_audit_lam0": p0[0], "E_u_audit_lam1": p0[1], "kin_audit_lam0": p0[3],
            "kin_audit_lam1": p0[4], "E_V_audit": p0[2]}
    # analytic-h vs registry-h on a localized member
    pa = fam.parts(0.03, 12.0, "exp2")
    pr = fam.parts(0.03, 12.0, "exp2", registry_u=True)
    ucheck = {"amp": 0.03, "R": 12.0, "kind": "exp2",
              "rel_dev_Eu1": abs(pa[1] - pr[1]) / abs(pr[1]),
              "rel_dev_kin1": abs(pa[4] - pr[4]) / abs(pr[4])}
    out = {"tag": tag, "n": n, "L": L, "h": h, "stencil": fam.cfg["stencil"],
           "amps": AMPS.tolist(), "R_grid": Rs, "certified_crosscheck": cert,
           "u_check": ucheck, "kinds": {}}
    for kind in kinds:
        log(f"{tag} kind {kind}: {len(Rs)} R x {len(AMPS)} amps")
        grid = {}
        for R in Rs:
            grid[R] = [(a, fam.parts(a, R, kind)) for a in AMPS]
            log(f"  R={R:g} done; E_stat(lam1) at amp 0.02: "
                f"{es_kin(grid[R][8][1], 1.0)[0]:.3f} kin {es_kin(grid[R][8][1], 1.0)[1]:.2f}")
        res = {"ladder": {}, "scan": {}}
        for R in Rs:
            res["ladder"][f"R_{R:g}"] = {
                "amp": [a for a, _ in grid[R]],
                "E_stat_lam0": [es_kin(p, 0.0)[0] for _, p in grid[R]],
                "E_stat_lam1": [es_kin(p, 1.0)[0] for _, p in grid[R]],
                "kin_lam0": [es_kin(p, 0.0)[1] for _, p in grid[R]],
                "kin_lam1": [es_kin(p, 1.0)[1] for _, p in grid[R]]}
        for lam in LAMS:
            for J in JS:
                for guard in ((0.0,) if lam > 0 else (0.0, 0.02)):
                    key = f"lam_{lam:g}_J_{J:g}" + ("" if guard == 0 else f"_guard_{guard:g}")
                    byR = {}
                    for R in Rs:
                        byR[f"R_{R:g}"] = min_over_amp(grid[R], J, lam, kin_guard=guard)
                    it = interior_test(Rs, [byR[f"R_{R:g}"]["E_J"] if byR[f"R_{R:g}"] else None for R in Rs])
                    res["scan"][key] = {"byR": byR, "R_test": it}
                    if (lam == 1.0 and J in (200.0, 800.0)) or (lam == 0.0 and J == 200.0 and guard == 0):
                        log(f"  {key}: R* {it['R_star']:g} ({'interior' if it['interior'] else 'wall' if it['at_wall'] else 'R_min'}), "
                            f"mono {it['monotone_decreasing_in_R']}, E_J(R) "
                            + " ".join(f"{e:.2f}" if e is not None else "inf" for e in it["E_J_of_R"])
                            + "  amp* " + " ".join(f"{byR[f'R_{R:g}']['amp_star']:.4f}" if byR[f"R_{R:g}"] else "-" for R in Rs)
                            + "  omega* " + " ".join(f"{byR[f'R_{R:g}']['omega_star']:.4f}" if byR[f"R_{R:g}"] else "-" for R in Rs))
        # pure static (J = 0): min over amp of E_stat per R
        res["scan"]["J_0_static"] = {}
        for lam in LAMS:
            byR = {}
            for R in Rs:
                es = np.array([es_kin(p, lam)[0] for _, p in grid[R]])
                i = int(np.argmin(es))
                byR[f"R_{R:g}"] = {"amp_star": float(AMPS[i]), "E_stat": float(es[i]),
                                   "E_stat_amp0": float(es[0]), "dressing_gain": float(es[i] - es[0])}
            res["scan"]["J_0_static"][f"lam_{lam:g}"] = byR
        out["kinds"][kind] = res
    out["runtime_s"] = round(time.time() - T0, 1)
    with open(os.path.join(DATA, f"{OUT_STEM}_{tag}.json"), "w") as f:
        json.dump(out, f, indent=1)
    return out


# ================= K2: alternative families =================
def mode_alt(n=32, L=48.0, kinds=("hard", "pow1", "pow2", "pow3"), tag="alt_n32_L48"):
    fam = LocFamily(n, L)
    Rs = [3.0, 4.5, 6.0, 9.0, 12.0, 18.0, L / 2.0]
    out = {"tag": tag, "n": n, "L": L, "h": fam.cfg["h"], "amps": AMPS.tolist(),
           "R_grid": Rs, "kinds": {}}
    for kind in kinds:
        log(f"{tag} kind {kind}")
        grid = {R: [(a, fam.parts(a, R, kind)) for a in AMPS] for R in Rs}
        res = {"scan": {}, "ladder_at_amp_0.02": {}}
        for R in Rs:
            p = grid[R][8][1]
            res["ladder_at_amp_0.02"][f"R_{R:g}"] = {"E_stat_lam1": es_kin(p, 1.0)[0],
                                                     "kin_lam1": es_kin(p, 1.0)[1]}
        for lam in (0.75, 1.0):
            for J in (50.0, 200.0, 800.0):
                key = f"lam_{lam:g}_J_{J:g}"
                byR = {f"R_{R:g}": min_over_amp(grid[R], J, lam) for R in Rs}
                it = interior_test(Rs, [byR[f"R_{R:g}"]["E_J"] for R in Rs])
                res["scan"][key] = {"byR": byR, "R_test": it}
                log(f"  {key}: R* {it['R_star']:g} ({'INTERIOR' if it['interior'] else 'wall' if it['at_wall'] else 'R_min'}), "
                    f"mono {it['monotone_decreasing_in_R']}, E_J(R) "
                    + " ".join(f"{e:.2f}" for e in it["E_J_of_R"]))
        out["kinds"][kind] = res
    out["runtime_s"] = round(time.time() - T0, 1)
    with open(os.path.join(DATA, f"{OUT_STEM}_{tag}.json"), "w") as f:
        json.dump(out, f, indent=1)
    return out


# ================= K2: derivative decomposition + dilation =================
def mode_decomp(n=32, L=48.0, tag="decomp_n32_L48"):
    fam = LocFamily(n, L)
    kind, lam = "exp2", 1.0
    out = {"tag": tag, "n": n, "L": L, "h": fam.cfg["h"], "kind": kind, "lam": lam, "rows": []}
    for J in (200.0, 800.0):
        for R in (12.0, 18.0):
            rows = [(a, fam.parts(a, R, kind)) for a in AMPS]
            m = min_over_amp(rows, J, lam)
            a_s = m["amp_star"]
            dR = 0.5
            vals = {}
            for Rq in (R - dR, R, R + dR):
                es, kk = es_kin(fam.parts(a_s, Rq, kind), lam)
                vals[Rq] = (es, kk, es + J * J / (4.0 * kk))
            dEs = (vals[R + dR][0] - vals[R - dR][0]) / (2 * dR)
            dEk = (J * J / (4.0 * vals[R + dR][1]) - J * J / (4.0 * vals[R - dR][1])) / (2 * dR)
            dEJ = (vals[R + dR][2] - vals[R - dR][2]) / (2 * dR)
            dkin = (vals[R + dR][1] - vals[R - dR][1]) / (2 * dR)
            # total derivative with re-minimization over amp
            mp = min_over_amp([(a, fam.parts(a, R + dR, kind)) for a in AMPS], J, lam)
            mm = min_over_amp([(a, fam.parts(a, R - dR, kind)) for a in AMPS], J, lam)
            dEJ_remin = (mp["E_J"] - mm["E_J"]) / (2 * dR)
            # dilation about mu = 1 at fixed amp*
            dmu = 0.02
            ep = es_kin(fam.parts(a_s, R, kind, mu=1 + dmu), lam)
            em = es_kin(fam.parts(a_s, R, kind, mu=1 - dmu), lam)
            EJp = ep[0] + J * J / (4.0 * ep[1]); EJm = em[0] + J * J / (4.0 * em[1])
            dEJ_dmu = (EJp - EJm) / (2 * dmu)
            dEs_dmu = (ep[0] - em[0]) / (2 * dmu)
            dEk_dmu = (J * J / (4.0 * ep[1]) - J * J / (4.0 * em[1])) / (2 * dmu)
            row = {"J": J, "R": R, "amp_star": a_s, "E_J": m["E_J"], "E_stat": vals[R][0],
                   "kin": vals[R][1], "J2_over_4kin": J * J / (4.0 * vals[R][1]),
                   "omega_star": J / (2.0 * vals[R][1]),
                   "dEstat_dR": dEs, "dJ2_4kin_dR": dEk, "dEJ_dR_fixed_amp": dEJ,
                   "dEJ_dR_reminimized": dEJ_remin, "dkin_dR": dkin, "dR": dR,
                   "dominant": "J2_over_4kin" if abs(dEk) > abs(dEs) else "E_stat",
                   "dilation_mu1": {"dEJ_dmu": dEJ_dmu, "dEstat_dmu": dEs_dmu, "dJ2_4kin_dmu": dEk_dmu,
                                    "dEJ_dmu_over_EJ": dEJ_dmu / m["E_J"], "dmu": dmu}}
            out["rows"].append(row)
            log(f"decomp J={J:g} R={R:g}: amp* {a_s:.4f} E_J {m['E_J']:.3f}  dE_stat/dR {dEs:+.4f}  "
                f"d(J2/4kin)/dR {dEk:+.4f}  total {dEJ:+.4f} (remin {dEJ_remin:+.4f})  "
                f"dEJ/dmu/EJ {dEJ_dmu / m['E_J']:+.4f}")
    # the vacuum-boost cost: E_stat(amp, R) - E_stat(0) vs amp at R = 12 (exponent)
    R = 12.0
    es0 = es_kin(fam.parts(0.0, R, kind), lam)[0]
    fit = []
    for a in (0.005, 0.01, 0.02, 0.04):
        es = es_kin(fam.parts(a, R, kind), lam)[0]
        fit.append((a, es - es0))
    fit = np.array(fit)
    out["static_cost_vs_amp_R12"] = {"amp": fit[:, 0].tolist(), "dE_stat": fit[:, 1].tolist(),
                                     "note": "E_stat(amp) - E_stat(0) at lam = 1, R = 12, exp2; the hedgehog background makes it non-quartic (cross terms with the base curvature)"}
    out["runtime_s"] = round(time.time() - T0, 1)
    with open(os.path.join(DATA, f"{OUT_STEM}_{tag}.json"), "w") as f:
        json.dump(out, f, indent=1)
    return out


# ================= K5: the plane-wave probe =================
def gen(kind):
    K = np.zeros((4, 4))
    if kind == "boost_x":
        K[0, 1] = K[1, 0] = 1.0
    elif kind == "boost_y":
        K[0, 2] = K[2, 0] = 1.0
    elif kind == "rot_z":
        K[1, 2], K[2, 1] = -1.0, 1.0
    elif kind == "rot_x":
        K[2, 3], K[3, 2] = -1.0, 1.0
    return K


def density_lam(A_sp, A_t, M, lam):
    """e_lam = 4[sum_{i<j} q_lam(F_ij) + sum_i q_lam(F_0i)] + V4 (omega = 1
    absorbed into A_t) at a single point; (I1 part, V4 part) returned."""
    h = TX.h_cov_np(M[None])[0]
    tot_e = tot_h = 0.0
    A = [A_t] + A_sp
    for mu in range(4):
        for nu in range(mu + 1, 4):
            F = B3.comm_eta(A[mu], A[nu])
            tot_e += float(B3.inner_eta(F, F))
            tot_h += float(R2.q_h(F, h))
    C = [(S * G) ** p + 1.0 + DELTA ** p for p in range(1, 5)]
    P = np.eye(4); v4 = 0.0
    for p in range(1, 5):
        P = P @ (M @ ETA)
        v4 += (np.trace(P) - C[p - 1]) ** 2
    v4 *= R2.W1
    return 4.0 * R2.mix(tot_e, tot_h, lam), v4


def wave_M(vac, gens, kvs, ws, eps, x, t, lorentz=True):
    """M and exact jets (A_t, A_x, A_y, A_z) of the superposition
    W = eps sum_a G_a cos(k_a.x - w_a t); Lorentz: M = e^W vac e^W^T
    (spectrum pinned, V4 = 0); non-Lorentz: M = vac + W (D_a symmetric)."""
    W = np.zeros((4, 4)); dW = [np.zeros((4, 4)) for _ in range(4)]
    for Ga, kv, w in zip(gens, kvs, ws):
        ph = kv @ x - w * t
        W += eps * Ga * np.cos(ph)
        kmu = np.array([-w, *kv])
        for mu in range(4):
            dW[mu] += -eps * Ga * np.sin(ph) * kmu[mu]
    if not lorentz:
        return vac + W, [B3.sym4(d) for d in dW]
    A = []
    E = expm(W)
    M = E @ vac @ E.T
    for mu in range(4):
        _, dE = expm_frechet(W, dW[mu])
        A.append(B3.sym4(dE @ vac @ E.T + E @ vac @ dE.T))
    return M, A


def mode_pw(tag="pw"):
    rng = np.random.default_rng(4444)
    vac = np.diag([-S * G, 1.0, DELTA, 0.0])
    out = {"tag": tag, "vac": vac.tolist(), "single_wave": [], "crossed_waves": [],
           "identity": "single-phase profile M(k.x - w t): A_mu = k_mu M'(phase) so "
                       "F_{mu nu} = k_mu k_nu (M' eta M' - M' eta M') = 0 identically; "
                       "I1 and I1_h vanish for ANY single plane wave of ANY amplitude and profile"}
    pairs = (("boost_x", "rot_z"), ("boost_x", "boost_y"), ("rot_x", "rot_z"))
    k1, k2 = np.array((0.5, 0.0, 0.0)), np.array((0.2, 0.3, -0.1))
    w1, w2 = 0.7, 0.4
    for (ga, gb) in pairs:
        Ga, Gb = gen(ga), gen(gb)
        for lam in (0.0, 1.0):
            # single wave carrying both generators (same phase): F == 0
            x0 = rng.uniform(-1, 1, 3); t0 = 0.3
            row = {"generators": [ga, gb], "k": k1.tolist(), "w": w1, "lam": lam, "eps": {}}
            for eps in (0.05, 0.1, 0.4, 1.0):
                M, A = wave_M(vac, [Ga, Gb], [k1, k1], [w1, w1], eps, x0, t0)
                # same phase, generator pair G_a cos + G_b sin: use a shifted phase for G_b
                Mb_, Ab_ = wave_M(vac, [Ga, Gb], [k1, k1], [w1, w1], eps, x0 + np.array([np.pi / (2 * k1[0]), 0, 0]), t0)
                d_I, d_V = density_lam(A[1:], A[0], M, lam)
                scale = sum(float(B3.inner_eta(a, a)) ** 2 for a in A) * 4.0 + 1e-300
                row["eps"][f"{eps:g}"] = {"I1part": d_I, "V4": d_V, "I1part_over_A4scale": d_I / scale}
            out["single_wave"].append(row)
            log(f"pw single {ga}+{gb} lam={lam:g}: I1-part/scale " + " ".join(f"{v['I1part_over_A4scale']:.1e}" for v in row["eps"].values()))
            # crossed waves (k1 != k2): F != 0 at O(eps^2), density O(eps^4)
            dens = {}
            EPS = (0.00125, 0.0025, 0.005, 0.01, 0.02, 0.04)
            for eps in EPS:
                M, A = wave_M(vac, [Ga, Gb], [k1, k2], [w1, w2], eps, x0, t0)
                dens[eps] = density_lam(A[1:], A[0], M, lam)
            rat = [dens[EPS[i + 1]][0] / dens[EPS[i]][0] for i in range(len(EPS) - 1)]
            out["crossed_waves"].append({"generators": [ga, gb], "k1": k1.tolist(), "k2": k2.tolist(),
                                         "w1": w1, "w2": w2, "lam": lam,
                                         "density_I1part": {f"{e:g}": dens[e][0] for e in dens},
                                         "V4_max": max(v[1] for v in dens.values()),
                                         "ratio_2x_amp_I1part": rat, "eps_ladder": list(EPS)})
            log(f"pw crossed {ga}+{gb} lam={lam:g}: I1-part(eps=0.01) {dens[0.01][0]:.3e} ratios " + " ".join(f"{r:.3f}" for r in rat) + f", V4 max {max(v[1] for v in dens.values()):.1e}")
    # non-Lorentz control: M = vac + eps D cos(k.x - w t), D = diag(0, 1, 0, 0)
    D = np.diag([0.0, 1.0, 0.0, 0.0]); D2 = np.zeros((4, 4)); D2[1, 2] = D2[2, 1] = 1.0; D2[0, 0] = 0.5
    x0 = rng.uniform(-1, 1, 3); t0 = 0.3
    single, crossed = {}, {}
    for eps in (0.01, 0.02, 0.04):
        M, A = wave_M(vac, [D], [k1], [w1], eps, x0, t0, lorentz=False)
        single[eps] = density_lam(A[1:], A[0], M, 1.0)
        M, A = wave_M(vac, [D, D2], [k1, k2], [w1, w2], eps, x0, t0, lorentz=False)
        crossed[eps] = density_lam(A[1:], A[0], M, 1.0)
    out["non_lorentz_control"] = {
        "single": {"D": "diag(0,1,0,0) cos(k1.x - w1 t) added to vac (spectrum NOT pinned)",
                   "I1part": {f"{e:g}": single[e][0] for e in single}, "V4": {f"{e:g}": single[e][1] for e in single},
                   "ratio_V4": [single[0.02][1] / single[0.01][1], single[0.04][1] / single[0.02][1]]},
        "crossed": {"D": "diag(0,1,0,0) cos(k1.x - w1 t) + (sym [1,2] entries 1, [0,0] 0.5) cos(k2.x - w2 t)",
                    "I1part": {f"{e:g}": crossed[e][0] for e in crossed}, "V4": {f"{e:g}": crossed[e][1] for e in crossed},
                    "ratio_I1part": [crossed[0.02][0] / crossed[0.01][0], crossed[0.04][0] / crossed[0.02][0]],
                    "ratio_V4": [crossed[0.02][1] / crossed[0.01][1], crossed[0.04][1] / crossed[0.02][1]]},
        "note": "V4 depends on M only (no jets): a mass-like quadratic term outside the Lorentz orbit, k-independent; no quadratic gradient kernel anywhere"}
    log(f"pw non-Lorentz: single I1part {[single[e][0] for e in single]} V4 ratios {out['non_lorentz_control']['single']['ratio_V4']}; "
        f"crossed I1 ratios {out['non_lorentz_control']['crossed']['ratio_I1part']} V4 ratios {out['non_lorentz_control']['crossed']['ratio_V4']}")
    out["runtime_s"] = round(time.time() - T0, 1)
    with open(os.path.join(DATA, f"{OUT_STEM}_{tag}.json"), "w") as f:
        json.dump(out, f, indent=1)
    return out


# ================= merge =================
def mode_merge(tags):
    out = {"candidate": "L_lambda = -4[(1-lambda) I1 + lambda I1_h] - V4, h_cov = eta + 2(eta u)(eta u)^T",
           "toy": {"g": G, "delta": DELTA, "s": S},
           "stencil": "sym (1/2 fwd + bwd, density per branch)",
           "conventions": {
               "energy": "E(omega) = E_stat + omega^2 kin; E_stat = h^3 sum 4 sum_{i<j} q_lam(F_ij) + E_V; kin = h^3 sum 4 sum_i q_lam(comm_eta(a0, A_i))",
               "fixed_J": "E_J = E_stat + J^2/(4 kin), omega* = J/(2 kin*); min over amp on the grid AMPS with parabolic refinement; R* = argmin over the R grid",
               "family": "M = Qb Mb Qb^T, Mb = m5_21_8 hedgehog (m = 0), a0 = Qb a0_unit Qb^T, h = Qb^-T Qb^-1 (analytic, checked against the registry eigen-solved u); rapidity b(r) = amp tanh(r/2) w(r; R)",
               "windows": {"exp2": "exp(-(r/R)^2)", "exp1": "exp(-r/R)", "hard": "1 for r < R, cos^2 ramp to 0 at R + 2",
                           "pow1/2/3": "(1 + r/R)^-q", "rigid": "1"},
               "kin_guard": "lam = 0 rows: kin <= guard excluded (guard 0 = the kin > 0 edge; guard 0.02 = the record's numeric guard)",
               "amps": AMPS.tolist(), "lams": list(LAMS), "Js": list(JS)},
           "parts": {}}
    for t in tags:
        p = os.path.join(DATA, f"{OUT_STEM}_{t}.json")
        if os.path.exists(p):
            out["parts"][t] = json.load(open(p))
    # K3 summary: omega* drift on the localized family and the R2 rigid numbers
    r2 = json.load(open(os.path.join(DATA, "m5_32_r2_audit_lattice.json")))
    out["K3_rigid_R2_audit"] = {k: {J: {"omega_star": v["byJ"][J]["omega_star"], "L": v["L"], "h": v["h"], "n": v["n"]}
                                    for J in v["byJ"]} for k, v in r2["L2"].items()}
    drift = {}
    for key in ("lam_1_J_200", "lam_1_J_800", "lam_0.75_J_200"):
        drift[key] = {}
        for t, part in out["parts"].items():
            if "kinds" in part and "exp2" in part["kinds"] and key in part["kinds"]["exp2"]["scan"]:
                sc = part["kinds"]["exp2"]["scan"][key]
                Rw = f"R_{part['L'] / 2:g}"
                drift[key][t] = {"L": part["L"], "n": part["n"], "h": part["h"],
                                 "omega_star_wall": sc["byR"][Rw]["omega_star"],
                                 "omega_star_R12": sc["byR"]["R_12"]["omega_star"],
                                 "R_star": sc["R_test"]["R_star"]}
    out["K3_localized_drift"] = drift
    with open(os.path.join(DATA, f"{OUT_STEM}.json"), "w") as f:
        json.dump(out, f, indent=1)
    log("merged")
    return out


if __name__ == "__main__":
    mode = sys.argv[1]
    kw = dict(a.split("=", 1) for a in sys.argv[2:])
    if mode == "box":
        n, L = int(kw.get("n", 32)), float(kw.get("L", 48))
        kinds = kw.get("kinds", "exp2").split(",")
        mode_box(n, L, kinds, kw.get("tag", f"box_n{n}_L{L:g}"))
    elif mode == "alt":
        mode_alt(int(kw.get("n", 32)), float(kw.get("L", 48)),
                 kw.get("kinds", "hard,pow1,pow2,pow3").split(","),
                 kw.get("tag", f"alt_n{kw.get('n', 32)}_L{float(kw.get('L', 48)):g}"))
    elif mode == "decomp":
        mode_decomp()
    elif mode == "pw":
        mode_pw()
    elif mode == "merge":
        mode_merge(kw.get("tags", "box_n32_L48,box_n48_L72,box_n64_L96,alt_n32_L48,alt_n48_L72,decomp_n32_L48,pw").split(","))
