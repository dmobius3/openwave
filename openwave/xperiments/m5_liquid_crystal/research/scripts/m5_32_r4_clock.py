"""M5.32 R4: the CLOCK INSTRUMENT on the audited R2 candidate, with a
LOCALIZED dressing family (the compact-arena false positive of the
rigid family is the thing under test).

EQUATIONS FIRST
---------------
Candidate (R2, covariant lambda-family; every symbol as in
m5_32_r2_b_bounded.py, imported here as RB, never modified):
    L_lambda = -4 [(1 - lambda) I1 + lambda I1_h] - V4
    I1   : q_eta(F) = tr(eta F eta F^T),  I1_h : q_h(F) = tr(h F h F^T)
    h(M) = eta + 2 (eta u)(eta u)^T,  u = timelike unit eigenvector of M eta
    E(M; omega) = E_stat_lambda(M) + omega^2 kin_lambda(M)      (exact in omega)
    E_stat = 4 h^3 sum_x sum_{i<j} q_lambda(F_ij) + V4
    kin    = 4 h^3 sum_x sum_i     q_lambda(comm_eta(a0, A_i))
Localized dressing family (NEW, the R4 object): on the certified 3x3
hedgehog Mb = Qh d4 Qh^T (m5_21_8 dressed(m = 0), d4 = diag(-g, 1, delta, 0))
with its clock flow a0_unit (the t-derivative of the twisting family):
    b(r; amp, R, p) = amp tanh(r/2) exp(-(r/R)^p)
    Qb(x) = exp(b(r) n.K)  (the radial boost about the time axis)
    M = Qb Mb Qb^T,   a0 = Qb a0_unit Qb^T
(the rigid family of R2 is the R -> infinity limit; every member is a
per-cell Lorentz orbit of the vacuum so V4 = 0 on the family.)
Fixed-J electron (I1: the fixed-Noether-charge frame, omega the multiplier):
    E_J(amp, R) = E_stat(amp, R) + J^2 / (4 kin(amp, R)),  J = 2 omega kin
    omega* = J / (2 kin*),   closure: dE_J*/dJ = omega*   (envelope theorem)
The compact-arena signature: R* -> L/2 (the minimizer runs to the box) or
kin* growing with L at fixed h (IR-extensive inertia, omega* ~ 1/L).
THE GATE (task text): omega* stable across L = 48 -> 72 at fixed h = 1.5
to <= 10 % to count as an electron property.
Relaxed read (stage relax): with a0 HELD at the optimum's clock flow,
    E_J[M] = E_stat[M] + J^2 / (4 kin[M]),
    grad E_J = grad E_stat - J^2 / (4 kin^2) grad kin
grad kin (a0 held) has the same two pieces as RB.energy_grad: (a) through
F_0i = a0 eta A_i - A_i eta a0 at fixed h (dq/dF = 2W, W = (1-lambda) eta
F eta + lambda h F h, chained with d1_adj), (b) through h(M) at fixed F
with S = sum_i F_0i h F_0i^T (the RB eigenvector-perturbation formula).
FD-gated (Richardson 4-point pair, <= 1e-6 relative) before use. Descent =
RB's energy-monotone backtracking FIRE with the pinned shell (Dirichlet at
the seed values), a fixed accepted-step budget.
Morse (I21): (i) the 2x2 Hessian of E_J over (amp, R) by central FD
(a0 moving with the family); (ii) basis-free: the second variation
d^2/dt^2 E_J[M* + t D] (a0 held) along 20 random smooth localized
symmetric directions D (5-point stencil, Richardson pair), plus the two
family tangents; the index = the count of negative second variations.
Radiation window (I2): the record (m5_21_16 s 3.5) has NO quadratic
vacuum kernel (the vacuum response is quartic). Measured here: (i) on the
exact vacuum the plane-wave energy per channel scales as eps^4 (doubling
ratio 16) for E_stat and for the omega^2 weight, so no linear dispersion
and no plane-wave gap exists on the vacuum; (ii) the channel dispersion
ON THE ELECTRON BACKGROUND in a far-field annulus window w(r): for
D = w(r) cos(k.x) a0_c(x), a0_c the channel flow (clock = the family's
a0; boost = K_z M + M K_z; rotation = J_z M - M J_z),
    omega^2(k) = S2(k) / (2 K_D(k)),  S2 = d^2 E_stat[M + eps D]/d eps^2,
    K_D = kin_lambda(M; D)   (harmonic balance, L = 1/2 m q'^2 - 1/2 s q^2)
    mu_c = omega(k -> 0),  c_c = d omega / dk at the smallest k
and the verdict compares omega* with mu_c.
Dilation (I9): the family is analytic, so x -> mu x is EXACT on the
family: b_mu(r) = amp tanh(r/(2 mu)) exp(-(r/(mu R))^p) on the same Mb
(direction-only, mu-invariant in the continuum, UV-cutoff dominated on
the lattice). Predictions: E_4 ~ mu^-1, kin ~ mu^+1, V4 ~ mu^3 for a
whole-field dilation; reported for the full field and for the dressing
excess (E - E(amp = 0)); the lattice map_coordinates dilation is the
second method; dE_J/dmu at mu = 1 by central FD = the stationarity read.
Inertia decomposition (task 6): in the frame Q(x) = Qb Qh (u = Q e0 to
1e-16, measured), F_0i = Q f_0i Q^T with f antisymmetric, and
    q_lambda(F_0i) = sum_{a<b} 2 w_lambda(a, b) f_ab^2,
    w_lambda(a, b) = (1 - lambda) eta_a eta_b + lambda
so kin splits exactly over the six eigen-pairs (a, b) of d4 (labels
-g, 1, delta, 0); the delta- and g-exponents of kin at fixed texture
(d4 varied, Q held) are the quadratic analog of the sextic claim.

STAGES (python3 m5_32_r4_clock.py STAGE [--n N --L L --lam X --J J]):
    gate      family limits + the kin-gradient FD gate (exit 1 on fail)
    fixedj    the (amp, R) scan + fixed-J minimization, one box
    relax     the relaxed fixed-J read (default lam 1, J 200, n 32, L 48)
    morse     Morse certification at the localized optimum
    radiation the vacuum quartic check + the background dispersion
    dilation  the Derrick probe at the optimum
    inertia   the eigen-pair decomposition of kin
    collect   merge checkpoints -> JSON + plots
Partials: ../checkpoints/m5_32_r4/*.json (gitignored)
Out: ../data/m5_32_r4_clock.json, ../plots/m5_32_r4_clock_*.png
"""
from __future__ import annotations

import argparse
import glob
import importlib.util
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
PLOTS = os.path.join(HERE, "..", "plots")
CKPT = os.path.join(HERE, "..", "checkpoints", "m5_32_r4")


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


RB = _load("m5_32_r2_b_bounded", "m5_32_r2_b_bounded.py")
B3, B8, LAG = RB.B3, RB.B8, RB.LAG
ETA = B3.ETA
G_MAIN, S_MAIN, DELTA = RB.G_MAIN, RB.S_MAIN, RB.DELTA
LAMBDAS = (0.0, 0.75, 1.0)
R_LADDER = (6.0, 9.0, 12.0, 18.0)
P_LADDER = (1, 2)
AMPS = (0.0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05, 0.07,
        0.1, 0.15, 0.2, 0.3)
J_16 = (50.0, 200.0, 800.0)
OMEGA_T15 = (0.1, 0.2, 0.5, 1.0, 1.5)
KIN_BASE_REF = 426.5070121483972     # n 32, L 48, h 1.5 undressed (R2 record)
J_15 = tuple(2.0 * KIN_BASE_REF * w for w in OMEGA_T15)
BOXES = ((32, 48.0), (48, 72.0), (48, 48.0))
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:8.1f}s] {msg}", flush=True)


def ck(name):
    os.makedirs(CKPT, exist_ok=True)
    return os.path.join(CKPT, name)


def dump(name, obj):
    with open(ck(name), "w") as f:
        json.dump(obj, f, indent=1)


def rd_ck(name):
    p = ck(name)
    if os.path.exists(p):
        with open(p) as f:
            return json.load(f)
    return None


# ================= the localized family =================
class Family:
    """cached lattice reads of the localized dressing family on one box."""

    def __init__(self, n, L):
        self.n, self.L = n, float(L)
        self.cfg = RB.cfg_of(n, L)
        self.R, self.K, self.K2 = RB.boost_geom(self.cfg)
        self.Mb = B8.dressed(self.cfg, 0.0)
        self.a0b = B8.a0_unit(self.cfg, 0.0)
        self.cache = {}
        self.n_reads = 0

    def b_of(self, amp, R, p, mu=1.0):
        r = self.R
        return amp * np.tanh(r / (2.0 * mu)) * np.exp(-(r / (mu * R)) ** p)

    def field(self, amp, R, p, mu=1.0):
        return RB.dress(self.Mb, self.a0b, self.b_of(amp, R, p, mu),
                        self.K, self.K2)

    def read(self, amp, R, p, mu=1.0):
        key = (round(float(amp), 10), round(float(R), 8), int(p),
               round(float(mu), 8))
        if key not in self.cache:
            Md, a0d = self.field(*key)
            rd = RB.reads(Md, self.cfg, a0d)
            rd.update(amp=key[0], R=key[1], p=key[2], mu=key[3])
            self.cache[key] = rd
            self.n_reads += 1
        return self.cache[key]

    def e_stat(self, rd, lam):
        return RB.e_stat_lam(rd, lam)

    def kin(self, rd, lam):
        return RB.kin_lam(rd, lam)

    def EJ(self, rd, lam, J):
        k = self.kin(rd, lam)
        return (self.e_stat(rd, lam) + J * J / (4.0 * k)) if k > 0 else np.inf


# ================= fixed-J minimization over (amp, R) =================
def nelder_mead(f, x0, step, bounds, maxfev=40, xtol=1e-4, ftol=1e-7):
    """2D Nelder-Mead with box bounds (clipped), returns (x, f, nfev, at_bound)."""
    lo, hi = np.array(bounds)[:, 0], np.array(bounds)[:, 1]

    def clip(x):
        return np.minimum(np.maximum(x, lo), hi)
    pts = [clip(np.array(x0, float))]
    for i in range(2):
        e = np.zeros(2); e[i] = step[i]
        pts.append(clip(pts[0] + e))
    vals = [f(p) for p in pts]
    nfev = 3
    while nfev < maxfev:
        order = np.argsort(vals)
        pts = [pts[i] for i in order]; vals = [vals[i] for i in order]
        if (np.max(np.abs(np.array(pts[1:]) - pts[0])) < xtol
                and abs(vals[-1] - vals[0]) < ftol * max(1.0, abs(vals[0]))):
            break
        c = (pts[0] + pts[1]) / 2.0
        xr = clip(c + (c - pts[2])); fr = f(xr); nfev += 1
        if fr < vals[0]:
            xe = clip(c + 2.0 * (c - pts[2])); fe = f(xe); nfev += 1
            if fe < fr:
                pts[2], vals[2] = xe, fe
            else:
                pts[2], vals[2] = xr, fr
        elif fr < vals[1]:
            pts[2], vals[2] = xr, fr
        else:
            xc = clip(c + 0.5 * (pts[2] - c)); fc = f(xc); nfev += 1
            if fc < vals[2]:
                pts[2], vals[2] = xc, fc
            else:
                for i in (1, 2):
                    pts[i] = clip(pts[0] + 0.5 * (pts[i] - pts[0]))
                    vals[i] = f(pts[i]); nfev += 1
    i = int(np.argmin(vals))
    x = pts[i]
    at_bound = bool(np.any(np.abs(x - lo) < 1e-9) or np.any(np.abs(x - hi) < 1e-9))
    return x, vals[i], nfev, at_bound


def minimize_EJ(fam, lam, J, p, grid_rows, x0=None, maxfev=40):
    """discrete grid minimum, then continuous (amp, R) refinement."""
    Lh = fam.L / 2.0
    if x0 is None:
        best = min(grid_rows, key=lambda r: fam.EJ(r, lam, J))
        x0 = (best["amp"], best["R"])
        E0 = fam.EJ(best, lam, J)
    else:
        E0 = None
    if not np.isfinite(fam.EJ(min(grid_rows, key=lambda r: fam.EJ(r, lam, J)), lam, J)):
        return {"unbounded": True, "lam": lam, "J": J, "p": p,
                "kin_nonpositive": [(r["amp"], r["R"]) for r in grid_rows
                                    if fam.kin(r, lam) <= 0]}

    def f(x):
        return fam.EJ(fam.read(x[0], x[1], p), lam, J)
    x, fx, nfev, atb = nelder_mead(
        f, x0, (0.004, 1.5), [(0.0, 0.6), (3.0, Lh)], maxfev=maxfev)
    r = fam.read(x[0], x[1], p)
    k = fam.kin(r, lam)
    # R-ladder profile at the optimal amp (the interior read)
    prof = {f"R_{R:g}": fam.EJ(fam.read(x[0], R, p), lam, J) for R in R_LADDER}
    pv = list(prof.values())
    return {"unbounded": False, "lam": lam, "J": J, "p": p,
            "amp_star": float(x[0]), "R_star": float(x[1]),
            "E_total": float(fx), "E_stat": float(fam.e_stat(r, lam)),
            "E_rot": float(J * J / (4 * k)), "kin_star": float(k),
            "omega_star": float(J / (2 * k)), "E_positive": bool(fx > 0),
            "nfev": nfev, "at_bound": atb, "grid_E": E0,
            "R_ladder_EJ_at_amp_star": prof,
            "R_ladder_monotone_falling": bool(all(pv[i + 1] < pv[i] for i in range(len(pv) - 1))),
            "interior_in_R": bool(3.0 < x[1] < 0.75 * Lh and not atb),
            "runs_to_box": bool(x[1] >= 0.75 * Lh)}


def stage_fixedj(n, L, ps=P_LADDER):
    fam = Family(n, L)
    tag = f"n{n}_L{L:g}"
    out = {"n": n, "L": fam.L, "h": fam.cfg["h"], "stencil": "sym",
           "g": G_MAIN, "s": S_MAIN, "delta": DELTA,
           "family": "b = amp tanh(r/2) exp(-(r/R)^p) on m5_21_8 dressed(m=0)",
           "amps": list(AMPS), "R_ladder": list(R_LADDER), "p_ladder": list(ps),
           "J_ladder": {"m5_21_16": list(J_16), "m5_21_15": list(J_15),
                        "J15_definition": "J = 2 kin_base(n32, L48) omega_t, "
                        "kin_base fixed at the R2 record value so J is the same "
                        "charge on every box"},
           "grid": {}, "byP": {}}
    base = fam.read(0.0, R_LADDER[0], 1)
    out["base"] = {"E_stat": base["A_I1"] * 4 + base["V4"], "kin": base["kin_I1"],
                   "note": "amp = 0: the undressed hedgehog (lambda-blind)"}
    log(f"FIXJ {tag} base E {out['base']['E_stat']:.4f} kin {out['base']['kin']:.4f}")
    for p in ps:
        rows = []
        for R in R_LADDER:
            for amp in AMPS:
                r = fam.read(amp, R, p)
                rows.append(r)
            log(f"FIXJ {tag} p{p} R{R:g} grid done: kin1 at amp .02 "
                f"{fam.kin(fam.read(0.02, R, p), 1.0):.3f} E1 "
                f"{fam.e_stat(fam.read(0.02, R, p), 1.0):.3f}")
        out["grid"][f"p{p}"] = [{k: r[k] for k in ("amp", "R", "A_I1", "A_I1h",
                                                    "V4", "kin_I1", "kin_I1h",
                                                    "min_gap", "timelike_eigvec_everywhere")}
                                for r in rows]
        dump(f"fixedj_{tag}.json", out)
        res = {}
        for lam in LAMBDAS:
            tab = {}
            for J in list(J_16) + list(J_15):
                t = minimize_EJ(fam, lam, J, p, rows,
                                maxfev=(40 if lam > 0 else 0))
                if not t["unbounded"] and lam > 0 and J in J_16:
                    x0 = (t["amp_star"], t["R_star"])
                    tp = minimize_EJ(fam, lam, 1.05 * J, p, rows, x0=x0, maxfev=25)
                    tm = minimize_EJ(fam, lam, 0.95 * J, p, rows, x0=x0, maxfev=25)
                    d = (tp["E_total"] - tm["E_total"]) / (0.10 * J)
                    t["dEdJ_numeric"] = float(d)
                    t["closure_rel"] = float(abs(d / t["omega_star"] - 1.0))
                    t["closure_le_3pct"] = bool(t["closure_rel"] <= 0.03)
                tab[f"J_{J:.6g}"] = t
                log(f"FIXJ {tag} p{p} lam {lam:g} J {J:.5g}: " + (
                    "UNBOUNDED (kin <= 0 on the grid)" if t["unbounded"] else
                    f"amp* {t['amp_star']:.4f} R* {t['R_star']:.2f} E {t['E_total']:.3f} "
                    f"w* {t['omega_star']:.4f} kin* {t['kin_star']:.1f} "
                    f"interior {t['interior_in_R']} box {t['runs_to_box']}"
                    + (f" clos {t['closure_rel']:.4f}" if "closure_rel" in t else "")))
            res[f"lam_{lam:g}"] = tab
        out["byP"][f"p{p}"] = res
        out["n_reads"] = fam.n_reads
        out["wall_s"] = round(time.time() - T0, 1)
        dump(f"fixedj_{tag}.json", out)
    return out


# ================= kin gradient at lambda (a0 held) =================
def kin_grad_lam(M, a0, cfg, lam):
    """gradient of kin_lambda(M; a0) wrt symmetric M, a0 held."""
    h3, h = cfg["h"] ** 3, cfg["h"]
    u0, V, lamv, sig, k0, ok, gap = RB.tl_eig(M)
    if not np.all(ok) and lam != 0.0:
        return None
    hh = RB.h_of(u0)
    G = np.zeros_like(M)
    S = np.zeros_like(M)
    a0E_T = (a0 @ ETA).swapaxes(-1, -2)
    Ea0_T = (ETA @ a0).swapaxes(-1, -2)
    for br, wt in B3.branches(cfg["stencil"]):
        A = [B3.d1(M, ax, h, br) for ax in range(3)]
        for i in range(3):
            F = B3.comm_eta(a0, A[i])
            W = 0.0
            if lam != 1.0:
                W = W + (1 - lam) * (ETA @ F @ ETA)
            if lam != 0.0:
                hF = hh @ F
                W = W + lam * (hF @ hh)
                S += wt * ((F @ hh) @ F.swapaxes(-1, -2))
            W = 2.0 * W
            dAi = a0E_T @ W - W @ Ea0_T
            G += wt * B3.d1_adj(dAi, i, h, br)
    if lam != 0.0:
        v = 8.0 * np.einsum("ab,...bc,cd,...d->...a", ETA, S, ETA, u0)
        vu = np.einsum("...a,...ak->...k", v, V)
        l0 = np.take_along_axis(lamv, k0[..., None], axis=-1)[..., 0]
        den = l0[..., None] - lamv
        mask = (np.arange(4)[None, :] != k0.reshape(-1, 1)).reshape(den.shape)
        c = np.where(mask, sig * vu / np.where(mask, den, 1.0), 0.0)
        w = np.einsum("...ak,...k->...a", V, c)
        G += lam * ((w @ ETA)[..., :, None] * (u0 @ ETA)[..., None, :])
    return 4.0 * h3 * B3.sym4(G)


def EJ_and_grad(M, a0, cfg, lam, J):
    E, G, info = RB.energy_grad(M, cfg, lam)
    if G is None:
        return np.nan, None, info, None
    rd = RB.reads(M, cfg, a0)
    k = RB.kin_lam(rd, lam)
    if k <= 0:
        return np.nan, None, dict(info, ok=False, kin=k), None
    Gk = kin_grad_lam(M, a0, cfg, lam)
    return E + J * J / (4 * k), G - J * J / (4 * k * k) * Gk, info, (E, k)


def stage_gate():
    fam = Family(16, 24.0)
    cfg = fam.cfg
    out = {"n": 16, "L": 24.0, "checks": {}}
    chk = out["checks"]
    # (a) the R -> infinity limit reproduces the rigid family
    Md_r, a0_r = RB.dress(fam.Mb, fam.a0b, 0.02 * np.tanh(fam.R / 2.0), fam.K, fam.K2)
    rr = RB.reads(Md_r, cfg, a0_r)
    rl = fam.read(0.02, 1e9, 1)
    chk["R_inf_reproduces_rigid_Estat1"] = float(abs(RB.e_stat_lam(rr, 1) - RB.e_stat_lam(rl, 1)) / abs(RB.e_stat_lam(rr, 1)))
    chk["R_inf_reproduces_rigid_kin1"] = float(abs(RB.kin_lam(rr, 1) - RB.kin_lam(rl, 1)) / abs(RB.kin_lam(rr, 1)))
    # (b) amp = 0 is the undressed hedgehog for every R, p
    r0 = fam.read(0.0, 6.0, 2)
    eu, ev = B3.e_parts(fam.Mb, cfg)
    chk["amp0_is_undressed_Estat"] = float(abs(RB.e_stat_lam(r0, 1) - (eu + ev)) / abs(eu + ev))
    chk["amp0_is_undressed_kin"] = float(abs(RB.kin_lam(r0, 1) - B3.kin_of(fam.Mb, fam.a0b, cfg)) / B3.kin_of(fam.Mb, fam.a0b, cfg))
    # (c) V4 = 0 on the family (Lorentz orbit per cell)
    chk["V4_on_family_max"] = float(max(abs(fam.read(a, R, p)["V4"]) for a in (0.02, 0.1) for R in (6.0, 12.0) for p in (1, 2)))
    # (d) the kin gradient FD gate (a0 held), lam 0 against the certified B3.kin_grad
    Md, a0d = fam.field(0.03, 9.0, 2)
    rng = np.random.default_rng(3204)
    Gk0 = kin_grad_lam(Md, a0d, cfg, 0.0)
    Gc = B3.kin_grad(Md, a0d, cfg)
    chk["kin_grad_lam0_vs_B3"] = float(np.max(np.abs(Gk0 - Gc)) / np.max(np.abs(Gc)))
    worst = 0.0
    fd = {}
    for lam in (0.75, 1.0):
        rows = []
        Gk = kin_grad_lam(Md, a0d, cfg, lam)
        for k in range(3):
            D = B3.sym4(rng.standard_normal(Md.shape))
            D /= np.sqrt(np.sum(D * D))
            gd = float(np.sum(Gk * D))

            def kat(t):
                return RB.kin_lam(RB.reads(Md + t * D, cfg, a0d), lam)

            def fd4(eps):
                return (8 * (kat(eps) - kat(-eps)) - (kat(2 * eps) - kat(-2 * eps))) / (12 * eps)
            f1, f2 = fd4(1e-3), fd4(5e-4)
            rich = (16 * f2 - f1) / 15.0
            rel = abs(rich - gd) / max(abs(gd), 1e-300)
            worst = max(worst, rel)
            rows.append({"g_dot_D": gd, "fd_richardson": rich, "rel_err": float(rel)})
        fd[f"lam{lam:g}"] = rows
    out["kin_grad_fd_gate"] = fd
    out["kin_grad_fd_worst_rel"] = float(worst)
    # (e) the E_J gradient along one direction (composite check)
    lam, J = 1.0, 200.0
    EJ0, GJ, _, _ = EJ_and_grad(Md, a0d, cfg, lam, J)
    D = B3.sym4(rng.standard_normal(Md.shape)); D /= np.sqrt(np.sum(D * D))

    def ej(t):
        return EJ_and_grad(Md + t * D, a0d, cfg, lam, J)[0]
    f1 = (8 * (ej(1e-3) - ej(-1e-3)) - (ej(2e-3) - ej(-2e-3))) / (12e-3)
    f2 = (8 * (ej(5e-4) - ej(-5e-4)) - (ej(1e-3) - ej(-1e-3))) / (6e-3)
    rich = (16 * f2 - f1) / 15.0
    chk["EJ_grad_fd_rel"] = float(abs(rich - np.sum(GJ * D)) / abs(np.sum(GJ * D)))
    out["gate_pass"] = bool(worst <= 1e-6 and all(v <= 1e-6 for v in chk.values()))
    dump("gate.json", out)
    log("GATE " + ", ".join(f"{k} {v:.2e}" for k, v in chk.items()) + f"; kin-grad FD worst {worst:.2e}")
    log(f"GATE PASS = {out['gate_pass']}")
    return out


# ================= the relaxed fixed-J read =================
def best_point(n, L, lam, J, p=None):
    fx = rd_ck(f"fixedj_n{n}_L{L:g}.json")
    if fx is None:
        raise SystemExit("run fixedj first")
    cands = []
    for pk, res in fx["byP"].items():
        if p is not None and pk != f"p{p}":
            continue
        t = res[f"lam_{lam:g}"][f"J_{J:.6g}"]
        if not t["unbounded"]:
            cands.append(t)
    return min(cands, key=lambda t: t["E_total"])


def stage_relax(lam=1.0, J=200.0, n=32, L=48.0, steps_acc=600, it_cap=1500,
                dt0=0.02, dt_max=0.2):
    fam = Family(n, L)
    cfg = fam.cfg
    t = best_point(n, L, lam, J)
    amp, R, p = t["amp_star"], t["R_star"], t["p"]
    M, a0 = fam.field(amp, R, p)
    free = (~B3.pin_shell(n, cfg["h"]))[..., None, None].astype(float)
    EJ0, G, info, parts = EJ_and_grad(M, a0, cfg, lam, J)
    out = {"lam": lam, "J": J, "n": n, "L": L, "h": cfg["h"], "p": p,
           "start": {"amp": amp, "R": R, "E_J": EJ0, "E_stat": parts[0], "kin": parts[1],
                     "omega": J / (2 * parts[1]), "E_J_parametric": t["E_total"]},
           "method": "TRUE fixed-J descent: E_J[M] = E_stat[M] + J^2/(4 kin[M]) with "
                     "a0 HELD at the parametric optimum's clock flow (frozen-a0 "
                     "convention, m5_21_9); grad kin exact (FD-gated); RB backtracking "
                     "FIRE; pinned shell depth 1.6 (Dirichlet at the seed)",
           "steps_acc_budget": steps_acc, "it_cap": it_cap, "trace": []}
    log(f"RELAX start E_J {EJ0:.5f} E_stat {parts[0]:.4f} kin {parts[1]:.3f} omega {out['start']['omega']:.5f}")
    v = np.zeros_like(M)
    dt, alpha, n_up = dt0, 0.1, 0
    F = -G * free
    E_prev, n_acc, n_rej, it = EJ0, 0, 0, 0
    stop = "budget"
    while it < it_cap and n_acc < steps_acc:
        it += 1
        P = float(np.sum(F * v))
        if P > 0.0:
            n_up += 1
            vn, fn = np.sqrt(np.sum(v * v)), np.sqrt(np.sum(F * F))
            v = (1 - alpha) * v + alpha * (F / max(fn, 1e-300)) * vn
            if n_up > 5:
                dt = min(dt * 1.1, dt_max); alpha *= 0.99
        else:
            v[:] = 0.0; alpha = 0.1; n_up = 0
        v_try = v + dt * F
        M_try = M + dt * v_try
        E, G, info, parts = EJ_and_grad(M_try, a0, cfg, lam, J)
        reject = (G is None) or not np.isfinite(E) or E > E_prev + 1e-12 * max(abs(E_prev), 1.0)
        if reject:
            n_rej += 1; dt *= 0.5; v[:] = 0.0; alpha, n_up = 0.1, 0
            if dt < 1e-7:
                stop = "STALLED (dt collapsed)" if G is not None else "LOCUS-HIT"
                break
            continue
        n_acc += 1
        M, v, E_prev = M_try, v_try, E
        F = -G * free
        if n_acc % 50 == 0 or n_acc == steps_acc:
            row = {"it": it, "acc": n_acc, "E_J": float(E), "E_stat": parts[0], "kin": parts[1],
                   "omega": J / (2 * parts[1]), "fmax": float(np.max(np.abs(F))), "dt": dt,
                   "min_gap": info["min_gap"], "rej": n_rej}
            out["trace"].append(row)
            log(f"RELAX acc {n_acc:4d} it {it:4d} E_J {E:.5f} E_stat {parts[0]:.4f} "
                f"kin {parts[1]:.3f} omega {row['omega']:.5f} fmax {row['fmax']:.2e} dt {dt:.1e} rej {n_rej}")
    rd = RB.reads(M, cfg, a0)
    k = RB.kin_lam(rd, lam)
    end = {"E_J": float(E_prev), "E_stat": float(RB.e_stat_lam(rd, lam)), "kin": float(k),
           "omega": float(J / (2 * k)), "V4": rd["V4"], "min_gap": rd["min_gap"],
           "max_abs_M0i": rd["max_abs_M0i"]}
    st = out["start"]
    out["end"] = end
    out["drift"] = {q: float((end[q] - st[q]) / abs(st[q])) for q in ("E_J", "E_stat", "kin", "omega")}
    out["stop"] = stop
    out["accepted"], out["rejected"], out["iterations"] = n_acc, n_rej, it
    tr = out["trace"]
    if len(tr) >= 4:
        q = [r for r in tr if r["acc"] >= 0.75 * n_acc]
        dE = q[-1]["E_J"] - q[0]["E_J"]
        out["last_quarter_rel"] = float(abs(dE) / max(abs(q[-1]["E_J"]), 1.0))
        out["verdict"] = "PLATEAU" if out["last_quarter_rel"] <= 1e-3 else "FALLING (still descending at the budget)"
    # the a0-refresh sensitivity: re-read kin with the family clock at the end field? not
    # defined off the family; report instead kin with a0 rescaled to the end u-frame: skipped
    out["wall_s"] = round(time.time() - T0, 1)
    np.savez_compressed(ck(f"relax_lam{lam:g}_J{J:g}_n{n}_L{L:g}.npz"), M=M, a0=a0)
    dump(f"relax_lam{lam:g}_J{J:g}_n{n}_L{L:g}.json", out)
    log(f"RELAX end {stop}: drift " + ", ".join(f"{k} {v:+.4%}" for k, v in out["drift"].items()))
    return out


# ================= Morse =================
def second_variation(fval, eps):
    """5-point second derivative with a Richardson pair (eps, eps/2)."""
    def d2(e):
        return (-fval(2 * e) + 16 * fval(e) - 30 * fval(0.0) + 16 * fval(-e) - fval(-2 * e)) / (12 * e * e)
    a, b = d2(eps), d2(eps / 2)
    return float((16 * b - a) / 15.0), float(abs(a - b) / max(abs(b), 1e-300))


def smooth_directions(fam, Rstar, n_dir, seed=3205):
    rng = np.random.default_rng(seed)
    X, Y, Z = B3.coords(fam.n, fam.cfg["h"])
    env = np.exp(-(fam.R / Rstar) ** 2)
    free = (~B3.pin_shell(fam.n, fam.cfg["h"])).astype(float)
    kmax = 3 * 2 * np.pi / fam.L
    dirs = []
    for _ in range(n_dir):
        D = np.zeros(X.shape + (4, 4))
        for m in range(6):
            k = rng.uniform(-kmax, kmax, 3)
            ph = rng.uniform(0, 2 * np.pi)
            S = B3.sym4(rng.standard_normal((4, 4)))
            D += np.cos(k[0] * X + k[1] * Y + k[2] * Z + ph)[..., None, None] * S
        D *= (env * free)[..., None, None]
        D /= np.sqrt(np.sum(D * D))
        dirs.append(D)
    return dirs


def stage_morse(n=32, L=48.0, lams=(0.75, 1.0), Js=(200.0, 800.0), n_dir=20):
    fam = Family(n, L)
    cfg = fam.cfg
    out = {"n": n, "L": L, "h": cfg["h"], "points": {}}
    for lam in lams:
        for J in Js:
            t = best_point(n, L, lam, J)
            amp, R, p = t["amp_star"], t["R_star"], t["p"]
            tag = f"lam{lam:g}_J{J:g}"
            # (i) 2x2 Hessian over (amp, R), a0 moving with the family
            da, dR = 0.05 * max(amp, 0.01), 0.05 * R

            def EJ(a, r):
                return fam.EJ(fam.read(a, r, p), lam, J)
            E0 = EJ(amp, R)
            Haa = (EJ(amp + da, R) - 2 * E0 + EJ(amp - da, R)) / da ** 2
            HRR = (EJ(amp, R + dR) - 2 * E0 + EJ(amp, R - dR)) / dR ** 2
            HaR = (EJ(amp + da, R + dR) - EJ(amp + da, R - dR) - EJ(amp - da, R + dR) + EJ(amp - da, R - dR)) / (4 * da * dR)
            H = np.array([[Haa, HaR], [HaR, HRR]])
            ev = np.linalg.eigvalsh(H)
            ga = (EJ(amp + da, R) - EJ(amp - da, R)) / (2 * da)
            gR = (EJ(amp, R + dR) - EJ(amp, R - dR)) / (2 * dR)
            # (ii) basis-free second variations, a0 held
            M0, a0 = fam.field(amp, R, p)
            dirs = smooth_directions(fam, R, n_dir)
            Mp, _ = fam.field(amp + da, R, p); Mm, _ = fam.field(amp - da, R, p)
            Da = (Mp - Mm) / (2 * da); Da /= np.sqrt(np.sum(Da * Da))
            Mp, _ = fam.field(amp, R + dR, p); Mm, _ = fam.field(amp, R - dR, p)
            DR = (Mp - Mm) / (2 * dR); DR /= np.sqrt(np.sum(DR * DR))
            rows = []
            for nm, D in [("family_amp", Da), ("family_R", DR)] + [(f"rand_{i}", d) for i, d in enumerate(dirs)]:
                def fval(tt):
                    rd = RB.reads(M0 + tt * D, cfg, a0)
                    k = RB.kin_lam(rd, lam)
                    return RB.e_stat_lam(rd, lam) + J * J / (4 * k)
                sv, rich_err = second_variation(fval, 0.5)
                # split: static vs rotational second variation
                svs, _ = second_variation(lambda tt: RB.e_stat_lam(RB.reads(M0 + tt * D, cfg, a0), lam), 0.5)
                rows.append({"dir": nm, "d2EJ": sv, "richardson_rel": rich_err, "d2Estat": svs,
                             "d2Erot": sv - svs})
                log(f"MORSE {tag} {nm}: d2E_J {sv:+.5g} (stat {svs:+.5g}) rich {rich_err:.1e}")
            neg = [r["dir"] for r in rows if r["d2EJ"] < 0]
            out["points"][tag] = {
                "amp_star": amp, "R_star": R, "p": p, "E_J": E0,
                "hessian_2x2": H.tolist(), "hessian_eigs": ev.tolist(),
                "grad_2x2": [ga, gR], "grad_rel": [ga * da / max(abs(E0), 1), gR * dR / max(abs(E0), 1)],
                "family_index": int(np.sum(ev < 0)),
                "second_variations": rows, "n_dir": len(rows),
                "n_negative": len(neg), "negative_dirs": neg,
                "basis_free_index": len(neg),
                "eps": 0.5, "eps_note": "D unit Frobenius over the grid; eps 0.5 = per-cell "
                                        "amplitude ~ 1e-3 (Richardson pair eps, eps/2)"}
            dump(f"morse_n{n}_L{L:g}.json", out)
    out["wall_s"] = round(time.time() - T0, 1)
    dump(f"morse_n{n}_L{L:g}.json", out)
    return out


# ================= radiation window =================
def channel_flows(M, a0_fam):
    Kz = np.zeros((4, 4)); Kz[0, 3] = Kz[3, 0] = 1.0
    Jz = np.zeros((4, 4)); Jz[1, 2], Jz[2, 1] = -1.0, 1.0
    return {"clock": a0_fam, "boost": Kz @ M + M @ Kz, "rotation": Jz @ M - M @ Jz}


def stage_radiation(n=32, L=48.0, lam_list=(0.75, 1.0), J=200.0):
    fam = Family(n, L)
    cfg = fam.cfg
    X, Y, Z = B3.coords(n, cfg["h"])
    out = {"n": n, "L": L, "h": cfg["h"], "vacuum": {}, "background": {}}
    # (i) the exact vacuum: plane-wave energy scaling per channel
    Mv = np.broadcast_to(B3.vac4(cfg), X.shape + (4, 4)).copy()
    k1 = 2 * np.pi / L * 2
    ch = channel_flows(Mv, None)
    ch["clock"] = np.broadcast_to(B8.G1 @ B3.vac4(cfg) - B3.vac4(cfg) @ B8.G1, Mv.shape).copy()
    # a single-generator plane wave has dD ~ D so F_ij = 0 and comm(Dt, dD) = 0
    # identically; the probe is a TWO-generator wave (channel c along z plus a
    # second channel c' along x) with the velocity pattern on the partner
    # generator, so both the curvature and the omega^2 weight are exercised
    partner = {"clock": "rotation", "boost": "clock", "rotation": "boost"}
    for lam in lam_list:
        vac = {}
        for nm, a0c in ch.items():
            rows = {}
            a0p = ch[partner[nm]]
            for eps in (0.01, 0.02):
                D = np.cos(k1 * Z)[..., None, None] * a0c + np.cos(k1 * X)[..., None, None] * a0p
                Dt = np.sin(k1 * Z)[..., None, None] * a0p
                rd = RB.reads(Mv + eps * D, cfg, eps * Dt)
                E0 = RB.e_stat_lam(RB.reads(Mv, cfg, Dt), lam)
                rows[f"eps_{eps:g}"] = {"E_stat_excess": RB.e_stat_lam(rd, lam) - E0,
                                        "curvature_excess": RB.e_stat_lam(rd, lam) - E0 - rd["V4"],
                                        "omega2_weight": RB.kin_lam(rd, lam), "V4": rd["V4"]}
            r1, r2 = rows["eps_0.01"], rows["eps_0.02"]
            vac[nm] = {"rows": rows, "partner_generator": partner[nm],
                       "doubling_ratio_curvature": float(r2["curvature_excess"] / r1["curvature_excess"]) if abs(r1["curvature_excess"]) > 1e-300 else None,
                       "doubling_ratio_Estat": float(r2["E_stat_excess"] / r1["E_stat_excess"]) if abs(r1["E_stat_excess"]) > 1e-300 else None,
                       "doubling_ratio_kin": float(r2["omega2_weight"] / r1["omega2_weight"]) if abs(r1["omega2_weight"]) > 1e-300 else None,
                       "doubling_ratio_V4": float(r2["V4"] / r1["V4"]) if abs(r1["V4"]) > 1e-300 else None}
            log(f"RAD vacuum lam {lam:g} {nm}: ratios curv {vac[nm]['doubling_ratio_curvature']} Estat {vac[nm]['doubling_ratio_Estat']} kin {vac[nm]['doubling_ratio_kin']} V4 {vac[nm]['doubling_ratio_V4']}")
        out["vacuum"][f"lam_{lam:g}"] = vac
    out["vacuum"]["reading"] = ("doubling ratio 16 = quartic: no quadratic kinetic or static kernel "
                                "on the exact vacuum, no linear plane wave, no plane-wave gap (the "
                                "m5_21_16 s 3.5 record reproduced for the lambda-family)")
    # (ii) the electron background: channel dispersion in the far-field annulus
    r = fam.R
    r_in, r_out = 8.0, min(0.42 * L, 20.0)
    w = 0.5 * (1 + np.tanh((r - r_in) / 1.5)) * 0.5 * (1 - np.tanh((r - r_out) / 1.5))
    ks = [0.0] + [2 * np.pi / L * m for m in (1, 2, 4)]
    for lam in lam_list:
        t = best_point(n, L, lam, J)
        M0, a0 = fam.field(t["amp_star"], t["R_star"], t["p"])
        ch = channel_flows(M0, a0)
        bg = {"point": {k: t[k] for k in ("amp_star", "R_star", "p", "omega_star")}, "window": [r_in, r_out], "channels": {}}
        for nm, a0c in ch.items():
            rows = []
            for kk in ks:
                for axis, C in (("z", Z), ("x", X)):
                    if kk == 0.0 and axis == "x":
                        continue
                    D = (w * np.cos(kk * C))[..., None, None] * a0c
                    D /= np.sqrt(np.sum(D * D))
                    KD = RB.kin_lam(RB.reads(M0, cfg, D), lam)
                    S2, rich = second_variation(lambda tt: RB.e_stat_lam(RB.reads(M0 + tt * D, cfg, a0), lam), 0.5)
                    om2 = S2 / (2 * KD) if KD > 0 else None
                    rows.append({"k": kk, "axis": axis, "S2": S2, "richardson_rel": rich, "K_D": KD,
                                 "omega2": om2, "omega": (np.sqrt(om2) if om2 is not None and om2 >= 0 else None)})
                    log(f"RAD bg lam {lam:g} {nm} k {kk:.3f}{axis}: S2 {S2:+.4g} K_D {KD:.4g} omega2 {om2}")
            mu2 = rows[0]["omega2"]
            z = [rw for rw in rows if rw["axis"] == "z" and rw["k"] > 0]
            slope = None
            if mu2 is not None and z and z[0]["omega2"] is not None:
                slope = (z[0]["omega2"] - mu2) / z[0]["k"] ** 2
            bg["channels"][nm] = {"rows": rows, "mu2": mu2, "mu": (np.sqrt(mu2) if mu2 is not None and mu2 >= 0 else None),
                                  "mu2_negative": bool(mu2 is not None and mu2 < 0),
                                  "c2_smallest_k": slope,
                                  "omega_star_below_mu": (bool(t["omega_star"] ** 2 < mu2) if mu2 is not None else None)}
        out["background"][f"lam_{lam:g}"] = bg
        dump(f"radiation_n{n}_L{L:g}.json", out)
    out["wall_s"] = round(time.time() - T0, 1)
    dump(f"radiation_n{n}_L{L:g}.json", out)
    return out


# ================= dilation =================
def stage_dilation(n=32, L=48.0, lam_list=(0.75, 1.0), J=200.0, mus=(0.8, 1.25)):
    from scipy.ndimage import map_coordinates
    fam = Family(n, L)
    cfg = fam.cfg
    out = {"n": n, "L": L, "h": cfg["h"], "points": {}}
    for lam in lam_list:
        t = best_point(n, L, lam, J)
        amp, R, p = t["amp_star"], t["R_star"], t["p"]
        rows = {}
        for mu in (0.95, 1.0, 1.05) + tuple(mus):
            rd = fam.read(amp, R, p, mu)
            r0 = fam.read(0.0, R, p, mu)
            k = RB.kin_lam(rd, lam)
            rows[f"mu_{mu:g}"] = {"E4": RB.e_stat_lam(rd, lam) - rd["V4"], "V4": rd["V4"], "kin": k,
                                  "E_J": RB.e_stat_lam(rd, lam) + J * J / (4 * k),
                                  "E4_dress_excess": RB.e_stat_lam(rd, lam) - RB.e_stat_lam(r0, lam),
                                  "kin_dress_excess": k - RB.kin_lam(r0, lam)}
        one = rows["mu_1"]
        scal = {}
        for mu in mus:
            q = rows[f"mu_{mu:g}"]
            scal[f"mu_{mu:g}"] = {
                "E4_ratio": q["E4"] / one["E4"], "E4_pred_mu^-1": 1 / mu,
                "E4_excess_ratio": q["E4_dress_excess"] / one["E4_dress_excess"],
                "kin_ratio": q["kin"] / one["kin"], "kin_pred_mu^+1": mu,
                "kin_excess_ratio": q["kin_dress_excess"] / one["kin_dress_excess"],
                "V4_ratio": (q["V4"] / one["V4"]) if abs(one["V4"]) > 1e-300 else None, "V4_pred_mu^3": mu ** 3}
        dEdmu = (rows["mu_1.05"]["E_J"] - rows["mu_0.95"]["E_J"]) / 0.10
        # lattice map_coordinates dilation (second method) at mu = 1.25
        M0, a0 = fam.field(amp, R, p)
        c = (n - 1) / 2.0
        idx = np.arange(n, dtype=float)
        lat = {}
        for mu in mus:
            co = (idx - c) / mu + c
            CI, CJ, CK = np.meshgrid(co, co, co, indexing="ij")
            Mm = np.zeros_like(M0); am = np.zeros_like(a0)
            for a in range(4):
                for b in range(4):
                    Mm[..., a, b] = map_coordinates(M0[..., a, b], [CI, CJ, CK], order=3, mode="nearest")
                    am[..., a, b] = map_coordinates(a0[..., a, b], [CI, CJ, CK], order=3, mode="nearest")
            rd = RB.reads(B3.sym4(Mm), cfg, B3.sym4(am))
            lat[f"mu_{mu:g}"] = {"E4_ratio": (RB.e_stat_lam(rd, lam) - rd["V4"]) / one["E4"],
                                 "kin_ratio": RB.kin_lam(rd, lam) / one["kin"], "V4": rd["V4"]}
        out["points"][f"lam_{lam:g}"] = {
            "amp_star": amp, "R_star": R, "p": p, "rows": rows, "scaling": scal,
            "lattice_interp_method": lat,
            "dEJ_dmu_at_1": float(dEdmu), "dEJ_dmu_rel": float(dEdmu / one["E_J"]),
            "stationary_under_dilation_1pct": bool(abs(dEdmu / one["E_J"]) < 0.01),
            "note": "the hedgehog core Mb is direction-only (dilation-invariant in the continuum, "
                    "UV-cutoff dominated on the lattice); Derrick's mu^-1 applies to the dressing "
                    "excess and the lattice-cutoff core does not scale"}
        log(f"DIL lam {lam:g}: E4 ratios {[(m, round(scal[f'mu_{m:g}']['E4_ratio'], 4)) for m in mus]} "
            f"excess {[(m, round(scal[f'mu_{m:g}']['E4_excess_ratio'], 4)) for m in mus]} kin {[(m, round(scal[f'mu_{m:g}']['kin_ratio'], 4)) for m in mus]} "
            f"dEJ/dmu {dEdmu:+.4g} ({dEdmu / one['E_J']:+.3%})")
    dump(f"dilation_n{n}_L{L:g}.json", out)
    return out


# ================= inertia decomposition =================
def stage_inertia(n=32, L=48.0, lam_list=(0.0, 0.75, 1.0), J=200.0):
    fam = Family(n, L)
    cfg = fam.cfg
    h3 = cfg["h"] ** 3
    X, Y, Z = B3.coords(n, cfg["h"])
    rho = np.sqrt(X * X + Y * Y); phi = np.arctan2(Y, X); th = -np.arctan2(Z, rho)
    Qh = np.einsum("...ab,...bc->...ac", B8.rot_field(B8.G3, phi), B8.rot_field(B8.G2, th))
    labels = ["-g", "1", "delta", "0"]
    out = {"n": n, "L": L, "h": cfg["h"], "points": {}}
    for lam in lam_list:
        t = best_point(n, L, 1.0 if lam == 0.0 else lam, J)
        amp, R, p = t["amp_star"], t["R_star"], t["p"]
        bl = fam.b_of(amp, R, p)
        Qb = (np.eye(4)[None, None, None] + np.sinh(bl)[..., None, None] * fam.K
              + (np.cosh(bl) - 1.0)[..., None, None] * fam.K2)
        Q = Qb @ Qh
        Qi = np.linalg.inv(Q)
        M0, a0 = fam.field(amp, R, p)
        u0 = RB.tl_eig(M0)[0]
        uchk = float(np.max(np.minimum(np.max(np.abs(u0 - Q[..., :, 0]), -1), np.max(np.abs(u0 + Q[..., :, 0]), -1))))
        pair = np.zeros((4, 4))
        tot = 0.0
        for br, wt in B3.branches(cfg["stencil"]):
            A = [B3.d1(M0, ax, cfg["h"], br) for ax in range(3)]
            for i in range(3):
                F = B3.comm_eta(a0, A[i])
                f = Qi @ F @ Qi.swapaxes(-1, -2)
                for a in range(4):
                    for b in range(a + 1, 4):
                        wl = (1 - lam) * ETA[a, a] * ETA[b, b] + lam
                        pair[a, b] += wt * 4 * h3 * 2 * wl * np.sum(f[..., a, b] ** 2)
                tot += wt * 4 * h3 * np.sum(RB.q_eta(F) * (1 - lam) + lam * RB.q_h(F, RB.h_of(u0)))
        kin_direct = RB.kin_lam(fam.read(amp, R, p), lam)
        table = {f"({labels[a]},{labels[b]})": float(pair[a, b]) for a in range(4) for b in range(a + 1, 4)}
        # eigenvalue exponents at fixed texture: vary d4 with Q held
        def kin_at(g, d):
            d4 = np.diag([-g, 1.0, d, 0.0])
            Mx = B3.sym4(np.einsum("...ab,bc,...dc->...ad", Q, d4, Q))
            # clock flow of the same texture: rotation G1 in the (delta, 0) plane
            Gm = B8.G1
            a0x = np.einsum("...ab,bc,...dc->...ad", Q, Gm @ d4 - d4 @ Gm, Q)
            # normalize to the family a0 convention (a0_unit = d/dt of the twisting family
            # = Q (G1 d4 - d4 G1) Q^T exactly at t = 0)
            return RB.kin_lam(RB.reads(Mx, cfg, B3.sym4(a0x)), lam)
        k0 = kin_at(G_MAIN, DELTA)
        a0chk = float(abs(k0 - kin_direct) / kin_direct)
        dd = [0.27, 0.33]; gg = [28.0, 36.0]
        kd = [kin_at(G_MAIN, d) for d in dd]
        kg = [kin_at(g, DELTA) for g in gg]
        exp_delta = float(np.log(kd[1] / kd[0]) / np.log(dd[1] / dd[0]))
        exp_g = float(np.log(kg[1] / kg[0]) / np.log(gg[1] / gg[0]))
        out["points"][f"lam_{lam:g}"] = {
            "amp_star": amp, "R_star": R, "p": p, "u_is_Qe0_maxdev": uchk,
            "kin_direct": kin_direct, "kin_pair_sum": float(pair.sum()), "kin_frame_total": float(tot),
            "pair_split": table,
            "pair_share": {k: v / pair.sum() for k, v in table.items()},
            "a0_analytic_vs_family_rel": a0chk,
            "delta_exponent": exp_delta, "g_exponent": exp_g,
            "kin_vs_delta": dict(zip(map(str, dd), kd)), "kin_vs_g": dict(zip(map(str, gg), kg)),
            "note": "quadratic-action analog: each pair (a,b) carries weight w_lambda(a,b) f_ab^2 "
                    "with f_ab built from eigenvalue differences; the exponents are the measured "
                    "kin ~ delta^x g^y at fixed texture (the clock rotates the (delta, 0) plane)"}
        log(f"INERTIA lam {lam:g}: kin {kin_direct:.3f} pair-sum {pair.sum():.3f}; shares "
            + ", ".join(f"{k} {v:.3f}" for k, v in out['points'][f'lam_{lam:g}']['pair_share'].items())
            + f"; delta-exp {exp_delta:.3f} g-exp {exp_g:.3f}; a0 check {a0chk:.1e}")
    dump(f"inertia_n{n}_L{L:g}.json", out)
    return out


# ================= collect =================
def stage_collect():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out = {"task": "M5.32 R4: the clock instrument on the localized dressing family",
           "candidate": "L_lambda = -4 [(1 - lambda) I1 + lambda I1_h] - V4, lambda in {0.75, 1} (+ lambda 0 control)",
           "point": {"g": G_MAIN, "s": S_MAIN, "delta": DELTA, "stencil": "sym"},
           "gate": rd_ck("gate.json"),
           "fixedj": {os.path.basename(p)[7:-5]: json.load(open(p)) for p in sorted(glob.glob(ck("fixedj_n*.json")))},
           "relax": {os.path.basename(p)[6:-5]: json.load(open(p)) for p in sorted(glob.glob(ck("relax_*.json")))},
           "morse": rd_ck("morse_n32_L48.json"), "radiation": rd_ck("radiation_n32_L48.json"),
           "dilation": rd_ck("dilation_n32_L48.json"), "inertia": rd_ck("inertia_n32_L48.json")}
    # the summary table + the L-drift gate
    summ = []
    for tag, fx in out["fixedj"].items():
        for pk, res in fx["byP"].items():
            for lk, tab in res.items():
                for Jk, t in tab.items():
                    row = {"box": tag, "n": fx["n"], "L": fx["L"], "h": fx["h"], "p": int(pk[1:]),
                           "lam": float(lk[4:]), "J": t["J"]}
                    if t["unbounded"]:
                        row.update(unbounded=True)
                    else:
                        row.update({k: t.get(k) for k in ("amp_star", "R_star", "omega_star", "E_total", "kin_star",
                                                          "E_stat", "closure_rel", "interior_in_R", "runs_to_box",
                                                          "at_bound", "R_ladder_monotone_falling")})
                    summ.append(row)
    out["summary_rows"] = summ
    drift = {}
    for p in P_LADDER:
        for lam in (0.75, 1.0):
            for J in J_16 + J_15:
                a = [r for r in summ if r["box"] == "n32_L48" and r["p"] == p and r["lam"] == lam and abs(r["J"] - J) < 1e-6 and not r.get("unbounded")]
                b = [r for r in summ if r["box"] == "n48_L72" and r["p"] == p and r["lam"] == lam and abs(r["J"] - J) < 1e-6 and not r.get("unbounded")]
                c = [r for r in summ if r["box"] == "n48_L48" and r["p"] == p and r["lam"] == lam and abs(r["J"] - J) < 1e-6 and not r.get("unbounded")]
                if a and b:
                    d = (b[0]["omega_star"] - a[0]["omega_star"]) / a[0]["omega_star"]
                    drift[f"p{p}_lam{lam:g}_J{J:.6g}"] = {
                        "omega_L48": a[0]["omega_star"], "omega_L72": b[0]["omega_star"], "drift_L": float(d),
                        "gate_le_10pct": bool(abs(d) <= 0.10),
                        "kin_L48": a[0]["kin_star"], "kin_L72": b[0]["kin_star"],
                        "R_L48": a[0]["R_star"], "R_L72": b[0]["R_star"],
                        "E_L48": a[0]["E_total"], "E_L72": b[0]["E_total"],
                        "omega_h1": c[0]["omega_star"] if c else None,
                        "drift_h": float((c[0]["omega_star"] - a[0]["omega_star"]) / a[0]["omega_star"]) if c else None}
    out["omega_drift_gate"] = drift
    passes = [v["gate_le_10pct"] for v in drift.values()]
    out["omega_drift_verdict"] = {"n_cases": len(passes), "n_pass": int(sum(passes)),
                                  "all_pass": bool(passes and all(passes))}
    os.makedirs(PLOTS, exist_ok=True)
    # plot 1: E_J(R) ladder profiles + omega* vs L
    if out["fixedj"]:
        fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
        for tag, fx in out["fixedj"].items():
            for pk, res in fx["byP"].items():
                for lk in ("lam_0.75", "lam_1"):
                    t = res[lk]["J_200"]
                    if t["unbounded"]:
                        continue
                    prof = t["R_ladder_EJ_at_amp_star"]
                    ax[0].plot(R_LADDER, list(prof.values()), "o-", label=f"{tag} {pk} {lk}")
        ax[0].set_xlabel("R (envelope radius)"); ax[0].set_ylabel("E_J at amp*, J = 200")
        ax[0].set_title("fixed-J energy along the R ladder"); ax[0].legend(fontsize=6)
        for key, d in drift.items():
            if "J200" in key or "J50" in key or "J800" in key:
                ax[1].plot([48, 72], [d["omega_L48"], d["omega_L72"]], "o-", label=key)
        ax[1].set_xlabel("L (h = 1.5)"); ax[1].set_ylabel("omega*"); ax[1].set_title("omega* across the L ladder (gate <= 10 %)")
        ax[1].legend(fontsize=6)
        for tag, fx in out["fixedj"].items():
            for pk, res in fx["byP"].items():
                tab = res["lam_1"]
                Js = [t["J"] for t in tab.values() if not t["unbounded"]]
                ws = [t["omega_star"] for t in tab.values() if not t["unbounded"]]
                o = np.argsort(Js)
                ax[2].plot(np.array(Js)[o], np.array(ws)[o], "s-", label=f"{tag} {pk} lam 1")
        ax[2].set_xscale("log"); ax[2].set_xlabel("J"); ax[2].set_ylabel("omega*"); ax[2].set_title("omega*(J), lambda = 1")
        ax[2].legend(fontsize=6)
        fig.suptitle("M5.32 R4: the localized fixed-J electron")
        fig.tight_layout(); fig.savefig(os.path.join(PLOTS, "m5_32_r4_clock_fixedj.png"), dpi=110)
    # plot 2: relax + morse + dispersion
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
    for tag, rl in out["relax"].items():
        tr = rl["trace"]
        if tr:
            ax[0].plot([r["acc"] for r in tr], [r["E_J"] for r in tr], "o-", label=f"E_J {tag}")
            ax[0].plot([r["acc"] for r in tr], [r["E_stat"] for r in tr], "s--", label=f"E_stat {tag}")
    ax[0].set_xlabel("accepted steps"); ax[0].set_title("relaxed fixed-J read (a0 held)"); ax[0].legend(fontsize=7)
    if out["morse"]:
        for tag, pt in out["morse"]["points"].items():
            sv = [r["d2EJ"] for r in pt["second_variations"]]
            ax[1].plot(range(len(sv)), sv, "o", label=f"{tag} index {pt['basis_free_index']}")
        ax[1].axhline(0, color="k", lw=0.5); ax[1].set_yscale("symlog", linthresh=1e-3)
        ax[1].set_xlabel("direction (0, 1 = family tangents)"); ax[1].set_title("second variation of E_J"); ax[1].legend(fontsize=7)
    if out["radiation"]:
        for lk, bg in out["radiation"]["background"].items():
            for nm, chn in bg["channels"].items():
                z = [r for r in chn["rows"] if r["axis"] == "z" and r["omega2"] is not None]
                ax[2].plot([r["k"] for r in z], [r["omega2"] for r in z], "o-", label=f"{lk} {nm}")
            ax[2].axhline(bg["point"]["omega_star"] ** 2, ls=":", label=f"{lk} omega*^2")
        ax[2].set_xlabel("k (along z)"); ax[2].set_ylabel("omega^2(k)"); ax[2].set_yscale("symlog", linthresh=1e-3)
        ax[2].set_title("channel dispersion on the electron background"); ax[2].legend(fontsize=6)
    fig.tight_layout(); fig.savefig(os.path.join(PLOTS, "m5_32_r4_clock_certify.png"), dpi=110)
    out["collected_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    with open(os.path.join(DATA, "m5_32_r4_clock.json"), "w") as f:
        json.dump(out, f, indent=1)
    log("COLLECT written")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["gate", "fixedj", "relax", "morse", "radiation", "dilation", "inertia", "collect"])
    ap.add_argument("--n", type=int, default=32)
    ap.add_argument("--L", type=float, default=48.0)
    ap.add_argument("--lam", type=float, default=1.0)
    ap.add_argument("--J", type=float, default=200.0)
    ap.add_argument("--steps", type=int, default=600)
    a = ap.parse_args()
    if a.stage == "gate":
        sys.exit(0 if stage_gate()["gate_pass"] else 1)
    elif a.stage == "fixedj":
        stage_fixedj(a.n, a.L)
    elif a.stage == "relax":
        stage_relax(a.lam, a.J, a.n, a.L, steps_acc=a.steps)
    elif a.stage == "morse":
        stage_morse(a.n, a.L)
    elif a.stage == "radiation":
        stage_radiation(a.n, a.L)
    elif a.stage == "dilation":
        stage_dilation(a.n, a.L)
    elif a.stage == "inertia":
        stage_inertia(a.n, a.L)
    else:
        stage_collect()
    log(f"done {a.stage}")


if __name__ == "__main__":
    main()
