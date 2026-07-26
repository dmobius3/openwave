"""M5.21.9a INDEPENDENT ADVERSARIAL AUDIT (second agent, own code).

Attempts to REFUTE the negative-delta claims CL1..CL8 of task
M5.21.9a with an independent implementation. Nothing is imported
from m5_21_9_a_negdelta.py (the script under audit). The pipeline
machinery reused here is the auditor's own prior work in
m5_21_8_audit_check.py (explicit-matrix ansatz build, plain central
FD, midpoint time average, composite Simpson), cross-gated against
a Rodrigues build in the M5.21.8 audit and re-gated here in the
negative-delta, negative-m regime.

Own additions for this audit:
  - exact polynomial algebra over integers (dict coefficients,
    factors multiplied out IN CODE, not hand-expanded) for the CL1
    factorization and the CL3 static identity, plus a tampered
    comparator negative control
  - exact rational (fractions.Fraction) evaluation of the omega^2
    coefficient across the delta < 0 half-plane and at the corner
    (g, delta) = (1e10, -+1e-10), killing float doubt, plus a
    naive-float64 fragility probe at the same corner
  - a t = 0 snapshot evaluator and a static clock-phase-average
    evaluator (alpha average, zero time derivatives), the
    omega -> 0 limit route for CL5
  - a generalized builder with arbitrary diagonal d for the CL7
    spectrum-shift test, plus the exact algebra
    M(d + c I) - M(d) = c Q Q^T checked directly
  - r-structure probe: E(r, om) = S / r^2 + K om^2 with K measured
    at r = 0.8, 1.0, 1.3, 2.0 (is the CL4 magnitude r-fortuitous?)
  - negative controls: omega probe-pair swaps, om choice in the
    m-scan, cone ladder extended to Delta = 0.01

Context: the author (reply 2026-07-19) suggested the region that
Minimize routes to -infinity for 0 < delta <= 2, g > 1 might be
bounded for delta <= 0, "not an issue, especially when shifting
spectrum". CL1..CL8 are the task's formalization of that check.
Verdicts go to data/m5_21_9_negdelta_audit.json.

Headless; prints a summary table.
"""
from __future__ import annotations

import json
import os
import sys
from fractions import Fraction as Fr

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
sys.path.insert(0, HERE)

from m5_21_8_audit_check import (  # noqa: E402
    XI, build_M, build_M_rod, fd_dM, plane_E, _simpson_w)


def mstar(g):
    return 0.5 * np.log((g + 1.0) / (g - 1.0))


# ================ exact polynomial algebra (own) ====================
def padd(a, b):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0) + v
    return {k: v for k, v in out.items() if v != 0}


def pmul(a, b):
    out = {}
    for (i, j), u in a.items():
        for (k, ll), v in b.items():
            key = (i + k, j + ll)
            out[key] = out.get(key, 0) + u * v
    return {k: v for k, v in out.items() if v != 0}


def pscale(a, s):
    return {k: s * v for k, v in a.items() if s * v != 0}


# keys are (power of g, power of d)
G = {(1, 0): 1}
D = {(0, 1): 1}
GM1 = {(1, 0): 1, (0, 0): -1}

# B = w-coefficient bracket of the transcribed Hm numerator:
# -2 g (d - 1) d + d^2 + g^2 (d - 2) d, multiplied out in code
B_POLY = padd(padd(
    pmul(pscale(G, -2), pmul(padd(D, {(0, 0): -1}), D)),
    pmul(D, D)),
    pmul(pmul(G, G), pmul(padd(D, {(0, 0): -2}), D)))


def coeff_transcribed(g, d):
    """omega^2 coefficient of the transcribed Hm; the r^2 in w
    cancels the 1/r^2 prefactor, so it is r-free. Works on floats
    and on Fractions."""
    return (8 * d * d * (-2 * g * (d - 1) * d + d * d
                         + g * g * (d - 2) * d) / ((g - 1) ** 2))


def coeff_factored(g, d):
    """the claimed factorization 8 d^3 ((g-1) d - 2g) / (g-1)."""
    return 8 * d ** 3 * ((g - 1) * d - 2 * g) / (g - 1)


def Hm_frac(g, d, r, om2):
    """the transcribed density, exact rational arithmetic."""
    w = -1 + 2 * r * r * om2
    return (4 * d * d * (-2 * g * (d - 1) * d * w + d * d * w
                         + g * g * (2 + (d - 2) * d * w))
            / ((g - 1) ** 2 * r * r))


def fin_frac(g, d, r0):
    """the transcribed finite-branch value, exact rational."""
    return ((8 * g * g * d * d + 8 * (g - 1) * g * d ** 3
             - 4 * (g - 1) ** 2 * d ** 4) / ((g - 1) ** 2 * r0))


# ================ claim CL1: exact factorization ====================
def cl1():
    # coefficient-level identity, denominators cleared:
    # [8 d^2 B] (g-1) == [8 d^3 ((g-1) d - 2g)] (g-1)^2 / (g-1)
    lhs_num = pscale(pmul(D, pmul(D, B_POLY)), 8)      # 8 d^2 B
    rhs_num = pscale(pmul(pmul(D, pmul(D, D)),
                          padd(pmul(GM1, D), pscale(G, -2))), 8)
    ident = pmul(lhs_num, GM1) == pmul(rhs_num, pmul(GM1, GM1))
    # comparator negative control: tamper 2g -> 2 in the factor
    bad_num = pscale(pmul(pmul(D, pmul(D, D)),
                          padd(pmul(GM1, D), {(0, 0): -2})), 8)
    control_detects = (pmul(lhs_num, GM1)
                       != pmul(bad_num, pmul(GM1, GM1)))
    # exact random cross-check of the two full rational forms
    rng = np.random.default_rng(21)
    n_exact, n_mismatch = 0, 0
    while n_exact < 400:
        g = Fr(int(rng.integers(-900, 900)), int(rng.integers(1, 97)))
        d = Fr(int(rng.integers(-900, 900)), int(rng.integers(1, 97)))
        if g == 1:
            continue
        if coeff_transcribed(g, d) != coeff_factored(g, d):
            n_mismatch += 1
        n_exact += 1
    # exact sign law at g > 1: coeff > 0 iff d < 0 or d > 2g/(g-1)
    law_ok = True
    for g in (Fr(1001, 1000), Fr(2), Fr(8), Fr(10 ** 6)):
        thr = 2 * g / (g - 1)
        cases = ((Fr(-5), 1), (Fr(-1, 10 ** 6), 1),
                 (thr / 2, -1), (Fr(1, 10), -1),
                 (thr * Fr(1001, 1000), 1), (thr, 0))
        for d, want in cases:
            c = coeff_factored(g, d)
            got = 0 if c == 0 else (1 if c > 0 else -1)
            if got != want:
                law_ok = False
    return {"poly_identity_exact": bool(ident),
            "tamper_control_detected": bool(control_detects),
            "n_exact_samples": n_exact,
            "n_mismatch": n_mismatch,
            "sign_law_exact_ok": bool(law_ok)}


# ================ claim CL2: delta < 0 half-plane ===================
def cl2():
    gs = [Fr(10 ** 9 + 1, 10 ** 9), Fr(3, 2), Fr(2), Fr(8),
          Fr(1000), Fr(10 ** 10)]
    ds = [Fr(-10 ** 6), Fr(-3), Fr(-2), Fr(-1), Fr(-3, 10),
          Fr(-1, 1000), Fr(-1, 10 ** 10), Fr(-1, 10 ** 30)]
    all_pos = True
    for g in gs:
        for d in ds:
            c = coeff_transcribed(g, d)
            if not (c > 0 and coeff_factored(g, d) == c):
                all_pos = False
    g = Fr(10 ** 10)
    cm = coeff_transcribed(g, Fr(-1, 10 ** 10))
    cp = coeff_transcribed(g, Fr(1, 10 ** 10))
    ok_corner = (cm > 0 and cp < 0
                 and cm == coeff_factored(g, Fr(-1, 10 ** 10))
                 and cp == coeff_factored(g, Fr(1, 10 ** 10)))
    # naive float64 fragility probe at the same corner
    fm = coeff_transcribed(1e10, -1e-10)
    fp = coeff_transcribed(1e10, 1e-10)
    return {"grid_all_positive_exact": bool(all_pos),
            "corner_neg_exact": str(cm),
            "corner_neg_float": float(cm),
            "corner_pos_float": float(cp),
            "corner_signs_exact_ok": bool(ok_corner),
            "float64_corner_neg": fm,
            "float64_corner_pos": fp,
            "float64_sign_agrees": bool(fm > 0 and fp < 0),
            "float64_rel_err_neg": abs(fm - float(cm))
            / abs(float(cm))}


# ================ claim CL3: linear in om^2, static min =============
def cl3():
    # exact linearity: second differences in om^2 vanish exactly
    rng = np.random.default_rng(22)
    lin_ok = True
    for _ in range(30):
        g = Fr(int(rng.integers(2, 900)), int(rng.integers(1, 47)))
        d = Fr(int(rng.integers(-900, 900)), int(rng.integers(1, 47)))
        r = Fr(int(rng.integers(1, 90)), int(rng.integers(1, 17)))
        if g == 1 or r == 0:
            continue
        a = Fr(int(rng.integers(0, 50)), 7)
        s = Fr(int(rng.integers(1, 50)), 11)
        h0 = Hm_frac(g, d, r, a)
        h1 = Hm_frac(g, d, r, a + s)
        h2 = Hm_frac(g, d, r, a + 2 * s)
        if h0 - 2 * h1 + h2 != 0:
            lin_ok = False
    # static identity as exact polynomials:
    # Hm numerator at om = 0 (w = -1) vs the fin numerator
    stat_num = pscale(pmul(D, pmul(D, padd(pscale(pmul(G, G), 2),
                                           pscale(B_POLY, -1)))), 4)
    fin_num = padd(padd(pscale(pmul(pmul(G, G), pmul(D, D)), 8),
                        pscale(pmul(pmul(GM1, G), pmul(D, pmul(D, D))),
                               8)),
                   pscale(pmul(pmul(GM1, GM1),
                               pmul(pmul(D, D), pmul(D, D))), -4))
    ident_static = stat_num == fin_num
    # the verified instance, exact
    h0 = Hm_frac(Fr(8), Fr(-3, 10), Fr(1), Fr(0))
    fi = fin_frac(Fr(8), Fr(-3, 10), Fr(1))
    # min over omega sits at 0 wherever the coefficient is positive:
    # affine in om^2 with slope coeff_transcribed (r-free), so for
    # g > 1, d < 0 the slope is positive (CL1 law) and argmin om = 0
    slope_pos = coeff_transcribed(Fr(8), Fr(-3, 10)) > 0
    grid = [Hm_frac(Fr(8), Fr(-3, 10), Fr(1), Fr(k, 10))
            for k in range(11)]
    monotone = all(grid[k + 1] > grid[k] for k in range(10))
    return {"exact_linear_om2": bool(lin_ok),
            "static_identity_poly": bool(ident_static),
            "instance_exact_equal": bool(h0 == fi),
            "instance_value": str(h0),
            "instance_float": float(h0),
            "slope_pos_at_8_m03": bool(slope_pos),
            "monotone_in_om2_grid": bool(monotone)}


# ================ pipeline probes (own machinery) ===================
def slope_E(g, dl, r=1.3, Delta=0.15, om_pair=(0.2, 0.4),
            nTH=61, nt=24):
    m = -mstar(g)
    o1, o2 = om_pair
    E1 = plane_E(m, g=g, dl=dl, om=o1, r=r, Delta=Delta,
                 nTH=nTH, nt=nt)
    E2 = plane_E(m, g=g, dl=dl, om=o2, r=r, Delta=Delta,
                 nTH=nTH, nt=nt)
    return (E2 - E1) / (o2 * o2 - o1 * o1), E1


def cl4():
    pts = [(8.0, -0.3), (8.0, -0.05), (3.0, -0.3), (100.0, -0.001)]
    claimed_pipe = [0.5522427721761849, 0.0023094867314125732,
                    0.7047960241503909, 1.598804894566053e-08]
    rows = []
    for (g, dl), cp in zip(pts, claimed_pipe):
        s, _ = slope_E(g, dl)
        cf = coeff_transcribed(g, dl)
        rows.append({"g": g, "dl": dl, "slope_mine": float(s),
                     "slope_claimed": cp,
                     "coeff_formula": float(cf),
                     "positive": bool(s > 0),
                     "rel_vs_claimed": abs(s - cp) / abs(cp)})
        print(f"  CL4 ({g:g},{dl:g}) slope {s:+.6e} "
              f"claimed {cp:+.6e}", flush=True)
    # negative controls at (8, -0.3): other omega probe pairs
    pairs = {}
    for pr in ((0.1, 0.5), (0.3, 0.6)):
        s, _ = slope_E(8.0, -0.3, om_pair=pr)
        pairs[f"{pr[0]:g}-{pr[1]:g}"] = float(s)
    # r-structure probe: E = S / r^2 + K om^2 ?
    Ks, Ss = {}, {}
    for r in (0.8, 1.0, 1.3, 2.0):
        s, E1 = slope_E(8.0, -0.3, r=r)
        Ks[f"{r:g}"] = float(s)
        Ss[f"{r:g}"] = float((E1 - 0.04 * s) * r * r)
    kv = list(Ks.values())
    sv = list(Ss.values())
    return {"rows": rows, "probe_pairs_slope": pairs,
            "K_of_r": Ks, "S_of_r": Ss,
            "K_rel_spread": (max(kv) - min(kv)) / abs(np.mean(kv)),
            "S_rel_spread": (max(sv) - min(sv)) / abs(np.mean(sv))}


def _simpE(Ha, TH, r):
    hh = TH[1] - TH[0]
    return 2.0 * float(hh / 3.0 * np.sum(_simpson_w(len(TH)) * Ha
                                         * np.cos(TH) * r * r))


def snapshot_E(m, g, dl, r, Delta, nTH=61, hfd=1e-5):
    """the naive om = 0 evaluator point: t = 0, om = 0."""
    TH = np.linspace(0.0, np.pi / 2 - Delta, nTH)
    X = r * np.cos(TH)
    Y = np.zeros_like(TH)
    Z = r * np.sin(TH)
    dM = fd_dM(np.zeros_like(TH), X, Y, Z, g, dl, m, 0.0, h=hfd)
    H = np.zeros_like(TH)
    for i in range(4):
        for j in range(4):
            F = dM[i] @ XI @ dM[j] - dM[j] @ XI @ dM[i]
            H += (F[..., 1, 2] ** 2 + F[..., 1, 3] ** 2
                  + F[..., 2, 3] ** 2 - F[..., 0, 1] ** 2
                  - F[..., 0, 2] ** 2 - F[..., 0, 3] ** 2)
    return _simpE(2.0 * H, TH, r)


def static_avg_E(m, g, dl, r, Delta, nTH=61, nal=24, hfd=1e-5):
    """own om -> 0 limit: clock-phase (alpha) average with zero
    time derivatives; the physically averaged static energy."""
    TH = np.linspace(0.0, np.pi / 2 - Delta, nTH)
    al = (np.arange(nal) + 0.5) * 2.0 * np.pi / nal
    A, THg = np.meshgrid(al, TH, indexing="ij")
    X = r * np.cos(THg)
    Y = np.zeros_like(THg)
    Z = r * np.sin(THg)
    dM = []
    for k in range(1, 4):
        v = [A, X, Y, Z]
        vp = list(v)
        vp[k] = vp[k] + hfd
        vm = list(v)
        vm[k] = vm[k] - hfd
        dM.append((build_M(*vp, g, dl, m, 1.0)
                   - build_M(*vm, g, dl, m, 1.0)) / (2.0 * hfd))
    H = np.zeros(A.shape)
    for i in range(3):
        for j in range(3):
            F = dM[i] @ XI @ dM[j] - dM[j] @ XI @ dM[i]
            H += (F[..., 1, 2] ** 2 + F[..., 1, 3] ** 2
                  + F[..., 2, 3] ** 2 - F[..., 0, 1] ** 2
                  - F[..., 0, 2] ** 2 - F[..., 0, 3] ** 2)
    Ha = 2.0 * H.mean(axis=0)
    return _simpE(Ha, TH, r)


def cl5():
    g, dl, r, Delta = 8.0, -0.3, 1.3, 0.15
    m = -mstar(g)
    oms = np.arange(1, 13) * 0.05
    Es = np.array([plane_E(m, g=g, dl=dl, om=o, r=r, Delta=Delta,
                           nTH=61, nt=24) for o in oms])
    Amat = np.vstack([np.ones_like(oms), oms * oms]).T
    (a0, a1), *_ = np.linalg.lstsq(Amat, Es, rcond=None)
    resid = np.max(np.abs(Es - Amat @ [a0, a1]) / Es)
    snap = snapshot_E(m, g, dl, r, Delta)
    stat = static_avg_E(m, g, dl, r, Delta)
    # small finite omega, full-period average: must sit ON the line
    E_tiny = plane_E(m, g=g, dl=dl, om=0.01, r=r, Delta=Delta,
                     nTH=61, nt=24)
    return {"slope_fit": float(a1), "intercept_fit": float(a0),
            "max_rel_resid": float(resid),
            "snapshot_t0": float(snap),
            "static_phase_avg": float(stat),
            "E_at_om_0p01": float(E_tiny),
            "snapshot_offset": float(snap - a0),
            "intercept_vs_static_rel": abs(a0 - stat)
            / abs(stat),
            "tiny_om_on_line_rel": abs(E_tiny
                                       - (a0 + a1 * 1e-4)) / E_tiny}


def _mscan(dl, om=0.3, r=1.0, Delta=0.15):
    def E(m):
        return plane_E(m, g=8.0, dl=dl, om=om, r=r, Delta=Delta,
                       nTH=41, nt=16)
    ms = np.linspace(-0.30, -0.02, 57)
    Es = np.array([E(m) for m in ms])
    i = int(np.argmin(Es))
    interior = bool(0 < i < len(ms) - 1)
    fine = np.linspace(ms[i] - 0.006, ms[i] + 0.006, 25)
    Ef = np.array([E(m) for m in fine])
    k = int(np.argmin(Ef))
    mf = float(fine[k])
    if 0 < k < len(fine) - 1:
        a, b, c = Ef[k - 1], Ef[k], Ef[k + 1]
        mf = float(fine[k] - 0.5 * (fine[1] - fine[0]) * (c - a)
                   / (c - 2.0 * b + a))
    return mf, interior, float(Es[0]), float(Es[-1]), float(np.min(Ef))


def cl6():
    m1, int1, eL1, eR1, _ = _mscan(-1e-4)
    m2, int2, eL2, eR2, emin2 = _mscan(-0.3)
    print(f"  CL6 dl=-1e-4 min {m1:+.5f}; dl=-0.3 min {m2:+.5f}",
          flush=True)
    # negative control: does the dl = -0.3 minimum depend on om?
    m2a, int2a, *_ = _mscan(-0.3, om=0.2)
    m2b, int2b, *_ = _mscan(-0.3, om=0.45)
    # evenness spot check at dl = -0.3
    Ep = plane_E(abs(m2), g=8.0, dl=-0.3, om=0.3, r=1.0,
                 Delta=0.15, nTH=41, nt=16)
    Em = plane_E(-abs(m2), g=8.0, dl=-0.3, om=0.3, r=1.0,
                 Delta=0.15, nTH=41, nt=16)
    return {"m_min_small_dl": m1, "interior_small": int1,
            "m_star_8": float(mstar(8.0)),
            "rel_dev_vs_mstar": abs(abs(m1) - mstar(8.0))
            / mstar(8.0),
            "m_min_big_dl": m2, "interior_big": int2,
            "E_ends_above_min_big": bool(eL2 > emin2
                                         and eR2 > emin2),
            "m_min_big_om0p2": m2a, "m_min_big_om0p45": m2b,
            "om_dependence_abs": max(abs(m2a - m2), abs(m2b - m2)),
            "evenness_rel": abs(Ep - Em) / abs(Ep)}


# ================ claim CL7: spectrum shift =========================
def build_M4(t, x, y, z, d4, m, om):
    """own generalized build, arbitrary diagonal d, pointwise."""
    G1 = np.zeros((4, 4))
    G1[2, 3] = -1.0
    G1[3, 2] = 1.0
    G2 = np.zeros((4, 4))
    G2[1, 3] = 1.0
    G2[3, 1] = -1.0
    G3 = np.zeros((4, 4))
    G3[1, 2] = -1.0
    G3[2, 1] = 1.0

    def rod(Gm, a):
        return (np.eye(4) + np.sin(a) * Gm
                + (1.0 - np.cos(a)) * (Gm @ Gm))

    phi = np.arctan2(y, x)
    th = -np.arctan2(z, np.sqrt(x * x + y * y))
    Qh = rod(G3, phi) @ rod(G2, th) @ rod(G1, om * t)
    R = np.sqrt(x * x + y * y + z * z)
    n = np.array([x, y, z]) / R
    K = np.zeros((4, 4))
    K[0, 1:] = n
    K[1:, 0] = n
    Qb = np.eye(4) + np.sinh(m) * K + (np.cosh(m) - 1.0) * (K @ K)
    Q = Qb @ Qh
    M = Q @ np.diag(d4) @ Q.T
    return 0.5 * (M + M.T), Q


def Hpt(t, x, y, z, d4, m, om, h=1e-5):
    """pointwise Hamiltonian density, own central FD."""
    v = [t, x, y, z]
    dM = []
    for k in range(4):
        vp = list(v)
        vp[k] += h
        vm = list(v)
        vm[k] -= h
        dM.append((build_M4(*vp, d4, m, om)[0]
                   - build_M4(*vm, d4, m, om)[0]) / (2.0 * h))
    H = 0.0
    for i in range(4):
        for j in range(4):
            F = dM[i] @ XI @ dM[j] - dM[j] @ XI @ dM[i]
            H += (F[1, 2] ** 2 + F[1, 3] ** 2 + F[2, 3] ** 2
                  - F[0, 1] ** 2 - F[0, 2] ** 2 - F[0, 3] ** 2)
    return float(H)


def cl7():
    # gate: build_M4 reproduces the gated Rodrigues build
    Mr, _ = build_M_rod(0.4, 1.1, 0.6, -0.8, 8.0, -0.3, 0.1, 0.3)
    M4, Q = build_M4(0.4, 1.1, 0.6, -0.8, (8.0, 1.0, -0.3, 0.0),
                     0.1, 0.3)
    gate = float(np.max(np.abs(Mr - M4)))
    # exact algebra: M(d + c I) - M(d) = c Q Q^T
    c = 0.7
    Ms, _ = build_M4(0.4, 1.1, 0.6, -0.8,
                     (8.0 + c, 1.0 + c, -0.3 + c, c), 0.1, 0.3)
    QQ = Q @ Q.T
    alg_dev = float(np.max(np.abs(Ms - M4 - c * QQ)))
    orth_defect = float(np.max(np.abs(QQ @ QQ.T - np.eye(4))))
    _, Q2 = build_M4(0.4, -0.5, 1.3, 0.9, (8.0, 1.0, -0.3, 0.0),
                     0.1, 0.3)
    pos_dep = float(np.max(np.abs(Q2 @ Q2.T - QQ)))
    # undressed branch: m = 0, random draws
    rng = np.random.default_rng(23)
    worst_m0 = 0.0
    for _ in range(5):
        g = rng.uniform(2.0, 10.0)
        dl = rng.uniform(-1.0, 1.0)
        cc = rng.uniform(-1.5, 1.5)
        p = rng.uniform(-1.5, 1.5, 4)
        if np.hypot(p[1], p[2]) < 0.4:
            continue
        h0 = Hpt(*p, (g, 1.0, dl, 0.0), 0.0, 0.3)
        h1 = Hpt(*p, (g + cc, 1.0 + cc, dl + cc, cc), 0.0, 0.3)
        worst_m0 = max(worst_m0, abs(h1 - h0)
                       / max(abs(h0), 1e-30))
    # boosted branch: fixed point and c, growing |m|
    grow = {}
    for m in (0.0013, 0.01, 0.05, 0.1):
        h0 = Hpt(0.4, 1.1, 0.6, -0.8, (8.0, 1.0, -0.3, 0.0), m, 0.3)
        h1 = Hpt(0.4, 1.1, 0.6, -0.8, (8.9, 1.9, 0.6, 0.9), m, 0.3)
        grow[f"{m:g}"] = abs(h1 - h0) / max(abs(h0), 1e-30)
    gv = list(grow.values())
    worst_boost = 0.0
    for _ in range(6):
        g = rng.uniform(2.0, 10.0)
        dl = rng.uniform(-1.0, 1.0)
        cc = rng.uniform(-1.5, 1.5)
        m = rng.uniform(0.01, 0.12) * rng.choice([-1.0, 1.0])
        p = rng.uniform(-1.5, 1.5, 4)
        if np.hypot(p[1], p[2]) < 0.4:
            continue
        h0 = Hpt(*p, (g, 1.0, dl, 0.0), m, 0.3)
        h1 = Hpt(*p, (g + cc, 1.0 + cc, dl + cc, cc), m, 0.3)
        worst_boost = max(worst_boost, abs(h1 - h0)
                          / max(abs(h0), 1e-30))
    return {"gate_M4_vs_rod": gate, "algebra_dev": alg_dev,
            "QQ_orth_defect_m0p1": orth_defect,
            "QQ_position_dep": pos_dep,
            "worst_rel_m0": worst_m0,
            "growth_fixed_point": grow,
            "growth_monotone": bool(all(gv[k + 1] > gv[k]
                                        for k in range(len(gv) - 1))),
            "worst_rel_boosted": worst_boost}


# ================ claim CL8: cone channel sign ======================
def cl8():
    m = -mstar(8.0)
    deltas = (0.4, 0.2, 0.1, 0.05, 0.02, 0.01)
    out = {}
    for dl in (0.3, -0.3):
        for r in (1.0, 1.3):
            lad = {}
            for Dd in deltas:
                lad[f"{Dd:g}"] = plane_E(m, g=8.0, dl=dl, om=0.3,
                                         r=r, Delta=Dd,
                                         nTH=161, nt=16)
            ls_main = (lad["0.05"] - lad["0.1"]) / np.log(2.0)
            ls_ext = (lad["0.01"] - lad["0.02"]) / np.log(2.0)
            out[f"dl{dl:+g}_r{r:g}"] = {
                "ladder": lad,
                "log_slope_0p1_0p05": float(ls_main),
                "log_slope_0p02_0p01": float(ls_ext),
                "rises_toward_small_Delta": bool(lad["0.05"]
                                                 > lad["0.4"])}
            print(f"  CL8 dl={dl:+g} r={r:g} E(0.4)="
                  f"{lad['0.4']:.4f} E(0.05)={lad['0.05']:.4f} "
                  f"E(0.01)={lad['0.01']:.4f}", flush=True)
    return out


# ================ main ==============================================
def main():
    os.makedirs(DATA, exist_ok=True)
    res = {}

    print("CL1 exact factorization ...", flush=True)
    c1 = cl1()
    v1 = (c1["poly_identity_exact"] and c1["tamper_control_detected"]
          and c1["n_mismatch"] == 0 and c1["sign_law_exact_ok"])
    res["CL1"] = {
        "claim": "om2_coeff factors exactly as "
                 "8 d^3 ((g-1) d - 2g)/(g-1); g>1 positive iff "
                 "d<0 or d>2g/(g-1)",
        "method": "own dict-polynomial algebra over exact integers, "
                  "factors multiplied out in code, denominators "
                  "cleared; 400 exact-Fraction random points; "
                  "tampered-coefficient comparator control; exact "
                  "sign-law cases incl the threshold d = 2g/(g-1)",
        "mine": c1,
        "claimed": {"factored": "8 d^3 ((g-1) d - 2g)/(g-1)",
                    "task_worst_rel": 1.1771868814335383e-12},
        "verdict": "CONFIRMED" if v1 else "REFUTED"}

    print("CL2 delta < 0 half-plane, exact corner ...", flush=True)
    c2 = cl2()
    v2 = (c2["grid_all_positive_exact"]
          and c2["corner_signs_exact_ok"]
          and abs(c2["corner_neg_float"] - 1.60000000024e-29)
          / 1.6e-29 < 1e-9)
    res["CL2"] = {
        "claim": "entire d<0 half-plane at g>1 omega-bounded; "
                 "corner (1e10, -1e-10) coeff +1.6e-29, "
                 "(1e10, +1e-10) -1.6e-29",
        "method": "exact fractions.Fraction on a 6x8 (g, d) grid "
                  "incl g = 1+1e-9 and d = -1e-30, plus the exact "
                  "corner; naive float64 route as fragility probe",
        "mine": c2,
        "claimed": {"corner_neg": 1.60000000024e-29,
                    "corner_pos": -1.6000000000800003e-29},
        "verdict": "CONFIRMED" if v2 else "REFUTED"}

    print("CL3 linear in om^2, static minimum ...", flush=True)
    c3 = cl3()
    v3 = (c3["exact_linear_om2"] and c3["static_identity_poly"]
          and c3["instance_exact_equal"]
          and abs(c3["instance_float"] - 0.6611510204081633) < 1e-15
          and c3["slope_pos_at_8_m03"]
          and c3["monotone_in_om2_grid"])
    res["CL3"] = {
        "claim": "Hm exactly linear in om^2; positive coeff gives "
                 "argmin om = 0 with value fin(g, d, r); instance "
                 "(8, -0.3, 1) both = 0.6611510204081633; "
                 "boundedness restored but state is STATIC",
        "method": "exact second differences in om^2 (30 random "
                  "Fraction draws); Hm(om=0, r=1) numerator == fin "
                  "numerator as exact polynomials; exact instance; "
                  "monotonicity of Hm on an om grid",
        "mine": c3,
        "claimed": {"both_equal": 0.6611510204081633,
                    "argmin_om": 0.0},
        "verdict": "CONFIRMED" if v3 else "REFUTED"}

    print("CL4 pipeline slope sign at negative delta ...", flush=True)
    c4 = cl4()
    v4 = (all(r["positive"] for r in c4["rows"])
          and all(r["rel_vs_claimed"] < 0.02 for r in c4["rows"])
          and all(s > 0 for s in c4["probe_pairs_slope"].values())
          and all(k > 0 for k in c4["K_of_r"].values()))
    res["CL4"] = {
        "claim": "ansatz-pipeline dE/d(om^2) positive at "
                 "(8,-0.3) ~ +0.552, (8,-0.05) ~ +0.00231, "
                 "(3,-0.3) ~ +0.705, (100,-0.001) ~ +1.60e-8 "
                 "(m = -|m*|, r = 1.3)",
        "method": "own pipeline (explicit-matrix build, plain "
                  "central FD, midpoint time avg nt=24, Simpson "
                  "nTH=61); probe pairs (0.2,0.4), (0.1,0.5), "
                  "(0.3,0.6); r-structure E = S/r^2 + K om^2 at "
                  "r = 0.8, 1.0, 1.3, 2.0",
        "mine": c4,
        "claimed": {"slopes": [0.5522427721761849,
                               0.0023094867314125732,
                               0.7047960241503909,
                               1.598804894566053e-08],
                    "coeffs": [0.5585142857142856,
                               0.002335714285714286,
                               0.7127999999999999,
                               1.616961616161616e-08]},
        "verdict": "CONFIRMED" if v4 else "REFUTED"}

    print("CL5 linearity + snapshot artifact ...", flush=True)
    c5 = cl5()
    v5 = (c5["max_rel_resid"] < 1e-8 and c5["slope_fit"] > 0
          and abs(c5["snapshot_t0"] - 0.5502061498648072) < 1e-3
          and c5["intercept_vs_static_rel"] < 1e-5
          and c5["tiny_om_on_line_rel"] < 1e-7
          and abs(c5["snapshot_offset"] + 0.264) < 2e-3)
    res["CL5"] = {
        "claim": "E linear in om^2 over (0.05..0.6] to ~1e-11 rel, "
                 "positive slope; om = 0 point off the line by "
                 "-0.264 = t = 0 snapshot artifact",
        "method": "own 12-point fit vs om^2; own t = 0 snapshot "
                  "evaluator; own static clock-phase-average "
                  "evaluator (alpha average, zero time derivs) as "
                  "the om -> 0 limit; E(om = 0.01) on-line check",
        "mine": c5,
        "claimed": {"resid": 8.968460342315011e-12,
                    "slope": 0.5522427722763334,
                    "intercept": 0.814528591438574,
                    "snapshot": 0.5502061498648072,
                    "offset": -0.26432244157376683},
        "verdict": "CONFIRMED" if v5 else "REFUTED"}

    print("CL6 m-sector survives negative delta ...", flush=True)
    c6 = cl6()
    v6 = (c6["interior_small"] and c6["rel_dev_vs_mstar"] < 0.01
          and c6["interior_big"]
          and abs(abs(c6["m_min_big_dl"]) - 0.1486) < 0.01
          and c6["E_ends_above_min_big"]
          and c6["evenness_rel"] < 1e-9)
    res["CL6"] = {
        "claim": "m-scan interior minimum survives d < 0: "
                 "|m| ~ 0.123 at (8, -1e-4, om 0.3) (resolution-"
                 "limited vs m* = 0.1257); |m| ~ 0.149 at (8, -0.3)",
        "method": "own scans, negative-m branch [-0.30, -0.02] step "
                  "0.005 + fine step 0.0005 + parabolic refine; om "
                  "controls 0.2 / 0.45 at dl = -0.3; evenness spot "
                  "check",
        "mine": c6,
        "claimed": {"m_min_small": -0.12312388633422496,
                    "m_star": 0.12565721414045308,
                    "m_min_big": -0.1486170209765923},
        "verdict": "CONFIRMED" if v6 else "REFUTED"}

    print("CL7 spectrum shift not exact when boosted ...", flush=True)
    c7 = cl7()
    v7 = (c7["gate_M4_vs_rod"] < 1e-12 and c7["algebra_dev"] < 1e-12
          and c7["QQ_orth_defect_m0p1"] > 1e-3
          and c7["QQ_position_dep"] > 1e-3
          and c7["worst_rel_m0"] < 1e-7
          and c7["growth_monotone"]
          and c7["worst_rel_boosted"] > 0.05)
    res["CL7"] = {
        "claim": "d -> d + c I shifts M by c Qb^2 (position-"
                 "dependent unless m = 0): H invariant at m = 0 "
                 "(~1e-10), broken at m != 0, violation grows with "
                 "|m| (up to ~0.41 at |m| ~ 0.1)",
        "method": "own arbitrary-diagonal builder (gated against "
                  "the Rodrigues build); exact algebra "
                  "M(d + cI) - M(d) = c Q Q^T; orthogonality defect "
                  "and position dependence of Q Q^T; pointwise H "
                  "with own FD, m = 0 draws vs growing |m|",
        "mine": c7,
        "claimed": {"undressed_worst": 4.275833557318291e-11,
                    "boosted_worst": 0.4102263682324696,
                    "small_m_rel": 7.510091779619473e-05},
        "verdict": "CONFIRMED" if v7 else "REFUTED"}

    print("CL8 cone channel flips sign ...", flush=True)
    c8 = cl8()
    neg = c8["dl-0.3_r1"]
    pos = c8["dl+0.3_r1"]
    lad_claim = {"0.4": 1.1009471878659454, "0.2": 1.3393166236184166,
                 "0.1": 1.544209968176653, "0.05": 1.7406159693389178}
    match = max(abs(neg["ladder"][k] - v) / v
                for k, v in lad_claim.items())
    v8 = (neg["rises_toward_small_Delta"]
          and neg["log_slope_0p1_0p05"] > 0.2
          and neg["log_slope_0p02_0p01"] > 0.2
          and not pos["rises_toward_small_Delta"]
          and pos["log_slope_0p1_0p05"] < 0
          and abs(pos["ladder"]["0.4"] - 0.686) < 0.01
          and abs(pos["ladder"]["0.05"] - 0.324) < 0.01
          and match < 5e-3
          and c8["dl-0.3_r1.3"]["rises_toward_small_Delta"]
          and not c8["dl+0.3_r1.3"]["rises_toward_small_Delta"])
    res["CL8"] = {
        "claim": "cone channel flips sign: dl = +0.3 E falls "
                 "0.686 -> 0.324 over Delta 0.4 -> 0.05 (collapse "
                 "channel); dl = -0.3 E rises 1.10 -> 1.74, log "
                 "slope ~ +0.28 (normal vortex core)",
        "method": "own Simpson quadrature nTH = 161, 6 Delta rungs "
                  "down to 0.01, both signs of dl, both r = 1.0 "
                  "and r = 1.3",
        "mine": c8,
        "claimed": {"neg_ladder": lad_claim,
                    "neg_log_slope": 0.2833539638776314,
                    "pos_ends": [0.686, 0.324]},
        "verdict": ("NUANCE" if v8 else "REFUTED"),
        "note": "numbers verified, but the claimed values sit at "
                "r = 1.0 while CL4/CL5 use r = 1.3 (hidden "
                "parameter switch; the sign flip itself holds at "
                "both r)"}

    with open(os.path.join(DATA, "m5_21_9_negdelta_audit.json"),
              "w") as f:
        json.dump(res, f, indent=1)

    print("\n=================== AUDIT SUMMARY ===================")
    for k in sorted(res):
        print(f"{k}: {res[k]['verdict']}")
    print("saved data/m5_21_9_negdelta_audit.json")
    return res


if __name__ == "__main__":
    main()
