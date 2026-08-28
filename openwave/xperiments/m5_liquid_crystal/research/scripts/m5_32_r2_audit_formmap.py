"""M5.32 R2 arm (a) INDEPENDENT ADVERSARIAL AUDIT of the form-map claims.

Built without reading the producer's m5_32_r2_a_formmap.py or
m5_32_r2_formmap.json. Imports only the audited registries
(m5_32_lagrangian.py: I1, I1_frob, V4, F_of_A; m5_32_terms_ext.py:
I1_h_np, J1_np, J2_np, Pgrad_np, timelike_eig_np) as ORACLES for the
densities; every form, split, channel, LP and path below is built here.

EQUATIONS FIRST
---------------
eta = diag(-1,1,1,1); M symmetric 4x4, M -> L M L^T; jets A_mu = d_mu M,
A_0 = omega a0 (clock direction a0, static background M).
    F_{0i} = a0 eta A_i - A_i eta a0            (linear in the spatial jet A_i)
    <F,G>_eta = tr(eta F eta G^T),  <F,G>_h = tr(h F h G^T),
    h = eta + 2 (eta u)(eta u)^T,  u = timelike unit eigenvector of M eta.
In the u-eigenframe (a Lorentz frame with u = e0) h = 1, so with the
entry split of the antisymmetric F_{0i} into its time-row part T_i
(entries (0,b), (b,0)) and its spatial-spatial part S_i (entries (a,b),
a,b >= 1):
    <F_0i,F_0i>_eta = S_i - T_i,   <F_0i,F_0i>_h = S_i + T_i,
    S_i = sum_{a,b>=1} F[a,b]^2 >= 0,   T_i = 2 sum_{b>=1} F[0,b]^2 >= 0.
Frame-free definitions used in the code (both Lorentz invariant):
    S := (<F,F>_h + <F,F>_eta)/2,   T := (<F,F>_h - <F,F>_eta)/2.
Lagrangian family and its Legendre energy (H = omega dL/domega - L, every
term exactly quadratic in omega, L = A + B omega + C omega^2 -> H = C
omega^2 - A):
    L_lambda = -4 [(1 - lambda) I1 + lambda I1_h] - V4,
    I1   = sum_{mu<nu} eta^mu eta^nu <F_mu nu, F_mu nu>_eta
         = -sum_i <F_0i,F_0i>_eta + sum_{i<j} <F_ij,F_ij>_eta,
    I1_h the same with <,>_h,
    H2(lambda) := omega^2 coefficient of the energy
               = 4 sum_i [(1 - lambda)(S_i - T_i) + lambda (S_i + T_i)]
               = 4 [S - (1 - 2 lambda) T],   S = sum_i S_i,  T = sum_i T_i.
    => affine in lambda, slope 8 T, H2(1/2) = 4 S >= 0.
Threshold per channel: H2(lambda) PSD  <=>  S(x) >= (1 - 2 lambda) T(x)
for every jet x  <=>  lambda >= (1 - mu)/2,  mu := inf_{T(x)>0} S(x)/T(x).
lambda* = 1/2 exactly on a channel iff mu = 0, i.e. iff a jet with S = 0,
T > 0 exists (pure time-row curvature).
Explicit pure-time-row jets at M_vac = diag(g, 1, delta, 0) (s = -1):
    boost_1 tangent a0 = (1+g)(e0 e1^T + e1 e0^T):
        F_0i = X - X^T, X = (1+g)(e0 e1^T - e1 e0^T) A_i, so
        spatial block = 0 iff A_i[0,2] = A_i[0,3] = 0 and
        F_0i[0,1] = (1+g)(A_i[1,1] + A_i[0,0]);  jet A_1 = e1 e1^T works.
    rot_3 tangent a0 = (delta - 1)(e1 e2^T + e2 e1^T):
        spatial block = 0 iff A[1,1] = A[2,2], A[1,3] = A[2,3] = 0 and
        F_0i[0,1] = (1 - delta) A_i[0,2], F_0i[0,2] = (1 - delta) A_i[0,1];
        jet A_1 = e0 e1^T + e1 e0^T works.
Degeneracy path (F3): M(t) = M_vac + t (e0 e1^T + e1 e0^T), N = M eta has
the block [[-g, t], [-t, 1]], eigenvalues collide at t* = (g + 1)/2 and
u (normalized u^T eta u = -1) has u^0 -> infinity there (the eigenvectors
tend to a null vector). V4(t) = w sum_p (tr N^p - C_p)^2, C_p = (-g)^p + 1
+ delta^p, w = registry W1. The minimal V4 on the 2-block degeneracy
locus: eigenvalue pair (m, m), V4_min = w min_m sum_p (2 m^p + rest_p - C_p)^2.
Pgrad (F4): d_0 u = omega du(a0), omega^2 Lagrangian coefficient
C_P = eta^00 q(dP, dP) = -q(dP, dP), dP = du (eta u)^T + u (eta du)^T,
q = sum eta_a eta_b X_ab^2; on a boost tangent du = e_k so q = -2, C_P = +2;
energy contribution with L containing c_P Pgrad is H2 containing c_P C_P = 2 c_P (c_P = -kappa).
LP (F5): H2_total(c) = -4 C_I1 + sum_k c_k C_k as 30x30 forms on every
channel; sampled constraints x^T H2 x >= 0 are NECESSARY, so an
infeasible LP is a rigorous certificate; feasible answers are refined
with cutting planes on the minimal eigenvector.

Conventions chosen here: s = -1 branch, M_vac = diag(g, 1, delta, 0);
jets = general symmetric A_i (30-dim, 10 entries per i); Q matrices by
exact polarization; u sign irrelevant (only u u^T enters).

Out: ../data/m5_32_r2_audit_formmap.json
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys
import time

import numpy as np
from scipy.linalg import expm
from scipy.optimize import linprog, minimize_scalar

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT = os.path.join(DATA, "m5_32_r2_audit_formmap.json")


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


L0 = _load("m5_32_lagrangian", "m5_32_lagrangian.py")
EXT = _load("m5_32_terms_ext", "m5_32_terms_ext.py")
ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
W1 = float(L0.W1)
T0 = time.time()
RES = {"conventions": {
    "eta": "diag(-1,1,1,1)", "branch": "s = -1, M_vac = diag(g,1,delta,0)",
    "jet": "general symmetric A_i, 30 real parameters",
    "H2": "omega^2 coefficient of the Legendre energy H = C omega^2 - A of "
          "L_lambda = -4[(1-lambda) I1 + lambda I1_h] - V4",
    "S_T": "S = (<F,F>_h + <F,F>_eta)/2, T = (<F,F>_h - <F,F>_eta)/2, summed over i",
    "w": W1, "toy": {"g": 32.0, "delta": 0.3}, "alt": {"g": 8.0, "delta": 0.3}},
    "claims": {}}
LOG = []


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.append(s)


# ---------------- building blocks (all my own) ----------------
def vac(g, delta):
    return np.diag([g, 1.0, delta, 0.0])


def gen_boost(k):
    G = np.zeros((4, 4)); G[0, k] = G[k, 0] = 1.0; return G


def gen_rot(i, j):
    G = np.zeros((4, 4)); G[i, j] = -1.0; G[j, i] = 1.0; return G


def gamma(t, r):
    Gm = np.zeros((4, 4))
    Gm[0, 1:] = t; Gm[1:, 0] = t
    Gm[1, 2], Gm[1, 3] = -r[2], r[1]
    Gm[2, 1], Gm[2, 3] = r[2], -r[0]
    Gm[3, 1], Gm[3, 2] = -r[1], r[0]
    return Gm


def coms(A, B):
    return A @ ETA @ B - B @ ETA @ A


def tangent(G, M):
    return G @ M + M @ G.T


def timelike_u(M):
    lam, V = np.linalg.eig(M @ ETA)
    if np.max(np.abs(lam.imag)) > 1e-9 * max(1.0, np.max(np.abs(lam))):
        return None, lam
    lam, V = lam.real, V.real
    n2 = np.einsum("ak,a,ak->k", V, np.diag(ETA), V)
    k = np.argmin(n2)
    if n2[k] >= 0:
        return None, lam
    return V[:, k] / np.sqrt(-n2[k]), lam


def h_of(M):
    u, _ = timelike_u(M)
    hu = ETA @ u
    return ETA + 2.0 * np.outer(hu, hu)


def F0(a0, Ai):
    return a0 @ ETA @ Ai - Ai @ ETA @ a0


def ip_eta(F):
    return np.trace(ETA @ F @ ETA @ F.T)


def ip_h(F, h):
    return np.trace(h @ F @ h @ F.T)


IU = np.triu_indices(4)


def jet30(x):
    A = np.zeros((3, 4, 4))
    for i in range(3):
        S = np.zeros((4, 4)); S[IU] = x[10 * i:10 * i + 10]
        A[i] = S + S.T - np.diag(np.diag(S))
    return A


def ST_of(a0, M, x):
    """(S, T) summed over i, frame-free definition."""
    h = h_of(M); A = jet30(x); S = T = 0.0
    for i in range(3):
        F = F0(a0, A[i]); e, hh = ip_eta(F), ip_h(F, h)
        S += 0.5 * (hh + e); T += 0.5 * (hh - e)
    return S, T


def ST_entries_vacframe(a0, M, x):
    """entry split (valid only when u = e0): S = spatial block, T = time row
    + time column (T = 2 sum_b F[0,b]^2 when F is antisymmetric, i.e. when
    a0 is symmetric; the coms(Gamma_0, M) clocks with a boost part are NOT
    symmetric, so F_0i is not antisymmetric there and both entries count)."""
    A = jet30(x); S = T = 0.0
    for i in range(3):
        F = F0(a0, A[i])
        # (0,0) carries eta_0 eta_0 = +1: it sits with S (zero when F is antisymmetric)
        S += np.sum(F[1:, 1:] ** 2) + F[0, 0] ** 2
        T += np.sum(F[0, 1:] ** 2) + np.sum(F[1:, 0] ** 2)
    return S, T


def polarize(fn, n=30):
    """quadratic form matrix of a scalar function fn(x) homogeneous quadratic."""
    Q = np.zeros((n, n)); d = np.zeros(n)
    for i in range(n):
        x = np.zeros(n); x[i] = 1.0; d[i] = fn(x); Q[i, i] = d[i]
    for i in range(n):
        for j in range(i + 1, n):
            x = np.zeros(n); x[i] = x[j] = 1.0
            Q[i, j] = Q[j, i] = 0.5 * (fn(x) - d[i] - d[j])
    return Q


def QST(a0, M):
    QS = polarize(lambda x: ST_of(a0, M, x)[0])
    QT = polarize(lambda x: ST_of(a0, M, x)[1])
    return QS, QT


def mineig(Q):
    return float(np.linalg.eigvalsh(0.5 * (Q + Q.T))[0])


def H2_form(QS, QT, lam):
    return 4.0 * (QS - (1.0 - 2.0 * lam) * QT)


def threshold(QS, QT, tol=1e-10):
    """smallest lambda in [0,1] with H2(lambda) PSD (bisection; monotone)."""
    scale = max(1.0, abs(np.linalg.eigvalsh(QT)[-1]), abs(np.linalg.eigvalsh(QS)[-1]))
    if mineig(H2_form(QS, QT, 0.0)) >= -tol * scale:
        return 0.0
    lo, hi = 0.0, 1.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if mineig(H2_form(QS, QT, mid)) >= -tol * scale:
            hi = mid
        else:
            lo = mid
    return hi


def mu_of(QS, QT, tol=1e-9):
    """inf S/T over T > 0: Schur complement over null(Q_T), then generalized eig."""
    w, V = np.linalg.eigh(QT)
    keep = w > tol * max(1.0, w[-1])
    Vr, Vn = V[:, keep], V[:, ~keep]
    Srr = Vr.T @ QS @ Vr; Srn = Vr.T @ QS @ Vn; Snn = Vn.T @ QS @ Vn
    Seff = Srr - Srn @ np.linalg.pinv(Snn, rcond=1e-10) @ Srn.T
    D = np.diag(1.0 / np.sqrt(w[keep]))
    return float(np.linalg.eigvalsh(D @ Seff @ D)[0]), int((~keep).sum())


def pure_time_row_jet(QS, QT, tol=1e-9):
    """max T on null(Q_S): returns (T_max, jet) or (0, None)."""
    w, V = np.linalg.eigh(QS)
    Vn = V[:, w < tol * max(1.0, w[-1])]
    if Vn.shape[1] == 0:
        return 0.0, None
    wt, Vt = np.linalg.eigh(Vn.T @ QT @ Vn)
    return float(wt[-1]), Vn @ Vt[:, -1]


def density_H2(term_density, a0, M, x, p):
    """omega^2 coefficient of a registry density (oracle) on the jet x."""
    def dens(w):
        A = np.zeros((4, 4, 4)); A[1:] = jet30(x); A[0] = w * a0
        return float(term_density(A, M, p))
    return 0.5 * (dens(1.0) + dens(-1.0)) - dens(0.0)


def lorentz(rng, scale=0.3):
    G = np.zeros((4, 4))
    for k in (1, 2, 3):
        G += rng.normal() * scale * gen_boost(k)
    G += rng.normal() * scale * gen_rot(1, 2) + rng.normal() * scale * gen_rot(2, 3) \
        + rng.normal() * scale * gen_rot(1, 3)
    return expm(G)


def channel_set(g, delta, rng, extended=True):
    M = vac(g, delta)
    ch = {}
    for k in (1, 2, 3):
        ch[f"boost_{k}"] = (tangent(gen_boost(k), M), M)
    for (i, j) in ((1, 2), (2, 3), (1, 3)):
        ch[f"rot_{i}{j}"] = (tangent(gen_rot(i, j), M), M)
    t = np.array([1.0, 0.7, -0.4]); r = np.array([0.5, -0.8, 0.3])
    ch["clock_tr"] = (coms(gamma(t, r), M), M)
    ch["clock_t"] = (coms(gamma(t, 0 * r), M), M)
    ch["clock_r"] = (coms(gamma(0 * t, r), M), M)
    if extended:
        # NEW 1: boost along a random direction (twisted about a different axis)
        n = rng.normal(size=3); n /= np.linalg.norm(n)
        Gb = sum(n[k - 1] * gen_boost(k) for k in (1, 2, 3))
        ch["boost_twisted"] = (tangent(Gb, M), M)
        # NEW 2: rotation-then-boost tangent (rotation axis != boost axis)
        R = expm(0.8 * gen_rot(2, 3))
        ch["boost_1_rot23"] = (tangent(R @ gen_boost(1) @ R.T, M), M)
        # NEW 3: boost tangent on a NON-vacuum static background (random symmetric)
        dM = rng.normal(size=(4, 4)); dM = 0.05 * (dM + dM.T) * g / 32.0
        Mb = M + dM
        ch["boost_1_nonvac_bg"] = (tangent(gen_boost(1), Mb), Mb)
        ch["clock_random_nonvac_bg"] = (coms(gamma(rng.normal(size=3), rng.normal(size=3)), Mb), Mb)
        # NEW 4: local clocks with random (t, r) at the vacuum
        for q in range(3):
            ch[f"clock_random_{q}"] = (coms(gamma(rng.normal(size=3), rng.normal(size=3)), M), M)
        # NEW 5: boost tangent on a Lorentz-dressed vacuum (u != e0)
        Ld = lorentz(rng, 0.5); Md = Ld @ M @ Ld.T
        ch["boost_1_dressed_bg"] = (tangent(gen_boost(1), Md), Md)
        # NEW 6: a generic symmetric a0 (not a tangent of anything)
        a = rng.normal(size=(4, 4)); ch["a0_random_symmetric"] = (a + a.T, M)
    return ch


# =============================== F1 ===============================
def audit_F1(rng):
    g, delta = 32.0, 0.3
    M = vac(g, delta)
    p = L0.default_params(s=-1.0, g=g, delta=delta)
    I1 = L0.REGISTRY["I1"]
    I1h = EXT.I1_h_np
    ch = channel_set(g, delta, rng, extended=False)
    names = ["boost_1", "boost_2", "rot_12", "rot_23", "clock_tr"]
    lams = np.array([0.0, 0.25, 0.5, 0.75, 1.0, 1.3])
    worst_affine = worst_slope = worst_half = worst_split = 0.0
    for nm in names:
        a0, Mb = ch[nm]
        for _ in range(20):
            x = rng.normal(size=30)
            S, T = ST_of(a0, Mb, x)
            Se, Te = ST_entries_vacframe(a0, Mb, x)
            worst_split = max(worst_split, abs(S - Se) + abs(T - Te))
            C1 = density_H2(I1.density, a0, Mb, x, p)
            Ch = density_H2(I1h, a0, Mb, x, p)
            H2 = -4.0 * ((1 - lams) * C1 + lams * Ch)      # V4 has no omega
            pred = 4.0 * (S - (1 - 2 * lams) * T)
            worst_affine = max(worst_affine, np.max(np.abs(H2 - pred)) / max(1.0, np.max(np.abs(pred))))
            slope = np.polyfit(lams, H2, 1)[0]
            worst_slope = max(worst_slope, abs(slope - 8 * T) / max(1.0, abs(8 * T)))
            H2half = -4.0 * (0.5 * C1 + 0.5 * Ch)
            worst_half = max(worst_half, abs(H2half - 4 * S) / max(1.0, abs(4 * S)))
    asym = {nm: float(np.linalg.norm(ch[nm][0] - ch[nm][0].T) / max(1e-300, np.linalg.norm(ch[nm][0])))
            for nm in ch}
    out = {"channels": names, "n_jets_per_channel": 20, "lambda_grid": lams.tolist(),
           "a0_asymmetry_|a0-a0T|/|a0|_per_channel": asym,
           "note_a0_symmetry": "coms(Gamma_0, M) with a boost part t is NOT symmetric (clock_t, clock_tr): "
                               "not a tangent of a symmetric M; the frame-free S/T split still holds there, "
                               "the antisymmetric-F entry formula T = 2 sum_b F[0,b]^2 does not",
           "max_rel_dev_H2_vs_4[S-(1-2lam)T]": worst_affine,
           "max_rel_dev_slope_vs_8T": worst_slope,
           "max_rel_dev_H2(1/2)_vs_4S": worst_half,
           "max_abs_dev_framefree_vs_entry_split": worst_split}
    # explicit pure time-row jets (analytic) verified numerically
    ex = {}
    x = np.zeros(30); x[[i for i in range(10) if (IU[0][i], IU[1][i]) == (1, 1)][0]] = 1.0
    S, T = ST_of(ch["boost_1"][0], M, x)
    Fb = F0(ch["boost_1"][0], jet30(x)[0])
    ex["boost_1: A_1 = e1 e1^T"] = {"S": S, "T": T, "F_01": Fb.tolist(),
                                    "F_01[0,1] predicted (1+g)": 1 + g}
    x = np.zeros(30); x[[i for i in range(10) if (IU[0][i], IU[1][i]) == (0, 1)][0]] = 1.0
    S, T = ST_of(ch["rot_12"][0], M, x)
    Fr = F0(ch["rot_12"][0], jet30(x)[0])
    ex["rot_12: A_1 = e0 e1^T + e1 e0^T"] = {"S": S, "T": T, "F_01": Fr.tolist(),
                                             "F_01[0,2] predicted (1-delta)": 1 - delta}
    # and the numerical null(Q_S) search on every R1 channel
    nulls = {}
    for nm, (a0, Mb) in ch.items():
        QS, QT = QST(a0, Mb)
        Tm, v = pure_time_row_jet(QS, QT)
        nulls[nm] = {"dim_null_QS": int(np.sum(np.linalg.eigvalsh(QS) < 1e-9 * max(1, np.linalg.eigvalsh(QS)[-1]))),
                     "max_T_on_null_QS": Tm,
                     "min_eig_QS": mineig(QS), "min_eig_QT": mineig(QT)}
    out["explicit_pure_time_row_jets"] = ex
    out["null_QS_search"] = nulls
    ok = (worst_affine < 1e-9 and worst_slope < 1e-9 and worst_half < 1e-9 and worst_split < 1e-9
          and ex["boost_1: A_1 = e1 e1^T"]["S"] < 1e-12 and ex["boost_1: A_1 = e1 e1^T"]["T"] > 0
          and ex["rot_12: A_1 = e0 e1^T + e1 e0^T"]["S"] < 1e-12 and ex["rot_12: A_1 = e0 e1^T + e1 e0^T"]["T"] > 0
          and all(v["max_T_on_null_QS"] > 1e-6 for v in nulls.values()))
    out["verdict"] = "CONFIRMED" if ok else "REFUTED"
    log(f"F1 a0 asymmetry per channel: { {k: round(v, 3) for k, v in asym.items()} }")
    log(f"F1 {out['verdict']}: affine dev {worst_affine:.1e}, slope dev {worst_slope:.1e}, "
        f"H2(1/2)=4S dev {worst_half:.1e}; explicit jets: boost T={ex['boost_1: A_1 = e1 e1^T']['T']:.4g}, "
        f"rot T={ex['rot_12: A_1 = e0 e1^T + e1 e0^T']['T']:.4g}")
    for nm, v in nulls.items():
        log(f"   {nm}: dim null(Q_S)={v['dim_null_QS']}, max T on it={v['max_T_on_null_QS']:.4g}")
    return out


# =============================== F2 ===============================
def audit_F2(rng):
    out = {}
    for (g, delta) in ((32.0, 0.3), (8.0, 0.3)):
        ch = channel_set(g, delta, rng, extended=True)
        rows = {}
        for nm, (a0, Mb) in ch.items():
            QS, QT = QST(a0, Mb)
            mu, dimnullT = mu_of(QS, QT)
            lam_star = threshold(QS, QT)
            probe = {str(l): mineig(H2_form(QS, QT, l)) for l in (0.45, 0.5, 0.55, 0.6, 1.0)}
            rows[nm] = {"mu": mu, "lambda_star": lam_star, "dim_null_QT": dimnullT,
                        "min_eig_QS": mineig(QS), "min_eig_QT": mineig(QT),
                        "min_eig_H2_at": probe,
                        "predicted_lambda_star_(1-mu)/2": 0.5 * (1 - max(mu, 0.0))}
            log(f"F2 g={g:g}: {nm:26s} mu={mu:+.3e} lambda*={lam_star:.6f} "
                f"mineig H2(0.55)={probe['0.55']:+.2e} H2(0.6)={probe['0.6']:+.2e}")
        out[f"g={g:g}"] = rows
    # homogeneity: the tangent a0 does not depend on the jet (it is a function of M only)
    g, delta = 32.0, 0.3; M = vac(g, delta)
    hom = {}
    for nm, (a0, Mb) in channel_set(g, delta, rng, extended=False).items():
        # a0 is built before any jet is drawn; check H2 is exactly homogeneous degree 2
        x = rng.normal(size=30); S1, T1 = ST_of(a0, Mb, x); S2, T2 = ST_of(a0, Mb, 3.0 * x)
        hom[nm] = {"S(3x)/S(x)": S2 / S1, "T(3x)/T(x)": T2 / T1}
    out["homogeneity_tangent_channels"] = hom
    # verdict
    allrows = [r for gg in ("g=32", "g=8") for r in out[gg].values()]
    lowered = [nm for gg in ("g=32", "g=8") for nm, r in out[gg].items() if r["lambda_star"] < 0.5 - 1e-6]
    raised = [nm for gg in ("g=32", "g=8") for nm, r in out[gg].items() if r["lambda_star"] > 0.5 + 1e-6]
    neg_above = [nm for gg in ("g=32", "g=8") for nm, r in out[gg].items()
                 if min(r["min_eig_H2_at"]["0.55"], r["min_eig_H2_at"]["0.6"]) < -1e-8]
    out["channels_lowering_threshold"] = lowered
    out["channels_raising_threshold"] = raised
    out["channels_negative_above_half"] = neg_above
    base = ["boost_1", "boost_2", "boost_3", "rot_12", "rot_23", "rot_13", "clock_tr", "clock_t", "clock_r"]
    base_ok = all(abs(out[gg][nm]["lambda_star"] - 0.5) < 1e-6 for gg in ("g=32", "g=8") for nm in base
                  if nm.startswith(("boost", "rot")))
    clocks_mu = {gg: {nm: out[gg][nm]["mu"] for nm in base if nm.startswith("clock")} for gg in ("g=32", "g=8")}
    out["clock_mu"] = clocks_mu
    hom_ok = all(abs(v["S(3x)/S(x)"] - 9) < 1e-9 and abs(v["T(3x)/T(x)"] - 9) < 1e-9 for v in hom.values())
    if base_ok and hom_ok and not raised and not neg_above:
        out["verdict"] = "CONFIRMED" if not lowered else "QUALIFIED"
    else:
        out["verdict"] = "REFUTED"
    log(f"F2 {out['verdict']}: lowered={lowered} raised={raised} negative above 1/2={neg_above}")
    return out


# =============================== F3 ===============================
def V4_density(M, g, delta):
    N = M @ ETA; P = np.eye(4); tot = 0.0
    for k in range(1, 5):
        P = P @ N
        Ck = (-g) ** k + 1.0 + delta ** k
        tot += (np.trace(P) - Ck) ** 2
    return W1 * tot


def audit_F3(rng):
    g, delta = 32.0, 0.3
    M = vac(g, delta)
    p = L0.default_params(s=-1.0, g=g, delta=delta)
    # (a) Lorentz dressings keep the spectrum + timelike u with gap g
    spec_dev = 0.0; gaps = []; u_ok = True
    vac_spec = np.sort(np.linalg.eigvals(M @ ETA).real)
    for _ in range(50):
        Ld = lorentz(rng, 0.7); Md = Ld @ M @ Ld.T
        u, lam = timelike_u(Md)
        if u is None:
            u_ok = False; continue
        lam = np.sort(lam.real)
        spec_dev = max(spec_dev, np.max(np.abs(lam - vac_spec)) / g)
        l0 = float((ETA @ u) @ (Md @ ETA @ u) / ((ETA @ u) @ u))
        gaps.append(float(np.min(np.abs(np.delete(lam, np.argmin(np.abs(lam - l0))) - l0))))
        # V4 exactly zero on the dressing
        spec_dev = max(spec_dev, V4_density(Md, g, delta))
    out = {"dressings": {"n": 50, "max_rel_spectrum_dev_and_V4": spec_dev,
                         "min_gap": float(np.min(gaps)), "max_gap": float(np.max(gaps)),
                         "u_exists_always": u_ok, "expected_gap": g}}
    # (b) the plain degeneracy path
    tstar = 0.5 * (g + 1)
    ts = [0.0, 4.0, 8.0, 12.0, 15.0, 16.0, 16.4, 16.49, 16.499, tstar, 17.0, 20.0]
    x = rng.normal(size=30); A = jet30(x)
    path = []
    for t in ts:
        Mt = M + t * (np.outer(E(0), E(1)) + np.outer(E(1), E(0)))
        u, lam = timelike_u(Mt)
        row = {"t": t, "V4": V4_density(Mt, g, delta), "spectrum_real": bool(u is not None)}
        if u is not None:
            h = h_of(Mt)
            row["u0_gamma"] = float(abs(u[0]))
            row["gap"] = float(np.sort(np.abs(np.diff(np.sort(lam.real))))[0])
            # static I1_h energy density on the FIXED random jet (h = I in the u frame)
            Ust = 0.0; Uet = 0.0
            for i in range(3):
                for j in range(i + 1, 3):
                    Fij = A[i] @ ETA @ A[j] - A[j] @ ETA @ A[i]
                    Ust += ip_h(Fij, h); Uet += ip_eta(Fij)
            row["4*I1_h_static_on_fixed_jet"] = 4 * Ust
            row["4*I1_static_on_fixed_jet"] = 4 * Uet
        path.append(row)
        log(f"F3 path t={t:7.3f}: V4={row['V4']:.4g} real={row['spectrum_real']} "
            + (f"gamma(u)={row['u0_gamma']:.4g} 4 I1_h={row['4*I1_h_static_on_fixed_jet']:.4g}" if u is not None else ""))
    V4star = V4_density(M + tstar * (np.outer(E(0), E(1)) + np.outer(E(1), E(0))), g, delta)
    # (c) minimal V4 on the 2-block degeneracy locus (eigenvalue pair (m, m))
    C = [(-g) ** k + 1.0 + delta ** k for k in range(1, 5)]
    fam = {}
    for nm, rest in (("(-g) merges with 1", [delta ** k + 0.0 for k in range(1, 5)]),
                     ("(-g) merges with delta", [1.0 + 0.0 for k in range(1, 5)]),
                     ("(-g) merges with 0", [1.0 + delta ** k for k in range(1, 5)])):
        f = lambda m, rest=rest: W1 * sum((2 * m ** k + rest[k - 1] - C[k - 1]) ** 2 for k in range(1, 5))
        best = None
        for m0 in np.linspace(-40, 40, 161):
            r = minimize_scalar(f, bracket=(m0 - 0.5, m0 + 0.5))
            if best is None or r.fun < best.fun:
                best = r
        fam[nm] = {"m": float(best.x), "V4_min": float(best.fun)}
    V4_locus_min = min(v["V4_min"] for v in fam.values())
    # independent energy scales
    row = json.load(open(os.path.join(DATA, L0.STORED3_ROW)))
    E_u_stored = row["E_u"]
    scale_producer = 1e4
    out.update({"path": path, "t_star": tstar, "V4_at_t_star": V4star,
                "V4_min_on_2block_degeneracy_locus": fam, "V4_locus_min": V4_locus_min,
                "scales": {"producer_stated_I1_h_dressed_electron": scale_producer,
                           "stored_electron_E_u_(g=8,n=32,undressed,integrated)": E_u_stored,
                           "note": "V4 here is a per-cell DENSITY; the 1e4 scale is an INTEGRATED energy: "
                                   "the ratio mixes units unless multiplied by the relaxing volume (h^3 per cell)"},
                "ratio_V4star_over_1e4": V4star / scale_producer,
                "ratio_V4locusmin_over_1e4": V4_locus_min / scale_producer})
    # verdict on the implicit protection claim
    grows = all(path[k]["4*I1_h_static_on_fixed_jet"] <= path[k + 1]["4*I1_h_static_on_fixed_jet"] + 1e-9
                for k in range(3, 8))
    out["I1_h_diverges_toward_locus"] = grows
    out["verdict"] = "QUALIFIED"
    log(f"F3 QUALIFIED: V4(t*={tstar}) = {V4star:.4g} (ratio to 1e4: {V4star / 1e4:.3g}); "
        f"min V4 on the 2-block locus = {V4_locus_min:.4g}; I1_h on a fixed jet grows toward the locus: {grows}")
    return out


def E(i):
    e = np.zeros(4); e[i] = 1.0; return e


# =============================== F4 ===============================
def audit_F4(rng):
    g, delta = 32.0, 0.3
    M = vac(g, delta)
    p = L0.default_params(s=-1.0, g=g, delta=delta)
    ch = channel_set(g, delta, rng, extended=False)
    out = {}
    for nm in ("boost_1", "boost_3", "rot_12", "rot_23", "clock_tr"):
        a0, Mb = ch[nm]
        Q = polarize(lambda x: density_H2(EXT.Pgrad_np, a0, Mb, x, p) - density_H2(EXT.Pgrad_np, a0, Mb, 0 * x, p))
        c0 = density_H2(EXT.Pgrad_np, a0, Mb, np.zeros(30), p)
        # my own analytic value: du = G u for a tangent; q = sum eta_a eta_b dP^2; C_P = -q
        out[nm] = {"C_P_Lagrangian_omega2_coeff_at_zero_jet": c0, "norm_of_jet_dependent_form": float(np.linalg.norm(Q)),
                   "energy_H2_with_L=-kappa*Pgrad": f"-kappa * {c0:.6g}"}
        log(f"F4 {nm}: C_P = {c0:+.6g} (constant), ||jet form|| = {np.linalg.norm(Q):.1e}")
    # explicit sign check on boost_1 by hand
    u = E(0); du = E(1)
    dP = np.outer(du, ETA @ u) + np.outer(u, ETA @ du)
    q = float(np.einsum("a,b,ab,ab->", np.diag(ETA), np.diag(ETA), dP, dP))
    out["hand_boost_1"] = {"q(dP,dP)": q, "C_P = eta^00 q": -q}
    ok = (abs(out["boost_1"]["C_P_Lagrangian_omega2_coeff_at_zero_jet"] - 2.0) < 1e-9
          and abs(out["rot_12"]["C_P_Lagrangian_omega2_coeff_at_zero_jet"]) < 1e-9
          and all(v["norm_of_jet_dependent_form"] < 1e-9 for k, v in out.items() if k != "hand_boost_1"))
    out["sign_reading"] = ("C_P = +2 on boosts is the LAGRANGIAN omega^2 coefficient; the energy is H2 = c_P C_P, "
                           "so with L = ... - kappa Pgrad (c_P = -kappa) the boost-channel energy gets -2 kappa omega^2: "
                           "kappa > 0 REMOVES omega^2 energy, kappa < 0 adds a constant that cannot lift a jet-direction "
                           "negative eigenvalue (the form is homogeneous in the jet, the constant is not)")
    out["verdict"] = "CONFIRMED" if ok else "REFUTED"
    log(f"F4 {out['verdict']}: {out['sign_reading']}")
    return out


# =============================== F5 ===============================
def lp_feasible(Qs, free, fixed, rng, nsamp=2000, objective=None):
    rows, rhs = [], []
    for nm, Q in Qs.items():
        n = next(iter(Q.values())).shape[0]
        X = rng.normal(size=(nsamp, n)); X /= np.linalg.norm(X, axis=1, keepdims=True)
        fx = np.einsum("si,ij,sj->s", X, sum(c * Q[k] for k, c in fixed.items()), X)
        A = np.stack([np.einsum("si,ij,sj->s", X, Q[k], X) for k in free], axis=1)
        rows.append(-A); rhs.append(fx)
    Aub = np.concatenate(rows); bub = np.concatenate(rhs)
    obj = np.zeros(len(free)) if objective is None else np.array(objective)
    for it in range(80):
        lp = linprog(obj, A_ub=Aub, b_ub=bub, bounds=[(-1e4, 1e4)] * len(free), method="highs")
        if lp.status != 0:
            return {"status": "INFEASIBLE" if lp.status == 2 else lp.message, "iterations": it}
        c = lp.x; cuts = []; worst = 0.0
        for nm, Q in Qs.items():
            H = sum(cc * Q[k] for k, cc in fixed.items()) + sum(ci * Q[k] for ci, k in zip(c, free))
            w, V = np.linalg.eigh(0.5 * (H + H.T))
            worst = min(worst, w[0])
            if w[0] < -1e-9 * max(1.0, abs(w[-1])):
                v = V[:, 0]
                cuts.append((-np.array([v @ Q[k] @ v for k in free]), v @ sum(cc * Q[k] for k, cc in fixed.items()) @ v))
        if not cuts:
            return {"status": "FEASIBLE", "c": dict(zip(free, map(float, c))), "min_eig": float(worst), "iterations": it}
        Aub = np.concatenate([Aub, np.stack([a for a, _ in cuts])]); bub = np.concatenate([bub, np.array([b for _, b in cuts])])
    return {"status": "UNDECIDED", "min_eig": float(worst)}


def audit_F5(rng):
    g, delta = 32.0, 0.3
    M = vac(g, delta)
    p = L0.default_params(s=-1.0, g=g, delta=delta)
    ch = channel_set(g, delta, rng, extended=True)
    terms = {"I1": L0.REGISTRY["I1"].density, "I1_h": EXT.I1_h_np, "J1": EXT.J1_np, "J2": EXT.J2_np,
             "I1_frob": L0.REGISTRY["I1_frob"].density}
    Qs = {}
    for nm, (a0, Mb) in ch.items():
        Qs[nm] = {k: polarize(lambda x, f=f: density_H2(f, a0, Mb, x, p)) for k, f in terms.items()}
    fixed = {"I1": -4.0}
    out = {"channels": list(Qs)}
    # single free I1_h: the feasible interval by LP objectives
    r_min = lp_feasible(Qs, ["I1_h"], fixed, rng, objective=[1.0])     # minimize c -> hits the box bound
    r_max = lp_feasible(Qs, ["I1_h"], fixed, rng, objective=[-1.0])    # maximize c
    out["I1_h_alone"] = {"minimize_c": r_min, "maximize_c": r_max}
    log(f"F5 c_I1h alone: max feasible c = {r_max.get('c')}, min run -> {r_min.get('c')} (box -1e4)")
    # scan c = -3.9, -4, -4.1 explicitly
    scan = {}
    for c in (-3.9, -4.0, -4.1):
        scan[str(c)] = min(mineig(-4 * Q["I1"] + c * Q["I1_h"]) for Q in Qs.values())
    out["I1_h_scan_min_eig"] = scan
    for fr in (["J1"], ["J2"], ["J1", "J2"]):
        r = lp_feasible(Qs, fr, fixed, rng)
        out["free_" + "+".join(fr)] = r
        log(f"F5 free {fr}: {r['status']}")
    # controls that MUST be feasible
    ctrl = {}
    ctrl["I1_h+J1+J2"] = lp_feasible(Qs, ["I1_h", "J1", "J2"], fixed, rng)
    ctrl["I1_frob"] = lp_feasible(Qs, ["I1_frob"], fixed, rng)
    ctrl["I1_h_vs_fixed_I1=+4 (sign-flipped control)"] = lp_feasible(Qs, ["I1_h"], {"I1": +4.0}, rng)
    out["controls"] = ctrl
    for k, v in ctrl.items():
        log(f"F5 control {k}: {v['status']} {v.get('c', '')}")
    # values of the omega^2 forms on the explicit pure-time-row boost jet (S = 0, T = 2178)
    v = np.zeros(30); v[[i for i in range(10) if (IU[0][i], IU[1][i]) == (1, 1)][0]] = 1.0
    out["form_values_on_pure_time_row_jet_boost_1_A1=e1e1"] = {
        k: float(v @ Qs["boost_1"][k] @ v) for k in ("I1", "I1_h", "J1", "J2")}
    # the LP restricted to null(Q_S) of boost_1 alone (a smaller, hand-checkable certificate)
    QS, QT = QST(*ch["boost_1"])
    w, V = np.linalg.eigh(QS); Vn = V[:, w < 1e-9 * w[-1]]
    Qn = {"boost_1_nullS": {k: Vn.T @ Qs["boost_1"][k] @ Vn for k in terms}}
    sub = {}
    for fr in (["J1"], ["J2"], ["J1", "J2"], ["I1_h"]):
        r = lp_feasible(Qn, fr, fixed, rng, nsamp=500)
        sub["+".join(fr)] = r["status"]
    out["LP_on_null_QS_of_boost_1_only"] = sub
    log(f"F5 LP restricted to null(Q_S) of boost_1: {sub}")
    ok = (r_max["status"] == "FEASIBLE" and abs(r_max["c"]["I1_h"] + 4.0) < 1e-6
          and all(out["free_" + k]["status"] == "INFEASIBLE" for k in ("J1", "J2", "J1+J2"))
          and ctrl["I1_h+J1+J2"]["status"] == "FEASIBLE" and ctrl["I1_frob"]["status"] == "FEASIBLE"
          and scan["-3.9"] < 0 and scan["-4.0"] > -1e-8 and scan["-4.1"] > -1e-8)
    out["verdict"] = "CONFIRMED" if ok else "REFUTED"
    log(f"F5 {out['verdict']}: scan min eig at c=-3.9/-4/-4.1: {scan}")
    return out


def main():
    rng = np.random.default_rng(20260827)
    RES["claims"]["F1"] = audit_F1(rng)
    RES["claims"]["F2"] = audit_F2(rng)
    RES["claims"]["F3"] = audit_F3(rng)
    RES["claims"]["F4"] = audit_F4(rng)
    RES["claims"]["F5"] = audit_F5(rng)
    RES["verdicts"] = {k: v["verdict"] for k, v in RES["claims"].items()}
    RES["runtime_s"] = time.time() - T0
    RES["log"] = LOG

    def conv(o):
        if isinstance(o, (np.floating, np.integer)):
            return o.item()
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, np.bool_):
            return bool(o)
        raise TypeError(type(o))
    with open(OUT, "w") as f:
        json.dump(RES, f, indent=1, default=conv)
    log(f"verdicts: {RES['verdicts']}  runtime {RES['runtime_s']:.1f}s -> {OUT}")


if __name__ == "__main__":
    main()
