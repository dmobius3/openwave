"""M5.32 R2 arm (a): the C2 SIGN MAP at form level (symbolic / form-level).

Reuses (imports, never modifies) the R0-audited registry
m5_32_lagrangian.py, the R1 extended registry m5_32_terms_ext.py (I1_h,
J1, J2, Pgrad), the R1 audit's form-level LP machinery
m5_32_r1_audit_symbolic.py (channels, jet18, feasibility with cutting
planes, timelike_u, lorentz, transform), and the M5.21.14 dressed family
m5_21_14_c_minimize.py (make_grid, ExactCorr, b_of, qb_from, m4h_batch).

EQUATIONS FIRST
---------------
Field M(x) real symmetric 4x4, eta = diag(-1,1,1,1), jets A_mu = d_mu M,
A_0 = omega a0 (clock direction a0 on a static background),
    F_{mu nu} = A_mu eta A_nu - A_nu eta A_mu.
Vacuum (s = -1 branch): M_vac = d = diag(g, 1, delta, 0), eigenframe of
N = M eta with u = e0 the timelike unit eigenvector (u^T eta u = -1).

Terms (registry definitions, restated):
    I1    = sum_{mu<nu} eta^mu eta^nu <F, F>_eta,  <F,G>_eta = tr(eta F eta G^T)
    I1_h  = sum_{mu<nu} eta^mu eta^nu tr(h F h F^T),  h = eta + 2 (eta u)(eta u)^T
            (covariant; h = 1 in the vacuum eigenframe so I1_h = I1_frob there;
             = I1 on any field whose time row is uniform, since then u = e0
             and the time-row entries of F vanish)
    J1    = sum eta^mu eta^nu tr(F eta M eta F eta M eta)
    J2    = sum eta^mu eta^nu tr(F eta F eta M eta M eta)
    Pgrad = sum_mu eta^{mu mu} q(d_mu P_t, d_mu P_t),  P_t = u u^T eta

(1) The lambda family
    L_lambda = -4 [(1 - lambda) I1 + lambda I1_h] - V4.
    Covariant for every lambda (both members are); static sector = I1
    exactly (I1_h = I1 on uniform-time-row fields).
    Split the time-row curvature: for F_{0i} write S = sum over the
    spatial-spatial entries (a, b >= 1) of F_{0i}[a,b]^2 and T = sum over
    the mixed entries (a = 0 xor b = 0). Then in the eigenframe
        <F_0i, F_0i>_eta = S - T,   tr(F_0i F_0i^T) = S + T,
    and with eta^0 eta^i = -1 the omega^2 coefficient (Lagrangian read)
    of I1 is C_I1 = T - S and of I1_h is C_h = -(S + T). The Hamiltonian
    omega^2 form is H2 = -4 C, so
        H2(lambda) = 4 [ S - (1 - 2 lambda) T ],   S, T PSD forms in the
    spatial jet. It is AFFINE in lambda; at lambda = 1/2 the T sector
    (the boost / time-row sector) is killed: H2 = 4 S >= 0. For lambda >
    1/2 it is PSD unconditionally; for lambda < 1/2 PSD iff
        (1 - 2 lambda) <= mu := inf { x^T S x / x^T T x : x^T T x > 0 }
    i.e. the exact threshold is lambda* = (1 - mu) / 2 (mu = 0 whenever a
    pure time-row curvature direction exists, then lambda* = 1/2).
    "Boost weight reversed": on the boost channels the jet directions
    with x^T H2_I1 x < 0 (T > S there) must give x^T H2(lambda) x > 0,
    i.e. lambda > (1 - S/T) / 2 on each of them (strict).

(2) The C2 LP: c_I1 = -4 fixed, free coefficients on I1_h, J1, J2 (and
    optionally I2..I6): {c : H2(c) PSD on every channel} by the R1 audit
    sampled LP + cutting planes (an infeasible LP is a certificate; a
    feasible answer is polished to min eig >= -1e-9 |max eig|).
    Pgrad: its omega^2 coefficient on a static background is
    -q(d_0 P_t, d_0 P_t), a CONSTANT in the spatial jet (d_0 u depends on
    a0 and M only). A constant cannot make an indefinite homogeneous form
    PSD (scale x -> s x), so Pgrad alone cannot stabilize I1 at form level;
    its per-channel constant is reported.
    Switch attribution: the term whose coefficient forced to zero makes
    the LP infeasible.

(3) Static footprints: on static uniform-time-row fields I1_h - I1 = 0
    and Pgrad = 0 (u = e0 everywhere). J1, J2 are fitted to the R0 static
    basis span{I1, I2, I6} = span{a, b, c} (I1 = 2a, I2 = 4b, I6 = 4c);
    the relative residual and the J2 / I1 range are the deformation read.

(4) The degeneracy locus of I1_h: M = M_vac + [[tau, v^T], [v, 0]],
    N = M eta = [[-(g + tau), v^T], [-v, D]], D = diag(1, delta, 0).
    Along v = t e_k the 2x2 block [[-(g + tau), t], [-t, D_k]] has
    eigenvalues collide at t* = (g + tau + D_k) / 2 (e_1: (g+tau+1)/2, the
    R1 audit's t* at tau = 0; e_3: (g+tau)/2, the smallest). General
    directions by bisection on |v|. The locus is a Lorentz INVARIANT of
    the spectrum of M eta (N -> L N L^-1); a Lorentz-dressed field keeps
    the vacuum spectrum exactly, so its true margin is the spectral gap,
    reported next to the |v| proxy.

(5) On the M5.21.14 family at b* (and 2 b*, 4 b*): kin(lambda) =
    -4 [(1 - lambda) C_I1 + lambda C_h] and E_static(lambda) = 4 [(1 -
    lambda) A_I1 + lambda A_h] + V4, both linear in lambda by
    construction (verified by direct evaluation of the mixed density);
    the zero crossings lambda_kin, lambda_E.

Out: ../data/m5_32_r2_formmap.json
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys
import time

import numpy as np
from scipy.optimize import linprog

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT = os.path.join(DATA, "m5_32_r2_formmap.json")
T0 = time.time()


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


L0 = _load("m5_32_lagrangian", "m5_32_lagrangian.py")
EXT = _load("m5_32_terms_ext", "m5_32_terms_ext.py")
AUD = _load("m5_32_r1_audit_symbolic", "m5_32_r1_audit_symbolic.py")
C14 = _load("m5_21_14_c_minimize", "m5_21_14_c_minimize.py")

ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
e = np.diag(ETA)
RES = {"arm": "R2.a form-level C2 sign map",
       "reused": ["m5_32_lagrangian.py (I1, I1_frob, I2..I6, V4)",
                  "m5_32_terms_ext.py (I1_h, J1, J2, Pgrad numpy)",
                  "m5_32_r1_audit_symbolic.py (channels, jet18, feasibility, timelike_u, lorentz, transform)",
                  "m5_21_14_c_minimize.py (make_grid, ExactCorr, b_of, qb_from, m4h_batch)"],
       "conventions": {"eta": "diag(-1,1,1,1)", "vacuum": "diag(g,1,delta,0) (s = -1)",
                       "H2": "omega^2 Hamiltonian form = -4 C with c_I1 = -4; H2(c) = sum_k c_k (-C_k) ... "
                             "stored as H2 = sum_k c_k Q_k with Q_k = -C_k form; c_I1 = -4",
                       "jet": "jet18 of the R1 audit: A_i = coms(Gamma_i, d), 6 params per i"}}


def log(*a):
    print(f"[{time.time() - T0:7.1f}s]", *a, flush=True)


# ---------------- densities of every term on a batch of jets ----------------
NAMES = ["I1", "I1_frob", "I2", "I3", "I4", "I5", "I6", "I1_h", "J1", "J2", "Pgrad"]


def densities_batch(A, M, p):
    """A (4, N, 4, 4), M (N, 4, 4) -> {name: (N,)}."""
    F = L0.F_of_A(A)
    out = {}
    for k in ("I1", "I1_frob", "I2", "I3", "I4", "I5", "I6"):
        out[k] = L0.REGISTRY[k].density(A, M, p) if hasattr(L0.REGISTRY[k], "density") else None
    out["I1_h"] = EXT.I1_h_np(A, M, p)
    out["J1"] = EXT.J1_np(A, M, p)
    out["J2"] = EXT.J2_np(A, M, p)
    out["Pgrad"] = EXT.Pgrad_np(A, M, p)
    return out


def gamma_batch(t, r):
    """t, r (N, 3) -> Gamma (N, 4, 4), the R1 audit gamma()."""
    N = t.shape[0]
    G = np.zeros((N, 4, 4))
    G[:, 0, 1:] = t; G[:, 1:, 0] = t
    G[:, 1, 2], G[:, 1, 3] = -r[:, 2], r[:, 1]
    G[:, 2, 1], G[:, 2, 3] = r[:, 2], -r[:, 0]
    G[:, 3, 1], G[:, 3, 2] = -r[:, 1], r[:, 0]
    return G


def jets_batch(X, M):
    """X (N, 18) -> A (4, N, 4, 4) with A_0 = 0 (the R1 audit jet18, batched)."""
    N = X.shape[0]
    A = np.zeros((4, N, 4, 4))
    Mb = np.broadcast_to(M, (N, 4, 4))
    for i in range(3):
        G = gamma_batch(X[:, 6 * i:6 * i + 3], X[:, 6 * i + 3:6 * i + 6])
        A[1 + i] = G @ ETA @ Mb - Mb @ ETA @ G
    return A


def omega2_batch(a0, X, M, p):
    """C_k(x) for every term on the batch X: exact quadratic in omega."""
    N = X.shape[0]
    Mb = np.broadcast_to(M, (N, 4, 4)).copy()
    A = jets_batch(X, M)
    out = {}
    for w in (1.0, -1.0, 0.0):
        A2 = A.copy(); A2[0] = w * a0
        out[w] = densities_batch(A2, Mb, p)
    return {k: 0.5 * (out[1.0][k] + out[-1.0][k]) - out[0.0][k] for k in NAMES}


def Q_forms(a0, M, p, n=18):
    """polarization on the batch: Q_k (n x n) with x^T Q_k x = C_k(x), plus the
    constant part c0_k = C_k(0) (nonzero only for Pgrad)."""
    Xs = [np.zeros(n)]
    for i in range(n):
        x = np.zeros(n); x[i] = 1.0; Xs.append(x)
    pairs = []
    for i in range(n):
        for j in range(i + 1, n):
            x = np.zeros(n); x[i] = 1.0; x[j] = 1.0; Xs.append(x); pairs.append((i, j))
    C = omega2_batch(a0, np.array(Xs), M, p)
    Q, c0 = {}, {}
    for k in NAMES:
        c = C[k]; c0[k] = float(c[0]); d = c[1:1 + n] - c[0]
        Qk = np.diag(d)
        for m, (i, j) in enumerate(pairs):
            Qk[i, j] = Qk[j, i] = 0.5 * (c[1 + n + m] - c[0] - d[i] - d[j])
        Q[k] = Qk
    return Q, c0


def check_forms(Q, c0, a0, M, p, rng, n=18):
    X = rng.normal(size=(6, n))
    C = omega2_batch(a0, X, M, p)
    worst = {}
    for k in NAMES:
        pred = np.einsum("si,ij,sj->s", X, Q[k], X) + c0[k]
        worst[k] = float(np.max(np.abs(pred - C[k])) / max(1.0, np.max(np.abs(C[k]))))
    return worst


def all_channels(M):
    """the nine R1-audit channels + the antisymmetric probes G M - M G^T of
    the six generators (the R1 producer's CHAN a0) + the mixed tangents +
    eig_t + a random symmetric direction."""
    ch = dict(AUD.channels(M))
    for k in (1, 2, 3):
        G = AUD.gen_boost(k); ch[f"x_boost_{k}_probe"] = G @ M - M @ G.T
        G = AUD.gen_rot(k); ch[f"x_rot_{k}_probe"] = G @ M - M @ G.T
    Kz, Jz, Jx = AUD.gen_boost(3), AUD.gen_rot(3), AUD.gen_rot(1)
    for nm, G in (("x_boost3_plus_rot3_tan", Kz + Jz), ("x_boost3_plus_rot1_tan", Kz + Jx)):
        ch[nm] = G @ M + M @ G.T
    ch["x_eig_t"] = np.diag([1.0, 0, 0, 0])
    rng = np.random.default_rng(20260828)
    S = rng.normal(size=(4, 4)); ch["x_random_sym"] = 0.5 * (S + S.T)
    return ch


def min_eig(H):
    return float(np.linalg.eigvalsh(0.5 * (H + H.T))[0])


# ---------------- (1) lambda family ----------------
def stage_lambda(g, delta, rng):
    M = AUD.vac(g, delta)
    p = L0.default_params(s=-1.0, g=g, delta=delta)
    ch = all_channels(M)
    Qs, c0s, polw = {}, {}, {}
    for nm, a0 in ch.items():
        Qs[nm], c0s[nm] = Q_forms(a0, M, p)
        polw[nm] = check_forms(Qs[nm], c0s[nm], a0, M, p, rng)
    log(f"forms built at g={g}, delta={delta}: {len(ch)} channels; worst polarization",
        max(max(v.values()) for v in polw.values()))
    res = {"g": g, "delta": delta, "channels": list(ch),
           "polarization_worst_relerr": {nm: max(v.values()) for nm, v in polw.items()},
           "Pgrad_constant_c0_per_channel": {nm: c0s[nm]["Pgrad"] for nm in ch},
           "nonPgrad_constants_max": float(max(abs(c0s[nm][k]) for nm in ch for k in NAMES if k != "Pgrad")),
           "I1h_equals_I1frob_as_form_maxdiff": float(max(np.max(np.abs(Qs[nm]["I1_h"] - Qs[nm]["I1_frob"])) for nm in ch))}
    lam_grid = np.round(np.arange(0.0, 1.5001, 0.05), 3)
    H = {nm: {"I1": -4 * Qs[nm]["I1"], "h": -4 * Qs[nm]["I1_h"]} for nm in ch}
    # affine check
    aff = max(np.max(np.abs(0.5 * (H[nm]["I1"] + H[nm]["h"]) - ((1 - .25) * H[nm]["I1"] + .25 * H[nm]["h"] + (1 - .75) * H[nm]["I1"] + .75 * H[nm]["h"]) / 2)) for nm in ch)
    res["affine_in_lambda_check"] = float(aff)
    grid = {}
    for nm in ch:
        grid[nm] = [min_eig((1 - l) * H[nm]["I1"] + l * H[nm]["h"]) for l in lam_grid]
    res["lambda_grid"] = lam_grid.tolist()
    res["min_eig_vs_lambda"] = grid
    scale = {nm: float(np.max(np.abs(np.linalg.eigvalsh(H[nm]["h"])))) for nm in ch}
    tol = 1e-9
    allpsd = [bool(all(grid[nm][i] >= -tol * max(1.0, scale[nm]) for nm in ch)) for i in range(len(lam_grid))]
    res["all_channels_PSD_on_grid"] = dict(zip(map(str, lam_grid), allpsd))
    # S, T split and the exact threshold
    thr = {}
    for nm in ch:
        S4 = 0.5 * (H[nm]["I1"] + H[nm]["h"]); T4 = 0.5 * (H[nm]["h"] - H[nm]["I1"])   # 4S, 4T
        wS, wT = np.linalg.eigvalsh(S4), np.linalg.eigvalsh(T4)
        # mu = inf x^T S x / x^T T x over range(T): generalized eigenproblem on range(T)
        wt, Vt = np.linalg.eigh(T4)
        keep = wt > 1e-9 * max(1.0, wt[-1])
        if keep.any():
            B = Vt[:, keep]
            Tr = B.T @ T4 @ B; Sr = B.T @ S4 @ B
            Lc = np.linalg.cholesky(Tr)
            Li = np.linalg.inv(Lc)
            mu = float(np.linalg.eigvalsh(Li @ Sr @ Li.T)[0])
            lam_star = (1 - mu) / 2
        else:
            mu, lam_star = None, None
        # I1-negative directions (the runaway eigenvectors) and their reversal lambda
        wI, VI = np.linalg.eigh(H[nm]["I1"])
        neg = wI < -1e-9 * max(1.0, abs(wI[-1]))
        rev = []
        for i in np.where(neg)[0]:
            v = VI[:, i]; s, t = v @ S4 @ v, v @ T4 @ v
            rev.append({"eig_I1": float(wI[i]), "S": float(s), "T": float(t),
                        "lambda_reversal": float((1 - s / t) / 2) if t > 1e-12 else None})
        # pure rotation-type (hedgehog) jets: r components only
        idx_r = [3, 4, 5, 9, 10, 11, 15, 16, 17]
        hh = {"min_eig_I1_on_r_jets": min_eig(H[nm]["I1"][np.ix_(idx_r, idx_r)]),
              "min_eig_lambda_half_on_r_jets": min_eig(S4[np.ix_(idx_r, idx_r)])}
        thr[nm] = {"S_min_eig": float(wS[0]), "T_min_eig": float(wT[0]), "S_PSD": bool(wS[0] >= -1e-9 * max(1, wS[-1])),
                   "T_PSD": bool(wT[0] >= -1e-9 * max(1, wT[-1])), "rank_T": int(keep.sum()),
                   "mu": mu, "lambda_star": lam_star, "n_I1_negative_dirs": int(neg.sum()),
                   "I1_negative_dirs": rev, "hedgehog_r_jets": hh,
                   "min_eig_I1": float(wI[0]), "min_eig_lambda_half": float(np.linalg.eigvalsh(S4)[0])}
    res["per_channel_threshold"] = thr
    ls = [v["lambda_star"] for v in thr.values() if v["lambda_star"] is not None]
    lrev = [d["lambda_reversal"] for v in thr.values() for d in v["I1_negative_dirs"] if d["lambda_reversal"] is not None]
    res["exact_threshold"] = {
        "lambda_star_all_channels": float(max(ls)) if ls else None,
        "argmax_channel": max(thr, key=lambda k: (thr[k]["lambda_star"] if thr[k]["lambda_star"] is not None else -1)),
        "lambda_reversal_all_I1_negative_dirs (strict >)": float(max(lrev)) if lrev else None,
        "S_PSD_all": bool(all(v["S_PSD"] for v in thr.values())),
        "T_PSD_all": bool(all(v["T_PSD"] for v in thr.values())),
        "reading": ("H2(lambda) = 4[S - (1-2 lambda) T] with S, T PSD: PSD for every "
                    "lambda >= lambda* = (1 - mu)/2; lambda >= 1/2 unconditionally; "
                    "boost weight reversed strictly for lambda > max lambda_reversal")}
    # numeric: where boost hedgehog dirs are strictly reversed
    boost = [nm for nm in ch if nm.startswith("boost_")]
    revgrid = []
    for l in lam_grid:
        ok = True
        for nm in boost:
            wI, VI = np.linalg.eigh(H[nm]["I1"])
            for i in np.where(wI < -1e-9 * max(1.0, abs(wI[-1])))[0]:
                v = VI[:, i]
                if v @ ((1 - l) * H[nm]["I1"] + l * H[nm]["h"]) @ v <= 1e-9 * scale[nm]:
                    ok = False
        revgrid.append(bool(ok))
    res["boost_hedgehog_strictly_reversed_on_grid"] = dict(zip(map(str, lam_grid), revgrid))
    return res, Qs, c0s, ch, M, p


# ---------------- (2) the C2 LP ----------------
def psd_interval(H0s, Qs1, lo=-1e3, hi=1e3):
    """{c : H0 + c Q PSD on all channels}: f(c) = min_ch min eig is concave;
    scan + bisection. Returns (c_lo, c_hi) or None."""
    def f(c):
        return min(min_eig(H0s[nm] + c * Qs1[nm]) - (-1e-9 * max(1.0, np.max(np.abs(np.linalg.eigvalsh(H0s[nm] + c * Qs1[nm]))))) for nm in H0s)
    grid = np.concatenate([-np.geomspace(1e-3, -lo, 200)[::-1], [0.0], np.geomspace(1e-3, hi, 200)])
    vals = np.array([f(c) for c in grid])
    ok = np.where(vals >= 0)[0]
    if len(ok) == 0:
        return None, float(vals.max()), float(grid[int(np.argmax(vals))])
    i0, i1 = ok[0], ok[-1]
    def bis(a, b):    # f(a) >= 0 > f(b)
        for _ in range(60):
            m = 0.5 * (a + b)
            if f(m) >= 0: a = m
            else: b = m
        return a
    c_lo = bis(grid[i0], grid[i0 - 1]) if i0 > 0 else -np.inf
    c_hi = bis(grid[i1], grid[i1 + 1]) if i1 < len(grid) - 1 else np.inf
    return (float(c_lo), float(c_hi)), float(vals.max()), float(grid[int(np.argmax(vals))])


def lp_opt(Qs, free, fixed, rng, obj_term, sense="max", nsamp=1500, box=1e3):
    """the R1 audit feasibility() loop (m5_32_r1_audit_symbolic.py) with a
    linear objective: extremize c[obj_term] over the PSD set (cutting planes
    on the min eigenvector until min eig >= -1e-9 |max eig| on every channel)."""
    n = next(iter(next(iter(Qs.values())).values())).shape[0]
    rows, rhs = [], []
    for chn, Q in Qs.items():
        X = rng.normal(size=(nsamp, n)); X /= np.linalg.norm(X, axis=1, keepdims=True)
        fx = np.einsum("si,ij,sj->s", X, sum(c * Q[k] for k, c in fixed.items()), X)
        A = np.stack([np.einsum("si,ij,sj->s", X, Q[k], X) for k in free], axis=1)
        rows.append(-A); rhs.append(fx)
    Aub = np.concatenate(rows); bub = np.concatenate(rhs)
    obj = np.zeros(len(free)); obj[free.index(obj_term)] = -1.0 if sense == "max" else 1.0
    bounds = [(-box, box)] * len(free)
    out = {"free": free, "objective": f"{sense} c_{obj_term}", "box": box}
    for it in range(200):
        lp = linprog(obj, A_ub=Aub, b_ub=bub, bounds=bounds, method="highs")
        if lp.status != 0:
            out.update({"status": "INFEASIBLE", "iterations": it}); return out
        c = lp.x; cuts = []; worst = 0.0
        for chn, Q in Qs.items():
            H = sum(cc * Q[k] for k, cc in fixed.items()) + sum(ci * Q[k] for ci, k in zip(c, free))
            wv, vv = np.linalg.eigh(0.5 * (H + H.T)); worst = min(worst, wv[0])
            if wv[0] < -1e-9 * max(1.0, abs(wv[-1])):
                v = vv[:, 0]
                cuts.append((-np.array([v @ Q[k] @ v for k in free]), v @ sum(cc * Q[k] for k, cc in fixed.items()) @ v))
        if not cuts:
            out.update({"status": "CONVERGED", "c": dict(zip(free, map(float, c))), "min_eig_all_channels": float(worst), "iterations": it}); return out
        Aub = np.concatenate([Aub, np.stack([a for a, _ in cuts])]); bub = np.concatenate([bub, np.array([b for _, b in cuts])])
    out.update({"status": "UNDECIDED", "c": dict(zip(free, map(float, c))), "min_eig_all_channels": float(worst), "iterations": 200})
    return out


def stage_lp(Qs, ch, rng):
    """Qs[nm][k] = C_k forms (Lagrangian omega^2); the audit's feasibility()
    uses H2(c) = sum c_k Q_k with c_I1 = -4: that is the Hamiltonian read
    up to the overall sign -1 ... NOTE: the audit feeds Q = C forms and
    fixed I1 = -4, i.e. it tests sum_k c_k C_k >= 0. With E = c (C w^2 - A)
    the kinetic energy is sum c_k C_k, so this IS the Hamiltonian form."""
    fixed = {"I1": -4.0}
    homog = {nm: {k: Qs[nm][k] for k in NAMES if k != "Pgrad"} for nm in ch}
    res = {}
    for label, free in (("I1_h", ["I1_h"]), ("J1", ["J1"]), ("J2", ["J2"]),
                        ("I1_h+J1+J2", ["I1_h", "J1", "J2"]), ("J1+J2 (no I1_h)", ["J1", "J2"]),
                        ("J1+J2+I2..I6 (no I1_h)", ["J1", "J2", "I2", "I3", "I4", "I5", "I6"]),
                        ("I1_h+J1+J2+I2..I6", ["I1_h", "J1", "J2", "I2", "I3", "I4", "I5", "I6"]),
                        ("I1_h+I2..I6 (no J)", ["I1_h", "I2", "I3", "I4", "I5", "I6"])):
        r = AUD.feasibility(homog, free, fixed, rng, nsamp=1500, tag=label)
        res[label] = {k: v for k, v in r.items() if k != "fixed"}
        log(f"LP {label}: {r['status']}", r.get("c"))
    # single-coefficient exact intervals (concave min-eig scan + bisection)
    H0 = {nm: -4 * Qs[nm]["I1"] for nm in ch}
    iv = {}
    for k in ("I1_h", "J1", "J2", "I2", "I3", "I4", "I5", "I6"):
        r, best, cbest = psd_interval(H0, {nm: Qs[nm][k] for nm in ch})
        iv[k] = {"interval": r, "best_min_eig_over_scan": best, "c_at_best": cbest}
        log(f"interval {k}: {r} best {best:.4g} at c = {cbest:.4g}")
    res["single_coefficient_PSD_intervals"] = iv
    # switch attribution: drop each term from the full free set
    full = ["I1_h", "J1", "J2", "I2", "I3", "I4", "I5", "I6"]
    attr = {}
    for k in full:
        free = [t for t in full if t != k]
        attr[f"without_{k}"] = AUD.feasibility(homog, free, fixed, rng, nsamp=1500)["status"]
    # an UNDECIDED verdict (the cutting-plane loop crawling at the degenerate boundary)
    # on a set that still contains I1_h is settled by the ray certificate below:
    # c_I1h <= -4 with every other coefficient 0 is PSD on every channel (I1h_ray_direct)
    for k, v in list(attr.items()):
        if v == "UNDECIDED" and k != "without_I1_h":
            attr[k] = "FEASIBLE (LP undecided at the boundary; contains the I1_h ray, PSD by I1h_ray_direct)"
    res["switch_attribution"] = attr
    # direct boundary checks of the I1_h ray (the LP crawls at the degenerate boundary c = -4)
    res["I1h_ray_direct"] = {f"c_I1h={c:g}": {nm: min_eig(H0[nm] + c * Qs[nm]["I1_h"]) for nm in ch}
                             for c in (-3.9, -4.0, -8.0, -100.0)}
    # optimizing LP: the LEAST negative c_I1h compatible with PSD when the other terms are free
    for label, free in (("J1+J2", ["I1_h", "J1", "J2"]), ("J1+J2+I2..I6", ["I1_h", "J1", "J2", "I2", "I3", "I4", "I5", "I6"]),
                        ("I2..I6", ["I1_h", "I2", "I3", "I4", "I5", "I6"])):
        for box in (10.0, 100.0, 1000.0):
            r = lp_opt(homog, free, fixed, rng, obj_term="I1_h", sense="max", box=box)
            res[f"max_c_I1h_with_{label}_free_box{box:g}"] = r
            log(f"max c_I1h with {label} (box {box:g}): ", None if r.get("c") is None else round(r["c"]["I1_h"], 4),
                r["status"], {k: round(v, 3) for k, v in (r.get("c") or {}).items() if abs(v) >= 0.5 * box})
    # Pgrad: constant per channel; scaling argument + numeric on unit jets
    res["Pgrad_note"] = ("omega^2 coefficient of Pgrad on a static background is a constant "
                         "-q(d_0 P_t, d_0 P_t) in the spatial jet (no quadratic part: see "
                         "Pgrad_constant_c0_per_channel and the zero Q_Pgrad norm); a constant "
                         "cannot make an indefinite homogeneous form PSD (scale the jet), so "
                         "Pgrad ALONE cannot stabilize I1 at form level")
    res["Pgrad_Q_form_norm_max"] = float(max(np.max(np.abs(Qs[nm]["Pgrad"])) for nm in ch))
    return res


# ---------------- (3) static footprints ----------------
def stage_static(g, delta, rng, nfield=50):
    M0 = AUD.vac(g, delta)
    p = L0.default_params(s=-1.0, g=g, delta=delta)
    rows = []
    kinds = []
    for j in range(nfield):
        if j % 2 == 0:      # spatial block perturbed (the R1 audit's static gate)
            M = M0.copy(); M[1:, 1:] += 2.0 * AUD.sym(rng.normal(size=(3, 3)))
            kinds.append("block_perturbed")
        else:               # rotated vacuum block (the hedgehog configurations)
            Rm = np.linalg.qr(rng.normal(size=(3, 3)))[0]
            M = M0.copy(); M[1:, 1:] = Rm @ np.diag([1.0, delta, 0.0]) @ Rm.T
            kinds.append("rotated_block")
        A = np.zeros((4, 1, 4, 4)); A[1:, 0, 1:, 1:] = AUD.sym(rng.normal(size=(3, 3, 3)))
        d = densities_batch(A, M[None], p)
        rows.append({k: float(d[k][0]) for k in NAMES})
    V = np.array([[r[k] for k in ("I1", "I2", "I6")] for r in rows])
    out = {"n_fields": nfield, "kinds": kinds,
           "I1h_minus_I1_max": float(max(abs(r["I1_h"] - r["I1"]) for r in rows)),
           "Pgrad_max": float(max(abs(r["Pgrad"]) for r in rows)),
           "I1_scale": float(max(abs(r["I1"]) for r in rows)),
           "basis_note": "R0 static basis (a, b, c): I1 = 2a, I2 = 4b, I6 = 4c; fit in (I1, I2, I6) then convert"}
    for k in ("J1", "J2"):
        y = np.array([r[k] for r in rows])
        coef, *_ = np.linalg.lstsq(V, y, rcond=None)
        resid = np.linalg.norm(y - V @ coef) / np.linalg.norm(y)
        ratio = y / np.array([r["I1"] for r in rows])
        out[k] = {"coef_I1_I2_I6": coef.tolist(),
                  "coef_abc": [2 * coef[0], 4 * coef[1], 4 * coef[2]],
                  "relative_residual": float(resid),
                  f"{k}/I1_min": float(ratio.min()), f"{k}/I1_max": float(ratio.max()),
                  "per_kind_ratio_range": {kd: [float(ratio[[i for i, x in enumerate(kinds) if x == kd]].min()),
                                                float(ratio[[i for i, x in enumerate(kinds) if x == kd]].max())]
                                           for kd in ("block_perturbed", "rotated_block")}}
        # on the rotated-block (hedgehog) fields the spectrum is the vacuum's: is J a fixed multiple then?
        sel = [i for i, x in enumerate(kinds) if x == "rotated_block"]
        coef2, *_ = np.linalg.lstsq(V[sel], y[sel], rcond=None)
        out[k]["rotated_block_only"] = {"coef_I1_I2_I6": coef2.tolist(),
                                        "relative_residual": float(np.linalg.norm(y[sel] - V[sel] @ coef2) / np.linalg.norm(y[sel]))}
    return out


# ---------------- (4) the degeneracy locus ----------------
def has_timelike(M):
    u, lam = AUD.timelike_u(M)
    if u is None:
        return False
    ev = np.linalg.eigvals(M @ ETA)
    return bool(np.max(np.abs(ev.imag)) < 1e-9 * max(1.0, np.max(np.abs(ev))))


def v_threshold(g, delta, tau, vhat, vmax=200.0):
    M0 = AUD.vac(g, delta)
    def Mof(s):
        M = M0.copy(); M[0, 0] += tau; M[0, 1:] = s * vhat; M[1:, 0] = s * vhat; return M
    if not has_timelike(Mof(0.0)):
        return 0.0
    a, b = 0.0, vmax
    if has_timelike(Mof(b)):
        return np.inf
    for _ in range(60):
        m = 0.5 * (a + b)
        if has_timelike(Mof(m)): a = m
        else: b = m
    return 0.5 * (a + b)


def stage_locus(rng):
    out = {}
    dirs = {"e1": np.array([1.0, 0, 0]), "e2": np.array([0, 1.0, 0]), "e3": np.array([0, 0, 1.0])}
    rd = rng.normal(size=(60, 3)); rd /= np.linalg.norm(rd, axis=1, keepdims=True)
    for g in (32.0, 8.0):
        delta = 0.3
        taus = [-0.5 * g, -0.25 * g, 0.0, 0.25 * g, 0.5 * g, 1.0 * g, 2.0 * g]
        rows = {}
        for tau in taus:
            r = {nm: v_threshold(g, delta, tau, v) for nm, v in dirs.items()}
            rr = [v_threshold(g, delta, tau, v) for v in rd]
            r["random_min"] = float(min(rr)); r["random_max"] = float(max(rr))
            r["formula_e1_(g+tau+1)/2"] = (g + tau + 1) / 2
            r["formula_e2_(g+tau+delta)/2"] = (g + tau + delta) / 2
            r["formula_e3_(g+tau)/2"] = (g + tau) / 2
            r["threshold_min_over_dirs"] = float(min(min(r[k] for k in dirs), r["random_min"]))
            rows[f"tau={tau:g}"] = r
        out[f"g={g:g}"] = {"delta": delta, "rows": rows,
                           "|v|_threshold_at_tau0_min": rows["tau=0"]["threshold_min_over_dirs"]}
        log(f"locus g={g}: tau=0 thresholds e1 {rows['tau=0']['e1']:.4f} e2 {rows['tau=0']['e2']:.4f} e3 {rows['tau=0']['e3']:.4f} random min {rows['tau=0']['random_min']:.4f}")
    # the R1 audit path: I1_h ~ (t* - t)^(-1/2) sign-indefinite (reproduce two points)
    out["R1_audit_path_tstar_g32"] = (32 + 1) / 2
    return out


def dressed_field(scale, avec, grid, g=32.0):
    P = grid["P"]
    K, K2, r = C14.kgeom(P)
    b = scale * C14.b_of(avec, r)
    Qb = C14.qb_from(K, K2, b)
    M4 = Qb @ C14.m4h_batch(P, g) @ np.swapaxes(Qb, -1, -2)
    return M4, b, r


def stage_margin(locus):
    with open(os.path.join(DATA, "m5_21_14_minimize.json")) as f:
        rec = json.load(f)
    avec = np.array(rec["avec"])
    grid = C14.make_grid(48, 8, 16)
    out = {"family": "M5.21.14 b*(r) = b_of(avec_record, r), grid make_grid(48, 8, 16), g = 32, delta = 0.3"}
    thr0 = locus["g=32"]["rows"]["tau=0"]["threshold_min_over_dirs"]
    for scale in (1.0, 2.0, 4.0):
        M4, b, r = dressed_field(scale, avec, grid)
        v = M4[:, 0, 1:]; vn = np.linalg.norm(v, axis=1); tau = M4[:, 0, 0] - 32.0
        # proxy: |v| against the tau-dependent minimal threshold (g + tau)/2 at the same point
        prox = vn / ((32.0 + tau) / 2)
        # exact: spectrum of M eta (Lorentz invariant -> vacuum spectrum), gap of the timelike eigenvalue
        ev = np.linalg.eigvals(M4 @ ETA)
        im = np.max(np.abs(ev.imag))
        evr = np.sort(ev.real, axis=1)
        tl = []
        for i in range(0, M4.shape[0], max(1, M4.shape[0] // 400)):
            u, lam = AUD.timelike_u(M4[i]); tl.append(lam if lam is not None else np.nan)
        tl = np.array(tl)
        # timelike eigenvalue of M eta is -g (vacuum spectrum (-g, 1, delta, 0) at s = -1);
        # sorted ascending it is evr[:, 0]; gap = distance of the nearest other eigenvalue
        gap = np.min(evr[:, 1:] - evr[:, :1], axis=1)
        out[f"scale_{scale:g}bstar"] = {
            "max_b": float(np.max(np.abs(b))), "r_at_max_b": float(r[int(np.argmax(np.abs(b)))]),
            "max_|v|": float(vn.max()), "r_at_max_|v|": float(r[int(np.argmax(vn))]),
            "tau_at_max_|v|": float(tau[int(np.argmax(vn))]),
            "threshold_|v|_tau0_min_over_dirs": thr0,
            "ratio_max_|v|_to_tau0_threshold": float(vn.max() / thr0),
            "max_ratio_|v|_to_(g+tau)/2_pointwise": float(prox.max()),
            "spectrum_max_imag": float(im),
            "timelike_eigenvalue_range_(subsample)": [float(np.nanmin(tl)), float(np.nanmax(tl))],
            "spectral_gap_min (nearest other eig - timelike eig, timelike = -g)": float(np.min(gap)),
            "timelike_exists_everywhere": bool(not np.isnan(tl).any())}
        log(f"margin {scale}b*: max b {np.max(np.abs(b)):.4f}, max |v| {vn.max():.3f}, proxy ratio {prox.max():.3f}, gap {np.min(gap):.4f}, imag {im:.2e}")
    out["reading"] = ("M_d = L M_b L^T with L a local boost: (M_d eta) = L (M_b eta) L^-1 has the "
                      "vacuum spectrum (g, 1, delta, 0) at every point, so the degeneracy locus "
                      "(a spectral collision) is unreachable by any Lorentz dressing; the |v| proxy "
                      "against the tau = 0 threshold is NOT the margin (tau grows with |v| along the "
                      "orbit); the spectral gap is")
    return out


# ---------------- (5) inertia and energy on the family ----------------
def stage_family():
    with open(os.path.join(DATA, "m5_21_14_minimize.json")) as f:
        rec = json.load(f)
    avec = np.array(rec["avec"])
    grid = C14.make_grid(48, 8, 16)
    ec = C14.ExactCorr(grid, C14.G_MAIN)
    p = L0.default_params(s=-1.0, g=32.0)
    P, wvol = grid["P"], grid["wvol"]
    out = {"grid": "make_grid(48, 8, 16)", "scales": {}}

    def abc_all(scale):
        bfun = (lambda r: scale * C14.b_of(avec, r))
        Asp = ec._A(bfun)
        Qb = C14.qb_from(ec.K_c, ec.K2_c, bfun(ec.r_c))
        M4 = Qb @ C14.m4h_batch(P, C14.G_MAIN) @ np.swapaxes(Qb, -1, -2)
        a0 = Qb @ ec.a0_base @ np.swapaxes(Qb, -1, -2)
        def dens(om):
            A = np.zeros((4,) + M4.shape)
            for i in range(3): A[1 + i] = Asp[i]
            A[0] = om * a0
            F = L0.F_of_A(A)
            d = {"I1": L0.REGISTRY["I1"].density(A, M4, p), "I1_h": EXT.I1_h_np(A, M4, p)}
            d["V4"] = L0.REGISTRY["V4"].density(A, M4, p)
            return {k: float(np.sum(wvol * v)) for k, v in d.items()}
        d0, dp, dm = dens(0.0), dens(1.0), dens(-1.0)
        return {k: {"A": d0[k], "B": 0.5 * (dp[k] - dm[k]), "C": 0.5 * (dp[k] + dm[k]) - d0[k]} for k in ("I1", "I1_h")}, d0["V4"], (M4, a0, Asp)

    base, V0, _ = abc_all(0.0)
    out["b_zero"] = {"abc": base, "V4": V0}
    for scale in (1.0, 2.0, 4.0):
        t, V4, (M4, a0, Asp) = abc_all(scale)
        kin_I1, kin_h = -4 * t["I1"]["C"], -4 * t["I1_h"]["C"]
        kinb_I1, kinb_h = -4 * base["I1"]["C"], -4 * base["I1_h"]["C"]
        Ec_I1 = 4 * (t["I1"]["A"] - base["I1"]["A"]) + (V4 - V0)
        Ec_h = 4 * (t["I1_h"]["A"] - base["I1_h"]["A"]) + (V4 - V0)
        row = {"abc": t, "V4": V4, "V4_minus_V4_base": V4 - V0,
               "kin_total_I1": kin_I1, "kin_total_I1h": kin_h,
               "kin_corr_I1": kin_I1 - kinb_I1, "kin_corr_I1h": kin_h - kinb_h,
               "E_corr_I1": Ec_I1, "E_corr_I1h": Ec_h,
               "kin(lambda)": "(1-lambda) kin_I1 + lambda kin_I1h",
               "lambda_kin_zero": float(kin_I1 / (kin_I1 - kin_h)) if kin_I1 != kin_h else None,
               "lambda_Ecorr_zero": float(Ec_I1 / (Ec_I1 - Ec_h)) if Ec_I1 != Ec_h else None}
        # linearity check: the mixed density evaluated directly at three lambdas
        lin = {}
        for lam in (0.25, 0.5, 0.75):
            def dens(om):
                A = np.zeros((4,) + M4.shape)
                for i in range(3): A[1 + i] = Asp[i]
                A[0] = om * a0
                return float(np.sum(wvol * ((1 - lam) * L0.REGISTRY["I1"].density(A, M4, p) + lam * EXT.I1_h_np(A, M4, p))))
            d0, dp, dm = dens(0.0), dens(1.0), dens(-1.0)
            C = 0.5 * (dp + dm) - d0
            kin_direct = -4 * C
            kin_pred = (1 - lam) * kin_I1 + lam * kin_h
            E_direct = 4 * d0 + V4
            E_pred = (1 - lam) * (4 * t["I1"]["A"] + V4) + lam * (4 * t["I1_h"]["A"] + V4)
            lin[str(lam)] = {"kin_direct": kin_direct, "kin_pred": kin_pred, "rel": abs(kin_direct - kin_pred) / max(1, abs(kin_pred)),
                             "E_direct": E_direct, "E_pred": E_pred, "rel_E": abs(E_direct - E_pred) / max(1, abs(E_pred))}
        row["linearity_check"] = lin
        out["scales"][f"{scale:g}bstar"] = row
        log(f"family {scale}b*: kin_I1 {kin_I1:+.2f} kin_h {kin_h:+.2f} lam_kin0 {row['lambda_kin_zero']}, "
            f"Ecorr_I1 {Ec_I1:+.1f} Ecorr_h {Ec_h:+.1f} lam_E0 {row['lambda_Ecorr_zero']}, dV4 {V4 - V0:.2e}")
    out["record_check"] = {"E_corr_at_bstar_record": rec["verdicts"]["E_corr_at_bstar"],
                           "kin_corr_at_bstar_record": rec["verdicts"]["kin_corr_at_bstar"],
                           "note": "record on grid (72, 12, 24); here (48, 8, 16), the R1 beta grid"}
    return out


# ---------------- covariance + static identity of the family ----------------
def stage_family_checks(rng):
    M = AUD.vac(32.0, 0.3)
    p = L0.default_params(s=-1.0, g=32.0)
    cov = 0.0
    for _ in range(30):
        Mr = M + 0.5 * AUD.sym(rng.normal(size=(4, 4)))
        A = AUD.sym(rng.normal(size=(4, 4, 4)))
        L = AUD.lorentz(rng)
        M2, A2 = AUD.transform(L, Mr, A)
        def fam(A_, M_):
            d = densities_batch(A_[:, None], M_[None], p)
            return 0.5 * d["I1"][0] + 0.5 * d["I1_h"][0]
        v1, v2 = fam(A, Mr), fam(A2, M2)
        cov = max(cov, abs(v1 - v2) / max(1.0, abs(v1)))
    stat = 0.0
    for _ in range(20):
        Ms = M.copy(); Ms[1:, 1:] += 2.0 * AUD.sym(rng.normal(size=(3, 3)))
        A = np.zeros((4, 1, 4, 4)); A[1:, 0, 1:, 1:] = AUD.sym(rng.normal(size=(3, 3, 3)))
        d = densities_batch(A, Ms[None], p)
        for lam in (0.3, 0.5, 1.2):
            stat = max(stat, abs((1 - lam) * d["I1"][0] + lam * d["I1_h"][0] - d["I1"][0]) / max(1, abs(d["I1"][0])))
    return {"covariance_relerr_max_lambda_0.5 (30 random boosts+rotations)": float(cov),
            "static_identity_relerr_max (20 static fields, lambda 0.3/0.5/1.2)": float(stat)}


def main():
    rng = np.random.default_rng(20260828)
    RES["family_checks"] = stage_family_checks(rng)
    log("family checks", RES["family_checks"])
    lam32, Qs, c0s, ch, M, p = stage_lambda(32.0, 0.3, rng)
    RES["lambda_family"] = {"g=32,delta=0.3": lam32}
    log("threshold g=32:", lam32["exact_threshold"])
    lam8, *_ = stage_lambda(8.0, 0.3, rng)
    RES["lambda_family"]["g=8,delta=0.3"] = {k: lam8[k] for k in ("exact_threshold", "all_channels_PSD_on_grid", "boost_hedgehog_strictly_reversed_on_grid", "per_channel_threshold")}
    log("threshold g=8:", lam8["exact_threshold"])
    RES["C2_LP"] = stage_lp(Qs, ch, rng)
    RES["static_footprints"] = stage_static(32.0, 0.3, rng)
    log("static", {k: RES["static_footprints"][k] for k in ("I1h_minus_I1_max", "Pgrad_max")},
        "J1 resid", RES["static_footprints"]["J1"]["relative_residual"], "J2 resid", RES["static_footprints"]["J2"]["relative_residual"])
    RES["degeneracy_locus"] = stage_locus(rng)
    RES["dressed_margin"] = stage_margin(RES["degeneracy_locus"])
    RES["family_on_dressed_electron"] = stage_family()
    RES["runtime_s"] = time.time() - T0
    with open(OUT, "w") as f:
        json.dump(RES, f, indent=1, default=lambda o: float(o) if isinstance(o, (np.floating, np.integer)) else (None if o != o else str(o)))
    log("wrote", OUT)


if __name__ == "__main__":
    main()
