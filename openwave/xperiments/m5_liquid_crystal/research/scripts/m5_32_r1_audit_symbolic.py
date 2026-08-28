"""M5.32 R1 arm (a) INDEPENDENT ADVERSARIAL AUDIT of the symbolic claims.

Built without reading the producer's m5_32_r1_a_symbolic.py,
m5_32_terms_ext.py or m5_32_r1_symbolic.json. Only the R0-audited
registry m5_32_lagrangian.py is imported (I1, I1_frob, I2..I6, V4 and
their conventions); every contraction below is re-implemented here with
explicit einsum metrics and cross-checked against the registry once.

EQUATIONS FIRST
---------------
Field M(x) real symmetric 4x4, eta = diag(-1,1,1,1). Jets A_mu = d_mu M,
A_0 = omega a0 (a clock direction a0 at a static background).
    F_{mu nu} = A_mu eta A_nu - A_nu eta A_mu,  F[mu,nu,a,b] (slots 0-1
    derivative, 2-3 internal), M -> L M L^T, d_mu covariant.
Contraction rule (registry docstring): DD eta, II eta, DI delta.

Registered even invariants (my einsum re-implementation, checked):
    I1 = 1/2 eta^m eta^n eta_a eta_b F[mnab] F[mnab]
    I1_frob = 1/2 eta^m eta^n F[mnab] F[mnab]
    I2 = F[mnab] F[abmn]
    I3 = eta^m eta_b F[mnab] F[manb]
    R[n,a] = sum_m F[m,n,a,m];  I4 = eta_n eta_a R[na]^2;
    I5 = R[na] R[an];  I6 = (sum_n R[nn])^2

A1. Parity-odd (one-epsilon) quadratic invariants: choose 4 of the 8
slots of F (x) F for eps^{....} (an epsilon slot on a derivative slot
pairs by delta, on an internal slot by eta, i.e. lower the internal index
first), pair the remaining 4 slots by the rule: 70 x 3 = 210 contractions.
Producer names E1 = (eps^{mnab} F_{mnab}) R, E2 = eps^{mn a r} F_{mnab}
R_r^b, E3 = eps^{mn a c} F_{mna}^r R_{rc}. Rank on 300 random jets.

A2. N1 = 4 I3 - I2 - 2 I1, N2 = 4 I4 - I6 - 2 I1, N3 = 4 I5 - I2 - I6
vanish on every STATIC UNIFORM-TIME-ROW field (my reading: A_0 = 0 and
A_i[0,:] = A_i[:,0] = 0, i.e. only the spatial 3x3 block of M varies).

A3. Channels: a0 = G M_vac + M_vac G^T for the six Lorentz generators
(boost_k: G[0k] = G[k0] = 1; rot_k: G[ij] = -eps_{kij}), plus the
notebook local clock a0 = coms(Gamma_0, d) (Gamma_0 with t = (1, .7, -.4),
r = (.5, -.8, .3), also the pure-t and pure-r variants). Spatial jet
family: A_i = coms(Gamma_i, d), 6 parameters per i (18 in all), Gamma_i
the notebook generator with time row t_i (symmetric) and spatial block
antisymmetric r_i (registry gamma_mu, boost_style real). Each term is
exactly quadratic in omega: I = A + B omega + C omega^2; C is a quadratic
form in the 18 jet entries, its 18x18 matrix Q_k obtained by exact
polarization. H2(c) = sum_k c_k Q_k with c_I1 = -4 fixed.
Feasibility: sampled LP (2000 jets per channel, all channels jointly):
{c : x^T H2(c) x >= 0 for every sample x}; an infeasible LP is a rigorous
certificate that no c exists (sampled constraints are necessary). If
feasible, a cutting-plane loop adds the minimal eigenvector of H2(c) as a
new constraint until min eig >= -1e-9 |max eig| (converges to the SDP answer) or the
LP becomes infeasible. Repeated at (g, delta) = (32, .3), (8, .3), (32, .1)
and with / without the epsilon terms, and also on the general 30-dim jet.

A4. u = timelike unit eigenvector of N = M eta (u^T eta u = -1),
h_cov = eta + 2 (eta u)(eta u)^T, I1_h = sum_{mu<nu} eta^mu eta^nu
tr(h_cov F h_cov F^T). d_mu u by first-order perturbation
(N - lam) du = -(dN - dlam) u, dlam = (eta u)^T dN u / ((eta u)^T u),
with u^T eta du = 0; checked against finite differences.
Degeneracy path: M(t) = M_vac + t (e0 e1^T + e1 e0^T), N = M eta has the
2x2 block [[-g, t], [-t, 1]] whose eigenvalues collide at t* = (g+1)/2.

A5. P_t = u u^T eta, Pgrad = sum_mu eta^{mu mu} q(d_mu P_t, d_mu P_t) with
q(X, X) = tr(X X^T) (Frobenius) and q_eta(X, X) = tr(eta X eta X^T); the
omega^2 coefficient on the bare vacuum for each clock.

Out: ../data/m5_32_r1_audit_symbolic.json
"""
from __future__ import annotations

import importlib.util
import itertools
import json
import os
import sys
import time

import numpy as np
import sympy as sp
from scipy.linalg import expm, null_space
from scipy.optimize import linprog

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT = os.path.join(DATA, "m5_32_r1_audit_symbolic.json")


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


REG = _load("m5_32_lagrangian", "m5_32_lagrangian.py")

ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
e = np.diag(ETA)  # metric signs
T0 = time.time()
RES = {"conventions": {
    "static_uniform_time_row": "A_0 = 0 and A_i[0,:] = A_i[:,0] = 0 (only the spatial 3x3 block of M varies)",
    "channel_a0": "a0 = G M_vac + M_vac G^T, G Lorentz generator (boost_k: G[0k]=G[k0]=1; rot_k: G[ij]=-eps_kij); local clock = coms(Gamma_0, d) notebook style",
    "spatial_jet": "A_i = coms(Gamma_i, d), Gamma_i notebook generator (t_i symmetric time row, r_i antisymmetric spatial block): 18 parameters; also the general symmetric 30-parameter jet",
    "H2": "H2(c) = sum_k c_k Q_k, Q_k = omega^2 Lagrangian coefficient matrix of I_k; c_I1 = -4 fixed; feasibility means H2(c) PSD on every channel",
    "epsilon_slot_metric": "epsilon slot on a derivative slot: delta; on an internal slot: eta (lower internal first)",
    "vacuum": "M_vac = diag(-s g, 1, delta, 0), s = -1 -> diag(g, 1, delta, 0)",
}}


def log(*a):
    print(f"[{time.time() - T0:7.1f}s]", *a, flush=True)


# ---------------- curvature + invariants ----------------
def F_of_A(A):
    """A (4,4,4) -> F[mu,nu,a,b]."""
    P = np.einsum("mab,bc,ncd->mnad", A, ETA, A)
    return P - P.transpose(1, 0, 2, 3)


def inv_even(F):
    """the six registered even invariants, explicit metrics."""
    I1 = 0.5 * np.einsum("m,n,a,b,mnab,mnab->", e, e, e, e, F, F)
    I1f = 0.5 * np.einsum("m,n,mnab,mnab->", e, e, F, F)
    I2 = np.einsum("mnab,abmn->", F, F)
    I3 = np.einsum("m,b,mnab,manb->", e, e, F, F)
    R = np.einsum("mnam->na", F)
    I4 = np.einsum("n,a,na,na->", e, e, R, R)
    I5 = np.einsum("na,an->", R, R)
    I6 = np.trace(R) ** 2
    return {"I1": I1, "I1_frob": I1f, "I2": I2, "I3": I3, "I4": I4, "I5": I5, "I6": I6}


def check_registry(rng):
    """my einsum invariants vs the registry K-matrices on random jets."""
    worst = 0.0
    for _ in range(20):
        A = sym(rng.normal(size=(4, 4, 4)))
        F = F_of_A(A)
        mine = inv_even(F)
        for k, v in mine.items():
            ref = float(REG.density_from_K(REG.F_of_A(A[:, None])[0], REG.REGISTRY[k]._K()))
            worst = max(worst, abs(v - ref) / max(1.0, abs(ref)))
    return worst


def sym(A):
    return 0.5 * (A + A.swapaxes(-1, -2))


# ---------------- epsilon contractions (A1) ----------------
EPS = np.zeros((4, 4, 4, 4))
for perm in itertools.permutations(range(4)):
    EPS[perm] = np.linalg.det(np.eye(4)[list(perm)])
SLOT_KIND = "ddiiddii"  # F (x) F


def eps_contractions():
    """all 210 one-epsilon contractions as (label, K 256x256 matrix)."""
    out = []
    for eps_slots in itertools.combinations(range(8), 4):
        rest = [s for s in range(8) if s not in eps_slots]
        pairings = [((rest[0], rest[1]), (rest[2], rest[3])),
                    ((rest[0], rest[2]), (rest[1], rest[3])),
                    ((rest[0], rest[3]), (rest[1], rest[2]))]
        for pr in pairings:
            K = np.zeros((4,) * 8)
            for idx in itertools.product(range(4), repeat=8):
                w = EPS[tuple(idx[s] for s in eps_slots)]
                if w == 0:
                    continue
                for s in eps_slots:
                    if SLOT_KIND[s] == "i":
                        w *= e[idx[s]]
                ok = True
                for (p, q) in pr:
                    if idx[p] != idx[q]:
                        ok = False
                        break
                    kinds = SLOT_KIND[p] + SLOT_KIND[q]
                    if kinds in ("dd", "ii"):
                        w *= e[idx[p]]
                if not ok:
                    continue
                K[idx] += w
            out.append((f"eps{eps_slots}_pair{pr}", K.reshape(256, 256)))
    return out


def E_named(F):
    """the producer's E1, E2, E3 in my conventions."""
    R = np.einsum("mnam->na", F)
    # R = R_a^a: the mixed (d,i) trace is a delta trace under the rule
    # (the producer's "eta-trace" wording: the eta-weighted trace is NOT
    # in the span of the 210 rule-contractions, checked; plain trace is)
    Rtr = np.trace(R)
    # E1 = eps^{mnab} F_{mnab} R: internal a,b lowered with eta
    E1 = np.einsum("mnab,a,b,mnab->", EPS, e, e, F) * Rtr
    # E2 = eps^{mn a r} F_{mnab} R_r^b : a lowered (eta), b (i,i) eta,
    # r = (eps upper, R derivative covariant) -> delta
    E2 = np.einsum("mnar,a,mnab,b,rb->", EPS, e, F, e, R)
    # E3 = eps^{mn a c} F_{mna}^r R_{rc}: a, c lowered (eta), r (i,d) delta
    E3 = np.einsum("mnac,a,c,mnar,rc->", EPS, e, e, F, R)
    return np.array([E1, E2, E3])


# ---------------- Lorentz + jets ----------------
def gen_boost(k):
    G = np.zeros((4, 4)); G[0, k] = G[k, 0] = 1.0; return G


def gen_rot(k):
    G = np.zeros((4, 4))
    i, j = [x for x in (1, 2, 3) if x != k]
    s = 1.0 if (k, i, j) in [(1, 2, 3), (2, 3, 1), (3, 1, 2)] else -1.0
    G[i, j] = -s; G[j, i] = s
    return G


def gamma(t, r):
    Gm = np.zeros((4, 4))
    Gm[0, 1:] = t; Gm[1:, 0] = t
    Gm[1, 2], Gm[1, 3] = -r[2], r[1]
    Gm[2, 1], Gm[2, 3] = r[2], -r[0]
    Gm[3, 1], Gm[3, 2] = -r[1], r[0]
    return Gm


def coms(A, B):
    return A @ ETA @ B - B @ ETA @ A


def vac(g, delta, s=-1.0):
    return np.diag([-s * g, 1.0, delta, 0.0])


def channels(M):
    ch = {}
    for k in (1, 2, 3):
        G = gen_boost(k); ch[f"boost_{k}"] = G @ M + M @ G.T
        G = gen_rot(k); ch[f"rot_{k}"] = G @ M + M @ G.T
    t = np.array([1.0, 0.7, -0.4]); r = np.array([0.5, -0.8, 0.3])
    ch["clock_tr"] = coms(gamma(t, r), M)
    ch["clock_t"] = coms(gamma(t, 0 * r), M)
    ch["clock_r"] = coms(gamma(0 * t, r), M)
    return ch


def jet18(x, M):
    """x (18,) -> A_i, i=1..3 via coms(Gamma_i, d)."""
    A = np.zeros((4, 4, 4))
    for i in range(3):
        A[1 + i] = coms(gamma(x[6 * i:6 * i + 3], x[6 * i + 3:6 * i + 6]), M)
    return A


IU = np.triu_indices(4)


def jet30(x, M):
    A = np.zeros((4, 4, 4))
    for i in range(3):
        S = np.zeros((4, 4)); S[IU] = x[10 * i:10 * i + 10]; A[1 + i] = sym(S) * 2 - np.diag(np.diag(S))
    return A


# ---------------- omega^2 forms (A3) ----------------
def densities(A, extra=None):
    F = F_of_A(A)
    d = inv_even(F)
    if extra is not None:
        d.update(extra(F))
    return d


def omega2_coeffs(a0, Aspace, extra=None):
    """C_k for every term at omega = +-1 (exact quadratic)."""
    def dens(w):
        A = Aspace.copy(); A[0] = w * a0
        return densities(A, extra)
    dp, dm, d0 = dens(1.0), dens(-1.0), dens(0.0)
    return {k: 0.5 * (dp[k] + dm[k]) - d0[k] for k in dp}


def Q_matrices(a0, M, jet, n, extra=None):
    """18x18 (or 30x30) matrices of the quadratic forms C_k by polarization."""
    names = None
    Cd = {}
    for i in range(n):
        x = np.zeros(n); x[i] = 1.0
        c = omega2_coeffs(a0, jet(x, M), extra)
        if names is None:
            names = list(c); Q = {k: np.zeros((n, n)) for k in names}
        for k in names:
            Q[k][i, i] = c[k]
        Cd[i] = c
    for i in range(n):
        for j in range(i + 1, n):
            x = np.zeros(n); x[i] = 1.0; x[j] = 1.0
            c = omega2_coeffs(a0, jet(x, M), extra)
            for k in names:
                Q[k][i, j] = Q[k][j, i] = 0.5 * (c[k] - Cd[i][k] - Cd[j][k])
    return Q


def check_polarization(Q, a0, M, jet, n, rng, extra=None):
    worst = 0.0
    for _ in range(5):
        x = rng.normal(size=n)
        c = omega2_coeffs(a0, jet(x, M), extra)
        for k in Q:
            worst = max(worst, abs(x @ Q[k] @ x - c[k]) / max(1.0, abs(c[k])))
    return worst


def feasibility(Qs, free, fixed, rng, nsamp=2000, tag=""):
    """Qs: {channel: {term: Q}}; free: list of free term names; fixed:
    {term: coef}. LP: for each sampled x and channel, sum_free c_k x^T Q_k x
    >= -sum_fixed c_k x^T Q_k x. Cutting planes on the min eigenvector."""
    n = next(iter(next(iter(Qs.values())).values())).shape[0]
    rows, rhs = [], []
    for ch, Q in Qs.items():
        X = rng.normal(size=(nsamp, n))
        X /= np.linalg.norm(X, axis=1, keepdims=True)
        fx = np.einsum("si,ij,sj->s", X, sum(c * Q[k] for k, c in fixed.items()), X)
        A = np.stack([np.einsum("si,ij,sj->s", X, Q[k], X) for k in free], axis=1)
        rows.append(-A); rhs.append(fx)      # -A c <= fx  <=> A c + fx >= 0
    Aub = np.concatenate(rows); bub = np.concatenate(rhs)
    bounds = [(-1e4, 1e4)] * len(free)
    out = {"free": free, "fixed": fixed, "nsamp": nsamp}
    for it in range(60):
        lp = linprog(np.zeros(len(free)), A_ub=Aub, b_ub=bub, bounds=bounds, method="highs")
        if lp.status != 0:
            out.update({"status": "INFEASIBLE", "lp_message": lp.message, "iterations": it,
                        "n_constraints": int(Aub.shape[0])})
            # also report the most-violated pair of channels: the farkas-like diagnostic
            return out
        c = lp.x
        worst = 0.0; cuts = []
        for ch, Q in Qs.items():
            H = sum(cc * Q[k] for k, cc in fixed.items()) + sum(ci * Q[k] for ci, k in zip(c, free))
            wv, vv = np.linalg.eigh(0.5 * (H + H.T))
            if wv[0] < worst:
                worst = wv[0]
            if wv[0] < -1e-9 * max(1.0, abs(wv[-1])):     # relative PSD tolerance
                v = vv[:, 0]
                fx = v @ sum(cc * Q[k] for k, cc in fixed.items()) @ v
                cuts.append((-np.array([v @ Q[k] @ v for k in free]), fx))
        if not cuts:
            out.update({"status": "FEASIBLE", "c": dict(zip(free, map(float, c))),
                        "min_eig_all_channels": float(worst), "iterations": it})
            return out
        Aub = np.concatenate([Aub, np.stack([a for a, _ in cuts])])
        bub = np.concatenate([bub, np.array([b for _, b in cuts])])
    out.update({"status": "UNDECIDED", "min_eig_all_channels": float(worst), "iterations": 60})
    return out


def compact_certificate(Qs, free, fixed, rng, nsamp=3000):
    """max t s.t. x^T H2(c) x >= t ||x||^2 on samples + cutting planes: t* < 0
    is the best achievable worst-case eigenvalue on the sampled jets (an
    UPPER bound on the true min over c of the min eigenvalue); the active
    constraints at the optimum are a small explicit certificate."""
    n = next(iter(next(iter(Qs.values())).values())).shape[0]
    rows, rhs, tags = [], [], []
    for ch, Q in Qs.items():
        X = rng.normal(size=(nsamp, n)); X /= np.linalg.norm(X, axis=1, keepdims=True)
        fx = np.einsum("si,ij,sj->s", X, sum(c * Q[k] for k, c in fixed.items()), X)
        A = np.stack([np.einsum("si,ij,sj->s", X, Q[k], X) for k in free], axis=1)
        rows.append(np.concatenate([-A, np.ones((nsamp, 1))], axis=1)); rhs.append(fx)
        tags += [(ch, x) for x in X]
    Aub = np.concatenate(rows); bub = np.concatenate(rhs)
    scale = np.max(np.abs(Aub[:, :-1])) + np.max(np.abs(bub))
    Aub[:, :-1] /= scale; bub = bub / scale
    obj = np.zeros(len(free) + 1); obj[-1] = -1.0
    bounds = [(-1e3, 1e3)] * len(free) + [(None, None)]
    for it in range(40):
        lp = linprog(obj, A_ub=Aub, b_ub=bub, bounds=bounds, method="highs")
        c = lp.x[:-1]; t = lp.x[-1]
        cuts = []
        for ch, Q in Qs.items():
            H = sum(cc * Q[k] for k, cc in fixed.items()) + sum(ci * Q[k] for ci, k in zip(c, free))
            wv, vv = np.linalg.eigh(0.5 * (H + H.T))
            if wv[0] < t * scale - 1e-6 * scale:
                v = vv[:, 0]
                fx = v @ sum(cc * Q[k] for k, cc in fixed.items()) @ v
                cuts.append((np.concatenate([-np.array([v @ Q[k] @ v for k in free]) / scale, [1.0]]), fx / scale, (ch, v)))
        if not cuts:
            break
        Aub = np.concatenate([Aub, np.stack([a for a, _, _ in cuts])])
        bub = np.concatenate([bub, np.array([b for _, b, _ in cuts])])
        tags += [tg for _, _, tg in cuts]
    slack = bub - Aub @ lp.x
    act = np.where(slack < 1e-9)[0]
    cert = []
    for i in act[:12]:
        ch, x = tags[i]
        Q = Qs[ch]
        row = {"channel": ch, "fixed_part": float(x @ sum(cc * Q[k] for k, cc in fixed.items()) @ x)}
        row.update({k: float(x @ Q[k] @ x) for k in free})
        row["jet_rounded"] = [round(float(v), 3) for v in x]
        cert.append(row)
    return {"best_worst_eigenvalue_t": float(t * scale), "c_at_optimum": dict(zip(free, map(float, c))),
            "iterations": it, "active_constraints": cert}


def channel_certificates(Qs, free, fixed):
    """per channel: the LP with samples restricted to one channel + the
    minimal-eigenvalue of -4 Q_I1 alone and of the pure-Frobenius choice."""
    out = {}
    for ch, Q in Qs.items():
        H1 = sum(c * Q[k] for k, c in fixed.items())
        out[ch] = {"min_eig_I1_only": float(np.linalg.eigvalsh(H1)[0]),
                   "max_eig_I1_only": float(np.linalg.eigvalsh(H1)[-1]),
                   "min_eig_I1frob_x_minus4": float(np.linalg.eigvalsh(-4 * Q["I1_frob"])[0]) if "I1_frob" in Q else None}
    return out


# ---------------- I1_h + eigenvector machinery (A4, A5) ----------------
def timelike_u(M):
    N = M @ ETA
    w, V = np.linalg.eig(N)
    best = None
    for i in range(4):
        if abs(w[i].imag) > 1e-9:
            continue
        v = V[:, i].real
        q = v @ ETA @ v
        if q < 0:
            if best is None or q < best[0]:
                best = (q, w[i].real, v)
    if best is None:
        return None, None
    q, lam, v = best
    return v / np.sqrt(-q), lam


def du_perturb(M, dM):
    u, lam = timelike_u(M)
    N = M @ ETA; dN = dM @ ETA
    wl = ETA @ u                       # left eigenvector
    dlam = (wl @ dN @ u) / (wl @ u)
    # solve (N - lam) du = -(dN - dlam) u with u^T eta du = 0 (least squares, augmented)
    Aug = np.vstack([N - lam * np.eye(4), (ETA @ u)[None]])
    b = np.concatenate([-(dN - dlam * np.eye(4)) @ u, [0.0]])
    du = np.linalg.lstsq(Aug, b, rcond=None)[0]
    return du, dlam


def du_fd(M, dM, h=1e-6):
    up, _ = timelike_u(M + h * dM); um, _ = timelike_u(M - h * dM)
    if up @ um < 0:
        um = -um
    return (up - um) / (2 * h)


def I1_h(F, M):
    u, _ = timelike_u(M)
    if u is None:
        return np.nan
    v = ETA @ u
    h = ETA + 2 * np.outer(v, v)
    tot = 0.0
    for mu in range(4):
        for nu in range(mu + 1, 4):
            Fm = F[mu, nu]
            tot += e[mu] * e[nu] * np.trace(h @ Fm @ h @ Fm.T)
    return tot


def lorentz(rng, scale=0.3):
    G = np.zeros((4, 4))
    v = scale * rng.normal(size=3); G[0, 1:] = v; G[1:, 0] = v
    r = scale * rng.normal(size=3)
    G[1, 2], G[2, 1] = -r[2], r[2]; G[1, 3], G[3, 1] = r[1], -r[1]; G[2, 3], G[3, 2] = -r[0], r[0]
    L = expm(G)
    assert np.allclose(L.T @ ETA @ L, ETA, atol=1e-10)
    return L


def transform(L, M, A):
    """M -> L M L^T, A_mu -> (L^{-T})_mu^nu L A_nu L^T (d covariant)."""
    Linv = np.linalg.inv(L)
    A2 = np.einsum("nm,nab->mab", Linv, A)         # A'_mu = (L^-1)^nu_mu A_nu
    A2 = np.einsum("ab,mbc,dc->mad", L, A2, L)
    return L @ M @ L.T, A2


def main():
    rng = np.random.default_rng(20260827)
    RES["registry_crosscheck_relerr"] = check_registry(rng)
    log("registry cross-check", RES["registry_crosscheck_relerr"])

    # ================= A1 =================
    log("A1: enumerating 210 one-epsilon contractions ...")
    Ks = eps_contractions()
    njet = 300
    vals = np.zeros((len(Ks), njet)); named = np.zeros((3, njet))
    for j in range(njet):
        A = sym(rng.normal(size=(4, 4, 4)))
        F = F_of_A(A); f = F.reshape(256)
        vals[:, j] = [f @ K @ f for _, K in Ks]
        named[:, j] = E_named(F)
    sv = np.linalg.svd(vals, compute_uv=False)
    rank = int(np.sum(sv > 1e-8 * sv[0]))
    svn = np.linalg.svd(named, compute_uv=False)
    rank_named = int(np.sum(svn > 1e-8 * svn[0]))
    # do E1..E3 span the full space? residual of the 210 rows on span(E)
    coef, *_ = np.linalg.lstsq(named.T, vals.T, rcond=None)
    resid = np.linalg.norm(vals.T - named.T @ coef) / np.linalg.norm(vals)
    # static uniform-time-row: do they vanish? and parity flip
    stat_max = 0.0; par_max = 0.0; nonzero_static = 0.0
    P = np.diag([1.0, -1.0, 1.0, 1.0])
    for j in range(50):
        A = np.zeros((4, 4, 4)); A[1:, 1:, 1:] = sym(rng.normal(size=(3, 3, 3)))
        F = F_of_A(A); f = F.reshape(256)
        stat_max = max(stat_max, max(abs(f @ K @ f) for _, K in Ks))
        nonzero_static = max(nonzero_static, abs(inv_even(F)["I1"]))
        A = sym(rng.normal(size=(4, 4, 4)))
        Ap = np.einsum("nm,nab->mab", P, A); Ap = np.einsum("ab,mbc,dc->mad", P, Ap, P)
        F = F_of_A(A); Fp = F_of_A(Ap)
        par_max = max(par_max, np.max(np.abs(E_named(F) + E_named(Fp))))
    RES["A1"] = {"n_contractions": len(Ks), "rank_210_on_300_jets": rank,
                 "top_singular_values": [float(x) for x in sv[:6]],
                 "rank_E1E2E3": rank_named, "sv_E": [float(x) for x in svn],
                 "relative_residual_of_210_on_span(E1,E2,E3)": float(resid),
                 "max_abs_on_static_uniform_time_row (50 fields)": float(stat_max),
                 "I1_scale_on_those_static_fields": float(nonzero_static),
                 "max_abs_E+E_parity (50 jets)": float(par_max)}
    fourth = None
    if rank > rank_named or resid > 1e-6:
        # name a contraction outside span(E)
        r2 = np.linalg.norm(vals.T - named.T @ coef, axis=0) / np.maximum(np.linalg.norm(vals, axis=1), 1e-300)
        i = int(np.argmax(r2)); fourth = {"label": Ks[i][0], "relative_residual": float(r2[i])}
    RES["A1"]["fourth_invariant"] = fourth
    log("A1 rank", rank, "rank(E)", rank_named, "resid", resid, "static", stat_max, "parity", par_max)

    # ================= A2 =================
    def Nvals(d):
        return {"N1": 4 * d["I3"] - d["I2"] - 2 * d["I1"],
                "N2": 4 * d["I4"] - d["I6"] - 2 * d["I1"],
                "N3": 4 * d["I5"] - d["I2"] - d["I6"]}
    res_static = {"N1": 0.0, "N2": 0.0, "N3": 0.0}; scale = 0.0
    for j in range(50):
        A = np.zeros((4, 4, 4)); A[1:, 1:, 1:] = sym(rng.normal(size=(3, 3, 3)))
        d = inv_even(F_of_A(A)); scale = max(scale, abs(d["I1"]))
        for k, v in Nvals(d).items():
            res_static[k] = max(res_static[k], abs(v))
    # also: a static field whose time row varies (NOT uniform) -> do they vanish? (probe the reading)
    res_static_full = {"N1": 0.0, "N2": 0.0, "N3": 0.0}
    for j in range(20):
        A = sym(rng.normal(size=(4, 4, 4))); A[0] = 0
        d = inv_even(F_of_A(A))
        for k, v in Nvals(d).items():
            res_static_full[k] = max(res_static_full[k], abs(v))
    # witnesses: vacuum + clock + random jet
    M = vac(32.0, 0.3)
    wit = {}
    for name, a0 in channels(M).items():
        A = jet18(rng.normal(size=18), M); A[0] = a0
        d = inv_even(F_of_A(A))
        wit[name] = {k: float(v) for k, v in Nvals(d).items()}
        wit[name]["I1"] = float(d["I1"])
    # is the static-preserving subspace exactly 3-dim mod I1? rank of the 6 invariants on static fields
    Vs = np.array([[inv_even(F_of_A(_A))[k] for k in ("I1", "I2", "I3", "I4", "I5", "I6")]
                   for _A in [np.pad(sym(rng.normal(size=(3, 3, 3))), ((1, 0), (1, 0), (1, 0))) for _ in range(200)]])
    svs = np.linalg.svd(Vs, compute_uv=False)
    ns = null_space(Vs, rcond=1e-9)
    RES["A2"] = {"max_residual_static_uniform_time_row": res_static, "I1_scale": float(scale),
                 "max_residual_static_but_time_row_varying": res_static_full,
                 "witness_nonvanishing": wit,
                 "rank_of_I1..I6_on_static_fields": int(np.sum(svs > 1e-9 * svs[0])),
                 "nullspace_dim_on_static_fields": int(ns.shape[1]),
                 "singular_values": [float(x) for x in svs]}
    log("A2", res_static, "rank on static", RES["A2"]["rank_of_I1..I6_on_static_fields"])

    # ================= A3 =================
    def extra_eps(F):
        E = E_named(F)
        return {"E1": E[0], "E2": E[1], "E3": E[2]}

    A3 = {}
    for (g, delta) in [(32.0, 0.3), (8.0, 0.3), (32.0, 0.1)]:
        M = vac(g, delta)
        key = f"g={g:g},delta={delta:g}"
        Qs = {}
        polw = 0.0
        for name, a0 in channels(M).items():
            Qs[name] = Q_matrices(a0, M, jet18, 18, extra_eps)
            polw = max(polw, check_polarization(Qs[name], a0, M, jet18, 18, rng, extra_eps))
        fixed = {"I1": -4.0}
        free = ["I2", "I3", "I4", "I5", "I6"]
        r = {"polarization_check": float(polw),
             "per_channel": channel_certificates(Qs, free, fixed)}
        r["LP_no_eps_all_channels"] = feasibility(Qs, free, fixed, rng)
        r["LP_with_eps_all_channels"] = feasibility(Qs, free + ["E1", "E2", "E3"], fixed, rng)
        # per-channel feasibility (which channels alone are feasible?)
        r["LP_no_eps_per_channel"] = {ch: feasibility({ch: Q}, free, fixed, rng, nsamp=1000)["status"] for ch, Q in Qs.items()}
        # channel subsets: clock alone, boosts alone, rotations alone
        r["LP_no_eps_subsets"] = {
            "boosts_only": feasibility({k: Qs[k] for k in Qs if k.startswith("boost")}, free, fixed, rng)["status"],
            "rots_only": feasibility({k: Qs[k] for k in Qs if k.startswith("rot")}, free, fixed, rng)["status"],
            "clock_tr+boost_3": feasibility({k: Qs[k] for k in ("clock_tr", "boost_3")}, free, fixed, rng)["status"],
            "six_lorentz_only": feasibility({k: Qs[k] for k in Qs if not k.startswith("clock")}, free, fixed, rng)["status"],
        }
        # the producer's stated certificates on the clock channel: probe my forms on those jet directions
        # (a) Gamma-tilde_11 on boost_3: x index t_1 of A_1 = x[0]; (Gamma-tilde_11 + Gamma-tilde_22) = x[0] + x[7]
        Qb = Qs["boost_3"]
        def lin(Q, x):
            return {k: float(x @ Q[k] @ x) for k in ("I1", "I2", "I3", "I4", "I5", "I6")}
        x = np.zeros(18); x[0] = 1.0
        r["boost_3_dir_Gt11"] = lin(Qb, x)
        x = np.zeros(18); x[0] = 1.0; x[7] = 1.0
        r["boost_3_dir_Gt11+Gt22"] = lin(Qb, x)
        if key == "g=32,delta=0.3":
            r["compact_certificate_all_channels"] = compact_certificate(Qs, free, fixed, rng)
            r["compact_certificate_boost_3"] = compact_certificate({"boost_3": Qs["boost_3"]}, free, fixed, rng)
            r["compact_certificate_clock_tr"] = compact_certificate({"clock_tr": Qs["clock_tr"]}, free, fixed, rng)
            # LP machinery sanity: with I1_frob fixed at -4 (PSD on every channel) the LP must be FEASIBLE
            r["LP_sanity_I1frob_fixed"] = feasibility(Qs, free, {"I1_frob": -4.0}, rng)["status"]
        # the I1_frob-only choice (variant A) as a form: PSD on every channel?
        r["minus4_I1frob_min_eig_per_channel"] = {ch: float(np.linalg.eigvalsh(-4 * Q["I1_frob"])[0]) for ch, Q in Qs.items()}
        # also I1_frob expressed inside span(I1..I6)? (is variant A reachable by c2..c6?)
        A3[key] = r
        log("A3", key, "no-eps:", r["LP_no_eps_all_channels"]["status"], "with-eps:", r["LP_with_eps_all_channels"]["status"],
            "subsets:", r["LP_no_eps_subsets"])
    # general 30-dim jet at the toy point
    M = vac(32.0, 0.3)
    Qs30 = {name: Q_matrices(a0, M, jet30, 30) for name, a0 in channels(M).items()}
    A3["g=32,delta=0.3,general_30_jet"] = {
        "LP_no_eps_all_channels": feasibility(Qs30, ["I2", "I3", "I4", "I5", "I6"], {"I1": -4.0}, rng)["status"],
        "minus4_I1frob_min_eig_per_channel": {ch: float(np.linalg.eigvalsh(-4 * Q["I1_frob"])[0]) for ch, Q in Qs30.items()}}
    # free c1 as well (is the no-go about the -4 or structural?)
    A3["g=32,delta=0.3,c1_plus4"] = feasibility(
        {n: Q_matrices(a0, M, jet18, 18) for n, a0 in channels(M).items()},
        ["I2", "I3", "I4", "I5", "I6"], {"I1": 4.0}, rng)
    RES["A3"] = A3
    log("A3 general30", A3["g=32,delta=0.3,general_30_jet"]["LP_no_eps_all_channels"], "c1=+4", A3["g=32,delta=0.3,c1_plus4"]["status"])

    # ================= A4 =================
    M = vac(32.0, 0.3)
    A4 = {}
    # (i) covariance on random jets, random M near vacuum, random boosts
    cov = 0.0
    for _ in range(50):
        Mr = M + 0.5 * sym(rng.normal(size=(4, 4)))
        A = sym(rng.normal(size=(4, 4, 4)))
        L = lorentz(rng)
        M2, A2 = transform(L, Mr, A)
        v1 = I1_h(F_of_A(A), Mr); v2 = I1_h(F_of_A(A2), M2)
        cov = max(cov, abs(v1 - v2) / max(1.0, abs(v1)))
    A4["covariance_relerr_max"] = float(cov)
    # sanity: I1 also invariant, I1_frob NOT
    inv1 = 0.0; invf = 0.0
    for _ in range(20):
        A = sym(rng.normal(size=(4, 4, 4))); L = lorentz(rng)
        M2, A2 = transform(L, M, A)
        d1 = inv_even(F_of_A(A)); d2 = inv_even(F_of_A(A2))
        inv1 = max(inv1, abs(d1["I1"] - d2["I1"]) / max(1, abs(d1["I1"])))
        invf = max(invf, abs(d1["I1_frob"] - d2["I1_frob"]) / max(1, abs(d1["I1_frob"])))
    A4["control_I1_invariance"] = float(inv1); A4["control_I1_frob_noninvariance"] = float(invf)
    # (ii) at vacuum equal I1_frob
    dd = 0.0
    for _ in range(20):
        A = sym(rng.normal(size=(4, 4, 4))); F = F_of_A(A)
        dd = max(dd, abs(I1_h(F, M) - inv_even(F)["I1_frob"]))
    A4["vacuum_I1h_minus_I1frob_max"] = float(dd)
    # (iii) static uniform-time-row: equal I1 (off shell, M varies in spatial block)
    dd = 0.0
    for _ in range(30):
        Ms = M.copy(); Ms[1:, 1:] += 2.0 * sym(rng.normal(size=(3, 3)))
        A = np.zeros((4, 4, 4)); A[1:, 1:, 1:] = sym(rng.normal(size=(3, 3, 3)))
        F = F_of_A(A)
        dd = max(dd, abs(I1_h(F, Ms) - inv_even(F)["I1"]))
    A4["static_uniform_time_row_I1h_minus_I1_max"] = float(dd)
    # (iv) H2 = -4 C form on every channel at vacuum (I1_h == I1_frob at the point)
    def extra_h(F):
        return {"I1_h": I1_h(F, M)}
    eigs = {}
    for name, a0 in channels(M).items():
        Q = Q_matrices(a0, M, jet18, 18, extra_h)
        w = np.linalg.eigvalsh(-4 * Q["I1_h"])
        eigs[name] = {"min_eig": float(w[0]), "max_eig": float(w[-1]),
                      "I1h_vs_I1frob_matrix_maxdiff": float(np.max(np.abs(Q["I1_h"] - Q["I1_frob"])))}
    A4["H2_minus4_I1h_per_channel"] = eigs
    # du: perturbation vs finite difference
    dmax = 0.0
    for _ in range(20):
        Mr = M + 0.5 * sym(rng.normal(size=(4, 4))); dM = sym(rng.normal(size=(4, 4)))
        du, _ = du_perturb(Mr, dM); dfd = du_fd(Mr, dM)
        dmax = max(dmax, np.max(np.abs(du - dfd)))
    A4["du_perturbation_vs_fd_maxdiff"] = float(dmax)
    # degeneracy locus: M(t) = vac + t (e0 e1^T + e1 e0^T); t* = (g+1)/2
    g = 32.0; tstar = (g + 1) / 2
    path = []
    A = sym(rng.normal(size=(4, 4, 4)))                 # fixed jet along the path
    Ab = np.zeros((4, 4, 4)); Ab[0] = channels(M)["boost_1"]  # boost clock alone
    for frac in [0.0, 0.5, 0.9, 0.99, 0.999, 0.9999, 1.0, 1.001, 1.01]:
        t = frac * tstar
        Mt = M.copy(); Mt[0, 1] = Mt[1, 0] = t
        u, lam = timelike_u(Mt)
        ev = np.linalg.eigvals(Mt @ ETA)
        row = {"t/t*": frac, "t": t, "eigs_Meta": [complex(x).real if abs(complex(x).imag) < 1e-12 else str(complex(x)) for x in ev]}
        if u is None:
            row["u"] = None; row["I1_h_random_jet"] = None; row["I1_h_boost_clock"] = None
        else:
            row["u_norm_euclid"] = float(np.linalg.norm(u)); row["lam"] = float(lam)
            row["I1_h_random_jet"] = float(I1_h(F_of_A(A), Mt))
            row["I1_random_jet"] = float(inv_even(F_of_A(A))["I1"])
            row["I1_h_boost_clock"] = float(I1_h(F_of_A(Ab), Mt))
        path.append(row)
    # sign of the divergent part near the locus over random jets: min and max of I1_h
    Mt = M.copy(); Mt[0, 1] = Mt[1, 0] = 0.9999 * tstar
    vals_near = [I1_h(F_of_A(sym(rng.normal(size=(4, 4, 4)))), Mt) for _ in range(200)]
    A4["degeneracy_path"] = path
    A4["degeneracy_t_star"] = tstar
    A4["I1h_near_locus_200_random_jets"] = {"min": float(min(vals_near)), "max": float(max(vals_near))}
    # a second locus: timelike eigenvalue crossing another real eigenvalue via a diagonal path (g -> 1)
    path2 = []
    for gg in [32.0, 4.0, 1.5, 1.01, 1.0, 0.99, 0.5]:
        Mt = vac(gg, 0.3)
        u, lam = timelike_u(Mt)
        path2.append({"g": gg, "lam": None if lam is None else float(lam), "u": None if u is None else [float(x) for x in u],
                      "I1_h_random_jet": None if u is None else float(I1_h(F_of_A(A), Mt))})
    A4["diagonal_crossing_path_g_to_1"] = path2
    RES["A4"] = A4
    log("A4 cov", cov, "vac", A4["vacuum_I1h_minus_I1frob_max"], "static", A4["static_uniform_time_row_I1h_minus_I1_max"],
        "du", dmax, "min eig", {k: v["min_eig"] for k, v in eigs.items()})

    # ================= A5 =================
    A5 = {}
    for name, a0 in channels(M).items():
        du, _ = du_perturb(M, a0)
        u, _ = timelike_u(M)
        dP = np.outer(du, u) @ ETA + np.outer(u, du) @ ETA    # d_0 P_t at omega = 1
        q_frob = float(np.sum(dP * dP)); q_eta = float(np.trace(ETA @ dP @ ETA @ dP.T))
        A5[name] = {"omega2_coeff_Pgrad_frob": -q_frob, "omega2_coeff_Pgrad_eta": -q_eta,
                    "du_norm": float(np.linalg.norm(du))}
    # cross-check with a finite difference along the clock orbit: u(omega tau) for boost_1
    L = expm(0.01 * gen_boost(1)); M2 = L @ M @ L.T
    u1, _ = timelike_u(M); u2, _ = timelike_u(M2)
    if u1 @ u2 < 0:
        u2 = -u2
    log("A5", {k: v["omega2_coeff_Pgrad_frob"] for k, v in A5.items()})
    A5["fd_check_boost_1_du_norm"] = float(np.linalg.norm((u2 - u1) / 0.01))
    RES["A5"] = A5

    RES["runtime_s"] = time.time() - T0
    with open(OUT, "w") as f:
        json.dump(RES, f, indent=1, default=float)
    log("wrote", OUT)


if __name__ == "__main__":
    main()
