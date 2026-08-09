"""M5.21.14 E: INDEPENDENT ADVERSARIAL AUDIT.

Attacks the M5.21.14 claims with the auditor's OWN implementations:
nothing here imports or reuses the task's a/b/c/d scripts. Allowed
loads: the certified instruments m5_21_3_a_4d.py (e_parts, kin_of,
base_cfg, coords, sym4) and m5_21_8_b_lattice.py (dressed, a0_unit),
which predate this task; the task's JSON outputs are read as the
claims under test.

Auditor routes (all different from the task's):
  C1  T1 forms: complex-step differentiation (h = 1e-20, exact to
      machine precision) of the auditor's own dressed-field builder on
      the auditor's own random smooth fields, then RICHARDSON
      EXTRAPOLATION IN 1/g (two-point, kills the 1/g term): the
      extrapolated dressed-minus-base density must hit the T1 formula
      to < 1e-6 relative. Both s signs; beta -> -beta evenness decay.
  C2  E_V invariance: sympy EXACT (symbolic b, two rational unit
      directions) + scipy expm as an independent Qb construction +
      mpmath 50-digit trace powers at full b (up to |b| = 2).
  C3  Unboundedness: sympy exact W = 0 for parallel v's; the negative
      channel nonvanishing on the hedgehog (auditor's algebraic
      hedgehog, no angle functions); descent ladder on the auditor's
      own oscillatory family beta_k = sin(kr)/sqrt(k)*window with an
      oscillation-resolved GL panel quadrature.
  C4  Sawtooth floor: own exact dressed evaluation at g = 32
      (complex-step + GL panels resolving the sawtooth), spot checks
      vs the BND record + own lam exponent.
  C5  Bulk kin flip: own radial-shell evaluation of base and dressed
      kin at the recorded avec profile, window r in (10, 20).
  C6  Realizable gain: certified lattice instruments, auditor's own
      Qb field assembly (K2 = K @ K, not hand-assembled), R2 profile
      vs rigid constant-b scan.
  C7  h-attribution arithmetic from the JSONs.
  C8  Retrodiction arithmetic from the stored V2 curves (own
      parabola interpolation).
  C9  wording review (structured notes).

The analytic family is rebuilt ALGEBRAICALLY (no atan2, complex-step
safe): M3 = n n^T + delta*phi phi^T with phi = (-y, x, 0)/rho, and
a0 = delta*(phi psi^T + psi phi^T), psi = (-xz/(r rho), -yz/(r rho),
rho/r); equality with the certified Rodrigues-angle builder is gated
at 1e-10 before use.

Run: python3 m5_21_14_e_audit.py [c1 c2 ... c9 | all]
Out: ../data/m5_21_14_audit.json (merged per stage)
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT = os.path.join(DATA, "m5_21_14_audit.json")

ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
DELTA = 0.3
CS_H = 1e-20                       # complex-step size (subtraction-free)


# ================= auditor eta algebra (own code) =================
def dens_u(A):
    """4 sum_{i<j} tr(eta F eta F^T), F = Ai eta Aj - Aj eta Ai."""
    tot = 0.0
    for i in range(3):
        for j in range(i + 1, 3):
            F = A[i] @ ETA @ A[j] - A[j] @ ETA @ A[i]
            EF = ETA @ F @ ETA
            tot = tot + 4.0 * np.einsum("...ab,...ab->...", EF, F)
    return tot


def dens_k(a0, A):
    tot = 0.0
    for i in range(3):
        F = a0 @ ETA @ A[i] - A[i] @ ETA @ a0
        EF = ETA @ F @ ETA
        tot = tot + 4.0 * np.einsum("...ab,...ab->...", EF, F)
    return tot


def qb_of(n, b):
    """Qb = I + sinh(b) K + (cosh(b)-1) K@K, K assembled from n; the
    square is COMPUTED (K @ K), not hand-assembled."""
    N = n.shape[0]
    K = np.zeros((N, 4, 4), dtype=n.dtype)
    K[:, 0, 1:] = n
    K[:, 1:, 0] = n
    K2 = K @ K
    return (np.eye(4, dtype=n.dtype)[None] + np.sinh(b)[:, None, None]
            * K + (np.cosh(b) - 1.0)[:, None, None] * K2)


def cs_grads(build, P):
    """complex-step spatial derivatives of a (N,4,4) or (N,3,3)
    matrix field: exact to machine precision, no subtraction."""
    A = []
    for ax in range(3):
        Pc = P.astype(complex).copy()
        Pc[:, ax] = Pc[:, ax] + 1j * CS_H
        A.append(build(Pc).imag / CS_H)
    return A


# ================= auditor T1 formula (own implementation) ==========
def t1_static(G, V):
    tot = 0.0
    for i in range(3):
        for j in range(i + 1, 3):
            C = G[i] @ G[j] - G[j] @ G[i]
            W = (V[j][..., :, None] * V[i][..., None, :]
                 - V[i][..., :, None] * V[j][..., None, :])
            w = (np.einsum("...ab,...b->...a", G[i], V[j])
                 - np.einsum("...ab,...b->...a", G[j], V[i]))
            tot = tot + 4.0 * (
                2.0 * np.einsum("...ab,...ab->...", C, W)
                + np.einsum("...ab,...ab->...", W, W)
                - 2.0 * np.einsum("...a,...a->...", w, w))
    return tot


def t1_static_parts(G, V):
    """(cross, quartic, negative) channel totals, for C3."""
    cr = qu = ng = 0.0
    for i in range(3):
        for j in range(i + 1, 3):
            C = G[i] @ G[j] - G[j] @ G[i]
            W = (V[j][..., :, None] * V[i][..., None, :]
                 - V[i][..., :, None] * V[j][..., None, :])
            w = (np.einsum("...ab,...b->...a", G[i], V[j])
                 - np.einsum("...ab,...b->...a", G[j], V[i]))
            cr = cr + 8.0 * np.einsum("...ab,...ab->...", C, W)
            qu = qu + 4.0 * np.einsum("...ab,...ab->...", W, W)
            ng = ng - 8.0 * np.einsum("...a,...a->...", w, w)
    return cr, qu, ng


def t1_kin(Md, V):
    tot = 0.0
    for i in range(3):
        u = np.einsum("...ab,...b->...a", Md, V[i])
        tot = tot - 8.0 * np.einsum("...a,...a->...", u, u)
    return tot


# ================= auditor algebraic analytic family ================
def rn_of(P):
    r = np.sqrt(P[:, 0] ** 2 + P[:, 1] ** 2 + P[:, 2] ** 2)
    return r, P / r[:, None]


def hh_m3(P, delta=DELTA):
    """spatial hedgehog block, algebraic: n n^T + delta phi phi^T."""
    r, n = rn_of(P)
    rho = np.sqrt(P[:, 0] ** 2 + P[:, 1] ** 2)
    phi = np.stack([-P[:, 1] / rho, P[:, 0] / rho,
                    np.zeros_like(rho)], axis=-1)
    return (n[:, :, None] * n[:, None, :]
            + delta * phi[:, :, None] * phi[:, None, :])


def hh_m4(P, g, delta=DELTA):
    """full base, s = -1 branch: time slot +g."""
    N = P.shape[0]
    M = np.zeros((N, 4, 4), dtype=P.dtype)
    M[:, 0, 0] = g
    M[:, 1:, 1:] = hh_m3(P, delta)
    return M


def hh_a0(P, delta=DELTA):
    """unit-omega clock flow at t = 0, algebraic:
    a0_spatial = delta*(phi psi^T + psi phi^T)."""
    r, n = rn_of(P)
    rho = np.sqrt(P[:, 0] ** 2 + P[:, 1] ** 2)
    phi = np.stack([-P[:, 1] / rho, P[:, 0] / rho,
                    np.zeros_like(rho)], axis=-1)
    psi = np.stack([-P[:, 0] * P[:, 2] / (r * rho),
                    -P[:, 1] * P[:, 2] / (r * rho), rho / r], axis=-1)
    N = P.shape[0]
    A = np.zeros((N, 4, 4), dtype=P.dtype)
    A[:, 1:, 1:] = delta * (phi[:, :, None] * psi[:, None, :]
                            + psi[:, :, None] * phi[:, None, :])
    return A


def hh_m4_dressed(P, g, bfun, delta=DELTA):
    r, n = rn_of(P)
    Q = qb_of(n, bfun(r))
    return Q @ hh_m4(P, g, delta) @ np.swapaxes(Q, -1, -2)


# ================= auditor quadrature (GL panels) =================
def gl_panels(r0, r1, width, ngl):
    x, w = np.polynomial.legendre.leggauss(ngl)
    edges = np.linspace(r0, r1, max(2, int(np.ceil((r1 - r0) / width))
                                    + 1))
    rs, ws = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        rs.append(0.5 * (b - a) * x + 0.5 * (a + b))
        ws.append(0.5 * (b - a) * w)
    return np.concatenate(rs), np.concatenate(ws)


def sphere_dirs(nu, nphi):
    u, wu = np.polynomial.legendre.leggauss(nu)
    phis = (np.arange(nphi) + 0.5) * 2 * np.pi / nphi
    st = np.sqrt(1 - u ** 2)
    dirs = np.stack([np.outer(st, np.cos(phis)).ravel(),
                     np.outer(st, np.sin(phis)).ravel(),
                     np.repeat(u, nphi)], axis=-1)
    wd = np.repeat(wu, nphi) * (2 * np.pi / nphi)
    return dirs, wd


def ball_grid(r0, r1, width, ngl, nu, nphi):
    rs, wr = gl_panels(r0, r1, width, ngl)
    dirs, wd = sphere_dirs(nu, nphi)
    P = (rs[:, None, None] * dirs[None]).reshape(-1, 3)
    W = (wr[:, None] * rs[:, None] ** 2 * wd[None]).ravel()
    return P, W


# ================= profile families (reimplemented) =================
RHOS9 = np.geomspace(0.5, 16.0, 9)
RES_RHOS = [4.3666, 6.7420, 10.4097, 16.0]   # d-script family constant


def b_of(avec, r):
    val = avec[0] * np.tanh(r / 2.0)
    for k, rho in enumerate(RHOS9):
        val = val + avec[k + 1] * (r / rho) * np.exp(-((r / rho) ** 2))
    return val


def taper(r):
    rr = np.real(r)
    t = np.ones_like(rr)
    t[rr >= 20.0] = 0.0
    band = (rr > 14.0) & (rr < 20.0)
    t = np.where(band, 0.5 * (1 + np.cos(np.pi
                                         * np.clip(rr - 14.0, 0, 6)
                                         / 6.0)), t)
    return t


def b_res_of(p, r):
    val = p[0] * np.tanh(r / 2.0)
    for k, rho in enumerate(RES_RHOS):
        val = val + p[k + 1] * (r / rho) * np.exp(-((r / rho) ** 2))
    return val * taper(r)


def saw_of(A, lam, r):
    mask = (np.real(r) <= 8.0)
    return A * np.sin(np.pi * r / lam) * mask


# ================= C1: the T1 forms =================
def stage_c1():
    rng = np.random.default_rng(977)
    KS = rng.normal(0.0, 0.4, size=(4, 3))
    PH = rng.uniform(0, 2 * np.pi, size=4)
    AM = rng.normal(0.0, 0.35, size=(4, 3, 3))
    AM = 0.5 * (AM + AM.transpose(0, 2, 1))
    PHD = rng.uniform(0, 2 * np.pi, size=4)
    AMD = rng.normal(0.0, 0.35, size=(4, 3, 3))
    AMD = 0.5 * (AMD + AMD.transpose(0, 2, 1))

    def m3f(P):
        M = np.zeros(P.shape[:1] + (3, 3), dtype=P.dtype)
        for k in range(4):
            M = M + AM[k] * np.cos(P @ KS[k] + PH[k])[:, None, None]
        return M

    def md3f(P):
        M = np.zeros(P.shape[:1] + (3, 3), dtype=P.dtype)
        for k in range(4):
            M = M + AMD[k] * np.cos(P @ KS[k] + PHD[k])[:, None, None]
        return M

    def betaf(r):
        return 0.7 * r * r * np.exp(-((r / 2.1) ** 2))

    def bnf(P):
        r, n = rn_of(P)
        return betaf(r)[:, None] * n

    def pad34(M3, P):
        M = np.zeros(P.shape[:1] + (4, 4), dtype=M3.dtype)
        M[:, 1:, 1:] = M3
        return M

    def m4d(P, g, s, sb):
        r, n = rn_of(P)
        M = pad34(m3f(P), P)
        M[:, 0, 0] = -s * g
        Q = qb_of(n, sb * betaf(r) / g)
        return Q @ M @ np.swapaxes(Q, -1, -2)

    def a0d(P, g, sb):
        r, n = rn_of(P)
        Q = qb_of(n, sb * betaf(r) / g)
        return Q @ pad34(md3f(P), P) @ np.swapaxes(Q, -1, -2)

    P = np.array([[1.1, -0.6, 0.8], [-1.7, 0.5, 1.2],
                  [0.9, 1.4, -1.1]])
    # the T1 reference (auditor formula, complex-step ingredients)
    G = [Gi[:, :, :] for Gi in cs_grads(m3f, P)]

    def bn_c(Pc):
        r, n = rn_of(Pc)
        return betaf(r)[:, None] * n
    Vfull = []
    for ax in range(3):
        Pc = P.astype(complex).copy()
        Pc[:, ax] += 1j * CS_H
        Vfull.append(bn_c(Pc).imag / CS_H)
    Md3 = md3f(P)
    T1u = t1_static(G, Vfull)
    T1k = t1_kin(Md3, Vfull)

    # base densities
    Ab = cs_grads(lambda Q: pad34(m3f(Q), Q), P)
    dub = dens_u(Ab)
    dkb = dens_k(pad34(md3f(P), P), Ab)

    gs = [1e3, 1e4, 1e5]
    rows = {}
    for s in (-1.0, 1.0):
        Du, Dk = {}, {}
        for g in gs:
            A = cs_grads(lambda Q, gg=g: m4d(Q, gg, s, +1), P)
            Du[g] = dens_u(A) - dub
            Dk[g] = dens_k(a0d(P, g, +1), A) - dkb
        eu = {g: np.max(np.abs(Du[g] - T1u) / np.abs(T1u))
              for g in gs}
        ek = {g: np.max(np.abs(Dk[g] - T1k) / np.abs(T1k))
              for g in gs}
        g1, g2 = 1e4, 1e5
        exu = (g2 * Du[g2] - g1 * Du[g1]) / (g2 - g1)
        exk = (g2 * Dk[g2] - g1 * Dk[g1]) / (g2 - g1)
        rows[f"s{int(s):+d}"] = {
            "err_u_by_g": {f"{g:g}": float(eu[g]) for g in gs},
            "err_k_by_g": {f"{g:g}": float(ek[g]) for g in gs},
            "slope_ratio_u": [float(eu[1e3] / eu[1e4]),
                              float(eu[1e4] / eu[1e5])],
            "extrap_res_u": float(np.max(np.abs(exu - T1u)
                                         / np.abs(T1u))),
            "extrap_res_k": float(np.max(np.abs(exk - T1k)
                                         / np.abs(T1k)))}
    # evenness: beta -> -beta at g = 1e3 and 1e4 (odd part ~ 1/g)
    ev = {}
    for g in (1e3, 1e4):
        Ap = cs_grads(lambda Q, gg=g: m4d(Q, gg, -1.0, +1), P)
        Am = cs_grads(lambda Q, gg=g: m4d(Q, gg, -1.0, -1), P)
        dp = dens_u(Ap) - dub
        dm = dens_u(Am) - dub
        ev[f"{g:g}"] = float(np.max(np.abs(dp - dm) / np.abs(T1u)))
    # s-difference decay
    sd = {}
    for g in (1e3, 1e4):
        Am1 = cs_grads(lambda Q, gg=g: m4d(Q, gg, -1.0, +1), P)
        Ap1 = cs_grads(lambda Q, gg=g: m4d(Q, gg, +1.0, +1), P)
        sd[f"{g:g}"] = float(np.max(np.abs(dens_u(Am1)
                                           - dens_u(Ap1))
                                    / np.abs(T1u)))
    # T1 formula evenness is exact by construction: verify anyway
    t1_even = float(np.max(np.abs(
        t1_static(G, [-v for v in Vfull]) - T1u)))
    ok = all(r["extrap_res_u"] < 1e-6 and r["extrap_res_k"] < 1e-6
             and all(8.0 < x < 12.5 for x in r["slope_ratio_u"])
             for r in rows.values())
    ok = ok and ev["10000"] < 0.2 * ev["1000"] + 1e-12 \
        and sd["10000"] < 0.2 * sd["1000"] + 1e-12 \
        and t1_even < 1e-12
    return {"verdict": "CONFIRMED" if ok else "REFUTED",
            "T1u_at_points": T1u.tolist(),
            "T1k_at_points": T1k.tolist(),
            "rows": rows, "evenness_odd_part_rel": ev,
            "s_sign_diff_rel": sd,
            "t1_formula_evenness_abs": t1_even,
            "sharpening": ("the measured odd part is 0.0 EXACTLY at "
                           "finite g: b -> -b is conjugation by "
                           "J = diag(1,-1,-1,-1) (JKJ = -K, JM4J = "
                           "M4 for block-diagonal M4), which the "
                           "trace densities cannot see; the evenness "
                           "is exact at ALL orders, stronger than "
                           "the note's leading-order claim"),
            "route": ("complex-step derivatives (exact), own fields, "
                      "own builder, 1/g Richardson extrapolation")}


# ================= C2: E_V invariance =================
def stage_c2():
    import sympy as sp
    from scipy.linalg import expm
    out = {}
    # exact symbolic, two rational unit directions, symbolic b
    b = sp.symbols("b", real=True)
    checks = []
    for nvec in ((sp.Rational(3, 7), sp.Rational(6, 7),
                  sp.Rational(2, 7)),
                 (sp.Rational(2, 11), sp.Rational(6, 11),
                  sp.Rational(9, 11))):
        K = sp.zeros(4, 4)
        for i in range(3):
            K[0, i + 1] = nvec[i]
            K[i + 1, 0] = nvec[i]
        Qb = (sp.eye(4) + sp.sinh(b) * K
              + (sp.cosh(b) - 1) * (K * K))
        eta = sp.diag(-1, 1, 1, 1)
        Z = sp.simplify(Qb * eta * Qb.T - eta)
        checks.append(Z == sp.zeros(4, 4))
    out["symbolic_eta_exact"] = all(checks)
    # expm route (independent construction of Qb)
    rng = np.random.default_rng(31415)
    worst_q, worst_eta = 0.0, 0.0
    for _ in range(20):
        n = rng.normal(size=3)
        n = n / np.linalg.norm(n)
        bb = float(rng.uniform(-2.0, 2.0))
        K = np.zeros((4, 4))
        K[0, 1:] = n
        K[1:, 0] = n
        Qe = expm(bb * K)
        Qc = qb_of(n[None], np.array([bb]))[0]
        worst_q = max(worst_q, float(np.max(np.abs(Qe - Qc))))
        worst_eta = max(worst_eta,
                        float(np.max(np.abs(Qe @ ETA @ Qe.T - ETA))))
    out["expm_vs_closed_form_max"] = worst_q
    out["expm_eta_residual_max"] = worst_eta
    # mpmath 50-digit trace powers at full b, random symmetric M
    import mpmath as mp
    mp.mp.dps = 50
    worst_tr = mp.mpf(0)
    rng = np.random.default_rng(2718)
    for _ in range(5):
        n = rng.normal(size=3)
        bb = float(rng.uniform(-2.0, 2.0))
        Ms = rng.normal(size=(4, 4))
        Ms = 0.5 * (Ms + Ms.T)
        # normalize n INSIDE mpmath so K^3 = K holds to 50 digits
        nm = [mp.mpf(v) for v in n]
        nrm = mp.sqrt(sum(v * v for v in nm))
        nm = [v / nrm for v in nm]
        Km = mp.zeros(4)
        for i in range(3):
            Km[0, i + 1] = nm[i]
            Km[i + 1, 0] = nm[i]
        Q = mp.eye(4) + mp.sinh(bb) * Km + (mp.cosh(bb) - 1) * Km * Km
        E = mp.diag([-1, 1, 1, 1])
        Mm = mp.matrix(Ms.tolist())
        Md = Q * Mm * Q.T
        A0, A1 = Mm * E, Md * E
        P0, P1 = mp.eye(4), mp.eye(4)
        for p in range(1, 5):
            P0, P1 = P0 * A0, P1 * A1
            t0 = sum(P0[i, i] for i in range(4))
            t1 = sum(P1[i, i] for i in range(4))
            worst_tr = max(worst_tr,
                           abs(t1 - t0) / max(abs(t0), mp.mpf(1)))
    out["trace_powers_rel_dev_50dig"] = float(worst_tr)
    ok = (out["symbolic_eta_exact"] and worst_q < 1e-12
          and worst_eta < 1e-12 and out["trace_powers_rel_dev_50dig"]
          < 1e-40)
    out["verdict"] = "CONFIRMED" if ok else "REFUTED"
    out["note"] = ("pointwise similarity => E_V invariant for ANY "
                   "b(r): the invariance is pointwise, so it holds "
                   "under any profile and any quadrature")
    return out


# ================= C3: unboundedness =================
def stage_c3():
    import sympy as sp
    out = {}
    # (a) W = 0 for v_i parallel to a common direction, symbolic
    c1, c2 = sp.symbols("c1 c2")
    nn = sp.Matrix(sp.symbols("q1 q2 q3"))
    v1, v2 = c1 * nn, c2 * nn
    W = v2 * v1.T - v1 * v2.T
    out["W_zero_parallel_symbolic"] = sp.expand(W) == sp.zeros(3, 3)
    # negative channel nonvanishing on the hedgehog with radial v's
    P = np.array([[0.9, 0.7, 1.1], [-1.3, 0.4, 0.6]])
    G = cs_grads(lambda Q: hh_m3(Q), P)
    r, n = rn_of(P)
    Vrad = [n[:, ax][:, None] * n for ax in range(3)]   # beta' = 1
    negsum = 0.0
    wq = 0.0
    for i in range(3):
        for j in range(i + 1, 3):
            w = (np.einsum("nab,nb->na", G[i], Vrad[j])
                 - np.einsum("nab,nb->na", G[j], Vrad[i]))
            negsum += 2.0 * np.einsum("na,na->n", w, w)
            W = (Vrad[j][:, :, None] * Vrad[i][:, None, :]
                 - Vrad[i][:, :, None] * Vrad[j][:, None, :])
            wq += np.einsum("nab,nab->n", W, W)
    out["neg_channel_hedgehog_radial_v"] = np.atleast_1d(
        negsum).tolist()
    out["quartic_W2_radial_v_max"] = float(np.max(np.abs(wq)))
    # (b) descent ladder, auditor family
    ks = [4.0, 8.0, 16.0, 32.0, 64.0]
    Pg, Wg = ball_grid(1.0, 7.0, 2 * np.pi / 64.0 / 6.0, 4, 10, 8)
    Gg = cs_grads(lambda Q: hh_m3(Q), Pg)
    rg, ng = rn_of(Pg)
    ladder = {}
    for k in ks:
        def betak(r, kk=k):
            return (np.sin(kk * r) / np.sqrt(kk)
                    * np.exp(-(((r - 4.0) / 1.5) ** 2)))

        def bnk(Pc, kk=k):
            rr, nn2 = rn_of(Pc)
            return betak(rr, kk)[:, None] * nn2
        V = []
        for ax in range(3):
            Pc = Pg.astype(complex).copy()
            Pc[:, ax] += 1j * CS_H
            V.append(bnk(Pc).imag / CS_H)
        cr, qu, ngc = t1_static_parts(Gg, V)
        ladder[f"k{k:g}"] = {
            "T1_total": float(np.sum(Wg * (cr + qu + ngc))),
            "cross": float(np.sum(Wg * cr)),
            "quartic": float(np.sum(Wg * qu)),
            "negative": float(np.sum(Wg * ngc))}
    Ts = [ladder[f"k{k:g}"]["T1_total"] for k in ks]
    out["descent_ladder"] = ladder
    out["descent_monotone"] = bool(all(Ts[i + 1] < Ts[i]
                                       for i in range(len(Ts) - 1)))
    out["descent_growth_T64_over_T4"] = float(Ts[-1] / Ts[0])
    ok = (out["W_zero_parallel_symbolic"]
          and out["quartic_W2_radial_v_max"] < 1e-25
          and min(np.atleast_1d(negsum)) > 1e-3
          and out["descent_monotone"] and Ts[-1] < 8 * Ts[0]
          and Ts[-1] < 0)
    out["verdict"] = "CONFIRMED" if ok else "REFUTED"
    return out


# ================= C4: sawtooth floor spot checks =================
def stage_c4():
    with open(os.path.join(DATA, "m5_21_14_minimize.json")) as f:
        MJ = json.load(f)
    bnd = {row["lam"]: dict(zip(row["amps"], row["E"]))
           for row in MJ["BND"]}
    g = 32.0
    pts = [(2.0, 0.02), (1.0, 0.02), (0.5, 0.02), (0.5, 0.1)]
    rows = []
    for lam, A in pts:
        Pg, Wg = ball_grid(0.05, 8.05, lam / 6.0, 5, 14, 8)

        def bfun(r, A=A, lam=lam):
            return saw_of(A, lam, r)
        Ad = cs_grads(lambda Q: hh_m4_dressed(Q, g, bfun), Pg)
        Ab = cs_grads(lambda Q: hh_m4(Q, g), Pg)
        E = float(np.sum(Wg * (dens_u(Ad) - dens_u(Ab))))
        claimed = bnd[lam][A]
        rows.append({"lam": lam, "A": A, "E_audit": E,
                     "E_claimed": claimed,
                     "rel_dev": float(abs(E - claimed)
                                      / max(abs(claimed), 1e-12)),
                     "same_sign": bool(np.sign(E)
                                       == np.sign(claimed))})
        print(json.dumps(rows[-1]), flush=True)
    # own lam exponent at A = 0.02 from audit numbers
    e02 = {r["lam"]: r["E_audit"] for r in rows if r["A"] == 0.02}
    lams = sorted(e02)
    slope = float(np.polyfit(np.log([lams[0], lams[-1]]),
                             np.log([abs(e02[lams[0]]),
                                     abs(e02[lams[-1]])]), 1)[0])
    exp_claim = float(np.polyfit(
        np.log([0.5, 2.0]),
        np.log([abs(bnd[0.5][0.02]), abs(bnd[2.0][0.02])]), 1)[0])
    ok_sign = all(r["same_sign"] for r in rows)
    ok_mag = all(0.1 < abs(r["E_audit"] / r["E_claimed"]) < 10
                 for r in rows)
    return {"verdict": "CONFIRMED" if (ok_sign and ok_mag)
            else ("PARTIAL" if ok_sign else "REFUTED"),
            "rows": rows,
            "lam_exponent_audit": slope,
            "lam_exponent_claimed_from_BND": exp_claim,
            "note": ("audit quadrature resolves the sawtooth "
                     "(GL panels lam/6); the task grid is geomspace "
                     "72-pt, under-resolved for lam = 0.5 at r > 3.5")}


# ================= C5: bulk kin slopes =================
def stage_c5():
    with open(os.path.join(DATA, "m5_21_14_minimize.json")) as f:
        MJ = json.load(f)
    avec = np.array(MJ["avec"])
    rs_task = np.array(MJ["rs"])
    g = 32.0
    win = rs_task[(rs_task > 10.0) & (rs_task < 20.0)]
    dirs, wd = sphere_dirs(16, 12)
    k_base, k_corr = [], []
    for rr in win:
        P = rr * dirs

        def bfun(r):
            return b_of(avec, r)
        Ad = cs_grads(lambda Q: hh_m4_dressed(Q, g, bfun), P)
        Ab = cs_grads(lambda Q: hh_m4(Q, g), P)
        a0b = hh_a0(P)
        r_, n_ = rn_of(P)
        Q0 = qb_of(n_, bfun(r_))
        a0d = Q0 @ a0b @ np.swapaxes(Q0, -1, -2)
        dkb = dens_k(a0b, Ab)
        dkd = dens_k(a0d, Ad)
        k_base.append(rr ** 2 * float(np.sum(wd * dkb)))
        k_corr.append(rr ** 2 * float(np.sum(wd * (dkd - dkb))))
    base = float(np.mean(k_base))
    corr = float(np.mean(k_corr))
    claimed = {"base": 13.481202393263455,
               "dressing": -19.604374521818748,
               "net": -6.123172128555293}
    ok = (abs(base - claimed["base"]) / abs(claimed["base"]) < 0.15
          and abs(corr - claimed["dressing"])
          / abs(claimed["dressing"]) < 0.15
          and (base + corr) < 0)
    return {"verdict": "CONFIRMED" if ok else "REFUTED",
            "audit_slopes": {"base": base, "dressing": corr,
                             "net": base + corr},
            "claimed_slopes": claimed,
            "window_radii": win.tolist(),
            "note": ("k(r) = r^2 * angular integral of the kin "
                     "density, averaged over the task's own window "
                     "radii; auditor fields + auditor quadrature")}


# ================= C6: realizable gain on the instrument ===========
def _load_certified():
    sp8 = importlib.util.spec_from_file_location(
        "b8x", os.path.join(HERE, "m5_21_8_b_lattice.py"))
    B8 = importlib.util.module_from_spec(sp8)
    sp8.loader.exec_module(B8)
    return B8.INS4, B8


def _lattice_qb_apply(INS4, B8, cfg, bfun):
    """auditor's own Qb field assembly on the certified grid."""
    n, h = cfg["n"], cfg["h"]
    X, Y, Z = INS4.coords(n, h)
    R = np.sqrt(X * X + Y * Y + Z * Z)
    bl = bfun(R.ravel()).reshape(R.shape)
    nv = np.stack([X / R, Y / R, Z / R], axis=-1)
    K = np.zeros(X.shape + (4, 4))
    K[..., 0, 1:] = nv
    K[..., 1:, 0] = nv
    K2 = K @ K
    Qb = (np.eye(4) + np.sinh(bl)[..., None, None] * K
          + (np.cosh(bl) - 1.0)[..., None, None] * K2)
    Mb = B8.dressed(cfg, 0.0)
    a0b = B8.a0_unit(cfg, 0.0)
    Md = INS4.sym4(Qb @ Mb @ np.swapaxes(Qb, -1, -2))
    a0d = INS4.sym4(Qb @ a0b @ np.swapaxes(Qb, -1, -2))
    return Mb, a0b, Md, a0d


def stage_c6():
    with open(os.path.join(DATA, "m5_21_14_resolution.json")) as f:
        RJ = json.load(f)
    INS4, B8 = _load_certified()
    cfg = INS4.base_cfg(s=-1.0, g=32.0, n=32, L=48.0)
    pres = np.array(RJ["R2"]["params"])
    Mb, a0b, Md, a0d = _lattice_qb_apply(
        INS4, B8, cfg, lambda r: b_res_of(pres, r))
    eu0, ev0 = INS4.e_parts(Mb, cfg)
    eud, evd = INS4.e_parts(Md, cfg)
    e_corr = float(eud + evd - eu0 - ev0)
    k_corr = float(INS4.kin_of(Md, a0d, cfg)
                   - INS4.kin_of(Mb, a0b, cfg))
    # rigid constant-b gain: fine scan + parabola
    mh = float(np.arctanh(1.0 / 32.0))
    ms = np.linspace(0.0, 1.4 * mh, 57)
    E0 = float(sum(INS4.e_parts(Mb, cfg)))
    Es = []
    for m in ms:
        if m == 0.0:
            Es.append(E0)
            continue
        _, _, Mdm, _ = _lattice_qb_apply(
            INS4, B8, cfg, lambda r, mm=m: np.full_like(r, mm))
        Es.append(float(sum(INS4.e_parts(Mdm, cfg))))
    Es = np.array(Es)
    i = int(np.argmin(Es))
    a, bq, c = Es[i - 1], Es[i], Es[i + 1]
    dm = ms[1] - ms[0]
    m_star = float(ms[i] - 0.5 * dm * (c - a) / (c - 2 * bq + a))
    E_star = float(bq - 0.125 * (c - a) ** 2 / (c - 2 * bq + a))
    gain_rigid = E_star - E0
    # cross-check: my constant-b dressing == B8.dressed(cfg, m)
    _, _, Mdx, _ = _lattice_qb_apply(
        INS4, B8, cfg, lambda r: np.full_like(r, 0.02))
    xdev = float(np.max(np.abs(Mdx - INS4.sym4(
        B8.dressed(cfg, 0.02)))))
    ratio = e_corr / gain_rigid
    ok = (abs(e_corr - RJ["R2"]["lattice_n32"]["E_corr"])
          / abs(RJ["R2"]["lattice_n32"]["E_corr"]) < 0.02
          and abs(gain_rigid - (-61.11)) / 61.11 < 0.05
          and 2.2 < ratio < 3.0 and xdev < 1e-9)
    return {"verdict": "CONFIRMED" if ok else "REFUTED",
            "E_corr_R2_audit": e_corr,
            "E_corr_R2_claimed": RJ["R2"]["lattice_n32"]["E_corr"],
            "kin_corr_R2_audit": k_corr,
            "rigid_gain_audit": gain_rigid,
            "rigid_m_star": m_star,
            "rigid_gain_claimed_V2_g32": -61.112549184653005,
            "variational_over_rigid": ratio,
            "const_b_equals_B8_dressed_max": xdev}


# ================= C7 + C8: JSON arithmetic =================
def stage_c7():
    with open(os.path.join(DATA, "m5_21_14_resolution.json")) as f:
        RJ = json.load(f)
    d32 = RJ["R4"]["E_rel_dev"]
    d64 = RJ["R4"]["E_rel_dev_n64"]
    k32 = RJ["R4"]["kin_rel_dev"]
    k64 = RJ["R4"]["kin_rel_dev_n64"]
    expo = float(np.log2(d32 / d64))
    ladder = RJ["R3"]["ladder"]
    return {"verdict": "PARTIAL",
            "masked_E_dev_n32_n64": [d32, d64],
            "ratio": float(d32 / d64),
            "implied_exponent": expo,
            "masked_kin_dev_n32_n64": [k32, k64],
            "kin_ratio": float(k32 / k64),
            "R3_full_ladder_E_corr": [r["E_corr"] for r in ladder],
            "caveats": [
                "two h points only: consistent with h^2 (exponent "
                "1.91) but cannot distinguish nearby exponents",
                "the kin deviation drops only ~2.1x (h^1), not h^2",
                "the FULL-profile R3 ladder E_corr is non-monotone "
                "and sign-flipping (+42.5 / -379.1 / +897.6 vs "
                "continuum -4309.5): the discretization attribution "
                "holds only for the r > 3 masked region"]}


def stage_c8():
    with open(os.path.join(DATA, "m5_21_14_verify.json")) as f:
        VJ = json.load(f)
    claimed_ratio = [0.8168, 0.8275, 0.8329]
    m5218_recorded = [0.8177, 0.8287, 0.8328]
    rows = []
    gains = []
    for c in VJ["V2"]["curves"]:
        m = np.array(c["m"])
        E = np.array(c["E"])
        i = int(np.argmin(E))
        a, b, cc = E[i - 1], E[i], E[i + 1]
        dm = m[1] - m[0]
        ms = float(m[i] - 0.5 * dm * (cc - a) / (cc - 2 * b + a))
        Es = float(b - 0.125 * (cc - a) ** 2 / (cc - 2 * b + a))
        E0 = float(E[int(np.argmin(np.abs(m)))])
        mh = float(np.arctanh(1.0 / c["g"]))
        # estimator-independence: quartic polyfit around the minimum
        sl = slice(max(0, i - 4), i + 5)
        co = np.polyfit(m[sl], E[sl], 4)
        rt = np.roots(np.polyder(co))
        rt = rt[np.isreal(rt)].real
        ms4 = float(rt[np.argmin(np.polyval(co, rt))])
        rows.append({"g": c["g"], "ratio_audit": abs(ms) / mh,
                     "ratio_quartic_fit": abs(ms4) / mh,
                     "gain_audit": Es - E0})
        gains.append(-(Es - E0))
    spread = max(gains) / min(gains)
    ok = (all(abs(rows[k]["ratio_audit"] - claimed_ratio[k]) < 2e-3
              for k in range(3))
          and abs(spread - 1.0134) < 0.005
          and all(abs(rows[k]["ratio_audit"] - m5218_recorded[k])
                  < 4e-3 for k in range(3)))
    return {"verdict": "CONFIRMED" if ok else "REFUTED",
            "rows": rows,
            "gain_spread_audit": float(spread),
            "gain_spread_claimed": 1.0134,
            "claimed_ratios": claimed_ratio,
            "m5_21_8_recorded": m5218_recorded,
            "nuance": ("the ratio's 3rd digit is estimator-dependent "
                       "at the +-0.002 level (parabola vs quartic "
                       "fit); the 0.82-0.84 band and the 1.3% gain "
                       "flatness are robust, but 'reproduced to 3 "
                       "digits' should read 'to ~0.002'")}


# ================= C9: wording =================
def stage_c9():
    return {"verdict": "PARTIAL",
            "fair": [
                "'MET at the term level' is a fair reading of "
                "T1_kin = -8 sum|Md3 v_i|^2 <= 0 (identity, both "
                "routes + this audit's C1)",
                "the boundedness row and the section 5 guard do "
                "state the T1_static runaway and the free-"
                "minimization ban where they are needed",
                "the C row states constant-omega oscillations do "
                "not self-start globally (box-limited ledger)"],
            "flags": [
                "section 4 axis B says 'negative ... for any "
                "nonzero dressing': the correct qualifier "
                "(section 1.3 has it) is any dressing with "
                "Md3 v_i != 0 for some i; T1_kin = 0 exactly when "
                "all Md3 v_i = 0, e.g. any dressing supported "
                "where Mdot3 vanishes",
                "negative kin_total means E(omega) descends as "
                "-|kin|*omega^2 without bound unless a stabilizer "
                "(fixed-J / constraint, the M5.21.8 section 6 "
                "convergence) is invoked; the note implies this "
                "via the C row but never states 'needs a "
                "stabilizer' explicitly for the flipped-sign "
                "regime: recommend one sentence in axis C or B"]}


STAGES = {"c1": stage_c1, "c2": stage_c2, "c3": stage_c3,
          "c4": stage_c4, "c5": stage_c5, "c6": stage_c6,
          "c7": stage_c7, "c8": stage_c8, "c9": stage_c9}


def main(argv):
    names = argv or list(STAGES)
    if names == ["all"]:
        names = list(STAGES)
    if os.path.exists(OUT):
        with open(OUT) as f:
            out = json.load(f)
    else:
        out = {"claims": {}}
    for nm in names:
        t0 = time.time()
        res = STAGES[nm]()
        res["stage_runtime_s"] = round(time.time() - t0, 1)
        out["claims"][nm[1]] = res
        print(json.dumps({nm: {"verdict": res["verdict"],
                               "runtime_s": res["stage_runtime_s"]}}),
              flush=True)
        with open(OUT, "w") as f:
            json.dump(out, f, indent=1)
    out["verdicts_summary"] = {k: v["verdict"]
                               for k, v in sorted(out["claims"].items())}
    with open(OUT, "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out.get("verdicts_summary", {})))


if __name__ == "__main__":
    main(sys.argv[1:])
