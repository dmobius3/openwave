"""M5.32 R7 ADVERSARIAL AUDIT: the time-row gradient penalty K_T (class C4).

An INDEPENDENT rebuild of the R7 claims A1-A6 (form level) and B1-B7
(lattice level). The producer's scripts (m5_32_r7_a_kt_form.py and
m5_32_r7_b_kt_lattice.py) were NOT read and are NOT imported; K_T is
re-implemented here from the written definition. Oracles: the certified
stack (m5_21_3_a_4d.py stencil / e_parts / kin_of / comm_eta, and the
m5_21_8_b_lattice.py hedgehog + clock flow) and the earlier auditors'
own modules (m5_32_r2_audit_lattice.py, m5_32_r4_audit_clock.py) which
were written by auditors, not by the R7 producer.

EQUATIONS FIRST
---------------
Signature eta = diag(-1, 1, 1, 1), index 0 = time (both as a derivative
index and as the internal row of M). M is a real symmetric 4x4 field with
RAW CONTRAVARIANT internal entries; jets A_mu = d_mu M; A_0 = omega a0.

    N        = M eta                       (the (1,1) endomorphism)
    u        = timelike unit eigenvector of N, u^T eta u = -1
    h        = eta + 2 (eta u)(eta u)^T     (h_cov, all indices low)
    K_T      = 1/2 sum_mu eta^{mu mu} [ tr(h A_mu h A_mu)
                                        - tr(eta A_mu eta A_mu) ]

Identity proved here (audit_A0, and used throughout): in the u-frame
(the frame where u = e0, so h = 1) and for a SYMMETRIC jet,

    K_T = 2 sum_mu eta^{mu mu} sum_j (A_mu)_{0j}^2
        = 2 [ sum_i sum_j (A_i)_{0j}^2 - omega^2 sum_j (a0)_{0j}^2 ]

because 1 - eta_a eta_b = 2 exactly on the mixed (time, space) index
pairs and 0 elsewhere. Reading it as a Lagrangian term L = ... - c2 K_T
with the Legendre rule H_I = C omega^2 - A gives

    E[-c2 K_T] = c2 ( 2 sum_i sum_j (A_i)_{0j}^2 )      static cost
               + c2 omega^2 ( 2 sum_j (a0)_{0j}^2 )     kinetic floor
both non-negative for c2 > 0 and for a symmetric a0. On an ANTISYMMETRIC
a0 the sign flips, since then (A)_{j0} = -(A)_{0j}.

Lattice: E = E_stat + omega^2 kin with
    E_stat = h^3 sum_x [4 sum_{i<j} q_lam(F_ij) + c2 sum_i k_i] + E_V
    kin    = h^3 sum_x [4 sum_i   q_lam([a0, A_i]_eta)] + c2 kappa
    q_lam  = (1-lam) tr(eta F eta F^T) + lam tr(h F h F^T)
    k_i    = 1/2[tr(h A_i h A_i) - tr(eta A_i eta A_i)]
    kappa  = h^3 sum_x 1/2[tr(h a0 h a0) - tr(eta a0 eta a0)]
Fixed J: E_J = E_stat + J^2/(4 kin), omega* = J/(2 kin).

Modes (each writes a checkpoint under ../checkpoints/m5_32_r7/):
    form   A1..A6
    lat    B1..B6 (own grid, own minimizer)
    b7     B7 (three-to-five box ladder, shell decomposition, delta ladder,
           cell-by-cell boundary-artifact test)
    merge  assembles ../data/m5_32_r7_audit.json
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys
import time

import numpy as np
from scipy.linalg import expm

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
CKPT = os.path.join(HERE, "..", "checkpoints", "m5_32_r7")
os.makedirs(CKPT, exist_ok=True)
T0 = time.time()


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


B3 = _load("m5_21_3_a_4d", "m5_21_3_a_4d.py")          # certified stack
B8 = _load("m5_21_8_b_lattice", "m5_21_8_b_lattice.py")  # hedgehog + clock

ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
G, DELTA, S = 32.0, 0.3, -1.0


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def dump(tag, obj):
    with open(os.path.join(CKPT, f"{tag}.json"), "w") as f:
        json.dump(obj, f, indent=1, default=float)
    log(f"checkpoint {tag}.json written")


def rel(a, b):
    d = max(abs(a), abs(b), 1e-300)
    return abs(a - b) / d


# ===================== my own K_T implementation =====================
def u_timelike(M):
    """MY eigen-solve: the timelike unit eigenvector u of N = M eta,
    u^T eta u = -1, contravariant. Returns u (..., 4) and a degeneracy
    measure (the gap between the timelike eigenvalue and the closest other)."""
    N = M @ ETA
    lam, V = np.linalg.eig(N)
    if np.max(np.abs(lam.imag)) > 1e-9 * max(np.max(np.abs(lam.real)), 1.0):
        raise ValueError("complex spectrum of M eta")
    lam, V = lam.real, V.real
    nrm = np.einsum("...ak,a,...ak->...k", V, np.diag(ETA), V)
    k0 = np.argmin(nrm, axis=-1)                      # most negative = timelike
    V = V / np.sqrt(np.abs(nrm))[..., None, :]
    u = np.take_along_axis(V, k0[..., None, None], axis=-1)[..., 0]
    l0 = np.take_along_axis(lam, k0[..., None], axis=-1)
    gap = np.min(np.where(np.abs(lam - l0) < 1e-14, np.inf, np.abs(lam - l0)),
                 axis=-1)
    return u, gap


def h_from_u(u):
    """h_cov = eta + 2 (eta u)(eta u)^T (all indices low)."""
    v = u @ ETA
    return ETA + 2.0 * v[..., :, None] * v[..., None, :]


def h_cov(M):
    return h_from_u(u_timelike(M)[0])


def tr_XAXA(X, A):
    """tr(X A X A) = X_ab A_bc X_cd A_da."""
    return np.einsum("...ab,...bc,...cd,...da->...", X, A, X, A, optimize=True)


def kt_density(A, h, weights=None, hmat=None):
    """K_T density from the DEFINITION. A: (4, ..., 4, 4) jets, h: h_cov.
    weights: the eta^{mu mu} derivative weights (mutable for mutants).
    hmat: override of the internal metric (mutable for mutants)."""
    w = np.diag(ETA) if weights is None else np.asarray(weights, float)
    hh = h if hmat is None else hmat
    tot = 0.0
    for mu in range(4):
        tot = tot + w[mu] * 0.5 * (tr_XAXA(hh, A[mu]) - tr_XAXA(ETA, A[mu]))
    return tot


def kt_uframe(A, u):
    """The claimed u-frame form 2 sum_mu eta^mumu sum_j (A_mu)_0j (A_mu)_j0,
    evaluated by boosting the jets into the u-frame with MY own boost."""
    Lam = boost_to_e0(u)
    tot = 0.0
    for mu in range(4):
        Ap = np.einsum("...ab,...bc,...dc->...ad", Lam, A[mu], Lam)
        tot = tot + np.diag(ETA)[mu] * 2.0 * np.sum(
            Ap[..., 0, 1:] * Ap[..., 1:, 0], axis=-1)
    return tot


def boost_to_e0(u):
    """MY Lambda in SO(1,3) with Lambda u = e0, built by the closed-form
    boost (no expm): batched over the leading axes."""
    u = np.asarray(u, float)
    sgn = np.sign(u[..., 0])
    sgn = np.where(sgn == 0, 1.0, sgn)
    u = u * sgn[..., None]
    ch = u[..., 0]                                   # cosh(rap)
    vec = u[..., 1:]
    sh = np.linalg.norm(vec, axis=-1)                # sinh(rap)
    nsafe = np.where(sh > 1e-300, sh, 1.0)
    nh = vec / nsafe[..., None]
    Lam = np.zeros(u.shape[:-1] + (4, 4))
    Lam[..., 0, 0] = ch
    Lam[..., 0, 1:] = -sh[..., None] * nh
    Lam[..., 1:, 0] = -sh[..., None] * nh
    eye3 = np.eye(3)
    Lam[..., 1:, 1:] = (eye3 + (ch - 1.0)[..., None, None]
                        * nh[..., :, None] * nh[..., None, :])
    return Lam


# ===================== lattice helpers (my own) =====================
def mydiff(f, ax, h, br):
    """MY forward / backward difference, zero-padded at the far edge."""
    out = np.zeros_like(f)
    sl = [slice(None)] * f.ndim
    def at(i):
        s = list(sl); s[ax] = i; return tuple(s)
    d = (f[at(slice(1, None))] - f[at(slice(0, -1))]) / h
    if br == "fwd":
        out[at(slice(0, -1))] = d
    elif br == "bwd":
        out[at(slice(1, None))] = d
    else:
        raise ValueError(br)
    return out


def q_eta(F):
    return np.einsum("...ab,...cd,ac,bd->...", F, F, ETA, ETA, optimize=True)


def q_h(F, h):
    return np.einsum("...ab,...bc,...cd,...ad->...", h, F, h, F, optimize=True)


def comm_eta(A, B):
    return A @ ETA @ B - B @ ETA @ A


BRANCHES = (("bwd", 0.5), ("fwd", 0.5))      # MY branch order (reversed)


def lattice_pieces(M, a0, h, cfg, want_density=False):
    """(E_u0, E_u1, E_V, kin0, kin1, E_KT, kappa) h^3-weighted, MY assembly.
    E_KT = h^3 sum_x sum_i k_i (the STATIC K_T cost, coefficient c2)
    kappa = h^3 sum_x 1/2[tr(h a0 h a0) - tr(eta a0 eta a0)] (the omega^2 floor)."""
    h3 = cfg["h"] ** 3
    eu0 = eu1 = k0 = k1 = kt = 0.0
    dens_kin = 0.0 if want_density else None
    for br, wt in BRANCHES:
        A = [mydiff(M, ax, cfg["h"], br) for ax in range(3)]
        for i in range(3):
            for j in range(i + 1, 3):
                F = comm_eta(A[i], A[j])
                eu0 += wt * 4.0 * np.sum(q_eta(F))
                eu1 += wt * 4.0 * np.sum(q_h(F, h))
            if a0 is not None:
                F = comm_eta(a0, A[i])
                dk = q_eta(F)
                k0 += wt * 4.0 * np.sum(dk)
                k1 += wt * 4.0 * np.sum(q_h(F, h))
                if want_density:
                    dens_kin = dens_kin + wt * 4.0 * dk
            kt += wt * 0.5 * np.sum(tr_XAXA(h, A[i]) - tr_XAXA(ETA, A[i]))
    _, ev = B3.e_parts(M, cfg)
    kap = 0.0
    if a0 is not None:
        kap = 0.5 * np.sum(tr_XAXA(h, a0) - tr_XAXA(ETA, a0))
    out = (h3 * eu0, h3 * eu1, float(ev), h3 * k0, h3 * k1, h3 * kt, h3 * kap)
    if want_density:
        return out, h3 * dens_kin
    return out


def mix(a, b, lam):
    return (1.0 - lam) * a + lam * b


def assemble(p, lam, c2):
    """(E_stat, kin) from the pieces."""
    eu0, eu1, ev, kk0, kk1, kt, kap = p
    return mix(eu0, eu1, lam) + ev + c2 * kt, mix(kk0, kk1, lam) + c2 * kap


# ===================== the dressing family (my own build) =====================
def qb_field(cfg, b):
    """MY radial boost field with per-cell rapidity b (n,n,n): closed form
    Qb = I + sinh(b) K + (cosh(b)-1) K^2, K = n.(E_0i + E_i0)."""
    X, Y, Z = B3.coords(cfg["n"], cfg["h"])
    r = np.sqrt(X * X + Y * Y + Z * Z)
    nh = np.stack([X / r, Y / r, Z / r], axis=-1)
    K = np.zeros(X.shape + (4, 4))
    K[..., 0, 1:] = nh
    K[..., 1:, 0] = nh
    K2 = np.zeros_like(K)
    K2[..., 0, 0] = 1.0
    K2[..., 1:, 1:] = nh[..., :, None] * nh[..., None, :]
    return (np.eye(4) + np.sinh(b)[..., None, None] * K
            + (np.cosh(b) - 1.0)[..., None, None] * K2)


def conj(Q, X):
    return np.einsum("...ab,...bc,...dc->...ad", Q, X, Q)


class Fam:
    """The R7 dressing family, rebuilt: b(r) = amp tanh(r/2) exp(-(r/R)^2)."""

    def __init__(self, n, L, delta=DELTA, g=G):
        self.cfg = B3.base_cfg(s=S, g=g, n=n, L=L, delta=delta)
        self.Mb = B8.dressed(self.cfg, 0.0)
        self.a0 = B8.a0_unit(self.cfg, 0.0)
        X, Y, Z = B3.coords(n, self.cfg["h"])
        self.r = np.sqrt(X * X + Y * Y + Z * Z)
        self.cache = {}

    def b_of(self, amp, R, mu=1.0):
        rr = self.r / mu
        return amp * np.tanh(rr / 2.0) * np.exp(-((rr / R) ** 2))

    def pieces(self, amp, R, mu=1.0, eig_h=False):
        key = (round(amp, 12), round(R, 9), round(mu, 9), eig_h)
        if key in self.cache:
            return self.cache[key]
        if amp == 0.0:
            Md, a0d, h = self.Mb, self.a0, np.broadcast_to(
                ETA, self.Mb.shape).copy()
            if eig_h:
                h = h_cov(Md)
        else:
            Qb = qb_field(self.cfg, self.b_of(amp, R, mu))
            Md = B3.sym4(conj(Qb, self.Mb))
            a0d = B3.sym4(conj(Qb, self.a0))
            h = h_cov(Md) if eig_h else h_from_u(
                np.einsum("...ab,b->...a", Qb, np.array([1.0, 0, 0, 0])))
        out = tuple(float(v) for v in lattice_pieces(Md, a0d, h, self.cfg))
        self.cache[key] = out
        return out


# ===================== Lorentz toolkit (my own) =====================
def gen_boost(k):
    K = np.zeros((4, 4)); K[0, k] = K[k, 0] = 1.0
    return K


def gen_rot(i, j):
    J = np.zeros((4, 4)); J[i, j] = -1.0; J[j, i] = 1.0
    return J


def lorentz(rap=(0, 0, 0), ang=(0, 0, 0)):
    Gm = sum(rap[k - 1] * gen_boost(k) for k in (1, 2, 3))
    Gm = Gm + ang[0] * gen_rot(1, 2) + ang[1] * gen_rot(2, 3) + ang[2] * gen_rot(1, 3)
    return expm(Gm)


def transform(Lm, M, A):
    """M -> L M L^T, A_mu -> (L^-1)^nu_mu L A_nu L^T (x -> L x)."""
    Mp = Lm @ M @ Lm.T
    Li = np.linalg.inv(Lm)
    inner = np.einsum("ab,nbc,dc->nad", Lm, A, Lm)
    Ap = np.einsum("nm,nab->mab", Li, inner)
    return Mp, Ap


def rand_sym(rng, scale=1.0):
    X = rng.normal(size=(4, 4)) * scale
    return X + X.T


def vac(g=G, delta=DELTA):
    return np.diag([-S * g, 1.0, delta, 0.0])


def tangent(Gm, M):
    """the symmetric tangent of M -> L M L^T along Gm."""
    return Gm @ M + M @ Gm.T


# =============================== A1 ===============================
def audit_A1(rng):
    out = {"transforms": {}, "mutants": {}, "scale_growth": {}}
    base = []
    for _ in range(6):
        M = vac() + 0.4 * rand_sym(rng)
        Aj = np.stack([rand_sym(rng) for _ in range(4)])
        base.append((M, Aj))
    tests = {}
    for rp in (0.5, 1.0, 2.0, 2.5, 3.0):
        tests[f"boost_x_rap{rp}"] = lorentz(rap=(rp, 0, 0))
        tests[f"boost_diag_rap{rp}"] = lorentz(rap=(rp / np.sqrt(3),) * 3)
    for an in (0.7, 2.0):
        tests[f"rot_rap{an}"] = lorentz(ang=(an, 0.3 * an, -0.5 * an))
    for rp in (1.0, 2.5, 3.0):
        tests[f"mixed_rap{rp}"] = lorentz(rap=(rp, 0.4 * rp, -0.2 * rp),
                                          ang=(0.9, -0.6, 0.4))
    worst = 0.0
    for nm, Lm in tests.items():
        w = 0.0
        for M, Aj in base:
            u, _ = u_timelike(M)
            k0 = kt_density(Aj, h_from_u(u))
            Mp, Ap = transform(Lm, M, Aj)
            up, _ = u_timelike(Mp)
            k1 = kt_density(Ap, h_from_u(up))
            w = max(w, rel(k0, k1))
        out["transforms"][nm] = w
        worst = max(worst, w)
    out["worst_drift"] = worst
    # NEGATIVE controls: each must FAIL
    Lm = lorentz(rap=(2.5, 1.0, -0.5), ang=(0.9, -0.6, 0.4))
    for mut in ("no_weights", "h_to_eye", "u_to_e1", "h_frozen"):
        w = 0.0
        for M, Aj in base:
            u, _ = u_timelike(M)
            h0 = h_from_u(u)
            Mp, Ap = transform(Lm, M, Aj)
            up, _ = u_timelike(Mp)
            h1 = h_from_u(up)
            if mut == "no_weights":
                a = kt_density(Aj, h0, weights=(1, 1, 1, 1))
                b = kt_density(Ap, h1, weights=(1, 1, 1, 1))
            elif mut == "h_to_eye":
                a = kt_density(Aj, h0, hmat=np.eye(4))
                b = kt_density(Ap, h1, hmat=np.eye(4))
            elif mut == "u_to_e1":
                hw = h_from_u(np.array([0.0, 1.0, 0.0, 0.0]))
                a = kt_density(Aj, h0, hmat=hw)
                b = kt_density(Ap, h1, hmat=hw)
            else:                       # h frozen at the untransformed value
                a = kt_density(Aj, h0)
                b = kt_density(Ap, h0)
            w = max(w, rel(a, b))
        out["mutants"][mut] = w
    # does the drift grow faster than the intermediate trace scale?
    Lm = lorentz(rap=(3.0, 0, 0))
    rows = []
    for sc in (0.1, 1.0, 10.0, 100.0):
        M = vac() + 0.4 * rand_sym(np.random.default_rng(11))
        Aj = sc * np.stack([rand_sym(np.random.default_rng(12 + q)) for q in range(4)])
        u, _ = u_timelike(M)
        k0 = kt_density(Aj, h_from_u(u))
        Mp, Ap = transform(Lm, M, Aj)
        up, _ = u_timelike(Mp)
        hp = h_from_u(up)
        k1 = kt_density(Ap, hp)
        inter = max(abs(tr_XAXA(hp, Ap[mu])) for mu in range(4))
        rows.append({"jet_scale": sc, "K_T": k0, "abs_drift": abs(k1 - k0),
                     "rel_drift": rel(k0, k1), "intermediate_scale": float(inter),
                     "drift_over_eps_times_intermediate":
                         abs(k1 - k0) / (2.22e-16 * float(inter))})
    out["scale_growth"]["rows"] = rows
    out["scale_growth"]["floor_multiple_max"] = max(
        r["drift_over_eps_times_intermediate"] for r in rows)
    dump("a1_covariance", out)
    return out


# =============================== A2 ===============================
def audit_A2(rng):
    out = {}
    # (i) analytic: for block-diagonal M, N = M eta is block diagonal, e0 is an
    # exact eigenvector with eta-norm -1 (the unique timelike one), so u = +-e0,
    # h = eta + 2 e0 e0^T = 1, and K_T = 2 sum_mu eta^mumu sum_j A_0j A_j0 = 0
    # since A_mu = d_mu M has A_0j = 0. EXACT, no cancellation.
    out["analytic"] = ("block-diagonal M => N = M eta block diagonal => u = e0 "
                       "exactly => h = identity => K_T = 2 sum_mu eta^mumu "
                       "sum_j (A_mu)_0j (A_mu)_j0 = 0 since (A_mu)_0j = 0")
    # (ii) my own random block-diagonal static fields on a small lattice
    n, L = 16, 48.0
    cfg = B3.base_cfg(s=S, g=G, n=n, L=L, delta=DELTA)
    X, Y, Z = B3.coords(n, cfg["h"])
    rr = np.sqrt(X * X + Y * Y + Z * Z)
    fields = []
    for q in range(20):
        Md = np.zeros(X.shape + (4, 4))
        Md[..., 0, 0] = G + rng.normal() * np.cos(rr / (2 + q % 5))
        B = np.zeros(X.shape + (3, 3))
        for a in range(3):
            for b in range(a, 3):
                v = (rng.normal() * np.sin((a + 1) * X / 4 + rng.normal())
                     + rng.normal() * np.cos((b + 1) * Y / 5) + rng.normal() * Z / 20)
                B[..., a, b] = v; B[..., b, a] = v
        Md[..., 1:, 1:] = B
        u, gap = u_timelike(Md)
        h = h_from_u(u)
        p = lattice_pieces(Md, None, h, cfg)
        fields.append({"E_u": p[0], "E_KT": p[5], "min_eig_gap": float(gap.min()),
                       "u_dev_from_e0": float(np.abs(np.abs(u[..., 0]) - 1).max())})
    out["random_block_diagonal"] = {
        "n_fields": len(fields), "n": n, "L": L,
        "max_abs_E_KT": max(abs(f["E_KT"]) for f in fields),
        "max_E_u": max(f["E_u"] for f in fields),
        "worst_rel": max(abs(f["E_KT"]) / max(f["E_u"], 1e-30) for f in fields),
        "min_eig_gap_seen": min(f["min_eig_gap"] for f in fields),
        "rows": fields[:5]}
    # (iii) a control WITH a time row: must be nonzero
    Md = np.zeros(X.shape + (4, 4))
    Md[..., 0, 0] = G
    Md[..., 1, 1] = 1.0; Md[..., 2, 2] = DELTA
    tr = 0.3 * np.sin(X / 3) * np.exp(-(rr / 12) ** 2)
    Md[..., 0, 1] = tr; Md[..., 1, 0] = tr
    u, _ = u_timelike(Md)
    p = lattice_pieces(Md, None, h_from_u(u), cfg)
    out["control_with_time_row"] = {"E_KT": p[5], "E_u": p[0]}
    # (iv) mutation: the wrong (spatial) u must break the identity
    hw = np.broadcast_to(h_from_u(np.array([0.0, 1.0, 0.0, 0.0])),
                         X.shape + (4, 4)).copy()
    Md2 = np.zeros(X.shape + (4, 4))
    Md2[..., 0, 0] = G + np.cos(rr / 4)
    Md2[..., 1, 1] = 1.0 + 0.2 * np.sin(X / 3)
    Md2[..., 2, 2] = DELTA; Md2[..., 3, 3] = 0.1 * np.cos(Y / 4)
    pm = lattice_pieces(Md2, None, hw, cfg)
    pc = lattice_pieces(Md2, None, h_from_u(u_timelike(Md2)[0]), cfg)
    out["mutation_u_to_e1"] = {"E_KT_mutant": pm[5], "E_KT_correct": pc[5]}
    # (v) the degeneracy locus: M(t) = M_vac + t (e0 e1^T + e1 e0^T), t* = (g+1)/2
    tstar = 0.5 * (G + 1.0)
    Aj = np.stack([rand_sym(np.random.default_rng(5)) for _ in range(4)])
    deg = []
    for f in (0.0, 0.5, 0.9, 0.99, 0.999, 0.9999, 1.0, 1.0001):
        t = f * tstar
        Mt = vac() + t * (np.outer([1, 0, 0, 0], [0, 1, 0, 0])
                          + np.outer([0, 1, 0, 0], [1, 0, 0, 0]))
        try:
            u, gap = u_timelike(Mt)
            deg.append({"t_over_tstar": f, "u0": float(u[0]),
                        "eig_gap": float(gap), "K_T": float(kt_density(Aj, h_from_u(u)))})
        except Exception as e:
            deg.append({"t_over_tstar": f, "error": str(e)})
    out["degeneracy_path"] = {"t_star": tstar, "rows": deg}
    dump("a2_static", out)
    return out


# ============ A2 scope: the campaign's own stored anchors ============
def audit_A2_anchors():
    import glob
    out = {"rows": [], "note": ""}
    files = sorted(glob.glob(os.path.join(DATA, "m5_32_r3_ii", "*.npz")))
    for fn in files:
        M = np.load(fn)["M"].astype(np.float64)
        n = M.shape[0]
        Lb = 48.0 if n == 32 else 72.0
        cfg = B3.base_cfg(s=S, g=G, n=n, L=Lb, delta=DELTA)
        u, gap = u_timelike(M)
        p = lattice_pieces(M, None, h_from_u(u), cfg)
        out["rows"].append({"file": os.path.basename(fn), "n": n,
                            "max_time_row": float(np.abs(M[..., 0, 1:]).max()),
                            "E_u": p[0], "E_V": p[2], "E_KT": p[5],
                            "min_eig_gap": float(gap.min())})
    d = {r["file"]: r for r in out["rows"]}
    # the Coulomb interaction energy carried by K_T: E_KT(pair) - 2 E_KT(single)
    single = d.get("lam0.75_dr1_single_d0_n32.npz", {}).get("E_KT", 0.0)
    ints = {}
    for kind in ("same", "anti"):
        for dd in (10, 14, 18, 24):
            k = f"lam0.75_dr1_{kind}_d{dd}_n32.npz"
            if k in d:
                ints[f"{kind}_d{dd}"] = d[k]["E_KT"] - 2.0 * single
    out["dressed_single_E_KT"] = single
    out["K_T_pair_interaction_per_c2"] = ints
    out["n_with_time_row"] = sum(1 for r in out["rows"] if r["max_time_row"] > 1e-12)
    out["n_total"] = len(out["rows"])
    out["max_E_KT_dressed"] = max(r["E_KT"] for r in out["rows"])
    out["max_E_KT_undressed"] = max(abs(r["E_KT"]) for r in out["rows"]
                                    if r["max_time_row"] <= 1e-12)
    dump("a2_anchors", out)
    return out


# =============================== A3 ===============================
def audit_A3():
    """Legendre fit of K_T(omega) on MY dressed clock, >= 7 omega samples."""
    out = {"rows": []}
    fam = Fam(32, 48.0)
    for (amp, R) in ((0.02, 6.0), (0.025, 12.0), (0.04, 18.0)):
        Qb = qb_field(fam.cfg, fam.b_of(amp, R))
        Md = B3.sym4(conj(Qb, fam.Mb))
        a0d = B3.sym4(conj(Qb, fam.a0))
        u, gap = u_timelike(Md)
        h = h_from_u(u)
        h3 = fam.cfg["h"] ** 3
        # spatial jets (sym stencil, my branch order)
        omegas = np.array([-1.5, -1.0, -0.5, -0.25, 0.0, 0.25, 0.5, 1.0, 1.5])
        vals = []
        for om in omegas:
            tot = 0.0
            for br, wt in BRANCHES:
                A = np.zeros((4,) + Md.shape)
                A[0] = om * a0d
                for ax in range(3):
                    A[ax + 1] = mydiff(Md, ax, fam.cfg["h"], br)
                tot += wt * h3 * np.sum(kt_density(A, h))
            vals.append(float(tot))
        vals = np.array(vals)
        c4 = np.polyfit(omegas, vals, 4)
        c2f = np.polyfit(omegas, vals, 2)
        # direct u-frame reading of the a0 time row
        Lam = boost_to_e0(u)
        a0p = np.einsum("...ab,...bc,...dc->...ad", Lam, a0d, Lam)
        kap_dens = 2.0 * np.sum(a0p[..., 0, 1:] ** 2, axis=-1)
        pieces = lattice_pieces(Md, a0d, h, fam.cfg)
        out["rows"].append({
            "amp": amp, "R": R,
            "quartic_coeff": float(c4[0]), "cubic_coeff": float(c4[1]),
            "omega2_coeff_deg4fit": float(c4[2]), "omega2_coeff_deg2fit": float(c2f[0]),
            "static_K_T": float(c2f[2]), "E_KT_static_from_pieces": pieces[5],
            "kappa_from_pieces": pieces[6],
            "max_a0_timerow_uframe": float(np.abs(a0p[..., 0, 1:]).max()),
            "max_a0_spatial_uframe": float(np.abs(a0p[..., 1:, 1:]).max()),
            "kappa_direct_uframe": float(h3 * np.sum(kap_dens)),
            "tr_h_a0_h_a0_scale": float(np.abs(tr_XAXA(h, a0d)).max()),
            "tr_eta_a0_eta_a0_scale": float(np.abs(tr_XAXA(ETA, a0d)).max()),
            "min_eig_gap": float(gap.min())})
        # MUTATION: the wrong u must give a nonzero omega^2 coefficient
        hw = h_from_u(np.broadcast_to(np.array([0.0, 1.0, 0.0, 0.0]),
                                      Md.shape[:-1]).copy())
        kap_mut = h3 * 0.5 * np.sum(tr_XAXA(hw, a0d) - tr_XAXA(ETA, a0d))
        out["rows"][-1]["kappa_mutant_u_to_e1"] = float(kap_mut)
    out["worst_abs_kappa"] = max(abs(r["kappa_from_pieces"]) for r in out["rows"])
    out["worst_rel_kappa_to_E_KT"] = max(
        abs(r["kappa_from_pieces"]) / abs(r["E_KT_static_from_pieces"])
        for r in out["rows"])
    out["worst_abs_quartic"] = max(abs(r["quartic_coeff"]) for r in out["rows"])
    dump("a3_clock_omega2", out)
    return out


# ============================= A4 / A5 =============================
def sym_jets30(x):
    """30-dim vector -> 3 symmetric spatial jets."""
    A = np.zeros((3, 4, 4))
    k = 0
    for i in range(3):
        for a in range(4):
            for b in range(a, 4):
                A[i, a, b] = x[k]; A[i, b, a] = x[k]; k += 1
    return A


def K0_of(a0, M):
    """the omega^2 ENERGY contribution of -c2 K_T, per unit c2:
    K0 = 2 sum_j (a0^uframe)_0j (a0^uframe)_j0."""
    u, _ = u_timelike(M)
    Lam = boost_to_e0(u)
    ap = Lam @ a0 @ Lam.T
    return float(2.0 * np.sum(ap[0, 1:] * ap[1:, 0]))


def H2_form(a0, M, lam, c2):
    """the 30x30 omega^2 form of the ENERGY plus the c2 constant."""
    u, _ = u_timelike(M)
    h = h_from_u(u)
    Q = np.zeros((30, 30))
    for p in range(30):
        for q in range(p, 30):
            def val(x):
                A = sym_jets30(x)
                t = 0.0
                for i in range(3):
                    F = comm_eta(a0, A[i])
                    t += 4.0 * ((1 - lam) * q_eta(F) + lam * q_h(F, h))
                return t
            ep = np.zeros(30); ep[p] = 1.0
            eq = np.zeros(30); eq[q] = 1.0
            v = 0.5 * (val(ep + eq) - val(ep) - val(eq))
            Q[p, q] = Q[q, p] = v
    return Q, K0_of(a0, M)


def audit_A4_A5(rng):
    out = {"channels": {}, "lambda_star": {}, "missing_channel": {}}
    M = vac()
    ch = {}
    for k in (1, 2, 3):
        ch[f"boost_{k}"] = (tangent(gen_boost(k), M), M)
    for (i, j) in ((1, 2), (2, 3), (1, 3)):
        ch[f"rot_{i}{j}"] = (tangent(gen_rot(i, j), M), M)
    # a boost tangent on a Lorentz-dressed background (u != e0)
    Ld = lorentz(rap=(0.6, -0.3, 0.2), ang=(0.4, 0.1, -0.2))
    Md = Ld @ M @ Ld.T
    ch["boost_1_dressed_bg"] = (tangent(gen_boost(1), Md), Md)
    # MY antisymmetric probes (NOT tangents of a symmetric field)
    ch["antisym_probe"] = (gen_boost(1) @ M - M @ gen_boost(1), M)
    # THE MISSING CHANNEL: the actual hedgehog clock flow at a lattice point
    fam = Fam(16, 48.0)
    idx = (11, 5, 9)
    a0_hh = fam.a0[idx]; M_hh = fam.Mb[idx]
    ch["hedgehog_clock_REAL"] = (a0_hh, M_hh)
    for nm, (a0, Mb) in ch.items():
        K0 = K0_of(a0, Mb)
        pred = None
        if nm.startswith("boost_") and nm[6:].isdigit():
            k = int(nm[6:])
            pred = 2.0 * (G + M[k, k]) ** 2
        sym_dev = float(np.abs(a0 - a0.T).max()) / max(np.abs(a0).max(), 1e-30)
        # jet-independence of the omega^2 coefficient
        devs = []
        for _ in range(6):
            x = rng.normal(size=30) * rng.choice([0.1, 1.0, 10.0])
            A = np.zeros((4, 4, 4)); A[1:] = sym_jets30(x)
            u, _ = u_timelike(Mb); h = h_from_u(u)
            def kt(w):
                A2 = A.copy(); A2[0] = w * a0
                return float(kt_density(A2, h))
            C = 0.5 * (kt(1.0) + kt(-1.0)) - kt(0.0)
            devs.append(abs(-C - K0))
        out["channels"][nm] = {
            "K0_per_c2": K0, "K0_predicted_2_g_plus_Mjj_sq": pred,
            "K0_rel_dev_vs_prediction": (None if pred is None else rel(K0, pred)),
            "a0_symmetric": bool(sym_dev < 1e-12),
            "a0_antisymmetry_rel": sym_dev,
            "jet_independence_worst_abs": float(max(devs))}
    # A5: the LP / crossing, MY way
    lams = np.array([0.0, 0.25, 0.4, 0.45, 0.49, 0.5, 0.51, 0.6, 0.75, 1.0])
    for nm in ("boost_1", "rot_12", "boost_1_dressed_bg", "hedgehog_clock_REAL"):
        a0, Mb = ch[nm]
        rows = {}
        for lam in lams:
            Q, K0 = H2_form(a0, Mb, lam, 0.0)
            w = np.linalg.eigvalsh(0.5 * (Q + Q.T))
            rows[f"{lam:g}"] = {"min_eig_H2_jetform": float(w.min()), "K0_per_c2": K0}
        # for each c2, the smallest lambda at which the TOTAL omega^2 energy
        # stays >= 0 for EVERY jet scale
        lstar = {}
        for c2 in (0.0, 0.1, 1.0, 10.0, 100.0, 1000.0):
            ok = [lam for lam in lams
                  if rows[f"{lam:g}"]["min_eig_H2_jetform"] >= -1e-6]
            lstar[f"c2_{c2:g}"] = (min(ok) if ok else None)
        # the explicit jet scale at which a given c2 stops helping
        K0 = rows["0"]["K0_per_c2"]
        fail = {}
        for c2 in (0.1, 1.0, 10.0, 100.0):
            me = rows["0"]["min_eig_H2_jetform"]
            fail[f"c2_{c2:g}"] = (None if me >= 0 or K0 <= 0
                                  else float(np.sqrt(c2 * K0 / abs(me))))
        out["lambda_star"][nm] = {"per_lambda": rows, "lambda_star_by_c2": lstar,
                                  "jet_scale_where_c2_fails_at_lam0": fail}
    out["missing_channel"] = {
        "name": "hedgehog_clock_REAL",
        "why": ("the a0 actually used by every lattice run (B8.a0_unit on the "
                "m5_21_8 hedgehog) is symmetric AND block diagonal, so its K0 "
                "is exactly 0: c2 gives ZERO help on the only clock channel "
                "the model realizes; the producer's 31-channel list contains "
                "no channel built from it"),
        "K0_per_c2": out["channels"]["hedgehog_clock_REAL"]["K0_per_c2"]}
    dump("a4_a5_channels", out)
    return out


# =============================== A6 ===============================
def audit_A6():
    out = {}
    fam = Fam(64, 96.0)
    amp = 0.02
    # (i) TRUE dilation: b_mu(r) = b(r/mu). The hedgehog background is exactly
    # scale free (it depends on the angles only), so this IS a full dilation of
    # the configuration and the continuum exponent must be exactly +1.
    for R in (6.0, 12.0):
        mus = np.array([0.6, 0.8, 1.0, 1.25, 1.6])
        e = [fam.pieces(amp, R, mu=m)[5] for m in mus]
        sl = np.polyfit(np.log(mus), np.log(e), 1)[0]
        pair = [float(np.log(e[k + 1] / e[k]) / np.log(mus[k + 1] / mus[k]))
                for k in range(len(mus) - 1)]
        out[f"dilation_R{R:g}"] = {"mu": mus.tolist(), "E_KT": e,
                                   "loglog_exponent": float(sl),
                                   "pair_exponents": pair}
    # (ii) the R ladder is NOT a dilation (the tanh(r/2) core is not rescaled)
    Rs = [3.0, 4.5, 6.0, 9.0, 12.0, 18.0, 24.0, 36.0, 48.0]
    e = [fam.pieces(amp, R)[5] for R in Rs]
    pair = [float(np.log(e[k + 1] / e[k]) / np.log(Rs[k + 1] / Rs[k]))
            for k in range(len(Rs) - 1)]
    out["R_ladder"] = {"R": Rs, "E_KT": e, "pair_exponents": pair,
                       "loglog_R6_to_48": float(np.polyfit(
                           np.log(Rs[2:]), np.log(e[2:]), 1)[0]),
                       "loglog_all": float(np.polyfit(np.log(Rs), np.log(e), 1)[0])}
    # (iii) resolution check: same L, finer h
    fine = Fam(64, 48.0)
    coarse = Fam(32, 48.0)
    out["resolution"] = {}
    for R in (6.0, 12.0):
        out["resolution"][f"R{R:g}"] = {
            "E_KT_h1.5": coarse.pieces(amp, R)[5],
            "E_KT_h0.75": fine.pieces(amp, R)[5]}
    mus = np.array([0.8, 1.0, 1.25])
    ef = [fine.pieces(amp, 6.0, mu=m)[5] for m in mus]
    out["resolution"]["dilation_exponent_h0.75_R6"] = float(
        np.polyfit(np.log(mus), np.log(ef), 1)[0])
    # (iv) positivity + the sign mutation
    out["positivity"] = {"all_positive": all(v > 0 for v in e),
                         "min_E_KT": min(e),
                         "sign_mutant_c2_negative_gives": [-v for v in e[:3]]}
    # (v) amp scaling (must be exactly 2 at small amp)
    amps = np.array([0.005, 0.01, 0.02, 0.04])
    ea = [fam.pieces(a, 12.0)[5] for a in amps]
    out["amp_exponent_R12"] = float(np.polyfit(np.log(amps), np.log(ea), 1)[0])
    dump("a6_derrick", out)
    return out


# ============================= B1..B6 =============================
# MY amp grid (finer and different from the producer's), MY R grid (with
# deliberate OFF-grid points to expose R* grid locking), MY minimizer.
MY_AMPS = np.array(sorted(set(np.round(np.concatenate([
    np.linspace(0.0, 0.032, 17), np.linspace(0.036, 0.06, 7),
    [0.07, 0.085, 0.1, 0.13, 0.17, 0.22, 0.3]]), 6))))


def my_R_grid(L):
    base = [3.0, 3.8, 4.5, 5.2, 6.0, 7.0, 8.0, 9.0, 10.5, 12.0,
            13.5, 15.0, 16.5, 18.0, 21.0]
    return [r for r in base if r < L / 2.0] + [L / 2.0]


def _loc_parab(xs, ys, i):
    """local quadratic through (i-1, i, i+1); returns a callable."""
    c = np.polyfit(xs[i - 1:i + 2], ys[i - 1:i + 2], 2)
    return lambda a: np.polyval(c, a)


def min_over_amp(rows, J, lam, c2):
    """MY minimizer: dense amp grid, then SEPARATE local quadratic models of
    E_stat(amp) and kin(amp) around the grid argmin, then golden section on
    E_J = E_stat + J^2/(4 kin) over that bracket. This differs structurally
    from a 3-point parabola fitted to E_J itself."""
    amps = np.array([a for a, _ in rows])
    es = np.array([assemble(p, lam, c2)[0] for _, p in rows])
    kk = np.array([assemble(p, lam, c2)[1] for _, p in rows])
    ok = kk > 0
    E = np.where(ok, es + J * J / (4.0 * np.where(ok, kk, 1.0)), np.inf)
    i = int(np.argmin(E))
    a_star, E_star, k_star = float(amps[i]), float(E[i]), float(kk[i])
    refined = False
    if 0 < i < len(amps) - 1 and ok[i - 1] and ok[i + 1]:
        fe = _loc_parab(amps, es, i)
        fk = _loc_parab(amps, kk, i)
        f = lambda a: fe(a) + J * J / (4.0 * max(fk(a), 1e-12))
        lo, hi = amps[i - 1], amps[i + 1]
        gr = 0.5 * (np.sqrt(5.0) - 1.0)
        x1, x2 = hi - gr * (hi - lo), lo + gr * (hi - lo)
        f1, f2 = f(x1), f(x2)
        for _ in range(60):
            if f1 < f2:
                hi, x2, f2 = x2, x1, f1
                x1 = hi - gr * (hi - lo); f1 = f(x1)
            else:
                lo, x1, f1 = x1, x2, f2
                x2 = lo + gr * (hi - lo); f2 = f(x2)
        a_star = float(0.5 * (lo + hi)); E_star = float(f(a_star))
        k_star = float(fk(a_star)); refined = True
    return {"amp_star": a_star, "E_J": E_star,
            "E_stat": float(E_star - J * J / (4.0 * k_star)), "kin": k_star,
            "omega_star": float(J / (2.0 * k_star)), "grid_index": i,
            "parabolic_refined": refined,
            "amp_edge": "interior" if 0 < i < len(amps) - 1
                        else ("amp_min" if i == 0 else "amp_max")}


def audit_B(n, L, lams=(0.75, 1.0), c2s=(0.0, 0.01, 0.03, 0.1, 0.3, -0.1),
            J=200.0, amps=None):
    fam = Fam(n, L)
    amps = MY_AMPS if amps is None else amps
    Rs = my_R_grid(L)
    log(f"B: box n={n} L={L}: {len(Rs)} R x {len(amps)} amps")
    grid = {}
    for R in Rs:
        grid[R] = [(float(a), fam.pieces(float(a), R)) for a in amps]
        log(f"  R={R:g} done")
    out = {"n": n, "L": L, "h": fam.cfg["h"], "J": J, "R_grid": Rs,
           "n_amps": len(amps), "amp_max": float(amps[-1]),
           "undressed": {}, "rows": {}}
    p0 = fam.pieces(0.0, Rs[0])
    for lam in lams:
        es, kk = assemble(p0, lam, 0.0)
        out["undressed"][f"lam_{lam:g}"] = {
            "E_stat": es, "kin": kk, "E_J": es + J * J / (4 * kk),
            "omega_star": J / (2 * kk), "E_KT": p0[5], "kappa": p0[6]}
    for lam in lams:
        for c2 in c2s:
            byR = {R: min_over_amp(grid[R], J, lam, c2) for R in Rs}
            EJ = np.array([byR[R]["E_J"] for R in Rs])
            i = int(np.argmin(EJ))
            und = out["undressed"][f"lam_{lam:g}"]["E_J"]
            best = byR[Rs[i]]
            out["rows"][f"lam_{lam:g}_c2_{c2:g}"] = {
                "R_star_grid": Rs[i], "R_index": i,
                "at_wall": i == len(Rs) - 1, "at_Rmin": i == 0,
                "interior": 0 < i < len(Rs) - 1,
                "amp_star": best["amp_star"], "omega_star": best["omega_star"],
                "kin": best["kin"], "E_J": best["E_J"], "E_stat": best["E_stat"],
                "E_J_undressed": und,
                "dressing_gain": best["E_J"] - und,
                "dressing_gain_frac": (best["E_J"] - und) / abs(und),
                "amp_edge": best["amp_edge"],
                "E_J_of_R": EJ.tolist(),
                "amp_of_R": [byR[R]["amp_star"] for R in Rs]}
    dump(f"b_box_n{n}_L{int(L)}", out)
    return out


# =============================== B7 ===============================
def kin_and_shells(n, L, delta=DELTA, stencil="sym", nshell=None):
    """undressed hedgehog: (E_u, E_V, kin) and the RADIAL SHELL profile of
    the kinetic density; also returns the per-cell kinetic density."""
    cfg = B3.base_cfg(s=S, g=G, n=n, L=L, delta=delta)
    Mb = B8.dressed(cfg, 0.0)
    a0 = B8.a0_unit(cfg, 0.0)
    h3 = cfg["h"] ** 3
    brs = BRANCHES if stencil == "sym" else ((stencil, 1.0),)
    dens = np.zeros((n, n, n))
    eu = 0.0
    for br, wt in brs:
        A = [mydiff(Mb, ax, cfg["h"], br) for ax in range(3)]
        for i in range(3):
            for j in range(i + 1, 3):
                eu += wt * 4.0 * np.sum(q_eta(comm_eta(A[i], A[j])))
            dens += wt * 4.0 * q_eta(comm_eta(a0, A[i]))
    dens *= h3
    _, ev = B3.e_parts(Mb, cfg)
    X, Y, Z = B3.coords(n, cfg["h"])
    r = np.sqrt(X * X + Y * Y + Z * Z)
    return {"n": n, "L": L, "h": cfg["h"], "delta": delta, "stencil": stencil,
            "E_u": h3 * eu, "E_V": float(ev), "kin": float(dens.sum()),
            "a0_norm_min": float(np.linalg.norm(a0, axis=(-2, -1)).min()),
            "a0_norm_max": float(np.linalg.norm(a0, axis=(-2, -1)).max())}, dens, r


def audit_B7(ns=(32, 48, 64, 80), h=1.5):
    out = {"boxes": [], "shells": {}, "stencils": {}, "delta_ladder": [],
           "cell_overlap": {}, "drift": {}}
    store = {}
    for n in ns:
        L = n * h
        t = time.time()
        row, dens, r = kin_and_shells(n, L)
        row["wall_s"] = round(time.time() - t, 2)
        out["boxes"].append(row)
        store[n] = (dens, r, L)
        log(f"B7 box n={n} L={L}: kin {row['kin']:.4f} E_u {row['E_u']:.4f} "
            f"({row['wall_s']}s)")
        dump("b7_partial", out)
    K = [b["kin"] for b in out["boxes"]]
    Ls = [b["L"] for b in out["boxes"]]
    out["linearity"] = {
        "L": Ls, "kin": K,
        "first_differences": [K[i + 1] - K[i] for i in range(len(K) - 1)],
        "slope_per_unit_L": [(K[i + 1] - K[i]) / (Ls[i + 1] - Ls[i])
                             for i in range(len(K) - 1)],
        "second_differences": [K[i + 2] - 2 * K[i + 1] + K[i]
                               for i in range(len(K) - 2)],
        "affine_fit": np.polyfit(Ls, K, 1).tolist(),
        "loglog_exponent_pairs": [
            float(np.log(K[i + 1] / K[i]) / np.log(Ls[i + 1] / Ls[i]))
            for i in range(len(K) - 1)],
        "loglog_exponent_all": float(np.polyfit(np.log(Ls), np.log(K), 1)[0]),
        "E_u": [b["E_u"] for b in out["boxes"]],
        "E_u_first_differences": [out["boxes"][i + 1]["E_u"] - out["boxes"][i]["E_u"]
                                  for i in range(len(K) - 1)]}
    # residual of the affine fit (a pure power law with exponent 1.07 would not fit)
    a, b = out["linearity"]["affine_fit"]
    out["linearity"]["affine_residual_rel"] = [
        float(abs(a * Ls[i] + b - K[i]) / K[i]) for i in range(len(K))]
    p = np.polyfit(np.log(Ls), np.log(K), 1)
    out["linearity"]["powerlaw_residual_rel"] = [
        float(abs(np.exp(np.polyval(p, np.log(Ls[i]))) - K[i]) / K[i])
        for i in range(len(K))]
    # shell decomposition in the biggest box
    nb = ns[-1]
    dens, r, L = store[nb]
    edges = np.arange(0.0, L / 2.0 * np.sqrt(3.0) + 1.5, 1.5)
    idx = np.digitize(r.ravel(), edges) - 1
    tot = np.bincount(idx, weights=dens.ravel(), minlength=len(edges))
    cnt = np.bincount(idx, minlength=len(edges))
    mid = 0.5 * (edges[:-1] + edges[1:])
    per_dr = tot[:len(mid)] / np.diff(edges)
    out["shells"] = {
        "n": nb, "L": L, "edges": edges.tolist(),
        "r_mid": mid.tolist(),
        "kin_per_shell": tot[:len(mid)].tolist(),
        "kin_per_unit_r": per_dr.tolist(),
        "cells_per_shell": cnt[:len(mid)].tolist(),
        "cumulative_fraction": (np.cumsum(tot[:len(mid)]) / tot.sum()).tolist(),
        "note": ("a 1/r^2 density gives kin_per_unit_r = constant while the "
                 "shell is fully inside the box; the fall past r = L/2 is the "
                 "cube corner geometry, not the field")}
    # is it core dominated? fraction inside r < 6 and r < 12
    frac = {}
    for rc in (3.0, 6.0, 12.0, 24.0, 36.0):
        frac[f"r_lt_{rc:g}"] = float(dens[r < rc].sum() / dens.sum())
    out["shells"]["fraction_inside"] = frac
    # the flat plateau level of kin_per_unit_r (core excluded, box wall excluded)
    m = (mid > 6.0) & (mid < L / 2.0 * 0.85)
    out["shells"]["plateau_kin_per_unit_r_mean"] = float(per_dr[m].mean())
    out["shells"]["plateau_kin_per_unit_r_std"] = float(per_dr[m].std())
    # BOUNDARY ARTIFACT test: cell by cell against the smallest box
    dsm, rsm, Lsm = store[ns[0]]
    o = (nb - ns[0]) // 2
    sub = dens[o:o + ns[0], o:o + ns[0], o:o + ns[0]]
    d = np.abs(sub - dsm)
    interior = np.zeros_like(d, dtype=bool)
    interior[1:-1, 1:-1, 1:-1] = True
    out["cell_overlap"] = {
        "small_box_n": ns[0], "big_box_n": nb,
        "max_abs_diff_all_cells": float(d.max()),
        "max_abs_diff_interior_cells": float(d[interior].max()),
        "sum_small": float(dsm.sum()), "sum_big_on_same_cells": float(sub.sum()),
        "rel_dev_sum_interior": float(
            abs(sub[interior].sum() - dsm[interior].sum())
            / abs(dsm[interior].sum())),
        "boundary_shell_sum_small": float(dsm[~interior].sum()),
        "boundary_shell_sum_big": float(sub[~interior].sum()),
        "note": ("h is identical and the grids align, so a cell-by-cell match "
                 "away from the small box face proves the L growth is pure "
                 "added volume, not a boundary or stencil artifact")}
    # stencil variants on the middle box
    nm = ns[len(ns) // 2]
    for st in ("sym", "fwd", "bwd"):
        r1, _, _ = kin_and_shells(ns[0], ns[0] * h, stencil=st)
        r2, _, _ = kin_and_shells(nm, nm * h, stencil=st)
        out["stencils"][st] = {
            "kin_small": r1["kin"], "kin_mid": r2["kin"],
            "slope_per_unit_L": (r2["kin"] - r1["kin"]) / (r2["L"] - r1["L"])}
    # delta ladder: a0 = delta * (direction field), so the extensive slope
    # should vanish at delta = 0 and grow like delta^2
    for dl in (0.0, 0.05, 0.1, 0.2, 0.3, 0.5):
        r1, _, _ = kin_and_shells(ns[0], ns[0] * h, delta=dl)
        r2, _, _ = kin_and_shells(nm, nm * h, delta=dl)
        out["delta_ladder"].append({
            "delta": dl, "kin_small": r1["kin"], "kin_mid": r2["kin"],
            "a0_norm": r1["a0_norm_max"],
            "slope_per_unit_L": (r2["kin"] - r1["kin"]) / (r2["L"] - r1["L"])})
    dl = [row for row in out["delta_ladder"] if row["delta"] > 0]
    out["delta_scaling_exponent"] = float(np.polyfit(
        np.log([row["delta"] for row in dl]),
        np.log([row["slope_per_unit_L"] for row in dl]), 1)[0])
    # omega* ladder at J = 200 and the drift
    J = 200.0
    om = [J / (2.0 * b["kin"]) for b in out["boxes"]]
    out["drift"] = {
        "J": J, "L": Ls, "omega_star": om,
        "drift_pairs": [om[i + 1] / om[i] - 1.0 for i in range(len(om) - 1)],
        "pure_1_over_L_prediction": [Ls[i] / Ls[i + 1] - 1.0
                                     for i in range(len(om) - 1)],
        "omega_star_times_L": [om[i] * Ls[i] for i in range(len(om))]}
    dump("b7_boxes", out)
    return out


# =============================== main ===============================
def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    rng = np.random.default_rng(20260828)
    if mode in ("a1", "form", "all"):
        audit_A1(rng)
    if mode in ("a2", "form", "all"):
        audit_A2(rng); audit_A2_anchors()
    if mode in ("a3", "form", "all"):
        audit_A3()
    if mode in ("a45", "form", "all"):
        audit_A4_A5(rng)
    if mode in ("a6", "form", "all"):
        audit_A6()
    if mode in ("lat", "all"):
        audit_B(32, 48.0)
        audit_B(48, 72.0)
    if mode in ("b7", "all"):
        audit_B7()
    if mode in ("merge", "all"):
        _load("merge", "m5_32_r7_audit_merge.py").main()
    log(f"mode {mode} done")


if __name__ == "__main__":
    main()
