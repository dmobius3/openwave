"""M5.32 R3 construction (iii): the fixed-J two-clock CROSS-INERTIA C(d) of
two boost-dressed electrons under the covariant lambda-family, and the
Newton sign per ensemble (fixed J vs fixed omega).

EQUATIONS FIRST (the two-clock Legendre algebra, our conventions)
-----------------------------------------------------------------
Field M(x) real symmetric 4x4 per cell, eta = diag(-1, 1, 1, 1), spatial
jets A_i = d_i M on the certified sym stencil, time jet A_0 = dM/dt.
Every registered term is exactly quadratic in A_0 (h = h(M) only), so the
kinetic energy of a configuration whose time derivative is a superposition
of two clock flows,
    A_0 = omega_1 a0_1 + omega_2 a0_2,
is the quadratic form (the FACTOR-4 BRIDGE of the registry, per lambda)
    T(omega_1, omega_2) = K(A_0) = 4 h^3 sum_br wt sum_cells sum_i
                          q_lambda(F_0i),   F_0i = A_0 eta A_i - A_i eta A_0,
    q_lambda(F, G)     = (1 - lambda) tr(eta F eta G^T) + lambda tr(h F h G^T)
so with the symmetric bilinear K(a, b) (K(a, a) = the record's kin_of at
lambda = 0)
    T = K11 omega_1^2 + K22 omega_2^2 + 2 K12 omega_1 omega_2,
    K_kl = K(a0_k, a0_l).
Convention used here (stated because the external formula's factors differ):
    T = 1/2 [I_1 omega_1^2 + I_2 omega_2^2 + 2 C omega_1 omega_2],
    I_k = 2 K_kk,  C = 2 K12                                         (1)
so that the single electron has T = 1/2 I_0 omega^2 = kin omega^2 and
J = dT/domega = I_0 omega = 2 kin omega (the record's J).  Conjugate momenta
    J_k = dT/domega_k,   J = I omega,   I = [[I_1, C], [C, I_2]]        (2)
Legendre transform at fixed (J_1, J_2) (the physical ensemble, each electron
carries its own J):
    E_J(d) = E_stat(d) + 1/2 J^T I(d)^-1 J                             (3)
    like clocks J_1 = J_2 = j:
    E_J = E_stat + j^2 (I_1 + I_2 - 2 C) / [2 (I_1 I_2 - C^2)]
        = E_stat + j^2 / (I_0 + C)          when I_1 = I_2 = I_0        (4)
    Delta E_J = E_J - E_J|_{C = 0} = - j^2 C / [I_0 (I_0 + C)]          (5)
Fixed omega (omega_1 = omega_2 = omega) instead:
    E_omega = E_stat + omega^2 (I_0 + C)                                (6)
so the clock-mediated interaction has OPPOSITE signs in the two ensembles:
at fixed J a cross-inertia C > 0 that grows as d shrinks LOWERS E_J toward
small d (ATTRACTION, force = -dE/dd < 0 means toward each other: E_J rises
with d); at fixed omega the same C > 0 RAISES E_omega toward small d
(REPULSION).  The external formula Delta E_J = -j^2 C / [2 I_0 (I_0 + C)]
is (5) in the convention T = I_0 omega^2 (no 1/2; I = kin, J = 2 I omega,
E = J^2 / (4 I)); the SIGN statement is convention-free:
    sign(Delta E_J) = -sign(C),  sign(Delta E_omega) = +sign(C).
The fixed-J energy is always read from the 2x2 inversion (3), never from
the raw fixed-omega value (import I5); both are reported.  Anti-aligned
clocks (omega_2 -> -omega_2) map C -> -C in (4) by the bilinearity of K:
a null test of the algebra, not a physical statement (control f).

THE CLOCK GENERATORS (the record's a0, localized per core)
----------------------------------------------------------
The record's unit-omega clock flow of the twisting hedgehog (m5_21_8_b
a0_unit = FD in t of Qh(t) d Qh(t)^T) is exactly the commutator with the
rotation generator about the LOCAL RADIAL direction n(x) = x / |x|:
    a0 = [G_n, M] = G_n M - M G_n,  G_n = local_rot(n)   (gate: <= 1e-8)
For two cores at x_k = (0, 0, +-d/2), n_k(x) = (x - x_k) / |x - x_k|,
    c_k = [G_{n_k}, M0]  (M0 the undressed composite pair, 4x4 embedded),
    a0_k = w_k(x) Qb c_k Qb^T     (the boost dressing outside, as the record)
with a scalar WINDOW w_k(x) localizing clock k to its own core. Windows
(all reported; the kinetic density is LOCAL in x, so K12 is supported only
where w_1 w_2 != 0):
    full   w_k = 1                        each clock stirs the whole field
    pou    w_1 = 1/2 (1 + tanh(z / s_w)), w_2 = 1 - w_1, s_w = 2.0
                                          the midplane partition of unity
    loc    w_k = exp(-(r_k / R_w)^2), R_w = 6.0   core-localized clocks
    hard   w_1 = [z > 0], w_2 = [z <= 0]  the Voronoi split: C == 0 EXACTLY
                                          (locality null, control)

THE PAIR FIELD
--------------
The certified two-center 3x3 composition of m5_21_4_a_pair (seed_pair
'same' = inverse-stereographic product with the +z escape tube; 'anti' =
the mirror-hedgehog sum) rebuilt WITHOUT the isotropic core blend (r_c = 4
in the record) so that the single-core limit is EXACTLY the record's
twisting hedgehog B8.dressed(cfg, 0) (gate: max dev 0; the blended seed
lowers kin by 14 % and fails the 1 % I_0 gate, measured); M = n n^T + delta
ph ph^T on the spatial block, M_00 = -s g = 32; each core boost-dressed on
the rigid family bl_k = amp tanh(r_k / 2) along n_k at the fixed-J optimum
amp* = opt_amp_refined of R2 arm b (m5_32_r2_bounded.json, n = 32, L = 48)
per (lambda, J), Qb = Qb_1 Qb_2 (order checked at d = 12).  NO relaxation:
this is the UN-RELAXED cross-inertia (the pair is not a stationary point of
L_lambda); the optional `relax` stage descends the d = 12 pair 300 steps
under L_lambda (R2 arm b energy-backtracking FIRE, pinned shell) and
re-reads C with a0_k = w_k [G_{n_k}, M_relaxed].

STAGES (python3 m5_32_r3_iii_crossinertia.py STAGE [--n --L --ds --lam])
    gate      a0 identity, single-limit identity, record kin_base, the
              polarization identity, locality + vacuum nulls, boost order
    scan      the (kind, d, amp, window) table of K11, K22, K12 (eta and h
              parts) + the static reads
    single    the single electron at each amp in the same box
    sens      window-width sensitivity of C (s_w 1/2/4, R_w 4/6/9)
    relax     d = 12 pair relaxed 300 steps under L_lambda, C re-read
    collect   fits, E_J / E_omega tables, force signs, plot, JSON
Partials: ../checkpoints/m5_32_r3_iii/*.json (gitignored)
Out: ../data/m5_32_r3_crossinertia.json, ../plots/m5_32_r3_crossinertia.png
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
CKPT = os.path.join(HERE, "..", "checkpoints", "m5_32_r3_iii")


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


R2B = _load("m5_32_r2_b_bounded", "m5_32_r2_b_bounded.py")
PAIR = _load("m5_21_4_a_pair", "m5_21_4_a_pair.py")
B3, B8 = R2B.B3, R2B.B8
ETA = B3.ETA
LAMBDAS = (0.0, 0.75, 1.0)
J_LADDER = (50.0, 200.0, 800.0)
DS = (8.0, 10.0, 12.0, 14.0, 18.0, 24.0)
WINDOWS = ("full", "pou", "loc", "hard")
S_W, R_W = 2.0, 6.0
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:8.1f}s] {msg}", flush=True)


def ck(name):
    os.makedirs(CKPT, exist_ok=True)
    return os.path.join(CKPT, name)


def dump(name, obj):
    with open(ck(name), "w") as f:
        json.dump(obj, f, indent=1)


# ================= amp* per (lambda, J) from R2 arm b =================
def amps_from_r2():
    with open(os.path.join(DATA, "m5_32_r2_bounded.json")) as f:
        fx = json.load(f)["fixedj"]["fixedj_n32"]
    out = {"lam_0": {f"J_{J:g}": 0.0 for J in J_LADDER}}
    for lam in (0.75, 1.0):
        tab = fx["byLambda"][f"lam_{lam:g}"]["table"]
        out[f"lam_{lam:g}"] = {f"J_{J:g}": float(tab[f"J_{J:g}"]["opt_amp_refined"])
                              for J in J_LADDER}
    return out


def amp_set():
    a = amps_from_r2()
    vals = sorted({round(v, 9) for d in a.values() for v in d.values()})
    return vals, a


# ================= geometry + generators =================
def local_rot(v):
    """G_n: rotation generator about the spatial unit vector n (the
    m5_21_3 gen_catalog convention; G_{e_x} = G1 of m5_21_8)."""
    W = np.zeros(v.shape[:-1] + (4, 4))
    n1, n2, n3 = v[..., 0], v[..., 1], v[..., 2]
    W[..., 1, 2], W[..., 1, 3] = -n3, n2
    W[..., 2, 1], W[..., 2, 3] = n3, -n1
    W[..., 3, 1], W[..., 3, 2] = -n2, n1
    return W


def core_geom(cfg, zc):
    """(r_k, n_k, K, K2) of the core at (0, 0, zc)."""
    X, Y, Z = B3.coords(cfg["n"], cfg["h"])
    Zc = Z - zc
    R = np.sqrt(X * X + Y * Y + Zc * Zc)
    nx, ny, nz = X / R, Y / R, Zc / R
    nh = np.stack([nx, ny, nz], -1)
    K = np.zeros(X.shape + (4, 4))
    K[..., 0, 1], K[..., 0, 2], K[..., 0, 3] = nx, ny, nz
    K[..., 1, 0], K[..., 2, 0], K[..., 3, 0] = nx, ny, nz
    K2 = np.zeros_like(K)
    K2[..., 0, 0] = 1.0
    for i, a in enumerate((nx, ny, nz)):
        for j, bb in enumerate((nx, ny, nz)):
            K2[..., 1 + i, 1 + j] = a * bb
    return R, nh, K, K2


def qb_of(bl, K, K2):
    return (np.eye(4)[None, None, None] + np.sinh(bl)[..., None, None] * K
            + (np.cosh(bl) - 1.0)[..., None, None] * K2)


def conj(Q, M):
    return B3.sym4(np.einsum("...ab,...bc,...dc->...ad", Q, M, Q))


# ================= the composite pair (unblended) =================
def composite_nhat(cfg, kind, d):
    """the director of m5_21_4_a_pair.seed_pair (verbatim formulas) for
    'same' / 'anti' / 'single'."""
    n, h = cfg["n"], cfg["h"]
    X, Y, Z = B3.coords(n, h)
    rho = np.sqrt(X * X + Y * Y)
    zt, zb = +d / 2.0, -d / 2.0
    tht = np.arctan2(rho, Z - zt)
    thb = np.arctan2(rho, Z - zb)
    if kind == "single":
        return PAIR._nhat_from_alpha(n, h, np.arctan2(rho, Z))
    if kind == "anti":
        return PAIR._nhat_from_alpha(n, h, tht + np.pi - thb)
    if kind == "same":
        tt = np.clip(np.tan(0.5 * tht), 0.0, 1e6)
        tb = np.clip(np.tan(0.5 * thb), 0.0, 1e6)
        umag = np.clip(tt * tb, 0.0, 1e12)
        rhos = np.where(rho < 1e-12, 1e-12, rho)
        c2p = (X * X - Y * Y) / (rhos * rhos)
        s2p = -(2.0 * X * Y) / (rhos * rhos)
        den = 1.0 + umag * umag
        nhat = np.stack([2.0 * umag * c2p / den, 2.0 * umag * s2p / den,
                         (1.0 - umag * umag) / den], axis=-1)
        r_t = 3.0
        zed = 0.5 * (np.tanh((Z - zb) / 2.0) - np.tanh((Z - zt) / 2.0))
        et = np.exp(-((rho / r_t) ** 2)) * zed
        nhat = nhat + et[..., None] * np.array([0.0, 0.0, 1.0])
        return nhat / np.maximum(np.linalg.norm(nhat, axis=-1)[..., None], 1e-300)
    raise ValueError(kind)


def tensor_unblended(cfg, nhat):
    """M3 = n n^T + delta ph ph^T (m5_21_4 _tensor_from_nhat with the core
    blend REMOVED: w = 1)."""
    n, h, delta = cfg["n"], cfg["h"], cfg["delta"]
    X, Y, Z = B3.coords(n, h)
    rho = np.sqrt(X * X + Y * Y)
    rhos = np.where(rho < 1e-12, 1e-12, rho)
    phihat = np.stack([-Y / rhos, X / rhos, np.zeros_like(Z)], axis=-1)
    dot = np.einsum("...a,...a->...", phihat, nhat)[..., None]
    ph = phihat - dot * nhat
    ph = ph / np.maximum(np.linalg.norm(ph, axis=-1)[..., None], 1e-300)
    return (nhat[..., :, None] * nhat[..., None, :]
            + delta * ph[..., :, None] * ph[..., None, :])


def pair_field(cfg, kind, d, amp, order="12", blended=False):
    """(M dressed, [a0_1, a0_2] base commutators dressed (unwindowed),
    geometry dict)."""
    if blended:
        M3 = PAIR.seed_pair(cfg, kind, d)
    else:
        M3 = tensor_unblended(cfg, composite_nhat(cfg, kind, d))
    M0 = B3.embed34(M3, cfg)
    centers = [+d / 2.0, -d / 2.0] if kind != "single" else [0.0]
    geo = [core_geom(cfg, zc) for zc in centers]
    cs = [local_rot(g[1]) @ M0 - M0 @ local_rot(g[1]) for g in geo]
    Qb = None
    seq = geo if order == "12" else geo[::-1]
    for R, nh, K, K2 in seq:
        Q = qb_of(amp * np.tanh(R / 2.0), K, K2)
        Qb = Q if Qb is None else Qb @ Q
    M = conj(Qb, M0)
    a0s = [conj(Qb, c) for c in cs]
    return M, a0s, geo, M0


def windows(cfg, geo):
    X, Y, Z = B3.coords(cfg["n"], cfg["h"])
    w1 = 0.5 * (1.0 + np.tanh(Z / S_W))
    out = {"full": (np.ones_like(Z), np.ones_like(Z)),
           "pou": (w1, 1.0 - w1),
           "loc": tuple(np.exp(-(g[0] / R_W) ** 2) for g in geo),
           "hard": ((Z > 0).astype(float), (Z <= 0).astype(float))}
    return out


# ================= the kinetic bilinear =================
class Field:
    """one eigen pass + one derivative pass per field; then any number of
    kinetic bilinears K(a, b) = (eta part, h part), kin = (1-lam) eta +
    lam h; plus the static reads (E_stat_lam = 4 A_lam + V4)."""

    def __init__(self, M, cfg):
        self.cfg, self.h3 = cfg, cfg["h"] ** 3
        u0, _, _, _, _, ok, gap = R2B.tl_eig(M)
        self.ok, self.gap = bool(np.all(ok)), float(np.min(gap))
        self.hh = R2B.h_of(u0)
        self.A = [(wt, [B3.d1(M, ax, cfg["h"], br) for ax in range(3)])
                  for br, wt in B3.branches(cfg["stencil"])]
        acc = np.zeros(2)
        for wt, A in self.A:
            for i in range(3):
                for j in range(i + 1, 3):
                    F = B3.comm_eta(A[i], A[j])
                    acc[0] += wt * np.sum(R2B.q_eta(F))
                    acc[1] += wt * np.sum(R2B.q_h(F, self.hh))
        _, ev = B3.e_parts(M, cfg)
        self.static = {"A_I1": float(self.h3 * acc[0]),
                       "A_I1h": float(self.h3 * acc[1]), "V4": float(ev),
                       "max_abs_M0i": float(np.max(np.abs(M[..., 0, 1:]))),
                       "min_gap": self.gap, "eigvec_ok": self.ok}

    def bilinear(self, a, b):
        acc = np.zeros(2)
        for wt, A in self.A:
            for i in range(3):
                Fa = B3.comm_eta(a, A[i])
                Fb = Fa if b is a else B3.comm_eta(b, A[i])
                acc[0] += wt * np.sum(B3.inner_eta(Fa, Fb))
                acc[1] += wt * np.sum(np.einsum(
                    "...ab,...bc,...cd,...ad->...", self.hh, Fa, self.hh, Fb,
                    optimize=True))
        return {"eta": float(4.0 * self.h3 * acc[0]),
                "h": float(4.0 * self.h3 * acc[1])}

    def e_stat(self, lam):
        s = self.static
        return 4.0 * ((1 - lam) * s["A_I1"] + lam * s["A_I1h"]) + s["V4"]


def kin_lam(kk, lam):
    return (1 - lam) * kk["eta"] + lam * kk["h"]


def kin_table(fld, a0s, win):
    """K11, K22, K12 per window (eta / h parts)."""
    out = {}
    for wn, (w1, w2) in win.items():
        a1 = w1[..., None, None] * a0s[0]
        a2 = w2[..., None, None] * a0s[1]
        out[wn] = {"K11": fld.bilinear(a1, a1), "K22": fld.bilinear(a2, a2),
                   "K12": fld.bilinear(a1, a2)}
    return out


# ================= GATE =================
def stage_gate(n=32, L=48.0):
    cfg = R2B.cfg_of(n, L)
    out = {"n": n, "L": L, "h": cfg["h"], "checks": {}}
    chk = out["checks"]
    # (1) the record's a0 == [G_n, M] on the undressed hedgehog
    Mb, a0b = B8.dressed(cfg, 0.0), B8.a0_unit(cfg, 0.0)
    R, nh, _, _ = core_geom(cfg, 0.0)
    a0c = local_rot(nh) @ Mb - Mb @ local_rot(nh)
    chk["a0_identity_rel"] = float(np.max(np.abs(a0c - a0b)) / np.max(np.abs(a0b)))
    # (2) the unblended single composite == B8.dressed(cfg, 0) exactly
    Ms, a0s, _, _ = pair_field(cfg, "single", 0.0, 0.0)
    chk["single_unblended_vs_record_max_dev"] = float(np.max(np.abs(Ms - Mb)))
    Msb, _, _, _ = pair_field(cfg, "single", 0.0, 0.0, blended=True)
    chk["single_blended_vs_record_max_dev"] = float(np.max(np.abs(Msb - Mb)))
    # (3) record kin_base + E_stat_base (control c) at lambda = 0
    with open(os.path.join(DATA, "m5_32_r2_bounded.json")) as f:
        fx = json.load(f)["fixedj"]["fixedj_n32"]
    fld = Field(Ms, cfg)
    kk = fld.bilinear(a0s[0], a0s[0])
    chk["kin_base_record"] = fx["kin_base"]["lam_0"]
    chk["kin_base_here_lam0"] = kk["eta"]
    # the record's a0 is a t-finite-difference (rel 6.7e-9 vs the exact
    # commutator), so the kin_base bar is 1e-6, not machine precision
    chk["kin_base_rel"] = float(abs(kk["eta"] / fx["kin_base"]["lam_0"] - 1.0))
    chk["kin_base_vs_B3_kin_of_rel"] = float(
        abs(kk["eta"] / B3.kin_of(Ms, a0s[0], cfg) - 1.0))
    chk["E_stat_base_record"] = fx["E_stat_base"]["lam_0"]
    chk["E_stat_base_here_lam0"] = fld.e_stat(0.0)
    chk["E_stat_base_rel"] = float(abs(fld.e_stat(0.0) / fx["E_stat_base"]["lam_0"] - 1.0))
    fldb = Field(Msb, cfg)
    kb = fldb.bilinear(local_rot(nh) @ Msb - Msb @ local_rot(nh),
                       local_rot(nh) @ Msb - Msb @ local_rot(nh))
    chk["kin_blended_seed_rel_to_record"] = float(kb["eta"] / fx["kin_base"]["lam_0"] - 1.0)
    # (3b) the dressed single at amp*(lam 1, J 200) vs the R2 fixed-J read
    amp1 = fx["byLambda"]["lam_1"]["table"]["J_200"]["opt_amp_refined"]
    Md, a0d, _, _ = pair_field(cfg, "single", 0.0, amp1)
    fd = Field(Md, cfg)
    kd = fd.bilinear(a0d[0], a0d[0])
    chk["dressed_kin_lam1_here"] = kd["h"]
    chk["dressed_kin_lam1_record"] = fx["byLambda"]["lam_1"]["table"]["J_200"]["kin_opt"]
    chk["dressed_kin_lam1_rel"] = float(abs(kd["h"] / chk["dressed_kin_lam1_record"] - 1.0))
    chk["dressed_Estat_lam1_here"] = fd.e_stat(1.0)
    chk["dressed_Estat_lam1_record"] = fx["byLambda"]["lam_1"]["table"]["J_200"]["E_stat_opt"]
    chk["dressed_Estat_lam1_rel"] = float(abs(fd.e_stat(1.0) / chk["dressed_Estat_lam1_record"] - 1.0))
    # (4) polarization identity vs the direct bilinear, at d = 12, amp1
    M, a0s, geo, _ = pair_field(cfg, "same", 12.0, amp1)
    fp = Field(M, cfg)
    win = windows(cfg, geo)
    pol = {}
    for wn in ("full", "pou"):
        w1, w2 = win[wn]
        a1, a2 = w1[..., None, None] * a0s[0], w2[..., None, None] * a0s[1]
        kp, km = fp.bilinear(a1 + a2, a1 + a2), fp.bilinear(a1 - a2, a1 - a2)
        k12 = fp.bilinear(a1, a2)
        pol[wn] = {"K12_direct": k12,
                   "K12_polarization": {k: (kp[k] - km[k]) / 4.0 for k in k12},
                   "rel": {k: float(abs((kp[k] - km[k]) / 4.0 - k12[k])
                                    / max(abs(k12[k]), 1e-300)) for k in k12},
                   "anti_aligned_K12": fp.bilinear(a1, -a2)}
    chk["polarization"] = pol
    # (5) locality null: hard window K12 == 0
    w1, w2 = win["hard"]
    chk["hard_window_K12"] = fp.bilinear(w1[..., None, None] * a0s[0],
                                         w2[..., None, None] * a0s[1])
    # (6) vacuum null: cores replaced by the vacuum, same two generators
    Mv = np.broadcast_to(B3.vac4(cfg), M.shape).copy()
    cv = [local_rot(g[1]) @ Mv - Mv @ local_rot(g[1]) for g in geo]
    fv = Field(Mv, cfg)
    chk["vacuum_null"] = {"K11": fv.bilinear(cv[0], cv[0]),
                          "K12": fv.bilinear(cv[0], cv[1]),
                          "a0_norm2": float(np.sum(cv[0] * cv[0]))}
    # (7) boost composition order at d = 12
    M21, a21, _, _ = pair_field(cfg, "same", 12.0, amp1, order="21")
    f21 = Field(M21, cfg)
    k12a = kin_table(fp, a0s, win)["full"]["K12"]
    k12b = kin_table(f21, a21, win)["full"]["K12"]
    chk["boost_order"] = {"K12_full_12": k12a, "K12_full_21": k12b,
                          "M_max_dev": float(np.max(np.abs(M - M21))),
                          "rel_h": float(abs(k12b["h"] / k12a["h"] - 1.0))}
    # (8) generator definition sensitivity: [G, M_dressed] vs Qb [G, M0] Qb^T
    alt = [local_rot(g[1]) @ M - M @ local_rot(g[1]) for g in geo]
    chk["generator_alt_commutator_on_dressed"] = {
        "K12_full_alt": fp.bilinear(alt[0], alt[1]),
        "K12_full_record": k12a}
    out["gate_pass"] = bool(chk["a0_identity_rel"] <= 1e-6
                            and chk["single_unblended_vs_record_max_dev"] <= 1e-12
                            and chk["kin_base_rel"] <= 1e-6
                            and chk["E_stat_base_rel"] <= 1e-10
                            and chk["dressed_kin_lam1_rel"] <= 0.01
                            and all(v <= 1e-9 for p in pol.values() for v in p["rel"].values())
                            and abs(chk["hard_window_K12"]["h"]) <= 1e-12
                            and abs(chk["vacuum_null"]["K12"]["h"]) <= 1e-12)
    dump("gate.json", out)
    log("GATE " + json.dumps({k: v for k, v in chk.items()
                              if not isinstance(v, dict)}, default=float))
    log(f"GATE polarization {json.dumps(pol, default=float)}")
    log(f"GATE hard {chk['hard_window_K12']} vac {chk['vacuum_null']} "
        f"order {chk['boost_order']} alt {chk['generator_alt_commutator_on_dressed']}")
    log(f"GATE PASS = {out['gate_pass']}")
    return out


# ================= SCAN =================
def stage_scan(n, L, ds, kinds=("same", "anti")):
    cfg = R2B.cfg_of(n, L)
    amps, amap = amp_set()
    out = {"n": n, "L": L, "h": cfg["h"], "amps": amps, "amp_map": amap,
           "windows": {"full": "w = 1", "pou": f"1/2(1 + tanh(z/{S_W}))",
                       "loc": f"exp(-(r_k/{R_W})^2)", "hard": "[z > 0]"},
           "rows": []}
    tag = f"n{n}_L{L:g}"
    for kind in kinds:
        for d in ds:
            for amp in amps:
                t0 = time.time()
                M, a0s, geo, _ = pair_field(cfg, kind, d, amp)
                fld = Field(M, cfg)
                row = {"kind": kind, "d": d, "amp": amp, "static": fld.static,
                       "K": kin_table(fld, a0s, windows(cfg, geo)),
                       "wall_s": round(time.time() - t0, 1)}
                out["rows"].append(row)
                kf, kp = row["K"]["full"], row["K"]["pou"]
                log(f"SCAN {tag} {kind} d {d:g} amp {amp:g}: E0 {fld.e_stat(0):.3f} "
                    f"E1 {fld.e_stat(1):.3f} | full K11 {kf['K11']['h']:.2f} "
                    f"K12 eta {kf['K12']['eta']:+.4f} h {kf['K12']['h']:+.4f} | "
                    f"pou K12 eta {kp['K12']['eta']:+.4f} h {kp['K12']['h']:+.4f} | "
                    f"loc K12 h {row['K']['loc']['K12']['h']:+.4f} "
                    f"hard {row['K']['hard']['K12']['h']:+.1e}")
                dump(f"scan_{tag}.json", out)
    return out


def stage_single(n, L):
    cfg = R2B.cfg_of(n, L)
    amps, amap = amp_set()
    out = {"n": n, "L": L, "h": cfg["h"], "rows": []}
    for amp in amps:
        M, a0s, geo, _ = pair_field(cfg, "single", 0.0, amp)
        fld = Field(M, cfg)
        row = {"amp": amp, "static": fld.static,
               "K": fld.bilinear(a0s[0], a0s[0])}
        out["rows"].append(row)
        log(f"SINGLE n{n} L{L:g} amp {amp:g}: E0 {fld.e_stat(0):.3f} E1 "
            f"{fld.e_stat(1):.3f} K eta {row['K']['eta']:.3f} h {row['K']['h']:.3f}")
    dump(f"single_n{n}_L{L:g}.json", out)
    return out


# ================= SENS: window-width sensitivity =================
def stage_sens(n=32, L=48.0, ds=(12.0, 24.0), lam=1.0, J=200.0):
    """C at lambda for the pou window s_w in {1, 2, 4} and the loc window
    R_w in {4, 6, 9}, like and opposite charge (a sensitivity report, the
    pre-registered windows stay s_w = 2, R_w = 6)."""
    cfg = R2B.cfg_of(n, L)
    _, amap = amp_set()
    amp = amap[f"lam_{lam:g}"][f"J_{J:g}"]
    X, Y, Z = B3.coords(n, cfg["h"])
    out = {"n": n, "L": L, "h": cfg["h"], "lam": lam, "J": J, "amp": amp, "rows": []}
    for kind in ("same", "anti"):
        for d in ds:
            M, a0s, geo, _ = pair_field(cfg, kind, d, amp)
            fld = Field(M, cfg)
            row = {"kind": kind, "d": d, "pou": {}, "loc": {}}
            for sw in (1.0, 2.0, 4.0):
                w1 = 0.5 * (1.0 + np.tanh(Z / sw))
                k = fld.bilinear(w1[..., None, None] * a0s[0], (1 - w1)[..., None, None] * a0s[1])
                k11 = fld.bilinear(w1[..., None, None] * a0s[0], w1[..., None, None] * a0s[0])
                row["pou"][f"s_w_{sw:g}"] = {"C": 2 * kin_lam(k, lam), "I1": 2 * kin_lam(k11, lam)}
            for rw in (4.0, 6.0, 9.0):
                w = [np.exp(-(g[0] / rw) ** 2) for g in geo]
                k = fld.bilinear(w[0][..., None, None] * a0s[0], w[1][..., None, None] * a0s[1])
                k11 = fld.bilinear(w[0][..., None, None] * a0s[0], w[0][..., None, None] * a0s[0])
                row["loc"][f"R_w_{rw:g}"] = {"C": 2 * kin_lam(k, lam), "I1": 2 * kin_lam(k11, lam)}
            out["rows"].append(row)
            log(f"SENS {kind} d {d:g}: " + json.dumps(row["pou"]) + " " + json.dumps(row["loc"]))
    dump(f"sens_n{n}_L{L:g}.json", out)
    return out


# ================= RELAX (optional) =================
def stage_relax(lam, n=32, L=48.0, d=12.0, steps=300, J=200.0):
    cfg = R2B.cfg_of(n, L)
    _, amap = amp_set()
    amp = amap[f"lam_{lam:g}"][f"J_{J:g}"]
    M0, a0s0, geo, _ = pair_field(cfg, "same", d, amp)
    win = windows(cfg, geo)
    free = (~B3.pin_shell(n, cfg["h"]))[..., None, None].astype(float)
    M = M0.copy()
    E, G, info = R2B.energy_grad(M, cfg, lam)
    v = np.zeros_like(M)
    dt, alpha, n_up, dt_min = 0.02, 0.1, 0, 1e-7
    F = -G * free
    E_prev, E0 = E, E
    n_acc = n_rej = 0
    stop = "budget"
    for it in range(1, steps + 1):
        P = float(np.sum(F * v))
        if P > 0.0:
            n_up += 1
            vn, fn = np.sqrt(np.sum(v * v)), np.sqrt(np.sum(F * F))
            v = (1 - alpha) * v + alpha * (F / max(fn, 1e-300)) * vn
            if n_up > 5:
                dt = min(dt * 1.1, 0.2)
                alpha *= 0.99
        else:
            v[:] = 0.0
            alpha, n_up = 0.1, 0
        v_try = v + dt * F
        M_try = M + dt * v_try
        E, G, info = R2B.energy_grad(M_try, cfg, lam)
        reject = (not info["ok"]) or not np.isfinite(E) or E > E_prev + 1e-12 * max(abs(E_prev), 1.0)
        if reject:
            n_rej += 1
            dt *= 0.5
            v[:] = 0.0
            alpha, n_up = 0.1, 0
            if dt < dt_min:
                stop = "STALLED"
                break
            continue
        n_acc += 1
        M, v, E_prev = M_try, v_try, E
        F = -G * free
        if it % 50 == 0:
            log(f"RELAX lam {lam:g} d {d:g} it {it} E {E:.4f} dt {dt:.2e} acc {n_acc} rej {n_rej}")
    fld = Field(M, cfg)
    a0_rec = a0s0                                   # the record generators
    a0_rel = [local_rot(g[1]) @ M - M @ local_rot(g[1]) for g in geo]
    out = {"lam": lam, "n": n, "L": L, "h": cfg["h"], "d": d, "amp": amp,
           "J": J, "steps": steps, "stop": stop, "E0": E0, "E_end": E_prev,
           "accepted": n_acc, "rejected": n_rej, "dt_final": dt,
           "rel_move": float(np.sqrt(np.sum((M - M0) ** 2)) / np.sqrt(np.sum((M0 - B3.vac4(cfg)) ** 2))),
           "static": fld.static,
           "K_generators_record": kin_table(fld, a0_rec, win),
           "K_generators_relaxed_commutator": kin_table(fld, a0_rel, win)}
    dump(f"relax_lam{lam:g}_d{d:g}_n{n}.json", out)
    log(f"RELAX lam {lam:g} done: E {E0:.4f} -> {E_prev:.4f} move {out['rel_move']:.3e} "
        f"K12 full h rec {out['K_generators_record']['full']['K12']['h']:+.4f} "
        f"rel {out['K_generators_relaxed_commutator']['full']['K12']['h']:+.4f} "
        f"pou h rec {out['K_generators_record']['pou']['K12']['h']:+.4f} "
        f"rel {out['K_generators_relaxed_commutator']['pou']['K12']['h']:+.4f}")
    return out


# ================= COLLECT =================
def fits(ds, C):
    ds, C = np.asarray(ds, float), np.asarray(C, float)
    out = {}
    for name, cols in (("A/d + B/d^2", [1 / ds, 1 / ds ** 2]),
                       ("A ln(d)/d + B/d", [np.log(ds) / ds, 1 / ds])):
        X = np.stack(cols, 1)
        coef, *_ = np.linalg.lstsq(X, C, rcond=None)
        res = C - X @ coef
        ss = float(np.sum((C - C.mean()) ** 2))
        out[name] = {"A": float(coef[0]), "B": float(coef[1]),
                     "rms_res": float(np.sqrt(np.mean(res ** 2))),
                     "R2": float(1 - np.sum(res ** 2) / ss) if ss > 0 else None}
    if np.all(C > 0) or np.all(C < 0):
        p = np.polyfit(np.log(ds), np.log(np.abs(C)), 1)
        out["power_law"] = {"exponent": float(p[0]), "prefactor": float(np.sign(C[0]) * np.exp(p[1]))}
    else:
        out["power_law"] = {"exponent": None, "note": "sign change in the window"}
    return out


def ensemble_rows(rows_by_d, single_K, lam, J, amp, wn):
    """per d: I1, I2, C, E_stat, E_J (2x2 inversion), E_omega (matching
    omega = J / I0_single), the clock parts."""
    I0 = 2.0 * kin_lam(single_K, lam)
    om = J / I0
    tab = []
    for d, row in rows_by_d:
        K = row["K"][wn]
        I1, I2, C = (2.0 * kin_lam(K[k], lam) for k in ("K11", "K22", "K12"))
        s = row["static"]
        Es = 4.0 * ((1 - lam) * s["A_I1"] + lam * s["A_I1h"]) + s["V4"]
        det = I1 * I2 - C * C
        EJ_clock = J * J * (I1 + I2 - 2 * C) / (2 * det) if det > 0 else float("nan")
        EJ_clock_noC = J * J * (I1 + I2) / (2 * I1 * I2)
        Eom_clock = om * om * (I1 + I2 + 2 * C) / 2.0
        tab.append({"d": d, "I1": I1, "I2": I2, "C": C, "C_over_I0": C / I0,
                    "det_pos": bool(det > 0), "E_stat": Es,
                    "E_J_clock": EJ_clock, "E_J_clock_noC": EJ_clock_noC,
                    "dE_J_from_C": EJ_clock - EJ_clock_noC,
                    "E_J_total": Es + EJ_clock,
                    "E_omega_clock": Eom_clock,
                    "dE_omega_from_C": om * om * C,
                    "E_omega_total": Es + Eom_clock,
                    "E_J_raw_fixed_omega_read_WRONG": Es + om * om * (I1 + I2 + 2 * C) / 2.0})
    ref = {"I0_single": I0, "omega_match": om,
           "E_stat_2single": None, "E_J_clock_2single": J * J / I0,
           "E_omega_clock_2single": om * om * I0}
    return tab, ref


def force_sign(tab, key, dmin=12.0):
    """sign of dE/dd over the outer window by finite differences; ATTRACTIVE
    iff E rises with d (dE/dd > 0) at every step."""
    pts = [(r["d"], r[key]) for r in tab if r["d"] >= dmin]
    sl = [(pts[i + 1][1] - pts[i][1]) / (pts[i + 1][0] - pts[i][0]) for i in range(len(pts) - 1)]
    if all(s > 0 for s in sl):
        v = "ATTRACTIVE (E rises with d)"
    elif all(s < 0 for s in sl):
        v = "REPULSIVE (E falls with d)"
    else:
        v = "NON-MONOTONE"
    return {"slopes_dE_dd": sl, "verdict": v, "window_d": [p[0] for p in pts]}


def stage_collect():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out = {"task": "M5.32 R3 construction (iii): fixed-J two-clock cross-inertia",
           "candidate": "L_lambda = -4 [(1 - lambda) I1 + lambda I1_h] - V4",
           "conventions": {"T": "1/2 [I1 w1^2 + I2 w2^2 + 2 C w1 w2], I_k = 2 K_kk, C = 2 K12",
                           "J": "J_k = dT/dw_k = (I w)_k; single: J = I0 w = 2 kin w",
                           "E_J": "E_stat + 1/2 J^T I^-1 J; like clocks j: E_stat + j^2 (I1+I2-2C)/(2(I1 I2 - C^2))",
                           "E_omega": "E_stat + w^2 (I1 + I2 + 2C)/2, w = j / I0_single",
                           "sign": "Delta E_J = -j^2 C/[I0 (I0+C)] (external: same with an extra 1/2 from T = I0 w^2)",
                           "attractive": "E(d) rises with d (force -dE/dd points toward the partner)"},
           "point": {"g": R2B.G_MAIN, "s": R2B.S_MAIN, "delta": R2B.DELTA, "stencil": "sym"},
           "limitation": "UN-RELAXED cross-inertia: the composite pair is the certified "
                         "seed, boost-dressed, not a stationary point of L_lambda; the "
                         "relax stage (if present) gives the 300-step descended read at d = 12"}

    def rd(name):
        p = ck(name)
        return json.load(open(p)) if os.path.exists(p) else None
    out["gate"] = rd("gate.json")
    scans = {os.path.basename(p)[5:-5]: json.load(open(p)) for p in sorted(glob.glob(ck("scan_*.json")))}
    singles = {os.path.basename(p)[7:-5]: json.load(open(p)) for p in sorted(glob.glob(ck("single_*.json")))}
    out["relax"] = {os.path.basename(p)[6:-5]: json.load(open(p)) for p in sorted(glob.glob(ck("relax_*.json")))}
    out["window_sensitivity"] = {os.path.basename(p)[5:-5]: json.load(open(p)) for p in sorted(glob.glob(ck("sens_*.json")))}
    out["boxes"] = {}
    for tag, sc in scans.items():
        sg = singles.get(tag)
        box = {"n": sc["n"], "L": sc["L"], "h": sc["h"], "windows": sc["windows"],
               "C_table": {}, "ensembles": {}, "controls": {}}
        amap = sc["amp_map"]
        single_by_amp = {round(r["amp"], 9): r for r in sg["rows"]} if sg else {}

        def rows_of(kind, amp):
            return sorted([(r["d"], r) for r in sc["rows"]
                           if r["kind"] == kind and abs(r["amp"] - amp) < 1e-9])
        # C(d) tables per lambda, window, kind at each lambda's amp*(J = 200)
        for lam in LAMBDAS:
            for J in J_LADDER:
                amp = amap[f"lam_{lam:g}"][f"J_{J:g}"]
                sK = single_by_amp.get(round(amp, 9), {}).get("K")
                for kind in ("same", "anti"):
                    rws = rows_of(kind, amp)
                    if not rws or sK is None:
                        continue
                    I0 = 2.0 * kin_lam(sK, lam)
                    for wn in WINDOWS:
                        ds = [d for d, _ in rws]
                        C = [2.0 * kin_lam(r["K"][wn]["K12"], lam) for _, r in rws]
                        I1 = [2.0 * kin_lam(r["K"][wn]["K11"], lam) for _, r in rws]
                        key = f"lam_{lam:g}/J_{J:g}/{kind}/{wn}"
                        ent = {"amp": amp, "d": ds, "C": C, "I1": I1,
                               "I0_single": I0, "C_over_I0": [c / I0 for c in C],
                               "sign_C": [int(np.sign(c)) for c in C],
                               "I1_over_I0_minus_1": [i / I0 - 1 for i in I1]}
                        if len(ds) == 2 and wn != "hard" and C[0] * C[1] > 0:
                            ent["two_point_exponent"] = float(np.log(abs(C[1]) / abs(C[0])) / np.log(ds[1] / ds[0]))
                        if len(ds) >= 3 and wn != "hard":
                            ent["fit_all"] = fits(ds, C)
                            outer = [(d, c) for d, c in zip(ds, C) if d >= 12]
                            if len(outer) >= 3:
                                ent["fit_outer_d_ge_12"] = fits([o[0] for o in outer], [o[1] for o in outer])
                        box["C_table"][key] = ent
                    if kind == "same":
                        for wn in ("full", "pou", "loc"):
                            tab, ref = ensemble_rows(rws, sK, lam, J, amp, wn)
                            ref["E_stat_2single"] = 2 * (4.0 * ((1 - lam) * single_by_amp[round(amp, 9)]["static"]["A_I1"]
                                                                + lam * single_by_amp[round(amp, 9)]["static"]["A_I1h"])
                                                         + single_by_amp[round(amp, 9)]["static"]["V4"])
                            box["ensembles"][f"lam_{lam:g}/J_{J:g}/{wn}"] = {
                                "rows": tab, "reference_2single": ref,
                                "force_dE_J_from_C_only": force_sign(tab, "dE_J_from_C"),
                                "force_dE_omega_from_C_only": force_sign(tab, "dE_omega_from_C"),
                                "force_E_J_clock": force_sign(tab, "E_J_clock"),
                                "force_E_omega_clock": force_sign(tab, "E_omega_clock"),
                                "force_E_stat": force_sign(tab, "E_stat"),
                                "force_E_J_total": force_sign(tab, "E_J_total"),
                                "force_E_omega_total": force_sign(tab, "E_omega_total")}
        # controls: (b) I0 recovery: single K vs record; (d) mutation; (e) anti
        if sg:
            rec = out["gate"]["checks"] if out["gate"] else {}
            box["controls"]["b_I0_single"] = {
                "kin_amp0_lam0": single_by_amp[0.0]["K"]["eta"],
                "record_kin_base_n32_L48": rec.get("kin_base_record"),
                "note": "the n32/L48 box reproduces the record; other boxes carry their own I0 (kin is IR-extensive, ~L)"}
        box["controls"]["d_mutation_lam1_vs_lam0"] = {
            k.replace("lam_1/", "MUT/"): {
                "C_lam1": v["C"],
                "C_lam0_same_amp": [2.0 * kin_lam(r["K"][k.split("/")[-1]]["K12"], 0.0)
                                    for _, r in rows_of(k.split("/")[2], v["amp"])],
                "sign_flip": [int(np.sign(a)) != int(np.sign(b)) for a, b in zip(
                    v["C"], [2.0 * kin_lam(r["K"][k.split("/")[-1]]["K12"], 0.0)
                             for _, r in rows_of(k.split("/")[2], v["amp"])])]}
            for k, v in box["C_table"].items() if k.startswith("lam_1/J_200/") and not k.endswith("hard")}
        box["controls"]["e_anti_vs_same"] = {
            wn: {"C_same": box["C_table"].get(f"lam_1/J_200/same/{wn}", {}).get("C"),
                 "C_anti": box["C_table"].get(f"lam_1/J_200/anti/{wn}", {}).get("C"),
                 "C_same_lam0": box["C_table"].get(f"lam_0/J_200/same/{wn}", {}).get("C"),
                 "C_anti_lam0": box["C_table"].get(f"lam_0/J_200/anti/{wn}", {}).get("C")}
            for wn in WINDOWS}
        out["boxes"][tag] = box
    # ---------- plot ----------
    os.makedirs(PLOTS, exist_ok=True)
    fig, ax = plt.subplots(2, 2, figsize=(13, 9))
    for tag, box in out["boxes"].items():
        for lam, ls in zip(LAMBDAS, ("--", "-.", "-")):
            for wn, mk in (("full", "o"), ("pou", "s"), ("loc", "^")):
                e = box["C_table"].get(f"lam_{lam:g}/J_200/same/{wn}")
                if e:
                    ax[0, 0].plot(e["d"], e["C_over_I0"], ls, marker=mk, ms=4,
                                  label=f"{tag} lam {lam:g} {wn}")
        for kind, c in (("same", "C0"), ("anti", "C3")):
            e = box["C_table"].get(f"lam_1/J_200/{kind}/full")
            if e:
                ax[0, 1].plot(e["d"], e["C"], "o-", color=c, label=f"{tag} {kind} full lam 1")
            e = box["C_table"].get(f"lam_1/J_200/{kind}/pou")
            if e:
                ax[0, 1].plot(e["d"], e["C"], "s--", color=c, label=f"{tag} {kind} pou lam 1")
        for wn, mk in (("full", "o"), ("pou", "s")):
            en = box["ensembles"].get(f"lam_1/J_200/{wn}")
            if en:
                dd = [r["d"] for r in en["rows"]]
                ej = np.array([r["E_J_clock"] for r in en["rows"]])
                eo = np.array([r["E_omega_clock"] for r in en["rows"]])
                ax[1, 0].plot(dd, ej - ej[-1], "-", marker=mk, label=f"{tag} fixed-J {wn}")
                ax[1, 0].plot(dd, eo - eo[-1], ":", marker=mk, label=f"{tag} fixed-omega {wn}")
                es = np.array([r["E_stat"] for r in en["rows"]])
                ax[1, 1].plot(dd, es - es[-1], "-", marker=mk, label=f"{tag} E_stat {wn}")
                et = np.array([r["E_J_total"] for r in en["rows"]])
                ax[1, 1].plot(dd, et - et[-1], "--", marker=mk, label=f"{tag} E_J total {wn}")
    ax[0, 0].axhline(0, color="k", lw=0.5); ax[0, 0].set_xlabel("d"); ax[0, 0].set_ylabel("C / I0")
    ax[0, 0].set_title("cross-inertia C(d) / I0, like charges, J = 200 amp*"); ax[0, 0].legend(fontsize=5)
    ax[0, 1].axhline(0, color="k", lw=0.5); ax[0, 1].set_xlabel("d"); ax[0, 1].set_ylabel("C")
    ax[0, 1].set_title("like vs opposite charge (Coulomb control), lambda = 1"); ax[0, 1].legend(fontsize=6)
    ax[1, 0].axhline(0, color="k", lw=0.5); ax[1, 0].set_xlabel("d"); ax[1, 0].set_ylabel("clock energy - value at d_max")
    ax[1, 0].set_title("clock-mediated energy per ensemble, lambda = 1, J = 200"); ax[1, 0].legend(fontsize=6)
    ax[1, 1].axhline(0, color="k", lw=0.5); ax[1, 1].set_xlabel("d"); ax[1, 1].set_ylabel("E - E(d_max)")
    ax[1, 1].set_title("static (un-relaxed) vs total fixed-J energy"); ax[1, 1].legend(fontsize=6)
    fig.suptitle("M5.32 R3 (iii): fixed-J two-clock cross-inertia on the lambda-family")
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS, "m5_32_r3_crossinertia.png"), dpi=110)
    out["collected_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    with open(os.path.join(DATA, "m5_32_r3_crossinertia.json"), "w") as f:
        json.dump(out, f, indent=1)
    log("COLLECT written")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["gate", "scan", "single", "sens", "relax", "collect"])
    ap.add_argument("--n", type=int, default=32)
    ap.add_argument("--L", type=float, default=48.0)
    ap.add_argument("--ds", type=str, default=",".join(f"{d:g}" for d in DS))
    ap.add_argument("--lam", type=float, default=1.0)
    ap.add_argument("--steps", type=int, default=300)
    a = ap.parse_args()
    ds = tuple(float(x) for x in a.ds.split(","))
    if a.stage == "gate":
        ok = stage_gate(a.n, a.L)["gate_pass"]
        sys.exit(0 if ok else 1)
    elif a.stage == "scan":
        stage_scan(a.n, a.L, ds)
    elif a.stage == "single":
        stage_single(a.n, a.L)
    elif a.stage == "sens":
        stage_sens(a.n, a.L)
    elif a.stage == "relax":
        stage_relax(a.lam, a.n, a.L, 12.0, a.steps)
    else:
        stage_collect()
    log(f"done {a.stage}")


if __name__ == "__main__":
    main()
