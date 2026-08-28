"""M5.32 R2 arm (b): lattice BOUNDEDNESS (G5) of the covariant lambda-family
and the clock preview (the fixed-J electron, the concavity certificate).

EQUATIONS FIRST
---------------
Field M(x): real symmetric 4x4 per cell; eta = diag(-1, 1, 1, 1); jets
A_i = d_i M on the certified sym stencil (1/2 (fwd + bwd), density per
branch); A_0 = omega a0 (a0 the clock flow of the configuration).
    F_{mu nu} = A_mu eta A_nu - A_nu eta A_mu
    <F, F>_eta = tr(eta F eta F^T)                       (the I1 quadratic form)
    u(M): the timelike unit eigenvector of N = M eta, u^T eta u = -1
    h(M) = eta + 2 (eta u)(eta u)^T                       (the covariant flip metric)
    <F, F>_h = tr(h F h F^T)                              (the I1_h quadratic form)
The lambda-family (lambda = 0 the certified action, lambda = 1 the covariant
flip I1_h; both terms are exactly quadratic in omega because h depends on M
only):
    L_lambda = -4 [(1 - lambda) I1 + lambda I1_h] - V4
    E_lambda(M; omega) = E_stat_lambda(M) + omega^2 kin_lambda(M)
    E_stat_lambda = 4 h^3 sum_br wt sum_cells sum_{i<j} q_lambda(F_ij) + V4
    kin_lambda    = 4 h^3 sum_br wt sum_cells sum_i     q_lambda(F_0i / omega)
    q_lambda(F)   = (1 - lambda) <F, F>_eta + lambda <F, F>_h
(the FACTOR-4 BRIDGE of m5_32_lagrangian.py; at lambda = 0 these are
exactly B3.e_total and B3.kin_of, gated to 1e-12 in the `gate` stage; at
lambda = 1 the omega^2 coefficient equals the registry's I1_h read, gated).
V4 = w sum_p (tr((M eta)^p) - C_p)^2 is conjugation-invariant, so V4 = 0
on every member of the boost-dressing families (they are Lorentz orbits of
the vacuum per cell).

Gradient of E_stat_lambda wrt M (the new instrument; the I1 part is the
certified B3.grad, the I1_h part has two pieces):
  (a) through F at fixed h: dq_lambda/dF = 2 W, W = (1-lambda) eta F eta
      + lambda h F h;  dF_ij = dA_i eta A_j + A_i eta dA_j - dA_j eta A_i
      - A_j eta dA_i, chained with the exact stencil adjoints d1_adj;
  (b) through h(M) at fixed F: with S = sum_{i<j} F_ij h F_ij^T (symmetric),
      d<F,F>_h = 2 tr(S dh),  dh = 2 [(eta du)(eta u)^T + (eta u)(eta du)^T]
      so d q_h = v . du with v = 8 eta S eta u, and the first-order
      eigenvector perturbation of N = M eta (m5_32_terms_ext du_np)
          du = sum_{k != 0} sigma_k (u_k^T eta dM eta u_0) / (lambda_0 - lambda_k) u_k
      gives  d q_h / dM = sym( (eta w)(eta u_0)^T ),
          w = sum_{k != 0} sigma_k (v . u_k) / (lambda_0 - lambda_k) u_k.
  GATE: 4-point central finite differences along random symmetric
  directions (Richardson pair eps, eps/2) must agree with g . D to <= 1e-6
  relative before any relaxation is run (the `gate` stage, exit code 1
  otherwise; the complex step is not available because the eigenvector
  normalization is not holomorphic).

The eigenvector-degeneracy locus (the R1 audit cost): where the timelike
eigenvalue of M eta collides with a spacelike one the spectrum turns
complex and u is undefined; I1_h ~ (t* - t)^(-1/2) there. On a Lorentz
orbit of the vacuum the spectrum of M eta is EXACTLY (-g, 1, delta, 0),
so the dressing families cannot reach the locus by construction; the
`family` stage measures this (per amplitude: max |M_0i| over the grid,
the minimum spectral gap min_cells min_k |lambda_0 - lambda_k|, and the
existence of exactly one real timelike eigenvector in every cell) and
the free descent is where the locus can be hit. A descent that loses
the eigenvector is reported LOCUS-HIT (step, max |M_0i|), one that
turns non-finite or falls below the dive floor DIVERGED; nothing is
fitted.

Fixed-J electron (the M5.21.16 FIXJ instrument, lattice version, NO guard):
on the rigid dressing ladder b = amp tanh(r/2) applied to the M5.21.8
twisting hedgehog with its clock flow a0,
    E(amp; J) = E_stat_lambda(amp) + J^2 / (4 kin_lambda(amp)),
    omega* = J / (2 kin_lambda(amp*)),   dE/dJ = omega*   (closure)
with the J ladders of the record (M5.21.16: 50, 200, 800; M5.21.15:
J = 2 kin_base omega_t, omega_t = 0.1 .. 1.5); a ladder member with
kin <= 0 makes the fixed-J energy unbounded below (guard needed).

Concavity certificate (lambda = 1): per configuration E is affine in
omega^2 (slope kin), so the FREE envelope E_free(omega) = min_amp E(amp;
omega) is a minimum of affine functions of omega^2, hence concave in
omega^2: no free interior minimum unless a slope is negative (the
M5.21.15 theorem); the slopes and the envelope second differences are
reported.

STAGES (python3 m5_32_r2_b_bounded.py STAGE [--lam] [--n] [--L] [--steps]):
    gate                    identities + the gradient gate (must pass)
    family --n N            wide + sawtooth families at L = 24 (all lambdas
                            from one pass: the reads are linear in lambda)
    descent --lam X --n N --L L --steps S   the free FIRE descent
    fixedj --n N            the fixed-J electron + envelope (all lambdas)
    collect                 merge the checkpoints into the JSON + plots
Partials: ../checkpoints/m5_32_r2_b/*.json (gitignored)
Out: ../data/m5_32_r2_bounded.json, ../plots/m5_32_r2_bounded_*.png
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
CKPT = os.path.join(HERE, "..", "checkpoints", "m5_32_r2_b")


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


LAG = _load("m5_32_lagrangian", "m5_32_lagrangian.py")
B3 = LAG.B3
B8 = _load("m5_21_8_b_lattice", "m5_21_8_b_lattice.py")
ETA = B3.ETA
ETA_D = np.diag(ETA)
G_MAIN, S_MAIN, DELTA = 32.0, -1.0, 0.3
LAMBDAS = (0.0, 0.5, 0.75, 1.0)
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:8.1f}s] {msg}", flush=True)


def cfg_of(n, L):
    return B3.base_cfg(s=S_MAIN, g=G_MAIN, n=n, L=float(L))


def params():
    return LAG.default_params(s=S_MAIN, g=G_MAIN)


# ================= eigenframe =================
def tl_eig(M):
    """u0 (...,4), V (...,4,4) columns eta-normalized, lam (...,4),
    sig (...,4), k0 (...), ok (...) bool, gap (...) = min_k |l0 - l_k|.
    Same construction as m5_32_terms_ext.timelike_eig_np, per-cell
    existence returned as a mask instead of raising."""
    N = M @ ETA
    lam, V = np.linalg.eig(N)
    scale = np.maximum(np.max(np.abs(lam.real), axis=-1), 1.0)
    real_ok = np.max(np.abs(lam.imag), axis=-1) <= 1e-9 * scale
    lam = lam.real
    V = V.real
    n2 = np.einsum("...ak,a,...ak->...k", V, ETA_D, V)
    n_time = np.sum(n2 < 0, axis=-1)
    nz_ok = np.min(np.abs(n2), axis=-1) > 1e-10 * np.max(
        np.einsum("...ak,...ak->...k", V, V), axis=-1)
    ok = real_ok & (n_time == 1) & nz_ok
    safe = np.where(np.abs(n2) > 1e-300, np.abs(n2), 1.0)
    V = V / np.sqrt(safe)[..., None, :]
    sig = np.sign(n2)
    k0 = np.argmin(n2, axis=-1)
    u0 = np.take_along_axis(V, k0[..., None, None], axis=-1)[..., 0]
    l0 = np.take_along_axis(lam, k0[..., None], axis=-1)[..., 0]
    d = np.abs(l0[..., None] - lam)
    np.put_along_axis(d, k0[..., None], np.inf, axis=-1)
    gap = np.min(d, axis=-1)
    return u0, V, lam, sig, k0, ok, gap


def h_of(u0):
    hu = u0 @ ETA
    return ETA + 2.0 * hu[..., :, None] * hu[..., None, :]


def locus_report(M):
    """the three locus numbers of a field."""
    _, _, _, _, _, ok, gap = tl_eig(M)
    m0i = np.max(np.abs(M[..., 0, 1:]))
    return {"max_abs_M0i": float(m0i), "min_gap": float(np.min(gap)),
            "timelike_eigvec_everywhere": bool(np.all(ok)),
            "n_bad_cells": int(np.sum(~ok))}


# ================= energies =================
def q_eta(F):
    return B3.inner_eta(F, F)


def q_h(F, h):
    return np.einsum("...ab,...bc,...cd,...ad->...", h, F, h, F,
                     optimize=True)


def reads(M, cfg, a0=None):
    """(A_I1, A_I1h, A_I1frob, V4, kin_I1, kin_I1h, kin_frob): the
    h^3-weighted lattice reads, E_stat = 4 A + V4, E_kin = omega^2 kin.
    One eigen pass; all lambdas combine linearly."""
    h3 = cfg["h"] ** 3
    u0, _, _, _, _, ok, gap = tl_eig(M)
    hh = h_of(u0)
    acc = np.zeros(6)
    for br, wt in B3.branches(cfg["stencil"]):
        A = [B3.d1(M, ax, cfg["h"], br) for ax in range(3)]
        for i in range(3):
            for j in range(i + 1, 3):
                F = B3.comm_eta(A[i], A[j])
                acc[0] += wt * np.sum(q_eta(F))
                acc[1] += wt * np.sum(q_h(F, hh))
                acc[2] += wt * np.sum(F * F)
            if a0 is not None:
                F = B3.comm_eta(a0, A[i])
                acc[3] += wt * np.sum(q_eta(F))
                acc[4] += wt * np.sum(q_h(F, hh))
                acc[5] += wt * np.sum(F * F)
    _, ev = B3.e_parts(M, cfg)
    out = {"A_I1": h3 * acc[0], "A_I1h": h3 * acc[1],
           "A_I1frob": h3 * acc[2], "V4": float(ev),
           "kin_I1": 4 * h3 * acc[3], "kin_I1h": 4 * h3 * acc[4],
           "kin_frob": 4 * h3 * acc[5]}
    out.update(locus_report(M))
    return {k: float(v) for k, v in out.items()}


def e_stat_lam(rd, lam):
    return 4.0 * ((1 - lam) * rd["A_I1"] + lam * rd["A_I1h"]) + rd["V4"]


def kin_lam(rd, lam):
    return (1 - lam) * rd["kin_I1"] + lam * rd["kin_I1h"]


def energy_grad(M, cfg, lam):
    """E_stat_lambda and its exact gradient wrt symmetric M."""
    h3, h = cfg["h"] ** 3, cfg["h"]
    n = M.shape[0]
    u0, V, lamv, sig, k0, ok, gap = tl_eig(M)
    if not np.all(ok) and lam != 0.0:
        # the flip metric is undefined here (the locus); the certified
        # action (lambda = 0) needs no eigenframe and is never aborted
        return np.nan, None, {"ok": False, "n_bad": int(np.sum(~ok)),
                              "min_gap": float(np.min(gap))}
    hh = h_of(u0)
    E = 0.0
    G = np.zeros_like(M)
    S = np.zeros_like(M)
    for br, wt in B3.branches(cfg["stencil"]):
        A = [B3.d1(M, ax, h, br) for ax in range(3)]
        EA_T = [(ETA @ a).swapaxes(-1, -2) for a in A]
        AE_T = [(a @ ETA).swapaxes(-1, -2) for a in A]
        dA = [np.zeros_like(M) for _ in range(3)]
        for i in range(3):
            for j in range(i + 1, 3):
                F = B3.comm_eta(A[i], A[j])
                W = 0.0
                if lam != 1.0:
                    E += wt * (1 - lam) * np.sum(q_eta(F))
                    W = W + (1 - lam) * (ETA @ F @ ETA)
                if lam != 0.0:
                    hF = hh @ F
                    E += wt * lam * np.sum(np.einsum(
                        "...ab,...ab->...", hF @ hh, F))
                    W = W + lam * (hF @ hh)
                    S += wt * ((F @ hh) @ F.swapaxes(-1, -2))   # F h F^T
                W = 2.0 * W
                dA[i] += W @ EA_T[j] - AE_T[j] @ W
                dA[j] += AE_T[i] @ W - W @ EA_T[i]
        for i in range(3):
            G += wt * B3.d1_adj(dA[i], i, h, br)
    if lam != 0.0:
        # piece (b): through h(M); S already carries the branch weights
        v = 8.0 * np.einsum("ab,...bc,cd,...d->...a", ETA, S, ETA, u0)
        vu = np.einsum("...a,...ak->...k", v, V)
        l0 = np.take_along_axis(lamv, k0[..., None], axis=-1)[..., 0]
        den = l0[..., None] - lamv
        mask = np.arange(4)[None, :] != k0.reshape(-1, 1)
        mask = mask.reshape(den.shape)
        c = np.where(mask, sig * vu / np.where(mask, den, 1.0), 0.0)
        w = np.einsum("...ak,...k->...a", V, c)
        ew = w @ ETA
        eu = u0 @ ETA
        G += lam * (ew[..., :, None] * eu[..., None, :])
    G = 4.0 * h3 * B3.sym4(G)
    _, ev = B3.e_parts(M, cfg)
    E = 4.0 * h3 * E + ev
    G = G + h3 * B3.sym4(_v4_grad(M, cfg))
    return float(E), G, {"ok": True, "min_gap": float(np.min(gap)),
                         "max_abs_M0i": float(np.max(np.abs(M[..., 0, 1:])))}


def _v4_grad(M, cfg):
    p = params()
    return LAG.v4_grad_np(M, p)


# ================= fields =================
def boost_geom(cfg):
    X, Y, Z = B3.coords(cfg["n"], cfg["h"])
    R = np.sqrt(X * X + Y * Y + Z * Z)
    nx, ny, nz = X / R, Y / R, Z / R
    K = np.zeros(X.shape + (4, 4))
    K[..., 0, 1], K[..., 0, 2], K[..., 0, 3] = nx, ny, nz
    K[..., 1, 0], K[..., 2, 0], K[..., 3, 0] = nx, ny, nz
    K2 = np.zeros_like(K)
    K2[..., 0, 0] = 1.0
    for i, a in enumerate((nx, ny, nz)):
        for j, bb in enumerate((nx, ny, nz)):
            K2[..., 1 + i, 1 + j] = a * bb
    return R, K, K2


def dress(Mb, a0b, bl, K, K2):
    Qb = (np.eye(4)[None, None, None] + np.sinh(bl)[..., None, None] * K
          + (np.cosh(bl) - 1.0)[..., None, None] * K2)
    Md = B3.sym4(np.einsum("...ab,...bc,...dc->...ad", Qb, Mb, Qb))
    a0d = None if a0b is None else B3.sym4(
        np.einsum("...ab,...bc,...dc->...ad", Qb, a0b, Qb))
    return Md, a0d


def saw_of(A, lam, r):
    """m5_21_14_c_minimize.saw_of, verbatim."""
    return A * np.sin(np.pi * r / lam) * (r <= 8.0)


def bstar_record():
    with open(os.path.join(DATA, "m5_21_14_minimize.json")) as f:
        rec = json.load(f)
    return np.array(rec["rs"]), np.array(rec["b_star"])


def dressed_electron(cfg, scale=1.0):
    """the b*-dressed electron of the M5.21.14 lattice cross-check."""
    rs, bstar = bstar_record()
    R, K, K2 = boost_geom(cfg)
    bl = scale * np.interp(R.ravel(), rs, bstar).reshape(R.shape)
    Mb = B8.dressed(cfg, 0.0)
    a0b = B8.a0_unit(cfg, 0.0)
    return dress(Mb, a0b, bl, K, K2)


def ck(name):
    os.makedirs(CKPT, exist_ok=True)
    return os.path.join(CKPT, name)


def dump(name, obj):
    with open(ck(name), "w") as f:
        json.dump(obj, f, indent=1)


# ================= GATE =================
def stage_gate():
    EXT = _load("m5_32_terms_ext", "m5_32_terms_ext.py")
    out = {"n": 16, "L": 24.0, "stencil": "sym", "checks": {}}
    cfg = cfg_of(16, 24.0)
    p = params()
    rng = np.random.default_rng(3202)
    Md, a0d = dressed_electron(cfg)
    Mr = Md + 0.3 * B3.sym4(rng.standard_normal(Md.shape)) \
        * np.exp(-(boost_geom(cfg)[0] / 6.0) ** 2)[..., None, None]
    chk = out["checks"]
    # (a) lambda = 0 identities against the certified stack
    for nm, M in (("dressed", Md), ("off_orbit", Mr)):
        E, G, _ = energy_grad(M, cfg, 0.0)
        Eb, Gb = B3.e_total(M, cfg), B3.grad(M, cfg)
        rd = reads(M, cfg, a0d)
        chk[f"lam0_energy_vs_B3_{nm}"] = float(abs(E - Eb) / abs(Eb))
        chk[f"lam0_grad_vs_B3_{nm}"] = float(
            np.max(np.abs(G - Gb)) / np.max(np.abs(Gb)))
        chk[f"lam0_reads_vs_B3_{nm}"] = float(
            abs(e_stat_lam(rd, 0.0) - Eb) / abs(Eb))
        chk[f"lam0_kin_vs_B3_{nm}"] = float(
            abs(kin_lam(rd, 0.0) - B3.kin_of(M, a0d, cfg))
            / abs(B3.kin_of(M, a0d, cfg)))
        # (b) I1_h density / kin against the registry
        A = LAG.lattice_jets(M, cfg)[0][0]
        d_mine = 0.0
        u0, *_ = tl_eig(M)
        hh = h_of(u0)
        F = LAG.F_of_A(A)
        for i in range(1, 4):
            for j in range(i + 1, 4):
                d_mine = d_mine + q_h(F[..., i, j, :, :], hh)
        d_reg = EXT.I1_h_np(A, M, p)
        chk[f"I1h_density_vs_registry_{nm}"] = float(
            np.max(np.abs(d_mine - d_reg)) / np.max(np.abs(d_reg)))
        Ah, Bh, Ch = LAG.omega_decompose(EXT.REGISTRY_EXT["I1_h"], M,
                                         cfg, p, a0d)
        chk[f"I1h_A_vs_registry_{nm}"] = float(
            abs(rd["A_I1h"] - Ah) / abs(Ah))
        chk[f"I1h_kin_vs_registry_C_{nm}"] = float(
            abs(rd["kin_I1h"] - (-4.0 * Ch)) / abs(4.0 * Ch))
        chk[f"I1h_B_omega_odd_{nm}"] = float(abs(Bh) / abs(Ah))
        # energy at lambda = 1 vs the registry read
        E1, _, _ = energy_grad(M, cfg, 1.0)
        chk[f"lam1_energy_vs_registry_{nm}"] = float(
            abs(E1 - (4.0 * Ah + rd["V4"])) / abs(E1))
    # vacuum-eigenframe identity: I1_h = I1_frob on the undressed hedgehog
    Mb = B8.dressed(cfg, 0.0)
    rd = reads(Mb, cfg, B8.a0_unit(cfg, 0.0))
    chk["I1h_eq_I1frob_undressed"] = float(
        abs(rd["A_I1h"] - rd["A_I1frob"]) / abs(rd["A_I1frob"]))
    chk["I1h_eq_I1frob_undressed_kin"] = float(
        abs(rd["kin_I1h"] - rd["kin_frob"]) / abs(rd["kin_frob"]))
    # (c) the gradient gate: 4-point central FD along random directions
    fd = {}
    worst = 0.0
    for lam in (0.5, 1.0):
        for nm, M in (("dressed", Md), ("off_orbit", Mr)):
            E0, G, info = energy_grad(M, cfg, lam)
            rows = []
            for k in range(3):
                D = B3.sym4(rng.standard_normal(M.shape))
                D /= np.sqrt(np.sum(D * D))
                gd = float(np.sum(G * D))

                def e_at(t):
                    return energy_grad(M + t * D, cfg, lam)[0]

                def fd4(eps):
                    return (8 * (e_at(eps) - e_at(-eps))
                            - (e_at(2 * eps) - e_at(-2 * eps))) / (12 * eps)
                eps = 1e-3
                f1, f2 = fd4(eps), fd4(eps / 2)
                rich = (16 * f2 - f1) / 15.0
                rel = abs(rich - gd) / max(abs(gd), 1e-300)
                rows.append({"g_dot_D": gd, "fd_eps": f1, "fd_eps_half": f2,
                             "fd_richardson": rich, "rel_err": float(rel)})
                worst = max(worst, rel)
            fd[f"lam{lam:g}_{nm}"] = {"E": E0, "min_gap": info["min_gap"],
                                      "dirs": rows}
    out["fd_gate"] = fd
    out["fd_worst_rel_err"] = float(worst)
    out["gate_pass"] = bool(worst <= 1e-6 and all(
        v <= 1e-9 for k, v in chk.items() if "omega_odd" not in k))
    dump("gate.json", out)
    log(f"GATE worst FD rel err {worst:.3e}; identities "
        + ", ".join(f"{k} {v:.1e}" for k, v in chk.items()))
    log(f"GATE PASS = {out['gate_pass']}")
    return out


# ================= FAMILY =================
WIDE_AMPS = (0.02, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0)
SAW = tuple((lam, A) for lam in (2.0, 1.0, 0.5) for A in (0.02, 0.05))


def stage_family(n, L=24.0):
    cfg = cfg_of(n, L)
    R, K, K2 = boost_geom(cfg)
    Mb = B8.dressed(cfg, 0.0)
    a0b = B8.a0_unit(cfg, 0.0)
    out = {"n": n, "L": L, "h": cfg["h"], "stencil": cfg["stencil"],
           "g": G_MAIN, "s": S_MAIN, "delta": DELTA, "members": {}}
    base = reads(Mb, cfg, a0b)
    out["members"]["base"] = base
    log(f"FAMILY n{n} base: A_I1 {base['A_I1']:.6g} A_I1h {base['A_I1h']:.6g} "
        f"V4 {base['V4']:.3g} gap {base['min_gap']:.4g}")
    members = [(f"wide_amp_{a:g}", a * np.tanh(R / 2.0)) for a in WIDE_AMPS]
    members += [(f"saw_lam_{l:g}_A_{A:g}", saw_of(A, l, R)) for l, A in SAW]
    for name, bl in members:
        Md, a0d = dress(Mb, a0b, bl, K, K2)
        rd = reads(Md, cfg, a0d)
        out["members"][name] = rd
        log(f"FAMILY n{n} {name}: dE(lam0) {e_stat_lam(rd, 0) - e_stat_lam(base, 0):+.6g} "
            f"dE(lam1) {e_stat_lam(rd, 1) - e_stat_lam(base, 1):+.6g} "
            f"max|M0i| {rd['max_abs_M0i']:.4g} gap {rd['min_gap']:.4g} "
            f"ok {rd['timelike_eigvec_everywhere']} V4 {rd['V4']:.2e}")
    # per-lambda summaries
    summ = {}
    for lam in LAMBDAS:
        e0 = e_stat_lam(base, lam)
        wide = {nm: e_stat_lam(out["members"][nm], lam) - e0
                for nm, _ in members if nm.startswith("wide")}
        ws = list(wide.values())
        imin = int(np.argmin(ws))
        row = {"wide_dE": wide,
               "wide_min": float(ws[imin]), "wide_argmin": list(wide)[imin],
               "wide_interior_min": bool(0 <= imin < len(ws) - 1 and
                                         (imin > 0 or ws[0] < 0)),
               "wide_monotone_falling": bool(all(
                   ws[k + 1] < ws[k] for k in range(len(ws) - 1))),
               "wide_goes_negative": bool(min(ws) < 0),
               "wide_last_positive": bool(ws[-1] > 0)}
        for l in (2.0, 1.0, 0.5):
            sw = {nm: e_stat_lam(out["members"][nm], lam) - e0
                  for nm, _ in members if nm.startswith(f"saw_lam_{l:g}_")}
            row[f"saw_lam_{l:g}_dE"] = sw
            row[f"saw_lam_{l:g}_negative"] = bool(min(sw.values()) < 0)
        summ[f"lam_{lam:g}"] = row
    out["per_lambda"] = summ
    ok_all = all(m["timelike_eigvec_everywhere"] for m in out["members"].values())
    out["locus_amplitude"] = None if ok_all else next(
        a for a in WIDE_AMPS
        if not out["members"][f"wide_amp_{a:g}"]["timelike_eigvec_everywhere"])
    out["locus_note"] = ("every family member is a per-cell Lorentz orbit "
                         "of the vacuum, so the spectrum of M eta is exactly "
                         "(-g, 1, delta, 0) and the locus is unreachable on "
                         "the family" if ok_all else "locus reached")
    dump(f"family_n{n}.json", out)
    return out


# ================= DESCENT =================
def stage_descent(lam, n, L, steps, log_every=50, dt0=0.02, dt_max=0.2,
                  dive_floor=-1e6):
    cfg = cfg_of(n, L)
    Md, a0d = dressed_electron(cfg)
    free = (~B3.pin_shell(n, cfg["h"]))[..., None, None].astype(float)
    M = Md.copy()
    E0, G, info = energy_grad(M, cfg, lam)
    rd0 = reads(M, cfg, a0d)
    out = {"lam": lam, "n": n, "L": L, "h": cfg["h"], "stencil": "sym",
           "steps_budget": steps, "pin": "B3.pin_shell depth 1.6 (Dirichlet "
           "at the initial b*-dressed values)", "E0": E0, "E0_lam0_control":
           float(B3.e_total(M, cfg)), "start": rd0, "trace": [],
           "fire": {"dt0": dt0, "dt_max": dt_max, "alpha0": 0.1}}
    # FIRE with ENERGY-MONOTONE BACKTRACKING: a trial step is accepted only
    # if E_lambda decreases (and, for lambda > 0, the eigenframe exists in
    # every cell); a rejected step halves dt and resets the velocity. The
    # certified fire (B3.fire) has no energy check and the record's own
    # relax from the lattice hedgehog went non-finite at g = 8 (M5.21.8
    # relax records, dt0 down to 1e-4): the g = 32 core is too stiff for a
    # fixed dt. LOCUS-HIT is declared only when dt collapses (< dt_min)
    # against locus-loss rejections: the descent cannot proceed without
    # leaving the domain of the flip metric.
    v = np.zeros_like(M)
    dt, alpha, n_up = dt0, 0.1, 0
    dt_min = 1e-7
    F = -G * free
    E_prev = E0
    stop = "budget"
    tag = f"lam{lam:g}_n{n}_L{L:g}" + ("" if steps == 2000 else f"_s{steps}")
    log(f"DESCENT {tag} E0 {E0:.6f} gap {info['min_gap']:.4g} "
        f"max|M0i| {info['max_abs_M0i']:.4g}")
    it_stop = steps
    n_rej, n_rej_locus, n_acc = 0, 0, 0
    last_locus = None
    fmax = float(np.max(np.abs(F)))
    for it in range(1, steps + 1):
        P = float(np.sum(F * v))
        if P > 0.0:
            n_up += 1
            vn = np.sqrt(np.sum(v * v))
            fn = np.sqrt(np.sum(F * F))
            v = (1 - alpha) * v + alpha * (F / max(fn, 1e-300)) * vn
            if n_up > 5:
                dt = min(dt * 1.1, dt_max)
                alpha *= 0.99
        else:
            v[:] = 0.0
            alpha = 0.1
            n_up = 0
        v_try = v + dt * F
        M_try = M + dt * v_try
        E, G, info = energy_grad(M_try, cfg, lam)
        locus_loss = not info["ok"]
        reject = locus_loss or not np.isfinite(E) or \
            E > E_prev + 1e-12 * max(abs(E_prev), 1.0)
        if reject:
            n_rej += 1
            if locus_loss:
                n_rej_locus += 1
                last_locus = {"step": it, "n_bad_cells": info["n_bad"],
                              "min_gap_trial": info["min_gap"],
                              "max_abs_M0i_trial":
                                  float(np.max(np.abs(M_try[..., 0, 1:]))),
                              "dt": dt}
            dt *= 0.5
            v[:] = 0.0
            alpha, n_up = 0.1, 0
            if dt < dt_min:
                stop = ("LOCUS-HIT" if locus_loss else
                        "STALLED (dt collapsed, no descent direction "
                        "accepted)")
                it_stop = it
                break
            continue
        n_acc += 1
        M, v, E_prev = M_try, v_try, E
        F = -G * free
        fmax = float(np.max(np.abs(F)))
        if it % log_every == 0 or it == steps:
            row = {"it": it, "E": E, "fmax": fmax, "dt": dt,
                   "min_gap": info["min_gap"],
                   "max_abs_M0i": info["max_abs_M0i"],
                   "n_accepted": n_acc, "n_rejected": n_rej,
                   "n_rejected_locus": n_rej_locus}
            out["trace"].append(row)
            log(f"DESCENT {tag} it {it:5d} E {E:14.4f} fmax {fmax:.3e} "
                f"gap {info['min_gap']:.4g} max|M0i| {info['max_abs_M0i']:.4g} "
                f"dt {dt:.2e} acc {n_acc} rej {n_rej} (locus {n_rej_locus})")
            if E < dive_floor:
                stop = "DIVERGED (dive floor)"
                it_stop = it
                break
    if stop == "LOCUS-HIT":
        out["locus_hit"] = dict(last_locus, E_last_accepted=E_prev,
                                min_gap_last_accepted=float(np.min(tl_eig(M)[6])),
                                max_abs_M0i_last_accepted=float(np.max(np.abs(M[..., 0, 1:]))))
        log(f"DESCENT {tag} LOCUS-HIT: {out['locus_hit']}")
    out["counts"] = {"accepted": n_acc, "rejected": n_rej,
                     "rejected_locus": n_rej_locus, "dt_final": dt}
    out["last_locus_rejection"] = last_locus
    out["E_last_accepted"] = float(E_prev)
    out["stop"] = stop
    out["steps_run"] = it_stop
    tr = out["trace"]
    if tr:
        last = tr[-1]
        out["E_end"] = last["E"]
        out["max_abs_M0i_end"] = last["max_abs_M0i"]
        out["min_gap_end"] = last["min_gap"]
        q = [r for r in tr if r["it"] >= 0.75 * it_stop]
        if len(q) >= 2 and stop == "budget":
            dE = q[-1]["E"] - q[0]["E"]
            out["last_quarter_dE"] = float(dE)
            out["last_quarter_rel"] = float(abs(dE) / max(abs(q[-1]["E"]), 1.0))
            out["verdict"] = ("PLATEAU" if out["last_quarter_rel"] <= 1e-3
                              else "FALLING (still descending at the budget)"
                              if dE < 0 else "RISING")
        elif stop == "budget":
            out["verdict"] = "budget (trace too short)"
    if stop != "budget":
        out["verdict"] = stop
    if stop == "budget":
        out["end"] = reads(M, cfg, a0d)
    out["wall_s"] = round(time.time() - T0, 1)
    dump(f"descent_{tag}.json", out)
    log(f"DESCENT {tag} verdict {out['verdict']} E_end "
        f"{out.get('E_end', float('nan')):.4f}")
    return out


# ================= FIXED-J =================
AMPS = (0.0, 0.0025, 0.005, 0.0075, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05,
        0.07, 0.1, 0.15, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0)
J_16 = (50.0, 200.0, 800.0)
OMEGA_T15 = (0.1, 0.2, 0.5, 1.0, 1.5)


def fixedj_table(rows, lam, J):
    """min over the amp ladder of E_stat + J^2/(4 kin), NO guard."""
    Es, kins, amps = [], [], []
    for r in rows:
        k = kin_lam(r, lam)
        amps.append(r["amp"])
        kins.append(k)
        Es.append(e_stat_lam(r, lam) + (J * J / (4.0 * k) if k > 0 else -np.inf))
    Es = np.array(Es)
    neg = [a for a, k in zip(amps, kins) if k <= 0]
    if neg:
        return {"J": J, "unbounded": True, "kin_nonpositive_at_amps": neg,
                "guard_needed": True}
    i = int(np.argmin(Es))
    return {"J": J, "unbounded": False, "opt_amp": amps[i],
            "E_total": float(Es[i]), "kin_opt": float(kins[i]),
            "omega_star": float(J / (2.0 * kins[i])),
            "E_positive": bool(Es[i] > 0),
            "interior": bool(0 < i < len(amps) - 1),
            "at_ladder_edge": bool(i == len(amps) - 1),
            "guard_needed": bool(i == len(amps) - 1 or Es[i] < 0),
            "E_ladder": Es.tolist()}


def refine_fixedj(read_at, rows, lam, J, iters=12):
    """golden-section search of E_stat + J^2/(4 kin) over amp between the
    ladder neighbors of the discrete optimum (real lattice reads, cached);
    returns the refined row or the discrete one if unbounded."""
    t = fixedj_table(rows, lam, J)
    if t["unbounded"]:
        return t
    amps = [r["amp"] for r in rows]
    i = amps.index(t["opt_amp"])
    lo = amps[i - 1] if i > 0 else 0.0
    hi = amps[i + 1] if i < len(amps) - 1 else amps[i]

    def obj(a):
        r = read_at(a)
        k = kin_lam(r, lam)
        return (e_stat_lam(r, lam) + J * J / (4.0 * k)) if k > 0 else np.inf, r
    gr = (np.sqrt(5.0) - 1.0) / 2.0
    a, b = lo, hi
    c, d = b - gr * (b - a), a + gr * (b - a)
    fc, fd = obj(c)[0], obj(d)[0]
    for _ in range(iters):
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - gr * (b - a)
            fc = obj(c)[0]
        else:
            a, c, fc = c, d, fd
            d = a + gr * (b - a)
            fd = obj(d)[0]
    best_a = c if fc < fd else d
    fbest, r = obj(best_a)
    # never worse than the discrete optimum
    if fbest > t["E_total"]:
        best_a, fbest, r = t["opt_amp"], t["E_total"], read_at(t["opt_amp"])
    k = kin_lam(r, lam)
    t.update({"opt_amp_refined": float(best_a), "E_total": float(fbest),
              "kin_opt": float(k), "omega_star": float(J / (2.0 * k)),
              "E_positive": bool(fbest > 0), "bracket": [lo, hi],
              "E_stat_opt": float(e_stat_lam(r, lam)),
              "E_rot_opt": float(J * J / (4.0 * k))})
    return t


def envelope_of(rows, kb):
    """the FREE omega envelope per lambda: E_free(omega) = min_amp
    [E_stat(amp) + omega^2 kin(amp)], a minimum of affine functions of
    x = omega^2 (hence concave in x); the second DIVIDED differences in x
    certify it (the omega grid is uniform, x is not)."""
    env = {}
    om = np.linspace(0.0, 3.0, 31)
    x = om ** 2
    for lam in LAMBDAS:
        Es = np.array([e_stat_lam(r, lam) for r in rows])
        ks = np.array([kin_lam(r, lam) for r in rows])
        Efree = np.array([np.min(Es + w * w * ks) for w in om])
        arg = [int(np.argmin(Es + w * w * ks)) for w in om]
        s1 = np.diff(Efree) / np.diff(x)
        d2 = np.diff(s1) / (x[2:] - x[:-2])
        i = int(np.argmin(Efree))
        env[f"lam_{lam:g}"] = {
            "omega_grid": om.tolist(), "E_free": Efree.tolist(),
            "argmin_amp": [rows[a]["amp"] for a in arg],
            "slopes_dE_domega2_per_amp": ks.tolist(),
            "envelope_slopes_along_omega": s1.tolist(),
            "all_slopes_positive": bool(np.all(ks > 0)),
            "envelope_concave_in_omega2": bool(np.all(
                d2 <= 1e-9 * max(1.0, np.max(np.abs(Efree))))),
            "second_divided_diff_max": float(np.max(d2)),
            "free_min_omega": float(om[i]),
            "free_interior_minimum": bool(0 < i < len(om) - 1),
            "E_free_min": float(Efree[i]),
            "E_free_unbounded_in_omega": bool(np.any(ks < 0))}
        Jm = 2.0 * kb[f"lam_{lam:g}"] * 0.5
        fom = []
        for w in om[1:]:
            kt = Jm / (2.0 * w)
            if ks.min() < kt < ks.max():
                j = int(np.argmin(np.abs(ks - kt)))
                fom.append({"omega": float(w), "kin_target": float(kt),
                            "amp": rows[j]["amp"], "kin": float(ks[j]),
                            "E": float(Es[j] + w * w * ks[j]),
                            "feasible": True})
            else:
                fom.append({"omega": float(w), "kin_target": float(kt),
                            "feasible": False})
        env[f"lam_{lam:g}"]["fixedJ_mid"] = {"J": Jm, "rows": fom}
    return env


def stage_fixedj(n, L=48.0):
    cfg = cfg_of(n, L)
    R, K, K2 = boost_geom(cfg)
    Mb = B8.dressed(cfg, 0.0)
    a0b = B8.a0_unit(cfg, 0.0)
    cache = {}

    def read_at(amp):
        key = round(float(amp), 9)
        if key not in cache:
            Md, a0d = dress(Mb, a0b, key * np.tanh(R / 2.0), K, K2)
            rd = reads(Md, cfg, a0d)
            rd["amp"] = key
            cache[key] = rd
        return cache[key]
    rows = []
    for amp in AMPS:
        rd = read_at(amp)
        rows.append(rd)
        log(f"FIXJ n{n} amp {amp:g}: E0 {e_stat_lam(rd, 0):.4f} E1 "
            f"{e_stat_lam(rd, 1):.4f} kin0 {kin_lam(rd, 0):.4f} kin1 "
            f"{kin_lam(rd, 1):.4f} kin_frob {rd['kin_frob']:.4f} "
            f"E_frob {4 * rd['A_I1frob'] + rd['V4']:.4f}")
    kb = {f"lam_{lam:g}": kin_lam(rows[0], lam) for lam in LAMBDAS}
    out = {"n": n, "L": L, "h": cfg["h"], "stencil": "sym", "amps": list(AMPS),
           "rows": rows, "kin_base": kb, "E_stat_base":
           {f"lam_{lam:g}": e_stat_lam(rows[0], lam) for lam in LAMBDAS},
           "J_hbar_half_note": ("the record defines no hbar/2 in program "
                                "units (plan E4: program units only); the "
                                "J ladders are the record's (M5.21.16: 50, "
                                "200, 800; M5.21.15: J = 2 kin_base omega_t)"),
           "byLambda": {}}
    for lam in LAMBDAS:
        Js = list(J_16) + [2.0 * kb[f"lam_{lam:g}"] * w for w in OMEGA_T15]
        tab = {}
        for J in Js:
            t = refine_fixedj(read_at, rows, lam, J)
            if not t["unbounded"]:
                # closure dE/dJ = omega* (central differences, +-5 %, each
                # end re-minimized over amp by the same golden search)
                tp = refine_fixedj(read_at, rows, lam, 1.05 * J)
                tm = refine_fixedj(read_at, rows, lam, 0.95 * J)
                d = (tp["E_total"] - tm["E_total"]) / (0.10 * J)
                t["dEdJ_numeric"] = float(d)
                t["closure_rel"] = float(abs(d / t["omega_star"] - 1.0))
                t["closure_le_3pct"] = bool(t["closure_rel"] <= 0.03)
            tab[f"J_{J:.6g}"] = t
        neg_amps = [r["amp"] for r in rows if kin_lam(r, lam) <= 0]
        out["byLambda"][f"lam_{lam:g}"] = {
            "J_ladder": Js, "kin_ladder": [kin_lam(r, lam) for r in rows],
            "E_stat_ladder": [e_stat_lam(r, lam) for r in rows],
            "kin_nonpositive_amps": neg_amps,
            "all_kin_positive": bool(not neg_amps), "table": tab}
        log(f"FIXJ n{n} lam {lam:g}: " + "; ".join(
            f"J {J:.5g}: " + ("UNBOUNDED" if t["unbounded"] else
                              f"amp* {t['opt_amp']:g} E {t['E_total']:.2f} "
                              f"w* {t['omega_star']:.4f} clos {t['closure_rel']:.3f}")
            for J, t in zip(Js, tab.values())))
    # variant-A (I1_frob) comparison at lambda = 1: same fields
    va = {"per_amp": []}
    for r in rows:
        va["per_amp"].append({
            "amp": r["amp"],
            "E_stat_h": e_stat_lam(r, 1.0),
            "E_stat_frob": 4 * r["A_I1frob"] + r["V4"],
            "kin_h": r["kin_I1h"], "kin_frob": r["kin_frob"],
            "dE_h_minus_frob": e_stat_lam(r, 1.0) - (4 * r["A_I1frob"] + r["V4"]),
            "dkin_h_minus_frob": r["kin_I1h"] - r["kin_frob"]})
    def as_frob(r):
        return dict(r, A_I1h=r["A_I1frob"], kin_I1h=r["kin_frob"])
    frob_rows = [as_frob(r) for r in rows]
    va["fixedJ_frob"] = {f"J_{J:g}": refine_fixedj(
        lambda a: as_frob(read_at(a)), frob_rows, 1.0, J) for J in J_16}
    va["fixedJ_h"] = {f"J_{J:g}": out["byLambda"]["lam_1"]["table"][f"J_{J:.6g}"]
                      for J in J_16}
    with open(os.path.join(DATA, "m5_21_16_fixedj.json")) as f:
        rec16 = json.load(f)["FIXJ"]
    va["record_m5_21_16_quadrature"] = {
        "e_stat_base": rec16["e_stat_base"], "kin_base_flip": rec16["kin_base_flip"],
        "byJ": {k: {kk: vv for kk, vv in v.items() if kk != "ladder"}
                for k, v in rec16["byJ"].items()},
        "note": "quadrature family units (a0 un-normalized, rmax 24); the "
                "lattice numbers above are h^3 sums on the n, L box"}
    out["variant_A_comparison"] = va
    out["envelope"] = envelope_of(rows, kb)
    dump(f"fixedj_n{n}.json", out)
    return out


# ================= COLLECT + PLOTS =================
def stage_collect():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out = {"task": "M5.32 R2 arm (b): lattice boundedness + clock preview",
           "candidate": "L_lambda = -4 [(1 - lambda) I1 + lambda I1_h] - V4",
           "lambdas": list(LAMBDAS), "point": {"g": G_MAIN, "s": S_MAIN,
                                                "delta": DELTA, "stencil": "sym"}}
    def rd(name):
        p = ck(name)
        if os.path.exists(p):
            with open(p) as f:
                return json.load(f)
        return None
    out["gate"] = rd("gate.json")
    out["family"] = {os.path.basename(p)[:-5]: json.load(open(p))
                     for p in sorted(glob.glob(ck("family_n*.json")))}
    out["descent"] = {os.path.basename(p)[8:-5]: json.load(open(p))
                      for p in sorted(glob.glob(ck("descent_*.json")))}
    out["fixedj"] = {os.path.basename(p)[:-5]: json.load(open(p))
                     for p in sorted(glob.glob(ck("fixedj_n*.json")))}
    for fx in out["fixedj"].values():
        fx["envelope"] = envelope_of(fx["rows"], fx["kin_base"])
    # G5 verdict table
    g5 = {}
    for lam in LAMBDAS:
        key = f"lam_{lam:g}"
        row = {}
        for fn, fam in out["family"].items():
            s = fam["per_lambda"][key]
            row[f"wide_{fn}"] = {
                "min_dE": s["wide_min"], "argmin": s["wide_argmin"],
                "monotone_falling": s["wide_monotone_falling"],
                "goes_negative": s["wide_goes_negative"],
                "verdict": ("UNBOUNDED (falls monotonically through the "
                            "ladder)" if s["wide_monotone_falling"] and
                            s["wide_goes_negative"] else
                            "bounded on the ladder (interior well)" if
                            s["wide_goes_negative"] else
                            "bounded (no negative member)")}
            for l in (2.0, 1.0, 0.5):
                row[f"saw_lam_{l:g}_{fn}"] = {
                    "dE": s[f"saw_lam_{l:g}_dE"],
                    "negative": s[f"saw_lam_{l:g}_negative"]}
        for dn, d in out["descent"].items():
            if abs(d["lam"] - lam) < 1e-12:
                row[f"descent_{dn}"] = {
                    "verdict": d["verdict"], "E0": d["E0"],
                    "E_end": d.get("E_end"), "steps_run": d["steps_run"],
                    "max_abs_M0i_end": d.get("max_abs_M0i_end"),
                    "min_gap_end": d.get("min_gap_end"),
                    "locus_hit": d.get("locus_hit")}
        g5[key] = row
    out["G5"] = g5
    out["locus_amplitude"] = {fn: fam["locus_amplitude"]
                              for fn, fam in out["family"].items()}
    os.makedirs(PLOTS, exist_ok=True)
    # plot 1: families
    if out["family"]:
        fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
        for fn, fam in out["family"].items():
            for lam in LAMBDAS:
                s = fam["per_lambda"][f"lam_{lam:g}"]
                ax[0].plot(WIDE_AMPS, list(s["wide_dE"].values()), "o-",
                           label=f"lam {lam:g} n{fam['n']}")
            m = fam["members"]
            ax[1].semilogy(WIDE_AMPS, [m[f"wide_amp_{a:g}"]["max_abs_M0i"]
                                       for a in WIDE_AMPS], "s-",
                           label=f"max|M0i| n{fam['n']}")
            ax[2].plot(WIDE_AMPS, [m[f"wide_amp_{a:g}"]["min_gap"]
                                   for a in WIDE_AMPS], "^-",
                       label=f"min gap n{fam['n']}")
        ax[0].set_xscale("log"); ax[0].set_yscale("symlog", linthresh=100)
        ax[0].axhline(0, color="k", lw=0.5)
        ax[0].set_xlabel("amp (b = amp tanh(r/2))"); ax[0].set_ylabel("E_stat - E_base")
        ax[0].set_title("wide boost-dressing family, L = 24, NO guard")
        ax[0].legend(fontsize=7)
        ax[1].set_xscale("log"); ax[1].set_xlabel("amp"); ax[1].legend()
        ax[1].set_title("time-mixing entries max |M_0i|")
        ax[2].set_xscale("log"); ax[2].set_xlabel("amp"); ax[2].legend()
        ax[2].set_title("spectral gap of M eta (locus at 0)")
        fig.suptitle("M5.32 R2.b: G5 family probes of the lambda-family")
        fig.tight_layout()
        fig.savefig(os.path.join(PLOTS, "m5_32_r2_bounded_family.png"), dpi=110)
    # plot 2: descents
    if out["descent"]:
        fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
        for dn, d in out["descent"].items():
            tr = d["trace"]
            if not tr:
                continue
            its = [r["it"] for r in tr]
            ax[0].plot(its, [r["E"] for r in tr], label=f"{dn}: {d['verdict'][:22]}")
            ax[1].semilogy(its, [r["max_abs_M0i"] for r in tr], label=dn)
            ax[2].semilogy(its, [r["min_gap"] for r in tr], label=dn)
        ax[0].set_yscale("symlog", linthresh=1e3)
        ax[0].set_xlabel("step"); ax[0].set_ylabel("E_stat_lambda"); ax[0].legend(fontsize=7)
        ax[0].set_title("free FIRE descent, pinned shell, from the b*-dressed electron")
        ax[1].set_xlabel("step"); ax[1].set_title("max |M_0i|"); ax[1].legend(fontsize=7)
        ax[2].set_xlabel("step"); ax[2].set_title("min spectral gap (locus at 0)")
        ax[2].legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(os.path.join(PLOTS, "m5_32_r2_bounded_descent.png"), dpi=110)
    # plot 3: fixed-J + envelope
    if out["fixedj"]:
        fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
        for fn, fx in out["fixedj"].items():
            amps = fx["amps"]
            for lam in LAMBDAS:
                b = fx["byLambda"][f"lam_{lam:g}"]
                ax[0].plot(amps, b["kin_ladder"], "o-", label=f"kin lam {lam:g} {fn}")
                t = b["table"]["J_200"]
                if not t["unbounded"]:
                    ax[1].plot(amps, t["E_ladder"], "o-", label=f"lam {lam:g} {fn}")
                e = fx["envelope"][f"lam_{lam:g}"]
                ax[2].plot(np.array(e["omega_grid"]) ** 2, e["E_free"], "-",
                           label=f"lam {lam:g} {fn}")
        ax[0].set_xscale("symlog", linthresh=0.01); ax[0].set_yscale("symlog", linthresh=100)
        ax[0].axhline(0, color="k", lw=0.5); ax[0].set_xlabel("amp"); ax[0].legend(fontsize=6)
        ax[0].set_title("clock inertia kin(amp)")
        ax[1].set_xscale("symlog", linthresh=0.01); ax[1].set_yscale("symlog", linthresh=1e3)
        ax[1].set_xlabel("amp"); ax[1].set_title("fixed-J energy at J = 200"); ax[1].legend(fontsize=6)
        ax[2].set_xlabel("omega^2"); ax[2].set_yscale("symlog", linthresh=1e3)
        ax[2].set_title("free envelope E_free(omega) (concave in omega^2)"); ax[2].legend(fontsize=6)
        fig.tight_layout()
        fig.savefig(os.path.join(PLOTS, "m5_32_r2_bounded_fixedj.png"), dpi=110)
    out["collected_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    with open(os.path.join(DATA, "m5_32_r2_bounded.json"), "w") as f:
        json.dump(out, f, indent=1)
    log("COLLECT written")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["gate", "family", "descent", "fixedj", "collect"])
    ap.add_argument("--lam", type=float, default=1.0)
    ap.add_argument("--n", type=int, default=32)
    ap.add_argument("--L", type=float, default=48.0)
    ap.add_argument("--steps", type=int, default=2000)
    a = ap.parse_args()
    if a.stage == "gate":
        ok = stage_gate()["gate_pass"]
        sys.exit(0 if ok else 1)
    elif a.stage == "family":
        stage_family(a.n, 24.0)
    elif a.stage == "descent":
        stage_descent(a.lam, a.n, a.L, a.steps)
    elif a.stage == "fixedj":
        stage_fixedj(a.n, a.L)
    else:
        stage_collect()
    log(f"done {a.stage}")


if __name__ == "__main__":
    main()
