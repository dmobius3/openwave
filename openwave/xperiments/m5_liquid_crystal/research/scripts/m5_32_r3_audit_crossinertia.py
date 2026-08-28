"""M5.32 R3 arm (iii) INDEPENDENT ADVERSARIAL AUDIT: the un-relaxed
two-clock cross-inertia C(d) of a like-charge electron pair, rebuilt from
the oracles with a DIFFERENT window family and a second clock generator.
The producer's script and JSON were NOT read.

EQUATIONS FIRST
---------------
Conventions (the registries m5_32_lagrangian.py + m5_32_terms_ext.py and
the certified stack m5_21_3_a_4d.py, imported, never modified): index 0 =
time, eta = diag(-1, 1, 1, 1), M(x) real symmetric 4x4 per cell, spatial
jets A_i = d_i M on the certified sym stencil (fwd / bwd branches, density
per branch, averaged), A_0 the clock jet,
    F_{0i}(a) = a eta A_i - A_i eta a                (linear in a = A_0)
    I1   = sum_{mu<nu} eta^mu eta^nu <F, F>_eta,  <F,G>_eta = tr(eta F eta G^T)
    I1_h = sum_{mu<nu} eta^mu eta^nu tr(h F h F^T),  h = eta + 2 (eta u)(eta u)^T,
           u the timelike unit eigenvector of M eta (m5_32_terms_ext.h_cov_np)
    L_lambda = -4 [(1 - lambda) I1 + lambda I1_h] - V4
Because eta^0 = -1, the time-row part of both invariants is
-sum_i q_lambda(F_0i), so the Legendre read of the kinetic energy is the
bilinear form (h^3-weighted, branch-averaged)
    K_lambda(a, b) = 4 h^3 sum_br wt sum_cells sum_i
                     [(1 - lambda) <F_0i(a), F_0i(b)>_eta
                      + lambda tr(h F_0i(a) h F_0i(b)^T)]
    E_kin(A_0 = omega a) = omega^2 K_lambda(a, a)
(gate: K_0(a, a) == B3.kin_of(M, a) and == -4 C_I1 of LAG.omega_decompose;
K_1(a, a) == -4 C_{I1_h} of the ext registry, all to roundoff).

Two clocks: A_0 = w1 a0_1 + w2 a0_2 gives the quadratic form
    T = K11 w1^2 + K22 w2^2 + 2 K12 w1 w2,   K_kl = K_lambda(a0_k, a0_l)
    self-inertia I_k = 2 K_kk,  cross-inertia C = 2 K12,  I0 = 2 K(single)
    fixed-J like clocks   E_J = E_stat + j^2 / (I0 + C)   (attraction iff C > 0
                                                            grows as d shrinks)
    fixed-omega           E_w = E_stat + w^2 (I0 + C)      (opposite sign)
Ratios reported: C/I0 and I_1/I0 (positive rescalings of a0 cancel).

Clock generators (per core k, center x_k on the z axis, n_k = (x - x_k)/r_k):
    a0_k = w_k(x) (G_k M - M G_k^T) = w_k [G_k, M]      (G_k antisymmetric)
    rad   : G_k = local_rot(n_k), the record generator (rotation about the
            local radial direction from core k; m5_21_3_a_4d.gen_catalog)
    z     : G_k = J_z for both cores (rotation about the pair axis; the
            two clocks then differ only by their windows)
    local : G_k = local_rot(leading spatial eigenvector of M) (gen_catalog
            clock_local; identical for both cores, window-only difference)
Un-normalized: omega is the physical rotation rate (unit-omega clock).

Windows (a DIFFERENT family from the producer's tanh midplane PoU / Gaussian):
    vor   : hard Voronoi split at the midplane z = 0 smoothed over one cell,
            w1 = 1/2 (1 + tanh(2 z / h)),  w2 = 1 - w1   (partition of unity)
    sph   : spherical window of radius d/2 about each core smoothed over one
            cell, w_k = 1/2 (1 + tanh((d/2 - r_k) / h))    (NOT a PoU)
    sph6  : fixed-radius sphere R_w = 6 about each core, same smoothing
    none  : w_k = 1 (the generator alone separates the clocks)
Single-electron reference for each window: the same window placed about
the single core (vor: plane at distance d/2 from the core; sph: radius d/2
about the core), so I_1(d)/I0_w(d) isolates the composite-director effect
from the truncation; I0_full = the un-windowed single (the producer's I0).

Fields (un-relaxed, the arm's seed construction):
    3x3 part: the certified two-center composition m5_21_4_a_pair.seed_pair
              ('same' = like charges, product ansatz + escape tube; 'anti' as
              the control; 'single'), isotropic-blended cores r_c = 4,
              embedded with M_00 = g (B3.embed34)
    4x4 part: per core the RIGID boost dressing of the R2 arm (b) ladder,
              rapidity b_k(x) = amp tanh(r_k / 2) along n_k,
              Q_k = 1 + sinh(b_k) K_k + (cosh(b_k) - 1) K_k^2,  Q = Q_1 Q_2,
              M = sym4(Q M_emb Q^T);  amp = 0.022 (the producer's lambda = 1,
              J = 200 point), amp = 0.10 as the amplitude-sensitivity control
Toy point: s = -1, g = 32, delta = 0.3 (M_vac = diag(32, 1, 0.3, 0)).

Reads per (box, kind, amp, window, generator, d):
    C/I0_full, I_1/I0_full, I_1/I0_w, and for the single I0_w/I0_full.
Power-law screen (X4): over d >= 12 with a constant sign of C, fit
log|C/I0_full| = a + p log d; PASS iff 0.8 <= |p| <= 1.2 and R^2 >= 0.95 in
BOTH boxes (n = 32, L = 48 and n = 48, L = 72, h = 1.5).

Out: ../data/m5_32_r3_audit_crossinertia.json
Run: nice -n 10 python3 m5_32_r3_audit_crossinertia.py [--quick]
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
OUT_JSON = os.path.join(DATA, "m5_32_r3_audit_crossinertia.json")


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


B3 = _load("m5_21_3_a_4d", "m5_21_3_a_4d.py")
LAG = _load("m5_32_lagrangian", "m5_32_lagrangian.py")
EXT = _load("m5_32_terms_ext", "m5_32_terms_ext.py")
PAIR = _load("m5_21_4_a_pair", "m5_21_4_a_pair.py")
ETA = B3.ETA

G_MAIN, S_MAIN, DELTA = 32.0, -1.0, 0.3
BOXES = ((32, 48.0), (48, 72.0))
DS_BOX = {32: (8.0, 10.0, 12.0, 14.0, 18.0, 24.0),
          48: (8.0, 12.0, 14.0, 18.0, 24.0, 30.0)}
AMPS = (0.022, 0.10)
WINDOWS = ("vor", "sph", "sph6", "none")
GENS = ("rad", "z", "local")
LAMBDAS = (1.0, 0.0)
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def cfg_of(n, L):
    return B3.base_cfg(s=S_MAIN, g=G_MAIN, n=n, L=float(L), delta=DELTA)


# ================= fields =================
def boost_rigid(cfg, zc, amp):
    X, Y, Z = B3.coords(cfg["n"], cfg["h"])
    Zc = Z - zc
    R = np.sqrt(X * X + Y * Y + Zc * Zc)
    nx, ny, nz = X / R, Y / R, Zc / R
    K = np.zeros(X.shape + (4, 4))
    K[..., 0, 1], K[..., 0, 2], K[..., 0, 3] = nx, ny, nz
    K[..., 1, 0], K[..., 2, 0], K[..., 3, 0] = nx, ny, nz
    K2 = np.zeros_like(K)
    K2[..., 0, 0] = 1.0
    for i, a in enumerate((nx, ny, nz)):
        for j, b in enumerate((nx, ny, nz)):
            K2[..., 1 + i, 1 + j] = a * b
    bl = amp * np.tanh(R / 2.0)
    return (np.eye(4)[None, None, None] + np.sinh(bl)[..., None, None] * K
            + (np.cosh(bl) - 1.0)[..., None, None] * K2)


def centers_of(kind, d):
    return [0.0] if kind == "single" else [+d / 2.0, -d / 2.0]


def build(cfg, kind, d, amp):
    M4 = B3.embed34(PAIR.seed_pair(cfg, kind, d), cfg)
    if amp == 0.0:
        return M4
    Qs = [boost_rigid(cfg, zc, amp) for zc in centers_of(kind, d)]
    Q = Qs[0] if len(Qs) == 1 else Qs[0] @ Qs[1]
    return B3.sym4(np.einsum("...ab,...bc,...dc->...ad", Q, M4, Q))


# ================= generators + windows =================
def local_rot(vhat):
    W = np.zeros(vhat.shape[:-1] + (4, 4))
    n1, n2, n3 = vhat[..., 0], vhat[..., 1], vhat[..., 2]
    W[..., 1, 2], W[..., 1, 3] = -n3, n2
    W[..., 2, 1], W[..., 2, 3] = n3, -n1
    W[..., 3, 1], W[..., 3, 2] = -n2, n1
    return W


def generator(cfg, M, gen, zc):
    X, Y, Z = B3.coords(cfg["n"], cfg["h"])
    if gen == "rad":
        Zc = Z - zc
        R = np.maximum(np.sqrt(X * X + Y * Y + Zc * Zc), 1e-12)
        return local_rot(np.stack([X / R, Y / R, Zc / R], axis=-1))
    if gen == "z":
        Jz = np.zeros((4, 4))
        Jz[1, 2], Jz[2, 1] = -1.0, 1.0
        return np.broadcast_to(Jz, M.shape)
    if gen == "local":
        _, V = np.linalg.eigh(M[..., 1:4, 1:4])
        return local_rot(V[..., :, 2])
    raise ValueError(gen)


def window(cfg, win, zc, d, plane_z=None):
    """w(x) for the core at z = zc; plane_z = the Voronoi plane."""
    X, Y, Z = B3.coords(cfg["n"], cfg["h"])
    h = cfg["h"]
    r = np.sqrt(X * X + Y * Y + (Z - zc) ** 2)
    if win == "none":
        return np.ones_like(r)
    if win == "vor":
        sgn = 1.0 if zc >= plane_z else -1.0
        return 0.5 * (1.0 + np.tanh(2.0 * sgn * (Z - plane_z) / h))
    if win == "sph":
        return 0.5 * (1.0 + np.tanh((d / 2.0 - r) / h))
    if win == "sph6":
        return 0.5 * (1.0 + np.tanh((6.0 - r) / h))
    raise ValueError(win)


def clock(cfg, M, gen, win, zc, d, plane_z):
    G = generator(cfg, M, gen, zc)
    w = window(cfg, win, zc, d, plane_z)
    return w[..., None, None] * (G @ M - M @ G.swapaxes(-1, -2))


# ================= the kinetic bilinear form =================
def jets(M, cfg):
    return [(wt, [B3.d1(M, ax, cfg["h"], br) for ax in range(3)])
            for br, wt in B3.branches(cfg["stencil"])]


def K_bilinear(a, b, jt, hf, lam, h3):
    tot = 0.0
    for wt, A in jt:
        for Ai in A:
            Fa = B3.comm_eta(a, Ai)
            Fb = B3.comm_eta(b, Ai)
            q = 0.0
            if lam != 1.0:
                q = q + (1.0 - lam) * np.sum(B3.inner_eta(Fa, Fb))
            if lam != 0.0:
                q = q + lam * np.sum(np.einsum(
                    "...ab,...bc,...cd,...ad->...", hf, Fa, hf, Fb,
                    optimize=True))
            tot = tot + wt * 4.0 * q
    return float(h3 * tot)


def gate(cfg, M, a):
    """K_lambda(a, a) against the oracles (kin_of, both registries)."""
    jt, hf, h3 = jets(M, cfg), EXT.h_cov_np(M), cfg["h"] ** 3
    p = LAG.default_params(s=S_MAIN, g=G_MAIN, delta=DELTA)
    k0 = K_bilinear(a, a, jt, hf, 0.0, h3)
    k1 = K_bilinear(a, a, jt, hf, 1.0, h3)
    kin = B3.kin_of(M, a, cfg)
    c_i1 = -4.0 * LAG.omega_decompose(LAG.REGISTRY["I1"], M, cfg, p, a)[2]
    c_h = -4.0 * LAG.omega_decompose(EXT.REGISTRY_EXT["I1_h"], M, cfg, p, a)[2]
    rel = lambda x, y: float(abs(x - y) / max(abs(y), 1e-300))  # noqa: E731
    return {"K0": k0, "kin_of": kin, "rel_K0_vs_kin_of": rel(k0, kin),
            "rel_K0_vs_registry_I1": rel(k0, c_i1),
            "K1": k1, "rel_K1_vs_registry_I1_h": rel(k1, c_h),
            "pass": bool(rel(k0, kin) < 1e-10 and rel(k0, c_i1) < 1e-10
                         and rel(k1, c_h) < 1e-10)}


# ================= the measurement =================
def measure(cfg, kind, d, amp):
    """all (window, generator, lambda) reads on one field."""
    M = build(cfg, kind, d, amp)
    jt, hf, h3 = jets(M, cfg), EXT.h_cov_np(M), cfg["h"] ** 3
    cs = centers_of(kind, d)
    out = {}
    for gen in GENS:
        for win in WINDOWS:
            if kind == "single":
                # the window placed about the single core as it sits in the pair
                a1 = clock(cfg, M, gen, win, 0.0, d, plane_z=-d / 2.0)
                afull = clock(cfg, M, gen, "none", 0.0, d, None)
                for lam in LAMBDAS:
                    out[f"{gen}|{win}|{lam:g}"] = {
                        "I0_w": 2.0 * K_bilinear(a1, a1, jt, hf, lam, h3),
                        "I0_full": 2.0 * K_bilinear(afull, afull, jt, hf, lam, h3)}
            else:
                a1 = clock(cfg, M, gen, win, cs[0], d, plane_z=0.0)
                a2 = clock(cfg, M, gen, win, cs[1], d, plane_z=0.0)
                for lam in LAMBDAS:
                    k11 = K_bilinear(a1, a1, jt, hf, lam, h3)
                    k22 = K_bilinear(a2, a2, jt, hf, lam, h3)
                    k12 = K_bilinear(a1, a2, jt, hf, lam, h3)
                    out[f"{gen}|{win}|{lam:g}"] = {
                        "I1": 2.0 * k11, "I2": 2.0 * k22, "C": 2.0 * k12}
    return out


def powerlaw(ds, cs):
    ds, cs = np.asarray(ds, float), np.asarray(cs, float)
    if len(ds) < 3 or np.any(cs == 0) or len(set(np.sign(cs))) != 1:
        return {"fit": None, "reason": "sign change or too few points"}
    x, y = np.log(ds), np.log(np.abs(cs))
    p, a = np.polyfit(x, y, 1)
    yh = a + p * x
    ss = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - float(np.sum((y - yh) ** 2)) / ss if ss > 0 else 0.0
    return {"fit": {"p": float(p), "a": float(a), "R2": float(r2),
                    "sign": int(np.sign(cs[0])), "npts": int(len(ds))}}


def main(quick=False):
    boxes = BOXES[:1] if quick else BOXES
    amps = AMPS[:1] if quick else AMPS
    kinds = ("same", "anti")
    res = {"task": "M5.32 R3 arm (iii) adversarial audit: cross-inertia",
           "toy_point": {"s": S_MAIN, "g": G_MAIN, "delta": DELTA},
           "conventions": __doc__.split("EQUATIONS FIRST")[1].split("Out:")[0],
           "gates": {}, "rows": [], "fits": {}, "python": sys.version}
    # ---- gates on the smallest pair field
    cfg = cfg_of(*BOXES[0])
    Mg = build(cfg, "same", 12.0, AMPS[0])
    ag = clock(cfg, Mg, "rad", "vor", 6.0, 12.0, 0.0)
    res["gates"]["K_vs_oracles_pair_same_d12"] = gate(cfg, Mg, ag)
    log(f"gate {res['gates']['K_vs_oracles_pair_same_d12']}")
    # V4 exactness of the rigid dressing (conjugation invariance)
    eu0, ev0 = B3.e_parts(B3.embed34(PAIR.seed_pair(cfg, "same", 12.0), cfg), cfg)
    eu1, ev1 = B3.e_parts(Mg, cfg)
    res["gates"]["V4_dressed_vs_undressed_d12"] = {
        "V4_undressed": float(ev0), "V4_dressed": float(ev1),
        "E_u_undressed": float(eu0), "E_u_dressed": float(eu1)}
    log(f"V4 gate {res['gates']['V4_dressed_vs_undressed_d12']}")
    # ---- the ladder
    for (n, L) in boxes:
        cfg = cfg_of(n, L)
        for amp in amps:
            singles = {}
            for d in DS_BOX[n]:
                singles[d] = measure(cfg, "single", d, amp)
                log(f"single n={n} amp={amp} d={d} done")
            for kind in kinds:
                for d in DS_BOX[n]:
                    mp = measure(cfg, kind, d, amp)
                    Mtmp = build(cfg, kind, d, amp)
                    m0i = float(np.max(np.sqrt(np.sum(Mtmp[..., 0, 1:] ** 2, -1))))
                    for key, v in mp.items():
                        s = singles[d][key]
                        gen, win, lam = key.split("|")
                        row = {"n": n, "L": L, "h": cfg["h"], "kind": kind,
                               "amp": amp, "max_norm_M0i": m0i, "d": d,
                               "gen": gen, "window": win, "lambda": float(lam),
                               "I1": v["I1"], "I2": v["I2"], "C": v["C"],
                               "I0_w": s["I0_w"], "I0_full": s["I0_full"],
                               "C_over_I0": v["C"] / s["I0_full"],
                               "I1_over_I0": v["I1"] / s["I0_full"],
                               "I1_over_I0w": v["I1"] / s["I0_w"],
                               "I0w_over_I0": s["I0_w"] / s["I0_full"]}
                        res["rows"].append(row)
                    log(f"{kind} n={n} amp={amp} d={d} maxM0i={m0i:.4f} "
                        + " ".join(f"{k.split('|')[0]}/{k.split('|')[1]}:"
                                   f"{v['C'] / singles[d][k]['I0_full']:+.3f}"
                                   for k, v in mp.items() if k.endswith("|1")))
            with open(OUT_JSON, "w") as f:
                json.dump(res, f, indent=1)
    # ---- power-law screen
    fits, best = {}, None
    for kind in kinds:
        for amp in amps:
            for gen in GENS:
                for win in WINDOWS:
                    for lam in LAMBDAS:
                        per_box, ok = {}, True
                        for (n, L) in boxes:
                            rows = [r for r in res["rows"]
                                    if r["kind"] == kind and r["amp"] == amp
                                    and r["gen"] == gen and r["window"] == win
                                    and r["lambda"] == lam and r["n"] == n
                                    and r["d"] >= 12.0]
                            rows.sort(key=lambda r: r["d"])
                            f = powerlaw([r["d"] for r in rows],
                                         [r["C_over_I0"] for r in rows])
                            per_box[f"n{n}"] = f
                            ft = f["fit"]
                            ok = ok and ft is not None and 0.8 <= abs(ft["p"]) <= 1.2 \
                                and ft["R2"] >= 0.95
                        key = f"{kind}|amp{amp:g}|{gen}|{win}|lam{lam:g}"
                        fits[key] = {"boxes": per_box, "clean_tail_both_boxes": bool(ok)}
                        sc = min((b["fit"]["R2"] if b["fit"] else -1.0)
                                 for b in per_box.values())
                        if best is None or sc > best[1]:
                            best = (key, sc)
    res["fits"] = fits
    res["best_min_R2"] = {"key": best[0], "min_R2": best[1],
                          "fits": fits[best[0]]} if best else None
    res["any_clean_tail"] = bool(any(v["clean_tail_both_boxes"] for v in fits.values()))
    res["runtime_s"] = time.time() - T0
    with open(OUT_JSON, "w") as f:
        json.dump(res, f, indent=1)
    log(f"any clean tail: {res['any_clean_tail']}; best min R2: {res['best_min_R2']}")


if __name__ == "__main__":
    main(quick="--quick" in sys.argv)
