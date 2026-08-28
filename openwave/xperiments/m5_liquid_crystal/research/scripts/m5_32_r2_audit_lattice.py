"""M5.32 R2 lattice-arms (b, c) ADVERSARIAL AUDIT: an independent
rebuild of every lattice claim of the candidate

    L_lambda = -4 [(1 - lambda) I1 + lambda I1_h] - V4,
    I1_h = I1 with eta -> h_cov = eta + 2 (eta u)(eta u)^T,
    u = the timelike unit eigenvector of M eta (u^T eta u = -1).

The producers' scripts (m5_32_r2_b_bounded.py, m5_32_r2_c_screen.py)
and their JSONs were NOT read; the oracles are the registries
(m5_32_lagrangian.py: I1, V4, F_of_A; m5_32_terms_ext.py: h_cov_np,
timelike_eig_np), the certified stack m5_21_3_a_4d.py (sym stencil,
e_parts, kin_of, embed34), the family builders m5_21_8_b_lattice.py
(dressed, a0_unit) and the M5.21.11 g-arm record (m5_21_11_garm.json).

EQUATIONS FIRST
---------------
Jets A_mu = d_mu M, curvature F_{mu nu} = A_mu eta A_nu - A_nu eta A_mu.
Per-pair quadratic forms on an antisymmetric F:
    q_eta(F) = <F,F>_eta = tr(eta F eta F^T) = sum_ab eta_a eta_b F_ab^2
    q_h(F)   = tr(h F h F^T),  h = h_cov(M)
    q_lam    = (1 - lam) q_eta + lam q_h          (LINEAR in lam)
Legendre energy density of L_lambda (A_0 = omega a0, every term
quadratic in omega, H = omega dL/domega - L):
    e_lam(x) = 4 [ sum_{i<j} q_lam(F_ij) + omega^2 sum_i q_lam(F_0i) ] + V4
    E_stat(lam) = h^3 sum_x 4 sum_{i<j} q_lam(F_ij) + E_V
    kin(lam)    = h^3 sum_x 4 sum_i q_lam(comm_eta(a0, A_i))
    E(omega)    = E_stat + omega^2 kin
(the certified E_cert = 4 (U + omega^2 T) + V4 is lam = 0: e_parts /
kin_of of m5_21_3_a_4d.py). Sym stencil: 1/2 (fwd + bwd) with the
density evaluated per branch; h_cov taken from M at the cell.

L1 (positivity). In the u-frame (Lambda u = e0, Lambda in SO(1,3))
h_cov = 1, q_eta = S - T, q_h = S + T with S = sum_{a,b >= 1} F_ab^2,
T = 2 sum_b F_0b^2 (both >= 0), so per pair
    q_lam = S - (1 - 2 lam) T  >= 0  for lam >= 1/2.
Both q_eta and q_h are invariant under internal Lorentz maps
F -> Lambda F Lambda^T (h_cov -> Lambda^-T h_cov Lambda^-1), so the
frame statement is pointwise for ANY symmetric M with a real timelike
eigenvector of M eta and ANY jets. The audit tests it on random
(M, jets), and measures the V4 cost of leaving the "u exists" region:
V4 depends on the spectrum only, V4 = w sum_p (tr((M eta)^p) - C_p)^2,
so min over spectra {m + i y, m - i y, r1, r2} (complex pair, y >= 0)
gives the floor of V4 on fields without u.

L2 (fixed-J electron). Rigid family: M(amp) = Qb Mb Qb^T with the
m5_21_8 hedgehog Mb (m = 0, rotation-only), Qb the radial boost with
b(r) = amp tanh(r/2) (the M5.21.16 family), a0 = Qb a0_base Qb^T.
    E(amp; J) = E_stat(amp) + J^2 / (4 kin(amp)),  omega* = J / (2 kin)
    closure: dE*/dJ = omega*(amp*)  (envelope theorem)
L3 (concavity). E_free(omega) = min_amp [E_stat + omega^2 kin] is a
min of affine functions of omega^2, hence concave in omega^2 for ANY
sign of kin; the substantive claim is kin(amp) > 0 for lam >= 1/2.
L4 (g-arm). E(m) of Qb(m) embed34(M3_end) Qb(m)^T, n = 48, L = 48,
m in +-3 artanh(1/g), 121 points, parabolic refinement (the M5.21.11
recipe); gain = E(m*) - E(0), q = slope of log(-gain) vs log artanh(1/g).
L5 (block identity). For M = embed34(M3): u = e0, h_cov = 1, F_ij is
3x3 block, so q_h(F) = q_eta(F) exactly; checked at the density level.

Out: ../data/m5_32_r2_audit_lattice.json
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys
import time

import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


L0 = _load("m5_32_lagrangian", "m5_32_lagrangian.py")
TX = _load("m5_32_terms_ext", "m5_32_terms_ext.py")
B3 = L0.B3
B8 = _load("m5_21_8_b_lattice", "m5_21_8_b_lattice.py")

ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
W1 = float(B3.W1)
G, DELTA, S = 32.0, 0.3, -1.0
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


# ================= densities =================
def q_eta(F):
    return B3.inner_eta(F, F)


def q_h(F, h):
    return np.einsum("...ab,...bc,...cd,...ad->...", h, F, h, F,
                     optimize=True)


def h_of(M):
    """registry h_cov (eigen-solved u)."""
    return TX.h_cov_np(M)


def lattice_parts(M, cfg, a0=None, h=None):
    """(E_u0, E_u1, E_V, kin0, kin1) h^3-weighted; the lam = 0 parts are
    the certified e_parts / kin_of by construction (same stencil)."""
    h3 = cfg["h"] ** 3
    if h is None:
        h = h_of(M)
    eu0 = eu1 = k0 = k1 = 0.0
    for br, (A, wt) in B3.a_fields(M, cfg).items():
        for i in range(3):
            for j in range(i + 1, 3):
                F = B3.comm_eta(A[i], A[j])
                eu0 += wt * 4.0 * np.sum(q_eta(F))
                eu1 += wt * 4.0 * np.sum(q_h(F, h))
            if a0 is not None:
                F = B3.comm_eta(a0, A[i])
                k0 += wt * 4.0 * np.sum(q_eta(F))
                k1 += wt * 4.0 * np.sum(q_h(F, h))
    _, ev = B3.e_parts(M, cfg)
    return (h3 * eu0, h3 * eu1, float(ev), h3 * k0, h3 * k1)


def mix(a, b, lam):
    return (1.0 - lam) * a + lam * b


def qb_field(cfg, bfield):
    """radial boost Qb with per-cell rapidity bfield (n,n,n)."""
    n, hh = cfg["n"], cfg["h"]
    X, Y, Z = B3.coords(n, hh)
    R = np.sqrt(X * X + Y * Y + Z * Z)
    nx, ny, nz = X / R, Y / R, Z / R
    K = np.zeros(X.shape + (4, 4))
    K[..., 0, 1], K[..., 0, 2], K[..., 0, 3] = nx, ny, nz
    K[..., 1, 0], K[..., 2, 0], K[..., 3, 0] = nx, ny, nz
    K2 = np.zeros_like(K)
    K2[..., 0, 0] = 1.0
    for i, a in enumerate((nx, ny, nz)):
        for j, b in enumerate((nx, ny, nz)):
            K2[..., 1 + i, 1 + j] = a * b
    return (np.eye(4)[None, None, None] + np.sinh(bfield)[..., None, None]
            * K + (np.cosh(bfield) - 1.0)[..., None, None] * K2)


def conj(Q, M):
    return np.einsum("...ab,...bc,...dc->...ad", Q, M, Q)


def h_from_Q(Q):
    """analytic h_cov for u = Q e0: h = Q^-T Q^-1 (Q in SO(1,3))."""
    Qi = np.linalg.inv(Q)
    return np.einsum("...ba,...bc->...ac", Qi, Qi)


# ================= L1: positivity on random (M, jets) =================
def boost_to_frame(u):
    """Lambda in SO(1,3) with Lambda u = e0 (u^T eta u = -1, u0 > 0)."""
    v = u[1:] / u[0]
    speed = np.linalg.norm(v)
    if speed < 1e-14:
        return np.eye(4)
    rap = np.arctanh(speed)
    nvec = v / speed
    K = np.zeros((4, 4))
    K[0, 1:] = nvec
    K[1:, 0] = nvec
    return expm(-rap * K)


def stage_L1(rng, n_per=2500, n_jets=4):
    vac = np.diag([-S * G, 1.0, DELTA, 0.0])
    buckets = {"near_0.5": 0.5, "mid_3": 3.0, "far_10": 10.0,
               "random_30": None}
    lams = [0.0, 0.4, 0.5, 0.75, 1.0]
    out = {"buckets": {}, "lams": lams}
    tot_neg = {f"{l:g}": 0 for l in lams}
    tot_dens = 0
    frame_err = 0.0
    inv_err = 0.0
    for nm, sig in buckets.items():
        R = rng.standard_normal((n_per, 4, 4))
        R = 0.5 * (R + R.swapaxes(-1, -2))
        M = vac[None] + sig * R if sig is not None else 30.0 * R
        N = M @ ETA
        lam_c = np.linalg.eigvals(N)
        cplx = np.max(np.abs(lam_c.imag), axis=-1) > 1e-9 * np.maximum(
            np.max(np.abs(lam_c.real), axis=-1), 1.0)
        Mu = M[~cplx]
        # V4 density of the samples (non-Lorentz deformation check)
        C = [(S * G) ** p + 1.0 + DELTA ** p for p in range(1, 5)]
        P = np.broadcast_to(np.eye(4), Mu.shape).copy()
        v4 = np.zeros(len(Mu))
        for p in range(1, 5):
            P = P @ (Mu @ ETA)
            v4 += (np.einsum("...kk->...", P) - C[p - 1]) ** 2
        v4 *= W1
        h = h_of(Mu)
        row = {"n_sampled": int(n_per), "n_complex_spectrum": int(cplx.sum()),
               "n_u_exists": int(len(Mu)),
               "V4_density_min": float(v4.min()), "V4_density_median": float(np.median(v4)),
               "frac_V4_positive": float(np.mean(v4 > 1e-12)),
               "neg_counts": {}, "min_density": {}}
        negs = {f"{l:g}": 0 for l in lams}
        mins = {f"{l:g}": np.inf for l in lams}
        ndens = 0
        for _ in range(n_jets):
            A = rng.standard_normal((4, len(Mu), 4, 4))
            A = 0.5 * (A + A.swapaxes(-1, -2))
            F = L0.F_of_A(A)                      # (N,4,4,4,4)
            # energy density = 4 sum_{mu<nu} q_lam(F_munu), NO eta^mu eta^nu
            qe = np.zeros(len(Mu))
            qh = np.zeros(len(Mu))
            for mu in range(4):
                for nu in range(mu + 1, 4):
                    Fm = F[:, mu, nu]
                    qe += q_eta(Fm)
                    qh += q_h(Fm, h)
            scale = 4.0 * (np.abs(qe) + np.abs(qh)) + 1e-300
            for l in lams:
                d = 4.0 * mix(qe, qh, l)
                negs[f"{l:g}"] += int(np.sum(d < -1e-10 * scale))
                mins[f"{l:g}"] = min(mins[f"{l:g}"], float(np.min(d / scale)))
            ndens += len(Mu)
            # u-frame decomposition check on 50 samples
            for k in range(min(50, len(Mu))):
                u = TX.timelike_eig_np(Mu[k])[0]
                if u[0] < 0:
                    u = -u
                Lb = boost_to_frame(u)
                Ff = np.einsum("ab,mnbc,dc->mnad", Lb, F[k], Lb)
                Sv = Tv = 0.0
                for mu in range(4):
                    for nu in range(mu + 1, 4):
                        Sv += np.sum(Ff[mu, nu, 1:, 1:] ** 2)
                        Tv += 2.0 * np.sum(Ff[mu, nu, 0, 1:] ** 2)
                for l in (0.4, 1.0):
                    pred = 4.0 * (Sv - (1.0 - 2.0 * l) * Tv)
                    d = 4.0 * mix(qe[k], qh[k], l)
                    frame_err = max(frame_err, abs(pred - d) / (abs(d) + abs(pred) + 1e-12))
                # invariance of q_h under a random internal Lorentz map
                Gm = rng.standard_normal((4, 4)) * 0.3
                Gm = ETA @ (Gm - Gm.T)            # eta-antisymmetric generator
                Lr = expm(Gm)
                ML = Lr @ Mu[k] @ Lr.T
                hL = h_of(ML[None])[0]
                FL = np.einsum("ab,mnbc,dc->mnad", Lr, F[k], Lr)
                qhL = sum(q_h(FL[mu, nu], hL) for mu in range(4) for nu in range(mu + 1, 4))
                inv_err = max(inv_err, abs(qhL - qh[k]) / (abs(qh[k]) + 1e-12))
        row["neg_counts"] = negs
        row["min_density_over_scale"] = {k: float(v) for k, v in mins.items()}
        row["n_densities"] = ndens
        for l in lams:
            tot_neg[f"{l:g}"] += negs[f"{l:g}"]
        tot_dens += ndens
        out["buckets"][nm] = row
        log(f"L1 {nm}: complex {cplx.sum()}/{n_per}, negs {negs} of {ndens}")
    out["total_negative_counts"] = tot_neg
    out["total_densities"] = tot_dens
    out["u_frame_decomposition_max_rel_err"] = float(frame_err)
    out["q_h_internal_lorentz_invariance_max_rel_err"] = float(inv_err)

    # V4 floor on fields WITHOUT u: spectrum {m + i y, m - i y, r1, r2}
    C = [(S * G) ** p + 1.0 + DELTA ** p for p in range(1, 5)]

    def v4_spec(x):
        m, y, r1, r2 = x
        z = complex(m, abs(y))
        return W1 * sum((2.0 * (z ** p).real + r1 ** p + r2 ** p - C[p - 1]) ** 2
                        for p in range(1, 5))

    best = None
    starts = []
    for m0 in np.linspace(-35, 35, 15):
        for r1 in (-32.0, 0.0, 1.0, 32.0):
            for r2 in (0.0, 0.3, 1.0):
                for y0 in (0.0, 1.0, 10.0):
                    starts.append([m0, y0, r1, r2])
    for x0 in starts:
        r = minimize(v4_spec, x0, method="Nelder-Mead",
                     options={"xatol": 1e-8, "fatol": 1e-10, "maxiter": 4000})
        if best is None or r.fun < best.fun:
            best = r
    # boundary (double real eigenvalue, y = 0) floor for comparison
    best0 = None
    for x0 in starts:
        r = minimize(lambda x: v4_spec([x[0], 0.0, x[1], x[2]]), [x0[0], x0[2], x0[3]],
                     method="Nelder-Mead", options={"xatol": 1e-8, "fatol": 1e-10, "maxiter": 4000})
        if best0 is None or r.fun < best0.fun:
            best0 = r
    # the WALL: the vacuum sector carries u on the sg eigenvalue (-32);
    # leaving the sector needs that eigenvalue to merge (double real or
    # complex pair) with another one: min V4 over pairs with mean m <= -8
    bestw = None
    for x0 in starts:
        r = minimize(lambda x: v4_spec([-8.0 - x[0] ** 2, x[1], x[2], x[3]]),
                     [np.sqrt(max(-8.0 - x0[0], 0.1)), x0[1], x0[2], x0[3]],
                     method="Nelder-Mead", options={"xatol": 1e-8, "fatol": 1e-10, "maxiter": 4000})
        if bestw is None or r.fun < bestw.fun:
            bestw = r
    # the ORDERED sector wall: in the vacuum sector the timelike eigenvalue
    # is the SMALLEST (-32); real eigenvalues cannot swap order without
    # meeting, so u can be lost only when the two smallest collide:
    # min V4 over {m, m, r1 >= m, r2 >= m}
    besto = None
    for m0 in np.linspace(-40, 5, 46):
        for d1 in (0.1, 1.0, 5.0, 10.0, 30.0):
            for d2 in (0.1, 1.0, 5.0, 10.0, 30.0):
                r = minimize(lambda x: v4_spec([x[0], 0.0, x[0] + x[1] ** 2, x[0] + x[2] ** 2]),
                             [m0, np.sqrt(d1), np.sqrt(d2)], method="Nelder-Mead",
                             options={"xatol": 1e-9, "fatol": 1e-12, "maxiter": 6000})
                if besto is None or r.fun < besto.fun:
                    besto = r
    mo = besto.x[0]
    m, y, r1, r2 = best.x
    # explicit 4x4 realization of the cheap no-u spectrum (y small > 0):
    # (t,x) block [[a,b],[b,c]] has M eta eigenvalues (c-a)/2 +- i sqrt(b^2 - ((a+c)/2)^2)
    yy = 0.5
    a, c = -m, m                      # (c - a)/2 = m
    b = np.sqrt(yy ** 2 + ((a + c) / 2) ** 2)
    Mreal = np.zeros((4, 4)); Mreal[0, 0], Mreal[1, 1], Mreal[0, 1], Mreal[1, 0] = a, c, b, b
    Mreal[2, 2], Mreal[3, 3] = r1, r2
    ev = np.linalg.eigvals(Mreal @ ETA)
    out["V4_floor_no_u"] = {
        "ordered_sector_wall_min": {"V4_density": float(besto.fun),
                                    "spectrum": [float(mo), float(mo), float(mo + besto.x[1] ** 2), float(mo + besto.x[2] ** 2)],
                                    "note": "the vacuum-sector wall: the timelike eigenvalue is the SMALLEST (-32) and can lose u only by colliding with the second smallest; min V4 over {m, m, r1 >= m, r2 >= m}"},
        "sg_merge_wall_min": {"V4_density": float(bestw.fun),
                              "spectrum": [-8.0 - bestw.x[0] ** 2, float(abs(bestw.x[1])), float(bestw.x[2]), float(bestw.x[3])],
                              "note": "min V4 with the sg-scale eigenvalue merged into a real-double/complex pair (mean <= -8)"},
        "explicit_no_u_field": {"M": Mreal.tolist(), "M_eta_eigenvalues": [[float(z.real), float(z.imag)] for z in ev],
                                "V4_density": float(v4_spec([m, yy, r1, r2])), "y": yy},
        "complex_pair_min": {"V4_density": float(best.fun), "m": float(m), "y": float(abs(y)),
                             "r1": float(r1), "r2": float(r2)},
        "double_real_locus_min": {"V4_density": float(best0.fun), "x": best0.x.tolist()},
        "h3_n32_L48": float(1.5 ** 3),
        "V4_per_cell_h3_weighted_n32_L48": float(best.fun * 1.5 ** 3),
        "note": "V4 depends on the spectrum of M eta only; the no-u region "
                "(complex pair) has this V4 floor per cell (density, W1 included)"}
    log(f"L1 V4 floor no-u: {best.fun:.4e} (double-real locus {best0.fun:.4e})")
    return out


# ================= L2 / L3: the rigid dressing family =================
class Family:
    def __init__(self, n, L):
        self.cfg = B3.base_cfg(s=S, g=G, n=n, L=L, delta=DELTA)
        self.Mb = B8.dressed(self.cfg, 0.0)
        self.a0 = B8.a0_unit(self.cfg, 0.0)
        X, Y, Z = B3.coords(n, self.cfg["h"])
        self.R = np.sqrt(X * X + Y * Y + Z * Z)
        self.cache = {}

    def parts(self, amp, use_registry_u=False):
        key = (round(amp, 12), use_registry_u)
        if key in self.cache:
            return self.cache[key]
        Qb = qb_field(self.cfg, amp * np.tanh(self.R / 2.0))
        Md = B3.sym4(conj(Qb, self.Mb))
        a0d = B3.sym4(conj(Qb, self.a0))
        h = h_of(Md) if use_registry_u else h_from_Q(Qb)
        r = lattice_parts(Md, self.cfg, a0d, h=h)
        self.cache[key] = r
        return r


def fixed_j(fam, amps, J, lam):
    rows = []
    for a in amps:
        eu0, eu1, ev, k0, k1 = fam.parts(a)
        es = mix(eu0, eu1, lam) + ev
        k = mix(k0, k1, lam)
        rows.append((a, es, k, es + J * J / (4.0 * k) if k > 0 else np.inf))
    rows = np.array(rows)
    i = int(np.argmin(rows[:, 3]))
    a_star, E_star = rows[i, 0], rows[i, 3]
    interior = 0 < i < len(rows) - 1
    if interior:
        a, b, c = rows[i - 1, 3], rows[i, 3], rows[i + 1, 3]
        da = rows[i + 1, 0] - rows[i, 0]
        if abs(c - 2 * b + a) > 1e-14 and abs(rows[i, 0] - rows[i - 1, 0] - da) < 1e-12:
            a_star = rows[i, 0] - 0.5 * da * (c - a) / (c - 2 * b + a)
            E_star = b - 0.125 * (c - a) ** 2 / (c - 2 * b + a)
    k_star = np.interp(a_star, rows[:, 0], rows[:, 2])
    es_star = np.interp(a_star, rows[:, 0], rows[:, 1])
    return {"J": J, "lam": lam, "amp_star": float(a_star), "E_total": float(E_star),
            "omega_star": float(J / (2.0 * k_star)), "kin_star": float(k_star),
            "E_stat_star": float(es_star), "interior": bool(interior),
            "grid_index": i, "n_amps": len(rows)}


def stage_L2(fam, tag):
    amps_fine = np.round(np.arange(0.0, 0.0601, 0.001), 6)
    amps_coarse = np.round(np.arange(0.0, 0.0601, 0.004), 6)
    amps_tail = np.array([0.08, 0.1, 0.15, 0.2, 0.3, 0.4])
    log(f"L2 {tag}: evaluating {len(amps_fine) + len(amps_tail)} amps")
    for a in np.concatenate([amps_fine, amps_tail]):
        fam.parts(a)
    # registry-u vs analytic-u check at one dressed amp
    p_an = fam.parts(0.02)
    p_reg = fam.parts(0.02, use_registry_u=True)
    ucheck = {"amp": 0.02, "rel_dev_Eu1": float(abs(p_an[1] - p_reg[1]) / abs(p_reg[1])),
              "rel_dev_kin1": float(abs(p_an[4] - p_reg[4]) / abs(p_reg[4]))}
    # certified cross-check at amp = 0 (lam = 0 parts == e_parts / kin_of)
    eu, ev = B3.e_parts(fam.Mb, fam.cfg)
    kk = B3.kin_of(fam.Mb, fam.a0, fam.cfg)
    p0 = fam.parts(0.0)
    cert = {"E_u_cert": float(eu), "E_u_audit_lam0": p0[0], "E_u_audit_lam1": p0[1],
            "kin_cert": float(kk), "kin_audit_lam0": p0[3], "kin_audit_lam1": p0[4],
            "E_V_base": p0[2]}
    ladder = {a: fam.parts(a) for a in np.concatenate([amps_fine, amps_tail])}
    table = {"amp": [], "E_stat_lam1": [], "kin_lam1": [], "E_stat_lam0": [], "kin_lam0": [], "E_V": []}
    for a, p in ladder.items():
        table["amp"].append(float(a)); table["E_stat_lam0"].append(p[0] + p[2])
        table["E_stat_lam1"].append(p[1] + p[2]); table["kin_lam0"].append(p[3])
        table["kin_lam1"].append(p[4]); table["E_V"].append(p[2])
    allamps = np.concatenate([amps_fine, amps_tail])
    res = {"tag": tag, "n": fam.cfg["n"], "L": fam.cfg["L"], "h": fam.cfg["h"],
           "stencil": fam.cfg["stencil"], "u_check": ucheck, "certified_crosscheck": cert,
           "ladder": table, "byJ": {}}
    for J in (50.0, 200.0, 800.0):
        fine = fixed_j(fam, allamps, J, 1.0)
        coarse = fixed_j(fam, np.concatenate([amps_coarse, amps_tail]), J, 1.0)
        # closure: dE*/dJ by central difference with re-minimization
        dJ = 0.01 * J
        Ep = fixed_j(fam, allamps, J + dJ, 1.0)["E_total"]
        Em = fixed_j(fam, allamps, J - dJ, 1.0)["E_total"]
        dEdJ = (Ep - Em) / (2 * dJ)
        fine["dEdJ_numeric"] = float(dEdJ)
        fine["closure_abs"] = float(abs(dEdJ - fine["omega_star"]))
        fine["coarse_step_0.004"] = coarse
        fine["E_dressing_gain"] = float(fine["E_stat_star"] - (p0[1] + p0[2]))
        # lam = 0 control at the same J (guard-free eta: kin can be negative)
        fine["lam0_control"] = fixed_j(fam, allamps, J, 0.0)
        res["byJ"][f"J_{J:g}"] = fine
        log(f"L2 {tag} J={J:g}: amp* {fine['amp_star']:.4f} E {fine['E_total']:.3f} "
            f"omega* {fine['omega_star']:.4f} closure {fine['closure_abs']:.2e}")
    return res


def stage_L3(fam):
    amps = np.round(np.concatenate([np.arange(0.0, 0.1, 0.005), np.arange(0.1, 3.01, 0.1)]), 6)
    out = {"amps": amps.tolist(), "kin": {}, "E_stat": {}, "neg_counts": {},
           "first_negative_amp": {}, "min_pointwise_density": {}}
    parts = [fam.parts(a) for a in amps]
    for lam in (0.0, 0.5, 0.75, 1.0):
        k = np.array([mix(p[3], p[4], lam) for p in parts])
        es = np.array([mix(p[0], p[1], lam) + p[2] for p in parts])
        out["kin"][f"{lam:g}"] = k.tolist()
        out["E_stat"][f"{lam:g}"] = es.tolist()
        out["neg_counts"][f"{lam:g}"] = int(np.sum(k < 0))
        neg = np.where(k < 0)[0]
        out["first_negative_amp"][f"{lam:g}"] = float(amps[neg[0]]) if len(neg) else None
        log(f"L3 lam={lam:g}: kin negatives {out['neg_counts'][f'{lam:g}']}/{len(amps)}, "
            f"first {out['first_negative_amp'][f'{lam:g}']}")
    # pointwise density minimum on the dressed fields (L1 on real fields)
    for amp in (0.02, 0.5, 3.0):
        Qb = qb_field(fam.cfg, amp * np.tanh(fam.R / 2.0))
        Md = B3.sym4(conj(Qb, fam.Mb))
        a0d = B3.sym4(conj(Qb, fam.a0))
        h = h_of(Md)
        dmin = {}
        for lam in (0.4, 0.5, 1.0):
            dm = np.inf
            for br, (A, wt) in B3.a_fields(Md, fam.cfg).items():
                d = 0.0
                for i in range(3):
                    for j in range(i + 1, 3):
                        F = B3.comm_eta(A[i], A[j])
                        d = d + 4.0 * mix(q_eta(F), q_h(F, h), lam)
                    F = B3.comm_eta(a0d, A[i])
                    d = d + 4.0 * mix(q_eta(F), q_h(F, h), lam)   # omega = 1
                dm = min(dm, float(np.min(d)))
            dmin[f"{lam:g}"] = dm
        out["min_pointwise_density"][f"amp_{amp:g}"] = dmin
    # concavity: E_free(omega^2) sampled, second differences
    om2 = np.linspace(0.0, 4.0, 21)
    conc = {}
    for lam in (0.5, 1.0):
        es = np.array(out["E_stat"][f"{lam:g}"]); k = np.array(out["kin"][f"{lam:g}"])
        Ef = np.array([np.min(es + w * k) for w in om2])
        d2 = np.diff(Ef, 2)
        conc[f"{lam:g}"] = {"E_free": Ef.tolist(), "max_second_diff": float(d2.max()),
                            "concave": bool(d2.max() <= 1e-9 * np.max(np.abs(Ef)))}
    out["om2"] = om2.tolist()
    out["concavity"] = conc
    return out


# ================= L4: the g-arm =================
def stage_L4():
    rec = json.load(open(os.path.join(DATA, "m5_21_11_garm.json")))
    rec_g = {(a["branch"], a["g"]): a for a in rec["arms"] if a["s"] == -1}
    out = {"n": 48, "L": 48.0, "h": 1.0, "stencil": "sym", "arms": []}
    for br in ("C", "A"):
        Z = np.load(os.path.join(DATA, f"m5_21_11_end_t11lad_{br}_n48_d0.3.npz"))
        M3 = Z["M"].astype(np.float64)
        for g in (8.0, 32.0):
            cfg = B3.base_cfg(s=S, g=g, n=48, L=48.0, delta=DELTA)
            M4 = B3.embed34(M3, cfg)
            m_his = float(np.arctanh(1.0 / g))
            ms = np.linspace(-3 * m_his, 3 * m_his, 121)
            E0s, E1s = [], []
            ucheck = []
            for k, m in enumerate(ms):
                Qb = qb_field(cfg, np.full(M3.shape[:3], m))
                Md = B3.sym4(conj(Qb, M4))
                h = h_from_Q(Qb)
                if k in (0, 60, 120):
                    hr = h_of(Md)
                    ucheck.append(float(np.max(np.abs(hr - h))))
                eu0, eu1, ev, _, _ = lattice_parts(Md, cfg, None, h=h)
                E0s.append(eu0 + ev); E1s.append(eu1 + ev)
            E0s, E1s = np.array(E0s), np.array(E1s)
            row = {"branch": br, "g": g, "m_his": m_his, "h_analytic_vs_registry_maxabs": max(ucheck),
                   "by_lam": {}}
            for lam in (0.0, 0.5, 0.75, 1.0):
                Es = mix(E0s, E1s, lam)
                i = int(np.argmin(Es)); ms_, Es_ = float(ms[i]), float(Es[i])
                if 0 < i < 120:
                    a, b, c = Es[i - 1], Es[i], Es[i + 1]
                    dm = ms[1] - ms[0]
                    if c - 2 * b + a > 0:
                        ms_ = float(ms[i] - 0.5 * dm * (c - a) / (c - 2 * b + a))
                        Es_ = float(b - 0.125 * (c - a) ** 2 / (c - 2 * b + a))
                E0 = float(Es[60])
                row["by_lam"][f"{lam:g}"] = {"m_star": ms_, "E0": E0, "E_star": Es_,
                                             "gain": Es_ - E0, "argmin_index": i}
            r = rec_g[(br, g)]
            row["record_lam0"] = {"m_star": r["m_star"], "gain": r["gain"], "E0": r["E0"]}
            row["match_record_gain_absdev"] = float(abs(row["by_lam"]["0"]["gain"] - r["gain"]))
            row["curve_lam0"] = [[float(m), float(e)] for m, e in zip(ms[::10], E0s[::10])]
            row["curve_lam1"] = [[float(m), float(e)] for m, e in zip(ms[::10], E1s[::10])]
            out["arms"].append(row)
            log(f"L4 {br} g={g:g}: lam0 gain {row['by_lam']['0']['gain']:.4f} (record {r['gain']:.4f}), "
                f"lam1 gain {row['by_lam']['1']['gain']:.4f} m* {row['by_lam']['1']['m_star']:.4f}")
    # q per branch and lam from the two-point slope g = 8 -> 32
    qs = {}
    for br in ("C", "A"):
        rows = {a["g"]: a for a in out["arms"] if a["branch"] == br}
        qs[br] = {}
        for lam in ("0", "0.5", "0.75", "1"):
            g8, g32 = rows[8.0]["by_lam"][lam]["gain"], rows[32.0]["by_lam"][lam]["gain"]
            x = np.log([np.arctanh(1 / 8), np.arctanh(1 / 32)])
            if g8 < 0 and g32 < 0:
                y = np.log([-g8, -g32])
                qs[br][lam] = float((y[1] - y[0]) / (x[1] - x[0]))
            else:
                qs[br][lam] = None
    out["q_two_point"] = qs
    return out


# ================= L5: the block identity =================
def density_gap(M4, cfg, a0=None):
    h = h_of(M4)
    gap_u = gap_k = 0.0
    scale_u = 0.0
    for br, (A, wt) in B3.a_fields(M4, cfg).items():
        for i in range(3):
            for j in range(i + 1, 3):
                F = B3.comm_eta(A[i], A[j])
                de, dh = q_eta(F), q_h(F, h)
                gap_u = max(gap_u, float(np.max(np.abs(de - dh))))
                scale_u = max(scale_u, float(np.max(np.abs(de))))
            if a0 is not None:
                F = B3.comm_eta(a0, A[i])
                gap_k = max(gap_k, float(np.max(np.abs(q_eta(F) - q_h(F, h)))))
    umax = float(np.max(np.abs(TX.timelike_eig_np(M4)[0][..., 1:])))
    return gap_u, gap_k, scale_u, umax


def stage_L5(rng):
    out = {}
    cfg = B3.base_cfg(s=S, g=8.0, n=32, L=48.0, delta=DELTA)
    M3 = np.load(os.path.join(DATA, L0.STORED3_NPZ))["M"].astype(np.float64)
    M4 = B3.embed34(M3, cfg)
    a0 = B3.gen_catalog(cfg, M4)["rot_z"]
    gu, gk, su, um = density_gap(M4, cfg, a0)
    out["end_state_A_n32_d0.3"] = {"max_abs_density_gap_u": gu, "max_abs_density_gap_kin": gk,
                                   "max_abs_density_u": su, "u_spatial_component_max": um,
                                   "n": 32, "L": 48.0, "g": 8.0}
    from scipy.ndimage import gaussian_filter
    cfg = B3.base_cfg(s=S, g=32.0, n=12, L=18.0, delta=DELTA)
    gaps = []
    for k in range(20):
        M3 = np.stack([[gaussian_filter(rng.normal(size=(12,) * 3), 1.5) for _ in range(3)]
                       for _ in range(3)], axis=-1).reshape(12, 12, 12, 3, 3)
        M3 = np.diag([1.0, DELTA, 0.0])[None, None, None] + 3.0 * 0.5 * (M3 + M3.swapaxes(-1, -2))
        M4 = B3.embed34(M3, cfg)
        a0 = B3.gen_catalog(cfg, M4)["clock_local"]
        gu, gk, su, um = density_gap(M4, cfg, a0)
        gaps.append({"gap_u": gu, "gap_kin": gk, "scale_u": su, "u_spatial": um})
    out["random_block_diag_20"] = {"max_gap_u": max(g["gap_u"] for g in gaps),
                                   "max_gap_kin": max(g["gap_kin"] for g in gaps),
                                   "max_scale_u": max(g["scale_u"] for g in gaps),
                                   "max_u_spatial": max(g["u_spatial"] for g in gaps),
                                   "n": 12, "L": 18.0, "rows": gaps}
    log(f"L5: end-state gap {out['end_state_A_n32_d0.3']['max_abs_density_gap_u']:.2e}, "
        f"random max gap {out['random_block_diag_20']['max_gap_u']:.2e}")
    return out


def main():
    rng = np.random.default_rng(3232)
    out = {"candidate": "L_lambda = -4[(1-lambda) I1 + lambda I1_h] - V4, h_cov = eta + 2(eta u)(eta u)^T",
           "toy": {"g": G, "delta": DELTA, "s": S}, "stencil": "sym (1/2 fwd + bwd, density per branch)",
           "conventions": {
               "energy_density": "4 [sum_{i<j} q_lam(F_ij) + omega^2 sum_i q_lam(F_0i)] + V4, q_lam = (1-lam) tr(eta F eta F^T) + lam tr(h F h F^T)",
               "h_cov": "registry m5_32_terms_ext.h_cov_np (eigen-solved u) for L1/L3-pointwise/L5 and checks; analytic h = Q^-T Q^-1 for the boost-dressed families in L2/L3/L4 (verified against the registry per stage)",
               "lambda_linearity": "E_lam = (1-lam) E_0 + lam E_1 exactly; every lam read from the (E_0, E_1) pair",
               "family_L2": "M(amp) = Qb Mb Qb^T, Mb = m5_21_8 hedgehog (m = 0), b(r) = amp tanh(r/2), a0 = Qb a0_unit Qb^T",
               "fixed_J": "E = E_stat + J^2/(4 kin), omega* = J/(2 kin(amp*)), parabolic refinement on the amp grid, closure by central difference in J with re-minimization"}}
    out["L1"] = stage_L1(rng)
    fam32 = Family(32, 48.0)
    out["L2"] = {"n32_L48": stage_L2(fam32, "n32_L48")}
    out["L3"] = stage_L3(fam32)
    del fam32
    out["L2"]["n48_L48"] = stage_L2(Family(48, 48.0), "n48_L48")
    out["L2"]["n48_L72"] = stage_L2(Family(48, 72.0), "n48_L72")
    out["L2"]["n64_L96"] = stage_L2(Family(64, 96.0), "n64_L96")
    out["L4"] = stage_L4()
    out["L5"] = stage_L5(rng)
    out["runtime_s"] = round(time.time() - T0, 1)
    with open(os.path.join(DATA, "m5_32_r2_audit_lattice.json"), "w") as f:
        json.dump(out, f, indent=1)
    log("done")


if __name__ == "__main__":
    main()
