#!/usr/bin/env python3
"""m8_1_eigensolve.py

Eigenvalues of -Delta on the parameter rectangle (y, w) in [0, pi*R] x [-W, W],
metric ds^2 = dy^2 + f(y)^2 dw^2, f(y) = cos(y/R), R = 1, area weight |f|.
Twisted seam: psi(pi, -w) = -psi(0, w)  (sign flip + w reversal).
Neumann at w = +-W.  Cone point at y = pi/2 where f = 0.

Separation of variables (spec "useful structure"):
  psi = u(y) Phi_n(w),  Phi_n(w) = cos(n pi (w + W)/(2W)),  mu_n = (n pi/(2W))^2.
  n even (Phi even in w)  -> u anti-periodic across seam: u(pi) = -u(0), u'(pi) = -u'(0)
  n odd  (Phi odd in w)   -> u periodic:                  u(pi) =  u(0), u'(pi) =  u'(0)
Reduced weighted Sturm-Liouville problem (as given in the spec):
  -( |f| u' )' + (mu/|f|) u = lambda |f| u   on (0, pi),  weight |f| = |cos y|.

Cone boundary data (spec):
  mu = 0 channel:  u(delta) = u_N ln(|delta|/(2R)) + u_D + o(1) on each side
                   (delta = y - pi/2; side '-' is delta<0, side '+' is delta>0).
  mu > 0 channels: regular Frobenius branch |delta|^sqrt(mu) on both sides.

Realization (A) REGULAR:  u_N^+ = 0 and u_N^- = 0 in the mu=0 channel (no condition
  links u_D^+ to u_D^-); regular branch in every mu>0 channel.
Realization (B) MATCHED(delta0):  utilde_D^+ = utilde_D^-  AND  u_N^+ + u_N^- = 0,
  with utilde_D = u_D + u_N ln(delta0/(2R)).  Other channels as in (A).

T1: realization (A): lowest 6 eigenvalues per W, first positive eigenvalue,
    P1 weighted FEM (mesh graded at the cone), 4 resolutions + Richardson.
T2: realization (B) at W = 1.0, delta0 in logspace(-3, 1, 15), implemented EXACTLY
    through the log-trace extraction via Frobenius series about the cone
    (variable s = 1 - sin y; both halves map to the same ODE in s, cone at s=0,
    seam end at s=1; series radius of convergence is 2, so s=1 is inside).

Outputs:
  ../data/m8_1_spectrum.json    (T1)
  ../data/m8_1_delta_scan.json  (T2)
"""

import json
import os
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))
R = 1.0
PI = np.pi

# ---------------------------------------------------------------------------
# Section 1: P1 weighted FEM for the reduced 1D Sturm-Liouville problems.
# Weak form:  int |f| u' v' dy + mu int (1/|f|) u v dy = lambda int |f| u v dy.
#  - conforming P1 => energy form finite => log branch (infinite energy) excluded
#    => converges to u_N = 0 on each side in the mu=0 channel  [realization (A)]
#  - mu>0: cone node DOFs removed (Dirichlet) = regular branch |delta|^sqrt(mu)
#    (the singular branch has infinite energy and infinite (mu/|f|) integral).
#  - cone node DOUBLED: no continuity is imposed between the two sides
#    (spec: "no condition links u_D^+ to u_D^-").
#  - seam handled by identifying the node at y=pi with the node at y=0 with a
#    sign: -1 for anti-periodic (n even), +1 for periodic (n odd).
# ---------------------------------------------------------------------------

GX, GW = np.polynomial.legendre.leggauss(6)  # Gauss-Legendre on [-1, 1]


def graded_nodes(nhalf, p=2):
    """Nodes on [0, pi/2], clustered at the cone y=pi/2 (grading exponent p)."""
    t = np.linspace(0.0, 1.0, nhalf + 1)
    return (PI / 2.0) * (1.0 - (1.0 - t) ** p)


def sector_eigs(W, n, nhalf, k=10, p=2):
    """Lowest k eigenvalues of the reduced problem for transverse index n."""
    mu = (n * PI / (2.0 * W)) ** 2
    seam_sign = -1.0 if n % 2 == 0 else 1.0  # even -> anti-periodic, odd -> periodic
    yl = graded_nodes(nhalf, p)              # 0 .. pi/2 (node nhalf = cone, left side)
    yr = (PI - yl)[::-1]                     # pi/2 .. pi (node 0 = cone, right side)
    yn = np.concatenate([yl, yr])            # cone node DOUBLED (indices nL-1, nL)
    nL = len(yl)
    nnode = len(yn)
    e0 = np.concatenate([np.arange(nL - 1), nL + np.arange(nL - 1)])
    e1 = e0 + 1

    ya, yb = yn[e0], yn[e1]
    h = yb - ya
    yg = ya[:, None] + h[:, None] * 0.5 * (1.0 + GX[None, :])   # Gauss points
    wg = h[:, None] * 0.5 * GW[None, :]                          # Gauss weights
    a = np.abs(np.cos(yg))                                       # weight |f|
    ph1 = (yb[:, None] - yg) / h[:, None]                        # hat at left node
    ph2 = (yg - ya[:, None]) / h[:, None]                        # hat at right node

    Ia = np.sum(a * wg, axis=1)
    k11 = Ia / h ** 2                                            # int a phi' phi'
    k12 = -k11
    k22 = k11
    m11 = np.sum(a * ph1 * ph1 * wg, axis=1)                     # int a phi phi
    m12 = np.sum(a * ph1 * ph2 * wg, axis=1)
    m22 = np.sum(a * ph2 * ph2 * wg, axis=1)
    if n >= 1:
        # potential term mu int (1/a) phi_i phi_j ; with the cone DOF removed all
        # remaining basis functions vanish linearly at the cone -> integrable.
        p11 = mu * np.sum(ph1 * ph1 / a * wg, axis=1)
        p12 = mu * np.sum(ph1 * ph2 / a * wg, axis=1)
        p22 = mu * np.sum(ph2 * ph2 / a * wg, axis=1)
        k11, k12, k22 = k11 + p11, k12 + p12, k22 + p22

    # DOF map: node -> (dof index, sign); -1 = removed (Dirichlet at cone, n>=1)
    dof = -np.ones(nnode, dtype=int)
    sgn = np.ones(nnode)
    idx = 0
    for inode in range(nnode - 1):           # last node (y=pi) is the seam duplicate
        if n >= 1 and inode in (nL - 1, nL):  # both cone copies removed
            continue
        dof[inode] = idx
        idx += 1
    dof[nnode - 1] = dof[0]                  # seam identification y=pi ~ y=0
    sgn[nnode - 1] = seam_sign
    ndof = idx

    rows = np.concatenate([e0, e0, e1, e1])
    cols = np.concatenate([e0, e1, e0, e1])
    kdat = np.concatenate([k11, k12, k12, k22])
    mdat = np.concatenate([m11, m12, m12, m22])
    r, c = dof[rows], dof[cols]
    s = sgn[rows] * sgn[cols]
    keep = (r >= 0) & (c >= 0)
    K = sp.coo_matrix((kdat[keep] * s[keep], (r[keep], c[keep])), shape=(ndof, ndof)).tocsc()
    M = sp.coo_matrix((mdat[keep] * s[keep], (r[keep], c[keep])), shape=(ndof, ndof)).tocsc()

    kk = min(k, ndof - 2)
    vals = eigsh(K, k=kk, M=M, sigma=-1.0, which="LM", return_eigenvectors=False)
    return np.sort(vals)


def richardson(vals):
    """Extrapolate a sequence of eigenvalues at h, h/2, h/4, h/8 (4 resolutions).

    Fit lambda_N = lambda_inf + C N^-q on the last 3 points; error estimate =
    |extrapolation(last 3) - extrapolation(first 3)| (fallback: last increment).
    """
    v = np.asarray(vals, float)

    def extrap3(a, b, c):
        d1, d2 = b - a, c - b
        if abs(d2) < 1e-13 or abs(d1) < 1e-13:
            return c, 2.0
        r = d1 / d2
        if r <= 1.0:                          # non-monotone / no convergence signal
            return c + d2, 2.0
        q = np.log2(r)
        return c + d2 / (r - 1.0) * 1.0, q    # lam_inf = c - e_c ; e_c = -d2/(r-1)

    lam_b, q_b = extrap3(v[-3], v[-2], v[-1])
    lam_a, _ = extrap3(v[-4], v[-3], v[-2])
    err = max(abs(lam_b - lam_a), 1e-12)
    return lam_b, err, q_b


# ---------------------------------------------------------------------------
# Section 2: EXACT mu=0 channel via Frobenius series about the cone.
#
# Substitution x = sin y, s = 1 - x maps BOTH halves onto s in (0, 1]
# (cone s=0, seam end s=1) with the SAME ODE:
#   s(2-s) u_ss + 2(1-s) u_s + lambda u = 0.
# Since |delta| = sqrt(2 s) (1 + O(s)), the spec's trace expansion
#   u = u_N ln(|delta|/(2R)) + u_D + o(1)  becomes  u = u_N * (1/2) ln(s/2) + u_D.
# Frobenius solutions (indicial root 0 double):
#   u_reg  = sum a_k s^k,  a_0 = 1,  a_{k+1} = (k(k+1) - lambda) a_k / (2 (k+1)^2)
#   u_log  = (1/2) ln(s/2) u_reg + sum b_k s^k,  b_0 = 0,
#     b_{k+1} = [ (k(k+1)-lambda) b_k - 2(k+1) a_{k+1} + (k+1/2) a_k ] / (2 (k+1)^2)
# so u_reg carries (u_N, u_D) = (0, 1) and u_log carries (u_N, u_D) = (1, 0):
# the log-trace extraction is EXACT (series coefficients), no approximation.
#
# Seam conditions in s (left half: ds/dy = -cos y, at y=0 ds/dy = -1;
# right half: at y=pi ds/dy = +1):
#   u(pi) = -u(0)   ->  uR(1) = -uL(1)
#   u'(pi) = -u'(0) ->  uR_s(1) = uL_s(1)
# With uL = A_L u_log + B_L u_reg, uR = A_R u_log + B_R u_reg and
# v_reg = u_reg(1), v_reg' = du_reg/ds(1), v_log' = du_log/ds(1):
#
# (A): A_L = A_R = 0  =>  (B_R + B_L) v_reg = 0 and (B_R - B_L) v_reg' = 0
#      => eigencondition  v_reg(lambda) = 0  OR  v_reg'(lambda) = 0.
# (B): A_R = -A_L =: -A (from u_N^+ + u_N^- = 0) and B_R - B_L = 2 A L,
#      L := ln(delta0/(2R))  (from utilde_D^+ = utilde_D^-). Seam value gives
#      (B_R + B_L) v_reg = 0 (v_log cancels); seam slope gives
#      A (v_log' - L v_reg') = 0.  => eigencondition
#      v_reg(lambda) = 0  (A = 0 family)   OR
#      F(lambda; delta0) := v_log'(lambda) - L v_reg'(lambda) = 0  (A != 0 family).
# For lambda < 0 all a_k > 0, so v_reg > 0 and v_reg' > 0: negative eigenvalues
# can only come from the F family.
# ---------------------------------------------------------------------------


def mu0_seam_series(lams, s=1.0):
    """Sum the two Frobenius series (and their s-derivatives) at s.

    Returns (v_reg, dv_reg, v_log, dv_log, esc): true value = returned * exp(esc)
    (per-lambda rescaling avoids overflow for deeply negative lambda; the common
    positive factor preserves signs and zeros of any fixed-lambda combination).
    """
    lam = np.atleast_1d(np.asarray(lams, float))
    kmax = int(3.5 * np.sqrt(np.max(np.abs(lam)) * s / 2.0)) + 400
    alpha = np.ones_like(lam)     # a_k s^k   (k = 0)
    beta = np.zeros_like(lam)     # b_k s^k
    Sa = alpha.copy()
    Sda = np.zeros_like(lam)
    Sb = np.zeros_like(lam)
    Sdb = np.zeros_like(lam)
    esc = np.zeros_like(lam)
    kmin = 2.0 * np.sqrt(np.max(np.abs(lam)) * s / 2.0) + 50.0
    for k in range(kmax):
        den = 2.0 * (k + 1) ** 2
        fac = (k * (k + 1) - lam)
        alpha_n = alpha * s * fac / den
        beta_n = (beta * s * fac - 2.0 * (k + 1) * alpha_n + (k + 0.5) * s * alpha) / den
        alpha, beta = alpha_n, beta_n
        kk = k + 1
        Sa += alpha
        Sda += kk * alpha / s
        Sb += beta
        Sdb += kk * beta / s
        big = np.abs(alpha) > 1e250
        if big.any():
            for arr in (alpha, beta, Sa, Sda, Sb, Sdb):
                arr[big] *= 1e-250
            esc[big] += np.log(1e250)
        if k > kmin and k % 50 == 0:
            if np.all(np.abs(alpha) + np.abs(beta) <=
                      1e-30 * (np.abs(Sa) + np.abs(Sb) + 1e-300)):
                break
    lnh = 0.5 * np.log(s / 2.0)
    v_reg = Sa
    dv_reg = Sda
    v_log = lnh * Sa + Sb
    dv_log = Sa / (2.0 * s) + lnh * Sda + Sdb
    return v_reg, dv_reg, v_log, dv_log, esc


def _vreg(lam):
    return mu0_seam_series(lam)[0]


def _dvreg(lam):
    return mu0_seam_series(lam)[1]


def F_matched(lam, L):
    """Matched-family eigencondition F = v_log' - L v_reg' (scaled, sign-exact)."""
    v_reg, dv_reg, v_log, dv_log, esc = mu0_seam_series(lam)
    return dv_log - L * dv_reg


def refine_zeros(fvec_scalar, grid, vals):
    """brentq at every sign change of vals over grid."""
    roots = []
    sgn = np.sign(vals)
    for i in range(len(grid) - 1):
        if sgn[i] == 0:
            roots.append(grid[i])
        elif sgn[i] * sgn[i + 1] < 0:
            roots.append(brentq(fvec_scalar, grid[i], grid[i + 1],
                                xtol=1e-13, rtol=8.9e-16, maxiter=200))
    return roots


# ---------------------------------------------------------------------------
# Drivers
# ---------------------------------------------------------------------------

W_LIST_T1 = [0.3, 0.6, 0.9, 1.2, 1.5, 1.5708, 1.65, 1.8, 2.2, 2.6, 3.0]
W_SUPPORT = [1.0]        # needed by T2 (and by the T3 cross-check)
NHALF_LIST = [128, 256, 512, 1024]   # elements per half; total N = 2*nhalf
MU_CUTOFF = 60.0         # skip sectors with mu_n > 60 (their spectrum >= mu_n)
NMAX_SECTOR = 12


def run_fem_all(W):
    """All sectors at all resolutions; per-sector Richardson for 10 eigenvalues."""
    sectors = [n for n in range(NMAX_SECTOR + 1)
               if (n * PI / (2.0 * W)) ** 2 <= MU_CUTOFF]
    out = {}
    for n in sectors:
        raw = {nh: sector_eigs(W, n, nh) for nh in NHALF_LIST}
        kmin = min(len(v) for v in raw.values())
        entries = []
        for idx in range(kmin):
            vals = [raw[nh][idx] for nh in NHALF_LIST]
            lam, err, q = richardson(vals)
            entries.append({"lambda": lam, "err": err, "q_fit": q,
                            "raw": [float(v) for v in vals]})
        out[n] = {"mu": (n * PI / (2.0 * W)) ** 2, "eigs": entries}
    return out


def run_T1():
    results = {}
    for W in W_LIST_T1 + W_SUPPORT:
        sec = run_fem_all(W)
        merged = []
        for n, d in sec.items():
            for e in d["eigs"]:
                merged.append({"lambda": e["lambda"], "err": e["err"],
                               "q_fit": e["q_fit"], "sector_n": n, "mu": d["mu"],
                               "raw_by_resolution": e["raw"]})
        merged.sort(key=lambda e: e["lambda"])
        lowest6 = merged[:6]
        lam1p = next(e for e in merged if e["lambda"] > 1e-6)
        n_excluded = max(sec.keys()) + 1      # sectors are contiguous 0..n_max
        mu_first_excluded = (n_excluded * PI / (2.0 * W)) ** 2
        results[f"{W:g}"] = {
            "W": W,
            "sectors_included": sorted(sec.keys()),
            "mu_by_sector": {str(n): sec[n]["mu"] for n in sec},
            "resolutions_elements_total": [2 * nh for nh in NHALF_LIST],
            "lowest6": lowest6,
            "lambda1plus": lam1p,
            "mu_first_excluded_sector": mu_first_excluded,
            "sector_coverage_ok": bool(mu_first_excluded > lowest6[-1]["lambda"] + 1.0),
        }
        print(f"W={W:8g}  lowest6=" +
              " ".join(f"{e['lambda']:.6f}(n={e['sector_n']})" for e in lowest6) +
              f"  lam1+={lam1p['lambda']:.8f}+-{lam1p['err']:.1e}")
    # exact mu=0 (A) eigenvalues (W-independent): zeros of v_reg and v_reg'
    grid = np.linspace(1e-6, 45.0, 4501)
    vr = mu0_seam_series(grid)[0]
    dvr = mu0_seam_series(grid)[1]
    zr = refine_zeros(lambda x: float(_vreg(np.array([x]))[0]), grid, vr)
    zdr = refine_zeros(lambda x: float(_dvreg(np.array([x]))[0]), grid, dvr)
    mu0_exact = sorted([0.0] + zr + zdr)   # lambda = 0: u = +-const halves (exact)
    print("mu0 (A) exact eigenvalues:", " ".join(f"{z:.10f}" for z in mu0_exact[:8]))
    return results, mu0_exact


def run_T2(t1_results, mu0_exact):
    W = 1.0
    # minima of the mu>0 sectors at W=1 (FEM + Richardson): candidates for the
    # first positive eigenvalue of realization (B) alongside the mu=0 channel
    sec = run_fem_all(W)
    nge1_entries = []
    for n, dd in sec.items():
        if n >= 1:
            e = dd["eigs"][0]
            nge1_entries.append({"sector_n": n, "mu": dd["mu"],
                                 "lambda": e["lambda"], "err": e["err"]})
    nge1_entries.sort(key=lambda e: e["lambda"])
    nge1_min = nge1_entries[0]

    # positive zeros of v_reg (B family with A=0; delta0-independent)
    grid = np.linspace(1e-6, 45.0, 4501)
    vr = mu0_seam_series(grid)[0]
    vreg_zeros = refine_zeros(lambda x: float(_vreg(np.array([x]))[0]), grid, vr)

    delta0s = np.logspace(-3, 1, 15)
    scan = []
    for d0 in delta0s:
        L = np.log(d0 / (2.0 * R))
        f_scalar = lambda x: float(F_matched(np.array([x]), L)[0])
        # negative-lambda scan of F in log(-lambda)
        lam_hi = max(1e4, 40.0 / d0 ** 2)
        ug = np.linspace(np.log(1e-8), np.log(lam_hi), 2000)
        lam_neg = -np.exp(ug)
        Fn = F_matched(lam_neg, L)
        neg_roots = []
        sgn = np.sign(Fn)
        for i in range(len(ug) - 1):
            if sgn[i] * sgn[i + 1] < 0:
                r = brentq(lambda u: f_scalar(-np.exp(u)), ug[i], ug[i + 1],
                           xtol=1e-14, rtol=8.9e-16, maxiter=200)
                neg_roots.append(-np.exp(r))
        neg_roots.sort()
        # positive-lambda F zeros
        gp = np.linspace(1e-6, 45.0, 4000)
        Fp = F_matched(gp, L)
        pos_F = refine_zeros(f_scalar, gp, Fp)
        # F at lambda = 0
        F0 = float(F_matched(np.array([0.0]), L)[0])
        mu0_pos = sorted([z for z in pos_F + vreg_zeros if z > 1e-6])
        first_pos_mu0 = mu0_pos[0] if mu0_pos else None
        cands = [(first_pos_mu0, "mu0-channel")] if first_pos_mu0 else []
        cands.append((nge1_min["lambda"], f"sector n={nge1_min['sector_n']}"))
        first_pos, src = min(cands, key=lambda t: t[0])
        allroots = neg_roots + mu0_pos
        zeroish = [z for z in allroots if abs(z) < 1e-6]
        scan.append({
            "delta0": float(d0), "L=ln(delta0/2R)": float(L),
            "negative_eigenvalues": neg_roots,
            "n_negative": len(neg_roots),
            "zeroish_eigenvalue_exists": bool(zeroish),
            "F_at_lambda0": F0,
            "first_positive": float(first_pos),
            "first_positive_source": src,
            "mu0_channel_positive_first3": mu0_pos[:3],
        })
        nr = " ".join(f"{z:.6g}" for z in neg_roots) or "-"
        print(f"delta0={d0:10.5g}  neg=[{nr}]  F(0)={F0:+.4f}  "
              f"first_pos={first_pos:.8f} ({src})")
    return {
        "W": W,
        "delta0_values": [float(x) for x in delta0s],
        "method": "exact log-trace Frobenius shooting (series about cone, s=1-sin y)",
        "eigencondition": "v_reg(lam)=0 OR F(lam)=v_log'(lam)-ln(delta0/2R)*v_reg'(lam)=0",
        "vreg_zeros_positive": vreg_zeros,
        "n_ge1_sector_minima_W1_fem": nge1_entries,
        "mu0_A_exact_eigs_for_reference": mu0_exact[:8],
        "per_delta0": scan,
    }


def main():
    os.makedirs(DATA, exist_ok=True)
    print("== T1: realization (A), P1 weighted FEM + Richardson ==")
    t1, mu0_exact = run_T1()
    with open(os.path.join(DATA, "m8_1_spectrum.json"), "w") as f:
        json.dump({"task": "T1 realization A",
                   "method": "P1 weighted FEM, graded mesh (p=2) at cone, "
                             "4 resolutions + Richardson",
                   "W_list": W_LIST_T1, "W_support": W_SUPPORT,
                   "mu0_A_exact_eigs_series": mu0_exact[:8],
                   "per_W": t1}, f, indent=1, default=float)
    print("== T2: realization (B) matched family at W=1.0 ==")
    t2 = run_T2(t1, mu0_exact)
    with open(os.path.join(DATA, "m8_1_delta_scan.json"), "w") as f:
        json.dump(t2, f, indent=1, default=float)
    print("JSONs written to", DATA)


if __name__ == "__main__":
    main()
