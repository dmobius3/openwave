#!/usr/bin/env python3
"""m8_1_audit_eigensolve.py

ADVERSARIAL AUDIT of m8_1_eigensolve.py / m8_1_eigensolve_xcheck.py.
Independent recomputation with methods sharing NO code and NO discretization
family with the audited scripts (they used P1 weighted FEM + a Frobenius-series
driver; this audit uses Chebyshev spectral collocation and Riccati shooting).

Problem (same spec): rectangle (y, w) in [0, pi] x [-W, W], metric
ds^2 = dy^2 + f^2 dw^2, f = cos y, weight |f|; twisted seam
psi(pi, -w) = -psi(0, w) (smooth continuation); Neumann at w = +-W; cone at
y = pi/2. Separated sector ODE: -(|f| u')' + (mu/|f|) u = lambda |f| u,
mu_n = (n pi / 2W)^2; n even -> anti-periodic u, n odd -> periodic u.

AUDIT METHOD 1 (realization A) - Chebyshev collocation, singular weight
handled EXPLICITLY:
  substitute t = sin y per half (both halves map to t in (0,1), cone t=1,
  seam t=0; weight |f| dy = dt). The sector ODE becomes the Legendre-type form
    -((1-t^2) u')' + (mu/(1-t^2)) u = lambda u .
  Factor the exact cone behavior u = (1-t^2)^{m/2} v (m = sqrt(mu)); v obeys
    -(1-t^2) v'' + 2(m+1) t v' + m(m+1) v = lambda v ,
  whose regular branch is smooth at t=1 (the singular branch (1-t^2)^{-m},
  resp. the log branch at m=0, is not polynomial-representable, so polynomial
  collocation selects EXACTLY the regular / u_N=0 realization). At t=1 the
  equation row degenerates to 2(m+1) v'(1) + m(m+1) v(1) = lambda v(1), the
  standard regularity row. Two-half coupled system with seam rows
    v2(0) - p v1(0) = 0,   v2'(0) + p v1'(0) = 0   (p=-1 anti-per., +1 per.)
  [u'(y=0) = +u_t(0) on half 1, u'(y=pi) = -u_t(0) on half 2].
  Generalized dense eigenproblem, two resolutions N=64/96 cross-matched.

CLOSED FORM (audit discovery, used as a cross-check only): with the
substitution above the regular solutions are Gegenbauer C_j^{(m+1/2)}(t),
so every sector-(A) eigenvalue is lambda = (m+j)(m+j+1), j = 0, 1, 2, ...
(j even <-> v'(0)=0 branch, j odd <-> v(0)=0 branch; the seam parity p only
fixes the relative sign of the two halves, not the eigenvalue).

AUDIT METHOD 2 (realization B, mu=0 channel) - my own Frobenius connection
(independently derived; different normalization from the audited script) plus
a Riccati log-derivative shooting cross-check:
  x = 1 - t = 1 - sin y (cone x=0, seam x=1):  (x(2-x) u')' + lambda u = 0.
  Basis at the cone: phi = sum a_k x^k (a_0=1), psi = phi ln x + h
  (h = sum b_k x^k, b_0 = 0), with
    a_{k+1} = a_k (k(k+1)-lambda) / (2(k+1)^2)
    b_{k+1} = ((1+2k) a_k - 4(k+1) a_{k+1} + (k(k+1)-lambda) b_k) / (2(k+1)^2).
  Trace dictionary (|delta|: 1-t = 1-cos(delta) => ln x = 2 ln|delta| - ln 2
  + O(delta^2)): u ~ a ln x + b  =>  u_N = 2a, u_D = b + a ln 2.
  With u_i = alpha_i phi + beta_i psi (i=1 left half '-', i=2 right half '+'),
  seam (anti-periodic): u2(0)=-u1(0), u2_t(0)=+u1_t(0), the (B) conditions
  u_N^+ + u_N^- = 0 and utilde_D^+ = utilde_D^- (utilde_D = u_D + u_N L,
  L = ln(delta0/2)) reduce to the determinant condition:
    branch 1 (delta0-independent):  phi(x=1; lambda) = 0
    branch 2 (matched family):      psi_x(1) = (2 ln delta0 - ln 2) phi_x(1).
  Branch-2 eigenfunctions are odd through the cone; equivalently they solve
  the HALF problem u_x(x=1) = 0 (seam) with cone data u ~ u_N ln(|delta|/delta0),
  i.e. b/a = ln 2 - 2 ln delta0. That half form is used for the independent
  Riccati shooting check: integrate V = x u_x / u from the seam (V(1)=0),
    dV/dxi = (-lambda x - (2-x) V^2 + x V)/(2-x),  xi = ln x,
  down to a small x_m, then extract b/a from V against the local (phi, psi)
  basis and solve b/a = ln 2 - 2 ln delta0 for lambda.
  Hand-derived anchors (verified numerically below):
    psi_x(1; lambda=2) = 2 - ln 2, phi_x(1; 2) = -1  =>  delta0* = 2/e;
    kappa->infinity Bessel matching  =>  lambda_neg ~ -4 e^{-2 gamma}/delta0^2.

Outputs: ../data/m8_1_audit.json ; comparison tables printed.
Standalone: needs only numpy + scipy (reads the audited JSONs if present,
skips those columns otherwise).
"""

import json
import math
import os

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import eig as dense_eig
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))
PI = math.pi
LN2 = math.log(2.0)
EULER_GAMMA = 0.5772156649015328606


# ---------------------------------------------------------------------------
# Realization (A): two-half Chebyshev collocation
# ---------------------------------------------------------------------------

def cheb(N):
    """Chebyshev differentiation matrix on x_j = cos(pi j / N) (Trefethen)."""
    x = np.cos(PI * np.arange(N + 1) / N)
    c = np.ones(N + 1)
    c[0] = 2.0
    c[N] = 2.0
    c *= (-1.0) ** np.arange(N + 1)
    X = np.tile(x, (N + 1, 1)).T
    dX = X - X.T
    D = np.outer(c, 1.0 / c) / (dX + np.eye(N + 1))
    D -= np.diag(D.sum(axis=1))
    return D, x


def sector_eigs_colloc(m, p, N, lam_cap):
    """Sector eigenvalues (one N).  m = sqrt(mu), p = seam parity."""
    D, xg = cheb(N)
    t = (1.0 - xg) / 2.0            # t[0]=0 (seam), t[N]=1 (cone)
    Dt = -2.0 * D
    D2 = Dt @ Dt
    n1 = N + 1
    size = 2 * n1
    A = np.zeros((size, size))
    B = np.zeros((size, size))
    for o in (0, n1):
        for j in range(1, n1):      # j=0 slots reserved for seam rows
            row = o + j
            tj = t[j]
            A[row, o:o + n1] = -(1.0 - tj ** 2) * D2[j, :] + 2.0 * (m + 1.0) * tj * Dt[j, :]
            A[row, row] += m * (m + 1.0)
            B[row, row] = 1.0       # at j=N this is the degenerate regularity row
    A[0, n1] = 1.0                  # v2(0) - p v1(0) = 0
    A[0, 0] = -p
    A[n1, n1:2 * n1] = Dt[0, :]     # v2'(0) + p v1'(0) = 0
    A[n1, 0:n1] = p * Dt[0, :]
    w = dense_eig(A, B, right=False)
    w = w[np.isfinite(w)]
    w = w[np.abs(w.imag) < 1e-6 * (1.0 + np.abs(w.real))].real
    w = w[(w > -5.0) & (w < lam_cap + 5.0)]
    return np.sort(w)


def sector_eigs(m, p, lam_cap):
    """Cross-matched eigenvalues from N=64 and N=96."""
    e1 = sector_eigs_colloc(m, p, 64, lam_cap)
    e2 = sector_eigs_colloc(m, p, 96, lam_cap)
    out = []
    for lam in e2:
        if len(e1) and np.min(np.abs(e1 - lam)) < 1e-7 * (1.0 + abs(lam)):
            if not out or lam - out[-1] > 1e-7 * (1.0 + abs(lam)):
                out.append(float(lam))
    return out


def spectrum_A(W, lam_cap=36.0):
    """Full realization-(A) spectrum below lam_cap, with sector labels."""
    ents = []
    n = 0
    min_seen = np.inf
    while True:
        m = n * PI / (2.0 * W)
        if m * m > lam_cap:         # sector minimum >= mu = m^2 (since 1/|f|>=|f|)
            break
        p = -1.0 if n % 2 == 0 else 1.0
        eigs = sector_eigs(m, p, lam_cap)
        if eigs:
            min_seen = min(min_seen, eigs[0])
        for lam in eigs:
            j_near = round((-(2 * m + 1) + math.sqrt((2 * m + 1) ** 2
                                                     + 4 * (lam - m * (m + 1)))) / 2.0) \
                if lam + 0.26 > m * (m + 1) else 0
            cf = (m + j_near) * (m + j_near + 1)
            ents.append({"lambda": lam, "sector_n": n, "m": m,
                         "closed_form_(m+j)(m+j+1)": cf,
                         "closed_form_dev": lam - cf})
        n += 1
    ents.sort(key=lambda d: d["lambda"])
    mu_excl = (n * PI / (2.0 * W)) ** 2
    return ents, mu_excl, min_seen


# ---------------------------------------------------------------------------
# Realization (B): Frobenius connection at x = 1  (my derivation)
# ---------------------------------------------------------------------------

def seam_series(lam):
    """phi, phi_x, psi-part sums at x=1 (scalar lam; joint rescale for deep
    negative lam).  Returns (P, dP, Q, dQ, esc):
      phi(1) = P * e^esc, phi_x(1) = dP * e^esc,
      psi(1) = Q * e^esc, psi_x(1) = (P + dQ) * e^esc."""
    lam = float(lam)
    a, b = 1.0, 0.0
    P, dP, Q, dQ = 1.0, 0.0, 0.0, 0.0
    esc = 0.0
    al = abs(lam)
    kturn = 2.0 * math.sqrt(al / 2.0) + 60.0
    kmax = int(3.2 * math.sqrt(al / 2.0)) + 900
    for k in range(kmax):
        den = 2.0 * (k + 1) * (k + 1)
        fac = k * (k + 1) - lam
        a1 = a * fac / den
        b1 = ((1.0 + 2 * k) * a - 4.0 * (k + 1) * a1 + fac * b) / den
        a, b = a1, b1
        kk = k + 1
        P += a
        dP += kk * a
        Q += b
        dQ += kk * b
        if abs(a) > 1e260 or abs(b) > 1e260:
            s = 1e-260
            a *= s; b *= s; P *= s; dP *= s; Q *= s; dQ *= s
            esc += 598.67342
        if k > kturn and abs(a) + abs(b) < 1e-26 * (abs(P) + abs(dP)
                                                    + abs(Q) + abs(dQ) + 1e-280):
            break
    return P, dP, Q, dQ, esc


def gB(lam, delta0):
    """Branch-2 eigencondition psi_x(1) - (2 ln d0 - ln 2) phi_x(1), scaled.
    Sign-exact for any lam (common positive scale factor e^esc dropped)."""
    P, dP, Q, dQ, esc = seam_series(lam)
    return (P + dQ) - (2.0 * math.log(delta0) - LN2) * dP


def phi_at_1(lam):
    return seam_series(lam)[0]


def bisect_sign(f, lo, hi, iters=200, rtol=1e-14):
    flo, fhi = f(lo), f(hi)
    if flo == 0.0:
        return lo
    if fhi == 0.0:
        return hi
    assert flo * fhi < 0.0, f"no bracket: f({lo})={flo}, f({hi})={fhi}"
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        if fm == 0.0:
            return mid
        if flo * fm < 0.0:
            hi = mid
        else:
            lo, flo = mid, fm
        if hi - lo < rtol * (1.0 + abs(mid)):
            break
    return 0.5 * (lo + hi)


def branch2_positive_roots(delta0, lam_max=45.0, ngrid=4000):
    grid = np.linspace(1e-6, lam_max, ngrid)
    vals = np.array([gB(x, delta0) for x in grid])
    roots = []
    for i in range(len(grid) - 1):
        if vals[i] == 0.0:
            roots.append(float(grid[i]))
        elif vals[i] * vals[i + 1] < 0.0:
            roots.append(brentq(lambda x: gB(x, delta0), grid[i], grid[i + 1],
                                xtol=1e-13, rtol=8.9e-16))
    return roots


def branch1_positive_roots(lam_max=45.0, ngrid=4000):
    grid = np.linspace(1e-6, lam_max, ngrid)
    vals = np.array([phi_at_1(x) for x in grid])
    roots = []
    for i in range(len(grid) - 1):
        if vals[i] * vals[i + 1] < 0.0:
            roots.append(brentq(phi_at_1, grid[i], grid[i + 1],
                                xtol=1e-13, rtol=8.9e-16))
    return roots


def negative_roots(delta0, kap_lo_fac=0.02, kap_hi_fac=6.0, ngrid=320):
    """All branch-2 roots with lambda<0 in kappa in [fac_lo, fac_hi]*kap_est.
    Also asserts phi(1)>0 on the scan (no branch-1 negatives)."""
    kap_est = 2.0 * math.exp(-EULER_GAMMA) / delta0
    klo = max(1e-4, kap_lo_fac * kap_est)
    khi = kap_hi_fac * kap_est
    kg = np.geomspace(klo, khi, ngrid)
    vals = []
    for kap in kg:
        lam = -kap * kap
        P, dP, Q, dQ, esc = seam_series(lam)
        assert P > 0.0 and dP > 0.0, "branch-1 negative root?!"
        vals.append((P + dQ) - (2.0 * math.log(delta0) - LN2) * dP)
    roots = []
    for i in range(len(kg) - 1):
        if vals[i] * vals[i + 1] < 0.0:
            kap = bisect_sign(lambda K: gB(-K * K, delta0), kg[i], kg[i + 1])
            roots.append(-kap * kap)
    roots.sort()
    return roots


# ---------------------------------------------------------------------------
# Riccati shooting cross-check for the branch-2 negative eigenvalue
# ---------------------------------------------------------------------------

def local_phi_psi(lam, x, nk=400):
    """phi, phi', psi, psi' at small x from the cone-side series.
    Accumulates the TERMS A_k = a_k x^k, B_k = b_k x^k (the bare coefficients
    a_k overflow for |lam| ~ 1e6 even though the terms are moderate):
      A_{k+1} = A_k x fac/den,
      B_{k+1} = ((1+2k) A_k x - 4(k+1) A_{k+1} + fac B_k x)/den."""
    A, B = 1.0, 0.0
    phi, dphi, h, dh = 1.0, 0.0, 0.0, 0.0
    kturn = 2.0 * math.sqrt(abs(lam) * x / 2.0) + 20.0
    for k in range(nk):
        den = 2.0 * (k + 1) * (k + 1)
        fac = k * (k + 1) - lam
        A1 = A * x * fac / den
        B1 = ((1.0 + 2 * k) * A * x - 4.0 * (k + 1) * A1 + fac * B * x) / den
        A, B = A1, B1
        kk = k + 1
        phi += A
        h += B
        dphi += kk * A / x
        dh += kk * B / x
        if k > kturn and abs(A) + abs(B) < 1e-24 * (abs(phi) + abs(h) + 1e-280):
            break
    lnx = math.log(x)
    psi = phi * lnx + h
    dpsi = dphi * lnx + phi / x + dh
    return phi, dphi, psi, dpsi


def riccati_b_over_a(lam):
    """Riccati shooting for lam < 0 in the Bessel variable z = kappa sqrt(2x)
    (R = u_z/u stays O(1)): integrate from the seam (R=0 at z0 = kappa sqrt 2)
    down to z_m, extract the cone ratio b/a against the local (phi, psi) basis.
    In z the sector ODE is u_zz + u_z (2-3x)/((2-x) z) - 2 u/(2-x) = 0."""
    kappa = math.sqrt(-lam)
    x_m = min(1e-3, 50.0 / max(1.0, abs(lam)))
    z0 = kappa * math.sqrt(2.0)
    z_m = kappa * math.sqrt(2.0 * x_m)

    def rhs(z, R):
        x = z * z / (2.0 * kappa * kappa)
        return [-R[0] * R[0] - R[0] * (2.0 - 3.0 * x) / ((2.0 - x) * z)
                + 2.0 / (2.0 - x)]

    sol = solve_ivp(rhs, (z0, z_m), [0.0], method="LSODA",
                    rtol=1e-11, atol=1e-12)
    if not sol.success:
        return math.nan
    S = (kappa * kappa / z_m) * sol.y[0, -1]        # u_x/u at x_m
    phi, dphi, psi, dpsi = local_phi_psi(lam, x_m)
    return (dpsi - S * psi) / (S * phi - dphi)


def riccati_negative_root(delta0, kap_bracket):
    target = LN2 - 2.0 * math.log(delta0)

    def f(kap):
        return riccati_b_over_a(-kap * kap) - target

    kap = brentq(f, *kap_bracket, xtol=1e-12, rtol=1e-13)
    return -kap * kap


# ---------------------------------------------------------------------------
# Anchors (hand-derived identities, asserted numerically)
# ---------------------------------------------------------------------------

def check_anchors():
    P, dP, Q, dQ, esc = seam_series(2.0)
    assert abs(P) < 1e-13, P                    # phi(1; 2) = 0  (phi = t)
    assert abs(dP + 1.0) < 1e-13, dP            # phi_x(1; 2) = -1
    assert abs((P + dQ) - (2.0 - LN2)) < 1e-12  # psi_x(1; 2) = 2 - ln 2
    P0, dP0, Q0, dQ0, e0 = seam_series(0.0)
    assert abs(P0 - 1.0) < 1e-14 and abs(dP0) < 1e-14
    assert abs((P0 + dQ0) - 2.0) < 1e-12        # psi_x(1; 0) = 2  (exact soln)
    # collocation anchor: W=2.2 sector n=1 ground state = (m)(m+1), m=pi/4.4
    m = PI / 4.4
    e = sector_eigs(m, +1.0, 10.0)
    assert abs(e[0] - m * (m + 1.0)) < 1e-9, e[0]


# ---------------------------------------------------------------------------
# Drivers
# ---------------------------------------------------------------------------

def load_json(path):
    try:
        with open(path) as f:
            return json.load(f)
    except Exception:
        return None


def run_A1():
    other = load_json(os.path.join(DATA, "m8_1_spectrum.json"))
    W_overlap = [0.6, 1.5, 1.5708, 2.2, 3.0]
    W_new = [0.45, 2.0]
    out = {}
    print("== A1: realization (A), Chebyshev collocation (N=64/96 matched) ==")
    for W in W_overlap + W_new + [1.0]:
        ents, mu_excl, min_seen = spectrum_A(W)
        lowest6 = ents[:6]
        lam1p = next(e for e in ents if e["lambda"] > 1e-6)
        rec = {"W": W,
               "lowest6_audit": lowest6,
               "lambda1plus_audit": lam1p["lambda"],
               "lambda1plus_sector_n": lam1p["sector_n"],
               "min_eigenvalue_seen": min_seen,
               "mu_first_excluded_sector": mu_excl,
               "max_closed_form_dev_lowest6": max(abs(e["closed_form_dev"])
                                                  for e in lowest6)}
        key = f"{W:g}"
        if other and key in other.get("per_W", {}):
            th = other["per_W"][key]
            rec["lowest6_other_agent"] = [e["lambda"] for e in th["lowest6"]]
            rec["lambda1plus_other_agent"] = th["lambda1plus"]["lambda"]
            rec["lowest6_absdiff"] = [abs(a["lambda"] - b["lambda"]) for a, b in
                                      zip(lowest6, th["lowest6"])]
        out[key] = rec
        print(f"W={W:8g} lam1+={lam1p['lambda']:.10f} (n={lam1p['sector_n']})  "
              "lowest6=" + " ".join(f"{e['lambda']:.8f}(n={e['sector_n']})"
                                    for e in lowest6))
    return out


def run_A2():
    other = load_json(os.path.join(DATA, "m8_1_delta_scan.json"))
    other_by_d0 = {}
    if other:
        for row in other["per_delta0"]:
            other_by_d0[row["delta0"]] = row
    print("== A2: realization (B) at W=1, independent trace matching ==")
    b1 = branch1_positive_roots()
    print("branch-1 (delta0-independent) roots:", [f"{r:.10f}" for r in b1[:4]])
    # mu>0 sector minimum at W=1 from the audit collocation
    m1 = PI / 2.0
    sec1 = sector_eigs(m1, +1.0, 20.0)
    nge1_min = sec1[0]
    print(f"W=1 n=1 sector minimum (collocation): {nge1_min:.10f} "
          f"[closed form {m1 * (m1 + 1):.10f}]")

    d0_spec = [0.001, 0.1, 0.71969, 1.3895]
    d0_grid_variants = {0.71969: 0.7196856730011514, 1.3895: 1.3894954943731375}
    per_d0 = {}
    for d0 in d0_spec + list(d0_grid_variants.values()):
        negs = negative_roots(d0)
        pos2 = branch2_positive_roots(d0)
        cands = sorted([r for r in (b1 + pos2) if r > 1e-6] + [nge1_min])
        first_pos = cands[0]
        rec = {"delta0": d0,
               "negative_eigenvalues_audit": negs,
               "n_negative_audit": len(negs),
               "first_positive_audit": first_pos,
               "branch2_positive_first3": pos2[:3]}
        d0_match = min(other_by_d0, key=lambda x: abs(x - d0)) if other_by_d0 else None
        if d0_match is not None and abs(d0_match - d0) < 5e-3 * d0:
            row = other_by_d0[d0_match]
            rec["other_agent_delta0"] = d0_match
            rec["other_agent_negative"] = row["negative_eigenvalues"]
            rec["other_agent_first_positive"] = row["first_positive"]
        per_d0[f"{d0:.10g}"] = rec
        print(f"delta0={d0:<12.8g} neg={negs}  first_pos={first_pos:.10f}")

    # (a) transition delta0*: smallest branch-2 positive root crosses 2
    def rho(d0):
        return branch2_positive_roots(d0, lam_max=5.9, ngrid=1200)[0]

    d0_star = bisect_sign(lambda d: rho(d) - 2.0, 0.4, 1.2, rtol=1e-12)
    d0_star_closed = 2.0 / math.e
    print(f"delta0* (bisection)   = {d0_star:.10f}")
    print(f"delta0* closed form 2/e = {d0_star_closed:.10f}")

    # (b) small-delta0 law of the negative eigenvalue
    Cvals = {}
    for d0 in (1e-3, 3e-4, 1e-4):
        kap_est = 2.0 * math.exp(-EULER_GAMMA) / d0
        kap = bisect_sign(lambda K: gB(-K * K, d0),
                          0.8 * kap_est, 1.25 * kap_est)
        lam_neg = -kap * kap
        Cvals[f"{d0:g}"] = {"lambda_neg": lam_neg, "C": lam_neg * d0 * d0}
        print(f"C({d0:g}) = {lam_neg * d0 * d0:.10f}   (lambda_neg={lam_neg:.6e})")
    C_34, C_14 = Cvals["0.0003"]["C"], Cvals["0.0001"]["C"]
    C0_extrap = (9.0 * C_14 - C_34) / 8.0          # Richardson in delta0^2
    C0_analytic = -4.0 * math.exp(-2.0 * EULER_GAMMA)
    print(f"C0 extrapolated = {C0_extrap:.10f}")
    print(f"C0 analytic  -4 e^-2gamma = {C0_analytic:.10f}")

    # Riccati shooting cross-checks (independent discretization family)
    ric = {}
    for d0 in (0.1, 0.001):
        kap_est = 2.0 * math.exp(-EULER_GAMMA) / d0
        lam_r = riccati_negative_root(d0, (0.85 * kap_est, 1.2 * kap_est))
        ric[f"{d0:g}"] = lam_r
        print(f"Riccati cross-check lambda_neg(delta0={d0:g}) = {lam_r:.8f}")

    return {"W": 1.0,
            "branch1_roots_delta0_independent": b1[:4],
            "n_ge1_sector_min_W1_audit": nge1_min,
            "per_delta0": per_d0,
            "delta0_star_bisection": d0_star,
            "delta0_star_closed_form_2_over_e": d0_star_closed,
            "C_of_delta0": Cvals,
            "C0_extrapolated": C0_extrap,
            "C0_analytic_minus4_exp_minus2gamma": C0_analytic,
            "riccati_crosscheck_lambda_neg": ric}


def run_A3_richardson_replay():
    """(vi): replay the audited Richardson rule on every raw ladder in the
    audited spectrum JSON and compare with its reported lambda/err."""
    other = load_json(os.path.join(DATA, "m8_1_spectrum.json"))
    if not other:
        return None

    def extrap3(a, b, c):
        d1, d2 = b - a, c - b
        if abs(d2) < 1e-13 or abs(d1) < 1e-13:
            return c
        r = d1 / d2
        if r <= 1.0:
            return c + d2
        return c + d2 / (r - 1.0)

    worst = 0.0
    nchecked = 0
    for key, blk in other["per_W"].items():
        for e in blk["lowest6"] + [blk["lambda1plus"]]:
            v = e["raw_by_resolution"]
            lam_b = extrap3(v[-3], v[-2], v[-1])
            lam_a = extrap3(v[-4], v[-3], v[-2])
            err = max(abs(lam_b - lam_a), 1e-12)
            worst = max(worst,
                        abs(lam_b - e["lambda"]) / (1.0 + abs(e["lambda"])),
                        abs(err - e["err"]) / (1.0 + e["err"]))
            nchecked += 1
    return {"ladders_replayed": nchecked,
            "max_rel_discrepancy_lambda_and_err": worst}


def main():
    os.makedirs(DATA, exist_ok=True)
    check_anchors()
    print("anchors OK: phi(1;2)=0, phi_x(1;2)=-1, psi_x(1;2)=2-ln2, "
          "psi_x(1;0)=2, colloc n=1 W=2.2 = m(m+1)")
    a1 = run_A1()
    a2 = run_A2()
    a3vi = run_A3_richardson_replay()
    out = {
        "audit_of": ["m8_1_eigensolve.py", "m8_1_eigensolve_xcheck.py",
                     "m8_1_spectrum.json", "m8_1_delta_scan.json",
                     "m8_1_xcheck.json"],
        "audit_methods": {
            "A_realization": "two-half Chebyshev collocation in t=sin y, "
                             "singular weight factored (1-t^2)^{m/2}, "
                             "N=64/96 cross-matched (no FEM, no cone series)",
            "closed_form_found": "sector eigenvalues lambda=(m+j)(m+j+1), "
                                 "m=n*pi/(2W), j>=0 (Gegenbauer)",
            "B_realization": "independent Frobenius connection at x=1 "
                             "(own normalization psi=phi ln x + h) + Riccati "
                             "log-derivative shooting cross-check",
        },
        "A1": a1,
        "A2": a2,
        "A3_vi_richardson_replay": a3vi,
    }
    path = os.path.join(DATA, "m8_1_audit.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=1, default=float)
    print("written", path)


if __name__ == "__main__":
    main()
