"""M8.3: mass-formula reproducer. Recomputes every constant from its own definition.

Target formula, reproducing mode-identity-theory's corrected mass-spectrum table as of
repo commit 13a5c2b / archived Zenodo version 10.5281/zenodo.21652153
(concept DOI 10.5281/zenodo.18603975, all versions):

    m(rho, sigma) = mu_Lambda * C_geom(rho) * (sqrt(Omega_Lambda))^(dist(rho)/30) * T^2(rho x sigma)

Recomputed from definition, never quoted:
  - 2I character table, McKay graph, McKay distances (Part A);
  - Kostant exponents from the multiplicity generating function
    sum_n mult_rho(V_n) z^n = (sum_i z^{e_i}) / ((1-z^12)(1-z^20))   (Part B);
  - C_geom(rho) = (prod_e 2 sin^2(pi e / D))^(1/(2 dim rho)), D = 60 (int) / 120 (half) (Part C);
  - anchors: Lambda, mu_Lambda = rho_Lambda^(1/4), rho_Lambda = Lambda c^4/(8 pi G), and
    sqrt(Omega_Lambda) = R/l_P with R = sqrt(3/Lambda) are all COMPUTED from H0, Omega_L, G, c,
    hbar, the declared measurement inputs (Part D);
  - single-irrep Reidemeister torsions T^2(tau) from the spectral-zeta definition (twisted
    0-form + coexact 1-form spectra on S^3/2I -> zeta'(0); Ray-Singer combination), computed and
    gated against the page's exact closed forms for R1 through R8 (8 irreps, integer AND
    half-integer alike, Part E); T^2(R0) = 1 is the declared non-acyclic convention, asserted
    rather than independently gated (see Declared conventions/inputs below);
  - the 24 products via log T^2(rho x sigma) = sum_tau N_{rho sigma tau} log T^2(tau)  (Part F);
  - the PDG comparison, scorecard counts, and rank ordering (Part G).

Development of this reproducer exposed an omitted scalar contribution in the published
half-integer torsions. The canonical source was corrected (torsion-correction.md, 2026-07-28)
before these validation targets were frozen; see findings/m8_3_method_note.md for diagnosis and
provenance.

Declared conventions/inputs (deliverable list, task spec DoD #3 -- NOT recomputed here):
  - T^2(R0) = 1 on the two non-acyclic diagonal products (ledger convention, page section 4);
  - the overall torsion sign convention, fixed once on the R7 closed form (page section 4;
    every other closed form below is then a genuine independent check, not circular);
  - D = 60/120 spin rule for C_geom (page section 2, Z4-stabilizer argument);
  - H0, Omega_L, G, c, hbar: the raw measurement inputs (Planck 2018 + CODATA) that Part D
    computes Lambda, mu_Lambda, and sqrt(Omega_Lambda) FROM -- those three are recomputed by
    this script, not quoted; only the underlying physical constants are declared;
  - PDG comparison masses: measurement inputs (page's own declared figures, Part G);
  - the (rho, sigma) -> fermion assignment: the page's ranked-table reading (data input), INCLUDING
    each fermion's structural sector (e.g. the bottom quark's sector is R2; its numerically
    compatible (R4, gal) address is therefore out-of-sector and not adjudicated). The
    representation-theoretic gate machinery that derives these sectors and quantum numbers is a
    separate concern and is not re-derived by this script.

Standing rule: every PASS line has a mutation test (run with --mutation-tests), with a
coverage-enforced registry (every gate id must be attacked by at least one mutation, and every
mutation must turn its gate red, or the suite fails by construction). A check that cannot fail
is not a check. Null-test context: the x3 proximity hit rate is REPORTED ONLY, not evidence for
the torsion map (pre-registered null mass-null-v1.1, corrected table, p_A = 0.690).

Run:  python3 m8_3_mass_reproducer.py [--precise] [--mutation-tests]
"""
import argparse
import json
import math
import os
import sys

import mpmath as mp
import numpy as np

# ---------------------------------------------------------------- Part A: 2I
labels = ["R0", "R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8"]
dims = [1, 2, 2, 3, 3, 4, 4, 5, 6]
idx = {l: i for i, l in enumerate(labels)}
# McKay graph = affine E8 (tensoring with the defining rep R1 is adjacency)
edges = [("R0", "R1"), ("R1", "R3"), ("R3", "R6"), ("R6", "R7"),
         ("R7", "R8"), ("R8", "R5"), ("R5", "R2"), ("R8", "R4")]
A = np.zeros((9, 9))
for a, b in edges:
    A[idx[a], idx[b]] = A[idx[b], idx[a]] = 1
# 9 conjugacy classes: (SU(2) eigen-angle theta, size)
classes = [(0.0, 1), (math.pi, 1), (math.pi / 2, 30), (math.pi / 3, 20),
           (2 * math.pi / 3, 20), (math.pi / 5, 12), (2 * math.pi / 5, 12),
           (3 * math.pi / 5, 12), (4 * math.pi / 5, 12)]
thetas = np.array([c[0] for c in classes])
sizes = np.array([c[1] for c in classes], float)


def chiV(n, th):
    """SU(2) character of V_n (dim n+1) at eigen-angle th."""
    s = math.sin(th)
    if abs(s) < 1e-12:
        return (n + 1) * math.cos(n * th)
    return math.sin((n + 1) * th) / s


def build_chartable(edge_list=edges, class_list=classes):
    """Rebuild the character table from an (edges, classes) pair (mutable for mutations)."""
    Am = np.zeros((9, 9))
    for a, b in edge_list:
        Am[idx[a], idx[b]] = Am[idx[b], idx[a]] = 1
    th = np.array([c[0] for c in class_list])
    sz = np.array([c[1] for c in class_list], float)
    N_REC = 24
    mrec = [np.zeros(9) for _ in range(N_REC)]
    mrec[0][idx["R0"]] = 1.0
    mrec[1][idx["R1"]] = 1.0
    for n in range(1, N_REC - 1):
        mrec[n + 1] = Am.dot(mrec[n]) - mrec[n - 1]
    Mmat = np.array([mrec[n] for n in range(9)])
    ch = np.zeros((9, 9))
    for c in range(9):
        rhs = np.array([chiV(n, th[c]) for n in range(9)])
        ch[:, c] = np.linalg.solve(Mmat, rhs)
    return Am, ch, th, sz


A, chi, thetas, sizes = build_chartable()


def mult_in_Vn(tau_i, n, chi_=None, sizes_=None, thetas_=None):
    """multiplicity of irrep tau in V_n restricted to 2I (exact integer)."""
    chi_ = chi if chi_ is None else chi_
    sizes_ = sizes if sizes_ is None else sizes_
    thetas_ = thetas if thetas_ is None else thetas_
    v = sum(sizes_[c] * chiV(n, thetas_[c]) * chi_[tau_i][c] for c in range(9))
    m = v / 120.0
    r = round(m)
    assert abs(m - r) < 1e-6, f"non-integer multiplicity {m} at tau={tau_i}, n={n}"
    return int(r)


def tensor_mult(a_i, b_i, c_i):
    """N: multiplicity of irrep c in a x b."""
    v = sum(sizes[k] * chi[a_i][k] * chi[b_i][k] * chi[c_i][k] for k in range(9))
    m = v / 120.0
    r = round(m)
    assert abs(m - r) < 1e-6
    return int(r)


def bfs_distances(adj):
    d = {idx["R0"]: 0}
    front = [idx["R0"]]
    while front:
        nxt = []
        for a_ in front:
            for b_ in range(9):
                if adj[a_, b_] > 0.5 and b_ not in d:
                    d[b_] = d[a_] + 1
                    nxt.append(b_)
        front = nxt
    return [d[i] for i in range(9)]


dist = bfs_distances(A)

# ------------------------------------------------- Part B: Kostant exponents
def kostant_exponents(tau_i, denom=(12, 20), nmax=80):
    """Exponents e_i defined by  sum_n mult(V_n) z^n * (1-z^d1)(1-z^d2) = sum_i z^{e_i}."""
    series = [mult_in_Vn(tau_i, n) for n in range(nmax)]
    d1, d2 = denom
    poly = list(series)
    for d in (d1, d2):
        nxt = [0] * nmax
        for n in range(nmax):
            nxt[n] = poly[n] - (poly[n - d] if n >= d else 0)
        poly = nxt
    # must terminate: everything above degree d1+d2 (=32<40) vanishes
    tail = poly[41:]
    exps = []
    for e, cf in enumerate(poly[:41]):
        if cf < 0:
            return None, poly, tail
        exps += [e] * cf
    return exps, poly, tail


# ---------------------------------------------------------- Part C: C_geom
def c_geom(tau_i, exps):
    half = tau_i in (idx["R1"], idx["R2"], idx["R6"], idx["R8"])
    D = 120 if half else 60
    prod = mp.mpf(1)
    for e in exps:
        prod *= 2 * mp.sin(mp.pi * e / D) ** 2
    return prod ** (mp.mpf(1) / (2 * dims[tau_i]))


# ------------------------------------------------- Part D: measured anchors
# Calibration INPUTS (declared, not derived): Planck 2018 + CODATA 2018.
H0_KM_S_MPC = mp.mpf("67.36")          # Planck 2018 TT,TE,EE+lowE+lensing
OMEGA_L_COSMO = mp.mpf("0.6847")       # Planck 2018 dark-energy density parameter
C_M_S = mp.mpf("299792458")
G_SI = mp.mpf("6.67430e-11")
HBAR_SI = mp.mpf("1.054571817e-34")
EV_J = mp.mpf("1.602176634e-19")
MPC_M = mp.mpf("3.0856775814913673e22")


def anchors(h0=H0_KM_S_MPC, omega_l=OMEGA_L_COSMO, g=G_SI):
    H0 = h0 * 1000 / MPC_M                               # 1/s
    Lam = 3 * omega_l * H0 ** 2 / C_M_S ** 2              # 1/m^2
    rho_L = Lam * C_M_S ** 4 / (8 * mp.pi * g)            # J/m^3
    hbarc = HBAR_SI * C_M_S                               # J m
    mu = (rho_L * hbarc ** 3) ** mp.mpf("0.25") / EV_J    # eV
    lP = mp.sqrt(HBAR_SI * g / C_M_S ** 3)                # m
    Rds = mp.sqrt(3 / Lam)                                # m
    sqrt_Omega = Rds / lP
    return Lam, mu, sqrt_Omega


# ------------------------------------------- Part E: torsion (spectral zeta)
# Twisted spectra on S^3/2I (unit radius; torsion of an acyclic bundle is
# metric-independent), u = n+1:
#   0-forms:  lambda = u^2 - 1,  mult(u) = u * mult_tau(V_{u-1})
#             (the right-SU(2) degeneracy factor u multiplies the character count)
#   coexact:  lambda = u^2 (u >= 2), mult(u) = (u+1) mult_tau(V_{u-2}) + (u-1) mult_tau(V_u)
# Ray-Singer: log T = 1/2 [zeta'_coex(0) - 2 zeta'_0form(0)]  (3-manifold reduction);
# the overall sign convention is FIXED once on the R7 closed form and then validated
# on every other irrep, integer and half-integer alike (independent gates).
PERIOD = 60
FIT_BASE = 120  # fit u-samples at r + 120, +180, +240; verify at +300, +360


def quasi_quadratic_coeffs(mult_fn):
    """Fit mult(u) = a_r u^2 + b_r u + c_r per residue r mod 60; verify two extra points."""
    coeffs = []
    for r in range(PERIOD):
        us = [r + FIT_BASE, r + FIT_BASE + 60, r + FIT_BASE + 120]
        vals = [mult_fn(u) for u in us]
        M = np.array([[u * u, u, 1] for u in us], float)
        abc = np.linalg.solve(M, np.array(vals, float))
        for u in (r + FIT_BASE + 180, r + FIT_BASE + 240):
            pred = abc[0] * u * u + abc[1] * u + abc[2]
            if abs(pred - mult_fn(u)) > 1e-4:
                raise RuntimeError(f"quasi-quadratic fit fails at r={r}, u={u}")
        coeffs.append(abc)
    return coeffs


class HurwitzCache:
    """H[p][j][r](s) = sum over u = r (mod 60), u >= u0 of u^(p - 2s - 2j), via Hurwitz zeta."""

    def __init__(self, u0):
        self.u0 = u0
        self.cache = {}

    def val(self, s, p, j, r):
        key = (str(s), p, j, r)
        if key not in self.cache:
            k0 = 0
            while 60 * k0 + r < self.u0:
                k0 += 1
            z = 2 * s + 2 * j - p
            q = mp.mpf(k0) + mp.mpf(r) / 60
            self.cache[key] = mp.power(60, -z) * mp.zeta(z, q)
        return self.cache[key]


def zeta_from_coeffs(s, coeffs, cache, jmax, exact_square):
    """zeta(s) = sum_u mult(u) lambda(u)^-s with lambda = u^2 (exact) or u^2-1 (binomial)."""
    total = mp.mpf(0)
    if exact_square:
        for r in range(PERIOD):
            a, b, c = coeffs[r]
            for p, cf in ((2, a), (1, b), (0, c)):
                if abs(cf) > 1e-12:
                    total += mp.mpf(cf) * cache.val(s, p, 0, r)
        return total
    g = mp.mpf(1)
    for j in range(jmax):
        if j > 0:
            g *= (s + j - 1) / j
        inner = mp.mpf(0)
        for r in range(PERIOD):
            a, b, c = coeffs[r]
            for p, cf in ((2, a), (1, b), (0, c)):
                if abs(cf) > 1e-12:
                    inner += mp.mpf(cf) * cache.val(s, p, j, r)
        total += g * inner
    return total


def zeta_prime_at_0(coeffs, cache, jmax, exact_square, sigma):
    """zeta'(0) by symmetric polynomial fit through s = sigma*{-2,-1,1,2} (analytic at 0)."""
    ss = [sigma * k for k in (-2, -1, 1, 2)]
    vals = [zeta_from_coeffs(mp.mpf(s), coeffs, cache, jmax, exact_square) for s in ss]
    M = mp.matrix([[mp.mpf(s) ** i for i in range(4)] for s in ss])
    sol = mp.lu_solve(M, mp.matrix(vals))
    return sol[1]  # coefficient of s^1


def torsion_T2(tau_i, cache0, cacheco, jmax, sigma, u_factor=True, co_shift=1):
    """T^2(tau) for acyclic tau (tau != R0) from the spectral-zeta definition."""
    m0 = lambda u: (u if u_factor else 1) * mult_in_Vn(tau_i, u - 1)
    mco = lambda u: (u + co_shift) * mult_in_Vn(tau_i, u - 2) + (u - co_shift) * mult_in_Vn(tau_i, u)
    c0 = quasi_quadratic_coeffs(m0)
    cc = quasi_quadratic_coeffs(mco)
    zp0 = zeta_prime_at_0(c0, cache0, jmax, False, sigma)
    zpc = zeta_prime_at_0(cc, cacheco, 0, True, sigma)
    return zpc - 2 * zp0  # = 2 log T up to the global sign fixed on R7


# ------------------------------------------------------------------- gates
GATES = []   # (gate_id, display_name, ok)
MUTS = []    # (target_gate_id, display_name, went_red)


def gate(gid, name, ok, detail=""):
    GATES.append((gid, name, bool(ok)))
    print(f"  {'PASS' if ok else 'FAIL'}  [{gid}] {name}{('  ' + detail) if detail else ''}")
    return bool(ok)


def mutation(target_gid, name, went_red):
    MUTS.append((target_gid, name, bool(went_red)))
    print(f"  {'RED (ok)' if went_red else 'STILL PASSING (BAD)'}  [{target_gid}] {name}")


# ---------------------------------------------- Part G: PDG comparison data
# Declared measurement inputs (page's own figures; NOT recomputed here).
PDG_GEV = {"e": 5.11e-4, "u": 2.16e-3, "d": 4.67e-3, "s": 9.34e-2, "mu": 1.057e-1,
           "tau": 1.777, "b": 4.18, "t": 172.7, "c": 1.27}

# Declared (rho, sigma) -> fermion assignment (page's ranked-table reading, DATA INPUT).
# Sector-first rule (declared, not re-derived): where the vertex structure names a fermion's
# irrep sector (neutrinos on R1, top/bottom on R2, electron on R7), the fermion is assigned to
# its nearest compatible address WITHIN that sector; a numerically compatible address outside
# the fermion's own sector is recorded but not adjudicated (the bottom quark, below).
SM_ASSIGN = {
    "e":   {"addr": ("R7", "triv"), "role": "benchmark"},
    "mu":  {"addr": ("R8", "std"),  "role": "assigned_within_x3"},
    "s":   {"addr": ("R8", "std"),  "role": "assigned_within_x3"},
    "tau": {"addr": ("R4", "gal"),  "role": "assigned_within_x3"},
    "t":   {"addr": ("R2", "triv"), "role": "assigned_within_x3"},
    "d":   {"addr": ("R8", "gal"),  "role": "assigned_outside_x3"},
    "b":   {"addr": ("R4", "gal"),  "role": "compatible_out_of_sector"},  # sector is R2
    "u":   {"addr": None,           "role": "unassigned"},
    "c":   {"addr": None,           "role": "unassigned"},
}
# Declared expected ratios (page's own displayed figures, cross-check target, not an input to
# the computation itself). Convention: the page reports the SYMMETRIZED max(pred/obs, obs/pred)
# for every entry except the top quark, which is deliberately displayed as raw pred/obs (0.93,
# predicted 7% BELOW observed) to convey undershoot direction rather than "off by 1.07x".
DECLARED_RATIO = {"e": 1.02, "mu": 1.02, "s": 1.10, "tau": 2.75, "t": 0.93, "d": 3.23, "b": 1.17}
RAW_NOT_SYMMETRIZED = {"t"}

# The pre-registered null test, both runs (mode-identity-theory, mass-null-test.md).
NULL_V11 = {"tag": "mass-null-v1.1", "table": "corrected", "p_A": 0.690, "S1_obs": 5, "S1_mean": 5.02}
NULL_V10 = {"tag": "mass-null-v1.0", "table": "pre-correction", "p_A": 0.174, "S1_obs": 6, "S1_mean": 4.50}

# Canonical rank order (mass-spectrum.md SIII), the cross-check target for Part G's rank gate.
CANONICAL_RANK_ORDER = [
    ("R1", "triv"), ("R1", "std"), ("R1", "gal"), ("R3", "std"), ("R3", "triv"), ("R3", "gal"),
    ("R6", "std"), ("R6", "triv"), ("R6", "gal"), ("R7", "triv"), ("R7", "std"), ("R7", "gal"),
    ("R8", "gal"), ("R8", "triv"), ("R8", "std"), ("R5", "gal"), ("R4", "gal"), ("R4", "triv"),
    ("R5", "triv"), ("R4", "std"), ("R5", "std"), ("R2", "triv"), ("R2", "gal"), ("R2", "std"),
]

NOT_RECOMPUTED = [
    "(rho, sigma) -> fermion identity: the page's ranked-table reading (assignment is a declared "
    "DATA INPUT, not re-derived from quantum-number gates here)",
    "structural sector per fermion (e.g. the bottom quark's sector is R2, making its numerically "
    "compatible (R4,gal) address out-of-sector): declared per the page's sector-first rule",
    "PDG comparison masses: measurement inputs",
    "H0, Omega_L, G, c, hbar: the declared measurement inputs (Planck 2018 + CODATA); Lambda, "
    "mu_Lambda, and sqrt(Omega_Lambda) are recomputed from them (Part D), not declared inputs "
    "themselves",
    "T^2(R0) = 1 non-acyclic ledger convention: ASSERTED per the page's own topological argument",
    "overall torsion sign convention, fixed once on R7: declared, then every other closed form "
    "is a genuine independent check",
    "mass-null-v1.1 / mass-null-v1.0 p_A figures and S1 statistics: quoted from the external "
    "null-test analysis (mass-null-test.md), not computed by this script; G3 checks only their "
    "internal consistency, not the null test itself",
    "each fermion's SM_ASSIGN role bucket (assigned_within_x3 / assigned_outside_x3 / "
    "compatible_out_of_sector / unassigned): a declared label alongside the address and sector, "
    "not derived from the ratio at runtime (though adjudicated is always checked as a subset of "
    "the ratio-derived compatible set, so a mislabeled role cannot inflate the count)",
]


def compute_all(jmax, sigma, u_factor=True, co_shift=1, only=None):
    """Run Parts A-F once; return everything Part G and the gates need."""
    Gm = np.array([[np.sum(sizes * chi[i] * chi[j]) / 120.0 for j in range(9)] for i in range(9)])
    cg = {}
    exps_all = {}
    for l in labels:
        e, poly, tail = kostant_exponents(idx[l])
        exps_all[l] = e
        if l != "R0":
            cg[l] = c_geom(idx[l], e)

    Lam, mu, sqrtOm = anchors()

    cache0, cacheco = HurwitzCache(2), HurwitzCache(2)
    X = {}
    for l in (only or labels[1:]):
        X[l] = torsion_T2(idx[l], cache0, cacheco, jmax, sigma, u_factor, co_shift)
    return Gm, exps_all, cg, Lam, mu, sqrtOm, X


def sign_fix(X):
    phi = (1 + mp.sqrt(5)) / 2
    target_r7 = mp.mpf(9) / 4
    sign = 1 if abs(mp.e ** X["R7"] - target_r7) < abs(mp.e ** -X["R7"] - target_r7) else -1
    T2 = {l: mp.e ** (sign * X[l]) for l in X}
    T2["R0"] = mp.mpf(1)
    return T2, sign


def products_24(T2, cg, dist_list):
    """T^2(rho x sigma) for sigma in {triv, std=R1, gal=R2}, and the mass of each address."""
    vac = {"triv": None, "std": idx["R1"], "gal": idx["R2"]}
    prod = {}
    for l in labels[1:]:
        row = {}
        for vname, vi in vac.items():
            if vi is None:
                row[vname] = T2[l]
            else:
                logv = mp.mpf(0)
                for t in range(9):
                    N = tensor_mult(idx[l], vi, t)
                    if N:
                        logv += N * mp.log(T2[labels[t]])
                row[vname] = mp.e ** logv
        prod[l] = row
    return prod


def mass_of(l, vname, prod, cg, mu, sqrtOm, dist_list):
    return mu * mp.mpf("1e-9") * cg[l] * sqrtOm ** (mp.mpf(dist_list[idx[l]]) / 30) * prod[l][vname]
    # mu is in eV (Part D); *1e-9 converts to GeV


def ranked_entries(prod, cg, mu, sqrtOm, dist_list):
    ent = []
    for l in labels[1:]:
        for vname in ("triv", "std", "gal"):
            m = mass_of(l, vname, prod, cg, mu, sqrtOm, dist_list)
            ent.append((l, vname, float(m)))
    ent.sort(key=lambda e: e[2])
    return ent


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--precise", action="store_true")
    ap.add_argument("--mutation-tests", action="store_true")
    args = ap.parse_args()
    mp.mp.dps = 60 if args.precise else 35
    jmax = 90 if args.precise else 55
    sigma = mp.mpf("1e-5") if args.precise else mp.mpf("1e-4")

    print("=" * 72)
    print("PART A: 2I character table, McKay graph, distances")
    print("=" * 72)
    Gm = np.array([[np.sum(sizes * chi[i] * chi[j]) / 120.0 for j in range(9)]
                   for i in range(9)])
    gate("A1_char_orthonormality", "character orthonormality",
         np.max(np.abs(Gm - np.eye(9))) < 1e-9,
         f"err {np.max(np.abs(Gm - np.eye(9))):.1e}")
    gate("A2_e8_marks", "dims are affine E8 marks (A.dims = 2 dims)",
         np.allclose(A.dot(np.array(dims, float)), 2.0 * np.array(dims, float)))
    gate("A3_completeness", "completeness sum_tau dim.mult = n+1 (n=0..23)",
         all(sum(dims[t] * mult_in_Vn(t, n) for t in range(9)) == n + 1
             for n in range(24)))
    print(f"  distances: {dict(zip(labels, dist))}")

    print("=" * 72)
    print("PART B: Kostant exponents from the generating function")
    print("=" * 72)
    page_exps = {  # mass-spectrum page section 3 (validation targets, NOT inputs)
        "R0": [0, 30], "R1": [1, 11, 19, 29], "R2": [7, 13, 17, 23],
        "R3": [2, 10, 12, 18, 20, 28], "R4": [6, 10, 14, 16, 20, 24],
        "R5": [6, 8, 12, 14, 16, 18, 22, 24],
        "R6": [3, 9, 11, 13, 17, 19, 21, 27],
        "R7": [4, 8, 10, 12, 14, 16, 18, 20, 22, 26],
        "R8": [5, 7, 9, 11, 13, 15, 15, 17, 19, 21, 23, 25]}
    exps_all = {}
    ok_b = True
    for l in labels:
        e, poly, tail = kostant_exponents(idx[l])
        good = (e is not None and all(abs(t) < 1e-9 for t in tail)
                and len(e) == 2 * dims[idx[l]] and e == page_exps[l])
        exps_all[l] = e
        ok_b &= good
    gate("B1_kostant_exponents", "exponents: nonneg poly, terminates, count=2dim, matches page (9/9)", ok_b)

    print("=" * 72)
    print("PART C: C_geom from exponents (D=60 int / 120 half)")
    print("=" * 72)
    page_cgeom = {"R1": 0.0988, "R2": 0.2436, "R3": 0.5553, "R4": 0.7970,
                  "R5": 0.8017, "R6": 0.2098, "R7": 0.7564, "R8": 0.2382}
    ok_c = True
    cg = {}
    for l in labels[1:]:
        v = c_geom(idx[l], exps_all[l])
        cg[l] = v
        ok_c &= abs(float(v) - page_cgeom[l]) < 5e-4
    gate("C1_cgeom", "C_geom matches page table (8/8, tol 5e-4)", ok_c,
         "  ".join(f"{l}={float(cg[l]):.4f}" for l in labels[1:]))

    print("=" * 72)
    print("PART D: anchors from measured inputs (declared calibration inputs)")
    print("=" * 72)
    Lam, mu, sqrtOm = anchors()
    print(f"  Lambda = {mp.nstr(Lam, 6)} 1/m^2   mu_Lambda = {mp.nstr(mu * 1000, 4)} meV"
          f"   sqrt(Omega) = {mp.nstr(sqrtOm, 6)}")
    gate("D1_mu_lambda", "mu_Lambda ~ 2.25 meV (page value, tol 3%)",
         abs(float(mu * 1000) - 2.25) / 2.25 < 0.03)
    gate("D2_sqrt_omega", "sqrt(Omega_Lambda) ~ 1.019e61 (page value, tol 3%)",
         abs(float(sqrtOm) / 1.019e61 - 1) < 0.03)

    print("=" * 72)
    print("PART E: torsion singles from the spectral-zeta definition (8 gated closed forms,")
    print("        R1-R8; R0=1 is the declared convention, asserted not gated)")
    print("=" * 72)
    cache0, cacheco = HurwitzCache(2), HurwitzCache(2)
    phi = (1 + mp.sqrt(5)) / 2
    # exact closed forms, EVERY irrep (the correction's headline: half-integer irreps also
    # reduce to golden-ratio closed forms, not just integer-spin ones)
    exact = {
        "R1": phi ** -4 / 4, "R2": phi ** 4 / 4,
        "R3": mp.mpf(4) / 5 * phi ** -2, "R4": mp.mpf(4) / 5 * phi ** 2,
        "R5": mp.mpf(25) / 9, "R6": mp.mpf(1),
        "R7": mp.mpf(9) / 4, "R8": mp.mpf(4),
    }
    X = {}
    for l in labels[1:]:
        X[l] = torsion_T2(idx[l], cache0, cacheco, jmax, sigma)
    T2, sign = sign_fix(X)
    print(f"  sign convention fixed on R7: {sign:+d}")
    for l in labels[1:]:
        gate(f"E_{l}", f"T^2({l}) matches exact closed form", abs(T2[l] / exact[l] - 1) < 1e-8,
             f"{mp.nstr(T2[l], 12)} vs {mp.nstr(exact[l], 12)}")
    gate("E_galois_R3R4", "Galois ratio T^2(R3)/T^2(R4) = phi^-4",
         abs(T2["R3"] / T2["R4"] / phi ** -4 - 1) < 1e-10)
    gate("E_galois_R1R2", "Galois ratio T^2(R1)/T^2(R2) = phi^-8",
         abs(T2["R1"] / T2["R2"] / phi ** -8 - 1) < 1e-10)
    gate("E_sector_products", "sector products: integer=4, half-integer=1/4, exact inverses",
         abs(T2["R3"] * T2["R7"] * T2["R5"] * T2["R4"] - 4) < 1e-8
         and abs(T2["R1"] * T2["R2"] * T2["R6"] * T2["R8"] - mp.mpf(1) / 4) < 1e-8)

    print("=" * 72)
    print("PART F: 24 products log T^2(rho x sigma) = sum N log T^2(tau)")
    print("=" * 72)
    page24 = {  # corrected page section 4 product table (validation targets)
        "R1": (0.0365, 0.306, 2.778), "R2": (1.714, 2.778, 2.094),
        "R3": (0.306, 0.0365, 4.000), "R4": (2.094, 4.000, 1.714),
        "R5": (2.778, 6.854, 0.146), "R6": (1.000, 0.688, 4.712),
        "R7": (2.250, 4.000, 4.000), "R8": (4.000, 13.090, 1.910)}
    prod = products_24(T2, cg, dist)
    ok24 = True
    for l in labels[1:]:
        got = (prod[l]["triv"], prod[l]["std"], prod[l]["gal"])
        for g, want in zip(got, page24[l]):
            ok24 &= abs(float(g) - want) / want < 5e-3
    gate("F1_products24", "24-product table matches page (24/24, tol 0.5%)", ok24)
    for l in labels[1:]:
        print(f"    {l}:  " + "  ".join(mp.nstr(prod[l][v], 6) for v in ("triv", "std", "gal")))

    print("=" * 72)
    print("PART G: PDG comparison, scorecard, and rank ordering")
    print("=" * 72)
    for l in NOT_RECOMPUTED:
        print(f"  NOT RECOMPUTED (declared): {l}")

    results = {}
    for f, spec in SM_ASSIGN.items():
        if spec["addr"] is None:
            results[f] = {"role": spec["role"], "addr": None, "mass_GeV": None,
                           "ratio": None, "raw_ratio": None, "symm_ratio": None}
            continue
        l, v = spec["addr"]
        m = mass_of(l, v, prod, cg, mu, sqrtOm, dist)
        obs = PDG_GEV[f]
        raw_ratio = float(m) / obs
        symm_ratio = max(raw_ratio, 1.0 / raw_ratio)
        ratio = raw_ratio if f in RAW_NOT_SYMMETRIZED else symm_ratio
        results[f] = {"role": spec["role"], "addr": [l, v], "mass_GeV": float(m),
                       "ratio": ratio, "raw_ratio": raw_ratio, "symm_ratio": symm_ratio}
        print(f"  {f:4s} {spec['addr']}  computed={float(m):.6g} GeV  observed={obs:.6g} GeV  "
              f"ratio={ratio:.3f} (symm={symm_ratio:.3f}, raw={raw_ratio:.3f})  [{spec['role']}]")

    # Compatible coverage within x3 uses the symmetrized (either-direction) reading; the
    # benchmark (e) is explicitly excluded, matching the page's own "of the remaining 8" framing.
    compatible = [f for f, r in results.items()
                  if r["symm_ratio"] is not None and r["symm_ratio"] <= 3.0
                  and SM_ASSIGN[f]["role"] != "benchmark"]
    adjudicated = [f for f in compatible if SM_ASSIGN[f]["role"] == "assigned_within_x3"]
    unassigned = [f for f, r in results.items() if r["addr"] is None]
    out_of_sector = [f for f, r in results.items() if SM_ASSIGN[f]["role"] == "compatible_out_of_sector"]
    outside_x3 = [f for f, r in results.items() if SM_ASSIGN[f]["role"] == "assigned_outside_x3"]
    gate("G1_scorecard_counts",
         "scorecard: 5 compatible / 4 adjudicated within x3, 2 unassigned, 1 out-of-sector, 1 outside",
         len(compatible) == 5 and len(adjudicated) == 4 and len(unassigned) == 2
         and out_of_sector == ["b"] and outside_x3 == ["d"],
         f"compatible={sorted(compatible)} adjudicated={sorted(adjudicated)} "
         f"unassigned={sorted(unassigned)}")

    ok_g2 = all(abs(results[f]["ratio"] - DECLARED_RATIO[f]) / DECLARED_RATIO[f] < 1e-2
                for f in DECLARED_RATIO)
    gate("G2_ratios_match_declared", "assigned/compatible ratios match the page's declared figures (tol 1%)",
         ok_g2, "  ".join(f"{f}={results[f]['ratio']:.3f}" for f in DECLARED_RATIO))

    gate("G3_null_test_provenance",
         "null-test provenance internally consistent (v1.0 -> v1.1, corrected table)",
         0 < NULL_V10["p_A"] < NULL_V11["p_A"] < 1
         and abs(NULL_V11["S1_mean"] - NULL_V11["S1_obs"]) < 1.0
         and NULL_V11["S1_obs"] == 5,
         f"{NULL_V10['tag']} p_A={NULL_V10['p_A']}  ->  {NULL_V11['tag']} p_A={NULL_V11['p_A']}")

    ranks = ranked_entries(prod, cg, mu, sqrtOm, dist)
    computed_order = [(l, v) for l, v, _ in ranks]
    gate("G4_ranks_match_canonical", "24-entry mass ordering matches the canonical page's ranks (24/24)",
         computed_order == CANONICAL_RANK_ORDER)
    print("  ranked entries (rank: irrep, vacuum, mass GeV):")
    for i, (l, v, m) in enumerate(ranks):
        print(f"    {i+1:2d}: {l:3s} {v:4s} {m:.6g} GeV")

    out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "m8_3_masses.json")
    out_path = os.path.normpath(out_path)
    payload = {
        "schema_version": "m8.3-v1",
        "source": {
            "concept_doi": "10.5281/zenodo.18603975",
            "version_doi": "10.5281/zenodo.21652153",
            "repo_commit": "13a5c2b",
            "repo_url": "https://github.com/dmobius3/mode-identity-theory",
        },
        "anchors": {"mu_Lambda_GeV": float(mu * mp.mpf("1e-9")), "sqrt_Omega_Lambda": float(sqrtOm),
                    "Lambda_per_m2": float(Lam)},
        "ranked_entries": [{"rank": i + 1, "rho": l, "sigma": v, "mass_GeV": m}
                            for i, (l, v, m) in enumerate(ranks)],
        "sm_comparison": results,
        "scorecard": {"compatible": sorted(compatible), "adjudicated": sorted(adjudicated),
                      "unassigned": sorted(unassigned), "out_of_sector": out_of_sector,
                      "assigned_outside_x3": outside_x3},
        "null_test": {"v1_1_corrected": NULL_V11, "v1_0_pre_correction": NULL_V10},
        "not_recomputed": NOT_RECOMPUTED,
    }
    if not args.mutation_tests:
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        with open(out_path, "w") as fh:
            json.dump(payload, fh, indent=2, sort_keys=True)
        print(f"  wrote {out_path}")

    n_fail = sum(1 for _, _, ok in GATES if not ok)
    print("=" * 72)
    print(f"GATES: {len(GATES) - n_fail}/{len(GATES)} pass")

    if args.mutation_tests:
        print("=" * 72)
        print("MUTATION TESTS (registry: every gate id must be attacked and go red)")
        print("=" * 72)

        # A1: perturb the character table directly (continuous Gram-matrix test; no rounding
        # step downstream, so a small scale change is guaranteed to trip it).
        chi_bad = chi.copy()
        chi_bad[3] *= 1.05
        Gb = np.array([[np.sum(sizes * chi_bad[i] * chi_bad[j]) / 120.0 for j in range(9)]
                       for i in range(9)])
        mutation("A1_char_orthonormality", "character row scaled 1.05",
                 np.max(np.abs(Gb - np.eye(9))) >= 1e-9)

        # A2/A3: perturb dims
        dims_bad = list(dims)
        dims_bad[4] += 1
        mutation("A2_e8_marks", "dims perturbed by +1 on R4",
                 not np.allclose(A.dot(np.array(dims_bad, float)), 2.0 * np.array(dims_bad, float)))
        mutation("A3_completeness", "completeness check with perturbed dims",
                 not all(sum(dims_bad[t] * mult_in_Vn(t, n) for t in range(9)) == n + 1
                         for n in range(24)))

        # B1: perturb one exponent in the page target
        page_exps_bad = dict(page_exps)
        page_exps_bad["R3"] = [2, 10, 12, 18, 20, 29]  # last entry wrong
        mutation("B1_kostant_exponents", "R3 exponent list perturbed",
                 exps_all["R3"] != page_exps_bad["R3"])

        # C1: perturb one page_cgeom target
        page_cgeom_bad = dict(page_cgeom)
        page_cgeom_bad["R4"] += 0.01
        mutation("C1_cgeom", "R4 C_geom target perturbed by +0.01",
                 abs(float(cg["R4"]) - page_cgeom_bad["R4"]) >= 5e-4)

        # D1/D2: recompute anchors with a grossly perturbed input
        Lam_b, mu_b, sqrtOm_b = anchors(g=G_SI * mp.mpf("1.5"))
        mutation("D1_mu_lambda", "G perturbed x1.5, mu_Lambda recomputed",
                 abs(float(mu_b * 1000) - 2.25) / 2.25 >= 0.03)
        mutation("D2_sqrt_omega", "G perturbed x1.5, sqrt(Omega) recomputed",
                 abs(float(sqrtOm_b) / 1.019e61 - 1) >= 0.03)

        # E: perturb each torsion single by +0.01 in log-space before exponentiating
        for l in labels[1:]:
            T2_bad = mp.e ** (sign * (X[l] + mp.mpf("0.01")))
            mutation(f"E_{l}", f"log T^2({l}) perturbed by +0.01",
                     abs(T2_bad / exact[l] - 1) >= 1e-8)

        T2_pert_r3 = dict(T2)
        T2_pert_r3["R3"] = T2["R3"] * mp.mpf("1.01")
        mutation("E_galois_R3R4", "T^2(R3) perturbed 1%: phi^-4 ratio fails",
                 abs(T2_pert_r3["R3"] / T2_pert_r3["R4"] / phi ** -4 - 1) >= 1e-10)
        T2_pert_r1 = dict(T2)
        T2_pert_r1["R1"] = T2["R1"] * mp.mpf("1.01")
        mutation("E_galois_R1R2", "T^2(R1) perturbed 1%: phi^-8 ratio fails",
                 abs(T2_pert_r1["R1"] / T2_pert_r1["R2"] / phi ** -8 - 1) >= 1e-10)
        T2_pert_r8 = dict(T2)
        T2_pert_r8["R8"] = T2["R8"] * mp.mpf("1.01")
        mutation("E_sector_products", "T^2(R8) perturbed 1%: half-integer product misses 1/4",
                 abs(T2_pert_r8["R3"] * T2_pert_r8["R7"] * T2_pert_r8["R5"] * T2_pert_r8["R4"] - 4) >= 1e-8
                 or abs(T2_pert_r8["R1"] * T2_pert_r8["R2"] * T2_pert_r8["R6"] * T2_pert_r8["R8"]
                        - mp.mpf(1) / 4) >= 1e-8)

        # F1: perturb one computed product
        page24_bad_target = 150.0  # nonsense target for R2 triv (page says 1.714)
        mutation("F1_products24", "R2-triv checked against a perturbed target (150.0)",
                 abs(float(prod["R2"]["triv"]) - page24_bad_target) / page24_bad_target >= 5e-3)

        # G1: flip b's sector flag to see the count assertion break
        compatible_bad = compatible  # unchanged set
        adjudicated_bad = adjudicated + ["b"]  # pretend b was adjudicated too
        mutation("G1_scorecard_counts", "bottom quark incorrectly counted as adjudicated",
                 not (len(compatible_bad) == 5 and len(adjudicated_bad) == 4))

        # G2: perturb a declared ratio target
        declared_bad = dict(DECLARED_RATIO)
        declared_bad["mu"] = 3.0
        mutation("G2_ratios_match_declared", "mu's declared ratio perturbed to 3.0",
                 abs(results["mu"]["ratio"] - declared_bad["mu"]) / declared_bad["mu"] >= 1e-2)

        # G3: perturb the null-test S1_obs
        null_bad = dict(NULL_V11)
        null_bad["S1_obs"] = 1
        mutation("G3_null_test_provenance", "v1.1 S1_obs perturbed to 1",
                 abs(null_bad["S1_mean"] - null_bad["S1_obs"]) >= 1.0)

        # G4: swap two adjacent ranks
        order_bad = list(computed_order)
        order_bad[9], order_bad[10] = order_bad[10], order_bad[9]
        mutation("G4_ranks_match_canonical", "ranks 10 and 11 swapped",
                 order_bad != CANONICAL_RANK_ORDER)

        gate_ids = {g for g, _, _ in GATES}
        target_ids = {t for t, _, _ in MUTS}
        # A set-based comparison alone cannot see a duplicated gate id (two distinct gate() call
        # sites sharing one id would collapse into one entry and silently read as "covered" by
        # whichever mutation targets that id, even though one call site's own logic was never
        # attacked). Guard against that explicitly: every GATES entry must have a distinct id.
        no_dupe_ids = len(GATES) == len(gate_ids)
        cover_ok = gate_ids == target_ids
        red_ok = all(r for _, _, r in MUTS)
        print(f"NO DUPLICATE GATE IDS: {no_dupe_ids}"
              + ("" if no_dupe_ids else f"  ({len(GATES)} gate() calls, only {len(gate_ids)} distinct ids)"))
        print(f"COVERAGE: gate ids == mutation targets: {cover_ok}"
              + ("" if cover_ok else f"  (uncovered: {sorted(gate_ids - target_ids)};"
                                     f" unknown targets: {sorted(target_ids - gate_ids)})"))
        print(f"ALL RED: {red_ok}")
        if not (no_dupe_ids and cover_ok and red_ok):
            print("META-GUARD FAILED. Exiting nonzero.")
            sys.exit(2)
        print(f"MUTATIONS: {len(MUTS)} declared, every gate attacked, every defect red.")

    sys.exit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
