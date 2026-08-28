"""M5.32 R6.b (symbolic, short): the C3 verdict. A Lorentz-invariant
potential V(M) with no derivatives is ORBIT-BLIND: constant along any
pointwise Lorentz dressing M(x) = Q(x) M_vac Q(x)^T, so no C3 member can
localize (or even see) the boost dressing that R4 found delocalized.

EQUATIONS FIRST
---------------
Setting. M(x) symmetric 4x4 (a (0,2) tensor), eta = diag(-1,1,1,1),
Q in O(1,3):  Q^T eta Q = eta.                                       (1)
Field transformation (the mandate's Lorentz action): M -> Q M Q^T.  (2)

LEMMA (conjugation).  From (1), Q^T eta = eta Q^-1, hence
    (Q M Q^T) eta = Q (M eta) Q^-1.                                  (3)
So the (1,1) tensor N := M eta transforms by CONJUGATION, and every
Lorentz-invariant algebraic function of M is a conjugation invariant
of N restricted to the eta-self-adjoint N (N^T = eta N eta).

THEOREM (orbit blindness). Let V(M) be any function of M alone (no
derivatives) with V(Q M Q^T) = V(M) for all Q in O(1,3), and let
M(x) = Q(x) M_vac Q(x)^T with Q(x) in O(1,3) pointwise. Then
    V(M(x)) = V(M_vac)   for every x,                                (4)
identically in the dressing (any amplitude, any profile, any window).
Proof: (4) is the invariance hypothesis evaluated at Q = Q(x). QED.
No spectral argument is needed for (4): the dressing family is BY
CONSTRUCTION a subset of one orbit, and V is constant on orbits.

COROLLARY (what V can act through). Under (3) the characteristic
polynomial of N = M eta is invariant, so the power sums
    t_p(M) = tr((M eta)^p),  p = 1..4,                               (5)
are invariant; by Newton's identities they generate the ring of
polynomial conjugation invariants of a 4x4 matrix (the coefficients
e_1..e_4 of det(x I - N)), and det(M) = det(N) det(eta)^-1 = -e_4(N)
(det eta = -1) lies in that ring: the det-hedge is not an escape.   (6)
Restricted to the vacuum sector (M symmetric, N diagonalizable over R
with one timelike eigenvector, the case of M_vac = diag(-s g, 1, delta,
0) and of the hedgehog Mb = Qh M_vac Qh^T), the eigenvalues of N are
real and the sorted spectrum is a continuous function of (5), so the
covariant T2 penalty
    T2_cov(M) = sum_i w_i (lambda_i(M eta) - v_i)^2,  sorted,         (7)
with ANY per-eigenvalue weights w_i is spectrum-only, hence orbit-blind.
Subtleties (stated, not needed for (4)):
  (a) O(1,3) is a proper subgroup of GL(4), so its conjugation
      invariants on eta-self-adjoint N are FINER than (5): the Segre
      type is invariant, and within the real-diagonalizable type the
      DISCRETE label "which eigenvalue owns the timelike eigenvector"
      is invariant but not a function of the unordered spectrum
      (diag(a,b,c,d) and diag(-b,-a,c,d) share the spectrum of N and
      lie on different orbits). A V may carry that label (it is what a
      per-eigenvalue weight keyed to the causal character does); it is
      still constant on every orbit, so (4) holds regardless.
  (b) The converse "spectrum determines the orbit" fails on the
      non-diagonalizable (null Jordan) types; irrelevant to (4).
  (c) tr(M^2) = tr(M M) is NOT of the form (5): it contracts with the
      Euclidean delta, is not Lorentz-invariant, and is the control
      that DOES see the dressing (it grows like cosh(2b)).

ESCAPE ROUTES (section 3) and their class:
  (i)  "potentials" containing derivatives (V(M, dM), the velocity-
       coupled V(M, d_0 M) named by the M5.21.15 theorem): these are
       kinetic terms; 2 derivatives = C4, 4 derivatives quartic in
       jets = C6 (non-commutator) or the curvature classes; not C3.
  (ii) V of the FIELD's u-frame data. If u = u(M) is the field's own
       timelike eigenvector, every Lorentz scalar built from (M, u(M))
       is again a Lorentz-invariant function of M alone (u transforms
       covariantly: u -> Q u), so it lies in the ring of section (5)-(6)
       and is orbit-blind; u^T eta M eta u is simply the timelike
       eigenvalue. If instead u is a FIXED external frame, V is not
       Lorentz-invariant: excluded by the mandate.
  (iii) the vacuum-vanishing coefficient c(M) = c0 h(V(M)) (I3, class
       C8): a function of a spectrum-only V, hence itself spectrum-only
       and orbit-blind as a MULTIPLIER; on the eigenvalue-pinned
       dressing families (E_V = 0 pointwise) it is the constant
       c0 h(0) over the whole family, so the multiplied term is seen
       only through its own derivative content.

DERRICK (section 4). Radial boost dressing b(r) dilated x -> mu x
(b_mu(r) = b(r/mu), mu > 1 = WIDER), 3 static dimensions, k first
derivatives per monomial: E_k(mu) = mu^(3-k) E_k(1):
    V (0 derivatives)          mu^3   but E_V == 0 on the orbit (4)
    C4  2-derivative           mu^1
    I1-type curvature^2, C6    mu^-1  (quartic in jets, 4 derivatives)
    C5  curvature^4            mu^-5  (8 derivatives)
Penalizing WIDE = dE/dmu > 0 with E bounded below on mu -> 0:
    C4 with coefficient > 0 (tr(h dM h dM) is PSD in the h-metric form):
        E = c mu, the ONLY bounded localizer of the four.
    C6 / I1-type with coefficient < 0: E = -|c|/mu rises with mu but
        is unbounded below as mu -> 0 (the certified lambda = 0 UV
        amplitude runaway pattern); with coefficient > 0 it decreases
        with mu (never penalizes wide).
    C5 with coefficient < 0: E = -|c|/mu^5, same UV runaway, worse;
        with coefficient > 0 it decreases with mu.
Cross terms with the hedgehog background (dQ ~ amp/R on r < R,
d Mb ~ 1/r): a 2n-derivative monomial with k derivatives on Q and
2n - k on Mb scales as amp^k R^(3-2n) when 2n - k < 3 (core-dominated
R^-k log R or R^-k otherwise), so the pure exponent 3 - 2n survives
the background: C4 still grows like R, quartics still fall like 1/R.
The R4 static gain is read from the audit JSON (no lattice run here)
and its R-exponent fitted, so the ladder can compare it with the mu^1
cost of C4: a gain rising faster than R^1 is not localized by C4.

CHECKS RUN
----------
S1  sympy: (1) for the exact rational boost (cosh, sinh) = (17/8, 15/8)
    and a rational rotation, and for a symbolic-rapidity boost; (3)
    exact; t_p(Q M Q^T) - t_p(M) == 0, p = 1..4, for a fully symbolic
    symmetric M (10 symbols); det(M) = -det(M eta) exact.
N1  50 random sector points M = Q0 D Q0^T (D random real spectrum with
    the vacuum's ordering, Q0 random boost*rotation) dressed by 50
    random pointwise Lorentz Q with rapidity up to 3: max |f(QMQ^T) -
    f(M)| for f in {V4, T2_cov (random weights), det, t_1..t_4} and the
    control tr(M^2) (plus T2_euclid = eigvalsh(M) penalty, non-covariant).
N2  the R4 dressing family (m5_32_r4_audit_clock.LocFamily, n = 32,
    L = 48, hedgehog Mb, radial boost b = amp tanh(r/2) w(r; R)): 50
    random members (amp uniform in (0, 3], R in {6, 9, 12, 18, 24},
    kinds exp2 / exp1 / hard / pow1 / pow2 / pow3): per member the max
    over cells of |f(Md)(x) - f(Mb)(x)|; expected 0 (float) for the
    covariant f, O(cosh(2 amp) g^2) for tr(M^2).
N3  the R4 gain exponent: from data/m5_32_r4_audit_clock_box_n64_L96.json
    the envelope static gain min_amp E_stat(amp, R) - E_stat(0) at
    lambda = 0.75, 1 versus R (6..48), log-log slope.

Output: ../data/m5_32_r6_orbitblind.json. Runtime ~1 min.
"""
from __future__ import annotations

import importlib.util
import json
import os
import time

import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT = os.path.join(DATA, "m5_32_r6_orbitblind.json")
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
G, DELTA, S = 32.0, 0.3, -1.0


# ================= S1: the symbolic proof =================
def stage_symbolic():
    eta = sp.diag(-1, 1, 1, 1)
    syms = sp.symbols("m00 m01 m02 m03 m11 m12 m13 m22 m23 m33")
    M = sp.zeros(4, 4)
    k = 0
    for i in range(4):
        for j in range(i, 4):
            M[i, j] = M[j, i] = syms[k]
            k += 1
    # exact rational boost along x (v = 3/5) and rational rotation (3-4-5) in yz
    ch, sh = sp.Rational(17, 8), sp.Rational(15, 8)
    Bx = sp.Matrix([[ch, sh, 0, 0], [sh, ch, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    c, s = sp.Rational(3, 5), sp.Rational(4, 5)
    Ryz = sp.Matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, c, -s], [0, 0, s, c]])
    Qrat = Bx * Ryz * Bx.T * Ryz.T          # a generic-looking rational element
    # symbolic-rapidity boost along a symbolic unit direction (nx, ny, nz)
    b, nx, ny, nz = sp.symbols("b n_x n_y n_z", real=True)
    K = sp.Matrix([[0, nx, ny, nz], [nx, 0, 0, 0], [ny, 0, 0, 0], [nz, 0, 0, 0]])
    K2 = sp.Matrix([[1, 0, 0, 0], [0, nx * nx, nx * ny, nx * nz],
                    [0, ny * nx, ny * ny, ny * nz], [0, nz * nx, nz * ny, nz * nz]])
    Qb = sp.eye(4) + sp.sinh(b) * K + (sp.cosh(b) - 1) * K2
    unit = {nz: sp.sqrt(1 - nx ** 2 - ny ** 2)}

    def zero(expr):
        e = sp.expand(expr.subs(unit))
        e = sp.simplify(e.rewrite(sp.exp)) if e.has(sp.cosh) or e.has(sp.sinh) else e
        return sp.simplify(e) == 0

    res = {}
    res["eq1_rational"] = all(zero(x) for x in (Qrat.T * eta * Qrat - eta))
    res["eq1_symbolic_boost"] = all(zero(x) for x in (Qb.T * eta * Qb - eta))
    # (3): (Q M Q^T) eta == Q (M eta) Q^-1 for the rational Q (exact)
    lhs = Qrat * M * Qrat.T * eta
    rhs = Qrat * (M * eta) * Qrat.inv()
    res["eq3_rational"] = all(sp.expand(x) == 0 for x in (lhs - rhs))
    # (5): trace invariants, rational Q (exact polynomial identity in 10 symbols)
    N = M * eta
    Nq = Qrat * M * Qrat.T * eta
    P, Pq = sp.eye(4), sp.eye(4)
    tr_ok = []
    for p in range(1, 5):
        P, Pq = P * N, Pq * Nq
        tr_ok.append(sp.expand(Pq.trace() - P.trace()) == 0)
    res["eq5_rational_p1to4"] = tr_ok
    # (5) with the symbolic-rapidity boost, p = 1, 2 (p = 3, 4 by the same
    # lemma; the p <= 2 symbolic check keeps the runtime short)
    Nb = Qb * M * Qb.T * eta
    P, Pb = sp.eye(4), sp.eye(4)
    tr_b = []
    for p in range(1, 3):
        P, Pb = P * N, Pb * Nb
        tr_b.append(zero(Pb.trace() - P.trace()))
    res["eq5_symbolic_boost_p1to2"] = tr_b
    # (6): det(M) = -det(M eta), and det invariant under Qrat
    res["eq6_det_sign"] = sp.expand(M.det() + (M * eta).det()) == 0
    res["eq6_det_invariant_rational"] = sp.expand((Qrat * M * Qrat.T).det() - M.det()) == 0
    # (a): the discrete label example, same spectrum, different orbit
    a_, b_, c_, d_ = sp.symbols("a b c d")
    N1 = sp.diag(a_, b_, c_, d_) * eta
    N2 = sp.diag(-b_, -a_, c_, d_) * eta
    x = sp.symbols("x")
    res["subtlety_a_same_charpoly"] = sp.expand(
        (x * sp.eye(4) - N1).det() - (x * sp.eye(4) - N2).det()) == 0
    # timelike eigenvalue (eigenvector e0) differs: -a vs b
    res["subtlety_a_timelike_eigenvalue"] = {"M1": str(N1[0, 0]), "M2": str(N2[0, 0])}
    log(f"S1 symbolic: {res}")
    return res


# ================= the potentials =================
def traces(M):
    Me = M @ ETA
    P = Me
    t = []
    for p in range(4):
        if p:
            P = P @ Me
        t.append(np.einsum("...kk->...", P))
    return t


C_P = tuple((S * G) ** k + 1.0 + DELTA ** k for k in range(1, 5))
W1 = 0.000724023879
V_VAC = np.sort(np.array([-S * G, 1.0, DELTA, 0.0]) * np.array([-1.0, 1.0, 1.0, 1.0]))
# spectrum of M_vac eta = (s g, 1, delta, 0) sorted


def V4(M):
    t = traces(M)
    return W1 * sum((t[k] - C_P[k]) ** 2 for k in range(4))


def spectrum_cov(M):
    lam = np.linalg.eigvals(M @ ETA)
    return np.sort(lam.real, axis=-1), np.abs(lam.imag).max()


def T2_cov(M, w):
    lam, im = spectrum_cov(M)
    return np.sum(w * (lam - V_VAC) ** 2, axis=-1), im


def T2_euclid(M, w):
    lam = np.linalg.eigvalsh(M)
    return np.sum(w * (lam - np.sort(V_VAC)) ** 2, axis=-1)


def det(M):
    return np.linalg.det(M)


def trM2(M):
    return np.einsum("...ij,...ji->...", M, M)


def boost(rap, nvec):
    K = np.zeros((4, 4))
    K[0, 1:] = nvec
    K[1:, 0] = nvec
    from scipy.linalg import expm
    return expm(rap * K)


def rotation(rng):
    A = rng.normal(size=(3, 3))
    Qr, _ = np.linalg.qr(A)
    if np.linalg.det(Qr) < 0:
        Qr[:, 0] *= -1
    R = np.eye(4)
    R[1:, 1:] = Qr
    return R


def random_lorentz(rng, rap_max):
    nvec = rng.normal(size=3)
    nvec /= np.linalg.norm(nvec)
    return boost(rng.uniform(0, rap_max), nvec) @ rotation(rng)


def conj(Q, M):
    return np.einsum("...ab,...bc,...dc->...ad", Q, M, Q)


def evaluate(M, w):
    t = traces(M)
    t2c, im = T2_cov(M, w)
    return {"V4": V4(M), "T2_cov": t2c, "det": det(M),
            "t1": t[0], "t2": t[1], "t3": t[2], "t4": t[3],
            "T2_euclid_control": T2_euclid(M, w), "trM2_control": trM2(M)}, im


# ================= N1: random sector points =================
def stage_random(rng, npts=50):
    w = rng.uniform(0.5, 2.0, size=4)
    rows = []
    for _ in range(npts):
        # a sector point: random real spectrum near the vacuum's, random frame
        d = np.array([-S * G, 1.0, DELTA, 0.0]) + rng.normal(scale=0.5, size=4)
        M0 = conj(random_lorentz(rng, 1.0), np.diag(d))
        Q = random_lorentz(rng, 3.0)
        f0, im0 = evaluate(M0, w)
        f1, im1 = evaluate(conj(Q, M0), w)
        rows.append({k: abs(f1[k] - f0[k]) / max(1.0, abs(f0[k])) for k in f0}
                    | {"max_imag": max(im0, im1), "rapidity": float(np.arccosh(Q[0, 0]))})
    out = {"npts": npts, "weights_T2": w.tolist(), "rapidity_max": 3.0,
           "max_rel_variation": {k: float(max(r[k] for r in rows)) for k in rows[0]}}
    log(f"N1 random sector: {out['max_rel_variation']}")
    return out


# ================= N2: the R4 dressing family =================
def stage_family(rng, nmem=50):
    R4 = _load("m5_32_r4_audit_clock", "m5_32_r4_audit_clock.py")
    R2 = R4.R2
    B3 = R4.B3
    fam = R4.LocFamily(32, 48.0)
    Mb = fam.Mb
    w = rng.uniform(0.5, 2.0, size=4)
    fb, imb = evaluate(Mb, w)
    kinds = ("exp2", "exp1", "hard", "pow1", "pow2", "pow3")
    Rs = (6.0, 9.0, 12.0, 18.0, 24.0)
    members, worst = [], {}
    for i in range(nmem):
        amp = float(rng.uniform(0.05, 3.0))
        R = float(rng.choice(Rs))
        kind = str(rng.choice(kinds))
        Qb = R2.qb_field(fam.cfg, fam.b_of(amp, R, kind))
        Md = B3.sym4(R2.conj(Qb, Mb))
        fd, imd = evaluate(Md, w)
        var = {k: float(np.max(np.abs(fd[k] - fb[k]))) for k in fd}
        scale = {k: float(np.max(np.abs(fb[k]))) for k in fb}
        rel = {k: var[k] / max(1.0, scale[k]) for k in var}
        members.append({"amp": amp, "R": R, "kind": kind, "max_abs_variation": var,
                        "max_rel_variation": rel, "max_imag_eig": float(max(imb, imd)),
                        "b_max": float(np.max(np.abs(fam.b_of(amp, R, kind))))})
        for k in var:
            worst[k] = max(worst.get(k, 0.0), rel[k])
        if i % 10 == 0:
            log(f"N2 member {i}: amp {amp:.2f} R {R:g} {kind}: V4 {var['V4']:.2e} "
                f"T2cov {var['T2_cov']:.2e} det {var['det']:.2e} trM2 {var['trM2_control']:.2e}")
    # E_V of the undressed hedgehog (eigenvalue-pinned: should be ~0)
    out = {"n": 32, "L": 48.0, "nmem": nmem, "amp_range": [0.05, 3.0], "R_grid": list(Rs),
           "kinds": list(kinds), "weights_T2": w.tolist(),
           "hedgehog_scale": {k: float(np.max(np.abs(fb[k]))) for k in fb},
           "hedgehog_V4_max_cell": float(np.max(fb["V4"])),
           "max_rel_variation_over_members": worst, "members": members}
    log(f"N2 family worst rel variation: {worst}")
    return out


# ================= N3: the R4 gain exponent (data only) =================
def stage_gain_fit():
    fn = os.path.join(DATA, "m5_32_r4_audit_clock_box_n64_L96.json")
    d = json.load(open(fn))
    lad = d["kinds"]["exp2"]["ladder"]
    out = {"source": os.path.basename(fn), "fits": {}}
    for lam in (0.75, 1.0):
        Rs, gains, amps = [], [], []
        for key, v in lad.items():
            R = float(key[2:])
            e = (1 - lam) * np.array(v["E_stat_lam0"]) + lam * np.array(v["E_stat_lam1"])
            i = int(np.argmin(e))
            Rs.append(R)
            gains.append(float(e[i] - e[0]))
            amps.append(float(v["amp"][i]))
        Rs, gains = np.array(Rs), np.array(gains)
        neg = gains < 0
        slope = None
        if neg.sum() >= 2:
            slope = float(np.polyfit(np.log(Rs[neg]), np.log(-gains[neg]), 1)[0])
        out["fits"][f"lam_{lam:g}"] = {"R": Rs.tolist(), "envelope_gain": gains.tolist(),
                                       "amp_star": amps, "loglog_slope_of_minus_gain": slope}
        log(f"N3 lam {lam:g}: R {Rs.tolist()} gain {np.round(gains, 3).tolist()} "
            f"amp* {amps} slope {slope}")
    return out


DERRICK = {
    "dilation": "b_mu(r) = b(r/mu), mu > 1 wider; E_k ~ mu^(3-k), k = first derivatives per monomial",
    "classes": {
        "V (C3)": {"k": 0, "exponent": 3, "on_orbit": "E_V == 0 identically (theorem), no Derrick weight",
                   "penalizes_wide": "never (blind)"},
        "C4 2-derivative (tr(h dM h dM), PSD)": {"k": 2, "exponent": 1,
                                                  "penalizes_wide": "coefficient > 0 (E = c mu, bounded below, the only bounded localizer)"},
        "I1-type curvature^2 / C6 non-commutator quartic": {"k": 4, "exponent": -1,
                                                              "penalizes_wide": "only with coefficient < 0, then unbounded below as mu -> 0 (UV runaway); coefficient > 0 favors wide"},
        "C5 curvature^4": {"k": 8, "exponent": -5,
                            "penalizes_wide": "only with coefficient < 0, UV runaway mu^-5; coefficient > 0 favors wide"},
    },
    "background_cross_terms": "k derivatives on Q (amp/R on r<R), 2n-k on the hedgehog (1/r): amp^k R^(3-2n) when 2n-k<3, else R^-k (log R): the pure exponent survives",
}

ESCAPES = {
    "(i) derivative-containing V (velocity-coupled V(M, d0 M))": "kinetic term: C4 (2 derivatives) or C6 / curvature classes (4+); not C3",
    "(ii) V of the field's u-frame data": "u = u(M) covariant: back in the invariant ring, orbit-blind (u^T eta M eta u = the timelike eigenvalue); fixed external u: not Lorentz-invariant, mandate-excluded",
    "(iii) c(M) = c0 h(V(M)) (I3, C8)": "spectrum-only multiplier, orbit-blind; constant c0 h(0) on every eigenvalue-pinned family; the product is seen only through its derivative factor",
}


def main():
    rng = np.random.default_rng(20260827)
    res = {"task": "M5.32 R6.b orbit-blind potential (C3 verdict)",
           "theorem": "V(Q M Q^T) = V(M) for Q in O(1,3) => V(Q(x) M_vac Q(x)^T) = V(M_vac) pointwise; "
                      "invariants = functions of the char polynomial of M eta (tr((M eta)^p), p=1..4, det(M) = -det(M eta)) "
                      "plus the discrete causal label; all constant on orbits",
           "S1_symbolic": stage_symbolic(),
           "N1_random_sector": stage_random(rng),
           "N2_r4_family": stage_family(rng),
           "N3_r4_gain_exponent": stage_gain_fit(),
           "escape_routes": ESCAPES, "derrick": DERRICK}
    res["verdict"] = ("C3 KILLED as a localizer: every Lorentz-invariant derivative-free V (V4, T2 per-eigenvalue, "
                      "det-hedge, LdG lift, c(M) multipliers) is constant along the boost dressing; localization "
                      "can only come from derivative terms, and by Derrick only a positive 2-derivative C4 term "
                      "penalizes wide dressings while staying bounded")
    res["runtime_s"] = time.time() - T0
    with open(OUT, "w") as f:
        json.dump(res, f, indent=1)
    log(f"wrote {OUT}")


if __name__ == "__main__":
    main()
