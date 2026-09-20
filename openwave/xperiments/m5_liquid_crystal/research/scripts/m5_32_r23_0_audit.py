"""M5.32 R23-0 / R23-3 adversarial audit: an independent attempt to refute the form-level
record (data/m5_32_r23_0_form.json) and the well record (data/m5_32_r23_3_well.json).

EQUATIONS FIRST (all re-derived here, nothing imported from the audited scripts)
---------------
eta = diag(-1, 1, 1, 1), N = M eta, t_p = tr N^p, C_p the vacuum power sums of
(-g, 1, delta, delta), r_p = t_p - C_p.
    V4 = w sum_p r_p^2,   L = sum_p a_p r_p,   V = V4 - c L,
    P(x) = sum_p p a_p x^(p-1) = (x + g)(x - 1)(x - delta).
Gradient on symmetric M:  d t_p / d M = p eta N^(p-1)  (eta N^k is symmetric), so
    dV/dM = eta Q(N),  Q(x) = sum_p p V_p x^(p-1):  critical iff Q(N) = 0 as a MATRIX,
i.e. iff the minimal polynomial of N divides Q. On the vacuum orbit N is diagonalizable with
minimal polynomial (x + g)(x - 1)(x - delta), which gives the one-dimensional family; on a
matrix with the same spectrum but a Jordan block at delta the minimal polynomial is quartic
and L is NOT critical (check A1c).
Hessian at the vacuum in the ten orthonormal coordinates E_k:
    H = 2 w sum_p g_p g_p^T - c sum_p a_p Hess t_p,   g_p[k] = p tr(N^(p-1) E_k eta),
    Hess t_p[k, l] = p sum_{j=0..p-2} tr(N^j E_k eta N^(p-2-j) E_l eta).
Orbit tangent at M0:  T_K = K^T M0 + M0 K,  K = eta Omega, Omega antisymmetric (6 generators).
Halo: eps'' = eps / r^2 + beta r^2 eps has the decaying closed form
    eps = sqrt(r) K_nu(sqrt(beta) r^2 / 2),  nu = sqrt(5) / 4,
(small r: r^(1/2 - 2 nu) = r^((1 - sqrt 5)/2)); an own RK4 inward integration cross-checks it.
Well: <F, F>_X = F_ab X_ac X_bd F_cd = || R F R^T ||_F^2 for X = R^T R positive definite, so
under ANY positive definite X both densities are non-negative for ANY field (check B3c).

Every check is one JSON entry and names the number that would make it fail.

Sectors. N is eta-self-adjoint, so a complex pair can only be born when the timelike
eigenvalue t collides with a spacelike one (opposite Krein sign). V depends on the spectrum
alone, hence for every path from the vacuum to a complex pair
    max V >= min_u max( max_{t in [-g, u]} f(t), fc(u) ),
f(t) the minimum of V over the three spacelike eigenvalues at fixed t, fc(u) the same with the
collision t = s = u imposed. Check A5iid evaluates this bound and an explicit polyline of
block-diagonal matrices that attains it (a real level crossing that hands the timelike label
from -g to delta, ending on the V = 0 twin diag(-delta, ...) of the vacuum).

Usage: python m5_32_r23_0_audit.py      (about 2 to 3 min, single process, 2 threads)
"""

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "2")

import json  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402
import sympy as sp  # noqa: E402
from scipy.linalg import expm  # noqa: E402
from scipy.optimize import brentq, differential_evolution, minimize  # noqa: E402
from scipy.special import kv  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_JSON = os.path.join(HERE, "..", "data", "m5_32_r23_0_audit.json")
EXT = os.path.join(
    HERE,
    "..",
    "..",
    "theory",
    "r22_mikulski_008",
    "reports",
    "008-i1-squared-clock",
    "results",
    "L_ladder",
)
T0 = time.time()
GG, DD = 8.0, 0.3
W1 = 0.000724023879
WW = 25.0 * W1
ALPHA = 1.0 / 137.035999
DPH = 1.0 - (ALPHA / (64.0 * np.pi)) ** 0.25
CS = (3e-5, 1e-4, 3e-4, 1e-3, 3e-3)
ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
RNG = np.random.default_rng(20260919)
CHECKS = []


def add(cid, claim, method, value, expect, fails_if, ok):
    CHECKS.append(
        {
            "id": cid,
            "claim": claim,
            "method": method,
            "value": value,
            "threshold or expectation": expect,
            "fails_if": fails_if,
            "verdict": "PASS" if ok else "FAIL",
        }
    )
    print(f"[{time.time() - T0:6.1f}s] {cid:8s} {'PASS' if ok else 'FAIL'}  {value}", flush=True)


# ---------------------------------------------------------------- basics
def a_vec(g, d):
    """a_p from integrating the cubic: P = x^3 + (g-1-d) x^2 + (d - g - g d) x + g d."""
    return np.array([g * d, (d - g - g * d) / 2.0, (g - 1.0 - d) / 3.0, 0.25])


def basis10():
    out = []
    for i in range(4):
        for j in range(i, 4):
            E = np.zeros((4, 4))
            if i == j:
                E[i, i] = 1.0
            else:
                E[i, j] = E[j, i] = 1.0 / np.sqrt(2.0)
            out.append(E)
    return out


E10 = basis10()
LABELS = ["00", "01", "02", "03", "11", "12", "13", "22", "23", "33"]


def to_mat(v):
    return sum(v[k] * E10[k] for k in range(10))


def to_vec(M):
    return np.array([np.sum(M * E) for E in E10])


def power_sums_of_matrix(M):
    N = M @ ETA
    P, t = np.eye(4), []
    for _ in range(4):
        P = P @ N
        t.append(np.trace(P))
    return np.array(t)


class Pot:
    def __init__(self, g, d, w, c):
        self.g, self.d, self.w, self.c = g, d, w, c
        self.spec0 = np.array([-g, 1.0, d, d])
        self.C = np.array([np.sum(self.spec0**p) for p in range(1, 5)])
        self.a = a_vec(g, d)
        self.M0 = np.diag([g, 1.0, d, d])

    def v_mat(self, M):
        r = power_sums_of_matrix(M) - self.C
        return self.w * r @ r - self.c * self.a @ r

    def v_and_grad_vec(self, v):
        M = to_mat(v)
        N = M @ ETA
        pw = [np.eye(4)]
        for _ in range(4):
            pw.append(pw[-1] @ N)
        r = np.array([np.trace(pw[p]) for p in range(1, 5)]) - self.C
        coef = 2.0 * self.w * r - self.c * self.a
        Gm = sum(coef[p - 1] * p * (ETA @ pw[p - 1]) for p in range(1, 5))
        Gm = 0.5 * (Gm + Gm.T)
        return self.w * r @ r - self.c * self.a @ r, to_vec(Gm)

    def v_spec(self, lam):
        lam = np.asarray(lam, dtype=float)
        r = np.array([np.sum(lam**p) for p in range(1, 5)]) - self.C
        return self.w * r @ r - self.c * self.a @ r

    def v_spec_grad(self, lam):
        lam = np.asarray(lam, dtype=float)
        r = np.array([np.sum(lam**p) for p in range(1, 5)]) - self.C
        coef = 2.0 * self.w * r - self.c * self.a
        return sum(coef[p - 1] * p * lam ** (p - 1) for p in range(1, 5))

    def v_pair(self, l1, l2, u, v):
        """spectrum (l1, l2, u + i v, u - i v)"""
        z = complex(u, v)
        t = np.array([l1**p + l2**p + 2.0 * (z**p).real for p in range(1, 5)])
        r = t - self.C
        return self.w * r @ r - self.c * self.a @ r

    def hessian(self):
        N = self.M0 @ ETA
        pw = [np.linalg.matrix_power(N, k) for k in range(4)]
        A = np.zeros((10, 10))
        B = np.zeros((10, 10))
        for p in range(1, 5):
            gp = np.array([p * np.trace(pw[p - 1] @ E @ ETA) for E in E10])
            A += 2.0 * np.outer(gp, gp)
            for k in range(10):
                for m in range(10):
                    s = 0.0
                    for j in range(p - 1):
                        s += np.trace(pw[j] @ E10[k] @ ETA @ pw[p - 2 - j] @ E10[m] @ ETA)
                    B[k, m] += self.a[p - 1] * p * s
        return A, B


def lorentz(rap=1.0):
    Om = np.zeros((4, 4))
    for i in range(4):
        for j in range(i + 1, 4):
            x = RNG.normal() * (rap if i == 0 else 1.0)
            Om[i, j], Om[j, i] = x, -x
    return expm(ETA @ Om)


# ---------------------------------------------------------------- A1
def check_a1():
    g, d = sp.symbols("g delta", positive=True)
    # (a) the gradient formula d t_p / dM = p eta N^(p-1), exact, generic symmetric M
    m = sp.symbols("m0:10")
    Ms = sp.zeros(4, 4)
    k = 0
    for i in range(4):
        for j in range(i, 4):
            Ms[i, j] = Ms[j, i] = m[k]
            k += 1
    et = sp.diag(-1, 1, 1, 1)
    Ns = Ms * et
    worst = 0
    for p in range(1, 5):
        tp = (Ns**p).trace()
        form = p * et * Ns ** (p - 1)
        k = 0
        for i in range(4):
            for j in range(i, 4):
                want = form[i, j] * (1 if i == j else 2)
                if sp.expand(sp.diff(tp, m[k]) - want) != 0:
                    worst += 1
                k += 1
    # (b) the full 10-gradient at the vacuum for a generic (V_1..V_4): J is 10 x 4
    V = sp.symbols("V1:5")
    N0 = sp.diag(-g, 1, d, d)
    Gm = sum((V[p - 1] * p * et * N0 ** (p - 1) for p in range(1, 5)), sp.zeros(4, 4))
    eqs = [sp.expand(Gm[i, j]) for i in range(4) for j in range(i, 4)]
    J = sp.Matrix([[sp.diff(e, v) for v in V] for e in eqs])
    ns = J.nullspace()
    x = sp.symbols("x")
    cub = sp.Poly(sp.expand((x + g) * (x - 1) * (x - d)), x)
    a_sym = sp.Matrix([cub.coeff_monomial(x ** (p - 1)) / p for p in range(1, 5)])
    prop = len(ns) == 1 and sp.simplify(ns[0] * a_sym[3] - a_sym * ns[0][3]) == sp.zeros(4, 1)
    num_ok = np.allclose(
        [float(v.subs({g: 8, d: sp.Rational(3, 10)})) for v in a_sym], a_vec(8.0, 0.3)
    )
    add(
        "A1a",
        "the first-derivative vectors keeping the vacuum critical on the FULL 10-entry gradient "
        "are one-dimensional, spanned by a_p",
        "sympy: d t_p/dM = p eta N^(p-1) verified entry by entry on a generic symmetric M "
        "(mismatches counted), then the nullspace of the 10 x 4 gradient Jacobian at M0, "
        "symbolic g and delta",
        {
            "gradient_formula_mismatches": worst,
            "jacobian_rank": int(J.rank()),
            "nullspace_dim": len(ns),
            "proportional_to_a_p": bool(prop),
            "a_p_numeric_g8_d0.3": [float(v) for v in a_vec(8.0, 0.3)],
        },
        "mismatches 0, rank 3, nullspace 1, proportional to a_p",
        "an off-diagonal gradient entry survived at the vacuum (rank 4 of J), or the nullspace "
        "had dimension 2 (degenerate pair adding freedom)",
        worst == 0 and J.rank() == 3 and bool(prop) and num_ok,
    )
    # (c) boosted vacua, numeric: rank of the 10 x 4 Jacobian and the null vector
    pot = Pot(GG, DD, WW, 0.0)
    worst_gap, worst_angle, worst_res = np.inf, 0.0, 0.0
    for _ in range(200):
        Lm = lorentz(1.0)
        M = Lm.T @ pot.M0 @ Lm
        N = M @ ETA
        Jn = np.array(
            [
                [p * np.sum((ETA @ np.linalg.matrix_power(N, p - 1)) * E) for p in range(1, 5)]
                for E in E10
            ]
        )
        sv = np.linalg.svd(Jn, compute_uv=False)
        worst_gap = min(worst_gap, sv[2] / sv[3] if sv[3] > 0 else np.inf)
        vt = np.linalg.svd(Jn)[2][-1]
        cosang = abs(vt @ pot.a) / np.linalg.norm(pot.a)
        worst_angle = max(worst_angle, 1.0 - cosang)
        worst_res = max(worst_res, np.linalg.norm(Jn @ pot.a) / np.linalg.norm(Jn))
    add(
        "A1b",
        "the same on Lorentz-transformed vacua M = Lambda^T M0 Lambda (random boosts and "
        "rotations): rank 3, null vector a_p",
        "200 random Lambda = expm(eta Omega); SVD of the numeric 10 x 4 Jacobian",
        {
            "min_ratio_sv3_over_sv4": float(worst_gap),
            "max_1_minus_cos_angle_to_a": float(worst_angle),
            "max_rel_residual_J_a": float(worst_res),
        },
        "sv3/sv4 > 1e6 (rank exactly 3), angle defect < 1e-8",
        "a boosted vacuum had a rank-4 Jacobian (no critical family) or rank 2 (extra freedom)",
        worst_gap > 1e6 and worst_angle < 1e-8,
    )
    # (c) wording: is criticality a property of the SPECTRUM? Jordan block at delta
    e = 0.2
    MJ = np.zeros((4, 4))
    MJ[0, 0], MJ[0, 1], MJ[1, 0], MJ[1, 1] = -DD + e, e, e, DD + e
    MJ[2, 2], MJ[3, 3] = -GG, 1.0
    lam = np.linalg.eigvals(MJ @ ETA)
    potc = Pot(GG, DD, WW, 1.0)
    r = power_sums_of_matrix(MJ) - potc.C
    _, gL = Pot(GG, DD, 0.0, 1.0).v_and_grad_vec(to_vec(MJ))
    add(
        "A1c",
        "WORDING: 'keep the spectrum (-g, 1, delta, delta) a critical point': criticality is a "
        "condition on eigenvalues alone",
        "a real symmetric M whose N has exactly that spectrum but a Jordan block at delta "
        "(the 2 x 2 block [[-delta+e, e], [e, delta+e]] in the (0,1) plane, -g and 1 spacelike); "
        "evaluate the power-sum residuals and |grad L|",
        {
            "eigenvalues_of_N": sorted(float(x.real) for x in lam),
            "max_abs_imag": float(np.abs(lam.imag).max()),
            "max_abs_power_sum_residual": float(np.abs(r).max()),
            "norm_grad_L": float(np.linalg.norm(gL)),
            "closed_form_norm_grad_L (g+delta)(1-delta) 2 e": (GG + DD) * (1 - DD) * 2 * e,
        },
        "if the claim's wording held, |grad L| = 0 at every M with that spectrum",
        "|grad L| > 1e-9 at a matrix with the vacuum spectrum: the condition is on the MINIMAL "
        "POLYNOMIAL (diagonalizable N), not on eigenvalues",
        float(np.linalg.norm(gL)) < 1e-9,
    )


# ---------------------------------------------------------------- A2
def check_a2():
    g, d, e, x = sp.symbols("g delta epsilon x")
    F = sp.integrate(sp.expand((x + g) * (x - 1) * (x - d)), x)  # L = sum_i F(l_i) - same at vac
    Ls = F.subs(x, d + e) + F.subs(x, d - e) - 2 * F.subs(x, d)
    diff = sp.simplify(sp.expand(Ls) - sp.expand((d + g) * (d - 1) * e**2 + e**4 / 2))
    # numeric, from an actual boosted and rotated matrix
    pot = Pot(GG, DD, 0.0, 1.0)
    worst = 0.0
    for eps in (1e-3, 0.05, 0.4, 1.3):
        Lm = lorentz(0.7)
        M = Lm.T @ np.diag([GG, 1.0, DD + eps, DD - eps]) @ Lm
        Lnum = -pot.v_mat(M)
        want = (DD + GG) * (DD - 1) * eps**2 + eps**4 / 2
        worst = max(worst, abs(Lnum - want))
    add(
        "A2",
        "L on the split (-g, 1, delta+eps, delta-eps) is (delta+g)(delta-1) eps^2 + eps^4/2",
        "own route: L = sum_i F(lambda_i) with F the antiderivative of the cubic (sympy), plus "
        "L from the traces of boosted and rotated split matrices, eps up to 1.3",
        {"symbolic_difference": str(diff), "max_abs_error_matrix_route": float(worst)},
        "difference 0; matrix route absolute error < 1e-9 (the traces are of order 4096, so "
        "1e-9 is their roundoff after a boost; a 1e-6 error in either coefficient would show as "
        "1.6e-7 at eps 0.4)",
        "a surviving eps^2 or eps^4 coefficient mismatch, or any odd term",
        diff == 0 and worst < 1e-9,
    )


# ---------------------------------------------------------------- A3, A4
def orbit_tangent(M0):
    rows, names = [], []
    for i in range(4):
        for j in range(i + 1, 4):
            Om = np.zeros((4, 4))
            Om[i, j], Om[j, i] = 1.0, -1.0
            K = ETA @ Om
            rows.append(to_vec(K.T @ M0 + M0 @ K))
            names.append(("boost" if i == 0 else "rot") + f"_{i}{j}")
    return np.array(rows), names


def check_a3():
    pot = Pot(GG, DD, WW, 0.0)
    A, B = pot.hessian()
    # finite-difference cross-check of the analytic Hessian at c = 1e-3
    c = 1e-3
    pc = Pot(GG, DD, WW, c)
    v0 = to_vec(pc.M0)
    Hfd = np.zeros((10, 10))
    for k in range(10):
        dd = []
        for h in (2e-3, 1e-3):
            ek = np.zeros(10)
            ek[k] = h
            dd.append((pc.v_and_grad_vec(v0 + ek)[1] - pc.v_and_grad_vec(v0 - ek)[1]) / (2 * h))
        Hfd[:, k] = (4 * dd[1] - dd[0]) / 3  # Richardson, O(h^4)
    H = WW * A - c * B
    # compare on the soft 9 x 9 part (drop M_00, whose 1.5e5 stiffness swamps the FD error)
    idx = [k for k in range(10) if k != 0]
    fd_err = float(np.abs(Hfd - H)[np.ix_(idx, idx)].max())
    T, names = orbit_tangent(pot.M0)
    svT = np.linalg.svd(T, compute_uv=False)
    rankT = int(np.sum(svT > 1e-12 * svT[0]))
    Q = np.linalg.svd(T.T)[0][:, :rankT]
    norms = {n: float(np.linalg.norm(t)) for n, t in zip(names, T)}
    e23 = np.zeros(10)
    e23[LABELS.index("23")] = 1.0
    res23 = float(np.linalg.norm(e23 - Q @ (Q.T @ e23)))
    rows, ok = [], True
    for cc in (0.0,) + CS:
        Hc = WW * A - cc * B
        ev, evec = np.linalg.eigh(Hc)
        scale = np.abs(ev).max()
        nzero = int(np.sum(np.abs(ev) < 1e-9 * scale))
        kill = float(np.abs(Hc @ Q).max())
        ray23 = float(e23 @ Hc @ e23)
        d22 = np.zeros(10)
        d22[LABELS.index("22")], d22[LABELS.index("33")] = np.sqrt(0.5), -np.sqrt(0.5)
        ray22 = float(d22 @ Hc @ d22)
        off = float(
            np.linalg.norm(Hc @ e23 - ray23 * e23) + np.linalg.norm(Hc @ d22 - ray22 * d22)
        )
        coef = (GG + DD) * (1 - DD) * cc
        nonzero = [float(x) for x in ev if abs(x) >= 1e-9 * scale]
        massive = [x for x in nonzero if not (cc > 0 and abs(x - coef) < 1e-6 * coef)]
        rows.append(
            {
                "c": cc,
                "zero_modes": nzero,
                "max_abs_H_times_orbit_tangent": kill,
                "M23_rayleigh": ray23,
                "M22m33_rayleigh": ray22,
                "doublet_eigenvector_residual": off,
                "closed_form_5.81c": coef,
                "massive": massive,
            }
        )
        ok &= nzero == (5 if cc > 0 else 7) and kill < 1e-9 * scale and off < 1e-9 * scale
        ok &= len(massive) == 3 and min(massive) > 0
        if cc > 0:
            ok &= abs(ray23 / coef - 1) < 1e-9 and abs(ray22 / coef - 1) < 1e-9
    add(
        "A3a",
        "Hessian of V4 - c L at the vacuum: 5 zero modes for c > 0 (7 at c = 0), exact doublet "
        "{M_23, (M_22 - M_33)/sqrt 2} at (g + delta)(1 - delta) c, three positive massive modes",
        "own analytic Hessian (trace formula), cross-checked against Richardson central "
        "differences of "
        "the analytic gradient; eigen-decomposition; M_23 and (M_22-M_33)/sqrt2 tested as exact "
        "eigenvectors",
        {"fd_cross_check_max_abs_error_soft_block": fd_err, "rows": rows},
        "zero modes 5 (7 at c 0); doublet rel error < 1e-9; 3 massive > 0; FD error < 1e-6",
        "a sixth zero mode for c > 0, a doublet eigenvalue off 5.81 c, or a massive mode <= 0",
        ok and fd_err < 1e-6,
    )
    add(
        "A3b",
        "the 5 zero modes are the Lorentz orbit (3 boosts + 2 tilts); M_23 is NOT an orbit "
        "direction; the collaborator's 'six zero modes, M_23 a stabilizer direction' is wrong",
        "tangent space of {Lambda^T M0 Lambda} from the six generators K = eta Omega, "
        "T_K = K^T M0 + M0 K; rank by SVD; the norm of each tangent vector; the residual of the "
        "unit M_23 vector after projection on the tangent space",
        {
            "tangent_rank": rankT,
            "tangent_vector_norms": norms,
            "M23_residual_after_projection_on_tangent": res23,
            "reading": "rot_23 generates the ZERO tangent vector (it is the stabilizer of the "
            "degenerate pair): a stabilizer removes a zero mode, it does not add one. The M_23 "
            "entry is the split doublet partner (a 45 degree rotation of (M_22 - M_33)/sqrt 2) "
            "and carries 5.81 c. Record right (5), collaborator wrong (6).",
        },
        "rank 5, |T(rot_23)| = 0, M_23 residual 1 (orthogonal to the orbit)",
        "tangent rank 6, or M_23 inside the tangent span (residual < 1)",
        rankT == 5 and norms["rot_23"] < 1e-14 and abs(res23 - 1.0) < 1e-12,
    )


def c_crit_of(d, w):
    pot = Pot(GG, d, w, 0.0)
    A, B = pot.hessian()
    T, _ = orbit_tangent(pot.M0)
    U = np.linalg.svd(T.T)[0]
    Qp = U[:, 5:]  # orthogonal complement of the orbit

    def low(c):
        return np.linalg.eigvalsh(Qp.T @ (w * A - c * B) @ Qp)[0]

    ev0 = np.linalg.eigvalsh(Qp.T @ (w * A) @ Qp)
    soft = float(min(x for x in ev0 if x > 1e-9 * ev0.max()))
    lo, hi = 1e-12, 1e-12
    while low(hi) >= 0:
        lo, hi = hi, hi * 1.5
    # low(c) is 5.81 c > 0 just above zero, so the first sign change is the massive crossing
    cc = brentq(low, lo, hi, xtol=1e-16, rtol=1e-12)
    ev, evec = np.linalg.eigh(Qp.T @ (w * A - cc * B) @ Qp)
    mode = Qp @ evec[:, 0]
    comp = {LABELS[k]: round(float(mode[k]), 4) for k in range(10) if abs(mode[k]) > 1e-6}
    return cc, soft, comp, (lambda c: float(np.linalg.eigvalsh(w * A - c * B)[0]))


def check_a4():
    kr = (0.7 / (1 - DPH)) ** 4
    c1, s1, m1, _ = c_crit_of(DD, WW)
    c2, s2, m2, _ = c_crit_of(DPH, WW)
    c3, s3, m3, low3 = c_crit_of(DPH, WW / kr)
    lowB = low3(3.094e-5)
    val = {
        "delta_phys": DPH,
        "w_ratio": kr,
        "c_crit_d0.3_w25": c1,
        "c_crit_over_3e-3": c1 / 3e-3,
        "c_crit_dphys_w25": c2,
        "c_crit_dphys_scaled_w": c3,
        "softest_massive_c0": {"d0.3": s1, "dphys_w25": s2, "dphys_scaled": s3},
        "lowest_eig_rowB_at_c_3.094e-05": lowB,
        "destabilizing_mode_components": {"d0.3": m1, "dphys": m2},
    }
    ok = (
        abs(c1 / 0.08633 - 1) < 5e-3
        and abs(c2 / 6.902e-3 - 1) < 5e-3
        and abs(c3 / 1.0433e-6 - 1) < 5e-3
        and abs(s3 / 1.5702e-7 - 1) < 1e-3
        and abs(lowB / -6.3229e-6 - 1) < 1e-3
        and abs(kr / 6615 - 1) < 1e-3
    )
    add(
        "A4",
        "c_crit about 0.086 (29 x 3e-3) at delta 0.3; about 6.9e-3 at delta phys with w 25 W1; "
        "with w / 6615 softest massive 1.57e-7, unstable at c 3.094e-5 (lowest -6.3e-6), stable "
        "below about 1e-6",
        "own Hessian projected on the complement of the orbit tangent; brentq root of the lowest "
        "projected eigenvalue (not a grid)",
        val,
        "each number within 0.5 percent of the record (the record used a 0.2 percent c grid)",
        "c_crit off by more than 0.5 percent, or the row B eigenvalue not negative at 3.094e-5",
        ok,
    )


# ---------------------------------------------------------------- A5
def star_state(pot):
    """the unconstrained minimum r = c a / 2w realized as a symmetric matrix.

    The pair enters through s = z + zbar, q = z zbar (regular at the double root, where the
    (u, v) chart is singular): p1 = s, p2 = s^2 - 2q, p3 = s^3 - 3 s q, p4 = s^4 - 4 s^2 q + 2 q^2.
    """
    from scipy.optimize import fsolve

    target = pot.C + pot.c * pot.a / (2 * pot.w)

    def eq(y):
        l1, l2, sm, q = y
        pair = [sm, sm**2 - 2 * q, sm**3 - 3 * sm * q, sm**4 - 4 * sm**2 * q + 2 * q**2]
        return [l1**p + l2**p + pair[p - 1] - target[p - 1] for p in range(1, 5)]

    y = fsolve(eq, [-GG, 1.0, 2 * pot.d, pot.d**2], xtol=1e-13)
    l1, l2, sm, q = y
    u, v2 = sm / 2, q - sm * sm / 4
    v = np.sqrt(v2) if v2 > 0 else float("nan")
    M = np.zeros((4, 4))
    M[0, 0], M[0, 3], M[3, 0], M[3, 3] = -u, v, v, u  # N block [[u, v], [-v, u]]
    M[1, 1], M[2, 2] = l2, l1  # slot 1 keeps 1, slot 2 (delta) takes -g
    return M, np.array([l1, l2, u, v]), float(np.abs(eq(y)).max())


def check_a5():
    # (i) real spectra
    rows, ok = {}, True
    for name, d, w, clist in (
        ("d0.3", DD, WW, CS),
        ("rowA", DPH, WW, (3.094e-5,)),
    ):
        for c in clist:
            pot = Pot(GG, d, w, c)
            best = np.inf
            for _ in range(250):
                l0 = RNG.uniform(-12, 12, size=4)
                res = minimize(
                    pot.v_spec, l0, jac=pot.v_spec_grad, method="BFGS", options={"gtol": 1e-13}
                )
                best = min(best, float(res.fun))
            de = differential_evolution(
                pot.v_spec,
                [(-12, 12)] * 4,
                seed=int(RNG.integers(1 << 30)),
                tol=1e-14,
                maxiter=300,
                popsize=25,
                polish=True,
                workers=1,
            )
            best = min(best, float(de.fun))
            # a fine local scan of the soft plane: split eps and the shift of the pair
            ee, ss = np.meshgrid(np.linspace(0, 1.5, 301), np.linspace(-0.5, 0.5, 201))
            loc = min(
                pot.v_spec([-GG, 1.0, d + s + e_, d + s - e_])
                for e_, s in zip(ee.ravel()[::7], ss.ravel()[::7])
            )
            rows[f"{name} c={c:g}"] = {"lowest_found": best, "soft_plane_scan_min": float(loc)}
            ok &= best > -1e-11 and loc > -1e-11
    add(
        "A5i",
        "on real spectra the vacuum is the global minimum of V4 - c L for every c used: no V < 0",
        "own search per c: 250 BFGS starts with analytic gradient on [-12, 12]^4, a differential "
        "evolution run, and a scan of the soft (split, pair-shift) plane",
        rows,
        "lowest value found > -1e-11",
        "any real spectrum with V < -1e-11",
        ok,
    )
    # (ii-a) closed form and the realizing matrix
    c = 3e-3
    pot = Pot(GG, DD, WW, c)
    closed = -c * c * float(pot.a @ pot.a) / (4 * WW)
    Ms, y, resid = star_state(pot)
    lam, vec = np.linalg.eig(Ms @ ETA)
    kre = []
    for i in range(4):
        if abs(lam[i].imag) < 1e-12:
            x = vec[:, i].real
            kre.append((float(lam[i].real), float(x @ ETA @ x)))
    pr = [i for i in range(4) if lam[i].imag > 1e-12][0]
    xr, xi = vec[:, pr].real, vec[:, pr].imag
    gram = np.array([[xr @ ETA @ xr, xr @ ETA @ xi], [xr @ ETA @ xi, xi @ ETA @ xi]])
    pair_tab = {}
    for cc in CS:
        yy = star_state(Pot(GG, DD, WW, cc))[1]
        pair_tab[f"{cc:g}"] = [float(yy[2]), float(yy[3]), float(abs(complex(yy[2], yy[3]) - DD))]
    add(
        "A5iia",
        "full symmetric matrices reach V = -c^2 sum a_p^2 / (4 w) (about -501 c^2), on a complex "
        "pair next to delta whose 2-plane has signature (1,1); -g and 1 sit on spacelike vectors",
        "solve r = c a / 2w for the spectrum (l1, l2, u +- i v), build the symmetric M with the "
        "block [[-u, v], [v, u]] in the (0,3) plane, evaluate V from the matrix, read the "
        "eta-norm "
        "of the real eigenvectors and the eta-Gram determinant of the pair's plane",
        {
            "sum_a2_over_4w": float(pot.a @ pot.a) / (4 * WW),
            "closed_form_c3e-3": closed,
            "V_of_matrix": float(pot.v_mat(Ms)),
            "spectrum_l1_l2_u_v": [float(x) for x in y],
            "solve_residual": resid,
            "real_eigenvalues_with_eta_norm": kre,
            "pair_plane_gram_det": float(np.linalg.det(gram)),
            "dist_from_M0_frobenius": float(np.linalg.norm(Ms - pot.M0)),
            "pair_u_v_and_distance_from_delta_by_c": pair_tab,
        },
        "V_of_matrix = closed form (rel 1e-8); Gram det < 0 (signature (1,1)); both real "
        "eigenvectors spacelike; coefficient 501.5",
        "V of the matrix differing from the closed form, or a positive Gram determinant",
        abs(pot.v_mat(Ms) / closed - 1) < 1e-6
        and np.linalg.det(gram) < 0
        and all(k[1] > 0 for k in kre)
        and abs(float(pot.a @ pot.a) / (4 * WW) - 501.54) < 0.1,
    )
    # (ii-b) can a perturbation of M0 (with M_0i) make the pair complex? nearest complex spectrum
    M0 = pot.M0
    worst_im, n_s = 0.0, 0
    for amp in (1e-3, 1e-2, 0.1, 0.5, 1.0, 2.0, 3.0, 4.0):
        v = RNG.normal(size=(4000, 10))
        v *= amp / np.linalg.norm(v, axis=1)[:, None]
        Mb = M0[None] + np.einsum("nk,kab->nab", v, np.array(E10))
        im = np.abs(np.linalg.eigvals(Mb @ ETA).imag).max()
        worst_im = max(worst_im, float(im))
        n_s += 4000

    def first_complex(u):
        u = u / np.linalg.norm(u)
        D = to_mat(u)
        s_lo, s_hi = 0.0, None
        for s in np.linspace(0.5, 40.0, 80):
            if np.abs(np.linalg.eigvals((M0 + s * D) @ ETA).imag).max() > 1e-9:
                s_hi = s
                break
            s_lo = s
        if s_hi is None:
            return 40.0
        for _ in range(40):
            sm = 0.5 * (s_lo + s_hi)
            if np.abs(np.linalg.eigvals((M0 + sm * D) @ ETA).imag).max() > 1e-9:
                s_hi = sm
            else:
                s_lo = sm
        return s_hi

    best_s, best_u = np.inf, None
    for k in range(24):
        u0 = RNG.normal(size=10)
        res = minimize(first_complex, u0, method="Nelder-Mead", options={"maxfev": 400})
        if res.fun < best_s:
            best_s, best_u = float(res.fun), res.x / np.linalg.norm(res.x)
    Mn = M0 + best_s * to_mat(best_u)
    add(
        "A5iib",
        "a real symmetric M near M0 (small perturbation, M_0i included) cannot give the pair "
        "(delta, delta) a complex split",
        "32000 random perturbations of Frobenius size 1e-3 to 4 (all ten entries): max |Im "
        "lambda|; then the nearest complex spectrum: Nelder-Mead over directions of the first "
        "scale s at which M0 + s D turns complex (24 starts)",
        {
            "max_abs_imag_over_random_perturbations_up_to_size_4": worst_im,
            "samples": n_s,
            "nearest_complex_onset_distance": best_s,
            "closed_form_(g+delta)/2": (GG + DD) / 2,
            "V4_at_that_onset": float(Pot(GG, DD, WW, 0.0).v_mat(Mn)),
            "direction_components": {
                LABELS[k]: round(float(best_u[k]), 3) for k in range(10) if abs(best_u[k]) > 0.02
            },
            "reading": "complex pairs need a Krein collision of the timelike eigenvalue with a "
            "spacelike one; the cheapest in Frobenius distance is -g meeting delta halfway",
        },
        "no complex eigenvalue below distance 4; onset at (g + delta)/2 = 4.15",
        "any complex eigenvalue at perturbation size <= 4, or an onset distance below 4.1",
        worst_im < 1e-9 and best_s > 4.1,
    )
    return pot, Ms


def check_a5_paths(pot, Ms):
    c, M0, v0 = pot.c, pot.M0, to_vec(pot.M0)
    closed = -c * c * float(pot.a @ pot.a) / (4 * WW)
    y_star = star_state(pot)[1]
    perms = []
    for slot_g in (1, 2):
        for slot_pair in (2, 3):
            if slot_pair == slot_g:
                continue
            Mp = np.zeros((4, 4))
            l1, l2, u, vv = y_star
            Mp[0, 0], Mp[0, slot_pair], Mp[slot_pair, 0], Mp[slot_pair, slot_pair] = -u, vv, vv, u
            rest = [q for q in (1, 2, 3) if q != slot_pair]
            Mp[slot_g, slot_g] = l1
            other = [q for q in rest if q != slot_g][0]
            Mp[other, other] = l2
            perms.append(Mp)
    starts = []
    for Mp in perms:
        starts.append(Mp)
        for _ in range(12):
            Lm = lorentz(0.6)
            starts.append(Lm.T @ Mp @ Lm)

    def near(v):
        dv = v - v0
        return dv @ dv, 2 * dv

    def nearest(margin):
        lim = margin * closed
        cons = {
            "type": "ineq",
            "fun": lambda v: lim - pot.v_and_grad_vec(v)[0],
            "jac": lambda v: -pot.v_and_grad_vec(v)[1],
        }
        bd, bM = np.inf, None
        for S in starts:
            # the start is nudged off the exact minimum, where the constraint gradient vanishes
            res = minimize(
                near,
                to_vec(S) + RNG.normal(size=10) * 1e-4,
                jac=True,
                constraints=[cons],
                method="SLSQP",
                options={"maxiter": 400, "ftol": 1e-14},
            )
            vfin = pot.v_and_grad_vec(res.x)[0]
            if vfin <= 0.5 * lim and res.fun < bd:
                bd, bM = float(res.fun), to_mat(res.x)
        return float(np.sqrt(bd)), bM

    R_first, best_M = nearest(0.2)
    R_edge, _ = nearest(0.01)
    eta_v = 0.2 * closed
    # Frobenius distance from M0 to the orbit of the V = 0 twin diag(-delta, -g, 1, delta)
    twin = np.diag([-DD, DD, 1.0, -GG])

    def twin_dist(q):
        Om = np.zeros((4, 4))
        k = 0
        for i in range(4):
            for j in range(i + 1, 4):
                Om[i, j], Om[j, i] = q[k], -q[k]
                k += 1
        Lm = expm(ETA @ Om)
        return np.sum((Lm.T @ twin @ Lm - M0) ** 2)

    d_twin = min(
        float(minimize(twin_dist, RNG.normal(size=6) * 0.3, method="BFGS").fun) for _ in range(12)
    )
    lamn = np.linalg.eigvals(best_M @ ETA)
    # straight path barrier
    ss = np.linspace(0, 1, 4001)
    vs = np.array([pot.v_mat(M0 + s * (best_M - M0)) for s in ss])
    # ball search: min V inside ||M - M0|| <= R
    ball = {}
    for R in (1.0, 2.0, 4.0, 6.0, 0.9 * R_first, 1.05 * R_first):
        con = {
            "type": "ineq",
            "fun": lambda v, R=R: R * R - (v - v0) @ (v - v0),
            "jac": lambda v, R=R: -2 * (v - v0),
        }
        lowest = np.inf
        seeds = [v0 + RNG.normal(size=10) * R / 4 for _ in range(30)]
        dirn = to_vec(best_M) - v0
        seeds += [v0 + dirn * min(1.0, 0.98 * R / np.linalg.norm(dirn))]
        for sd in seeds:
            res = minimize(
                lambda v: pot.v_and_grad_vec(v),
                sd,
                jac=True,
                constraints=[con],
                method="SLSQP",
                options={"maxiter": 300, "ftol": 1e-15},
            )
            if (res.x - v0) @ (res.x - v0) <= R * R * (1 + 1e-9):
                lowest = min(lowest, float(res.fun))
        ball[f"R={R:.3f}"] = lowest
    below = [k for k, x in ball.items() if x < 0.5 * eta_v]
    ok_ball = all(x > -1e-9 for k, x in ball.items() if float(k[2:]) < 0.95 * R_first) and any(
        float(k[2:]) > R_first for k in below
    )
    add(
        "A5iic",
        "V < 0 is not reachable from the vacuum without a large barrier (the record: 'a V4 "
        "barrier of order w * 8^8')",
        "c = 3e-3. (1) SLSQP: minimize |M - M0|_F subject to V <= 0.2 V_min (and 0.01 V_min), "
        "52 starts (the "
        "closed-form state in 4 slot assignments, each with 12 random Lorentz transforms); "
        "(2) V along the straight line M0 -> that state; (3) min V inside balls |M - M0| <= R",
        {
            "first_R_with_V_below_0.2_Vmin": R_first,
            "first_R_with_V_below_0.01_Vmin": R_edge,
            "distance_M0_to_the_orbit_of_the_V0_twin": float(np.sqrt(d_twin)),
            "spectrum_of_that_state": [[float(x.real), float(x.imag)] for x in lamn],
            "straight_path_barrier_max_V": float(vs.max()),
            "straight_path_barrier_over_abs_Vmin": float(vs.max() / abs(closed)),
            "record_order_w_8^8": WW * 8.0**8,
            "min_V_in_ball": ball,
        },
        "no V < 0 in any ball of radius below 0.95 R_first; V < 0 found just above R_first; "
        "the straight-path barrier is many orders above |V_min| = 4.5e-3",
        "V < -1e-9 inside a ball of radius below 0.95 R_first (a nearer state than the search "
        "found), or a straight-path barrier below 1e3 |V_min|",
        ok_ball and vs.max() > 1e3 * abs(closed),
    )

    # the barrier of ANY path (Krein collision argument) and an explicit path that attains it
    p0 = pot
    ts = np.concatenate(
        [
            [-GG],
            -GG + np.geomspace(1e-6, 1.0, 900)[:-1],
            np.arange(-7.0, -6.0, 0.0005),
            np.arange(-6.0, DD + 1e-9, 0.002),
        ]
    )
    ts[-1] = DD
    branch, fvals = [], []
    x = np.array([1.0, DD, DD])
    for k, t in enumerate(ts):
        # a small symmetry-breaking nudge lets the continuation leave the l2 = l3 branch at its
        # pitchfork (near t = -g + 0.003) instead of riding the saddle
        x0 = x + (np.array([0.0, 1e-3, -1e-3]) if t < -7.9 else 0.0)
        res = minimize(lambda q: p0.v_spec([t, *q]), x0, method="BFGS", options={"gtol": 1e-12})
        x = res.x
        branch.append(x.copy())
        fvals.append(float(res.fun))
    branch, fvals = np.array(branch), np.array(fvals)
    # the LOWER bound needs the GLOBAL f(t): multi-start on a sub-grid (every 8th knot between
    # -7 and -6, every 40th elsewhere); the continuation value is kept only as one more start
    sub = [k for k in range(len(ts)) if (k % 8 == 0 and -7.0 <= ts[k] <= -6.0) or k % 40 == 0]
    f_glob, glob_gap, gap_near_peak = {}, 0.0, 0.0
    for k in sub:
        bestg = min(
            float(minimize(lambda q: p0.v_spec([ts[k], *q]), st, method="BFGS").fun)
            for st in [RNG.uniform(-9, 5, 3) for _ in range(8)] + [[-5.0, 0.0, 3.0]]
        )
        f_glob[k] = min(bestg, fvals[k])
        glob_gap = max(glob_gap, fvals[k] - f_glob[k])
        if -7.0 <= ts[k] <= -6.0:
            gap_near_peak = max(gap_near_peak, fvals[k] - f_glob[k])
    jumps = np.abs(np.diff(branch, axis=0)).max(axis=1)
    k_jump = int(np.argmax(jumps)) + 1
    gap = np.abs(branch - ts[:, None]).min(axis=1)
    k_cross = int(np.argmin(gap[: int(np.searchsorted(ts, -1.0))]))
    kmax = int(np.argmax(fvals))
    # fc(u) near the crossing: the collision constraint t = s = u
    fcs = {}
    for k in range(max(0, k_cross - 400), min(len(ts), k_cross + 400), 8):
        u = ts[k]
        rest = [branch[k][j] for j in range(3) if j != int(np.argmin(np.abs(branch[k] - u)))]
        fcs[k] = min(
            float(
                minimize(
                    lambda q: p0.v_spec([u, u, *q]), st, method="BFGS", options={"gtol": 1e-12}
                ).fun
            )
            for st in [rest, [0.0, 4.0], [1.0, DD]] + [RNG.uniform(-9, 5, 2) for _ in range(4)]
        )
    sub_k = np.array(sorted(f_glob))
    run_glob = np.maximum.accumulate(np.array([f_glob[k] for k in sub_k]))

    def run_at(k):
        j = int(np.searchsorted(sub_k, k, side="right")) - 1
        return run_glob[j] if j >= 0 else 0.0

    lower = min(max(run_at(k), fcs[k]) for k in fcs)
    fc_below = min(
        float(minimize(lambda q, u=u: p0.v_spec([u, u, *q]), [1.0, DD], method="BFGS").fun)
        for u in (-8.0, -8.5, -9.0, -10.0)
    )
    # explicit path, leg 1: block-diagonal matrices diag(-t, l1, l2, l3) along the branch
    # V is evaluated ALONG the straight segments between consecutive knots (sub-sampled by the
    # knot distance), so the polyline is a continuous path whatever the knot spacing
    leg1 = [np.diag([-t, *b]) for t, b in zip(ts, branch)]
    v_leg1, n_eval = [], 0
    for Ma, Mb in zip(leg1[:-1], leg1[1:]):
        nsub = max(3, int(np.ceil(np.linalg.norm(Mb - Ma) / 0.004)))
        for q in np.linspace(0.0, 1.0, nsub, endpoint=False):
            v_leg1.append(pot.v_mat(Ma + q * (Mb - Ma)))
        n_eval += nsub
    v_leg1 = np.array(v_leg1 + [pot.v_mat(leg1[-1])])
    end = leg1[-1]
    # leg 2: from the twin to the V < 0 minimum, the pair in the (0, slot of delta) plane
    sl_d = 1 + int(np.argmin(np.abs(branch[-1] - DD)))
    sl_g = 1 + int(np.argmin(branch[-1]))
    sl_1 = [q for q in (1, 2, 3) if q not in (sl_d, sl_g)][0]
    l1, l2, u, vv = y_star
    Mst = np.zeros((4, 4))
    Mst[0, 0], Mst[0, sl_d], Mst[sl_d, 0], Mst[sl_d, sl_d] = -u, vv, vv, u
    Mst[sl_g, sl_g], Mst[sl_1, sl_1] = l1, l2
    v_leg2 = np.array([pot.v_mat(end + q * (Mst - end)) for q in np.linspace(0, 1, 2001)])
    upper = float(max(v_leg1.max(), v_leg2.max()))
    w88 = WW * 8.0**8
    add(
        "A5iid",
        "WORDING: the V < 0 sector 'cannot be reached without moving the timelike eigenvalue "
        "from -g to about delta (a V4 barrier of order w * 8^8)' = 3.0e5",
        "lower bound for EVERY path: a complex pair needs the timelike eigenvalue t to collide "
        "with a spacelike one at some u with a real spectrum up to there, so max V >= min_u max( "
        "max_{t <= u} f(t), fc(u) ), f(t) = min over the three spacelike eigenvalues "
        "(continuation "
        "in t, geometric steps off the vacuum, 5e-4 near the crossing, global multi-start on a "
        "sub-grid), fc(u) the same with the collision "
        "imposed. Upper bound: an explicit continuous path of MATRICES, leg 1 block-diagonal "
        "diag(-t, l(t)) along the branch through a REAL crossing t = l1 to the V = 0 twin, leg 2 "
        "straight from the twin to the V < 0 minimum; V evaluated from the matrices",
        {
            "barrier_lower_bound": float(lower),
            "barrier_upper_bound_explicit_path": upper,
            "upper_bound_in_units_of_w": upper / WW,
            "record_w_8^8": w88,
            "record_over_actual": w88 / upper,
            "f_max_at_t": float(ts[kmax]),
            "real_crossing_at_t": float(ts[k_cross]),
            "spacelike_eigenvalues_at_crossing": [float(q) for q in branch[k_cross]],
            "largest_knot_to_knot_move": float(jumps.max()),
            "largest_knot_to_knot_move_at_t": float(ts[k_jump]),
            "polyline_V_evaluations_leg1": n_eval,
            "continuation_minus_global_f_max_gap_anywhere": float(glob_gap),
            "continuation_minus_global_f_max_gap_between_-7_and_-6": float(gap_near_peak),
            "min_fc_for_collisions_at_or_below_-g": fc_below,
            "max_V_on_leg2_twin_to_minimum": float(v_leg2.max()),
            "V_at_leg2_end": float(v_leg2[-1]),
            "leg1_end_matrix_diag": [float(q) for q in np.diag(end)],
            "leg1_max_abs_M0i": 0.0,
            "reading": "the timelike label is carried from -g to delta by a level crossing with a "
            "spacelike eigenvalue that comes down to meet it while a third goes to about +4; the "
            "multiset never has to give up its fourth power sum, so the cost is of order 1.5e3 w, "
            "not 8^8 w. Leg 1 keeps M_0i = 0 exactly (it lies inside the instrument's sector).",
        },
        "if the record's order of magnitude were right, no such path would exist: the explicit "
        "upper bound would sit within a factor 10 of 3.0e5",
        "an explicit continuous path (a polyline of matrices, V sampled every 0.004 along it) "
        "with max V more than 10 x below w 8^8",
        w88 / upper < 10.0,
    )
    add(
        "A5iie",
        "the qualitative reading: the V < 0 sector cannot be reached from the vacuum sector at "
        "small V (a barrier many orders above |V_min|)",
        "the path-independent lower bound of A5iid against |V_min| at the largest c, and against "
        "the largest V_0 of the record's conversion table (0.110, zero spatial block)",
        {
            "barrier_lower_bound": float(lower),
            "abs_Vmin_c3e-3": abs(closed),
            "ratio_to_abs_Vmin": float(lower / abs(closed)),
            "ratio_to_V0_zero_block_0.1103": float(lower / 0.11034084891072946),
        },
        "ratio to |V_min| > 1e3",
        "a collision configuration with V below 1e3 |V_min| reachable with a real spectrum from "
        "the vacuum",
        float(lower / abs(closed)) > 1e3,
    )

    # Frobenius distance to the ORBIT is ill-posed (non-compact orbit): demonstration + search
    k4 = np.array([1.0, 1.0, 0.0, 0.0])
    Apar = -3.0
    Mrest = M0 + Apar * np.outer(k4, k4)
    chi = 9.0
    Om = np.zeros((4, 4))
    Om[0, 1], Om[1, 0] = chi, -chi
    Lb = expm(ETA @ Om)
    for sign in (1, -1):
        Lb = expm(sign * (ETA @ Om))
        dist = np.linalg.norm(Lb.T @ (Mrest - M0) @ Lb)
        if dist < 1e-3:
            break
    lam_b = np.linalg.eigvals(Mrest @ ETA)

    def v_null(q):
        th, ph, A_, b1, b2 = q
        n = np.array([np.cos(th), np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph)])
        e1 = np.array([-np.sin(th), np.cos(th) * np.cos(ph), np.cos(th) * np.sin(ph)])
        e2 = np.cross(n, e1)
        kk = np.concatenate([[1.0], n])
        b = np.concatenate([[0.0], b1 * e1 + b2 * e2])
        return pot.v_mat(M0 + A_ * np.outer(kk, kk) + np.outer(kk, b) + np.outer(b, kk))

    lowest = np.inf
    for _ in range(150):
        q0 = np.concatenate(
            [RNG.uniform(0, np.pi, 2), RNG.normal(size=3) * RNG.choice([0.3, 3, 10])]
        )
        res = minimize(v_null, q0, method="BFGS", options={"gtol": 1e-12})
        lowest = min(lowest, float(res.fun))
    add(
        "A5iif",
        "the requested 'distance R to the vacuum orbit' search is well posed and V < 0 stays at "
        "finite distance from the orbit",
        "(1) demonstration: the boosted vacuum plus e^(-2 chi) A k k^T (k null, A = -3, chi = 9) "
        "is at Frobenius distance ~1e-7 of the orbit yet has a complex pair and V of order 50; "
        "(2) every state at ZERO inf-distance from the orbit through contracting components is "
        "M0 + A k k^T + k b^T + b k^T in some frame: minimize V over (n, A, b), 150 BFGS starts",
        {
            "demo_frobenius_distance_to_orbit_point": float(dist),
            "demo_spectrum": [[float(x.real), float(x.imag)] for x in lam_b],
            "demo_V": float(pot.v_mat(Mrest)),
            "min_V_over_null_family": lowest,
            "reading": "Frobenius distance is not Lorentz invariant and the orbit is not "
            "compact: states with a complex spectrum sit arbitrarily close to highly boosted "
            "vacua, at large V. V itself is invariant, so the invariant statement is the "
            "collision barrier of A5iid. No V < 0 state was found in the zero-distance family.",
        },
        "min V over the null family >= -1e-10 (no V < 0 at zero orbit distance)",
        "V < -1e-10 somewhere in the null family: V < 0 arbitrarily close (Frobenius) to the "
        "orbit",
        lowest > -1e-10,
    )

    # (iii) block-diagonal fields
    n = 200000
    S = RNG.normal(size=(n, 3, 3)) * RNG.choice([0.1, 1.0, 10.0], size=(n, 1, 1))
    S = S + np.swapaxes(S, 1, 2)
    Mb = np.zeros((n, 4, 4))
    Mb[:, 1:, 1:] = S
    Mb[:, 0, 0] = RNG.normal(size=n) * 10
    im = float(np.abs(np.linalg.eigvals(Mb @ ETA).imag).max())
    Malt = np.diag([-DD, -GG, 1.0, DD])
    v_alt = float(pot.v_mat(Malt))
    Halt = np.zeros((10, 10))
    va = to_vec(Malt)
    for k in range(10):
        ek = np.zeros(10)
        ek[k] = 1e-5
        Halt[:, k] = (pot.v_and_grad_vec(va + ek)[1] - pot.v_and_grad_vec(va - ek)[1]) / 2e-5
    Halt = 0.5 * (Halt + Halt.T)
    blk = [k for k in range(10) if LABELS[k][0] != "0" or LABELS[k] == "00"]
    add(
        "A5iii",
        "a block-diagonal field (M_0i = 0) always has a real spectrum, so the instrument cannot "
        "reach the V < 0 states",
        "200000 random block-diagonal symmetric M at three scales: max |Im lambda|; plus the "
        "block-diagonal twin diag(-delta, -g, 1, delta) (timelike eigenvalue delta): its V, the "
        "lowest Hessian eigenvalue inside the block-diagonal subspace and on all ten entries",
        {
            "max_abs_imag": im,
            "V_of_block_diagonal_twin": v_alt,
            "twin_lowest_hessian_eig_block_diagonal_subspace": float(
                np.linalg.eigvalsh(Halt[np.ix_(blk, blk)])[0]
            ),
            "twin_lowest_hessian_eig_all_entries": float(np.linalg.eigvalsh(Halt)[0]),
            "closed_form_-(g+delta)(1-delta)c": -(GG + DD) * (1 - DD) * c,
            "reading": "true as algebra (N = diag(-M_00, S), S symmetric). Not said in the "
            "record: the block-diagonal sector CONTAINS a V = 0 twin of the vacuum for every c, "
            "which is a saddle of the full problem (unstable along M_0i at -5.81 c) and stable "
            "inside the block-diagonal subspace; the vacuum is a global minimum there, not the "
            "unique one.",
        },
        "max |Im| = 0",
        "a block-diagonal M with a complex eigenvalue",
        im < 1e-12,
    )


# ---------------------------------------------------------------- A6
def halo_closed(r, beta):
    nu = np.sqrt(5.0) / 4.0
    return np.sqrt(r) * kv(nu, 0.5 * np.sqrt(beta) * r * r)


def rk4_inward(beta, r_out, r_in, n=40000):
    """own fixed-step RK4 on u = ln eps: state (eps, eps') rescaled each step to avoid overflow"""
    h = (r_in - r_out) / n
    y = np.array([1.0, -np.sqrt(beta) * r_out - 0.5 / r_out])  # WKB tail incl. the prefactor
    r, logs = r_out, 0.0
    rs, ys, ls = [r], [y.copy()], [0.0]

    def f(r_, y_):
        return np.array([y_[1], (1.0 / r_**2 + beta * r_**2) * y_[0]])

    for _ in range(n):
        k1 = f(r, y)
        k2 = f(r + h / 2, y + h / 2 * k1)
        k3 = f(r + h / 2, y + h / 2 * k2)
        k4 = f(r + h, y + h * k3)
        y = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        r += h
        s = abs(y[0])
        y, logs = y / s, logs + np.log(s)
        rs.append(r)
        ys.append(y.copy())
        ls.append(logs)
    return np.array(rs), np.array(ys), np.array(ls)


def check_a6():
    K = 8 * (1 - DD) ** 2
    p = (np.sqrt(5.0) - 1) / 2
    rec = {3e-5: 18.514450722536125, 1e-4: 13.982649132456624, 3e-4: 10.977173858692936}
    rec.update({1e-3: 8.632181609080455, 3e-3: 7.155807790389519})
    rows, radii, worst, worst_rk = [], [], 0.0, 0.0
    for c in CS:
        beta = (GG + DD) * (1 - DD) * c / (2 * K)

        def fun(r, beta=beta):
            return halo_closed(r, beta) * r**p - 0.5 * halo_closed(4.5, beta) * 4.5**p

        rr = brentq(fun, 4.5001, 80.0, xtol=1e-12)
        # RK4 route
        rs, ys, ls = rk4_inward(beta, 7.5 * beta**-0.25, 4.0)
        logf = ls + np.log(np.abs(ys[:, 0])) + p * np.log(rs)
        ref = np.interp(4.5, rs[::-1], logf[::-1])
        ri, li = rs[::-1], logf[::-1]
        k = np.where((li < ref + np.log(0.5)) & (ri > 4.5))[0][0]
        r_rk = float(np.interp(-(ref + np.log(0.5)), -li[k - 1 : k + 1], ri[k - 1 : k + 1]))
        rows.append(
            {"c": c, "beta": beta, "reader_closed_form": rr, "reader_rk4": r_rk, "record": rec[c]}
        )
        radii.append(rr)
        worst = max(worst, abs(rr - rec[c]))
        worst_rk = max(worst_rk, abs(r_rk - rr))
    s5 = float(np.polyfit(np.log(CS), np.log(radii), 1)[0])
    s3 = float(np.polyfit(np.log(CS[:3]), np.log(radii[:3]), 1)[0])
    add(
        "A6a",
        "the reader 'eps r^0.618 falls to half of its r = 4.5 value' returns 18.5, 14.0, 11.0, "
        "8.6, 7.2; log slope -0.21 (five) and -0.23 (lowest three)",
        "closed form eps = sqrt(r) K_nu(sqrt(beta) r^2 / 2), nu = sqrt(5)/4, root by brentq; "
        "cross-checked by an own fixed-step RK4 inward integration in log amplitude started on "
        "the WKB tail at 7.5 beta^(-1/4)",
        {
            "rows": rows,
            "slope_all_five": s5,
            "slope_lowest_three": s3,
            "max_abs_diff_to_record": worst,
            "max_abs_diff_rk4_vs_closed": worst_rk,
        },
        "each radius within 0.01 of the record (its grid step is 0.003); slopes within 0.003 of "
        "-0.2071 and -0.2271; RK4 within 0.01 of the closed form",
        "a radius off by more than 0.01, or a slope off by more than 0.003",
        worst < 0.01 and abs(s5 + 0.2071) < 0.003 and abs(s3 + 0.2271) < 0.003 and worst_rk < 0.01,
    )
    # intrinsic radius scaling and the power law
    betas = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2]
    rad = []
    for b in betas:

        def slope(r, b=b):
            e = 1e-6 * r
            return (
                r * (np.log(halo_closed(r + e, b)) - np.log(halo_closed(r - e, b))) / (2 * e) + 2
            )

        rad.append(brentq(slope, 0.05 * b**-0.25, 5 * b**-0.25))
    expo = float(np.polyfit(np.log(betas), np.log(rad), 1)[0])
    # own RK4 check of the same radius at one beta
    rs, ys, ls = rk4_inward(1e-4, 75.0, 1.0)
    ri, si = rs[::-1], (rs * ys[:, 1] / ys[:, 0] + 2)[::-1]
    kk = np.where(si < 0)[0][0]
    r_rk = float(np.interp(0.0, -si[kk - 1 : kk + 1], ri[kk - 1 : kk + 1]))
    # small-r power law: the local log slope of the closed form, against r^((1 - sqrt 5)/2)
    small = [float(np.log(halo_closed(2e-3, 1e-4) / halo_closed(1e-3, 1e-4)) / np.log(2.0))]
    add(
        "A6b",
        "every intrinsic radius of the decaying solution scales exactly as beta^(-1/4); without "
        "the potential term the decaying power law is r^((1 - sqrt 5)/2) = r^(-0.618)",
        "radius where the log slope of the closed form reaches -2, five decades of beta, "
        "fitted exponent; the same radius from the RK4 route at beta 1e-4; the local log slope "
        "of the closed form at r = 1e-3 to 2e-3; sympy check that r^((1 - sqrt 5)/2) solves "
        "eps'' = eps / r^2",
        {
            "exponent": expo,
            "radius_closed_beta1e-4": rad[2],
            "radius_rk4_beta1e-4": r_rk,
            "small_r_log_slope": small[0],
            "(1-sqrt5)/2": (1 - np.sqrt(5)) / 2,
            "indicial_residual": str(
                sp.simplify(((1 - sp.sqrt(5)) / 2) * ((1 - sp.sqrt(5)) / 2 - 1) - 1)
            ),
        },
        "exponent -0.25 within 1e-6; RK4 radius within 1e-3 relative; small-r slope -0.618 "
        "within 1e-4; indicial residual 0",
        "an exponent off -1/4 or a different small-r power",
        abs(expo + 0.25) < 1e-6
        and abs(r_rk / rad[2] - 1) < 1e-3
        and abs(small[0] - (1 - np.sqrt(5)) / 2) < 1e-4,
    )


def check_a6c():
    r, K, m2 = sp.symbols("r K m2", positive=True)
    e = sp.Function("e")
    # energy per unit solid angle: r^2 [ (K / r^2)(2 e'^2 + 2 e^2 / r^2) + m2 e^2 ]
    dens = r**2 * ((K / r**2) * (2 * e(r).diff(r) ** 2 + 2 * e(r) ** 2 / r**2) + m2 * e(r) ** 2)
    el = sp.diff(dens, e(r)) - sp.diff(sp.diff(dens, e(r).diff(r)), r)
    rhs = sp.solve(sp.Eq(el, 0), e(r).diff(r, 2))[0]
    want = e(r) / r**2 + (m2 / (2 * K)) * r**2 * e(r)
    diff = sp.simplify(rhs - want)
    add(
        "A6c",
        "beta = m2 / (2 K) with m2 = (g + delta)(1 - delta) c: the halo equation follows from "
        "e_2 = (K / r^2)(2 eps_r^2 + 2 eps^2 / r^2) plus the potential m2 eps^2 (the -c L "
        "quadratic term of A2), radial measure r^2 dr",
        "sympy Euler-Lagrange equation of r^2 (e_2 + m2 eps^2), solved for eps''",
        {
            "eps''": str(sp.simplify(rhs)),
            "difference_from_claimed": str(diff),
            "note": "K = 8 (1 - delta)^2 and the angular factor of the potential term are taken "
            "from the record (not audited here); a potential (m2 / 2) eps^2 would halve beta",
        },
        "difference 0",
        "a factor 2 in beta (potential normalized as m2 eps^2 / 2) or a different centrifugal "
        "coefficient",
        diff == 0,
    )


# ---------------------------------------------------------------- A7
def check_a7():
    pot = Pot(GG, DD, WW, 0.0)
    a_iso = (1 + 2 * DD) / 3
    v1 = pot.v_mat(np.diag([GG, a_iso, a_iso, a_iso]))
    from scipy.optimize import minimize_scalar

    res = minimize_scalar(
        lambda m: pot.v_mat(np.diag([m, a_iso, a_iso, a_iso])),
        bounds=(7.0, 9.0),
        method="bounded",
        options={"xatol": 1e-12},
    )
    v3 = pot.v_mat(np.diag([GG, 0.0, 0.0, 0.0]))
    add(
        "A7",
        "V4 at the isotropic block a I, a = (1 + 2 delta)/3: 0.019253 (M_00 = g), 0.009955 "
        "(M_00 re-solved), 0.11034 (zero spatial block)",
        "V4 from the traces of the explicit 4 x 4 matrices; M_00 by bounded scalar minimization",
        {
            "M00_g": float(v1),
            "M00_solved": float(res.fun),
            "M00_at_min": float(res.x),
            "zero_block": float(v3),
            "record_M00_solved": 0.009954763483326615,
            "note": "the record's 0.00995476 comes from a 1e-5 grid in M_00 on a curvature of "
            "1.5e5, good to about 2e-6 only; the minimum is 0.00995457 at M_00 = 8.00035. The "
            "quoted 0.009955 holds.",
        },
        "0.019253, 0.009955, 0.11034 at the quoted precision (5, 4 and 5 significant figures)",
        "any of the three differing at the quoted precision",
        f"{v1:.5g}" == "0.019253" and f"{res.fun:.4g}" == "0.009955" and f"{v3:.5g}" == "0.11034",
    )


# ---------------------------------------------------------------- B
def one_sided(f, ax, fwd, h):
    d = np.diff(f, axis=ax) / h
    pad = [(0, 0)] * f.ndim
    pad[ax] = (0, 1) if fwd else (1, 0)
    return np.pad(d, pad)


def bracket(A, B):
    return A @ ETA @ B - B @ ETA @ A


def contract(F, X):
    Y = np.swapaxes(X, -1, -2) @ F @ X if X.ndim > 2 else X.T @ F @ X
    return np.sum(Y * F, axis=(-1, -2))


def metric_G(M, delta):
    x = ETA @ M
    I4 = np.eye(4)
    q = x @ (x - I4) @ (x - delta * I4) / (GG * (GG - 1) * (GG - delta))
    return ETA - 2.0 * q @ ETA


def load_field():
    Mraw = np.load(os.path.join(EXT, "M_G_polished_N32.npz"))["M"]
    seed = np.load(os.path.join(EXT, "seeds", "m5_21_2b_end_A_T2_sym_e0_n32_d0.3_pinned.npz"))[
        "M"
    ].astype(np.float64)
    n = Mraw.shape[0]
    wc = max(1, int(np.ceil(1.6 / 1.5)))
    idx = np.arange(n)
    edge = (idx < wc) | (idx >= n - wc)
    shell = edge[:, None, None] | edge[None, :, None] | edge[None, None, :]
    full = np.zeros((n, n, n, 4, 4))
    full[..., 1:, 1:] = seed
    full[..., 0, 0] = -GG
    Msym = 0.5 * (Mraw + np.swapaxes(Mraw, -1, -2))
    Msym[shell] = full[shell]
    return Msym


def tangent_a0(M, h, renv=10.0):
    n = M.shape[0]
    x = (np.arange(n) - (n - 1) / 2.0) * h
    r2 = x[:, None, None] ** 2 + x[None, :, None] ** 2 + x[None, None, :] ** 2
    env = np.exp(-((np.sqrt(r2) / renv) ** 4))
    Wb = np.zeros((4, 4))
    Wb[0, 1] = Wb[1, 0] = 1.0
    a = env[..., None, None] * (Wb @ M + M @ Wb.T)
    return a / np.sqrt(np.sum(a * a))


def well_numbers(M, h, X):
    a0 = tangent_a0(M, h)
    i1s = np.zeros(M.shape[:3])
    k1 = np.zeros(M.shape[:3])
    for fwd in (True, False):
        A = [one_sided(M, ax, fwd, h) for ax in range(3)]
        for i in range(3):
            k1 += 2.0 * contract(bracket(a0, A[i]), X)
            for j in range(i + 1, 3):
                i1s += 2.0 * contract(bracket(A[i], A[j]), X)
    C1 = h**3 * float(np.sum(i1s * k1))
    C2 = h**3 * float(np.sum(k1 * k1))
    return C1, C2, i1s, k1


def check_b():
    h = 1.5
    M = load_field()
    Gm = metric_G(M, DD)
    dev = float(np.abs(Gm - np.eye(4)).max())
    m0i = float(np.abs(M[..., 0, 1:]).max())
    gmin = float(np.linalg.eigvalsh(0.5 * (Gm + np.swapaxes(Gm, -1, -2))).min())
    asym = float(np.abs(Gm - np.swapaxes(Gm, -1, -2)).max())
    add(
        "B1",
        "on the committed N = 32 field with its pinned shell, G = eta - 2 q eta differs from the "
        "identity by at most 3.9e-4 and M_0i = 0 exactly: the report's contraction is the "
        "positive one",
        "own rebuild of the pinned shell (boolean edge mask, width ceil(1.6/1.5) = 2 cells) and "
        "of G; max |G - 1|, max |M_0i|, the lowest eigenvalue of G over the grid",
        {
            "max_abs_G_minus_identity": dev,
            "max_abs_M0i": m0i,
            "min_eig_G": gmin,
            "max_abs_G_asymmetry": asym,
        },
        "max |G - 1| = 3.9220e-4 (record) within 1e-9; M_0i = 0; min eig G > 0",
        "a cell with |G - 1| above 3.93e-4, a non-zero M_0i, or a non-positive G eigenvalue",
        abs(dev - 0.00039219775454779615) < 1e-9 and m0i == 0.0 and gmin > 0,
    )
    rec = {
        "G": (9.075435624180153e-06, 8.564773350793515e-05),
        "identity": (9.072930577981556e-06, 8.561973197903893e-05),
        "eta": (-9.072930577981556e-06, 8.561973197903893e-05),
    }
    vals, ok, signs = {}, True, {}
    for name, X in (("G", Gm), ("identity", np.eye(4)), ("eta", ETA)):
        C1, C2, i1s, k1 = well_numbers(M, h, X)
        vals[name] = {"C1": C1, "C2": C2, "omega_E": float(np.sqrt(C1 / C2)) if C1 > 0 else None}
        signs[name] = {
            "min_i1s": float(i1s.min()),
            "min_k1": float(k1.min()),
            "max_k1": float(k1.max()),
        }
        ok &= abs(C1 / rec[name][0] - 1) < 1e-9 and abs(C2 / rec[name][1] - 1) < 1e-9
    ok &= abs(vals["G"]["omega_E"] - 0.32551859539640887) < 1e-10
    add(
        "B2",
        "C1 = +9.0754e-6 (G), +9.0729e-6 (identity), -9.0729e-6 (eta); C2 = 8.5648e-5 (G); "
        "omega_E = 0.32551859539640887 (G)",
        "own stencils (np.diff padded one-sided differences, fwd/bwd averaged), own bracket, the "
        "contraction as sum((X^T F X) * F), a0 normalized over the grid",
        {"values": vals, "pointwise": signs},
        "C1 and C2 within 1e-9 relative of the record under all three contractions; omega_E "
        "within 1e-10",
        "any C1 or C2 off by more than 1e-9 relative",
        ok,
    )
    # B3a: algebra
    al = sp.symbols("al0:3")
    Ss = [
        sp.Matrix(3, 3, lambda i, j, k=k: sp.Symbol(f"s{k}_{min(i, j)}{max(i, j)}"))
        for k in range(3)
    ]
    As = [sp.diag(al[k], Ss[k]) for k in range(3)]
    vv = sp.symbols("v1:4")
    a0s = sp.zeros(4, 4)
    for i in range(3):
        a0s[0, i + 1] = a0s[i + 1, 0] = vv[i]
    et = sp.diag(-1, 1, 1, 1)
    leak_rot, leak_boost, antis = 0, 0, 0
    for i in range(3):
        F = (a0s * et * As[i] - As[i] * et * a0s).expand()
        leak_boost += sum(1 for p in range(1, 4) for q in range(1, 4) if F[p, q] != 0)
        leak_boost += 1 if F[0, 0] != 0 else 0
        antis += 0 if (F + F.T).expand() == sp.zeros(4, 4) else 1
        for j in range(i + 1, 3):
            Fr = (As[i] * et * As[j] - As[j] * et * As[i]).expand()
            leak_rot += sum(1 for q in range(4) if Fr[0, q] != 0 or Fr[q, 0] != 0)
            antis += 0 if (Fr + Fr.T).expand() == sp.zeros(4, 4) else 1
    # a0 of a block-diagonal M has only 0i entries
    mm = sp.symbols("mm")
    Mb = sp.diag(mm, Ss[0])
    Wb = sp.zeros(4, 4)
    Wb[0, 1] = Wb[1, 0] = 1
    a0b = Wb * Mb + Mb * Wb.T
    a0_leak = sum(1 for p in range(1, 4) for q in range(1, 4) if a0b[p, q] != 0) + (a0b[0, 0] != 0)
    # sign under eta for a pure boost-block F: <F, F>_eta = -2 sum F_0i^2
    f = sp.symbols("f1:4")
    Fb = sp.zeros(4, 4)
    for i in range(3):
        Fb[0, i + 1], Fb[i + 1, 0] = f[i], -f[i]
    Yb = et.T * Fb * et
    eta_val = sp.expand(sum(Yb[p, q] * Fb[p, q] for p in range(4) for q in range(4)))
    add(
        "B3a",
        "on a block-diagonal field i1s lives entirely in the rotation block and k1 entirely in "
        "the boost block; under eta k1 = -2 |F_0i|^2 <= 0",
        "sympy with generic block-diagonal symmetric A_i = diag(alpha_i, S_i) and a0 with only "
        "0i entries: count the non-zero entries outside the claimed block; a0 = W M + M W^T on "
        "a block-diagonal M; <F, F>_eta on a pure boost-block F",
        {
            "rotation_bracket_entries_in_boost_block": leak_rot,
            "boost_bracket_entries_outside_boost_block": leak_boost,
            "non_antisymmetric_brackets": antis,
            "a0_entries_outside_0i": int(a0_leak),
            "<F,F>_eta_on_boost_block": str(eta_val),
        },
        "all leak counts 0; <F,F>_eta = -2 (f1^2 + f2^2 + f3^2)",
        "any non-zero entry outside the claimed block",
        leak_rot == 0
        and leak_boost == 0
        and antis == 0
        and a0_leak == 0
        and sp.expand(eta_val + 2 * sum(x**2 for x in f)) == 0,
    )
    # B3b: positivity needs no block structure at all
    worst = np.inf
    for _ in range(20000):
        Fm = RNG.normal(size=(4, 4))
        Fm = Fm - Fm.T
        Rm = RNG.normal(size=(4, 4))
        Xp = Rm.T @ Rm + 1e-3 * np.eye(4)
        val = float(contract(Fm, Xp))
        chol = np.linalg.cholesky(Xp).T  # Xp = chol^T chol
        alt = float(np.sum((chol @ Fm @ chol.T) ** 2))
        worst = min(worst, val)
        if abs(val - alt) > 1e-8 * max(1.0, abs(alt)):
            worst = -np.inf
    add(
        "B3b",
        "WORDING: 'on a block-diagonal field ... so under any positive contraction i1s >= 0 and "
        "k1 >= 0': the block structure is presented as the reason for positivity",
        "identity <F, F>_X = || R F R^T ||_F^2 for X = R^T R: 20000 random antisymmetric F (all "
        "six entries, boost and rotation mixed) against random positive definite X, both routes",
        {"min_<F,F>_X_over_samples": worst},
        "min >= 0 and the two routes agree: positivity holds for ANY F, block structure is not "
        "needed (it is needed only for the sign under eta)",
        "a negative <F, F>_X for a positive definite X, or the two routes disagreeing",
        worst >= 0,
    )
    # B3c: perturb M_0i
    rows = {}
    ok = True
    n = M.shape[0]
    xg = (np.arange(n) - (n - 1) / 2.0) * h
    Xg, Yg, Zg = np.meshgrid(xg, xg, xg, indexing="ij")
    envs = np.exp(-((Xg**2 + Yg**2 + Zg**2) / 15.0**2))
    smooth = np.zeros(M.shape[:3] + (3,))
    for i in range(3):
        for _ in range(4):
            kx, ky, kz = RNG.integers(-2, 3, size=3) * 2 * np.pi / (n * h)
            smooth[..., i] += RNG.normal() * np.cos(
                kx * Xg + ky * Yg + kz * Zg + RNG.uniform(0, 6.28)
            )
    smooth *= envs[..., None] / np.abs(smooth).max()
    for kind, amp in [("noise", a_) for a_ in (1e-3, 1e-2, 1e-1, 1.0)] + [
        ("smooth", a_) for a_ in (1e-2, 1e-1, 1.0)
    ]:
        Mp = M.copy()
        pert = (RNG.normal(size=M.shape[:3] + (3,)) if kind == "noise" else smooth) * amp
        Mp[..., 0, 1:] += pert
        Mp[..., 1:, 0] += pert
        Gp = metric_G(Mp, DD)
        gmin_p = float(np.linalg.eigvalsh(0.5 * (Gp + np.swapaxes(Gp, -1, -2))).min())
        row = {"min_eig_G": gmin_p}
        for name, X in (("G", Gp), ("identity", np.eye(4)), ("eta", ETA)):
            C1, C2, i1s, k1 = well_numbers(Mp, h, X)
            row[name] = {
                "C1": C1,
                "C2": C2,
                "frac_i1s_negative": float(np.mean(i1s < 0)),
                "frac_k1_negative": float(np.mean(k1 < 0)),
            }
        rows[f"{kind} amp={amp:g}"] = row
        ok &= row["identity"]["C1"] > 0 and row["identity"]["frac_k1_negative"] == 0.0
        if gmin_p > 0:
            ok &= row["G"]["C1"] > 0
    add(
        "B3c",
        "C1 >= 0 under a positive contraction for ANY field (tested off the block-diagonal "
        "sector: a symmetric M_0i added to the committed field, white noise cell by cell at "
        "amplitude 1e-3 to 1 and a smooth low-k packet at 1e-2 to 1)",
        "rebuild G, a0, i1s, k1 on the perturbed field; C1 and the pointwise sign fractions under "
        "G, identity and eta; the lowest eigenvalue of G (is the report's G still positive?)",
        rows,
        "identity: C1 > 0 and no negative cell at every amplitude; G: C1 > 0 wherever G stays "
        "positive definite; eta: signs become mixed once M_0i != 0",
        "C1 <= 0 or a negative k1 cell under the identity contraction, or C1 <= 0 under a "
        "positive definite G",
        ok,
    )
    # B3d: is the well informative? random fields with no soliton at all
    rnd = {}
    ok_r = True
    for k in range(3):
        n = 16
        Mr = np.zeros((n, n, n, 4, 4))
        S = RNG.normal(size=(n, n, n, 3, 3)) * 0.2
        Mr[..., 1:, 1:] = S + np.swapaxes(S, -1, -2) + np.diag([1.0, DD, DD])
        Mr[..., 0, 0] = -GG + RNG.normal(size=(n, n, n)) * 0.2
        C1, C2, _, _ = well_numbers(Mr, h, np.eye(4))
        rnd[f"noise_field_{k}"] = {"C1": C1, "C2": C2, "omega_E": float(np.sqrt(C1 / C2))}
        ok_r &= C1 > 0
    Mu = np.zeros((8, 8, 8, 4, 4))
    Mu[...] = np.diag([-GG, 1.0, DD, DD])
    C1u, C2u, _, _ = well_numbers(Mu, h, np.eye(4))
    add(
        "B3d",
        "the well (C1 > 0, a real omega_E) exists for any field under a positive contraction",
        "C1 and omega_E on three white-noise block-diagonal fields with no soliton, and on the "
        "uniform vacuum (where every derivative vanishes)",
        {
            "noise": rnd,
            "uniform_vacuum": {"C1": C1u, "C2": C2u},
            "reading": "C1 = h^3 sum i1s k1 is a sum of products of non-negative densities, so "
            "C1 > 0 for every non-uniform field, noise included: the EXISTENCE of the "
            "frozen-profile "
            "well carries no information about the field or its dynamics. It is fixed by the form "
            "(i1s - omega^2 k1)^2 and the sign of the contraction on the boost block. Only the "
            "magnitude of omega_E could inform, and that is normalization dependent (a0 norm, "
            "envelope, h^(-3/2)).",
        },
        "C1 > 0 on noise; C1 = 0 on the uniform vacuum",
        "a non-uniform field with C1 <= 0 under the identity contraction",
        ok_r and C1u == 0.0,
    )


def main():
    check_a1()
    check_a2()
    check_a3()
    check_a4()
    pot, Ms = check_a5()
    check_a5_paths(pot, Ms)
    check_a6()
    check_a6c()
    check_a7()
    check_b()
    n_pass = sum(1 for c in CHECKS if c["verdict"] == "PASS")
    out = {
        "audit": "M5.32 R23-0 form record and R23-3 well record, independent refutation attempt",
        "checks": CHECKS,
        "runtime_s": time.time() - T0,
        "counts": {"total": len(CHECKS), "PASS": n_pass, "FAIL": len(CHECKS) - n_pass},
    }
    with open(OUT_JSON, "w") as fjs:
        json.dump(out, fjs, indent=1, default=str)
    print("counts:", out["counts"], f"runtime {out['runtime_s']:.0f} s")


if __name__ == "__main__":
    main()
