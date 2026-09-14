"""M5.32 R21 plan check on the author's 2026-09-14 14:45 UTC reply (discussioncomment-18435636).

Two form-level claims, each one function, the audited versions in R21-0 (e) and (f):

(1) The light-trajectory slope nu' = -V_{nu Delta} / V_{nu nu} on the spectrum (nu + Delta, nu, nu, nu),
    nu(Delta) the valley (argmin over nu at fixed Delta). The author: -1/4 at Delta = 0 for every
    separable V = sum_i f(lambda_i), evaded only by rank weighting a f(lambda_max) + b sum_rest f
    (nu' = -a / (a + 3 b)). Here: V4 = W1 sum_p (tr N^p - C_p)^2 (symmetric in the eigenvalues of N,
    not separable) at Delta = 0 and at the hierarchy Delta = -9, with the separable and the
    rank-weighted controls, and the symmetric-point identity V_nD / V_nn = (A + 3 B) / (4 A + 12 B).

(2) The boost dressing M -> Q M Q^T, Q = exp(s K), K = E_01 + E_10 (Q^T eta Q = eta): the spectrum of
    N = M eta is conjugated (silent to all orders for every eigenvalue function of N, and det M with
    det Q = 1); the Euclidean norm tr(M^2) has first derivative 4 tr(M^2 K) (zero on block-diagonal
    fields) and second derivative 8 tr(K^2 M^2) + 8 tr(K M K M).

Run from the repo root: python3 openwave/xperiments/m5_liquid_crystal/research/scripts/m5_32_r21_plan_check.py
Writes data/m5_32_r21_plan_check.json next to the other M5.32 records.
"""

import json
import os

import numpy as np
import sympy as sp
from scipy.linalg import expm

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "data", "m5_32_r21_plan_check.json")

W1 = sp.Rational(724023879, 10**12)  # B3.W1 to the digits used in the stack
SPEC_VAC = [-8, 1, sp.Rational(3, 10), 0]  # the N-spectrum of the certified branch (s = -1)


def valley_slopes(V, nu, D, at_D):
    """(nu0, nu') on every valley branch (V_nn > 0) of V(nu, Delta) at Delta = at_D."""
    Vn, Vnn, VnD = sp.diff(V, nu), sp.diff(V, nu, 2), sp.diff(V, nu, D)
    out = []
    for r in sp.Poly(Vn.subs(D, at_D), nu).nroots():
        if abs(sp.im(r)) > 1e-12:
            continue
        r = sp.re(r)
        if Vnn.subs({nu: r, D: at_D}) > 0:
            out.append((float(r), float(-VnD.subs({nu: r, D: at_D}) / Vnn.subs({nu: r, D: at_D}))))
    return out


def part1():
    nu, D = sp.symbols("nu Delta", real=True)
    lam = [nu + D, nu, nu, nu]
    C = [sum(sp.Integer(1) * l**p for l in SPEC_VAC) for p in range(1, 5)]
    V4 = W1 * sum((sum(l**p for l in lam) - C[p - 1]) ** 2 for p in range(1, 5))
    f = lambda l: (l - 2) ** 2 + l**4 / 7  # a non-quadratic separable control
    Vsep = sum(f(l) for l in lam)
    a, b = sp.Integer(3), sp.Integer(1)  # rank weighting: lambda_max = nu + D for D > 0
    Vrank = a * f(nu + D) + b * sum(f(l) for l in lam[1:])
    # the symmetric-point identity for any symmetric V: V_nD / V_nn = (A + 3 B) / (4 A + 12 B) = 1 / 4
    l1, l2, l3, l4 = sp.symbols("l1:5", real=True)
    S = W1 * sum((sum(l**p for l in (l1, l2, l3, l4)) - C[p - 1]) ** 2 for p in range(1, 5))
    x = sp.symbols("x", real=True)
    sub = {l1: x, l2: x, l3: x, l4: x}
    A = sp.simplify(sp.diff(S, l1, 2).subs(sub))
    B = sp.simplify(sp.diff(S, l1, l2).subs(sub))
    ratio = sp.simplify((A + 3 * B) / (4 * A + 12 * B))
    res = {
        "separable_control_D0": valley_slopes(Vsep, nu, D, 0),
        "rank_weighted_control_D0_a3_b1": valley_slopes(Vrank, nu, D, sp.Rational(1, 10**6)),
        "rank_weighted_prediction_-a/(a+3b)": float(-a / (a + 3 * b)),
        "V4_D0": valley_slopes(V4, nu, D, 0),
        "V4_D-9_hierarchy": valley_slopes(V4, nu, D, -9),
        "V4_D8": valley_slopes(V4, nu, D, 8),
        "symmetric_point_identity_(A+3B)/(4A+12B)": str(ratio),
    }
    for k, v in res.items():
        print(f"(1) {k}: {v}")
    return res


def part2():
    eta = np.diag([-1.0, 1, 1, 1])
    K = np.zeros((4, 4))
    K[0, 1] = K[1, 0] = 1.0
    Q = lambda s: expm(s * K)
    assert np.allclose(Q(0.3).T @ eta @ Q(0.3), eta)
    rng = np.random.default_rng(0)
    Mv = np.diag([8.0, 1, 0.3, 0])
    Bm = rng.standard_normal((3, 3))
    Mb = np.zeros((4, 4))
    Mb[0, 0], Mb[1:, 1:] = 8.0, Bm + Bm.T  # block diagonal, M_0i = 0 (R20's static fields)
    Mg = rng.standard_normal((4, 4))
    Mg = Mg + Mg.T
    res = {}
    for name, M in (("vacuum", Mv), ("block_diagonal", Mb), ("general", Mg)):
        Ms = lambda s: Q(s) @ M @ Q(s).T
        spec = lambda X: np.sort(np.linalg.eigvals(X @ eta).real)
        drift = max(np.max(np.abs(spec(Ms(s)) - spec(M))) for s in (0.02, 0.05, 0.3))
        n = lambda s: np.trace(Ms(s) @ Ms(s))
        h = 1e-3
        d1, d2 = (n(h) - n(-h)) / (2 * h), (n(h) - 2 * n(0) + n(-h)) / h**2
        closed1 = 4 * np.trace(M @ M @ K)
        closed2 = 8 * np.trace(K @ K @ M @ M) + 8 * np.trace(K @ M @ K @ M)
        ddet = np.linalg.det(Ms(0.05)) - np.linalg.det(M)
        res[name] = {
            "N_spectrum_drift": float(drift),
            "d1_trM2": float(d1),
            "d1_closed": float(closed1),
            "d2_trM2": float(d2),
            "d2_closed": float(closed2),
            "detM_drift_s005": float(ddet),
        }
        print(
            f"(2) {name:15s} spectrum drift {drift:.1e}; d tr(M^2)/ds {d1:+.4f} (closed {closed1:+.4f}); "
            f"d2/ds2 {d2:.3f} (closed {closed2:.3f}); det M drift {ddet:+.1e}"
        )
    return res


if __name__ == "__main__":
    out = {"part1_valley_slope": part1(), "part2_norm_along_boost": part2()}
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    print("wrote", os.path.relpath(OUT))
