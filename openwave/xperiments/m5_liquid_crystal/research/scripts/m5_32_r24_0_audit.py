"""M5.32 R24-0 adversarial audit: the form-level claims C1 to C6.

This is the adversarial audit of the R24-0 form-level claims. It was written by a fresh agent
that did not read the audited script (m5_32_r24_0_form.py): the agent received the claims as
statements and built its own code for each one. Nothing is imported from the audited script.
Nine independent audit scripts were consolidated here into one file, one function per check
group; the methods are the auditor's own and were kept as written.

CLAIMS UNDER ATTACK
-------------------
C1  V4 on the diagonal split (-g, 1, delta + eps, delta - eps) equals
    4 eps^4 (1 + 9 d^2 + 36 d^4 + 12 d^2 eps^2 + eps^4): no eps^2 term.
C2  The decaying solution of eps'' = (A / r^2 + beta r^2) eps is sqrt(r) K_nu(sqrt(beta) r^2 / 2)
    with nu = sqrt(1 + 4 A) / 4, and the form sqrt(rho) K_nu is off by 31 to 39 percent.
C3  The Riesz cluster projector gives a doublet block B that vanishes on every Lorentz image of
    the vacuum, has tr B^2 = 2 eps^2 on every image of the split, and is smooth through the
    coalescence; the Euclidean construction is not invariant.
C4  On a boost ripple exp(theta(x) K_a) of a split background tr(dB dB) = -2 |B a|^2 theta'^2,
    while the Frobenius and the h-metric contractions equal +2 |B a|^2 theta'^2.
C5  The same ripple on the stored relaxed field lowers the production energy, through the
    curvature term alone (dE = -3.283, -10.947, -68.457 at the three claimed (A, k)).
C6  Second variation of the curvature density around the hedgehog: (a) no linear term on the
    doublet jets, (b) the listed jet table, (c) lowest angular eigenvalue 48 (A = 3, ladder
    48, 96, 160), (iii) the coupling to the director does not change it.

VERDICTS
--------
Each check records PASS or FAIL against the one-line claim stored beside it, where FAIL means
the claim was refuted or not reproduced. kind = "claim" attacks a statement of C1 to C6;
kind = "scope" attacks an extension of a claim beyond the conditions it was stated for (a
FAIL there marks a limit of validity, not an error in the claim as stated). could_fail = false
marks a check that is an identity once another check holds. The auditor's scripts printed
numbers and carried no thresholds; the thresholds here were attached at consolidation, after
the auditor's logs existed, and sit orders of magnitude away from the observed values.

INDEPENDENCE
------------
Own code throughout: sympy and mpmath for C1 and C2, a Riesz projector from the eigenvectors
with eta-norms for C3 and C4, a Cartesian sympy jet calculus for the pointwise C6, a Cartesian
tensor Galerkin method on the sphere for the angular problem (no theta grid, no spin-weighted
harmonic formula), a brute-force finite-difference second variation on the exact nonlinear
family, and an own lattice energy for C5. The C5 part also imports the production module
m5_32_r23_1_cscan.py (same folder) to evaluate the production energy on the stored field; if
the stored field file is absent that part reports SKIPPED.

Modes:
    run     C1, C2, C3, C4, the four C6 groups and the extras (no production code)
    stack   C5 and the vacuum artifact (the parts that go with the production stack)
Run: python3 m5_32_r24_0_audit.py run | stack
Runtime: run about 1 min, stack about 30 s (2 threads, no pools).
Writes: ../data/m5_32_r24_0_audit.json (the two modes merge into one record).
"""

import os

for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_k] = "2"

import importlib.util  # noqa: E402
import itertools  # noqa: E402
import json  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from types import SimpleNamespace  # noqa: E402

import mpmath as mp  # noqa: E402
import numpy as np  # noqa: E402
import sympy as sp  # noqa: E402
from scipy.integrate import solve_ivp  # noqa: E402
from scipy.linalg import eigh, expm  # noqa: E402
from scipy.special import kv  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT_JSON = os.path.join(DATA, "m5_32_r24_0_audit.json")
PRODUCTION = "m5_32_r23_1_cscan.py"

G8, DELTA = 8.0, 0.3
ETA = np.diag([-1.0, 1, 1, 1])
PAIRS = [(0, 1), (0, 2), (1, 2)]
# the claimed C5 record: E0 and dE at (A, k h / pi)
CLAIM_E0 = 5.416820972916
CLAIM_DE = {
    (0.02, 1 / 4): -3.2829791043,
    (0.02, 1 / 2): -10.9469360784,
    (0.05, 1 / 2): -68.4567721656,
}
T0 = time.time()
CHECKS = {}
_CACHE = {}


# ================= the record =================
def _clean(x):
    if isinstance(x, dict):
        return {str(k): _clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_clean(v) for v in x]
    if isinstance(x, np.ndarray):
        return _clean(x.tolist())
    if isinstance(x, (bool, np.bool_)):
        return bool(x)
    if isinstance(x, (int, np.integer)):
        return int(x)
    if isinstance(x, (float, np.floating)):
        return float(x)
    if x is None or isinstance(x, str):
        return x
    return str(x)


def record(cid, claim, passed, numbers, could_fail=True, kind="claim", note=None):
    verdict = passed if isinstance(passed, str) else ("PASS" if passed else "FAIL")
    CHECKS[cid] = {
        "claim": claim,
        "verdict": verdict,
        "numbers": _clean(numbers),
        "could_fail": bool(could_fail),
        "kind": kind,
    }
    if note:
        CHECKS[cid]["note"] = note
    print("  [%s] %s (%s): %s" % (verdict, cid, kind, claim), flush=True)


def _counts(checks):
    v = [c["verdict"] for c in checks.values()]
    return len(v), v.count("PASS"), v.count("FAIL"), v.count("SKIPPED")


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


# ================= C1 =================
def check_c1():
    """V4 on the diagonal split, symbolic and with the matrix itself."""
    print("\n=== C1: V4 on the diagonal split ===")
    g, d, e = sp.symbols("g delta eps", real=True)
    ev_vac = [-g, 1, d, d]
    ev_spl = [-g, 1, d + e, d - e]

    def V(pmax):
        tot = 0
        for p in range(1, pmax + 1):
            tot += (sum(x**p for x in ev_spl) - sum(x**p for x in ev_vac)) ** 2
        return sp.expand(tot)

    V4 = V(4)
    claim4 = sp.expand(4 * e**4 * (1 + 9 * d**2 + 36 * d**4 + 12 * d**2 * e**2 + e**4))
    r4 = sp.simplify(V4 - claim4)
    print("C1 V4 - claim =", r4)
    print("C1 V4 factored  =", sp.factor(V4))
    V3 = V(3)
    print("C1 V3 factored  =", sp.factor(V3))
    claim3 = sp.expand(4 * e**4 * (1 + 9 * d**2))
    r3 = sp.simplify(V3 - claim3)
    print("C1 V3 - 4e^4(1+9d^2) =", r3)
    lowest = min(m[0] for m in sp.Poly(V4, e).monoms())
    # numeric spot check with the matrix itself (not eigenvalue shortcut)
    rng = np.random.default_rng(1)
    rows = []
    for _ in range(3):
        gg, dd, ee = 8.0, 0.3, rng.uniform(-0.2, 0.2)
        N = np.diag([-gg, 1, dd + ee, dd - ee])
        N0 = np.diag([-gg, 1, dd, dd])
        v = sum(
            (np.trace(np.linalg.matrix_power(N, p)) - np.trace(np.linalg.matrix_power(N0, p))) ** 2
            for p in range(1, 5)
        )
        c = 4 * ee**4 * (1 + 9 * dd**2 + 36 * dd**4 + 12 * dd**2 * ee**2 + ee**8 / ee**4)
        rel = abs(v - c) / abs(v)
        rows.append({"eps": ee, "V": v, "claim": c, "rel": rel})
        print("  numeric eps=%.4f  V=%.12e  claim=%.12e  rel=%.2e" % (ee, v, c, rel))
    record(
        "C1.v4_closed_form",
        "V4 on the split equals 4 eps^4 (1 + 9 d^2 + 36 d^4 + 12 d^2 eps^2 + eps^4), symbolically",
        r4 == 0,
        {"V4_minus_claim": str(r4), "V4_factored": str(sp.factor(V4))},
    )
    record(
        "C1.v3_closed_form",
        "the p <= 3 member equals 4 eps^4 (1 + 9 d^2), symbolically",
        r3 == 0,
        {"V3_minus_claim": str(r3), "V3_factored": str(sp.factor(V3))},
    )
    record(
        "C1.no_eps2_term",
        "the lowest power of eps in V4 on the split is 4",
        lowest == 4,
        {"lowest_power_of_eps": lowest},
    )
    record(
        "C1.numeric_matrix_traces",
        "the closed form matches V4 from matrix-power traces at three random eps (rel < 1e-6)",
        max(r["rel"] for r in rows) < 1e-6,
        {"rows": rows, "threshold": 1e-6},
    )


# ================= C2 =================
def _matches_31_39(lo, hi):
    a, b = sorted((abs(lo), abs(hi)))
    return abs(a - 0.31) < 0.03 and abs(b - 0.39) < 0.03


def check_c2():
    """Bessel form of the decaying solution, and the '31 to 39 percent' figure."""
    print("\n=== C2: Bessel form of the decaying solution ===")
    # numeric residual test of the ODE at high precision with mpmath
    mp.mp.dps = 40

    def sol(rr, Av, bv):
        n = mp.sqrt(1 + 4 * Av) / 4
        return mp.sqrt(rr) * mp.besselk(n, mp.sqrt(bv) * rr**2 / 2)

    residuals = {}
    for Av in (1, 3, 0.37):
        for bv in (7.41e-4, 0.5):
            worst = 0
            for rr in (0.7, 3.0, 9.0, 18.0):

                def f(x, Av=Av, bv=bv):
                    return sol(x, Av, bv)

                lhs = mp.diff(f, rr, 2)
                rhs = (Av / rr**2 + bv * rr**2) * f(rr)
                worst = max(worst, abs(lhs - rhs) / abs(rhs))
            residuals["A=%s beta=%s" % (Av, bv)] = float(worst)
            print("C2 ODE residual A=%s beta=%s: max rel = %s" % (Av, bv, mp.nstr(worst, 3)))

    # scaling collapse over the window r 6 to 18, and the deviation of the sqrt(rho) K_nu form
    c = 1e-3
    beta = 5.81 * c / (2 * 3.92)
    print("beta =", beta)
    collapse = {}
    readings = []
    for Av in (1.0, 3.0):
        n = np.sqrt(1 + 4 * Av) / 4

        def rho(rr):
            return np.sqrt(beta) * rr**2 / 2

        rs = np.linspace(6, 18, 13)
        eps = np.sqrt(rs) * kv(n, rho(rs))
        good = eps * rs ** (2 * n - 0.5) / (rho(rs) ** n * kv(n, rho(rs)))
        bad = eps * rs ** (2 * n - 0.5) / (np.sqrt(rho(rs)) * kv(n, rho(rs)))
        print("A=%g nu=%.6f" % (Av, n))
        print(
            "  eps r^(2nu-1/2) / (rho^nu K_nu): min %.6f max %.6f (const => spread %.2e)"
            % (good.min(), good.max(), good.max() / good.min() - 1)
        )
        print(
            "  eps r^(2nu-1/2) / (sqrt(rho) K_nu): min %.6f max %.6f ratio max/min %.4f"
            % (bad.min(), bad.max(), bad.max() / bad.min())
        )
        row = {
            "nu": n,
            "rho_nu_K_nu_spread": good.max() / good.min() - 1,
            "sqrt_rho_K_nu_max_over_min": bad.max() / bad.min(),
            "deviation_ranges": {},
        }
        readings.append(("A=%g spread" % Av, 0.0, bad.max() / bad.min() - 1))
        # 'off by' relative to a reference: the deviation range normalized at the window's
        # middle and at each end
        for ref in (0, len(rs) // 2, -1):
            dev = bad / bad[ref] - 1
            row["deviation_ranges"]["normalized at r=%.1f" % rs[ref]] = [dev.min(), dev.max()]
            readings.append(("A=%g normalized at r=%.1f" % (Av, rs[ref]), dev.min(), dev.max()))
            print(
                "    normalized at r=%.1f: deviation range [%.1f%%, %.1f%%]"
                % (rs[ref], 100 * dev.min(), 100 * dev.max())
            )
        row["rho_window"] = [rho(6), rho(18)]
        row["rho_pow_nu_minus_half_at_ends"] = [rho(6) ** (n - 0.5), rho(18) ** (n - 0.5)]
        print(
            "    rho window: %.4f .. %.4f ; rho^(nu-1/2) at ends: %.4f %.4f"
            % (rho(6), rho(18), rho(6) ** (n - 0.5), rho(18) ** (n - 0.5))
        )
        collapse["A=%g" % Av] = row

    # independent ODE integration check of the decaying branch
    Av = 1.0
    n = np.sqrt(1 + 4 * Av) / 4

    def f(rr):
        return np.sqrt(rr) * kv(n, np.sqrt(beta) * rr**2 / 2)

    h = 1e-5
    R1 = 40.0
    y0 = [f(R1), (f(R1 + h) - f(R1 - h)) / (2 * h)]
    s = solve_ivp(
        lambda t, y: [y[1], (Av / t**2 + beta * t**2) * y[0]],
        [R1, 6.0],
        y0,
        rtol=1e-11,
        atol=1e-30,
        dense_output=True,
    )
    backward = {}
    for rr in (30, 18, 12, 6):
        num = s.sol(rr)[0]
        backward["r=%g" % rr] = {"num": num, "claim": f(rr), "rel": abs(num - f(rr)) / abs(f(rr))}
        print("  backward ODE integration r=%g: num %.8e  claim %.8e" % (rr, num, f(rr)))

    record(
        "C2.ode_residual",
        "sqrt(r) K_nu(sqrt(beta) r^2 / 2) solves eps'' = (A / r^2 + beta r^2) eps (rel < 1e-12)",
        max(residuals.values()) < 1e-12,
        {"max_rel_residual": residuals, "threshold": 1e-12, "mpmath_dps": 40},
    )
    record(
        "C2.rho_nu_K_nu_collapse",
        "eps r^(2 nu - 1/2) / (rho^nu K_nu) is one constant over r 6 to 18 (spread < 1e-10)",
        max(v["rho_nu_K_nu_spread"] for v in collapse.values()) < 1e-10,
        {k: v["rho_nu_K_nu_spread"] for k, v in collapse.items()},
        could_fail=False,
        note="an identity once C2.ode_residual holds; the ODE checks carry the test",
    )
    record(
        "C2.backward_integration",
        "a backward integration from r = 40 lands on the closed form at r 30, 18, 12, 6",
        max(v["rel"] for v in backward.values()) < 1e-6,
        {"rows": backward, "threshold": 1e-6},
    )
    hit = [
        nm
        for nm, lo, hi in readings
        if _matches_31_39(lo, hi) or (lo == 0.0 and 0.31 <= hi <= 0.39)
    ]
    record(
        "C2.off_by_31_to_39_percent",
        "the form sqrt(rho) K_nu is off by 31 to 39 percent over the window r 6 to 18",
        len(hit) > 0,
        {"collapse": collapse, "readings_that_reproduce_the_figure": hit},
        note="not reproduced: A = 1 gives max / min - 1 = 13.8 percent, A = 3 gives 141.6 "
        "percent, and no normalization of the deviation range gives 31 to 39; only the window "
        "r 6 to 18 was tested",
    )


# ================= C3 and C4: the auditor's own Lorentz and projector tools =================
def rand_lorentz(rng, scale_boost=1.0, scale_rot=2.0):
    A = np.zeros((4, 4))
    for i in range(1, 4):
        v = rng.normal() * scale_boost
        A[0, i], A[i, 0] = v, -v
    for i in range(1, 4):
        for j in range(i + 1, 4):
            v = rng.normal() * scale_rot
            A[i, j], A[j, i] = v, -v
    X = ETA @ A  # X^T eta + eta X = 0
    return expm(X)


def riesz_B(N, return_all=False):
    w, V = np.linalg.eig(N)
    assert np.max(np.abs(w.imag)) < 1e-7, w
    w, V = w.real, V.real
    norms = np.array([V[:, i] @ ETA @ V[:, i] for i in range(4)])
    # timelike = most negative eta-norm (normalized)
    ig = int(np.argmin(norms / np.sum(V**2, axis=0)))
    assert norms[ig] < 0
    rest = [i for i in range(4) if i != ig]
    i1 = rest[int(np.argmax(w[rest]))]

    def P(i):
        return np.outer(V[:, i], ETA @ V[:, i]) / norms[i]

    Pg, P1 = P(ig), P(i1)
    P23 = np.eye(4) - Pg - P1
    B = P23 @ N @ P23 - 0.5 * np.trace(P23 @ N) * P23
    if return_all:
        return B, Pg, P1, P23, V[:, ig] / np.sqrt(-norms[ig]), w
    return B


def eucl_B(N):
    w, V = np.linalg.eig(N)
    w, V = w.real, V.real
    norms = np.array([V[:, i] @ ETA @ V[:, i] for i in range(4)])
    ig = int(np.argmin(norms / np.sum(V**2, axis=0)))
    rest = [i for i in range(4) if i != ig]
    i1 = rest[int(np.argmax(w[rest]))]
    u = V[:, ig] / np.linalg.norm(V[:, ig])
    n = V[:, i1] / np.linalg.norm(V[:, i1])
    P = np.eye(4) - np.outer(u, u) - np.outer(n, n)
    return P @ N @ P - 0.5 * np.trace(P @ N) * P


def single_projector(N):
    w, V = np.linalg.eig(N)
    w, V = w.real, V.real
    norms = np.array([V[:, i] @ ETA @ V[:, i] for i in range(4)])
    ig = int(np.argmin(norms / np.sum(V**2, axis=0)))
    rest = [i for i in range(4) if i != ig]
    rest = sorted(rest, key=lambda i: -w[i])
    i2 = rest[1]  # upper member of the pair
    return np.outer(V[:, i2], ETA @ V[:, i2]) / norms[i2]


def boostK(a):
    K = np.zeros((4, 4))
    K[0, 1:] = a
    K[1:, 0] = a
    return K


def run_case(a, eps, theta_fun, x, label, frame=None, show=True, info=None):
    a = a / np.linalg.norm(a)
    K = boostK(a)
    Ns = np.diag([-G8, 1, DELTA + eps, DELTA - eps])
    if frame is not None:  # optional global Lorentz frame for the split background
        Ns = frame @ Ns @ np.linalg.inv(frame)
    hh = 1e-5

    def fields(xx):
        L = expm(theta_fun(xx) * K)
        N = L @ Ns @ np.linalg.inv(L)
        return N

    N = fields(x)
    M = N @ ETA
    sym = np.abs(M - M.T).max()
    trs = [
        np.trace(np.linalg.matrix_power(N, p)) - np.trace(np.linalg.matrix_power(Ns, p))
        for p in range(1, 5)
    ]

    # 4th-order FD of B
    def Bf(xx):
        return riesz_B(fields(xx), True)

    B0, Pg, P1, P23, u, w = Bf(x)
    dB = (-Bf(x + 2 * hh)[0] + 8 * Bf(x + hh)[0] - 8 * Bf(x - hh)[0] + Bf(x - 2 * hh)[0]) / (
        12 * hh
    )
    thp = (
        -theta_fun(x + 2 * hh)
        + 8 * theta_fun(x + hh)
        - 8 * theta_fun(x - hh)
        + theta_fun(x - 2 * hh)
    ) / (12 * hh)
    a4 = np.concatenate([[0.0], a])
    # |B a|^2 : evaluate with the UNBOOSTED split B (frame of the background)
    Bs = riesz_B(Ns)
    Ba = Bs @ a4
    Ba2_eucl = Ba @ Ba
    Ba2_eta = Ba @ ETA @ Ba
    t_mixed = np.trace(dB @ dB)
    t_frob = np.trace(dB @ dB.T)
    Hlo = ETA + 2 * np.outer(ETA @ u, ETA @ u)
    Hup = ETA + 2 * np.outer(u, u)
    t_h = np.trace(dB.T @ Hlo @ dB @ Hup)
    if info is not None:
        info.update(
            {
                "theta": theta_fun(x),
                "theta_prime": thp,
                "M_asym": sym,
                "max_trace_drift": max(abs(t) for t in trs),
                "covariant_trace": t_mixed,
                "frobenius": t_frob,
                "h_metric": t_h,
                "claim": 2 * Ba2_eucl * thp**2,
            }
        )
    if show:
        print(
            "%s: a=%s eps=%.3f theta=%.3f theta'=%.3f"
            % (label, np.round(a, 3), eps, theta_fun(x), thp)
        )
        print("    M asym %.1e | max trace drift %.1e" % (sym, max(abs(t) for t in trs)))
        print(
            "    tr(dBdB)=%+.8e  claim -2|Ba|^2 th'^2=%+.8e (eucl |Ba|^2) / %+.8e (eta |Ba|^2)"
            % (t_mixed, -2 * Ba2_eucl * thp**2, -2 * Ba2_eta * thp**2)
        )
        print(
            "    Frobenius=%+.8e  h-metric=%+.8e  claim=+%.8e"
            % (t_frob, t_h, 2 * Ba2_eucl * thp**2)
        )
    return t_mixed, t_frob, t_h, 2 * Ba2_eucl * thp**2


def check_c3(rng):
    """Riesz cluster projector. rng is shared with check_c4 and check_extras (one stream)."""
    print("\n=== C3a: B on Lorentz-transformed vacua ===")
    Nvac = np.diag([-G8, 1, DELTA, DELTA])
    worst = 0
    worst_sym = 0
    rap = []
    for _ in range(400):
        L = rand_lorentz(rng)
        N = L @ Nvac @ np.linalg.inv(L)
        M = N @ ETA
        worst_sym = max(worst_sym, np.abs(M - M.T).max() / np.abs(M).max())
        B = riesz_B(N)
        worst = max(worst, np.abs(B).max() / np.abs(N).max())
        rap.append(np.arccosh(abs(L[0, 0])))
    print(
        "max |B|/|N| over 400 random SO(3,1): %.3e ; max rapidity %.2f ; M symmetry defect %.2e"
        % (worst, max(rap), worst_sym)
    )
    record(
        "C3a.vacuum_orbit",
        "the Riesz B vanishes on 400 random Lorentz images of the vacuum (|B| / |N| < 1e-8)",
        worst < 1e-8,
        {
            "max_B_over_N": worst,
            "max_rapidity": max(rap),
            "M_symmetry_defect": worst_sym,
            "threshold": 1e-8,
        },
    )

    print("=== C3d: Euclidean construction on the same ===")
    vals = []
    for _ in range(400):
        L = rand_lorentz(rng)
        N = L @ Nvac @ np.linalg.inv(L)
        vals.append(np.abs(eucl_B(N)).max())
    print(
        "Euclidean |B|_max: min %.3e median %.3e max %.3e"
        % (min(vals), np.median(vals), max(vals))
    )
    # pure rotations as control: Euclidean must vanish there
    vals_rot = []
    for _ in range(50):
        L = rand_lorentz(rng, scale_boost=0.0)
        N = L @ Nvac @ np.linalg.inv(L)
        vals_rot.append(np.abs(eucl_B(N)).max())
    print("control, pure rotations, Euclidean |B|_max: %.3e" % max(vals_rot))
    # small boost scaling
    small = {}
    for s in (1e-1, 1e-2, 1e-3):
        K = np.zeros((4, 4))
        K[0, 2] = K[2, 0] = 1
        L = expm(s * K)
        N = L @ Nvac @ np.linalg.inv(L)
        small["rapidity %.0e" % s] = {
            "euclid": np.abs(eucl_B(N)).max(),
            "riesz": np.abs(riesz_B(N)).max(),
        }
        print(
            "  boost rapidity %.0e along doublet axis: Euclid |B| = %.3e ; Riesz |B| = %.3e"
            % (s, np.abs(eucl_B(N)).max(), np.abs(riesz_B(N)).max())
        )
    record(
        "C3d.euclidean_not_invariant",
        "the Euclidean construction gives B != 0 on boosted vacua (median |B| > 1e-3) and "
        "B = 0 on rotated ones (< 1e-10)",
        np.median(vals) > 1e-3 and max(vals_rot) < 1e-10,
        {
            "euclid_B_min_median_max": [min(vals), np.median(vals), max(vals)],
            "pure_rotation_control_max": max(vals_rot),
            "boost_along_doublet_axis": small,
        },
        note="a boost along a doublet axis leaves the Euclidean B at zero: see the extras",
    )

    print("=== C3b: tr B^2 = 2 eps^2 on transformed splits ===")
    worst = 0
    for _ in range(400):
        eps = rng.uniform(-0.25, 0.25)
        L = rand_lorentz(rng)
        N = L @ np.diag([-G8, 1, DELTA + eps, DELTA - eps]) @ np.linalg.inv(L)
        B = riesz_B(N)
        worst = max(worst, abs(np.trace(B @ B) - 2 * eps**2) / (2 * eps**2))
    print("max rel error of tr B^2 vs 2 eps^2: %.3e" % worst)
    record(
        "C3b.trB2_on_split_orbit",
        "tr B^2 = 2 eps^2 on 400 random Lorentz images of the split, |eps| < 0.25 (rel < 1e-5)",
        worst < 1e-5,
        {"max_rel_error": worst, "threshold": 1e-5},
    )
    # a check that could fail: eps large enough that delta+eps > 1 (ordering swap)
    eps = 0.9
    N = np.diag([-G8, 1, DELTA + eps, DELTA - eps])
    B = riesz_B(N)
    trB2 = np.trace(B @ B)
    print(
        "ordering-swap case eps=0.9 (delta+eps=1.2 > 1): tr B^2 = %.4f vs 2eps^2 = %.4f"
        "  <-- 'largest spacelike' picks the wrong member" % (trB2, 2 * eps**2)
    )
    record(
        "C3b.ordering_swap",
        "tr B^2 = 2 eps^2 persists at eps = 0.9, where delta + eps passes the eigenvalue 1",
        abs(trB2 - 2 * eps**2) / (2 * eps**2) < 1e-5,
        {"eps": eps, "trB2": trB2, "two_eps2": 2 * eps**2},
        kind="scope",
        note="the 'largest spacelike eigenvalue' selection rule picks the wrong member once "
        "delta + eps > 1; the cluster read is valid for delta + |eps| < 1 only",
    )

    print("=== C3c: smoothness through coalescence on a generic line ===")
    L = rand_lorentz(rng)
    M0 = L @ Nvac @ np.linalg.inv(L) @ ETA
    M0 = 0.5 * (M0 + M0.T)
    X = rng.normal(size=(4, 4))
    X = X + X.T
    print("   t        |B(t)+B(-t)|/|B(t)|   |B(t)|/t     |P2(t)-P2(-t)|")
    line = {}
    for t in (1e-2, 1e-3, 1e-4, 1e-5, 1e-6):
        Bp, Bm = riesz_B((M0 + t * X) @ ETA), riesz_B((M0 - t * X) @ ETA)
        Pp, Pm = single_projector((M0 + t * X) @ ETA), single_projector((M0 - t * X) @ ETA)
        line["%.0e" % t] = {
            "even_part_over_B": np.abs(Bp + Bm).max() / np.abs(Bp).max(),
            "B_over_t": np.abs(Bp).max() / t,
            "single_projector_jump": np.abs(Pp - Pm).max(),
        }
        print(
            "  %.0e   %.3e            %.6f   %.4f"
            % (
                t,
                np.abs(Bp + Bm).max() / np.abs(Bp).max(),
                np.abs(Bp).max() / t,
                np.abs(Pp - Pm).max(),
            )
        )
    # analyticity: fit B(t) entries to a polynomial across t=0 and check residual
    ts = np.linspace(-2e-2, 2e-2, 41)
    Bs = np.array([riesz_B((M0 + t * X) @ ETA) for t in ts])
    res_B = 0
    for i in range(4):
        for j in range(4):
            cfit = np.polyfit(ts, Bs[:, i, j], 5)
            res_B = max(res_B, np.abs(np.polyval(cfit, ts) - Bs[:, i, j]).max())
    print(
        "degree-5 polynomial fit of B(t) across t=0: max residual %.3e (|B| scale %.3e)"
        % (res_B, np.abs(Bs).max())
    )
    Ps = np.array([single_projector((M0 + t * X) @ ETA) for t in ts])
    res_P = 0
    for i in range(4):
        for j in range(4):
            cfit = np.polyfit(ts, Ps[:, i, j], 5)
            res_P = max(res_P, np.abs(np.polyval(cfit, ts) - Ps[:, i, j]).max())
    print("same fit for the single projector: max residual %.3e (O(1) => jump)" % res_P)
    even6, even5 = line["1e-06"]["even_part_over_B"], line["1e-05"]["even_part_over_B"]
    record(
        "C3c.smooth_through_coalescence",
        "B(t) is smooth across the coalescence on a generic line (even part falls linearly in "
        "t, degree-5 fit residual < 1e-2 of scale) while a single projector jumps by O(1)",
        even6 < 1e-3
        and 5 < even5 / even6 < 20
        and res_B / np.abs(Bs).max() < 1e-2
        and res_P > 0.1
        and line["1e-06"]["single_projector_jump"] > 0.5,
        {
            "line": line,
            "fit_residual_B": res_B,
            "B_scale": np.abs(Bs).max(),
            "fit_residual_single_projector": res_P,
        },
    )


def check_c4(rng):
    """Boost ripple witness. Continues the random stream of check_c3."""
    print("\n=== C4: boost ripple witness ===")

    def th1(x):
        return 0.3 * np.sin(1.7 * x) + 0.1 * x

    def big(x):
        return 1.5 * np.sin(x) + 0.8

    named = {}
    for a, eps, fun, x, label in (
        (np.array([0, 1.0, 0]), 0.1, th1, 0.4, "doublet axis 2"),
        (np.array([0, 0, 1.0]), 0.1, th1, 0.4, "doublet axis 3"),
        (np.array([1.0, 0, 0]), 0.1, th1, 0.4, "the '1' axis (not in doublet plane)"),
        (np.array([0, 1.0, 1.0]), 0.1, th1, 0.4, "doublet diagonal (B a rotates a)"),
        (np.array([0, 1.0, 0]), 0.1, big, 0.3, "LARGE theta, doublet axis"),
    ):
        info = {}
        run_case(a, eps, fun, x, label, info=info)
        named[label] = info
    print("--- random sweep ---")
    worst_d, worst_h, worst_f, anypos = 0, 0, 0, 0
    for k in range(300):
        a = rng.normal(size=3)
        eps = rng.uniform(-0.25, 0.25)
        c1, c2, c3 = rng.normal(size=3)

        def th(x, c1=c1, c2=c2, c3=c3):
            return 0.5 * c1 * np.sin(c2 * x + c3)

        x = rng.uniform(-2, 2)
        tm, tf, thm, cl = run_case(a, eps, th, x, "rnd", show=False)
        if cl > 1e-14:
            worst_d = max(worst_d, abs(tm + cl) / cl)
            worst_h = max(worst_h, abs(thm - cl) / cl)
            worst_f = max(worst_f, abs(tf - cl) / cl)
        if tm > 1e-12:
            anypos += 1
    print(
        "random (a, eps, theta): max rel err (d) %.2e ; h-metric %.2e ; Frobenius %.2e ; "
        "cases with tr(dBdB)>0: %d" % (worst_d, worst_h, worst_f, anypos)
    )
    print(
        "--- split background in a random global Lorentz frame "
        "(a is then not the frame's own axis) ---"
    )
    worst = 0
    pos = 0
    for k in range(200):
        F = rand_lorentz(rng)
        a = rng.normal(size=3)
        eps = rng.uniform(-0.25, 0.25)
        tm, tf, thm, cl = run_case(
            a, eps, lambda x: 0.4 * np.sin(x), 0.2, "frame", frame=F, show=False
        )
        pos += tm > 1e-12
        worst = max(worst, abs(tm + cl) / max(cl, 1e-14))
    print(
        "generic-frame background: tr(dBdB) > 0 in %d of 200 cases; formula -2|Ba|^2 th'^2 "
        "(B of the background, Euclidean norm) max rel err %.2e" % (pos, worst)
    )

    live = [v for v in named.values() if v["claim"] > 1e-14]
    dead = [v for v in named.values() if v["claim"] <= 1e-14]
    record(
        "C4.covariant_trace_named_cases",
        "tr(dB dB) = -2 |B a|^2 theta'^2 on the five named cases (rel < 1e-5; zero off the "
        "doublet plane)",
        max(abs(v["covariant_trace"] + v["claim"]) / v["claim"] for v in live) < 1e-5
        and max(abs(v["covariant_trace"]) for v in dead) < 1e-12,
        {"cases": named, "threshold": 1e-5},
    )
    record(
        "C4.covariant_trace_sweep",
        "the same closed form on 300 random (a, eps, theta) with the background at rest "
        "(rel < 1e-4), never positive",
        worst_d < 1e-4 and anypos == 0,
        {"max_rel_error": worst_d, "positive_cases": anypos, "threshold": 1e-4},
    )
    record(
        "C4.h_metric_sweep",
        "the h-metric contraction equals +2 |B a|^2 theta'^2 on the same sweep (rel < 1e-4)",
        worst_h < 1e-4,
        {"max_rel_error": worst_h, "threshold": 1e-4},
    )
    big_case = named["LARGE theta, doublet axis"]
    record(
        "C4.frobenius_finite_theta",
        "the Frobenius contraction equals +2 |B a|^2 theta'^2 at finite theta (rel < 1e-4)",
        worst_f < 1e-4,
        {
            "max_rel_error_sweep": worst_f,
            "large_theta_case_frobenius_over_claim": big_case["frobenius"] / big_case["claim"],
            "threshold": 1e-4,
        },
        note="refuted at finite theta: the Frobenius contraction is frame dependent and grows "
        "as cosh(4 theta) (extras, EX.c4e); it equals the claim at theta = 0 only",
    )
    record(
        "C4.generic_frame_background",
        "the closed form, with B of the background, holds when the split background sits in "
        "a random global Lorentz frame",
        pos == 0 and worst < 1e-4,
        {"positive_cases_of_200": int(pos), "max_rel_error": worst},
        kind="scope",
        note="with a boosted background the lab-frame boost generator is a boost plus a "
        "rotation in the rest frame, and the rotation part contributes with the other sign",
    )


# ================= C5: the ripple on the stored relaxed field =================
def check_c5():
    """Production energy and the auditor's own energy on the stored field. Returns the vacuum
    control for check_vacuum_artifact, or None when skipped."""
    print("\n=== C5: boost ripple on the stored relaxed field ===")
    n, L, c = 32, 48.0, 1e-3
    if not os.path.exists(os.path.join(HERE, PRODUCTION)):
        record("C5.stored_field", "C5 needs the production module", "SKIPPED",
               {"reason": "production module absent", "module": PRODUCTION})  # fmt: skip
        return None
    CS = _load("m5_32_r23_1_cscan", PRODUCTION)
    tag = CS.job_tag(
        {"seed": "rad", "delta": 0.3, "w1s": 25.0, "c": 1e-3, "n": 32, "L": 48.0, "src": None}
    )
    rel_path = "../data/m5_32_r23_1/" + tag + ".npz"
    path = os.path.join(CS.OUT_NPZ, tag + ".npz")
    if not os.path.exists(path):
        record("C5.stored_field", "C5 needs the stored relaxed field", "SKIPPED",
               {"reason": "stored field file absent", "file": rel_path})  # fmt: skip
        return None
    cfg = CS.R21.cfg_of(32, 48.0, 8.0, 0.3)
    p = CS.R21.params_of(8.0, 0.3)
    roots = CS.R0.roots_of(cfg, degenerate=True)
    pot = ("v4", roots, CS.W1 * 25.0)
    h = cfg["h"]
    print("field:", rel_path, "h =", h, "roots =", roots, "W1 =", CS.W1)
    M0 = np.load(path)["M"].astype(np.float64)
    print(
        "shape", M0.shape,
        "max |M_0a| (time-space entries of the stored field) =", np.abs(M0[..., 0, 1:]).max(),
    )  # fmt: skip
    X, Y, Z = CS.B3.coords(n, h)
    r = np.sqrt(X**2 + Y**2 + Z**2)
    WW = np.outer(np.diag(ETA), np.diag(ETA))

    def E_prod(M):
        return CS.energy_grad(M, cfg, p, pot, c, need_grad=False)[0]

    # ---------- own energy (fwd/bwd averaged first differences, eta commutator, V4, -cL) ------
    def diff(f, ax, kind):
        out = np.zeros_like(f)
        a = [slice(None)] * f.ndim
        b = [slice(None)] * f.ndim
        a[ax] = slice(1, None)
        b[ax] = slice(0, -1)
        d = (f[tuple(a)] - f[tuple(b)]) / h
        if kind == "fwd":
            out[tuple(b)] = d
        else:
            out[tuple(a)] = d
        return out

    def my_parts(M):
        ec = 0.0
        for kind in ("fwd", "bwd"):
            A = [diff(M, ax, kind) for ax in range(3)]
            for i in range(3):
                for j in range(i + 1, 3):
                    F = A[i] @ ETA @ A[j] - A[j] @ ETA @ A[i]
                    ec += 0.5 * 4.0 * np.sum(F * F * WW)
        N = M @ ETA
        P = np.broadcast_to(np.eye(4), N.shape)
        v = 0.0
        lin = 0.0
        distinct = sorted(set(float(q) for q in roots))
        poly = np.poly(distinct)[::-1]
        for k in range(1, 5):
            P = P @ N
            t = np.einsum("...kk->...", P) - sum(q**k for q in roots)
            v = v + t**2
            lin = lin + poly[k - 1] / k * t
        return h**3 * ec, h**3 * CS.W1 * 25.0 * v.sum(), -c * h**3 * lin.sum()

    def boost_field(theta, axis=2):
        """Lam = exp(theta K), K mixes index 0 with spatial index `axis` (1=x,2=y,3=z)."""
        ch, sh = np.cosh(theta), np.sinh(theta)
        Lam = np.broadcast_to(np.eye(4), theta.shape + (4, 4)).copy()
        Lam[..., 0, 0] = ch
        Lam[..., axis, axis] = ch
        Lam[..., 0, axis] = sh
        Lam[..., axis, 0] = sh
        return Lam

    def rot_field(theta, i=1, j=2):
        cs, sn = np.cos(theta), np.sin(theta)
        R = np.broadcast_to(np.eye(4), theta.shape + (4, 4)).copy()
        R[..., i, i] = cs
        R[..., j, j] = cs
        R[..., i, j] = -sn
        R[..., j, i] = sn
        return R

    mask = r < (L / 2 - 3 * h)

    def theta_of(A, k, renv=12.0, phase="sin"):
        car = np.sin(k * X) if phase == "sin" else np.cos(k * X)
        return A * car * np.exp(-((r / renv) ** 4)) * mask

    def apply(Lam, M):
        return Lam @ M @ np.swapaxes(Lam, -1, -2)

    def predicted(theta, axis=2):
        """small-theta continuum-form prediction on the lattice stencil:
        D_i = A_i + (d_i theta) Q, Q = K M + M K; dE = E_curv[D] - E_curv[A]."""
        K = np.zeros((4, 4))
        K[0, axis] = K[axis, 0] = 1
        Q = K @ M0 + M0 @ K
        tot = 0.0
        cross = 0.0
        for kind in ("fwd", "bwd"):
            A = [diff(M0, ax, kind) for ax in range(3)]
            T = [diff(theta, ax, kind) for ax in range(3)]
            for i in range(3):
                for j in range(i + 1, 3):
                    F = A[i] @ ETA @ A[j] - A[j] @ ETA @ A[i]
                    Ci = A[i] @ ETA @ Q - Q @ ETA @ A[i]
                    Cj = A[j] @ ETA @ Q - Q @ ETA @ A[j]
                    Gm = T[j][..., None, None] * Ci - T[i][..., None, None] * Cj
                    tot += 0.5 * 4 * np.sum(Gm * Gm * WW)
                    cross += 0.5 * 4 * 2 * np.sum(F * Gm * WW)
        return h**3 * tot, h**3 * cross

    E0 = E_prod(M0)
    parts0 = my_parts(M0)
    pp = CS.R20.energy_parts(M0, cfg, p, pot)
    print(
        "E0 production = %.6f ; production parts: E_curv %.6f V %.6f" % (E0, pp["E_curv"], pp["V"])
    )
    print("E0 mine       = %.6f ; my parts: curv %.6f V4 %.6f -cL %.6f" % (sum(parts0), *parts0))

    rows = []
    print("\n--- claimed cases + scaling scan (boost along y) ---")
    print(
        "   A        k*h/pi    dE_prod       dE_mine      d_curv       d_V4        d_lin"
        "      pred_quad    pred_cross   dE/(A^2)"
    )
    cases = [(0.02, 1 / 4), (0.02, 1 / 2), (0.05, 1 / 2)]
    for kk in (1 / 2, 1 / 4, 1 / 8, 1 / 16, 1 / 32, 0.0):
        for A in (0.05, 0.02, 0.005, 0.001, 1e-4):
            if (A, kk) not in cases:
                cases.append((A, kk))
    for A, kk in cases:
        k = kk * np.pi / h
        th = theta_of(A, k, phase="sin" if kk > 0 else "cos")
        Mr = apply(boost_field(th), M0)
        asym = np.abs(Mr - np.swapaxes(Mr, -1, -2)).max()
        dE = E_prod(Mr) - E0
        pm = my_parts(Mr)
        pq, pc = predicted(th)
        print(
            " %.0e   %.5f   %+.6e  %+.6e  %+.4e  %+.2e  %+.2e  %+.6e  %+.2e  %+.5e"
            % (
                A,
                kk,
                dE,
                sum(pm) - sum(parts0),
                pm[0] - parts0[0],
                pm[1] - parts0[1],
                pm[2] - parts0[2],
                pq,
                pc,
                dE / A**2,
            )
        )
        rows.append(
            {
                "A": A,
                "kh_over_pi": kk,
                "dE_prod": dE,
                "dE_mine": sum(pm) - sum(parts0),
                "d_curv": pm[0] - parts0[0],
                "d_V4": pm[1] - parts0[1],
                "d_lin": pm[2] - parts0[2],
                "pred_quad": pq,
                "pred_cross": pc,
                "asym": float(asym),
            }
        )

    print("\n--- controls ---")
    th = np.full_like(r, 0.3)
    d_uniform = E_prod(apply(boost_field(th), M0)) - E0
    print(
        "uniform boost theta=0.3 everywhere (incl. pinned shell): dE = %+.3e  "
        "(invariance control; must be ~0)" % d_uniform
    )
    rotation = {}
    for A, kk in ((0.02, 1 / 4), (0.02, 1 / 16), (0.005, 1 / 4)):
        th = theta_of(A, kk * np.pi / h)
        dEr = E_prod(apply(rot_field(th, 1, 2), M0)) - E0
        rotation["A=%g kh/pi=%g" % (A, kk)] = dEr
        print(
            "ROTATION ripple (x-y plane) A=%.3f kh/pi=%.4f: dE = %+.6e   "
            "(relaxed minimum => expected >= 0)" % (A, kk, dEr)
        )
    variants = {}
    for axis, nm in ((1, "x (along the carrier)"), (3, "z")):
        th = theta_of(0.02, np.pi / (4 * h))
        variants["boost along " + nm] = E_prod(apply(boost_field(th, axis), M0)) - E0
        print("boost along %s, A=0.02 kh/pi=1/4: dE = %+.6e" % (nm, variants["boost along " + nm]))
    # envelope width scan at k=0 (pure long-wavelength bump)
    for renv in (6.0, 9.0, 12.0):
        th = theta_of(0.02, 0.0, renv=renv, phase="cos")
        variants["k=0 bump renv=%.0f" % renv] = E_prod(apply(boost_field(th), M0)) - E0
        print(
            "k=0 bump, renv=%.0f, A=0.02: dE = %+.6e"
            % (renv, variants["k=0 bump renv=%.0f" % renv])
        )
    # far-field ripple: confined to a shell where the background gradient is small
    for rc, wd in ((6.0, 3.0), (12.0, 3.0), (16.0, 2.5)):
        th = 0.02 * np.sin(np.pi / (4 * h) * X) * np.exp(-(((r - rc) / wd) ** 2)) * mask
        key = "shell r=%.0f width %.1f" % (rc, wd)
        variants[key] = E_prod(apply(boost_field(th), M0)) - E0
        print(
            "shell ripple at r=%.0f width %.1f, A=0.02 kh/pi=1/4: dE = %+.6e"
            % (rc, wd, variants[key])
        )
    # vacuum background control: uniform degenerate vacuum -> C_i = 0 -> dE must vanish
    Mv = np.broadcast_to(np.diag([8.0, 1, 0.3, 0.3]), M0.shape).copy()
    Ev = E_prod(Mv)
    th = theta_of(0.05, np.pi / (2 * h))
    d_vac = E_prod(apply(boost_field(th), Mv)) - Ev
    print(
        "uniform VACUUM diag(8,1,.3,.3) + ripple A=0.05 kh/pi=1/2: E_vac %.3e  dE = %+.3e"
        % (Ev, d_vac)
    )
    Ms = np.broadcast_to(np.diag([8.0, 1, 0.4, 0.2]), M0.shape).copy()
    Es = E_prod(Ms)
    d_split = E_prod(apply(boost_field(th), Ms)) - Es
    print("uniform SPLIT diag(8,1,.4,.2) + same ripple: dE = %+.3e" % d_split)

    by = {(row["A"], row["kh_over_pi"]): row for row in rows}
    record(
        "C5.baseline_energy",
        "own lattice energy equals the production energy on the stored field, and both equal "
        "the claimed E0 (rel < 1e-9)",
        abs(sum(parts0) - E0) / E0 < 1e-9 and abs(E0 - CLAIM_E0) / CLAIM_E0 < 1e-9,
        {
            "field": rel_path,
            "E0_production": E0,
            "E0_own": sum(parts0),
            "E0_claimed": CLAIM_E0,
            "own_parts_curv_V4_lin": list(parts0),
            "production_parts": {"E_curv": pp["E_curv"], "V": pp["V"]},
            "max_abs_time_space_entries": np.abs(M0[..., 0, 1:]).max(),
        },
    )
    record(
        "C5.claimed_cases",
        "the three claimed dE reproduce with the production energy (rel < 1e-8)",
        max(abs(by[k]["dE_prod"] - v) / abs(v) for k, v in CLAIM_DE.items()) < 1e-8,
        {
            "rows": [
                {"A": k[0], "kh_over_pi": k[1], "dE_production": by[k]["dE_prod"], "claimed": v}
                for k, v in CLAIM_DE.items()
            ],
            "threshold": 1e-8,
        },
    )
    record(
        "C5.own_energy_agrees",
        "own energy gives the same dE as the production energy on all 33 rows (rel < 1e-6)",
        max(abs(row["dE_mine"] - row["dE_prod"]) / abs(row["dE_prod"]) for row in rows) < 1e-6,
        {"rows": rows, "threshold": 1e-6},
    )
    record(
        "C5.negative_on_every_row",
        "dE < 0 on all 33 (A, k) rows, with the symmetry of M kept (asymmetry < 1e-12)",
        all(row["dE_prod"] < 0 for row in rows) and max(row["asym"] for row in rows) < 1e-12,
        {
            "max_dE": max(row["dE_prod"] for row in rows),
            "min_dE": min(row["dE_prod"] for row in rows),
            "max_asym": max(row["asym"] for row in rows),
        },
    )
    record(
        "C5.curvature_term_alone",
        "the change sits in the curvature term: V4 and the linear member move by < 1e-10",
        max(max(abs(row["d_V4"]), abs(row["d_lin"])) for row in rows) < 1e-10,
        {
            "max_abs_d_V4": max(abs(row["d_V4"]) for row in rows),
            "max_abs_d_lin": max(abs(row["d_lin"]) for row in rows),
        },
    )
    quad = {}
    pred = {}
    for kk in (1 / 2, 1 / 4, 1 / 8, 1 / 16, 1 / 32, 0.0):
        a3, a4 = by[(0.001, kk)], by[(1e-4, kk)]
        quad["kh/pi=%g" % kk] = abs(a3["dE_prod"] / 1e-6 - a4["dE_prod"] / 1e-8) / abs(
            a4["dE_prod"] / 1e-8
        )
        pred["kh/pi=%g" % kk] = abs(a4["pred_quad"] - a4["dE_prod"]) / abs(a4["dE_prod"])
    record(
        "C5.quadratic_in_A",
        "dE / A^2 is one constant at small A for every k (A = 1e-3 against 1e-4, rel < 1e-3): "
        "a second-order instability, not a finite-amplitude effect",
        max(quad.values()) < 1e-3,
        {"rel_change_of_dE_over_A2": quad, "threshold": 1e-3},
    )
    record(
        "C5.second_order_prediction",
        "own small-theta second-order form on the lattice stencil predicts dE at A = 1e-4 "
        "(rel < 2e-2)",
        max(pred.values()) < 2e-2,
        {
            "rel_error": pred,
            "max_abs_cross_term": max(abs(row["pred_cross"]) for row in rows),
            "threshold": 2e-2,
        },
    )
    record(
        "C5.uniform_boost_invariance",
        "a uniform boost of the whole field leaves the production energy unchanged (< 1e-9)",
        abs(d_uniform) < 1e-9,
        {"dE": d_uniform, "threshold": 1e-9},
    )
    record(
        "C5.rotation_ripple_control",
        "a rotation ripple of the same shape raises the energy (the field is a minimum "
        "against rotations)",
        all(v > 0 for v in rotation.values()),
        {"dE": rotation},
    )
    record(
        "C5.axis_envelope_shell_variants",
        "dE stays negative for the boost along x and z, for a k = 0 bump at three widths and "
        "for shell-confined ripples at r = 6, 12, 16",
        all(v < 0 for v in variants.values()),
        {"dE": variants},
    )
    record(
        "C5.uniform_background_control",
        "on a uniform vacuum, where the continuum form gives zero, the lattice residue is "
        "below 1e-3 of the field's dE for the same ripple",
        abs(d_vac) / abs(by[(0.05, 1 / 2)]["dE_prod"]) < 1e-3,
        {
            "dE_uniform_vacuum": d_vac,
            "dE_uniform_split": d_split,
            "dE_field_same_ripple": by[(0.05, 1 / 2)]["dE_prod"],
            "E_vacuum": Ev,
        },
    )
    return {"dE_uniform_vacuum_A0.05_k0.5": d_vac}


# ================= vacuum artifact =================
def check_vacuum_artifact(c5=None):
    """Size and scaling of the pure LATTICE artifact (ripple on a UNIFORM vacuum, where the
    continuum form gives exactly zero because d_i M = theta_i Lam Q Lam^T are all parallel).
    Own curvature energy only; the production stack is not used here."""
    print("\n=== vacuum artifact: the ripple on a uniform vacuum ===")
    n, L = 32, 48.0
    h = L / n
    x = (np.arange(n) - (n - 1) / 2) * h
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    r = np.sqrt(X**2 + Y**2 + Z**2)
    eta = np.array([-1.0, 1, 1, 1])
    WW = np.outer(eta, eta)

    def diff(f, ax, kind):
        out = np.zeros_like(f)
        a = [slice(None)] * f.ndim
        b = [slice(None)] * f.ndim
        a[ax] = slice(1, None)
        b[ax] = slice(0, -1)
        d = (f[tuple(a)] - f[tuple(b)]) / h
        out[tuple(b) if kind == "fwd" else tuple(a)] = d
        return out

    def ecurv(M):
        e = 0.0
        for kind in ("fwd", "bwd"):
            A = [diff(M, ax, kind) for ax in range(3)]
            for i in range(3):
                for j in range(i + 1, 3):
                    F = A[i] @ ETA @ A[j] - A[j] @ ETA @ A[i]
                    e += 2.0 * np.sum(F * F * WW)
        return h**3 * e

    def boost(theta):
        Lam = np.broadcast_to(np.eye(4), theta.shape + (4, 4)).copy()
        Lam[..., 0, 0] = Lam[..., 2, 2] = np.cosh(theta)
        Lam[..., 0, 2] = Lam[..., 2, 0] = np.sinh(theta)
        return Lam

    Mv = np.broadcast_to(np.diag([8.0, 1, 0.3, 0.3]), (n, n, n, 4, 4)).copy()
    mask = r < L / 2 - 3 * h
    rows = {}
    print("   A      kh/pi   dE(vacuum)      dE/A^4        dE/A^6")
    for kk in (0.5, 0.25, 0.125):
        for A in (0.1, 0.05, 0.02, 0.01):
            th = A * np.sin(kk * np.pi / h * X) * np.exp(-((r / 12) ** 4)) * mask
            Lm = boost(th)
            dE = ecurv(Lm @ Mv @ np.swapaxes(Lm, -1, -2))
            rows[(A, kk)] = dE
            print(" %.3f  %.4f  %+.4e  %+.4e  %+.4e" % (A, kk, dE, dE / A**4, dE / A**6))
    a6 = {}
    a4 = {}
    for kk in (0.5, 0.25, 0.125):
        hi, lo = rows[(0.02, kk)], rows[(0.01, kk)]
        a6["kh/pi=%g" % kk] = abs(hi / 0.02**6 - lo / 0.01**6) / abs(lo / 0.01**6)
        a4["kh/pi=%g" % kk] = (hi / 0.02**4) / (lo / 0.01**4)
    table = [
        {"A": k[0], "kh_over_pi": k[1], "dE": v, "dE_over_A6": v / k[0] ** 6}
        for k, v in rows.items()
    ]
    record(
        "VA.scales_as_A6",
        "the ripple energy on the uniform vacuum is a lattice artifact of order A^6 (dE / A^6 "
        "constant to 1 percent between A = 0.02 and 0.01), so it cannot feed an A^2 signal",
        max(a6.values()) < 1e-2,
        {"rel_change_of_dE_over_A6": a6, "ratio_of_dE_over_A4": a4, "rows": table},
    )
    ref = abs(CLAIM_DE[(0.05, 1 / 2)])
    record(
        "VA.negligible_against_c5",
        "at A = 0.05, k h / pi = 1/2 the artifact is below 1e-3 of the claimed C5 dE",
        abs(rows[(0.05, 0.5)]) / ref < 1e-3,
        {
            "dE_vacuum": rows[(0.05, 0.5)],
            "claimed_dE_field": -ref,
            "ratio": abs(rows[(0.05, 0.5)]) / ref,
        },
    )
    if c5 is not None:
        prod = c5["dE_uniform_vacuum_A0.05_k0.5"]
        record(
            "VA.matches_production_energy",
            "own curvature energy and the production energy agree on the vacuum artifact "
            "(rel < 1e-6)",
            abs(prod - rows[(0.05, 0.5)]) / abs(prod) < 1e-6,
            {"own": rows[(0.05, 0.5)], "production": prod},
        )


# ================= C6 pointwise: Cartesian jet calculus =================
def pointwise_machinery():
    """Pointwise first and second variation of the curvature density
    e = 4 sum_{i<j} tr(F_ij F_ij^T), F_ij = [d_i M, d_j M], around the hedgehog, derived in
    CARTESIAN coordinates (sympy differentiates the Cartesian expressions of the frame
    tensors; no spherical frame-derivative table is used).

    Fields: psi = a T1 + b T2 (doublet) + d1 D1 + d2 D2 (director) + s Srr + u Utt (massive)
      D1 = rhat that^T + that rhat^T, D2 = rhat phat^T + phat rhat^T,
      Srr = rhat rhat^T, Utt = that that^T + phat phat^T.
    Jets per field f: f, f_r, f_t, f_p
      (d_k f = rhat_k f_r + that_k f_t / r + phat_k f_p / (r sin))."""
    if "pointwise" in _CACHE:
        return _CACHE["pointwise"]
    x, y, z = sp.symbols("x y z", real=True)
    X = sp.Matrix([x, y, z])
    r = sp.sqrt(x**2 + y**2 + z**2)
    rho = sp.sqrt(x**2 + y**2)
    rh = X / r
    ph = sp.Matrix([-y, x, 0]) / rho
    th = sp.Matrix([x * z, y * z, -(rho**2)]) / (r * rho)  # d rhat / d theta
    delta = sp.Rational(3, 10)

    T = {
        "a": th * th.T - ph * ph.T,
        "b": th * ph.T + ph * th.T,
        "d1": rh * th.T + th * rh.T,
        "d2": rh * ph.T + ph * rh.T,
        "s": rh * rh.T,
        "u": th * th.T + ph * ph.T,
    }
    M0 = delta * sp.eye(3) + (1 - delta) * rh * rh.T
    coords = (x, y, z)
    dM0 = [M0.diff(c) for c in coords]
    dT = {k: [v.diff(c) for c in coords] for k, v in T.items()}

    lam = {k: sp.lambdify((x, y, z), v, "numpy") for k, v in T.items()}
    lam_d = {k: [sp.lambdify((x, y, z), m, "numpy") for m in v] for k, v in dT.items()}
    lam_dM0 = [sp.lambdify((x, y, z), m, "numpy") for m in dM0]
    lam_fr = [sp.lambdify((x, y, z), v, "numpy") for v in (rh, th, ph)]

    FIELDS = ["a", "b", "d1", "d2", "s", "u"]
    JETS = ["", "_r", "_t", "_p"]
    NAMES = [f + j for f in FIELDS for j in JETS]

    def jets_at(R, TH, PH):
        """list of 24 'gradient tensors' G[k] (3x3 each, k = x,y,z): the d_k psi produced by
        a unit jet."""
        pt = (R * np.sin(TH) * np.cos(PH), R * np.sin(TH) * np.sin(PH), R * np.cos(TH))
        rhat, that, phat = [np.array(f(*pt), float).ravel() for f in lam_fr]
        out = []
        for f in FIELDS:
            Tm = np.array(lam[f](*pt), float)
            dTm = [np.array(g(*pt), float) for g in lam_d[f]]
            out.append(np.array(dTm))  # jet f: d_k psi = d_k T
            for vec, scale in ((rhat, 1.0), (that, 1.0 / R), (phat, 1.0 / (R * np.sin(TH)))):
                out.append(np.array([vec[k] * scale * Tm for k in range(3)]))
        A0 = np.array([np.array(g(*pt), float) for g in lam_dM0])
        return A0, out

    def comm(A, B):
        return A @ B - B @ A

    def variations(A0, G):
        nj = len(G)
        F0 = [comm(A0[i], A0[j]) for i, j in PAIRS]
        e0 = 4 * sum(np.sum(F * F) for F in F0)
        F1 = [[comm(A0[i], g[j]) + comm(g[i], A0[j]) for i, j in PAIRS] for g in G]
        lin = np.array([8 * sum(np.sum(F0[q] * F1[a][q]) for q in range(3)) for a in range(nj)])
        H = np.zeros((nj, nj))
        for a in range(nj):
            for b in range(nj):
                v = 0.0
                for q, (i, j) in enumerate(PAIRS):
                    v += np.sum(F1[a][q] * F1[b][q])
                    F2 = comm(G[a][i], G[b][j]) + comm(G[b][i], G[a][j])
                    v += np.sum(F0[q] * F2)
                H[a, b] = 8 * v  # H_ab = d^2 e / d v_a d v_b  (e2 = 1/2 v^T H v)
        return e0, lin, H

    def claimed_H(R, TH):
        """the claimant's table, scaled back by (1-delta)^2 / r^4, keyed by jet-name pairs."""
        s, c = np.sin(TH), np.cos(TH)
        t = {}
        t[("a_r", "a_r")] = t[("b_r", "b_r")] = 16 * R**2
        t[("a_t", "a_t")] = t[("b_t", "b_t")] = 8
        t[("a_p", "a_p")] = t[("b_p", "b_p")] = 8 / s**2
        t[("a", "a")] = t[("b", "b")] = 32 / s**2 - 48
        t[("a", "a_t")] = t[("b", "b_t")] = 48 * c / s
        t[("a", "b_p")] = 16 * c / s**2
        t[("a_p", "b")] = -16 * c / s**2
        t[("a_t", "b_p")] = 24 / s
        t[("a_p", "b_t")] = -24 / s
        return t

    def fd_check(A0, G, H, lin, rng):
        """independent of the bilinear algebra: evaluate the full quartic density on random
        jets and take a numerical second directional derivative."""

        def dens(v):
            D = A0 + sum(v[a] * G[a] for a in range(len(G)))
            return 4 * sum(np.sum(comm(D[i], D[j]) ** 2) for i, j in PAIRS)

        worst = 0
        for _ in range(5):
            v = rng.normal(size=len(G))
            hh = 1e-3
            # e(t v) is a quartic in t: five-point stencils are exact up to rounding
            f = [dens(k * hh * v) for k in (-2, -1, 0, 1, 2)]
            d1 = (f[0] - 8 * f[1] + 8 * f[3] - f[4]) / (12 * hh)
            d2 = (-f[0] + 16 * f[1] - 30 * f[2] + 16 * f[3] - f[4]) / (12 * hh**2)
            worst = max(
                worst, abs(d1 - lin @ v) / (abs(d1) + 1e-12), abs(d2 - v @ H @ v) / abs(d2)
            )
        return worst

    m = SimpleNamespace(
        delta=float(delta),
        NAMES=NAMES,
        jets_at=jets_at,
        variations=variations,
        claimed_H=claimed_H,
        fd_check=fd_check,
    )
    _CACHE["pointwise"] = m
    return m


def check_c6_pointwise():
    """C6 (a), (b): linear terms and the jet Hessian at random points."""
    print("\n=== C6 pointwise: first and second variation in Cartesian jets ===")
    m = pointwise_machinery()
    NAMES = m.NAMES
    rng = np.random.default_rng(7)
    dl = m.delta
    doublet = [k for k, n in enumerate(NAMES) if n[0] in "ab"]
    direc = [k for k, n in enumerate(NAMES) if n.startswith("d")]
    massive = [k for k, n in enumerate(NAMES) if n[0] in "su"]
    worst_lin_doublet = 0
    worst_claim = 0
    worst_claim2 = 0
    worstfd = 0
    first = {}
    with np.printoptions(linewidth=200, precision=4, suppress=True):
        for trial in range(8):
            R, TH, PH = rng.uniform(0.5, 9), rng.uniform(0.2, 2.9), rng.uniform(0, 6.28)
            A0, G = m.jets_at(R, TH, PH)
            e0, lin, H = m.variations(A0, G)
            worstfd = max(worstfd, m.fd_check(A0, G, H, lin, rng))
            sc = R**4 / (1 - dl) ** 2
            Hs = H * sc
            worst_lin_doublet = max(worst_lin_doublet, np.abs(lin[doublet]).max() * sc)
            cl = m.claimed_H(R, TH)
            Hc = np.zeros((8, 8))
            dn = [NAMES[k] for k in doublet]
            for (p, q), v in cl.items():
                i, j = dn.index(p), dn.index(q)
                Hc[i, j] = Hc[j, i] = v
            Hd = Hs[np.ix_(doublet, doublet)]
            worst_claim = max(worst_claim, np.abs(Hd - Hc).max() / np.abs(Hc).max())
            worst_claim2 = max(worst_claim2, np.abs(Hd - 2 * Hc).max() / np.abs(2 * Hc).max())
            if trial == 0:
                lin_scaled = {n: float(v) for n, v in zip(NAMES, np.round(lin * sc, 6))}
                first = {
                    "point_r_theta_phi": [R, TH, PH],
                    "e0_scaled": e0 * sc,
                    "max_abs_doublet_x_director": np.abs(Hs[np.ix_(doublet, direc)]).max(),
                    "max_abs_doublet_x_massive": np.abs(Hs[np.ix_(doublet, massive)]).max(),
                    "nonzero_linear_terms_scaled": {
                        n: v for n, v in lin_scaled.items() if abs(v) > 1e-9
                    },
                }
                print(
                    "point r=%.3f theta=%.3f phi=%.3f ; e0 r^4/(1-d)^2 = %.6f"
                    % (R, TH, PH, e0 * sc)
                )
                print("doublet jets:", dn)
                print("MY Hessian d2e/dv dv * r^4/(1-d)^2 (doublet block):\n", Hd)
                print("CLAIMED table read as symmetric Hessian entries:\n", Hc)
                print("linear terms * r^4/(1-d)^2, all 24 jets:\n", lin_scaled)
                print("MIXED block doublet (rows) x director (cols):", [NAMES[k] for k in direc])
                print(Hs[np.ix_(doublet, direc)])
                print("MIXED block doublet (rows) x massive (cols):", [NAMES[k] for k in massive])
                print(Hs[np.ix_(doublet, massive)])
                print("director x director block:\n", Hs[np.ix_(direc, direc)])
        print()
        print("over 8 random points:")
        print("  max |linear term| on doublet jets (scaled)        : %.3e" % worst_lin_doublet)
        print("  max rel deviation of doublet Hessian from claim   : %.3e" % worst_claim)
        print("  max rel deviation of doublet Hessian from 2 x claim: %.3e" % worst_claim2)
        print("  FD self-check of my (lin, H) vs the quartic density: %.3e" % worstfd)

        # analytic-form fit of the mixed block: print at a clean point to read off closed forms
        for TH in (np.pi / 3, np.pi / 2, 1.0):
            R, PH = 1.0, 0.7
            A0, G = m.jets_at(R, TH, PH)
            e0, lin, H = m.variations(A0, G)
            Hs = H * R**4 / (1 - dl) ** 2
            print(
                "\ntheta = %.4f  (sin %.4f cos %.4f, 1/sin %.4f, cot %.4f)"
                % (TH, np.sin(TH), np.cos(TH), 1 / np.sin(TH), 1 / np.tan(TH))
            )
            print("doublet x director:\n", Hs[np.ix_(doublet, direc)])
            print("doublet x massive:\n", Hs[np.ix_(doublet, massive)])
            lin_scaled = {
                n: float(v) for n, v in zip(NAMES, np.round(lin * R**4 / (1 - dl) ** 2, 6))
            }
            print("linear terms (scaled):", lin_scaled)

    record(
        "C6a.no_linear_term_on_doublet",
        "the first variation vanishes on all 8 doublet jets at 8 random points (scaled, < 1e-10)",
        worst_lin_doublet < 1e-10,
        {"max_abs_linear_term_scaled": worst_lin_doublet, "threshold": 1e-10},
    )
    record(
        "C6.fd_self_check",
        "own (linear, Hessian) pair matches five-point stencils of the full quartic density "
        "(rel < 1e-8)",
        worstfd < 1e-8,
        {"max_rel_error": worstfd, "threshold": 1e-8},
    )
    record(
        "C6b.table_as_quadratic_form",
        "the listed table t gives the second variation as e2 = v^T t v, that is "
        "d2e / dv dv = 2 t, on the doublet block (rel < 1e-10)",
        worst_claim2 < 1e-10,
        {"max_rel_deviation_from_2t": worst_claim2, "threshold": 1e-10},
    )
    record(
        "C6b.table_as_hessian_entries",
        "the listed table equals the Hessian entries d2e / dv dv as labeled (rel < 1e-10)",
        worst_claim < 1e-10,
        {"max_rel_deviation": worst_claim, "threshold": 1e-10},
        note="a uniform factor 2 on every entry, a labeling convention: the table holds the "
        "coefficients of the quadratic form, not second derivatives; the angular eigenvalue "
        "is a ratio inside the table and does not move",
    )
    record(
        "C6.mixed_blocks_nonzero_pointwise",
        "the doublet x director and the doublet x massive blocks are nonzero pointwise, so "
        "any decoupling has to come from the angular integration or the constraint",
        first["max_abs_doublet_x_director"] > 1 and first["max_abs_doublet_x_massive"] > 1,
        first,
    )


# ================= C6 Galerkin on the sphere =================
def galerkin_machinery(LMAX):
    """Angular eigenproblem by a CARTESIAN TENSOR GALERKIN method on the sphere.
    No theta grid, no spin-weighted-harmonic formula: basis = symmetric 3x3 tensors with
    polynomial entries in xhat (degree <= LMAX), projected pointwise onto a sector:
       doublet  : P S P - 1/2 tr(P S P) P          (P = I - xhat xhat^T)
       director : P S xx^T + xx^T S P
       rr       : (x^T S x) xx^T ;  tt : tr(PSP) P
    Energy density e = 4 sum_{i<j} |[d_i M, d_j M]|^2 depends on gradients only. For
    psi = r^p Psi(xhat):
       d_k psi |_{r=1} = p xhat_k Psi + grad_tan_k Psi .
    Second-variation bilinear form B (e2 = B(psi,psi)), in units (1-delta)^2 / r^4:
       R = B(radial jet, radial jet), C = B(radial, tangential), K = B(tangential, tangential).
    Power-law solutions satisfy [K + C - p(p-1) R + p (C^T - C)] v = 0; if C = 0:
    K v = A R v, A = p(p-1), i.e. the radial equation eps'' = A eps / r^2.
    Constrained (eigenvalues frozen) director block: + L(psi_2) = -48 tau_a.tau_b (units
    above)."""
    if ("galerkin", LMAX) in _CACHE:
        return _CACHE[("galerkin", LMAX)]
    NT, NP = 36, 72

    # quadrature
    ct, wt = np.polynomial.legendre.leggauss(NT)
    phi = (np.arange(NP) + 0.5) * 2 * np.pi / NP
    CT, PH = np.meshgrid(ct, phi, indexing="ij")
    ST = np.sqrt(1 - CT**2)
    XH = np.stack([ST * np.cos(PH), ST * np.sin(PH), CT], -1).reshape(-1, 3)
    WQ = (wt[:, None] * np.ones(NP)[None, :] * 2 * np.pi / NP).ravel()
    NPT = len(WQ)
    assert abs(WQ.sum() - 4 * np.pi) < 1e-12

    # monomials and free gradients
    exps = [e for e in itertools.product(range(LMAX + 1), repeat=3) if sum(e) <= LMAX]

    def mono(e):
        return XH[:, 0] ** e[0] * XH[:, 1] ** e[1] * XH[:, 2] ** e[2]

    def dmono(e):
        g = np.zeros((NPT, 3))
        for m in range(3):
            if e[m] > 0:
                ee = list(e)
                ee[m] -= 1
                g[:, m] = e[m] * mono(ee)
        return g

    units = []
    for i in range(3):
        for j in range(i, 3):
            E = np.zeros((3, 3))
            E[i, j] = E[j, i] = 1
            units.append(E)
    S_val = (
        np.array([[mono(e)[:, None, None] * U for U in units] for e in exps])
        .reshape(-1, NPT, 3, 3)
        .transpose(1, 0, 2, 3)
    )
    S_grd = (
        np.array([[dmono(e)[:, :, None, None] * U for U in units] for e in exps])
        .reshape(-1, NPT, 3, 3, 3)
        .transpose(1, 0, 2, 3, 4)
    )
    NB = S_val.shape[1]
    I3 = np.eye(3)
    xx = np.einsum("pi,pj->pij", XH, XH)
    Pv = I3[None] - xx  # (p,3,3)
    dxx = np.einsum("mi,pj->pmij", I3, XH) + np.einsum("pi,mj->pmij", XH, I3)
    Pg = -dxx  # (p,m,3,3)

    class TF:  # tensor field with free gradient; shapes (p, nb, 3,3) and (p, nb, m, 3,3)
        def __init__(s, v, g):
            s.v, s.g = v, g

        def __matmul__(s, o):
            return TF(s.v @ o.v, s.g @ o.v[:, :, None] + s.v[:, :, None] @ o.g)

        def __add__(s, o):
            return TF(s.v + o.v, s.g + o.g)

        def __sub__(s, o):
            return TF(s.v - o.v, s.g - o.g)

        def scal(s, f, df):  # multiply by scalar field f (p,nb), df (p,nb,m)
            return TF(
                f[..., None, None] * s.v,
                df[..., None, None] * s.v[:, :, None] + f[:, :, None, None, None] * s.g,
            )

        def tr(s):
            return np.einsum("pbii->pb", s.v), np.einsum("pbmii->pbm", s.g)

    def bc(v, g):
        return TF(
            np.broadcast_to(v[:, None], (NPT, NB, 3, 3)),
            np.broadcast_to(g[:, None], (NPT, NB, 3, 3, 3)),
        )

    S = TF(S_val, S_grd)
    P = bc(Pv, Pg)
    XXf = bc(xx, dxx)
    PSP = P @ S @ P
    t, dt = PSP.tr()
    sectors = {
        "doublet": PSP - P.scal(0.5 * t, 0.5 * dt),
        "director": (P @ S @ XXf) + (XXf @ S @ P),
    }
    xSx, dxSx = (XXf @ S @ XXf).tr()
    sectors["rr"] = XXf.scal(xSx, dxSx)
    sectors["tt"] = P.scal(t, dt)

    # background
    A0 = (1 - DELTA) * np.einsum("pkm,pmij->pkij", Pv, dxx)  # d_k M0 at r=1

    def comm(A, B):
        return A @ B - B @ A

    F0 = [comm(A0[:, i], A0[:, j]) for i, j in PAIRS]
    e0 = 4 * sum(np.sum(F * F, axis=(1, 2)) for F in F0) / (1 - DELTA) ** 2
    print("check e0 r^4/(1-d)^2 =", e0[:3], "(pointwise check: 3.92)")

    def jets(tf):
        Xt = np.einsum("pkm,pbmij->pbkij", Pv, tf.g)  # tangential gradient
        Xr = np.einsum("pk,pbij->pbkij", XH, tf.v)  # radial jet (times p)
        return Xr, Xt

    def bilinear(Xa, Xb):
        """B_ab = int dOmega 4 sum_q [<F1a,F1b> + <F0, [Xa_i,Xb_j] + [Xb_i,Xa_j]>],
        / (1-delta)^2."""
        na, nb = Xa.shape[1], Xb.shape[1]
        out = np.zeros((na, nb))
        sw = np.sqrt(WQ)[:, None, None, None]
        for q, (i, j) in enumerate(PAIRS):
            Ai, Aj = A0[:, None, i], A0[:, None, j]
            F1a = comm(Ai, Xa[:, :, j]) + comm(Xa[:, :, i], Aj)
            F1b = comm(Ai, Xb[:, :, j]) + comm(Xb[:, :, i], Aj)
            out += 4 * np.tensordot(
                (sw * F1a).transpose(1, 0, 2, 3).reshape(na, -1),
                (sw * F1b).transpose(1, 0, 2, 3).reshape(nb, -1).T,
                1,
            )
            G = F0[q][:, None]  # (p,1,3,3)
            GT = np.swapaxes(G, -1, -2)

            # <G, Xa_i Xb_j - Xb_j Xa_i + Xb_i Xa_j - Xa_j Xb_i> ;
            # tr(G^T U V) = sum (G^T U)_{ln} V_{nl}
            def tr3(U, V, GT=GT):
                L_ = (sw**2 * (GT @ U)).transpose(1, 0, 2, 3).reshape(U.shape[1], -1)
                R_ = np.swapaxes(V, -1, -2).transpose(1, 0, 2, 3).reshape(V.shape[1], -1)
                return L_ @ R_.T  # [u, v]

            out += 4 * (
                tr3(Xa[:, :, i], Xb[:, :, j])
                - tr3(Xb[:, :, j], Xa[:, :, i]).T
                + tr3(Xb[:, :, i], Xa[:, :, j]).T
                - tr3(Xa[:, :, j], Xb[:, :, i])
            )
        return out / (1 - DELTA) ** 2

    def gram(tfa, tfb):
        return np.einsum("p,paij,pbij->ab", WQ, tfa.v, tfb.v, optimize=True)

    def reduce_basis(G, tol=1e-9):
        w, V = np.linalg.eigh(G)
        keep = w > tol * w.max()
        return V[:, keep] / np.sqrt(w[keep])

    def solve(names, constrained=True, label=""):
        Xr = np.concatenate([jets(sectors[n])[0] for n in names], 1)
        Xt = np.concatenate([jets(sectors[n])[1] for n in names], 1)
        allv = TF(np.concatenate([sectors[n].v for n in names], 1), None)
        G = np.einsum("p,paij,pbij->ab", WQ, allv.v, allv.v, optimize=True)
        K = bilinear(Xt, Xt)
        R = bilinear(Xr, Xr)
        C = bilinear(Xr, Xt)
        if constrained and "director" in names:
            o = sum(NB for n in names[: names.index("director")])
            d = sectors["director"].v
            # tau_a . tau_b = <psi_a, psi_b>_F / 2 for psi = r tau^T + tau r^T
            K[o : o + NB, o : o + NB] += (
                -48.0 * 0.5 * np.einsum("p,paij,pbij->ab", WQ, d, d, optimize=True)
            )
        T = reduce_basis(G)
        Kr, Rr, Cr = T.T @ K @ T, T.T @ R @ T, T.T @ C @ T
        print("\n=== %s : sectors %s, LMAX %d, rank %d ===" % (label, names, LMAX, T.shape[1]))
        print(
            "  |K - K^T| %.2e  |C| %.2e (vs |K| %.2e)  |R - 8*Gram(=I)| %.2e"
            % (
                np.abs(Kr - Kr.T).max(),
                np.abs(Cr).max(),
                np.abs(Kr).max(),
                np.abs(Rr - 8 * np.eye(len(Rr))).max(),
            )
        )
        wR = np.linalg.eigvalsh(0.5 * (Rr + Rr.T))
        print("  R eigen range: %.4f .. %.4f" % (wR.min(), wR.max()))
        Ks = 0.5 * (Kr + Kr.T) + 0.5 * (Cr + Cr.T)
        w, V = eigh(Ks, 0.5 * (Rr + Rr.T))
        lam = w * 16  # angular eigenvalue in the claimant's units (16 r^2 radial coefficient)
        # group multiplicities
        groups = []
        for val in lam:
            if groups and abs(val - groups[-1][0]) < 1e-6 * max(1, abs(val)):
                groups[-1][1] += 1
            else:
                groups.append([val, 1])
        print("  lowest angular eigenvalues (claimant units; A = value/16) with multiplicity:")
        for val, mult in groups[:14]:
            print("     lambda = %+10.5f   A = %+8.5f   mult %d" % (val, val / 16, mult))
        info = {
            "rank": T.shape[1],
            "K_asymmetry": np.abs(Kr - Kr.T).max(),
            "C_max": np.abs(Cr).max(),
            "K_max": np.abs(Kr).max(),
            "R_eigen_range": [wR.min(), wR.max()],
            "groups": [[float(v), int(mu)] for v, mu in groups[:14]],
        }
        return lam, V, T, names, info

    m = SimpleNamespace(
        LMAX=LMAX, NB=NB, WQ=WQ, XH=XH, A0=A0, F0=F0, sectors=sectors, comm=comm, jets=jets,
        bilinear=bilinear, gram=gram, reduce_basis=reduce_basis, solve=solve,
    )  # fmt: skip
    _CACHE[("galerkin", LMAX)] = m
    return m


def check_c6_galerkin():
    """C6 (c) and (iii) at LMAX 2 and LMAX 3."""
    for LMAX in (2, 3):
        print("\n=== C6 Galerkin on the sphere, LMAX %d ===" % LMAX)
        m = galerkin_machinery(LMAX)
        sectors, NB, WQ = m.sectors, m.NB, m.WQ
        tag = "_L%d" % LMAX
        lam_d, _, _, _, info_d = m.solve(["doublet"], label="(c) doublet alone")
        _, _, _, _, info_dir = m.solve(
            ["director"],
            constrained=True,
            label="director alone, CONSTRAINED (n unit, eigenvalues frozen)",
        )
        _, _, _, _, info_dir_lin = m.solve(
            ["director"],
            constrained=False,
            label="director alone, linear parameterization (no -48 term)",
        )
        lam_c, V, T, names, info_c = m.solve(
            ["doublet", "director"],
            constrained=True,
            label="(iii) COUPLED doublet + director, constrained",
        )
        # doublet weight of each coupled eigenvector
        allv = np.concatenate([sectors[n].v for n in names], 1)
        Gfull = np.einsum("p,paij,pbij->ab", WQ, allv, allv, optimize=True)
        Pd = np.zeros_like(Gfull)
        Pd[:NB, :NB] = Gfull[:NB, :NB]
        coef = T @ V
        wd = np.einsum("ai,ab,bi->i", coef, Pd, coef) / np.einsum("ai,ab,bi->i", coef, Gfull, coef)
        print("  doublet weight of the lowest coupled modes:")
        for k in range(min(40, len(lam_c))):
            print(
                "     lambda %+10.5f  A %+8.5f  doublet weight %.4f"
                % (lam_c[k], lam_c[k] / 16, wd[k])
            )
        _, _, _, _, info_c_lin = m.solve(
            ["doublet", "director"],
            constrained=False,
            label="COUPLED doublet + director, LINEAR parameterization (for reference)",
        )

        print(
            "\n=== integrated MIXED blocks between sectors "
            "(each sector in its own Frobenius-orthonormal basis) ==="
        )
        red = {}
        for n in sectors:
            red[n] = m.reduce_basis(m.gram(sectors[n], sectors[n]))
        J = {n: m.jets(sectors[n]) for n in sectors}
        names4 = list(sectors)
        mixed = {}
        for i, na in enumerate(names4):
            for nb_ in names4[i + 1 :]:
                Kab = red[na].T @ m.bilinear(J[na][1], J[nb_][1]) @ red[nb_]
                Rab = red[na].T @ m.bilinear(J[na][0], J[nb_][0]) @ red[nb_]
                Cab = red[na].T @ m.bilinear(J[na][0], J[nb_][1]) @ red[nb_]
                Cba = red[na].T @ m.bilinear(J[na][1], J[nb_][0]) @ red[nb_]
                sv = np.linalg.svd(Kab, compute_uv=False)
                mixed["%s x %s" % (na, nb_)] = {
                    "K_max": np.abs(Kab).max(),
                    "K_largest_singular_value": sv[0],
                    "R_max": np.abs(Rab).max(),
                    "C_rt_max": np.abs(Cab).max(),
                    "C_tr_max": np.abs(Cba).max(),
                }
                print(
                    "  %-8s x %-8s : max|K| %.3e (largest singular value %.4f)  max|R| %.3e  "
                    "max|C_rt| %.3e  max|C_tr| %.3e"
                    % (
                        na,
                        nb_,
                        np.abs(Kab).max(),
                        sv[0],
                        np.abs(Rab).max(),
                        np.abs(Cab).max(),
                        np.abs(Cba).max(),
                    )
                )
        # pointwise (non-integrated) size of the doublet x director mixed density, to show it
        # is not zero before integration
        Xa, Xb = J["doublet"][1][:, :6], J["director"][1][:, :6]
        dens = 0
        for q, (i, j) in enumerate(PAIRS):
            Ai, Aj = m.A0[:, None, i], m.A0[:, None, j]
            F1a = m.comm(Ai, Xa[:, :, j]) + m.comm(Xa[:, :, i], Aj)
            F1b = m.comm(Ai, Xb[:, :, j]) + m.comm(Xb[:, :, i], Aj)
            dens = dens + 4 * np.einsum("paij,pbij->pab", F1a, F1b)
            F0T = np.swapaxes(m.F0[q], -1, -2)

            def U(A_, B_, F0T=F0T):
                return np.einsum("pij,pajk,pbki->pab", F0T, A_, B_)

            dens = dens + 4 * (
                U(Xa[:, :, i], Xb[:, :, j])
                - np.swapaxes(U(Xb[:, :, j], Xa[:, :, i]), 1, 2)
                + np.swapaxes(U(Xb[:, :, i], Xa[:, :, j]), 1, 2)
                - U(Xa[:, :, j], Xb[:, :, i])
            )
        dens_max = np.abs(dens).max() / (1 - DELTA) ** 2
        dens_int = np.abs(np.einsum("p,pab->ab", WQ, dens)).max() / (1 - DELTA) ** 2
        print(
            "  pointwise doublet x director mixed density (first 6x6 basis pairs): "
            "max |density| = %.4f ; max |integral| = %.3e" % (dens_max, dens_int)
        )

        gd = info_d["groups"]
        record(
            "C6c.doublet_lowest_48" + tag,
            "the lowest angular eigenvalue of the doublet sector is 48 (A = 3), multiplicity "
            "10 = 2 x 5",
            abs(lam_d[0] - 48) < 1e-6 and gd[0][1] == 10,
            {"lowest": lam_d[0], "multiplicity": gd[0][1], "rank": info_d["rank"]},
        )
        record(
            "C6c.doublet_ladder" + tag,
            "the doublet ladder starts 48, 96, 160 = 8 (l (l + 1) - 4) + 32 for l = 2, 3, 4",
            all(abs(gd[k][0] - v) < 1e-6 for k, v in enumerate((48.0, 96.0, 160.0))),
            {"groups_value_multiplicity": gd},
        )
        record(
            "C6.radial_tangential_block_vanishes" + tag,
            "the radial x tangential block C vanishes on the doublet (|C| / |K| < 1e-10), so "
            "power laws separate and A = p (p - 1)",
            info_d["C_max"] / info_d["K_max"] < 1e-10,
            {"C_max": info_d["C_max"], "K_max": info_d["K_max"]},
        )
        in48 = [k for k in range(len(lam_c)) if abs(lam_c[k] - 48) < 1e-6]
        below = [[float(v), int(mu)] for v, mu in info_c["groups"] if v < 48 - 1e-6]
        record(
            "C6iii.coupled_48_survives" + tag,
            "with the constrained director coupled in, the eigenvalue 48 keeps multiplicity 10 "
            "and doublet weight 1 (> 0.999)",
            len(in48) == 10 and min(wd[in48]) > 0.999,
            {
                "multiplicity_at_48": len(in48),
                "min_doublet_weight_at_48": min(wd[in48]) if in48 else None,
                "groups_below_48_value_multiplicity": below,
                "max_doublet_weight_below_48": max(
                    [wd[k] for k in range(len(lam_c)) if lam_c[k] < 48 - 1e-6], default=None
                ),
                "coupled_groups": info_c["groups"],
                "coupled_groups_linear_parameterization": info_c_lin["groups"],
                "director_alone_constrained": info_dir["groups"],
                "director_alone_linear": info_dir_lin["groups"],
            },
            note="the modes below 48 (A = 0 and A = 2) are pure director modes, doublet weight 0",
        )
        dd = mixed["doublet x director"]
        record(
            "C6iii.doublet_director_block_integrates_to_zero" + tag,
            "the doublet x director block is nonzero pointwise and integrates to zero over "
            "the sphere (< 1e-9)",
            dd["K_max"] < 1e-9 and dens_int < 1e-9 and dens_max > 1,
            {
                "integrated": dd,
                "pointwise_density_max": dens_max,
                "pointwise_density_integral_max": dens_int,
            },
        )
        record(
            "C6.doublet_massive_block" + tag,
            "the doublet decouples from the massive rr and tt sectors in the curvature term "
            "(integrated K block < 1e-9)",
            mixed["doublet x rr"]["K_max"] < 1e-9 and mixed["doublet x tt"]["K_max"] < 1e-9,
            {k: v for k, v in mixed.items() if k != "doublet x director"},
            kind="scope",
            note="the integrated doublet x massive block is of the size of K itself; the "
            "massive modes carry a potential mass, and the Schur check sizes the effect",
        )
        if LMAX != 2:
            _CACHE.pop(("galerkin", LMAX), None)


# ================= C6 massive Schur complement =================
def check_c6_massive_schur():
    """The doublet x massive block is nonzero. Size of its effect: adiabatic Schur complement
    of the massive (rr, tt) modes, with the V4 - cL Hessian as their mass (M_00 slaved), at
    several r. Reuses the Galerkin machinery at LMAX 2."""
    print("\n=== C6 massive Schur complement (LMAX 2) ===")
    m = galerkin_machinery(2)
    sectors, jets, bilinear, reduce_basis = m.sectors, m.jets, m.bilinear, m.reduce_basis
    WQ, XH, NB = m.WQ, m.XH, m.NB
    W = 25.0 * 0.000724023879
    C_LIN = 1e-3

    # potential Hessian in (m00, s, u): N eigenvalues (-m00, 1+s, delta+u, delta+u)
    roots = [-G8, 1.0, DELTA, DELTA]
    poly = np.poly(sorted(set(roots)))[::-1]
    a_p = [poly[p - 1] / p for p in range(1, 5)]
    Cp = [sum(q**p for q in roots) for p in range(1, 5)]

    def V(v):
        m_, s, u = v
        ev = [-(G8 + m_), 1 + s, DELTA + u, DELTA + u]
        t = [sum(x**p for x in ev) - Cp[p - 1] for p in range(1, 5)]
        return W * sum(x * x for x in t) - C_LIN * sum(a_p[k] * t[k] for k in range(4))

    hh = 1e-4
    HV = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            ei, ej = np.eye(3)[i] * hh, np.eye(3)[j] * hh
            HV[i, j] = (V(ei + ej) - V(ei - ej) - V(-ei + ej) + V(-ei - ej)) / (4 * hh * hh)
    print("potential Hessian (m00, s, u):\n", HV)
    Hsu = HV[1:, 1:] - np.outer(HV[1:, 0], HV[0, 1:]) / HV[0, 0]
    print("after slaving M_00, (s,u) Hessian:\n", Hsu, " eig:", np.linalg.eigvalsh(Hsu))

    names_x = ["doublet", "director"]
    names_m = ["rr", "tt"]
    J = {n: jets(sectors[n]) for n in sectors}

    def cat(names, k):
        return np.concatenate([J[n][k] for n in names], 1)

    vx = np.concatenate([sectors[n].v for n in names_x], 1)
    vm = np.concatenate([sectors[n].v for n in names_m], 1)
    Tx = reduce_basis(np.einsum("p,paij,pbij->ab", WQ, vx, vx, optimize=True))
    Tm = reduce_basis(np.einsum("p,paij,pbij->ab", WQ, vm, vm, optimize=True))
    Kxx = bilinear(cat(names_x, 1), cat(names_x, 1))
    d = sectors["director"].v
    # constrained director block
    Kxx[NB:, NB:] += -24.0 * np.einsum("p,paij,pbij->ab", WQ, d, d, optimize=True)
    Rxx = bilinear(cat(names_x, 0), cat(names_x, 0))
    Kxm = bilinear(cat(names_x, 1), cat(names_m, 1))
    Kmm = bilinear(cat(names_m, 1), cat(names_m, 1))
    xx = np.einsum("pi,pj->pij", XH, XH)
    s_co = np.einsum("paij,pij->pa", vm, xx)  # coefficient of rhat rhat^T
    u_co = 0.5 * (np.einsum("paii->pa", vm) - s_co)  # coefficient of P
    Pot = 0.5 * (
        Hsu[0, 0] * np.einsum("p,pa,pb->ab", WQ, s_co, s_co)
        + Hsu[1, 1] * np.einsum("p,pa,pb->ab", WQ, u_co, u_co)
        + Hsu[0, 1]
        * (np.einsum("p,pa,pb->ab", WQ, s_co, u_co) + np.einsum("p,pa,pb->ab", WQ, u_co, s_co))
    )
    Kxx, Rxx, Kxm, Kmm, Pot = (
        Tx.T @ Kxx @ Tx,
        Tx.T @ Rxx @ Tx,
        Tx.T @ Kxm @ Tm,
        Tm.T @ Kmm @ Tm,
        Tm.T @ Pot @ Tm,
    )

    w0, V0 = eigh(0.5 * (Kxx + Kxx.T), Rxx)
    print("\nno massive coupling: lowest A:", np.round(w0[:30], 4))
    # doublet projector for weights
    Gx = Tx.T @ np.einsum("p,paij,pbij->ab", WQ, vx, vx, optimize=True) @ Tx
    Gd = np.zeros((2 * NB, 2 * NB))
    Gd[:NB, :NB] = np.einsum("p,paij,pbij->ab", WQ, vx[:, :NB], vx[:, :NB], optimize=True)
    Gd = Tx.T @ Gd @ Tx
    print(
        "\n   r     doublet-dominated eigenvalues A_eff (adiabatic, massive modes eliminated "
        "at that r); weight>0.5 and A<4"
    )
    rows = {}
    for r in (4.0, 6.0, 9.0, 12.0, 18.0, 1e3):
        mass = Pot * r**4 / (1 - DELTA) ** 2
        Keff = Kxx - Kxm @ np.linalg.solve(Kmm + mass, Kxm.T)
        w, Vv = eigh(0.5 * (Keff + Keff.T), Rxx)
        wd = np.einsum("ai,ab,bi->i", Vv, Gd, Vv) / np.einsum("ai,ab,bi->i", Vv, Gx, Vv)
        sel = [
            (round(float(a), 4), round(float(x), 3)) for a, x in zip(w, wd) if x > 0.5 and a < 4.5
        ]
        rows[r] = {
            "min_eig_Kmm_plus_mass": np.linalg.eigvalsh(Kmm + mass).min(),
            "overall_lowest_A": w[0],
            "lowest_doublet_like_A": min(a for a, x in sel),
            "its_doublet_weight": min(sel)[1],
            "doublet_like_A_weight": sel,
        }
        print(
            "  %6.1f  min eig(Kmm+mass) %.2f  overall lowest A %.4f ; doublet-like "
            "(A, doublet weight): %s" % (r, rows[r]["min_eig_Kmm_plus_mass"], w[0], sel)
        )

    def dev(r):
        return abs(rows[r]["lowest_doublet_like_A"] - 3) / 3

    numbers = {
        "potential_hessian_m00_s_u": HV,
        "s_u_hessian_after_slaving_M00": Hsu,
        "its_eigenvalues": np.linalg.eigvalsh(Hsu),
        "lowest_A_without_massive_coupling": np.round(w0[:30], 4),
        "rows_by_r": rows,
    }
    record(
        "C6m.A3_asymptotic",
        "with the massive modes eliminated adiabatically, the lowest doublet-like A returns "
        "to 3 at large r (r = 1000, rel < 1e-3)",
        dev(1e3) < 1e-3,
        numbers,
    )
    record(
        "C6m.A3_outer_window",
        "A_eff stays within 5 percent of 3 at r = 12 and r = 18",
        dev(12.0) < 0.05 and dev(18.0) < 0.05,
        {"rel_deviation": {"r=12": dev(12.0), "r=18": dev(18.0)}, "threshold": 0.05},
    )
    record(
        "C6m.A3_inner_window",
        "A_eff stays within 5 percent of 3 at r = 6 and r = 9",
        dev(6.0) < 0.05 and dev(9.0) < 0.05,
        {
            "A_eff": {
                "r=6": rows[6.0]["lowest_doublet_like_A"],
                "r=9": rows[9.0]["lowest_doublet_like_A"],
            },
            "rel_deviation": {"r=6": dev(6.0), "r=9": dev(9.0)},
            "overall_lowest_A_at_r4": rows[4.0]["overall_lowest_A"],
            "threshold": 0.05,
        },
        kind="scope",
        note="an adiabatic estimate, not a solved radial problem: the coupling to the massive "
        "modes pulls A_eff below 3 inside r of about 10, where the mode is also only partly "
        "doublet (weight 0.55 at r = 6)",
    )


# ================= C6 constrained finite differences =================
def check_c6_constrained_fd():
    """C6 (i) second route + check of the constrained director block: brute-force numerical
    second variation of the curvature density on the EXACT nonlinear family
       M = delta I + (1-delta) n n^T + a (e1 e1^T - e2 e2^T) + b (e1 e2^T + e2 e1^T),
       n = normalize(rhat + (d1 that + d2 phat)/(1-delta)),  e1,e2 = Gram-Schmidt(that, phat | n),
    with a, b, d1, d2 affine in (r, theta, phi) around the point. Spatial derivatives by
    6th-order central differences in CARTESIAN x,y,z; jet derivatives by 5-point stencils.
    Nothing shared with the pointwise check except the final comparison."""
    print("\n=== C6 constrained finite differences on the exact nonlinear family ===")
    rng = np.random.default_rng(11)
    a4 = pointwise_machinery()

    def sph(X):
        r = np.linalg.norm(X)
        th = np.arccos(X[2] / r)
        ph = np.arctan2(X[1], X[0])
        return r, th, ph

    def frame(th, ph):
        rh = np.array([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)])
        tt = np.array([np.cos(th) * np.cos(ph), np.cos(th) * np.sin(ph), -np.sin(th)])
        pp = np.array([-np.sin(ph), np.cos(ph), 0.0])
        return rh, tt, pp

    def M_of(X, v, P0):
        r, th, ph = sph(X)
        dphi = (ph - P0[2] + np.pi) % (2 * np.pi) - np.pi
        basis = np.array([1.0, r - P0[0], th - P0[1], dphi])
        a, b, d1, d2 = [v[4 * k : 4 * k + 4] @ basis for k in range(4)]
        rh, tt, pp = frame(th, ph)
        n = rh + (d1 * tt + d2 * pp) / (1 - DELTA)
        n /= np.linalg.norm(n)
        e1 = tt - (tt @ n) * n
        e1 /= np.linalg.norm(e1)
        e2 = pp - (pp @ n) * n - (pp @ e1) * e1
        e2 /= np.linalg.norm(e2)
        return (
            DELTA * np.eye(3)
            + (1 - DELTA) * np.outer(n, n)
            + a * (np.outer(e1, e1) - np.outer(e2, e2))
            + b * (np.outer(e1, e2) + np.outer(e2, e1))
        )

    def density(v, X0, P0, h):
        D = []
        for k in range(3):
            e = np.zeros(3)
            e[k] = h
            D.append(
                (
                    45 * (M_of(X0 + e, v, P0) - M_of(X0 - e, v, P0))
                    - 9 * (M_of(X0 + 2 * e, v, P0) - M_of(X0 - 2 * e, v, P0))
                    + (M_of(X0 + 3 * e, v, P0) - M_of(X0 - 3 * e, v, P0))
                )
                / (60 * h)
            )
        tot = 0.0
        for i in range(3):
            for j in range(i + 1, 3):
                F = D[i] @ D[j] - D[j] @ D[i]
                tot += 4 * np.sum(F * F)
        return tot

    def hessian(X0, P0, s=2e-2):
        h = 2e-2 * P0[0]

        def f(v):
            return density(v, X0, P0, h)

        n = 16
        H = np.zeros((n, n))
        g = np.zeros(n)
        f0 = f(np.zeros(n))
        E = np.eye(n) * s
        for i in range(n):
            fm2, fm1, fp1, fp2 = f(-2 * E[i]), f(-E[i]), f(E[i]), f(2 * E[i])
            g[i] = (fm2 - 8 * fm1 + 8 * fp1 - fp2) / (12 * s)
            H[i, i] = (-fm2 + 16 * fm1 - 30 * f0 + 16 * fp1 - fp2) / (12 * s * s)
        for i in range(n):
            for j in range(i + 1, n):
                # 4th-order mixed derivative
                def c(k, i=i, j=j):
                    return (
                        f(k * E[i] + k * E[j])
                        - f(k * E[i] - k * E[j])
                        - f(-k * E[i] + k * E[j])
                        + f(-k * E[i] - k * E[j])
                    )

                H[i, j] = H[j, i] = (16 * c(1) - c(2)) / (48 * s * s)
        return g, H

    points = []
    with np.printoptions(linewidth=220, precision=3, suppress=True):
        for trial in range(2):
            P0 = np.array([rng.uniform(1.0, 5.0), rng.uniform(0.5, 2.6), rng.uniform(0.3, 5.9)])
            X0 = P0[0] * frame(P0[1], P0[2])[0]
            g, H = hessian(X0, P0)
            sc = P0[0] ** 4 / (1 - DELTA) ** 2
            H *= sc
            g *= sc
            A0, G = a4.jets_at(*P0)
            e0, lin, Hl = a4.variations(A0, G)
            Hl = Hl[:16, :16] * sc
            names = a4.NAMES[:16]
            sn, cs = np.sin(P0[1]), np.cos(P0[1])
            print("\npoint r=%.3f th=%.3f ph=%.3f" % tuple(P0))
            d_pointwise = np.abs(H[:8, :8] - Hl[:8, :8]).max()
            print("  doublet block: max |H_bruteforce - H_pointwise| = %.2e" % d_pointwise)
            cl = a4.claimed_H(P0[0], P0[1])
            Hc = np.zeros((8, 8))
            for (p, q), val in cl.items():
                i, j = names.index(p), names.index(q)
                Hc[i, j] = Hc[j, i] = val
            d_table = np.abs(H[:8, :8] - 2 * Hc).max()
            print(
                "  doublet block: max |H_bruteforce - 2 * claimed| = %.2e  (scale %.1f)"
                % (d_table, np.abs(Hc).max())
            )
            Ddd = H[8:, 8:] - Hl[8:, 8:]
            print("  director block: H_constrained - H_linear =\n", Ddd)
            pred = np.zeros((8, 8))
            pred[0, 0] = pred[4, 4] = -96.0
            d_dir = np.abs(Ddd - pred).max()
            print("  predicted -96 on (d1,d1),(d2,d2): max dev %.2e" % d_dir)
            Dmix = H[:8, 8:] - Hl[:8, 8:]
            # predicted mixed correction: L_dir applied to tau' = -(psi tau)/(1-delta):
            # tau'_1 = -(a d1 + b d2)/(1-d), tau'_2 = -(b d1 - a d2)/(1-d);
            # L = 11.2 [ d1'_t + cot d1' + d2'_p / sin ]
            pm = np.zeros((8, 8))
            k = -11.2 / (1 - DELTA)
            ia, ia_t, ia_p, ib, ib_t, ib_p = 0, 2, 3, 4, 6, 7
            jd1, jd1_t, jd1_p, jd2, jd2_t, jd2_p = 0, 2, 3, 4, 6, 7
            # d/dt (a d1 + b d2) + cot (a d1 + b d2) + (1/sin) d/dp (b d1 - a d2)
            for i, j, val in (
                (ia_t, jd1, 1),
                (ia, jd1_t, 1),
                (ib_t, jd2, 1),
                (ib, jd2_t, 1),
                (ia, jd1, cs / sn),
                (ib, jd2, cs / sn),
                (ib_p, jd1, 1 / sn),
                (ib, jd1_p, 1 / sn),
                (ia_p, jd2, -1 / sn),
                (ia, jd2_p, -1 / sn),
            ):
                pm[i, j] += k * val
            d_mix = np.abs(Dmix - pm).max()
            print(
                "  mixed block: H_constrained - H_linear, max dev from the divergence-form "
                "prediction: %.2e (size of correction %.2f)" % (d_mix, np.abs(Dmix).max())
            )
            print("  constrained MIXED block doublet x director (pointwise, scaled):\n", H[:8, 8:])
            print("  linear terms brute force (scaled):", np.round(g, 4))
            points.append(
                {
                    "point_r_theta_phi": P0,
                    "doublet_block_vs_pointwise": d_pointwise,
                    "doublet_block_vs_2x_table": d_table,
                    "table_scale": np.abs(Hc).max(),
                    "director_block_dev_from_minus_96": d_dir,
                    "mixed_block_dev_from_divergence_form": d_mix,
                    "mixed_correction_size": np.abs(Dmix).max(),
                    "constrained_mixed_block_max_abs": np.abs(H[:8, 8:]).max(),
                    "max_abs_linear_term_on_doublet_scaled": np.abs(g[:8]).max(),
                    "linear_terms_scaled": np.round(g, 4),
                }
            )

    def worst(key):
        return max(pt[key] for pt in points)

    record(
        "C6fd.doublet_block_vs_pointwise",
        "a brute-force second variation on the exact nonlinear family reproduces the "
        "pointwise doublet Hessian (abs < 1e-4 on a scale of 40)",
        worst("doublet_block_vs_pointwise") < 1e-4,
        {"points": points, "threshold": 1e-4},
    )
    record(
        "C6fd.doublet_block_vs_2x_table",
        "the same brute-force Hessian equals 2 x the listed table (abs < 1e-4)",
        worst("doublet_block_vs_2x_table") < 1e-4,
        {"max_abs_deviation": worst("doublet_block_vs_2x_table"), "threshold": 1e-4},
    )
    record(
        "C6fd.director_constraint_term",
        "freezing the eigenvalues (unit n) adds -96 to the (d1, d1) and (d2, d2) Hessian "
        "entries and nothing else (abs < 1e-2): the -48 of the Galerkin director block",
        worst("director_block_dev_from_minus_96") < 1e-2,
        {"max_abs_deviation": worst("director_block_dev_from_minus_96"), "threshold": 1e-2},
    )
    record(
        "C6fd.mixed_block_divergence_form",
        "the constraint correction of the doublet x director block has the divergence form "
        "(abs < 1e-2 on a correction of size 16 to 22)",
        worst("mixed_block_dev_from_divergence_form") < 1e-2,
        {
            "max_abs_deviation": worst("mixed_block_dev_from_divergence_form"),
            "threshold": 1e-2,
        },
    )
    record(
        "C6fd.constrained_mixed_block_vanishes",
        "on the constrained family the doublet x director block vanishes pointwise "
        "(abs < 1e-2), and no linear term survives on the doublet jets",
        worst("constrained_mixed_block_max_abs") < 1e-2
        and worst("max_abs_linear_term_on_doublet_scaled") < 1e-2,
        {
            "max_abs_mixed_block": worst("constrained_mixed_block_max_abs"),
            "max_abs_linear_term_on_doublet": worst("max_abs_linear_term_on_doublet_scaled"),
            "threshold": 1e-2,
        },
    )


# ================= extras =================
def check_extras(rng_shared):
    """C3d nuance (which boosts the Euclidean construction survives), C4d on rotated / boosted
    backgrounds, C4e Frobenius cosh(4 theta) law, C2 candidate readings of the '31 to 39
    percent'. rng_shared is the stream left by check_c3 and check_c4."""
    print("\n=== extras ===")
    Nvac = np.diag([-G8, 1, DELTA, DELTA])
    print("C3d: Euclidean |B| on single boosts of the vacuum (rapidity 0.5):")
    single = {}
    for nm, a in (
        ("along e1 (the '1' axis)", [1, 0, 0]),
        ("along e2 (doublet axis)", [0, 1, 0]),
        ("along (e2+e3)/sqrt2", [0, 1, 1]),
        ("along (e1+e2)/sqrt2", [1, 1, 0]),
        ("along (1,1,1)/sqrt3", [1, 1, 1]),
    ):
        a = np.array(a, float)
        a /= np.linalg.norm(a)
        L = expm(0.5 * boostK(a))
        N = L @ Nvac @ np.linalg.inv(L)
        single[nm] = {
            "euclid": np.abs(eucl_B(N)).max(),
            "riesz": np.abs(riesz_B(N)).max(),
            "has_e1_component": bool(a[0] != 0),
        }
        print(
            "   %-26s Euclid %.3e   Riesz %.3e" % (nm, single[nm]["euclid"], single[nm]["riesz"])
        )
    record(
        "EX.c3d_single_boosts",
        "on single boosts at rapidity 0.5 the Riesz B vanishes (< 1e-10) for every axis, and "
        "the Euclidean B is nonzero exactly when the boost has a component along e1",
        all(v["riesz"] < 1e-10 for v in single.values())
        and all((v["euclid"] > 1e-3) == v["has_e1_component"] for v in single.values()),
        single,
    )

    print(
        "\nC4d on a ROTATED split background (spatial rotation only), random a: "
        "formula -2|B a|^2 th'^2 with B of that background"
    )
    rng = np.random.default_rng(3)
    worst = 0
    pos = 0
    for _ in range(100):
        F = rand_lorentz(rng_shared, scale_boost=0.0)
        tm, tf, thm, cl = run_case(
            rng.normal(size=3),
            rng.uniform(-0.25, 0.25),
            lambda x: 0.4 * np.sin(x),
            0.2,
            "rot",
            frame=F,
            show=False,
        )
        worst = max(worst, abs(tm + cl) / cl)
        pos += tm > 0
    print("   max rel err %.2e ; positive cases %d / 100" % (worst, pos))
    record(
        "EX.c4d_rotated_background",
        "the closed form holds on a spatially rotated split background (rel < 1e-4, never "
        "positive)",
        worst < 1e-4 and pos == 0,
        {"max_rel_error": worst, "positive_cases_of_100": int(pos)},
    )
    print(
        "C4d on a BOOSTED split background (lab-frame boost generator is then boost+rotation "
        "in the rest frame):"
    )
    boosted = {}
    for rap in (0.1, 0.5, 1.0):
        pos = 0
        vals = []
        for _ in range(100):
            b = rng.normal(size=3)
            b /= np.linalg.norm(b)
            F = expm(rap * boostK(b))
            tm, tf, thm, cl = run_case(
                rng.normal(size=3), 0.1, lambda x: 0.4 * np.sin(x), 0.2, "b", frame=F, show=False
            )
            pos += tm > 0
            vals.append(tm)
        boosted["rapidity %.1f" % rap] = {
            "positive_cases_of_100": int(pos),
            "range": [min(vals), max(vals)],
        }
        print(
            "   background rapidity %.1f: tr(dBdB) > 0 in %d / 100 ; range [%.2e, %.2e]"
            % (rap, pos, min(vals), max(vals))
        )
    record(
        "EX.c4d_boosted_background",
        "tr(dB dB) stays negative on a lab-frame boost ripple when the split background is "
        "itself boosted (rapidity 0.1, 0.5, 1.0)",
        all(v["positive_cases_of_100"] == 0 for v in boosted.values()),
        boosted,
        kind="scope",
        note="the sign is a statement about boosts in the rest frame of the background; a "
        "lab-frame boost of a moving background carries a rotation part of the other sign",
    )

    print("\nC4e: Frobenius / (2 eps^2 th'^2) against cosh(4 theta), a along a doublet axis")
    law = {}
    for th0 in (0.0, 0.02, 0.05, 0.2, 0.5, 1.0):
        tm, tf, thm, cl = run_case(
            np.array([0, 1.0, 0]), 0.1, lambda x, t=th0: t + 0.3 * x, 0.0, "f", show=False
        )
        law["theta %.2f" % th0] = {
            "frobenius_over_claim": tf / cl,
            "cosh_4_theta": np.cosh(4 * th0),
            "h_metric_over_claim": thm / cl,
            "covariant_over_claim": tm / cl,
        }
        print(
            "   theta %.2f: Frobenius/claim %.6f  cosh(4theta) %.6f  h-metric/claim %.6f  "
            "mixed/claim %.6f" % (th0, tf / cl, np.cosh(4 * th0), thm / cl, tm / cl)
        )
    record(
        "EX.c4e_frobenius_equals_claim",
        "Frobenius / (2 |B a|^2 theta'^2) = 1 for theta up to 1 (rel < 1e-4)",
        max(abs(v["frobenius_over_claim"] - 1) for v in law.values()) < 1e-4,
        law,
        note="refuted at finite theta: the ratio is 1.0032 at theta 0.02, 3.76 at 0.5, 27.3 "
        "at 1.0; the claim holds at theta = 0 only",
    )
    record(
        "EX.c4e_cosh_law",
        "the Frobenius ratio follows cosh(4 theta) while the h-metric ratio stays +1 and the "
        "covariant ratio stays -1 (rel < 1e-5)",
        max(
            max(
                abs(v["frobenius_over_claim"] / v["cosh_4_theta"] - 1),
                abs(v["h_metric_over_claim"] - 1),
                abs(v["covariant_over_claim"] + 1),
            )
            for v in law.values()
        )
        < 1e-5,
        {"threshold": 1e-5},
    )

    print(
        "\nC2: candidate readings of 'off by 31 to 39 percent' "
        "(A=1, r 6..18, beta=5.81e-3/(2*3.92))"
    )
    beta = 5.81e-3 / (2 * 3.92)
    nu = np.sqrt(5) / 4
    r = np.linspace(6, 18, 121)
    rho = np.sqrt(beta) * r**2 / 2
    X = r ** (2 * nu) * kv(nu, rho)
    cands = {}
    for nm, Y in (
        ("sqrt(rho)K_nu, const=1", np.sqrt(rho) * kv(nu, rho)),
        (
            "sqrt(rho)K_nu x (2/sqrt beta)^nu",
            np.sqrt(rho) * kv(nu, rho) * (2 / np.sqrt(beta)) ** nu,
        ),
        (
            "sqrt(rho)K_nu x (2/sqrt beta)^(1/2)",
            np.sqrt(rho) * kv(nu, rho) * (2 / np.sqrt(beta)) ** 0.5,
        ),
        (
            "sqrt(rho)K_nu normalized at rho->0 limit of rho^nu K_nu",
            np.sqrt(rho) * kv(nu, rho) * (2 / np.sqrt(beta)) ** nu,
        ),
    ):
        q = Y / X
        cands[nm] = {
            "Y_over_X": [q.min(), q.max()],
            "X_over_Y": [(1 / q).min(), (1 / q).max()],
        }
        print(
            "   %-58s Y/X in [%.4f, %.4f] ; X/Y in [%.4f, %.4f]"
            % (nm, q.min(), q.max(), (1 / q).min(), (1 / q).max())
        )
    slopes = {}
    for A_ in (1.0, 3.0):
        n_ = np.sqrt(1 + 4 * A_) / 4
        slopes["A=%g" % A_] = {
            "rho_pow_nu_minus_half_window": [rho.min() ** (n_ - 0.5), rho.max() ** (n_ - 0.5)],
            "log_slope_mismatch": n_ - 0.5,
        }
        print(
            "   A=%g: rho^(nu-1/2) over the window: [%.4f, %.4f]; log-slope mismatch "
            "d ln/d ln rho = %.4f"
            % (A_, rho.min() ** (n_ - 0.5), rho.max() ** (n_ - 0.5), n_ - 0.5)
        )
    hit = []
    for nm, v in cands.items():
        for key in ("Y_over_X", "X_over_Y"):
            lo, hi = v[key][0] - 1, v[key][1] - 1
            if _matches_31_39(lo, hi):
                hit.append(nm + " " + key)
    record(
        "EX.c2_candidate_readings",
        "some normalization of sqrt(rho) K_nu against r^(2 nu) K_nu reproduces 'off by 31 to "
        "39 percent' on r 6 to 18 (both ends within 3 points)",
        len(hit) > 0,
        {"candidates": cands, "log_slopes": slopes, "readings_that_reproduce_the_figure": hit},
        note="not reproduced: the closest reading, the (2 / sqrt beta)^(1/2) normalization, "
        "gives 24 to 41 percent; supplements C2.off_by_31_to_39_percent",
    )


# ================= main =================
def write_record(mode, runtime):
    rec = {}
    if os.path.exists(OUT_JSON):
        with open(OUT_JSON) as f:
            rec = json.load(f)
    rec["audit"] = "M5.32 R24-0 adversarial audit of the form-level claims C1 to C6"
    rec["verdict_convention"] = (
        "PASS or FAIL against the claim line of each check; FAIL = refuted or not reproduced; "
        "kind scope = an extension beyond the stated conditions of a claim"
    )
    checks = {k: v for k, v in rec.get("checks", {}).items() if v.get("mode") != mode}
    for cid, c in CHECKS.items():
        checks[cid] = dict(c, mode=mode)
    rec["checks"] = checks
    n, p, f, s = _counts(CHECKS)
    rec.setdefault("modes", {})[mode] = {
        "runtime_s": round(runtime, 1),
        "checks": n,
        "PASS": p,
        "FAIL": f,
        "SKIPPED": s,
        "fail_ids": [k for k, v in CHECKS.items() if v["verdict"] == "FAIL"],
    }
    n, p, f, s = _counts(checks)
    rec["summary"] = {
        "checks": n,
        "PASS": p,
        "FAIL": f,
        "SKIPPED": s,
        "fail_ids": [k for k, v in checks.items() if v["verdict"] == "FAIL"],
    }
    os.makedirs(DATA, exist_ok=True)
    with open(OUT_JSON, "w") as fh:
        json.dump(rec, fh, indent=1)
    return rec


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else ""
    if mode not in ("run", "stack"):
        print("usage: python3 m5_32_r24_0_audit.py run | stack")
        return 2
    if mode == "run":
        rng_shared = np.random.default_rng(20260920)  # one stream: C3, then C4, then extras
        steps = (
            ("C1", check_c1),
            ("C2", check_c2),
            ("C3", lambda: check_c3(rng_shared)),
            ("C4", lambda: check_c4(rng_shared)),
            ("C6 pointwise", check_c6_pointwise),
            ("C6 Galerkin", check_c6_galerkin),
            ("C6 massive Schur", check_c6_massive_schur),
            ("C6 constrained finite differences", check_c6_constrained_fd),
            ("extras", lambda: check_extras(rng_shared)),
        )
        for name, fun in steps:
            t = time.time()
            fun()
            print("[%7.1fs] %s done in %.1f s" % (time.time() - T0, name, time.time() - t))
    else:
        t = time.time()
        c5 = check_c5()
        print("[%7.1fs] C5 done in %.1f s" % (time.time() - T0, time.time() - t))
        t = time.time()
        check_vacuum_artifact(c5)
        print("[%7.1fs] vacuum artifact done in %.1f s" % (time.time() - T0, time.time() - t))
    rec = write_record(mode, time.time() - T0)
    print("\nFAIL checks of this mode:")
    if not any(c["verdict"] == "FAIL" for c in CHECKS.values()):
        print("  none")
    for cid, c in CHECKS.items():
        if c["verdict"] == "FAIL":
            print("  %s (%s): %s" % (cid, c["kind"], c["claim"]))
    s = rec["summary"]
    print(
        "record, all modes merged: %d checks, %d PASS, %d FAIL, %d SKIPPED"
        % (s["checks"], s["PASS"], s["FAIL"], s["SKIPPED"])
    )
    print(
        "runtime %.1f s, mode %s, wrote ../data/%s"
        % (time.time() - T0, mode, os.path.basename(OUT_JSON))
    )
    n, p, f, sk = _counts(CHECKS)
    tail = ", %d SKIPPED" % sk if sk else ""
    print("R24-0 audit: %d checks, %d PASS, %d FAIL%s" % (n, p, f, tail))
    return 0


if __name__ == "__main__":
    sys.exit(main())
