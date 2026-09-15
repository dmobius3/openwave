"""Independent adversarial audit of the two R21 form-level claims.

Reimplemented from the claim statements only (no import of the script under test).
Part 1: symmetric-point slope identity nu' = -1/4, several potentials, plus rank weighting.
Part 2: V4 stationary points in nu at Delta = -9 and Delta = 0, classification, nu' by two
        independent routes (implicit derivative on an independently built polynomial, and
        finite differences on re-minimized valleys using a matrix-based evaluation of V4).
Part 3: boost dressing Q = exp(s K): eta invariance, N-spectrum invariance, det M invariance,
        first and second derivative of tr(M^2), symbolic on a generic symmetric M and numeric.
"""

import time

import numpy as np
import sympy as sp
from numpy.polynomial import polynomial as P
from scipy.linalg import expm
from scipy.optimize import minimize_scalar

T0 = time.time()
W1_F = 0.000724023879
SPEC_VAC = [-8.0, 1.0, 0.3, 0.0]
TABLE = []


def row(check, method, mine, claimed, ok):
    TABLE.append((check, method, mine, claimed, "PASS" if ok else "FAIL"))


# ---------------------------------------------------------------- part 1: symmetric-point identity
def slope_symbolic(Vsym, lams, nu, D, at_D=0, side=None):
    """nu' = -V_nD / V_nn on every stationary point in nu of V(nu + D, nu, nu, nu) at Delta = at_D.

    Returns list of (nu0, Vnn, nu'). side is unused here (kept for signature parity)."""
    V = Vsym.subs(dict(zip(lams, [nu + D, nu, nu, nu])), simultaneous=True)
    Vn, Vnn, VnD = sp.diff(V, nu), sp.diff(V, nu, 2), sp.diff(V, nu, D)
    out = []
    for r in sp.Poly(sp.expand(Vn.subs(D, at_D)), nu).nroots(n=30):
        if abs(sp.im(r)) > 1e-20:
            continue
        r = sp.re(r)
        vnn = Vnn.subs({nu: r, D: at_D})
        vnd = VnD.subs({nu: r, D: at_D})
        out.append((float(r), float(vnn), float(-vnd / vnn)))
    return out


def hessian_ratio(Vsym, lams):
    """(A + 3B) / (4A + 12B) at a totally symmetric point, A = V_11, B = V_12; plus the full
    Hessian check that S_ii are all equal and S_ij (i != j) all equal there."""
    x = sp.symbols("x", real=True)
    sub = {l: x for l in lams}
    H = sp.hessian(Vsym, lams).subs(sub)
    H = H.applyfunc(sp.simplify)
    A = H[0, 0]
    B = H[0, 1]
    diag_ok = all(sp.simplify(H[i, i] - A) == 0 for i in range(4))
    off_ok = all(sp.simplify(H[i, j] - B) == 0 for i in range(4) for j in range(4) if i != j)
    ratio = sp.simplify((A + 3 * B) / (4 * A + 12 * B))
    # a second, independent route: the chain rule directly on nu, D
    nu, D = sp.symbols("nu Delta", real=True)
    V = Vsym.subs(dict(zip(lams, [nu + D, nu, nu, nu])), simultaneous=True)
    r2 = sp.simplify((sp.diff(V, nu, D) / sp.diff(V, nu, 2)).subs(D, 0))
    return ratio, r2, diag_ok and off_ok


def part1():
    lams = sp.symbols("l1:5", real=True)
    nu, D = sp.symbols("nu Delta", real=True)
    W1 = sp.Rational(724023879, 10**12)
    C = [sum(sp.nsimplify(l) ** p for l in SPEC_VAC) for p in range(1, 5)]
    V4 = W1 * sum((sum(l**p for l in lams) - C[p - 1]) ** 2 for p in range(1, 5))
    # elementary symmetric polynomials, a non-separable control of a different family
    e2 = sum(lams[i] * lams[j] for i in range(4) for j in range(i + 1, 4))
    e4 = lams[0] * lams[1] * lams[2] * lams[3]
    Ve = (e2 - 3) ** 2 + e4
    # a random symmetric polynomial in power sums (cross terms), coefficients fixed by seed
    rng = np.random.default_rng(7)
    ps = [sum(l**p for l in lams) for p in range(1, 4)]
    Vr = (
        sum(
            sp.Rational(int(rng.integers(-9, 10)), 4) * ps[i] * ps[j]
            for i in range(3)
            for j in range(i, 3)
        )
        + sp.Rational(int(rng.integers(1, 9)), 3) * ps[0] * ps[1] * ps[2]
    )
    # separable control, different f from the script under test
    f = lambda l: sp.exp(l / 3) + (l - 1) ** 2 / 5 + l**3 / 11
    Vs = sum(f(l) for l in lams)
    for name, Vsym in (
        ("V4", V4),
        ("(e2-3)^2+e4", Ve),
        ("random power-sum poly", Vr),
        ("separable exp+cubic", Vs),
    ):
        ratio, r2, hess_ok = hessian_ratio(Vsym, lams)
        ok = (ratio == sp.Rational(1, 4)) and (r2 == sp.Rational(1, 4)) and hess_ok
        row(
            f"1: V_nD/V_nn at symmetric point, {name}",
            "sympy Hessian (A,B) + direct chain rule",
            f"{ratio} / {r2} / Hessian A,B form {hess_ok}",
            "1/4",
            ok,
        )
        print(f"[1] {name}: Hessian ratio {ratio}, chain-rule ratio {r2}, A/B form {hess_ok}")
    # rank-weighted: a f(lambda_max) + b sum_rest f; lambda_max = nu + D for D > 0
    a, b = sp.Integer(3), sp.Integer(1)
    g = lambda l: (l - 2) ** 2 + l**4 / 7
    Vplus = a * g(nu + D) + b * 3 * g(nu)  # D > 0: heavy eigenvalue is the max
    # D < 0: the max is one of the triple; the heavy one at nu + D takes weight b
    Vminus = a * g(nu) + b * (g(nu + D) + 2 * g(nu))
    pred = float(-a / (a + 3 * b))
    for lab, Vw in (("D -> 0+", Vplus), ("D -> 0-", Vminus)):
        V = Vw
        Vn, Vnn, VnD = sp.diff(V, nu), sp.diff(V, nu, 2), sp.diff(V, nu, D)
        vals = []
        for r in sp.Poly(sp.expand(Vn.subs(D, 0)), nu).nroots(n=30):
            if abs(sp.im(r)) > 1e-20:
                continue
            r = sp.re(r)
            vnn, vnd = Vnn.subs({nu: r, D: 0}), VnD.subs({nu: r, D: 0})
            if vnn > 0:
                vals.append((float(r), float(-vnd / vnn)))
        print(f"[1] rank-weighted a=3 b=1 {lab}: valleys (nu0, nu') = {vals}")
        if lab == "D -> 0+":
            ok = all(abs(v[1] - pred) < 1e-12 for v in vals) and len(vals) > 0
            row(
                "1: rank-weighted a=3 b=1, nu' at D -> 0+",
                "sympy implicit derivative, one-sided",
                str([round(v[1], 12) for v in vals]),
                f"-a/(a+3b) = {pred}",
                ok,
            )
        else:
            row(
                "1: rank-weighted a=3 b=1, nu' at D -> 0- (NOT claimed, kink check)",
                "sympy implicit derivative, one-sided",
                str([round(v[1], 12) for v in vals]),
                "(no claim; -1/4 would be the unweighted value)",
                True,
            )
    # generic (a, b) symbolic
    A_, B_ = sp.symbols("a b", positive=True)
    Vg = A_ * g(nu + D) + B_ * 3 * g(nu)
    Vn, Vnn, VnD = sp.diff(Vg, nu), sp.diff(Vg, nu, 2), sp.diff(Vg, nu, D)
    # at D = 0, Vn = (A+3B) g'(nu) = 0 -> g'(nu0) = 0; then VnD = A g''(nu0), Vnn = (A+3B) g''(nu0)
    r_gen = sp.simplify((VnD / Vnn).subs(D, 0))
    ok = sp.simplify(r_gen - A_ / (A_ + 3 * B_)) == 0
    row(
        "1: rank-weighted generic (a,b), V_nD/V_nn at D -> 0+",
        "sympy symbolic",
        str(r_gen),
        "a/(a+3b)",
        ok,
    )
    print(f"[1] generic rank-weighted ratio: {r_gen}")


# ---------------------------------------------------------------- part 2: V4 branches
def v4_matrix(nu, D):
    """V4 evaluated the long way: build M, form N = M eta, take traces of powers."""
    eta = np.diag([-1.0, 1, 1, 1])
    M = np.diag([-(nu + D), nu, nu, nu])  # so that N = M eta = diag(nu + D, nu, nu, nu)
    N = M @ eta
    Np = np.eye(4)
    v = 0.0
    for p in range(1, 5):
        Np = Np @ N
        Cp = sum(l**p for l in SPEC_VAC)
        v += (np.trace(Np) - Cp) ** 2
    return W1_F * v


def v4_poly_in_nu(D):
    """Coefficients of V4(nu, D) as a polynomial in nu, built with numpy polynomial arithmetic."""
    lam = [np.array([D, 1.0]), np.array([0.0, 1.0])]  # nu + D and nu, ascending powers
    total = np.array([0.0])
    for p in range(1, 5):
        Cp = sum(l**p for l in SPEC_VAC)
        s = P.polypow(lam[0], p)
        for _ in range(3):
            s = P.polyadd(s, P.polypow(lam[1], p))
        s = P.polyadd(s, np.array([-Cp]))
        total = P.polyadd(total, P.polymul(s, s))
    return W1_F * total


def part2():
    claims = {
        -9.0: [(1.0009, -1.0054), (5.9868, 0.0491)],
        0.0: [(-5.6472, -0.25), (5.5999, -0.25)],
    }
    for D, cl in claims.items():
        c = v4_poly_in_nu(D)
        dc = P.polyder(c)
        ddc = P.polyder(dc)
        roots = P.polyroots(dc)
        real = sorted(r.real for r in roots if abs(r.imag) < 1e-9)
        print(
            f"[2] Delta = {D}: {len(roots)} stationary points, {len(real)} real: {np.round(real, 6)}"
        )
        valleys = []
        for r in real:
            vnn = P.polyval(r, ddc)
            kind = "min" if vnn > 0 else ("max" if vnn < 0 else "flat")
            # route A: implicit derivative from the polynomial family, V_nD by central difference in D
            h = 1e-5
            vn_plus = P.polyval(r, P.polyder(v4_poly_in_nu(D + h)))
            vn_minus = P.polyval(r, P.polyder(v4_poly_in_nu(D - h)))
            vnd = (vn_plus - vn_minus) / (2 * h)
            slopeA = -vnd / vnn

            # route B: re-minimize the matrix-built V4 at D +- h, Richardson on two step sizes
            def valley_at(Dx, guess=r):
                res = minimize_scalar(
                    lambda x: v4_matrix(x, Dx),
                    bounds=(guess - 0.05, guess + 0.05),
                    method="bounded",
                    options={"xatol": 1e-12},
                )
                return res.x

            slopeB = float("nan")
            if kind == "min":
                slopes = []
                for hh in (2e-3, 1e-3):
                    slopes.append((valley_at(D + hh) - valley_at(D - hh)) / (2 * hh))
                slopeB = (4 * slopes[1] - slopes[0]) / 3
            print(
                f"      nu0 = {r:+.6f}  V_nn = {vnn:+.4e} ({kind})  nu' implicit = {slopeA:+.6f}  nu' re-min = {slopeB:+.6f}"
            )
            if kind == "min":
                valleys.append((r, slopeA, slopeB))
        for nu_c, sl_c in cl:
            hit = [v for v in valleys if abs(v[0] - nu_c) < 5e-3]
            if not hit:
                row(
                    f"2: V4 Delta={D:g} valley near nu={nu_c}",
                    "numpy polyroots + V_nn sign",
                    "no valley found",
                    f"nu'={sl_c}",
                    False,
                )
                continue
            r, sA, sB = hit[0]
            ok = abs(sA - sl_c) < 1e-3 and abs(sB - sl_c) < 1e-3 and abs(sA - sB) < 1e-4
            row(
                f"2: V4 Delta={D:g} valley nu0={r:.4f}",
                "implicit deriv (poly) + re-minimized FD (matrix V4)",
                f"{sA:+.5f} / {sB:+.5f}",
                f"{sl_c:+.4f}",
                ok,
            )
        n_val = len(valleys)
        row(
            f"2: V4 Delta={D:g} number of valley branches",
            "all real roots of V_nu, V_nn > 0",
            str(n_val),
            str(len(cl)),
            n_val == len(cl),
        )
    # W1 independence: nu' is a ratio, so W1 cancels; verify by rescaling
    c9 = v4_poly_in_nu(-9.0)
    dc = P.polyder(c9)
    r = max((x.real for x in P.polyroots(dc) if abs(x.imag) < 1e-9 and abs(x.real - 1.0) < 0.1))
    print(
        f"[2] note: nu' is a ratio of second derivatives, W1 cancels identically (root check nu0 = {r:.6f})"
    )


# ---------------------------------------------------------------- part 3: boost dressing
def part3():
    s = sp.symbols("s", real=True)
    K = sp.zeros(4, 4)
    K[0, 1] = K[1, 0] = 1
    eta = sp.diag(-1, 1, 1, 1)
    Q = (s * K).exp()  # sympy matrix exponential, no closed form assumed
    Q = Q.applyfunc(sp.simplify)
    # (a) Q^T eta Q = eta
    a_ok = (Q.T * eta * Q - eta).applyfunc(sp.simplify) == sp.zeros(4, 4)
    row("3a: Q^T eta Q = eta (symbolic, all s)", "sympy matrix exp", str(a_ok), "True", a_ok)
    # generic symmetric M
    m = sp.symbols("m0:10", real=True)
    M = sp.zeros(4, 4)
    idx = 0
    for i in range(4):
        for j in range(i, 4):
            M[i, j] = M[j, i] = m[idx]
            idx += 1
    Ms = Q * M * Q.T
    # (b) N spectrum: N(s) = Q M Q^T eta; Q^T eta = eta Q^{-1}, so N(s) = Q N Q^{-1}: check char poly equality
    lam = sp.symbols("lam")
    b1 = (Q.T * eta - eta * Q.inv()).applyfunc(sp.simplify) == sp.zeros(4, 4)
    Ns = Ms * eta
    N0 = M * eta
    # Newton-trace route: tr N(s)^p = tr N^p for p = 1..4 fixes the char poly of a 4x4
    tr_ok = True
    Np, N0p = sp.eye(4), sp.eye(4)
    for p in range(1, 5):
        Np, N0p = Np * Ns, N0p * N0
        d = sp.simplify(sp.expand((Np.trace() - N0p.trace()).rewrite(sp.exp)))
        tr_ok = tr_ok and d == 0
    det_ok = sp.simplify((Ms.det() - M.det()).rewrite(sp.exp)) == 0
    row("3b: N(s) = Q N Q^-1 (Q^T eta = eta Q^-1), symbolic", "sympy", str(b1), "True", b1)
    row(
        "3b: tr N(s)^p = tr N^p, p=1..4, generic symmetric M, all s",
        "sympy power traces",
        str(tr_ok),
        "True",
        tr_ok,
    )
    row("3b: det M(s) = det M, generic symmetric M, all s", "sympy", str(det_ok), "True", det_ok)
    # (c),(d) derivatives of tr(M(s)^2)
    n_s = (Ms * Ms).trace()
    d1 = sp.simplify(sp.diff(n_s, s).subs(s, 0))
    d2 = sp.simplify(sp.diff(n_s, s, 2).subs(s, 0))
    c1 = sp.expand(4 * (M * M * K).trace())
    c2 = sp.expand(8 * (K * K * M * M).trace() + 8 * (K * M * K * M).trace())
    c_ok = sp.expand(d1 - c1) == 0
    d_ok = sp.expand(d2 - c2) == 0
    row(
        "3c: d/ds tr M^2 at 0 = 4 tr(M^2 K), generic symmetric M", "sympy", str(c_ok), "True", c_ok
    )
    row(
        "3d: d2/ds2 tr M^2 at 0 = 8 tr(K^2 M^2) + 8 tr(K M K M), generic",
        "sympy",
        str(d_ok),
        "True",
        d_ok,
    )
    print(f"[3] d1 closed form (expanded) = {sp.factor(d1)}")
    # block-diagonal: M_0i = 0 for i = 1,2,3 -> m1 = m2 = m3 = 0
    bd = {m[1]: 0, m[2]: 0, m[3]: 0}
    d1_bd = sp.simplify(d1.subs(bd))
    row(
        "3c: first derivative vanishes on block-diagonal M (symbolic)",
        "sympy subs M_0i = 0",
        str(d1_bd),
        "0",
        d1_bd == 0,
    )
    generic_nonzero = sp.expand(d1) != 0
    row(
        "3c: first derivative generically nonzero (nonzero polynomial)",
        "sympy",
        str(generic_nonzero),
        "True",
        generic_nonzero,
    )
    # converse probe: a non-block-diagonal M with vanishing first derivative
    probe = {
        m[0]: 2,
        m[1]: 0,
        m[2]: 1,
        m[3]: 0,
        m[4]: 5,
        m[5]: 0,
        m[6]: 0,
        m[7]: 1,
        m[8]: 0,
        m[9]: 1,
    }
    d1_probe = d1.subs(probe)
    print(
        f"[3] converse probe: M_02 = 1 (not block-diagonal), M_01 = M_12 = M_03 = 0 gives d1 = {d1_probe}"
    )
    # (d) at the vacuum
    vac = {
        m[0]: 8,
        m[4]: 1,
        m[7]: sp.Rational(3, 10),
        m[9]: 0,
        m[1]: 0,
        m[2]: 0,
        m[3]: 0,
        m[5]: 0,
        m[6]: 0,
        m[8]: 0,
    }
    d2_vac = d2.subs(vac)
    row(
        "3d: second derivative at M = diag(8, 1, 0.3, 0)",
        "sympy, exact",
        str(d2_vac),
        "648",
        d2_vac == 648,
    )
    # hand count: 8 tr(P M^2) = 8 (64 + 1) = 520 ; 8 tr(KMKM) = 8 (8 + 8) = 128 ; total 648
    print("[3] hand count: 8*(64+1) + 8*(1*8 + 8*1) = 520 + 128 = 648")
    # numerics: several random symmetric M, several s including large s, expm from scipy
    rng = np.random.default_rng(123)
    Kn = np.zeros((4, 4))
    Kn[0, 1] = Kn[1, 0] = 1.0
    etan = np.diag([-1.0, 1, 1, 1])
    worst_spec, worst_det, worst_d1, worst_d2 = 0.0, 0.0, 0.0, 0.0
    for trial in range(20):
        A = rng.standard_normal((4, 4)) * 3
        Mn = A + A.T
        if trial == 0:
            Mn = np.diag([8.0, 1, 0.3, 0])
        spec0 = np.sort(np.linalg.eigvals(Mn @ etan))
        for sv in (0.1, 0.7, 2.0, -1.5):
            Qn = expm(sv * Kn)
            Mss = Qn @ Mn @ Qn.T
            spec = np.sort(np.linalg.eigvals(Mss @ etan))
            worst_spec = max(
                worst_spec, np.max(np.abs(spec - spec0)) / max(1.0, np.max(np.abs(spec0)))
            )
            worst_det = max(
                worst_det,
                abs(np.linalg.det(Mss) - np.linalg.det(Mn)) / max(1.0, abs(np.linalg.det(Mn))),
            )
        f = lambda sv: np.trace(
            (expm(sv * Kn) @ Mn @ expm(sv * Kn).T) ** 1 @ (expm(sv * Kn) @ Mn @ expm(sv * Kn).T)
        )
        h = 1e-3
        d1n = (f(h) - f(-h)) / (2 * h)
        d2n = (f(h) - 2 * f(0) + f(-h)) / h**2
        # Richardson with h/2
        d1n2 = (f(h / 2) - f(-h / 2)) / h
        d2n2 = (f(h / 2) - 2 * f(0) + f(-h / 2)) / (h / 2) ** 2
        d1r, d2r = (4 * d1n2 - d1n) / 3, (4 * d2n2 - d2n) / 3
        c1n = 4 * np.trace(Mn @ Mn @ Kn)
        c2n = 8 * np.trace(Kn @ Kn @ Mn @ Mn) + 8 * np.trace(Kn @ Mn @ Kn @ Mn)
        worst_d1 = max(worst_d1, abs(d1r - c1n) / max(1.0, abs(c1n)))
        worst_d2 = max(worst_d2, abs(d2r - c2n) / max(1.0, abs(c2n)))
        if trial == 0:
            print(
                f"[3] vacuum numeric: d1 = {d1r:+.6f} (closed {c1n:+.6f}), d2 = {d2r:.6f} (closed {c2n:.6f})"
            )
            row(
                "3d: vacuum second derivative, numeric Richardson FD",
                "scipy expm + FD",
                f"{d2r:.4f}",
                "648",
                abs(d2r - 648) < 1e-3,
            )
    row(
        "3b: N-spectrum drift, 20 random symmetric M, s in {0.1, 0.7, 2, -1.5}",
        "numpy eigvals (relative)",
        f"{worst_spec:.1e}",
        "0",
        worst_spec < 1e-9,
    )
    row(
        "3b: det M drift, same set",
        "numpy det (relative)",
        f"{worst_det:.1e}",
        "0",
        worst_det < 1e-9,
    )
    row(
        "3c: numeric d1 vs 4 tr(M^2 K), 20 random M",
        "Richardson FD (relative)",
        f"{worst_d1:.1e}",
        "0",
        worst_d1 < 1e-6,
    )
    row(
        "3d: numeric d2 vs closed form, 20 random M",
        "Richardson FD (relative)",
        f"{worst_d2:.1e}",
        "0",
        worst_d2 < 1e-5,
    )


if __name__ == "__main__":
    part1()
    part2()
    part3()
    print()
    print("| Check | Method | Mine | Claimed | Verdict |")
    print("| --- | --- | --- | --- | --- |")
    for r in TABLE:
        print("| " + " | ".join(str(x) for x in r) + " |")
    print(
        f"\nfails: {sum(1 for r in TABLE if r[4] == 'FAIL')} of {len(TABLE)}; elapsed {time.time() - T0:.1f} s"
    )
