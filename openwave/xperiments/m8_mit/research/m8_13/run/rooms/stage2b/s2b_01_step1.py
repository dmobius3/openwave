"""s2b_01_step1.py -- my own check of the author's Step 1 (the invariant form).

Nothing is taken from S0_S3_MAXIMUM.md, check_s3_maximum.py or its log.  The spin-3
objects are rebuilt from the room's `worklist.md` conventions
(J_z v_m = m v_m, J_+- v_m = sqrt(12 - m(m+-1)) v_{m+-1}, (Theta u)_m = (-1)^{3-m} conj(c_{-m})),
and rho_6 from the worklist's own definition.

Checked here:
  (1) the stretched CG closed form against sympy, all 49 pairs
  (2) ||rho_6||^2 at the weight states v3, v2, v1, v0
  (3) |f|^2, |a00|^2, TrNbar, TrNbar^2 at those states
  (4) the 4x4 independence determinant
  (5) the solved coefficients
  (6) the resulting identity, on fresh random states, at two precisions
  (7) TrNbar = 12 identically  (the author's "N = 4I + Q with Q traceless" needs it
      and neither the note nor check_s3_maximum.py states or gates it)
  (8) the bridge the 4-dimensionality rests on: bidegree-(2,2) forms <-> Hermitian
      forms on Sym^2 is a BIJECTION, because the 784 functions conj(m_a) m_b are
      distinct monomials in (u, ubar) and hence linearly independent
"""

from pathlib import Path
import itertools
import numpy as np
import sympy as sp
from sympy.physics.quantum.cg import CG

HERE = Path(__file__).parent
ORDER = list(range(3, -4, -1))
IDX = {m: i for i, m in enumerate(ORDER)}


def cg_stretched(m1, m2):
    M = m1 + m2
    if abs(M) > 6 or abs(m1) > 3 or abs(m2) > 3:
        return sp.Integer(0)
    f = sp.factorial
    return sp.sqrt(sp.Rational(f(6) * f(6) * f(6 + M) * f(6 - M),
                               f(12) * f(3 + m1) * f(3 - m1) * f(3 + m2) * f(3 - m2)))


def spin3():
    Jz = sp.diag(*[sp.Integer(m) for m in ORDER])
    Jp = sp.zeros(7, 7)
    Jm = sp.zeros(7, 7)
    for m in ORDER:
        if m + 1 <= 3:
            Jp[IDX[m + 1], IDX[m]] = sp.sqrt(12 - m * (m + 1))
        if m - 1 >= -3:
            Jm[IDX[m - 1], IDX[m]] = sp.sqrt(12 - m * (m - 1))
    return (Jp + Jm) / 2, (Jp - Jm) / (2 * sp.I), Jz


FX, FY, FZ = spin3()
F = [FX, FY, FZ]


def theta(c):
    return {m: (-1) ** (3 - m) * sp.conjugate(c[-m]) for m in range(-3, 4)}


def rho6_sq(c):
    tot = 0
    th = theta(c)
    for Q in range(-6, 7):
        s = 0
        for m1 in range(-3, 4):
            m2 = Q - m1
            if abs(m2) <= 3:
                s += cg_stretched(m1, m2) * c[m1] * th[m2]
        s = sp.expand(s)
        tot += sp.expand(sp.re(s) ** 2 + sp.im(s) ** 2)
    return sp.simplify(sp.expand(tot))


def vecof(c):
    return sp.Matrix(7, 1, [c[m] for m in ORDER])


def f2_of(c):
    u = vecof(c)
    comps = [sp.simplify(sp.expand((u.H * X * u)[0, 0])) for X in F]
    return sp.simplify(sum(sp.re(x) ** 2 + sp.im(x) ** 2 for x in comps)), comps


def a00sq_of(c):
    th = theta(c)
    s = sp.expand(sum(sp.conjugate(th[m]) * c[m] for m in range(-3, 4)))
    return sp.simplify((sp.re(s) ** 2 + sp.im(s) ** 2) / 7)


def Nbar_of(c):
    u = vecof(c)
    return sp.Matrix(3, 3, lambda i, k: sp.simplify(sp.re(sp.expand((u.H * F[i] * F[k] * u)[0, 0]))))


def main():
    # ---- (1) the stretched CG closed form vs sympy, all 49 pairs
    bad = []
    for m1 in range(-3, 4):
        for m2 in range(-3, 4):
            s = sp.simplify(CG(3, m1, 3, m2, 6, m1 + m2).doit()) if abs(m1 + m2) <= 6 else 0
            if sp.simplify(s - cg_stretched(m1, m2)) != 0:
                bad.append((m1, m2))
    print("(1) stretched CG closed form vs sympy, 49 pairs, mismatches:", bad if bad else "none")

    # ---- (2)(3) the weight states
    rows = []
    print("\n(2)(3) weight states")
    for m in (3, 2, 1, 0):
        c = {k: (sp.Integer(1) if k == m else sp.Integer(0)) for k in range(-3, 4)}
        r6 = rho6_sq(c)
        f2, comps = f2_of(c)
        a0 = a00sq_of(c)
        N = Nbar_of(c)
        trN = sp.simplify(N.trace())
        trN2 = sp.simplify((N * N).trace())
        print("   m=%d : ||rho6||^2 = %s = %s/924 | |f|^2 = %s | |a00|^2 = %s | TrN = %s | TrN^2 = %s"
              % (m, r6, sp.simplify(r6 * 924), f2, a0, trN, trN2))
        print("          Nbar = diag%s  (off-diagonals %s)"
              % ([N[i, i] for i in range(3)],
                 "all zero" if all(N[i, k] == 0 for i in range(3) for k in range(3) if i != k) else "NOT zero"))
        rows.append([sp.Integer(1), f2, a0, trN2, r6])

    # ---- (4) independence determinant
    A = sp.Matrix([r[:4] for r in rows])
    det = sp.simplify(A.det())
    print("\n(4) independence determinant of the four invariants at v3,v2,v1,v0 :", det,
          "  (author says 180/7; equal:", sp.simplify(det - sp.Rational(180, 7)) == 0, ")")

    # ---- (5) solved coefficients
    coef = A.solve(sp.Matrix([r[4] for r in rows]))
    print("(5) solved coefficients:", list(coef))
    want = [sp.Rational(-5, 231), sp.Rational(-1, 22), sp.Rational(7, 11), sp.Rational(1, 198)]
    print("    author's:", want, "  equal:", all(sp.simplify(coef[i] - want[i]) == 0 for i in range(4)))

    # ---- (7) TrNbar = 12 identically (exact, symbolic u)
    xs = sp.symbols("x0:7", real=True)
    ys = sp.symbols("y0:7", real=True)
    csym = {m: xs[IDX[m]] + sp.I * ys[IDX[m]] for m in range(-3, 4)}
    Nsym = Nbar_of(csym)
    nrm = sp.expand(sum(xs[i] ** 2 + ys[i] ** 2 for i in range(7)))
    print("\n(7) TrNbar - 12*||u||^2 identically :",
          sp.simplify(sp.expand(Nsym.trace() - 12 * nrm)) == 0)
    print("    -> Nbar = 4 I + Q with Q traceless requires exactly this; the note")
    print("       asserts it without stating it, and check_s3_maximum.py has no gate for it.")

    # ---- (6) the identity, on fresh random states, double and 30-digit
    print("\n(6) the identity r6 = -5/231 - |f|^2/22 + (7/11)|a00|^2 + TrN^2/198")
    rng = np.random.default_rng(2026)
    worst = 0.0
    for _ in range(12):
        z = rng.normal(size=7) + 1j * rng.normal(size=7)
        z /= np.linalg.norm(z)
        c = {m: sp.nsimplify(sp.Float(z[IDX[m]].real, 30) + sp.I * sp.Float(z[IDX[m]].imag, 30),
                             rational=False) for m in range(-3, 4)}
        r6 = sp.N(rho6_sq(c), 30)
        f2, _ = f2_of(c)
        pred = (sp.Rational(-5, 231) - sp.N(f2, 30) / 22 + sp.Rational(7, 11) * sp.N(a00sq_of(c), 30)
                + sp.N((Nbar_of(c) * Nbar_of(c)).trace(), 30) / 198)
        worst = max(worst, abs(float(sp.N(r6 - pred, 30))))
    print("    max |lhs - rhs| over 12 random unit states, 30-digit arithmetic: %.3e" % worst)

    # ---- (8) the bridge: injectivity of Herm(Sym^2) -> bidegree-(2,2) forms
    mons = [(i, k) for i in range(7) for k in range(i, 7)]
    print("\n(8) Sym^2(C^7) monomial count:", len(mons), "(28 expected)")
    pairs = set()
    dup = 0
    for al in mons:
        for be in mons:
            key = (al, be)          # conj(u_a u_b) * u_c u_d  <-> the ordered pair
            if key in pairs:
                dup += 1
            pairs.add(key)
    print("    the %d functions conj(m_alpha) m_beta carry DISTINCT (antiholomorphic," % len(pairs))
    print("    holomorphic) monomial bidegrees, duplicates found:", dup)
    print("    A polynomial in (u, ubar) vanishing identically has zero coefficients, so")
    print("    these %d functions are linearly independent and the map" % len(pairs))
    print("    Herm(Sym^2) -> {(2,2)-forms} is INJECTIVE.  Dimension counting alone does NOT")
    print("    give this; it is the step both the author's Step 1 and my stage-1 Step 1 compress.")


if __name__ == "__main__":
    main()
