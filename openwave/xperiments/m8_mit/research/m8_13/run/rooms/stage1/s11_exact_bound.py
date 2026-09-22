"""s11_exact_bound.py -- the exact inequality on the elliptope.

Goal: show  R = 20 N /(77 Den^2) <= 463/924  on E = {detG = 1 + 2p - q1 >= 0},
i.e.  Phi := 463 Den^2 - 240 N >= 0  there, and find the equality locus.
"""

from pathlib import Path
import sympy as sp
import itertools
import numpy as np

HERE = Path(__file__).parent
a, b, c, L = sp.symbols("a b c L", real=True)

q1 = a**2 + b**2 + c**2
q2 = a**2 * b**2 + b**2 * c**2 + c**2 * a**2
p = a * b * c
Den = 3 * q1 - 2 * p + 5
Nn = 15 * q1**2 + 7 * q2 - 20 * p * q1 + 9 * p**2 + 71 * q1 - 108 * p + 30
detG = 1 + 2 * p - q1
R = sp.Rational(20, 77) * Nn / Den**2


def main():
    Phi = sp.expand(463 * Den**2 - 240 * Nn)
    print("Phi = 463 Den^2 - 240 N =")
    print("   ", sp.factor(Phi))
    print("\n  expanded:", Phi)
    print("\n  Phi at the hexagon (1/2,-1/2,1/2):",
          sp.simplify(Phi.subs({a: sp.Rational(1, 2), b: sp.Rational(-1, 2), c: sp.Rational(1, 2)})))
    print("  Phi at (0,0,0):", Phi.subs({a: 0, b: 0, c: 0}))
    print("  Phi at (1,1,1):", Phi.subs({a: 1, b: 1, c: 1}))

    # is Phi >= 0 on the whole cube, or only on the elliptope?
    print("\n  scan Phi on a grid of the CUBE [-1,1]^3 and of the ELLIPTOPE:")
    Pf = sp.lambdify((a, b, c), Phi, "numpy")
    Df = sp.lambdify((a, b, c), detG, "numpy")
    g = np.linspace(-1, 1, 61)
    A, B, C = np.meshgrid(g, g, g, indexing="ij")
    PV = Pf(A, B, C)
    DV = Df(A, B, C)
    print("    min Phi over cube grid      : %.6f  at %s" %
          (PV.min(), np.round([A.flat[PV.argmin()], B.flat[PV.argmin()], C.flat[PV.argmin()]], 4)))
    mask = DV >= 0
    print("    min Phi over elliptope grid : %.6f  at %s" %
          (PV[mask].min(), np.round([A[mask][PV[mask].argmin()], B[mask][PV[mask].argmin()],
                                     C[mask][PV[mask].argmin()]], 4)))

    # try: is Phi - lambda * detG * (something) a sum of squares?
    print("\n  Phi as a polynomial in detG:")
    dd = sp.symbols("d")
    # substitute q1 = 1 + 2p - d  to eliminate q1 in favour of detG
    print("   writing q1 = 1 + 2p - d, Den and N become:")
    Den_d = sp.expand(Den.subs(q1, 1 + 2 * p - dd))
    print("     Den =", sp.simplify(3 * (1 + 2 * p - dd) - 2 * p + 5))
    # careful: q2 is not a function of q1 and p alone -> keep symbolic
    print("     (q2 is independent of q1,p, so keep it)")

    Q1, Q2, P = sp.symbols("Q1 Q2 P", real=True)
    Den_s = 3 * Q1 - 2 * P + 5
    N_s = 15 * Q1**2 + 7 * Q2 - 20 * P * Q1 + 9 * P**2 + 71 * Q1 - 108 * P + 30
    Phi_s = sp.expand(463 * Den_s**2 - 240 * N_s)
    print("\n  Phi in terms of Q1,Q2,P:")
    print("   ", sp.collect(Phi_s, [Q2]))
    print("   coefficient of Q2:", sp.diff(Phi_s, Q2))

    # so Phi = -1680 Q2 + (stuff in Q1,P).  Need an upper bound for Q2.
    print("\n  Phi = -1680*Q2 + G(Q1,P) with G =", sp.expand(Phi_s + 1680 * Q2))

    # Newton-type inequality: Q2 = e2 of (a^2,b^2,c^2), Q1 = e1, P^2 = e3
    print("\n  With x=a^2,y=b^2,z=c^2 >= 0: Q1 = e1, Q2 = e2, P^2 = e3.")
    print("  Maclaurin/Newton: e2 <= e1^2/3, and e2 >= ... need the RIGHT bound.")
    print("  Since Phi has -1680 Q2, Phi >= 0 needs an UPPER bound on Q2 at given Q1,P.")
    x, y, z = sp.symbols("x y z", nonnegative=True)
    print("   e2 <= e1^2/3 gives Phi >= -560 Q1^2 + G(Q1,P) =",
          sp.expand(sp.simplify(-560 * Q1**2 + (Phi_s + 1680 * Q2))))

    # evaluate that weaker bound on the elliptope numerically
    expr = sp.expand(-560 * Q1**2 + (Phi_s + 1680 * Q2))
    f = sp.lambdify((Q1, P), expr, "numpy")
    print("\n   scan of the weakened bound over feasible (Q1,P):")
    worst = (1e9, None)
    for _ in range(200000):
        V = np.random.normal(size=(3, 3))
        V /= np.linalg.norm(V, axis=1)[:, None]
        aa, bb, cc = V[0] @ V[1], V[0] @ V[2], V[1] @ V[2]
        val = float(f(aa**2 + bb**2 + cc**2, aa * bb * cc))
        if val < worst[0]:
            worst = (val, (aa, bb, cc))
    print("    min of weakened bound over sampled elliptope: %.6f at %s"
          % (worst[0], np.round(worst[1], 6)))


if __name__ == "__main__":
    main()
