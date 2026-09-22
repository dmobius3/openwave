"""s12_chain_verify.py -- verification of the exact proof chain for the MAXIMUM,
and of the equality analysis for the MINIMUM.

MAXIMUM chain (all steps exact, verified symbolically here):

  Phi := 463*Den^2 - 240*N  =  -1680*Q2 + G(Q1,P),
         G = 567 Q1^2 - 308 P^2 - 756 P Q1 + 16660 P - 3150 Q1 + 4375
  (S1) Maclaurin on x=a^2,y=b^2,z=c^2 >= 0 :   Q2 = e2 <= e1^2/3 = Q1^2/3
       =>  Phi >= Psi(Q1,P) := 7 Q1^2 - 308 P^2 - 756 P Q1 + 16660 P - 3150 Q1 + 4375
  (S2) dPsi/dQ1 = 14 Q1 - 756 P - 3150 < 0 on the feasible box (Q1 <= 3, P >= -1/8),
       and Q1 <= 1 + 2P (Gram PSD)  =>  Psi(Q1,P) >= Psi(1+2P, P) =: psi(P)
  (S3) psi(P) = -1792 P^2 + 9632 P + 1232 = 1792 (P + 1/8)(11/2 - P)
  (S4) P = abc >= -1/8 for any three unit vectors:
       AM-GM gives Q1 >= 3 |P|^(2/3); with Q1 <= 1+2P and P = -sigma^3 < 0 this is
       1 - 2 sigma^3 >= 3 sigma^2, i.e. (2 sigma - 1)(sigma+1)^2 <= 0, i.e. sigma <= 1/2.
  => Phi >= 0, i.e. rhat_6 <= 463/924.

EQUALITY forces: Q2 = Q1^2/3 (=> a^2=b^2=c^2), Q1 = 1+2P (coplanar), P = -1/8,
hence |a|=|b|=|c|=1/2 and abc = -1/8: three coplanar axes at 60 degrees.
"""

from pathlib import Path
import sympy as sp
import numpy as np

HERE = Path(__file__).parent
a, b, c = sp.symbols("a b c", real=True)
Q1, Q2, P = sp.symbols("Q1 Q2 P", real=True)
sig = sp.symbols("sigma", nonnegative=True)

q1e = a**2 + b**2 + c**2
q2e = a**2 * b**2 + b**2 * c**2 + c**2 * a**2
pe = a * b * c
Den = 3 * q1e - 2 * pe + 5
Nn = 15 * q1e**2 + 7 * q2e - 20 * pe * q1e + 9 * pe**2 + 71 * q1e - 108 * pe + 30
detG = 1 + 2 * pe - q1e
R = sp.Rational(20, 77) * Nn / Den**2


def main():
    print("=== MAXIMUM: verification of each step ===\n")

    Phi = sp.expand(463 * Den**2 - 240 * Nn)
    Phi_s = sp.expand(463 * (3 * Q1 - 2 * P + 5) ** 2
                      - 240 * (15 * Q1**2 + 7 * Q2 - 20 * P * Q1 + 9 * P**2 + 71 * Q1 - 108 * P + 30))
    print("S0  Phi in (Q1,Q2,P):", sp.expand(Phi_s))
    print("    matches Phi(a,b,c):",
          sp.simplify(sp.expand(Phi - Phi_s.subs({Q1: q1e, Q2: q2e, P: pe}))) == 0)

    G = sp.expand(Phi_s + 1680 * Q2)
    print("    coefficient of Q2 is", sp.diff(Phi_s, Q2), "; G =", G)

    Psi = sp.expand(G - 560 * Q1**2)
    print("\nS1  Psi = G - 560 Q1^2 =", Psi)
    print("    Phi - Psi = 1680*(Q1^2/3 - Q2) >= 0 by Maclaurin:",
          sp.simplify(Phi_s - Psi - 1680 * (Q1**2 / 3 - Q2)) == 0)
    # Maclaurin e2 <= e1^2/3 for nonneg x,y,z, with equality iff x=y=z
    x, y, z = sp.symbols("x y z", nonnegative=True)
    e1, e2 = x + y + z, x * y + y * z + z * x
    print("    e1^2/3 - e2 = (1/3)*",
          sp.factor(sp.expand(3 * (e1**2 / 3 - e2))),
          " = (1/3)*((x-y)^2+(y-z)^2+(z-x)^2)/2 :",
          sp.simplify(sp.expand(e1**2 / 3 - e2) - sp.expand(((x - y)**2 + (y - z)**2 + (z - x)**2) / 6)) == 0)

    print("\nS2  dPsi/dQ1 =", sp.diff(Psi, Q1))
    print("    on Q1 in [0,3], P in [-1/8,1]: max value =",
          max(sp.diff(Psi, Q1).subs({Q1: qq, P: pp})
              for qq in (0, 3) for pp in (sp.Rational(-1, 8), 1)))

    psi = sp.expand(Psi.subs(Q1, 1 + 2 * P))
    print("\nS3  psi(P) = Psi(1+2P,P) =", psi)
    print("    factored:", sp.factor(psi))
    print("    roots:", sp.solve(psi, P))

    print("\nS4  P >= -1/8 for three unit vectors:")
    expr = sp.expand(1 - 2 * sig**3 - 3 * sig**2)
    print("    1 - 2 sigma^3 - 3 sigma^2 =", sp.factor(expr), " (P = -sigma^3)")
    print("    so feasibility 1+2P >= Q1 >= 3|P|^(2/3) forces (2 sigma -1)(sigma+1)^2 <= 0 -> sigma <= 1/2")
    # numeric corroboration of P >= -1/8
    rng = np.random.default_rng(77)
    V = rng.normal(size=(600000, 3, 3))
    V /= np.linalg.norm(V, axis=2)[:, :, None]
    aa = np.einsum("ni,ni->n", V[:, 0], V[:, 1])
    bb = np.einsum("ni,ni->n", V[:, 0], V[:, 2])
    cc = np.einsum("ni,ni->n", V[:, 1], V[:, 2])
    print("    numeric: min abc over 6e5 random triples = %.10f   (-1/8 = %.10f)"
          % ((aa * bb * cc).min(), -1 / 8))

    print("\nEQUALITY analysis:")
    print("    Q2 = Q1^2/3  <=>  a^2 = b^2 = c^2")
    print("    Q1 = 1+2P    <=>  det Gram = 0 (coplanar axes)")
    print("    psi(P)=0     <=>  P = -1/8 (P = 11/2 impossible, |a|,|b|,|c| <= 1)")
    print("    => a^2=b^2=c^2=t^2, 3t^2 = 1+2(-1/8) = 3/4 => t = 1/2, abc = -1/8")
    sols = []
    for sa in (1, -1):
        for sb in (1, -1):
            for sc in (1, -1):
                if sa * sb * sc == -1:
                    sols.append((sp.Rational(sa, 2), sp.Rational(sb, 2), sp.Rational(sc, 2)))
    print("    sign patterns with abc = -1/8 :", sols)
    for s in sols:
        v = sp.simplify(R.subs({a: s[0], b: s[1], c: s[2]}))
        d = sp.simplify(detG.subs({a: s[0], b: s[1], c: s[2]}))
        print("      (a,b,c)=%s  R=%s  detG=%s" % (s, v, d))
    print("    all four are the SAME axis configuration: flipping the sign of one n_k")
    print("    flips exactly two of (a,b,c).  Gram = [[1,-1/2,-1/2],[-1/2,1,-1/2],[-1/2,-1/2,1]]")
    Gm = sp.Matrix([[1, sp.Rational(-1, 2), sp.Rational(-1, 2)],
                    [sp.Rational(-1, 2), 1, sp.Rational(-1, 2)],
                    [sp.Rational(-1, 2), sp.Rational(-1, 2), 1]])
    print("    det =", Gm.det(), " rank =", Gm.rank(), " eigenvalues =", Gm.eigenvals())

    # direct global numeric maximisation over axes (corroboration)
    print("\n  numeric corroboration: max over 6e5 random axis triples = %.15f"
          % float(sp.lambdify((a, b, c), R, "numpy")(aa, bb, cc).max()))
    print("  463/924 = %.15f" % (463 / 924))

    print("\n=== MINIMUM: equality analysis for p_6 = 1 ===")
    # Vandermonde identity used in the Cauchy-Schwarz step
    n = sp.symbols("n", integer=True)
    ok = all(sp.simplify(sum(sp.binomial(6, k) * sp.binomial(6, nn - k) for k in range(0, 7))
                         - sp.binomial(12, nn)) == 0 for nn in range(13))
    print("  Vandermonde sum_k C(6,k)C(6,n-k) = C(12,n) for n=0..12 :", ok)
    print("  n_J constants (s04b): n_0=13/7 > n_2=65/84 > n_4=13/154 > n_6=1/924, strict:",
          sp.Rational(13, 7) > sp.Rational(65, 84) > sp.Rational(13, 154) > sp.Rational(1, 924))
    print("  924*n_J =", [924 * x for x in (sp.Rational(13, 7), sp.Rational(65, 84),
                                            sp.Rational(13, 154), sp.Rational(1, 924))],
          " = C(13,6),C(13,4),C(13,2),C(13,0) :",
          [sp.binomial(13, 6), sp.binomial(13, 4), sp.binomial(13, 2), sp.binomial(13, 0)])


if __name__ == "__main__":
    main()
