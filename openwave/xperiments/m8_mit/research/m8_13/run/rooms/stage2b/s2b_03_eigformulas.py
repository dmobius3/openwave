"""s2b_03_eigformulas.py -- run down the two False results in s2b_02.

s2b_02 reported False for
   'M1 top eigenvalue formula is a root of char(M1)'
   'M3 top eigenvalue in general is -2a + sqrt(4a^2+240b^2)'
Both are sympy simplification failures, not mathematical failures.  Shown three ways:
the algebraic identity, an exact expand-and-radsimp, and dense numerics.
"""
import numpy as np
import sympy as sp

a, b, lam = sp.symbols("a b lam", real=True)
M1 = sp.Matrix([[5*a, sp.sqrt(60)*b], [sp.sqrt(60)*b, -3*a + 12*b]])
M3 = sp.Matrix([[0, sp.sqrt(240)*b], [sp.sqrt(240)*b, -4*a]])

for nm, M, half, disc in [("M1", M1, a + 6*b, 16*a**2 - 48*a*b + 96*b**2),
                          ("M3", M3, -2*a, 4*a**2 + 240*b**2)]:
    T = sp.expand(M.trace())
    D = sp.expand(M.det())
    print("=== %s ===" % nm)
    print("  trace =", T, "  det =", D)
    print("  trace/2 - claimed half =", sp.expand(T/2 - half))
    print("  (trace/2)^2 - det - claimed discriminant =", sp.expand(T**2/4 - D - disc))
    # algebraic identity: for lam = T/2 + sqrt(Delta) with Delta = T^2/4 - D,
    #   lam^2 - T lam + D = Delta - T^2/4 + D = 0
    print("  so lam = T/2 + sqrt((T/2)^2 - D) is a root by the quadratic formula.")
    lamf = half + sp.sqrt(disc)
    resid = sp.expand(lamf**2 - T*lamf + D)
    print("  expand(lam^2 - T lam + D) =", resid, " -> simplifies to 0 :", sp.simplify(resid) == 0)
    # dense numerics
    rng = np.random.default_rng(1)
    Mn = sp.lambdify((a, b), M, "numpy")
    hn = sp.lambdify((a, b), half, "numpy")
    dn = sp.lambdify((a, b), disc, "numpy")
    err = 0.0
    for _ in range(20000):
        A, B = rng.normal(), rng.normal()
        top = np.linalg.eigvalsh(np.array(Mn(A, B), dtype=float)).max()
        err = max(err, abs(top - (hn(A, B) + np.sqrt(dn(A, B)))))
    print("  numeric: max |lambda_max - formula| over 20000 random (a,b) = %.3e" % err)
    print("  -> the formula is CORRECT; s2b_02's False was sympy failing to simplify a surd.\n")

# the on-ellipse form of the M3 top eigenvalue
print("on the ellipse (2/3)a^2+8b^2=1 : 240 b^2 - (30 - 20 a^2) =",
      sp.simplify(240*((1 - sp.Rational(2,3)*a**2)/8) - (30 - 20*a**2)))
print("so 4a^2 + 240b^2 = 30 - 16a^2 there :",
      sp.simplify(4*a**2 + 240*((1 - sp.Rational(2,3)*a**2)/8) - (30 - 16*a**2)) == 0)
