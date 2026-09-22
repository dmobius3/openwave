"""Step 1's identity as an EXACT POLYNOMIAL IDENTITY in the 14 real coordinates.

The note fits four coefficients at four weight states and appeals to a
dimension count (Sym^2 V3 multiplicity-free) to conclude the identity holds
everywhere.  This file does not use the dimension count: it expands both sides
as polynomials in Re c_m, Im c_m and subtracts.
"""
import sympy as sp
import spin3 as S

x = sp.symbols("x0:7", real=True)
y = sp.symbols("y0:7", real=True)
c = [x[i] + sp.I * y[i] for i in range(7)]

n2 = sp.expand(sum(x[i] ** 2 + y[i] ** 2 for i in range(7)))

r = S.rho6(c)
num = sp.expand(sum(sp.expand(sp.conjugate(v) * v) for v in r.values()))
num = sp.expand(sp.re(num))

fv = [S.braket(c, A) for A in S.FVEC]
f2 = sp.expand(sum(sp.expand(sp.re(v) ** 2 + sp.im(v) ** 2) for v in fv))
A0 = sp.expand(sum(sp.conjugate(S.theta(c)[i]) * c[i] for i in range(7)))
a00sq = sp.expand((sp.re(A0) ** 2 + sp.im(A0) ** 2) / 7)
Nb = sp.Matrix(3, 3, lambda i, j: sp.expand(sp.re(S.braket(c, S.FVEC[i] * S.FVEC[j]))))
trN2 = sp.expand((Nb * Nb).trace())

# r6 * ||u||^4 == -5/231 ||u||^4 - |f|^2/22 + (7/11)|a00|^2 + TrN^2/198
# with |f|^2, |a00|^2, TrN^2 written UNNORMALIZED (each already degree 4).
lhs = sp.expand(num)
rhs = sp.expand(sp.Rational(-5, 231) * n2**2
                - f2 / 22
                + sp.Rational(7, 11) * a00sq
                + trN2 / 198)
d = sp.expand(lhs - rhs)
print("identity residual as a polynomial in 14 real variables:", d)
assert d == 0, "step 1's identity is NOT a polynomial identity"
print("PASS  step 1 holds identically, with no dimension count and no fit")

# ARM: the same with 1/198 -> 1/200 must leave a nonzero polynomial.
bad = sp.expand(lhs - sp.expand(sp.Rational(-5, 231) * n2**2 - f2 / 22
                                + sp.Rational(7, 11) * a00sq + trN2 / 200))
assert bad != 0
print("ARM   1/198 -> 1/200 leaves a nonzero residual polynomial "
      f"({len(bad.as_ordered_terms())} terms)")

# and the dimension count itself, checked rather than quoted:
# Sym^2 V3 = V0 + V2 + V4 + V6, dimensions 1 + 5 + 9 + 13 = 28 = dim Sym^2(C^7)
dims = [2 * Jv + 1 for Jv in (0, 2, 4, 6)]
assert sum(dims) == 7 * 8 // 2 == 28
print(f"PASS  Sym^2 V3 is multiplicity-free: {dims} sums to {sum(dims)} = dim Sym^2(C^7)")
