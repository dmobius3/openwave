"""s15_affine_relations.py -- exact check of the affine relations among p_2,p_4,p_6
on the real form V_R, which s06 only identified NUMERICALLY (rank test).

Claimed (read off numerically in s06 from two sample points plus the hexagon):
    p_2 = (11/7) p_6 - 24/49 ,    p_4 = 66/49 - (18/7) p_6      (||u|| = 1)
Homogeneous form, with n = ||u||^2 :
    49 p_2 - 77 p_6 + 24 n^2 = 0 ,   49 p_4 + 126 p_6 - 66 n^2 = 0
These play NO role in the final proof; this is an independent cross-check.
"""
from pathlib import Path
import sympy as sp
from sympy.physics.quantum.cg import CG

ORDER = list(range(3, -4, -1)); IDX = {m: i for i, m in enumerate(ORDER)}

def cgx(m1, m2, Jt, Q):
    if m1 + m2 != Q or abs(m1) > 3 or abs(m2) > 3 or abs(Q) > Jt:
        return sp.Integer(0)
    return sp.nsimplify(sp.simplify(CG(3, m1, 3, m2, Jt, Q).doit()))

t = sp.symbols("t0:7", real=True)
cc = {0: sp.I * t[0]}
for m, (pp, qq) in zip([1, 2, 3], [(t[1], t[2]), (t[3], t[4]), (t[5], t[6])]):
    cc[m] = pp + sp.I * qq
    cc[-m] = (-1) ** (3 + m) * sp.conjugate(cc[m])
u = {m: sp.expand(cc[m]) for m in range(-3, 4)}
nrm = sp.expand(sum(sp.Abs(u[m]) ** 2 for m in range(-3, 4)))
nrm = sp.expand(sp.simplify(nrm))
print("||u||^2 =", nrm)

def pJ(J):
    tot = 0
    for Q in range(-J, J + 1):
        e = sp.expand(sum(cgx(m1, Q - m1, J, Q) * u[m1] * u[Q - m1]
                          for m1 in range(-3, 4) if abs(Q - m1) <= 3))
        tot += sp.expand(sp.re(e) ** 2 + sp.im(e) ** 2)
    return sp.expand(sp.simplify(tot))

P = {J: pJ(J) for J in (0, 2, 4, 6)}
print("p_0 - n^2/7 =", sp.simplify(sp.expand(P[0] - nrm**2 / 7)))
print("sum_J p_J - n^2 =", sp.simplify(sp.expand(sum(P.values()) - nrm**2)))
print("49 p_2 - 77 p_6 + 24 n^2 =", sp.simplify(sp.expand(49*P[2] - 77*P[6] + 24*nrm**2)))
print("49 p_4 + 126 p_6 - 66 n^2 =", sp.simplify(sp.expand(49*P[4] + 126*P[6] - 66*nrm**2)))
