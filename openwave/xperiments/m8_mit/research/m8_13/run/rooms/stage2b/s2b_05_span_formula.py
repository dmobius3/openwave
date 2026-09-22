"""s2b_05_span_formula.py -- r6 on span{|3,3>_n, |3,-3>_n}, cleanly.

s2b_04 printed r6(t) in a form sympy cluttered with Abs(); here it is done on the
real parametrisation u = cos(s)|3,3> + sin(s)|3,-3>, and independently from the
author's Step-1 identity, to pin the strictness the note's Step 4 needs.
"""
import numpy as np
import sympy as sp
from sympy.physics.quantum.cg import CG

ORDER = list(range(3, -4, -1)); IDX = {m: i for i, m in enumerate(ORDER)}
s = sp.Symbol("s", real=True)

def cg_str(m1, m2):
    M = m1 + m2
    if abs(M) > 6 or abs(m1) > 3 or abs(m2) > 3:
        return sp.Integer(0)
    f = sp.factorial
    return sp.sqrt(sp.Rational(f(6)*f(6)*f(6+M)*f(6-M),
                               f(12)*f(3+m1)*f(3-m1)*f(3+m2)*f(3-m2)))

cv = {m: sp.Integer(0) for m in range(-3,4)}
cv[3], cv[-3] = sp.cos(s), sp.sin(s)
th = {m: (-1)**(3-m)*sp.conjugate(cv[-m]) for m in range(-3,4)}
tot = 0
for Q in range(-6,7):
    e = sp.expand(sum(cg_str(m1,Q-m1)*cv[m1]*th[Q-m1] for m1 in range(-3,4) if abs(Q-m1)<=3))
    tot += sp.expand(sp.re(e)**2 + sp.im(e)**2)
r6 = sp.simplify(sp.expand(sp.trigsimp(tot)))
print("r6(s) on span{|3,3>,|3,-3>} =", sp.nsimplify(sp.simplify(r6)))
target = sp.Rational(463,924) - sp.cos(2*s)**2/2
print("  equals 463/924 - cos^2(2s)/2 :", sp.simplify(sp.expand(sp.trigsimp(r6 - target))) == 0)
for val in [0, sp.pi/8, sp.pi/6, sp.pi/4, sp.pi/3]:
    print("   s = %-10s r6 = %-14s = %.12f" % (val, sp.nsimplify(sp.simplify(r6.subs(s,val))),
                                               float(sp.N(r6.subs(s,val)))))
print("  so r6 < 463/924 STRICTLY unless cos(2s) = 0, i.e. |alpha| = |beta|.")
print("  and with |alpha|^2 = t: r6 = 463/924 - (2t-1)^2/2 ; at t=3/4: ",
      sp.nsimplify(sp.Rational(463,924) - sp.Rational(1,8)), "= 695/1848")

# TrNbar^2 is constant 171/2 on the whole span (the note's equality set for that bound)
def spin3():
    Jz = sp.diag(*[sp.Integer(m) for m in ORDER]); Jp = sp.zeros(7,7); Jm = sp.zeros(7,7)
    for m in ORDER:
        if m+1 <= 3: Jp[IDX[m+1], IDX[m]] = sp.sqrt(12-m*(m+1))
        if m-1 >= -3: Jm[IDX[m-1], IDX[m]] = sp.sqrt(12-m*(m-1))
    return (Jp+Jm)/2, (Jp-Jm)/(2*sp.I), Jz
FX,FY,FZ = spin3(); F=[FX,FY,FZ]
p1,p2 = sp.symbols("phi1 phi2", real=True)
cw = {m: sp.Integer(0) for m in range(-3,4)}
cw[3] = sp.cos(s)*sp.exp(sp.I*p1); cw[-3] = sp.sin(s)*sp.exp(sp.I*p2)
uu = sp.Matrix(7,1,[cw[m] for m in ORDER])
Nb = sp.Matrix(3,3, lambda i,k: sp.simplify(sp.re(sp.expand((uu.H*F[i]*F[k]*uu)[0,0]))))
print("  TrNbar^2 on the span, all phases:", sp.nsimplify(sp.simplify(sp.trigsimp(sp.expand((Nb*Nb).trace())))))
print("  |f|^2 on the span:", sp.nsimplify(sp.simplify(sp.trigsimp(sp.expand(
    sum(sp.re((uu.H*X*uu)[0,0])**2 + sp.im((uu.H*X*uu)[0,0])**2 for X in F))))))
