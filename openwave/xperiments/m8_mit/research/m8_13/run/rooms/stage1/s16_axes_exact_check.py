"""s16_axes_exact_check.py -- EXACT (not Monte-Carlo) confirmation of the
three-axes formula on asymmetric rational configurations.

s07 checked the formula exactly only at three symmetric configurations and s08
checked it to 1e-16 numerically.  Here the state is built from the axes in exact
arithmetic and rhat_6 is computed by the CG route, then compared with
R(a,b,c) = 20 N / (77 Den^2).
"""
import sympy as sp
from sympy.physics.quantum.cg import CG

ORDER = list(range(3, -4, -1)); IDX = {m: i for i, m in enumerate(ORDER)}
z1, z2 = sp.symbols("z1 z2")

def T6():
    out = []
    for Q in range(-6, 7):
        M = sp.zeros(7, 7)
        for i, mm in enumerate(ORDER):
            for k, m1 in enumerate(ORDER):
                if m1 - mm == Q:
                    M[i, k] = sp.nsimplify(sp.simplify(CG(3, m1, 3, -mm, 6, Q).doit())) * (-1)**(3+mm)
        out.append(M)
    return out
T6L = T6()

def state_from_axes(ns):
    P = sp.Integer(1)
    for n in ns:
        w = n[0] + sp.I*n[1]
        P *= (w/2)*z1**2 - n[2]*z1*z2 - (sp.conjugate(w)/2)*z2**2
    P = sp.expand(P)
    u = sp.zeros(7, 1)
    for k in range(7):
        ak = P.coeff(z1, k).coeff(z2, 6-k)
        u[IDX[k-3]] = ak*sp.sqrt(sp.factorial(k)*sp.factorial(6-k))
    nn = sp.sqrt(sp.simplify(sum(sp.Abs(x)**2 for x in u)))
    return sp.simplify(u/nn)

def rhat6(u):
    return sp.simplify(sum(sp.Abs(sp.simplify((u.H*M*u)[0,0]))**2 for M in T6L))

a,b,c = sp.symbols("a b c", real=True)
q1=a**2+b**2+c**2; q2=a**2*b**2+b**2*c**2+c**2*a**2; p=a*b*c
R = sp.Rational(20,77)*(15*q1**2+7*q2-20*p*q1+9*p**2+71*q1-108*p+30)/(3*q1-2*p+5)**2

cfgs = {
 "orthogonal (octahedron)": [(0,0,1),(1,0,0),(0,1,0)],
 "asym 1": [(0,0,1),(1,0,0),(sp.Rational(3,5),sp.Rational(4,5),0)],
 "asym 2": [(0,0,1),(sp.Rational(3,5),0,sp.Rational(4,5)),(sp.Rational(3,5),sp.Rational(4,5),0)],
 "two coincident": [(0,0,1),(0,0,1),(1,0,0)],
}
for name, ns in cfgs.items():
    ns = [tuple(sp.nsimplify(x) for x in n) for n in ns]
    u = state_from_axes(ns)
    direct = rhat6(u)
    g = (sum(ns[0][i]*ns[1][i] for i in range(3)),
         sum(ns[0][i]*ns[2][i] for i in range(3)),
         sum(ns[1][i]*ns[2][i] for i in range(3)))
    form = sp.simplify(R.subs({a:g[0], b:g[1], c:g[2]}))
    print("%-24s Gram=%s  direct=%s  formula=%s  equal=%s"
          % (name, g, sp.nsimplify(direct), sp.nsimplify(form),
             sp.simplify(direct-form)==0))
    print("     state c_m =", [sp.simplify(x) for x in u])
