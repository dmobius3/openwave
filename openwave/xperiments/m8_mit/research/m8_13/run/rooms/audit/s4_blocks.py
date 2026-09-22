"""C4, C5, C6, C7: the operator A(e) = e1 fx^2 + e2 fy^2 + e3 fz^2 on the ellipse.

Everything here is exact sympy; the only numerics is a corroborating scan.
"""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import numpy as np
import sympy as sp
import core as C

HERE = pathlib.Path(__file__).parent
out = []
P = out.append

a, b = sp.symbols('a b', real=True)
Jx, Jy, Jz = C.JOPS_EXACT
e1 = -a / 3 + 2 * b
e2 = -a / 3 - 2 * b
e3 = 2 * a / 3
P("C5 parametrization check (from the claim's own words)")
P(f"  e3 = 2a/3, e1-e2 = 4b, trace e1+e2+e3 = {sp.simplify(e1+e2+e3)}")
P(f"  e1-e2 = {sp.simplify(e1-e2)}")
P(f"  |e|^2 = {sp.expand(e1**2+e2**2+e3**2)}   (claim says (2/3)a^2 + 8b^2)")
P(f"  |e|^2 - (2a^2/3 + 8b^2) = {sp.simplify(e1**2+e2**2+e3**2 - (sp.Rational(2,3)*a**2+8*b**2))}")

A = sp.expand(e1 * Jx * Jx + e2 * Jy * Jy + e3 * Jz * Jz)
A = sp.Matrix(7, 7, lambda i, j: sp.nsimplify(sp.simplify(sp.expand(A[i, j]))))
P("")
P("A(e) in the m-basis (rows/cols m = 3,2,1,0,-1,-2,-3):")
for i in range(7):
    P("   " + "  ".join(f"{sp.simplify(A[i,j])}" for j in range(7)))

# --- block basis: P: v_m -> v_{-m} symmetric / antisymmetric ---------------
r2 = 1 / sp.sqrt(2)


def vec(d):
    v = sp.zeros(7, 1)
    for m, x in d.items():
        v[C.IDX[m]] = x
    return v


bas = {
    "sym(v2,v-2)": vec({2: r2, -2: r2}), "v0": vec({0: 1}),
    "sym(v3,v-3)": vec({3: r2, -3: r2}), "sym(v1,v-1)": vec({1: r2, -1: r2}),
    "asym(v3,v-3)": vec({3: r2, -3: -r2}), "asym(v1,v-1)": vec({1: r2, -1: -r2}),
    "asym(v2,v-2)": vec({2: r2, -2: -r2}),
}
names = list(bas)
Bm = sp.Matrix.hstack(*[bas[n] for n in names])
Ab = sp.simplify(Bm.T * A * Bm)
P("")
P("A in the parity-adapted basis " + str(names) + ":")
for i in range(7):
    P("   " + "  ".join(f"{sp.simplify(Ab[i,j])}" for j in range(7)))
P(f"  is block diagonal (0-2 | 2-4 | 4-6 | 6): "
  f"{all(sp.simplify(Ab[i,j])==0 for i in range(7) for j in range(7) if (i//2 != j//2 or i==6 or j==6) and not (i==6 and j==6) and not (i//2==j//2 and i<6 and j<6))}")
P(f"  entry for asym(v2,v-2) (claim: a zero eigenvalue) = {sp.simplify(Ab[6,6])}")
P(f"  its coupling to everything else = "
  f"{[sp.simplify(Ab[6,j]) for j in range(7)]}")

M3 = sp.simplify(Ab[0:2, 0:2])
M1 = sp.simplify(Ab[2:4, 2:4])
M2 = sp.simplify(Ab[4:6, 4:6])
P("")
P(f"block on (sym(v2,v-2), v0)   = {M3.tolist()}   claim M3 = [[0,sqrt240 b],[sqrt240 b,-4a]]")
P(f"   difference from claim      = {sp.simplify(M3 - sp.Matrix([[0, sp.sqrt(240)*b],[sp.sqrt(240)*b, -4*a]]))}")
P(f"block on (sym(v3,v-3), sym(v1,v-1)) = {M1.tolist()}   claim M1 = [[5a,sqrt60 b],[sqrt60 b,-3a+12b]]")
P(f"   difference from claim      = {sp.simplify(M1 - sp.Matrix([[5*a, sp.sqrt(60)*b],[sp.sqrt(60)*b, -3*a+12*b]]))}")
P(f"block on (asym3, asym1)       = {M2.tolist()}")
M1_mb = M1.subs(b, -b)
P(f"claim: M2 = M1 at -b -> M1(-b) = {sp.simplify(M1_mb).tolist()}")
P(f"   M2 - M1(-b) (as matrices)  = {sp.simplify(M2 - M1_mb).tolist()}")
P(f"   charpolys equal?           = "
  f"{sp.simplify(sp.expand(M2.charpoly(sp.Symbol('x')).as_expr() - M1_mb.charpoly(sp.Symbol('x')).as_expr()))==0}")
P("   (sign of the off-diagonal is a basis-orientation choice; it does not move the spectrum)")

# --- C6/C7: where does lambda_max reach 15/sqrt6 on the ellipse? -----------
lam = 15 / sp.sqrt(6)
P("")
P(f"lambda* = 15/sqrt6 = {sp.nsimplify(lam)} = {sp.N(lam,20)};  lambda*^2 = {sp.simplify(lam**2)}")
al, be = sp.symbols('alpha beta', real=True)
sub = {a: sp.sqrt(6) * al, b: sp.sqrt(6) * be}
ell = sp.simplify((sp.Rational(2, 3) * a**2 + 8 * b**2 - 1).subs(sub))   # 4al^2+48be^2-1
P(f"ellipse in (alpha,beta) with a=sqrt6*alpha, b=sqrt6*beta : {sp.expand(ell)} = 0")

blocks = {"M1": M1, "M2": M2, "M3": M3}
zeros = {}
for nm, Mb in blocks.items():
    Mb = sp.simplify(Mb.subs(sub))
    g = sp.expand(sp.simplify((lam - Mb[0, 0]) * (lam - Mb[1, 1]) - Mb[0, 1] * Mb[1, 0]))
    t = sp.expand(sp.simplify(2 * lam - Mb[0, 0] - Mb[1, 1]))
    # reduce g modulo the ellipse: replace beta^2 by (1-4 alpha^2)/48
    gred = sp.simplify(sp.expand(g).subs(be**2, (1 - 4 * al**2) / 48))
    gred = sp.simplify(sp.expand(gred.subs(be**2, (1 - 4 * al**2) / 48)))
    P("")
    P(f"{nm}:  det(lambda* I - M) on the ellipse = {sp.factor(sp.simplify(gred))}")
    P(f"      trace(lambda* I - M)             = {t}")
    zeros[nm] = sp.factor(sp.simplify(gred))

P("")
P("Sign analysis of each factor on the ellipse 4a^2+48b^2=1 (so |alpha|<=1/2, |beta|<=1/(4sqrt3)):")
P(f"  max of (u*alpha+v*beta) on the ellipse = sqrt(u^2/4 + v^2/48)  [Cauchy-Schwarz]")
for u_, v_, desc in [(-1, 6, "6beta - alpha  (appears in M1 factor alpha+1-6beta)"),
                     (-1, -6, "-6beta - alpha (appears in M2 factor alpha+1+6beta)"),
                     (2, 0, "2alpha         (appears in 1-2alpha)"),
                     (2, 12, "2alpha+12beta  (trace test M1)"),
                     (4, 0, "4alpha         (trace test M3)")]:
    mx = sp.sqrt(sp.Rational(u_, 1)**2 / 4 + sp.Rational(v_, 1)**2 / 48)
    P(f"  max({desc}) = {sp.nsimplify(mx)} = {sp.N(mx,12)}")

P("")
P("Exact solution of {g_k = 0} & ellipse:")
sols_all = set()
for nm in blocks:
    sols = sp.solve([zeros[nm], ell], [al, be], dict=True)
    P(f"  {nm}: {sols}")
    for s in sols:
        sols_all.add((sp.nsimplify(s[al]), sp.nsimplify(s[be])))
P(f"  union of contact points (alpha,beta): {sorted(sols_all, key=lambda t: (float(t[0]), float(t[1])))}")
P("  in (a,b) = sqrt6*(alpha,beta):")
for s in sorted(sols_all, key=lambda t: (float(t[0]), float(t[1]))):
    P(f"     a = {sp.nsimplify(sp.sqrt(6)*s[0])}, b = {sp.nsimplify(sp.sqrt(6)*s[1])}   "
      f"e = ({sp.nsimplify(sp.simplify(e1.subs({a:sp.sqrt(6)*s[0], b:sp.sqrt(6)*s[1]})))}, "
      f"{sp.nsimplify(sp.simplify(e2.subs({a:sp.sqrt(6)*s[0], b:sp.sqrt(6)*s[1]})))}, "
      f"{sp.nsimplify(sp.simplify(e3.subs({a:sp.sqrt(6)*s[0], b:sp.sqrt(6)*s[1]})))})"
      f"   sqrt6*e = ({sp.nsimplify(sp.simplify(sp.sqrt(6)*e1.subs({a:sp.sqrt(6)*s[0], b:sp.sqrt(6)*s[1]})))},"
      f"{sp.nsimplify(sp.simplify(sp.sqrt(6)*e2.subs({a:sp.sqrt(6)*s[0], b:sp.sqrt(6)*s[1]})))},"
      f"{sp.nsimplify(sp.simplify(sp.sqrt(6)*e3.subs({a:sp.sqrt(6)*s[0], b:sp.sqrt(6)*s[1]})))})")

# --- C7: at each contact point, full 7x7 spectrum --------------------------
P("")
P("C7: full spectrum of A at each contact point, and which blocks reach lambda*")
pts = [(sp.sqrt(6) / 2, 0), (-sp.sqrt(6) / 4, sp.sqrt(6) / 8), (-sp.sqrt(6) / 4, -sp.sqrt(6) / 8)]
for (av, bv) in pts:
    Ap = sp.Matrix(7, 7, lambda i, j: sp.nsimplify(sp.simplify(A[i, j].subs({a: av, b: bv}))))
    ev = Ap.eigenvals()
    evs = {sp.nsimplify(sp.simplify(k)): v for k, v in ev.items()}
    P(f"  (a,b) = ({sp.nsimplify(av)}, {sp.nsimplify(bv)})")
    P(f"    spectrum = { {str(sp.nsimplify(k)): v for k, v in evs.items()} }")
    top = max(evs, key=lambda k: float(sp.N(k)))
    P(f"    top = {sp.nsimplify(top)} = {sp.N(top,18)};  equals lambda*? "
      f"{sp.simplify(top-lam)==0};  multiplicity = {evs[top]}")
    for nm, Mb in blocks.items():
        Mp = sp.simplify(Mb.subs({a: av, b: bv}))
        bev = {sp.nsimplify(sp.simplify(k)): v for k, v in Mp.eigenvals().items()}
        hit = sum(v for k, v in bev.items() if sp.simplify(k - lam) == 0)
        P(f"    {nm} eigs = { {str(k): v for k,v in bev.items()} }  -> reaches lambda* with mult {hit}")
    z = sp.simplify(Ap * bas["asym(v2,v-2)"])
    P(f"    A * asym(v2,v-2) = {z.T.tolist()}  (claim: zero eigenvector)")

txt = "\n".join(map(str, out))
print(txt)
(HERE / "out_s4.txt").write_text(txt + "\n")
