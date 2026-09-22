"""s2a_02b_equality_set.py -- the equality set of step 3c, exhaustively and fast.

Same content as s2a_02_equality_set.py (which was killed: sympy's solve over the
radical system ran past the budget).  Here the surds are cleared first by the
rational substitution  a = sqrt(6) u ,  b = sqrt(6) v , which turns

   ellipse   (2/3)a^2 + 8b^2 = 1        into    4u^2 + 48v^2 = 1
   L0        15/sqrt(6) = 5 sqrt(6)/2   into    an exactly rational system

so every polynomial below has rational coefficients and the zero-dimensional
system is solved exactly.

For a 2x2 real symmetric M,   lambda_max(M) = L0   <=>   char_M(L0) = 0  AND  2 L0 >= tr M.
"""
import sympy as sp

u, v, L = sp.symbols("u v L", real=True)
S6 = sp.sqrt(6)
a_, b_ = S6*u, S6*v
L0 = 15/sp.sqrt(6)
ELL = sp.expand(4*u**2 + 48*v**2 - 1)

M1 = sp.Matrix([[5*a_, sp.sqrt(60)*b_], [sp.sqrt(60)*b_, -3*a_ + 12*b_]])
M2 = sp.Matrix([[5*a_, sp.sqrt(60)*b_], [sp.sqrt(60)*b_, -3*a_ - 12*b_]])
M3 = sp.Matrix([[0, sp.sqrt(240)*b_], [sp.sqrt(240)*b_, -4*a_]])
BLOCKS = {"M1": M1, "M2": M2, "M3": M3}

def lam_max(M):
    T, D = sp.simplify(M.trace()), sp.simplify(M.det())
    return sp.simplify(T/2 + sp.sqrt(sp.simplify(T**2/4 - D)))

def e_of(av, bv):
    return tuple(sp.radsimp(sp.simplify(x)) for x in (-av/3+2*bv, -av/3-2*bv, 2*av/3))

print("L0 = 15/sqrt(6) =", sp.nsimplify(sp.radsimp(L0)), "=", sp.N(L0, 25))
print("seventh eigenvalue is identically 0; 0 != L0 since L0 =", sp.N(L0,10), "> 0 -> no equality point there\n")

allpts = {}
for name, M in BLOCKS.items():
    cp = sp.expand(sp.radsimp(sp.simplify(M.charpoly(L).as_expr().subs(L, L0))))
    cp = sp.simplify(sp.nsimplify(cp))
    print("=== %s ===" % name)
    print("  det(M - L0 I) as a polynomial in (u,v):", sp.expand(cp))
    # eliminate v^2 using the ellipse where possible, then solve the system exactly
    sols = sp.solve([sp.expand(cp*6), ELL], [u, v], dict=True)
    pts = []
    for s in sols:
        uv, vv = sp.simplify(s[u]), sp.simplify(s[v])
        if sp.im(sp.N(uv)) != 0 or sp.im(sp.N(vv)) != 0:
            continue
        Tr = sp.simplify(M.trace().subs({u: uv, v: vv}))
        larger = sp.N(2*L0 - Tr) >= -sp.Float("1e-25")
        dd = sp.simplify(sp.radsimp(lam_max(M).subs({u: uv, v: vv}) - L0))
        pts.append((uv, vv, larger, dd))
    print("  real intersections with the ellipse: %d" % len(pts))
    for uv, vv, larger, dd in pts:
        print("     (u,v)=(%s,%s)  ->  (a,b)=(%s,%s)   2L0>=trM: %s   lam_max-L0 = %s"
              % (uv, vv, sp.radsimp(S6*uv), sp.radsimp(S6*vv), larger, sp.simplify(dd)))
    keep = [(uv, vv) for (uv, vv, larger, dd) in pts if larger and sp.simplify(dd) == 0]
    allpts[name] = keep
    print("  EQUALITY POINTS of %s (in (u,v)): %s\n" % (name, keep))

union = []
for name in allpts:
    for p in allpts[name]:
        if not any(sp.simplify(p[0]-q[0])==0 and sp.simplify(p[1]-q[1])==0 for q in union):
            union.append(p)
perms = {"(2,-1,-1)/sqrt6": (2,-1,-1), "(-1,2,-1)/sqrt6": (-1,2,-1), "(-1,-1,2)/sqrt6": (-1,-1,2)}
print("=== UNION over the three blocks:", len(union), "distinct points ===")
for (uv, vv) in union:
    av, bv = sp.radsimp(S6*uv), sp.radsimp(S6*vv)
    ev = e_of(av, bv)
    lab = [k for k,t in perms.items()
           if all(sp.simplify(ev[i] - sp.Rational(t[i],1)/S6)==0 for i in range(3))]
    print("  (a,b)=(%s,%s)  e=%s  -> %s" % (av, bv, ev, lab[0] if lab else "NOT a permutation"))

print("\n=== which block attains L0 where ===")
print("  %-26s %-5s %-5s %-5s" % ("(u,v)", "M1", "M2", "M3"))
for (uv, vv) in union:
    row = ["yes" if any(sp.simplify(uv-x)==0 and sp.simplify(vv-y)==0 for (x,y) in allpts[n]) else "no"
           for n in ("M1","M2","M3")]
    print("  %-26s %-5s %-5s %-5s" % ("(%s, %s)"%(uv,vv), *row))

print("\n=== the two loci NAMED in the text, intersected with the ellipse ===")
named = [("b = 0 with a > 0", v, lambda U,V: sp.N(U)>0),
         ("a = -2b with b > 0", u + 2*v, lambda U,V: sp.N(V)>0)]
count = 0
for label, cond, sgn in named:
    for s in sp.solve([cond, ELL], [u, v], dict=True):
        U, V = sp.simplify(s[u]), sp.simplify(s[v])
        if not sgn(U, V):
            continue
        count += 1
        av, bv = sp.radsimp(S6*U), sp.radsimp(S6*V)
        print("  %-20s -> (u,v)=(%s,%s)  e=%s" % (label, U, V, e_of(av, bv)))
print("  points delivered by the two named loci:", count)

print("\n=== the locus a = 2b with b < 0 (NOT named in the text) ===")
for s in sp.solve([u - 2*v, ELL], [u, v], dict=True):
    U, V = sp.simplify(s[u]), sp.simplify(s[v])
    if sp.N(V) >= 0:
        continue
    av, bv = sp.radsimp(S6*U), sp.radsimp(S6*V)
    print("  (u,v)=(%s,%s)  e=%s" % (U, V, e_of(av, bv)))

print("\n=== lambda_max at the DISCARDED signs (why a>0 and b>0 are needed) ===")
for label, U, V in [("b=0, a<0", sp.Rational(-1,2), sp.Integer(0)),
                    ("a=-2b, b<0", sp.Rational(1,4), sp.Rational(-1,8))]:
    for name, M in BLOCKS.items():
        lm = sp.radsimp(sp.simplify(lam_max(M).subs({u:U, v:V})))
        print("  %-12s %s: lam_max = %-18s = %-18s (L0 = %s)"
              % (label, name, sp.nsimplify(lm), sp.N(lm,14), sp.N(L0,14)))
