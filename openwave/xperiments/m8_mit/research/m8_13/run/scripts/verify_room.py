"""Maintainer's check of the stage-1 return, by the maintainer's own code.

The room's route is not the author's and not mine, so nothing here is a
re-run of its scripts: each claim is recomputed from the worklist's own
definitions in spin3.py.
"""
import sympy as sp
from sympy.physics.quantum.cg import CG
import spin3 as S

P, F = [], []
def ck(name, cond, detail=""):
    (P if cond else F).append(name)
    print(f"{'PASS' if cond else 'FAIL'}  {name}" + (f"   {detail}" if detail else ""))

# --- the FULL multipole spectrum, r_K for K = 0..6 -------------------------
def cg_tab(K):
    t = {}
    for Q in range(-K, K + 1):
        for m1 in S.MS:
            m2 = Q - m1
            if m2 in S.IDX:
                t[(m1, m2, Q)] = sp.nsimplify(sp.simplify(CG(3, m1, 3, m2, K, Q).doit()))
    return t

TAB = {K: cg_tab(K) for K in range(7)}

def rK(c, K, holo=False):
    """K-th multipole weight.  holo=False uses Theta u (the worklist's rho_K);
    holo=True uses u itself (the room's p_J on Sym^2)."""
    tu = c if holo else S.theta(c)
    tot = sp.Integer(0)
    for Q in range(-K, K + 1):
        s = sp.Integer(0)
        for m1 in S.MS:
            m2 = Q - m1
            if m2 in S.IDX:
                s += TAB[K][(m1, m2, Q)] * c[S.IDX[m1]] * tu[S.IDX[m2]]
        s = sp.expand(s)
        tot += sp.expand(sp.conjugate(s) * s)
    return sp.simplify(sp.expand(tot) / sp.expand(S.norm2(c) ** 2))

HEX = S.evec(p3=1, m3=1)
V3 = S.evec(p3=1)

spec_hex = [rK(HEX, K) for K in range(7)]
spec_v3 = [rK(V3, K) for K in range(7)]
ck("1 r_K(hexagon) = (1/7, 0, 25/84, 0, 9/154, 0, 463/924)",
   spec_hex == [sp.Rational(1,7), 0, sp.Rational(25,84), 0, sp.Rational(9,154), 0,
                sp.Rational(463,924)], str([str(x) for x in spec_hex]))
ck("2 r_K(v3) = (1/7, 9/28, 25/84, 1/6, 9/154, 1/84, 1/924)",
   spec_v3 == [sp.Rational(1,7), sp.Rational(9,28), sp.Rational(25,84), sp.Rational(1,6),
               sp.Rational(9,154), sp.Rational(1,84), sp.Rational(1,924)],
   str([str(x) for x in spec_v3]))
ck("3 both spectra sum to 1", sum(spec_hex) == 1 and sum(spec_v3) == 1)

# --- identically: sum_K r_K = 1 and r_0 = 1/7, as polynomials --------------
x = sp.symbols("x0:7", real=True); y = sp.symbols("y0:7", real=True)
cs = [x[i] + sp.I*y[i] for i in range(7)]
n2 = sp.expand(sum(x[i]**2 + y[i]**2 for i in range(7)))

def rK_num(c, K, holo=False):
    tu = c if holo else S.theta(c)
    tot = sp.Integer(0)
    for Q in range(-K, K + 1):
        s = sp.Integer(0)
        for m1 in S.MS:
            m2 = Q - m1
            if m2 in S.IDX:
                s += TAB[K][(m1, m2, Q)] * c[S.IDX[m1]] * tu[S.IDX[m2]]
        s = sp.expand(s)
        tot += sp.expand(sp.re(sp.expand(sp.conjugate(s)*s)))
    return sp.expand(tot)

tot_all = sp.expand(sum(rK_num(cs, K) for K in range(7)))
ck("4 sum_K ||rho_K||^2 = ||u||^4 IDENTICALLY (14 real variables)",
   sp.expand(tot_all - n2**2) == 0)
ck("5 r_0 = 1/7 identically, so the crude bound r_6 <= 6/7 is real",
   sp.expand(rK_num(cs, 0) - n2**2/7) == 0)

# --- the room's engine: r_6 = sum_J n_J p_J, p_J from the HOLOMORPHIC square
NJ = {0: sp.Rational(13,7), 2: sp.Rational(65,84), 4: sp.Rational(13,154),
      6: sp.Rational(1,924)}
lhs = rK_num(cs, 6)                                   # ||rho_6||^2, uses Theta
rhs = sp.expand(sum(NJ[J]*rK_num(cs, J, holo=True) for J in (0,2,4,6)))
ck("6 r_6 = (13/7)p_0 + (65/84)p_2 + (13/154)p_4 + (1/924)p_6, IDENTICALLY",
   sp.expand(lhs - rhs) == 0)
bad = sp.expand(sum((NJ[J] + (sp.Rational(1,1000) if J==2 else 0))*rK_num(cs, J, holo=True)
                    for J in (0,2,4,6)))
ck("6-arm perturbing n_2 by 1/1000 breaks it", sp.expand(lhs - bad) != 0)
ck("7 n_J = C(13, 6-J)/924 exactly",
   [NJ[J] for J in (0,2,4,6)] ==
   [sp.Rational(sp.binomial(13, 6-J), 924) for J in (0,2,4,6)],
   str([str(sp.binomial(13,6-J)) for J in (0,2,4,6)]))
ck("8 sum_J p_J = 1 identically (J even only; odd J vanish on Sym^2)",
   sp.expand(sum(rK_num(cs, J, holo=True) for J in (0,2,4,6)) - n2**2) == 0
   and all(sp.expand(rK_num(cs, J, holo=True)) == 0 for J in (1,3,5)))

# --- rho_6 of the hexagon, entry by entry ---------------------------------
r6h = S.rho6(HEX)
nz = {Q: sp.simplify(v) for Q, v in r6h.items() if sp.simplify(v) != 0}
ck("9 rho_6(hexagon) is non-zero only at Q = -6, 0, +6",
   set(nz) == {-6, 0, 6}, str({k: str(v) for k, v in nz.items()}))
# rho_6 is bilinear in (c, Theta u), so it scales with ||u||^2: my HEX above is
# UNNORMALIZED (c = 1, 1), which is why the raw entries came out 2x.  Normalize.
HEXN = S.evec(p3=1/sp.sqrt(2), m3=1/sp.sqrt(2))
nzn = {Q: sp.simplify(v) for Q, v in S.rho6(HEXN).items() if sp.simplify(v) != 0}
ck("10 on the UNIT hexagon, entries are 1/2 at Q = +-6 and 1/sqrt(924) at Q = 0",
   sp.simplify(nzn.get(6) - sp.Rational(1,2)) == 0
   and sp.simplify(nzn.get(-6) - sp.Rational(1,2)) == 0
   and sp.simplify(nzn.get(0) - 1/sp.sqrt(924)) == 0,
   str({k: str(v) for k, v in nzn.items()}))
ck("11 so ||rho_6||^2 = 1/4 + 1/4 + 1/924 = 463/924",
   sp.simplify(sum(v**2 for v in nzn.values()) - sp.Rational(463,924)) == 0)

# --- the two other constellations the room names --------------------------
ck("12 r_6(v0) = 100/231, the (1,1,1) constellation", S.r6(S.evec(p0=1)) == sp.Rational(100,231))
ck("13 r_6((v1 - v-1)/sqrt2) = 145/308, the (1,0,0) constellation",
   S.r6(S.evec(p1=1, m1=-1)) == sp.Rational(145,308),
   str(S.r6(S.evec(p1=1, m1=-1))))

# --- the room's final algebra, independently --------------------------------
q1, p_, q2 = sp.symbols("q1 p q2", real=True)
Den = 3*q1 - 2*p_ + 5
N_  = 15*q1**2 + 7*q2 - 20*p_*q1 + 9*p_**2 + 71*q1 - 108*p_ + 30
Phi = sp.expand(463*Den**2 - 240*N_)
G   = 567*q1**2 - 308*p_**2 - 756*p_*q1 + 16660*p_ - 3150*q1 + 4375
ck("14 Phi = 463 Den^2 - 240 N equals -1680 q2 + G, exactly",
   sp.expand(Phi - (-1680*q2 + G)) == 0)
Psi = 7*q1**2 - 308*p_**2 - 756*p_*q1 + 16660*p_ - 3150*q1 + 4375
ck("15 Psi = G with q1^2 reduced from 567 to 7, i.e. G - 1680*(q1^2/3) ",
   sp.expand(G - 1680*q1**2/3 - Psi) == 0)
psi = sp.expand(Psi.subs(q1, 1 + 2*p_))
ck("16 psi(p) = Psi(1+2p, p) = -1792 p^2 + 9632 p + 1232 = 1792(p+1/8)(11/2-p)",
   sp.expand(psi - 1792*(p_ + sp.Rational(1,8))*(sp.Rational(11,2) - p_)) == 0,
   f"psi = {sp.factor(psi)}")
ck("17 dPsi/dq1 = 14 q1 - 756 p - 3150 < 0 on q1 in [0,3], p in [-1/8,1]",
   sp.diff(Psi, q1).subs({q1: 3, p_: -sp.Rational(1,8)}) < 0,
   f"max value {sp.diff(Psi,q1).subs({q1:3,p_:-sp.Rational(1,8)})}")
sig = sp.Symbol("sigma", positive=True)
ck("18 abc >= -1/8: 1 - 2s^3 - 3s^2 = -(s+1)^2(2s-1), so s <= 1/2",
   sp.expand(1 - 2*sig**3 - 3*sig**2 + (sig+1)**2*(2*sig-1)) == 0)
ck("19 Maclaurin e2 <= e1^2/3 as a sum of squares",
   sp.expand((sp.Symbol('X')+sp.Symbol('Y')+sp.Symbol('Z'))**2/3
             - (sp.Symbol('X')*sp.Symbol('Y')+sp.Symbol('Y')*sp.Symbol('Z')
                +sp.Symbol('Z')*sp.Symbol('X'))
             - ((sp.Symbol('X')-sp.Symbol('Y'))**2+(sp.Symbol('Y')-sp.Symbol('Z'))**2
                +(sp.Symbol('Z')-sp.Symbol('X'))**2)/6) == 0)
ck("20 the room's S2 counterexample: Psi(q1=3, p=-1/8) < 0, so the Gram bound is load-bearing",
   Psi.subs({q1: 3, p_: -sp.Rational(1,8)}) < 0,
   f"Psi = {sp.nsimplify(Psi.subs({q1:3,p_:-sp.Rational(1,8)}))} = {float(Psi.subs({q1:3,p_:-sp.Rational(1,8)})):.1f}")
# a^2=b^2=c^2=t^2 gives q1 = 3t^2; with q1 = 1+2p and p = abc = +-t^3, psi(p)=0
# forces p = -1/8, hence t = 1/2.
sol = sp.solve([sp.Eq(3*sp.Symbol('t2'), 1 + 2*p_), sp.Eq(psi, 0)],
               [sp.Symbol('t2'), p_], dict=True)
real = [d for d in sol if d[p_].is_real and 0 <= d[sp.Symbol('t2')] <= 1]
ck("21 q2 = q1^2/3 with q1 = 1+2p and psi(p) = 0 forces t^2 = 1/4 and p = -1/8",
   len(real) == 1 and real[0][p_] == sp.Rational(-1,8)
   and real[0][sp.Symbol('t2')] == sp.Rational(1,4),
   str([{str(k): str(v) for k, v in d.items()} for d in sol]))
Gram = sp.Matrix([[1,-sp.Rational(1,2),-sp.Rational(1,2)],
                  [-sp.Rational(1,2),1,-sp.Rational(1,2)],
                  [-sp.Rational(1,2),-sp.Rational(1,2),1]])
ck("22 that Gram has eigenvalues 3/2, 3/2, 0: rank 2, three coplanar axes at 60 degrees",
   sorted(Gram.eigenvals().keys()) == [0, sp.Rational(3,2)] and Gram.det() == 0)

# --- W2: the three-axes formula, attacked at points the room did not use ----
# r_6 = (20/77) N / Den^2 on Theta-invariant states.  The room verified the
# spinor identity behind it at seven configurations.  Here the formula is
# tested against r_6 values computed independently from the worklist's own
# definition, at three constellations whose states are known exactly.
print()
print("--- W2: the three-axes formula against independently computed values ---")

def axes_formula(a, b, c):
    q1 = a**2 + b**2 + c**2
    q2 = a**2*b**2 + b**2*c**2 + c**2*a**2
    pp = a*b*c
    Dn = 3*q1 - 2*pp + 5
    Nn = 15*q1**2 + 7*q2 - 20*pp*q1 + 9*pp**2 + 71*q1 - 108*pp + 30
    return sp.simplify(sp.Rational(20,77) * Nn / Dn**2)

CASES = [
    ("hexagon, 3 coplanar axes at 60 deg", (sp.Rational(-1,2),)*3, HEXN, sp.Rational(463,924)),
    ("v0, three coincident axes (1,1,1)", (sp.Integer(1),)*3, S.evec(p0=1), sp.Rational(100,231)),
    ("(v1 - v-1)/sqrt2, axes (1,0,0)", (sp.Integer(1), sp.Integer(0), sp.Integer(0)),
     S.evec(p1=1, m1=-1), sp.Rational(145,308)),
]
allok = True
for nm, gram, st, known in CASES:
    f = axes_formula(*gram)
    direct = S.r6(st)
    ok = (f == known == direct)
    allok &= ok
    print(f"  {'PASS' if ok else 'FAIL'}  {nm}: formula {f}, direct {direct}, expected {known}")
ck("23 the three-axes formula reproduces three independently computed exact values",
   allok)
ck("24 and it is discriminating: the three values differ",
   len({axes_formula(*g) for _, g, _, _ in CASES}) == 3)

# --- W1 is closed by check 6 above -----------------------------------------
print()
print("--- W1: the room flagged off-block vanishing of N as resting on Schur")
print("    plus a 1e-16 numerical check.  Check 6 above proves the consequence")
print("    (r_6 = sum_J n_J p_J) as an EXACT polynomial identity in 14 real")
print("    variables, which does not go through N at all.  Link closed.")

print()
print(f"FINAL {len(P)} passed, {len(F)} failed")
if F: print("FAILED: " + "; ".join(F))
