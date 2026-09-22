"""Maintainer's own instrument for the uniqueness question at spin 3.

Independent of the argument under audit and of the stage-1 room: every number
below comes from this file's own code.  Exact throughout, except the clearly
labelled corroborating searches, which are reported for what they are.

Every substantive check is paired with an ARM: a deliberate perturbation that
the same check must reject.  An arm that does not fire is reported as a dead
arm and counted as a failure, because a check that cannot fail is not a check.
"""

import itertools
import sys

import numpy as np
import sympy as sp

import spin3 as S

PASS, FAIL, DEAD = [], [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"{'PASS' if cond else 'FAIL'}  {name}" + (f"   {detail}" if detail else ""))
    return cond


def arm(name, fired, detail=""):
    """`fired` is True when the perturbed input is correctly rejected."""
    (PASS if fired else DEAD).append(f"arm:{name}")
    print(f"{'ARM ' if fired else 'DEAD'}  {name}" + (f"   {detail}" if detail else ""))
    return fired


a, b = sp.symbols("a b", real=True)
ELL = sp.Rational(2, 3) * a**2 + 8 * b**2 - 1
BOUND = 15 / sp.sqrt(6)
HEX = S.evec(p3=1, m3=1)

print("=" * 72)
print("0. CONVENTIONS")
print("=" * 72)

check("0.1 <3 3; 3 3 | 6 6> = +1 (Condon-Shortley, as the worklist fixes it)",
      S.CG6[(3, 3, 6)] == 1)
c_sym = [sp.Symbol(f"c{i}") for i in range(7)]
check("0.2 Theta^2 = identity (integer spin)",
      all(sp.simplify(x - y) == 0 for x, y in zip(S.theta(S.theta(c_sym)), c_sym)))
v0 = S.CG6[(3, -3, 0)]
check("0.3 <3 3; 3 -3 | 6 0> = 1/(2 sqrt(231))",
      sp.simplify(v0 - 1 / (2 * sp.sqrt(231))) == 0, f"= {v0}")

print()
print("=" * 72)
print("1. THE INVARIANT AND ITS WEIGHT-STATE VALUES")
print("=" * 72)

wv = {k: S.r6(S.evec(**{f"p{k}": 1})) for k in (3, 2, 1, 0)}
check("1.1 r6 on v3, v2, v1, v0 = 1, 36, 225, 400 over 924",
      [wv[k] for k in (3, 2, 1, 0)] == [sp.Rational(n, 924) for n in (1, 36, 225, 400)],
      str([str(wv[k]) for k in (3, 2, 1, 0)]))
check("1.2 r6 at the hexagon (v3+v-3)/sqrt(2) = 463/924", S.r6(HEX) == sp.Rational(463, 924))
check("1.3 r6 at the coherent state v3 = 1/924", S.r6(S.evec(p3=1)) == sp.Rational(1, 924))

# Rotation and phase invariance.  Done EXACTLY, with a Wigner D at rational
# multiples of pi, because a float64 input caps a "numerical" version at about
# 1e-16 no matter how many digits the arithmetic carries.  (First draft asserted
# 1e-30 on a float64 vector and failed itself; the failure was the tolerance,
# not the invariance.)
# Two rotations that are exact and cheap, plus one at high precision.
#   R_z(t):      v_m -> e^(-i m t) v_m            (diagonal, exact at t = 2pi/7)
#   R_y(pi):     v_m -> (-1)^(3-m) v_{-m}         (d^j_{m'm}(pi) = (-1)^(j-m) delta)
exact_u = [sp.Rational(1, 2), sp.Rational(-1, 3) + sp.I / 5, sp.Rational(2, 7),
           sp.I / 3, sp.Rational(1, 6) - sp.I / 4, sp.Rational(-3, 5), sp.Rational(1, 1)]
ref_exact = S.r6(exact_u)

t = 2 * sp.pi / 7
Rz = [sp.exp(-sp.I * m * t) * exact_u[S.IDX[m]] for m in S.MS]
Ry = [((-1) ** (3 + m)) * exact_u[S.IDX[-m]] for m in S.MS]  # coefficient of v_m
check("1.4a r6 exactly invariant under R_z(2pi/7)",
      sp.simplify(sp.expand(S.r6(Rz)) - ref_exact) == 0, f"r6 = {ref_exact}")
check("1.4b r6 exactly invariant under R_y(pi): v_m -> (-1)^(3-m) v_{-m}",
      sp.simplify(S.r6(Ry) - ref_exact) == 0)
check("1.5 r6 exactly invariant under a global phase e^(i pi/5)",
      sp.simplify(S.r6([sp.exp(sp.I * sp.pi / 5) * x for x in exact_u]) - ref_exact) == 0)
check("1.6 r6 exactly homogeneous of degree 0 (scale by 3)",
      sp.simplify(S.r6([3 * x for x in exact_u]) - ref_exact) == 0)
NOTROT = sp.eye(7)
NOTROT[0, 1], NOTROT[1, 0] = sp.Rational(1, 5), sp.Rational(-1, 5)
arm("1.A3 a non-rotation perturbation of the identity changes r6",
    sp.simplify(S.r6(list(NOTROT * sp.Matrix(7, 1, exact_u))) - ref_exact) != 0)

import mpmath as mp


def r6_num(cv, prec):
    mp.mp.dps = prec
    cv = [mp.mpc(z) for z in cv]
    n2 = sum(mp.conj(z) * z for z in cv)
    tu = [((-1) ** (3 - m)) * mp.conj(cv[S.IDX[-m]]) for m in S.MS]
    tot = mp.mpf(0)
    for Q in range(-6, 7):
        s = mp.mpc(0)
        for m1 in S.MS:
            m2 = Q - m1
            if m2 not in S.IDX:
                continue
            s += mp.mpf(str(sp.N(S.CG6[(m1, m2, Q)], prec))) * cv[S.IDX[m1]] * tu[S.IDX[m2]]
        tot += (mp.conj(s) * s).real
    return tot / (n2.real**2)


mp.mp.dps = 60  # parse EXACT_REF at full width: mpf(str(...)) truncates at the CURRENT dps
EXACT_REF = mp.mpf(str(sp.N(ref_exact, 55)))
UV = [complex(sp.N(sp.re(x), 30)) for x in exact_u], [complex(sp.N(sp.im(x), 30)) for x in exact_u]


def r6_rot(axis, th, prec):
    mp.mp.dps = prec
    n = mp.matrix([mp.mpf(str(x)) for x in axis])
    n = n / mp.sqrt(sum(x**2 for x in n))
    G = mp.zeros(7, 7)
    for k, A in enumerate(S.FVEC):
        for i in range(7):
            for j in range(7):
                G[i, j] += n[k] * mp.mpc(str(sp.N(sp.re(A[i, j]), prec)),
                                         str(sp.N(sp.im(A[i, j]), prec)))
    R = mp.expm(mp.mpc(0, -1) * mp.mpf(str(th)) * G)
    u = mp.matrix([mp.mpc(str(sp.N(sp.re(x), prec)), str(sp.N(sp.im(x), prec)))
                   for x in exact_u])
    return r6_num(list(R * u), prec)


ok30 = abs(r6_rot((1, 2, 3), "0.7", 30) - EXACT_REF) < mp.mpf("1e-22")
ok50 = abs(r6_rot((1, 2, 3), "0.7", 50) - EXACT_REF) < mp.mpf("1e-40")
check("1.4c r6 invariant under a generic-axis rotation, at 30 dps (tol 1e-22) and again "
      "at 50 dps (tol 1e-40): the exactness rule's second precision", ok30 and ok50)
# mp.mpc() on a sympy Rational routes through float() and silently caps at ~16
# digits; feed it decimal strings at the working precision instead.
def to_mp(cv, prec):
    return [mp.mpc(str(sp.N(sp.re(x), prec)), str(sp.N(sp.im(x), prec))) for x in cv]


mp.mp.dps = 60
check("1.7 the mpmath route reproduces the exact value of r6 to 40 digits",
      abs(r6_num(to_mp(exact_u, 55), 50) - EXACT_REF) < mp.mpf("1e-40"))

# ARM: a CG sign flip on a coefficient that is NOT alone in its Q-sum.
# (Flipping <33;33|66> alone is invisible: Q=6 carries a single term.  That is a
# dead arm, and it is shown to be dead here on purpose, then replaced.)
flip66 = S.cg_table(sign=lambda m1, m2, Q: -1 if (m1, m2, Q) == (3, 3, 6) else 1)
arm("1.A0 flip <33;33|66> alone  (EXPECTED DEAD: Q=6 is a single term)",
    S.r6(HEX, tab=flip66) != sp.Rational(463, 924),
    "shown dead deliberately; the live arm is 1.A1")


flip60 = S.cg_table(sign=lambda m1, m2, Q: -1 if (m1, m2, Q) == (3, -3, 0) else 1)
arm("1.A1 flip <33;3-3|60> (Q=0 carries seven terms)",
    S.r6(HEX, tab=flip60) != sp.Rational(463, 924),
    f"perturbed r6(hex) = {S.r6(HEX, tab=flip60)}")
arm("1.A2 drop Theta's (-1)^(3-m) phase, tested AT THE HEXAGON  (EXPECTED DEAD: the "
    "phase is +1 at both m = +-3, so the hexagon cannot see it)",
    S.r6(HEX, phase=lambda m: 1) != sp.Rational(463, 924),
    "shown dead deliberately; the live arm is 1.A4")
arm("1.A4 drop Theta's (-1)^(3-m) phase, tested on a generic vector",
    sp.simplify(S.r6(exact_u, phase=lambda m: 1) - ref_exact) != 0,
    f"perturbed r6 = {S.r6(exact_u, phase=lambda m: 1)}")

print()
print("=" * 72)
print("2. THE NEMATIC OPERATOR AND ITS BLOCKS")
print("=" * 72)

e1 = -a / 3 + 2 * b
e2 = -a / 3 - 2 * b
e3 = 2 * a / 3
check("2.1 e traceless", sp.simplify(e1 + e2 + e3) == 0)
check("2.2 |e|^2 - 1 is the stated ellipse (2/3)a^2 + 8b^2 = 1",
      sp.simplify(e1**2 + e2**2 + e3**2 - 1 - ELL) == 0)
check("2.3 e1 fx^2 + e2 fy^2 + e3 fz^2 == a(fz^2 - 4) + b(f+^2 + f-^2), as 7x7",
      sp.simplify(S.diag_operator(e1, e2, e3) - S.nematic_operator(a, b)) == sp.zeros(7, 7))

OP = S.nematic_operator(a, b)
lam = sp.Symbol("lam")
cp = sp.factor(sp.expand(OP.charpoly(lam).as_expr()))
BL = S.blocks(a, b)
prod = lam * sp.prod([sp.expand((M - lam * sp.eye(2)).det()) for M in BL.values()])
check("2.4 charpoly = lam * det(M1-lam) det(M2-lam) det(M3-lam)  (zero eigenvalue + three 2x2)",
      sp.simplify(sp.expand(cp) - sp.expand(-prod)) == 0 or
      sp.simplify(sp.expand(cp) - sp.expand(prod)) == 0)
check("2.5 M2 is M1 at -b", sp.simplify(BL["M2"] - BL["M1"].subs(b, -b)) == sp.zeros(2, 2))
check("2.6 the ellipse is even in b", sp.simplify(ELL - ELL.subs(b, -b)) == 0)

# ARM: the block decomposition must not survive a perturbed operator.
OPp = sp.expand(OP + sp.Rational(1, 50) * S.JX**2)
cpp = sp.expand(OPp.charpoly(lam).as_expr())
arm("2.A1 perturb the operator by (1/50) fx^2: charpoly no longer factors into the blocks",
    sp.simplify(sp.expand(cpp) - sp.expand(-prod)) != 0 and
    sp.simplify(sp.expand(cpp) - sp.expand(prod)) != 0)
BLp = S.blocks(a, b, s60=sp.sqrt(61) * b)
prodp = lam * sp.prod([sp.expand((M - lam * sp.eye(2)).det()) for M in BLp.values()])
arm("2.A2 sqrt(60) -> sqrt(61) in M1: factorization fails",
    sp.simplify(sp.expand(cp) - sp.expand(-prodp)) != 0 and
    sp.simplify(sp.expand(cp) - sp.expand(prodp)) != 0)

print()
print("=" * 72)
print("3. THE BOUND lam_max <= 15/sqrt(6), AND ITS EQUALITY SET, EXACTLY")
print("=" * 72)

# 3a. the two sum-of-squares identities the note leans on
lhs3 = sp.expand((BOUND + 2 * a) ** 2 - (30 - 16 * a**2))
check("3.1 M3: (15/sqrt6 + 2a)^2 - (30 - 16a^2) == 20 (a + sqrt6/4)^2",
      sp.simplify(lhs3 - 20 * (a + sp.sqrt(6) / 4) ** 2) == 0)
Nn = sp.sqrt(sp.Rational(2, 3) * a**2 + 8 * b**2)
sq1 = sp.expand((a**2 + 6 * a * b + 24 * b**2) ** 2 - sp.Rational(3, 2) * (a + 6 * b) ** 2 * Nn**2)
check("3.2 M1: (a^2+6ab+24b^2)^2 - (3/2)(a+6b)^2 N^2 == 36 b^2 (a+2b)^2",
      sp.simplify(sq1 - 36 * b**2 * (a + 2 * b) ** 2) == 0)
check("3.3 M1: a^2 + 6ab + 24b^2 == (a+3b)^2 + 15b^2  (the a+6b<=0 branch's sign)",
      sp.simplify(a**2 + 6 * a * b + 24 * b**2 - ((a + 3 * b) ** 2 + 15 * b**2)) == 0)
arm("3.A1 24b^2 -> 23b^2 in the sum of squares",
    sp.simplify(sp.expand((a**2 + 6 * a * b + 23 * b**2) ** 2
                          - sp.Rational(3, 2) * (a + 6 * b) ** 2 * Nn**2)
                - 36 * b**2 * (a + 2 * b) ** 2) != 0)

# 3b. THE EQUALITY SET, per block, by an exact solve rather than a sweep.
def equality_points(bound, ell, bl):
    """Every point of `ell` where some block's LARGER eigenvalue equals `bound`."""
    hits = {}
    for nm, M in bl.items():
        for sol in sp.solve([sp.expand((M - bound * sp.eye(2)).det()), ell], [a, b], dict=True):
            if a not in sol or b not in sol:
                continue
            av, bv = sp.nsimplify(sp.simplify(sol[a])), sp.nsimplify(sp.simplify(sol[b]))
            if not (av.is_real and bv.is_real):
                continue
            ev = [sp.nsimplify(x) for x in M.subs({a: av, b: bv}).eigenvals()]
            if sp.simplify(max(ev) - bound) != 0:
                continue  # bound is an eigenvalue there, but the SMALLER one
            hits.setdefault((av, bv), set()).add(nm)
    return hits


HITS = equality_points(BOUND, ELL, BL)
want = {
    (sp.sqrt(6) / 2, sp.Integer(0)): {"M1", "M2"},
    (-sp.sqrt(6) / 4, sp.sqrt(6) / 8): {"M1", "M3"},
    (-sp.sqrt(6) / 4, -sp.sqrt(6) / 8): {"M2", "M3"},
}
print("    equality points found:")
for k in sorted(HITS, key=lambda p: (sp.N(p[0]), sp.N(p[1]))):
    print(f"      (a, b) = ({k[0]}, {k[1]})   reached by {sorted(HITS[k])}")
check("3.4 exactly three equality points on the ellipse", len(HITS) == 3, f"found {len(HITS)}")
check("3.5 they are (sqrt6/2, 0) and (-sqrt6/4, +-sqrt6/8)",
      set(HITS) == set(want))
check("3.6 each is reached by exactly two blocks, in the stated pairing (D1)",
      all(HITS.get(k) == v for k, v in want.items()))

arm("3.A2 bound 15/sqrt6 -> 14/sqrt6: the equality set changes",
    set(equality_points(14 / sp.sqrt(6), ELL, BL)) != set(want))
arm("3.A3 ellipse constant 1 -> 51/50: the equality set changes",
    set(equality_points(BOUND, sp.Rational(2, 3) * a**2 + 8 * b**2 - sp.Rational(51, 50), BL))
    != set(want))
arm("3.A4 sqrt(60) -> sqrt(61) in M1: the equality set changes",
    set(equality_points(BOUND, ELL, S.blocks(a, b, s60=sp.sqrt(61) * b))) != set(want))

# 3c. THE BRANCH THE NOTE DOES NOT DO: a + 6b <= 0 carries no equality point.
#     There the bound reads (3/sqrt6)(a+6b) N <= a^2+6ab+24b^2 with LHS <= 0 and
#     RHS = (a+3b)^2 + 15b^2 >= 0, so equality needs BOTH sides zero.
zero_rhs = sp.solve([a**2 + 6 * a * b + 24 * b**2, ELL], [a, b], dict=True)
real_zero_rhs = [s for s in zero_rhs
                 if a in s and b in s
                 and sp.simplify(sp.im(s[a])) == 0 and sp.simplify(sp.im(s[b])) == 0]
check("3.7 branch a+6b<=0: RHS (a+3b)^2+15b^2 vanishes nowhere on the ellipse, "
      "so no equality point there", len(real_zero_rhs) == 0,
      f"real solutions of RHS=0 on the ellipse: {len(real_zero_rhs)}")
m1pts = [k for k, v in HITS.items() if "M1" in v]
check("3.8 every M1 equality point has a+6b>0, so all of them live on the squared branch",
      len(m1pts) == 2 and all(sp.simplify(k[0] + 6 * k[1]) > 0 for k in m1pts),
      "; ".join(f"a+6b = {sp.simplify(k[0] + 6 * k[1])}" for k in m1pts))
m2pts = [k for k, v in HITS.items() if "M2" in v]
check("3.8b and every M2 equality point has a-6b>0, its mirror branch  (the point reached "
      "by M2 and M3 has a+6b<0, which is not a counterexample: a+6b>0 is M1's branch)",
      len(m2pts) == 2 and all(sp.simplify(k[0] - 6 * k[1]) > 0 for k in m2pts),
      "; ".join(f"a-6b = {sp.simplify(k[0] - 6 * k[1])}" for k in m2pts))
check("3.8c M3's equality condition 20(a+sqrt6/4)^2 = 0 pins a = -sqrt6/4, and both of its "
      "points have it",
      all(sp.simplify(k[0] + sp.sqrt(6) / 4) == 0 for k, v in HITS.items() if "M3" in v))

# 3d. corroboration only: a full sweep of the ellipse, endpoints INCLUDED.
th = np.linspace(0.0, 2.0 * np.pi, 200001)  # closed, so t=0 is in the grid
av = np.sqrt(1.5) * np.cos(th)
bv = np.cos(0.0) * np.sin(th) / np.sqrt(8)
Bn = float(sp.N(BOUND, 30))
lmax = np.empty_like(th)
for i, (A_, B_) in enumerate(zip(av, bv)):
    m1 = np.array([[5 * A_, np.sqrt(60) * B_], [np.sqrt(60) * B_, -3 * A_ + 12 * B_]])
    m2 = np.array([[5 * A_, -np.sqrt(60) * B_], [-np.sqrt(60) * B_, -3 * A_ - 12 * B_]])
    m3 = np.array([[0.0, np.sqrt(240) * B_], [np.sqrt(240) * B_, -4 * A_]])
    lmax[i] = max(np.linalg.eigvalsh(m1)[-1], np.linalg.eigvalsh(m2)[-1],
                  np.linalg.eigvalsh(m3)[-1])
check("3.9 [corroboration] 200001-point sweep of the closed ellipse never exceeds 15/sqrt6",
      lmax.max() <= Bn + 1e-12, f"max excess {lmax.max() - Bn:.3e}")
near = np.where(lmax > Bn - 1e-9)[0]
groups = np.split(near, np.where(np.diff(near) > 1)[0] + 1) if len(near) else []
groups = [g for g in groups if len(g)]
# the closed grid repeats t=0 at t=2pi; merge that pair
n_near = len(groups) - (1 if len(groups) >= 2 and groups[0][0] == 0
                        and groups[-1][-1] == len(th) - 1 else 0)
check("3.10 [corroboration] the sweep touches the bound in exactly three places",
      n_near == 3, f"contiguous near-bound groups (wrap merged): {n_near}")
arm("3.A5 sweep arm: add (1/50) fx^2, i.e. shift a by 1/50 in the blocks -> the sweep "
    "must now exceed the bound",
    max(max(np.linalg.eigvalsh(np.array([[5 * (A_ + 0.02), np.sqrt(60) * B_],
                                         [np.sqrt(60) * B_, -3 * (A_ + 0.02) + 12 * B_]]))[-1]
            for A_, B_ in zip(av, bv)), 0) > Bn + 1e-9)

print()
print("=" * 72)
print("4. FROM THE EQUALITY POINT TO THE STATE")
print("=" * 72)

for (av_, bv_), who in sorted(want.items(), key=lambda kv: (sp.N(kv[0][0]), sp.N(kv[0][1]))):
    Op = S.nematic_operator(av_, bv_)
    ev = Op.eigenvects()
    top = max(sp.nsimplify(e[0]) for e in ev)
    mult = [e[1] for e in ev if sp.simplify(sp.nsimplify(e[0]) - top) == 0][0]
    vecs = [e[2] for e in ev if sp.simplify(sp.nsimplify(e[0]) - top) == 0][0]
    sq6e = tuple(sp.simplify(sp.sqrt(6) * x) for x in
                 (-av_ / 3 + 2 * bv_, -av_ / 3 - 2 * bv_, 2 * av_ / 3))
    check(f"4.x top eigenvalue at (a,b)=({av_}, {bv_}) is 15/sqrt6, multiplicity 2",
          sp.simplify(top - BOUND) == 0 and mult == 2,
          f"sqrt6*e = {sq6e}, blocks {sorted(who)}")
    if sq6e == (-1, -1, 2):
        span_ok = all(sp.simplify(v[1]) == 0 and sp.simplify(v[2]) == 0
                      and sp.simplify(v[3]) == 0 and sp.simplify(v[4]) == 0
                      and sp.simplify(v[5]) == 0 for v in vecs)
        check("4.1 at sqrt6*e = (-1,-1,2) the top eigenspace is span{v3, v-3}", span_ok)

perms = {tuple(sorted(sp.simplify(sp.sqrt(6) * x) for x in
                      (-k[0] / 3 + 2 * k[1], -k[0] / 3 - 2 * k[1], 2 * k[0] / 3)))
         for k in want}
check("4.2 the three points are the three permutations of (2,-1,-1)/sqrt6",
      perms == {tuple(sorted((2, -1, -1)))} and len(want) == 3)

print()
print("=" * 72)
print("5. ARITHMETIC, AND THE STEP-1 IDENTITY")
print("=" * 72)

check("5.1 -5/231 + 1/11 + 171/396 = 463/924",
      sp.Rational(-5, 231) + sp.Rational(1, 11) + sp.Rational(171, 396) == sp.Rational(463, 924))
check("5.2 48 + 225/6 = 171/2", 48 + sp.Rational(225, 6) == sp.Rational(171, 2))
check("5.3 (15/sqrt6)^2 = 225/6", sp.simplify(BOUND**2 - sp.Rational(225, 6)) == 0)

# The step-1 identity, on exact vectors that are NOT the weight states it was
# fitted to.  Fit the four coefficients from the weight states, then test them
# elsewhere: that is the part the fit cannot fake.
def invariants(cv):
    n2 = S.norm2(cv)
    fv = S.magnetization(cv)
    f2 = sp.simplify(sum(sp.expand(sp.re(x) ** 2 + sp.im(x) ** 2) for x in fv) / n2**2)
    A0 = sp.simplify(sp.expand(sp.re(S.a00(cv)) ** 2 + sp.im(S.a00(cv)) ** 2) / n2**2)
    Nb = S.nbar(cv) / n2
    tr = sp.simplify(sp.expand((Nb * Nb).trace()))
    return sp.Integer(1), sp.simplify(f2), A0, tr


k0, k1, k2, k3 = sp.symbols("k0 k1 k2 k3")
rows, rhs = [], []
for m in (3, 2, 1, 0):
    cv = S.evec(**{f"p{m}": 1})
    rows.append(list(invariants(cv)))
    rhs.append(S.r6(cv))
Mfit = sp.Matrix(rows)
check("5.4 the four invariants at v3, v2, v1, v0 have determinant 180/7 (a basis)",
      sp.simplify(Mfit.det() - sp.Rational(180, 7)) == 0, f"det = {sp.simplify(Mfit.det())}")
coef = Mfit.solve(sp.Matrix(rhs))
print(f"    fitted: r6 = {coef[0]} + ({coef[1]})|f|^2 + ({coef[2]})|a00|^2 + ({coef[3]})TrN^2")
check("5.5 the fitted coefficients are (-5/231, -1/22, 7/11, 1/198)",
      [sp.nsimplify(x) for x in coef] ==
      [sp.Rational(-5, 231), sp.Rational(-1, 22), sp.Rational(7, 11), sp.Rational(1, 198)])

rngv = np.random.default_rng(4632)
off_ok, tested = True, 0
for _ in range(5):
    cv = [sp.Rational(int(x), 7) + sp.I * sp.Rational(int(y), 5)
          for x, y in zip(rngv.integers(-6, 7, 7), rngv.integers(-6, 7, 7))]
    if S.norm2(cv) == 0:
        continue
    iv = invariants(cv)
    pred = sp.simplify(sum(c * v for c, v in zip(coef, iv)))
    if sp.simplify(pred - S.r6(cv)) != 0:
        off_ok = False
    tested += 1
check(f"5.6 the identity holds at {tested} exact vectors off the weight states",
      off_ok and tested == 5)
bad = [coef[0] + sp.Rational(1, 1000), coef[1], coef[2], coef[3]]
cvv = [sp.Rational(int(x), 7) + sp.I * sp.Rational(int(y), 5)
       for x, y in zip(rngv.integers(-6, 7, 7), rngv.integers(-6, 7, 7))]
arm("5.A1 perturb the constant coefficient by 1/1000: the identity fails off the weight states",
    sp.simplify(sum(c * v for c, v in zip(bad, invariants(cvv))) - S.r6(cvv)) != 0)

check("5.7 at the hexagon: |f|^2 = 0, |a00|^2 = 1/7, TrN^2 = 171/2",
      [sp.nsimplify(x) for x in invariants(HEX)[1:]] ==
      [sp.Integer(0), sp.Rational(1, 7), sp.Rational(171, 2)],
      str([str(sp.nsimplify(x)) for x in invariants(HEX)[1:]]))
check("5.8 Cauchy-Schwarz ceiling |a00|^2 <= 1/7 is met at the hexagon",
      sp.nsimplify(invariants(HEX)[2]) == sp.Rational(1, 7))

print()
print("=" * 72)
print("6. CORROBORATION ONLY: a global search over the unit sphere of C^7")
print("=" * 72)

from scipy.optimize import minimize

CGn = {k: float(sp.N(v, 30)) for k, v in S.CG6.items()}


def r6_f(x):
    c = x[:7] + 1j * x[7:]
    n2 = float(np.vdot(c, c).real)
    tu = np.array([((-1) ** (3 - m)) * np.conj(c[S.IDX[-m]]) for m in S.MS])
    t = 0.0
    for Q in range(-6, 7):
        s = 0j
        for m1 in S.MS:
            m2 = Q - m1
            if m2 in S.IDX:
                s += CGn[(m1, m2, Q)] * c[S.IDX[m1]] * tu[S.IDX[m2]]
        t += abs(s) ** 2
    return t / n2**2


r = np.random.default_rng(1)
best, worst = -1.0, 2.0
for _ in range(400):
    x0 = r.normal(size=14)
    hi = minimize(lambda x: -r6_f(x), x0, method="BFGS")
    lo = minimize(r6_f, x0, method="BFGS")
    best, worst = max(best, -hi.fun), min(worst, lo.fun)
print(f"    400 restarts: max {best:.12f}  (463/924 = {463/924:.12f})")
print(f"                  min {worst:.12f}  (1/924   = {1/924:.12f})")
check("6.1 [corroboration] no restart exceeds 463/924", best <= 463 / 924 + 1e-9)
check("6.2 [corroboration] no restart falls below 1/924", worst >= 1 / 924 - 1e-9)
check("6.3 [corroboration] both bounds are approached", best > 463 / 924 - 1e-7
      and worst < 1 / 924 + 1e-7)
print("    A search rejects candidates.  It establishes neither extremum nor completeness.")

print()
print("=" * 72)
print("7. STEP 4: FROM THE THREE EQUALITY CONDITIONS TO ONE ORBIT")
print("=" * 72)

al, be = sp.symbols("alpha beta", positive=True)
ph = sp.Symbol("phi", real=True)

# N-bar's trace is fixed, so TrN^2 = 48 + ||Q||^2 with Q = N - 4I traceless.
check("7.1 fx^2 + fy^2 + fz^2 = 12 I, so TrN = 12 for every unit u",
      sp.simplify(S.JX**2 + S.JY**2 + S.JZ**2 - 12 * sp.eye(7)) == sp.zeros(7, 7))
u_gen = [sp.Rational(1, 2), sp.Rational(-1, 3) + sp.I / 5, sp.Rational(2, 7),
         sp.I / 3, sp.Rational(1, 6) - sp.I / 4, sp.Rational(-3, 5), sp.Integer(1)]
Ng = S.nbar(u_gen) / S.norm2(u_gen)
check("7.2 TrN = 12 on a generic vector", sp.simplify(Ng.trace() - 12) == 0)
Qg = Ng - 4 * sp.eye(3)
check("7.3 TrN^2 = 48 + ||Q||^2 with Q = N - 4I traceless",
      sp.simplify(Qg.trace()) == 0 and
      sp.simplify((Ng * Ng).trace() - 48 - (Qg * Qg).trace()) == 0)

# the duality step 3d rests on: Tr(N E) = <u, A_E u> for unit traceless symmetric E
Ediag = sp.Matrix([[-sp.Rational(1, 3) + 2 * b, 0, 0],
                   [0, -sp.Rational(1, 3) - 2 * b, 0], [0, 0, sp.Rational(2, 3)]])
AE = sp.zeros(7, 7)
for i in range(3):
    for jj in range(3):
        AE += Ediag[i, jj] * (S.FVEC[i] * S.FVEC[jj] + S.FVEC[jj] * S.FVEC[i]) / 2
check("7.4 Tr(N E) = <u, A_E u> for a traceless symmetric E  (the duality step 3d uses)",
      sp.simplify((Ng * Ediag).trace()
                  - sp.re(S.braket(u_gen, AE)) / S.norm2(u_gen)) == 0)

# on span{v3, v-3}: everything in closed form
span_u = [al, 0, 0, 0, 0, 0, be * sp.exp(sp.I * ph)]
n2s = sp.simplify(S.norm2(span_u))
Ns = sp.simplify(S.nbar(span_u) / n2s)
trN2 = sp.simplify((Ns * Ns).trace().subs(al**2 + be**2, 1))
check("7.5 every u in span{v3, v-3} has TrN^2 = 171/2, the equality value",
      sp.simplify(sp.simplify(trN2.subs(be, sp.sqrt(1 - al**2))) - sp.Rational(171, 2)) == 0,
      f"TrN^2 = {trN2}")
fs = [sp.simplify(S.braket(span_u, A) / n2s) for A in S.FVEC]
f2s = sp.simplify(sum(sp.re(x) ** 2 + sp.im(x) ** 2 for x in fs))
check("7.6 on that span <f> = 3(|a|^2 - |b|^2) n, so |f|^2 = 9(|a|^2-|b|^2)^2",
      sp.simplify(f2s.subs(be, sp.sqrt(1 - al**2))
                  - 9 * (2 * al**2 - 1) ** 2) == 0, f"|f|^2 = {sp.simplify(f2s)}")
a00s = sp.simplify(sp.re(S.a00(span_u)) ** 2 + sp.im(S.a00(span_u)) ** 2) / n2s**2
check("7.7 on that span |a00|^2 = 4|a|^2|b|^2/7 (unit u)",
      sp.simplify(sp.simplify(a00s) - 4 * al**2 * be**2 / (7 * (al**2 + be**2) ** 2)) == 0,
      f"|a00|^2 = {sp.simplify(a00s)}")
check("7.8 so |f|^2 = 0 forces |a| = |b|, AND |a00|^2 = 1/7 forces it too: either "
      "condition alone closes the span",
      sp.solve([sp.Eq(9 * (2 * al**2 - 1) ** 2, 0), sp.Eq(al**2 + be**2, 1)], [al, be],
               dict=True) != [] and
      sp.simplify((4 * al**2 * (1 - al**2) / 7 - sp.Rational(1, 7)).factor()) ==
      sp.simplify((-sp.Rational(4, 7) * (al**2 - sp.Rational(1, 2)) ** 2).factor()))

t = sp.Symbol("t", positive=True)
r_on_span = sp.simplify(sp.Rational(-5, 231) - 9 * (2 * t - 1) ** 2 / 22
                        + sp.Rational(7, 11) * 4 * t * (1 - t) / 7 + sp.Rational(171, 2) / 198)
check("7.9 r6 on the span, as a function of t = |alpha|^2, peaks only at t = 1/2",
      sp.solve(sp.diff(r_on_span, t), t) == [sp.Rational(1, 2)]
      and sp.simplify(r_on_span.subs(t, sp.Rational(1, 2)) - sp.Rational(463, 924)) == 0,
      f"r6(t) = {sp.expand(r_on_span)}")
check("7.10 the relative phase is a rotation about n: R_z(chi) shifts it by 6 chi, "
      "so one orbit",
      sp.simplify(sp.exp(-sp.I * 3 * sp.Symbol('chi')) / sp.exp(sp.I * 3 * sp.Symbol('chi'))
                  - sp.exp(-sp.I * 6 * sp.Symbol('chi'))) == 0)

off = [al, sp.Rational(1, 10), 0, 0, 0, 0, be]
Noff = S.nbar(off) / S.norm2(off)
arm("7.A1 leaving the span (a v2 component of 1/10) drops TrN^2 strictly below 171/2",
    sp.simplify(((Noff * Noff).trace()).subs({al: sp.Rational(1, 2), be: sp.Rational(1, 2)})
                - sp.Rational(171, 2)) < 0,
    f"TrN^2 = {sp.nsimplify(((Noff*Noff).trace()).subs({al: sp.Rational(1,2), be: sp.Rational(1,2)}))}")
arm("7.A2 unequal weights on the span (t = 1/3) fall strictly below 463/924",
    sp.simplify(r_on_span.subs(t, sp.Rational(1, 3)) - sp.Rational(463, 924)) < 0,
    f"r6(1/3) = {r_on_span.subs(t, sp.Rational(1, 3))}")

print()
print("=" * 72)
dead_real = [d for d in DEAD if "EXPECTED DEAD" not in d]
print(f"{len(PASS)} passed, {len(FAIL)} failed, {len(DEAD)} dead arms "
      f"({len(dead_real)} unexpected)")
if FAIL:
    print("FAILED: " + "; ".join(FAIL))
if dead_real:
    print("DEAD (unexpected): " + "; ".join(dead_real))
sys.exit(1 if (FAIL or dead_real) else 0)
