"""Items 0 and 1: group order, perfectness, CG value, sector multiplicities at levels n <= 18.
Everything group-theoretic is DERIVED from the two generators (exact arithmetic in Q(sqrt5))."""
import json
import math
from fractions import Fraction as F
import mpmath
from common import (HERE, Q5, Q1, Q2, ONE, qmul, qconj, qnorm2, close_group, generate_gamma,
                    cg_exact, quat_to_su2, Dj_generic, q5num)

out = {}
# ---------------- item 0(a)
n1, n2 = qnorm2(Q1), qnorm2(Q2)
print("||q1||^2 =", n1, " ||q2||^2 =", n2)
out["norm2_q1"] = str(n1)
out["norm2_q2"] = str(n2)
G = generate_gamma()
order = len(G)
print("order of Gamma =", order)
out["order"] = order
# closure check that can fail: every product lies in G, every inverse lies in G
Gs = set(G)
assert all(qmul(a, b) in Gs for a in G for b in G), "not closed"
assert all(qconj(a) in Gs for a in G), "inverse missing"
assert all(qnorm2(a) == Q5(1) for a in G)
comms = {qmul(qmul(a, b), qmul(qconj(a), qconj(b))) for a in G for b in G}
Dsub = close_group(list(comms))
print("number of distinct commutators =", len(comms), " derived subgroup order =", len(Dsub))
out["n_commutators"] = len(comms)
out["derived_order"] = len(Dsub)
out["perfect"] = (len(Dsub) == order)

# conjugacy classes
classes = []
seen = set()
for x in G:
    if x in seen:
        continue
    cl = {qmul(qmul(g, x), qconj(g)) for g in G}
    seen |= cl
    classes.append(sorted(cl, key=lambda q: tuple((c.a, c.b) for c in q)))
classes.sort(key=lambda c: (len(c), -float(c[0][0].a) - float(c[0][0].b) * 5 ** 0.5))
print("number of classes =", len(classes), " sizes =", [len(c) for c in classes])
reps = [c[0] for c in classes]
sizes = [len(c) for c in classes]
traces = [r[0] for r in reps]   # Re(q) = cos(theta); SU(2) trace = 2 Re q
# element orders
def order_of(q):
    k, y = 1, q
    while y != ONE:
        y = qmul(y, q)
        k += 1
    return k
orders = [order_of(r) for r in reps]
print("class Re(q):", traces)
print("class element orders:", orders)
out["classes"] = [{"size": s, "Re_q": str(t), "elem_order": o} for s, t, o in zip(sizes, traces, orders)]
# distinct Re(q) per class => character of any SU(2) rep is determined per class by Re(q)
assert len(set(traces)) == len(traces)


# ---------------- exact characters of V_j restricted to Gamma: chi_j = U_{2j}(w), w = Re q
def cheb_U(n, w):
    if n == 0:
        return Q5(1)
    u0, u1 = Q5(1), Q5(2) * w
    for _ in range(n - 1):
        u0, u1 = u1, Q5(2) * w * u1 - u0
    return u1


chiV = {n: [cheb_U(n, w) for w in traces] for n in range(0, 19)}   # level n = 2j

# numeric check that can fail: trace of D^{n/2}(U(q)) equals U_n(Re q) for each class rep and n <= 18
mp = mpmath.mp
mp.dps = 40
maxdev = 0
for r in reps:
    qn = [q5num(c, mp) for c in r]
    U = quat_to_su2(qn, mp.mpf(1), mp.mpc(0, 1))
    for n in range(0, 19):
        D = Dj_generic(n, U, lambda k: mp.factorial(k), mp.sqrt, mp.mpc(0))
        tr = sum(D[i][i] for i in range(n + 1))
        maxdev = max(maxdev, abs(tr - q5num(chiV[n][reps.index(r)], mp)))
print("max |tr D^j - U_2j(Re q)| over classes, n<=18:", mpmath.nstr(maxdev, 5))
out["check_trace_vs_chebyshev_maxdev"] = mpmath.nstr(maxdev, 5)
assert maxdev < mp.mpf(10) ** (-30)

# ---------------- decompose V_3 numerically into Gamma-irreducibles (invariant Hermitian averaging)
Dall = []
for q in G:
    qn = [q5num(c, mp) for c in q]
    U = quat_to_su2(qn, mp.mpf(1), mp.mpc(0, 1))
    Dall.append(mp.matrix(Dj_generic(6, U, lambda k: mp.factorial(k), mp.sqrt, mp.mpc(0))))
H = mp.matrix(7, 7)
for i in range(7):
    for k in range(7):
        H[i, k] = mp.mpc(mp.mpf(1 + i * 7 + 2 * k * k) / (13 + i + k), mp.mpf(i - k) / (5 + i + 2 * k))
H = (H + H.H) / 2
A = mp.matrix(7, 7)
for D in Dall:
    A += D * H * D.H
A /= order
E, Qm = mp.eighe(A)
print("eigenvalues of averaged Hermitian:", [mpmath.nstr(e, 12) for e in E])
# cluster
ev = sorted([(E[i], i) for i in range(7)], key=lambda t: t[0])
clusters = []
for e, i in ev:
    if clusters and abs(e - clusters[-1][0][0]) < mp.mpf(10) ** (-25):
        clusters[-1].append((e, i))
    else:
        clusters.append([(e, i)])
dims = sorted(len(c) for c in clusters)
print("cluster dimensions:", dims)
assert dims == [3, 4]
sectors = {}
for c in clusters:
    d = len(c)
    eta = mp.matrix(7, d)
    for a, (e, i) in enumerate(c):
        for r in range(7):
            eta[r, a] = Qm[r, i]
    sectors[d] = eta


def z5_identify(x):
    """x = p + q*phi with integers |p|,|q| <= 10 (character values are algebraic integers of Q(sqrt5))."""
    phi = (1 + mp.sqrt(5)) / 2
    best = None
    for p in range(-10, 11):
        for q in range(-10, 11):
            err = abs(x - (p + q * phi))
            if best is None or err < best[0]:
                best = (err, p, q)
    err, p, q = best
    assert err < mp.mpf(10) ** (-25), (x, best)
    return Q5(F(p) + F(q, 2), F(q, 2))


chis = {}
for d, eta in sectors.items():
    # invariance residual D eta = eta sigma, sigma = eta^dag D eta
    res = 0
    vals = []
    for r in reps:
        idx = G.index(r)
        D = Dall[idx]
        sig = eta.H * D * eta
        res = max(res, mp.mnorm(D * eta - eta * sig, 1))
        vals.append(sum(sig[a, a] for a in range(d)))
    print(f"sector d={d}: max invariance residual {mpmath.nstr(res, 5)}; imag parts max",
          mpmath.nstr(max(abs(mp.im(v)) for v in vals), 5))
    chis[d] = [z5_identify(mp.re(v)) for v in vals]
    print(f"  character (exact, per class):", chis[d])


def inner(c1, c2):
    s = Q5(0)
    for sz, a, b in zip(sizes, c1, c2):
        s = s + Q5(sz) * a * b      # characters real (values in Q(sqrt5) subset of R)
    return s * Q5(F(1, order))


# exact checks that can fail
for d in (3, 4):
    nrm = inner(chis[d], chis[d])
    print(f"<chi_{d},chi_{d}> =", nrm, " chi(1) =", chis[d][0])
    assert nrm == Q5(1) and chis[d][0] == Q5(d)
sumok = all(chis[3][i] + chis[4][i] == chiV[6][i] for i in range(len(reps)))
print("chi_3 + chi_4 == chi_{V_3} exactly on every class:", sumok)
assert sumok
out["character_3dim"] = [str(x) for x in chis[3]]
out["character_4dim"] = [str(x) for x in chis[4]]
out["character_V3"] = [str(x) for x in chiV[6]]

mult = {}
for d in (3, 4):
    row = []
    for n in range(0, 19):
        m = inner(chis[d], chiV[n])
        assert m.b == 0 and m.a.denominator == 1 and m.a >= 0, m
        row.append(int(m.a))
    mult[d] = row
triv = []
for n in range(0, 19):
    m = inner([Q5(1)] * len(reps), chiV[n])
    assert m.b == 0 and m.a.denominator == 1
    triv.append(int(m.a))
print("level n:          ", list(range(19)))
print("dim Hom(3-dim, V):", mult[3])
print("dim Hom(4-dim, V):", mult[4])
print("dim of invariants:", triv)
out["mult_3dim"] = {str(n): mult[3][n] for n in range(19)}
out["mult_4dim"] = {str(n): mult[4][n] for n in range(19)}
out["mult_trivial"] = {str(n): triv[n] for n in range(19)}

# ---------------- item 0(b)
S, Aq = cg_exact(3, 3, 3, -3, 6, 0)
print("<3 3; 3 -3 | 6 0> = ", S, "* sqrt(", Aq, ")  =", float(S) * math.sqrt(float(Aq)))
val2 = S * S * Aq
out["cg_33_3m3_60_squared"] = str(val2)
out["cg_33_3m3_60_sign"] = "+" if S > 0 else "-"
out["cg_33_3m3_60_float"] = float(S) * math.sqrt(float(Aq))
# cross-check with sympy (independent implementation) on a grid of coefficients (can fail)
from sympy.physics.quantum.cg import CG
from sympy import Rational, nsimplify
bad = 0
cnt = 0
for (j1, j2, J) in [(3, 3, 6), (3, 3, 0), (3, 3, 3), (6, 3, 3), (6, 3, 7), (3, 3, 5), (9, 3, 6), (4, 6, 3)]:
    for m1 in range(-j1, j1 + 1):
        for m2 in range(-j2, j2 + 1):
            M = m1 + m2
            if abs(M) > J:
                continue
            s = CG(j1, m1, j2, m2, J, M).doit()
            S_, A_ = cg_exact(j1, m1, j2, m2, J, M)
            mine2 = S_ * S_ * A_
            sy2 = Rational(s ** 2)
            cnt += 1
            if F(int(sy2.p), int(sy2.q)) != mine2 or ((s > 0) != (S_ > 0) and S_ != 0):
                bad += 1
print(f"Racah CG vs sympy CG: {cnt} coefficients compared, {bad} mismatches")
out["cg_vs_sympy"] = {"compared": cnt, "mismatches": bad}
assert bad == 0
(HERE / "out_item01.json").write_text(json.dumps(out, indent=1))
