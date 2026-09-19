"""Items 0 and 1 (exact part): group, derived subgroup, isotypic projectors of V_3, Hom dimensions."""
import pickle, os, json
from fractions import Fraction
from mq import MQ, ZERO, ONE, mq
from su2 import qmul, qinv, qkey, q_to_su2, D_exact, cg, cg_table_hw
from checks import check, save_log

os.makedirs("out", exist_ok=True)
R = {}
s5 = MQ.sqrt(5)
half = Fraction(1, 2)
phi = (ONE + s5) * half
phiinv = (s5 - ONE) * half
q1 = (mq(half), mq(half), mq(half), mq(half))
q2 = (phi * half, phiinv * half, mq(half), ZERO)
E = (ONE, ZERO, ZERO, ZERO)

def qnorm2(q):
    return sum((x * x for x in q), ZERO)

# ---- 0(a): norms
n1, n2 = qnorm2(q1), qnorm2(q2)
check("||q1||^2 == 1 exactly", lambda v: v == ONE, n1, n1 + Fraction(1, 10**6), "add 1e-6")
check("||q2||^2 == 1 exactly", lambda v: v == ONE, n2, qnorm2((phi * half, phi * half, mq(half), ZERO)), "q2 with phi in place of 1/phi")
R["norm_q1_sq"] = n1.pretty(); R["norm_q2_sq"] = n2.pretty()

# ---- closure
def closure(gens):
    elems = {qkey(E): E}
    frontier = [E]
    while frontier:
        nf = []
        for a in frontier:
            for g in gens:
                b = qmul(a, g)
                k = qkey(b)
                if k not in elems:
                    elems[k] = b
                    nf.append(b)
        frontier = nf
    return list(elems.values())

G = closure([q1, q2])
order = len(G)
print("order of Gamma:", order)
R["order"] = order
# sanity: closure really closed (all products inside), and inverses inside
keys = {qkey(g) for g in G}
closed = all(qkey(qmul(a, b)) in keys for a in G[:30] for b in G)
check("group closed under products (30 x all spot check) and every element unit", lambda c: c[0] and c[1],
      (closed, all(qnorm2(g) == ONE for g in G)), (closed, False), "pretend one norm != 1")

# ---- derived subgroup
comms = []
ck = set()
for a in G:
    for b in G:
        c = qmul(qmul(a, b), qmul(qinv(a), qinv(b)))
        k = qkey(c)
        if k not in ck:
            ck.add(k); comms.append(c)
Dg = closure(comms)
print("order of commutator subgroup:", len(Dg), " #distinct commutators:", len(comms))
R["derived_order"] = len(Dg)
R["perfect"] = len(Dg) == order
check("derived subgroup order == |Gamma| (Gamma perfect)", lambda x: x[0] == x[1], (len(Dg), order), (len(Dg) - 1, order), "derived order - 1")

# ---- conjugacy classes
cls = []
seen = set()
for g in G:
    if qkey(g) in seen:
        continue
    c = {}
    for h in G:
        x = qmul(qmul(h, g), qinv(h))
        c[qkey(x)] = x
    seen.update(c.keys())
    cls.append(list(c.values()))
cls_info = [(len(c), c[0][0].pretty()) for c in cls]
print("classes (size, Re):", cls_info)
R["classes"] = cls_info
R["class_sum_check"] = sum(len(c) for c in cls)

# ---- 0(b)
v = cg(3, 3, 3, -3, 6, 0)
v2 = cg_table_hw(3, 3, 6)[(3, -3, 0)]
print("<3 3; 3 -3 | 6 0> =", v.pretty(), " (HW construction:", v2.pretty(), ")")
R["cg_33_3m3_60"] = v.pretty()
check("CG <33;3-3|60>: Racah == highest-weight construction", lambda p: p[0] == p[1], (v, v2), (v, v2 * 2), "double HW value")
check("CG <33;3-3|60> > 0 (Condon-Shortley)", lambda x: x.sign() > 0, v, -v, "negate")

# ---- D^3 of all elements, exact
D3 = [D_exact(3, q_to_su2(g)) for g in G]
def matmul(A, B):
    n, m, p = len(A), len(B), len(B[0])
    return [[sum((A[i][k] * B[k][j] for k in range(m) if not A[i][k].is_zero() and not B[k][j].is_zero()), ZERO) for j in range(p)] for i in range(n)]
def madd(A, B):
    return [[a + b for a, b in zip(ra, rb)] for ra, rb in zip(A, B)]
def mscale(A, s):
    return [[a * s for a in r] for r in A]
def meq(A, B):
    return all((a - b).is_zero() for ra, rb in zip(A, B) for a, b in zip(ra, rb))
def ident(n):
    return [[ONE if i == j else ZERO for j in range(n)] for i in range(n)]
def dag(A):
    return [[A[j][i].conj() for j in range(len(A))] for i in range(len(A[0]))]

# homomorphism spot check: D(q1 q2) = D(q1) D(q2)
iq1 = [qkey(g) for g in G].index(qkey(q1)); iq2 = [qkey(g) for g in G].index(qkey(q2))
Dq1q2 = D_exact(3, q_to_su2(qmul(q1, q2)))
check("D^3(q1 q2) == D^3(q1) D^3(q2) exactly", lambda M: meq(M, matmul(D3[iq1], D3[iq2])), Dq1q2, matmul(D3[iq2], D3[iq1]), "use D(q2)D(q1)")
check("D^3(q2) unitary exactly", lambda M: meq(matmul(dag(M), M), ident(7)), D3[iq2], madd(D3[iq2], mscale(ident(7), Fraction(1, 1000))), "add 1e-3 I")

# class of q1
c1 = [c for c in cls if any(qkey(x) == qkey(q1) for x in c)][0]
Z = [[ZERO] * 7 for _ in range(7)]
keyidx = {qkey(g): i for i, g in enumerate(G)}
for x in c1:
    Z = madd(Z, D3[keyidx[qkey(x)]])
Z2 = matmul(Z, Z)
# find c with Z^2 = c Z
cval = None
for i in range(7):
    for j in range(7):
        if not Z[i][j].is_zero():
            cval = Z2[i][j] / Z[i][j]; break
    if cval is not None:
        break
print("class size of q1:", len(c1), " Z^2 = c Z with c =", cval.pretty())
check("Z^2 == c Z exactly (class sum has eigenvalues 0 and c)", lambda cc: meq(Z2, mscale(Z, cc)), cval, cval + 1, "c+1")
Pz = mscale(Z, cval.inv())
Pc = madd(ident(7), mscale(Pz, -1))
tr = lambda A: sum((A[i][i] for i in range(len(A))), ZERO)
print("rank(Z/c) =", tr(Pz).pretty(), " rank(I - Z/c) =", tr(Pc).pretty())
P = {}
for M in (Pz, Pc):
    P[int(tr(M).t[1][0])] = M
R["class_q1_size"] = len(c1); R["Zc"] = cval.pretty()
for d, M in P.items():
    check(f"P_{d}^2 == P_{d}", lambda A: meq(matmul(A, A), A), M, mscale(M, 2), "2P")
    check(f"P_{d} commutes with D3(q1), D3(q2)", lambda A: meq(matmul(A, D3[iq1]), matmul(D3[iq1], A)) and meq(matmul(A, D3[iq2]), matmul(D3[iq2], A)),
          M, madd(M, [[ONE if (i, j) == (0, 1) else ZERO for j in range(7)] for i in range(7)]), "add E_01")
    check(f"P_{d} Hermitian", lambda A: meq(dag(A), A), M, mscale(M, MQ.I), "multiply by i")

# characters and irreducibility
chi = {d: [tr(matmul(P[d], D3[i])) for i in range(order)] for d in P}
for d in P:
    nrm = sum((x * x.conj() for x in chi[d]), ZERO) * Fraction(1, order)
    print(f"sector {d}: <chi,chi> =", nrm.pretty(), " chi real:", all(x.is_real() for x in chi[d]))
    check(f"sector {d} irreducible: <chi,chi> == 1", lambda v: v == ONE, nrm, nrm * 2, "double")
    check(f"sector {d} character real-valued", lambda cs: all(x.is_real() for x in cs), chi[d], [chi[d][0] * MQ.I] + chi[d][1:], "multiply one value by i")
    R[f"chi_{d}_by_class"] = [(len(c), c[0][0].pretty(), chi[d][keyidx[qkey(c[0])]].pretty()) for c in cls]

# chi_j(h) = U_{2j}(w)
def cheb_U(n, w):
    if n == 0:
        return ONE
    a, b = ONE, w * 2
    for _ in range(n - 1):
        a, b = b, w * 2 * b - a
    return b

dims = {}
for d in P:
    dims[d] = {}
    for twoj in range(0, 19):
        s = ZERO
        for i, g in enumerate(G):
            s = s + chi[d][i].conj() * cheb_U(twoj, g[0])
        s = s * Fraction(1, order)
        dims[d][twoj] = s
    print(f"sector {d}: dim Hom(sigma, V_(n/2)) for n=0..18:", [dims[d][n].pretty() for n in range(19)])
R["dims"] = {d: {n: dims[d][n].pretty() for n in dims[d]} for d in dims}
# check: sum over sectors weighted... check trace formula of chi_j vs actual D^j trace at q2 for j=3
chk = tr(D3[iq2]) - cheb_U(6, q2[0])
check("tr D^3(q2) == U_6(Re q2)", lambda x: x.is_zero(), chk, chk + 1, "+1")
check("dim Hom(sigma, V_3) == 1 for both sectors", lambda dd: all(dd[k][6] == ONE for k in dd), dims, {3: {6: mq(2)}, 4: dims[4]}, "set 2")

with open("out/group.pkl", "wb") as f:
    pickle.dump({"G": G, "q1": q1, "q2": q2, "P": P, "dims": {d: {n: dims[d][n] for n in dims[d]} for d in dims},
                 "classes": cls, "chi": chi}, f)
with open("out/item01.json", "w") as f:
    json.dump(R, f, indent=1)
save_log("s01")
