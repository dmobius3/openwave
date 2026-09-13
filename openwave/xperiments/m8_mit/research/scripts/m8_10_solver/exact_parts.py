"""Exact (symbolic) route for items 2, 3, 4, 5, 7, 8.

1. The sector projector P = (d/|G|) sum_h chi(h) D^3(h) computed EXACTLY: D^3 = S D' S^-1 with D' the monomial-basis matrix
   (entries in Q(sqrt5)(i), exact rational arithmetic) and S = diag(sqrt((3+m)!(3-m)!)).
2. T_L = sum_b [eta_b (x) Theta eta_b]_L = sum_{k1+k2=M} CG (-1)^{k2} P_{k1,-k2} (depends on P only), exact in sympy.
3. r_L = tr(Op3(T_L) P)/d, F_J = tr(OpJ(T_6)^dag OpJ(T_6) P), OpJ(T)_{Kk} = <6 K-k; 3 k|J K> T_{K-k}.
4. Structure (derived in RETURN.md): for J != 3, Pi_J N(Phi_{sigma,u}) has coefficients M^{6,J}(u)_M (R^{6,J})_{Ka},
   so ||Pi_J N||^2 = |M^{6,J}(u)|^2 F_J/(2J+1); kappa = 1 + r_6 (7/d) <v, M_6(v)> for unit v, u = sqrt(7/d) v.
Every symbolic quantity is checked to be an exact rational (or exact zero) and compared to the high-precision engine."""
import json
from fractions import Fraction as Fr
import math
import sympy as sp
from common import HERE, Q5, generate_gamma, cg_exact

items01 = json.loads((HERE / "out_item01.json").read_text())


def parse_q5(s):
    s = s.strip("()")
    a, b = s.split(" + ")
    return Q5(Fr(a), Fr(b.replace("*sqrt5", "")))


class C5:
    """re + i im, re, im in Q(sqrt5)."""
    __slots__ = ("re", "im")

    def __init__(self, re, im=None):
        self.re = re if isinstance(re, Q5) else Q5(re)
        self.im = im if isinstance(im, Q5) else Q5(0 if im is None else im)

    def __add__(s, o):
        return C5(s.re + o.re, s.im + o.im)

    def __sub__(s, o):
        return C5(s.re - o.re, s.im - o.im)

    def __mul__(s, o):
        if not isinstance(o, C5):
            o = C5(o)
        return C5(s.re * o.re - s.im * o.im, s.re * o.im + s.im * o.re)

    __rmul__ = __mul__

    def conj(s):
        return C5(s.re, -s.im)

    def __eq__(s, o):
        return s.re == o.re and s.im == o.im

    def is0(s):
        return s.re == Q5(0) and s.im == Q5(0)

    def __pow__(s, n):
        r = C5(1)
        for _ in range(n):
            r = r * s
        return r


def mono_D(n, q):
    w, x, y, z = q
    U00, U01 = C5(w, x), C5(y, z)
    U10, U11 = C5(-y, z), C5(w, -x)
    D = [[C5(0) for _ in range(n + 1)] for _ in range(n + 1)]
    for kk in range(n + 1):
        p, qq = kk, n - kk
        for s in range(p + 1):
            for t in range(qq + 1):
                c = math.comb(p, s) * math.comb(qq, t)
                D[s + t][kk] = D[s + t][kk] + (U00 ** s) * (U10 ** (p - s)) * (U01 ** t) * (U11 ** (qq - t)) * C5(c)
    return D


def matmul(A, B):
    n, m, p = len(A), len(B), len(B[0])
    return [[sum((A[i][k] * B[k][j] for k in range(m)), C5(0)) for j in range(p)] for i in range(n)]


G = generate_gamma()
cls_re = [parse_q5(c["Re_q"]) for c in items01["classes"]]
fac = [math.factorial(i) * math.factorial(6 - i) for i in range(7)]      # f_i, i = m + 3
gens = [G[0]]  # placeholder
from common import Q1, Q2
Dq1, Dq2 = mono_D(6, Q1), mono_D(6, Q2)
out = {}
Pi = {}
for d, key in ((3, "character_3dim"), (4, "character_4dim")):
    chi = [parse_q5(x) for x in items01[key]]
    acc = [[C5(0) for _ in range(7)] for _ in range(7)]
    for q in G:
        Dm = mono_D(6, q)
        c = chi[cls_re.index(q[0])]
        for i in range(7):
            for k in range(7):
                acc[i][k] = acc[i][k] + Dm[i][k] * C5(c)
    P = [[acc[i][k] * C5(Q5(Fr(d, len(G)))) for k in range(7)] for i in range(7)]
    # exact checks that can fail
    P2 = matmul(P, P)
    idem = all((P2[i][k] - P[i][k]).is0() for i in range(7) for k in range(7))
    comm = all((matmul(P, Dg)[i][k] - matmul(Dg, P)[i][k]).is0() for Dg in (Dq1, Dq2) for i in range(7) for k in range(7))
    tr = sum((P[i][i] for i in range(7)), C5(0))
    herm = all((P[k][i].conj() * C5(fac[k]) - P[i][k] * C5(fac[i])).is0() for i in range(7) for k in range(7))
    print(f"sector {d}: Pi^2 = Pi {idem}, [Pi, D'(q1)] = [Pi, D'(q2)] = 0 {comm}, tr = {tr.re}+i{tr.im}, P hermitian {herm}")
    assert idem and comm and tr == C5(d) and herm
    Pi[d] = P
# P3 + P4 = identity (exact)
assert all((Pi[3][i][k] + Pi[4][i][k] - C5(1 if i == k else 0)).is0() for i in range(7) for k in range(7))
print("P3 + P4 = I exactly: True")

s5 = sp.sqrt(5)


def to_sp(c):
    return (sp.Rational(c.re.a) + sp.Rational(c.re.b) * s5) + sp.I * (sp.Rational(c.im.a) + sp.Rational(c.im.b) * s5)


def cgs(j1, m1, j2, m2, J, M):
    S, A = cg_exact(j1, m1, j2, m2, J, M)
    if S == 0:
        return sp.Integer(0)
    return sp.Rational(S) * sp.sqrt(sp.Rational(A))


def simp(x):
    return sp.nsimplify(sp.expand(sp.radsimp(sp.expand(x))))


def exact_rational(x, label):
    y = sp.expand(sp.radsimp(sp.expand(x)))
    if not y.is_Rational:
        y2 = sp.simplify(y)
        if not y2.is_Rational:
            raise ValueError(f"{label}: not reduced to a rational: {y}")
        y = y2
    return sp.Rational(y)


Psp = {}
for d in (3, 4):
    Psp[d] = sp.Matrix(7, 7, lambda i, k: to_sp(Pi[d][i][k]) * sp.sqrt(sp.Rational(fac[i], fac[k])))


def Tvec(P, L):
    T = []
    for M in range(-L, L + 1):
        t = sp.Integer(0)
        for k1 in range(-3, 4):
            k2 = M - k1
            if abs(k2) > 3:
                continue
            t += cgs(3, k1, 3, k2, L, M) * (-1) ** (k2 % 2) * P[k1 + 3, -k2 + 3]
        T.append(sp.expand(t))
    return T


def Op(T, L, J):
    return sp.Matrix(2 * J + 1, 7, lambda iK, ik: (cgs(L, (iK - J) - (ik - 3), 3, ik - 3, J, iK - J) * T[(iK - J) - (ik - 3) + L]
                                                  if abs((iK - J) - (ik - 3)) <= L else sp.Integer(0)))


sector_consts = {}
for d in (3, 4):
    P = Psp[d]
    rec = {}
    TL = {}
    for L in range(0, 7):
        T = Tvec(P, L)
        n2 = exact_rational(sum(sp.expand(t * sp.conjugate(t)) for t in T), f"|T_{L}|^2")
        TL[L] = T
        rec[f"T{L}_norm2"] = str(n2)
    print(f"sector {d}: |T_L|^2 for L=0..6:", [rec[f'T{L}_norm2'] for L in range(7)])
    r = {}
    for L in (0, 6):
        O = Op(TL[L], L, 3)
        rv = sp.expand(sp.radsimp(sp.expand((O * P).trace() / d)))
        r2 = exact_rational(rv * sp.conjugate(rv), f"r_{L}^2")
        # real and sign
        rnum = complex(sp.N(rv, 50))
        assert abs(rnum.imag) < 1e-40
        r[L] = (r2, 1 if rnum.real > 0 else -1)
        rec[f"r{L}_squared"] = str(r2)
        rec[f"r{L}_sign"] = "+" if rnum.real > 0 else "-"
        print(f"  r_{L} = {'+' if rnum.real > 0 else '-'}sqrt({r2}) = {rnum.real:.15f}")
    F = {}
    for J in range(3, 10):
        O = Op(TL[6], 6, J)
        val = exact_rational((O.H * O * P).trace(), f"F_{J}")
        F[J] = val
        rec[f"F_{J}"] = str(val)
    print("  F_J (J=3..9):", {J: str(F[J]) for J in F})
    sector_consts[d] = {"r": r, "F": F}
    out[f"sector_{d}"] = rec


# ---------------------------------------------------------------- ray parts, exact
def theta(u):
    return [(-1) ** ((i - 3) % 2) * sp.conjugate(u[6 - i]) for i in range(7)]


def couple(x, j1, y, j2, J):
    res = []
    for M in range(-J, J + 1):
        t = sp.Integer(0)
        for m1 in range(-j1, j1 + 1):
            m2 = M - m1
            if abs(m2) <= j2:
                c = cgs(j1, m1, j2, m2, J, M)
                if c != 0:
                    t += c * x[m1 + j1] * y[m2 + j2]
        res.append(sp.expand(t))
    return res


def M6J(v, J):
    return couple(couple(v, 3, theta(v), 3, 6), 6, v, 3, J)


e = lambda m: [sp.Integer(1) if i - 3 == m else sp.Integer(0) for i in range(7)]
RAYS = {
    "R1": e(3), "R2": e(0),
    "R3": [(x + y) / sp.sqrt(2) for x, y in zip(e(2), e(-2))],
    "R4": [(x + y) / sp.sqrt(2) for x, y in zip(e(3), e(-3))],
    "R5": [sp.sqrt(13) / 5 * x + 2 * sp.sqrt(3) / 5 * y for x, y in zip(e(2), e(-3))],
}
for name, v in RAYS.items():
    assert exact_rational(sum(sp.expand(x * sp.conjugate(x)) for x in v), "norm") == 1

engine = json.loads((HERE / "out_engine_110.json").read_text())
mu = {J: 4 * J * (J + 1) - 48 for J in range(10)}
sector_levels = {3: [J for J in range(10) if items01["mult_3dim"][str(2 * J)] > 0],
                 4: [J for J in range(10) if items01["mult_4dim"][str(2 * J)] > 0]}
out["rays"] = {}
raypart = {}
for name, v in RAYS.items():
    M6 = M6J(v, 3)
    proj = sp.radsimp(sp.expand(sum(sp.expand(sp.conjugate(a) * b) for a, b in zip(v, M6))))
    # must be exactly (rational) * sqrt(91): check that proj/sqrt(91) is rational
    exact_rational(proj / sp.sqrt(91), "<v,M6 v>/sqrt(91)")
    perp = [sp.expand(b - proj * a) for a, b in zip(v, M6)]
    perp2 = exact_rational(sum(sp.expand(x * sp.conjugate(x)) for x in perp), "|M6 perp|^2")
    mJ = {}
    zero_comps = {}
    for J in range(3, 10):
        comps = M6J(v, J)
        mJ[J] = exact_rational(sum(sp.expand(x * sp.conjugate(x)) for x in comps), f"|M^(6,{J})|^2")
        zero_comps[J] = all(sp.expand(sp.radsimp(x)) == 0 for x in comps)
    raypart[name] = (proj, perp2, mJ)
    print(f"{name}: <v,M6(v)> = {proj}, |M6(v) - <v,M6 v> v|^2 = {perp2}, |M^(6,J)(v)|^2 J=3..9: {[str(mJ[J]) for J in range(3, 10)]}")
    out["rays"][name] = {"v_M6v": str(proj), "M6_perp_norm2": str(perp2),
                         "M6J_norm2": {str(2 * J): str(mJ[J]) for J in mJ},
                         "M6J_identically_zero": {str(2 * J): zero_comps[J] for J in zero_comps}}

# assemble exact values and compare to the 110-digit engine
import mpmath
mpmath.mp.dps = 110
assembled = {}
maxdev = 0
for d in (3, 4):
    r6 = sector_consts[d]["r"][6]
    F = sector_consts[d]["F"]
    for name in RAYS:
        proj, perp2, mJ = raypart[name]
        s = sp.Rational(7, d)
        # kappa = 1 + r6 * (7/d) * <v, M6 v>; r6 = sign*sqrt(r6^2): exact value kappa - 1 = sign*sqrt(r6^2)*(7/d)*proj
        km1 = r6[1] * sp.sqrt(r6[0]) * s * proj
        km1 = exact_rational(km1, "kappa - 1")
        kappa = 1 + km1
        # |Pi_6 N - kappa Phi|^2 = (d/7) |w(u) - kappa u|^2 = (d/7) r6^2 (7/d)^3 |M6(v) - <v,M6 v> v|^2
        par2 = sp.Rational(d, 7) * r6[0] * s ** 3 * perp2
        Nn = {J: s ** 3 * mJ[J] * F[J] / (2 * J + 1) for J in range(4, 10)}
        xin = {2 * J: Nn[J] / mu[J] ** 2 for J in Nn}
        lam4 = -3 * sum(Nn[J] / mu[J] for J in Nn)
        e_ = engine["rays"][f"{d}:{name}"]
        devs = [abs(mpmath.mpf(e_["kappa_re"]) - mpmath.mpf(sp.N(kappa, 120))),
                abs(mpmath.mpf(e_["lambda4_over_g2"]) - mpmath.mpf(sp.N(lam4, 120)))]
        for J in range(4, 10):
            devs.append(abs(mpmath.mpf(e_["xi_level_norm2_over_g2"][str(2 * J)]) - mpmath.mpf(sp.N(xin[2 * J], 120))))
        dv = max(devs)
        maxdev = max(maxdev, dv)
        assembled[f"{d}:{name}"] = {"kappa": str(kappa), "kappa_minus_1": str(km1), "parallel_resid2": str(par2),
                                    "xi_norm2_over_g2": {str(n): str(xin[n]) for n in xin},
                                    "sector_levels_n": [2 * J for J in sector_levels[d]],
                                    "lambda4_over_g2": str(lam4), "max_dev_vs_engine110": mpmath.nstr(dv, 5)}
        print(f"[d={d} {name}] kappa = {kappa}, lambda4/g^2 = {lam4}, ||Pi_n xi||^2/g^2 = "
              f"{ {n: str(xin[n]) for n in xin} }, dev vs engine {mpmath.nstr(dv, 3)}")
print("max deviation exact vs 110-digit engine:", mpmath.nstr(maxdev, 5))
assert maxdev < mpmath.mpf(10) ** (-95)
out["assembled"] = assembled

# item 8 ratios
ratios = {}
for name in RAYS:
    l3 = sp.Rational(assembled[f"3:{name}"]["lambda4_over_g2"])
    l4 = sp.Rational(assembled[f"4:{name}"]["lambda4_over_g2"])
    k3 = sp.nsimplify(assembled[f"3:{name}"]["kappa_minus_1"])
    k4 = sp.nsimplify(assembled[f"4:{name}"]["kappa_minus_1"])
    ratios[name] = {"lambda4_ratio_3_over_4": str(l3 / l4), "kappa_minus_1_ratio_3_over_4": str(sp.nsimplify(k3 / k4))}
    print(f"{name}: lambda4(3)/lambda4(4) = {l3 / l4} = {float(l3 / l4):.10f};  (kappa-1)(3)/(kappa-1)(4) = {sp.nsimplify(k3 / k4)}")
out["item8"] = ratios
(HERE / "out_exact.json").write_text(json.dumps(out, indent=1))
print("wrote out_exact.json")
