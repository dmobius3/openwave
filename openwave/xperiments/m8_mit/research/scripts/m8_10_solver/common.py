"""Shared exact and high-precision helpers. Paths relative to this file only."""
import pathlib
from fractions import Fraction as F
from functools import lru_cache
import math

HERE = pathlib.Path(__file__).parent


# ---------------------------------------------------------------- Q(sqrt5)
class Q5:
    """a + b*sqrt(5), a, b rational (exact)."""
    __slots__ = ("a", "b")

    def __init__(self, a, b=0):
        self.a = F(a)
        self.b = F(b)

    def __add__(s, o):
        o = o if isinstance(o, Q5) else Q5(o)
        return Q5(s.a + o.a, s.b + o.b)

    __radd__ = __add__

    def __neg__(s):
        return Q5(-s.a, -s.b)

    def __sub__(s, o):
        return s + (-(o if isinstance(o, Q5) else Q5(o)))

    def __rsub__(s, o):
        return Q5(o) - s

    def __mul__(s, o):
        o = o if isinstance(o, Q5) else Q5(o)
        return Q5(s.a * o.a + 5 * s.b * o.b, s.a * o.b + s.b * o.a)

    __rmul__ = __mul__

    def __eq__(s, o):
        o = o if isinstance(o, Q5) else Q5(o)
        return s.a == o.a and s.b == o.b

    def __hash__(s):
        return hash((s.a, s.b))

    def __repr__(s):
        return f"({s.a} + {s.b}*sqrt5)"


def q5num(x, mp):
    return mp.mpf(x.a.numerator) / x.a.denominator + mp.mpf(x.b.numerator) / x.b.denominator * mp.sqrt(5)


def q5float(x):
    return float(x.a) + float(x.b) * math.sqrt(5)


HALF = F(1, 2)
PHI = Q5(HALF, HALF)          # (1+sqrt5)/2
PHIINV = Q5(-HALF, HALF)      # phi - 1 = 1/phi


# ---------------------------------------------------------------- quaternions (w, x, y, z) over Q5
def qmul(p, q):
    a1, b1, c1, d1 = p
    a2, b2, c2, d2 = q
    return (a1 * a2 - b1 * b2 - c1 * c2 - d1 * d2,
            a1 * b2 + b1 * a2 + c1 * d2 - d1 * c2,
            a1 * c2 - b1 * d2 + c1 * a2 + d1 * b2,
            a1 * d2 + b1 * c2 - c1 * b2 + d1 * a2)


def qconj(q):
    return (q[0], -q[1], -q[2], -q[3])


def qnorm2(q):
    return q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]


ONE = (Q5(1), Q5(0), Q5(0), Q5(0))
Q1 = (Q5(HALF), Q5(HALF), Q5(HALF), Q5(HALF))
Q2 = (PHI * Q5(HALF), PHIINV * Q5(HALF), Q5(HALF), Q5(0))


def close_group(gens):
    """Closure under multiplication of a finite set of unit quaternions (exact)."""
    elems = {ONE}
    frontier = [ONE]
    while frontier:
        new = []
        for x in frontier:
            for g in gens:
                y = qmul(x, g)
                if y not in elems:
                    elems.add(y)
                    new.append(y)
        frontier = new
    return elems


def generate_gamma():
    return sorted(close_group([Q1, Q2]), key=lambda q: tuple((c.a, c.b) for c in q))


# ---------------------------------------------------------------- Clebsch-Gordan (Racah formula, integer or half-integer via doubled args)
@lru_cache(maxsize=None)
def cg_exact(j1, m1, j2, m2, J, M):
    """Condon-Shortley CG for INTEGER spins. Returns (S, A) with CG = S*sqrt(A), S, A Fractions, A>0 or CG=0 -> (0,1)."""
    if m1 + m2 != M or abs(m1) > j1 or abs(m2) > j2 or abs(M) > J:
        return (F(0), F(1))
    if J < abs(j1 - j2) or J > j1 + j2:
        return (F(0), F(1))
    f = math.factorial
    A = F((2 * J + 1) * f(J + j1 - j2) * f(J - j1 + j2) * f(j1 + j2 - J), f(j1 + j2 + J + 1))
    A *= f(J + M) * f(J - M) * f(j1 - m1) * f(j1 + m1) * f(j2 - m2) * f(j2 + m2)
    S = F(0)
    for k in range(0, j1 + j2 + J + 2):
        args = [k, j1 + j2 - J - k, j1 - m1 - k, j2 + m2 - k, J - j2 + m1 + k, J - j1 - m2 + k]
        if min(args) < 0:
            continue
        den = 1
        for x in args:
            den *= f(x)
        S += F((-1) ** k, den)
    if S == 0:
        return (F(0), F(1))
    return (S, A)


def cg_mp(j1, m1, j2, m2, J, M, mp):
    S, A = cg_exact(j1, m1, j2, m2, J, M)
    if S == 0:
        return mp.mpf(0)
    return mp.mpf(S.numerator) / S.denominator * mp.sqrt(mp.mpf(A.numerator) / A.denominator)


def cg_float(j1, m1, j2, m2, J, M):
    S, A = cg_exact(j1, m1, j2, m2, J, M)
    return float(S) * math.sqrt(float(A)) if S != 0 else 0.0


# ---------------------------------------------------------------- SU(2) matrices and D^j (symmetric power realization, CS basis)
def quat_to_su2(q, one, I):
    """q = w + x i + y j + z k  ->  [[w + i x, y + i z], [-y + i z, w - i x]] (i->diag(i,-i), j->[[0,1],[-1,0]], k->[[0,i],[i,0]])."""
    w, x, y, z = q
    return [[w + I * x, y + I * z], [-y + I * z, w - I * x]]


def Dj_generic(j2, U, fact, sqrt, zero):
    """D^j(U) with j = j2/2 in the basis e_m = x^{j+m} y^{j-m}/sqrt((j+m)!(j-m)!), indices i = m + j (0..2j).
    g x = U00 x + U10 y, g y = U01 x + U11 y; D_{mk} = coefficient of e_m in g e_k.
    Only integer j used by callers of the mp/float versions (j2 even), but formula works for any j2."""
    n = j2
    U00, U01 = U[0][0], U[0][1]
    U10, U11 = U[1][0], U[1][1]
    D = [[zero for _ in range(n + 1)] for _ in range(n + 1)]
    for kk in range(n + 1):          # kk = j + k : power of x in e_k
        p, qq = kk, n - kk           # e_k = x^p y^qq
        for s in range(p + 1):
            for t in range(qq + 1):
                mm = s + t           # power of x in result = j + m
                c = math.comb(p, s) * math.comb(qq, t)
                term = c * (U00 ** s) * (U10 ** (p - s)) * (U01 ** t) * (U11 ** (qq - t))
                D[mm][kk] = D[mm][kk] + term * sqrt(fact(mm) * fact(n - mm)) / sqrt(fact(kk) * fact(n - kk))
    return D


def spin_matrices_generic(J, sqrt, zero, one):
    """Jz, J+, J- for integer J, index i = m + J, CS convention (J+ has positive real entries)."""
    n = 2 * J + 1
    Jz = [[zero] * n for _ in range(n)]
    Jp = [[zero] * n for _ in range(n)]
    Jm = [[zero] * n for _ in range(n)]
    for i in range(n):
        m = i - J
        Jz[i][i] = one * m
        if m + 1 <= J:
            Jp[i + 1][i] = sqrt((J - m) * (J + m + 1))
        if m - 1 >= -J:
            Jm[i - 1][i] = sqrt((J + m) * (J - m + 1))
    return Jz, Jp, Jm


# ---------------------------------------------------------------- identification
def identify_rational(x, mp, bound, tol):
    """Best rational with denominator <= bound; accepted if |x - p/q| < tol. Returns Fraction or None."""
    s = mp.nstr(x, mp.dps, strip_zeros=False)
    fr = F(s).limit_denominator(bound)
    err = abs(x - mp.mpf(fr.numerator) / fr.denominator)
    return (fr if err < tol else None), err


def identify_q5(x, mp, maxcoeff, tol):
    """Find integers (p, q, r), r>0, with x = (p + q sqrt5)/r via PSLQ on [x, 1, sqrt5]."""
    rel = mp.pslq([x, mp.mpf(1), mp.sqrt(5)], maxcoeff=maxcoeff, maxsteps=10 ** 6)
    if rel is None or rel[0] == 0:
        return None, None
    r, p, q = rel[0], -rel[1], -rel[2]
    if r < 0:
        r, p, q = -r, -p, -q
    val = (mp.mpf(p) + q * mp.sqrt(5)) / r
    err = abs(x - val)
    return ((F(p, r), F(q, r)) if err < tol else None), err
