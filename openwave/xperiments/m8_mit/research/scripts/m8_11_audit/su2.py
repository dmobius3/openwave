"""SU(2) machinery built from scratch: CG coefficients (two constructions), D^j matrices, the finite group."""
from fractions import Fraction
from math import factorial, comb
from functools import lru_cache
from mq import MQ, ZERO, ONE, mq

# ---------------------------------------------------------------- CG, construction 1: Racah formula
@lru_cache(maxsize=None)
def cg(j1, m1, j2, m2, J, M):
    """<j1 m1; j2 m2 | J M>, Condon-Shortley, integer spins, exact MQ."""
    if m1 + m2 != M or abs(m1) > j1 or abs(m2) > j2 or abs(M) > J:
        return ZERO
    if J < abs(j1 - j2) or J > j1 + j2:
        return ZERO
    f = factorial
    pref = Fraction((2 * J + 1) * f(J + j1 - j2) * f(J - j1 + j2) * f(j1 + j2 - J), f(j1 + j2 + J + 1))
    pref *= f(J + M) * f(J - M) * f(j1 - m1) * f(j1 + m1) * f(j2 - m2) * f(j2 + m2)
    s = Fraction(0)
    for k in range(0, j1 + j2 - J + 1):
        dens = [k, j1 + j2 - J - k, j1 - m1 - k, j2 + m2 - k, J - j2 + m1 + k, J - j1 - m2 + k]
        if min(dens) < 0:
            continue
        d = 1
        for x in dens:
            d *= f(x)
        s += Fraction((-1) ** k, d)
    if s == 0:
        return ZERO
    return MQ.sqrt(pref) * s


# ---------------------------------------------------------------- CG, construction 2: highest weight + lowering
def w_mono(j, m):
    """||x^{j+m} y^{j-m}||^2 in the SU(2)-invariant inner product normalized so that v_m are orthonormal."""
    return Fraction(factorial(j + m) * factorial(j - m), factorial(2 * j))


def _nullspace_rational(A):
    """A: list of rows (Fractions). Return a basis of the nullspace (list of vectors)."""
    rows = [list(r) for r in A]
    ncol = len(rows[0]) if rows else 0
    piv = []
    r = 0
    for c in range(ncol):
        p = None
        for i in range(r, len(rows)):
            if rows[i][c] != 0:
                p = i
                break
        if p is None:
            continue
        rows[r], rows[p] = rows[p], rows[r]
        pv = rows[r][c]
        rows[r] = [x / pv for x in rows[r]]
        for i in range(len(rows)):
            if i != r and rows[i][c] != 0:
                fct = rows[i][c]
                rows[i] = [a - fct * b for a, b in zip(rows[i], rows[r])]
        piv.append(c)
        r += 1
    free = [c for c in range(ncol) if c not in piv]
    basis = []
    for fc in free:
        v = [Fraction(0)] * ncol
        v[fc] = Fraction(1)
        for i, pc in enumerate(piv):
            v[pc] = -rows[i][fc]
        basis.append(v)
    return basis


def cg_table_hw(j1, j2, J):
    """All <j1 m1; j2 m2|J M> by the highest-weight construction. Returns dict (m1,m2,M)->MQ."""
    # monomial basis e_m = x^{j+m} y^{j-m}; J+ = x d/dy : e_m -> (j-m) e_{m+1}; J- = y d/dx : e_m -> (j+m) e_{m-1}
    pairs = [(m1, J - m1) for m1 in range(-j1, j1 + 1) if abs(J - m1) <= j2]
    targets = [(m1, J + 1 - m1) for m1 in range(-j1, j1 + 1) if abs(J + 1 - m1) <= j2]
    A = []
    for (t1, t2) in targets:
        row = []
        for (m1, m2) in pairs:
            c = Fraction(0)
            if (m1 + 1, m2) == (t1, t2):
                c += j1 - m1
            if (m1, m2 + 1) == (t1, t2):
                c += j2 - m2
            row.append(c)
        A.append(row)
    if A:
        ns = _nullspace_rational(A)
    else:
        ns = [[Fraction(1)] * len(pairs)]
    assert len(ns) == 1, (j1, j2, J, len(ns))
    hw = {p: c for p, c in zip(pairs, ns[0])}
    # lower rationally in monomial basis: vec at weight M -> weight M-1
    vecs = {J: hw}
    cur = hw
    for M in range(J, -J, -1):
        nxt = {}
        for (m1, m2), c in cur.items():
            if m1 - 1 >= -j1:
                k = (m1 - 1, m2)
                nxt[k] = nxt.get(k, Fraction(0)) + c * (j1 + m1)
            if m2 - 1 >= -j2:
                k = (m1, m2 - 1)
                nxt[k] = nxt.get(k, Fraction(0)) + c * (j2 + m2)
        vecs[M - 1] = nxt
        cur = nxt
    # convert to orthonormal coordinates: coefficient on v_m1 v_m2 is c * sqrt(w1 w2); normalize each M
    out = {}
    # phase: coefficient of (j1, J-j1) at M=J positive
    ref = hw.get((j1, J - j1), Fraction(0))
    assert ref != 0
    sgn = 1 if ref > 0 else -1
    for M, vec in vecs.items():
        n2 = sum(c * c * w_mono(j1, m1) * w_mono(j2, m2) for (m1, m2), c in vec.items())
        inv_norm = MQ.sqrt(1 / n2)
        for (m1, m2), c in vec.items():
            if c != 0:
                out[(m1, m2, M)] = MQ.sqrt(w_mono(j1, m1) * w_mono(j2, m2)) * inv_norm * (sgn * c)
    return out


def couple(x, j1, y, j2, L):
    """[x (x) y]_L components (index M+L), x, y lists of MQ (index m+j)."""
    out = [ZERO] * (2 * L + 1)
    for m1 in range(-j1, j1 + 1):
        a = x[m1 + j1]
        if a.is_zero():
            continue
        for m2 in range(max(-j2, -L - m1), min(j2, L - m1) + 1):
            b = y[m2 + j2]
            if b.is_zero():
                continue
            c = cg(j1, m1, j2, m2, L, m1 + m2)
            if c.is_zero():
                continue
            out[m1 + m2 + L] = out[m1 + m2 + L] + c * a * b
    return out


def theta(x, j):
    """(Theta x)_m = (-1)^m conj(x_{-m})."""
    return [(x[-m + j].conj() * ((-1) ** (m % 2))) for m in range(-j, j + 1)]


# ---------------------------------------------------------------- quaternions and SU(2)
def qmul(p, q):
    a1, b1, c1, d1 = p
    a2, b2, c2, d2 = q
    return (a1 * a2 - b1 * b2 - c1 * c2 - d1 * d2,
            a1 * b2 + b1 * a2 + c1 * d2 - d1 * c2,
            a1 * c2 - b1 * d2 + c1 * a2 + d1 * b2,
            a1 * d2 + b1 * c2 - c1 * b2 + d1 * a2)


def qinv(q):
    return (q[0], -q[1], -q[2], -q[3])


def qkey(q):
    return tuple(hash(x) for x in q)


def q_to_su2(q):
    """w+xi+yj+zk -> [[w - i z, -y - i x], [y - i x, w + i z]]  (i->-i sx, j->-i sy, k->-i sz)."""
    w, x, y, z = q
    I = MQ.I
    return [[w - I * z, -y - I * x], [y - I * x, w + I * z]]


def D_exact(j, g):
    """D^j(g) in the orthonormal CS basis, exact; g = 2x2 list of MQ; rows/cols index m+j."""
    (g11, g12), (g21, g22) = g
    n = 2 * j
    # powers
    def poly_pow(a, b, p):
        # (a x + b y)^p -> coeffs by power of x
        return {k: mq(comb(p, k)) * (a ** k) * (b ** (p - k)) for k in range(p + 1)}
    Dm = [[ZERO] * (n + 1) for _ in range(n + 1)]
    for m in range(-j, j + 1):
        P1 = poly_pow(g11, g21, j + m)   # g.x = g11 x + g21 y
        P2 = poly_pow(g12, g22, j - m)   # g.y = g12 x + g22 y
        for k1, c1 in P1.items():
            if c1.is_zero():
                continue
            for k2, c2 in P2.items():
                if c2.is_zero():
                    continue
                mp = k1 + k2 - j   # x power = j + m'
                Dm[mp + j][m + j] = Dm[mp + j][m + j] + c1 * c2
    D = [[ZERO] * (n + 1) for _ in range(n + 1)]
    for mp in range(-j, j + 1):
        for m in range(-j, j + 1):
            if not Dm[mp + j][m + j].is_zero():
                D[mp + j][m + j] = Dm[mp + j][m + j] * MQ.sqrt(w_mono(j, mp) / w_mono(j, m))
    return D


def spin_matrices(j):
    """Jz, J+, J- as float lists in the CS orthonormal basis (index m+j)."""
    import numpy as np
    n = 2 * j + 1
    Jz = np.diag([float(m) for m in range(-j, j + 1)]).astype(complex)
    Jp = np.zeros((n, n), complex)
    for m in range(-j, j):
        Jp[m + 1 + j, m + j] = ((j - m) * (j + m + 1)) ** 0.5
    Jm = Jp.T.copy()
    return Jz, Jp, Jm
