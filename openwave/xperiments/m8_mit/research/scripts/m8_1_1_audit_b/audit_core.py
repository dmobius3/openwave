#!/usr/bin/env python3
"""Audit core: exact arithmetic engine + group construction, built by a
deliberately different route from the audited script.

Differences from the audited implementation, on purpose:
  * cyclotomic field: reduction by a PRECOMPUTED power table, inverse by an
    exact linear solve against the multiplication matrix (not extended Euclid);
    the cyclotomic polynomial itself is built by dividing x^N - 1 by the
    lower Phi_d, not taken from sympy.
  * groups: C_n and BD_n are ABSTRACT (index arithmetic), with the SU(2)
    matrix realisation verified separately; 2T/2O/2I are UNIT QUATERNIONS over
    a real quadratic field Q(sqrt d), closed by multiplication.
  * conjugacy classes: brute force conjugation by EVERY element.
  * character tables: induced characters from every cyclic subgroup + Sym^n
    restrictions + orthogonality-driven peeling (no Burnside-Dixon, no mod p).
"""

from fractions import Fraction as Fr
from math import gcd
import itertools

# ---------------------------------------------------------------- arithmetic


def divisors(n):
    return sorted(d for d in range(1, n + 1) if n % d == 0)


def mobius(n):
    if n == 1:
        return 1
    res = 1
    d = 2
    m = n
    while d * d <= m:
        if m % d == 0:
            m //= d
            if m % d == 0:
                return 0
            res = -res
        d += 1
    if m > 1:
        res = -res
    return res


def ramanujan(m, k):
    """c_m(k) = sum_{gcd(a,m)=1} zeta_m^{ka}  =  sum_{d | gcd(k,m)} d*mu(m/d)."""
    g = gcd(k % m, m)
    return sum(d * mobius(m // d) for d in divisors(g))


def Tsum(m, j):
    """T(m,j) = sum_{a=1}^{m-1} zeta_m^{ja} / (2 - 2 cos(2 pi a/m))
              = (m^2-1)/12 - j(m-j)/2   for 0 <= j <= m-1.  Exact rational.
    Proved in the report: second difference and total sum both match."""
    j %= m
    return Fr(m * m - 1, 12) - Fr(j * (m - j), 2)


def prim_T(m, j):
    """Same sum restricted to a coprime to m (primitive m-th roots)."""
    j %= m
    tot = Tsum(m, j)
    for d in divisors(m):
        if 1 < d < m:
            tot -= prim_T(d, j % d)
    return tot


# ----------------------------------------------------------- integer polys
def poly_mul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    out[i + j] += x * y
    return out


def poly_divexact(a, b):
    """Exact division of integer polynomials (b monic-ish, leading divides)."""
    a = a[:]
    q = [0] * (len(a) - len(b) + 1)
    for i in range(len(a) - len(b), -1, -1):
        c = a[i + len(b) - 1]
        if c == 0:
            continue
        assert c % b[-1] == 0
        c //= b[-1]
        q[i] = c
        for k in range(len(b)):
            a[i + k] -= c * b[k]
    assert all(x == 0 for x in a), 'inexact polynomial division'
    return q


_CYC = {}


def cyclotomic(n):
    """Phi_n(x) as an integer coefficient list, low degree first.
    Built as (x^n - 1) / prod_{d|n, d<n} Phi_d.  No library table."""
    if n in _CYC:
        return _CYC[n]
    num = [-1] + [0] * (n - 1) + [1]
    den = [1]
    for d in divisors(n):
        if d < n:
            den = poly_mul(den, cyclotomic(d))
    res = poly_divexact(num, den)
    _CYC[n] = res
    return res


def poly_mod(a, m):
    """Reduce integer poly a modulo monic integer poly m."""
    a = a[:]
    dm = len(m) - 1
    for i in range(len(a) - 1, dm - 1, -1):
        c = a[i]
        if c:
            a[i] = 0
            base = i - dm
            for k in range(dm):
                a[base + k] -= c * m[k]
    a = a[:dm]
    while len(a) < dm:
        a.append(0)
    return a


class CycloField:
    """Q(zeta_N) in the power basis, table-driven reduction, inverse by solve."""

    def __init__(self, N):
        self.N = N
        self.phi = cyclotomic(N)
        self.deg = len(self.phi) - 1
        d = self.deg
        # x^k reduced, for k = 0 .. 2d-2  (needed by convolution reduction)
        self.red = []
        for k in range(max(2 * d - 1, N)):
            e = [0] * (k + 1)
            e[k] = 1
            self.red.append(poly_mod(e, self.phi))
        self.ZERO = tuple(Fr(0) for _ in range(d))
        self.ONE = self.rat(1)
        self._inv = {}

    def rat(self, q):
        q = Fr(q)
        v = [Fr(0)] * self.deg
        v[0] = q
        return tuple(v)

    def zeta(self, k):
        k %= self.N
        return tuple(Fr(c) for c in self.red[k])

    def add(self, a, b):
        return tuple(x + y for x, y in zip(a, b))

    def sub(self, a, b):
        return tuple(x - y for x, y in zip(a, b))

    def neg(self, a):
        return tuple(-x for x in a)

    def scal(self, a, q):
        q = Fr(q)
        return tuple(x * q for x in a)

    def mul(self, a, b):
        d = self.deg
        conv = [Fr(0)] * (2 * d - 1)
        for i in range(d):
            ai = a[i]
            if ai:
                for j in range(d):
                    bj = b[j]
                    if bj:
                        conv[i + j] += ai * bj
        out = [Fr(0)] * d
        for k in range(2 * d - 1):
            c = conv[k]
            if c:
                rk = self.red[k]
                for t in range(d):
                    if rk[t]:
                        out[t] += c * rk[t]
        return tuple(out)

    def inv(self, a):
        if a in self._inv:
            return self._inv[a]
        d = self.deg
        cols = []
        for i in range(d):
            e = [Fr(0)] * d
            e[i] = Fr(1)
            cols.append(self.mul(a, tuple(e)))
        M = [[cols[i][r] for i in range(d)] for r in range(d)]
        rhs = [Fr(1) if r == 0 else Fr(0) for r in range(d)]
        x = lin_solve(M, rhs)
        assert x is not None, 'not invertible'
        out = tuple(x)
        assert self.mul(a, out) == self.ONE
        self._inv[a] = out
        return out

    def galois(self, a, t):
        out = [Fr(0)] * self.deg
        for j in range(self.deg):
            if a[j]:
                zk = self.red[(j * t) % self.N]
                for s in range(self.deg):
                    if zk[s]:
                        out[s] += a[j] * zk[s]
        return tuple(out)

    def conj(self, a):
        return self.galois(a, self.N - 1)

    def is_rat(self, a):
        return all(x == 0 for x in a[1:])

    def to_rat(self, a):
        assert self.is_rat(a), f'not rational: {a}'
        return a[0]


def lin_solve(M, b):
    n = len(M)
    m = len(M[0])
    A = [[Fr(M[i][j]) for j in range(m)] + [Fr(b[i])] for i in range(n)]
    piv = []
    r = 0
    for c in range(m):
        pr = None
        for rr in range(r, n):
            if A[rr][c] != 0:
                pr = rr
                break
        if pr is None:
            continue
        A[r], A[pr] = A[pr], A[r]
        pv = A[r][c]
        A[r] = [v / pv for v in A[r]]
        for rr in range(n):
            if rr != r and A[rr][c] != 0:
                f = A[rr][c]
                A[rr] = [A[rr][k] - f * A[r][k] for k in range(m + 1)]
        piv.append(c)
        r += 1
        if r == n:
            break
    for row in A:
        if all(row[j] == 0 for j in range(m)) and row[m] != 0:
            return None
    x = [Fr(0)] * m
    for i, c in enumerate(piv):
        x[c] = A[i][m]
    return x


def null_space(M, m):
    n = len(M)
    A = [[Fr(M[i][j]) for j in range(m)] for i in range(n)]
    piv = []
    r = 0
    for c in range(m):
        pr = None
        for rr in range(r, n):
            if A[rr][c] != 0:
                pr = rr
                break
        if pr is None:
            continue
        A[r], A[pr] = A[pr], A[r]
        pv = A[r][c]
        A[r] = [v / pv for v in A[r]]
        for rr in range(n):
            if rr != r and A[rr][c] != 0:
                f = A[rr][c]
                A[rr] = [A[rr][k] - f * A[r][k] for k in range(m)]
        piv.append(c)
        r += 1
        if r == n:
            break
    free = [c for c in range(m) if c not in piv]
    out = []
    for fc in free:
        v = [Fr(0)] * m
        v[fc] = Fr(1)
        for i, c in enumerate(piv):
            v[c] = -A[i][fc]
        out.append(v)
    return out


def mat_rank(M, m):
    if not M:
        return 0
    return m - len(null_space(M, m))


# --------------------------------------------------- real quadratic field
class Quad:
    """a + b*sqrt(d), a,b in Q.  d squarefree positive (d=1 -> plain Q)."""

    def __init__(self, d):
        self.d = d
        self.ZERO = (Fr(0), Fr(0))
        self.ONE = (Fr(1), Fr(0))

    def rat(self, q):
        return (Fr(q), Fr(0))

    def add(self, x, y):
        return (x[0] + y[0], x[1] + y[1])

    def sub(self, x, y):
        return (x[0] - y[0], x[1] - y[1])

    def neg(self, x):
        return (-x[0], -x[1])

    def mul(self, x, y):
        return (x[0] * y[0] + self.d * x[1] * y[1], x[0] * y[1] + x[1] * y[0])

    def scal(self, x, q):
        q = Fr(q)
        return (x[0] * q, x[1] * q)


def qmul(Q, p, q):
    """Hamilton product of quaternions with Quad coefficients."""
    a1, b1, c1, d1 = p
    a2, b2, c2, d2 = q
    A = Q.sub(Q.sub(Q.sub(Q.mul(a1, a2), Q.mul(b1, b2)), Q.mul(c1, c2)), Q.mul(d1, d2))
    B = Q.sub(Q.add(Q.add(Q.mul(a1, b2), Q.mul(b1, a2)), Q.mul(c1, d2)), Q.mul(d1, c2))
    C = Q.add(Q.add(Q.sub(Q.mul(a1, c2), Q.mul(b1, d2)), Q.mul(c1, a2)), Q.mul(d1, b2))
    D = Q.add(Q.sub(Q.add(Q.mul(a1, d2), Q.mul(b1, c2)), Q.mul(c1, b2)), Q.mul(d1, a2))
    return (A, B, C, D)


def qnorm(Q, p):
    tot = Q.ZERO
    for x in p:
        tot = Q.add(tot, Q.mul(x, x))
    return tot
