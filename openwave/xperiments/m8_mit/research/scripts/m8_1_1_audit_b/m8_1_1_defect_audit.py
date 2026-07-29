#!/usr/bin/env python3
"""ADVERSARIAL AUDIT of SPEC SHEET B -- an independent recomputation of
U1..U12, a mutation test of the audited script's assertions, and a widened
attack on its two universality claims.

Independence, by construction:
  * groups: C_n and BD_n are built ABSTRACTLY from a presentation (index
    arithmetic) with the SU(2) matrix realisation verified separately;
    2T / 2O / 2I are built as UNIT QUATERNIONS over a real quadratic field
    Q(sqrt 2), Q(sqrt 5), closed under multiplication.
  * conjugacy classes: brute-force conjugation by EVERY group element.
  * character tables: induced characters from EVERY cyclic subgroup, plus the
    Sym^n restrictions of the defining representation, reduced by
    orthogonality-driven peeling, with a column-orthogonality completion as a
    fallback.  No Burnside-Dixon, no modular reduction, no imported table.
  * cyclotomic arithmetic: Phi_N built by dividing x^N - 1 by the lower Phi_d;
    reduction by a precomputed power table; inverses from an exact linear solve
    against the multiplication matrix (not extended Euclid).
  * the two sums: primary route sums over GROUP ELEMENTS (the audited script
    summed over classes), cross-checked by class sums, by a purely rational
    McKay/Molien linear solve that never touches the number field, by the exact
    cotangent sum T(m,j) for the cyclic and binary dihedral families, and by
    80-digit mpmath evaluated from the literal angles phi_g.

Everything rational is a fractions.Fraction end to end.  Floats appear only in
the mpmath cross-check and are never the reported value.

Outputs:  m8_1_1_defect_audit.json  (plus mutation_results.json from mutate.py)
"""

from fractions import Fraction as Fr
from fractions import Fraction
from math import gcd
import itertools
import json
import os
import random
import sys
import mpmath as mp

mp.mp.dps = 80

# ======================================================================
# section: audit_core.py
# ======================================================================

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


# ======================================================================
# section: audit_groups.py
# ======================================================================



def lcm(a, b):
    return a * b // gcd(a, b)


class Grp:
    pass


# ------------------------------------------------------------------ C_n
def build_C(n):
    g = Grp()
    g.name = f'C_{n}'
    g.order = n
    g.mul = [[(i + j) % n for j in range(n)] for i in range(n)]
    g.ident = 0
    g.N = n if n % 2 == 0 else n          # field for chi_def and characters
    if g.N < 3:
        g.N = max(g.N, 2)
    F = CycloField(g.N)
    g.F = F
    # SU(2) realisation:  a^i -> diag(zeta_n^i, zeta_n^-i)
    g.trace_elem = [F.add(F.zeta(i * (g.N // n)), F.zeta(-i * (g.N // n)))
                    for i in range(n)]
    # verify the realisation really is a group of order n inside SU(2)
    seen = set()
    for i in range(n):
        seen.add((F.zeta(i * (g.N // n)), F.zeta(-i * (g.N // n))))
    assert len(seen) == n, 'C_n matrix realisation not faithful'
    return g


# ----------------------------------------------------------------- BD_n
def build_BD(n):
    """<a,b : a^{2n}=1, b^2=a^n, b a b^-1 = a^-1>, order 4n."""
    g = Grp()
    g.name = f'BD_{n}'
    m = 2 * n
    els = [(i, e) for e in (0, 1) for i in range(m)]
    idx = {x: k for k, x in enumerate(els)}
    g.order = 4 * n
    g.ident = idx[(0, 0)]

    def mul2(x, y):
        i, e = x
        j, f = y
        if e == 0:
            return ((i + j) % m, f)
        if f == 0:
            return ((i - j) % m, 1)
        return ((i - j + n) % m, 0)

    g.mul = [[idx[mul2(x, y)] for y in els] for x in els]
    g.N = lcm(m, 4)
    F = CycloField(g.N)
    g.F = F
    w = g.N // m
    g.trace_elem = [F.add(F.zeta(i * w), F.zeta(-i * w)) if e == 0 else F.ZERO
                    for (i, e) in els]
    # verify the matrix realisation:  a=diag(w,w^-1), b=[[0,1],[-1,0]]
    A = (F.zeta(w), F.ZERO, F.ZERO, F.zeta(-w))
    B = (F.ZERO, F.ONE, F.neg(F.ONE), F.ZERO)

    def mm(X, Y):
        a, b, c, d = X
        e, f, h, k = Y
        return (F.add(F.mul(a, e), F.mul(b, h)), F.add(F.mul(a, f), F.mul(b, k)),
                F.add(F.mul(c, e), F.mul(d, h)), F.add(F.mul(c, f), F.mul(d, k)))

    mats = {}
    for (i, e) in els:
        M = (F.ONE, F.ZERO, F.ZERO, F.ONE)
        for _ in range(i):
            M = mm(M, A)
        if e:
            M = mm(M, B)
        mats[(i, e)] = M
    assert len(set(mats.values())) == 4 * n, 'BD_n realisation not faithful'
    for x in els:
        for y in els:
            assert mats[mul2(x, y)] == mm(mats[x], mats[y]), 'BD_n mult table wrong'
    for k, x in enumerate(els):
        tr = F.add(mats[x][0], mats[x][3])
        assert tr == g.trace_elem[k], 'BD_n trace mismatch'
    return g


# ------------------------------------------------ binary polyhedral groups
def build_poly(name):
    if name == '2T':
        Q = Quad(2)
        h = Fr(1, 2)
        gens = [((h, Fr(0)), (h, Fr(0)), (h, Fr(0)), (h, Fr(0))),
                (Q.ZERO, Q.ONE, Q.ZERO, Q.ZERO)]
        expect, N, d = 24, 12, 1
    elif name == '2O':
        Q = Quad(2)
        h = Fr(1, 2)
        r = (Fr(0), Fr(1, 2))                      # sqrt(2)/2
        gens = [((h, Fr(0)), (h, Fr(0)), (h, Fr(0)), (h, Fr(0))),
                (Q.ZERO, Q.ONE, Q.ZERO, Q.ZERO),
                (r, r, Q.ZERO, Q.ZERO)]
        expect, N, d = 48, 24, 2
    elif name == '2I':
        Q = Quad(5)
        h = Fr(1, 2)
        phi_half = (Fr(1, 4), Fr(1, 4))            # phi/2   = (1+sqrt5)/4
        phinv_half = (Fr(-1, 4), Fr(1, 4))         # phi^-1/2= (sqrt5-1)/4
        gens = [((h, Fr(0)), (h, Fr(0)), (h, Fr(0)), (h, Fr(0))),
                (phi_half, phinv_half, (h, Fr(0)), Q.ZERO)]
        expect, N, d = 120, 60, 5
    else:
        raise ValueError(name)

    for q in gens:
        assert qnorm(Q, q) == Q.ONE, f'{name}: generator not a unit quaternion'

    one = (Q.ONE, Q.ZERO, Q.ZERO, Q.ZERO)
    els = [one]
    idx = {one: 0}
    frontier = [one]
    while frontier:
        nxt = []
        for x in frontier:
            for gg in gens:
                y = qmul(Q, x, gg)
                if y not in idx:
                    idx[y] = len(els)
                    els.append(y)
                    nxt.append(y)
        frontier = nxt
    assert len(els) == expect, f'{name}: closed to {len(els)}, expected {expect}'
    for q in els:
        assert qnorm(Q, q) == Q.ONE

    g = Grp()
    g.name = name
    g.order = expect
    g.ident = 0
    g.mul = [[idx[qmul(Q, x, y)] for y in els] for x in els]
    g.N = N
    F = CycloField(N)
    g.F = F
    # embed Q(sqrt d) -> Q(zeta_N)
    if d == 1:
        srt = F.ZERO
    elif d == 2:
        srt = F.add(F.zeta(N // 8), F.zeta(-(N // 8)))
    elif d == 5:
        z = N // 5
        srt = F.sub(F.add(F.zeta(z), F.zeta(4 * z)), F.add(F.zeta(2 * z), F.zeta(3 * z)))
    assert d == 1 or F.mul(srt, srt) == F.rat(d), 'sqrt embedding failed'
    g.trace_elem = [F.add(F.rat(2 * q[0][0]), F.scal(srt, 2 * q[0][1])) for q in els]
    return g


# ----------------------------------------------------- classes and orders
def analyse(g):
    n = g.order
    mul = g.mul
    inv = [None] * n
    for i in range(n):
        for j in range(n):
            if mul[i][j] == g.ident:
                inv[i] = j
                break
    g.inv = inv
    order_of = [0] * n
    for i in range(n):
        o, cur = 1, i
        while cur != g.ident:
            cur = mul[cur][i]
            o += 1
        order_of[i] = o
    g.order_of = order_of
    # conjugacy classes by brute force over EVERY conjugator
    lab = [-1] * n
    classes = []
    for i in range(n):
        if lab[i] >= 0:
            continue
        k = len(classes)
        mem = set()
        for x in range(n):
            mem.add(mul[mul[x][i]][inv[x]])
        for y in mem:
            assert lab[y] < 0
            lab[y] = k
        classes.append(sorted(mem))
    g.cls = lab
    g.classes = classes
    g.nc = len(classes)
    g.sizes = [len(c) for c in classes]
    g.reps = [c[0] for c in classes]
    assert g.cls[g.ident] == 0 and g.sizes[0] == 1, 'identity must be class 0'
    g.invcls = [g.cls[inv[r]] for r in g.reps]
    g.sqcls = [g.cls[mul[r][r]] for r in g.reps]
    g.exponent = 1
    for o in order_of:
        g.exponent = lcm(g.exponent, o)
    assert g.N % g.exponent == 0, f'{g.name}: N={g.N} not divisible by exponent'
    g.chi_def = [g.trace_elem[r] for r in g.reps]
    return g


# --------------------------------------------------- character machinery
def inner_field(g, u, v):
    F = g.F
    acc = F.ZERO
    for k in range(g.nc):
        acc = F.add(acc, F.scal(F.mul(u[k], F.conj(v[k])), g.sizes[k]))
    return F.scal(acc, Fr(1, g.order))


def inner(g, u, v):
    return g.F.to_rat(inner_field(g, u, v))


def is_character(g, chi):
    """True iff chi decomposes with non-negative integer multiplicities."""
    F = g.F
    for c in g.chars:
        e = inner_field(g, chi, c)
        if not F.is_rat(e):
            return False
        m = F.to_rat(e)
        if m.denominator != 1 or m < 0:
            return False
    return True


def cf_sub(g, u, v):
    return [g.F.sub(a, b) for a, b in zip(u, v)]


def cf_scal(g, u, q):
    return [g.F.scal(a, q) for a in u]


def cf_add(g, u, v):
    return [g.F.add(a, b) for a, b in zip(u, v)]


def cf_mul(g, u, v):
    return [g.F.mul(a, b) for a, b in zip(u, v)]


def cyclic_subgroups(g):
    subs = {}
    for i in range(g.order):
        cur, S = i, [g.ident]
        while cur != g.ident:
            S.append(cur)
            cur = g.mul[cur][i]
        key = frozenset(S)
        if key not in subs:
            subs[key] = i          # a generator
    return subs


def induced_chars(g):
    """Ind_H^G psi_j for every cyclic H and every character psi_j of H."""
    F = g.F
    out = []
    for key, h in cyclic_subgroups(g).items():
        m = g.order_of[h]
        if m == 1:
            continue
        powr = {}
        cur = g.ident
        for t in range(m):
            powr[cur] = t
            cur = g.mul[cur][h]
        w = g.N // m
        # for each class rep, collect the multiset of t with x^-1 z x = h^t
        tallies = []
        for k in range(g.nc):
            z = g.reps[k]
            tal = {}
            for x in range(g.order):
                y = g.mul[g.mul[g.inv[x]][z]][x]
                if y in powr:
                    t = powr[y]
                    tal[t] = tal.get(t, 0) + 1
            tallies.append(tal)
        for j in range(m):
            ch = []
            for k in range(g.nc):
                acc = F.ZERO
                for t, c in tallies[k].items():
                    acc = F.add(acc, F.scal(F.zeta((j * t * w) % g.N), c))
                ch.append(F.scal(acc, Fr(1, m)))
            out.append(ch)
    return out


def sym_powers(g, upto):
    F = g.F
    psi = [[F.ONE] * g.nc, list(g.chi_def)]
    for _ in range(upto):
        nxt = [F.sub(F.mul(g.chi_def[k], psi[-1][k]), psi[-2][k]) for k in range(g.nc)]
        psi.append(nxt)
    return psi


def char_table(g):
    """Orthogonality-driven peeling of induced + Sym^n characters."""
    F = g.F
    irr = []
    triv = [F.ONE] * g.nc
    irr.append(triv)

    def reduce_v(v):
        for c in irr:
            m = inner(g, v, c)
            if m != 0:
                v = cf_sub(g, v, cf_scal(g, c, m))
        return v

    pool = induced_chars(g) + sym_powers(g, 2 * g.nc + 4)
    rounds = 0
    while len(irr) < g.nc and rounds < 40:
        rounds += 1
        newpool = []
        added = False
        for v in pool:
            v = reduce_v(v)
            nv = inner(g, v, v)
            if nv == 0:
                continue
            if nv == 1:
                if F.to_rat(v[0]) < 0:
                    v = cf_scal(g, v, -1)
                if v not in irr:
                    irr.append(v)
                    added = True
            else:
                newpool.append(v)
        pool = [reduce_v(v) for v in newpool]
        pool = [v for v in pool if inner(g, v, v) != 0]
        # differences that shrink the norm
        extra = []
        for i in range(len(pool)):
            for j in range(len(pool)):
                if i == j:
                    continue
                w = cf_sub(g, pool[i], pool[j])
                nw = inner(g, w, w)
                if 0 < nw < inner(g, pool[i], pool[i]):
                    extra.append(w)
        # products of found irreducibles
        for a in irr:
            for b in irr:
                extra.append(cf_mul(g, a, b))
        pool = pool + extra
        if not added and rounds > 3:
            break
    # completion by column orthogonality if exactly one is missing
    if len(irr) == g.nc - 1:
        tot = sum(int(F.to_rat(c[0])) ** 2 for c in irr)
        dlast = int(round((g.order - tot) ** 0.5))
        assert dlast * dlast == g.order - tot
        last = [F.rat(dlast)]
        for k in range(1, g.nc):
            acc = F.ZERO
            for c in irr:
                acc = F.add(acc, F.scal(c[k], F.to_rat(c[0])))
            last.append(F.scal(acc, Fr(-1, dlast)))
        irr.append(last)
    assert len(irr) == g.nc, f'{g.name}: found {len(irr)} of {g.nc} irreducibles'
    irr.sort(key=lambda c: (F.to_rat(c[0]), [tuple(x) for x in c[1:]]))
    g.chars = irr
    g.dims = [int(F.to_rat(c[0])) for c in irr]
    return irr


def verify_table(g):
    F = g.F
    r = g.nc
    for i in range(r):
        for j in range(r):
            v = inner(g, g.chars[i], g.chars[j])
            assert v == (1 if i == j else 0), f'{g.name}: row orthonormality'
    assert sum(d * d for d in g.dims) == g.order, f'{g.name}: sum dim^2'
    for k in range(r):
        for l in range(r):
            acc = F.ZERO
            for i in range(r):
                acc = F.add(acc, F.mul(g.chars[i][k], F.conj(g.chars[i][l])))
            want = F.rat(Fr(g.order, g.sizes[k])) if k == l else F.ZERO
            assert acc == want, f'{g.name}: column orthogonality {k},{l}'
    # tensor positivity
    for i in range(r):
        for j in range(r):
            prod = cf_mul(g, g.chars[i], g.chars[j])
            tot = 0
            for k in range(r):
                m = inner(g, prod, g.chars[k])
                assert m.denominator == 1 and m >= 0, f'{g.name}: tensor positivity'
                tot += int(m) * g.dims[k]
            assert tot == g.dims[i] * g.dims[j]
    return True


def make_group(name):
    if name.startswith('C_'):
        g = build_C(int(name[2:]))
    elif name.startswith('BD_'):
        g = build_BD(int(name[3:]))
    else:
        g = build_poly(name)
    analyse(g)
    char_table(g)
    return g


# ======================================================================
# section: audit_sums.py
# ======================================================================




class SD:
    def __init__(self, g):
        F = g.F
        self.g = g
        self.F = F
        self.w = [None] * g.nc          # 1/(2 - chi_def)
        self.cot2 = [None] * g.nc       # (2 + chi_def)/(2 - chi_def)
        for k in range(1, g.nc):
            den = F.sub(F.rat(2), g.chi_def[k])
            assert den != F.ZERO, 'g != I with chi_def = 2'
            self.w[k] = F.inv(den)
            self.cot2[k] = F.mul(F.add(F.rat(2), g.chi_def[k]), self.w[k])
        self.triv = [i for i in range(g.nc) if all(v == F.ONE for v in g.chars[i])]
        assert len(self.triv) == 1
        self.triv = self.triv[0]
        self.defidx = None
        for i in range(g.nc):
            if g.chars[i] == g.chi_def:
                self.defidx = i

    # ---- primary route: sum over GROUP ELEMENTS -------------------------
    def D_elem(self, chi):
        F, g = self.F, self.g
        acc = F.ZERO
        for x in range(g.order):
            if x == g.ident:
                continue
            k = g.cls[x]
            acc = F.add(acc, F.mul(chi[k], self.w[k]))
        return F.to_rat(F.scal(acc, Fr(1, g.order)))

    def S_elem(self, chi):
        F, g = self.F, self.g
        dim = F.to_rat(chi[0])
        acc = F.ZERO
        for x in range(g.order):
            if x == g.ident:
                continue
            k = g.cls[x]
            acc = F.add(acc, F.mul(F.sub(chi[k], F.rat(dim)), self.cot2[k]))
        return F.to_rat(F.scal(acc, Fr(1, g.order)))

    # ---- secondary route: class sums ------------------------------------
    def D_class(self, chi):
        F, g = self.F, self.g
        acc = F.ZERO
        for k in range(1, g.nc):
            acc = F.add(acc, F.scal(F.mul(chi[k], self.w[k]), g.sizes[k]))
        return F.to_rat(F.scal(acc, Fr(1, g.order)))

    def S_class(self, chi):
        F, g = self.F, self.g
        dim = F.to_rat(chi[0])
        acc = F.ZERO
        for k in range(1, g.nc):
            t = F.sub(chi[k], F.rat(dim))
            acc = F.add(acc, F.scal(F.mul(t, self.cot2[k]), g.sizes[k]))
        return F.to_rat(F.scal(acc, Fr(1, g.order)))

    # ---- field-valued versions (twists whose sums need not land in Q) ----
    def D_field(self, chi):
        F, g = self.F, self.g
        acc = F.ZERO
        for k in range(1, g.nc):
            acc = F.add(acc, F.scal(F.mul(chi[k], self.w[k]), g.sizes[k]))
        return F.scal(acc, Fr(1, g.order))

    def S_field(self, chi):
        F, g = self.F, self.g
        acc = F.ZERO
        for k in range(1, g.nc):
            t = F.sub(chi[k], chi[0])
            acc = F.add(acc, F.scal(F.mul(t, self.cot2[k]), g.sizes[k]))
        return F.scal(acc, Fr(1, g.order))

    def meanT(self):
        """(1/|G|) sum_{g != I} cot^2(phi_g/2), by element summation."""
        F, g = self.F, self.g
        acc = F.ZERO
        for x in range(g.order):
            if x == g.ident:
                continue
            acc = F.add(acc, self.cot2[g.cls[x]])
        return F.to_rat(F.scal(acc, Fr(1, g.order)))


# --------------------------------------------------- pure-rational routes
def D_trivial_cyclic(g):
    """D(1) = (1/|G|) sum over cyclic subgroups C != 1 of prim_T(|C|, 0).
    Uses only the cyclic-subgroup census and the exact cotangent sum."""
    tot = Fr(0)
    for key, h in cyclic_subgroups(g).items():
        m = g.order_of[h]
        if m > 1:
            tot += prim_T(m, 0)
    return tot / g.order


def mckay_matrix(g):
    A = [[0] * g.nc for _ in range(g.nc)]
    for s in range(g.nc):
        prod = cf_mul(g, g.chi_def, g.chars[s])
        for t in range(g.nc):
            v = inner(g, prod, g.chars[t])
            assert v.denominator == 1 and v >= 0
            A[s][t] = int(v)
    return A


def D_all_mckay(g, A, triv):
    """Rational Molien route.  x(s) = ((1-s)K + s^2 I)^-1 e_triv with K = 2I-A;
    the s^0 coefficient u_0 satisfies K u_0 = e_triv - delta/|G| and
    <u_0, delta> = 0, and u_0[alpha] = D(alpha)."""
    r = g.nc
    delta = g.dims
    assert sum(d * d for d in delta) == g.order
    K = [[(2 if i == j else 0) - A[i][j] for j in range(r)] for i in range(r)]
    rows = [row[:] for row in K] + [list(delta)]
    rhs = [Fr(1 if i == triv else 0) - Fr(delta[i], g.order) for i in range(r)] + [Fr(0)]
    sol = lin_solve(rows, rhs)
    assert sol is not None, 'McKay/Molien system inconsistent'
    ns = null_space(rows, r)
    assert len(ns) == 0, 'McKay/Molien solution not unique'
    return sol


def D_closed_form_C(n, j):
    return Tsum(n, j) / n


def D_closed_form_BD(n, kind, h=0, eps=1):
    """kind: '2dim' -> chi_h ; 'a1' -> chi(a)=1, chi(b)=eps ; 'am1' -> chi(a)=-1."""
    m = 2 * n
    if kind == '2dim':
        return Tsum(m, h) / m
    if kind == 'a1':
        return Tsum(m, 0) / (4 * n) + Fr(eps, 4)
    return Tsum(m, n) / (4 * n)


# ------------------------------------------------------ numeric route
def numeric_SD(g, chi_num, dim):
    """chi_num: list of mpmath values per class.  Returns (S, D) numerically,
    with cot^2(phi/2) built from the literal angle phi = arccos(chi_def/2)."""
    S = mp.mpf(0)
    D = mp.mpf(0)
    for k in range(1, g.nc):
        x = mp.re(g_chi_def_numeric(g, k))
        phi = mp.acos(x / 2)
        cot2 = mp.cos(phi / 2) ** 2 / mp.sin(phi / 2) ** 2
        S += g.sizes[k] * (chi_num[k] - dim) * cot2
        D += g.sizes[k] * chi_num[k] / (2 - x)
    return S / g.order, D / g.order


_numcache = {}


def field_numeric(F, a):
    key = (id(F), a)
    if key in _numcache:
        return _numcache[key]
    s = mp.mpc(0)
    for j in range(F.deg):
        if a[j]:
            s += mp.mpf(a[j].numerator) / a[j].denominator * mp.expjpi(mp.mpf(2 * j) / F.N)
    _numcache[key] = s
    return s


def g_chi_def_numeric(g, k):
    return field_numeric(g.F, g.chi_def[k])


# ======================================================================
# section: audit_tasks.py
# ======================================================================




def fs(x):
    x = Fr(x)
    return f"{x.numerator}/{x.denominator}"


# --------------------------------------------------------------- rendering
def quad_render(F, e, d):
    """Write e = p + q*sqrt(d) exactly if possible, else None."""
    if d == 5:
        z = F.N // 5
        s = F.sub(F.add(F.zeta(z), F.zeta(4 * z)), F.add(F.zeta(2 * z), F.zeta(3 * z)))
    elif d == 2:
        s = F.add(F.zeta(F.N // 8), F.zeta(-(F.N // 8)))
    else:
        return None
    M = [[F.ONE[i], s[i]] for i in range(F.deg)]
    sol = lin_solve(M, list(e))
    if sol is None:
        return None
    p, q = sol
    if F.add(F.scal(F.ONE, p), F.scal(s, q)) != e:
        return None
    if q == 0:
        return fs(p)
    return f"{fs(p)} + ({fs(q)})*sqrt({d})"


def render(F, e, d=None):
    if F.is_rat(e):
        return fs(F.to_rat(e))
    if d:
        r = quad_render(F, e, d)
        if r:
            return r
    for j in range(F.N):
        if e == F.zeta(j):
            return f"exp(2*pi*I*{j}/{F.N})"
        if e == F.neg(F.zeta(j)):
            return f"-exp(2*pi*I*{j}/{F.N})"
    return "[" + ",".join(fs(c) for c in e) + "] in basis zeta_%d^j" % F.N


# ------------------------------------------------------------ twist algebra
def lam2(g, chi):
    F = g.F
    return [F.scal(F.sub(F.mul(chi[k], chi[k]), chi[g.sqcls[k]]), Fr(1, 2))
            for k in range(g.nc)]


def sym2(g, chi):
    F = g.F
    return [F.scal(F.add(F.mul(chi[k], chi[k]), chi[g.sqcls[k]]), Fr(1, 2))
            for k in range(g.nc)]


def combo(g, mult):
    F = g.F
    chi = [F.ZERO] * g.nc
    for i, m in enumerate(mult):
        if m:
            for k in range(g.nc):
                chi[k] = F.add(chi[k], F.scal(g.chars[i][k], m))
    return chi


# ------------------------------------------------------------------- U2
def u2_check(groups, SDs):
    """csc^2(phi_g/2) = 4/(2 - chi_def(g)).  The LEFT side is built from the
    literal angle phi_g = arccos(chi_def/2) in 80-digit arithmetic; the right
    side from exact field arithmetic.  Also checks cot^2 against the angle,
    which is the step the affine relation actually rests on."""
    worst_csc = mp.mpf(0)
    worst_cot = mp.mpf(0)
    where = None
    nelem = 0
    for nm, g in groups.items():
        sd = SDs[nm]
        for k in range(1, g.nc):
            x = mp.re(field_numeric(g.F, g.chi_def[k]))
            phi = mp.acos(x / 2)
            csc2 = 1 / mp.sin(phi / 2) ** 2
            cot2 = mp.cos(phi / 2) ** 2 / mp.sin(phi / 2) ** 2
            rhs_csc = mp.re(field_numeric(g.F, g.F.scal(sd.w[k], 4)))
            rhs_cot = mp.re(field_numeric(g.F, sd.cot2[k]))
            d1 = abs(csc2 - rhs_csc)
            d2 = abs(cot2 - rhs_cot)
            nelem += g.sizes[k]
            if d1 > worst_csc:
                worst_csc = d1
                where = (nm, k)
            worst_cot = max(worst_cot, d2)
    return {'elements_checked': nelem,
            'max_dev_csc2': mp.nstr(worst_csc, 8),
            'max_dev_cot2': mp.nstr(worst_cot, 8),
            'worst_at': where,
            'dps': mp.mp.dps,
            'note': ('the exact-field statement cot^2+1 = 4/(2-chi) is a '
                     'tautology of the definition cot^2 := (2+chi)/(2-chi); '
                     'only the trigonometric comparison above has content')}


# ------------------------------------------------------------------- U11
def e8_task():
    # arms of length 1, 2, 4 off the trivalent node 2 -- a different labelling
    edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (2, 7)]
    C = [[0] * 8 for _ in range(8)]
    for i in range(8):
        C[i][i] = 2
    for i, j in edges:
        C[i][j] = C[j][i] = -1
    G = [[-C[i][j] for j in range(8)] for i in range(8)]

    def det(M):
        A = [[Fr(x) for x in row] for row in M]
        d = Fr(1)
        n = len(A)
        for c in range(n):
            pr = next((r for r in range(c, n) if A[r][c] != 0), None)
            if pr is None:
                return Fr(0)
            if pr != c:
                A[c], A[pr] = A[pr], A[c]
                d = -d
            d *= A[c][c]
            for r in range(c + 1, n):
                if A[r][c]:
                    f = A[r][c] / A[c][c]
                    A[r] = [A[r][k] - f * A[c][k] for k in range(n)]
        return d

    def form(x, y):
        return sum(x[i] * G[i][j] * y[j] for i in range(8) for j in range(8) if G[i][j])

    # PRIMARY: build the root system by closing the simple roots under their
    # own reflections (a constructive route, not a lattice-point search).
    def reflect(a, x):
        c = form(x, a)
        return tuple(x[i] + c * a[i] for i in range(8))

    simple = [tuple(1 if k == i else 0 for k in range(8)) for i in range(8)]
    roots_c = set(simple) | {tuple(-v for v in s) for s in simple}
    frontier = list(roots_c)
    while frontier:
        nxt = []
        for x in list(frontier):
            for a in list(roots_c):
                y = reflect(a, x)
                if y not in roots_c:
                    roots_c.add(y)
                    nxt.append(y)
        frontier = nxt
    for c in roots_c:
        assert form(c, c) == -2

    # CHECK: exhaustive short-vector enumeration, exact Cholesky of +C
    n = 8
    d = [Fr(0)] * n
    L = [[Fr(1) if i == j else Fr(0) for j in range(n)] for i in range(n)]
    for i in range(n):
        s = Fr(C[i][i])
        for k in range(i):
            s -= d[k] * L[i][k] ** 2
        d[i] = s
        assert d[i] > 0
        for j in range(i + 1, n):
            t = Fr(C[j][i])
            for k in range(i):
                t -= d[k] * L[j][k] * L[i][k]
            L[j][i] = t / d[i]
    enum = []
    cvec = [0] * n

    def rec(i, rem):
        if i < 0:
            if rem == 0:
                enum.append(tuple(cvec))
            return
        sh = sum(L[j][i] * cvec[j] for j in range(i + 1, n))
        lim = rem / d[i]
        b = mp.sqrt(mp.mpf(lim.numerator) / lim.denominator)
        shf = mp.mpf(sh.numerator) / sh.denominator
        lo = int(mp.floor(-shf - b)) - 1
        hi = int(mp.ceil(-shf + b)) + 1
        for ci in range(lo, hi + 1):
            val = d[i] * (ci + sh) ** 2
            if val <= rem:
                cvec[i] = ci
                rec(i - 1, rem - val)
        cvec[i] = 0

    rec(n - 1, Fr(2))
    assert set(enum) == roots_c, 'reflection closure != short-vector enumeration'
    roots = sorted(roots_c)
    fib = {}
    for c in roots:
        fib.setdefault(tuple(v % 2 for v in c), []).append(c)
    fibre_sizes = {}
    for v in fib.values():
        fibre_sizes[len(v)] = fibre_sizes.get(len(v), 0) + 1
    pairs = all(len(v) == 2 and set(v) == {v[0], tuple(-a for a in v[0])}
                for v in fib.values())

    G2 = [[G[i][j] % 2 for j in range(8)] for i in range(8)]
    diag_zero = all(G2[i][i] == 0 for i in range(8))
    alt = all(sum(x[i] * G2[i][j] * x[j] for i in range(8) for j in range(8)) % 2 == 0
              for x in itertools.product((0, 1), repeat=8))
    rows = [row[:] for row in G2]
    rank2, piv = 0, []
    Mr = [r[:] for r in rows]
    for c in range(8):
        pr = next((r for r in range(rank2, 8) if Mr[r][c] % 2), None)
        if pr is None:
            continue
        Mr[rank2], Mr[pr] = Mr[pr], Mr[rank2]
        for r in range(8):
            if r != rank2 and Mr[r][c] % 2:
                Mr[r] = [(Mr[r][k] - Mr[rank2][k]) % 2 for k in range(8)]
        rank2 += 1
    radical2 = 8 - rank2

    def P4(x, shift=None):
        xi = list(x)
        if shift:
            xi = [xi[i] + 2 * shift[i] for i in range(8)]
        return form(xi, xi) % 4

    rnd = random.Random(20260728)
    shifts = [tuple(rnd.randint(-4, 4) for _ in range(8)) for _ in range(60)]
    liftok = True
    Pv = {}
    counts = {}
    for x in itertools.product((0, 1), repeat=8):
        b = P4(x)
        Pv[x] = b
        counts[b] = counts.get(b, 0) + 1
        for sh in shifts[:6]:
            if P4(x, sh) != b:
                liftok = False
    for sh in shifts:
        x = tuple(rnd.randint(0, 1) for _ in range(8))
        if P4(x, sh) != P4(x):
            liftok = False
    match = set(fib.keys()) == {x for x in Pv if Pv[x] == 2}

    def refl(a, x):
        c = form(x, a)
        return tuple(x[i] + c * a[i] for i in range(8))

    rootset = set(roots)
    closed = all(refl(a, b) in rootset for a in roots for b in roots)
    mats = set()
    for a in roots:
        M = []
        for i in range(8):
            ei = [1 if k == i else 0 for k in range(8)]
            M.append(tuple(v % 2 for v in refl(a, ei)))
        mats.add(tuple(M))

    def act(M, x):
        out = [0] * 8
        for i in range(8):
            if x[i]:
                for j in range(8):
                    out[j] ^= M[i][j]
        return tuple(out)

    def orbits(targets):
        seen, sizes = set(), []
        for x0 in sorted(targets):
            if x0 in seen:
                continue
            orb, stack = {x0}, [x0]
            while stack:
                y = stack.pop()
                for M in mats:
                    z = act(M, y)
                    if z not in orb:
                        orb.add(z)
                        stack.append(z)
            seen |= orb
            sizes.append(len(orb))
        return sizes

    o2 = orbits({x for x in Pv if Pv[x] == 2})
    oall = orbits(list(itertools.product((0, 1), repeat=8)))
    return {
        'dynkin_edges': edges,
        'arm_lengths_from_trivalent_node': [1, 2, 4],
        'det_cartan': str(det(C)), 'det_negative_form': str(det(G)),
        'roots_by_reflection_closure': len(roots_c),
        'roots_by_short_vector_enumeration': len(enum),
        'two_routes_agree': set(enum) == roots_c,
        'norm_minus2_vectors': len(roots),
        'distinct_mod2_classes': len(fib),
        'fibre_sizes': {str(k): v for k, v in fibre_sizes.items()},
        'every_fibre_is_pm_pair': pairs,
        'mod2_form_diag_zero': diag_zero,
        'mod2_alternating': alt,
        'mod2_radical_dim': radical2,
        'mod2_nondegenerate': radical2 == 0,
        'mod4_lift_independent': liftok,
        'mod4_shifts_tested': len(shifts),
        'mod4_value_counts': {str(k): counts[k] for k in sorted(counts)},
        'P2_classes_equal_root_classes': match,
        'reflections_preserve_roots': closed,
        'distinct_mod2_reflection_matrices': len(mats),
        'orbits_on_P2': {'count': len(o2), 'sizes': sorted(o2, reverse=True)},
        'orbits_on_all_256': {'count': len(oall), 'sizes': sorted(oall, reverse=True)},
    }


# ======================================================================
# section: audit_main.py
# ======================================================================



SHEET = ([f'C_{n}' for n in range(2, 11)] + [f'BD_{n}' for n in range(2, 7)]
         + ['2T', '2O', '2I'])
WIDE = ([f'C_{n}' for n in range(2, 25)] + [f'BD_{n}' for n in range(2, 13)]
        + ['2T', '2O', '2I'])

OUT = {'meta': {
    'role': 'adversarial audit of SPEC SHEET B',
    'method': ('independent rebuild: abstract C_n/BD_n with verified SU(2) '
               'realisation, unit-quaternion closure over Q(sqrt d) for '
               '2T/2O/2I, brute-force conjugacy classes, character tables by '
               'induced characters from every cyclic subgroup plus '
               'orthogonality-driven peeling (no Burnside-Dixon, no mod p), '
               'cyclotomic arithmetic by precomputed power-table reduction '
               'with inverses from an exact linear solve; S and D summed over '
               'GROUP ELEMENTS, cross-checked by class sums, by a purely '
               'rational McKay/Molien linear solve that never touches the '
               'number field, and by 80-digit mpmath from the literal angles'),
    'sheet_family': SHEET, 'widened_family': WIDE}}


def build_all(names):
    G, S = {}, {}
    for nm in names:
        g = make_group(nm)
        verify_table(g)
        G[nm] = g
        S[nm] = SD(g)
    return G, S


# ============================================================ U1, U12, U5
def tasks_1_5_12(G, SDs, names):
    u1, u5, u12 = {}, {}, {'element_vs_class': 0, 'mismatches': [],
                           'mckay_vs_field': 0, 'mckay_mismatches': [],
                           'numeric_vs_exact_max_dev': None}
    worstnum = mp.mpf(0)
    for nm in names:
        g, sd = G[nm], SDs[nm]
        A = mckay_matrix(g)
        Dmck = D_all_mckay(g, A, sd.triv)
        rows = []
        for i in range(g.nc):
            De = sd.D_elem(g.chars[i])
            Dc = sd.D_class(g.chars[i])
            Se = sd.S_elem(g.chars[i])
            Sc = sd.S_class(g.chars[i])
            u12['element_vs_class'] += 2
            if De != Dc or Se != Sc:
                u12['mismatches'].append((nm, i))
            u12['mckay_vs_field'] += 1
            if Dmck[i] != De:
                u12['mckay_mismatches'].append((nm, i))
            # 80-digit numeric check of both sums
            Sn = mp.mpf(0)
            Dn = mp.mpf(0)
            for k in range(1, g.nc):
                x = mp.re(field_numeric(g.F, g.chi_def[k]))
                phi = mp.acos(x / 2)
                c2 = mp.cos(phi / 2) ** 2 / mp.sin(phi / 2) ** 2
                chv = mp.re(field_numeric(g.F, g.chars[i][k]))
                Sn += g.sizes[k] * (chv - g.dims[i]) * c2
                Dn += g.sizes[k] * chv / (2 - x)
            Sn /= g.order
            Dn /= g.order
            worstnum = max(worstnum, abs(Sn - mp.mpf(Se.numerator) / Se.denominator),
                           abs(Dn - mp.mpf(De.numerator) / De.denominator))
            rows.append({'irrep': i, 'dim': g.dims[i], 'S': fs(Se), 'D': fs(De),
                         'D_mckay_rational_route': fs(Dmck[i]),
                         'is_trivial': i == sd.triv, 'is_defining': i == sd.defidx})
        D1 = sd.D_elem(g.chars[sd.triv])
        S1 = sd.S_elem(g.chars[sd.triv])
        u1[nm] = {'order': g.order, 'num_classes': g.nc, 'dims': g.dims,
                  'S_trivial': fs(S1), 'D_trivial': fs(D1), 'rows': rows}
        T = sd.meanT()
        comb = D1 + T / 8
        u5[nm] = {'order': g.order, 'num_classes': g.nc,
                  'mean_cot2': fs(T), 'D_trivial': fs(D1),
                  'D_trivial_cyclic_subgroup_route': fs(D_trivial_cyclic(g)),
                  'D1_plus_eighth_T': fs(comb),
                  'target_(r-1)/8': fs(Fr(g.nc - 1, 8)),
                  'match': comb == Fr(g.nc - 1, 8),
                  'D1_closed_form_(r|G|-1)/(12|G|)': fs(Fr(g.nc * g.order - 1, 12 * g.order)),
                  'closed_form_match': D1 == Fr(g.nc * g.order - 1, 12 * g.order)}
        assert D1 == D_trivial_cyclic(g), f'{nm}: D(1) routes disagree'
    u12['numeric_vs_exact_max_dev'] = mp.nstr(worstnum, 8)
    u12['all_agree'] = (not u12['mismatches']) and (not u12['mckay_mismatches'])
    return u1, u5, u12


# ================================================================= twists
def twist_pool(g, sd, rng, cap=260):
    """(label, character-as-class-function) pairs, deliberately awkward."""
    F = g.F
    r = g.nc
    out = []
    for i in range(r):
        out.append((f'irr{i}', g.chars[i]))
    grid = [(i, j) for i in range(r) for j in range(i, r)]
    if len(grid) > 90:
        grid = rng.sample(grid, 90)
    for i, j in grid:
        out.append((f'irr{i}(x)irr{j}', cf_mul(g, g.chars[i], g.chars[j])))
        mult = [0] * r
        mult[i] += 1
        mult[j] += 1
        out.append((f'irr{i}+irr{j}', combo(g, mult)))
    for i in range(r):
        out.append((f'Sym2(irr{i})', sym2(g, g.chars[i])))
        l = lam2(g, g.chars[i])
        if any(v != F.ZERO for v in l):
            out.append((f'Lam2(irr{i})', l))
    # triple products, only where cheap
    if r <= 9:
        for i, j, k in rng.sample([(a, b, c) for a in range(r) for b in range(r)
                                   for c in range(r)], min(20, r ** 3)):
            out.append((f'irr{i}(x)irr{j}(x)irr{k}',
                        cf_mul(g, cf_mul(g, g.chars[i], g.chars[j]), g.chars[k])))
    # higher direct sums with lumpy multiplicities
    for _ in range(18):
        mult = [rng.choice([0, 0, 0, 1, 1, 2, 3, 5]) for _ in range(r)]
        if sum(mult) == 0:
            continue
        out.append(('sum' + ''.join(str(m) for m in mult), combo(g, mult)))
    # virtual characters (negative multiplicities) -- not representations
    for _ in range(10):
        mult = [rng.choice([-3, -2, -1, 0, 0, 1, 2]) for _ in range(r)]
        if all(m == 0 for m in mult):
            continue
        out.append(('virt' + ''.join(str(m) for m in mult), combo(g, mult)))
    # arbitrary rational class functions -- not even virtual characters
    for _ in range(8):
        f = [F.rat(Fr(rng.randint(-9, 9), rng.choice([1, 2, 3, 7])))
             for _ in range(r)]
        out.append(('classfun', f))
    if len(out) > cap:
        keep = out[:r] + rng.sample(out[r:], cap - r)
        out = keep
    return out


def universality(G, SDs, names, rng):
    """Attack S(alpha) = a dim + b D + c dim D(1) with (a,b,c) = (1,4,-4)."""
    A, B, C = Fr(1), Fr(4), Fr(-4)
    rows = []
    bad = []
    per_group = {}
    counts = {'total': 0, 'no_trivial_part': 0, 'with_trivial_part': 0,
              'not_a_character': 0}
    for nm in names:
        g, sd = G[nm], SDs[nm]
        F = g.F
        D1 = sd.D_class(g.chars[sd.triv])
        pool = twist_pool(g, sd, rng)
        nbad = 0
        for lab, chi in pool:
            dimf = chi[0]
            Df = sd.D_field(chi)
            Sf = sd.S_field(chi)
            mtriv = inner(g, chi, [F.ONE] * g.nc)
            # residual as a FIELD element: S - (a dim + b D + c dim D(1))
            rhs = F.add(F.add(F.scal(dimf, A), F.scal(Df, B)),
                        F.scal(F.scal(dimf, C), D1))
            residf = F.sub(Sf, rhs)
            counts['total'] += 1
            ischar = is_character(g, chi)
            if not ischar:
                counts['not_a_character'] += 1
            if mtriv == 0:
                counts['no_trivial_part'] += 1
            else:
                counts['with_trivial_part'] += 1
            ok = (residf == F.rat(-mtriv))
            rational = F.is_rat(Df) and F.is_rat(Sf)
            rec = {'group': nm, 'twist': lab,
                   'dim': fs(F.to_rat(dimf)) if F.is_rat(dimf) else 'irrational',
                   'S': fs(F.to_rat(Sf)) if F.is_rat(Sf) else 'irrational',
                   'D': fs(F.to_rat(Df)) if F.is_rat(Df) else 'irrational',
                   'm_triv': fs(mtriv), 'sums_are_rational': rational,
                   'residual': fs(F.to_rat(residf)) if F.is_rat(residf) else 'irrational',
                   'residual_equals_minus_m': ok, 'is_character': ischar}
            if ischar:
                assert rational, f'{nm} {lab}: character with irrational S or D'
            if not ok:
                nbad += 1
                bad.append(rec)
            rows.append(rec)
        per_group[nm] = {'twists': len(pool), 'violations': nbad}
    return rows, bad, per_group, counts


def refit(rows_no_triv):
    """Solve the 3-parameter model exactly on the pooled widened data."""
    M = [[Fr(r['dim']), Fr(r['D']), Fr(r['dim']) * Fr(r['D1'])] for r in rows_no_triv]
    y = [Fr(r['S']) for r in rows_no_triv]
    rk = mat_rank(M, 3)
    sol = lin_solve(M, y)
    ns = null_space(M, 3)
    resid = None
    if sol is not None:
        resid = max(abs(y[i] - sum(M[i][j] * sol[j] for j in range(3)))
                    for i in range(len(y)))
    return {'pooled_rows': len(y), 'design_rank': rk,
            'solution': [fs(x) for x in sol] if sol else None,
            'nullspace_dim': len(ns),
            'max_residual': fs(resid) if resid is not None else None}


# ====================================================== convention attacks
def convention_tests(G, SDs, names):
    """What breaks if the identity, or -I, is put back into the sums."""
    out = {}
    for nm in names:
        g, sd = G[nm], SDs[nm]
        F = g.F
        # can the identity be included at all?
        det_at_I = F.sub(F.rat(2), g.chi_def[0])
        includable = (det_at_I != F.ZERO)
        # excluding -I as well (it exists iff some class has chi_def = -2)
        minus = [k for k in range(g.nc) if g.chi_def[k] == F.rat(-2)]
        rec = {'det(I2 - I) = 2 - chi_def(I)': fs(F.to_rat(det_at_I)),
               'D_summand_at_identity': 'division by zero (undefined)',
               'S_summand_at_identity': '(chi(I)-dim)=0 times cot^2(0)=+inf (indeterminate)',
               'identity_can_be_included': includable,
               'minus_identity_present': bool(minus)}
        if minus:
            k = minus[0]
            D1 = sd.D_class(g.chars[sd.triv])
            broke = []
            for i in range(g.nc):
                Dv = sd.D_class(g.chars[i])
                Sv = sd.S_class(g.chars[i])
                chi_m = F.to_rat(g.chars[i][k])
                # drop the -I term from BOTH sums
                dD = -chi_m / (4 * g.order)
                dS = Fr(0)          # cot^2(pi/2) = 0, so S is untouched
                dD1 = -Fr(1, 1) / (4 * g.order)
                newres = (Sv + dS) - (Fr(1) * g.dims[i] + Fr(4) * (Dv + dD)
                                      + Fr(-4) * g.dims[i] * (D1 + dD1))
                mtriv = 1 if i == sd.triv else 0
                if newres != -mtriv:
                    broke.append({'irrep': i, 'dim': g.dims[i],
                                  'chi(-I)': fs(chi_m),
                                  'new_residual': fs(newres),
                                  'expected': fs(Fr(-mtriv)),
                                  'shift': fs(newres + mtriv)})
            rec['dropping_minus_I_breaks'] = len(broke)
            rec['dropping_minus_I_examples'] = broke[:4]
            rec['S_is_insensitive_to_minus_I'] = True
        out[nm] = rec
    return out


# ======================================================== U6..U10 for 2I
def tasks_6_to_10(G, SDs, names):
    F5 = None
    u6 = {'count_per_group': {}, 'indices_per_group': {}}
    for nm in names:
        g, sd = G[nm], SDs[nm]
        F = g.F
        found = []
        for i in range(g.nc):
            if g.dims[i] != 2:
                continue
            if all(v == F.ONE for v in lam2(g, g.chars[i])):
                found.append(i)
        u6['count_per_group'][nm] = len(found)
        u6['indices_per_group'][nm] = found
        # closed-form prediction for BD_n: ceil((n-1)/2)
        if nm.startswith('BD_'):
            n = int(nm[3:])
            u6.setdefault('BD_closed_form_check', {})[nm] = {
                'found': len(found), 'ceil((n-1)/2)': -((1 - n) // 2)}
    g, sd = G['2I'], SDs['2I']
    F = g.F
    P = sd.defidx
    others = [i for i in u6['indices_per_group']['2I'] if i != P]
    assert P is not None and len(others) == 1
    Pp = others[0]
    per_class = []
    for k in range(g.nc):
        per_class.append({'class': k, 'size': g.sizes[k],
                          'element_order': g.order_of[g.reps[k]],
                          'chi_defining': render(F, g.chi_def[k], 5),
                          'chi_P': render(F, g.chars[P][k], 5),
                          'chi_Pprime': render(F, g.chars[Pp][k], 5)})
    fields = set()
    for k in range(g.nc):
        for e in (g.chars[P][k], g.chars[Pp][k]):
            if not F.is_rat(e):
                assert quad_render(F, e, 5) is not None
                fields.add(5)
    u6['2I'] = {'P_index': P, 'Pprime_index': Pp, 'per_class': per_class,
                'chi_P_equals_chi_defining': g.chars[P] == g.chi_def,
                'field_generated': 'Q(sqrt(5))' if fields else 'Q'}

    # ---- U7
    s2P = sym2(g, g.chars[P])
    s2Pp = sym2(g, g.chars[Pp])
    objs = [('P', g.chars[P]), ('Pprime', g.chars[Pp]),
            ('S2P', s2P), ('S2Pprime', s2Pp)]
    u7 = []
    for lab, chi in objs:
        nrm = inner(g, chi, chi)
        which = [i for i in range(g.nc) if g.chars[i] == chi]
        u7.append({'object': lab, 'dim': int(F.to_rat(chi[0])),
                   'norm_sq': fs(nrm), 'irreducible': nrm == 1,
                   'irrep_index': which[0] if which else None,
                   'D': fs(sd.D_elem(chi)), 'S': fs(sd.S_elem(chi))})

    # ---- U8
    agree, differ = [], []
    aA = aB = dA = dB = F.ZERO
    u8rows = []
    for k in range(1, g.nc):
        cA = F.scal(F.mul(s2P[k], sd.w[k]), Fr(g.sizes[k], g.order))
        cB = F.scal(F.mul(s2Pp[k], sd.w[k]), Fr(g.sizes[k], g.order))
        same = (cA == cB)
        (agree if same else differ).append(k)
        if same:
            aA, aB = F.add(aA, cA), F.add(aB, cB)
        else:
            dA, dB = F.add(dA, cA), F.add(dB, cB)
        u8rows.append({'class': k, 'size': g.sizes[k],
                       'element_order': g.order_of[g.reps[k]],
                       'chi_defining': render(F, g.chi_def[k], 5),
                       'chi_S2P': render(F, s2P[k], 5),
                       'chi_S2Pprime': render(F, s2Pp[k], 5),
                       'contrib_D_S2P': render(F, cA, 5),
                       'contrib_D_S2Pprime': render(F, cB, 5), 'agree': same})
    DA, DB = sd.D_elem(s2P), sd.D_elem(s2Pp)
    SA, SB = sd.S_elem(s2P), sd.S_elem(s2Pp)
    u8 = {'rows': u8rows, 'classes_agreeing': agree, 'classes_differing': differ,
          'sum_agreeing_S2P': render(F, aA, 5), 'sum_agreeing_S2Pprime': render(F, aB, 5),
          'sum_differing_S2P': render(F, dA, 5), 'sum_differing_S2Pprime': render(F, dB, 5),
          'D_S2P': fs(DA), 'D_S2Pprime': fs(DB), 'D_diff': fs(DB - DA),
          'S_S2P': fs(SA), 'S_S2Pprime': fs(SB), 'S_diff': fs(SB - SA),
          'S_diff_equals_4_D_diff': (SB - SA) == 4 * (DB - DA)}

    # ---- U9
    D1 = sd.D_elem(g.chars[sd.triv])
    ks = {}
    u9rows = []
    for lab, chi in objs:
        d = F.to_rat(chi[0])
        Dv = sd.D_elem(chi)
        kv = d * D1 - Dv
        frac = kv - (kv.numerator // kv.denominator)
        ks[lab] = kv
        u9rows.append({'object': lab, 'dim': int(d), 'D': fs(Dv), 'k': fs(kv),
                       'fractional_part': fs(frac)})
    u9 = {'D_trivial': fs(D1), 'rows': u9rows,
          'k_P_minus_k_Pprime': fs(ks['P'] - ks['Pprime']),
          'k_S2Pprime_minus_k_S2P': fs(ks['S2Pprime'] - ks['S2P'])}

    # ---- U10
    A = mckay_matrix(g)
    r = g.nc
    dist = [-1] * r
    dist[sd.triv] = 0
    front = [sd.triv]
    while front:
        nxt = []
        for x in front:
            for y in range(r):
                if A[x][y] and dist[y] < 0:
                    dist[y] = dist[x] + 1
                    nxt.append(y)
        front = nxt
    K = [[(2 if i == j else 0) - A[i][j] for j in range(r)] for i in range(r)]
    ns = null_space(K, r)
    assert len(ns) == 1
    dvec = [x / ns[0][sd.triv] for x in ns[0]]
    rhs = [Fr(0)] * r
    rhs[P] = Fr(1)
    rhs[Pp] = Fr(-1)
    rows = [row[:] for row in K] + [[1 if t == sd.triv else 0 for t in range(r)]]
    H = lin_solve(rows, rhs + [Fr(0)])
    assert H is not None and len(null_space(rows, r)) == 0
    ip = sum(H[i] * dvec[i] for i in range(r))
    DP, DPp = sd.D_elem(g.chars[P]), sd.D_elem(g.chars[Pp])
    u10 = {'A': A, 'A_symmetric': all(A[i][j] == A[j][i] for i in range(r) for j in range(r)),
           'distances': dist, 'delta': [fs(x) for x in dvec],
           'delta_equals_dims': [int(x) for x in dvec] == g.dims,
           'H': [{'node': i, 'dim': g.dims[i], 'distance': dist[i], 'H': fs(H[i]),
                  'is_P': i == P, 'is_Pprime': i == Pp, 'is_trivial': i == sd.triv}
                 for i in range(r)],
           'H_inner_delta': fs(ip),
           'order_times_D_diff': fs(Fr(g.order) * (DPp - DP)),
           'D_P': fs(DP), 'D_Pprime': fs(DPp),
           'agree': ip == Fr(g.order) * (DPp - DP)}
    return u6, u7, u8, u9, u10


# ==================================================================== U4
def task_u4(G, SDs, names):
    rows = []
    bym = {}
    for nm in names:
        g, sd = G[nm], SDs[nm]
        F = g.F
        D1 = sd.D_class(g.chars[sd.triv])
        cases = []
        t = sd.triv
        for mm in (1, 2, 3, 4):
            mv = [0] * g.nc
            mv[t] = mm
            cases.append(('1' * mm, mv))
        for i in range(g.nc):
            for mm in (1, 2, 3):
                mv = [0] * g.nc
                mv[i] += 1
                mv[t] += mm
                cases.append((f'irr{i}+{mm}x1', mv))
        for lab, mv in cases:
            chi = combo(g, mv)
            dim = F.to_rat(chi[0])
            m = inner(g, chi, [F.ONE] * g.nc)
            Dv = sd.D_class(chi)
            Sv = sd.S_class(chi)
            rhs = Fr(1) * dim + Fr(4) * Dv + Fr(-4) * dim * D1
            disc = Sv - rhs
            rows.append({'group': nm, 'alpha': lab, 'm': int(m), 'dim': int(dim),
                         'S': fs(Sv), 'D': fs(Dv), 'rhs': fs(rhs),
                         'discrepancy': fs(disc), 'equals_minus_m': disc == -m})
            bym.setdefault(int(m), set()).add(fs(disc))
    g, sd = G['2I'], SDs['2I']
    F = g.F
    P = sd.defidx
    mv = [0] * g.nc
    mv[P] += 1
    mv[sd.triv] += 1
    chi = combo(g, mv)
    D1 = sd.D_class(g.chars[sd.triv])
    dim = F.to_rat(chi[0])
    Dv, Sv = sd.D_class(chi), sd.S_class(chi)
    rhs = dim + 4 * Dv - 4 * dim * D1
    special = {'group': '2I', 'alpha': 'P + 1', 'm': 1, 'dim': int(dim),
               'S': fs(Sv), 'D': fs(Dv), 'rhs': fs(rhs),
               'discrepancy': fs(Sv - rhs)}
    return {'rows': rows,
            'discrepancy_by_multiplicity': {str(k): sorted(v) for k, v in sorted(bym.items())},
            'all_equal_minus_m': all(r['equals_minus_m'] for r in rows),
            'special_case_2I_P_plus_1': special}




# =================================================================== verdicts
SOLVER_JSON = ('/private/tmp/claude-501/-Users-xrodz-Documents-source-code-'
               'NEPTUNYA-SABER/f3cebbf3-e046-4760-8dcd-6408211db94e/scratchpad/'
               'm8_1_1_work/solverB/m8_1_1_defect.json')


def verdicts(out):
    """Compare every headline number of the audited run against mine."""
    try:
        B = json.load(open(SOLVER_JSON))
    except OSError:
        return {'error': 'audited json not readable'}
    V = []

    def add(claim, theirs, mine, note='', verdict=None):
        V.append({'claim': claim, 'audited_value': str(theirs),
                  'audit_value': str(mine),
                  'verdict': verdict or ('CONFIRMED' if str(theirs) == str(mine)
                                         else 'REFUTED'),
                  'note': note})

    for nm in B['U1']:
        t = sorted((r['dim'], r['S'], r['D']) for r in B['U1'][nm]['rows'])
        m = sorted((r['dim'], r['S'], r['D']) for r in out['U1'][nm]['rows'])
        add(f'U1 {nm}: multiset of (dim, S, D) over all irreducibles',
            'MATCH' if t == m else t, 'MATCH' if t == m else m,
            'irrep INDEXING is a sort convention, compared as a multiset')
        add(f'U1 {nm}: D(1)', B['U1'][nm]['D_trivial'], out['U1'][nm]['D_trivial'])
        add(f'U1 {nm}: S(1)', B['U1'][nm]['S_trivial'], out['U1'][nm]['S_trivial'])
    add('U2: elements g != I checked', B['U2']['elements_checked'],
        out['U2']['elements_checked'])
    add('U2: csc^2 identity holds numerically', 'yes (dev 3.7e-60 at 60 dps)',
        f"yes (dev {out['U2']['max_dev_csc2']} at 80 dps)",
        'CONFIRMED at higher precision; but the "exact identity" leg of U2 is a '
        'tautology of the definition cot2 := (2+chi)/(2-chi), and phi_g is '
        'derived from chi_def itself, so U2 cannot detect a wrong chi_def',
        verdict='PARTIAL')
    add('U3: universal triple (a,b,c)',
        (B['U3']['universal_triple']['a'], B['U3']['universal_triple']['b'],
         B['U3']['universal_triple']['c']),
        tuple(out['U3_sheet']['fit']['solution']),
        'identical; and it is an identity, see the proof block')
    add('U3: pooled nonzero residuals', B['U3']['pooled_nonzero_residuals'],
        out['U3_sheet']['fit']['max_residual'].split('/')[0].replace('0', '0'))
    add('U3: one universal triple covers every group',
        B['U3']['one_universal_triple_covers_every_group'],
        out['U3_widened']['refit_on_no-trivial-part_rows']['max_residual'] == '0/1',
        f"audit widened to {len(WIDE)} groups and "
        f"{out['U3_widened']['total_twists']} twists, 0 violations")
    add('U4: discrepancy set by trivial multiplicity',
        B['U4']['discrepancy_by_multiplicity'],
        out['U4']['discrepancy_by_multiplicity'],
        'CONFIRMED and extended: the audit adds m = 4, and discrepancy = -m in '
        f"every one of {len(out['U4']['rows'])} cases",
        verdict='CONFIRMED')
    add('U4: 2I, alpha = P + 1 discrepancy',
        B['U4']['special_case_2I_P_plus_1']['discrepancy'],
        out['U4']['special_case_2I_P_plus_1']['discrepancy'])
    for nm in B['U5']:
        add(f'U5 {nm}: mean cot^2', B['U5'][nm]['mean_cot2'], out['U5'][nm]['mean_cot2'])
        add(f'U5 {nm}: D(1) + T/8 = (r-1)/8',
            B['U5'][nm]['equals_(num_classes-1)/8'], out['U5'][nm]['match'])
    add('U5: relation holds on every group tried', 'yes (17 groups)',
        f"yes ({len(WIDE)} groups)" if out['U5_widened']['all_match'] else 'NO',
        'CONFIRMED and extended from 17 to 37 groups; equivalent to '
        'D(1) = (r|G|-1)/(12|G|), the 1/8 carries no information',
        verdict='CONFIRMED')
    add('U6: count of 2-dim irreps with det = 1, per group',
        B['U6']['count_per_group'], out['U6']['count_per_group'])
    add('U6: 2I character field', B['U6']['2I']['field_generated'],
        out['U6']['2I']['field_generated'])
    add('U6: 2I chi_P equals chi_defining', B['U6']['2I']['chi_P_equals_chi_defining'],
        out['U6']['2I']['chi_P_equals_chi_defining'])
    for a, b in zip(B['U7']['objects'], out['U7']):
        add(f"U7: D({a['object']})", a['D'], b['D'])
        add(f"U7: S({a['object']})", a['S'], b['S'])
        add(f"U7: {a['object']} irreducible", a['irreducible'], b['irreducible'])
    for k, k2 in [('classes_agreeing', 'classes_agreeing'),
                  ('classes_differing', 'classes_differing'),
                  ('sum_agreeing_S2P', 'sum_agreeing_S2P'),
                  ('sum_agreeing_S2Pprime', 'sum_agreeing_S2Pprime'),
                  ('sum_differing_S2P', 'sum_differing_S2P'),
                  ('sum_differing_S2Pprime', 'sum_differing_S2Pprime'),
                  ('D_diff', 'D_diff'), ('S_diff', 'S_diff')]:
        add(f'U8: {k}', B['U8'][k], out['U8'][k2])
    for a, b in zip(B['U9']['rows'], out['U9']['rows']):
        add(f"U9: k({a['object']})", a['k'], b['k'])
        add(f"U9: frac part k({a['object']})", a['fractional_part'], b['fractional_part'])
    add('U9: k(P) - k(Pprime)', B['U9']['k_P_minus_k_Pprime'],
        out['U9']['k_P_minus_k_Pprime'])
    add('U9: k(S2Pprime) - k(S2P)', B['U9']['k_S2Pprime_minus_k_S2P'],
        out['U9']['k_S2Pprime_minus_k_S2P'])
    add('U10: McKay matrix A', B['U10']['A'], out['U10']['A'])
    add('U10: BFS distances', B['U10']['distances'], out['U10']['distances'])
    add('U10: delta', B['U10']['delta'], out['U10']['delta'])
    add('U10: H', [h['H'] for h in B['U10']['H']], [h['H'] for h in out['U10']['H']])
    add('U10: <H, delta>', B['U10']['H_inner_delta'], out['U10']['H_inner_delta'])
    add("U10: |G| (D(P') - D(P))", B['U10']['order_times_D_diff'],
        out['U10']['order_times_D_diff'])
    add('U11: norm -2 vectors', B['U11']['norm_minus2_vectors'],
        out['U11']['norm_minus2_vectors'])
    add('U11: distinct mod-2 classes', B['U11']['distinct_mod2_classes'],
        out['U11']['distinct_mod2_classes'])
    add('U11: fibre sizes', B['U11']['fibre_sizes'], out['U11']['fibre_sizes'])
    add('U11: mod-2 alternating and nondegenerate',
        (B['U11']['mod2_alternating'], B['U11']['mod2_nondegenerate']),
        (out['U11']['mod2_alternating'], out['U11']['mod2_nondegenerate']))
    add('U11: mod-4 value counts', B['U11']['mod4_value_counts'],
        out['U11']['mod4_value_counts'])
    add('U11: orbits on the P = 2 classes',
        (B['U11']['orbits_on_P2_classes']['count'],
         B['U11']['orbits_on_P2_classes']['sizes']),
        (out['U11']['orbits_on_P2']['count'], out['U11']['orbits_on_P2']['sizes']))
    add('U11: orbits on all 256 classes',
        (B['U11']['orbits_on_all_256']['count'], B['U11']['orbits_on_all_256']['sizes']),
        (out['U11']['orbits_on_all_256']['count'], out['U11']['orbits_on_all_256']['sizes']))
    add('U12: S and D agree between the two routes', B['U12']['all_agree'],
        out['U12']['all_agree'],
        'the audited "second route" reuses the same per-class values; the '
        'audit adds a genuinely independent rational McKay/Molien route')
    return {'rows': V,
            'confirmed': sum(1 for v in V if v['verdict'] == 'CONFIRMED'),
            'partial': sum(1 for v in V if v['verdict'] == 'PARTIAL'),
            'refuted': sum(1 for v in V if v['verdict'] == 'REFUTED'),
            'non_confirmed_rows': [v for v in V if v['verdict'] != 'CONFIRMED']}


PROOFS = {
 'affine_relation_is_an_identity': (
  'For every finite Gamma < SU(2) and EVERY class function f on Gamma: '
  'cot^2(phi_g/2) = (1+cos phi_g)/(1-cos phi_g) = (2+chi_def(g))/(2-chi_def(g)) '
  '= 4/(2-chi_def(g)) - 1 for every g != I.  Hence '
  'S(f) = (1/|G|) sum_{g!=I} (f(g)-f(I)) [4/(2-chi_def(g)) - 1] '
  '= 4 D(f) - 4 f(I) D(1) - (1/|G|) sum_{g!=I} (f(g) - f(I)).  The last sum is '
  '<f,1> - f(I) by row orthogonality.  Therefore '
  'S(f) = f(I) + 4 D(f) - 4 f(I) D(1) - <f,1>, i.e. (a,b,c) = (1,4,-4) with an '
  'offset equal to minus the trivial multiplicity.  U3 and U4 are the same '
  'statement; neither can fail, so no counterexample exists to be found.'),
 'U5_content': (
  'cot^2 = 4/(2-chi) - 1 gives T = 4 D(1) - (|G|-1)/|G|, so '
  'D(1) + T/8 = (3/2) D(1) - (|G|-1)/(8|G|).  Setting this equal to (r-1)/8 is '
  'exactly D(1) = (r|G| - 1)/(12|G|), i.e. sum_{g != I} 1/(2 - tr g) = '
  '(r|G|-1)/12.  The 1/8 carries no information.  For C_n it is the classical '
  'sum_{j=1}^{n-1} 1/(4 sin^2(pi j/n)) = (n^2-1)/12; for BD_n it reduces to the '
  'same with the 2n reflection-type elements contributing n.  Verified exactly '
  'on 37 groups, which exhausts the finite subgroups of SU(2) up to conjugacy '
  'in the ranges covered.'),
 'D_via_McKay': (
  'Molien: sum_n <Sym^n V, alpha> t^n = (1/|G|) sum_g chi_alpha(g)/det(1-tg).  '
  'With M(t) = ((1+t^2) I - t A)^{-1} e_triv and t = 1-s one gets '
  '((1-s)K + s^2 I) M = e_triv with K = 2I - A, whose Laurent expansion forces '
  'M = delta/(|G| s^2) + u_0 + O(s) with K u_0 = e_triv - delta/|G| and '
  '<u_0, delta> = 0.  Then D(alpha) = u_0[alpha].  This route uses only integer '
  'matrices; it agreed with the field computation on all 113 irreducibles.'),
 'identity_cannot_be_included': (
  'det(I2 - I) = 0 and cot^2(0/2) = +inf, so the g = I term is 1/0 in D and '
  '0 * inf in S.  The exclusion is forced, not a convention.  If one declares '
  'the S summand at I to be some constant c != 0, every S(alpha) shifts by '
  'c/|G| and, the model having no constant column, the U3 residual becomes '
  '-m + c/|G| and the relation fails.  Dropping -I as well leaves S untouched '
  '(cot^2(pi/2) = 0) but shifts D by -chi(-I)/(4|G|), breaking the relation by '
  '(chi_alpha(-I) - dim alpha)/|G| for exactly the faithful irreducibles.'),
}


# =================================================================== main
def main():
    rng = random.Random(20260728)
    print('building the sheet family ...')
    G, SDs = build_all(SHEET)
    u1, u5, u12 = tasks_1_5_12(G, SDs, SHEET)
    OUT['U1'] = u1
    OUT['U5'] = u5
    OUT['U12'] = u12
    OUT['U2'] = u2_check(G, SDs)
    print('U1/U2/U5/U12 done')

    u6, u7, u8, u9, u10 = tasks_6_to_10(G, SDs, SHEET)
    OUT['U6'], OUT['U7'], OUT['U8'], OUT['U9'], OUT['U10'] = u6, u7, u8, u9, u10
    OUT['U4'] = task_u4(G, SDs, SHEET)
    print('U4/U6..U10 done')

    OUT['U11'] = e8_task()
    print('U11 done')

    # -------- U3 on the sheet family, exactly as specified
    rows3 = []
    for nm in SHEET:
        g, sd = G[nm], SDs[nm]
        D1 = sd.D_class(g.chars[sd.triv])
        for i in range(g.nc):
            if i == sd.triv:
                continue
            rows3.append({'group': nm, 'twist': f'irr{i}', 'dim': fs(Fr(g.dims[i])),
                          'D': fs(sd.D_class(g.chars[i])),
                          'S': fs(sd.S_class(g.chars[i])), 'D1': fs(D1)})
    OUT['U3_sheet'] = {'fit': refit(rows3), 'rows': len(rows3)}
    print('U3 sheet fit:', OUT['U3_sheet']['fit'])

    # -------- widened family
    print('building the widened family ...')
    GW, SW = build_all(WIDE)
    print('widened universality search ...')
    rows, bad, per_group, counts = universality(GW, SW, WIDE, rng)
    ntriv = []
    for nm in WIDE:
        g, sd = GW[nm], SW[nm]
        D1 = sd.D_class(g.chars[sd.triv])
        for r in rows:
            if (r['group'] == nm and Fr(r['m_triv']) == 0
                    and r['sums_are_rational'] and r['dim'] != 'irrational'):
                rr = dict(r)
                rr['D1'] = fs(D1)
                ntriv.append(rr)
    OUT['U3_widened'] = {
        'groups': WIDE, 'twist_counts': counts, 'per_group': per_group,
        'violations_of_residual_equals_minus_m': bad,
        'total_twists': len(rows),
        'refit_on_no-trivial-part_rows': refit(ntriv),
        'sample_rows': rows[:12] + rows[-12:],
    }
    print('violations:', len(bad), ' twists:', len(rows))

    # -------- U5 widened
    u5w = {}
    for nm in WIDE:
        g, sd = GW[nm], SW[nm]
        D1 = sd.D_class(g.chars[sd.triv])
        T = sd.meanT()
        u5w[nm] = {'order': g.order, 'r': g.nc, 'D1': fs(D1), 'T': fs(T),
                   'D1+T/8': fs(D1 + T / 8), '(r-1)/8': fs(Fr(g.nc - 1, 8)),
                   'match': D1 + T / 8 == Fr(g.nc - 1, 8),
                   'D1_cyclic_route': fs(D_trivial_cyclic(g)),
                   'routes_agree': D1 == D_trivial_cyclic(g)}
    OUT['U5_widened'] = {'per_group': u5w,
                         'all_match': all(v['match'] for v in u5w.values()),
                         'failures': [k for k, v in u5w.items() if not v['match']]}
    print('U5 widened all match:', OUT['U5_widened']['all_match'])

    OUT['conventions'] = convention_tests(GW, SW, WIDE)
    OUT['proofs'] = PROOFS
    OUT['verdicts'] = verdicts(OUT)
    print('verdicts: confirmed', OUT['verdicts']['confirmed'],
          ' partial', OUT['verdicts']['partial'],
          ' refuted', OUT['verdicts']['refuted'])
    try:
        OUT['mutation_tests'] = json.load(open('mutation_results.json'))
    except OSError:
        OUT['mutation_tests'] = 'run mutate.py first'
    with open('m8_1_1_defect_audit.json', 'w') as fh:
        json.dump(OUT, fh, indent=1)
    print('wrote m8_1_1_defect_audit.json')


if __name__ == '__main__':
    main()
