#!/usr/bin/env python3
"""SPEC SHEET B -- two class-function sums on finite subgroups of SU(2),
and an integral lattice count.  Tasks U1 .. U12.  Self-contained.

Everything is exact.  Group elements are 2x2 matrices over Q(zeta_N) held in a
canonical power-basis representation with integer numerators over a common
denominator, so equality is structural and elements are hashable.  Character
tables are computed by the Burnside-Dixon modular method (class multiplication
matrices, simultaneous eigenvectors over F_p with p = 1 mod exponent, then the
eigenvalue-multiplicity lift back to Q(zeta_N)); nothing is imported.  S(alpha)
and D(alpha) land as exact Fractions; character values that live in a real
quadratic field are carried in the cyclotomic field and rendered exactly.

Outputs:  m8_1_1_defect.json  plus a full transcript on stdout.
"""

from fractions import Fraction
from math import gcd, isqrt
import itertools
import json
import sys

import mpmath as mp
import sympy

mp.mp.dps = 60


# =====================================================================
# PART 1 -- exact arithmetic in the cyclotomic field Q(zeta_N)
# =====================================================================
_X = sympy.Symbol('x')


class Cyc:
    def __init__(self, N):
        self.N = N
        p = sympy.Poly(sympy.cyclotomic_poly(N, _X), _X)
        cs = [int(c) for c in p.all_coeffs()]          # highest degree first
        self.deg = len(cs) - 1
        self.mod = list(reversed(cs))                  # mod[i] = coeff of x^i
        assert self.mod[self.deg] == 1
        self.ZERO = ((0,) * self.deg, 1)
        self.ONE = self.rat(1)
        self._inv_cache = {}

    # ---------- construction ----------
    def rat(self, q):
        q = Fraction(q)
        nums = [0] * self.deg
        nums[0] = q.numerator
        return self._norm(nums, q.denominator)

    def zeta(self, k):
        """zeta_N^k as a field element."""
        k %= self.N
        L = [0] * max(self.deg, k + 1)
        L[k] = 1
        return self._norm(self._reduce(L), 1)

    # ---------- internals ----------
    def _norm(self, nums, den):
        if den < 0:
            nums = [-n for n in nums]
            den = -den
        g = den
        for n in nums:
            g = gcd(g, n if n >= 0 else -n)
        if g > 1:
            nums = [n // g for n in nums]
            den //= g
        return (tuple(nums), den)

    def _reduce(self, L):
        d = len(L) - 1
        while d >= self.deg:
            c = L[d]
            if c:
                L[d] = 0
                base = d - self.deg
                for i in range(self.deg):
                    L[base + i] -= c * self.mod[i]
            d -= 1
        out = L[:self.deg]
        while len(out) < self.deg:
            out.append(0)
        return out

    # ---------- ring operations ----------
    def add(self, a, b):
        an, ad = a
        bn, bd = b
        g = gcd(ad, bd)
        m1 = bd // g
        m2 = ad // g
        nums = [an[i] * m1 + bn[i] * m2 for i in range(self.deg)]
        return self._norm(nums, ad * m1)

    def neg(self, a):
        return (tuple(-x for x in a[0]), a[1])

    def sub(self, a, b):
        return self.add(a, self.neg(b))

    def mul(self, a, b):
        an, ad = a
        bn, bd = b
        d = self.deg
        conv = [0] * (2 * d - 1)
        for i in range(d):
            ai = an[i]
            if ai:
                for j in range(d):
                    bj = bn[j]
                    if bj:
                        conv[i + j] += ai * bj
        return self._norm(self._reduce(conv), ad * bd)

    def scal(self, a, q):
        q = Fraction(q)
        return self._norm([x * q.numerator for x in a[0]], a[1] * q.denominator)

    def inv(self, a):
        if a in self._inv_cache:
            return self._inv_cache[a]
        if a == self.ZERO:
            raise ZeroDivisionError('inverse of 0 in Q(zeta_N)')
        # extended Euclid in Q[x] between a(x) and Phi_N(x)
        A = [Fraction(n, a[1]) for n in a[0]]
        while A and A[-1] == 0:
            A.pop()
        M = [Fraction(c) for c in self.mod]
        r0, r1 = M[:], A[:]
        s0, s1 = [Fraction(0)], [Fraction(1)]
        while any(c != 0 for c in r1):
            q, r = _pdivmod(r0, r1)
            r0, r1 = r1, r
            s0, s1 = s1, _psub(s0, _pmul(q, s1))
        # r0 = gcd (a nonzero constant since Phi_N is irreducible)
        c = r0[-1]
        res = [x / c for x in s0]
        res = res[:self.deg] + [Fraction(0)] * max(0, self.deg - len(res))
        den = 1
        for x in res:
            den = den * x.denominator // gcd(den, x.denominator)
        out = self._norm([int(x * den) for x in res], den)
        chk = self.mul(a, out)
        assert chk == self.ONE, 'inverse failed'
        self._inv_cache[a] = out
        return out

    def div(self, a, b):
        return self.mul(a, self.inv(b))

    def galois(self, a, t):
        """Apply zeta -> zeta^t (t coprime to N)."""
        t %= self.N
        size = max(self.deg, (self.deg - 1) * t + 1)
        L = [0] * size
        for j in range(self.deg):
            if a[0][j]:
                L[j * t] += a[0][j]
        return self._norm(self._reduce(L), a[1])

    def conj(self, a):
        return self.galois(a, self.N - 1)

    # ---------- predicates / conversion ----------
    def is_rational(self, a):
        return all(x == 0 for x in a[0][1:])

    def to_fraction(self, a):
        assert self.is_rational(a), 'element is not rational'
        return Fraction(a[0][0], a[1])

    def numeric(self, a):
        s = mp.mpc(0)
        for j in range(self.deg):
            if a[0][j]:
                s += mp.mpf(a[0][j]) * mp.expjpi(mp.mpf(2 * j) / self.N)
        return s / a[1]


# ---------- helper polynomial ops over Q ----------
def _ptrim(p):
    p = list(p)
    while p and p[-1] == 0:
        p.pop()
    return p


def _psub(a, b):
    n = max(len(a), len(b))
    out = [Fraction(0)] * n
    for i, x in enumerate(a):
        out[i] += x
    for i, x in enumerate(b):
        out[i] -= x
    return _ptrim(out)


def _pmul(a, b):
    if not a or not b:
        return []
    out = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    out[i + j] += x * y
    return _ptrim(out)


def _pdivmod(a, b):
    a = _ptrim(a)
    b = _ptrim(b)
    q = [Fraction(0)] * max(0, len(a) - len(b) + 1)
    r = a[:]
    while r and len(r) >= len(b):
        c = r[-1] / b[-1]
        d = len(r) - len(b)
        q[d] = c
        for i in range(len(b)):
            r[d + i] -= c * b[i]
        r = _ptrim(r)
    return _ptrim(q), r

# =====================================================================
# PART 2 -- exact rendering of field elements and rational linear algebra
# =====================================================================


def fs(x):
    x = Fraction(x)
    return f"{x.numerator}/{x.denominator}"


def rref_frac(M, ncols):
    rows = [[Fraction(v) for v in r] for r in M]
    piv = []
    r = 0
    for c in range(ncols):
        pr = None
        for rr in range(r, len(rows)):
            if rows[rr][c] != 0:
                pr = rr
                break
        if pr is None:
            continue
        rows[r], rows[pr] = rows[pr], rows[r]
        pv = rows[r][c]
        rows[r] = [v / pv for v in rows[r]]
        for rr in range(len(rows)):
            if rr != r and rows[rr][c] != 0:
                f = rows[rr][c]
                rows[rr] = [rows[rr][k] - f * rows[r][k] for k in range(len(rows[rr]))]
        piv.append(c)
        r += 1
        if r == len(rows):
            break
    return rows, piv


def solve_exact(A, b):
    """Solve A x = b exactly.  Returns (particular, nullbasis) or None."""
    n = len(A)
    m = len(A[0]) if n else 0
    aug = [[Fraction(A[i][j]) for j in range(m)] + [Fraction(b[i])] for i in range(n)]
    rows, piv = rref_frac(aug, m)
    for row in rows:
        if all(row[j] == 0 for j in range(m)) and row[m] != 0:
            return None
    x = [Fraction(0)] * m
    for i, c in enumerate(piv):
        x[c] = rows[i][m]
    free = [c for c in range(m) if c not in piv]
    null = []
    for fc in free:
        v = [Fraction(0)] * m
        v[fc] = Fraction(1)
        for i, c in enumerate(piv):
            v[c] = -rows[i][fc]
        null.append(v)
    return x, null


def minpoly(F, e):
    """Monic minimal polynomial coefficients [c0, c1, ..., 1] over Q."""
    pows = [F.ONE]
    cur = F.ONE
    for _ in range(F.deg):
        cur = F.mul(cur, e)
        pows.append(cur)
    for k in range(1, F.deg + 1):
        A = [[Fraction(pows[i][0][row], pows[i][1]) for i in range(k)]
             for row in range(F.deg)]
        b = [-Fraction(pows[k][0][row], pows[k][1]) for row in range(F.deg)]
        sol = solve_exact(A, b)
        if sol is not None:
            return list(sol[0]) + [Fraction(1)]
    raise RuntimeError('minimal polynomial not found')


def squarefree(n):
    """n = f^2 * d with d squarefree; return (f, d) for n > 0."""
    f = 1
    d = n
    i = 2
    while i * i <= d:
        while d % (i * i) == 0:
            d //= i * i
            f *= i
        i += 1
    return f, d


def render(F, e):
    """Return (kind, string) exact representation."""
    if F.is_rational(e):
        return 'rational', fs(F.to_fraction(e))
    for j in range(F.N):
        if e == F.zeta(j):
            return 'root_of_unity', f"exp(2*I*pi*{j}/{F.N})"
        if e == F.neg(F.zeta(j)):
            return 'root_of_unity', f"-exp(2*I*pi*{j}/{F.N})"
    mp_ = minpoly(F, e)
    if len(mp_) == 3:
        c0, c1 = mp_[0], mp_[1]
        disc = c1 * c1 - 4 * c0
        a = -c1 / 2
        val = F.numeric(e)
        if disc > 0:
            f, d = squarefree(disc.numerator * disc.denominator)
            bpos = Fraction(f, disc.denominator) / 2
            cand = [mp.mpf(a.numerator) / a.denominator
                    + s * mp.mpf(bpos.numerator) / bpos.denominator * mp.sqrt(d)
                    for s in (1, -1)]
            sign = 1 if abs(val.real - cand[0]) < abs(val.real - cand[1]) else -1
            b = sign * bpos
            return 'quadratic', f"{fs(a)} + ({fs(b)})*sqrt({d})"
        else:
            f, d = squarefree((-disc).numerator * (-disc).denominator)
            bpos = Fraction(f, (-disc).denominator) / 2
            sign = 1 if val.imag > 0 else -1
            b = sign * bpos
            return 'quadratic_imag', f"{fs(a)} + ({fs(b)})*I*sqrt({d})"
    poly = " + ".join(f"({fs(c)})*x^{i}" for i, c in enumerate(mp_) if c != 0)
    coords = "[" + ", ".join(fs(Fraction(n, e[1])) for n in e[0]) + "]"
    return 'algebraic', f"root of {poly}; coords in basis zeta_{F.N}^j: {coords}"


def numstr(F, e, dps=25):
    v = F.numeric(e)
    if abs(v.imag) < mp.mpf(10) ** (-30):
        return mp.nstr(v.real, 20)
    return mp.nstr(v.real, 20) + (" + " if v.imag >= 0 else " - ") + mp.nstr(abs(v.imag), 20) + "*I"

# =====================================================================
# PART 3 -- groups, conjugacy classes, Burnside-Dixon character tables
# =====================================================================


# ------------------------------------------------------------------ matrices
def mmul(F, X, Y):
    a, b, c, d = X
    e, f, g, h = Y
    return (F.add(F.mul(a, e), F.mul(b, g)),
            F.add(F.mul(a, f), F.mul(b, h)),
            F.add(F.mul(c, e), F.mul(d, g)),
            F.add(F.mul(c, f), F.mul(d, h)))


def minv(F, X):
    """Inverse of a determinant-one 2x2 matrix."""
    a, b, c, d = X
    return (d, F.neg(b), F.neg(c), a)


def mdet(F, X):
    a, b, c, d = X
    return F.sub(F.mul(a, d), F.mul(b, c))


def mtr(F, X):
    return F.add(X[0], X[3])


def quat(F, a, b, c, d, I):
    """a + b i + c j + d k  ->  [[a+bI, c+dI], [-c+dI, a-bI]]."""
    return (F.add(a, F.mul(b, I)), F.add(c, F.mul(d, I)),
            F.add(F.neg(c), F.mul(d, I)), F.sub(a, F.mul(b, I)))


def lcm(a, b):
    return a * b // gcd(a, b)


# ------------------------------------------------------------- constructions
def make_group(name):
    """Return (F, generators, expected_order, description)."""
    if name.startswith('C_'):
        n = int(name[2:])
        N = lcm(n, 4)
        F = Cyc(N)
        z = F.zeta(N // n)
        zi = F.zeta(-(N // n))
        g = (z, F.ZERO, F.ZERO, zi)
        return F, [g], n
    if name.startswith('BD_'):
        n = int(name[3:])
        N = lcm(2 * n, 4)
        F = Cyc(N)
        w = F.zeta(N // (2 * n))
        wi = F.zeta(-(N // (2 * n)))
        a = (w, F.ZERO, F.ZERO, wi)
        one = F.ONE
        b = (F.ZERO, one, F.neg(one), F.ZERO)
        return F, [a, b], 4 * n
    if name == '2T':
        F = Cyc(12)
        I = F.zeta(3)                      # zeta_12^3 = i
        h = F.rat(Fraction(1, 2))
        g1 = quat(F, h, h, h, h, I)        # (1+i+j+k)/2
        g2 = quat(F, F.ZERO, F.ONE, F.ZERO, F.ZERO, I)   # i
        return F, [g1, g2], 24
    if name == '2O':
        F = Cyc(24)
        I = F.zeta(6)                      # zeta_24^6 = i
        h = F.rat(Fraction(1, 2))
        g1 = quat(F, h, h, h, h, I)
        g2 = quat(F, F.ZERO, F.ONE, F.ZERO, F.ZERO, I)
        z8 = F.zeta(3)                     # zeta_24^3 = zeta_8 = (1+i)/sqrt2
        g3 = (z8, F.ZERO, F.ZERO, F.zeta(-3))
        return F, [g1, g2, g3], 48
    if name == '2I':
        F = Cyc(60)
        I = F.zeta(15)                     # zeta_60^15 = i
        h = F.rat(Fraction(1, 2))
        g1 = quat(F, h, h, h, h, I)
        z5 = F.zeta(12)                    # zeta_60^12 = zeta_5
        # sqrt5 = z5 + z5^4 - z5^2 - z5^3
        s5 = F.sub(F.add(F.zeta(12), F.zeta(48)), F.add(F.zeta(24), F.zeta(36)))
        assert F.mul(s5, s5) == F.rat(5)
        phi = F.scal(F.add(F.ONE, s5), Fraction(1, 2))          # (1+sqrt5)/2
        phinv = F.scal(F.sub(s5, F.ONE), Fraction(1, 2))        # (sqrt5-1)/2
        assert F.mul(phi, phinv) == F.ONE
        g2 = quat(F, F.scal(phi, Fraction(1, 2)),
                  F.scal(phinv, Fraction(1, 2)),
                  F.rat(Fraction(1, 2)), F.ZERO, I)
        return F, [g1, g2], 120
    raise ValueError(name)


def close_group(F, gens):
    ident = (F.ONE, F.ZERO, F.ZERO, F.ONE)
    seen = {ident: 0}
    order = [ident]
    frontier = [ident]
    while frontier:
        nxt = []
        for x in frontier:
            for g in gens:
                y = mmul(F, x, g)
                if y not in seen:
                    seen[y] = len(order)
                    order.append(y)
                    nxt.append(y)
        frontier = nxt
    return order, seen


def conj_classes(F, elems, index, gens):
    n = len(elems)
    cls = [-1] * n
    classes = []
    for i in range(n):
        if cls[i] >= 0:
            continue
        k = len(classes)
        members = [i]
        cls[i] = k
        stack = [i]
        while stack:
            a = stack.pop()
            xa = elems[a]
            for g in gens:
                y = mmul(F, mmul(F, g, xa), minv(F, g))
                j = index[y]
                if cls[j] < 0:
                    cls[j] = k
                    members.append(j)
                    stack.append(j)
        classes.append(sorted(members))
    return classes, cls


def elem_order(F, elems, index, i):
    ident = 0
    x = elems[i]
    o = 1
    cur = x
    while index[cur] != ident:
        cur = mmul(F, cur, x)
        o += 1
    return o


# ------------------------------------------------------- modular linear algebra
def nullspace_mod(M, p, nc):
    rows = [r[:] for r in M]
    nr = len(rows)
    piv = []
    r = 0
    for c in range(nc):
        pr = None
        for rr in range(r, nr):
            if rows[rr][c] % p:
                pr = rr
                break
        if pr is None:
            continue
        rows[r], rows[pr] = rows[pr], rows[r]
        iv = pow(rows[r][c], p - 2, p)
        rows[r] = [(x * iv) % p for x in rows[r]]
        for rr in range(nr):
            if rr != r and rows[rr][c] % p:
                f = rows[rr][c]
                rows[rr] = [(rows[rr][cc] - f * rows[r][cc]) % p for cc in range(nc)]
        piv.append(c)
        r += 1
        if r == nr:
            break
    free = [c for c in range(nc) if c not in piv]
    basis = []
    for fc in free:
        v = [0] * nc
        v[fc] = 1
        for i, pc in enumerate(piv):
            v[pc] = (-rows[i][fc]) % p
        basis.append(v)
    return basis


def is_prime(n):
    if n < 2:
        return False
    for q in range(2, isqrt(n) + 1):
        if n % q == 0:
            return False
    return True


def primitive_root(p):
    fac = []
    m = p - 1
    d = 2
    while d * d <= m:
        if m % d == 0:
            fac.append(d)
            while m % d == 0:
                m //= d
        d += 1
    if m > 1:
        fac.append(m)
    for g in range(2, p):
        if all(pow(g, (p - 1) // f, p) != 1 for f in fac):
            return g
    raise RuntimeError('no primitive root')


# ------------------------------------------------------------- Dixon's method
class GroupData:
    pass


def build(name):
    F, gens, expected = make_group(name)
    elems, index = close_group(F, gens)
    assert len(elems) == expected, f'{name}: closed to {len(elems)}, expected {expected}'
    G = len(elems)
    ident = (F.ONE, F.ZERO, F.ZERO, F.ONE)
    assert index[ident] == 0

    # every element is in SU(2): det = 1 and g^{-1} = conjugate transpose
    for x in elems:
        assert mdet(F, x) == F.ONE
        a, b, c, d = x
        assert F.conj(a) == d and F.conj(b) == F.neg(c)

    classes, cls = conj_classes(F, elems, index, gens)
    r = len(classes)
    reps = [c[0] for c in classes]
    sizes = [len(c) for c in classes]
    assert cls[0] == 0 and sizes[0] == 1

    orders = [elem_order(F, elems, index, reps[k]) for k in range(r)]
    e = 1
    for o in orders:
        e = lcm(e, o)
    assert F.N % e == 0, f'{name}: exponent {e} does not divide N={F.N}'

    # class of g^j for each class rep
    powmap = []
    for k in range(r):
        x = elems[reps[k]]
        row = [0]
        cur = x
        for _ in range(1, orders[k]):
            row.append(cls[index[cur]])
            cur = mmul(F, cur, x)
        # row[j] = class of rep^j for j = 0..order-1 (row[0] -> identity class)
        row[0] = 0
        powmap.append(row)

    # inverse-class map
    invcls = [cls[index[minv(F, elems[reps[k]])]] for k in range(r)]

    # class multiplication coefficients a[i][j][k]
    a = [[[0] * r for _ in range(r)] for _ in range(r)]
    for k in range(r):
        z = elems[reps[k]]
        for xi, x in enumerate(elems):
            y = mmul(F, minv(F, x), z)
            a[cls[xi]][cls[index[y]]][k] += 1

    # choose the prime
    p = None
    cand = e + 1
    bound = 2 * isqrt(G) + 2
    while True:
        if cand > bound and is_prime(cand):
            p = cand
            break
        cand += e
    z = pow(primitive_root(p), (p - 1) // e, p)
    assert pow(z, e, p) == 1
    for q in set(_prime_factors(e)):
        assert pow(z, e // q, p) != 1

    Ms = [[[a[i][j][k] % p for k in range(r)] for j in range(r)] for i in range(r)]

    # simultaneous eigenvectors over F_p
    spaces = [[[1 if t == s else 0 for t in range(r)] for s in range(r)]]
    for i in range(r):
        newsp = []
        for B in spaces:
            k = len(B)
            if k == 1:
                newsp.append(B)
                continue
            Mi = Ms[i]
            MB = [[sum(Mi[rr][cc] * B[t][cc] for cc in range(r)) % p
                   for t in range(k)] for rr in range(r)]
            got = 0
            for lam in range(p):
                Cm = [[(MB[rr][t] - lam * B[t][rr]) % p for t in range(k)]
                      for rr in range(r)]
                ns = nullspace_mod(Cm, p, k)
                if ns:
                    nb = [[sum(B[t][rr] * v[t] for t in range(k)) % p
                           for rr in range(r)] for v in ns]
                    newsp.append(nb)
                    got += len(ns)
                    if got == k:
                        break
            assert got == k, 'eigenspace split incomplete'
        spaces = newsp
    assert all(len(B) == 1 for B in spaces) and len(spaces) == r, 'split failed'

    chars = []
    for B in spaces:
        v = B[0]
        iv = pow(v[0], p - 2, p)
        om = [(x * iv) % p for x in v]
        t = 0
        for kk in range(r):
            t = (t + om[kk] * om[invcls[kk]] * pow(sizes[kk], p - 2, p)) % p
        d2 = (G % p) * pow(t, p - 2, p) % p
        deg = None
        for cand_d in range(1, isqrt(G) + 1):
            if cand_d * cand_d % p == d2:
                deg = cand_d
                break
        assert deg is not None, 'degree recovery failed'
        chi_mod = [(deg * om[kk] * pow(sizes[kk], p - 2, p)) % p for kk in range(r)]
        assert chi_mod[0] == deg % p
        # lift class by class
        vals = []
        for kk in range(r):
            m = orders[kk]
            zm = pow(z, e // m, p)
            mults = []
            for l in range(m):
                s = 0
                for j in range(m):
                    s += chi_mod[powmap[kk][j]] * pow(zm, (-l * j) % m, p)
                s = s % p * pow(m, p - 2, p) % p
                mults.append(s)
            assert sum(mults) == deg, 'multiplicity lift failed'
            val = F.ZERO
            for l in range(m):
                if mults[l]:
                    val = F.add(val, F.scal(F.zeta((F.N // m) * l), mults[l]))
            vals.append(val)
        assert vals[0] == F.rat(deg)
        chars.append(vals)

    chars.sort(key=lambda v: (F.to_fraction(v[0]), [x for x in v[1:]]))

    gd = GroupData()
    gd.name = name
    gd.F = F
    gd.gens = gens
    gd.elems = elems
    gd.index = index
    gd.order = G
    gd.classes = classes
    gd.cls = cls
    gd.reps = reps
    gd.sizes = sizes
    gd.orders = orders
    gd.exponent = e
    gd.prime = p
    gd.powmap = powmap
    gd.invcls = invcls
    gd.chars = chars
    gd.nirr = r
    gd.chi_def = [mtr(F, elems[reps[k]]) for k in range(r)]
    gd.sq_class = [cls[index[mmul(F, elems[reps[k]], elems[reps[k]])]] for k in range(r)]
    return gd


def _prime_factors(n):
    out = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            out.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        out.append(n)
    return out


def inner(gd, u, v):
    """(1/|G|) sum_classes |C| u(g) conj(v(g)) -- exact."""
    F = gd.F
    acc = F.ZERO
    for k in range(gd.nirr):
        acc = F.add(acc, F.scal(F.mul(u[k], F.conj(v[k])), gd.sizes[k]))
    return F.scal(acc, Fraction(1, gd.order))


def verify_table(gd):
    F = gd.F
    r = gd.nirr
    ok = {'row_orthonormal': True, 'col_orthogonal': True, 'sum_dim_sq': None}
    for i in range(r):
        for j in range(r):
            v = inner(gd, gd.chars[i], gd.chars[j])
            want = F.ONE if i == j else F.ZERO
            if v != want:
                ok['row_orthonormal'] = False
    tot = 0
    for i in range(r):
        d = F.to_fraction(gd.chars[i][0])
        assert d.denominator == 1
        tot += int(d) ** 2
    ok['sum_dim_sq'] = tot
    assert tot == gd.order, f'{gd.name}: sum dim^2 = {tot} != {gd.order}'
    # column orthogonality
    for k in range(r):
        for l in range(r):
            acc = F.ZERO
            for i in range(r):
                acc = F.add(acc, F.mul(gd.chars[i][k], F.conj(gd.chars[i][l])))
            want = F.rat(Fraction(gd.order, gd.sizes[k])) if k == l else F.ZERO
            if acc != want:
                ok['col_orthogonal'] = False
    return ok

# =====================================================================
# PART 4 -- the two sums S, D and the tasks U1 .. U12
# =====================================================================


GROUP_NAMES = ([f'C_{n}' for n in range(2, 11)]
               + [f'BD_{n}' for n in range(2, 7)]
               + ['2T', '2O', '2I'])

OUT = {'meta': {
    'sheet': 'SPEC SHEET B',
    'groups': None,
    'arithmetic': ('exact: 2x2 matrices over Q(zeta_N) in a canonical power-basis '
                   'representation with Fraction coefficients; character tables by '
                   'Burnside-Dixon (class multiplication matrices, simultaneous '
                   'eigenvectors over F_p with p = 1 mod exponent, eigenvalue-'
                   'multiplicity lift), verified by row and column orthogonality '
                   'and sum of squared dimensions'),
    'no_imported_character_tables': True,
    'rational_format': 'numerator/denominator',
    'quadratic_format': 'a/b + (c/d)*sqrt(n)',
}}


def banner(t):
    print('\n' + '=' * 78)
    print(t)
    print('=' * 78)


# --------------------------------------------------------------- S and D
class GroupSD:
    def __init__(self, gd):
        F = gd.F
        self.gd = gd
        self.F = F
        r = gd.nirr
        self.dinv = [None] * r          # 1 / (2 - chi_def)
        self.cot2 = [None] * r          # cot^2(phi_g/2) = (2+chi)/(2-chi)
        for k in range(1, r):
            den = F.sub(F.rat(2), gd.chi_def[k])
            assert den != F.ZERO
            self.dinv[k] = F.inv(den)
            self.cot2[k] = F.mul(F.add(F.rat(2), gd.chi_def[k]), self.dinv[k])
        self.dims = [int(F.to_fraction(c[0])) for c in gd.chars]
        self.triv = [i for i in range(r) if all(v == F.ONE for v in gd.chars[i])]
        assert len(self.triv) == 1
        self.triv = self.triv[0]
        # index of the defining character among irreducibles (may be absent)
        self.defidx = None
        for i in range(r):
            if gd.chars[i] == gd.chi_def:
                self.defidx = i

    # class-sum route ---------------------------------------------------
    def D_class(self, chi):
        F, gd = self.F, self.gd
        acc = F.ZERO
        for k in range(1, gd.nirr):
            acc = F.add(acc, F.scal(F.mul(chi[k], self.dinv[k]), gd.sizes[k]))
        acc = F.scal(acc, Fraction(1, gd.order))
        return F.to_fraction(acc)

    def S_class(self, chi, dim):
        F, gd = self.F, self.gd
        acc = F.ZERO
        for k in range(1, gd.nirr):
            t = F.sub(chi[k], F.rat(dim))
            acc = F.add(acc, F.scal(F.mul(t, self.cot2[k]), gd.sizes[k]))
        acc = F.scal(acc, Fraction(1, gd.order))
        return F.to_fraction(acc)

    # element-sum route (independent check) ------------------------------
    def D_elem(self, chi):
        F, gd = self.F, self.gd
        acc = F.ZERO
        for i in range(1, gd.order):
            k = gd.cls[i]
            acc = F.add(acc, F.mul(chi[k], self.dinv[k]))
        acc = F.scal(acc, Fraction(1, gd.order))
        return F.to_fraction(acc)

    def S_elem(self, chi, dim):
        F, gd = self.F, self.gd
        acc = F.ZERO
        for i in range(1, gd.order):
            k = gd.cls[i]
            acc = F.add(acc, F.mul(F.sub(chi[k], F.rat(dim)), self.cot2[k]))
        acc = F.scal(acc, Fraction(1, gd.order))
        return F.to_fraction(acc)


def tensor_positivity(gd):
    """Strong check: every product chi_i * chi_j decomposes over the computed
    table with non-negative integer multiplicities of the right total degree.
    Also checks the defining character itself."""
    F = gd.F
    r = gd.nirr
    dims = [int(F.to_fraction(c[0])) for c in gd.chars]
    for i in range(r):
        for j in range(r):
            prod = [F.mul(gd.chars[i][k], gd.chars[j][k]) for k in range(r)]
            tot = 0
            for k in range(r):
                m = F.to_fraction(inner(gd, prod, gd.chars[k]))
                if m.denominator != 1 or m < 0:
                    return False
                tot += int(m) * dims[k]
            if tot != dims[i] * dims[j]:
                return False
    tot = 0
    for k in range(r):
        m = F.to_fraction(inner(gd, gd.chi_def, gd.chars[k]))
        if m.denominator != 1 or m < 0:
            return False
        tot += int(m) * dims[k]
    return tot == 2


def phi_over_pi(gd, k):
    """phi_g / pi as an exact rational, from 2 cos(phi) = chi_def."""
    F = gd.F
    tr = gd.chi_def[k]
    N = F.N
    target = F.numeric(tr).real
    best = None
    for j in range(0, N // 2 + 1):
        c = 2 * mp.cos(2 * mp.pi * j / N)
        if abs(c - target) < mp.mpf(10) ** (-40):
            best = Fraction(2 * j, N)
            break
    assert best is not None, 'phi not a rational multiple of pi'
    return best


def main():
    # ============================================================ build all
    banner('GROUP CONSTRUCTION AND CHARACTER TABLES')
    G = {}
    SD = {}
    OUT['groups'] = {}
    for nm in GROUP_NAMES:
        gd = build(nm)
        chk = verify_table(gd)
        assert chk['row_orthonormal'] and chk['col_orthogonal']
        tp = tensor_positivity(gd)
        assert tp, f'{nm}: tensor product positivity failed'
        G[nm] = gd
        SD[nm] = GroupSD(gd)
        dims = SD[nm].dims
        print(f'{nm:6s} |G|={gd.order:4d}  classes={gd.nirr:2d}  exponent={gd.exponent:3d}  '
              f'field=Q(zeta_{gd.F.N}) deg={gd.F.deg:2d}  Dixon p={gd.prime:3d}  '
              f'dims={dims}  sum d^2={chk["sum_dim_sq"]}  tensor-positivity={tp}')
        OUT['groups'][nm] = {
            'order': gd.order,
            'num_classes': gd.nirr,
            'exponent': gd.exponent,
            'field_N': gd.F.N,
            'field_degree': gd.F.deg,
            'dixon_prime': gd.prime,
            'irrep_dims': dims,
            'trivial_index': SD[nm].triv,
            'defining_char_index': SD[nm].defidx,
            'orthonormality_verified': True,
            'tensor_product_positivity_verified': tp,
            'sum_dim_squared': chk['sum_dim_sq'],
            'classes': [
                {'index': k, 'size': gd.sizes[k], 'element_order': gd.orders[k],
                 'chi_defining': render(gd.F, gd.chi_def[k])[1],
                 'phi_over_pi': fs(phi_over_pi(gd, k)) if k > 0 else '0/1'}
                for k in range(gd.nirr)],
            'character_table': [[render(gd.F, v)[1] for v in ch] for ch in gd.chars],
        }

    # ================================================================== U1
    banner('U1  --  S(alpha) and D(alpha) for every irreducible, every group')
    OUT['U1'] = {}
    for nm in GROUP_NAMES:
        gd, sd = G[nm], SD[nm]
        rows = []
        for i in range(gd.nirr):
            d = sd.dims[i]
            Dv = sd.D_class(gd.chars[i])
            Sv = sd.S_class(gd.chars[i], d)
            rows.append({'irrep': i, 'dim': d, 'S': fs(Sv), 'D': fs(Dv),
                         'is_trivial': i == sd.triv,
                         'is_defining': i == sd.defidx})
        D1 = Fraction(rows[sd.triv]['D'].split('/')[0]) / int(rows[sd.triv]['D'].split('/')[1])
        OUT['U1'][nm] = {'order': gd.order, 'rows': rows,
                         'S_trivial': rows[sd.triv]['S'],
                         'D_trivial': rows[sd.triv]['D']}
        print(f'\n{nm}  |G|={gd.order}   S(1)={rows[sd.triv]["S"]}   D(1)={rows[sd.triv]["D"]}')
        print('  ' + 'irrep'.ljust(6) + 'dim'.ljust(5) + 'S(alpha)'.ljust(22) + 'D(alpha)')
        for rw in rows:
            tag = ' <- trivial' if rw['is_trivial'] else (' <- defining' if rw['is_defining'] else '')
            print('  ' + str(rw['irrep']).ljust(6) + str(rw['dim']).ljust(5)
                  + rw['S'].ljust(22) + rw['D'] + tag)

    # ================================================================== U2
    banner('U2  --  csc^2(phi_g/2) = 4/(2 - chi_defining(g)) elementwise')
    maxdev = mp.mpf(0)
    worst = None
    exact_ok = True
    total_elems = 0
    for nm in GROUP_NAMES:
        gd, sd = G[nm], SD[nm]
        F = gd.F
        for k in range(1, gd.nirr):
            # exact algebraic side: csc^2 = cot^2 + 1
            lhs_exact = F.add(sd.cot2[k], F.ONE)
            rhs_exact = F.scal(sd.dinv[k], 4)
            if lhs_exact != rhs_exact:
                exact_ok = False
            # numeric trigonometric side
            ph = phi_over_pi(gd, k) * mp.pi
            csc2 = 1 / mp.sin(ph / 2) ** 2
            rhs = F.numeric(rhs_exact).real
            dev = abs(csc2 - rhs)
            total_elems += gd.sizes[k]
            if dev > maxdev:
                maxdev = dev
                worst = (nm, k, gd.sizes[k])
    print(f'exact identity csc^2 = cot^2 + 1 = 4/(2-chi) holds in the field for every class: {exact_ok}')
    print(f'elements g != I checked: {total_elems}')
    print(f'largest numeric deviation (mpmath, 60 dps): {mp.nstr(maxdev, 8)}   at {worst}')
    OUT['U2'] = {'exact_identity_holds': exact_ok,
                 'elements_checked': total_elems,
                 'largest_deviation_decimal': mp.nstr(maxdev, 12),
                 'largest_deviation_at': {'group': worst[0], 'class_index': worst[1],
                                          'class_size': worst[2]} if worst else None,
                 'working_precision_dps': mp.mp.dps}

    # ================================================================== U5
    banner('U5  --  (1/|G|) sum_{g!=I} cot^2(phi_g/2),  D(1),  D(1) + (1/8) sum')
    OUT['U5'] = {}
    print('  ' + 'group'.ljust(7) + '|G|'.ljust(6) + 'T = mean cot^2'.ljust(22)
          + 'D(1)'.ljust(16) + 'D(1)+T/8'.ljust(10) + 'comparison')
    for nm in GROUP_NAMES:
        gd, sd = G[nm], SD[nm]
        F = gd.F
        acc = F.ZERO
        for k in range(1, gd.nirr):
            acc = F.add(acc, F.scal(sd.cot2[k], gd.sizes[k]))
        T = F.to_fraction(F.scal(acc, Fraction(1, gd.order)))
        D1 = sd.D_class(gd.chars[sd.triv])
        comb = D1 + T / 8
        OUT['U5'][nm] = {'order': gd.order, 'mean_cot2': fs(T), 'D_trivial': fs(D1),
                         'D1_plus_eighth': fs(comb),
                         'num_classes': gd.nirr,
                         'equals_(num_classes-1)/8': comb == Fraction(gd.nirr - 1, 8),
                         'decimal': {'mean_cot2': float(T), 'D_trivial': float(D1),
                                     'D1_plus_eighth': float(comb)}}
        print('  ' + nm.ljust(7) + str(gd.order).ljust(6) + fs(T).ljust(22)
              + fs(D1).ljust(16) + fs(comb).ljust(10)
              + f'(#classes-1)/8 = {fs(Fraction(gd.nirr - 1, 8))}'
              + ('  MATCH' if comb == Fraction(gd.nirr - 1, 8) else '  NO MATCH'))

    # ================================================================= U12
    banner('U12  --  second route: element sum vs class sum')
    agree = True
    ndisagree = 0
    checks = 0
    per_group = {}
    for nm in GROUP_NAMES:
        gd, sd = G[nm], SD[nm]
        ok = True
        for i in range(gd.nirr):
            d = sd.dims[i]
            if (sd.D_class(gd.chars[i]) != sd.D_elem(gd.chars[i])
                    or sd.S_class(gd.chars[i], d) != sd.S_elem(gd.chars[i], d)):
                ok = False
                ndisagree += 1
            checks += 2
        per_group[nm] = ok
        agree = agree and ok
    print(f'{checks} exact comparisons (S and D, every irrep, every group).')
    print(f'all agree: {agree}   disagreements: {ndisagree}')
    OUT['U12'] = {'comparisons': checks, 'all_agree': agree,
                  'disagreements': ndisagree, 'per_group': per_group}

    return G, SD



# ------------------------------------------------------------------- helpers
def combo(gd, mult):
    """Character of the direct sum with multiplicity vector mult."""
    F = gd.F
    chi = [F.ZERO] * gd.nirr
    for i, m in enumerate(mult):
        if m:
            for k in range(gd.nirr):
                chi[k] = F.add(chi[k], F.scal(gd.chars[i][k], m))
    return chi


def task_U3(G, SD):
    banner('U3  --  is S(alpha) affine in (dim, D(alpha), dim*D(1)) for alpha with no trivial part?')
    data = []          # (group, label, dim, D, dim*D1, S)
    for nm in GROUP_NAMES:
        gd, sd = G[nm], SD[nm]
        D1 = sd.D_class(gd.chars[sd.triv])
        for i in range(gd.nirr):
            if i == sd.triv:
                continue
            d = sd.dims[i]
            Dv = sd.D_class(gd.chars[i])
            Sv = sd.S_class(gd.chars[i], d)
            data.append([nm, f'irr{i}', Fraction(d), Dv, Fraction(d) * D1, Sv])
        # extra check rows: sums of two non-trivial irreps (still no trivial part)
        nont = [i for i in range(gd.nirr) if i != sd.triv]
        for a, b in list(itertools.combinations_with_replacement(nont, 2))[:6]:
            mult = [0] * gd.nirr
            mult[a] += 1
            mult[b] += 1
            chi = combo(gd, mult)
            d = sd.dims[a] + sd.dims[b]
            Dv = sd.D_class(chi)
            Sv = sd.S_class(chi, d)
            data.append([nm, f'irr{a}+irr{b}', Fraction(d), Dv, Fraction(d) * D1, Sv])

    # pooled exact fit on a 3-row subset, then verify on the rest
    idx3 = None
    for trip in itertools.combinations(range(len(data)), 3):
        M = [[data[t][2], data[t][3], data[t][4]] for t in trip]
        det = (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
               - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
               + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))
        if det != 0:
            idx3 = trip
            break
    assert idx3 is not None, 'pooled design matrix has rank < 3'
    A3 = [[data[t][2], data[t][3], data[t][4]] for t in idx3]
    b3 = [data[t][5] for t in idx3]
    sol, null = solve_exact(A3, b3)
    a_u, b_u, c_u = sol
    print(f'fit subset rows {idx3}  ->  '
          f'({data[idx3[0]][0]},{data[idx3[0]][1]}), '
          f'({data[idx3[1]][0]},{data[idx3[1]][1]}), '
          f'({data[idx3[2]][0]},{data[idx3[2]][1]})')
    print(f'UNIVERSAL TRIPLE   a = {fs(a_u)}   b = {fs(b_u)}   c = {fs(c_u)}')

    resid = []
    nz = 0
    for row in data:
        nm, lab, d, Dv, dD1, Sv = row
        r = Sv - (a_u * d + b_u * Dv + c_u * dD1)
        resid.append({'group': nm, 'alpha': lab, 'dim': int(d), 'D': fs(Dv),
                      'S': fs(Sv), 'residual': fs(r)})
        if r != 0:
            nz += 1
    print(f'rows in the no-trivial-constituent pool: {len(data)}  '
          f'(irreps + 2-term sums)')
    print(f'nonzero residuals against the universal triple: {nz}')
    print(f'max |residual| = {fs(max((abs(Fraction(x["residual"])) for x in resid), default=Fraction(0)))}')

    # per-group analysis: the 3-column design is rank-deficient inside one group
    pergroup = {}
    print('\n  per-group fits (columns dim and dim*D(1) are proportional within a group):')
    print('  ' + 'group'.ljust(7) + 'rows'.ljust(6) + 'rank3'.ljust(7) + 'rank2'.ljust(7)
          + 'A = a + c*D(1)'.ljust(20) + 'b'.ljust(8) + 'univ. resid')
    for nm in GROUP_NAMES:
        rows = [r for r in data if r[0] == nm]
        M3 = [[r[2], r[3], r[4]] for r in rows]
        y = [r[5] for r in rows]
        rk3 = len(rref_rank(M3, 3))
        M2 = [[r[2], r[3]] for r in rows]
        rk2 = len(rref_rank(M2, 2))
        s2 = solve_exact(M2, y)
        D1 = SD[nm].D_class(G[nm].chars[SD[nm].triv])
        univ_res = max(abs(y[i] - (a_u * rows[i][2] + b_u * rows[i][3] + c_u * rows[i][4]))
                       for i in range(len(rows)))
        if s2 is None:
            Astr, bstr = 'INCONSISTENT', 'INCONSISTENT'
            part, nul = None, None
        else:
            part, nul = s2
            Astr = fs(part[0]) + ('' if not nul else ' (+free)')
            bstr = fs(part[1]) + ('' if not nul else ' (+free)')
        pergroup[nm] = {
            'rows': len(rows), 'rank_3col': rk3, 'rank_2col': rk2,
            'A_reduced': Astr, 'b_reduced': bstr,
            'nullspace_dim_2col': (len(nul) if nul is not None else None),
            'D_trivial': fs(D1),
            'a_plus_c_D1_from_universal': fs(a_u + c_u * D1),
            'universal_triple_max_residual': fs(univ_res),
        }
        print('  ' + nm.ljust(7) + str(len(rows)).ljust(6) + str(rk3).ljust(7)
              + str(rk2).ljust(7) + Astr.ljust(20) + bstr.ljust(8) + fs(univ_res))

    OUT['U3'] = {
        'universal_triple': {'a': fs(a_u), 'b': fs(b_u), 'c': fs(c_u)},
        'fit_subset_rows': [{'group': data[t][0], 'alpha': data[t][1]} for t in idx3],
        'pooled_rows': len(data),
        'pooled_nonzero_residuals': nz,
        'one_universal_triple_covers_every_group': nz == 0,
        'pooled_design_rank': 3,
        'residuals': resid,
        'per_group': pergroup,
        'note': ('within a single group the columns dim and dim*D(1) are '
                 'proportional, so the 3-column design has rank <= 2 there; '
                 'only the combination a + c*D(1) and b are determined per group.'),
    }
    return (a_u, b_u, c_u), data


def rref_rank(M, ncols):
    rows, piv = rref_frac([list(r) for r in M], ncols)
    return piv


def task_U4(G, SD, coeffs):
    a_u, b_u, c_u = coeffs
    banner('U4  --  alpha containing the trivial rep with multiplicity m: signed discrepancy')
    rows = []
    for nm in GROUP_NAMES:
        gd, sd = G[nm], SD[nm]
        D1 = sd.D_class(gd.chars[sd.triv])
        cases = []
        t = sd.triv
        mv = [0] * gd.nirr; mv[t] = 1
        cases.append(('1', list(mv)))
        mv = [0] * gd.nirr; mv[t] = 2
        cases.append(('1+1', list(mv)))
        for i in range(gd.nirr):
            mv = [0] * gd.nirr; mv[i] += 1; mv[t] += 1
            cases.append((f'irr{i}+1', list(mv)))
            mv = [0] * gd.nirr; mv[i] += 1; mv[t] += 2
            cases.append((f'irr{i}+1+1', list(mv)))
        for lab, mv in cases:
            chi = combo(gd, mv)
            dim = sum(mv[i] * sd.dims[i] for i in range(gd.nirr))
            m = mv[t]
            Dv = sd.D_class(chi)
            Sv = sd.S_class(chi, dim)
            rhs = a_u * Fraction(dim) + b_u * Dv + c_u * Fraction(dim) * D1
            rows.append({'group': nm, 'alpha': lab, 'm': m, 'dim': dim,
                         'S': fs(Sv), 'D': fs(Dv), 'rhs_U3': fs(rhs),
                         'discrepancy': fs(Sv - rhs)})
    bym = {}
    for r in rows:
        bym.setdefault(r['m'], set()).add(r['discrepancy'])
    print('  multiplicity m -> set of signed discrepancies S(alpha) - RHS_U3 observed:')
    for m in sorted(bym):
        print(f'    m = {m}: {sorted(bym[m])}   ({sum(1 for r in rows if r["m"] == m)} cases)')
    # the specifically requested case
    gd, sd = G['2I'], SD['2I']
    Pidx = sd.defidx
    mv = [0] * gd.nirr; mv[Pidx] += 1; mv[sd.triv] += 1
    chi = combo(gd, mv)
    dim = sd.dims[Pidx] + 1
    D1 = sd.D_class(gd.chars[sd.triv])
    Dv = sd.D_class(chi); Sv = sd.S_class(chi, dim)
    rhs = a_u * dim + b_u * Dv + c_u * dim * D1
    special = {'group': '2I', 'alpha': 'P + 1', 'm': 1, 'dim': dim,
               'S': fs(Sv), 'D': fs(Dv), 'rhs_U3': fs(rhs),
               'discrepancy': fs(Sv - rhs)}
    print(f'\n  requested case 2I, alpha = P + 1 : dim={dim}  S={fs(Sv)}  D={fs(Dv)}  '
          f'RHS={fs(rhs)}  discrepancy={fs(Sv - rhs)}  (m = 1)')
    OUT['U4'] = {'rows': rows, 'discrepancy_by_multiplicity':
                 {str(m): sorted(bym[m]) for m in sorted(bym)},
                 'special_case_2I_P_plus_1': special}



def lam2_char(gd, chi):
    """chi_{Lambda^2 rho}(g) = (chi(g)^2 - chi(g^2))/2."""
    F = gd.F
    return [F.scal(F.sub(F.mul(chi[k], chi[k]), chi[gd.sq_class[k]]), Fraction(1, 2))
            for k in range(gd.nirr)]


def sym2_char(gd, chi):
    """chi_{S^2 rho}(g) = (chi(g)^2 + chi(g^2))/2."""
    F = gd.F
    return [F.scal(F.add(F.mul(chi[k], chi[k]), chi[gd.sq_class[k]]), Fraction(1, 2))
            for k in range(gd.nirr)]


def task_U6(G, SD):
    banner('U6  --  every 2-dimensional irreducible rho with det rho = 1')
    res = {}
    for nm in GROUP_NAMES:
        gd, sd = G[nm], SD[nm]
        F = gd.F
        found = []
        for i in range(gd.nirr):
            if sd.dims[i] != 2:
                continue
            l2 = lam2_char(gd, gd.chars[i])
            if all(v == F.ONE for v in l2):
                found.append(i)
        res[nm] = found
        print(f'  {nm:6s}  2-dim irreps: {[i for i in range(gd.nirr) if sd.dims[i] == 2]}'
              f'   with det = 1: {found}   count = {len(found)}')
    # detail for 2I
    gd, sd = G['2I'], SD['2I']
    F = gd.F
    P = sd.defidx
    others = [i for i in res['2I'] if i != P]
    assert P in res['2I'] and len(others) == 1
    Pp = others[0]
    print(f'\n  2I:  P = irrep index {P} (character = trace of the defining matrices), '
          f'P-prime = irrep index {Pp}')
    print('  ' + 'class'.ljust(7) + 'size'.ljust(6) + 'order'.ljust(7)
          + 'chi_defining'.ljust(26) + 'chi_P'.ljust(26) + 'chi_P-prime')
    detail = []
    fields = set()
    for k in range(gd.nirr):
        a = render(F, gd.chars[P][k])
        b = render(F, gd.chars[Pp][k])
        c = render(F, gd.chi_def[k])
        for kind, s in (a, b):
            if kind.startswith('quadratic'):
                fields.add(s.split('sqrt(')[1].split(')')[0])
        detail.append({'class': k, 'size': gd.sizes[k], 'element_order': gd.orders[k],
                       'chi_defining': c[1], 'chi_P': a[1], 'chi_Pprime': b[1]})
        print('  ' + str(k).ljust(7) + str(gd.sizes[k]).ljust(6) + str(gd.orders[k]).ljust(7)
              + c[1].ljust(26) + a[1].ljust(26) + b[1])
    fieldstr = 'Q(sqrt(' + ', '.join(sorted(fields)) + '))' if fields else 'Q'
    print(f'  character values of P and P-prime generate: {fieldstr}   '
          f'(chi_P equals chi_defining on every class: '
          f'{gd.chars[P] == gd.chi_def})')
    OUT['U6'] = {'count_per_group': {nm: len(res[nm]) for nm in GROUP_NAMES},
                 'indices_per_group': {nm: res[nm] for nm in GROUP_NAMES},
                 '2I': {'P_index': P, 'Pprime_index': Pp, 'per_class': detail,
                        'field_generated': fieldstr,
                        'chi_P_equals_chi_defining': gd.chars[P] == gd.chi_def}}
    return P, Pp


def task_U7_U8_U9(G, SD, P, Pp):
    gd, sd = G['2I'], SD['2I']
    F = gd.F
    banner('U7  --  2I: S^2 P and S^2 P-prime; D and S for P, P-prime, S^2P, S^2P-prime')
    s2P = sym2_char(gd, gd.chars[P])
    s2Pp = sym2_char(gd, gd.chars[Pp])
    objs = [('P', gd.chars[P], 2), ('P-prime', gd.chars[Pp], 2),
            ('S^2P', s2P, 3), ('S^2P-prime', s2Pp, 3)]
    u7 = []
    for lab, chi, dim in objs:
        nrm = inner(gd, chi, chi)
        irr = nrm == F.ONE
        which = [i for i in range(gd.nirr) if gd.chars[i] == chi]
        d = F.to_fraction(chi[0])
        Dv = sd.D_class(chi)
        Sv = sd.S_class(chi, int(d))
        u7.append({'object': lab, 'dim': int(d), 'norm_sq': render(F, nrm)[1],
                   'irreducible': irr, 'irrep_index': (which[0] if which else None),
                   'D': fs(Dv), 'S': fs(Sv)})
        print(f'  {lab:11s} dim={int(d)}  <chi,chi>={render(F, nrm)[1]}  '
              f'irreducible={irr}  irrep index={which[0] if which else None}  '
              f'D={fs(Dv):12s} S={fs(Sv)}')
    OUT['U7'] = {'objects': u7}

    banner('U8  --  2I: class-by-class contributions to D for S^2P and S^2P-prime')
    rows = []
    agreeA = F.ZERO; agreeB = F.ZERO
    diffA = F.ZERO; diffB = F.ZERO
    diffclasses = []
    print('  ' + 'class'.ljust(6) + 'size'.ljust(6) + 'ord'.ljust(5)
          + 'chi_S2P'.ljust(24) + 'contrib D(S^2P)'.ljust(30) + 'contrib D(S^2P-prime)')
    for k in range(1, gd.nirr):
        cA = F.scal(F.mul(s2P[k], sd.dinv[k]), Fraction(gd.sizes[k], gd.order))
        cB = F.scal(F.mul(s2Pp[k], sd.dinv[k]), Fraction(gd.sizes[k], gd.order))
        same = (cA == cB)
        if same:
            agreeA = F.add(agreeA, cA); agreeB = F.add(agreeB, cB)
        else:
            diffA = F.add(diffA, cA); diffB = F.add(diffB, cB)
            diffclasses.append(k)
        rows.append({'class': k, 'size': gd.sizes[k], 'element_order': gd.orders[k],
                     'chi_defining': render(F, gd.chi_def[k])[1],
                     'chi_S2P': render(F, s2P[k])[1],
                     'chi_S2Pprime': render(F, s2Pp[k])[1],
                     'contrib_D_S2P': render(F, cA)[1],
                     'contrib_D_S2Pprime': render(F, cB)[1],
                     'agree': same})
        print('  ' + str(k).ljust(6) + str(gd.sizes[k]).ljust(6) + str(gd.orders[k]).ljust(5)
              + render(F, s2P[k])[1].ljust(24) + render(F, cA)[1].ljust(30)
              + render(F, cB)[1] + ('' if same else '   <-- differs'))
    DA = sd.D_class(s2P); DB = sd.D_class(s2Pp)
    SA = sd.S_class(s2P, 3); SB = sd.S_class(s2Pp, 3)
    print(f'\n  classes that AGREE: {[k for k in range(1, gd.nirr) if k not in diffclasses]}')
    print(f'  classes that DIFFER: {diffclasses}')
    print(f'  summed agreeing contribution   S^2P: {render(F, agreeA)[1]}   '
          f'S^2P-prime: {render(F, agreeB)[1]}')
    print(f'  summed differing contribution  S^2P: {render(F, diffA)[1]}   '
          f'S^2P-prime: {render(F, diffB)[1]}')
    print(f'  D(S^2P) = {fs(DA)}   D(S^2P-prime) = {fs(DB)}   '
          f'difference D(S^2P-prime) - D(S^2P) = {fs(DB - DA)}')
    print(f'  S(S^2P) = {fs(SA)}   S(S^2P-prime) = {fs(SB)}   '
          f'difference S(S^2P-prime) - S(S^2P) = {fs(SB - SA)}')
    OUT['U8'] = {'rows': rows,
                 'classes_agreeing': [k for k in range(1, gd.nirr) if k not in diffclasses],
                 'classes_differing': diffclasses,
                 'sum_agreeing_S2P': render(F, agreeA)[1],
                 'sum_agreeing_S2Pprime': render(F, agreeB)[1],
                 'sum_differing_S2P': render(F, diffA)[1],
                 'sum_differing_S2Pprime': render(F, diffB)[1],
                 'D_S2P': fs(DA), 'D_S2Pprime': fs(DB),
                 'D_diff': fs(DB - DA),
                 'S_S2P': fs(SA), 'S_S2Pprime': fs(SB),
                 'S_diff': fs(SB - SA)}

    banner('U9  --  2I: k(alpha) = dim(alpha)*D(1) - D(alpha)')
    D1 = sd.D_class(gd.chars[sd.triv])
    ks = {}
    print('  ' + 'object'.ljust(12) + 'dim'.ljust(5) + 'D(alpha)'.ljust(14)
          + 'k(alpha)'.ljust(16) + 'fractional part'.ljust(16) + 'decimal')
    u9 = []
    for lab, chi, _ in objs:
        d = int(F.to_fraction(chi[0]))
        Dv = sd.D_class(chi)
        kv = Fraction(d) * D1 - Dv
        frac = kv - (kv.numerator // kv.denominator)
        ks[lab] = kv
        u9.append({'object': lab, 'dim': d, 'D': fs(Dv), 'k': fs(kv),
                   'fractional_part': fs(frac), 'k_decimal': float(kv),
                   'frac_decimal': float(frac)})
        print('  ' + lab.ljust(12) + str(d).ljust(5) + fs(Dv).ljust(14)
              + fs(kv).ljust(16) + fs(frac).ljust(16) + f'{float(kv):.10f}')
    dPP = ks['P'] - ks['P-prime']
    dSS = ks['S^2P-prime'] - ks['S^2P']
    print(f'  k(P) - k(P-prime)          = {fs(dPP)}   = {float(dPP):.10f}')
    print(f'  k(S^2P-prime) - k(S^2P)    = {fs(dSS)}   = {float(dSS):.10f}')
    OUT['U9'] = {'D_trivial': fs(D1), 'rows': u9,
                 'k_P_minus_k_Pprime': fs(dPP),
                 'k_S2Pprime_minus_k_S2P': fs(dSS)}
    return s2P, s2Pp


def task_U10(G, SD, P, Pp):
    gd, sd = G['2I'], SD['2I']
    F = gd.F
    r = gd.nirr
    banner('U10  --  2I: McKay matrix A, null vector delta, and (2Id - A) H = e_P - e_Pprime')
    A = [[0] * r for _ in range(r)]
    for s in range(r):
        prod = [F.mul(gd.chi_def[k], gd.chars[s][k]) for k in range(r)]
        for t in range(r):
            v = inner(gd, gd.chars[t], prod)
            q = F.to_fraction(v)
            assert q.denominator == 1
            A[s][t] = int(q)
    print('  A (rows = sigma, cols = sigma-prime), irreps ordered by index:')
    print('        ' + ''.join(f'{t:4d}' for t in range(r)) + '     dim')
    for s in range(r):
        print(f'   {s:2d}   ' + ''.join(f'{A[s][t]:4d}' for t in range(r))
              + f'     {sd.dims[s]}')
    sym = all(A[s][t] == A[t][s] for s in range(r) for t in range(r))
    # BFS distances from the trivial node
    dist = [-1] * r
    dist[sd.triv] = 0
    frontier = [sd.triv]
    while frontier:
        nxt = []
        for x in frontier:
            for y in range(r):
                if A[x][y] and dist[y] < 0:
                    dist[y] = dist[x] + 1
                    nxt.append(y)
        frontier = nxt
    print(f'  A symmetric: {sym}    BFS distances from trivial node {sd.triv}: {dist}')
    M = [[(2 if s == t else 0) - A[s][t] for t in range(r)] for s in range(r)]
    # null vector, normalized at the trivial node
    sol = solve_exact(M, [0] * r)
    part, null = sol
    assert len(null) == 1, f'kernel dimension {len(null)}'
    dvec = null[0]
    dvec = [x / dvec[sd.triv] for x in dvec]
    print(f'  delta (null vector of 2Id - A, delta[trivial] = 1): '
          f'{[fs(x) for x in dvec]}')
    print(f'  delta equals the vector of irrep dimensions: '
          f'{[int(x) for x in dvec] == sd.dims}')
    rhs = [0] * r
    rhs[P] = 1
    rhs[Pp] = -1
    Aug = [row[:] for row in M] + [[1 if t == sd.triv else 0 for t in range(r)]]
    bb = list(rhs) + [0]
    solH = solve_exact(Aug, bb)
    assert solH is not None, 'system inconsistent'
    H, nullH = solH
    assert len(nullH) == 0, 'H not unique after normalization'
    ip = sum(H[i] * dvec[i] for i in range(r))
    print('  H (solution of (2Id - A) H = e_P - e_Pprime with H[trivial] = 0):')
    print('  ' + 'node'.ljust(6) + 'dim'.ljust(5) + 'dist'.ljust(6)
          + 'H entry'.ljust(14) + 'decimal')
    hrows = []
    for i in range(r):
        hrows.append({'node': i, 'dim': sd.dims[i], 'distance': dist[i],
                      'H': fs(H[i]), 'decimal': float(H[i]),
                      'is_P': i == P, 'is_Pprime': i == Pp,
                      'is_trivial': i == sd.triv})
        tag = ''
        if i == P: tag = '  <- P'
        if i == Pp: tag = '  <- P-prime'
        if i == sd.triv: tag = '  <- trivial'
        print('  ' + str(i).ljust(6) + str(sd.dims[i]).ljust(5) + str(dist[i]).ljust(6)
              + fs(H[i]).ljust(14) + f'{float(H[i]):.10f}' + tag)
    DP = sd.D_class(gd.chars[P]); DPp = sd.D_class(gd.chars[Pp])
    cmpv = Fraction(gd.order) * (DPp - DP)
    print(f'  <H, delta> = {fs(ip)}   = {float(ip):.10f}')
    print(f'  |Gamma| * (D(P-prime) - D(P)) = 120 * ({fs(DPp)} - {fs(DP)}) = {fs(cmpv)}')
    OUT['U10'] = {'A': A, 'A_symmetric': sym, 'distances': dist,
                  'delta': [fs(x) for x in dvec],
                  'delta_equals_dims': [int(x) for x in dvec] == sd.dims,
                  'H': hrows, 'H_inner_delta': fs(ip),
                  'order_times_D_diff': fs(cmpv),
                  'D_P': fs(DP), 'D_Pprime': fs(DPp)}



# ------------------------------------------------------------------ U11 / E8
E8_EDGES = [(0, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (1, 3)]


def e8_cartan_positive():
    A = [[0] * 8 for _ in range(8)]
    for i in range(8):
        A[i][i] = 2
    for i, j in E8_EDGES:
        A[i][j] = A[j][i] = -1
    return A


def det_int(M):
    n = len(M)
    A = [[Fraction(x) for x in row] for row in M]
    det = Fraction(1)
    for c in range(n):
        pr = None
        for r in range(c, n):
            if A[r][c] != 0:
                pr = r
                break
        if pr is None:
            return Fraction(0)
        if pr != c:
            A[c], A[pr] = A[pr], A[c]
            det = -det
        det *= A[c][c]
        pv = A[c][c]
        for r in range(c + 1, n):
            if A[r][c] != 0:
                f = A[r][c] / pv
                A[r] = [A[r][k] - f * A[c][k] for k in range(n)]
    return det


def enumerate_norm(Apos, target):
    """All integer c with c^T Apos c == target, by LDL^T (Fincke-Pohst)."""
    n = len(Apos)
    d = [Fraction(0)] * n
    L = [[Fraction(1) if i == j else Fraction(0) for j in range(n)] for i in range(n)]
    for i in range(n):
        s = Fraction(Apos[i][i])
        for k in range(i):
            s -= d[k] * L[i][k] ** 2
        d[i] = s
        assert d[i] > 0, 'form is not positive definite'
        for j in range(i + 1, n):
            t = Fraction(Apos[j][i])
            for k in range(i):
                t -= d[k] * L[j][k] * L[i][k]
            L[j][i] = t / d[i]
    sols = []
    c = [0] * n

    def rec(i, rem):
        if i < 0:
            if rem == 0:
                sols.append(tuple(c))
            return
        s = sum(L[j][i] * c[j] for j in range(i + 1, n))
        bnd = rem / d[i]
        b = float(bnd) ** 0.5
        lo = int(-float(s) - b) - 2
        hi = int(-float(s) + b) + 2
        for ci in range(lo, hi + 1):
            val = d[i] * (ci + s) ** 2
            if val <= rem:
                c[i] = ci
                rec(i - 1, rem - val)
        c[i] = 0

    rec(n - 1, Fraction(target))
    return sols


def task_U11():
    banner('U11  --  the E8 root lattice from its (negative definite) Cartan matrix')
    Apos = e8_cartan_positive()
    Gm = [[-Apos[i][j] for j in range(8)] for i in range(8)]
    dpos = det_int(Apos)
    dneg = det_int(Gm)
    print(f'  E8 Dynkin edges (0-indexed Bourbaki nodes 1..8): {E8_EDGES}')
    print(f'  det(positive Cartan) = {dpos}    det(negative form G = -Cartan) = {dneg}')
    print('  G (the negative definite Gram matrix used throughout):')
    for row in Gm:
        print('    ' + ' '.join(f'{v:3d}' for v in row))

    def form(x, y):
        return sum(x[i] * Gm[i][j] * y[j] for i in range(8) for j in range(8))

    roots = [tuple(c) for c in enumerate_norm(Apos, 2)]
    for c in roots:
        assert form(c, c) == -2
    print(f'\n  lattice vectors of norm -2 : {len(roots)}')
    normcount = {}
    for c in roots:
        normcount[form(c, c)] = normcount.get(form(c, c), 0) + 1
    # mod 2 fibres
    fib = {}
    for c in roots:
        fib.setdefault(tuple(v % 2 for v in c), []).append(c)
    sizes = sorted(set(len(v) for v in fib.values()))
    print(f'  distinct mod-2 classes among them : {len(fib)}')
    print(f'  fibre sizes : {sorted({s: sum(1 for v in fib.values() if len(v) == s) for s in sizes}.items())}'
          f'  (size -> how many classes)')
    negpairs = all(len(v) == 2 and set(v) == {v[0], tuple(-a for a in v[0])}
                   for v in fib.values())
    print(f'  every fibre is exactly a pair {{alpha, -alpha}} : {negpairs}')

    # mod 2 intersection form
    G2 = [[Gm[i][j] % 2 for j in range(8)] for i in range(8)]
    alternating = all(sum(x[i] * G2[i][j] * x[j] for i in range(8) for j in range(8)) % 2 == 0
                      for x in itertools.product((0, 1), repeat=8))
    diag_zero = all(G2[i][i] == 0 for i in range(8))
    rank2 = len(nullspace_mod([row[:] for row in G2], 2, 8))
    nondeg = (rank2 == 0)
    print(f'\n  mod-2 intersection form G mod 2:')
    for row in G2:
        print('    ' + ' '.join(str(v) for v in row))
    print(f'  alternating (x.x = 0 for all x in F_2^8, equivalently zero diagonal): '
          f'{alternating and diag_zero}')
    print(f'  nondegenerate (dim of radical = {rank2}, det mod 2 = {int(dneg) % 2}): {nondeg}')

    # mod 4 refinement
    def P4(x, lift_shift=None):
        xi = list(x)
        if lift_shift is not None:
            xi = [xi[i] + 2 * lift_shift[i] for i in range(8)]
        return form(xi, xi) % 4

    shifts = [tuple((7 * i + 3 * j) % 5 - 2 for j in range(8)) for i in range(6)]
    liftok = True
    for x in itertools.product((0, 1), repeat=8):
        base = P4(x)
        for sh in shifts:
            if P4(x, sh) != base:
                liftok = False
    counts = {}
    Pvals = {}
    for x in itertools.product((0, 1), repeat=8):
        v = P4(x)
        Pvals[x] = v
        counts[v] = counts.get(v, 0) + 1
    print(f'\n  mod-4 refinement P(x) = xi.xi mod 4, lift-independence verified on '
          f'{len(shifts)} shifts for all 256 classes: {liftok}')
    print(f'  multiset of P over the 256 classes: '
          f'{ {k: counts[k] for k in sorted(counts)} }')
    rootclasses = set(fib.keys())
    match = rootclasses == {x for x in Pvals if Pvals[x] == 2}
    print(f'  the P = 2 classes are exactly the 120 root classes: {match}')

    # reflections
    def refl(alpha, x):
        c = form(x, alpha)
        return tuple(x[i] + c * alpha[i] for i in range(8))

    rootset = set(roots)
    closed = all(refl(a, b) in rootset for a in roots for b in roots)
    print(f'\n  reflections s_alpha(x) = x + (x.alpha) alpha map roots to roots '
          f'(all {len(roots)}x{len(roots)} = {len(roots) ** 2} pairs): {closed}')
    gens = []
    for a in roots:
        M = []
        for i in range(8):
            ei = [1 if k == i else 0 for k in range(8)]
            M.append(tuple(v % 2 for v in refl(a, ei)))
        gens.append(M)
    uniq = {tuple(M) for M in gens}
    print(f'  distinct mod-2 reflection matrices from the 240 roots: {len(uniq)}')

    def act(M, x):
        out = [0] * 8
        for i in range(8):
            if x[i]:
                for j in range(8):
                    out[j] ^= M[i][j]
        return tuple(out)

    targets = {x for x in Pvals if Pvals[x] == 2}
    seen = set()
    orbits = []
    for x0 in sorted(targets):
        if x0 in seen:
            continue
        orb = {x0}
        stack = [x0]
        while stack:
            y = stack.pop()
            for M in uniq:
                z = act(M, y)
                if z not in orb:
                    orb.add(z)
                    stack.append(z)
        seen |= orb
        orbits.append(len(orb))
    # also the full action on all 256 classes, for context
    seen2 = set()
    orb_all = []
    for x0 in itertools.product((0, 1), repeat=8):
        if x0 in seen2:
            continue
        orb = {x0}
        stack = [x0]
        while stack:
            y = stack.pop()
            for M in uniq:
                z = act(M, y)
                if z not in orb:
                    orb.add(z)
                    stack.append(z)
        seen2 |= orb
        orb_all.append(len(orb))
    print(f'  orbits of the reflection group on the {len(targets)} classes with P = 2: '
          f'{len(orbits)} orbit(s), sizes {sorted(orbits, reverse=True)}')
    print(f'  (for context) orbits on all 256 classes: {len(orb_all)}, '
          f'sizes {sorted(orb_all, reverse=True)}')

    OUT['U11'] = {
        'dynkin_edges_0indexed': E8_EDGES,
        'det_positive_cartan': str(dpos),
        'det_negative_form': str(dneg),
        'gram_negative': Gm,
        'norm_minus2_vectors': len(roots),
        'distinct_mod2_classes': len(fib),
        'fibre_sizes': {str(s): sum(1 for v in fib.values() if len(v) == s) for s in sizes},
        'every_fibre_is_plus_minus_pair': negpairs,
        'reflection_closure_verified_pairs': len(roots) ** 2,
        'reflections_preserve_root_set': closed,
        'mod2_form': G2,
        'mod2_alternating': bool(alternating and diag_zero),
        'mod2_radical_dim': rank2,
        'mod2_nondegenerate': nondeg,
        'mod4_lift_independent': liftok,
        'mod4_value_counts': {str(k): counts[k] for k in sorted(counts)},
        'P2_classes_equal_root_classes': match,
        'reflection_matrices_mod2_distinct': len(uniq),
        'orbits_on_P2_classes': {'count': len(orbits), 'sizes': sorted(orbits, reverse=True)},
        'orbits_on_all_256': {'count': len(orb_all), 'sizes': sorted(orb_all, reverse=True)},
    }


if __name__ == '__main__':
    G, SD = main()
    COEF, U3DATA = task_U3(G, SD)
    task_U4(G, SD, COEF)
    P, Pp = task_U6(G, SD)
    task_U7_U8_U9(G, SD, P, Pp)
    task_U10(G, SD, P, Pp)
    task_U11()
    OUT['meta']['groups'] = GROUP_NAMES
    with open('m8_1_1_defect.json', 'w') as fh:
        json.dump(OUT, fh, indent=1, sort_keys=False)
    print('\nwrote m8_1_1_defect.json')
