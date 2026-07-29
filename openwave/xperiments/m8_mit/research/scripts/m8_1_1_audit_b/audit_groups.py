#!/usr/bin/env python3
"""Group construction + conjugacy classes + character tables (induced
characters from cyclic subgroups, orthogonality-driven peeling)."""

from fractions import Fraction as Fr
from math import gcd
import itertools

from audit_core import (CycloField, Quad, qmul, qnorm, lin_solve, null_space,
                        mat_rank, divisors, Tsum, prim_T, ramanujan)


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
