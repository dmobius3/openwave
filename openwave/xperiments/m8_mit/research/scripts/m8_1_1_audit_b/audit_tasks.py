#!/usr/bin/env python3
"""Task-by-task independent recomputation (U1..U12) plus the widened
universality search."""

from fractions import Fraction as Fr
from math import gcd
import itertools
import random
import mpmath as mp

from audit_core import (lin_solve, null_space, mat_rank, divisors, Tsum,
                        prim_T, ramanujan, CycloField)
from audit_groups import (make_group, verify_table, inner, cf_mul, cf_add,
                          cf_scal, cf_sub, cyclic_subgroups)
from audit_sums import (SD, D_trivial_cyclic, mckay_matrix, D_all_mckay,
                        D_closed_form_C, D_closed_form_BD, field_numeric)

mp.mp.dps = 80


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
