#!/usr/bin/env python3
"""
m8_1_1_coexact_solver.py

Self-contained execution of SPEC SHEET A.

Everything (groups, conjugacy classes, irreducible characters, branchings,
invariant dimensions) is built from the explicit 2x2 generator matrices given
in the spec.  Nothing is imported from a stored character table or group
library.  numpy / mpmath / fractions / stdlib only.

Numerical policy
----------------
* group combinatorics (multiplication table, classes, power maps) uses
  float64 matrices; element identification is exact in practice because the
  minimum pairwise distance between distinct group elements is O(0.1) while
  float64 round-off is O(1e-15).  Both facts are measured and reported.
* character theory uses mpmath at 50 decimal digits.  Every quantity that
  theory forces to be an integer is rounded, and the maximum deviation from
  integrality over the whole run is reported.
* an INDEPENDENT EXACT check of the branching multiplicities is done with
  integer arithmetic in the cyclotomic ring Z[zeta_e] (e = group exponent):
  the Gram matrix G[a][b] = <chi_a, chi_b> is computed exactly and compared
  with M M^T, M[a][sigma] = multiplicity of sigma in V_a.
* a third, fully independent method (explicit averaging projectors + SVD
  rank) is used for the invariant dimensions of task T3.
"""

import json
import math
import os
import time
from fractions import Fraction
from itertools import product

import numpy as np
import mpmath as mp

mp.mp.dps = 50

TOL_ELEM = 1e-8          # tolerance for identifying two group elements
SVD_TOL = 1e-8           # stated tolerance for the numerical-rank decision
JMAX = 14                # symmetric powers V_0 .. V_14
MMIN, MMAX = 2, 12       # E_m for m = 2..12

AUDIT = {
    "max_integer_rounding_deviation": 0.0,
    "max_integer_rounding_where": None,
    "min_pairwise_element_distance": None,
    "element_lookup_dict_misses": 0,
}


# ----------------------------------------------------------------------------
# small helpers
# ----------------------------------------------------------------------------
def note_int(x, where):
    """Round a high-precision value that theory forces to be an integer."""
    v = float(mp.re(x)) if isinstance(x, (mp.mpf, mp.mpc)) else float(np.real(x))
    n = int(round(v))
    dev = abs(v - n)
    if isinstance(x, (mp.mpf, mp.mpc)):
        dev = max(dev, abs(float(mp.im(x))))
    else:
        dev = max(dev, abs(float(np.imag(x))))
    if dev > AUDIT["max_integer_rounding_deviation"]:
        AUDIT["max_integer_rounding_deviation"] = dev
        AUDIT["max_integer_rounding_where"] = where
    if dev > 1e-6:
        raise RuntimeError("value not integral at %s: %r (dev %g)" % (where, x, dev))
    return n


def fmt(x, digits=30):
    return mp.nstr(x, digits, strip_zeros=False)


def as_rational_or_none(x, tol=mp.mpf(10) ** (-35), maxden=10 ** 6):
    """Return 'p/q' if x is (numerically) an exact rational, else None."""
    if abs(mp.im(x)) > tol:
        return None
    r = mp.re(x)
    fr = Fraction(float(r)).limit_denominator(maxden)
    if abs(mp.mpf(fr.numerator) / fr.denominator - r) < tol:
        return "%d/%d" % (fr.numerator, fr.denominator)
    return None


# ----------------------------------------------------------------------------
# 1. groups
# ----------------------------------------------------------------------------
def quat_to_mat(a, b, c, d):
    """q = a + b i + c j + d k  ->  [[a+bi, c+di], [-c+di, a-bi]]."""
    return np.array([[a + b * 1j, c + d * 1j],
                     [-c + d * 1j, a - b * 1j]], dtype=complex)


def close_group(gens, cap=4000):
    ident = np.eye(2, dtype=complex)
    elems = [ident]
    frontier = [ident]
    while frontier:
        new = []
        for M in frontier:
            for g in gens:
                P = M @ g
                if any(np.max(np.abs(P - E)) < TOL_ELEM for E in elems):
                    continue
                if any(np.max(np.abs(P - E)) < TOL_ELEM for E in new):
                    continue
                new.append(P)
        elems.extend(new)
        frontier = new
        if len(elems) > cap:
            raise RuntimeError("closure blew past cap")
    return elems


PHI = (1.0 + math.sqrt(5.0)) / 2.0


def explicit_2I_quaternions():
    """The fallback explicit 120-element list from the spec sheet."""
    qs = []
    for s in (1, -1):
        for pos in range(4):
            v = [0.0, 0.0, 0.0, 0.0]
            v[pos] = float(s)
            qs.append(tuple(v))
    for s0 in (1, -1):
        for s1 in (1, -1):
            for s2 in (1, -1):
                for s3 in (1, -1):
                    qs.append((s0 * 0.5, s1 * 0.5, s2 * 0.5, s3 * 0.5))
    # even permutations of (0, +-1, +-phi^-1, +-phi)/2
    even_perms = [p for p in _perms4() if _parity(p) == 0]
    base = [0.0, 1.0, 1.0 / PHI, PHI]
    for p in even_perms:
        for s1 in (1, -1):
            for s2 in (1, -1):
                for s3 in (1, -1):
                    vals = [0.0, s1 * base[1], s2 * base[2], s3 * base[3]]
                    q = [0.0] * 4
                    for src in range(4):
                        q[p[src]] = vals[src] / 2.0
                    qs.append(tuple(q))
    # dedupe
    out = []
    for q in qs:
        if not any(max(abs(q[i] - r[i]) for i in range(4)) < 1e-9 for r in out):
            out.append(q)
    return out


def _perms4():
    from itertools import permutations
    return list(permutations(range(4)))


def _parity(p):
    n = 0
    for i in range(4):
        for j in range(i + 1, 4):
            if p[i] > p[j]:
                n += 1
    return n % 2


def build_groups():
    groups = []
    for n in range(1, 11):
        z = np.exp(2j * math.pi / n)
        gens = [np.array([[z, 0], [0, np.conj(z)]], dtype=complex)]
        groups.append(("C_%d" % n, gens, n))
    for n in range(2, 7):
        z = np.exp(1j * math.pi / n)
        gens = [np.array([[z, 0], [0, np.conj(z)]], dtype=complex),
                np.array([[0, 1], [-1, 0]], dtype=complex)]
        groups.append(("BD_%d" % n, gens, 4 * n))
    g2T = [quat_to_mat(0.5, 0.5, 0.5, 0.5), quat_to_mat(0, 1, 0, 0)]
    groups.append(("2T", g2T, 24))
    g2O = g2T + [quat_to_mat(1 / math.sqrt(2), 1 / math.sqrt(2), 0, 0)]
    groups.append(("2O", g2O, 48))
    g2I = [quat_to_mat(0.5, 0.5, 0.5, 0.5),
           quat_to_mat(PHI / 2, (1 / PHI) / 2, 0.5, 0.0)]
    groups.append(("2I", g2I, 120))
    return groups


# ----------------------------------------------------------------------------
# 2. symmetric powers
# ----------------------------------------------------------------------------
def sym_power_monomial(g, a):
    """Sym^a(g) in the monomial basis x^(a-k) y^k (columns indexed by k)."""
    al, be, ga, de = g[0, 0], g[0, 1], g[1, 0], g[1, 1]
    S = np.zeros((a + 1, a + 1), dtype=complex)
    for k in range(a + 1):
        # (al x + ga y)^(a-k) * (be x + de y)^k
        for p in range(a - k + 1):
            cp = math.comb(a - k, p) * al ** (a - k - p) * ga ** p
            for qq in range(k + 1):
                cq = math.comb(k, qq) * be ** (k - qq) * de ** qq
                S[p + qq, k] += cp * cq
    return S


def sym_power_unitary(g, a):
    """Sym^a(g) in the orthonormal basis sqrt(C(a,k)) x^(a-k) y^k."""
    S = sym_power_monomial(g, a)
    d = np.array([math.sqrt(math.comb(a, k)) for k in range(a + 1)])
    return (S * d[None, :]) / d[:, None]


# ----------------------------------------------------------------------------
# 3. exact cyclotomic arithmetic in Z[x]/(x^e - 1), reduce mod Phi_e to test
# ----------------------------------------------------------------------------
_PHI_CACHE = {}


def poly_divmod_monic(num, den):
    """Exact division of integer polys, den monic; returns (quot, rem)."""
    num = list(num)
    dn = len(den) - 1
    q = [0] * max(0, len(num) - dn)
    for i in range(len(num) - 1, dn - 1, -1):
        c = num[i]
        if c:
            q[i - dn] = c
            for j in range(dn + 1):
                num[i - dn + j] -= c * den[j]
    while len(num) > 1 and num[-1] == 0:
        num.pop()
    return q, num


def cyclotomic_poly(n):
    if n in _PHI_CACHE:
        return _PHI_CACHE[n]
    num = [-1] + [0] * (n - 1) + [1]          # x^n - 1
    for d in range(1, n):
        if n % d == 0:
            num, r = poly_divmod_monic(num, cyclotomic_poly(d))
            assert all(c == 0 for c in r), (n, d)
    while len(num) > 1 and num[-1] == 0:
        num.pop()
    _PHI_CACHE[n] = num
    return num


class Cyc(object):
    """Element of Z[x]/(x^e - 1) with integer coefficients."""
    __slots__ = ("e", "c")

    def __init__(self, e, c=None):
        self.e = e
        self.c = [0] * e if c is None else list(c)

    @staticmethod
    def zeta(e, k):
        z = Cyc(e)
        z.c[k % e] = 1
        return z

    def __add__(self, o):
        return Cyc(self.e, [a + b for a, b in zip(self.c, o.c)])

    def __mul__(self, o):
        e = self.e
        out = [0] * e
        for i, a in enumerate(self.c):
            if a:
                for j, b in enumerate(o.c):
                    if b:
                        out[(i + j) % e] += a * b
        return Cyc(e, out)

    def scale(self, n):
        return Cyc(self.e, [n * a for a in self.c])

    def conj(self):
        e = self.e
        out = [0] * e
        for i, a in enumerate(self.c):
            out[(-i) % e] += a
        return Cyc(e, out)

    def to_rational(self):
        """Reduce mod Phi_e; return Fraction if the result is rational."""
        e = self.e
        r = list(self.c)
        while len(r) > 1 and r[-1] == 0:
            r.pop()
        _, rem = poly_divmod_monic(r, cyclotomic_poly(e))
        # rem is the representative in Z[x]/Phi_e; rational iff degree 0
        while len(rem) > 1 and rem[-1] == 0:
            rem.pop()
        if len(rem) == 1:
            return Fraction(rem[0], 1)
        return None


# ----------------------------------------------------------------------------
# 4. per-group analysis
# ----------------------------------------------------------------------------
class GroupData(object):
    pass


def analyse_group(name, gens, order_expected, log):
    G = GroupData()
    G.name = name
    t0 = time.time()

    elems = close_group(gens)
    fallback_used = False
    if len(elems) != order_expected and name == "2I":
        log.append("2I closure gave %d, falling back to explicit list" % len(elems))
        qs = explicit_2I_quaternions()
        elems = [quat_to_mat(*q) for q in qs]
        fallback_used = True
    G.explicit_list_agrees = None
    if name == "2I":
        ex = [quat_to_mat(*q) for q in explicit_2I_quaternions()]
        matched = 0
        for M in ex:
            if any(np.max(np.abs(M - X)) < TOL_ELEM for X in elems):
                matched += 1
        G.explicit_list_agrees = (matched == len(ex) == len(elems))
        log.append("2I: closure has %d elements, explicit spec list has %d, "
                   "%d of the explicit elements found in the closure"
                   % (len(elems), len(ex), matched))
    G.order = len(elems)
    G.order_expected = order_expected
    G.order_ok = (G.order == order_expected)
    G.fallback_used = fallback_used
    if not G.order_ok:
        log.append("%s: ORDER MISMATCH built %d expected %d"
                   % (name, G.order, order_expected))
    E = np.array(elems)                        # (N,2,2)
    N = G.order

    # unitarity / determinant / minimum separation
    G.max_unitarity_err = float(max(
        np.max(np.abs(E[i].conj().T @ E[i] - np.eye(2))) for i in range(N)))
    G.max_det_err = float(max(
        abs(E[i][0, 0] * E[i][1, 1] - E[i][0, 1] * E[i][1, 0] - 1.0)
        for i in range(N)))
    mind = min(float(np.max(np.abs(E[i] - E[j])))
               for i in range(N) for j in range(i + 1, N)) if N > 1 else float("inf")
    G.min_sep = mind
    if AUDIT["min_pairwise_element_distance"] is None or mind < AUDIT["min_pairwise_element_distance"]:
        AUDIT["min_pairwise_element_distance"] = mind

    # index lookup
    def key(M):
        return tuple(int(round(v * 1e6)) for v in
                     np.concatenate([M.real.ravel(), M.imag.ravel()]))

    kmap = {}
    for i in range(N):
        kmap.setdefault(key(E[i]), i)
    G.key_injective = (len(kmap) == N)

    def idx(M):
        k = key(M)
        if k in kmap:
            return kmap[k]
        AUDIT["element_lookup_dict_misses"] += 1
        best, bd = -1, 1e9
        for i in range(N):
            d = float(np.max(np.abs(E[i] - M)))
            if d < bd:
                bd, best = d, i
        if bd > TOL_ELEM:
            raise RuntimeError("element not found, dist %g" % bd)
        return best

    # multiplication table + closure verification
    mul = np.zeros((N, N), dtype=int)
    for i in range(N):
        for j in range(N):
            mul[i, j] = idx(E[i] @ E[j])
    G.closed = True
    # closure is verified by idx() succeeding for every product
    inv = np.zeros(N, dtype=int)
    for i in range(N):
        r = np.where(mul[i] == 0)[0]
        assert len(r) == 1
        inv[i] = r[0]

    # element orders
    eorder = np.zeros(N, dtype=int)
    for i in range(N):
        c, o = i, 1
        while c != 0:
            c = mul[c, i]
            o += 1
        eorder[i] = o

    # -I present?
    minusI = -np.eye(2, dtype=complex)
    G.has_minus_I = any(np.max(np.abs(E[i] - minusI)) < TOL_ELEM for i in range(N))

    # ---- conjugacy classes (3.1) -------------------------------------------
    seen = [False] * N
    classes = []
    for i in range(N):
        if seen[i]:
            continue
        orb = set()
        for g in range(N):
            orb.add(mul[mul[g, i], inv[g]])
        for x in orb:
            seen[x] = True
        classes.append(sorted(orb))
    # class sort key: (element order, -trace, size, min index)
    def cls_key(c):
        r = c[0]
        return (int(eorder[r]), -round(float(np.real(np.trace(E[r]))), 9), len(c), min(c))
    classes.sort(key=cls_key)
    # identity class first
    for t, c in enumerate(classes):
        if 0 in c:
            classes.insert(0, classes.pop(t))
            break
    G.classes = classes
    nc = len(classes)
    G.nc = nc
    cls_of = np.zeros(N, dtype=int)
    for ci, c in enumerate(classes):
        for x in c:
            cls_of[x] = ci
    G.cls_of = cls_of
    G.class_sizes = [len(c) for c in classes]
    G.class_reps = [c[0] for c in classes]
    G.class_orders = [int(eorder[c[0]]) for c in classes]
    G.class_traces = [float(np.real(np.trace(E[c[0]]))) for c in classes]

    # eigenvalue lambda = e^{2 pi i m / n} for each class (exact data)
    G.class_lambda = []
    for ci, c in enumerate(classes):
        n = G.class_orders[ci]
        tr = G.class_traces[ci]
        best, bd = None, 1e9
        for m in range(n):
            if math.gcd(m, n) != 1 and not (n == 1 and m == 0):
                continue
            d = abs(2 * math.cos(2 * math.pi * m / n) - tr)
            if d < bd:
                bd, best = d, m
        G.class_lambda.append((best, n))
        if bd > 1e-7:
            log.append("%s: lambda id failure class %d dev %g" % (name, ci, bd))
    G.exponent = 1
    for n in G.class_orders:
        G.exponent = G.exponent * n // math.gcd(G.exponent, n)

    # power map: class of g^p
    def class_of_power(ci, p):
        g = G.class_reps[ci]
        c = 0
        for _ in range(p):
            c = mul[c, g]
        return int(cls_of[c])
    G.class_sq = [class_of_power(ci, 2) for ci in range(nc)]

    # ---- symmetric-power characters chi_a on each class ---------------------
    # exact via lambda; high precision via mpmath
    chi = mp.zeros(JMAX + 1, nc)
    for ci in range(nc):
        m, n = G.class_lambda[ci]
        th = 2 * mp.pi * mp.mpf(m) / n
        lam = mp.exp(mp.mpc(0, 1) * th)
        for a in range(JMAX + 1):
            s = mp.mpc(0)
            for k in range(a + 1):
                s += lam ** (a - 2 * k)
            chi[a, ci] = s
    G.chi_sym = chi

    # ---- irreducible characters (3.2) --------------------------------------
    # class multiplication matrices, exact integers.
    # a_ijk = #{(x,y) : x in C_i, y in C_j, x*y = g_k} for the FIXED
    # representative g_k, computed as #{x in C_i : x^{-1} g_k in C_j}.
    # (M_i)_{k,j} = a_ijk represents multiplication by the class sum C_i^.
    Ms = [np.zeros((nc, nc), dtype=object) for _ in range(nc)]
    for i in range(nc):
        for k in range(nc):
            gk = G.class_reps[k]
            for x in classes[i]:
                j = cls_of[mul[inv[x], gk]]
                Ms[i][k, j] += 1
    G.class_matrices = [[[int(Ms[i][k, j]) for j in range(nc)] for k in range(nc)]
                        for i in range(nc)]

    # simultaneous eigenvectors via a random integer combination
    rng = np.random.RandomState(12345 + len(name))
    for attempt in range(60):
        r = rng.randint(1, 40, size=nc)
        Mc = sum(int(r[i]) * Ms[i].astype(object) for i in range(nc))
        Mm = mp.matrix([[mp.mpf(int(Mc[a, b])) for b in range(nc)] for a in range(nc)])
        Ev, ER = mp.eig(Mm)
        ok = True
        for a in range(nc):
            for b in range(a + 1, nc):
                if abs(Ev[a] - Ev[b]) < mp.mpf(10) ** (-20):
                    ok = False
        if ok:
            break
    if not ok:
        raise RuntimeError("could not separate eigenvalues for %s" % name)
    G.eig_attempts = attempt + 1

    Mmats = [mp.matrix([[mp.mpf(int(Ms[i][a, b])) for b in range(nc)]
                        for a in range(nc)]) for i in range(nc)]
    rows = []
    omega_consistency = mp.mpf(0)
    for t in range(nc):
        v = mp.matrix([ER[i, t] for i in range(nc)])
        kbest = max(range(nc), key=lambda i: abs(v[i]))
        vmax = abs(v[kbest])
        omega = []
        for i in range(nc):
            w = Mmats[i] * v
            om = w[kbest] / v[kbest]
            for k2 in range(nc):
                if abs(v[k2]) > vmax / 10:
                    omega_consistency = max(omega_consistency,
                                            abs(w[k2] / v[k2] - om))
            omega.append(om)
        ssum = mp.mpc(0)
        for i in range(nc):
            ssum += (omega[i] * mp.conj(omega[i])) / G.class_sizes[i]
        deg = mp.sqrt(mp.mpf(N) / mp.re(ssum))
        degi = note_int(deg, "%s degree" % name)
        row = [omega[i] * degi / G.class_sizes[i] for i in range(nc)]
        rows.append((degi, row))

    # canonical ordering: trivial first, then by (dim, numeric key)
    def is_trivial(row):
        return all(abs(x - 1) < mp.mpf(10) ** (-20) for x in row)
    triv = [t for t, (d, row) in enumerate(rows) if is_trivial(row)]
    assert len(triv) == 1, (name, len(triv))
    ti = triv[0]
    others = [t for t in range(nc) if t != ti]
    others.sort(key=lambda t: (rows[t][0],
                               tuple(round(float(mp.re(x)), 12) for x in rows[t][1]),
                               tuple(round(float(mp.im(x)), 12) for x in rows[t][1])))
    orderidx = [ti] + others
    G.dims = [rows[t][0] for t in orderidx]
    G.chartab = [rows[t][1] for t in orderidx]

    # verifications
    err_row = mp.mpf(0)
    for a in range(nc):
        for b in range(nc):
            s = mp.mpc(0)
            for i in range(nc):
                s += G.class_sizes[i] * G.chartab[a][i] * mp.conj(G.chartab[b][i])
            tgt = mp.mpf(N) if a == b else mp.mpf(0)
            err_row = max(err_row, abs(s - tgt))
    err_col = mp.mpf(0)
    for i in range(nc):
        for j in range(nc):
            s = mp.mpc(0)
            for a in range(nc):
                s += G.chartab[a][i] * mp.conj(G.chartab[a][j])
            tgt = mp.mpf(N) / G.class_sizes[i] if i == j else mp.mpf(0)
            err_col = max(err_col, abs(s - tgt))
    G.omega_consistency = float(omega_consistency)
    G.orth_row_err = float(err_row)
    G.orth_col_err = float(err_col)
    G.sum_dim_sq = sum(d * d for d in G.dims)
    G.sum_dim_sq_ok = (G.sum_dim_sq == N)

    # ---- inner products -----------------------------------------------------
    def ip(f, g_):
        """<f, g> = (1/|G|) sum_gamma f(gamma) conj(g(gamma)) over classes."""
        s = mp.mpc(0)
        for i in range(nc):
            s += G.class_sizes[i] * f[i] * mp.conj(g_[i])
        return s / N
    G.ip = ip

    # multiplicity matrix M[a][sigma] = <chi_sigma, chi_a>
    Mult = [[0] * nc for _ in range(JMAX + 1)]
    for a in range(JMAX + 1):
        fa = [chi[a, i] for i in range(nc)]
        for s in range(nc):
            Mult[a][s] = note_int(ip(fa, G.chartab[s]),
                                  "%s mult V_%d sigma%d" % (name, a, s))
    G.Mult = Mult
    G.dim_check = all(sum(Mult[a][s] * G.dims[s] for s in range(nc)) == a + 1
                      for a in range(JMAX + 1))

    # ---- adjacency (3.3) ----------------------------------------------------
    chi1 = [chi[1, i] for i in range(nc)]
    A = [[0] * nc for _ in range(nc)]
    for s in range(nc):
        prod = [chi1[i] * G.chartab[s][i] for i in range(nc)]
        for s2 in range(nc):
            A[s][s2] = note_int(ip(prod, G.chartab[s2]), "%s A[%d][%d]" % (name, s, s2))
    G.A = A
    G.A_symmetric = all(A[i][j] == A[j][i] for i in range(nc) for j in range(nc))
    G.A_rowsum_ok = all(sum(A[s][s2] * G.dims[s2] for s2 in range(nc)) == 2 * G.dims[s]
                        for s in range(nc))

    # ---- distances (3.4) ----------------------------------------------------
    dist = [-1] * nc
    dist[0] = 0
    frontier = [0]
    while frontier:
        nxt = []
        for u in frontier:
            for v2 in range(nc):
                if A[u][v2] != 0 and dist[v2] < 0:
                    dist[v2] = dist[u] + 1
                    nxt.append(v2)
        frontier = nxt
    G.dist = dist
    G.connected = all(d >= 0 for d in dist)
    # graph diameter (all-pairs BFS)
    diam = 0
    for src in range(nc):
        dd = [-1] * nc
        dd[src] = 0
        fr = [src]
        while fr:
            nx = []
            for u in fr:
                for v2 in range(nc):
                    if A[u][v2] != 0 and dd[v2] < 0:
                        dd[v2] = dd[u] + 1
                        nx.append(v2)
            fr = nx
        diam = max(diam, max(dd))
    G.diameter = diam

    # ---- exact cyclotomic cross-check of the Gram matrix --------------------
    e = G.exponent
    chi_cyc = []
    for a in range(JMAX + 1):
        row = []
        for ci in range(nc):
            m, n = G.class_lambda[ci]
            step = m * (e // n)
            z = Cyc(e)
            for k in range(a + 1):
                z.c[((a - 2 * k) * step) % e] += 1
            row.append(z)
        chi_cyc.append(row)
    Gram_exact = [[None] * (JMAX + 1) for _ in range(JMAX + 1)]
    gram_ok = True
    for a in range(JMAX + 1):
        for b in range(a, JMAX + 1):
            acc = Cyc(e)
            for ci in range(nc):
                acc = acc + (chi_cyc[a][ci] * chi_cyc[b][ci].conj()).scale(G.class_sizes[ci])
            fr = acc.to_rational()
            if fr is None:
                gram_ok = False
                Gram_exact[a][b] = Gram_exact[b][a] = None
                continue
            val = fr / N
            if val.denominator != 1:
                gram_ok = False
            Gram_exact[a][b] = Gram_exact[b][a] = int(val)
    G.Gram_exact = Gram_exact
    G.gram_rational_ok = gram_ok
    MMt = [[sum(Mult[a][s] * Mult[b][s] for s in range(nc)) for b in range(JMAX + 1)]
           for a in range(JMAX + 1)]
    G.Gram_from_Mult = MMt
    G.gram_matches = (Gram_exact == MMt)

    # ---- explicit irreducible representation matrices -----------------------
    U = {}
    for a in range(0, max(JMAX, MMAX) + 1):
        U[a] = np.array([sym_power_unitary(E[i], a) for i in range(N)])
    G.max_sym_unitarity_err = float(max(
        np.max(np.abs(U[a][i].conj().T @ U[a][i] - np.eye(a + 1)))
        for a in range(0, JMAX + 1) for i in range(0, N, max(1, N // 8))))
    # homomorphism spot check
    hom = 0.0
    for a in (2, 3, 5):
        if a > JMAX:
            continue
        for i in range(0, N, max(1, N // 5)):
            for j in range(0, N, max(1, N // 5)):
                hom = max(hom, float(np.max(np.abs(
                    U[a][i] @ U[a][j] - U[a][mul[i, j]]))))
    G.max_sym_hom_err = hom

    irrep_mats = {}
    irrep_source = {}
    for s in range(nc):
        ds = G.dims[s]
        if ds == 1:
            mats = np.array([[[complex(G.chartab[s][cls_of[i]])]] for i in range(N)])
            irrep_mats[s] = mats
            irrep_source[s] = "1-dim from character"
            continue
        found = None
        for a in range(JMAX + 1):
            if Mult[a][s] == 1:
                found = a
                break
        if found is None:
            raise RuntimeError("%s: no multiplicity-one host for sigma %d" % (name, s))
        cs = np.array([complex(mp.conj(G.chartab[s][cls_of[i]])) for i in range(N)])
        P = np.zeros((found + 1, found + 1), dtype=complex)
        for i in range(N):
            P += cs[i] * U[found][i]
        P *= ds / N
        uu, sv, vh = np.linalg.svd(P)
        rk = int(np.sum(sv > SVD_TOL))
        assert rk == ds, (name, s, rk, ds, sv)
        B = uu[:, :ds]
        mats = np.array([B.conj().T @ U[found][i] @ B for i in range(N)])
        irrep_mats[s] = mats
        irrep_source[s] = "isotypic block of V_%d" % found
    G.irrep_source = irrep_source

    # verify explicit irreps
    cerr = 0.0
    uerr = 0.0
    herr = 0.0
    for s in range(nc):
        ms = irrep_mats[s]
        for i in range(N):
            cerr = max(cerr, abs(np.trace(ms[i]) - complex(G.chartab[s][cls_of[i]])))
            uerr = max(uerr, float(np.max(np.abs(
                ms[i].conj().T @ ms[i] - np.eye(G.dims[s])))))
        step = max(1, N // 6)
        for i in range(0, N, step):
            for j in range(0, N, step):
                herr = max(herr, float(np.max(np.abs(ms[i] @ ms[j] - ms[mul[i, j]]))))
    G.irrep_char_err = cerr
    G.irrep_unitary_err = uerr
    G.irrep_hom_err = herr
    G.irrep_mats = irrep_mats

    # ---- T1 -----------------------------------------------------------------
    T1 = []
    for s in range(nc):
        occ = [(a, Mult[a][s]) for a in range(JMAX + 1) if Mult[a][s] != 0]
        least = occ[0][0] if occ else None
        T1.append({"sigma": s, "dim": G.dims[s], "d": dist[s],
                   "least_a": least,
                   "occurrences": [{"a": a, "mult": m_} for a, m_ in occ]})
    G.T1 = T1

    # ---- T2 -----------------------------------------------------------------
    # the parity test is run for every group; the spec only ASKS for it when
    # -I is present, so 'applicable' records that, but the measurement is
    # reported either way.
    viol = []
    for s in range(nc):
        for a in range(JMAX + 1):
            if Mult[a][s] != 0 and (a - dist[s]) % 2 != 0:
                viol.append({"sigma": s, "dim": G.dims[s], "d": dist[s],
                             "a": a, "mult": Mult[a][s]})
    G.T2 = {"minus_I_in_group": G.has_minus_I,
            "rule_applicable": G.has_minus_I,
            "parity_rule_holds": (not viol),
            "n_violations": len(viol),
            "violations": viol}

    # ---- twists tau: all irreducibles, plus S(rho) from T4 ------------------
    taus = []          # (label, character-on-classes, matrices)
    for s in range(nc):
        taus.append(("sigma%d" % s, [G.chartab[s][i] for i in range(nc)],
                     irrep_mats[s], G.dims[s]))

    # T4: 2-dim irreps with det = 1
    T4 = []
    for s in range(nc):
        if G.dims[s] != 2:
            continue
        lam2 = []
        for i in range(nc):
            v = (G.chartab[s][i] ** 2 - G.chartab[s][G.class_sq[i]]) / 2
            lam2.append(v)
        det1 = all(abs(v - 1) < mp.mpf(10) ** (-20) for v in lam2)
        if not det1:
            continue
        Schar = [(G.chartab[s][i] ** 2 + G.chartab[s][G.class_sq[i]]) / 2
                 for i in range(nc)]
        cons = []
        for s2 in range(nc):
            m_ = note_int(ip(Schar, G.chartab[s2]), "%s S(rho%d) mult %d" % (name, s, s2))
            if m_:
                cons.append({"sigma": s2, "dim": G.dims[s2], "d": dist[s2], "mult": m_})
        assert sum(c["dim"] * c["mult"] for c in cons) == 3, (name, s, cons)
        Smats = np.array([sym_power_unitary(irrep_mats[s][i], 2) for i in range(N)])
        taus.append(("S(rho%d)" % s, Schar, Smats, 3))
        T4.append({"rho_sigma": s, "constituents": cons,
                   "lambda2_is_one": det1})
    G.T4_raw = T4

    # ---- T3 / T5 / T6: mu_tau(m) -------------------------------------------
    mu_tables = {}
    for (lab, ch, mats, dtau) in taus:
        rows2 = []
        for m_ in range(MMIN, MMAX + 1):
            chiE = [(m_ - 1) * chi[m_, i] + (m_ + 1) * chi[m_ - 2, i] for i in range(nc)]
            # convention (a): sum chi_E * chi_tau  (no conjugate)
            sa = mp.mpc(0)
            sb = mp.mpc(0)
            for i in range(nc):
                sa += G.class_sizes[i] * chiE[i] * ch[i]
                sb += G.class_sizes[i] * chiE[i] * mp.conj(ch[i])
            va = note_int(sa / N, "%s mu %s m=%d convA" % (name, lab, m_))
            vb = note_int(sb / N, "%s mu %s m=%d convB" % (name, lab, m_))
            # method (ii): explicit averaging projector, SVD rank
            # E_m(g) = Sym^m(g) (x) I_{m-1}  +  Sym^{m-2}(g) (x) I_{m+1}
            # rank(P_{E_m (x) tau}) = (m-1) rank(P_{V_m (x) tau})
            #                        + (m+1) rank(P_{V_{m-2} (x) tau})
            rk = 0
            svinfo = []
            for (a, rep) in ((m_, m_ - 1), (m_ - 2, m_ + 1)):
                dimP = (a + 1) * dtau
                P = np.zeros((dimP, dimP), dtype=complex)
                for i in range(N):
                    P += np.kron(U[a][i], mats[i])
                P /= N
                sv = np.linalg.svd(P, compute_uv=False)
                r_ = int(np.sum(sv > SVD_TOL))
                rk += rep * r_
                above = float(sv[r_ - 1]) if r_ > 0 else None
                below = float(sv[r_]) if r_ < len(sv) else None
                svinfo.append({"a": a, "copies": rep, "rank": r_,
                               "smallest_sv_above_tol": above,
                               "largest_sv_below_tol": below,
                               "projector_trace_real": float(np.real(np.trace(P))),
                               "projector_trace_imag": float(np.imag(np.trace(P)))})
            rows2.append({"m": m_, "char_tau": va, "char_tau_conj": vb,
                          "svd_rank": rk, "blocks": svinfo,
                          "agree": (va == vb == rk)})
        mu_tables[lab] = {"dim_tau": dtau, "rows": rows2}
    G.mu = mu_tables

    def first_nonzero_m(lab):
        for r_ in mu_tables[lab]["rows"]:
            if r_["char_tau"] != 0:
                return r_["m"]
        return None

    for rec in T4:
        lab = "S(rho%d)" % rec["rho_sigma"]
        qq = first_nonzero_m(lab)
        rec["q"] = qq
        rec["q_squared"] = None if qq is None else qq * qq
        # T7: least k with <chi_tau, chi_k> != 0
        ch = dict((l_, c_) for (l_, c_, _m, _d) in taus)[lab]
        leastk = None
        for k in range(JMAX + 1):
            v = note_int(ip([chi[k, i] for i in range(nc)], ch),
                         "%s T7 %s k=%d" % (name, lab, k))
            if v != 0:
                leastk = k
                break
        rec["T7_least_k"] = leastk
        rec["T7_k_k_plus_2"] = None if leastk is None else leastk * (leastk + 2)
    G.T4 = T4

    # T5
    q_triv = first_nonzero_m("sigma0")
    leastk_triv = None
    for k in range(JMAX + 1):
        v = note_int(ip([chi[k, i] for i in range(nc)], G.chartab[0]),
                     "%s T7 trivial k=%d" % (name, k))
        if v != 0:
            leastk_triv = k
            break
    G.T5 = {"tau": "trivial", "q": q_triv,
            "q_squared": None if q_triv is None else q_triv * q_triv,
            "T7_least_k": leastk_triv,
            "T7_k_k_plus_2": None if leastk_triv is None else leastk_triv * (leastk_triv + 2)}

    # T6
    T6 = []
    for s in range(nc):
        lab = "sigma%d" % s
        e_sigma = first_nonzero_m(lab)
        T6.append({"sigma": s, "dim": G.dims[s], "d": dist[s],
                   "least_sym_level": T1[s]["least_a"], "e": e_sigma})
    G.T6 = T6

    G.elapsed = time.time() - t0
    G.N = N
    G.E = E
    G.mul = mul
    G.U = U
    return G


# ----------------------------------------------------------------------------
# T9: full explicit projector on E_m (x) V_tau, no block reduction
# ----------------------------------------------------------------------------
def t9_full_projector(G, m_, s, log):
    N = G.N
    U = G.U
    tau = G.irrep_mats[s]
    dtau = G.dims[s]
    a1, r1 = m_, m_ - 1
    a2, r2 = m_ - 2, m_ + 1
    dimE = (a1 + 1) * r1 + (a2 + 1) * r2
    dim = dimE * dtau
    P = np.zeros((dim, dim), dtype=complex)
    for i in range(N):
        Eg = np.zeros((dimE, dimE), dtype=complex)
        b1 = (a1 + 1) * r1
        Eg[:b1, :b1] = np.kron(U[a1][i], np.eye(r1))
        Eg[b1:, b1:] = np.kron(U[a2][i], np.eye(r2))
        P += np.kron(Eg, tau[i])
    P /= N
    sv = np.linalg.svd(P, compute_uv=False)
    rk = int(np.sum(sv > SVD_TOL))
    chiE = [(m_ - 1) * G.chi_sym[m_, i] + (m_ + 1) * G.chi_sym[m_ - 2, i]
            for i in range(G.nc)]
    sa = mp.mpc(0)
    for i in range(G.nc):
        sa += G.class_sizes[i] * chiE[i] * G.chartab[s][i]
    csum = note_int(sa / N, "T9 %s" % G.name)
    return {"group": G.name, "m": m_, "tau_sigma": s, "dim_tau": dtau,
            "dim_E_m": dimE, "matrix_dim": dim,
            "svd_tolerance": SVD_TOL,
            "svd_rank": rk, "character_sum": csum, "agree": rk == csum,
            "projector_idempotency_err": float(np.max(np.abs(P @ P - P))),
            "projector_trace_real": float(np.real(np.trace(P))),
            "projector_trace_imag": float(np.imag(np.trace(P))),
            "singular_values": [float(x) for x in sv]}


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------
def main():
    log = []
    out = {"spec": "SPEC SHEET A", "conventions": {}, "groups": {}}
    out["conventions"] = {
        "quaternion_map": "a+bi+cj+dk -> [[a+bi, c+di], [-c+di, a-bi]]",
        "sym_power": ("Sym^a(g) acts on monomials x^(a-k) y^k by "
                      "x -> alpha x + gamma y, y -> beta x + delta y for "
                      "g=[[alpha,beta],[gamma,delta]]; matrices are then "
                      "conjugated to the orthonormal basis sqrt(C(a,k)) x^(a-k) y^k"),
        "inner_product": "<f,g> = (1/|Gamma|) sum_gamma f(gamma) conj(g(gamma))",
        "E_m": "chi_(E_m)(g) = (m-1) chi_m(g) + (m+1) chi_(m-2)(g)",
        "mu_convention_A": "(1/|Gamma|) sum chi_(E_m) * chi_tau   (no conjugate)",
        "mu_convention_B": "(1/|Gamma|) sum chi_(E_m) * conj(chi_tau)",
        "svd_tolerance": SVD_TOL,
        "mpmath_dps": mp.mp.dps,
        "j_range": [0, JMAX],
        "m_range": [MMIN, MMAX],
    }

    Gs = {}
    for (name, gens, order) in build_groups():
        G = analyse_group(name, gens, order, log)
        Gs[name] = G
        print("%-6s |G|=%3d classes=%2d dims=%-22s diam=%d  (%.1fs)"
              % (name, G.order, G.nc, str(G.dims), G.diameter, G.elapsed))

    # T9b: systematic sweep in which the FULL unreduced projector on
    # E_m (x) V_tau is built for every tau and every m, for whole groups.
    # This tests the block reduction used by method (ii) in the main table.
    t9b = []
    for gname in ("BD_2", "BD_3", "2T"):
        G = Gs[gname]
        agree = 0
        total = 0
        worst_kept = 1e9
        worst_dropped = 0.0
        for lab, tab in G.mu.items():
            s = None
            if lab.startswith("sigma"):
                s = int(lab[5:])
            for row in tab["rows"]:
                m_ = row["m"]
                if s is not None:
                    r = t9_full_projector(G, m_, s, log)
                    rk = r["svd_rank"]
                    sv = r["singular_values"]
                else:
                    continue
                total += 1
                if rk == row["char_tau"] == row["svd_rank"]:
                    agree += 1
                if rk > 0:
                    worst_kept = min(worst_kept, sv[rk - 1])
                if rk < len(sv):
                    worst_dropped = max(worst_dropped, sv[rk])
        t9b.append({"group": gname, "cases": total, "agreements": agree,
                    "all_agree": agree == total,
                    "smallest_singular_value_kept": worst_kept,
                    "largest_singular_value_dropped": worst_dropped})
        print("T9b %-5s full unreduced projector, every irreducible tau and "
              "every m: %d/%d agree" % (gname, agree, total))
    out["T9b_full_unreduced_sweep"] = t9b

    # T9 instrumentation on nonabelian groups
    t9 = []
    t9.append(t9_full_projector(Gs["BD_3"], 4, 0, log))
    t9.append(t9_full_projector(Gs["2T"], 4, 0, log))
    t9.append(t9_full_projector(Gs["2T"], 6, 6, log))     # tau = the 3-dim irrep
    t9.append(t9_full_projector(Gs["2I"], 6, 0, log))     # negative control, mu = 0
    t9.append(t9_full_projector(Gs["2O"], 11, 7, log))    # tau = the 4-dim irrep
    t9.append(t9_full_projector(Gs["2I"], 11, 8, log))    # tau = the 6-dim irrep
    t9.append(t9_full_projector(Gs["2I"], 12, 7, log))    # tau = the 5-dim irrep
    out["T9"] = t9
    for r in t9:
        print("T9 %-5s m=%-3d tau=sigma%d(dim %d)  matrix %dx%d  "
              "rank=%d charsum=%d agree=%s"
              % (r["group"], r["m"], r["tau_sigma"], r["dim_tau"],
                 r["matrix_dim"], r["matrix_dim"], r["svd_rank"],
                 r["character_sum"], r["agree"]))

    # serialize
    for name, G in Gs.items():
        blk = {}
        blk["order_expected"] = G.order_expected
        blk["order_verified"] = G.order
        blk["order_matches"] = G.order_ok
        blk["fallback_generator_list_used"] = G.fallback_used
        blk["closure_equals_explicit_spec_list"] = G.explicit_list_agrees
        blk["exponent"] = G.exponent
        blk["minus_I_in_group"] = G.has_minus_I
        blk["checks"] = {
            "max_unitarity_error": G.max_unitarity_err,
            "max_det_minus_one": G.max_det_err,
            "min_pairwise_element_distance": G.min_sep,
            "element_key_injective": G.key_injective,
            "closed_under_multiplication": G.closed,
            "row_orthogonality_max_error": G.orth_row_err,
            "col_orthogonality_max_error": G.orth_col_err,
            "sum_dim_squared": G.sum_dim_sq,
            "sum_dim_squared_equals_order": G.sum_dim_sq_ok,
            "adjacency_symmetric": G.A_symmetric,
            "adjacency_row_weight_ok": G.A_rowsum_ok,
            "branching_dimension_check": G.dim_check,
            "graph_connected": G.connected,
            "sym_power_max_unitarity_error": G.max_sym_unitarity_err,
            "sym_power_max_homomorphism_error": G.max_sym_hom_err,
            "irrep_matrix_character_error": G.irrep_char_err,
            "irrep_matrix_unitarity_error": G.irrep_unitary_err,
            "irrep_matrix_homomorphism_error": G.irrep_hom_err,
            "exact_cyclotomic_gram_rational": G.gram_rational_ok,
            "exact_cyclotomic_gram_matches_MMt": G.gram_matches,
            "eig_randomisation_attempts": G.eig_attempts,
            "joint_eigenvector_consistency_error": G.omega_consistency,
        }
        blk["classes"] = [
            {"index": i, "size": G.class_sizes[i], "element_order": G.class_orders[i],
             "trace": fmt(mp.mpf(G.class_traces[i]), 20),
             "lambda": "exp(2*pi*i*%d/%d)" % G.class_lambda[i],
             "class_of_square": G.class_sq[i],
             "representative": [[str(G.E[G.class_reps[i]][r, c])
                                 for c in range(2)] for r in range(2)]}
            for i in range(G.nc)]
        blk["class_multiplication_matrices"] = G.class_matrices
        blk["irreducible_characters"] = [
            {"sigma": s, "dim": G.dims[s], "d": G.dist[s],
             "source_of_explicit_matrices": G.irrep_source[s],
             "values": [{"class": i,
                         "re": fmt(mp.re(G.chartab[s][i]), 30),
                         "im": fmt(mp.im(G.chartab[s][i]), 30),
                         "exact_rational": as_rational_or_none(G.chartab[s][i])}
                        for i in range(G.nc)]}
            for s in range(G.nc)]
        blk["sym_power_characters"] = [
            {"a": a, "values": [{"class": i, "re": fmt(mp.re(G.chi_sym[a, i]), 30),
                                 "exact_rational": as_rational_or_none(G.chi_sym[a, i])}
                                for i in range(G.nc)]}
            for a in range(JMAX + 1)]
        blk["adjacency_A"] = G.A
        blk["distance_vector"] = G.dist
        blk["graph_diameter"] = G.diameter
        blk["branching_multiplicities_V_a"] = G.Mult
        blk["gram_exact_cyclotomic"] = G.Gram_exact
        blk["gram_from_multiplicities"] = G.Gram_from_Mult
        blk["T1"] = G.T1
        blk["T2"] = G.T2
        blk["T3_mu_tables"] = G.mu
        blk["T4"] = G.T4
        blk["T5"] = G.T5
        blk["T6"] = G.T6
        if name == "2I":
            t8 = []
            for a in range(0, 9):
                cons = [{"sigma": s, "dim": G.dims[s], "d": G.dist[s],
                         "mult": G.Mult[a][s]}
                        for s in range(G.nc) if G.Mult[a][s]]
                t8.append({"a": a, "dim_V_a": a + 1, "constituents": cons,
                           "dimension_check": sum(c["dim"] * c["mult"] for c in cons)})
            blk["T8_branching_2I"] = t8
        out["groups"][name] = blk

    out["audit"] = {
        "max_integer_rounding_deviation": AUDIT["max_integer_rounding_deviation"],
        "max_integer_rounding_where": AUDIT["max_integer_rounding_where"],
        "min_pairwise_element_distance_over_all_groups":
            AUDIT["min_pairwise_element_distance"],
        "element_lookup_dict_misses": AUDIT["element_lookup_dict_misses"],
        "log": log,
    }

    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "m8_1_1_coexact.json"), "w") as f:
        json.dump(out, f, indent=1, sort_keys=False)
    print("\nwrote m8_1_1_coexact.json")
    print("audit:", json.dumps(out["audit"], indent=1)[:800])
    return out, Gs


if __name__ == "__main__":
    main()
