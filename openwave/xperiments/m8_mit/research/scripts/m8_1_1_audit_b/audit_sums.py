#!/usr/bin/env python3
"""The two sums S and D, computed by ELEMENT summation (primary) plus three
mutually independent cross-routes:
  * closed form for C_n / BD_n from the exact cotangent-sum T(m,j);
  * a purely rational McKay/Molien linear solve  (2I - A) u = e_triv - delta/|G|,
    <u, delta> = 0, whose alpha-entry is D(alpha)   -- no field arithmetic at all;
  * high precision mpmath evaluation from the literal angles phi_g.
"""

from fractions import Fraction as Fr
from math import gcd
import mpmath as mp

from audit_core import lin_solve, null_space, divisors, Tsum, prim_T
from audit_groups import inner, cf_mul, cf_add, cf_scal, cf_sub, cyclic_subgroups

mp.mp.dps = 80


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
