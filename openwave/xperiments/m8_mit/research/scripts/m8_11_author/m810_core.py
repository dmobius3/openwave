"""M8.10 fresh derivation: core SU(2) and 2I machinery, written from scratch (2026-09-12).

Independence: imports numpy, scipy and the standard library only. Nothing from
M8_DYNAMICS/build, the openwave repository, or any C/C2 room.

Conventions: Condon-Shortley; basis |j m>, m = j, j-1, ..., -j; normalized Haar measure
on SU(2); R = 1, so level n = 2j has -Delta eigenvalue n(n+2). A quaternion
q = (a, b, c, d) acts as U(q) = a I - i(b sx + c sy + d sz).
"""
import itertools
from fractions import Fraction
from math import factorial, prod, sqrt

import numpy as np
from scipy.linalg import expm

PHI = (1 + sqrt(5)) / 2


def mvals(j):
    return [j - k for k in range(2 * j + 1)]


def spin_matrices(j):
    ms, n = mvals(j), 2 * j + 1
    jp = np.zeros((n, n), complex)
    for c in range(1, n):  # J+ |j m> = sqrt(j(j+1) - m(m+1)) |j m+1>
        m = ms[c]
        jp[c - 1, c] = sqrt(j * (j + 1) - m * (m + 1))
    jm = jp.conj().T
    return (jp + jm) / 2, (jp - jm) / (2 * 1j), np.diag(ms).astype(complex)


def qmul(p, q):
    a1, b1, c1, d1 = p
    a2, b2, c2, d2 = q
    return np.array([a1 * a2 - b1 * b2 - c1 * c2 - d1 * d2,
                     a1 * b2 + b1 * a2 + c1 * d2 - d1 * c2,
                     a1 * c2 - b1 * d2 + c1 * a2 + d1 * b2,
                     a1 * d2 + b1 * c2 - c1 * b2 + d1 * a2])


def su2(q):
    a, b, c, d = q
    return np.array([[a - 1j * d, -1j * b - c], [-1j * b + c, a + 1j * d]])


def wigner_D(j, q):
    """Spin-j matrix of U(q) = exp(-i theta n.sigma/2), computed as exp(-i theta n.J)."""
    a, v = q[0], np.asarray(q[1:], float)
    s = np.linalg.norm(v)
    if s < 1e-14:
        return np.eye(2 * j + 1, dtype=complex) * (1.0 if a > 0 else (-1.0) ** (2 * j))
    theta, n = 2 * np.arctan2(s, a), v / s
    jx, jy, jz = spin_matrices(j)
    return expm(-1j * theta * (n[0] * jx + n[1] * jy + n[2] * jz))


def conj_matrix(j):
    """W with conj(D^j) = W D^j W^T (integer j)."""
    ms = mvals(j)
    w = np.zeros((2 * j + 1, 2 * j + 1))
    for i, m in enumerate(ms):
        w[i, ms.index(-m)] = (-1) ** (j - m)
    return w


def cg_exact(j1, m1, j2, m2, J, M):
    """<j1 m1; j2 m2 | J M> = s * sqrt(pre), with s and pre exact Fractions (Racah formula)."""
    if m1 + m2 != M or not abs(j1 - j2) <= J <= j1 + j2 or abs(m1) > j1 or abs(m2) > j2 or abs(M) > J:
        return Fraction(0), Fraction(0)
    f = factorial
    pre = Fraction((2 * J + 1) * f(J + j1 - j2) * f(J - j1 + j2) * f(j1 + j2 - J)
                   * f(J + M) * f(J - M) * f(j1 - m1) * f(j1 + m1) * f(j2 - m2) * f(j2 + m2),
                   f(j1 + j2 + J + 1))
    s = Fraction(0)
    for k in range(j1 + j2 - J + 1):
        den = [k, j1 + j2 - J - k, j1 - m1 - k, j2 + m2 - k, J - j2 + m1 + k, J - j1 - m2 + k]
        if min(den) >= 0:
            s += Fraction((-1) ** k, prod(f(x) for x in den))
    return s, pre


def cg(j1, m1, j2, m2, J, M):
    s, pre = cg_exact(j1, m1, j2, m2, J, M)
    return float(s) * sqrt(float(pre))


def coupling(j1, j2, J):
    """C_J : V_j1 (x) V_j2 -> V_J, rows |J M>, columns i1*(2 j2 + 1) + i2."""
    n2 = 2 * j2 + 1
    C = np.zeros((2 * J + 1, (2 * j1 + 1) * n2))
    for a, M in enumerate(mvals(J)):
        for i1, m1 in enumerate(mvals(j1)):
            for i2, m2 in enumerate(mvals(j2)):
                C[a, i1 * n2 + i2] = cg(j1, m1, j2, m2, J, M)
    return C


def _parity(p):
    return sum(1 for i in range(len(p)) for k in range(i + 1, len(p)) if p[i] > p[k]) % 2


def binary_icosahedral():
    els = []
    for i in range(4):
        for s in (1.0, -1.0):
            q = [0.0] * 4
            q[i] = s
            els.append(q)
    els += [list(s) for s in itertools.product((0.5, -0.5), repeat=4)]
    base = (0.0, 0.5, 0.5 / PHI, 0.5 * PHI)
    for p in itertools.permutations(range(4)):
        if _parity(p) == 0:
            for sb, sc, sd in itertools.product((1, -1), repeat=3):
                v = (base[0], sb * base[1], sc * base[2], sd * base[3])
                els.append([v[p[k]] for k in range(4)])
    return np.array(els)


# Rotation-angle classes of I = 2I/{+-1}, keyed by cos(theta) = 2a^2 - 1, with (chi_3', chi_4).
CLASSES = [(1.0, 'e', 3, 4), (-1.0, 'C2', -1, 0), (-0.5, 'C3', 0, 1),
           ((PHI - 1) / 2, 'C5 (2pi/5)', 1 - PHI, -1), (-PHI / 2, 'C5 (4pi/5)', PHI, -1)]
CHI_3 = {'e': 3, 'C2': -1, 'C3': 0, 'C5 (2pi/5)': PHI, 'C5 (4pi/5)': 1 - PHI}   # for the planted control
DIM = {'3p': 3, '4': 4, '3': 3}


def rot_class(q):
    c = 2 * q[0] ** 2 - 1
    hits = [k for k in CLASSES if abs(k[0] - c) < 1e-9]
    assert len(hits) == 1, (q, c)
    return hits[0]


def character(sector, q):
    k = rot_class(q)
    return {'3p': k[2], '4': k[3], '3': CHI_3[k[1]]}[sector]


def isotypic_projector(j, sector, G):
    P = sum(character(sector, q) * wigner_D(j, q) for q in G)
    return P * DIM[sector] / len(G)
