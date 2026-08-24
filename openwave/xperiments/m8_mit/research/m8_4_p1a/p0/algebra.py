"""Quaternion algebra, SU(2) matrices, and symmetric-power representations.

Self-contained: no imports from outside this room.
"""

import numpy as np
from math import comb


def qmul(a, b):
    """Quaternion product a*b, Hamilton convention (w, x, y, z)."""
    w1, x1, y1, z1 = a
    w2, x2, y2, z2 = b
    return np.array([
        w1*w2 - x1*x2 - y1*y2 - z1*z2,
        w1*x2 + x1*w2 + y1*z2 - z1*y2,
        w1*y2 - x1*z2 + y1*w2 + z1*x2,
        w1*z2 + x1*y2 - y1*x2 + z1*w2,
    ])


def qconj(q):
    """Quaternion conjugate."""
    return np.array([q[0], -q[1], -q[2], -q[3]])


def quat_to_su2(q):
    """Unit quaternion to 2x2 SU(2) matrix."""
    w, x, y, z = q
    return np.array([[w + 1j*x, y + 1j*z],
                     [-y + 1j*z, w - 1j*x]])


def su2_character(q, n):
    """Character of V_n at unit quaternion q: sin((n+1)α)/sin(α), α = arccos(w)."""
    w = np.clip(float(q[0]), -1.0, 1.0)
    a = np.arccos(w)
    s = np.sin(a)
    if abs(s) < 1e-12:
        return float(n + 1) * (1.0 if abs(a) < 1e-12 else (-1.0)**n)
    return np.sin((n + 1) * a) / s


def sym_power(U, n):
    """(n+1)x(n+1) matrix of the n-th symmetric power of a 2x2 matrix U.

    Convention: (g·f)(z) = f(g⁻¹z), monomial basis e_k = z1^k z2^(n-k),
    k = 0, ..., n.  For U in SU(2), U⁻¹ = U†.
    """
    if n == 0:
        return np.array([[1.0 + 0j]])
    d = n + 1
    Ui = np.conj(U).T
    a, b = Ui[0, 0], Ui[0, 1]
    c, dd = Ui[1, 0], Ui[1, 1]
    R = np.zeros((d, d), dtype=complex)
    for k in range(d):
        for i in range(k + 1):
            ci = comb(k, i) * (a**i) * (b**(k - i))
            for ll in range(n - k + 1):
                cl = comb(n - k, ll) * (c**ll) * (dd**(n - k - ll))
                R[i + ll, k] += ci * cl
    return R


def invariant_gram(n):
    """The SU(2)-invariant Gram matrix on the monomial basis of V_n.

    H_{jk} = δ_{jk} / C(n, j).  The representation matrices satisfy
    ρ(g)† H ρ(g) = H for all g in SU(2).
    """
    return np.diag([1.0 / comb(n, k) for k in range(n + 1)])


def Ad(v):
    """Adjoint action Ad(v) : su(2) → su(2) ≅ R³, for unit quaternion v.

    Ad(v)(ξ) = v ξ v⁻¹ where ξ is a pure imaginary quaternion.
    Returns a 3x3 real matrix.
    """
    IMAG = [np.array([0.0, 1, 0, 0]),
            np.array([0.0, 0, 1, 0]),
            np.array([0.0, 0, 0, 1])]
    vc = qconj(v)
    A = np.zeros((3, 3))
    for j in range(3):
        img = qmul(qmul(v, IMAG[j]), vc)
        A[:, j] = img[1:4]
    return A
