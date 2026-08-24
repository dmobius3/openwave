"""RBF-FD weight computation on S³.

Self-contained reimplementation of the polyharmonic-spline + polynomial
augmentation scheme used by the production backend.
"""

import numpy as np
from itertools import product as iproduct
from scipy.spatial import cKDTree


def monomials(p):
    """All 4D multiindices with total degree ≤ p, sorted by degree then lex."""
    monos = []
    for deg in range(p + 1):
        for a in range(deg + 1):
            for b in range(deg - a + 1):
                for c in range(deg - a - b + 1):
                    d = deg - a - b - c
                    monos.append((a, b, c, d))
    return monos


def eval_monos(P, monos):
    """Evaluate monomials at an array of points. Returns (len(P), len(monos))."""
    Q = np.zeros((len(P), len(monos)))
    for j, alpha in enumerate(monos):
        col = np.ones(len(P))
        for dim in range(4):
            if alpha[dim] > 0:
                col *= P[:, dim] ** alpha[dim]
        Q[:, j] = col
    return Q


def _flat_laplacian_mono(center, alpha):
    """Δ_{R⁴} of monomial x^alpha evaluated at center."""
    val = 0.0
    for dim in range(4):
        e = alpha[dim]
        if e >= 2:
            coeff = e * (e - 1)
            term = coeff
            for d2 in range(4):
                if d2 == dim:
                    term *= center[d2] ** (e - 2)
                else:
                    term *= center[d2] ** alpha[d2]
            val += term
    return val


def _degree(alpha):
    return sum(alpha)


def surf_lap_monos(center, monos):
    """Δ_{S³} of each monomial at center.

    For homogeneous polynomial of degree d on S³:
        Δ_{S³} p = Δ_{R⁴} p − d(d+2) p
    """
    result = np.zeros(len(monos))
    for j, alpha in enumerate(monos):
        d = _degree(alpha)
        mono_val = 1.0
        for dim in range(4):
            mono_val *= center[dim] ** alpha[dim]
        flat_lap = _flat_laplacian_mono(center, alpha)
        result[j] = flat_lap - d * (d + 2) * mono_val
    return result


def rbf_row(center, P, m=7, p=4):
    """One RBF-FD row for the surface Laplacian on S³.

    Polyharmonic spline |x-x_j|^m with polynomial augmentation up to degree p.
    Returns weights w such that Δ_{S³}f(center) ≈ Σ w_j f(x_j).
    """
    monos = monomials(p)
    npoly = len(monos)
    k = len(P)
    dist = np.linalg.norm(P[:, None, :] - P[None, :, :], axis=-1)
    Mx = np.zeros((k + npoly, k + npoly))
    Mx[:k, :k] = dist ** m
    Q = eval_monos(P, monos)
    Mx[:k, k:] = Q
    Mx[k:, :k] = Q.T

    rhs = np.zeros(k + npoly)
    u0 = np.linalg.norm(P - center, axis=1)
    rhs[:k] = (m * (m + 1) * np.where(u0 > 0, u0 ** (m - 2), 0.0)
               - (m * (m + 4) / 4.0) * u0 ** m)
    rhs[k:] = surf_lap_monos(center, monos)

    return np.linalg.lstsq(Mx, rhs, rcond=None)[0][:k]
