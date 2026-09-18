"""Float D^j(g) (same construction as the exact one), vectorized over many g."""
import numpy as np
from math import comb, factorial


def D_np(j, g11, g12, g21, g22):
    """g entries: arrays of shape (P,). Returns (P, 2j+1, 2j+1) orthonormal CS basis."""
    g11, g12, g21, g22 = map(np.atleast_1d, (g11, g12, g21, g22))
    P = g11.shape[0]
    n = 2 * j
    D = np.zeros((P, n + 1, n + 1), complex)
    pw = {}
    def pow_(a, k, key):
        if (key, k) not in pw:
            pw[(key, k)] = a ** k
        return pw[(key, k)]
    for m in range(-j, j + 1):
        for k1 in range(j + m + 1):
            c1 = comb(j + m, k1) * pow_(g11, k1, 0) * pow_(g21, j + m - k1, 1)
            for k2 in range(j - m + 1):
                c2 = comb(j - m, k2) * pow_(g12, k2, 2) * pow_(g22, j - m - k2, 3)
                mp_ = k1 + k2 - j
                D[:, mp_ + j, m + j] += c1 * c2
    w = np.array([factorial(j + m) * factorial(j - m) / factorial(2 * j) for m in range(-j, j + 1)])
    D *= np.sqrt(w)[None, :, None] / np.sqrt(w)[None, None, :]
    return D


def D_of(j, g):
    g = np.asarray(g)
    return D_np(j, g[0, 0], g[0, 1], g[1, 0], g[1, 1])[0]


def su2_axis_angle(n, theta):
    n = np.asarray(n, float)
    n = n / np.linalg.norm(n)
    sx = np.array([[0, 1], [1, 0]], complex)
    sy = np.array([[0, -1j], [1j, 0]], complex)
    sz = np.array([[1, 0], [0, -1]], complex)
    H = n[0] * sx + n[1] * sy + n[2] * sz
    return np.cos(theta / 2) * np.eye(2) - 1j * np.sin(theta / 2) * H
