"""Shared definitions for the stage-1 V2 recomputation (R = 1 throughout)."""
import json
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
R = 1.0

# Width labels exactly as required in stage1.json, with W/R values.
WIDTHS = [
    ("1/4", 0.25),
    ("1/2", 0.5),
    ("1", 1.0),
    ("7/5", 1.4),
    ("3/2", 1.5),
    ("pi/2", math.pi / 2),
    ("17/10", 1.7),
]
NEIG = 6


def path(*p):
    return os.path.join(HERE, *p)


def save_json(name, obj):
    with open(path(name), "w") as f:
        json.dump(obj, f, indent=1, sort_keys=False)


def load_json(name):
    with open(path(name)) as f:
        return json.load(f)


def richardson(seq):
    """Order from the last three levels (ratio-2 refinement), extrapolation from
    the last two, error estimate |extrap - finest|.  Returns dict (all floats;
    p is None when the differences do not shrink monotonically with the same
    sign, in which case extrap = finest and err = |finest - previous|)."""
    a, b, c = seq[-3], seq[-2], seq[-1]
    d1, d2 = b - a, c - b
    if d1 != 0 and d2 != 0 and d1 * d2 > 0 and abs(d2) < abs(d1):
        p = math.log2(d1 / d2)
        ext = c + d2 / (2.0 ** p - 1.0)
        err = abs(ext - c)
        return {"p": p, "extrap": ext, "err": err}
    return {"p": None, "extrap": c, "err": abs(d2)}


def lagrange_1d(nodes):
    """Return functions (vals, ders) evaluating the Lagrange basis on the
    reference nodes at points x (arrays shape (len(nodes), len(x)))."""
    nodes = np.asarray(nodes, float)
    n = len(nodes)

    def vals(x):
        x = np.asarray(x, float)
        V = np.ones((n, len(x)))
        for i in range(n):
            for j in range(n):
                if j != i:
                    V[i] *= (x - nodes[j]) / (nodes[i] - nodes[j])
        return V

    def ders(x):
        x = np.asarray(x, float)
        D = np.zeros((n, len(x)))
        for i in range(n):
            for m in range(n):
                if m == i:
                    continue
                t = np.ones(len(x)) / (nodes[i] - nodes[m])
                for j in range(n):
                    if j != i and j != m:
                        t = t * (x - nodes[j]) / (nodes[i] - nodes[j])
                D[i] += t
        return D

    return vals, ders


def gll_nodes(q):
    """Gauss-Lobatto-Legendre nodes of degree q on [-1, 1]."""
    if q == 1:
        return np.array([-1.0, 1.0])
    P = np.polynomial.legendre.Legendre.basis(q).deriv()
    inner = np.sort(np.real(P.roots()))
    return np.concatenate([[-1.0], inner, [1.0]])


def graded_s(N, g):
    """Distances from the collapsed fiber on one side: s_k = (pi/2)(k/N)^g."""
    k = np.arange(N + 1)
    return (math.pi / 2) * (k / N) ** g


def element_matrices_s(sa, sb, q, nquad, mu):
    """Element matrices on [sa, sb] (s = distance to the fiber, |cos y| = sin s)
    for degree-q Lagrange on GLL nodes.  Returns (A, B, Mm):
      A = int sin(s) N_i' N_j' ds,  B = int N_i N_j / sin(s) ds,
      Mm = int sin(s) N_i N_j ds.   mu is unused here (kept for clarity)."""
    xi, wq = np.polynomial.legendre.leggauss(nquad)
    vals, ders = lagrange_1d(gll_nodes(q))
    V, D = vals(xi), ders(xi)
    h = sb - sa
    s = sa + (xi + 1) * h / 2
    jac = h / 2
    sn = np.sin(s)
    A = (D * (wq * sn / jac)) @ D.T
    B = (V * (wq / sn * jac)) @ V.T
    Mm = (V * (wq * sn * jac)) @ V.T
    return A, B, Mm
