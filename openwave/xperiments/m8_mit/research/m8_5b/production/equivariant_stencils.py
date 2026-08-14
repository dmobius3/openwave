"""Deck-equivariant scalar operator assembly by orbit transport.

WHY THIS EXISTS, stated as the demonstrated reason rather than a guess:

    Independent nearest-neighbour selection is not an admissible quotient-stencil
    construction, because distance ties at the stencil boundary can break
    deck-equivariance even on an exactly invariant point cloud.  Orbit-related
    rows are therefore constructed by transport from a common representative
    stencil.

The three seams were tested separately before this was written, and only the
middle one failed:

    cloud covariance   A_gamma x_i = x_pi(i) under the frozen u x v law:  0.000e+00
    stencil-set cov.   pi(N(i)) = N(pi(i)):  FAILS, 11 of 840 on L(7;1,2),
                       ALL of them at the k-th distance tie
    weight covariance  same transported stencil, weights recomputed:  9.36e-11

So the RBF-FD formulas, the distance metric, the stencil size and the seed
ladder are NOT implicated and are left untouched.  Only selection changes.

A deterministic tie-break would be reproducible but still not equivariant, so it
is NOT a substitute: reproducibility is a G8 property, equivariance is the
structural property the quotient construction requires.

SCALAR ONLY.  One-form rows need the frozen coefficient action carried through
the transport as well, which is a separate certification with its own
manufactured pointwise check.  Do not generalize this rule to blocks by copying.
"""

import numpy as np
from scipy.spatial import cKDTree

import sys
sys.path.insert(0, "../pilot")
from route_a_tuned import monomials, eval_monos, surf_lap_monos   # noqa: E402
from route_a_twosided import qmul                                  # noqa: E402

__all__ = ["group_multiplication_table", "rbf_row", "orbit_stencils",
           "build_L_equivariant", "STENCIL_COVARIANCE", "RBF_FD_WEIGHT_COVARIANCE"]


# TWO DIFFERENT PROPERTIES UNDER TWO DIFFERENT NAMES, and the second does NOT
# follow from the first.  Frozen here so a later reader cannot collapse them.
STENCIL_COVARIANCE = {
    "statement": "pi(N(i)) = N(pi(i)) for every node and every deck element",
    "how_established": "EXACT BY CONSTRUCTION under orbit transport",
    "evidence": "identically zero; the transported index set IS the image set",
    "warning": ("this is partly tautological once rows are transported, and on its "
                "own it certifies NOTHING about the weights"),
}
RBF_FD_WEIGHT_COVARIANCE = {
    "statement": ("on a COMMON transported stencil, independently recomputed RBF-FD "
                  "weights are covariant"),
    "how_established": "MEASURED, independently, before the transport was written",
    "evidence": "seam 3 on L(7;1,2): 9.36e-11",
    "warning": ("this is what LICENSES orbit transport.  It does not follow from "
                "STENCIL_COVARIANCE by any argument; seam 3 earned it separately, "
                "and removing seam 3 would leave the transport unjustified"),
}


def group_multiplication_table(pairs, tol=1e-9):
    def same(a, b):
        return abs(a[0] - b[0]).max() < tol and abs(a[1] - b[1]).max() < tol

    def idx_of(p):
        for i, q in enumerate(pairs):
            if same(p, q):
                return i
        raise ValueError("product left the supplied element list")

    return [[idx_of((qmul(pairs[i][0], pairs[j][0]),
                     qmul(pairs[i][1], pairs[j][1])))
             for j in range(len(pairs))] for i in range(len(pairs))]


def rbf_row(center, P, m=7, p=4):
    """One RBF-FD row, transcribed unchanged from the pilot's `build_L`."""
    monos = monomials(p)
    npoly = len(monos)
    k = len(P)
    d = np.linalg.norm(P[:, None, :] - P[None, :, :], axis=-1)
    Mx = np.zeros((k + npoly, k + npoly))
    Mx[:k, :k] = d ** m
    Q = eval_monos(P, monos)
    Mx[:k, k:] = Q
    Mx[k:, :k] = Q.T
    rhs = np.zeros(k + npoly)
    u0 = np.linalg.norm(P - center, axis=1)
    rhs[:k] = (m * (m + 1) * np.where(u0 > 0, u0 ** (m - 2), 0.0)
               - (m * (m + 4) / 4.0) * u0 ** m)
    rhs[k:] = surf_lap_monos(center, monos)
    return np.linalg.lstsq(Mx, rhs, rcond=None)[0][:k]


def orbit_stencils(X, oid, gid, pairs, k=110):
    """Stencil SELECTION only, shared by every operator built on this cloud.

    Returns `(node_of, mult, plan)` where `plan[orbit] = (rep, idx, moved_by_g)`:
    one stencil chosen per ORBIT at the canonical representative (the seed itself,
    `gid == 0`, which is deterministic from the cloud construction rather than
    from an incidental node ordering), plus that stencil transported to every
    other orbit member.

    This is GEOMETRY, not a transformation law, which is why the scalar and
    one-form assemblies share it: they must sit on identical stencils.  The LAWS
    they impose on the weights are different and are NOT shared.
    """
    N = len(X)
    G = len(pairs)
    mult = group_multiplication_table(pairs)
    node_of = {(int(oid[i]), int(gid[i])): i for i in range(N)}
    tree = cKDTree(X)

    plan = {}
    for orbit in sorted({int(o) for o in oid}):
        rep = node_of[(orbit, 0)]
        _, idx = tree.query(X[rep], k=k)
        idx = np.asarray(idx)
        moved_by_g = [np.array([node_of[(int(oid[s]), mult[g][int(gid[s])])]
                                for s in idx]) for g in range(G)]
        plan[orbit] = (rep, idx, moved_by_g)
    return node_of, mult, plan


def build_L_equivariant(X, oid, gid, pairs, k=110, m=7, p=4):
    """Scalar Laplacian assembled so that P_gamma^-1 L P_gamma = L by construction.

    Each orbit member receives the representative's row transported by the frozen
    node permutation: neighbour indices mapped, WEIGHTS IDENTICAL.  A scalar is a
    trivial fiber, so the transport carries no fiber map.  That is the whole
    difference from the one-form case, and it is why the one-form transport is a
    separate certification rather than a copy of this one.

    Equivariance is then a property of the construction, not an emergent property
    that appears once distance ties happen to stop occurring.
    """
    node_of, mult, plan = orbit_stencils(X, oid, gid, pairs, k)
    G = len(pairs)
    L = np.zeros((len(X), len(X)))

    for orbit, (rep, idx, moved_by_g) in plan.items():
        w = rbf_row(X[rep], X[idx], m, p)
        for g in range(G):
            tgt = node_of[(orbit, g)]
            L[tgt, moved_by_g[g]] = w   # weights are invariant; only indices move
    return -L                           # sign convention matches `operators`
