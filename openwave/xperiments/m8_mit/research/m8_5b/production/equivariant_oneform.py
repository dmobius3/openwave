"""Deck-equivariant one-form assembly by orbit transport with the fiber map.

THE LAW WAS MEASURED, NOT ADOPTED BY ANALOGY.

The scalar repair transports a row by permuting neighbour indices and reusing the
weights unchanged.  A one-form row additionally carries a FIBER index, and the
placement of the fiber map is not obvious.  The tempting move is to assume the
sandwich `B -> A B A^-1`, which looks natural and is what an operator on one-forms
obeys.  It is WRONG for the object actually being built.

`D_a[i, j]` is the RBF-FD weight for the directional derivative along `e_a` at node
`i`, acting on a SCALAR sample at node `j`.  One frame index, one scalar index.
Only the codomain index rotates.  Evaluating an actual one-form under the frozen
section 2 pullback gives, for `A_gamma = Ad(v_gamma^-1)`:

    d_{gamma i, gamma j} = A_gamma d_{ij}          ONE-SIDED, a vector law

and the sandwich is the DERIVED consequence one level up, on curl and Delta_H,
which carry two frame indices:

    B_{gamma i, gamma j} = A_gamma B_{ij} A_gamma^-1

Measured on L(7;1,2), 60 seeds, 420 nodes, k = 110, by recomputing the image rows
INDEPENDENTLY on the common transported stencil (weight scale 3.374e+00):

    A . w          [derived]              3.58e-10   ADMISSIBLE
    A^-1 . w       [inverse/transpose]    6.58e+00
    w              [omitted]              6.60e+00
    Ad(u^-1) . w   [wrong pair entry]     6.13e+00

The three rejected placements fail at full amplitude, not marginally.  The
wrong-pair-entry row is only meaningful because L(7;1,2) is genuinely two-sided;
a `v = 1` configuration cannot see it, which is the same trap section 8 records
for the manufactured pullback gate.

The sandwich is then CONFIRMED rather than asserted: with these rows in place the
assembled Hodge operator commutes with `T_gamma = kron(Ad(v^-1), P_gamma)` at
3.90e-16.

WHAT IS SHARED WITH THE SCALAR MODULE AND WHAT IS NOT.  Stencil SELECTION is
shared, through `orbit_stencils`, because the two operators must sit on identical
stencils and selection is geometry.  The transformation LAWS are not shared and
must not be: the scalar transport carries no fiber map at all.
"""

import numpy as np

from equivariant_stencils import orbit_stencils
from route_a_tuned import monomials, eval_monos                       # noqa: E402
from route_a_twosided import qmul, qconj, Ad                          # noqa: E402
from route_a_hodge import hodge_matrix                                # noqa: E402

__all__ = ["dirderiv_row", "fiber_map", "transport_matrix",
           "build_D_equivariant", "build_operators_equivariant",
           "ROW_TRANSPORT_MODES"]

IMAG = [np.array([0.0, 1, 0, 0]), np.array([0.0, 0, 1, 0]), np.array([0.0, 0, 0, 1])]

# The mutation vocabulary, named so a harness cannot silently run the right one.
ROW_TRANSPORT_MODES = {
    "correct": "A . w, the measured law, A = Ad(v^-1)",
    "inverse": "A^-1 . w, the inverse/transpose variant",
    "omitted": "w, node indices transported but the fiber map dropped",
    "wrong_entry": "Ad(u^-1) . w, the frozen action's other pair entry",
}


def fiber_map(pair, mode="correct"):
    """The 3x3 coefficient action for one deck element, under one transport mode."""
    u, v = pair
    if mode == "correct":
        return Ad(qconj(v))
    if mode == "inverse":
        return Ad(qconj(v)).T
    if mode == "omitted":
        return np.eye(3)
    if mode == "wrong_entry":
        return Ad(qconj(u))
    raise ValueError(f"unknown transport mode {mode!r}")


def dirderiv_row(center, P, direction, m=7, p=4):
    """One RBF-FD directional-derivative row, transcribed from `build_dirderiv`.

    Transcribed rather than imported because `build_dirderiv` selects its own
    stencils internally, which is precisely the layer being replaced.  The
    weight FORMULA is unchanged, and that is the part seam 3 certified.
    """
    monos = monomials(p)
    npoly = len(monos)
    k = len(P)
    d = np.linalg.norm(P[:, None, :] - P[None, :, :], axis=-1)
    Mx = np.zeros((k + npoly, k + npoly))
    Mx[:k, :k] = d ** m
    Q = eval_monos(P, monos)
    Mx[:k, k:] = Q
    Mx[k:, :k] = Q.T
    diff = center - P
    u0 = np.linalg.norm(diff, axis=1)
    rhs = np.zeros(k + npoly)
    with np.errstate(divide="ignore", invalid="ignore"):
        coef = m * np.where(u0 > 0, u0 ** (m - 2), 0.0)
    rhs[:k] = coef * (diff @ direction)
    for jj, a in enumerate(monos):
        g = 0.0
        for c, e in enumerate(a):
            if e:
                t = e * center[c] ** (e - 1)
                for cc, ee in enumerate(a):
                    if cc != c and ee:
                        t *= center[cc] ** ee
                g += t * direction[c]
        rhs[k + jj] = g
    return np.linalg.lstsq(Mx, rhs, rcond=None)[0][:k]


def build_D_equivariant(X, oid, gid, pairs, k=110, m=7, p=4, mode="correct"):
    """The three directional derivatives, assembled by orbit transport.

    `mode` selects the fiber placement and exists so the law can be MUTATED.  Only
    `"correct"` is admissible; the others are the required negative controls and
    must never be reachable from production code by default.
    """
    node_of, mult, plan = orbit_stencils(X, oid, gid, pairs, k)
    N, G = len(X), len(pairs)
    D = [np.zeros((N, N)) for _ in range(3)]

    for orbit, (rep, idx, moved_by_g) in plan.items():
        # the representative's three rows, in the left-invariant frame AT rep
        w = np.array([dirderiv_row(X[rep], X[idx], qmul(X[rep], IMAG[a]), m, p)
                      for a in range(3)])                   # (3, k)
        for g in range(G):
            tgt = node_of[(orbit, g)]
            wt = fiber_map(pairs[g], mode) @ w              # the measured law
            for a in range(3):
                D[a][tgt, moved_by_g[g]] = wt[a]
    return D


def build_operators_equivariant(X, oid, gid, pairs, k=110, m=7, p=4, mode="correct"):
    """Drop-in replacement for `route_a_twosided.operators` on a Gamma-closed cloud.

    Returns `(Lpos, H)` with the same meaning and the same sign convention.  The
    scalar and one-form layers sit on the SAME transported stencils.
    """
    from equivariant_stencils import build_L_equivariant
    Lpos = build_L_equivariant(X, oid, gid, pairs, k, m, p)
    D = build_D_equivariant(X, oid, gid, pairs, k, m, p, mode)
    return Lpos, hodge_matrix(Lpos, D, len(X))


def transport_matrix(X, oid, gid, pairs, g, mult=None, node_of=None, mode="correct"):
    """`T_gamma = kron(A_gamma, P_gamma)` in the component-major 3N layout.

    The component index is OUTER and the node index INNER, matching
    `hodge_matrix`'s `a*N + i` blocking.  Getting that backwards produces a
    plausible-looking operator that measures nothing.
    """
    from equivariant_stencils import group_multiplication_table
    N = len(X)
    if mult is None:
        mult = group_multiplication_table(pairs)
    if node_of is None:
        node_of = {(int(oid[i]), int(gid[i])): i for i in range(N)}
    perm = np.array([node_of[(int(oid[i]), mult[g][int(gid[i])])] for i in range(N)])
    P = np.zeros((N, N))
    P[perm, np.arange(N)] = 1.0
    return np.kron(fiber_map(pairs[g], mode), P)
