"""Route (a)'s two-sided scalar harmonic representation: counts AND subspaces.

ROUTE-A ONLY.  Section 2 method disjointness: route (a)'s primary multiplicity
measurement is CHARACTER-FREE, by commutant or intertwiner rank.  Nothing here
computes or consumes a character, and no route (b) module may import it.

REALIZATION.  Frozen in `TWOSIDED_REALIZATION.md` and derived there from
section 2's pullback law, not chosen:

    column-major vec:   M_n(u, v) = rho_n(v)^T  kron  rho_n(u)
    row-major  vec:     M_n(u, v) = rho_n(u)    kron  rho_n(v)^T

the SAME operator.  The transpose on the right factor is load-bearing and is not
an artifact of the vec convention.

WHY THIS MODULE COMPUTES BOTH COUNTS AND BASES FROM ONE OPERATOR.  The
no-transpose form `rho_n(u) kron rho_n(v)` is character-equivalent, so it
reproduces every invariant DIMENSION exactly while fixing a maximally different
SUBSPACE at n >= 3.  Splitting counts and bases across two realizations would let
that divergence live inside a single route undetected.  One operator serves both.

    Count-only safety is not subspace correctness.  Any use of this machinery to
    build bases, projectors or symmetry-adapted subspaces must pass the pointwise
    realization test in `verify_realization`.

THE n = 2 TRAP.  At the first nontrivial level the correct and no-transpose
subspaces COINCIDE.  A regression that stops there passes under the wrong
realization.  `verify_realization` therefore sweeps n >= 3 by default.
"""

import sys

import numpy as np

sys.path.insert(0, "../pilot")
from route_a_nonabelian import quat_to_su2, sym_power   # noqa: E402

__all__ = ["level_operator", "invariant_dim_and_basis", "verify_realization",
           "REALIZATION", "VEC_ORDER"]

REALIZATION = "rho_n(v)^T kron rho_n(u), column-major vec (TWOSIDED_REALIZATION.md)"
VEC_ORDER = "F"


def level_operator(u, v, n, no_transpose=False):
    """M_n(u, v) on the (n+1)^2-dimensional level space, column-major vec.

    `no_transpose=True` selects the WRONG character-equivalent form and exists
    only so the mutation harness can exhibit its signature.  It is never used in
    production.
    """
    U = sym_power(quat_to_su2(np.asarray(u, dtype=float)), n)
    V = sym_power(quat_to_su2(np.asarray(v, dtype=float)), n)
    if no_transpose:
        return np.kron(U, V)
    return np.kron(V.T, U)


def coefficient_operator(u, v, n, no_transpose=False):
    """The action on COEFFICIENT VECTORS: M_n(u, v)^T.

    `level_operator` acts on BASIS VALUES; it is verified as
    `M_n vec(rho_n(x)) = vec(rho_n(u x v))`.  A function is
    `f(x) = c . vec(rho_n(x))`, so

        f(u x v) = c . (M_n vec(rho_n(x))) = (M_n^T c) . vec(rho_n(x))

    and the coefficient action is the TRANSPOSE.  Verified pointwise at 1.1e-15
    against 4.9e+00 for the untransposed form.

    RECORDED NARROWLY: `M_n` is the wrong coefficient action, even though on the
    tuning groups tested, replacing it by `M_n^T` changes neither the invariant
    dimensions nor the resulting invariant subspace (sin theta_max = 0 at
    n = 2..6).  Whether that agreement is a coincidence or a theorem about these
    groups is NOT settled here and is not relied on: the pointwise test is the
    evidence, and the tuning agreement is only a regression.
    """
    return level_operator(u, v, n, no_transpose).T


def invariant_dim_and_basis(pairs, n, tol_rel=1e-8, no_transpose=False):
    """Character-free invariant dimension AND an orthonormal invariant basis.

    Nullity of the stacked `(M_n(u,v)^T - I)` over the supplied effective group,
    from singular values, with the section 6.3 projector-rank tolerance of 1e-8
    relative to the largest singular value.

    Uses the COEFFICIENT action, since invariance is a statement about functions:
    `f(gamma . x) = f(x)` for all x.

    Returns (dimension, basis, conditioning_gap).  The basis spans the fixed
    space of the SAME operator the dimension came from.
    """
    dim = (n + 1) ** 2
    blocks = [coefficient_operator(u, v, n, no_transpose) - np.eye(dim)
              for (u, v) in pairs]
    A = np.vstack(blocks)
    _, s, Vh = np.linalg.svd(A)
    smax = s[0] if s.size else 0.0
    cutoff = max(tol_rel * smax, 1e-12)
    k = int(np.sum(s < cutoff))
    below, above = s[s < cutoff], s[s >= cutoff]

    # Conditioning diagnostic: how cleanly the nullspace separates from the rest.
    # The zero cases are branched EXPLICITLY rather than suppressed, because each
    # has a distinct mathematical meaning that `inf` alone would conflate:
    #   below empty      no nullspace at all; separation is vacuous
    #   above empty      everything is null; separation is vacuous
    #   below.max() == 0 the nullspace is exactly zero to machine precision,
    #                    which is PERFECT separation, genuinely infinite gap
    if below.size == 0 or above.size == 0:
        gap = float("inf")
        gap_state = "vacuous: one side of the cut is empty"
    elif below.max() == 0.0:
        gap = float("inf")
        gap_state = "exact: nullspace is identically zero, separation is perfect"
    else:
        gap = float(above.min() / below.max())
        gap_state = "measured"

    basis = Vh.conj().T[:, dim - k:] if k else np.zeros((dim, 0), dtype=complex)
    return k, basis, {"gap": gap, "state": gap_state}


def verify_realization(rng=None, levels=(3, 4, 5), trials=12, tol=1e-10):
    """The manufactured pointwise test that fixes the realization.

    Evaluates `rho_n` at the transformed point `u x v` directly, with no
    representation-theoretic assumption, and requires the operator to reproduce
    it.  This is the only instrument that can separate the candidate right-factor
    conventions, since all of them share a character.

    Sweeps n >= 3 by default: at n = 2 the correct and no-transpose subspaces
    coincide and the test cannot discriminate.
    """
    rng = rng or np.random.default_rng(20260811)
    worst_correct, worst_wrong = 0.0, np.inf
    for n in levels:
        for _ in range(trials):
            q = [rng.normal(size=4) for _ in range(3)]
            u, v, x = [a / np.linalg.norm(a) for a in q]
            truth = sym_power(quat_to_su2(u) @ quat_to_su2(x) @ quat_to_su2(v), n)
            src = sym_power(quat_to_su2(x), n).reshape(-1, order=VEC_ORDER)
            tgt = truth.reshape(-1, order=VEC_ORDER)
            worst_correct = max(worst_correct,
                                float(np.abs(level_operator(u, v, n) @ src - tgt).max()))
            worst_wrong = min(worst_wrong,
                              float(np.abs(level_operator(u, v, n, True) @ src - tgt).max()))
    return {
        "realization": REALIZATION,
        "levels_swept": list(levels),
        "worst_residual_correct": worst_correct,
        "best_residual_no_transpose": worst_wrong,
        "pass": bool(worst_correct < tol and worst_wrong > tol),
        "note": ("the no-transpose form is character-equivalent; only this "
                 "pointwise evaluation separates them"),
    }


# --- n_max, section 6.2, route-A-owned ---------------------------------------

def derive_n_max(pairs, search_to=40, tol_rel=1e-8):
    """Section 6.2's harmonic ceiling, derived character-free by route (a).

    FROZEN RULE, reading A:

        scan harmonic levels n > 0 in increasing order;
        retain those with invariant dimension > 0;
        n_max = the SECOND retained n.

    "Nonzero" qualifies the LEVEL, not the invariant dimension.  Section 6.2's
    own worked example calls lambda = 48 the first nonzero scalar invariant level
    and lambda = 80 the second, and section 6.3 contrasts "a NONZERO level" with
    "The ZERO level" throughout, so the quotient's zero level is a distinct
    spectral object rather than the first member of this count.

    THE DISCARDED READING is a required regression, not a matter of taste:
    counting n = 0 (whose invariant dimension is always 1) returns n_max = 2
    mechanically on every case where n = 1 is absent, which would destroy the
    case-dependence section 6.2 states the ceiling has, and would shrink the
    required band on every quotient.  `derive_n_max_wrong_reading` exists solely
    so a test can reject it.

    Derived from route (a)'s own two-sided character-free rank machinery.  It is
    never taken from route (b); section 6.2 requires both routes to certify the
    SAME n_max, which is only a check if each derives it independently.
    """
    positive = []
    for n in range(1, search_to + 1):
        k, _, _ = invariant_dim_and_basis(pairs, n, tol_rel=tol_rel)
        if k > 0:
            positive.append(n)
            if len(positive) == 2:
                return {"n_max": positive[1], "first_two_positive_levels": positive,
                        "rule": "second n > 0 with positive scalar invariant dimension",
                        "derived_by": "route (a) two-sided character-free rank"}
    raise ValueError(f"fewer than two positive invariant levels below n={search_to}")


def derive_n_max_wrong_reading(pairs, search_to=40, tol_rel=1e-8):
    """The DISCARDED reading: counts n = 0.  Present only as a mutation target."""
    positive = []
    for n in range(0, search_to + 1):
        k, _, _ = invariant_dim_and_basis(pairs, n, tol_rel=tol_rel)
        if k > 0:
            positive.append(n)
            if len(positive) == 2:
                return positive[1]
    raise ValueError("fewer than two positive invariant levels")
