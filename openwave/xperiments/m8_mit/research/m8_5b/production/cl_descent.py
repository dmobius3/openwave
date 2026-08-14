"""C_L descent and branch-rank extraction on the repaired quotient backend.

WHAT C_L IS, derived rather than inherited.

Under the frozen section 2 action `x -> u x v` the left-invariant frame is carried
by `Ad(v^-1)` alone, so the SU(2)_L factor acts on a one-form by moving the base
point and leaving the frame components alone.  The left Casimir on the one-form
carrier is therefore the scalar Casimir applied COMPONENTWISE:

    C_L = kron(I_3, L_pos)                in the component-major 3N layout

and on `Omega^1 = (+)_n V_n(u) (x) [V_n (x) V_2](v)` it acts as `n(n+2)` on the
`V_n(u)` factor, for every right-factor summand `m` in `{n+2, n, n-2}`.  That `n`
is exactly section 6.1's harmonic level, which "labels the SCALAR factor V_n(u)
throughout, never the resulting representation and never the eigenvalue index".
So `C_L` is the operator that separates the two coexact branches sharing one
eigenvalue, which is the whole reason it is needed.

DESCENT IS NOW STRUCTURAL, NOT NUMERICAL.

    [C_L, T_gamma] = kron(I_3, L) kron(A, P) - kron(A, P) kron(I_3, L)
                   = kron(A, L P) - kron(A, P L)
                   = kron(A, [L, P])

so descent holds EXACTLY whenever the scalar operator commutes exactly with the
node permutation.  On the repaired backend `[L, P] = 0` identically, so this is
zero by construction rather than at a numerical floor.  On the old independent-kNN
backend it was `3e-03`, which is why this work was blocked: the previous failure
was never a property of `C_L`, it was the scalar stencil defect seen through
`C_L`.  Nothing is recalibrated from that number.

THE TWO FIELD SEMANTICS, which the current artifact conflates.

    cluster_degeneracy_count      dim E_lambda, the numerical cluster dimension
    measured_rank_(n, branch)     dim ker( C_L|E_lambda - n(n+2) I )

The producer currently emits the first under the second's name.  That is the
defect being corrected, and the old behaviour is retained here as a REQUIRED
NEGATIVE MUTATION, `mode="cluster_count"`.  It is a real historical defect rather
than a hypothetical one.

`quotient_multiplicity` is not read until after the branch ranks are produced.
The numerical route must reach its own answer first.
"""

import sys

import numpy as np

sys.path.insert(0, "../pilot")

import equivariant_oneform as ofm                                  # noqa: E402
import equivariant_stencils as eqs                                 # noqa: E402
from route_a_twosided import (cloud, relax, reduce_scalar,          # noqa: E402
                              reduce_oneform)

CLUSTER_WINDOW = 0.35          # frozen, section 6.3
BRANCH_WINDOW = 0.35           # the same frozen window, applied to C_L eigenvalues

# FROZEN VACUOUS-CASE CONVENTION, not a tolerance.
#
# `band_leak` is the ratio `|(I - P) C_L P| / |C_L P|`.  Where `C_L|E` is the zero
# operator the denominator is roundoff and the ratio is 0/0.  A zero operator
# preserves every subspace, so the leak is vacuously zero and the ratio is simply
# undefined, not small.  `LEAK_DENOM_FLOOR` substitutes the operator's own scale
# in that case, and the substitution is REPORTED in the record it produces.
#
# It is a definition covering an undefined ratio, NOT a threshold fitted to
# observed residuals.  Nothing about its value was chosen by looking at the
# measured leaks: every well-scaled cluster here sits near 1e-13, thirteen orders
# below anything this constant could reach, so no measurement depends on it.  It
# is frozen now, before any gate consumes it.
LEAK_DENOM_FLOOR = 1e-8

__all__ = ["build_CL", "descent_residual", "band_leak", "branch_ranks",
           "reduced_system"]


def build_CL(Lpos):
    """`kron(I_3, L_pos)`: component index OUTER, node index INNER."""
    return np.kron(np.eye(3), Lpos)


def descent_residual(CL, X, oid, gid, pairs, trials=6, rng_seed=5):
    """Does `C_L` map Gamma-invariant fields to Gamma-invariant fields?

    The invariant projector is built from the FROZEN transport, so this asks the
    question `measured_rank` actually needs rather than a global commutator.
    Both the input and the output invariance are reported: an input that is not
    invariant to begin with makes the output number meaningless.
    """
    N = len(X)
    mult = eqs.group_multiplication_table(pairs)
    node_of = {(int(oid[i]), int(gid[i])): i for i in range(N)}
    Ts = [ofm.transport_matrix(X, oid, gid, pairs, g, mult, node_of, mode="correct")
          for g in range(len(pairs))]

    rng = np.random.default_rng(rng_seed)
    worst_in = worst_out = 0.0
    for _ in range(trials):
        v = rng.normal(size=3 * N)
        w = sum(T @ v for T in Ts) / len(Ts)            # project onto H_Gamma
        nw = np.linalg.norm(w)
        worst_in = max(worst_in,
                       max(np.linalg.norm(T @ w - w) for T in Ts) / nw)
        Cw = CL @ w
        nCw = np.linalg.norm(Cw)
        worst_out = max(worst_out,
                        max(np.linalg.norm(T @ Cw - Cw) for T in Ts) / nCw)
    return worst_in, worst_out


def reduced_system(X, oid, gid, pairs, M, k=110, m=7, p=4):
    """The quotient operators, all three through the SAME production reduction."""
    Lpos, H = ofm.build_operators_equivariant(X, oid, gid, pairs, k, m, p)
    CL = build_CL(Lpos)
    N = len(X)
    return {
        "Lpos": Lpos, "H": H, "CL": CL,
        "scal_red": reduce_scalar(Lpos, oid, M),
        "H_red": reduce_oneform(H, oid, gid, pairs, M, N),
        "CL_red": reduce_oneform(CL, oid, gid, pairs, M, N),
    }


def _cluster_basis(ev, V, lam, window=CLUSTER_WINDOW, rtol=1e-8):
    """Rank-revealing real basis of a numerical eigenvalue cluster (6.1 rule 2)."""
    idx = np.where(np.abs(np.real(ev) - lam) <= window)[0]
    if len(idx) == 0:
        return np.zeros((V.shape[0], 0)), 0
    B = np.hstack([np.real(V[:, idx]), np.imag(V[:, idx])])
    U, s, _ = np.linalg.svd(B, full_matrices=False)
    r = int(np.sum(s > rtol * s[0])) if s.size else 0
    return U[:, :r], len(idx)


def band_leak(CL_red, ev, V, lam):
    """Does `C_L` PRESERVE the numerical cluster?  If not, restricting it is void.

    `|(I - P) C_L P| / |C_L P|` on the cluster's own orthonormal basis.  This is
    the quotient-space counterpart of the certified-band invariance earned on the
    unquotiented sphere; it is EARNED here, not inherited from there.

    THE DENOMINATOR NEEDS A FLOOR, and leaving it out produced a real false alarm.
    At `lambda = 0` the restriction `C_L|E` is the zero operator, so `|C_L P|` is
    roundoff and the ratio is 0/0, which came out as `1.0e+00`: an apparent total
    leak at the one coordinate that is most obviously fine.  A zero operator
    preserves every subspace, so the denominator is floored at the operator's own
    scale and the substitution is REPORTED rather than silent.
    """
    U, _ = _cluster_basis(ev, V, lam)
    if U.shape[1] == 0:
        return None
    Y = CL_red @ U
    out = float(np.linalg.norm(Y - U @ (U.T @ Y)))
    den = float(np.linalg.norm(Y))
    floor = LEAK_DENOM_FLOOR * float(np.linalg.norm(CL_red, 2))
    if den > floor:
        return {"value": out / den, "denominator": "|C_L P|",
                "denominator_state": "well scaled"}
    return {"value": out / floor if floor > 0 else 0.0,
            "denominator": f"{LEAK_DENOM_FLOOR:g} |C_L|_2",
            "denominator_state": ("C_L|E is numerically the zero operator, which "
                                  "preserves every subspace; |C_L P| is roundoff "
                                  "and would give a meaningless 0/0 ratio")}


def branch_ranks(CL_red, ev, V, lam, mode="branch"):
    """Ranks of `C_L|E_lambda` at each predicted left-Casimir value `n(n+2)`.

    `mode="cluster_count"` is the REQUIRED NEGATIVE MUTATION: it reports the
    cluster dimension for every branch, which is the producer's current defect.

    Two independent readings of the same restriction are returned, and they must
    agree.  The eigenvalue reading uses the frozen 0.35 window; the nullity
    reading is a rank-revealing SVD of `C_L|E - c I`.  Agreement is what makes
    the branch rank a measurement rather than a rounding choice.
    """
    U, n_raw = _cluster_basis(ev, V, lam)
    dim = U.shape[1]
    if dim == 0:
        return {"cluster_dim": 0, "raw_cluster_indices": n_raw, "branches": {}}

    R = U.T @ CL_red @ U
    w = np.linalg.eigvals(R)
    wr = np.sort(np.real(w))
    spread = float(np.abs(np.imag(w)).max())

    branches = {}
    for n in range(0, 12):
        c = float(n * (n + 2))
        by_eig = int(np.sum(np.abs(wr - c) <= BRANCH_WINDOW))
        s = np.linalg.svd(R - c * np.eye(dim), compute_uv=False)
        by_null = int(np.sum(s <= 1e-6 * max(float(s[0]), 1.0)))
        if by_eig == 0 and by_null == 0:
            continue
        branches[n] = {
            "left_casimir": c,
            "rank_by_eigenvalue_window": by_eig,
            "rank_by_svd_nullity": by_null,
            "agree": by_eig == by_null,
            "reported": dim if mode == "cluster_count" else by_eig,
        }
    return {"cluster_dim": dim, "raw_cluster_indices": n_raw,
            "cl_eigenvalues": [float(x) for x in wr],
            "imaginary_spread": spread, "branches": branches}
