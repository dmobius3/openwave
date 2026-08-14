"""Route (a)'s two-sided ONE-FORM action, derived pointwise.

CARRIER.  A one-form is `omega = f_j sigma^j` in the LEFT-INVARIANT coframe, so
its data is a coefficient function `f : S^3 -> R^3`.  The carrier at harmonic
level n is therefore

    R^3  tensor  (level-n scalar coefficient space)

DERIVED ACTION (section 2's pullback law, `(F* omega)(x) = Ad(v) f(u x v)`):

    T_n(u, v)  =  Ad(v)  kron  M_n(u, v)^T          [frame SLOW]
               =  M_n(u, v)^T  kron  Ad(v)          [frame FAST]

the same operator under the two index layouts, matching the geometric pullback
at 3.2e-15.  Neither form was assumed: `derive_oneform_action.py` measured all
five candidates against the geometric object.

TWO INDEPENDENT HAZARDS, KEPT AS SEPARATE PREDICATES.  They are different
mistakes with different fixes, and a single "pointwise action failed" boolean
would say the composite operator is wrong without saying which convention was
under test:

    HARMONIC_COEFFICIENT_ACTION_MATCHES_PULLBACK
        M_n^T is the action on COEFFICIENTS; M_n acts on BASIS VALUES.
        Substituting M_n misses by ~1.8e+01.

    FRAME_ACTION_MATCHES_AD_V
        the frame factor is Ad(v), not Ad(v^-1) = Ad(v)^T.
        Substituting Ad(v^-1) misses by ~8.5e+00.

BOTH are invisible to characters (`chi(Ad v) = chi(Ad v^-1)` exactly), invisible
to invariant dimensions, and invisible to cross-route multiplicity agreement.
Only the pointwise evaluation separates them.

NOT YET SECTORED.  This module gives the action and its invariants on the FULL
carrier.  The frozen records need exact and coexact separately, and section 6.1
fixes that decomposition as `V_n tensor V_2 = V_{n+2} + V_n + V_{n-2}` on the
RIGHT factor.  A correct total invariant dimension can still hide a wrong
allocation between sectors, so the sector maps must be built from the frozen
convention rather than by splitting a total afterwards.  That is the next step
and `quotient_multiplicity` for one-form records stays unpopulated until it
exists.
"""

import sys

import numpy as np

sys.path.insert(0, "../pilot")
from route_a_nonabelian import quat_to_su2, sym_power     # noqa: E402
from route_a_twosided import Ad, qmul                     # noqa: E402

import route_a_repn as repn                               # noqa: E402

__all__ = ["oneform_operator", "oneform_invariant_dim_and_basis",
           "verify_harmonic_coefficient_action", "verify_frame_action",
           "FRAME_LAYOUT"]

FRAME_LAYOUT = "slow"          # Ad(v) kron M^T


def oneform_operator(u, v, n, harmonic_mutation=False, frame_mutation=False):
    """T_n(u, v) on the one-form carrier, frame slow.

    The two mutation flags exist ONLY for the harness and are never set in
    production.  They are separate arguments precisely so a test can redden one
    predicate without disturbing the other.
    """
    M = (repn.level_operator(u, v, n) if harmonic_mutation
         else repn.coefficient_operator(u, v, n))
    A = np.asarray(Ad(v), dtype=float)
    if frame_mutation:
        A = np.linalg.inv(A)
    return np.kron(A, M)


def _basis_values(x, n):
    return sym_power(quat_to_su2(np.asarray(x, dtype=float)),
                     n).reshape(-1, order=repn.VEC_ORDER)


def _pullback_pointwise(coeffs, n, u, v, x):
    """(F* omega)(x) = Ad(v) f(u x v), computed geometrically."""
    uxv = qmul(qmul(u, x), v)
    vals = np.array([coeffs[j] @ _basis_values(uxv, n) for j in range(3)])
    return np.asarray(Ad(v), dtype=float) @ vals


def _residual(n, trials, rng, harmonic_mutation=False, frame_mutation=False):
    dim = (n + 1) ** 2
    worst = 0.0
    for _ in range(trials):
        q = [rng.normal(size=4) for _ in range(3)]
        u, v, x = [a / np.linalg.norm(a) for a in q]
        c = rng.normal(size=(3, dim)) + 1j * rng.normal(size=(3, dim))
        truth = _pullback_pointwise(c, n, u, v, x)
        T = oneform_operator(u, v, n, harmonic_mutation, frame_mutation)
        got = (T @ c.reshape(-1)).reshape(3, dim) @ _basis_values(x, n)
        worst = max(worst, float(np.abs(got - truth).max()))
    return worst


def verify_harmonic_coefficient_action(levels=(2, 3, 4), trials=20, seed=20260811,
                                       tol=1e-10):
    """HARMONIC_COEFFICIENT_ACTION_MATCHES_PULLBACK, with its own mutation."""
    rng = np.random.default_rng(seed)
    ok = max(_residual(n, trials, rng) for n in levels)
    rng = np.random.default_rng(seed)
    mut = min(_residual(n, trials, rng, harmonic_mutation=True) for n in levels)
    return {"predicate": "HARMONIC_COEFFICIENT_ACTION_MATCHES_PULLBACK",
            "shipped_residual": ok, "mutation_residual": mut,
            "mutation": "M_n^T -> M_n (basis-value action)",
            "pass": bool(ok < tol and mut > tol)}


def verify_frame_action(levels=(2, 3, 4), trials=20, seed=20260811, tol=1e-10):
    """FRAME_ACTION_MATCHES_AD_V, with its own mutation."""
    rng = np.random.default_rng(seed)
    ok = max(_residual(n, trials, rng) for n in levels)
    rng = np.random.default_rng(seed)
    mut = min(_residual(n, trials, rng, frame_mutation=True) for n in levels)
    return {"predicate": "FRAME_ACTION_MATCHES_AD_V",
            "shipped_residual": ok, "mutation_residual": mut,
            "mutation": "Ad(v) -> Ad(v^-1) = Ad(v)^T",
            "pass": bool(ok < tol and mut > tol)}


def oneform_invariant_dim_and_basis(pairs, n, tol_rel=1e-8, **mut):
    """Invariant dimension and basis of the FULL one-form carrier at level n.

    NOT a per-sector quantity.  See the module docstring: the exact/coexact
    split must come from the frozen section 6.1 decomposition, not from
    apportioning this total.
    """
    dim = 3 * (n + 1) ** 2
    blocks = [oneform_operator(u, v, n, **mut) - np.eye(dim) for (u, v) in pairs]
    A = np.vstack(blocks)
    _, s, Vh = np.linalg.svd(A)
    cutoff = max(tol_rel * (s[0] if s.size else 0.0), 1e-12)
    k = int(np.sum(s < cutoff))
    basis = Vh.conj().T[:, dim - k:] if k else np.zeros((dim, 0), dtype=complex)
    return k, basis


# --- certified sector projectors and sector-specific invariant ranks ----------

from math import comb                                        # noqa: E402

_AX = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
_EPS = np.zeros((3, 3, 3))
for _i, _j, _k in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
    _EPS[_i, _j, _k] = 1.0
    _EPS[_i, _k, _j] = -1.0


def H_scalar(n):
    """Invariant Gram on the monomial basis: diag(1/C(n,k)).  Verified 1.8e-15."""
    return np.diag([1.0 / comb(n, k) for k in range(n + 1)])


def H_oneform_carrier(n):
    """I_3 kron (H_n kron H_n^-1).  The transposed u-factor carries the INVERSE."""
    return np.kron(np.eye(3), np.kron(H_scalar(n), np.linalg.inv(H_scalar(n))))


def _axis_quat(a, t):
    a = np.asarray(a, float) / np.linalg.norm(a)
    return np.concatenate([[np.cos(t / 2)], np.sin(t / 2) * a])


def _casimir(n, h=1e-5):
    def RA(a, t):
        v = _axis_quat(a, t)
        return np.kron(np.asarray(Ad(v), float), sym_power(quat_to_su2(v), n))
    J = [(RA(a, h) - RA(a, -h)) / (2 * h) for a in _AX]
    return sum(Ja @ Ja for Ja in J)


def _casimir_of_Vm(m, h=1e-5):
    def R(a, t):
        return sym_power(quat_to_su2(_axis_quat(a, t)), m)
    J = [(R(a, h) - R(a, -h)) / (2 * h) for a in _AX]
    return float(np.real(np.linalg.eigvals(sum(Ja @ Ja for Ja in J)).mean()))


def sector_projector(n, m):
    """Algebraic spectral projector onto the V_m summand of (frame x V_n(v)).

    Labels come from the CALIBRATED Casimir eigenvalue, never from dimension.
    Spectral form V E V^-1, not Q Q^H: the monomial basis is non-orthonormal, so
    the Euclidean form fails P_i P_j = 0 and equivariance at n >= 2 while still
    having the right ranks.
    """
    C = _casimir(n)
    w, V = np.linalg.eig(C)
    Vi = np.linalg.inv(V)
    E = np.zeros_like(C)
    for i, e in enumerate(w.real):
        if min(range(0, n + 4), key=lambda mm: abs(_casimir_of_Vm(mm) - e)) == m:
            E[i, i] = 1.0
    return V @ E @ Vi


def sector_invariant_rank(pairs, n, m, tol_rel=1e-8):
    """Character-free invariant rank INSIDE the certified sector V_m.

    Builds an H-orthonormal basis of ran(P_m) lifted to the full carrier, then
    restricts every group element to that basis and takes the stacked nullity
    there.  "Rank inside the sector" is therefore literal: the computation never
    forms the full-carrier fixed space and apportions it.
    """
    P = np.kron(sector_projector(n, m), np.eye(n + 1))
    H = H_oneform_carrier(n)
    S = np.linalg.cholesky(H).conj().T                 # H = S^* S

    w, V = np.linalg.eig(P)
    idx = [i for i, x in enumerate(w.real) if abs(x - 1) < 1e-6]
    if not idx:
        return 0
    B = V[:, idx]
    Q, _ = np.linalg.qr(S @ B)                          # orthonormal in the H-metric
    B = np.linalg.solve(S, Q)                           # back to carrier coordinates

    blocks = []
    for (u, v) in pairs:
        R = oneform_operator(u, v, n)
        Rm = B.conj().T @ H @ R @ B                     # restriction, H-orthonormal basis
        blocks.append(Rm - np.eye(Rm.shape[0]))
    A = np.vstack(blocks)
    s = np.linalg.svd(A, compute_uv=False)
    cutoff = max(tol_rel * (s[0] if s.size else 0.0), 1e-12)
    return int(np.sum(s < cutoff))


def admissible_sectors(n):
    """The frozen section 6.1 summands that EXIST at level n, by level.

    Constructed explicitly rather than by a general formula, so no phantom
    V_{n-2} sector is manufactured at the bottom of the tower.
    """
    out = []
    if n >= 1:
        out.append(("oneform_exact", n, n * (n + 2)))
    out.append(("oneform_coexact_up", n + 2, (n + 2) ** 2))
    if n >= 2:
        out.append(("oneform_coexact_down", n - 2, n * n))
    return out
