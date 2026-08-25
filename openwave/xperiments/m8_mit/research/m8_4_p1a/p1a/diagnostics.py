"""P1A.2: Non-normality and self-adjointness diagnostics.

Measures the SUBSPACE observable's conditioning, not just individual
eigenvectors. Two projectors are in play:
  P_Q = Q Q^H M_h  — M_h-orthogonal projector (unit norm in M_h)
  P_spec = Riesz    — spectral projector (generally oblique)
"""

import numpy as np
from .mass_matrix import Mh_matvec, Mh_inner, Mh_norm_F


def henrici_departure(L):
    """Henrici departure from normality: sqrt(||L||_F^2 - sum|lambda|^2) / ||L||_F."""
    evals = np.linalg.eigvals(L)
    fro2 = float(np.sum(np.abs(L)**2))
    spec2 = float(np.sum(np.abs(evals)**2))
    if fro2 < 1e-30:
        return 0.0
    return float(np.sqrt(max(0, fro2 - spec2)) / np.sqrt(fro2))


def commutator_departure(L):
    """||[L, L^H]|| / ||L||^2."""
    Lh = L.conj().T
    comm = L @ Lh - Lh @ L
    return float(np.linalg.norm(comm) / np.linalg.norm(L)**2)


def asymmetry_ratio(L):
    """||L - L^H|| / ||L||."""
    return float(np.linalg.norm(L - L.conj().T) / np.linalg.norm(L))


def self_adjointness_residual(L, Mh_diag):
    """Dimensionless self-adjointness residual:

    eps_SA = ||L^H M_h - M_h L||_F / (||L^H M_h||_F + ||M_h L||_F)

    Same formula at every resolution and every sector.
    """
    n = L.shape[0]
    Mh = np.diag(Mh_diag)

    LhM = L.conj().T @ Mh
    ML = Mh @ L

    num = float(np.linalg.norm(LhM - ML, 'fro'))
    den = float(np.linalg.norm(LhM, 'fro') + np.linalg.norm(ML, 'fro'))

    if den < 1e-30:
        return 0.0
    return num / den


def spectral_separation(evals, target_indices, complement_indices):
    """Spectral separation between target cluster and complement.

    Returns min|lambda_target - lambda_complement|.
    """
    if len(target_indices) == 0 or len(complement_indices) == 0:
        return float('inf')
    target = evals[target_indices]
    complement = evals[complement_indices]
    min_sep = float('inf')
    for t in target:
        for c in complement:
            min_sep = min(min_sep, abs(t - c))
    return min_sep


def schur_sep(T, k):
    """Schur sep-type quantity for the top-left k×k block.

    For T in Schur form with blocks T11 (k×k) and T22 ((n-k)×(n-k)),
    sep(T11, T22) = min singular value of the Sylvester operator
    X -> T11 X - X T22.
    """
    n = T.shape[0]
    T11 = T[:k, :k]
    T22 = T[k:, k:]
    m = n - k

    def sylvester_op(X_flat):
        X = X_flat.reshape(k, m)
        return (T11 @ X - X @ T22).ravel()

    from scipy.sparse.linalg import LinearOperator, svds
    op = LinearOperator((k * m, k * m), matvec=sylvester_op)
    try:
        s = svds(op, k=1, which='SM', return_singular_vectors=False)
        return float(s[0])
    except Exception:
        A = np.zeros((k * m, k * m), dtype=complex)
        for j in range(k * m):
            e = np.zeros(k * m, dtype=complex)
            e[j] = 1.0
            A[:, j] = sylvester_op(e)
        return float(np.min(np.linalg.svd(A, compute_uv=False)))


def riesz_projector(L, centre, radius, n_quad=64):
    """Riesz spectral projector via contour integral.

    P_spec = (1/2pi i) oint_Gamma (zI - L)^{-1} dz

    Contour is a circle of given centre and radius.
    n_quad points on the circle, trapezoidal rule.
    """
    n = L.shape[0]
    I = np.eye(n, dtype=complex)
    P = np.zeros((n, n), dtype=complex)
    for j in range(n_quad):
        theta = 2 * np.pi * j / n_quad
        z = centre + radius * np.exp(1j * theta)
        dz = 1j * radius * np.exp(1j * theta) * (2 * np.pi / n_quad)
        resolvent = np.linalg.solve(z * I - L, I)
        P += resolvent * dz
    P /= (2 * np.pi * 1j)
    return P


def riesz_projector_norm(P_spec):
    """||P_spec||_2 — spectral projector norm in Euclidean metric."""
    return float(np.linalg.norm(P_spec, 2))


def invariance_residual(L, Q, Mh_diag):
    """||(I - P_Q) L Q||_{M_h} where P_Q = Q Q^H M_h.

    Measures how well the subspace spanned by Q is invariant under L.
    """
    n = L.shape[0]
    k = Q.shape[1]
    LQ = L @ Q
    QhM = Q.conj().T * Mh_diag[None, :]
    proj_LQ = Q @ (QhM @ LQ)
    residual = LQ - proj_LQ
    return float(Mh_norm_F(Mh_diag, residual))


def eigenvector_condition(L):
    """Eigenvector condition number: cond(V) where L = V Lambda V^{-1}."""
    evals, V = np.linalg.eig(L)
    return float(np.linalg.cond(V))


def run_nonnormality_diagnostics(L, Mh_diag, label=""):
    """Full P1A.2 diagnostic suite for one operator."""
    results = {
        "label": label,
        "henrici": henrici_departure(L),
        "commutator": commutator_departure(L),
        "asymmetry": asymmetry_ratio(L),
        "eps_SA": self_adjointness_residual(L, Mh_diag),
        "eigvec_cond": eigenvector_condition(L),
    }
    return results
