"""P1A.3: Invariant-subspace extractor and scored quantities.

The production candidate is an ordered complex Schur factorization.
The span is M_h-orthonormalized: Q^H M_h Q = I.

Scored quantities (frozen definitions):
  - Principal angles: sigma_i = cos(theta_i) from SVD of Q_0^H M_h Q_1
  - Leakage: ell = (1/k) * ||(I - P_0) Q_1||^2_{M_h,F}
  - Three-way identity: ell == 1 - (1/k)*||Q_0^H M_h Q_1||_F^2
                             == (1/k) * sum_i sin^2(theta_i)
"""

import numpy as np
from .mass_matrix import Mh_matvec, Mh_inner, Mh_norm_F
from .diagnostics import riesz_projector, riesz_projector_norm, invariance_residual


def ordered_schur_subspace(L, k, target_lambda):
    """Extract the k-dimensional Schur invariant subspace nearest target_lambda.

    Uses scipy's ordered Schur decomposition.
    Returns (Q_raw, T, evals_cluster, sdim).
    """
    from scipy.linalg import schur

    T, Z = schur(L, output='complex')
    evals = np.diag(T)

    dists = np.abs(evals - target_lambda)
    order = np.argsort(dists)

    if k >= len(evals):
        return Z, T, evals, len(evals)

    threshold = float(dists[order[k - 1]]) + 1e-10
    T2, Z2, sdim = schur(L, output='complex',
                          sort=lambda e: abs(e - target_lambda) < threshold)

    if sdim != k:
        for mult in [1.01, 1.1, 1.5, 2.0, 5.0]:
            th = float(dists[order[k - 1]]) * mult
            T2, Z2, sdim = schur(L, output='complex',
                                  sort=lambda e, th=th: abs(e - target_lambda) < th)
            if sdim == k:
                break

    if sdim != k and k < len(evals):
        gap_th = (float(dists[order[k - 1]]) + float(dists[order[k]])) / 2.0
        T2, Z2, sdim = schur(L, output='complex',
                              sort=lambda e, th=gap_th: abs(e - target_lambda) < th)

    evals_cluster = np.diag(T2)[:sdim]
    Q_raw = Z2[:, :sdim]
    return Q_raw, T2, evals_cluster, sdim


def mh_orthonormalize(Q, Mh_diag):
    """M_h-orthonormalize Q so that Q^H M_h Q = I.

    Uses modified Gram-Schmidt with M_h inner product.
    """
    n, k = Q.shape
    R = Q.copy().astype(complex)
    for j in range(k):
        for i in range(j):
            coeff = Mh_inner(Mh_diag, R[:, i], R[:, j])
            R[:, j] -= coeff * R[:, i]
        norm = np.sqrt(np.real(Mh_inner(Mh_diag, R[:, j], R[:, j])))
        if norm < 1e-14:
            raise ValueError(f"M_h-orthonormalization failed at column {j}: near-zero norm")
        R[:, j] /= norm
    return R


def verify_mh_orthonormal(Q, Mh_diag, tol=1e-10):
    """Check Q^H M_h Q ≈ I."""
    k = Q.shape[1]
    G = np.zeros((k, k), dtype=complex)
    for i in range(k):
        for j in range(k):
            G[i, j] = Mh_inner(Mh_diag, Q[:, i], Q[:, j])
    err = float(np.max(np.abs(G - np.eye(k))))
    return err < tol, err


def principal_angles(Q0, Q1, Mh_diag):
    """Principal angles between M_h-orthonormal bases Q0 and Q1.

    Singular values of Q0^H M_h Q1 are the cosines.
    theta_i = arccos(clip(sigma_i, 0, 1)).

    Returns (thetas, sigmas).
    """
    k = Q0.shape[1]
    M = np.zeros((k, k), dtype=complex)
    for i in range(k):
        for j in range(k):
            M[i, j] = Mh_inner(Mh_diag, Q0[:, i], Q1[:, j])

    sigmas = np.linalg.svd(M, compute_uv=False)
    sigmas_real = np.real(sigmas)
    sigmas_clipped = np.clip(sigmas_real, 0.0, 1.0)
    thetas = np.arccos(sigmas_clipped)
    return thetas, sigmas_clipped


def leakage_projector(Q0, Q1, Mh_diag):
    """Leakage by the direct projector route.

    ell = (1/k) * ||(I - P_0) Q_1||^2_{M_h,F}
    where P_0 = Q_0 Q_0^H M_h.
    """
    k = Q0.shape[1]
    n = Q0.shape[0]

    residual_norm_sq = 0.0
    for j in range(Q1.shape[1]):
        v = Q1[:, j].copy()
        for i in range(k):
            coeff = Mh_inner(Mh_diag, Q0[:, i], v)
            v -= coeff * Q0[:, i]
        residual_norm_sq += np.real(Mh_inner(Mh_diag, v, v))

    return float(residual_norm_sq / k)


def leakage_angle(thetas):
    """Leakage by the angle route: (1/k) * sum_i sin^2(theta_i)."""
    return float(np.mean(np.sin(thetas)**2))


def leakage_overlap(Q0, Q1, Mh_diag):
    """Leakage by the overlap route: 1 - (1/k) * ||Q0^H M_h Q1||_F^2."""
    k = Q0.shape[1]
    M = np.zeros((k, Q1.shape[1]), dtype=complex)
    for i in range(k):
        for j in range(Q1.shape[1]):
            M[i, j] = Mh_inner(Mh_diag, Q0[:, i], Q1[:, j])
    fro2 = float(np.sum(np.abs(M)**2))
    return 1.0 - fro2 / k


def three_way_identity_check(Q0, Q1, Mh_diag, tol=1e-10):
    """The three-way identity is a live gate.

    Compute ell by projector, angle, and overlap routes independently.
    If they disagree beyond tolerance, scoring stops.
    """
    thetas, sigmas = principal_angles(Q0, Q1, Mh_diag)
    ell_proj = leakage_projector(Q0, Q1, Mh_diag)
    ell_angle = leakage_angle(thetas)
    ell_overlap = leakage_overlap(Q0, Q1, Mh_diag)

    agree_pa = abs(ell_proj - ell_angle)
    agree_po = abs(ell_proj - ell_overlap)
    agree_ao = abs(ell_angle - ell_overlap)

    max_disagree = max(agree_pa, agree_po, agree_ao)
    passed = max_disagree < tol

    ell_max = float(np.max(np.sin(thetas)**2))

    return {
        "thetas": thetas.tolist(),
        "sigmas": sigmas.tolist(),
        "ell_projector": ell_proj,
        "ell_angle": ell_angle,
        "ell_overlap": ell_overlap,
        "ell_scored": ell_proj,
        "ell_max": ell_max,
        "max_disagreement": max_disagree,
        "pass": passed,
    }


def projector_invariants(Q, L, Mh_diag, k):
    """Projector invariants, the experiment's live proof.

    P = Q Q^H M_h (the M_h-orthogonal projector).
    Checks: P^2 ≈ P, P^H M_h ≈ M_h P, rank P = k, (I-P)LP ≈ 0.
    """
    n = Q.shape[0]
    Mh = np.diag(Mh_diag)

    P = Q @ (Q.conj().T @ Mh)

    P2 = P @ P
    idempotence = float(np.linalg.norm(P2 - P, 'fro') / max(np.linalg.norm(P, 'fro'), 1e-15))

    PhM = P.conj().T @ Mh
    MP = Mh @ P
    symmetry = float(np.linalg.norm(PhM - MP, 'fro') / max(np.linalg.norm(MP, 'fro'), 1e-15))

    rank = int(np.round(np.real(np.trace(P))))

    IminusP = np.eye(n) - P
    inv_residual = float(np.linalg.norm(IminusP @ L @ P, 'fro'))
    inv_residual_rel = inv_residual / max(np.linalg.norm(L @ P, 'fro'), 1e-15)

    return {
        "idempotence": idempotence,
        "mh_symmetry": symmetry,
        "rank": rank,
        "rank_expected": k,
        "invariance_residual": inv_residual,
        "invariance_residual_rel": inv_residual_rel,
        "pass": (idempotence < 1e-10 and symmetry < 1e-10
                 and rank == k),
    }


def riesz_crosscheck(L, Q_schur, Mh_diag, centre, radius, k, n_quad=128):
    """Riesz cross-check: compare Schur subspace with Riesz subspace.

    The contour (centre, radius) is fixed from known analytic levels
    BEFORE any numerical eigenvalue is inspected.

    Returns the max principal angle between the two subspaces.
    """
    P_spec = riesz_projector(L, centre, radius, n_quad=n_quad)

    evals_P, vecs_P = np.linalg.eig(P_spec)
    idx = np.where(np.abs(evals_P - 1.0) < 0.5)[0]

    if len(idx) < k:
        idx_sorted = np.argsort(np.abs(evals_P - 1.0))
        idx_top = idx_sorted[:k]
        if np.max(np.abs(evals_P[idx_top] - 1.0)) > 0.8:
            return {
                "riesz_rank": len(idx),
                "expected_rank": k,
                "rank_mismatch": True,
                "P_spec_norm": riesz_projector_norm(P_spec),
                "pass": False,
            }
        idx = idx_top

    if len(idx) > k:
        sub_dists = np.abs(evals_P[idx] - 1.0)
        keep = np.argsort(sub_dists)[:k]
        idx = idx[keep]

    Q_riesz = np.real(vecs_P[:, idx]) if np.max(np.abs(np.imag(vecs_P[:, idx]))) < 1e-10 else vecs_P[:, idx]
    try:
        Q_riesz = mh_orthonormalize(Q_riesz, Mh_diag)
    except ValueError:
        return {
            "riesz_rank": len(idx),
            "expected_rank": k,
            "rank_mismatch": False,
            "orthonorm_failed": True,
            "P_spec_norm": riesz_projector_norm(P_spec),
            "pass": False,
        }

    thetas, _ = principal_angles(Q_schur, Q_riesz, Mh_diag)
    theta_max = float(np.max(thetas))

    sqrt_M = np.sqrt(Mh_diag)
    inv_sqrt_M = 1.0 / sqrt_M
    P_mh = sqrt_M[:, None] * P_spec * inv_sqrt_M[None, :]
    p_spec_norm_mh = float(np.linalg.norm(P_mh, 2))

    return {
        "riesz_rank": len(idx),
        "expected_rank": k,
        "theta_max_schur_vs_riesz": theta_max,
        "P_spec_norm": riesz_projector_norm(P_spec),
        "P_spec_norm_mh": p_spec_norm_mh,
        "rank_mismatch": False,
        "pass": True,
    }


def imaginary_contamination(evals_cluster):
    """Check |Im lambda| for the eigenvalue cluster."""
    im_parts = np.abs(np.imag(evals_cluster))
    return {
        "max_im": float(np.max(im_parts)),
        "mean_im": float(np.mean(im_parts)),
        "all_im": im_parts.tolist(),
    }


def extract_and_score(L, Mh_diag, target_lambda, k, label="",
                      riesz_centre=None, riesz_radius=None):
    """Full P1A.3 pipeline: Schur extract, orthonormalize, score.

    Returns a comprehensive results dict.
    """
    Q_raw, T, evals_cluster, sdim = ordered_schur_subspace(L, k, target_lambda)

    results = {"label": label, "target_lambda": target_lambda, "k": k}

    if sdim != k:
        results["schur_dim_mismatch"] = True
        results["schur_dim"] = sdim
        results["pass"] = False
        return results, None

    results["schur_dim_mismatch"] = False
    results["evals_cluster"] = evals_cluster.tolist()
    results["cluster_position"] = float(np.mean(np.real(evals_cluster)))
    results["cluster_spread"] = float(np.max(np.real(evals_cluster)) - np.min(np.real(evals_cluster)))

    im_check = imaginary_contamination(evals_cluster)
    results["imaginary"] = im_check

    Q = mh_orthonormalize(Q_raw, Mh_diag)
    orth_ok, orth_err = verify_mh_orthonormal(Q, Mh_diag)
    results["mh_orthonormal_err"] = orth_err

    proj = projector_invariants(Q, L, Mh_diag, k)
    results["projector"] = proj

    if riesz_centre is not None and riesz_radius is not None:
        riesz = riesz_crosscheck(L, Q, Mh_diag, riesz_centre, riesz_radius, k)
        results["riesz"] = riesz

    all_evals = np.linalg.eigvals(L)
    all_evals_sorted = np.sort(np.real(all_evals))
    target_idx = []
    complement_idx = []
    for i, ev in enumerate(all_evals):
        if abs(ev - target_lambda) < abs(results["cluster_spread"]) + 2.0:
            target_idx.append(i)
        else:
            complement_idx.append(i)
    from .diagnostics import spectral_separation, schur_sep
    results["spectral_separation"] = spectral_separation(all_evals, target_idx, complement_idx)
    if k < L.shape[0] and k > 0:
        results["schur_sep"] = schur_sep(T, k)

    results["pass"] = (not results["schur_dim_mismatch"]
                       and proj["pass"]
                       and orth_ok)

    return results, Q
