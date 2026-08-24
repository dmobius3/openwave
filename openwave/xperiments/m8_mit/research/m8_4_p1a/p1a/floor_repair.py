"""P1A.4b: Corrected eigensolver-uncertainty floor.

The P1A.4a floor F = eps_mach * ||L|| * max(||P_spec||, 1) is a one-matvec
backward-error scale (~1e-12).  Measured roundoff-limited contamination sits
near 1e-10, so genuinely clean observations are classified ABOVE_FLOOR and
enter the pooled regression nine orders below the converging group.

The defect: backward error is not eigenvalue uncertainty for a non-normal
matrix.  It must be amplified by spectral sensitivity.

DERIVATION
----------
The bundle Laplacian L is M_h-self-adjoint up to discretization error.
Its true eigenvalues are real; computed imaginary parts arise entirely
from numerical perturbation E where L_computed = L_true + E.

For a cluster C of eigenvalues with spectral projector P_C, first-order
perturbation theory (Stewart, Kato) gives:

    |delta lambda_j| <= ||P_C||_2 * ||E||_2

for each eigenvalue lambda_j in C.

The perturbation norm is bounded by the backward error:

    ||E||_2 <= max(max_j beta_j, eps_mach * ||L||_2)

where beta_j = ||L v_j - lambda_j v_j||_2 / ||v_j||_2 is the eigenpair
residual from the computed eigendecomposition.

Therefore the eigensolver-uncertainty floor is:

    F = ||P_spec||_2 * max(max_j beta_j, eps_mach * ||L||_2)

The spectral separation sep(C, sigma(L)\\C) validates the bound:
||P_spec|| * backward_error << sep must hold for the first-order
approximation to be trustworthy.

No observed |Im lambda| enters the derivation.  The floor is computed
from the eigensolver's own residuals, the operator norm, and the
spectral projector norm.
"""

import hashlib
import numpy as np

FLOOR_FORMULA = (
    "F = ||P_spec||_2 * max(max_j(||L*v_j - lambda_j*v_j||_2 / ||v_j||_2),"
    " eps_mach * ||L||_2)"
)
FLOOR_FORMULA_HASH = hashlib.sha256(FLOOR_FORMULA.encode()).hexdigest()

EPS_MACH = np.finfo(float).eps

POWER_GATE_FACTOR = 10.0
DISCRIMINATION_FACTOR = 10.0


def cluster_residuals(L, target_lambda, k):
    """Eigenpair residuals and spectral separation for the target cluster.

    Returns dict with max_beta, betas, cluster_evals, complement_evals, sep.
    """
    evals, vecs = np.linalg.eig(L)
    dists = np.abs(evals - target_lambda)
    order = np.argsort(dists)

    cluster_idx = order[:k]
    complement_idx = order[k:]

    betas = []
    for j in cluster_idx:
        v = vecs[:, j]
        r = L @ v - evals[j] * v
        beta = float(np.linalg.norm(r) / max(np.linalg.norm(v), 1e-300))
        betas.append(beta)

    cluster_evals = evals[cluster_idx]
    complement_evals = evals[complement_idx]

    if len(complement_evals) > 0:
        sep = float(np.min(np.abs(
            cluster_evals[:, None] - complement_evals[None, :])))
    else:
        sep = float('inf')

    return {
        "max_beta": float(np.max(betas)),
        "betas": betas,
        "cluster_evals": cluster_evals,
        "complement_evals": complement_evals,
        "sep": sep,
    }


def eigensolver_floor(L, target_lambda, k, P_spec_norm):
    """Compute the eigensolver-uncertainty floor.

    F = ||P_spec|| * max(max_j beta_j, eps_mach * ||L||)
    """
    cr = cluster_residuals(L, target_lambda, k)
    L_norm = float(np.linalg.norm(L, 2))
    backward_error = max(cr["max_beta"], EPS_MACH * L_norm)
    F = P_spec_norm * backward_error
    perturbation_over_sep = (
        P_spec_norm * backward_error / cr["sep"]
        if cr["sep"] > 0 else float('inf'))

    return F, {
        "max_beta": cr["max_beta"],
        "L_norm": L_norm,
        "eps_L": EPS_MACH * L_norm,
        "backward_error": backward_error,
        "P_spec_norm": P_spec_norm,
        "sep": cr["sep"],
        "perturbation_over_sep": perturbation_over_sep,
        "floor": F,
    }


def power_gate(fit):
    """E(h) < 10 * T(h) requires exp(r_max) < 10."""
    if fit.get("status") != "ok":
        return False, "no valid fit"
    if fit.get("p", 0) <= 0:
        return False, "p <= 0"
    exp_rmax = np.exp(fit["r_max"])
    passed = exp_rmax < POWER_GATE_FACTOR
    reason = (f"exp(r_max) = {exp_rmax:.2f} < {POWER_GATE_FACTOR}"
              if passed else
              f"exp(r_max) = {exp_rmax:.2f} >= {POWER_GATE_FACTOR}")
    return passed, reason


def trend_value(s_rho, A, p, h):
    """T(h) = s_rho * A * h^p — central trend without margin."""
    return s_rho * A * h ** p


def mutate_cloud_for_gate_test(X, oid, n_orbits_to_cluster=5,
                                cluster_fraction=0.95):
    """Deterministic cloud mutation: cluster several orbits toward centroids.

    Moves nodes in the specified seed orbits toward their orbit centroid
    on S^3, collapsing separation radius and violating quasi-uniformity.
    Multiple orbits at high fraction ensures the gate trips.
    """
    X_mut = X.copy()
    unique_orbits = np.unique(oid)
    targets = unique_orbits[:n_orbits_to_cluster]

    for seed_id in targets:
        mask = (oid == seed_id)
        if not np.any(mask):
            continue
        pts = X_mut[mask]
        centroid = pts.mean(axis=0)
        centroid /= np.linalg.norm(centroid)
        X_mut[mask] = ((1.0 - cluster_fraction) * pts
                       + cluster_fraction * centroid)
        norms = np.linalg.norm(X_mut[mask], axis=1, keepdims=True)
        X_mut[mask] /= norms
    return X_mut
