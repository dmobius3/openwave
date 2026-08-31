"""Galerkin projector: cubic nonlinearity P_N[<psi,psi>psi] via the transform route.

Uses the fast Hopf-FFT route (fast_transform.py) for production-rule projections
at N >= 8, falling back to direct quadrature for small N and monitor-rule projections.
"""
import numpy as np
from build.quadrature import (
    production_rule, monitor_rule, synthesis, analysis,
)
from build.fast_transform import fast_project_cubic, fast_synthesis, fast_analysis
from build.sections import build_basis_object, total_modes
from build.group import multiplicity, DIMS


_FAST_THRESHOLD = 8


def project_cubic(coeffs, basis_obj, N, rho='R0', rule='production'):
    """Compute P_N[|psi|^2 psi] via the transform route.

    Parameters
    ----------
    coeffs : complex 1-D array, length total_modes(rho, N)
    basis_obj : dict from build_basis_object
    N : int, spectral cutoff
    rho : str
    rule : 'production' (4N) or 'monitor' (6N)

    Returns
    -------
    projected_coeffs : complex 1-D array
    """
    if rule == 'production' and N >= _FAST_THRESHOLD:
        return fast_project_cubic(coeffs, basis_obj, N, rho=rho)

    if rule == 'production':
        nodes, weights = production_rule(N)
    elif rule == 'monitor':
        nodes, weights = monitor_rule(N)
    else:
        raise ValueError(f"unknown rule: {rule}")

    vals = synthesis(coeffs, basis_obj, nodes, rho=rho)
    abs2 = np.sum(np.abs(vals) ** 2, axis=1, keepdims=True)
    cubic_vals = abs2 * vals
    return analysis(cubic_vals, basis_obj, nodes, weights)


def project_cubic_direct(coeffs, basis_obj, N, rho='R0', rule='production'):
    """Direct quadrature route (no FFT). For dual-route comparison."""
    if rule == 'production':
        nodes, weights = production_rule(N)
    elif rule == 'monitor':
        nodes, weights = monitor_rule(N)
    else:
        raise ValueError(f"unknown rule: {rule}")

    vals = synthesis(coeffs, basis_obj, nodes, rho=rho)
    abs2 = np.sum(np.abs(vals) ** 2, axis=1, keepdims=True)
    cubic_vals = abs2 * vals
    return analysis(cubic_vals, basis_obj, nodes, weights)


def cascade_reading(coeffs, basis_obj_production, basis_obj_monitor, N, rho='R0'):
    """Compute C_N = ||P_{N<n<=3N} f(psi)||_2 / ||P_{n<=3N} f(psi)||_2.

    Uses the fast FFT transform with 6N degree for the monitor band.
    """
    D_monitor = 6 * N
    field, u, wu, K = fast_synthesis(coeffs, basis_obj_production, N,
                                     rho=rho, degree=D_monitor)
    abs2 = np.sum(np.abs(field) ** 2, axis=-1, keepdims=True)
    cubic_field = abs2 * field

    full_coeffs = fast_analysis(cubic_field, basis_obj_monitor, N,
                                u, wu, K, rho=rho)
    prod_coeffs = fast_analysis(cubic_field, basis_obj_production, N,
                                u, wu, K, rho=rho)

    n_prod = len(prod_coeffs)
    high_band = full_coeffs.copy()
    high_band[:n_prod] -= prod_coeffs

    norm_total = np.linalg.norm(full_coeffs)
    norm_high = np.linalg.norm(high_band)
    if norm_total < 1e-30:
        return 0.0
    return float(norm_high / norm_total)
