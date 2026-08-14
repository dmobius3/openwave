"""Field-level M5/Faber curvature scale observable.

This module evaluates the regularized hedgehog used by M5.6.4b,

    q0 = r0 / sqrt(r^2 + r0^2),       q = x / sqrt(r^2 + r0^2),
    Gamma_i = q0 d_i q - (d_i q0) q + q x d_i q,
    R_ij = Gamma_i x Gamma_j.

The dimensionless, convention-free shell observable is

    C(rho) = < r^2 sqrt(sum_{i<j} |R_ij|^2) >_shell,   rho = r/r0.

For the analytic profile, before shell averaging,

    |R|(r) = sqrt(r^2 + 3 r0^2) / (r^2 + r0^2)^(3/2),
    C(rho) = rho^2 sqrt(rho^2 + 3) / (rho^2 + 1)^(3/2),
    d log(C) / d log(mu) = -6 / ((rho^2 + 1)(rho^2 + 3)),

where mu/mu0 = 1/rho.  ``C`` is the measured field-amplitude/form-factor
curve.  A source/action normalization is still required to identify it with a
renormalized coupling.  ``coupling_interpretations`` reports two preregistered
normalizations without selecting either from agreement with a comparator.

All finite differences below are nonperiodic centered differences.  Boundary
cells are marked NaN and excluded rather than wrapped with ``np.roll``.
"""

from __future__ import annotations

import numpy as np


def regularized_hedgehog(r0: float, half_extent: float, n: int) -> dict:
    """Return the M5.6.4b q field on a cubic grid."""
    if r0 <= 0 or half_extent <= 0 or n < 5:
        raise ValueError("r0 and half_extent must be positive; n must be >= 5")
    xs = np.linspace(-half_extent, half_extent, n, dtype=float)
    h = float(xs[1] - xs[0])
    x, y, z = np.meshgrid(xs, xs, xs, indexing="ij")
    r2 = x * x + y * y + z * z
    rn = np.sqrt(r2 + r0 * r0)
    q0 = r0 / rn
    q = np.stack((x / rn, y / rn, z / rn), axis=-1)
    return {"h": h, "r": np.sqrt(r2), "q0": q0, "q": q}


def centered_difference(field: np.ndarray, axis: int, h: float) -> np.ndarray:
    """Second-order centered derivative with NaNs on two boundary planes."""
    out = np.full(field.shape, np.nan, dtype=float)
    middle = [slice(None)] * field.ndim
    plus = [slice(None)] * field.ndim
    minus = [slice(None)] * field.ndim
    middle[axis] = slice(1, -1)
    plus[axis] = slice(2, None)
    minus[axis] = slice(None, -2)
    out[tuple(middle)] = (field[tuple(plus)] - field[tuple(minus)]) / (2.0 * h)
    return out


def connection(
    q0: np.ndarray,
    q: np.ndarray,
    h: float,
    *,
    include_cross: bool = True,
    scale: float = 1.0,
) -> list[np.ndarray]:
    """Evaluate the three Gamma_i fields from the displayed equation."""
    gamma = []
    for axis in range(3):
        dq = centered_difference(q, axis, h)
        dq0 = centered_difference(q0, axis, h)
        value = q0[..., None] * dq - dq0[..., None] * q
        if include_cross:
            value = value + np.cross(q, dq)
        gamma.append(scale * value)
    return gamma


def curvature_magnitude(gamma: list[np.ndarray]) -> np.ndarray:
    """Return sqrt(sum_{i<j} |Gamma_i x Gamma_j|^2)."""
    total = np.zeros(gamma[0].shape[:-1], dtype=float)
    for i, j in ((0, 1), (0, 2), (1, 2)):
        rij = np.cross(gamma[i], gamma[j])
        total += np.einsum("...a,...a->...", rij, rij)
    return np.sqrt(total)


def analytic_curvature_magnitude(r: np.ndarray, r0: float) -> np.ndarray:
    """Closed form of the continuum curvature norm for the same q field."""
    return np.sqrt(r * r + 3.0 * r0 * r0) / (r * r + r0 * r0) ** 1.5


def analytic_proxy(rho: np.ndarray) -> np.ndarray:
    """Pointwise continuum C(rho), before finite-width shell averaging."""
    rho = np.asarray(rho, dtype=float)
    return rho**2 * np.sqrt(rho**2 + 3.0) / (rho**2 + 1.0) ** 1.5


def analytic_log_slope(rho: np.ndarray) -> np.ndarray:
    """Exact d log(C) / d log(mu), with mu/mu0 = 1/rho."""
    rho = np.asarray(rho, dtype=float)
    return -6.0 / ((rho**2 + 1.0) * (rho**2 + 3.0))


def shell_profile(
    r: np.ndarray,
    rmag: np.ndarray,
    r0: float,
    rho_centers: np.ndarray,
    half_width: float,
) -> dict:
    """Average r^2 |R| over fixed dimensionless spherical shells."""
    rho = r / r0
    values = []
    counts = []
    for center in np.asarray(rho_centers, dtype=float):
        mask = np.isfinite(rmag) & (rho >= center - half_width) & (rho < center + half_width)
        count = int(np.count_nonzero(mask))
        if count == 0:
            raise ValueError(f"empty shell at rho={center}")
        values.append(float(np.mean(r[mask] ** 2 * rmag[mask])))
        counts.append(count)
    return {"rho": np.asarray(rho_centers), "C": np.asarray(values), "counts": counts}


def numerical_profile(
    *,
    r0: float,
    half_extent_over_r0: float,
    n: int,
    rho_centers: np.ndarray,
    shell_half_width: float,
    include_cross: bool = True,
    connection_scale: float = 1.0,
) -> dict:
    """Build the q field, differentiate it, and return numerical/exact shell curves."""
    grid = regularized_hedgehog(r0, half_extent_over_r0 * r0, n)
    gamma = connection(
        grid["q0"],
        grid["q"],
        grid["h"],
        include_cross=include_cross,
        scale=connection_scale,
    )
    rmag = curvature_magnitude(gamma)
    measured = shell_profile(grid["r"], rmag, r0, rho_centers, shell_half_width)
    exact_rmag = analytic_curvature_magnitude(grid["r"], r0)
    exact_shell = shell_profile(grid["r"], exact_rmag, r0, rho_centers, shell_half_width)
    return {
        "n": n,
        "h_over_r0": grid["h"] / r0,
        "half_extent_over_r0": half_extent_over_r0,
        "rho": measured["rho"],
        "C": measured["C"],
        "C_exact_shell": exact_shell["C"],
        "counts": measured["counts"],
    }


def local_polynomial_derivative(x: np.ndarray, y: np.ndarray, width: int = 5) -> np.ndarray:
    """Differentiate y(x) with overlapping local cubic fits."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.ndim != 1 or y.shape != x.shape or width < 5 or width % 2 == 0:
        raise ValueError("x/y must be 1-D peers and width an odd integer >= 5")
    result = np.empty_like(y)
    half = width // 2
    for index in range(len(x)):
        lo = max(0, min(index - half, len(x) - width))
        hi = lo + width
        shifted = x[lo:hi] - x[index]
        degree = min(3, width - 1)
        coeff = np.polynomial.polynomial.polyfit(shifted, y[lo:hi], degree)
        result[index] = coeff[1]
    return result


def coupling_interpretations(rho: np.ndarray, c_values: np.ndarray) -> dict:
    """Return both preregistered inverse-coupling curves and slope estimators.

    The farthest shell is the reference.  ``energy`` means g^2/g_ref^2=C/C_ref;
    ``amplitude`` means g/g_ref=C/C_ref.  Neither is selected here.
    """
    rho = np.asarray(rho, dtype=float)
    c_values = np.asarray(c_values, dtype=float)
    log_mu = -np.log(rho)
    c_ref = float(c_values[-1])
    outputs = {}
    for name, power in (("energy", 1.0), ("amplitude", 2.0)):
        inverse_g2 = (c_ref / c_values) ** power
        outputs[name] = {
            "inverse_g2": inverse_g2,
            "slope_gradient": np.gradient(inverse_g2, log_mu, edge_order=2),
            "slope_local_cubic": local_polynomial_derivative(log_mu, inverse_g2),
        }
    return {"C_ref": c_ref, "log_mu_over_mu0": log_mu, "schemes": outputs}
