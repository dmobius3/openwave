"""Replacement imaginary-contamination rule.

The original per-sector power-law envelope (C_rho * n^(-alpha_rho)) was
nine independent calibrations.  The pathology: roundoff-dominated sectors
had their noise fitted (exponents from -2.7 to +10.7, prefactors 5e-15
to 4e+10), while genuinely contaminated sectors got plausible exponents
and thus wide envelopes.  Fit quality replaced contamination magnitude as
the effective test.  All nine rulings were withdrawn before any target
execution.

Replacement — two arms:
  Arm A: Numerical floor from eps_mach, ||L||_2, ||P_spec||_2.
         If max|Im lambda| <= floor, the sector is FLOOR_LIMITED.
  Arm B: One pooled convergence law y(h) ~ A h^p across all above-floor
         sectors, with analytic sector normalization s_rho and geometry-
         derived fill distance h.  Margin = exp(r_max) where r_max is the
         worst positive calibration residual.
"""

import numpy as np
from scipy.spatial import cKDTree


EPS_MACH = np.finfo(float).eps


def compute_fill_distance(X):
    """Max nearest-neighbor geodesic distance on S^3 — proxy for fill distance.

    d_geo = 2 arcsin(chord / 2) converts chord to geodesic.
    """
    tree = cKDTree(X)
    dists, _ = tree.query(X, k=2)
    max_chord = float(np.max(dists[:, 1]))
    return 2.0 * np.arcsin(np.clip(max_chord / 2.0, 0.0, 1.0))


def numerical_floor(L_norm_2, P_spec_norm):
    """Arm A: precision floor for imaginary eigenvalue contamination.

    F = eps_mach * ||L||_2 * max(||P_spec||_2, 1)

    Inputs measure the numerical calculation, not the observed Im lambda.
    Below F, imaginary parts are indistinguishable from roundoff in the
    Schur / eigenvalue computation.

    ||P_spec|| = 0 is EXTRACTION_FAIL (Section 4 of LADDER_RULE_TASK):
    a nonzero projector has spectral radius >= 1, so ||P|| = 0 means
    extraction failed catastrophically.  Returns None to signal that
    the floor cannot be computed.
    """
    if P_spec_norm is not None and P_spec_norm < 0.5:
        return None
    amp = max(float(P_spec_norm) if P_spec_norm and P_spec_norm > 0 else 1.0,
              1.0)
    return EPS_MACH * float(L_norm_2) * amp


def analytic_scale(d_rho):
    """s_rho = (1 + d_rho(d_rho + 2)) / R^2, R = 1 for unit S^3."""
    return 1.0 + d_rho * (d_rho + 2)


def pooled_convergence_fit(points):
    """Fit ONE pooled power law log y = log A + p log h.

    points: list of (h, y, tag_str) tuples, all above-floor.
    Returns dict with A, p, r_max, per-point residuals.
    """
    if len(points) < 3:
        return {"status": "insufficient_data", "n_points": len(points)}

    h_arr = np.array([pt[0] for pt in points], dtype=float)
    y_arr = np.array([pt[1] for pt in points], dtype=float)
    tags = [pt[2] for pt in points]

    log_h = np.log(h_arr)
    log_y = np.log(y_arr)

    coeffs = np.polyfit(log_h, log_y, 1)
    p = float(coeffs[0])
    log_A = float(coeffs[1])
    A = float(np.exp(log_A))

    predicted = log_A + p * log_h
    residuals = (log_y - predicted).tolist()
    r_max = float(np.max(residuals))

    return {
        "status": "ok",
        "A": A, "p": p, "r_max": r_max,
        "n_points": len(points),
        "residuals": residuals,
        "tags": tags,
    }


def envelope_value(s_rho, A, p, r_max, h):
    """E_rho(h) = s_rho * A * h^p * exp(r_max)."""
    return s_rho * A * h ** p * np.exp(r_max)


def evaluate_entry(im_max, L_norm_2, P_spec_norm, s_rho, h, fit):
    """Evaluate one (sector, resolution) entry against the frozen rule.

    Returns dict with regime, floor, envelope, eligibility.
    If floor is None (extraction failed), returns EXTRACTION_FAIL.
    """
    floor = numerical_floor(L_norm_2, P_spec_norm)
    if floor is None:
        return {
            "floor": None, "im_max": im_max,
            "L_norm_2": L_norm_2, "P_spec_norm_used": P_spec_norm,
            "s_rho": s_rho, "h": h,
            "regime": "EXTRACTION_FAIL",
            "convergence": "N/A",
            "I_over_E": None,
            "eligibility": "EXTRACTION_FAIL",
        }
    is_floor = im_max <= floor

    result = {
        "floor": floor,
        "im_max": im_max,
        "L_norm_2": L_norm_2,
        "P_spec_norm_used": P_spec_norm,
        "s_rho": s_rho,
        "h": h,
    }

    if is_floor:
        result.update({
            "regime": "FLOOR_LIMITED",
            "convergence": "N/A_FLOOR",
            "I_over_E": None,
            "eligibility": "QUALIFIED",
        })
        return result

    if fit is None or fit.get("status") != "ok":
        result.update({
            "regime": "ABOVE_FLOOR",
            "convergence": "MODEL_FAILED",
            "I_over_E": None,
            "eligibility": "NO_LABEL",
        })
        return result

    A, p, r_max = fit["A"], fit["p"], fit["r_max"]

    if p <= 0:
        result.update({
            "regime": "ABOVE_FLOOR",
            "convergence": "FAIL_p<=0",
            "I_over_E": None,
            "eligibility": "NO_LABEL",
        })
        return result

    E = envelope_value(s_rho, A, p, r_max, h)
    within = im_max <= E

    result.update({
        "regime": "ABOVE_FLOOR",
        "convergence": "PASS" if within else "FAIL",
        "envelope": E,
        "I_over_E": im_max / E if E > 0 else float('inf'),
        "eligibility": "QUALIFIED" if within else "NO_LABEL",
    })
    return result
