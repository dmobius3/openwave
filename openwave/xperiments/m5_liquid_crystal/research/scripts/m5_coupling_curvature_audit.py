"""Independent adversarial audit of the M5/Faber curvature scale scan.

This script deliberately does not import either contribution module.  It reads
the tracked scan JSON as a set of claims and reconstructs the load-bearing
algebra with SymPy, a fourth-order point derivative on fresh grids, and a
CubicSpline derivative estimator.

Run from the repository root:

    python3 openwave/xperiments/m5_liquid_crystal/research/scripts/\
        m5_coupling_curvature_audit.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.interpolate import CubicSpline


HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parent
PRIMARY_JSON = RESEARCH / "data" / "m5_coupling_curvature_scan.json"
OUT_JSON = RESEARCH / "data" / "m5_coupling_curvature_audit.json"


def relative_l2(actual: np.ndarray, expected: np.ndarray) -> float:
    """Return ||actual - expected||_2 / ||expected||_2."""
    return float(np.linalg.norm(actual - expected) / np.linalg.norm(expected))


def symbolic_derivation() -> dict:
    """Re-derive Gamma, |R|^2, and d(log C)/d(log mu) exactly."""
    x, y, z, r0 = sp.symbols("x y z r0", real=True, positive=True)
    rho = sp.symbols("rho", real=True, positive=True)
    coordinates = (x, y, z)
    radius_squared = x**2 + y**2 + z**2
    denominator = sp.sqrt(radius_squared + r0**2)
    q0 = r0 / denominator
    q = sp.Matrix([x, y, z]) / denominator

    gamma = []
    gamma_without_cross = []
    for coordinate in coordinates:
        dq0 = sp.diff(q0, coordinate)
        dq = q.diff(coordinate)
        abelian_part = q0 * dq - dq0 * q
        gamma_without_cross.append(sp.simplify(abelian_part))
        gamma.append(sp.simplify(abelian_part + q.cross(dq)))

    expected_gamma = [
        sp.Matrix([r0, z, -y]),
        sp.Matrix([-z, r0, x]),
        sp.Matrix([y, -x, r0]),
    ]
    expected_gamma = [value / (radius_squared + r0**2) for value in expected_gamma]
    gamma_identity = all(
        all(sp.simplify(component) == 0 for component in (actual - expected))
        for actual, expected in zip(gamma, expected_gamma)
    )

    def curvature_squared(connection: list[sp.Matrix]) -> sp.Expr:
        total = 0
        for i, j in ((0, 1), (0, 2), (1, 2)):
            curvature = connection[i].cross(connection[j])
            total += curvature.dot(curvature)
        return sp.factor(sp.simplify(total))

    curvature_sq = curvature_squared(gamma)
    no_cross_curvature_sq = curvature_squared(gamma_without_cross)
    expected_curvature_sq = (radius_squared + 3 * r0**2) / (
        radius_squared + r0**2
    ) ** 3
    expected_no_cross_sq = 3 * r0**4 / (radius_squared + r0**2) ** 4

    proxy = rho**2 * sp.sqrt(rho**2 + 3) / (rho**2 + 1) ** sp.Rational(3, 2)
    log_slope = sp.factor(-rho * sp.diff(sp.log(proxy), rho))
    expected_log_slope = -6 / ((rho**2 + 1) * (rho**2 + 3))

    return {
        "gamma_identity": bool(gamma_identity),
        "curvature_squared": sp.sstr(curvature_sq),
        "curvature_squared_identity": bool(
            sp.simplify(curvature_sq - expected_curvature_sq) == 0
        ),
        "no_cross_curvature_squared": sp.sstr(no_cross_curvature_sq),
        "no_cross_curvature_squared_identity": bool(
            sp.simplify(no_cross_curvature_sq - expected_no_cross_sq) == 0
        ),
        "dlogC_dlogmu": sp.sstr(log_slope),
        "log_slope_identity": bool(sp.simplify(log_slope - expected_log_slope) == 0),
    }


def hedgehog_values(points: np.ndarray, r0: float) -> tuple[np.ndarray, np.ndarray]:
    """Evaluate q0 and q directly at arbitrary points."""
    denominator = np.sqrt(np.sum(points * points, axis=-1) + r0 * r0)
    return r0 / denominator, points / denominator[..., None]


def fourth_order_derivative(
    points: np.ndarray, axis: int, spacing: float, r0: float
) -> tuple[np.ndarray, np.ndarray]:
    """Differentiate q0 and q with a five-point, fourth-order point stencil."""
    samples = []
    for offset in (-2.0, -1.0, 1.0, 2.0):
        shifted = points.copy()
        shifted[..., axis] += offset * spacing
        samples.append(hedgehog_values(shifted, r0))
    dq0 = (samples[0][0] - 8 * samples[1][0] + 8 * samples[2][0] - samples[3][0]) / (
        12 * spacing
    )
    dq = (samples[0][1] - 8 * samples[1][1] + 8 * samples[2][1] - samples[3][1]) / (
        12 * spacing
    )
    return dq0, dq


def fourth_order_shell_profile(
    *, n: int, half_extent: float, rho: np.ndarray, half_width: float, r0: float = 1.0
) -> dict:
    """Build an independent shell curve with fourth-order point derivatives."""
    coordinates = np.linspace(-half_extent, half_extent, n)
    spacing = float(coordinates[1] - coordinates[0])
    points = np.stack(np.meshgrid(coordinates, coordinates, coordinates, indexing="ij"), -1)
    q0, q = hedgehog_values(points, r0)

    gamma = []
    for axis in range(3):
        dq0, dq = fourth_order_derivative(points, axis, spacing, r0)
        gamma.append(q0[..., None] * dq - dq0[..., None] * q + np.cross(q, dq))

    curvature_squared = np.zeros(points.shape[:-1])
    for i, j in ((0, 1), (0, 2), (1, 2)):
        curvature = np.cross(gamma[i], gamma[j])
        curvature_squared += np.sum(curvature * curvature, axis=-1)

    radius = np.linalg.norm(points, axis=-1)
    numerical_curvature = np.sqrt(curvature_squared)
    exact_curvature = np.sqrt(radius * radius + 3 * r0 * r0) / (
        radius * radius + r0 * r0
    ) ** 1.5
    numerical = []
    exact = []
    counts = []
    for center in rho:
        shell = (radius / r0 >= center - half_width) & (
            radius / r0 < center + half_width
        )
        counts.append(int(np.count_nonzero(shell)))
        numerical.append(float(np.mean(radius[shell] ** 2 * numerical_curvature[shell])))
        exact.append(float(np.mean(radius[shell] ** 2 * exact_curvature[shell])))
    return {
        "n": n,
        "h_over_r0": spacing / r0,
        "C_fourth_order": np.asarray(numerical),
        "C_exact_shell": np.asarray(exact),
        "counts": counts,
    }


def reconstruct_exact_shell(primary_row: dict, rho: np.ndarray, half_width: float) -> dict:
    """Rebuild a stored exact shell without using contribution code."""
    n = int(primary_row["n"])
    extent = float(primary_row["half_extent_over_r0"])
    coordinates = np.linspace(-extent, extent, n)
    x, y, z = np.meshgrid(coordinates, coordinates, coordinates, indexing="ij")
    radius_squared = x * x + y * y + z * z
    radius = np.sqrt(radius_squared)
    proxy = radius_squared * np.sqrt(radius_squared + 3.0) / (radius_squared + 1.0) ** 1.5
    values = []
    counts = []
    for center in rho:
        shell = (radius >= center - half_width) & (radius < center + half_width)
        values.append(float(np.mean(proxy[shell])))
        counts.append(int(np.count_nonzero(shell)))
    return {"C_exact_shell": np.asarray(values), "counts": counts}


def audit_coupling_conventions(primary: dict, rho: np.ndarray) -> dict:
    """Check both conditional convention formulas and slopes independently."""
    point_c = rho**2 * np.sqrt(rho**2 + 3.0) / (rho**2 + 1.0) ** 1.5
    point_c_ref = float(point_c[-1])
    log_mu = -np.log(rho)
    denominator = (rho**2 + 1.0) * (rho**2 + 3.0)
    schemes = {}
    for name, power in (("energy", 1.0), ("amplitude", 2.0)):
        inverse_g2 = (point_c_ref / point_c) ** power
        exact_slope = 6.0 * power * inverse_g2 / denominator
        spline = CubicSpline(log_mu[::-1], inverse_g2[::-1])
        spline_slope = spline(log_mu, 1)
        stored = primary["analytic_point_interpretations"]["schemes"][name]
        stored_inverse = np.asarray(stored["inverse_g2"], dtype=float)
        stored_local_slope = np.asarray(stored["slope_local_cubic"], dtype=float)
        schemes[name] = {
            "formula_power": power,
            "stored_inverse_g2_max_abs_error": float(
                np.max(np.abs(stored_inverse - inverse_g2))
            ),
            "cubic_spline_vs_exact_slope_rel_l2_all": relative_l2(
                spline_slope, exact_slope
            ),
            "cubic_spline_vs_exact_slope_rel_l2_interior": relative_l2(
                spline_slope[2:-2], exact_slope[2:-2]
            ),
            "stored_local_cubic_vs_exact_slope_rel_l2_interior": relative_l2(
                stored_local_slope[2:-2], exact_slope[2:-2]
            ),
            "exact_slope_first_last": [float(exact_slope[0]), float(exact_slope[-1])],
            "spline_slope_first_last": [float(spline_slope[0]), float(spline_slope[-1])],
            "wrong_mu_equals_r_signs_are_negative": bool(np.all(-exact_slope < 0.0)),
        }

    measured_endpoint_errors = {}
    for name in ("energy", "amplitude"):
        measured = np.asarray(
            primary["measured_interpretations"]["schemes"][name]["slope_local_cubic"]
        )
        exact_shell = np.asarray(
            primary["analytic_shell_interpretations"]["schemes"][name][
                "slope_local_cubic"
            ]
        )
        measured_endpoint_errors[name] = (
            np.abs(measured[[0, -1]] / exact_shell[[0, -1]] - 1.0).tolist()
        )

    return {
        "units": {
            "q": "dimensionless",
            "Gamma": "length^-1",
            "R": "length^-2",
            "C": "dimensionless",
            "d_inverse_g2_d_log_mu": "dimensionless, conditional on a dimensionless g_ref",
        },
        "schemes": schemes,
        "measured_endpoint_relative_errors_vs_same_shell_estimator": measured_endpoint_errors,
    }


def audit_mutations(primary: dict, rho: np.ndarray) -> dict:
    """Attack the cross term, scale direction, and normalization sensitivity."""
    far = rho >= 2.5
    full_curvature = np.sqrt(rho**2 + 3.0) / (rho**2 + 1.0) ** 1.5
    no_cross_curvature = np.sqrt(3.0) / (rho**2 + 1.0) ** 2
    full_exponent = float(np.polyfit(np.log(rho[far]), np.log(full_curvature[far]), 1)[0])
    no_cross_exponent = float(
        np.polyfit(np.log(rho[far]), np.log(no_cross_curvature[far]), 1)[0]
    )
    stored_mutation = primary["mutations"]["remove_q_cross_dq"]

    scales = (0.7, 1.0, 1.9)
    dimensionless_curves = []
    for r0 in scales:
        radius = r0 * rho
        curvature = np.sqrt(radius**2 + 3 * r0**2) / (radius**2 + r0**2) ** 1.5
        dimensionless_curves.append(radius**2 * curvature)
    covariance_residual = max(
        relative_l2(curve, dimensionless_curves[1]) for curve in dimensionless_curves
    )

    lambda_gamma = 1.2
    return {
        "full_exact_far_exponent": full_exponent,
        "no_cross_exact_far_exponent": no_cross_exponent,
        "exponent_change": no_cross_exponent - full_exponent,
        "stored_full_exponent_abs_error": abs(
            full_exponent - stored_mutation["full_connection_R_falloff_exponent"]
        ),
        "stored_no_cross_exponent_abs_error": abs(
            no_cross_exponent - stored_mutation["no_cross_R_falloff_exponent"]
        ),
        "scale_covariance_rel_l2": covariance_residual,
        "connection_scale_identity": {
            "lambda": lambda_gamma,
            "raw_C_ratio": lambda_gamma**2,
            "far_normalized_shape_residual": 0.0,
            "audit_assessment": (
                "algebraically true but tautological; it is not an independent shape mutation"
            ),
        },
    }


def run() -> dict:
    primary = json.loads(PRIMARY_JSON.read_text())
    rho = np.asarray(primary["rho"], dtype=float)
    half_width = 0.16

    symbolic = symbolic_derivation()

    exact_shell_reconstruction = []
    recomputed_primary_errors = []
    for row in primary["spatial_refinement"]:
        reconstructed = reconstruct_exact_shell(row, rho, half_width)
        stored_exact = np.asarray(row["C_exact_shell"], dtype=float)
        stored_field = np.asarray(row["C"], dtype=float)
        exact_shell_reconstruction.append(
            {
                "n": row["n"],
                "max_abs_C_exact_shell_error": float(
                    np.max(np.abs(reconstructed["C_exact_shell"] - stored_exact))
                ),
                "counts_match": reconstructed["counts"] == row["counts"],
            }
        )
        recomputed_primary_errors.append(relative_l2(stored_field, stored_exact))

    fourth_order_rows = []
    for n in (41, 61):
        independent = fourth_order_shell_profile(
            n=n, half_extent=8.0, rho=rho, half_width=half_width
        )
        primary_row = next(row for row in primary["spatial_refinement"] if row["n"] == n)
        primary_c = np.asarray(primary_row["C"], dtype=float)
        fourth_order_rows.append(
            {
                "n": n,
                "h_over_r0": independent["h_over_r0"],
                "fourth_order_vs_exact_shell_rel_l2": relative_l2(
                    independent["C_fourth_order"], independent["C_exact_shell"]
                ),
                "primary_second_order_vs_fourth_order_rel_l2": relative_l2(
                    primary_c, independent["C_fourth_order"]
                ),
                "minimum_shell_count": min(independent["counts"]),
            }
        )

    conventions = audit_coupling_conventions(primary, rho)
    mutations = audit_mutations(primary, rho)

    stored_spatial_errors = np.asarray(primary["spatial_rel_l2_errors"], dtype=float)
    recomputed_spatial_errors = np.asarray(recomputed_primary_errors)
    gates = {
        "symbolic_gamma_identity": symbolic["gamma_identity"],
        "symbolic_curvature_identity": symbolic["curvature_squared_identity"],
        "symbolic_no_cross_identity": symbolic["no_cross_curvature_squared_identity"],
        "symbolic_log_slope_identity": symbolic["log_slope_identity"],
        "stored_exact_shells_reconstruct_to_1e_minus_14": max(
            row["max_abs_C_exact_shell_error"] for row in exact_shell_reconstruction
        )
        < 1e-14,
        "stored_shell_counts_reconstruct": all(
            row["counts_match"] for row in exact_shell_reconstruction
        ),
        "stored_spatial_errors_recompute_to_1e_minus_14": bool(
            np.max(np.abs(stored_spatial_errors - recomputed_spatial_errors)) < 1e-14
        ),
        "independent_fourth_order_converges": (
            fourth_order_rows[1]["fourth_order_vs_exact_shell_rel_l2"]
            < fourth_order_rows[0]["fourth_order_vs_exact_shell_rel_l2"]
        ),
        "independent_n61_fourth_order_vs_exact_below_0p1pct": (
            fourth_order_rows[1]["fourth_order_vs_exact_shell_rel_l2"] < 0.001
        ),
        "primary_n61_agrees_with_independent_route_below_2pct": (
            fourth_order_rows[1]["primary_second_order_vs_fourth_order_rel_l2"] < 0.02
        ),
        "both_inverse_coupling_formulas_match": all(
            values["stored_inverse_g2_max_abs_error"] < 1e-14
            for values in conventions["schemes"].values()
        ),
        "independent_spline_slopes_match_exact_interior_below_0p1pct": all(
            values["cubic_spline_vs_exact_slope_rel_l2_interior"] < 0.001
            for values in conventions["schemes"].values()
        ),
        "wrong_mu_direction_flips_slope_sign": all(
            values["wrong_mu_equals_r_signs_are_negative"]
            for values in conventions["schemes"].values()
        ),
        "cross_term_mutation_changes_far_exponent_by_more_than_1": abs(
            mutations["exponent_change"]
        )
        > 1.0,
        "stored_mutation_exponents_match_exact_route_below_0p01": (
            mutations["stored_full_exponent_abs_error"] < 0.01
            and mutations["stored_no_cross_exponent_abs_error"] < 0.01
        ),
        "dimensionless_profile_is_scale_covariant": mutations["scale_covariance_rel_l2"]
        < 1e-14,
    }

    verdicts = [
        {
            "claim": "Displayed hedgehog implies the stated Gamma, curvature norm, C(rho), and logarithmic slope.",
            "verdict": "CONFIRMED",
            "reason": "Independent exact SymPy differentiation proves all four identities.",
        },
        {
            "claim": "The tracked second-order shell curve resolves the continuum ansatz profile and converges spatially.",
            "verdict": "CONFIRMED",
            "reason": (
                "All stored shell oracles and errors reconstruct, while an independent fourth-order grid route "
                "converges and agrees with the n=61 stored curve within 2%."
            ),
        },
        {
            "claim": "The fixed-h box ladder is domain-refinement evidence.",
            "verdict": "PARTIAL",
            "reason": (
                "The reported near-zero spread is real, but all three odd grids share the same h and identical "
                "sample points on the measured shells; this local analytic stencil has no boundary solve, so the "
                "check is box invariance, not a nontrivial domain-convergence test."
            ),
        },
        {
            "claim": "The energy/action and field-amplitude formulas are the two stated conditional inverse-coupling readings.",
            "verdict": "CONFIRMED",
            "reason": (
                "Independent algebra and CubicSpline differentiation recover both powers, positive slopes for "
                "mu/mu0=1/rho, and the exact pointwise derivatives."
            ),
        },
        {
            "claim": "The reported measured slope endpoints are validated by the frozen derivative gates.",
            "verdict": "PARTIAL",
            "reason": (
                "The gates deliberately exclude two samples at each end.  The far-end measured local-cubic "
                "slopes differ from the same-shell estimator by about 9.6%, so endpoints should be labeled "
                "diagnostic or the quoted range should use the validated interior."
            ),
        },
        {
            "claim": "Removing q cross dq destroys the Coulomb-like far-field plateau.",
            "verdict": "CONFIRMED",
            "reason": (
                "The independent exact mutation gives |R|=sqrt(3) r0^2/(r^2+r0^2)^2 and reproduces the "
                "stored exponent change from approximately -1.97 to -3.73."
            ),
        },
        {
            "claim": "Gamma rescaling is an independent load-bearing mutation.",
            "verdict": "PARTIAL",
            "reason": (
                "R scales as lambda^2 and far normalization cancels it identically; the numbers are correct, "
                "but this is a dimensional identity rather than an independent sensitivity test."
            ),
        },
        {
            "claim": "This scan measures an M5 renormalized running coupling or beta-function coefficient.",
            "verdict": "REFUTED",
            "reason": (
                "Only a classical single-ansatz curvature form factor is fixed.  The source/action dictionary, "
                "stationary two-core scheme, field map, b0, and flavour fit are all absent; the note's explicit "
                "not-computed boundary correctly acknowledges this."
            ),
        },
        {
            "claim": "The note's explicit not-computed list keeps the contribution within instrument scope.",
            "verdict": "CONFIRMED",
            "reason": (
                "It does not select a coupling convention or claim b0, n_f, QCD matching, or a two-core force."
            ),
        },
    ]

    return {
        "status": "independent adversarial audit complete; mathematical core confirmed with scope qualifications",
        "independence": {
            "primary_modules_imported": False,
            "primary_json_treated_as_claims": str(PRIMARY_JSON.relative_to(RESEARCH)),
            "routes": [
                "exact SymPy differentiation",
                "five-point fourth-order point derivatives on independently constructed grids",
                "SciPy CubicSpline differentiation",
                "exact analytic load-bearing mutations",
            ],
        },
        "symbolic_derivation": symbolic,
        "exact_shell_reconstruction": exact_shell_reconstruction,
        "recomputed_primary_spatial_rel_l2_errors": recomputed_primary_errors,
        "independent_fourth_order_profiles": fourth_order_rows,
        "coupling_conventions": conventions,
        "mutations": mutations,
        "scope_fixes_needed": [
            "Call the fixed-h box ladder a box-invariance check, not substantive domain refinement.",
            "Label the two endpoint slope values diagnostic or quote the interior range covered by the derivative gates.",
            "Call gradient/local-cubic two estimators, not independent derivations; the audit CubicSpline is the independent estimator.",
            "Keep 'running coupling' conditional: the deliverable is a classical ansatz form factor until an M5 source/action or stationary two-core dictionary selects g_R.",
        ],
        "verdicts": verdicts,
        "gates": gates,
        "all_mathematical_gates_pass": all(gates.values()),
    }


def main() -> int:
    result = run()
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"gates": result["gates"]}, indent=2))
    print(f"wrote {OUT_JSON}")
    return 0 if result["all_mathematical_gates_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
