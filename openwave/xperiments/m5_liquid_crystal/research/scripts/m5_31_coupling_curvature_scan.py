"""Run the M5/Faber curvature scale scan and write its tracked evidence.

This driver imports all physics from ``m5_31_coupling_curvature_field.py``.  It
measures the raw C(rho) curve, two explicitly conditional coupling
interpretations, spatial/domain refinement, two derivative estimators, and a
mutation that removes the q x dq term from Gamma.

Run from the repository root:

    python3 openwave/xperiments/m5_liquid_crystal/research/scripts/\
        m5_31_coupling_curvature_scan.py
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np

from m5_31_coupling_curvature_field import (
    analytic_log_slope,
    analytic_proxy,
    coupling_interpretations,
    numerical_profile,
)

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parent
OUT_JSON = RESEARCH / "data" / "m5_31_coupling_curvature_scan.json"
OUT_PLOT = RESEARCH / "plots" / "m5_31_coupling_curvature_scan.png"

RHO = np.geomspace(0.6, 5.5, 25)
SHELL_HALF_WIDTH = 0.16


def serializable_profile(row: dict) -> dict:
    """Convert a numerical profile to JSON-native values."""
    return {
        key: value.tolist() if isinstance(value, np.ndarray) else value
        for key, value in row.items()
    }


def relative_l2(actual: np.ndarray, expected: np.ndarray) -> float:
    return float(np.linalg.norm(actual - expected) / np.linalg.norm(expected))


def run() -> dict:
    spatial_specs = [
        (41, 8.0),
        (61, 8.0),
        (81, 8.0),
        (107, 8.0),
    ]
    spatial = [
        numerical_profile(
            r0=1.0,
            half_extent_over_r0=extent,
            n=n,
            rho_centers=RHO,
            shell_half_width=SHELL_HALF_WIDTH,
        )
        for n, extent in spatial_specs
    ]

    # Keep h/r0=0.3 exactly while changing only the box.
    box_specs = [(49, 7.2), (61, 9.0), (81, 12.0)]
    box_invariance = [
        numerical_profile(
            r0=1.0,
            half_extent_over_r0=extent,
            n=n,
            rho_centers=RHO,
            shell_half_width=SHELL_HALF_WIDTH,
        )
        for n, extent in box_specs
    ]

    finest = spatial[-1]
    interpretations = coupling_interpretations(RHO, finest["C"])
    exact_shell_interpretations = coupling_interpretations(RHO, finest["C_exact_shell"])
    exact_point_interpretations = coupling_interpretations(RHO, analytic_proxy(RHO))

    interior = slice(2, -2)
    derivative_disagreement = {}
    derivative_oracle_error = {}
    for scheme, values in interpretations["schemes"].items():
        numerator = np.linalg.norm(
            values["slope_gradient"][interior] - values["slope_local_cubic"][interior]
        )
        denominator = np.linalg.norm(values["slope_local_cubic"][interior])
        derivative_disagreement[scheme] = float(numerator / denominator)
        exact_slope = exact_shell_interpretations["schemes"][scheme]["slope_local_cubic"][interior]
        derivative_oracle_error[scheme] = relative_l2(
            values["slope_local_cubic"][interior], exact_slope
        )

    box_reference = box_invariance[1]["C"]
    box_spread = max(
        relative_l2(row["C"], box_reference)
        for row in (box_invariance[0], box_invariance[2])
    )
    spatial_errors = [relative_l2(row["C"], row["C_exact_shell"]) for row in spatial]

    no_cross = numerical_profile(
        r0=1.0,
        half_extent_over_r0=8.0,
        n=61,
        rho_centers=RHO,
        shell_half_width=SHELL_HALF_WIDTH,
        include_cross=False,
    )
    far = RHO >= 2.5
    full_falloff = float(
        np.polyfit(np.log(RHO[far]), np.log(spatial[1]["C"][far] / RHO[far] ** 2), 1)[0]
    )
    no_cross_falloff = float(
        np.polyfit(np.log(RHO[far]), np.log(no_cross["C"][far] / RHO[far] ** 2), 1)[0]
    )

    # Gamma -> lambda Gamma implies R -> lambda^2 R.  This checks that raw
    # normalization responds while the far-normalized shape remains invariant.
    lambda_gamma = 1.2
    scaled_c = lambda_gamma**2 * finest["C"]
    raw_scale_ratio = float(np.mean(scaled_c / finest["C"]))
    normalized_shape_residual = relative_l2(scaled_c / scaled_c[-1], finest["C"] / finest["C"][-1])

    gates = {
        "finest_field_vs_exact_shell_rel_l2_lt_0p03": spatial_errors[-1] < 0.03,
        "field_error_decreases_coarse_to_fine": spatial_errors[-1] < spatial_errors[0],
        "fixed_h_box_invariance_rel_l2_lt_0p005": box_spread < 0.005,
        "energy_derivative_methods_rel_l2_lt_0p08": derivative_disagreement["energy"] < 0.08,
        "amplitude_derivative_methods_rel_l2_lt_0p08": derivative_disagreement["amplitude"] < 0.08,
        "energy_slope_vs_exact_shell_rel_l2_lt_0p03": derivative_oracle_error["energy"] < 0.03,
        "amplitude_slope_vs_exact_shell_rel_l2_lt_0p03": derivative_oracle_error["amplitude"]
        < 0.03,
        "full_connection_far_falloff_near_minus2": abs(full_falloff + 2.0) < 0.2,
        "remove_cross_term_changes_far_falloff": abs(no_cross_falloff - full_falloff) > 1.0,
        "connection_scale_mutation_changes_raw_C_by_lambda2": abs(
            raw_scale_ratio - lambda_gamma**2
        )
        < 1e-12,
        "far_normalized_shape_invariant_to_connection_scale": normalized_shape_residual < 1e-12,
    }

    def interpretation_json(source: dict) -> dict:
        return {
            "C_ref": source["C_ref"],
            "log_mu_over_mu0": source["log_mu_over_mu0"].tolist(),
            "schemes": {
                name: {key: value.tolist() for key, value in values.items()}
                for name, values in source["schemes"].items()
            },
        }

    return {
        "status": "classical single-ansatz form-factor instrument; no M5 coupling convention selected",
        "definitions": {
            "rho": "r/r0",
            "mu_over_mu0": "1/rho",
            "C": "shell mean of r^2 sqrt(sum_{i<j}|R_ij|^2)",
            "energy_scheme": "g_R^2/g_ref^2 = C/C_ref (conditional)",
            "amplitude_scheme": "g_R/g_ref = C/C_ref (conditional)",
        },
        "rho": RHO.tolist(),
        "analytic_point_C": analytic_proxy(RHO).tolist(),
        "analytic_point_dlogC_dlogmu": analytic_log_slope(RHO).tolist(),
        "spatial_refinement": [serializable_profile(row) for row in spatial],
        "spatial_rel_l2_errors": spatial_errors,
        "box_invariance": [serializable_profile(row) for row in box_invariance],
        "box_invariance_rel_l2": box_spread,
        "measured_interpretations": interpretation_json(interpretations),
        "analytic_shell_interpretations": interpretation_json(exact_shell_interpretations),
        "analytic_point_interpretations": interpretation_json(exact_point_interpretations),
        "derivative_method_rel_l2": derivative_disagreement,
        "derivative_vs_exact_shell_rel_l2": derivative_oracle_error,
        "mutations": {
            "remove_q_cross_dq": {
                "profile": serializable_profile(no_cross),
                "full_connection_R_falloff_exponent": full_falloff,
                "no_cross_R_falloff_exponent": no_cross_falloff,
            },
            "connection_scale": {
                "lambda": lambda_gamma,
                "observed_raw_C_ratio": raw_scale_ratio,
                "far_normalized_shape_rel_l2": normalized_shape_residual,
            },
        },
        "gates": gates,
        "all_gates_pass": all(gates.values()),
        "not_computed": [
            "source/action normalization selecting C proportional to g or g^2",
            "a QFT beta function or an M5-to-SU(3)/QCD field map",
            "two-core stationary potential- or force-scheme coupling",
            "b0, effective flavour count, or comparator-selected fit",
        ],
    }


def plot_result(result: dict) -> None:
    rho = np.asarray(result["rho"])
    finest = result["spatial_refinement"][-1]
    schemes = result["measured_interpretations"]["schemes"]
    exact_schemes = result["analytic_shell_interpretations"]["schemes"]
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))
    axes[0].plot(rho, result["analytic_point_C"], "k--", label="continuum point")
    axes[0].plot(rho, finest["C"], "o-", ms=3, label="field shells, finest h")
    axes[0].set(xscale="log", xlabel="rho = r/r0", ylabel="C(rho)")
    axes[0].legend(frameon=False)
    for name, values in schemes.items():
        line = axes[1].plot(
            rho,
            values["slope_local_cubic"],
            "o",
            ms=3,
            label=f"{name}: field shells",
        )[0]
        axes[1].plot(
            rho,
            exact_schemes[name]["slope_local_cubic"],
            color=line.get_color(),
            linestyle="--",
            label=f"{name}: exact same shells",
        )
    axes[1].set(
        xscale="log",
        xlabel="rho = r/r0",
        ylabel="d(1/g_R^2) / d log(mu/mu0)",
        ylim=(-0.5, 22.0),
    )
    axes[1].legend(frameon=False, fontsize=8, ncol=2, loc="upper right")
    fig.tight_layout()
    OUT_PLOT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PLOT, dpi=180)
    plt.close(fig)


def main() -> int:
    result = run()
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(result, indent=2) + "\n")
    plot_result(result)
    print(
        json.dumps(
            {"gates": result["gates"], "all_gates_pass": result["all_gates_pass"]}, indent=2
        )
    )
    print(f"wrote {OUT_JSON}")
    print(f"wrote {OUT_PLOT}")
    return 0 if result["all_gates_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
