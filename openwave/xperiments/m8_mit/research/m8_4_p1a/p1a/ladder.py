"""P1A.4: Resolution ladder over all nine bundles.

One frozen cloud rule, one stencil rule, one resolution sequence.
No sector deciding it deserves more nodes after being looked at.

Records per sector and resolution:
  - cluster position error
  - raw cluster spread
  - principal-angle error (vs free subspace)
  - leakage
  - imaginary eigenvalue contamination
  - eps_SA (self-adjointness residual)
  - spectral separation
  - Riesz cross-check
"""

import numpy as np
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from p0.group import LABELS, DIMS, MCKAY_DIST, build_icosians, build_character_table
from p0.representations import build_all_representations
from p0.cloud import fibonacci_seeds_s3, build_orbit_cloud
from p0.bundle_operator import build_L_bundle

from .mass_matrix import build_Mh_base, build_Mh_rho
from .diagnostics import (run_nonnormality_diagnostics, self_adjointness_residual,
                          riesz_projector_norm)
from .subspace import (extract_and_score, mh_orthonormalize, three_way_identity_check,
                       principal_angles, imaginary_contamination)


RESOLUTION_SEQUENCE = [8, 12, 16, 20, 30, 40, 60]

K_STENCIL_RULE = lambda n_seeds: min(110, max(20, int(n_seeds * 120 * 0.015)))

K_DENSITY_DEFAULT = None


def analytic_riesz_params(label):
    """Analytically determined Riesz contour parameters.

    Centre and radius fixed from known analytic eigenvalue and its
    analytic neighbours BEFORE any numerical eigenvalue is inspected.

    For E_rho at McKay distance d, eigenvalue = d(d+2).
    Next lower level: (d-1)(d-1+2) = d^2-1
    Next higher level: (d+1)(d+1+2) = d^2+4d+3

    Centre at the target eigenvalue, radius = half gap to nearest neighbour.
    """
    idx = LABELS.index(label)
    d = MCKAY_DIST[idx]
    lam = d * (d + 2)

    if d == 0:
        lam_next = 3
        centre = lam
        radius = 1.5
    elif d == 1:
        lam_prev = 0
        lam_next = 8
        gap = min(lam - lam_prev, lam_next - lam)
        centre = lam
        radius = gap * 0.4
    else:
        lam_prev = (d - 1) * (d + 1)
        lam_next = (d + 1) * (d + 3)
        gap = min(lam - lam_prev, lam_next - lam)
        centre = lam
        radius = gap * 0.4

    return float(centre), float(radius)


def imaginary_envelope(im_values_free, n_seeds_list):
    """Derive numerical envelope for imaginary contamination.

    ONE frozen formula. The envelope is max(|Im lambda|) at each resolution,
    fitted as C * n_seeds^(-alpha).
    """
    if len(im_values_free) < 2:
        return {"formula": "insufficient_data"}

    ns = np.array(n_seeds_list, dtype=float)
    ims = np.array(im_values_free, dtype=float)

    mask = ims > 1e-16
    if np.sum(mask) < 2:
        return {"formula": "all_zero", "max_im": float(np.max(ims))}

    log_ns = np.log(ns[mask])
    log_ims = np.log(ims[mask])
    coeffs = np.polyfit(log_ns, log_ims, 1)
    alpha = -coeffs[0]
    C = np.exp(coeffs[1])

    def envelope(n):
        return C * n ** (-alpha)

    return {
        "formula": f"{C:.4e} * n_seeds^(-{alpha:.2f})",
        "C": C, "alpha": alpha,
        "envelope_fn": envelope,
    }


def run_single_resolution(n_seeds, elems, chi, reps, label, W_base, k_stencil):
    """Run one sector at one resolution. Returns a results dict."""
    idx = LABELS.index(label)
    d_rho = MCKAY_DIST[idx]
    d_fiber = DIMS[idx]
    expected_lambda = d_rho * (d_rho + 2)
    expected_dim = d_rho + 1

    seeds = fibonacci_seeds_s3(n_seeds)
    X, oid, gid = build_orbit_cloud(seeds, elems)

    L, seed_orbits = build_L_bundle(X, oid, gid, elems, reps[label],
                                     k=k_stencil)

    Mh_diag = build_Mh_rho(W_base, d_fiber)

    diag_nonnorm = run_nonnormality_diagnostics(L, Mh_diag, label=label)

    eps_sa = self_adjointness_residual(L, Mh_diag)

    centre, radius = analytic_riesz_params(label)

    extract_result, Q = extract_and_score(
        L, Mh_diag, expected_lambda, expected_dim, label=label,
        riesz_centre=centre, riesz_radius=radius)

    result = {
        "n_seeds": n_seeds,
        "k_stencil": k_stencil,
        "label": label,
        "d_rho": d_rho,
        "d_fiber": d_fiber,
        "expected_lambda": expected_lambda,
        "expected_dim": expected_dim,
        "nonnormality": diag_nonnorm,
        "eps_SA": eps_sa,
        "extraction": extract_result,
    }

    return result, Q, L, Mh_diag


def run_ladder(print_fn=print):
    """Run the full resolution ladder over all nine bundles.

    Returns the complete results structure.
    """
    print_fn("Building group and representations...")
    elems = build_icosians()
    chi = build_character_table(elems)
    _, _ = build_all_representations(elems, chi)
    reps_dict, _ = build_all_representations(elems, chi)

    all_results = {}

    for n_seeds in RESOLUTION_SEQUENCE:
        print_fn(f"\n=== Resolution: {n_seeds} seeds ({n_seeds * 120} nodes) ===")

        seeds = fibonacci_seeds_s3(n_seeds)
        X, oid, gid = build_orbit_cloud(seeds, elems)
        W_base = build_Mh_base(X, oid, n_seeds)

        k_stencil = K_STENCIL_RULE(n_seeds)
        print_fn(f"  k_stencil = {k_stencil}")

        pos_check = bool(np.all(W_base > 0))
        print_fn(f"  M_h positivity: {pos_check}")

        for label in LABELS:
            idx = LABELS.index(label)
            d_fiber = DIMS[idx]

            try:
                result, Q, L, Mh_diag = run_single_resolution(
                    n_seeds, elems, chi, reps_dict, label, W_base, k_stencil)

                key = (n_seeds, label)
                all_results[key] = result

                ex = result["extraction"]
                pos = ex.get("cluster_position", "N/A")
                spread = ex.get("cluster_spread", "N/A")
                eps = result["eps_SA"]

                if isinstance(pos, float):
                    print_fn(f"  {label}: pos={pos:.4f}, spread={spread:.4e}, "
                            f"eps_SA={eps:.4e}")
                else:
                    print_fn(f"  {label}: extraction failed")

            except Exception as e:
                print_fn(f"  {label}: ERROR - {e}")
                all_results[(n_seeds, label)] = {"error": str(e)}

    return all_results
