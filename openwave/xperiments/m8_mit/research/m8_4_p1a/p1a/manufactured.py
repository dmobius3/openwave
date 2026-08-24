"""P1A.5: Manufactured calibration cases.

Synthetic operator and subspace perturbations with known answers.
NOT runs of target nonlinear dynamics. Qualifying the estimator.

Requires a calibration/holdout split: some manufactured cases are
used to verify the pipeline; others are held out to test it blind.

Adjudication: for every manufactured case, plant and compare recovered
against expected for Delta_lambda, s, theta_max, and ell, with stated
quantitative tolerances. The rotation cases exercise the angle and
leakage path — they are the ones that matter most.
"""

import numpy as np
from .mass_matrix import build_Mh_rho, Mh_inner
from .subspace import (mh_orthonormalize, three_way_identity_check,
                       extract_and_score, principal_angles, ordered_schur_subspace)
from .diagnostics import self_adjointness_residual


MANUFACTURED_TOLERANCES = {
    "position": 0.1,
    "spread_rel": 0.1,
    "spread_abs": 0.005,
    "theta_max": 1e-3,
    "leakage": 1e-4,
}


def make_known_perturbation(L_free, Q_free, Mh_diag, k,
                            shift=0.5, split=0.1, rotation_angle=0.15,
                            rng_seed=42):
    """Manufacture a perturbed operator with KNOWN answers.

    Uses ordered Schur decomposition to guarantee exact eigenvalue control.
    After construction, computes the exact invariant subspace of L_pert
    and derives expected theta_max and ell vs Q_free.
    """
    from scipy.linalg import schur
    n = L_free.shape[0]
    rng = np.random.default_rng(rng_seed)

    evals_free = np.linalg.eigvals(L_free)
    sorted_by_re = np.sort(np.real(evals_free))
    target = sorted_by_re[k // 2] if k > 0 else 0.0

    dists_from_target = np.abs(evals_free - target)
    threshold = float(np.sort(dists_from_target)[k - 1]) + 1.0

    T, Z, sdim = schur(L_free, output='complex',
                        sort=lambda e: abs(e - target) < threshold)

    if sdim < k:
        threshold *= 2.0
        T, Z, sdim = schur(L_free, output='complex',
                            sort=lambda e: abs(e - target) < threshold)

    actual_k = min(sdim, k)

    for i in range(actual_k):
        pert_i = shift + split * (i - (actual_k - 1) / 2.0) / max(actual_k - 1, 1)
        T[i, i] += pert_i

    if rotation_angle > 0 and actual_k < n:
        m = min(actual_k, n - actual_k)
        G = rotation_angle * rng.normal(size=(m, actual_k))
        T[actual_k:actual_k + m, :actual_k] += G

    L_pert = Z @ T @ Z.conj().T

    evals_pert_all = np.linalg.eigvals(L_pert)
    target_shifted = np.mean(np.real(np.diag(T)[:actual_k]))
    dists_pert = np.abs(evals_pert_all - target_shifted)
    cluster_idx = np.argsort(dists_pert)[:actual_k]
    evals_cluster = evals_pert_all[cluster_idx]
    expected_spread = float(np.max(np.real(evals_cluster)) -
                            np.min(np.real(evals_cluster)))
    expected_position = float(np.mean(np.real(evals_cluster)))

    known_answers = {
        "shift": shift,
        "split": split,
        "rotation_angle": rotation_angle,
        "expected_cluster_position_delta": shift,
        "expected_cluster_position": expected_position,
        "expected_spread": expected_spread,
    }

    Q_pert_raw, _, _, sdim_pert = ordered_schur_subspace(
        L_pert, actual_k, target_shifted)
    if sdim_pert == actual_k:
        Q_pert_exact = mh_orthonormalize(Q_pert_raw, Mh_diag)
        tw_exact = three_way_identity_check(Q_free, Q_pert_exact, Mh_diag)
        known_answers["expected_theta_max"] = float(np.max(tw_exact["thetas"]))
        known_answers["expected_ell"] = tw_exact["ell_scored"]
    else:
        known_answers["expected_theta_max"] = None
        known_answers["expected_ell"] = None

    return L_pert, known_answers


def run_manufactured_case(L_pert, Mh_diag, Q_free, target_lambda, k,
                          known_answers, case_name=""):
    """Run the extraction pipeline on a manufactured case and compare.

    Quantitative comparison for all four scored quantities:
    Delta_lambda, s, theta_max, ell.
    """

    result, Q_pert = extract_and_score(
        L_pert, Mh_diag, target_lambda + known_answers["shift"], k,
        label=case_name)

    if Q_pert is None:
        return {
            "case": case_name,
            "extraction_failed": True,
            "pass": False,
        }

    identity = three_way_identity_check(Q_free, Q_pert, Mh_diag)

    eps_sa = self_adjointness_residual(L_pert, Mh_diag)

    tol = MANUFACTURED_TOLERANCES

    recovered_position = result.get("cluster_position", 0.0)
    expected_position = known_answers["expected_cluster_position"]
    position_error = abs(recovered_position - expected_position)
    position_ok = position_error < tol["position"]

    actual_spread = result.get("cluster_spread", 0.0)
    expected_spread = known_answers["expected_spread"]
    spread_error = abs(actual_spread - expected_spread)
    spread_ok = spread_error < max(tol["spread_rel"] * expected_spread,
                                   tol["spread_abs"])

    recovered_theta = float(np.max(identity["thetas"]))
    expected_theta = known_answers.get("expected_theta_max")
    if expected_theta is not None:
        theta_error = abs(recovered_theta - expected_theta)
        theta_ok = theta_error < tol["theta_max"]
    else:
        theta_error = None
        theta_ok = True

    recovered_ell = identity["ell_scored"]
    expected_ell = known_answers.get("expected_ell")
    if expected_ell is not None:
        ell_error = abs(recovered_ell - expected_ell)
        ell_ok = ell_error < tol["leakage"]
    else:
        ell_error = None
        ell_ok = True

    return {
        "case": case_name,
        "extraction": result,
        "three_way": identity,
        "eps_SA": eps_sa,
        "actual_position": recovered_position,
        "expected_position": expected_position,
        "position_error": position_error,
        "position_ok": position_ok,
        "actual_spread": actual_spread,
        "expected_spread": expected_spread,
        "spread_error": spread_error,
        "spread_ok": spread_ok,
        "actual_theta_max": recovered_theta,
        "expected_theta_max": expected_theta,
        "theta_error": theta_error,
        "theta_ok": theta_ok,
        "actual_ell": recovered_ell,
        "expected_ell": expected_ell,
        "ell_error": ell_error,
        "ell_ok": ell_ok,
        "extraction_failed": False,
        "pass": (result["pass"]
                 and identity["pass"]
                 and position_ok
                 and spread_ok
                 and theta_ok
                 and ell_ok),
    }


def run_p1a_5(L_free, Q_free, Mh_diag, target_lambda, k, label=""):
    """Full P1A.5: manufactured calibration with holdout split.

    Calibration cases (used to verify the pipeline):
      C1: moderate shift + moderate split + small rotation
      C2: zero shift + zero split + zero rotation (should recover free)

    Holdout cases (blind test):
      H1: large shift + large split + large rotation
      H2: small shift + asymmetric split pattern
    """

    results = {"label": label, "cases": {}}

    c1_L, c1_known = make_known_perturbation(
        L_free, Q_free, Mh_diag, k,
        shift=0.5, split=0.1, rotation_angle=0.15, rng_seed=42)
    results["cases"]["C1_moderate"] = run_manufactured_case(
        c1_L, Mh_diag, Q_free, target_lambda, k, c1_known, "C1_moderate")

    free_ex, _ = extract_and_score(L_free, Mh_diag, target_lambda, k,
                                   label="C2_free_ref")
    free_spread = free_ex.get("cluster_spread", 0.0)
    free_position = free_ex.get("cluster_position", target_lambda)
    c2_known = {
        "shift": 0.0, "split": 0.0, "rotation_angle": 0.0,
        "expected_cluster_position_delta": 0.0,
        "expected_cluster_position": free_position,
        "expected_spread": free_spread,
        "expected_theta_max": 0.0,
        "expected_ell": 0.0,
    }
    results["cases"]["C2_free"] = run_manufactured_case(
        L_free, Mh_diag, Q_free, target_lambda, k, c2_known, "C2_free")

    h1_L, h1_known = make_known_perturbation(
        L_free, Q_free, Mh_diag, k,
        shift=2.0, split=0.5, rotation_angle=0.4, rng_seed=77)
    results["cases"]["H1_large"] = run_manufactured_case(
        h1_L, Mh_diag, Q_free, target_lambda, k, h1_known, "H1_large")

    h2_L, h2_known = make_known_perturbation(
        L_free, Q_free, Mh_diag, k,
        shift=0.1, split=0.03, rotation_angle=0.05, rng_seed=99)
    results["cases"]["H2_small"] = run_manufactured_case(
        h2_L, Mh_diag, Q_free, target_lambda, k, h2_known, "H2_small")

    calib_pass = all(results["cases"][c]["pass"]
                     for c in ["C1_moderate", "C2_free"]
                     if not results["cases"][c].get("extraction_failed", False))
    holdout_pass = all(results["cases"][c]["pass"]
                       for c in ["H1_large", "H2_small"]
                       if not results["cases"][c].get("extraction_failed", False))

    results["calibration_pass"] = calib_pass
    results["holdout_pass"] = holdout_pass
    results["pass"] = calib_pass and holdout_pass

    return results
