"""P1A.4a: Cloud and operator admissibility gate.

Qualifies the resolution family as a one-parameter refinement family
before any contamination refit.  All gate checks use geometry-only
quantities or the R0 (scalar/trivial) bundle operator.  No sector-
specific spectral information enters the gate.

Measurement definitions (frozen):
  h  = max nearest-neighbour geodesic distance on the full S^3 orbit cloud
  q  = (1/2) * min pairwise geodesic distance on the full orbit cloud
  geodesic distance: d = 2 arcsin(chord / 2)  where chord is Euclidean R^4
  Fill distance and separation are measured on the FULL 120n-point cloud,
  not on quotient representatives.

The stencil parameter k = 110 is the P0-qualified value and is fixed
across the entire ladder.  K_STENCIL_RULE (which varied k with n) changed
the discretization SCHEME along the ladder; it is excluded from P1A.4a.
"""

import numpy as np
from scipy.spatial import cKDTree

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from p0.cloud import fibonacci_seeds_s3, build_orbit_cloud
from p0.bundle_operator import build_L_bundle, orbit_stencils
from p0.rbffd import rbf_row, monomials, eval_monos, surf_lap_monos


# ── Frozen parameters ──────────────────────────────────────────────────
K_FIXED = 110
RBF_M = 7
RBF_P = 4

CALIBRATION_CANDIDATES = [
    8, 12, 16, 20, 24, 28, 30, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72
]
HOLDOUT_FIXED = [80, 120]
FINER_HOLDOUT_SEARCH = [140, 160, 180, 200, 220, 240]
FINER_HOLDOUT_C = 0.90

GATE = {
    "mesh_ratio_max": 20.0,
    "stencil_cond_max": 1e15,
    "const_residual_max": 1e-6,
    "reprod_max": 0.1,
}

N_STABILITY_REF = 10000


# ── Geometry ───────────────────────────────────────────────────────────

def compute_fill_distance(X):
    """h = max nearest-neighbour geodesic on the full orbit cloud."""
    tree = cKDTree(X)
    dists, _ = tree.query(X, k=2)
    max_chord = float(np.max(dists[:, 1]))
    return 2.0 * np.arcsin(np.clip(max_chord / 2.0, 0.0, 1.0))


def compute_separation_radius(X):
    """q = (1/2) min pairwise geodesic distance."""
    tree = cKDTree(X)
    dists, _ = tree.query(X, k=2)
    min_chord = float(np.min(dists[:, 1]))
    d_geo = 2.0 * np.arcsin(np.clip(min_chord / 2.0, 0.0, 1.0))
    return d_geo / 2.0


def fill_distance_reference(X, n_ref=None):
    """Fill distance from a dense independent reference set on S^3.

    Uses a Fibonacci lattice of n_ref points on S^3 (no orbit expansion)
    as a reference set.  The fill distance is max over reference points of
    min distance to the cloud.  This is an upper bound on the true fill
    distance, whereas the cloud-based measure is a lower bound.
    """
    if n_ref is None:
        n_ref = N_STABILITY_REF
    ref_seeds = fibonacci_seeds_s3(n_ref)
    ref = np.array(ref_seeds)
    tree = cKDTree(X)
    dists, _ = tree.query(ref, k=1)
    max_chord = float(np.max(dists))
    return 2.0 * np.arcsin(np.clip(max_chord / 2.0, 0.0, 1.0))


# ── Stencil quality ───────────────────────────────────────────────────

def stencil_conditioning_stats(X, oid, gid, elems,
                                k=K_FIXED, m=RBF_M, p=RBF_P):
    """Effective condition numbers of the augmented RBF-FD system.

    One value per orbit representative (equivariant copies share geometry).
    The augmented system uses extrinsic R^4 polynomials on S^3, so
    polynomial dependencies (x1^2+x2^2+x3^2+x4^2=1 etc.) create
    ~15 inherently near-zero singular values.  The effective condition
    number excludes those rank-deficient directions.
    """
    node_of, mult, plan = orbit_stencils(X, oid, gid, elems, k)
    monos = monomials(p)
    npoly = len(monos)
    conds = []
    for orbit in sorted(plan.keys()):
        rep, idx, moved = plan[orbit]
        P = X[idx]
        center = X[rep]
        dist = np.linalg.norm(P[:, None, :] - P[None, :, :], axis=-1)
        Mx = np.zeros((k + npoly, k + npoly))
        Mx[:k, :k] = dist ** m
        Q = eval_monos(P, monos)
        Mx[:k, k:] = Q
        Mx[k:, :k] = Q.T
        sv = np.linalg.svd(Mx, compute_uv=False)
        tol = max(Mx.shape) * np.finfo(float).eps * sv[0]
        eff_rank = int(np.sum(sv > tol))
        if eff_rank > 0:
            conds.append(float(sv[0] / sv[eff_rank - 1]))
        else:
            conds.append(float('inf'))
    arr = np.array(conds)
    return {
        "max": float(np.max(arr)),
        "median": float(np.median(arr)),
        "p95": float(np.percentile(arr, 95)),
        "p99": float(np.percentile(arr, 99)),
        "n_stencils": len(arr),
    }


def polynomial_reproduction_residual(X, oid, gid, elems,
                                      k=K_FIXED, m=RBF_M, p=RBF_P):
    """Max relative error of RBF-FD weights reproducing surface Laplacian
    on degree-1 monomials (x_0 through x_3).

    Exact: Delta_{S^3}(x_i) = -3 x_i  on the unit 3-sphere.
    """
    node_of, mult, plan = orbit_stencils(X, oid, gid, elems, k)
    max_rel = 0.0
    for orbit in sorted(plan.keys()):
        rep, idx, moved = plan[orbit]
        P = X[idx]
        center = X[rep]
        w = rbf_row(center, P, m, p)
        for dim in range(4):
            approx = float(np.dot(w, P[:, dim]))
            exact = -3.0 * center[dim]
            err = abs(approx - exact)
            ref = 3.0
            max_rel = max(max_rel, err / ref)
    return max_rel


def constant_function_residual(L_R0):
    """||L_R0 * 1||_inf — Laplacian of the constant function."""
    ones = np.ones(L_R0.shape[1])
    return float(np.max(np.abs(L_R0 @ ones)))


# ── Combined gate ──────────────────────────────────────────────────────

def cloud_gate(X, oid, gid, elems, reps_R0,
               k=K_FIXED, m=RBF_M, p=RBF_P, print_fn=print):
    """Full cloud admissibility gate on the R0 (scalar) operator.

    Returns (passed, diagnostics_dict).
    """
    n_seeds = int(np.max(oid)) + 1

    h = compute_fill_distance(X)
    q = compute_separation_radius(X)
    h_ref = fill_distance_reference(X)
    mesh_ratio = h / q if q > 1e-30 else float('inf')

    print_fn(f"    geometry: h={h:.4f}, q={q:.4f}, h/q={mesh_ratio:.1f}, "
             f"h_ref={h_ref:.4f}")

    sc = stencil_conditioning_stats(X, oid, gid, elems, k, m, p)
    print_fn(f"    stencil cond: max={sc['max']:.2e}, "
             f"median={sc['median']:.2e}, p99={sc['p99']:.2e}")

    reprod = polynomial_reproduction_residual(X, oid, gid, elems, k, m, p)
    print_fn(f"    poly reproduction residual: {reprod:.2e}")

    L_R0, _ = build_L_bundle(X, oid, gid, elems, reps_R0, k=k, m=m, p=p)
    L_R0 = np.real(L_R0)
    L_norm = float(np.linalg.norm(L_R0, 2))
    h2L = h ** 2 * L_norm
    const_res = constant_function_residual(L_R0)
    print_fn(f"    R0: ||L||={L_norm:.1f}, h^2||L||={h2L:.2f}, "
             f"const_res={const_res:.2e}")

    diag = {
        "n_seeds": n_seeds, "N": len(X),
        "h": h, "q": q, "h_ref": h_ref,
        "h_stability": abs(h - h_ref) / max(h_ref, 1e-15),
        "mesh_ratio": mesh_ratio,
        "stencil_cond_max": sc["max"],
        "stencil_cond_median": sc["median"],
        "stencil_cond_p95": sc["p95"],
        "stencil_cond_p99": sc["p99"],
        "L_R0_norm": L_norm,
        "h2_L": h2L,
        "const_residual": const_res,
        "reprod_residual": reprod,
    }

    reasons = []
    if mesh_ratio > GATE["mesh_ratio_max"]:
        reasons.append(f"mesh_ratio={mesh_ratio:.1f}")
    if sc["max"] > GATE["stencil_cond_max"]:
        reasons.append(f"stencil_cond_max={sc['max']:.2e}")
    if const_res > GATE["const_residual_max"]:
        reasons.append(f"const_residual={const_res:.2e}")
    if reprod > GATE["reprod_max"]:
        reasons.append(f"reprod_residual={reprod:.2e}")

    diag["pass"] = len(reasons) == 0
    diag["fail_reasons"] = reasons
    return diag["pass"], diag


# ── Extraction failure detection ───────────────────────────────────────

def is_extraction_ok(ex_result):
    """Check if spectral extraction is valid for contamination analysis.

    ||P_spec|| = 0 is a hard failure (Section 4 of LADDER_RULE_TASK).
    A nonzero projector has spectral radius 1, so ||P|| >= 1 under any
    induced norm.  ||P_spec|| = 0 means extraction failed catastrophically.
    """
    if not ex_result.get("pass"):
        return False
    riesz = ex_result.get("riesz", {})
    if riesz.get("rank_mismatch"):
        return False
    P_spec_norm = riesz.get("P_spec_norm")
    if P_spec_norm is None or P_spec_norm < 0.5:
        return False
    return True
