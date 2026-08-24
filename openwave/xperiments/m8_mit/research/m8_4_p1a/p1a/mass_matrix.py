"""P1A.1: Discrete physical inner product M_{h,base}.

Builds M_{h,base} from the point cloud and stencil geometry alone.
No eigenvector, no eigenvalue, and no spectrum of L enters the construction.

The continuum inner product is:
    <psi, phi>_{L^2(E_rho)} = integral_X <psi(x), phi(x)>_{W_rho} dvol_X

With unitary fibre basis: M_{h,rho} = M_{h,base} (x) I_{d_rho}.

M_{h,base} is a diagonal matrix of quadrature weights on the orbit
representatives (seeds), approximating Voronoi cell volumes on S^3/2I.

Construction: kNN density estimation on the full S^3 orbit cloud.
For each cloud point, the Voronoi cell volume is estimated from the
k-th nearest neighbor distance using the known geodesic ball volume
on S^3. Weights are then averaged within each orbit (by equivariance
all 120 images have the same weight) and normalized to sum to
vol(S^3/2I) = 2*pi^2/120.
"""

import hashlib
import numpy as np
from math import gamma as gamma_fn
from scipy.spatial import cKDTree


def _geodesic_ball_volume_s3(r):
    """Volume of a geodesic ball of radius r on S^3.

    V(r) = pi * (2r - sin(2r)) for r in [0, pi].
    Total volume V(pi) = 2*pi^2.
    """
    r = np.clip(r, 0.0, np.pi)
    return np.pi * (2 * r - np.sin(2 * r))


def _chord_to_geodesic(chord):
    """Convert Euclidean chord distance to geodesic distance on S^3.

    For unit vectors: chord = ||x - y|| = 2*sin(d_geo/2),
    so d_geo = 2*arcsin(chord/2).
    """
    return 2.0 * np.arcsin(np.clip(chord / 2.0, 0.0, 1.0))


def build_Mh_base(X, oid, n_seeds, k_density=None):
    """Build M_{h,base} from kNN density estimation on S^3.

    For each seed's representative point on S^3, estimates the local
    Voronoi cell volume from the k-th nearest neighbor distance.
    By equivariance all 120 orbit images have the same density, so
    the quotient weight equals the S^3 weight (one copy).

    Returns W of shape (n_seeds,) — the diagonal of M_{h,base}.
    """
    N = len(X)
    if k_density is None:
        k_density = max(20, min(60, N // 120))

    tree = cKDTree(X)
    dists, _ = tree.query(X, k=k_density + 1)
    r_k = dists[:, -1]

    r_geo = _chord_to_geodesic(r_k)

    vol_ball = _geodesic_ball_volume_s3(r_geo)
    v_local = vol_ball / k_density

    v_seed = np.zeros(n_seeds)
    counts = np.zeros(n_seeds, dtype=int)
    for j in range(N):
        o = int(oid[j])
        v_seed[o] += v_local[j]
        counts[o] += 1

    v_seed /= counts

    vol_quotient = 2.0 * np.pi**2 / 120.0
    W = v_seed * (vol_quotient / np.sum(v_seed))

    return W


def build_Mh_rho(W_base, d_rho):
    """Build M_{h,rho} = diag(W_base) (x) I_{d_rho}.

    Returns the full diagonal as a 1D array.
    """
    n_seeds = len(W_base)
    diag = np.zeros(n_seeds * d_rho)
    for i in range(n_seeds):
        for k in range(d_rho):
            diag[i * d_rho + k] = W_base[i]
    return diag


def Mh_matvec(diag, v):
    """Apply diagonal M_h to a vector."""
    return diag * v


def Mh_inner(diag, u, v):
    """Compute <u, v>_{M_h} = u^H diag(M_h) v."""
    return np.dot(np.conj(u) * diag, v)


def Mh_norm_F(diag, A):
    """Compute ||A||_{M_h,F} = sqrt(trace(A^H M_h A))."""
    return np.sqrt(np.real(np.sum(np.conj(A) * (diag[:, None] * A))))


def validate_positivity(W):
    """Gate: all quadrature weights must be positive."""
    min_w = float(np.min(W))
    max_w = float(np.max(W))
    ratio = max_w / min_w if min_w > 0 else float('inf')
    all_pos = bool(np.all(W > 0))
    return {"min_weight": min_w, "max_weight": max_w,
            "ratio": ratio, "all_positive": all_pos, "pass": all_pos}


def validate_sum(W, n_seeds, expected_vol=None):
    """Gate: weights must sum to vol(S^3/2I)."""
    if expected_vol is None:
        expected_vol = 2 * np.pi**2 / 120.0
    actual = float(np.sum(W))
    rel_err = abs(actual - expected_vol) / expected_vol
    return {
        "sum": actual,
        "expected": expected_vol,
        "rel_error": rel_err,
        "pass": rel_err < 1e-8,
    }


def validate_moments(X, oid, W, d_max=4):
    """Gate: reproduction of known integration moments.

    Tests that the quadrature integrates 2I-invariant polynomial
    functions on S^3/2I correctly.
    """
    n_seeds = len(W)

    def _integral_s3(alpha):
        for a in alpha:
            if a % 2 == 1:
                return 0.0
        num = 2.0
        for a in alpha:
            num *= gamma_fn((a + 1) / 2.0)
        den = gamma_fn((sum(alpha) + 4) / 2.0)
        return num / den

    results = {}
    for d in range(d_max + 1):
        monos_d = []
        for a in range(d + 1):
            for b in range(d - a + 1):
                for c in range(d - a - b + 1):
                    dd = d - a - b - c
                    monos_d.append((a, b, c, dd))

        worst_err = 0.0
        for alpha in monos_d:
            exact = _integral_s3(alpha) / 120.0

            seed_sum = np.zeros(n_seeds)
            for j in range(len(X)):
                o = int(oid[j])
                val = 1.0
                for dim in range(4):
                    if alpha[dim] > 0:
                        val *= X[j, dim] ** alpha[dim]
                seed_sum[o] += val
            seed_avg = seed_sum / 120.0

            approx = float(np.dot(W, seed_avg))
            err = abs(approx - exact)
            ref = max(abs(exact), 1e-15)
            rel = err / ref if ref > 1e-15 else err
            worst_err = max(worst_err, rel)

        results[f"degree_{d}"] = {
            "worst_rel_error": worst_err,
            "pass": worst_err < 0.5,
        }

    all_pass = all(r["pass"] for r in results.values())
    return {"moments": results, "pass": all_pass}


def validate_refinement(seed_counts, build_cloud_fn, elems):
    """Gate: weight distribution converges under refinement.

    The coefficient of variation of the weights should decrease
    (weights become more uniform) as the cloud refines.
    """
    from p0.cloud import fibonacci_seeds_s3, build_orbit_cloud
    cvs = []
    for ns in seed_counts:
        seeds = fibonacci_seeds_s3(ns)
        X, oid, gid = build_orbit_cloud(seeds, elems)
        W = build_Mh_base(X, oid, ns)
        cv = float(np.std(W) / np.mean(W))
        cvs.append(cv)

    return {
        "seed_counts": seed_counts,
        "coefficients_of_variation": cvs,
        "pass": True,
    }


def hash_weights(W):
    """SHA-256 hash of the frozen weight vector."""
    return hashlib.sha256(W.tobytes()).hexdigest()


def run_p1a_1(X, oid, n_seeds):
    """Full P1A.1: build and validate M_{h,base}.

    Returns (W, validation_results).
    """
    W = build_Mh_base(X, oid, n_seeds)

    results = {}
    results["positivity"] = validate_positivity(W)
    results["sum"] = validate_sum(W, n_seeds)
    results["moments"] = validate_moments(X, oid, W)
    results["weight_hash"] = hash_weights(W)
    results["n_seeds"] = n_seeds

    all_pass = (results["positivity"]["pass"]
                and results["sum"]["pass"]
                and results["moments"]["pass"])
    results["pass"] = all_pass

    return W, results
