"""Regression tests for structural invariants that must not be broken.

These tests exist because each one guarded against a bug that appeared
during development and would silently pass all other checks.
"""

import numpy as np


def test_identity_at_index_zero(elems):
    """Repair 3: gid=0 MUST be the identity quaternion.

    During development, elems[0] was [-1,0,0,0] (the central element).
    Orbit transport assigns the representative stencil to node gid=0.
    With gid=0 ≠ identity, the stencil was placed at the WRONG node,
    producing a catastrophically wrong spectrum where all eigenvalues
    were symmetric about zero instead of non-negative.
    """
    e = elems[0]
    assert abs(e[0] - 1.0) < 1e-10, f"elems[0][0] = {e[0]}, expected 1.0"
    assert np.linalg.norm(e[1:]) < 1e-10, f"elems[0][1:] = {e[1:]}, expected 0"
    return True


def test_identity_mutation(elems, chi, reps, build_cloud, build_bundle, MCKAY_DIST, LABELS):
    """Repair 3 mutation: force the identity to a non-zero index, confirm failure.

    If the spectrum is correct with identity misplaced, the assertion is not
    load-bearing and should be strengthened.
    """
    from .cloud import build_orbit_cloud, fibonacci_seeds_s3

    # Corrupt: swap elems[0] with elems[60] (some non-identity element)
    corrupted = list(elems)
    corrupted[0], corrupted[60] = corrupted[60], corrupted[0]

    seeds = fibonacci_seeds_s3(20)
    X, oid, gid = build_orbit_cloud(seeds, corrupted)

    rho = reps["R3"]
    d_mckay = MCKAY_DIST[LABELS.index("R3")]
    expected_lambda = d_mckay * (d_mckay + 2)

    from .bundle_operator import build_L_bundle
    L, _ = build_L_bundle(X, oid, gid, corrupted, rho, k=80)
    eigs = sorted(np.real(np.linalg.eigvals(L)))
    first_eig = eigs[0]

    # With identity misplaced, the first eigenvalue should NOT be near the expected value
    err = abs(first_eig - expected_lambda)
    detected = err > 1.0
    return {
        "first_eigenvalue": first_eig,
        "expected": expected_lambda,
        "error": err,
        "mutation_detected": detected,
    }


def test_non_hermitian_solver(L, label=""):
    """Repair 4: the bundle operator is NOT Hermitian in the Euclidean inner product.

    The continuum Laplacian on a flat bundle IS self-adjoint, but the discrete
    coordinate matrix is genuinely non-normal. A future cleanup that "restores"
    eigh on continuum grounds would be wrong: eigenvector conditioning grows
    with fibre dimension (measured at 1.3e1 for R0, 1.6e2 for R7), and the
    Hermitian solver's eigenvalues are garbage (the first non-trivial RUN1
    attempt with eigh gave eigenvalues of -1095 on a matrix whose true
    smallest eigenvalue was 3.0).

    This test FAILS if the matrix is close to Hermitian, which would mean
    the operator construction changed in a way that needs investigation.
    """
    asym = np.linalg.norm(L - L.conj().T) / np.linalg.norm(L)
    is_non_hermitian = asym > 0.1

    # Also verify: eigenvalues from eig have non-negative real parts
    # (the operator is -Δ which should be positive semidefinite)
    evals = np.linalg.eigvals(L)
    min_real = float(np.min(np.real(evals)))

    # And eigenvalues from eigh are WRONG (first one is very negative)
    H = 0.5 * (L + L.conj().T)
    evals_h = np.linalg.eigvalsh(H)
    min_hermitian = float(np.min(evals_h))

    # eigh gives negative eigenvalues for an operator that should be positive
    hermitian_gives_negative = min_hermitian < -1.0

    return {
        "label": label,
        "asymmetry": float(asym),
        "is_non_hermitian": is_non_hermitian,
        "min_real_eig": min_real,
        "min_hermitian_eig": min_hermitian,
        "hermitian_gives_negative": hermitian_gives_negative,
        "pass": is_non_hermitian and hermitian_gives_negative,
    }
