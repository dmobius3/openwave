"""P0.0: Construct and certify the nine irreducible representations of 2I.

Six irreps (R0, R1, R3, R6, R7, R8) are direct restrictions of V_n (n=0..5).
Three (R4, R5, R2) are extracted from V_6 and V_7 by character projection.

Every representation is unitary with respect to a stated Hermitian metric.
"""

import numpy as np
from .algebra import quat_to_su2, sym_power, invariant_gram, su2_character
from .group import (LABELS, DIMS, MCKAY_DIST, DIRECT_LEVEL, PROJECTED,
                    build_character_table, class_of)


def _isotypic_projector(elems, chi_table, irrep_idx, ambient_reps):
    """Character-based isotypic projector for irrep sigma in ambient V_n.

    P_σ = (d_σ / |G|) Σ_g conj(χ_σ(g)) ρ_ambient(g)

    Returns a d_ambient × d_ambient matrix.
    """
    d_sigma = DIMS[irrep_idx]
    d_amb = ambient_reps[0].shape[0]
    P = np.zeros((d_amb, d_amb), dtype=complex)
    for gi, g in enumerate(elems):
        c = class_of(g)
        chi_val = chi_table[irrep_idx, c]
        P += np.conj(chi_val) * ambient_reps[gi]
    P *= d_sigma / len(elems)
    return P


def _extract_subrep(elems, ambient_reps, projector, expected_dim, gram=None):
    """Extract an irreducible sub-representation from the projector image.

    If gram is provided, orthonormalize with respect to it.
    Returns (basis, restricted_reps) where basis is d_amb × d_sigma.
    """
    if gram is not None:
        S = np.linalg.cholesky(gram)
        P_orth = S @ projector @ np.linalg.inv(S)
    else:
        S = np.eye(projector.shape[0])
        P_orth = projector

    w, V = np.linalg.eigh(0.5 * (P_orth + P_orth.conj().T))
    idx = np.where(np.abs(w - 1.0) < 0.1)[0]
    assert len(idx) == expected_dim, \
        f"projector rank {len(idx)} != expected {expected_dim}"

    B_orth = V[:, idx]
    if gram is not None:
        B = np.linalg.solve(S, B_orth)
    else:
        B = B_orth

    B, _ = np.linalg.qr(B)

    reps = []
    for gi in range(len(elems)):
        R = B.conj().T @ ambient_reps[gi] @ B
        reps.append(R)
    return B, reps


def _orthonormalize_ambient(raw_reps, n):
    """Transform V_n representations to be unitary w.r.t. standard inner product.

    S = sqrt(invariant_gram(n)), then ρ̃(g) = S ρ(g) S⁻¹.
    """
    H = invariant_gram(n)
    S = np.diag(np.sqrt(np.diag(H)))
    Si = np.diag(1.0 / np.sqrt(np.diag(H)))
    return [S @ R @ Si for R in raw_reps]


def build_all_representations(elems, chi_table):
    """Construct explicit matrix representations for all nine 2I irreps.

    All representations are unitary with respect to the standard inner
    product (identity Gram matrix).

    Returns a dict: label -> list of d_ρ × d_ρ complex matrices,
    one per group element, in the same order as `elems`.
    """
    reps = {}
    bases = {}

    for label, n in DIRECT_LEVEL.items():
        raw = []
        for g in elems:
            U = quat_to_su2(g)
            raw.append(sym_power(U, n))
        reps[label] = _orthonormalize_ambient(raw, n)
        bases[label] = np.eye(n + 1, dtype=complex)

    for label, n in PROJECTED.items():
        raw = []
        for g in elems:
            U = quat_to_su2(g)
            raw.append(sym_power(U, n))
        ambient = _orthonormalize_ambient(raw, n)

        irrep_idx = LABELS.index(label)
        P = _isotypic_projector(elems, chi_table, irrep_idx, ambient)
        B, sub_reps = _extract_subrep(
            elems, ambient, P, DIMS[irrep_idx])
        reps[label] = sub_reps
        bases[label] = B

    return reps, bases


def certify_representations(elems, chi_table, reps, mult_table):
    """P0.0 certification: all six required checks.

    Returns a dict with results for each check.
    """
    results = {}
    G = len(elems)

    # 1. Correct complex dimension
    dim_check = {}
    for i, label in enumerate(LABELS):
        d_actual = reps[label][0].shape[0]
        d_expected = DIMS[i]
        dim_check[label] = {"actual": d_actual, "expected": d_expected,
                            "pass": d_actual == d_expected}
    results["dimensions"] = dim_check

    # 2. ρ(e) = I
    identity_idx = None
    for gi, g in enumerate(elems):
        if abs(g[0] - 1.0) < 1e-10 and np.linalg.norm(g[1:]) < 1e-10:
            identity_idx = gi
            break
    assert identity_idx is not None
    id_check = {}
    for label in LABELS:
        d = DIMS[LABELS.index(label)]
        err = float(np.max(np.abs(reps[label][identity_idx] - np.eye(d))))
        id_check[label] = {"residual": err, "pass": err < 1e-12}
    results["identity"] = id_check

    # 3. ρ(gh) = ρ(g)ρ(h) across the FULL multiplication table
    homomorphism_check = {}
    for label in LABELS:
        worst = 0.0
        count = 0
        for i in range(G):
            for j in range(G):
                k = mult_table[i][j]
                prod = reps[label][i] @ reps[label][j]
                err = float(np.max(np.abs(prod - reps[label][k])))
                worst = max(worst, err)
                count += 1
        homomorphism_check[label] = {
            "pairs_checked": count,
            "worst_residual": worst,
            "pass": worst < 1e-8}
    results["homomorphism"] = homomorphism_check

    # 4. Unitarity: ρ(g)† ρ(g) = I (all reps orthonormalized)
    unitarity_check = {}
    for label in LABELS:
        d = DIMS[LABELS.index(label)]
        worst = 0.0
        for gi in range(G):
            R = reps[label][gi]
            err = float(np.max(np.abs(R.conj().T @ R - np.eye(d))))
            worst = max(worst, err)
        unitarity_check[label] = {
            "gram_matrix": "I (orthonormalized basis)",
            "worst_residual": worst,
            "pass": worst < 1e-8}
    results["unitarity"] = unitarity_check

    # 5. Traces match the character table
    trace_check = {}
    for i, label in enumerate(LABELS):
        worst = 0.0
        for gi, g in enumerate(elems):
            tr = float(np.real(np.trace(reps[label][gi])))
            c = class_of(g)
            expected = chi_table[i, c]
            err = abs(tr - expected)
            worst = max(worst, err)
        trace_check[label] = {
            "worst_residual": worst,
            "pass": worst < 1e-6}
    results["traces"] = trace_check

    # 6. Label correspondence stated explicitly
    label_map = {}
    for i, label in enumerate(LABELS):
        label_map[label] = {
            "dimension": DIMS[i],
            "mckay_distance": MCKAY_DIST[i],
            "source": f"V_{DIRECT_LEVEL[label]}|_2I" if label in DIRECT_LEVEL
                      else f"projection from V_{PROJECTED[label]}|_2I",
        }
    results["label_correspondence"] = label_map

    return results


def mutation_test(elems, chi_table, reps, mult_table):
    """P0.0 mutation: swap two representations and corrupt one.

    Mutation A: swap R3 and R4 (both dim 3) — traces must go red.
    Mutation B: replace R1 with a random unitary — homomorphism must go red.
    """
    results = {}

    # Mutation A: swap R3 and R4
    swapped = dict(reps)
    swapped["R3"], swapped["R4"] = reps["R4"], reps["R3"]
    for label in ["R3", "R4"]:
        i = LABELS.index(label)
        worst = 0.0
        for gi, g in enumerate(elems):
            tr = float(np.real(np.trace(swapped[label][gi])))
            c = class_of(g)
            expected = chi_table[i, c]
            err = abs(tr - expected)
            worst = max(worst, err)
        results[f"swap_{label}_trace_err"] = worst
    results["swap_detected"] = (results["swap_R3_trace_err"] > 0.1 or
                                results["swap_R4_trace_err"] > 0.1)

    # Mutation B: corrupt R1 with a random unitary
    rng = np.random.default_rng(42)
    H = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
    Q, _ = np.linalg.qr(H)
    corrupted_R1 = [Q @ R @ Q.conj().T for R in reps["R1"]]
    corrupted_R1[0] = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
    worst_hom = 0.0
    for i in range(min(20, len(elems))):
        for j in range(min(20, len(elems))):
            k = mult_table[i][j]
            prod = corrupted_R1[i] @ corrupted_R1[j]
            err = float(np.max(np.abs(prod - corrupted_R1[k])))
            worst_hom = max(worst_hom, err)
    results["corrupt_homomorphism_err"] = worst_hom
    results["corrupt_detected"] = worst_hom > 0.1

    results["pass"] = results["swap_detected"] and results["corrupt_detected"]
    return results


def covariance_check(elems, reps, mult_table):
    """Verify that ρ(g) → U ρ(g) U† leaves spectral conclusions unchanged.

    Pick a random unitary U for each representation, compute the
    conjugated representation, and verify that the invariant dimensions
    from the full-group averaging projector are identical.
    """
    rng = np.random.default_rng(99)
    results = {}
    for label in LABELS:
        d = DIMS[LABELS.index(label)]
        H = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
        U, _ = np.linalg.qr(H)

        orig_mats = reps[label]
        conj_mats = [U @ R @ U.conj().T for R in orig_mats]

        G = len(elems)
        orig_proj = sum(orig_mats) / G
        conj_proj = sum(conj_mats) / G

        orig_eigs = sorted(np.abs(np.linalg.eigvals(orig_proj)))
        conj_eigs = sorted(np.abs(np.linalg.eigvals(conj_proj)))

        err = float(max(abs(a - b) for a, b in zip(orig_eigs, conj_eigs)))
        results[label] = {"basis_change_eigenvalue_err": err,
                          "pass": err < 1e-10}

    results["all_pass"] = all(r["pass"] for l, r in results.items()
                              if l in LABELS)
    return results
