"""Section bases for each sector and cutoff, per protocol section 3.

Construction: stacked-constraint SVD in the frozen element order, Lowdin
symmetric orthonormalization, sign rule (first nonzero component positive).
"""
import numpy as np
import hashlib
import scipy.linalg
from m85c_group import (
    G120, KLASS, IRREP_NAMES, D_RHO, DIM_RHO, CHARTAB,
    FIBRES, GEN_PAIR, pi_unitary, compute_intertwiners, multiplicity,
    qmul, quat_to_su2, sym_power, unitarized_sym_power, binom_sqrt_diag
)


def intertwiner_basis_deterministic(rho_name, n):
    """Compute an orthonormal basis for Hom_2I(V_n, W_rho) using the deterministic
    stacked-constraint SVD construction. Returns list of d_rho x (n+1) matrices.
    Uses raw (monomial-basis) sym_power for the constraint to avoid conditioning
    issues at high n, then converts to the unitarized basis."""
    d_rho = DIM_RHO[rho_name]
    d_n = n + 1
    rho_mats = FIBRES[rho_name]
    T = binom_sqrt_diag(n)

    m_expected = multiplicity(rho_name, n)
    if m_expected == 0:
        return [], 0

    rows = []
    for idx in GEN_PAIR:
        rho_g = rho_mats[idx]
        U_g = quat_to_su2(np.asarray(G120[idx], float))
        pi_g_raw = sym_power(U_g, n)
        constraint = np.kron(rho_g, np.eye(d_n)) - np.kron(np.eye(d_rho), pi_g_raw.T)
        rows.append(constraint)
    stacked = np.vstack(rows)

    U, s, Vh = np.linalg.svd(stacked, full_matrices=True)
    ncols = stacked.shape[1]
    null_count = m_expected
    boundary = ncols - null_count
    if boundary > 0 and null_count > 0:
        gap = s[boundary - 1] / max(s[boundary] if boundary < len(s) else 1e-16, 1e-300)
        assert gap > 1e2, (f"{rho_name} n={n}: SVD gap ratio {gap:.2e} too small "
                           f"(s[{boundary-1}]={s[boundary-1]:.4e}, "
                           f"s[{boundary}]={s[boundary] if boundary < len(s) else 0:.4e})")
    null_vecs = Vh[ncols - null_count:]

    intertwiners = []
    for i in range(null_count):
        B = null_vecs[i].conj().reshape(d_rho, d_n)
        A = B * T[None, :]
        intertwiners.append(A)

    if null_count > 1:
        gram = np.array([[np.trace(A.conj().T @ B) for B in intertwiners] for A in intertwiners])
        eigvals, eigvecs = np.linalg.eigh(gram)
        G_sqrt_inv = eigvecs @ np.diag(1.0 / np.sqrt(eigvals)) @ eigvecs.conj().T
        new_ints = []
        for i in range(null_count):
            A = sum(G_sqrt_inv[j, i] * intertwiners[j] for j in range(null_count))
            new_ints.append(A)
        intertwiners = new_ints

    signed = []
    for A in intertwiners:
        flat = A.flatten()
        for c in flat:
            if abs(c) > 1e-14:
                phase = c / abs(c)
                A = A / phase
                break
        signed.append(A)

    return signed, null_count


def level_gram(rho_name, n, N_quad):
    """Compute the Gram matrix for the section space at level n, sector rho,
    using the production quadrature at degree N_quad."""
    from m85c_quadrature import hopf_rule, eval_section_level

    X, W = hopf_rule(N_quad)
    ints, m = intertwiner_basis_deterministic(rho_name, n)
    if m == 0:
        return np.zeros((0, 0))

    d_rho = DIM_RHO[rho_name]
    d_n = n + 1
    basis_dim = m * d_n * d_rho

    vals = eval_section_level(X, rho_name, n, ints)
    gram = vals.conj().T @ (W[:, None] * vals)
    return gram


def lowdin_orthonormalize(gram, vecs):
    """Lowdin symmetric orthonormalization: S^{-1/2} applied to the basis."""
    if gram.shape[0] == 0:
        return vecs, gram
    eigvals, eigvecs = np.linalg.eigh(gram)
    assert np.all(eigvals > 1e-14), f"Gram matrix nearly singular: min eigval {eigvals.min()}"
    S_inv_sqrt = eigvecs @ np.diag(1.0 / np.sqrt(eigvals)) @ eigvecs.conj().T
    new_vecs = S_inv_sqrt @ vecs
    return new_vecs, S_inv_sqrt


def build_sector_basis(rho_name, N_cutoff):
    """Build the complete orthonormalized basis for sector rho at cutoff N.
    Returns: dict with level-keyed intertwiner arrays, total mode count,
    and the basis hash."""
    levels = {}
    total_modes = 0

    for n in range(N_cutoff + 1):
        m = multiplicity(rho_name, n)
        if m == 0:
            continue
        ints, nc = intertwiner_basis_deterministic(rho_name, n)
        assert nc == m, f"{rho_name} n={n}: SVD null dim {nc} != character mult {m}"
        levels[n] = ints
        total_modes += m * (n + 1) * DIM_RHO[rho_name]

    return levels, total_modes


def hash_basis_object(levels, rho_name):
    """SHA-256 of the basis object: level-major, intertwiner-index, then flat matrix entries
    as little-endian float64 (real, imag) pairs."""
    h = hashlib.sha256()
    h.update(rho_name.encode())
    for n in sorted(levels.keys()):
        h.update(n.to_bytes(4, 'little'))
        for A in levels[n]:
            flat = A.flatten()
            for c in flat:
                h.update(np.float64(c.real).tobytes())
                h.update(np.float64(c.imag).tobytes())
    return h.hexdigest()


def build_r0_basis(N_cutoff):
    """Build the R0 (2I-invariant) basis at cutoff N.
    Uses the Reynolds projector to find invariant subspaces per level."""
    levels = {}
    total_modes = 0

    for n in range(N_cutoff + 1):
        m = multiplicity("R0", n)
        if m == 0:
            continue
        reynolds = np.zeros((n + 1, n + 1), dtype=complex)
        for q in G120:
            reynolds += pi_unitary(n, q)
        reynolds /= 120.0

        eigvals, eigvecs = np.linalg.eigh(reynolds)
        mask = eigvals > 0.5
        assert mask.sum() == m, f"R0 n={n}: Reynolds rank {mask.sum()} != mult {m}"
        B = eigvecs[:, mask]

        signed_cols = []
        for j in range(m):
            col = B[:, j].copy()
            for c in col:
                if abs(c) > 1e-14:
                    phase = c / abs(c)
                    col = col / phase
                    break
            signed_cols.append(col)
        B = np.column_stack(signed_cols)

        levels[n] = B
        total_modes += m * (n + 1)

    return levels, total_modes
