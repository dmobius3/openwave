"""Sector bases: intertwiner spaces Hom_2I(V_n, W_rho) and section bases for
all sectors at all required levels.

The deterministic construction per § 3: stacked-constraint SVD over the
generating pair, Lowdin orthonormalization against the per-level Gram, sign
fixed by the first nonzero component.

Each section f_{i,j}(x) = A_i * pi_n(x) * e_j takes values in W_rho,
where A_i is d_rho x (n+1), pi_n(x) is the unitarized Sym^n at point x,
and e_j is the j-th standard basis vector.
"""
import numpy as np
import hashlib
from scipy.linalg import sqrtm, svd as scipy_svd
from build.group import (G120, GEN_INDICES, IRREP_NAMES, DIMS, D_RHO,
                          multiplicity, pi_n_unitary, pi_n_mono, _sym_gram,
                          rep_rho, ALL_REPS, qmul, quat_to_su2, sym_power)


def compute_intertwiners(rho, n):
    """Compute Hom_2I(V_n, W_rho) via stacked equivariance constraints
    rho(g) A = A pi_n(g) over the generating pair, in the UNITARIZED basis.

    Uses the CG-recurrence-based pi_n_unitary for numerical stability at all n.

    Returns array of shape (mult, d_rho, n+1) in the UNITARIZED bases.
    """
    d_rho = DIMS[rho]
    dim_n = n + 1
    if d_rho == 0 or dim_n == 0:
        return np.zeros((0, 0, 0), dtype=complex)

    m = multiplicity(rho, n)
    if m == 0:
        return np.zeros((0, d_rho, dim_n), dtype=complex)

    vec_dim = d_rho * dim_n

    rows = []
    for gen_idx in GEN_INDICES:
        rho_g = ALL_REPS[rho][gen_idx]
        pi_g = pi_n_unitary(n, G120[gen_idx])
        C = np.kron(np.eye(dim_n), rho_g) - np.kron(pi_g.T, np.eye(d_rho))
        rows.append(C)
    M = np.vstack(rows)

    U, S, Vh = scipy_svd(M, full_matrices=True, lapack_driver='gesvd')
    tol = max(1e-10 * S[0], 1e-12) if len(S) > 0 else 1e-12
    null_mask = S < tol
    n_null_from_sv = null_mask.sum()
    n_trailing = max(0, Vh.shape[0] - len(S))
    n_null = n_null_from_sv + n_trailing

    assert n_null == m, f"expected mult {m} for {rho} at n={n}, got nullity {n_null}"

    null_vecs = Vh[-n_null:].conj()
    # kron uses column-major vectorization: reshape accordingly
    intertwiners = np.array([v.reshape(dim_n, d_rho).T for v in null_vecs])

    # sign fix: first nonzero component has positive real part
    for i in range(m):
        flat = intertwiners[i].ravel()
        for c in flat:
            if abs(c) > 1e-14:
                phase = c / abs(c)
                intertwiners[i] /= phase
                break

    return intertwiners


def level_gram_analytic(n):
    """Per-level Gram matrix for sections in the UNITARIZED basis.

    For the unitarized representation, Schur orthogonality gives:
    integral pi_n_u(x)^H pi_n_u(x) dx = I / (n+1)

    So for sections f_{i,j}(x) = A_i pi_n_u(x) e_j:
    <f_{i,j}, f_{k,l}> = sum_alpha,b,c [A_i^H]_{b,alpha} [A_k]_{alpha,c} delta_{jl} / (n+1)
    = delta_{jl} [A_i^H A_k]_{j...} ...

    Actually, more carefully:
    <f_{i,j}, f_{k,l}> = integral sum_a f_{i,j,a}(x)* f_{k,l,a}(x) dx
    = integral sum_a (sum_b [A_i]_{a,b} [pi_n(x)]_{b,j})* (sum_c [A_k]_{a,c} [pi_n(x)]_{c,l}) dx
    = sum_{a,b,c} [A_i]_{a,b}* [A_k]_{a,c} integral [pi_n(x)]_{b,j}* [pi_n(x)]_{c,l} dx

    For UNITARIZED pi_n: integral [pi_n_u]_{bj}* [pi_n_u]_{cl} dx = delta_{bc} delta_{jl} / (n+1)

    So <f_{i,j}, f_{k,l}> = delta_{jl} / (n+1) * sum_a [A_i]_{a,j}* [A_k]_{a,j}
    Wait, that's not right. Let me redo:
    = sum_{a,b,c} [A_i*]_{a,b} [A_k]_{a,c} * delta_{bc} delta_{jl} / (n+1)
    = delta_{jl} / (n+1) * sum_{a,b} [A_i*]_{a,b} [A_k]_{a,b}
    = delta_{jl} / (n+1) * tr(A_i^H A_k)

    Hmm, but that means the Gram in the (i,j) indices is:
    G_{(i,j),(k,l)} = delta_{jl} / (n+1) * tr(A_i^H A_k)

    This is block-diagonal in j,l with the intertwiner overlap matrix
    H_{ik} = tr(A_i^H A_k) / (n+1) replicated (n+1) times.

    Lowdin orthonormalization on this Gram makes the basis orthonormal.
    """
    pass  # we compute it numerically via quadrature for safety


def hopf_rule(D):
    """Hopf-coordinate product rule exact for polynomials of degree <= D on S^3.
    Returns (nodes, weights) with weights normalized to sum to 1."""
    K = D + 1
    nu = D // 2 + 1
    xs, ws = np.polynomial.legendre.leggauss(nu)
    u = (xs + 1) / 2
    wu = ws / 2
    xi = 2 * np.pi * np.arange(K) / K
    ce, se = np.sqrt(1 - u), np.sqrt(u)
    X, W = [], []
    for cu, su, w in zip(ce, se, wu):
        for a in xi:
            for b in xi:
                X.append([cu * np.cos(a), cu * np.sin(a),
                          su * np.cos(b), su * np.sin(b)])
                W.append(w)
    X = np.array(X)
    W = np.array(W)
    W /= W.sum()
    return X, W


def eval_pi_n_unitary_at_nodes(n, nodes):
    """Evaluate pi_n_unitary at each node. Returns array of shape (n_nodes, n+1, n+1)."""
    g = _sym_gram(n)
    Dh = np.diag(np.sqrt(g))
    Dih = np.diag(1.0 / np.sqrt(g))
    result = np.zeros((len(nodes), n + 1, n + 1), dtype=complex)
    for i, x in enumerate(nodes):
        M = quat_to_su2(np.asarray(x, float))
        result[i] = Dh @ sym_power(M, n) @ Dih
    return result


def build_basis_object(rho, N_max):
    """Build the § 3 basis object for sector rho at all levels n = 0..N_max.

    Returns a dict: level n -> array of shape (total_modes_at_n, d_rho)
    where total_modes_at_n = mult(rho, n) * (n+1), representing the
    section coefficients in the unitarized harmonic basis.

    The basis is Lowdin-orthonormalized and sign-fixed per § 3's construction.
    """
    d_rho = DIMS[rho]
    basis = {}

    for n in range(N_max + 1):
        m = multiplicity(rho, n)
        if m == 0:
            continue
        dim_n = n + 1
        total = m * dim_n

        # compute intertwiners (shape: m x d_rho x dim_n)
        A = compute_intertwiners(rho, n)
        assert A.shape == (m, d_rho, dim_n)

        # Gram matrix for the section basis:
        # G_{(i,j),(k,l)} = delta_{jl} * tr(A_i^H A_k) / (n+1)
        # This is block diagonal: for each j, the m x m block is H / (n+1)
        # where H_{ik} = tr(A_i^H A_k)
        H = np.zeros((m, m), dtype=complex)
        for i in range(m):
            for k in range(m):
                H[i, k] = np.trace(A[i].conj().T @ A[k])
        H /= (n + 1)

        # Lowdin: B = H^{-1/2}, applied to each j-block independently
        # Since H is the same for all j, this simplifies
        H_half = sqrtm(H)
        try:
            B = np.linalg.inv(H_half)
        except np.linalg.LinAlgError:
            raise RuntimeError(f"singular Gram at {rho}, n={n}")

        # orthonormalized intertwiners: A_orth = sum_k B_{ik} A_k
        A_orth = np.einsum('ik,kab->iab', B, A)

        # verify orthonormality
        H_check = np.zeros((m, m), dtype=complex)
        for i in range(m):
            for k in range(m):
                H_check[i, k] = np.trace(A_orth[i].conj().T @ A_orth[k])
        H_check /= (n + 1)
        assert np.abs(H_check - np.eye(m)).max() < 1e-10, \
            f"Lowdin failed at {rho}, n={n}: {np.abs(H_check - np.eye(m)).max():.2e}"

        # sign fix each intertwiner
        for i in range(m):
            flat = A_orth[i].ravel()
            for c in flat:
                if abs(c) > 1e-14:
                    phase = c / abs(c)
                    A_orth[i] /= phase
                    break

        basis[n] = A_orth  # shape (m, d_rho, dim_n)

    return basis


def basis_hash(basis_obj, rho):
    """SHA-256 of the basis object's canonical serialization.
    Coefficients as little-endian float64 (real, imag) pairs in level-major,
    then intertwiner-index, then row, then column order."""
    h = hashlib.sha256()
    for n in sorted(basis_obj.keys()):
        A = basis_obj[n]  # (m, d_rho, dim_n)
        for i in range(A.shape[0]):
            for r in range(A.shape[1]):
                for c in range(A.shape[2]):
                    z = A[i, r, c]
                    h.update(np.float64(z.real).tobytes())
                    h.update(np.float64(z.imag).tobytes())
    return h.hexdigest()


def total_modes(rho, N):
    """Total complex Galerkin modes in sector rho at cutoff N."""
    return sum(multiplicity(rho, n) * (n + 1) for n in range(N + 1))
