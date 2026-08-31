"""Hopf-coordinate quadrature on S^3 for spectral Galerkin methods.

Product rules on S^3 using Hopf coordinates:
    x = (cos eta cos xi1, cos eta sin xi1, sin eta cos xi2, sin eta sin xi2)
with measure proportional to cos(eta) sin(eta) d(eta) d(xi1) d(xi2).
The substitution u = sin^2(eta) linearizes the radial part.

Quadrature nodes and weights, plus direct (non-FFT) synthesis and analysis
transforms for equivariant section bases.
"""

import numpy as np
from build.group import pi_n_unitary, quat_to_su2, _sym_gram, DIMS, multiplicity


# ---------------------------------------------------------------------------
# Quadrature rules
# ---------------------------------------------------------------------------

def hopf_rule(D):
    """Hopf-coordinate product rule exact for polynomials of degree <= D on S^3.

    Parameters
    ----------
    D : int
        Maximum polynomial degree for which the rule is exact.

    Returns
    -------
    nodes : ndarray, shape (N_nodes, 4)
        Points on S^3 in Cartesian coordinates.
    weights : ndarray, shape (N_nodes,)
        Quadrature weights, normalized so that weights.sum() == 1.

    Notes
    -----
    The rule uses:
    - K = D+1 uniform points in each of xi1, xi2 in [0, 2*pi)
      (exact for trigonometric polynomials of degree <= D via the DFT).
    - nu = D//2 + 1 Gauss-Legendre points in u = sin^2(eta) in [0, 1]
      (exact for polynomials of degree <= 2*nu - 1 >= D in u;
       the change of variables from eta to u absorbs the cos(eta)*sin(eta)
       Jacobian, so the u-integrand is a polynomial of degree <= D//2 in u
       when the angular part has been eliminated by the trapezoidal rule).

    Node count = nu * K^2.
    """
    K = D + 1
    nu = D // 2 + 1

    # Gauss-Legendre on [-1, 1], mapped to u in [0, 1]
    xs, ws = np.polynomial.legendre.leggauss(nu)
    u = (xs + 1) / 2
    wu = ws / 2

    # Uniform grid in each angle
    xi = 2 * np.pi * np.arange(K) / K

    # Precompute cos(eta), sin(eta) from u = sin^2(eta)
    ce = np.sqrt(1 - u)   # cos(eta)
    se = np.sqrt(u)        # sin(eta)

    # Build node array and weight array
    n_nodes = nu * K * K
    nodes = np.empty((n_nodes, 4))
    weights = np.empty(n_nodes)

    idx = 0
    for i_u in range(nu):
        c, s, w = ce[i_u], se[i_u], wu[i_u]
        for i_a in range(K):
            ca, sa = np.cos(xi[i_a]), np.sin(xi[i_a])
            for i_b in range(K):
                cb, sb = np.cos(xi[i_b]), np.sin(xi[i_b])
                nodes[idx, 0] = c * ca
                nodes[idx, 1] = c * sa
                nodes[idx, 2] = s * cb
                nodes[idx, 3] = s * sb
                weights[idx] = w
                idx += 1

    weights /= weights.sum()
    return nodes, weights


def production_rule(N):
    """4N-exact Hopf rule for Galerkin production integrands at band limit N.

    The cubic nonlinearity |psi|^2 psi with psi band-limited to level N
    produces integrands of polynomial degree <= 4N when projected onto
    modes of degree <= N.

    Node count: nu * K^2 where K = 4N+1, nu = 2N+1.
    So node count = (2N+1) * (4N+1)^2.

    Parameters
    ----------
    N : int
        Band limit of the field.

    Returns
    -------
    nodes : ndarray, shape (N_nodes, 4)
    weights : ndarray, shape (N_nodes,)
    """
    return hopf_rule(4 * N)


def monitor_rule(N):
    """6N-exact Hopf rule for cascade-monitor integrands at band limit N.

    The cascade monitor measures the band N < n <= 3N of |psi|^2 psi.
    Projecting onto modes of degree up to 3N against a degree-3N integrand
    requires exactness to degree 6N.

    Node count: nu * K^2 where K = 6N+1, nu = 3N+1.
    So node count = (3N+1) * (6N+1)^2.

    Parameters
    ----------
    N : int
        Band limit of the field.

    Returns
    -------
    nodes : ndarray, shape (N_nodes, 4)
    weights : ndarray, shape (N_nodes,)
    """
    return hopf_rule(6 * N)


# ---------------------------------------------------------------------------
# Representation evaluation at nodes
# ---------------------------------------------------------------------------

def eval_pi_n_unitary_at_nodes(n, nodes):
    """Evaluate the unitarized n-th symmetric power at each quadrature node.

    Parameters
    ----------
    n : int
        Representation level (dimension = n+1).
    nodes : ndarray, shape (N_nodes, 4)
        Points on S^3 (unit quaternions).

    Returns
    -------
    pi_vals : ndarray, shape (N_nodes, n+1, n+1), dtype complex
        pi_vals[i] = pi_n_unitary(n, nodes[i]).
    """
    if n == 0:
        return np.ones((len(nodes), 1, 1), dtype=complex)

    n_nodes = len(nodes)
    dim = n + 1
    result = np.empty((n_nodes, dim, dim), dtype=complex)

    for i in range(n_nodes):
        result[i] = pi_n_unitary(n, nodes[i])

    return result


# ---------------------------------------------------------------------------
# Section evaluation
# ---------------------------------------------------------------------------

def eval_section_at_nodes(intertwiner, n, nodes, pi_cache=None):
    """Evaluate one intertwiner's section basis functions at quadrature nodes.

    For intertwiner A of shape (d_rho, n+1), the section basis functions are:
        f_{A,j}(x) = A . pi_n_unitary(x) . e_j    for j = 0, ..., n

    giving a d_rho-valued function for each j.

    Parameters
    ----------
    intertwiner : ndarray, shape (d_rho, n+1)
        The intertwiner matrix A.
    n : int
        Representation level.
    nodes : ndarray, shape (N_nodes, 4)
        Quadrature nodes on S^3.
    pi_cache : ndarray, shape (N_nodes, n+1, n+1) or None
        Precomputed pi_n_unitary values. If None, they are computed here.

    Returns
    -------
    values : ndarray, shape (N_nodes, d_rho, n+1)
        values[i, a, j] = sum_b A[a,b] * pi_n_unitary(nodes[i])[b,j].
    """
    if pi_cache is None:
        pi_cache = eval_pi_n_unitary_at_nodes(n, nodes)

    # intertwiner is (d_rho, dim_n), pi_cache is (N_nodes, dim_n, dim_n)
    # result[i] = A @ pi_cache[i], shape (d_rho, dim_n)
    # Use einsum: result[i, a, j] = sum_b A[a,b] * pi_cache[i,b,j]
    values = np.einsum('ab,ibj->iaj', intertwiner, pi_cache)
    return values


# ---------------------------------------------------------------------------
# Direct synthesis (coefficients -> field values at nodes)
# ---------------------------------------------------------------------------

def synthesis(coefficients, basis_obj, nodes, rho=None):
    """Evaluate an equivariant field at quadrature nodes from Galerkin coefficients.

    The field is:
        psi(x) = sum_{n} sum_{i=0}^{m_n-1} sum_{j=0}^{n}
                     c_{n,i,j} * A_{n,i} . pi_n(x) . e_j

    where A_{n,i} are the orthonormalized intertwiners from basis_obj.

    Parameters
    ----------
    coefficients : ndarray, shape (N_total,) complex
        Galerkin coefficients in level-major, intertwiner-index, then j order.
        N_total = sum over active levels of mult(rho, n) * (n+1).
    basis_obj : dict
        level n -> ndarray of shape (mult, d_rho, n+1), the orthonormalized
        intertwiners. Produced by sections.build_basis_object().
    nodes : ndarray, shape (N_nodes, 4)
        Quadrature nodes on S^3.
    rho : str or None
        Sector name (for extracting d_rho). If None, inferred from basis_obj.

    Returns
    -------
    field : ndarray, shape (N_nodes, d_rho), dtype complex
        The field evaluated at each node.
    """
    # Determine d_rho from the first available level
    levels = sorted(basis_obj.keys())
    if len(levels) == 0:
        raise ValueError("empty basis_obj")
    d_rho = basis_obj[levels[0]].shape[1]

    n_nodes = len(nodes)
    field = np.zeros((n_nodes, d_rho), dtype=complex)

    # Walk through coefficients level by level
    offset = 0
    for n in levels:
        A_all = basis_obj[n]           # shape (mult, d_rho, n+1)
        m = A_all.shape[0]
        dim_n = n + 1
        block_size = m * dim_n

        # Extract this level's coefficients: shape (m, dim_n)
        c_block = coefficients[offset:offset + block_size].reshape(m, dim_n)
        offset += block_size

        # Compute pi_n at all nodes: shape (N_nodes, dim_n, dim_n)
        pi_vals = eval_pi_n_unitary_at_nodes(n, nodes)

        # For each intertwiner i and column j:
        #   contribution = c_{n,i,j} * A_{n,i} @ pi_n(x) @ e_j
        # Vectorized: for intertwiner i,
        #   A_i @ pi_vals[node] has shape (d_rho, dim_n)
        #   sum over j of c[i,j] * (A_i @ pi_vals)[node, :, j]
        #     = A_i @ pi_vals[node] @ c[i, :]

        # section_vals[i, node, :] = A_i @ pi_vals[node] @ c[i, :]
        for i in range(m):
            # A_all[i] @ pi_vals: shape (N_nodes, d_rho, dim_n)
            Apn = np.einsum('ab,nbj->naj', A_all[i], pi_vals)
            # Contract with c_block[i]: (N_nodes, d_rho)
            field += np.einsum('naj,j->na', Apn, c_block[i])

    return field


# ---------------------------------------------------------------------------
# Direct analysis (field values at nodes -> coefficients)
# ---------------------------------------------------------------------------

def analysis(field_values, basis_obj, nodes, weights):
    """Compute Galerkin coefficients from field values at quadrature nodes.

    This is the quadrature-weighted adjoint of synthesis:
        c_{n,i,j} = sum_nodes w_node * f_{n,i,j}(x_node)^H . psi(x_node)

    where f_{n,i,j}(x) = A_{n,i} . pi_n(x) . e_j is the section basis function.

    For orthonormal sections and an exact-enough rule, analysis(synthesis(c)) = c.

    Parameters
    ----------
    field_values : ndarray, shape (N_nodes, d_rho), dtype complex
        The field evaluated at each node.
    basis_obj : dict
        level n -> ndarray of shape (mult, d_rho, n+1).
    nodes : ndarray, shape (N_nodes, 4)
        Quadrature nodes on S^3.
    weights : ndarray, shape (N_nodes,)
        Quadrature weights (summing to 1).

    Returns
    -------
    coefficients : ndarray, shape (N_total,), dtype complex
        Galerkin coefficients in level-major order.
    """
    levels = sorted(basis_obj.keys())
    if len(levels) == 0:
        raise ValueError("empty basis_obj")

    # Weighted field values: shape (N_nodes, d_rho)
    weighted_field = weights[:, None] * field_values

    coeffs_list = []

    for n in levels:
        A_all = basis_obj[n]           # shape (mult, d_rho, n+1)
        m = A_all.shape[0]
        dim_n = n + 1

        # pi_n at all nodes: shape (N_nodes, dim_n, dim_n)
        pi_vals = eval_pi_n_unitary_at_nodes(n, nodes)

        c_block = np.zeros((m, dim_n), dtype=complex)

        for i in range(m):
            # Section basis values: f_{i,j}(x)[a] = sum_b A[a,b] pi[b,j]
            # f_vals shape: (N_nodes, d_rho, dim_n)
            f_vals = np.einsum('ab,nbj->naj', A_all[i], pi_vals)

            # Inner product: c[i,j] = sum_node sum_a w[node] * f[node,a,j]^* * psi[node,a]
            # = sum_node sum_a f_vals[node,a,j]^* * weighted_field[node,a]
            c_block[i] = np.einsum('naj,na->j', f_vals.conj(), weighted_field)

        # No normalization factor needed: the Lowdin-orthonormalized
        # intertwiners satisfy tr(A_i^H A_k) / (n+1) = delta_{ik},
        # which combined with Schur orthogonality (integral of
        # pi[b,j]^* pi[c,l] = delta_{bc} delta_{jl} / (n+1)) gives
        # orthonormal sections: <f_{i,j}, f_{k,l}> = delta_{ik} delta_{jl}.

        coeffs_list.append(c_block.ravel())

    return np.concatenate(coeffs_list)
