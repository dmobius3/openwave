"""Hopf-coordinate product rule on S^3 and mode evaluation.

Protocol section 4.1: x = (cos(eta)cos(xi1), cos(eta)sin(xi1), sin(eta)cos(xi2), sin(eta)sin(xi2)),
uniform 4N+1 points in each of xi1, xi2, Gauss-Legendre 2N+1 points in u = sin^2(eta),
weights normalized to unit total measure; node count (2N+1)(4N+1)^2.
"""
import numpy as np
from m85c_group import pi_unitary, DIM_RHO, FIBRES, unitarized_sym_power, quat_to_su2


def hopf_rule(D):
    """Product rule exact for polynomials of degree <= D on S^3 (unit total measure)."""
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


def eval_pi_unitary_all(x, n_max):
    """Evaluate pi_unitary(n, x) for n = 0..n_max. Returns list of matrices."""
    U = quat_to_su2(np.asarray(x, float))
    result = [np.array([[1.0 + 0j]])]
    for n in range(1, n_max + 1):
        result.append(unitarized_sym_power(U, n))
    return result


def eval_scalar_modes(X, n_max):
    """Evaluate all scalar harmonic modes at nodes X, levels 0..n_max.
    Returns: node x mode complex matrix. Mode ordering: level-major, then
    (n+1)^2 matrix elements in row-major order."""
    npts = len(X)
    mode_list = []
    for n in range(n_max + 1):
        mode_list.extend(range((n + 1)**2))
    nmodes = sum((n + 1)**2 for n in range(n_max + 1))
    Y = np.zeros((npts, nmodes), dtype=complex)

    offset = 0
    for i, x in enumerate(X):
        mats = eval_pi_unitary_all(x, n_max)
        off = 0
        for n in range(n_max + 1):
            d = (n + 1)**2
            Y[i, off:off + d] = mats[n].flatten()
            off += d

    return Y


def eval_section_level(X, rho_name, n, intertwiners):
    """Evaluate the W_rho-valued section basis functions at nodes X for level n.
    Section f_{A_i,j}(x) = A_i pi_n(x) e_j, where A_i are the intertwiners.
    Returns: npts x (m * d_n * d_rho) complex matrix, with column ordering:
    intertwiner-index i, multiplet-index j, component-index c."""
    d_rho = DIM_RHO[rho_name]
    d_n = n + 1
    m = len(intertwiners)
    npts = len(X)
    ncols = m * d_n * d_rho

    vals = np.zeros((npts, ncols), dtype=complex)

    for p, x in enumerate(X):
        pin = pi_unitary(n, x)
        col = 0
        for i, A in enumerate(intertwiners):
            for j in range(d_n):
                vec = A @ pin @ np.eye(d_n)[:, j]
                vals[p, col:col + d_rho] = vec
                col += d_rho

    return vals


def eval_r0_modes(X, r0_levels, n_max):
    """Evaluate R0-invariant modes at nodes X.
    r0_levels: dict {n: B_n} where B_n is (n+1) x m_n matrix of invariant vectors.
    Returns: npts x total_modes complex matrix."""
    npts = len(X)
    mode_list = []
    for n in sorted(r0_levels.keys()):
        if n > n_max:
            break
        B = r0_levels[n]
        m = B.shape[1]
        mode_list.append((n, B, m))

    total = sum(entry[2] * (entry[0] + 1) for entry in mode_list)
    Y = np.zeros((npts, total), dtype=complex)

    for p, x in enumerate(X):
        col = 0
        for n, B, m in mode_list:
            pin = pi_unitary(n, x)
            for i in range(m):
                a = B[:, i]
                for j in range(n + 1):
                    Y[p, col] = np.dot(a.conj(), pin[:, j])
                    col += 1

    return Y
