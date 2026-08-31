"""Fast separable Hopf-coordinate synthesis/analysis via FFT.

The unitarized representation pi_n_u factors in Hopf coordinates as:
  pi_u[a,b](eta, xi1, xi2) = d_u[a,b](eta) * exp(i*(n-a-b)*xi1) * exp(i*(b-a)*xi2)

where d_u[a,b](eta) = pi_n_unitary(n, [cos(eta), 0, sin(eta), 0])[a,b] is the
representation evaluated at a real SU(2) element.

The angular part is a 2D Fourier series handled by FFT, leaving only a 1D
transform in u = sin^2(eta) for each Fourier mode pair.

Complexity: O(nu * n^2 + K^2 log K) per level instead of O(n_nodes * n^2).
"""
import numpy as np
from build.group import pi_n_unitary, multiplicity, DIMS
from build.sections import total_modes


def _d_function_at_u(n, u_array):
    """Evaluate d_u[a,b](eta) = pi_n_unitary at real SU(2) elements for each u.

    Parameters
    ----------
    n : int
        Level (dimension = n+1).
    u_array : ndarray of shape (nu,)
        u = sin^2(eta) values.

    Returns
    -------
    d_all : ndarray of shape (nu, n+1, n+1)
    """
    nu = len(u_array)
    dim = n + 1
    d_all = np.empty((nu, dim, dim))
    for i in range(nu):
        ce = np.sqrt(1.0 - u_array[i])
        se = np.sqrt(u_array[i])
        q = np.array([ce, 0.0, se, 0.0])
        d_all[i] = pi_n_unitary(n, q).real
    return d_all


def fast_synthesis(coefficients, basis_obj, N, rho='R0', degree=None):
    """Fast synthesis using separable Hopf-coordinate FFT.

    Parameters
    ----------
    degree : int or None
        Override the quadrature degree (default: 4*N for production).

    Returns
    -------
    field_3d : ndarray of shape (nu, K, K, d_rho), complex
    u_array, wu : quadrature points and weights in u
    K : angular grid size
    """
    d_rho = DIMS[rho]
    D = degree if degree is not None else 4 * N
    K = D + 1
    nu_pts = D // 2 + 1

    xs, ws = np.polynomial.legendre.leggauss(nu_pts)
    u_array = (xs + 1) / 2
    wu = ws / 2

    levels = sorted(basis_obj.keys())
    fourier_3d = np.zeros((nu_pts, K, K, d_rho), dtype=complex)

    offset = 0
    for n in levels:
        A_all = basis_obj[n]
        mult = A_all.shape[0]
        dim_n = n + 1
        block_size = mult * dim_n

        c_block = coefficients[offset:offset + block_size].reshape(mult, dim_n)
        offset += block_size

        d_all = _d_function_at_u(n, u_array)

        # V[alpha, a, b] = sum_i A[i, alpha, a] * c[i, b]
        V = np.einsum('ira,ib->rab', A_all, c_block)

        for a in range(dim_n):
            for b in range(dim_n):
                p = (n - a - b) % K
                q = (b - a) % K
                for alpha in range(d_rho):
                    fourier_3d[:, p, q, alpha] += V[alpha, a, b] * d_all[:, a, b]

    field_3d = np.zeros((nu_pts, K, K, d_rho), dtype=complex)
    for alpha in range(d_rho):
        for ui in range(nu_pts):
            field_3d[ui, :, :, alpha] = np.fft.ifft2(fourier_3d[ui, :, :, alpha]) * K * K

    return field_3d, u_array, wu, K


def fast_analysis(field_3d, basis_obj, N, u_array, wu, K, rho='R0'):
    """Fast analysis: field values on Hopf grid -> Galerkin coefficients."""
    d_rho = DIMS[rho]
    nu_pts = len(u_array)
    levels = sorted(basis_obj.keys())

    fourier_3d = np.zeros_like(field_3d)
    for alpha in range(d_rho):
        for ui in range(nu_pts):
            fourier_3d[ui, :, :, alpha] = np.fft.fft2(field_3d[ui, :, :, alpha]) / (K * K)

    coeffs_list = []

    for n in levels:
        A_all = basis_obj[n]
        mult = A_all.shape[0]
        dim_n = n + 1

        d_all = _d_function_at_u(n, u_array)

        c_block = np.zeros((mult, dim_n), dtype=complex)

        for a in range(dim_n):
            for b in range(dim_n):
                p = (n - a - b) % K
                q = (b - a) % K

                for alpha in range(d_rho):
                    # w[u] * d[u,a,b] * F[u,p,q,alpha]  integrated over u
                    integrand = d_all[:, a, b] * fourier_3d[:, p, q, alpha]
                    integral = np.sum(wu * integrand)

                    # c[i,b] += conj(A[i,alpha,a]) * integral
                    for i in range(mult):
                        c_block[i, b] += np.conj(A_all[i, alpha, a]) * integral

        coeffs_list.append(c_block.ravel())

    return np.concatenate(coeffs_list)


def fast_project_cubic(coefficients, basis_obj, N, rho='R0'):
    """Compute P_N[|psi|^2 psi] via the fast Hopf-FFT transform route."""
    field_3d, u_array, wu, K = fast_synthesis(coefficients, basis_obj, N, rho=rho)

    abs2 = np.sum(np.abs(field_3d) ** 2, axis=-1, keepdims=True)
    cubic_3d = abs2 * field_3d

    return fast_analysis(cubic_3d, basis_obj, N, u_array, wu, K, rho=rho)
