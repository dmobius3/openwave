"""Linear operators for the spectral Galerkin system on S^3 with binary
icosahedral (2I) symmetry.

Assembles:
  - Laplacian eigenvalue array for any sector rho at cutoff N
  - Diagonal linear operator (Laplacian in the Galerkin basis)
  - Real Jacobian of the cubic residual R(phi; omega) = Delta phi - c1 <phi,phi> phi + omega^2 phi
  - Diagnostics: self-adjointness defect K and imaginary-spectrum defect J
"""
import numpy as np
from build.group import IRREP_NAMES, DIMS, multiplicity
from build.sections import total_modes


def laplacian_eigenvalues(rho, N):
    """Laplacian eigenvalue array for sector rho at cutoff N.

    Returns a 1-D real array of length total_modes(rho, N).  Entry i carries
    the eigenvalue n(n+2) for the level n to which basis function i belongs.
    Modes are ordered level-major, then intertwiner-index, then multiplet-index.
    """
    eigs = []
    for n in range(N + 1):
        m = multiplicity(rho, n)
        if m == 0:
            continue
        count = m * (n + 1)          # modes at this level
        eigs.extend([n * (n + 2)] * count)
    eigs = np.array(eigs, dtype=float)
    assert len(eigs) == total_modes(rho, N), (
        f"eigenvalue count {len(eigs)} != total_modes({rho}, {N}) = {total_modes(rho, N)}"
    )
    return eigs


def linear_operator(rho, N):
    """Diagonal linear operator for sector rho at cutoff N.

    Returns a diagonal matrix of shape (m, m) where m = total_modes(rho, N),
    with diagonal entries -n(n+2) for each basis function at level n.
    This is the Galerkin representation of the Laplacian on S^3.
    """
    eigs = laplacian_eigenvalues(rho, N)
    return np.diag(-eigs)


def real_jacobian(phi, omega, rho, N, c1=1.0):
    """Real Jacobian of the cubic residual R(phi; omega) at state (phi, omega).

    The residual is:
        R(phi; omega) = Delta phi - c1 <phi,phi> phi + omega^2 phi

    where <phi,phi> = phi^H phi (Hermitian inner product, real-valued) and
    Delta is diagonal with entries -n(n+2) in the Galerkin basis.

    Because phi -> <phi,phi> phi is NOT complex-differentiable, the Jacobian
    is real-linear on R^{2m} (m = number of complex modes).  Writing
    phi = x + iy and v = [x; y] in R^{2m}:

        J = diag([l; l]) + omega^2 I_{2m} - c1 <phi,phi> I_{2m} - 2 c1 v v^T

    where l is the vector of Laplacian eigenvalues -n(n+2).

    Parameters
    ----------
    phi : complex array of length m
        Current complex state vector.
    omega : float
        Frequency parameter (the Jacobian depends on omega^2).
    rho : str
        Sector name (e.g. 'R0').
    N : int
        Spectral cutoff.
    c1 : float
        Nonlinear coupling constant (default 1.0).

    Returns
    -------
    J : ndarray of shape (2m, 2m)
        The real Jacobian matrix.
    """
    phi = np.asarray(phi, dtype=complex).ravel()
    m = total_modes(rho, N)
    assert len(phi) == m, f"phi has {len(phi)} entries, expected {m}"

    l = -laplacian_eigenvalues(rho, N)          # diagonal of Laplacian: -n(n+2)
    phi_norm_sq = np.vdot(phi, phi).real        # <phi, phi>
    v = np.concatenate([phi.real, phi.imag])     # real representation

    # Diagonal part: (L + omega^2 I - c1 <phi,phi> I) in real coords
    diag_val = l + omega**2 - c1 * phi_norm_sq
    diag_block = np.concatenate([diag_val, diag_val])
    J = np.diag(diag_block)

    # Rank-1 update: -2 c1 v v^T
    J -= 2.0 * c1 * np.outer(v, v)

    return J


def K_defect(A):
    """Self-adjointness defect: ||A - A^H||_2 (spectral norm of the skew part).

    For a real matrix, A^H = A^T, so this measures ||A - A^T||_2.
    """
    A = np.asarray(A)
    return np.linalg.norm(A - A.conj().T, 2)


def J_defect(A):
    """Imaginary-spectrum defect: max |Im lambda(A)| over all eigenvalues.

    A self-adjoint operator has purely real spectrum, so this measures
    departure from self-adjointness in spectral terms.
    """
    A = np.asarray(A)
    eigvals = np.linalg.eigvals(A)
    return float(np.max(np.abs(eigvals.imag)))
