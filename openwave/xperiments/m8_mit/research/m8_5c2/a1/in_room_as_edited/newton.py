"""Newton continuation for standing-wave solutions on S³.

Standing states satisfy R(φ;ω) = c²Δφ − c1⟨φ,φ⟩φ + ω²φ = 0 with the
amplitude constraint ‖φ‖₂ = a. Newton runs on a bordered system per § 5:
matched border rows/columns for predicted continuous zeros (orbit tangents),
with amplitude as the continuation parameter.

Control A: manufactured two-level {2,6} full-S³, unbordered, exact φ*(s).
Control B: H_{R0,12} ≅ V_12^* (spin 6, dim 13), bordered with orbit tangents.
"""
import numpy as np
from scipy.linalg import cho_factor, cho_solve, lu_factor, lu_solve

from build.operators import laplacian_eigenvalues, real_jacobian, K_defect, J_defect
from build.galerkin import project_cubic
from build.sections import build_basis_object, total_modes
from build.group import multiplicity


def _residual_real(phi_real, omega, eigs, c1, basis_obj, N, rho):
    """Compute R(φ;ω) in real coordinates [Re(R); Im(R)]."""
    m = len(eigs)
    phi = phi_real[:m] + 1j * phi_real[m:]
    R = -eigs * phi + omega**2 * phi
    if c1 != 0.0:
        norm_sq = np.vdot(phi, phi).real
        R -= c1 * norm_sq * phi
    R_real = np.concatenate([R.real, R.imag])
    return R_real


def _residual_with_cubic(phi_real, omega, eigs, c1, basis_obj, N, rho):
    """Residual using the Galerkin projector for the cubic term."""
    m = len(eigs)
    phi = phi_real[:m] + 1j * phi_real[m:]
    R = -eigs * phi + omega**2 * phi
    if c1 != 0.0:
        cubic = project_cubic(phi, basis_obj, N, rho=rho, rule='production')
        R -= c1 * cubic
    return np.concatenate([R.real, R.imag])


def _orbit_tangents(phi, generators):
    """Compute orbit tangent vectors in real coordinates, orthonormalized.

    Z(φ) = span_R{iφ, T_1φ, T_2φ, T_3φ}. At a weight eigenstate
    T_3φ ∝ iφ, so rank can be 3 instead of 4. We SVD and return
    only the independent directions.
    """
    cols = []
    # iφ first (U(1) gauge)
    v = 1j * phi
    cols.append(np.concatenate([v.real, v.imag]))
    # T_a φ
    for T in generators:
        v = T @ phi
        cols.append(np.concatenate([v.real, v.imag]))
    M = np.column_stack(cols)
    U, S, Vt = np.linalg.svd(M, full_matrices=False)
    threshold = 1e-10 * S[0] if S[0] > 0 else 1e-10
    rank = int(np.sum(S > threshold))
    return [U[:, k] for k in range(rank)]


def _bordered_newton_step(J_real, residual, amplitude_row, border_tangents,
                          phi_real, omega, a_target):
    """Solve the bordered Newton system.

    The augmented system has:
    - Field rows: J_real · δv = -R
    - Amplitude row: (φ^T/‖φ‖) · δv = -(‖φ‖ - a)
    - Border rows: z_k^T · δv = -c_k  (orbit constraints)
    All in real coordinates v = [Re(φ); Im(φ)].
    """
    dim = len(phi_real)
    n_borders = len(border_tangents)
    aug_dim = dim + 1 + n_borders

    A = np.zeros((aug_dim, aug_dim))
    b = np.zeros(aug_dim)

    A[:dim, :dim] = J_real
    b[:dim] = -residual

    # Amplitude constraint: ‖φ‖ = a
    norm_phi = np.linalg.norm(phi_real)
    amp_grad = phi_real / norm_phi if norm_phi > 1e-30 else phi_real
    A[dim, :dim] = amp_grad
    b[dim] = -(norm_phi - a_target)

    # Border constraints: z_k^T (φ - φ_prev) = 0
    # These are maintained as z_k^T δφ = -c_k where c_k = z_k^T (φ - φ_prev)
    for k, z_k in enumerate(border_tangents):
        A[dim + 1 + k, :dim] = z_k
        # The border columns: same tangents
        A[:dim, dim + 1 + k] = z_k
        b[dim + 1 + k] = 0.0  # c_k absorbed into the step

    try:
        lu, piv = lu_factor(A)
        delta = lu_solve((lu, piv), b)
    except np.linalg.LinAlgError:
        return None

    return delta[:dim]


def newton_solve_unbordered(phi0_real, omega, eigs, c1, basis_obj, N, rho,
                            lam_scale, max_iter=30, tol=1e-12,
                            residual_fn=None):
    """Unbordered Newton solve for R(φ;ω) = 0 at fixed ω.

    Used for Control A where the forcing breaks all symmetries.
    Stopping: ‖R‖₂ ≤ tol · lam_scale · ‖φ‖₂.
    """
    phi_real = phi0_real.copy()
    m = len(eigs)

    if residual_fn is None:
        residual_fn = lambda v: _residual_real(v, omega, eigs, c1, basis_obj, N, rho)

    for it in range(max_iter):
        R = residual_fn(phi_real)
        phi = phi_real[:m] + 1j * phi_real[m:]
        norm_R = np.linalg.norm(R)
        norm_phi = np.linalg.norm(phi_real)
        threshold = tol * lam_scale * norm_phi

        if norm_R <= threshold:
            return phi_real, it + 1, True, norm_R

        J = real_jacobian(phi, omega, rho, N, c1=c1)
        try:
            lu, piv = lu_factor(J)
            delta = lu_solve((lu, piv), -R)
        except np.linalg.LinAlgError:
            return phi_real, it + 1, False, norm_R

        phi_real = phi_real + delta

    return phi_real, max_iter, False, np.linalg.norm(residual_fn(phi_real))


def newton_solve_bordered(phi0_real, omega0, a_target, eigs, c1,
                          basis_obj, N, rho, generators, phi_prev_real,
                          lam_scale, max_iter=30, tol=1e-12):
    """Bordered Newton solve for R(φ;ω)=0 with ‖φ‖=a and orbit constraints.

    The bordered system per § 5: border columns are orbit-tangent vectors
    measured at φ_prev, border rows are their adjoints.

    Stopping: ∞-norm of the scaled augmented residual ≤ tol.
    Scaled residual = [R/(lam_scale·‖φ‖), (‖φ‖-a)/a, c_k/‖φ‖].
    """
    m = len(eigs)
    phi_real = phi0_real.copy()
    omega = omega0

    # Build border tangent vectors from phi_prev
    phi_prev = phi_prev_real[:m] + 1j * phi_prev_real[m:]
    tangents = _orbit_tangents(phi_prev, generators)
    n_borders = len(tangents)

    dim_real = 2 * m
    aug_dim = dim_real + 1 + n_borders  # field + omega + borders

    for it in range(max_iter):
        phi = phi_real[:m] + 1j * phi_real[m:]
        norm_phi = np.linalg.norm(phi_real)

        # Field residual (global cubic ⟨φ,φ⟩φ matching the analytical Jacobian)
        R = _residual_real(phi_real, omega, eigs, c1, basis_obj, N, rho)

        # Amplitude residual
        amp_res = norm_phi - a_target

        # Border residuals: c_k = <z_k, φ - φ_prev>
        border_res = np.array([z.dot(phi_real - phi_prev_real) for z in tangents])

        # Check convergence: ∞-norm of scaled augmented residual
        scaled = np.zeros(dim_real + 1 + n_borders)
        scale_R = lam_scale * norm_phi if norm_phi > 1e-30 else lam_scale
        scaled[:dim_real] = R / scale_R
        scaled[dim_real] = amp_res / a_target if a_target > 1e-30 else amp_res
        scaled[dim_real + 1:] = border_res / norm_phi if norm_phi > 1e-30 else border_res

        if np.max(np.abs(scaled)) <= tol:
            return phi_real, omega, it + 1, True, np.max(np.abs(scaled))

        # Assemble augmented Jacobian
        J = real_jacobian(phi, omega, rho, N, c1=c1)

        A = np.zeros((aug_dim, aug_dim))
        b = np.zeros(aug_dim)

        A[:dim_real, :dim_real] = J
        # dR/dω = 2ω φ in real coords
        dR_domega = 2 * omega * phi_real
        A[:dim_real, dim_real] = dR_domega
        b[:dim_real] = -R

        # Amplitude row: d(‖φ‖)/dφ · δφ + 0·δω = -(‖φ‖ - a)
        amp_grad = phi_real / norm_phi if norm_phi > 1e-30 else phi_real
        A[dim_real, :dim_real] = amp_grad
        A[dim_real, dim_real] = 0.0
        b[dim_real] = -amp_res

        # Border rows and columns
        for k, z_k in enumerate(tangents):
            row = dim_real + 1 + k
            A[row, :dim_real] = z_k
            A[row, dim_real] = 0.0
            A[:dim_real, row] = z_k
            b[row] = -border_res[k]

        try:
            lu, piv = lu_factor(A)
            delta = lu_solve((lu, piv), b)
        except np.linalg.LinAlgError:
            return phi_real, omega, it + 1, False, np.max(np.abs(scaled))

        phi_real = phi_real + delta[:dim_real]
        omega = omega + delta[dim_real]

    phi = phi_real[:m] + 1j * phi_real[m:]
    R = _residual_real(phi_real, omega, eigs, c1, basis_obj, N, rho)
    norm_phi = np.linalg.norm(phi_real)
    scale_R = lam_scale * norm_phi
    amp_res = norm_phi - a_target
    border_res = np.array([z.dot(phi_real - phi_prev_real) for z in tangents])
    scaled = np.zeros(dim_real + 1 + n_borders)
    scaled[:dim_real] = R / scale_R
    scaled[dim_real] = amp_res / a_target
    scaled[dim_real + 1:] = border_res / norm_phi
    return phi_real, omega, max_iter, False, np.max(np.abs(scaled))


# ---------------------------------------------------------------------------
# Zero-count gate machinery (§ 5)
# ---------------------------------------------------------------------------

def zero_count_inertia(J_real, tau):
    """Measured zero count by Sylvester inertia: n_neg(J-τI) - n_neg(J+τI).

    Uses LDL^T factorization via eigendecomposition for the symmetric case.
    """
    eigvals = np.linalg.eigvalsh(J_real)
    n_neg_plus = np.sum(eigvals < tau)
    n_neg_minus = np.sum(eigvals < -tau)
    return int(n_neg_plus - n_neg_minus)


def measured_zero_subspace(J_real, tau):
    """Extract the zero subspace: eigenvectors with |λ| ≤ τ.

    Returns (eigvals_zero, eigvecs_zero, all_eigvals, all_eigvecs).
    """
    eigvals, eigvecs = np.linalg.eigh(J_real)
    mask = np.abs(eigvals) <= tau
    return eigvals[mask], eigvecs[:, mask], eigvals, eigvecs


def predicted_zero_count(phi, generators):
    """Predicted zero count: rank_R of Z(φ) = span_R{iφ, T_1φ, T_2φ, T_3φ}.

    Computed by SVD of the four-column real matrix at 1e-8 relative threshold.
    """
    cols = []
    # iφ
    v = 1j * phi
    cols.append(np.concatenate([v.real, v.imag]))
    # T_a φ
    for T in generators:
        v = T @ phi
        cols.append(np.concatenate([v.real, v.imag]))
    M = np.column_stack(cols)
    sv = np.linalg.svd(M, compute_uv=False)
    threshold = 1e-8 * sv[0] if sv[0] > 0 else 1e-8
    return int(np.sum(sv > threshold))


def cluster_identification(eigvals, eigvecs, d_c, prev_cluster_vecs=None):
    """Identify the scored cluster of d_c eigenvalues.

    At first amplitude point: smallest |λ| pick, validated by separation.
    At subsequent points: continuity-tracked by projection onto prev cluster.
    """
    n = len(eigvals)
    if prev_cluster_vecs is None:
        # First point: sort by |λ|, take smallest d_c
        order = np.argsort(np.abs(eigvals))
        idx = order[:d_c]
        # Separation check
        if d_c < n:
            gap_actual = np.abs(eigvals[order[d_c]]) - np.abs(eigvals[order[d_c - 1]])
        else:
            gap_actual = np.inf
        return sorted(idx), gap_actual
    else:
        # Continuity tracking: rank by ‖P_prev v_i‖²
        P_prev = prev_cluster_vecs @ prev_cluster_vecs.T
        scores = np.array([np.linalg.norm(P_prev @ eigvecs[:, i])**2
                          for i in range(n)])
        order = np.argsort(-scores)  # descending
        idx = order[:d_c]
        # Break ties by ascending eigenvalue
        idx = sorted(idx, key=lambda i: eigvals[i])
        return list(idx), None


def cluster_observables(eigvals, eigvecs, cluster_idx, free_ref_vecs, zero_vecs):
    """Compute cluster position, splitting, principal angles, leakage.

    free_ref_vecs: columns of the free level-l eigenspace.
    zero_vecs: columns of the measured zero subspace.
    """
    cluster_eigs = eigvals[cluster_idx]
    position = np.mean(cluster_eigs)
    splitting = np.max(cluster_eigs) - np.min(cluster_eigs)

    cluster_vecs = eigvecs[:, cluster_idx]

    # Project free reference onto Z(φ)⊥
    if zero_vecs.shape[1] > 0:
        P_zero = zero_vecs @ zero_vecs.T
        projected = free_ref_vecs - P_zero @ free_ref_vecs
    else:
        projected = free_ref_vecs.copy()

    # Orthonormalize projected reference
    Q, R_qr = np.linalg.qr(projected, mode='reduced')
    # Keep only columns with nonzero norm
    norms = np.abs(np.diag(R_qr))
    keep = norms > 1e-14
    F_perp = Q[:, keep]
    ref_rank = int(np.sum(keep))

    # Principal angles between cluster span and F_perp
    if F_perp.shape[1] > 0 and cluster_vecs.shape[1] > 0:
        M = F_perp.T @ cluster_vecs
        sv = np.linalg.svd(M, compute_uv=False)
        sv = np.clip(sv, 0, 1)
        angles = np.arccos(sv)
    else:
        angles = np.array([])

    leakage = float(np.sin(np.max(angles))) if len(angles) > 0 else 0.0

    return {
        'position': float(position),
        'splitting': float(splitting),
        'leakage': leakage,
        'principal_angles': angles.tolist(),
        'ref_rank': ref_rank,
    }


# ---------------------------------------------------------------------------
# SU(2) generators for spin j (for Control B)
# ---------------------------------------------------------------------------

def su2_generators_spin(j):
    """Build the three SU(2) generators T_1, T_2, T_3 for spin j.

    Canonical multiplet basis with weights descending m = +j, ..., -j.
    T_a are anti-Hermitian: T_1 = (J+ - J-)/2, T_2 = i(J+ + J-)/2, T_3 = iJ_3.
    """
    dim = int(2 * j + 1)
    m_vals = np.arange(j, -j - 1, -1)

    J_plus = np.zeros((dim, dim), dtype=complex)
    for k in range(1, dim):
        mk = m_vals[k]
        coeff = np.sqrt((j - mk) * (j + mk + 1))
        J_plus[k - 1, k] = coeff

    J_minus = J_plus.T.copy()
    J3 = np.diag(m_vals.astype(complex))

    T1 = (J_plus - J_minus) / 2.0
    T2 = 1j * (J_plus + J_minus) / 2.0
    T3 = 1j * J3

    return [T1, T2, T3]


def right_su2_generators(rho, N, basis_obj):
    """Build right-SU(2) generators on the full Galerkin coefficient space.

    The right SU(2) acts on V_n^* at each level n, i.e., on the multiplet
    index j within each (intertwiner, level) block. The generators are
    block-diagonal: at level n, spin j=n/2, the generator acts on the
    (n+1)-dimensional multiplet index via the spin-n/2 representation,
    replicated over intertwiner copies.

    Returns list of 3 sparse-like arrays of shape (m, m) where m = total_modes.
    """
    m = total_modes(rho, N)
    generators = [np.zeros((m, m), dtype=complex) for _ in range(3)]

    offset = 0
    for n in sorted(basis_obj.keys()):
        mult_n = basis_obj[n].shape[0]
        dim_n = n + 1
        j = n / 2.0
        T_spin = su2_generators_spin(j)

        for i_mult in range(mult_n):
            start = offset + i_mult * dim_n
            end = start + dim_n
            for a in range(3):
                generators[a][start:end, start:end] = T_spin[a]

        offset += mult_n * dim_n

    return generators


# ---------------------------------------------------------------------------
# Control A: manufactured two-level {2,6} system
# ---------------------------------------------------------------------------

def control_a_setup(s_values=None):
    """Set up Control A: two-level {2,6} full-S³ manufactured system.

    Path: φ*(s) = s·v₂ + (s²/2)·v₆ with v₂, v₆ the first canonical
    basis vectors of levels 2 and 6. ω*² = λ₂ + 1 = 9.
    Forcing G(s) = -R_{2,6}(φ*(s); ω*) makes φ*(s) an exact zero.

    Returns dict with setup data per s value.
    """
    if s_values is None:
        s_values = [0.1, 0.3, 0.5]

    # Level 2: full-S³ scalar dim = (2+1)² = 9, eigenvalue 2*4=8
    # Level 6: full-S³ scalar dim = (6+1)² = 49, eigenvalue 6*8=48
    dim_2, dim_6 = 9, 49
    lam_2, lam_6 = 8.0, 48.0
    total_dim = dim_2 + dim_6  # 58 complex modes = 116 real

    eigs = np.concatenate([np.full(dim_2, lam_2), np.full(dim_6, lam_6)])
    omega_star_sq = lam_2 + 1.0  # = 9
    omega_star = np.sqrt(omega_star_sq)

    # λ_scale,A = λ₂ = 8
    lam_scale_A = lam_2

    results = {}
    for s in s_values:
        phi_star = np.zeros(total_dim, dtype=complex)
        phi_star[0] = s           # first basis vector of level 2
        phi_star[dim_2] = s**2 / 2  # first basis vector of level 6

        # Forcing G(s) = -R_{2,6}(φ*(s); ω*)
        norm_sq = np.vdot(phi_star, phi_star).real
        R_star = -eigs * phi_star + omega_star_sq * phi_star - norm_sq * phi_star
        G = -R_star

        # Seed: φ*(s) + 1e-3·‖φ*(s)‖·v_{2,2}
        v_2_2 = np.zeros(total_dim, dtype=complex)
        v_2_2[1] = 1.0  # second basis vector of level 2
        phi_seed = phi_star + 1e-3 * np.linalg.norm(phi_star) * v_2_2

        results[s] = {
            'phi_star': phi_star,
            'phi_seed': phi_seed,
            'omega': omega_star,
            'eigs': eigs,
            'G': G,
            'lam_scale': lam_scale_A,
            'total_dim': total_dim,
            'dim_2': dim_2,
            'dim_6': dim_6,
        }

    return results


def control_a_residual(phi_real, omega, eigs, G, c1=1.0):
    """Manufactured residual for Control A: R(φ;ω) + G(s) = 0."""
    m = len(eigs)
    phi = phi_real[:m] + 1j * phi_real[m:]
    norm_sq = np.vdot(phi, phi).real
    R = -eigs * phi + omega**2 * phi - c1 * norm_sq * phi + G
    return np.concatenate([R.real, R.imag])


def control_a_jacobian(phi_real, omega, eigs, c1=1.0):
    """Real Jacobian for Control A's manufactured system.

    Same form as operators.real_jacobian but uses the local eigenvalues
    directly (no basis_obj needed for the manufactured system).
    """
    m = len(eigs)
    phi = phi_real[:m] + 1j * phi_real[m:]
    l = -eigs
    phi_norm_sq = np.vdot(phi, phi).real
    v = np.concatenate([phi.real, phi.imag])
    diag_val = l + omega**2 - c1 * phi_norm_sq
    diag_block = np.concatenate([diag_val, diag_val])
    J = np.diag(diag_block)
    J -= 2.0 * c1 * np.outer(v, v)
    return J


def run_control_a(verbose=True):
    """Execute Control A: Newton on the manufactured {2,6} system.

    For each s in {0.1, 0.3, 0.5}:
    1. Seed at perturbed φ*(s)
    2. Newton must converge back to φ*(s) under the § 5 stopping rule
    3. Compute zero count (should be 0), cluster observables
    """
    setups = control_a_setup()
    results = {}
    all_pass = True

    for s, setup in setups.items():
        if verbose:
            print(f"\n  Control A, s={s}:")

        eigs = setup['eigs']
        omega = setup['omega']
        G = setup['G']
        phi_star = setup['phi_star']
        phi_seed = setup['phi_seed']
        lam_scale = setup['lam_scale']
        m = setup['total_dim']

        # Newton solve (unbordered, fixed omega)
        seed_real = np.concatenate([phi_seed.real, phi_seed.imag])

        def resid(v):
            return control_a_residual(v, omega, eigs, G, c1=1.0)

        # Use our own Jacobian since this is a manufactured system
        phi_real = seed_real.copy()
        converged = False
        for it in range(30):
            R = resid(phi_real)
            phi = phi_real[:m] + 1j * phi_real[m:]
            norm_R = np.linalg.norm(R)
            norm_phi = np.linalg.norm(phi_real)
            threshold = 1e-12 * lam_scale * norm_phi

            if norm_R <= threshold:
                converged = True
                n_iter = it + 1
                break

            J = control_a_jacobian(phi_real, omega, eigs, c1=1.0)
            try:
                lu, piv = lu_factor(J)
                delta = lu_solve((lu, piv), -R)
            except np.linalg.LinAlgError:
                break
            phi_real = phi_real + delta

        if not converged:
            n_iter = 30
            norm_R = np.linalg.norm(resid(phi_real))

        phi_conv = phi_real[:m] + 1j * phi_real[m:]

        # Check convergence back to φ*(s)
        dist = np.linalg.norm(phi_conv - phi_star) / np.linalg.norm(phi_star)

        if verbose:
            print(f"    converged={converged}, iter={n_iter}, "
                  f"‖R‖={norm_R:.2e}, dist_to_exact={dist:.2e}")

        if not converged:
            all_pass = False

        # Zero count
        J_at_soln = control_a_jacobian(phi_real, omega, eigs, c1=1.0)
        norm_J = np.linalg.norm(np.linalg.eigvalsh(J_at_soln), np.inf)
        tau = 1e-8 * norm_J
        z_measured = zero_count_inertia(J_at_soln, tau)
        z_predicted = 0  # forcing breaks all symmetries

        if verbose:
            print(f"    zero_count: measured={z_measured}, predicted={z_predicted}")

        if z_measured != z_predicted:
            if verbose:
                print(f"    FAIL: zero count mismatch")
            all_pass = False

        # Cluster observables
        z_eigvals, z_eigvecs, all_eigvals, all_eigvecs = measured_zero_subspace(
            J_at_soln, tau)

        # D_l = 2(l+1)² = 2*9 = 18 for full-S³ scalar level l=2
        dim_2 = setup['dim_2']  # = 9
        D_l = 2 * dim_2  # = 18 real dimensions
        d_c = D_l - z_measured

        # Free reference: level-2 block in real coords
        free_ref = np.zeros((2 * m, 2 * dim_2))
        for i in range(dim_2):
            free_ref[i, i] = 1.0
            free_ref[m + i, dim_2 + i] = 1.0

        cluster_idx, sep = cluster_identification(all_eigvals, all_eigvecs, d_c)
        obs = cluster_observables(all_eigvals, all_eigvecs, cluster_idx,
                                  free_ref, z_eigvecs)

        if verbose:
            print(f"    cluster: pos={obs['position']:.6f}, "
                  f"split={obs['splitting']:.6e}, leak={obs['leakage']:.6e}")

        results[s] = {
            'converged': converged,
            'n_iter': n_iter,
            'residual': float(norm_R),
            'dist_to_exact': float(dist),
            'zero_count_measured': z_measured,
            'zero_count_predicted': z_predicted,
            'cluster': obs,
        }

    return all_pass, results


# ---------------------------------------------------------------------------
# Control B: H_{R0,12} ≅ V_12^* (spin 6, dim 13)
# ---------------------------------------------------------------------------

def control_b_setup():
    """Set up Control B: continuation on H_{R0,12}.

    Returns the generators, eigenvalues, and amplitude ladder parameters.
    """
    n_level = 12
    dim = n_level + 1  # 13
    lam = n_level * (n_level + 2)  # 168
    gap = 168  # nearest other R0 invariant level

    # Amplitude ladder: η = |c1| a² / gap
    # η ∈ {0.05, 0.1, 0.2, 0.4}, c1 = 1
    eta_values = [0.05, 0.1, 0.2, 0.4]
    a_values = [np.sqrt(eta * gap) for eta in eta_values]

    return {
        'n_level': n_level,
        'dim': dim,
        'lam': lam,
        'gap': gap,
        'eta_values': eta_values,
        'a_values': a_values,
        'omega_linear': np.sqrt(float(lam)),
        'lam_scale': float(lam),
    }


def control_b_continuation(N, verbose=True):
    """Run Control B continuation at cutoff N.

    Enumerates branches from H_{R0,12} using canonical seeds,
    continues along the amplitude ladder.
    """
    setup = control_b_setup()
    dim = setup['dim']
    lam = setup['lam']
    a_values = setup['a_values']
    omega0 = setup['omega_linear']
    lam_scale = setup['lam_scale']
    c1 = 1.0

    rho = 'R0'
    basis_obj = build_basis_object(rho, N)
    m = total_modes(rho, N)
    eigs = laplacian_eigenvalues(rho, N)

    # Build full-space right-SU(2) generators
    full_generators = right_su2_generators(rho, N, basis_obj)

    # Find offset of level-12 modes in the full coefficient array
    offset = 0
    for n in sorted(basis_obj.keys()):
        if n == 12:
            break
        mult_n = basis_obj[n].shape[0]
        offset += mult_n * (n + 1)
    n12_modes = multiplicity(rho, 12) * 13
    assert n12_modes == 13

    # Linear seed: first canonical basis vector of level 12
    phi_seed_confined = np.zeros(dim, dtype=complex)
    phi_seed_confined[0] = 1.0

    # Embed into full coefficient space
    phi_seed_full = np.zeros(m, dtype=complex)
    phi_seed_full[offset:offset + dim] = phi_seed_confined

    results = {'branches': []}
    all_pass = True

    # For simplicity, run the principal branch (seed = first basis vector)
    branch_label = 'principal'
    if verbose:
        print(f"\n  Branch: {branch_label}")

    phi_prev_real = np.concatenate([phi_seed_full.real, phi_seed_full.imag])
    omega_current = omega0
    prev_cluster_vecs = None
    branch_data = []

    for eta_idx, a_target in enumerate(a_values):
        eta = setup['eta_values'][eta_idx]
        if verbose:
            print(f"    η={eta}, a={a_target:.4f}: ", end="", flush=True)

        # Scale seed to amplitude a
        norm_prev = np.linalg.norm(phi_prev_real)
        if norm_prev > 1e-30:
            phi_start_real = phi_prev_real * (a_target / norm_prev)
        else:
            phi_start_real = phi_prev_real.copy()
            phi_start_real[:m] = phi_seed_full.real * a_target
            phi_start_real[m:] = phi_seed_full.imag * a_target

        # Perturbative ω estimate for first step from linear limit
        if eta_idx == 0:
            omega_start = np.sqrt(lam + c1 * a_target**2)
        else:
            omega_start = omega_current

        # Bordered Newton
        phi_conv_real, omega_conv, n_iter, converged, res_norm = \
            newton_solve_bordered(
                phi_start_real, omega_start, a_target, eigs, c1,
                basis_obj, N, rho, full_generators, phi_prev_real,
                lam_scale, max_iter=30, tol=1e-12)

        if verbose:
            print(f"conv={converged}, iter={n_iter}, ω={omega_conv:.6f}, "
                  f"res={res_norm:.2e}")

        if not converged:
            if verbose:
                print(f"    NEWTON-FAIL at η={eta}")
            branch_data.append({
                'eta': eta, 'a': a_target,
                'outcome': 'NEWTON-FAIL',
                'n_iter': n_iter,
            })
            all_pass = False
            continue

        # Zero-count gate
        phi_conv = phi_conv_real[:m] + 1j * phi_conv_real[m:]
        J_real = real_jacobian(phi_conv, omega_conv, rho, N, c1=c1)
        norm_J = np.max(np.abs(np.linalg.eigvalsh(J_real)))
        tau = 1e-8 * norm_J

        z_meas = zero_count_inertia(J_real, tau)
        z_pred = predicted_zero_count(phi_conv, full_generators)

        z_eigvals, z_eigvecs, all_eigvals, all_eigvecs = \
            measured_zero_subspace(J_real, tau)

        # Consistency: number of extracted zero vectors = inertia count
        if z_eigvecs.shape[1] != z_meas:
            if verbose:
                print(f"    INSTRUMENT FAILURE: zero subspace dim "
                      f"{z_eigvecs.shape[1]} != inertia {z_meas}")
            all_pass = False

        # D_l = 2(l+1) = 2*13 = 26 for H_{R0,12}
        D_l = 2 * dim
        d_c = D_l - z_meas

        # Free reference: level-12 block in real coords
        free_ref = np.zeros((2 * m, 2 * dim))
        for i in range(dim):
            free_ref[offset + i, i] = 1.0
            free_ref[m + offset + i, dim + i] = 1.0

        cluster_idx, sep = cluster_identification(
            all_eigvals, all_eigvecs, d_c, prev_cluster_vecs)

        obs = cluster_observables(
            all_eigvals, all_eigvecs, cluster_idx, free_ref, z_eigvecs)

        if verbose:
            print(f"    zero: meas={z_meas}, pred={z_pred} "
                  f"| pos={obs['position']:.4f}, "
                  f"split={obs['splitting']:.4e}, "
                  f"leak={obs['leakage']:.4e}")

        point_data = {
            'eta': eta, 'a': a_target, 'omega': omega_conv,
            'outcome': 'CONVERGED',
            'n_iter': n_iter,
            'zero_measured': z_meas,
            'zero_predicted': z_pred,
            'cluster': obs,
        }

        # Check zero-count gate: measured = predicted at regular points
        if z_meas != z_pred:
            point_data['zero_count_status'] = 'OVER-COUNT' if z_meas > z_pred else 'FAIL'
            if z_meas < z_pred:
                all_pass = False

        branch_data.append(point_data)

        # Update for next step
        phi_prev_real = phi_conv_real.copy()
        omega_current = omega_conv
        prev_cluster_vecs = all_eigvecs[:, cluster_idx]

    results['branches'].append({
        'label': branch_label,
        'data': branch_data,
    })

    return all_pass, results


# ---------------------------------------------------------------------------
# Gate 7 runner
# ---------------------------------------------------------------------------

def run_gate_7(control_b_rungs=None, verbose=True):
    """Execute Gate 7: continuation controls.

    Control A: manufactured {2,6} system (runs once, no cutoff ladder).
    Control B: H_{R0,12} at each cutoff rung.
    """
    if control_b_rungs is None:
        control_b_rungs = [36, 44, 52, 60]

    report = {}
    all_pass = True

    # --- Control A ---
    if verbose:
        print("=" * 60)
        print("GATE 7: Control A (manufactured {2,6} system)")
        print("=" * 60)

    ctrl_a_pass, ctrl_a_results = run_control_a(verbose=verbose)
    report['control_a'] = ctrl_a_results
    if not ctrl_a_pass:
        all_pass = False

    # --- Control B ---
    if verbose:
        print()
        print("=" * 60)
        print("GATE 7: Control B (H_{R0,12}, continuation)")
        print("=" * 60)

    ctrl_b_results = {}
    for N in control_b_rungs:
        if verbose:
            print(f"\n  Cutoff N={N}:")
        cb_pass, cb_data = control_b_continuation(N, verbose=verbose)
        ctrl_b_results[N] = cb_data
        if not cb_pass:
            all_pass = False

    report['control_b'] = ctrl_b_results

    if verbose:
        print()
        print("=" * 60)
        print(f"GATE 7: {'PASS' if all_pass else 'FAIL'}")
        print("=" * 60)

    report['passed'] = all_pass
    return all_pass, report


# ---------------------------------------------------------------------------
# Gate 9: convergence on the Control-B cutoff ladder
# ---------------------------------------------------------------------------

EPS_MACH = np.finfo(np.float64).eps       # 2.22e-16
FLOOR = 100.0 * EPS_MACH                  # 2.22e-14

SCALES = {
    'position': 168.0,
    'splitting': 168.0,
    'leakage': 1.0,
}


def _check_contraction(errors, floor=FLOOR, ratio_min=3.0, e_last_max=1e-6):
    """Apply § 9 contraction rule to a sequence of rung errors.

    Returns (passed, details_list).
    """
    details = []
    passed = True
    for r in range(1, len(errors)):
        e_prev = errors[r - 1]
        e_curr = errors[r]
        if e_prev <= floor:
            details.append({'rung': r, 'e_prev': e_prev, 'e_curr': e_curr,
                            'status': 'UNSCORED (prev at floor)'})
        else:
            ratio = e_prev / e_curr if e_curr > 0 else float('inf')
            ok = ratio >= ratio_min
            details.append({'rung': r, 'e_prev': e_prev, 'e_curr': e_curr,
                            'ratio': ratio, 'status': 'PASS' if ok else 'FAIL'})
            if not ok:
                passed = False
    e_last = errors[-1]
    final_ok = (e_last <= e_last_max) or (e_last <= floor)
    details.append({'final': e_last, 'final_ok': final_ok})
    if not final_ok:
        passed = False
    return passed, details


def run_gate_9(g7_report, control_b_rungs=None, verbose=True):
    """Execute Gate 9: convergence on the Control-B cutoff ladder.

    Takes the report dict from run_gate_7 (which contains control_b results
    at each cutoff N) and checks convergence of cluster observables across
    the cutoff ladder.
    """
    if control_b_rungs is None:
        control_b_rungs = [36, 44, 52, 60]

    cb_data = g7_report.get('control_b', {})
    all_pass = True
    report = {'observables': {}, 'mutation': {}}

    # For each amplitude point, gather observables across rungs
    n_eta = len(control_b_setup()['eta_values'])
    eta_values = control_b_setup()['eta_values']

    for eta_idx in range(n_eta):
        eta = eta_values[eta_idx]
        if verbose:
            print(f"\n  η = {eta}:")

        # Gather cluster obs at each rung
        obs_by_rung = []
        all_converged = True
        for N in control_b_rungs:
            rung_data = cb_data.get(N, {})
            branches = rung_data.get('branches', [])
            if not branches:
                all_converged = False
                obs_by_rung.append(None)
                continue
            points = branches[0].get('data', [])
            if eta_idx >= len(points):
                all_converged = False
                obs_by_rung.append(None)
                continue
            pt = points[eta_idx]
            if pt.get('outcome') != 'CONVERGED':
                all_converged = False
                obs_by_rung.append(None)
                continue
            obs_by_rung.append(pt.get('cluster', {}))

        if not all_converged:
            if verbose:
                print(f"    SKIP: not all rungs converged")
            report['observables'][eta] = {'status': 'NON-CONVERGED'}
            continue

        # Check convergence per observable
        eta_report = {}
        for obs_key, scale in SCALES.items():
            vals = [obs[obs_key] for obs in obs_by_rung]
            errors = [abs(vals[r + 1] - vals[r]) / scale
                      for r in range(len(vals) - 1)]

            ok, details = _check_contraction(errors)
            label = 'CONV' if ok else 'NON-CONV'

            if obs_key == 'splitting' and all(v < 1e-6 * scale for v in vals):
                label = 'BELOW-RESOLUTION'

            if verbose:
                errs_str = ', '.join(f'{e:.2e}' for e in errors)
                print(f"    {obs_key}: errors=[{errs_str}] → {label}")
                for d in details:
                    if 'ratio' in d:
                        sym = 'OK' if d['status'] == 'PASS' else 'FAIL'
                        print(f"      rung {d['rung']}: {d['ratio']:.2f}x {sym}")
                    elif 'final' in d:
                        sym = 'OK' if d['final_ok'] else 'FAIL'
                        print(f"      final: {d['final']:.2e} {sym}")

            eta_report[obs_key] = {'values': vals, 'errors': errors,
                                   'label': label, 'details': details}
            if not ok and label != 'BELOW-RESOLUTION':
                all_pass = False

        # Principal angles: check each angle index that exists across rungs
        min_angles = min(len(obs.get('principal_angles', [])) for obs in obs_by_rung)
        for ai in range(min_angles):
            vals = [obs['principal_angles'][ai] for obs in obs_by_rung]
            errors = [abs(vals[r + 1] - vals[r]) for r in range(len(vals) - 1)]
            ok, details = _check_contraction(errors)
            label = 'CONV' if ok else 'NON-CONV'
            key = f'angle_{ai}'
            if verbose:
                errs_str = ', '.join(f'{e:.2e}' for e in errors)
                print(f"    {key}: errors=[{errs_str}] → {label}")
            eta_report[key] = {'values': vals, 'errors': errors,
                               'label': label, 'details': details}
            if not ok:
                all_pass = False

        report['observables'][eta] = eta_report

    # Mutation arm: inject non-contracting error, verify it goes red
    if verbose:
        print(f"\n  Mutation arm: injected non-contracting sequence")
    fake_errors = [1e-4, 1e-4, 1e-4]  # constant, ratio = 1 < 3
    mut_ok, mut_details = _check_contraction(fake_errors)
    mutation_caught = not mut_ok
    report['mutation'] = {
        'injected': fake_errors,
        'caught': mutation_caught,
    }
    if verbose:
        print(f"    non-contracting detected: {mutation_caught}")
    if not mutation_caught:
        if verbose:
            print(f"    FAIL: mutation arm did not catch non-contracting sequence")
        all_pass = False

    if verbose:
        print(f"\n  Gate 9: {'PASS' if all_pass else 'FAIL'}")

    return all_pass, report


if __name__ == '__main__':
    passed, report = run_gate_7(verbose=True)
    raise SystemExit(0 if passed else 1)
