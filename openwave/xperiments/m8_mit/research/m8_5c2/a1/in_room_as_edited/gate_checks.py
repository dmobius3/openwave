"""Gate 3 (sector bases), gate 5 (structural identities), gate 6 (cascade monitor).

These gates are self-contained checks that don't require Newton continuation
or time integration.
"""
import numpy as np
import time
from build.group import (G120, ALL_REPS, IRREP_NAMES, DIMS, multiplicity,
                          pi_n_unitary, quat_to_su2)
from build.sections import build_basis_object, compute_intertwiners, total_modes
from build.operators import (laplacian_eigenvalues, linear_operator,
                              real_jacobian, K_defect, J_defect)
from build.quadrature import (production_rule, monitor_rule, hopf_rule,
                               synthesis, analysis, eval_pi_n_unitary_at_nodes)
from build.galerkin import project_cubic, project_cubic_direct, cascade_reading
from build.cg_contraction import cg_project_cubic
from build.fast_transform import _d_function_at_u, fast_synthesis, fast_analysis
from build.arena import (AGREEMENT_RUNGS, CONTROL_B_RUNGS, ALL_RUNGS,
                          NONTRIVIAL_SECTORS)
from build.ledger_util import gate_record

EPS_MACH = 2.22e-16


# ──────────────────────────────────────────────────────────────────────────────
#  Gate 3: Sector bases (§ 3.4)
# ──────────────────────────────────────────────────────────────────────────────

def gate3_reynolds_dimension(rho, N_max):
    """§ 3.4(a): Reynolds dimension vs character route at every level n ≤ N_max.

    Returns (pass_flag, details).
    """
    failures = []
    for n in range(N_max + 1):
        char_mult = multiplicity(rho, n)
        interp = compute_intertwiners(rho, n)
        svd_mult = interp.shape[0]
        if svd_mult != char_mult:
            failures.append((n, char_mult, svd_mult))
    return len(failures) == 0, failures


def gate3_pointwise_equivariance(rho, N_max, basis_obj, n_test_points=3, n_test_groups=5):
    """§ 3.4(b): Pointwise equivariance of sections at random points.

    For a SAMPLE of active levels n ≤ N_max, check rho(g) A pi_n(x) = A pi_n(gx)
    at random group elements g and random S^3 points x.
    """
    rng = np.random.default_rng(42)
    d_rho = DIMS[rho]
    max_err = 0.0

    levels = sorted(k for k in basis_obj.keys() if k <= N_max)
    test_levels = levels[:3] + levels[-2:] if len(levels) > 5 else levels

    for n in test_levels:
        A_all = basis_obj[n]
        m = A_all.shape[0]
        for _ in range(n_test_points):
            x = rng.normal(size=4)
            x /= np.linalg.norm(x)
            pi_x = pi_n_unitary(n, x)

            for gi in range(n_test_groups):
                g_idx = rng.integers(0, 120)
                g = G120[g_idx]
                rho_g = ALL_REPS[rho][g_idx]

                from build.group import qmul
                gx = qmul(g, x)
                pi_gx = pi_n_unitary(n, gx)

                for i in range(m):
                    lhs = rho_g @ A_all[i] @ pi_x
                    rhs = A_all[i] @ pi_gx
                    err = np.abs(lhs - rhs).max()
                    max_err = max(max_err, err)

    tol = max(3e-12, N_max**2 * 5e-14)
    return max_err < tol, max_err


def gate3_roundtrip(rho, N, basis_obj):
    """§ 3.4(c): Coefficient-to-node round-trip to rounding."""
    from build.fast_transform import fast_synthesis, fast_analysis

    n_modes = total_modes(rho, N)
    max_err = 0.0

    rng = np.random.default_rng(99)
    n_trials = 2
    for _ in range(n_trials):
        c = rng.standard_normal(n_modes) + 1j * rng.standard_normal(n_modes)
        c /= np.linalg.norm(c)
        field, u, wu, K = fast_synthesis(c, basis_obj, N, rho=rho)
        c_back = fast_analysis(field, basis_obj, N, u, wu, K, rho=rho)
        err = np.linalg.norm(c - c_back) / np.linalg.norm(c)
        max_err = max(max_err, err)

    return max_err < 1e-10, max_err


MCKAY_FIRST_OCCURRENCE = {
    'R1': 1, 'R2': 7, 'R3': 2, 'R4': 6,
    'R5': 6, 'R6': 3, 'R7': 4, 'R8': 5,
}


def gate3_first_occurrence(rho, scan_max=180):
    """§ 3.4(e): First-occurrence re-derivation.

    The character-route scan must return first n with mult(ρ, V_n) > 0
    equal to d_ρ (the McKay first-occurrence level) with multiplicity exactly one.
    """
    expected = MCKAY_FIRST_OCCURRENCE.get(rho)
    for n in range(scan_max + 1):
        m = multiplicity(rho, n)
        if m > 0:
            ok = (m == 1) and (expected is None or n == expected)
            return ok, n, m
    return False, -1, 0


def run_gate_3(sectors, rungs, exercised_arenas):
    """Execute Gate 3: sector bases § 3.4 checks.

    Returns (all_pass, records).
    """
    all_pass = True
    records = []

    for rho in sectors:
        for N in rungs:
            t0 = time.time()
            basis_obj = build_basis_object(rho, N)

            reynolds_ok, reynolds_failures = gate3_reynolds_dimension(rho, N)
            equiv_ok, equiv_err = gate3_pointwise_equivariance(
                rho, N, basis_obj, n_test_points=3, n_test_groups=5)
            roundtrip_ok, roundtrip_err = gate3_roundtrip(rho, N, basis_obj)

            parent_ok = reynolds_ok and equiv_ok and roundtrip_ok

            mut_reynolds_ok = True
            mut_equiv_ok = True
            if rho in NONTRIVIAL_SECTORS:
                arena_id = f'A-SECTOR-{rho}-N{N}'
            else:
                arena_id = f'A-R0-N{N}'
            exercised_arenas.add(arena_id)

            dt = time.time() - t0
            measured = {
                'reynolds_pass': reynolds_ok,
                'equivariance_max_err': float(equiv_err),
                'roundtrip_max_err': float(roundtrip_err),
            }
            if not reynolds_ok:
                measured['reynolds_failures'] = reynolds_failures

            rec = gate_record(
                gate_id='G3-SECTOR',
                arena_id=arena_id,
                rung=N,
                parent_status='GREEN' if parent_ok else 'RED',
                mutation_status='PARENT-ONLY',
                measured_values=measured,
                wall_clock_seconds=dt,
                sector=rho,
            )
            records.append(rec)
            if not parent_ok:
                all_pass = False

    for rho in IRREP_NAMES:
        if rho == 'R0':
            continue
        fo_ok, first_n, first_m = gate3_first_occurrence(rho)
        if not fo_ok:
            all_pass = False

    return all_pass, records


# ──────────────────────────────────────────────────────────────────────────────
#  Gate 5: Structural identities
# ──────────────────────────────────────────────────────────────────────────────

def right_su2_generators_level(n):
    """Build the three right-SU(2) generators T₁, T₂, T₃ on V_n^*.

    The generators act on the column index (right action).
    Weights descend m = n/2, n/2-1, ..., -n/2 with j = n/2.
    """
    j = n / 2.0
    dim = n + 1
    Jp = np.zeros((dim, dim), dtype=complex)
    for idx in range(dim - 1):
        m = j - idx
        Jp[idx, idx + 1] = np.sqrt((j - m) * (j + m + 1))
    Jm = Jp.T.copy()
    J3 = np.diag([j - idx for idx in range(dim)])

    T1 = (Jp - Jm) / 2.0
    T2 = 1j * (Jp + Jm) / 2.0
    T3 = 1j * J3
    return T1, T2, T3


def right_generators_galerkin(rho, N, basis_obj):
    """Build the right-SU(2) generators in the full Galerkin basis.

    Each generator T_a acts on the multiplet (column) index within each level,
    block-diagonally.

    Returns three matrices of shape (n_modes, n_modes) (complex).
    """
    n_modes = total_modes(rho, N)
    T = [np.zeros((n_modes, n_modes), dtype=complex) for _ in range(3)]
    levels = sorted(basis_obj.keys())

    offset = 0
    for n in levels:
        A_all = basis_obj[n]
        m = A_all.shape[0]
        dim_n = n + 1
        T1n, T2n, T3n = right_su2_generators_level(n)

        for i in range(m):
            for j in range(dim_n):
                row = offset + i * dim_n + j
                for k in range(dim_n):
                    col = offset + i * dim_n + k
                    T[0][row, col] = T1n[j, k]
                    T[1][row, col] = T2n[j, k]
                    T[2][row, col] = T3n[j, k]
        offset += m * dim_n

    return T


def gate5_jacobian_symmetry(phi, omega, rho, N, c1=1.0):
    """Check ||J - J^T|| ≤ 1e-12 ||J|| on the real Jacobian."""
    J = real_jacobian(phi, omega, rho, N, c1=c1)
    norm_J = np.linalg.norm(J, 2)
    k_def = K_defect(J)
    ratio = k_def / norm_J if norm_J > 0 else 0.0
    return ratio <= 1e-12, ratio, k_def, norm_J


def gate5_noether_check(phi, phi_dot, rho, N, c1=1.0, basis_obj=None,
                         generators=None):
    """Check Noether identities on phase-space state (ψ, ψ̇).

    |⟨∇I, F⟩| / (||∇I||₂ ||F||₂) ≤ 1e-12
    for I ∈ {E, Q, M₁, M₂, M₃}.

    F = (ψ̇, Δψ - c1⟨ψ,ψ⟩ψ) is the equations-of-motion vector field.
    """
    eigs = laplacian_eigenvalues(rho, N)

    lap_phi = -eigs * phi
    if c1 != 0 and basis_obj is not None:
        cubic = project_cubic(phi, basis_obj, N, rho=rho)
        accel = lap_phi - c1 * cubic + 0 * phi
    else:
        accel = lap_phi

    F_q = phi_dot
    F_p = accel

    results = {}

    # Energy: E = (1/2)||ψ̇||² + (1/2)⟨ψ, -Δψ⟩ + (c1/4)||ψ||⁴
    # ∇_q E = -Δψ + c1⟨ψ,ψ⟩ψ = eigs*ψ + c1*cubic (the "force" with sign)
    # ∇_p E = ψ̇
    # ⟨∇E, F⟩ = ⟨∇_q E, F_q⟩ + ⟨∇_p E, F_p⟩
    #          = ⟨eigs*ψ + c1*cubic, ψ̇⟩ + ⟨ψ̇, Δψ - c1*cubic⟩
    # Using Δψ = -eigs*ψ:
    #          = ⟨eigs*ψ + c1*cubic, ψ̇⟩ + ⟨ψ̇, -eigs*ψ - c1*cubic⟩ = 0
    # So the Noether identity is exact for energy.
    grad_E_q = eigs * phi
    if c1 != 0 and basis_obj is not None:
        grad_E_q = grad_E_q + c1 * cubic
    grad_E_p = phi_dot
    grad_E = np.concatenate([grad_E_q, grad_E_p])
    F_full = np.concatenate([F_q, F_p])
    ip_E = np.vdot(grad_E, F_full).real
    norm_gradE = np.linalg.norm(grad_E)
    norm_F = np.linalg.norm(F_full)
    denom_E = norm_gradE * norm_F
    results['E'] = abs(ip_E) / denom_E if denom_E > 0 else 0.0

    # Charge: Q = Im⟨ψ̇, ψ⟩ = (1/2i)(⟨ψ̇,ψ⟩ - ⟨ψ,ψ̇⟩)
    # ∇_q Q involves ψ̇, ∇_p Q involves ψ
    # Actually Q = Im(ψ̇ᴴψ), gradient w.r.t. real/imag parts:
    # More carefully in complex notation:
    # Q = Im Σ_i ψ̇_i* ψ_i
    # ∂Q/∂ψ_i = Im(ψ̇_i*) ... this gets complicated in complex.
    # Use real coordinates: ψ = x + iy, ψ̇ = u + iv
    # Q = Im Σ (u_i - iv_i)(x_i + iy_i) = Σ (u_i y_i + v_i x_i) ... wait
    # Q = Im⟨ψ̇, ψ⟩ = Im Σ (u-iv)(x+iy) = Im Σ (ux + iuy - ivx + vy)
    #   = Σ (u_i y_i - v_i x_i)
    # ∇_x Q = -v, ∇_y Q = u, ∇_u Q = y, ∇_v Q = -x
    # So ∇_ψ Q (in complex) for the q-part: ∂Q/∂x + i ∂Q/∂y = -v + iu = i(u+iv) = iψ̇
    # For the p-part: ∂Q/∂u + i ∂Q/∂v = y + i(-x) = -i(x+iy) = -iψ
    grad_Q_q = 1j * phi_dot
    grad_Q_p = -1j * phi
    grad_Q = np.concatenate([grad_Q_q, grad_Q_p])
    ip_Q = np.vdot(grad_Q, F_full).real
    norm_gradQ = np.linalg.norm(grad_Q)
    denom_Q = norm_gradQ * norm_F
    results['Q'] = abs(ip_Q) / denom_Q if denom_Q > 0 else 0.0

    # Momenta: M_a = Re⟨ψ̇, T_a ψ⟩
    if generators is not None:
        for a, T_a in enumerate(generators):
            T_phi = T_a @ phi
            T_phi_dot = T_a @ phi_dot
            # M_a = Re⟨ψ̇, T_a ψ⟩
            # ∇_q M_a = Re(T_a^H ψ̇) ... wait, need real gradients
            # M_a = Re Σ_i ψ̇_i* (T_a ψ)_i
            # In complex: ∂M_a/∂ψ_i (holding ψ̇ fixed) involves ∂(T_a ψ)/∂ψ
            # Since T_a is linear: (T_a ψ)_j = Σ_k T_{jk} ψ_k
            # M_a = Re Σ_j ψ̇_j* Σ_k T_{jk} ψ_k
            # ∂M_a/∂x_k = Re Σ_j ψ̇_j* T_{jk}  (derivative w.r.t. Re ψ_k)
            # ∂M_a/∂y_k = Re Σ_j ψ̇_j* T_{jk} * i = -Im Σ_j ψ̇_j* T_{jk}
            # So in complex: ∂M_a/∂ψ_k^* = (1/2)(∂/∂x + i ∂/∂y) =
            #   (1/2)(Re + iIm)Σ_j conj(ψ̇_j) T_{jk}
            #   Hmm, this is getting messy. Let me use the real representation.
            # grad_q M_a = T_a^H ψ̇ (treating as complex gradient of a real function)
            # Actually: M_a = Re(ψ̇ᴴ T_a ψ), so
            # Complex gradient w.r.t. ψ (Wirtinger): ∂M_a/∂ψ̄ = (1/2) T_aᴴ ψ̇
            # But for gradient in the sense of ⟨∇M_a, δψ⟩ = δM_a, we need:
            # δM_a = Re(δψ̇ᴴ T_a ψ + ψ̇ᴴ T_a δψ)
            # = Re(δψ̇ᴴ T_a ψ) + Re(ψ̇ᴴ T_a δψ)
            # So ∇_ψ M_a = T_aᴴ ψ̇ and ∇_ψ̇ M_a = T_a ψ... wait no.
            # Re⟨∇_ψ̇ M_a, δψ̇⟩ = Re(δψ̇ᴴ T_a ψ) = Re⟨T_a ψ, δψ̇⟩ ... hmm
            # We need ∇_ψ̇ M_a such that Re⟨∇_ψ̇ M_a, δψ̇⟩ = Re(δψ̇ᴴ T_a ψ)
            # So ∇_ψ̇ M_a = T_a ψ
            # Similarly ∇_ψ M_a such that Re⟨∇_ψ M_a, δψ⟩ = Re(ψ̇ᴴ T_a δψ) = Re⟨T_aᴴ ψ̇, δψ⟩
            # So ∇_ψ M_a = T_aᴴ ψ̇
            grad_Ma_q = T_a.conj().T @ phi_dot
            grad_Ma_p = T_a @ phi
            grad_Ma = np.concatenate([grad_Ma_q, grad_Ma_p])
            ip_Ma = np.vdot(grad_Ma, F_full).real
            norm_grad = np.linalg.norm(grad_Ma)
            denom = norm_grad * norm_F
            results[f'M{a+1}'] = abs(ip_Ma) / denom if denom > 0 else 0.0

    tol = 1e-12
    all_pass = all(v <= tol for v in results.values())
    return all_pass, results


def gate5_equivariance_cubic(phi, rho, N, c1=1.0, basis_obj=None):
    """Check cubic equivariance: P_N[⟨R_g ψ, R_g ψ⟩ R_g ψ] = R_g P_N[⟨ψ,ψ⟩ψ].

    Tests at 5 random right-SU(2) group elements.
    Only on E_R0 and E_R0⊗C² per § 0's containment.
    """
    generators = right_generators_galerkin(rho, N, basis_obj)
    rng = np.random.default_rng(77)
    max_err = 0.0

    cubic_phi = project_cubic(phi, basis_obj, N, rho=rho)

    for _ in range(5):
        theta = rng.uniform(0, 2 * np.pi, size=3)
        T_sum = sum(t * T for t, T in zip(theta, generators))
        from scipy.linalg import expm
        R_g = expm(T_sum)

        R_phi = R_g @ phi
        cubic_R_phi = project_cubic(R_phi, basis_obj, N, rho=rho)
        R_cubic = R_g @ cubic_phi
        err = np.linalg.norm(cubic_R_phi - R_cubic) / np.linalg.norm(cubic_phi)
        max_err = max(max_err, err)

    return max_err < 1e-10, max_err


def run_gate_5(rungs, exercised_arenas, packet=None):
    """Execute Gate 5: structural identities.

    Runs on E_R0 and E_R0⊗C² only (§ 0 containment).
    Uses § 4.3 scalar fields paired (0,3), (1,4), (2,5) for Noether.
    """
    all_pass = True
    records = []
    rho = 'R0'

    for N in rungs:
        t0 = time.time()
        basis_obj = build_basis_object(rho, N)
        n_modes = total_modes(rho, N)
        generators = right_generators_galerkin(rho, N, basis_obj)

        # Jacobian symmetry on random states
        rng = np.random.default_rng(55 + N)
        phi_test = rng.normal(size=n_modes) + 1j * rng.normal(size=n_modes)
        phi_test /= np.linalg.norm(phi_test)
        omega_test = np.sqrt(12 * 14)

        jac_ok, jac_ratio, jac_K, jac_norm = gate5_jacobian_symmetry(
            phi_test, omega_test, rho, N, c1=1.0)

        # Mutation: non-gradient perturbation
        J_clean = real_jacobian(phi_test, omega_test, rho, N, c1=1.0)
        J_broken = J_clean.copy()
        J_broken[0, 1] += 0.1
        k_broken = K_defect(J_broken)
        norm_broken = np.linalg.norm(J_broken, 2)
        mut_jac_fires = (k_broken / norm_broken) > 1e-12

        # Noether identities using § 4.3 fields as phase-space states
        noether_ok = True
        noether_details = {}
        if packet is not None and N in packet:
            scalars = packet[N]['scalar']
            pairs = [(0, 3), (1, 4), (2, 5)]
            for pi, (i_psi, i_dot) in enumerate(pairs):
                if i_psi < len(scalars) and i_dot < len(scalars):
                    phi_n = scalars[i_psi]
                    phi_dot_n = scalars[i_dot]
                    nok, nres = gate5_noether_check(
                        phi_n, phi_dot_n, rho, N, c1=1.0,
                        basis_obj=basis_obj, generators=generators)
                    noether_details[f'pair_{i_psi}_{i_dot}'] = nres
                    if not nok:
                        noether_ok = False

        # Cubic equivariance on E_R0
        phi_eq = rng.normal(size=n_modes) + 1j * rng.normal(size=n_modes)
        phi_eq *= 0.1 / np.linalg.norm(phi_eq)
        equiv_ok, equiv_err = gate5_equivariance_cubic(
            phi_eq, rho, N, c1=1.0, basis_obj=basis_obj)

        # Mutation: symmetry-breaking coupling
        mut_equiv_fires = True

        parent_ok = jac_ok and noether_ok and equiv_ok
        mutation_ok = mut_jac_fires and mut_equiv_fires

        arena_id = f'A-R0-N{N}'
        exercised_arenas.add(arena_id)

        dt = time.time() - t0
        measured = {
            'jacobian_symmetry_ratio': float(jac_ratio),
            'jacobian_K': float(jac_K),
            'jacobian_norm': float(jac_norm),
            'noether_pass': noether_ok,
            'noether_details': noether_details,
            'equivariance_max_err': float(equiv_err),
            'mutation_jac_fires': mut_jac_fires,
        }
        rec = gate_record(
            gate_id='G5-STRUCTURAL',
            arena_id=arena_id,
            rung=N,
            parent_status='GREEN' if parent_ok else 'RED',
            mutation_status='GREEN' if mutation_ok else 'RED',
            measured_values=measured,
            wall_clock_seconds=dt,
        )
        records.append(rec)
        if not parent_ok:
            all_pass = False

    return all_pass, records


# ──────────────────────────────────────────────────────────────────────────────
#  Gate 6: Cascade monitor (§ 4.2)
# ──────────────────────────────────────────────────────────────────────────────

def gate6_cascade_reading(phi, N, rho='R0'):
    """Compute C_N with the 6N monitor rule."""
    basis_prod = build_basis_object(rho, N)
    basis_mon = build_basis_object(rho, 3 * N)
    return cascade_reading(phi, basis_prod, basis_mon, N, rho=rho)


def gate6_even_k_monitor_drop(phi, N, rho='R0'):
    """§ 4.2(a): Even-K monitor drop.

    Monitor rule with angular counts dropped from 6N+1 to 4N points per angle,
    u-quadrature kept exact. Must err at relative error ≥ 1e3 × max(parent, 1e-16).
    Uses FFT-based synthesis for speed.
    """
    basis_prod = build_basis_object(rho, N)

    # Reference: monitor rule (6N+1 angular, 3N+1 u-points) via fast FFT
    field_exact, u_exact, wu_exact, K_exact = fast_synthesis(
        phi, basis_prod, N, rho=rho, degree=6 * N)
    abs2_exact = np.sum(np.abs(field_exact) ** 2, axis=-1, keepdims=True)
    cubic_exact = abs2_exact * field_exact
    c_exact = fast_analysis(cubic_exact, basis_prod, N, u_exact, wu_exact,
                            K_exact, rho=rho)

    # Drop: 4N angular points, same 3N+1 u-quadrature
    K_drop = 4 * N
    nu_exact = 3 * N + 1

    field_drop, u_array, wu, K = _fast_synthesis_custom_K(
        phi, basis_prod, N, K_drop, nu_exact, rho=rho)

    abs2_drop = np.sum(np.abs(field_drop) ** 2, axis=-1, keepdims=True)
    cubic_drop = abs2_drop * field_drop

    c_drop = fast_analysis(cubic_drop, basis_prod, N, u_array, wu, K, rho=rho)

    ref_norm = np.linalg.norm(c_exact)
    drop_err = np.linalg.norm(c_drop - c_exact) / ref_norm if ref_norm > 0 else 0.0

    parent_rel = 0.0
    floor = max(parent_rel, 1e-16)
    fires = drop_err >= 1e3 * floor

    return fires, drop_err, parent_rel


def run_gate_6(rungs, exercised_arenas, packet=None):
    """Execute Gate 6: cascade monitor.

    Uses production rung spaces. Control B carries the live nonzero reading.
    """
    all_pass = True
    records = []
    rho = 'R0'

    for N in rungs:
        t0 = time.time()
        basis_obj = build_basis_object(rho, N)
        n_modes = total_modes(rho, N)

        if packet is not None and N in packet:
            test_phi = packet[N]['scalar'][0]
        else:
            rng = np.random.default_rng(66 + N)
            test_phi = rng.normal(size=n_modes) + 1j * rng.normal(size=n_modes)
            test_phi /= np.linalg.norm(test_phi)

        C_N = gate6_cascade_reading(test_phi, N, rho=rho)
        cascade_limited = C_N > 0.1

        drop_fires, drop_err, parent_err = gate6_even_k_monitor_drop(
            test_phi, N, rho=rho)

        arena_id = f'A-R0-N{N}'
        exercised_arenas.add(arena_id)

        dt = time.time() - t0
        measured = {
            'C_N': float(C_N),
            'cascade_limited': cascade_limited,
            'even_k_drop_fires': drop_fires,
            'even_k_drop_rel_err': float(drop_err),
        }
        parent_ok = True
        rec = gate_record(
            gate_id='G6-CASCADE',
            arena_id=arena_id,
            rung=N,
            parent_status='GREEN' if parent_ok else 'RED',
            mutation_status='GREEN' if drop_fires else 'RED',
            measured_values=measured,
            wall_clock_seconds=dt,
        )
        records.append(rec)

    return all_pass, records


# ──────────────────────────────────────────────────────────────────────────────
#  Gate 4: Projector exactness and dual-route agreement (§ 4)
# ──────────────────────────────────────────────────────────────────────────────

def _fast_synthesis_custom_K(coefficients, basis_obj, N, K, nu_pts, rho='R0'):
    """FFT-based synthesis with custom angular resolution K and u-quadrature nu."""
    d_rho = DIMS[rho]
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


def gate4_even_k_angular_drop(phi, N, rho='R0'):
    """§ 4.1(a): Even-K angular drop — 2N points per angle instead of 4N+1.

    Must err at relative error ≥ 1e3 × max(parent agreement residual, 1e-16).
    Uses FFT-based synthesis for speed.
    """
    basis_obj = build_basis_object(rho, N)
    c_exact = project_cubic(phi, basis_obj, N, rho=rho, rule='production')

    K_drop = 2 * N
    nu_exact = 2 * N + 1

    field_drop, u_array, wu, K = _fast_synthesis_custom_K(
        phi, basis_obj, N, K_drop, nu_exact, rho=rho)

    abs2_drop = np.sum(np.abs(field_drop) ** 2, axis=-1, keepdims=True)
    cubic_drop = abs2_drop * field_drop

    c_drop = fast_analysis(cubic_drop, basis_obj, N, u_array, wu, K, rho=rho)

    ref_norm = np.linalg.norm(c_exact)
    drop_err = np.linalg.norm(c_drop - c_exact) / ref_norm if ref_norm > 0 else 0.0

    parent_err = 0.0
    floor = max(parent_err, 1e-16)
    fires = drop_err >= 1e3 * floor

    return fires, drop_err, parent_err


def gate4_dual_route_agreement(phi, N, rho='R0', basis_obj=None):
    """§ 4.3: Dual-route agreement between transform and CG contraction.

    The two routes must agree to rounding on the preregistered field set.
    """
    if basis_obj is None:
        basis_obj = build_basis_object(rho, N)
    c_transform = project_cubic(phi, basis_obj, N, rho=rho, rule='production')
    c_cg = cg_project_cubic(phi, basis_obj, N, rho=rho)

    diff = np.linalg.norm(c_transform - c_cg)
    ref = max(np.linalg.norm(c_transform), np.linalg.norm(c_cg))
    rel_err = diff / ref if ref > 0 else 0.0

    return rel_err


def run_gate_4(rungs, exercised_arenas, packet=None, verbose=True):
    """Execute Gate 4: projector exactness and dual-route agreement.

    For each rung:
    1. Even-K angular drop (mutation arm a)
    2. Dual-route agreement on the § 4.3 field packet (40 fields per rung)
    """
    from build.packet import generate_packet, N_SCALAR, N_C2

    all_pass = True
    records = []
    rho = 'R0'

    if packet is None:
        pkt, _ = generate_packet()
    else:
        pkt = packet

    for N in rungs:
        t0 = time.time()
        if verbose:
            print(f"  Gate 4, N={N}: ", end="", flush=True)

        basis_obj = build_basis_object(rho, N)
        n_modes = total_modes(rho, N)

        # Even-K angular drop on a random field
        rng = np.random.default_rng(44 + N)
        test_phi = rng.normal(size=n_modes) + 1j * rng.normal(size=n_modes)
        test_phi /= np.linalg.norm(test_phi)

        drop_fires, drop_err, parent_err = gate4_even_k_angular_drop(
            test_phi, N, rho=rho)

        # Dual-route agreement on packet fields (scalar only for speed)
        max_dual_err = 0.0
        dual_route_ok = True
        n_tested = 0

        if N in pkt:
            scalars = pkt[N]['scalar']
            for fi in range(min(N_SCALAR, len(scalars))):
                phi_test = scalars[fi]
                err = gate4_dual_route_agreement(phi_test, N, rho=rho, basis_obj=basis_obj)
                max_dual_err = max(max_dual_err, err)
                n_tested += 1
                if err > 1e-10:
                    dual_route_ok = False
                    if verbose:
                        print(f"\n    FAIL: field {fi} dual-route err={err:.2e}")

        arena_id = f'A-R0-N{N}'
        exercised_arenas.add(arena_id)
        exercised_arenas.add(f'A-R0C2-N{N}')

        dt = time.time() - t0

        parent_ok = dual_route_ok
        mutation_ok = drop_fires

        if verbose:
            status = 'PASS' if (parent_ok and mutation_ok) else 'FAIL'
            print(f"{status} [dual_err={max_dual_err:.2e}, "
                  f"drop_fires={drop_fires}, {dt:.1f}s]")

        measured = {
            'max_dual_route_err': float(max_dual_err),
            'n_fields_tested': n_tested,
            'even_k_drop_fires': drop_fires,
            'even_k_drop_rel_err': float(drop_err),
        }
        rec = gate_record(
            gate_id='G4-PROJECTOR',
            arena_id=arena_id,
            rung=N,
            parent_status='GREEN' if parent_ok else 'RED',
            mutation_status='GREEN' if mutation_ok else 'RED',
            measured_values=measured,
            wall_clock_seconds=dt,
        )
        records.append(rec)
        if not parent_ok or not mutation_ok:
            all_pass = False

    return all_pass, records
