"""Störmer-Verlet leapfrog time integrator for gate 8 of the M8.5-C2
qualification protocol.

Implements:
  - Störmer-Verlet (second order symplectic) and symplectic Euler (first order)
  - RK45 adaptive reference integrator
  - Two controls: (i) linear standing wave on H_{R0,12}, (ii) nonlinear
    constant section on E_R0
  - Conservation checks (energy, charge, momenta)
  - Full-width wiring check
  - Gate 8 runner
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import ellipj, ellipk

from build.operators import laplacian_eigenvalues
from build.galerkin import project_cubic
from build.sections import build_basis_object, total_modes
from build.group import multiplicity


# ---------------------------------------------------------------------------
# Force computation
# ---------------------------------------------------------------------------

def _force(q, eigs, c1, basis_obj, N, rho, nodes=None, weights=None):
    """Compute F(q) = -diag(eigs)*q - c1 * P_N[|psi|^2 psi].

    Parameters
    ----------
    q : complex 1-D array  (Galerkin coefficients)
    eigs : real 1-D array  (eigenvalues n(n+2), NOT negated)
    c1 : float             (coupling constant)
    basis_obj : dict       (from build_basis_object)
    N : int                (spectral cutoff)
    rho : str              (sector name)

    Returns
    -------
    F : complex 1-D array
    """
    F = -eigs * q
    if c1 != 0.0:
        cubic = project_cubic(q, basis_obj, N, rho=rho, rule='production')
        F -= c1 * cubic
    return F


def _force_confined_i(q, eigs):
    """Force for control (i): purely linear, c1=0.
    F(q) = -diag(eigs) * q
    """
    return -eigs * q


def _force_confined_ii(q, c1):
    """Force for control (ii): single constant mode, lambda_0=0.
    F(q) = -c1 |q|^2 q
    """
    return -c1 * np.abs(q)**2 * q


# ---------------------------------------------------------------------------
# Integrators
# ---------------------------------------------------------------------------

def stormer_verlet_step(q, p, h, eigs, c1, basis_obj, N, rho,
                        nodes=None, weights=None, force_fn=None):
    """One Störmer-Verlet (leapfrog) step.

    Kick-drift-kick:
        p_{n+1/2} = p_n + (h/2) F(q_n)
        q_{n+1}   = q_n + h p_{n+1/2}
        p_{n+1}   = p_{n+1/2} + (h/2) F(q_{n+1})

    Parameters
    ----------
    q : complex array  – position (Galerkin coefficients)
    p : complex array  – velocity
    h : float          – timestep
    eigs : real array  – Laplacian eigenvalues n(n+2)
    c1 : float         – coupling constant
    basis_obj : dict   – basis object
    N : int            – spectral cutoff
    rho : str          – sector name
    force_fn : callable or None – override force; signature force_fn(q)

    Returns
    -------
    (q_new, p_new)
    """
    if force_fn is None:
        def force_fn(qq):
            return _force(qq, eigs, c1, basis_obj, N, rho, nodes, weights)

    F0 = force_fn(q)
    p_half = p + (h / 2.0) * F0          # kick
    q_new = q + h * p_half               # drift
    F1 = force_fn(q_new)
    p_new = p_half + (h / 2.0) * F1      # kick
    return q_new, p_new


def symplectic_euler_step(q, p, h, eigs, c1, basis_obj, N, rho,
                          nodes=None, weights=None, force_fn=None):
    """One kick-drift symplectic Euler step (first order).

    p_{n+1} = p_n + h F(q_n)     # kick
    q_{n+1} = q_n + h p_{n+1}    # drift
    """
    if force_fn is None:
        def force_fn(qq):
            return _force(qq, eigs, c1, basis_obj, N, rho, nodes, weights)

    F0 = force_fn(q)
    p_new = p + h * F0
    q_new = q + h * p_new
    return q_new, p_new


# ---------------------------------------------------------------------------
# RK45 reference integrator
# ---------------------------------------------------------------------------

def rk45_reference(q0, p0, T_horizon, eigs, c1, basis_obj, N, rho,
                   t_eval=None, rtol=1e-10, atol=1e-13, force_fn=None):
    """Adaptive RK45 reference integration.

    Packs (q, p) into a real state vector and uses scipy solve_ivp.

    Returns
    -------
    sol : bunch with .t, .y; also adds .q_traj, .p_traj as complex arrays.
    """
    m = len(q0)

    def rhs(t, y):
        q_re = y[:m]
        q_im = y[m:2*m]
        p_re = y[2*m:3*m]
        p_im = y[3*m:4*m]
        q = q_re + 1j * q_im
        p = p_re + 1j * p_im
        if force_fn is not None:
            F = force_fn(q)
        else:
            F = _force(q, eigs, c1, basis_obj, N, rho)
        return np.concatenate([p.real, p.imag, F.real, F.imag])

    y0 = np.concatenate([q0.real, q0.imag, p0.real, p0.imag])

    sol = solve_ivp(rhs, [0, T_horizon], y0, method='RK45',
                    t_eval=t_eval, rtol=rtol, atol=atol,
                    max_step=T_horizon / 100)

    # Unpack trajectories
    sol.q_traj = sol.y[:m] + 1j * sol.y[m:2*m]       # shape (m, n_times)
    sol.p_traj = sol.y[2*m:3*m] + 1j * sol.y[3*m:4*m]
    return sol


# ---------------------------------------------------------------------------
# Control (i): linear standing wave on H_{R0,12}
# ---------------------------------------------------------------------------

def control_i_setup():
    """Set up control (i): linear standing wave on H_{R0,12}.

    c1 = 0 (manufactured LINEAR operator)
    eigenvalue: 12*14 = 168, so omega = sqrt(168)
    The level-12 eigenspace of R0 has multiplicity 1, dimension 13.

    phi = normalized version of sum_{j=1}^{13} (1 + j/13) e^{i j/3} v_j

    Exact solution: psi(t) = cos(omega t) phi + sin(omega t)/omega phi_dot(0)
    Since p(0) = i omega phi:
        psi(t) = e^{i omega t} phi   (complex standing wave)

    Returns
    -------
    phi : complex array of length 13
    omega : float
    T_ctrl : float  (period)
    generators_T : list of 3 arrays, each (13,13)
    initial_momenta : array of length 3
    """
    n_level = 12
    dim = n_level + 1  # = 13
    lam = n_level * (n_level + 2)  # = 168
    omega = np.sqrt(lam)

    # Check multiplicity
    m = multiplicity('R0', n_level)
    assert m == 1, f"expected multiplicity 1 at level 12, got {m}"

    # Build initial data: phi = normalized sum of (1+j/13) e^{ij/3} v_j
    # j runs from 1 to 13, indexing basis vectors v_1 ... v_13
    # In our 0-based array: index j-1 for j=1..13
    phi_raw = np.zeros(dim, dtype=complex)
    for j in range(1, dim + 1):
        phi_raw[j - 1] = (1.0 + j / 13.0) * np.exp(1j * j / 3.0)

    # Normalize to unit L2
    phi = phi_raw / np.linalg.norm(phi_raw)

    # Period
    T_ctrl = 2.0 * np.pi / omega

    # Initial velocity for complex standing wave: p(0) = i omega phi
    p0 = 1j * omega * phi

    # Generators for SU(2) spin-6 representation (dim 13)
    # Canonical multiplet basis carries weights DESCENDING: m = +6, +5, ..., -6
    j_spin = 6.0
    m_vals = np.arange(j_spin, -j_spin - 1, -1)  # [+6, +5, ..., -6]

    # J_+ raises m by 1: J_+ |m> = sqrt((j-m)(j+m+1)) |m+1>
    # In our descending basis, |m+1> is one index BEFORE |m>
    # So J_+[k-1, k] = sqrt((j - m_k)(j + m_k + 1)) where m_k = m_vals[k]
    J_plus = np.zeros((dim, dim), dtype=complex)
    for k in range(1, dim):
        mk = m_vals[k]
        coeff = np.sqrt((j_spin - mk) * (j_spin + mk + 1))
        J_plus[k - 1, k] = coeff

    J_minus = J_plus.T.copy()  # J_- = J_+^T (real matrix transposed)
    J3 = np.diag(m_vals)

    # Generators T_a (anti-Hermitian):
    # T_1 = (J_+ - J_-) / 2
    # T_2 = i(J_+ + J_-) / 2
    # T_3 = i J_3
    T1 = (J_plus - J_minus) / 2.0
    T2 = 1j * (J_plus + J_minus) / 2.0
    T3 = 1j * J3

    generators_T = [T1, T2, T3]

    # Momenta: M_a(psi, psi_dot) = Re<psi_dot, T_a psi>
    initial_momenta = np.zeros(3)
    for a in range(3):
        initial_momenta[a] = np.real(np.vdot(p0, generators_T[a] @ phi))

    # Eigenvalues for this confined system (all the same: 168)
    eigs_confined = np.full(dim, float(lam))

    # --- Setup assertions ---

    # Scale for conservation
    norm_phi = np.linalg.norm(phi)
    norm_p0 = np.linalg.norm(p0)
    S_cons = norm_phi * norm_p0 + norm_phi**2

    # All three |M_a| >= 1e-6 * S_cons
    for a in range(3):
        assert abs(initial_momenta[a]) >= 1e-6 * S_cons, (
            f"|M_{a+1}| = {abs(initial_momenta[a]):.2e} < 1e-6 * S_cons = {1e-6 * S_cons:.2e}"
        )

    # |M_1/M_2| = tan(1/3) to 1e-12
    ratio = abs(initial_momenta[0] / initial_momenta[1])
    expected_ratio = abs(np.tan(1.0 / 3.0))
    assert abs(ratio - expected_ratio) < 1e-12, (
        f"|M_1/M_2| = {ratio:.15e}, expected tan(1/3) = {expected_ratio:.15e}, "
        f"diff = {abs(ratio - expected_ratio):.2e}"
    )

    # M_3 sign check: should be negative (see docstring reasoning)
    # M_3 = Re<i omega phi, i J_3 phi> = omega * Re(phi^H J_3 phi)
    # = omega * sum_j m_j |c_j|^2
    # Descending m_j, ascending |c_j| -> sum is negative -> M_3 < 0
    phiH_J3_phi = np.vdot(phi, J3 @ phi).real
    assert phiH_J3_phi < 0, f"phi^H J_3 phi = {phiH_J3_phi}, expected negative"
    expected_M3 = omega * phiH_J3_phi
    assert abs(initial_momenta[2] - expected_M3) < 1e-12 * abs(expected_M3), (
        f"M_3 = {initial_momenta[2]:.15e}, expected {expected_M3:.15e}"
    )

    return phi, omega, T_ctrl, generators_T, initial_momenta, eigs_confined, p0


# ---------------------------------------------------------------------------
# Control (ii): nonlinear constant section on E_R0
# ---------------------------------------------------------------------------

def control_ii_setup():
    """Set up control (ii): nonlinear constant section on E_R0.

    c1 = +1, real initial data q(0) = 1, q_dot(0) = 0.
    The constant mode (level 0) has eigenvalue lambda_0 = 0.

    Exact solution: q(t) = q0 * cn(omega*t, k) where
      omega = sqrt(c1) * q0 = 1
      k^2 = c1 * q0^2 / (2 omega^2) = 1/2
      k = 1/sqrt(2)

    Period: T = 4 K(k^2) / omega = 4 K(0.5)

    Returns
    -------
    q0_val : complex array of length 1
    p0_val : complex array of length 1
    T_ctrl : float
    omega : float
    k_sq : float  (= 0.5, the parameter m for ellipj)
    eigs_confined : real array of length 1
    """
    c1 = 1.0
    q0_scalar = 1.0
    omega = np.sqrt(c1) * q0_scalar  # = 1.0
    k_sq = 0.5  # m = k^2

    K_val = ellipk(k_sq)
    T_ctrl = 4.0 * K_val / omega

    q0_val = np.array([q0_scalar + 0j])
    p0_val = np.array([0.0 + 0j])
    eigs_confined = np.array([0.0])  # level 0: lambda = 0*2 = 0

    return q0_val, p0_val, T_ctrl, omega, k_sq, eigs_confined, c1


def control_ii_exact(t, omega, k_sq, q0_scalar=1.0):
    """Exact solution for control (ii): q(t) = q0 * cn(omega*t, k^2)."""
    sn, cn, dn, ph = ellipj(omega * np.atleast_1d(t), k_sq)
    return q0_scalar * cn


# ---------------------------------------------------------------------------
# Error metric
# ---------------------------------------------------------------------------

def error_metric(q_h, q_ref, t_eval=None):
    """Compute relative supremum error.

    error = max_t ||q_h(t) - q_ref(t)||_2 / max_t ||q_ref(t)||_2

    Parameters
    ----------
    q_h : array, shape (m, n_times) or (n_times,) for scalar
    q_ref : array, shape (m, n_times) or (n_times,)
    """
    q_h = np.atleast_2d(q_h)
    q_ref = np.atleast_2d(q_ref)
    diff_norms = np.sqrt(np.sum(np.abs(q_h - q_ref)**2, axis=0))
    ref_norms = np.sqrt(np.sum(np.abs(q_ref)**2, axis=0))
    max_ref = np.max(ref_norms)
    if max_ref < 1e-30:
        return 0.0
    return float(np.max(diff_norms) / max_ref)


# ---------------------------------------------------------------------------
# Conservation checks
# ---------------------------------------------------------------------------

def _energy(q, p, eigs, c1):
    """Hamiltonian: H = (1/2)||p||^2 + (1/2) sum_i lambda_i |q_i|^2
       + (c1/4) ||q||^4  (for the confined constant-mode case with lambda=0,
       this simplifies to (1/2)|p|^2 + (c1/4)|q|^4).

    More precisely, the energy for the second-order system
    q_tt = -lambda q - c1 |q|^2 q is:
    H = (1/2)||p||^2 + (1/2) q^H diag(lambda) q + (c1/4)(q^H q)^2
    """
    kinetic = 0.5 * np.vdot(p, p).real
    potential_linear = 0.5 * np.sum(eigs * np.abs(q)**2)
    norm_sq = np.vdot(q, q).real
    potential_nonlinear = (c1 / 4.0) * norm_sq**2
    return kinetic + potential_linear + potential_nonlinear


def _charge(q, p):
    """Charge (U(1) Noether): Q = Im<q, p>."""
    return np.imag(np.vdot(q, p))


def conservation_check(q_traj, p_traj, dt, T_ctrl, c1, eigs,
                        generators=None, omega_ctrl=None):
    """Check conservation of energy, charge, and momenta.

    Parameters
    ----------
    q_traj : list of complex arrays (one per timestep)
    p_traj : list of complex arrays
    dt : float
    T_ctrl : float  (characteristic period)
    c1 : float
    eigs : real array
    generators : list of (dim, dim) arrays or None
    omega_ctrl : float or None (for energy envelope bound)

    Returns
    -------
    results : dict with conservation data
    """
    n_steps = len(q_traj)

    # Scale
    norms_q = [np.linalg.norm(q) for q in q_traj]
    norms_p = [np.linalg.norm(p) for p in p_traj]
    S_cons = max(nq * np_val + nq**2 for nq, np_val in zip(norms_q, norms_p))

    # Energy
    energies = np.array([_energy(q, p, eigs, c1)
                         for q, p in zip(q_traj, p_traj)])
    E_mean = np.mean(energies)
    E_dev = energies - E_mean
    E_envelope = np.max(np.abs(E_dev))

    # Charge
    charges = np.array([_charge(q, p) for q, p in zip(q_traj, p_traj)])
    charge_drift = np.max(np.abs(charges - charges[0]))

    # Momenta
    momenta_drift = None
    if generators is not None:
        n_gen = len(generators)
        momenta = np.zeros((n_steps, n_gen))
        for t_idx in range(n_steps):
            for a in range(n_gen):
                momenta[t_idx, a] = np.real(
                    np.vdot(p_traj[t_idx], generators[a] @ q_traj[t_idx])
                )
        momenta_drift = np.max(np.abs(momenta - momenta[0:1, :]))

    results = {
        'S_cons': S_cons,
        'energy_mean': E_mean,
        'energy_envelope': E_envelope,
        'energy_relative_envelope': E_envelope / abs(E_mean) if abs(E_mean) > 1e-30 else 0.0,
        'charge_drift': charge_drift,
        'charge_relative': charge_drift / S_cons if S_cons > 1e-30 else 0.0,
    }

    if momenta_drift is not None:
        results['momenta_drift'] = momenta_drift
        results['momenta_relative'] = momenta_drift / S_cons if S_cons > 1e-30 else 0.0

    # Energy envelope bound: 5(omega h)^2 / 8 relative to mean
    if omega_ctrl is not None and abs(E_mean) > 1e-30:
        h_r = dt  # timestep used
        expected_envelope = 5.0 * (omega_ctrl * h_r)**2 / 8.0
        results['energy_envelope_bound'] = expected_envelope
        results['energy_envelope_ratio'] = E_envelope / abs(E_mean)

        # Secular test: linear fit slope * horizon <= 10% of envelope
        t_arr = np.arange(n_steps) * dt
        if n_steps > 2:
            slope = np.polyfit(t_arr, energies, 1)[0]
            horizon = t_arr[-1]
            secular = abs(slope) * horizon
            results['energy_secular'] = secular
            results['energy_secular_fraction'] = (
                secular / E_envelope if E_envelope > 1e-30 else 0.0
            )

    return results


# ---------------------------------------------------------------------------
# Full-width wiring check
# ---------------------------------------------------------------------------

def wiring_check(control_id, N=36):
    """Full-width wiring check: 1000 steps at N=36.

    Embeds the confined control into the full E_{R0} space at cutoff N,
    runs 1000 Störmer-Verlet steps, and checks that the full-width
    evolution matches the confined evolution at every step.

    Parameters
    ----------
    control_id : int (1 or 2)
    N : int (default 36)

    Returns
    -------
    max_relative_error : float
    """
    rho = 'R0'
    M = total_modes(rho, N)

    # Build basis object for the full system
    basis_obj = build_basis_object(rho, N)

    # Full eigenvalue array
    eigs_full = laplacian_eigenvalues(rho, N)

    # Determine T_min from the highest active eigenvalue
    # lambda_max at level N: need highest level with nonzero multiplicity
    max_level = 0
    for n in range(N + 1):
        if multiplicity(rho, n) > 0:
            max_level = n
    lam_max = max_level * (max_level + 2)
    omega_max = np.sqrt(lam_max)
    T_min = 2.0 * np.pi / omega_max
    h_wire = T_min / 100.0
    n_steps = 1000

    if control_id == 1:
        # Control (i): linear standing wave at level 12
        phi, omega, T_ctrl, generators_T, _, eigs_confined, p0_confined = control_i_setup()
        c1 = 0.0
        q_confined = phi.copy()
        p_confined = p0_confined.copy()

        # Find the offset of level-12 modes in the full array
        offset = 0
        for n in sorted(basis_obj.keys()):
            if n == 12:
                break
            m_n = basis_obj[n].shape[0]
            offset += m_n * (n + 1)
        n12_modes = multiplicity(rho, 12) * 13

        # Embed into full space
        q_full = np.zeros(M, dtype=complex)
        p_full = np.zeros(M, dtype=complex)
        q_full[offset:offset + n12_modes] = q_confined
        p_full[offset:offset + n12_modes] = p_confined

        def force_confined(q):
            return _force_confined_i(q, eigs_confined)

        def force_full(q):
            return _force(q, eigs_full, c1, basis_obj, N, rho)

    elif control_id == 2:
        # Control (ii): constant mode at level 0
        q0_val, p0_val, T_ctrl, omega_ctrl, k_sq, eigs_confined, c1 = control_ii_setup()
        q_confined = q0_val.copy()
        p_confined = p0_val.copy()

        # Level 0 is the first mode
        q_full = np.zeros(M, dtype=complex)
        p_full = np.zeros(M, dtype=complex)
        q_full[0] = q_confined[0]
        p_full[0] = p_confined[0]
        offset = 0
        n_confined = 1

        def force_confined(q):
            return _force_confined_ii(q, c1)

        def force_full(q):
            return _force(q, eigs_full, c1, basis_obj, N, rho)

    else:
        raise ValueError(f"unknown control_id: {control_id}")

    # Run both confined and full-width integrations
    max_rel_err = 0.0
    for step in range(n_steps):
        # Confined step
        q_confined, p_confined = stormer_verlet_step(
            q_confined, p_confined, h_wire, eigs_confined, c1,
            None, None, rho, force_fn=force_confined
        )
        # Full-width step
        q_full, p_full = stormer_verlet_step(
            q_full, p_full, h_wire, eigs_full, c1,
            basis_obj, N, rho, force_fn=force_full
        )

        # Extract the confined subspace from full
        if control_id == 1:
            q_sub = q_full[offset:offset + n12_modes]
        else:
            q_sub = q_full[offset:offset + 1]

        # Check match
        diff = np.linalg.norm(q_sub - q_confined)
        ref_norm = np.linalg.norm(q_confined)
        if ref_norm > 1e-30:
            rel_err = diff / ref_norm
        else:
            rel_err = diff
        max_rel_err = max(max_rel_err, rel_err)

        # Also check that modes outside the confined subspace stay zero
        q_outside = q_full.copy()
        if control_id == 1:
            q_outside[offset:offset + n12_modes] = 0.0
        else:
            q_outside[0] = 0.0
        outside_norm = np.linalg.norm(q_outside)
        if ref_norm > 1e-30:
            max_rel_err = max(max_rel_err, outside_norm / ref_norm)

    return max_rel_err


# ---------------------------------------------------------------------------
# Run a timestep ladder
# ---------------------------------------------------------------------------

def _run_ladder(q0, p0, T_horizon, dt_base, eigs, c1, force_fn,
                t_eval_ref, integrator='sv', n_rungs=4):
    """Run a timestep ladder {h, h/2, h/4, h/8} and collect trajectories.

    Returns list of (dt, q_at_eval) for each rung, where q_at_eval has
    shape (m, len(t_eval_ref)) sampled at the reference times.
    """
    results = []
    for rung in range(n_rungs):
        dt = dt_base / (2**rung)
        n_steps = int(round(T_horizon / dt))

        # Storage at eval times
        # We sample at multiples of the base dt that are also eval points
        eval_set = set(np.round(t_eval_ref / dt).astype(int))
        eval_set = {k for k in eval_set if 0 <= k <= n_steps}

        q = q0.copy()
        p = p0.copy()
        m = len(q0)

        # Collect trajectory at every step for simplicity
        q_list = [q.copy()]
        for step in range(1, n_steps + 1):
            if integrator == 'sv':
                q, p = stormer_verlet_step(
                    q, p, dt, eigs, c1, None, None, 'R0', force_fn=force_fn
                )
            elif integrator == 'euler':
                q, p = symplectic_euler_step(
                    q, p, dt, eigs, c1, None, None, 'R0', force_fn=force_fn
                )
            q_list.append(q.copy())

        # Interpolate to t_eval_ref times
        t_all = np.arange(n_steps + 1) * dt
        q_arr = np.array(q_list)  # shape (n_steps+1, m)

        # Find closest step for each eval time
        q_eval = np.zeros((m, len(t_eval_ref)), dtype=complex)
        for i, te in enumerate(t_eval_ref):
            idx = int(round(te / dt))
            idx = min(idx, n_steps)
            q_eval[:, i] = q_arr[idx]

        results.append((dt, q_eval))

    return results


# ---------------------------------------------------------------------------
# Gate 8 runner
# ---------------------------------------------------------------------------

def run_gate_8(verbose=True):
    """Execute gate 8: the supportive time arm.

    For each control:
      1. Set up the control
      2. Compute RK45 reference
      3. Run Störmer-Verlet at dt ladder, horizon 10*T_ctrl
      4. Compute error metric at each rung
      5. Check error <= 1e-3 at finest rung
      6. Check contraction >= 3x per rung
      7. Check conservation
      8. Run full-width wiring check

    Mutation arm: symplectic Euler on control (ii) must FAIL 3x contraction.

    Returns
    -------
    passed : bool
    report : dict
    """
    report = {}
    all_pass = True

    # ======================================================================
    # CONTROL (i): linear standing wave
    # ======================================================================
    if verbose:
        print("=" * 60)
        print("CONTROL (i): linear standing wave on H_{R0,12}")
        print("=" * 60)

    phi, omega_i, T_i, generators_T, momenta_i, eigs_i, p0_i = control_i_setup()
    if verbose:
        print(f"  omega = sqrt(168) = {omega_i:.6f}")
        print(f"  T = {T_i:.6f}")
        print(f"  momenta = {momenta_i}")

    c1_i = 0.0
    dim_i = len(phi)
    T_horizon_i = 10.0 * T_i

    def force_i(q):
        return _force_confined_i(q, eigs_i)

    # Exact solution: psi(t) = e^{i omega t} phi
    def exact_i(t_arr):
        """Returns shape (dim, n_times)"""
        return phi[:, None] * np.exp(1j * omega_i * t_arr[None, :])

    # Base timestep: T/100 per § 7
    dt_base_i = T_i / 100.0
    n_eval = 201
    t_eval_i = np.linspace(0, T_horizon_i, n_eval)

    # Reference: exact solution (no need for RK45 on linear system)
    q_ref_i = exact_i(t_eval_i)

    # Ladder
    ladder_i = _run_ladder(phi, p0_i, T_horizon_i, dt_base_i, eigs_i,
                           c1_i, force_i, t_eval_i, integrator='sv', n_rungs=4)

    errors_i = []
    for dt_r, q_r in ladder_i:
        err = error_metric(q_r, q_ref_i)
        errors_i.append(err)
        if verbose:
            print(f"  dt={dt_r:.6e}  error={err:.6e}")

    # Check finest rung error
    if errors_i[-1] > 1e-3:
        if verbose:
            print(f"  FAIL: finest rung error {errors_i[-1]:.2e} > 1e-3")
        all_pass = False
    else:
        if verbose:
            print(f"  PASS: finest rung error {errors_i[-1]:.2e} <= 1e-3")

    # Check contraction >= 3x
    contractions_i = []
    for k in range(1, len(errors_i)):
        if errors_i[k] > 1e-30:
            c = errors_i[k - 1] / errors_i[k]
        else:
            c = float('inf')
        contractions_i.append(c)
        if verbose:
            print(f"  contraction rung {k}: {c:.2f}x")
        if c < 3.0:
            if verbose:
                print(f"  FAIL: contraction {c:.2f} < 3.0")
            all_pass = False

    # Conservation check
    # Run a trajectory storing all steps at finest dt
    dt_finest_i = dt_base_i / 8.0
    n_steps_i = int(round(T_horizon_i / dt_finest_i))
    q_traj_i = [phi.copy()]
    p_traj_i = [p0_i.copy()]
    q_c, p_c = phi.copy(), p0_i.copy()
    for _ in range(n_steps_i):
        q_c, p_c = stormer_verlet_step(
            q_c, p_c, dt_finest_i, eigs_i, c1_i, None, None, 'R0',
            force_fn=force_i
        )
        q_traj_i.append(q_c.copy())
        p_traj_i.append(p_c.copy())

    cons_i = conservation_check(q_traj_i, p_traj_i, dt_finest_i, T_i,
                                 c1_i, eigs_i, generators=generators_T,
                                 omega_ctrl=omega_i)
    if verbose:
        print(f"  conservation: charge_rel={cons_i['charge_relative']:.2e}, "
              f"momenta_rel={cons_i.get('momenta_relative', 'N/A')}")

    if cons_i['charge_relative'] > 1e-12:
        if verbose:
            print(f"  FAIL: charge drift {cons_i['charge_relative']:.2e} > 1e-12")
        all_pass = False

    if cons_i.get('momenta_relative', 0) > 1e-12:
        if verbose:
            print(f"  FAIL: momenta drift {cons_i['momenta_relative']:.2e} > 1e-12")
        all_pass = False

    report['control_i'] = {
        'errors': errors_i,
        'contractions': contractions_i,
        'conservation': cons_i,
    }

    # ======================================================================
    # CONTROL (ii): nonlinear constant section
    # ======================================================================
    if verbose:
        print()
        print("=" * 60)
        print("CONTROL (ii): nonlinear constant section on E_R0")
        print("=" * 60)

    q0_ii, p0_ii, T_ii, omega_ii, k_sq_ii, eigs_ii, c1_ii = control_ii_setup()
    if verbose:
        print(f"  omega = {omega_ii:.6f}, k^2 = {k_sq_ii}")
        print(f"  T_cn = {T_ii:.6f}")

    T_horizon_ii = 10.0 * T_ii

    def force_ii(q):
        return _force_confined_ii(q, c1_ii)

    dt_base_ii = T_ii / 100.0
    n_eval_ii = 201
    t_eval_ii = np.linspace(0, T_horizon_ii, n_eval_ii)

    # RK45 reference
    sol_ii = rk45_reference(q0_ii, p0_ii, T_horizon_ii, eigs_ii, c1_ii,
                            None, None, 'R0', t_eval=t_eval_ii,
                            force_fn=force_ii)
    q_ref_ii = sol_ii.q_traj  # shape (1, n_eval)

    # Also check against exact cn solution
    q_exact_ii = control_ii_exact(t_eval_ii, omega_ii, k_sq_ii)
    exact_err = np.max(np.abs(q_ref_ii[0] - q_exact_ii))
    if verbose:
        print(f"  RK45 vs exact cn: max|diff| = {exact_err:.2e}")
    assert exact_err < 1e-8, f"RK45 vs cn mismatch: {exact_err}"

    # Ladder
    ladder_ii = _run_ladder(q0_ii, p0_ii, T_horizon_ii, dt_base_ii, eigs_ii,
                            c1_ii, force_ii, t_eval_ii, integrator='sv', n_rungs=4)

    errors_ii = []
    for dt_r, q_r in ladder_ii:
        err = error_metric(q_r, q_ref_ii)
        errors_ii.append(err)
        if verbose:
            print(f"  dt={dt_r:.6e}  error={err:.6e}")

    if errors_ii[-1] > 1e-3:
        if verbose:
            print(f"  FAIL: finest rung error {errors_ii[-1]:.2e} > 1e-3")
        all_pass = False
    else:
        if verbose:
            print(f"  PASS: finest rung error {errors_ii[-1]:.2e} <= 1e-3")

    contractions_ii = []
    for k in range(1, len(errors_ii)):
        if errors_ii[k] > 1e-30:
            c = errors_ii[k - 1] / errors_ii[k]
        else:
            c = float('inf')
        contractions_ii.append(c)
        if verbose:
            print(f"  contraction rung {k}: {c:.2f}x")
        if c < 3.0:
            if verbose:
                print(f"  FAIL: contraction {c:.2f} < 3.0")
            all_pass = False

    # Conservation for control (ii)
    dt_finest_ii = dt_base_ii / 8.0
    n_steps_ii = int(round(T_horizon_ii / dt_finest_ii))
    q_traj_ii = [q0_ii.copy()]
    p_traj_ii = [p0_ii.copy()]
    q_c, p_c = q0_ii.copy(), p0_ii.copy()
    for _ in range(n_steps_ii):
        q_c, p_c = stormer_verlet_step(
            q_c, p_c, dt_finest_ii, eigs_ii, c1_ii, None, None, 'R0',
            force_fn=force_ii
        )
        q_traj_ii.append(q_c.copy())
        p_traj_ii.append(p_c.copy())

    cons_ii = conservation_check(q_traj_ii, p_traj_ii, dt_finest_ii, T_ii,
                                  c1_ii, eigs_ii, omega_ctrl=omega_ii)
    if verbose:
        print(f"  conservation: charge_rel={cons_ii['charge_relative']:.2e}, "
              f"energy_envelope={cons_ii['energy_envelope']:.2e}")

    report['control_ii'] = {
        'errors': errors_ii,
        'contractions': contractions_ii,
        'conservation': cons_ii,
    }

    # ======================================================================
    # MUTATION ARM: symplectic Euler on control (ii) must FAIL 3x contraction
    # ======================================================================
    if verbose:
        print()
        print("=" * 60)
        print("MUTATION ARM: symplectic Euler on control (ii)")
        print("=" * 60)

    ladder_euler = _run_ladder(q0_ii, p0_ii, T_horizon_ii, dt_base_ii, eigs_ii,
                               c1_ii, force_ii, t_eval_ii, integrator='euler',
                               n_rungs=4)

    errors_euler = []
    for dt_r, q_r in ladder_euler:
        err = error_metric(q_r, q_ref_ii)
        errors_euler.append(err)
        if verbose:
            print(f"  dt={dt_r:.6e}  error={err:.6e}")

    # Check that contraction FAILS to reach 3x at some rung
    euler_has_bad_contraction = False
    for k in range(1, len(errors_euler)):
        if errors_euler[k] > 1e-30:
            c = errors_euler[k - 1] / errors_euler[k]
        else:
            c = float('inf')
        if verbose:
            print(f"  contraction rung {k}: {c:.2f}x")
        if c < 3.0:
            euler_has_bad_contraction = True

    if not euler_has_bad_contraction:
        if verbose:
            print("  FAIL: symplectic Euler achieved 3x contraction (mutation not caught)")
        all_pass = False
    else:
        if verbose:
            print("  PASS: symplectic Euler fails 3x contraction (mutation detected)")

    report['mutation_arm'] = {
        'errors': errors_euler,
        'euler_fails_contraction': euler_has_bad_contraction,
    }

    # ======================================================================
    # WIRING CHECKS (skip if verbose=False for speed; always run in gate)
    # ======================================================================
    if verbose:
        print()
        print("=" * 60)
        print("WIRING CHECKS (N=36)")
        print("=" * 60)

    for ctrl in [1, 2]:
        if verbose:
            print(f"  Control ({ctrl}): ", end="", flush=True)
        try:
            wire_err = wiring_check(ctrl, N=36)
            if verbose:
                print(f"max_rel_error = {wire_err:.2e}")
            if wire_err > 1e-10:
                if verbose:
                    print(f"  FAIL: wiring error {wire_err:.2e} > 1e-10")
                all_pass = False
            else:
                if verbose:
                    print(f"  PASS")
        except Exception as e:
            if verbose:
                print(f"ERROR: {e}")
            all_pass = False
            report[f'wiring_{ctrl}'] = str(e)

    # ======================================================================
    # FINAL VERDICT
    # ======================================================================
    if verbose:
        print()
        print("=" * 60)
        if all_pass:
            print("GATE 8: PASS")
        else:
            print("GATE 8: FAIL")
        print("=" * 60)

    report['passed'] = all_pass
    return all_pass, report


# ---------------------------------------------------------------------------
# Standalone entry point
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    passed, report = run_gate_8(verbose=True)
    raise SystemExit(0 if passed else 1)
