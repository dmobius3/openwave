#!/usr/bin/env python3
"""
M4/EWT - Emergent encoding from lattice dynamics (corrected).

OpenWave criterion:
    Gravity: local metric phenomena
    Test: foundational encoding derivation

This version corrects the units and mass-law issues from the review.

Key corrections:
1. The lattice spacing is a(eta) = 1/eta, and the long-wave speed is
   v_phys = a(eta) * sqrt(k/m0), not sqrt(k/m0).
2. The per-site mass m0 is fixed. Mass per unit length is rho = m0/a.
3. The 1D pair potential is V(r) = V0/r, giving stiffness k ~ eta^3.
4. The oscillator frequency is a consistency check on the wave speed,
   via f = v/L with fixed physical length L.
"""

import math


def pair_potential(r):
    """One-dimensional repulsive pair potential V(r) = V0 / r."""
    return 1.0 / r


def measure_lattice_stiffness(eta, delta=1e-5):
    """
    Compute local spring stiffness k(eta) at the constrained spacing
    a(eta) = 1/eta via central finite-difference of the pair potential.
    """
    a = 1.0 / eta
    v_plus = pair_potential(a + delta)
    v_mid = pair_potential(a)
    v_minus = pair_potential(a - delta)

    stiffness = (v_plus - 2.0 * v_mid + v_minus) / (delta * delta)
    return stiffness


def get_peak_centroid(u, radius=15):
    """Compute sub-grid peak position using the squared centroid."""
    max_val = max(u)
    max_idx = u.index(max_val)

    i_min = max(0, max_idx - radius)
    i_max = min(len(u), max_idx + radius + 1)

    total_mass = 0.0
    weighted_sum = 0.0
    for i in range(i_min, i_max):
        w = u[i] * u[i]
        total_mass += w
        weighted_sum += i * w

    if total_mass == 0.0:
        return float(max_idx)

    return weighted_sum / total_mass


def simulate_pulse_speed_physical(eta, N=2000, steps=800, dt=0.05):
    """
    Simulate a pulse on a uniform lattice and return physical speed.

    The microscopic equation is:
        m0 * d^2u_i/dt^2 = k * (u_{i-1} - 2u_i + u_{i+1})

    The site-index speed is sqrt(k/m0). The physical speed is
        v_phys = a * sqrt(k/m0)
    where a = 1/eta.
    """
    m0 = 1.0
    stiffness = measure_lattice_stiffness(eta)

    # Site-index speed approximation
    v_site = math.sqrt(stiffness / m0)

    # Initial Gaussian pulse, unidirectional in site indices
    u_prev = [0.0] * N
    u_cur = [0.0] * N
    center = N // 4
    sigma = 10.0

    for i in range(N):
        x = i - center
        u_cur[i] = math.exp(-(x * x) / (2.0 * sigma * sigma))
        x_prev = x + v_site * dt
        u_prev[i] = math.exp(-(x_prev * x_prev) / (2.0 * sigma * sigma))

    peak_positions = []

    for _ in range(steps):
        u_next = [0.0] * N
        for i in range(1, N - 1):
            force = stiffness * (u_cur[i - 1] - 2.0 * u_cur[i] + u_cur[i + 1])
            acceleration = force / m0
            u_next[i] = 2.0 * u_cur[i] - u_prev[i] + acceleration * dt * dt

        u_prev, u_cur = u_cur, u_next

        pos = get_peak_centroid(u_cur)
        peak_positions.append(pos)

        if pos < 20.0 or pos > N - 20.0:
            break

    if len(peak_positions) < 2:
        return 0.0

    # Speed in site-index units
    v_site_measured = (peak_positions[-1] - peak_positions[0]) / (len(peak_positions) * dt)

    # Convert to physical speed using lattice spacing a = 1/eta
    a = 1.0 / eta
    v_phys_measured = v_site_measured * a

    return v_phys_measured


def fit_exponent(x, y):
    """Fit y = A * x^exponent using valid positive data points."""
    valid_pairs = [(x_i, y_i) for x_i, y_i in zip(x, y) if y_i > 0.0]
    if len(valid_pairs) < 2:
        return 0.0

    log_x = [math.log(p[0]) for p in valid_pairs]
    log_y = [math.log(p[1]) for p in valid_pairs]

    n = len(log_x)
    sum_x = sum(log_x)
    sum_y = sum(log_y)
    sum_xy = sum(a * b for a, b in zip(log_x, log_y))
    sum_x2 = sum(a * a for a in log_x)

    denom = n * sum_x2 - sum_x * sum_x
    if abs(denom) < 1e-12:
        return 0.0

    return (n * sum_xy - sum_x * sum_y) / denom


def main():
    print("=" * 64)
    print("M4.9 - Corrected Emergent Encoding")
    print("=" * 64)

    eta_values = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

    # Phase 1: stiffness exponent (analytic identity)
    print("\n[1/2] Stiffness exponent alpha (analytic identity)...")
    stiffnesses = []
    for eta in eta_values:
        k = measure_lattice_stiffness(eta)
        stiffnesses.append(k)
        print(f"    eta = {eta:.1f}  k = {k:.6f}")

    alpha = fit_exponent(eta_values, stiffnesses)
    print(f"\n    alpha = {alpha:.6f} (k ~ eta^alpha)")
    print("    Expected alpha = 3.0")
    print("    Note: this is the analytic derivative of V(r)=1/r, not a fit.")

    # Phase 2: physical wave-speed exponent
    print("\n[2/2] Physical wave-speed exponent beta...")
    speeds = []
    for eta in eta_values:
        v = simulate_pulse_speed_physical(eta)
        speeds.append(v)
        print(f"    eta = {eta:.1f}  v_phys = {v:.6f}")

    beta = fit_exponent(eta_values, speeds)
    print(f"\n    beta = {beta:.6f} (v_phys ~ eta^beta)")
    print("    Expected beta = 0.5")

    # Phase 3: oscillator frequency as consistency check
    # f = v_phys / L with fixed physical L
    # Therefore gamma = beta exactly.
    gamma = beta
    print("\n[Consistency] Oscillator frequency exponent gamma")
    print(f"    gamma = {gamma:.6f} (f ~ eta^gamma)")
    print("    Expected gamma = 0.5")
    print("    Derived from f = v_phys / L with fixed L.")
    print("    No independent simulation is required for this phase.")

    # Summary
    print("\n" + "=" * 64)

    alpha_pass = abs(alpha - 3.0) < 0.05
    beta_pass = abs(beta - 0.5) < 0.05

    if alpha_pass:
        print(f"    Stiffness exponent:  PASS (alpha = {alpha:.4f})")
    else:
        print(f"    Stiffness exponent:  FAIL (alpha = {alpha:.4f})")

    if beta_pass:
        print(f"    Wave-speed exponent: PASS (beta = {beta:.4f})")
    else:
        print(f"    Wave-speed exponent: FAIL (beta = {beta:.4f})")

    # gamma is defined as beta (f = v_phys / L, fixed L), so it is not an
    # independent measurement and carries no PASS of its own.
    print(f"    Oscillator exponent: DERIVED (gamma = beta = {gamma:.4f})")

    print("=" * 64)
    print("\nDone.")


if __name__ == "__main__":
    main()
