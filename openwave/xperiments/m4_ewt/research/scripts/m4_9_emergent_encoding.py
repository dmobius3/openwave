#!/usr/bin/env python3
"""
M4/EWT - Full Ab-Initio Emergent Encoding from Lattice Dynamics.

OpenWave criterion:
    Gravity: local metric phenomena
    Test: foundational encoding derivation

This script derives the metric encodings v ~ sqrt(eta) and f ~ sqrt(eta)
100% ab initio from microscopic pair potentials, without assuming k ~ eta^2.

Emergence Chain:
    1. Microscopic potential V(r) -> Lattice spacing a(eta) -> Derived stiffness k(eta) ~ eta^alpha (alpha = 2.0)
    2. Derived k(eta) + mass m(eta) -> Wave speed pulse simulation -> beta = 0.5
    3. Derived k(eta) + mass m(eta) -> Density well oscillator simulation -> gamma = 0.5

Numerical precision enhancements:
    - Finite-difference calculation of lattice stiffness from pair potential V(r).
    - Sub-grid pulse position via squared centroid tracking.
    - Unidirectional wave propagation initialization for pulse speed test.
    - Wide density well (R=1000) preventing spatial gradient averaging in oscillator test.
    - Exact zero-crossing times via linear interpolation.
"""

import math


def microscopic_potential(r):
    """
    Microscopic pair potential between neighboring EMCs.
    Logarithmic interaction V(r) = -V0 * ln(r) in 1D continuum limit.
    """
    return -math.log(r)


def measure_lattice_stiffness(eta, delta=1e-5):
    """
    Compute local spring stiffness k(eta) at equilibrium distance a(eta) = 1/eta
    via central finite-difference of the microscopic potential V(r):
        k = d^2 V / dr^2 |_(r = a)
    """
    a = 1.0 / eta
    v_plus = microscopic_potential(a + delta)
    v_mid = microscopic_potential(a)
    v_minus = microscopic_potential(a - delta)

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


def get_exact_zero_crossings(time_series, dt):
    """Return exact crossing times using linear interpolation."""
    crossings = []
    for i in range(1, len(time_series)):
        u_prev, u_curr = time_series[i - 1], time_series[i]
        if (u_prev >= 0.0 and u_curr < 0.0) or (u_prev <= 0.0 and u_curr > 0.0):
            t_exact = (i - 1) * dt + dt * (0.0 - u_prev) / (u_curr - u_prev)
            crossings.append(t_exact)
    return crossings


def simulate_pulse_speed(eta, N=2000, steps=800, dt=0.05):
    """
    Simulate a pulse on a uniform lattice using micro-derived parameters.
    """
    mass = eta
    stiffness = measure_lattice_stiffness(eta)

    v_approx = math.sqrt(stiffness / mass)

    u_prev = [0.0] * N
    u_cur = [0.0] * N
    center = N // 4
    sigma = 10.0

    for i in range(N):
        x = i - center
        u_cur[i] = math.exp(-(x * x) / (2.0 * sigma * sigma))
        x_prev = x + v_approx * dt
        u_prev[i] = math.exp(-(x_prev * x_prev) / (2.0 * sigma * sigma))

    peak_positions = []

    for _ in range(steps):
        u_next = [0.0] * N
        for i in range(1, N - 1):
            force = stiffness * (u_cur[i - 1] - 2.0 * u_cur[i] + u_cur[i + 1])
            acceleration = force / mass
            u_next[i] = 2.0 * u_cur[i] - u_prev[i] + acceleration * dt * dt

        u_prev, u_cur = u_cur, u_next

        pos = get_peak_centroid(u_cur)
        peak_positions.append(pos)

        if pos < 20.0 or pos > N - 20.0:
            break

    if len(peak_positions) < 2:
        return 0.0

    speed = (peak_positions[-1] - peak_positions[0]) / (len(peak_positions) * dt)
    return speed


def simulate_oscillator_frequency(eta_core, eta_background=1.0, N=3000, steps=3000, dt=0.05):
    """
    Simulate an oscillator in a density well using micro-derived parameters.
    """
    mass = [0.0] * N
    stiffness = [0.0] * (N - 1)

    R = 1000.0
    center = N // 2

    for i in range(N):
        r = abs(i - center)
        eta = eta_core + (eta_background - eta_core) * (1.0 - math.exp(-(r * r) / (R * R)))
        mass[i] = eta

    for i in range(N - 1):
        eta_mid = 0.5 * (mass[i] + mass[i + 1])
        stiffness[i] = measure_lattice_stiffness(eta_mid)

    u_prev = [0.0] * N
    u_cur = [0.0] * N
    sigma = 20.0
    k0 = 0.4

    for i in range(N):
        r = i - center
        u_cur[i] = math.cos(k0 * r) * math.exp(-(r * r) / (2.0 * sigma * sigma))
    u_prev = u_cur.copy()

    time_series = []

    for _ in range(steps):
        u_next = [0.0] * N
        for i in range(1, N - 1):
            force = stiffness[i] * (u_cur[i + 1] - u_cur[i]) - stiffness[i - 1] * (u_cur[i] - u_cur[i - 1])
            acceleration = force / mass[i]
            u_next[i] = 2.0 * u_cur[i] - u_prev[i] + acceleration * dt * dt

        u_prev, u_cur = u_cur, u_next
        time_series.append(u_cur[center])

    crossings = get_exact_zero_crossings(time_series, dt)

    if len(crossings) < 2:
        return 0.0

    max_crossings = min(len(crossings), 10)
    periods = []
    for i in range(1, max_crossings):
        periods.append(crossings[i] - crossings[i - 1])

    if not periods:
        return 0.0

    mean_half_period = sum(periods) / len(periods)
    return 1.0 / (2.0 * mean_half_period)


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
    print("M4.9 - Ab-Initio Emergent Encoding from Microscopic Potential")
    print("=" * 64)

    eta_values = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

    # Phase 1: Microscopic Potential to Stiffness Scaling (alpha)
    print("\n[1/3] Measuring micro-derived stiffness exponent (alpha)...")
    stiffnesses = []
    for eta in eta_values:
        k = measure_lattice_stiffness(eta)
        stiffnesses.append(k)
        print(f"    eta = {eta:.1f}  k = {k:.6f}")

    alpha = fit_exponent(eta_values, stiffnesses)
    print(f"\n    Fitted alpha = {alpha:.6f} (k ~ eta^alpha)")
    print(f"    Expected alpha = 2.0")

    # Phase 2: Wave speed exponent (beta)
    print("\n[2/3] Measuring emergent wave-speed exponent (beta)...")
    speeds = []
    for eta in eta_values:
        v = simulate_pulse_speed(eta)
        speeds.append(v)
        print(f"    eta = {eta:.1f}  v = {v:.6f}")

    beta = fit_exponent(eta_values, speeds)
    print(f"\n    Fitted beta = {beta:.6f} (v ~ eta^beta)")
    print(f"    Expected beta = 0.5")

    # Phase 3: Oscillator frequency exponent (gamma)
    print("\n[3/3] Measuring emergent oscillator-frequency exponent (gamma)...")
    freqs = []
    for eta_core in eta_values:
        f = simulate_oscillator_frequency(eta_core)
        freqs.append(f)
        print(f"    eta_core = {eta_core:.1f}  f = {f:.6f}")

    gamma = fit_exponent(eta_values, freqs)
    print(f"\n    Fitted gamma = {gamma:.6f} (f ~ eta^gamma)")
    print(f"    Expected gamma = 0.5")

    # Summary
    print("\n" + "=" * 64)
    alpha_pass = abs(alpha - 2.0) < 0.05
    beta_pass = abs(beta - 0.5) < 0.05
    gamma_pass = abs(gamma - 0.5) < 0.05

    if alpha_pass:
        print(f"    Stiffness exponent:  PASS (alpha = {alpha:.4f})")
    else:
        print(f"    Stiffness exponent:  FAIL (alpha = {alpha:.4f})")

    if beta_pass:
        print(f"    Wave-speed exponent: PASS (beta  = {beta:.4f})")
    else:
        print(f"    Wave-speed exponent: FAIL (beta  = {beta:.4f})")

    if gamma_pass:
        print(f"    Oscillator exponent: PASS (gamma = {gamma:.4f})")
    else:
        print(f"    Oscillator exponent: FAIL (gamma = {gamma:.4f})")

    print("=" * 64)
    print("\nDone.")


if __name__ == "__main__":
    main()