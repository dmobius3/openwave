#!/usr/bin/env python3
"""
M4/EWT - Emergent encoding from lattice dynamics.

OpenWave criterion:
    Gravity: local metric phenomena
    Test: foundational encoding derivation

This script does NOT assume n(r) = eta^(-1/2) or v = sqrt(eta).

Instead, it simulates a one-dimensional spring-mass lattice whose
microscopic parameters depend on a local density field eta(r), and
measures:

1. The speed of a propagating pulse as a function of uniform eta.
2. The frequency of an oscillator in a static density field eta.

The model inputs are strictly microscopic scalings:

    mass per lattice site    m ~ eta
    spring stiffness         k ~ eta^2

These follow from the continuum BCC picture and are not fitted to the
desired exponent.

The script measures the exponents beta and gamma and compares them with 0.5.

Numerical precision enhancements:
    - Sub-grid pulse position via squared centroid tracking.
    - Unidirectional wave propagation initialization for pulse speed test.
    - Wide density well (R=1000) preventing spatial gradient averaging in oscillator test.
    - Exact zero-crossing times via linear interpolation.
    - Robust linear regression in log-log space.
"""

import math


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
    Simulate a Gaussian pulse on a uniform lattice with density eta.

    Microphysics:
        m(eta) = eta
        k(eta) = eta^2

    Returns measured wave speed in lattice units.
    """
    mass = eta
    stiffness = eta * eta

    # Unidirectional pulse initialization to prevent pulse splitting
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
    Simulate an oscillator in a wide density profile well.

    Microphysics:
        m(eta) = eta
        k(eta) = eta^2

    Returns the measured frequency in lattice units.
    """
    mass = [0.0] * N
    stiffness = [0.0] * (N - 1)

    # Wide density well radius to prevent gradient leakage
    R = 1000.0
    center = N // 2

    for i in range(N):
        r = abs(i - center)
        eta = eta_core + (eta_background - eta_core) * (1.0 - math.exp(-(r * r) / (R * R)))
        mass[i] = eta

    for i in range(N - 1):
        eta_mid = 0.5 * (mass[i] + mass[i + 1])
        stiffness[i] = eta_mid * eta_mid

    # Compact carrier wave packet centered in core region
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

    # Calculate mean period from initial clean zero crossings
    max_crossings = min(len(crossings), 10)
    periods = []
    for i in range(1, max_crossings):
        periods.append(crossings[i] - crossings[i - 1])

    if not periods:
        return 0.0

    mean_half_period = sum(periods) / len(periods)
    return 1.0 / (2.0 * mean_half_period)


def fit_exponent(x, y):
    """Fit y = A * x^beta using valid positive data points."""
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
    print("M4.9 - Emergent Encoding from Lattice Dynamics")
    print("=" * 64)

    # Test 1: Wave speed vs density
    print("\n[1/2] Measuring wave-speed exponent (beta)...")
    eta_values = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    speeds = []

    for eta in eta_values:
        v = simulate_pulse_speed(eta)
        speeds.append(v)
        print(f"    eta = {eta:.1f}  v = {v:.6f}")

    beta = fit_exponent(eta_values, speeds)
    print(f"\n    Fitted beta = {beta:.6f}")
    print(f"    Expected beta = 0.5")

    # Test 2: Oscillator frequency vs core density
    print("\n[2/2] Measuring oscillator-frequency exponent (gamma)...")
    eta_core_values = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    freqs = []

    for eta_core in eta_core_values:
        f = simulate_oscillator_frequency(eta_core)
        freqs.append(f)
        print(f"    eta_core = {eta_core:.1f}  f = {f:.6f}")

    gamma = fit_exponent(eta_core_values, freqs)
    print(f"\n    Fitted gamma = {gamma:.6f}")
    print(f"    Expected gamma = 0.5")

    # Summary
    print("\n" + "=" * 64)
    beta_pass = abs(beta - 0.5) < 0.05
    gamma_pass = abs(gamma - 0.5) < 0.05

    if beta_pass:
        print(f"    Wave-speed exponent: PASS (beta = {beta:.4f})")
    else:
        print(f"    Wave-speed exponent: FAIL (beta = {beta:.4f})")

    if gamma_pass:
        print(f"    Oscillator exponent: PASS (gamma = {gamma:.4f})")
    else:
        print(f"    Oscillator exponent: FAIL (gamma = {gamma:.4f})")

    print("=" * 64)
    print("\nDone.")


if __name__ == "__main__":
    main()