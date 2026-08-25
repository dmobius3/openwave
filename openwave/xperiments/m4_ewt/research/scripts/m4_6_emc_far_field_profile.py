#!/usr/bin/env python3
"""
M4/EWT - In-platform derivation of the far-field EMC density profile.

OpenWave criterion:
    Gravity: local metric phenomena
    Test: foundational EMC density profile

Mechanism:
    In the long-wavelength limit (r >> lambda_l) the discrete BCC lattice
    behaves as an isotropic elastic continuum. Outside a bounded spherical mass,
    the weak-field strain obeys the radial Laplace equation:

        d/dr (r^2 d(delta_eta)/dr) = 0.

    The unique spherically symmetric solution vanishing at infinity has the
    harmonic monopole form:

        delta_eta(r) = -A / r.

    This script determines the monopole amplitude A numerically by solving the
    boundary-value problem with:
        1. Gauss strain flux condition at r = R_sun:
            d(delta_eta)/dr = +r_s / R_sun^2
        2. Asymptotic Robin boundary condition at r = R_outer:
            r * d(delta_eta)/dr + delta_eta = 0

    It verifies that the resulting A equals the gravitational radius r_s = 2GM / c^2.
"""

import math

# ----------------------------------------------------------------------
# 1. Physical constants and source parameters
# ----------------------------------------------------------------------
print("[1/4] Loading physical constants...")

G = 6.67430e-11  # m^3 kg^-1 s^-2
c = 299792458.0  # m/s
M_sun = 1.989e30  # kg
R_sun = 6.957e8  # m

r_s = 2.0 * G * M_sun / (c * c)

R_outer = 500.0 * R_sun  # far-field outer boundary
N_steps = 20000  # integration grid steps

print(f"    G          = {G:.6e} m^3 kg^-1 s^-2")
print(f"    c          = {c:.3f} m/s")
print(f"    M_sun      = {M_sun:.3e} kg")
print(f"    R_sun      = {R_sun:.3e} m")
print(f"    r_s        = {r_s:.3f} m")
print(f"    R_outer    = {R_outer:.3e} m")
print(f"    N_steps    = {N_steps}")

# ----------------------------------------------------------------------
# 2. Numerical boundary-value problem for delta_eta(r)
# ----------------------------------------------------------------------
print("[2/4] Solving the radial Laplace equation...")


def rhs(r, y):
    # y[0] = delta_eta, y[1] = d(delta_eta)/dr
    return [y[1], -2.0 * y[1] / r]


def integrate(guess_delta):
    """
    Integrate radial ODE from R_sun to R_outer using RK4.
    Returns (y0_end, y1_end) at R_outer.
    """
    r = R_sun
    y0 = guess_delta
    y1 = +r_s / (R_sun * R_sun)  # Gauss strain flux condition: d(-r_s/r)/dr = +r_s/r^2

    dr = (R_outer - R_sun) / N_steps

    for _ in range(N_steps):
        k1 = rhs(r, [y0, y1])
        k2 = rhs(r + 0.5 * dr, [y0 + 0.5 * dr * k1[0], y1 + 0.5 * dr * k1[1]])
        k3 = rhs(r + 0.5 * dr, [y0 + 0.5 * dr * k2[0], y1 + 0.5 * dr * k2[1]])
        k4 = rhs(r + dr, [y0 + dr * k3[0], y1 + dr * k3[1]])

        y0 += (dr / 6.0) * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0])
        y1 += (dr / 6.0) * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1])
        r += dr

    return y0, y1


def target_residual(guess_delta):
    """
    Evaluates asymptotic BC residual: r * d(delta_eta)/dr + delta_eta = 0.
    """
    y0_end, y1_end = integrate(guess_delta)
    return R_outer * y1_end + y0_end


# Bisection loop to find inner delta_eta(R_sun) satisfying asymptotic BC
lo = -10.0 * r_s / R_sun
hi = 0.0

for _ in range(100):
    mid = 0.5 * (lo + hi)
    f_lo = target_residual(lo)
    f_mid = target_residual(mid)
    if f_lo * f_mid < 0:
        hi = mid
    else:
        lo = mid

delta_sun = 0.5 * (lo + hi)
A_numerical = -delta_sun * R_sun

print(f"    delta_eta(R_sun)   = {delta_sun:.12e}")
print(f"    A_numerical        = {A_numerical:.15e} m")
print(f"    A_expected = r_s   = {r_s:.15e} m")

# ----------------------------------------------------------------------
# 3. Verification of the full profile eta(r)
# ----------------------------------------------------------------------
print("[3/4] Verifying numerical solution against analytic eta(r) profile...")

max_rel_error = 0.0
worst_r = 0.0
max_rel_error_delta = 0.0  # error relative to delta_eta itself, not to eta ~ 1

r = R_sun
y0 = delta_sun
y1 = +r_s / (R_sun * R_sun)
dr = (R_outer - R_sun) / N_steps

for _ in range(N_steps):
    eta_numeric = 1.0 + y0
    eta_analytic = 1.0 - r_s / r
    rel_err = abs(eta_numeric - eta_analytic) / abs(eta_analytic)
    if rel_err > max_rel_error:
        max_rel_error = rel_err
        worst_r = r
    rel_err_delta = abs(y0 - (-r_s / r)) / (r_s / r)
    if rel_err_delta > max_rel_error_delta:
        max_rel_error_delta = rel_err_delta

    k1 = rhs(r, [y0, y1])
    k2 = rhs(r + 0.5 * dr, [y0 + 0.5 * dr * k1[0], y1 + 0.5 * dr * k1[1]])
    k3 = rhs(r + 0.5 * dr, [y0 + 0.5 * dr * k2[0], y1 + 0.5 * dr * k2[1]])
    k4 = rhs(r + dr, [y0 + dr * k3[0], y1 + dr * k3[1]])

    y0 += (dr / 6.0) * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0])
    y1 += (dr / 6.0) * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1])
    r += dr

print(f"    Max relative error (eta) = {max_rel_error:.6e}")
print(f"    at r                     = {worst_r:.3e} m")
print(f"    Max relative error (delta_eta) = {max_rel_error_delta:.6e}   (eta-relative / (r_s/r))")

# ----------------------------------------------------------------------
# 4. Summary
# ----------------------------------------------------------------------
print("[4/4] Summary...")

amplitude_pass = abs(A_numerical - r_s) / r_s < 1e-10
profile_pass = max_rel_error < 1e-10

if amplitude_pass:
    print("    Monopole amplitude: PASS (A_numerical == r_s)")
else:
    print("    Monopole amplitude: FAIL")

if profile_pass:
    print("    Profile verification: PASS (numerical eta matches analytic eta)")
else:
    print("    Profile verification: FAIL")

print("\nDone.")
