#!/usr/bin/env python3
"""
M4/EWT - Gravitational time dilation from the internal EMC soliton clock.

OpenWave criterion:
    Gravity: local metric phenomena
    Test: gravitational time dilation / redshift

Mechanism:
    A clock is a standing-wave soliton. Its period is the time for an
    internal longitudinal EMC lattice wave to travel from one side of
    the soliton to the opposite side and back.

    In natural lattice units (c=1), the internal wave speed is

        v_clock(r) = sqrt(eta(r)),

    where eta(r) = N_nu(r) / N_stat = 1 - r_s/r.

    The fractional frequency shift of a static clock is then

        df/f = v_clock/c - 1 = sqrt(1 - r_s/r) - 1,

    which is identically the exact Schwarzschild redshift factor for a
    static clock once eta = 1 - r_s/r is granted. Its weak-field
    expansion is -Phi_N + O(Phi_N^2), with Phi_N = GM/(c^2 r).

This script computes the redshift at the solar limb, checks the
identity against the exact GR factor, and prints the first-order
value -Phi_N as a labelled sanity line (the difference from it is the
truncation error of the first-order reference, x/4 with x = r_s/r).
"""

import math

# ----------------------------------------------------------------------
# 1. Physical constants
# ----------------------------------------------------------------------
print("[1/3] Loading physical constants...")

G = 6.67430e-11  # m^3 kg^-1 s^-2
c = 299792458.0  # m/s
M_sun = 1.989e30  # kg
R_sun = 6.957e8  # m

r_s = 2.0 * G * M_sun / (c * c)

print(f"    G          = {G:.6e} m^3 kg^-1 s^-2")
print(f"    c          = {c:.3f} m/s")
print(f"    M_sun      = {M_sun:.3e} kg")
print(f"    R_sun      = {R_sun:.3e} m")
print(f"    r_s        = {r_s:.3f} m")

# ----------------------------------------------------------------------
# 2. EMC density ratio and internal clock speed at the solar limb
# ----------------------------------------------------------------------
print("[2/3] Computing EMC density ratio and internal clock speed...")

eta = 1.0 - r_s / R_sun
v_clock_over_c = math.sqrt(eta)
Phi_N = G * M_sun / (c * c * R_sun)

print(f"    eta = N_nu/N_stat = {eta:.12f}")
print(f"    v_clock / c       = {v_clock_over_c:.12f}")
print(f"    Phi_N at limb     = {Phi_N:.6e}")

# ----------------------------------------------------------------------
# 3. Gravitational redshift
# ----------------------------------------------------------------------
print("[3/3] Computing gravitational redshift and comparing with GR...")

# EWT prediction directly from the internal clock speed
df_f_ewt = v_clock_over_c - 1.0

# Exact GR value for a static clock in the Schwarzschild field, written in
# the GR variable Phi_N (computed from G*M directly, not from r_s)
df_f_gr_exact = math.sqrt(1.0 - 2.0 * Phi_N) - 1.0

# First-order (weak-field) reference, for orientation only
df_f_first_order = -Phi_N

abs_diff_exact = abs(df_f_ewt - df_f_gr_exact)
rel_diff_first_order = abs(df_f_ewt - df_f_first_order) / abs(df_f_first_order) * 100.0
truncation_ratio = (r_s / R_sun) / 4.0 * 100.0  # x/4, the first-order truncation

print(f"    EWT predicted df/f        = {df_f_ewt:.6e}")
print(f"    GR exact df/f             = {df_f_gr_exact:.6e}")
print(f"    |EWT - GR exact|          = {abs_diff_exact:.3e}  (identity check)")
print(f"    First-order -Phi_N        = {df_f_first_order:.6e}  (sanity line)")
print(f"    Rel. diff. vs first-order = {rel_diff_first_order:.6f}%")
print(f"    Expected truncation x/4   = {truncation_ratio:.6f}%")

# The identity check: the encoding reproduces the exact factor at machine
# precision. The sanity check: the departure from the first-order reference
# is the reference's own truncation error, not a model discrepancy.
if abs_diff_exact < 1e-12 and abs(rel_diff_first_order - truncation_ratio) < 1e-6:
    print("    RESULT: PASS (exact GR factor reproduced; first-order residue = x/4)")
else:
    print("    RESULT: FAIL")

print("\nDone.")
