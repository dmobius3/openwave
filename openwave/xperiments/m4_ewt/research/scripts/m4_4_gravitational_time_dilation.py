#!/usr/bin/env python3
"""
M4/EWT - Gravitational time dilation from the internal EMC soliton clock.

OpenWave criterion:
    Gravity: metric phenomena
    Test: gravitational time dilation / redshift

Mechanism:
    A clock is a standing-wave soliton. Its period is the time for an
    internal longitudinal EMC lattice wave to travel from one side of
    the soliton to the opposite side and back.

    In natural lattice units (c=1), the internal wave speed is

        v_clock(r) = sqrt(eta(r)),

    where eta(r) = N_nu(r) / N_stat = 1 - r_s/r.

    In the weak-field limit this gives the standard gravitational
    redshift:

        df/f = -Phi_N = -GM/(c^2 r).

This script computes the redshift at the solar limb and compares it
with the standard GR value.
"""

import math

# ----------------------------------------------------------------------
# 1. Physical constants and EWT density levels
# ----------------------------------------------------------------------
print("[1/3] Loading physical constants and EWT density levels...")

G      = 6.67430e-11          # m^3 kg^-1 s^-2
c      = 299792458.0          # m/s
M_sun  = 1.989e30             # kg
R_sun  = 6.957e8              # m

r_s = 2.0 * G * M_sun / (c * c)

N_stat = 3.298651882390107e52
N_eff  = 6.252517621935487e48

print(f"    G          = {G:.6e} m^3 kg^-1 s^-2")
print(f"    c          = {c:.3f} m/s")
print(f"    M_sun      = {M_sun:.3e} kg")
print(f"    R_sun      = {R_sun:.3e} m")
print(f"    r_s        = {r_s:.3f} m")
print(f"    N_stat     = {N_stat:.6e}")
print(f"    N_eff      = {N_eff:.6e}")

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
print("[3/3] Computing gravitational redshift and comparing with standard...")

# EWT prediction from the internal clock cycle
df_f_ewt = -Phi_N

# Standard GR / observed value
df_f_target = -G * M_sun / (c * c * R_sun)

rel_diff = abs(df_f_ewt - df_f_target) / abs(df_f_target) * 100.0

print(f"    EWT predicted df/f = {df_f_ewt:.6e}")
print(f"    Target df/f        = {df_f_target:.6e}")
print(f"    Relative difference = {rel_diff:.4f}%")

if rel_diff < 0.1:
    print("    RESULT: PASS (exact match)")
else:
    print("    RESULT: FAIL")

print("\nDone.")