#!/usr/bin/env python3
"""
M4/EWT - Monopole flux from lattice geometry with derived zeta.

OpenWave criteria:
    Gravity: local metric phenomena (the flux amplitude of M4.6)
    Gravity: Newton limit (GEM), the strength (G) clause
    Test: foundational flux condition, no criterion claimed

Mechanism:
    The far-field EMC density profile is

        eta(r) = 1 - A/r,

    where A is the monopole amplitude with dimensions of length.

    This script performs two tests:

    1. Pure BCC geometry:
       N_ideal = 8*pi^4
       G_EWT_pure = f(N_ideal, L_p_geom)

    2. BCC geometry corrected by the sphere-packing impedance:
       zeta_est = (1 - eta_BCC) / (eta_BCC * N_ideal)
       N_zeta = N_ideal * (1 - zeta_est)
       G_EWT_zeta = f(N_zeta, L_p_geom)

    The second test is the physical one. It uses no alpha-calibrated
    N_final. The lattice parameters are the BCC coordination number, the
    ideal packing fraction, and the geometric lattice projection factor;
    the lattice length lambda_l is the Planck length sqrt(hbar G / c^3),
    so the measured G enters the chain once, as G_EWT ~ lambda_l^(3/2)
    ~ G^(3/4) (maintainer note in the finding).

    The monopole amplitude is then computed as

        A = 2 * G_EWT_zeta * M / c^2

    and used as the boundary condition for the far-field Laplace
    equation.
"""

import math

# ----------------------------------------------------------------------
# 1. Physical constants and EWT lattice parameters
# ----------------------------------------------------------------------
print("[1/4] Loading physical constants and EWT lattice parameters...")

c = 299792458.0  # m/s
m_e = 9.1093837015e-31  # kg
r_e = 2.8179403262e-15  # m

M_sun = 1.989e30  # kg
R_sun = 6.957e8  # m

pi = math.pi

# BCC geometry
N_ideal = 8.0 * pi**4
A_pi = 4.0 * pi**3 + pi**2 + pi
eps_M_geo = 1.0 / (8.0 * pi**7)
alpha_geo = 1.0 / (A_pi - eps_M_geo)

# Ideal BCC projection factor
L_p_geom = 2.0 / math.sqrt(3.0)

# BCC sphere packing fraction
eta_BCC = math.sqrt(3.0) * pi / 8.0

# Estimated packing impedance from the occupied-sublattice correction
zeta_est = (1.0 - eta_BCC) / (eta_BCC * N_ideal)

# Lattice stiffness corrected by the packing impedance
N_zeta = N_ideal * (1.0 - zeta_est)

# Wave center count for the electron
K_neutrinos = 10

# Background EMC density parameters
r_nu_val = 2.81794e-17
lambda_l = 1.6162e-35  # Planck length sqrt(hbar G / c^3), 5-digit value
e_euler = math.e

N_nu_statutory = (r_nu_val / (2.0 * lambda_l * e_euler)) ** 3

# Unified coupling for the geometric lattice projection
C_Unif_geo = (1.0 / K_neutrinos) + 1.0 + (alpha_geo / (pi * L_p_geom))

# Effective volume deficit
X_eff_geo = (A_pi * 3.0 * K_neutrinos * math.sqrt(2.0)) / C_Unif_geo
N_nu_effective_geo = N_nu_statutory / X_eff_geo

# Base gravitational scale from electron soliton
G_Base = (c * c * r_e) / m_e


def compute_G(N_stiffness):
    return (
        (G_Base / A_pi)
        * (1.0 / (N_stiffness * A_pi) ** 3)
        * (1.0 / (K_neutrinos * math.sqrt(N_nu_effective_geo)))
    )


G_EWT_pure = compute_G(N_ideal)
G_EWT_zeta = compute_G(N_zeta)

G_CODATA = 6.67430e-11

print(f"    c                    = {c:.3f} m/s")
print(f"    m_e                  = {m_e:.15e} kg")
print(f"    r_e                  = {r_e:.15e} m")
print(f"    M_sun                = {M_sun:.3e} kg")
print(f"    R_sun                = {R_sun:.3e} m")
print(f"    N_ideal              = {N_ideal:.15f}")
print(f"    zeta_est             = {zeta_est:.15e}")
print(f"    N_zeta               = {N_zeta:.15f}")
print(f"    alpha_geo            = {alpha_geo:.15f}")
print(f"    L_p_geom             = {L_p_geom:.15f}")

# ----------------------------------------------------------------------
# 2. Test 1: pure geometry
# ----------------------------------------------------------------------
print("\n[2/4] Test 1: pure BCC geometry...")

A_pure = 2.0 * G_EWT_pure * M_sun / (c * c)
ppm_pure = abs(G_EWT_pure - G_CODATA) / G_CODATA * 1e6

print(f"    G_EWT_pure           = {G_EWT_pure:.15e} m^3 kg^-1 s^-2")
print(f"    A_pure               = {A_pure:.15e} m")
print(f"    Difference in ppm    = {ppm_pure:.2f} ppm")

# ----------------------------------------------------------------------
# 3. Test 2: corrected geometry with zeta_est
# ----------------------------------------------------------------------
print("\n[3/4] Test 2: BCC geometry corrected by packing impedance...")

A_zeta = 2.0 * G_EWT_zeta * M_sun / (c * c)
ppm_zeta = abs(G_EWT_zeta - G_CODATA) / G_CODATA * 1e6

print(f"    G_EWT_zeta           = {G_EWT_zeta:.15e} m^3 kg^-1 s^-2")
print(f"    A_zeta               = {A_zeta:.15e} m")
print(f"    Difference in ppm    = {ppm_zeta:.2f} ppm")

# ----------------------------------------------------------------------
# 4. Solve radial Laplace equation using A_zeta
# ----------------------------------------------------------------------
print("\n[4/4] Solving radial Laplace equation with corrected geometric flux...")

R_outer = 500.0 * R_sun
N_steps = 20000


def rhs(r, y):
    return [y[1], -2.0 * y[1] / r]


dr = (R_outer - R_sun) / N_steps

y0 = -A_zeta / R_sun
y1 = +A_zeta / (R_sun * R_sun)

r = R_sun
max_rel_error_eta = 0.0
worst_r = 0.0

for _ in range(N_steps):
    eta_numeric = 1.0 + y0
    eta_analytic = 1.0 - A_zeta / r

    rel_err = abs(eta_numeric - eta_analytic) / abs(eta_analytic)
    if rel_err > max_rel_error_eta:
        max_rel_error_eta = rel_err
        worst_r = r

    k1 = rhs(r, [y0, y1])
    k2 = rhs(r + 0.5 * dr, [y0 + 0.5 * dr * k1[0], y1 + 0.5 * dr * k1[1]])
    k3 = rhs(r + 0.5 * dr, [y0 + 0.5 * dr * k2[0], y1 + 0.5 * dr * k2[1]])
    k4 = rhs(r + dr, [y0 + dr * k3[0], y1 + dr * k3[1]])

    y0 += (dr / 6.0) * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0])
    y1 += (dr / 6.0) * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1])
    r += dr

scale_factor = A_zeta / R_outer
robin_abs_residual = R_outer * y1 + y0
robin_rel_residual = abs(robin_abs_residual) / scale_factor

print(f"    Max relative error of eta  = {max_rel_error_eta:.6e}")
print(f"    at r                       = {worst_r:.3e} m")
print(f"    Robin relative residual    = {robin_rel_residual:.6e}")

# Summary
print("\nSummary...")

pure_pass = ppm_pure < 2000.0
zeta_pass = ppm_zeta < 100.0
profile_pass = max_rel_error_eta < 1e-10
robin_pass = robin_rel_residual < 1e-8

if pure_pass:
    print(f"    Pure geometry test: PASS ({ppm_pure:.2f} ppm, expected scale)")
else:
    print("    Pure geometry test: FAIL")

if zeta_pass:
    print(f"    Zeta-corrected test: PASS ({ppm_zeta:.2f} ppm)")
else:
    print("    Zeta-corrected test: FAIL")

if profile_pass:
    print("    Profile verification: PASS")
else:
    print("    Profile verification: FAIL")

if robin_pass:
    print("    Robin boundary condition: PASS")
else:
    print("    Robin boundary condition: FAIL")

print("\nDone.")
