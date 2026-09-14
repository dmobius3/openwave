#!/usr/bin/env python3
"""M4/EWT - Newtonian Force from EMC Push-Out with Fully Geometric G

OpenWave criterion:
    Gravity: Newton limit (GEM)

Purpose:
    Verify that the EMC field-overlap geometry recovers Newton's 1/r^2
    force law when the monopole amplitude A and the coupling K_emc are
    both built from the self-consistent geometric trinity G_geom,
    hbar_geom, lambda_l derived in the M4.7 emergence engine.

    The strength clause:

      circularity : discharged. G_geom is derived from BCC geometry
                    without re-entering through an input.
      accuracy    : stated and not met. The G_geom residual is 21.9x the
                    CODATA uncertainty on G.

Dimensional anchors:
    r_e, m_e, c are the measured anchors used to build the dimensionless
    geometric ratio. G_geom is a dimensionless ratio against c^2 r_e / m_e,
    not G from nothing. The anchors are named at the point of use.

Physics & Structural Identification:
1. Exact Domain Integration: mapped radial integration r = R/(1-t) over
   t in [0, 1) evaluates the pre-computed angular result J(r) = 4pi/r^2
   for r >= R.
2. Algebraic Identity: A = 2 G M / c^2 and K_emc = c^4 / (16 pi G) use
   the same G_geom, so G_geom cancels identically in F_EMC = F_Newton.
   The gate passes for any value of G_geom: it checks the normalization
   of A, K_emc and 4pi, not G_geom.
3. Residual Floor: the reported ~1.2e-5 % difference is accumulated
   floating-point roundoff of the midpoint summation over a constant
   integrand, not a physical residual.
"""

import math
import sys

try:
    from m4_7_ewt_emergence_engine import (
        PI,
        C0,
        M_E,
        R_E,
        E_CHARGE_CODATA,
        G_CODATA,
        BCC_IDEAL_PROJECTION_LP,
        compute_alpha_geometric,
        derive_eps_M_from_BCC,
        derive_planck_charge_from_e,
        derive_neutrino_radius,
        derive_lambda_l_geometric,
        gravity_sector,
    )
except ImportError:
    raise ImportError("This module requires m4_7_ewt_emergence_engine.py in the same directory.")


# ----------------------------------------------------------------------
# 1. Geometric G from the self-consistent trinity
# ----------------------------------------------------------------------


def derive_G_geom():
    """Derive G_geom from the M4.7 self-consistent trinity.

    Returns
    -------
    dict with:
        G_geom        : derived gravitational constant [m^3 kg^-1 s^-2]
        lambda_l      : derived EMC lattice spacing [m]
        alpha_geom    : geometric fine-structure constant
        r_nu          : neutrino radius [m]
        N_geom        : effective geometric stiffness
        rel_err_vs_CODATA        : relative residual vs CODATA value
        rel_err_over_uncertainty : residual as a multiple of CODATA uncertainty
    """
    bcc = derive_eps_M_from_BCC(8.0 * PI**4)
    N_geom = bcc["N_geom"]
    eps_M = bcc["eps_M"]

    alpha_inv = compute_alpha_geometric(eps_M)
    alpha_geom = 1.0 / alpha_inv

    # Dimensional anchors r_e, m_e, c are named here at the point of use.
    q_P = derive_planck_charge_from_e(alpha_geom, E_CHARGE_CODATA)
    nu_res = derive_neutrino_radius(alpha_geom, q_P)
    r_nu = nu_res["r_nu"]

    lambda_l = derive_lambda_l_geometric(
        alpha_geom=alpha_geom,
        r_e=R_E,
        r_nu=r_nu,
        N_geom=N_geom,
        L_p_geom=BCC_IDEAL_PROJECTION_LP,
        K_WC=10,
    )

    res = gravity_sector(
        alpha_geom=alpha_geom,
        r_nu=r_nu,
        N_geom=N_geom,
        L_p_geom=BCC_IDEAL_PROJECTION_LP,
        K_WC=10,
        lambda_l=lambda_l,
        r_e=R_E,
        m_e=M_E,
        c0=C0,
    )
    G_geom = res["G_EWT"]

    rel_err = abs(G_geom - G_CODATA) / G_CODATA
    codata_uncertainty = 2.2e-5  # CODATA 2022 relative uncertainty on G
    rel_err_over_uncertainty = rel_err / codata_uncertainty

    return {
        "G_geom": G_geom,
        "lambda_l": lambda_l,
        "alpha_geom": alpha_geom,
        "r_nu": r_nu,
        "N_geom": N_geom,
        "rel_err_vs_CODATA": rel_err,
        "rel_err_over_uncertainty": rel_err_over_uncertainty,
    }


# ----------------------------------------------------------------------
# 2. Exact-domain overlap integral
# ----------------------------------------------------------------------


def compute_overlap_integral_exact_domain(A1, A2, R, num_pts=10000):
    """Evaluate the spatial overlap integral I(R) over r in [R, inf)
    using the coordinate transformation t in [0, 1).
    """

    def overlap_integral_at(dist):
        dt = 1.0 / num_pts
        i_sum = 0.0

        for i in range(num_pts):
            t = (i + 0.5) * dt
            r = dist / (1.0 - t)
            dr_dt = dist / ((1.0 - t) * (1.0 - t))

            angular_integral = 4.0 * math.pi / (r * r)
            i_sum += angular_integral * dr_dt * dt

        return A1 * A2 * i_sum

    dR = R * 1e-6
    i_plus = overlap_integral_at(R + dR)
    i_minus = overlap_integral_at(R - dR)

    f_numeric = -(i_plus - i_minus) / (2.0 * dR)
    return f_numeric


# ----------------------------------------------------------------------
# 3. Main
# ----------------------------------------------------------------------


def main():
    print("[1/5] Deriving G_geom from the M4.7 self-consistent trinity...")
    geom = derive_G_geom()
    G_geom = geom["G_geom"]

    print(f"    G_geom             = {G_geom:.15e} m^3 kg^-1 s^-2")
    print(f"    lambda_l (derived) = {geom['lambda_l']:.15e} m")
    print(f"    alpha_geom         = {geom['alpha_geom']:.15f}")
    print(f"    r_nu               = {geom['r_nu']:.15e} m")
    print(f"    N_geom             = {geom['N_geom']:.10f}")

    print("\n[2/5] Accuracy statement (strength clause, accuracy half)...")
    print(f"    G_CODATA            = {G_CODATA:.15e} m^3 kg^-1 s^-2")
    print(f"    relative residual   = {geom['rel_err_vs_CODATA']*100:.6f} %")
    print(f"    CODATA rel. unc. on G = 22 ppm (2.2e-5)")
    print(f"    residual / uncertainty = " f"{geom['rel_err_over_uncertainty']:.1f}x")
    print("    Note: r_e, m_e, c are the dimensional anchors. The derivation")
    print("          produces G_geom as a dimensionless ratio against")
    print("          c^2 r_e / m_e, not G from nothing.")

    print("\n[3/5] Loading physical parameters for the force test...")
    c = C0
    M1, M2 = 1.989e30, 5.972e24
    R = 1.495978707e11

    A1 = 2.0 * G_geom * M1 / (c * c)
    A2 = 2.0 * G_geom * M2 / (c * c)

    print(f"    c          = {c:.3f} m/s")
    print(f"    G_geom     = {G_geom:.15e} m^3 kg^-1 s^-2")
    print(f"    A1 (Sun)   = {A1:.3f} m")
    print(f"    A2 (Earth) = {A2:.3f} m")
    print(f"    R (1 AU)   = {R:.3e} m")

    print("\n[4/5] Performing mapped 3D spatial integration over [R, inf)...")
    F_numeric = compute_overlap_integral_exact_domain(A1, A2, R)

    print("\n[5/5] Converting field integral to physical force...")
    K_emc = c**4 / (16.0 * math.pi * G_geom)
    F_emc = K_emc * F_numeric

    F_newton = G_geom * M1 * M2 / (R * R)
    rel_diff = abs(F_emc - F_newton) / F_newton * 100.0

    print(f"    |F_EMC| (numerical) = {F_emc:.12e} N")
    print(f"    |F_Newton|          = {F_newton:.12e} N")
    print(f"    Rel. diff.          = {rel_diff:.12e}%")

    ok = rel_diff < 1e-3
    if ok:
        print("    RESULT: PASS (normalization of A, K_emc and 4pi consistent; G_geom cancels)")
    else:
        print("    RESULT: FAIL")
    return ok


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
