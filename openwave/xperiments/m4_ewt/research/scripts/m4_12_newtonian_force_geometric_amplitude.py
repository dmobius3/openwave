#!/usr/bin/env python3
"""M4/EWT - Newtonian Force with the Amplitude in M4.7 Chain Factors

OpenWave criterion:
    Gravity: Newton limit (GEM)

Purpose:
    Write the monopole amplitude A in the M4.7 chain factors, using the
    ratio N_nu_stat / N_nu_eff explicitly, without the symbol G_geom.
    The product is the engine's G_EWT formula times 2 M / c^2, term for
    term, so A = 2 G_geom M / c^2 by value and G_geom still cancels in
    the force gate. The gate remains a normalization-consistency check,
    as in M4.10.

    The mutation arm changes G_geom in the coupling only, one copy of
    the formula and not the other, so it confirms the identity and does
    not test the value of G_geom.

Dimensional anchors:
    r_e, m_e, c are the measured anchors used to build the dimensionless
    geometric ratio. The amplitude A is a length built from these anchors
    and the BCC geometry.

Structure:
    [1] Derive the M4.7 geometric trinity.
    [2] Write A in chain factors (no G symbol; equals 2 G_geom M / c^2).
    [3] State the accuracy half of the strength clause.
    [4] Numerical overlap integral and force baseline gate.
    [5] One-copy mutation: G_geom in the coupling only.
"""

import math
import sys

try:
    from m4_7_ewt_emergence_engine import (
        PI, C0, M_E, R_E, E_CHARGE_CODATA,
        G_CODATA,
        BCC_IDEAL_PROJECTION_LP,
        compute_alpha_core,
        compute_alpha_geometric,
        derive_eps_M_from_BCC,
        derive_planck_charge_from_e,
        derive_neutrino_radius,
        derive_lambda_l_geometric,
        gravity_sector,
    )
except ImportError as e:
    raise ImportError(
        "This module requires m4_7_ewt_emergence_engine.py in the same directory."
    ) from e


K_WC = 10


# ----------------------------------------------------------------------
# 1. Geometry from the M4.7 trinity
# ----------------------------------------------------------------------

def derive_geometry():
    """Return the M4.7 geometric quantities used by the amplitude test.

    Returns
    -------
    dict with G_geom, N_geom, X_eff, N_nu_stat, N_nu_eff, alpha_geom,
    lambda_l, r_nu, A_pi.
    """
    bcc = derive_eps_M_from_BCC(8.0 * PI**4)
    N_geom = bcc["N_geom"]
    eps_M = bcc["eps_M"]

    alpha_inv = compute_alpha_geometric(eps_M)
    alpha_geom = 1.0 / alpha_inv

    q_P = derive_planck_charge_from_e(alpha_geom, E_CHARGE_CODATA)
    nu_res = derive_neutrino_radius(alpha_geom, q_P)
    r_nu = nu_res["r_nu"]

    lambda_l = derive_lambda_l_geometric(
        alpha_geom=alpha_geom,
        r_e=R_E,
        r_nu=r_nu,
        N_geom=N_geom,
        L_p_geom=BCC_IDEAL_PROJECTION_LP,
        K_WC=K_WC,
    )

    res = gravity_sector(
        alpha_geom=alpha_geom,
        r_nu=r_nu,
        N_geom=N_geom,
        L_p_geom=BCC_IDEAL_PROJECTION_LP,
        K_WC=K_WC,
        lambda_l=lambda_l,
        r_e=R_E,
        m_e=M_E,
        c0=C0,
    )

    return {
        "G_geom": res["G_EWT"],
        "N_geom": N_geom,
        "N_nu_stat": res["N_nu_statutory"],
        "N_nu_eff": res["N_nu_eff"],
        "X_eff": res["X_eff"],
        "alpha_geom": alpha_geom,
        "lambda_l": lambda_l,
        "r_nu": r_nu,
        "A_pi": compute_alpha_core(),
    }


# ----------------------------------------------------------------------
# 2. Geometric amplitude from the EMC deficit
# ----------------------------------------------------------------------

def geometric_amplitude(M, geom):
    """Derive the monopole amplitude A from the EMC deficit ratio.

    A = 2 M r_e / m_e * sqrt(X_eff) / (A_pi^4 N_geom^3 K_WC sqrt(N_nu_stat))

    with X_eff = N_nu_stat / N_nu_eff. No G_geom symbol enters, but the
    product equals 2 G_geom M / c^2 (engine gravity_sector).
    """
    r_e = R_E
    m_e = M_E
    A_pi = geom["A_pi"]
    N_geom = geom["N_geom"]
    X_eff = geom["X_eff"]
    N_nu_stat = geom["N_nu_stat"]

    numerator = 2.0 * M * r_e * math.sqrt(X_eff)
    denominator = m_e * A_pi**4 * N_geom**3 * K_WC * math.sqrt(N_nu_stat)
    return numerator / denominator


# ----------------------------------------------------------------------
# 3. Exact-domain overlap integral
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
    return -(i_plus - i_minus) / (2.0 * dR)


# ----------------------------------------------------------------------
# 4. Force with a given G in the coupling only
# ----------------------------------------------------------------------

def force_ratio(A1, A2, M1, M2, R, G_for_coupling, F_numeric):
    """Return the EMC and Newton forces with G_for_coupling in K_emc and
    F_Newton, while A1 and A2 are kept as given (geometric)."""
    c = C0
    K_emc = c**4 / (16.0 * math.pi * G_for_coupling)
    F_emc = K_emc * F_numeric
    F_newton = G_for_coupling * M1 * M2 / (R * R)
    return F_emc, F_newton


# ----------------------------------------------------------------------
# 5. Main
# ----------------------------------------------------------------------

def main():
    print("[1/6] Deriving geometry from the M4.7 trinity...")
    geom = derive_geometry()
    G_geom = geom["G_geom"]

    print(f"    G_geom       = {G_geom:.15e} m^3 kg^-1 s^-2")
    print(f"    N_geom       = {geom['N_geom']:.10f}")
    print(f"    N_nu_stat    = {geom['N_nu_stat']:.6e}")
    print(f"    N_nu_eff     = {geom['N_nu_eff']:.6e}")
    print(f"    X_eff        = {geom['X_eff']:.6f}")
    print(f"    A_pi         = {geom['A_pi']:.10f}")

    print("\n[2/6] Deriving geometric amplitude from N_nu_stat / N_nu_eff...")
    M1, M2 = 1.989e30, 5.972e24
    R = 1.495978707e11
    c = C0

    A1 = geometric_amplitude(M1, geom)
    A2 = geometric_amplitude(M2, geom)

    # Cross-check against the 2 G_geom M / c^2 form (must agree)
    A1_alt = 2.0 * G_geom * M1 / (c * c)
    A2_alt = 2.0 * G_geom * M2 / (c * c)

    print(f"    A_1 (geometric) = {A1:.6e} m")
    print(f"    A_1 (2GM/c^2)   = {A1_alt:.6e} m")
    print(f"    rel. diff A_1   = {abs(A1 - A1_alt)/A1_alt:.3e}")
    print(f"    A_2 (geometric) = {A2:.6e} m")
    print(f"    A_2 (2GM/c^2)   = {A2_alt:.6e} m")
    print(f"    rel. diff A_2   = {abs(A2 - A2_alt)/A2_alt:.3e}")
    print("    The two forms are the same expression: the amplitude carries")
    print("    G_geom by value, so G_geom cancels in the force gate.")

    print("\n[3/6] Accuracy statement (strength clause, accuracy half)...")
    rel_err = abs(G_geom - G_CODATA) / G_CODATA
    codata_uncertainty = 2.2e-5
    rel_err_over_uncertainty = rel_err / codata_uncertainty
    print(f"    G_CODATA            = {G_CODATA:.15e} m^3 kg^-1 s^-2")
    print(f"    relative residual   = {rel_err*100:.6f} %")
    print(f"    residual / uncertainty = {rel_err_over_uncertainty:.1f}x")
    print("    Note: r_e, m_e, c are the dimensional anchors.")

    print("\n[4/6] Performing mapped 3D spatial integration over [R, inf)...")
    F_numeric = compute_overlap_integral_exact_domain(A1, A2, R)

    print("\n[5/6] Baseline gate with the derived G_geom in the coupling...")
    F_emc, F_newton = force_ratio(A1, A2, M1, M2, R, G_geom, F_numeric)
    rel_diff = abs(F_emc - F_newton) / F_newton * 100.0
    print(f"    |F_EMC|      = {F_emc:.12e} N")
    print(f"    |F_Newton|   = {F_newton:.12e} N")
    print(f"    Rel. diff.   = {rel_diff:.12e}%")
    baseline_pass = rel_diff < 1e-3
    if baseline_pass:
        print("    RESULT: PASS (normalization consistency; G_geom cancels)")
    else:
        print("    RESULT: FAIL")

    print("\n[6/6] Mutation test: replace G_geom with 1.0 in K_emc and F_Newton...")
    G_mut = 1.0
    F_emc_mut, F_newton_mut = force_ratio(A1, A2, M1, M2, R, G_mut, F_numeric)
    rel_diff_mut = abs(F_emc_mut - F_newton_mut) / F_newton_mut * 100.0
    print(f"    |F_EMC|_mut  = {F_emc_mut:.12e} N")
    print(f"    |F_Newton|_mut = {F_newton_mut:.12e} N")
    print(f"    Rel. diff._mut = {rel_diff_mut:.6f}%")
    mutation_fails = rel_diff_mut > 10.0
    if mutation_fails:
        print("    RESULT: FAIL as expected (one copy mutated; confirms the identity)")
    else:
        print("    RESULT: PASS (unexpected; the gate ignores the coupling)")

    overall = baseline_pass and mutation_fails
    print()
    if overall:
        print("OVERALL: PASS (normalization-consistency gate; G_geom cancels)")
    else:
        print("OVERALL: FAIL")
    return overall


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
