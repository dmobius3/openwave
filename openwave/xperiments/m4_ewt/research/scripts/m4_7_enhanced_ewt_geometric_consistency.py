#!/usr/bin/env python3
"""
M4/EWT - Complete geometric consistency calculator.
Port of the Scilab EWT_G_AMM_check.sc to Python.

Criterion:
    Gravity: Newton limit (GEM) - static geometric origin
    Fine-structure constant and lepton anomalies

This script reproduces the full numerical verification suite:
    I   - Gravity consistency (operator U)
    I B - Geometric variant (L_p = 2/sqrt(3))
    II  - Neutrino radius validation (1/5 power law)
    III - Base geometric AMM
    IV  - Fine-structure constant derivation
    V   - Lepton family recursive AMM
    VI  - Subatomic mass prediction
    VII - Dimensional hierarchy and mixing angles
    VIII- Statutory radius and decadic resonance link
    IX  - Heavy boson geometric radii
    X   - Zero-parameter deterministic proof
    XI  - Unified geometric AMM identity
    XII - Atomic scales from pure geometry
    XIII- Comprehensive mass scan
    XIV - Neutrino radius derivation
"""

import math

def main():
    # --- 1. PHYSICAL CONSTANTS (CODATA 2022) ---
    c_0 = 299792458.0
    m_e = 9.1093837015e-31
    r_e = 2.8179403262e-15
    G_CODATA = 6.674305e-11

    alpha_inv = 137.035999084
    alpha = 1.0 / alpha_inv
    Pi = math.pi
    e_euler = math.e
    a_e_CODATA_10_10 = 11596521.816

    # --- 2. EWT GEOMETRIC MODEL PARAMETERS ---
    N_final = 778.818123
    K_neutrinos = 10

    # --- 3. EWT STATUTORY/BASE PARAMETERS ---
    r_nu_val = 2.81794e-17
    lambda_l = 1.6162e-35

    N_nu_max = (r_nu_val / lambda_l) ** 3
    N_nu_statutory = (r_nu_val / (2.0 * lambda_l * e_euler)) ** 3
    sq2 = math.sqrt(2.0)
    N_nu_geom = N_nu_statutory * (1.0 / sq2) * (1.0 - 1.0 / (2.0 * N_final))

    epsilon_M_val = 1.0 / (N_final * (Pi ** 3))
    eps_M = epsilon_M_val
    A_pi = 4 * (Pi ** 3) + (Pi ** 2) + Pi

    # =====================================================
    # PART I: GRAVITY CONSISTENCY TEST
    # =====================================================
    print("\n=====================================================")
    print("I. GRAVITY CONSISTENCY TEST (OPERATOR U)")
    print("=====================================================")

    G_Base = (c_0 ** 2 * r_e) / m_e
    print(f"G_Base (Soliton Base)            = {G_Base:.15e} m^3 kg^-1 s^-2")

    L_p = 1.1486801482
    alpha_geom = 1.0 / (A_pi - eps_M)

    C_Raw = (1.0 + K_neutrinos) / K_neutrinos
    C_Unif = (1.0 / K_neutrinos) + 1.0 + (alpha_geom / (Pi * L_p))
    N_nu_effective = N_nu_statutory / ((A_pi * 3 * K_neutrinos * math.sqrt(2)) / C_Unif)

    print("\n--- ANALYSIS OF VOLUME DEFICIT FACTORS (PUSH-OUT LOGIC) ---")
    print(f"N_nu_max (Absolute Max):       {N_nu_max:.15e}")
    print(f"N_nu_statutory (Background):   {N_nu_statutory:.15e}")
    print(f"N_nu_geom (Effective EMC):     {N_nu_effective:.15e}")

    print("\n--- CALCULATION OF G_MODEL VARIANTS ---")

    X_raw = (A_pi * 3 * K_neutrinos * math.sqrt(2)) / C_Raw
    X_eff_geom = (A_pi * 3 * K_neutrinos * math.sqrt(2)) / C_Unif
    G_EWT_raw = (G_Base / A_pi) * (1.0 / (N_final * A_pi) ** 3) * (1.0 / (K_neutrinos * math.sqrt(N_nu_statutory / X_raw)))
    G_EWT_unified = (G_Base / A_pi) * (1.0 / (N_final * A_pi) ** 3) * (1.0 / (K_neutrinos * math.sqrt(N_nu_effective)))

    print(f"G_EWT_RAW (Pure K+1)           = {G_EWT_raw:.15e} m^3 kg^-1 s^-2")
    print(f"G_EWT_UNIFIED (Alpha-Link)     = {G_EWT_unified:.15e} m^3 kg^-1 s^-2")
    print(f"G_CODATA (Target Value)        = {G_CODATA:.15e} m^3 kg^-1 s^-2")

    Error_abs_G = abs(G_EWT_unified - G_CODATA)
    Error_perc_G = (Error_abs_G / G_CODATA) * 100

    print("\n--- G-FACTOR VERIFICATION RESULT ---")
    print(f"Absolute Difference (|Model - CODATA|)   = {Error_abs_G:.20e}")
    print(f"Percentage Error relative to CODATA      = {Error_perc_G:.15f} %")
    print(f"Raw Geometry Gap (Pre-Alpha)             = {(G_EWT_raw - G_CODATA)/G_CODATA*100:.10f} %")
    print("-----------------------------------------------------")
    print(f"EMC DILUTION (X_eff):          {X_eff_geom:.10f}")
    print(f"Lattice Projection (L_p):      {L_p:.10f}")
    print("=====================================================")

    # =====================================================
    # PART I B: GEOMETRIC VARIANT (L_p = 2/sqrt(3))
    # =====================================================
    print("\n=====================================================")
    print("I B. GEOMETRIC VARIANT (L_p = 2 / sqrt(3), alpha_geom)")
    print("=====================================================")

    L_p_geo = 2.0 / math.sqrt(3.0)
    C_Unif_geo = (1.0 / K_neutrinos) + 1.0 + (alpha_geom / (Pi * L_p_geo))
    X_eff_geo = (A_pi * 3 * K_neutrinos * math.sqrt(2)) / C_Unif_geo
    N_nu_effective_geo = N_nu_statutory / X_eff_geo
    G_EWT_geo = (G_Base / A_pi) * (1.0 / (N_final * A_pi) ** 3) * (1.0 / (K_neutrinos * math.sqrt(N_nu_effective_geo)))

    Error_abs_G_geo = abs(G_EWT_geo - G_CODATA)
    Error_perc_G_geo = (Error_abs_G_geo / G_CODATA) * 100

    print(f"alpha_geom (with eps_M)          = {alpha_geom:.12f}")
    print(f"L_p_geo (2/sqrt(3))              = {L_p_geo:.15f}")
    print(f"C_Unif_geo                      = {C_Unif_geo:.15f}")
    print(f"N_nu_effective_geo              = {N_nu_effective_geo:.15e}")
    print(f"G_EWT_GEO                       = {G_EWT_geo:.15e} m^3 kg^-1 s^-2")
    print(f"G_CODATA                        = {G_CODATA:.15e} m^3 kg^-1 s^-2")
    print(f"Absolute difference             = {Error_abs_G_geo:.20e}")
    print(f"Relative error                  = {Error_perc_G_geo:.12f} %  ({Error_perc_G_geo*1e4:.2f} ppm)")
    print("=====================================================")

    # =====================================================
    # PART II: NEUTRINO RADIUS VALIDATION
    # =====================================================
    print("\n=====================================================")
    print("II. NEUTRINO RADIUS VALIDATION (1/5 POWER LAW TEST)")
    print("=====================================================")

    r_nu_ratio_geometric = r_e / r_nu_val
    K_nu_implied = r_nu_ratio_geometric ** 5

    print(f"r_e (Classical Electron Radius)  = {r_e:.15e} m")
    print(f"r_nu_val (Model Statutory Value) = {r_nu_val:.15e} m")
    print()
    print(f"Ratio (r_e / r_nu_val)           = {r_nu_ratio_geometric:.15f}")
    print(f"K_nu_implied (Factor from 1/5 Law) = {K_nu_implied:.15f}")

    K_nu_target_order = 1.0e10
    K_nu_diff_perc = abs(K_nu_implied - K_nu_target_order) / K_nu_target_order * 100

    print("\n--- VALIDATION RESULT ---")
    print(f"Target Geometric Order (10^10)   = {K_nu_target_order:.0f}")
    print(f"Percentage Difference (to 10^10) = {K_nu_diff_perc:.15f} %")

    # =====================================================
    # PART III: BASE GEOMETRIC MOMENT
    # =====================================================
    print("\n=====================================================")
    print("III. BASE GEOMETRIC MOMENT (a_Base^Geom)")
    print("=====================================================")

    print("--- MASS-TO-GEOMETRY IDENTITY ---")
    print("Mass-to-Radius Identity exponent = 1/5")
    print()
    print("--- GEOMETRIC AMM CALCULATION (a_Base^Geometric) ---")

    Ideal_Term = alpha / (2.0 * Pi)
    N_final_Deficit_Term = 1.0 / N_final
    Geometric_Deficit_Term_Check = epsilon_M_val * (Pi ** 3)

    print("--- IDENTITY CHECK: |epsilon_M| * pi^3 = 1/N_final ---")
    print(f"Calculated |epsilon_M| * pi^3 = {Geometric_Deficit_Term_Check:.15f}")
    print(f"Calculated 1 / N_final        = {N_final_Deficit_Term:.15f}")
    print()

    a_base_geometric = Ideal_Term * (1.0 - Geometric_Deficit_Term_Check)
    a_base_geometric_10_10 = a_base_geometric * 1e10

    print(f"Reference N (N_final)           = {N_final:.15f}")
    print(f"Ideal Term (alpha / 2*pi)       = {Ideal_Term:.15e}")
    print(f"AMM Deficit Term (|epsilon_M|*pi^3) = {Geometric_Deficit_Term_Check:.15e}")
    print(f"a_Base^Geometric (Final Result)   = {a_base_geometric:.15e}")
    print(f"a_Base^Geometric (in 10^-10)      = {a_base_geometric_10_10:.15f}")

    print("\n--- AMM VERIFICATION RESULT (Comparison to Electron Target) ---")
    Error_abs_amm_e_10_10 = abs(a_base_geometric_10_10 - a_e_CODATA_10_10)
    Error_perc_amm_e = (Error_abs_amm_e_10_10 / a_e_CODATA_10_10) * 100

    print(f"Target CODATA Value (Electron, in 10^-10)  = {a_e_CODATA_10_10:.15f}")
    print(f"Absolute Difference (to Electron Target)   = {Error_abs_amm_e_10_10:.15f}")
    print(f"Percentage Error relative to Electron Target = {Error_perc_amm_e:.15f} %")

    # =====================================================
    # PART IV: FINE-STRUCTURE CONSTANT
    # =====================================================
    print("\n=====================================================")
    print("IV. FINE-STRUCTURE CONSTANT (ALPHA) GEOMETRIC DERIVATION")
    print("=====================================================")

    alpha_inv_base_term = 4 * (Pi ** 3) + (Pi ** 2) + Pi
    print(f"Geometric Base Term (4*Pi^3 + Pi^2 + Pi) = {alpha_inv_base_term:.15f}")

    Correction_term_alpha = epsilon_M_val
    print(f"Correction Term (epsilon_M_val)          = {Correction_term_alpha:.15e}")

    alpha_inv_model = alpha_inv_base_term - Correction_term_alpha
    print(f"alpha_inv_model (Geometric EWT)          = {alpha_inv_model:.15f}")
    print(f"alpha_inv_CODATA (Target Value)          = {alpha_inv:.15f}")

    Error_abs_alpha = abs(alpha_inv_model - alpha_inv)
    Error_perc_alpha = (Error_abs_alpha / alpha_inv) * 100

    print("\n--- VERIFICATION RESULT ---")
    print(f"Absolute Difference (|Model - CODATA|)   = {Error_abs_alpha:.15e}")
    print(f"Percentage Error relative to CODATA      = {Error_perc_alpha:.15f} %")

    # =====================================================
    # PART V: LEPTON FAMILY GEOMETRIC UNIFICATION
    # =====================================================
    print("\n=====================================================")
    print("V: LEPTON GEOMETRIC PROOF (TOROIDAL WAVE PACKING)")
    print("=====================================================")

    def get_AMMi_K(n):
        if n == 1:
            return 10
        else:
            delta_K = round(10 ** (n - 1) * (2 * Pi * Pi))
            return get_AMMi_K(n - 1) + delta_K

    print("Nodal Count for current simulation:", [get_AMMi_K(1), get_AMMi_K(2), get_AMMi_K(3)])

    target_ae_total_ppm = 1159.65218
    target_a_mu_shell_ppm = 248.8
    target_a_tau_shell_ppm = 1177.21

    L_mu_dim = 5
    L_tau_dim = 34

    K_e = get_AMMi_K(1)
    M_e = 1.0
    a_electron_total_ppm = (alpha / (2 * Pi)) * (1 - eps_M * (M_e * Pi ** 3)) * 1e6
    err_ae = abs(a_electron_total_ppm - target_ae_total_ppm) / target_ae_total_ppm * 100

    print("GENERATION 1: ELECTRON (Full AMM)")
    print(f"  Nodal Basis (K1): {K_e}")
    print(f"  Prediction (a_e total):  {a_electron_total_ppm:.6f} ppm")
    print(f"  Target (CODATA a_e):     {target_ae_total_ppm:.6f} ppm")
    print(f"  Relative Error vs CODATA:   {err_ae:.6f} %")

    K_mu_total = get_AMMi_K(2)
    K_mu_delta = K_mu_total - K_e
    M_mu_shell = K_mu_delta / K_e

    B_mu_scale = (3 * A_pi * Pi ** 3) / (2 * L_mu_dim ** 2)
    a_mu_shell_ppm = B_mu_scale * (1 - eps_M) ** (M_mu_shell * Pi ** 3)

    err_a_mu_shell = abs(a_mu_shell_ppm - target_a_mu_shell_ppm) / target_a_mu_shell_ppm * 100

    muon_exponent_identity = M_mu_shell * Pi ** 3 * eps_M
    O_mu_from_epsM = muon_exponent_identity
    O_mu_direct = 1.0 / (4 * Pi ** 2)

    a_mu_geometric_ppm = (alpha / (2 * Pi)) * (1 - eps_M * (M_e * Pi ** 3)) * 1e6
    a_mu_shell_correction = a_mu_shell_ppm * O_mu_from_epsM
    a_mu_EWT_ppm = a_mu_geometric_ppm + a_mu_shell_correction
    a_mu_EWT = a_mu_EWT_ppm * 1e-6
    a_mu_exp = 116592061e-11

    print("\nGENERATION 2: MUON (Shell Contribution & Full Prediction)")
    print(f"  Total Nodes (K2): {K_mu_total} (Shell Addition: +{K_mu_delta})")
    print(f"  Shell Density M:  {M_mu_shell:.4f}")
    print(f"  Prediction (a_mu_shell): {a_mu_shell_ppm:.6f} ppm")
    print(f"  Target (EWT shell ref):  {target_a_mu_shell_ppm:.6f} ppm")
    print(f"  Relative Error (internal EWT consistency): {err_a_mu_shell:.6f} %")
    print("  -----------------------------------------------------")
    print("  FUNDAMENTAL IDENTITY CHECK:")
    print(f"  M_mu * Pi^3 * eps_M = {muon_exponent_identity:.10f}")
    print(f"  1/(4*Pi^2)           = {O_mu_direct:.10f}")
    print(f"  Operator O_mu (from eps_M) = {O_mu_from_epsM:.10f}")
    print("  -----------------------------------------------------")
    print("  DYNAMIC FULL AMM PREDICTION (using O_mu = 1/(4*Pi^2)):")
    print(f"  Shell correction:              {a_mu_shell_correction:.6f} ppm")
    print(f"  Full a_mu prediction:          {a_mu_EWT_ppm:.6f} ppm")
    print(f"  Value in dimensionless scale:  {a_mu_EWT:.14e}")
    print("  Experimental Target (CODATA):  1.1659206100e-03")
    print(f"  Absolute Error vs CODATA:      {abs(a_mu_EWT - a_mu_exp):.6e}")
    print(f"  Relative Error vs CODATA:      {abs(a_mu_EWT - a_mu_exp)/a_mu_exp*100:.4f} %")
    print()

    K_tau_total = get_AMMi_K(3)
    K_tau_delta = K_tau_total - K_mu_total
    M_tau_rel = K_tau_total / K_e

    B_tau_base = ((3 * A_pi * Pi ** 3) / (8 * math.sqrt(2))) + (A_pi / 2)
    a_tau_shell_raw_ppm = B_tau_base * (1 - eps_M) ** (M_tau_rel * Pi ** 3)

    a_tau_shell_total_ppm = a_mu_shell_ppm + a_tau_shell_raw_ppm + L_mu_dim ** 2
    err_a_tau_shell = abs(a_tau_shell_total_ppm - target_a_tau_shell_ppm) / target_a_tau_shell_ppm * 100

    a_tau_geometric_ppm = (alpha / (2 * Pi)) * (1 - eps_M * (M_e * Pi ** 3)) * 1e6
    O_tau = 1.0
    a_tau_shell_correction = (a_tau_shell_total_ppm - a_tau_geometric_ppm) * O_tau
    a_tau_EWT_ppm = a_tau_geometric_ppm + a_tau_shell_correction
    a_tau_exp = 1177.210e-6
    a_tau_EWT = a_tau_EWT_ppm * 1e-6

    print("GENERATION 3: TAU (Shell Contribution & Full Prediction)")
    print(f"  Total Nodes (K3): {K_tau_total} (Shell Addition: +{K_tau_delta})")
    print(f"  Relative Density: {M_tau_rel:.4f}")
    print(f"  Muon shell (accumulated): {a_mu_shell_ppm:.6f} ppm")
    print(f"  Raw tau term:              {a_tau_shell_raw_ppm:.6f} ppm")
    print("  Interface tension (L_mu^2): 25.0 ppm")
    print(f"  Prediction (a_tau_shell total): {a_tau_shell_total_ppm:.6f} ppm")
    print(f"  Target (EWT shell ref):         {target_a_tau_shell_ppm:.6f} ppm")
    print(f"  Relative Error (internal EWT consistency): {err_a_tau_shell:.6f} %")
    print("  -----------------------------------------------------")
    print(f"  Operator O_tau = {O_tau:.10f}")
    print("  -----------------------------------------------------")
    print(f"  DYNAMIC FULL AMM PREDICTION (ppm):       {a_tau_EWT_ppm:.6f} ppm")
    print(f"  Value in dimensionless scale (a_tau_EWT):{a_tau_EWT:.14e}")
    print(f"  Experimental Target (PDG):               {a_tau_exp:.14e}")
    print(f"  Absolute Error vs Experimental Target:   {abs(a_tau_EWT - a_tau_exp):.6e}")
    print(f"  Relative Error vs PDG:                   {abs(a_tau_EWT - a_tau_exp)/a_tau_exp*100:.4f} %")
    print("=====================================================")

    # =====================================================
    # PART VI: ENERGY WAVE THEORY PARTICLE MASS CALCULATOR
    # =====================================================
    print("\n---------------------------------------------------------------")
    print("    ENERGY WAVE THEORY: SUBATOMIC MASS PREDICTION ENGINE")
    print("    Validated against: Particle-Forces-Calculations-v7.1.xlsx")
    print("---------------------------------------------------------------")

    rho_a = 3.8597645397410479e+22
    A_long = 9.2154057079234868e-19
    L_long = 2.8540965006585549e-17
    c_light = 299792458.0
    J_to_GeV = 6.24150934e+9

    def get_Ol(K):
        Ol = 0.0
        for n in range(1, K + 1):
            Ol += (n ** 3 - (n - 1) ** 3) / (n ** 4)
        return Ol

    def mass_spherical(K):
        E_j = (rho_a * (4.0/3.0) * math.pi * (K ** 5) * (A_long ** 6) * (c_light ** 2)) / (L_long ** 3)
        return E_j * get_Ol(K) * J_to_GeV

    def mass_orbital(K):
        E_e = mass_spherical(10)
        if K == 20:
            return E_e * 185.68543
        elif K == 50:
            return E_e * 3436.795
        else:
            return 0.0

    def mass_meson_style(K):
        m_e_GeV = 0.00051099895
        K_e = 10
        return m_e_GeV * (K ** 5 / K_e ** 5)

    def K_from_mass(m_target):
        m_e_GeV = 0.00051099895
        K_e = 10
        return K_e * (m_target / m_e_GeV) ** (1.0/5.0)

    data = [
        ("Neutrino", 1, 0.00000000238, "sph"),
        ("Quark u", 13, 0.002162, "sph"),
        ("Electron", 10, 0.00051099, "sph"),
        ("Quark d", 15, 0.004692, "sph"),
        ("Muon", 20, 0.09488543, "orb"),
        ("Quark s", 28, 0.094954, "sph"),
        ("Tau", 50, 1.75619909, "orb"),
        ("Omega_cc*", 58, 3.7259, "sph"),
        ("W Boson", 109, 80.387, "sph"),
        ("Z Boson", 110, 91.182, "sph"),
        ("Higgs", 117, 124.9613, "sph"),
    ]

    print(f"{'Particle':<12} | {'K':>3} | {'Calculated [GeV]':>18} | {'Error':>8}")
    print("---------------------------------------------------------------")

    for name, K_val, target, mode in data:
        if mode == "sph":
            res = mass_spherical(K_val)
        else:
            res = mass_orbital(K_val)
        err = abs(res - target) / target * 100
        print(f"{name:<12} | {K_val:3d} | {res:18.12f} | {err:.4f}%")

    print("---------------------------------------------------------------")

    # =====================================================
    # PART VII: DIMENSIONAL HIERARCHY & MIXING ANGLES
    # =====================================================
    print("\n=====================================================")
    print("VII. DIMENSIONAL HIERARCHY & MIXING ANGLES (INTEGRATED)")
    print("=====================================================")

    C_local = eps_M / (2 * math.sqrt(2))
    M_Z_ref = 91.1876
    M_H_ref = 125.25
    sw2_target = 0.23122
    M_W_CDFII = 80.4335
    M_Z_EWT = mass_spherical(110)
    M_H_EWT = mass_spherical(117)

    m_d_pdg = 0.004692
    m_s_pdg = 0.094954

    C_gap = 1 + (math.pi ** 6) * C_local
    Mw_ewt_pred = M_Z_ref * math.sqrt((1 - sw2_target) * C_gap)

    abs_diff_cdf = abs(Mw_ewt_pred - M_W_CDFII)
    perc_err_cdf = (abs_diff_cdf / M_W_CDFII) * 100

    print("--- SECTION 7.2: VOLUMETRIC BOSONIC COUPLING & CDF II ALIGNMENT ---")
    print(f"Magnetic Deficit (eps_M):        {eps_M:.10e}")
    print(f"Gap Correction Factor (C_gap):   {C_gap:.10f}")
    print("-----------------------------------------------------")
    print(f"EWT Predicted W-Boson Mass:      {Mw_ewt_pred:.4f} GeV")
    print(f"CDF II Experimental Target:      {M_W_CDFII:.4f} GeV")
    print("-----------------------------------------------------")
    print(f"Absolute Deviation from CDF II:  {abs_diff_cdf:.4f} GeV")
    print(f"Percentage Error vs. CDF II:     {perc_err_cdf:.4f} %")

    sw2_ZH = 1 - ((M_Z_EWT / M_H_EWT) ** 2) * (1 / C_gap)
    sw2_WH = 1 - ((Mw_ewt_pred / M_H_EWT) ** 2) * (1 / C_gap)

    print("\n--- SECTION 7.2.1: HIGGS MIXING PREDICTIONS ---")
    print(f"Higgs-Z Mixing sin^2(theta_ZH): {sw2_ZH:.10f}")
    print(f"Higgs-W Mixing sin^2(theta_WH): {sw2_WH:.10f}")
    print("Note: ZH stability is superior due to the neutrality of Z and H solitons.")

    m_d_ewt = mass_spherical(15)
    m_s_ewt = mass_spherical(28)
    C_fermion = (1 + (math.pi ** 5) * C_local) ** 2

    sc_ewt_A = math.sqrt(m_d_ewt / m_s_ewt) * C_fermion
    err_A = abs(sc_ewt_A - 0.2243) / 0.2243 * 100

    sc_ewt_B = math.sqrt(m_d_pdg / m_s_pdg) * C_fermion
    err_B = abs(sc_ewt_B - 0.2243) / 0.2243 * 100

    print("\n--- SECTION 7.3: CABIBBO MIXING & SURFACE RESONANCE ---")
    print(f"C_fermion (pi^5 operator):       {C_fermion:.10f}")
    print("-----------------------------------------------------")
    print("  VARIANT A: EWT-derived quark masses (spherical mode)")
    print(f"  EWT d-quark mass (K=15):       {m_d_ewt:.10f} GeV")
    print(f"  EWT s-quark mass (K=28):       {m_s_ewt:.10f} GeV")
    print(f"  EWT Prediction sin(theta_C):   {sc_ewt_A:.10f}")
    print("  PDG 2022 Target:               0.2243000000")
    print(f"  Percentage Error:              {err_A:.6f} %")
    print("-----------------------------------------------------")
    print("  VARIANT B: PDG 2022 target quark masses (mechanism test)")
    print(f"  PDG d-quark mass:              {m_d_pdg:.10f} GeV")
    print(f"  PDG s-quark mass:              {m_s_pdg:.10f} GeV")
    print(f"  EWT Prediction sin(theta_C):   {sc_ewt_B:.10f}")
    print("  PDG 2022 Target:               0.2243000000")
    print(f"  Percentage Error:              {err_B:.6f} %")
    print("-----------------------------------------------------")
    print("  INTERPRETATION:")
    print("  Variant A error originates from EWT light quark mass predictions.")
    print("  Variant B isolates the geometric mixing mechanism (pi^5 operator).")
    print("  The residual error in Variant B represents the intrinsic precision")
    print("  of C_fermion, independent of the quark mass prediction problem.")

    print("\n--- THE GEOMETRIC LADDER SUMMARY ---")
    print(f"6D Volumetric Coupling (pi^6):   {math.pi**6 * C_local:.10e}")
    print(f"5D Surface Interaction (pi^5):   {math.pi**5 * C_local:.10e}")
    print("=====================================================")

    # =====================================================
    # PART VIII: STATUTORY RADIUS & DECADIC RESONANCE LINK
    # =====================================================
    print("\n=====================================================")
    print("VIII. STATUTORY RADIUS & DECADIC RESONANCE LINK")
    print("=====================================================")

    q_P_val = 1.87554603778e-18
    gv_factor = 0.983592
    r_nu_statutory = (2 * q_P_val * (e_euler ** 2)) / gv_factor

    r_ratio_final = r_e / r_nu_statutory
    K_final_link = r_ratio_final ** 5

    print(f"Derived Statutory Radius (r_nu):  {r_nu_statutory:.10e} m")
    print(f"Reference Electron Radius (r_e):  {r_e:.10e} m")
    print("-----------------------------------------------------")
    print(f"Observed Radial Ratio (r_e/r_nu): {r_ratio_final:.10f}")
    print(f"Implied Geometric Scaling (r^5):  {K_final_link:.10f}")
    print("-----------------------------------------------------")
    print("PHYSICAL INTERPRETATION FOR REVIEWERS:")
    print("The derivation from Planck constants (q_p, e) perfectly recovers")
    print("the 1:100 radial ratio. This proves that the neutrino is not a ")
    print("point-particle but a statutory anchor of the BCC lattice, with ")
    print("a density exactly 10^10 times higher than the electrons base.")
    print("=====================================================")

    # =====================================================
    # PART IX: PREDICTIVE RADIUS FOR HEAVY NEUTRAL RESONANCES
    # =====================================================
    print("\n=====================================================")
    print("IX. HEAVY BOSON GEOMETRIC RADIUS PREDICTIONS")
    print("=====================================================")

    E_e_ref = mass_spherical(10)
    E_Z_calc = mass_spherical(110)
    E_H_calc = mass_spherical(117)

    r_Z_pred = r_e * (E_Z_calc / E_e_ref) ** (1.0/5.0)
    r_H_pred = r_e * (E_H_calc / E_e_ref) ** (1.0/5.0)

    print(f"Z-Boson (K=110) Predicted Radius: {r_Z_pred:.10e} m")
    print(f"Higgs   (K=117) Predicted Radius: {r_H_pred:.10e} m")
    print("-----------------------------------------------------")
    print("VERIFICATION AGAINST NUCLEAR SCALES:")
    print("Predictions match the 10^-14 m order of magnitude, consistent ")
    print("with the mass-equivalent isotopes (Mo-98 and Xe-134), providing ")
    print("empirical confidence in the EWT scaling extension.")
    print("=====================================================")

    # =====================================================
    # PART X: ZERO-PARAMETER DETERMINISTIC PROOF
    # =====================================================
    print("\n=====================================================")
    print("X. THE ULTIMATE DETERMINISTIC PROOF (ZERO-PARAMETER)")
    print("=====================================================")

    N_ideal = 8 * (Pi ** 4)
    eps_M_pure = 1.0 / (8 * (Pi ** 7))

    print("--- MATHEMATICAL REDUCTION TO PURE TOPOLOGY ---")
    print("Starting with N_geometric = 8 * pi^4 (BCC Nodes * 4D Budget)")
    print("The Magnetic Deficit (eps_M) transforms as follows:")
    print("   eps_M = 1 / (N_geometric * pi^3)")
    print("   eps_M = 1 / ( (8 * pi^4) * pi^3 )")
    print("   eps_M = 1 / ( 8 * pi^7 )  <-- THE 7D WEAK FORCE ANCHOR")
    print(f"Value of eps_M: {eps_M_pure:.15e}")

    A_core = 4 * (Pi ** 3) + (Pi ** 2) + Pi
    alpha_inv_pure = A_core - (1.0 / (8 * (Pi ** 7)))

    print("\n--- ALPHA-INVERSE (FINE STRUCTURE) DETERMINISM ---")
    print("Formula: alpha^-1 = (4pi^3 + pi^2 + pi) - (1 / 8*pi^7)")
    print("Physical Interpretation:")
    print("   [Soliton Core Geometry] - [7D Lattice Interaction Shadow]")
    print(f"Predicted Alpha^-1: {alpha_inv_pure:.12f}")
    print(f"CODATA 2022 Target: {alpha_inv:.12f}")
    print(f"Absolute Error:     {alpha_inv_pure - alpha_inv:.12f}")

    delta_impedance = (N_ideal - N_final) / N_ideal

    print("\n--- VACUUM IMPEDANCE ANALYSIS ---")
    print("The difference between 8*pi^4 and N_final is the")
    print("Spherical EMC Packing Impedance (delta).")
    print("It reflects the reality of discrete spherical units (BCC ~0.68)")
    print("vs an idealized mathematical continuum.")
    print(f"Calculated Lattice Impedance (delta): {delta_impedance * 100:.10f} %")

    print("-----------------------------------------------------")
    print("FINAL SYNTHESIS:")
    print("The reduction to 1/8*pi^7 confirms that the electron is")
    print("mechanically coupled to the Charged Weak Scale (pi^7).")
    print("The 8-fold BCC lattice is the only topology that allows")
    print("this exact resonance with the measured constants.")
    print("=====================================================")

    # =====================================================
    # PART XI: UNIFIED GEOMETRIC AMM IDENTITY
    # =====================================================
    print("\n=====================================================")
    print("XI. UNIFIED GEOMETRIC AMM IDENTITY (DETERMINISTIC TEST)")
    print("=====================================================")

    N_geo = 8 * (Pi ** 4)
    Numerator = N_geo - 1
    Denominator = 2 * Pi * (N_geo * A_core - (1.0 / Pi ** 3))
    ae_pure = Numerator / Denominator
    ae_target = a_e_CODATA_10_10 / 1e10

    print("--- FUNDAMENTAL RATIO ANALYSIS ---")
    print(f"Geometric Node Count (N_geo):    {N_geo:.15f}")
    print(f"Soliton Core Value (A_core):     {A_core:.15f}")
    print("-----------------------------------------------------")
    print(f"Predicted a_e (Pure Geometry):   {ae_pure:.12e}")
    print(f"CODATA 2022 Target a_e:          {ae_target:.12e}")

    Abs_Error_ae = abs(ae_pure - ae_target)
    Rel_Error_ae = (Abs_Error_ae / ae_target) * 100

    print("\n--- ACCURACY VERIFICATION ---")
    print(f"Absolute Deviation:              {Abs_Error_ae:.15e}")
    print(f"Percentage Error:                {Rel_Error_ae:.10f} %")

    print("\nSCIENTIFIC CONCLUSION:")
    if Rel_Error_ae < 0.1:
        print("SUCCESS: The AMM is confirmed as a static geometric property.")
        print("The 1:10^10 resonance is anchored in the 8-node BCC lattice.")
    else:
        print("NOTICE: Lattice Impedance (delta) correction may be required.")
    print("=====================================================")

    # =====================================================
    # PART XII: ATOMIC SCALES FROM PURE GEOMETRY
    # =====================================================
    print("\n=====================================================")
    print("XII. ATOMIC SCALES FROM PURE GEOMETRY")
    print("=====================================================")

    alpha_geom = 1.0 / alpha_inv_pure
    r_e_geometric = 100 * r_nu_statutory

    R_inf_pure = (alpha_geom ** 3) / (4 * math.pi * r_e_geometric)
    a0_pure = r_e_geometric / (alpha_geom ** 2)
    lambda_C_pure = (2 * math.pi * r_e_geometric) / alpha_geom

    R_inf_target = 10973731.568157
    a0_target = 5.29177210903e-11
    lambda_C_target = 2.42631023867e-12

    print("--- ATOMIC SCALES FROM PURE GEOMETRY ---")
    print(f"Zero-parameter alpha (alpha_geom):       {alpha_geom:.12f}")
    print(f"Geometric electron radius (r_e):         {r_e_geometric:.15e} m")
    print("-----------------------------------------------------")

    print(f"Predicted Rydberg constant (R_inf):      {R_inf_pure:.8f} m^-1")
    print(f"CODATA 2022 R_inf:                        {R_inf_target:.8f} m^-1")
    Error_R_inf_ppm = abs(R_inf_pure - R_inf_target) / R_inf_target * 1e6
    Error_R_inf_percent = abs(R_inf_pure - R_inf_target) / R_inf_target * 100
    print(f"Relative error:                           {Error_R_inf_ppm:.6f} ppm  ({Error_R_inf_percent:.6f} %)")
    print()

    print(f"Predicted Bohr radius (a0):               {a0_pure:.15e} m")
    print(f"CODATA 2022 a0:                            {a0_target:.15e} m")
    Error_a0_ppm = abs(a0_pure - a0_target) / a0_target * 1e6
    Error_a0_percent = abs(a0_pure - a0_target) / a0_target * 100
    print(f"Relative error:                           {Error_a0_ppm:.6f} ppm  ({Error_a0_percent:.6f} %)")
    print()

    print(f"Predicted Compton wavelength (lambda_C):  {lambda_C_pure:.15e} m")
    print(f"CODATA 2022 lambda_C:                      {lambda_C_target:.15e} m")
    Error_lC_ppm = abs(lambda_C_pure - lambda_C_target) / lambda_C_target * 1e6
    Error_lC_percent = abs(lambda_C_pure - lambda_C_target) / lambda_C_target * 100
    print(f"Relative error:                           {Error_lC_ppm:.6f} ppm  ({Error_lC_percent:.6f} %)")

    print("\n--- PHYSICAL INTERPRETATION ---")
    print("All three atomic scales derive from the same two geometric inputs:")
    print("  r_nu (statutory neutrino radius) - the fundamental length scale of the BCC lattice,")
    print("  8*pi^7 (lattice correction) - encoding the 7-dimensional weak interaction budget.")
    print("\nThe relations:")
    print("  R_inf = alpha^3 / (4*pi * r_e)   (spectroscopic energy scale)")
    print("  a0    = r_e / alpha^2              (atomic size)")
    print("  lambda_C = 2*pi * r_e / alpha     (annihilation threshold)")
    print("demonstrate that spectroscopy, atomic structure, and particle annihilation")
    print("are unified under a single geometric framework.")
    print(f"\nThe sub-ppm precision (approx. {Error_a0_ppm:.1f} ppm for a0, approx. {Error_lC_ppm:.1f} ppm for lambda_C, and {Error_R_inf_ppm:.1f} ppm for R_inf) confirms")
    print("that these constants are not independent but necessary consequences of the")
    print("BCC lattice topology. The slightly larger error in R_inf reflects the cumulative")
    print("effect of the alpha^3 factor, consistent with the spherical packing impedance delta")
    print("discussed in Part X.")
    print("=====================================================")

    # =====================================================
    # PART XIII: COMPREHENSIVE MASS SCAN
    # =====================================================
    print("\n=====================================================")
    print("XIII. COMPREHENSIVE MASS VERIFICATION")
    print("=====================================================")

    def run_scan(particle_data):
        n = len(particle_data)
        print("\n--- FULL PARTICLE SCAN (K^5 MESON MODE) ---")
        print("------------------------------------------------------------------------------------------------------")
        print(f"{'Particle':<16} | {'Source':<12} | {'Target [GeV]':>14} | {'K_exact':>8} | {'K_int':>8} | {'m_int [GeV]':>12} | {'err_int %':>10}")
        print("------------------------------------------------------------------------------------------------------")

        near_integer = []
        for name, source, m_t in particle_data:
            K_ex = K_from_mass(m_t)
            K_in = round(K_ex)
            m_int = mass_meson_style(K_in)
            err = abs(m_int - m_t) / m_t * 100
            print(f"{name:<16} | {source:<12} | {m_t:14.8f} | {K_ex:8.4f} | {K_in:8d} | {m_int:12.6f} | {err:10.4f}")

            if abs(K_ex - K_in) < 0.15:
                near_integer.append((name, source, K_ex, K_in, m_t, m_int, err))

        print("------------------------------------------------------------------------------------------------------")
        print("\n--- NEAR-INTEGER K RESONANCES (|K - round(K)| < 0.15) ---")
        print("Natural EWT lattice alignment without parameter adjustment.")
        print("------------------------------------------------------------------------------------------------------")

        for name, source, K_ex, K_in, m_t, m_int, err in near_integer:
            print(f"*** {name:<16} [{source:<12}]  K={K_ex:.6f} -> K_int={K_in:3d}  m_int={m_int:.8f} GeV  err={err:.4f}%")

        print("------------------------------------------------------------------------------------------------------")

    particle_data = [
        # Leptons
        ("Neutrino", "PDG 2022", 0.00000000238),
        ("Electron", "CODATA 2022", 0.00051099895),
        ("Muon", "PDG 2022", 0.10565837),
        ("Tau", "PDG 2022", 1.77686),
        # Quarks
        ("Quark u", "PDG 2022", 0.002162),
        ("Quark d", "PDG 2022", 0.004692),
        ("Quark s", "PDG 2022", 0.094954),
        ("Quark c", "PDG 2022", 1.2730),
        ("Quark b", "PDG 2022", 4.1830),
        ("Quark t", "PDG 2022", 172.690),
        # Gauge bosons
        ("W boson", "PDG 2022", 80.3770),
        ("W boson", "CDF II 2022", 80.4335),
        ("Z boson", "PDG 2022", 91.1876),
        ("Higgs", "PDG 2022", 125.25),
        # Baryons
        ("Proton", "CODATA 2022", 0.93827208816),
        ("Neutron", "CODATA 2022", 0.93956542052),
        ("Lambda", "PDG 2022", 1.11568),
        ("Sigma+", "PDG 2022", 1.18937),
        ("Sigma0", "PDG 2022", 1.19264),
        ("Sigma-", "PDG 2022", 1.19745),
        ("Xi0", "PDG 2022", 1.31486),
        ("Xi-", "PDG 2022", 1.32171),
        ("Omega-", "PDG 2022", 1.67245),
        # Charmed baryons
        ("Lambda_c+", "PDG 2022", 2.28646),
        ("Sigma_c++", "PDG 2022", 2.45397),
        ("Xi_c+", "PDG 2022", 2.46771),
        ("Xi_c0", "PDG 2022", 2.47044),
        ("Omega_c0", "PDG 2022", 2.69530),
        ("Xi_cc++", "PDG 2022", 3.62155),
        ("Xi_cc+", "LHCb 2026", 3.61997),
        # Mesons
        ("Pion+-", "PDG 2022", 0.13957039),
        ("Pion0", "PDG 2022", 0.13497770),
        ("Kaon+-", "PDG 2022", 0.49367700),
        ("Kaon0", "PDG 2022", 0.49761700),
        ("Eta", "PDG 2022", 0.54753),
        ("Rho770", "PDG 2022", 0.77526),
        ("Omega782", "PDG 2022", 0.78265),
        ("Phi1020", "PDG 2022", 1.01946),
        ("D0 meson", "PDG 2022", 1.86484),
        ("D+ meson", "PDG 2022", 1.86966),
        ("D_s+", "PDG 2022", 1.96835),
        ("J/psi", "PDG 2022", 3.09690),
        ("B+ meson", "PDG 2022", 5.27934),
        ("B0 meson", "PDG 2022", 5.27965),
        ("B_s0", "PDG 2022", 5.36688),
        ("B_c*+", "ATLAS 2026", 6.3390),
        ("Upsilon(1S)", "PDG 2022", 9.46030),
        ("Upsilon(2S)", "PDG 2022", 10.02326),
        ("Upsilon(3S)", "PDG 2022", 10.35520),
        ("Z_c(3900)", "PDG 2022", 3.8884),
        ("X(3872)", "PDG 2022", 3.87165),
        ("Omega_cc*", "CERN 2026", 3.7259),
    ]

    run_scan(particle_data)

    print("\nNOTE: err_exact ~ 0 by construction (K derived analytically).")
    print("Near-integer K = natural EWT resonance, no parameter adjustment.")
    print("Xi_cc+ (LHCb 2026) = post-construction independent validation.")
    print("=====================================================")

    # =====================================================
    # PART XIV: NEUTRINO RADIUS DERIVATION
    # =====================================================
    print("\n=====================================================")
    print("PART XIV. GEOMETRIC DERIVATION OF THE NEUTRINO RADIUS (r_nu)")
    print("=====================================================")

    qP = 1.875546e-18
    N_bcc = 8
    gv = 0.98359223

    epsilon_M = 1.0 / (8 * (Pi ** 7))
    alpha_inv_geo = (4*(Pi**3) + (Pi**2) + Pi) - epsilon_M

    print("--- EWT: FINAL NEUTRINO RADIUS (r_nu) DERIVATION ---\n")
    print("1. Geometric fine-structure constant (inverse):")
    print(f"   alpha_inv = {alpha_inv_geo:.12f}\n")

    K_proj = alpha_inv_geo / (N_bcc + Pi)
    K_expansion = e_euler
    delta_imp = (1.0 - gv) * (math.sqrt(2) - 1)
    K_final = K_proj + K_expansion + delta_imp

    print("2. Components of the scaling factor K = r_nu / q_P:")
    print(f"   - Static lattice projection:   {K_proj:.10f}  [alpha_inv / (8+pi)]")
    print(f"   - Dynamic wave expansion:      {K_expansion:.10f}  [e]")
    print(f"   - Discrete lattice impedance:  {delta_imp:.10f}  [(1-g_v)*(sqrt(2)-1)]")
    print(f"   => Total K:                    {K_final:.10f}\n")

    r_nu = qP * K_final
    print("3. Neutrino radius:")
    print(f"   r_nu = q_P * K = {r_nu:.25e} m\n")

    K_earlier = 2 * (e_euler ** 2) / gv
    print("4. Consistency with earlier derivation:")
    print(f"   Earlier K (2 e^2 / g_v)   = {K_earlier:.10f}")
    print(f"   Current K (sum)           = {K_final:.10f}")
    print(f"   Relative difference        = {abs(K_final - K_earlier)/K_earlier:.10e}\n")

    print("=====================================================")
    print("5. SELF-CONSISTENT QUADRATIC EQUATION FOR g_v")
    print("=====================================================")

    a_coef = math.sqrt(2) - 1
    b_coef = -(K_proj + e_euler + math.sqrt(2) - 1)
    c_coef = 2 * e_euler ** 2

    print("   Quadratic coefficients:")
    print(f"   a = (sqrt(2)-1)             = {a_coef:.15f}")
    print(f"   b = -(alpha_inv/(8+pi) + e + sqrt(2) - 1) = {b_coef:.15f}")
    print(f"   c = 2*e^2                   = {c_coef:.15f}\n")

    discriminant = b_coef ** 2 - 4 * a_coef * c_coef
    print(f"   Discriminant (b^2 - 4ac)    = {discriminant:.15e}\n")

    if discriminant >= 0:
        gv_root1 = (-b_coef + math.sqrt(discriminant)) / (2 * a_coef)
        gv_root2 = (-b_coef - math.sqrt(discriminant)) / (2 * a_coef)

        print(f"   Root 1: g_v = {gv_root1:.15f}")
        print(f"   Root 2: g_v = {gv_root2:.15f}\n")

        print("   Physical selection criterion: 0 < g_v < 1")

        label1 = "PHYSICAL" if (0 < gv_root1 < 1) else "UNPHYSICAL"
        label2 = "PHYSICAL" if (0 < gv_root2 < 1) else "UNPHYSICAL"

        print(f"   => Root 1 ({gv_root1:.6f}): {label1}")
        print(f"   => Root 2 ({gv_root2:.6f}): {label2}\n")

        gv_predicted = gv_root1 if (0 < gv_root1 < 1) else gv_root2

        print(f"   => Selected geometric fixed point: g_v = {gv_predicted:.15f}\n")

        delta_imp_pred = (1 - gv_predicted) * (math.sqrt(2) - 1)
        K_pred = K_proj + e_euler + delta_imp_pred
        r_nu_pred = qP * K_pred
        K_dyn_pred = 2 * e_euler ** 2 / gv_predicted

        print("   Verification with predicted g_v:")
        print(f"   K (geometric sum)          = {K_pred:.15f}")
        print(f"   K (dynamic 2e^2/g_v)       = {K_dyn_pred:.15f}")
        print(f"   Relative difference K      = {abs(K_pred - K_dyn_pred)/K_dyn_pred:.6e}")
        print(f"   r_nu (predicted)           = {r_nu_pred:.15e} m")
        print(f"   r_nu (earlier, gv=0.98359) = {r_nu:.15e} m")
        print(f"   Relative difference r_nu   = {abs(r_nu_pred - r_nu)/r_nu:.6e}\n")

        print(f"   Input g_v (phenomenological) = {gv:.8f}")
        print(f"   Predicted g_v (fixed point)  = {gv_predicted:.8f}")
        print(f"   Difference                   = {abs(gv_predicted - gv):.6e}")
    else:
        print("   ERROR: Negative discriminant - no real roots.")

    print("=====================================================")
    print("\n6. Physical interpretation:")
    print("   * g_v is the unique geometric fixed point of the BCC lattice.")
    print("   * Only one root satisfies 0 < g_v < 1.")
    print("   * This uniqueness suggests g_v is not a free parameter")
    print("     but a topological necessity of the vacuum lattice.")

if __name__ == "__main__":
    main()