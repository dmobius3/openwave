#!/usr/bin/env python3
"""M4/EWT - Newtonian Force Mutual Consistency Check via EMC Push-Out

OpenWave criterion:
    Gravity: Newton limit (GEM)
    Test: Mutual consistency of monopole amplitude A = 2GM/c², push-out
          coupling K_emc = c⁴/(16πG), and 4π spatial overlap coefficient.

Physics & Structural Identification:
1. Exact Domain Integration: Mapped radial integration r = R/(1-t) over
   t in [0, 1) evaluates the pre-computed angular result J(r) = 4π/r² for r >= R.
2. Algebraic Identity & Consistency Gate: G, c, M, and R cancel identically
   in F_emc ≡ F_Newton. This script serves as a consistency check discriminating
   against incorrect normalization parameters (e.g., angular factor 4π, amplitude
   scaling, or stress-tensor coupling 16πG).
3. Residual Floor: The reported ~1.2e-5 % difference represents accumulated
   floating-point roundoff of the midpoint summation over a constant integrand,
   not a physical residual.
"""

import math


def compute_overlap_integral_exact_domain(A1, A2, R, num_pts=10000):
    """Evaluates the spatial overlap integral I(R) over r in [R, inf)

    using coordinate transformation t in [0, 1).
    """

    def overlap_integral_at(dist):
        dt = 1.0 / num_pts
        i_sum = 0.0

        for i in range(num_pts):
            t = (i + 0.5) * dt
            r = dist / (1.0 - t)
            dr_dt = dist / ((1.0 - t) * (1.0 - t))

            # Angularly integrated overlap for r >= R
            angular_integral = 4.0 * math.pi / (r * r)
            i_sum += angular_integral * dr_dt * dt

        return A1 * A2 * i_sum

    dR = R * 1e-6
    i_plus = overlap_integral_at(R + dR)
    i_minus = overlap_integral_at(R - dR)

    f_numeric = -(i_plus - i_minus) / (2.0 * dR)
    return f_numeric


def main():
    print("[1/4] Loading physical parameters...")
    c = 299792458.0
    G_EWT_geo = 6.674336927110799e-11
    M1, M2 = 1.989e30, 5.972e24
    R = 1.495978707e11

    A1 = 2.0 * G_EWT_geo * M1 / (c * c)
    A2 = 2.0 * G_EWT_geo * M2 / (c * c)

    print(f"    c         = {c:.3f} m/s")
    print(f"    G_EWT_geo = {G_EWT_geo:.15e} m^3 kg^-1 s^-2")
    print(f"    A1 (Sun)  = {A1:.3f} m")
    print(f"    A2 (Earth)= {A2:.3f} m")
    print(f"    R (1 AU)  = {R:.3e} m")

    print("\n[2/4] Performing mapped 3D spatial integration over [R, inf)...")
    F_numeric = compute_overlap_integral_exact_domain(A1, A2, R)

    print("\n[3/4] Converting field integral to physical force...")
    K_emc = c**4 / (16.0 * math.pi * G_EWT_geo)
    F_emc = K_emc * F_numeric

    print("\n[4/4] Comparing with Newton's law...")
    F_newton = G_EWT_geo * M1 * M2 / (R * R)
    rel_diff = abs(F_emc - F_newton) / F_newton * 100.0

    print(f"    F_EMC (numerical) = {F_emc:.12e} N")
    print(f"    F_Newton          = {F_newton:.12e} N")
    print(f"    Rel. diff.        = {rel_diff:.12e}%")

    # Gate threshold set to 1e-3% to verify mutual algebraic consistency
    if rel_diff < 1e-3:
        print("    RESULT: PASS (Mutual consistency of A, K_emc, and 4π overlap verified)")
    else:
        print("    RESULT: FAIL")


if __name__ == "__main__":
    main()
