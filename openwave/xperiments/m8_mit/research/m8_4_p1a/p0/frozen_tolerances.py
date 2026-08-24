"""Frozen qualification tolerances for P0.

These values were derived during the first engineering run (RUN1) from:
  - M8.5-B one-form precedent: 3.58e-10 admissible, O(1) rejection
  - R_0 collapse residual: 2.84e-14 real, 0.0 imaginary
  - RBF-FD weight covariance on S³/2I: ~1e-10 (60 seeds, k=110)
  - P0.4 oracle floors: ~1e-9 for d_ρ ≤ 4, ~5e-4 for d_ρ = 7

They are frozen HERE, in this file, BEFORE the qualification rerun.
The file's SHA-256 hash is recorded in the qualification note to prove
that no tolerance was adjusted after seeing a rerun residual.
"""

TOLERANCES = {
    "p0_0_homomorphism": 1e-8,
    "p0_0_unitarity": 1e-8,
    "p0_0_trace": 1e-6,
    "p0_0_identity": 1e-12,

    "p0_2_collapse_real": 1e-10,
    "p0_2_collapse_imag": 1e-14,

    "p0_3_eigenvalue_window": 2.0,

    "p0_4_admissible_floor": 1e-3,
    "p0_4_rejection_floor": 1.0,
    "p0_4_min_separation": 1e3,
}

GRID_PARAMS = {
    "n_seeds": 60,
    "k_stencil": 110,
    "rbf_m": 7,
    "rbf_p": 4,
}
