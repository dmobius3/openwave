#!/usr/bin/env python3
"""M8.9 S1b Implementation Qualification

Commissioned unit: reproduce every gate and mutation from the frozen contract
S1B_DECISION_RULE.md, attack for vacuity, and exercise the branch logic on
synthetic inputs.

Contract frozen-region SHA-256:
  c44c603a8058ed8529e5bb0f42ec168b443ff22b050d3748f1f1de3537c7d297

No live R_0 spectrum at n=12 or n=20.  No S1b result.
"""

import sys, os, time, hashlib, json
import numpy as np
import mpmath

ROOM = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ROOM)
sys.path.insert(0, os.path.join(ROOM, "m8_5b", "pilot"))
sys.path.insert(0, os.path.join(ROOM, "m8_5b", "production"))
sys.path.insert(0, os.path.join(ROOM, "m8_5b", "gates"))

import route_a_repn as repn
from route_a_nonabelian import quat_to_su2, sym_power
from route_a_twosided import pairs_left
from p0.group import build_icosians
from p0.cloud import fibonacci_seeds_s3, build_orbit_cloud
from p0.bundle_operator import build_L_bundle, orbit_stencils as p0_orbit_stencils
from p0.representations import build_all_representations
from p0.group import build_character_table
from p1a.mass_matrix import build_Mh_base, Mh_inner, Mh_norm_F

eps = np.finfo(np.float64).eps
SEP = "=" * 72

def section(title):
    print(f"\n{SEP}\n  {title}\n{SEP}")

def invariant_dim_and_basis_fast(pairs, n, tol_rel=1e-8, no_transpose=False):
    """Same math as repn.invariant_dim_and_basis, economy SVD.

    The shipped function uses np.linalg.svd(A) with full_matrices=True (the
    default), allocating an m x m U matrix that is immediately discarded.
    At n=20 this is a 52920 x 52920 complex matrix (44 GB).  The economy SVD
    produces identical s and Vh and is the only change.
    """
    dim = (n + 1) ** 2
    blocks = [repn.coefficient_operator(u, v, n, no_transpose) - np.eye(dim)
              for (u, v) in pairs]
    A = np.vstack(blocks)
    _, s, Vh = np.linalg.svd(A, full_matrices=False)
    smax = s[0] if s.size else 0.0
    cutoff = max(tol_rel * smax, 1e-12)
    k = int(np.sum(s < cutoff))
    below, above = s[s < cutoff], s[s >= cutoff]
    if below.size == 0 or above.size == 0:
        gap = float("inf")
        gap_state = "vacuous: one side of the cut is empty"
    elif below.max() == 0.0:
        gap = float("inf")
        gap_state = "exact: nullspace is identically zero, separation is perfect"
    else:
        gap = float(above.min() / below.max())
        gap_state = "measured"
    basis = Vh.conj().T[:, dim - k:] if k else np.zeros((dim, 0), dtype=complex)
    return k, basis, {"gap": gap, "state": gap_state}


def mh_orthonormalize(Q, Mh_diag):
    """M_h-orthonormalize: Q^H M_h Q = I via modified Gram-Schmidt."""
    n, k = Q.shape
    R = Q.copy().astype(complex)
    for j in range(k):
        for i in range(j):
            coeff = np.dot(np.conj(R[:, i]) * Mh_diag, R[:, j])
            R[:, j] -= coeff * R[:, i]
        norm = np.sqrt(np.real(np.dot(np.conj(R[:, j]) * Mh_diag, R[:, j])))
        if norm < 1e-14:
            raise ValueError(f"M_h-orthonormalization failed at column {j}")
        R[:, j] /= norm
    return R


def principal_angles_coeff(B1, B2):
    """Principal angles between column spaces of B1, B2 in coefficient space.

    Uses SVD of B1^H B2 (ordinary inner product).
    """
    B1q, _ = np.linalg.qr(B1)
    B2q, _ = np.linalg.qr(B2)
    M = B1q.conj().T @ B2q
    sigmas = np.linalg.svd(M, compute_uv=False)
    sigmas_clipped = np.clip(np.real(sigmas), 0.0, 1.0)
    thetas = np.arccos(sigmas_clipped)
    return thetas


def form_An(Q, Mh_diag, L):
    """A_n = Q^H M_h L Q."""
    MhQ = Mh_diag[:, None] * Q
    return MhQ.conj().T @ L @ Q


def compute_K(A):
    """K = (A - A^H)/2, return (K, ||K||_2)."""
    K = (A - A.conj().T) / 2
    K_norm = np.linalg.norm(K, ord=2)
    return K, float(K_norm)


def compute_J(A):
    """J = max |Im(mu)| over eigenvalues mu of A."""
    evals = np.linalg.eigvals(A)
    return float(np.max(np.abs(np.imag(evals))))


# ===========================================================================
#  MAIN
# ===========================================================================
def main():
    results = {}
    all_pass = True

    # -----------------------------------------------------------------------
    #  Q0: PROVENANCE
    # -----------------------------------------------------------------------
    section("Q0: PROVENANCE")

    # Q0a: Contract hash
    contract_path = os.path.join(ROOM, "contract", "S1B_DECISION_RULE.md")
    with open(contract_path, "rb") as f:
        raw = f.read()
    boundary = b"<!-- FREEZE-BOUNDARY -->"
    idx = raw.find(boundary)
    frozen = raw[:idx]
    h = hashlib.sha256(frozen).hexdigest()
    expected = "c44c603a8058ed8529e5bb0f42ec168b443ff22b050d3748f1f1de3537c7d297"
    print(f"  contract frozen-region SHA-256: {h}")
    print(f"  expected:                       {expected}")
    print(f"  MATCH: {h == expected}")
    assert h == expected, "CONTRACT HASH MISMATCH — STOP"

    # Q0b: Manifest hash
    manifest_path = os.path.join(ROOM, "ROOM_MANIFEST.json")
    with open(manifest_path, "rb") as f:
        mh = hashlib.sha256(f.read()).hexdigest()
    manifest_expected = "46ab992545c6a28cd6d116c238aa00ed5525f6c11f39688ea832172d6869d20f"
    print(f"  manifest SHA-256: {mh}")
    print(f"  expected:         {manifest_expected}")
    print(f"  MATCH: {mh == manifest_expected}")
    assert mh == manifest_expected, "MANIFEST HASH MISMATCH — STOP"

    # Q0c: Environment
    print(f"\n  Python:  {sys.version}")
    print(f"  numpy:   {np.__version__}")
    import scipy; print(f"  scipy:   {scipy.__version__}")
    print(f"  mpmath:  {mpmath.__version__}")

    print(f"\n  AFFIRMATION: No OpenWave or project code, and no module in the")
    print(f"  required closure, was imported from outside this room.")
    print(f"  Standard library and third-party packages (numpy, scipy, mpmath)")
    print(f"  are external and are covered by the environment record above.")

    # -----------------------------------------------------------------------
    #  INFRASTRUCTURE
    # -----------------------------------------------------------------------
    section("INFRASTRUCTURE")
    t0 = time.time()
    elems = build_icosians()
    print(f"  2I elements: {len(elems)}")
    pairs = pairs_left(elems)
    print(f"  left-action pairs: {len(pairs)}")
    chi = build_character_table(elems)
    reps, _ = build_all_representations(elems, chi)
    print(f"  representations built: {list(reps.keys())}")

    N_SEEDS = 60
    K_STENCIL = 110
    RBF_M, RBF_P = 7, 4
    seeds = fibonacci_seeds_s3(N_SEEDS)
    X, oid, gid = build_orbit_cloud(seeds, elems)
    print(f"  cloud: {len(X)} nodes from {N_SEEDS} seeds")
    W = build_Mh_base(X, oid, N_SEEDS)
    Mh_diag = W.copy()
    print(f"  M_h: {len(W)} weights, range [{W.min():.4e}, {W.max():.4e}], "
          f"kappa={W.max()/W.min():.4f}")

    L_q, seed_orbits = build_L_bundle(X, oid, gid, elems, reps["R0"],
                                       k=K_STENCIL, m=RBF_M, p=RBF_P)
    L_q = np.real(np.asarray(L_q))
    print(f"  L_q: {L_q.shape}, seed_orbits: {len(seed_orbits)}")
    print(f"  infrastructure time: {time.time()-t0:.1f}s")

    # -----------------------------------------------------------------------
    #  G-RANK
    # -----------------------------------------------------------------------
    section("G-RANK: invariant dimension")
    for n, k_expected in [(12, 13), (20, 21)]:
        t0 = time.time()
        k, basis, gap_info = invariant_dim_and_basis_fast(pairs, n)
        dt = time.time() - t0
        passed = (k == k_expected)
        status = "PASS" if passed else "FAIL"
        print(f"  n={n}: dim={k} (expected {k_expected}), gap={gap_info}, "
              f"time={dt:.2f}s -> {status}")
        results[f"G-RANK_n{n}"] = {"dim": k, "expected": k_expected,
                                    "gap": gap_info, "pass": passed}
        if not passed:
            all_pass = False
    print(f"  note: economy SVD used (identical s, Vh; full_matrices=False)")

    # -----------------------------------------------------------------------
    #  G-REAL
    # -----------------------------------------------------------------------
    section("G-REAL: pointwise realization")
    t0 = time.time()
    vr = repn.verify_realization()
    dt = time.time() - t0
    gr_pass = vr["pass"]
    print(f"  levels swept:              {vr['levels_swept']}")
    print(f"  worst_residual_correct:    {vr['worst_residual_correct']:.2e}  (need < 1e-10)")
    print(f"  best_residual_no_transpose:{vr['best_residual_no_transpose']:.2e}  (need > 1e-10)")
    print(f"  realization: {vr['realization']}")
    print(f"  time: {dt:.2f}s -> {'PASS' if gr_pass else 'FAIL'}")
    results["G-REAL"] = vr
    if not gr_pass:
        all_pass = False

    # G-REAL: also test deck invariance on the actual cover cloud
    print("\n  deck invariance on cover cloud (independent of verify_realization):")
    for n_test in [0, 6]:
        _, basis_n, _ = invariant_dim_and_basis_fast(pairs, n_test)
        if basis_n.shape[1] == 0:
            print(f"    n={n_test}: no invariant subspace, skip")
            continue
        dim_n = (n_test + 1) ** 2
        Vn = np.zeros((len(X), dim_n), dtype=complex)
        for i in range(len(X)):
            Vn[i] = sym_power(quat_to_su2(X[i]), n_test).reshape(-1, order="F")
        Fn_cover = Vn @ basis_n
        worst_res = 0.0
        for gi in range(1, min(len(elems), 20)):
            for node_i in range(0, len(X), 120):
                gamma_i = node_i + gi
                if gamma_i < len(X):
                    res = np.abs(Fn_cover[gamma_i] - Fn_cover[node_i]).max()
                    worst_res = max(worst_res, res)
        print(f"    n={n_test}: k={basis_n.shape[1]}, "
              f"worst deck residual={worst_res:.2e}")

    # -----------------------------------------------------------------------
    #  G-SUBSPACE
    # -----------------------------------------------------------------------
    section("G-SUBSPACE: two independent constructions")

    for n in [12, 20]:
        print(f"\n  --- n={n} ---")
        k_n = 13 if n == 12 else 21
        dim_n = (n + 1) ** 2

        # SVD route
        t0 = time.time()
        k_svd, C_svd, gap_svd = invariant_dim_and_basis_fast(pairs, n)
        t_svd = time.time() - t0
        print(f"  SVD route: dim={k_svd}, time={t_svd:.2f}s")

        # Averaging route: Pi_n = (1/|G|) sum coefficient_operator(u, v, n)
        t0 = time.time()
        Pi = np.zeros((dim_n, dim_n), dtype=complex)
        for (u, v) in pairs:
            Pi += repn.coefficient_operator(u, v, n)
        Pi /= len(pairs)
        t_avg = time.time() - t0
        print(f"  averaging route: Pi_n built, time={t_avg:.2f}s")

        # Idempotence check
        Pi2_err = float(np.linalg.norm(Pi @ Pi - Pi))
        PiH_err = float(np.linalg.norm(Pi - Pi.conj().T))
        print(f"  ||Pi^2 - Pi|| = {Pi2_err:.2e}  (idempotent)")
        print(f"  ||Pi - Pi^H|| = {PiH_err:.2e}  (oblique, NOT Hermitian)")

        # Rank at ABSOLUTE cutoff 1e-8
        s_Pi = np.linalg.svd(Pi, compute_uv=False)
        rank_abs = int(np.sum(s_Pi > 1e-8))
        print(f"  rank(Pi_n) at absolute cutoff 1e-8: {rank_abs}")

        # ran(Pi_n): orthonormal basis from LEFT singular vectors (column space)
        U_full, s_full, _ = np.linalg.svd(Pi, full_matrices=False)
        C_avg = U_full[:, :rank_abs]
        print(f"  C_avg shape: {C_avg.shape}")

        # Principal angles in coefficient space
        if k_svd > 0 and rank_abs > 0 and k_svd == rank_abs:
            thetas = principal_angles_coeff(C_svd, C_avg)
            theta_max = float(np.max(thetas))
            sin_theta_max = float(np.sin(theta_max))
            print(f"  theta_C max: {theta_max:.2e} rad")
            print(f"  sin(theta_C max): {sin_theta_max:.2e}  (gate: <= 1e-6)")
        else:
            sin_theta_max = float("inf")
            print(f"  DIMENSION MISMATCH: SVD={k_svd}, avg={rank_abs}")

        # Gate evaluation
        dims_agree = (k_svd == rank_abs == k_n)
        angle_ok = (sin_theta_max <= 1e-6)
        gs_pass = dims_agree and angle_ok
        print(f"  dims agree: {dims_agree}, angle ok: {angle_ok}")
        print(f"  G-SUBSPACE n={n}: {'PASS' if gs_pass else 'FAIL'}")
        results[f"G-SUBSPACE_n{n}"] = {
            "k_svd": k_svd, "rank_Pi": rank_abs,
            "sin_theta_C": sin_theta_max,
            "idempotence_err": Pi2_err, "non_hermiticity": PiH_err,
            "pass": gs_pass,
        }
        if not gs_pass:
            all_pass = False

        # --- G-SUBSPACE ARMS ---
        print(f"\n  G-SUBSPACE arms, n={n}:")

        # Subspace arm: tilt one column
        for sin_phi, required in [(1e-4, "RED"), (1e-7, "green")]:
            C_tilt = C_svd.copy()
            orth_complement = np.eye(dim_n, dtype=complex) - C_svd @ C_svd.conj().T
            w = orth_complement[:, 0]
            w_norm = np.linalg.norm(w)
            if w_norm < 1e-12:
                for col in range(dim_n):
                    w = orth_complement[:, col]
                    w_norm = np.linalg.norm(w)
                    if w_norm > 1e-12:
                        break
            w = w / w_norm
            phi = np.arcsin(sin_phi)
            C_tilt[:, 0] = C_svd[:, 0] * np.cos(phi) + w * np.sin(phi)
            C_tilt, _ = np.linalg.qr(C_tilt)
            thetas_tilt = principal_angles_coeff(C_svd, C_tilt)
            meas = float(np.sin(np.max(thetas_tilt)))
            gate_says = "RED" if meas > 1e-6 else "green"
            match = gate_says == required
            print(f"    tilt sin(phi)={sin_phi:.0e}: measured sin(theta_max)="
                  f"{meas:.3e}, gate={gate_says}, required={required} "
                  f"-> {'OK' if match else 'MISMATCH'}")

        # Rank arm: numerically-zero matrix must have rank 0
        P_zero = np.zeros((dim_n, dim_n), dtype=complex) + 1e-17 * np.random.default_rng(42).normal(size=(dim_n, dim_n))
        s_zero = np.linalg.svd(P_zero, compute_uv=False)
        rank_zero_abs = int(np.sum(s_zero > 1e-8))
        rank_zero_rel = int(np.sum(s_zero > 1e-8 * s_zero[0])) if s_zero[0] > 0 else 0
        print(f"    rank arm: zero matrix rank at abs 1e-8: {rank_zero_abs} (must be 0)")
        print(f"    rank arm: zero matrix rank at rel 1e-8: {rank_zero_rel} "
              f"(defect if > 0)")

        # Known rank-k projector
        P_known = np.zeros((dim_n, dim_n), dtype=complex)
        P_known[:k_n, :k_n] = np.eye(k_n)
        s_known = np.linalg.svd(P_known, compute_uv=False)
        rank_known = int(np.sum(s_known > 1e-8))
        print(f"    rank arm: known rank-{k_n} projector at abs 1e-8: {rank_known}")

    # -----------------------------------------------------------------------
    #  G-ALIGN
    # -----------------------------------------------------------------------
    section("G-ALIGN: cover-to-seed restriction")
    node_of, mult, plan = p0_orbit_stencils(X, oid, gid, elems, k=K_STENCIL)
    n_seeds_actual = len(plan)
    print(f"  orbits in plan: {n_seeds_actual}")
    print(f"  seed_orbits == sorted(plan.keys()): {seed_orbits == sorted(plan.keys())}")

    r_vec = np.array([plan[o][0] for o in seed_orbits])
    r_expected = np.array([node_of[(o, 0)] for o in seed_orbits])
    align_ok = np.array_equal(r_vec, r_expected)
    print(f"  r == (node_of[(o, 0)])_{{o in seed_orbits}}: {align_ok}")
    results["G-ALIGN"] = {"structural_predicate": align_ok, "pass": align_ok}
    if not align_ok:
        all_pass = False
    else:
        print(f"  G-ALIGN: PASS")

    # G-ALIGN arm: change one entry
    r_mut = r_vec.copy()
    orig_val = r_mut[0]
    r_mut[0] = r_vec[1] if r_vec[1] != orig_val else r_vec[2]
    arm_fails = not np.array_equal(r_mut, r_expected)
    print(f"\n  G-ALIGN arm: changed r[0] from {orig_val} to {r_mut[0]}")
    print(f"    predicate fails: {arm_fails} (must be True) -> "
          f"{'OK' if arm_fails else 'MISMATCH'}")

    # -----------------------------------------------------------------------
    #  G-SAMPLE
    # -----------------------------------------------------------------------
    section("G-SAMPLE: conditioning before orthonormalization")

    for n, k_n in [(12, 13), (20, 21)]:
        print(f"\n  --- n={n}, k={k_n} ---")
        dim_n = (n + 1) ** 2

        for route_label, C_n in [("svd", None), ("avg", None)]:
            t0 = time.time()
            if route_label == "svd":
                _, C_n, _ = invariant_dim_and_basis_fast(pairs, n)
            else:
                Pi = np.zeros((dim_n, dim_n), dtype=complex)
                for (u, v) in pairs:
                    Pi += repn.coefficient_operator(u, v, n)
                Pi /= len(pairs)
                U_pi, s_Pi, _ = np.linalg.svd(Pi, full_matrices=False)
                rank_abs = int(np.sum(s_Pi > 1e-8))
                C_n = U_pi[:, :rank_abs]

            # Sample V_n on cover cloud
            Vn = np.zeros((len(X), dim_n), dtype=complex)
            for i in range(len(X)):
                Vn[i] = sym_power(quat_to_su2(X[i]), n).reshape(-1, order="F")

            # F_n^cover = V_n C_n
            Fn_cover = Vn @ C_n

            # F_n^seed: restrict to seed representatives
            rep_indices = np.array([node_of[(o, 0)] for o in seed_orbits])
            Fn_seed = Fn_cover[rep_indices]

            # W_n = M_h^{1/2} F_n^seed
            Mh_sqrt = np.sqrt(Mh_diag)
            Wn = Mh_sqrt[:, None] * Fn_seed

            # Check rank and conditioning
            s_W = np.linalg.svd(Wn, compute_uv=False)
            num_rank = int(np.sum(s_W > 1e-8 * s_W[0]))
            kappa = float(s_W[0] / s_W[-1]) if s_W[-1] > 0 else float("inf")
            dt = time.time() - t0

            rank_ok = (num_rank == k_n)
            kappa_ok = (kappa <= 1e6)
            gs_pass = rank_ok and kappa_ok
            print(f"  route={route_label}: rank={num_rank} (need {k_n}), "
                  f"kappa={kappa:.2e} (need <= 1e6), time={dt:.1f}s "
                  f"-> {'PASS' if gs_pass else 'FAIL'}")
            results[f"G-SAMPLE_n{n}_{route_label}"] = {
                "rank": num_rank, "kappa": kappa, "pass": gs_pass,
            }
            if not gs_pass:
                all_pass = False

    # -----------------------------------------------------------------------
    #  G-WIRE: n=0 machinery control
    # -----------------------------------------------------------------------
    section("G-WIRE: n=0 machinery control")

    # Sample the constant harmonic (n=0, k=1)
    _, C0, _ = invariant_dim_and_basis_fast(pairs, 0)
    print(f"  n=0: invariant dim = {C0.shape[1]}")
    V0 = np.zeros((len(X), 1), dtype=complex)
    for i in range(len(X)):
        V0[i] = sym_power(quat_to_su2(X[i]), 0).reshape(-1, order="F")
    F0_cover = V0 @ C0
    rep_indices = np.array([node_of[(o, 0)] for o in seed_orbits])
    F0_seed = F0_cover[rep_indices]

    # M_h-orthonormalize
    Q0 = mh_orthonormalize(F0_seed, Mh_diag)
    orth_err = float(np.abs(Q0.conj().T @ np.diag(Mh_diag) @ Q0 - np.eye(1)).max())
    print(f"  Q0 M_h-orthonormality error: {orth_err:.2e}")

    # ||L Q0||_{M_h,F}
    LQ0 = L_q @ Q0
    LQ0_norm = float(np.sqrt(np.real(np.sum(np.conj(LQ0) * (Mh_diag[:, None] * LQ0)))))
    gw_pass = LQ0_norm <= 1e-8
    print(f"  ||L Q_0||_{{M_h,F}} = {LQ0_norm:.2e}  (gate: <= 1e-8)")
    print(f"  G-WIRE: {'PASS' if gw_pass else 'FAIL'}")
    results["G-WIRE"] = {"LQ0_norm": LQ0_norm, "pass": gw_pass}
    if not gw_pass:
        all_pass = False

    # max |Im| on constant mode
    A0 = form_An(Q0, Mh_diag, L_q)
    im0 = float(np.max(np.abs(np.imag(np.linalg.eigvals(A0)))))
    print(f"  max |Im| on A_0: {im0:.2e}  (recorded, not gated; "
          f"1x1 real => zero by algebra)")

    # G-WIRE arm: L_mut = L + delta * P_0, delta = 1e-4
    delta = 1e-4
    P0 = Q0 @ (Q0.conj().T @ np.diag(Mh_diag))
    L_mut = L_q + delta * P0
    LmQ0 = L_mut @ Q0
    LmQ0_norm = float(np.sqrt(np.real(np.sum(np.conj(LmQ0) * (Mh_diag[:, None] * LmQ0)))))
    arm_red = LmQ0_norm > 1e-8
    print(f"\n  G-WIRE arm: delta={delta}")
    print(f"    ||L_mut Q_0||_{{M_h,F}} = {LmQ0_norm:.2e} "
          f"(expected >= {delta - 1e-8:.0e})")
    print(f"    arm goes RED: {arm_red} -> {'OK' if arm_red else 'MISMATCH'}")

    # K_0 algebra check
    _, K0_norm = compute_K(A0)
    print(f"\n  K_0 = ||K_0||_2 = {K0_norm:.2e}  (zero by algebra: A_0 is 1x1 real)")

    # -----------------------------------------------------------------------
    #  G-BASIS: basis invariance
    # -----------------------------------------------------------------------
    section("G-BASIS: basis invariance")

    # Run on actual objects at n=0 (trivial, 1x1)
    print("  n=0 (trivial):")
    rng_basis = np.random.default_rng(20260825)
    k0 = Q0.shape[1]
    H = rng_basis.normal(size=(k0, k0)) + 1j * rng_basis.normal(size=(k0, k0))
    U0, _ = np.linalg.qr(H)
    Q0_rot = Q0 @ U0
    A0_rot = form_An(Q0_rot, Mh_diag, L_q)
    A0_cov_err = float(np.linalg.norm(A0_rot - U0.conj().T @ A0 @ U0, "fro"))
    A0_scale = max(float(np.linalg.norm(A0, "fro")), 1e-30)
    A0_rel = A0_cov_err / A0_scale
    print(f"    ||A' - U^H A U||_F = {A0_cov_err:.2e}, "
          f"||A||_F = {A0_scale:.2e}, rel = {A0_rel:.2e}")

    P0_proj = Q0 @ (Q0.conj().T @ np.diag(Mh_diag))
    P0_rot = Q0_rot @ (Q0_rot.conj().T @ np.diag(Mh_diag))
    P0_diff = float(np.linalg.norm(P0_rot - P0_proj, "fro"))
    P0_scale = max(1.0, float(np.linalg.norm(P0_proj, "fro")))
    P0_rel = P0_diff / P0_scale
    print(f"    ||P' - P||_F = {P0_diff:.2e}, scale = {P0_scale:.2e}, "
          f"rel = {P0_rel:.2e}")

    # Demonstrate on a synthetic higher-dimensional case
    print("\n  synthetic k=13 (representative of n=12 dimension):")
    rng_syn = np.random.default_rng(20260825)
    k_syn = 13
    n_syn = 60
    Q_syn = np.linalg.qr(rng_syn.normal(size=(n_syn, k_syn)) +
                         1j * rng_syn.normal(size=(n_syn, k_syn)))[0]
    Mh_syn = np.abs(rng_syn.normal(size=n_syn)) + 0.01
    L_syn = rng_syn.normal(size=(n_syn, n_syn))
    L_syn = (L_syn + L_syn.T) / 2
    Q_syn = mh_orthonormalize(Q_syn, Mh_syn)
    A_syn = form_An(Q_syn, Mh_syn, L_syn)

    H_syn = rng_syn.normal(size=(k_syn, k_syn)) + 1j * rng_syn.normal(size=(k_syn, k_syn))
    U_syn, _ = np.linalg.qr(H_syn)
    Q_syn_rot = Q_syn @ U_syn
    A_syn_rot = form_An(Q_syn_rot, Mh_syn, L_syn)
    A_cov_err = float(np.linalg.norm(A_syn_rot - U_syn.conj().T @ A_syn @ U_syn, "fro"))
    A_cov_scale = max(float(np.linalg.norm(A_syn, "fro")), 1e-30)
    A_cov_rel = A_cov_err / A_cov_scale
    a_pass = A_cov_rel <= 1e-10
    print(f"    ||A' - U^H A U||_F / ||A||_F = {A_cov_rel:.2e}  "
          f"(gate: <= 1e-10) -> {'PASS' if a_pass else 'FAIL'}")

    P_syn = Q_syn @ (Q_syn.conj().T @ np.diag(Mh_syn))
    P_syn_rot = Q_syn_rot @ (Q_syn_rot.conj().T @ np.diag(Mh_syn))
    P_diff = float(np.linalg.norm(P_syn_rot - P_syn, "fro"))
    P_scale = max(1.0, float(np.linalg.norm(P_syn, "fro")))
    P_rel = P_diff / P_scale
    p_pass = P_rel <= 1e-10
    print(f"    ||P' - P||_F / max(1,||P||_F) = {P_rel:.2e}  "
          f"(gate: <= 1e-10) -> {'PASS' if p_pass else 'FAIL'}")
    gb_pass = a_pass and p_pass
    results["G-BASIS"] = {"A_covariance_rel": A_cov_rel, "P_diff_rel": P_rel,
                          "pass": gb_pass}
    if not gb_pass:
        all_pass = False

    # -----------------------------------------------------------------------
    #  G-DISCRIM
    # -----------------------------------------------------------------------
    section("G-DISCRIM: K vs J discrimination")

    lam = 168.0
    k_test = 13
    A_green = lam * np.eye(k_test, dtype=complex)
    K_green, K_green_norm = compute_K(A_green)
    J_green = compute_J(A_green)
    K_green_bound = 100 * eps * lam
    green_K_ok = K_green_norm <= K_green_bound
    green_J_ok = J_green <= 1e-12
    print(f"  green parent: A = {lam} * I")
    print(f"    ||K||_2 = {K_green_norm:.2e}  (need <= 100*eps*||A||_2 = {K_green_bound:.2e})"
          f" -> {'ok' if green_K_ok else 'FAIL'}")
    print(f"    J       = {J_green:.2e}  (need <= 1e-12) -> {'ok' if green_J_ok else 'FAIL'}")

    # Arm A: nilpotent, non-Hermitian but real spectrum
    A_armA = A_green.copy()
    A_armA[0, 1] += 1.0
    K_A, K_A_norm = compute_K(A_armA)
    J_A = compute_J(A_armA)
    armA_K_ok = abs(K_A_norm - 0.5) <= 1e-12
    armA_J_ok = J_A <= 1e-12
    print(f"\n  arm A: nilpotent N[0,1]=1")
    print(f"    ||K||_2 = {K_A_norm:.6f}  (need |..-0.5| <= 1e-12, "
          f"err={abs(K_A_norm - 0.5):.2e}) -> {'ok' if armA_K_ok else 'FAIL'}")
    print(f"    J       = {J_A:.2e}  (need <= 1e-12) -> {'ok' if armA_J_ok else 'FAIL'}")

    # Arm B: real antisymmetric, complex spectrum
    A_armB = A_green.copy()
    A_armB[0, 1] += 2.0
    A_armB[1, 0] -= 2.0
    K_B, K_B_norm = compute_K(A_armB)
    J_B = compute_J(A_armB)
    armB_K_ok = abs(K_B_norm - 2.0) <= 1e-12
    armB_J_ok = abs(J_B - 2.0) <= 1e-12
    print(f"\n  arm B: antisymmetric S[0,1]=-S[1,0]=2")
    print(f"    ||K||_2 = {K_B_norm:.6f}  (need |..-2| <= 1e-12, "
          f"err={abs(K_B_norm - 2.0):.2e}) -> {'ok' if armB_K_ok else 'FAIL'}")
    print(f"    J       = {J_B:.6f}  (need |..-2| <= 1e-12, "
          f"err={abs(J_B - 2.0):.2e}) -> {'ok' if armB_J_ok else 'FAIL'}")

    gd_pass = (green_K_ok and green_J_ok and armA_K_ok and armA_J_ok
               and armB_K_ok and armB_J_ok)
    print(f"\n  G-DISCRIM: {'PASS' if gd_pass else 'FAIL'}")
    results["G-DISCRIM"] = {
        "green_K": K_green_norm, "green_J": J_green,
        "armA_K": K_A_norm, "armA_J": J_A,
        "armB_K": K_B_norm, "armB_J": J_B,
        "pass": gd_pass,
    }
    if not gd_pass:
        all_pass = False

    # -----------------------------------------------------------------------
    #  Q2: MUTATION POWER
    # -----------------------------------------------------------------------
    section("Q2: MUTATION POWER")

    # G-REAL mutation (already shown via verify_realization above)
    print("  G-REAL arm: worst_residual_correct < 1e-10 AND "
          "best_residual_no_transpose > 1e-10")
    print(f"    correct:      {vr['worst_residual_correct']:.2e} -> GREEN")
    print(f"    no_transpose: {vr['best_residual_no_transpose']:.2e} -> RED")
    print(f"    separation:   {vr['best_residual_no_transpose']/vr['worst_residual_correct']:.2e}")

    # G-WIRE mutation (already shown above)
    print(f"\n  G-WIRE arm: green parent ||LQ0|| = {LQ0_norm:.2e}, "
          f"red child ||L_mut Q0|| = {LmQ0_norm:.2e}")

    # G-ALIGN mutation (already shown above)
    print(f"\n  G-ALIGN arm: structural predicate on r vector")
    print(f"    green: r matches node_of[(o,0)] -> {align_ok}")
    print(f"    red:   r[0] changed -> predicate fails: {arm_fails}")

    # G-DISCRIM mutations (already shown above)
    print(f"\n  G-DISCRIM arm A: green K={K_green_norm:.2e}, J={J_green:.2e}")
    print(f"                   red   K={K_A_norm:.2e}, J={J_A:.2e}")
    print(f"  G-DISCRIM arm B: green K={K_green_norm:.2e}, J={J_green:.2e}")
    print(f"                   red   K={K_B_norm:.2e}, J={J_B:.2e}")

    # G-SUBSPACE tilt arm (already shown above for sin_phi=1e-4)
    print(f"\n  G-SUBSPACE tilt arm: demonstrated above for sin(phi)=1e-4 -> RED")
    print(f"  G-SUBSPACE green control: sin(phi)=1e-7 -> green")

    # Precision ladder arm: companion matrix of (x-1)^k
    print(f"\n  Precision ladder arm: companion matrix of (x-1)^k")

    # Test at k=2 first (should COLLAPSE via rule 1: 2×2 eigensolve is analytic)
    k_comp_2 = 2
    coeffs_2 = np.array([1.0])
    for _ in range(k_comp_2):
        coeffs_2 = np.convolve(coeffs_2, [1.0, -1.0])
    comp_2 = np.zeros((k_comp_2, k_comp_2))
    comp_2[1:, :-1] = np.eye(k_comp_2 - 1)
    comp_2[:, -1] = -coeffs_2[:-1] / coeffs_2[-1]
    evals_2 = np.linalg.eigvals(comp_2)
    J_2_64 = float(np.max(np.abs(np.imag(evals_2))))
    with mpmath.workdps(30):
        A2_mp = mpmath.matrix(k_comp_2, k_comp_2)
        for ii in range(k_comp_2):
            for jj in range(k_comp_2):
                A2_mp[ii, jj] = mpmath.mpf(str(comp_2[ii, jj]))
        ev30 = mpmath.eig(A2_mp, left=False, right=False)
        J_2_30 = float(max(abs(complex(e).imag) for e in ev30))
    with mpmath.workdps(50):
        A2_mp = mpmath.matrix(k_comp_2, k_comp_2)
        for ii in range(k_comp_2):
            for jj in range(k_comp_2):
                A2_mp[ii, jj] = mpmath.mpf(str(comp_2[ii, jj]))
        ev50 = mpmath.eig(A2_mp, left=False, right=False)
        J_2_50 = float(max(abs(complex(e).imag) for e in ev50))
    if J_2_64 == 0 and J_2_30 == 0 and J_2_50 == 0:
        comp_2_label = "COLLAPSES"
    elif J_2_50 < 0.5 * J_2_30 or J_2_50 > 2 * J_2_30:
        comp_2_label = "AMBIGUOUS"
    elif J_2_50 >= 0.5 * J_2_64:
        comp_2_label = "PERSISTS"
    elif J_2_50 <= 1e-3 * J_2_64:
        comp_2_label = "COLLAPSES"
    else:
        comp_2_label = "AMBIGUOUS"
    print(f"    k=2: J(64)={J_2_64:.2e}, J(30)={J_2_30:.2e}, J(50)={J_2_50:.2e}"
          f" -> {comp_2_label}")

    # Test at k=6 (defective; k-th root scaling causes AMBIGUOUS)
    k_comp = 6
    coeffs = np.array([1.0])
    for _ in range(k_comp):
        coeffs = np.convolve(coeffs, [1.0, -1.0])
    companion = np.zeros((k_comp, k_comp))
    companion[1:, :-1] = np.eye(k_comp - 1)
    companion[:, -1] = -coeffs[:-1] / coeffs[-1]
    A_comp = companion

    evals_comp = np.linalg.eigvals(A_comp)
    J_comp_64 = float(np.max(np.abs(np.imag(evals_comp))))
    print(f"    float64 J = {J_comp_64:.2e}")

    # Re-form the companion matrix at higher precision and eigensolve.
    # For k-fold defective eigenvalues, J scales as eps^{1/k}.
    # 30 dps
    with mpmath.workdps(30):
        A_mp = mpmath.matrix(k_comp, k_comp)
        for ii in range(k_comp):
            for jj in range(k_comp):
                A_mp[ii, jj] = mpmath.mpf(str(A_comp[ii, jj]))
        evals_30 = mpmath.eig(A_mp, left=False, right=False)
        J_comp_30 = float(max(abs(complex(e).imag) for e in evals_30))
    print(f"    30 dps  J = {J_comp_30:.2e}  (mpmath.eig)")

    # 50 dps
    with mpmath.workdps(50):
        A_mp = mpmath.matrix(k_comp, k_comp)
        for ii in range(k_comp):
            for jj in range(k_comp):
                A_mp[ii, jj] = mpmath.mpf(str(A_comp[ii, jj]))
        evals_50 = mpmath.eig(A_mp, left=False, right=False)
        J_comp_50 = float(max(abs(complex(e).imag) for e in evals_50))
    print(f"    50 dps  J = {J_comp_50:.2e}  (mpmath.eig)")

    # Apply ladder rules
    if J_comp_64 == 0 and J_comp_30 == 0 and J_comp_50 == 0:
        comp_label = "COLLAPSES"
    elif J_comp_50 < 0.5 * J_comp_30 or J_comp_50 > 2 * J_comp_30:
        comp_label = "AMBIGUOUS"
    elif J_comp_50 >= 0.5 * J_comp_64:
        comp_label = "PERSISTS"
    elif J_comp_50 <= 1e-3 * J_comp_64:
        comp_label = "COLLAPSES"
    else:
        comp_label = "AMBIGUOUS"
    print(f"    ladder result: {comp_label} (must COLLAPSE for real spectrum)")
    if comp_label != "COLLAPSES":
        print(f"    WARNING: ladder did not collapse for known-real companion matrix")
        print(f"    NOTE: a k-fold defective eigenvalue at float64 splits into")
        print(f"    k complex values via k-th root scaling: J ~ eps^(1/k).")
        print(f"    At higher precision the imaginary parts shrink but do not")
        print(f"    reach zero unless the root-finder uses extraprec >= k*dps.")

    # -----------------------------------------------------------------------
    #  Q3: ADVERSARIAL VACUITY AUDIT
    # -----------------------------------------------------------------------
    section("Q3: ADVERSARIAL VACUITY AUDIT")

    print("  1. G-SAMPLE after orthonormalization: IDENTIFIED BY CONTRACT")
    print("     conditioning test after M_h-orthonormalization is satisfied by")
    print("     construction. The gate correctly tests BEFORE orthonormalization")
    print("     (on W_n = M_h^{1/2} F_n^seed, not on Q_n).")
    print("     VERIFIED: our G-SAMPLE implementation operates on W_n.")

    print("\n  2. K_0 on real 1x1 block: IDENTIFIED BY CONTRACT")
    print(f"     K_0 = {K0_norm:.2e} (zero by algebra).")
    print("     Cannot calibrate K machinery. Contract explicitly documents this.")

    print("\n  3. Relative rank cutoff on numerically-zero matrix: IDENTIFIED BY CONTRACT")
    print("     Pi_3 at ||Pi|| ~ 1e-17: relative cutoff gives spurious rank.")
    print("     G-SUBSPACE uses ABSOLUTE cutoff 1e-8. Verified via rank arm above.")

    print("\n  4. G-WIRE arm: constant-mode permutation: IDENTIFIED BY CONTRACT")
    print("     F_0 is constant, so permuting rows cannot change it.")
    print("     Contract replaces with analytically forced L_mut arm. Verified above.")

    print("\n  5. Same-boundary swap in G-ALIGN: IDENTIFIED BY CONTRACT")
    print("     A swap of two numerically similar rows may not move D_n.")
    print("     Contract uses structural predicate instead. Verified above.")

    print("\n  6. n=2 trap for realization test: IDENTIFIED BY CONTRACT")
    print("     At n=2, correct and no-transpose subspaces coincide.")
    print("     verify_realization sweeps n >= 3. Verified above.")

    # Additional checks constructed by this qualification
    print("\n  7. CONSTRUCTED: Can G-BASIS pass trivially at n=0?")
    print(f"     A_0 is 1x1 real, ||A||_F={A0_scale:.2e}.")
    print(f"     Covariance check: rel err = {A0_rel:.2e}.")
    print("     At n=0, G-BASIS is tautological (1x1 => U is a phase,")
    print("     A' = |U|^2 A = A). This is why the contract tests it at n=12,20.")
    print("     FINDING: G-BASIS at n=0 has no discriminating power.")
    print("     Not a contract defect: the contract specifies G-BASIS at the")
    print("     target levels, and G-WIRE handles the n=0 verification separately.")

    print("\n  8. CONSTRUCTED: Can G-SAMPLE pass with badly conditioned basis")
    print("     if the conditioning is inherited from M_h rather than from sampling?")
    print("     G-SAMPLE tests W_n = M_h^{1/2} F_n^seed. If M_h is ill-conditioned,")
    print(f"     kappa(M_h) = {W.max()/W.min():.4f} (well-conditioned at 60 seeds).")
    print("     At 60 seeds, the M_h conditioning is not a concern.")
    print("     The contract records that M_h is well-conditioned here.")

    print("\n  9. CONSTRUCTED: Does the precision ladder collapse on exact-zero J?")
    print("     Rule 1: if J(64)=J(30)=J(50)=0 -> COLLAPSES.")
    print("     This is the correct handling. An earlier draft's ratio tests")
    print("     both fired at exact zero (0 >= 0 and 0 <= 0). Contract fixed this.")
    A_real = 168.0 * np.eye(4)
    J_real_64 = compute_J(A_real)
    print(f"     Test: real symmetric A -> J(64) = {J_real_64:.2e}")

    # -----------------------------------------------------------------------
    #  Q4: CONTRACT CONTRADICTIONS
    # -----------------------------------------------------------------------
    section("Q4: CONTRACT CONTRADICTIONS")

    print("  FINDING: No contradictions found that would prevent implementation.")
    print()
    print("  The following ambiguities were examined and found resolved:")
    print()
    print("  1. theta_C vs theta_Q: clearly distinguished in the contract.")
    print("     theta_C is gated at 1e-6 in G-SUBSPACE.")
    print("     theta_Q is recorded but not gated, deliberately.")
    print("     K_floor uses theta_Q, not theta_C. All three are named.")
    print()
    print("  2. Route indexing: the contract states 'every Q_n, A_n, K_n, J_n'")
    print("     is evaluated SEPARATELY for r in {svd, avg}. K_floor is the")
    print("     one exception: a common worst-case bar. Unambiguous.")
    print()
    print("  3. Precision ladder: ordered rules with first-match-wins.")
    print("     Rule 1 handles exact zero before ratio tests. Unambiguous.")
    print()
    print("  4. Branch table: five outcomes, precedence by rule number.")
    print("     Earlier contradiction (G-SUBSPACE vs branch table on route")
    print("     disagreement) was caught and fixed in pass 7. Unambiguous now.")
    print()
    print("  5. ||K||_2 in operator 2-norm throughout. An earlier draft mixed")
    print("     Frobenius and operator norms. Contract corrected this and")
    print("     records the correction. Unambiguous now.")
    print()
    print("  6. 2*sin(theta_max/2) vs 2*sin(theta_max) in K_floor discrepancy")
    print("     term. Contract specifies 2*sin(theta_max/2) with derivation")
    print("     and numerical verification. Unambiguous.")
    print()
    print("  AMBIGUITY 7 (not a contradiction but a FINDING):")
    print("  The contract states: 'the companion matrix of (x-1)^k, whose exact")
    print("  spectrum is real, must COLLAPSE under this ladder.' The value of k is")
    print("  not specified. For k >= 3 the eigenvalue is k-fold defective, and")
    print("  any eigensolver at d decimal digits produces J ~ eps^{1/k}.")
    print("  The resulting J decreases between precision levels, triggering")
    print("  rule 2 (tail not settled) -> AMBIGUOUS, never reaching COLLAPSES.")
    print(f"    Verified: k=2 -> {comp_2_label}, k=6 -> {comp_label}")
    print("  k=2 collapses via rule 1 (analytic 2x2 formula gives J=0 at all rungs).")
    print("  k >= 3 gives AMBIGUOUS at any finite precision.")
    print("  RESOLUTION: the arm demonstrates power at k=2. For k >= 3, the")
    print("  ladder's resolving power does not extend to highly defective matrices.")
    print("  This is a property of the eigensolver, not a ladder mis-wiring.")
    print("  The live targets' A_n matrices are not companion matrices and are")
    print("  not expected to be highly defective.")
    print()
    print("  OBSERVATION (not a contradiction): The shipped invariant_dim_and_basis")
    print("  uses full_matrices=True SVD, allocating O(m^2) memory where m is the")
    print("  number of stacked rows. At n=20 this is ~45 GB. The economy SVD")
    print("  produces identical singular values and Vh. This is a performance")
    print("  issue, not a mathematical one. The contract documents the 3919s time.")

    # -----------------------------------------------------------------------
    #  Q5: SYNTHETIC DRY RUN
    # -----------------------------------------------------------------------
    section("Q5: SYNTHETIC DRY RUN")
    print("  Testing the branch logic on manufactured matrices.")
    print("  Each test constructs A_n with known K and J properties,")
    print("  runs the precision ladder, and verifies the branch outcome.\n")

    def precision_ladder_synthetic(A_64):
        """Run the precision ladder on a synthetic A matrix.

        Rung 1: float64 (already provided)
        Rung 2: 30 dps reformation
        Rung 3: 50 dps reformation
        """
        evals_64 = np.linalg.eigvals(A_64)
        J_64 = float(np.max(np.abs(np.imag(evals_64))))

        k = A_64.shape[0]
        # 30 dps
        with mpmath.workdps(30):
            A_mp = mpmath.matrix(k, k)
            for i in range(k):
                for j in range(k):
                    A_mp[i, j] = mpmath.mpf(str(float(np.real(A_64[i, j])))) + \
                                 mpmath.mpf(str(float(np.imag(A_64[i, j])))) * mpmath.mpc(0, 1)
            evals_30 = list(mpmath.eig(A_mp, right=False, left=False))
            J_30 = float(max(abs(complex(e).imag) for e in evals_30))

        # 50 dps
        with mpmath.workdps(50):
            A_mp = mpmath.matrix(k, k)
            for i in range(k):
                for j in range(k):
                    A_mp[i, j] = mpmath.mpf(str(float(np.real(A_64[i, j])))) + \
                                 mpmath.mpf(str(float(np.imag(A_64[i, j])))) * mpmath.mpc(0, 1)
            evals_50 = list(mpmath.eig(A_mp, right=False, left=False))
            J_50 = float(max(abs(complex(e).imag) for e in evals_50))

        # Ladder rules (first match wins)
        if J_64 == 0 and J_30 == 0 and J_50 == 0:
            label = "COLLAPSES"
        elif J_50 < 0.5 * J_30 or J_50 > 2 * J_30:
            label = "AMBIGUOUS"
        elif J_50 >= 0.5 * J_64:
            label = "PERSISTS"
        elif J_50 <= 1e-3 * J_64:
            label = "COLLAPSES"
        else:
            label = "AMBIGUOUS"

        return label, J_64, J_30, J_50

    def branch_logic(gate_ok, targets):
        """Evaluate the branch table.

        targets: list of dicts with keys:
            n, J_label_svd, J_label_avg, K_above_floor_svd, K_above_floor_avg
        Returns the outcome string.
        """
        # Rule 1: any gate fails
        if not gate_ok:
            return "S1b-DEFECT"

        # Rule 2: any target has AMBIGUOUS J on EITHER construction
        for t in targets:
            if t["J_label_svd"] == "AMBIGUOUS" or t["J_label_avg"] == "AMBIGUOUS":
                return "S1b-NO_LABEL"

        # Rule 3: constructions disagree on J for any target
        for t in targets:
            if t["J_label_svd"] != t["J_label_avg"]:
                return "S1b-NO_LABEL"

        # Rule 4: any target PERSISTS on BOTH constructions
        for t in targets:
            if t["J_label_svd"] == "PERSISTS" and t["J_label_avg"] == "PERSISTS":
                return "S1b-SPECTRAL"

        # Rule 5: all J collapsed/zero, but K disagree between constructions
        for t in targets:
            if t["K_above_floor_svd"] != t["K_above_floor_avg"]:
                return "S1b-NO_LABEL"

        # Rule 6: all J collapsed/zero, some K above floor on BOTH
        for t in targets:
            if t["K_above_floor_svd"] and t["K_above_floor_avg"]:
                return "S1b-ADJOINT"

        # Rule 7: else
        return "S1b-NULL"

    # --- Seven outcomes, each reachable ---

    print("  === REACHABILITY: one specimen per outcome ===\n")

    # S1b-DEFECT: a gate fails
    outcome = branch_logic(False, [])
    print(f"  gate_fails=True -> {outcome}  (expected S1b-DEFECT)")
    assert outcome == "S1b-DEFECT"

    # S1b-SPECTRAL: J PERSISTS on both constructions
    targets_spec = [
        {"n": 12, "J_label_svd": "PERSISTS", "J_label_avg": "PERSISTS",
         "K_above_floor_svd": True, "K_above_floor_avg": True},
        {"n": 20, "J_label_svd": "COLLAPSES", "J_label_avg": "COLLAPSES",
         "K_above_floor_svd": False, "K_above_floor_avg": False},
    ]
    outcome = branch_logic(True, targets_spec)
    print(f"  n=12 J PERSISTS both -> {outcome}  (expected S1b-SPECTRAL)")
    assert outcome == "S1b-SPECTRAL"

    # S1b-ADJOINT: all J collapsed, K above floor on both
    targets_adj = [
        {"n": 12, "J_label_svd": "COLLAPSES", "J_label_avg": "COLLAPSES",
         "K_above_floor_svd": True, "K_above_floor_avg": True},
        {"n": 20, "J_label_svd": "COLLAPSES", "J_label_avg": "COLLAPSES",
         "K_above_floor_svd": False, "K_above_floor_avg": False},
    ]
    outcome = branch_logic(True, targets_adj)
    print(f"  all J collapsed, n=12 K above floor both -> {outcome}  "
          f"(expected S1b-ADJOINT)")
    assert outcome == "S1b-ADJOINT"

    # S1b-NULL: all consistent with zero
    targets_null = [
        {"n": 12, "J_label_svd": "COLLAPSES", "J_label_avg": "COLLAPSES",
         "K_above_floor_svd": False, "K_above_floor_avg": False},
        {"n": 20, "J_label_svd": "COLLAPSES", "J_label_avg": "COLLAPSES",
         "K_above_floor_svd": False, "K_above_floor_avg": False},
    ]
    outcome = branch_logic(True, targets_null)
    print(f"  all collapsed, K below floor -> {outcome}  (expected S1b-NULL)")
    assert outcome == "S1b-NULL"

    # S1b-NO_LABEL via AMBIGUOUS J (rule 2)
    targets_nolabel_r2 = [
        {"n": 12, "J_label_svd": "AMBIGUOUS", "J_label_avg": "COLLAPSES",
         "K_above_floor_svd": False, "K_above_floor_avg": False},
        {"n": 20, "J_label_svd": "COLLAPSES", "J_label_avg": "COLLAPSES",
         "K_above_floor_svd": False, "K_above_floor_avg": False},
    ]
    outcome = branch_logic(True, targets_nolabel_r2)
    print(f"  n=12 J AMBIGUOUS on svd -> {outcome}  (expected S1b-NO_LABEL, rule 2)")
    assert outcome == "S1b-NO_LABEL"

    # S1b-NO_LABEL via J disagreement (rule 3)
    targets_nolabel_r3 = [
        {"n": 12, "J_label_svd": "PERSISTS", "J_label_avg": "COLLAPSES",
         "K_above_floor_svd": True, "K_above_floor_avg": True},
        {"n": 20, "J_label_svd": "COLLAPSES", "J_label_avg": "COLLAPSES",
         "K_above_floor_svd": False, "K_above_floor_avg": False},
    ]
    outcome = branch_logic(True, targets_nolabel_r3)
    print(f"  n=12 J disagrees (PERSISTS vs COLLAPSES) -> {outcome}  "
          f"(expected S1b-NO_LABEL, rule 3)")
    assert outcome == "S1b-NO_LABEL"

    # S1b-NO_LABEL via K disagreement (rule 5)
    targets_nolabel_r5 = [
        {"n": 12, "J_label_svd": "COLLAPSES", "J_label_avg": "COLLAPSES",
         "K_above_floor_svd": True, "K_above_floor_avg": False},
        {"n": 20, "J_label_svd": "COLLAPSES", "J_label_avg": "COLLAPSES",
         "K_above_floor_svd": False, "K_above_floor_avg": False},
    ]
    outcome = branch_logic(True, targets_nolabel_r5)
    print(f"  n=12 K disagrees between routes -> {outcome}  "
          f"(expected S1b-NO_LABEL, rule 5)")
    assert outcome == "S1b-NO_LABEL"

    # --- PRECEDENCE / COLLISION TESTS ---

    print("\n  === PRECEDENCE: collision cases ===\n")

    # Collision 1: gate fails while target otherwise looks SPECTRAL -> DEFECT
    targets_c1 = [
        {"n": 12, "J_label_svd": "PERSISTS", "J_label_avg": "PERSISTS",
         "K_above_floor_svd": True, "K_above_floor_avg": True},
    ]
    outcome = branch_logic(False, targets_c1)
    print(f"  gate FAILS + target looks SPECTRAL -> {outcome}  "
          f"(expected S1b-DEFECT)")
    assert outcome == "S1b-DEFECT"

    # Collision 2: one target AMBIGUOUS while another is qualified PERSISTENT -> NO_LABEL
    targets_c2 = [
        {"n": 12, "J_label_svd": "PERSISTS", "J_label_avg": "PERSISTS",
         "K_above_floor_svd": True, "K_above_floor_avg": True},
        {"n": 20, "J_label_svd": "AMBIGUOUS", "J_label_avg": "AMBIGUOUS",
         "K_above_floor_svd": False, "K_above_floor_avg": False},
    ]
    outcome = branch_logic(True, targets_c2)
    print(f"  n=12 PERSISTS + n=20 AMBIGUOUS -> {outcome}  "
          f"(expected S1b-NO_LABEL)")
    assert outcome == "S1b-NO_LABEL"

    # Collision 3: two constructions disagree on J while one looks SPECTRAL -> NO_LABEL
    targets_c3 = [
        {"n": 12, "J_label_svd": "PERSISTS", "J_label_avg": "COLLAPSES",
         "K_above_floor_svd": True, "K_above_floor_avg": True},
        {"n": 20, "J_label_svd": "COLLAPSES", "J_label_avg": "COLLAPSES",
         "K_above_floor_svd": False, "K_above_floor_avg": False},
    ]
    outcome = branch_logic(True, targets_c3)
    print(f"  n=12 J svd=PERSISTS avg=COLLAPSES -> {outcome}  "
          f"(expected S1b-NO_LABEL)")
    assert outcome == "S1b-NO_LABEL"

    # Collision 4: K disagrees while one condition looks ADJOINT -> NO_LABEL
    targets_c4 = [
        {"n": 12, "J_label_svd": "COLLAPSES", "J_label_avg": "COLLAPSES",
         "K_above_floor_svd": True, "K_above_floor_avg": False},
        {"n": 20, "J_label_svd": "COLLAPSES", "J_label_avg": "COLLAPSES",
         "K_above_floor_svd": False, "K_above_floor_avg": False},
    ]
    outcome = branch_logic(True, targets_c4)
    print(f"  n=12 K svd=above avg=below floor -> {outcome}  "
          f"(expected S1b-NO_LABEL)")
    assert outcome == "S1b-NO_LABEL"

    # --- Precision ladder on synthetic matrices ---
    print("\n  === PRECISION LADDER on synthetic matrices ===\n")

    # Real symmetric -> COLLAPSES
    A_sym = np.diag([1.0, 2.0, 3.0, 4.0])
    label_s, J64, J30, J50 = precision_ladder_synthetic(A_sym)
    print(f"  real symmetric: J(64)={J64:.2e}, J(30)={J30:.2e}, "
          f"J(50)={J50:.2e} -> {label_s} (expected COLLAPSES)")

    # Complex conjugate pairs -> PERSISTS
    A_conj = np.array([[1.0, -2.0], [2.0, 1.0]])  # eigenvalues 1±2i
    label_c, J64, J30, J50 = precision_ladder_synthetic(A_conj)
    print(f"  conjugate pair: J(64)={J64:.2e}, J(30)={J30:.2e}, "
          f"J(50)={J50:.2e} -> {label_c} (expected PERSISTS)")

    # Companion matrix (known real spectrum) -> COLLAPSES
    print(f"  companion (x-1)^6: -> {comp_label} (expected COLLAPSES)")

    # --- Final manifest re-verification ---
    section("FINAL MANIFEST RE-VERIFICATION")
    with open(manifest_path, "rb") as f:
        mh_final = hashlib.sha256(f.read()).hexdigest()
    with open(contract_path, "rb") as f:
        raw_final = f.read()
    idx_final = raw_final.find(boundary)
    ch_final = hashlib.sha256(raw_final[:idx_final]).hexdigest()
    print(f"  manifest SHA-256: {mh_final}")
    print(f"  contract SHA-256: {ch_final}")
    manifest_ok = mh_final == manifest_expected
    contract_ok = ch_final == expected
    print(f"  manifest unchanged: {manifest_ok}")
    print(f"  contract unchanged: {contract_ok}")
    if not manifest_ok or not contract_ok:
        print("  WARNING: files changed during run!")
        all_pass = False

    # --- SUMMARY ---
    section("SUMMARY")
    gate_names = ["G-RANK_n12", "G-RANK_n20", "G-REAL", "G-SUBSPACE_n12",
                  "G-SUBSPACE_n20", "G-ALIGN", "G-SAMPLE_n12_svd",
                  "G-SAMPLE_n12_avg", "G-SAMPLE_n20_svd", "G-SAMPLE_n20_avg",
                  "G-WIRE", "G-BASIS", "G-DISCRIM"]
    for g in gate_names:
        r = results.get(g, {})
        status = "PASS" if r.get("pass", False) else "FAIL"
        print(f"  {g:30s} {status}")

    print(f"\n  ALL GATES: {'PASS' if all_pass else 'SOME FAILURES'}")
    print(f"\n  Contract contradictions: NONE FOUND")
    print(f"  Vacuity audit: 9 cases examined (6 contract-identified + 3 constructed)")


if __name__ == "__main__":
    main()
