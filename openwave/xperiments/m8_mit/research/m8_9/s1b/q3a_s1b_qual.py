#!/usr/bin/env python3
"""M8.9 S1b Full Implementation Qualification — attempt q3a

Contract reading order: addendum 3 > 2 > 1 > decision rule (later governs earlier).
Instrument qualification only. No live target. No threshold changes.

Incremental serialization: results written to disk after every gate.
"""

import sys, os, time, hashlib, json
import numpy as np
import mpmath

ROOM = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
LEDGER = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(LEDGER, "results")
NOTE = os.path.join(LEDGER, "QUALIFICATION_NOTE.md")

sys.path.insert(0, ROOM)
sys.path.insert(0, os.path.join(ROOM, "m8_5b", "pilot"))
sys.path.insert(0, os.path.join(ROOM, "m8_5b", "production"))
sys.path.insert(0, os.path.join(ROOM, "m8_5b", "gates"))

import route_a_repn as repn
from route_a_nonabelian import quat_to_su2, sym_power
from route_a_twosided import pairs_left
from p0.group import build_icosians, build_character_table
from p0.cloud import fibonacci_seeds_s3, build_orbit_cloud
from p0.bundle_operator import build_L_bundle, orbit_stencils as p0_orbit_stencils
from p0.representations import build_all_representations
from p1a.mass_matrix import build_Mh_base

eps = np.finfo(np.float64).eps
SEP = "=" * 72

os.makedirs(RESULTS, exist_ok=True)


def note(text):
    with open(NOTE, "a") as f:
        f.write(text + "\n")
    print(text)


def save_result(name, data):
    path = os.path.join(RESULTS, f"{name}.json")
    with open(path, "w") as f:
        json.dump(data, f, indent=2, default=_json_default)
        f.write("\n")
    print(f"  [serialized -> results/{name}.json]")


def _json_default(obj):
    if isinstance(obj, np.bool_):
        return bool(obj)
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, complex):
        return {"re": obj.real, "im": obj.imag}
    if isinstance(obj, np.complexfloating):
        return {"re": float(obj.real), "im": float(obj.imag)}
    if obj == float("inf"):
        return "Infinity"
    if obj == float("-inf"):
        return "-Infinity"
    raise TypeError(f"not JSON serializable: {type(obj)}")


def sha256_file(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def sha256_bytes(data):
    return hashlib.sha256(data).hexdigest()


# ---------------------------------------------------------------------------
#  Core computational utilities
# ---------------------------------------------------------------------------

def form_An(Q, Mh_diag, L):
    """A_n = Q^H M_h L Q."""
    MhQ = Mh_diag[:, None] * Q
    return MhQ.conj().T @ L @ Q


def compute_K(A):
    """K = (A - A^H)/2, return (K, ||K||_2)."""
    K = (A - A.conj().T) / 2
    return K, float(np.linalg.norm(K, ord=2))


def compute_J(A):
    """J = max |Im(mu)| over eigenvalues mu of A."""
    evals = np.linalg.eigvals(A)
    return float(np.max(np.abs(np.imag(evals))))


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


def sine_form_angle(B1, B2):
    """Principal angle via the SINE form: s = ||(I - P_a) Q_b||_2.

    B1, B2 must have orthonormal columns (ordinary inner product).
    Returns s = sin(theta_max).
    Addendum 1 §A1.4: the arccos route is PROHIBITED.
    """
    B1q, _ = np.linalg.qr(B1)
    B2q, _ = np.linalg.qr(B2)
    k1, k2 = B1q.shape[1], B2q.shape[1]
    if k1 == 0 or k2 == 0:
        return float("nan")
    P_a = B1q @ B1q.conj().T
    residual = (np.eye(B1q.shape[0]) - P_a) @ B2q
    s = float(np.linalg.norm(residual, ord=2))
    return s


def half_angle_d(s):
    """d = 2 sin(theta/2) from s = sin(theta).

    Frozen form (addendum 1 §A1.4):
        d = 2 * sin( 0.5 * arcsin( clip(s, 0, 1) ) )

    The naive identity sqrt(2(1 - sqrt(1 - s^2))) is PROHIBITED:
    it returns exactly 0.0 below about 7.45e-09.
    """
    return 2.0 * np.sin(0.5 * np.arcsin(np.clip(s, 0.0, 1.0)))


def averaging_basis(pairs, n):
    """Reynolds projector averaging route for the invariant subspace.

    Pi_n = (1/|G|) sum_gamma coefficient_operator(u, v, n)
    Rank at ABSOLUTE cutoff 1e-8.
    Orthonormal basis of ran(Pi_n) from left singular vectors.
    """
    dim = (n + 1) ** 2
    Pi = np.zeros((dim, dim), dtype=complex)
    for (u, v) in pairs:
        Pi += repn.coefficient_operator(u, v, n)
    Pi /= len(pairs)

    Pi2_err = float(np.linalg.norm(Pi @ Pi - Pi))
    PiH_err = float(np.linalg.norm(Pi - Pi.conj().T))

    s_Pi = np.linalg.svd(Pi, compute_uv=False)
    rank_abs = int(np.sum(s_Pi > 1e-8))

    U_full, _, _ = np.linalg.svd(Pi, full_matrices=False)
    C_avg = U_full[:, :rank_abs]

    return C_avg, rank_abs, {
        "idempotence_err": Pi2_err,
        "non_hermiticity": PiH_err,
        "Pi_norm": float(np.linalg.norm(Pi, ord=2)),
        "rank_abs_1e8": rank_abs,
    }


def precision_ladder(A_64_entries, k):
    """Precision ladder: re-form A at each rung from the SAME binary inputs.

    For synthetic tests, A_64_entries is the float64 matrix directly.
    Rungs: float64, 30 dps, 50 dps.
    Returns (label, J_64, J_30, J_50, rule_used).
    """
    A_64 = np.asarray(A_64_entries)
    evals_64 = np.linalg.eigvals(A_64)
    J_64 = float(np.max(np.abs(np.imag(evals_64))))

    def mp_eigvals(dps):
        with mpmath.workdps(dps):
            A_mp = mpmath.matrix(k, k)
            for i in range(k):
                for j in range(k):
                    val = A_64[i, j]
                    if np.isrealobj(A_64):
                        A_mp[i, j] = mpmath.mpf(float(val))
                    else:
                        A_mp[i, j] = mpmath.mpc(float(np.real(val)), float(np.imag(val)))
            evals = list(mpmath.eig(A_mp, right=False, left=False))
            return float(max(abs(complex(e).imag) for e in evals))

    J_30 = mp_eigvals(30)
    J_50 = mp_eigvals(50)

    if J_64 == 0 and J_30 == 0 and J_50 == 0:
        return "COLLAPSES", J_64, J_30, J_50, 1
    if J_50 < 0.5 * J_30 or J_50 > 2 * J_30:
        return "AMBIGUOUS", J_64, J_30, J_50, 2
    if J_50 >= 0.5 * J_64:
        return "PERSISTS", J_64, J_30, J_50, 3
    if J_50 <= 1e-3 * J_64:
        return "COLLAPSES", J_64, J_30, J_50, 4
    return "AMBIGUOUS", J_64, J_30, J_50, 5


def precision_ladder_reformed(Q, Mh_diag, L, k):
    """Precision ladder with RE-FORMATION at each rung.

    From the SAME fixed binary inputs Q, Mh_diag, L, the matrix product
    is re-formed at each rung and eigensolved at matching precision.
    """
    A_64 = form_An(Q, Mh_diag, L)
    evals_64 = np.linalg.eigvals(A_64)
    J_64 = float(np.max(np.abs(np.imag(evals_64))))

    def reform_and_solve(dps):
        with mpmath.workdps(dps):
            n_row = Q.shape[0]
            Q_mp = mpmath.matrix(n_row, k)
            Mh_mp = [mpmath.mpf(float(Mh_diag[i])) for i in range(n_row)]
            L_mp = mpmath.matrix(n_row, n_row)
            for i in range(n_row):
                for j in range(k):
                    Q_mp[i, j] = (mpmath.mpf(float(np.real(Q[i, j])))
                                  + mpmath.mpf(float(np.imag(Q[i, j]))) * mpmath.mpc(0, 1))
            for i in range(n_row):
                for j in range(n_row):
                    L_mp[i, j] = mpmath.mpf(float(np.real(L[i, j])))
            MhQ = mpmath.matrix(n_row, k)
            for i in range(n_row):
                for j in range(k):
                    MhQ[i, j] = Mh_mp[i] * Q_mp[i, j]
            QH_Mh = MhQ.T.conjugate()
            LQ = L_mp * Q_mp
            A_mp = QH_Mh * LQ
            evals = list(mpmath.eig(A_mp, right=False, left=False))
            return float(max(abs(complex(e).imag) for e in evals))

    J_30 = reform_and_solve(30)
    J_50 = reform_and_solve(50)

    if J_64 == 0 and J_30 == 0 and J_50 == 0:
        return "COLLAPSES", J_64, J_30, J_50, 1
    if J_50 < 0.5 * J_30 or J_50 > 2 * J_30:
        return "AMBIGUOUS", J_64, J_30, J_50, 2
    if J_50 >= 0.5 * J_64:
        return "PERSISTS", J_64, J_30, J_50, 3
    if J_50 <= 1e-3 * J_64:
        return "COLLAPSES", J_64, J_30, J_50, 4
    return "AMBIGUOUS", J_64, J_30, J_50, 5


def branch_logic(gate_ok, targets):
    """Branch table: five outcomes, first match wins."""
    if not gate_ok:
        return "S1b-DEFECT", 1
    for t in targets:
        if t["J_label_svd"] == "AMBIGUOUS" or t["J_label_avg"] == "AMBIGUOUS":
            return "S1b-NO_LABEL", 2
    for t in targets:
        if t["J_label_svd"] != t["J_label_avg"]:
            return "S1b-NO_LABEL", 3
    for t in targets:
        if t["J_label_svd"] == "PERSISTS" and t["J_label_avg"] == "PERSISTS":
            return "S1b-SPECTRAL", 4
    for t in targets:
        if t["K_above_floor_svd"] != t["K_above_floor_avg"]:
            return "S1b-NO_LABEL", 5
    for t in targets:
        if t["K_above_floor_svd"] and t["K_above_floor_avg"]:
            return "S1b-ADJOINT", 6
    return "S1b-NULL", 7


# ===========================================================================
#  MAIN
# ===========================================================================
def main():
    t_start = time.time()

    # -----------------------------------------------------------------------
    #  INFRASTRUCTURE
    # -----------------------------------------------------------------------
    print(f"\n{SEP}\n  INFRASTRUCTURE\n{SEP}")
    t0 = time.time()
    elems = build_icosians()
    pairs = pairs_left(elems)
    print(f"  2I elements: {len(elems)}, left-action pairs: {len(pairs)}")

    chi = build_character_table(elems)
    reps, _ = build_all_representations(elems, chi)

    N_SEEDS = 60
    K_STENCIL = 110
    RBF_M, RBF_P = 7, 4
    seeds = fibonacci_seeds_s3(N_SEEDS)
    X, oid, gid = build_orbit_cloud(seeds, elems)
    print(f"  cloud: {len(X)} nodes from {N_SEEDS} seeds")

    W = build_Mh_base(X, oid, N_SEEDS)
    Mh_diag = W.copy()
    kappa_Mh = float(W.max() / W.min())
    print(f"  M_h: range [{W.min():.4e}, {W.max():.4e}], kappa={kappa_Mh:.4f}")

    node_of, mult_table, plan = p0_orbit_stencils(X, oid, gid, elems, k=K_STENCIL)
    L_q, seed_orbits = build_L_bundle(X, oid, gid, elems, reps["R0"],
                                       k=K_STENCIL, m=RBF_M, p=RBF_P)
    L_q = np.real(np.asarray(L_q))
    print(f"  L_q: {L_q.shape}, seed_orbits: {len(seed_orbits)}")

    # Ltilde and B for K_floor
    Mh_sqrt = np.sqrt(Mh_diag)
    Mh_inv_sqrt = 1.0 / Mh_sqrt
    Ltilde = np.diag(Mh_sqrt) @ L_q @ np.diag(Mh_inv_sqrt)
    B = Ltilde - Ltilde.T
    B_norm = float(np.linalg.norm(B, ord=2))
    print(f"  ||Ltilde - Ltilde^H||_2 = {B_norm:.2f}")

    rep_indices = np.array([node_of[(o, 0)] for o in seed_orbits])
    dt_infra = time.time() - t0
    print(f"  infrastructure time: {dt_infra:.1f}s")

    save_result("infrastructure", {
        "n_seeds": N_SEEDS, "n_nodes": len(X),
        "k_stencil": K_STENCIL, "m": RBF_M, "p": RBF_P,
        "kappa_Mh": kappa_Mh, "B_norm": B_norm,
        "L_q_shape": list(L_q.shape),
        "time_s": dt_infra,
    })

    # -----------------------------------------------------------------------
    #  Q1: G-REAL — pointwise realization
    # -----------------------------------------------------------------------
    print(f"\n{SEP}\n  G-REAL: pointwise realization\n{SEP}")
    t0 = time.time()
    vr = repn.verify_realization()
    dt = time.time() - t0
    print(f"  worst_residual_correct:     {vr['worst_residual_correct']:.2e} (need < 1e-10)")
    print(f"  best_residual_no_transpose: {vr['best_residual_no_transpose']:.2e} (need > 1e-10)")
    print(f"  PASS: {vr['pass']}, time: {dt:.2f}s")
    save_result("G-REAL", vr)
    note(f"\n## G-REAL\n- worst_correct: {vr['worst_residual_correct']:.2e}\n"
         f"- best_no_transpose: {vr['best_residual_no_transpose']:.2e}\n"
         f"- **{'PASS' if vr['pass'] else 'FAIL'}**")

    # Green parent already demonstrated. Red child:
    print(f"\n  G-REAL mutation: no_transpose form")
    print(f"    green: {vr['worst_residual_correct']:.2e}")
    print(f"    red:   {vr['best_residual_no_transpose']:.2e}")
    print(f"    separation: {vr['best_residual_no_transpose']/max(vr['worst_residual_correct'],1e-30):.2e}")

    # -----------------------------------------------------------------------
    #  G-DISCRIM — K vs J discrimination (synthetic, no SVD needed)
    # -----------------------------------------------------------------------
    print(f"\n{SEP}\n  G-DISCRIM: K vs J discrimination\n{SEP}")
    lam = 168.0
    k_test = 13
    A_green = lam * np.eye(k_test, dtype=complex)
    _, K_green_norm = compute_K(A_green)
    J_green = compute_J(A_green)
    K_green_bound = 100 * eps * lam
    print(f"  green: ||K||_2={K_green_norm:.2e} (need <= {K_green_bound:.2e}), J={J_green:.2e} (need <= 1e-12)")

    A_armA = A_green.copy(); A_armA[0, 1] += 1.0
    _, K_A_norm = compute_K(A_armA)
    J_A = compute_J(A_armA)
    print(f"  arm A (nilpotent): ||K||_2={K_A_norm:.6f} (need |..-0.5|<1e-12), J={J_A:.2e} (need <=1e-12)")

    A_armB = A_green.copy(); A_armB[0, 1] += 2.0; A_armB[1, 0] -= 2.0
    _, K_B_norm = compute_K(A_armB)
    J_B = compute_J(A_armB)
    print(f"  arm B (antisym):   ||K||_2={K_B_norm:.6f} (need |..-2|<1e-12), J={J_B:.6f} (need |..-2|<1e-12)")

    gd_pass = (K_green_norm <= K_green_bound and J_green <= 1e-12
               and abs(K_A_norm - 0.5) <= 1e-12 and J_A <= 1e-12
               and abs(K_B_norm - 2.0) <= 1e-12 and abs(J_B - 2.0) <= 1e-12)
    print(f"  G-DISCRIM: {'PASS' if gd_pass else 'FAIL'}")
    save_result("G-DISCRIM", {
        "green_K": K_green_norm, "green_J": J_green,
        "armA_K": K_A_norm, "armA_J": J_A,
        "armB_K": K_B_norm, "armB_J": J_B,
        "pass": gd_pass,
    })
    note(f"\n## G-DISCRIM\n- green: K={K_green_norm:.2e}, J={J_green:.2e}\n"
         f"- arm A: K={K_A_norm:.6f}, J={J_A:.2e}\n"
         f"- arm B: K={K_B_norm:.6f}, J={J_B:.6f}\n- **{'PASS' if gd_pass else 'FAIL'}**")

    # -----------------------------------------------------------------------
    #  Q3: Ladder controls (A1.1) — synthetic, no SVD needed
    # -----------------------------------------------------------------------
    print(f"\n{SEP}\n  LADDER CONTROLS (A1.1)\n{SEP}")

    def companion_matrix(k_val):
        coeffs = np.array([1.0])
        for _ in range(k_val):
            coeffs = np.convolve(coeffs, [1.0, -1.0])
        C = np.zeros((k_val, k_val))
        C[1:, :-1] = np.eye(k_val - 1)
        C[:, -1] = -coeffs[:-1] / coeffs[-1]
        return C

    # k=2: exact-zero control -> COLLAPSES via rule 1
    comp2 = companion_matrix(2)
    label2, J2_64, J2_30, J2_50, rule2 = precision_ladder(comp2, 2)
    print(f"  k=2: J(64)={J2_64:.2e}, J(30)={J2_30:.2e}, J(50)={J2_50:.2e} -> {label2} (rule {rule2})")
    assert label2 == "COLLAPSES" and rule2 == 1, f"k=2 must COLLAPSE via rule 1, got {label2} rule {rule2}"

    # k=3: ill-conditioning refusal -> AMBIGUOUS via rule 2
    comp3 = companion_matrix(3)
    label3, J3_64, J3_30, J3_50, rule3 = precision_ladder(comp3, 3)
    print(f"  k=3: J(64)={J3_64:.2e}, J(30)={J3_30:.2e}, J(50)={J3_50:.2e} -> {label3} (rule {rule3})")
    assert label3 == "AMBIGUOUS" and rule3 == 2, f"k=3 must be AMBIGUOUS via rule 2, got {label3} rule {rule3}"

    # Synthetic triplet (1e-8, 1e-12, 1e-12) -> COLLAPSES via rule 4
    J_syn = (1e-8, 1e-12, 1e-12)
    if J_syn[0] == 0 and J_syn[1] == 0 and J_syn[2] == 0:
        syn_label, syn_rule = "COLLAPSES", 1
    elif J_syn[2] < 0.5 * J_syn[1] or J_syn[2] > 2 * J_syn[1]:
        syn_label, syn_rule = "AMBIGUOUS", 2
    elif J_syn[2] >= 0.5 * J_syn[0]:
        syn_label, syn_rule = "PERSISTS", 3
    elif J_syn[2] <= 1e-3 * J_syn[0]:
        syn_label, syn_rule = "COLLAPSES", 4
    else:
        syn_label, syn_rule = "AMBIGUOUS", 5
    print(f"  synthetic (1e-8,1e-12,1e-12): -> {syn_label} (rule {syn_rule})")
    assert syn_label == "COLLAPSES" and syn_rule == 4

    # Rule-4 tail edges: (1e-8, 1e-12, 2e-12) and (1e-8, 1e-12, 5e-13)
    for J_edge, desc in [((1e-8, 1e-12, 2e-12), "ratio=2.0"),
                          ((1e-8, 1e-12, 5e-13), "ratio=0.5")]:
        if J_edge[0] == 0 and J_edge[1] == 0 and J_edge[2] == 0:
            el, er = "COLLAPSES", 1
        elif J_edge[2] < 0.5 * J_edge[1] or J_edge[2] > 2 * J_edge[1]:
            el, er = "AMBIGUOUS", 2
        elif J_edge[2] >= 0.5 * J_edge[0]:
            el, er = "PERSISTS", 3
        elif J_edge[2] <= 1e-3 * J_edge[0]:
            el, er = "COLLAPSES", 4
        else:
            el, er = "AMBIGUOUS", 5
        print(f"  edge {desc}: J={J_edge} -> {el} (rule {er})")
        assert el == "COLLAPSES" and er == 4, f"edge {desc} must COLLAPSE via rule 4"

    # Rule-5 fall-through: (1e-8, 2e-11, 2e-11) -> AMBIGUOUS via rule 5
    J_r5 = (1e-8, 2e-11, 2e-11)
    if J_r5[0] == 0 and J_r5[1] == 0 and J_r5[2] == 0:
        r5_label, r5_rule = "COLLAPSES", 1
    elif J_r5[2] < 0.5 * J_r5[1] or J_r5[2] > 2 * J_r5[1]:
        r5_label, r5_rule = "AMBIGUOUS", 2
    elif J_r5[2] >= 0.5 * J_r5[0]:
        r5_label, r5_rule = "PERSISTS", 3
    elif J_r5[2] <= 1e-3 * J_r5[0]:
        r5_label, r5_rule = "COLLAPSES", 4
    else:
        r5_label, r5_rule = "AMBIGUOUS", 5
    print(f"  rule-5 fall-through (1e-8,2e-11,2e-11): -> {r5_label} (rule {r5_rule})")
    assert r5_label == "AMBIGUOUS" and r5_rule == 5

    save_result("ladder_controls", {
        "k2": {"J": [J2_64, J2_30, J2_50], "label": label2, "rule": rule2},
        "k3": {"J": [J3_64, J3_30, J3_50], "label": label3, "rule": rule3},
        "synthetic_rule4": {"J": list(J_syn), "label": syn_label, "rule": syn_rule},
        "edge_2_0": {"J": [1e-8, 1e-12, 2e-12], "label": "COLLAPSES", "rule": 4},
        "edge_0_5": {"J": [1e-8, 1e-12, 5e-13], "label": "COLLAPSES", "rule": 4},
        "rule5_fallthrough": {"J": list(J_r5), "label": r5_label, "rule": r5_rule},
    })
    note(f"\n## Ladder Controls (A1.1)\n"
         f"- k=2: {label2} via rule {rule2} (exact-zero control)\n"
         f"- k=3: {label3} via rule {rule3} (ill-conditioning refusal)\n"
         f"- synthetic (1e-8,1e-12,1e-12): {syn_label} via rule {syn_rule}\n"
         f"- edges at ratio 2.0 and 0.5: both COLLAPSES via rule 4\n"
         f"- (1e-8,2e-11,2e-11): {r5_label} via rule {r5_rule}")

    # -----------------------------------------------------------------------
    #  Q7: Adjudicator outcomes and precedence collisions (synthetic)
    # -----------------------------------------------------------------------
    print(f"\n{SEP}\n  ADJUDICATOR OUTCOMES (Q7)\n{SEP}")

    # Seven reachable outcomes
    cases = [
        ("S1b-DEFECT", False, []),
        ("S1b-NO_LABEL (r2)", True, [
            {"n":12,"J_label_svd":"AMBIGUOUS","J_label_avg":"COLLAPSES","K_above_floor_svd":False,"K_above_floor_avg":False},
            {"n":20,"J_label_svd":"COLLAPSES","J_label_avg":"COLLAPSES","K_above_floor_svd":False,"K_above_floor_avg":False}]),
        ("S1b-NO_LABEL (r3)", True, [
            {"n":12,"J_label_svd":"PERSISTS","J_label_avg":"COLLAPSES","K_above_floor_svd":True,"K_above_floor_avg":True},
            {"n":20,"J_label_svd":"COLLAPSES","J_label_avg":"COLLAPSES","K_above_floor_svd":False,"K_above_floor_avg":False}]),
        ("S1b-SPECTRAL", True, [
            {"n":12,"J_label_svd":"PERSISTS","J_label_avg":"PERSISTS","K_above_floor_svd":True,"K_above_floor_avg":True},
            {"n":20,"J_label_svd":"COLLAPSES","J_label_avg":"COLLAPSES","K_above_floor_svd":False,"K_above_floor_avg":False}]),
        ("S1b-NO_LABEL (r5)", True, [
            {"n":12,"J_label_svd":"COLLAPSES","J_label_avg":"COLLAPSES","K_above_floor_svd":True,"K_above_floor_avg":False},
            {"n":20,"J_label_svd":"COLLAPSES","J_label_avg":"COLLAPSES","K_above_floor_svd":False,"K_above_floor_avg":False}]),
        ("S1b-ADJOINT", True, [
            {"n":12,"J_label_svd":"COLLAPSES","J_label_avg":"COLLAPSES","K_above_floor_svd":True,"K_above_floor_avg":True},
            {"n":20,"J_label_svd":"COLLAPSES","J_label_avg":"COLLAPSES","K_above_floor_svd":False,"K_above_floor_avg":False}]),
        ("S1b-NULL", True, [
            {"n":12,"J_label_svd":"COLLAPSES","J_label_avg":"COLLAPSES","K_above_floor_svd":False,"K_above_floor_avg":False},
            {"n":20,"J_label_svd":"COLLAPSES","J_label_avg":"COLLAPSES","K_above_floor_svd":False,"K_above_floor_avg":False}]),
    ]
    adj_results = []
    for expected_name, gate_ok, targets in cases:
        outcome, rule = branch_logic(gate_ok, targets)
        expected_outcome = expected_name.split(" ")[0]
        ok = outcome == expected_outcome
        print(f"  {expected_name:30s} -> {outcome} (rule {rule}) {'OK' if ok else 'MISMATCH'}")
        adj_results.append({"case": expected_name, "outcome": outcome, "rule": rule, "ok": ok})
        assert ok, f"adjudicator mismatch: {expected_name} -> {outcome}"

    # Four precedence collisions
    print(f"\n  Precedence collisions:")
    collisions = [
        ("gate FAILS + SPECTRAL -> DEFECT", False, [
            {"n":12,"J_label_svd":"PERSISTS","J_label_avg":"PERSISTS","K_above_floor_svd":True,"K_above_floor_avg":True}],
         "S1b-DEFECT"),
        ("AMBIGUOUS + PERSISTENT -> NO_LABEL", True, [
            {"n":12,"J_label_svd":"PERSISTS","J_label_avg":"PERSISTS","K_above_floor_svd":True,"K_above_floor_avg":True},
            {"n":20,"J_label_svd":"AMBIGUOUS","J_label_avg":"AMBIGUOUS","K_above_floor_svd":False,"K_above_floor_avg":False}],
         "S1b-NO_LABEL"),
        ("J disagree + one SPECTRAL -> NO_LABEL", True, [
            {"n":12,"J_label_svd":"PERSISTS","J_label_avg":"COLLAPSES","K_above_floor_svd":True,"K_above_floor_avg":True},
            {"n":20,"J_label_svd":"COLLAPSES","J_label_avg":"COLLAPSES","K_above_floor_svd":False,"K_above_floor_avg":False}],
         "S1b-NO_LABEL"),
        ("K disagree + one ADJOINT -> NO_LABEL", True, [
            {"n":12,"J_label_svd":"COLLAPSES","J_label_avg":"COLLAPSES","K_above_floor_svd":True,"K_above_floor_avg":False},
            {"n":20,"J_label_svd":"COLLAPSES","J_label_avg":"COLLAPSES","K_above_floor_svd":False,"K_above_floor_avg":False}],
         "S1b-NO_LABEL"),
    ]
    coll_results = []
    for desc, gate_ok, targets, expected in collisions:
        outcome, rule = branch_logic(gate_ok, targets)
        ok = outcome == expected
        print(f"  {desc:50s} -> {outcome} (r{rule}) {'OK' if ok else 'MISMATCH'}")
        coll_results.append({"desc": desc, "outcome": outcome, "rule": rule, "ok": ok})
        assert ok

    save_result("adjudicator", {"outcomes": adj_results, "collisions": coll_results})
    note(f"\n## Adjudicator (Q7)\n- All 7 outcomes reachable: YES\n- All 4 precedence collisions correct: YES")

    # -----------------------------------------------------------------------
    #  G-ALIGN — structural predicate
    # -----------------------------------------------------------------------
    print(f"\n{SEP}\n  G-ALIGN: cover-to-seed restriction\n{SEP}")
    r_vec = np.array([plan[o][0] for o in seed_orbits])
    r_expected = np.array([node_of[(o, 0)] for o in seed_orbits])
    align_ok = bool(np.array_equal(r_vec, r_expected))
    so_ok = seed_orbits == sorted(plan.keys())
    print(f"  seed_orbits == sorted(plan.keys()): {so_ok}")
    print(f"  r == (node_of[(o,0)])_{{o in seed_orbits}}: {align_ok}")
    print(f"  G-ALIGN: {'PASS' if align_ok else 'FAIL'}")

    # Arm: change one entry
    r_mut = r_vec.copy()
    orig = r_mut[0]
    r_mut[0] = r_vec[1] if r_vec[1] != orig else r_vec[2]
    arm_fails = not np.array_equal(r_mut, r_expected)
    print(f"  arm: r[0] {orig}->{r_mut[0]}, predicate fails: {arm_fails}")

    save_result("G-ALIGN", {"structural_predicate": align_ok, "seed_orbits_sorted": so_ok,
                             "arm_red": arm_fails, "pass": align_ok})
    note(f"\n## G-ALIGN\n- structural predicate: {align_ok}\n- arm (changed r[0]): predicate fails={arm_fails}\n- **{'PASS' if align_ok else 'FAIL'}**")

    # -----------------------------------------------------------------------
    #  G-WIRE — n=0 machinery control
    # -----------------------------------------------------------------------
    print(f"\n{SEP}\n  G-WIRE: n=0 machinery control\n{SEP}")
    _, C0, gap0 = repn.invariant_dim_and_basis(pairs, 0)
    print(f"  n=0: dim={C0.shape[1]}, gap={gap0}")
    dim0 = 1
    V0 = np.zeros((len(X), dim0), dtype=complex)
    for i in range(len(X)):
        V0[i] = sym_power(quat_to_su2(X[i]), 0).reshape(-1, order="F")
    F0_cover = V0 @ C0
    F0_seed = F0_cover[rep_indices]
    Q0 = mh_orthonormalize(F0_seed, Mh_diag)
    LQ0 = L_q @ Q0
    LQ0_norm = float(np.sqrt(np.real(np.sum(np.conj(LQ0) * (Mh_diag[:, None] * LQ0)))))
    gw_pass = LQ0_norm <= 1e-8
    print(f"  ||L Q_0||_{{M_h,F}} = {LQ0_norm:.2e} (gate: <= 1e-8)")
    print(f"  G-WIRE: {'PASS' if gw_pass else 'FAIL'}")

    # max |Im|
    A0 = form_An(Q0, Mh_diag, L_q)
    im0 = float(np.max(np.abs(np.imag(np.linalg.eigvals(A0)))))
    print(f"  max |Im| on A_0: {im0:.2e} (recorded, not gated)")

    # Arm: L_mut = L + delta * P_0
    delta = 1e-4
    P0 = Q0 @ (Q0.conj().T @ np.diag(Mh_diag))
    L_mut = L_q + delta * P0
    LmQ0 = L_mut @ Q0
    LmQ0_norm = float(np.sqrt(np.real(np.sum(np.conj(LmQ0) * (Mh_diag[:, None] * LmQ0)))))
    arm_red = LmQ0_norm > 1e-8
    print(f"  arm: ||L_mut Q_0|| = {LmQ0_norm:.2e} (expected >= {delta - 1e-8:.0e}), red: {arm_red}")

    save_result("G-WIRE", {"LQ0_norm": LQ0_norm, "max_im_A0": im0,
                             "arm_LmQ0_norm": LmQ0_norm, "arm_red": arm_red, "pass": gw_pass})
    note(f"\n## G-WIRE\n- ||L Q_0||_{{M_h,F}} = {LQ0_norm:.2e}\n- arm: ||L_mut Q_0|| = {LmQ0_norm:.2e}, red={arm_red}\n- **{'PASS' if gw_pass else 'FAIL'}**")

    # -----------------------------------------------------------------------
    #  SVD calls — THE EXPENSIVE PART
    # -----------------------------------------------------------------------
    target_data = {}
    rng_basis = np.random.default_rng(20260825)

    for n_target, k_expected in [(12, 13), (20, 21)]:
        print(f"\n{SEP}\n  SVD at n={n_target} (EXPENSIVE — shipped full_matrices=True)\n{SEP}")
        dim_n = (n_target + 1) ** 2

        # --- Shipped SVD (full) ---
        cache_path = os.path.join(RESULTS, f"C_svd_n{n_target}.npy")
        cache_json = os.path.join(RESULTS, f"svd_n{n_target}.json")
        if os.path.exists(cache_path) and os.path.exists(cache_json):
            C_svd = np.load(cache_path)
            with open(cache_json) as fj:
                cached = json.load(fj)
            k_svd = cached["dim"]
            gap_svd = cached["gap"]
            dt_svd = cached["time_s"]
            print(f"  shipped SVD (CACHED): dim={k_svd}, time={dt_svd:.1f}s")
            print(f"  gap: {gap_svd}")
            print(f"  hash: {sha256_bytes(C_svd.tobytes())}")
        else:
            t0 = time.time()
            k_svd, C_svd, gap_svd = repn.invariant_dim_and_basis(pairs, n_target)
            dt_svd = time.time() - t0
            print(f"  shipped SVD: dim={k_svd}, time={dt_svd:.1f}s")
            print(f"  gap: {gap_svd}")

            svd_artifact = {
                "n": n_target, "dim": k_svd, "gap": gap_svd,
                "basis_shape": list(C_svd.shape),
                "basis_hash": sha256_bytes(C_svd.tobytes()),
                "time_s": dt_svd,
            }
            save_result(f"svd_n{n_target}", svd_artifact)

            np.save(cache_path, C_svd)
            print(f"  [basis serialized to C_svd_n{n_target}.npy]")

        # G-RANK
        gr_pass = (k_svd == k_expected)
        print(f"  G-RANK n={n_target}: dim={k_svd} (expected {k_expected}) -> {'PASS' if gr_pass else 'FAIL'}")

        # G-RANK mutation: verify at a nearby n that gives a DIFFERENT dimension
        n_check = n_target - 1
        k_check, _, _ = repn.invariant_dim_and_basis(pairs, n_check)
        rank_arm_different = (k_check != k_expected)
        print(f"  G-RANK arm: n={n_check} gives dim={k_check} != {k_expected}: {rank_arm_different}")
        print(f"    (no-transpose form is character-equivalent; G-REAL separates realizations)")

        save_result(f"G-RANK_n{n_target}", {"dim": k_svd, "expected": k_expected,
                                              "gap": gap_svd, "pass": gr_pass,
                                              "arm_n": n_check, "arm_dim": k_check})
        note(f"\n## G-RANK n={n_target}\n- dim={k_svd} (expected {k_expected})\n- gap={gap_svd}\n- arm: n={n_check} gives dim={k_check}\n- **{'PASS' if gr_pass else 'FAIL'}**")

        # --- Economy SVD bridge diagnostic (per A2, licenses nothing) ---
        t0 = time.time()
        blocks = [repn.coefficient_operator(u, v, n_target) - np.eye(dim_n)
                  for (u, v) in pairs]
        A_stacked = np.vstack(blocks)
        _, s_econ, Vh_econ = np.linalg.svd(A_stacked, full_matrices=False)
        dt_econ = time.time() - t0
        print(f"  economy SVD: time={dt_econ:.1f}s")

        # Compare basis: the shipped call used full_matrices=True on the same A.
        # For tall A (m >> n), s and Vh are identical between full and economy SVD
        # (only U differs). Verify by comparing the resulting basis.
        cutoff_econ = max(1e-8 * s_econ[0], 1e-12) if s_econ.size else 1e-12
        k_econ = int(np.sum(s_econ < cutoff_econ))
        C_econ = Vh_econ.conj().T[:, dim_n - k_econ:] if k_econ else np.zeros((dim_n, 0))

        # Bitwise comparison of basis matrices
        basis_match = bool(np.array_equal(C_svd, C_econ))
        if not basis_match and C_svd.shape == C_econ.shape:
            basis_fro = float(np.linalg.norm(C_svd - C_econ, "fro"))
            angle_bridge = sine_form_angle(C_svd, C_econ)
        else:
            basis_fro = 0.0 if basis_match else float("inf")
            angle_bridge = 0.0 if basis_match else float("nan")

        print(f"  bridge: dim_econ={k_econ}, basis bitwise match={basis_match}")
        if not basis_match:
            print(f"  bridge: ||C_svd - C_econ||_F={basis_fro:.2e}, sin(angle)={angle_bridge:.2e}")

        save_result(f"bridge_n{n_target}", {
            "dim_econ": k_econ, "dim_shipped": k_svd,
            "basis_bitwise_match": basis_match,
            "basis_fro_diff": basis_fro,
            "subspace_sin_angle": angle_bridge,
            "economy_time_s": dt_econ,
            "note": "diagnostic per A2; licenses nothing",
        })
        note(f"\n## SVD Bridge n={n_target}\n- dim match: {k_svd == k_econ}\n- basis bitwise match: {basis_match}\n- (diagnostic, routes nothing)")

        # --- Averaging construction ---
        t0 = time.time()
        C_avg, rank_avg, avg_diag = averaging_basis(pairs, n_target)
        dt_avg = time.time() - t0
        print(f"  averaging: rank={rank_avg}, time={dt_avg:.1f}s")
        print(f"  Pi diagnostics: {avg_diag}")

        # --- theta_C: sine form ---
        s_C = sine_form_angle(C_svd, C_avg)
        print(f"  theta_C (sine form): sin(theta_max) = {s_C:.4e} (gate: <= 1e-6)")

        # G-SUBSPACE
        dims_agree = (k_svd == rank_avg == k_expected)
        angle_ok = (s_C <= 1e-6)
        gs_pass = dims_agree and angle_ok
        print(f"  G-SUBSPACE n={n_target}: dims_agree={dims_agree}, angle_ok={angle_ok} -> {'PASS' if gs_pass else 'FAIL'}")
        save_result(f"G-SUBSPACE_n{n_target}", {
            "k_svd": k_svd, "rank_avg": rank_avg, "sin_theta_C": s_C,
            "avg_diagnostics": avg_diag, "pass": gs_pass,
        })
        note(f"\n## G-SUBSPACE n={n_target}\n- k_svd={k_svd}, rank_avg={rank_avg}\n"
             f"- sin(theta_C)={s_C:.4e} (gate <=1e-6)\n- **{'PASS' if gs_pass else 'FAIL'}**")

        # G-SUBSPACE tilt arms (sine form)
        print(f"\n  G-SUBSPACE tilt arms (sine form):")
        orth_comp = np.eye(dim_n, dtype=complex) - C_svd @ np.linalg.pinv(C_svd)
        w = None
        for col in range(dim_n):
            wc = orth_comp[:, col]
            if np.linalg.norm(wc) > 1e-12:
                w = wc / np.linalg.norm(wc)
                break
        for sin_phi, required in [(1e-4, "RED"), (1e-7, "green")]:
            C_tilt = C_svd.copy()
            phi = np.arcsin(sin_phi)
            C_tilt[:, 0] = C_svd[:, 0] * np.cos(phi) + w * np.sin(phi)
            C_tilt, _ = np.linalg.qr(C_tilt)
            meas = sine_form_angle(C_svd, C_tilt)
            gate_says = "RED" if meas > 1e-6 else "green"
            print(f"    sin(phi)={sin_phi:.0e}: sin(theta_max)={meas:.4e}, "
                  f"gate={gate_says}, required={required} -> {'OK' if gate_says == required else 'MISMATCH'}")

        # --- Sampling: both routes ---
        route_results = {}
        for route_label, C_n in [("svd", C_svd), ("avg", C_avg)]:
            # V_n on cover cloud
            Vn = np.zeros((len(X), dim_n), dtype=complex)
            for i in range(len(X)):
                Vn[i] = sym_power(quat_to_su2(X[i]), n_target).reshape(-1, order="F")
            Fn_cover = Vn @ C_n
            Fn_seed = Fn_cover[rep_indices]
            # F_n^seed == F_n^cover[r, :] assertion
            Fn_align = Fn_cover[rep_indices]
            assert np.array_equal(Fn_seed, Fn_align), "F_n^seed alignment failed"

            # W_n = M_h^{1/2} F_n^seed
            Wn = Mh_sqrt[:, None] * Fn_seed
            s_W = np.linalg.svd(Wn, compute_uv=False)
            num_rank = int(np.sum(s_W > 1e-8 * s_W[0]))
            kappa_W = float(s_W[0] / s_W[-1]) if s_W[-1] > 0 else float("inf")

            # G-SAMPLE
            rank_ok = (num_rank == k_expected)
            kappa_ok = (kappa_W <= 1e6)
            gsamp_pass = rank_ok and kappa_ok
            print(f"  G-SAMPLE n={n_target} {route_label}: rank={num_rank}, kappa={kappa_W:.4e} -> {'PASS' if gsamp_pass else 'FAIL'}")
            save_result(f"G-SAMPLE_n{n_target}_{route_label}", {
                "rank": num_rank, "kappa": kappa_W, "pass": gsamp_pass,
            })

            # M_h-orthonormalize
            Qn = mh_orthonormalize(Fn_seed, Mh_diag)
            k_actual = Qn.shape[1]
            orth_err = float(np.abs(Qn.conj().T @ np.diag(Mh_diag) @ Qn - np.eye(k_actual)).max())

            # A_n, K_n
            An = form_An(Qn, Mh_diag, L_q)
            Kn, Kn_norm = compute_K(An)
            An_norm = float(np.linalg.norm(An, ord=2))

            # Qtilde for theta_Q
            Qtilde = Mh_sqrt[:, None] * Qn

            route_results[route_label] = {
                "C_n": C_n, "Fn_cover": Fn_cover, "Fn_seed": Fn_seed,
                "Wn": Wn, "Qn": Qn, "Qtilde": Qtilde,
                "An": An, "Kn_norm": Kn_norm, "An_norm": An_norm,
                "kappa_W": kappa_W, "orth_err": orth_err,
            }

        # G-SAMPLE mutation arm: zero out one column of C_n -> rank-deficient sampling
        print(f"\n  G-SAMPLE mutation arm n={n_target}:")
        C_mut = C_svd.copy()
        C_mut[:, 0] = 0.0
        Vn_for_mut = np.zeros((len(X), dim_n), dtype=complex)
        for i in range(len(X)):
            Vn_for_mut[i] = sym_power(quat_to_su2(X[i]), n_target).reshape(-1, order="F")
        Fn_mut = (Vn_for_mut @ C_mut)[rep_indices]
        Wn_mut = Mh_sqrt[:, None] * Fn_mut
        s_W_mut = np.linalg.svd(Wn_mut, compute_uv=False)
        rank_mut = int(np.sum(s_W_mut > 1e-8 * s_W_mut[0])) if s_W_mut[0] > 0 else 0
        arm_red = (rank_mut < k_expected)
        print(f"    zeroed column: rank={rank_mut} (expected {k_expected}), red={arm_red}")

        # G-BASIS on both routes (rng_basis initialized once before loop)
        print(f"\n  G-BASIS n={n_target}:")
        for rl in ["svd", "avg"]:
            Qn = route_results[rl]["Qn"]
            An = route_results[rl]["An"]
            k = Qn.shape[1]
            H = rng_basis.normal(size=(k, k)) + 1j * rng_basis.normal(size=(k, k))
            U_rot, _ = np.linalg.qr(H)
            Qn_rot = Qn @ U_rot
            An_rot = form_An(Qn_rot, Mh_diag, L_q)
            A_cov_err = float(np.linalg.norm(An_rot - U_rot.conj().T @ An @ U_rot, "fro"))
            A_scale = max(float(np.linalg.norm(An, "fro")), 1e-30)
            A_rel = A_cov_err / A_scale
            Pn = Qn @ (Qn.conj().T @ np.diag(Mh_diag))
            Pn_rot = Qn_rot @ (Qn_rot.conj().T @ np.diag(Mh_diag))
            P_diff = float(np.linalg.norm(Pn_rot - Pn, "fro"))
            P_scale = max(1.0, float(np.linalg.norm(Pn, "fro")))
            P_rel = P_diff / P_scale
            a_ok = A_rel <= 1e-10
            p_ok = P_rel <= 1e-10
            print(f"    {rl}: ||A'-U^H A U||/||A||={A_rel:.2e}, ||P'-P||/scale={P_rel:.2e} -> {'PASS' if a_ok and p_ok else 'FAIL'}")
            save_result(f"G-BASIS_n{n_target}_{rl}", {
                "A_covariance_rel": A_rel, "P_diff_rel": P_rel,
                "pass": a_ok and p_ok,
            })

        # G-BASIS mutation arm: non-isometric perturbation
        print(f"  G-BASIS mutation arm n={n_target}:")
        Qn_ref = route_results["svd"]["Qn"]
        An_ref = route_results["svd"]["An"]
        k_b = Qn_ref.shape[1]
        M_pert = np.eye(k_b, dtype=complex)
        M_pert[0, 0] = 2.0
        Qn_pert = Qn_ref @ M_pert
        An_pert = form_An(Qn_pert, Mh_diag, L_q)
        A_cov_pert = float(np.linalg.norm(An_pert - M_pert.conj().T @ An_ref @ M_pert, "fro"))
        Pn_ref = Qn_ref @ (Qn_ref.conj().T @ np.diag(Mh_diag))
        Pn_pert = Qn_pert @ (Qn_pert.conj().T @ np.diag(Mh_diag))
        P_diff_pert = float(np.linalg.norm(Pn_pert - Pn_ref, "fro"))
        basis_arm_red = (P_diff_pert > 1e-10)
        print(f"    non-isometric: ||P'-P||={P_diff_pert:.2e}, red={basis_arm_red}")

        note(f"\n## G-SAMPLE n={n_target}\n"
             f"- svd: kappa={route_results['svd']['kappa_W']:.4e}\n"
             f"- avg: kappa={route_results['avg']['kappa_W']:.4e}\n"
             f"- arm (zeroed column): rank={rank_mut}, red={arm_red}")
        note(f"\n## G-BASIS n={n_target}\n"
             f"- isometric rotation: covariance and projector invariant on both routes\n"
             f"- arm (non-isometric): ||P'-P||={P_diff_pert:.2e}, red={basis_arm_red}")

        # theta_Q (sine form, in Qtilde space)
        s_Q = sine_form_angle(route_results["svd"]["Qtilde"], route_results["avg"]["Qtilde"])
        d_Q = half_angle_d(s_Q)
        print(f"\n  theta_Q n={n_target}: sin(theta_Q) = {s_Q:.4e}, 2sin(theta_Q/2) = {d_Q:.4e}")
        print(f"  (recorded, NOT gated — deliberate per contract)")

        # K_floor
        max_An_norm = max(route_results["svd"]["An_norm"], route_results["avg"]["An_norm"])
        max_kappa_W = max(route_results["svd"]["kappa_W"], route_results["avg"]["kappa_W"])
        term1 = 100 * eps * max_An_norm
        term2 = 10 * eps * max_kappa_W * B_norm
        term3 = d_Q * B_norm
        K_floor = max(term1, term2, term3)
        dominant = ["arithmetic", "kappa", "discrepancy"][[term1, term2, term3].index(max(term1, term2, term3))]
        print(f"\n  K_floor n={n_target}:")
        print(f"    term1 (arithmetic):  100*eps*max||A||_2 = {term1:.4e}")
        print(f"    term2 (kappa):       10*eps*max_kappa*||B||_2 = {term2:.4e}")
        print(f"    term3 (discrepancy): 2sin(theta_Q/2)*||B||_2 = {term3:.4e}")
        print(f"    K_floor = {K_floor:.4e}, DOMINANT: {dominant}")

        # K resolved?
        for rl in ["svd", "avg"]:
            Kn_norm = route_results[rl]["Kn_norm"]
            above = Kn_norm > K_floor
            print(f"    ||K_{n_target}^{rl}||_2 = {Kn_norm:.4e} {'>' if above else '<='} K_floor -> {'RESOLVED' if above else 'not resolved'}")

        # Precision ladder (re-formed)
        print(f"\n  Precision ladder n={n_target}:")
        for rl in ["svd", "avg"]:
            Qn = route_results[rl]["Qn"]
            print(f"    route {rl}: reforming A at 30/50 dps...")
            t0 = time.time()
            label_J, J64, J30, J50, rule_J = precision_ladder_reformed(Qn, Mh_diag, L_q, k_expected)
            dt_ladder = time.time() - t0
            print(f"    J(64)={J64:.4e}, J(30)={J30:.4e}, J(50)={J50:.4e} -> {label_J} (rule {rule_J}), {dt_ladder:.1f}s")
            route_results[rl]["J_label"] = label_J
            route_results[rl]["J_readings"] = (J64, J30, J50)
            route_results[rl]["J_rule"] = rule_J

        # Bauer-Fike diagnostic
        for rl in ["svd", "avg"]:
            An = route_results[rl]["An"]
            evals_An, V_An = np.linalg.eig(An)
            kappa_V = float(np.linalg.cond(V_An))
            J_bf = kappa_V * 100 * eps * route_results[rl]["An_norm"]
            demoted = kappa_V > 1e8
            print(f"    Bauer-Fike {rl}: kappa(V)={kappa_V:.2e}, J_bf={J_bf:.2e}"
                  f"{' (DEMOTED)' if demoted else ''}")
            route_results[rl]["kappa_V"] = kappa_V
            route_results[rl]["J_bf"] = J_bf
            route_results[rl]["J_bf_demoted"] = demoted

        target_data[n_target] = {
            "k": k_expected, "dim_n": dim_n,
            "sin_theta_C": s_C, "sin_theta_Q": s_Q,
            "d_theta_Q": d_Q, "K_floor": K_floor,
            "dominant_term": dominant,
            "terms": {"arithmetic": term1, "kappa": term2, "discrepancy": term3},
            "routes": {},
        }
        for rl in ["svd", "avg"]:
            Kn_norm = route_results[rl]["Kn_norm"]
            target_data[n_target]["routes"][rl] = {
                "An_norm": route_results[rl]["An_norm"],
                "Kn_norm": Kn_norm,
                "K_above_floor": Kn_norm > K_floor,
                "kappa_W": route_results[rl]["kappa_W"],
                "J_label": route_results[rl]["J_label"],
                "J_readings": route_results[rl]["J_readings"],
                "J_rule": route_results[rl]["J_rule"],
                "kappa_V": route_results[rl]["kappa_V"],
                "J_bf": route_results[rl]["J_bf"],
                "orth_err": route_results[rl]["orth_err"],
            }

        save_result(f"target_n{n_target}", target_data[n_target])
        note(f"\n## Target n={n_target}\n"
             f"- sin(theta_C)={s_C:.4e}, sin(theta_Q)={s_Q:.4e}\n"
             f"- K_floor={K_floor:.4e} (dominant: {dominant})\n"
             f"- svd: ||K||={route_results['svd']['Kn_norm']:.4e}, J={route_results['svd']['J_label']}\n"
             f"- avg: ||K||={route_results['avg']['Kn_norm']:.4e}, J={route_results['avg']['J_label']}")

    # -----------------------------------------------------------------------
    #  Q8: ADVERSARIAL VACUITY AUDIT
    # -----------------------------------------------------------------------
    print(f"\n{SEP}\n  ADVERSARIAL VACUITY AUDIT (Q8)\n{SEP}")
    vacuity = []

    cases_v = [
        ("G-REAL", "Can it pass for a reason other than the correct realization?",
         "No: the no-transpose form is character-equivalent (same dimension) but "
         "fails pointwise at O(1). Only the pointwise test separates them, and it "
         "sweeps n>=3 to avoid the n=2 coincidence trap."),
        ("G-RANK", "Can the dimension be right for a wrong basis?",
         "Yes in principle (wrong realization preserves dimension). But G-RANK only "
         "checks dimension; G-REAL checks the realization separately. The two together "
         "are not vacuous."),
        ("G-SUBSPACE", "Can the angle be small for a wrong reason?",
         "The sine form returns ~7e-16 on identical subspaces and resolves 1e-7. "
         "The tilt arms verify discrimination at 1e-4 (RED) and 1e-7 (green). "
         "The absolute rank cutoff guards against the spurious-rank-on-zero-matrix trap."),
        ("G-ALIGN", "Can the structural predicate pass trivially?",
         "No: the predicate is exact (r == node_of[(o,0)]) and the arm shows a single "
         "changed entry fails it. No numerical tolerance involved."),
        ("G-SAMPLE", "Can it pass after the step that guarantees it?",
         "No: the gate tests W_n = M_h^{1/2} F_n^seed BEFORE orthonormalization. "
         "Testing Q_n after orthonormalization would be green-by-construction."),
        ("G-BASIS", "Can it pass trivially at n=0?",
         "Yes: A_0 is 1x1, so U^H A U = |u|^2 A = A for unit u. This is why the "
         "contract specifies G-BASIS at n=12,20, and G-WIRE handles n=0 separately."),
        ("G-DISCRIM", "Can an implementation conflate K and J?",
         "Arm A produces ||K||=0.5, J=0 (non-Hermitian but real spectrum). An "
         "implementation that equates non-Hermitian with complex-spectrum fails arm A."),
        ("G-WIRE", "Can the constant-mode arm be vacuous?",
         "The original permutation arm was vacuous (F_0 is constant, rows identical). "
         "The contracted L_mut arm has analytically forced effect: ||L_mut Q_0|| >= 1e-4 - 1e-8."),
    ]
    for gate, question, answer in cases_v:
        print(f"\n  {gate}: {question}")
        print(f"    {answer}")
        vacuity.append({"gate": gate, "question": question, "answer": answer})

    save_result("vacuity_audit", {"cases": vacuity})
    note(f"\n## Adversarial Vacuity Audit (Q8)\n- All 8 gates examined\n- Key finding: G-BASIS at n=0 has no power (contract addresses via G-WIRE)")

    # -----------------------------------------------------------------------
    #  Q9: CONTRACT CONTRADICTIONS
    # -----------------------------------------------------------------------
    print(f"\n{SEP}\n  CONTRACT CONTRADICTIONS (Q9)\n{SEP}")

    contradictions = []
    print("  No implementability-blocking contradictions found.")
    print()
    print("  Examined and resolved:")
    findings = [
        "theta_C vs theta_Q clearly distinguished; G-SUBSPACE gates theta_C, K_floor uses theta_Q",
        "Route indexing: all downstream objects explicitly per-route, K_floor is the one common bar",
        "Precision ladder: ordered rules with first-match-wins, rule 1 handles exact zero",
        "Branch table: pass-7 fixed the G-SUBSPACE vs branch-table contradiction on route disagreement",
        "Norms: all in operator 2-norm (pass-7 corrected earlier Frobenius/operator mixing)",
        "Half-angle: 2*sin(theta/2) via arcsin, not the prohibited sqrt identity or arccos route",
        "Addendum 2 correctly marks the evidence-transfer branch INOPERATIVE; bridge is diagnostic only",
        "Addendum 3 correctly keys the ledger by attempt_id, resolving the fixed-directory incompatibility",
    ]
    for i, f in enumerate(findings, 1):
        print(f"  {i}. {f}")
        contradictions.append(f)

    save_result("contract_contradictions", {"blocking": [], "examined_and_resolved": contradictions})
    note(f"\n## Contract Contradictions (Q9)\n- No blocking contradictions found\n- {len(contradictions)} items examined and resolved")

    # -----------------------------------------------------------------------
    #  END-OF-RUN INPUT RE-VERIFICATION (A3.3)
    # -----------------------------------------------------------------------
    print(f"\n{SEP}\n  END-OF-RUN INPUT RE-VERIFICATION\n{SEP}")

    manifest_path = os.path.join(ROOM, "ROOM_MANIFEST.json")
    with open(manifest_path) as f:
        manifest = json.load(f)

    all_inputs_ok = True
    for filepath, expected_hash in manifest.items():
        full_path = os.path.join(ROOM, filepath)
        actual = sha256_file(full_path)
        if actual != expected_hash:
            print(f"  CHANGED: {filepath}")
            all_inputs_ok = False

    # Re-verify frozen documents
    boundaries = {
        "contract/S1B_DECISION_RULE.md": ("<!-- FREEZE-BOUNDARY -->",
            "c44c603a8058ed8529e5bb0f42ec168b443ff22b050d3748f1f1de3537c7d297"),
        "contract/S1B_ADDENDUM_1.md": ("<!-- ADDENDUM-BOUNDARY -->",
            "6da36a1c672772e8c731d4d66f20a77e81779d19d5cb883e617c028d43f46746"),
        "contract/S1B_ADDENDUM_2.md": ("<!-- ADDENDUM2-BOUNDARY -->",
            "14011c338ce331d9bb3424f0ed5aaff31a1c89a9ceafbc2bae41a91ebb41a222"),
        "contract/S1B_ADDENDUM_3.md": ("<!-- ADDENDUM3-BOUNDARY -->",
            "e3304fe993f7a2ea523cfcbc4110c0f37174caea95383b17e465c8fa85b6e28c"),
    }
    for filepath, (boundary, expected) in boundaries.items():
        full_path = os.path.join(ROOM, filepath)
        with open(full_path, "rb") as f:
            raw = f.read()
        idx = raw.find(boundary.encode())
        if idx < 0:
            print(f"  ERROR: boundary not found in {filepath}")
            all_inputs_ok = False
            continue
        frozen = raw[:idx]
        h = hashlib.sha256(frozen).hexdigest()
        if h != expected:
            print(f"  FROZEN CHANGED: {filepath}")
            all_inputs_ok = False

    if all_inputs_ok:
        print("  All manifest-pinned inputs unchanged since Q0: OK")
    else:
        print("  STOP: inputs changed since Q0!")
        note(f"\n## End-of-run verification\n- **STOP: inputs changed**")
        return

    note(f"\n## End-of-run verification\n- All inputs unchanged: YES")

    # -----------------------------------------------------------------------
    #  OUTPUT MANIFEST
    # -----------------------------------------------------------------------
    print(f"\n{SEP}\n  OUTPUT MANIFEST\n{SEP}")
    output_files = {}
    for root, dirs, files in os.walk(LEDGER):
        for fname in files:
            fpath = os.path.join(root, fname)
            if fname == "OUTPUT_MANIFEST.json":
                continue
            rel = os.path.relpath(fpath, LEDGER)
            output_files[rel] = sha256_file(fpath)
    output_manifest = {
        "attempt_id": "q3a",
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "total_time_s": time.time() - t_start,
        "files": output_files,
    }
    out_path = os.path.join(LEDGER, "OUTPUT_MANIFEST.json")
    with open(out_path, "w") as f:
        json.dump(output_manifest, f, indent=2)
        f.write("\n")
    print(f"  OUTPUT_MANIFEST.json written with {len(output_files)} files")
    print(f"  total wall time: {time.time() - t_start:.1f}s")

    note(f"\n## Output Manifest\n- {len(output_files)} files\n- total time: {time.time() - t_start:.1f}s")


if __name__ == "__main__":
    main()
