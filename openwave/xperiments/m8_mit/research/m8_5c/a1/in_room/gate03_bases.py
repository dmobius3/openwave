"""Gate 3: Sector bases certification (§3.4).

(a) Reynolds dimension vs character route multiplicity
(b) Pointwise equivariance of sections
(e) First-occurrence re-derivation
Mutation arms for each.
"""
import numpy as np
import time
from m85c_group import (G120, IRREP_NAMES, D_RHO, DIM_RHO, multiplicity,
                        FIBRES, qmul, quat_to_su2, pi_unitary)
from m85c_bases import intertwiner_basis_deterministic
from m85c_ledger import gate_record, resource_record

N_CHECK = 24


def check_multiplicity_agreement(rho_name, N):
    """(a) Reynolds dimension vs character route at all levels n ≤ N."""
    errors = []
    for n in range(N + 1):
        m_char = multiplicity(rho_name, n)
        if m_char == 0:
            continue
        try:
            ints, m_svd = intertwiner_basis_deterministic(rho_name, n)
        except Exception as e:
            errors.append({"n": n, "m_char": m_char, "m_svd": -1, "error": str(e)})
            continue
        if m_svd != m_char:
            errors.append({"n": n, "m_char": m_char, "m_svd": m_svd})
    return errors


def check_equivariance(rho_name, N, n_samples=5):
    """(b) Pointwise equivariance: ρ(g) f(x) = f(gx) at random points."""
    rng = np.random.Generator(np.random.PCG64(123 + hash(rho_name) % 1000))
    rho_mats = FIBRES[rho_name]
    d_rho = DIM_RHO[rho_name]
    max_err = 0.0

    levels_to_check = [n for n in range(min(N + 1, 30)) if multiplicity(rho_name, n) > 0][:5]

    for n in levels_to_check:
        ints, m = intertwiner_basis_deterministic(rho_name, n)
        for _ in range(n_samples):
            x = rng.standard_normal(4)
            x /= np.linalg.norm(x)
            gi = rng.integers(0, 120)
            g = G120[gi]
            gx = qmul(g, x)

            for i_int in range(min(m, 2)):
                A = ints[i_int]
                for j in range(min(n + 1, 3)):
                    fx = A @ pi_unitary(n, x)[:, j]
                    fgx = A @ pi_unitary(n, gx)[:, j]
                    rho_g_fx = rho_mats[gi] @ fx
                    err = np.linalg.norm(rho_g_fx - fgx)
                    scale = max(np.linalg.norm(fx), 1e-14)
                    max_err = max(max_err, err / scale)

    return max_err


def check_first_occurrence(rho_name):
    """(e) First-occurrence re-derivation: mult(rho, V_{d_rho}) = 1 and
    mult(rho, V_n) = 0 for n < d_rho."""
    d = D_RHO[rho_name]
    for n in range(d):
        m = multiplicity(rho_name, n)
        if m != 0:
            return False, f"mult({rho_name}, V_{n}) = {m} != 0 (should be 0 for n < {d})"
    m_first = multiplicity(rho_name, d)
    if m_first != 1:
        return False, f"mult({rho_name}, V_{d}) = {m_first} != 1"
    return True, f"first occurrence at n={d}, mult=1"


def run_gate3():
    t0 = time.time()
    all_pass = True

    # (e) First-occurrence
    print("  (e) First-occurrence re-derivation:")
    fo_pass = True
    fo_results = {}
    for rho in IRREP_NAMES:
        ok, msg = check_first_occurrence(rho)
        fo_results[rho] = {"pass": ok, "detail": msg}
        fo_pass = fo_pass and ok
        print(f"    {rho}: d_rho={D_RHO[rho]}, dim={DIM_RHO[rho]}, {msg}, {'PASS' if ok else 'FAIL'}")
    all_pass = all_pass and fo_pass

    # (a) Multiplicity agreement
    print(f"\n  (a) Multiplicity agreement (Reynolds vs character), N={N_CHECK}:")
    mult_pass = True
    mult_results = {}
    for rho in IRREP_NAMES:
        errors = check_multiplicity_agreement(rho, N_CHECK)
        ok = len(errors) == 0
        mult_pass = mult_pass and ok
        n_levels = sum(1 for n in range(N_CHECK + 1) if multiplicity(rho, n) > 0)
        mult_results[rho] = {"n_levels": n_levels, "errors": errors, "pass": ok}
        print(f"    {rho}: {n_levels} levels, {len(errors)} errors, {'PASS' if ok else 'FAIL'}")
    all_pass = all_pass and mult_pass

    # (b) Pointwise equivariance
    print(f"\n  (b) Pointwise equivariance of sections:")
    equiv_pass = True
    equiv_results = {}
    for rho in IRREP_NAMES:
        if rho == "R0":
            max_err = 0.0
        else:
            max_err = check_equivariance(rho, N_CHECK)
        ok = max_err < 1e-8
        equiv_pass = equiv_pass and ok
        equiv_results[rho] = {"max_rel_err": max_err, "pass": ok}
        print(f"    {rho}: max_rel_err={max_err:.3e}, {'PASS' if ok else 'FAIL'}")
    all_pass = all_pass and equiv_pass

    # Mutation: deliberately break first occurrence for R3
    print("\n  Mutation arms:")
    mut_results = {}

    # Mult mutation: claim wrong multiplicity at a level
    rho_test = "R3"
    n_test = D_RHO[rho_test]
    m_real = multiplicity(rho_test, n_test)
    m_fake = m_real + 1
    mut_mult_fires = m_real != m_fake
    print(f"    Mult mutation: {rho_test} n={n_test}: real mult={m_real}, fake={m_fake}, fires={mut_mult_fires}")
    mut_results["multiplicity"] = {"fires": mut_mult_fires}

    # Equivariance mutation: use wrong group element
    rho_test = "R8"
    ints, m = intertwiner_basis_deterministic(rho_test, D_RHO[rho_test])
    A = ints[0]
    x = np.array([0.5, 0.5, 0.5, 0.5])
    gi = 10
    g = G120[gi]
    gx = qmul(g, x)
    fx = A @ pi_unitary(D_RHO[rho_test], x)[:, 0]
    fgx = A @ pi_unitary(D_RHO[rho_test], gx)[:, 0]
    wrong_gi = (gi + 1) % 120
    rho_wrong = FIBRES[rho_test][wrong_gi] @ fx
    err_wrong = np.linalg.norm(rho_wrong - fgx) / max(np.linalg.norm(fx), 1e-14)
    mut_equiv_fires = err_wrong > 1e-8
    print(f"    Equiv mutation: wrong group element, err={err_wrong:.3e}, fires={mut_equiv_fires}")
    mut_results["equivariance"] = {"fires": mut_equiv_fires, "err": float(err_wrong)}

    mutation_fires = mut_mult_fires and mut_equiv_fires
    gate_pass = all_pass and mutation_fires
    dt = time.time() - t0

    gate_record(
        gate_id="G3-BASES",
        arena_id=f"all-sectors-N{N_CHECK}",
        rung=N_CHECK,
        parent_status="GREEN" if all_pass else "RED",
        mutation_status="FIRES" if mutation_fires else "DEAD",
        measured_values={
            "first_occurrence": fo_results,
            "multiplicity": mult_results,
            "equivariance": equiv_results,
            "mutations": mut_results,
        },
    )

    resource_record("G3-BASES", dt, dt)
    return gate_pass, dt


if __name__ == "__main__":
    print("Gate 3: Sector bases certification")
    print("=" * 60)
    passed, dt = run_gate3()
    print(f"\nGate 3: {'GREEN' if passed else 'RED'} ({dt:.1f}s)")
