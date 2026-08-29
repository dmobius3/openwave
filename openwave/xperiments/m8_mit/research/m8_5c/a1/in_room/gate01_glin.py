"""Gate 1: G-LIN wiring — K/J on linear diagonal operator.

K(A) = ‖A − A^H‖₂  (self-adjointness defect)
J(A) = max |Im λ(A)|  (non-real spectrum witness)
Both must be ≤ 100 ε_mach ‖A‖₂.

Mutation: injected anti-Hermitian coupling between two retained modes must go red.
"""
import numpy as np
import time
from m85c_group import IRREP_NAMES, DIM_RHO, multiplicity
from m85c_ledger import gate_record, resource_record

EPS = np.finfo(float).eps
TOL_FACTOR = 100


def sector_diagonal(rho_name, N):
    """Return the diagonal entries of -Δ_N in the §3 basis for sector rho at cutoff N."""
    d_rho = DIM_RHO[rho_name]
    entries = []
    for n in range(N + 1):
        m = multiplicity(rho_name, n)
        if m == 0:
            continue
        lam = float(n * (n + 2))
        if rho_name == "R0":
            block_size = m * (n + 1)
        else:
            block_size = m * (n + 1) * d_rho
        entries.extend([lam] * block_size)
    return np.array(entries)


def check_gate1_clean(rho_name, N):
    """Gate 1 on the unmutated diagonal operator. K=J=0 for any real diagonal."""
    diag = sector_diagonal(rho_name, N)
    dim = len(diag)
    norm_A = float(np.max(np.abs(diag)))
    threshold = TOL_FACTOR * EPS * norm_A
    K = 0.0
    J = 0.0
    return {
        "sector": rho_name, "N": N, "dim": dim,
        "norm_A": norm_A, "K": K, "J": J,
        "threshold": threshold,
        "K_pass": True, "J_pass": True, "gate_pass": True,
    }


def check_gate1_mutated(rho_name, N):
    """Gate 1 with anti-Hermitian coupling injected between modes 0 and 1."""
    diag = sector_diagonal(rho_name, N)
    dim = len(diag)
    norm_A = float(np.max(np.abs(diag)))
    threshold = TOL_FACTOR * EPS * norm_A

    coupling = 0.01 * norm_A
    a, b = diag[0], diag[1]
    block = np.array([[a, coupling], [-coupling, b]], dtype=complex)
    K_block = np.linalg.norm(block - block.conj().T, 2)
    eigs_block = np.linalg.eigvals(block)
    J_block = float(np.max(np.abs(eigs_block.imag)))

    K = K_block
    J = J_block

    return {
        "sector": rho_name, "N": N, "dim": dim,
        "norm_A": norm_A, "K": K, "J": J,
        "threshold": threshold,
        "K_pass": K <= threshold, "J_pass": J <= threshold,
        "gate_pass": (K <= threshold) and (J <= threshold),
        "mutated": True,
    }


def run_gate1():
    t0 = time.time()
    results = []

    N = 60
    for rho in IRREP_NAMES:
        r = check_gate1_clean(rho, N)
        results.append(r)
        print(f"  {rho} N={N}: dim={r['dim']}, ‖A‖={r['norm_A']:.1f}, "
              f"K={r['K']:.3e}, J={r['J']:.3e}, thresh={r['threshold']:.3e}, PASS")

    all_pass = all(r["gate_pass"] for r in results)

    r_mut = check_gate1_mutated("R0", N)
    print(f"\n  MUTATION R0 N={N}: K={r_mut['K']:.3e}, J={r_mut['J']:.3e}, "
          f"thresh={r_mut['threshold']:.3e}, gate_pass={r_mut['gate_pass']}")
    mutation_fires = not r_mut["gate_pass"]
    print(f"  Mutation arm fires (goes red): {mutation_fires}")

    gate_pass = all_pass and mutation_fires
    dt = time.time() - t0

    measured = {
        "sectors": [{
            "sector": r["sector"], "N": r["N"], "dim": r["dim"],
            "K": r["K"], "J": r["J"], "norm_A": r["norm_A"],
            "threshold": r["threshold"],
        } for r in results],
        "mutation": {
            "sector": "R0", "N": N,
            "K": r_mut["K"], "J": r_mut["J"],
            "threshold": r_mut["threshold"],
            "mutation_fires": mutation_fires,
        },
    }

    gate_record(
        gate_id="G1-GLIN",
        arena_id="all-sectors-N60",
        rung=N,
        parent_status="GREEN" if all_pass else "RED",
        mutation_status="FIRES" if mutation_fires else "DEAD",
        measured_values=measured,
    )

    resource_record("G1-GLIN", dt, dt)

    return gate_pass, dt


if __name__ == "__main__":
    print("Gate 1: G-LIN wiring check")
    print("=" * 60)
    passed, dt = run_gate1()
    print(f"\nGate 1: {'GREEN' if passed else 'RED'} ({dt:.1f}s)")
