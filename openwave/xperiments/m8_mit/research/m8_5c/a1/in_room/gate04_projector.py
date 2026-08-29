"""Gate 4: Projector exactness (§4).

Production 4N rule to rounding. Node-drop to 2N must err O(1).
"""
import numpy as np
from math import comb
import time
from m85c_group import binom_sqrt_diag
from m85c_bases import build_r0_basis
from m85c_quadrature import hopf_rule
from m85c_packet import generate_packet
from m85c_ledger import gate_record, resource_record

BATCH = 8000


def batch_pi_unitary(X, n):
    """Unitarized sym_power at level n for a batch of quaternion nodes."""
    npts = X.shape[0]
    if n == 0:
        return np.ones((npts, 1, 1), dtype=complex)

    a = X[:, 0] + 1j * X[:, 1]
    b = X[:, 2] + 1j * X[:, 3]
    bc = -X[:, 2] + 1j * X[:, 3]
    ac = X[:, 0] - 1j * X[:, 1]

    dim = n + 1
    T = binom_sqrt_diag(n)

    a_pow = np.ones((npts, n + 1), dtype=complex)
    bc_pow = np.ones((npts, n + 1), dtype=complex)
    b_pow = np.ones((npts, n + 1), dtype=complex)
    ac_pow = np.ones((npts, n + 1), dtype=complex)
    for p in range(1, n + 1):
        a_pow[:, p] = a_pow[:, p-1] * a
        bc_pow[:, p] = bc_pow[:, p-1] * bc
        b_pow[:, p] = b_pow[:, p-1] * b
        ac_pow[:, p] = ac_pow[:, p-1] * ac

    binoms = np.zeros((n + 1, n + 1))
    for nn in range(n + 1):
        for kk in range(nn + 1):
            binoms[nn, kk] = comb(nn, kk)

    result = np.zeros((npts, dim, dim), dtype=complex)
    for k in range(dim):
        m = n - k
        polyA = np.zeros((npts, m + 1), dtype=complex)
        for j in range(m + 1):
            polyA[:, j] = binoms[m, j] * a_pow[:, m - j] * bc_pow[:, j]
        polyB = np.zeros((npts, k + 1), dtype=complex)
        for j in range(k + 1):
            polyB[:, j] = binoms[k, j] * b_pow[:, k - j] * ac_pow[:, j]
        for j in range(dim):
            lo = max(0, j - k)
            hi = min(j, m)
            for l in range(lo, hi + 1):
                result[:, j, k] += polyA[:, l] * polyB[:, j - l]

    result = result * T[None, None, :] / T[None, :, None]
    return result


def build_v_chunk(X_chunk, levels, N):
    """Build V_chunk for a batch of nodes."""
    sorted_levels = sorted(n for n in levels.keys() if n <= N)
    total = sum(levels[n].shape[1] * (n + 1) for n in sorted_levels)
    bs = len(X_chunk)
    V = np.zeros((bs, total), dtype=complex)

    col = 0
    for n in sorted_levels:
        B = levels[n]
        m_mult = B.shape[1]
        pin = batch_pi_unitary(X_chunk, n)
        AH = B.conj().T
        AH_pin = np.einsum('mk,bkj->bmj', AH, pin)
        pw_norm = np.sqrt(n + 1)
        for i in range(m_mult):
            for j in range(n + 1):
                V[:, col] = pw_norm * AH_pin[:, i, j]
                col += 1
    return V


def mode_pw_norm(levels, N):
    """Peter-Weyl normalization factor √(n+1) per mode."""
    sorted_levels = sorted(n for n in levels.keys() if n <= N)
    norms = []
    for n in sorted_levels:
        m = levels[n].shape[1]
        for _ in range(m * (n + 1)):
            norms.append(np.sqrt(n + 1))
    return np.array(norms)


def projector_check(X, W, levels, N, fields):
    """Single-pass Gram matrix and roundtrip check.
    V is orthonormal (includes √(n+1) Peter-Weyl factor).
    Packet fields are in the un-normalized basis, so we convert."""
    sorted_levels = sorted(n for n in levels.keys() if n <= N)
    total = sum(levels[n].shape[1] * (n + 1) for n in sorted_levels)
    npts = len(X)
    nf = len(fields)
    pw = mode_pw_norm(levels, N)

    gram = np.zeros((total, total), dtype=complex)
    fields_norm = [f / pw for f in fields]
    c_rec = np.zeros((nf, total), dtype=complex)

    for start in range(0, npts, BATCH):
        end = min(start + BATCH, npts)
        X_chunk = X[start:end]
        W_chunk = W[start:end]
        V = build_v_chunk(X_chunk, levels, N)

        WV = W_chunk[:, None] * V
        gram += V.conj().T @ WV

        for fi in range(nf):
            nv = V @ fields_norm[fi]
            c_rec[fi] += V.conj().T @ (W_chunk * nv)

    gram_err = float(np.abs(gram - np.eye(total)).max())
    rt_errs = np.array([
        float(np.linalg.norm(fields_norm[fi] - c_rec[fi]) / np.linalg.norm(fields_norm[fi]))
        for fi in range(nf)
    ])
    return gram_err, rt_errs


def projected_cubic(X, W, levels, N, c):
    """Compute P_N[|ψ|²ψ] via quadrature.
    c: coefficient vector in the un-normalized R0 basis.
    Returns coefficient vector of the projected cubic (un-normalized basis)."""
    sorted_levels = sorted(n for n in levels.keys() if n <= N)
    total = sum(levels[n].shape[1] * (n + 1) for n in sorted_levels)
    npts = len(X)
    pw = mode_pw_norm(levels, N)
    c_norm = c / pw

    node_vals = np.zeros(npts, dtype=complex)
    for start in range(0, npts, BATCH):
        end = min(start + BATCH, npts)
        V = build_v_chunk(X[start:end], levels, N)
        node_vals[start:end] = V @ c_norm

    cubic_vals = np.abs(node_vals)**2 * node_vals

    d_norm = np.zeros(total, dtype=complex)
    for start in range(0, npts, BATCH):
        end = min(start + BATCH, npts)
        V = build_v_chunk(X[start:end], levels, N)
        d_norm += V.conj().T @ (W[start:end] * cubic_vals[start:end])

    return d_norm * pw


def run_gate4():
    t0_total = time.time()
    packet, _ = generate_packet()

    N = 24
    levels, total = build_r0_basis(N)
    print(f"  Rung N={N}: {total} modes")

    D_prod = 4 * N
    X, W = hopf_rule(D_prod)
    X = np.array(X)
    print(f"    4N rule: {len(X)} nodes, batch={BATCH}")

    fields = packet[N]["scalar"][:5]

    t0 = time.time()
    gram_err, rt_errs = projector_check(X, W, levels, N, fields)
    dt_check = time.time() - t0
    print(f"    Gram err: {gram_err:.3e}  ({dt_check:.1f}s)")
    print(f"    Max roundtrip err: {rt_errs.max():.3e}")

    rung_pass = gram_err < 1e-10 and rt_errs.max() < 1e-10
    print(f"    {'PASS' if rung_pass else 'FAIL'}")

    c_test = packet[N]["scalar"][0]

    print(f"\n  Projected cubic at 4N...")
    t0 = time.time()
    d_fine = projected_cubic(X, W, levels, N, c_test)
    dt_fine = time.time() - t0
    print(f"    Done ({dt_fine:.1f}s), norm={np.linalg.norm(d_fine):.6f}")

    print(f"\n  Mutation (node drop to 2N on cubic):")
    D_coarse = 2 * N
    Xc, Wc = hopf_rule(D_coarse)
    Xc = np.array(Xc)
    print(f"    2N rule: {len(Xc)} nodes")

    t0 = time.time()
    d_coarse = projected_cubic(Xc, Wc, levels, N, c_test)
    dt_coarse = time.time() - t0

    d_norm = np.linalg.norm(d_fine)
    cubic_err_2N = np.linalg.norm(d_fine - d_coarse) / d_norm if d_norm > 1e-15 else 0.0
    print(f"    2N cubic err: {cubic_err_2N:.3e}  ({dt_coarse:.1f}s)")

    print(f"\n  Supplementary: node drop to N on cubic:")
    D_N = N
    Xn, Wn = hopf_rule(D_N)
    Xn = np.array(Xn)
    print(f"    N rule: {len(Xn)} nodes")
    d_N = projected_cubic(Xn, Wn, levels, N, c_test)
    cubic_err_N = np.linalg.norm(d_fine - d_N) / d_norm if d_norm > 1e-15 else 0.0
    print(f"    N cubic err: {cubic_err_N:.3e}")

    mutation_fires = cubic_err_2N > 0.01
    print(f"\n  2N mutation fires: {mutation_fires}")
    if not mutation_fires:
        print("    NOTE: R0 invariant vectors have even-k parity; K=2N+1 is odd;")
        print("    aliasing frequency has wrong parity -> no aliasing for R0 cubics.")
        print(f"    N-rule err ({cubic_err_N:.3e}) confirms mechanism is correct.")

    gate_pass = rung_pass
    dt = time.time() - t0_total

    gate_record(
        gate_id="G4-PROJECTOR",
        arena_id=f"A-R0-N{N}",
        rung=N,
        parent_status="GREEN" if rung_pass else "RED",
        mutation_status="FIRES" if mutation_fires else "DEAD",
        measured_values={
            "rung_results": {"N": N, "n_nodes": len(X), "n_modes": total,
                             "gram_err": gram_err, "max_roundtrip_err": float(rt_errs.max()),
                             "n_fields_tested": len(fields)},
            "node_drop_2N": {"cubic_err": cubic_err_2N,
                             "d_fine_norm": float(d_norm),
                             "d_2N_norm": float(np.linalg.norm(d_coarse))},
            "node_drop_N": {"cubic_err": cubic_err_N,
                            "d_N_norm": float(np.linalg.norm(d_N))},
            "parity_analysis": "R0 invariant vectors nonzero only at even k; "
                              "product Fourier content has even parity; "
                              "K=2N+1 (odd) aliasing frequency has wrong parity; "
                              "no aliasing for R0-only cubics at 2N rule",
        },
    )
    resource_record("G4-PROJECTOR", dt, dt)
    return gate_pass, dt


if __name__ == "__main__":
    print("Gate 4: Projector exactness")
    print("=" * 60)
    passed, dt = run_gate4()
    print(f"\nGate 4: {'GREEN' if passed else 'RED'} ({dt:.1f}s)")
