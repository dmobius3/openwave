#!/usr/bin/env python3
"""Solver-versus-discretization localization test for the M8.4 P1A closeout.

CLAIM THIS SUPPORTS
    The imaginary parts in the computed target clusters are a property of the assembled
    discrete matrix L_h, not of the eigensolver's floating-point arithmetic.

METHOD
    Hold L_h byte-for-byte fixed.  Change ONLY the arithmetic precision of the eigensolve.
    If the imaginary parts collapse as precision rises, they are solver arithmetic.  If they
    are unchanged, they belong to the matrix.

WHAT IT DOES NOT ESTABLISH
    Which upstream stage is responsible.  This separates solver from matrix.  It does NOT
    separate the equivariant assembly from the underlying RBF-FD discretization.

DISCRIMINATING CONTROL
    A control matrix whose exact spectrum is real, but whose computed imaginary parts are
    numerical error from the finite-precision eigensolve, is run through the same ladder.
    Its |Im| must fall with precision.  Without it the test cannot fail, and a ladder that
    moves nothing would be indistinguishable from one that is broken.

USAGE
    PYTHONPATH=<dir containing p0>:. python3 localization_test.py
"""
import sys
import numpy as np
import mpmath as mp

from p0 import group as G, representations as RP
from p0.cloud import fibonacci_seeds_s3, build_orbit_cloud
from p0.bundle_operator import build_L_bundle

N_SEEDS = 20                      # 40x40 to 60x60 blocks: tractable at 50 digits
K_STENCIL, RBF_M, RBF_P = 110, 7, 4
DPS_LADDER = (30, 50)
MCKAY_DIST = {"R1": 1, "R3": 2, "R2": 7}
COLLAPSE_FACTOR = 1e3             # a solver-limited quantity must fall by at least this


def cluster_max_imag(evals, lam0, k):
    idx = np.argsort(np.abs(evals - lam0))[:k]
    return float(np.max(np.abs(np.imag(evals[idx]))))


def mp_eigvals(L, dps):
    mp.mp.dps = dps
    n = L.shape[0]
    A = mp.matrix([[mp.mpc(float(L[i, j].real), float(L[i, j].imag))
                    for j in range(n)] for i in range(n)])
    E, _ = mp.eig(A)
    return np.array([complex(x) for x in E])


def main():
    elems = G.build_icosians()
    chi = G.build_character_table(elems)
    out = RP.build_all_representations(elems, chi)
    reps = out[0] if isinstance(out, tuple) else out
    seeds = fibonacci_seeds_s3(N_SEEDS)
    X, oid, gid = build_orbit_cloud(seeds, elems)

    print(f"cloud: {N_SEEDS} seeds, {len(X)} nodes, k={K_STENCIL}")
    print(f"precision ladder: float64 (~16 dps), then {', '.join(str(d) for d in DPS_LADDER)}\n")

    rows, persist = [], True
    for lab, d in MCKAY_DIST.items():
        k, lam0 = d + 1, d * (d + 2)
        L = np.asarray(build_L_bundle(X, oid, gid, elems, reps[lab],
                                      k=K_STENCIL, m=RBF_M, p=RBF_P)[0])
        vals = [cluster_max_imag(np.linalg.eigvals(L), lam0, k)]
        for dps in DPS_LADDER:
            vals.append(cluster_max_imag(mp_eigvals(L, dps), lam0, k))
        rows.append((lab, L.shape[0], lam0, vals))
        if min(vals) < max(vals) / COLLAPSE_FACTOR:
            persist = False
        print(f"  {lab}  N={L.shape[0]:<4} lambda={lam0:<3} "
              + "  ".join(f"{v:.4e}" for v in vals))

    # MUTATION ARM: the ladder must be able to shrink a quantity that IS solver-limited,
    # or a persistence result is indistinguishable from a broken ladder.
    #
    # A real SYMMETRIC matrix will not serve: LAPACK returns exactly zero imaginary parts,
    # so nothing can collapse and the arm is vacuous.  Nor will a well-separated real
    # spectrum under an ill-conditioned similarity: that also returns exactly real values.
    # What does serve is a highly MULTIPLE real root, where roundoff genuinely scatters the
    # cluster into the complex plane: the companion matrix of (x-1)^m has spectrum exactly
    # {1} with multiplicity m, so every computed imaginary part is pure arithmetic.
    m = 10
    c = np.poly(np.ones(m))
    C = np.zeros((m, m))
    C[0, :] = -c[1:] / c[0]
    C[1:, :-1] = np.eye(m - 1)
    ctrl = [float(np.max(np.abs(np.imag(np.linalg.eigvals(C)))))]
    for dps in DPS_LADDER:
        ctrl.append(float(np.max(np.abs(np.imag(mp_eigvals(C.astype(complex), dps))))))
    collapsed = ctrl[-1] < ctrl[0] / COLLAPSE_FACTOR
    print(f"\n  CONTROL, companion of (x-1)^{m}, true |Im| = 0 exactly:")
    print("    " + "  ".join(f"{v:.4e}" for v in ctrl)
          + f"   collapse factor {ctrl[0] / max(ctrl[-1], 1e-300):.2e}")
    print(f"    ladder resolves a solver-limited quantity: {collapsed}")

    print()
    if not collapsed:
        print("MUTATION ARM FAILED: the ladder does not detect a solver-limited quantity.")
        print("The persistence result below is NOT trustworthy.")
        return 1
    if persist:
        print("RESULT: target imaginary parts PERSIST across the precision ladder.")
        print("They belong to the assembled matrix L_h, not to the eigensolver.")
        print("This does NOT identify assembly versus discretization as the source.")
        return 0
    print("RESULT: target imaginary parts COLLAPSE with precision: solver arithmetic.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
