#!/usr/bin/env python3
"""
M8.5-B engineering pilot, route (a) prototype 8: Hodge Laplacian via Weitzenbock.

NON-EVIDENTIARY pilot code (frozen pre-registration section 6). Deck group 2T,
order 24, nonabelian and NOT the target group. No 2I-specific quantity here.

WHY THIS REPLACES THE CURL BUILD
  Prototype 7 built curl directly. Its operators were verified correct and it recovered
  the exactly-known lambda = 4 level to 1e-13, but only 18 to 21 percent of the reduced
  spectrum was effectively real and that fraction did not improve with resolution. A
  first-order operator on a scattered cloud is not dissipative and needs stabilization.
  The Hodge Laplacian is second order, like the scalar operator that behaved well
  throughout, so this builds it directly and never forms curl.

THE DERIVATION, DONE HERE RATHER THAN LOOKED UP
  Weitzenbock on one-forms, from the pinned reference sheet:

      Delta_H  =  nabla* nabla  +  Ric,      Ric = 2 g on the unit S^3

  For a bi-invariant metric the Levi-Civita connection on left-invariant fields is
  nabla_{e_i} e_j = (1/2)[e_i, e_j] = eps_ijk e_k, so for omega = f_j sigma^j

      (nabla_i omega)_j = e_i(f_j) - eps_ijk f_k

  Tracing the second covariant derivative, and using the two contractions
  sum_{i,a} eps_ija eps_iab = -2 delta_jb  and  sum_{i,k} eps_ijk eps_ikl = -2 delta_jl
  (both verified by enumeration before use), gives

      (Delta_H omega)_j  =  L_pos f_j  +  2 eps_ijk e_i(f_k)  +  4 f_j

  with L_pos the POSITIVE scalar Laplacian. The dominant term is the second-order scalar
  operator; the coupling is first order and subordinate, which is the structural reason
  to expect this to behave where curl did not.

  Two independent consequences are checked before any spectrum is read:
    constant components (the left-invariant one-forms) must give exactly 4;
    omega = df must give exactly the scalar eigenvalue n(n+2).
"""

import argparse
import itertools
import json
import time

import numpy as np

from route_a_tuned import relax_seeds, gamma_cloud_from_seeds, build_L
from route_a_oneform import left_frame, build_dirderiv
from route_a_nonabelian import close_group

EPS = np.zeros((3, 3, 3))
for _i, _j, _k in itertools.permutations(range(3)):
    EPS[_i, _j, _k] = np.sign(np.linalg.det(np.eye(3)[[_i, _j, _k]]))


def hodge_matrix(Lpos, D, N):
    """(Delta_H)_j = L_pos f_j + 2 eps_ijk D_i f_k + 4 f_j, as 3N x 3N."""
    H = np.zeros((3 * N, 3 * N))
    for j in range(3):
        H[j*N:(j+1)*N, j*N:(j+1)*N] += Lpos + 4.0 * np.eye(N)
        for i in range(3):
            for k in range(3):
                if EPS[i, j, k]:
                    H[j*N:(j+1)*N, k*N:(k+1)*N] += 2.0 * EPS[i, j, k] * D[i]
    return H


def orbit_reduce_block(A, oid, M, N):
    reps = [np.where(oid == mm)[0][0] for mm in range(M)]
    R = np.zeros((3 * M, 3 * M))
    for a in range(3):
        for b in range(3):
            blk = A[a*N:(a+1)*N, b*N:(b+1)*N]
            for mm, i in enumerate(reps):
                row = blk[i]
                for l in range(M):
                    R[a*M + mm, b*M + l] = row[oid == l].sum()
    return R


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, nargs="+", default=[60, 120])
    ap.add_argument("--k", type=int, default=60)
    ap.add_argument("--m", type=int, default=5)
    ap.add_argument("--p", type=int, default=3)
    ap.add_argument("--json")
    args = ap.parse_args()

    G = close_group([np.array([0.0, 1, 0, 0]), np.array([0.5, 0.5, 0.5, 0.5])])
    rows = []
    for ns in args.seeds:
        rng = np.random.default_rng(20260731)
        S = rng.normal(size=(ns, 4)); S /= np.linalg.norm(S, axis=1, keepdims=True)
        S = relax_seeds(S, G, 40)
        X, oid, M = gamma_cloud_from_seeds(S, G)
        N = len(X)

        t0 = time.time()
        Lpos = -build_L(X, args.k, args.m, args.p)
        D = build_dirderiv(X, left_frame(X), args.k, args.m, args.p)
        H = hodge_matrix(Lpos, D, N)
        t_build = time.time() - t0

        # check 1: left-invariant one-forms must give exactly 4
        v = np.zeros(3 * N); v[0:N] = 1.0
        c1 = float(np.abs((H @ v)[0:N] - 4.0).max())

        # check 2: omega = df must give the scalar eigenvalue.  Use a degree-1 harmonic,
        # whose scalar eigenvalue is 3, so Delta_H (df) must equal 3 df.
        a = np.array([0.3, -0.5, 0.7, 0.2]); a /= np.linalg.norm(a)
        f = X @ a
        df = np.concatenate([D[j] @ f for j in range(3)])
        c2 = float(np.abs(H @ df - 3.0 * df).max() / np.abs(df).max())

        Rh = orbit_reduce_block(H, oid, M, N)
        ev = np.linalg.eigvals(Rh)
        real_frac = float((np.abs(ev.imag) < 1e-6).mean())
        r = np.sort(ev.real)

        print(f"  nodes={N:>5} orbits={M:>4}  build {t_build:5.1f}s")
        print(f"    check 1, left-invariant forms -> 4 : max err {c1:.2e}")
        print(f"    check 2, Delta_H(df) = 3 df        : rel err {c2:.2e}")
        print(f"    effectively real modes             : {100*real_frac:5.1f}%  "
              f"(curl build gave 18-21%)")
        print(f"    lowest eigenvalues                 : "
              f"{np.array2string(r[:8], precision=4)}")
        rows.append({"nodes": N, "orbits": M, "t_build": round(t_build, 2),
                     "check_leftinv": c1, "check_exact": c2,
                     "real_fraction": real_frac,
                     "lowest": [float(x) for x in r[:8]]})
        print()

    if args.json:
        with open(args.json, "w") as fh:
            json.dump({"group": "2T", "note": "pilot, non-evidentiary",
                       "operator": "Hodge via Weitzenbock, no curl formed",
                       "rows": rows}, fh, indent=2, sort_keys=True)
            fh.write("\n")
        print(f"wrote {args.json}")


if __name__ == "__main__":
    main()
