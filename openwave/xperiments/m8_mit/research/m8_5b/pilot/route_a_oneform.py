#!/usr/bin/env python3
"""
M8.5-B engineering pilot, route (a) prototype 7: one-forms on the grid.

NON-EVIDENTIARY pilot code (frozen pre-registration section 6). Deck group 2T,
order 24, nonabelian and NOT the target group. No 2I-specific quantity here.

THE REPRESENTATION, AND WHY IT MAKES THE QUOTIENT EASY
  A one-form is carried by its three components in the LEFT-invariant coframe,
  omega = f_1 sigma^1 + f_2 sigma^2 + f_3 sigma^3, with the left-invariant fields
  e_j(q) = q * u_j for u_j in {i, j, k} under quaternion multiplication.

  A deck group acting by LEFT multiplication moves the base point but leaves the frame
  alone, so it acts on the three coefficient functions independently. The orbit
  reduction that worked for scalars therefore applies COMPONENTWISE, with no new
  identification machinery. That is the payoff of the frame choice.

THE OPERATOR
  With the unit-curvature normalization the left-invariant fields satisfy
  [e_i, e_j] = 2 eps_ijk e_k, and curl acts on components as

      (curl omega)_i = eps_ijk e_j(f_k) - 2 f_i

  On coexact one-forms the Hodge Laplacian is curl^2, and curl annihilates exact forms.
  So the two sectors separate without ever forming the Hodge Laplacian directly: the
  kernel of curl is the exact sector, and the nonzero curl eigenvalues give the coexact
  spectrum as lambda = (curl eigenvalue)^2.

  The constant in the curl is CHECKED before use: constant components represent the
  left-invariant one-forms themselves, and the formula must return exactly -2 on them.

  Directional derivatives are built by RBF-FD on the same Gamma-closed cloud, with the
  polyharmonic derivative d/dv |x-y|^m = m |x-y|^(m-2) (x-y).v.
"""

import argparse
import itertools
import json
import time

import numpy as np
from scipy.spatial import cKDTree

from route_a_tuned import relax_seeds, gamma_cloud_from_seeds, monomials, eval_monos
from route_a_nonabelian import close_group, qmul

IMAG = [np.array([0.0, 1, 0, 0]), np.array([0.0, 0, 1, 0]), np.array([0.0, 0, 0, 1])]


def left_frame(X):
    """e_j(q) = q * u_j, the left-invariant orthonormal frame, as tangent vectors."""
    return [np.array([qmul(q, u) for q in X]) for u in IMAG]


def build_dirderiv(X, frame, k, m, p):
    """RBF-FD directional derivatives along each left-invariant field."""
    monos = monomials(p)
    npoly = len(monos)
    N = len(X)
    tree = cKDTree(X)
    _, idx = tree.query(X, k=k)
    D = [np.zeros((N, N)) for _ in range(3)]
    for i in range(N):
        s = idx[i]
        P = X[s]
        d = np.linalg.norm(P[:, None, :] - P[None, :, :], axis=-1)
        A = d ** m
        Q = eval_monos(P, monos)
        M = np.zeros((k + npoly, k + npoly))
        M[:k, :k] = A
        M[:k, k:] = Q
        M[k:, :k] = Q.T
        diff = X[i] - P                       # (k,4)
        u0 = np.linalg.norm(diff, axis=1)
        for j in range(3):
            v = frame[j][i]
            rhs = np.zeros(k + npoly)
            with np.errstate(divide="ignore", invalid="ignore"):
                coef = m * np.where(u0 > 0, u0 ** (m - 2), 0.0)
            rhs[:k] = coef * (diff @ v)
            # directional derivative of each monomial at X[i]
            for jj, a in enumerate(monos):
                g = 0.0
                for c, e in enumerate(a):
                    if e:
                        t = e * X[i][c] ** (e - 1)
                        for cc, ee in enumerate(a):
                            if cc != c and ee:
                                t *= X[i][cc] ** ee
                        g += t * v[c]
                rhs[k + jj] = g
            w = np.linalg.lstsq(M, rhs, rcond=None)[0][:k]
            D[j][i, s] = w
    return D


def curl_matrix(D, N):
    """(curl omega)_i = eps_ijk D_j f_k - 2 f_i, as a 3N x 3N operator."""
    C = np.zeros((3 * N, 3 * N))
    eps = np.zeros((3, 3, 3))
    for i, j, k in itertools.permutations(range(3)):
        eps[i, j, k] = np.sign(np.linalg.det(np.eye(3)[[i, j, k]]))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                if eps[i, j, k]:
                    C[i*N:(i+1)*N, k*N:(k+1)*N] += eps[i, j, k] * D[j]
        C[i*N:(i+1)*N, i*N:(i+1)*N] -= 2.0 * np.eye(N)
    return C


def orbit_reduce_block(C, oid, M, N):
    """Componentwise orbit reduction: valid because the frame is left-invariant."""
    reps = [np.where(oid == mm)[0][0] for mm in range(M)]
    R = np.zeros((3 * M, 3 * M))
    for a in range(3):
        for b in range(3):
            blk = C[a*N:(a+1)*N, b*N:(b+1)*N]
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
        frame = left_frame(X)

        # frame sanity: orthonormal and tangent, checked not assumed
        tang = max(abs(float(frame[j][i] @ X[i])) for j in range(3) for i in range(0, N, 97))
        norm = max(abs(np.linalg.norm(frame[j][i]) - 1) for j in range(3) for i in range(0, N, 97))

        t0 = time.time()
        D = build_dirderiv(X, frame, args.k, args.m, args.p)
        C = curl_matrix(D, N)
        t_build = time.time() - t0

        # constant-component check: curl of a left-invariant one-form must be exactly -2
        v = np.zeros(3 * N); v[0:N] = 1.0
        got = C @ v
        const_err = float(np.abs(got[0:N] + 2.0).max())

        Rc = orbit_reduce_block(C, oid, M, N)
        ev = np.linalg.eigvals(Rc)
        ev = ev[np.argsort(np.abs(ev.real))]
        nz = ev[np.abs(ev.real) > 0.5]
        lam = np.sort((nz.real ** 2))

        print(f"  nodes={N:>5} orbits={M:>4}  build {t_build:5.1f}s")
        print(f"    frame tangency max |e.q| = {tang:.1e}, unit-norm dev = {norm:.1e}")
        print(f"    curl on left-invariant forms: max |curl + 2| = {const_err:.2e}")
        print(f"    lowest coexact lambda = curl^2: "
              f"{np.array2string(lam[:8], precision=3)}")
        rows.append({"nodes": N, "orbits": M, "t_build": round(t_build, 2),
                     "frame_tangency": tang, "const_check": const_err,
                     "lowest_coexact": [float(x) for x in lam[:8]]})
        print()

    if args.json:
        with open(args.json, "w") as fh:
            json.dump({"group": "2T", "note": "pilot, non-evidentiary",
                       "rows": rows}, fh, indent=2, sort_keys=True)
            fh.write("\n")
        print(f"wrote {args.json}")


if __name__ == "__main__":
    main()
