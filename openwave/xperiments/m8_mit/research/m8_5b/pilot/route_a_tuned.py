#!/usr/bin/env python3
"""
M8.5-B engineering pilot, route (a) prototype 6: tuning for a converging ladder.

NON-EVIDENTIARY pilot code (frozen pre-registration section 6). Deck group 2T,
order 24, nonabelian and NOT the target group. No 2I-specific quantity here.

THE QUESTION
  Prototype 5 reached the required band (first invariant level near 48) but the error
  there plateaued near 1 percent and did not improve with node count. The frozen
  contract requires a converging resolution ladder and a convergence statistic, so a
  plateau is not good enough. This asks whether stencil order, stencil size, and node
  quality can be chosen to turn that plateau into a demonstrably converging sequence.

WHAT PROTOTYPE 5 GOT WRONG
  For polyharmonic-spline RBF-FD the convergence order is governed by the POLYNOMIAL
  augmentation degree, not by the spline exponent. Prototype 5 used phi = u^5 with
  degree-1 augmentation, so a high-order kernel was paired with a nearly order-free
  polynomial space. That is the plateau's most likely cause, and it is a tuning error
  rather than a limit of the method.

WHAT IS SWEPT HERE
  p : polynomial augmentation degree, with the exact surface Laplacian of each
      monomial, Delta_{S^3} x^a = sum_i a_i(a_i-1) x^(a-2e_i) - n(n+2) x^a  on |x|=1
  k : stencil size, held at roughly twice the polynomial space dimension
  node quality : seeds relaxed by Riesz repulsion BEFORE orbits are generated, so the
      cloud stays exactly Gamma-closed while becoming quasi-uniform

  Monomials of degree <= p are linearly dependent on the sphere, so the local saddle
  system is solved in least-squares form rather than by direct inversion.
"""

import argparse
import itertools
import json
import time

import numpy as np
from scipy.spatial import cKDTree

from route_a_nonabelian import close_group, qmul


def monomials(p):
    """Exponent tuples of total degree <= p in 4 variables."""
    out = []
    for total in range(p + 1):
        for a in itertools.product(range(total + 1), repeat=4):
            if sum(a) == total:
                out.append(a)
    return out


def eval_monos(X, monos):
    M = np.ones((len(X), len(monos)))
    for j, a in enumerate(monos):
        for i, e in enumerate(a):
            if e:
                M[:, j] *= X[:, i] ** e
    return M


def surf_lap_monos(x, monos):
    """Delta_{S^3} of each monomial at a single point x on the unit sphere."""
    idx = {a: j for j, a in enumerate(monos)}
    out = np.zeros(len(monos))
    for j, a in enumerate(monos):
        n = sum(a)
        val = 1.0
        for i, e in enumerate(a):
            if e:
                val *= x[i] ** e
        out[j] -= n * (n + 2) * val
        for i, e in enumerate(a):
            if e >= 2:
                b = list(a)
                b[i] -= 2
                b = tuple(b)
                if b in idx:
                    v = 1.0
                    for ii, ee in enumerate(b):
                        if ee:
                            v *= x[ii] ** ee
                    out[idx[b]] += 0.0      # placeholder, handled below
                    out[j] += 0.0
                # accumulate the ambient term directly into the value vector
    # ambient Laplacian term, computed as a separate pass into monomial values
    amb = np.zeros(len(monos))
    for j, a in enumerate(monos):
        for i, e in enumerate(a):
            if e >= 2:
                b = list(a); b[i] -= 2; b = tuple(b)
                v = e * (e - 1)
                for ii, ee in enumerate(b):
                    if ee:
                        v *= x[ii] ** ee
                amb[j] += v
    return out + amb


def relax_seeds(seeds, G, n_iter, step=0.05):
    """Riesz repulsion on the FULL orbit cloud, moving only the seeds."""
    S = seeds.copy()
    for _ in range(n_iter):
        cloud, owner = [], []
        for si, v in enumerate(S):
            for g in G:
                cloud.append(qmul(g, v)); owner.append(si)
        cloud = np.array(cloud); owner = np.array(owner)
        tree = cKDTree(cloud)
        d, idx = tree.query(cloud, k=8)
        F = np.zeros_like(S)
        for i in range(len(cloud)):
            if owner[i] != -1:
                diff = cloud[i] - cloud[idx[i, 1:]]
                w = 1.0 / np.maximum(d[i, 1:], 1e-6) ** 3
                F[owner[i]] += (diff * w[:, None]).sum(axis=0)
        S = S + step * F / (np.linalg.norm(F, axis=1, keepdims=True) + 1e-12)
        S /= np.linalg.norm(S, axis=1, keepdims=True)
    return S


def gamma_cloud_from_seeds(S, G):
    nodes, oid = [], []
    for si, v in enumerate(S):
        orbit = []
        for g in G:
            w = qmul(g, v)
            if not any(np.linalg.norm(w - u) < 1e-9 for u in orbit):
                orbit.append(w)
        for w in orbit:
            nodes.append(w); oid.append(si)
    return np.array(nodes), np.array(oid), len(S)


def build_L(X, k, m, p):
    monos = monomials(p)
    npoly = len(monos)
    N = len(X)
    tree = cKDTree(X)
    _, idx = tree.query(X, k=k)
    L = np.zeros((N, N))
    for i in range(N):
        s = idx[i]
        P = X[s]
        d = np.linalg.norm(P[:, None, :] - P[None, :, :], axis=-1)
        A = d ** m
        Q = eval_monos(P, monos)
        Mx = np.zeros((k + npoly, k + npoly))
        Mx[:k, :k] = A
        Mx[:k, k:] = Q
        Mx[k:, :k] = Q.T
        rhs = np.zeros(k + npoly)
        u0 = np.linalg.norm(P - X[i], axis=1)
        rhs[:k] = m * (m + 1) * np.where(u0 > 0, u0 ** (m - 2), 0.0) \
            - (m * (m + 4) / 4.0) * u0 ** m
        rhs[k:] = surf_lap_monos(X[i], monos)
        w = np.linalg.lstsq(Mx, rhs, rcond=None)[0][:k]
        L[i, s] = w
    return L


def reduced_spectrum(L, oid, M, n_ev=4):
    reps = [np.where(oid == mm)[0][0] for mm in range(M)]
    A = -L
    R = np.zeros((M, M))
    for mm, i in enumerate(reps):
        row = A[i]
        for l in range(M):
            R[mm, l] = row[oid == l].sum()
    ev = np.linalg.eigvals(R)
    return np.sort(ev.real)[:n_ev]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, nargs="+", default=[60, 120, 240])
    ap.add_argument("--configs", type=str, default="2:30:5,3:60:5,3:70:7")
    ap.add_argument("--relax", type=int, default=40)
    ap.add_argument("--json")
    args = ap.parse_args()

    G = close_group([np.array([0.0, 1, 0, 0]), np.array([0.5, 0.5, 0.5, 0.5])])
    TARGET = 48.0
    cfgs = []
    for c in args.configs.split(","):
        p, k, m = (int(v) for v in c.split(":"))
        cfgs.append((p, k, m))

    print(f"deck group 2T, order {len(G)}; target first invariant level = {TARGET}")
    print(f"seed relaxation iterations: {args.relax}")
    print(f"configs (p=poly degree, k=stencil, m=PHS): {cfgs}")
    print()

    out = {}
    for (p, k, m) in cfgs:
        npoly = len(monomials(p))
        print(f"  p={p} (poly dim {npoly})  k={k}  m={m}")
        errs, rows = [], []
        for ns in args.seeds:
            rng = np.random.default_rng(20260731)
            S = rng.normal(size=(ns, 4)); S /= np.linalg.norm(S, axis=1, keepdims=True)
            S = relax_seeds(S, G, args.relax)
            X, oid, M = gamma_cloud_from_seeds(S, G)
            t0 = time.time()
            L = build_L(X, k, m, p)
            ev = reduced_spectrum(L, oid, M)
            dt = time.time() - t0
            lam1 = ev[1]
            err = abs(lam1 - TARGET) / TARGET
            errs.append(err)
            rows.append({"nodes": len(X), "orbits": M, "lambda1": float(lam1),
                         "rel_err": float(err), "seconds": round(dt, 2)})
            print(f"    nodes={len(X):>5}  lambda1={lam1:9.4f}  relerr={err:8.2e}  {dt:5.1f}s")
        if len(errs) >= 2:
            orders = [np.log2(a / b) for a, b in zip(errs[:-1], errs[1:])]
            print(f"    observed order between rungs: "
                  + ", ".join(f"{o:.2f}" for o in orders))
            converging = all(b < a for a, b in zip(errs[:-1], errs[1:]))
            print(f"    monotone decreasing: {converging}")
            out[f"p{p}_k{k}_m{m}"] = {"rows": rows, "orders": orders,
                                      "monotone": bool(converging)}
        print()

    if args.json:
        with open(args.json, "w") as fh:
            json.dump({"group": "2T", "target": TARGET,
                       "relax_iters": args.relax,
                       "note": "pilot, non-evidentiary; not the target group",
                       "configs": out}, fh, indent=2, sort_keys=True)
            fh.write("\n")
        print(f"wrote {args.json}")


if __name__ == "__main__":
    main()
