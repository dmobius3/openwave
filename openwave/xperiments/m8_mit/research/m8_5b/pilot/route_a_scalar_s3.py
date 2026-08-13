#!/usr/bin/env python3
"""
M8.5-B engineering pilot, route (a) prototype 1: scalar Laplacian on the unit S^3.

NON-EVIDENTIARY. This is pilot code under the frozen pre-registration's section 6:
its only purpose is to select [PILOT] values. It produces no certification result,
touches neither sealed packet, and examines no 2I-specific quantity.

SCHEME UNDER TEST
  Hopf (toroidal) coordinates on S^3:
      (eta, xi1, xi2)  ->  (cos(eta) e^{i xi1}, sin(eta) e^{i xi2})
      eta in [0, pi/2],  xi1, xi2 in [0, 2pi)
      ds^2 = d eta^2 + cos^2(eta) d xi1^2 + sin^2(eta) d xi2^2
      sqrt(g) = cos(eta) sin(eta)

  Both xi are periodic, so Fourier-diagonalize them exactly:
      f = e^{i m1 xi1} e^{i m2 xi2} g(eta)
  leaving, for each (m1, m2), a 1D Sturm-Liouville problem in eta:

      -(1/w) d/d eta ( w  dg/d eta ) + ( m1^2/cos^2 eta + m2^2/sin^2 eta ) g = lambda g
      w(eta) = cos(eta) sin(eta)

  Discretized with a staggered (cell-centred) finite-difference grid, which keeps
  the sample points off the coordinate-singular endpoints where w -> 0.

  This is the "intrinsic charts" half of the route (a) menu. The eta direction is
  the only one carrying a grid; the group action will later act on the two periodic
  fibres, which is why this coordinate choice is worth measuring first.

KNOWN ANSWER (rung 1, unit radius, pinned in the reference sheet)
  eigenvalues  n(n+2),  n >= 0  ->  0, 3, 8, 15, 24, 35, ...
  multiplicity (n+1)^2          ->  1, 4, 9, 16, 25, 36, ...

WHAT IS MEASURED
  accuracy   : max relative error on the resolved levels
  cost       : wall time, and the size of the largest matrix factorized
  degeneracy : the spread within each eigenvalue cluster, which is the quantity
               the frozen G7-analogue for route (a) will have to report
"""

import argparse
import json
import time

import numpy as np
from scipy.linalg import eigh_tridiagonal


def radial_spectrum(m1, m2, N, n_want):
    """
    Eigenvalues of the 1D operator for one (m1, m2) Fourier sector.

    Cell-centred grid on (0, pi/2): eta_j = (j + 1/2) h, j = 0..N-1, h = (pi/2)/N.
    The second-order term is discretized in self-adjoint form using w at the cell
    faces, which keeps the matrix symmetric tridiagonal.
    """
    h = (np.pi / 2) / N
    j = np.arange(N)
    eta = (j + 0.5) * h
    w = np.cos(eta) * np.sin(eta)                    # weight at centres
    eta_f = (j[:-1] + 1.0) * h                       # interior faces
    w_f = np.cos(eta_f) * np.sin(eta_f)              # weight at faces

    # -(1/w) d/deta ( w dg/deta )  ->  symmetric tridiagonal after scaling by 1/w
    off = -w_f / (h * h * np.sqrt(w[:-1] * w[1:]))
    diag = np.zeros(N)
    diag[:-1] += w_f / (h * h * w[:-1])
    diag[1:] += w_f / (h * h * w[1:])

    # centrifugal terms
    diag += (m1 * m1) / np.cos(eta) ** 2 + (m2 * m2) / np.sin(eta) ** 2

    k = min(n_want, N)
    vals = eigh_tridiagonal(diag, off, select="i", select_range=(0, k - 1),
                            eigvals_only=True)
    return vals


def run(N, m_max, lam_max):
    """Assemble the full scalar spectrum up to lam_max by sweeping Fourier sectors."""
    t0 = time.time()
    found = []
    sectors = 0
    for m1 in range(-m_max, m_max + 1):
        for m2 in range(-m_max, m_max + 1):
            # a sector cannot contribute below its centrifugal floor
            floor = m1 * m1 + m2 * m2
            if floor > lam_max + 1:
                continue
            sectors += 1
            vals = radial_spectrum(m1, m2, N, n_want=max(4, int(np.sqrt(lam_max)) + 2))
            found.extend(v for v in vals if v <= lam_max + 1e-6)
    elapsed = time.time() - t0
    return np.sort(np.array(found)), elapsed, sectors


def cluster(vals, tol):
    """Group eigenvalues into degeneracy clusters with an absolute tolerance."""
    out = []
    for v in vals:
        if out and abs(v - out[-1][-1]) <= tol:
            out[-1].append(v)
        else:
            out.append([v])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ladder", type=int, nargs="+", default=[40, 80, 160, 320])
    ap.add_argument("--lam-max", type=float, default=25.0)
    ap.add_argument("--m-max", type=int, default=6)
    ap.add_argument("--cluster-tol", type=float, default=0.35)
    ap.add_argument("--json")
    args = ap.parse_args()

    # pinned rung-1 reference, unit radius
    ref = [(n, n * (n + 2), (n + 1) ** 2) for n in range(0, 12)]
    ref = [(n, lam, mult) for (n, lam, mult) in ref if lam <= args.lam_max]

    print("route (a) prototype 1: Hopf coordinates, Fourier x staggered FD")
    print(f"target levels through lambda <= {args.lam_max}: "
          + ", ".join(f"n={n}:{lam}(x{mult})" for n, lam, mult in ref))
    print()

    rows = []
    for N in args.ladder:
        vals, elapsed, sectors = run(N, args.m_max, args.lam_max)
        clusters = cluster(vals, args.cluster_tol)
        got = [(float(np.mean(c)), len(c), float(np.ptp(c))) for c in clusters]

        print(f"N = {N:>4}   sectors={sectors:>4}   {elapsed:6.2f}s   "
              f"clusters found: {len(got)}")
        worst_rel = 0.0
        mult_ok = True
        for (n, lam, mult), (gv, gm, spread) in zip(ref, got):
            rel = abs(gv - lam) / lam if lam else abs(gv)
            worst_rel = max(worst_rel, rel)
            flag = "" if gm == mult else f"  <-- multiplicity {gm}, expected {mult}"
            mult_ok &= (gm == mult)
            print(f"    n={n}  lambda={lam:<3} got={gv:12.8f}  relerr={rel:9.2e}  "
                  f"mult={gm:<3} spread={spread:9.2e}{flag}")
        rows.append({
            "N": N, "seconds": round(elapsed, 3), "sectors": sectors,
            "clusters": len(got), "worst_rel_err": worst_rel,
            "multiplicities_correct": bool(mult_ok),
            "max_cluster_spread": max((s for _, _, s in got), default=0.0),
        })
        print()

    print("convergence summary")
    print(f"  {'N':>5} {'sec':>8} {'worst rel err':>15} {'max spread':>12}  mults")
    for r in rows:
        print(f"  {r['N']:>5} {r['seconds']:>8.2f} {r['worst_rel_err']:>15.3e} "
              f"{r['max_cluster_spread']:>12.3e}  {'ok' if r['multiplicities_correct'] else 'WRONG'}")

    if len(rows) >= 2:
        print("\n  observed order (successive halving of h):")
        for a, b in zip(rows[:-1], rows[1:]):
            if a["worst_rel_err"] > 0 and b["worst_rel_err"] > 0:
                p = np.log2(a["worst_rel_err"] / b["worst_rel_err"])
                print(f"    N {a['N']:>4} -> {b['N']:<4}   p ~ {p:5.2f}")

    if args.json:
        with open(args.json, "w") as fh:
            json.dump({"scheme": "hopf_fourier_staggered_fd",
                       "lam_max": args.lam_max, "m_max": args.m_max,
                       "cluster_tol": args.cluster_tol, "ladder": rows},
                      fh, indent=2, sort_keys=True)
            fh.write("\n")
        print(f"\nwrote {args.json}")


if __name__ == "__main__":
    main()
