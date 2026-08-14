#!/usr/bin/env python3
"""
M8.5-B engineering pilot, route (a) prototype 2: the identification map on a
cyclic quotient of S^3.

NON-EVIDENTIARY pilot code (frozen pre-registration section 6). Tuning set only.
The sealed rung-3a/3b adjudication case is NOT touched here and is excluded from
the tuning set by name.

THE IDENTIFICATION, AND WHY THIS COORDINATE CHOICE
  In Hopf coordinates the generator of a cyclic deck group acts as a pure shift of
  the two periodic fibres, leaving eta alone:

      (eta, xi1, xi2)  ->  (eta, xi1 + 2 pi / p, xi2 + 2 pi q / p)

  A Fourier mode e^{i m1 xi1} e^{i m2 xi2} picks up the phase
  exp(2 pi i (m1 + q m2) / p), so the equivariance constraint is exactly

      m1 + q * m2  ==  0   (mod p)

  The radial problem in eta is untouched by the quotient. So the identification map
  costs nothing here: it is a sublattice selection on the Fourier sweep, not a
  ghost-cell bookkeeping exercise. That is the measurement this prototype exists to
  make, since "the identification map is fiddly" is the named route (a) risk.

CHARACTER-FREE, as the frozen section 2 rule requires
  Multiplicities are obtained by applying the numerical group action to the computed
  modes and counting those it fixes. No character table, no irreducible labels, and
  no representation-theoretic input of any kind appears in this file. Route (b) owns
  the character-averaging prediction.
"""

import argparse
import json
import time

import numpy as np

from route_a_scalar_s3 import radial_spectrum, cluster

# Tuning set: small cyclic quotients, preregistered for the pilot.
# The sealed adjudication case is deliberately NOT in this list and is never run here.
TUNING_SET = {
    "L(2,1)": (2, 1),
    "L(3,1)": (3, 1),
    "L(4,1)": (4, 1),
}

# Published scalar multiplicities for the tuning cases, k = 0..9, unit radius.
# Reference values for the tuning set only. Used to confirm the identification map
# is right before any parameter is frozen.
TUNING_REFERENCE = {
    "L(2,1)": [1, 0, 9, 0, 25, 0, 49, 0, 81, 0],
    "L(3,1)": [1, 0, 3, 8, 5, 12, 21, 16, 27, 40],
    "L(4,1)": [1, 0, 3, 0, 15, 0, 21, 0, 45, 0],
}


def invariant_modes(p, q, m_max, lam_max):
    """The (m1, m2) sectors surviving the equivariance constraint."""
    keep = []
    for m1 in range(-m_max, m_max + 1):
        for m2 in range(-m_max, m_max + 1):
            if (m1 + q * m2) % p != 0:
                continue
            if m1 * m1 + m2 * m2 > lam_max + 1:
                continue
            keep.append((m1, m2))
    return keep


def quotient_spectrum(p, q, N, m_max, lam_max):
    t0 = time.time()
    sectors = invariant_modes(p, q, m_max, lam_max)
    found = []
    for (m1, m2) in sectors:
        vals = radial_spectrum(m1, m2, N, n_want=max(4, int(np.sqrt(lam_max)) + 2))
        found.extend(v for v in vals if v <= lam_max + 1e-6)
    return np.sort(np.array(found)), time.time() - t0, len(sectors)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=320)
    ap.add_argument("--k-max", type=int, default=9)
    ap.add_argument("--m-max", type=int, default=12)
    ap.add_argument("--cluster-tol", type=float, default=0.35)
    ap.add_argument("--json")
    args = ap.parse_args()

    lam_max = args.k_max * (args.k_max + 2)
    print("route (a) prototype 2: cyclic quotient via Fourier selection rule")
    print(f"N = {args.N}, levels k = 0..{args.k_max} (lambda <= {lam_max})")
    print("tuning set:", ", ".join(TUNING_SET), " (sealed case excluded by construction)")
    print()

    results = {}
    for name, (p, q) in TUNING_SET.items():
        vals, elapsed, nsec = quotient_spectrum(p, q, args.N, args.m_max, lam_max)
        clusters = cluster(vals, args.cluster_tol)

        # bin clusters onto the k(k+2) ladder
        got = {}
        for c in clusters:
            mean = float(np.mean(c))
            k = int(round((np.sqrt(mean + 1) - 1)))
            if abs(k * (k + 2) - mean) < 0.35:
                got[k] = got.get(k, 0) + len(c)

        ref = TUNING_REFERENCE[name]
        line = [got.get(k, 0) for k in range(args.k_max + 1)]
        match = line == ref[:args.k_max + 1]
        print(f"  {name}  p={p} q={q}   {elapsed:5.2f}s   sectors={nsec}")
        print(f"    computed : {line}")
        print(f"    reference: {ref[:args.k_max + 1]}")
        print(f"    {'MATCH' if match else 'MISMATCH'}")
        print()
        results[name] = {"p": p, "q": q, "seconds": round(elapsed, 3),
                         "sectors": nsec, "computed": line,
                         "reference": ref[:args.k_max + 1], "match": bool(match)}

    allok = all(r["match"] for r in results.values())
    print("VERDICT:", "identification map reproduces every tuning case"
          if allok else "IDENTIFICATION MAP FAILED on at least one case")

    if args.json:
        with open(args.json, "w") as fh:
            json.dump({"N": args.N, "k_max": args.k_max,
                       "scheme": "hopf_fourier_selection_rule",
                       "note": "pilot, non-evidentiary; sealed case not present",
                       "cases": results}, fh, indent=2, sort_keys=True)
            fh.write("\n")
        print(f"wrote {args.json}")


if __name__ == "__main__":
    main()
