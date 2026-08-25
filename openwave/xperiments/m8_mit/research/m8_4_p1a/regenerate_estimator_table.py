#!/usr/bin/env python3
"""Regenerate the estimator-qualification block quoted in the P1A closeout.

The P1A.4a and P1A.4b notes record P1A.1 to P1A.3 as "PASS, unchanged, not rerun", so the
numbers in the closeout's "What qualified" table have no standing record in those notes.
This script recomputes them from the shipped p0/ and p1a/ packages, at the production cloud,
so the claim is reproducible from the evidence package rather than resting on an overwritten
earlier note.

    PYTHONPATH=.:../m8_5b/pilot:../m8_5b/production python3 regenerate_estimator_table.py

Runtime is a few minutes: nine bundles at 60 seeds, each with a Schur extraction and a
128-point Riesz contour.
"""
import sys
import numpy as np

from p0 import group as G, representations as RP
from p0.cloud import fibonacci_seeds_s3, build_orbit_cloud
from p0.bundle_operator import build_L_bundle
from p1a.mass_matrix import build_Mh_base, build_Mh_rho
from p1a.diagnostics import henrici_departure, riesz_projector, riesz_projector_norm
from p1a.subspace import (extract_and_score, projector_invariants,
                          three_way_identity_check, riesz_crosscheck)

N_SEEDS, K, RBF_M, RBF_P = 60, 110, 7, 4
LABELS = ["R0", "R1", "R3", "R6", "R7", "R8", "R4", "R5", "R2"]
DIST = {"R0": 0, "R1": 1, "R3": 2, "R6": 3, "R7": 4, "R8": 5, "R4": 6, "R5": 6, "R2": 7}


def main():
    elems = G.build_icosians()
    chi = G.build_character_table(elems)
    out = RP.build_all_representations(elems, chi)
    reps = out[0] if isinstance(out, tuple) else out
    X, oid, gid = build_orbit_cloud(fibonacci_seeds_s3(N_SEEDS), elems)
    W = build_Mh_base(X, oid, N_SEEDS)

    print(f"cloud {len(X)} nodes, k={K}\n")
    print(f"{'sector':7s} {'idem':>10} {'Mh-sym':>10} {'inv resid':>11} "
          f"{'inv rel':>11} {'Riesz th':>10} {'3-way':>10} {'Henrici':>8} {'|P_spec|':>9}")
    agg = {k: [] for k in ("idem", "sym", "abs", "inv", "riesz", "three", "hen", "pspec")}
    for lab in LABELS:
        d = DIST[lab]; k = d + 1; lam0 = d * (d + 2)
        L = np.asarray(build_L_bundle(X, oid, gid, elems, reps[lab], k=K, m=RBF_M, p=RBF_P)[0])
        Mh = build_Mh_rho(W, reps[lab][0].shape[0])
        ex, Q = extract_and_score(L, Mh, lam0, k, label=lab)
        pv = projector_invariants(Q, L, Mh, k)
        tw = three_way_identity_check(Q, Q, Mh)
        radius = max(1.0, 0.25 * max(lam0, 1))
        rz = riesz_crosscheck(L, Q, Mh, lam0, radius, k)
        pspec = riesz_projector_norm(riesz_projector(L, lam0, radius))
        hen = henrici_departure(L)
        row = (pv["idempotence"], pv["mh_symmetry"], pv["invariance_residual"],
               pv["invariance_residual_rel"],
               rz["theta_max_schur_vs_riesz"], tw["max_disagreement"], hen, pspec)
        for key, v in zip(agg, row):
            agg[key].append(float(v))
        print(f"{lab:7s} " + " ".join(f"{v:>10.2e}" for v in row[:6])
              + f" {row[6]:>8.3f} {row[7]:>9.2f}")
    print("\nranges across the nine bundles:")
    for key, name in (("idem", "projector idempotence"), ("sym", "M_h-symmetry"),
                      ("abs", "invariance residual (absolute)"),
                      ("inv", "invariance residual (relative)"), ("riesz", "Schur vs Riesz theta_max"),
                      ("three", "three-way leakage identity"), ("hen", "Henrici departure"),
                      ("pspec", "||P_spec|| on target clusters")):
        v = agg[key]
        print(f"  {name:32s} {min(v):.3e} to {max(v):.3e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
