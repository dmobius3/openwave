# N4b design , close the N4 pending tasks for a complete Duda round (no loose ends)

The maintainer (2026-06-22): nail items 2, 3, 6 + refine 4 as ONE bundle "N4b" (extension of N4), so Duda gets a
complete round. DEFER N5 (article) + N6 (masses) behind it.

| N4b sub-task | Source item | Question | Script |
| --- | --- | --- | --- |
| **N4b-1** | item 2 | does the TBM gate + delta_CP survive a real LdG TENSOR potential scan `(a,b,c)`? (Duda's #1 concern: "the potential is crucial") | `m5_11_n4b_potential.py` |
| **N4b-2** | item 3 | the theta12 (35.26 vs 33.68) + theta23 residuals: can a minimal natural breaking absorb them, and at what cost to the predictions? | `m5_11_n4b_residual.py` |
| **N4b-3** | item 6 | where does `g_chiral` come from + its natural scale? (the achiral substrate gives NO CP; the chiral Lifshitz term is required; g_chiral~O(1) -> theta13~O(10deg)) | `m5_11_n4b_chiral_origin.py` |
| **N4b-4** | item 4 (refine) | crisp write-up of the TWO delta scales (mixing-delta ~0.1 vs quantum-phase delta ~1e-10) | findings doc + small numeric illustration |

## Method notes

- **N4b-1.** Replace the crude on-site `kappa*P` with the REAL LdG potential Hessian projected onto flavour
  displacements: `H^V_ab = INT [2a Tr(dA dB) - 6b Tr(Mvac dA dB) + c(8 Tr(Mvac dA)Tr(Mvac dB) + 4 Tr(Mvac^2)
  Tr(dA dB))]` (spatial 3x3 block). M_mass = K(kinetic) + H^V. Re-find magic alpha* per (a,b,c); check the 3
  TBM angles + (with the chiral term) delta_CP. Expect: robust by symmetry (Mvac is mu-tau-symmetric -> H^V
  preserves mu-tau + magic), confirmed numerically.
- **N4b-2.** Quantify residuals; add ONE small real mu-tau / magic breaking; map theta12 toward 33.68; report
  whether it spoils theta23/theta13/delta_CP (the predictions-vs-fit trade-off). Honest, not forced.
- **N4b-3.** Show C (chiral overlap) = 0 for achiral loops (chi=0) -> achiral substrate = CP-conserving;
  chiral screw sources C -> CP. g_chiral* ~ 0.94 (O(1)) -> theta13 ~ O(10 deg) NATURAL. The open substrate
  question for Duda: does the M5 LdG carry a chiral (Lifshitz) term to favour the screw?
- **N4b-4.** Document: the mixing sees an O(0.1) effective biaxiality, NOT the 1e-10 quantum-phase delta;
  the two are distinct scales (N3 finding), reconciled by the chiral mechanism.

## Discipline
Headless numpy f64, 16-core. Checkpoint each. #236 HELD, git untouched. Total invisibility. Cross-link
n4b_findings <-> n4_findings <-> 10a <-> #199. Output: `n4b_findings.md` + scripts + a robustness figure.
