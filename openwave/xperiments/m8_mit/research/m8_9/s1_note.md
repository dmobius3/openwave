# M8.9 S1: INSTRUMENT DEFECT. The `R_0` block does not carry its analytic content

> **Run 2026-08-25 under the filed decision rule**
> ([`../findings/m8_9_s1_decision_rule.md`](../findings/m8_9_s1_decision_rule.md)), frozen region
> SHA-256 `68df11a0…4938934`, verified before any spectrum was computed. G-MULT failed. **No
> S1-A/B/C classification is issued and `I_star` was not formed**, because the filed rule computes
> it only on a passing gate. S1 did not answer M8.9's question. It established that this control
> cannot answer it.

## Provenance

| Gate | Result |
| --- | --- |
| `main` contains merge `9ebb0ffc` | ✅ ancestor of the run commit |
| frozen-rule SHA-256 recomputed from the filed artifact | ✅ `68df11a02ee5097d23712b4eace9b65220d2f540ea147660ce3cd8ae4a938934` |
| import closure over `p0/`, `p1a/`, `m8_5b/pilot`, `m8_5b/production` | ✅ clean, no workaround |
| operator route identical to P1A: `build_orbit_cloud(fibonacci_seeds_s3(60))`, `build_L_bundle(..., reps["R0"], k=110, m=7, p=4)` | ✅ 7200-node cloud, 60x60 block |
| spectrum route identical to P1A: `np.linalg.eigvals(L)`, `M_h` used for the inner product only | ✅ `extract_and_score` takes the same route |
| independent tie to the merged P1A record | ✅ the level-0 mode reads `5.954e-14` against the shipped `R_0` value `5.95e-14` |

## Two runs, and the first one is on the record

**Run 1: INVALID ENGINEERING REHEARSAL.** Mutation discrimination was not established before the
spectral read. The implementation read the live spectrum first, found G-MULT already red, and then
applied the filed mutation to that already-failing spectrum. Observing failure there proves
nothing: a gate that already fails cannot demonstrate discrimination by failing again. The arm was
vacuous and the run does not qualify.

**Run 2: QUALIFICATION RUN.** Governing rule unchanged; the filed mutation implemented correctly,
on a known-green parent, before the live spectral read. This note reports run 2.

The repair was mechanically dictated by the already-filed rule, which specifies the mutation and
places the arm among the pre-execution steps. No adjudication criterion was touched: windows,
multiplicities, thresholds and branches were all publicly frozen and independently verified before
either run. But the spectrum HAD been seen when the repair was made, and hiding that ordering
failure would make this look cleaner than it was.

## G-MULT

**Mutation arm, run BEFORE the real spectrum.** Arming against the real spectrum would have been
vacuous, because a gate that already fails cannot demonstrate anything by failing again. The arm
therefore ran on the exact analytic spectrum, which passes by construction: moving ONE eigenvalue
from 168 to 305, across the frozen 304 midpoint, with every other eigenvalue unchanged, moved
`C_12` by `-1` and `C_20` by `+1` and failed the gate. **The gate discriminates.**

**On the real spectrum: FAIL.**

| cluster | window | required | found | status |
| --- | --- | --- | --- | --- |
| `C_0` | `Re < 84` | 1 | 1 | ok |
| `C_12` | `84 <= Re < 304` | 13 | 13 | ok |
| `C_20` | `304 <= Re < 532` | 21 | **20** | mismatch |
| `C_24` | `Re >= 532` | 25 | **26** | mismatch |

## The record

Reported because the commission requires it, not because any of it licenses a branch.

| cluster | center | Re spread | `max\|Im lambda\|` | Schur residual |
| --- | --- | --- | --- | --- |
| `C_0` | -0.0000 | 0.000e+00 | 5.954e-14 | 2.641e-12 |
| `C_12` | 168.7638 | 1.109e+01 | 3.035e+00 † | 7.001e-12 |
| `C_20` | 445.8225 | 1.826e+02 | 2.241e+01 | 7.777e-12 |
| `C_24` | 813.3837 | 7.693e+02 | 3.126e+01 | 8.474e-12 |

† **Non-adjudicative observation from an instrument-failed run.** `C_12`'s 3.035 is recorded
prominently because it is a reason a successor experiment is worthwhile, and for no other reason.
It is not evidence for any branch.

## What the defect is

**Language note.** G-MULT's failure is exactly the statement that the analytic level identities
are NOT established for this block. Everything below is therefore phrased in terms of the frozen
assignment WINDOWS, never in terms of which continuum level an eigenvalue belongs to. Saying "the
twenty-first member of the level-20 band" would presuppose the thing the gate could not certify.

The count deficit in `C_20` and the excess in `C_24` are produced at the 532 boundary: the nearest
eigenvalue on the `C_24` side is `Re lambda = 549.739`. That is the proximate failure, and it is a
symptom of a larger one.

**The block does not reproduce the analytic decomposition the rule predicted.** Eigenvalues
assigned to the `C_12` window span 164.1 to 175.2 about that window's analytic center 168. Those
assigned to `C_20` span 347.7 to 530.3 about 440. Those assigned to `C_24` span 549.7 to **1319.0**,
while the block's analytic content stops at 624. Four tight clusters were predicted, and the upper
windows are not clustered at all.

**The contamination is not confined to the top of the block.** 28 of the 60 eigenvalues carry
`|Im| > 1e-06`, in 14 complex-conjugate pairs, and such pairs already appear WITHIN the `C_12`
assignment window, the lowest of them at `Re lambda = 166.32`, rather than only near the
finite-dimensional cutoff. The lowest mode is clean at 5.954e-14. Per window: `C_0` 0 of 1 dirty,
`C_12` 8 of 13, `C_20` 8 of 20, `C_24` 12 of 26.

## Why `C_12`'s exact count does not rescue a classification

`C_12` returned exactly its required 13 and its `max|Im|` of 3.035 is some 228 times `I_high_min`.
Reading that as S1-A is the obvious temptation and it is barred, for a reason this column already
had in writing before S1 was designed. `route_a_repn.py`: **count-only safety is not subspace
correctness.** A window holding the right NUMBER of eigenvalues is not evidence that its contents
are the level-12 eigenspace, and with the `C_20` and `C_24` windows
demonstrably failing their counts, there is no ground for assuming the `C_12` window's contents are
faithful either. The filed rule requires all four
counts, and it requires them precisely so that a partially-correct spectrum cannot be mined for
the branch its cleanest window happens to support.

Ordering G-MULT ahead of classification is what made this decidable rather than arguable.

## Disposition

**INSTRUMENT DEFECT.** No branch. No licensed sentence. `I_star` not formed. **No high-level
localization conclusion is licensed, including from the apparently dirty `C_12` window.**

The filed rule routes to S2 from S1-B and S1-C, and this is neither, so **S1 does not trigger S2**.
M8.9's question stands exactly where it stood before the run: unanswered.

What S1 does establish, and it is a real result about the instrument rather than about the physics,
is that the `R_0` production block at 60 seeds cannot serve as a high-level trivial-fibre control.
The design assumed a 60-dimensional block spanned by four analytic levels would resolve all four.
It resolves the first two and not the rest. Any successor control must either raise the resolution
until the analytic decomposition is reproduced, or find a trivial-fibre probe that does not depend
on the top of a block being faithful.

**Designing that successor is not S1's job and is not authorized here.** S1 stops.
