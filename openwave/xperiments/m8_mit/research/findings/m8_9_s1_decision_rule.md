# M8.9 S1 decision rule, frozen before the run

> **Status: FROZEN and FILED. S1 is not commissioned and this document does not authorize the
> run.** It is filed now, rather than kept author-side, so that when S1 reports, its governing rule
> is already public and timestamped rather than cited from a document outside the evidence package.
> Its only purpose is to fix the reading of S1's outcome BEFORE the 60x60 matrix is diagonalized,
> so a suggestive number cannot pick its own interpretation afterwards. Nothing below is contingent
> on any S1 output, because no S1 high-level `R_0` spectrum has been computed.
>
> Claim boundary, inherited from the #468 merge commitment: locate the source of the unphysical
> complex spectrum. No repair, no nonlinear dynamics, no reopening of P1A. P1A's outcome is final
> and this task cannot change it.

## The question, and why the existing table cannot answer it

M8.4 P1A established that the imaginary parts in the computed bundle spectra belong to the
assembled matrix rather than to the eigensolver: with `L_h` held byte-for-byte fixed, they are
unchanged from float64 through 50-digit arithmetic. What P1A did NOT establish is which upstream
stage produces them, the equivariant quotient assembly or the RBF-FD discretization underneath it.

The shipped contamination table cannot separate them, because the two candidate causes are
perfectly confounded. Each nontrivial bundle was inspected only at its own first-section level, so

    rho  <->  d_rho  <->  lambda = d_rho (d_rho + 2)

move together across every row. The rise from 1.4e-10 at lambda = 3 to 1.8e-01 at lambda = 63 is
equally consistent with "high harmonic level is hard for this discretization" and with "nontrivial
fibre transport manufactures the complex part".

## What S1 measures

The trivial bundle `R_0` supplies the missing cell: a TRIVIAL fibre at HIGH harmonic level. Its
production block is exactly 60x60 at 60 seeds, and the 2I-invariant harmonics fill it exactly:

| n | lambda = n(n+2) | multiplicity | cumulative |
| --- | --- | --- | --- |
| 0 | 0 | 1 | 1 |
| 12 | 168 | 13 | 14 |
| 20 | 440 | 21 | 35 |
| 24 | 624 | 25 | **60** |

The next invariant level, n = 30, would need 31 more dimensions and does not fit. So the block's
entire analytic content is four clusters, and their multiplicities are an exact prediction rather
than a fitted expectation. P1A examined only the lambda = 0 cluster.

Per cluster, S1 records: center, real spread, `max|Im lambda|`, Schur extraction residual, and
multiplicity against the table above. Same 60-seed cloud, same `k = 110`, same RBF parameters,
same shipped `R_0` matrix. Nothing is rebuilt differently.

**G-MULT, and the assignment is frozen so the gate can actually fail.** "Resolves into four
clusters" is not self-defining. An implementation that sorts the 60 eigenvalues and slices them
1/13/21/25 would pass by construction, which is a tautological gate, not a test. Assignment is
therefore by nearest analytic level in `Re lambda`, with the Voronoi boundaries fixed at the
midpoints between exact centers, 84, 304 and 532:

| cluster | assignment window | required count |
| --- | --- | --- |
| `C_0` | `Re lambda < 84` | 1 |
| `C_12` | `84 <= Re lambda < 304` | 13 |
| `C_20` | `304 <= Re lambda < 532` | 21 |
| `C_24` | `Re lambda >= 532` | 25 |

G-MULT passes only when all four counts match. A single eigenvalue drifting across a frozen
midpoint changes a count and fails it, which is what makes it falsifiable. The gate is
mutation-armed under the standing rule, and the mutation is specified rather than left to choice:
MOVE EXACTLY ONE eigenvalue across ONE frozen midpoint, leaving every other eigenvalue unchanged;
confirm that one adjacent count decreases by one, the other increases by one, and G-MULT fails.
Record both. A swap of two eigenvalues across the same boundary is explicitly NOT an acceptable
mutation: it preserves both counts, so G-MULT stays green legitimately and the arm would prove
nothing. Failing G-MULT means the block is not carrying the analytic content
this test assumes, and S1 reports an instrument defect rather than a localization result.

## The decision rule, and it is deliberately ASYMMETRIC

The two outcomes do NOT carry equal weight, and the classification is three-way rather than
binary so that an unforeseen intermediate result cannot be forced into either branch.

**The thresholds are the two regimes already observed in closed P1A data, not new inventions.**
At the production cloud the low-contamination group runs up to `I_low_max = 1.27e-08` (`R_7`, the
dirtiest of `R_0`, `R_1`, `R_3`, `R_6`, `R_7`), and the high-contamination group starts at
`I_high_min = 1.33e-02` (`R_8`, the cleanest of `R_8`, `R_5`, `R_4`, `R_2`). Those two regimes are
separated by a factor of 1.05e+06 and are exactly what motivated S1 in the first place. Nothing
here is selected from S1's own output.

Write `I` for `max|Im lambda|` on a cluster. Only `C_12` (lambda = 168) and `C_20`
(lambda = 440) adjudicate, and they adjudicate through a SINGLE statistic, frozen here:

    I_star = max(I_12, I_20)

Routing on one number rather than on two separately-quantified clauses is what makes the three
branches provably disjoint and exhaustive: `I_star` is a single real, the two thresholds satisfy
`I_low_max < I_high_min`, and the three conditions below partition the line. An earlier draft
quantified each branch over "either cluster" independently, which was NOT a partition: a run
returning `I_12 = 2e-02` and `I_20 = 1e-05` fired both S1-A and S1-C at once.

**S1-A, base discretization strongly implicated.** `I_star >= I_high_min`. Equivalently, at least
one diagnostic trivial-fibre cluster reaches the previously observed dirty regime. Licensed
conclusion: nontrivial fibre transport is NOT NECESSARY for high-scale contamination, and the base
discretization is strongly implicated. This branch would substantially end the search.

**S1-B, the simple high-lambda explanation is rejected.** `I_star <= I_low_max`. Equivalently,
BOTH diagnostic clusters stay inside the prior clean regime. Licensed conclusion, and no more than
this: high harmonic level by itself, in the trivial quotient sector, is INSUFFICIENT to reproduce
the high-contamination regime. This does NOT exonerate RBF-FD and does NOT convict the fibre
transport. The defect may be an INTERACTION between the non-self-adjoint stencil, the equivariant
reduction, and a nontrivial fibre, with no single ingredient defective alone; such an interaction
vanishes in `R_0` by construction. S2 is required.

**S1-C, indeterminate.** `I_low_max < I_star < I_high_min`. Equivalently, neither diagnostic
cluster reaches the dirty regime, but at least one has left the clean one. S1 yields NO
localization conclusion and S2 is required. This branch exists so that a result at, say, 1e-05 is
reported as indeterminate rather than argued into a regime it does not belong to.

**lambda = 624 carries no S1 decision weight.** `C_24` is the edge-of-block cluster, consuming the
final 25 dimensions with no represented invariant level above it, so its numerical error is
confounded with resolution and truncation exhaustion. Its values are RECORDED and its multiplicity
remains part of G-MULT, but its contamination does not choose an S1 branch. A dirty `C_24` is not
evidence of cleanliness either: it simply does not affect S1 classification.

Recording that before the run is what stops a dirty edge cluster from being read as
"discretization implicated" when it is only the block running out of room.

## S2, named here so S1 is not asked to carry it

S2 compares two mechanisms for imposing equivariance on the SAME cloud and the SAME RBF weights:
the production quotient transport assembly against an independent subspace restriction of the
cover-space operator by projector or intertwiner. Three sectors suffice: `R_0` as control, one
clean low-level nontrivial sector, and `R_2` as the dirty extreme. Both dirty implicates the
collocation operator; independent-restriction clean with production dirty localizes an assembly
defect; both dirty but quantitatively different measures each layer's contribution to an
interaction, which is still an answer.

**Pre-commitment, made now to remove a later degree of freedom.** The "clean low-level nontrivial
sector" is named `R_1`: it has the lowest nontrivial first-section level, `d_rho = 1` at
`lambda = 3`. Naming it before any S2 output exists costs nothing now and cannot be chosen to suit
a result later.

S1 is a control, not a verdict. Its job is to decide whether S2 is needed.
<!-- FREEZE-BOUNDARY -->

**Freeze record.** SHA-256 covers every byte ABOVE the `FREEZE-BOUNDARY` marker: `68df11a02ee5097d23712b4eace9b65220d2f540ea147660ce3cd8ae4a938934`

```bash
sed '/FREEZE-BOUNDARY/,$d' m8_9_s1_decision_rule.md | shasum -a 256
```
