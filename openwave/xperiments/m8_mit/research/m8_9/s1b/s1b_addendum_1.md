# S1b addendum 1: three defects the independent qualification exposed

> **APPEND-ONLY with respect to the parent.** The frozen region of `S1B_DECISION_RULE.md`, SHA-256
> `c44c603a8058ed8529e5bb0f42ec168b443ff22b050d3748f1f1de3537c7d297`, is UNCHANGED and still
> verifies. This addendum adds to it and has its own freeze record.
>
> **A first draft of this addendum, SHA-256 `0847f3b5652e11cc5fd481be464449655557da435a026582bcdd1f4f6a402acb`,
> is WITHDRAWN and was never adopted.** It mislabelled the `k = 2` control and under-scoped the SVD
> rerun. It is retained as `S1B_ADDENDUM_1_WITHDRAWN_0847f3b5.md` so the record shows the
> correction rather than a clean text.
>
> **PRECEDENCE.** For the ladder controls, the SVD evidence transfer, the principal-angle
> measurement method, and the resulting `K_floor` quantities covered below, this addendum SUPERSEDES
> the corresponding parent statements while leaving the parent record intact. Where both carry a
> number for the same quantity, the addendum's governs and the parent's is historical.
>
> Everything here was found BEFORE any live `R_0` spectrum at `n = 12` or `n = 20` was computed.
> Nothing is derived from a target observation.

## A1.1 The collapse control's missing parameter, and what `k = 2` actually exercises

**The defect.** The parent arms the ladder with "the companion matrix of `(x-1)^k` ... must
COLLAPSE" and never fixes `k`. That omission is OUTCOME-BEARING:

| `k` | `J(64)` | `J(30)` | `J(50)` | verdict | by |
| --- | --- | --- | --- | --- | --- |
| 2 | 0 | 0 | 0 | COLLAPSES | rule 1, exact zero |
| 3 | 6.469e-06 | 9.101e-11 | 2.193e-17 | AMBIGUOUS | rule 2, tail |
| 4 | 2.390e-04 | 3.052e-08 | 4.694e-13 | AMBIGUOUS | rule 2, tail |
| 10 | 5.048e-02 | 1.610e-03 | 1.652e-05 | AMBIGUOUS | rule 2, tail |

**`k = 2` is the RULE-1 EXACT-ZERO control and nothing more.** A first draft of this addendum
described it as showing "an arithmetic artifact resolving." It does not: its three readings are
identically zero, so it never reaches the ratio tests at all. What it demonstrates is that a
minimally defective real-spectrum matrix lands on the exact-zero branch in this environment.
Calling it an artifact-collapse demonstration would be green-for-the-wrong-reason, which is the
failure family this program has now paid for seven times.

**Frozen ladder controls, each naming the rule it exercises:**

| control | input | required |
| --- | --- | --- |
| exact-zero control | companion of `(x-1)^2` | `COLLAPSES`, reached via rule 1 |
| ill-conditioning refusal control | companion of `(x-1)^3` | `AMBIGUOUS`, reached via rule 2 |
| rule-4 reachability control | the synthetic triplet `J = (1e-8, 1e-12, 1e-12)` | `COLLAPSES`, reached via rule 4 |

The third is adjudicator-only: it feeds the ladder a manufactured triple with a stable tail and
`J(50) <= 1e-3 J(64)`, and exists solely to show rule 4 is reachable without pretending `k = 2`
exercises it. Verified, together with its edges: `(1e-8, 1e-12, 2e-12)` and `(1e-8, 1e-12, 5e-13)`
both COLLAPSE
by rule 4, at tail ratios 2.0 and 0.5 exactly, which rule 2 admits because its tests are strict.

**A fourth control, and a correction.** A previous draft offered `(1e-8, 1e-11, 1e-11)` as a
fall-through specimen returning AMBIGUOUS by rule 5. It does not. Rule 4 tests
`J(50) <= 1e-3 * J(64)`, and `1e-3 * 1e-8` evaluates to `1.0000000000000001e-11` in float64, so the
comparison holds and the triplet COLLAPSES by rule 4. In exact arithmetic it holds at equality.
Either way the specimen's verdict turns on rounding at a boundary, which makes it unusable as a
control regardless of which side it lands on. The verified specimen behind that draft was in fact
`(1e-8, 1e-11, 1.5e-11)`, and it was transcribed into the addendum with the wrong third entry.

| rule-5 fall-through control | `J = (1e-8, 2e-11, 2e-11)` | `AMBIGUOUS`, reached via rule 5 |

Verified: stable tail, rule 3 false, rule 4 false at a margin of 2.0x rather than at equality.

**Why higher `k` behaves that way, and why rule 2 is NOT changed.** The eigenvalue 1 of the
companion of `(x-1)^k` has algebraic multiplicity `k` and geometric multiplicity 1, so a
perturbation `delta` splits it into roots at distance of order `delta^(1/k)`. Raising precision
shrinks `delta`, so `J` keeps falling and never settles between rungs. Refusing to classify a
maximally defective spectrum is the correct conservative response.

## A1.2 A standing open question, recorded and deliberately NOT resolved

`k >= 3` collapsing monotonically into AMBIGUOUS raises a question broader than this control: a `J`
that is a PURE arithmetic artifact shrinks at every rung, so it can trip rule 2 before reaching
rule 4. Rule 4 then licenses `COLLAPSES` only for a `J` with a genuine nonzero floor beneath the
float64 noise. Since `S1b-NULL` and `S1b-ADJOINT` both require every `J` collapsed or zero, that
would narrow the reachable outcomes.

**Not fixed here.** Rule 1 catches exact zero, and for a REAL `A_n` LAPACK returns exactly-zero
imaginary parts for real eigenvalues rather than a smooth artifact, so the degenerate regime may
not arise for the actual target at all. Whether it does is unknown and must stay unknown until the
instrument is authorized. Changing rule 2 on speculation would be the redesign this addendum
exists to avoid.

## A1.3 The SVD: the evidence chain is repaired by a bridge, not by assertion

**What happened.** The parent names `invariant_dim_and_basis(pairs, n)[1]` for `C_n^svd`. The
qualification substituted a local reimplementation differing only in `full_matrices=False` and used
it to build the qualifying `C_svd`. Verified by inspection: the shipped function is called ZERO
times in `s1b_qualification.py`; its sole occurrence is in the replacement's docstring. Every
qualification result downstream of `C_n^svd` was therefore obtained from a non-contract
construction, and dimension agreement alone does not transfer that evidence.

**Runnability, settled as a binary rather than left ambiguous.** With `full_matrices=True` the SVD
materializes an `m x m` `U` that is discarded: 20280 x 20280 at `n = 12`, 6.58 GB, and
52920 x 52920 at `n = 20`, 44.8 GB nominal against 25.8 GB of physical RAM on this machine. It
nonetheless RAN TO COMPLETION at `n = 20`, once, in 3918.5 s, presumably through compressed memory.
**The shipped call is therefore expensive and near the machine's limit, but demonstrated runnable,
and this addendum requires it.** It is NOT declared impracticable. Reproducibility at that size is
not guaranteed, which is exactly why the bridge below requires it ONCE per target rather than for
every downstream surface.

**G-SVD-BRIDGE, and its transfer rule is bitwise identity, not a perturbation bound.** Run the
shipped `invariant_dim_and_basis` once at each of `n = 12, 20`. Let `C^ship` be its basis and
`C^econ` the economy construction already used.

**Why not a perturbation bridge.** A previous draft licensed transfer on
`2 sin(theta_bridge/2) ‖B‖_2 <= 100 eps ‖A_n‖_2` with `B = Ltilde - Ltilde^H`. That bound is valid
for `K(Q) = (1/2) Q^H B Q` and for nothing else. It does NOT bound the change in the full
compression `A(Q) = Q^H Ltilde Q`, whose relevant scale is `‖Ltilde‖_2` rather than `‖B‖_2`, and it
therefore cannot transfer `J`-bearing evidence at all: `A` may be nonnormal, so an arbitrarily small
perturbation can move its eigenvalues disproportionately, which is the reason the precision ladder
exists. Building a second uncertainty estimator to rescue the cheap part of a rerun, when the
expensive part must be paid regardless, would only create another estimator needing its own
qualification. Withdrawn.

**The transfer rule, frozen:**

    bitwise-identical `s` AND `Vh` between the shipped and economy calls  ->  evidence TRANSFERS
    anything else                                                        ->  RERUN the downstream
                                                                             SVD-route surfaces

Bitwise identity means `C^ship == C^econ` byte for byte, so the same deterministic sampling pipeline
receives the same coefficients and there is no evidence-transfer question to answer. Note what that
does and does not imply numerically: it makes the underlying SUBSPACES identical, and the measured
sine residual between them should then fall at the numerical projector floor, about 7e-16, NOT at
bitwise zero. Claiming otherwise would repeat this addendum's own central lesson, that a numerical
angle formula can manufacture apparent structure where none exists. Measured at
`n = 12`: `‖s_ship - s_econ‖_inf = 0.000e+00` and `Vh` bitwise identical, at 804.9 s against 1.2 s.
Whether it holds at `n = 20` is unknown and is not assumed.

On the fallback branch, rerun every surface downstream of `C_n^svd` on `C^ship`: G-RANK, the SVD
half of G-SUBSPACE, the SVD half of G-SAMPLE, and any check that consumed `Q_n^svd`. The shipped
`n = 20` call is required on both branches, so the fallback costs only the cheap downstream work.

**Recorded as DIAGNOSTICS, routing nothing.** Both sine-form comparisons are measured and reported
because they are informative about the two constructions, and neither licenses transfer on its own:

- coefficient space, `‖(I - P_ship) Q_econ‖_2`, reported against an engineering reference of
  `1e-12`. **That number is a declared reporting threshold, not a frozen gate**, and it is not
  asked to certify anything;
- weighted sampled space, the same quantity between `ran Qtilde^ship` and `ran Qtilde^econ`,
  reported against the same reference. A coefficient-space agreement does not imply a sampled-space
  one, since `C -> F^seed -> W -> Qtilde` can amplify by up to `kappa(W)`; that is the `theta_C`
  versus `theta_Q` lesson of A1.4 one level down, and it is why both are reported.

**Arms, because a gate and its diagnostics must both be shown to have power.** None needs the
65-minute SVD:

| arm | construction | required |
| --- | --- | --- |
| bitwise, green | `C^econ` against itself | identical, transfer path taken |
| bitwise, red | one entry of `Vh` perturbed by 1 ulp | NOT identical, fallback path taken |
| diagnostic, green | untilted basis | residual at the projector floor, about 7e-16 |
| diagnostic, red | one column tilted by a frozen orthogonal component, `sin(phi) = 1e-10` | residual 1e-10, two orders above the `1e-12` reference |

The sampled-space diagnostic takes the same tilt, propagated through the fixed sampling map. Its
resulting residual is RECORDED, not required to remain at 1e-10: the map may amplify or contract
the perturbation by up to `kappa(W)`. Only the coefficient-space arm carries the analytic 1e-10
target; the sampled-space version is a propagated diagnostic and is read as such.

**The economy SVD is NOT adopted into the live path.** The bridge licenses the historical evidence
only. Adoption would need its own addendum.

## A1.4 The principal angles in the parent were measured with an ill-conditioned formula

**This was found after the first draft of this addendum was frozen, and it corrects the parent's
central uncertainty claim.** `theta_C` was computed as `arccos` of the singular values of
`Q_a^H Q_b`. Near `sigma = 1` that is catastrophically ill-conditioned: for IDENTICAL subspaces it
returns 2.98e-08, and it cannot distinguish a true angle of 0 from 1e-8 or 1e-9, all of which read
about 3e-08. The stable form `‖(I - P_a) Q_b‖_2` resolves 1e-9 correctly and returns 7e-16 on
identical subspaces.

Measured on the real constructions both ways:

| target | arccos route, as recorded in the parent | sine route, correct | overstated by |
| --- | --- | --- | --- |
| `n = 12` | 2.581e-08 | 7.289e-11 | 354x |
| `n = 20` | 2.107e-08 | 8.086e-11 | 261x |

In exact mathematics the Reynolds-projector range and the invariant nullspace ARE the same
invariant subspace, so this is a NUMERICAL discrepancy between two constructions of one space and
not evidence that the exact spaces differ. That discrepancy is resolved at 7e-11, well above the
7e-16 identical floor, but roughly 350 times smaller than the parent records.

**Consequence, frozen.** `theta_C` and `theta_Q` are measured by the SINE form throughout, with

    s = sin(theta_max) = ‖(I - P_a) Q_b‖_2

the measured quantity. The arccos implementation may not be substituted anywhere in S1b.

**The half-angle recovery is frozen too, because prohibiting arccos leaves an implicit route back
to it.** `K_floor` needs `2 sin(theta/2)`, not `sin(theta)`. Frozen:

    d = 2 * sin( 0.5 * arcsin( clip(s, 0, 1) ) )

with the clip for roundoff safety only. `arcsin` is well conditioned at 0, which is where S1b's
angles live, exactly as `arccos` is not at 1. **The naive identity
`2 sin(theta/2) = sqrt(2 (1 - sqrt(1 - s^2)))` is PROHIBITED**: it cancels catastrophically at the
scales relevant here. It returns EXACTLY 0.0 below about 7.45e-09, which includes BOTH measured
`theta_C` values, 7.289e-11 and 8.086e-11, so it would silently annihilate the discrepancy term.
Just above that boundary it is already severely quantized rather than accurate: at `s = 1e-8` it
returns 1.4901161193847656e-08, an error of 49 percent, because `1 - s^2/2` stops resolving near
`sqrt(eps) = 1.49e-08`. An earlier draft of this addendum claimed exact zero "for every `s` at or
below 1e-8", which is false at 1e-8 itself; the boundary is 7.45e-09 and the regime just above it
is quantized, not correct. The frozen form is exact to machine precision
across `s` from 1e-16 to 0.5, verified. This is the same defect class as the arccos finding: a
formula that manufactures or destroys structure at the scales S1b operates at.

**The parent's claim that "subspace-construction ambiguity is the DOMINANT uncertainty, not the
arithmetic and not `kappa`" is WITHDRAWN. It is NOT replaced by a new dominance claim**, because
neither competing term is yet known:

- **the `kappa` term uses the MEASURED `max_r kappa(W_n^r)`, not the `1e6` admissibility ceiling.**
  A previous draft used the ceiling and got 7.18e-06 at both targets. The qualification measured
  `kappa(W_12) = 37.1` and `kappa(W_20) = 628`, giving `10 eps kappa ‖B‖_2` of 2.66e-10 and
  4.51e-09. The ceiling overstated that term by four to five orders of magnitude;
- **the discrepancy term uses `theta_Q`, not `theta_C`**, and `theta_Q` has not been measured under
  the sine form. The corrected coefficient-space values are not the operative input.

| target | arithmetic `100 eps ‖A‖_2` | `kappa` term, measured | `theta_C` arccos, withdrawn | `theta_C` sine, correct |
| --- | --- | --- | --- | --- |
| `n = 12` | 3.73e-12 | 2.66e-10 | 2.581e-08 | 7.289e-11 |
| `n = 20` | 9.77e-12 | 4.51e-09 | 2.107e-08 | 8.086e-11 |

**Dominance among the three `K_floor` terms is recomputed during the narrow rerun and is not
prejudged here.** For scale only, adjudicating nothing: a discrepancy term built from the corrected
`theta_C` would be 2.36e-07 and 2.61e-07, exceeding the measured `kappa` terms by roughly 885x and
58x. That is the opposite of a previous draft's conclusion, and it is exactly why no dominance
claim is frozen until `theta_Q` is in hand.

No gate outcome changes under any of this: G-SUBSPACE passes at `1e-6` on either angle measurement,
and the tilt arms discriminate either way.

The G-SUBSPACE gate threshold stays at `sin(theta_C) <= 1e-6`, and its arms are re-run under the
sine form: the `1e-4` tilt must go RED and the `1e-7` tilt must stay green. Under the arccos form
the `1e-7` control read 9.541e-08 against an identical-subspace floor of about 3e-08, so it was
numerically distinguishable from zero but badly conditioned and biased: it did not provide a
clean near-zero calibration. Under the sine form it reads 1.000e-07 and does.

## A1.5 What must be rerun

- ALL FOUR ladder controls of A1.1: the `k = 2` exact-zero control, the `k = 3` refusal control,
  the synthetic rule-4 triplet, and the rule-5 fall-through specimen;
- any synthetic routing case that consumes a ladder verdict;
- G-SVD-BRIDGE at both targets, which requires the shipped call once each, plus its four arms;
  on the fallback branch, the downstream SVD-route surfaces named in A1.3;
- every `theta_C` and `theta_Q` measurement, under the sine form, and the G-SUBSPACE tilt arms
with it;
- `K_floor` recomputed with the corrected discrepancy term.

Everything else stands on the hash-pinned room and the prior outputs, CONDITIONAL on G-SVD-BRIDGE
passing. If it fails, the `C^svd`-dependent surfaces named in A1.3 are rerun as well.

Still no live target. Nothing here authorizes the S1b run.
<!-- ADDENDUM-BOUNDARY -->

**Freeze record, addendum 1.** SHA-256 covers every byte ABOVE the boundary comment: `6da36a1c672772e8c731d4d66f20a77e81779d19d5cb883e617c028d43f46746`

```bash
sed '/^<!-- ADDENDUM-BOUNDARY -->$/,$d' S1B_ADDENDUM_1.md | shasum -a 256
```

The parent rule's freeze record is untouched and verifies independently.
