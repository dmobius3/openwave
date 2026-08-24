# M8.4 P1A closeout: the subspace estimator qualifies, this discrete operator does not supply usable spectra

> **Phase closed 2026-08-24.** P1A was chartered to answer one question before any threshold,
> amplitude ladder, or coupling was chosen: what do "angle" and "leakage" physically mean for the
> discrete operator, and can they be measured reliably? It answered both halves. The machinery
> qualifies. The operator does not supply spectra controlled enough for M8.4 to use it.
>
> No nonlinear target configuration was ever run. No target sector was spent.

## Status

| Component | Final P1A status |
| --- | --- |
| P0 substrate qualification | **PASS** |
| P1A.0 to P1A.3, invariant-subspace estimator | **PASS** |
| `M_h` self-adjointness convergence | **NEGATIVE** |
| fixed-`k = 110` refinement family | **PASS** |
| cloud admissibility and mutation discrimination | **PASS** |
| imaginary-contamination qualification | **FAIL globally** |
| per-sector scientific eligibility | **NO VERDICTS ISSUED** |
| P1A.5 manufactured label calibration | **UNINSPECTED, NOT REACHED** |
| nonlinear target execution | **NOT AUTHORIZED** |

## What qualified, and it is the hard part

The invariant-subspace machinery works, and survived adversarial checking.

| Check | Result |
| --- | --- |
| projector idempotence `‖P² − P‖` | 2.4e-16 to 6.0e-16, all nine sectors |
| `M_h`-symmetry `‖P†M_h − M_hP‖` | 1.2e-16 to 1.8e-16 |
| invariance residual `‖(I − P_Q)LQ‖_{M_h}` | 3.9e-14 to 9.7e-14 |
| Schur against Riesz, `θ_max`, contour preregistered from analytic levels | 0 to 3.9e-08 |
| three-way leakage identity, projector against angle against overlap | max disagreement 8.1e-16 |
| `‖P_spec‖` on the target clusters | **2.2 to 3.1** |

That last row matters most. Despite a Henrici departure of 0.63, the invariant subspaces M8.4
scores are well conditioned. Strong non-normality of the operator did NOT translate into
sensitivity of the subspaces, which was the central worry the phase existed to settle.

The supporting apparatus also qualified: the `M_h` metric passes its geometry-only gates with
moment reproduction to 1e-11; the fixed-`k` refinement family is a genuine one-parameter ladder;
and the cloud admissibility gate is demonstrably discriminating, rejecting a manufactured cloud
whose mesh ratio is driven to 83.4 against a frozen ceiling of 20.0.

## What failed, and why it is substantive

Imaginary contamination. The continuum operator is self-adjoint, so any imaginary part in the
computed spectrum is error. Measured on the target clusters at the production cloud, it ranges
from `9e-15` on `R_0` to `1.8e-01` on `R_2`, and the pooled envelope needed to cover that spread
is about 58,688 times its own fitted trend. Against the power gate frozen BEFORE the fit,
`E(h) < 10·T(h)`, that fails, and a rule with no power issues no verdicts.

**The failure is substantive, not another instrument defect.** That distinction was established
by a standalone localization test, run commissioner-side rather than by the qualification unit,
because it selects no threshold and its outcome cannot alter any P1A gate. With the assembled
`L_h` held byte-for-byte fixed and only the arithmetic precision varied, the target imaginary
parts are unchanged from float64 through 30- and 50-digit arithmetic.

| sector | float64 | 30 dps | 50 dps |
| --- | --- | --- | --- |
| `R_1` | 1.021e-10 | 1.021e-10 | 1.021e-10 |
| `R_3` | 1.942e-10 | 1.942e-10 | 1.942e-10 |
| `R_2` | 2.222e-01 | 2.222e-01 | 2.222e-01 |

Identical to four significant figures across sixteen to fifty digits. **The imaginary spectrum is
a property of the assembled discrete matrix; its source lies upstream of the eigensolver, in the
assembly or the discretization.** The test separates solver from matrix. It does NOT separate
equivariant assembly from the underlying RBF-FD discretization, and this closeout makes no claim
about which is responsible.

Two consequences follow, neither discretionary. The observed complex parts cannot be reclassified
as solver noise by any numerical floor, however closely an empirical formula happens to fit them.
And a further Arm-A repair would therefore be chasing something that is not there.

## Three qualification rules failed before this one, and each failure was informative

Recorded so the sequence is legible rather than looking like thrashing. Each was caught by
adversarial adjudication, and every one of them before a target sector was spent.

| Attempt | Defect | What it taught |
| --- | --- | --- |
| per-sector power-law envelope | nine independent fits; the dirtiest sector passed and two of the cleanest failed | fit quality had replaced contamination magnitude as the effective test |
| pooled envelope on the original ladder | seed count did not control geometric resolution, and the stencil rule changed the discrete scheme across most of the sequence | it was not a qualified ONE-PARAMETER refinement family; fixing `k = 110` repaired it |
| eigensolver-residual floor | the solver is accurate to 3e-12; the floor moved 3x and reclassified nothing | the contamination is not solver error, which the localization test then confirmed causally |

## What this does and does not mean

**It does not erase the positive result.** The invariant-subspace machinery genuinely qualified,
under adversarial checks it could have failed and did not. What failed is a separate question:
whether THIS discrete operator supplies spectra controlled enough for M8.4 to use that machinery
scientifically. The answer is no, on the frozen rules, and the localization test establishes that
the answer cannot be changed by treating the observed complex parts as solver noise.

**It is a preregistered outcome, not a breakdown.** The contract's § 10 already provided that
numerical inability to resolve the cluster yields estimator-qualification failure and NO
scientific label, never `destroy`. That rule was frozen before any manufactured case was built,
precisely so solver pathology could not become a physics result. It is now doing its job.

**M8.4's claim ceiling is untouched**, because nothing was claimed. No sector verdicts issued, no
labels assigned, P1A.5 never opened.

## The boundary

Any continuation would require a separately commissioned, independently justified structural
discretization or assembly correction, whose purpose is to recover the continuum self-adjoint
structure to a precision adequate for the frozen M8.4 observables. **That work is not a repair to
P1A and does not reopen its outcome.**

Establishing which upstream stage carries the defect, equivariant assembly or the underlying
RBF-FD discretization, is the natural first step of any such effort, and is not answered here.

## Artefacts, in this repository

| Path | Contents |
| --- | --- |
| [`../m8_4_p1a/p0/`](../m8_4_p1a/p0/) | the P0 substrate-qualification package, with `frozen_tolerances.py` and `regression_tests.py` |
| [`../m8_4_p1a/p1a/`](../m8_4_p1a/p1a/) | the P1A packages: metric, diagnostics, subspace extractor, ladder, cloud gate, contamination rules, floor repair |
| [`../m8_4_p1a/notes/`](../m8_4_p1a/notes/) | the three qualification notes: P0, P1A.4a, P1A.4b, each carrying its frozen constants and hashes |
| [`../m8_4_p1a/localization/`](../m8_4_p1a/localization/) | the localization finding and its runnable reproducer |

**Reproducing the localization test**, which is the claim-bearing piece of § "What failed":

```bash
cd openwave/xperiments/m8_mit/research/m8_4_p1a
PYTHONPATH=.:../m8_5b/pilot:../m8_5b/production python3 localization/localization_test.py
```

It carries a deliberately precision-sensitive control: the companion matrix of `(x-1)^10`, whose
exact spectrum is real, so any computed imaginary component there is numerical error from the
finite-precision eigensolve. The ladder must shrink it, and does, from `4.79e-02` to `1.17e-05`. Without that arm a persistence
result would be indistinguishable from a broken ladder. The script exits nonzero if the arm
fails.

Commissioner adjudication history, the verification notes behind each of the three failed
qualification rules, is an external author-side record and is deliberately NOT part of this
evidence package. Nothing in this closeout depends on it.

P0's bytes were verified intact on every P1A run, 10 of 10, and the ladder constants frozen
before the first contamination value was computed were unchanged at every subsequent run,
including across an interrupted session.
