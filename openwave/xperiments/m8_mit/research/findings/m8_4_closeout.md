# M8.4 closeout: the dynamics program closes UNRESOLVED, on two instrument failures

> **Closed 2026-08-26, unresolved.** P1A, the pre-target qualification phase, RAN and CLOSED. The
> nonlinear pilot is BLOCKED: the one dynamics substrate that was built did not qualify, and the
> alternate the M8.5 charter contemplated was never built as a dynamics substrate at all.
> **No target configuration was ever executed. No nontrivial sector was spent; all eight remain
> unspent.**
> Nothing here is a verdict on MIT dynamics, and nothing here attributes the numerical defect to
> a specific stage. The successor chassis is elected separately as M8.5-C, target-free, with its
> protocol pending.

## What this closeout does and does not say

**Says:** the M8.4 preregistration was FILED 2026-08-23, naming `M4L_Erho` as the target-bearing
family: NINE flat bundles `E_ρ`, one per irrep of 2I, of which eight are target-bearing and `E_R0`
is the mandatory null control. Its § 2 gives `M4_int`'s eight-slot physical-section comparison a
structural N/A, and the filing executes no `M4_int` dynamics; § 0 preserves `M4_int` as a scope and
control object, and says a later three-sector dynamical study, if still worth doing, gets its own
compact appendix. This program then failed to produce a numerical substrate trustworthy enough to
run `M4L_Erho`, so no target configuration was ever executed.

**Does NOT say**, and no reader may infer:

- that MIT dynamics is falsified, or that no field equation on S³/2I realizes the McKay structure.
  OQ1's TWISTED branch is untested; its NATIVE branch was already closed by the kinematic close,
  and neither of those is a dynamics verdict. "Untested" is not a promise: even a successful
  successor answers only the narrower question, per the ceiling at the end of this document;
- that any sector was scored, spent, or partially spent. **No scored observable of the
  preregistration's § 7 was estimated on any candidate.** P1A did ship per-sector numbers, the
  `‖P_spec‖` of 2.2 to 3.1 and the invariance residuals across all nine sectors, but those are
  free-limit estimator qualification, zero-credit and expressly permitted by § 8. The `M4_int`
  structural N/A and the kinematic close stand as they were; the eight target bundles are unspent;
- that the numerical defect is attributed. The S1b adjudication left base RBF-FD assembly versus
  scalar quotient reduction explicitly UNSEPARATED, and this closeout does not separate them
  either;
- that the structural results are affected. The kinematic close, the M8.1/M8.2/M8.3 certifications,
  M8.5-A/B, M8.8 and M8.9 stand unchanged.

**Closed, not abandoned.** The reopening path is pre-committed and named below.

## The history, in two chassis

**Chassis 1, RBF-FD on the M8.5-B quotient backend.** P1A
([closeout](m8_4_p1a_closeout.md), #468) ran as the pilot's chartered qualification phase and split
its verdict across nine rows and six dispositions. PASSED: P0 substrate qualification, the
invariant-subspace estimator P1A.0 to P1A.3, the fixed-`k = 110` refinement family, and cloud
admissibility with mutation discrimination. NEGATIVE: `M_h` self-adjointness convergence, the
disposition most directly about the property this closure turns on. FAILED GLOBALLY:
imaginary-contamination qualification. NO VERDICTS ISSUED: per-sector scientific eligibility.
UNINSPECTED, NOT REACHED: the P1A.5 manufactured label calibration. NOT AUTHORIZED: nonlinear
target execution, which is the row that says in the phase's own table that no target ever ran. P1A's own adjudication is that this failure is SUBSTANTIVE rather than
another defect of the measuring machinery. Imaginary contamination ran from `5.95e-14` on the
`R_0` control to `1.8e-01` on `R_2`, and a precision ladder holding `L_h` byte-fixed showed it
unchanged from float64 through 50-digit arithmetic: a property of the assembled matrix, not the
eigensolver. That is where the nonlinear pilot stopped, with no target spent.

**The localization, M8.9.** S1 ([note](../m8_9/s1_note.md), #472) closed INSTRUMENT DEFECT: the
60-seed `R_0` block does not reproduce its own analytic decomposition, so the global spectrum
cannot identify its continuum levels. S1b ([closeout](m8_9_s1b_closeout.md), #484) removed that
assumption by supplying the invariant subspaces independently of the operator, and adjudicated
`S1b-SPECTRAL`: the TRIVIAL fibre alone gives `J = 2.201` at `n = 12` and `19.254` at `n = 20`.
Nontrivial fibre transport is not necessary for the non-real action; the defect lies upstream, in
the base discretization or the scalar quotient reduction, not separated.

**Chassis 2 was never built.** The M8.5 charter reads "prototype both, choose on evidence" for
(a) a 2I-equivariant grid and (b) a spectral method in 2I-symmetric harmonics. Only (a) was built
as a dynamics substrate. Route (b) existed solely as character-averaging certification machinery
(`route_b_core.py`), never as a simulation backend. **The evidence for route (b) as dynamics was
never gathered.**

## Why the program closes here rather than repairing

P1A tried the floor and the refinement family; S1 tried the spectrum; S1b removed the last
confound. The `M_h`-adjoint defect is O(1), the same scale as the operator, and survives 50-digit
arithmetic. No qualified refinement demonstrated recovery of the self-adjoint substrate the M8.4
observables require. What the two closeouts do support is a partial localization: the defect sits
in the assembled matrix rather than in the eigensolver or the estimator, and within that matrix it
is unseparated between the base discretization and the scalar quotient reduction. No localized
repair target inside the retired RBF-FD lineage was identified; that is narrower than saying no
reformulation could exist. Another rung of seeds or precision is not the answer, and the chassis fork
memo records the evidence bar under which the RBF-FD route could ever return.

**M8.9's deferred COMMISSIONING DECISION is answered here, since it is retired rather than left
dangling.** The S1b closeout said whether the remaining base-versus-reduction distinction was worth
a separately commissioned task was "a decision this closeout does not make." It is made now: **it
is not commissioned.** The technical disjunction itself remains UNANSWERED, and this closeout does
not answer it. That localization carries decision value only while the RBF-FD lineage is a
candidate for M8.4 dynamics, which it is not unless the chassis memo's two-condition readmission
bar (structure AND consistency, separately frozen, pre-target) is met. The door is left ajar on
those terms; nothing is pending behind it.

## What was learned, and it is not nothing

The negatives are narrowing rather than empty:

- the invariant-subspace estimator qualified and is reusable: `M_h`-orthogonal projectors with
  idempotence 2.4e-16 to 6.0e-16 across all nine sectors, a verified three-way leakage identity, and `‖P_spec‖` of 2.2 to 3.1 on the
  target clusters, so strong non-normality did NOT make the scored subspaces ill-conditioned,
  which was the central worry;
- the S1b instrument, its `K`/`J` witnesses, precision ladder and dual-construction agreement, is
  operator-agnostic and carries forward to any successor chassis;
- three numerical-formula and control-specification defects were found and fixed before any
  instrument met a target, each recorded in the S1b addenda;
- the kinematic close stands: native single-valued fields on S³/2I carry no nontrivial McKay slot
  at any level, closing OQ1's native branch by theorem with no dynamics run.

## Status changes this closeout makes

| Item | Before | After |
| --- | --- | --- |
| M8.4 row | BACKLOG, nonlinear pilot blocked | DONE, closed unresolved (author's call) |
| M8.9's deferred COMMISSIONING decision | "a decision this closeout does not make" | not commissioned; the technical base-vs-reduction question stays open and unanswered |
| M8.7 gate CONDITION | a validated in-platform dynamics | UNCHANGED, not touched by this closeout |
| M8.7 gate ROUTING PROSE | pointed at an M8.4-lineage result | repointed: M8.5-C, then a fresh preregistration, then a validated dynamics |
| MODELS.md dynamics rows | verdicts `not yet tested`, prose naming M8.4 as the live program | verdicts UNCHANGED; prose repointed, since the program named there has closed |
| The eight target sectors | unspent | unspent |

MODELS.md's verdict cells are untouched: every affected row already read `🚧 not yet tested` and
still does, because nothing was claimed and nothing was tested. What changes is prose that
described a runnable task, which is the second question the #375 propagation rule asks of every
document and the one that caught two docs last time.

## The reopening path, pre-committed

The chassis decision is recorded in the M8.4 chassis fork memo, an author-side working document
held outside the repository, frozen at

    44c664d1ceb17949da78e55dcc5fe322cd447375bfe91ab3d33256775f386f4b

which covers every byte above the memo's boundary comment and is checked with

    sed '/^<!-- MEMO-BOUNDARY -->$/,$d' CHASSIS_DECISION_MEMO.md | shasum -a 256

The boundary convention is load-bearing: a range that includes the marker line yields a different
digest. Text may be appended below the boundary without disturbing the freeze.

**THREE errata are filed below that boundary, and this closeout inherits the corrected text, not
the frozen text.** (1) The aperture table closes projector tolerance and oversampling by "gate 4,
gate 11a", but the adopted decision is 11(b), which makes 11(a) inoperative; gate 4 is unaffected
and still closes the aperture. (2) The Option C rationale reads "Two chassis failed", which
contradicts the memo's own inputs table; corrected to the proposition used throughout this
document, that the one chassis built failed qualification, the alternate was never built as a
dynamics substrate, and no target sector was spent. (3) The execution sequence reads "MODELS.md
untouched", and it was not: two prose lines named M8.4 as the live program and both were repointed.
The narrower true statement is that its VERDICT CELLS are untouched, both already reading `🚧 not
yet tested` and still doing so.

None of the three changes the decision, none touches the hashed region, and the digest above
recomputes identical with all three present. All three are filed where the claims live rather than
handled as commentary elsewhere, which is the same rule the #375 sweep applies to this column's own
documents.

**Nothing checkable in this closeout rests on the memo.** It is cited for three things only: the
decision it records, the readmission bar it fixes, and its own hash as a provenance pointer. Every
measured quantity in this document is sourced to a shipped artifact in this repository. The memo
does carry verified computations of its own, the degree-`4N` exactness result with its mutation arm
and the per-sector mode-count table; those are design inputs for M8.5-C and ship with the M8.5-C
protocol, not from here. This follows P1A's practice of citing an author-side commission for a
procedural disposition while stating plainly that nothing numerical rests on it.

**The pre-redline draft is preserved, not discarded.** This branch was rebuilt before the PR
opened, so its history reads as the corrected record rather than as a false statement followed by a
correction living somewhere else. A commit subject cannot carry an erratum the way this memo can:
there is no place beside the claim to put one, and when co-location is impossible the honest option
is removal. The four-commit draft, including the wrong propositions the redlines caught, is kept
unmerged at `m8_4-closeout-predraft`, `bbdc1c83b01556f242e3467f07c8db8bb851102e`. That is the shape
M8.8 used for its attempt 1 and M8.5-B for `M85B-ADJ-04`: retire the artifact on the record, keep
it locatable, do not let it into the official run.

**Option C, cubic only.** The spectral/Galerkin chassis is elected as **M8.5-C**, target-free
simulation engineering under the charter that already names route (b). Its protocol is a separate
artifact and is NOT filed as of this closeout, so M8.5-C is an elected successor with its
preregistration pending, not a running task; the roadmap row says so. Its room contains no M8.4
target object and its qualification executable cannot load one; work at zero amplitude on any
`E_ρ`, all work on `E_R0`, and manufactured operators are qualification. Two terminal outcomes,
no repair round after execution begins:

- **`M8.5-C-QUALIFIED`** plus a frozen cost estimate under the named ceiling files a fresh
  successor preregistration under a NEW identity, inheriting a chassis that was never in a room
  with the target;
- **`M8.5-C-FAILED`** publishes the failure and closes the spectral route; the M8.5 row then
  records the grid backend retained for what `M85B-ADJ-07` certified and excluded from
  M8.4-lineage dynamics, and the reopening path closes with it.

**The chain to M8.7 has three links, and M8.5-C is only the first.** M8.5-C is target-free
engineering and cannot by itself validate a dynamics. The precommitted route is: M8.5-C reaches
`M8.5-C-QUALIFIED`; THEN a fresh preregistration under a new identity runs a target; THEN that
validated in-platform dynamics opens M8.7's gate. A chassis qualification alone does not open it.

A future M8.4-lineage question is therefore a fresh preregistration, not a resumption of this row.
That is the point: the instrument gets built and judged somewhere the experiment cannot reach.

## What a successor would inherit, and what it would not

**Inherits:** the M8.2 lock and the M8.4 preregistration's frozen § 7 observable and § 9
anti-tuning clauses; the certified fibre realizations; the invariant-basis machinery; the S1b
qualification instrument; and the pinned equation of motion, read from `m4_ewt/wave_engine.py`
at `c9dc3796`: second order in time, leapfrog, `∂²ψ/∂t² = c²∇²ψ − dV(ψ)`, with `cubic_nls` being
`v_mode 1`, `V = (c1/4)u²` on `u = ⟨ψ,ψ⟩`.

**Does not inherit:** the RBF-FD dynamics substrate; the `saturating` configuration, dropped from
the successor's scope by the frozen decision, which narrows the experiment matrix relative to
M8.4 § 9 and must be reported as such; and any claim to have tested OQ1's dynamical branch.

**And a ceiling that survives the closure.** Should a successor ever run, `M4L_Erho` installs one
flat bundle per McKay slot BY CONSTRUCTION. No result from it can be reported as "the dynamics
selected the McKay slots," the free eigenvalues earn zero credit, and small-amplitude continuation
is pre-labelled "nonlinear persistence of the installed free structure" by the preregistration's
own § 1. That ceiling is stated here so it survives the gap between this closeout and any
reopening.
