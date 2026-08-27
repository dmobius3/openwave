# M8.4 closeout: the dynamics program closes UNRESOLVED, on two instrument failures

> **Closed 2026-08-26, unresolved.** Two numerical substrates were built and neither qualified.
> **No M8.4 dynamics were ever run. No nontrivial sector was spent; all eight remain unspent.**
> Nothing here is a verdict on MIT dynamics, and nothing here attributes the numerical defect to
> a specific stage. The successor chassis is commissioned separately as M8.5-C, target-free.

## What this closeout does and does not say

**Says:** this program did not produce a numerical substrate trustworthy enough to run the
`M4L_Erho` experiment, so the experiment was never run.

**Does NOT say**, and no reader may infer:

- that MIT dynamics is falsified, or that no field equation on S³/2I realizes the McKay structure.
  OQ1 is untouched; it was never tested;
- that any sector was scored, spent, or partially spent. The `M4_int` structural N/A and the
  kinematic close stand as they were; the eight target bundles are unspent;
- that the numerical defect is attributed. The S1b adjudication left base RBF-FD assembly versus
  scalar quotient reduction explicitly UNSEPARATED, and this closeout does not separate them
  either;
- that the structural results are affected. The kinematic close, the M8.1/M8.2/M8.3 certifications,
  M8.5-A/B, M8.8 and M8.9 stand unchanged.

**Closed, not abandoned.** The reopening path is pre-committed and named below.

## The history, in two chassis

**Chassis 1, RBF-FD on the M8.5-B quotient backend.** P1A
([closeout](m8_4_p1a_closeout.md), #468) qualified the invariant-subspace estimator and FAILED the
substrate's spectral qualification globally. Imaginary contamination ran from `5.95e-14` on the
`R_0` control to `1.8e-01` on `R_2`, and a precision ladder holding `L_h` byte-fixed showed it
unchanged from float64 through 50-digit arithmetic: a property of the assembled matrix, not the
eigensolver. The pilot was blocked with no target spent.

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
observables require. Another rung of seeds or precision is not the answer, and the chassis fork
memo records the evidence bar under which the RBF-FD route could ever return.

**M8.9's deferred question is answered here, since it is retired rather than left dangling.** The
S1b closeout said whether the remaining base-versus-reduction distinction was worth a separately
commissioned task was "a decision this closeout does not make." It is made now: **it is not
commissioned.** That localization has decision value only if the RBF-FD lineage is retained for
M8.4 dynamics, and it is not. The disjunction stays open on the record as an unanswered question,
not as pending work.

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
| M8.4 row | BACKLOG, pilot blocked | DONE, closed unresolved |
| M8.9's deferred separation question | "a decision this closeout does not make" | not commissioned; retired with the lineage |
| M8.7 gate | gated on an M8.4-lineage result | gated still, now via the M8.5-C reopening path rather than a closed lineage |
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

**Option C, cubic only.** The spectral/Galerkin chassis is commissioned as **M8.5-C**, target-free
simulation engineering under the charter that already names route (b). Its room contains no M8.4
target object and its qualification executable cannot load one; work at zero amplitude on any
`E_ρ`, all work on `E_R0`, and manufactured operators are qualification. Two terminal outcomes,
no repair round after execution begins:

- **`M8.5-C-QUALIFIED`** plus a frozen cost estimate under the named ceiling files a fresh
  successor preregistration under a NEW identity, inheriting a chassis that was never in a room
  with the target;
- **`M8.5-C-FAILED`** publishes the failure and closes the spectral route; the M8.5 row then
  records the grid backend retained for what `M85B-ADJ-07` certified and excluded from
  M8.4-lineage dynamics, and the reopening path closes with it.

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
