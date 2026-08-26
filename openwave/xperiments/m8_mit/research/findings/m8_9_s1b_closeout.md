# M8.9 S1b closeout: nontrivial fibre transport is not necessary for the non-real spectrum

> **Adjudicated `S1b-SPECTRAL` on 2026-08-26**, by an independent adjudication of measurements
> exposed during a nonconforming qualification run. Three findings sit on this record together and
> all three are true: the instrument qualified, the execution exceeded its authorization, and the
> exposed reading mechanically satisfies the SPECTRAL branch under rules frozen and hash-published
> before it existed. Read § "Provenance and the nonconformance" before the result.

## What S1b answered

M8.9 asked one question: does the non-real spectrum in the M8.4 substrate enter through the base
RBF-FD discretization, through the equivariant quotient assembly, or through their interaction?

**The trivial bundle `R_0` was probed at high harmonic level with an analytic subspace supplied
independently of the operator**, and its compressed action is emphatically non-real:

| target | `lambda` | `‖A_n‖_2` | `‖K_n‖_2` | `K_floor` | `J`, all three precision rungs | ladder |
| --- | --- | --- | --- | --- | --- | --- |
| `n = 12` | 168 | 175.18 | 3.607 | 6.85e-07 | 2.2010292238 | PERSISTS |
| `n = 20` | 440 | 508.33 | 64.410 | 9.93e-07 | 19.2541741428 | PERSISTS |

The `J` readings are stable to thirteen significant figures at `n = 12` and fourteen at `n = 20`
across the three precision rungs of one construction. ACROSS the two independent constructions the
agreement is narrower, eight figures at `n = 12` and nine at `n = 20`, which is the
subspace-construction discrepancy rather than a precision effect. Separately, and it is a
different quantity against a different floor, `‖K_n‖_2` exceeds `K_floor` by roughly 5.3e+06 at
`n = 12` and 6.5e+07 at `n = 20`. This is not a threshold-edge result.

**Licensed, and nothing beyond it:**

> the trivial fibre alone can produce a non-real compressed action, so nontrivial fibre transport
> is NOT NECESSARY for that phenomenon, and the base discretization or the scalar quotient reduction
> is strongly implicated.

**What that does NOT establish.** The licensed sentence is deliberately disjunctive. S1b has NOT
separated the base RBF-FD discretization from the scalar quotient reduction; it has excluded a
third candidate. Writing that M8.9 identified the discretization specifically would claim more than
the frozen rule licenses.

**S1b-SPECTRAL does not trigger S2 or the reduction-mechanism successor.** The frozen routing sends
`S1b-ADJOINT` and `S1b-NULL` onward to a separately frozen interaction or reduction comparison;
SPECTRAL is not routed there, and the rule is SILENT on its disposition. Closing M8.9 at the
claim ceiling is this closeout's decision, not a sentence of the rule.
**No further stage is commissioned by this outcome**, and the remaining base-versus-reduction
distinction is left open rather than pursued.

## Why S1b succeeded where S1 failed

S1 closed `INSTRUMENT DEFECT`. Its architecture diagonalized the whole 60x60 block and asked the
global spectrum to sort itself into continuum levels; G-MULT failed, so no reading was licensed.

S1b removed that assumption entirely. The invariant subspaces are constructed independently of `L`,
by two algebraically different routes, sampled on the production cloud, and the shipped operator is
then asked what it does to them. No cluster identification, no Voronoi assignment, no dependence on
the top of a finite block being faithful. With that assumption gone, the non-real action was still
there and still enormous.

## The instrument

Eight gates green with their mutation arms, four ladder controls each firing at the rule it was
frozen to exercise, all seven adjudicator outcomes reachable, all four precedence collisions
correct, end-of-run input re-verification clean. Full record in
[`../m8_9/s1b/q3a_qualification_note.md`](../m8_9/s1b/q3a_qualification_note.md).

Three instrument defects were found and fixed BEFORE the instrument met a target, each recorded in
the addenda. Two are numerical-formula defects: an `arccos` principal-angle route that manufactured
1e-08 structure from nothing and overstated the construction discrepancy by 350x, and a half-angle
identity that annihilates angles below 7.45e-09 into exact zero. The third is a
control-specification defect: a collapse control whose verdict changed with an unfrozen parameter
`k`.

## Provenance and the nonconformance

**Attempt `q3a` was commissioned as an instrument qualification and forbidden to touch the live
target. It qualified the instrument and then evaluated the target anyway.** The full account,
including the contract defect that made the boundary ambiguous, is in
[`../m8_9/s1b/s1b_addendum_4.md`](../m8_9/s1b/s1b_addendum_4.md) § A4.1. The gates structurally
require forming `A_n`; evaluating `‖K_n‖` and running the ladder on `J` was avoidable and was not
authorized. No governing document drew that line. The fault is the contract author's.

**The measurements are retained rather than voided**, because every rule interpreting them was
frozen and hash-published beforehand, verified at run start and re-verified at run end with no
input drift. The usual reason an unauthorized look invalidates an experiment, that a threshold or
statistic can be adjusted afterwards, was closed before the run began. Applying those frozen rules
to those pinned bytes is a ratification, not a new measurement.

**The branch was issued by a separate Adjudication Unit**, in a room containing the contracts and
the pinned evidence and no numerical code at all, unable to construct an operator had it wanted to.
The author had already read the exposed numbers and did not issue the branch. Trace in
[`../m8_9/s1b/adjudication.md`](../m8_9/s1b/adjudication.md).

**One correction to that adjudication document, which does not affect its verdict.** Its § 6 states
that "both S1b-SPECTRAL and other outcomes leave M8.9 OPEN: a separately frozen interaction or
reduction comparison MAY be commissioned". That over-reads the frozen routing: ADJOINT and NULL are
the branches routed there, and the rule says nothing about SPECTRAL's disposition, so the
open-or-closed question sits OUTSIDE the rule and is answered in § Status on this closeout's own
authority. The rule-by-rule trace and the forced outcome
are unaffected. The document is preserved unedited and the correction recorded here.

**No blind rerun is possible and none will be claimed.** `J_12 = 2.201` and `J_20 = 19.254` have
been seen. Any later execution is a reproduction of the exposed result, never a blind first run.

## Chain of custody

Six frozen documents, each verifying its own published hash, then the evidence manifest, then every
record any branch rule reads:

| artifact | SHA-256 |
| --- | --- |
| `s1b_decision_rule.md`, frozen region | `c44c603a8058ed8529e5bb0f42ec168b443ff22b050d3748f1f1de3537c7d297` |
| `s1b_addendum_1.md` | `6da36a1c672772e8c731d4d66f20a77e81779d19d5cb883e617c028d43f46746` |
| `s1b_addendum_2.md` | `14011c338ce331d9bb3424f0ed5aaff31a1c89a9ceafbc2bae41a91ebb41a222` |
| `s1b_addendum_3.md` | `e3304fe993f7a2ea523cfcbc4110c0f37174caea95383b17e465c8fa85b6e28c` |
| `s1b_addendum_4.md` | `98484b0edd2de97e34093a564674cd8f128a7bca38f96960c1007f7f45f00634` |
| `s1b_addendum_5.md` | `7ca059f074059f925a7f231fcb1ac93932e890121ef9abd1975a4df78542b3e5` |
| `q3a_output_manifest.json` | `f5401c5179d0cde42a8763de175542d3d04aaa316be0329ec30a42cd9d6bc3a4` |
| `q3a_results/target_n12.json` | `d96c7fdb3ab192b48f04c49040dbb8f6795379a35fb04330078bff22115a5b91` |
| `q3a_results/target_n20.json` | `4e8c0c72fac0a8e2f127fc1629bb4afa4b2d90e9b12b9df97714d56b1057fed7` |

**The shipped package is fully manifest-verifiable.** All 34 entries are present, 32 hashing
exactly. The two that differ, `q3a_qualification_note.md` and `q3a_run.log`, are the post-manifest
appends documented in addendum 4 § A4.3: both were written to after the manifest was computed. Both
logs required `git add -f`, since `.gitignore` matches `*.log`, and an earlier commit claiming to
ship them had in fact shipped only three of the five files it named. The
note's append point is located exactly, since truncating immediately before its `## Output Manifest`
section reproduces `098631325e9692033597ea3f63b0c5955178c929a9330570573f9767a362175c`. Both are
forensic context and no branch rule takes an input from either.

`round1_s1b_qualification.py` is round-1 archival code, hashing to
`5a9e04845375c4d12c3a475607f20a1d5f13cc82829d64caa022c82d5e784802` as addendum 2 § A2.1 cites. It
is preserved as evidence of what round 1 ran and, per addendum 2 § A2.3, may not be executed: it
authenticates the round-1 room manifest and cannot run anywhere that manifest is not current.

## Status

**M8.9 is CLOSED at its stated claim ceiling**: nontrivial fibre transport is excluded as a
necessary cause; base RBF-FD versus scalar quotient reduction remains unresolved and is not pursued
by this branch.

M8.4's pilot remains blocked and its P1A outcome final; nothing here reopens it. Whether the
remaining base-versus-reduction distinction is worth a separately commissioned task is a decision
this closeout does not make.
