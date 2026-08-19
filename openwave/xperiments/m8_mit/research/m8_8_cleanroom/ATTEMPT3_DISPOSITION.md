# M8.8 clean-room attempt 3: RETIRED PRE-REVEAL

> **Disposition, 2026-08-19.** The run completed, is computationally the strongest of the
> three, and is retired before § 8 step 6 on two evidentiary defects. No answer material was
> opened, the firewall held, and the outputs are preserved unaltered at `f82ef464`.

## What this attempt got right, and it is most of it

| Item | Evidence |
| --- | --- |
| Manifest before production, with the freeze as an explicit event | validation artifact at 11:44, manifest declared final at 11:52 carrying the required status line verbatim, production implementation at 12:27 |
| The validation boundary held in substance | the validation script performs no file writes at all and states in its own header that it does not compute the result or populate the raw output; all outputs are written by the production implementation alone |
| Neither earlier failure mode recurred | the manifest was never amended after production began, and no departure from it arose |
| Reproducibility | both scripts rerun from the committed bytes at exit 0, and all four output files regenerate byte-identical |
| Firewall | all four seeded inputs byte-unchanged at the end of the run |

Nothing here suggests the computed values are wrong.

## Defect 1: the validation phase did not validate the manifest's own construction

The revised instruction requires every manifest-defined construction or convention checkable
from the permitted inputs to undergo a target-blind internal-consistency check before the
manifest may be declared final.

The manifest defines a substantial construction of the objects the whole computation runs on.
Those definitions are checkable before any result is computed, by dimension, character,
irreducibility and adjacency arguments, and the previous retirement is the direct evidence
that they need to be, since that run's failure was exactly a defect of that kind surviving
into implementation.

The recorded validation covers packet hashes, group enumeration, relators, chain relations,
the augmentation, free ranks and one boundary correspondence. Verified directly: the
validation script never opens the manifest at all and contains no checking of the
construction class in question. The manifest was then declared final.

So the causal sequence was recorded impeccably while the phase it recorded did not reach the
load-bearing content it exists to cover. This cannot be repaired after the fact without
moving a manifest check from before the freeze to after it, which is the one move the
architecture forbids.

## Defect 2: the frozen mutation-coverage contract was not discharged

Protocol § 9 states it in bold: every gate carries a runnable mutation that must redden it,
with enforced coverage and a nonzero exit, and manual attestation is not accepted. § 5.5
requires the raw output to carry every pre-reveal gate by identifier with its outcome.

The manifest instantiates fifteen gates. Verified by enumerating the raw output's own gate
identifiers: seven carry a runnable mutation arm and the remainder do not. Positive checks
exist for several of the uncovered ones, but a positive check is not the mutation arm the
contract requires, and the run nonetheless completed and reported success.

The deeper problem is the second half of that sentence in § 9. Coverage was not enforced at
all: nothing in the run could fail on incomplete coverage, so the requirement was discharged
by a reader noticing rather than by the artifact refusing. That is the same defect class the
protocol names elsewhere, a check that cannot fail, applied one level up to the suite.

## The commissioner-side miss, recorded rather than smoothed over

The landing record for this attempt reported the mutation arms that existed and flagged gate
sufficiency as a review question rather than answering it. The contract is explicit and
machine-checkable, the artifact was in hand, and enumerating instantiated gates against
mutation arms takes one command. Flagging a checkable question as somebody else's is the
same error the protocol warns about in its own domain, and it is recorded here so the next
landing record does not repeat it.

## Where this sits in the failure tree

| Attempt | Failure mode |
| --- | --- |
| 1 | the manifest was changed to fit the implementation |
| 2 | the manifest stayed frozen and the implementation departed from it instead of stopping |
| 3 | ordering and freeze finally correct, but the validation phase did not cover the manifest's own construction, and the frozen all-gates mutation-coverage requirement was not instantiated or enforced |

The progression is upstream, which is the useful direction: the architecture is no longer in
question, only the completeness of what it requires.

## The repair for attempt 4: enforce coverage in the run, not in the review

Two generic requirements are added, and the existing architecture is untouched.

The manifest must carry a coverage table with one row per construction or convention it
defines, each row naming either the validation artifact and specific check that tests it or
stating that it is not checkable before implementation and why. The final-status line may not
be written while any row is uncovered, and a row whose check exists only as prose counts as
uncovered.

Before production output counts as complete, an automated coverage check must establish, for
every pre-reveal gate the manifest instantiates, that a gate result exists, that the declared
runnable mutation exists, and that its red outcome is recorded. Incomplete coverage must make
the run exit nonzero, and the run may not print an all-pass summary while any instantiated
gate lacks mutation evidence.

The point is that the run becomes incapable of reporting success while its declared coverage
is incomplete. No reviewer has to notice a missing row after a three-hour run.

## What attempt 4 does not receive

No mention of any earlier commission, and nothing from this attempt's mathematics: no
construction, no object names, no route, no gate identifiers, and no statement of which gates
were uncovered. The fresh implementer proves completeness from its own manifest.

## Preserved

| Object | State |
| --- | --- |
| Attempt-3 outputs | `f82ef464`, pushed, no pull request, no comparison ever run |
| Attempt-3 room | `~/Desktop/OpenWave/M8.8_ATTEMPT3_RETIRED/`, intact |
| Canonical answer packet | SEALED, never opened. The case remains unconsumed |
