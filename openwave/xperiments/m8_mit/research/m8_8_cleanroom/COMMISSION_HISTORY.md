# M8.8 clean room: the commission history

> **What this file is.** The canonical record of the five clean-room commissions run under
> § 8 of the reproduction protocol. It exists so that one file answers three questions
> without spelunking through commits: why five runs were made, why the first four were
> retired, and why the fifth is held as the Phase A candidate.
>
> **The answer packet has never been opened.** No comparison has run against any attempt, no
> orientation has been selected, and no verdict exists. The four retirements below are all
> pre-reveal, and the fifth commission is pre-reveal on qualification hold, so the case is
> unconsumed and a conforming run remains fully available.
>
> Per-file hashes are re-derivable from the tree at each landed commit and are not restated
> here. The one hash inventory that is operative rather than descriptive, Phase A's thirteen
> frozen artifacts, lives in
> [Addendum 1](../findings/m8_8_reproduction_protocol.md) and is not duplicated.

## Why there were five

Each commission was a genuinely fresh context, launched into a rebuilt room holding the same
frozen technical inputs: the protocol at the content commit, the group packet, and the
construction packet, all three verified against the § 8 lock manifest before every launch. No
attempt was told that an earlier one existed, and no attempt received any artifact from a
predecessor.

Attempts 2 through 4 each received a revised `TASK.md` in response to the preceding
retirement, and those revisions are recorded per attempt below; each closed an ambiguity the
previous run had exposed. Attempt 5 deliberately reused attempt 4's instruction unchanged, to
test whether a fresh implementer could comply with the mature contract as written rather than
whether the contract could be tuned until something passed.

## The five commissions

### Attempt 1, retired: the manifest was back-fitted to the implementation

| Field | Value |
| --- | --- |
| Instruction received | `TASK.md` `c6132075…` |
| Outputs landed, unread | `e6c9a44b` |
| Disposition | RETIRED PRE-REVEAL, manifest-order nonconformance |

The method-and-gate manifest was written before any implementation code, as § 8 step 4
requires. It was then amended, after the implementation existed, to correct a stated `SU(2)`
embedding so that it matched what the code did. An embedding convention is the class of
methodological statement the manifest exists to commit in advance, so the record could no
longer distinguish a mis-transcription from a choice settled during implementation and
back-filled.

Explicitly not the defect: the firewall held, all four seeded inputs byte-unchanged; the
computation was never suspected wrong; a run-order dependency in the output was disclosed and
reproducible, with a leaf-level comparison finding 320 of 320 shared values equal.

Room modification times showed the manifest's final write postdating the implementation by
thirteen seconds, which is how the amendment surfaced. The lesson carried forward is that
mtimes alone are not sufficient ordering evidence: the causal order rested on the session
record, and later instructions required the freeze to be an explicit declared event rather
than something inferred from timestamps.

### Attempt 2, retired: the implementation departed from a frozen manifest instead of stopping

| Field | Value |
| --- | --- |
| Instruction received | `TASK.md` `08820050…`, adding: once implementation begins the manifest is immutable, and a mismatch is resolved by conforming the implementation or stopping, never by amending the manifest |
| Outputs landed, unread | `d3c7859a` |
| Disposition | RETIRED PRE-REVEAL, manifest-departure nonconformance |

The new clause worked: the implementer found a genuine mathematical error in its own frozen
manifest, did not amend it, and disclosed the mismatch in a dedicated artifact that quoted the
frozen statement and named what the code did instead. The manifest had declared an irreducible
as the Galois conjugate of a representation whose character is rational, so that conjugate
coincides with the original rather than being distinct.

Conforming was therefore unavailable, which left stopping as the only authorized branch. The
implementer reported and then continued with a replacement construction.

A commissioner-side argument that the repair was forced rather than chosen, since the group
has exactly two irreducibles of the relevant dimension, was verified as mathematics and
rejected as procedure: accepting it would have installed an unregistered rule that a departure
is permissible whenever the commissioner later judges it constrained enough, which is the
post-hoc adjudication a precommitted manifest exists to abolish.

Explicitly not the defect: the manifest was not amended, the firewall held, the diagnosis was
independently verified as correct, and the disclosure was exemplary.

### Attempt 3, retired: the validation phase did not cover the manifest's own construction

| Field | Value |
| --- | --- |
| Instruction received | `TASK.md` `75cef74d…`, adding pre-implementation manifest validation, the validation-artifact boundary, and an explicit machine-readable freeze line |
| Outputs landed, unread | `f82ef464` |
| Disposition | RETIRED PRE-REVEAL, two evidentiary defects |

Ordering and the freeze were impeccable, and the boundary held in substance: the validation
script wrote no files at all and the outputs came from production alone.

First defect: the validation phase never reached the manifest's own construction of the
objects the computation runs on. The validation script never opened the manifest and contained
no checking of that class, which is exactly the class whose defect had retired attempt 2, and
the manifest was declared final anyway.

Second defect: the frozen mutation-coverage contract was not discharged. Fifteen gates were
instantiated and seven carried mutation arms, and more importantly coverage was not enforced
at all, so nothing in the run could fail on incomplete coverage.

Commissioner-side process finding: the landing record for this attempt flagged gate
sufficiency as a review question instead of answering it, when the contract is explicit and
machine-checkable and the enumeration is one command. Later landings answer it directly.

### Attempt 4, retired: gate-contract violations against already-explicit requirements

| Field | Value |
| --- | --- |
| Instruction received | `TASK.md` `e3d9b908…`, adding the coverage table, ID-bearing registries, registry-to-coverage set equality, and a production coverage checker that must take the frozen registry as its authoritative expected set |
| Outputs landed, unread | `9aa740ea` |
| Disposition | RETIRED PRE-REVEAL, multiple gate-contract nonconformances |

The registry architecture worked: the manifest carried construction, convention and gate
registries, the validator genuinely read the manifest and proved exact set equality, the
freeze was explicit, and the production coverage check blocked with a nonzero exit.

The violations were of requirements the instruction and protocol already stated explicitly.
Pre-implementation exemptions were claimed for constructions the instruction expressly permits
a validation script to check. The production coverage checker hard-coded its expected gate set
instead of reading the frozen manifest, which is the duplicated-list implementation the
instruction prohibits by name. One mutation verdict was assigned as a literal with a comment
reading "by design", and one gate verdict was an attestation with no check behind it. The
exact saturation requirement was replaced with per-irrep acyclicity plus Euler-characteristic
reasoning on the stated grounds that the exact route is expensive, which is the substitution
§ 9 forbids and supplies a counterexample against, after the frozen manifest had itself
promised the exact route. Two gates evaluated one predicate, one mutation arm asserted that
adding one to a value changes it, and the environment record overstated its own scope.

This attempt produced the finding that the instruction was not amended afterwards: the
coverage enforcement added for it verifies that a mutation key exists and reads `RED`, and a
hard-coded literal satisfies exactly that, so the mechanism built to guarantee mutations exist
cannot distinguish an executed mutation from an asserted one.

### Attempt 5, held: the mathematical reproduction phase completed, qualification did not

| Field | Value |
| --- | --- |
| Instruction received | `TASK.md` `e3d9b908…`, unchanged from attempt 4, deliberately |
| Outputs landed, unread | `ba38f2e5` |
| Disposition | PRE-REVEAL INCOMPLETE, QUALIFICATION HOLD. Frozen as the Phase A candidate under Addendum 1 |

Running on an unamended instruction, this attempt independently repaired essentially every
substantive defect its predecessors had shown. The representation constructions were validated
before the manifest froze, including an explicit record that the self-conjugacy fact whose
absence retired attempt 2 was checked. The exact integral saturation route this protocol
demands was supplied, with unimodular minors obtained for all three nonzero boundary images.
No verdict anywhere across nine scripts is a hard-coded literal. Only two coverage rows claim
not-checkable, both because the torsion output is itself the target and cannot be checked
target-blind. Ordering, the firewall, byte-identical reproduction and the intent scan all
hold.

It falls short on one requirement. Fifteen of nineteen declared gate mutations were never
executed: the `mutation` member of each gate result is a description rather than a record of
an execution observed to fail. Of the seven validation artifacts, only the fixture executes
mutations, four of them, and the production script contains no gate identifiers, no coverage
logic and no mutation execution. There is no coverage checker at any level and the run exits
zero.

That is materially attempt 4's defect expressed as prose rather than as a literal, and § 9
refuses both in the same sentence.

## What the sequence shows

Read in order, the five outcomes move upstream through the same structure:

```text
back-fit the manifest  ->  depart from it instead of stopping  ->  validate incompletely and
enforce nothing  ->  violate explicit gate-contract requirements  ->  complete the
mathematical phase with its internal mathematical gates satisfied while mutation
qualification remained incomplete
```

Every run was a fresh context receiving the same three frozen technical inputs. Across the
later commissions, and still after the earlier methodological defects had been independently
repaired, the recurring unresolved requirement was § 9's executable-mutation and
enforced-coverage contract. Attempts 1 and 2 were retired on manifest-order and
manifest-departure grounds, and this record does not establish a § 9 finding against them;
attempts 3, 4 and 5 each failed § 9 in a different form.

Attempt 5 is the strongest evidence: on the unchanged mature instruction it repaired the
substantive defects while leaving fifteen of nineteen declared mutations unexecuted.

The inference drawn from that is deliberately narrow. The mutation requirement is not
impractical; **for this task, the monolithic single-context commissioning architecture has
not reliably discharged it.** That is an operational diagnosis about this workflow, not a
general law derived from five heterogeneous runs.

The standard is therefore preserved in full and only the division of labor changes. Addendum 1
freezes attempt 5's artifacts as Phase A and commissions a separate pre-reveal Phase B whose
sole task is to build and execute the mutation machinery against those immutable bytes. A
second argument for that split does not depend on the failures at all: an implementer
certifying the adequacy of its own complete test suite is weaker evidence than a separate
context attacking the frozen implementation.

## Rooms, preserved off-tree

Each retired room is kept intact rather than deleted, in the author's out-of-git working
area alongside the commissioning record, one preserved directory per attempt. The rooms are
not in the repository because they are working state rather than evidence. The artifacts they
produced are all here, per attempt, verbatim. At the end of every run the three frozen
technical inputs were reverified byte-unchanged against the § 8 lock manifest, and that
attempt's commissioned `TASK.md` was separately reverified against its own recorded
instruction hash.
