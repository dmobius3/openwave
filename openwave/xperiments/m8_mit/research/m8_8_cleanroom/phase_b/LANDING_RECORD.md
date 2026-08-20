# M8.8 Phase B: the qualification output, landed unread against the answer packet

> **Claims no result.** The answer packet has not been opened, no comparison has run, and no
> verdict exists. This commit discharges the ordering requirement that Phase B's deliverables
> land before § 8 step 6.

## What landed

| File | SHA-256 | Role |
| --- | --- | --- |
| `run_phase_b.py` | `cae7993897d483ab3ced79df1ade4f3e371569ed0d822955a674b742945ff772` | the qualifier: executes all nineteen declared mutations against the frozen machinery |
| `MUTATION_RESULTS.json` | `a22e6da390ab31c477a97a0de29a3734ca2f7fda8c61d5dd08854860c4e59a7a` | the execution record, one per gate, with object mutated, predicate, baseline, mutated result and observed red |
| `QUALIFICATION_RECORD.md` | `99fae82dc8a4369780d19ee842ebe4e7ee090bab5b88f592c06dea4834f3001a` | the qualifier's own summary |

## Coverage, enumerated by the commissioner rather than reported

Parsing the frozen manifest's § 4 registry independently and comparing against the execution
record:

| Check | Result |
| --- | --- |
| Registry gates, my own parse | 19 |
| Executed mutation records | 19 |
| Exact set equality | yes, no missing, no duplicate, no unregistered gate |
| Records carrying all required fields | 19 of 19: object mutated, gate predicate, baseline result, mutated result, observed red |
| Records with `red_outcome` true | 19 of 19 |
| Records showing a PASS baseline moving to a FAIL under mutation | 19 of 19 |
| Stub or empty mutation objects | none; each names a concrete mutation |
| Phase A integrity | all 13 hashes verified before and after by the qualifier, and independently by me at this tip; no Phase A byte moved |
| Reproduction | the qualifier reruns from the committed bytes at exit 0 |
| A9 | standard library only, no network, no subprocess, no eval |

## The qualifier's own gate, mutation-tested by the commissioner

The qualifier is the object that certifies coverage, so its accept side was attacked
directly. Three mutations, each applied to a scratch copy:

| Mutation | Result |
| --- | --- |
| Append a byte to a frozen Phase A artifact | exit 1, `HASH MISMATCH: compute_torsion.py` |
| Alter one gate row in the frozen manifest's registry | exit 1, `HASH MISMATCH: METHOD_AND_GATE_MANIFEST.md` |
| Drop one executed record before the coverage comparison | exit 1, `MISSING: {'G-M03'}`, exact set equality false |
| Force one mutation's `red_outcome` to false | exit 1, `G-M01: mutation did not redden` |

The clean run returns to exit 0 after each. The coverage gate discriminates in every
direction it is required to.

## One finding, stated rather than deferred

**The expected gate set is a hard-coded Python list, not parsed from the manifest.**
`run_phase_b.py` carries `MANIFEST_GATES` as a literal of nineteen identifiers. Addendum 1
requires that the qualifier "parses the nineteen gate identifiers and their declared
mutations from `METHOD_AND_GATE_MANIFEST.md` § 4", and `QUALIFICATION_RECORD.md` asserts that
it did so. The mechanism is transcription.

What is verified about it, so the finding is not overstated: the transcription is exactly
correct, confirmed against an independent parse at nineteen of nineteen; and it is bound to
the manifest by hash, since altering a single registry row in the manifest aborts the run
with a hash mismatch before the list is used. The failure mode the rule exists to prevent, a
duplicated list omitting the same gate the runtime omits, is therefore closed in this
instance.

What remains true regardless: the addendum specifies parsing, this transcribes, and the
qualification record claims a mechanism the code does not use. An earlier attempt was retired
in part for hard-coding an expected gate set, and a "the risk is closed here" argument of
exactly this shape was raised and rejected once before in this task, on the grounds that it
requires the commissioner to adjudicate after the fact whether a deviation was harmless
enough. That adjudication is the adjudicator's, and this record does not make it.
