# M8.8 clean-room attempt 1: RETIRED PRE-REVEAL

> **Disposition, 2026-08-19.** The first commissioned clean-room run completed and produced a
> full deliverable set. It is retired before § 8 step 6 for manifest-order nonconformance. No
> answer material was opened at any point, the firewall held, and the run's outputs are
> preserved unaltered. This record exists outside the room and outside the repository until it
> lands with the eventual result.

## The retiring defect, stated exactly

`TASK.md` § 3 and protocol § 8 step 4 require the method-and-gate manifest to be written
first and completely, before any implementation code, and the manifest is the object that
freezes the selected route, the route-native gates, **the conventions consumed**, and the
declared native orientation. `TASK.md` states in terms that a manifest written after the fact
is a different and weaker object.

The manifest was written first. It was then amended, after the implementation existed, to
correct the stated `SU(2)` embedding so that it matched the embedding the code actually used.
An embedding convention is not a typo or a formatting repair; it is precisely the class of
methodological statement the manifest exists to commit in advance.

The causal record therefore reads:

```text
manifest_0  ->  implementation exists  ->  mismatch exposed by the implementation
            ->  manifest_1 corrected to match the implementation
```

rather than the required:

```text
complete manifest  ->  implementation
```

The record under-determines which of two things happened: that the implementer always
intended the code's embedding and mis-transcribed it into the manifest, or that the
convention settled during implementation and was back-filled. Intent is not recoverable from
the artifacts, and the manifest exists precisely so that this question never has to be asked.

**What this is not.** It is not a finding that the computation is wrong. The two embeddings
appear to be relabelings of the same quaternion algebra, under which the torsion values would
be unchanged, and every gate the run reported passed. The defect is to the evidentiary
object, not demonstrably to the number.

## Explicitly NOT the retiring defect

| Item | Why it does not retire the run |
| --- | --- |
| The run-order dependency in `RAW_OUTPUT.json` | disclosed, reproducible, and affects only whether an embedded gate block is present; a leaf-level comparison found 320 of 320 shared leaves equal |
| The hand-transcribed hash error in the landing record | commissioner-side, never reached a commit, and was caught by a check that was functioning as designed |
| Anything about the firewall | it held: all four seeded inputs were byte-unchanged at the end of the run, and no quarantined object was opened |

## What is preserved, and where

| Object | State |
| --- | --- |
| Attempt-1 outputs | commit `e6c9a44b` on branch `m8_8-cleanroom-output`, pushed, UNTOUCHED, no pull request opened, no comparison ever run against it |
| Attempt-1 room | `~/Desktop/OpenWave/M8.8_ATTEMPT1_RETIRED/`, preserved intact rather than deleted |
| Canonical answer packet | SEALED and unopened, in author custody, exactly as before the commission |
| Maintainer-side construction-audit artifact | still quarantined outside the tree, publication remains a § 8 step 9 event |

Nothing about this retirement consumes the case. The answer has never been visible to anyone
in the loop, so a conforming run remains fully available.

## Recommission terms

A second implementer context is commissioned from a FRESHLY REBUILT four-file room extracted
again from the frozen lock commit, never from the populated attempt-1 room.

**The second implementer is not given attempt 1's manifest, code, or outputs.** Handing over
the corrected manifest would convert the first implementer's methodological choice into an
input to the second, which would destroy the independence the recommission exists to buy. It
generates its own complete manifest before touching implementation.

**The commission conditions are not identical, and this record does not pretend they are.**
The three frozen technical inputs are byte-identical, re-extracted from the frozen commits
rather than copied from the retired room. The operational task file differs by exactly one
generic process clause, added because attempt 1 demonstrated that the instruction was
under-specified at one process fork:

| Object | Attempt 1 | Attempt 2 |
| --- | --- | --- |
| `PROTOCOL.md` | `23f313d0cd47e4ff644ae7c9730cbc9eb380cb6f7c64e9daddb7b98cd6885d87` | identical |
| `m8_5a_packet.json` | `e3b0c945bbbb15b4549fa641234c9461062c2337b3d1e372af621b614d4883a9` | identical |
| `m8_8_construction_packet.json` | `df00c0222f98c481eb56b882cd867a6c3a4f8604b8633e81dec0cce1f8460a06` | identical |
| `TASK.md` | `c61320750c082dc7b2927bdb7d8c7d3fe4bf7d230b31462a5509ee1c7a871363` | `088200506c0a8eb3b81b4de501e4147cee23853c77eb83c168f5da93a5bedc43` |

The added clause states that once implementation begins the manifest is immutable, and that a
mismatch exposed by implementation is resolved either by conforming the implementation to the
still-intended manifest or by stopping and reporting, never by amending the manifest to match
the implementation. It is deliberately phrased to cover the third case as well: where the
manifest itself proves substantively wrong, the correct response is to stop rather than to
force the implementation toward a convention now believed incorrect.

**The revision carries no technical content from attempt 1.** It names no convention, no
embedding, no orientation, no route, no gate, no value, and no outcome, and it does not
disclose that an earlier commission existed. It is a process repair derived from the fact of
the retirement rather than from anything the retired run computed, and it is legitimate
precisely because the answer is still sealed.

That makes the recommission a genuine test rather than a retry. Three outcomes are all
informative: the second implementer independently selects the same convention and reproduces
the same values; it selects a different permitted route or convention, which is equally
legitimate; or it meets the same ambiguity and stops before implementing, which would say
something real about the sufficiency of the frozen inputs rather than about either
implementer.

## Why this record makes the eventual result stronger

A reproduction whose commission history shows a run retired for a process defect, with the
answer still sealed and the defect named precisely, is worth more than one with no visible
failures. The alternative was to proceed to comparison carrying an amended manifest and to
discover the objection after the answer was visible, at which point no clean run on this case
would have been possible again.
