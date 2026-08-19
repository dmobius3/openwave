# M8.8 clean-room attempt 2: RETIRED PRE-REVEAL

> **Disposition, 2026-08-19.** The recommissioned run completed and produced a full deliverable
> set. It is retired before § 8 step 6 for manifest-departure nonconformance. No answer
> material was opened, the firewall held, and the outputs are preserved unaltered at
> `d3c7859a`.

## The retiring defect, stated exactly

The task file offered two paths when implementation exposes a mismatch with the frozen
manifest: conform the implementation to the still-intended manifest, or STOP and report.

The implementer established that the manifest's construction for one object was
mathematically wrong and therefore could not remain the intended method. Conforming was
genuinely unavailable. At that point the only authorized branch was STOP. Instead it reported
the mismatch and continued, substituting a different construction and carrying the run to
completion.

**That is the whole defect.** It is narrower than attempt 1's and it is procedural rather than
mathematical.

## What is explicitly NOT the defect

| Item | Status |
| --- | --- |
| The mathematics | not suspected wrong. The implementer's diagnosis was verified independently on the commissioner side, with machinery sharing nothing with the implementation, and it is correct: the manifest's named object coincides with one already present rather than supplying a distinct one |
| The manifest | NOT amended. It was written first and its bytes were never modified again. On the failure mode that retired attempt 1, this run did exactly the right thing |
| The firewall | held. All four seeded inputs were byte-unchanged at the end of the run and no quarantined object was opened |
| The disclosure | exemplary. `MANIFEST_MISMATCH.md` quotes the frozen statement, explains why it fails, names the replacement used, and states plainly that the manifest was not amended |
| Reproducibility | the implementation reruns from its committed bytes at exit 0 and regenerates its raw output byte-identical |

## Why the "forced repair" argument does not save it

The commissioner-side analysis was that the repair introduced no discretion, since the group
has exactly two objects of the relevant kind, so once the manifest's named one is shown to
coincide with the other, the replacement is uniquely determined. That analysis is correct as
mathematics and was verified.

It is still the wrong basis for acceptance. Accepting here would install a rule that was never
pre-registered:

> an implementer may depart from the frozen manifest after implementation starts, provided the
> commissioner later determines the correction was mathematically forced.

That rule requires adjudicating, after seeing each implementation, whether a departure was
forced enough. Post-hoc judgment of exactly that kind is what a precommitted manifest exists
to abolish, so an exception granted on those grounds dissolves the thing it is an exception
to.

The point is sharper still because the instruction was not ambiguous this time. The clause was
added specifically to cover the case where the manifest itself proves substantively wrong, and
the attempt-1 disposition says so in terms. Waiving it on the first occasion it fired would
leave it with no force in any later run.

## Two distinct failure modes, now both on record

| Attempt | Failure mode |
| --- | --- |
| 1 | the manifest was changed to fit the implementation |
| 2 | the manifest stayed frozen, and the implementation departed from it instead of stopping |

These are different, and neither is a repetition of the other. Both share one upstream cause:
the manifest was required to be written before implementation, but never required to be
CHECKED before it froze.

## The repair for attempt 3: validate before freezing, not loosen after

Attempt 3's task file adds one clause immediately before the immutability sentence, requiring
that every manifest-defined construction or convention checkable from the permitted inputs
undergo a target-blind internal-consistency check, recorded in the manifest, before the
manifest is complete, and permitting defects found in that phase to be corrected and recorded
before implementation begins.

The immutability sentence is retained unchanged. Together they give a three-part causal
structure with corrections allowed only on the left of the freeze:

```text
draft manifest  ->  manifest self-validation  ->  FINAL MANIFEST  ->  implementation
```

This is a process improvement rather than a retry. Commissioning fresh contexts until one
happens to author a flawless manifest would be testing authoring luck; requiring the method to
be validated before it becomes immutable is a reusable gate that addresses both recorded
failure modes at once.

## What attempt 3 does not receive

No mention of any earlier commission, and nothing from attempt 2's mathematics: no object
names, no construction, no representation, no convention, and no result. The added clause is
written generically and names none of it. Attempt 3 receives only the stronger requirement
that its own manifest be checked before code exists.

## Preserved

| Object | State |
| --- | --- |
| Attempt-2 outputs | `d3c7859a`, pushed, no pull request, no comparison ever run |
| Attempt-2 room | `~/Desktop/OpenWave/M8.8_ATTEMPT2_RETIRED/`, intact |
| `MANIFEST_MISMATCH.md` | kept deliberately; it belongs in the eventual process history as the artifact that proved the immutability clause works |
| Canonical answer packet | SEALED, never opened. The case remains unconsumed |
