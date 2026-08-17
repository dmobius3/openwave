# M8.5-B Phase A: freeze record

    Declared     2026-08-13, describing the corrected build produced
                 in response to the second review round. Acceptance
                 has not occurred; the evidence of it will be the
                 merge and its commit ancestry, not this record. An
                 earlier record dated 2026-08-12 described the pre-review
                 build and is superseded by this one
    Scope        the scalar and one-form adjudication machinery for rungs 3a
                 and 3b, as published in this directory
    Manifest     `qualification/MANIFEST.json`
                 SHA-256 9e6dcd7b5d177040cf25f1eda0bb4614f4002cd85f25660ca9ec75dc9734d934
    Reproduce    python3 run_qualification.py

## What is frozen

The implementation in this tree, identified by `MANIFEST.json`, which carries
the SHA-256 of every shipped file except three. Two are excluded because they
would close a hash cycle: the manifest itself, which cannot contain its own
hash, and this record, which carries the manifest's. The third, `.gitignore`,
is excluded for a different reason: it is repository hygiene rather than
qualification evidence, no battery reads it, and repository policy marks it
`export-ignore`, so pinning it made a correct exported tree fail this gate.
The chain runs freeze -> manifest -> every other file, and this record's own
anchor is the commit that carries it. Changing any shipped file changes the
manifest hash and therefore breaks this pin, which is the point.

## What the freeze claims

The adjudication machinery that Addendum 12.1 specifies exists here, is
qualified, and was qualified **before any sealed case it will consume was
created**. That ordering is the whole claim: the comparator cannot have been
shaped by an answer it had not seen.

The qualification is reproducible, from a clean clone and from an exported
source archive alike, by the single command above, which establishes, in
order: no first-party module imported by the controlling process resolves
outside this tree; the route (a) group-closure battery added by Addendum
12.3, whose mutation arm proves the central-equivalence repair is
load-bearing; the confinement those batteries run under, observed in a
child and a grandchild launched through the same helper they use and probed
under a deliberately poisoned parent environment; the environment; the
manifest; the Packet I and Packet II gate batteries; the structural battery
target-scored at 8 of 8 with non-vacuity proved by suppressing each predicate
in turn; the Q4 integration rehearsal including the route-(b) deletion limb;
and the integrated battery.

## What the freeze does NOT claim

**No qualification run in this tree is adjudication evidence.** `M85B-ADJ-04`
invoked route (a) before this repair and stopped structurally before producing
an output; it is retired under Addendum 12.3 and its Packet II remains sealed.
Every case exercised by this qualification command is synthetic or a frozen
tuning case, no external reference value has been compared to any route
output, and no certification language from § 1 of the preregistration is
engaged. The M8.5-A gate condition and every § 0 claim ceiling stand
unchanged.
