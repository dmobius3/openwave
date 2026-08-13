# M8.5-B Phase A: freeze record

    Declared     2026-08-13, describing the build accepted after the
                 second review round. An earlier freeze record dated
                 2026-08-12 described the pre-review build and is
                 superseded by this one
    Scope        the scalar and one-form adjudication machinery for rungs 3a
                 and 3b, as published in this directory
    Manifest     `qualification/MANIFEST.json`
                 SHA-256 f79c4c9d8a0fc8bcab5244d6422baa1f3955fe925f49a676e740e0a4ba9605ae
    Reproduce    python3 run_qualification.py

## What is frozen

The implementation in this tree, identified by `MANIFEST.json`, which carries
the SHA-256 of every shipped file except two: itself, which cannot contain its
own hash, and this record, which contains the manifest's hash and would
otherwise close a cycle. The chain runs freeze -> manifest -> every other file,
and this record's own anchor is the commit that carries it. Changing any
shipped file changes the manifest
hash and therefore breaks this pin, which is the point.

## What the freeze claims

The adjudication machinery that Addendum 12.1 specifies exists here, is
qualified, and was qualified **before any sealed case it will consume was
created**. That ordering is the whole claim: the comparator cannot have been
shaped by an answer it had not seen.

The qualification is reproducible from a clean checkout by the single command
above, which establishes, in order: no first-party module imported by the
controlling process resolves outside this tree, and every battery runs under
subprocess path confinement; the environment; the manifest; the Packet I and
Packet II gate batteries; the structural battery target-scored at 8 of 8 with
non-vacuity proved by suppressing each predicate in turn; the Q4 integration
rehearsal including the
route-(b) deletion limb; and the integrated battery.

## What the freeze does NOT claim

**No rung has run on a real adjudication case.** Every case exercised here is
synthetic or a frozen tuning case, so nothing in this directory is adjudication
evidence, no external reference value has been compared to any route output, and
no certification language from § 1 of the preregistration is engaged. The
M8.5-A gate condition and every § 0 claim ceiling stand unchanged.
