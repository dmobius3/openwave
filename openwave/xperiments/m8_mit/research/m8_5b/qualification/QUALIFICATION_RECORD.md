# M8.5-B Phase A: Qualification Record and Method Note

**This document serves as the § 9 method note for the Phase A qualification
build.** It states what this machinery is, what was tested and what was not, how
to reproduce the result in one command, and what each record here means.

## What this is

The scalar and one-form adjudication machinery for rungs 3a and 3b: the packet
sealing gate, the two comparison surfaces that Addendum 12.1 specifies, the
Step-3 record surface and its full-lattice validator, both route producers, and
the standalone rung-3b theorem evaluator.

## Qualification is not adjudication

**No rung has run on a real adjudication case.** Every case exercised in this
directory is synthetic or a frozen § 6.1 tuning case. Nothing here compares an
external published reference value against a route output, so nothing here is
adjudication evidence and no § 1 certification language is engaged. What is
established is narrower and is the point: the machinery that will consume a
future sealed answer packet exists, is qualified, and was qualified before that
packet was created.

## Reproducing it

    python3 run_qualification.py

One command, from a clean checkout. It establishes, in order:

| step | what it shows |
| --- | --- |
| provenance | every first-party module resolves inside this tree, so a green run cannot be borrowing code from elsewhere |
| environment | the interpreter and library versions, compared against the recorded ones |
| manifest | every shipped file matches its recorded SHA-256 |
| schema | the Packet I and Packet II gate batteries pass, and the suite demonstrates it can fail |
| structural | the target-scored battery covers 8 of 8 predicates, and suppressing any one of them reddens its own item |
| rehearsal | the Q4 integration rehearsal is GREEN for 3a and 3b on both routes, and the deletion limb completes the adjudication from committed artifacts alone with the route-(b) recompute physically removed |
| integrated | the battery returns 31 of 31, every failure asserted in its required layer |
| records | the fresh run reproduces the shipped records |

Expected verdict:

    PHASE A QUALIFICATION: PASS - structural 8/8; integrated 31/31;
    Q4 GREEN; deletion GREEN

## What the records mean

    qualification/QUALIFY_RECORD.json    the integrated battery, per item, with
                                         the layer each failure was required to
                                         occur in. 31 items, 0 failed
    qualification/REHEARSAL_RECORD.json  the Q4 rehearsal: packet hashes, the
                                         committed artifacts, and the deletion
                                         limb's absent-file set
    qualification/ENVIRONMENT.json       the environment qualification was
                                         recorded under
    qualification/MANIFEST.json          SHA-256 of every shipped file except
                                         itself, which cannot contain its own
                                         hash
    qualification/PHASE_A_FREEZE.md      pins the manifest, and therefore the
                                         bytes; states the ordering claim

A fresh run regenerates its outputs under `rehearsal/`, leaving the shipped
records untouched, so the two can be compared rather than one overwriting the
other.

## On `pilot/`

Those modules are numerical support code retained because the qualified
producers import them. They are dependencies, not evidence, and the directory
name is historical.

## The ordering claim

This machinery is frozen, and its bytes are pinned by a manifest that is itself
pinned by the freeze record, **before any answer packet it will consume
exists**. A comparator cannot be shaped by an answer it has not seen, and that
is the property the freeze exists to make checkable rather than asserted.
