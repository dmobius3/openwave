# Historical commission records, NOT governing

**Nothing in this directory is governing, adjudicative, or part of the frozen evidence chain.**
These are the briefs the commissioned units were given and the room provenance that surrounded
them. No branch rule reads any of it, no freeze hash covers it, and the `S1b-SPECTRAL` adjudication
does not depend on a single byte here. The governing documents are the six frozen files in
[`../../`](../../), and the adjudicative evidence is `q3a_output_manifest.json` and what it pins.

They are kept because the S1b closeout turns on a commission-boundary failure, and a reader
evaluating that failure should be able to see what the units were actually told rather than take
the closeout's account of it.

| file | what it is |
| --- | --- |
| `qualification_TASK.md` | the final full-requalification brief, the one attempt `q3a` executed |
| `round1_TASK.md` | the first brief, whose run produced code but no note |
| `adjudication_TASK.md` | the read-only brief for the unit that issued the branch |
| `qualification_room_import_gate.py` | the launch precondition for the qualification room |
| `qualification_ROOM_MANIFEST.json` | that room's INPUT manifest, 82 entries |
| `adjudication_ROOM_MANIFEST.json` | the adjudication room's manifest, 29 entries |

## Why the import gate is here

Several provenance defects the addenda describe were found BY this gate rather than by reading:
a missing `route_gates` dependency on the room's first build, a further `packet_schema` dependency
under it, and a false-green in which modules calling `sys.exit(0)` at import terminated the checker
with status 0 and no output, so the gate reported nothing and looked fine. It is archived as the
instrument that caught them.

## Provenance of these copies

Both room manifests were verified against their rooms immediately before archiving: 82 of 82 and
29 of 29 input files unchanged since their runs. The Desktop rooms were then deleted, so these
manifests describe directories that no longer exist. The adjudication room's `evidence/` copies are
deliberately NOT archived, being duplicates of artifacts already pinned in
`q3a_output_manifest.json`, and a second copy of pinned bytes would only raise the question of which
is canonical.

**These are not roots of trust.** The two files here named `ROOM_MANIFEST.json` authenticate the
HISTORICAL ROOM STATES described on this page and nothing else. They are not part of the chain of
custody for the adjudicated `S1b-SPECTRAL` result, and no verification of that result should route
through them. The canonical result provenance is the six frozen S1b contracts, `q3a`'s
`OUTPUT_MANIFEST.json`, and the pinned evidence artifacts, all in the closeout package at
[`../../`](../../) and listed in the closeout's own chain-of-custody table. The naming similarity is
the only thing these manifests share with that one.

## The round-1 brief is superseded, not authoritative

`round1_TASK.md` describes a commission whose results receive no evidentiary credit, per
[`../../s1b_addendum_2.md`](../../s1b_addendum_2.md) § A2.1. It is here as the record of what was
asked, not as a live instruction. The same applies to the qualification brief's prohibition on the
live target, which attempt `q3a` crossed: addendum 4 § A4.1 records that the brief and the gates it
required were mutually inconsistent, and that the fault is the contract author's.
