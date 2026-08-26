# S1b addendum 2: round-1 evidence is not available, and the transfer branch is inoperative

> **APPEND-ONLY.** The frozen regions of `S1B_DECISION_RULE.md`
> (`c44c603a8058ed8529e5bb0f42ec168b443ff22b050d3748f1f1de3537c7d297`) and `S1B_ADDENDUM_1.md`
> (`6da36a1c672772e8c731d4d66f20a77e81779d19d5cb883e617c028d43f46746`) are UNCHANGED and both still
> verify. This addendum adds one provenance fact and its consequence. It has its own freeze record.
>
> **PRECEDENCE.** Where addendum 1 licenses reuse of historical round-1 evidence, this addendum
> supersedes it for the present qualification. Addendum 1's mathematics is untouched; only the
> availability of the evidence it would have transferred has changed.

## A2.1 The fact

Round 1 delivered code but **no qualification note and no result records**. Its results exist only
as a terminal transcript, which is not a hash-pinned artifact and receives no evidentiary weight.

The pinned script `s1b_qualification.py` proves what code existed in round 1 but is NOT a replayable
result artifact, for two independent reasons, both verified:

- it resolves its room as `ROOM = os.path.dirname(os.path.abspath(__file__))` at line 18, so it
  reads `contract/`, `p0/` and the manifest relative to its own location;
- at line 148 it asserts `ROOM_MANIFEST.json` hashes to `46ab992545c6a28cd6d116c238aa00ed5525f6c11f39688ea832172d6869d20f`,
  the ROUND-1 input manifest, and halts otherwise.

Any round-2 room differs from the round-1 room by construction, since it contains addendum 1. That
assertion can therefore never pass, and the script cannot execute in any room in which it would be
useful. Reconstructing a bit-exact round-1 replica is also not possible: the round-1
`room_import_gate.py` was edited in place afterwards and its bytes are not recoverable, so the
manifest value it authenticates cannot be reproduced.

**This is a defect in how round 2 was staged, not in the script and not in addendum 1.** The
commissioned unit encountered it and correctly moved toward a STOP; its attempted workarounds,
rewriting `__file__` and bypassing the manifest assertion, would have destroyed the provenance
discipline the assertion exists to enforce, and are prohibited below.

## A2.2 The consequence

**No round-1 result receives evidentiary credit.** Addendum 1's historical-evidence transfer branch
is INOPERATIVE for this qualification. There is nothing to transfer, so the branch has no subject.

- The shipped `invariant_dim_and_basis` output feeds the qualifying path DIRECTLY. It is run once
  per target, its returned dimension, basis and diagnostics are serialized immediately as generated
  qualification output, that output is hashed, and every downstream gate consumes that exact
  serialized artifact.
- The bitwise shipped-versus-economy comparison of `s` and `Vh` **remains required and remains
  recorded**, as an equivalence DIAGNOSTIC. It no longer licenses anything, because a license to
  reuse nonexistent evidence is empty. A bitwise match is a fact about two constructions, not a
  transfer.
- **Selective credit is prohibited.** `G-REAL`, `G-ALIGN`, `G-WIRE`, `G-BASIS`, `G-DISCRIM` and the
  adjudicator cases are in exactly the same evidentiary position as the SVD-route results: recalled
  from a transcript, not preserved. Crediting the ones that are remembered favourably while
  rerunning the ones that are not would reproduce the failure this addendum records. The
  qualification is therefore run in FULL.

## A2.3 The round-1 artifact is archival

`prior/s1b_qualification.py` is archival evidence of what round 1 ran. Its pinned hash is verified
and it is READ-ONLY in every sense: it may not be executed, relocated, patched, run with an altered
`__file__`, have its manifest assertion bypassed, or have any output of it used for qualification
credit. Today's discovered impossibility becomes a guard rather than another trap.

## A2.4 Evidence must be produced incrementally, not at the end

Two independent failures have now shown that a terminal transcript is not a durable scientific
artifact: round 1 finished and wrote no note, and round 2 died at minute 16 with everything it had
established lost. Both would have survived an artifact written as it went.

**Frozen for this and every subsequent S1b qualification.** After Q0 provenance succeeds, and
BEFORE any further work, the unit creates an output ledger and writes Q0 into it immediately, then
appends after every completed gate, control and mutation:

    qualification_round2/QUALIFICATION_NOTE.md      appended as facts are established
    qualification_round2/results/<gate>.json        one machine-readable record per gate
    qualification_round2/OUTPUT_MANIFEST.json       written last, pinning note, code and records

The output ledger is a separate provenance layer from the room. The INPUT manifest pins the frozen
contract and the supplied code and is never modified by the unit. The OUTPUT manifest pins the
qualification's own code, note and measured results. A run that dies partway then leaves everything
it had established on disk and hashable.

**Interruption does not license resumption of scientific credit.** A later unit may INSPECT a
partial ledger, but unless a checkpoint-continuation protocol is separately commissioned it reruns
the qualification independently rather than inheriting another run's partial results. That is the
same rule this addendum applies to round 1, applied to itself.
<!-- ADDENDUM2-BOUNDARY -->

**Freeze record, addendum 2.** SHA-256 covers every byte ABOVE the boundary comment: `14011c338ce331d9bb3424f0ed5aaff31a1c89a9ceafbc2bae41a91ebb41a222`

```bash
sed '/^<!-- ADDENDUM2-BOUNDARY -->$/,$d' S1B_ADDENDUM_2.md | shasum -a 256
```

The parent rule and addendum 1 are untouched and verify independently.
