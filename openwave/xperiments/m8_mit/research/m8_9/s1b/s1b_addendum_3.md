# S1b addendum 3: the output ledger is per attempt, not per round

> **APPEND-ONLY.** The frozen regions of `S1B_DECISION_RULE.md`
> (`c44c603a8058ed8529e5bb0f42ec168b443ff22b050d3748f1f1de3537c7d297`), `S1B_ADDENDUM_1.md`
> (`6da36a1c672772e8c731d4d66f20a77e81779d19d5cb883e617c028d43f46746`) and `S1B_ADDENDUM_2.md`
> (`14011c338ce331d9bb3424f0ed5aaff31a1c89a9ceafbc2bae41a91ebb41a222`) are UNCHANGED and all three
> still verify.
>
> **PRECEDENCE.** This addendum supersedes addendum 2 § A2.4's ledger PATHS and nothing else.
> A2.4's requirements, that the ledger is created before any other work, appended to after every
> completed gate rather than batched, kept separate from the input manifest, and never used to
> inherit another run's scientific credit, all stand unchanged.
>
> This is procedural provenance, not scientific design. No gate, threshold, control or measurement
> is touched.

## A3.1 The defect

A2.4 froze a single fixed directory, `qualification_round2/`, for "this and every subsequent S1b
qualification", while simultaneously providing that a run may be interrupted and that a LATER unit
reruns independently rather than inheriting the partial ledger.

Those two provisions are incompatible. A later unit rerunning into the same fixed directory has only
bad options: append to a dead run's note, overwrite it, delete it, or invent an unauthorized
pathname. Each defeats the provenance the ledger exists to create, and the last is precisely the
kind of unilateral choice the contract forbids.

This has already cost this program three artifacts: a round-1 note that was never written, a round-2
run whose evidence died with its terminal, and a round-2 commission overwritten in place while
staging the fix for the first two.

## A3.2 The ledger is keyed by attempt

    qualification_runs/<attempt_id>/QUALIFICATION_NOTE.md
    qualification_runs/<attempt_id>/results/<gate>.json
    qualification_runs/<attempt_id>/OUTPUT_MANIFEST.json

`<attempt_id>` is **supplied in the handoff prompt before execution and is never chosen, derived or
altered by the unit.** A later attempt after an interruption receives a new id, so evidence from
distinct attempts can never mix.

**Collision guard, and it is a STOP.** Before writing anything, the unit checks whether
`qualification_runs/<attempt_id>/` already exists. If it does, the unit STOPS and reports, without
reading, writing, appending to or deleting anything inside it. A pre-existing directory means either
the id was reissued in error or an earlier attempt used it; both require a decision the unit is not
authorized to make.

**A unit writes only inside its own attempt directory.** Earlier attempts' ledgers are inspectable
by a human or a commissioner, never by a running unit, and never modified.

## A3.3 The input layer stays immutable, and is re-verified at the end

A2.4's two-layer separation stands: the INPUT manifest pins the frozen contracts and supplied code
and is never modified by the unit; the OUTPUT manifest pins that attempt's own note, code and result
records.

**Frozen addition.** Immediately before writing `OUTPUT_MANIFEST.json`, the unit re-hashes every
file listed in `ROOM_MANIFEST.json` and re-verifies all frozen-document hashes. **Any change to a
manifest-pinned input since Q0 is a STOP**, and the partial ledger is left in place as evidence.

This is a re-verification of INPUTS only, not a rerun of the import gate: by that point the room
legitimately contains generated qualification output that the input manifest does not list, so the
gate's own accounting would no longer be the right instrument. Given how many failures in this
program have involved a room changing underneath its own report, an end-of-run input check is worth
its cost.
<!-- ADDENDUM3-BOUNDARY -->

**Freeze record, addendum 3.** SHA-256 covers every byte ABOVE the boundary comment: `e3304fe993f7a2ea523cfcbc4110c0f37174caea95383b17e465c8fa85b6e28c`

```bash
sed '/^<!-- ADDENDUM3-BOUNDARY -->$/,$d' S1B_ADDENDUM_3.md | shasum -a 256
```

The parent rule and addenda 1 and 2 are untouched and verify independently.
