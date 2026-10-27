# T10: A validation receipt per reproduction run

> Platform task **T10**. Status: ✅ **CLOSED 2026-10-27** (go 11:29 EDT, batched with T9) · Roadmap:
> [`../platform_roadmap.md`](../platform_roadmap.md) · Deliverable lands beside each
> reproduction record under `research/data/<id>_receipt.json` with its rule in
> [`METHOD_NOTE.md`](../METHOD_NOTE.md).

## PLANNING

**Scope.** One small JSON per reproduction or verification run, in a fixed shape: the base commit of this repo, the commit or DOI of the package reproduced, the interpreter and the pinned library versions, the exact command per job, the exit code and wall time, the comparator's output per file, and the verdict words. The M5.33 to M5.35 records already hold every item in prose and in `*_reproduction.json`; the receipt fixes the shape so a reader, or a later run, can diff two receipts without reading two task docs. The pattern is the `validation-receipt.yaml` of the repository that hosts the M5 discussion (base, head, declared impact surface, command, result), reduced to what a reproduction needs.

**Definition of done.** A 60-line writer `dev_docs/utils/write_receipt.py` that takes a logs folder and a comparator output and writes the receipt; the shape documented in `METHOD_NOTE.md` (one table); the M5.35 run converted as the first receipt from its existing logs; the next reproduction run writes one at FINISH by rule.

**Gating.** None beyond the user's go.

**Blindspot pass.** (i) A receipt is not evidence of correctness, only of what ran: the verdict words stay in the record with their gates. (ii) Machine-specific numbers (wall times) are in the receipt because they are a run's facts, with the machine named by class only (the receipt never carries a hostname or a path). (iii) The receipt is tracked (small JSON), never the logs.

## FINDINGS (2026-10-27)

| Deliverable | Where | State |
| --- | --- | --- |
| The writer | [`write_receipt.py`](../utils/write_receipt.py): logs folder plus comparator text in, one JSON out; the interpreter recorded by name, never by path; exit 1 when any job's exit is non-zero (105 lines after `black`, against the plan's 60: the argument parsing and the docstring are the difference, logged as a deviation) | ✅ |
| The shape | [`METHOD_NOTE.md` § The validation receipt](../METHOD_NOTE.md#the-validation-receipt): nine fields, the machine by class only | ✅ |
| The first receipt | [`m5_35_receipt.json`](../../openwave/xperiments/m5_liquid_crystal/research/data/m5_35_receipt.json) from the M5.35 logs: 60 jobs, 0 non-zero, 24 files compared, base `fea758b6`, the package pinned by DOI | ✅ |
| The mutation check | a copy of two logs with one `exit=0` flipped to `exit=1` turns the writer's exit to 1 and lists the job under `exit_nonzero` (the AI_HYGIENE sentence of 2026-10-27, applied to the writer itself) | ✅ |
| The rule | the pre-send checklist of `METHOD_NOTE.md`: a reproduction or verification run writes its receipt at FINISH | ✅ |

**Scope**: a record format and one converted run; a receipt says what ran and never that it ran correctly, and the M5.35 verdict words it carries are the task record's, not the receipt's.
