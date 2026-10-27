# T10: A validation receipt per reproduction run

> Platform task **T10**. Status: 🚧 **STAGED 2026-10-27**, no go yet · Roadmap:
> [`../platform_roadmap.md`](../platform_roadmap.md) · Deliverable lands beside each
> reproduction record under `research/data/<id>_receipt.json` with its rule in
> [`METHOD_NOTE.md`](../METHOD_NOTE.md).

## PLANNING

**Scope.** One small JSON per reproduction or verification run, in a fixed shape: the base commit of this repo, the commit or DOI of the package reproduced, the interpreter and the pinned library versions, the exact command per job, the exit code and wall time, the comparator's output per file, and the verdict words. The M5.33 to M5.35 records already hold every item in prose and in `*_reproduction.json`; the receipt fixes the shape so a reader, or a later run, can diff two receipts without reading two task docs. The pattern is the `validation-receipt.yaml` of the repository that hosts the M5 discussion (base, head, declared impact surface, command, result), reduced to what a reproduction needs.

**Definition of done.** A 60-line writer `dev_docs/utils/write_receipt.py` that takes a logs folder and a comparator output and writes the receipt; the shape documented in `METHOD_NOTE.md` (one table); the M5.35 run converted as the first receipt from its existing logs; the next reproduction run writes one at FINISH by rule.

**Gating.** None beyond the user's go.

**Blindspot pass.** (i) A receipt is not evidence of correctness, only of what ran: the verdict words stay in the record with their gates. (ii) Machine-specific numbers (wall times) are in the receipt because they are a run's facts, with the machine named by class only (the receipt never carries a hostname or a path). (iii) The receipt is tracked (small JSON), never the logs.
