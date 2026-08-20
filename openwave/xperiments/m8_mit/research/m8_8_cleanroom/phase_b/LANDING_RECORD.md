# M8.8 Phase B: the repaired qualification, landed unread against the answer packet

> **Claims no result.** The answer packet has not been opened, no comparison has run, and no
> verdict exists. Lands before § 8 step 6. The first Phase B landing, `f441e0ec`, is preserved
> in history and is not rewritten.

## What landed

| File | SHA-256 | Role |
| --- | --- | --- |
| `run_phase_b.py` | `3a621d0ebaa5a6b7819e4ad3f57825f1bda064684d36f49cdd4eb99d913bfba7` | the repaired qualifier: parses the § 4 registry at runtime, then executes all nineteen declared mutations |
| `MUTATION_RESULTS.json` | `1540d2a53aaa63f5f6850544fc366df781773d3420425b2f30a0a5eb6836f460` | the execution record, schema 2, carrying the parsed, implemented and executed identifier sets |
| `QUALIFICATION_RECORD.md` | `2fef2b8aacb35e1f9b0765c62f06ab37a999f010f2e2d6f9e264cd0d4d58f55f` | the qualifier's summary, its parsing statement now literally true |

## The mechanism nonconformance is cured

The prior landing hard-coded its expected gate set while its record claimed runtime parsing.
The qualifier now parses the frozen manifest's § 4 registry at runtime, extracting each gate
identifier and its declared mutation, and rejects malformed rows, malformed identifiers,
empty mutation declarations and duplicate identifiers. No hard-coded gate list remains
anywhere in the file. The qualification record's parsing statement is now literally true.

## Commissioner enumeration, run before this record was written

| Check | Result |
| --- | --- |
| Phase A artifacts | 13 of 13 exact; no generated files under `phase_a_frozen/` |
| Registry parsed at runtime | yes; the recorded `manifest_sha256` matches the pinned manifest |
| Their parsed set against my own independent parse | identical, 19 and 19 |
| parsed == implemented handlers == executed records | true, proven before execution and again after |
| All nineteen mutations executed and red | yes; every record carries object mutated, predicate, a passing baseline, a failing mutated result and an observed red |
| Reproduction | clean run exits 0 from a scratch copy |

## The linkage property, tested directly rather than inferred

The previous landing's self-attacks proved that the integrity layer fires. They did not prove
that the parsed registry drives coverage once integrity has passed, which are separate
properties. That property was tested here, below the integrity layer, by altering the parsed
set after hash verification and after the qualifier's own self-test:

| Injection | Result |
| --- | --- |
| Remove one registered identifier from the parsed set | exit 1, `COVERAGE FAILURE`, naming `G-M03` as implemented but not parsed, `FATAL: Pre-execution coverage mismatch` |
| Add an identifier no handler implements | exit 1, same fatal coverage mismatch |
| Parser fed scratch text with a gate removed | returns 18 rather than 19 |
| Parser fed a duplicated identifier | raises, `Duplicate gate identifier` |
| Parser fed an empty mutation declaration | raises, `Empty mutation declaration` |

Neither failing injection is a hash failure. The run stops on coverage, which is what
establishes that the manifest parser, not a literal, is what qualification rests on.

The qualifier's own parser-to-coverage self-test is retained, and no frozen byte was touched
by any of this testing: `phase_a_frozen/` holds thirteen files and no bytecode afterward.
