# M8.8 clean-room attempt 3: the output, landed unread against the answer packet

> **Claims no result.** The answer packet has not been opened, its hash has not been verified,
> no comparison has run, and no verdict exists. This commit discharges § 3's requirement that
> the implementation, environment record, derivation artifacts and raw output be committed
> BEFORE any quarantined object is unsealed, with the ordering carried by commit ancestry.

## What landed

| File | SHA-256 | Role |
| --- | --- | --- |
| `TASK.md` | `75cef74d3c325690c77d36f48d20977860505caca81c8b0d4f48e86a47beb899` | the operational instruction, committed BEFORE this run started |
| `METHOD_AND_GATE_MANIFEST.md` | `22a7ec59b1be17a9960d5c3e66a14faa530e950d80511676921c5663214b8d38` | § 8 step 4, finalized 11:52 carrying the explicit freeze declaration |
| `validate_manifest.py` | `83a6ac1a0967d0461739edcfcac049dd014c71eb7f569af3c70a4e7ecc7e9810` | pre-implementation validation artifact, 11:44, writes nothing |
| `reproduce.py` | `1e4471b3390b4b07b5cf83d436c18fb58a65ea0257b6dfe837716a0c68fcb28e` | production implementation, 12:27, begun after the freeze |
| `ENVIRONMENT.md` | `f401169784ee75f097856590416bbb61c48289eb3c86198c252ff83f7cbec2d2` | environment record |
| `DERIVATION_ARTIFACTS.json` | `60b63ec577ca40ad0740c402c35bd06dddcc9a77d7f6631b83fa995097e33bad` | the § 7 derivation artifacts |
| `RAW_OUTPUT.json` | `600dc6b9756bd0a276233e25a93d29e314494b9a077ab4acbaea616e85331fb1` | the § 5.5 raw output |
| `CONSULTED_FILES.md` | `964a557f1d1e85fc5a920a659bf30e4a24f627906a2a83a070cb0cf174502880` | consulted-files manifest, generated in the room |

## Both recorded failure modes are absent, and the new gate is visibly what prevented them

The causal sequence the revised instruction demands is present in the artifacts themselves,
in order:

```text
validate_manifest.py  11:44  ->  MANIFEST FINAL  11:52  ->  reproduce.py  12:27  ->  outputs  12:27
```

| Requirement | Evidence |
| --- | --- |
| Manifest validated before freezing | `validate_manifest.py` predates the manifest's finalization by eight minutes, and the manifest's § 6 lists it as a validation artifact with its checks and their results |
| The freeze is an explicit event, not an inference | the manifest carries the required line `MANIFEST STATUS: FINAL; pre-implementation validation complete.` verbatim |
| Validation stayed inside its boundary | `validate_manifest.py` performs no file writes at all, and its own header states that it does not compute torsion values or populate the raw output. The four output files are written by the production implementation alone |
| Production implementation began only after the freeze | `reproduce.py` postdates the finalized manifest by 35 minutes |
| The manifest was not amended after implementation began | its bytes predate every production artifact; the attempt-1 failure mode did not recur |
| No departure from the frozen manifest | no mismatch or stop report exists, so the attempt-2 failure mode did not arise either |

## Author-side verification of the landed bytes

Run outside the room on a scratch copy, so no manifest-listed file could be overwritten by the
audit itself:

| Check | Result |
| --- | --- |
| `validate_manifest.py` from the committed bytes | exit 0, `Overall: ALL PASS` |
| `reproduce.py` from the committed bytes | exit 0, about seven seconds |
| `RAW_OUTPUT.json` regenerated | byte-identical |
| `DERIVATION_ARTIFACTS.json` regenerated | byte-identical |
| `ENVIRONMENT.md` regenerated | byte-identical |
| `CONSULTED_FILES.md` regenerated | byte-identical |
| A9 intent scan | clean: standard library only, no network, no subprocess, no eval, no pickle |
| Firewall | all four seeded inputs byte-unchanged at the end of the run |

The run reports 55 gate results including a mutation arm for each of M1, M2, M6, M7, D1, D2
and T2, plus a convention fixture with module-side, vector-convention and boundary-direction
mutation instances. Whether those gates are individually sufficient is a review question and
is not asserted here.

## What is deliberately absent

No adjudication record, no comparison output, no orientation selection, and no answer packet.
Those are § 8 steps 6 through 9 and they land after this commit is in the history.
