# M8.8 § 8 step 5: the clean-room output, landed unread against the answer packet

> **What this commit is.** The implementer's complete deliverable set from the § 8 step 3
> clean room, copied out byte-identical and committed BEFORE the canonical answer packet was
> opened or its hash verified. This is the ordering record § 3 requires, and it is carried by
> commit ancestry rather than by anyone's word: this commit exists in the history before any
> commit that unseals a quarantined object.
>
> **This commit claims no result.** No comparison has run, no answer packet has been opened,
> and no verdict of any kind is recorded here. The § 8 steps 6 through 9 sequence is entirely
> ahead of this point.

## What landed

| File | SHA-256 | Role |
| --- | --- | --- |
| `METHOD_AND_GATE_MANIFEST.md` | `e7655a397445ac689202a43f1c2cf8f1582f39828a5037131136b8d6e3353e64` | § 8 step 4: route, conventions, gates, declared native orientation |
| `m88_torsion.py` | `171e346ae0259700cbe2ff1ddf58f37ecd2f080763645d2e40304cf656d0cd61` | the implementation |
| `m88_gates.py` | `a8613bcdd22f0f0ac25114578c09baa511baa4c6ae860693bff220f310016835` | the gate and mutation battery |
| `ENVIRONMENT.md` | `f894579a7b797544aa9516d910689b84ef4a6e84b5fa2ecb4bea98d767e71a10` | environment record |
| `RAW_OUTPUT.json` | `466852df9e66e857198186d050564b3276c05eb5cc8143127a389a9ef65a3c72` | the § 5.5 raw output |
| `CONSULTED_FILES.md` | `b95eb2cfaa2a58aee9d6473ac5e6d28363794cd4a11c6ec2eac257e78f07ae9e` | consulted-files manifest, generated inside the room |
| `gate_results.json` | `cd4f833d8f4466d08e145f5ef76101c75fab551b1935e887b94fd94a44761007` | machine-readable gate outcomes |

Copied verbatim. No reformatting was applied, deliberately: the bytes are the object, and a
formatter pass would move every hash above. The room's compiled-bytecode directory is not an
artifact and was excluded.

## The room, and that it held

The room opened from the § 8 step 2 lock commit with exactly four files. All four were
verified byte-unchanged at the end of the run against the hashes frozen in the lock manifest:
the protocol at the content commit, the group packet, the construction packet, and the
operational task file. The implementer declared its § 0 eligibility before reading anything
and recorded an affirmative statement that no external references were consulted.

## Disclosure: one post-implementation edit to the step 4 manifest

The manifest was written first, before any implementation code, as § 8 step 4 requires. It was
then edited once, after the implementation existed, because its stated `SU(2)` embedding did
not match the embedding the code actually used; the implementer found the mismatch itself and
resolved it by correcting the manifest to match the code. The edit is disclosed here rather
than left to be inferred, and two consequences follow.

First, filesystem modification times inside the room are NOT the ordering evidence and must
not be read as such: the manifest's final mtime postdates the implementation file by thirteen
seconds for exactly this reason, while the true order is manifest-then-code. Second, the
amendment is a convention statement rather than a result, and the two embeddings are
relabelings of the same quaternion algebra, so the torsion values are insensitive to the
choice. That last point is stated as a claim to be checked at adjudication, not as a finding
already established.

## Author-side verification of the landed bytes

Run outside the room, on a scratch copy of the committed files plus the two public packets,
so the room itself was not touched and no manifest-listed file could be overwritten by the
audit:

| Check | Result |
| --- | --- |
| `m88_torsion.py` from the committed bytes | exit 0 |
| `m88_gates.py` from the committed bytes | exit 0, every gate reported PASS |
| `gate_results.json` regenerated | byte-identical to the committed copy |
| `RAW_OUTPUT.json` regenerated | byte-identical to the committed copy |
| Determinism | no timestamps or environment-dependent fields; a leaf-level comparison of the parsed output found 320 of 320 shared leaves equal |

**One reproduction note a rerunner needs.** `m88_torsion.py` embeds `gate_results.json` into
its output when that file is present and prints a warning when it is not, so a first run in an
empty directory produces a `RAW_OUTPUT.json` that is correct in every computed value but
carries no `gate_results` block. Run the gate battery first, or run the computation a second
time after it. Byte-identity against the committed output holds under the documented order,
and the 24 keys that differ in the other order are the embedded gate block alone, never a
computed value.

This verifies that the artifact reproduces itself from its own committed source. It is not a
comparison against anything, and it establishes nothing about agreement.

## What is deliberately absent

No adjudication record, no comparison harness output, no orientation selection, and no answer
packet. Those are § 8 steps 6 through 9 and they land later, after this commit is in the
history. The maintainer-side construction-audit artifact also remains outside the tree until
step 9, per § 4.2.
