# M8.5-C attempt A1: the archived record

The single commissioned Build Unit run of the M8.5-C qualification protocol. TERMINATED
2026-08-29 without adjudication: neither of the protocol's two terminal sentences issued.
The canonical account of why, and the ruling that dispositions this record, is the
protocol's dated Addendum 1 ([`../../findings/m8_5c_protocol.md`](../../findings/m8_5c_protocol.md),
below its freeze boundary), which carries the maintainer's ruling from
[#501](https://github.com/openwave-labs/openwave/issues/501) verbatim. This README is a
map, not a second account.

| file | what it is |
| --- | --- |
| `ledger/OUTPUT_LEDGER.jsonl` | the unit's append-only run ledger: 4 GATE, 4 RESOURCE, 1 STOP record; no gate-4 verdict of either color |
| `ledger/COMMISSIONING.md` | preflight: provenance checks, launch-gate output, environment, input-manifest closing hash |
| `ledger/INPUT_MANIFEST.json` | the § 12 input manifest, closed and hashed before the first GATE record |
| `A1_AUTHOR_ANNOTATION.md` | dated post-stop author annotation: two ledger-semantics corrections a mechanical reader needs, archive hashes, the denominator, post-stop provenance labels |
| `00_COMMISSION.md` | the commission (NON-GOVERNING room mechanics, archived per § 12) |
| `ROOM_MANIFEST.json` | the supplied-input whitelist with hashes, generated at commissioning |
| `room_import_gate.py` | the § 12 launch gate: per-module subprocess with sentinel, self-arming |
| `in_room/` | every line the unit wrote in-room, as left at the stop (the room's `build/`; renamed here only because `build/` is gitignored platform-wide, bytes identical) |

The unit's files are byte-identical to the stopped room. The attempt never resumes, and
nothing here carries execution credit into any successor.
