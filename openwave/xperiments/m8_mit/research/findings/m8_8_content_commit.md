# M8.8 § 8 step 1: the content commit record

> The commit that adds this file is the CONTENT COMMIT required by
> [`m8_8_reproduction_protocol.md` § 8 step 1](m8_8_reproduction_protocol.md). The protocol,
> the group-packet re-audit and the construction packet are recorded by containment: all
> three are in this commit's tree, at the hashes below. The two objects that must NOT enter
> the tree, the maintainer-side construction-audit artifact (§ 4.2 item 6) and the canonical
> answer packet (§ 4.4), are frozen here by hash only.

## In this commit's tree

| Object | Path | SHA-256 |
| --- | --- | --- |
| The protocol (settled text per [#430](https://github.com/openwave-labs/openwave/pull/430)) | [`m8_8_reproduction_protocol.md`](m8_8_reproduction_protocol.md) | `23f313d0cd47e4ff644ae7c9730cbc9eb380cb6f7c64e9daddb7b98cd6885d87` |
| Construction packet (equals the § 11 pin) | [`../data/m8_8_construction_packet.json`](../data/m8_8_construction_packet.json) | `df00c0222f98c481eb56b882cd867a6c3a4f8604b8633e81dec0cce1f8460a06` |
| Group packet (equals the § 11 pin) | [`../data/m8_5a_packet.json`](../data/m8_5a_packet.json) | `e3b0c945bbbb15b4549fa641234c9461062c2337b3d1e372af621b614d4883a9` |
| Group-packet re-audit, run against M8.8's own § 4.3 forbidden list per the [task record](../tasks/m8_8_task_details.md), not inherited from M8.5-A | [`../data/m8_8_group_packet_reaudit.json`](../data/m8_8_group_packet_reaudit.json) | `6e039277fbba121a81d1ee2bad60784640432964e9f21d16ce0089c0ab57be2e` |

## Frozen by hash, outside the tree

| Object | SHA-256 | Custody until publication |
| --- | --- | --- |
| Maintainer-side construction-audit artifact `m8_8_packet_audit.json`, 9463 bytes: 12 of 12 checks pass, 12 of 12 mutations detected, produced by [`m8_8_packet_audit.py`](../scripts/m8_8_packet_audit.py) at `e3165904` against the construction packet above | `d5bb04b9c747d3780a3e931d33d8ed9c7ab79c759a9626ccf664d88b478ef0bb` | maintainer, on an orphan commit outside every base the clean room can open from; PUBLISHED with the adjudication evidence at § 8 step 9 and verified against this hash |
| Canonical answer packet (the § 11 pin, recorded by the author; issued and frozen 2026-08-10) | `744c7f25e2312d90fc356b11510da685328f05e80ae62721d0a0f418dcf9697e` | author, QUARANTINED per § 3 until § 8 step 6 |

The construction packet's provenance class is `derived`. Its provenance record is the
archive adversarially verified on the [#408](https://github.com/openwave-labs/openwave/pull/408)
record (2026-08-05): ciphertext SHA-256
`2ba72660c74b69d458141e9b0842e4da289408558654021c75d6782133059765` (81522 bytes), anchored by
the tag `m8.8-provenance-02`; the plaintext archive stays unpublished until commitment, its
hash on the same thread record.

The LOCK COMMIT of § 8 step 2 follows as a separate commit; its manifest, never this file,
carries the commit identifiers, because a commit cannot contain its own hash.
