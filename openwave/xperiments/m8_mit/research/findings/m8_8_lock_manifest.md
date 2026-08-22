# M8.8 § 8 step 2: the LOCK MANIFEST

> The commit that adds this file is the LOCK COMMIT required by
> [`m8_8_reproduction_protocol.md` § 8 step 2](m8_8_reproduction_protocol.md). Per § 10.3 it
> is the SOLE freeze boundary: from this commit the addenda-only rule of § 12 binds, whether
> or not anything has merged. Any rebase, force-push or other history rewrite after it VOIDS
> the freeze and requires a new lock record and a clean-room restart; both § 8 commits are
> preserved through annotated tags (`m8.8-content`, `m8.8-lock`) so no merge strategy can
> orphan them.

## The content commit

| Field | Value |
| --- | --- |
| § 8 step 1 CONTENT COMMIT | `b1b6ce48ac54a2b861782dc67487c7b5cf06b745` |
| Its record | [`m8_8_content_commit.md`](m8_8_content_commit.md), in that commit's tree |

## Every packet and audit hash (SHA-256)

| Object | SHA-256 |
| --- | --- |
| Group packet [`../data/m8_5a_packet.json`](../data/m8_5a_packet.json) | `e3b0c945bbbb15b4549fa641234c9461062c2337b3d1e372af621b614d4883a9` |
| Group-packet re-audit [`../data/m8_8_group_packet_reaudit.json`](../data/m8_8_group_packet_reaudit.json) | `6e039277fbba121a81d1ee2bad60784640432964e9f21d16ce0089c0ab57be2e` |
| Construction packet [`../data/m8_8_construction_packet.json`](../data/m8_8_construction_packet.json) | `df00c0222f98c481eb56b882cd867a6c3a4f8604b8633e81dec0cce1f8460a06` |
| Maintainer-side construction-audit artifact (outside the tree until § 8 step 9) | `d5bb04b9c747d3780a3e931d33d8ed9c7ab79c759a9626ccf664d88b478ef0bb` |
| Construction-provenance archive, ciphertext (tag `m8.8-provenance-02`; plaintext unpublished until commitment, its hash on the [#408](https://github.com/openwave-labs/openwave/pull/408) record) | `2ba72660c74b69d458141e9b0842e4da289408558654021c75d6782133059765` |
| Canonical answer packet (quarantined with the author until § 8 step 6) | `744c7f25e2312d90fc356b11510da685328f05e80ae62721d0a0f418dcf9697e` |
| M8.3 method note, the § 1 analytic-side control (§ 11 pin) | `3e0c1901d4089991a1de7cff0b1cde453257a29891793249ce63955f144ef06d` |
| The protocol text frozen by this lock (the version at the content commit) | `23f313d0cd47e4ff644ae7c9730cbc9eb380cb6f7c64e9daddb7b98cd6885d87` |

## The clean-room base

| Field | Value |
| --- | --- |
| Base the clean room opens from | EXACTLY the lock commit, the commit adding this file (§ 8 step 3); its identifier is carried by the annotated tag `m8.8-lock`, since a commit cannot contain its own hash |
| Clean-room standard | the M8.5-A clean-room standard at commit `e53a64d493c33324318c0b5b3007f566f4d82f5d` (§ 11 pin) |

## State after this commit

The maintainer-side sequence of § 8 steps 1 and 2 is complete. Next in the frozen sequence:
the author commissions the fresh implementer (§ 4 makes the author permanently ineligible to
serve as one), the clean room opens from the lock commit, and steps 4 and 5 are the
implementer's commits. The § 12 addenda-only rule binds from here.
