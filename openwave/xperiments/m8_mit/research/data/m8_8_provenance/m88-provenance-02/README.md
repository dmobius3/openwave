# M8.8 construction-packet provenance archive, `M88-CONSTR-02`

This archive is the derivation source and environment record for the M8.8 construction
packet, delivered to the maintainer before the § 8 lock commit and published in full at
commitment. If you are reading it after publication, this is the material that lets you
rerun the derivation of the packet and check that nothing in it could have been steered.

| | |
| --- | --- |
| construction packet SHA-256 | `df00c0222f98c481eb56b882cd867a6c3a4f8604b8633e81dec0cce1f8460a06` |
| derivation source SHA-256 | `8a3a1c87f54372a446356a5c2a5ece4d9b4ba7a32367ef129b8baf18b44733f6` |
| M8.5-A group packet SHA-256 | `e3b0c945bbbb15b4549fa641234c9461062c2337b3d1e372af621b614d4883a9` |
| provenance class | `derived` (§ 4.2 of the enclosed protocol) |
| supersedes | `M88-CONSTR-01`, packet `2b51ce55…`, REJECTED by the maintainer audit |

The packet hash is recorded in the review record of the PR that landed the packet, and is
pinned in the protocol's § 11 at the lock commit; before that commit the § 11 row reads
`[PIN at landing]` by design. The group-packet hash is published in the platform repository
as `authoritative_sha256` in `research/data/m8_5a_packet_audit.json`, which is a repository
path, not a file in this archive.

## Where to look

- **`MANIFEST.md`**: every file with its SHA-256 and byte count, in three classes; only the
  six FROZEN files define `source_content_sha256`. Recompute rather than read.
- **`ENVIRONMENT.md`**: the concatenation rule, tested verification and rebuild recipes, the
  maintainer-audit invocation, the exact reproducibility statement, and why `certify.py`
  travels with the archive. Every command in it was executed as written before being
  documented.
- **`build/audit_A1_A11_output.txt`**: the maintainer audit's full transcript, A1 through
  A11 plus the mutation suite, regenerated from inside this archive; its header carries the
  packet hash it ran against.
- **`build/author_side/construction_audit.md`**: the author-side account, including the
  candidate-1 rejection, the corrections made along the way, and what this archive does NOT
  establish.
- **`recover.py`**: an author-supplied self-check. Its verdict is informational only; the
  independent recipes in `ENVIRONMENT.md` are the audit path.

## The frozen six

`source_content_sha256` covers exactly six files, concatenated in this order and nothing
else: `build/qphi.py`, `build/complex.py`, `build/kernel.py`, `build/sat.py`,
`build/search_sym.py`, `build/build_packet_v2.py`. Editing any of them, including whitespace
or a comment, invalidates the published hash. Everything else here is execution, recovery,
or context, classified file by file in `MANIFEST.md`.

The frozen sources contain no absolute filesystem paths, so the derivation runs from any
location; paths inside them are relative to their working directory, and the recipes in
`ENVIRONMENT.md` state the working directory each command needs. (`recover.py` scans every text
file in the archive, itself included, for machine-specific paths, meaning user and home
directories; deliberate portable paths in the recipes, like `/tmp` scratch directories and
the container mount, are visible to a human read and are not defects.)

## Candidate 1, for the record

The first construction packet passed all 31 author-side gates and its own 8-mutation battery
and was wrong: `im ∂₃` sat at finite index inside `ker ∂₂`, caught by the maintainer audit's
universal-cover check (A7). Its bytes are retained by the author outside this archive, with
their own rejection record; nothing in THIS archive verifies the candidate-1 figures, which
are recorded in `construction_audit.md` as history. The rejection is part of why this archive
exists in the form it does: the audit that caught it could not have been run by the author's
own gate set.

## Clean-room status

`build/author_side/` stays outside the clean room until commitment. (§ 4.2 check 6 governs
the maintainer-side construction-audit artifact; this author-side account is delivered as
input to that held record, and both stay out of the room on the same schedule.) The
construction packet itself is public. The `m8_8_reproduction_protocol.md` here is a SNAPSHOT
of the protocol as of delivery, not a pin; the live text is the PR until the lock commit
freezes it.
