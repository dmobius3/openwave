# M8.8 construction-packet provenance delivery, `M88-CONSTR-02`

Encrypted provenance archive for the M8.8 reproduction protocol, delivered publicly per the
agreement in [PR #408](https://github.com/openwave-labs/openwave/pull/408).

This is an ORPHAN branch. It shares no history with `main` and is not part of any base the
clean room opens from. It holds exactly two files.

| | |
| --- | --- |
| file | `m88-provenance-02.tar.gz.age` |
| form | ASCII armored `age`, 81522 bytes |
| ciphertext SHA-256 | `2ba72660c74b69d458141e9b0842e4da289408558654021c75d6782133059765` |
| recipient | `age1s8zd9xr73fkm857n44d59hyrn6f0aus56dl2luf6pgh7dmfjuuxqnfvsl5` |
| `age` version | v1.3.1, both ends |

## Revision history

| commit | plaintext | ciphertext | why |
| --- | --- | --- | --- |
| `46c30841` | `f14e68c2…` | `a5665d3e…` | first delivery, SUPERSEDED |
| this one | `4fa0228b…` | `2ba72660…` | certificate defect (maintainer-found) fixed, then the archive rebuilt under three independent cold-read rounds |

`46c30841` is retained rather than rewritten: it is the delivery the maintainer verified on
the record. Its archive shipped a `saturation_certificate.json` that asserted its containment
premises as literals and pinned a packet not in the archive. That fix and the full account of
everything the subsequent adversarial reads found are in the archive's own
`construction_audit.md` and `ENVIRONMENT.md`. `source_content_sha256` did NOT move
(`8a3a1c87…`), and the construction packet did not change (`df00c022…`).

## What it is, and what it is not

The **plaintext** archive inside is the canonical provenance object for construction packet
`df00c022…` under provenance ID `M88-CONSTR-02`. It stays unpublished until commitment, at
which point it is published and checked against the plaintext hash recorded in the PR.

The **ciphertext** published here is a transport object carrying a public timestamp. `age`
encryption is randomized, so re-encrypting the same plaintext yields different bytes; this
hash identifies this one object and nothing more.

Confidentiality now rests entirely on custody of the recipient secret key through
commitment. That is a deliberate transfer of risk, made because an independent timestamp is
worth more than a private channel, and recorded here rather than left implicit.

## Verification, after decryption

1. plaintext SHA-256 `4fa0228bc7c99bca0770399c82bc9981f5f3c934c608773c77ddb798e5ad0913`
2. all 19 manifest hashes, recomputed independently rather than via the archive's own
   `recover.py`, whose verdict is informational only
3. the frozen-six concatenation reproducing
   `8a3a1c87f54372a446356a5c2a5ece4d9b4ba7a32367ef129b8baf18b44733f6` from inside the
   extracted tree
