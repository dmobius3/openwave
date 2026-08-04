# M8.8 construction-packet provenance delivery, `M88-CONSTR-02`

Encrypted provenance archive for the M8.8 reproduction protocol, delivered publicly per the
agreement in [PR #408](https://github.com/openwave-labs/openwave/pull/408).

This is an ORPHAN branch. It shares no history with `main` and is not part of any base the
clean room opens from. It holds exactly two files.

| | |
| --- | --- |
| file | `m88-provenance-02.tar.gz.age` |
| form | ASCII armored `age`, 79962 bytes |
| ciphertext SHA-256 | `a5665d3e1f677a7bf89f492544e81d0b66a7a92e9ab97375c39c4158f0d768eb` |
| recipient | `age1s8zd9xr73fkm857n44d59hyrn6f0aus56dl2luf6pgh7dmfjuuxqnfvsl5` |
| `age` version | v1.3.1, both ends |

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

1. plaintext SHA-256 `f14e68c2d068043b8c4c5b11297b5b9d37cf0638f015cc9b7a7014043be77518`
2. all 21 manifest hashes, recomputed independently rather than via the archive's own
   `recover.py`, whose verdict is informational only
3. the frozen-six concatenation reproducing
   `8a3a1c87f54372a446356a5c2a5ece4d9b4ba7a32367ef129b8baf18b44733f6` from inside the
   extracted tree
