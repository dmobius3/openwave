# M8.8 construction-packet provenance archive, M88-CONSTR-02, plaintext

> The canonical provenance object for the `derived` construction packet. Its ciphertext was
> delivered on 2026-08-05 and is anchored by the tag `m8.8-provenance-02`; the plaintext was
> committed by hash on [#408](https://github.com/openwave-labs/openwave/pull/408) and stayed
> unpublished until commitment, per the content-commit record. Published here after § 8
> step 9 ([#458](https://github.com/openwave-labs/openwave/pull/458)).

| Object | SHA-256 | Bytes |
| --- | --- | --- |
| `m88-provenance-02.tar.gz`, the committed plaintext | `4fa0228bc7c99bca0770399c82bc9981f5f3c934c608773c77ddb798e5ad0913` | 59948 |
| its ciphertext, tag `m8.8-provenance-02`, commit `65189763` | `2ba72660c74b69d458141e9b0842e4da289408558654021c75d6782133059765` | 81522 |

The tarball is the commitment: the plaintext hash above is the value posted on #408 before
delivery, and the maintainer verified the decrypted bytes against it on 2026-08-05. The
directory `m88-provenance-02/` beside it is extracted from those exact bytes for reading;
verify against the tarball, not the tree.

What the archive establishes is stated in its own `MANIFEST` and was fixed on the thread:
REPRODUCIBLE DERIVATION AND STRUCTURAL VERIFICATION of the construction packet, under the
protocol's `derived` provenance class (§ 4.2). It does not establish faithful extraction
from literature, independent discovery of the topological model, or uniqueness of the
derived complex.
