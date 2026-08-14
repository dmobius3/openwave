"""Staged, auditable ingestion of a canonical hash-pinned reference.

The stages are kept separate because the gates downstream must be creditable
independently of the ones upstream:

    1. canonical bytes are hash-verified
    2. bytes are parsed
    3. the map / reference is read from the PARSED object
    4. the comparison indexes THROUGH that map
    5. negative arms mutate an IN-MEMORY COPY, never the bytes

Without stage separation, "the transcription mutation was caught" can silently
mean "the SHA-256 changed", which tests nothing about the comparison layer.
Protocol section 4.1 requires the step-5 mutation to be "made on an in-memory
copy downstream of the completed hash check, so it tests the comparison layer
rather than SHA-256 (G5a)".

Every derived object carries the provenance of the bytes it came from, so a
reader can see that a mutated object still descends from a successful hash
verification.

Canonical form, per section 11.7: keys sorted, two-space indent, ASCII, LF,
single trailing newline.
"""

import copy
import hashlib
import json

__all__ = ["canonical_bytes", "IngestionError", "Ingested", "ingest"]


def canonical_bytes(obj):
    """The frozen canonical serialization."""
    return (json.dumps(obj, sort_keys=True, indent=2, ensure_ascii=True)
            + "\n").encode("ascii")


class IngestionError(Exception):
    """Raised when hash verification or parsing fails.  Never swallowed."""


class Ingested:
    """A parsed object that remembers how it was obtained.

    `hash_verified` stays True through mutation, and `mutations` records what was
    done to the in-memory copy.  That pair is the audit marker: a reader can see
    the source bytes verified AND that the perturbation touched only the parsed
    representation.
    """

    def __init__(self, data, sha256, byte_length, canonical_confirmed):
        self.data = data
        self.sha256 = sha256
        self.byte_length = byte_length
        self.canonical_confirmed = canonical_confirmed
        self.hash_verified = True
        self.mutations = []

    def provenance(self):
        return {
            "source_sha256": self.sha256,
            "byte_length": self.byte_length,
            "hash_verified": self.hash_verified,
            "bytes_were_canonical": self.canonical_confirmed,
            "mutations_applied_to_parsed_copy": list(self.mutations),
            "note": ("source bytes hash-verified successfully; any mutation "
                     "listed above was applied only to the parsed in-memory "
                     "representation and never to the bytes"),
        }

    def mutated_copy(self, description, mutate_fn):
        """An in-memory copy with one perturbation, downstream of the hash check.

        The copy inherits `hash_verified = True` deliberately: the bytes DID
        verify, and that is precisely what makes a red result downstream
        attributable to the comparison layer rather than to SHA-256.
        """
        clone = Ingested(copy.deepcopy(self.data), self.sha256, self.byte_length,
                         self.canonical_confirmed)
        clone.mutations = list(self.mutations) + [description]
        mutate_fn(clone.data)
        return clone


def ingest(raw, expected_sha256, require_canonical=True):
    """Stage 1 then stage 2.  Raises rather than returning a degraded object.

    A failure here must never be reported as a downstream gate catching
    something; the caller is expected to let this propagate.
    """
    if isinstance(raw, str):
        raw = raw.encode("ascii")
    got = hashlib.sha256(raw).hexdigest()
    if got != expected_sha256:
        raise IngestionError(
            f"hash mismatch: expected {expected_sha256}, computed {got}")

    try:
        data = json.loads(raw.decode("ascii"))
    except Exception as exc:
        raise IngestionError(f"canonical bytes did not parse: {exc}") from exc

    canonical_confirmed = (canonical_bytes(data) == raw)
    if require_canonical and not canonical_confirmed:
        raise IngestionError("bytes are not their own canonical serialization")

    return Ingested(data, got, len(raw), canonical_confirmed)
