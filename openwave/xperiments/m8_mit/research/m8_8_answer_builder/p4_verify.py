"""Hash verification, in a module that cannot parse JSON.

WHY THIS IS A SEPARATE FILE. § 4.4 requires the consumer to verify the hash
BEFORE parsing. Six successive attempts to demonstrate that by measurement each
failed against a parse route the instrumentation did not name, the last one by
a `JSONDecoder` instantiated at import time: wrapping the class, the scanner
factories and the module's own default decoder cannot retroactively replace the
`scan_once` already bound on an instance that existed first. Enumerating parser
callables will always admit one more.

So the ordering is no longer observed. It is STRUCTURAL. Verification lives here,
this module imports `hashlib` and `pathlib` and nothing else, and it holds no
reference to any parser. `p4.py` runs it in a subprocess where `json` and `_json`
are unimportable: if the bytes verify and a wrong hash is refused in an
interpreter that has no JSON implementation at all, no JSON parse can have
happened during verification, whatever any instance anywhere was holding.

The loader in `p4_reference.py` then composes the two in one expression, and the
control checks that shape on the syntax tree rather than trusting it. Nothing may
run before the verification call, because there is nowhere to put it.

Keep the import list here at `hashlib` and `pathlib`. Adding a third import is
not a style question: it is the guarantee.
"""
import hashlib
import pathlib


class HashMismatch(Exception):
    """Raised before any caller gets bytes back, which is the whole point."""


def verified_bytes(path, expected_sha256):
    """Return the file's bytes, or raise. No parser is reachable from here.

    ONE ENTRY SHAPE. This took a `bytes` path too, as a convenience, and a
    reviewer raised it twice. The arm was never unsafe: if `raw` hashes to
    `expected_sha256` then `raw` IS the canonical bytes whichever arm produced
    it, and the protocol rests on that resistance everywhere else. The reason it
    is gone is that it was unreachable, and that it cost an allowlist entry.
    `isinstance` sat in the caller's permitted-call set solely to admit it, so
    removing the arm turns that permission into a prohibition on the one file
    whose entire value is being pinnable."""
    raw = pathlib.Path(path).read_bytes()
    got = hashlib.sha256(raw).hexdigest()
    if got != expected_sha256:
        raise HashMismatch(f"{got[:16]} is not the expected {expected_sha256[:16]}; "
                           "refused before parsing")
    return raw
