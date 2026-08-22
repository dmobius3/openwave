"""P6: issue the IMMUTABLE pre-lock hash record.

The candidate record is NOT replaced or deleted. It stays as evidence of the
pre-issuance state P2 through P5 actually reviewed. The issued record is a new
artifact, and from it the packet and the record do not move: a change after this
is an explicit reissue, not a revision.

DERIVED, NEVER TRANSCRIBED. Every value below is computed from the packet bytes
this script reads. The one thing this build has caught in itself more than any
other is a value typed in one place and compared against the typing, so the hash
and the byte length are taken from the file and nothing is copied from the
candidate record except by explicit comparison against what was recomputed.

The status transition is made VISIBLE rather than done in place. The issued
record carries its own `record_format_version`, its own `status`, and names the
candidate it supersedes by hash. Reusing a field whose meaning changes under
issuance would hide exactly the transition this step exists to mark.
"""
import argparse
import hashlib
import json
import pathlib
import sys

import p2_schema as S

RECORD_FORMAT_VERSION = "m8_8-issuance-record-1"

BANNER = (
    "ISSUED, ANSWER-BEARING. IMMUTABLE from this record onward.\n"
    "This record and the packet it names do not move. A change to either is an\n"
    "explicit REISSUE, not a revision, and requires a new issuance record.\n"
    "It supersedes the CANDIDATE pre-lock hash record for every operational\n"
    "purpose; that candidate is retained as history and governs nothing.\n"
    "Published at the § 8 step 9 adjudication record, with its bytes verified\n"
    "against the frozen hash.\n"
)


def build(packet_path, candidate_record_path):
    """The issued record, every field derived from the bytes on disk."""
    raw = pathlib.Path(packet_path).read_bytes()
    digest = hashlib.sha256(raw).hexdigest()

    # § 10.2 is re-verified HERE rather than assumed from an earlier gate. An
    # issuance that froze bytes it had not itself checked would be freezing
    # someone else's verdict.
    reserialized = S.canonical_bytes(json.loads(raw.decode("ascii")))
    if reserialized != raw:
        raise ValueError("the packet is not its own canonical serialization; "
                         "§ 10.2 does not hold and nothing is issued")

    packet = json.loads(raw.decode("ascii"))
    cand_raw = pathlib.Path(candidate_record_path).read_bytes()
    cand = json.loads(cand_raw.decode("ascii"))

    # The candidate is COMPARED against what was recomputed, never copied from.
    for field, computed in (("expected_canonical_plaintext_sha256", digest),
                            ("canonical_plaintext_byte_length", len(raw))):
        if cand.get(field) != computed:
            raise ValueError(f"the candidate record's {field} is "
                             f"{cand.get(field)!r}, recomputed {computed!r}; "
                             "issue nothing until that is understood")
    if cand.get("status") != "CANDIDATE":
        raise ValueError(f"the record being superseded says {cand.get('status')!r}, "
                         "not CANDIDATE")

    return {
        "_banner": BANNER,
        "record_format_version": RECORD_FORMAT_VERSION,
        "status": "ISSUED",
        "answer_packet_format_version": packet["format_version"],
        "expected_canonical_plaintext_sha256": digest,
        "canonical_plaintext_byte_length": len(raw),
        "canonicalization_rule": cand["canonicalization_rule"],
        "supersedes": {
            "record_format_version": cand["record_format_version"],
            "status": cand["status"],
            "sha256": hashlib.sha256(cand_raw).hexdigest(),
        },
    }


def verify(record_bytes, packet_path):
    """Read the issued record back and check it against the packet."""
    rec = json.loads(record_bytes.decode("ascii"))
    raw = pathlib.Path(packet_path).read_bytes()
    problems = []
    if rec["expected_canonical_plaintext_sha256"] != hashlib.sha256(raw).hexdigest():
        problems.append("the issued hash is not the packet's hash")
    if rec["canonical_plaintext_byte_length"] != len(raw):
        problems.append("the issued byte length is not the packet's length")
    if rec["status"] != "ISSUED":
        problems.append(f"status is {rec['status']!r}")
    if S.canonical_bytes(rec) != record_bytes:
        problems.append("the record is not its own canonical serialization")
    return problems


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--packet", required=True)
    ap.add_argument("--candidate-record", required=True)
    ap.add_argument("--emit", metavar="PATH")
    a = ap.parse_args()

    rec = build(a.packet, a.candidate_record)
    out = S.canonical_bytes(rec)
    problems = verify(out, a.packet)

    print("P6, issuance\n")
    print(f"  packet          {rec['expected_canonical_plaintext_sha256']}")
    print(f"  byte length     {rec['canonical_plaintext_byte_length']}")
    print(f"  record          {hashlib.sha256(out).hexdigest()}")
    print(f"  supersedes      {rec['supersedes']['sha256']}  ({rec['supersedes']['status']})")
    print()
    for p in problems:
        print(f"  PROBLEM: {p}")
    if problems:
        return 1
    print("  the issued record is canonical and binds the packet it names.")

    if a.emit:
        path = pathlib.Path(a.emit)
        if path.exists():
            print(f"\n  REFUSED: {path} exists. An issued record is immutable; "
                  "overwriting one is a reissue and is not this script's job.")
            return 1
        path.write_bytes(out)
        print(f"\n  wrote {path}")
        print("  the candidate record is untouched and retained as history.")
    else:
        print("\n  dry run; nothing written. Pass --emit to issue.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
