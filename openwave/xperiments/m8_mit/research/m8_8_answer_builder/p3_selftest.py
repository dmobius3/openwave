"""P3 ENGINEERING SELF-TEST. Not protocol evidence.

WHAT THIS IS. A development prerequisite for trusting the P3 gate: proof that
each of its twelve checks is reachable and can reject something, before anyone
relies on a green result. A check that cannot fail is not a check, and this file
exists so that P3's green run means something when it is first read.

WHAT THIS IS NOT. It is NOT the plan's mutation battery. `ANSWER_PACKET_PLAN.md`
assigns that to P4: "every P3 check reddened by a targeted mutation", alongside
the § 4.4 ingestion controls. That is the formal adversarial evidence and it is
P4's deliverable, deliberately produced by a different step than the one that
wrote the gate. This file and P4's content-mutation half will cover overlapping
ground, and they have different evidentiary standing on purpose.

The two are named differently for a reason this build learned the hard way: when
one object quietly acquires two evidentiary jobs, a claim about one of them ends
up resting on work done for the other.

Everything here runs through P3's INJECTED mode. The production hash pin does not
apply and is not consulted, which is the point of having the two modes: the pin
is input integrity, and this is about whether the checks work.
"""

import copy
import hashlib
import json
import pathlib
import sys

import p3

CASES = []


def case(target, note):
    def deco(fn):
        CASES.append((fn.__name__, target, fn, note))
        return fn
    return deco


def _sidecar_for(pb):
    """A sidecar committing to exactly these packet bytes, so a packet mutation
    does not incidentally redden the sidecar checks."""
    return {
        "_banner": SIDECAR["_banner"],
        "record_format_version": SIDECAR["record_format_version"],
        "status": SIDECAR["status"],
        "status_note": SIDECAR["status_note"],
        "answer_packet_format_version": SIDECAR["answer_packet_format_version"],
        "canonicalization_rule": SIDECAR["canonicalization_rule"],
        "canonical_plaintext_byte_length": len(pb),
        "expected_canonical_plaintext_sha256": hashlib.sha256(pb).hexdigest(),
    }


def mutate_packet(fn):
    """Apply a mutation to the packet object, re-canonicalize, rebuild the
    sidecar to match. Returns (packet_bytes, sidecar_bytes)."""
    pkt = copy.deepcopy(PACKET)
    fn(pkt)
    pb = p3.canonical_bytes(pkt)
    return pb, p3.canonical_bytes(_sidecar_for(pb))


# ------------------------------------------------------------------- cases

@case("KEYS_EXACTLY_SEVEN", "an eighth top-level key")
def extra_key():
    return mutate_packet(lambda p: p.update({"notes": "extra"}))


@case("PACKET_CANONICAL_BYTES", "four-space indent: parses the same, § 10.2 says two")
def packet_indent_wrong():
    bad = (json.dumps(PACKET, sort_keys=True, indent=4, ensure_ascii=True) + "\n").encode()
    return bad, p3.canonical_bytes(_sidecar_for(bad))


@case("SIDECAR_CANONICAL_BYTES", "CR in the sidecar; § 10.2 requires LF")
def sidecar_has_cr():
    pb = p3.canonical_bytes(PACKET)
    return pb, p3.canonical_bytes(_sidecar_for(pb)).replace(b"\n", b"\r\n", 1)


@case("SIDECAR_COMMITS_TO_PACKET", "sidecar naming a different hash")
def sidecar_wrong_hash():
    pb = p3.canonical_bytes(PACKET)
    rec = _sidecar_for(pb)
    h = rec["expected_canonical_plaintext_sha256"]
    rec["expected_canonical_plaintext_sha256"] = h[:-1] + ("0" if h[-1] != "0" else "1")
    return pb, p3.canonical_bytes(rec)


@case("TRIPLES_NORMALIZED", "a triple scaled by 2: same number, not the frozen encoding")
def unnormalized():
    def f(p):
        a, b, c = p["rows"][0]["value"]
        p["rows"][0]["value"] = [a * 2, b * 2, c * 2]
    return mutate_packet(f)


@case("CLASS_CENSUS", "the declared class moved onto a free row")
def census_wrong():
    def f(p):
        for r in p["rows"]:
            if r["evidentiary_class"] == p3.CLASS_FREE:
                r["evidentiary_class"] = p3.CLASS_DECLARED
                return
    return mutate_packet(f)


@case("ADJUDICATES_MATCH_PUBLIC_PACKETS", "one hex digit of the group packet hash")
def adjudicates_wrong():
    def f(p):
        h = p["adjudicates"]["group_packet_sha256"]
        p["adjudicates"]["group_packet_sha256"] = h[:-1] + ("0" if h[-1] != "0" else "1")
    return mutate_packet(f)


@case("SIGNATURES_DISTINCT", "two rows sharing a signature: § 5.5 must separate all nine")
def duplicate_signature():
    def f(p):
        p["rows"][1]["signature"] = copy.deepcopy(p["rows"][0]["signature"])
    return mutate_packet(f)


@case("DECLARED_ROW_CONVENTION", "R0 no longer declaring T^2 = 1")
def declared_row_wrong():
    def f(p):
        for r in p["rows"]:
            if r["label"] == p3.DECLARED_ROW:
                r["value"] = [2, 0, 1]
    return mutate_packet(f)


@case("IDENTITIES_RECOMPUTE_POSITION_WISE", "one slot's expected value taken from another")
def identity_drifts():
    def f(p):
        p["identities"][0]["expected_value"] = list(p["identities"][2]["expected_value"])
    return mutate_packet(f)


@case("VALUES_MATCH_BOTH_SOURCE_FAMILIES", "two row values exchanged: both still normalized")
def values_swapped():
    def f(p):
        p["rows"][0]["value"], p["rows"][8]["value"] = \
            p["rows"][8]["value"], p["rows"][0]["value"]
    return mutate_packet(f)


@case("SELECTOR_NOT_SELF_INVERSE", "a selector equal to 1: § 5.4's invalid-anchor branch")
def selector_self_inverse():
    def f(p):
        for r in p["rows"]:
            if r["label"] == p3.SELECTOR_ROW:
                r["value"] = [1, 0, 1]
    return mutate_packet(f)


# -------------------------------------------------------------------- main

def main():
    if len(sys.argv) != 6:
        print("usage: p3_selftest.py <p1> <packet> <sidecar> <group> <construction>",
              file=sys.stderr)
        return 2
    global PACKET, SIDECAR
    p1 = json.loads(pathlib.Path(sys.argv[1]).read_bytes())
    PACKET = json.loads(pathlib.Path(sys.argv[2]).read_bytes())
    SIDECAR = json.loads(pathlib.Path(sys.argv[3]).read_bytes())
    gb = pathlib.Path(sys.argv[4]).read_bytes()
    cb = pathlib.Path(sys.argv[5]).read_bytes()

    print("P3 ENGINEERING SELF-TEST. Not protocol evidence; the formal\n"
          "'every P3 check reddens' claim is P4's, per the frozen plan.\n")

    base_pb = p3.canonical_bytes(PACKET)
    base_sb = p3.canonical_bytes(SIDECAR)
    checks, stops = p3.gate(base_pb, base_sb, p1, gb, cb)
    if stops:
        print(f"  BASELINE IS NOT GREEN: {stops[:2]}")
        return 1
    print(f"  baseline: {len(checks)} checks, all green\n")

    observed, failures = set(), []
    for name, target, fn, note in CASES:
        pb, sb = fn()
        ch, _ = p3.gate(pb, sb, p1, gb, cb)
        red = {c for c, ok, _ in ch if not ok}
        detail = {c: d for c, ok, d in ch if not ok}
        hit = target in red and not str(detail.get(target, "")).startswith("BLOCKED:")
        if hit:
            observed.add(target)
        else:
            failures.append(name)
        extra = sorted(red - {target})
        print(f"  {'RED ' if hit else 'MISS'}  {name}")
        print(f"          {note}")
        print(f"          target {target}" + (f"; also red: {extra}" if extra else ""))

    never = [c for c in p3.CHECK_IDS if c not in observed]
    print(f"\n  {len(CASES) - len(failures)}/{len(CASES)} cases reddened their target")
    print(f"  reachability: {len(observed)}/{len(p3.CHECK_IDS)} P3 checks observed rejecting")
    if never:
        print(f"  NEVER OBSERVED REJECTING: {never}")
    ok = not failures and not never
    print(f"  {'every P3 check is reachable' if ok else 'REACHABILITY INCOMPLETE'}")
    print("\n  Reminder: this is an engineering prerequisite, not the protocol's\n"
          "  mutation evidence. P4 owns that claim.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
