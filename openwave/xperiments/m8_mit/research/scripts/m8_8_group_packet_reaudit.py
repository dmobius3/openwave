"""M8.8 GROUP-PACKET RE-AUDIT: the M8.5-A packet re-audited against M8.8's forbidden list.

WHAT THIS IS.  The [M8.8 reproduction protocol](../findings/m8_8_reproduction_protocol.md)
section 4.1 carries the M8.5-A group packet as PUBLIC and already audited, and section 4.3
requires it RE-AUDITED against M8.8's own forbidden-input list rather than inherited on its
[M8.5-A audit record](../data/m8_5a_packet_audit.json), because "the two tasks have different
forbidden sets, and a packet clean for one is not automatically clean for the other".  Section
8 step 1 names this re-audit as part of the CONTENT COMMIT, so it runs at freeze.

THE TWO LISTS DIFFER IN BOTH DIRECTIONS, which is the whole argument for re-running:

  1  M8.5-A forbade the MCKAY-side target: dimensions, distances, character values, class
     labels.  Its scan does not carry the word `torsion` at all, so a torsion value planted
     in this packet's declared prose passes the M8.5-A audit green.  The mutation
     `torsion_word_in_declared_prose` below demonstrates exactly that, and it is the single
     clearest reason this file exists.
  2  M8.8's CONSTRUCTION-packet scan forbids `phi`, which is correct there (that packet
     carries integers and element IDs only) and wrong here (this packet's coefficient field
     IS Q(phi)).  The token is therefore ABSENT from the list below by decision rather than
     by oversight, the record states how many times `phi` legitimately occurs, and the
     invariance test `extra_phi_occurrence_in_prose` holds the decision in place.

  The exemption is expressed as an omission from the declared list, NOT as a subtraction
  applied to the scanned text.  A subtraction that silently matches nothing is how the
  construction audit's A9 exemption came to be inert for three revisions.

WHAT THIS VERIFIES.  Five checks, each mutation-tested:

  G1  identity: the delivered bytes reproduce the protocol section 11 pin, and the declared
      canonicalization is a no-op on them, so the pin names a GROUP and not a transcription
  G2  the M8.5-A structural audit re-runs GREEN on these exact bytes, A1 through A8, by
      importing that audit rather than restating it, so the two records cannot drift
  G3  no section 4.3 answer-bearing vocabulary, under M8.8's list as declared below
  G4  no evaluated numeric content anywhere: no JSON numeric literal, no decimal or
      scientific rendering.  This is the check a vocabulary scan cannot do, since a planted
      torsion value need not be spelled with any forbidden word
  G5  the canonical element enumeration is the one the construction packet's element IDs
      address: 120 elements, ranked by lexicographic order of the canonical coordinate
      tuple per the section 4.2 encoding rule, pinned by its own digest under the
      serialization `ENUMERATION_ENCODING` declares and the artifact carries

WHY G5 IS THIS RE-AUDIT'S ADDITION.  M8.5-A consumed the group as a group; nothing in it
depended on WHICH element received which index.  M8.8's construction packet addresses
elements by canonical ID, so the enumeration is load-bearing between the two public packets:
the same 120 quaternions in a different order would silently change what every boundary-map
entry means.  The enumeration is a function of the packet alone, which G5 records as a
frozen digest and the invariance test `generator_order_reversed` shows is independent of the
order the packet happens to list its generators in.

WHY THE ENCODING IS DECLARED AND EMITTED.  The first revision of this file pinned the digest
while stating only the ORDERING, which left the number computable by one side alone: the
four-plus plausible renderings of "120 tuples in rank order" (compact JSON, spaced JSON,
trailing newline, flat integer list, component strings) each give a different SHA-256, and
`27ff780d...` names exactly one of them.  A pin only its author can evaluate is an assertion
wearing a hash, the same defect class as naming a file "at main" with no revision.
`ENUMERATION_ENCODING` below therefore states the rule in section 4.2 register, the emitted
artifact carries it beside the digest, and the byte length is recorded so a mismatching third
party can localize rather than guess.

WHAT THIS DOES NOT VERIFY.

  1  That the supplied group is the one the M8 column needs.  As in M8.5-A, this audit
     establishes what the packet IS; the pre-registration answers whether it is the right
     input.
  2  Absence of leakage by construction.  G3 is a scan against a named list, and a scan is a
     floor rather than a proof.  G4 exists because that floor is low.
  3  G5 in isolation.  Every closure change reachable by mutating this packet also reddens
     G2, because a different 120-element unit-quaternion group fails the M8.5-A order census
     and a second embedding of 2I inside this coefficient lattice was not found.  G5 earns
     its place by PINNING the enumeration, which no other check records, not by an isolated
     red.
  4  The binding between this packet and the construction packet.  That the construction
     packet's `group_packet_sha256` reproduces from these bytes is check A3 of
     [`m8_8_packet_audit.py`](m8_8_packet_audit.py), and it stays there: this file runs from
     the group packet alone.

USAGE.
    python3 m8_8_group_packet_reaudit.py                   re-audit, write JSON, exit 0 on green
    python3 m8_8_group_packet_reaudit.py --mutation-tests  additionally run the mutation suite
    python3 m8_8_group_packet_reaudit.py --packet PATH     re-audit a packet elsewhere

Exit code is nonzero if any check fails or any mutation goes undetected, so a red re-audit
cannot be mistaken for a green one.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path
from typing import Any, Iterator

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
DEFAULT_PACKET = DATA / "m8_5a_packet.json"
DEFAULT_OUT = DATA / "m8_8_group_packet_reaudit.json"

sys.path.insert(0, str(HERE))

import m8_5a_packet_audit as m85  # noqa: E402
import m8_8_packet_audit as m88  # noqa: E402

# Declared here so every check compares against a STATED reference rather than against
# whatever the packet happens to produce.
EXPECTED_PACKET_SHA256 = "e3b0c945bbbb15b4549fa641234c9461062c2337b3d1e372af621b614d4883a9"
EXPECTED_ENUMERATION_SHA256 = "27ff780d28d5d854d464ead87e8fc20244fac8334bda9f0600c6ee1b3c89561e"
EXPECTED_ORDER = 120

# The encoding the digest above is taken over, stated so both sides can compute it.  Emitted
# into the artifact verbatim: the rule travels with the number rather than with this file.
ENUMERATION_ENCODING = (
    "CANONICAL TUPLE: the element's four quaternion components in the packet's "
    "`quaternion_basis` order (1, i, j, k), each contributing the integer pair (A, B) of "
    "its `(A + B*phi)/2` numerator, rational part before phi part, the denominator fixed "
    "at 2 and dropped, giving 8 signed integers.  RANK: the 120 tuples sorted ascending in "
    "lexicographic order, compared entrywise as signed integers with the first entry most "
    "significant; rank is the 0-based position and is the element ID of section 4.2.  "
    "DIGEST INPUT: one JSON array of the 120 rank-ordered 8-integer arrays, rendered with "
    "no whitespace at all (JSON separators ',' and ':'), integers as bare decimal with '-' "
    "for negative, no '+', no leading zero, no exponent; ASCII; no trailing newline.  "
    "DIGEST: SHA-256 of those bytes, lowercase hex."
)

# Section 4.3's forbidden construction inputs, reduced to the words their objects are made
# of: the M8.3 artifacts, the spectral side they were computed from, and the torsion-side
# quantities the run must not receive.  Matched as a prefix at a word boundary, the same
# matcher both prior audits use.  `phi` is deliberately absent; see the header.
FORBIDDEN_VOCABULARY = (
    "torsion",
    "reidemeister",
    "singer",
    "whitehead",
    "zeta",
    "spectral",
    "spectrum",
    "eigen",
    "laplac",
    "heat",
    "determinant",
    "character",
    "irrep",
    "irreducible",
    "ratio",
    "sector",
    "sqrt",
    "mass",
    "answer",
    "target",
)

DECIMAL_RE = re.compile(r"\d+\.\d+|\d+[eE][-+]?\d+")


# --------------------------------------------------------------------------------------
# The canonical enumeration, built the cheap way and tied once to the audit's own Group
# --------------------------------------------------------------------------------------


def enumeration(generators: list, cap: int = 5000) -> list:
    """The section 4.2 enumeration: canonical coordinate tuples in lexicographic order.

    `m8_8_packet_audit.Group` produces the same list on its way to a full Cayley table.
    This builds only the closure, so the mutation suite can afford one per mutation; main()
    checks the two against each other once on the delivered packet.

    The cap is the M8.5-A `close()` cap and it is load-bearing rather than defensive: a
    mutated generator of norm 3 generates an INFINITE group, and the first draft of this
    file, which had no cap, hung on its own `norm_not_unit` mutation instead of reporting
    it.  A non-finite generating set is a loud failure here, never a wait.
    """
    seen = {m88.canonical_key(m88.H_ONE)}
    frontier = [m88.H_ONE]
    while frontier:
        nxt = []
        for x in frontier:
            for g in generators:
                y = m88.h_mul(x, g)
                key = m88.canonical_key(y)
                if key not in seen:
                    seen.add(key)
                    nxt.append(y)
                    if len(seen) > cap:
                        raise RuntimeError(f"closure exceeded {cap} elements; not finite")
        frontier = nxt
    return sorted(seen)


def enumeration_blob(keys: list) -> bytes:
    """The exact bytes the digest is taken over, per `ENUMERATION_ENCODING`.

    Separated from the hashing so the byte length reaches the artifact and so this one
    expression, rather than a sentence about it, is what the declared encoding describes.
    On the delivered packet the blob is 2389 bytes and runs
    `[[-2,0,0,0,0,0,0,0],[-1,0,-1,0,-1,0,-1,0],` ... `,[1,0,1,0,1,0,1,0],[2,0,0,0,0,0,0,0]]`.
    """
    return json.dumps([list(k) for k in keys], separators=(",", ":")).encode("ascii")


def enumeration_digest(keys: list) -> str:
    return hashlib.sha256(enumeration_blob(keys)).hexdigest()


def numeric_leaves(node: Any, path: str = "$") -> Iterator[tuple]:
    """Every int, float or bool leaf in the parsed packet, with its path."""
    if isinstance(node, dict):
        for key, value in node.items():
            yield from numeric_leaves(value, f"{path}.{key}")
    elif isinstance(node, list):
        for i, value in enumerate(node):
            yield from numeric_leaves(value, f"{path}[{i}]")
    elif isinstance(node, (int, float)):  # bool included on purpose, it is an int subclass
        yield (path, node)


# --------------------------------------------------------------------------------------
# The checks
# --------------------------------------------------------------------------------------


def run_checks(raw: bytes, packet: dict) -> tuple:
    """Returns (checks, detail).  Each check is (id, description, pass, observed)."""
    checks: list = []
    detail: dict = {}

    def add(cid, desc, passed, observed):
        checks.append({"id": cid, "check": desc, "pass": bool(passed), "observed": observed})

    text = raw.decode("utf-8", errors="replace").lower()

    # G1 identity.  Every mutation below also reddens this one, since it pins the exact
    # bytes; each mutation row names the check it TARGETS, and the suite is a subset test.
    incoming = hashlib.sha256(raw).hexdigest()
    canon = m85.canonical_bytes(packet)
    authoritative = hashlib.sha256(canon).hexdigest()
    g1 = incoming == authoritative == EXPECTED_PACKET_SHA256
    add(
        "G1",
        "delivered bytes reproduce the protocol section 11 pin and are already canonical",
        g1,
        {
            "incoming_sha256": incoming,
            "authoritative_sha256": authoritative,
            "expected_sha256": EXPECTED_PACKET_SHA256,
            "canonicalization_changed_bytes": incoming != authoritative,
        },
    )

    # G2 the M8.5-A structural audit, imported rather than restated
    m85_checks, m85_detail = m85.run_checks(raw, packet)
    m85_reds = sorted(c["id"] for c in m85_checks if not c["pass"])
    add(
        "G2",
        "the M8.5-A structural audit re-runs green on these bytes, A1 through A8",
        not m85_reds,
        {"m8_5a_checks_red": m85_reds, "m8_5a_checks_run": len(m85_checks)},
    )
    detail["m8_5a_detail"] = m85_detail

    # G3 section 4.3 vocabulary, M8.8's list
    hits = sorted({t for t in FORBIDDEN_VOCABULARY if re.search(rf"\b{re.escape(t)}", text)})
    add(
        "G3",
        "no section 4.3 answer-bearing vocabulary, under M8.8's list rather than M8.5-A's",
        not hits,
        {
            "token_hits": hits,
            "tokens_scanned": len(FORBIDDEN_VOCABULARY),
            # recorded so the deliberate absence of `phi` from the list is visible in the
            # artifact rather than only in this file's header
            "phi_occurrences_permitted": len(re.findall(r"\bphi", text)),
        },
    )

    # G4 no evaluated numeric content
    literals = list(numeric_leaves(packet))
    decimals = sorted(set(DECIMAL_RE.findall(text)))
    add(
        "G4",
        "no JSON numeric literal and no decimal or scientific rendering anywhere",
        not literals and not decimals,
        {
            "numeric_literals": [{"path": p, "value": v} for p, v in literals],
            "decimal_renderings": decimals,
        },
    )

    # G5 the canonical enumeration the construction packet's element IDs address
    try:
        gens = [tuple(m88.parse_component(c) for c in g) for g in packet["generators"]]
        keys = enumeration(gens)
        blob = enumeration_blob(keys)
        digest = hashlib.sha256(blob).hexdigest()
        error = None
    except (ValueError, KeyError, TypeError, RuntimeError) as exc:
        keys, blob, digest, error = [], b"", None, f"{type(exc).__name__}: {exc}"
    ordered = keys == sorted(keys)
    g5 = len(keys) == EXPECTED_ORDER and ordered and digest == EXPECTED_ENUMERATION_SHA256
    add(
        "G5",
        "canonical enumeration: 120 elements in lexicographic rank order, digest as pinned",
        g5,
        {
            "size": len(keys),
            "lexicographically_ordered": ordered,
            "enumeration_sha256": digest,
            "expected_sha256": EXPECTED_ENUMERATION_SHA256,
            # the rule and the length travel WITH the number, so the digest is checkable by
            # a party that has only this artifact and the group packet
            "encoding": ENUMERATION_ENCODING,
            "serialized_bytes": len(blob),
            "error": error,
        },
    )
    detail["enumeration_sha256"] = digest
    detail["group_order"] = len(keys)
    return checks, detail


# --------------------------------------------------------------------------------------
# Mutation suite: every check above is shown to go red, and the two decisions to stay green
# --------------------------------------------------------------------------------------


def mutations(packet: dict) -> list:
    """Each entry is (name, mutated packet, the check ids that MUST go red)."""
    out = []

    m = json.loads(json.dumps(packet))
    m["format_version"] = "m8_5a-packet-v2"
    out.append(("format_version_string_changed", m, ["G1"]))

    # The isolation demonstration, and the reason this file exists: M8.5-A's scan does not
    # carry `torsion`, so G2 stays GREEN here while G3 reddens.
    m = json.loads(json.dumps(packet))
    m["coefficient_form"] = packet["coefficient_form"] + "; torsion of the R7 sector"
    out.append(("torsion_word_in_declared_prose", m, ["G3"]))

    # A planted value that no vocabulary scan can see, in either list.
    m = json.loads(json.dumps(packet))
    m["coefficient_form"] = packet["coefficient_form"] + "; 1.6180339887"
    out.append(("decimal_rendering_in_declared_prose", m, ["G4"]))

    # Numeric leaf under a declared key, so the key set stays clean and only G4 sees it.
    m = json.loads(json.dumps(packet))
    m["format_version"] = 1
    out.append(("numeric_literal_under_declared_key", m, ["G4"]))

    # a generator of norm 3 generates an infinite group, so both closures refuse at their cap
    m = json.loads(json.dumps(packet))
    m["generators"][0][1] = "(3 + 0*phi)/2"
    out.append(("norm_not_unit", m, ["G2", "G5"]))

    m = json.loads(json.dumps(packet))
    m["generators"] = [packet["generators"][0], packet["generators"][0]]
    out.append(("generators_collapsed_to_one", m, ["G2", "G5"]))

    m = json.loads(json.dumps(packet))
    m["mckay_distance"] = 3
    out.append(("stray_key_added", m, ["G2"]))

    return out


def invariances(packet: dict) -> list:
    """Each entry is (name, mutated packet, the check ids that MUST STAY green).

    These hold two decisions in place: that `phi` is permitted vocabulary here, and that the
    enumeration does not depend on the order the packet lists its generators in.
    """
    out = []

    m = json.loads(json.dumps(packet))
    # deliberately not "the golden ratio", which would redden G3 on `ratio` and prove the
    # opposite of what this row is for
    m["coefficient_form"] = packet["coefficient_form"] + "; phi generates the field"
    out.append(("extra_phi_occurrence_in_prose", m, ["G3"]))

    m = json.loads(json.dumps(packet))
    m["generators"] = list(reversed(packet["generators"]))
    out.append(("generator_order_reversed", m, ["G5"]))

    return out


def run_mutation_suite(packet: dict) -> list:
    results = []
    for name, mutated, expected_red in mutations(packet):
        raw = m85.canonical_bytes(mutated)
        try:
            checks, _ = run_checks(raw, mutated)
            reds = sorted(c["id"] for c in checks if not c["pass"])
            error = None
        except Exception as exc:  # a mutation that breaks the run counts as red, loudly
            reds, error = ["EXCEPTION"], f"{type(exc).__name__}: {exc}"
        ok = all(cid in reds or "EXCEPTION" in reds for cid in expected_red)
        results.append(
            {
                "mutation": name,
                "expected_red": expected_red,
                "observed_red": reds,
                "detected": ok,
                "error": error,
            }
        )
    return results


def run_invariance_suite(packet: dict) -> list:
    results = []
    for name, mutated, expected_green in invariances(packet):
        raw = m85.canonical_bytes(mutated)
        try:
            checks, _ = run_checks(raw, mutated)
            greens = sorted(c["id"] for c in checks if c["pass"])
            error = None
        except Exception as exc:
            greens, error = [], f"{type(exc).__name__}: {exc}"
        ok = all(cid in greens for cid in expected_green)
        results.append(
            {
                "invariance": name,
                "expected_green": expected_green,
                "observed_green": greens,
                "held": ok,
                "error": error,
            }
        )
    return results


# --------------------------------------------------------------------------------------


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--packet", type=Path, default=DEFAULT_PACKET)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--mutation-tests", action="store_true")
    args = ap.parse_args()

    raw = args.packet.read_bytes()
    packet = json.loads(raw)

    checks, detail = run_checks(raw, packet)
    all_pass = all(c["pass"] for c in checks)

    # Tie the cheap closure used above to the enumeration the construction audit actually
    # consumes, once, on the delivered packet.  If these ever disagree, G5 pins nothing.
    gens = [tuple(m88.parse_component(c) for c in g) for g in packet["generators"]]
    group = m88.Group(gens)
    via_group = [m88.canonical_key(e) for e in group.element]
    agrees = enumeration_digest(via_group) == detail.get("enumeration_sha256")
    detail["agrees_with_construction_audit_group_class"] = agrees
    detail["identity_rank"] = group.identity

    mut = run_mutation_suite(packet) if args.mutation_tests else None
    mut_ok = all(m["detected"] for m in mut) if mut else None
    inv = run_invariance_suite(packet) if args.mutation_tests else None
    inv_ok = all(i["held"] for i in inv) if inv else None
    m85_mut = m85.run_mutation_suite(packet) if args.mutation_tests else None
    m85_mut_ok = all(m["detected"] for m in m85_mut) if m85_mut else None

    result = {
        "what": "M8.8 group-packet re-audit, run from the group packet alone",
        "why": (
            "protocol section 4.3 requires the M8.5-A packet re-audited against M8.8's "
            "forbidden list rather than inherited on its M8.5-A audit; section 8 step 1 "
            "carries this record in the content commit"
        ),
        "packet_file": args.packet.name,
        "inherited_record_not_relied_on": "m8_5a_packet_audit.json",
        "forbidden_vocabulary": list(FORBIDDEN_VOCABULARY),
        "vocabulary_note": (
            "`phi` is absent from this list by decision: it is forbidden in the M8.8 "
            "CONSTRUCTION packet, which carries integers and element IDs only, and is the "
            "coefficient-field generator here.  Held by the extra_phi_occurrence_in_prose "
            "invariance rather than by a text subtraction."
        ),
        "checks": checks,
        "detail": detail,
        "all_checks_pass": all_pass,
        "enumeration_agrees_with_construction_audit": agrees,
        "mutation_suite": mut,
        "mutation_suite_all_detected": mut_ok,
        "invariance_suite": inv,
        "invariance_suite_all_held": inv_ok,
        "m8_5a_mutation_suite_all_detected": m85_mut_ok,
    }
    args.out.write_text(json.dumps(result, indent=2) + "\n")

    print(f"packet          {args.packet}")
    print(f"pin             {EXPECTED_PACKET_SHA256}")
    print(f"enumeration     {detail.get('enumeration_sha256')}")
    print()
    for c in checks:
        print(f"  {'PASS' if c['pass'] else 'FAIL'}  {c['id']}  {c['check']}")
        # the declared encoding is a paragraph and belongs in the artifact, not in a
        # terminal line; everything else in `observed` prints as it stands
        terse = {k: v for k, v in c["observed"].items() if k != "encoding"}
        print(f"        {json.dumps(terse)}")
    print()
    print(f"  enumeration agrees with the construction audit's Group class: {agrees}")
    if mut is not None:
        print()
        for m in mut:
            print(
                f"  {'DETECTED' if m['detected'] else 'MISSED  '}  {m['mutation']:36s} "
                f"expected red {m['expected_red']}, observed {m['observed_red']}"
            )
        for i in inv:
            print(
                f"  {'HELD    ' if i['held'] else 'BROKEN  '}  {i['invariance']:36s} "
                f"expected green {i['expected_green']}"
            )
        print(f"  M8.5-A mutation suite all detected: {m85_mut_ok}")
    print()
    print(f"ALL CHECKS PASS: {all_pass}", end="")
    if mut is not None:
        print(f" | MUTATIONS ALL DETECTED: {mut_ok} | INVARIANCES ALL HELD: {inv_ok}", end="")
    print()
    print(f"written: {args.out}")

    green = all_pass and agrees
    return 0 if green and mut_ok is not False and inv_ok is not False else 1


if __name__ == "__main__":
    sys.exit(main())
