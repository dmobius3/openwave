"""P2 mutation battery: every declared P2 check must be shown able to reject,
BY THE BRANCH UNDER TEST.

Target-free in the § 2 sense: not answer-bearing. Mutations are structural, and
no value is transcribed here FROM the record. Stated exactly, because two
earlier versions of this sentence were false. It once claimed every substitute
value was read out of the packet, which the literals below contradict. It then
claimed no target value appears here, while listing [1,0,1] among them, and
[1,0,1] is verbatim the committed value of two rows. The truth is narrower and
worth saying plainly: the triple literals are [1,1,1], [0,0,1], [1,0,1] and
[1,0,0], each chosen for a structural property (a wrong signature component, a
value with no inverse, the multiplicative unit, and a zero denominator), and the
unit coincides with a committed value because the multiplicative unit is what
that row's value IS. Nothing here was copied from the record, and one literal
nonetheless equals a target. This is the FOURTH wording of this sentence; the
first three were each false in a different way, twice by overclaiming and once by
enumerating incompletely, which is why it now says what it says.

The rules P1 ended at, carried forward whole:

  1. A case runs `p2.assess`, the same predicate the driver runs.
  2. EVERY named target must reject, not any one of them.
  3. A case declares the REASON its target's detail must contain. Verifying only
     THAT a check rejected let a fix in one place silently hollow out coverage
     in another, and nothing could see it.
  4. `require_green` names checks that must stay green, so a case cannot pass by
     tripping an earlier gate instead of the one it claims to test.
  5. Coverage is computed from checks OBSERVED rejecting, never from declared
     targets.
"""

import contextlib
import copy
import json
import pathlib
import shutil
import sys
import tempfile

import p2
import p2_emit
import p2_schema as S

CASES = []


def case(*targets, reason=None, require_green=(), allow_blocked=False):
    def deco(fn):
        CASES.append((fn.__name__, tuple(targets), fn, reason,
                      tuple(require_green), allow_blocked))
        return fn
    return deco


def failed_with_details(**kw):
    try:
        checks, stops, _ = p2.assess(**kw)
    except Exception as e:                                    # noqa: BLE001
        return {"__raised__": f"{type(e).__name__}: {e}"}, {}
    # EVERY case asserts a COMPLETE verdict set, not only its own targets.
    # Nothing did: the fourth audit's "a failed pin blocks EVERYTHING" repair
    # could be reverted with the battery still fully green, because no case
    # cared how many verdicts came back. The contract says 38, so 38 is checked
    # on every single run rather than in one case about pins.
    emitted = [c for c, _, _ in checks]
    if sorted(emitted) != sorted(p2.CHECK_IDS):
        return {"__raised__": f"INCOMPLETE REGISTRY: {len(emitted)} of "
                              f"{len(p2.CHECK_IDS)} verdicts emitted"}, {}
    bad = {c for c, ok, _ in checks if not ok}
    details = {c: d for c, ok, d in checks if not ok}
    for s in stops:
        k = s.split(":", 1)[0]
        if k in EXTRA_TARGETS:
            bad.add(k)
            details[k] = s
    return bad, details


EXTRA_TARGETS = {"REGISTRY"}


def case_passed(targets, bad, raised, reason=None, details=None,
                require_green=(), allow_blocked=False):
    """A BLOCKED target does not count as rejection by the branch under test.

    One case deleted a required nested field, which made the nested-schema pass
    return early and BLOCK the check it named. The target appeared among the
    failures, so the case passed, while the branch it claimed to exercise never
    ran. Same shape as the P1 lesson one level up: rejection is not coverage
    unless it is rejection by the branch under test."""
    if raised or not set(targets) <= set(bad):
        return False
    if any(g in bad for g in require_green):
        return False
    if not allow_blocked:
        d = details or {}
        if any(str(d.get(x, "")).startswith("BLOCKED:") for x in targets):
            return False
    if reason is not None:
        d = details or {}
        want = reason if isinstance(reason, dict) else {t: reason for t in targets}
        return all(sub in str(d.get(t, "")) for t, sub in want.items())
    return True


# ------------------------------------------------------------------ input

@case("P1_ARTIFACT_PINNED", reason="not the pinned")
def p1_artifact_not_pinned(base, tmp):
    """P2 consumes exactly one P1 version. Anything else is refused."""
    src = pathlib.Path(base["p1_artifact"])
    dst = tmp / "p1.json"
    dst.write_bytes(src.read_bytes() + b"\n")
    return {"p1_artifact": str(dst)}


# ------------------------------------------------------------ top-level

def _mutate(base, fn):
    """Build the real packet, apply a structural mutation, hand it to assess."""
    p1 = p2_emit.load_p1(base["p1_artifact"])
    pkt = copy.deepcopy(p2_emit.build_packet(p1))
    fn(pkt)          # mutates in place; the return value is DISCARDED, because
                     # `.pop()` returns the popped object and merging that into
                     # the kwargs turned three cases into TypeErrors instead of
                     # rejections
    return {"p1_artifact": base["p1_artifact"], "packet": pkt}


@case("KEYS_EXACTLY_SEVEN", reason="want", require_green=("P1_ARTIFACT_PINNED",))
def extra_top_level_key(base, tmp):
    return _mutate(base, lambda p: p.update({"notes": "extra"}))


@case("FORMAT_VERSION", require_green=("KEYS_EXACTLY_SEVEN",))
def wrong_format_version(base, tmp):
    return _mutate(base, lambda p: p.update({"format_version": "m8_8-answer-packet-0"}))


@case("TARGET_ID", require_green=("KEYS_EXACTLY_SEVEN",))
def wrong_target_id(base, tmp):
    return _mutate(base, lambda p: p.update({"target_id": "M88-ADJ-99"}))


@case("ADJUDICATES_MATCH_P1", require_green=("KEYS_EXACTLY_SEVEN",))
def wrong_adjudicates(base, tmp):
    def f(p):
        h = p["adjudicates"]["group_packet_sha256"]
        p["adjudicates"]["group_packet_sha256"] = h[:-1] + ("0" if h[-1] != "0" else "1")
    return _mutate(base, f)


# ------------------------------------------------------------------ rows

@case("ROWS_COUNT", require_green=("ROWS_IN_CANONICAL_ORDER",))
def row_dropped(base, tmp):
    return _mutate(base, lambda p: p["rows"].pop())


@case("ROWS_IN_CANONICAL_ORDER", reason="differs from the frozen sort key",
      require_green=("ROWS_COUNT", "ROWS_SIGNATURES_DISTINCT"))
def rows_out_of_order(base, tmp):
    """Swap two adjacent rows. Emitting the array already ordered is what makes
    the order verifiable rather than inferred."""
    def f(p):
        p["rows"][0], p["rows"][1] = p["rows"][1], p["rows"][0]
    return _mutate(base, f)


@case("ROWS_SIGNATURES_DISTINCT", require_green=("ROWS_COUNT",))
def duplicate_signature(base, tmp):
    def f(p):
        p["rows"][1]["signature"] = copy.deepcopy(p["rows"][0]["signature"])
    return _mutate(base, f)


@case("ROWS_TRIPLES_NORMALIZED", reason="not normalized",
      require_green=("ROWS_NO_FLOATS",))
def unnormalized_triple(base, tmp):
    """Scale a value by 2 in both numerator and denominator: same number, no
    longer the frozen encoding."""
    def f(p):
        a, b, c = p["rows"][0]["value"]
        p["rows"][0]["value"] = [a * 2, b * 2, c * 2]
    return _mutate(base, f)


@case("ROWS_NO_FLOATS", reason="floats at")
def decimal_literal(base, tmp):
    """§ 4.4: never a decimal. A float anywhere in the packet fails."""
    def f(p):
        p["rows"][0]["value"] = [float(x) for x in p["rows"][0]["value"]]
    return _mutate(base, f)


@case("ROWS_VALUES_MATCH_P1", require_green=("ROWS_TRIPLES_NORMALIZED",))
def value_swapped_between_rows(base, tmp):
    """Give one row another row's value. Both stay normalized, both stay
    integers; only the tie back to P1 catches it."""
    def f(p):
        p["rows"][0]["value"], p["rows"][8]["value"] = \
            p["rows"][8]["value"], p["rows"][0]["value"]
    return _mutate(base, f)


@case("ROWS_SIGNATURES_MATCH_P1", require_green=("ROWS_SIGNATURES_DISTINCT",))
def signature_component_altered(base, tmp):
    def f(p):
        p["rows"][3]["signature"]["t"] = [1, 1, 1]
    return _mutate(base, f)


@case("CLASS_CENSUS", require_green=("ROWS_COUNT",))
def wrong_class_census(base, tmp):
    def f(p):
        for r in p["rows"]:
            if r["evidentiary_class"] == S.CLASS_FREE:
                r["evidentiary_class"] = S.CLASS_DECLARED
                break
    return _mutate(base, f)


@case("CLASS_SELECTOR_NOT_SELF_INVERSE", reason="SELF-INVERSE",
      require_green=("CLASS_CENSUS",))
def selector_self_inverse(base, tmp):
    """A selector whose value is its own inverse cannot select an orientation:
    the § 8 invalid-anchor branch."""
    def f(p):
        for r in p["rows"]:
            if r["evidentiary_class"] == S.CLASS_SELECTOR:
                r["value"] = [1, 0, 1]
    return _mutate(base, f)


# ------------------------------------------------------------ identities

@case("IDY_COUNT", require_green=("IDY_FACTOR_ORDER_CANONICAL",))
def identity_dropped(base, tmp):
    return _mutate(base, lambda p: p["identities"].pop())


@case("IDY_POSITIONS_CONTIGUOUS", require_green=("IDY_COUNT",))
def identity_position_gap(base, tmp):
    def f(p):
        p["identities"][2]["position"] = 7
    return _mutate(base, f)


@case("IDY_FACTORS_RESOLVE_BY_SIGNATURE", reason="do not resolve uniquely",
      require_green=("ROWS_SIGNATURES_DISTINCT",))
def factor_signature_matches_no_row(base, tmp):
    def f(p):
        p["identities"][0]["factors"][0]["signature"]["dimension"] = 99
    return _mutate(base, f)


@case("IDY_RECOMPUTES_FROM_ROWS", reason="recompute, want",
      require_green=("IDY_FACTORS_RESOLVE_BY_SIGNATURE",))
def identity_expected_value_drifts(base, tmp):
    """Point one slot's expected value at another slot's. The identity no longer
    follows from the rows the packet itself commits."""
    def f(p):
        p["identities"][0]["expected_value"] = list(p["identities"][2]["expected_value"])
    return _mutate(base, f)


# ---------------------------------------------------------- indexing map

@case("NESTED_SCHEMA_EXACT", reason="missing field",
      require_green=("TOP_LEVEL_TYPES", "KEYS_EXACTLY_SEVEN"))
def indexing_domain_missing(base, tmp):
    """Retargeted. This deletes a required nested field, so the nested-schema
    pass is the branch that rejects it; the domain check is BLOCKED and never
    runs. It previously claimed IDXMAP_DOMAINS_DECLARED as its target and passed
    on the block. `indexing_arbitrary_domain_name` exercises that check
    directly."""
    return _mutate(base, lambda p: p["indexing_map"].pop("destination_domain"))


@case("IDXMAP_ZERO_BASED_CONTIGUOUS", reason="not zero-based contiguous",
      require_green=("IDXMAP_BIJECTIVE",))
def indexing_one_based(base, tmp):
    def f(p):
        for e in p["indexing_map"]["entries"]:
            e["source_position"] += 1
            e["destination_position"] += 1
    return _mutate(base, f)


@case("IDXMAP_BIJECTIVE", reason="not a bijection",
      require_green=("IDXMAP_DOMAINS_DECLARED",))
def indexing_not_bijective(base, tmp):
    def f(p):
        p["indexing_map"]["entries"][1]["destination_position"] = \
            p["indexing_map"]["entries"][0]["destination_position"]
    return _mutate(base, f)


@case("IDXMAP_DISCRIMINATING", reason="cannot discriminate",
      require_green=("IDXMAP_BIJECTIVE",))
def indexing_is_identity(base, tmp):
    """An identity map cannot distinguish a harness that applies it from one
    that ignores it, which is the whole point of the § 4.4 control."""
    def f(p):
        for e in p["indexing_map"]["entries"]:
            e["destination_position"] = e["source_position"]
    return _mutate(base, f)


# --------------------------------------------------------- convention map

@case("CONV_MAP_COMPLETE")
def convention_section_missing(base, tmp):
    return _mutate(base, lambda p: p["convention_map"].pop("bridge"))


@case("CONV_BASING_MATCHES_ADJUDICATES", reason="pins a DIFFERENT construction packet",
      require_green=("CONV_MAP_COMPLETE",))
def basing_pins_a_different_packet(base, tmp):
    def f(p):
        h = p["convention_map"]["basing_reference"]["construction_packet_sha256"]
        p["convention_map"]["basing_reference"]["construction_packet_sha256"] = \
            h[:-1] + ("0" if h[-1] != "0" else "1")
    return _mutate(base, f)


# ----------------------------------------------------------------- labels

@case("LABELS_ARE_NOT_JOIN_KEYS", "ROWS_SIGNATURES_DISTINCT",
      reason={"LABELS_ARE_NOT_JOIN_KEYS": "cannot be resolved without labels",
              "ROWS_SIGNATURES_DISTINCT": "distinct of"})
def label_would_be_needed_to_resolve(base, tmp):
    """Two rows sharing a signature: a factor can then only be resolved by label.

    HONEST DEPENDENCY, recorded rather than hidden. Once the nested schema is
    exact, a stray label anywhere outside `rows` is an undeclared field and is
    rejected earlier, so this check can no longer redden alone. Signatures being
    the only join key is now a consequence of the frozen schema rather than an
    independent property. It is kept because a future schema edit could make it
    independent again, and its case names both checks that necessarily fire."""
    def f(p):
        p["rows"][1]["signature"] = copy.deepcopy(p["rows"][0]["signature"])
    return _mutate(base, f)


# --------------------------------------------------------- bytes and record

@case("CANONICAL_BYTES_FIXED_POINT", reason="canonical serialization",
      require_green=("KEYS_EXACTLY_SEVEN",))
def non_canonical_bytes(base, tmp):
    """Four-space indent instead of two: parses identically, different bytes."""
    p1 = p2_emit.load_p1(base["p1_artifact"])
    pkt = p2_emit.build_packet(p1)
    bad = (json.dumps(pkt, sort_keys=True, indent=4, ensure_ascii=True) + "\n").encode()
    return {"p1_artifact": base["p1_artifact"], "packet": pkt, "packet_bytes": bad,
            "hash_record": p2_emit.build_hash_record(bad)}


@case("HASH_RECORD_MATCHES_PACKET", reason="sha256",
      require_green=("CANONICAL_BYTES_FIXED_POINT",))
def hash_record_wrong_hash(base, tmp):
    p1 = p2_emit.load_p1(base["p1_artifact"])
    pkt = p2_emit.build_packet(p1)
    rec = p2_emit.build_hash_record(S.canonical_bytes(pkt))
    h = rec["expected_canonical_plaintext_sha256"]
    rec["expected_canonical_plaintext_sha256"] = h[:-1] + ("0" if h[-1] != "0" else "1")
    return {"p1_artifact": base["p1_artifact"], "packet": pkt, "hash_record": rec}


@case("HASH_RECORD_IS_CANDIDATE", require_green=("HASH_RECORD_MATCHES_PACKET",))
def hash_record_claims_issued(base, tmp):
    """Immutability attaches at P6. A P2 record claiming ISSUED would mean a
    later gate failure forced a reissue of something never issued."""
    p1 = p2_emit.load_p1(base["p1_artifact"])
    pkt = p2_emit.build_packet(p1)
    rec = p2_emit.build_hash_record(S.canonical_bytes(pkt))
    rec["status"] = "ISSUED"
    return {"p1_artifact": base["p1_artifact"], "packet": pkt, "hash_record": rec}


@case("REGISTRY", reason="REGISTRY")
def registry_check_dropped(base, tmp):
    orig = p2.validate_registry
    p2.validate_registry = lambda emitted: orig(list(emitted)[:-1])
    try:
        return {"p1_artifact": base["p1_artifact"], "__restore__": orig}
    finally:
        pass


# ------------------------------- the audit's five internally-consistent fakes
# Each of these produced a packet that passed every check while being a
# materially different object. They are the reason P2 now checks fidelity to P1
# and to a frozen template, not merely self-consistency.

@case("NESTED_SCHEMA_EXACT", reason="undeclared field",
      require_green=("KEYS_EXACTLY_SEVEN",))
def undeclared_nested_field(base, tmp):
    """An extra field on a row. Freezing only the top-level keys let this pass."""
    def f(p):
        p["rows"][0]["note"] = "undeclared"
    return _mutate(base, f)


@case("NESTED_SCHEMA_EXACT", reason="missing field",
      require_green=("KEYS_EXACTLY_SEVEN",))
def missing_nested_field(base, tmp):
    """A missing nested field used to make a later stage RAISE rather than
    return a verdict."""
    def f(p):
        del p["rows"][0]["signature"]["t"]
    return _mutate(base, f)


@case("IDY_DEFINITIONS_MATCH_P1", reason="differs from P1",
      require_green=("IDY_FACTORS_RESOLVE_BY_SIGNATURE", "IDY_RECOMPUTES_FROM_ROWS",
                     "IDY_COUNT", "IDY_POSITIONS_CONTIGUOUS"))
def identity_redefined_consistently(base, tmp):
    """THE decisive one. Point a factor at a different resolving row signature
    and update the slot's expected value to the recomputed result. The packet
    stays internally consistent and is no longer P1's identity."""
    from qphi_exact import Phi

    def f(p):
        slot = p["identities"][0]
        used = {S.signature_of(x) for x in slot["factors"]}
        donor = next(r for r in p["rows"] if S.signature_of(r) not in used)
        slot["factors"][0]["signature"] = copy.deepcopy(donor["signature"])
        slot["factors"].sort(key=S.factor_sort_key)
        by_sig = {S.signature_of(r): r for r in p["rows"]}
        acc = Phi(1)
        for fac in slot["factors"]:
            a, b, c = by_sig[S.signature_of(fac)]["value"]
            acc = acc * ((Phi(a, b) / Phi(c)) ** fac["exponent"])
        slot["expected_value"] = list(acc.triple())
    return _mutate(base, f)


@case("IDY_FACTOR_ORDER_CANONICAL", reason="unordered",
      require_green=("IDY_DEFINITIONS_MATCH_P1", "IDY_RECOMPUTES_FROM_ROWS"))
def factors_out_of_frozen_order(base, tmp):
    """Equivalent product, different serialization. The bytes get hashed."""
    def f(p):
        p["identities"][2]["factors"].reverse()
    return _mutate(base, f)


@case("IDXMAP_MATCHES_P1", reason="differs from the P1-derived",
      require_green=("IDXMAP_BIJECTIVE", "IDXMAP_DISCRIMINATING",
                     "IDXMAP_ZERO_BASED_CONTIGUOUS"))
def indexing_valid_but_wrong(base, tmp):
    """Exchange two destinations. Still contiguous, still bijective, still
    non-identity, and no longer the correspondence anyone declared."""
    def f(p):
        e = p["indexing_map"]["entries"]
        e[0]["destination_position"], e[3]["destination_position"] = \
            e[3]["destination_position"], e[0]["destination_position"]
    return _mutate(base, f)


@case("IDXMAP_DOMAINS_DECLARED", reason="destination_domain",
      require_green=("IDXMAP_BIJECTIVE",))
def indexing_arbitrary_domain_name(base, tmp):
    """A plausible but undeclared domain string."""
    def f(p):
        p["indexing_map"]["destination_domain"] = "some_other_ordering"
    return _mutate(base, f)


@case("CONV_MAP_EXACT", reason="differs from the frozen template",
      require_green=("CONV_MAP_COMPLETE", "CONV_BASING_MATCHES_ADJUDICATES"))
def bridge_identification_wrong(base, tmp):
    """Replace the bridge identification. All four subsections still present."""
    def f(p):
        p["convention_map"]["bridge"]["identification"] = "WRONG"
    return _mutate(base, f)


@case("CONV_MAP_EXACT", reason="differs from the frozen template",
      require_green=("CONV_MAP_COMPLETE",))
def anchor_rule_wrong(base, tmp):
    def f(p):
        p["convention_map"]["orientation_anchor_rule"]["rule"] = \
            "selected by the implementation"
    return _mutate(base, f)


@case("CONV_MAP_EXACT", reason="differs from the frozen template",
      require_green=("CONV_MAP_COMPLETE",))
def value_encoding_wrong(base, tmp):
    def f(p):
        p["convention_map"]["value_encoding"]["floating_point"] = \
            "accepted to 12 significant figures"
    return _mutate(base, f)


@case("HASH_RECORD_SCHEMA_EXACT", reason="missing field",
      require_green=("HASH_RECORD_MATCHES_PACKET", "HASH_RECORD_IS_CANDIDATE"))
def hash_record_missing_core_fields(base, tmp):
    """Drop the canonicalization rule and the packet format version. The record
    still hashed the right bytes, which was all the checker required."""
    p1 = p2_emit.load_p1(base["p1_artifact"])
    pkt = p2_emit.build_packet(p1)
    rec = p2_emit.build_hash_record(S.canonical_bytes(pkt))
    rec.pop("canonicalization_rule")
    rec.pop("answer_packet_format_version")
    return {"p1_artifact": base["p1_artifact"], "packet": pkt, "hash_record": rec}


# ------------------------- the second audit: five more internally-consistent
# fakes, each of which passed all 34 checks before these were added.

@case("CLASS_ASSIGNMENT_MATCHES_P1", reason="wrong class",
      require_green=("CLASS_CENSUS", "ROWS_COUNT"))
def declared_class_moved_census_preserved(base, tmp):
    """Swap the declared class onto a free row. 1/7/1 is preserved, so the
    census check is satisfied and the packet is a different object."""
    def f(p):
        a = next(r for r in p["rows"] if r["evidentiary_class"] == S.CLASS_DECLARED)
        b = next(r for r in p["rows"] if r["evidentiary_class"] == S.CLASS_FREE)
        a["evidentiary_class"], b["evidentiary_class"] = \
            b["evidentiary_class"], a["evidentiary_class"]
    return _mutate(base, f)


@case("CLASS_ASSIGNMENT_MATCHES_P1", "CONV_MAP_EXACT",
      reason={"CLASS_ASSIGNMENT_MATCHES_P1": "wrong class",
              "CONV_MAP_EXACT": "differs from the frozen template"},
      require_green=("CLASS_CENSUS", "CLASS_SELECTOR_NOT_SELF_INVERSE"))
def selector_reassigned_with_convention_following(base, tmp):
    """Move the selector class to another row AND update the convention map's
    anchor signature to follow it. The map used to derive its expectation from
    whichever row claimed the class, so the two agreed with each other and with
    nothing else."""
    def f(p):
        a = next(r for r in p["rows"] if r["evidentiary_class"] == S.CLASS_SELECTOR)
        b = next(r for r in p["rows"]
                 if r["evidentiary_class"] == S.CLASS_FREE
                 and as_phi_local(r["value"]) != (1, 0, 1))
        a["evidentiary_class"], b["evidentiary_class"] = \
            b["evidentiary_class"], a["evidentiary_class"]
        p["convention_map"]["orientation_anchor_rule"]["anchor_row_signature"] = \
            copy.deepcopy(b["signature"])
    return _mutate(base, f)


def as_phi_local(triple):
    from qphi_exact import Phi
    a, b, c = triple
    v = Phi(a, b) / Phi(c)
    return (v * v).triple()


@case("BYTES_BIND_TO_PACKET", reason="not exactly the canonical serialization",
      require_green=("CANONICAL_BYTES_FIXED_POINT", "HASH_RECORD_MATCHES_PACKET",
                     "TARGET_ID"))
def bytes_are_a_different_packet(base, tmp):
    """A valid packet object, canonical bytes for a DIFFERENT packet, and a
    sidecar matching those bytes. Semantic checks read the object, byte checks
    read the bytes, and nothing tied them together."""
    p1 = p2_emit.load_p1(base["p1_artifact"])
    good = p2_emit.build_packet(p1)
    other = copy.deepcopy(good)
    other["target_id"] = "M88-ADJ-99"
    ob = S.canonical_bytes(other)
    return {"p1_artifact": base["p1_artifact"], "packet": good,
            "packet_bytes": ob, "hash_record": p2_emit.build_hash_record(ob)}


@case("HASH_RECORD_EXACT", reason="differs from the frozen template",
      require_green=("HASH_RECORD_SCHEMA_EXACT", "HASH_RECORD_MATCHES_PACKET",
                     "HASH_RECORD_IS_CANDIDATE"))
def sidecar_canonicalization_rule_falsified(base, tmp):
    """A false canonicalization rule. The record still hashed the right bytes and
    carried every declared field, which was all the checker required."""
    p1 = p2_emit.load_p1(base["p1_artifact"])
    pkt = p2_emit.build_packet(p1)
    rec = p2_emit.build_hash_record(S.canonical_bytes(pkt))
    rec["canonicalization_rule"] = "JSON, four-space indent, UTF-8, CRLF"
    return {"p1_artifact": base["p1_artifact"], "packet": pkt, "hash_record": rec}


@case("HASH_RECORD_EXACT", reason="differs from the frozen template",
      require_green=("HASH_RECORD_SCHEMA_EXACT", "HASH_RECORD_IS_CANDIDATE"))
def sidecar_banner_weakened(base, tmp):
    """A weakened quarantine banner is the one that actually costs something."""
    p1 = p2_emit.load_p1(base["p1_artifact"])
    pkt = p2_emit.build_packet(p1)
    rec = p2_emit.build_hash_record(S.canonical_bytes(pkt))
    rec["_banner"] = "Internal draft. Share as needed.\n"
    return {"p1_artifact": base["p1_artifact"], "packet": pkt, "hash_record": rec}


@case("TOP_LEVEL_TYPES", reason="want list", require_green=("KEYS_EXACTLY_SEVEN",))
def rows_is_not_a_list(base, tmp):
    """`rows: None` used to raise before any verdict existed, because the nested
    pass iterated it before its container type was established."""
    def f(p):
        p["rows"] = None
    return _mutate(base, f)


@case("TOP_LEVEL_TYPES", reason="want dict", require_green=("KEYS_EXACTLY_SEVEN",))
def convention_map_is_not_a_dict(base, tmp):
    def f(p):
        p["convention_map"] = "not an object"
    return _mutate(base, f)


@case("TOP_LEVEL_TYPES", reason="hash_record", require_green=("KEYS_EXACTLY_SEVEN",))
def hash_record_is_not_a_dict(base, tmp):
    """A list where the sidecar belongs used to raise after the shape check."""
    p1 = p2_emit.load_p1(base["p1_artifact"])
    pkt = p2_emit.build_packet(p1)
    return {"p1_artifact": base["p1_artifact"], "packet": pkt, "hash_record": []}


# ------------------------- the third audit: model-level holes found by reading
# the predicate rather than by running it.

@case("BYTES_BIND_TO_PACKET", reason="not exactly the canonical serialization",
      require_green=("CANONICAL_BYTES_FIXED_POINT", "HASH_RECORD_MATCHES_PACKET",
                     "ROWS_NO_FLOATS"))
def bytes_use_float_lexeme_for_integer(base, tmp):
    """Python says 1 == 1.0; the protocol does not.

    Object carries the integer, bytes carry the float lexeme, sidecar hashes
    those bytes. The semantic float check reads the OBJECT and sees none, the
    alternate bytes are themselves canonical, and decoded-object equality was
    True. Only byte-exact binding catches it."""
    p1 = p2_emit.load_p1(base["p1_artifact"])
    good = p2_emit.build_packet(p1)
    other = copy.deepcopy(good)
    other["rows"][0]["value"][0] = float(other["rows"][0]["value"][0])
    ob = S.canonical_bytes(other)
    return {"p1_artifact": base["p1_artifact"], "packet": good,
            "packet_bytes": ob, "hash_record": p2_emit.build_hash_record(ob)}


@case("NESTED_SCHEMA_EXACT", reason="not an object",
      require_green=("TOP_LEVEL_TYPES", "KEYS_EXACTLY_SEVEN"))
def row_entry_is_not_an_object(base, tmp):
    """A list where a row belongs. exact_shape reported it correctly and the
    next line then called .get() on it and raised."""
    return _mutate(base, lambda p: p["rows"].__setitem__(0, []))


@case("NESTED_SCHEMA_EXACT", reason="not an object",
      require_green=("TOP_LEVEL_TYPES", "KEYS_EXACTLY_SEVEN"))
def factor_entry_is_not_an_object(base, tmp):
    def f(p):
        p["identities"][0]["factors"][0] = []
    return _mutate(base, f)


@case("TOP_LEVEL_TYPES", reason="packet")
def packet_is_not_an_object(base, tmp):
    """`sorted()` succeeds on a list, so even the key check could not run."""
    return {"p1_artifact": base["p1_artifact"], "packet": []}


# ------------------------ the fourth audit: leaf types, arithmetic on invalid
# input, and an expectation still drawn from the object under test.

def _reparsed(base, fn):
    """Mutate through a JSON round trip, so leaf types are what a real consumer
    would see rather than what Python happens to hold."""
    p1 = p2_emit.load_p1(base["p1_artifact"])
    pkt = json.loads(S.canonical_bytes(p2_emit.build_packet(p1)).decode())
    fn(pkt)
    b = S.canonical_bytes(pkt)
    return {"p1_artifact": base["p1_artifact"], "packet": pkt, "packet_bytes": b,
            "hash_record": p2_emit.build_hash_record(b)}


@case("NESTED_SCHEMA_EXACT", reason="three integer members",
      require_green=("TOP_LEVEL_TYPES", "KEYS_EXACTLY_SEVEN", "ROWS_NO_FLOATS"))
def boolean_leaf_in_expected_value(base, tmp):
    """JSON `true` where an integer belongs. Python says True == 1, so every
    arithmetic comparison downstream accepted it and the container type check
    saw a list, as declared."""
    def f(p):
        p["identities"][0]["expected_value"][2] = True
    return _reparsed(base, f)


@case("NESTED_SCHEMA_EXACT", reason="three integer members",
      require_green=("TOP_LEVEL_TYPES", "KEYS_EXACTLY_SEVEN"))
def boolean_leaf_in_factor_signature(base, tmp):
    def f(p):
        p["identities"][0]["factors"][0]["signature"]["s"][1] = True
    return _reparsed(base, f)


@case("NESTED_SCHEMA_EXACT", reason="three integer members",
      require_green=("TOP_LEVEL_TYPES", "KEYS_EXACTLY_SEVEN"))
def short_signature_triple(base, tmp):
    """An empty list where a triple belongs used to raise IndexError inside the
    canonical sort key, before any verdict existed."""
    def f(p):
        p["rows"][0]["signature"]["s"] = []
    return _reparsed(base, f)


@case("ROWS_TRIPLES_NORMALIZED", "CLASS_SELECTOR_NOT_SELF_INVERSE",
      reason={"ROWS_TRIPLES_NORMALIZED": "not normalized",
              "CLASS_SELECTOR_NOT_SELF_INVERSE": "BLOCKED"},
      allow_blocked=True,
      require_green=("NESTED_SCHEMA_EXACT", "ROWS_COUNT"))
def zero_denominator_triple(base, tmp):
    """A zero denominator is shape-valid and arithmetically fatal. The exact
    inverse raised ZeroDivisionError instead of the run returning verdicts.
    Blocking the dependent check is the CORRECT outcome here, so this case
    declares that explicitly rather than counting a block as a rejection."""
    def f(p):
        sel = next(r for r in p["rows"]
                   if r["evidentiary_class"] == S.CLASS_SELECTOR)
        sel["value"] = [1, 0, 0]
    return _reparsed(base, f)


@case("ROWS_SIGNATURES_MATCH_P1", "IDY_DEFINITIONS_MATCH_P1",
      reason={"ROWS_SIGNATURES_MATCH_P1": "R", "IDY_DEFINITIONS_MATCH_P1": "differs from P1"},
      require_green=("ROWS_SIGNATURES_DISTINCT", "IDY_FACTORS_RESOLVE_BY_SIGNATURE",
                     "IDY_RECOMPUTES_FROM_ROWS"))
def row_signature_and_its_factors_moved_together(base, tmp):
    """Change a row's signature AND every factor referencing it, coherently.

    The identity check used to build its P1 expectation from the candidate's own
    rows, so both sides moved together and it stayed green while only the row
    check caught the edit. Both must reject."""
    def f(p):
        r = p["rows"][1]
        old = copy.deepcopy(r["signature"])
        new = copy.deepcopy(old)
        new["st"] = [new["st"][0], new["st"][1] + 1, new["st"][2]]
        r["signature"] = new
        for i in p["identities"]:
            for fac in i["factors"]:
                if fac["signature"] == old:
                    fac["signature"] = copy.deepcopy(new)
            i["factors"].sort(key=S.factor_sort_key)
    return _mutate(base, f)


@case("P1_ARTIFACT_PINNED", reason="not the pinned")
def bad_pin_with_valid_packet_supplied(base, tmp):
    """A wrong P1 pin alongside a valid packet used to continue without the P1
    object and emit a SHORT registry. Fail-closed either way, but a truncated
    verdict set is worse to diagnose than a complete one."""
    src = pathlib.Path(base["p1_artifact"])
    dst = tmp / "p1.json"
    dst.write_bytes(src.read_bytes() + b"\n")
    p1 = p2_emit.load_p1(base["p1_artifact"])
    return {"p1_artifact": str(dst), "packet": p2_emit.build_packet(p1)}


# ------------------------ the fifth audit: two JSON-reachable raises, three
# expectations the candidate helped define, and the § 5.4 clauses the frozen
# convention map did not carry.

@case("IDY_RECOMPUTES_FROM_ROWS",
      reason="recompute, want",
      require_green=("ROWS_TRIPLES_NORMALIZED", "NESTED_SCHEMA_EXACT",
                     "IDY_FACTORS_RESOLVE_BY_SIGNATURE"))
def zero_row_value_under_a_negative_exponent(base, tmp):
    """A ZERO value on a row some identity inverts.

    Normalized (c > 0, gcd = 1) and shape-valid, so ROWS_TRIPLES_NORMALIZED
    passes it and the `norm_ok` gate licenses the exact recompute; Q(phi) then
    inverts zero and the run RAISED instead of returning verdicts. The gate was
    independent of the check it guarded but not sufficient for it. The row is
    found structurally, by which signatures carry a negative exponent."""
    def f(p):
        neg = {S.signature_of(fac) for i in p["identities"]
               for fac in i["factors"] if fac["exponent"] < 0}
        for r in p["rows"]:
            if S.signature_of(r) in neg:
                r["value"] = [0, 0, 1]
                return
    return _reparsed(base, f)


@case("ROWS_VALUES_MATCH_P1", "ROWS_SIGNATURES_MATCH_P1",
      reason={"ROWS_VALUES_MATCH_P1": "UNKNOWN_ROW",
              "ROWS_SIGNATURES_MATCH_P1": "UNKNOWN_ROW"},
      require_green=("IDXMAP_MATCHES_P1", "NESTED_SCHEMA_EXACT",
                     "LABELS_ARE_NOT_JOIN_KEYS"))
def row_label_absent_from_p1(base, tmp):
    """A label P1 does not carry. A string, so the nested schema accepts it.

    This RAISED KeyError out of the predicate while the indexing expectation
    joined through the candidate's labels. `require_green` on IDXMAP_MATCHES_P1
    is the assertion that matters now: a display label can no longer move the
    map's expectation at all, so only the two checks that bind labels to P1
    reject."""
    def f(p):
        p["rows"][0]["label"] = "UNKNOWN_ROW"
    return _mutate(base, f)


@case("IDXMAP_MATCHES_P1", "ROWS_COUNT", "IDXMAP_ZERO_BASED_CONTIGUOUS",
      reason={"IDXMAP_MATCHES_P1": "differs from the P1-derived",
              "IDXMAP_ZERO_BASED_CONTIGUOUS": "not zero-based contiguous"},
      require_green=("IDXMAP_BIJECTIVE", "IDXMAP_DOMAINS_DECLARED"))
def row_and_its_indexing_entry_dropped_together(base, tmp):
    """Eight rows and eight consistent entries. Contiguity, bijectivity and
    non-identity all still hold ACROSS THE SHRUNKEN SET, so every structural
    indexing check passed while the map described eight of nine rows. The
    expectation is nine entries because P1 has nine rows, not because the
    candidate does."""
    def f(p):
        p["rows"].pop()
        p["indexing_map"]["entries"] = [
            e for e in p["indexing_map"]["entries"]
            if e["source_position"] != len(p["rows"])]
    return _mutate(base, f)


@case("IDXMAP_MATCHES_P1", "ROWS_SIGNATURES_MATCH_P1",
      reason={"IDXMAP_MATCHES_P1": "differs from the P1-derived"},
      require_green=("IDXMAP_BIJECTIVE", "IDXMAP_ZERO_BASED_CONTIGUOUS",
                     "ROWS_VALUES_MATCH_P1", "IDXMAP_DISCRIMINATING"))
def labels_swapped_between_rows_sharing_a_value(base, tmp):
    """Two rows that carry the SAME value, labels exchanged, destinations
    following. Values still match P1 label-wise because the values are equal,
    so only the signature tie caught it while the map itself, joined through
    those labels, endorsed a correspondence that was wrong. Found structurally:
    any two rows with equal values."""
    def f(p):
        pairs = [(i, j) for i in range(len(p["rows"]))
                 for j in range(i + 1, len(p["rows"]))
                 if p["rows"][i]["value"] == p["rows"][j]["value"]]
        if not pairs:
            raise RuntimeError("no two rows share a value; case not applicable")
        i, j = pairs[0]
        p["rows"][i]["label"], p["rows"][j]["label"] = \
            p["rows"][j]["label"], p["rows"][i]["label"]
        ent = p["indexing_map"]["entries"]
        ent[i]["destination_position"], ent[j]["destination_position"] = \
            ent[j]["destination_position"], ent[i]["destination_position"]
    return _mutate(base, f)


@case("ADJUDICATES_MATCH_P1", "CONV_MAP_EXACT",
      reason={"CONV_MAP_EXACT": "differs from the frozen template"},
      require_green=("CONV_BASING_MATCHES_ADJUDICATES", "CONV_MAP_COMPLETE"))
def adjudicates_and_basing_moved_together(base, tmp):
    """A fabricated construction hash in BOTH `adjudicates` and the basing
    reference. The two agree with each other, so the basing check is satisfied;
    while the template's expectation was read off the candidate's own
    `adjudicates`, CONV_MAP_EXACT was satisfied too, and only one check stood
    between the packet and a construction packet nobody declared."""
    def f(p):
        h = p["adjudicates"]["construction_packet_sha256"]
        wrong = h[:-1] + ("0" if h[-1] != "0" else "1")
        p["adjudicates"]["construction_packet_sha256"] = wrong
        p["convention_map"]["basing_reference"]["construction_packet_sha256"] = wrong
    return _mutate(base, f)


@case("CONV_MAP_EXACT", reason="differs from the frozen template",
      require_green=("CONV_MAP_COMPLETE", "CONV_BASING_MATCHES_ADJUDICATES"))
def anchor_uniqueness_gate_dropped(base, tmp):
    """§ 5.4 makes uniqueness a GATE rather than an assumption, and names the
    INVALID ANCHOR outcome when both branches hold. The frozen template carried
    the anchor rule without its gate until the fifth audit."""
    def f(p):
        p["convention_map"]["orientation_anchor_rule"].pop("uniqueness_gate")
    return _mutate(base, f)


@case("CONV_MAP_EXACT", reason="differs from the frozen template",
      require_green=("CONV_MAP_COMPLETE",))
def bridge_involution_ordering_dropped(base, tmp):
    """The last clause of § 5.4's enumerated identification: the involution acts
    on T^2_target AFTER it, never before. An implementation that inverted first
    satisfies every other clause in the statement."""
    def f(p):
        s = p["convention_map"]["bridge"]["statement"]
        p["convention_map"]["bridge"]["statement"] = s.split("; and the involution")[0]
    return _mutate(base, f)


@case("HASH_RECORD_EXACT", "HASH_RECORD_SCHEMA_EXACT",
      reason={"HASH_RECORD_EXACT": "differs from the frozen template",
              "HASH_RECORD_SCHEMA_EXACT": "canonical_plaintext_byte_length"},
      require_green=("HASH_RECORD_MATCHES_PACKET", "HASH_RECORD_IS_CANDIDATE",
                     "CANONICAL_BYTES_FIXED_POINT"))
def sidecar_byte_length_is_a_float(base, tmp):
    """Python says 10974.0 == 10974, so a float byte length passed a check whose
    name is EXACT, and passed HASH_RECORD_MATCHES_PACKET too. Only the byte
    comparison sees it. Same class as the packet's 1 vs 1.0, one artifact over."""
    p1 = p2_emit.load_p1(base["p1_artifact"])
    pkt = p2_emit.build_packet(p1)
    rec = p2_emit.build_hash_record(S.canonical_bytes(pkt))
    rec["canonical_plaintext_byte_length"] = float(rec["canonical_plaintext_byte_length"])
    return {"p1_artifact": base["p1_artifact"], "packet": pkt, "hash_record": rec}


# ------------------------ the sixth audit: the § 5.4 branch-to-orientation
# mapping, and the last JSON-reachable raise.

def _rec_reparsed(base, fn):
    """Mutate the SIDECAR through a JSON round trip, so its leaf types are what
    a consumer parsing the file would hold rather than what Python happens to."""
    p1 = p2_emit.load_p1(base["p1_artifact"])
    pkt = p2_emit.build_packet(p1)
    pb = S.canonical_bytes(pkt)
    rec = json.loads(json.dumps(p2_emit.build_hash_record(pb)))
    fn(rec)
    return {"p1_artifact": base["p1_artifact"], "packet": pkt,
            "packet_bytes": pb, "hash_record": rec}


@case("HASH_RECORD_MATCHES_PACKET", reason="record says 0",
      require_green=("CANONICAL_BYTES_FIXED_POINT", "HASH_RECORD_IS_CANDIDATE",
                     "BYTES_BIND_TO_PACKET"))
def sidecar_hash_field_is_not_a_string(base, tmp):
    """An INTEGER where the sidecar's hash belongs. JSON-reachable, and it
    RAISED rather than returning verdicts: the detail sliced the field, and
    `.get(k, "")` defends only against a missing key, never a present one of the
    wrong type. HASH_RECORD_SCHEMA_EXACT saw the type and did not gate what
    followed, which is the not-sufficient-validation class again. null, true,
    a float and an object all reached the same line."""
    def f(rec):
        rec["expected_canonical_plaintext_sha256"] = 0
    return _rec_reparsed(base, f)


@case("CONV_MAP_EXACT", reason="differs from the frozen template",
      require_green=("CONV_MAP_COMPLETE", "CONV_BASING_MATCHES_ADJUDICATES"))
def anchor_selection_mapping_dropped(base, tmp):
    """§ 5.4: "whichever of the two orientations agrees with the analytic
    convention at R7 is the selected one". The rule named who selects and
    between what, the gate named when a selection is valid, and nothing named
    which branch chooses which orientation."""
    def f(p):
        p["convention_map"]["orientation_anchor_rule"].pop("selection_mapping")
    return _mutate(base, f)


@case("CONV_MAP_EXACT", reason="differs from the frozen template",
      require_green=("CONV_MAP_COMPLETE", "CONV_BASING_MATCHES_ADJUDICATES"))
def anchor_selection_mapping_reversed(base, tmp):
    """The mapping present but BACKWARDS: x = r selecting the global inverse.

    The decisive case of the pair. Dropping a field is caught by anything that
    compares field sets; a reversed mapping is caught only by comparing content,
    and it is the error an implementer could actually make."""
    def f(p):
        rule = p["convention_map"]["orientation_anchor_rule"]
        s = rule["selection_mapping"]
        rule["selection_mapping"] = (s.replace("COMMITTED table", "\x00")
                                      .replace("GLOBAL INVERSE", "COMMITTED table")
                                      .replace("\x00", "GLOBAL INVERSE"))
    return _mutate(base, f)


# ------------------------ the ninth audit: § 5.4's closing paragraph, `null` as
# a supplied value, details that asserted what had just failed, and an unbounded
# exponent.

@case("CONV_MAP_EXACT", reason="differs from the frozen template",
      require_green=("CONV_MAP_COMPLETE", "CONV_BASING_MATCHES_ADJUDICATES"))
def bridge_exclusivity_dropped(base, tmp):
    """§ 5.4: "The bridge is one involution, and nothing else". The clause that
    DEFINES the global inversion the anchor rule invokes by name, and forecloses
    every other transformation."""
    def f(p):
        p["convention_map"]["bridge"].pop("exclusivity")
    return _mutate(base, f)


@case("CONV_MAP_EXACT", reason="differs from the frozen template",
      require_green=("CONV_MAP_COMPLETE", "CONV_BASING_MATCHES_ADJUDICATES"))
def bridge_exclusivity_weakened(base, tmp):
    """Present but WEAKENED: rows inverted independently rather than all at once.

    The decisive one of the pair, and the error the clause exists to forbid: a
    per-row inversion satisfies every other clause in the packet."""
    def f(p):
        b = p["convention_map"]["bridge"]
        b["exclusivity"] = b["exclusivity"].replace(
            "applied to EVERY ROW AT ONCE", "applied to each row independently")
    return _mutate(base, f)


@case("TOP_LEVEL_TYPES", reason="packet", require_green=("P1_ARTIFACT_PINNED",))
def whole_document_null_packet(base, tmp):
    """A packet FILE whose entire content is `null`.

    It parses to None, and while None was the "not supplied" default this was
    silently replaced by the emitter's own packet and returned 38 GREEN verdicts
    about an object the caller never supplied. Every other whole-document scalar
    was rejected correctly, so the hole was the sentinel, not the type check."""
    p1 = p2_emit.load_p1(base["p1_artifact"])
    del p1
    return {"p1_artifact": base["p1_artifact"], "packet": json.loads("null")}


@case("TOP_LEVEL_TYPES", reason="hash_record",
      require_green=("KEYS_EXACTLY_SEVEN", "P1_ARTIFACT_PINNED"))
def whole_document_null_sidecar(base, tmp):
    """The same at the sidecar. `null` is a JSON document, and the contract
    covers any state either file can be in."""
    p1 = p2_emit.load_p1(base["p1_artifact"])
    return {"p1_artifact": base["p1_artifact"], "packet": p2_emit.build_packet(p1),
            "hash_record": json.loads("null")}


@case("NESTED_SCHEMA_EXACT", reason="magnitude exceeds",
      require_green=("TOP_LEVEL_TYPES", "KEYS_EXACTLY_SEVEN", "ROWS_NO_FLOATS"))
def factor_exponent_magnitude_absurd(base, tmp):
    """A 13-byte integer literal buying an unbounded computation.

    The exact recompute raised the row value to the power before anything
    considered magnitude, so this ended in MemoryError: a JSON-reachable raise
    from a shape-valid document. The bound is on cost, not on mathematics;
    IDY_DEFINITIONS_MATCH_P1 is what pins the exponents to P1's."""
    def f(p):
        p["identities"][0]["factors"][0]["exponent"] = 10 ** 6
    return _reparsed(base, f)


# ------------------------ the tenth audit: a diagnostic that pointed at the
# healthy half of what it was reporting.

@case("IDXMAP_ZERO_BASED_CONTIGUOUS", reason="destinations",
      require_green=("IDXMAP_BIJECTIVE", "IDXMAP_DOMAINS_DECLARED",
                     "IDXMAP_DISCRIMINATING"))
def indexing_destinations_not_contiguous(base, tmp):
    """A DESTINATION-side contiguity break, sources untouched.

    Still nine distinct destinations, so bijectivity holds and only contiguity
    reddens. The detail used to print the sources, which were fine, under a red
    verdict about the destinations. This case exists to keep the detail pointed
    at the half that actually failed."""
    def f(p):
        for i, e in enumerate(p["indexing_map"]["entries"]):
            e["destination_position"] = i if i < len(p["rows"]) - 1 else len(p["rows"])
    return _mutate(base, f)


# ------------------------ the eleventh audit: details naming a strict subset of
# their own predicate, and the two cost operands the exponent bound missed.

@case("IDXMAP_DOMAINS_DECLARED", reason="zero_based",
      require_green=("IDXMAP_BIJECTIVE", "IDXMAP_ZERO_BASED_CONTIGUOUS",
                     "IDXMAP_MATCHES_P1", "NESTED_SCHEMA_EXACT"))
def indexing_not_zero_based(base, tmp):
    """`zero_based: false`, schema-valid at a bool field.

    The third clause of a three-clause predicate, and the one that had no case:
    the detail printed the two domain names, both correct, and never named it."""
    def f(p):
        p["indexing_map"]["zero_based"] = False
    return _mutate(base, f)


@case("HASH_RECORD_MATCHES_PACKET", reason="byte length",
      require_green=("HASH_RECORD_SCHEMA_EXACT", "HASH_RECORD_IS_CANDIDATE",
                     "CANONICAL_BYTES_FIXED_POINT"))
def sidecar_byte_length_wrong_hash_right(base, tmp):
    """The length wrong, the hash right. The check tests both and used to print
    only the hash, so a red verdict about the record matching the packet was
    elaborated with the half that matched."""
    p1 = p2_emit.load_p1(base["p1_artifact"])
    pkt = p2_emit.build_packet(p1)
    rec = p2_emit.build_hash_record(S.canonical_bytes(pkt))
    rec["canonical_plaintext_byte_length"] += 1
    return {"p1_artifact": base["p1_artifact"], "packet": pkt, "hash_record": rec}


@case("NESTED_SCHEMA_EXACT", reason="exceeds",
      require_green=("TOP_LEVEL_TYPES", "KEYS_EXACTLY_SEVEN", "ROWS_NO_FLOATS"))
def factor_count_absurd(base, tmp):
    """Repeated factors: shape-valid, each resolving to exactly one row.

    Cost is exponent x factor count x value magnitude, and bounding the exponent
    alone left the product unbounded. A packet of a few tens of KB drove assess()
    past 400 seconds toward the same MemoryError the exponent bound was added to
    stop."""
    def f(p):
        slot = p["identities"][0]
        f0 = copy.deepcopy(slot["factors"][0])
        slot["factors"] = sorted(
            [copy.deepcopy(f0) for _ in range(S.FACTORS_PER_SLOT_MAX + 1)],
            key=S.factor_sort_key)
    return _reparsed(base, f)


# ------------------------ the twelfth audit: clauses and shape calls that could
# each be deleted with the battery still fully green. Coverage of a CHECK is not
# coverage of its CLAUSES.

@case("IDXMAP_DOMAINS_DECLARED", reason="source_domain",
      require_green=("IDXMAP_BIJECTIVE", "IDXMAP_ZERO_BASED_CONTIGUOUS",
                     "NESTED_SCHEMA_EXACT"))
def indexing_arbitrary_source_domain(base, tmp):
    """The source_domain clause had no case at all: only the destination side
    and `zero_based` were exercised, so that clause could be deleted outright
    with the battery still reporting full coverage."""
    def f(p):
        p["indexing_map"]["source_domain"] = "some_other_ordering"
    return _mutate(base, f)


@case("NESTED_SCHEMA_EXACT", reason="bool where int required",
      require_green=("TOP_LEVEL_TYPES", "KEYS_EXACTLY_SEVEN", "ROWS_NO_FLOATS"))
def boolean_at_a_scalar_int_field(base, tmp):
    """JSON `true` at a scalar integer field, as opposed to inside a triple.

    `exact_shape`'s bool-where-int rule is the sole defense at every scalar int
    in the schema, and nothing exercised it: the leaf cases all target triples.
    Python says True == 1, so without this rule the value is invisible
    downstream."""
    def f(p):
        p["identities"][0]["position"] = True
    return _reparsed(base, f)


@case("NESTED_SCHEMA_EXACT", reason="identities[0]",
      require_green=("TOP_LEVEL_TYPES", "KEYS_EXACTLY_SEVEN"))
def undeclared_field_on_an_identity_slot(base, tmp):
    """`exact_shape` on an identity slot had no case; only rows and factors did,
    so that call could be removed and an undeclared field on a slot would pass."""
    def f(p):
        p["identities"][0]["note"] = "undeclared"
    return _mutate(base, f)


@case("NESTED_SCHEMA_EXACT", reason="indexing_map.entries[0]",
      require_green=("TOP_LEVEL_TYPES", "KEYS_EXACTLY_SEVEN"))
def malformed_indexing_entry(base, tmp):
    """`exact_shape` on an indexing entry had no case. It is also what stops a
    missing entry field from making a later stage dereference it."""
    def f(p):
        del p["indexing_map"]["entries"][0]["destination_position"]
    return _mutate(base, f)


@case("REGISTRY", reason="more than once")
def registry_duplicate_id(base, tmp):
    """`validate_registry` has three branches and only `missing` had a case, and
    that case's reason matched all three stop strings. Duplicated ids."""
    orig = p2.validate_registry
    p2.validate_registry = lambda emitted: orig(list(emitted) + [list(emitted)[0]])
    return {"p1_artifact": base["p1_artifact"], "__restore__": orig}


@case("REGISTRY", reason="undeclared check ids")
def registry_undeclared_id(base, tmp):
    orig = p2.validate_registry
    p2.validate_registry = lambda emitted: orig(list(emitted) + ["NOT_A_DECLARED_CHECK"])
    return {"p1_artifact": base["p1_artifact"], "__restore__": orig}


@case("NESTED_SCHEMA_EXACT", reason="digits",
      require_green=("TOP_LEVEL_TYPES", "KEYS_EXACTLY_SEVEN", "ROWS_NO_FLOATS"))
def row_value_magnitude_absurd(base, tmp):
    """The THIRD cost operand. Exponent and factor count were bounded; the base
    was not, so a shape-valid packet no larger than the real one still drove the
    recompute to hundreds of seconds."""
    def f(p):
        big = int("9" * (S.VALUE_ABS_DIGITS_MAX + 1))
        p["rows"][0]["value"] = [big, 0, big - 1]
    return _reparsed(base, f)


@case("CLASS_SELECTOR_NOT_SELF_INVERSE", reason="ZERO",
      require_green=("CLASS_CENSUS", "ROWS_TRIPLES_NORMALIZED"))
def selector_value_is_zero(base, tmp):
    """Zero is a normalized triple with no inverse, so the reciprocal class
    {x, x^-1} § 5.2 requires does not exist. Testing only v^2 = 1 called that
    "not self-inverse", which is true and beside the point."""
    def f(p):
        for r in p["rows"]:
            if r["evidentiary_class"] == S.CLASS_SELECTOR:
                r["value"] = [0, 0, 1]
    return _reparsed(base, f)


# ------------------------ the thirteenth audit: clauses that no case pinned.

# The frozen leaf inventory of the convention map, owned HERE, on the CHECKER
# side, and never read off the template.
#
# The first version of this generated its cases FROM `S.convention_map_template`,
# which reproduced, inside the fix for a coverage gap, the exact
# expectation-ownership defect this build has been closing for eight rounds:
# deleting a leaf from the template simply generated one case fewer, so the
# battery stayed fully green on a packet that had lost a clause. Only hollowing
# the template caught it.
#
# Owned on this side, deleting a leaf from the emitter's template makes the
# matching case's `.pop()` raise, and adding one trips the inventory check
# below. Either way it takes an edit in two files that disagree.
EXPECTED_CONVENTION_LEAVES = (
    ("basing_reference", "construction_packet_sha256"),
    ("basing_reference", "protocol_section"),
    ("basing_reference", "source"),
    ("bridge", "exclusivity"),
    ("bridge", "identification"),
    ("bridge", "protocol_section"),
    ("bridge", "statement"),
    ("orientation_anchor_rule", "anchor_row_label"),
    ("orientation_anchor_rule", "anchor_row_signature"),
    ("orientation_anchor_rule", "protocol_section"),
    ("orientation_anchor_rule", "rule"),
    ("orientation_anchor_rule", "selection_mapping"),
    ("orientation_anchor_rule", "uniqueness_gate"),
)


# The frozen CONTENT of both templates, digested, owned on the CHECKER side.
#
# The leaf inventory above pins the convention map's leaf NAMES. It does not pin
# what those leaves SAY, and the sidecar template was pinned by nothing at all:
# both live in the same file the emitter builds from, so a single-file edit moved
# the emitted artifact and its expectation together, which is the defect this
# build has closed eight times elsewhere. Gutting a leaf's content, or rewriting
# the sidecar's quarantine banner, left the driver green and the battery fully
# green.
#
# A digest pins content exhaustively and cheaply, in both directions, for both
# templates, and does not need a case per string. Regenerate deliberately when a
# template legitimately changes: that is the point, since it forces the edit to
# be made twice, in two files, by someone who means it.
TEMPLATE_DIGESTS = {
    "convention_map_template":
        "72e5728756cf61c46f41d74bc39ab1cec34d4cba210a16e4209ba82bab9b4b96",
    "hash_record_template":
        "1d70ea7ac27a97d1e412c057c128eee08175291962b118e8a6dc1f311df131d5",
}


def template_digests_match():
    """Both frozen templates against their frozen digests, dynamics normalised."""
    import hashlib
    got = {
        "convention_map_template": S.convention_map_template(
            "0" * 64, S.SELECTOR_LABEL,
            {"dimension": 1, "s": [1, 0, 1], "t": [1, 0, 1], "st": [1, 0, 1]}),
        "hash_record_template": S.hash_record_template("0" * 64, 0),
    }
    return [f"{k}: {hashlib.sha256(S.canonical_bytes(v)).hexdigest()[:16]} != "
            f"{TEMPLATE_DIGESTS[k][:16]}"
            for k, v in got.items()
            if hashlib.sha256(S.canonical_bytes(v)).hexdigest() != TEMPLATE_DIGESTS[k]]


def convention_leaf_inventory_matches():
    """The template's actual leaf set against the frozen inventory above.

    Catches a leaf ADDED to the template, which no per-leaf deletion case can
    see. Reported as a REGISTRY-class stop, because an unpinned clause is a
    coverage defect in the battery rather than a defect in the packet."""
    probe = S.convention_map_template(
        "0" * 64, S.SELECTOR_LABEL,
        {"dimension": 1, "s": [1, 0, 1], "t": [1, 0, 1], "st": [1, 0, 1]})
    actual = {(sec, leaf) for sec in probe for leaf in probe[sec]}
    return sorted(actual - set(EXPECTED_CONVENTION_LEAVES)), \
        sorted(set(EXPECTED_CONVENTION_LEAVES) - actual)


def _register_convention_leaf_cases():
    """One DELETION case per FROZEN leaf. Deletion, not replacement: with the
    key absent from the template a replacement merely adds an undeclared key,
    CONV_MAP_EXACT still reddens for the wrong reason, and the case passes,
    which is why three leaves that had cases were nonetheless unpinned."""
    for section, leaf in EXPECTED_CONVENTION_LEAVES:
        def make(sec, lf):
            def fn(base, tmp):
                return _mutate(base, lambda p: p["convention_map"][sec].pop(lf))
            fn.__name__ = f"convention_leaf_dropped__{sec}__{lf}"
            fn.__doc__ = (f"§ 4.4 requires the convention map to carry "
                          f"{sec}.{lf}; dropping it must redden CONV_MAP_EXACT.")
            return fn
        case("CONV_MAP_EXACT", reason="differs from the frozen template",
             require_green=("CONV_MAP_COMPLETE",))(make(section, leaf))


_register_convention_leaf_cases()


@case("TOP_LEVEL_TYPES", reason="identities", require_green=("KEYS_EXACTLY_SEVEN",))
def identities_is_not_a_list(base, tmp):
    """`identities: null`.

    Its entry in the frozen top-level type map could be deleted with the battery
    fully green, which is why this case exists. An earlier version of this
    docstring, and of P2_REPORT § 4l, said it was "the one of three containers
    with no case"; that was false, since `rows` and `convention_map` have cases
    and so does this one now. What was true is the deletability, and that is all
    this claims."""
    return _mutate(base, lambda p: p.update({"identities": None}))


@case("NESTED_SCHEMA_EXACT", reason="factors[0].signature",
      require_green=("TOP_LEVEL_TYPES", "KEYS_EXACTLY_SEVEN"))
def undeclared_field_in_a_factor_signature(base, tmp):
    """The `exact_shape` call on a FACTOR's signature had no case, though it is
    the sole defence against an undeclared field there."""
    def f(p):
        p["identities"][0]["factors"][0]["signature"]["note"] = "undeclared"
    return _mutate(base, f)


@case("CONV_BASING_MATCHES_ADJUDICATES", reason="DIFFERENT construction packet",
      require_green=("CONV_MAP_COMPLETE", "TOP_LEVEL_TYPES"))
def basing_reference_is_not_an_object(base, tmp):
    """A non-object `basing_reference`. The isinstance guard that makes this a
    verdict rather than an AttributeError had no case."""
    def f(p):
        p["convention_map"]["basing_reference"] = "not an object"
    return _mutate(base, f)


@case("IDY_POSITIONS_CONTIGUOUS", "IDY_COUNT",
      reason={"IDY_POSITIONS_CONTIGUOUS": "0, 1, 2, 3, 4"},
      require_green=("NESTED_SCHEMA_EXACT",))
def identity_slot_added_with_contiguous_positions(base, tmp):
    """A FIFTH slot, positions still contiguous from zero.

    This is the case that pins the twelfth audit's repair: while the expected
    range came from the candidate's own len(idents), a five-slot packet with
    positions 0..4 satisfied the check. Against the frozen IDENTITY_SLOTS it
    cannot. Without this case the repair was invisible and revertible."""
    def f(p):
        extra = copy.deepcopy(p["identities"][0])
        extra["position"] = len(p["identities"])
        p["identities"].append(extra)
    return _mutate(base, f)


@case("CONV_MAP_EXACT", reason="differs from the frozen template",
      require_green=("CONV_MAP_COMPLETE", "TOP_LEVEL_TYPES", "KEYS_EXACTLY_SEVEN"))
def deeply_nested_value_in_the_convention_map(base, tmp):
    """2000 nested empty lists where a convention-map leaf belongs.

    This pins the twelfth audit's iterative `_find_floats`. That walk runs ahead
    of every gate and outside the exception boundary, so while it recursed this
    input escaped assess() with RecursionError and no registry at all; the
    battery could not see the repair and could not see it reverted. The nested
    schema does not descend into `bridge`, so the packet reaches the float scan
    intact."""
    def f(p):
        node = []
        for _ in range(2000):
            node = [node]
        p["convention_map"]["bridge"]["statement"] = node
    return _mutate(base, f)


# ------------------------ the fourteenth audit: repairs nothing was pinning.

@case("CANONICAL_BYTES_FIXED_POINT", reason="hash record",
      require_green=("BYTES_BIND_TO_PACKET", "HASH_RECORD_IS_CANDIDATE"))
def sidecar_carries_a_non_json_constant(base, tmp):
    """A sidecar FILE carrying `NaN`, which json.loads accepts and JSON does not.

    The record half of this check was described in a comment as "a fixed point
    by construction", and that stopped being true the moment `parse_constant`
    was added: this is a live JSON-reachable failure branch, and deleting
    `and rec_ok` left the battery fully green."""
    p1 = p2_emit.load_p1(base["p1_artifact"])
    pkt = p2_emit.build_packet(p1)
    pb = S.canonical_bytes(pkt)
    rec = p2_emit.build_hash_record(pb)
    rec["canonical_plaintext_byte_length"] = float("nan")
    return {"p1_artifact": base["p1_artifact"], "packet": pkt,
            "packet_bytes": pb, "hash_record": rec}


@case("CLASS_SELECTOR_NOT_SELF_INVERSE", reason="designated selector signature",
      require_green=("CLASS_CENSUS", "ROWS_TRIPLES_NORMALIZED"))
def selector_signature_absent_from_the_rows(base, tmp):
    """No row carries P1's designated selector signature.

    This pins the thirteenth audit's repair. While the subject was whichever row
    CLAIMED the class, altering the designated row's signature left this check
    examining some other row's value and passing; against P1's signature there is
    no subject at all, and the check says so."""
    def f(p):
        for r in p["rows"]:
            if r["evidentiary_class"] == S.CLASS_SELECTOR:
                r["signature"]["st"] = [r["signature"]["st"][0] + 7,
                                        r["signature"]["st"][1],
                                        r["signature"]["st"][2]]
    return _mutate(base, f)


@case("IDXMAP_ZERO_BASED_CONTIGUOUS", reason="destinations",
      require_green=("IDXMAP_BIJECTIVE", "IDXMAP_DOMAINS_DECLARED"))
def indexing_entries_contiguous_but_short(base, tmp):
    """Eight entries over 0..7 on a nine-row packet.

    Pins the eleventh audit's frozen range: against `len(rows)` this was green
    whenever the entry array and the row array shrank together."""
    def f(p):
        p["indexing_map"]["entries"] = [e for e in p["indexing_map"]["entries"]
                                        if e["source_position"] != 8]
    return _mutate(base, f)


@case("CANONICAL_BYTES_FIXED_POINT", reason="not ASCII",
      require_green=("HASH_RECORD_MATCHES_PACKET", "HASH_RECORD_IS_CANDIDATE"))
def packet_bytes_not_ascii(base, tmp):
    """Non-ASCII bytes. § 10.2 requires ASCII, and that branch had no case: only
    the final byte-equality comparison of six branches was pinned."""
    p1 = p2_emit.load_p1(base["p1_artifact"])
    pkt = p2_emit.build_packet(p1)
    bad = S.canonical_bytes(pkt)[:-1] + "\u00e9".encode("utf-8") + b"\n"
    return {"p1_artifact": base["p1_artifact"], "packet": pkt, "packet_bytes": bad,
            "hash_record": p2_emit.build_hash_record(bad)}


@case("CANONICAL_BYTES_FIXED_POINT", reason="contains CR",
      require_green=("HASH_RECORD_MATCHES_PACKET",))
def packet_bytes_contain_cr(base, tmp):
    """§ 10.2 requires LF. That branch had no case either."""
    p1 = p2_emit.load_p1(base["p1_artifact"])
    pkt = p2_emit.build_packet(p1)
    bad = S.canonical_bytes(pkt).replace(b"\n", b"\r\n", 1)
    return {"p1_artifact": base["p1_artifact"], "packet": pkt, "packet_bytes": bad,
            "hash_record": p2_emit.build_hash_record(bad)}


@case("CANONICAL_BYTES_FIXED_POINT", reason="trailing newline",
      require_green=("HASH_RECORD_MATCHES_PACKET",))
def packet_bytes_double_trailing_newline(base, tmp):
    """§ 10.2 requires exactly one."""
    p1 = p2_emit.load_p1(base["p1_artifact"])
    pkt = p2_emit.build_packet(p1)
    bad = S.canonical_bytes(pkt) + b"\n"
    return {"p1_artifact": base["p1_artifact"], "packet": pkt, "packet_bytes": bad,
            "hash_record": p2_emit.build_hash_record(bad)}


@case("IDY_RECOMPUTES_FROM_ROWS", "IDY_FACTORS_RESOLVE_BY_SIGNATURE",
      reason={"IDY_RECOMPUTES_FROM_ROWS": "BLOCKED: factors do not resolve"},
      allow_blocked=True, require_green=("NESTED_SCHEMA_EXACT", "IDY_COUNT"))
def recompute_blocked_when_factors_do_not_resolve(base, tmp):
    """Pins the BLOCKED branch's WORDING, not merely that some verdict appears.

    Reverting it to the old summary line reported an ordinary red reading
    "0/0 recompute, want 4 of 4" when nothing had been compared at all."""
    def f(p):
        p["identities"][0]["factors"][0]["signature"]["dimension"] = 99
    return _mutate(base, f)


def main():
    if len(sys.argv) == 2 and sys.argv[1] == "--counts":
        named = {t for _n, ts, _f, _r, _g, _b in CASES for t in ts}
        print(f"declared_checks={len(p2.CHECK_IDS)}")
        print(f"mutation_cases={len(CASES)}")
        print(f"cases_with_causal_reason={sum(1 for c in CASES if c[3] is not None)}")
        print(f"checks_named_by_a_case={len(named & set(p2.CHECK_IDS))}")
        return 0
    if len(sys.argv) != 2:
        print("usage: p2_mutations.py <p1_artifact>", file=sys.stderr)
        print("       p2_mutations.py --counts", file=sys.stderr)
        return 2
    base = {"p1_artifact": sys.argv[1]}

    covered = {t for _n, ts, _f, _r, _g, _b in CASES for t in ts}
    unknown = sorted(covered - set(p2.CHECK_IDS) - EXTRA_TARGETS)

    tmp = pathlib.Path(tempfile.mkdtemp(prefix="p2mut-"))
    print("P2 mutation battery: each case must make its TARGET check reject,\n"
          "                     by the branch it names\n")
    results = []
    try:
        for name, targets, fn, reason, req_green, allow_blocked in CASES:
            d = tmp / name
            d.mkdir()
            kw = fn(base, d)
            restore = kw.pop("__restore__", None)
            try:
                bad, details = failed_with_details(**kw)
            finally:
                if restore is not None:
                    p2.validate_registry = restore
            raised = isinstance(bad, dict) and "__raised__" in bad
            ok = case_passed(set(targets), bad if not raised else set(), raised,
                             reason=reason, details=details,
                             require_green=req_green, allow_blocked=allow_blocked)
            results.append((name, targets, ok, bad))
            print(f"  {'RED ' if ok else 'MISS'}  {name}")
            print(f"         targets: {', '.join(targets)}")
            if raised:
                print(f"         RAISED instead of rejecting: {bad['__raised__']}")
            elif ok:
                extra = sorted(set(bad) - set(targets))
                why = f' via "{reason}"' if reason else ""
                print(f"         rejected: {sorted(targets)}{why}"
                      + (f"  (also: {extra})" if extra else ""))
            else:
                still = sorted(set(targets) - set(bad))
                greenv = [g for g in req_green if g in bad]
                if still:
                    print(f"         STILL GREEN: {still}; failures were {sorted(bad)}")
                elif greenv:
                    print(f"         TRIPPED AN EARLIER GATE instead: {greenv}")
                else:
                    print(f"         WRONG REASON. Actual: "
                          + "; ".join(f"{t}={details.get(t)!r}" for t in targets))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    extra_leaves, missing_leaves = convention_leaf_inventory_matches()
    if extra_leaves or missing_leaves:
        print(f"\n  CONVENTION TEMPLATE LEAF INVENTORY MISMATCH: "
              f"unpinned {extra_leaves}, stale {missing_leaves}")
    digest_bad = template_digests_match()
    if digest_bad:
        print(f"\n  FROZEN TEMPLATE CONTENT CHANGED: {digest_bad}")

    observed = set()
    for _n, targets, ok, _b in results:
        if ok:
            observed |= set(targets)
    never = [c for c in p2.CHECK_IDS if c not in observed]
    red = sum(1 for r in results if r[2])
    print(f"\n  {red}/{len(results)} cases made EVERY named target reject")
    print(f"  coverage: {len(observed & set(p2.CHECK_IDS))}/{len(p2.CHECK_IDS)} "
          f"declared checks OBSERVED rejecting")
    if unknown:
        print(f"  CASES TARGET UNDECLARED IDS: {unknown}")
    if never:
        print(f"  CHECKS NEVER OBSERVED REJECTING: {never}")
    ok = (red == len(results) and not never and not unknown
          and not extra_leaves and not missing_leaves and not digest_bad)
    print(f"  {'EVERY DECLARED P2 CHECK IS OPERATIVE' if ok else 'COVERAGE INCOMPLETE'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
