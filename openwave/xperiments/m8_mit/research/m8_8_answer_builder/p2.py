"""P2 driver: emit the CANDIDATE answer packet and check it structurally.

Same architecture as P1, deliberately, because every part of it was paid for:

  * every check has a STABLE ID, and `assess()` is the single acceptance
    predicate the mutation battery calls;
  * registry integrity is FATAL, not a printed note;
  * every stage sits behind an exception boundary, because an exception is not
    a verdict;
  * publication is atomic and invalidates prior artifacts BEFORE the checks can
    stop, so a failed run cannot leave a stale green packet standing;
  * the input is PINNED and the build fails closed if it is not the pinned one.

The CONTRACT, scoped deliberately. `assess()` returns a complete verdict set,
never a raise, for any state THE ANSWER PACKET OR ITS SIDECAR can actually be
in: anything reachable by parsing either JSON file, whatever its content. Both
artifacts, because both are assessed here and the sidecar is where the last
raise lived; scoping this to the packet alone would have put a real defect
outside the promised surface on a technicality.

It does NOT promise verdicts for arbitrary Python objects handed to its keyword
surface in process, which is not a state either artifact can be in. TWO such
inputs are known to raise, a non-string dict key and an unserializable value
inside `convention_map`, and are left raising on purpose: hardening them would
make the predicate's obligation the whole of Python, and that requirement has no
end. A non-`bytes` `packet_bytes` was the third until the `null`-sentinel repair
required guarding it; it now returns BLOCKED verdicts. The count is stated here
because it was three, and a fix that changes it must change this sentence.

Fixed because they ARE reachable by parsing a file: a zero row value, an unknown
row label, and a sidecar hash field of the wrong type.

The output is a CANDIDATE. Nothing here is issued.
"""

import argparse
import hashlib
import json
import pathlib
import sys

import p2_emit
import p2_schema as S
from qphi_exact import Phi

CHECK_IDS = [
    "P1_ARTIFACT_PINNED",
    "KEYS_EXACTLY_SEVEN", "NESTED_SCHEMA_EXACT", "FORMAT_VERSION", "TARGET_ID",
    "ADJUDICATES_MATCH_P1",
    "ROWS_COUNT", "ROWS_IN_CANONICAL_ORDER", "ROWS_SIGNATURES_DISTINCT",
    "ROWS_TRIPLES_NORMALIZED", "ROWS_NO_FLOATS", "ROWS_VALUES_MATCH_P1",
    "ROWS_SIGNATURES_MATCH_P1", "CLASS_CENSUS", "CLASS_ASSIGNMENT_MATCHES_P1",
    "CLASS_SELECTOR_NOT_SELF_INVERSE",
    "IDY_COUNT", "IDY_POSITIONS_CONTIGUOUS", "IDY_FACTORS_RESOLVE_BY_SIGNATURE",
    "IDY_RECOMPUTES_FROM_ROWS", "IDY_FACTOR_ORDER_CANONICAL",
    "IDY_DEFINITIONS_MATCH_P1",
    "IDXMAP_DOMAINS_DECLARED", "IDXMAP_ZERO_BASED_CONTIGUOUS", "IDXMAP_BIJECTIVE",
    "IDXMAP_DISCRIMINATING", "IDXMAP_MATCHES_P1",
    "CONV_MAP_COMPLETE", "CONV_MAP_EXACT",
    "CONV_REPRESENTATION_EVALUATION_MATCHES_4_2",
    "CONV_BASING_MATCHES_ADJUDICATES",
    "LABELS_ARE_NOT_JOIN_KEYS",
    "TOP_LEVEL_TYPES", "BYTES_BIND_TO_PACKET",
    "CANONICAL_BYTES_FIXED_POINT", "HASH_RECORD_SCHEMA_EXACT", "HASH_RECORD_EXACT",
    "HASH_RECORD_MATCHES_PACKET", "HASH_RECORD_IS_CANDIDATE",
]


def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()


def environment():
    import platform
    out = {"python": platform.python_version(), "platform": platform.platform()}
    for mod in ("numpy", "mpmath"):
        try:
            out[mod] = __import__(mod).__version__
        except Exception:                                     # noqa: BLE001
            out[mod] = "absent"
    return out


def validate_registry(emitted):
    """Registry integrity as a function, so the battery can drive its input."""
    missing = [c for c in CHECK_IDS if c not in emitted]
    dupes = sorted({c for c in emitted if emitted.count(c) > 1})
    undeclared = sorted(set(emitted) - set(CHECK_IDS))
    stops = []
    if missing:
        stops.append(f"REGISTRY: declared checks not executed: {missing}")
    if dupes:
        stops.append(f"REGISTRY: check ids emitted more than once: {dupes}")
    if undeclared:
        stops.append(f"REGISTRY: undeclared check ids emitted: {undeclared}")
    return stops, {"declared": len(CHECK_IDS), "emitted": len(emitted),
                   "missing": missing, "duplicated": dupes, "undeclared": undeclared}


def as_phi(t):
    a, b, c = t
    return Phi(a, b) / Phi(c)


PACKET_CHECKS = [c for c in CHECK_IDS if c != "P1_ARTIFACT_PINNED"]


_UNSET = object()


def assess(p1_artifact, packet=_UNSET, hash_record=_UNSET, packet_bytes=_UNSET):
    """The acceptance predicate. Returns (checks, stops, data).

    packet/hash_record may be supplied to check a candidate that was NOT built
    by this call, which is how the mutation battery reaches structural checks
    without rebuilding a valid packet first.

    `None` is a SUPPLIED VALUE, not an absent one. While these defaulted to
    None, a packet or sidecar file whose entire content was `null` parsed to
    None, was read as "nothing supplied", and was silently replaced by this
    module's own freshly built artifact: the run then returned a fully green
    verdict set describing an object the caller never provided. That is "the
    object assessed was not the object under review" recurring at the parameter
    surface rather than inside the objects, which is why the content sweeps
    could not see it. A private sentinel is the only default that distinguishes
    absent from null."""
    checks, stops, data = [], [], {}

    def add(cid, ok, detail):
        checks.append((cid, ok, detail))
        if not ok:
            stops.append(f"{cid}: {detail}")
        return ok

    def blocked(cids, e):
        for c in cids:
            add(c, False, f"BLOCKED: {type(e).__name__}: {e}")

    # the sole input version, fail closed --------------------------------
    p1 = None
    try:
        p1 = p2_emit.load_p1(p1_artifact)
        add("P1_ARTIFACT_PINNED", True,
            f"{p2_emit.P1_ARTIFACT_SHA256[:16]} as pinned")
    except Exception as e:                                    # noqa: BLE001
        add("P1_ARTIFACT_PINNED", False, f"{type(e).__name__}: {e}")

    # A failed pin blocks EVERYTHING, supplied packet or not. Continuing without
    # the P1 object produced a short registry: fail-closed, but the contract says
    # a complete verdict set, and a truncated one is worse to diagnose.
    if p1 is None:
        blocked(PACKET_CHECKS,
                RuntimeError("P1 pin failed; no packet check was evaluated"))
        reg, data["registry"] = validate_registry([c for c, _, _ in checks])
        stops += reg
        data["environment"] = environment()
        return checks, stops, data

    try:
        if packet is _UNSET:
            packet = p2_emit.build_packet(p1)
        pbytes = S.canonical_bytes(packet) if packet_bytes is _UNSET else packet_bytes
        # Required BY the sentinel change, not a reopening of the declined
        # keyword-surface hardening: making None a supplied value means a
        # supplied None reaches _canonical_conformance, which would raise. One
        # isinstance keeps the change verdict-neutral instead of trading one
        # raise for another.
        if not isinstance(pbytes, bytes):
            raise TypeError(f"packet_bytes: {type(pbytes).__name__}, want bytes")
        if hash_record is _UNSET:
            hash_record = p2_emit.build_hash_record(pbytes)
        data["packet"] = packet
        data["packet_bytes"] = pbytes
        data["hash_record"] = hash_record
        data["hash_record_bytes"] = S.canonical_bytes(hash_record)
    except Exception as e:                                    # noqa: BLE001
        blocked(PACKET_CHECKS, e)
        reg, data["registry"] = validate_registry([c for c, _, _ in checks])
        stops += reg
        data["environment"] = environment()
        return checks, stops, data

    # top-level shape -----------------------------------------------------
    # The packet's OWN type first: `sorted()` succeeds on a list and `.get()`
    # then raises, so even KEYS_EXACTLY_SEVEN cannot run until this holds.
    if not isinstance(packet, dict):
        add("TOP_LEVEL_TYPES", False,
            f"packet: {type(packet).__name__}, want dict")
        blocked([c for c in PACKET_CHECKS if c != "TOP_LEVEL_TYPES"],
                RuntimeError("packet is not an object"))
        reg, data["registry"] = validate_registry([c for c, _, _ in checks])
        stops += reg
        data["environment"] = environment()
        return checks, stops, data

    # A whole-packet scan, not a descent: safe on any structure, so it runs
    # before the nested pass can block it. § 4.4 says never a decimal, anywhere.
    floats = _find_floats(packet)
    add("ROWS_NO_FLOATS", not floats,
        "no floating-point anywhere in the packet" if not floats
        else f"floats at {floats[:3]}")

    keys = tuple(sorted(packet))
    add("KEYS_EXACTLY_SEVEN", keys == S.PACKET_KEYS,
        f"{len(keys)} keys" if keys == S.PACKET_KEYS
        else f"got {keys}, want {S.PACKET_KEYS}")
    add("FORMAT_VERSION", packet.get("format_version") == S.FORMAT_VERSION,
        packet.get("format_version"))
    add("TARGET_ID", packet.get("target_id") == S.TARGET_ID, packet.get("target_id"))

    # TOP-LEVEL CONTAINER TYPES FIRST. The nested pass iterates `rows`,
    # `identities` and `indexing_map`, so `rows: None` raised before any verdict
    # existed. Nothing may be iterated or `.get()`-chained until its container
    # type is established.
    tl = [f"{k}: {type(packet.get(k)).__name__}, want {ty.__name__}"
          for k, ty in S.TOP_LEVEL_TYPES.items()
          if not isinstance(packet.get(k), ty)]
    if not isinstance(hash_record, dict):
        tl.append(f"hash_record: {type(hash_record).__name__}, want dict")
    tl_ok = add("TOP_LEVEL_TYPES", not tl,
                "every top-level container has its declared type" if not tl
                else f"{len(tl)} wrong: {tl[:2]}")
    if not tl_ok:
        blocked([c for c in PACKET_CHECKS
                 if c not in ("KEYS_EXACTLY_SEVEN", "TOP_LEVEL_TYPES",
                              "FORMAT_VERSION", "TARGET_ID", "ROWS_NO_FLOATS")],
                RuntimeError("top-level types failed; dependent checks not evaluated"))
        reg, data["registry"] = validate_registry([c for c, _, _ in checks])
        stops += reg
        data["environment"] = environment()
        return checks, stops, data

    # EXACT NESTED SHAPE, before anything dereferences a nested field. Without
    # this an undeclared field on a row passed every check, and a MISSING nested
    # field made a later stage raise instead of returning a verdict.
    # NON-DEREFERENCING. exact_shape correctly reports "not an object", and the
    # NEXT line then called .get() on that same non-object and raised. Every
    # descent is now guarded by the type check that precedes it.
    shape = []
    shape += S.exact_shape(packet.get("adjudicates"), S.ADJUDICATES_KEYS, "adjudicates")

    for n, r in enumerate(packet.get("rows", [])):
        shape += S.exact_shape(r, S.ROW_KEYS, f"rows[{n}]")
        if isinstance(r, dict):
            shape += S.exact_shape(r.get("signature"), S.SIGNATURE_KEYS,
                                   f"rows[{n}].signature")
            shape += S.triple_shape_problems(r.get("value"), f"rows[{n}].value")
            shape += S.signature_leaf_problems(r.get("signature"),
                                               f"rows[{n}].signature")

    for n, i in enumerate(packet.get("identities", [])):
        shape += S.exact_shape(i, S.IDENTITY_KEYS, f"identities[{n}]")
        if not isinstance(i, dict):
            continue
        shape += S.triple_shape_problems(i.get("expected_value"),
                                         f"identities[{n}].expected_value")
        factors = i.get("factors")
        if isinstance(factors, list):
            if len(factors) > S.FACTORS_PER_SLOT_MAX:
                shape.append(f"identities[{n}].factors: {len(factors)} factors, "
                             f"exceeds {S.FACTORS_PER_SLOT_MAX}")
            for m, f in enumerate(factors):
                shape += S.exact_shape(f, S.FACTOR_KEYS,
                                       f"identities[{n}].factors[{m}]")
                if isinstance(f, dict):
                    shape += S.exact_shape(f.get("signature"), S.SIGNATURE_KEYS,
                                           f"identities[{n}].factors[{m}].signature")
                    shape += S.signature_leaf_problems(
                        f.get("signature"), f"identities[{n}].factors[{m}].signature")
                    if type(f.get("exponent")) is not int:
                        shape.append(f"identities[{n}].factors[{m}].exponent: "
                                     "not a true int")
                    elif abs(f["exponent"]) > S.EXPONENT_ABS_MAX:
                        shape.append(f"identities[{n}].factors[{m}].exponent: "
                                     f"magnitude exceeds {S.EXPONENT_ABS_MAX}")

    car = packet.get("convention_map")
    if isinstance(car, dict):
        anch = car.get("orientation_anchor_rule")
        if isinstance(anch, dict):
            shape += S.signature_leaf_problems(
                anch.get("anchor_row_signature"),
                "convention_map.orientation_anchor_rule.anchor_row_signature")

    im0 = packet.get("indexing_map")
    shape += S.exact_shape(im0, S.INDEXING_KEYS, "indexing_map")
    if isinstance(im0, dict):
        entries = im0.get("entries")
        if isinstance(entries, list):
            for n, e in enumerate(entries):
                shape += S.exact_shape(e, S.INDEXING_ENTRY_KEYS,
                                       f"indexing_map.entries[{n}]")

    shape_ok = add("NESTED_SCHEMA_EXACT", not shape,
                   "every nested object has exactly its declared fields and types"
                   if not shape else f"{len(shape)} problem(s): {shape[:2]}")

    if not shape_ok:
        blocked([c for c in PACKET_CHECKS
                 if c not in ("KEYS_EXACTLY_SEVEN", "NESTED_SCHEMA_EXACT",
                              "TOP_LEVEL_TYPES", "FORMAT_VERSION", "TARGET_ID",
                              "ROWS_NO_FLOATS")],
                RuntimeError("nested schema failed; dependent checks not evaluated"))
        reg, data["registry"] = validate_registry([c for c, _, _ in checks])
        stops += reg
        data["environment"] = environment()
        return checks, stops, data

    adj = packet.get("adjudicates", {})
    want_adj = ({"group_packet_sha256": p1["adjudicates"]["group_packet"],
                 "construction_packet_sha256": p1["adjudicates"]["construction_packet"]}
                if p1 else None)
    add("ADJUDICATES_MATCH_P1", want_adj is not None and adj == want_adj,
        "both packet hashes match the P1 artifact" if adj == want_adj else adj)

    # rows ----------------------------------------------------------------
    rows = packet.get("rows", [])
    add("ROWS_COUNT", len(rows) == S.ROW_COUNT, len(rows))

    ordered = sorted(rows, key=S.row_sort_key)
    add("ROWS_IN_CANONICAL_ORDER", rows == ordered,
        "already in row_sort_key order" if rows == ordered
        else "array order differs from the frozen sort key")

    sigs = [S.signature_of(r) for r in rows]
    distinct_ok = len(set(sigs)) == len(sigs) == S.ROW_COUNT
    add("ROWS_SIGNATURES_DISTINCT", distinct_ok,
        f"{len(set(sigs))} distinct of {len(sigs)}" if distinct_ok
        else f"{len(set(sigs))} distinct of {len(sigs)}, want "
             f"{S.ROW_COUNT} of {S.ROW_COUNT}")

    trips = [r["value"] for r in rows] + [c for r in rows
                                          for c in (r["signature"]["s"],
                                                    r["signature"]["t"],
                                                    r["signature"]["st"])]
    bad_norm = [t for t in trips if not S.is_normalized_triple(t)]
    norm_ok = add("ROWS_TRIPLES_NORMALIZED", not bad_norm,
                  f"{len(trips)} triples normalized" if not bad_norm
                  else f"{len(bad_norm)} not normalized")

    if p1:
        p1_vals = {v["label"]: list(v["triple"]) for v in p1["values"]}
        mismatch = [r["label"] for r in rows if r["value"] != p1_vals.get(r["label"])]
        add("ROWS_VALUES_MATCH_P1", not mismatch,
            mismatch or f"{len(rows)}/{S.ROW_COUNT} values match P1")
        p1_sig = {k: (v["dimension"], tuple(v["s"]), tuple(v["t"]), tuple(v["st"]))
                  for k, v in p1["signatures"].items()}
        smis = [r["label"] for r in rows if S.signature_of(r) != p1_sig.get(r["label"])]
        add("ROWS_SIGNATURES_MATCH_P1", not smis,
            smis or f"{len(rows)}/{S.ROW_COUNT} signatures match P1")

    census = {}
    for r in rows:
        census[r["evidentiary_class"]] = census.get(r["evidentiary_class"], 0) + 1
    add("CLASS_CENSUS", census == S.CLASS_CENSUS, census)

    # FIDELITY, not census. Swapping the declared class onto a free row, or
    # moving the selector class and updating the convention map to follow it,
    # preserved 1/7/1 and passed everything. The designated rows are named in the
    # PUBLIC protocol (§ 5.2, § 5.4) and their signatures come from P1, so the
    # binding is to a signature rather than to whichever row claims the class.
    if p1:
        want_cls = {}
        for lab, cls in ((S.DECLARED_LABEL, S.CLASS_DECLARED),
                         (S.SELECTOR_LABEL, S.CLASS_SELECTOR)):
            sg = p1["signatures"].get(lab)
            if sg:
                want_cls[(sg["dimension"], tuple(sg["s"]), tuple(sg["t"]),
                          tuple(sg["st"]))] = cls
        misassigned = []
        for r in rows:
            k = S.signature_of(r)
            expect = want_cls.get(k, S.CLASS_FREE)
            if r["evidentiary_class"] != expect:
                misassigned.append(expect)
        add("CLASS_ASSIGNMENT_MATCHES_P1", not misassigned,
            "declared and selector classes sit on the protocol-designated rows"
            if not misassigned else f"{len(misassigned)} row(s) carry the wrong class")

    # The subject is P1's DESIGNATED selector row, found by signature, not
    # whichever row claims the class. Picking the claimant meant a packet that
    # moved the class also moved which value this check examined, so the check
    # and the packet agreed with each other; CLASS_ASSIGNMENT_MATCHES_P1 caught
    # the move, but this check was individually satisfiable by the wrong row.
    p1_sel_sig = None
    if p1 and p1["signatures"].get(S.SELECTOR_LABEL):
        _sg = p1["signatures"][S.SELECTOR_LABEL]
        p1_sel_sig = (_sg["dimension"], tuple(_sg["s"]), tuple(_sg["t"]), tuple(_sg["st"]))
    sel = [r for r in rows if S.signature_of(r) == p1_sel_sig]
    # Q(phi) arithmetic on a triple with a zero denominator raises. Normalization
    # having FAILED is not a licence to keep computing: the check that noticed it
    # must gate the ones that would divide by it.
    if len(sel) == 1 and not norm_ok:
        add("CLASS_SELECTOR_NOT_SELF_INVERSE", False,
            "BLOCKED: triples are not normalized; exact arithmetic not attempted")
    elif len(sel) == 1:
        v = as_phi(sel[0]["value"])
        if v == Phi(0):
            # § 5.2's reciprocal class {x, x^-1} and § 5.4's "exactly one of
            # x = r and x^-1 = r" are both undefined at zero, which has no
            # inverse. Testing only v^2 = 1 passed a zero selector as "not
            # self-inverse", which is true and beside the point.
            add("CLASS_SELECTOR_NOT_SELF_INVERSE", False,
                "ZERO: the selector has no inverse, so its reciprocal class "
                "does not exist")
        else:
            si = (v * v).triple() == (1, 0, 1)
            add("CLASS_SELECTOR_NOT_SELF_INVERSE", not si,
                "selector value is not self-inverse" if not si else "SELF-INVERSE")
    else:
        # The detail must describe the SUBJECT THIS CHECK ACTUALLY USES. The
        # thirteenth audit moved that from "whichever row claims the selector
        # class" to "the row whose signature is P1's designated selector" and
        # left this sentence describing the old one, so a red verdict named a
        # condition the predicate no longer tests. Found by four lenses.
        add("CLASS_SELECTOR_NOT_SELF_INVERSE", False,
            f"{len(sel)} rows carry P1's designated selector signature, "
            f"expected exactly 1" if p1_sel_sig is not None
            else "P1 declares no signature for the designated selector row")

    # identities -----------------------------------------------------------
    idents = packet.get("identities", [])
    add("IDY_COUNT", len(idents) == S.IDENTITY_SLOTS, len(idents))

    positions = [i.get("position") for i in idents]
    # The FROZEN slot count, not the candidate's own length. Building the range
    # from `len(idents)` made this green for an identity array of any length,
    # the identity-side twin of the contiguity defect repaired at
    # IDXMAP_ZERO_BASED_CONTIGUOUS. Repairing one and not the other is what
    # "closed the instance, not the class" looks like from the inside.
    add("IDY_POSITIONS_CONTIGUOUS", positions == list(range(S.IDENTITY_SLOTS)),
        positions)

    by_sig = {}
    for r in rows:
        by_sig.setdefault(S.signature_of(r), []).append(r)
    unresolved = []
    for i in idents:
        for f in i.get("factors", []):
            key = S.signature_of(f)
            if len(by_sig.get(key, [])) != 1:
                unresolved.append((i.get("position"), len(by_sig.get(key, []))))
    add("IDY_FACTORS_RESOLVE_BY_SIGNATURE", not unresolved,
        "every factor resolves to exactly one row by signature"
        if not unresolved else f"{len(unresolved)} factors do not resolve uniquely")

    recomputed = []
    if not norm_ok:
        add("IDY_RECOMPUTES_FROM_ROWS", False,
            "BLOCKED: triples are not normalized; exact arithmetic not attempted")
    elif unresolved:
        # Previously this fell through to the summary line and reported an
        # ordinary red reading "0/0 recompute", as though four slots had been
        # compared and none matched. Nothing was compared. A check that did not
        # run says so.
        add("IDY_RECOMPUTES_FROM_ROWS", False,
            "BLOCKED: factors do not resolve; exact arithmetic not attempted")
    else:
        for i in idents:
            try:
                acc = Phi(1)
                for f in i["factors"]:
                    acc = acc * (as_phi(by_sig[S.signature_of(f)][0]["value"])
                                 ** f["exponent"])
                got = list(acc.triple())
            except ZeroDivisionError:
                # A ZERO row value is normalized: c > 0 and gcd(0, 0, 1) = 1, so
                # ROWS_TRIPLES_NORMALIZED passes it and `norm_ok` licenses
                # arithmetic it never validated. A negative exponent on such a
                # row then inverts zero. The gate was independent of this check
                # but not sufficient for it, and a raise is not a verdict.
                recomputed.append((i["position"], False))
                continue
            want = i.get("expected_value")
            # Compared elementwise on TRUE ints. Plain list equality accepts
            # True for 1, which is how a boolean leaf survived to here.
            same = (S.is_triple_shape(want)
                    and all(a == b and type(b) is int for a, b in zip(got, want)))
            recomputed.append((i["position"], same))
    if norm_ok and not unresolved:
        # Three clauses, and the detail used to report only the third: a red
        # verdict caused by the slot COUNT printed a sentence saying every slot
        # recomputed correctly. The eleventh audit claimed to have swept this
        # class; it had not, and this was the survivor.
        n_ok = sum(ok for _, ok in recomputed)
        rec_ok = (bool(recomputed) and len(recomputed) == S.IDENTITY_SLOTS
                  and all(ok for _, ok in recomputed))
        add("IDY_RECOMPUTES_FROM_ROWS", rec_ok,
            f"{n_ok}/{len(recomputed)} recompute from the committed rows, "
            "position-wise" if rec_ok else
            f"{n_ok}/{len(recomputed)} recompute, want {S.IDENTITY_SLOTS} of "
            f"{S.IDENTITY_SLOTS}")

    # Canonical factor order, declared rather than inherited from prose.
    order_bad = [i["position"] for i in idents
                 if i["factors"] != sorted(i["factors"], key=S.factor_sort_key)]
    add("IDY_FACTOR_ORDER_CANONICAL", not order_bad,
        "factors in frozen order" if not order_bad else f"slots {order_bad} unordered")

    # FIDELITY, not self-consistency. Every check above is satisfied by a packet
    # whose factors were swapped for other resolving signatures with the expected
    # value updated to match: internally coherent, and no longer P1's identity.
    if p1:
        # From P1's OWN signatures, not from the candidate's rows. Deriving the
        # expectation from the object under test let a coherent edit to a row
        # signature and every factor referencing it move both sides together, so
        # this check individually accepted a different identity object.
        sig_by_label = {lab: (sg["dimension"], tuple(sg["s"]), tuple(sg["t"]),
                              tuple(sg["st"]))
                        for lab, sg in (p1.get("signatures") or {}).items()}
        missing_rows = sorted({lab for i in p1.get("identities", [])
                               for lab, _ in i["factors"]} - set(sig_by_label))
        if missing_rows:
            add("IDY_DEFINITIONS_MATCH_P1", False,
                f"P1 identities reference rows absent from P1's signatures: {missing_rows}")
        else:
            # Compared as MULTISETS, position-wise across slots. Factor ORDER is
            # a separate frozen property with its own check; folding the two
            # together meant an ordering mutation tripped this one too and
            # neither could be tested alone.
            want_slots = [sorted((sig_by_label[lab], int(exp))
                                 for lab, exp in i["factors"])
                          for i in p1.get("identities", [])]
            got_slots = [sorted((S.signature_of(f), f["exponent"])
                                for f in i["factors"]) for i in idents]
            same = got_slots == want_slots
            add("IDY_DEFINITIONS_MATCH_P1", same,
                f"{len(want_slots)} slots match P1's definitions position-wise"
                if same else "an identity definition differs from P1's")

    # indexing map ---------------------------------------------------------
    im = packet.get("indexing_map", {})
    # Domains compared EXACTLY. Requiring only non-empty strings accepted any
    # domain name at all, so the map could describe a correspondence nobody
    # declared.
    # A detail must cover its OWN predicate. This tested three clauses and
    # printed two, so `zero_based: false` reddened while the diagnostic showed
    # two correct domain names and never named the failing clause.
    dom_bad = [n for n, ok in
               (("source_domain", im.get("source_domain") == S.INDEXING_SOURCE_DOMAIN),
                ("destination_domain",
                 im.get("destination_domain") == S.INDEXING_DESTINATION_DOMAIN),
                ("zero_based", im.get("zero_based") is True)) if not ok]
    add("IDXMAP_DOMAINS_DECLARED", not dom_bad,
        f"{im.get('source_domain')} -> {im.get('destination_domain')}" if not dom_bad
        else f"wrong: {dom_bad}; "
             f"{im.get('source_domain')!r} -> {im.get('destination_domain')!r}, "
             f"zero_based={im.get('zero_based')!r}")

    ent = im.get("entries", [])
    src = [e.get("source_position") for e in ent]
    dst = [e.get("destination_position") for e in ent]
    # Details are CONDITIONAL. These two stated the property as their detail on
    # failure as well as on success, so a red verdict printed a line asserting
    # the thing that had just failed, and the stops list contradicted itself.
    # It also made battery rule 3 structurally unavailable for them: a reason
    # substring cannot distinguish branches when both branches say the same
    # sentence, which is why their cases sat among the reason-less ones.
    # The detail names the FAILING domain. The predicate covers both, and an
    # earlier conditional detail printed the sources unconditionally, so a
    # destination-side break showed a healthy prefix of the domain that was
    # fine: a red verdict elaborated with innocent evidence. Misdirecting rather
    # than asserting, but the same family as the details that claimed the
    # property they had just failed.
    # The range is the FROZEN row count, not the candidate's. Reading it off the
    # object under test left an eight-row packet with eight consistent entries
    # green here, sixth instance of the class and the one indexing check the
    # earlier repair did not reach.
    want_range = list(range(S.ROW_COUNT))
    src_ok = sorted(src) == want_range
    dst_ok = sorted(dst) == want_range
    contig = src_ok and dst_ok
    add("IDXMAP_ZERO_BASED_CONTIGUOUS", contig,
        f"{len(ent)} entries over 0..{S.ROW_COUNT - 1}" if contig
        else f"not zero-based contiguous over 0..{S.ROW_COUNT - 1}, {len(ent)} entries; "
             + "; ".join(f"{nm} {sorted(v)[:5]}" for nm, v, ok in
                         (("sources", src, src_ok), ("destinations", dst, dst_ok))
                         if not ok))
    bij = len(set(src)) == len(src) == len(set(dst)) == len(dst)
    add("IDXMAP_BIJECTIVE", bij, "one-to-one and onto" if bij
        else f"not a bijection: {len(set(src))}/{len(src)} distinct sources, "
             f"{len(set(dst))}/{len(dst)} distinct destinations")
    # A map that is the identity cannot distinguish a harness that applies it
    # from one that ignores it. Stated precisely: § 4.4 requires a SYNTHETIC
    # nonidentity fixture to prove the map-processing path is operative, and
    # does not require every legitimate live map to be nonidentity. This
    # packet's is, because the canonical order and the record order disagree in
    # one transposition, so the live map does that work for free. The check is a
    # candidate-specific strengthening, not a § 4.4 schema requirement, and a
    # future geometry whose two orderings agreed would redden a correct packet.
    disc = src != dst and bool(ent)
    add("IDXMAP_DISCRIMINATING", disc, "non-identity permutation" if disc
        else ("no entries at all" if not ent else "IDENTITY: cannot discriminate"))

    # Non-identity is a PROPERTY, not correctness: exchanging two destinations
    # keeps it a valid contiguous non-identity bijection and changes what the map
    # says. The exact expected entries are derived from the canonical row
    # signatures and P1's declared record ordering.
    # DERIVED FROM P1 ALONE, both ends. Enumerating the candidate's rows for the
    # source positions and joining through the candidate's labels for the
    # destinations let the object under test help define its own expectation:
    # an eight-row packet whose entries dropped the same row passed, and
    # swapping the labels of two rows that share a value left this check
    # endorsing a correspondence that was simply wrong. The same raw lookup
    # `ref_pos[r["label"]]` also raised KeyError on any label P1 does not carry,
    # which is a verdict the predicate owes rather than an exception. Both are
    # gone: nothing the candidate says enters the expectation, so the check now
    # compares the candidate to P1 instead of to itself.
    if p1:
        ref_pos = {v["label"]: n for n, v in enumerate(p1["values"])}
        p1_rows = [{"label": lab,
                    "signature": {k: sg[k] for k in ("dimension", "s", "t", "st")}}
                   for lab, sg in p1["signatures"].items()]
        p1_rows.sort(key=S.row_sort_key)
        orphan = sorted(r["label"] for r in p1_rows if r["label"] not in ref_pos)
        if orphan:
            add("IDXMAP_MATCHES_P1", False,
                f"P1 signatures name rows absent from P1's values: {orphan}")
        else:
            # § 5.5 / p2_schema.INDEXING_ENTRY_ORDER: ascending source_position.
            want_entries = [{"source_position": n,
                             "destination_position": ref_pos[r["label"]]}
                            for n, r in enumerate(p1_rows)]
            add("IDXMAP_MATCHES_P1", ent == want_entries,
                "entries match the correspondence derived from P1"
                if ent == want_entries
                else "map differs from the P1-derived correspondence")

    # convention map -------------------------------------------------------
    cm = packet.get("convention_map", {})
    # FIVE now. § 4.4 asks this map for the § 5.4 bridge and anchor rule, the
    # § 4.2 basing reference, and the evaluation convention. The single field
    # once named `evaluation_convention` carried the Q(phi) encoding, which is
    # the encoding of the VALUES and is fixed by § 4.4's `rows` row anyway; the
    # convention § 4.2 defines, `g -> rho(g)`, was reachable only by following
    # the basing reference, which is a different item on the same list.
    want_cm = {"bridge", "orientation_anchor_rule", "basing_reference",
               "representation_evaluation"}
    add("CONV_MAP_COMPLETE", set(cm) == want_cm,
        sorted(cm) if set(cm) != want_cm
        else "bridge, anchor, basing, representation evaluation")
    # Derived from P1's designated selector row, NOT from whichever row claims
    # the class. Taking it from the claimant let the map follow a reassigned
    # selector and agree with itself.
    p1_sel = p1["signatures"].get(S.SELECTOR_LABEL) if p1 else None
    sel_sig = ({"dimension": p1_sel["dimension"], "s": list(p1_sel["s"]),
                "t": list(p1_sel["t"]), "st": list(p1_sel["st"])}
               if p1_sel else None)
    # The construction hash comes from P1, NOT from the candidate's own
    # `adjudicates`. While the template's one dynamic field was read off the
    # object under test, moving a fabricated hash into `adjudicates` and into
    # `basing_reference` together left this check and
    # CONV_BASING_MATCHES_ADJUDICATES both green: the two agreed with each other
    # and with nothing else, and the ensemble was one unrecorded dependency on
    # ADJUDICATES_MATCH_P1 deep.
    # The § 4.2 convention comes from the CONSTRUCTION PACKET, hash-pinned,
    # not from a literal here. Until P5 read the protocol rather than the packet,
    # the template was both what the emitter built from and what this compared
    # against, so P2's reading of § 4.4 could not fail any check in the build.
    try:
        want_repr_eval = S.representation_evaluation_from_construction(
            p1["adjudicates"]["construction_packet"])
    except Exception as e:                                    # noqa: BLE001
        want_repr_eval = None
        repr_eval_why = f"{type(e).__name__}: {e}"
    want_cm_full = (S.convention_map_template(
        p1["adjudicates"]["construction_packet"],
        S.SELECTOR_LABEL, sel_sig, want_repr_eval)
        if sel_sig and want_repr_eval is not None else None)
    add("CONV_MAP_EXACT", want_cm_full is not None and cm == want_cm_full,
        "convention map matches the frozen template exactly" if cm == want_cm_full
        else "convention map content differs from the frozen template")

    # Named on its own so the property can redden by itself rather than only as
    # part of the whole-map comparison. § 4.4 requires this map to CARRY § 4.2's
    # evaluation convention; reaching it through `basing_reference` is a
    # different item on § 4.4's list and does not discharge this one.
    re_block = cm.get("representation_evaluation") if isinstance(cm, dict) else None
    declared_re = re_block.get("convention") if isinstance(re_block, dict) else None
    re_ok = want_repr_eval is not None and declared_re == want_repr_eval
    add("CONV_REPRESENTATION_EVALUATION_MATCHES_4_2", re_ok,
        "the declared representation-evaluation convention is the construction "
        "packet's own `basing.evaluation`, read from the pinned bytes"
        if re_ok else
        (f"declared {declared_re!r}, construction packet says {want_repr_eval!r}"
         if want_repr_eval is not None else
         f"could not read the construction packet: {repr_eval_why}"))

    br = cm.get("basing_reference")
    br_sha = br.get("construction_packet_sha256") if isinstance(br, dict) else None
    basing_ok = br_sha == adj.get("construction_packet_sha256")
    add("CONV_BASING_MATCHES_ADJUDICATES", basing_ok,
        "basing reference pins the same construction packet as `adjudicates`"
        if basing_ok else
        f"basing reference pins a DIFFERENT construction packet: "
        f"{br_sha!r} vs adjudicates "
        f"{adj.get('construction_packet_sha256')!r}")

    # labels must not be load-bearing --------------------------------------
    # Strip every label, then require the packet to remain fully resolvable.
    #
    # Stated exactly, because the removal of the unreachable clause made the
    # stripping itself inert: `_resolves_without_labels` no longer reads a label
    # at all, so the strip changes nothing it looks at. It is kept because the
    # invariant it expresses is the one that matters, and because a future
    # clause that DOES read a label must find them already gone. The comment
    # used to say "if anything joins by label, this is where it shows", which
    # after that removal was a claim about a mechanism no longer present.
    stripped = json.loads(json.dumps(packet))
    for r in stripped["rows"]:
        r.pop("label", None)
    lab_ok, lab_detail = _resolves_without_labels(stripped)
    add("LABELS_ARE_NOT_JOIN_KEYS", lab_ok, lab_detail)

    # bytes and hash record -------------------------------------------------
    # Validate the BYTES, not a round trip of an object we just serialized.
    # Re-canonicalizing our own output is a fixed point by construction and could
    # not fail; this inspects the actual bytes against § 10.2 and is reddened by
    # supplying non-canonical bytes.
    # The object and the bytes must be the same artifact, compared as BYTES.
    # Decoded-object equality was defeatable because Python says 1 == 1.0: a
    # packet with an integer could be paired with canonical bytes carrying a
    # float lexeme, the semantic float check read the object and saw none, and
    # the sidecar hashed the alternate bytes correctly.
    try:
        bound = pbytes == S.canonical_bytes(packet)
        bind_why = ("bytes are exactly the canonical serialization of the assessed "
                    "object" if bound else
                    "supplied bytes are not exactly the canonical serialization "
                    "of the supplied packet")
    except Exception as e:                                    # noqa: BLE001
        bound, bind_why = False, f"packet not canonically serializable: {type(e).__name__}"
    add("BYTES_BIND_TO_PACKET", bound, bind_why)

    # ONE emission covering both artifacts. The record half used to be a fixed
    # point by construction and this comment said so; adding `parse_constant`
    # made that FALSE. `json.loads` admits NaN, Infinity and -Infinity, and
    # `json.dumps` re-emits them, so a sidecar FILE carrying one of those
    # lexemes canonicalizes to bytes this function now refuses: a live,
    # JSON-reachable failure branch. It was previously a second `add()` of this
    # same id, which would have emitted a duplicate check id had it ever fired.
    ok_bytes, why = _canonical_conformance(pbytes)
    rec_ok, rec_conf_why = _canonical_conformance(data["hash_record_bytes"])
    conf_bad = [f"packet: {why}"] if not ok_bytes else []
    conf_bad += [f"hash record: {rec_conf_why}"] if not rec_ok else []
    add("CANONICAL_BYTES_FIXED_POINT", ok_bytes and rec_ok,
        f"{len(pbytes)} bytes, § 10.2 conformant; hash record conformant"
        if not conf_bad else "; ".join(conf_bad))
    hr_shape = S.exact_shape(hash_record, S.HASH_RECORD_KEYS, "hash_record")
    add("HASH_RECORD_SCHEMA_EXACT", not hr_shape,
        "record has exactly its declared fields" if not hr_shape else hr_shape[:2])

    # Compared as BYTES, like the packet. Dict equality is not type-sensitive,
    # so a float byte length passed a check whose name is EXACT: Python says
    # 10974.0 == 10974. Same class as the packet's 1 vs 1.0, at the sidecar.
    # Guarded the way BYTES_BIND_TO_PACKET is, so the repair cannot introduce
    # the raise it was written to remove.
    want_rec = S.hash_record_template(sha256_bytes(pbytes), len(pbytes))
    try:
        rec_exact = S.canonical_bytes(hash_record) == S.canonical_bytes(want_rec)
        rec_why = ("sidecar matches the frozen template exactly, byte for byte"
                   if rec_exact else "sidecar content differs from the frozen template")
    except Exception as e:                                    # noqa: BLE001
        rec_exact = False
        rec_why = f"sidecar not canonically serializable: {type(e).__name__}"
    add("HASH_RECORD_EXACT", rec_exact, rec_why)

    # The DETAIL, not the predicate, was the raise. `.get(k, "")[:16]` defends
    # against a MISSING key and not against a present one of the wrong type, so
    # a sidecar carrying an integer, null, boolean, float or object at the hash
    # field sliced a non-string and escaped assess(). The comparison above was
    # always safe, since `==` across types is merely False.
    #
    # Fixed HERE rather than by gating these checks on HASH_RECORD_SCHEMA_EXACT,
    # which is the shape this build usually reaches for. That gate would have
    # been fix-induced damage: `hash_record_missing_core_fields` requires
    # HASH_RECORD_MATCHES_PACKET and HASH_RECORD_IS_CANDIDATE to stay GREEN
    # while the schema check reddens, which is the whole point of that case, and
    # blocking them would have turned a passing battery case into a MISS.
    got_sha = hash_record.get("expected_canonical_plaintext_sha256")
    got_len = hash_record.get("canonical_plaintext_byte_length")
    rec_bad = [n for n, ok in (("sha256", got_sha == sha256_bytes(pbytes)),
                               ("byte length", got_len == len(pbytes))) if not ok]
    add("HASH_RECORD_MATCHES_PACKET", not rec_bad,
        (got_sha[:16] if isinstance(got_sha, str)
         else f"{type(got_sha).__name__}: {got_sha!r}") if not rec_bad
        else f"wrong: {rec_bad}; record says {got_sha!r} at {got_len!r} bytes, "
             f"packet is {sha256_bytes(pbytes)[:16]}… at {len(pbytes)} bytes")
    add("HASH_RECORD_IS_CANDIDATE", hash_record.get("status") == "CANDIDATE",
        hash_record.get("status"))

    reg, data["registry"] = validate_registry([c for c, _, _ in checks])
    stops += reg
    data["environment"] = environment()
    return checks, stops, data


def _canonical_conformance(b):
    """§ 10.2: JSON, keys sorted, two-space indent, ASCII, LF, one trailing
    newline. Checked against the bytes themselves."""
    if b"\r" in b:
        return False, "contains CR; § 10.2 requires LF"
    try:
        text = b.decode("ascii")
    except UnicodeDecodeError:
        return False, "not ASCII"
    if not text.endswith("\n") or text.endswith("\n\n"):
        return False, "not exactly one trailing newline"
    try:
        # `parse_constant` fires on NaN, Infinity and -Infinity, which Python's
        # parser accepts and JSON does not. Without it, bytes carrying those
        # lexemes round-tripped and were certified § 10.2 conformant.
        obj = json.loads(text, parse_constant=_reject_constant)
    except Exception as e:                                    # noqa: BLE001
        return False, f"not parseable JSON: {e}"
    try:
        canon = S.canonical_bytes(obj)
    except Exception as e:                                    # noqa: BLE001
        # The two structurally identical canonicalizations beside this one are
        # guarded; this one was not.
        return False, f"content not canonically serializable: {type(e).__name__}"
    if canon != b:
        return False, "bytes differ from the canonical serialization of their own content"
    return True, "conformant"


def _reject_constant(name):
    raise ValueError(f"{name} is not JSON; § 10.2 admits no such literal")


def _find_floats(obj, path="$"):
    """ITERATIVE, and holding only a SEGMENT per stacked node.

    Two defects, one per revision. It began recursive, so a packet nested deeper
    than Python's frame limit escaped with RecursionError; `json.loads` accepts
    roughly ten times the depth the interpreter can recurse over, and this is
    called ahead of every gate and outside the exception boundary, so the run
    returned no registry at all.

    The iterative repair then held a fully materialised path string for every
    node on the stack at once. Depth-first recursion had freed each branch's
    paths as it unwound; a stack holds every sibling's simultaneously, so peak
    memory went from O(depth^2) to O(nodes x depth) and a wide, deep packet
    exhausted memory instead. Each entry now carries only its own short segment
    plus a parent index, and a path is assembled ONLY for a node that is
    actually a float, which is the rare case. Emission order is preserved by
    pushing children reversed."""
    out, stack, trail, built = [], [(obj, None, None)], [], {}
    while stack:
        o, parent, seg = stack.pop()
        me = len(trail)
        trail.append((parent, seg))
        if isinstance(o, float):
            # Memoised. The previous revision walked the parent chain per
            # float, which moved the O(depth) cost off every node and onto
            # every float rather than removing it: a deep structure holding many
            # floats was still quadratic. Each node's path is now built once,
            # from its parent's already-built path, and only along ancestries
            # that actually lead to a float.
            if me not in built:
                chain = []
                node = me
                while node is not None and node not in built:
                    chain.append(node)
                    node = trail[node][0]
                base = built[node] if node is not None else path
                for n in reversed(chain):
                    seg = trail[n][1]
                    base = base + seg if seg is not None else base
                    built[n] = base
            out.append(built[me])
        elif isinstance(o, dict):
            stack.extend(reversed([(v, me, f".{k}") for k, v in o.items()]))
        elif isinstance(o, list):
            stack.extend(reversed([(v, me, f"[{i}]") for i, v in enumerate(o)]))
    return out


def _resolves_without_labels(packet):
    """Labels must never be needed to resolve a factor to its row.

    ONE requirement now, stated as one. This once had two: a stray-label scan and
    a resolution invariant. The scan was subsumed when the nested schema became
    exact, and the residual `label in row` test is unreachable because the only
    caller strips every row label before calling. The docstring and the success
    detail both went on claiming two requirements after only one remained."""
    # The stray-label scan that used to live here is now subsumed: a label on a
    # factor or an indexing entry is an undeclared field, and NESTED_SCHEMA_EXACT
    # rejects it first. What remains is the resolution invariant.
    by_sig = {}
    for r in packet.get("rows", []):
        by_sig.setdefault(S.signature_of(r), []).append(r)
    for i in packet.get("identities", []):
        for f in i.get("factors", []):
            if len(by_sig.get(S.signature_of(f), [])) != 1:
                return False, "an identity factor cannot be resolved without labels"
    # The `label in row` test that used to sit here was unreachable: the caller
    # pops every row label before calling, so it could never fire. An
    # unreachable clause reported as a requirement is a claim the code does not
    # keep, which is the class this build has been closing in documentation.
    return True, "labels are never needed to resolve a factor to its row"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--p1-artifact", required=True)
    ap.add_argument("--emit", metavar="DIR")
    a = ap.parse_args()

    if a.emit:
        prior = pathlib.Path(a.emit)
        for name in ("m8_8_answer_packet.CANDIDATE.json",
                     "m8_8_prelock_hash_record.CANDIDATE.json",
                     "p2_runtime.log.json"):
            f = prior / name
            if f.exists():
                f.rename(f.with_suffix(f.suffix + ".superseded"))

    checks, stops, data = assess(a.p1_artifact)

    print("P2 candidate emission\n")
    w = max(len(c) for c, _, _ in checks)
    for cid, ok, detail in checks:
        print(f"  {'ok  ' if ok else 'FAIL'} {cid:{w}s}  {detail}")
    env = data["environment"]
    print(f"\n  environment: python {env['python']}, numpy {env['numpy']}, "
          f"mpmath {env['mpmath']}")

    print()
    if stops:
        print("  STOP. No candidate emitted.")
        for s in stops:
            print(f"    - {s}")
        if a.emit:
            e = pathlib.Path(a.emit)
            e.mkdir(parents=True, exist_ok=True)
            for f in e.glob("*.superseded"):
                f.unlink()
            (e / "P2_FAILED.txt").write_text(
                "P2 STOPPED. No candidate was produced, and any candidate from an\n"
                "earlier run has been removed rather than left to be mistaken for\n"
                "current output.\n\n" + "\n".join(stops) + "\n")
        return 1
    print("  P2 GREEN on all declared checks. CANDIDATE, not issued.")

    if a.emit:
        out = pathlib.Path(a.emit)
        out.mkdir(parents=True, exist_ok=True)
        pbytes = data["packet_bytes"]

        pk = out / "m8_8_answer_packet.CANDIDATE.json"
        tmp = pk.with_suffix(".partial")
        tmp.write_bytes(pbytes)
        tmp.replace(pk)

        # EXACTLY the assessed object. main() used to copy the record, add the
        # banner, and serialize that changed object, so the bytes on disk were
        # not the bytes assess() validated. The banner is now part of the record
        # as built, and nothing is modified after the verdict.
        rec = out / "m8_8_prelock_hash_record.CANDIDATE.json"
        tmp = rec.with_suffix(".partial")
        tmp.write_bytes(data["hash_record_bytes"])
        tmp.replace(rec)

        (out / "p2_runtime.log.json").write_text(json.dumps({
            "_note": "NONCANONICAL. Machine-local detail; not part of any artifact.",
            "p1_artifact_path": str(pathlib.Path(a.p1_artifact).resolve()),
            "environment": data["environment"],
        }, indent=2, sort_keys=True) + "\n")

        for f in out.glob("*.superseded"):
            f.unlink()
        (out / "P2_FAILED.txt").unlink(missing_ok=True)
        print(f"\n  wrote {pk}")
        print(f"  wrote {rec}")
        print(f"  candidate packet sha256: {sha256_bytes(pbytes)}")
        print(f"  canonical byte length:   {len(pbytes)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
