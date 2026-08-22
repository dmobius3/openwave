"""P3: the mechanical schema gate.

ONE question. Are these exact bytes a valid seven-key answer packet and pre-lock
hash record under the frozen schema and the § 10.2 canonical form?

INDEPENDENCE IS THE POINT, and it is a hard rule rather than a preference:

    P3 may consume the P2 packet, the sidecar, the pinned P1 artifact and the
    governing protocol sections. It may NOT import or execute P2 emitter,
    schema or checker code.

Every expectation below is encoded here from the protocol text, or derived from
the pinned P1 artifact and the two public packets, and never read out of
`p2_schema`. A gate that validated the packet against the emitter's own frozen
constants would establish only that P2 agreed with its own declaration, which P2
established already and which cannot fail. Arithmetic comes from `p3_qphi`,
written independently for the same reason.

TWO MODES, and the distinction matters more than it looks:

    production   the three input hashes must match, and then every check runs
    injected     a candidate is supplied directly; the hash pin does not apply,
                 and the checks run against the same independent expectations

The pin is an input-integrity PRECONDITION. It is never evidence that the schema
is correct, and it is reported separately so it cannot be mistaken for one. If
the pin were the only way in, a hash mismatch would stop everything and no
individual check could ever be exercised, which is precisely what P4 has to do.
This is the same separation P1 and P2 each had to make, between trusting the
input version and testing the acceptance predicate.

DELIBERATELY NARROW. P3 checks what the protocol requires of the packet and
nothing else. P2 carries 38 checks, many of them defensive additions paid for by
its own audit history; importing that breadth here would trade P3's independence
for coverage it was not asked for. The plan's checklist is the scope.

`p3_selftest.py` is an ENGINEERING self-test, not protocol evidence. The formal
claim that every P3 check can be reddened belongs to P4's content-mutation half,
where the frozen plan puts it.
"""

import argparse
import hashlib
import json
import pathlib
import sys

import p3_forms as FORMS
import p3_qphi as Q

# ---------------------------------------------------------------- expectations
# Encoded HERE from the protocol text. Not imported.

# § 4.4, the answer packet's key table: exactly these seven.
PACKET_KEYS = ("adjudicates", "convention_map", "format_version", "identities",
               "indexing_map", "rows", "target_id")

# § 5.2, the evidentiary classes and their census: one declared, seven free, one
# free orientation selector.
CLASS_DECLARED = "declared_convention"
CLASS_FREE = "free"
CLASS_SELECTOR = "free_orientation_selector"
CLASS_CENSUS = {CLASS_DECLARED: 1, CLASS_FREE: 7, CLASS_SELECTOR: 1}

# § 5.2 names both designated rows publicly, and § 5.4 names the anchor.
DECLARED_ROW = "R0"
SELECTOR_ROW = "R7"

# § 5.2: the declared ledger convention is T^2(R0) = 1.
DECLARED_ROW_VALUE = [1, 0, 1]

# § 5.5: nine irreps.
ROW_COUNT = 9
# § 5.3: four identity slots.
IDENTITY_SLOTS = 4
# § 5.5: the row signature is dimension plus exact characters at s, t and st.
SIGNATURE_FIELDS = ("dimension", "s", "t", "st")

# The GROUP packet hash, § 11's pin, verbatim and typed here because § 11 is
# where it lives.
#
# The CONSTRUCTION packet hash is deliberately NOT typed here. § 11 still reads
# "[PIN at landing]" for it, so there is no protocol pin to quote, and adding a
# second hand-typed literal would create a second trust root for P3 to be wrong
# about. It is read instead from the P1 artifact P3 already pins by hash, which
# makes the chain explicit and checkable end to end:
#
#     P3's own P1 hash pin  ->  P1.adjudicates.construction_packet
#                           ->  the supplied construction packet bytes
#
# Both ends are enforced below: the candidate must name the same hash, and the
# file handed to P3 must actually hash to it.
PROTOCOL_GROUP_PIN = \
    "e3b0c945bbbb15b4549fa641234c9461062c2337b3d1e372af621b614d4883a9"

CHECK_IDS = [
    "KEYS_EXACTLY_SEVEN",
    "PACKET_CANONICAL_BYTES",
    "SIDECAR_CANONICAL_BYTES",
    "SIDECAR_COMMITS_TO_PACKET",
    "TRIPLES_NORMALIZED",
    "CLASS_CENSUS",
    "ADJUDICATES_MATCH_PUBLIC_PACKETS",
    "SIGNATURES_DISTINCT",
    "DECLARED_ROW_CONVENTION",
    "IDENTITIES_RECOMPUTE_POSITION_WISE",
    "VALUES_MATCH_BOTH_SOURCE_FAMILIES",
    "SELECTOR_NOT_SELF_INVERSE",
]


def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()


def canonical_bytes(obj):
    """§ 10.2: JSON, keys sorted, two-space indent, ASCII, LF, one trailing
    newline. P3's own implementation, so a defect in P2's would be visible."""
    text = json.dumps(obj, sort_keys=True, indent=2, ensure_ascii=True)
    return (text + "\n").replace("\r\n", "\n").encode("ascii")


def canonical_problems(b, what):
    """§ 10.2 conformance, checked against the BYTES rather than by round trip."""
    if not isinstance(b, bytes):
        return [f"{what}: not bytes"]
    if b"\r" in b:
        return [f"{what}: contains CR, § 10.2 requires LF"]
    try:
        text = b.decode("ascii")
    except UnicodeDecodeError:
        return [f"{what}: not ASCII"]
    if not text.endswith("\n") or text.endswith("\n\n"):
        return [f"{what}: not exactly one trailing newline"]
    try:
        obj = json.loads(text, parse_constant=_not_json)
    except Exception as e:                                    # noqa: BLE001
        return [f"{what}: not parseable JSON: {e}"]
    try:
        if canonical_bytes(obj) != b:
            return [f"{what}: not the canonical serialization of its own content"]
    except Exception as e:                                    # noqa: BLE001
        return [f"{what}: content not canonically serializable: {type(e).__name__}"]
    return []


def _not_json(name):
    raise ValueError(f"{name} is not JSON; § 10.2 admits no such literal")


def signature_of(row):
    """§ 5.5's label-free identity, as a hashable tuple."""
    s = row["signature"]
    return (s["dimension"], tuple(s["s"]), tuple(s["t"]), tuple(s["st"]))


def validate_registry(emitted):
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


def gate(packet_bytes, sidecar_bytes, p1, group_bytes, construction_bytes):
    """The acceptance predicate. Returns (checks, stops).

    Takes BYTES for both artifacts, because § 10.2 conformance is a property of
    bytes and validating a re-serialization of an object we just parsed would be
    a fixed point. Everything else is parsed from those same bytes, so the object
    checked and the bytes checked cannot diverge.
    """
    checks, stops = [], []

    def add(cid, ok, detail):
        checks.append((cid, ok, detail))
        if not ok:
            stops.append(f"{cid}: {detail}")
        return ok

    def blocked(cids, why):
        for c in cids:
            add(c, False, f"BLOCKED: {why}")

    # § 10.2, both artifacts, before anything is parsed for meaning ----------
    pk_bad = canonical_problems(packet_bytes, "packet")
    add("PACKET_CANONICAL_BYTES", not pk_bad,
        f"{len(packet_bytes)} bytes, § 10.2 conformant" if not pk_bad else pk_bad[0])
    sc_bad = canonical_problems(sidecar_bytes, "sidecar")
    add("SIDECAR_CANONICAL_BYTES", not sc_bad,
        f"{len(sidecar_bytes)} bytes, § 10.2 conformant" if not sc_bad else sc_bad[0])

    remaining = [c for c in CHECK_IDS
                 if c not in ("PACKET_CANONICAL_BYTES", "SIDECAR_CANONICAL_BYTES")]
    if pk_bad or sc_bad:
        blocked(remaining, "artifact bytes are not § 10.2 conformant")
        s, reg = validate_registry([c for c, _, _ in checks])
        return checks, stops + s

    packet = json.loads(packet_bytes.decode("ascii"))
    sidecar = json.loads(sidecar_bytes.decode("ascii"))

    # the sidecar commits to THESE bytes --------------------------------------
    want_sha, want_len = sha256_bytes(packet_bytes), len(packet_bytes)
    got_sha = sidecar.get("expected_canonical_plaintext_sha256")
    got_len = sidecar.get("canonical_plaintext_byte_length")
    bad = [n for n, ok in (("sha256", got_sha == want_sha),
                           ("byte length", got_len == want_len and type(got_len) is int))
           if not ok]
    add("SIDECAR_COMMITS_TO_PACKET", not bad,
        f"{want_sha[:16]} at {want_len} bytes" if not bad
        else f"wrong: {bad}; record says {got_sha!r} at {got_len!r}")

    # § 4.4, exactly seven keys ----------------------------------------------
    keys = tuple(sorted(packet)) if isinstance(packet, dict) else ()
    keys_ok = add("KEYS_EXACTLY_SEVEN", keys == PACKET_KEYS,
                  f"{len(keys)} keys as declared" if keys == PACKET_KEYS
                  else f"got {keys}, want {PACKET_KEYS}")
    if not keys_ok:
        blocked([c for c in remaining
                 if c not in ("KEYS_EXACTLY_SEVEN", "SIDECAR_COMMITS_TO_PACKET")],
                "top-level key set is not the declared seven")
        s, reg = validate_registry([c for c, _, _ in checks])
        return checks, stops + s

    rows = packet["rows"]
    idents = packet["identities"]

    # § 4.4's encoding: every triple exact and normalized ---------------------
    trips = []
    shape_bad = []
    if isinstance(rows, list):
        for n, r in enumerate(rows):
            if not isinstance(r, dict) or not isinstance(r.get("signature"), dict):
                shape_bad.append(f"rows[{n}]: not a row object")
                continue
            trips.append((f"rows[{n}].value", r.get("value")))
            for f in ("s", "t", "st"):
                trips.append((f"rows[{n}].signature.{f}", r["signature"].get(f)))
            if type(r["signature"].get("dimension")) is not int:
                shape_bad.append(f"rows[{n}].signature.dimension: not a true int")
    if isinstance(idents, list):
        for n, i in enumerate(idents):
            if not isinstance(i, dict):
                shape_bad.append(f"identities[{n}]: not an object")
                continue
            trips.append((f"identities[{n}].expected_value", i.get("expected_value")))
    unnorm = [w for w, t in trips if not Q.is_normalized(t)]
    norm_ok = add("TRIPLES_NORMALIZED", not unnorm and not shape_bad,
                  f"{len(trips)} triples exact and normalized"
                  if not unnorm and not shape_bad
                  else f"{len(unnorm) + len(shape_bad)} problem(s): "
                       f"{(unnorm + shape_bad)[:2]}")

    # § 5.2's census ----------------------------------------------------------
    census = {}
    for r in rows if isinstance(rows, list) else []:
        if isinstance(r, dict):
            census[r.get("evidentiary_class")] = census.get(r.get("evidentiary_class"), 0) + 1
    add("CLASS_CENSUS", census == CLASS_CENSUS,
        f"{CLASS_CENSUS} as § 5.2 declares" if census == CLASS_CENSUS else str(census))

    # `adjudicates` against the PUBLIC PACKETS, hashed here -------------------
    # Not against P1's record of them and not against P2's: P3 hashes the two
    # public files itself, which is the only form of this check that is
    # independent of both earlier steps.
    adj = packet["adjudicates"]
    p1_construction = (p1.get("adjudicates") or {}).get("construction_packet")
    why = []
    if adj.get("group_packet_sha256") != PROTOCOL_GROUP_PIN:
        why.append("`adjudicates` group hash is not § 11's pin")
    if sha256_bytes(group_bytes) != PROTOCOL_GROUP_PIN:
        why.append("the supplied group packet is not the § 11 pinned object")
    if not p1_construction:
        why.append("the pinned P1 artifact declares no construction packet hash")
    else:
        if adj.get("construction_packet_sha256") != p1_construction:
            why.append("`adjudicates` construction hash is not the one the pinned "
                       "P1 adjudication record declares")
        if sha256_bytes(construction_bytes) != p1_construction:
            why.append("the supplied construction packet does not hash to the one "
                       "the pinned P1 adjudication record declares")
    add("ADJUDICATES_MATCH_PUBLIC_PACKETS", not why,
        "group hash is § 11's pin; construction hash is the pinned P1's "
        "declaration; both supplied files hash to them" if not why
        else "; ".join(why))

    # § 5.5: the signature must separate all nine ----------------------------
    try:
        sigs = [signature_of(r) for r in rows]
        add("SIGNATURES_DISTINCT",
            len(set(sigs)) == len(sigs) == ROW_COUNT,
            f"{len(set(sigs))} distinct of {len(sigs)}"
            + ("" if len(sigs) == ROW_COUNT else f", want {ROW_COUNT}"))
    except Exception as e:                                    # noqa: BLE001
        add("SIGNATURES_DISTINCT", False, f"unreadable signatures: {type(e).__name__}")

    by_label = {r.get("label"): r for r in rows if isinstance(r, dict)}

    # § 5.2: R0 carries the declared convention, value 1 ---------------------
    r0 = by_label.get(DECLARED_ROW)
    if r0 is None:
        add("DECLARED_ROW_CONVENTION", False, f"no row labelled {DECLARED_ROW}")
    else:
        why = []
        if r0.get("value") != DECLARED_ROW_VALUE:
            why.append(f"value {r0.get('value')}, want {DECLARED_ROW_VALUE}")
        if r0.get("evidentiary_class") != CLASS_DECLARED:
            why.append(f"class {r0.get('evidentiary_class')!r}, want {CLASS_DECLARED!r}")
        add("DECLARED_ROW_CONVENTION", not why,
            f"{DECLARED_ROW} declares T^2 = 1 as § 5.2 requires" if not why
            else "; ".join(why))

    # § 5.3: recompute every identity from the committed rows, POSITION-WISE --
    if not norm_ok:
        add("IDENTITIES_RECOMPUTE_POSITION_WISE", False,
            "BLOCKED: triples are not normalized; exact arithmetic not attempted")
    else:
        by_sig = {}
        for r in rows:
            by_sig.setdefault(signature_of(r), []).append(r)
        results, why = [], None
        try:
            for slot, i in enumerate(idents):
                if i.get("position") != slot:
                    why = (f"slot {slot} declares position {i.get('position')!r}; "
                           "§ 5.3 compares position-wise")
                    break
                acc = Q.ONE
                for f in i["factors"]:
                    hits = by_sig.get(signature_of(f), [])
                    if len(hits) != 1:
                        raise KeyError(f"factor resolves to {len(hits)} rows")
                    acc = Q.mul(acc, Q.power(tuple(hits[0]["value"]), f["exponent"]))
                results.append(list(acc) == i.get("expected_value"))
        except Exception as e:                                # noqa: BLE001
            why = f"{type(e).__name__}: {e}"
        ok = (why is None and len(results) == IDENTITY_SLOTS and all(results))
        add("IDENTITIES_RECOMPUTE_POSITION_WISE", ok,
            f"{IDENTITY_SLOTS}/{IDENTITY_SLOTS} recompute from the committed rows"
            if ok else (why or f"{sum(results)}/{len(results)} recompute, "
                               f"want {IDENTITY_SLOTS} of {IDENTITY_SLOTS}"))

    # BOTH SOURCE FAMILIES, established by evaluating what each source WROTE.
    #
    # An earlier version of this check compared the candidate to P1's adjudicated
    # triple and then confirmed P1 had recorded both families and had marked its
    # own agreement gates green. That is P1's assertion re-reported, not an
    # independent check, and it left the substantive half of the plan's claim
    # resting on P1's verdict bits.
    #
    # P3 now parses each stated closed form with its own grammar, evaluates it
    # with its own exact arithmetic, and requires the result to equal the value
    # the packet commits. That is source verification, not rederivation: no
    # program is run and no spectral quantity is recomputed. The family of a
    # rendering is derived from the source key rather than taken from P1's
    # `families` field, and the two are cross-checked, so neither can be trusted
    # alone. An unparseable form is a STOP, not a skip.
    try:
        p1v = {v["label"]: v for v in p1["values"]}
        bad_rows = []
        for r in rows:
            lab, val = r.get("label"), r.get("value")
            rec = p1v.get(lab)
            if rec is None:
                bad_rows.append(f"{lab!r}: absent from P1")
                continue
            if list(rec["triple"]) != val:
                bad_rows.append(f"{lab}: packet value differs from P1's")
            fams = set()
            for src, form in sorted((rec.get("as_written") or {}).items()):
                try:
                    got = list(FORMS.evaluate(form))
                except Exception as e:                        # noqa: BLE001
                    bad_rows.append(f"{lab}/{src.split()[0]}: unparseable stated "
                                    f"form ({e})")
                    continue
                if got != val:
                    bad_rows.append(f"{lab}/{src.split()[0]}: stated form "
                                    f"evaluates to {got}, packet commits {val}")
                parts = src.split()
                fams.add(parts[1].split("/")[0] if len(parts) > 1 else "?")
            if not {"theory", "platform"} <= fams:
                bad_rows.append(f"{lab}: renderings cover {sorted(fams)}, "
                                "want both families")
            if fams != set(rec.get("families") or []):
                bad_rows.append(f"{lab}: P1 lists families {rec.get('families')}, "
                                f"source keys say {sorted(fams)}")
        add("VALUES_MATCH_BOTH_SOURCE_FAMILIES", not bad_rows,
            f"{len(rows)} rows: every stated form in both families evaluates, "
            "in P3's own arithmetic, to the value the packet commits"
            if not bad_rows else f"{len(bad_rows)}: {bad_rows[:2]}")
    except Exception as e:                                    # noqa: BLE001
        add("VALUES_MATCH_BOTH_SOURCE_FAMILIES", False,
            f"P1 artifact not in the expected shape: {type(e).__name__}: {e}")

    # § 5.4: the anchor must be able to discriminate -------------------------
    # Binds the CLASS as well as the value, mirroring the R0 check. Checking only
    # the value let the selector class be moved off R7 onto a free row with the
    # 1/7/1 census preserved and every P3 check still green: the census counted
    # the classes and nothing tied the selector class to the protocol-designated
    # row. That exact subject-ownership defect occurred in P2 as well.
    r7 = by_label.get(SELECTOR_ROW)
    if r7 is None:
        add("SELECTOR_NOT_SELF_INVERSE", False, f"no row labelled {SELECTOR_ROW}")
    elif not norm_ok:
        add("SELECTOR_NOT_SELF_INVERSE", False,
            "BLOCKED: triples are not normalized; exact arithmetic not attempted")
    else:
        try:
            why7 = []
            if r7.get("evidentiary_class") != CLASS_SELECTOR:
                why7.append(f"{SELECTOR_ROW} carries "
                            f"{r7.get('evidentiary_class')!r}, want {CLASS_SELECTOR!r}")
            v = tuple(r7["value"])
            if Q.mul(v, v) == Q.ONE:
                why7.append("SELF-INVERSE: cannot select an orientation")
            add("SELECTOR_NOT_SELF_INVERSE", not why7,
                f"{SELECTOR_ROW} carries the selector class and is not "
                "self-inverse, so it can discriminate" if not why7
                else "; ".join(why7))
        except Exception as e:                                # noqa: BLE001
            add("SELECTOR_NOT_SELF_INVERSE", False,
                f"selector value unusable: {type(e).__name__}: {e}")

    s, _reg = validate_registry([c for c, _, _ in checks])
    return checks, stops + s


# --------------------------------------------------------------- production
# The pinned inputs. An INPUT-INTEGRITY PRECONDITION and never evidence about
# the schema: reported separately, below the checks, so the two cannot be
# confused. P3 does not regenerate anything on a mismatch; it stops.
PINS = {
    "p1_artifact": "93acd8376da92687626cb6715aae7a5cd35c8adbb8c9d3eb7a0fd2ee006b3df4",
    "packet": "744c7f25e2312d90fc356b11510da685328f05e80ae62721d0a0f418dcf9697e",
    "sidecar": "76705591ee67390ce6572fc05cc9fcf40686a0c6ac109997b056b978ee622ba4",
}


def check_pins(paths):
    out = []
    for name, path in paths.items():
        got = sha256_bytes(pathlib.Path(path).read_bytes())
        out.append((name, got == PINS[name], got))
    return out


def main():
    ap = argparse.ArgumentParser(description="P3 mechanical schema gate")
    ap.add_argument("--p1-artifact", required=True)
    ap.add_argument("--packet", required=True)
    ap.add_argument("--sidecar", required=True)
    ap.add_argument("--group-packet", required=True)
    ap.add_argument("--construction-packet", required=True)
    a = ap.parse_args()

    print("P3 schema gate\n")
    pins = check_pins({"p1_artifact": a.p1_artifact, "packet": a.packet,
                       "sidecar": a.sidecar})
    print("  input integrity, a PRECONDITION and not evidence about the schema:")
    for name, ok, got in pins:
        print(f"    {'ok  ' if ok else 'FAIL'} {name:12s} {got[:16]}")
    if not all(ok for _, ok, _ in pins):
        print("\n  STOP. An input is not the pinned one. P3 does not regenerate it;\n"
              "  something upstream moved and that is what needs finding.")
        return 1

    pb = pathlib.Path(a.packet).read_bytes()
    sb = pathlib.Path(a.sidecar).read_bytes()
    p1 = json.loads(pathlib.Path(a.p1_artifact).read_bytes())
    gb = pathlib.Path(a.group_packet).read_bytes()
    cb = pathlib.Path(a.construction_packet).read_bytes()

    checks, stops = gate(pb, sb, p1, gb, cb)
    print()
    w = max(len(c) for c, _, _ in checks)
    for cid, ok, detail in checks:
        print(f"  {'ok  ' if ok else 'FAIL'} {cid:{w}s}  {detail}")
    print()
    if stops:
        print("  STOP. The candidate does not pass the schema gate.")
        for s in stops:
            print(f"    - {s}")
        return 1
    print(f"  P3 GREEN on all {len(checks)} declared checks.")
    print("  The candidate is a valid seven-key packet and sidecar under the")
    print("  frozen schema and § 10.2. This says nothing about whether the")
    print("  values are correct, which is P5's question.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
