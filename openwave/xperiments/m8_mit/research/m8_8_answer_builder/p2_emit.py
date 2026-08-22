"""P2 emission: build the CANDIDATE answer packet from the P1 artifact.

Target-free. Every source value, signature and identity factor is read out of the
P1 cross-check table. Identity expected values are the ONE derived field: they
are recomputed exactly from the committed rows rather than copied from P1's
stated right-hand sides, so a row and its identity cannot drift apart. P1 is
closed and its artifact hash is the sole input version: re-deriving any of it
during emission would create a second, unaudited path to the same numbers.

CANDIDATE, not the issued record. Immutability attaches at P6 issuance, after
the schema gate, the mutation battery and the cold read are all green. P3, P4
and P5 may reject and revise freely until then; calling this immutable would
mean a later gate failure forced a formal reissue of something never issued.
"""

import json
import pathlib

import p2_schema as S

# The banner and the P1 pin live in the schema, which is their single owner.
# Two dead copies of the banner text sat here and in the driver; both were
# unreferenced, and an edit to either would have looked effective while
# changing nothing. Aliases kept only so existing references resolve.
BANNER = S.HASH_RECORD_BANNER
P1_ARTIFACT_SHA256 = S.P1_ARTIFACT_SHA256
SELECTOR_LABEL = S.SELECTOR_LABEL
DECLARED_LABEL = S.DECLARED_LABEL


def load_p1(path):
    """Read the P1 artifact and refuse anything but the pinned bytes."""
    import hashlib
    raw = pathlib.Path(path).read_bytes()
    got = hashlib.sha256(raw).hexdigest()
    if got != S.P1_ARTIFACT_SHA256:
        raise ValueError(
            f"P1 artifact hash {got[:16]} is not the pinned "
            f"{S.P1_ARTIFACT_SHA256[:16]}. P2 consumes exactly one P1 version; "
            "regenerating P1 means re-running P1.")
    return json.loads(raw)


def build_rows(p1):
    """One row per irrep, in canonical sort-key order.

    The array is emitted ALREADY ordered, so a consumer can verify the order
    rather than infer it."""
    values = {v["label"]: v["triple"] for v in p1["values"]}
    rows = []
    for label, sig in p1["signatures"].items():
        if label == S.DECLARED_LABEL:
            cls = S.CLASS_DECLARED
        elif label == S.SELECTOR_LABEL:
            cls = S.CLASS_SELECTOR
        else:
            cls = S.CLASS_FREE
        rows.append({
            "label": label,                       # display metadata only
            "signature": {"dimension": sig["dimension"], "s": list(sig["s"]),
                          "t": list(sig["t"]), "st": list(sig["st"])},
            "evidentiary_class": cls,
            "value": list(values[label]),
        })
    rows.sort(key=S.row_sort_key)
    return rows


def build_identities(p1, rows):
    """Identity slots with explicit zero-based contiguous positions.

    Factors carry the row SIGNATURE and an integer exponent, never a label:
    § 5.3 wants the formula defined label-free, and a label join is exactly the
    path by which a mislabelled row would silently take another row's value."""
    by_label = {r["label"]: r for r in rows}
    out = []
    for pos, ident in enumerate(p1["identities"]):
        factors = []
        for label, exponent in ident["factors"]:
            factors.append({
                "signature": dict(by_label[label]["signature"]),
                "exponent": int(exponent),
            })
        factors.sort(key=S.factor_sort_key)   # frozen order, not source order
        out.append({
            "position": pos,
            "factors": factors,
            "expected_value": list(_expected(ident, by_label)),
        })
    return out


def _expected(ident, by_label):
    """The identity's expected value, recomputed from the committed rows.

    Not read from the record's stated right-hand side: § 5.3 requires the
    expected value to be recomputable from the rows the packet itself commits,
    so a row and its identity cannot drift apart unnoticed. P1 already checked
    the two agree (IDY_HOLDS)."""
    from qphi_exact import Phi
    acc = Phi(1)
    for label, exponent in ident["factors"]:
        a, b, c = by_label[label]["value"]
        acc = acc * ((Phi(a, b) / Phi(c)) ** int(exponent))
    return acc.triple()


def build_indexing_map(rows, p1):
    """The correspondence between the packet's canonical row order and the
    reference order used by the M8.3 record.

    Both domains are declared, both are zero-based, and both are required
    contiguous. It is a genuine NON-IDENTITY permutation here, because sorting
    by the ten-integer signature key does not agree with the record's own row
    order, which is what makes it able to discriminate a harness that applies
    the map from one that ignores it."""
    ref_order = [v["label"] for v in p1["values"]]
    pos_of_ref = {lab: i for i, lab in enumerate(ref_order)}
    return {
        "source_domain": "packet_rows_canonical_position",
        "destination_domain": "m8_3_record_row_position",
        "zero_based": True,
        "entries": [{"source_position": i,
                     "destination_position": pos_of_ref[r["label"]]}
                    for i, r in enumerate(rows)],
    }


def build_convention_map(p1, rows):
    """§ 5.4 bridge and anchor rule, § 4.2 basing, and the two conventions.

    The complete expected content lives in the schema and is compared exactly.
    The construction hash, the selector's own signature and § 4.2's
    representation-evaluation convention are the dynamic fields; the last is READ
    from the construction packet rather than composed here."""
    sel = [r for r in rows if r["evidentiary_class"] == S.CLASS_SELECTOR]
    if len(sel) != 1:
        raise ValueError(f"{len(sel)} rows carry the selector class, expected 1")
    return S.convention_map_template(
        p1["adjudicates"]["construction_packet"],
        S.SELECTOR_LABEL,
        dict(sel[0]["signature"]),
        S.representation_evaluation_from_construction(
            p1["adjudicates"]["construction_packet"]))


def build_packet(p1):
    rows = build_rows(p1)
    return {
        "format_version": S.FORMAT_VERSION,
        "target_id": S.TARGET_ID,
        "adjudicates": {
            "group_packet_sha256": p1["adjudicates"]["group_packet"],
            "construction_packet_sha256": p1["adjudicates"]["construction_packet"],
        },
        "rows": rows,
        "identities": build_identities(p1, rows),
        "indexing_map": build_indexing_map(rows, p1),
        "convention_map": build_convention_map(p1, rows),
    }


def build_hash_record(packet_bytes):
    """The CANDIDATE pre-lock hash record: sidecar object (a), not yet issued.

    Ships the canonicalization rule WITH the hash, because a quarantined object
    offers no other checkable property and a hash over an uncanonical rendering
    pins a transcription rather than an object.

    It does NOT carry § 10.2's incoming-bytes hash or its statement of whether
    canonicalization changed the bytes, and this docstring used to claim it did.
    Those attach to the ISSUED authoritative hash; this record is explicitly
    pre-issuance, and the emitter writes the canonical bytes directly, so there
    is no delivered rendering for either field to describe. They belong to the
    P6 record, and pulling P6 semantics backward to satisfy a docstring would
    have been the wrong repair."""
    import hashlib
    canon = hashlib.sha256(packet_bytes).hexdigest()
    rec = S.hash_record_template(canon, len(packet_bytes))
    return rec
