"""The AUTHOR-SIDE REFERENCE COMPARISON FIXTURE.

WHAT THIS IS FOR, stated before anything else, because the claim is narrow.

The § 4.4 ingestion controls are properties of a harness that CONSUMES the answer
packet. We do not own that harness: it is the maintainer's, and they must
independently bind and mutation-test the corresponding controls there. What P4
can establish is that the controls are IMPLEMENTABLE against the frozen packet,
which is to say the packet is CONSUMABLE in the way § 4.4 requires. This fixture
is the minimum consumer needed to demonstrate that, and nothing here claims his
harness is correct.

The comparison it performs is the one § 5.3 and § 5.4 describe:

  1. load the canonical bytes and hash-verify BEFORE parsing;
  2. JOIN the committed rows to the packet's rows BY SIGNATURE, never by label;
  3. APPLY the indexing map and require it to agree with that join;
  4. locate the anchor by the SIGNATURE the packet carries for it, not by label;
  5. select the orientation, per § 5.4's rule;
  6. recompute the identities position-wise from the SELECTED committed values;
  7. compare.

Steps 2 and 4 were absent from the first version of this fixture, which took the
committed output as a bare list of values, aligned it positionally, and found the
anchor through a label in the P1 artifact. It therefore never read a committed
signature at all: swapping two signatures on the committed side left it agreeing.
It claimed to resolve by signature and did not. § 5.5 makes the signature the
join key and says the raw output carries it, so a consumer that aligns
positionally is not the consumer the protocol describes.

TWO MECHANISMS, CERTIFYING TWO DIFFERENT THINGS, and both must stay.

  MAP CORRECTNESS. The signature join independently derives what the
  correspondence OUGHT to be, and the map is required to state exactly that.

  MAP APPLICATION. The comparison then walks THE MAP, not the join. The map is
  the mechanism the comparison path actually follows, so a consumer that
  validated it and then ignored it produces a different answer.

An intermediate version collapsed the two: it joined by signature and compared
through the join, leaving the map a cross-check only. That silently gave up the
§ 4.4 property the synthetic nonidentity fixture exists to demonstrate.
`corr_mode="join"` models that harness exactly, validation included.

And it cannot be caught by comparing VALUES. A correct map and the join agree by
construction, since validation requires it, so both modes reach the same answer.
The ingestion control therefore measures whether the map's destinations are READ
during the comparison, which is the only thing that separates them.

NO TRANSCRIPTION. There is not a single value, triple or expected result written
into this file. Everything compared comes from the bytes handed in. That is
control (5), and it is checkable by reading the file.

Arithmetic is `p3_qphi`, which shares no code with the emitter.
"""

import hashlib
import json

import p4_verify
import pathlib

import p3_qphi as Q


# Raised BEFORE the parser is reached. Defined in `p4_verify`, and re-exported
# here so callers keep one name: a second exception class would let the loader
# raise one the verifier never raises, which is the sort of seam this control has
# already failed to hold through.
HashMismatch = p4_verify.HashMismatch


def load_hash_verified(path, expected_sha256, parse=json.loads):
    """Verify, and only then parse. ONE expression, deliberately.

    The body is a single composition and must stay one. `p4_verify` cannot parse
    JSON, and `p4.py` proves that by running it where `json` does not exist; this
    call site is where that guarantee turns into an ORDERING, because a parse
    that ran before the hash check would need a statement to live in and there is
    no statement here. The control checks that on the syntax tree.

    `parse` is injectable for one reason: the ingestion tests pass a sentinel, so
    the bytes handed to the parser can be compared to the bytes that were hashed.
    Hashing the canonical file and parsing a re-rendering of it did not hold an
    earlier version. `json.loads` takes bytes, so there is nothing to re-render.
    """
    return parse(p4_verify.verified_bytes(path, expected_sha256))


class _IntLike:
    """What counts as a position besides a real `int`.

    The ingestion control substitutes objects that record being used as an index;
    they are deliberately not `int` subclasses, so the integer check has to admit
    them by capability rather than by type. `bool` is excluded explicitly above,
    since it IS an `int` subclass and that is exactly the hole."""

    def __instancecheck__(self, obj):                         # pragma: no cover
        return hasattr(type(obj), "__index__")


_IntLike = _IntLike()


def signature_of(row):
    s = row["signature"]
    return (s["dimension"], tuple(s["s"]), tuple(s["t"]), tuple(s["st"]))


def signature_key(obj):
    """The § 5.5 join key, from either a packet row or a committed row."""
    s = obj["signature"]
    return (s["dimension"], tuple(s["s"]), tuple(s["t"]), tuple(s["st"]))


def join_by_signature(packet, committed_rows):
    """The label-free correspondence. Raises with what went wrong."""
    prows = packet["rows"]
    psigs = [signature_key(r) for r in prows]
    csigs = [signature_key(r) for r in committed_rows]
    if len(set(psigs)) != len(psigs):
        raise ValueError("packet rows do not have distinct signatures")
    if len(set(csigs)) != len(csigs):
        raise ValueError("committed rows do not have distinct signatures")
    if set(psigs) != set(csigs):
        raise ValueError(
            f"signature sets differ: {len(set(psigs) - set(csigs))} packet row(s) "
            f"unmatched, {len(set(csigs) - set(psigs))} committed row(s) unmatched")
    where = {s: i for i, s in enumerate(csigs)}
    return {s: (prows[i], where[s]) for i, s in enumerate(psigs)}


# The indexing map's HEADER, as this fixture's own code interprets it. The
# consumer reads `entries` and treats each one as a zero-based packet-row
# position mapped to a zero-based committed-row position; those three
# declarations say exactly that, and it never checked any of them. Setting
# `zero_based` false, or renaming either domain, was accepted with all nine rows
# and four identities agreeing. Same class as the convention-map leaves: relying
# on a declaration you do not read.
INDEXING_MAP_HEADER = {
    "source_domain": "packet_rows_canonical_position",
    "destination_domain": "m8_3_record_row_position",
    "zero_based": True,
}


def check_indexing_map(packet, committed_rows, join):
    """APPLY the map and require it to agree with the signature join.

    The map says packet canonical position -> the committed output's own
    position; the join says the same thing by signature. If they disagree, the
    packet describes a correspondence its own rows contradict.

    First, though, the map must MEAN that. The header declarations are checked
    against what this code implements, and positions must be true integers:
    `False == 0` in Python and a bool indexes a list, so a position written as
    JSON `false` was silently accepted as position zero."""
    imap = packet["indexing_map"]
    for key, expected in INDEXING_MAP_HEADER.items():
        # A declared BOOLEAN must be a bool, by identity. `1 == True` in Python,
        # so `!= True` accepted a `zero_based` written as JSON `1`: the mirror
        # image of the `False == 0` hole the position type-check exists for.
        # Strings still compare by value, which is what they should do.
        got = imap.get(key)
        wrong = (got is not expected) if isinstance(expected, bool) else (got != expected)
        if wrong:
            raise ValueError(f"indexing_map.{key} is {imap.get(key)!r}, not the "
                             f"{expected!r} this consumer implements")
    rows = packet["rows"]
    seen = set()
    for entry in packet["indexing_map"]["entries"]:
        for field in ("source_position", "destination_position"):
            v = entry.get(field)
            # A tagged destination is the ingestion control's instrument and is
            # an integer by contract; a bool is not, and never should be.
            if isinstance(v, bool) or not isinstance(v, (int, _IntLike)):
                raise ValueError(f"indexing map {field} {v!r} is not an integer")
        src, dst = entry["source_position"], entry["destination_position"]
        # `dst` may be a tagged object rather than a bare int: the ingestion
        # control substitutes one to measure whether the comparison actually
        # USES map destinations as indices, as opposed to merely reading them.
        # Only ordering and equality are needed here, never indexing.
        if not (0 <= src < len(rows) and 0 <= dst < len(committed_rows)):
            raise IndexError(f"indexing entry {src} -> {dst} is out of range")
        if dst in seen:
            raise ValueError(f"two sources map to destination {dst}")
        seen.add(dst)
        _row, joined = join[signature_key(rows[src])]
        if joined != dst:
            raise ValueError(f"indexing map sends canonical row {src} to {dst}, but "
                             f"its signature matches committed row {joined}")
    if len(seen) != len(committed_rows):
        raise ValueError("indexing map does not cover every destination position")


# The convention-map semantics THIS FIXTURE'S OWN CODE implements. A consumer
# validates the declarations it relies on; that is what "loads and validates its
# convention map" means, and it is the line that keeps this from becoming a
# second P3. Each entry below corresponds to behaviour in this file: the
# orientation transform inverts every row at once, the selection follows x = r
# and x^-1 = r, and `p3_qphi` assumes exactly the declared encoding.
# The convention-map declarations THIS FIXTURE'S OWN CODE relies on, held as
# COMPLETE EXPECTED VALUES rather than as required substrings.
#
# Substring presence is not validation. A leaf can carry the required phrase and
# explicitly assert its negation: "DO NOT apply the inverse to EVERY ROW AT ONCE;
# apply it only to R7" contains the phrase and means the opposite, and the
# substring version accepted it. What the consumer needs to know is that the
# declaration it is relying on says what it thinks it says, and only an exact
# comparison establishes that.
#
# These are protocol prose, not answer values; the transcription screen looks
# for value-shaped literals and is unaffected.
# The schema string this fixture is written against. Every field access below
# is specific to it, and it was neither read nor validated: setting it to an
# incompatible version, or deleting it outright, left the comparison agreeing on
# all nine rows and four identities. Seventh occurrence of relying on a
# declaration the code never reads, one level above the indexing map's header.
FORMAT_VERSION = "m8_8-answer-packet-1"

RELIED_ON_DECLARATIONS = {
    ("bridge", "identification"):
        'T^2_target(rho) := |tau_rho|^2',
    ("bridge", "exclusivity"):
        "the SOLE admitted convention bridge between the routes is the preregistered global sign of log T^2, equivalently the inversion T^2 <-> (T^2)^-1 applied to EVERY ROW AT ONCE and anchored once at R7; this is what 'global inverse' means wherever it appears here. Outputs not related by exactly that transformation are a DISAGREEMENT, never a convention difference",
    ("orientation_anchor_rule", "selection_mapping"):
        "writing x for the committed native R7 value and r for this packet's R7 reference value, whichever orientation agrees at R7 is the selected one: x = r selects the COMMITTED table as committed, x^-1 = r selects its GLOBAL INVERSE, and the full comparison proceeds under the selected orientation",
    ("orientation_anchor_rule", "uniqueness_gate"):
        'the selection must be unique and that is a gate, not an assumption: exactly one of x = r and x^-1 = r may hold. Both holding means r = r^-1, a self-inverse reference that cannot discriminate orientation, and the adjudicator records an INVALID ANCHOR and selects no orientation. Neither holding is a § 8 disagreement, never a convention difference',
}

# THE Q(phi) ENCODING IS NOT IN THAT LIST ANY MORE, and the reason matters more
# than the four entries do. This fixture's arithmetic depends on it, so for most
# of the build it was listed here and validated against the packet's own
# declaration. An independent protocol-first read then adjudicated that
# declaration out of the packet: § 4.4 assigns the encoding to `rows` and
# enumerates a different content for the convention map, so carrying it in both
# places created two statements that could drift.
#
# What the fixture relies on has not changed. WHERE THAT RELIANCE IS ANSWERED
# has. The encoding is a PROTOCOL CONSTANT, fixed by § 4.4's prose, and this code
# implements it directly rather than reading it back from the object under test.
# There is nothing left here to validate, which is the correct end state: a
# declaration you validate against the thing that declared it was never
# independent evidence in the first place.

def validate_convention_map(packet, construction_bytes):
    """Require the packet to DECLARE the semantics this fixture implements.

    Without this the fixture hardcoded the semantics and never read the
    declarations: replacing `bridge.identification` with nonsense, or the
    declared field with a different number field, left the comparison green. A
    consumer that assumes semantics it never checked is not validating its
    convention map."""
    if packet.get("format_version") != FORMAT_VERSION:
        raise ValueError(f"format_version is {packet.get('format_version')!r}, "
                         f"not the {FORMAT_VERSION!r} this consumer implements")
    cm = packet.get("convention_map")
    if not isinstance(cm, dict):
        raise ValueError("convention_map is not an object")
    for (section, leaf), expected in RELIED_ON_DECLARATIONS.items():
        block = cm.get(section)
        if not isinstance(block, dict) or not isinstance(block.get(leaf), str):
            raise ValueError(f"convention_map.{section}.{leaf} is missing")
        if block[leaf] != expected:
            raise ValueError(f"convention_map.{section}.{leaf} is not the declaration "
                             "this consumer relies on")
    basing = cm.get("basing_reference")
    if not isinstance(basing, dict):
        raise ValueError("convention_map.basing_reference is missing")
    if basing.get("construction_packet_sha256") != hashlib.sha256(
            construction_bytes).hexdigest():
        raise ValueError("the basing reference does not pin the construction packet "
                         "actually supplied")


def check_adjudicates(packet, group_bytes, construction_bytes):
    """§ 4.4: `adjudicates` must match the run's ACTUAL packets.

    A consumer-side control, and it was absent: the fixture took no public
    packets at all, so zeroing an `adjudicates` hash after a successful
    hash-verified load left the comparison green."""
    adj = packet.get("adjudicates")
    if not isinstance(adj, dict):
        raise ValueError("adjudicates is not an object")
    want = {"group_packet_sha256": hashlib.sha256(group_bytes).hexdigest(),
            "construction_packet_sha256": hashlib.sha256(construction_bytes).hexdigest()}
    if adj != want:
        wrong = sorted(k for k in want if adj.get(k) != want[k])
        raise ValueError(f"`adjudicates` does not match the packets supplied to this "
                         f"run: {wrong}")


# § 5.4 does not leave the anchor open: "only after the answer packet opens may
# the adjudicator use `R7` to select between the committed table and its global
# inverse", and § 5.2 gives R7 its standing as the orientation selector. So the
# protocol DESIGNATES the anchor by evidentiary class, and a packet that declares
# some other row's signature as its anchor is not exercising a freedom.
ANCHOR_CLASS = "free_orientation_selector"


def anchor_signature(packet):
    """The anchor's signature, taken from the declaration and then CHECKED.

    § 5.4 carries the anchor's SIGNATURE in the convention map precisely so the
    reference resolves without the public label. The fixture consumed that
    declaration and never validated it: substituting any other row's signature
    was accepted for six of the eight alternatives, and the § 5.4 selection then
    ran on the wrong row while convention-map validation stayed green. The two
    that failed did so only because they happen to be self-inverse.

    This is the claimant-versus-designated-subject shape again. The rule comes
    from the protocol; which row carries the class comes from the packet; the
    consumer requires the two to agree."""
    rule = packet["convention_map"]["orientation_anchor_rule"]
    declared = signature_key({"signature": rule["anchor_row_signature"]})
    owning = [r for r in packet["rows"] if r.get("evidentiary_class") == ANCHOR_CLASS]
    if len(owning) != 1:
        raise ValueError(f"{len(owning)} rows carry the {ANCHOR_CLASS} class; § 5.2 "
                         "gives it to exactly one")
    if declared != signature_key(owning[0]):
        raise ValueError("the declared anchor is not the row § 5.4 designates: "
                         f"the {ANCHOR_CLASS} row has a different signature")
    if rule.get("anchor_row_label") != owning[0].get("label"):
        raise ValueError("the anchor rule's label and the designated row disagree")
    return declared


def select_orientation(anchor_committed, anchor_reference):
    """§ 5.4's one-bit selection, and its gate.

    Returns "committed", "inverse", "INVALID_ANCHOR" or "DISAGREEMENT". Both
    branches holding means a self-inverse reference, and the protocol says the
    adjudicator then selects nothing rather than silently picking one."""
    x, r = tuple(anchor_committed), tuple(anchor_reference)
    direct, inverted = x == r, Q.inv(x) == r
    if direct and inverted:
        return "INVALID_ANCHOR"
    if direct:
        return "committed"
    if inverted:
        return "inverse"
    return "DISAGREEMENT"


def compare(packet, committed_rows, group_bytes, construction_bytes,
            corr_mode="map"):
    """`committed_rows` are the implementation's own records, each carrying the
    § 5.5 signature and its value, in the implementation's own order.

    `corr_mode` selects what the comparison FOLLOWS. "map" is production.
    "join" models the documented regression: a consumer that validates the map
    and then aligns through the signature join instead.

    A correct map and the join agree BY CONSTRUCTION, since validation requires
    it, so no comparison of VALUES can distinguish the two modes. An earlier
    control tried to, using positional identity as the shortcut, which is a
    different and easier failure than the one that actually occurred. Telling
    them apart requires instrumentation on whether the map's destinations are
    read during the comparison, and that is what the ingestion control measures.

    Returns (ok, detail). Pure: nothing is written anywhere."""
    rows = packet["rows"]
    try:
        check_adjudicates(packet, group_bytes, construction_bytes)
        validate_convention_map(packet, construction_bytes)
        join = join_by_signature(packet, committed_rows)
        check_indexing_map(packet, committed_rows, join)
        if corr_mode == "map":
            corr = {e["source_position"]: e["destination_position"]
                    for e in packet["indexing_map"]["entries"]}
        elif corr_mode == "join":
            corr = {i: join[signature_key(r)][1] for i, r in enumerate(rows)}
        else:
            raise ValueError(f"unknown corr_mode {corr_mode!r}")
        asig = anchor_signature(packet)
        aidx = next((i for i, r in enumerate(rows) if signature_key(r) == asig), None)
        if aidx is None:
            return False, "the anchor signature the packet declares matches no row"
        orientation = select_orientation(committed_rows[corr[aidx]]["value"],
                                         rows[aidx]["value"])
    except Exception as e:                                    # noqa: BLE001
        return False, f"{type(e).__name__}: {e}"

    if orientation in ("INVALID_ANCHOR", "DISAGREEMENT"):
        return False, orientation

    def selected(v):
        v = tuple(v)
        return v if orientation == "committed" else Q.inv(v)

    differ = [i for i in range(len(rows))
              if selected(committed_rows[corr[i]]["value"]) != tuple(rows[i]["value"])]
    if differ:
        return False, (f"orientation {orientation}; {len(differ)} row(s) differ after "
                       f"selection, following the {corr_mode}")

    for slot, ident in enumerate(packet["identities"]):
        if ident["position"] != slot:
            return False, f"identity slot {slot} declares position {ident['position']}"
        acc = Q.ONE
        for factor in ident["factors"]:
            fidx = next((i for i, r in enumerate(rows)
                         if signature_key(r) == signature_key(factor)), None)
            if fidx is None:
                return False, f"slot {slot}: a factor resolves to no row"
            acc = Q.mul(acc, Q.power(selected(committed_rows[corr[fidx]]["value"]),
                                     factor["exponent"]))
        if list(acc) != ident["expected_value"]:
            return False, (f"identity slot {slot} recomputes to {list(acc)}, "
                           f"packet declares {ident['expected_value']}")

    return True, (f"orientation {orientation}; {len(join)} rows joined by signature, "
                  f"compared through the {corr_mode}, and "
                  f"{len(packet['identities'])} identities agree")
