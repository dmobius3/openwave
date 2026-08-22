"""The frozen nested schema for the § 4.4 answer packet.

Target-free: this file declares STRUCTURE, never a value or a membership. Every
value, signature and identity factor comes from the P1 artifact at build time.

Everything structural is frozen HERE rather than decided at emission time, for
the reason the plan gives: the mutation suite needs an unambiguous target, and
"sorted by the signature" is exactly the unstated-ordering defect one layer up.
A signature is ten integers and naming it does not order them.
"""

# § 4.4: exactly these keys, no more and no fewer.
PACKET_KEYS = ("adjudicates", "convention_map", "format_version", "identities",
               "indexing_map", "rows", "target_id")

FORMAT_VERSION = "m8_8-answer-packet-1"

# The P1 artifact this build consumes, pinned. Held HERE rather than in the
# emitter, for the reason the designators below give: while the module that
# builds the packet from P1 also owned the hash that says WHICH P1 is
# legitimate, a single-file edit moved the artifact and its expectation
# together and no check could notice.
#
# Recorded plainly, because the move is organization rather than a new trust
# boundary: an edit to this literal is still accepted, exactly as an edit to
# DECLARED_LABEL would be. What it buys is that the emitter can no longer
# relicense its own input. A separately frozen input-version record consumed by
# both sides is the stronger form and is not built here.
P1_ARTIFACT_SHA256 = \
    "93acd8376da92687626cb6715aae7a5cd35c8adbb8c9d3eb7a0fd2ee006b3df4"

# The PUBLIC protocol designators (§ 5.2, § 5.4) and the row count, frozen HERE
# so the emitter and the checker cannot move together. While the emitter owned
# them, a coherent change to SELECTOR_LABEL would have moved the emitted class,
# the expected assignment, the anchor label and the anchor signature in step:
# they would agree with each other and with nothing else.
DECLARED_LABEL = "R0"
SELECTOR_LABEL = "R7"
ROW_COUNT = 9

# The sidecar banner, frozen here for the same reason. While the emitter owned
# it, "only the hash and the byte length are dynamic" was not literally true.
#
# It used to bar the P5 reader, and that was WRONG. What the quarantine protects
# is the reproducing context: the party computing these values by the other
# route must not receive the answer key. P5 is not that party. P5 recomputes
# faithful capture FROM the public sources and needs this record to do it, and
# REDLINE_PROTOCOL § 2.1 hands it over by name, "the immutable pre-lock hash
# record". A banner barring P5 sat inside the file § 2.1 gives P5, so whichever
# document a reader happened to hold decided the answer. That is the same defect
# P1 recorded fixing between its report and the bundle manifest, recurring
# between this build's artifact and the governing document. § 2.1's withhold
# list names "our cross-check table", which is P1's artifact and not this one.
HASH_RECORD_BANNER = (
    "QUARANTINED, ANSWER-BEARING. CANDIDATE, not issued.\n"
    "BARRED from the reproducing context and every route into it: the PR,\n"
    "commit messages, the protocol, the public ciphertext README, and any\n"
    "§ 4.3 permitted input.\n"
    "TWO authorized destinations and no others: an author-side informed redline\n"
    "session, which is permanently author-side once it receives this; and the P5\n"
    "independent reader, whose declared input set includes this record per\n"
    "REDLINE_PROTOCOL § 2.1.\n"
)
TARGET_ID = "M88-ADJ-01"

# § 5.2, verbatim. Assigned per row at build time from the P1 artifact's labels;
# labels are display metadata in the EMITTED packet and are never a join key
# there, which is a different question from how the author-side build indexes
# its own inputs.
CLASS_DECLARED = "declared_convention"
CLASS_FREE = "free"
CLASS_SELECTOR = "free_orientation_selector"
CLASS_CENSUS = {CLASS_DECLARED: 1, CLASS_FREE: 7, CLASS_SELECTOR: 1}

IDENTITY_SLOTS = 4

# A bound on identity factor exponents, to bound COMPUTATION rather than to
# constrain mathematics. The exponent is a JSON integer literal and the exact
# recompute raises it before anything considers magnitude, so a 13-byte literal
# bought an effectively unbounded computation ending in MemoryError: a
# JSON-reachable raise. Every § 5.3 exponent in the record is +/-1, so this is
# generous by any measure; what actually pins the exponents to P1's is
# IDY_DEFINITIONS_MATCH_P1, and this only stops the arithmetic running away
# before that check is reached.
EXPONENT_ABS_MAX = 64

# A bound on FACTORS PER SLOT, for the same reason and because the exponent
# bound alone was not enough. Repeated factors are shape-valid and each resolves
# to exactly one row, so a packet of a few tens of KB drove assess() past 400
# seconds toward the same MemoryError. Every P1 slot carries 2 or 4.
FACTORS_PER_SLOT_MAX = 16

# The THIRD operand. The comment here used to name all three, "exponent x factor
# count x value magnitude", bound two of them, and then claim the amplification
# was a fixed constant. That claim was false as written and the symptom survived
# it: a shape-valid packet no larger than the real one still drove assess() to
# hundreds of seconds, because a row value is a JSON integer literal of any
# length and the recompute raises it to a bounded power a bounded number of
# times, which is still unbounded when the base is not.
#
# Bounding the digit count of every integer inside a triple closes it. With all
# three bounded, each slot costs a fixed constant and the total is linear in the
# file's own size, which is what the earlier comment claimed prematurely. The
# largest component anywhere in the real packet is two digits.
VALUE_ABS_DIGITS_MAX = 64


def row_sort_key(row):
    """The frozen ten-integer key, `dimension` most significant, compared
    lexicographically as SIGNED integers. Never by label."""
    s = row["signature"]
    return (s["dimension"],
            s["s"][0], s["s"][1], s["s"][2],
            s["t"][0], s["t"][1], s["t"][2],
            s["st"][0], s["st"][1], s["st"][2])


def signature_of(row):
    """The label-free identity, as a hashable tuple. This is the ONLY join key."""
    s = row["signature"]
    return (s["dimension"], tuple(s["s"]), tuple(s["t"]), tuple(s["st"]))


def is_triple_shape(t):
    """Exactly three members, each a true int.

    `type(x) is int` rather than isinstance, because bool subclasses int and
    Python says True == 1. Declaring a field as `list` froze the container and
    left every leaf open: JSON `true` in a triple position survived a round trip
    and compared equal to 1 in every arithmetic check downstream. Same class as
    the 1 vs 1.0 finding, one level deeper."""
    return (isinstance(t, list) and len(t) == 3
            and all(type(x) is int for x in t))


def triple_shape_problems(obj, where):
    """Shape AND magnitude. Magnitude here bounds COST, not mathematics: these
    integers are the base of the identity recompute's exponentiation, and an
    unbounded base defeats the exponent and factor-count bounds. What pins the
    values to P1's is ROWS_VALUES_MATCH_P1."""
    if not is_triple_shape(obj):
        return [f"{where}: not exactly three integer members (bools excluded)"]
    over = [i for i, x in enumerate(obj) if len(str(abs(x))) > VALUE_ABS_DIGITS_MAX]
    return [] if not over else [
        f"{where}: component(s) {over} exceed {VALUE_ABS_DIGITS_MAX} digits"]


def signature_leaf_problems(sig, where):
    """Every triple inside a signature object, checked at the leaf."""
    if not isinstance(sig, dict):
        return []
    out = []
    if type(sig.get("dimension")) is not int:
        out.append(f"{where}.dimension: not a true int")
    for k in ("s", "t", "st"):
        out += triple_shape_problems(sig.get(k), f"{where}.{k}")
    return out


def is_normalized_triple(t):
    """(a, b, c) with c > 0 and gcd(a, b, c) = 1. Integers only, never a float."""
    from math import gcd
    if not is_triple_shape(t):
        return False
    a, b, c = t
    return c > 0 and gcd(gcd(abs(a), abs(b)), c) == 1


def canonical_bytes(obj):
    """§ 10.2: JSON, keys sorted, two-space indent, ASCII, LF, one trailing
    newline. Returned as bytes so the caller hashes what it writes."""
    import json
    text = json.dumps(obj, sort_keys=True, indent=2, ensure_ascii=True)
    return (text + "\n").replace("\r\n", "\n").encode("ascii")


# ---------------------------------------------------------------------------
# EXACT NESTED SHAPE. Freezing only the seven top-level keys left every nested
# object open: an undeclared field could be added to a row and every check still
# passed, while a missing nested field made the checker dereference it and raise.
# Structure is frozen all the way down, or it is not frozen.

ROW_KEYS = {"label": str, "signature": dict, "evidentiary_class": str, "value": list}
SIGNATURE_KEYS = {"dimension": int, "s": list, "t": list, "st": list}
IDENTITY_KEYS = {"position": int, "factors": list, "expected_value": list}
FACTOR_KEYS = {"signature": dict, "exponent": int}
INDEXING_KEYS = {"source_domain": str, "destination_domain": str,
                 "zero_based": bool, "entries": list}
INDEXING_ENTRY_KEYS = {"source_position": int, "destination_position": int}
ADJUDICATES_KEYS = {"group_packet_sha256": str, "construction_packet_sha256": str}

HASH_RECORD_KEYS = {
    "_banner": str,
    "record_format_version": str,
    "status": str,
    "status_note": str,
    "answer_packet_format_version": str,
    "canonicalization_rule": str,
    "canonical_plaintext_byte_length": int,
    "expected_canonical_plaintext_sha256": str,
}

# The declared domains of the indexing map, compared EXACTLY. Requiring only a
# non-empty string accepted arbitrary domain names.
INDEXING_SOURCE_DOMAIN = "packet_rows_canonical_position"
INDEXING_DESTINATION_DOMAIN = "m8_3_record_row_position"

# The ORDER of the `entries` array, declared here rather than left implicit.
# The checker compares the array to a P1-derived list built in this order, so
# the order was already load-bearing while being named nowhere, which is the
# unstated-ordering defect this file exists to prevent. The bytes get hashed.
INDEXING_ENTRY_ORDER = "ascending source_position, contiguous from zero"


def factor_sort_key(factor):
    """The frozen order of identity factors: by the factor's row sort key, then
    exponent. Declared because leaving it to source order means two equivalent
    products can serialize differently, and the bytes are about to be hashed."""
    s = factor["signature"]
    return (s["dimension"], s["s"][0], s["s"][1], s["s"][2],
            s["t"][0], s["t"][1], s["t"][2],
            s["st"][0], s["st"][1], s["st"][2], factor["exponent"])


def representation_evaluation_from_construction(construction_sha256):
    """§ 4.2's evaluation convention, read from the construction packet itself.

    Hash-checked against the pin before parsing, so the string cannot come from
    an object other than the one `adjudicates` names."""
    import hashlib
    import json as _json
    import pathlib
    path = pathlib.Path(__file__).resolve().parent / "bundle" / \
        "construction__m8_8_construction_packet.json"
    raw = path.read_bytes()
    got = hashlib.sha256(raw).hexdigest()
    if got != construction_sha256:
        raise ValueError(f"construction packet is {got[:16]}, not the pinned "
                         f"{construction_sha256[:16]}")
    basing = _json.loads(raw)["basing"]
    return basing["evaluation"]


def convention_map_template(construction_sha256, selector_label, selector_signature,
                            representation_evaluation):
    """The COMPLETE expected convention map, compared exactly.

    Checking only that the subsection names are present let the bridge
    identification be replaced with anything at all. The construction hash, the
    selector's own signature and the § 4.2 representation-evaluation convention
    are the dynamic fields.

    `representation_evaluation` IS NOT TYPED HERE. It is read from the
    construction packet's `basing` block by the caller, which is the whole point:
    § 4.4 requires this map to carry § 4.2's evaluation convention, and for six
    build steps the field of that name carried the Q(phi) encoding instead. The
    reason nothing caught it is that this template is what the emitter builds
    from AND what the checker compares against, so P2's reading of § 4.4 could
    not be wrong in any way P2, P3 or P4 could see. An independent protocol-first
    read found it in one pass. The expectation for this field now comes from a
    pinned external object rather than from the packet's own design."""
    return {
        "bridge": {
            "protocol_section": "5.4",
            "identification": "T^2_target(rho) := |tau_rho|^2",
            # The LAST clause of § 5.4's enumerated identification was missing:
            # "and the involution below acting on T^2_target AFTER this
            # identification". It is not commentary, it is an ordering
            # requirement on the one admitted bridge, and a harness that
            # inverted before identifying would satisfy every other clause here.
            "statement": ("the squared modulus of the Reidemeister torsion of the "
                          "based acyclic complex in the declared bases, determinants "
                          "over C in the evaluated representation; no field norm to a "
                          "subfield, no absolute value elsewhere, no squaring beyond "
                          "the displayed one; and the involution acting on T^2_target "
                          "AFTER this identification, never before it"),
            # § 5.4's closing normative paragraph, "The bridge is one involution,
            # and nothing else". It does two things nothing else here did: it
            # DEFINES the global inversion that `rule` and `selection_mapping`
            # both invoke by name, and it forecloses every other transformation.
            # Without it a consumer could satisfy every clause present while
            # inverting rows individually, or admit a second transformation as a
            # convention difference. The term was used but never defined, which
            # is the defect that was closed at x and r one audit earlier.
            "exclusivity": (
                "the SOLE admitted convention bridge between the routes is the "
                "preregistered global sign of log T^2, equivalently the "
                "inversion T^2 <-> (T^2)^-1 applied to EVERY ROW AT ONCE and "
                "anchored once at R7; this is what 'global inverse' means "
                "wherever it appears here. Outputs not related by exactly that "
                "transformation are a DISAGREEMENT, never a convention "
                "difference"),
        },
        "orientation_anchor_rule": {
            "protocol_section": "5.4",
            # The label is the PUBLIC protocol reference (§ 5.2, § 5.4 name it).
            # The signature is carried alongside so the reference resolves
            # without the label, which is what keeps labels non-load-bearing.
            "anchor_row_label": selector_label,
            "anchor_row_signature": selector_signature,
            "rule": ("the one-bit selection between the committed table and its "
                     "global inverse is made by the ADJUDICATOR after this packet "
                     "opens, never by the implementation"),
            # § 5.4 makes uniqueness "a gate rather than an assumption", and the
            # packet is the harness's machine-readable input for the selection
            # it gates. Carrying the rule without its gate left the packet
            # silent on the case the protocol wrote the gate for, so a harness
            # could pick an orientation the protocol says is not selectable.
            "uniqueness_gate": (
                "the selection must be unique and that is a gate, not an "
                "assumption: exactly one of x = r and x^-1 = r may hold. Both "
                "holding means r = r^-1, a self-inverse reference that cannot "
                "discriminate orientation, and the adjudicator records an "
                "INVALID ANCHOR and selects no orientation. Neither holding is "
                "a § 8 disagreement, never a convention difference"),
            # WHICH branch selects WHICH orientation, and what x and r mean.
            # The rule said who selects and between what; the gate said when a
            # selection is valid; nothing said that the branch which HOLDS is
            # the branch that is CHOSEN. § 5.4: "whichever of the two
            # orientations agrees with the analytic convention at R7 is the
            # selected one, and the full comparison proceeds under it". A
            # consumer reading only this packet could satisfy every other clause
            # here and still select the inverse when the committed table is the
            # one that agrees. The notation was also used without ever being
            # defined, which § 5.4 does define.
            "selection_mapping": (
                "writing x for the committed native R7 value and r for this "
                "packet's R7 reference value, whichever orientation agrees at "
                "R7 is the selected one: x = r selects the COMMITTED table as "
                "committed, x^-1 = r selects its GLOBAL INVERSE, and the full "
                "comparison proceeds under the selected orientation"),
        },
        "basing_reference": {
            "protocol_section": "4.2",
            "source": "the construction packet's `basing` block",
            "construction_packet_sha256": construction_sha256,
        },
        # The Q(phi) encoding is NOT here. It was, under two names, and the
        # independent protocol-first read adjudicated it out: § 4.4 assigns the
        # encoding to `rows`, and enumerates a different content for this map.
        # Carrying it in both places created two declarations that could drift,
        # which is the failure this build has met more than any other. The
        # encoding is a protocol constant, not a packet declaration.
        # NEW. § 4.2's evaluation convention for group-ring entries, which the
        # protocol singles out because a `rho(g)` versus `rho(g^-1)^T` error is
        # CHARACTER-INVISIBLE: no character check can see it, so it has to be
        # declared. Quoted from the construction packet, not composed here.
        "representation_evaluation": {
            "protocol_section": "4.2",
            "convention": representation_evaluation,
            "source": "quoted from the construction packet's `basing.evaluation`",
        },
    }


def exact_shape(obj, spec, where):
    """Exact key set and value types. Returns a list of problems."""
    if not isinstance(obj, dict):
        return [f"{where}: not an object"]
    problems = []
    extra = sorted(set(obj) - set(spec))
    missing = sorted(set(spec) - set(obj))
    if extra:
        problems.append(f"{where}: undeclared field(s) {extra}")
    if missing:
        problems.append(f"{where}: missing field(s) {missing}")
    for k, ty in spec.items():
        if k in obj and not isinstance(obj[k], ty):
            problems.append(f"{where}.{k}: {type(obj[k]).__name__}, want {ty.__name__}")
        if k in obj and ty is int and isinstance(obj[k], bool):
            problems.append(f"{where}.{k}: bool where int required")
    return problems


# The COMPLETE expected sidecar, compared exactly. Freezing only its field TYPES
# let a false canonicalization rule or a weakened quarantine banner pass green.
# Only the hash and the byte length are dynamic.
def hash_record_template(canonical_sha256, byte_length):
    return {
        "_banner": HASH_RECORD_BANNER,
        "record_format_version": "m8_8-prelock-hash-record-1",
        "status": "CANDIDATE",
        "status_note": ("NOT ISSUED. Immutability attaches at P6, after the schema "
                        "gate, the mutation battery and the cold read are green. "
                        "Until then this record may be revised freely."),
        "answer_packet_format_version": FORMAT_VERSION,
        "canonicalization_rule": ("JSON, keys sorted, two-space indent, ASCII, LF, "
                                  "single trailing newline (protocol § 10.2)"),
        "canonical_plaintext_byte_length": byte_length,
        "expected_canonical_plaintext_sha256": canonical_sha256,
    }


# Top-level container types, established BEFORE anything iterates or chains
# `.get()`. Passing None for `rows` used to raise rather than return a verdict.
TOP_LEVEL_TYPES = {
    "format_version": str, "target_id": str, "adjudicates": dict,
    "rows": list, "identities": list, "indexing_map": dict, "convention_map": dict,
}
