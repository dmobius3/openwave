"""P4: the mutation battery, in the two halves the frozen plan specifies.

HALF ONE, CONTENT. Every P3 check reddened by a targeted mutation, plus the named
cases the plan lists: unnormalized triple, decimal literal, swapped rows, broken
slot, wrong `adjudicates`, duplicated signature, self-inverse anchor substitute,
non-canonical bytes.

It also formalizes a class rather than a bug. CLAIMANT VERSUS DESIGNATED SUBJECT:
wherever a check says "the row with property X", the protocol usually says
"specific row Y must have property X", and those differ whenever a packet can
move the property. This defect appeared in P2 and then appeared AGAIN,
independently, in P3, in an implementation written from the protocol rather than
from P2's code. Two independent occurrences make it a shape the design invites,
so it gets paired cases here instead of a footnote.

HALF TWO, INGESTION. Six § 4.4 ingestion controls plus two § 5.3 and § 5.4
consumer controls, exercised against the author-side
reference comparison fixture in `p4_reference.py`.

    CLAIM SCOPE, stated exactly, because it is easy to overstate.
    These tests demonstrate that the § 4.4 INGESTION CONTROLS ASSIGNED TO P4 are
    IMPLEMENTABLE against the frozen packet: that the packet is CONSUMABLE in the
    way § 4.4 requires. The real adjudication harness is maintainer-owned, and the
    maintainer must independently bind and mutation-test the corresponding
    controls there. We are not claiming that harness is correct.

    NOT ALL of § 4.4's controls are P4's. Its controls paragraph also requires the
    canonical packet to be PUBLISHED once the adjudication record is committed,
    with its bytes verified against the frozen hash. That is issuance, it happens
    at P6, and P4 does not exercise it. An earlier wording of this scope sentence
    said "all § 4.4 controls" and was wrong by exactly that one.

EIGHT CONTROLS IN TWO CATEGORIES, and the categories are not interchangeable.

  SIX § 4.4 INGESTION CONTROLS ASSIGNED TO P4. § 4.4's controls paragraph names
  SIX obligations. P4 exercises five of them, listed below; the sixth is
  publication, which is issuance and belongs to P6. Two of the six controls split
  one obligation into its ordering half and its identity half, which is why five
  obligations give six controls:
    (1) SAME_PARSER_PATH          the bytes enter the fixture's parser, and THAT
                                  OBJECT is what the comparison consumes
    (2) WRONG_HASH_BEFORE_PARSE   hash-verified BEFORE parsing
    (3) SYNTHETIC_NONIDENTITY_MAP the map-processing path is operative
    (4) DOWNSTREAM_CELL_MUTATION  one loaded cell mutated downstream of the
                                  completed hash check reddens the comparison
    (5) NO_MANUAL_TRANSCRIPTION   no manual transcription anywhere in the path
    (6) ADJUDICATES_BINDING       `adjudicates` matches the run's actual packets
  § 4.4's sixth obligation, PUBLICATION of the canonical packet once the
  adjudication record is committed, is issuance and belongs to P6.

  TWO ADDITIONAL CONTROLS, not § 4.4 controls. These cover § 5.3 and § 5.4
  behaviour the reference consumer must have for the six above to mean anything,
  and they exist because redlines showed the fixture lacked it:
    (7) CONVENTION_MAP_VALIDATED       the consumer validates the declarations
                                       its own code relies on
    (8) ORIENTATION_SELECTION_EXERCISED  § 5.4's decision procedure on all three
                                       of its outcomes, not only the native one
  Counting eight is right internally. Saying "§ 4.4 assigns eight" would not be.

TWO CONTROLS ARE EXPLICITLY BOUNDED. Six successive measurements of parse
ordering did not hold, and five of map application, each by a route around the
instrumentation chosen. Narrowing the claims was not enough on its own: an
external redline showed both did not hold again inside their stated bounds. Both were then
rebuilt so the property holds by CONSTRUCTION rather than by observation, and the
bounds below say what that construction covers:

    WRONG_HASH_BEFORE_PARSE and SYNTHETIC_NONIDENTITY_MAP each carry one. Both
    are held once, in the `BOUNDS` dict below, and printed from there with the
    result and in the scope block. They are not restated in prose anywhere in
    this file, because four hand-maintained copies is what let a narrower
    exclusion survive in one docstring while the stronger claim stayed in the
    printed output. P4_REPORT.md and the review brief quote them, and the
    self-check compares those quotations to this output.

THIS BATTERY WRITES NOTHING. Not into the tree, not anywhere. That is the
`--emit` lesson from P2, where a battery that wrote could have left a stale green
artifact standing. Its transcript stays author-side.
"""

import ast
import copy
import hashlib
import json
import os
import pathlib
import re
import subprocess
import sys
import textwrap

# WRITES NOTHING, and that has to include what the interpreter writes on our
# behalf. The documented invocation created `build/__pycache__/*.pyc` on import,
# so the claim was false as stated even though no application code wrote a byte.
# Set before the local imports, which is the only point at which it can matter.
sys.dont_write_bytecode = True

import p3
import p3_qphi as Q
import p4_reference as REF

PACKET_SHA = "744c7f25e2312d90fc356b11510da685328f05e80ae62721d0a0f418dcf9697e"

CASES = []


def case(target, *, reason=None, require_green=(), note="", precondition=None):
    """`precondition` asserts the MUTATED INPUT really has the property the case
    is about, before any P3 rejection is counted.

    Generic reason substrings were the weak point. `reason="wrong"` is satisfied
    by either branch of a two-branch check, so a case could silently drift onto a
    different branch of the same check and still be credited. Where a detail
    string cannot separate branches at all, as with three different ways to make
    a triple invalid, only a predicate on the input can."""
    def deco(fn):
        CASES.append((fn.__name__, target, fn, reason, tuple(require_green), note,
                      precondition))
        return fn
    return deco


def _row0_value(pb):
    return json.loads(pb.decode("ascii"))["rows"][0]["value"]


# ---------------------------------------------------------------- half one

def _mutate(fn):
    """Mutate the packet, recanonicalize, and rebuild the sidecar so a packet
    mutation does not incidentally redden the sidecar checks."""
    pkt = copy.deepcopy(PACKET)
    fn(pkt)
    pb = p3.canonical_bytes(pkt)
    rec = dict(SIDECAR)
    rec["canonical_plaintext_byte_length"] = len(pb)
    rec["expected_canonical_plaintext_sha256"] = hashlib.sha256(pb).hexdigest()
    return pb, p3.canonical_bytes(rec), P1, GROUP, CONSTRUCTION


def _p1_mutate(fn):
    p1 = copy.deepcopy(P1)
    fn(p1)
    return (p3.canonical_bytes(PACKET), p3.canonical_bytes(SIDECAR), p1,
            GROUP, CONSTRUCTION)


@case("KEYS_EXACTLY_SEVEN", reason="want", note="an eighth top-level key")
def extra_top_level_key():
    return _mutate(lambda p: p.update({"notes": "extra"}))


@case("PACKET_CANONICAL_BYTES", reason="canonical serialization",
      note="NON-CANONICAL BYTES: four-space indent, parses identically")
def non_canonical_bytes():
    bad = (json.dumps(PACKET, sort_keys=True, indent=4, ensure_ascii=True) + "\n").encode()
    rec = dict(SIDECAR)
    rec["canonical_plaintext_byte_length"] = len(bad)
    rec["expected_canonical_plaintext_sha256"] = hashlib.sha256(bad).hexdigest()
    return bad, p3.canonical_bytes(rec), P1, GROUP, CONSTRUCTION


@case("SIDECAR_CANONICAL_BYTES", reason="contains CR", note="CR in the sidecar")
def sidecar_not_canonical():
    pb = p3.canonical_bytes(PACKET)
    return pb, p3.canonical_bytes(SIDECAR).replace(b"\n", b"\r\n", 1), P1, GROUP, CONSTRUCTION


@case("SIDECAR_COMMITS_TO_PACKET", reason="wrong: ['sha256']",
      note="the record naming a different hash")
def sidecar_wrong_commitment():
    pb = p3.canonical_bytes(PACKET)
    rec = dict(SIDECAR)
    h = hashlib.sha256(pb).hexdigest()
    rec["canonical_plaintext_byte_length"] = len(pb)
    rec["expected_canonical_plaintext_sha256"] = h[:-1] + ("0" if h[-1] != "0" else "1")
    return pb, p3.canonical_bytes(rec), P1, GROUP, CONSTRUCTION


@case("TRIPLES_NORMALIZED", reason="problem",
      note="UNNORMALIZED TRIPLE: scaled by 2, same number, wrong encoding",
      precondition=lambda pb, sb: (all(type(x) is int for x in _row0_value(pb)),
                                   "every member is still a true int"))
def unnormalized_triple():
    def f(p):
        a, b, c = p["rows"][0]["value"]
        p["rows"][0]["value"] = [a * 2, b * 2, c * 2]
    return _mutate(f)


@case("TRIPLES_NORMALIZED", reason="problem",
      note="DECIMAL LITERAL: § 4.4 admits no decimal anywhere",
      precondition=lambda pb, sb: (any(isinstance(x, float) for x in _row0_value(pb)),
                                   "a float is present in the mutated triple"))
def decimal_literal():
    def f(p):
        p["rows"][0]["value"] = [float(x) for x in p["rows"][0]["value"]]
    return _mutate(f)


@case("TRIPLES_NORMALIZED", reason="problem",
      note="JSON true where an integer belongs; Python says True == 1",
      precondition=lambda pb, sb: (any(isinstance(x, bool) for x in _row0_value(pb)),
                                   "a bool is present in the mutated triple"))
def boolean_leaf():
    def f(p):
        p["rows"][0]["value"][1] = True
    return _mutate(f)


# `reason="{"` is generic, and deliberately left so: CLASS_CENSUS has a single
# failing branch that renders the observed census, so there is no sibling branch
# a substring could select by mistake. SIGNATURES_DISTINCT below is the opposite
# case and was changed, because it has two branches that both render as
# "N distinct of M" and its substring matched either.
@case("CLASS_CENSUS", reason="{", note="the 1/7/1 census broken")
def census_broken():
    def f(p):
        for r in p["rows"]:
            if r["evidentiary_class"] == p3.CLASS_FREE:
                r["evidentiary_class"] = p3.CLASS_DECLARED
                return
    return _mutate(f)


@case("ADJUDICATES_MATCH_PUBLIC_PACKETS", reason="group hash is not § 11's pin",
      note="WRONG ADJUDICATES: one hex digit of the group hash")
def wrong_adjudicates():
    def f(p):
        h = p["adjudicates"]["group_packet_sha256"]
        p["adjudicates"]["group_packet_sha256"] = h[:-1] + ("0" if h[-1] != "0" else "1")
    return _mutate(f)


@case("ADJUDICATES_MATCH_PUBLIC_PACKETS", reason="supplied group packet",
      note="a substituted public packet WITH `adjudicates` updated to match it")
def paired_public_packet_substitution():
    gb = GROUP + b" "
    pkt = copy.deepcopy(PACKET)
    pkt["adjudicates"]["group_packet_sha256"] = hashlib.sha256(gb).hexdigest()
    pb = p3.canonical_bytes(pkt)
    rec = dict(SIDECAR)
    rec["canonical_plaintext_byte_length"] = len(pb)
    rec["expected_canonical_plaintext_sha256"] = hashlib.sha256(pb).hexdigest()
    return pb, p3.canonical_bytes(rec), P1, gb, CONSTRUCTION


@case("ADJUDICATES_MATCH_PUBLIC_PACKETS", reason="construction hash is not the one",
      note="the CONSTRUCTION hash wrong, with the group hash left correct")
def wrong_adjudicates_construction():
    """The other half of the same check, and it was never exercised.

    Both existing cases touch only the group packet, so either construction
    branch could be disabled with P4 still at 20/20 and 12/12. `control_6` had
    already recorded fixing exactly this on the ingestion side; the repair never
    reached the content half, which is why the reason substring here names the
    construction branch specifically rather than matching any rejection."""
    def f(p):
        h = p["adjudicates"]["construction_packet_sha256"]
        p["adjudicates"]["construction_packet_sha256"] = \
            h[:-1] + ("0" if h[-1] != "0" else "1")
    return _mutate(f)


@case("ADJUDICATES_MATCH_PUBLIC_PACKETS",
      reason="supplied construction packet does not hash",
      note="a substituted CONSTRUCTION packet, `adjudicates` left correct")
def substituted_construction_packet():
    return (p3.canonical_bytes(PACKET), p3.canonical_bytes(SIDECAR), P1,
            GROUP, CONSTRUCTION + b" ")


@case("SIGNATURES_DISTINCT", reason="8 distinct of 9", require_green=("CLASS_CENSUS",),
      note="DUPLICATED SIGNATURE: § 5.5 must separate all nine")
def duplicated_signature():
    def f(p):
        p["rows"][1]["signature"] = copy.deepcopy(p["rows"][0]["signature"])
    return _mutate(f)


@case("VALUES_MATCH_BOTH_SOURCE_FAMILIES", reason="want both families",
      note="a whole FAMILY removed: only theory renderings left")
def one_family_missing():
    """The branch the check is NAMED after, and it was not exercised.

    Three cases targeted this check and all three altered a stated form, so the
    `{theory, platform} <= fams` branch could be disabled with P4 fully green:
    the check's own headline property was the one going untested."""
    def f(p1):
        for v in p1["values"]:
            if v["label"] == "R5":
                for k in [k for k in v["as_written"] if k.startswith("B")]:
                    del v["as_written"][k]
                v["families"] = ["theory"]
    return _p1_mutate(f)


@case("VALUES_MATCH_BOTH_SOURCE_FAMILIES", reason="source keys say",
      note="`families` desynchronised from the source keys it summarises")
def families_field_desynchronised():
    """The cross-check that keeps P1's summary honest about its own record."""
    def f(p1):
        for v in p1["values"]:
            if v["label"] == "R5":
                v["families"] = ["theory"]
    return _p1_mutate(f)


@case("VALUES_MATCH_BOTH_SOURCE_FAMILIES", reason="unparseable stated form",
      note="a stated closed form in notation the grammar does not model")
def unparseable_stated_form():
    """P3 calls this a deliberate design decision and nothing exercised it.

    `p3_forms` is a restricted grammar over the syntaxes the frozen sources
    actually use, and its comment says an unparseable form is A STOP, NOT A
    SKIP. A decision that consequential should have a case; without one, the
    `except` arm could be turned into `continue` and P4 would not notice."""
    def f(p1):
        for v in p1["values"]:
            if v["label"] == "R5":
                k = next(iter(v["as_written"]))
                # The numerator here is the established DECOY, one off the
                # committed one, as the two family cases above already use. The
                # first version of this case wrote the row's real value into a
                # target-free source and the leak scan caught it; naming that
                # value in this comment would have done the same thing again.
                v["as_written"][k] = "sqrt[5]{26/9}"
    return _p1_mutate(f)


@case("PACKET_CANONICAL_BYTES", reason="not ASCII",
      note="§ 10.2's ASCII requirement, a non-ASCII byte inside a declaration")
def packet_not_ascii():
    """Injected at the BYTE level, because the canonicalizer cannot emit this.

    `canonical_bytes` serializes with `ensure_ascii`, so a non-ASCII character
    placed in a value comes back escaped and the bytes stay ASCII: the first
    version of this case mutated the object and never reddened anything. The
    requirement is on the BYTES, so the case has to be too."""
    b = p3.canonical_bytes(PACKET)
    assert b.count(b"admitted") == 1
    return (b.replace(b"admitted", b"admitt\xc3\xa9d"),
            p3.canonical_bytes(SIDECAR), P1, GROUP, CONSTRUCTION)


@case("PACKET_CANONICAL_BYTES", reason="not exactly one trailing newline",
      note="§ 10.2's single-trailing-newline requirement")
def packet_double_trailing_newline():
    return (p3.canonical_bytes(PACKET) + b"\n", p3.canonical_bytes(SIDECAR),
            P1, GROUP, CONSTRUCTION)


@case("TRIPLES_NORMALIZED", reason="dimension: not a true int",
      note="a signature DIMENSION written as JSON true, the bool-is-an-int trap")
def signature_dimension_boolean():
    """`True == 1`, so a dimension written as `true` indexes and compares as 1.

    The same trap the value fields already have a case for, at the signature
    field, which had none."""
    def f(p):
        p["rows"][0]["signature"]["dimension"] = True
    return _mutate(f)


@case("DECLARED_ROW_CONVENTION", reason="value",
      note="R0 no longer declaring T^2 = 1")
def declared_row_value_wrong():
    def f(p):
        for r in p["rows"]:
            if r["label"] == p3.DECLARED_ROW:
                r["value"] = [2, 0, 1]
    return _mutate(f)


@case("DECLARED_ROW_CONVENTION", reason="class",
      require_green=("CLASS_CENSUS",),
      note="CLAIMANT VS DESIGNATED: declared class moved off R0, census preserved")
def declared_class_moved_off_designated_row():
    def f(p):
        a = next(r for r in p["rows"] if r["label"] == p3.DECLARED_ROW)
        b = next(r for r in p["rows"] if r["evidentiary_class"] == p3.CLASS_FREE)
        a["evidentiary_class"], b["evidentiary_class"] = \
            b["evidentiary_class"], a["evidentiary_class"]
    return _mutate(f)


@case("IDENTITIES_RECOMPUTE_POSITION_WISE", reason="recompute",
      note="BROKEN SLOT: one expected value taken from another slot")
def broken_identity_slot():
    def f(p):
        p["identities"][0]["expected_value"] = list(p["identities"][2]["expected_value"])
    return _mutate(f)


@case("IDENTITIES_RECOMPUTE_POSITION_WISE", reason="position",
      note="a slot declaring a position that is not its own index")
def identity_position_wrong():
    def f(p):
        p["identities"][2]["position"] = 7
    return _mutate(f)


@case("VALUES_MATCH_BOTH_SOURCE_FAMILIES", reason="differs",
      note="SWAPPED ROWS: two values exchanged, both still normalized")
def swapped_rows():
    def f(p):
        p["rows"][0]["value"], p["rows"][8]["value"] = \
            p["rows"][8]["value"], p["rows"][0]["value"]
    return _mutate(f)


@case("VALUES_MATCH_BOTH_SOURCE_FAMILIES", reason="/A: stated form",
      note="the THEORY family's stated closed form altered")
def theory_family_form_altered():
    def f(p1):
        for v in p1["values"]:
            if v["label"] == "R5":
                k = next(k for k in v["as_written"] if k.startswith("A "))
                v["as_written"][k] = "mp.mpf(26) / 9"
    return _p1_mutate(f)


@case("VALUES_MATCH_BOTH_SOURCE_FAMILIES", reason="/B1: stated form",
      note="the PLATFORM family's stated closed forms altered")
def platform_family_form_altered():
    def f(p1):
        for v in p1["values"]:
            if v["label"] == "R5":
                for k in list(v["as_written"]):
                    if k.startswith("B"):
                        v["as_written"][k] = "26/9"
    return _p1_mutate(f)


@case("SELECTOR_NOT_SELF_INVERSE", reason="SELF-INVERSE",
      note="SELF-INVERSE ANCHOR SUBSTITUTE: R7 set to 1")
def self_inverse_anchor():
    def f(p):
        for r in p["rows"]:
            if r["label"] == p3.SELECTOR_ROW:
                r["value"] = [1, 0, 1]
    return _mutate(f)


@case("SELECTOR_NOT_SELF_INVERSE", reason="carries",
      require_green=("CLASS_CENSUS",),
      note="CLAIMANT VS DESIGNATED: selector class moved off R7, census preserved")
def selector_class_moved_off_designated_row():
    def f(p):
        a = next(r for r in p["rows"] if r["label"] == p3.SELECTOR_ROW)
        b = next(r for r in p["rows"] if r["evidentiary_class"] == p3.CLASS_FREE)
        a["evidentiary_class"], b["evidentiary_class"] = \
            b["evidentiary_class"], a["evidentiary_class"]
    return _mutate(f)


# ---------------------------------------------------------------- half two

class ParseSentinel:
    """Counts entries into the parser, retains what it returned, AND what it ate.

    Retaining the result matters: `SAME_PARSER_PATH` used to prove only that the
    sentinel was called, not that the object handed onward CAME FROM that call. A
    loader could parse once through the sentinel, discard it, and return a second
    parser's result. Identity is checked with `is`.

    Retaining the INPUT matters for the other half of the same claim, and it was
    missing. The control is mapped to "loads the canonical bytes directly", but
    nothing compared the parser's input to the bytes that were hashed. A loader
    could hash the canonical file, re-render it, and parse the re-rendering:
    `return parse(raw.decode("ascii") + " ")` kept object identity, kept the
    control green, and meant the object compared was parsed from something other
    than the hashed bytes.
    """

    def __init__(self):
        self.calls = 0
        self.result = None
        self.inputs = []

    def __call__(self, payload):
        self.calls += 1
        self.inputs.append(payload)
        self.result = _WatchedPacket(json.loads(payload))
        return self.result


# Type-correct stand-ins for the five top-level keys the comparison needs. A
# consumer starved of these fails cleanly rather than crashing, so the control
# gets a verdict instead of a traceback.
_HOLLOW = {"rows": [], "identities": [], "indexing_map": {"entries": []},
           "convention_map": {}, "adjudicates": {}}


class _WatchedPacket(dict):
    """The parsed packet, whose real content is reachable ONLY through itself.

    `pkt is sentinel.result` is checked before `compare` runs and says nothing
    about what `compare` then does. Rebinding `packet = dict(packet)` on entry
    gave the comparison a different object while the control still printed that
    the object compared IS the object the parser returned. Counting reads on this
    object caught that, and did not hold either, once reading the five keys and
    THEN copying: a dead warm-up read satisfied the counter while every
    substantive access went to the copy. Counting proves touching, not use.

    So the content is not in the dict. The five keys § 5.3 and § 5.4 need are
    held in a side table and the concrete storage carries `_HOLLOW` stand-ins;
    `__getitem__` and `get` serve the real values from the side table. CPython's
    `dict(other)` copies through the concrete storage without calling
    `__getitem__`, so a shallow copy inherits the stand-ins and any consumer
    working from it compares nothing. Warm-up reads no longer help, because what
    they warm is not where the copy looks.

    The criterion is therefore causal rather than observational, as it is for
    `_MapOnlyRows`: a comparison that AGREES has taken its data through this
    object, whatever route it took, because no other route carries the data."""

    def __init__(self, real):
        super().__init__({k: (list(v) if isinstance(v, list) else dict(v))
                          for k, v in _HOLLOW.items()} | 
                         {k: v for k, v in real.items() if k not in _HOLLOW})
        self._real = real
        self.read = set()

    def __getitem__(self, key):
        if key in self._real:
            self.read.add(key)
            return self._real[key]
        return super().__getitem__(key)

    def get(self, key, default=None):
        if key in self._real:
            self.read.add(key)
            return self._real[key]
        return super().get(key, default)


# The two bounded claims, held ONCE, in code. Every place that states a bound
# reads it from here: the two controls' result lines and the scope block. The
# report and the review brief quote them, and the self-check compares those
# quotations to this output.
#
# They used to be prose, restated in a module docstring, a function docstring, a
# printed line and a report. An external redline found what that invites: the
# function docstring carried an exclusion the printed bound did not, so the
# implementation had been narrowed while the stronger claim survived in prose,
# and the self-check compared the module docstring and missed it. Drift between
# four hand-maintained copies is not a thing to check for. It is a thing to make
# impossible, so no docstring below restates a bound; they refer to this dict.
BOUNDS = {
    "WRONG_HASH_BEFORE_PARSE":
        "demonstrates ordering for THIS loader, structurally: verification runs "
        "where no JSON implementation can be imported, and the loader's only "
        "parse call takes its result as its sole argument. It does not constrain "
        "a consumer that does not use this loader",
    "SYNTHETIC_NONIDENTITY_MAP":
        "demonstrates that, in the reference fixture, all committed-row "
        "addressing used by comparison is map-derived; it is an implementation "
        "test, not a proof against a consumer written to misreport its own behaviour",
}


def control_1_same_parser_path():
    """The canonical bytes enter the fixture's parser, and THAT OBJECT is what
    the comparison consumes. The identity check is the part that was missing."""
    canonical = pathlib.Path(PACKET_PATH).read_bytes()
    sentinel = ParseSentinel()
    pkt = REF.load_hash_verified(PACKET_PATH, PACKET_SHA, parse=sentinel)
    same = pkt is sentinel.result
    # The bytes that were hashed, compared to what the parser was actually given.
    fed = sentinel.inputs[:1]
    verbatim = fed == [canonical]
    ok, why = REF.compare(pkt, COMMITTED, GROUP, CONSTRUCTION)
    # Read ON THE PARSED OBJECT during the comparison, not merely present in it.
    # And again: `not missed` is vacuous on an empty inventory, which printed
    # "0/0 needed keys" and stayed green. The count is frozen here.
    NEEDED = {"rows", "identities", "indexing_map", "convention_map", "adjudicates"}
    NEEDED_COUNT = 5
    missed = sorted(NEEDED - pkt.read)
    used = same and not missed and len(NEEDED) == NEEDED_COUNT
    return (sentinel.calls == 1 and used and verbatim and ok,
            f"parser entered {sentinel.calls} time(s) and was fed the hashed bytes "
            f"{'VERBATIM' if verbatim else 'ALTERED: ' + repr(fed)[:60]}; the object "
            f"compared {'IS' if same else 'IS NOT'} the object that parser returned, "
            f"and its content was reachable ONLY through it: {len(NEEDED & pkt.read)}"
            f"/{len(NEEDED)} needed keys served from that object, a shallow copy "
            f"carrying stand-ins"
            + ("" if not missed else f" -- NOT served from it: {missed}, so something "
               "else was compared")
            + f"; comparison {'agrees' if ok else 'FAILED: ' + why}")


_NO_JSON_PROBE = r'''
import sys

class _RefuseJson:
    """Refuse `json` and its C accelerator at import, before anything binds one."""
    def find_module(self, name, path=None):
        return self if name.split(".")[0] in ("json", "_json") else None
    def find_spec(self, name, path=None, target=None):
        if name.split(".")[0] in ("json", "_json"):
            raise ImportError(f"{name} is unavailable in this interpreter")
        return None

for _m in [m for m in sys.modules if m.split(".")[0] in ("json", "_json")]:
    del sys.modules[_m]
sys.meta_path.insert(0, _RefuseJson())

# POSITIVE CONTROL. If this import succeeds the import guard is not working and
# everything below would be worthless.
try:
    import json
    print("BLOCKER-FAILED")
    raise SystemExit(0)
except ImportError:
    pass

sys.path.insert(0, sys.argv[1])
import p4_verify

# Nothing in the verifier's namespace came from a JSON implementation.
leaked = sorted(k for k, v in vars(p4_verify).items()
                if getattr(getattr(v, "__module__", "") or "", "split", None)
                and str(getattr(v, "__module__", "")).split(".")[0] in ("json", "_json"))

raw = p4_verify.verified_bytes(sys.argv[2], sys.argv[3])
good = len(raw)
try:
    p4_verify.verified_bytes(sys.argv[2], "0" * 64)
    refused = False
except p4_verify.HashMismatch:
    refused = True
print(f"OK {good} {refused} {leaked}")
'''


def control_2_wrong_hash_before_parse():
    """A wrong hash fails before any parsing, and that is now STRUCTURAL.

    SIX narrower measurements failed before this one, each looking total until
    the next route was found: an exit code; an injectable sentinel counting only
    the parser passed in; a patch of `json.loads`, which `JSONDecoder().decode`
    does not pass through; a patch of `raw_decode`, which `scan_once` does not pass through; those plus
    the scanner FACTORY, which the module's own pre-built `_default_decoder`
    does not pass through; and those plus that one instance, which ANY OTHER decoder built
    before the wrappers go on does not pass through. `EARLY_DECODER = json.JSONDecoder()` at
    module scope, called before the hash check, kept the control printing zero.

    The sixth failure is the one that settles the method. `scan_once` is bound
    when a decoder is constructed, so no later patch can reach a decoder that
    already exists, and an enumeration of parser callables will always admit one
    more. Measurement is the wrong instrument for this property.

    So the ordering is not observed here. It is a consequence of structure, in
    FOUR parts, each of which can fail. The fourth was added after a reviewer
    defeated all three of the others, and this sentence still said "three" a
    round later: a docstring understating its own implementation is the drift the
    bound machinery exists to prevent, at a site the paraphrase check cannot see,
    since that check keys on the opening of a bound and this quotes none.

      ISOLATION.  `p4_verify` imports `hashlib` and `pathlib` and nothing else,
                  checked on its syntax tree rather than by reading it.
      EXECUTION.  That module verifies the real bytes, and refuses a wrong hash,
                  in a subprocess where `json` and `_json` cannot be imported at
                  all. A positive control proves the import guard works: `import json`
                  must raise there, or the whole probe is discarded. Verification
                  completing in an interpreter with no JSON implementation is
                  what rules out every pre-bound instance at once, since none can
                  exist.
      COMPOSITION. `load_hash_verified` is one `return parse(verified_bytes(...))`
                  and the syntax tree is required to be exactly that: one
                  statement, the parse call taking the verification call as its
                  only argument. A parse before the hash check would need a
                  statement to live in, and there is nowhere to put one.

    The claim is bounded, and the bound is BOUNDS["WRONG_HASH_BEFORE_PARSE"],
    printed with every result. It is deliberately not restated here. A previous
    version of this docstring carried an exclusion the printed bound did not,
    which is exactly the drift the review brief warns about, and it survived
    because the self-check read the module docstring and not this one.

    The old enumeration is kept below as corroboration, with its negative
    control, because a measurement that agrees with a structural argument is
    worth having. It is no longer what the claim rests on.
    """
    import json.decoder
    import json.scanner

    # ---- ISOLATION: the verifier's imports AND its body, from its syntax tree
    #
    # Imports alone were not enough. Python's built-ins need no import, and
    # `eval(raw.decode("ascii").replace("true", "True")...)` parses the packet
    # with nothing imported at all, runs happily in the no-JSON subprocess, and
    # left this limb green. "Cannot reach a JSON module" is a weaker property
    # than "does nothing but hash and compare", so the body is pinned too: a
    # fixed sequence of statements, each shape-checked, with no call permitted
    # anywhere except the four the verification actually needs.
    # THE SOURCE PATH IS FIXED HERE, not asked of the object under test.
    # `REF.p4_verify.__file__` let the object nominate its own evidence: binding
    # the name `p4_verify` to a delegate that pre-parses and then calls the real
    # verifier, carrying the real `__file__`, made this limb read the honest file
    # while a different object did the work. That is the FIFTH occurrence in this
    # build of a check taking its expectation from what it is checking.
    here = pathlib.Path(__file__).resolve().parent
    vpath = here / "p4_verify.py"
    vsrc = vpath.read_text()
    vtree = ast.parse(vsrc)
    imported = set()
    for node in ast.walk(vtree):
        if isinstance(node, ast.Import):
            imported.update(a.name.split(".")[0] for a in node.names)
        elif isinstance(node, ast.ImportFrom):
            imported.add((node.module or "").split(".")[0])

    # FIVE, not six. `isinstance` was here only to admit `verified_bytes`'s
    # second entry shape; with that arm gone the allowlist is a prohibition
    # rather than a permission, and an `isinstance` appearing in the verifier
    # now reddens this limb.
    ALLOWED_CALLS = {"read_bytes", "Path", "sha256", "hexdigest", "HashMismatch"}
    # The FIRST `def verified_bytes` in the file, which is the one whose body the
    # shape check below pins. The identity limb requires the RUNNING function to
    # be this one, not a later rebinding of the name.
    vfn = vfn_def = next((n for n in ast.walk(vtree)
                          if isinstance(n, ast.FunctionDef)
                          and n.name == "verified_bytes"), None)
    vbody = [s for s in (vfn.body if vfn else [])
             if not (isinstance(s, ast.Expr) and isinstance(s.value, ast.Constant))]
    called = set()
    for node in ast.walk(vfn) if vfn else ():
        if isinstance(node, ast.Call):
            f = node.func
            called.add(f.attr if isinstance(f, ast.Attribute)
                       else f.id if isinstance(f, ast.Name) else "?")
    body_shape = (len(vbody) == 4
                  and isinstance(vbody[0], ast.Assign)       # raw = ...read_bytes()
                  and isinstance(vbody[1], ast.Assign)       # got = sha256(raw)...
                  and isinstance(vbody[2], ast.If)           # if got != expected: raise
                  and len(vbody[2].body) == 1
                  and isinstance(vbody[2].body[0], ast.Raise)
                  and not vbody[2].orelse
                  and isinstance(vbody[3], ast.Return)       # return raw
                  and isinstance(vbody[3].value, ast.Name)
                  and vbody[3].value.id == vbody[0].targets[0].id)
    isolated = (imported == {"hashlib", "pathlib"} and body_shape
                and called <= ALLOWED_CALLS)

    # ---- EXECUTION: verification where no JSON implementation exists --------
    # `-B`, because this battery claims to write nothing and that includes what
    # the interpreter writes on its behalf. Without it the subprocess imports
    # `p4_verify` and leaves `build/__pycache__/p4_verify.*.pyc` behind, and the
    # claim would have been false again in a new way.
    probe = subprocess.run(
        [sys.executable, "-B", "-c", _NO_JSON_PROBE,
         str(pathlib.Path(__file__).parent), PACKET_PATH, PACKET_SHA],
        capture_output=True, text=True, timeout=60,
        env=dict(os.environ, PYTHONDONTWRITEBYTECODE="1"))
    line = probe.stdout.strip().splitlines()[-1] if probe.stdout.strip() else ""
    if line == "BLOCKER-FAILED":
        return False, "the no-JSON subprocess check could still import json; its result proves nothing"
    if not line.startswith("OK "):
        return False, f"the no-JSON subprocess check did not complete: {(probe.stderr or probe.stdout)[-160:]}"
    _, nbytes, refused_there, leaked = line.split(" ", 3)
    no_json_ok = (int(nbytes) == len(pathlib.Path(PACKET_PATH).read_bytes())
                  and refused_there == "True" and leaked == "[]")

    # ---- COMPOSITION: the loader is one expression -------------------------
    lpath = here / "p4_reference.py"
    lsrc = lpath.read_text()
    loader = next((n for n in ast.walk(ast.parse(lsrc))
                   if isinstance(n, ast.FunctionDef) and n.name == "load_hash_verified"),
                  None)
    body = [s for s in (loader.body if loader else [])
            if not (isinstance(s, ast.Expr) and isinstance(s.value, ast.Constant))]

    def bare(node, name):
        return isinstance(node, ast.Name) and node.id == name

    # Pinning the SHAPE was not enough. Requiring an outer `parse(...)` around an
    # attribute call named `verified_bytes` left the inner call's ARGUMENTS
    # unconstrained, and Python evaluates those before the call runs, so a parse
    # placed in an argument expression happened strictly before the hash
    # check while this test still passed. "No statement to live in" is only true
    # if there is no expression to live in either. Both arguments must therefore
    # be bare names, which admits no evaluation of anything.
    composed = False
    if len(body) == 1 and isinstance(body[0], ast.Return) and isinstance(body[0].value, ast.Call):
        outer = body[0].value
        inner = outer.args[0] if len(outer.args) == 1 and not outer.keywords else None
        composed = (bare(outer.func, "parse")
                    and isinstance(inner, ast.Call)
                    and isinstance(inner.func, ast.Attribute)
                    and bare(inner.func.value, "p4_verify")
                    and inner.func.attr == "verified_bytes"
                    and not inner.keywords
                    and len(inner.args) == 2
                    and bare(inner.args[0], "path")
                    and bare(inner.args[1], "expected_sha256"))

    # ---- IDENTITY: the names resolve to the objects that were inspected ----
    #
    # Everything above is a fact about two source FILES. None of it constrains
    # what `p4_verify` and `load_hash_verified` are BOUND TO at runtime, and two
    # routes used exactly that gap: a delegate object under the name, and a
    # decorator that leaves the body byte-identical and runs before it. A
    # structural argument about text is worth nothing if the text is not what
    # executes.
    vmod = sys.modules.get("p4_verify")
    vfn = getattr(REF, "p4_verify", None)
    loader_fn = getattr(REF, "load_hash_verified", None)
    identity = []
    if vmod is None or vfn is not vmod:
        identity.append("`p4_verify` in the fixture is not the imported module")
    elif pathlib.Path(getattr(vmod.verified_bytes, "__code__",
                              vmod).co_filename).resolve() != vpath:
        identity.append("`verified_bytes` was not defined in the file inspected")
    # A FILE IS NOT A FUNCTION. Pinning `co_filename` alone left a second `def`
    # in the SAME file, bound to the name after the honest one, satisfying every
    # limb at once: the AST still found and pinned the honest def, `eval` needs
    # no import so the import set held, `called` walked only the honest def, and
    # the subprocess ran the impostor happily with no JSON implementation
    # anywhere. The loader beside this already pinned its first line; the
    # verifier did not, and nothing made "which function do I pin" and "which
    # function runs" the same question.
    elif vmod.verified_bytes.__code__.co_firstlineno != (
            vfn_def.lineno if vfn_def else -1):
        identity.append(
            f"`verified_bytes` at runtime starts at line "
            f"{vmod.verified_bytes.__code__.co_firstlineno}, but the `def` whose "
            f"body was pinned is at {vfn_def.lineno if vfn_def else '?'}; a "
            "second definition rebound to the name does that")
    if loader_fn is None or getattr(loader_fn, "__wrapped__", None) is not None:
        identity.append("the loader is wrapped, so something runs before its body")
    elif loader_fn.__code__.co_firstlineno != (loader.lineno if loader else -1):
        identity.append(
            f"the loader at runtime starts at line "
            f"{loader_fn.__code__.co_firstlineno}, the `def` inspected is at "
            f"{loader.lineno if loader else '?'}; a decorator does that")
    elif pathlib.Path(loader_fn.__code__.co_filename).resolve() != lpath:
        identity.append("the loader was not defined in the file inspected")
    # THE GLOBALS THE PINNED BODY CALLS INTO. Three bindings were pinned and the
    # two the body actually resolves were not: replacing `p4_verify.pathlib`
    # with a shim whose `Path(...).read_bytes()` parses before returning left the
    # file on disk untouched, the running function the honest one at the honest
    # line, the loader the honest composition, and the subprocess importing the
    # file fresh with no JSON present. A pre-built decoder is invisible to the
    # enumeration, so the measurement still read zero. Pinning a function is not
    # pinning what it calls.
    if vmod is not None:
        g = getattr(vmod.verified_bytes, "__globals__", {})
        if g.get("pathlib") is not pathlib or g.get("hashlib") is not hashlib:
            identity.append("`verified_bytes` resolves `pathlib`/`hashlib` to "
                            "something other than the stdlib modules inspected")

    # ---- CORROBORATION: the old enumerated measurement ---------------------
    counted = {"n": 0}

    def wrap(fn):
        def counting(*a, **k):
            counted["n"] += 1
            return fn(*a, **k)
        return counting

    def wrap_factory(fn):
        def counting_factory(ctx):
            return wrap(fn(ctx))
        return counting_factory

    targets = [(json, "loads", wrap),
               (json.JSONDecoder, "raw_decode", wrap),
               (json.JSONDecoder, "decode", wrap),
               (json.decoder.scanner, "make_scanner", wrap_factory),
               (json.scanner, "py_make_scanner", wrap_factory),
               (json._default_decoder, "scan_once", wrap)]
    if hasattr(json.scanner, "c_make_scanner"):
        targets.append((json.scanner, "c_make_scanner", wrap_factory))
    saved = [(o, a, getattr(o, a)) for o, a, _ in targets if hasattr(o, a)]

    wrong = PACKET_SHA[:-1] + ("0" if PACKET_SHA[-1] != "0" else "1")
    for obj, attr, how in targets:
        if hasattr(obj, attr):
            setattr(obj, attr, how(getattr(obj, attr)))
    try:
        try:
            REF.load_hash_verified(PACKET_PATH, wrong)
            return False, "the wrong hash was ACCEPTED"
        except REF.HashMismatch:
            production = counted["n"]
        counted["n"] = 0
        raw = pathlib.Path(PACKET_PATH).read_bytes()
        json._default_decoder.scan_once(raw.decode("ascii"), 0)
        detector = counted["n"]
    finally:
        for obj, attr, real in saved:
            setattr(obj, attr, real)

    ok = (isolated and no_json_ok and composed and not identity
          and production == 0 and detector > 0)
    return (ok,
            f"STRUCTURAL: the verifier imports {sorted(imported)}"
            f"{'' if imported == {'hashlib', 'pathlib'} else ' AND SHOULD NOT'}, "
            f"has the pinned four-statement body (read, hash, raise on mismatch, "
            f"return the same bytes)"
            f"{'' if body_shape else ' -- NO IT DOES NOT'}, and calls only "
            f"{sorted(called)}"
            f"{'' if called <= ALLOWED_CALLS else ' INCLUDING SOMETHING UNPERMITTED'}; it "
            f"verified the packet and refused a wrong hash in a subprocess where "
            f"`json` cannot be imported at all"
            f"{'' if no_json_ok else ' -- BUT THE PROBE DISAGREES: ' + line}"
            f", with a positive control confirming the import guard; and the loader is "
            f"{'one `return parse(verified_bytes(...))` and nothing else'
             if composed else 'NOT the single composition this depends on'}; "
            f"and the names resolve to the objects inspected"
            f"{'' if not identity else ' -- NO: ' + '; '.join(identity)}. "
            f"CORROBORATED by measurement: {production} parse(s) across "
            f"{len(saved)} enumerated stdlib entry points, negative control seen "
            f"({detector}).")


class _TaggedDestination:
    """A destination that records being used as an index into the committed rows.

    Not an `int` subclass on purpose: CPython indexes a list with an int subclass
    without calling `__index__`, so a subclass could be used without the tag
    firing."""

    __slots__ = ("_v", "_c")

    def __init__(self, value, counter):
        self._v, self._c = value, counter

    def __index__(self):
        self._c["map"] += 1
        return self._v

    def __eq__(self, o):
        return self._v == (o._v if isinstance(o, _TaggedDestination) else o)

    def __hash__(self):
        return hash(self._v)

    def __lt__(self, o):
        return self._v < (o._v if isinstance(o, _TaggedDestination) else o)

    def __le__(self, o):
        return self._v <= (o._v if isinstance(o, _TaggedDestination) else o)

    def __gt__(self, o):
        return self._v > (o._v if isinstance(o, _TaggedDestination) else o)

    def __ge__(self, o):
        return self._v >= (o._v if isinstance(o, _TaggedDestination) else o)

    def __repr__(self):
        return f"dst({self._v})"


# A value no row holds, well-formed so that reaching it is a DISAGREEMENT rather
# than a crash: normalized, invertible, and nowhere in the packet.
_STAND_IN = [999983, 0, 1]


class _MapOnlyRows:
    """Committed rows that yield their real VALUE only through a map destination.

    Five versions of this control did not hold, each by doing the instrumented
    thing and then not using it. Counting index USE failed against a consumer
    that indexed and discarded. Counting plain versus tagged access failed
    against one that did not index at all: `list.__iter__` never calls `__getitem__`, so
    `next(r for j, r in enumerate(rows) if j == join[...])` reached every row
    while the counters still read 1 tagged, 0 plain and the control passed.

    Counting is the wrong instrument. Each successive measurement was a claim
    about what a consumer HAPPENED to do, and there is always one more way to do
    the same thing differently. This is a claim about what a consumer CAN reach:
    the real values are served only in response to a tagged destination, and
    every other route, plain index or iteration, is served `_STAND_IN` instead. A
    comparison that agrees has therefore read the real values through the map,
    whatever route it took to get there, because no other route carries them.

    Not a `list` subclass, deliberately: `list.__iter__` and the other C-level
    paths would shortcut the wrapper, which is exactly the defect this replaces.

    The join is unaffected and needs to be. It reads signatures, never values,
    and it legitimately walks the whole table before any orientation is selected;
    masked rows carry their real signatures for that reason."""

    def __init__(self, rows, counter):
        self._rows = [dict(r) for r in rows]
        self._c = counter

    def _stand_in(self, i):
        r = dict(self._rows[i])
        r["value"] = list(_STAND_IN)
        return r

    def __len__(self):
        return len(self._rows)

    def __iter__(self):
        self._c["iter"] += 1
        return iter([self._stand_in(i) for i in range(len(self._rows))])

    def __getitem__(self, i):
        if isinstance(i, _TaggedDestination):
            self._c["tagged"] += 1
            return self._rows[int(i)]
        self._c["plain"] += 1
        return self._stand_in(int(i))


def control_3_synthetic_nonidentity_map():
    """A SYNTHETIC nonidentity map. What is demonstrable here, and what is not.

    § 4.4 asks a synthetic nonidentity fixture to prove the map-processing path
    is operative, "since a live map may be too simple to discriminate applying
    from ignoring it". That control assumes a consumer for which the map IS the
    alignment mechanism. This fixture also joins by signature, per § 5.5, and
    validation requires the two to agree, so ignoring the map yields the same
    ANSWER. No comparison of values can separate them, and four instrumentations
    each failed to hold against a consumer that performed the instrumented action and
    discarded its result.

    So this no longer measures what the consumer does. It changes what the
    consumer can REACH. The committed table is served by `_MapOnlyRows`, which
    returns real values in response to a tagged map destination and `_STAND_IN` to
    every other route. Map mode must then AGREE, which it can only do by reading
    through the map, and join mode must FAIL, which is what proves the mask is
    live rather than a wrapper that never fired. Agreement in map mode without
    that paired failure would mean nothing.

    The bound is BOUNDS["SYNTHETIC_NONIDENTITY_MAP"], printed with the result
    and not restated here. In particular a sufficiently determined consumer could
    construct `_TaggedDestination` objects of its own; the fixture is our own
    code, and the question is whether it consumes the map.
    """
    n = 4
    sigs = [{"dimension": i + 1, "s": [i, 0, 1], "t": [0, i, 1], "st": [1, i, 1]}
            for i in range(n)]
    # S2 and S3 deliberately COMMIT THE SAME VALUE. Every earlier synthetic had
    # all-distinct values, so any malformed map disagreed on values and the
    # "malformed map rejected" subcase passed without `check_indexing_map`
    # existing at all: removing that function outright left P4 green, because the
    # wrong map merely produced a DISAGREEMENT further downstream. With a
    # duplicated value, swapping those two destinations is invisible to
    # comparison, so only the validator can catch it.
    vals = [[i + 2, 0, 1] for i in range(n)]
    vals[3] = list(vals[2])
    # The anchor is S0, so S0 carries the class § 5.4 designates. A synthetic
    # fixture that declared an anchor no row owned would now be rejected by the
    # consumer, and rightly: that is the defect this synthetic must not model.
    rows = [{"label": f"S{i}",
             "evidentiary_class": p3.CLASS_SELECTOR if i == 0 else p3.CLASS_FREE,
             "value": list(vals[i]), "signature": sigs[i]} for i in range(n)]
    perm = {i: (i + 1) % n for i in range(n)}
    plain = [None] * n
    for s, d in perm.items():
        plain[d] = {"signature": sigs[s], "value": list(vals[s])}

    def build(counter):
        entries = [{"source_position": s,
                    "destination_position": _TaggedDestination(d, counter)}
                   for s, d in sorted(perm.items())]
        # ONE IDENTITY SLOT, not zero. § 5.3's stage-two recomputation is a
        # committed-row addressing site, and with no slots this control never
        # reached it: the bound said all such addressing is map-derived while
        # that site went unvisited, so the identity stage could be rewritten
        # to align through the signature join and P4 stayed green. A synthetic
        # that skips a stage cannot bound what the stage does.
        # The synthetic declares the version it is written against, as the
        # real packet does; a fixture that omitted it would now be rejected.
        return ({"format_version": REF.FORMAT_VERSION,
                 "rows": rows,
                 "identities": [{"position": 0,
                                 "expected_value": list(vals[1]),
                                 "factors": [{"signature": sigs[1],
                                              "exponent": 1}]}],
                 # The synthetic map carries the same HEADER the real one does.
                 # A synthetic that declared different semantics from the ones
                 # the consumer implements would be rejected, and rightly.
                 "indexing_map": dict(REF.INDEXING_MAP_HEADER, entries=entries),
                 "adjudicates": copy.deepcopy(PACKET["adjudicates"]),
                 "convention_map": copy.deepcopy(PACKET["convention_map"]) |
                 {"orientation_anchor_rule":
                  dict(PACKET["convention_map"]["orientation_anchor_rule"],
                       anchor_row_signature=sigs[0], anchor_row_label="S0")}},
                _MapOnlyRows(plain, counter))

    def counter():
        return {"map": 0, "tagged": 0, "plain": 0, "iter": 0}

    cm, cj = counter(), counter()
    pk_m, rows_m = build(cm)
    map_ok, map_why = REF.compare(pk_m, rows_m, GROUP, CONSTRUCTION, corr_mode="map")
    pk_j, rows_j = build(cj)
    # The paired negative. A consumer that validates the map and then aligns
    # through the join reaches only stand-in values, so it must now DISAGREE.
    join_ok, join_why = REF.compare(pk_j, rows_j, GROUP, CONSTRUCTION, corr_mode="join")

    # MAP CORRECTNESS, one subcase per branch of `check_indexing_map`, each
    # required to reject BY THE BRANCH IT NAMES. Requiring only "rejected" was
    # not enough: three of these branches could be hollowed and P4 stayed green.
    def malform(fn):
        """Apply one malformation and report how the consumer answered.

        An UNCAUGHT exception is NOT a rejection. `check_indexing_map` raises and
        `compare` catches it, so a caught ValueError is the branch doing its job;
        a KeyError arriving from the comparison stage means the branch is gone
        and the malformed map got that far. Hollowing a branch used to produce a
        traceback here rather than a verdict, which reddens P4 without saying
        why."""
        pk, rows_ = build(counter())
        # The whole map, not only its entries: three of the malformations below
        # are header declarations, which the consumer read and never checked.
        fn(pk["indexing_map"])
        try:
            accepted, reason = REF.compare(pk, rows_, GROUP, CONSTRUCTION)
        except Exception as e:                                # noqa: BLE001
            return False, f"UNCAUGHT {type(e).__name__}: {e}"
        return accepted, str(reason)

    def swap_equal_valued(m):
        # S2 and S3 hold the same value, so this is undetectable by comparison.
        e = m["entries"]
        e[2]["destination_position"], e[3]["destination_position"] = \
            e[3]["destination_position"], e[2]["destination_position"]

    def boolean_position(m):
        # `False == 0` and a bool indexes a list, so a position written as JSON
        # `false` was accepted as position zero and nothing noticed.
        e = m["entries"]
        for i, entry in enumerate(e):
            if int(entry["source_position"]) == 0:
                e[i] = dict(entry, source_position=False)
                return

    BRANCHES = [
        ("states a correspondence its own rows contradict", swap_equal_valued,
         "indexing map sends canonical row"),
        ("two sources to one destination",
         lambda m: m["entries"].__setitem__(
             1, dict(m["entries"][1],
                     destination_position=m["entries"][0]["destination_position"])),
         "two sources map to destination"),
        ("a destination left uncovered", lambda m: m["entries"].pop(),
         "does not cover every destination position"),
        ("a destination out of range",
         lambda m: m["entries"].__setitem__(
             0, dict(m["entries"][0],
                     destination_position=_TaggedDestination(99, counter()))),
         "is out of range"),
        # The HEADER, which the consumer relied on and never read. Each of these
        # left the comparison agreeing on all nine rows and four identities.
        ("declares itself not zero-based", lambda m: m.__setitem__("zero_based", False),
         "indexing_map.zero_based"),
        ("renames the source domain",
         lambda m: m.__setitem__("source_domain", "bogus-source"),
         "indexing_map.source_domain"),
        ("renames the destination domain",
         lambda m: m.__setitem__("destination_domain", "bogus-destination"),
         "indexing_map.destination_domain"),
        ("a source position written as JSON false", boolean_position,
         "is not an integer"),
        # The guard loops over BOTH fields and only the source arm was reached,
        # so narrowing it to source alone left P4 green while a destination
        # written as JSON `false` or `true` compared successfully.
        ("a destination position written as JSON false",
         lambda m: m["entries"].__setitem__(
             0, dict(m["entries"][0], destination_position=False)),
         "is not an integer"),
    ]
    map_checked, wrong_branch = [], []
    for label, fn, want in BRANCHES:
        accepted, reason = malform(fn)
        if accepted:
            wrong_branch.append(f"{label}: ACCEPTED")
        elif want not in reason:
            wrong_branch.append(f"{label}: rejected but by {reason[:48]!r}")
        else:
            map_checked.append(label)
    branch_failures = bool(wrong_branch)   # named for what it reports

    # Same repair: `len(map_checked) == BRANCH_COUNT` rather than "no branch
    # failed", which an empty BRANCHES satisfies while printing "0/8 rejected".
    BRANCH_COUNT = 9
    # `join_why == "DISAGREEMENT"`, not `not join_ok`: any failure of join mode
    # satisfied the old form, including one with nothing to do with the
    # stand-in. Breaking `corr_mode` dispatch made join mode raise, and the
    # control printed "DISAGREES, as ValueError: unknown corr_mode" while green.
    ok = (map_ok and join_why == "DISAGREEMENT" and cm["tagged"] > 0
          and not branch_failures
          and len(BRANCHES) == BRANCH_COUNT and len(map_checked) == BRANCH_COUNT)
    return (ok, f"4-cycle derangement, committed values served ONLY to a tagged map "
                f"destination: map mode {'agrees' if map_ok else 'FAILS: ' + map_why}, "
                f"so every value it compared came through the map "
                f"({cm['tagged']} tagged, {cm['plain']} plain, {cm['iter']} iteration"
                f"(s) for the signature join); the same comparison aligned through the "
                f"join instead {'DISAGREES' if not join_ok else 'STILL AGREES, so the '
                'mask never fired'}"
                f"{', as ' + str(join_why) if not join_ok else ''}; malformed map "
                f"{len(map_checked)}/{len(BRANCHES)} rejected BY "
                f"`check_indexing_map` itself"
                + ("" if not wrong_branch else f" -- FAILURES: {wrong_branch}"))


def control_4_downstream_cell_mutation():
    """One reference cell mutated on an in-memory copy AFTER the hash check
    completed, so the test exercises the comparison layer rather than SHA-256.

    TWO cells, at two stages. A row value alone left the § 5.3 stage-two identity
    recomputation entirely unexercised: the loop could be replaced by one over an
    empty list and every control stayed green. An identity's expected value is
    mutated too, with every row left untouched, so the rejection has to come from
    the identity stage or not at all."""
    pkt = REF.load_hash_verified(PACKET_PATH, PACKET_SHA)

    row_bad = copy.deepcopy(pkt)
    a, b, c = row_bad["rows"][3]["value"]
    row_bad["rows"][3]["value"] = [a + 1, b, c]
    row_ok, row_why = REF.compare(row_bad, COMMITTED, GROUP, CONSTRUCTION)

    idy_bad = copy.deepcopy(pkt)
    ea, eb, ec = idy_bad["identities"][1]["expected_value"]
    idy_bad["identities"][1]["expected_value"] = [ea + 1, eb, ec]
    idy_ok, idy_why = REF.compare(idy_bad, COMMITTED, GROUP, CONSTRUCTION)
    from_identity = "identity slot" in idy_why

    return (not row_ok and not idy_ok and from_identity,
            f"row cell: {'reddened' if not row_ok else 'ACCEPTED'}; identity expected "
            f"value with all rows untouched: "
            f"{'reddened' if not idy_ok else 'ACCEPTED'}"
            f"{' at the identity stage' if from_identity else ' but NOT at the identity stage'}")


def _triple_shaped_literals(src):
    """Every 3-element list/tuple whose members are all numeric EXPRESSIONS.

    Syntax-aware rather than textual. A regex over the source missed
    `(1 + 0, 0, 1)`, because it required a bare integer followed by a comma, and
    a value written as an expression is still a transcribed value."""
    def numeric(node):
        if isinstance(node, ast.Constant):
            return isinstance(node.value, (int, float)) and not isinstance(node.value, bool)
        if isinstance(node, ast.UnaryOp):
            return numeric(node.operand)
        if isinstance(node, ast.BinOp):
            return numeric(node.left) and numeric(node.right)
        return False
    out = []
    for node in ast.walk(ast.parse(src)):
        if isinstance(node, (ast.List, ast.Tuple)) and len(node.elts) == 3 \
                and all(numeric(e) for e in node.elts):
            out.append(ast.unparse(node))
    return out


def _committed_assignments(src):
    """Numeric literals in any statement that assigns to COMMITTED.

    The comparison-input construction is the path a transcribed reference value
    would actually be inserted into, and an earlier version of this control
    could not see it: it located the block by searching for the text
    `COMMITTED = `, which occurs FIRST inside this control's own source, so the
    scan sliced itself and read its own code. Located by syntax now."""
    out = []
    for node in ast.walk(ast.parse(src)):
        targets = (node.targets if isinstance(node, ast.Assign)
                   else [node.target] if isinstance(node, ast.AugAssign) else [])
        if not targets:
            continue
        names = {n.id for tgt in targets for n in ast.walk(tgt)
                 if isinstance(n, ast.Name)}
        if "COMMITTED" not in names:
            continue
        for n in ast.walk(node.value):
            if isinstance(n, ast.Constant) and isinstance(n.value, (int, float)) \
                    and not isinstance(n.value, bool):
                out.append(ast.unparse(node)[:60])
                break
    return out


def control_5_no_manual_transcription():
    """A SCREEN, and named one deliberately.

    It scans the reference fixture for triple-shaped literals and this file's
    COMMITTED assignments for numeric constants, both by syntax rather than by
    text. What it CANNOT do is establish "no manual transcription anywhere":
    that property rests on provenance, namely that every value compared is
    derived from the pinned artifacts handed in and the fixture carries no table.

    Two earlier versions overclaimed. The first scanned only the fixture, so a
    value inserted into the comparison-input path here was invisible. The second
    tried to widen that and located the input block by text search, matching its
    own source instead. Both are why this is now described as a screen.
    """
    # FIXED PATHS, not cwd-relative. Reading "p4_reference.py" relative to the
    # working directory let the SUBJECT be chosen by where the battery was
    # launched from: a transcribed triple in the fixture is caught from `build/`
    # and invisible from a directory holding a one-line decoy of the same name.
    # Same class as the delegate that defeated the isolation limb, at the two
    # sites that repair did not reach.
    here = pathlib.Path(__file__).resolve().parent
    hits = [f"fixture: {h}" for h in
            _triple_shaped_literals((here / "p4_reference.py").read_text())]
    hits += [f"input path: {h}" for h in
             _committed_assignments((here / "p4.py").read_text())]
    # A POSITIVE CONTROL. "The scan found nothing" is satisfied by a scanner
    # that cannot find anything: making both return `[]` left this green. Each
    # is now run on a planted sample first, and a scanner that comes up empty on
    # its own bait fails the control rather than clearing it.
    BAIT_FIXTURE = "_planted = [7, -3, 11]\n"
    BAIT_INPUT = "def _f():\n    COMMITTED = [{'value': [7, -3, 11]}]\n"
    blind = []
    if not _triple_shaped_literals(BAIT_FIXTURE):
        blind.append("the literal scanner finds nothing in a planted sample")
    if not _committed_assignments(BAIT_INPUT):
        blind.append("the input-path scanner finds nothing in a planted sample")
    if blind:
        return False, "SCANNER BLIND: " + "; ".join(blind)

    return (not hits,
            "screen clear on a scan proven to fire: no triple-shaped literal in "
            "the fixture, no numeric "
            "constant in the comparison-input construction; the property itself "
            "rests on provenance, not on this scan"
            if not hits else f"transcription found: {hits[:3]}")


def control_6_adjudicates_binding():
    """§ 4.4: `adjudicates` must match the run's ACTUAL packets, checked by the
    CONSUMER. The fixture took no public packets at all, so zeroing a hash after
    a successful hash-verified load left the comparison green."""
    pkt = REF.load_hash_verified(PACKET_PATH, PACKET_SHA)
    base_ok, _ = REF.compare(pkt, COMMITTED, GROUP, CONSTRUCTION)
    bad = copy.deepcopy(pkt)
    bad["adjudicates"]["group_packet_sha256"] = "0" * 64
    zero_ok, why = REF.compare(bad, COMMITTED, GROUP, CONSTRUCTION)
    # BOTH halves, separately. Mutating only the group hash left the
    # construction branch removable with this control still green.
    badc = copy.deepcopy(pkt)
    badc["adjudicates"]["construction_packet_sha256"] = "0" * 64
    c_ok, c_why = REF.compare(badc, COMMITTED, GROUP, CONSTRUCTION)
    # and the rejection must come from the adjudicates binding, not from the
    # basing reference, which pins the same packet by a different route.
    from_adj = "`adjudicates` does not match" in c_why
    sub_g, _ = REF.compare(pkt, COMMITTED, GROUP + b" ", CONSTRUCTION)
    sub_c, _ = REF.compare(pkt, COMMITTED, GROUP, CONSTRUCTION + b" ")
    ok = base_ok and not zero_ok and not c_ok and from_adj and not sub_g and not sub_c
    return (ok, f"real packet {'binds' if base_ok else 'FAILS'}; zeroed group hash "
                f"{'rejected' if not zero_ok else 'ACCEPTED'}; zeroed construction hash "
                f"{'rejected' if not c_ok else 'ACCEPTED'}"
                f"{' by the adjudicates binding itself' if from_adj else ' but NOT by the adjudicates binding'}; "
                f"substituted group packet {'rejected' if not sub_g else 'ACCEPTED'}; "
                f"substituted construction packet {'rejected' if not sub_c else 'ACCEPTED'}")


# The convention-map requirements the consumer MUST validate, frozen HERE on the
# control side and never read from the fixture.
#
# The first version of this control iterated `REF.IMPLEMENTED_SEMANTICS`, so
# deleting a requirement from the fixture simply made the control test one fewer
# thing and P4 stayed green: the control derived its expectation from the object
# under test. That is the same defect as the P2 convention-leaf generator, at a
# new site, which is why it is written out here instead.
# The convention-map declarations the consumer MUST rely on and validate, frozen
# HERE on the control side as KEYS and never read from the fixture.
#
# An earlier version iterated the fixture's own list, so deleting a requirement
# made the control test one fewer thing and P4 stayed green: the control derived
# its expectation from the object under test. That was the fourth occurrence of
# that class in this build.
EXPECTED_DECLARATIONS = (
    ("bridge", "identification"),
    ("bridge", "exclusivity"),
    ("orientation_anchor_rule", "selection_mapping"),
    ("orientation_anchor_rule", "uniqueness_gate"),
)


def control_7_convention_map_validated():
    """The consumer LOADS AND VALIDATES its convention map.

    NOT one of the six § 4.4 controls, despite what this line used to say. § 5.3
    and § 5.4 are what require a consumer to know the semantics it is applying;
    § 4.4's own list stops short of this. The category split is a scope claim the
    report makes, so mislabelling it here would have contradicted the report from
    inside the code it describes.

    Three things.

    INVENTORY. The fixture's relied-on set must equal the frozen set above, so a
    declaration dropped from the fixture is a failure rather than one less test.

    CONTRADICTION. Each declaration is rewritten to assert the OPPOSITE while
    retaining its distinctive wording. This is the case that broke the previous
    version: it validated by substring presence, so "DO NOT apply the inverse to
    EVERY ROW AT ONCE; apply it only to R7" contained the required phrase and was
    accepted. Substring presence is not validation, and the fixture compares
    complete declared values now.

    REMOVAL. Each declaration is also deleted outright, which a substring check
    would have caught but which still has to hold.
    """
    # format_version rides alongside the convention-map leaves: same class,
    # same control. Two subcases, wrong and removed.
    fv_bad = []
    for label, mutate_fv in (("format_version wrong",
                              lambda p: p.__setitem__("format_version", "m8_9-INCOMPATIBLE")),
                             ("format_version removed",
                              lambda p: p.pop("format_version", None))):
        bad = copy.deepcopy(PACKET)
        mutate_fv(bad)
        acc, _ = REF.compare(bad, COMMITTED, GROUP, CONSTRUCTION)
        if acc:
            fv_bad.append(label)

    fixture = set(REF.RELIED_ON_DECLARATIONS)
    frozen = set(EXPECTED_DECLARATIONS)
    if fixture != frozen:
        return False, (f"the consumer's relied-on set does not match the frozen "
                       f"inventory: missing {sorted(frozen - fixture)}, "
                       f"unexpected {sorted(fixture - frozen)}")

    pkt = REF.load_hash_verified(PACKET_PATH, PACKET_SHA)
    accepted = []
    for section, leaf in EXPECTED_DECLARATIONS:
        real = pkt["convention_map"][section][leaf]
        contradiction = f"DO NOT rely on the following, it does not hold: {real}"
        for label, value in (("contradicted", contradiction), ("removed", "")):
            bad = copy.deepcopy(pkt)
            bad["convention_map"][section][leaf] = value
            if REF.compare(bad, COMMITTED, GROUP, CONSTRUCTION)[0]:
                accepted.append(f"{section}.{leaf} {label}")
    basing = copy.deepcopy(pkt)
    basing["convention_map"]["basing_reference"]["construction_packet_sha256"] = "0" * 64
    if REF.compare(basing, COMMITTED, GROUP, CONSTRUCTION)[0]:
        accepted.append("basing_reference")
    accepted.extend(fv_bad)
    total = len(EXPECTED_DECLARATIONS) * 2 + 1 + 2
    return (not accepted,
            f"inventory matches the frozen set; {total} mutations, each declaration "
            "both CONTRADICTED while keeping its wording and removed outright, plus "
            "`format_version` wrong and removed: all rejected" if not accepted
            else f"still accepted: {accepted}")


def control_8_orientation_selection_exercised():
    """§ 5.4's selection, exercised on all three of its outcomes.

    Every other fixture selects the COMMITTED branch, so the whole selector could
    be replaced by a constant, and the global inversion disabled, with every
    control still green. The report calls this fixture the minimum consumer of
    § 5.3 and § 5.4; a consumer that never inverts has not shown it implements
    § 5.4.

    The inverse arm uses the real packet against globally inverted committed
    values, which is exactly the case § 5.4's one-bit selection exists for. The
    invalid-anchor and disagreement arms need a self-inverse and a mismatched
    anchor, neither of which the real packet can provide, so they are small
    synthetics."""
    pkt = REF.load_hash_verified(PACKET_PATH, PACKET_SHA)
    inverted = [{"signature": r["signature"], "value": list(Q.inv(tuple(r["value"])))}
                for r in COMMITTED]
    inv_ok, inv_why = REF.compare(pkt, inverted, GROUP, CONSTRUCTION)
    selected_inverse = "orientation inverse" in inv_why

    def synth(anchor_committed, anchor_reference):
        sigs = [{"dimension": i + 1, "s": [i, 0, 1], "t": [0, i, 1], "st": [1, i, 1]}
                for i in range(2)]
        rows = [{"label": f"S{i}",
                 "evidentiary_class": p3.CLASS_SELECTOR if i == 0 else p3.CLASS_FREE,
                 "value": list(v), "signature": sigs[i]}
                for i, v in enumerate((anchor_reference, (3, 0, 1)))]
        com = [{"signature": sigs[i], "value": list(v)}
               for i, v in enumerate((anchor_committed, (3, 0, 1)))]
        # The synthetic declares the version it is written against, as the
        # real packet does; a fixture that omitted it would now be rejected.
        return ({"format_version": REF.FORMAT_VERSION,
                 "rows": rows, "identities": [],
                 "indexing_map": dict(REF.INDEXING_MAP_HEADER,
                                      entries=[{"source_position": i,
                                                "destination_position": i}
                                               for i in range(2)]),
                 "adjudicates": copy.deepcopy(pkt["adjudicates"]),
                 "convention_map": copy.deepcopy(pkt["convention_map"]) |
                 {"orientation_anchor_rule":
                  dict(pkt["convention_map"]["orientation_anchor_rule"],
                       anchor_row_signature=sigs[0], anchor_row_label="S0")}}, com)

    p_i, c_i = synth((1, 0, 1), (1, 0, 1))        # self-inverse reference
    _, invalid_why = REF.compare(p_i, c_i, GROUP, CONSTRUCTION)
    p_d, c_d = synth((5, 0, 1), (7, 0, 1))        # neither branch holds
    _, disagree_why = REF.compare(p_d, c_d, GROUP, CONSTRUCTION)

    # The anchor's DESIGNATION, which the fixture consumed and never checked.
    # Substituting any other row's signature was accepted for six of the eight
    # alternatives, so § 5.4's selection ran on the wrong row. Every alternative
    # must be rejected now, and the two that used to fail must fail for the
    # designation rather than for happening to be self-inverse.
    # ALTERNATES is frozen here, and `tried` is counted, because requiring only
    # that nothing failed is satisfied perfectly by trying nothing: emptying this
    # loop left the control green and still printing "all 8 alternate signatures
    # rejected". The 8 was a literal in the format string and the loop's own
    # result was never checked against it.
    ALTERNATES = 8
    swapped, tried = [], 0
    real_sig = pkt["convention_map"]["orientation_anchor_rule"]["anchor_row_signature"]
    for row in pkt["rows"]:
        if row["signature"] == real_sig:
            continue
        bad = copy.deepcopy(pkt)
        bad["convention_map"]["orientation_anchor_rule"]["anchor_row_signature"] = \
            row["signature"]
        acc, why = REF.compare(bad, COMMITTED, GROUP, CONSTRUCTION)
        tried += 1
        if acc or "designates" not in str(why):
            swapped.append(f"{row['label']}:{'ACCEPTED' if acc else why}")
    bad_label = copy.deepcopy(pkt)
    bad_label["convention_map"]["orientation_anchor_rule"]["anchor_row_label"] = "R1"
    label_acc, _ = REF.compare(bad_label, COMMITTED, GROUP, CONSTRUCTION)

    ok = (inv_ok and selected_inverse and invalid_why == "INVALID_ANCHOR"
          and disagree_why == "DISAGREEMENT" and tried == ALTERNATES
          and not swapped and not label_acc)
    return (ok, f"globally inverted committed table: "
                f"{'agrees under the INVERSE orientation' if inv_ok and selected_inverse else 'FAILED: ' + inv_why}; "
                f"self-inverse anchor -> {invalid_why}; "
                f"neither branch -> {disagree_why}; anchor DESIGNATION: "
                f"{tried}/{ALTERNATES} alternate signatures TRIED, each rejected as "
                f"not the row § 5.4 designates"
                + ("" if not swapped else f" -- EXCEPT {swapped}")
                + f", and a mismatched anchor label "
                f"{'rejected' if not label_acc else 'ACCEPTED'}")


CONTROLS = [
    ("SAME_PARSER_PATH", control_1_same_parser_path),
    ("WRONG_HASH_BEFORE_PARSE", control_2_wrong_hash_before_parse),
    ("SYNTHETIC_NONIDENTITY_MAP", control_3_synthetic_nonidentity_map),
    ("DOWNSTREAM_CELL_MUTATION", control_4_downstream_cell_mutation),
    ("NO_MANUAL_TRANSCRIPTION", control_5_no_manual_transcription),
    ("ADJUDICATES_BINDING", control_6_adjudicates_binding),
    ("CONVENTION_MAP_VALIDATED", control_7_convention_map_validated),
    ("ORIENTATION_SELECTION_EXERCISED", control_8_orientation_selection_exercised),
]


# -------------------------------------------------------------------- main

def main():
    global PACKET, SIDECAR, P1, GROUP, CONSTRUCTION, PACKET_PATH, COMMITTED
    if len(sys.argv) != 6:
        print("usage: p4.py <p1> <packet> <sidecar> <group> <construction>",
              file=sys.stderr)
        return 2
    p1_path, PACKET_PATH, sc_path, g_path, c_path = sys.argv[1:6]
    P1 = json.loads(pathlib.Path(p1_path).read_bytes())
    PACKET = json.loads(pathlib.Path(PACKET_PATH).read_bytes())
    SIDECAR = json.loads(pathlib.Path(sc_path).read_bytes())
    GROUP = pathlib.Path(g_path).read_bytes()
    CONSTRUCTION = pathlib.Path(c_path).read_bytes()
    # The committed raw output as § 5.5 describes it: each row carries its own
    # SIGNATURE and its value. The first version passed bare values, which is what
    # let the fixture align positionally and never read a committed signature.
    COMMITTED = [{"signature": P1["signatures"][v["label"]],
                  "value": list(v["triple"])} for v in P1["values"]]

    print("P4 mutation battery, two halves\n")
    base = p3.gate(p3.canonical_bytes(PACKET), p3.canonical_bytes(SIDECAR),
                   P1, GROUP, CONSTRUCTION)
    if base[1]:
        print(f"  BASELINE P3 IS NOT GREEN: {base[1][:2]}")
        return 1
    print(f"  baseline: P3 green on {len(base[0])} checks\n")

    print("  HALF ONE, CONTENT: every P3 check reddened by a targeted mutation\n")
    observed, misses = set(), []
    for name, target, fn, reason, req_green, note, pre in CASES:
        pb, sb, p1, gb, cb = fn()
        pre_ok, pre_why = pre(pb, sb) if pre else (True, "")
        # `stops` as well as `checks`. `add()`'s `if not ok: stops.append(...)`
        # is the ONLY link between a red check and P3 actually rejecting, and no
        # case touched it: setting that line to `if False:` left every per-check
        # boolean correct, `stops` permanently empty, P3 printing "GREEN on all
        # 12 declared checks" on a packet with a red check, and P4 at 28/28.
        # Reading only the booleans tested the verdicts and not the rejection.
        checks, stops = p3.gate(pb, sb, p1, gb, cb)
        red = {c for c, ok, _ in checks if not ok}
        detail = {c: d for c, ok, d in checks if not ok}
        blocked = str(detail.get(target, "")).startswith("BLOCKED:")
        why_ok = reason is None or reason in str(detail.get(target, ""))
        green_ok = not any(g in red for g in req_green)
        stopped = any(s.startswith(target + ":") for s in stops)
        hit = (pre_ok and target in red and stopped and not blocked
               and why_ok and green_ok)
        observed.add(target) if hit else misses.append(name)
        print(f"    {'RED ' if hit else 'MISS'}  {name}")
        print(f"            {note}" + (f"  [{pre_why}]" if pre and pre_ok else ""))
        if not hit:
            if not pre_ok:
                print(f"            PRECONDITION FAILED: {pre_why}")
            print(f"            target {target}: "
                  f"{'not red' if target not in red else 'BLOCKED' if blocked else 'wrong reason' if not why_ok else 'tripped ' + str([g for g in req_green if g in red])}")

    never = [c for c in p3.CHECK_IDS if c not in observed]
    print(f"\n    {len(CASES) - len(misses)}/{len(CASES)} cases reddened their target")
    print(f"    coverage: {len(observed)}/{len(p3.CHECK_IDS)} P3 checks observed rejecting")
    if never:
        print(f"    NEVER OBSERVED REJECTING: {never}")

    print("\n  HALF TWO, INGESTION: six § 4.4 ingestion controls plus two § 5.3 "
          "and § 5.4\n  consumer controls, against the reference fixture\n")
    ctl_bad = []
    for cid, fn in CONTROLS:
        ok, why = fn()
        if not ok:
            ctl_bad.append(cid)
        # The bound comes from BOUNDS, never from the control's own text, so a
        # control cannot state one the scope block does not.
        if cid in BOUNDS:
            why = f"{why} BOUNDED: {BOUNDS[cid]}"
        print(f"    {'ok  ' if ok else 'FAIL'} {cid:26s} {why}")

    ok = not misses and not never and not ctl_bad
    print(f"\n  {'P4 GREEN. ' if ok else 'P4 INCOMPLETE. '}"
          f"Content: {len(CASES) - len(misses)}/{len(CASES)}, "
          f"coverage {len(observed)}/{len(p3.CHECK_IDS)}. "
          f"Ingestion: {len(CONTROLS) - len(ctl_bad)}/{len(CONTROLS)}.")
    print("\n  Scope: SIX of the eight are the § 4.4 ingestion controls ASSIGNED TO\n"
          "  P4, shown IMPLEMENTABLE against the frozen packet, so it is consumable\n"
          "  as § 4.4 requires; § 4.4's publication control is issuance and belongs\n"
          "  to P6. The other two cover § 5.3 and § 5.4 consumer behaviour and are\n"
          "  NOT § 4.4 controls. Two are explicitly bounded, and each bound is held\n"
          "  once in this file and printed from there, here and with its own result,\n"
          "  so neither can be widened in prose while the implementation stays\n"
          "  narrow:")
    for cid, bound in BOUNDS.items():
        # break_on_hyphens=False: wrapping "map-derived" across a line break
        # would leave the printed bound no longer a verbatim match for the one
        # in the report, which is the exact drift this is meant to prevent.
        print("\n".join(textwrap.wrap(f"{cid} {bound}.", 72,
                                      initial_indent="    ",
                                      subsequent_indent="    ",
                                      break_on_hyphens=False,
                                      break_long_words=False)))
    print("  The real adjudication harness is maintainer-owned, and the maintainer\n"
          "  must independently bind and mutation-test the corresponding controls\n"
          "  there. Nothing here claims that harness is correct. This battery wrote\n"
          "  nothing.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
