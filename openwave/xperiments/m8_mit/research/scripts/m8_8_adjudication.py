"""M8.8 ADJUDICATION HARNESS: protocol section 8 steps 6 to 8, written and committed PRE-REVEAL.

WHAT THIS IS.  The maintainer-side comparison harness that
[`m8_8_reproduction_protocol.md`](../findings/m8_8_reproduction_protocol.md) section 8 step 7
names and section 4.4 controls.  It is committed BEFORE the canonical answer packet is
decrypted, against the section 4.4 schema and the committed Phase A raw output only, so
the harness cannot have been shaped by the reference values.  The ordering is auditable
from the commit graph: this file lands, then the packet opens.

THE INPUTS, ALL PINNED.
  raw output        m8_8_cleanroom/attempt5/RAW_OUTPUT.json      1a9b56ce... (Addendum 1)
  manifest          m8_8_cleanroom/attempt5/METHOD_AND_GATE_MANIFEST.md   8aa140e3...
  group packet      data/m8_5a_packet.json                        e3b0c945... (section 11)
  construction      data/m8_8_construction_packet.json            df00c022... (section 11)
  answer packet     delivered out of band, plaintext              744c7f25... (section 11)
  phase B record    m8_8_cleanroom/phase_b/MUTATION_RESULTS.json  (lands by PR before step 6)

THE SEQUENCE IMPLEMENTED, section 8 steps 6 to 8, in this order and no other:
  step 6   the delivered answer-packet bytes are SHA-256 hashed against the section 11 pin
           BEFORE anything is consumed.  The pin is over the section 10.2 CANONICAL form,
           so a delivery in another rendering is parsed only to re-serialize canonically,
           and is consumed only if the canonical hash equals the pin, with the record
           stating that canonicalization changed the bytes; any other outcome refuses.
  step 7   `adjudicates` is checked against the raw output's consumed packet hashes and
           both against the section 11 pins; the raw output's own bytes are checked against
           the Addendum 1 pin; the POSITIONAL indexing map is APPLIED to pair packet rows
           with committed rows; the convention map is LOADED AND VALIDATED (four object
           items, the `representation_evaluation` item carrying the construction packet's
           `basing.evaluation` byte-equal) without applying any orientation.
  step 8   the section 5.4 selection at R7: exactly one of x = r, x^-1 = r must hold
           (both: INVALID ANCHOR; neither: DISAGREEMENT; no row comparison in either case).
           The selected identity or global inverse is applied to the COMMITTED rows, the
           seven free forms and R7 compared exactly in Q(phi), and the four identities
           recomputed from the SELECTED rows under the packet's slot definitions and
           compared POSITION-WISE (section 5.3).  Outcome recorded by section 8 category.

THE PACKET SHAPE THIS HARNESS READS.  Section 4.4 fixes the key set and describes each
value; it does not print a JSON skeleton.  The first committed reading of this harness
(PR #452) guessed a signature-dispatched map; the author then extracted the packet's
key-level shape mechanically (#451, 2026-08-21, keys and types only, no values) and it
differs in five places, the substantive one being a POSITIONAL indexing map.  The first
adapter (PR #453) bound that shape with three readings left PENDING and failing closed;
the author answered all three on #451 (2026-08-21 14:28), and this second adapter pins
them.  Both are visible pre-reveal diffs landed before any ciphertext exists.
  A THIRD adapter (2026-08-22, POST-reveal, visible) changed two things and no semantics:
  the source-domain STRING is the author's builder's own (`packet_rows_canonical_position`,
  read off the delivered packet; the #453 placeholder `answer_packet.rows` is refused as
  superseded), and the terminal summary no longer crashes on the `_sector_coverage` record
  entry (the record was already written; only the exit code was lost).  C3k and C17 cover
  both, and the committed harness's step-6 refusal record is published with the rerun.

  {
    "format_version": str,
    "target_id": str,
    "adjudicates": {"group_packet_sha256": hex, "construction_packet_sha256": hex},
    "rows": [ {"label": str,
               "signature": {"dimension": int, "s": [a,b,c], "t": [a,b,c], "st": [a,b,c]},
               "evidentiary_class": "declared_convention" | "free" | "free_orientation_selector",
               "value": [a,b,c]}                x 9 ],
    "identities": [ {"position": int,
                     "factors": [{"signature": {...}, "exponent": int}, ...]  (EXACT keys),
                     "expected_value": [a,b,c]}  x 4 ],
    "indexing_map": {"source_domain": SOURCE_DOMAIN,
                     "destination_domain": "m8_3_record_row_position", "zero_based": bool,
                     "entries": [{"source_position": int, "destination_position": int}] x 9},
    "convention_map": {
      "bridge": {exclusivity, identification, protocol_section, statement},
      "orientation_anchor_rule": {anchor_row_label, anchor_row_signature{dimension,s,t,st},
                                  protocol_section, rule, selection_mapping, uniqueness_gate},
      "basing_reference": {construction_packet_sha256, protocol_section, source},
      "representation_evaluation": {convention, protocol_section, source} }
  }

  Every Q(phi) value is the normalized triple (a, b, c) for (a + b*phi)/c, c > 0,
  gcd(a, b, c) = 1, as section 4.2 fixes and the raw output already uses.  The committed
  RAW_OUTPUT.json keeps its own field names (`row_signature` with `dim`, `chi_s`, `chi_t`,
  `chi_st`); both spell the same section 5.5 signature and reduce to one key.

  THE POSITIONAL MAP (author, #451).  `source_position` indexes the packet's CANONICAL
  row order, its rows sorted by the ten-integer signature key, never the file order;
  `destination_position` indexes the M8.3 record's own row order R0..R8
  (`m8_3_record_row_position`); `zero_based` fixes the base of both.  The destination
  order is DERIVED FROM THE TWO PUBLIC PACKETS ALONE (`m8_3_record_order`: the labels are
  fixed publicly by dimension and McKay distance up to the Galois flip, each Galois pair
  is split through the group packet's embedding as section 5.5 states, R1 being the
  embedding itself, and the signatures follow by the symmetric-power tower), never from the
  realized RAW_OUTPUT.json array position, which the packet's 2026-08-10 freeze predates
  and which the author ruled out as a reading.  The derived nine keys must be exactly the
  committed signatures or the run refuses.  Source positions must cover the nine packet
  rows exactly once, destination positions must be distinct and in range (injectivity),
  and the pairing is DRIVEN BY THE MAP: the packet's own signatures are checked to be
  distinct and present among the committed rows, and a LIVE row whose own signature
  differs from the committed row the map pairs it with is REFUSED (section 5.5 makes the
  signature the identity), while the self-test arms record the difference so ignoring the
  map fails AT THE COMPARISON.  On the committed table the canonical order and the record
  order differ (R5 and R6 exchange), so the honest map IS section 4.4's nonidentity
  fixture: applying it versus ignoring it changes the outcome at the comparison (C3).

  RESOLVED BY THE AUTHOR (#451, 2026-08-21 14:28), each pinned at its point of use:
  (1) the domains, above; (2) the inner keys of the four convention objects, exact
  (CONVENTION_ITEM_KEYS), with `representation_evaluation.convention` the field byte-
  compared to the construction packet's `basing.evaluation` and every other field
  recorded, never gated; (3) the remaining keys: rows `evidentiary_class` and `value`,
  identities `expected_value`, factors EXACTLY `signature` and `exponent`: `conjugate` is
  outside the frozen factor schema and, since present-and-true would change the computed
  identity, it is REFUSED at the schema layer rather than admitted as optional (the #453
  reading), and the evaluator carries no conjugate branch.

THE CONTROLS, section 4.4, run by --self-test and REQUIRED green before any live run.
An independent adversarial audit of the first version found C1, C2 and C3 tautological
(a comparison skipping six of eight rows, or a set-wise identity comparison, passed them)
and found the anchor label-driven rather than R7-driven; the controls below are the
repaired set, each phrased so the named defect would fail it:
  C1  EACH of the eight nontrivial reference cells mutated in turn, downstream of the
      completed hash check: exactly that one label mismatches (R7: `disagreement`);
  C2  EACH committed raw cell mutated the same way, from the other side;
  C3  the NONIDENTITY positional fixture is the real relation: the honest map (canonical
      packet order -> M8.3 record order) is a non-identity permutation on the committed
      table, moving exactly the rows where the two orders differ; the supplied map is
      GREEN; ignoring the map (identity arm) and a preregistered wrong permutation (three
      free rows cycled) each fail AT THE COMPARISON on exactly the moved rows, not by an
      earlier refusal; a map moving R0 onto an acyclic row refuses on class (C3d); a LIVE
      packet whose own signatures contradict the map is refused naming both signatures,
      and the same fixture as a self-test arm is green with the difference RECORDED
      (C3e); a one-based map is honored (C3f), a misdeclared base refuses on range (C3g),
      a source position listed twice (C3h), a destination past the committed rows (C3i)
      and a JSON boolean position (C3j) refuse; an unimplemented domain string, including
      the ruled-out `RAW_OUTPUT.rows`, is refused before any pairing, the map FAILING
      CLOSED (C3k); a second domain string maps through ITS OWN order (C3l); registered
      orders are pairwise distinct and none is the committed file order (C3m); the
      public derivation of R0..R8 reproduces the committed signature set and refuses the
      binding when a committed signature is altered (C3n); the Galois convention has an
      oracle, R1 = the embedding's trace recomputed by a second route, R2/R4 its
      conjugates, R3 on R1's side (C3p: the (dim, distance) table alone is blind to the
      flip, an audit showed); exchanged label strings are recorded, never gated (C3q);
      the packet's file order is irrelevant to the pairing (C3o);
  C4  the global inverse of the committed table resolves to `convention difference`; the
      two sector-product expected values swapped redden BOTH slots on the identity and the
      inverted fixture, which is the position-wise rule of section 5.3 exercised on the
      harness rather than on the fixture arithmetic;
  C5  a self-inverse R7 reference resolves to `invalid anchor` with no row comparison;
  C6  an R7 reference equal to neither x nor x^-1 resolves to `disagreement`;
  C7  a tampered identity `expected` reddens the identity layer alone;
  C8  tampered `adjudicates`, content-tampered bytes, a duplicate key in the delivered
      bytes, a `representation_evaluation.convention` not byte-equal to the declaration,
      a convention object with a wrong inner key set (one per item), and a convention
      item delivered as a bare string each REFUSE; the verbatim field is recorded by
      name; an uncanonical rendering of the pinned object is accepted with the change
      recorded, and pin-matching bytes in another canonicalizer's rendering are accepted
      with the form recorded; the bridge, anchor-rule and source prose are RECORDED,
      never gated (a second audit showed a substring check refusing section 5.4's own
      wording);
  C9  the selector class on any row but the dim-5 irrep refuses, even with a correct table;
  C10 a non-injective map refuses (run as a self-test arm: on the live path the C3e
      self-contradiction check is reached first, so a live non-injective map refuses there);
  C11 structurally trivial identities refuse; two identities at one position, or a
      position that is not an int, refuse;
  C12 a malformed but pinned packet yields a RECORDED structural failure, not a traceback;
  C13 EXACTNESS: negation, Galois conjugation and a 1e-20 near-miss of one cell are each
      exactly one mismatch, so an approximate, sign-blind or conjugation-blind equality
      would fail here;
  C14 one control per identity rule (non-Galois ratio, overlapping products, R0 factor,
      zero exponent, repeated factor, same pair in both ratios, a `conjugate` key on a
      factor as true and as a string, exponent beyond the bound), plus a conjugate factor
      refused at the schema layer even when its expected value was recomputed with it;
  C15 the Phase B checker, each rule in isolation on a synthetic record;
  C16 empty `target_id`, an own signature absent from the committed output, a
      declared-convention value other than 1, an extra key on a row, and the superseded
      #452 signature spelling each refuse (rows, identities and factors carry exact key
      sets).
  The fixtures are built from the committed raw table and public data only; they contain
  no reference value and prove nothing about the reproduction, only about the harness.
  `hypothesis failure` is unreachable with this raw output (R0 is its only non-acyclic row
  and must pair with `declared_convention`), so no control exercises it.

WHAT THIS FILE DOES NOT DO.  It does not re-derive torsion, does not read the sealed
packet until step 6, does not assign evidentiary classes (those are answer-side metadata,
read from the packet after signature matching, section 5.5), and does not upgrade any
claim: section 9's explicitly-not-verified list stands whatever the category.

EXIT CODES.  0 = adjudication completed, `reproduced` or `convention difference` (both
successes under section 8).  2 = adjudication completed with a finding: `partial
disagreement`, `disagreement`, `invalid anchor`, `hypothesis failure`.  1 = refusal or
`structural failure`, or a self-test control that did not behave.

Run:
  python3 m8_8_adjudication.py --self-test
  python3 m8_8_adjudication.py --packet PATH [--phase-b PATH] [--json PATH]
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import sys
from fractions import Fraction
from math import gcd
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parent
RAW_OUTPUT_PATH = RESEARCH / "m8_8_cleanroom" / "attempt5" / "RAW_OUTPUT.json"
MANIFEST_PATH = RESEARCH / "m8_8_cleanroom" / "attempt5" / "METHOD_AND_GATE_MANIFEST.md"
GROUP_PACKET_PATH = RESEARCH / "data" / "m8_5a_packet.json"
CONSTRUCTION_PACKET_PATH = RESEARCH / "data" / "m8_8_construction_packet.json"
PHASE_B_DEFAULT = RESEARCH / "m8_8_cleanroom" / "phase_b" / "MUTATION_RESULTS.json"
DEFAULT_OUT = RESEARCH / "data" / "m8_8_adjudication.json"

# Section 11 and Addendum 1 pins, every one recomputed at the #451 review.
PIN_ANSWER_PACKET = "744c7f25e2312d90fc356b11510da685328f05e80ae62721d0a0f418dcf9697e"
PIN_GROUP_PACKET = "e3b0c945bbbb15b4549fa641234c9461062c2337b3d1e372af621b614d4883a9"
PIN_CONSTRUCTION_PACKET = "df00c0222f98c481eb56b882cd867a6c3a4f8604b8633e81dec0cce1f8460a06"
PIN_RAW_OUTPUT = "1a9b56ce70bae73e5cf8c4ef00f6e43bf76937afb9075801605f6bf5047d1002"
PIN_MANIFEST = "8aa140e3978366ca38f7c1d5926d1a2972305733be434595f2905e9df512f838"

RAW_SCHEMA = "m8_8-raw-output-1"
PACKET_KEYS = {
    "format_version",
    "target_id",
    "adjudicates",
    "rows",
    "identities",
    "indexing_map",
    "convention_map",
}
CONVENTION_KEYS = {
    "bridge",
    "orientation_anchor_rule",
    "basing_reference",
    "representation_evaluation",
}
INDEXING_KEYS = {"source_domain", "destination_domain", "zero_based", "entries"}
ENTRY_KEYS = {"source_position", "destination_position"}
ROW_KEYS = {"label", "signature", "evidentiary_class", "value"}
IDENTITY_KEYS = {"position", "factors", "expected_value"}
FACTOR_KEYS = {"signature", "exponent"}  # EXACT: `conjugate` is outside the frozen schema
# The four convention objects' inner key sets (author, #451 2026-08-21 14:28), exact.  Only
# `representation_evaluation.convention` is gated (byte-equal to the construction packet's
# `basing.evaluation`); every other field is recorded for the adjudication record.
CONVENTION_ITEM_KEYS = {
    "bridge": {"exclusivity", "identification", "protocol_section", "statement"},
    "orientation_anchor_rule": {
        "anchor_row_label",
        "anchor_row_signature",
        "protocol_section",
        "rule",
        "selection_mapping",
        "uniqueness_gate",
    },
    "basing_reference": {"construction_packet_sha256", "protocol_section", "source"},
    "representation_evaluation": {"convention", "protocol_section", "source"},
}
EVALUATION_FIELD = "convention"
# The domains (author, #451 2026-08-21 14:28).  SOURCE: the packet's canonical row order,
# its rows sorted by the ten-integer signature key (NOT file order).  DESTINATION:
# `m8_3_record_row_position`, the M8.3 record's own row order R0..R8, derived below from the
# two PUBLIC packets alone (`m8_3_record_order`), never from the committed array position.
# The map FAILS CLOSED: any other domain string is refused before pairing.
# The SOURCE STRING is the author's builder's own, read off the delivered packet at step 6
# (2026-08-22): `packet_rows_canonical_position`.  The #453 spelling `answer_packet.rows`
# was a maintainer placeholder (the author's 2026-08-21 shape post gave keys and types
# only, no values) and is refused as superseded; the semantics did not move.
SOURCE_DOMAIN = "packet_rows_canonical_position"
SOURCE_DOMAINS = {SOURCE_DOMAIN}
DESTINATION_DOMAIN = "m8_3_record_row_position"
# destination domain string -> the function producing that ordering of committed signature
# keys; the order is looked up BY the string, so a vocabulary entry without its ordering
# cannot silently reuse another domain's order (filled next to `m8_3_record_order` below).
DESTINATION_ORDERS: dict = {}
CLASSES = {"declared_convention", "free", "free_orientation_selector"}
SIG_FIELDS = ("dimension", "s", "t", "st")  # the answer packet's `signature` (author-confirmed)
RAW_SIG_FIELDS = ("dim", "chi_s", "chi_t", "chi_st")  # the committed raw output's `row_signature`
N_ROWS = 9
N_IDENTITIES = 4
R7_DIM = 5  # the orientation anchor is the dim-5 irrep (M8.5-A label map); public
MAX_EXPONENT = 8
PHASE_B_SCHEMA_PREFIX = "m8_8-phase-b-mutation-results-"

SUCCESS = {"reproduced", "convention difference"}
FINDINGS = {"partial disagreement", "disagreement", "invalid anchor", "hypothesis failure"}


class Refusal(Exception):
    """A precondition failed before any comparison; recorded as `structural failure`."""


# ----------------------------------------------------------------------------------------
# Exact Q(phi) arithmetic.  x + y*phi with x, y in Q and phi^2 = phi + 1.
# ----------------------------------------------------------------------------------------
class QPhi:
    __slots__ = ("x", "y")

    def __init__(self, x, y=0):
        self.x = Fraction(x)
        self.y = Fraction(y)

    @classmethod
    def from_triple(cls, t) -> "QPhi":
        if (
            not isinstance(t, list)
            or len(t) != 3
            or not all(isinstance(v, int) and not isinstance(v, bool) for v in t)
        ):
            raise Refusal(f"Q(phi) value is not an integer triple: {t!r}")
        a, b, c = t
        if c <= 0:
            raise Refusal(f"Q(phi) triple with c <= 0: {t!r}")
        if gcd(gcd(abs(a), abs(b)), c) != 1:
            raise Refusal(f"Q(phi) triple not normalized (gcd != 1): {t!r}")
        return cls(Fraction(a, c), Fraction(b, c))

    def to_triple(self) -> list:
        c = self.x.denominator * self.y.denominator // gcd(self.x.denominator, self.y.denominator)
        a = self.x.numerator * (c // self.x.denominator)
        b = self.y.numerator * (c // self.y.denominator)
        g = gcd(gcd(abs(a), abs(b)), c)
        return [a // g, b // g, c // g]

    def __eq__(self, other) -> bool:
        return isinstance(other, QPhi) and self.x == other.x and self.y == other.y

    def __hash__(self):
        return hash((self.x, self.y))

    def __add__(self, other: "QPhi") -> "QPhi":
        return QPhi(self.x + other.x, self.y + other.y)

    def __sub__(self, other: "QPhi") -> "QPhi":
        return QPhi(self.x - other.x, self.y - other.y)

    def __mul__(self, other: "QPhi") -> "QPhi":
        # (x1 + y1 p)(x2 + y2 p) = x1x2 + (x1y2 + y1x2) p + y1y2 (p + 1)
        return QPhi(
            self.x * other.x + self.y * other.y,
            self.x * other.y + self.y * other.x + self.y * other.y,
        )

    def conjugate(self) -> "QPhi":
        # phi -> 1 - phi
        return QPhi(self.x + self.y, -self.y)

    def norm(self) -> Fraction:
        return self.x * self.x + self.x * self.y - self.y * self.y

    def inverse(self) -> "QPhi":
        n = self.norm()
        if n == 0:
            raise Refusal("inverse of zero in Q(phi)")
        c = self.conjugate()
        return QPhi(c.x / n, c.y / n)

    def __pow__(self, e: int) -> "QPhi":
        if abs(e) > MAX_EXPONENT:
            raise Refusal(f"exponent {e} exceeds the bound {MAX_EXPONENT}")
        base = self if e >= 0 else self.inverse()
        out = QPhi(1, 0)
        for _ in range(abs(e)):
            out = out * base
        return out

    def __repr__(self):
        return f"QPhi{tuple(self.to_triple())}"


ONE = QPhi(1, 0)


# ----------------------------------------------------------------------------------------
# Loading, with the section 4.4 / 10.2 preconditions.
# ----------------------------------------------------------------------------------------
def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def canonical_bytes(obj) -> bytes:
    """Section 10.2: keys sorted, two-space indent, ASCII, LF, single trailing newline."""
    return (json.dumps(obj, sort_keys=True, indent=2, ensure_ascii=True) + "\n").encode("ascii")


def load_pinned_json(path: Path, pin: str, what: str) -> dict:
    raw = path.read_bytes()
    got = sha256_bytes(raw)
    if got != pin:
        raise Refusal(f"{what}: SHA-256 {got} != pinned {pin}; not parsed")
    return json.loads(raw)


def need(obj, key: str, ctx: str):
    """Checked access on packet data: a missing key is a Refusal, never a traceback."""
    if not isinstance(obj, dict):
        raise Refusal(f"{ctx}: expected an object, got {type(obj).__name__}")
    if key not in obj:
        raise Refusal(f"{ctx}: missing key {key!r}")
    return obj[key]


def load_answer_packet_bytes(raw: bytes, pin: str = PIN_ANSWER_PACKET) -> tuple:
    """Step 6.  Returns (packet, load_record).

    The pin is the section 10.2 hash of the CANONICAL form.  The delivered bytes are hashed
    first; if they equal the pin nothing further is needed.  If they do not, section 10.2
    allows canonicalization of the delivered rendering, so the bytes are parsed ONLY to
    re-serialize them canonically, and nothing is consumed unless the canonical hash equals
    the pin; the record states that canonicalization changed the bytes.  Any other outcome
    is a refusal before the object reaches the comparison.
    """
    delivered = sha256_bytes(raw)
    rec = {"delivered_sha256": delivered, "delivered_length": len(raw)}

    def parse(b: bytes):
        # Duplicate keys are refused outright: json.loads would keep the last occurrence,
        # so a decoy rendering could otherwise be "canonicalized" into the pinned object.
        def no_dupes(pairs):
            d = {}
            for k, v in pairs:
                if k in d:
                    raise Refusal(f"answer packet: duplicate key {k!r} in delivered bytes")
                d[k] = v
            return d

        try:
            return json.loads(b, object_pairs_hook=no_dupes)
        except Refusal:
            raise
        except Exception as e:  # ValueError, RecursionError, UnicodeDecodeError, ...
            raise Refusal(f"answer packet: bytes do not parse: {type(e).__name__}: {e}")

    if delivered == pin:
        # The pin proves the object.  A rendering that is not byte-identical to this
        # harness's canonicalizer is RECORDED, not refused: the hash already matched.
        packet = parse(raw)
        rec["canonicalization_changed_bytes"] = False
        rec["canonical_form_ok"] = canonical_bytes(packet) == raw
        rec["canonical_sha256"] = pin
    else:
        packet = parse(raw)
        canon = canonical_bytes(packet)
        canon_hash = sha256_bytes(canon)
        if canon_hash != pin:
            raise Refusal(
                f"answer packet: delivered SHA-256 {delivered} != pin {pin}, and the canonical "
                f"form hashes to {canon_hash}, also != pin; not consumed"
            )
        rec["canonicalization_changed_bytes"] = True
        rec["canonical_form_ok"] = True
        rec["canonical_sha256"] = canon_hash
        rec["canonical_length"] = len(canon)
    if not isinstance(packet, dict) or set(packet) != PACKET_KEYS:
        raise Refusal(
            f"answer packet: key set "
            f"{sorted(packet) if isinstance(packet, dict) else type(packet).__name__} "
            f"!= {sorted(PACKET_KEYS)}"
        )
    for k in ("format_version", "target_id"):
        if not isinstance(packet[k], str) or not packet[k]:
            raise Refusal(f"answer packet: {k} must be a non-empty string")
    return packet, rec


def signature_key(sig: dict, fields: tuple = SIG_FIELDS, what: str = "signature") -> tuple:
    """The section 5.5 row signature as one hashable key, whichever spelling carries it.

    The answer packet spells it `signature{dimension, s, t, st}` and the committed raw
    output `row_signature{dim, chi_s, chi_t, chi_st}`; both reduce to
    (dim, chi_s, chi_t, chi_st) with the characters as normalized triples.
    """
    if not isinstance(sig, dict) or set(sig) != set(fields):
        raise Refusal(
            f"{what} fields {sorted(sig) if isinstance(sig, dict) else sig!r} != {sorted(fields)}"
        )
    dim = sig[fields[0]]
    if not isinstance(dim, int) or isinstance(dim, bool) or dim <= 0:
        raise Refusal(f"{what}.{fields[0]} not a positive int: {dim!r}")
    return (dim,) + tuple(tuple(QPhi.from_triple(sig[f]).to_triple()) for f in fields[1:])


def raw_key(r: dict) -> tuple:
    return signature_key(r["row_signature"], RAW_SIG_FIELDS, "row_signature")


def load_raw_rows(raw_output: dict) -> dict:
    """Committed rows keyed by signature -> {'acyclic': bool, 'value': QPhi|None}."""
    rows = {}
    for r in raw_output["rows"]:
        key = raw_key(r)
        if key in rows:
            raise Refusal(f"raw output: duplicate row signature {key}")
        acyclic = bool(r.get("acyclic"))
        val = QPhi.from_triple(r["T_squared_native"]) if acyclic else None
        if acyclic and val == QPhi(0, 0):
            raise Refusal(f"raw output: acyclic row {key} carries zero torsion")
        rows[key] = {"acyclic": acyclic, "value": val}
    if len(rows) != N_ROWS:
        raise Refusal(f"raw output: {len(rows)} rows != {N_ROWS}")
    return rows


def committed_order(raw_output: dict) -> list:
    """The committed rows' signature keys in FILE order.  NOT a destination domain: the
    author's map indexes the M8.3 record order, and reading it against the realized array
    position is the reading #451 (2026-08-21 14:28) ruled out.  Kept for the self-test
    arms (C3l, C3m) only."""
    return [raw_key(r) for r in raw_output["rows"]]


# ---- the M8.3 record order, from the two PUBLIC packets alone ------------------------
# The record's rows are R0..R8, labels fixed PUBLICLY before M8.8 by dimension and McKay
# distance (m8_2_first_occurrence.py: dims 1,2,2,3,3,4,4,5,6; distances 0,1,7,2,6,6,3,4,5).
# That table is invariant under the Galois flip (R1 <-> R2, R3 <-> R4), so it fixes the
# labels only up to the outer automorphism; what splits each Galois pair is section 5.5's
# own rule, "through the group packet's embedding": R1 IS the embedding (character = the
# quaternion trace), R2 its Galois conjugate.  With that rule the signatures follow with no
# further input: the symmetric-power tower V_n (dim n+1) has
# chi_n = chi_1 chi_{n-1} - chi_{n-2}, and R0..R8 read
#   R0 = V_0, R1 = V_1, R2 = sigma(V_1), R3 = V_2, R4 = sigma(V_2), R5 = V_1 (x) sigma(V_1),
#   R6 = V_3, R7 = V_4, R8 = V_5,
# sigma the Galois map phi -> 1 - phi on characters (the two dim-4 irreps are Galois-fixed
# and are told apart by distance: R6 is the tower's V_3, R5 the remaining one).  Element
# IDs follow the protocol's enumeration (section 4.1): rank of the eight numerator integers
# with the fixed denominator 2 dropped, lexicographic.  The derived keys must each be a
# committed signature, or the live run refuses: the public derivation and the committed
# output disagreeing is a structural failure, never a silent pairing.
LABEL_ORDER = ("R0", "R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8")
_PUBLIC_CACHE: dict = {}


def _parse_packet_coordinate(text: str) -> QPhi:
    """`(A + B*phi)/2`, the group packet's fixed-denominator form, as a Q(phi) value."""
    s = text.strip()
    if not (s.startswith("(") and s.endswith(")/2")):
        raise Refusal(f"group packet coordinate not in (A + B*phi)/2 form: {text!r}")
    inner = s[1:-3].replace(" ", "")
    head, _, tail = inner.partition("*phi")
    if tail or not head:
        raise Refusal(f"group packet coordinate not in (A + B*phi)/2 form: {text!r}")
    cut = max(head.rfind("+"), head.rfind("-"))
    if cut <= 0:
        raise Refusal(f"group packet coordinate not in (A + B*phi)/2 form: {text!r}")
    return QPhi(Fraction(int(head[:cut]), 2), Fraction(int(head[cut:]), 2))


def _qmul(p: tuple, q: tuple) -> tuple:
    w1, x1, y1, z1 = p
    w2, x2, y2, z2 = q
    return (
        w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
        w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
        w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
        w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2,
    )


def _enumeration_key(q: tuple) -> tuple:
    """Section 4.1 ID sort key: the eight numerators over the fixed denominator 2."""
    return tuple(int(v * 2) for c in q for v in (c.x, c.y))


def public_label_signatures(group_packet: dict, construction_packet: dict) -> dict:
    """label -> signature key (dim, chi_s, chi_t, chi_st), from the public packets only."""
    if group_packet.get("quaternion_basis") != ["1", "i", "j", "k"]:
        raise Refusal("group packet: quaternion_basis is not [1, i, j, k]")
    gens = [tuple(_parse_packet_coordinate(c) for c in g) for g in group_packet["generators"]]
    if any(len(g) != 4 for g in gens):
        raise Refusal("group packet: a generator is not a 4-component quaternion")
    zero = QPhi(0, 0)
    one = tuple([QPhi(1, 0), zero, zero, zero])
    elems = {one, *gens}
    negs = [tuple(zero - c for c in g) for g in gens]
    while True:
        fresh = set()
        for a in list(elems):
            for b in gens + negs:
                for c in (_qmul(a, b), _qmul(b, a)):
                    if c not in elems:
                        fresh.add(c)
        if not fresh:
            break
        elems |= fresh
    if len(elems) != 120:
        raise Refusal(f"group packet closure has {len(elems)} elements, not 120")
    by_id = sorted(elems, key=_enumeration_key)
    ag = construction_packet["abstract_generators"]
    s, t = by_id[ag["s"]], by_id[ag["t"]]
    st = _qmul(s, t)
    chi1 = {name: g[0] * QPhi(2, 0) for name, g in (("s", s), ("t", t), ("st", st))}
    _PUBLIC_CACHE["embedding_trace"] = {k: tuple(v.to_triple()) for k, v in chi1.items()}
    tower = [{k: QPhi(1, 0) for k in chi1}, chi1]
    for n in range(2, 6):
        tower.append({k: chi1[k] * tower[n - 1][k] - tower[n - 2][k] for k in chi1})

    def sigma(ch):
        return {k: v.conjugate() for k, v in ch.items()}

    chars = {
        "R0": (1, tower[0]),
        "R1": (2, tower[1]),
        "R2": (2, sigma(tower[1])),
        "R3": (3, tower[2]),
        "R4": (3, sigma(tower[2])),
        "R5": (4, {k: tower[1][k] * tower[1][k].conjugate() for k in chi1}),
        "R6": (4, tower[3]),
        "R7": (5, tower[4]),
        "R8": (6, tower[5]),
    }
    return {
        lab: (dim,) + tuple(tuple(ch[k].to_triple()) for k in ("s", "t", "st"))
        for lab, (dim, ch) in chars.items()
    }


def m8_3_record_order(raw_output: dict) -> list:
    """`m8_3_record_row_position`: the committed signature keys in the M8.3 record's row
    order R0..R8, derived from the public packets (cached) and checked against the
    committed rows.  `raw_output` is consulted only for that presence check."""
    if "labels" not in _PUBLIC_CACHE:
        gp = load_pinned_json(GROUP_PACKET_PATH, PIN_GROUP_PACKET, "group packet")
        cp = load_pinned_json(CONSTRUCTION_PACKET_PATH, PIN_CONSTRUCTION_PACKET, "construction")
        _PUBLIC_CACHE["labels"] = public_label_signatures(gp, cp)
    labels = _PUBLIC_CACHE["labels"]
    committed = set(committed_order(raw_output))
    order = [labels[lab] for lab in LABEL_ORDER]
    missing = [lab for lab, key in zip(LABEL_ORDER, order) if key not in committed]
    if missing:
        raise Refusal(
            f"m8_3_record_row_position: the public derivation of {missing} is not among the "
            f"committed signatures; the destination order cannot be bound"
        )
    if len(set(order)) != N_ROWS:
        raise Refusal("m8_3_record_row_position: the public derivation does not separate R0..R8")
    return order


DESTINATION_ORDERS[DESTINATION_DOMAIN] = m8_3_record_order


def canonical_source_order(rows: list, own: dict) -> list:
    """The packet's canonical row order (the SOURCE domain): its rows sorted by the
    ten-integer signature key, independent of file order."""
    return sorted(range(len(rows)), key=lambda i: own[rows[i]["label"]])


def check_raw_output_structure(raw_output: dict) -> list:
    """Section 8 `structural failure` conditions visible from the committed raw output."""
    problems = []
    if raw_output.get("schema_version") != RAW_SCHEMA:
        problems.append(f"raw schema {raw_output.get('schema_version')!r} != {RAW_SCHEMA!r}")
    if raw_output.get("group_packet_sha256") != PIN_GROUP_PACKET:
        problems.append("raw output consumed a group packet other than the section 11 pin")
    if raw_output.get("construction_packet_sha256") != PIN_CONSTRUCTION_PACKET:
        problems.append("raw output consumed a construction packet other than the pin")
    if raw_output.get("manifest_sha256") != PIN_MANIFEST:
        problems.append("raw output names a manifest other than the Addendum 1 pin")
    gates = raw_output.get("gate_results", {})
    for gid, g in sorted(gates.items()):
        if g.get("outcome") != "PASS":
            problems.append(f"pre-reveal gate {gid} recorded {g.get('outcome')!r}")
    non_acyclic = [r for r in raw_output["rows"] if not r.get("acyclic")]
    if len(non_acyclic) != 1 or non_acyclic[0]["row_signature"]["dim"] != 1:
        problems.append("expected exactly one non-acyclic row, the trivial irrep R0")
    return problems


def check_phase_b(phase_b: dict, manifest_text: str) -> list:
    """Addendum 1: every declared mutation executed and red, exact registry coverage."""
    problems = []
    if not isinstance(phase_b, dict):
        return ["phase B record: not an object"]
    if not str(phase_b.get("schema_version", "")).startswith(PHASE_B_SCHEMA_PREFIX):
        problems.append(f"phase B: schema_version {phase_b.get('schema_version')!r} unexpected")
    if phase_b.get("manifest_sha256") != PIN_MANIFEST:
        problems.append("phase B: manifest_sha256 != the Addendum 1 manifest pin")
    for flag in (
        "phase_a_hashes_verified_pre",
        "phase_a_hashes_verified_post",
        "manifest_parsed_at_runtime",
        "parser_self_test_passed",
        "all_mutations_reddened",
    ):
        if phase_b.get(flag) is not True:
            problems.append(f"phase B: {flag} is not true")
    records = phase_b.get("results")
    if not isinstance(records, list) or not records:
        return problems + ["phase B record: `results` is not a non-empty list"]
    executed = {}
    for rec in records:
        if not isinstance(rec, dict):
            problems.append("phase B: a result record is not an object")
            continue
        gid = rec.get("gate_id")
        if gid in executed:
            problems.append(f"phase B: duplicate record for {gid}")
        executed[gid] = rec
        if rec.get("red_outcome") is not True:
            problems.append(f"phase B: {gid} red_outcome is not true")
        if not str(rec.get("baseline_result", "")).startswith("PASS"):
            problems.append(f"phase B: {gid} baseline_result is not a PASS")
        if not str(rec.get("mutated_result", "")).startswith("FAIL"):
            problems.append(f"phase B: {gid} mutated_result is not a FAIL (red not evidenced)")
        for field in ("object_mutated", "implemented_mutation", "declared_mutation"):
            if not isinstance(rec.get(field), str) or not rec.get(field).strip():
                problems.append(f"phase B: {gid} {field} is empty")
    cov = phase_b.get("registry_coverage")
    if not isinstance(cov, dict):
        problems.append("phase B: registry_coverage missing")
    else:
        if sorted(cov.get("executed_gate_ids") or []) != sorted(executed):
            problems.append("phase B: registry_coverage.executed_gate_ids != results")
        if cov.get("count") != len(records):
            problems.append("phase B: registry_coverage.count != number of results")
        for flag in ("pre_execution_set_equality", "post_execution_set_equality"):
            if cov.get(flag) is not True:
                problems.append(f"phase B: registry_coverage.{flag} is not true")
    declared = set()
    for line in manifest_text.splitlines():
        if line.startswith("| G-") and "|" in line[1:]:
            declared.add(line.split("|")[1].strip().split()[0])
    if declared != set(executed):
        problems.append(
            f"phase B: executed set != manifest registry; "
            f"missing {sorted(declared - set(executed))}, "
            f"extra {sorted(set(executed) - declared)}"
        )
    return problems


# ----------------------------------------------------------------------------------------
# Step 7: adjudicates, indexing map, convention map.
# ----------------------------------------------------------------------------------------
def check_adjudicates(packet: dict, raw_output: dict) -> None:
    adj = packet["adjudicates"]
    if not isinstance(adj, dict) or set(adj) != {
        "group_packet_sha256",
        "construction_packet_sha256",
    }:
        raise Refusal("adjudicates: wrong key set")
    if adj["group_packet_sha256"] != PIN_GROUP_PACKET:
        raise Refusal("adjudicates.group_packet_sha256 != section 11 pin")
    if adj["construction_packet_sha256"] != PIN_CONSTRUCTION_PACKET:
        raise Refusal("adjudicates.construction_packet_sha256 != section 11 pin")
    if adj["group_packet_sha256"] != raw_output["group_packet_sha256"]:
        raise Refusal("adjudicates: group packet != the one the run consumed")
    if adj["construction_packet_sha256"] != raw_output["construction_packet_sha256"]:
        raise Refusal("adjudicates: construction packet != the one the run consumed")


def validate_convention_map(packet: dict, construction_packet: dict) -> dict:
    cm = packet["convention_map"]
    if not isinstance(cm, dict) or set(cm) != CONVENTION_KEYS:
        raise Refusal(
            f"convention_map: key set {sorted(cm) if isinstance(cm, dict) else cm!r}"
            f" != {sorted(CONVENTION_KEYS)}"
        )
    for k, v in cm.items():
        if not isinstance(v, dict) or not v:
            raise Refusal(f"convention_map.{k}: not a non-empty object")
        if set(v) != CONVENTION_ITEM_KEYS[k]:
            raise Refusal(
                f"convention_map.{k}: key set {sorted(v)} != {sorted(CONVENTION_ITEM_KEYS[k])}"
            )
    # Section 4.4: the fourth item reproduces the construction packet's `basing.evaluation`
    # VERBATIM and nothing else may be filled in against it.  The field is `convention`
    # (author-confirmed, #451): byte-equal to the declaration, or refused.
    declared = construction_packet["basing"]["evaluation"]
    ev = cm["representation_evaluation"]
    if ev[EVALUATION_FIELD] != declared:
        raise Refusal(
            f"convention_map.representation_evaluation.{EVALUATION_FIELD} is not byte-equal "
            "to the construction packet's basing.evaluation verbatim"
        )
    # The other fields and the other three items are prose the protocol does not pin
    # verbatim.  A free-text check cannot tell a wrong bridge from a paraphrase (an earlier
    # heuristic here refused section 5.4's own wording), so they are RECORDED for the
    # adjudication record and judged by the adjudicator there, never gated by substring.
    return {
        "evaluation_verbatim_gated": True,
        "other_items_gated": False,
        "orientation_applied": False,
        "representation_evaluation": ev,
        "evaluation_verbatim_fields": [EVALUATION_FIELD],
        "bridge": cm["bridge"],
        "orientation_anchor_rule": cm["orientation_anchor_rule"],
        "basing_reference": cm["basing_reference"],
    }


def identity_indexing_map(zero_based: bool = True) -> dict:
    """The trivial positional map, used ONLY by self-test arms that simulate ignoring the
    packet's map (pairing canonical packet row i with M8.3 record row i)."""
    base = 0 if zero_based else 1
    return {
        "source_domain": SOURCE_DOMAIN,
        "destination_domain": DESTINATION_DOMAIN,
        "zero_based": zero_based,
        "entries": [
            {"source_position": i + base, "destination_position": i + base} for i in range(N_ROWS)
        ],
    }


def apply_indexing_map(
    packet: dict, raw_rows: dict, raw_output: dict, map_override: dict | None = None
) -> dict:
    """Pair each packet row with ONE committed row BY THE POSITIONAL MAP.  Returns
    label -> paired record.

    `map_override` is used ONLY by the self-test arms (identity arm, wrong-map arm); the
    live run always consumes the packet's own map.
    """
    imap = map_override if map_override is not None else packet["indexing_map"]
    if not isinstance(imap, dict) or set(imap) != INDEXING_KEYS:
        raise Refusal(
            f"indexing_map: key set {sorted(imap) if isinstance(imap, dict) else imap!r}"
            f" != {sorted(INDEXING_KEYS)}"
        )
    # FAIL CLOSED on the domains: an unrecognized domain string is refused, so the map is
    # never applied under a guessed reading.
    for k, allowed in (
        ("source_domain", SOURCE_DOMAINS),
        ("destination_domain", set(DESTINATION_ORDERS)),
    ):
        if not isinstance(imap[k], str) or imap[k] not in allowed:
            raise Refusal(
                f"indexing_map.{k} {imap[k]!r} is not an implemented domain {sorted(allowed)}; "
                f"the map is not applied"
            )
    zero_based = imap["zero_based"]
    if not isinstance(zero_based, bool):
        raise Refusal("indexing_map.zero_based must be a JSON boolean")
    base = 0 if zero_based else 1
    raw_order = DESTINATION_ORDERS[imap["destination_domain"]](raw_output)

    rows = packet["rows"]
    if not isinstance(rows, list) or len(rows) != N_ROWS:
        raise Refusal(
            f"packet rows: {len(rows) if isinstance(rows, list) else rows!r} != {N_ROWS}"
        )
    for r in rows:
        if not isinstance(r, dict) or set(r) != ROW_KEYS:
            raise Refusal(
                f"packet row: key set {sorted(r) if isinstance(r, dict) else r!r} != {sorted(ROW_KEYS)}"
            )
    labels = [need(r, "label", "packet row") for r in rows]
    if len(set(labels)) != N_ROWS or any(not isinstance(lab, str) for lab in labels):
        raise Refusal("packet rows: labels not nine distinct strings")
    # The packet's own signatures: pairwise distinct and each present among the committed
    # signatures, so a packet cannot carry rows that describe nothing.  They do NOT pair.
    own = {r["label"]: signature_key(need(r, "signature", f"row {r['label']}")) for r in rows}
    if len(set(own.values())) != N_ROWS:
        raise Refusal("packet rows: signatures not pairwise distinct")
    for lab, key in own.items():
        if key not in raw_rows:
            raise Refusal(f"row {lab}: own signature {key} absent from the committed output")

    entries = imap["entries"]
    if not isinstance(entries, list) or len(entries) != N_ROWS:
        raise Refusal("indexing_map.entries: not nine entries")
    src_to_dst = {}
    for e in entries:
        if not isinstance(e, dict) or set(e) != ENTRY_KEYS:
            raise Refusal(f"indexing_map entry: key set != {sorted(ENTRY_KEYS)}")
        for name in ENTRY_KEYS:
            if not isinstance(e[name], int) or isinstance(e[name], bool):
                raise Refusal(f"indexing_map entry: {name} must be an int")
        s, d = e["source_position"] - base, e["destination_position"] - base
        if not 0 <= s < N_ROWS:
            raise Refusal(
                f"indexing_map entry: source_position {e['source_position']} outside the "
                f"packet rows (zero_based={zero_based})"
            )
        if not 0 <= d < N_ROWS:
            raise Refusal(
                f"indexing_map entry: destination_position {e['destination_position']} "
                f"outside the committed rows (zero_based={zero_based})"
            )
        if s in src_to_dst:
            raise Refusal(
                f"indexing_map.entries: source_position {e['source_position']} listed twice"
            )
        src_to_dst[s] = d
    # nine distinct in-range source positions cover the nine packet rows exactly once

    # SOURCE positions index the packet's CANONICAL row order (rows sorted by the
    # ten-integer signature key), never the file order (author, #451).
    source_order = canonical_source_order(rows, own)
    paired, used = {}, set()
    for pos, file_pos in enumerate(source_order):
        r = rows[file_pos]
        dst = src_to_dst[pos]
        key = raw_order[dst]
        if key in used:
            raise Refusal(f"indexing map is not injective at committed position {dst + base}")
        used.add(key)
        if map_override is None and own[r["label"]] != key:
            # Section 5.5 makes the signature the row's identity.  A LIVE packet whose map
            # pairs a row with a committed row of another signature contradicts itself and
            # is refused, never silently categorized.  The self-test arms (map_override)
            # keep the difference as a record so ignoring the map fails AT THE COMPARISON.
            raise Refusal(
                f"row {r['label']}: own signature {own[r['label']]} differs from the "
                f"committed row the map pairs it with, {key}"
            )
        cls = need(r, "evidentiary_class", f"row {r['label']}")
        if cls not in CLASSES:
            raise Refusal(f"row {r['label']}: evidentiary_class {cls!r} not in {sorted(CLASSES)}")
        paired[r["label"]] = {
            "signature": key,
            "position": pos + base,
            "file_position": file_pos + base,
            "committed_position": dst + base,
            "own_signature_differs_from_map": own[r["label"]] != key,
            "label_is_record_label": r["label"] == LABEL_ORDER[dst],
            "class": cls,
            "reference": QPhi.from_triple(need(r, "value", f"row {r['label']}")),
            "committed": raw_rows[key]["value"],
            "acyclic": raw_rows[key]["acyclic"],
        }
    classes = [p["class"] for p in paired.values()]
    if classes.count("free_orientation_selector") != 1:
        raise Refusal("packet rows: exactly one free_orientation_selector (R7) required")
    if classes.count("declared_convention") != 1:
        raise Refusal("packet rows: exactly one declared_convention (R0) required")
    for label, p in paired.items():
        if p["class"] == "free_orientation_selector" and p["signature"][0] != R7_DIM:
            raise Refusal(
                f"{label}: the orientation selector is paired with a dim-{p['signature'][0]} "
                f"row; section 5.4 anchors at R7, the dim-{R7_DIM} irrep"
            )
        if p["class"] == "declared_convention":
            if p["acyclic"]:
                raise Refusal(
                    f"{label}: declared_convention row is acyclic in the committed output"
                )
            if p["reference"] != ONE:
                raise Refusal(f"{label}: declared_convention value is not 1 (section 9 item 2)")
        elif not p["acyclic"]:
            # A NONTRIVIAL row the run found non-acyclic: section 8 hypothesis failure.
            p["hypothesis_failure"] = True
    return paired


# ----------------------------------------------------------------------------------------
# Step 8: the R7 selection, the row comparison, the identities.
# ----------------------------------------------------------------------------------------
def select_orientation(paired: dict) -> dict:
    ((label, p),) = [
        (lab, p) for lab, p in paired.items() if p["class"] == "free_orientation_selector"
    ]
    x, r = p["committed"], p["reference"]
    if x is None:
        return {
            "anchor": label,
            "category": "hypothesis failure",
            "reason": "R7 is non-acyclic in the committed output",
        }
    direct = x == r
    inverse = x.inverse() == r
    if direct and inverse:
        return {
            "anchor": label,
            "category": "invalid anchor",
            "orientation": None,
            "reason": "both x = r and x^-1 = r hold; the reference is self-inverse",
        }
    if not direct and not inverse:
        return {
            "anchor": label,
            "category": "disagreement",
            "orientation": None,
            "reason": "neither x = r nor x^-1 = r holds at R7",
        }
    return {
        "anchor": label,
        "category": None,
        "orientation": "identity" if direct else "global inverse",
    }


def selected_rows(paired: dict, orientation: str) -> dict:
    """The committed table under the selected orientation, by signature.  R0 carries the
    declared convention value 1 (section 9 item 2), never a computed one."""
    out = {}
    for p in paired.values():
        if p["class"] == "declared_convention":
            out[p["signature"]] = ONE
        elif p["committed"] is not None:
            v = p["committed"]
            out[p["signature"]] = v if orientation == "identity" else v.inverse()
    return out


def compare_rows(paired: dict, sel: dict) -> list:
    mismatches = []
    for label, p in sorted(paired.items()):
        if p["class"] == "declared_convention":
            continue
        got = sel.get(p["signature"])
        if got is None or got != p["reference"]:
            mismatches.append(
                {
                    "label": label,
                    "signature": list(p["signature"]),
                    "reference": p["reference"].to_triple(),
                    "selected_committed": got.to_triple() if got else None,
                }
            )
    return mismatches


def compute_identities(packet: dict, sel: dict) -> list:
    ids = packet["identities"]
    if not isinstance(ids, list) or len(ids) != N_IDENTITIES:
        raise Refusal(
            f"identities: {len(ids) if isinstance(ids, list) else ids!r} != {N_IDENTITIES}"
        )
    trivial_keys = {k for k, v in sel.items() if k[0] == 1}
    out, kinds, product_members, ratio_pairs = [], [], [], set()
    for slot_entry in ids:
        if not isinstance(slot_entry, dict) or set(slot_entry) != IDENTITY_KEYS:
            raise Refusal(f"identity entry: key set != {sorted(IDENTITY_KEYS)}")
        position = need(slot_entry, "position", "identity")
        if not isinstance(position, int) or isinstance(position, bool) or position < 0:
            raise Refusal(f"identities: position must be a non-negative int, got {position!r}")
        slot = f"position {position}"
        factors = need(slot_entry, "factors", f"identity {slot}")
        if not isinstance(factors, list) or len(factors) < 2:
            raise Refusal(f"identities: {slot} malformed or with fewer than two factors")
        acc, keys, exps = ONE, [], []
        for f in factors:
            # EXACT key set: `conjugate` is outside the frozen factor schema (author, #451),
            # and since it would change the computed identity it is refused, not ignored.
            if not isinstance(f, dict) or set(f) != FACTOR_KEYS:
                raise Refusal(
                    f"identity {slot} factor: key set "
                    f"{sorted(f) if isinstance(f, dict) else f!r} != {sorted(FACTOR_KEYS)}"
                )
            key = signature_key(need(f, "signature", f"identity {slot} factor"))
            if key not in sel:
                raise Refusal(f"identity {slot}: factor signature {key} not among selected rows")
            if key in trivial_keys:
                raise Refusal(f"identity {slot}: the declared-convention row is not a factor")
            e = need(f, "exponent", f"identity {slot} factor")
            if not isinstance(e, int) or isinstance(e, bool) or e == 0:
                raise Refusal(f"identity {slot}: exponent must be a nonzero int")
            acc = acc * (sel[key] ** e)
            keys.append(key)
            exps.append(e)
        if len(set(keys)) != len(keys):
            raise Refusal(f"identity {slot}: repeated factor")
        # Section 5.1 names the four: two Galois ratios, two sector products.  A ratio is
        # two factors with exponents {+1, -1} on a Galois pair, same dim and conjugate
        # characters; a product is two or more factors with every exponent +1.
        if len(keys) == 2 and sorted(exps) == [-1, 1]:
            a, b = keys
            conj = (a[0],) + tuple(
                tuple(QPhi.from_triple(list(t)).conjugate().to_triple()) for t in a[1:]
            )
            if b != conj:
                raise Refusal(f"identity {slot}: ratio factors are not a Galois pair")
            pair = frozenset(keys)
            if pair in ratio_pairs:
                raise Refusal(f"identity {slot}: the same Galois pair fills two ratio slots")
            ratio_pairs.add(pair)
            kinds.append("ratio")
        elif all(e == 1 for e in exps):
            kinds.append("product")
            product_members.append(set(keys))
        else:
            raise Refusal(f"identity {slot}: neither a Galois ratio nor a sector product")
        expected = QPhi.from_triple(need(slot_entry, "expected_value", f"identity {slot}"))
        out.append(
            {
                "position": position,
                "slot": slot,
                "kind": kinds[-1],
                "recomputed": acc.to_triple(),
                "expected": expected.to_triple(),
                "equal": acc == expected,
            }
        )
    if len({i["position"] for i in out}) != N_IDENTITIES:
        raise Refusal("identities: positions not distinct")
    if kinds.count("ratio") != 2 or kinds.count("product") != 2:
        raise Refusal(f"identities: need two Galois ratios and two sector products, got {kinds}")
    if product_members[0] & product_members[1]:
        raise Refusal("identities: the two sector products are not disjoint")
    # Coverage of the eight nontrivial rows by the two sector products is RECORDED, not
    # gated: the sector partition is an M8.3 structural assignment the protocol places out
    # of scope (section 5.3), so the packet declares it and the record shows it.
    covered = product_members[0] | product_members[1]
    nontrivial = {k for k in sel if k[0] != 1}
    out.append(
        {
            "position": None,
            "slot": "_sector_coverage",
            "kind": "record",
            "rows_in_products": len(covered),
            "nontrivial_rows": len(nontrivial),
            "uncovered": [list(k) for k in sorted(nontrivial - covered)],
            "equal": True,
        }
    )
    return out


def rows_record(paired: dict) -> dict:
    """The pairing as adjudicated, per label, for the written record."""
    return {
        lab: {
            "signature": list(p["signature"]),
            "position": p["position"],
            "file_position": p["file_position"],
            "committed_position": p["committed_position"],
            "class": p["class"],
            "own_signature_differs_from_map": p["own_signature_differs_from_map"],
            "label_is_record_label": p["label_is_record_label"],
        }
        for lab, p in sorted(paired.items())
    }


def adjudicate(
    packet: dict, raw_output: dict, construction_packet: dict, map_override: dict | None = None
) -> dict:
    """Steps 7 and 8 on already-loaded, already-hash-verified objects."""
    check_adjudicates(packet, raw_output)
    raw_rows = load_raw_rows(raw_output)
    cm = validate_convention_map(packet, construction_packet)
    paired = apply_indexing_map(packet, raw_rows, raw_output, map_override)

    hyp = [lab for lab, p in paired.items() if p.get("hypothesis_failure")]
    if hyp:
        return {
            "category": "hypothesis failure",
            "convention_map": cm,
            "rows": rows_record(paired),
            "non_acyclic_nontrivial": hyp,
            "orientation": None,
            "row_mismatches": None,
            "identities": None,
        }

    selection = select_orientation(paired)
    if selection["category"]:
        return {
            "category": selection["category"],
            "convention_map": cm,
            "selection": selection,
            "rows": rows_record(paired),
            "orientation": None,
            "row_mismatches": None,
            "identities": None,
        }

    sel = selected_rows(paired, selection["orientation"])
    mismatches = compare_rows(paired, sel)
    identities = compute_identities(packet, sel)
    ids_ok = all(i["equal"] for i in identities)
    if not mismatches and ids_ok:
        category = (
            "reproduced" if selection["orientation"] == "identity" else "convention difference"
        )
    else:
        category = "partial disagreement"
    return {
        "category": category,
        "convention_map": cm,
        "selection": selection,
        "rows": rows_record(paired),
        "orientation": selection["orientation"],
        "row_mismatches": mismatches,
        "identities": identities,
        "selected_table": {str(list(k)): v.to_triple() for k, v in sorted(sel.items())},
    }


# ----------------------------------------------------------------------------------------
# The self-test: section 4.4 controls on synthetic fixtures built from public data only.
# ----------------------------------------------------------------------------------------
def _sig_dict(key: tuple) -> dict:
    return {"dimension": key[0], "s": list(key[1]), "t": list(key[2]), "st": list(key[3])}


def build_fixture(
    raw_output: dict,
    construction_packet: dict,
    *,
    invert: bool = False,
    zero_based: bool = True,
    dishonest_signatures: bool = False,
    file_order: str = "record",
) -> dict:
    """A synthetic answer packet IN THE AUTHOR-CONFIRMED SHAPE whose reference table IS the
    committed table (or its global inverse).  It contains no reference value and proves
    only harness behavior.

    Labels are the public M8.3 labels (`m8_3_record_order`).  The positional map pairs
    the packet's CANONICAL row order (rows sorted by signature key) with the M8.3 record
    order; on the committed table the two orders differ, so the honest map is a genuine
    non-identity permutation (section 4.4's nonidentity fixture is the real relation, not
    a synthetic shuffle) and a harness that ignores the map fails AT THE COMPARISON on
    exactly the rows the permutation moves.  `file_order` lists the rows in record order
    (`record`), reversed (`reversed`) or canonical order (`canonical`): the map must not
    depend on it.  With dishonest_signatures=True two free rows exchange their own
    `signature` while the map stays true: a LIVE packet of that shape is refused, and the
    same fixture as a self-test arm records `own_signature_differs_from_map` on those rows.
    """
    raw_rows = load_raw_rows(raw_output)
    order = m8_3_record_order(raw_output)  # record position -> committed signature key
    label_of = dict(zip(order, LABEL_ORDER))
    trivial = order[0]
    r7 = order[LABEL_ORDER.index("R7")]
    free_keys = [k for k in order if k not in (trivial, r7)]
    base = 0 if zero_based else 1

    def value_of(k):
        if k == trivial:
            return ONE
        v = raw_rows[k]["value"]
        return v.inverse() if invert else v

    def cls_of(lab):
        if lab == "R0":
            return "declared_convention"
        return "free_orientation_selector" if lab == "R7" else "free"

    # own signature shown per record position; dishonest = the two dim-4 rows exchanged
    shown = dict(zip(order, order))
    if dishonest_signatures:
        a, b = [k for k in free_keys if k[0] == 4]
        shown[a], shown[b] = b, a
    file_positions = list(range(N_ROWS))
    if file_order == "reversed":
        file_positions.reverse()
    elif file_order == "canonical":
        file_positions.sort(key=lambda i: shown[order[i]])
    elif file_order != "record":
        raise ValueError(file_order)
    rows = []
    for rec_pos in file_positions:
        key = order[rec_pos]
        rows.append(
            {
                "label": label_of[key],
                "signature": _sig_dict(shown[key]),
                "evidentiary_class": cls_of(label_of[key]),
                "value": value_of(key).to_triple(),
            }
        )
    # the map: canonical rank of each row (by its SHOWN signature) -> record position
    canon = sorted(range(N_ROWS), key=lambda i: shown[order[i]])
    entries = [
        {"source_position": rank + base, "destination_position": rec_pos + base}
        for rank, rec_pos in enumerate(canon)
    ]

    # Four identity positions: two Galois-pair ratios (dims 2 and 3 pairs), two sector
    # products over disjoint halves.  Membership here is synthetic; the live packet
    # declares its own.  Expected values are recomputed from the same table, so the
    # fixture is self-consistent by construction.
    def ident(position, spec):
        acc = ONE
        for k, e in spec:
            acc = acc * (value_of(k) ** e)
        return {
            "position": position,
            "factors": [{"signature": _sig_dict(k), "exponent": e} for k, e in spec],
            "expected_value": acc.to_triple(),
        }

    dim2 = [k for k in free_keys if k[0] == 2]
    dim3 = [k for k in free_keys if k[0] == 3]
    half_a, half_b = free_keys[:3], free_keys[3:6]
    identities = [
        ident(base + 0, [(dim2[0], 1), (dim2[1], -1)]),
        ident(base + 1, [(dim3[0], 1), (dim3[1], -1)]),
        ident(base + 2, [(k, 1) for k in half_a]),
        ident(base + 3, [(k, 1) for k in half_b]),
    ]
    return {
        "format_version": "m8_8-answer-FIXTURE-3",
        "target_id": "SELF-TEST FIXTURE, NOT THE CANONICAL PACKET",
        "adjudicates": {
            "group_packet_sha256": PIN_GROUP_PACKET,
            "construction_packet_sha256": PIN_CONSTRUCTION_PACKET,
        },
        "rows": rows,
        "identities": identities,
        "indexing_map": {
            "source_domain": SOURCE_DOMAIN,
            "destination_domain": DESTINATION_DOMAIN,
            "zero_based": zero_based,
            "entries": entries,
        },
        "convention_map": {
            "bridge": {
                "statement": "T^2_target(rho) := |tau_rho|^2; involution T^2 <-> 1/T^2",
                "identification": "FIXTURE",
                "exclusivity": "FIXTURE",
                "protocol_section": "5.4",
            },
            "orientation_anchor_rule": {
                "rule": "select at R7 between the committed table and its global inverse",
                "anchor_row_label": "R7",
                "anchor_row_signature": _sig_dict(r7),
                "selection_mapping": "FIXTURE",
                "uniqueness_gate": "FIXTURE",
                "protocol_section": "5.4",
            },
            "basing_reference": {
                "construction_packet_sha256": PIN_CONSTRUCTION_PACKET,
                "source": "construction packet `basing`",
                "protocol_section": "4.2",
            },
            "representation_evaluation": {
                EVALUATION_FIELD: construction_packet["basing"]["evaluation"],
                "source": "construction packet `basing.evaluation`",
                "protocol_section": "4.2",
            },
        },
    }


def _run_fixture(packet, raw_output, construction_packet, **kw):
    try:
        return adjudicate(packet, raw_output, construction_packet, **kw)
    except Refusal as e:
        return {"category": "structural failure", "refusal": str(e)}
    except Exception as e:  # a crash is a structural failure with a record, never silence
        return {"category": "structural failure", "exception": f"{type(e).__name__}: {e}"}


def self_test(raw_output: dict, construction_packet: dict) -> list:
    results = []

    def record(name, expect, observed, detail=None):
        results.append(
            {
                "control": name,
                "expected": expect,
                "observed": observed,
                "ok": observed == expect,
                **({"detail": detail} if detail else {}),
            }
        )

    def sel_by_class(pkt, cls):
        return [r for r in pkt["rows"] if r["evidentiary_class"] == cls][0]

    def by_label(pkt, lab):
        return [r for r in pkt["rows"] if r["label"] == lab][0]

    def unequal_positions(out):
        return sorted(
            i["position"]
            for i in (out.get("identities") or [])
            if not i["equal"] and i["position"] is not None
        )

    # Baseline: the committed table against itself must be `reproduced`.
    fx = build_fixture(raw_output, construction_packet)
    base = _run_fixture(fx, raw_output, construction_packet)
    record("C0 baseline identity fixture", "reproduced", base["category"])

    # C1: EACH nontrivial reference cell mutated in turn on an in-memory copy, downstream
    # of hashing; the comparison must go red with exactly that one label mismatching, so a
    # comparison that silently skipped any row would be caught here.
    nontrivial = [
        r["label"] for r in fx["rows"] if r["evidentiary_class"] != "declared_convention"
    ]
    for lab in nontrivial:
        m = copy.deepcopy(fx)
        row = by_label(m, lab)
        row["value"] = (QPhi.from_triple(row["value"]) * QPhi(2, 0)).to_triple()
        out = _run_fixture(m, raw_output, construction_packet)
        expected_cat = "disagreement" if lab == "R7" else "partial disagreement"
        flagged = sorted(x["label"] for x in (out.get("row_mismatches") or []))
        record(
            f"C1 reference cell {lab} mutated",
            (expected_cat, [] if lab == "R7" else [lab]),
            (out["category"], flagged),
        )

    # C2: EACH committed cell mutated the same way, from the other side.
    for r0 in [r for r in raw_output["rows"] if r.get("acyclic")]:
        mraw = copy.deepcopy(raw_output)
        row = [r for r in mraw["rows"] if r["row_signature"] == r0["row_signature"]][0]
        row["T_squared_native"] = (
            QPhi.from_triple(row["T_squared_native"]) * QPhi(3, 0)
        ).to_triple()
        out = _run_fixture(fx, mraw, construction_packet)
        is_r7 = r0["row_signature"]["dim"] == R7_DIM
        n_flagged = len(out.get("row_mismatches") or [])
        record(
            f"C2 committed cell dim-{r0['row_signature']['dim']} mutated",
            ("disagreement", 0) if is_r7 else ("partial disagreement", 1),
            (out["category"], n_flagged),
        )

    # C3: the NONIDENTITY positional fixture (section 4.4).  The author's map pairs the
    # packet's CANONICAL row order (sorted by signature key) with the M8.3 record order;
    # on the committed table the two orders differ, so the honest map itself is the
    # genuine non-identity permutation and no synthetic shuffle is needed.  Applying it
    # must be GREEN; ignoring it (the identity arm) and a preregistered wrong permutation
    # must each reach the comparison and fail THERE, on the rows the permutation moves.
    dfx = build_fixture(raw_output, construction_packet)
    moved_entries = [
        e
        for e in dfx["indexing_map"]["entries"]
        if e["source_position"] != e["destination_position"]
    ]
    n_moved = len(moved_entries)
    rec_order = m8_3_record_order(raw_output)
    canon_order = sorted(rec_order)
    n_differ = sum(1 for a, b in zip(rec_order, canon_order) if a != b)
    record(
        "C3 the honest map is a NON-IDENTITY permutation (canonical order != record order)",
        (True, n_differ),
        (n_moved > 0, n_moved),
        {"moved_source_positions": sorted(e["source_position"] for e in moved_entries)},
    )
    supplied = _run_fixture(dfx, raw_output, construction_packet)
    identity_arm = _run_fixture(
        dfx, raw_output, construction_packet, map_override=identity_indexing_map()
    )
    wrong = copy.deepcopy(dfx["indexing_map"])
    free_dst = [
        e
        for e in wrong["entries"]
        if rec_order[e["destination_position"]] not in (rec_order[0], rec_order[7])
    ][:3]
    dsts = [e["destination_position"] for e in free_dst]
    for e, d in zip(free_dst, dsts[1:] + dsts[:1]):
        e["destination_position"] = d
    wrong_arm = _run_fixture(dfx, raw_output, construction_packet, map_override=wrong)
    record("C3a the honest map, supplied (applied)", "reproduced", supplied["category"])
    record(
        "C3b the map IGNORED (identity arm): red AT THE COMPARISON on the moved rows",
        ("partial disagreement", n_moved),
        (identity_arm["category"], len(identity_arm.get("row_mismatches") or [])),
    )
    n_wrong = sum(1 for e, d in zip(free_dst, dsts) if e["destination_position"] != d)
    record(
        "C3c a preregistered wrong permutation (three free rows cycled): red AT THE COMPARISON",
        ("partial disagreement", n_wrong),
        (wrong_arm["category"], len(wrong_arm.get("row_mismatches") or [])),
    )
    # C3d: a map that moves R0 onto an acyclic row is refused on class, separately.
    r0move = copy.deepcopy(dfx["indexing_map"])
    e0 = [e for e in r0move["entries"] if e["destination_position"] == 0][0]
    e1 = [e for e in r0move["entries"] if e["destination_position"] == 1][0]
    e0["destination_position"], e1["destination_position"] = (
        e1["destination_position"],
        e0["destination_position"],
    )
    out = _run_fixture(dfx, raw_output, construction_packet, map_override=r0move)
    record(
        "C3d map moving R0 onto an acyclic row refused BY THE CLASS CHECK",
        ("structural failure", True),
        (out["category"], "is acyclic" in out.get("refusal", "")),
    )
    # C3e: a LIVE packet whose own signatures contradict the map (section 5.5 makes the
    # signature the identity) is REFUSED naming both signatures; the same fixture run as a
    # self-test arm (map_override) is green with the difference RECORDED on the two rows,
    # so the live refusal and the comparison-level arms are both exercised.
    sfx3 = build_fixture(raw_output, construction_packet, dishonest_signatures=True)
    out = _run_fixture(sfx3, raw_output, construction_packet)
    record(
        "C3e live packet with own signatures contradicting the map refused BY THAT CHECK",
        ("structural failure", True),
        (
            out["category"],
            "differs from the committed row the map pairs" in out.get("refusal", ""),
        ),
    )
    out = _run_fixture(sfx3, raw_output, construction_packet, map_override=sfx3["indexing_map"])
    record(
        "C3e' the same fixture as a self-test arm: green, differences recorded on two rows",
        ("reproduced", 2),
        (
            out["category"],
            sum(
                1 for r in (out.get("rows") or {}).values() if r["own_signature_differs_from_map"]
            ),
        ),
    )
    # C3f: zero_based = false with both positions shifted by one is the same map.
    ofx = build_fixture(raw_output, construction_packet, zero_based=False)
    out = _run_fixture(ofx, raw_output, construction_packet)
    record("C3f one-based map honored", "reproduced", out["category"])
    # C3g: the base misdeclared (one-based positions under zero_based = true) refuses on
    # range, never silently re-pairs.
    mis = copy.deepcopy(ofx)
    mis["indexing_map"]["zero_based"] = True
    out = _run_fixture(mis, raw_output, construction_packet)
    record(
        "C3g misdeclared base refused ON RANGE",
        ("structural failure", True),
        (out["category"], "outside the" in out.get("refusal", "")),
    )
    # C3h: a source position listed twice (the entries no longer cover every packet row).
    twice = copy.deepcopy(dfx)
    twice["indexing_map"]["entries"][1]["source_position"] = twice["indexing_map"]["entries"][2][
        "source_position"
    ]
    out = _run_fixture(twice, raw_output, construction_packet)
    record(
        "C3h source position listed twice refused BY THAT CHECK",
        ("structural failure", True),
        (out["category"], "listed twice" in out.get("refusal", "")),
    )

    # C3i: a destination position past the last committed row refuses ON RANGE (not an
    # IndexError); C3j: a JSON true in an entry is not an int; C3k: an unimplemented domain
    # string is refused BEFORE any pairing (the map fails closed).
    oor = copy.deepcopy(dfx)
    oor["indexing_map"]["entries"][0]["destination_position"] = N_ROWS
    out = _run_fixture(oor, raw_output, construction_packet)
    record(
        "C3i destination position past the committed rows refused ON RANGE",
        ("structural failure", True),
        (out["category"], "outside the committed rows" in out.get("refusal", "")),
    )
    boo = copy.deepcopy(dfx)
    boo["indexing_map"]["entries"][0]["destination_position"] = True
    out = _run_fixture(boo, raw_output, construction_packet)
    record(
        "C3j JSON boolean as a position refused",
        ("structural failure", True),
        (out["category"], "must be an int" in out.get("refusal", "")),
    )
    for k, bad in (
        ("source_domain", "m8_5a.label_order"),
        ("source_domain", "answer_packet.rows"),
        ("destination_domain", "RAW_OUTPUT.rows"),
    ):
        dom = copy.deepcopy(dfx)
        dom["indexing_map"][k] = bad
        out = _run_fixture(dom, raw_output, construction_packet)
        record(
            f"C3k unimplemented {k} {bad!r} refused, map not applied (fails closed)",
            ("structural failure", True),
            (out["category"], "not an implemented domain" in out.get("refusal", "")),
        )

    # C3l: the destination order is looked up BY the domain string: binding a second
    # string to an order with two free rows exchanged re-pairs exactly those two rows and
    # fails at the comparison, so a vocabulary entry cannot silently reuse another order.
    def _two_swapped(ro):
        order = m8_3_record_order(ro)
        a, b = [i for i, k in enumerate(order) if k[0] in (2, 3)][:2]
        order[a], order[b] = order[b], order[a]
        return order

    DESTINATION_ORDERS["_selftest.two_swapped"] = _two_swapped
    try:
        rev = copy.deepcopy(dfx["indexing_map"])
        rev["destination_domain"] = "_selftest.two_swapped"
        out = _run_fixture(dfx, raw_output, construction_packet, map_override=rev)
    finally:
        del DESTINATION_ORDERS["_selftest.two_swapped"]
    record(
        "C3l a second domain string maps through ITS OWN order: red at the comparison",
        ("partial disagreement", 2),
        (out["category"], len(out.get("row_mismatches") or [])),
    )

    # C3m: no two registered destination domains may yield the same order on the committed
    # table, and none may equal the committed file order, so a future entry that merely
    # aliases the realized array position fails here.
    orders = {name: tuple(fn(raw_output)) for name, fn in DESTINATION_ORDERS.items()}
    record(
        "C3m registered destination domains are pairwise distinct orders, none the file order",
        (len(orders), False),
        (len(set(orders.values())), tuple(committed_order(raw_output)) in orders.values()),
        {"domains": sorted(orders)},
    )
    # C3n: the destination order is PUBLIC-DERIVED: the nine keys are exactly the committed
    # signatures (set equality), and a committed signature altered in a scratch copy makes
    # the binding refuse rather than pair.
    record(
        "C3n the public derivation of R0..R8 reproduces the committed signature set",
        True,
        set(rec_order) == set(committed_order(raw_output)),
        {"record_order_as_file_index": [committed_order(raw_output).index(k) for k in rec_order]},
    )
    alt = copy.deepcopy(raw_output)
    alt["rows"][1]["row_signature"]["chi_st"] = [1, 0, 1]
    out = _run_fixture(dfx, alt, construction_packet)
    record(
        "C3n' a committed signature the public derivation cannot reproduce refuses the binding",
        ("structural failure", True),
        (out["category"], "cannot be bound" in out.get("refusal", "")),
    )
    # C3p: the Galois convention has an oracle.  The (dim, distance) table is blind to the
    # flip R1 <-> R2, R3 <-> R4 (an audit showed the self-test green under both swaps), so
    # the binding is checked against section 5.5's rule directly: R1's character at
    # (s, t, st) is the embedding's trace, recomputed here from the generator quaternions
    # by a second route (2 x the real part of the elements at the declared IDs), R2 is its
    # conjugate, and R3/R4 follow the same orientation (R3 = chi_1^2 - 1 on R1's side).
    gp_ = load_pinned_json(GROUP_PACKET_PATH, PIN_GROUP_PACKET, "group packet")
    gens_ = [tuple(_parse_packet_coordinate(c) for c in g) for g in gp_["generators"]]
    z_ = QPhi(0, 0)
    el_ = {tuple([QPhi(1, 0), z_, z_, z_]), *gens_}
    while True:
        new_ = {
            c
            for a in list(el_)
            for b in gens_ + [tuple(z_ - x for x in g) for g in gens_]
            for c in (_qmul(a, b), _qmul(b, a))
            if c not in el_
        }
        if not new_:
            break
        el_ |= new_
    ids_ = sorted(el_, key=_enumeration_key)
    ag_ = construction_packet["abstract_generators"]
    s_, t_ = ids_[ag_["s"]], ids_[ag_["t"]]
    trace_ = {
        "s": (s_[0] + s_[0]).to_triple(),
        "t": (t_[0] + t_[0]).to_triple(),
        "st": (_qmul(s_, t_)[0] + _qmul(s_, t_)[0]).to_triple(),
    }
    key_of_label = dict(zip(LABEL_ORDER, rec_order))
    r1, r2, r3, r4 = (key_of_label[k] for k in ("R1", "R2", "R3", "R4"))
    r1_is_embedding = [list(x) for x in r1[1:]] == [trace_["s"], trace_["t"], trace_["st"]]
    r2_is_conjugate = r2[1:] == tuple(
        tuple(QPhi.from_triple(list(x)).conjugate().to_triple()) for x in r1[1:]
    )
    r3_on_r1_side = r3[1:] == tuple(
        tuple((QPhi.from_triple(list(x)) ** 2 - ONE).to_triple()) for x in r1[1:]
    )
    r4_is_conjugate = r4[1:] == tuple(
        tuple(QPhi.from_triple(list(x)).conjugate().to_triple()) for x in r3[1:]
    )
    record(
        "C3p Galois convention oracle: R1 = the embedding's trace (section 5.5), R2/R4 conjugate, R3 on R1's side",
        (True, True, True, True),
        (r1_is_embedding, r2_is_conjugate, r3_on_r1_side, r4_is_conjugate),
        {"embedding_trace": trace_},
    )
    # C3q: the packet's label strings are RECORDED against the record position they map to,
    # never gated (labels are answer-side metadata, section 5.5): a packet with two labels
    # exchanged adjudicates the same and the record shows the two rows flagged.
    swapl = copy.deepcopy(dfx)
    a_, b_ = [r for r in swapl["rows"] if r["label"] in ("R3", "R4")]
    a_["label"], b_["label"] = b_["label"], a_["label"]
    out = _run_fixture(swapl, raw_output, construction_packet)
    record(
        "C3q exchanged label strings: same category, flagged on exactly those rows",
        ("reproduced", ["R3", "R4"]),
        (
            out["category"],
            sorted(
                l for l, r in (out.get("rows") or {}).items() if not r["label_is_record_label"]
            ),
        ),
    )
    # C3o: the packet's FILE order is irrelevant: the same rows listed reversed or in
    # canonical order pair identically and stay green.
    for fo in ("reversed", "canonical"):
        ffx = build_fixture(raw_output, construction_packet, file_order=fo)
        out = _run_fixture(ffx, raw_output, construction_packet)
        strip = lambda rows: {  # noqa: E731
            lab: {k: v for k, v in r.items() if k != "file_position"}
            for lab, r in (rows or {}).items()
        }
        record(
            f"C3o rows listed in {fo} file order: same pairing, green",
            ("reproduced", strip(supplied["rows"])),
            (out["category"], strip(out.get("rows"))),
        )

    # C4: the global inverse resolves to convention difference; the sector-product
    # positions EXCHANGE values under inversion, and the harness compares them
    # position-wise: swapping the two expected values reddens both positions, on the
    # identity and the inverted fixture.
    ifx = build_fixture(raw_output, construction_packet, invert=True)
    out = _run_fixture(ifx, raw_output, construction_packet)
    record("C4a inverted fixture", "convention difference", out["category"])
    for name, base_fx in (("identity", fx), ("inverted", ifx)):
        sw = copy.deepcopy(base_fx)
        ia = [i for i in sw["identities"] if i["position"] == 2][0]
        ib = [i for i in sw["identities"] if i["position"] == 3][0]
        ia["expected_value"], ib["expected_value"] = ib["expected_value"], ia["expected_value"]
        out = _run_fixture(sw, raw_output, construction_packet)
        record(
            f"C4b sector expected values swapped on the {name} fixture: position-wise red",
            ("partial disagreement", [2, 3], 0),
            (out["category"], unequal_positions(out), len(out.get("row_mismatches") or [])),
        )

    # C5: self-inverse anchor -> invalid anchor.  Force R7's reference AND committed to 1.
    sfx = copy.deepcopy(fx)
    sraw = copy.deepcopy(raw_output)
    sel_by_class(sfx, "free_orientation_selector")["value"] = [1, 0, 1]
    for r in sraw["rows"]:
        if r.get("acyclic") and r["row_signature"]["dim"] == R7_DIM:
            r["T_squared_native"] = [1, 0, 1]
    out = _run_fixture(sfx, sraw, construction_packet)
    record("C5 self-inverse anchor", "invalid anchor", out["category"])

    # C6: anchor equal to neither x nor x^-1 -> disagreement, no row comparison.
    nfx = copy.deepcopy(fx)
    r = sel_by_class(nfx, "free_orientation_selector")
    r["value"] = (QPhi.from_triple(r["value"]) * QPhi(7, 0)).to_triple()
    out = _run_fixture(nfx, raw_output, construction_packet)
    record(
        "C6 anchor matches neither orientation",
        "disagreement",
        out["category"],
        {"row_mismatches_issued": out.get("row_mismatches") is not None},
    )

    # C7: identity expected tampered, rows untouched -> identity layer alone goes red.
    tfx = copy.deepcopy(fx)
    tfx["identities"][2]["expected_value"] = (
        QPhi.from_triple(tfx["identities"][2]["expected_value"]) * QPhi(5, 0)
    ).to_triple()
    out = _run_fixture(tfx, raw_output, construction_packet)
    record(
        "C7 identity expected tampered",
        ("partial disagreement", 0, [2]),
        (out["category"], len(out.get("row_mismatches") or []), unequal_positions(out)),
    )

    # C8: refusals before comparison.
    afx = copy.deepcopy(fx)
    afx["adjudicates"]["construction_packet_sha256"] = "0" * 64
    out = _run_fixture(afx, raw_output, construction_packet)
    record("C8a adjudicates tampered", "structural failure", out["category"])
    good = canonical_bytes(fx)
    pin = sha256_bytes(good)
    tampered = copy.deepcopy(fx)
    tampered["rows"][1]["value"] = (
        QPhi.from_triple(tampered["rows"][1]["value"]) * QPhi(2, 0)
    ).to_triple()
    try:
        load_answer_packet_bytes(canonical_bytes(tampered), pin)
        record(
            "C8b content-tampered bytes refused (neither delivered nor canonical hash)",
            "refused",
            "accepted",
        )
    except Refusal as e:
        record(
            "C8b content-tampered bytes refused (neither delivered nor canonical hash)",
            "refused",
            "refused",
            {"message": str(e)[:80]},
        )
    noncanon = json.dumps(fx).encode()  # same object, uncanonical rendering
    try:
        _, rec = load_answer_packet_bytes(noncanon, pin)
        record(
            "C8c uncanonical rendering of the pinned object accepted WITH the change recorded",
            True,
            rec.get("canonicalization_changed_bytes") is True,
        )
    except Refusal:
        record(
            "C8c uncanonical rendering of the pinned object accepted WITH the change recorded",
            True,
            False,
        )
    _, rec = load_answer_packet_bytes(noncanon, sha256_bytes(noncanon))
    record(
        "C8d pin over an uncanonical rendering: accepted (the pin proves the object), form recorded",
        (False, False),
        (rec.get("canonical_form_ok"), rec.get("canonicalization_changed_bytes")),
    )
    bad_cm = copy.deepcopy(fx)
    bad_cm["convention_map"]["representation_evaluation"][EVALUATION_FIELD] = "g |-> rho(g^-1)"
    out = _run_fixture(bad_cm, raw_output, construction_packet)
    record(
        "C8e representation_evaluation.convention not byte-equal refused BY THAT CHECK",
        ("structural failure", True),
        (out["category"], "byte-equal" in out.get("refusal", "")),
    )
    for item, field, text in (
        ("bridge", "statement", "|\\tau_\\rho|^2 with no field norm to a subfield (section 5.4)"),
        ("orientation_anchor_rule", "rule", "anchor: R7, compare x and 1/x (section 5.4)"),
        ("representation_evaluation", "source", "a paraphrased source note"),
    ):
        para = copy.deepcopy(fx)
        para["convention_map"][item][field] = text
        out = _run_fixture(para, raw_output, construction_packet)
        record(
            f"C8f {item}.{field} prose is recorded, never gated (a paraphrase passes)",
            ("reproduced", para["convention_map"][item]),
            (out["category"], (out.get("convention_map") or {}).get(item)),
        )
    # C8e': each convention object is held to its EXACT inner key set (author, #451): an
    # extra field beside the verbatim one, and a missing field, each refuse.
    for item, mut in (
        ("representation_evaluation", lambda o: o.update(note="g |-> rho(g^-1)")),
        ("bridge", lambda o: o.pop("exclusivity")),
        ("orientation_anchor_rule", lambda o: o.update(extra=1)),
        ("basing_reference", lambda o: o.pop("source")),
    ):
        two = copy.deepcopy(fx)
        mut(two["convention_map"][item])
        out = _run_fixture(two, raw_output, construction_packet)
        record(
            f"C8e' {item} with a wrong inner key set refused BY THE KEY-SET CHECK",
            ("structural failure", True),
            (out["category"], f"convention_map.{item}: key set" in out.get("refusal", "")),
        )
    # An IDENTICAL duplicate: without the hook the object would still equal the pinned one
    # and be accepted through the fallback, so this control is diagnostic of the hook.
    dup = canonical_bytes(fx)[:-2] + b',\n  "format_version": "m8_8-answer-FIXTURE-2"\n}\n'
    try:
        load_answer_packet_bytes(dup, sha256_bytes(canonical_bytes(fx)))
        record("C8g duplicate key in delivered bytes refused", "refused", "accepted")
    except Refusal:
        record("C8g duplicate key in delivered bytes refused", "refused", "refused")
    greek = copy.deepcopy(fx)
    greek["convention_map"]["bridge"] = {"statement": "T²_target(ρ) := |τ_ρ|², section 5.4"}
    utf8 = (json.dumps(greek, sort_keys=True, indent=2, ensure_ascii=False) + "\n").encode("utf-8")
    _, rec = load_answer_packet_bytes(utf8, sha256_bytes(utf8))
    record(
        "C8h pin-matching bytes in another canonicalizer's rendering accepted, form recorded",
        (False, False),
        (rec.get("canonical_form_ok"), rec.get("canonicalization_changed_bytes")),
    )
    # C8i: a convention item delivered as a bare string (the #452 reading) is refused: the
    # author-confirmed items are objects, and a silent string fallback would hide a shape
    # mismatch at the reveal.
    strcm = copy.deepcopy(fx)
    strcm["convention_map"]["basing_reference"] = "construction packet `basing`"
    out = _run_fixture(strcm, raw_output, construction_packet)
    record(
        "C8i convention item as a bare string refused",
        ("structural failure", True),
        (out["category"], "non-empty object" in out.get("refusal", "")),
    )
    # C8j: the verbatim field is reported by name, and the evaluation text is the
    # construction packet's declaration, so the record shows WHICH field carried it.
    record(
        "C8j verbatim evaluation field recorded by name",
        [EVALUATION_FIELD],
        (base.get("convention_map") or {}).get("evaluation_verbatim_fields"),
    )

    # C9: the selector class on a non-R7 row is refused, even with a correct table.
    sel9 = copy.deepcopy(fx)
    r7 = sel_by_class(sel9, "free_orientation_selector")
    r4 = [r for r in sel9["rows"] if r["signature"]["dimension"] == 4][0]
    r7["evidentiary_class"], r4["evidentiary_class"] = "free", "free_orientation_selector"
    out = _run_fixture(sel9, raw_output, construction_packet)
    record("C9 selector class on a non-R7 row refused", "structural failure", out["category"])

    # C10: two packet rows sent to one committed position is refused (injectivity).
    inj = copy.deepcopy(dfx)
    ents = inj["indexing_map"]["entries"]
    pos_r1 = inj["rows"].index(by_label(inj, "R1"))
    pos_r2 = inj["rows"].index(by_label(inj, "R2"))
    e1 = [e for e in ents if e["source_position"] == pos_r1][0]
    e2 = [e for e in ents if e["source_position"] == pos_r2][0]
    e2["destination_position"] = e1["destination_position"]
    # run as an arm: on the live path the self-contradiction check (C3e) fires first
    out = _run_fixture(dfx, raw_output, construction_packet, map_override=inj["indexing_map"])
    record(
        "C10 non-injective map refused BY THE INJECTIVITY CHECK",
        ("structural failure", True),
        (out["category"], "not injective" in out.get("refusal", "")),
    )

    # C11: structurally trivial identities are refused (R0 factor; zero exponent; wrong kinds).
    triv = copy.deepcopy(fx)
    r0sig = sel_by_class(triv, "declared_convention")["signature"]
    triv["identities"] = [
        {
            "position": i,
            "factors": [{"signature": r0sig, "exponent": 0}] * 2,
            "expected": [1, 0, 1],
        }
        for i in range(4)
    ]
    out = _run_fixture(triv, raw_output, construction_packet)
    record("C11a trivial identities (R0^0) refused", "structural failure", out["category"])
    kinds = copy.deepcopy(fx)
    kinds["identities"][0] = copy.deepcopy(kinds["identities"][2])
    kinds["identities"][0]["position"] = 0
    out = _run_fixture(kinds, raw_output, construction_packet)
    record("C11b three products and one ratio refused", "structural failure", out["category"])
    dpos = copy.deepcopy(fx)
    dpos["identities"][1]["position"] = 0
    out = _run_fixture(dpos, raw_output, construction_packet)
    record(
        "C11c two identities at one position refused BY THAT CHECK",
        ("structural failure", True),
        (out["category"], "positions not distinct" in out.get("refusal", "")),
    )
    spos = copy.deepcopy(fx)
    spos["identities"][1]["position"] = "1"
    out = _run_fixture(spos, raw_output, construction_packet)
    record(
        "C11d identity position as a string refused",
        ("structural failure", True),
        (out["category"], "non-negative int" in out.get("refusal", "")),
    )

    # C12: a malformed but otherwise pinned packet yields a recorded structural failure,
    # never a traceback.
    mal = copy.deepcopy(fx)
    del mal["rows"][3]["value"]
    out = _run_fixture(mal, raw_output, construction_packet)
    record(
        "C12 malformed row (missing T_squared) is a recorded refusal",
        ("structural failure", True),
        (out["category"], "refusal" in out),
    )

    # C13: EXACTNESS.  A comparison that were approximate, sign-blind or conjugation-blind
    # would pass C1; these three mutations of one free cell each require exactly one
    # mismatch: negation, Galois conjugation, and a near-miss rational.
    # The target is a free row with a nonzero phi-part, so conjugation is a real mutation.
    target = [
        r
        for r in fx["rows"]
        if r["evidentiary_class"] == "free" and QPhi.from_triple(r["value"]).y != 0
    ][0]
    v = QPhi.from_triple(target["value"])
    near = QPhi(v.x + Fraction(1, 10**20), v.y)  # beyond any float tolerance
    for name, mutant in (
        ("negation", QPhi(-v.x, -v.y)),
        ("Galois conjugation", v.conjugate()),
        ("near-miss rational (1e-20 off)", near),
    ):
        if mutant == v:
            record(f"C13 exactness: {name} is a real mutation", True, False)
            continue
        ex = copy.deepcopy(fx)
        by_label(ex, target["label"])["value"] = mutant.to_triple()
        out = _run_fixture(ex, raw_output, construction_packet)
        record(
            f"C13 exactness: {name} of one cell is exactly one mismatch",
            ("partial disagreement", [target["label"]]),
            (out["category"], sorted(x["label"] for x in (out.get("row_mismatches") or []))),
        )

    # C14: one control per identity-structure rule, each violating exactly that rule on an
    # otherwise valid fixture (C11a hit several at once and so proved none individually).
    # Each rule is identified by its refusal MESSAGE, not only by the category, so a rule
    # masked by a sibling refusal cannot pass its own control.
    C14_MSG = {
        "non-Galois ratio": "not a Galois pair",
        "overlapping products": "not disjoint",
        "R0 factor": "declared-convention row is not a factor",
        "zero exponent": "nonzero int",
        "repeated factor": "repeated factor",
        "same pair in both ratios": "same Galois pair",
        "conjugate key on a factor (true)": "factor: key set",
        "conjugate key on a factor (string)": "factor: key set",
        "exponent beyond the bound": "exceeds the bound",
    }

    def with_ids(mutator, msg):
        z = copy.deepcopy(fx)
        mutator(z["identities"])
        out = _run_fixture(z, raw_output, construction_packet)
        return out["category"], msg in out.get("refusal", "")

    dim3 = [r["signature"] for r in fx["rows"] if r["signature"]["dimension"] == 3]
    r0sig = sel_by_class(fx, "declared_convention")["signature"]

    def non_galois(ids):  # dim-2 over dim-3: same exponents, not a pair
        ids[0]["factors"][1]["signature"] = dim3[0]

    def overlap(ids):  # product b re-uses a row of product a
        ids[3]["factors"][0]["signature"] = ids[2]["factors"][0]["signature"]

    def r0_factor(ids):
        ids[2]["factors"].append({"signature": r0sig, "exponent": 1})

    def zero_exp(ids):
        ids[2]["factors"][0]["exponent"] = 0

    def repeated(ids):
        ids[2]["factors"].append(copy.deepcopy(ids[2]["factors"][0]))

    def dup_pair(ids):  # both ratio slots on the same Galois pair
        ids[1]["factors"] = copy.deepcopy(ids[0]["factors"])

    def conj_true(ids):  # outside the frozen factor schema: refused, never applied
        ids[2]["factors"][0]["conjugate"] = True

    def conj_string(ids):  # the string "false" is refused the same way
        ids[2]["factors"][0]["conjugate"] = "false"

    def exp_nine(ids):
        ids[2]["factors"][0]["exponent"] = 9

    for name, mut in (
        ("non-Galois ratio", non_galois),
        ("overlapping products", overlap),
        ("R0 factor", r0_factor),
        ("zero exponent", zero_exp),
        ("repeated factor", repeated),
        ("same pair in both ratios", dup_pair),
        ("conjugate key on a factor (true)", conj_true),
        ("conjugate key on a factor (string)", conj_string),
        ("exponent beyond the bound", exp_nine),
    ):
        record(
            f"C14 identity rule: {name} refused BY ITS OWN CHECK",
            ("structural failure", True),
            with_ids(mut, C14_MSG[name]),
        )
    # A `conjugate: true` factor whose expected value was recomputed WITH the conjugation
    # is still refused at the schema layer: the key is outside the frozen factor schema and
    # would change the computed identity, so it is never honored silently (author, #451).
    cj = copy.deepcopy(fx)
    f0 = cj["identities"][2]["factors"][0]
    f0["conjugate"] = True
    key0 = signature_key(f0["signature"])
    tbl = {signature_key(r["signature"]): QPhi.from_triple(r["value"]) for r in cj["rows"]}
    acc = ONE
    for f in cj["identities"][2]["factors"]:
        val = tbl[signature_key(f["signature"])]
        acc = acc * (val.conjugate() if signature_key(f["signature"]) == key0 else val)
    cj["identities"][2]["expected_value"] = acc.to_triple()
    out = _run_fixture(cj, raw_output, construction_packet)
    record(
        "C14 a conjugate factor is refused at the schema layer even with a consistent expected value",
        ("structural failure", True),
        (out["category"], "factor: key set" in out.get("refusal", "")),
    )

    # C15: the Phase B checker, each rule in isolation, on a synthetic record shaped like
    # the author's.  (The live run consumes the committed record; this proves the checker.)
    man = MANIFEST_PATH.read_text()
    gids = sorted(
        line.split("|")[1].strip().split()[0]
        for line in man.splitlines()
        if line.startswith("| G-")
    )
    good_pb = {
        "schema_version": PHASE_B_SCHEMA_PREFIX + "0",
        "manifest_sha256": PIN_MANIFEST,
        "phase_a_hashes_verified_pre": True,
        "phase_a_hashes_verified_post": True,
        "manifest_parsed_at_runtime": True,
        "parser_self_test_passed": True,
        "all_mutations_reddened": True,
        "registry_coverage": {
            "executed_gate_ids": gids,
            "count": len(gids),
            "pre_execution_set_equality": True,
            "post_execution_set_equality": True,
        },
        "results": [
            {
                "gate_id": g,
                "red_outcome": True,
                "baseline_result": "PASS (x)",
                "mutated_result": "FAIL (y)",
                "object_mutated": "o",
                "implemented_mutation": "i",
                "declared_mutation": "d",
            }
            for g in gids
        ],
    }
    record("C15 phase B: well-formed record passes", [], check_phase_b(good_pb, man))

    def pb(mutator):
        z = copy.deepcopy(good_pb)
        mutator(z)
        return len(check_phase_b(z, man)) > 0

    for name, mut in (
        ("missing gate", lambda z: z["results"].pop()),
        (
            "duplicate gate",
            lambda z: (
                z["results"].append(copy.deepcopy(z["results"][0])),
                z["registry_coverage"].update(count=len(z["results"])),
            ),
        ),
        ("non-red gate", lambda z: z["results"][0].update(red_outcome=False)),
        ("baseline not PASS", lambda z: z["results"][0].update(baseline_result="FAIL")),
        ("mutated not FAIL", lambda z: z["results"][0].update(mutated_result="PASS")),
        ("manifest pin mismatch", lambda z: z.update(manifest_sha256="0" * 64)),
        ("post-hash flag false", lambda z: z.update(phase_a_hashes_verified_post=False)),
        ("coverage contradicts results", lambda z: z["registry_coverage"].update(count=1)),
        ("empty implemented_mutation", lambda z: z["results"][0].update(implemented_mutation="")),
    ):
        record(f"C15 phase B: {name} flagged", True, pb(mut))

    # C16: remaining refusals, one each.
    tid = copy.deepcopy(fx)
    tid["target_id"] = ""
    try:
        load_answer_packet_bytes(canonical_bytes(tid), sha256_bytes(canonical_bytes(tid)))
        record("C16a empty target_id refused", "refused", "accepted")
    except Refusal:
        record("C16a empty target_id refused", "refused", "refused")
    own99 = copy.deepcopy(fx)
    own99["rows"][2]["signature"]["dimension"] = 99
    out = _run_fixture(own99, raw_output, construction_packet)
    record(
        "C16b own signature absent from the committed output refused BY THAT CHECK",
        ("structural failure", True),
        (out["category"], "own signature" in out.get("refusal", "")),
    )
    # C16c: T^2(R0) = 1 is a declared convention (section 9 item 2); any other value refuses.
    r0v = copy.deepcopy(fx)
    sel_by_class(r0v, "declared_convention")["value"] = [7, 0, 1]
    out = _run_fixture(r0v, raw_output, construction_packet)
    record(
        "C16c declared_convention row carrying a value other than 1 refused",
        ("structural failure", True),
        (out["category"], "is not 1" in out.get("refusal", "")),
    )
    # C16e: an extra key on a row (a stale `row_signature` beside `signature`) refuses.
    extra = copy.deepcopy(fx)
    extra["rows"][2]["row_signature"] = {"dim": 1}
    out = _run_fixture(extra, raw_output, construction_packet)
    record(
        "C16e extra key on a packet row refused BY THE KEY-SET CHECK",
        ("structural failure", True),
        (out["category"], "packet row: key set" in out.get("refusal", "")),
    )
    # C16d: the #452 spelling of a signature (`row_signature` with `dim`/`chi_*`) is
    # refused by the field check, so an old-shape packet cannot pass through by accident.
    old = copy.deepcopy(fx)
    sig = old["rows"][2].pop("signature")
    old["rows"][2]["row_signature"] = {
        "dim": sig["dimension"],
        "chi_s": sig["s"],
        "chi_t": sig["t"],
        "chi_st": sig["st"],
    }
    out = _run_fixture(old, raw_output, construction_packet)
    record(
        "C16d the superseded #452 signature spelling refused",
        ("structural failure", True),
        (out["category"], "packet row: key set" in out.get("refusal", "")),
    )

    # C17: the terminal summary survives a SUCCESS-shaped result.  The honest fixture is
    # the only self-test arm that reaches a success category, and its identities list ends
    # with the `_sector_coverage` record entry, which carries no `recomputed`/`expected`;
    # the 2026-08-22 step-6 diagnostic crashed exactly there, after the record was written
    # and before the exit code.  The printer is pure, so it is called on the real result.
    honest = build_fixture(raw_output, construction_packet)
    out = _run_fixture(honest, raw_output, construction_packet)
    try:
        text = summary_lines(out)
        printed = (out["category"] in SUCCESS, "_sector_coverage" in text, "identity" in text)
    except Exception as e:  # noqa: BLE001
        printed = (f"{type(e).__name__}: {e}", False, False)
    record(
        "C17 summary printer handles the coverage record on a success result",
        (True, True, True),
        printed,
    )
    return results


# ----------------------------------------------------------------------------------------
def summary_lines(result: dict) -> str:
    """The terminal summary of an adjudication result.  Pure: the record is already written
    before this runs, so a failure here could only lose the exit code, which is why it is
    exercised by the self-test (C17) on a success-shaped result that carries the
    `_sector_coverage` record entry the identities list ends with."""
    lines = [f"category     {result['category'].upper()}", f"orientation  {result.get('orientation')}"]
    for m in result.get("row_mismatches") or []:
        lines.append(
            f"  mismatch {m['label']} {m['signature']}: reference {m['reference']} "
            f"vs selected committed {m['selected_committed']}"
        )
    for i in result.get("identities") or []:
        if i.get("kind") == "record":
            lines.append(
                f"  {i['slot']}: {i['rows_in_products']}/{i['nontrivial_rows']} nontrivial rows "
                f"in the two sector products, uncovered {i['uncovered']}"
            )
            continue
        lines.append(
            f"  identity {i['slot']}: {'equal' if i['equal'] else 'UNEQUAL'} "
            f"recomputed {i['recomputed']} expected {i['expected']}"
        )
    return "\n".join(lines)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--packet", type=Path, help="decrypted canonical answer packet (step 6)")
    ap.add_argument("--phase-b", type=Path, default=PHASE_B_DEFAULT)
    ap.add_argument("--json", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()

    try:
        raw_output = load_pinned_json(RAW_OUTPUT_PATH, PIN_RAW_OUTPUT, "raw output")
        construction = load_pinned_json(
            CONSTRUCTION_PACKET_PATH, PIN_CONSTRUCTION_PACKET, "construction packet"
        )
        load_pinned_json(GROUP_PACKET_PATH, PIN_GROUP_PACKET, "group packet")
        manifest_bytes = MANIFEST_PATH.read_bytes()
        if sha256_bytes(manifest_bytes) != PIN_MANIFEST:
            raise Refusal("manifest: SHA-256 != Addendum 1 pin")
    except Refusal as e:
        print(f"REFUSED  {e}")
        return 1

    if args.self_test:
        results = self_test(raw_output, construction)
        all_ok = all(r["ok"] for r in results)
        for r in results:
            print(
                f"  {'OK  ' if r['ok'] else 'FAIL'}  {r['control']}: expected "
                f"{r['expected']!r}, observed {r['observed']!r}"
            )
        print(
            f"self-test    {'ALL CONTROLS BEHAVED' if all_ok else 'A CONTROL DID NOT BEHAVE'}"
            f" ({sum(r['ok'] for r in results)}/{len(results)})"
        )
        return 0 if all_ok else 1

    if args.packet is None:
        print("nothing to do: pass --self-test, or --packet PATH at step 6")
        return 1

    # Live run.  Every precondition is a refusal, recorded, before any comparison.
    record = {
        "what": "M8.8 section 8 steps 6-8 adjudication",
        "pins": {
            "answer_packet": PIN_ANSWER_PACKET,
            "raw_output": PIN_RAW_OUTPUT,
            "manifest": PIN_MANIFEST,
            "group_packet": PIN_GROUP_PACKET,
            "construction_packet": PIN_CONSTRUCTION_PACKET,
        },
    }
    try:
        structural = check_raw_output_structure(raw_output)
        if not args.phase_b.exists():
            raise Refusal(
                f"phase B record absent at {args.phase_b}; Addendum 1 requires it " "before step 6"
            )
        structural += check_phase_b(
            json.loads(args.phase_b.read_text()), manifest_bytes.decode("utf-8")
        )
        if structural:
            record["structural_problems"] = structural
            raise Refusal("; ".join(structural))
        packet_bytes = args.packet.read_bytes()
        record["delivered_bytes"] = {
            "sha256": sha256_bytes(packet_bytes),
            "length": len(packet_bytes),
        }
        packet, load_rec = load_answer_packet_bytes(packet_bytes)
        record["load"] = load_rec
        record["target_id"] = packet["target_id"]
        record["format_version"] = packet["format_version"]
        result = adjudicate(packet, raw_output, construction)
    except Refusal as e:
        record["category"] = "structural failure"
        record["refusal"] = str(e)
        args.json.write_text(json.dumps(record, indent=2) + "\n")
        print(f"STRUCTURAL FAILURE  {e}")
        print(f"written      {args.json}")
        return 1
    except Exception as e:  # malformed-but-pinned input: a record is still owed
        record["category"] = "structural failure"
        record["exception"] = f"{type(e).__name__}: {e}"
        args.json.write_text(json.dumps(record, indent=2) + "\n")
        print(f"STRUCTURAL FAILURE  uncaught {type(e).__name__}: {e}")
        print(f"written      {args.json}")
        return 1

    record.update(result)
    args.json.write_text(json.dumps(record, indent=2, default=str) + "\n")
    print(summary_lines(result))
    print(f"written      {args.json}")
    if result["category"] in SUCCESS:
        return 0
    if result["category"] in FINDINGS:
        return 2
    return 1


if __name__ == "__main__":
    sys.exit(main())
