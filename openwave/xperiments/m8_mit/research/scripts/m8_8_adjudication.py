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
           the Addendum 1 pin; the indexing map is APPLIED to pair packet rows with
           committed rows; the convention map is LOADED AND VALIDATED (four items, the
           `basing_evaluation` item byte-equal to the construction packet's
           `basing.evaluation`) without applying any orientation.
  step 8   the section 5.4 selection at R7: exactly one of x = r, x^-1 = r must hold
           (both: INVALID ANCHOR; neither: DISAGREEMENT; no row comparison in either case).
           The selected identity or global inverse is applied to the COMMITTED rows, the
           seven free forms and R7 compared exactly in Q(phi), and the four identities
           recomputed from the SELECTED rows under the packet's slot definitions and
           compared POSITION-WISE (section 5.3).  Outcome recorded by section 8 category.

THE PACKET SHAPE THIS HARNESS READS.  Section 4.4 fixes the key set and describes each
value; it does not print a JSON skeleton.  The key-level shape below is this harness's
reading of section 4.4, stated here so it can be confirmed by the author BEFORE reveal
(keys and types only, no values) and so any post-reveal adapter change is visible as a
diff against a committed file rather than a silent accommodation.

  {
    "format_version": str,
    "target_id": str,
    "adjudicates": {"group_packet_sha256": hex, "construction_packet_sha256": hex},
    "rows": [ {"label": str,
               "row_signature": {"dim": int, "chi_s": [a,b,c], "chi_t": [a,b,c],
                                 "chi_st": [a,b,c]},
               "class": "declared_convention" | "free" | "free_orientation_selector",
               "T_squared": [a,b,c]}            x 9 ],
    "identities": [ {"slot": str,
                     "factors": [{"row_signature": {...}, "exponent": int,
                                  "conjugate": bool (optional, default false)}, ...],
                     "expected": [a,b,c]}        x 4 ],
    "indexing_map": {"mode": "signature_identity"}
                  | {"mode": "signature_table",
                     "entries": [{"label": str, "row_signature": {...}}, ...]},
    "convention_map": {"bridge": str, "anchor_rule": str,
                       "basing_reference": str, "basing_evaluation": str}
  }

  Every Q(phi) value is the normalized triple (a, b, c) for (a + b*phi)/c, c > 0,
  gcd(a, b, c) = 1, as section 4.2 fixes and the raw output already uses.  In
  `signature_identity` mode a packet row is paired with the committed row carrying the
  same signature.  In `signature_table` mode the packet's own row signatures are NOT
  used for pairing: each packet row's label is looked up in `entries` and the signature
  recorded THERE selects the committed row.  The synthetic fixture of the self-test uses
  table mode with a derangement, which is how the map-processing path is proved operative
  (section 4.4: a live map may be too simple to discriminate applying from ignoring it).

THE CONTROLS, section 4.4, run by --self-test and REQUIRED green before any live run.
An independent adversarial audit of the first version found C1, C2 and C3 tautological
(a comparison skipping six of eight rows, or a set-wise identity comparison, passed them)
and found the anchor label-driven rather than R7-driven; the controls below are the
repaired set, each phrased so the named defect would fail it:
  C1  EACH of the eight nontrivial reference cells mutated in turn, downstream of the
      completed hash check: exactly that one label mismatches (R7: `disagreement`);
  C2  EACH committed raw cell mutated the same way, from the other side;
  C3  synthetic NONIDENTITY indexing fixture (a derangement over the seven FREE labels,
      R0 and R7 fixed): the supplied map is GREEN; ignoring the map and a preregistered
      wrong map each fail AT THE COMPARISON with seven mismatches, not by an earlier
      refusal; a map moving R0 onto an acyclic row refuses on class (C3d);
  C4  the global inverse of the committed table resolves to `convention difference`; the
      two sector-product expected values swapped redden BOTH slots on the identity and the
      inverted fixture, which is the position-wise rule of section 5.3 exercised on the
      harness rather than on the fixture arithmetic;
  C5  a self-inverse R7 reference resolves to `invalid anchor` with no row comparison;
  C6  an R7 reference equal to neither x nor x^-1 resolves to `disagreement`;
  C7  a tampered identity `expected` reddens the identity layer alone;
  C8  tampered `adjudicates`, content-tampered bytes, a duplicate key in the delivered
      bytes and a non-verbatim `basing_evaluation` each REFUSE; an uncanonical rendering
      of the pinned object is accepted with the change recorded, and pin-matching bytes in
      another canonicalizer's rendering are accepted with the form recorded; the bridge
      and anchor prose are RECORDED, never gated (a second audit showed a substring check
      refusing section 5.4's own wording);
  C9  the selector class on any row but the dim-5 irrep refuses, even with a correct table;
  C10 a non-injective map refuses;
  C11 structurally trivial identities refuse;
  C12 a malformed but pinned packet yields a RECORDED structural failure, not a traceback;
  C13 EXACTNESS: negation, Galois conjugation and a 1e-9 near-miss of one cell are each
      exactly one mismatch, so an approximate, sign-blind or conjugation-blind equality
      would fail here;
  C14 one control per identity rule (non-Galois ratio, overlapping products, R0 factor,
      zero exponent, repeated factor, same pair in both ratios, conjugate inside a ratio,
      conjugate as a string, exponent beyond the bound), plus a genuine conjugate factor
      honored;
  C15 the Phase B checker, each rule in isolation on a synthetic record;
  C16 empty `target_id` and an own signature absent from the committed output refuse.
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
CONVENTION_KEYS = {"bridge", "anchor_rule", "basing_reference", "basing_evaluation"}
CLASSES = {"declared_convention", "free", "free_orientation_selector"}
SIG_FIELDS = ("dim", "chi_s", "chi_t", "chi_st")
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


def signature_key(sig: dict) -> tuple:
    if not isinstance(sig, dict) or set(sig) != set(SIG_FIELDS):
        raise Refusal(
            f"row_signature fields {sorted(sig) if isinstance(sig, dict) else sig!r}"
            f" != {sorted(SIG_FIELDS)}"
        )
    dim = sig["dim"]
    if not isinstance(dim, int) or isinstance(dim, bool) or dim <= 0:
        raise Refusal(f"row_signature.dim not a positive int: {dim!r}")
    return (dim,) + tuple(tuple(QPhi.from_triple(sig[f]).to_triple()) for f in SIG_FIELDS[1:])


def load_raw_rows(raw_output: dict) -> dict:
    """Committed rows keyed by signature -> {'acyclic': bool, 'value': QPhi|None}."""
    rows = {}
    for r in raw_output["rows"]:
        key = signature_key(r["row_signature"])
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
        if not isinstance(v, str) or not v.strip():
            raise Refusal(f"convention_map.{k}: empty or not a string")
    declared = construction_packet["basing"]["evaluation"]
    if cm["basing_evaluation"] != declared:
        raise Refusal(
            "convention_map.basing_evaluation is not the construction packet's "
            "basing.evaluation verbatim"
        )
    # The other three items are prose the protocol does not pin verbatim.  A free-text
    # check cannot tell a wrong bridge from a paraphrase (an earlier heuristic here refused
    # section 5.4's own wording), so they are RECORDED for the adjudication record and
    # judged by the adjudicator there, never gated by substring.
    return {
        "validated": True,
        "orientation_applied": False,
        "basing_evaluation": cm["basing_evaluation"],
        "bridge": cm["bridge"],
        "anchor_rule": cm["anchor_rule"],
        "basing_reference": cm["basing_reference"],
    }


def apply_indexing_map(packet: dict, raw_rows: dict, map_override: dict | None = None) -> dict:
    """Pair each packet row with ONE committed row.  Returns label -> paired record.

    `map_override` is used ONLY by the self-test arms (identity arm, wrong-map arm); the
    live run always consumes the packet's own map.
    """
    imap = map_override if map_override is not None else packet["indexing_map"]
    if not isinstance(imap, dict) or "mode" not in imap:
        raise Refusal("indexing_map: missing mode")
    rows = packet["rows"]
    if not isinstance(rows, list) or len(rows) != N_ROWS:
        raise Refusal(
            f"packet rows: {len(rows) if isinstance(rows, list) else rows!r} != {N_ROWS}"
        )
    labels = [need(r, "label", "packet row") for r in rows]
    if len(set(labels)) != N_ROWS or any(not isinstance(lab, str) for lab in labels):
        raise Refusal("packet rows: labels not nine distinct strings")
    # The packet's own signatures: pairwise distinct and each present among the committed
    # signatures, in every mode, so a packet cannot carry rows that describe nothing.
    own = {r["label"]: signature_key(need(r, "row_signature", f"row {r['label']}")) for r in rows}
    if len(set(own.values())) != N_ROWS:
        raise Refusal("packet rows: row signatures not pairwise distinct")
    for lab, key in own.items():
        if key not in raw_rows:
            raise Refusal(f"row {lab}: own signature {key} absent from the committed output")
    if imap["mode"] == "signature_identity":
        lookup = dict(own)
    elif imap["mode"] == "signature_table":
        entries = imap.get("entries")
        if not isinstance(entries, list) or len(entries) != N_ROWS:
            raise Refusal("indexing_map.entries: not nine entries")
        lookup = {}
        for e in entries:
            lab = need(e, "label", "indexing_map entry")
            if lab in lookup:
                raise Refusal(f"indexing_map.entries: duplicate label {lab}")
            lookup[lab] = signature_key(need(e, "row_signature", f"indexing_map entry {lab}"))
        if set(lookup) != set(labels):
            raise Refusal("indexing_map.entries: label set != packet row labels")
    else:
        raise Refusal(f"indexing_map.mode {imap['mode']!r} unknown")

    paired, used = {}, set()
    for r in rows:
        key = lookup[r["label"]]
        if key not in raw_rows:
            raise Refusal(
                f"indexing map sends {r['label']} to a signature absent from the "
                f"committed output: {key}"
            )
        if key in used:
            raise Refusal(f"indexing map is not injective at {key}")
        used.add(key)
        cls = need(r, "class", f"row {r['label']}")
        if cls not in CLASSES:
            raise Refusal(f"row {r['label']}: class {cls!r} not in {sorted(CLASSES)}")
        paired[r["label"]] = {
            "signature": key,
            "own_signature_differs_from_map": own[r["label"]] != key,
            "class": cls,
            "reference": QPhi.from_triple(need(r, "T_squared", f"row {r['label']}")),
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
        slot = need(slot_entry, "slot", "identity")
        factors = need(slot_entry, "factors", f"identity {slot}")
        if not isinstance(slot, str) or not isinstance(factors, list) or len(factors) < 2:
            raise Refusal(f"identities: slot {slot!r} malformed or with fewer than two factors")
        acc, keys, exps, conjs = ONE, [], [], []
        for f in factors:
            key = signature_key(need(f, "row_signature", f"identity {slot} factor"))
            if key not in sel:
                raise Refusal(f"identity {slot}: factor signature {key} not among selected rows")
            if key in trivial_keys:
                raise Refusal(f"identity {slot}: the declared-convention row is not a factor")
            e = need(f, "exponent", f"identity {slot} factor")
            if not isinstance(e, int) or isinstance(e, bool) or e == 0:
                raise Refusal(f"identity {slot}: exponent must be a nonzero int")
            v = sel[key]
            conj = f.get("conjugate", False)
            if not isinstance(conj, bool):
                raise Refusal(f"identity {slot}: conjugate flag must be a JSON boolean")
            if conj:
                v = v.conjugate()
            acc = acc * (v**e)
            keys.append(key)
            exps.append(e)
            conjs.append(conj)
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
            if any(conjs):
                # conj(T(rho')) = T(rho) on a Galois pair, so a conjugate flag inside a
                # ratio makes it identically 1 and carries no information.
                raise Refusal(f"identity {slot}: conjugate flag inside a Galois ratio")
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
        expected = QPhi.from_triple(need(slot_entry, "expected", f"identity {slot}"))
        out.append(
            {
                "slot": slot,
                "kind": kinds[-1],
                "recomputed": acc.to_triple(),
                "expected": expected.to_triple(),
                "equal": acc == expected,
            }
        )
    if len({i["slot"] for i in out}) != N_IDENTITIES:
        raise Refusal("identities: slots not distinct")
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
            "class": p["class"],
            "own_signature_differs_from_map": p["own_signature_differs_from_map"],
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
    paired = apply_indexing_map(packet, raw_rows, map_override)

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
    return {"dim": key[0], "chi_s": list(key[1]), "chi_t": list(key[2]), "chi_st": list(key[3])}


def build_fixture(
    raw_output: dict, construction_packet: dict, *, invert: bool = False, derange: bool = False
) -> dict:
    """A synthetic answer packet whose reference table IS the committed table (or its
    global inverse).  It contains no reference value and proves only harness behavior.

    Labels follow M8.5-A's label-to-dimension map for readability; R0 is the trivial row,
    R7 the dim-5 row.  With derange=True the packet's rows are labeled by a cyclic shift
    and an explicit signature table restores the true pairing, so a harness that ignores
    the table pairs every row wrongly.
    """
    raw_rows = load_raw_rows(raw_output)
    keys = sorted(raw_rows)  # deterministic order
    trivial = [k for k in keys if not raw_rows[k]["acyclic"]][0]
    acyclic = [k for k in keys if raw_rows[k]["acyclic"]]
    r7 = [k for k in acyclic if k[0] == 5][0]
    labels = ["R0"] + [f"R{i}" for i in range(1, 9)]
    order = (
        [trivial]
        + [k for k in acyclic if k != r7][:6]
        + [r7]
        + [k for k in acyclic if k != r7][6:]
    )
    # order: R0 trivial, R1..R6 six free rows, R7 the selector, R8 the seventh free row
    assert len(order) == N_ROWS

    def value_of(k):
        if k == trivial:
            return ONE
        v = raw_rows[k]["value"]
        return v.inverse() if invert else v

    true_pairs = dict(zip(labels, order))
    # Derangement over the seven FREE labels only (R0 and R7 keep their signatures), so a
    # harness that ignores the table reaches the comparison and fails THERE, on values,
    # rather than being stopped earlier by a class refusal.
    free_labels = [lab for lab in labels if lab not in ("R0", "R7")]
    shifted = dict(zip(free_labels, free_labels[1:] + free_labels[:1]))
    rows, entries = [], []
    for lab in labels:
        true_key = true_pairs[lab]
        shown_key = true_pairs[shifted[lab]] if (derange and lab in shifted) else true_key
        cls = (
            "declared_convention"
            if lab == "R0"
            else "free_orientation_selector" if lab == "R7" else "free"
        )
        rows.append(
            {
                "label": lab,
                "row_signature": _sig_dict(shown_key),
                "class": cls,
                "T_squared": value_of(true_key).to_triple(),
            }
        )
        entries.append({"label": lab, "row_signature": _sig_dict(true_key)})

    # Four identity slots: two Galois-pair ratios (dims 2 and 3 pairs), two sector
    # products over disjoint halves.  Slot membership here is synthetic; the live packet
    # declares its own.  Expected values are recomputed from the same table, so the
    # fixture is self-consistent by construction.
    def ident(slot, spec):
        acc = ONE
        for k, e in spec:
            acc = acc * (value_of(k) ** e)
        return {
            "slot": slot,
            "factors": [{"row_signature": _sig_dict(k), "exponent": e} for k, e in spec],
            "expected": acc.to_triple(),
        }

    dim2 = [k for k in acyclic if k[0] == 2]
    dim3 = [k for k in acyclic if k[0] == 3]
    half_a = [order[1], order[2], order[3]]
    half_b = [order[4], order[5], order[8]]
    identities = [
        ident("galois_ratio_dim2", [(dim2[0], 1), (dim2[1], -1)]),
        ident("galois_ratio_dim3", [(dim3[0], 1), (dim3[1], -1)]),
        ident("sector_product_a", [(k, 1) for k in half_a]),
        ident("sector_product_b", [(k, 1) for k in half_b]),
    ]
    imap = (
        {"mode": "signature_table", "entries": entries}
        if derange
        else {"mode": "signature_identity"}
    )
    return {
        "format_version": "m8_8-answer-FIXTURE-1",
        "target_id": "SELF-TEST FIXTURE, NOT THE CANONICAL PACKET",
        "adjudicates": {
            "group_packet_sha256": PIN_GROUP_PACKET,
            "construction_packet_sha256": PIN_CONSTRUCTION_PACKET,
        },
        "rows": rows,
        "identities": identities,
        "indexing_map": imap,
        "convention_map": {
            "bridge": "T^2_target(rho) := |tau_rho|^2, section 5.4; involution T^2 <-> 1/T^2",
            "anchor_rule": "select at R7 between the committed table and its global inverse",
            "basing_reference": "construction packet `basing`, protocol section 4.2",
            "basing_evaluation": construction_packet["basing"]["evaluation"],
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

    # Baseline: the committed table against itself must be `reproduced`.
    fx = build_fixture(raw_output, construction_packet)
    base = _run_fixture(fx, raw_output, construction_packet)
    record("C0 baseline identity fixture", "reproduced", base["category"])

    # C1: EACH nontrivial reference cell mutated in turn on an in-memory copy, downstream
    # of hashing; the comparison must go red with exactly that one label mismatching, so a
    # comparison that silently skipped any row would be caught here.
    nontrivial = [r["label"] for r in fx["rows"] if r["class"] != "declared_convention"]
    for lab in nontrivial:
        m = copy.deepcopy(fx)
        row = [r for r in m["rows"] if r["label"] == lab][0]
        row["T_squared"] = (QPhi.from_triple(row["T_squared"]) * QPhi(2, 0)).to_triple()
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

    # C3: nonidentity (derangement over the seven free labels) indexing fixture, three
    # arms.  The identity arm and the wrong arm must reach the comparison and fail there.
    dfx = build_fixture(raw_output, construction_packet, derange=True)
    supplied = _run_fixture(dfx, raw_output, construction_packet)
    identity_arm = _run_fixture(
        dfx, raw_output, construction_packet, map_override={"mode": "signature_identity"}
    )
    wrong = copy.deepcopy(dfx["indexing_map"])
    free_idx = [i for i, e in enumerate(wrong["entries"]) if e["label"] not in ("R0", "R7")]
    sigs = [wrong["entries"][i]["row_signature"] for i in free_idx]
    for i, sig in zip(free_idx, sigs[2:] + sigs[:2]):
        wrong["entries"][i]["row_signature"] = sig
    wrong_arm = _run_fixture(dfx, raw_output, construction_packet, map_override=wrong)
    record("C3a derangement map, supplied", "reproduced", supplied["category"])
    record(
        "C3b derangement map, identity arm (map ignored): red AT THE COMPARISON",
        ("partial disagreement", 7),
        (identity_arm["category"], len(identity_arm.get("row_mismatches") or [])),
    )
    record(
        "C3c derangement map, preregistered wrong map: red AT THE COMPARISON",
        ("partial disagreement", 7),
        (wrong_arm["category"], len(wrong_arm.get("row_mismatches") or [])),
    )
    # C3d: a map that moves R0 onto an acyclic row is refused on class, separately.
    r0move = copy.deepcopy(dfx["indexing_map"])
    e0 = [e for e in r0move["entries"] if e["label"] == "R0"][0]
    e1 = [e for e in r0move["entries"] if e["label"] == "R1"][0]
    e0["row_signature"], e1["row_signature"] = e1["row_signature"], e0["row_signature"]
    out = _run_fixture(dfx, raw_output, construction_packet, map_override=r0move)
    record("C3d map moving R0 onto an acyclic row refused", "structural failure", out["category"])

    # C4: the global inverse resolves to convention difference; the sector-product slots
    # EXCHANGE values under inversion, and the harness compares them position-wise: swapping
    # the two expected values reddens both slots, on the identity and the inverted fixture.
    ifx = build_fixture(raw_output, construction_packet, invert=True)
    out = _run_fixture(ifx, raw_output, construction_packet)
    record("C4a inverted fixture", "convention difference", out["category"])
    for name, base_fx in (("identity", fx), ("inverted", ifx)):
        sw = copy.deepcopy(base_fx)
        ia = [i for i in sw["identities"] if i["slot"] == "sector_product_a"][0]
        ib = [i for i in sw["identities"] if i["slot"] == "sector_product_b"][0]
        ia["expected"], ib["expected"] = ib["expected"], ia["expected"]
        out = _run_fixture(sw, raw_output, construction_packet)
        record(
            f"C4b sector expected values swapped on the {name} fixture: position-wise red",
            ("partial disagreement", ["sector_product_a", "sector_product_b"], 0),
            (
                out["category"],
                sorted(i["slot"] for i in (out.get("identities") or []) if not i["equal"]),
                len(out.get("row_mismatches") or []),
            ),
        )

    # C5: self-inverse anchor -> invalid anchor.  Force R7's reference AND committed to 1.
    sfx = copy.deepcopy(fx)
    sraw = copy.deepcopy(raw_output)
    for r in sfx["rows"]:
        if r["class"] == "free_orientation_selector":
            r["T_squared"] = [1, 0, 1]
    for r in sraw["rows"]:
        if r.get("acyclic") and r["row_signature"]["dim"] == 5:
            r["T_squared_native"] = [1, 0, 1]
    out = _run_fixture(sfx, sraw, construction_packet)
    record("C5 self-inverse anchor", "invalid anchor", out["category"])

    # C6: anchor equal to neither x nor x^-1 -> disagreement, no row comparison.
    nfx = copy.deepcopy(fx)
    for r in nfx["rows"]:
        if r["class"] == "free_orientation_selector":
            r["T_squared"] = (QPhi.from_triple(r["T_squared"]) * QPhi(7, 0)).to_triple()
    out = _run_fixture(nfx, raw_output, construction_packet)
    record(
        "C6 anchor matches neither orientation",
        "disagreement",
        out["category"],
        {"row_mismatches_issued": out.get("row_mismatches") is not None},
    )

    # C7: identity expected tampered, rows untouched -> identity layer alone goes red.
    tfx = copy.deepcopy(fx)
    tfx["identities"][2]["expected"] = (
        QPhi.from_triple(tfx["identities"][2]["expected"]) * QPhi(5, 0)
    ).to_triple()
    out = _run_fixture(tfx, raw_output, construction_packet)
    record(
        "C7 identity expected tampered",
        "partial disagreement",
        out["category"],
        {
            "row_mismatches": len(out.get("row_mismatches") or []),
            "identities_unequal": [
                i["slot"] for i in (out.get("identities") or []) if not i["equal"]
            ],
        },
    )

    # C8: refusals before comparison.
    afx = copy.deepcopy(fx)
    afx["adjudicates"]["construction_packet_sha256"] = "0" * 64
    out = _run_fixture(afx, raw_output, construction_packet)
    record("C8a adjudicates tampered", "structural failure", out["category"])
    good = canonical_bytes(fx)
    pin = sha256_bytes(good)
    tampered = copy.deepcopy(fx)
    tampered["rows"][1]["T_squared"] = (
        QPhi.from_triple(tampered["rows"][1]["T_squared"]) * QPhi(2, 0)
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
    bad_cm["convention_map"]["basing_evaluation"] = "g |-> rho(g^-1)"
    out = _run_fixture(bad_cm, raw_output, construction_packet)
    record("C8e basing_evaluation not verbatim refused", "structural failure", out["category"])
    para = copy.deepcopy(fx)
    para["convention_map"][
        "bridge"
    ] = "|\\tau_\\rho|^2 with no field norm to a subfield (section 5.4)"
    out = _run_fixture(para, raw_output, construction_packet)
    record(
        "C8f bridge prose is recorded, never gated (a paraphrase of section 5.4 passes)",
        ("reproduced", para["convention_map"]["bridge"]),
        (out["category"], (out.get("convention_map") or {}).get("bridge")),
    )
    # An IDENTICAL duplicate: without the hook the object would still equal the pinned one
    # and be accepted through the fallback, so this control is diagnostic of the hook.
    dup = canonical_bytes(fx)[:-2] + b',\n  "format_version": "m8_8-answer-FIXTURE-1"\n}\n'
    try:
        load_answer_packet_bytes(dup, sha256_bytes(canonical_bytes(fx)))
        record("C8g duplicate key in delivered bytes refused", "refused", "accepted")
    except Refusal:
        record("C8g duplicate key in delivered bytes refused", "refused", "refused")
    greek = copy.deepcopy(fx)
    greek["convention_map"]["bridge"] = "T²_target(ρ) := |τ_ρ|², section 5.4"
    utf8 = (json.dumps(greek, sort_keys=True, indent=2, ensure_ascii=False) + "\n").encode("utf-8")
    _, rec = load_answer_packet_bytes(utf8, sha256_bytes(utf8))
    record(
        "C8h pin-matching bytes in another canonicalizer's rendering accepted, form recorded",
        (False, False),
        (rec.get("canonical_form_ok"), rec.get("canonicalization_changed_bytes")),
    )

    # C9: the selector class on a non-R7 row is refused, even with a correct table.
    sel9 = copy.deepcopy(fx)
    r7 = [r for r in sel9["rows"] if r["class"] == "free_orientation_selector"][0]
    r4 = [r for r in sel9["rows"] if r["row_signature"]["dim"] == 4][0]
    r7["class"], r4["class"] = "free", "free_orientation_selector"
    out = _run_fixture(sel9, raw_output, construction_packet)
    record("C9 selector class on a non-R7 row refused", "structural failure", out["category"])

    # C10: two labels sent to one committed signature is refused (injectivity).
    inj = copy.deepcopy(dfx)
    ents = inj["indexing_map"]["entries"]
    e1 = [e for e in ents if e["label"] == "R1"][0]
    e2 = [e for e in ents if e["label"] == "R2"][0]
    e2["row_signature"] = copy.deepcopy(e1["row_signature"])
    out = _run_fixture(inj, raw_output, construction_packet)
    record(
        "C10 non-injective map refused BY THE INJECTIVITY CHECK",
        ("structural failure", True),
        (out["category"], "not injective" in out.get("refusal", "")),
    )

    # C11: structurally trivial identities are refused (R0 factor; zero exponent; wrong kinds).
    triv = copy.deepcopy(fx)
    r0sig = [r for r in triv["rows"] if r["class"] == "declared_convention"][0]["row_signature"]
    triv["identities"] = [
        {
            "slot": f"x{i}",
            "factors": [{"row_signature": r0sig, "exponent": 0}] * 2,
            "expected": [1, 0, 1],
        }
        for i in range(4)
    ]
    out = _run_fixture(triv, raw_output, construction_packet)
    record("C11a trivial identities (R0^0) refused", "structural failure", out["category"])
    kinds = copy.deepcopy(fx)
    kinds["identities"][0] = copy.deepcopy(kinds["identities"][2])
    kinds["identities"][0]["slot"] = "third_product"
    out = _run_fixture(kinds, raw_output, construction_packet)
    record("C11b three products and one ratio refused", "structural failure", out["category"])

    # C12: a malformed but otherwise pinned packet yields a recorded structural failure,
    # never a traceback.
    mal = copy.deepcopy(fx)
    del mal["rows"][3]["T_squared"]
    out = _run_fixture(mal, raw_output, construction_packet)
    record(
        "C12 malformed row (missing T_squared) is a recorded refusal",
        ("structural failure", True),
        (out["category"], "refusal" in out),
    )

    # C13: EXACTNESS.  A comparison that were approximate, sign-blind or conjugation-blind
    # would pass C1; these three mutations of one free cell each require exactly one
    # mismatch: negation, Galois conjugation, and a near-miss rational.
    target = [r for r in fx["rows"] if r["class"] == "free"][2]
    v = QPhi.from_triple(target["T_squared"])
    near = QPhi(v.x + Fraction(1, 10**20), v.y)  # beyond any float tolerance
    for name, mutant in (
        ("negation", QPhi(-v.x, -v.y)),
        ("Galois conjugation", v.conjugate()),
        ("near-miss rational (1e-9 off)", near),
    ):
        if mutant == v:
            record(f"C13 exactness: {name} is a real mutation", True, False)
            continue
        ex = copy.deepcopy(fx)
        [r for r in ex["rows"] if r["label"] == target["label"]][0][
            "T_squared"
        ] = mutant.to_triple()
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
        "conjugate flag inside a ratio": "conjugate flag inside",
        "conjugate flag as a string": "JSON boolean",
        "exponent beyond the bound": "exceeds the bound",
    }

    def with_ids(mutator, msg):
        z = copy.deepcopy(fx)
        mutator(z["identities"])
        out = _run_fixture(z, raw_output, construction_packet)
        return out["category"], msg in out.get("refusal", "")

    dim3 = [r["row_signature"] for r in fx["rows"] if r["row_signature"]["dim"] == 3]
    r0sig = [r for r in fx["rows"] if r["class"] == "declared_convention"][0]["row_signature"]

    def non_galois(ids):  # dim-2 over dim-3: same exponents, not a pair
        ids[0]["factors"][1]["row_signature"] = dim3[0]

    def overlap(ids):  # product b re-uses a row of product a
        ids[3]["factors"][0]["row_signature"] = ids[2]["factors"][0]["row_signature"]

    def r0_factor(ids):
        ids[2]["factors"].append({"row_signature": r0sig, "exponent": 1})

    def zero_exp(ids):
        ids[2]["factors"][0]["exponent"] = 0

    def repeated(ids):
        ids[2]["factors"].append(copy.deepcopy(ids[2]["factors"][0]))

    def dup_pair(ids):  # both ratio slots on the same Galois pair
        ids[1]["factors"] = copy.deepcopy(ids[0]["factors"])

    def conj_in_ratio(ids):
        ids[0]["factors"][0]["conjugate"] = True

    def conj_string(ids):  # truthiness trap: the string "false"
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
        ("conjugate flag inside a ratio", conj_in_ratio),
        ("conjugate flag as a string", conj_string),
        ("exponent beyond the bound", exp_nine),
    ):
        record(
            f"C14 identity rule: {name} refused BY ITS OWN CHECK",
            ("structural failure", True),
            with_ids(mut, C14_MSG[name]),
        )
    # A genuine conjugate factor in a product, expected recomputed accordingly, is GREEN.
    cj = copy.deepcopy(fx)
    f0 = cj["identities"][2]["factors"][0]
    f0["conjugate"] = True
    key0 = signature_key(f0["row_signature"])
    tbl = {signature_key(r["row_signature"]): QPhi.from_triple(r["T_squared"]) for r in cj["rows"]}
    acc = ONE
    for f in cj["identities"][2]["factors"]:
        val = tbl[signature_key(f["row_signature"])]
        acc = acc * (val.conjugate() if signature_key(f["row_signature"]) == key0 else val)
    cj["identities"][2]["expected"] = acc.to_triple()
    out = _run_fixture(cj, raw_output, construction_packet)
    record("C14 genuine conjugate factor in a product is honored", "reproduced", out["category"])

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
    own99["rows"][2]["row_signature"]["dim"] = 99
    out = _run_fixture(own99, raw_output, construction_packet)
    record(
        "C16b own signature absent from the committed output refused BY THAT CHECK",
        ("structural failure", True),
        (out["category"], "own signature" in out.get("refusal", "")),
    )
    # C16c: T^2(R0) = 1 is a declared convention (section 9 item 2); any other value refuses.
    r0v = copy.deepcopy(fx)
    [r for r in r0v["rows"] if r["class"] == "declared_convention"][0]["T_squared"] = [7, 0, 1]
    out = _run_fixture(r0v, raw_output, construction_packet)
    record(
        "C16c declared_convention row carrying a value other than 1 refused",
        ("structural failure", True),
        (out["category"], "is not 1" in out.get("refusal", "")),
    )
    return results


# ----------------------------------------------------------------------------------------
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
    print(f"category     {result['category'].upper()}")
    print(f"orientation  {result.get('orientation')}")
    for m in result.get("row_mismatches") or []:
        print(
            f"  mismatch {m['label']} {m['signature']}: reference {m['reference']} "
            f"vs selected committed {m['selected_committed']}"
        )
    for i in result.get("identities") or []:
        print(
            f"  identity {i['slot']}: {'equal' if i['equal'] else 'UNEQUAL'} "
            f"recomputed {i['recomputed']} expected {i['expected']}"
        )
    print(f"written      {args.json}")
    if result["category"] in SUCCESS:
        return 0
    if result["category"] in FINDINGS:
        return 2
    return 1


if __name__ == "__main__":
    sys.exit(main())
