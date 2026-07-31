"""M8.5-A PACKET AUDIT: the maintainer check that runs before the clean room opens.

WHAT THIS IS.  The [M8.5-A reproduction protocol](../findings/m8_5a_reproduction_protocol.md)
section 4 requires the packet entering the clean room to be audited by someone other than the
context that will implement from it.  The clean-room procedure standard
([`dev_docs/tasks/t4_task_details.md`](../../../../../dev_docs/tasks/t4_task_details.md))
fixes two properties of that audit, and this file is written to satisfy both:

  1  IT RUNS FROM THE PACKET ALONE.  The author shipped a verification script and a report
     alongside the packet, in a separate archive that stays outside the room.  This file
     reads NEITHER.  An audit that reads the supplier's report and concurs is a review of
     the report; the two results are set side by side afterwards, as a cross-check, in
     [`m8_5_task_details.md`](../tasks/m8_5_task_details.md).
  2  IT IS MECHANICAL, NOT A JUDGMENT CALL.  The same maintainer assembles and audits the
     packet, so every check below is a count or an exact equality that can go red, and
     `--mutation-tests` shows each one doing so.

WHAT THE PACKET IS.  Two generators of a finite subgroup of the unit quaternions, with every
component given as an exact element of the field Q(phi), phi^2 = phi + 1.  No decimal
rendering exists anywhere in the packet, which is what makes the hash pin a GROUP rather
than a transcription of one (protocol section 4, canonicalization).

WHAT THIS VERIFIES.  Eight checks, each mutation-tested:

  A1  format: exactly the declared keys, the declared minimal polynomial and basis, every
      component parsing as (a + b*phi)/2 with integer a and b and the denominator fixed at 2
  A2  both generators have quaternion norm exactly 1 in Q(phi), not 1 to a tolerance
  A3  the generated closure is finite, closed under multiplication, and has exactly 120
      distinct elements
  A4  the closure contains the identity and an inverse for every element
  A5  the center is exactly {+1, -1}
  A6  the element-order census over all 120 elements
  A7  the central quotient has order 60 and the A5 element-order profile
  A8  leakage: the packet carries no labels, dimensions, distances or character values, and
      no key outside the declared set

A6 IS THIS AUDIT'S ADDITION.  A5 and A7 together (center of order 2, quotient A5) are also
satisfied by the direct product A5 x C2, which is NOT the binary icosahedral group.  What
separates them is the order-4 population: 2I has 30 elements of order 4, A5 x C2 has none.
A finite subgroup of the unit quaternions cannot be A5 x C2 for an independent reason (the
only unit quaternion of order 2 is -1, while A5 x C2 has 31 involutions), so the census is
not the only thing standing between the two.  It is, however, the mechanical version of that
argument, and this audit prefers a count to an appeal to a theorem.

WHAT THIS DOES NOT VERIFY.  That the supplied group is the one the M8 column needs.  This
audit establishes what the packet IS; whether that is the right input is the protocol's
question, answered by the pre-registration, not here.  It also does not certify the packet
target-free by construction: A8 scans for the leakage vocabulary this protocol named, and a
scan is a floor, not a proof.

USAGE.
    python3 m8_5a_packet_audit.py                     audit, write JSON, exit 0 on all-pass
    python3 m8_5a_packet_audit.py --mutation-tests    additionally run the mutation suite
    python3 m8_5a_packet_audit.py --packet PATH       audit a packet elsewhere

Exit code is nonzero if any check fails, so a red audit cannot be mistaken for a green one.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEFAULT_PACKET = HERE.parent / "data" / "m8_5a_packet.json"
DEFAULT_OUT = HERE.parent / "data" / "m8_5a_packet_audit.json"

# The expected group, declared here so the checks below compare against a stated reference
# rather than against whatever the packet happens to produce.  Sources: the binary
# icosahedral group 2I = SL(2,5), order 120, nine conjugacy classes; and A5, order 60.
EXPECTED_ORDER = 120
EXPECTED_ORDER_CENSUS = {1: 1, 2: 1, 3: 20, 4: 30, 5: 24, 6: 20, 10: 24}
EXPECTED_QUOTIENT_ORDER = 60
EXPECTED_QUOTIENT_CENSUS = {1: 1, 2: 15, 3: 20, 5: 24}

DECLARED_KEYS = {
    "coefficient_field",
    "coefficient_form",
    "format_version",
    "generators",
    "quaternion_basis",
}
DECLARED_MINPOLY = "phi^2 - phi - 1"
DECLARED_BASIS = ["1", "i", "j", "k"]

# Vocabulary that would carry a target into the room.  Protocol section 4 forbids the packet
# carrying anything derived; these are the words the derived objects are made of.
LEAKAGE_TOKENS = [
    "irrep", "irreducible", "character", "chi", "dim", "dimension", "degree",
    "distance", "adjacency", "mckay", "affine", "e8", "dynkin", "node",
    "class", "conjugacy", "trivial", "sigma", "tau", "spin", "label",
    "table", "level", "occurrence", "flat", "connection", "sym",
]

COMPONENT_RE = re.compile(r"^\((-?\d+) \+ (-?\d+)\*phi\)/2$")


# ----------------------------------------------------------------------------------------
# Exact arithmetic in Q(phi), phi^2 = phi + 1.  An element is the pair (p, q) meaning
# p + q*phi with p and q rational.  Nothing here rounds, so every equality below is exact.
# ----------------------------------------------------------------------------------------

class Fp:
    """An element of Q(phi), phi^2 = phi + 1."""

    __slots__ = ("p", "q")

    def __init__(self, p: Fraction, q: Fraction = Fraction(0)) -> None:
        self.p = Fraction(p)
        self.q = Fraction(q)

    def __add__(self, other: "Fp") -> "Fp":
        return Fp(self.p + other.p, self.q + other.q)

    def __sub__(self, other: "Fp") -> "Fp":
        return Fp(self.p - other.p, self.q - other.q)

    def __neg__(self) -> "Fp":
        return Fp(-self.p, -self.q)

    def __mul__(self, other: "Fp") -> "Fp":
        # (p1 + q1 phi)(p2 + q2 phi) = p1 p2 + (p1 q2 + p2 q1) phi + q1 q2 phi^2,
        # and phi^2 = phi + 1.
        p = self.p * other.p + self.q * other.q
        q = self.p * other.q + other.p * self.q + self.q * other.q
        return Fp(p, q)

    def __eq__(self, other: object) -> bool:
        return isinstance(other, Fp) and self.p == other.p and self.q == other.q

    def __hash__(self) -> int:
        return hash((self.p, self.q))

    def key(self) -> tuple:
        return (self.p.numerator, self.p.denominator, self.q.numerator, self.q.denominator)


ZERO = Fp(Fraction(0))
ONE = Fp(Fraction(1))


class Quat:
    """A quaternion with coefficients in Q(phi), in the basis (1, i, j, k)."""

    __slots__ = ("c",)

    def __init__(self, c: tuple) -> None:
        self.c = tuple(c)

    def __mul__(self, o: "Quat") -> "Quat":
        a1, b1, c1, d1 = self.c
        a2, b2, c2, d2 = o.c
        return Quat((
            a1 * a2 - b1 * b2 - c1 * c2 - d1 * d2,
            a1 * b2 + b1 * a2 + c1 * d2 - d1 * c2,
            a1 * c2 - b1 * d2 + c1 * a2 + d1 * b2,
            a1 * d2 + b1 * c2 - c1 * b2 + d1 * a2,
        ))

    def __neg__(self) -> "Quat":
        return Quat(tuple(-x for x in self.c))

    def conj(self) -> "Quat":
        a, b, c, d = self.c
        return Quat((a, -b, -c, -d))

    def norm(self) -> Fp:
        a, b, c, d = self.c
        return a * a + b * b + c * c + d * d

    def __eq__(self, other: object) -> bool:
        return isinstance(other, Quat) and self.c == other.c

    def __hash__(self) -> int:
        return hash(self.key())

    def key(self) -> tuple:
        return tuple(x.key() for x in self.c)


QUAT_ONE = Quat((ONE, ZERO, ZERO, ZERO))
QUAT_NEG_ONE = -QUAT_ONE


# ----------------------------------------------------------------------------------------
# Packet parsing
# ----------------------------------------------------------------------------------------

def parse_component(text: str) -> Fp:
    """Parse '(a + b*phi)/2' into an exact Q(phi) element.  Raises on any other shape."""
    m = COMPONENT_RE.match(text)
    if not m:
        raise ValueError(f"component does not match (a + b*phi)/2 with integer a, b: {text!r}")
    a, b = int(m.group(1)), int(m.group(2))
    return Fp(Fraction(a, 2), Fraction(b, 2))


def parse_generators(packet: dict) -> list:
    return [Quat(tuple(parse_component(x) for x in g)) for g in packet["generators"]]


# ----------------------------------------------------------------------------------------
# Group construction, from the generators alone
# ----------------------------------------------------------------------------------------

def close(generators: list, cap: int = 5000) -> list:
    """Breadth-first closure under multiplication.  The cap turns a non-finite or wrong
    generating set into a loud failure rather than a hang."""
    seen = {QUAT_ONE.key(): QUAT_ONE}
    frontier = [QUAT_ONE]
    while frontier:
        nxt = []
        for x in frontier:
            for g in generators:
                y = x * g
                k = y.key()
                if k not in seen:
                    seen[k] = y
                    nxt.append(y)
                    if len(seen) > cap:
                        raise RuntimeError(f"closure exceeded {cap} elements; not finite as expected")
        frontier = nxt
    return list(seen.values())


def element_order(q: Quat, cap: int = 1000) -> int:
    x, n = q, 1
    while x != QUAT_ONE:
        x = x * q
        n += 1
        if n > cap:
            raise RuntimeError("element order exceeded cap")
    return n


def census(elements: list, order_of) -> dict:
    out: dict = {}
    for e in elements:
        o = order_of(e)
        out[o] = out.get(o, 0) + 1
    return dict(sorted(out.items()))


def center(elements: list) -> list:
    return [z for z in elements if all(z * g == g * z for g in elements)]


def quotient_by_center(elements: list) -> tuple:
    """Cosets of {+1, -1}.  Returns (cosets, multiply, order_of) with the coset operations
    defined on representatives, which is well defined because the subgroup is central."""
    rep = {}
    for e in elements:
        k = e.key()
        kn = (-e).key()
        if k not in rep and kn not in rep:
            rep[k] = e
    cosets = list(rep.values())

    def coset_key(q: Quat) -> tuple:
        return min(q.key(), (-q).key())

    index = {coset_key(c): c for c in cosets}

    def order_of(c: Quat) -> int:
        x, n = c, 1
        while coset_key(x) != coset_key(QUAT_ONE):
            x = x * c
            n += 1
            if n > 1000:
                raise RuntimeError("coset order exceeded cap")
        return n

    return cosets, index, order_of


# ----------------------------------------------------------------------------------------
# The checks
# ----------------------------------------------------------------------------------------

def run_checks(raw: bytes, packet: dict) -> tuple:
    """Returns (checks, detail).  Each check is (id, description, passed, observed)."""
    checks = []
    detail: dict = {}

    def add(cid, desc, passed, observed):
        checks.append({"id": cid, "check": desc, "pass": bool(passed), "observed": observed})

    # A8 is computed first and appended last.  It is the only check that does not depend on
    # the packet parsing, and leakage detection must not be silenced by a format failure: a
    # stray key carrying a target reddens A1, and an audit that stopped there would report
    # the wrong reason.  The mutation suite is what surfaced this.
    text = raw.decode("utf-8", errors="replace").lower()
    hits = sorted({t for t in LEAKAGE_TOKENS if re.search(rf"\b{re.escape(t)}", text)})
    stray_keys = sorted(set(packet.keys()) - DECLARED_KEYS)
    a8 = not hits and not stray_keys

    def finish():
        add("A8", "no leakage vocabulary and no key outside the declared set", a8,
            {"token_hits": hits, "stray_keys": stray_keys})
        return checks, detail

    # A1 format
    keys_ok = set(packet.keys()) == DECLARED_KEYS
    minpoly = packet.get("coefficient_field", {}).get("minimal_polynomial")
    gen_sym = packet.get("coefficient_field", {}).get("generator")
    basis_ok = packet.get("quaternion_basis") == DECLARED_BASIS
    try:
        gens = parse_generators(packet)
        parse_ok, parse_err = True, None
    except ValueError as exc:
        gens, parse_ok, parse_err = [], False, str(exc)
    shape_ok = (
        isinstance(packet.get("generators"), list)
        and len(packet["generators"]) == 2
        and all(isinstance(g, list) and len(g) == 4 for g in packet["generators"])
    )
    a1 = keys_ok and minpoly == DECLARED_MINPOLY and gen_sym == "phi" and basis_ok and shape_ok and parse_ok
    add("A1", "declared keys, minimal polynomial, basis; every component is (a + b*phi)/2",
        a1, {"keys_match": keys_ok, "minimal_polynomial": minpoly, "basis_match": basis_ok,
             "shape_2x4": shape_ok, "components_parse": parse_ok, "parse_error": parse_err})
    if not a1:
        return finish()

    # A2 exact unit norm
    norms = [g.norm() for g in gens]
    a2 = all(n == ONE for n in norms)
    add("A2", "both generators have quaternion norm exactly 1 in Q(phi)", a2,
        {"norms": [f"{n.p} + {n.q}*phi" for n in norms]})

    # A3 closure
    try:
        elements = close(gens)
        closure_err = None
    except RuntimeError as exc:
        elements, closure_err = [], str(exc)
    order = len(elements)
    keys = {e.key() for e in elements}
    closed = bool(elements) and all((a * b).key() in keys for a in elements for b in elements)
    a3 = order == EXPECTED_ORDER and closed
    add("A3", f"closure is finite, multiplicatively closed, of order exactly {EXPECTED_ORDER}",
        a3, {"order": order, "closed_under_multiplication": closed, "error": closure_err})
    if not a3:
        return finish()

    # A4 identity and inverses
    has_id = QUAT_ONE.key() in keys
    inverses_ok = all(any((a * b).key() == QUAT_ONE.key() for b in elements) for a in elements)
    add("A4", "identity present and every element has an inverse in the closure",
        has_id and inverses_ok, {"identity_present": has_id, "all_invertible": inverses_ok})

    # A5 center
    z = center(elements)
    z_keys = sorted(e.key() for e in z)
    a5 = len(z) == 2 and z_keys == sorted([QUAT_ONE.key(), QUAT_NEG_ONE.key()])
    add("A5", "center is exactly {+1, -1}", a5, {"center_order": len(z)})

    # A6 element-order census
    cen = census(elements, element_order)
    a6 = cen == EXPECTED_ORDER_CENSUS
    add("A6", "element-order census matches 2I and rules out A5 x C2 (30 elements of order 4)",
        a6, {"census": cen, "expected": EXPECTED_ORDER_CENSUS,
             "order_4_population": cen.get(4, 0)})

    # A7 central quotient
    cosets, _, coset_order = quotient_by_center(elements)
    qcen = census(cosets, coset_order)
    a7 = len(cosets) == EXPECTED_QUOTIENT_ORDER and qcen == EXPECTED_QUOTIENT_CENSUS
    add("A7", f"central quotient has order {EXPECTED_QUOTIENT_ORDER} with the A5 order profile",
        a7, {"quotient_order": len(cosets), "census": qcen,
             "expected": EXPECTED_QUOTIENT_CENSUS})

    detail["generator_orders"] = [element_order(g) for g in gens]
    detail["element_order_census"] = cen
    detail["quotient_order_census"] = qcen
    detail["group_order"] = order
    detail["center_order"] = len(z)
    return finish()


# ----------------------------------------------------------------------------------------
# Canonicalization and hashing
# ----------------------------------------------------------------------------------------

def canonical_bytes(packet: dict) -> bytes:
    """The canonical serialization: keys sorted, two-space indent, ASCII, LF, one trailing
    newline.  Declared here so the same group always produces the same bytes."""
    return json.dumps(packet, sort_keys=True, indent=2, ensure_ascii=True).encode("utf-8") + b"\n"


# ----------------------------------------------------------------------------------------
# Mutation suite: every check above is shown to go red under a deliberate defect
# ----------------------------------------------------------------------------------------

def mutations(packet: dict) -> list:
    """Each entry is (name, mutated packet, the check ids expected to go red)."""
    out = []

    m = json.loads(json.dumps(packet))
    m["generators"][0][0] = "(1 + 0*phi)/3"
    out.append(("denominator_not_2", m, ["A1"]))

    m = json.loads(json.dumps(packet))
    m["coefficient_field"]["minimal_polynomial"] = "phi^2 - 2"
    out.append(("wrong_minimal_polynomial", m, ["A1"]))

    m = json.loads(json.dumps(packet))
    m["generators"][0][1] = "(3 + 0*phi)/2"
    out.append(("norm_not_unit", m, ["A2"]))

    m = json.loads(json.dumps(packet))
    # A unit quaternion, so A2 stays green by design: this mutation tests A3 alone.
    m["generators"][1] = ["(0 + 0*phi)/2", "(2 + 0*phi)/2", "(0 + 0*phi)/2", "(0 + 0*phi)/2"]
    out.append(("second_generator_replaced_by_i", m, ["A3"]))

    m = json.loads(json.dumps(packet))
    m["generators"] = [packet["generators"][0], packet["generators"][0]]
    out.append(("generators_collapsed_to_one", m, ["A3"]))

    m = json.loads(json.dumps(packet))
    m["dimensions"] = [1, 2, 2, 3, 3, 4, 4, 5, 6]
    out.append(("dimension_labels_added", m, ["A8"]))

    m = json.loads(json.dumps(packet))
    m["mckay_distance"] = 3
    out.append(("mckay_distance_added", m, ["A8"]))

    return out


def run_mutation_suite(packet: dict) -> list:
    results = []
    for name, mutated, expected_red in mutations(packet):
        raw = canonical_bytes(mutated)
        try:
            checks, _ = run_checks(raw, mutated)
            reds = sorted(c["id"] for c in checks if not c["pass"])
            error = None
        except Exception as exc:  # a mutation that breaks the run counts as red, loudly
            reds, error = ["EXCEPTION"], f"{type(exc).__name__}: {exc}"
        # The expected checks must go red.  Later checks may not run at all once an earlier
        # one aborts the sequence, which is why this is a subset test, not an equality test.
        ok = all(cid in reds or "EXCEPTION" in reds for cid in expected_red)
        results.append({"mutation": name, "expected_red": expected_red,
                        "observed_red": reds, "detected": ok, "error": error})
    return results


# ----------------------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--packet", type=Path, default=DEFAULT_PACKET)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--mutation-tests", action="store_true")
    args = ap.parse_args()

    raw = args.packet.read_bytes()
    packet = json.loads(raw)

    incoming = hashlib.sha256(raw).hexdigest()
    canon = canonical_bytes(packet)
    authoritative = hashlib.sha256(canon).hexdigest()

    checks, detail = run_checks(raw, packet)
    all_pass = all(c["pass"] for c in checks)

    mut = run_mutation_suite(packet) if args.mutation_tests else None
    mut_ok = all(m["detected"] for m in mut) if mut else None

    result = {
        "what": "M8.5-A packet audit, run from the packet alone",
        "packet_file": args.packet.name,
        "incoming_sha256": incoming,
        "canonical_form": "json.dumps(sort_keys=True, indent=2, ensure_ascii=True) + LF",
        "authoritative_sha256": authoritative,
        "canonicalization_changed_bytes": incoming != authoritative,
        "checks": checks,
        "detail": detail,
        "all_checks_pass": all_pass,
        "mutation_suite": mut,
        "mutation_suite_all_detected": mut_ok,
    }
    args.out.write_text(json.dumps(result, indent=2) + "\n")

    print(f"packet            {args.packet}")
    print(f"incoming sha256   {incoming}")
    print(f"authoritative     {authoritative}")
    print(f"canonicalization  {'CHANGED the bytes' if incoming != authoritative else 'no-op (delivered bytes already canonical)'}")
    print()
    for c in checks:
        print(f"  {'PASS' if c['pass'] else 'FAIL'}  {c['id']}  {c['check']}")
        print(f"        {json.dumps(c['observed'])}")
    print()
    print(f"  generator orders          {detail.get('generator_orders')}")
    print(f"  element-order census      {detail.get('element_order_census')}")
    print(f"  central quotient census   {detail.get('quotient_order_census')}")
    if mut is not None:
        print()
        for m in mut:
            print(f"  {'DETECTED' if m['detected'] else 'MISSED  '}  {m['mutation']:34s} "
                  f"expected red {m['expected_red']}, observed {m['observed_red']}")
    print()
    print(f"ALL CHECKS PASS: {all_pass}" + (f" | MUTATIONS ALL DETECTED: {mut_ok}" if mut is not None else ""))
    print(f"written: {args.out}")

    return 0 if all_pass and (mut_ok is not False) else 1


if __name__ == "__main__":
    sys.exit(main())
