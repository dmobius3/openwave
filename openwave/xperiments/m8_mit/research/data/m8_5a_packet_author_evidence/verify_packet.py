#!/usr/bin/env python3
"""
M8.5-A clean-room packet: exact verification of the supplied group input.

Author-side verification. Confirms that the generators in the packet are unit
quaternions over Q(sqrt5) and that they close to the binary icosahedral group
2I, and NOT merely to some order-120 subgroup of SU(2).

Exact arithmetic throughout: coefficients live in Q(phi) with phi^2 = phi + 1,
represented as Fraction pairs (p, q) meaning p + q*phi. No floating point is
used anywhere, so there is no tolerance, no dedup ambiguity, and no closure
ambiguity.

GATES (each printed PASS line is mutation-tested; see --mutation-tests)

  P1  every supplied generator has quaternion norm exactly 1
  P2  the generators close under exact multiplication to exactly 120 distinct
      elements
  P3  the identity is present, every element's inverse is in the set, and the
      center is exactly {+1, -1}
  P4  the quotient by the center has order 60 and the element-order profile of
      A_5 (1 of order 1, 15 of order 2, 20 of order 3, 24 of order 5, none of
      order 4 or 6)

P4 is the gate that distinguishes 2I from another order-120 subgroup: order and
closure alone do not, which is why P2 is not sufficient on its own.

WHAT THIS DOES NOT DO: it computes no characters, no irreducible dimensions, no
McKay data, and no first-occurrence anything. It certifies the packet's group
input and nothing downstream of it. The verification report carries gate
outcomes and hashes only.

Usage:
  python3 verify_packet.py --packet m8_5a_packet.json [--report report.json]
  python3 verify_packet.py --packet m8_5a_packet.json --mutation-tests
"""

import argparse
import hashlib
import json
import platform
import sys
from fractions import Fraction
from itertools import permutations

# --------------------------------------------------------------------------
# exact arithmetic in Q(phi), phi^2 = phi + 1
# an element is a pair (p, q) of Fractions meaning p + q*phi
# --------------------------------------------------------------------------

ZERO = (Fraction(0), Fraction(0))
ONE = (Fraction(1), Fraction(0))


def qadd(a, b):
    return (a[0] + b[0], a[1] + b[1])


def qsub(a, b):
    return (a[0] - b[0], a[1] - b[1])


def qneg(a):
    return (-a[0], -a[1])


def qmul(a, b):
    # (p1 + q1 phi)(p2 + q2 phi) = p1p2 + (p1q2 + q1p2) phi + q1q2 phi^2
    # phi^2 = phi + 1
    p = a[0] * b[0] + a[1] * b[1]
    q = a[0] * b[1] + a[1] * b[0] + a[1] * b[1]
    return (p, q)


def parse_coeff(s):
    """Parse the canonical form '(a + b*phi)/d' into a Q(phi) element."""
    txt = s.strip()
    if not (txt.startswith("(") and "/" in txt):
        raise ValueError(f"coefficient not in canonical form: {s!r}")
    body, denom = txt.rsplit("/", 1)
    body = body.strip()
    if not (body.startswith("(") and body.endswith(")")):
        raise ValueError(f"coefficient not in canonical form: {s!r}")
    body = body[1:-1]
    d = int(denom.strip())
    if d <= 0:
        raise ValueError(f"denominator must be positive: {s!r}")
    # body is 'a + b*phi' or 'a - b*phi'
    if "+" in body:
        a_txt, b_txt = body.split("+", 1)
        sign = 1
    elif "-" in body[1:]:
        idx = body.index("-", 1)
        a_txt, b_txt = body[:idx], body[idx + 1:]
        sign = -1
    else:
        raise ValueError(f"coefficient not in canonical form: {s!r}")
    a = int(a_txt.strip())
    b_txt = b_txt.strip()
    if not b_txt.endswith("*phi"):
        raise ValueError(f"coefficient not in canonical form: {s!r}")
    b = sign * int(b_txt[:-4].strip())
    return (Fraction(a, d), Fraction(b, d))


def format_coeff(x):
    """Canonical '(a + b*phi)/d': integers, d > 0, gcd(|a|,|b|,d) = 1."""
    p, q = x
    d = _lcm(p.denominator, q.denominator)
    a = int(p * d)
    b = int(q * d)
    g = _gcd3(abs(a), abs(b), d)
    if g:
        a, b, d = a // g, b // g, d // g
    sign = "+" if b >= 0 else "-"
    return f"({a} {sign} {abs(b)}*phi)/{d}"


def _gcd(a, b):
    while b:
        a, b = b, a % b
    return a


def _gcd3(a, b, c):
    return _gcd(_gcd(a, b), c)


def _lcm(a, b):
    return a * b // _gcd(a, b)


# --------------------------------------------------------------------------
# quaternions with Q(phi) coefficients, basis (1, i, j, k)
# --------------------------------------------------------------------------

QONE = (ONE, ZERO, ZERO, ZERO)


def hmul(x, y):
    w1, x1, y1, z1 = x
    w2, x2, y2, z2 = y
    w = qsub(qsub(qsub(qmul(w1, w2), qmul(x1, x2)), qmul(y1, y2)), qmul(z1, z2))
    i = qadd(qadd(qmul(w1, x2), qmul(x1, w2)), qsub(qmul(y1, z2), qmul(z1, y2)))
    j = qadd(qadd(qmul(w1, y2), qmul(y1, w2)), qsub(qmul(z1, x2), qmul(x1, z2)))
    k = qadd(qadd(qmul(w1, z2), qmul(z1, w2)), qsub(qmul(x1, y2), qmul(y1, x2)))
    return (w, i, j, k)


def hconj(x):
    w, i, j, k = x
    return (w, qneg(i), qneg(j), qneg(k))


def hnorm(x):
    """Quaternion norm w^2 + x^2 + y^2 + z^2, exact."""
    total = ZERO
    for c in x:
        total = qadd(total, qmul(c, c))
    return total


def hneg(x):
    return tuple(qneg(c) for c in x)


def horder(x, cap=600):
    """Multiplicative order of a unit quaternion, exact."""
    cur = x
    for n in range(1, cap + 1):
        if cur == QONE:
            return n
        cur = hmul(cur, x)
    raise ValueError("order exceeded cap")


def close_group(gens, cap=100000):
    """Exact closure under multiplication. Returns the element set."""
    elems = {QONE}
    frontier = [QONE]
    while frontier:
        nxt = []
        for a in frontier:
            for g in gens:
                b = hmul(a, g)
                if b not in elems:
                    elems.add(b)
                    nxt.append(b)
                    if len(elems) > cap:
                        raise ValueError("closure exceeded cap")
        frontier = nxt
    return elems


# --------------------------------------------------------------------------
# gates
# --------------------------------------------------------------------------

A5_PROFILE = {1: 1, 2: 15, 3: 20, 5: 24}


def run_gates(gens, mut=None):
    res = {}

    # P1: unit norms
    norms_ok = all(hnorm(g) == ONE for g in gens)
    res["P1_generators_unit_norm"] = norms_ok

    # P2: exact closure to 120
    try:
        G = close_group(gens)
    except ValueError:
        res["P2_closure_order_120"] = False
        res["P3_identity_inverses_center"] = False
        res["P4_central_quotient_is_A5"] = False
        return res, None
    order = len(G)
    res["P2_closure_order_120"] = (order == 120)

    # P3: identity, inverses, center exactly {+1, -1}
    has_id = QONE in G
    inv_ok = all(hconj(g) in G for g in G)  # unit quaternion inverse = conjugate
    center = {g for g in G if all(hmul(g, h) == hmul(h, g) for h in G)}
    minus_one = hneg(QONE)
    center_ok = (center == {QONE, minus_one})
    res["P3_identity_inverses_center"] = bool(has_id and inv_ok and center_ok)

    # P4: quotient by center has A_5 order profile
    if mut == "quotient_by_trivial":
        cosets = [frozenset({g}) for g in G]
    else:
        cosets = []
        seen = set()
        for g in G:
            if g in seen:
                continue
            pair = frozenset({g, hneg(g)})
            cosets.append(pair)
            seen |= pair
    profile = {}
    if len(cosets) == 60:
        rep_of = {}
        for c in cosets:
            for e in c:
                rep_of[e] = c
        for c in cosets:
            g = next(iter(c))
            n = 1
            cur = g
            while rep_of[cur] != rep_of[QONE] or n == 0:
                cur = hmul(cur, g)
                n += 1
                if n > 120:
                    break
            profile[n] = profile.get(n, 0) + 1
    res["P4_central_quotient_is_A5"] = (len(cosets) == 60 and profile == A5_PROFILE)
    return res, order


# --------------------------------------------------------------------------
# packet io
# --------------------------------------------------------------------------

def load_packet(path):
    raw = open(path, "rb").read()
    digest = hashlib.sha256(raw).hexdigest()
    doc = json.loads(raw.decode("utf-8"))
    gens = []
    for row in doc["generators"]:
        if len(row) != 4:
            raise ValueError("generator is not a 4-tuple")
        gens.append(tuple(parse_coeff(c) for c in row))
    return doc, gens, digest


MUTATIONS = {
    "perturb_generator": "add 1 to the scalar numerator of generator 0",
    "drop_golden_generator": "keep only the rational generator",
    "duplicate_generator": "replace generator 1 by generator 0",
    "quotient_by_trivial": "quotient by the trivial subgroup instead of the center",
}


def mutate(gens, name):
    if name == "perturb_generator":
        g = list(gens[0])
        p, q = g[0]
        g[0] = (p + 1, q)
        return [tuple(g)] + list(gens[1:])
    if name == "drop_golden_generator":
        return [gens[0]]
    if name == "duplicate_generator":
        return [gens[0], gens[0]]
    if name == "quotient_by_trivial":
        return list(gens)
    raise ValueError(name)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--packet", required=True)
    ap.add_argument("--report")
    ap.add_argument("--mutation-tests", action="store_true")
    args = ap.parse_args()

    doc, gens, digest = load_packet(args.packet)
    res, order = run_gates(gens)

    print(f"packet          : {args.packet}")
    print(f"packet sha256   : {digest}")
    print(f"generators      : {len(gens)}")
    print(f"closure order   : {order}")
    print("")
    for k, v in res.items():
        print(f"  {'PASS' if v else 'FAIL'}  {k}")
    baseline_ok = all(res.values())
    print("")
    print("VERDICT:", "all gates PASS" if baseline_ok else "GATE FAILURE")

    if args.mutation_tests:
        print("\nmutation tests (each must redden at least one gate):")
        if not baseline_ok:
            print("  baseline is not green; aborting")
            return 1
        covered = set()
        allred = True
        for name, desc in MUTATIONS.items():
            mg = mutate(gens, name)
            mres, _ = run_gates(mg, mut=name if name == "quotient_by_trivial" else None)
            red = [k for k, v in mres.items() if not v]
            covered |= set(red)
            ok = bool(red)
            allred &= ok
            print(f"  {'OK  ' if ok else 'MISS'} {name:24s} reddened: {','.join(red) if red else 'NOTHING'}  ({desc})")
        uncovered = set(res) - covered
        if uncovered:
            print(f"  UNCOVERED GATES: {sorted(uncovered)}")
        if not allred or uncovered:
            print("\nMUTATION SUITE FAILED")
            return 1
        print("\nmutation suite: every gate reddened by at least one mutation")

    if args.report:
        report = {
            "packet_sha256": digest,
            "packet_format_version": doc.get("format_version"),
            "generator_count": len(gens),
            "closure_order": order,
            "gates": res,
            "environment": {
                "python": platform.python_version(),
                "platform": platform.platform(),
                "exact_arithmetic": "fractions.Fraction over Q(phi), no floating point",
            },
        }
        with open(args.report, "w") as fh:
            json.dump(report, fh, indent=2, sort_keys=True)
            fh.write("\n")
        print(f"\nreport written: {args.report}")

    return 0 if baseline_ok else 1


if __name__ == "__main__":
    sys.exit(main())
