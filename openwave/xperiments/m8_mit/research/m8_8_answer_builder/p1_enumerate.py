"""P1 step 1: rebuild the 120-element enumeration from the group packet alone,
and check it against the digest pinned in protocol § 4.2 before anything is built
on top of it.

Target-free. It handles group elements and their IDs, never a torsion value.

The § 4.2 warning this implements: the ID sort key fixes the denominator at 2 and
does NOT reduce. Applying the § 4.4 value normalization here instead renames 26
of the 120 IDs, and the misread is quiet. So the digest is checked FIRST, and a
mismatch is fatal rather than advisory.
"""

import hashlib
import json
import pathlib
import re
import sys

from qphi_exact import Phi, Quat

PACKET_DIGEST = "27ff780d28d5d854d464ead87e8fc20244fac8334bda9f0600c6ee1b3c89561e"
PACKET_BYTES = 2389
PINNED_RANKS = {0: [-2, 0, 0, 0, 0, 0, 0, 0],
                118: [1, 0, 1, 0, 1, 0, 1, 0],
                119: [2, 0, 0, 0, 0, 0, 0, 0]}

_COMP = re.compile(r"^\(\s*(-?\d+)\s*\+\s*(-?\d+)\s*\*\s*phi\s*\)\s*/\s*2$")


def parse_component(s):
    """'(A + B*phi)/2' -> Phi. The group packet's declared coefficient_form."""
    m = _COMP.match(s.strip())
    if not m:
        raise ValueError(f"component does not match the declared form: {s!r}")
    from fractions import Fraction as F
    return Phi(F(int(m.group(1)), 2), F(int(m.group(2)), 2))


def load_generators(packet):
    if packet.get("coefficient_form", "").split(" with ")[0] != "(a + b*phi)/2":
        raise ValueError("coefficient_form is not the expected fixed-denominator form")
    if packet.get("quaternion_basis") != ["1", "i", "j", "k"]:
        raise ValueError("quaternion_basis is not (1, i, j, k); the sort key assumes it")
    return [Quat([parse_component(c) for c in g]) for g in packet["generators"]]


def close_group(gens, cap=200):
    """Multiplicative closure. Cap is a runaway guard, not an expected size."""
    ident = Quat([Phi(1), Phi(0), Phi(0), Phi(0)])
    seen = {ident}
    frontier = [ident]
    while frontier:
        nxt = []
        for x in frontier:
            for g in gens:
                y = x * g
                if y not in seen:
                    seen.add(y)
                    nxt.append(y)
                    if len(seen) > cap:
                        raise RuntimeError(f"closure exceeded {cap}; not a finite group here")
        frontier = nxt
    return seen


def enumerate_group(packet):
    """Returns the elements in § 4.2 rank order. Index == element ID."""
    elems = close_group(load_generators(packet))
    return sorted(elems, key=lambda q: q.key())


def digest_of(ordered):
    """SHA-256 over one JSON array of the 120 rank-ordered eight-integer arrays,
    no whitespace, bare decimal integers, ASCII, no trailing newline (§ 4.2)."""
    payload = json.dumps([list(q.key()) for q in ordered],
                         separators=(",", ":"), ensure_ascii=True)
    b = payload.encode("ascii")
    return hashlib.sha256(b).hexdigest(), len(b)


def check(packet_path):
    packet = json.loads(pathlib.Path(packet_path).read_text())
    ordered = enumerate_group(packet)
    got, nbytes = digest_of(ordered)

    results = []
    results.append(("order 120", 120, len(ordered)))
    results.append(("digest", PACKET_DIGEST, got))
    results.append(("byte length", PACKET_BYTES, nbytes))
    for rank, key in PINNED_RANKS.items():
        results.append((f"rank {rank} key", key,
                        list(ordered[rank].key()) if rank < len(ordered) else None))

    ok = all(exp == act for _, exp, act in results)
    return ordered, results, ok


def main():
    if len(sys.argv) != 2:
        print("usage: p1_enumerate.py <group_packet.json>", file=sys.stderr)
        return 2
    ordered, results, ok = check(sys.argv[1])
    print("P1 step 1: enumeration rebuilt from the group packet alone")
    for name, exp, act in results:
        mark = "ok  " if exp == act else "FAIL"
        shown = f"{act}" if exp == act else f"{act}   expected {exp}"
        print(f"  {mark} {name:16s} {shown}")
    print(f"\n  {'ENUMERATION REPRODUCED' if ok else 'ENUMERATION DOES NOT MATCH; STOP'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
