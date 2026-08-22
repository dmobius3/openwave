"""Exact Q(phi) and quaternion arithmetic, and the 2I closure with canonical IDs."""
from fractions import Fraction as F
import json, pathlib, re

# Q(phi): a + b*phi, phi^2 = phi + 1
def qadd(x, y):  return (x[0]+y[0], x[1]+y[1])
def qsub(x, y):  return (x[0]-y[0], x[1]-y[1])
def qmulf(x, y):
    a1,b1 = x; a2,b2 = y
    return (a1*a2 + b1*b2, a1*b2 + a2*b1 + b1*b2)
def qneg(x): return (-x[0], -x[1])
ZERO, ONE = (F(0),F(0)), (F(1),F(0))

# quaternions over Q(phi): 4-tuples of Q(phi) elements
def hmul(p, q):
    a1,b1,c1,d1 = p; a2,b2,c2,d2 = q
    return (
      qsub(qsub(qsub(qmulf(a1,a2), qmulf(b1,b2)), qmulf(c1,c2)), qmulf(d1,d2)),
      qadd(qadd(qadd(qmulf(a1,b2), qmulf(b1,a2)), qmulf(c1,d2)), qneg(qmulf(d1,c2))),
      qadd(qadd(qsub(qmulf(a1,c2), qmulf(b1,d2)), qmulf(c1,a2)), qmulf(d1,b2)),
      qadd(qadd(qadd(qmulf(a1,d2), qmulf(b1,c2)), qneg(qmulf(c1,b2))), qmulf(d1,a2)))
def hconj(p): return (p[0], qneg(p[1]), qneg(p[2]), qneg(p[3]))
HONE = (ONE, ZERO, ZERO, ZERO)

def parse_packet(path):
    """A-packet components are (a + b*phi)/2 with integer a, b."""
    d = json.loads(pathlib.Path(path).read_text())
    out = []
    for g in d["generators"]:
        comps = []
        for c in g:
            m = re.fullmatch(r"\((-?\d+) \+ (-?\d+)\*phi\)/(\d+)", c.strip())
            if not m:
                raise ValueError(f"unparsable component {c!r}")
            a, b, den = int(m.group(1)), int(m.group(2)), int(m.group(3))
            comps.append((F(a, den), F(b, den)))
        out.append(tuple(comps))
    return out, d

def close_group(gens):
    elems = [HONE]
    frontier = [HONE]
    while frontier:
        nxt = []
        for x in frontier:
            for g in gens:
                y = hmul(x, g)
                if y not in elems:
                    elems.append(y); nxt.append(y)
        frontier = nxt
    return elems

def canonical_ids(elems):
    """IDs by lexicographic rank of the exact coordinate tuple, per the protocol."""
    key = lambda e: tuple((c[0], c[1]) for c in e)
    order = sorted(elems, key=key)
    return order, {key(e): i for i, e in enumerate(order)}
