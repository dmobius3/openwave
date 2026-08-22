"""Mechanical parsing of the STATED closed forms recorded in the P1 artifact.

This is source verification, not rederivation. P3 does not run the M8.3 program,
does not touch a spectral computation, and does not reconcile anything: it reads
the closed form each source WROTE DOWN, evaluates it in exact Q(phi) with P3's
own arithmetic, and requires the result to equal the value the packet commits.

Why it has to exist. Without it, P3's cross-family claim collapses onto P1's own
verdict bits: "the packet equals P1's consensus, and P1 says both families
agreed". That is P1's assertion re-reported, not an independent check, and the
frozen plan asks P3 to establish that packet values match BOTH source families.

The two families write their forms differently, which is what makes the
comparison worth doing:

    theory     `mp.mpf(4) / 5 * phi ** -2`
    platform   `(4/5)phi^-2`      and     `mp.mpf(4) / 5 * phi ** -2`

One grammar covers both. It is deliberately tiny and STRICT: integers, `mp.mpf`
wrappers, parentheses, `phi` with `**` or `^` powers, explicit `*` and `/`, and
juxtaposition as multiplication. Anything else raises, and a raise is a STOP. If
a source changes its syntax, P3 must halt and be looked at rather than quietly
skip a row.
"""

import re

import p3_qphi as Q

_NUM = re.compile(r"-?\d+")
PHI = (0, 1, 1)


def evaluate(text):
    """The stated closed form as a normalized Q(phi) triple. Raises if it is not
    exactly one of the shapes this grammar admits."""
    if not isinstance(text, str) or not text.strip():
        raise ValueError("empty or non-string closed form")
    s = text
    pos = 0

    def skip():
        nonlocal pos
        while pos < len(s) and s[pos].isspace():
            pos += 1

    def starts(tok):
        return s.startswith(tok, pos)

    def factor():
        nonlocal pos
        skip()
        if pos >= len(s):
            raise ValueError("expected a factor, found end of form")
        if starts("mp.mpf("):
            pos += len("mp.mpf(")
            skip()
            m = _NUM.match(s, pos)
            if not m:
                raise ValueError("mp.mpf without an integer")
            pos = m.end()
            skip()
            if pos >= len(s) or s[pos] != ")":
                raise ValueError("unclosed mp.mpf(")
            pos += 1
            return (int(m.group()), 0, 1)
        if s[pos] == "(":
            pos += 1
            v = expr()
            skip()
            if pos >= len(s) or s[pos] != ")":
                raise ValueError("unclosed (")
            pos += 1
            return v
        if starts("phi"):
            pos += 3
            skip()
            if starts("**"):
                pos += 2
            elif pos < len(s) and s[pos] == "^":
                pos += 1
            else:
                return PHI
            skip()
            m = _NUM.match(s, pos)
            if not m:
                raise ValueError("phi power without an integer")
            pos = m.end()
            return Q.power(PHI, int(m.group()))
        m = _NUM.match(s, pos)
        if not m:
            raise ValueError(f"unparseable at {s[pos:pos + 12]!r}")
        pos = m.end()
        return (int(m.group()), 0, 1)

    def more():
        """True when a factor follows with no operator: juxtaposition, which
        these sources use for multiplication, as in `(4/5)phi^-2`."""
        return pos < len(s) and (s[pos] == "(" or s[pos].isdigit()
                                 or starts("phi") or starts("mp.mpf("))

    def expr():
        nonlocal pos
        acc = factor()
        while True:
            skip()
            if pos < len(s) and s[pos] in "*/" and not starts("**"):
                op = s[pos]
                pos += 1
                acc = Q.mul(acc, factor()) if op == "*" else Q.mul(acc, Q.inv(factor()))
            elif more():
                acc = Q.mul(acc, factor())
            else:
                return acc

    value = expr()
    skip()
    if pos != len(s):
        raise ValueError(f"trailing {s[pos:]!r}")
    return value
