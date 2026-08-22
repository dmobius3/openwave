"""Exact Q(phi) for P3, written INDEPENDENTLY of P2's.

P3's evidentiary value comes entirely from not inheriting P2's expectations, and
that has to include the arithmetic. `qphi_exact.Phi` represents an element as a
pair of Fractions and normalizes only when asked for a triple. This module
represents it as the integer triple `(a, b, c)` meaning `(a + b*phi)/c`
throughout, normalizing after every operation, and derives inverses from the
field norm rather than from a Fraction division.

That is a different representation reached by a different route, so agreement
between the two is evidence rather than tautology. If P3 imported P2's
arithmetic, a defect in it would be invisible to exactly the gate meant to catch
it.

Nothing here is answer-bearing: this is the number system, not a value in it.
`phi^2 = phi + 1`, so `(a1 + b1 phi)(a2 + b2 phi) = a1a2 + b1b2 + (a1b2 + b1a2 +
b1b2) phi`.
"""

from math import gcd


def normalize(a, b, c):
    """`(a, b, c)` with `c > 0` and `gcd(a, b, c) = 1`, per § 4.4's encoding."""
    if c == 0:
        raise ZeroDivisionError("Q(phi) element with zero denominator")
    if c < 0:
        a, b, c = -a, -b, -c
    g = gcd(gcd(abs(a), abs(b)), c)
    if g:
        a, b, c = a // g, b // g, c // g
    return (a, b, c)


def is_normalized(t):
    """Exactly three true ints, `c > 0`, `gcd = 1`.

    `type(x) is int` rather than isinstance: `bool` subclasses `int`, and JSON
    `true` in a numeric position would otherwise pass and then compare equal to 1
    in every arithmetic step after it.
    """
    if not (isinstance(t, list) and len(t) == 3 and all(type(x) is int for x in t)):
        return False
    a, b, c = t
    return c > 0 and gcd(gcd(abs(a), abs(b)), c) == 1


def mul(x, y):
    a1, b1, c1 = x
    a2, b2, c2 = y
    return normalize(a1 * a2 + b1 * b2,
                     a1 * b2 + b1 * a2 + b1 * b2,
                     c1 * c2)


def inv(x):
    """From the field norm. `N(a + b phi) = a^2 + ab - b^2`, and the conjugate of
    `a + b phi` under `phi -> 1 - phi` is `(a + b) - b phi`, so
    `c / (a + b phi) = c((a + b) - b phi) / N`."""
    a, b, c = x
    n = a * a + a * b - b * b
    if n == 0:
        raise ZeroDivisionError("Q(phi) inverse of zero")
    return normalize(c * (a + b), -c * b, n)


def power(x, e):
    if e < 0:
        return power(inv(x), -e)
    r = ONE
    for _ in range(e):
        r = mul(r, x)
    return r


ONE = (1, 0, 1)
