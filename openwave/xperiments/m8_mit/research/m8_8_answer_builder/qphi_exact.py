"""Exact Q(phi) arithmetic for the P1 data assembly.

Elements are `p + q*phi` with p, q rational and `phi^2 = phi + 1`. Nothing here
is answer-bearing: this is the number system, not any value in it.

The normalized triple `(a, b, c)` means `(a + b*phi)/c` with `c > 0` and
`gcd(a, b, c) = 1`, which is the protocol's frozen value encoding (§ 4.4).
Note that this is NOT the element-ID sort key of § 4.2; that one fixes the
denominator at 2 and does not reduce. The two canonical forms sit four
sentences apart in the protocol and conflating them renames 26 of 120 IDs.
"""

from fractions import Fraction as F
from math import gcd


class Phi:
    """p + q*phi, exact."""

    __slots__ = ("p", "q")

    def __init__(self, p=0, q=0):
        self.p = F(p)
        self.q = F(q)

    def __add__(self, o):
        o = _co(o)
        return Phi(self.p + o.p, self.q + o.q)

    def __radd__(self, o):
        return _co(o) + self

    def __neg__(self):
        return Phi(-self.p, -self.q)

    def __sub__(self, o):
        return self + (-_co(o))

    def __rsub__(self, o):
        return _co(o) + (-self)

    def __mul__(self, o):
        o = _co(o)
        # (p1 + q1 phi)(p2 + q2 phi) = p1p2 + q1q2  +  (p1q2 + q1p2 + q1q2) phi
        return Phi(self.p * o.p + self.q * o.q,
                   self.p * o.q + self.q * o.p + self.q * o.q)

    def __rmul__(self, o):
        return _co(o) * self

    def conj(self):
        """The nontrivial Galois automorphism, phi -> 1 - phi."""
        return Phi(self.p + self.q, -self.q)

    def norm(self):
        """p^2 + pq - q^2, rational. Zero only for the zero element."""
        n = self * self.conj()
        assert n.q == 0
        return n.p

    def inv(self):
        n = self.norm()
        if n == 0:
            raise ZeroDivisionError("Q(phi) inverse of zero")
        c = self.conj()
        return Phi(c.p / n, c.q / n)

    def __truediv__(self, o):
        return self * _co(o).inv()

    def __rtruediv__(self, o):
        return _co(o) * self.inv()

    def __pow__(self, n):
        if n < 0:
            return self.inv() ** (-n)
        r, b = Phi(1), self
        while n:
            if n & 1:
                r = r * b
            b, n = b * b, n >> 1
        return r

    def __eq__(self, o):
        o = _co(o)
        return self.p == o.p and self.q == o.q

    def __hash__(self):
        return hash((self.p, self.q))

    def __repr__(self):
        return f"({self.p} + {self.q}phi)"

    def triple(self):
        """Normalized (a, b, c) with c > 0 and gcd(a, b, c) = 1, per § 4.4."""
        c = self.p.denominator * self.q.denominator // gcd(
            self.p.denominator, self.q.denominator)
        a = int(self.p * c)
        b = int(self.q * c)
        g = gcd(gcd(abs(a), abs(b)), c)
        if g:
            a, b, c = a // g, b // g, c // g
        if c < 0:
            a, b, c = -a, -b, -c
        return (a, b, c)

    def approx(self):
        return float(self.p) + float(self.q) * (1 + 5 ** 0.5) / 2


def _co(o):
    return o if isinstance(o, Phi) else Phi(o, 0)


PHI = Phi(0, 1)
ONE = Phi(1, 0)


class Quat:
    """Quaternion with Q(phi) components, basis order (1, i, j, k)."""

    __slots__ = ("c",)

    def __init__(self, c):
        self.c = tuple(c)

    def __mul__(self, o):
        a0, a1, a2, a3 = self.c
        b0, b1, b2, b3 = o.c
        return Quat((
            a0 * b0 - a1 * b1 - a2 * b2 - a3 * b3,
            a0 * b1 + a1 * b0 + a2 * b3 - a3 * b2,
            a0 * b2 - a1 * b3 + a2 * b0 + a3 * b1,
            a0 * b3 + a1 * b2 - a2 * b1 + a3 * b0,
        ))

    def conj(self):
        a0, a1, a2, a3 = self.c
        return Quat((a0, -a1, -a2, -a3))

    def key(self):
        """§ 4.2 ID sort key: eight signed integers, the numerator pair of each
        component over the FIXED denominator 2, in quaternion_basis order.
        Not reduced. Raises if a component is not in (1/2)Z[phi]."""
        out = []
        for x in self.c:
            for r in (x.p, x.q):
                v = r * 2
                if v.denominator != 1:
                    raise ValueError(f"component not over denominator 2: {r}")
                out.append(int(v))
        return tuple(out)

    def __eq__(self, o):
        return self.c == o.c

    def __hash__(self):
        return hash(self.c)
