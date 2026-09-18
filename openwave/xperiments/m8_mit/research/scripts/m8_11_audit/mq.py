"""Exact multi-quadratic numbers: sum_s c_s sqrt(s), s squarefree positive int, c_s in Q(i).

Zero test is exact (square roots of distinct squarefree integers are linearly independent over Q(i)).
"""
from fractions import Fraction
from math import gcd, isqrt
import mpmath

_sqf_cache = {}


def _factor(n):
    f = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            f[p] = f.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def _sqf(n):
    """n = k^2 * s with s squarefree; return (k, s)."""
    if n in _sqf_cache:
        return _sqf_cache[n]
    k, s = 1, 1
    for p, e in _factor(n).items():
        k *= p ** (e // 2)
        if e % 2:
            s *= p
    _sqf_cache[n] = (k, s)
    return k, s


def _primes_of(s):
    return list(_factor(s).keys())


class MQ:
    __slots__ = ("t",)

    def __init__(self, t=None):
        # t: dict s -> (re, im) Fractions, no zero entries
        self.t = t if t is not None else {}

    # ---- constructors
    @staticmethod
    def rat(q, qi=0):
        q = Fraction(q)
        qi = Fraction(qi)
        if q == 0 and qi == 0:
            return MQ()
        return MQ({1: (q, qi)})

    @staticmethod
    def sqrt(q):
        """sqrt of a nonnegative rational."""
        q = Fraction(q)
        if q < 0:
            raise ValueError("negative")
        if q == 0:
            return MQ()
        n, d = q.numerator, q.denominator
        k, s = _sqf(n * d)
        return MQ({s: (Fraction(k, d), Fraction(0))})

    I = None  # set below

    # ---- arithmetic
    def __add__(self, o):
        if not isinstance(o, MQ):
            o = MQ.rat(o)
        t = dict(self.t)
        for s, (a, b) in o.t.items():
            if s in t:
                c, d = t[s]
                c, d = c + a, d + b
                if c == 0 and d == 0:
                    del t[s]
                else:
                    t[s] = (c, d)
            else:
                t[s] = (a, b)
        return MQ(t)

    __radd__ = __add__

    def __neg__(self):
        return MQ({s: (-a, -b) for s, (a, b) in self.t.items()})

    def __sub__(self, o):
        if not isinstance(o, MQ):
            o = MQ.rat(o)
        return self + (-o)

    def __rsub__(self, o):
        return (-self) + o

    def __mul__(self, o):
        if not isinstance(o, MQ):
            if isinstance(o, (int, Fraction)):
                o = Fraction(o)
                if o == 0:
                    return MQ()
                return MQ({s: (a * o, b * o) for s, (a, b) in self.t.items()})
            o = MQ.rat(o)
        t = {}
        for s, (a, b) in self.t.items():
            for u, (c, d) in o.t.items():
                g = gcd(s, u)
                key = (s // g) * (u // g)
                re = (a * c - b * d) * g
                im = (a * d + b * c) * g
                if key in t:
                    x, y = t[key]
                    t[key] = (x + re, y + im)
                else:
                    t[key] = (re, im)
        return MQ({k: v for k, v in t.items() if v[0] != 0 or v[1] != 0})

    __rmul__ = __mul__

    def conj(self):
        return MQ({s: (a, -b) for s, (a, b) in self.t.items()})

    def re(self):
        return MQ({s: (a, Fraction(0)) for s, (a, b) in self.t.items() if a != 0})

    def im(self):
        return MQ({s: (b, Fraction(0)) for s, (a, b) in self.t.items() if b != 0})

    def is_zero(self):
        return not self.t

    def __eq__(self, o):
        if not isinstance(o, MQ):
            o = MQ.rat(o)
        return (self - o).is_zero()

    def __hash__(self):
        return hash(tuple(sorted(self.t.items())))

    def _galois(self, p):
        return MQ({s: ((-a, -b) if s % p == 0 else (a, b)) for s, (a, b) in self.t.items()})

    def inv(self):
        if self.is_zero():
            raise ZeroDivisionError("MQ zero")
        primes = set()
        for s in self.t:
            primes.update(_primes_of(s))
        if not primes:
            a, b = self.t[1]
            n = a * a + b * b
            return MQ({1: (a / n, -b / n)})
        p = max(primes)
        c = self._galois(p)
        return c * (self * c).inv()

    def __truediv__(self, o):
        if not isinstance(o, MQ):
            o = MQ.rat(o)
        return self * o.inv()

    def __rtruediv__(self, o):
        return MQ.rat(o) * self.inv() if not isinstance(o, MQ) else o * self.inv()

    def __pow__(self, n):
        r = MQ.rat(1)
        for _ in range(n):
            r = r * self
        return r

    # ---- conversions
    def to_mpc(self):
        s = mpmath.mpc(0)
        for k, (a, b) in self.t.items():
            r = mpmath.sqrt(k)
            s += r * mpmath.mpc(mpmath.mpf(a.numerator) / a.denominator, mpmath.mpf(b.numerator) / b.denominator)
        return s

    def to_complex(self):
        return complex(self.to_mpc())

    def is_real(self):
        return all(b == 0 for (a, b) in self.t.values())

    def sign(self):
        """sign of a real MQ (exact decision helped by high-precision evaluation, then certified by zero test)."""
        if self.is_zero():
            return 0
        with mpmath.workdps(60):
            v = self.to_mpc()
            assert abs(v.imag) < mpmath.mpf(10) ** -50
            # nonzero exactly; value is an algebraic number of small height, 60 digits suffice to decide sign
            return 1 if v.real > 0 else -1

    def __repr__(self):
        return self.pretty()

    def pretty(self):
        if not self.t:
            return "0"
        parts = []
        for s in sorted(self.t):
            a, b = self.t[s]
            if b == 0:
                c = str(a)
            elif a == 0:
                c = f"{b}*I"
            else:
                c = f"({a}+{b}*I)"
            parts.append(c if s == 1 else f"{c}*sqrt({s})")
        return " + ".join(parts)


MQ.I = MQ({1: (Fraction(0), Fraction(1))})
ZERO = MQ()
ONE = MQ.rat(1)


def mq(x):
    return x if isinstance(x, MQ) else MQ.rat(x)
