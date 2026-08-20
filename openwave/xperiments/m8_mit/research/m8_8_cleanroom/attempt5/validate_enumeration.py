"""
Pre-implementation validation: enumerate the 120 elements of 2I,
verify the canonical ID sort key and SHA-256 hash.
"""
from fractions import Fraction
import hashlib
import json

# Q(phi) arithmetic: represent a + b*phi where phi = (1+sqrt(5))/2
# phi^2 = phi + 1

class QGold:
    """Element of Q(phi) = Q(golden ratio)."""
    __slots__ = ('a', 'b')

    def __init__(self, a, b=0):
        self.a = Fraction(a)
        self.b = Fraction(b)

    def __repr__(self):
        return f"QGold({self.a}, {self.b})"

    def __eq__(self, other):
        if isinstance(other, (int, Fraction)):
            return self.a == other and self.b == 0
        return self.a == other.a and self.b == other.b

    def __hash__(self):
        return hash((self.a, self.b))

    def __add__(self, other):
        if isinstance(other, (int, Fraction)):
            return QGold(self.a + other, self.b)
        return QGold(self.a + other.a, self.b + other.b)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if isinstance(other, (int, Fraction)):
            return QGold(self.a - other, self.b)
        return QGold(self.a - other.a, self.b - other.b)

    def __rsub__(self, other):
        if isinstance(other, (int, Fraction)):
            return QGold(other - self.a, -self.b)
        return QGold(other.a - self.a, other.b - self.b)

    def __neg__(self):
        return QGold(-self.a, -self.b)

    def __mul__(self, other):
        if isinstance(other, (int, Fraction)):
            return QGold(self.a * other, self.b * other)
        # (a1 + b1*phi)(a2 + b2*phi) = a1*a2 + (a1*b2+a2*b1)*phi + b1*b2*phi^2
        # = a1*a2 + b1*b2 + (a1*b2 + a2*b1 + b1*b2)*phi
        return QGold(self.a * other.a + self.b * other.b,
                     self.a * other.b + self.b * other.a + self.b * other.b)

    def __rmul__(self, other):
        return self.__mul__(other)

    def norm(self):
        """Norm from Q(phi) to Q: N(a+b*phi) = a^2 + a*b - b^2."""
        return self.a * self.a + self.a * self.b - self.b * self.b

    def conjugate(self):
        """Galois conjugate: phi -> 1-phi."""
        return QGold(self.a + self.b, -self.b)

    def __truediv__(self, other):
        if isinstance(other, (int, Fraction)):
            return QGold(self.a / other, self.b / other)
        n = other.norm()
        conj = other.conjugate()
        num = self * conj
        return QGold(num.a / n, num.b / n)

    def is_zero(self):
        return self.a == 0 and self.b == 0


# Quaternion with Q(phi) components
class Quat:
    """Quaternion with components in Q(phi)."""
    __slots__ = ('w', 'x', 'y', 'z')

    def __init__(self, w, x, y, z):
        self.w = w if isinstance(w, QGold) else QGold(w)
        self.x = x if isinstance(x, QGold) else QGold(x)
        self.y = y if isinstance(y, QGold) else QGold(y)
        self.z = z if isinstance(z, QGold) else QGold(z)

    def __repr__(self):
        return f"Quat({self.w}, {self.x}, {self.y}, {self.z})"

    def __eq__(self, other):
        return (self.w == other.w and self.x == other.x and
                self.y == other.y and self.z == other.z)

    def __hash__(self):
        return hash((self.w, self.x, self.y, self.z))

    def __mul__(self, other):
        # (w1 + x1*i + y1*j + z1*k)(w2 + x2*i + y2*j + z2*k)
        w = self.w * other.w - self.x * other.x - self.y * other.y - self.z * other.z
        x = self.w * other.x + self.x * other.w + self.y * other.z - self.z * other.y
        y = self.w * other.y - self.x * other.z + self.y * other.w + self.z * other.x
        z = self.w * other.z + self.x * other.y - self.y * other.x + self.z * other.w
        return Quat(w, x, y, z)

    def __neg__(self):
        return Quat(-self.w, -self.x, -self.y, -self.z)

    def sort_key(self):
        """8-integer sort key per protocol § 4.2.
        Components in quaternion_basis order (1, i, j, k).
        Each component (A + B*phi)/2 gives the pair (A, B) as signed integers.
        """
        key = []
        for comp in [self.w, self.x, self.y, self.z]:
            # comp = a + b*phi = (A + B*phi)/2 where A, B are integers
            A = comp.a * 2
            B = comp.b * 2
            assert A.denominator == 1, f"Non-integer A: {A}"
            assert B.denominator == 1, f"Non-integer B: {B}"
            key.append(int(A))
            key.append(int(B))
        return tuple(key)


def parse_generator(gen_strs):
    """Parse a generator from the group packet format."""
    components = []
    for s in gen_strs:
        # format: "(a + b*phi)/2"
        s = s.strip()
        assert s.startswith('(') and s.endswith(')/2'), f"Unexpected format: {s}"
        inner = s[1:-3]  # remove ( and )/2
        # Parse "a + b*phi" or "a - b*phi"
        parts = inner.replace(' ', '').replace('*phi', 'p')
        # Now parts is like "1+0p" or "-1+1p" or "0+1p"
        a_val = 0
        b_val = 0
        # Split by + or - (keeping the sign)
        tokens = []
        current = ''
        for ch in parts:
            if ch in '+-' and current:
                tokens.append(current)
                current = ch
            else:
                current += ch
        if current:
            tokens.append(current)

        for tok in tokens:
            if 'p' in tok:
                coeff = tok.replace('p', '')
                if coeff == '' or coeff == '+':
                    b_val = 1
                elif coeff == '-':
                    b_val = -1
                else:
                    b_val = int(coeff)
            else:
                a_val = int(tok)

        components.append(QGold(Fraction(a_val, 2), Fraction(b_val, 2)))

    return Quat(*components)


def main():
    # Load group packet
    with open('m8_5a_packet.json', 'r') as f:
        gp = json.load(f)

    # Parse generators
    g1 = parse_generator(gp['generators'][0])
    g2 = parse_generator(gp['generators'][1])

    print(f"Generator 1: w={g1.w}, x={g1.x}, y={g1.y}, z={g1.z}")
    print(f"Generator 2: w={g2.w}, x={g2.x}, y={g2.y}, z={g2.z}")
    print(f"g1 sort key: {g1.sort_key()}")
    print(f"g2 sort key: {g2.sort_key()}")

    # Generate the group by closure under multiplication
    identity = Quat(QGold(1), QGold(0), QGold(0), QGold(0))
    elements = {identity}
    frontier = [identity, g1, g2]
    elements.add(g1)
    elements.add(g2)

    while True:
        new_elements = set()
        for a in list(elements):
            for b in [g1, g2, -g1, -g2]:
                ab = a * b
                if ab not in elements and ab not in new_elements:
                    new_elements.add(ab)
                ba = b * a
                if ba not in elements and ba not in new_elements:
                    new_elements.add(ba)
        if not new_elements:
            break
        elements.update(new_elements)

    print(f"\nGroup order: {len(elements)}")
    assert len(elements) == 120, f"Expected 120 elements, got {len(elements)}"

    # Sort by canonical key
    sorted_elems = sorted(elements, key=lambda q: q.sort_key())

    # Print first few and last few
    print(f"\nRank 0: {sorted_elems[0].sort_key()}")
    print(f"Rank 118: {sorted_elems[118].sort_key()}")
    print(f"Rank 119: {sorted_elems[119].sort_key()}")

    # Verify expected values from protocol
    assert sorted_elems[0].sort_key() == (-2, 0, 0, 0, 0, 0, 0, 0), f"Rank 0 mismatch: {sorted_elems[0].sort_key()}"
    assert sorted_elems[118].sort_key() == (1, 0, 1, 0, 1, 0, 1, 0), f"Rank 118 mismatch: {sorted_elems[118].sort_key()}"
    assert sorted_elems[119].sort_key() == (2, 0, 0, 0, 0, 0, 0, 0), f"Rank 119 mismatch: {sorted_elems[119].sort_key()}"

    # Check identity at rank 119
    assert sorted_elems[119] == identity, "Identity not at rank 119"

    # Compute SHA-256
    key_lists = [list(q.sort_key()) for q in sorted_elems]
    json_str = json.dumps(key_lists, separators=(',', ':'))
    sha = hashlib.sha256(json_str.encode('ascii')).hexdigest()
    print(f"\nEnumeration JSON length: {len(json_str)} bytes")
    print(f"SHA-256: {sha}")

    expected_sha = "27ff780d28d5d854d464ead87e8fc20244fac8334bda9f0600c6ee1b3c89561e"
    expected_len = 2389
    assert len(json_str) == expected_len, f"Length mismatch: {len(json_str)} vs {expected_len}"
    assert sha == expected_sha, f"SHA-256 mismatch:\n  got:      {sha}\n  expected: {expected_sha}"
    print("SHA-256 VERIFIED ✓")

    # Check abstract generators from construction packet
    with open('m8_8_construction_packet.json', 'r') as f:
        cp = json.load(f)

    s_id = cp['abstract_generators']['s']  # 118
    t_id = cp['abstract_generators']['t']  # 80

    s = sorted_elems[s_id]
    t = sorted_elems[t_id]

    print(f"\ns = element {s_id}: sort_key = {s.sort_key()}")
    print(f"t = element {t_id}: sort_key = {t.sort_key()}")

    # Check relators: s^3 (st)^-2 = 1 and t^5 (st)^-2 = 1
    # Equivalently: s^3 = (st)^2 and t^5 = (st)^2
    s3 = s * s * s
    st = s * t
    st2 = st * st
    t5 = t * t * t * t * t

    print(f"\ns^3 sort_key: {s3.sort_key()}")
    print(f"(st)^2 sort_key: {st2.sort_key()}")
    print(f"t^5 sort_key: {t5.sort_key()}")

    assert s3 == st2, "Relator s^3 = (st)^2 FAILED"
    print("Relator s^3 = (st)^2: VERIFIED ✓")

    assert t5 == st2, "Relator t^5 = (st)^2 FAILED"
    print("Relator t^5 = (st)^2: VERIFIED ✓")

    # Check orders
    def element_order(q):
        curr = q
        for n in range(1, 121):
            if curr == identity:
                return n
            curr = curr * q
        return None

    s_order = element_order(s)
    t_order = element_order(t)
    st_order = element_order(st)

    print(f"\norder(s) = {s_order}")
    print(f"order(t) = {t_order}")
    print(f"order(st) = {st_order}")

    # s^3 = -1 => order(s) = 6; t^5 = -1 => order(t) = 10; (st)^2 = -1 => order(st) = 4
    assert s_order == 6, f"Expected order(s) = 6, got {s_order}"
    assert t_order == 10, f"Expected order(t) = 10, got {t_order}"
    assert st_order == 4, f"Expected order(st) = 4, got {st_order}"
    print("Generator orders: VERIFIED ✓")

    # Verify s and t generate the full group
    generated = {identity}
    frontier = [identity]
    while frontier:
        new_frontier = []
        for a in frontier:
            for b in [s, t]:
                for c in [a * b, b * a, a * (Quat(-b.w, -b.x, -b.y, -b.z))]:
                    if c not in generated:
                        generated.add(c)
                        new_frontier.append(c)
        frontier = new_frontier

    print(f"\n<s,t> generates {len(generated)} elements")
    assert len(generated) == 120, "s and t do not generate the full group!"
    print("<s,t> = 2I: VERIFIED ✓")

    # Check element 119 is identity
    print(f"\nElement 119 == identity: {sorted_elems[119] == identity}")

    # Check that s = element 118 = g1 (first group packet generator)
    print(f"s (element 118) == g1: {s == g1}")

    print("\n=== ALL ENUMERATION CHECKS PASSED ===")


if __name__ == '__main__':
    main()
