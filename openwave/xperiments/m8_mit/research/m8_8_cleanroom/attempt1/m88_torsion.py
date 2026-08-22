#!/usr/bin/env python3
"""M8.8 Clean Room: Reidemeister torsion of S^3/2I from the based chain complex."""

import json
import hashlib
import sys
from fractions import Fraction
from itertools import combinations
from math import comb

# ============================================================
# Section 1: Exact arithmetic over Q(phi) and Q(phi, i)
# ============================================================

class QPhi:
    """Element of Q(phi) where phi = (1+sqrt(5))/2, phi^2 = phi + 1."""
    __slots__ = ('a', 'b')

    def __init__(self, a=0, b=0):
        self.a = a if isinstance(a, Fraction) else Fraction(a)
        self.b = b if isinstance(b, Fraction) else Fraction(b)

    def __add__(self, other):
        if isinstance(other, (int, Fraction)): other = QPhi(other)
        return QPhi(self.a + other.a, self.b + other.b)
    def __radd__(self, other): return self.__add__(other)
    def __sub__(self, other):
        if isinstance(other, (int, Fraction)): other = QPhi(other)
        return QPhi(self.a - other.a, self.b - other.b)
    def __rsub__(self, other): return QPhi(other).__sub__(self)
    def __neg__(self): return QPhi(-self.a, -self.b)
    def __mul__(self, other):
        if isinstance(other, (int, Fraction)): other = QPhi(other)
        if isinstance(other, QPhiI): return QPhiI(self) * other
        # (a+b*phi)(c+d*phi) = ac + bd + (ad+bc+bd)*phi  [phi^2=phi+1]
        return QPhi(self.a*other.a + self.b*other.b,
                     self.a*other.b + self.b*other.a + self.b*other.b)
    def __rmul__(self, other): return self.__mul__(other)
    def __eq__(self, other):
        if isinstance(other, (int, Fraction)): other = QPhi(other)
        if not isinstance(other, QPhi): return NotImplemented
        return self.a == other.a and self.b == other.b
    def __hash__(self): return hash((self.a, self.b))
    def __bool__(self): return self.a != 0 or self.b != 0
    def norm_to_Q(self):
        # N(a+b*phi) = (a+b*phi)(a+b*(1-phi)) = a^2+ab-b^2
        return self.a*self.a + self.a*self.b - self.b*self.b
    def galois_conj(self):
        # sigma(phi) = 1-phi, so sigma(a+b*phi) = a+b-b*phi = (a+b) + (-b)*phi
        return QPhi(self.a + self.b, -self.b)
    def inv(self):
        n = self.norm_to_Q()
        if n == 0: raise ZeroDivisionError("QPhi division by zero")
        return QPhi((self.a + self.b) / n, -self.b / n)
    def __truediv__(self, other):
        if isinstance(other, (int, Fraction)): other = QPhi(other)
        return self * other.inv()
    def __repr__(self):
        return f"QPhi({self.a}, {self.b})"
    def to_triple(self):
        """Normalized triple (a, b, c) for (a + b*phi)/c with c>0, gcd(|a|,|b|,c)=1."""
        from math import gcd
        if self.a == 0 and self.b == 0:
            return (0, 0, 1)
        p, q = self.a, self.b  # value = p + q*phi
        # Express as (num_a + num_b*phi) / den
        dp = p.denominator
        dq = q.denominator
        from math import lcm
        den = lcm(dp, dq)
        num_a = int(p * den)
        num_b = int(q * den)
        g = gcd(gcd(abs(num_a), abs(num_b)), den)
        return (num_a // g, num_b // g, den // g)


class QPhiI:
    """Element of Q(phi, i) where i^2 = -1. Represents re + im*i with re,im in Q(phi)."""
    __slots__ = ('re', 'im')

    def __init__(self, re=None, im=None):
        if re is None: re = QPhi()
        elif isinstance(re, (int, Fraction)): re = QPhi(re)
        elif isinstance(re, QPhiI):
            if im is not None: raise ValueError
            self.re = re.re; self.im = re.im; return
        if im is None: im = QPhi()
        elif isinstance(im, (int, Fraction)): im = QPhi(im)
        self.re = re
        self.im = im

    def __add__(self, other):
        other = self._c(other)
        return QPhiI(self.re + other.re, self.im + other.im)
    def __radd__(self, other): return self.__add__(other)
    def __sub__(self, other):
        other = self._c(other)
        return QPhiI(self.re - other.re, self.im - other.im)
    def __rsub__(self, other): return self._c(other).__sub__(self)
    def __neg__(self): return QPhiI(-self.re, -self.im)
    def __mul__(self, other):
        other = self._c(other)
        return QPhiI(self.re*other.re - self.im*other.im,
                      self.re*other.im + self.im*other.re)
    def __rmul__(self, other): return self._c(other).__mul__(self)
    def __eq__(self, other):
        other = self._c(other)
        return self.re == other.re and self.im == other.im
    def __hash__(self): return hash((self.re, self.im))
    def __bool__(self): return bool(self.re) or bool(self.im)
    def conj(self):
        """Complex conjugate: i -> -i."""
        return QPhiI(self.re, -self.im)
    def abs_sq(self):
        """Return |z|^2 = z * conj(z) in Q(phi)."""
        return self.re * self.re + self.im * self.im
    def inv(self):
        n = self.abs_sq()
        return QPhiI(self.re / n, QPhi() - self.im / n)
    def __truediv__(self, other):
        other = self._c(other)
        return self * other.inv()
    def is_zero(self):
        return not bool(self.re) and not bool(self.im)
    def galois_conj(self):
        """Galois conjugate: phi -> 1-phi (does NOT touch i)."""
        return QPhiI(self.re.galois_conj(), self.im.galois_conj())
    def _c(self, other):
        if isinstance(other, QPhiI): return other
        if isinstance(other, QPhi): return QPhiI(other)
        if isinstance(other, (int, Fraction)): return QPhiI(QPhi(other))
        raise TypeError(f"Cannot coerce {type(other)} to QPhiI")
    def __repr__(self):
        return f"QPhiI({self.re}, {self.im})"


ZERO = QPhiI(QPhi(0), QPhi(0))
ONE = QPhiI(QPhi(1), QPhi(0))
PHI = QPhi(0, 1)  # phi as a QPhi element
IMAG = QPhiI(QPhi(0), QPhi(1))  # the complex i

# ============================================================
# Section 2: Matrix operations over QPhiI
# ============================================================

def mat_zeros(r, c):
    return [[QPhiI() for _ in range(c)] for _ in range(r)]

def mat_id(n):
    M = mat_zeros(n, n)
    for i in range(n): M[i][i] = ONE
    return M

def mat_scale(M, s):
    s = QPhiI(s) if not isinstance(s, QPhiI) else s
    return [[M[i][j] * s for j in range(len(M[0]))] for i in range(len(M))]

def mat_add(A, B):
    return [[A[i][j] + B[i][j] for j in range(len(A[0]))] for i in range(len(A))]

def mat_sub(A, B):
    return [[A[i][j] - B[i][j] for j in range(len(A[0]))] for i in range(len(A))]

def mat_mul(A, B):
    rA, cA, cB = len(A), len(A[0]), len(B[0])
    C = mat_zeros(rA, cB)
    for i in range(rA):
        for k in range(cA):
            if A[i][k].is_zero(): continue
            for j in range(cB):
                C[i][j] = C[i][j] + A[i][k] * B[k][j]
    return C

def mat_transpose(M):
    r, c = len(M), len(M[0])
    return [[M[j][i] for j in range(r)] for i in range(c)]

def mat_conj_transpose(M):
    """Conjugate transpose (dagger)."""
    r, c = len(M), len(M[0])
    return [[M[j][i].conj() for j in range(r)] for i in range(c)]

def mat_trace(M):
    return sum((M[i][i] for i in range(len(M))), QPhiI())

def mat_submatrix(M, rows, cols):
    return [[M[r][c] for c in cols] for r in rows]

def mat_det(M):
    """Determinant via Gaussian elimination over QPhiI."""
    n = len(M)
    A = [[M[i][j] for j in range(n)] for i in range(n)]
    det = QPhiI(QPhi(1))
    for col in range(n):
        pivot = None
        for row in range(col, n):
            if not A[row][col].is_zero():
                pivot = row
                break
        if pivot is None:
            return QPhiI()
        if pivot != col:
            A[col], A[pivot] = A[pivot], A[col]
            det = -det
        det = det * A[col][col]
        inv_p = A[col][col].inv()
        for row in range(col + 1, n):
            if A[row][col].is_zero(): continue
            factor = A[row][col] * inv_p
            for j in range(col + 1, n):
                A[row][j] = A[row][j] - factor * A[col][j]
            A[row][col] = QPhiI()
    return det

def mat_rank(M):
    """Rank via Gaussian elimination."""
    r, c = len(M), len(M[0])
    A = [[M[i][j] for j in range(c)] for i in range(r)]
    rank = 0
    for col in range(c):
        pivot = None
        for row in range(rank, r):
            if not A[row][col].is_zero():
                pivot = row
                break
        if pivot is None: continue
        A[rank], A[pivot] = A[pivot], A[rank]
        inv_p = A[rank][col].inv()
        for row in range(rank + 1, r):
            if A[row][col].is_zero(): continue
            factor = A[row][col] * inv_p
            for j in range(col, c):
                A[row][j] = A[row][j] - factor * A[rank][j]
        rank += 1
    return rank

def mat_is_pos_def_hermitian(H):
    """Check positive definiteness via Cholesky (exact)."""
    n = len(H)
    L = mat_zeros(n, n)
    for j in range(n):
        s = QPhi()
        for k in range(j):
            s = s + L[j][k].abs_sq()
        diag = H[j][j].re - s
        if diag.norm_to_Q() <= 0 or diag.b != 0:
            if diag == QPhi(0): return False
            # Check sign: need diag > 0. diag = a + b*phi. phi ~ 1.618.
            # diag > 0 iff a + b*1.618... > 0
            # Use: diag > 0 iff (a+b > 0 and a > 0) or ... just check numerically
            approx = float(diag.a) + float(diag.b) * 1.6180339887
            if approx <= 0: return False
        L[j][j] = QPhiI(diag)  # store diag^2, not diag itself; just checking positivity
        for i in range(j + 1, n):
            s2 = QPhiI()
            for k in range(j):
                s2 = s2 + L[i][k] * L[j][k].conj()
            L[i][j] = (H[i][j] - s2)  # / L[j][j] but we only need sign of diag
    return True

# ============================================================
# Section 3: Quaternion arithmetic over Q(phi)
# ============================================================

def quat_mul(q1, q2):
    """Multiply two quaternions (w,x,y,z) over Q(phi)."""
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    return (
        w1*w2 - x1*x2 - y1*y2 - z1*z2,
        w1*x2 + x1*w2 + y1*z2 - z1*y2,
        w1*y2 - x1*z2 + y1*w2 + z1*x2,
        w1*z2 + x1*y2 - y1*x2 + z1*w2,
    )

def quat_inv(q):
    """Inverse of a unit quaternion = conjugate."""
    w, x, y, z = q
    return (w, -x, -y, -z)

def quat_to_8int(q):
    """Convert quaternion (w,x,y,z) in Q(phi) to the 8-integer key.
    Each component is (A + B*phi)/2; key is [A_w, B_w, A_x, B_x, A_y, B_y, A_z, B_z]."""
    key = []
    for comp in q:
        # comp = a + b*phi. Express as (A + B*phi)/2: A = 2a, B = 2b
        A = int(comp.a * 2)
        B = int(comp.b * 2)
        if Fraction(A, 2) != comp.a or Fraction(B, 2) != comp.b:
            raise ValueError(f"Component {comp} not in (A+B*phi)/2 form")
        key.extend([A, B])
    return key

# ============================================================
# Section 4: Group construction and canonical enumeration
# ============================================================

def build_group(gen_quats):
    """Build the 120-element group 2I by closure under multiplication."""
    elements = {}
    def add(q):
        key = tuple(quat_to_8int(q))
        if key not in elements:
            elements[key] = q
            return True
        return False

    for g in gen_quats:
        add(g)
    # Also add inverses and identity
    identity = (QPhi(1), QPhi(0), QPhi(0), QPhi(0))
    add(identity)
    neg_id = (QPhi(-1), QPhi(0), QPhi(0), QPhi(0))
    add(neg_id)

    changed = True
    while changed:
        changed = False
        keys = list(elements.keys())
        for k1 in keys:
            for k2 in keys:
                p = quat_mul(elements[k1], elements[k2])
                if add(p):
                    changed = True
        if len(elements) > 130:
            raise RuntimeError(f"Group too large: {len(elements)} elements")

    if len(elements) != 120:
        raise RuntimeError(f"Expected 120 elements, got {len(elements)}")
    return elements

def canonical_enumeration(elements_dict):
    """Sort elements by 8-integer key, return ordered list and verify SHA-256."""
    keys_sorted = sorted(elements_dict.keys())
    ordered = []
    for k in keys_sorted:
        ordered.append(elements_dict[k])

    # Build the JSON array for SHA-256
    arr = [list(k) for k in keys_sorted]
    json_str = json.dumps(arr, separators=(',', ':'))
    sha = hashlib.sha256(json_str.encode('ascii')).hexdigest()

    expected_sha = "27ff780d28d5d854d464ead87e8fc20244fac8334bda9f0600c6ee1b3c89561e"
    print(f"Enumeration SHA-256: {sha}")
    print(f"Expected:           {expected_sha}")
    print(f"Match: {sha == expected_sha}")

    # Verify specific ranks
    assert list(keys_sorted[0]) == [-2,0,0,0,0,0,0,0], f"Rank 0: {keys_sorted[0]}"
    assert list(keys_sorted[118]) == [1,0,1,0,1,0,1,0], f"Rank 118: {keys_sorted[118]}"
    assert list(keys_sorted[119]) == [2,0,0,0,0,0,0,0], f"Rank 119: {keys_sorted[119]}"

    return ordered, keys_sorted, sha == expected_sha

def build_mult_table(ordered, keys_sorted):
    """Build multiplication table: mult[i][j] = ID of product of element i and element j."""
    key_to_id = {k: i for i, k in enumerate(keys_sorted)}
    n = len(ordered)
    mult = [[0]*n for _ in range(n)]
    inv_table = [0]*n
    for i in range(n):
        for j in range(n):
            p = quat_mul(ordered[i], ordered[j])
            pk = tuple(quat_to_8int(p))
            mult[i][j] = key_to_id[pk]
        # Inverse
        qi = quat_inv(ordered[i])
        qk = tuple(quat_to_8int(qi))
        inv_table[i] = key_to_id[qk]
    return mult, inv_table

# ============================================================
# Section 5: Representation construction
# ============================================================

def quat_to_su2(q):
    """Map quaternion (w,x,y,z) to 2x2 matrix over Q(phi,i).
    M(q) = [[w+z*i, x+y*i], [-x+y*i, w-z*i]]  (i = complex imaginary unit)."""
    w, x, y, z = q
    return [
        [QPhiI(w, z), QPhiI(x, y)],
        [QPhiI(QPhi()-x, y), QPhiI(w, QPhi()-z)]
    ]

def sym_power_matrix(M2x2, n):
    """Compute Sym^n of a 2x2 matrix. Returns (n+1)x(n+1) matrix.
    Basis: {e1^{n-k} e2^k : k = 0, ..., n}.
    Entry (j,k) = sum over valid r of C(n-k,r)*C(k,j-r)*a^{n-k-r}*c^r*b^{k-j+r}*d^{j-r}
    where M = [[a,b],[c,d]]."""
    if n == 0:
        return [[ONE]]
    if n == 1:
        return [[M2x2[i][j] for j in range(2)] for i in range(2)]
    a, b = M2x2[0][0], M2x2[0][1]
    c, d = M2x2[1][0], M2x2[1][1]
    dim = n + 1
    S = mat_zeros(dim, dim)
    # Precompute powers
    pow_a = [ONE]
    pow_b = [ONE]
    pow_c = [ONE]
    pow_d = [ONE]
    for _ in range(n):
        pow_a.append(pow_a[-1] * a)
        pow_b.append(pow_b[-1] * b)
        pow_c.append(pow_c[-1] * c)
        pow_d.append(pow_d[-1] * d)
    for j in range(dim):
        for k in range(dim):
            val = QPhiI()
            r_lo = max(0, j - k)
            r_hi = min(j, n - k)
            for r in range(r_lo, r_hi + 1):
                coeff = comb(n - k, r) * comb(k, j - r)
                term = pow_a[n - k - r] * pow_c[r] * pow_b[k - j + r] * pow_d[j - r]
                val = val + term * coeff
            S[j][k] = val
    return S

def build_all_irreps(ordered, inv_table):
    """Build all 9 irreducible representations of 2I.
    Returns list of (label, dim, rep_matrices) where rep_matrices[element_id] = matrix."""
    n_elts = len(ordered)

    # Natural 2-dim rep V1 via SU(2) embedding
    v1_mats = [quat_to_su2(q) for q in ordered]

    # Galois conjugate V7: apply phi -> 1-phi to all QPhi components
    def galois_conj_mat(M):
        return [[e.galois_conj() for e in row] for row in M]
    v7_mats = [galois_conj_mat(m) for m in v1_mats]

    # Symmetric powers of V1
    print("  Building Sym^2(V1)...")
    v2_mats = [sym_power_matrix(m, 2) for m in v1_mats]
    print("  Building Sym^3(V1)...")
    v3_mats = [sym_power_matrix(m, 3) for m in v1_mats]
    print("  Building Sym^4(V1)...")
    v4_mats = [sym_power_matrix(m, 4) for m in v1_mats]
    print("  Building Sym^5(V1)...")
    v5_mats = [sym_power_matrix(m, 5) for m in v1_mats]

    # Galois conjugates of symmetric powers
    v8_mats = [galois_conj_mat(m) for m in v2_mats]  # Sym^2(V7) = sigma(Sym^2(V1))
    v6_mats = [galois_conj_mat(m) for m in v3_mats]  # Sym^3(V7) = sigma(Sym^3(V1))

    # Trivial rep
    v0_mats = [[[ONE]] for _ in range(n_elts)]

    irreps = [
        ("V0", 1, v0_mats),
        ("V1", 2, v1_mats),
        ("V7", 2, v7_mats),
        ("V2", 3, v2_mats),
        ("V8", 3, v8_mats),
        ("V3", 4, v3_mats),
        ("V6", 4, v6_mats),
        ("V4", 5, v4_mats),
        ("V5", 6, v5_mats),
    ]
    return irreps

# ============================================================
# Section 6: Boundary map parsing and evaluation
# ============================================================

def parse_boundary_maps(packet):
    """Parse d1, d2, d3 from the construction packet.
    Returns each as a list of lists of group-ring elements.
    A group-ring element is a list of (coeff, element_id) pairs."""
    d1 = packet["boundary_maps"]["d1"]
    d2 = packet["boundary_maps"]["d2"]
    d3 = packet["boundary_maps"]["d3"]
    return d1, d2, d3

def eval_group_ring_element(gr_elem, rep_mats, dim):
    """Evaluate a group-ring element at a representation.
    gr_elem: list of [coeff, element_id] pairs.
    rep_mats: list of matrices indexed by element_id.
    Returns a dim x dim matrix."""
    result = mat_zeros(dim, dim)
    for coeff, eid in gr_elem:
        M = rep_mats[eid]
        for i in range(dim):
            for j in range(dim):
                result[i][j] = result[i][j] + M[i][j] * coeff
    return result

def eval_boundary_map(d_matrix, rep_mats, dim):
    """Evaluate a boundary map (matrix of group-ring elements) at a representation.
    d_matrix[i][j] is a group-ring element (list of [coeff, eid] pairs).
    Returns a (rows*dim) x (cols*dim) block matrix."""
    rows = len(d_matrix)
    cols = len(d_matrix[0])
    big = mat_zeros(rows * dim, cols * dim)
    for i in range(rows):
        for j in range(cols):
            block = eval_group_ring_element(d_matrix[i][j], rep_mats, dim)
            for bi in range(dim):
                for bj in range(dim):
                    big[i*dim + bi][j*dim + bj] = block[bi][bj]
    return big

# ============================================================
# Section 7: Torsion computation
# ============================================================

def find_nonsingular_rows(M, nrows, ncols, count):
    """Find 'count' rows forming a nonsingular submatrix via pivot selection."""
    A = [[M[i][j] for j in range(ncols)] for i in range(nrows)]
    used = [False] * nrows
    pivot_rows = []
    for col in range(ncols):
        best = None
        for row in range(nrows):
            if not used[row] and not A[row][col].is_zero():
                best = row
                break
        if best is None: continue
        used[best] = True
        pivot_rows.append(best)
        inv_p = A[best][col].inv()
        for row in range(nrows):
            if row != best and not A[row][col].is_zero():
                factor = A[row][col] * inv_p
                for j in range(col, ncols):
                    A[row][j] = A[row][j] - factor * A[best][j]
        if len(pivot_rows) == count: break
    return sorted(pivot_rows) if len(pivot_rows) == count else None

def find_nonsingular_cols(M, nrows, ncols, count):
    """Find 'count' columns forming a nonsingular submatrix via pivot selection."""
    MT = mat_transpose(M)
    return find_nonsingular_rows(MT, ncols, nrows, count)

def compute_torsion_sq(D1, D2, D3, d, return_factors=False):
    """Compute T^2 = |tau|^2 for an acyclic 4-term complex.
    D1: 2d x d, D2: 2d x 2d, D3: d x 2d.
    Returns T^2 as a QPhi value.
    If return_factors=True, also returns (Ip, I, J, Jp, det_D1, det_D2, det_D3)."""

    all_idx = list(range(2 * d))

    Ip = find_nonsingular_rows(D1, 2*d, d, d)
    if Ip is None:
        raise ValueError("Cannot find nonsingular rows in D1")
    I = sorted(set(all_idx) - set(Ip))

    J = find_nonsingular_cols(D3, d, 2*d, d)
    if J is None:
        raise ValueError("Cannot find nonsingular cols in D3")
    Jp = sorted(set(all_idx) - set(J))

    det_D1 = mat_det(mat_submatrix(D1, Ip, list(range(d))))
    det_D3 = mat_det(mat_submatrix(D3, list(range(d)), J))
    det_D2 = mat_det(mat_submatrix(D2, Jp, I))

    if det_D2.is_zero():
        raise ValueError("D2 submatrix is singular")

    abs2_D1 = det_D1.abs_sq()
    abs2_D2 = det_D2.abs_sq()
    abs2_D3 = det_D3.abs_sq()

    T2 = abs2_D2 / (abs2_D1 * abs2_D3)
    if return_factors:
        return T2, (Ip, I, J, Jp, det_D1, det_D2, det_D3)
    return T2

# ============================================================
# Section 8: Row signature computation
# ============================================================

def compute_row_signature(rep_mats, dim, s_id, t_id, mult_table):
    """Compute the row signature: (dim, chi(s), chi(t), chi(st))."""
    chi_s = mat_trace(rep_mats[s_id])
    chi_t = mat_trace(rep_mats[t_id])
    st_id = mult_table[s_id][t_id]
    chi_st = mat_trace(rep_mats[st_id])
    # Convert to QPhi triples
    def to_triple(z):
        # Character should be real (in Q(phi)), verify im = 0
        assert z.im == QPhi(0), f"Character has nonzero imaginary part: {z}"
        return z.re.to_triple()
    return (dim, to_triple(chi_s), to_triple(chi_t), to_triple(chi_st))

# ============================================================
# Section 9: Main computation
# ============================================================

def main():
    # Load packets
    print("=" * 60)
    print("M8.8 Clean Room: Reidemeister Torsion Computation")
    print("=" * 60)

    with open("m8_5a_packet.json") as f:
        group_packet = json.load(f)
    with open("m8_8_construction_packet.json") as f:
        constr_packet = json.load(f)

    # Verify packet hashes
    with open("m8_5a_packet.json", "rb") as f:
        gp_hash = hashlib.sha256(f.read()).hexdigest()
    with open("m8_8_construction_packet.json", "rb") as f:
        cp_hash = hashlib.sha256(f.read()).hexdigest()
    print(f"\nGroup packet SHA-256:        {gp_hash}")
    print(f"Construction packet SHA-256: {cp_hash}")
    print(f"Construction packet declares: {constr_packet['group_packet_sha256']}")
    print(f"Match: {gp_hash == constr_packet['group_packet_sha256']}")

    # Parse generators
    print("\n--- Building group 2I ---")
    gens_raw = group_packet["generators"]
    phi = QPhi(0, 1)
    def parse_qphi(s):
        # Parse "(a + b*phi)/2"
        import re
        m = re.match(r'\((-?\d+)\s*\+\s*(-?\d+)\*phi\)/2', s)
        a, b = int(m.group(1)), int(m.group(2))
        return QPhi(Fraction(a, 2), Fraction(b, 2))

    gen_quats = []
    for g in gens_raw:
        q = tuple(parse_qphi(s) for s in g)
        gen_quats.append(q)

    print(f"Generator 0: {[quat_to_8int(gen_quats[0])]}")
    print(f"Generator 1: {[quat_to_8int(gen_quats[1])]}")

    # Build group
    elements_dict = build_group(gen_quats)
    print(f"Group size: {len(elements_dict)}")

    # Canonical enumeration
    print("\n--- Canonical enumeration ---")
    ordered, keys_sorted, enum_ok = canonical_enumeration(elements_dict)
    if not enum_ok:
        print("FATAL: Enumeration SHA-256 mismatch!")
        sys.exit(1)
    print("Enumeration gate E1: PASS")

    # Build multiplication table
    print("\n--- Building multiplication table ---")
    mult_table, inv_table = build_mult_table(ordered, keys_sorted)
    print("Done.")

    # Verify abstract generators
    s_id = constr_packet["abstract_generators"]["s"]
    t_id = constr_packet["abstract_generators"]["t"]
    print(f"\nAbstract generators: s = element {s_id}, t = element {t_id}")
    print(f"s 8-int key: {list(keys_sorted[s_id])}")
    print(f"t 8-int key: {list(keys_sorted[t_id])}")

    # Check relators: s^3 (st)^-2 = 1 and t^5 (st)^-2 = 1
    identity_id = 119  # element [2,0,0,0,0,0,0,0]
    st_id = mult_table[s_id][t_id]

    def group_power(eid, n):
        if n == 0: return identity_id
        result = eid
        for _ in range(n - 1):
            result = mult_table[result][eid]
        return result

    s3 = group_power(s_id, 3)
    t5 = group_power(t_id, 5)
    st2 = group_power(st_id, 2)
    print(f"s^3 = element {s3} (should be central: -1 = element 0)")
    print(f"t^5 = element {t5}")
    print(f"(st)^2 = element {st2}")
    print(f"s^3 = (st)^2: {s3 == st2}")
    print(f"t^5 = (st)^2: {t5 == st2}")

    # Check orders
    s_order = 1
    eid = s_id
    while eid != identity_id:
        eid = mult_table[eid][s_id]
        s_order += 1
        if s_order > 120: break
    t_order = 1
    eid = t_id
    while eid != identity_id:
        eid = mult_table[eid][t_id]
        t_order += 1
        if t_order > 120: break
    print(f"order(s) = {s_order}, order(t) = {t_order}")

    # Parse boundary maps
    print("\n--- Parsing boundary maps ---")
    d1, d2, d3 = parse_boundary_maps(constr_packet)
    print(f"d1: {len(d1)} x {len(d1[0])} over Z[2I]")
    print(f"d2: {len(d2)} x {len(d2[0])} over Z[2I]")
    print(f"d3: {len(d3)} x {len(d3[0])} over Z[2I]")

    # Model gate M6: check eps(d1) = 0 and generator correspondence
    print("\n--- Model gate M6: augmentation and generator correspondence ---")
    def augment(gr_elem):
        return sum(c for c, _ in gr_elem)
    for i, row in enumerate(d1):
        for j, entry in enumerate(row):
            aug = augment(entry)
            print(f"  eps(d1[{i}][{j}]) = {aug}", end="")
            if aug != 0:
                print(" FAIL")
            else:
                print(" OK")
    # d1 should be [[s-1], [t-1]]
    assert d1[0][0] == [[1, s_id], [-1, identity_id]], f"d1[0][0] mismatch: {d1[0][0]}"
    assert d1[1][0] == [[1, t_id], [-1, identity_id]], f"d1[1][0] mismatch: {d1[1][0]}"
    print("  Generator correspondence: PASS")

    # Model gate M1: d_compose check (d3.d2 = 0 and d2.d1 = 0 over Z[2I])
    print("\n--- Model gate M1: d3*d2 = 0 and d2*d1 = 0 over Z[2I] ---")
    def gr_mul(a, b, mt):
        """Multiply two group-ring elements using multiplication table."""
        result = {}
        for ca, ea in a:
            for cb, eb in b:
                eid = mt[ea][eb]
                result[eid] = result.get(eid, 0) + ca * cb
        return [[c, e] for e, c in result.items() if c != 0]

    def gr_is_zero(elem):
        return len(elem) == 0

    # d3 (1x2) * d2 (2x2) = 1x2 matrix
    print("  Checking d3*d2:")
    for j in range(2):
        prod = []
        for k in range(2):
            p = gr_mul(d3[0][k], d2[k][j], mult_table)
            # Add to prod
            combined = {}
            for c, e in prod:
                combined[e] = combined.get(e, 0) + c
            for c, e in p:
                combined[e] = combined.get(e, 0) + c
            prod = [[c, e] for e, c in combined.items() if c != 0]
        status = "PASS" if gr_is_zero(prod) else f"FAIL: {prod}"
        print(f"    (d3*d2)[0][{j}] = 0: {status}")

    # d2 (2x2) * d1 (2x1) = 2x1 matrix
    print("  Checking d2*d1:")
    for i in range(2):
        prod = []
        for k in range(2):
            p = gr_mul(d2[i][k], d1[k][0], mult_table)
            combined = {}
            for c, e in prod:
                combined[e] = combined.get(e, 0) + c
            for c, e in p:
                combined[e] = combined.get(e, 0) + c
            prod = [[c, e] for e, c in combined.items() if c != 0]
        status = "PASS" if gr_is_zero(prod) else f"FAIL: {prod}"
        print(f"    (d2*d1)[{i}][0] = 0: {status}")

    # Build representations
    print("\n--- Building irreducible representations ---")
    irreps = build_all_irreps(ordered, inv_table)
    print(f"Built {len(irreps)} irreps.")

    # Verify representations are homomorphisms (spot check)
    print("\n--- Verifying representations (spot check) ---")
    for label, dim, mats in irreps:
        # Check rho(s)*rho(t) = rho(st)
        prod = mat_mul(mats[s_id], mats[t_id])
        expected = mats[st_id]
        ok = all(prod[i][j] == expected[i][j] for i in range(dim) for j in range(dim))
        print(f"  {label} (dim {dim}): rho(s)*rho(t) = rho(st): {'PASS' if ok else 'FAIL'}")

    # Compute row signatures
    print("\n--- Computing row signatures ---")
    signatures = []
    for label, dim, mats in irreps:
        sig = compute_row_signature(mats, dim, s_id, t_id, mult_table)
        signatures.append(sig)
        print(f"  {label}: dim={sig[0]}, chi(s)={sig[1]}, chi(t)={sig[2]}, chi(st)={sig[3]}")

    # Per-irrep computation: evaluate boundaries, check acyclicity, compute torsion
    print("\n" + "=" * 60)
    print("Per-irrep torsion computation")
    print("=" * 60)

    results = []
    for idx, (label, dim, mats) in enumerate(irreps):
        print(f"\n--- {label} (dim {dim}) ---")

        # Evaluate boundary maps
        print("  Evaluating boundary maps...")
        D1 = eval_boundary_map(d1, mats, dim)
        D2 = eval_boundary_map(d2, mats, dim)
        D3 = eval_boundary_map(d3, mats, dim)
        print(f"  D1: {len(D1)}x{len(D1[0])}, D2: {len(D2)}x{len(D2[0])}, D3: {len(D3)}x{len(D3[0])}")

        # Gate D1: check D3*D2 = 0 and D2*D1 = 0
        print("  Checking D3*D2 = 0...", end=" ")
        prod32 = mat_mul(D3, D2)
        ok32 = all(prod32[i][j].is_zero() for i in range(len(prod32)) for j in range(len(prod32[0])))
        print("PASS" if ok32 else "FAIL")

        print("  Checking D2*D1 = 0...", end=" ")
        prod21 = mat_mul(D2, D1)
        ok21 = all(prod21[i][j].is_zero() for i in range(len(prod21)) for j in range(len(prod21[0])))
        print("PASS" if ok21 else "FAIL")

        # Check ranks
        r3 = mat_rank(D3)
        r2 = mat_rank(D2)
        r1 = mat_rank(D1)
        print(f"  Ranks: D3={r3}, D2={r2}, D1={r1}")

        acyclic = (r3 == dim) and (r2 == dim) and (r1 == dim) and (r3 + r2 == 2*dim)
        if label == "V0":
            # Trivial rep: expected non-acyclic
            print(f"  Acyclicity: {'non-acyclic (expected for R0)' if not acyclic else 'UNEXPECTED: acyclic'}")
            results.append((label, dim, signatures[idx], None, "non-acyclic",
                            D1, D2, D3, None, None, None, None, None, None, None))
            continue

        if not acyclic:
            print(f"  Acyclicity: FAIL (nontrivial irrep should be acyclic)")
            results.append((label, dim, signatures[idx], None, "non-acyclic-unexpected",
                            D1, D2, D3, None, None, None, None, None, None, None))
            continue

        print(f"  Acyclicity: PASS")

        # Compute torsion
        print("  Computing torsion T^2...")
        try:
            T2, factors = compute_torsion_sq(D1, D2, D3, dim, return_factors=True)
            Ip, I, J, Jp, det_D1, det_D2, det_D3 = factors
            triple = T2.to_triple()
            print(f"  T^2 = {T2} = ({triple[0]} + {triple[1]}*phi) / {triple[2]}")
            results.append((label, dim, signatures[idx], triple, "acyclic",
                            D1, D2, D3, Ip, I, J, Jp, det_D1, det_D2, det_D3))
        except Exception as e:
            print(f"  FAILED: {e}")
            results.append((label, dim, signatures[idx], None, f"error: {e}",
                            None, None, None, None, None, None, None, None, None, None))

    # Print summary
    print("\n" + "=" * 60)
    print("RESULTS SUMMARY")
    print("=" * 60)
    for r in results:
        label, dim, sig, T2_triple, status = r[0], r[1], r[2], r[3], r[4]
        if T2_triple is not None:
            a, b, c = T2_triple
            print(f"  {label} (dim {dim}): T^2 = ({a} + {b}*phi) / {c}   [{status}]")
        else:
            print(f"  {label} (dim {dim}): {status}")

    # Check Galois consistency (gate D5)
    print("\n--- Gate D5: Galois consistency ---")
    galois_pairs = [("V1", "V7"), ("V2", "V8"), ("V3", "V6")]
    result_map = {r[0]: r for r in results}
    for l1, l2 in galois_pairs:
        r1 = result_map[l1]
        r2 = result_map[l2]
        if r1[3] is None or r2[3] is None:
            print(f"  {l1}/{l2}: skipped (missing data)")
            continue
        a1, b1, c1 = r1[3]
        a2, b2, c2 = r2[3]
        ga, gb, gc = a1 + b1, -b1, c1
        from math import gcd
        g = gcd(gcd(abs(ga), abs(gb)), gc)
        ga, gb, gc = ga // g, gb // g, gc // g
        match = (ga == a2 and gb == b2 and gc == c2)
        print(f"  {l1} -> sigma -> ({ga} + {gb}*phi)/{gc}, {l2} = ({a2} + {b2}*phi)/{c2}: {'PASS' if match else 'FAIL'}")

    for lbl in ["V4", "V5"]:
        r = result_map[lbl]
        if r[3] is None: continue
        a, b, c = r[3]
        self_conj = (b == 0)
        print(f"  {lbl} self-conjugate (b=0): {'PASS' if self_conj else f'FAIL (b={b})'}")

    # Build derivation artifacts
    print("\n--- Building derivation artifacts ---")
    derivation_artifacts = {}
    for r in results:
        label = r[0]
        D1_ev, D2_ev, D3_ev = r[5], r[6], r[7]
        if D1_ev is None:
            continue

        def mat_to_str(M):
            rows = []
            for row in M:
                rows.append("[" + ", ".join(str(x) for x in row) + "]")
            return "[" + ", ".join(rows) + "]"

        artifact = {"label": label}
        D1_str = mat_to_str(D1_ev)
        D2_str = mat_to_str(D2_ev)
        D3_str = mat_to_str(D3_ev)
        artifact["D1_sha256"] = hashlib.sha256(D1_str.encode()).hexdigest()
        artifact["D2_sha256"] = hashlib.sha256(D2_str.encode()).hexdigest()
        artifact["D3_sha256"] = hashlib.sha256(D3_str.encode()).hexdigest()

        if r[3] is not None:  # acyclic
            Ip, I_set, J, Jp = r[8], r[9], r[10], r[11]
            det_D1, det_D2, det_D3 = r[12], r[13], r[14]
            artifact["index_sets"] = {
                "Ip": list(Ip), "I": list(I_set), "J": list(J), "Jp": list(Jp)
            }
            det_str = f"det_D1={det_D1}, det_D2={det_D2}, det_D3={det_D3}"
            artifact["det_factors_sha256"] = hashlib.sha256(det_str.encode()).hexdigest()
            artifact["det_D1_Ip"] = str(det_D1)
            artifact["det_D2_JpI"] = str(det_D2)
            artifact["det_D3_J"] = str(det_D3)

        derivation_artifacts[label] = artifact
        print(f"  {label}: D1={artifact['D1_sha256'][:16]}... D2={artifact['D2_sha256'][:16]}... D3={artifact['D3_sha256'][:16]}...")

    # Load gate results
    gate_results = {}
    try:
        with open("gate_results.json") as f:
            gate_results = json.load(f)
    except FileNotFoundError:
        print("  WARNING: gate_results.json not found")

    # Compute manifest SHA-256
    with open("METHOD_AND_GATE_MANIFEST.md", "rb") as f:
        manifest_sha256 = hashlib.sha256(f.read()).hexdigest()

    # Build output
    print("\n--- Building RAW_OUTPUT.json ---")
    output_rows = []
    for r in results:
        label, dim, sig, T2_triple, status = r[0], r[1], r[2], r[3], r[4]
        row = {
            "dimension": dim,
            "row_signature": {
                "dimension": sig[0],
                "chi_s": list(sig[1]),
                "chi_t": list(sig[2]),
                "chi_st": list(sig[3]),
            },
        }
        if status == "acyclic" and T2_triple is not None:
            row["acyclic"] = True
            row["T2_target"] = list(T2_triple)
        elif status == "non-acyclic":
            row["acyclic"] = False
        else:
            row["acyclic"] = False
            row["error"] = status
        output_rows.append(row)

    raw_output = {
        "schema_version": "m8_8-raw-output-v1",
        "group_packet_sha256": gp_hash,
        "construction_packet_sha256": cp_hash,
        "manifest_sha256": manifest_sha256,
        "rows": output_rows,
        "derivation_artifacts": derivation_artifacts,
        "gate_results": gate_results,
    }

    with open("RAW_OUTPUT.json", "w") as f:
        json.dump(raw_output, f, indent=2, sort_keys=True)
    print("RAW_OUTPUT.json written.")

    print("\nDone.")

if __name__ == '__main__':
    main()
