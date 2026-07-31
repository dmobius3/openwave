#!/usr/bin/env python3
"""M8.5-A context-isolated independent-method reproduction (scalar first-occurrence table).

Implements PROTOCOL.md (frozen 2026-07-30) against the group input in m8_5a_packet.json.

Conventions (stated per PROTOCOL.md section 5.4):
  - Inner product <a,b>_G = (1/|G|) sum_g a(g) * conj(b(g)); multiplicities over C.
  - All characters computed here are exactly real elements of Q(phi); conjugation is
    therefore the identity and cannot flip any result. Realness is asserted, not assumed.
  - ALL arithmetic is exact over Q(phi) (phi^2 = phi + 1) using fractions.Fraction.
    No floating point enters the pipeline. Every tolerance is therefore fixed at
    exactly zero (see TOLERANCES below); the section 8 stability-under-increased-
    precision requirement is discharged by exactness itself.
  - The occurrence test "multiplicity > 1/2" (section 5.4) is applied literally,
    after the integer-nearness gate (G8) passes for that multiplicity.

Exit codes: 0 = success; 2 = gate/structural failure (fail loud);
            3 = mutation-harness failure (uncovered gate or unreddened mutation).
"""

import argparse
import hashlib
import json
import os
import platform
import sys
from fractions import Fraction

# --------------------------------------------------------------------------
# Frozen run parameters
# --------------------------------------------------------------------------

SCHEMA_VERSION = "m8_5a-v1"
PACKET_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "m8_5a_packet.json")
PACKET_SHA256 = "e3b0c945bbbb15b4549fa641234c9461062c2337b3d1e372af621b614d4883a9"
EXPECTED_GROUP_ORDER = 120       # permitted construction input (PROTOCOL.md section 4);
                                 # checked by G1, never used to build anything.
N_MAX = 24                       # scalar search bound (section 5.5)
M_MAX = 24                       # coexact search bound (section 6)
PEEL_N_MAX = 24                  # bound for character-extraction peeling; fail loud beyond

# Every tolerance in this run, fixed in source before any comparison (TASK req. 3):
# all comparisons are exact equalities in Q(phi); the tolerance is exactly 0.
TOLERANCES = {
    "G3_orthonormality": "0 (exact arithmetic in Q(phi))",
    "G5_Q_identification": "0 (exact)",
    "G8_integer_nearness": "0 (exact: value must BE a nonnegative rational integer)",
    "occurrence_test": "> 1/2 applied to an exact integer after G8",
}


class GateFailure(Exception):
    """A gate failed; carries the gate name. Fail loud (PROTOCOL section 5.5, TASK req. 4)."""

    def __init__(self, gate, detail):
        self.gate = gate
        self.detail = detail
        super().__init__("%s: %s" % (gate, detail))


# --------------------------------------------------------------------------
# Exact arithmetic in Q(phi), phi^2 = phi + 1
# --------------------------------------------------------------------------

class QPhi(object):
    """Element a + b*phi of Q(phi), a and b exact rationals. Immutable, hashable."""

    __slots__ = ("a", "b")

    def __init__(self, a, b=0):
        object.__setattr__(self, "a", Fraction(a))
        object.__setattr__(self, "b", Fraction(b))

    def __setattr__(self, *args):
        raise AttributeError("QPhi is immutable")

    def __add__(self, o):
        o = _coerce(o)
        return QPhi(self.a + o.a, self.b + o.b)

    def __sub__(self, o):
        o = _coerce(o)
        return QPhi(self.a - o.a, self.b - o.b)

    def __neg__(self):
        return QPhi(-self.a, -self.b)

    def __mul__(self, o):
        o = _coerce(o)
        # (a1 + b1 phi)(a2 + b2 phi) = a1a2 + b1b2  +  (a1b2 + a2b1 + b1b2) phi
        return QPhi(self.a * o.a + self.b * o.b,
                    self.a * o.b + o.a * self.b + self.b * o.b)

    __radd__ = __add__
    __rmul__ = __mul__

    def __rsub__(self, o):
        return _coerce(o) - self

    def __eq__(self, o):
        o = _coerce(o)
        return self.a == o.a and self.b == o.b

    def __ne__(self, o):
        return not self.__eq__(o)

    def __hash__(self):
        return hash((self.a, self.b))

    def galois(self):
        """Field automorphism phi -> 1 - phi:  a + b*phi -> (a+b) - b*phi."""
        return QPhi(self.a + self.b, -self.b)

    def is_zero(self):
        return self.a == 0 and self.b == 0

    def is_rational(self):
        return self.b == 0

    def is_nonneg_integer(self):
        """Exact integer-nearness (tolerance 0): a nonnegative rational integer."""
        return self.b == 0 and self.a.denominator == 1 and self.a >= 0

    def as_int(self):
        if not (self.b == 0 and self.a.denominator == 1):
            raise ValueError("not an exact rational integer: %s" % self)
        return int(self.a)

    def sort_key(self):
        return (self.a, self.b)

    def __repr__(self):
        if self.b == 0:
            return str(self.a)
        return "(%s + %s*phi)" % (self.a, self.b)


def _coerce(x):
    if isinstance(x, QPhi):
        return x
    return QPhi(x)


ZERO = QPhi(0)
ONE = QPhi(1)


# --------------------------------------------------------------------------
# Quaternions over Q(phi): tuples (w, x, y, z) of QPhi, basis (1, i, j, k)
# --------------------------------------------------------------------------

def qmul(p, q):
    w1, x1, y1, z1 = p
    w2, x2, y2, z2 = q
    return (w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
            w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
            w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
            w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2)


def qconj(q):
    w, x, y, z = q
    return (w, -x, -y, -z)


def qnorm2(q):
    w, x, y, z = q
    return w * w + x * x + y * y + z * z


QUAT_ONE = (ONE, ZERO, ZERO, ZERO)


def quat_galois(q):
    """Componentwise Galois conjugation: a group isomorphism onto the conjugate embedding."""
    return tuple(c.galois() for c in q)


# --------------------------------------------------------------------------
# Packet parsing
# --------------------------------------------------------------------------

def parse_component(s):
    """Parse '(a + b*phi)/2' with integer a, b into QPhi (a + b*phi)/2."""
    s = s.strip()
    if not (s.startswith("(") and s.endswith(")/2")):
        raise ValueError("bad component format: %r" % s)
    inner = s[1:-3]
    left, right = inner.split("+")
    a = int(left.strip())
    rb = right.strip()
    if not rb.endswith("*phi"):
        raise ValueError("bad component format: %r" % s)
    b = int(rb[:-4].strip())
    return QPhi(Fraction(a, 2), Fraction(b, 2))


def load_packet():
    with open(PACKET_FILE, "rb") as f:
        raw = f.read()
    digest = hashlib.sha256(raw).hexdigest()
    if digest != PACKET_SHA256:
        raise GateFailure("packet", "SHA-256 mismatch: %s" % digest)
    packet = json.loads(raw.decode("utf-8"))
    if packet.get("format_version") != "m8_5a-packet-v1":
        raise GateFailure("packet", "unexpected format_version")
    gens = []
    for g in packet["generators"]:
        quat = tuple(parse_component(c) for c in g)
        if qnorm2(quat) != ONE:
            raise GateFailure("packet", "generator is not a unit quaternion: %s" % (quat,))
        gens.append(quat)
    return gens


# --------------------------------------------------------------------------
# Group construction: closure, conjugacy classes, power map
# --------------------------------------------------------------------------

def generate_group(gens, mutation=None):
    """Breadth-first closure of the generator set under quaternion multiplication."""
    elements = [QUAT_ONE]
    seen = {QUAT_ONE}
    frontier = [QUAT_ONE]
    while frontier:
        nxt = []
        for q in frontier:
            for g in gens:
                p = qmul(q, g)
                if p not in seen:
                    seen.add(p)
                    elements.append(p)
                    nxt.append(p)
        frontier = nxt
        if len(elements) > 10000:
            raise GateFailure("G1", "closure exceeded 10000 elements; not a small finite group")
    if mutation == "G1_drop_element":
        elements = elements[:-1]
        seen = set(elements)
    return elements, seen


def check_closure(elements, seen, gens):
    """True iff the element set is closed under right multiplication by every generator."""
    for q in elements:
        for g in gens:
            if qmul(q, g) not in seen:
                return False
    return True


def conjugacy_classes(elements, seen, gens, mutation=None):
    """Partition into conjugacy classes: orbits under conjugation by the generators."""
    inv_gens = [qconj(g) for g in gens]  # unit quaternions: inverse = conjugate
    class_of = {}
    classes = []
    for e in elements:
        if e in class_of:
            continue
        orbit = [e]
        oseen = {e}
        stack = [e]
        while stack:
            x = stack.pop()
            for g, gi in zip(gens, inv_gens):
                y = qmul(qmul(g, x), gi)
                if y not in oseen:
                    oseen.add(y)
                    orbit.append(y)
                    stack.append(y)
        idx = len(classes)
        classes.append(orbit)
        for x in orbit:
            class_of[x] = idx
    # Deterministic class order: by (size, trace sort key); trace = 2*Re is a class invariant.
    order = sorted(range(len(classes)),
                   key=lambda i: (len(classes[i]),
                                  (classes[i][0][0] + classes[i][0][0]).sort_key()))
    classes = [classes[i] for i in order]
    class_of = {}
    for i, cl in enumerate(classes):
        for x in cl:
            class_of[x] = i
    if mutation == "G2_swap_elements":
        # Swap one element between the two largest classes: breaks conjugation-closure.
        a, b = len(classes) - 1, len(classes) - 2
        classes[a][0], classes[b][0] = classes[b][0], classes[a][0]
        class_of[classes[a][0]] = a
        class_of[classes[b][0]] = b
    return classes, class_of


# --------------------------------------------------------------------------
# Complex numbers over Q(phi) and explicit matrix representations
# --------------------------------------------------------------------------
# A complex value is a pair (re, im) of QPhi. Matrices are tuples of row-tuples.

CZERO = (ZERO, ZERO)
CONE = (ONE, ZERO)


def cadd(u, v):
    return (u[0] + v[0], u[1] + v[1])


def cmul(u, v):
    return (u[0] * v[0] - u[1] * v[1], u[0] * v[1] + u[1] * v[0])


def su2_matrix(q):
    """Packet's frozen embedding: quaternion w + xi + yj + zk as the SU(2) matrix
    [[w + x i,  y + z i], [-y + z i,  w - x i]]  (trace = 2w = 2 Re q)."""
    w, x, y, z = q
    return (((w, x), (y, z)),
            ((-y, z), (w, -x)))


def mat_trace_real(M):
    """Trace of a complex matrix, asserted exactly real; returns the QPhi real part."""
    re, im = ZERO, ZERO
    for i in range(len(M)):
        re = re + M[i][i][0]
        im = im + M[i][i][1]
    if not im.is_zero():
        raise GateFailure("realness", "matrix trace has nonzero imaginary part: %s" % im)
    return re


def sym2_matrix(M):
    """Explicit symmetric square: the induced matrix on the monomial basis
    {e_i e_j : i <= j}. This is the constructed-object route of TASK requirement 2:
    the coefficient character is the TRACE of this matrix, not the character identity."""
    d = len(M)
    basis = [(i, j) for i in range(d) for j in range(i, d)]
    S = []
    for (k, l) in basis:
        row = []
        for (i, j) in basis:
            if i == j:
                if k == l:
                    c = cmul(M[k][i], M[k][i])
                else:
                    c = cmul((QPhi(2), ZERO), cmul(M[k][i], M[l][i]))
            else:
                if k == l:
                    c = cmul(M[k][i], M[k][j])
                else:
                    c = cadd(cmul(M[k][i], M[l][j]), cmul(M[l][i], M[k][j]))
            row.append(c)
        S.append(tuple(row))
    return tuple(S)


def build_sigma_reps(class_reps):
    """Explicit matrices, per conjugacy-class representative, for the three declared
    flat-connection classes. 'standard' is the packet's frozen embedding itself;
    'galois' is the Galois-conjugated embedding (componentwise phi -> 1 - phi, a
    quaternion-algebra automorphism, hence a genuine 2-dim representation of the
    same abstract group). 'trivial' is the 1-dimensional identity representation."""
    reps = {}
    reps["trivial"] = [((CONE,),) for _ in class_reps]
    reps["standard"] = [su2_matrix(g) for g in class_reps]
    reps["galois"] = [su2_matrix(quat_galois(g)) for g in class_reps]
    return reps


def build_tau_chars(sigma_reps, mutation=None):
    """Coefficient representation tau_sigma = Sym^2(sigma), constructed explicitly.
    Returns {name: tuple of QPhi class-function values} - THE object consumed by the
    first-occurrence calculation and by gate G11 alike (Addendum 1 discipline)."""
    tau = {}
    for name, mats in sigma_reps.items():
        vals = []
        for M in mats:
            if mutation == "G11_chi_squared":
                t = mat_trace_real(M)
                vals.append(t * t)          # the Addendum-1 mutation: chi^2, not Sym^2
            else:
                vals.append(mat_trace_real(sym2_matrix(M)))
        tau[name] = tuple(vals)
    return tau


def build_sigma_chars(sigma_reps):
    """chi_sigma from the same constructed matrices (used by G5 and G11)."""
    return {name: tuple(mat_trace_real(M) for M in mats)
            for name, mats in sigma_reps.items()}


# --------------------------------------------------------------------------
# The scalar tower V_n restricted to the group: exact characters
# --------------------------------------------------------------------------

def build_tower(class_reps, n_max, mutation=None):
    """chi_{V_n} on each class, built from the weight sum: the SU(2) element with
    eigenvalues exp(+-i alpha) (2 cos alpha = 2 Re q =: t) gives
    chi_{V_n} = sum of the n+1 weights exp(i(n-2k) alpha), k = 0..n, evaluated
    exactly via the power sums s_m = exp(im alpha) + exp(-im alpha):
    s_0 = 2, s_1 = t, s_m = t s_{m-1} - s_{m-2}.
    The G12 mutation builds each level with n weights instead of n+1."""
    tower = []
    per_class_s = []
    for g in class_reps:
        t = g[0] + g[0]                      # 2 Re q, exact in Q(phi)
        s = [QPhi(2), t]
        for m in range(2, n_max + 1):
            s.append(t * s[m - 1] - s[m - 2])
        per_class_s.append(s)
    for n in range(0, n_max + 1):
        count = n if mutation == "G12_n_weights" else n + 1
        vals = []
        for s in per_class_s:
            top = count - 1
            total = ZERO
            m = top
            while m > 0:
                total = total + s[m]
                m -= 2
            if m == 0 and top >= 0:
                total = total + ONE
            vals.append(total)
        tower.append(tuple(vals))
    if mutation == "G8_perturb_tower":
        v = list(tower[2])
        v[0] = v[0] + QPhi(Fraction(1, 3))   # non-integer perturbation
        tower[2] = tuple(v)
    return tower


# --------------------------------------------------------------------------
# Class-function inner product and integer-nearness (G8)
# --------------------------------------------------------------------------

def make_ip(class_sizes, group_order):
    inv = Fraction(1, group_order)
    def ip(f, g):
        """<f,g>_G = (1/|G|) sum_c |C_c| f(c) conj(g(c)); all values here are exactly
        real (asserted at construction), so conjugation is the identity."""
        acc = ZERO
        for size, fv, gv in zip(class_sizes, f, g):
            acc = acc + QPhi(size) * fv * gv
        return QPhi(inv) * acc
    return ip


_G8_COUNT = [0]   # multiplicities checked in the current pipeline run (reporting only)


def multiplicity(ip, f, g, context, g8_log):
    """Every character inner product used as a multiplicity passes through here:
    the section 5.4 integer-nearness rule at tolerance exactly 0."""
    v = ip(f, g)
    _G8_COUNT[0] += 1
    if not v.is_nonneg_integer():
        g8_log.append((context, repr(v)))
        raise GateFailure("G8", "non-integer multiplicity in %s: %s" % (context, v))
    return v.as_int()


# --------------------------------------------------------------------------
# Irreducible characters by exact peeling of the restricted tower
# --------------------------------------------------------------------------

def extract_irreducibles(tower, ip, id_class, group_order, g8_log):
    """Peel irreducible characters out of the exact restrictions chi_{V_n}|_G.

    Invariant: after subtracting all known irreducible content from a restriction,
    the remainder is a genuine character (a nonnegative integer combination of
    not-yet-found irreducibles). A remainder of norm 1 IS a new irreducible
    character, exactly. Norm >= 2 remainders are pooled and re-reduced whenever a
    new irreducible lands. Completes when the sum of squared dimensions equals the
    group order; fails loud if the tower bound is exhausted first."""
    known = []
    pool = []

    def reduce_fn(f, ctx):
        r = list(f)
        for chi in known:
            m = multiplicity(ip, f, chi, "restriction multiplicity (%s)" % ctx, g8_log)
            if m:
                r = [rv - QPhi(m) * cv for rv, cv in zip(r, chi)]
        return tuple(r)

    def norm2(r, ctx):
        v = ip(r, r)
        if not v.is_nonneg_integer():
            g8_log.append((ctx, repr(v)))
            raise GateFailure("G8", "non-integer remainder norm in %s: %s" % (ctx, v))
        return v.as_int()

    def dims_done():
        s = 0
        for chi in known:
            d = chi[id_class]
            if not d.is_nonneg_integer():
                raise GateFailure("character-extraction", "non-integer dimension %s" % d)
            s += d.as_int() ** 2
        return s == group_order

    def drain_pool():
        progress = True
        while progress:
            progress = False
            for i in range(len(pool) - 1, -1, -1):
                r = reduce_fn(pool[i], "pool re-reduction")
                n2 = norm2(r, "pool remainder norm")
                if n2 == 0:
                    pool.pop(i)
                    progress = True
                elif n2 == 1:
                    pool.pop(i)
                    known.append(r)
                    progress = True

    for n in range(len(tower)):
        if dims_done():
            break
        r = reduce_fn(tower[n], "V_%d restriction" % n)
        n2 = norm2(r, "V_%d remainder norm" % n)
        if n2 == 1:
            known.append(r)
            drain_pool()
        elif n2 > 1:
            pool.append(r)
            drain_pool()

    if not dims_done():
        raise GateFailure("character-extraction",
                          "peeling did not complete within the tower bound "
                          "(sum of squared dims != group order); structural failure")
    if pool:
        raise GateFailure("character-extraction", "unresolved remainder pool nonempty")
    return known


# --------------------------------------------------------------------------
# McKay matrix, graph distances
# --------------------------------------------------------------------------

def build_mckay(irreps, chi_q, ip, g8_log, mutation=None):
    """A_{rho,tau} = <chi_rho * chi_Q, chi_tau>, derived in-implementation
    (PROTOCOL section 5.7); entries pass the integer-nearness rule before rounding."""
    k = len(irreps)
    A = []
    for r in range(k):
        prod = tuple(a * b for a, b in zip(irreps[r], chi_q))
        row = []
        for t in range(k):
            row.append(multiplicity(ip, prod, irreps[t],
                                    "McKay adjacency A[%d][%d]" % (r, t), g8_log))
        A.append(row)
    if mutation == "G9_asymmetric":
        A[0][1] += 1
    if mutation == "G6_zero_entry":
        for r in range(k):
            for t in range(k):
                if A[r][t] > 0:
                    A[r][t] = 0
                    return A
    if mutation == "G7_disconnect":
        last = k - 1
        for t in range(k):
            A[last][t] = 0
            A[t][last] = 0
    return A


def bfs_distances(A, start):
    """Graph distance from `start` over the support A > 0. None = unreachable."""
    k = len(A)
    dist = [None] * k
    dist[start] = 0
    frontier = [start]
    while frontier:
        nxt = []
        for u in frontier:
            for v in range(k):
                if A[u][v] > 0 and dist[v] is None:
                    dist[v] = dist[u] + 1
                    nxt.append(v)
        frontier = nxt
    return dist


def bipartition_ok(A, dist):
    """True iff every edge of the support graph joins opposite BFS parities.
    Computational witness for the parity lemma used by the coexact derivation."""
    k = len(A)
    for u in range(k):
        for v in range(k):
            if A[u][v] > 0 and dist[u] is not None and dist[v] is not None:
                if (dist[u] + dist[v]) % 2 == 0:
                    return False
    return True


# --------------------------------------------------------------------------
# First occurrences: scalar tower (section 5.4) and coexact tower (section 6)
# --------------------------------------------------------------------------

def occurrence_multiplicity(ip, tower_n, tau_chi, rho_chi, ctx, g8_log):
    prod = tuple(a * b for a, b in zip(tower_n, tau_chi))
    return multiplicity(ip, prod, rho_chi, ctx, g8_log)


def scalar_first_occurrences(tower, tau_chars, irreps, ip, g8_log):
    """n_first(rho, sigma) = min{ n >= 0 : <chi_{V_n} chi_{tau_sigma}, chi_rho> > 0 },
    n <= N_MAX, occurrence tested as multiplicity > 1/2 AFTER the G8 integer gate.
    Also returns the full multiplicity table (reused by the coexact module, which
    is defined over the same tower restrictions)."""
    mult_table = {}
    n_first = {}
    not_found = []
    for r, rho in enumerate(irreps):
        for name, tau in tau_chars.items():
            first = None
            for n in range(0, N_MAX + 1):
                m = occurrence_multiplicity(
                    ip, tower[n], tau, rho,
                    "first-occurrence <V_%d x tau_%s, rho_%d>" % (n, name, r), g8_log)
                mult_table[(r, name, n)] = m
                if first is None and Fraction(m) > Fraction(1, 2):
                    first = n
            if first is None:
                not_found.append((r, name))
            n_first[(r, name)] = first
    return n_first, mult_table, not_found


def coexact_first_occurrences(mult_table, irreps, tau_names):
    """Coexact module (section 6). Own harmonic analysis (method note section C):
    on S^3 = SU(2), Peter-Weyl gives Omega^1 = sum_n V_n (x) (V_n (x) V_2) in the
    left-invariant coframe; exact forms carry the (V_n, V_n) blocks, so the coexact
    eigenspace at level m (Hodge eigenvalue m^2/R^2, m >= 2) is
    (V_m (x) V_{m-2}) + (V_{m-2} (x) V_m). The deck group acts on the LEFT factor,
    so level m contributes V_m|_G + V_{m-2}|_G (right-factor dimensions are nonzero
    multiplicities, which cannot change a first occurrence). Under the section 5.4
    pattern: occurrence at level m iff either bracket multiplicity is positive.
    m_first(rho, sigma) = min{ m >= 2 : mult(m) > 0 or mult(m-2) > 0 }, m <= M_MAX."""
    m_first = {}
    not_found = []
    for r in range(len(irreps)):
        for name in tau_names:
            first = None
            for m in range(2, M_MAX + 1):
                if (mult_table[(r, name, m)] > 0
                        or mult_table[(r, name, m - 2)] > 0):
                    first = m
                    break
            if first is None:
                not_found.append((r, name))
            m_first[(r, name)] = first
    return m_first, not_found


def coexact_rule_prediction(d):
    """The ASSERTED entry rule (section 6): level d for d >= 2, 2 for d = 0, 3 for d = 1."""
    if d == 0:
        return 2
    if d == 1:
        return 3
    return d


# --------------------------------------------------------------------------
# Section 7 comparison harness (used at adjudication time; self-tested by G10)
# --------------------------------------------------------------------------

def compare_tables(rows_a, rows_b, sabotage=False):
    """Label-free comparison of two scalar tables by (dim, mckay_distance) signature.
    Returns (category, details). Exact integer equality, no tolerances."""
    if sabotage:
        return ("reproduced", [])          # G10 mutation target: a harness that cannot fail
    sig_a = [(r["dim"], r["mckay_distance"]) for r in rows_a]
    sig_b = [(r["dim"], r["mckay_distance"]) for r in rows_b]
    if len(set(sig_a)) != len(sig_a) or len(set(sig_b)) != len(sig_b):
        return ("structural failure", ["signatures not pairwise distinct"])
    if len(rows_a) != len(rows_b):
        return ("structural failure", ["row counts differ: %d vs %d"
                                       % (len(rows_a), len(rows_b))])
    if set(sig_a) != set(sig_b):
        return ("structural failure", ["signature sets differ"])
    by_sig = {(r["dim"], r["mckay_distance"]): r for r in rows_b}
    mismatches = []
    for r in rows_a:
        other = by_sig[(r["dim"], r["mckay_distance"])]
        for col in ("trivial", "standard", "galois"):
            va, vb = r["n_first"][col], other["n_first"][col]
            if va != vb:
                mismatches.append("(dim=%d,d=%d) %s: %s vs %s"
                                  % (r["dim"], r["mckay_distance"], col, va, vb))
    if mismatches:
        return ("partial disagreement", mismatches)
    return ("reproduced", [])


def perturb_rows(rows):
    """The doc_typo transcription mutation: one transcribed cell perturbed by +1."""
    import copy
    out = copy.deepcopy(rows)
    out[0]["n_first"]["trivial"] += 1
    return out


# --------------------------------------------------------------------------
# Pipeline: construct everything, evaluate every gate on the consumed objects
# --------------------------------------------------------------------------

CONSULTED_FILES = ["PROTOCOL.md", "TASK.md", "GROUP_INPUT.md", "m8_5a_packet.json"]
GATE_NAMES = ["G1", "G2", "G3", "G4", "G5", "G6",
              "G7", "G8", "G9", "G10", "G11", "G12"]


def run_pipeline(mutation=None):
    gates = []
    g8_log = []

    def gate(name, ok, detail):
        gates.append({"gate": name, "pass": bool(ok), "detail": detail})

    def enforce():
        for g in gates:
            if not g["pass"]:
                e = GateFailure(g["gate"], g["detail"])
                e.gates = gates
                raise e

    try:
        # -- group ---------------------------------------------------------
        gens = load_packet()
        elements, seen = generate_group(gens, mutation)
        closed = check_closure(elements, seen, gens)
        order = len(elements)
        gate("G1", order == EXPECTED_GROUP_ORDER and closed,
             "constructed order %d (required %d); closed under generators: %s"
             % (order, EXPECTED_GROUP_ORDER, closed))
        enforce()

        classes, class_of = conjugacy_classes(elements, seen, gens, mutation)
        sizes = [len(c) for c in classes]
        g2_ok = (sum(sizes) == order
                 and all(order % s == 0 for s in sizes)
                 and len(classes[class_of[QUAT_ONE]]) == 1)
        if g2_ok:
            for i, cl in enumerate(classes):
                for x in cl:
                    for g in gens:
                        if class_of.get(qmul(qmul(g, x), qconj(g))) != i:
                            g2_ok = False
        gate("G2", g2_ok,
             "%d classes, sizes %s: sum to %d, each divides %d, identity is a "
             "singleton, each class closed under conjugation" % (len(classes), sizes,
                                                                 sum(sizes), order))
        enforce()

        class_reps = [c[0] for c in classes]
        id_class = class_of[QUAT_ONE]
        ip = make_ip(sizes, order)

        # -- constructed objects consumed downstream ------------------------
        tower = build_tower(class_reps, max(N_MAX, PEEL_N_MAX), mutation)
        sigma_reps = build_sigma_reps(class_reps)
        if mutation == "G5_swap_embeddings":
            sigma_reps["standard"], sigma_reps["galois"] = (
                sigma_reps["galois"], sigma_reps["standard"])
        sigma_chars = build_sigma_chars(sigma_reps)
        tau_chars = build_tau_chars(sigma_reps, mutation)

        irreps = extract_irreducibles(tower, ip, id_class, order, g8_log)
        if mutation == "G3_corrupt_char":
            c = (id_class + 1) % len(classes)
            v = list(irreps[1])
            v[c] = v[c] + ONE
            irreps[1] = tuple(v)
        if mutation == "G4_drop_irrep":
            irreps.pop()

        dims = []
        for chi in irreps:
            d = chi[id_class]
            if not d.is_nonneg_integer() or d.as_int() == 0:
                raise GateFailure("character-extraction", "bad dimension %s" % d)
            dims.append(d.as_int())

        # -- character-table gates ------------------------------------------
        ortho_ok = True
        for i in range(len(irreps)):
            for j in range(i, len(irreps)):
                want = ONE if i == j else ZERO
                if ip(irreps[i], irreps[j]) != want:
                    ortho_ok = False
        gate("G3", ortho_ok,
             "character rows exactly orthonormal in <.,.>_G (tolerance 0, exact "
             "Q(phi) arithmetic); %d irreducibles" % len(irreps))

        gate("G4", sum(d * d for d in dims) == order,
             "sum of squared dimensions %s = %d (group order %d)"
             % (dims, sum(d * d for d in dims), order))

        expected_q = tuple(g[0] + g[0] for g in class_reps)  # 2 cos theta = 2 Re q, raw
        two_dims = [i for i, d in enumerate(dims) if d == 2]
        q_matches = [i for i in two_dims if irreps[i] == expected_q]
        g5_ok = (sigma_chars["standard"] == expected_q
                 and len(q_matches) == 1
                 and len(two_dims) == 2)
        q_idx = q_matches[0] if q_matches else None
        qp_idx = None
        if g5_ok:
            qp_idx = [i for i in two_dims if i != q_idx][0]
            g5_ok = sigma_chars["galois"] == irreps[qp_idx]
        gate("G5", g5_ok,
             "exactly one 2-dim irreducible has chi(g) = 2 cos theta(g) on every class "
             "under the packet embedding (that is Q); the constructed 'standard' rep "
             "has that character; the constructed 'galois' rep has the character of "
             "the unique other 2-dim irreducible (that is Q')")
        enforce()

        # -- McKay graph ------------------------------------------------------
        A = build_mckay(irreps, irreps[q_idx], ip, g8_log, mutation)
        triv_candidates = [i for i, chi in enumerate(irreps)
                           if all(v == ONE for v in chi)]
        if len(triv_candidates) != 1:
            raise GateFailure("character-extraction", "trivial character not unique")
        triv_idx = triv_candidates[0]
        dist = bfs_distances(A, triv_idx)

        g6_ok = all(sum(A[r][t] * dims[t] for t in range(len(dims))) == 2 * dims[r]
                    for r in range(len(dims)))
        gate("G6", g6_ok, "McKay mark condition A.dims = 2.dims holds exactly")

        sigs = [(dims[i], dist[i]) for i in range(len(dims))]
        g7_ok = all(d is not None for d in dist) and len(set(sigs)) == len(sigs)
        gate("G7", g7_ok,
             "graph distance defined for every irreducible; (dim, distance) "
             "signatures pairwise distinct: %s" % sorted(
                 s for s in sigs if s[1] is not None))

        g9_ok = all(A[r][t] == A[t][r]
                    for r in range(len(dims)) for t in range(len(dims)))
        gate("G9", g9_ok, "McKay tensor-multiplicity matrix is symmetric")
        enforce()

        bipartite = bipartition_ok(A, dist)

        # -- first occurrences (scalar) --------------------------------------
        n_first, mult_table, nf = scalar_first_occurrences(
            tower, tau_chars, irreps, ip, g8_log)
        if nf:
            gate("NOT-FOUND", False,
                 "no occurrence within n <= %d for %s" % (N_MAX, nf))
            enforce()

        gate("G8", not g8_log,
             "integer-nearness held for all %d multiplicities computed "
             "(restriction, adjacency, first-occurrence), at tolerance exactly 0"
             % _G8_COUNT[0])

        g12_ok = all(tower[n][id_class] == QPhi(n + 1) for n in range(N_MAX + 1))
        gate("G12", g12_ok,
             "the consumed V_n character object satisfies chi_{V_n}(e) = n+1 "
             "for 0 <= n <= %d" % N_MAX)

        half = QPhi(Fraction(1, 2))
        g11_ok = True
        g11_detail = []
        for name in ("trivial", "standard", "galois"):
            for c, rep in enumerate(class_reps):
                c2 = class_of[qmul(rep, rep)]     # raw group multiplication
                lhs = tau_chars[name][c]          # the consumed coefficient character
                rhs = (sigma_chars[name][c] * sigma_chars[name][c]
                       + sigma_chars[name][c2]) * half
                if lhs != rhs:
                    g11_ok = False
                    g11_detail.append("%s class %d: %s vs %s" % (name, c, lhs, rhs))
        gate("G11", g11_ok,
             "consumed tau_sigma character equals (chi_sigma(g)^2 + chi_sigma(g^2))/2 "
             "on every class for all three connection classes"
             + ("" if g11_ok else "; mismatches: %s" % g11_detail[:4]))

        # -- output rows ------------------------------------------------------
        row_order = sorted(range(len(irreps)), key=lambda i: (dist[i], dims[i]))
        rows = []
        for i in row_order:
            rows.append({
                "dim": dims[i],
                "mckay_distance": dist[i],
                "n_first": {name: n_first[(i, name)]
                            for name in ("trivial", "standard", "galois")},
            })

        sab = (mutation == "G10_sabotage_compare")
        cat_same, _ = compare_tables(rows, rows, sabotage=sab)
        cat_mut, mut_detail = compare_tables(rows, perturb_rows(rows), sabotage=sab)
        gate("G10", cat_same == "reproduced" and cat_mut == "partial disagreement",
             "comparison harness: identity compares 'reproduced'; doc_typo "
             "transcription mutation compares '%s' (%s)" % (cat_mut, mut_detail[:1]))
        enforce()

        # -- coexact module (section 6; pre-declared to RUN) -------------------
        m_first, nf2 = coexact_first_occurrences(
            mult_table, irreps, ("trivial", "standard", "galois"))
        witness_nd = all(n_first[(i, "trivial")] == dist[i]
                         for i in range(len(irreps)))
        rule_cells = []
        rule_all_match = True
        for i in row_order:
            pred = coexact_rule_prediction(dist[i])
            comp = m_first[(i, "trivial")]
            match = (pred == comp)
            rule_all_match = rule_all_match and match
            rule_cells.append({"dim": dims[i], "mckay_distance": dist[i],
                               "predicted": pred, "computed": comp, "match": match})
        if nf2:
            verdict = "not resolved"
            verdict_reason = ("NOT-FOUND cells within m <= %d: %s" % (M_MAX, nf2))
        elif not rule_all_match:
            verdict = "contradicted"
            verdict_reason = ("at least one trivial-column cell disagrees under "
                              "exact arithmetic and the stated convention map")
        elif bipartite and witness_nd:
            verdict = "structurally derived and reproduced"
            verdict_reason = (
                "general argument supplied (method_note_draft.md section C): "
                "n_first(rho, trivial) = McKay distance via the Clebsch-Gordan "
                "recursion, plus bipartite parity; coexact level-m content is "
                "V_m + V_{m-2} with m >= 2; computed table agrees cell-by-cell")
        else:
            verdict = "numerically consistent, not derived"
            verdict_reason = ("computed first occurrences agree across the declared "
                              "range, but a supporting witness of the general "
                              "argument failed (bipartite=%s, n_first=distance=%s)"
                              % (bipartite, witness_nd))

        module = {
            "ran": True,
            "declared_before_unsealing": True,
            "own_convention": ("coexact level m: Hodge-Laplacian eigenvalue m^2/R^2 on "
                               "(V_m x V_{m-2}) + (V_{m-2} x V_m), m >= 2"),
            "convention_map": ("identity on levels: this implementation's level m maps "
                               "to the source convention's coexact level m; both print "
                               "eigenvalue m^2/R^2"),
            "rule_scope": ("the asserted rule references only the McKay distance d, so "
                           "it is adjudicated against the trivial-connection column; "
                           "all three columns are reported as data"),
            "m_first_rows": [{
                "dim": dims[i],
                "mckay_distance": dist[i],
                "m_first": {name: m_first[(i, name)]
                            for name in ("trivial", "standard", "galois")},
            } for i in row_order],
            "rule_adjudication": rule_cells,
            "supporting_checks": {
                "graph_bipartite": bipartite,
                "n_first_trivial_equals_mckay_distance": witness_nd,
            },
            "verdict": verdict,
            "verdict_reason": verdict_reason,
            "standing_note": ("per PROTOCOL section 6, numerical agreement never "
                              "upgrades the rule's standing; only the derivation can, "
                              "and it is examined in the adversarial audit"),
        }

        return {
            "gates": gates,
            "rows": rows,
            "module": module,
            "order": order,
            "num_classes": len(classes),
            "class_sizes": sizes,
            "dims": dims,
            "distances": dist,
            "q_index": q_idx,
            "qprime_index": qp_idx,
            "g8_checked": _G8_COUNT[0],
        }
    except GateFailure as e:
        if not hasattr(e, "gates"):
            e.gates = gates
        raise


# --------------------------------------------------------------------------
# Mutation harness (PROTOCOL section 8): every PASS line must be shown to fail
# --------------------------------------------------------------------------

MUTATIONS = [
    ("G1_drop_element", "G1", "one element removed from the closure result"),
    ("G2_swap_elements", "G2", "two elements swapped between conjugacy classes"),
    ("G3_corrupt_char", "G3", "one character value shifted by +1"),
    ("G4_drop_irrep", "G4", "one irreducible removed from the extracted set"),
    ("G5_swap_embeddings", "G5", "standard and galois embeddings interchanged"),
    ("G6_zero_entry", "G6", "one positive McKay adjacency entry zeroed"),
    ("G7_disconnect", "G7", "one irreducible's adjacency row and column zeroed"),
    ("G8_perturb_tower", "G8", "a tower character value perturbed by +1/3"),
    ("G9_asymmetric", "G9", "A[0][1] incremented, breaking symmetry"),
    ("G10_sabotage_compare", "G10", "comparison harness forced to report agreement"),
    ("G11_chi_squared", "G11", "Sym^2(sigma) replaced by chi_sigma^2 (Addendum 1)"),
    ("G12_n_weights", "G12", "tower built with n rather than n+1 weights (Addendum 1)"),
]


def run_mutation_tests():
    """For each mutation: rerun the full pipeline and demand the TARGET gate goes red.
    Exits nonzero if any mutation fails to redden its target or any gate is uncovered."""
    print("== mutation harness ==")
    failures = 0
    covered = set()
    for name, target, desc in MUTATIONS:
        _G8_COUNT[0] = 0
        try:
            run_pipeline(mutation=name)
            reddened = set()
        except GateFailure as e:
            reddened = {g["gate"] for g in e.gates if not g["pass"]}
            reddened.add(e.gate)
        ok = target in reddened
        if ok:
            covered.add(target)
        else:
            failures += 1
        print("%-24s target %-4s -> reddened %-28s %s   (%s)"
              % (name, target, sorted(reddened) or "NOTHING",
                 "OK" if ok else "MUTATION NOT DETECTED", desc))
    uncovered = [g for g in GATE_NAMES if g not in covered]
    _G8_COUNT[0] = 0
    try:
        run_pipeline(mutation=None)
        clean_ok = True
    except GateFailure as e:
        clean_ok = False
        print("clean pipeline FAILED under no mutation: %s" % e)
    print("gates covered and reddened: %s" % sorted(covered))
    if uncovered:
        print("UNCOVERED GATES: %s" % uncovered)
    if failures or uncovered or not clean_ok:
        print("MUTATION HARNESS: FAIL (%d unreddened, %d uncovered, clean=%s)"
              % (failures, len(uncovered), clean_ok))
        return 3
    print("MUTATION HARNESS: PASS (%d mutations, all %d gates covered, "
          "clean run green)" % (len(MUTATIONS), len(GATE_NAMES)))
    return 0


# --------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------

def environment_record():
    return {
        "interpreter": "python %s" % platform.python_version(),
        "implementation": platform.python_implementation(),
        "libraries": "standard library only (fractions, json, hashlib); no numpy",
        "os": "%s %s" % (platform.system(), platform.release()),
        "hardware": platform.machine(),
        "seeds": "none: the pipeline is deterministic and fully exact over Q(phi)",
    }


def print_report(res):
    print("M8.5-A context-isolated independent-method reproduction")
    print("schema_version: %s" % SCHEMA_VERSION)
    print("packet: %s sha256 %s" % (os.path.basename(PACKET_FILE), PACKET_SHA256))
    env = environment_record()
    for k in sorted(env):
        print("env.%s: %s" % (k, env[k]))
    print("")
    print("group: order %d, %d conjugacy classes, sizes %s"
          % (res["order"], res["num_classes"], res["class_sizes"]))
    print("irreducible dimensions: %s" % sorted(res["dims"]))
    print("")
    print("gates (every tolerance fixed at exactly 0; %d multiplicities checked):"
          % res["g8_checked"])
    for g in res["gates"]:
        print("  %-4s %s  %s" % (g["gate"], "PASS" if g["pass"] else "FAIL",
                                 g["detail"]))
    print("")
    print("scalar first-occurrence table (rows label-free, by (dim, McKay distance)):")
    print("  dim  distance  trivial  standard  galois")
    for r in res["rows"]:
        print("  %3d  %8d  %7d  %8d  %6d"
              % (r["dim"], r["mckay_distance"], r["n_first"]["trivial"],
                 r["n_first"]["standard"], r["n_first"]["galois"]))
    print("")
    mod = res["module"]
    print("coexact module (section 6, ASSERTED rule; ran: %s):" % mod["ran"])
    print("  convention map: %s" % mod["convention_map"])
    print("  m_first table:")
    print("  dim  distance  trivial  standard  galois")
    for r in mod["m_first_rows"]:
        print("  %3d  %8d  %7d  %8d  %6d"
              % (r["dim"], r["mckay_distance"], r["m_first"]["trivial"],
                 r["m_first"]["standard"], r["m_first"]["galois"]))
    print("  rule adjudication (trivial column): %s"
          % ("all cells match" if all(c["match"] for c in mod["rule_adjudication"])
             else [c for c in mod["rule_adjudication"] if not c["match"]]))
    print("  supporting checks: %s" % mod["supporting_checks"])
    print("  verdict: %s" % mod["verdict"])
    print("  reason: %s" % mod["verdict_reason"])
    print("  standing: %s" % mod["standing_note"])
    print("")
    print("claim ceiling (section 2): context-isolated independent-method "
          "reproduction at most; adjudication against the pinned section 6.1 "
          "table is the maintainers' step and is NOT performed here.")
    print("run status: completed; all gates PASS; awaiting section 7 adjudication")


def write_result_json(res, path):
    out = {
        "schema_version": SCHEMA_VERSION,
        "rows": res["rows"],
        "environment": environment_record(),
        "gates": res["gates"],
        "tolerances": TOLERANCES,
        "consulted_files": CONSULTED_FILES + [
            "see manifest.md for automatically loaded context"],
        "coexact_module": res["module"],
        "group": {"order": res["order"], "num_classes": res["num_classes"],
                  "class_sizes": res["class_sizes"],
                  "irreducible_dimensions": sorted(res["dims"])},
        "claim_ceiling": ("context-isolated independent-method reproduction "
                          "(PROTOCOL section 2); no stronger label available"),
        "packet_sha256": PACKET_SHA256,
    }
    with open(path, "w") as f:
        json.dump(out, f, indent=2, sort_keys=True)
        f.write("\n")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--mutation-tests", action="store_true",
                    help="run the section 8 mutation harness and exit")
    ap.add_argument("--json-out", default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "result.json"))
    args = ap.parse_args(argv)

    if args.mutation_tests:
        return run_mutation_tests()

    _G8_COUNT[0] = 0
    try:
        res = run_pipeline(mutation=None)
    except GateFailure as e:
        print("STRUCTURAL FAILURE (fail loud): gate %s: %s" % (e.gate, e.detail))
        for g in getattr(e, "gates", []):
            print("  %-4s %s  %s" % (g["gate"], "PASS" if g["pass"] else "FAIL",
                                     g["detail"]))
        return 2
    print_report(res)
    write_result_json(res, args.json_out)
    print("result.json written: %s" % args.json_out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
