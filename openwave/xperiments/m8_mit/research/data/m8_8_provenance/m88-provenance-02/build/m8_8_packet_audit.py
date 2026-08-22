"""M8.8 CONSTRUCTION-PACKET AUDIT: the maintainer check that runs before the clean room opens.

WHAT THIS IS.  The [M8.8 reproduction protocol](../findings/m8_8_reproduction_protocol.md)
section 4.2 requires the construction packet entering the clean room to be audited by
someone other than the context that will implement from it, and names the audit load-bearing
for section 1 branch four: the internal model gates cannot identify the intended manifold
complex on their own.  This file is the maintainer side of that requirement, and it is built
to the same two properties as the [M8.5-A packet audit](m8_5a_packet_audit.py):

  1  IT RUNS FROM THE TWO PUBLIC PACKETS ALONE.  It reads the construction packet and the
     M8.5-A group packet.  It reads NEITHER the author-side verification script nor the
     author-side audit artifact, which stay outside the room and are set beside this result
     afterwards.  An audit that reads the supplier's report and concurs is a review of the
     report.
  2  IT IS MECHANICAL.  Every check is a count or an exact equality that can go red, and
     `--mutation-tests` shows each one doing so.

WHAT THE PACKET IS.  A finite based 3-dimensional chain complex over Z[2I] in degrees 0..3,
free ranks [1, 2, 2, 1], boundary maps as matrices whose entries are group-ring elements
encoded as (coefficient, element_id) terms against the 120-element closure of the group
packet.  Integers and element IDs only: no evaluated matrix, no decimal rendering.

WHAT THIS VERIFIES.  Eleven checks, each mutation-tested:

  A1  format: exactly the declared key set, degree_range [0, 3], truncation_rule consistent
      with model_kind, matrix dimensions matching free_ranks, every numeric leaf an integer
  A2  canonicalization: the delivered bytes are the declared canonical form, ASCII with LF
      and a single trailing newline, and the authoritative hash is issued over that form
  A3  group binding: the recorded group-packet hash reproduces from the group packet's own
      bytes; the abstract generators resolve inside the independently rebuilt closure; the
      presentation relations hold; and the two generators generate all 120 elements
  A4  the complex is a complex: d_(n+1) d_n = 0 over Z[2I]
  A5  census: free ranks as declared and chi = 0, a closed 3-manifold's complex rather than
      a presentation 2-complex
  A6  trivial sector: eps(d1) = 0, eps(d2) unimodular, eps(d3) = 0, hence
      H_*(Z (x) C_*) = (Z, 0, 0, Z) through the DECLARED augmentation.  d1 is separately
      confirmed to be the (g - 1) 1-cell correspondence and NOT the augmentation
  A7  universal cover: H_*(C_* (x) F_p) = (1, 0, 0, 1) for p in {2, 3, 5, 7, 10007}
  A8  per-irrep acyclicity, against nine independently constructed irreducibles: the trivial
      representation NON-acyclic (expected, and a pass), every nontrivial one acyclic
  A9  leakage: no answer-bearing vocabulary outside the fields the schema requires prose in,
      no decimal literal, no key outside the declared set
  A10 relator provenance: d1 is the (g - 1) correspondence, and each d2 row is the Fox
      jacobian of exactly one word in the two-syllable family over {s, t, st, ts}, whose
      rendering must appear in the packet's own declared basis_order
  A11 exact top boundary: im(d3) = ker(d2) as integral lattices, certified with no prime
      set on the accept side

A7 AND A11 ARE THIS AUDIT'S ADDITIONS, AND THE PAIR IS THE POINT.  Every other model gate in
the protocol is blind to an overall change of d3 inside ker(d2): eps(d3) = 0 either way, the
trivial-sector homology is computed from eps(d2) alone, and per-irrep acyclicity is a rank
statement.  So `d3 -> 2 d3` passes A1 through A6 and A8 while multiplying every torsion
value by 2^(+-dim).  Viewing C_* as a Z-complex restores the missing information: it is then
the cellular chain complex of the universal cover S^3, whose homology is (Z, 0, 0, Z) at
every prime.  A rank drop mod p means im(d3) sits at finite index inside ker(d2), which is
exactly the case a rational-rank acceptance predicate cannot see.

A7 IS A REJECT SCREEN AND NOT AN ACCEPT CRITERION, WHICH IS WHY A11 EXISTS.  A drop at any
tested prime proves a finite index, but passing at every tested prime proves nothing: no
finite prime set can exclude an index supported elsewhere.  `d3 -> 11 d3` leaves all five of
A7's Betti profiles correct while multiplying every torsion value by 11^(+-dim), and the
mutation suite runs exactly that case.  A11 is the accept side and is prime-set independent:
containment comes from d3 d2 = 0, the ranks are pinned by mod-p lower bounds meeting
ceilings that the chain relations supply, and im(d3) is shown saturated by unit pivots
alone.  Two saturated lattices of equal rank, one inside the other, are equal.

WHAT THIS DOES NOT VERIFY.  That the packet came from where its provenance record says it
did: the provenance material is held maintainer-side until commitment (section 4.2 check 6),
and this audit checks the object, not its origin.  It does not verify any torsion value, and
it evaluates none.  Nor does it establish that the complex is THE intended model rather than
A model with the right invariants; A7 and A11 narrow that gap considerably but do not close
it, which is why section 1 branch four stays open and bounded rather than excluded.

The residue is worth naming precisely, because A11 makes it sharper rather than smaller.
C_3 has rank 1, so any two generators of ker(d2) differ by d3' = u d3 for a unit u of
Z[2I]/(N), and for every nontrivial irrep rho, which kills N, the torsion moves by
det rho(u).  The protocol's basis-invariance argument covers u = +-g, where det is a root of
unity and the modulus is exactly 1.  It does not cover the rest: 2I has 9 complex irreps in
7 Galois orbits with every character real, so the Whitehead group has rank 9 - 7 = 2 and
nontrivial units exist.  A11 certifies that the packet's d3 GENERATES ker(d2); it does not
certify WHICH generator, and the reproduced modulus depends on that choice.  Closing this
needs provenance, not another gate, which is exactly the division of labour section 4.2 sets
up between the model gates and the maintainer-side construction audit.

USAGE.
    python3 m8_8_packet_audit.py                     audit, write JSON, exit 0 on all-pass
    python3 m8_8_packet_audit.py --mutation-tests    additionally run the mutation suite
    python3 m8_8_packet_audit.py --packet PATH       audit a packet elsewhere
    python3 m8_8_packet_audit.py --no-write          print only, write no artifact

Exit code is nonzero if any check fails, so a red audit cannot be mistaken for a green one.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import re
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any, Callable, Iterator, Sequence

import numpy as np

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
DEFAULT_PACKET = DATA / "m8_8_construction_packet.json"
DEFAULT_GROUP_PACKET = DATA / "m8_5a_packet.json"
DEFAULT_OUT = DATA / "m8_8_packet_audit.json"

DECLARED_KEYS = {
    "format_version",
    "group_packet_sha256",
    "abstract_generators",
    "model_kind",
    "degree_range",
    "free_ranks",
    "boundary_maps",
    "top_closure",
    "truncation_rule",
    "basing",
    "provenance_id",
}
MODEL_KINDS = {"finite_cellular", "resolution_derived"}
COVER_PRIMES = (2, 3, 5, 7, 10007)
MUTATION_PRIMES = (2, 3)
# largest p for which the int64 elimination in rank_mod_p is exact: (p-1)^2 < 2^63
MAX_INT64_PRIME = 3_037_000_499
# the rank lower bound for A11; any prime in range works, a large one is simply less
# likely to need a second opinion
BOUND_PRIME = 2_147_483_647
# entries in the A11 elimination stay tiny on a real packet; a runaway means the greedy
# pivot order met something pathological, and the honest answer there is "inconclusive"
A11_GROWTH_CAP = 1 << 40
LEAKAGE_VOCABULARY = (
    "torsion",
    "character",
    "determinant",
    "eigen",
    "zeta",
    "irrep",
    "ratio",
    "sqrt",
    "phi",
)

# --------------------------------------------------------------------------------------
# Q(phi), quaternions over it, and the group
# --------------------------------------------------------------------------------------

QElem = tuple[Fraction, Fraction]  # (a, b) means a + b*phi, phi^2 = phi + 1
Quaternion = tuple[QElem, QElem, QElem, QElem]
QZERO: QElem = (Fraction(0), Fraction(0))
QONE: QElem = (Fraction(1), Fraction(0))
COMPONENT = re.compile(r"^\((-?\d+) \+ (-?\d+)\*phi\)/2$")


def q_add(x: QElem, y: QElem) -> QElem:
    return (x[0] + y[0], x[1] + y[1])


def q_sub(x: QElem, y: QElem) -> QElem:
    return (x[0] - y[0], x[1] - y[1])


def q_mul(x: QElem, y: QElem) -> QElem:
    a, b = x
    c, d = y
    return (a * c + b * d, a * d + b * c + b * d)


def q_galois(x: QElem) -> QElem:
    """The nontrivial automorphism of Q(phi): phi -> 1 - phi."""
    return (x[0] + x[1], -x[1])


def h_mul(p: Quaternion, q: Quaternion) -> Quaternion:
    a1, b1, c1, d1 = p
    a2, b2, c2, d2 = q
    return (
        q_sub(q_sub(q_sub(q_mul(a1, a2), q_mul(b1, b2)), q_mul(c1, c2)), q_mul(d1, d2)),
        q_add(q_sub(q_add(q_mul(a1, b2), q_mul(b1, a2)), q_mul(d1, c2)), q_mul(c1, d2)),
        q_add(q_add(q_sub(q_mul(a1, c2), q_mul(b1, d2)), q_mul(c1, a2)), q_mul(d1, b2)),
        q_add(q_add(q_sub(q_mul(a1, d2), q_mul(c1, b2)), q_mul(b1, c2)), q_mul(d1, a2)),
    )


H_ONE: Quaternion = (QONE, QZERO, QZERO, QZERO)


def parse_component(text: str) -> QElem:
    match = COMPONENT.match(text)
    if not match:
        raise ValueError(f"component does not parse as (a + b*phi)/2: {text!r}")
    return (Fraction(int(match.group(1)), 2), Fraction(int(match.group(2)), 2))


def canonical_key(h: Quaternion) -> tuple[int, ...]:
    """The packet's canonical exact coordinate tuple, in (a + b*phi)/2 form."""
    out: list[int] = []
    for a, b in h:
        doubled_a, doubled_b = 2 * a, 2 * b
        if doubled_a.denominator != 1 or doubled_b.denominator != 1:
            raise ValueError("component is not of the form (a + b*phi)/2 with integer a, b")
        out.extend((int(doubled_a), int(doubled_b)))
    return tuple(out)


class Group:
    """The 120-element closure, its canonical enumeration, and its Cayley table.

    Element IDs are the rank of each element's canonical coordinate tuple in lexicographic
    order, exactly as the protocol section 4.2 encoding rule defines them.  Both sides derive
    the same enumeration from the group packet alone, which is what makes the packet's
    element IDs meaningful without any further agreement.
    """

    def __init__(self, generators: Sequence[Quaternion]) -> None:
        seen: dict[tuple[int, ...], Quaternion] = {canonical_key(H_ONE): H_ONE}
        frontier = [H_ONE]
        while frontier:
            nxt: list[Quaternion] = []
            for x in frontier:
                for g in generators:
                    y = h_mul(x, g)
                    key = canonical_key(y)
                    if key not in seen:
                        seen[key] = y
                        nxt.append(y)
            frontier = nxt

        self.size = len(seen)
        keys = sorted(seen)
        self.element = [seen[k] for k in keys]
        index = {k: i for i, k in enumerate(keys)}
        self.mul = [
            [
                index[canonical_key(h_mul(self.element[i], self.element[j]))]
                for j in range(self.size)
            ]
            for i in range(self.size)
        ]
        self.identity = index[canonical_key(H_ONE)]
        self.inv = [
            next(j for j in range(self.size) if self.mul[i][j] == self.identity)
            for i in range(self.size)
        ]

    def power(self, i: int, n: int) -> int:
        if n < 0:
            return self.power(self.inv[i], -n)
        result = self.identity
        for _ in range(n):
            result = self.mul[result][i]
        return result

    def order(self, i: int) -> int:
        n, y = 1, i
        while y != self.identity:
            y = self.mul[y][i]
            n += 1
        return n

    def generated_by(self, seeds: Sequence[int]) -> int:
        have = {self.identity}
        frontier = [self.identity]
        while frontier:
            nxt: list[int] = []
            for x in frontier:
                for g in seeds:
                    y = self.mul[x][g]
                    if y not in have:
                        have.add(y)
                        nxt.append(y)
            frontier = nxt
        return len(have)


# --------------------------------------------------------------------------------------
# group-ring arithmetic
# --------------------------------------------------------------------------------------

RingElement = dict[int, int]
RingMatrix = list[list[RingElement]]


def ring(terms: Sequence[Sequence[int]]) -> RingElement:
    acc: RingElement = {}
    for coefficient, element in terms:
        acc[element] = acc.get(element, 0) + coefficient
    return {g: c for g, c in acc.items() if c}


def ring_add(x: RingElement, y: RingElement) -> RingElement:
    acc = dict(x)
    for g, c in y.items():
        acc[g] = acc.get(g, 0) + c
    return {g: c for g, c in acc.items() if c}


def ring_mul(x: RingElement, y: RingElement, group: Group) -> RingElement:
    acc: RingElement = {}
    for g, a in x.items():
        row = group.mul[g]
        for h, b in y.items():
            k = row[h]
            acc[k] = acc.get(k, 0) + a * b
    return {g: c for g, c in acc.items() if c}


def ring_matrix(raw: Sequence[Sequence[Sequence[Sequence[int]]]]) -> RingMatrix:
    return [[ring(entry) for entry in row] for row in raw]


def ring_matmul(a: RingMatrix, b: RingMatrix, group: Group) -> RingMatrix:
    inner = len(b)
    out: RingMatrix = []
    for row in a:
        new_row: list[RingElement] = []
        for j in range(len(b[0])):
            acc: RingElement = {}
            for k in range(inner):
                acc = ring_add(acc, ring_mul(row[k], b[k][j], group))
            new_row.append(acc)
        out.append(new_row)
    return out


def is_zero_matrix(m: RingMatrix) -> bool:
    return all(not entry for row in m for entry in row)


def augment(m: RingMatrix) -> list[list[int]]:
    """Apply the declared augmentation eps: every group element maps to 1."""
    return [[sum(entry.values()) for entry in row] for row in m]


# --------------------------------------------------------------------------------------
# the universal-cover complex, over F_p
# --------------------------------------------------------------------------------------


def regular_evaluate(m: RingMatrix, group: Group) -> np.ndarray:
    """Evaluate a group-ring matrix in the regular representation.

    The result is the corresponding boundary map of the cellular chain complex of the
    universal cover, as an integer matrix.  Blocks are permutation sums, so this is built
    entry by entry rather than by multiplying 120 x 120 matrices.
    """
    n = group.size
    rows, cols = len(m), len(m[0])
    out = np.zeros((rows * n, cols * n), dtype=np.int64)
    for i in range(rows):
        for j in range(cols):
            for g, coefficient in m[i][j].items():
                row_of = group.mul[group.inv[g]]
                for a in range(n):
                    out[i * n + a, j * n + row_of[a]] += coefficient
    return out


def rank_mod_p(matrix: np.ndarray, p: int) -> int:
    """Exact rank over F_p by Gaussian elimination on an int64 array.

    The int64 path is exact only while the elimination's products fit: the update forms
    `factor * pivot_row` with both below p, so p must satisfy (p-1)^2 < 2^63.  Above that
    the products wrap SILENTLY and the function returns a wrong rank rather than raising,
    which is the worst failure mode a gate can have.  Overflow is refused rather than
    absorbed: every caller here is well inside the bound, and a caller that is not has made
    a mistake worth hearing about.
    """
    if p > MAX_INT64_PRIME:
        raise ValueError(
            f"p={p} exceeds the int64-safe bound {MAX_INT64_PRIME}; "
            "the elimination would overflow and return a wrong rank"
        )
    a = np.mod(matrix, p).astype(np.int64)
    rows, cols = a.shape
    rank = 0
    for col in range(cols):
        pivot = None
        column = a[rank:, col]
        nonzero = np.nonzero(column)[0]
        if nonzero.size:
            pivot = rank + int(nonzero[0])
        if pivot is None:
            continue
        if pivot != rank:
            a[[rank, pivot]] = a[[pivot, rank]]
        a[rank] = (a[rank] * pow(int(a[rank, col]), p - 2, p)) % p
        factors = a[:, col].copy()
        factors[rank] = 0
        nz = np.nonzero(factors)[0]
        if nz.size:
            a[nz] = (a[nz] - np.outer(factors[nz], a[rank])) % p
        rank += 1
        if rank == rows:
            break
    return rank


def cover_betti(
    matrices: Sequence[np.ndarray], ranks: Sequence[int], size: int, p: int
) -> tuple[int, ...]:
    r1, r2, r3 = (rank_mod_p(m, p) for m in matrices)
    c0, c1, c2, c3 = (r * size for r in ranks)
    return (c0 - r1, c1 - r1 - r2, c2 - r2 - r3, c3 - r3)


def saturation_certificate(matrix: np.ndarray) -> tuple[bool | None, dict[str, Any]]:
    """Is the row image of an integer matrix SATURATED, that is, a direct summand?

    Equivalently: are all its elementary divisors 1?  Row and column operations used here
    are unimodular, so the product of the pivots this elimination produces equals the
    product of the elementary divisors up to sign.  Every pivot being +-1 therefore settles
    the question, with no Smith form and no determinant of a 119 x 119 minor.

    Three outcomes, and the third is deliberate.  True: the elimination cleared the matrix
    using only +-1 pivots.  False: it stalled on a nonzero remainder whose entries share a
    common factor, which IS an elementary divisor above 1.  None: it stalled with coprime
    entries and no unit pivot, where this routine has not established either answer, and
    section 4.2's own rule for a one-sided predicate applies: report inconclusive rather
    than issue a verdict.  A None is treated as a red check by the caller, so the audit
    fails closed.
    """
    a = np.array(matrix, dtype=np.int64)
    rows, cols = a.shape
    live_rows = list(range(rows))
    live_cols = list(range(cols))
    pivots = 0
    while live_rows and live_cols:
        block = a[np.ix_(live_rows, live_cols)]
        if not block.any():
            return True, {"pivots": pivots, "outcome": "cleared with unit pivots"}
        units = np.argwhere(np.abs(block) == 1)
        if units.size == 0:
            common = int(np.gcd.reduce(np.abs(block[block != 0])))
            if common > 1:
                return False, {
                    "pivots": pivots,
                    "outcome": "elementary divisor above 1",
                    "common_factor": common,
                }
            return None, {"pivots": pivots, "outcome": "inconclusive: no unit pivot"}
        # prefer a unit sitting in a sparse row: the pivot row is added to every row it
        # has to clear, so a dense one is what makes entries grow
        weights = np.count_nonzero(block, axis=1)
        bi, bj = min(((int(i), int(j)) for i, j in units), key=lambda ij: weights[ij[0]])
        pivot_row, pivot_col = live_rows[bi], live_cols[bj]
        pivot = int(a[pivot_row, pivot_col])
        # only live rows matter: a retired row is a previous pivot row, and the column
        # operations that would clear it touch nothing still in play
        others = np.array([r for r in live_rows if r != pivot_row], dtype=np.intp)
        factors = a[others, pivot_col]
        touched = others[factors != 0]
        if touched.size:
            a[touched] -= np.outer(a[touched, pivot_col] * pivot, a[pivot_row])
        if int(np.abs(a).max()) > A11_GROWTH_CAP:
            return None, {"pivots": pivots, "outcome": "inconclusive: entry growth"}
        live_rows.remove(pivot_row)
        live_cols.remove(pivot_col)
        pivots += 1
    return True, {"pivots": pivots, "outcome": "cleared with unit pivots"}


def exact_top_boundary_check(
    matrices: Sequence[np.ndarray], ranks: Sequence[int], size: int, augmentation_is_zero: bool
) -> tuple[bool, dict[str, Any]]:
    """Certify `im d3 = ker d2` as integral lattices, with no prime set on the accept side.

    A mod-p battery is a REJECT screen.  A rank drop at p proves a finite index, but no
    finite prime set can exclude one: `d3 -> 11 d3` leaves every Betti number in A7 correct
    at 2, 3, 5, 7 and 10007 while multiplying every torsion value by 11^(+-dim).  This check
    is the accept side, and it is prime-set independent:

      containment   from d3 d2 = 0, already exact over Z[2I] in A4
      ranks         each mod-p lower bound is met by a ceiling from the chain relations:
                    eps d1 = 0 caps rank d1, then d2 d1 = 0 caps rank d2, then d3 d2 = 0
                    caps rank d3.  Meeting a ceiling turns a bound into an equality
      saturation    of im d3, by unit pivots alone (see saturation_certificate)
      ker d2        saturated automatically, being the kernel of an integer matrix

    Two saturated lattices of equal rank, one inside the other, are equal.
    """
    d1, d2, d3 = matrices
    c0, c1, c2 = (r * size for r in ranks[:3])
    r1, r2, r3 = (rank_mod_p(m, BOUND_PRIME) for m in matrices)
    ceilings = (c0 - 1 if augmentation_is_zero else c0, c1 - r1, c2 - r2)
    closed = (r1, r2, r3) == ceilings
    saturated, certificate = saturation_certificate(d3)
    kernel_rank = c2 - r2
    observed = {
        "ranks": [r1, r2, r3],
        "ceilings": list(ceilings),
        "bounds_close": closed,
        "rank_ker_d2": kernel_rank,
        "im_d3_saturated": saturated,
        "certificate": certificate,
    }
    return bool(closed and saturated is True and r3 == kernel_rank), observed


# --------------------------------------------------------------------------------------
# the nine irreducible representations, constructed independently
# --------------------------------------------------------------------------------------

CElem = tuple[Fraction, Fraction, Fraction, Fraction]  # a + b*phi + (c + d*phi)*i
CMatrix = list[list[CElem]]
C_ZERO: CElem = (Fraction(0),) * 4
C_ONE: CElem = (Fraction(1), Fraction(0), Fraction(0), Fraction(0))


def c_add(x: CElem, y: CElem) -> CElem:
    return (x[0] + y[0], x[1] + y[1], x[2] + y[2], x[3] + y[3])


def c_neg(x: CElem) -> CElem:
    return (-x[0], -x[1], -x[2], -x[3])


def c_mul(x: CElem, y: CElem) -> CElem:
    real = q_sub(q_mul((x[0], x[1]), (y[0], y[1])), q_mul((x[2], x[3]), (y[2], y[3])))
    imag = q_add(q_mul((x[0], x[1]), (y[2], y[3])), q_mul((x[2], x[3]), (y[0], y[1])))
    return (real[0], real[1], imag[0], imag[1])


def c_conj(x: CElem) -> CElem:
    return (x[0], x[1], -x[2], -x[3])


def c_galois(x: CElem) -> CElem:
    real, imag = q_galois((x[0], x[1])), q_galois((x[2], x[3]))
    return (real[0], real[1], imag[0], imag[1])


def c_inv(x: CElem) -> CElem:
    conjugate = c_conj(x)
    norm = c_mul(x, conjugate)
    a, b = norm[0], norm[1]
    field_conj = (a + b, -b)
    denominator = q_mul((a, b), field_conj)[0]
    return c_mul(
        conjugate,
        (field_conj[0] / denominator, field_conj[1] / denominator, Fraction(0), Fraction(0)),
    )


def c_sum(values: Sequence[CElem]) -> CElem:
    acc = C_ZERO
    for v in values:
        acc = c_add(acc, v)
    return acc


def c_matmul(a: CMatrix, b: CMatrix) -> CMatrix:
    inner = len(b)
    return [
        [c_sum([c_mul(a[i][k], b[k][j]) for k in range(inner)]) for j in range(len(b[0]))]
        for i in range(len(a))
    ]


def su2(h: Quaternion) -> CMatrix:
    """w + xi + yj + zk  ->  [[w + xI, y + zI], [-y + zI, w - xI]], a homomorphism."""
    w, x, y, z = h
    embed = lambda v: (v[0], v[1], Fraction(0), Fraction(0))  # noqa: E731
    times_i = lambda v: (Fraction(0), Fraction(0), v[0], v[1])  # noqa: E731
    return [
        [c_add(embed(w), times_i(x)), c_add(embed(y), times_i(z))],
        [c_add(c_neg(embed(y)), times_i(z)), c_add(embed(w), c_neg(times_i(x)))],
    ]


def symmetric_power(m: CMatrix, k: int) -> CMatrix:
    """Sym^k of a 2 x 2 matrix, monomial basis e1^(k-j) e2^j, columns are images.

    The column convention is load-bearing: storing the expansion by rows instead produces
    the transpose, which is an ANTI-homomorphism that still passes character-norm and rank
    checks.  `--mutation-tests` includes that mutation for exactly that reason.
    """
    n = k + 1
    out: CMatrix = [[C_ZERO] * n for _ in range(n)]
    a, b = m[0][0], m[0][1]
    c, d = m[1][0], m[1][1]
    for j in range(n):
        poly: dict[int, CElem] = {0: C_ONE}
        for factor_left, factor_right, repeat in ((a, c, k - j), (b, d, j)):
            for _ in range(repeat):
                new: dict[int, CElem] = {}
                for exponent, coefficient in poly.items():
                    new[exponent] = c_add(
                        new.get(exponent, C_ZERO), c_mul(coefficient, factor_left)
                    )
                    new[exponent + 1] = c_add(
                        new.get(exponent + 1, C_ZERO), c_mul(coefficient, factor_right)
                    )
                poly = new
        for exponent, coefficient in poly.items():
            out[exponent][j] = coefficient
    return out


def kronecker(a: CMatrix, b: CMatrix) -> CMatrix:
    p, q = len(b), len(b[0])
    return [
        [c_mul(a[i // p][j // q], b[i % p][j % q]) for j in range(len(a[0]) * q)]
        for i in range(len(a) * p)
    ]


def build_irreps(group: Group) -> dict[str, list[CMatrix]]:
    """Nine irreducibles of 2I, from the quaternion embedding and its Galois twist.

    Sym^k of the fundamental gives dimensions 1..6; the Galois twist gives the partners of
    the 2 and the 3; and the remaining 4 is 2 (x) 2', which factors through A5.  Sym^3 is
    Galois-self-conjugate, so twisting it produces no ninth representation.
    """
    fundamental = [su2(group.element[g]) for g in range(group.size)]
    twisted = [[[c_galois(x) for x in row] for row in m] for m in fundamental]
    reps: dict[str, list[CMatrix]] = {
        "R0": [[[C_ONE]] for _ in range(group.size)],
        "4a": [kronecker(fundamental[g], twisted[g]) for g in range(group.size)],
    }
    for k, name in ((1, "2a"), (2, "3a"), (3, "4b"), (4, "5"), (5, "6")):
        reps[name] = [symmetric_power(fundamental[g], k) for g in range(group.size)]
    for k, name in ((1, "2b"), (2, "3b")):
        reps[name] = [symmetric_power(twisted[g], k) for g in range(group.size)]
    return reps


def evaluate(m: RingMatrix, rep: Sequence[CMatrix]) -> CMatrix:
    dim = len(rep[0])
    rows, cols = len(m), len(m[0])
    out: CMatrix = [[C_ZERO] * (cols * dim) for _ in range(rows * dim)]
    for i in range(rows):
        for j in range(cols):
            for g, coefficient in m[i][j].items():
                scalar = (Fraction(coefficient), Fraction(0), Fraction(0), Fraction(0))
                block = rep[g]
                for a in range(dim):
                    for b in range(dim):
                        out[i * dim + a][j * dim + b] = c_add(
                            out[i * dim + a][j * dim + b], c_mul(block[a][b], scalar)
                        )
    return out


def exact_rank(m: CMatrix) -> int:
    a = [row[:] for row in m]
    rows, cols = len(a), len(a[0])
    rank = 0
    for col in range(cols):
        pivot = next((i for i in range(rank, rows) if a[i][col] != C_ZERO), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        scale = c_inv(a[rank][col])
        a[rank] = [c_mul(x, scale) for x in a[rank]]
        for i in range(rows):
            if i != rank and a[i][col] != C_ZERO:
                factor = a[i][col]
                a[i] = [c_add(a[i][k], c_neg(c_mul(factor, a[rank][k]))) for k in range(cols)]
        rank += 1
        if rank == rows:
            break
    return rank


# --------------------------------------------------------------------------------------
# Fox calculus, for the relator-provenance check
# --------------------------------------------------------------------------------------

Word = list[tuple[int, int]]


def fox_derivative(word: Word, letter: int, group: Group) -> RingElement:
    """Left Fox derivative, evaluated through the quotient map onto the group.

    d(uv) = d(u) + u.d(v),  d(x) = 1,  d(y) = 0 for y != x,  d(x^-1) = -x^-1.
    """
    acc: RingElement = {}
    prefix: RingElement = {group.identity: 1}
    for generator, exponent in word:
        if exponent == 1:
            if generator == letter:
                acc = ring_add(acc, prefix)
            prefix = ring_mul(prefix, {generator: 1}, group)
        else:
            prefix = ring_mul(prefix, {group.inv[generator]: 1}, group)
            if generator == letter:
                acc = ring_add(acc, {g: -c for g, c in prefix.items()})
    return acc


def syllable(name: str, exponent: int, s: int, t: int) -> Word:
    base = {"s": [(s, 1)], "t": [(t, 1)], "st": [(s, 1), (t, 1)], "ts": [(t, 1), (s, 1)]}[name]
    if exponent >= 0:
        return base * exponent
    return [(g, -e) for g, e in reversed(base)] * (-exponent)


def render(a: str, m: int, b: str, n: int) -> str:
    wrap = lambda name: name if len(name) == 1 else f"({name})"  # noqa: E731
    return f"{wrap(a)}^{m} {wrap(b)}^{n}"


def identify_relators(
    d2: RingMatrix, s: int, t: int, group: Group, span: int = 6
) -> list[list[str]]:
    """For each d2 row, every two-syllable word whose Fox jacobian equals it."""
    names = ("s", "t", "st", "ts")
    matches: list[list[str]] = [[] for _ in d2]
    for a in names:
        for b in names:
            if a == b:
                continue
            for m in range(-span, span + 1):
                for n in range(-span, span + 1):
                    if m == 0 or n == 0:
                        continue
                    word = syllable(a, m, s, t) + syllable(b, n, s, t)
                    jacobian = [fox_derivative(word, s, group), fox_derivative(word, t, group)]
                    for row, declared in enumerate(d2):
                        if jacobian[0] == declared[0] and jacobian[1] == declared[1]:
                            matches[row].append(render(a, m, b, n))
    return matches


# --------------------------------------------------------------------------------------
# the checks
# --------------------------------------------------------------------------------------


def numeric_leaves(value: Any) -> Iterator[Any]:
    if isinstance(value, dict):
        for item in value.values():
            yield from numeric_leaves(item)
    elif isinstance(value, list):
        for item in value:
            yield from numeric_leaves(item)
    elif value is not None and not isinstance(value, str):
        yield value


def canonical_bytes(packet: dict[str, Any]) -> str:
    return json.dumps(packet, sort_keys=True, indent=2, ensure_ascii=True) + "\n"


def sha256(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def run_checks(
    packet: dict[str, Any],
    raw: str,
    group_bytes: str,
    group: Group,
    reps: dict[str, list[CMatrix]],
    cover_primes: Sequence[int] = COVER_PRIMES,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    checks: list[dict[str, Any]] = []
    detail: dict[str, Any] = {}

    def record(identifier: str, description: str, passed: bool, observed: Any) -> None:
        checks.append(
            {"id": identifier, "check": description, "pass": bool(passed), "observed": observed}
        )

    ranks = packet.get("free_ranks", [])
    maps = packet.get("boundary_maps", {})
    try:
        d1, d2, d3 = (ring_matrix(maps[k]) for k in ("d1", "d2", "d3"))
    except (KeyError, TypeError, ValueError):
        record(
            "A1",
            "format: declared key set and well-formed boundary maps",
            False,
            "boundary_maps missing or malformed",
        )
        return checks, detail

    shapes = [(len(d1), len(d1[0])), (len(d2), len(d2[0])), (len(d3), len(d3[0]))]
    expected_shapes = (
        [(ranks[1], ranks[0]), (ranks[2], ranks[1]), (ranks[3], ranks[2])]
        if len(ranks) == 4
        else []
    )
    format_ok = (
        set(packet) == DECLARED_KEYS
        and packet.get("degree_range") == [0, 3]
        and packet.get("model_kind") in MODEL_KINDS
        and (packet.get("truncation_rule") is None)
        == (packet.get("model_kind") == "finite_cellular")
        and shapes == expected_shapes
        and all(isinstance(v, int) for v in numeric_leaves(packet))
    )
    record(
        "A1",
        "format: key set, degree range, truncation rule, shapes, integer leaves",
        format_ok,
        {
            "keys_ok": set(packet) == DECLARED_KEYS,
            "shapes": shapes,
            "expected_shapes": expected_shapes,
            "model_kind": packet.get("model_kind"),
        },
    )

    canonical = canonical_bytes(packet)
    canonical_ok = (
        canonical == raw
        and raw.isascii()
        and "\r" not in raw
        and raw.endswith("\n")
        and not raw.endswith("\n\n")
    )
    record(
        "A2",
        "canonicalization: delivered bytes are the declared canonical form",
        canonical_ok,
        {
            "canonical_matches_delivered": canonical == raw,
            "ascii": raw.isascii(),
            "crlf": "\r" in raw,
        },
    )

    s = packet["abstract_generators"].get("s")
    t = packet["abstract_generators"].get("t")
    in_range = (
        isinstance(s, int) and isinstance(t, int) and 0 <= s < group.size and 0 <= t < group.size
    )
    relations_ok = False
    generated = 0
    if in_range:
        central = group.power(s, 3)
        relations_ok = (
            central == group.power(t, 5) == group.power(group.mul[s][t], 2)
            and central != group.identity
            and group.power(central, 2) == group.identity
        )
        generated = group.generated_by([s, t])
    binding_ok = (
        packet.get("group_packet_sha256") == sha256(group_bytes)
        and in_range
        and relations_ok
        and generated == group.size
    )
    record(
        "A3",
        "group binding: hash, generator correspondence, relations, generation",
        binding_ok,
        {
            "group_hash_reproduces": packet.get("group_packet_sha256") == sha256(group_bytes),
            "s": s,
            "t": t,
            "orders": (
                [group.order(s), group.order(t), group.order(group.mul[s][t])]
                if in_range
                else None
            ),
            "relations_hold": relations_ok,
            "elements_generated": generated,
        },
    )
    detail["identity_element_id"] = group.identity

    complex_ok = is_zero_matrix(ring_matmul(d2, d1, group)) and is_zero_matrix(
        ring_matmul(d3, d2, group)
    )
    record(
        "A4",
        "the complex is a complex: d_(n+1) d_n = 0 over Z[2I]",
        complex_ok,
        {
            "d2_d1_zero": is_zero_matrix(ring_matmul(d2, d1, group)),
            "d3_d2_zero": is_zero_matrix(ring_matmul(d3, d2, group)),
        },
    )

    chi = ranks[0] - ranks[1] + ranks[2] - ranks[3] if len(ranks) == 4 else None
    record(
        "A5",
        "census: declared free ranks and chi = 0",
        chi == 0,
        {"free_ranks": ranks, "chi": chi},
    )

    e1, e2, e3 = augment(d1), augment(d2), augment(d3)
    determinant = e2[0][0] * e2[1][1] - e2[0][1] * e2[1][0] if len(e2) == 2 else None
    correspondence = [
        [ring([[1, s], [-1, group.identity]])],
        [ring([[1, t], [-1, group.identity]])],
    ]
    trivial_ok = (
        all(v == 0 for row in e1 for v in row)
        and determinant is not None
        and abs(determinant) == 1
        and all(v == 0 for row in e3 for v in row)
        and d1 == correspondence
    )
    record(
        "A6",
        "trivial sector: H_*(Z (x) C_*) = (Z, 0, 0, Z) through the declared eps",
        trivial_ok,
        {
            "eps_d1": e1,
            "eps_d2": e2,
            "eps_d2_det": determinant,
            "eps_d3": e3,
            "d1_is_the_1_cell_correspondence": d1 == correspondence,
        },
    )
    detail["augmented_homology"] = "(Z, 0, 0, Z)" if trivial_ok else "not (Z, 0, 0, Z)"

    cover_matrices = [regular_evaluate(m, group) for m in (d1, d2, d3)]
    betti = {p: cover_betti(cover_matrices, ranks, group.size, p) for p in cover_primes}
    cover_ok = all(b == (1, 0, 0, 1) for b in betti.values())
    record(
        "A7",
        "universal cover: H_*(C_* (x) F_p) = (1, 0, 0, 1) at every tested prime",
        cover_ok,
        {str(p): list(b) for p, b in betti.items()},
    )
    detail["cover_betti"] = {str(p): list(b) for p, b in betti.items()}
    detail["cover_ranks"] = {
        str(p): [rank_mod_p(m, p) for m in cover_matrices] for p in cover_primes
    }

    acyclicity: dict[str, bool] = {}
    for name, rep in reps.items():
        dim = len(rep[0])
        acyclicity[name] = all(exact_rank(evaluate(m, rep)) == dim for m in (d1, d2, d3))
    acyclic_ok = (not acyclicity["R0"]) and all(v for k, v in acyclicity.items() if k != "R0")
    record(
        "A8",
        "per-irrep acyclicity: trivial non-acyclic (expected), nontrivial acyclic",
        acyclic_ok,
        acyclicity,
    )

    prose = json.dumps(packet.get("basing", {})) + json.dumps(packet.get("top_closure", {}))
    scanned = raw.replace(prose.strip("{}"), "")
    hits = [word for word in LEAKAGE_VOCABULARY if word in scanned.lower()]
    leakage_ok = not hits and not re.search(r"\d\.\d", raw) and set(packet) <= DECLARED_KEYS
    record(
        "A9",
        "leakage: no answer-bearing vocabulary, no decimal literal, no extra key",
        leakage_ok,
        {
            "vocabulary_hits": hits,
            "decimal_literal": bool(re.search(r"\d\.\d", raw)),
            "extra_keys": sorted(set(packet) - DECLARED_KEYS),
        },
    )

    if in_range:
        matches = identify_relators(d2, s, t, group)
        basis_order = str(packet.get("basing", {}).get("basis_order", ""))
        unique = all(len(m) == 1 for m in matches)
        declared_ok = unique and all(m[0] in basis_order for m in matches)
    else:
        matches, unique, declared_ok = [], False, False
    record(
        "A10",
        "relator provenance: each d2 row is the Fox jacobian of its declared relator",
        declared_ok,
        {"unique_match_per_row": unique, "matched_relators": matches},
    )
    detail["matched_relators"] = matches

    if len(ranks) == 4:
        exact_ok, exact_observed = exact_top_boundary_check(
            cover_matrices, ranks, group.size, all(v == 0 for row in e1 for v in row)
        )
    else:
        exact_ok, exact_observed = False, {"outcome": "free_ranks malformed"}
    record(
        "A11",
        "exact top boundary: im d3 = ker d2 as lattices, certified without a prime set",
        exact_ok,
        exact_observed,
    )
    detail["exact_top_boundary"] = exact_observed

    return checks, detail


# --------------------------------------------------------------------------------------
# mutation suite
# --------------------------------------------------------------------------------------


def _set_coefficient(packet: dict[str, Any]) -> None:
    packet["boundary_maps"]["d3"][0][0][0][0] += 1


def _swap_entries(packet: dict[str, Any]) -> None:
    row = packet["boundary_maps"]["d3"][0]
    row[0], row[1] = row[1], row[0]


def _shift_element_id(packet: dict[str, Any]) -> None:
    term = packet["boundary_maps"]["d2"][0][0][0]
    term[1] = (term[1] + 1) % 120


def _flip_sign(packet: dict[str, Any]) -> None:
    term = packet["boundary_maps"]["d2"][1][1][2]
    term[0] = -term[0]


def _drop_augmentation_term(packet: dict[str, Any]) -> None:
    packet["boundary_maps"]["d1"][0][0].pop()


def _remap_generator(packet: dict[str, Any]) -> None:
    packet["abstract_generators"]["s"] = 20


def _inflate_rank(packet: dict[str, Any]) -> None:
    packet["free_ranks"] = [1, 2, 2, 2]


def _scale_top_boundary(packet: dict[str, Any]) -> None:
    packet["boundary_maps"]["d3"][0] = [
        [[2 * c, g] for c, g in entry] for entry in packet["boundary_maps"]["d3"][0]
    ]


def _relabel_relator(packet: dict[str, Any]) -> None:
    """Mislabel the FIRST relator, which is the one currently agreeing with its matrix.

    Mutating the second relator's label would not do: swapping `(st)` for `(ts)` there is
    the correction, not a defect, and A10 goes green under it.  That asymmetry is itself the
    causal demonstration that A10 reads the matrices rather than the prose.
    """
    packet["basing"]["basis_order"] = packet["basing"]["basis_order"].replace("s^3", "s^4")


def _scale_top_boundary_offlist(packet: dict[str, Any]) -> None:
    """Scale d3 by a prime OUTSIDE the A7 list, which is what makes A11 load-bearing.

    A7 sees `2 d3` only because 2 is in its prime list.  Multiply by 11 instead and every
    Betti number A7 computes is still (1, 0, 0, 1), at 2, 3, 5, 7 and 10007 alike, while
    every torsion value is multiplied by 11^(+-dim).  Only A11 reddens here, and that is
    the point of it: no finite prime set is a sufficient accept criterion.
    """
    packet["boundary_maps"]["d3"][0] = [
        [[11 * c, g] for c, g in entry] for entry in packet["boundary_maps"]["d3"][0]
    ]


MUTATIONS: tuple[tuple[str, Callable[[dict[str, Any]], None], str], ...] = (
    ("d3 coefficient off by one", _set_coefficient, "A4"),
    ("d3 two entries swapped", _swap_entries, "A4"),
    ("d2 element id shifted", _shift_element_id, "A4"),
    ("d2 sign flipped on one term", _flip_sign, "A4"),
    ("d1 augmentation term dropped", _drop_augmentation_term, "A6"),
    ("abstract generator s remapped", _remap_generator, "A3"),
    ("free_ranks top rank inflated", _inflate_rank, "A1"),
    ("d3 scaled to 2 d3 inside ker d2", _scale_top_boundary, "A7"),
    ("declared relator relabelled", _relabel_relator, "A10"),
    ("d3 scaled to 11 d3, invisible to A7", _scale_top_boundary_offlist, "A11"),
)


def run_mutation_suite(
    packet: dict[str, Any], group_bytes: str, group: Group, reps: dict[str, list[CMatrix]]
) -> list[dict[str, Any]]:
    results: list[dict[str, Any]] = []
    for label, mutate, expected in MUTATIONS:
        mutated = copy.deepcopy(packet)
        mutate(mutated)
        raw = canonical_bytes(mutated)
        checks, _ = run_checks(
            mutated, raw, group_bytes, group, reps, cover_primes=MUTATION_PRIMES
        )
        red = sorted(c["id"] for c in checks if not c["pass"])
        results.append(
            {
                "mutation": label,
                "expected_red": expected,
                "observed_red": red,
                "detected": expected in red,
            }
        )
    return results


# --------------------------------------------------------------------------------------
# entry point
# --------------------------------------------------------------------------------------


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--packet", type=Path, default=DEFAULT_PACKET)
    parser.add_argument("--group-packet", type=Path, default=DEFAULT_GROUP_PACKET)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--mutation-tests", action="store_true")
    parser.add_argument("--no-write", action="store_true")
    args = parser.parse_args()

    for path in (args.packet, args.group_packet):
        if not path.exists():
            print(f"packet not found: {path}")
            print("the construction packet lands with the M8.8 protocol; see the task record")
            return 2

    raw = args.packet.read_text()
    packet = json.loads(raw)
    group_bytes = args.group_packet.read_text()
    group_packet = json.loads(group_bytes)

    generators = [tuple(parse_component(c) for c in g) for g in group_packet["generators"]]
    group = Group(generators)
    if group.size != 120:
        print(f"group packet closure is {group.size} elements, expected 120")
        return 2
    reps = build_irreps(group)

    checks, detail = run_checks(packet, raw, group_bytes, group, reps)
    all_pass = all(c["pass"] for c in checks)
    mutations = (
        run_mutation_suite(packet, group_bytes, group, reps) if args.mutation_tests else None
    )
    mutations_ok = all(m["detected"] for m in mutations) if mutations else None

    incoming = sha256(raw)
    authoritative = sha256(canonical_bytes(packet))
    result = {
        "what": "M8.8 construction-packet audit, run from the two public packets alone",
        "packet_file": args.packet.name,
        "group_packet_file": args.group_packet.name,
        "incoming_sha256": incoming,
        "canonical_form": "json.dumps(sort_keys=True, indent=2, ensure_ascii=True) + LF",
        "authoritative_sha256": authoritative,
        "canonicalization_changed_bytes": incoming != authoritative,
        "checks": checks,
        "detail": detail,
        "all_checks_pass": all_pass,
        "mutation_suite": mutations,
        "mutation_suite_all_detected": mutations_ok,
    }
    if not args.no_write:
        args.out.write_text(json.dumps(result, indent=2) + "\n")

    print(f"packet            {args.packet}")
    print(f"incoming sha256   {incoming}")
    print(f"authoritative     {authoritative}")
    changed = "CHANGED the bytes" if incoming != authoritative else "no-op (already canonical)"
    print(f"canonicalization  {changed}")
    print()
    for check in checks:
        print(f"  {'PASS' if check['pass'] else 'FAIL'}  {check['id']:4s} {check['check']}")
        print(f"        {json.dumps(check['observed'])}")
    print()
    print(f"  universal-cover Betti     {json.dumps(detail.get('cover_betti'))}")
    print(f"  matched relators          {json.dumps(detail.get('matched_relators'))}")
    if mutations is not None:
        print()
        for mutation in mutations:
            flag = "DETECTED" if mutation["detected"] else "MISSED  "
            print(
                f"  {flag}  {mutation['mutation']:34s} expected red {mutation['expected_red']},"
                f" observed {mutation['observed_red']}"
            )
    print()
    print(
        f"ALL CHECKS PASS: {all_pass}"
        + (f" | MUTATIONS ALL DETECTED: {mutations_ok}" if mutations is not None else "")
    )
    if not args.no_write:
        print(f"written: {args.out}")

    return 0 if all_pass and (mutations_ok is not False) else 1


if __name__ == "__main__":
    sys.exit(main())
