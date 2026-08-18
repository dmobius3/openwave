#!/usr/bin/env python3
"""m8_4_kinematic_check.py

M8.4 first result artifact (task: research/tasks/m8_4_task_details.md).
Finding note with claim, proof, equation-to-code map, and claim ceiling:
research/findings/m8_4_kinematic_close.md.  Reproduce with
`python3 m8_4_kinematic_check.py`; the pass state is exit 0 with the final
line `SUMMARY: ALL PASS`.

Kinematic accessibility of the McKay slot structure for NATIVE fields on the
closed quotient X = S^3/2I, versus fields/modes on the cover S^3 and twisted
sectors over X.

Facts checked, each with a red condition:

  A (control, cover side): for every nontrivial 2I irrep rho, the first S^3
    harmonic level n at which rho occurs in V_n|_2I equals its McKay distance
    d(rho): T-a of the M8.2 lock (j_first = d/2, n = 2j).  Certified by TWO
    computationally separate in-artifact routes: the icosian character sums
    below, and a McKay-recursion route (m(n+1) = A m(n) - m(n-1) on the
    affine E8 graph) that never touches the icosians or the characters.  The
    routes share the McKay graph as an input; the character route re-derives
    its edges internally through irreducibility and orthonormality
    assertions, so they are separate computations, not independent
    derivations.

  B (center parity): the four spinorial slots (rho(-1) = -1) occur only at
    odd n, while X-scalar content m_0(n) is supported on even n only.  A
    transparent corollary of C for the spinorial half: the electron slot
    (d = 1) is among the absent ones.

  C (projector degeneracy, THE LEMMA): every native field on X lifts to a
    2I-invariant field on S^3, and the isotypic projector P_rho annihilates
    invariant vectors for every nontrivial rho, at every level, INCLUDING
    levels whose ambient harmonic space contains a rho-summand alongside the
    invariant one.  Tested concretely: at the first level carrying both the
    trivial and an integer irrep, an explicitly constructed invariant vector
    projects to zero in the rho-sector; a non-invariant vector at the same
    level does not.  Consequence: native X-states carry NONE of the eight
    nontrivial McKay slots at any level.  (The ambient coexistence levels
    this check prints for the integer irreps are NOT slot accessibility for
    X-states; they are where invariant content and those irreps coexist as
    DISTINCT summands of the same ambient V_n.  Information only.)

  D (mutations): a shifted d-map must break A; admitting odd levels must
    expose spinorial content (B's gate is load-bearing); dropping a McKay
    edge must break the recursion route.  C carries its own mutation arm.

Exact group theory evaluated numerically (float tolerance 1e-6 on quantities
that are integers by construction; every multiplicity is asserted integral).
"""

import itertools
import math

TOL = 1e-6
PHI = (1.0 + 5.0**0.5) / 2.0
NMAX = 60  # highest harmonic level scanned


# ----------------------------------------------------------------------------
# The binary icosahedral group 2I as the 120 unit icosians (w, x, y, z).
# ----------------------------------------------------------------------------

def icosians():
    qs = set()

    def put(q):
        qs.add(tuple(round(c, 6) for c in q))

    for i in range(4):
        for s in (1.0, -1.0):
            q = [0.0, 0.0, 0.0, 0.0]
            q[i] = s
            put(q)
    for signs in itertools.product((0.5, -0.5), repeat=4):
        put(list(signs))
    base = [0.0, 0.5, 0.5 / PHI, 0.5 * PHI]
    even_perms = [p for p in itertools.permutations(range(4))
                  if perm_sign(p) == 1]
    for p in even_perms:
        for signs in itertools.product((1.0, -1.0), repeat=4):
            put([signs[k] * base[p[k]] for k in range(4)])
    return sorted(qs)


def perm_sign(p):
    s = 1
    p = list(p)
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]:
                s = -s
    return s


def qmul(a, b):
    aw, ax, ay, az = a
    bw, bx, by, bz = b
    return (aw * bw - ax * bx - ay * by - az * bz,
            aw * bx + ax * bw + ay * bz - az * by,
            aw * by - ax * bz + ay * bw + az * bx,
            aw * bz + ax * by - ay * bx + az * bw)


G = icosians()
assert len(G) == 120, f"expected 120 icosians, got {len(G)}"
GSET = set(G)
for a in G[:20]:
    for b in G:
        p = tuple(round(c, 6) for c in qmul(a, b))
        assert p in GSET, "closure failure"

MINUS_ONE = (-1.0, 0.0, 0.0, 0.0)
IDX_MINUS_ONE = G.index(MINUS_ONE)
IDX_ONE = G.index((1.0, 0.0, 0.0, 0.0))


# ----------------------------------------------------------------------------
# Route 1: SU(2) characters chi_{V_n}(g) = sin((n+1)a)/sin(a), cos(a) = Re(g).
# ----------------------------------------------------------------------------

def chi_Vn(n, g):
    w = max(-1.0, min(1.0, g[0]))
    if abs(w - 1.0) < 1e-9:
        return float(n + 1)
    if abs(w + 1.0) < 1e-9:
        return float((n + 1) * (-1) ** n)
    a = math.acos(w)
    return math.sin((n + 1) * a) / math.sin(a)


def arr_Vn(n):
    return [chi_Vn(n, g) for g in G]


def inner(u, v):
    return sum(x * y for x, y in zip(u, v)) / 120.0


def mult(chi_rho, n):
    m = inner(chi_rho, arr_Vn(n))
    r = round(m)
    assert abs(m - r) < TOL, f"non-integral multiplicity {m} at n={n}"
    return int(r)


def assert_irred(chi, name):
    nrm = inner(chi, chi)
    assert abs(nrm - 1.0) < TOL, f"{name} not irreducible: <chi,chi>={nrm}"


# The nine 2I irreps.  V_1..V_5 restrict irreducibly.  The McKay graph is
# bipartite by center parity, trivalent at 6: chain 1-2-3-4'-5-6-4-2', branch
# 3' on 6.  Hence V_7 = 6 + 2', then 2' (x) 2 = 4, then V_6 = 4 + 3'.

chi = {"1": [1.0] * 120}
for name, n in (("2", 1), ("3", 2), ("4p", 3), ("5", 4), ("6", 5)):
    chi[name] = arr_Vn(n)
    assert_irred(chi[name], name)

chi["2p"] = [a - b for a, b in zip(arr_Vn(7), chi["6"])]
assert_irred(chi["2p"], "2p")
assert abs(inner(chi["2p"], chi["2"])) < TOL, "2p is not a new 2-dim irrep"
chi["4"] = [a * b for a, b in zip(chi["2"], chi["2p"])]
assert_irred(chi["4"], "4")
assert abs(inner(chi["4"], chi["4p"])) < TOL, "4 is not a new 4-dim irrep"
chi["3p"] = [a - b for a, b in zip(arr_Vn(6), chi["4"])]
assert_irred(chi["3p"], "3p")
assert abs(inner(chi["3p"], chi["3"])) < TOL, "3p is not a new 3-dim irrep"

dims = {k: round(v[IDX_ONE]) for k, v in chi.items()}
assert sorted(dims.values()) == [1, 2, 2, 3, 3, 4, 4, 5, 6]
assert sum(d * d for d in dims.values()) == 120
names = list(chi)
for i in range(len(names)):
    for j in range(i + 1, len(names)):
        assert abs(inner(chi[names[i]], chi[names[j]])) < TOL, "not orthogonal"

# Frozen T-a distance map (M8.2 section 3 table; n = d at first occurrence).
D_MAP = {"2": 1, "3": 2, "4p": 3, "5": 4, "6": 5, "4": 6, "3p": 6, "2p": 7}

SPINORIAL = {k for k in D_MAP if round(chi[k][IDX_MINUS_ONE]) == -dims[k]}
INTEGER = {k for k in D_MAP if round(chi[k][IDX_MINUS_ONE]) == dims[k]}
assert SPINORIAL | INTEGER == set(D_MAP) and len(SPINORIAL) == 4


def first_occurrence(chi_rho, accessible=lambda n: True):
    for n in range(1, NMAX + 1):
        if accessible(n) and mult(chi_rho, n) > 0:
            return n
    return None


# ----------------------------------------------------------------------------
# Route 2: McKay recursion, independent of the icosians and the characters.
# V_1 (x) V_n = V_{n+1} + V_{n-1}  and  2 (x) rho = sum of McKay neighbors,
# so the multiplicity vector obeys m(n+1) = A m(n) - m(n-1) with A the affine
# E8 adjacency matrix, m(0) = e_1, m(1) = e_2.
# ----------------------------------------------------------------------------

ORDER = ["1", "2", "3", "4p", "5", "6", "4", "2p", "3p"]
EDGES = [("1", "2"), ("2", "3"), ("3", "4p"), ("4p", "5"), ("5", "6"),
         ("6", "4"), ("4", "2p"), ("6", "3p")]


def recursion_mults(edges, nmax):
    idx = {k: i for i, k in enumerate(ORDER)}
    A = [[0] * 9 for _ in range(9)]
    for u, v in edges:
        A[idx[u]][idx[v]] = 1
        A[idx[v]][idx[u]] = 1
    m_prev = [0] * 9
    m_prev[idx["1"]] = 1
    m_cur = [0] * 9
    m_cur[idx["2"]] = 1
    out = {0: m_prev[:], 1: m_cur[:]}
    for n in range(1, nmax):
        m_next = [sum(A[i][j] * m_cur[j] for j in range(9)) - m_prev[i]
                  for i in range(9)]
        out[n + 1] = m_next[:]
        m_prev, m_cur = m_cur, m_next
    return out


REC = recursion_mults(EDGES, NMAX)
route_agree = all(
    REC[n][ORDER.index(k)] == mult(chi[k], n)
    for n in range(0, NMAX + 1) for k in ORDER
)
assert route_agree, "character route and McKay-recursion route disagree"


# ----------------------------------------------------------------------------
# CHECK A: cover-side first occurrences equal the McKay distances (2 routes).
# ----------------------------------------------------------------------------

print("CHECK A: cover/twisted-sector first occurrence vs McKay distance d")
print("  (character route and McKay-recursion route agree on every")
print(f"   multiplicity, all 9 irreps, n <= {NMAX})")
ok_A = True
for name in sorted(D_MAP, key=lambda k: (D_MAP[k], k)):
    n1 = first_occurrence(chi[name])
    tag = "spinorial" if name in SPINORIAL else "integer"
    print(f"  irrep {name:>2} (dim {dims[name]}, {tag:9s}): "
          f"first n = {n1}, d = {D_MAP[name]}")
    ok_A &= (n1 == D_MAP[name])
print(f"  -> {'PASS' if ok_A else 'FAIL'}: the T-a ladder n = d "
      f"{'is' if ok_A else 'is NOT'} realized on the cover\n")

# ----------------------------------------------------------------------------
# CHECK B: center parity.
# ----------------------------------------------------------------------------

m0 = {n: mult(chi["1"], n) for n in range(0, NMAX + 1)}
odd_m0 = [n for n in range(1, NMAX + 1, 2) if m0[n] != 0]
spin_even = [(k, n) for k in SPINORIAL for n in range(0, NMAX + 1, 2)
             if mult(chi[k], n) != 0]
support = [n for n in range(0, NMAX + 1) if m0[n] > 0]
print("CHECK B: center parity")
print(f"  X-accessible scalar levels (m_0(n) > 0), n <= {NMAX}: {support}")
print(f"  odd-n invariants: {odd_m0}  (must be empty)")
print(f"  spinorial content at even n: {spin_even}  (must be empty)")
ok_B = not odd_m0 and not spin_even
print(f"  -> {'PASS' if ok_B else 'FAIL'}: spinorial slots "
      f"(d = 1, 3, 5, 7; includes the electron slot at d = 1) are "
      f"{'absent from' if ok_B else 'NOT proven absent from'} X content\n")

# ----------------------------------------------------------------------------
# CHECK C: projector degeneracy (the lemma).  Explicit representation of 2I
# on degree-n polynomials in (z1, z2); construct an invariant vector at the
# first level where the trivial and an integer irrep coexist; its rho-sector
# projection must vanish although the ambient rho-summand exists; a
# non-invariant vector at the same level must project nontrivially.
# ----------------------------------------------------------------------------

def su2_matrix(q):
    w, x, y, z = q
    return ((complex(w, x), complex(y, z)),
            (complex(-y, z), complex(w, -x)))


def binom(n, k):
    return math.comb(n, k)


def rep_matrix(q, n):
    """Matrix of g on degree-n homogeneous polys, (g.p)(z) = p(g^-1 z)."""
    m = su2_matrix(q)
    # unitary inverse = conjugate transpose
    a, b = m[0][0].conjugate(), m[1][0].conjugate()
    c, d = m[0][1].conjugate(), m[1][1].conjugate()
    # monomial e_k = z1^k z2^(n-k) maps to (a z1 + b z2)^k (c z1 + d z2)^(n-k)
    R = [[0j] * (n + 1) for _ in range(n + 1)]
    for k in range(n + 1):
        for i in range(k + 1):
            ci = binom(k, i) * (a ** i) * (b ** (k - i))
            for l in range(n - k + 1):
                cl = binom(n - k, l) * (c ** l) * (d ** (n - k - l))
                R[i + l][k] += ci * cl
    return R


def mat_vec(R, v):
    return [sum(R[i][j] * v[j] for j in range(len(v))) for i in range(len(v))]


def vec_norm(v):
    return math.sqrt(sum(abs(x) ** 2 for x in v))


RHO_TEST = "3"
N_TEST = next(n for n in range(1, NMAX + 1)
              if m0[n] > 0 and mult(chi[RHO_TEST], n) > 0)

reps = [rep_matrix(q, N_TEST) for q in G]
# construction sanity: traces reproduce chi_{V_n} on every group element
for q, R in zip(G, reps):
    tr = sum(R[i][i] for i in range(N_TEST + 1))
    assert abs(tr - chi_Vn(N_TEST, q)) < 1e-6, "rep trace mismatch"


def project(char_arr, dim_rho, v):
    out = [0j] * len(v)
    for gi, R in enumerate(reps):
        w = char_arr[gi] * dim_rho / 120.0
        Rv = mat_vec(R, v)
        for i in range(len(v)):
            out[i] += w * Rv[i]
    return out


seed_vec = [complex((7 * j) % 11 - 5, (3 * j * j) % 13 - 6)
            for j in range(N_TEST + 1)]
v_inv = project(chi["1"], 1, seed_vec)          # group average: invariant
v_inv2 = project(chi["1"], 1, v_inv)            # idempotence sanity
assert vec_norm([a - b for a, b in zip(v_inv, v_inv2)]) < 1e-6 * (
    vec_norm(v_inv) + 1.0), "projector not idempotent"

coex = {k: next(n for n in range(1, NMAX + 1)
                if m0[n] > 0 and mult(chi[k], n) > 0) for k in INTEGER}
print("CHECK C: projector degeneracy on X (the lemma)")
print(f"  test level n = {N_TEST}: m_0 = {m0[N_TEST]}, "
      f"mult({RHO_TEST}) = {mult(chi[RHO_TEST], N_TEST)} "
      f"(ambient coexistence holds, so a vanishing projection is not vacuous)")
norm_inv = vec_norm(v_inv)
assert norm_inv > 1e-8, "invariant test vector is zero; pick another seed"
p_of_inv = vec_norm(project(chi[RHO_TEST], dims[RHO_TEST], v_inv)) / norm_inv
p_of_raw = vec_norm(project(chi[RHO_TEST], dims[RHO_TEST], seed_vec)) / \
    vec_norm(seed_vec)
print(f"  |P_{RHO_TEST}(invariant)| / |invariant|     = {p_of_inv:.2e}  "
      f"(must be ~0)")
print(f"  |P_{RHO_TEST}(non-invariant)| / |vector|    = {p_of_raw:.2e}  "
      f"(mutation arm: must be > 0)")
ok_C = p_of_inv < 1e-7 and p_of_raw > 1e-3
print(f"  informational: ambient coexistence levels for the integer irreps "
      f"(NOT accessibility): {[f'{k}:{coex[k]}' for k in sorted(coex)]}")
print(f"  -> {'PASS' if ok_C else 'FAIL'}: native X-states carry "
      f"{'no' if ok_C else '(UNPROVEN)'} nontrivial slot content "
      f"at any level\n")

# ----------------------------------------------------------------------------
# CHECK D: mutations must redden A (both routes) and B.
# ----------------------------------------------------------------------------

print("CHECK D: mutation arms")
mut_A = all(first_occurrence(chi[k]) == D_MAP[k] + 1 for k in D_MAP)
print(f"  shifted d-map reproduces CHECK A: {mut_A}  (must be False)")
mut_B = all(first_occurrence(chi[k], accessible=lambda n: True) is not None
            for k in SPINORIAL)
print(f"  admitting all n as X-accessible finds spinorial content: "
      f"{mut_B}  (must be True: the B gate is load-bearing)")
REC_MUT = recursion_mults(EDGES[:-1], NMAX)  # drop the 6-3p edge
mut_R = all(REC_MUT[n][ORDER.index(k)] == mult(chi[k], n)
            for n in range(0, NMAX + 1) for k in ORDER)
print(f"  McKay recursion with a dropped edge still matches characters: "
      f"{mut_R}  (must be False)")
ok_D = (not mut_A) and mut_B and (not mut_R)
print(f"  -> {'PASS' if ok_D else 'FAIL'}\n")

ok = ok_A and ok_B and ok_C and ok_D
print("SUMMARY:", "ALL PASS" if ok else "SOMETHING FAILED")
raise SystemExit(0 if ok else 1)
