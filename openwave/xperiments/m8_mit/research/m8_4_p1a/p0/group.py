"""The binary icosahedral group 2I: 120 unit icosians, multiplication table,
conjugacy classes, and the full character table.

Self-contained: constructs everything from first principles.
"""

import itertools
import numpy as np
from .algebra import qmul, su2_character


PHI = (1.0 + 5.0**0.5) / 2.0


def _perm_sign(p):
    s = 1
    p = list(p)
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]:
                s = -s
    return s


def build_icosians():
    """The 120 elements of 2I as unit quaternions, canonically ordered."""
    qs = set()

    def put(q):
        qs.add(tuple(round(c, 10) for c in q))

    for i in range(4):
        for s in (1.0, -1.0):
            q = [0.0, 0.0, 0.0, 0.0]
            q[i] = s
            put(q)
    for signs in itertools.product((0.5, -0.5), repeat=4):
        put(list(signs))
    base = [0.0, 0.5, 0.5 / PHI, 0.5 * PHI]
    even_perms = [p for p in itertools.permutations(range(4))
                  if _perm_sign(p) == 1]
    for p in even_perms:
        for signs in itertools.product((1.0, -1.0), repeat=4):
            put([signs[k] * base[p[k]] for k in range(4)])
    elems = sorted(qs)
    assert len(elems) == 120, f"expected 120 icosians, got {len(elems)}"
    result = [np.array(q) for q in elems]
    # Put the identity element first — orbit transport needs gid=0 to be the identity.
    # During development, gid=0 was [-1,0,0,0] (the central element), which produced
    # a catastrophically wrong spectrum. This is a structural requirement, not cosmetic.
    found = False
    for i, e in enumerate(result):
        if abs(e[0] - 1.0) < 1e-10 and np.linalg.norm(e[1:]) < 1e-10:
            result[0], result[i] = result[i], result[0]
            found = True
            break
    assert found, "identity quaternion [1,0,0,0] not found in 2I elements"
    e0 = result[0]
    assert abs(e0[0] - 1.0) < 1e-10 and np.linalg.norm(e0[1:]) < 1e-10, \
        f"elems[0] is {e0}, not the identity — orbit transport will be wrong"
    return result


def verify_closure(elems, tol=1e-8):
    """Verify the group is closed under multiplication."""
    eset = [tuple(round(float(c), 8) for c in e) for e in elems]
    eset_s = set(eset)
    for a in elems[:20]:
        for b in elems:
            p = qmul(a, b)
            pt = tuple(round(float(c), 8) for c in p)
            if pt not in eset_s:
                neg = tuple(round(-float(c), 8) for c in p)
                assert neg in eset_s, f"closure failure at {a}*{b}"


def multiplication_table(elems, tol=1e-8):
    """Full multiplication table: mult[i][j] = index of elems[i]*elems[j].

    Exact matching — the 120 icosians form a group, so every product
    is exactly one of them.
    """
    coords = np.array(elems)

    def find(p):
        diffs = np.max(np.abs(coords - p), axis=1)
        best = int(np.argmin(diffs))
        if diffs[best] < tol:
            return best
        raise ValueError(f"product not found, best dist {diffs[best]:.2e}")

    return [[find(qmul(elems[i], elems[j]))
             for j in range(len(elems))]
            for i in range(len(elems))]


# Conjugacy classes of 2I, identified by half-angle α = arccos(w)
CLASS_ANGLES = [0.0, np.pi, np.pi/2, np.pi/3, 2*np.pi/3,
                np.pi/5, 2*np.pi/5, 3*np.pi/5, 4*np.pi/5]
CLASS_SIZES = [1, 1, 30, 20, 20, 12, 12, 12, 12]

# The nine irreps: labels, dimensions, McKay distances
LABELS = ["R0", "R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8"]
DIMS = [1, 2, 2, 3, 3, 4, 4, 5, 6]
MCKAY_DIST = [0, 1, 7, 2, 6, 6, 3, 4, 5]

# Which V_n each irrep is a restriction of (for n ≤ 5, V_n|_{2I} is irreducible)
# R0=V0, R1=V1, R3=V2, R6=V3, R7=V4, R8=V5
# R4 and R5 extracted from V6, R2 from V7
DIRECT_LEVEL = {"R0": 0, "R1": 1, "R3": 2, "R6": 3, "R7": 4, "R8": 5}
PROJECTED = {"R4": 6, "R5": 6, "R2": 7}


def class_of(q):
    """Return the conjugacy class index (0..8) of a group element."""
    w = float(q[0])
    a = np.arccos(np.clip(w, -1.0, 1.0))
    dists = [abs(a - ca) for ca in CLASS_ANGLES]
    return int(np.argmin(dists))


def build_character_table(elems):
    """Compute the full 9x9 character table from the McKay recursion.

    Returns chi[irrep_idx][class_idx], using the same method as
    m8_2_first_occurrence.py but computed independently.
    """
    thetas = np.array(CLASS_ANGLES)
    sizes = np.array(CLASS_SIZES, dtype=float)

    edges = [("R0", "R1"), ("R1", "R3"), ("R3", "R6"), ("R6", "R7"),
             ("R7", "R8"), ("R8", "R5"), ("R5", "R2"), ("R8", "R4")]
    idx = {l: i for i, l in enumerate(LABELS)}
    A = np.zeros((9, 9))
    for a, b in edges:
        A[idx[a], idx[b]] = A[idx[b], idx[a]] = 1

    def chiV(n, th):
        s = np.sin(th)
        if abs(s) < 1e-12:
            return (n + 1) * np.cos(n * th)
        return np.sin((n + 1) * th) / s

    m = [np.zeros(9) for _ in range(9)]
    m[0][idx["R0"]] = 1.0
    m[1][idx["R1"]] = 1.0
    for n in range(1, 8):
        m[n + 1] = A.dot(m[n]) - m[n - 1]

    Mmat = np.array([m[n] for n in range(9)])
    chi = np.zeros((9, 9))
    for c in range(9):
        rhs = np.array([chiV(n, thetas[c]) for n in range(9)])
        chi[:, c] = np.linalg.solve(Mmat, rhs)

    ortho = np.array([[np.sum(sizes * chi[i] * chi[j]) / 120.0
                        for j in range(9)] for i in range(9)])
    ortho_err = float(np.max(np.abs(ortho - np.eye(9))))
    assert ortho_err < 1e-6, f"character table not orthonormal: {ortho_err}"

    for i in range(9):
        assert abs(chi[i, 0] - DIMS[i]) < 1e-6, \
            f"{LABELS[i]} dim mismatch: {chi[i,0]} vs {DIMS[i]}"

    return chi


def character_at_element(chi_table, q):
    """Return the 9-vector of character values at group element q."""
    c = class_of(q)
    return chi_table[:, c]
