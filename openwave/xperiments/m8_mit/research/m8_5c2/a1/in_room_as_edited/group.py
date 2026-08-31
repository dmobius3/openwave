"""Binary icosahedral group 2I: 120 unit quaternions, explicit construction,
character table, generating pair, and unitary representations for all 9 sectors.

Representations R1, R3, R6, R7, R8 are unitarized symmetric powers of the
SU(2) fundamental. R2 is extracted from V_7 via the Reynolds projector (first
occurrence of R2 is at level 7). R4 = Sym^2(R2) unitarized. R5 = R1 x R2.

All constructions match right_translation_check.py and mode_count.py, verified
at import time.
"""
import numpy as np
import itertools
from math import factorial
from scipy.linalg import sqrtm

TAU = (1 + np.sqrt(5)) / 2
TAU_INV = (np.sqrt(5) - 1) / 2

# ---- quaternion arithmetic ----
def qmul(a, b):
    w1, x1, y1, z1 = a
    w2, x2, y2, z2 = b
    return np.array([
        w1*w2 - x1*x2 - y1*y2 - z1*z2,
        w1*x2 + x1*w2 + y1*z2 - z1*y2,
        w1*y2 - x1*z2 + y1*w2 + z1*x2,
        w1*z2 + x1*y2 - y1*x2 + z1*w2,
    ])

def quat_to_su2(q):
    w, x, y, z = q
    return np.array([[w + 1j*x, y + 1j*z],
                     [-y + 1j*z, w - 1j*x]], dtype=complex)

def sym_power(M, n):
    dim = n + 1
    out = np.zeros((dim, dim), dtype=complex)
    for k in range(dim):
        polyA = np.array([1.0 + 0j])
        for _ in range(n - k):
            polyA = np.convolve(polyA, np.array([M[0, 0], M[1, 0]]))
        polyB = np.array([1.0 + 0j])
        for _ in range(k):
            polyB = np.convolve(polyB, np.array([M[0, 1], M[1, 1]]))
        out[:, k] = np.convolve(polyA, polyB)
    return out

# ---- group construction (matches right_translation_check.py) ----
def _build_icosians():
    Q = []
    for i in range(4):
        for s in (1.0, -1.0):
            q = np.zeros(4); q[i] = s; Q.append(q)
    for signs in range(16):
        Q.append(np.array([0.5 if (signs >> k) & 1 == 0 else -0.5 for k in range(4)]))
    base = np.array([0.0, 1.0, 1/TAU, TAU]) / 2
    evens = [p for p in itertools.permutations(range(4))
             if sum(1 for a in range(4) for b in range(a+1, 4) if p[a] > p[b]) % 2 == 0]
    for p in evens:
        for signs in range(8):
            sb = np.array([1.0, (-1)**(signs & 1), (-1)**((signs >> 1) & 1),
                           (-1)**((signs >> 2) & 1)])
            v = base * sb
            Q.append(np.array([v[p.index(k)] for k in range(4)]))
    return np.array(Q)

G120 = _build_icosians()
assert G120.shape == (120, 4)

def order_of(q):
    p = q.copy()
    for k in range(1, 12):
        if np.abs(p - np.array([1, 0, 0, 0])).sum() < 1e-9:
            return k
        p = qmul(p, q)
    return -1

# ---- verification ----
_unique_count = len({tuple(np.round(q, 9)) for q in G120})
assert _unique_count == 120
assert np.allclose(np.einsum('ij,ij->i', G120, G120), 1.0)
_order_census = {}
for q in G120:
    o = order_of(q)
    _order_census[o] = _order_census.get(o, 0) + 1
assert _order_census == {1: 1, 2: 1, 3: 20, 4: 30, 5: 24, 6: 20, 10: 24}

def _closure_residual(G):
    maxres = 0.0
    for i in range(len(G)):
        for j in range(len(G)):
            p = qmul(G[i], G[j])
            dists = np.sqrt(((G - p)**2).sum(axis=1))
            maxres = max(maxres, dists.min())
    return maxres

_cres = _closure_residual(G120)
assert _cres < 1e-9

# ---- generating pair ----
def _find_generators():
    ord10s = [i for i, q in enumerate(G120) if order_of(q) == 10]
    ord4s = [i for i, q in enumerate(G120) if order_of(q) == 4]
    def generate(seed_idx):
        S = {tuple(np.round(G120[i], 9)) for i in seed_idx}
        frontier = [G120[i] for i in seed_idx]
        while frontier:
            nxt = []
            for a in frontier:
                for i in seed_idx:
                    p = qmul(a, G120[i])
                    tp = tuple(np.round(p, 9))
                    if tp not in S:
                        S.add(tp)
                        nxt.append(p)
            frontier = nxt
        return S
    for a in ord10s[:3]:
        for b in ord4s:
            if len(generate([a, b])) == 120:
                return [a, b]
    raise RuntimeError("no generating pair found")

GEN_INDICES = _find_generators()

# ---- conjugacy class snapping ----
CLASS_ANGLES = np.array([0, np.pi, np.pi/2, np.pi/3, 2*np.pi/3,
                         np.pi/5, 2*np.pi/5, 3*np.pi/5, 4*np.pi/5])
CLASS_SIZES = np.array([1, 1, 30, 20, 20, 12, 12, 12, 12])

def snap_class(q):
    th = np.arccos(np.clip(q[0], -1, 1))
    k = int(np.argmin(np.abs(CLASS_ANGLES - th)))
    assert abs(CLASS_ANGLES[k] - th) < 1e-9, f"angle {th} off every 2I class"
    return k

ELEMENT_CLASSES = np.array([snap_class(q) for q in G120])

# ---- character table ----
GALOIS_CLASS = {5: 7, 7: 5, 6: 8, 8: 6}

def chiV(n, k):
    th = CLASS_ANGLES[k]
    if k == 0:
        return float(n + 1)
    if k == 1:
        return float((n + 1) * (-1)**n)
    return np.sin((n + 1) * th) / np.sin(th)

IRREP_NAMES = ["R0", "R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8"]

def _build_chartab():
    IRREPS8 = {
        "R0": lambda k: 1.0,
        "R1": lambda k: chiV(1, k),
        "R2": lambda k: chiV(1, GALOIS_CLASS.get(k, k)),
        "R3": lambda k: chiV(2, k),
        "R4": lambda k: chiV(2, GALOIS_CLASS.get(k, k)),
        "R6": lambda k: chiV(3, k),
        "R7": lambda k: chiV(4, k),
        "R8": lambda k: chiV(5, k),
    }
    others_vec = {r: np.array([IRREPS8[r](k) for k in ELEMENT_CLASSES]) for r in IRREPS8}
    import collections
    csize = collections.Counter(ELEMENT_CLASSES)
    cent = {k: 120 // csize[k] for k in csize}
    absq = {k: cent[k] - sum(IRREPS8[r](k)**2 for r in IRREPS8) for k in csize}
    amb = [k for k in sorted(csize) if absq[k] > 0.5 and k != 0]
    hits = []
    for bits in range(2 ** len(amb)):
        row = {k: float(np.sqrt(max(absq[k], 0.0))) for k in csize}
        for j, k in enumerate(amb):
            if (bits >> j) & 1:
                row[k] = -row[k]
        vec = np.array([row[k] for k in ELEMENT_CLASSES])
        if (all(abs(np.mean(vec * others_vec[r])) < 1e-9 for r in IRREPS8)
                and abs(np.mean(vec * vec) - 1.0) < 1e-9):
            hits.append(vec)
    assert len(hits) == 1, f"R5 sign solution not unique: {len(hits)}"
    tab = {r: others_vec[r] for r in IRREPS8}
    tab["R5"] = hits[0]
    return tab

CHARTAB = _build_chartab()

_orth = np.array([[np.mean(CHARTAB[a] * CHARTAB[b]) for b in IRREP_NAMES]
                  for a in IRREP_NAMES])
assert np.abs(_orth - np.eye(9)).max() < 1e-9

DIMS = {name: int(round(CHARTAB[name][0])) for name in IRREP_NAMES}
D_RHO = {"R0": 0, "R1": 1, "R3": 2, "R6": 3, "R7": 4, "R8": 5, "R4": 6, "R5": 6, "R2": 7}

def multiplicity(rho, n):
    return int(round(np.mean(np.array([chiV(n, k) for k in ELEMENT_CLASSES]) * CHARTAB[rho])))

# ---- unitarization of symmetric powers ----
def _sym_gram(n):
    """Diagonal Gram for Sym^n in monomial basis: G_kk = k!(n-k)!/n!"""
    fn = factorial(n)
    return np.array([factorial(k) * factorial(n - k) / fn for k in range(n + 1)])

def pi_n_mono(n, q):
    if n == 0:
        return np.array([[1.0 + 0j]])
    return sym_power(quat_to_su2(np.asarray(q, float)), n)

def _pi_n_cg_recurrence(n, U):
    """Sym^n in unitarized (weight) basis via CG recurrence. Stable at all n."""
    D = U.copy()
    for step in range(1, n):
        dim_old = step + 1
        dim_new = step + 2
        D_pad = np.zeros((dim_new, dim_new), dtype=complex)
        D_pad[:dim_old, :dim_old] = D
        D_pm = np.zeros_like(D_pad)
        D_pm[:, 1:] = D_pad[:, :dim_new - 1]
        D_mp = np.zeros_like(D_pad)
        D_mp[1:, :] = D_pad[:dim_new - 1, :]
        D_mm = np.zeros_like(D_pad)
        D_mm[1:, 1:] = D_pad[:dim_new - 1, :dim_new - 1]
        i_arr = np.arange(dim_new, dtype=float)
        cg_op = np.sqrt(np.clip((dim_old - i_arr) / dim_old, 0, None))
        cg_om = np.sqrt(np.clip(i_arr / dim_old, 0, None))
        CG_pp = cg_op[:, None] * cg_op[None, :]
        CG_pm = cg_op[:, None] * cg_om[None, :]
        CG_mp = cg_om[:, None] * cg_op[None, :]
        CG_mm = cg_om[:, None] * cg_om[None, :]
        D = (CG_pp * D_pad * U[0, 0] + CG_pm * D_pm * U[0, 1] +
             CG_mp * D_mp * U[1, 0] + CG_mm * D_mm * U[1, 1])
    return D

def pi_n_unitary(n, q):
    if n == 0:
        return np.array([[1.0 + 0j]])
    U = quat_to_su2(np.asarray(q, float))
    if n <= 40:
        g = _sym_gram(n)
        Dh = np.diag(np.sqrt(g))
        Dih = np.diag(1.0 / np.sqrt(g))
        return Dh @ sym_power(U, n) @ Dih
    return _pi_n_cg_recurrence(n, U)

# ---- extract R2 from V_7 via Reynolds projector ----
def _build_R2_matrices():
    """Extract R2 (dim 2) from V_7|_{2I} using the character-weighted Reynolds projector."""
    d = 8  # dim V_7
    all_V7 = np.zeros((120, d, d), dtype=complex)
    for i, q in enumerate(G120):
        all_V7[i] = pi_n_unitary(7, q)
    chi_R2 = CHARTAB["R2"]
    P = np.zeros((d, d), dtype=complex)
    for i in range(120):
        P += chi_R2[i].conj() * all_V7[i]
    P *= 2.0 / 120.0
    eigvals, eigvecs = np.linalg.eigh(P)
    mask = eigvals > 0.5
    assert mask.sum() == 2, f"expected rank 2, got {mask.sum()}"
    U = eigvecs[:, mask]  # 8x2
    R2_mats = np.zeros((120, 2, 2), dtype=complex)
    for i in range(120):
        R2_mats[i] = U.conj().T @ all_V7[i] @ U
    return R2_mats

_R2_ALL = _build_R2_matrices()

# ---- build R4 = Sym^2(R2) ----
def _sym2_of_mat(M):
    """Sym^2 of a 2x2 matrix in the monomial basis {u^2, uv, v^2}."""
    return sym_power(M, 2)

def _build_R4_matrices():
    """R4 = unitarized Sym^2(R2).
    R2 is already unitary (det 1), so Sym^2(R2) uses the same Gram as Sym^2(SU2)."""
    g = _sym_gram(2)
    Dh = np.diag(np.sqrt(g))
    Dih = np.diag(1.0 / np.sqrt(g))
    R4_mats = np.zeros((120, 3, 3), dtype=complex)
    for i in range(120):
        R4_mats[i] = Dh @ _sym2_of_mat(_R2_ALL[i]) @ Dih
    return R4_mats

_R4_ALL = _build_R4_matrices()

# ---- build R5 = R1 tensor R2 ----
def _build_R5_matrices():
    R5_mats = np.zeros((120, 4, 4), dtype=complex)
    for i, q in enumerate(G120):
        R1 = quat_to_su2(np.asarray(q, float))
        R5_mats[i] = np.kron(R1, _R2_ALL[i])
    return R5_mats

_R5_ALL = _build_R5_matrices()

# ---- precompute all rep matrices ----
def _precompute_all():
    """Precompute all 120 matrices for every sector."""
    all_reps = {}
    for rho in IRREP_NAMES:
        d = DIMS[rho]
        mats = np.zeros((120, d, d), dtype=complex)
        for i, q in enumerate(G120):
            mats[i] = _rep_rho_single(rho, i)
        all_reps[rho] = mats
    return all_reps

def _rep_rho_single(rho, idx):
    q = G120[idx]
    if rho == "R0":
        return np.array([[1.0 + 0j]])
    elif rho == "R1":
        return quat_to_su2(np.asarray(q, float))
    elif rho == "R2":
        return _R2_ALL[idx]
    elif rho == "R3":
        return pi_n_unitary(2, q)
    elif rho == "R4":
        return _R4_ALL[idx]
    elif rho == "R5":
        return _R5_ALL[idx]
    elif rho == "R6":
        return pi_n_unitary(3, q)
    elif rho == "R7":
        return pi_n_unitary(4, q)
    elif rho == "R8":
        return pi_n_unitary(5, q)
    raise ValueError(f"unknown sector: {rho}")

ALL_REPS = _precompute_all()

def rep_rho(rho, idx):
    """Unitary representation matrix of G120[idx] in sector rho."""
    return ALL_REPS[rho][idx]

def rep_rho_quat(rho, q):
    """Unitary representation matrix of quaternion q in sector rho.
    q must be in G120 (found by nearest-neighbor lookup)."""
    dists = np.sqrt(((G120 - q)**2).sum(axis=1))
    idx = int(np.argmin(dists))
    assert dists[idx] < 1e-9, f"quaternion not in 2I"
    return ALL_REPS[rho][idx]

# ---- verification ----
def _verify_reps():
    for rho in IRREP_NAMES:
        d = DIMS[rho]
        mats = ALL_REPS[rho]
        # unitarity
        for i in range(120):
            M = mats[i]
            assert M.shape == (d, d)
            assert np.abs(M @ M.conj().T - np.eye(d)).max() < 1e-9, \
                f"{rho}[{i}]: not unitary"
        # character
        for i in range(120):
            tr = np.trace(mats[i]).real
            expected = CHARTAB[rho][i]
            assert abs(tr - expected) < 1e-7, \
                f"{rho}[{i}]: trace {tr} != {expected}"
        # homomorphism on full group
        for _ in range(50):
            a, b = np.random.randint(0, 120, 2)
            p = qmul(G120[a], G120[b])
            pidx = int(np.argmin(np.sqrt(((G120 - p)**2).sum(axis=1))))
            err = np.abs(mats[a] @ mats[b] - mats[pidx]).max()
            assert err < 1e-9, f"{rho}: homomorphism error {err}"
        # rho(-1) = +-I
        minus1_idx = int(np.argmin(np.sqrt(((G120 - np.array([-1, 0, 0, 0]))**2).sum(axis=1))))
        Mm = mats[minus1_idx]
        assert np.abs(np.abs(Mm) - np.eye(d)).max() < 1e-9, f"{rho}: rho(-1) not +-I"

_verify_reps()
