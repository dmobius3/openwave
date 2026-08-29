"""Binary icosahedral group 2I: 120 unit icosians, character table, fibre realizations.

All 9 irreps of 2I. Fibre realizations are UNITARY matrix representations built
from the explicit 120-icosian set via Reynolds projection into Sym^n(SU(2)).
"""
import numpy as np
import itertools
import hashlib
from scipy.linalg import sqrtm

TAU = (1 + np.sqrt(5)) / 2

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
        coeffs = np.convolve(polyA, polyB)
        out[:, k] = coeffs
    return out

def icosians():
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

def order_of(q):
    p = q.copy()
    for k in range(1, 12):
        if np.abs(p - np.array([1, 0, 0, 0])).sum() < 1e-9:
            return k
        p = qmul(p, q)
    return -1

def closure_residual(G):
    prods = np.einsum("iab,jb->ija",
                      np.array([[[q[0], -q[1], -q[2], -q[3]],
                                 [q[1],  q[0], -q[3],  q[2]],
                                 [q[2],  q[3],  q[0], -q[1]],
                                 [q[3], -q[2],  q[1],  q[0]]] for q in G]), G)
    return np.abs(prods[:, :, None, :] - G[None, None, :, :]).sum(axis=3).min(axis=2).max()


# ---- build the group and its structure ----

G120 = icosians()

CLASS_ANGLES = np.array([0, np.pi, np.pi/2, np.pi/3, 2*np.pi/3,
                         np.pi/5, 2*np.pi/5, 3*np.pi/5, 4*np.pi/5])
GALOIS_IDX = {5: 7, 7: 5, 6: 8, 8: 6}

def snap_class(q):
    th = np.arccos(np.clip(q[0], -1, 1))
    k = int(np.argmin(np.abs(CLASS_ANGLES - th)))
    assert abs(CLASS_ANGLES[k] - th) < 1e-9, f"angle {th} off every 2I class"
    return k

KLASS = np.array([snap_class(q) for q in G120])

def chiV(n, k):
    th = CLASS_ANGLES[k]
    if k == 0:
        return float(n + 1)
    if k == 1:
        return float((n + 1) * (-1)**n)
    return np.sin((n + 1) * th) / np.sin(th)

# ---- character table: 8 from SU(2) formulas, R5 derived in-room ----

IRREP_NAMES = ["R0", "R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8"]
D_RHO = {"R0": 0, "R1": 1, "R2": 7, "R3": 2, "R4": 6, "R5": 6,
          "R6": 3, "R7": 4, "R8": 5}
DIM_RHO = {"R0": 1, "R1": 2, "R2": 2, "R3": 3, "R4": 3, "R5": 4,
            "R6": 4, "R7": 5, "R8": 6}

IRREPS8_FORMULAS = {
    "R0": lambda k: 1.0,
    "R1": lambda k: chiV(1, k),
    "R2": lambda k: chiV(1, GALOIS_IDX.get(k, k)),
    "R3": lambda k: chiV(2, k),
    "R4": lambda k: chiV(2, GALOIS_IDX.get(k, k)),
    "R6": lambda k: chiV(3, k),
    "R7": lambda k: chiV(4, k),
    "R8": lambda k: chiV(5, k),
}

def build_character_table():
    import collections
    csize = collections.Counter(KLASS)
    cent = {k: 120 // csize[k] for k in csize}
    others_vec = {r: np.array([IRREPS8_FORMULAS[r](k) for k in KLASS]) for r in IRREPS8_FORMULAS}
    absq = {k: cent[k] - sum(IRREPS8_FORMULAS[r](k)**2 for r in IRREPS8_FORMULAS) for k in csize}
    amb = [k for k in sorted(csize) if absq[k] > 0.5 and k != 0]
    hits = []
    for bits in range(2 ** len(amb)):
        row = {k: float(np.sqrt(max(absq[k], 0.0))) for k in csize}
        for j, k in enumerate(amb):
            if (bits >> j) & 1:
                row[k] = -row[k]
        vec = np.array([row[k] for k in KLASS])
        if (all(abs(np.mean(vec * others_vec[r])) < 1e-9 for r in IRREPS8_FORMULAS)
                and abs(np.mean(vec * vec) - 1.0) < 1e-9):
            hits.append(vec)
    assert len(hits) == 1, f"R5 sign solution not unique: {len(hits)} candidates"
    chartab = {r: others_vec[r] for r in IRREPS8_FORMULAS}
    chartab["R5"] = hits[0]
    dev = np.abs(np.array([[np.mean(chartab[a] * chartab[b]) for b in IRREP_NAMES]
                           for a in IRREP_NAMES]) - np.eye(9)).max()
    assert dev < 1e-9, f"character orthonormality failed: {dev}"
    return chartab

CHARTAB = build_character_table()

def multiplicity(rho, n):
    return round(float(np.mean(np.array([chiV(n, k) for k in KLASS]) * CHARTAB[rho])))


# ---- unitarized Wigner D-matrices ----

def binom_sqrt_diag(n):
    from math import comb
    return np.array([np.sqrt(comb(n, k)) for k in range(n + 1)])

def unitarized_sym_power(M, n):
    T = binom_sqrt_diag(n)
    raw = sym_power(M, n)
    return (raw * T[None, :]) / T[:, None]

def pi_unitary(n, q):
    if n == 0:
        return np.array([[1.0 + 0j]])
    return unitarized_sym_power(quat_to_su2(np.asarray(q, float)), n)


# ---- fibre realizations for all 9 irreps ----

def _find_generating_pair():
    ord10s = [i for i, q in enumerate(G120) if order_of(q) == 10]
    ord4s = [i for i, q in enumerate(G120) if order_of(q) == 4]
    for a in ord10s[:5]:
        for b in ord4s:
            S = {tuple(np.round(G120[i], 9)) for i in [a, b]}
            frontier = [G120[a], G120[b]]
            while frontier:
                nxt = []
                for x in frontier:
                    for idx in [a, b]:
                        p = qmul(x, G120[idx])
                        tp = tuple(np.round(p, 9))
                        if tp not in S:
                            S.add(tp)
                            nxt.append(p)
                frontier = nxt
            if len(S) == 120:
                return [a, b]
    raise RuntimeError("no generating pair found")

GEN_PAIR = _find_generating_pair()

def build_fibre_realization(rho_name):
    d = DIM_RHO[rho_name]
    if rho_name == "R0":
        return [np.array([[1.0 + 0j]]) for _ in G120]

    n0 = D_RHO[rho_name]
    chi_vec = CHARTAB[rho_name]
    reynolds = np.zeros((n0 + 1, n0 + 1), dtype=complex)
    for i, q in enumerate(G120):
        reynolds += chi_vec[i].conj() * pi_unitary(n0, q)
    reynolds *= d / 120.0

    eigvals, eigvecs = np.linalg.eigh(reynolds)
    mask = eigvals > 0.5
    assert mask.sum() == d, f"{rho_name}: Reynolds rank {mask.sum()} != dim {d}"
    B = eigvecs[:, mask]

    rho_mats = []
    for q in G120:
        M = B.conj().T @ pi_unitary(n0, q) @ B
        rho_mats.append(M)
    return rho_mats


def build_all_fibre_realizations():
    fibres = {}
    for rho in IRREP_NAMES:
        fibres[rho] = build_fibre_realization(rho)
    return fibres

FIBRES = build_all_fibre_realizations()


def verify_fibre_realization(rho_name, rho_mats):
    d = DIM_RHO[rho_name]
    n = len(rho_mats)
    assert n == 120

    hom_err = 0.0
    for i in range(120):
        for j in range(120):
            prod_q = qmul(G120[i], G120[j])
            k = int(np.argmin(np.array([np.linalg.norm(prod_q - G120[m]) for m in range(120)])))
            err = np.abs(rho_mats[i] @ rho_mats[j] - rho_mats[k]).max()
            hom_err = max(hom_err, err)

    unit_err = max(np.abs(M @ M.conj().T - np.eye(d)).max() for M in rho_mats)

    char_err = max(abs(np.trace(rho_mats[i]) - CHARTAB[rho_name][i]) for i in range(120))

    minus1_idx = int(np.argmin(np.array([np.linalg.norm(G120[i] - np.array([-1, 0, 0, 0]))
                                         for i in range(120)])))
    neg1_mat = rho_mats[minus1_idx]
    parity = np.trace(neg1_mat).real / d
    spinorial = D_RHO[rho_name] % 2 == 1
    expected_parity = -1.0 if spinorial else 1.0
    parity_ok = abs(parity - expected_parity) < 1e-9

    return {
        "homomorphism_err": hom_err,
        "unitarity_err": unit_err,
        "character_err": char_err,
        "parity_ok": parity_ok,
        "parity_value": parity,
        "expected_parity": expected_parity,
        "all_pass": hom_err < 1e-9 and unit_err < 1e-9 and char_err < 1e-9 and parity_ok,
    }


# ---- intertwiner spaces ----

def compute_intertwiners(rho_name, n):
    d_rho = DIM_RHO[rho_name]
    d_n = n + 1
    rho_mats = FIBRES[rho_name]

    rows = []
    for idx in GEN_PAIR:
        rho_g = rho_mats[idx]
        pi_g = pi_unitary(n, G120[idx])
        constraint = np.kron(rho_g, np.eye(d_n)) - np.kron(np.eye(d_rho), pi_g.T)
        rows.append(constraint)
    stacked = np.vstack(rows)

    U, s, Vh = np.linalg.svd(stacked)
    tol = 1e-10 * s[0] if s.size > 0 else 1e-10
    null_mask = s < tol
    extra_null = max(0, Vh.shape[0] - s.size)
    null_count = int(null_mask.sum()) + extra_null
    null_vecs = Vh[-(null_count):] if null_count > 0 else np.zeros((0, Vh.shape[1]))

    intertwiners = []
    for i in range(null_count):
        A = null_vecs[i].reshape(d_rho, d_n)
        intertwiners.append(A)

    return intertwiners, null_count


def compute_multiplicity_character(rho_name, n):
    return multiplicity(rho_name, n)
