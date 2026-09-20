"""M5.32 R23-1 / R23-2 adversarial audit (independent re-implementation).

EQUATIONS FIRST
---------------
Field: real symmetric 4 x 4 M(x), block-diagonal (M_0i = 0), spatial block S (3 x 3),
time entry M_00, u = -M_00. Lattice n^3, h = L / n, x_i = (i - (n - 1) / 2) h.
    E      = E_curv + V4 + E_lin
    E_curv = 4 h^3 sum_cells sum_{branch in fwd, bwd} (1/2) sum_{i<j} || [A_i, A_j] ||_F^2,
             A_i = one-sided difference of S along axis i (zero on the missing layer)
    V4     = w h^3 sum_cells sum_p r_p^2,   r_p = u^p + tr S^p - C_p,  C_p = sum q^p,
             q = (-g, 1, delta, delta),  w = w1s * 0.000724023879
    E_lin  = -c h^3 sum_cells sum_p a_p r_p,  sum_p p a_p x^(p-1) = (x + g)(x - 1)(x - delta)
Gradient (all nine entries of S treated as independent, the symmetric result):
    dE_curv / dA_i = 2 [F_ij, A_j] coef,  dE_curv / dA_j = -2 [F_ij, A_i] coef,
    chained through the adjoint of the one-sided difference;
    dV / dS = h^3 sum_p (2 w r_p - c a_p) p S^(p-1).
The record's fmax is max |G_ab| in that convention (an off-diagonal independent variable has
derivative 2 G_ab).

Nothing here is imported from, or copied out of, any other script in this folder. The audited
definitions were read (not imported) from m5_32_r23_1_cscan.py and the modules it loads.

Run: python3 m5_32_r23_1_audit.py      (about 45 s, 3 threads, no pools)
Writes: ../data/m5_32_r23_1_audit.json
"""

import os

for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_k] = "3"

import itertools  # noqa: E402
import json  # noqa: E402
import time  # noqa: E402
from collections import deque  # noqa: E402

import numpy as np  # noqa: E402
from scipy.ndimage import map_coordinates  # noqa: E402
from scipy.optimize import least_squares, minimize  # noqa: E402
from scipy.spatial.transform import Rotation  # noqa: E402
from scipy.special import kve  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
ROWS_JSON = os.path.join(DATA, "m5_32_r23_1_cscan.json")
EXTRA_JSON = os.path.join(DATA, "m5_32_r23_1_cscan_extra.json")
COLLECT_JSON = os.path.join(DATA, "m5_32_r23_1_cscan_collect.json")
FIELDS = os.path.join(DATA, "m5_32_r23_1")
R22S_FIELD = os.path.join(DATA, "m5_32_r22_1s", "rad_pin_d0.3_w25_n32_L48_x9000_s.npz")
OUT_JSON = os.path.join(DATA, "m5_32_r23_1_audit.json")

G8 = 8.0
W_UNIT = 0.000724023879
PIN_DEPTH = 1.6
R_CORE = 4.0
SHELLS = (1.5, 3.0, 4.5, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0)
PROFILE_R = np.arange(3.0, 22.6, 0.75)
NU = np.sqrt(5.0) / 4.0
T0 = time.time()
CHECKS = []
FROM_EXTRA = []


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def add(cid, claim, method, value, expectation, fails_if, ok):
    CHECKS.append(
        {
            "id": cid,
            "claim": claim,
            "method": method,
            "value": value,
            "expectation": expectation,
            "fails_if": fails_if,
            "verdict": "PASS" if ok else "FAIL",
        }
    )
    log(f"{cid}: {'PASS' if ok else 'FAIL'}")


# ======================= the model, re-implemented =======================
class Par:
    def __init__(self, delta, w1s, c, n, L):
        self.delta, self.c, self.n, self.L = float(delta), float(c), int(n), float(L)
        self.h = self.L / self.n
        self.w = float(w1s) * W_UNIT
        q = np.array([-G8, 1.0, self.delta, self.delta])
        self.Cp = [float(np.sum(q**p)) for p in range(1, 5)]
        poly = np.polymul(np.polymul([1.0, G8], [1.0, -1.0]), [1.0, -self.delta])
        asc = poly[::-1]
        self.a = [float(asc[p - 1]) / p for p in range(1, 5)]
        x = (np.arange(self.n) - (self.n - 1) / 2.0) * self.h
        self.X = np.stack(np.meshgrid(x, x, x, indexing="ij"), axis=-1)
        self.r = np.sqrt(np.sum(self.X**2, axis=-1))
        k = max(1, int(np.ceil(PIN_DEPTH / self.h)))
        free = np.zeros((self.n,) * 3, dtype=bool)
        free[k:-k, k:-k, k:-k] = True
        self.free = free
        self.gate = 1e-3 * ((1.0 - self.delta) / 0.7) ** 4


def one_sided(S, ax, h, fwd):
    out = np.zeros_like(S)
    d = np.diff(S, axis=ax) / h
    idx = [slice(None)] * S.ndim
    idx[ax] = slice(0, -1) if fwd else slice(1, None)
    out[tuple(idx)] = d
    return out


def curv_density(S, h, ax0=0):
    """per-cell curvature energy (h^3 weighted); the spatial axes start at ax0."""
    e = np.zeros(S.shape[:-2])
    for fwd in (True, False):
        A = [one_sided(S, ax0 + i, h, fwd) for i in range(3)]
        for i in range(3):
            for j in range(i + 1, 3):
                P = A[i] @ A[j]
                F = P - np.swapaxes(P, -1, -2)
                e += 0.5 * 4.0 * np.sum(F * F, axis=(-1, -2))
    return h**3 * e


def curv_energy_4x4(M, h):
    """the literal 4 x 4 eta form, used only to validate the 3 x 3 reduction."""
    eta = np.diag([-1.0, 1.0, 1.0, 1.0])
    tot = 0.0
    for fwd in (True, False):
        A = [one_sided(M, i, h, fwd) for i in range(3)]
        for i in range(3):
            for j in range(i + 1, 3):
                F = A[i] @ eta @ A[j] - A[j] @ eta @ A[i]
                tot += 0.5 * np.sum(np.einsum("...ab,ac,bd,...cd->...", F, eta, eta, F))
    return 4.0 * h**3 * float(tot)


def traces(S):
    S2 = S @ S
    s1 = np.einsum("...aa->...", S)
    s2 = np.einsum("...aa->...", S2)
    s3 = np.sum(S2 * S, axis=(-1, -2))
    s4 = np.sum(S2 * S2, axis=(-1, -2))
    return [s1, s2, s3, s4]


def residuals(S, m00, P):
    u = -m00
    s = traces(S)
    return [u ** (p + 1) + s[p] - P.Cp[p] for p in range(4)], u


def pot_density(S, m00, P):
    """(V4, E_lin) per cell, h^3 weighted."""
    r, _ = residuals(S, m00, P)
    v4 = P.w * sum(q * q for q in r)
    lin = -P.c * sum(P.a[p] * r[p] for p in range(4))
    return P.h**3 * v4, P.h**3 * lin


def solve_u(S, P, m_start, iters=40):
    """my own per-cell minimizer over M_00 (Newton with a convexity guard)."""
    s = traces(S)
    u = -np.array(m_start, dtype=float)
    for _ in range(iters):
        r = [u ** (p + 1) + s[p] - P.Cp[p] for p in range(4)]
        k = [2.0 * P.w * r[p] - P.c * P.a[p] for p in range(4)]
        d = [(p + 1) * u**p for p in range(4)]
        dd = [(p + 1) * p * u ** (p - 1) if p > 0 else 0.0 * u for p in range(4)]
        g = sum(k[p] * d[p] for p in range(4))
        H = sum(2.0 * P.w * d[p] ** 2 + k[p] * dd[p] for p in range(4))
        u = u - g / np.where(H > 0, H, np.inf)
    return -u


def m00_stationarity(S, m00, P):
    """g = d(potential density)/du, its Hessian, the sum of |terms|, and the roundoff floor of g
    (the p = 4 residual is a difference of numbers of size u^4, so g cannot be resolved below
    2 w 4 |u|^3 times a few ulp of u^4 + tr S^4 + C_4)."""
    r, u = residuals(S, m00, P)
    k = [2.0 * P.w * r[p] - P.c * P.a[p] for p in range(4)]
    d = [(p + 1) * u**p for p in range(4)]
    dd = [(p + 1) * p * u ** (p - 1) if p > 0 else 0.0 * u for p in range(4)]
    g = sum(k[p] * d[p] for p in range(4))
    H = sum(2.0 * P.w * d[p] ** 2 + k[p] * dd[p] for p in range(4))
    scale = sum((np.abs(2.0 * P.w * r[p]) + abs(P.c * P.a[p])) * np.abs(d[p]) for p in range(4))
    s4 = traces(S)[3]
    floor = 2.0 * P.w * np.abs(d[3]) * 8.0 * np.finfo(float).eps * (u**4 + s4 + P.Cp[3])
    return g, H, scale, u, floor


def energy_parts(S, m00, P):
    ec = curv_density(S, P.h)
    v4, lin = pot_density(S, m00, P)
    return ec, v4, lin


def gradient_S(S, m00, P):
    """analytic dE/dS at fixed M_00, nine-entry convention (symmetric)."""
    h = P.h
    Gs = np.zeros_like(S)
    coef = 0.5 * 4.0 * h**3
    for fwd in (True, False):
        A = [one_sided(S, i, h, fwd) for i in range(3)]
        GA = [np.zeros_like(S) for _ in range(3)]
        for i in range(3):
            for j in range(i + 1, 3):
                Pm = A[i] @ A[j]
                F = Pm - np.swapaxes(Pm, -1, -2)
                GA[i] += coef * 2.0 * (F @ A[j] - A[j] @ F)
                GA[j] -= coef * 2.0 * (F @ A[i] - A[i] @ F)
        for ax in range(3):
            lo = [slice(None)] * S.ndim
            hi = [slice(None)] * S.ndim
            lo[ax], hi[ax] = slice(0, -1), slice(1, None)
            src = tuple(lo) if fwd else tuple(hi)
            Gs[tuple(hi)] += GA[ax][src] / h
            Gs[tuple(lo)] -= GA[ax][src] / h
    r, _ = residuals(S, m00, P)
    k = [2.0 * P.w * r[p] - P.c * P.a[p] for p in range(4)]
    I3 = np.broadcast_to(np.eye(3), S.shape)
    S2 = S @ S
    S3 = S2 @ S
    pw = [I3, S, S2, S3]
    for p in range(4):
        Gs += P.h**3 * ((p + 1) * k[p])[..., None, None] * pw[p]
    return Gs


def fd_gradient(S, m00, P, cells, ents, t=5e-4, resolve=True):
    """five-point central differences of my own energy on 5^3 patches (the curvature part is a
    quartic polynomial in S, so this stencil is exact for it); M_00 re-solved when resolve."""
    B = len(cells)
    pat = np.empty((B, 5, 5, 5, 3, 3))
    for b, (i, j, k) in enumerate(cells):
        pat[b] = S[i - 2 : i + 3, j - 2 : j + 3, k - 2 : k + 3]
    m0 = np.array([m00[i, j, k] for i, j, k in cells])
    out = {}
    for mult in (1.0, -1.0, 2.0, -2.0):
        q = pat.copy()
        for b, (a, c) in enumerate(ents):
            q[b, 2, 2, 2, a, c] += mult * t
            if a != c:
                q[b, 2, 2, 2, c, a] += mult * t
        ec = curv_density(q, P.h, ax0=1)[:, 1:4, 1:4, 1:4].sum(axis=(1, 2, 3))
        Sc = q[:, 2, 2, 2]
        mm = solve_u(Sc, P, m0) if resolve else m0
        v4, lin = pot_density(Sc, mm, P)
        out[mult] = ec + v4 + lin
    return (8.0 * (out[1.0] - out[-1.0]) - (out[2.0] - out[-2.0])) / (12.0 * t)


# ======================= the reads, re-implemented =======================
def r_half_of(e, r):
    """the smallest cell radius R with E(r <= R) >= E / 2, plus a linearly interpolated version."""
    rr, inv = np.unique(np.round(r.ravel(), 9), return_inverse=True)
    es = np.bincount(inv, weights=e.ravel())
    cum = np.cumsum(es)
    half = 0.5 * cum[-1]
    i = int(np.searchsorted(cum, half))
    lo_c = cum[i - 1] if i > 0 else 0.0
    lo_r = rr[i - 1] if i > 0 else 0.0
    interp = lo_r + (half - lo_c) * (rr[i] - lo_r) / (cum[i] - lo_c)
    step = rr[i] - lo_r
    return float(rr[i]), float(interp), float(step), bool(np.all(np.diff(cum) >= 0))


def shell_means(lam, r, R, half_width):
    sh = np.abs(r - R) < half_width
    if not sh.any():
        return None
    v = lam[sh]
    return {
        "mean": [float(v[:, a].mean()) for a in range(3)],
        "std": [float(v[:, a].std()) for a in range(3)],
        "cells": int(sh.sum()),
    }


def eps_profile(lam, r, h, half_width):
    eps = 0.5 * (lam[..., 1] - lam[..., 0])
    out = []
    for R in PROFILE_R:
        sh = np.abs(r - R) < half_width
        if sh.any():
            out.append([float(R), float(eps[sh].mean())])
    return np.array(out)


def make_seed(P, kind):
    rh = P.X / np.maximum(P.r, 1e-300)[..., None]
    I3 = np.eye(3)
    ext = P.delta * I3 + (1.0 - P.delta) * rh[..., :, None] * rh[..., None, :]
    if kind == "rad":
        C = (1.0 + 2.0 * P.delta) / 3.0 * I3
    else:
        e = min(P.delta, 0.3)
        C = np.diag([1.0, P.delta + e, P.delta - e])
    u = (1.0 - np.exp(-((P.r / R_CORE) ** 2)))[..., None, None]
    return u * ext + (1.0 - u) * C, ext


# ======================= the degree reader, re-implemented =======================
def cube_surface(lo, hi):
    """faces of the cube [lo, hi]^3 as (m, m, 3) index grids."""
    rng = np.arange(lo, hi + 1)
    faces = []
    for ax in range(3):
        for side in (lo, hi):
            b, c = np.meshgrid(rng, rng, indexing="ij")
            idx = np.zeros(b.shape + (3,), dtype=int)
            others = [a for a in range(3) if a != ax]
            idx[..., ax] = side
            idx[..., others[0]] = b
            idx[..., others[1]] = c
            faces.append(idx)
    return faces


def solid_angle(a, b, c):
    num = np.einsum("...i,...i->...", a, np.cross(b, c))
    den = (
        1.0
        + np.einsum("...i,...i->...", a, b)
        + np.einsum("...i,...i->...", b, c)
        + np.einsum("...i,...i->...", c, a)
    )
    return 2.0 * np.arctan2(num, den)


def degree_read(vec_top, gap_tm, P, lo, hi):
    n = P.n
    faces = cube_surface(lo, hi)
    center = 0.5 * (lo + hi)
    lin = lambda ix: (ix[..., 0] * n + ix[..., 1]) * n + ix[..., 2]  # noqa: E731
    v = vec_top.reshape(-1, 3)
    gap = gap_tm.ravel()
    links = []
    frustrated, quad_gap = 0, []
    quads = []
    for f in faces:
        L = lin(f)
        links.append(np.stack([L[:-1, :].ravel(), L[1:, :].ravel()], axis=1))
        links.append(np.stack([L[:, :-1].ravel(), L[:, 1:].ravel()], axis=1))
        q = [L[:-1, :-1], L[1:, :-1], L[1:, 1:], L[:-1, 1:]]
        prod = np.ones(q[0].shape)
        for s in range(4):
            prod *= np.sign(np.sum(v[q[s]] * v[q[(s + 1) % 4]], axis=-1))
        bad = prod < 0
        frustrated += int(bad.sum())
        gq = np.minimum.reduce([gap[q[s]] for s in range(4)])
        quad_gap += list(gq[bad])
        quads.append(np.stack([x.ravel() for x in q], axis=1))
    links = np.unique(np.sort(np.concatenate(links), axis=1), axis=0)
    quads = np.concatenate(quads)
    verts = np.unique(links)
    # orientation 1: breadth-first flood fill over the surface graph
    adj = {int(a): [] for a in verts}
    for a, b in links:
        adj[int(a)].append(int(b))
        adj[int(b)].append(int(a))
    sign = {int(verts[0]): 1.0}
    dq = deque([int(verts[0])])
    while dq:
        a = dq.popleft()
        for b in adj[a]:
            if b not in sign:
                sign[b] = sign[a] * (1.0 if float(v[a] @ v[b]) >= 0 else -1.0)
                dq.append(b)
    sg = np.zeros(len(v))
    for a, s_ in sign.items():
        sg[a] = s_
    vo = v * sg[:, None]
    conf_bfs = int(np.sum(np.sum(vo[links[:, 0]] * vo[links[:, 1]], axis=1) < 0))
    # orientation 2: along r-hat from the box center
    pos = P.X.reshape(-1, 3)
    sr = np.where(np.sum(v * pos, axis=1) < 0, -1.0, 1.0)
    vr = v * sr[:, None]
    conf_rhat = int(np.sum(np.sum(vr[links[:, 0]] * vr[links[:, 1]], axis=1) < 0))
    # the degree, both quad splittings averaged, triangles ordered outward geometrically
    ijk = np.stack(np.unravel_index(np.arange(len(v)), (n, n, n)), axis=1).astype(float)
    ijk -= center

    def deg_of(field):
        tot = 0.0
        for tri in ((0, 1, 2), (0, 2, 3), (0, 1, 3), (1, 2, 3)):
            p, q, r_ = quads[:, tri[0]], quads[:, tri[1]], quads[:, tri[2]]
            orient = np.sign(np.einsum("ni,ni->n", ijk[p], np.cross(ijk[q], ijk[r_])))
            om = solid_angle(field[p], field[q], field[r_]) * orient
            tot += 0.5 * float(om.sum())
        return tot / (4.0 * np.pi)

    if np.mean(np.sum(vo[verts] * pos[verts], axis=1)) < 0:
        vo = -vo
    return {
        "frustrated_plaquettes": frustrated,
        "conflicts_bfs": conf_bfs,
        "conflicts_rhat": conf_rhat,
        "degree_bfs": deg_of(vo) if frustrated == 0 else None,
        "degree_rhat": deg_of(vr) if conf_rhat == 0 else None,
        "gap_top_mid_at_frustrated": [float(x) for x in sorted(quad_gap)[:6]],
        "gap_top_mid_surface_median": float(np.median(gap[verts])),
        "gap_top_mid_surface_min": float(np.min(gap[verts])),
    }


def volume_disclinations(vec_top, P):
    """frustrated elementary squares of the top-eigenvector line field in the whole volume."""
    v = vec_top
    n = P.n
    rc = []
    total = 0
    for a, b in ((0, 1), (0, 2), (1, 2)):

        def sh(da, db, a=a, b=b):
            idx = [slice(None)] * 3
            idx[a] = slice(da, n - 1 + da)
            idx[b] = slice(db, n - 1 + db)
            return tuple(idx)

        q = [v[sh(0, 0)], v[sh(1, 0)], v[sh(1, 1)], v[sh(0, 1)]]
        prod = np.ones(q[0].shape[:-1])
        for s in range(4):
            prod *= np.sign(np.sum(q[s] * q[(s + 1) % 4], axis=-1))
        rq = 0.25 * (P.r[sh(0, 0)] + P.r[sh(1, 0)] + P.r[sh(1, 1)] + P.r[sh(0, 1)])
        fq = P.free[sh(0, 0)] & P.free[sh(1, 0)] & P.free[sh(1, 1)] & P.free[sh(0, 1)]
        bad = (prod < 0) & fq
        total += int(bad.sum())
        rc += list(rq[bad])
    rc = np.array(rc)
    return {
        "frustrated_squares_free_region": total,
        "r_min": float(rc.min()) if total else None,
        "r_median": float(np.median(rc)) if total else None,
        "r_max": float(rc.max()) if total else None,
    }


# ======================= symmetry images =======================
def images():
    out = []
    for perm in itertools.permutations(range(3)):
        for sg in itertools.product((1.0, -1.0), repeat=3):
            out.append((perm, sg))
    return out


def apply_image(S, m00, perm, sg):
    T = np.transpose(S, perm + (3, 4))
    m = np.transpose(m00, perm)
    for a in range(3):
        if sg[a] < 0:
            T = np.flip(T, axis=a)
            m = np.flip(m, axis=a)
    T = T[..., list(perm), :][..., :, list(perm)]
    s = np.array(sg)
    T = T * s[:, None] * s[None, :]
    return np.ascontiguousarray(T), np.ascontiguousarray(m)


def is_lattice_exact(perm, sg):
    return len(set(sg)) == 1


def image_matrix(perm, sg):
    Rm = np.zeros((3, 3))
    for a in range(3):
        Rm[a, perm[a]] = sg[a]
    return Rm


def rotated_copy(S, Rm, P):
    """S'(x) = R S(R^T x) R^T by trilinear interpolation (exact for lattice rotations)."""
    n = P.n
    idx = (P.X.reshape(-1, 3) @ Rm) / P.h + (n - 1) / 2.0
    out = np.empty((n**3, 3, 3))
    for a in range(3):
        for b in range(a, 3):
            v = map_coordinates(S[..., a, b], idx.T, order=1, mode="nearest")
            out[:, a, b] = v
            out[:, b, a] = v
    out = Rm @ out @ Rm.T
    return out.reshape(n, n, n, 3, 3)


def best_rotation(SA, SB, P, seed, n_random=150, n_refine=4):
    """min over SO(3) (and SO(3) times inversion) of the RMS of SA - rotated SB on r < 18."""
    mask = P.r < 18.0
    rng_r = np.random.default_rng(seed)

    def rms_of(Rm, SBx):
        D_ = SA[mask] - rotated_copy(SBx, Rm, P)[mask]
        return float(np.sqrt(np.mean(np.sum(D_**2, axis=(-1, -2)))))

    out = {}
    for name, SBx in (
        ("proper", SB),
        ("with_inversion", np.ascontiguousarray(SB[::-1, ::-1, ::-1])),
    ):
        cands = [np.eye(3)]
        cands += [image_matrix(p_, s_) for p_, s_ in images()]
        cands = [c_ for c_ in cands if np.linalg.det(c_) > 0]
        cands += list(Rotation.random(n_random, random_state=rng_r).as_matrix())
        vals = np.array([rms_of(c_, SBx) for c_ in cands])
        order = np.argsort(vals)[:n_refine]
        best = (float(vals[order[0]]), cands[order[0]])
        for i in order:
            rv0 = Rotation.from_matrix(cands[i]).as_rotvec()
            res = minimize(
                lambda rv: rms_of(Rotation.from_rotvec(rv).as_matrix(), SBx),
                rv0,
                method="Nelder-Mead",
                options={"maxiter": 60, "xatol": 1e-3, "fatol": 1e-5},
            )
            if res.fun < best[0]:
                best = (float(res.fun), Rotation.from_rotvec(res.x).as_matrix())
        out[name] = {
            "rms_identity": float(vals[0]),
            "rms_random_median": float(np.median(vals[-n_random:])),
            "rms_min": best[0],
            "angle_deg_of_best": float(
                np.degrees(np.linalg.norm(Rotation.from_matrix(best[1]).as_rotvec()))
            ),
        }
    return out


# ======================= fits =======================
def fit_candidates(R, y):
    """three two-parameter candidates, least squares in log space and in linear space."""
    out = {}

    def knu(par, r):
        z = 0.5 * np.exp(par[1]) * r**2
        # sqrt(r) K_nu(z) = sqrt(r) kve(nu, z) exp(-z)
        return par[0] + 0.5 * np.log(r) + np.log(kve(NU, z)) - z

    def yuk(par, r):
        return par[0] - r / np.exp(par[1]) - np.log(r)

    def powl(par, r):
        return par[0] - par[1] * np.log(r)

    starts = {
        "knu_cutoff": (knu, [[0.0, np.log(b)] for b in (1e-3, 1e-2, 3e-2, 1e-1)]),
        "yukawa": (yuk, [[0.0, np.log(s)] for s in (1.0, 3.0, 10.0)]),
        "power": (powl, [[0.0, p] for p in (0.6, 2.0, 5.0)]),
    }
    ly = np.log(y)
    for name, (fn, st) in starts.items():
        best = None
        for p0 in st:
            p0 = [ly[0] - fn([0.0, p0[1]], R)[0], p0[1]]
            res = least_squares(lambda p, fn=fn: fn(p, R) - ly, p0)
            rl = least_squares(lambda p, fn=fn: np.exp(fn(p, R)) - y, res.x)
            rec = {
                "rms_log": float(np.sqrt(np.mean(res.fun**2))),
                "rms_lin_rel_to_max": float(np.sqrt(np.mean(rl.fun**2)) / y.max()),
                "par_log_fit": [float(q) for q in res.x],
                "par_lin_fit": [float(q) for q in rl.x],
            }
            if best is None or rec["rms_log"] < best["rms_log"]:
                best = rec
        if name == "knu_cutoff":
            best["beta_log_fit"] = float(np.exp(best["par_log_fit"][1]) ** 2)
            best["R_cut_log_fit"] = float(np.exp(best["par_log_fit"][1]) ** -0.5)
            best["R_cut_lin_fit"] = float(np.exp(best["par_lin_fit"][1]) ** -0.5)
        elif name == "yukawa":
            best["R_log_fit"] = float(np.exp(best["par_log_fit"][1]))
            best["R_lin_fit"] = float(np.exp(best["par_lin_fit"][1]))
        else:
            best["exponent_log_fit"] = best["par_log_fit"][1]
            best["exponent_lin_fit"] = best["par_lin_fit"][1]
        out[name] = best
    return out


# ======================= main =======================
def main():
    with open(ROWS_JSON) as f:
        rows = json.load(f)["rows"]
    extra_rows = {}
    if os.path.exists(EXTRA_JSON):
        with open(EXTRA_JSON) as f:
            extra_rows = json.load(f)["rows"]
    tags = [
        t
        for t, r in rows.items()
        if r.get("status") == "OK"
        and r.get("n") == 32
        and os.path.exists(os.path.join(FIELDS, t + ".npz"))
    ]
    tags.sort()
    log(f"{len(tags)} finished n32 rows")
    D = {}
    for t in tags:
        row = rows[t]
        P = Par(row["delta"], row["w1s"], row["c"], row["n"], row["L"])
        M = np.load(os.path.join(FIELDS, t + ".npz"))["M"]
        last = row["chunks"][-1]
        S = np.ascontiguousarray(M[..., 1:, 1:])
        m00 = np.ascontiguousarray(M[..., 0, 0])
        ec, v4, lin = energy_parts(S, m00, P)
        if abs(ec.sum() - last["E_curv"]) > 1e-6 * abs(last["E_curv"]) and t in extra_rows:
            row = extra_rows[t]
            last = row["chunks"][-1]
            FROM_EXTRA.append(t)
        lam, vec = np.linalg.eigh(S)
        D[t] = {
            "row": row, "P": P, "M": M, "S": S, "m00": m00, "last": last,
            "ec": ec, "v4": v4, "lin": lin, "e": ec + v4 + lin, "lam": lam, "vec": vec,
        }  # fmt: skip

    def tag_of(seed, delta, w1s, c, src=""):
        return f"{seed}_pin_d{delta:.4g}_w{w1s:.4g}_c{c:g}_n32_L48{src}"

    c_scan = (0.0, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3)
    T_RAD = {c: tag_of("rad", 0.3, 25.0, c) for c in c_scan}
    T_TWO = tag_of("rad", 0.3, 25.0, 3e-4, "_from_r22s")
    T_BIA = {c: tag_of("bia", 0.3, 25.0, c) for c in (1e-4, 1e-3)}
    d_phys = 1.0 - ((1.0 / 137.035999) / (64.0 * np.pi)) ** 0.25
    T_A0 = tag_of("rad", d_phys, 25.0, 0.0)
    T_AC = tag_of("rad", d_phys, 25.0, 3.094e-05)
    T_B0 = tag_of("rad", d_phys, 25.0 / (0.7 / (1.0 - d_phys)) ** 4, 0.0)
    expected = list(T_RAD.values()) + [T_TWO] + list(T_BIA.values()) + [T_A0, T_AC, T_B0]
    missing = [t for t in expected if t not in D]
    add(
        "C0.rows",
        "all 12 finished n32 rows named in the claims exist as end fields",
        "tag reconstruction from the claim text against the JSON and the npz folder",
        {
            "found": len(tags),
            "missing": missing,
            "rows_whose_end_field_matches_the_continuation_pool_json": FROM_EXTRA,
        },
        "no missing tag",
        "a claimed row has no status OK row or no end field",
        not missing,
    )

    # ---------------- C1: the instrument ----------------
    blk, sym = {}, {}
    for t, d in D.items():
        M = d["M"]
        blk[t] = float(max(np.abs(M[..., 0, 1:]).max(), np.abs(M[..., 1:, 0]).max()))
        sym[t] = float(np.abs(d["S"] - np.swapaxes(d["S"], -1, -2)).max())
    add(
        "C1.block_diagonal",
        "every end field is block-diagonal (M_0i = 0) with a symmetric spatial block",
        "max |M_0i| and max |S - S^T| over all cells",
        {"max_M0i": max(blk.values()), "max_asym": max(sym.values())},
        "both exactly zero",
        "any nonzero M_0i or asymmetry above 1e-15",
        max(blk.values()) == 0.0 and max(sym.values()) < 1e-15,
    )
    two = [T_RAD[3e-3], T_B0]
    v44 = {}
    for t in two:
        d = D[t]
        e4 = curv_energy_4x4(d["M"], d["P"].h)
        v44[t] = {
            "mine_4x4_eta": e4,
            "mine_3x3": float(d["ec"].sum()),
            "stored": d["last"]["E_curv"],
            "rel_4x4_vs_stored": abs(e4 - d["last"]["E_curv"]) / abs(d["last"]["E_curv"]),
            "rel_3x3_vs_4x4": abs(e4 - d["ec"].sum()) / abs(e4),
        }
    add(
        "C1.curv_two_rows",
        "my curvature energy (literal 4 x 4 eta form and the 3 x 3 reduction) reproduces the "
        "stored E_curv of two rows before it is used for anything else",
        "own one-sided stencils (zero on the missing layer), fwd and bwd each weighted 1/2",
        v44,
        "relative difference under 1e-8 for both rows and both forms",
        "any relative difference at or above 1e-8",
        all(v["rel_4x4_vs_stored"] < 1e-8 and v["rel_3x3_vs_4x4"] < 1e-8 for v in v44.values()),
    )
    rep, worst = {}, 0.0
    for t, d in D.items():
        mine = {
            "E": float(d["e"].sum()),
            "E_curv": float(d["ec"].sum()),
            "V4": float(d["v4"].sum()),
            "E_lin": float(d["lin"].sum()),
        }
        rel = {}
        for k, v in mine.items():
            ref = d["last"][k]
            rel[k] = abs(v - ref) / abs(ref) if ref != 0 else abs(v)
        rep[t] = {"mine": mine, "rel_diff": rel}
        worst = max(worst, max(rel.values()))
    add(
        "C1.energies_all_rows",
        "stored E, E_curv, V4, E_lin of every finished n32 row reproduce from its end field",
        "own energy code on the end field against the last chunk of the row",
        {"worst_rel_diff": worst, "rows": rep},
        "every relative difference under 1e-8 (absolute for a stored zero)",
        "any part of any row off by 1e-8 relative or more",
        worst < 1e-8,
    )
    pins = {}
    for t, d in D.items():
        seed, _ = make_seed(d["P"], d["row"]["seed"])
        pin = ~d["P"].free
        pins[t] = {
            "max_abs_S_minus_seed_on_pin": float(np.abs(d["S"][pin] - seed[pin]).max()),
            "max_abs_M00_minus_g_on_pin": float(np.abs(d["m00"][pin] - G8).max()),
        }
    wp = max(max(v.values()) for v in pins.values())
    add(
        "C1.pinned_shell",
        "the outer shell of depth 1.6 (2 cells at h 1.5) is pinned at the seed in every end field",
        "own seed construction (rad and bia cores, r_c 4) compared on the shell cells",
        {"worst": wp, "rows": pins},
        "difference under 1e-13 (the from_r22s row is pinned at the same analytic hedgehog)",
        "any pinned cell differs from my seed by 1e-13 or more",
        wp < 1e-13,
    )
    m00r = {}
    for t, d in D.items():
        P = d["P"]
        g, H, scale, u, floor = m00_stationarity(d["S"], d["m00"], P)
        fr = P.free
        ratio = np.abs(g[fr]) / scale[fr]
        resolved = np.abs(g[fr]) > floor[fr]
        step = np.abs(g[fr] / H[fr]) / np.abs(u[fr])
        mine = solve_u(d["S"][fr], P, d["m00"][fr])
        m00r[t] = {
            "max_abs_g_over_sum_abs_terms_raw": float(ratio.max()),
            "max_ratio_on_cells_where_g_is_above_its_roundoff_floor": float(
                ratio[resolved].max() if resolved.any() else 0.0
            ),
            "fraction_of_free_cells_with_g_above_roundoff_floor": float(resolved.mean()),
            "max_abs_g_over_roundoff_floor": float((np.abs(g[fr]) / floor[fr]).max()),
            "max_newton_step_rel": float(step.max()),
            "min_hessian_over_w": float(H[fr].min() / P.w),
            "max_abs_M00_minus_my_minimizer": float(np.abs(mine - d["m00"][fr]).max()),
        }
    wr = max(v["max_ratio_on_cells_where_g_is_above_its_roundoff_floor"] for v in m00r.values())
    ws = max(v["max_newton_step_rel"] for v in m00r.values())
    wm = max(v["max_abs_M00_minus_my_minimizer"] for v in m00r.values())
    hpos = min(v["min_hessian_over_w"] for v in m00r.values())
    add(
        "C1.m00_minimizer",
        "M_00 in the end fields is the per-cell minimizer of the potential density (its M_00 "
        "derivative is zero to 1e-8 of its scale)",
        "g = d/dM_00 of w sum r_p^2 - c sum a_p r_p on free cells against the sum of the "
        "absolute values of its terms, on the cells where g is above its double-precision "
        "floor (2 w 4 |u|^3 times 8 ulp of u^4 + tr S^4 + C_4); the relative Newton step g / "
        "(H u); the Hessian sign; the difference from my own Newton solve. The raw ratio is "
        "reported too: it reaches 5e-4 only in near-vacuum cells where the sum of |terms| is "
        "itself of roundoff size.",
        {
            "worst_resolved_ratio": wr,
            "worst_newton_step_rel": ws,
            "worst_diff_from_my_minimizer": wm,
            "min_hessian_over_w": hpos,
            "rows": m00r,
        },
        "resolved ratio under 1e-8, Newton step under 1e-12 of |M_00|, my minimizer within "
        "1e-12, a positive Hessian",
        "a free cell with a resolvable M_00 derivative at or above 1e-8 of its scale, a Newton "
        "step at or above 1e-12, or a non-positive Hessian",
        wr < 1e-8 and ws < 1e-12 and wm < 1e-12 and hpos > 0,
    )
    rng = np.random.default_rng(20260919)
    gr = {}
    ok_gate, ok_fd, ok_label = True, True, True
    for t, d in D.items():
        P = d["P"]
        Gs = gradient_S(d["S"], d["m00"], P)
        Gf = np.abs(Gs) * P.free[..., None, None]
        fmax = float(Gf.max())
        free_idx = np.argwhere(P.free)
        pick = free_idx[rng.choice(len(free_idx), size=260, replace=False)]
        ents = [tuple(sorted(int(x) for x in rng.choice(3, size=2))) for _ in range(260)]
        top = np.argsort(Gf.ravel())[-20:]
        for q in top:
            i, j, k, a, b = np.unravel_index(q, Gf.shape)
            pick = np.vstack([pick, [i, j, k]])
            ents.append((min(a, b), max(a, b)))
        cells = [tuple(int(x) for x in p) for p in pick]
        ana = np.array([Gs[c][a, b] * (1.0 if a == b else 2.0) for c, (a, b) in zip(cells, ents)])
        fd_res = fd_gradient(d["S"], d["m00"], P, cells, ents, resolve=True)
        fd_fix = fd_gradient(d["S"], d["m00"], P, cells, ents, resolve=False)
        ref = max(np.abs(ana).max(), 1e-3 * P.gate)
        off = Gf.copy()
        off[..., 0, 0] = off[..., 1, 1] = off[..., 2, 2] = 0.0
        gr[t] = {
            "label": d["row"]["label"],
            "gate": P.gate,
            "fmax_mine_record_convention": fmax,
            "fmax_stored": d["last"].get("fmax_spatial"),
            "fmax_independent_variable_convention": float(
                max(Gf[..., [0, 1, 2], [0, 1, 2]].max(), 2.0 * off.max())
            ),
            "rms_G_free": float(np.sqrt(np.mean(Gs[P.free] ** 2))),
            "fd_samples": len(cells),
            "fd_resolved_vs_analytic_max_abs_over_max_G": float(np.abs(fd_res - ana).max() / ref),
            "fd_fixed_m00_vs_resolved_max_abs_over_max_G": float(
                np.abs(fd_fix - fd_res).max() / ref
            ),
            "fd_max_abs_record_convention": float(
                max(abs(v) / (1.0 if a == b else 2.0) for v, (a, b) in zip(fd_res, ents))
            ),
        }
        if d["row"]["label"] == "AT_GATE":
            ok_gate &= gr[t]["fd_max_abs_record_convention"] < P.gate and fmax < P.gate
        ok_fd &= gr[t]["fd_resolved_vs_analytic_max_abs_over_max_G"] < 1e-4
        ok_label &= abs(fmax - d["last"]["fmax_spatial"]) < 1e-6 * fmax + 1e-16
        d["G"] = Gs
    add(
        "C1.gradient_fd",
        "my analytic gradient equals central differences of my own reduced energy (M_00 "
        "re-solved per perturbed cell), so the envelope argument holds",
        "280 (cell, entry) samples per row: 260 random free cells plus the 20 largest |G| "
        "entries; 5^3 patches, five-point stencil, step 5e-4; M_00 re-solved by my own Newton "
        "(the envelope argument is tested by also holding M_00: the two must agree)",
        {t: {k: v for k, v in g.items() if k.startswith("fd_")} for t, g in gr.items()},
        "max |FD - analytic| under 1e-4 of max(max |G|, gate / 1000) on every row",
        "a row where the finite differences disagree with the analytic gradient by 1e-4 of max "
        "|G| or more",
        ok_fd,
    )
    add(
        "C1.at_gate_stationary",
        "rows labeled AT_GATE have max |dE/dS| over free cells below 1e-3 "
        "((1 - delta) / 0.7)^4 in the record's convention",
        "own analytic gradient over ALL free cells, its largest entries confirmed by finite "
        "differences; the independent-variable convention (2 G_ab off the diagonal) is reported",
        {t: g for t, g in gr.items()},
        "every AT_GATE row under its gate, and my fmax equal to the stored fmax",
        "an AT_GATE row whose max |G| (mine, full field, or FD) is at or above the gate, or a "
        "stored fmax I cannot reproduce to 1e-6",
        ok_gate and ok_label,
    )
    conv = {
        t: g["fmax_independent_variable_convention"] / g["gate"]
        for t, g in gr.items()
        if g["label"] == "AT_GATE"
    }
    add(
        "C1.gate_convention",
        "WORDING: 'max |dE/dS| below the gate' holds when dE/dS means the derivative with "
        "respect to the independent symmetric entries",
        "max over free cells of |dE/dS_aa| and 2 |G_ab| against the gate",
        conv,
        "ratio under 1 for every AT_GATE row",
        "an AT_GATE row passes only in the halved (sym4) convention",
        all(v < 1.0 for v in conv.values()),
    )

    # ---------------- C2: the half-energy radius ----------------
    stored_rh = {0.0: 17.02, 3e-5: 16.48, 1e-4: 15.93, 3e-4: 13.97, 1e-3: 11.30, 3e-3: 7.46}
    rh = {}
    ok = True
    for c, t in T_RAD.items():
        d = D[t]
        a, b, step, mono = r_half_of(d["e"], d["P"].r)
        rh[f"{c:g}"] = {
            "mine": a,
            "stored_json": d["last"]["r_half"],
            "claimed": stored_rh[c],
            "interpolated": b,
            "reader_step_at_r_half": step,
            "reader_step_percent": 100.0 * step / a,
        }
        ok &= abs(a - d["last"]["r_half"]) < 1e-9 and abs(a - stored_rh[c]) < 0.006
    P03 = D[T_RAD[0.0]]["P"]
    seed_rad, ext03 = make_seed(P03, "rad")
    es = sum(energy_parts(seed_rad, np.full(P03.r.shape, G8), P03))
    rh["seed"] = {"mine": r_half_of(es, P03.r)[0], "claimed": 4.44}
    ok &= abs(rh["seed"]["mine"] - 4.44) < 0.006
    add(
        "C2.r_half_values",
        "r_half end values 17.02, 16.48, 15.93, 13.97, 11.30, 7.46 for c = 0 to 3e-3 and 4.44 "
        "for the seed",
        "own per-cell density, cells grouped by radius, first radius whose ball holds half",
        rh,
        "each within 0.006 of the claim and equal to the stored JSON value",
        "any recomputed radius off by more than 0.006",
        ok,
    )
    neg = {}
    for t, d in D.items():
        neg[t] = {
            "min_cell_total": float(d["e"].min()),
            "min_cell_lin_part": float(d["lin"].min()),
            "cells_with_negative_lin": int((d["lin"] < 0).sum()),
            "cumulative_monotone": r_half_of(d["e"], d["P"].r)[3],
        }
    add(
        "C2.negative_density",
        "the per-cell density including -c L is never negative in the end fields, so the "
        "cumulative read is monotone",
        "min over cells of my density and of the -c L part alone; monotonicity of the shell "
        "cumulative sum",
        neg,
        "min total density >= 0 and a monotone cumulative sum on every row",
        "a negative cell or a non-monotone cumulative sum on any row",
        all(v["min_cell_total"] >= 0 and v["cumulative_monotone"] for v in neg.values()),
    )
    uni = r_half_of(np.ones(P03.r.shape), P03.r)
    uni_free = r_half_of(P03.free.astype(float), P03.r)
    ceil = {
        "uniform_whole_cube": uni[0],
        "uniform_whole_cube_continuum": (3.0 * 48.0**3 / (8.0 * np.pi)) ** (1.0 / 3.0),
        "uniform_free_region_only": uni_free[0],
        "uniform_free_region_continuum": (3.0 * 42.0**3 / (8.0 * np.pi)) ** (1.0 / 3.0),
    }
    rho2 = r_half_of(1.0 / P03.r**2, P03.r)[0]
    ceil["density_r_minus_2"] = rho2
    add(
        "C2.reader_ceiling",
        "the reader has a finite ceiling in this cube: a uniform density reads about 23.6 "
        "(whole cube) and about 20.7 (free region only)",
        "r_half of a constant density on the same 32^3 cells, and on the free cells alone",
        ceil,
        "lattice values within 0.3 of the continuum 23.63 and 20.68",
        "the lattice ceiling differs from the continuum estimate by more than 0.3",
        abs(uni[0] - ceil["uniform_whole_cube_continuum"]) < 0.3
        and abs(uni_free[0] - ceil["uniform_free_region_continuum"]) < 0.3,
    )
    hh = energy_parts(ext03, np.full(P03.r.shape, G8), P03)[0]
    dec = {}
    for c in (0.0, 3e-3):
        d = D[T_RAD[c]]
        r = d["P"].r
        inside = r < 6.0
        pin = ~d["P"].free
        tot = float(d["e"].sum())
        ext_only = r_half_of(np.where(inside, 0.0, d["e"]), r)[0]
        hh_only = r_half_of(np.where(inside, 0.0, hh), r)[0]
        dec[f"{c:g}"] = {
            "E": tot,
            "frac_inside_r6": float(d["e"][inside].sum() / tot),
            "frac_outside_r6": float(d["e"][~inside].sum() / tot),
            "frac_pinned_shell": float(d["e"][pin].sum() / tot),
            "inside_r6_total": float(d["e"][inside].sum()),
            "seed_inside_r6_total": float(es[inside].sum()),
            "outside_r6_total": float(d["e"][~inside].sum()),
            "outside_r6_curv": float(d["ec"][~inside].sum()),
            "outside_r6_V4": float(d["v4"][~inside].sum()),
            "outside_r6_lin": float(d["lin"][~inside].sum()),
            "unrelaxed_uniaxial_hedgehog_curv_outside_r6": float(hh[~inside].sum()),
            "r_half_of_the_exterior_energy_alone": ext_only,
            "r_half_of_the_hedgehog_tail_alone": hh_only,
        }
    add(
        "C2.decomposition",
        "at c = 0 r_half reads the exterior (halo), not the core: most of E sits outside r = 6 "
        "and it is not distributed like the Coulomb tail of the uniaxial hedgehog",
        "E inside / outside r = 6 by part; the half-energy radius of the exterior energy alone "
        "against that of the curvature of the unrelaxed uniaxial hedgehog on the same cells",
        dec,
        "outside fraction above 0.5 at c = 0 and an exterior-only radius more than 20 percent "
        "above the hedgehog tail's",
        "the core holds half of E or more at c = 0, or the exterior energy is distributed like "
        "the 1 / r^4 tail",
        dec["0"]["frac_outside_r6"] > 0.5
        and dec["0"]["r_half_of_the_exterior_energy_alone"]
        > 1.2 * dec["0"]["r_half_of_the_hedgehog_tail_alone"],
    )
    floor = {}
    for c, t in T_RAD.items():
        d = D[t]
        inside = d["P"].r < 6.0
        ref = r_half_of(np.where(inside, d["e"], hh), d["P"].r)
        act = r_half_of(d["e"], d["P"].r)
        floor[f"{c:g}"] = {
            "r_half": act[0],
            "r_half_interpolated": act[1],
            "r_half_with_exterior_replaced_by_bare_hedgehog_tail": ref[0],
            "same_interpolated": ref[1],
            "excess_over_halo_free_reference_percent": 100.0 * (act[1] - ref[1]) / ref[1],
        }
    low = [k for k, v in floor.items() if abs(v["excess_over_halo_free_reference_percent"]) < 10]
    add(
        "C2.coulomb_floor",
        "r_half responds to the halo at every c of the scan (so a slope of r_half against c is "
        "a slope of the halo radius)",
        "reference read: the row's own density inside r = 6 plus the curvature density of the "
        "unrelaxed uniaxial hedgehog outside (no halo at all); the row's r_half must exceed it",
        {"rows": floor, "rows_within_10_percent_of_the_halo_free_reference": low},
        "every row more than 10 percent above its halo-free reference",
        "a row whose r_half is within 10 percent of the halo-free reference: there the read is "
        "set by the core and the Coulomb tail, and it floors instead of following the halo",
        not low,
    )

    # ---------------- C3: the spectrum and the halo ----------------
    d0 = D[T_RAD[0.0]]
    h = d0["P"].h
    sh0 = {f"{R:g}": shell_means(d0["lam"], d0["P"].r, R, 0.75 * h) for R in SHELLS}
    inner = [sh0[f"{R:g}"]["mean"] for R in SHELLS if R <= 15.0]
    rngs = [[min(x[a] for x in inner), max(x[a] for x in inner)] for a in range(3)]
    bands = [(0.14, 0.18), (0.44, 0.50), (0.96, 0.98)]
    in_band = all(
        bands[a][0] - 0.005 <= rngs[a][0] and rngs[a][1] <= bands[a][1] + 0.005 for a in range(3)
    )
    add(
        "C3.c0_spectrum_bands",
        "at c = 0 the shell-mean spectrum is about (0.14 to 0.18, 0.44 to 0.50, 0.96 to 0.98) "
        "and flat in r from 1.5 to 15",
        "own shells |r - R| < 0.75 h, R in 1.5, 3, 4.5, 6, 9, 12, 15; min and max over R",
        {"range_small_mid_top": rngs, "shells": sh0},
        "every shell mean inside the quoted band (0.005 rounding allowance)",
        "a shell mean between r 1.5 and r 15 outside its band",
        in_band,
    )
    spread = {}
    eps0 = 0.5 * (d0["lam"][..., 1] - d0["lam"][..., 0])
    gap0 = d0["lam"][..., 2] - d0["lam"][..., 1]
    for R in SHELLS:
        if R > 15.0:
            continue
        m = np.abs(d0["P"].r - R) < 0.75 * h
        sm = sh0[f"{R:g}"]
        spread[f"{R:g}"] = {
            "std_small_mid_top": sm["std"],
            "std_mid_over_gap_mid_small": sm["std"][1] / (sm["mean"][1] - sm["mean"][0]),
            "eps_percentiles_5_50_95": [float(x) for x in np.percentile(eps0[m], [5, 50, 95])],
            "gap_top_mid_min_in_shell": float(gap0[m].min()),
        }
    wsp = max(v["std_mid_over_gap_mid_small"] for v in spread.values())
    add(
        "C3.c0_angular_spread",
        "WORDING: 'biaxial with a spectrum about (...) flat in r' describes the field, not only "
        "its shell average",
        "standard deviation of each eigenvalue over the cells of a shell at c = 0 against the "
        "mid - small gap; the 5 / 50 / 95 percentiles of eps in the shell; the smallest top - "
        "mid gap in the shell",
        spread,
        "within-shell std of the mid eigenvalue under 25 percent of the mid - small gap",
        "the within-shell scatter is 25 percent of the gap or more (the shell mean then hides "
        "an angular structure, here a disclination network of the top eigenvector, see C5)",
        wsp < 0.25,
    )
    r22, d_out, d_core, near_out, near_core = None, None, None, None, None
    if os.path.exists(R22S_FIELD):
        M22 = np.load(R22S_FIELD)["M"]
        l22 = np.linalg.eigvalsh(M22[..., 1:, 1:])
        r22 = {f"{R:g}": shell_means(l22, d0["P"].r, R, 0.75 * h)["mean"] for R in SHELLS}
        q22 = (0.17, 0.45, 0.98)

        def dmax_of(Rs, other):
            return max(abs(r22[f"{R:g}"][a] - other(R, a)) for R in Rs for a in range(3))

        d_out = dmax_of((6.0, 9.0, 12.0, 15.0), lambda R, a: sh0[f"{R:g}"]["mean"][a])
        d_core = dmax_of((1.5, 3.0, 4.5), lambda R, a: sh0[f"{R:g}"]["mean"][a])
        near_out = dmax_of((6.0, 9.0, 12.0, 15.0), lambda R, a: q22[a])
        near_core = dmax_of((1.5, 3.0, 4.5), lambda R, a: q22[a])
    r22_val = {
        "r22_1s_shells": r22,
        "c0_row_vs_r22_max_abs_diff_r6_to_15": d_out,
        "c0_row_vs_r22_max_abs_diff_r1.5_to_4.5": d_core,
        "r22_vs_quoted_triple_r6_to_15": near_out,
        "r22_vs_quoted_triple_r1.5_to_4.5": near_core,
    }
    add(
        "C3.reproduces_r22_halo",
        "the c = 0 row reproduces the earlier R22 interior (0.17, 0.45, 0.98) from a fresh seed, "
        "on the shells r 6 to 15",
        "shell means (0.75 h) of the R22-1s end field (read only) against the c = 0 row and "
        "against the quoted triple",
        r22_val,
        "both within 0.02 on r 6 to 15",
        "a shell mean off by more than 0.02 on r 6 to 15, or the R22 field is missing",
        d_out is not None and d_out < 0.02 and near_out < 0.02,
    )
    add(
        "C3.reproduces_r22_core",
        "WORDING: the same reproduction holds for the whole 'interior', the shells r 1.5 to 4.5 "
        "included",
        "as above on the shells r 1.5, 3, 4.5",
        r22_val,
        "both within 0.02 on r 1.5 to 4.5",
        "the two fields, or the R22 field and the quoted triple, differ by more than 0.02 inside "
        "r 4.5 (the triple is then a halo value, not an interior value)",
        d_core is not None and d_core < 0.02 and near_core < 0.02,
    )
    quoted = {
        6.0: (0.289, 0.456, 0.944),
        9.0: (0.301, 0.364, 0.979),
        12.0: (0.3015, 0.317, 0.994),
        15.0: (0.301, 0.304, 0.9985),
    }
    d3 = D[T_RAD[3e-3]]
    q3, worst3, sens3 = {}, 0.0, 0.0
    for R, qv in quoted.items():
        rec = {}
        for nm, hw in (("0.5h", 0.5 * h), ("0.75h", 0.75 * h), ("1.0h", 1.0 * h)):
            rec[nm] = shell_means(d3["lam"], d3["P"].r, R, hw)["mean"]
        rec["quoted"] = list(qv)
        q3[f"{R:g}"] = rec
        worst3 = max(worst3, max(abs(rec["0.75h"][a] - qv[a]) for a in range(3)))
        sens3 = max(sens3, max(abs(rec["0.5h"][a] - rec["1.0h"][a]) for a in range(3)))
    add(
        "C3.c3e-3_shell_means",
        "at c = 3e-3 the shell means are (0.289, 0.456, 0.944) at r 6, (0.301, 0.364, 0.979) at "
        "r 9, (0.3015, 0.317, 0.994) at r 12, (0.301, 0.304, 0.9985) at r 15",
        "own shells of half-width 0.5 h, 0.75 h (the record's) and 1.0 h",
        {"shells": q3, "worst_abs_diff_0.75h": worst3, "width_sensitivity": sens3},
        "within 0.002 of the quoted values at half-width 0.75 h",
        "any quoted shell mean off by more than 0.002",
        worst3 < 0.002,
    )
    mono, cut, wig = {}, {}, {}
    ok_mono, ok_wig = True, True
    for c, t in T_RAD.items():
        d = D[t]
        pr = eps_profile(d["lam"], d["P"].r, h, 0.5 * h)
        sel = pr[(pr[:, 0] >= 4.5) & (pr[:, 0] <= 21.0)]
        up = np.diff(sel[:, 1]) / sel[:-1, 1]
        wig[f"{c:g}"] = {
            "largest_upward_step_rel": float(up.max()),
            "upward_steps": int((up > 0).sum()),
            "upward_steps_over_5_percent": int((up > 0.05).sum()),
            "r_of_largest_upward_step": float(sel[int(np.argmax(up)), 0]),
        }
        eps_c = 0.5 * (d["lam"][..., 1] - d["lam"][..., 0])
        cont = np.array(
            [
                [R, float(eps_c[np.abs(d["P"].r - R) < 0.5 * h].mean())]
                for R in np.arange(4.5, 21.1, 1.5)
            ]
        )
        upc = np.diff(cont[:, 1]) / cont[:-1, 1]
        mono[f"{c:g}"] = {
            "eps_on_contiguous_shells": [[float(a_), float(b_)] for a_, b_ in cont],
            "largest_upward_step_rel": float(upc.max()),
            "r_of_largest_upward_step": float(cont[int(np.argmax(upc)), 0]),
        }
        e6 = float(np.interp(6.0, pr[:, 0], pr[:, 1]))
        e18 = float(np.interp(18.0, pr[:, 0], pr[:, 1]))
        cut[f"{c:g}"] = {
            "eps_r6": e6,
            "eps_r18": e18,
            "ratio": e18 / e6,
            "visibly_cut": bool(e18 / e6 < 0.25),
        }
        if c > 0:
            ok_mono &= mono[f"{c:g}"]["largest_upward_step_rel"] <= 0.05
            ok_wig &= wig[f"{c:g}"]["upward_steps_over_5_percent"] == 0
        d["prof"] = pr
    add(
        "C3.eps_monotone",
        "with the stiffness on, the split eps(r) decays with r (every c > 0, rad seed)",
        "own shells: contiguous, non-overlapping, width h = 1.5, centers 4.5 to 21",
        mono,
        "no upward step above 5 percent on any c > 0 row",
        "a c > 0 row with an upward step of eps(r) above 5 percent between neighboring shells",
        ok_mono,
    )
    outer = {
        k: max(
            (b_[1] - a_[1]) / a_[1]
            for a_, b_ in zip(
                v["eps_on_contiguous_shells"][:-1], v["eps_on_contiguous_shells"][1:]
            )
            if a_[0] >= 10.5
        )
        for k, v in mono.items()
    }
    add(
        "C3.eps_decays_outer",
        "the decay of eps(r) holds in the outer halo (r >= 10.5) on every row",
        "the same contiguous shells restricted to centers 10.5 to 21",
        {"largest_upward_step_rel_r_ge_10.5": outer},
        "no upward step above 5 percent",
        "an upward step above 5 percent beyond r = 10.5",
        all(v <= 0.05 for v in outer.values()),
    )
    add(
        "C3.eps_profile_smooth",
        "WORDING: the stored profile (shells of width h on a 0.75 grid, overlapping) is a "
        "smooth decaying curve that a reader can interpolate (the r_eps reader and any fit do)",
        "the same profile rebuilt with my shells, 4.5 <= r <= 21; relative upward steps",
        wig,
        "no upward step above 5 percent on any c > 0 row",
        "an upward step above 5 percent (lattice sampling of a strongly angle-dependent eps)",
        ok_wig,
    )
    vis = [k for k, v in cut.items() if v["visibly_cut"]]
    add(
        "C3.halo_cut_in_box",
        "the halo is cut inside the L 48 box only at the larger c (the record treats c = 1e-3 "
        "and 3e-3 as cut and c <= 3e-4 as box-limited)",
        "eps(18) / eps(6) against 0.25 (half of the free r^-0.618 value 0.507)",
        {"rows": cut, "visibly_cut": vis},
        "exactly the rows c = 1e-3 and 3e-3 fall under 0.25",
        "the set of visibly cut rows is not {0.001, 0.003}",
        sorted(vis) == ["0.001", "0.003"],
    )

    # ---------------- C4: the halo shape ----------------
    rr = np.linspace(2.0, 12.0, 4001)
    beta_t = 3e-3
    yk = (
        np.sqrt(rr)
        * kve(NU, 0.5 * np.sqrt(beta_t) * rr**2)
        * np.exp(-0.5 * np.sqrt(beta_t) * rr**2)
    )
    ypp = np.gradient(np.gradient(yk, rr), rr)
    res_ode = np.abs(ypp - (1.0 / rr**2 + beta_t * rr**2) * yk)[50:-50] / np.abs(yk[50:-50])
    add(
        "C4.knu_solves_ode",
        "sqrt(r) K_nu(sqrt(beta) r^2 / 2) with nu = sqrt(5) / 4 solves eps'' = eps / r^2 + "
        "beta r^2 eps",
        "second finite difference of the closed form on a fine grid, beta = 3e-3",
        {"max_rel_residual": float(res_ode.max())},
        "relative ODE residual under 1e-4",
        "the closed form does not satisfy the ODE",
        float(res_ode.max()) < 1e-4,
    )
    fits = {}
    for c in (1e-3, 3e-3):
        pr = D[T_RAD[c]]["prof"]
        sel = pr[(pr[:, 0] >= 5.0) & (pr[:, 0] <= 18.0) & (pr[:, 1] > 0)]
        fits[f"{c:g}"] = fit_candidates(sel[:, 0], sel[:, 1])
        fits[f"{c:g}"]["points"] = len(sel)
        fits[f"{c:g}"]["eps_first_last"] = [float(sel[0, 1]), float(sel[-1, 1])]
    rc = [fits[k]["knu_cutoff"]["R_cut_log_fit"] for k in ("0.001", "0.003")]
    rl = [fits[k]["knu_cutoff"]["R_cut_lin_fit"] for k in ("0.001", "0.003")]
    fits["knu_scale_against_c"] = {
        "R_cut_log_fit_c1e-3_c3e-3": rc,
        "implied_exponent_log_fit": float(np.log(rc[1] / rc[0]) / np.log(3.0)),
        "R_cut_lin_fit_c1e-3_c3e-3": rl,
        "implied_exponent_lin_fit": float(np.log(rl[1] / rl[0]) / np.log(3.0)),
        "note": "two points and a poor K_nu fit at c = 1e-3: informs, settles nothing",
    }
    better = all(
        f["knu_cutoff"]["rms_log"] < f["yukawa"]["rms_log"]
        and f["knu_cutoff"]["rms_lin_rel_to_max"] < f["yukawa"]["rms_lin_rel_to_max"]
        for k, f in fits.items()
        if k != "knu_scale_against_c"
    )
    add(
        "C4.knu_vs_yukawa",
        "the measured profiles at c = 1e-3 and 3e-3 are better described by the K_nu cutoff "
        "than by exp(-r / R) / r",
        "two-parameter least squares (amplitude and one scale) on 5 <= r <= 18, residuals in "
        "log space and in linear space; a pure power law as the third candidate",
        fits,
        "K_nu residual below the Yukawa residual on both rows in both metrics",
        "the Yukawa form fits as well or better on either row in either metric",
        better,
    )

    good = {k: f["knu_cutoff"]["rms_log"] for k, f in fits.items() if k != "knu_scale_against_c"}
    add(
        "C4.knu_describes_profile",
        "WORDING: 'better described by the K_nu cutoff' also means 'described by it'",
        "RMS log residual of the two-parameter K_nu fit on 5 <= r <= 18",
        good,
        "RMS log residual under 0.15 (a 15 percent scatter) on both rows",
        "a row where the best of the three candidates still misses by more than 15 percent RMS",
        all(v < 0.15 for v in good.values()),
    )

    # ---------------- C5: the degree ----------------
    deg = {}
    ok_deg = True
    for t, d in D.items():
        P = d["P"]
        top = d["vec"][..., :, 2]
        gap = d["lam"][..., 2] - d["lam"][..., 1]
        rec = {}
        for hv in (6.0, 9.0, 12.0):
            k = int(round(hv / P.h))
            c0 = P.n // 2
            mine = degree_read(top, gap, P, c0 - k, c0 + k)
            cen = degree_read(top, gap, P, c0 - k - 1, c0 + k)
            st = d["last"]["degree_top"][f"{hv:g}"]
            mine["stored_degree_conflicts"] = st
            mine["centered_cube_frustrated"] = cen["frustrated_plaquettes"]
            mine["centered_cube_degree"] = cen["degree_bfs"]
            rec[f"{hv:g}"] = mine
            same_class = (st[1] == 0) == (mine["frustrated_plaquettes"] == 0)
            ok_deg &= same_class
            if st[1] == 0 and mine["frustrated_plaquettes"] == 0:
                ok_deg &= abs(abs(mine["degree_bfs"]) - abs(st[0])) < 1e-6
                ok_deg &= mine["conflicts_bfs"] == 0
        rec["volume"] = volume_disclinations(top, P)
        deg[t] = rec
    add(
        "C5.degree_all_rows",
        "the stored surface-oriented degree reads hold: zero conflicts exactly where the top "
        "eigenvector is orientable on the cube, and |degree| = 1 there",
        "own reader: frustrated elementary squares on the cube surface (gauge invariant), "
        "breadth-first flood-fill orientation, r-hat orientation, solid-angle degree with both "
        "quad splittings and geometric outward ordering; also a cube centered on the box",
        deg,
        "stored conflicts are zero if and only if my frustrated-square count is zero, and "
        "|degree| agrees to 1e-6 where defined",
        "a cube the record calls clean has a frustrated square, a cube the record calls "
        "conflicted has none, or the degrees differ",
        ok_deg,
    )
    claim5 = {
        T_RAD[3e-3]: (True, True, True),
        T_RAD[1e-3]: (True, True, True),
        T_RAD[0.0]: (False, False, True),
        T_RAD[3e-5]: (False, None, None),
        T_RAD[1e-4]: (False, None, None),
    }
    got5, ok5 = {}, True
    for t, want in claim5.items():
        g = [deg[t][k]["frustrated_plaquettes"] == 0 for k in ("6", "9", "12")]
        dg = [deg[t][k]["degree_bfs"] for k in ("6", "9", "12")]
        got5[t] = {"orientable_r6_r9_r12": g, "degree_signed_outward": dg}
        for a in range(3):
            if want[a] is not None:
                ok5 &= g[a] == want[a]
                if want[a]:
                    ok5 &= abs(abs(dg[a]) - 1.0) < 1e-6
    add(
        "C5.claimed_pattern",
        "clean unit degree on all three cubes at c = 3e-3 and 1e-3; at c = 0 conflicts on r 6 "
        "and r 9 and a unit degree on r 12; at c = 3e-5 and 1e-4 conflicts on r 6",
        "own reader as above, rad seed rows",
        got5,
        "the orientable / non-orientable pattern and the unit degrees as claimed",
        "any cube of the claimed pattern reads differently with my reader",
        ok5,
    )
    cnt, flips = {}, {}
    for t in D:
        for k in ("6", "9", "12"):
            q = deg[t][k]
            if q["stored_degree_conflicts"][1] > 0:
                cnt[f"{t} r{k}"] = {
                    "stored_conflicts": q["stored_degree_conflicts"][1],
                    "my_flood_fill_conflicts": q["conflicts_bfs"],
                    "my_rhat_conflicts": q["conflicts_rhat"],
                    "frustrated_squares_gauge_invariant": q["frustrated_plaquettes"],
                }
            a_ = q["frustrated_plaquettes"] == 0
            b_ = q["centered_cube_frustrated"] == 0
            if a_ != b_:
                flips[f"{t} r{k}"] = {
                    "record_cube_frustrated": q["frustrated_plaquettes"],
                    "box_centered_cube_frustrated": q["centered_cube_frustrated"],
                    "box_centered_cube_degree": q["centered_cube_degree"],
                }
    add(
        "C5.conflict_counts",
        "WORDING: the stored conflict counts (28 and 6 at c = 0; 35, 63, 89 on row A) are "
        "numbers of the field",
        "the same cubes oriented by a breadth-first flood fill and by r-hat, against the "
        "gauge-invariant count of frustrated elementary squares (disclination piercings)",
        cnt,
        "my flood-fill count equals the stored count on every conflicted cube",
        "a conflicted cube where another spanning tree gives another count (the count is then "
        "a property of the tree; only zero / nonzero and the piercing count are invariant)",
        all(v["stored_conflicts"] == v["my_flood_fill_conflicts"] for v in cnt.values()),
    )
    add(
        "C5.cube_centering",
        "WORDING: the clean / conflicted reading of a cube does not depend on the reader's cube "
        "being centered on cell n // 2 (half a cell, 0.75, off the box center)",
        "the same read on the cube centered on the box center (one cell wider on the low side)",
        flips,
        "no cube changes class",
        "a cube that is conflicted on the record's cube and clean on the box-centered one, or "
        "the reverse",
        not flips,
    )
    ext = {f"{c:g}": deg[T_RAD[c]]["volume"] for c in c_scan}
    add(
        "C5.disclination_extent",
        "the clean cubes at c = 1e-3 and 3e-3 mean the line defects of the top eigenvector are "
        "confined inside r = 6 there, and reach beyond r = 6 for c <= 3e-4",
        "frustrated elementary squares of the top-eigenvector line field over the whole free "
        "region (all three plane orientations), with their radii",
        ext,
        "r_max under 6 at c = 1e-3 and 3e-3, above 6 for c <= 3e-4",
        "frustrated squares beyond r = 6 at the two large c, or none beyond r = 6 at the small c",
        all((ext[f"{c:g}"]["r_max"] or 0) < 6.0 for c in (1e-3, 3e-3))
        and all((ext[f"{c:g}"]["r_max"] or 0) > 6.0 for c in (0.0, 3e-5, 1e-4, 3e-4)),
    )
    signs = {t: [deg[t][k]["degree_bfs"] for k in ("6", "9", "12")] for t in D}
    ok_sign = all(v is None or v > 0 for vs in signs.values() for v in vs)
    add(
        "C5.degree_sign",
        "the mixed signs of the stored degrees (for example -1, +1, +1 at c = 1e-3) are the "
        "arbitrary root sign of the reader and not a physical sign change between cubes",
        "degree with the global sign fixed by mean (v . r-hat) > 0 on each cube",
        signs,
        "every defined degree is +1 once oriented outward",
        "a cube with degree -1 relative to the outward orientation",
        ok_sign,
    )

    # ---------------- C6: convergence and uniqueness ----------------
    imgs = images()
    da = D[T_RAD[3e-4]]
    inv = {}
    base = float(da["e"].sum())
    for perm, sg in imgs:
        Si, mi = apply_image(da["S"], da["m00"], perm, sg)
        Ei = float(sum(energy_parts(Si, mi, da["P"])).sum())
        inv[str((perm, sg))] = abs(Ei - base) / base
    exact = [k for k, v in inv.items() if v < 1e-12]
    n_exact_expected = sum(is_lattice_exact(p, s) for p, s in imgs)
    rest = [v for k, v in inv.items() if v >= 1e-12]
    seed_inv = max(
        float(np.abs(apply_image(seed_rad, np.full(P03.r.shape, G8), p, s)[0] - seed_rad).max())
        for p, s in imgs
    )
    add(
        "C6.twelve_symmetries",
        "the one-sided fwd / bwd averaged functional has only 12 exact lattice symmetries "
        "(axis permutations times full inversion)",
        "energy of the 48 cubic images of the c = 3e-4 end field; the analytic hedgehog seed "
        "must map to itself under all 48 (a test of my image code)",
        {
            "exact_images": len(exact),
            "min_rel_change_of_the_other_36": float(min(rest)) if rest else None,
            "max_rel_change_of_the_other_36": float(max(rest)) if rest else None,
            "seed_invariance_max_abs": seed_inv,
        },
        "exactly 12 images leave E unchanged to 1e-12 and they are the predicted ones",
        "the count of exact images is not 12, or my image code does not fix the seed",
        len(exact) == n_exact_expected == 12
        and seed_inv < 1e-13
        and all((inv[str((p, s))] < 1e-12) == is_lattice_exact(p, s) for p, s in imgs),
    )
    pairs = {
        "two_sided_c3e-4": (T_RAD[3e-4], T_TWO),
        "seed_c1e-4": (T_RAD[1e-4], T_BIA[1e-4]),
        "seed_c1e-3": (T_RAD[1e-3], T_BIA[1e-3]),
    }
    pr_out = {}
    for nm, (ta, tb) in pairs.items():
        A, Bf = D[ta], D[tb]
        P = A["P"]
        r = P.r
        regions = {
            "all_free": P.free,
            "r_lt_6": r < 6.0,
            "r_6_12": (r >= 6.0) & (r < 12.0),
            "r_lt_12": r < 12.0,
            "r_ge_12_free": (r >= 12.0) & P.free,
        }

        def rms(X, m):
            return float(np.sqrt(np.mean(np.sum(X[m] ** 2, axis=(-1, -2)))))

        best12, best48 = None, None
        for perm, sg in imgs:
            Si, mi = apply_image(Bf["S"], Bf["m00"], perm, sg)
            v = rms(A["S"] - Si, P.free)
            if best48 is None or v < best48[0]:
                best48 = (v, (perm, sg), Si, mi)
            if is_lattice_exact(perm, sg) and (best12 is None or v < best12[0]):
                best12 = (v, (perm, sg), Si, mi)
        ident = rms(A["S"] - Bf["S"], P.free)
        halo_scale = rms(A["S"] - make_seed(P, "rad")[1], P.free)
        lam_d = np.sqrt(np.mean(np.sum((A["lam"] - Bf["lam"]) ** 2, axis=-1)[P.free]))
        Sb = best12[2]
        lam_b = np.linalg.eigvalsh(Sb)
        reg = {
            k: {
                "rms_S_diff_best12": rms(A["S"] - Sb, m),
                "rms_sorted_eigenvalue_diff_best12": float(
                    np.sqrt(np.mean(np.sum((A["lam"] - lam_b) ** 2, axis=-1)[m]))
                ),
            }
            for k, m in regions.items()
        }
        # a straight path between the two end fields, M_00 re-solved
        path = []
        for tt in np.linspace(0.0, 1.0, 11):
            St = (1.0 - tt) * A["S"] + tt * Sb
            mt = A["m00"].copy()
            mt[P.free] = solve_u(St[P.free], P, A["m00"][P.free])
            path.append(float(sum(energy_parts(St, mt, P)).sum()))
        gdot = float(np.sum((A["G"] * (Sb - A["S"]))[P.free]))
        shell_cmp = {
            f"{R:g}": [
                shell_means(A["lam"], r, R, 0.75 * P.h)["mean"],
                shell_means(Bf["lam"], r, R, 0.75 * P.h)["mean"],
            ]
            for R in (3.0, 6.0, 9.0, 12.0, 15.0)
        }
        rot = best_rotation(A["S"], Bf["S"], P, seed=len(pr_out) + 7)
        log(f"rotation search {nm}: {rot}")
        rha, rhb = r_half_of(A["e"], r), r_half_of(Bf["e"], r)
        pr_out[nm] = {
            "labels": [A["row"]["label"], Bf["row"]["label"]],
            "E": [float(A["e"].sum()), float(Bf["e"].sum())],
            "E_rel_diff_percent": 100.0 * abs(A["e"].sum() - Bf["e"].sum()) / Bf["e"].sum(),
            "r_half": [rha[0], rhb[0]],
            "r_half_interpolated": [rha[1], rhb[1]],
            "r_half_rel_diff_percent": 100.0 * abs(rha[0] - rhb[0]) / rhb[0],
            "r_half_interp_rel_diff_percent": 100.0 * abs(rha[1] - rhb[1]) / rhb[1],
            "continuous_rotation_search_r_lt_18": rot,
            "rms_S_diff_identity_free": ident,
            "rms_S_diff_min_over_12": best12[0],
            "image_12": str(best12[1]),
            "rms_S_diff_min_over_48": best48[0],
            "image_48": str(best48[1]),
            "rms_of_field_minus_uniaxial_hedgehog": halo_scale,
            "rms_sorted_eigenvalue_diff_identity": float(lam_d),
            "regions": reg,
            "shell_spectra_a_b": shell_cmp,
            "straight_path_E": path,
            "straight_path_barrier_above_higher_end": max(path) - max(path[0], path[-1]),
            "dE_dt_at_a_from_my_gradient": gdot,
        }
    e_pc = [v["E_rel_diff_percent"] for v in pr_out.values()]
    r_pc = [v["r_half_rel_diff_percent"] for v in pr_out.values()]
    add(
        "C6.three_percent",
        "the radius is reproducible to about 3 percent across start fields and the energy to "
        "about 2 to 3 percent",
        "own energies and own r_half on the three pairs of end fields",
        {"E_percent": e_pc, "r_half_percent": r_pc, "pairs": pr_out},
        "every radius and every energy difference under 3.5 percent",
        "a pair differing by 3.5 percent or more in the radius or in the energy",
        max(r_pc) < 3.5 and max(e_pc) < 3.5,
    )
    ri_pc = [v["r_half_interp_rel_diff_percent"] for v in pr_out.values()]
    add(
        "C6.three_percent_reader_free",
        "WORDING: the same 3 percent holds when the one-cell quantization of the r_half reader "
        "is removed, and 'energy to about 2 to 3 percent' describes the AT_GATE pairs",
        "r_half linearly interpolated between the bracketing cell radii; the energy differences "
        "of the two pairs where both rows are AT_GATE",
        {
            "r_half_interp_percent": ri_pc,
            "E_percent_two_sided_and_seed_c1e-3": [e_pc[0], e_pc[2]],
            "E_percent_pair_with_a_FALLING_row": e_pc[1],
        },
        "interpolated radius differences under 3.5 percent and both AT_GATE energy differences "
        "inside 2 to 3 percent",
        "an interpolated radius difference of 3.5 percent or more, or an AT_GATE energy "
        "difference outside 2 to 3 percent",
        max(ri_pc) < 3.5 and all(2.0 <= x <= 3.0 for x in (e_pc[0], e_pc[2])),
    )
    rot_chk = {k: v["continuous_rotation_search_r_lt_18"] for k, v in pr_out.items()}
    lat = max(
        float(
            np.abs(
                rotated_copy(da["S"], image_matrix(p_, s_), da["P"])
                - apply_image(da["S"], da["m00"], p_, s_)[0]
            ).max()
        )
        for p_, s_ in imgs
        if np.linalg.det(image_matrix(p_, s_)) > 0
    )
    add(
        "C6.not_a_rotated_copy",
        "the paired end fields are not one texture in two orientations (a continuous rotation, "
        "which the lattice and the cubic pin break only weakly, does not map one onto the other)",
        "min over SO(3), with and without the inversion, of the RMS of S_a - R S_b(R^T x) R^T on "
        "r < 18: 24 lattice rotations plus 150 random ones, the best 4 refined by Nelder-Mead; "
        "my interpolated rotation reproduces my exact lattice images to the value shown",
        {"pairs": rot_chk, "lattice_rotation_code_check_max_abs": lat},
        "the best rotation leaves more than half of the identity RMS on every pair",
        "a rotation that brings the RMS under half of its identity value (then the fields are "
        "the same object turned), or my rotation code failing its lattice check",
        lat < 1e-12
        and all(
            min(v["proper"]["rms_min"], v["with_inversion"]["rms_min"])
            > 0.5 * v["proper"]["rms_identity"]
            for v in rot_chk.values()
        ),
    )
    not_same = all(
        v["rms_S_diff_min_over_48"] > 0.1 * v["rms_of_field_minus_uniaxial_hedgehog"]
        for v in pr_out.values()
    )
    add(
        "C6.not_the_same_field",
        "the paired end fields are NOT the same field, also after the symmetry images",
        "RMS Frobenius difference of the spatial blocks over free cells, minimized over the 12 "
        "exact and all 48 cubic images, against the RMS departure of the field from the "
        "uniaxial hedgehog",
        {
            k: {
                "min12": v["rms_S_diff_min_over_12"],
                "min48": v["rms_S_diff_min_over_48"],
                "scale": v["rms_of_field_minus_uniaxial_hedgehog"],
            }
            for k, v in pr_out.items()
        },
        "min over 48 images above 10 percent of the halo scale on every pair",
        "a pair that coincides (under 10 percent of the halo scale) after some image",
        not_same,
    )
    where = {
        k: {
            "core_r_lt_6": v["regions"]["r_lt_6"]["rms_S_diff_best12"],
            "r_6_12": v["regions"]["r_6_12"]["rms_S_diff_best12"],
            "r_ge_12": v["regions"]["r_ge_12_free"]["rms_S_diff_best12"],
            "eig_core": v["regions"]["r_lt_6"]["rms_sorted_eigenvalue_diff_best12"],
            "eig_r_ge_12": v["regions"]["r_ge_12_free"]["rms_sorted_eigenvalue_diff_best12"],
        }
        for k, v in pr_out.items()
    }
    add(
        "C6.where_the_difference_lives",
        "the differences between paired end fields are a frame (orientation) difference more "
        "than a spectrum difference",
        "per-region RMS of S_a - S_b against the RMS of the sorted-eigenvalue difference",
        where,
        "in every region the eigenvalue RMS is under half of the full-matrix RMS",
        "the eigenvalue difference carries half or more of the matrix difference in a region",
        all(
            v["eig_core"] < 0.5 * v["core_r_lt_6"] and v["eig_r_ge_12"] < 0.5 * v["r_ge_12"]
            for v in where.values()
        ),
    )
    barrier = {
        k: {
            "path_E": v["straight_path_E"],
            "barrier": v["straight_path_barrier_above_higher_end"],
            "dE_dt_at_a": v["dE_dt_at_a_from_my_gradient"],
        }
        for k, v in pr_out.items()
    }
    add(
        "C6.distinct_stationary_points",
        "WORDING: 'the gate does not select a unique stationary point' needs the two AT_GATE "
        "fields to be separated, not two stops along one slow valley",
        "energy along the straight path between the pair (best exact image, M_00 re-solved); a "
        "straight path between two frame textures is a poor path, so a barrier here shows the "
        "fields are far apart and does not prove two separate minima",
        barrier,
        "a barrier above the higher end on the two pairs where both rows are AT_GATE",
        "the path energy falls monotonically from the higher field to the lower one (then the "
        "higher AT_GATE field is a slow-valley stop and the gate is too loose, not non-unique)",
        all(
            pr_out[k]["straight_path_barrier_above_higher_end"] > 1e-6
            for k in ("two_sided_c3e-4", "seed_c1e-3")
        ),
    )

    # ---------------- C7: the c 3e-3 row and the drift ----------------
    drift = {}
    for t, d in D.items():
        ch = [q for q in d["row"]["chunks"] if q["iters"] > 0]
        last_it = ch[-1]["iters"]
        out = {"label": rows[t]["label"], "label_of_the_end_field_row": d["row"]["label"]}
        out["iters"] = last_it
        for win in (500, 1000):
            w = [q["r_half"] for q in ch if q["iters"] >= last_it - win]
            out[f"drift_{win}"] = float((max(w) - min(w)) / np.mean(w)) if len(w) > 1 else None
        rr_ = r_half_of(d["e"], d["P"].r)
        out["reader_step_percent"] = 100.0 * rr_[2] / rr_[0]
        out["E_drop_last_chunk_over_E"] = ch[-1].get("drop", 0.0) / abs(ch[-1]["E"])
        drift[t] = out
    t3 = T_RAD[3e-3]
    main_ch = [q for q in rows[t3]["chunks"] if 0 < q["iters"] <= 1000]
    w = [q["r_half"] for q in main_ch if q["iters"] >= 500]
    d500 = float((max(w) - min(w)) / np.mean(w))
    add(
        "C7.c3e-3_drift",
        "the c = 3e-3 row met the gate at 1000 iterations with r_half still rising: the drift "
        "over the last 500 iterations exceeds 2 percent",
        "stored chunk reads at iterations 500, 750, 1000 of the main-pool row",
        {"r_half_500_to_1000": w, "drift": d500, "reader_step_percent": drift[t3]},
        "drift above 0.02",
        "drift at or under 0.02",
        d500 > 0.02,
    )
    cont = None
    try:
        Zs = np.load(os.path.join(FIELDS, t3 + "_stage.npz"), allow_pickle=True)
        cont = [
            [q["iters"], q["r_half"], q["E"]]
            for q in json.loads(str(Zs["chunks"]))
            if q["iters"] > 1000
        ]
    except Exception as e:  # noqa: BLE001
        cont = repr(e)
    moved = None
    if isinstance(cont, list) and cont:
        moved = float(max(abs(q[1] - w[-1]) for q in cont) / w[-1])
    rr3 = r_half_of(D[t3]["e"], D[t3]["P"].r)
    add(
        "C7.c3e-3_continuation",
        "WORDING: at the gate (1000 iterations) the c = 3e-3 radius was 'still rising', which "
        "is why the row is held out of a slope fit",
        "chunks past iteration 1000 in the row's stage file (read only, a snapshot if the "
        "continuation pool is still writing); the interpolated r_half of the present end field. "
        "This audit's first pass, run before the continuation pool replaced the 1000-iteration "
        "end field, read an interpolated r_half of 7.3297 on it (not reproducible from the "
        "files any more)",
        {
            "chunks_iters_r_half_E": cont,
            "r_half_at_1000": w[-1],
            "largest_relative_move_after_1000": moved,
            "end_field_r_half": rr3[0],
            "end_field_r_half_interpolated": rr3[1],
            "interpolated_r_half_of_the_1000_iteration_field_first_pass": 7.329710910570718,
            "end_field_row_label": D[t3]["row"]["label"],
            "end_field_row_iters": D[t3]["row"]["iters"],
        },
        "r_half moves by more than 2 percent somewhere in the continuation (or no continued "
        "chunk exists yet)",
        "the continuation leaves r_half within 2 percent of its 1000-iteration value: the "
        "radius was then final at the gate and the 4.2 percent 'drift' was one step of the "
        "reader (see C7.drift_rule_resolution)",
        moved is None or moved > 0.02,
    )
    proj = {}
    for t, d in D.items():
        if rows[t]["label"] != "AT_GATE" or t in FROM_EXTRA:
            continue
        ch = [q for q in d["row"]["chunks"] if q["iters"] > 0 and q.get("chunk_iters", 250) > 5]
        dr = [q["drop"] for q in ch][-5:]
        rec = {"last_drops": dr, "E": ch[-1]["E"]}
        if len(dr) >= 3 and dr[-1] > 0 and dr[0] > 0:
            rho = (dr[-1] / dr[0]) ** (1.0 / (len(dr) - 1))
            rec["drop_ratio_per_chunk"] = float(rho)
            rec["projected_remaining_fall_percent_of_E"] = (
                float(100.0 * dr[-1] * rho / (1.0 - rho) / abs(ch[-1]["E"])) if rho < 1 else None
            )
        proj[t] = rec
    bad = {
        t: v.get("projected_remaining_fall_percent_of_E")
        for t, v in proj.items()
        if v.get("projected_remaining_fall_percent_of_E", 0.0) is None
        or v.get("projected_remaining_fall_percent_of_E", 0.0) > 0.5
    }
    add(
        "C7.energy_convergence_at_gate",
        "WORDING: AT_GATE energies are converged well inside the 1.5 to 3 percent start-field "
        "differences they are compared at",
        "geometric extrapolation of the last five chunk drops of each AT_GATE row (ratio per "
        "chunk, remaining fall as a percent of E)",
        {"rows": proj, "rows_with_projected_fall_above_0.5_percent": bad},
        "projected remaining fall under 0.5 percent of E on every AT_GATE row",
        "an AT_GATE row whose drops are not shrinking (ratio >= 1) or whose projected remaining "
        "fall exceeds 0.5 percent of E",
        not bad,
    )
    others = {
        t: v
        for t, v in drift.items()
        if v["label"] == "AT_GATE" and t != t3 and (v["drift_500"] or 0) > 0.02
    }
    add(
        "C7.other_at_gate_rows",
        "no other AT_GATE row stopped with r_half drifting more than 2 percent over its last "
        "500 iterations",
        "the same drift read on every AT_GATE row (500 and 1000 iteration windows reported)",
        {"offenders_500": list(others), "all": drift},
        "no other AT_GATE row above 0.02 on the 500-iteration window",
        "another AT_GATE row above 0.02",
        not others,
    )
    coarse = {
        t: v["reader_step_percent"]
        for t, v in drift.items()
        if v["label"] == "AT_GATE" and v["reader_step_percent"] > 2.0
    }
    add(
        "C7.drift_rule_resolution",
        "WORDING: a 2 percent drift rule is only meaningful where one step of the r_half reader "
        "(the gap between neighboring cell radii) is under 2 percent",
        "gap between the distinct cell radii bracketing r_half, in percent of r_half",
        {"rows_with_reader_step_above_2_percent": coarse},
        "no AT_GATE row whose reader step exceeds the 2 percent rule",
        "an AT_GATE row where a single reader step is larger than the drift tolerance (a "
        "one-step move then reads as drift and a sub-step move as none)",
        not coarse,
    )

    # ---------------- C8: the delta 0.9224 rows ----------------
    out8 = {}
    for nm, t in (("A_c0", T_A0), ("A_c", T_AC), ("B_c0", T_B0)):
        d = D[t]
        P = d["P"]
        lam = d["lam"]
        vac4 = (1.0 - P.delta) ** 4
        g1 = ((lam[..., 2] - lam[..., 1]) * (lam[..., 2] - lam[..., 0])) ** 2
        seed, _ = make_seed(P, "rad")
        ls = np.linalg.eigvalsh(seed)
        g1s = ((ls[..., 2] - ls[..., 1]) * (ls[..., 2] - ls[..., 0])) ** 2
        ecs = curv_density(seed, P.h)
        sh = {}
        for R in SHELLS:
            m = np.abs(P.r - R) < 0.75 * P.h
            sh[f"{R:g}"] = {
                "spectrum": [float(lam[..., a][m].mean()) for a in range(3)],
                "G1sq_over_vacuum": float(g1[m].mean() / vac4),
                "seed_G1sq_over_vacuum": float(g1s[m].mean() / vac4),
                "curv_density_over_seed": float(d["ec"][m].sum() / ecs[m].sum()),
                "trace": float(lam[m].sum(axis=-1).mean()),
            }
        idx = np.indices(P.r.shape)
        cheb = np.max(np.abs(idx - (P.n - 1) / 2.0), axis=0)
        layers = {
            "pinned_shell": cheb > P.n / 2.0 - 2.0,
            "first_3_free_layers": (cheb <= P.n / 2.0 - 2.0) & (cheb > P.n / 2.0 - 5.0),
            "interior": cheb <= P.n / 2.0 - 5.0,
        }
        tot = float(d["e"].sum())
        rr_ = r_half_of(d["e"], P.r)
        top = d["vec"][..., :, 2]
        rh_ = P.X / np.maximum(P.r, 1e-300)[..., None]
        align = np.abs(np.sum(top * rh_, axis=-1))
        aniso = lam[..., 2] - 0.5 * (lam[..., 1] + lam[..., 0])
        rs = np.arange(1.5, 21.1, 1.5)
        an = np.array([aniso[np.abs(P.r - R) < 0.5 * P.h].mean() for R in rs])
        k_lin = float(np.sum(an * rs) / np.sum(rs * rs))
        ramp = {
            "r": [float(x) for x in rs],
            "anisotropy_over_vacuum": [float(x / (1.0 - P.delta)) for x in an],
            "slope_of_a_through_origin_line_fit_over_vacuum": k_lin / (1.0 - P.delta),
            "rms_residual_of_the_line_over_vacuum": float(
                np.sqrt(np.mean((an - k_lin * rs) ** 2)) / (1.0 - P.delta)
            ),
            "log_log_slope_4.5_to_18": float(
                np.polyfit(np.log(rs[2:12]), np.log(np.maximum(an[2:12], 1e-300)), 1)[0]
            ),
        }
        out8[nm] = {
            "anisotropy_ramp": ramp,
            "label": d["row"]["label"],
            "iters": d["row"]["iters"],
            "E": tot,
            "E_seed": float((ecs + sum(pot_density(seed, np.full(P.r.shape, G8), P))).sum()),
            "E_curv": float(d["ec"].sum()),
            "V4": float(d["v4"].sum()),
            "V4_of_a_fully_isotropic_free_region": float(
                pot_density(
                    np.broadcast_to((1 + 2 * P.delta) / 3 * np.eye(3), d["S"].shape),
                    solve_u(
                        np.broadcast_to((1 + 2 * P.delta) / 3 * np.eye(3), d["S"].shape),
                        P,
                        np.full(P.r.shape, G8),
                    ),
                    P,
                )[0][P.free].sum()
            ),
            "r_half": rr_[0],
            "r_half_interpolated": rr_[1],
            "fmax_mine": gr[t]["fmax_mine_record_convention"],
            "gate": P.gate,
            "shells": sh,
            "energy_fraction_by_layer": {
                k: float(d["e"][m].sum() / tot) for k, m in layers.items()
            },
            "volume_fraction_by_layer": {k: float(m.mean()) for k, m in layers.items()},
            "top_eigvec_alignment_with_rhat_mean_free": float(align[P.free].mean()),
            "min_gap_top_mid_free": float((lam[..., 2] - lam[..., 1])[P.free].min()),
            "min_gap_mid_small_free": float((lam[..., 1] - lam[..., 0])[P.free].min()),
            "degree": {k: deg[t][k] for k in ("6", "9", "12")},
            "volume_disclinations": deg[t]["volume"],
        }
    a0 = out8["A_c0"]
    sp = [a0["shells"][f"{R:g}"]["spectrum"] for R in SHELLS if R <= 15.0]
    okA0 = (
        all(abs(x[0] - 0.905) < 0.01 and abs(x[1] - 0.949) < 0.01 and abs(x[2] - 0.992) < 0.01
            for x in sp)
        and abs(a0["r_half"] - 19.9) < 0.1
        and all(a0["degree"][k]["frustrated_plaquettes"] > 0 for k in ("6", "9", "12"))
    )  # fmt: skip
    add(
        "C8.rowA_c0",
        "row A at c = 0: spectrum about (0.905, 0.949, 0.992) flat in r from 1.5 to 15, r_half "
        "19.9, the top eigenvector not orientable on all three cubes",
        "own shells (0.75 h), own r_half, own degree reader",
        {k: a0[k] for k in ("shells", "r_half", "degree", "label", "fmax_mine", "gate")},
        "every shell mean within 0.01 of the quoted triple, r_half within 0.1, frustrated "
        "squares on all three cubes",
        "a shell off by more than 0.01, another radius, or an orientable cube",
        okA0,
    )
    ac = out8["A_c"]
    sp = [ac["shells"][f"{R:g}"]["spectrum"] for R in SHELLS if R <= 15.0]
    okAc = (
        all(abs(x[0] - 0.910) < 0.01 and abs(x[1] - 0.940) < 0.01 and abs(x[2] - 0.995) < 0.01
            for x in sp)
        and abs(ac["r_half"] - 15.2) < 0.1
        and ac["degree"]["6"]["frustrated_plaquettes"] > 0
        and ac["degree"]["9"]["frustrated_plaquettes"] == 0
        and ac["degree"]["12"]["frustrated_plaquettes"] == 0
    )  # fmt: skip
    add(
        "C8.rowA_c3.094e-5",
        "row A at c = 3.094e-5: spectrum about (0.910, 0.940, 0.995) flat in r, r_half 15.2, "
        "unit degree on r 9 and r 12, conflicts on r 6",
        "own shells (0.75 h), own r_half, own degree reader",
        {k: ac[k] for k in ("shells", "r_half", "degree", "label", "fmax_mine", "gate")},
        "every shell mean r 1.5 to 15 within 0.01 of the quoted triple, r_half within 0.1, the "
        "claimed orientability pattern",
        "a shell off by more than 0.01 (the spectrum is then not flat), another radius, or "
        "another pattern",
        okAc,
    )
    b0 = out8["B_c0"]
    sB = b0["shells"]
    okB = (
        max(abs(sB["1.5"]["spectrum"][a] - (0.946, 0.948, 0.954)[a]) for a in range(3)) < 0.003
        and max(abs(sB["15"]["spectrum"][a] - (0.930, 0.937, 0.980)[a]) for a in range(3)) < 0.003
        and b0["fmax_mine"] < b0["gate"]
        and abs(b0["E"] - 1.545e-4) < 1e-7
        and abs(b0["E_seed"] - 8.66e-4) < 1e-6
    )
    add(
        "C8.rowB_numbers",
        "row B (w1s 0.003779, c = 0): fmax below the gate, E 1.545e-4 from 8.66e-4, spectrum "
        "about (0.946, 0.948, 0.954) at r 1.5 to (0.930, 0.937, 0.980) at r 15",
        "own energy, own gradient over all free cells, own shells, own seed energy",
        {k: b0[k] for k in ("E", "E_seed", "fmax_mine", "gate", "label", "iters", "shells")},
        "shell means within 0.003, energies as quoted, my fmax under the gate",
        "any of those numbers off",
        okB,
    )
    g_in = sB["1.5"]["G1sq_over_vacuum"]
    g_out = sB["15"]["G1sq_over_vacuum"]
    add(
        "C8.rowB_gap_closing",
        "row B lowers the curvature energy by closing the eigenvalue gaps (the quartic is "
        "weighted by products of gaps): 'the potential at this weight does not hold the "
        "spectrum'",
        "shell mean of G_1^2 = ((top - mid)(top - small))^2 over its vacuum value (1 - delta)^4, "
        "the curvature density over the seed's, the V4 paid against the curvature saved, and "
        "the V4 cost of a fully isotropic free region",
        {
            "G1sq_over_vacuum_by_shell": {k: v["G1sq_over_vacuum"] for k, v in sB.items()},
            "seed_G1sq_over_vacuum_by_shell": {
                k: v["seed_G1sq_over_vacuum"] for k, v in sB.items()
            },
            "curv_density_over_seed_by_shell": {
                k: v["curv_density_over_seed"] for k, v in sB.items()
            },
            "anisotropy_ramp_top_minus_mean_of_pair": b0["anisotropy_ramp"],
            "V4_paid": b0["V4"],
            "curvature_saved": b0["E_seed"] - b0["E_curv"],
            "V4_fully_isotropic": b0["V4_of_a_fully_isotropic_free_region"],
            "energy_fraction_by_layer": b0["energy_fraction_by_layer"],
            "volume_fraction_by_layer": b0["volume_fraction_by_layer"],
            "degree": {k: v["degree_bfs"] for k, v in b0["degree"].items()},
        },
        "G_1^2 under 10 percent of vacuum in the interior (r 1.5), rising toward the pin, and "
        "V4 paid under 10 percent of the curvature saved",
        "the gaps are not closed in the interior, or the V4 paid is comparable to the saving",
        g_in < 0.1 and g_out > g_in and b0["V4"] < 0.1 * (b0["E_seed"] - b0["E_curv"]),
    )
    ga = a0["shells"]
    tear = {
        "G1sq_over_vacuum_by_shell": {k: v["G1sq_over_vacuum"] for k, v in ga.items()},
        "min_gap_top_mid_free": a0["min_gap_top_mid_free"],
        "vacuum_gap": 1.0 - D[T_A0]["P"].delta,
        "volume_disclinations": a0["volume_disclinations"],
        "frustrated_by_cube": {
            k: {
                "count": v["frustrated_plaquettes"],
                "gap_at_frustrated": v["gap_top_mid_at_frustrated"],
                "gap_surface_median": v["gap_top_mid_surface_median"],
            }
            for k, v in a0["degree"].items()
        },
        "top_alignment_with_rhat": a0["top_eigvec_alignment_with_rhat_mean_free"],
        "rowB_top_alignment_with_rhat": b0["top_eigvec_alignment_with_rhat_mean_free"],
    }
    n_fr = [v["count"] for v in tear["frustrated_by_cube"].values()]
    gap_ratio = max(
        (v["gap_at_frustrated"][0] / v["gap_surface_median"])
        for v in tear["frustrated_by_cube"].values()
        if v["gap_at_frustrated"]
    )
    add(
        "C8.rowA_tearing",
        "row A at c = 0 lost the hedgehog by tearing (disclination lines of the top eigenvector "
        "piercing the cubes where top and mid cross locally), not by a uniform gap closing",
        "gauge-invariant frustrated squares on each cube and in the volume, and the top - mid "
        "gap at the frustrated squares against the surface median",
        tear,
        "a small even number of piercings per cube with the gap at the piercings under 30 "
        "percent of the surface median, while G_1^2 stays above 10 percent of vacuum",
        "no localized gap collapse at the piercings (then it is not tearing), or G_1^2 under 10 "
        "percent everywhere (then it is the row B melting)",
        all(x > 0 and x % 2 == 0 for x in n_fr)
        and gap_ratio < 0.3
        and min(tear["G1sq_over_vacuum_by_shell"].values()) > 0.1,
    )

    # ---------------- C9: the labels ----------------
    ref = D[T_RAD[0.0]]["last"]["r_half"]
    lab = {
        "rowB_r_half": b0["r_half"],
        "rowA_c0_r_half": a0["r_half"],
        "rowA_c_r_half": ac["r_half"],
        "delta03_c0_r_half": ref,
        "rowB_over_reference": b0["r_half"] / ref,
        "TWIN_SAME_OBJECT": bool(abs(b0["r_half"] / ref - 1.0) < 0.10),
        "DELTA_COMPACT_any": bool(min(b0["r_half"], a0["r_half"], ac["r_half"]) < 8.0),
        "dip_confined_to_2_cells_rowA": bool(
            abs(a0["shells"]["6"]["spectrum"][1] - D[T_A0]["P"].delta) < 0.005
        ),
        "ceiling_uniform_cube": uni[0],
        "ceiling_uniform_free_region": uni_free[0],
        "rowB_minus_free_ceiling": b0["r_half"] - uni_free[0],
        "rowA_c0_minus_free_ceiling": a0["r_half"] - uni_free[0],
    }
    add(
        "C9.labels",
        "DELTA_BOX_FILLING for both rows: not TWIN_SAME_OBJECT (row B radius not within 10 "
        "percent of 17.02), not CORE_LATTICE_PINNED (the eigenvalue dip is not confined to 2 "
        "cells), not DELTA_COMPACT (no radius under 8)",
        "own radii and own shell spectra against the pre-registered label definitions",
        lab,
        "all three alternative labels fail on my numbers",
        "any alternative label is met",
        (not lab["TWIN_SAME_OBJECT"])
        and (not lab["DELTA_COMPACT_any"])
        and (not lab["dip_confined_to_2_cells_rowA"]),
    )
    add(
        "C9.ceiling_distinguishable",
        "WORDING: an r_half of about 20 is a measured radius (distinguishable from the reader's "
        "uniform-density ceiling)",
        "row A (c = 0) and row B r_half against the r_half of a uniform density on the free "
        "region (20.7) and on the whole cube (23.6)",
        lab,
        "both radii more than 10 percent below the free-region ceiling",
        "a radius within 10 percent of the uniform free-region reading (it then carries no "
        "size information beyond 'the box')",
        b0["r_half"] < 0.9 * uni_free[0] and a0["r_half"] < 0.9 * uni_free[0],
    )

    n_pass = sum(1 for c in CHECKS if c["verdict"] == "PASS")
    out = {
        "task": "M5.32 R23-1 / R23-2 adversarial audit",
        "independence": "own energy, gradient, M_00 solve, radius readers, shells, degree "
        "reader and symmetry images; no import from any script in this folder",
        "checks": CHECKS,
        "counts": {"PASS": n_pass, "FAIL": len(CHECKS) - n_pass, "total": len(CHECKS)},
        "wall_s": round(time.time() - T0, 1),
    }
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=1, default=float)
    log(f"PASS {n_pass} FAIL {len(CHECKS) - n_pass} -> {os.path.relpath(OUT_JSON, HERE)}")
    for c in CHECKS:
        if c["verdict"] == "FAIL":
            log(f"FAIL {c['id']}")


if __name__ == "__main__":
    main()
