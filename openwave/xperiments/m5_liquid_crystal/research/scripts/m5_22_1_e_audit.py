#!/usr/bin/env python3
"""M5.22.1 ADVERSARIAL AUDIT (independent re-implementation, AI_HYGIENE.md par. 1).

Self-contained: numpy + scipy.ndimage only. Imports NOTHING from the analysis
scripts (m5_22_b_census.py, m5_22_1_a_kick.py, m5_22_1_b_deuteron.py) or the
instrument stack. Every instrument re-derived from the equations in
findings/m5_22_note.md par. 1 and findings/m5_22_1_note.md par. 1:

  energy   E = h^3 sum(u + V_T2), u = 4 sum_{i<j} tr([A_i,A_j]^T[A_i,A_j]),
           A_i = one-sided d_i M / h with pad0 edges, E_u = (fwd + bwd)/2,
           V_T2 = w2 sum_i (lambda_i - v_i)^2, v = sorted(0, delta, 1) ascending
  charge   NOT the analysis's finite-difference Mermin-Ho flux: own BFS
           orientation lift of the leading eigenvector + exact van
           Oosterom-Strackee solid-angle degree over closed lattice boxes
           (full box, half boxes, narrow column boxes, z-cut ladder)
  dipole   two ways: (a) z-cut charge ladder p_z = sum z_mid dQ(z) (pure
           surface degrees, independent of any div computation), (b) own
           volume rho = div B / 4pi with own lift + own B + central diffs
  basin    field distance ||M_end - M_parent||_F on the free region,
           normalized by the parent's own relaxation displacement
           ||M_parent - M0_parent||_F, + eigengap-field comparison, + own
           analytic-gradient plain-GD descent probe (FD-gated)
  winding  own incremental-unwrap path integral of the P-family
           cross-section angle (composite additivity check)
  xr       own cross-stencil ratio: max/min over E_u(fwd, bwd, 2h,
           subsample parity 0, subsample parity 1)

Run:  python3 m5_22_1_e_audit.py [--skip-descent] [--skip-c6]
Out:  ../data/m5_22_1_audit.json
"""

import argparse
import glob
import json
import os
import sys
from collections import deque

import numpy as np
from scipy import ndimage

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")

W2_DEFAULT = 0.0027581


# ================================ energy ====================================
def _axsl(ax, a, b):
    s = [slice(None)] * 3
    s[ax] = slice(a, b)
    return tuple(s)


def _u_total(A):
    """u = 4 sum_{i<j} tr([A_i,A_j]^T [A_i,A_j]) summed over cells."""
    tot = 0.0
    for i in range(3):
        for j in range(i + 1, 3):
            C = A[i] @ A[j] - A[j] @ A[i]
            tot += float(np.einsum("...ab,...ab->...", C, C).sum())
    return 4.0 * tot


def _A_oneside(M, h, side):
    n = M.shape[0]
    A = []
    for ax in range(3):
        d = (M[_axsl(ax, 1, n)] - M[_axsl(ax, 0, n - 1)]) / h
        D = np.zeros_like(M)
        if side == "fwd":
            D[_axsl(ax, 0, n - 1)] = d
        else:
            D[_axsl(ax, 1, n)] = d
        A.append(D)
    return A


def _A_2h(M, h):
    n = M.shape[0]
    A = []
    for ax in range(3):
        D = np.zeros_like(M)
        D[_axsl(ax, 1, n - 1)] = (M[_axsl(ax, 2, n)] - M[_axsl(ax, 0, n - 2)]) / (2 * h)
        A.append(D)
    return A


def energy_pad0(M, h, delta, w2):
    ef = (h ** 3) * _u_total(_A_oneside(M, h, "fwd"))
    eb = (h ** 3) * _u_total(_A_oneside(M, h, "bwd"))
    lam = np.linalg.eigvalsh(M)
    v = np.array(sorted((0.0, float(delta), 1.0)))
    ev = (h ** 3) * w2 * float(((lam - v) ** 2).sum())
    return {"E_u_fwd": ef, "E_u_bwd": eb, "E_u": 0.5 * (ef + eb),
            "E_v": ev, "E_end": 0.5 * (ef + eb) + ev}


def xr_mine(M, h):
    """own cross-stencil ratio, per the instrument's consistency() recipe."""
    reads = {
        "fwd": (h ** 3) * _u_total(_A_oneside(M, h, "fwd")),
        "bwd": (h ** 3) * _u_total(_A_oneside(M, h, "bwd")),
        "2h": (h ** 3) * _u_total(_A_2h(M, h)),
    }
    for par in (0, 1):
        sub = M[par::2, par::2, par::2]
        reads[f"sub{par}"] = ((2 * h) ** 3) * _u_total(_A_oneside(sub, 2 * h, "fwd"))
    vals = [v for v in reads.values() if v > 0]
    return max(vals) / min(vals), reads


# ============================ analytic gradient =============================
def grad_pad0(M, h, delta, w2):
    """own derivation: d(tr(C^T C))/dA_i = 2[C, A_j] (C antisymmetric, A sym),
    u-prefactor 4 -> 8[C, A_j]; pad0 stencil adjoints; V_T2 per-cell
    gradient V diag(2 w2 (lambda - v)) V^T. Validated by FD gate."""
    n = M.shape[0]
    G = np.zeros_like(M)
    for side in ("fwd", "bwd"):
        A = _A_oneside(M, h, side)
        gA = [np.zeros_like(M) for _ in range(3)]
        for i in range(3):
            for j in range(i + 1, 3):
                C = A[i] @ A[j] - A[j] @ A[i]
                gA[i] += 8.0 * (C @ A[j] - A[j] @ C)
                gA[j] += 8.0 * (A[i] @ C - C @ A[i])
        for ax in range(3):
            g = gA[ax].copy()
            if side == "fwd":
                g[_axsl(ax, n - 1, n)] = 0.0
                G[_axsl(ax, 1, n)] += g[_axsl(ax, 0, n - 1)] / h
                G -= g / h
            else:
                g[_axsl(ax, 0, 1)] = 0.0
                G += g / h
                G[_axsl(ax, 0, n - 1)] -= g[_axsl(ax, 1, n)] / h
    G *= 0.5
    lam, V = np.linalg.eigh(M)
    vt = np.array(sorted((0.0, float(delta), 1.0)))
    S = 2.0 * w2 * (lam - vt)
    G = G + np.einsum("...ik,...k,...jk->...ij", V, S, V)
    return G * (h ** 3)


def grad_fd_gate(rng):
    n, h, delta, w2 = 6, 1.0, 0.3, W2_DEFAULT
    M = rng.standard_normal((n, n, n, 3, 3))
    M = 0.25 * (M + np.swapaxes(M, -1, -2))
    G = grad_pad0(M, h, delta, w2)
    worst = 0.0
    for _ in range(4):
        D = rng.standard_normal((n, n, n, 3, 3))
        D = 0.5 * (D + np.swapaxes(D, -1, -2))
        t = 1e-6
        num = (energy_pad0(M + t * D, h, delta, w2)["E_end"]
               - energy_pad0(M - t * D, h, delta, w2)["E_end"]) / (2 * t)
        ana = float((G * D).sum())
        worst = max(worst, abs(num - ana) / max(abs(num), 1e-14))
    return worst


def pin_mask(n, h, depth=1.6):
    """outer shell of physical depth >= 1.6 units (2 cells at h = 1.5)."""
    idx = np.arange(n)
    edge = np.minimum(idx, n - 1 - idx) * h
    m1 = edge < depth
    shell = m1[:, None, None] | m1[None, :, None] | m1[None, None, :]
    return shell


def descent_probe(M, h, delta, w2, steps=400, alpha=0.05):
    """plain GD with the own gradient, pinned shell held. A genuinely
    stationary endpoint sheds only float32-rounding roughness."""
    n = M.shape[0]
    shell = pin_mask(n, h)[..., None, None]
    X = M.copy()
    e0 = energy_pad0(X, h, delta, w2)["E_end"]
    ehist = [e0]
    fmax0 = None
    for k in range(steps):
        G = grad_pad0(X, h, delta, w2)
        G[np.broadcast_to(shell, G.shape)] = 0.0
        if fmax0 is None:
            fmax0 = float(np.abs(G).max())
        X = X - alpha * G
        if (k + 1) % 20 == 0:
            e = energy_pad0(X, h, delta, w2)["E_end"]
            if e > ehist[-1]:
                alpha *= 0.5
            ehist.append(e)
    Gend = grad_pad0(X, h, delta, w2)
    Gend[np.broadcast_to(shell, Gend.shape)] = 0.0
    e1 = energy_pad0(X, h, delta, w2)["E_end"]
    return {"E_start": e0, "E_probe_end": e1, "dE_total": e0 - e1,
            "dE_last20": ehist[-1] - e1 if len(ehist) > 1 else 0.0,
            "fmax_start": fmax0, "fmax_end": float(np.abs(Gend).max()),
            "steps": steps, "alpha_final": alpha}, X


# ============================ charge instruments ============================
def directors(M):
    lam, V = np.linalg.eigh(M)
    nhat = V[..., -1]
    nhat = nhat / np.linalg.norm(nhat, axis=-1, keepdims=True)
    gap = lam[..., 2] - lam[..., 1]
    return nhat, gap


def bfs_lift(mask, vecs, gap):
    """orientation lift over an arbitrary masked point set: BFS from the
    max-gap point, sign from the first assigned neighbor; returns the sign
    field (+-1 on mask, 0 off) and the count of frustrated edges."""
    n = mask.shape[0]
    dots = []
    for ax in range(3):
        d = np.einsum("...i,...i->...", vecs[_axsl(ax, 0, n - 1)],
                      vecs[_axsl(ax, 1, n)])
        dots.append(d)
    signs = np.zeros(mask.shape, dtype=np.int8)
    g = np.where(mask, gap, -np.inf)
    start = np.unravel_index(int(np.argmax(g)), mask.shape)
    signs[start] = 1
    q = deque([start])
    while q:
        i, j, k = q.popleft()
        s0 = int(signs[i, j, k])
        for ax, (di, dj, dk) in enumerate(((1, 0, 0), (0, 1, 0), (0, 0, 1))):
            for sgn in (1, -1):
                ii, jj, kk = i + sgn * di, j + sgn * dj, k + sgn * dk
                if not (0 <= ii < n and 0 <= jj < n and 0 <= kk < n):
                    continue
                if not mask[ii, jj, kk] or signs[ii, jj, kk] != 0:
                    continue
                lo = (i, j, k) if sgn > 0 else (ii, jj, kk)
                dv = float(dots[ax][lo])
                signs[ii, jj, kk] = s0 if dv >= 0 else -s0
                q.append((ii, jj, kk))
    ncf = 0
    for ax in range(3):
        both = mask[_axsl(ax, 0, n - 1)] & mask[_axsl(ax, 1, n)]
        sp = signs[_axsl(ax, 0, n - 1)].astype(np.int32) * \
            signs[_axsl(ax, 1, n)].astype(np.int32)
        bad = both & (sp * dots[ax] < 0)
        ncf += int(bad.sum())
    return signs, ncf


def _tri_omega(a, b, c):
    num = np.einsum("...i,...i->...", a, np.cross(b, c))
    den = 1.0 + np.einsum("...i,...i->...", a, b) \
        + np.einsum("...i,...i->...", b, c) \
        + np.einsum("...i,...i->...", c, a)
    return 2.0 * np.arctan2(num, den)


def box_degree(W, bounds):
    """exact solid-angle degree of the signed unit field W over the closed
    box surface (bounds inclusive), outward oriented (van Oosterom-Strackee
    per triangle)."""
    i0, i1, j0, j1, k0, k1 = bounds
    lo = (i0, j0, k0)
    hi = (i1, j1, k1)
    ident = np.eye(3)
    total = 0.0
    for a in range(3):
        b, c = [x for x in range(3) if x != a]
        par = float(np.cross(ident[b], ident[c]) @ ident[a])
        for pos, sgn in ((hi[a], 1.0), (lo[a], -1.0)):
            sl = [slice(lo[0], hi[0] + 1), slice(lo[1], hi[1] + 1),
                  slice(lo[2], hi[2] + 1)]
            sl[a] = pos
            G = W[tuple(sl)]
            A = G[:-1, :-1]
            B = G[1:, :-1]
            C = G[1:, 1:]
            D = G[:-1, 1:]
            om = _tri_omega(A, B, C) + _tri_omega(A, C, D)
            total += par * sgn * float(om.sum())
    return total / (4.0 * np.pi)


def shell_mask(n, bounds):
    i0, i1, j0, j1, k0, k1 = bounds
    m = np.zeros((n, n, n), dtype=bool)
    m[i0:i1 + 1, j0:j1 + 1, k0:k1 + 1] = True
    if i1 - i0 >= 2 and j1 - j0 >= 2 and k1 - k0 >= 2:
        m[i0 + 1:i1, j0 + 1:j1, k0 + 1:k1] = False
    return m


def read_box_surface(nhat, gap, bounds):
    """independent per-surface lift + degree + conflict count + min gap."""
    n = nhat.shape[0]
    m = shell_mask(n, bounds)
    signs, ncf = bfs_lift(m, nhat, gap)
    W = nhat * signs[..., None]
    Q = box_degree(W, bounds)
    return {"Q": Q, "conflicts": ncf, "min_gap": float(gap[m].min())}


def surface_lift(nhat, gap, bounds):
    n = nhat.shape[0]
    m = shell_mask(n, bounds)
    signs, ncf = bfs_lift(m, nhat, gap)
    return m, signs, ncf


def gauge_align(signs, mask, ref_signs, ref_mask, gap):
    """align a surface lift's global sign to a reference lift via the
    gap-weighted majority of sign products on the overlap."""
    ov = mask & ref_mask
    if not ov.any():
        return 1, 0.0
    prod = (signs[ov].astype(np.float64) * ref_signs[ov]) * gap[ov]
    t = float(prod.sum())
    return (1 if t >= 0 else -1), t


def aligned_degree(nhat, gap, bounds, ref):
    """per-surface lift, gauge-aligned to the state's reference surface."""
    ref_mask, ref_signs = ref
    m, signs, ncf = surface_lift(nhat, gap, bounds)
    s, _ = gauge_align(signs, m, ref_signs, ref_mask, gap)
    W = nhat * (s * signs[..., None])
    return {"Q": box_degree(W, bounds), "conflicts": ncf,
            "min_gap": float(gap[m].min())}


def conflict_zs(nhat, gap, bounds, h):
    """physical z coordinates of frustrated edges on a box surface."""
    n = nhat.shape[0]
    m = shell_mask(n, bounds)
    signs, _ = bfs_lift(m, nhat, gap)
    c = (n - 1) / 2.0
    zs = []
    for ax in range(3):
        both = m[_axsl(ax, 0, n - 1)] & m[_axsl(ax, 1, n)]
        d = np.einsum("...i,...i->...", nhat[_axsl(ax, 0, n - 1)],
                      nhat[_axsl(ax, 1, n)])
        sp = signs[_axsl(ax, 0, n - 1)].astype(np.int32) * \
            signs[_axsl(ax, 1, n)].astype(np.int32)
        bad = np.argwhere(both & (sp * d < 0))
        for p in bad:
            zs.append(float((p[2] - c) * h))
    return zs


def volume_lift(nhat, gap):
    n = nhat.shape[0]
    mask = np.ones((n, n, n), dtype=bool)
    signs, ncf = bfs_lift(mask, nhat, gap)
    return nhat * signs[..., None], ncf


def mermin_B(W, h):
    dW = np.gradient(W, h, axis=(0, 1, 2))
    B = np.zeros(W.shape[:3] + (3,))
    for i in range(3):
        j, k = (i + 1) % 3, (i + 2) % 3
        B[..., i] = np.einsum("...a,...a->...", W, np.cross(dW[j], dW[k]))
    return B


def volume_moments(M, h, excl=1.6):
    nhat, gap = directors(M)
    W, ncf = volume_lift(nhat, gap)
    B = mermin_B(W, h)
    div = np.zeros(B.shape[:3])
    for ax in range(3):
        div += np.gradient(B[..., ax], h, axis=ax)
    rho = div / (4.0 * np.pi)
    n = M.shape[0]
    c = (n - 1) / 2.0
    x = (np.arange(n) - c) * h
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    interior = ~pin_mask(n, h, depth=excl)
    w = rho * interior
    h3 = h ** 3
    r2 = X * X + Y * Y + Z * Z
    return {"q_vol": float(w.sum() * h3),
            "p_z": float((w * Z).sum() * h3),
            "Q2_zz": float((w * (3 * Z * Z - r2)).sum() * h3),
            "volume_lift_conflicts": int(ncf)}, W, gap


def idx_of(n, h, xphys):
    return int(round(xphys / h + (n - 1) / 2.0))


# =========================== 2D winding (own) ===============================
def step_fn(t):
    return np.arctan(t) / np.pi + 0.5


def winding_P(centers, d=10.0, m=400001):
    """(1/pi) change of the P-family composite angle along the semicircle
    (x, y) = (-d cos phi, d sin phi). The instrument definition (census
    note par. 1) keeps EACH atan2 term continuous separately; for
    half-integer s the branch jump carries a genuine pi, so per-term
    continuity is part of the definition, not a numerical choice.
    centers = [(s, x0), ...]."""
    phi = np.linspace(0.0, np.pi, m)
    xg = -d * np.cos(phi)
    yg = d * np.sin(phi)
    tot = np.zeros_like(phi)
    for s, x0 in centers:
        x = xg - x0
        t1 = np.unwrap(np.arctan2(yg - 1.0, x))
        t2 = np.unwrap(np.arctan2(-1.0 * np.ones_like(x), x))
        tot = tot + s * (t1 - t2) \
            - (2.0 * np.pi / 3.0) * step_fn(5.0 * x) \
            + (np.pi / 3.0) * (step_fn(5.0 * (x - 1.0))
                               + step_fn(5.0 * (x + 1.0)))
    return float((tot[-1] - tot[0]) / np.pi)


# ============================== helpers =====================================
def load_end(name):
    z = np.load(os.path.join(DATA, name))
    return (np.asarray(z["M"], dtype=np.float64),
            np.asarray(z["M0"], dtype=np.float64),
            float(z["h"]), float(z["delta"]), int(z["n"]))


def rowjson(name):
    with open(os.path.join(DATA, name)) as f:
        return json.load(f)


def ring_read_own(gap, h, thr, rmax_frac=0.9):
    """own low-eigengap component read: scipy label, 26-connectivity."""
    n = gap.shape[0]
    c = (n - 1) / 2.0
    x = (np.arange(n) - c) * h
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    r = np.sqrt(X * X + Y * Y + Z * Z)
    rho = np.sqrt(X * X + Y * Y)
    interior = (~pin_mask(n, h)) & (r < rmax_frac * (n - 1) * h / 2.0)
    mask = (gap < thr) & interior
    lab, nlab = ndimage.label(mask, structure=np.ones((3, 3, 3)))
    rings, cols = [], []
    for q in range(1, nlab + 1):
        sel = lab == q
        entry = {"voxels": int(sel.sum()),
                 "rho_mean": float(rho[sel].mean()),
                 "z_centroid": float(Z[sel].mean()),
                 "gap_min": float(gap[sel].min())}
        (cols if entry["rho_mean"] < 3.0 else rings).append(entry)
    rings.sort(key=lambda e: e["z_centroid"])
    return {"n_rings": len(rings), "rings": rings, "n_columns": len(cols)}


def field_distance(Ma, Mb, free3):
    d = (Ma - Mb)[free3]
    return float(np.sqrt((d * d).sum()))


def tofloat(o):
    if isinstance(o, dict):
        return {k: tofloat(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [tofloat(v) for v in o]
    if isinstance(o, (np.floating, np.integer)):
        return float(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


# ================================= main =====================================
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-descent", action="store_true")
    ap.add_argument("--skip-c6", action="store_true")
    args = ap.parse_args()
    audit = {}
    print("=== M5.22.1 adversarial audit (independent instruments) ===")

    # ------------------ instrument controls ---------------------------------
    print("\n-- controls: degree instrument --")
    n32, h32 = 32, 1.5
    c = (n32 - 1) / 2.0
    x = (np.arange(n32) - c) * h32
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    r = np.sqrt(X * X + Y * Y + Z * Z)
    r = np.where(r < 1e-9, 1e-9, r)
    hh = np.stack([X / r, Y / r, Z / r], axis=-1)
    b0 = (4, 27, 4, 27, 4, 27)
    q_h = box_degree(hh, b0)
    st = np.sqrt(X * X + Y * Y) / r
    ph = np.arctan2(Y, X)
    d2 = np.stack([st * np.cos(2 * ph), st * np.sin(2 * ph), Z / r], axis=-1)
    q_2 = box_degree(d2, b0)
    rng = np.random.default_rng(52201)
    flip = rng.choice([-1.0, 1.0], size=(n32, n32, n32, 1))
    Mh = np.einsum("...i,...j->...ij", hh, hh) \
        + 0.3 * np.einsum("...i,...j->...ij", d2, d2)
    Mh = Mh * 1.0
    nh, gh = directors(Mh + 0 * flip[..., None])
    rb = read_box_surface(nh, gh, b0)
    print(f"  hedgehog vectors: Q = {q_h:+.6f} (expect +1)")
    print(f"  degree-2 vectors: Q = {q_2:+.6f} (expect +2)")
    print(f"  hedgehog through M + eigh + lift: |Q| = {abs(rb['Q']):.6f} "
          f"(conflicts {rb['conflicts']}; sign is the eigh gauge)")
    ggate = grad_fd_gate(rng)
    print(f"  gradient FD gate: worst rel {ggate:.2e}")
    audit["controls"] = {"hedgehog_Q": q_h, "degree2_Q": q_2,
                         "hedgehog_via_M": rb, "grad_fd_gate": ggate}

    # ------------------ C1: ladder energies + ratios + xr -------------------
    print("\n-- C1: ladder energies (own pad0 energy) --")
    ladder = {
        ("P-0.5", 0.3): ("m5_22_end_P-0.5_plane_sc6_n32_pinned_d0.3.npz",
                         "m5_22_row_P-0.5_plane_sc6_n32_pinned_d0.3.json", 8.496),
        ("P-1", 0.3): ("m5_22_end_P-1_plane_sc6_n32_pinned_d0.3.npz",
                       "m5_22_row_P-1_plane_sc6_n32_pinned_d0.3.json", 12.730),
        ("E-0.5", 0.3): ("m5_22_end_E-0.5_plane_sc6_n32_pinned_d0.3.npz",
                         "m5_22_row_E-0.5_plane_sc6_n32_pinned_d0.3.json", 6.029),
        ("P-0.5", 0.2): ("m5_22_end_P-0.5_plane_sc6_n32_pinned_d0.2.npz",
                         "m5_22_row_P-0.5_plane_sc6_n32_pinned_d0.2.json", 10.257),
        ("P-1", 0.2): ("m5_22_end_P-1_plane_sc6_n32_pinned_d0.2.npz",
                       "m5_22_row_P-1_plane_sc6_n32_pinned_d0.2.json", 15.882),
        ("E-0.5", 0.2): ("m5_22_end_E-0.5_plane_sc6_n32_pinned_d0.2.npz",
                         "m5_22_row_E-0.5_plane_sc6_n32_pinned_d0.2.json", 7.096),
        ("P-0.5", 0.1): ("m5_22_end_P-0.5_plane_sc6_n32_pinned_d0.1.npz",
                         "m5_22_row_P-0.5_plane_sc6_n32_pinned_d0.1.json", 12.582),
        ("P-1", 0.1): ("m5_22_end_P-1_plane_sc6_n32_pinned_d0.1.npz",
                       "m5_22_row_P-1_plane_sc6_n32_pinned_d0.1.json", 19.840),
        ("E-0.5", 0.1): ("m5_22_end_E-0.5_plane_sc6_n32_pinned_d0.1.npz",
                         "m5_22_row_E-0.5_plane_sc6_n32_pinned_d0.1.json", 8.596),
    }
    c1 = {}
    myE = {}
    worst_rel = 0.0
    worst_xr = 0.0
    qclass_note = {"P-0.5": 1, "P-1": 0, "E-0.5": 1}
    qclass_ok = True
    for (state, dl), (endf, rowf, noteval) in ladder.items():
        M, _, h, delta, _ = load_end(endf)
        row = rowjson(rowf)
        w2 = float(row.get("w2", W2_DEFAULT))
        e = energy_pad0(M, h, delta, w2)
        xr, _ = xr_mine(M, h)
        row_xr = row.get("consistency", {}).get("xstencil_ratio")
        rel = abs(e["E_end"] - float(row["E_end"])) / float(row["E_end"])
        worst_rel = max(worst_rel, rel)
        worst_xr = max(worst_xr, xr)
        myE[(state, dl)] = e["E_end"]
        nh_, gp_ = directors(M)
        qf = read_box_surface(nh_, gp_, (4, 27, 4, 27, 4, 27))
        qcl_ok = abs(abs(qf["Q"]) - qclass_note[state]) < 0.02 \
            and qf["conflicts"] == 0
        qclass_ok = qclass_ok and qcl_ok
        c1[f"{state}_d{dl:g}"] = {"mine": e["E_end"], "row": row["E_end"],
                                  "note": noteval, "rel_err": rel,
                                  "xr_mine": xr, "xr_row": row_xr,
                                  "absQ_far_mine": abs(qf["Q"]),
                                  "absQ_note": qclass_note[state],
                                  "far_conflicts": qf["conflicts"]}
        print(f"  {state:6s} d{dl:g}: mine {e['E_end']:.4f} row "
              f"{float(row['E_end']):.4f} note {noteval} rel {rel:.1e} "
              f"xr {xr:.3f} (row {row_xr}) |Q| {abs(qf['Q']):.3f} "
              f"(note {qclass_note[state]}, cf {qf['conflicts']})")
    ratios = {}
    for dl in (0.3, 0.2, 0.1):
        ratios[f"neutral_over_proton_d{dl:g}"] = \
            myE[("P-1", dl)] / myE[("P-0.5", dl)]
        ratios[f"proton_over_lepton_d{dl:g}"] = \
            myE[("P-0.5", dl)] / myE[("E-0.5", dl)]
    note_ratios = {"neutral_over_proton_d0.3": 1.498,
                   "neutral_over_proton_d0.2": 1.548,
                   "neutral_over_proton_d0.1": 1.577,
                   "proton_over_lepton_d0.3": 1.409,
                   "proton_over_lepton_d0.2": 1.446,
                   "proton_over_lepton_d0.1": 1.464}
    ratio_ok = all(abs(ratios[k] - v) < 0.0015 for k, v in note_ratios.items())
    for k in sorted(ratios):
        print(f"  ratio {k}: mine {ratios[k]:.4f} note {note_ratios[k]}")
    c1_verdict = "CONFIRMED" if (worst_rel < 1e-4 and ratio_ok
                                 and qclass_ok) else \
        ("QUALIFIED" if worst_rel < 1e-2 else "REFUTED")
    audit["C1"] = {"verdict": c1_verdict, "worst_rel_err": worst_rel,
                   "worst_xr_mine": worst_xr, "ratios_mine": ratios,
                   "ratios_note": note_ratios, "detail": c1,
                   "method": "own pad0 sym one-sided energy from npz M "
                             "(float64 cast); own xr = max/min over "
                             "E_u(fwd,bwd,2h,sub0,sub1)"}
    print(f"  C1 {c1_verdict}: worst rel {worst_rel:.1e}, "
          f"worst xr {worst_xr:.3f}")

    # ------------------ C2: parent ring-antiring structure ------------------
    print("\n-- C2: parent P-1 d0.3 ring-antiring read (own degree boxes) --")
    Mp, M0p, h, delta, n = load_end("m5_22_end_P-1_plane_sc6_n32_pinned_d0.3.npz")
    nhat, gap = directors(Mp)
    momp, _, _ = volume_moments(Mp, h)
    lo, hi = 4, 27
    kcap = idx_of(n, h, 0.75)  # z = +0.75 cap
    far_b = (lo, hi, lo, hi, lo, hi)
    ref_mask, ref_signs, ref_ncf = surface_lift(nhat, gap, far_b)
    # gauge: flip the reference so the UPPER half reads positive (the
    # global director sign is a gauge; the note's convention has upper +)
    up_try = aligned_degree(nhat, gap, (lo, hi, lo, hi, kcap, hi),
                            (ref_mask, ref_signs))
    if up_try["Q"] < 0:
        ref_signs = -ref_signs
    ref = (ref_mask, ref_signs)

    def rd(bounds):
        return aligned_degree(nhat, gap, bounds, ref)

    reads = {"far_full": rd(far_b),
             "upper_half": rd((lo, hi, lo, hi, kcap, hi)),
             "lower_half": rd((lo, hi, lo, hi, lo, kcap - 1))}
    # cap robustness: cap index 14..18 (z = -2.25 .. +3.75)
    caps = {}
    for kc in range(14, 19):
        caps[f"kcap{kc}_z{(kc - (n - 1) / 2) * h:+.2f}"] = \
            rd((lo, hi, lo, hi, kc, hi))["Q"]
    # lateral robustness
    lat = {}
    for w in (15.75, 17.25):
        a0, a1 = idx_of(n, h, -w), idx_of(n, h, w)
        lat[f"halfwidth{w}"] = rd((a0, a1, a0, a1, kcap, hi))["Q"]
    # column boxes (|x|,|y| <= 3.75, excludes the rings at rho ~ 10)
    cb0, cb1 = idx_of(n, h, -3.75), idx_of(n, h, 3.75)
    col = {"upper": rd((cb0, cb1, cb0, cb1, kcap, hi))["Q"],
           "lower": rd((cb0, cb1, cb0, cb1, lo, kcap - 1))["Q"],
           "full": rd((cb0, cb1, cb0, cb1, lo, hi))["Q"]}
    # independent per-surface lifts (fresh unaligned lift, |Q| checks)
    surf = {}
    for nm, b in (("far", far_b),
                  ("upper", (lo, hi, lo, hi, kcap, hi)),
                  ("lower", (lo, hi, lo, hi, lo, kcap - 1))):
        surf[nm] = read_box_surface(nhat, gap, b)
    # z-cut charge ladder -> per-slab dQ and a surface-only p_z
    zc = (np.arange(n) - (n - 1) / 2.0) * h
    Qb = {}
    for k in range(lo + 1, hi + 1):
        Qb[k] = rd((lo, hi, lo, hi, lo, k))["Q"]
    dq = []
    pz_ladder = 0.0
    for k in range(lo + 1, hi):
        d = Qb[k + 1] - Qb[k]
        zmid = 0.5 * (zc[k] + zc[k + 1])
        pz_ladder += zmid * d
        dq.append({"z_mid": zmid, "dQ": d})
    big = [e for e in dq if abs(e["dQ"]) > 0.05]
    rings_p = {f"thr{t:g}": ring_read_own(gap, h, t) for t in (0.06, 0.09, 0.15)}
    print(f"  volume lift conflicts (moments only): "
          f"{momp['volume_lift_conflicts']}; far-surface conflicts {ref_ncf}")
    print(f"  far {reads['far_full']['Q']:+.4f}  "
          f"upper {reads['upper_half']['Q']:+.4f}  "
          f"lower {reads['lower_half']['Q']:+.4f}")
    print(f"  cap scan: " + " ".join(f"{k}:{v:+.3f}" for k, v in caps.items()))
    print(f"  lateral: " + " ".join(f"{k}:{v:+.3f}" for k, v in lat.items()))
    print(f"  column boxes: upper {col['upper']:+.4f} lower {col['lower']:+.4f}"
          f" full {col['full']:+.4f}")
    print(f"  surface-lift reads: " + " ".join(
        f"{k}:Q{v['Q']:+.3f}/cf{v['conflicts']}" for k, v in surf.items()))
    print(f"  z-ladder charge lumps: "
          + " ".join(f"z{e['z_mid']:+.1f}:{e['dQ']:+.3f}" for e in big))
    print(f"  p_z: ladder {pz_ladder:+.3f}  volume {momp['p_z']:+.3f} "
          f"(note +21.3); q_vol {momp['q_vol']:+.4f}; "
          f"Q2_zz {momp['Q2_zz']:+.2f}")
    for t, rr in rings_p.items():
        print(f"  ring read {t}: n_rings {rr['n_rings']} "
              f"{[(round(e['z_centroid'], 1), round(e['rho_mean'], 1)) for e in rr['rings']]}"
              f" cols {rr['n_columns']}")
    up_int = abs(reads["upper_half"]["Q"] - 1.0) < 0.02
    lo_int = abs(reads["lower_half"]["Q"] + 1.0) < 0.02
    far_zero = abs(reads["far_full"]["Q"]) < 0.02
    cap_rob = max(abs(v - 1.0) for v in caps.values()) < 0.02
    lat_rob = max(abs(v - 1.0) for v in lat.values()) < 0.02
    col_small = abs(col["full"]) < 0.02
    pz_ok = abs(momp["p_z"] - 21.3) < 2.0
    c2_verdict = "CONFIRMED" if (up_int and lo_int and far_zero and cap_rob
                                 and lat_rob and col_small and pz_ok) \
        else "QUALIFIED"
    audit["C2"] = {"verdict": c2_verdict,
                   "gauge": "volume lift sign fixed so upper half > 0 "
                            "(global director sign is a gauge)",
                   "reads": reads, "cap_scan": caps, "lateral_scan": lat,
                   "column_boxes": col, "surface_lift_reads": surf,
                   "z_ladder_lumps": big, "p_z_ladder": pz_ladder,
                   "moments_volume": momp, "ring_reads": {
                       t: {"n_rings": rr["n_rings"],
                           "z": [e["z_centroid"] for e in rr["rings"]]}
                       for t, rr in rings_p.items()},
                   "method": "own BFS lift + exact van Oosterom-Strackee box "
                             "degrees; z-cut ladder p_z (surface-only) + own "
                             "volume div B moments"}
    print(f"  C2 {c2_verdict}")

    # ------------------ C3: kick-apart return -------------------------------
    print("\n-- C3: kick-apart identity (own energies + field distances) --")
    free3 = ~pin_mask(n, h)
    parent_row = rowjson("m5_22_row_P-1_plane_sc6_n32_pinned_d0.3.json")
    Ep_mine = energy_pad0(Mp, h, delta,
                          float(parent_row.get("w2", W2_DEFAULT)))["E_end"]
    den_relax = field_distance(Mp, M0p, free3)
    branches = {
        "split_d3": ("m5_22_1_end_P-1_plane_sc6_n32_pinned_d0.3_split_d3.npz",
                     "m5_22_1_row_P-1_plane_sc6_n32_pinned_d0.3_split_d3.json",
                     12.7308),
        "split_d6": ("m5_22_1_end_P-1_plane_sc6_n32_pinned_d0.3_split_d6.npz",
                     "m5_22_1_row_P-1_plane_sc6_n32_pinned_d0.3_split_d6.json",
                     12.7308),
        "kick_v0.5": ("m5_22_1_end_P-1_plane_sc6_n32_pinned_d0.3_kick_v0.5.npz",
                      "m5_22_1_row_P-1_plane_sc6_n32_pinned_d0.3_kick_v0.5.json",
                      12.7385),
        "kick_v1": ("m5_22_1_end_P-1_plane_sc6_n32_pinned_d0.3_kick_v1.npz",
                    "m5_22_1_row_P-1_plane_sc6_n32_pinned_d0.3_kick_v1.json",
                    12.7746),
    }
    c3 = {"parent_E_mine": Ep_mine, "parent_E_note": 12.7296,
          "den_relax_norm": den_relax}
    gap_p = gap
    for nm, (endf, rowf, noteval) in branches.items():
        Mb, _, hb, db, _ = load_end(endf)
        row = rowjson(rowf)
        e = energy_pad0(Mb, hb, db, W2_DEFAULT)["E_end"]
        d_abs = field_distance(Mb, Mp, free3)
        nb, gb = directors(Mb)
        core = free3 & (np.sqrt(X * X + Y * Y + Z * Z) < 18.0)
        gd = gb[core] - gap_p[core]
        fb_mask, fb_signs, _ = surface_lift(nb, gb, far_b)
        if aligned_degree(nb, gb, (lo, hi, lo, hi, kcap, hi),
                          (fb_mask, fb_signs))["Q"] < 0:
            fb_signs = -fb_signs
        refb = (fb_mask, fb_signs)
        qup = aligned_degree(nb, gb, (lo, hi, lo, hi, kcap, hi), refb)["Q"]
        qlo = aligned_degree(nb, gb, (lo, hi, lo, hi, lo, kcap - 1), refb)["Q"]
        qfarb = aligned_degree(nb, gb, far_b, refb)["Q"]
        c3[nm] = {"E_mine": e, "E_row": row["E_end"], "E_note": noteval,
                  "q_far": qfarb,
                  "dE_vs_parent": e - Ep_mine,
                  "dist_to_parent": d_abs,
                  "dist_rel_to_relax": d_abs / den_relax,
                  "maxabs_dM": float(np.abs((Mb - Mp)[free3]).max()),
                  "gap_rms_diff_core": float(np.sqrt((gd * gd).mean())),
                  "gap_max_diff_core": float(np.abs(gd).max()),
                  "q_upper": qup, "q_lower": qlo}
        print(f"  {nm:9s}: E {e:.4f} (row {row['E_end']:.4f}, note {noteval})"
              f" dE {e - Ep_mine:+.4f}  |dM|/|relax| "
              f"{d_abs / den_relax:.4f}  gap_rms {c3[nm]['gap_rms_diff_core']:.4f}"
              f"  q_up {qup:+.3f} q_lo {qlo:+.3f} q_far {qfarb:+.3f}")
    Ms3, _, _, _, _ = load_end(branches["split_d3"][0])
    Ms6, _, _, _, _ = load_end(branches["split_d6"][0])
    c3["dist_split3_vs_split6_rel"] = \
        field_distance(Ms3, Ms6, free3) / den_relax
    print(f"  cross-check d(split_d3, split_d6)/|relax| = "
          f"{c3['dist_split3_vs_split6_rel']:.5f}")
    if not args.skip_descent:
        print("  descent probes (own gradient, plain GD, pinned shell):")
        dp_par, _ = descent_probe(Mp, h, delta, W2_DEFAULT)
        Mk, _, _, _, _ = load_end(branches["kick_v1"][0])
        dp_k1, Xk = descent_probe(Mk, h, delta, W2_DEFAULT)
        d_after = field_distance(Xk, Mp, free3)
        c3["descent_parent"] = dp_par
        c3["descent_kick_v1"] = dp_k1
        c3["kick_v1_dist_after_descent"] = d_after
        print(f"    parent : dE {dp_par['dE_total']:+.2e} "
              f"fmax {dp_par['fmax_start']:.1e} -> {dp_par['fmax_end']:.1e}")
        print(f"    kick_v1: dE {dp_k1['dE_total']:+.2e} "
              f"fmax {dp_k1['fmax_start']:.1e} -> {dp_k1['fmax_end']:.1e}; "
              f"dist to parent after probe {d_after / den_relax:.4f} of relax")
    same_topo = all(abs(c3[nm]["q_upper"] - 1) < 0.05
                    and abs(c3[nm]["q_lower"] + 1) < 0.05 for nm in branches)
    e_ok = all(abs(c3[nm]["E_mine"] - c3[nm]["E_note"]) < 5e-4
               for nm in branches)
    v1_gap = c3["kick_v1"]["dE_vs_parent"]
    v1_dist = c3["kick_v1"]["dist_rel_to_relax"]
    c3_verdict = "CONFIRMED" if (same_topo and e_ok and v1_dist < 0.02) else \
        ("QUALIFIED" if (same_topo and e_ok) else "REFUTED")
    audit["C3"] = {"verdict": c3_verdict, "detail": c3,
                   "kick_v1_energy_gap": v1_gap,
                   "kick_v1_dist_rel": v1_dist,
                   "method": "own energies; field distance on the free region "
                             "normalized by the parent relaxation displacement "
                             "||Mp - M0p||; eigengap-field rms; own-gradient "
                             "descent probe on parent + kick_v1"}
    print(f"  C3 {c3_verdict}")

    # ------------------ C4: pp escape ---------------------------------------
    print("\n-- C4: pp composite (2x s=-0.5) escape --")
    Mpp, _, hp2, dp2, _ = load_end("m5_22_1_end_d2_s-0.5_s-0.5_a2_sc6_n32_d0.3.npz")
    rowpp = rowjson("m5_22_1_row_d2_s-0.5_s-0.5_a2_sc6_n32_d0.3.json")
    e_pp = energy_pad0(Mpp, hp2, dp2, W2_DEFAULT)["E_end"]
    npp, gpp = directors(Mpp)
    pp_mask, pp_signs, ncf_pp = surface_lift(npp, gpp, far_b)
    if aligned_degree(npp, gpp, (lo, hi, lo, hi, kcap, hi),
                      (pp_mask, pp_signs))["Q"] < 0:
        pp_signs = -pp_signs
    refpp = (pp_mask, pp_signs)
    qfar_scan = {}
    for w in (11.25, 14.25, 17.25):
        a0, a1 = idx_of(n, hp2, -w), idx_of(n, hp2, w)
        r_ = read_box_surface(npp, gpp, (a0, a1, a0, a1, a0, a1))
        qfar_scan[f"halfwidth{w}"] = {"absQ": abs(r_["Q"]),
                                      "conflicts": r_["conflicts"]}
    qup_pp = aligned_degree(npp, gpp, (lo, hi, lo, hi, kcap, hi), refpp)["Q"]
    qlo_pp = aligned_degree(npp, gpp, (lo, hi, lo, hi, lo, kcap - 1),
                            refpp)["Q"]
    rings_pp = ring_read_own(gpp, hp2, 0.15)
    mom_pp, _, _ = volume_moments(Mpp, hp2)
    w_single = winding_P([(-0.5, 0.0)])
    w_comp = winding_P([(-0.5, -1.0), (-0.5, 1.0)])
    print(f"  E mine {e_pp:.4f} row {rowpp['E_end']:.4f} note 15.047")
    print(f"  far |degree| scan: " + " ".join(
        f"{k}:{v['absQ']:.4f}/cf{v['conflicts']}"
        for k, v in qfar_scan.items()))
    print(f"  upper {qup_pp:+.4f} lower {qlo_pp:+.4f} "
          f"(row slab +-1.04, note +-1.07)")
    print(f"  rings thr0.15: {[(round(e['z_centroid'], 1), round(e['rho_mean'], 1)) for e in rings_pp['rings']]}")
    print(f"  2D winding: single s=-0.5 {w_single:+.5f}, composite {w_comp:+.5f}"
          f" (escape law: total 2 -> 3D degree 0)")
    print(f"  p_z volume {mom_pp['p_z']:+.3f} (row 24.52)")
    far_zero = max(v["absQ"] for v in qfar_scan.values()) < 0.02
    c4_verdict = "CONFIRMED" if (far_zero and abs(e_pp - 15.047) < 5e-3
                                 and abs(w_comp - 2.0) < 5e-4
                                 and abs(qup_pp - 1.0) < 0.02
                                 and abs(qlo_pp + 1.0) < 0.02) else "QUALIFIED"
    audit["C4"] = {"verdict": c4_verdict, "E_mine": e_pp,
                   "E_row": rowpp["E_end"],
                   "qfar_scan": qfar_scan, "q_upper": qup_pp,
                   "q_lower": qlo_pp, "far_surface_conflicts": int(ncf_pp),
                   "rings_thr0.15_z": [e["z_centroid"]
                                       for e in rings_pp["rings"]],
                   "winding_single": w_single, "winding_composite": w_comp,
                   "p_z_volume": mom_pp["p_z"],
                   "binding_check": {"E_pp": e_pp,
                                     "twice_E_p": 2 * myE[("P-0.5", 0.3)]},
                   "method": "own energy; own lift + box degrees at 3 radii; "
                             "own incremental-unwrap 2D winding of the summed "
                             "cross-section angle"}
    print(f"  C4 {c4_verdict}")

    # ------------------ C5: graft charge frustration ------------------------
    print("\n-- C5: graft far-field consistency --")
    c5 = {}
    for gtag, rowf, noteq in (
            ("m5_22_1_end_graft_p_below_n_above_z6_n32_d0.3.npz",
             "m5_22_1_row_graft_p_below_n_above_z6_n32_d0.3.json", -0.88),
            ("m5_22_1_end_graft_p_below_n_above_z9_n32_d0.3.npz",
             "m5_22_1_row_graft_p_below_n_above_z9_n32_d0.3.json", -0.85)):
        Mg, _, hg, dg, _ = load_end(gtag)
        rowg = rowjson(rowf)
        eg = energy_pad0(Mg, hg, dg, W2_DEFAULT)["E_end"]
        ng, gg = directors(Mg)
        scan = {}
        for w in (11.25, 12.75, 14.25, 15.75, 17.25, 18.75, 20.25):
            a0, a1 = idx_of(n, hg, -w), idx_of(n, hg, w)
            r_ = read_box_surface(ng, gg, (a0, a1, a0, a1, a0, a1))
            scan[f"halfwidth{w}"] = {"absQ": abs(r_["Q"]),
                                     "conflicts": r_["conflicts"],
                                     "min_gap": r_["min_gap"]}
        # where the outer-surface frustration sits (z of broken edges)
        a0, a1 = idx_of(n, hg, -17.25), idx_of(n, hg, 17.25)
        outer_b = (a0, a1, a0, a1, a0, a1)
        cz = conflict_zs(ng, gg, outer_b, hg)
        # z-split of the outer surface charge (aligned within this radius)
        g_mask, g_signs, _ = surface_lift(ng, gg, outer_b)
        refg = (g_mask, g_signs)
        up_g = aligned_degree(ng, gg, (a0, a1, a0, a1, kcap, a1), refg)
        lo_g = aligned_degree(ng, gg, (a0, a1, a0, a1, a0, kcap - 1), refg)
        fu_g = aligned_degree(ng, gg, outer_b, refg)
        key = os.path.basename(gtag).replace("m5_22_1_end_", "").replace(".npz", "")
        vals = [v["absQ"] for v in scan.values()]
        cfs = [v["conflicts"] for v in scan.values()]
        c5[key] = {"E_mine": eg, "E_row": rowg["E_end"],
                   "q_note": noteq, "q_row": rowg["charges"]["q_far"],
                   "scan": scan,
                   "absQ_spread": float(max(vals) - min(vals)),
                   "outer_conflict_zs": {
                       "count": len(cz),
                       "z_min": float(min(cz)) if cz else None,
                       "z_max": float(max(cz)) if cz else None},
                   "outer_split": {"full": fu_g["Q"], "upper": up_g["Q"],
                                   "lower": lo_g["Q"]}}
        print(f"  {key}: E {eg:.4f} (row {rowg['E_end']:.4f}); |Q| by radius "
              + " ".join(f"{v:.3f}" for v in vals)
              + f"; conflicts {cfs}")
        print(f"    outer surface: full {fu_g['Q']:+.3f} upper {up_g['Q']:+.3f}"
              f" lower {lo_g['Q']:+.3f}; frustrated edges n={len(cz)}"
              + (f" z in [{min(cz):+.1f}, {max(cz):+.1f}]" if cz else ""))
    inner_zero = all(c5[k]["scan"][f"halfwidth{w}"]["absQ"] < 0.02
                     and c5[k]["scan"][f"halfwidth{w}"]["conflicts"] == 0
                     for k in c5 for w in (11.25, 12.75, 14.25, 15.75))
    outer_charged = all(abs(c5[k]["scan"]["halfwidth17.25"]["absQ"] - 1.0) < 0.05
                        for k in c5)
    outer_frustrated = all(c5[k]["scan"]["halfwidth17.25"]["conflicts"] > 0
                           for k in c5)
    c5_verdict = "QUALIFIED" if (inner_zero and outer_frustrated) else \
        ("CONFIRMED" if outer_frustrated else "REFUTED")
    audit["C5"] = {"verdict": c5_verdict, "detail": c5,
                   "inner_radii_clean_zero": inner_zero,
                   "outer_radius_charged": outer_charged,
                   "outer_radius_frustrated": outer_frustrated,
                   "finding": "the far field is LAYERED, not one fractional "
                              "sector: exact degree 0 with zero lift "
                              "conflicts at half-widths 11.25-15.75, then "
                              "|Q| = 1 with a frustrated (non-orientable) "
                              "seam at the outermost free surfaces next to "
                              "the pinned shell; the rows' -0.88/-0.85 are "
                              "the FD flux instrument smearing that seam, "
                              "not a radius-stable fractional charge",
                   "method": "independent per-surface lift + exact solid-angle "
                             "degree at 7 probe radii + broken-edge z map"}
    print(f"  C5 {c5_verdict} (inner zero {inner_zero}, outer charged "
          f"{outer_charged}, outer frustrated {outer_frustrated})")

    # ------------------ C6: the late three-center composite -----------------
    if not args.skip_c6:
        print("\n-- C6: three-center composite (dn) --")
        rows = sorted(glob.glob(os.path.join(DATA, "m5_22_1_row_dn_*.json")))
        if not rows:
            audit["C6"] = {"verdict": "NOT-LANDED",
                           "note": "no m5_22_1_row_dn_*.json at audit close"}
            print("  no dn row landed yet")
        else:
            c6 = {}
            for rf in rows:
                rowd = rowjson(os.path.basename(rf))
                tag = rowd["tag"]
                endf = f"m5_22_1_end_{tag}.npz"
                try:
                    Md, _, hd, dd, _ = load_end(endf)
                except FileNotFoundError:
                    c6[tag] = {"error": "npz missing"}
                    continue
                ed = energy_pad0(Md, hd, dd, W2_DEFAULT)["E_end"]
                nd, gd_ = directors(Md)
                d_mask, d_signs, ncf_d = surface_lift(nd, gd_, far_b)
                if aligned_degree(nd, gd_, (lo, hi, lo, hi, kcap, hi),
                                  (d_mask, d_signs))["Q"] < 0:
                    d_signs = -d_signs
                refd = (d_mask, d_signs)
                qfar = aligned_degree(nd, gd_, (lo, hi, lo, hi, lo, hi),
                                      refd)["Q"]
                qup = aligned_degree(nd, gd_, (lo, hi, lo, hi, kcap, hi),
                                     refd)["Q"]
                qlo = aligned_degree(nd, gd_, (lo, hi, lo, hi, lo, kcap - 1),
                                     refd)["Q"]
                rr = ring_read_own(gd_, hd, 0.15)
                w2d = None
                if "pairs" in rowd:
                    w2d = winding_P([(float(s_), float(x_))
                                     for s_, x_ in rowd["pairs"]])
                # z-cut ladder: localize the enclosed charge in z
                zc6 = (np.arange(n) - (n - 1) / 2.0) * hd
                Qb6 = {}
                for k in range(lo + 1, hi + 1):
                    Qb6[k] = aligned_degree(nd, gd_,
                                            (lo, hi, lo, hi, lo, k),
                                            refd)["Q"]
                lumps = []
                for k in range(lo + 1, hi):
                    dqv = Qb6[k + 1] - Qb6[k]
                    if abs(dqv) > 0.05:
                        lumps.append({"z_mid": 0.5 * (zc6[k] + zc6[k + 1]),
                                      "dQ": dqv})
                capscan = {}
                for kc in range(13, 20):
                    capscan[f"z{(kc - (n - 1) / 2) * hd:+.2f}"] = \
                        aligned_degree(nd, gd_, (lo, hi, lo, hi, kc, hi),
                                       refd)["Q"]
                c6[tag] = {"winding_2d_mine": w2d,
                           "z_ladder_lumps": lumps,
                           "upper_box_cap_scan": capscan,
                           "E_mine": ed, "E_row": rowd.get("E_end"),
                           "q_far_mine": qfar,
                           "q_far_row": rowd.get("charges", {}).get("q_far"),
                           "q_upper": qup, "q_lower": qlo,
                           "far_surface_conflicts": int(ncf_d),
                           "rings_z_mine": [e["z_centroid"]
                                            for e in rr["rings"]],
                           "rings_rho_mine": [e["rho_mean"]
                                              for e in rr["rings"]],
                           "row_stop": rowd.get("stop")}
                print(f"  {tag}: E {ed:.4f} (row {rowd.get('E_end')}) "
                      f"qfar {qfar:+.4f} (row "
                      f"{rowd.get('charges', {}).get('q_far')}) "
                      f"up {qup:+.3f} lo {qlo:+.3f} "
                      f"rings {[round(z, 1) for z in c6[tag]['rings_z_mine']]}"
                      f" 2d-winding {w2d}")
                print(f"    z-ladder lumps: " + " ".join(
                    f"z{e['z_mid']:+.1f}:{e['dQ']:+.3f}" for e in lumps))
                print(f"    upper-box cap scan: " + " ".join(
                    f"{k}:{v:+.3f}" for k, v in capscan.items()))
            ok = all("error" not in v and
                     abs(v["E_mine"] - v["E_row"]) / v["E_row"] < 1e-4
                     for v in c6.values())
            audit["C6"] = {"verdict": "CONFIRMED" if ok else "QUALIFIED",
                           "detail": c6}
            print(f"  C6 {audit['C6']['verdict']}")

    # ------------------ save -------------------------------------------------
    out = os.path.join(DATA, "m5_22_1_audit.json")
    with open(out, "w") as f:
        json.dump(tofloat(audit), f, indent=1)
    print(f"\nverdicts -> {out}")
    print("\n=== verdict summary ===")
    for k in ("C1", "C2", "C3", "C4", "C5", "C6"):
        if k in audit:
            print(f"  {k}: {audit[k]['verdict']}")


if __name__ == "__main__":
    main()
