"""M5.21.4 independent adversarial audit, rounds 1 + 2.

Round 1 (claims 1, 2, 3, 6, 7, 8): static claims.
Round 2 (R2-1..R2-5): the capture-run claims (note section 5). CLI:
    python3 m5_21_4_audit_check.py            round 1 (unchanged)
    python3 m5_21_4_audit_check.py round2     compute round 2 -> scratch
    python3 m5_21_4_audit_check.py merge2     merge verdicts into
                                              ../data/m5_21_4_audit.json


Independent implementations throughout (none of the task's functions are
used for audited quantities):
  energy      own stencil code (fwd/bwd slicing, own commutator sum, own
              eigvalsh potential), density variant for regional splits
  charge      own oriented lift (BFS flood fill, not sweep+majority), own
              Mermin-Ho B (np.gradient), own trapezoid-weighted cube-face
              flux (different halves than the task's), own spherical
              solid-angle degree on (a) interpolated lattice directors and
              (b) EXACT off-lattice evaluation of the analytic seed maps
  core scan   deviator-norm dip along near-axis columns (not the task's
              eigen-gap tracker) + full-grid low-gap locus
Seed construction is reused via P.seed_pair per the audit brief; local
parameterized copies (r_t tube radius, r_c blend radius) exist ONLY for
the artifact probes and are verified equal to P.seed_pair at defaults.

Writes ../data/m5_21_4_audit.json. Read-only on all task artifacts.
"""

import importlib.util
import json
import os
import sys
import time
from collections import deque

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")

_s = importlib.util.spec_from_file_location(
    "pair", os.path.join(HERE, "m5_21_4_a_pair.py"))
P = importlib.util.module_from_spec(_s)
_s.loader.exec_module(P)
INS = P.INS

W2 = 0.000724023879        # the 2b w2 (note quotes 7.2402e-4)
DELTA = 0.3
VAC = np.array([0.0, DELTA, 1.0])   # ascending


# ================= my energy =================
def my_energy(M, h, w2=W2, split=False):
    """own implementation: E = h^3 [ 4 sum_{i<j} ||[D_i M, D_j M]||_F^2
    (avg of fwd/bwd branch densities) + w2 sum_k (lam_k - vac_k)^2 ]."""
    M = np.asarray(M, dtype=np.float64)
    e_u = 0.0
    for br in ("f", "b"):
        A = []
        for ax in range(3):
            D = np.zeros_like(M)
            hi = [slice(None)] * 5
            lo = [slice(None)] * 5
            hi[ax], lo[ax] = slice(1, None), slice(0, -1)
            if br == "f":
                D[tuple(lo)] = (M[tuple(hi)] - M[tuple(lo)]) / h
            else:
                D[tuple(hi)] = (M[tuple(hi)] - M[tuple(lo)]) / h
            A.append(D)
        for i in range(3):
            for j in range(i + 1, 3):
                C = A[i] @ A[j] - A[j] @ A[i]
                e_u += 0.5 * 4.0 * float((C * C).sum())
        del A
    lam = np.linalg.eigvalsh(M)
    e_v = w2 * float(((lam - VAC) ** 2).sum())
    if split:
        return h ** 3 * e_u, h ** 3 * e_v
    return h ** 3 * (e_u + e_v)


def my_edens(M, h, w2=W2):
    """per-voxel density (same definition), for regional splits."""
    M = np.asarray(M, dtype=np.float64)
    dens = np.zeros(M.shape[:3])
    for br in ("f", "b"):
        A = []
        for ax in range(3):
            D = np.zeros_like(M)
            hi = [slice(None)] * 5
            lo = [slice(None)] * 5
            hi[ax], lo[ax] = slice(1, None), slice(0, -1)
            if br == "f":
                D[tuple(lo)] = (M[tuple(hi)] - M[tuple(lo)]) / h
            else:
                D[tuple(hi)] = (M[tuple(hi)] - M[tuple(lo)]) / h
            A.append(D)
        for i in range(3):
            for j in range(i + 1, 3):
                C = A[i] @ A[j] - A[j] @ A[i]
                dens += 0.5 * 4.0 * np.einsum("...ab,...ab->...", C, C)
        del A
    lam = np.linalg.eigvalsh(M)
    dens += w2 * ((lam - VAC) ** 2).sum(axis=-1)
    return dens


# ================= my oriented lift (BFS flood fill) =================
def my_orient(M):
    lam, vec = np.linalg.eigh(np.asarray(M, dtype=np.float64))
    v = vec[..., :, 2].copy()
    n = v.shape[0]
    sgn = np.zeros((n, n, n), dtype=np.int8)
    start = (1, 1, 1)
    sgn[start] = 1
    dq = deque([start])
    nbrs = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0),
            (0, 0, 1), (0, 0, -1))
    while dq:
        i, j, k = dq.popleft()
        vi = v[i, j, k] * sgn[i, j, k]
        for di, dj, dk in nbrs:
            a, b, c = i + di, j + dj, k + dk
            if 0 <= a < n and 0 <= b < n and 0 <= c < n and sgn[a, b, c] == 0:
                d = vi[0] * v[a, b, c, 0] + vi[1] * v[a, b, c, 1] \
                    + vi[2] * v[a, b, c, 2]
                sgn[a, b, c] = 1 if d >= 0.0 else -1
                dq.append((a, b, c))
    vo = v * sgn[..., None].astype(np.float64)
    confs = []
    for ax in range(3):
        sl_hi = [slice(None)] * 4
        sl_lo = [slice(None)] * 4
        sl_hi[ax], sl_lo[ax] = slice(1, None), slice(0, -1)
        dts = np.einsum("...a,...a->...", vo[tuple(sl_hi)], vo[tuple(sl_lo)])
        for idx in np.argwhere(dts < -1e-9):
            p = idx.astype(float)
            p[ax] += 0.5
            confs.append(p)
    return vo, confs


def conflict_stats(confs, h, n):
    if not confs:
        return {"count": 0}
    arr = (np.array(confs) - (n - 1) / 2.0) * h
    rho = np.hypot(arr[:, 0], arr[:, 1])
    az = np.abs(arr[:, 2])
    in_tube = float(np.mean((rho < 5.5) & (az < 10.5)))
    return {"count": len(confs), "rho_min": float(rho.min()),
            "rho_max": float(rho.max()), "absz_max": float(az.max()),
            "frac_in_tube_rho5.5_z10.5": in_tube}


# ================= my flux instruments =================
def my_B(nv, h):
    g = [np.gradient(nv, h, axis=ax) for ax in range(3)]
    Fyz = np.einsum("...a,...a->...", nv, np.cross(g[1], g[2]))
    Fzx = np.einsum("...a,...a->...", nv, np.cross(g[2], g[0]))
    Fxy = np.einsum("...a,...a->...", nv, np.cross(g[0], g[1]))
    return np.stack([Fyz, Fzx, Fxy], axis=-1)


def my_cube_flux(B, h, center, half):
    """trapezoid-weighted closed-cube flux, own face bookkeeping."""
    n = B.shape[0]
    c = (n - 1) / 2.0
    ic = [int(round(center[a] / h + c)) for a in range(3)]
    k = int(round(half / h))
    lo = [ic[a] - k for a in range(3)]
    hi = [ic[a] + k for a in range(3)]
    if min(lo) < 1 or max(hi) > n - 2:
        return float("nan")
    m = 2 * k + 1
    w1 = np.ones(m)
    w1[0] = w1[-1] = 0.5
    W = np.outer(w1, w1)
    s = 0.0
    for ax in range(3):
        sl_hi = [slice(lo[a], hi[a] + 1) for a in range(3)] + [ax]
        sl_lo = [slice(lo[a], hi[a] + 1) for a in range(3)] + [ax]
        sl_hi[ax], sl_lo[ax] = hi[ax], lo[ax]
        s += float((B[tuple(sl_hi)] * W).sum())
        s -= float((B[tuple(sl_lo)] * W).sum())
    return s * h * h / (4.0 * np.pi)


def _trilin(F, idx):
    n = F.shape[0]
    idx = np.clip(idx, 0.0, n - 1 - 1e-9)
    I = np.floor(idx).astype(int)
    f = idx - I
    out = np.zeros(idx.shape[:-1] + (3,))
    for dx in (0, 1):
        wx = f[..., 0] if dx else 1.0 - f[..., 0]
        for dy in (0, 1):
            wy = f[..., 1] if dy else 1.0 - f[..., 1]
            for dz in (0, 1):
                wz = f[..., 2] if dz else 1.0 - f[..., 2]
                out += (wx * wy * wz)[..., None] * F[
                    np.minimum(I[..., 0] + dx, n - 1),
                    np.minimum(I[..., 1] + dy, n - 1),
                    np.minimum(I[..., 2] + dz, n - 1)]
    return out


def sphere_degree_lattice(vf, h, center, R, nth=180, nph=360):
    """solid-angle degree of the (oriented) lattice director on a sphere,
    via trilinear interpolation; my own surface shape (not a cube)."""
    n = vf.shape[0]
    th = (np.arange(nth) + 0.5) * np.pi / nth
    ph = (np.arange(nph) + 0.5) * 2.0 * np.pi / nph
    TH, PH = np.meshgrid(th, ph, indexing="ij")
    pts = np.stack([center[0] + R * np.sin(TH) * np.cos(PH),
                    center[1] + R * np.sin(TH) * np.sin(PH),
                    center[2] + R * np.cos(TH)], axis=-1)
    nv = _trilin(vf, pts / h + (n - 1) / 2.0)
    nv /= np.maximum(np.linalg.norm(nv, axis=-1, keepdims=True), 1e-30)
    return _deg_from_grid(nv, th, ph)


def _deg_from_grid(nv, th, ph):
    dth = np.gradient(nv, th, axis=0)
    dpv = 2.0 * np.pi / len(ph)
    dph = (np.roll(nv, -1, axis=1) - np.roll(nv, 1, axis=1)) / (2.0 * dpv)
    J = np.einsum("...a,...a->...", nv, np.cross(dth, dph))
    return float(J.sum() * (np.pi / len(th)) * dpv / (4.0 * np.pi))


# ---- exact off-lattice evaluation of the analytic seed director maps ----
def nhat_anti_exact(x, y, z, d):
    rho = np.hypot(x, y)
    a = np.arctan2(rho, z - d / 2.0) + np.pi - np.arctan2(rho, z + d / 2.0)
    rhos = np.where(rho < 1e-12, 1e-12, rho)
    cph, sph = x / rhos, y / rhos
    return np.stack([np.sin(a) * cph, np.sin(a) * sph, np.cos(a)], axis=-1)


def nhat_single_exact(x, y, z, mirror=False):
    rho = np.hypot(x, y)
    a = np.arctan2(rho, z)
    if mirror:
        a = np.pi - a
    rhos = np.where(rho < 1e-12, 1e-12, rho)
    cph, sph = x / rhos, y / rhos
    return np.stack([np.sin(a) * cph, np.sin(a) * sph, np.cos(a)], axis=-1)


def nhat_same_exact(x, y, z, d, r_t=3.0):
    rho = np.hypot(x, y)
    zt, zb = d / 2.0, -d / 2.0
    tht = np.arctan2(rho, z - zt)
    thb = np.arctan2(rho, z - zb)
    tt = np.clip(np.tan(0.5 * tht), 0.0, 1e6)
    tb = np.clip(np.tan(0.5 * thb), 0.0, 1e6)
    umag = np.clip(tt * tb, 0.0, 1e12)
    rhos = np.where(rho < 1e-12, 1e-12, rho)
    c2p = (x * x - y * y) / (rhos * rhos)
    s2p = -(2.0 * x * y) / (rhos * rhos)
    den = 1.0 + umag * umag
    nv = np.stack([2.0 * umag * c2p / den, 2.0 * umag * s2p / den,
                   (1.0 - umag * umag) / den], axis=-1)
    zed = 0.5 * (np.tanh((z - zb) / 2.0) - np.tanh((z - zt) / 2.0))
    eta = np.exp(-((rho / r_t) ** 2)) * zed
    nv = nv + eta[..., None] * np.array([0.0, 0.0, 1.0])
    return nv / np.maximum(np.linalg.norm(nv, axis=-1, keepdims=True),
                           1e-300)


def sphere_degree_exact(fn, center, R, nth=400, nph=800):
    th = (np.arange(nth) + 0.5) * np.pi / nth
    ph = (np.arange(nph) + 0.5) * 2.0 * np.pi / nph
    TH, PH = np.meshgrid(th, ph, indexing="ij")
    x = center[0] + R * np.sin(TH) * np.cos(PH)
    y = center[1] + R * np.sin(TH) * np.sin(PH)
    z = center[2] + R * np.cos(TH)
    nv = fn(x, y, z)
    return _deg_from_grid(nv, th, ph)


# ================= my core scan =================
def axis_deviator(M, cfg):
    """deviator-norm profile averaged over the 4 near-axis columns."""
    M = np.asarray(M, dtype=np.float64)
    n, h = cfg["n"], cfg["h"]
    c = (n - 1) // 2
    prof = np.zeros(n)
    for (i, j) in ((c, c), (c, c + 1), (c + 1, c), (c + 1, c + 1)):
        col = M[i, j]
        dev = col - (np.trace(col, axis1=-2, axis2=-1) / 3.0)[..., None,
                                                             None] * np.eye(3)
        prof += np.sqrt(np.einsum("...ab,...ab->...", dev, dev))
    prof /= 4.0
    zax = (np.arange(n) - (n - 1) / 2.0) * h
    return zax, prof


def find_dips(zax, prof, bulk=0.72, frac=0.5, margin=4):
    """z locations where the deviator norm dips below frac*bulk (core
    signature: isotropic blend), boundary margin excluded."""
    sel = slice(margin, len(prof) - margin)
    lo = prof[sel] < frac * bulk
    zs = zax[sel][lo]
    return [float(z) for z in zs], float(prof[sel].min())


def lowgap_locus(M, cfg, thr=0.12):
    M = np.asarray(M, dtype=np.float64)
    lam = np.linalg.eigvalsh(M)
    gap = np.minimum(lam[..., 2] - lam[..., 1], lam[..., 1] - lam[..., 0])
    X, Y, Z = INS.coords(cfg["n"], cfg["h"])
    m = gap < thr
    # exclude the pinned/edge shell
    r_inf = np.maximum(np.maximum(np.abs(X), np.abs(Y)), np.abs(Z))
    m &= r_inf < 0.5 * cfg["L"] - 3 * cfg["h"]
    if not m.any():
        return {"count": 0}
    rho = np.hypot(X[m], Y[m])
    return {"count": int(m.sum()), "rho_med": float(np.median(rho)),
            "rho_min": float(rho.min()), "rho_max": float(rho.max()),
            "z_med": float(np.median(Z[m])),
            "absz_max": float(np.abs(Z[m]).max())}


# ================= parameterized construction copies (probes only) ====
def my_tensor_from_nhat(cfg, nhat, centers, r_c=4.0):
    n, h, delta = cfg["n"], cfg["h"], cfg["delta"]
    X, Y, Z = INS.coords(n, h)
    rho = np.sqrt(X * X + Y * Y)
    rhos = np.where(rho < 1e-12, 1e-12, rho)
    phihat = np.stack([-Y / rhos, X / rhos, np.zeros_like(Z)], axis=-1)
    near = rho < 1e-9
    if np.any(near):
        phihat[near] = np.array([0.0, 1.0, 0.0])
    dot = np.einsum("...a,...a->...", phihat, nhat)[..., None]
    ph = phihat - dot * nhat
    ph = ph / np.maximum(np.linalg.norm(ph, axis=-1)[..., None], 1e-300)
    S = (nhat[..., :, None] * nhat[..., None, :]
         + delta * ph[..., :, None] * ph[..., None, :])
    a = (1.0 + delta) / 3.0
    w = np.ones_like(rho)
    for zc in centers:
        di = np.sqrt(X * X + Y * Y + (Z - zc) ** 2)
        w = w * (1.0 - np.exp(-((di / r_c) ** 2)))
    return w[..., None, None] * S \
        + (1.0 - w[..., None, None]) * (a * np.eye(3))


def my_seed(cfg, kind, d, r_t=3.0, r_c=4.0):
    n, h = cfg["n"], cfg["h"]
    X, Y, Z = INS.coords(n, h)
    if kind == "same":
        nv = nhat_same_exact(X, Y, Z, d, r_t=r_t)
        return my_tensor_from_nhat(cfg, nv, [d / 2.0, -d / 2.0], r_c=r_c)
    if kind == "anti":
        nv = nhat_anti_exact(X, Y, Z, d)
        return my_tensor_from_nhat(cfg, nv, [d / 2.0, -d / 2.0], r_c=r_c)
    if kind in ("single", "mirror"):
        nv = nhat_single_exact(X, Y, Z, mirror=(kind == "mirror"))
        return my_tensor_from_nhat(cfg, nv, [0.0], r_c=r_c)
    raise ValueError(kind)


def jload(name):
    with open(os.path.join(DATA, name)) as f:
        return json.load(f)


def npz_state(name):
    with np.load(os.path.join(DATA, name)) as z:
        return np.asarray(z["M"], dtype=np.float64)


# ================= main =================
def main():
    t00 = time.time()
    R = {"meta": {"date": "2026-07-21", "round": 1,
                  "scope": [1, 2, 3, 6, 7, 8],
                  "skipped": {"4": "m5_21_4_ev_cap.npz not on disk",
                              "5": "m5_21_4_ev_cap.npz not on disk"}}}
    cfg32 = INS.base_cfg(term="T2", n=32, L=48.0, bc="pinned",
                         stencil="sym")
    cfg48 = INS.base_cfg(term="T2", n=48, L=72.0, bc="pinned",
                         stencil="sym")
    h = cfg32["h"]

    # ---- self tests -------------------------------------------------
    st = {}
    n = 32
    X, Y, Z = INS.coords(n, h)
    rr = np.sqrt(X * X + Y * Y + Z * Z)
    rs = np.where(rr < 1e-12, 1e-12, rr)
    vhedge = np.stack([X / rs, Y / rs, Z / rs], axis=-1)
    st["sphere_deg_hedgehog_R10"] = sphere_degree_lattice(
        vhedge, h, (0, 0, 0), 10.0)
    Bh = my_B(vhedge, h)
    st["cube_flux_hedgehog_half9"] = my_cube_flux(Bh, h, (0, 0, 0), 9.0)
    st["exact_sphere_deg_single"] = sphere_degree_exact(
        lambda x, y, z: nhat_single_exact(x, y, z), (0, 0, 0), 8.0)
    # construction copies == task construction at defaults
    for kind, d in (("same", 18.0), ("anti", 18.0), ("single", 0.0)):
        dev = float(np.abs(my_seed(cfg32, kind, d)
                           - P.seed_pair(cfg32, kind, d)).max())
        st[f"copy_vs_task_{kind}"] = dev
    R["self_tests"] = st
    print("self tests", json.dumps(st, indent=1), flush=True)

    # ---- claim 1: signed-charge instrument --------------------------
    c1 = {}
    seeds = {k: P.seed_pair(cfg32, k, d) for k, d in
             (("single", 0.0), ("mirror", 0.0), ("anti", 18.0),
              ("same", 18.0))}
    # exact analytic degrees (off-lattice, integer-quantized ground truth)
    ex = {}
    ex["single_far_R8"] = sphere_degree_exact(
        lambda x, y, z: nhat_single_exact(x, y, z), (0, 0, 0), 8.0)
    ex["single_far_R16"] = sphere_degree_exact(
        lambda x, y, z: nhat_single_exact(x, y, z), (0, 0, 0), 16.0)
    ex["mirror_far_R8"] = sphere_degree_exact(
        lambda x, y, z: nhat_single_exact(x, y, z, mirror=True),
        (0, 0, 0), 8.0)
    for lab, c, rad in (("top", (0, 0, 9.0), 5.0), ("bot", (0, 0, -9.0),
                                                    5.0),
                        ("far", (0, 0, 0), 20.0)):
        ex[f"anti_{lab}_R{rad:g}"] = sphere_degree_exact(
            lambda x, y, z: nhat_anti_exact(x, y, z, 18.0), c, rad)
        ex[f"same_{lab}_R{rad:g}"] = sphere_degree_exact(
            lambda x, y, z: nhat_same_exact(x, y, z, 18.0), c, rad)
    ex["same_top_R4"] = sphere_degree_exact(
        lambda x, y, z: nhat_same_exact(x, y, z, 18.0), (0, 0, 9.0), 4.0)
    ex["same_far_R16"] = sphere_degree_exact(
        lambda x, y, z: nhat_same_exact(x, y, z, 18.0), (0, 0, 0), 16.0)
    c1["exact_degrees"] = ex
    # lattice instrument: my BFS lift + my B + my cubes + my spheres
    latt = {}
    for kind, M in seeds.items():
        vo, confs = my_orient(M)
        # global sign: align far field with the analytic map
        Xg, Yg, Zg = INS.coords(32, h)
        if kind == "same":
            ref = nhat_same_exact(Xg, Yg, Zg, 18.0)
        elif kind == "anti":
            ref = nhat_anti_exact(Xg, Yg, Zg, 18.0)
        else:
            ref = nhat_single_exact(Xg, Yg, Zg, mirror=(kind == "mirror"))
        far = rr > 12.0
        if float(np.sum(np.einsum("...a,...a->...", vo, ref)[far])) < 0:
            vo = -vo
        B = my_B(vo, h)
        row = {"conflicts": conflict_stats(confs, h, 32)}
        d = 18.0
        row["cube_top_h4.5"] = my_cube_flux(B, h, (0, 0, d / 2), 4.5)
        row["cube_top_h6"] = my_cube_flux(B, h, (0, 0, d / 2), 6.0)
        row["cube_bot_h4.5"] = my_cube_flux(B, h, (0, 0, -d / 2), 4.5)
        row["cube_bot_h6"] = my_cube_flux(B, h, (0, 0, -d / 2), 6.0)
        row["cube_far_h16.5"] = my_cube_flux(B, h, (0, 0, 0), 16.5)
        row["cube_far_h19.5"] = my_cube_flux(B, h, (0, 0, 0), 19.5)
        row["sphere_far_R16"] = sphere_degree_lattice(vo, h, (0, 0, 0),
                                                      16.0)
        row["sphere_far_R20"] = sphere_degree_lattice(vo, h, (0, 0, 0),
                                                      20.0)
        if kind in ("anti", "same"):
            row["sphere_top_R5"] = sphere_degree_lattice(
                vo, h, (0, 0, d / 2), 5.0)
            row["sphere_bot_R5"] = sphere_degree_lattice(
                vo, h, (0, 0, -d / 2), 5.0)
        latt[kind] = row
        print("claim1", kind, json.dumps(row), flush=True)
    c1["lattice_instrument"] = latt
    R["claim1"] = c1

    # ---- claim 2: it400 saved-state energy re-reads ------------------
    c2 = {"rows": []}
    j400 = jload("m5_21_4_ladder_it400.json")
    for kind, d in (("anti", 18.0), ("anti", 24.0), ("same", 18.0)):
        M = npz_state(f"m5_21_4_lad_{kind}_d{d:g}_n32_it400.npz")
        eu, ev = my_energy(M, h, split=True)
        row = next(r for r in j400["rows"]
                   if r["kind"] == kind and r["d"] == d)
        c2["rows"].append({
            "kind": kind, "d": d, "E_claimed": row["E"],
            "E_mine": eu + ev, "E_u_claimed": row["E_u"], "E_u_mine": eu,
            "rel_dev": abs(eu + ev - row["E"]) / abs(row["E"])})
        print("claim2", json.dumps(c2["rows"][-1]), flush=True)
    R["claim2"] = c2

    # ---- claim 3: slope signs + core stationarity --------------------
    c3 = {}
    sl = jload("m5_21_4_seed_ladder.json")

    def eint(sec, box):
        ref = sl[box]["E_single"] + sl[box]["E_mirror"]
        return {r["d"]: r["E"] - ref for r in sl[box]["rows"]
                if r["kind"] == sec}
    a32, s32 = eint("anti", "n32"), eint("same", "n32")
    ds = sorted(a32)
    c3["anti_Eint_n32"] = {str(d): a32[d] for d in ds}
    c3["anti_dEdd_all_positive"] = bool(
        all(a32[ds[i + 1]] > a32[ds[i]] for i in range(len(ds) - 1)))
    dss = sorted(s32)
    c3["same_Eint_n32"] = {str(d): s32[d] for d in dss}
    c3["same_dEdd_all_positive"] = bool(
        all(s32[dss[i + 1]] > s32[dss[i]] for i in range(len(dss) - 1)))
    # my core scans: seed control, it120 half-gone, it400 gone
    scans = {}
    Mc = P.seed_pair(cfg32, "anti", 18.0)
    zax, prof = axis_deviator(Mc, cfg32)
    zs, mn = find_dips(zax, prof)
    scans["seed_anti_d18"] = {"dip_zs": zs, "min_dev": mn}
    for tag, fn in (("anti_d18_it400", "m5_21_4_lad_anti_d18_n32_it400.npz"),
                    ("anti_d24_it400", "m5_21_4_lad_anti_d24_n32_it400.npz"),
                    ("anti_d24_it120",
                     "m5_21_4_lad_anti_d24_n32_it120.npz")):
        M = npz_state(fn)
        zax, prof = axis_deviator(M, cfg32)
        zs, mn = find_dips(zax, prof)
        scans[tag] = {"dip_zs": zs, "min_dev": mn}
    # it400 charge re-reads (my instruments) on anti d24 + same d18
    M = npz_state("m5_21_4_lad_anti_d24_n32_it400.npz")
    vo, confs = my_orient(M)
    B = my_B(vo, h)
    scans["anti_d24_it400_farflux_mine"] = my_cube_flux(B, h, (0, 0, 0),
                                                        16.5)
    scans["anti_d24_it400_conflicts"] = conflict_stats(confs, h, 32)
    Ms = npz_state("m5_21_4_lad_same_d18_n32_it400.npz")
    vo, confs = my_orient(Ms)
    B = my_B(vo, h)
    scans["same_d18_it400_cube_half9_mine"] = my_cube_flux(B, h, (0, 0, 0),
                                                           9.0)
    scans["same_d18_it400_cube_half15_mine"] = my_cube_flux(
        B, h, (0, 0, 0), 15.0)
    scans["same_d18_it400_cube_half16.5_mine"] = my_cube_flux(
        B, h, (0, 0, 0), 16.5)
    scans["same_d18_it400_cube_half19.5_mine"] = my_cube_flux(
        B, h, (0, 0, 0), 19.5)
    scans["same_d18_it400_sphere_R20_mine"] = sphere_degree_lattice(
        vo, h, (0, 0, 0), 20.0)
    scans["same_d18_it400_conflicts"] = conflict_stats(confs, h, 32)
    scans["same_d18_it400_lowgap"] = lowgap_locus(Ms, cfg32)
    c3["scans"] = scans
    R["claim3"] = c3
    print("claim3", json.dumps(c3, indent=1), flush=True)

    # ---- claim 7: seed-ladder energy re-reads (own functional) -------
    c7 = {"rows": []}
    checks32 = [("single", 0.0, sl["n32"]["E_single"]),
                ("mirror", 0.0, sl["n32"]["E_mirror"]),
                ("anti", 12.0, None), ("anti", 24.0, None),
                ("same", 12.0, None), ("same", 24.0, None)]
    dens_cache = {}
    for kind, d, claimed in checks32:
        if claimed is None:
            claimed = next(r["E"] for r in sl["n32"]["rows"]
                           if r["kind"] == kind and r["d"] == d)
        M = P.seed_pair(cfg32, kind, d)
        if kind in ("single", "mirror") or (kind, d) in (
                ("same", 12.0), ("same", 24.0)):
            dens_cache[(kind, d)] = my_edens(M, h)
        e = my_energy(M, h)
        c7["rows"].append({"box": "n32", "kind": kind, "d": d,
                           "E_claimed": claimed, "E_mine": e,
                           "rel_dev": abs(e - claimed) / abs(claimed)})
        print("claim7", json.dumps(c7["rows"][-1]), flush=True)
    # n48 spot check (one E_int point: pair + single)
    e48s = my_energy(P.seed_pair(cfg48, "single", 0.0), cfg48["h"])
    e48a = my_energy(P.seed_pair(cfg48, "anti", 18.0), cfg48["h"])
    c7["rows"].append({"box": "n48", "kind": "single", "d": 0.0,
                       "E_claimed": sl["n48"]["E_single"], "E_mine": e48s,
                       "rel_dev": abs(e48s - sl["n48"]["E_single"])
                       / sl["n48"]["E_single"]})
    cl = next(r["E"] for r in sl["n48"]["rows"]
              if r["kind"] == "anti" and r["d"] == 18.0)
    c7["rows"].append({"box": "n48", "kind": "anti", "d": 18.0,
                       "E_claimed": cl, "E_mine": e48a,
                       "rel_dev": abs(e48a - cl) / cl})
    print("claim7 n48", json.dumps(c7["rows"][-2:]), flush=True)
    # slope + linearity + ratio columns from the JSON (my arithmetic)
    for box in ("n32", "n48"):
        s_e = eint("same", box)
        dd = sorted(s_e)
        loc = [(s_e[dd[i + 1]] - s_e[dd[i]]) / (dd[i + 1] - dd[i])
               for i in range(len(dd) - 1)]
        A = np.vstack([dd, np.ones(len(dd))]).T
        coef, res, *_ = np.linalg.lstsq(A, np.array([s_e[d] for d in dd]),
                                        rcond=None)
        yv = np.array([s_e[d] for d in dd])
        r2 = 1.0 - float(res[0]) / float(((yv - yv.mean()) ** 2).sum())
        c7[f"same_slope_{box}"] = {
            "avg_slope": (s_e[dd[-1]] - s_e[dd[0]]) / (dd[-1] - dd[0]),
            "fit_slope": float(coef[0]), "local_slopes": loc,
            "R2": r2,
            "local_slope_spread": (max(loc) - min(loc)) / float(coef[0])}
    a48 = eint("anti", "n48")
    pred48 = 64.0 * np.pi * sl["n48"]["c2_seed"]
    ratios = {str(d): a48[d] / (-pred48 / d) for d in sorted(a48)}
    c7["anti_n48_ratio_mine"] = ratios
    rv = [ratios[k] for k in sorted(ratios, key=float)]
    c7["anti_n48_ratio_monotone"] = bool(
        all(rv[i + 1] < rv[i] for i in range(len(rv) - 1)))
    c7["anti_all_negative"] = bool(
        all(v < 0 for v in list(a32.values()) + list(a48.values())))
    c7["same_all_positive"] = bool(
        all(v > 0 for v in list(s32.values()) + list(eint("same",
                                                          "n48").values())))
    R["claim7"] = c7
    print("claim7 slopes", json.dumps({k: v for k, v in c7.items()
                                       if k != "rows"}, indent=1),
          flush=True)

    # ---- claim 6: Coulomb-form honesty + c2 self-cal ----------------
    c6 = {}
    dens32 = dens_cache[("single", 0.0)]
    rr32 = np.sqrt(X * X + Y * Y + Z * Z)
    sel = (rr32 > 9.0) & (rr32 < 0.5 * cfg32["L"] - 4 * h)
    c6["c2_seed_n32_mine"] = float(np.median(dens32[sel] * rr32[sel] ** 4)
                                   / 8.0)
    c6["c2_seed_n32_claimed"] = sl["n32"]["c2_seed"]
    for wlab, lo, hi_ in (("w10_16", 10.0, 16.0), ("w12_18", 12.0, 18.0)):
        s2 = (rr32 > lo) & (rr32 < hi_)
        c6[f"c2_n32_{wlab}"] = float(np.median(dens32[s2] * rr32[s2] ** 4)
                                     / 8.0)
    # density local exponent in the window (r^-4 form check)
    lr = np.log(rr32[sel])
    ld = np.log(np.maximum(dens32[sel], 1e-300))
    Aq = np.vstack([lr, np.ones_like(lr)]).T
    coef, *_ = np.linalg.lstsq(Aq, ld, rcond=None)
    c6["dens_exponent_n32_seed"] = float(coef[0])
    dens48 = my_edens(P.seed_pair(cfg48, "single", 0.0), cfg48["h"])
    X8, Y8, Z8 = INS.coords(48, cfg48["h"])
    rr48 = np.sqrt(X8 * X8 + Y8 * Y8 + Z8 * Z8)
    s8 = (rr48 > 9.0) & (rr48 < 0.5 * cfg48["L"] - 4 * cfg48["h"])
    c6["c2_seed_n48_mine"] = float(np.median(dens48[s8] * rr48[s8] ** 4)
                                   / 8.0)
    c6["c2_seed_n48_claimed"] = sl["n48"]["c2_seed"]
    lr = np.log(rr48[s8])
    ld = np.log(np.maximum(dens48[s8], 1e-300))
    Aq = np.vstack([lr, np.ones_like(lr)]).T
    coef, *_ = np.linalg.lstsq(Aq, ld, rcond=None)
    c6["dens_exponent_n48_seed"] = float(coef[0])
    # E_int local exponent, n48 anti (note admits ~1.5)
    dd = sorted(a48)
    c6["Eint_local_exponent_n48"] = [
        float(np.log(a48[dd[i]] / a48[dd[i + 1]])
              / np.log(dd[i + 1] / dd[i])) for i in range(len(dd) - 1)]
    # fixed-d box trend of the ratio (does box growth approach 1?)
    pred32 = 64.0 * np.pi * sl["n32"]["c2_seed"]
    c6["ratio_fixed_d"] = {
        str(d): {"n32": a32[d] / (-pred32 / d),
                 "n48": a48[d] / (-pred48 / d)}
        for d in (12.0, 18.0, 24.0)}
    R["claim6"] = c6
    print("claim6", json.dumps(c6, indent=1), flush=True)

    # ---- claim 8: string vs tube artifact ---------------------------
    c8 = {}
    ref32 = sl["n32"]["E_single"] + sl["n32"]["E_mirror"]
    rt_rows = {}
    for r_t in (2.0, 3.0, 4.5):
        es = {}
        for d in (12.0, 18.0, 24.0):
            es[d] = my_energy(my_seed(cfg32, "same", d, r_t=r_t), h) - ref32
        rt_rows[str(r_t)] = {
            "Eint": {str(d): es[d] for d in es},
            "avg_slope_12_24": (es[24.0] - es[12.0]) / 12.0,
            "slope_18_24": (es[24.0] - es[18.0]) / 6.0}
        print("claim8 r_t", r_t, json.dumps(rt_rows[str(r_t)]), flush=True)
    c8["rt_probe"] = rt_rows
    # regional decomposition of the same-sector E_int (default r_t)
    dens_ref = dens_cache[("single", 0.0)] + my_edens(
        P.seed_pair(cfg32, "mirror", 0.0), h)
    rho32 = np.hypot(X, Y)
    dec = {}
    for d in (12.0, 24.0):
        dp = dens_cache[("same", d)]
        dloc = (dp - dens_ref) * h ** 3
        tube = (rho32 < 6.0) & (np.abs(Z) < d / 2.0)
        cores = (np.sqrt(X ** 2 + Y ** 2 + (Z - d / 2) ** 2) < 6.0) | \
                (np.sqrt(X ** 2 + Y ** 2 + (Z + d / 2) ** 2) < 6.0)
        cores &= ~tube
        far = (rr32 > 18.0) & ~tube & ~cores
        rest = ~(tube | cores | far)
        dec[str(d)] = {"tube": float(dloc[tube].sum()),
                       "cores": float(dloc[cores].sum()),
                       "far": float(dloc[far].sum()),
                       "rest": float(dloc[rest].sum()),
                       "total": float(dloc.sum())}
    dec["slope_by_region_12_24"] = {
        k: (dec["24.0"][k] - dec["12.0"][k]) / 12.0
        for k in ("tube", "cores", "far", "rest", "total")}
    c8["region_decomposition"] = dec
    # r_c blend probe on the anti near zone (unclaimed hazard)
    rc_rows = {}
    for r_c in (2.5, 4.0):
        e_s = my_energy(my_seed(cfg32, "single", 0.0, r_c=r_c), h)
        e_m = my_energy(my_seed(cfg32, "mirror", 0.0, r_c=r_c), h)
        row = {}
        for d in (10.0, 12.0):
            ep = my_energy(my_seed(cfg32, "anti", d, r_c=r_c), h)
            row[str(d)] = ep - (e_s + e_m)
        rc_rows[str(r_c)] = row
        print("claim8 r_c", r_c, json.dumps(row), flush=True)
    c8["rc_probe_anti"] = rc_rows
    # singular-line probe: the same-seed inter-core axis is NOT smoothly
    # escaped; in-plane director magnitude as rho -> 0 (smooth escape
    # would require it -> 0)
    sing = {}
    for zc in (0.0, 4.0):
        row = {}
        for rho in (1e-3, 0.1, 0.5, 1.0):
            nv = nhat_same_exact(np.array([rho]), np.array([0.0]),
                                 np.array([zc]), 18.0)
            row[str(rho)] = float(np.hypot(nv[0, 0], nv[0, 1]))
        sing[f"z{zc:g}"] = row
    c8["same_seed_axis_inplane_mag"] = sing
    R["claim8"] = c8

    # ---- verdicts ---------------------------------------------------
    R["verdicts"] = [
        {"claim": "1 signed-charge instrument (+1/-1, antipair (+1,-1) "
                  "far 0, same far +2, conflict localization)",
         "claimed": "single +1.050, mirror -1.050, anti (+1.11, -1.12) "
                    "far 0.001, same far +2.078/+2.067, 20 conflicts at "
                    "the tube",
         "audit": "exact spheres: single +0.99998, mirror -0.99998, anti "
                  "per-core -/+0.99992 far 0.0 (1e-16), same far "
                  "magnitude 1.9998; my cubes 0.995-0.997; conflicts "
                  "0/0/0 and exactly 20, all at rho 0.75, |z|<=6.75",
         "verdict": "CONFIRMED",
         "note": "their cube +1.05 is quadrature overshoot (true +1.000); "
                 "the anti per-core SIGN assignment is lift-convention "
                 "dependent (my analytic-aligned lift reads top -1, bot "
                 "+1); same-sector per-core non-integers are a real field "
                 "property (singular line), exact spheres -0.067/-1.35"},
        {"claim": "2 E(d) 400-it table energies",
         "claimed": "anti d18 1.302179, anti d24 2.328807, same d18 "
                    "13.938505",
         "audit": "1.302179 / 2.328807 / 13.938505 (own functional), rel "
                  "dev <= 4.8e-8 (f32 storage floor); E_u splits match",
         "verdict": "CONFIRMED", "note": ""},
        {"claim": "3 attraction/repulsion signs + core stationarity",
         "claimed": "antipair BOUND (E_int<0 all d, attraction); same = "
                    "string not repulsion; no static separated-antipair "
                    "regime (heal protocol failed)",
         "audit": "anti dE/dd>0 at all 6 d (attraction as d shrinks); "
                  "same dE/dd>0 too (net inward, NOT repulsion); my "
                  "deviator core scan: seed dips at +/-9 present, it400 "
                  "dips ABSENT (min_dev 0.43-0.51 vs bulk 0.726), it120 "
                  "d24 mid-annihilation; anti it400 far flux -4e-5, 0 "
                  "conflicts",
         "verdict": "CONFIRMED",
         "note": "the original DoD stationarity criterion is indeed "
                 "unachievable; the annihilation-scan re-label is "
                 "honest and necessary"},
        {"claim": "6 Coulomb-form comparison honesty + c2 self-cal",
         "claimed": "c2_seed 0.5309 (n32) / 0.5251 (n48); form test "
                    "within stated caveats; ratio approaching 1 as d and "
                    "box grow",
         "audit": "c2 reproduced exactly, window-stable +/-0.5%, density "
                  "exponent -4.02 (both boxes); E_int local exponent "
                  "1.66/1.56/1.38 (matches admitted ~1.5); BUT at fixed "
                  "d the ratio moves AWAY from 1 with box growth (d18 "
                  "1.01->1.57, d24 0.68->1.34)",
         "verdict": "PARTIAL",
         "note": "caveats are stated and repeated; the one overreach is "
                 "'approaching ... as the box grows': only the within-"
                 "n48 d-trend approaches; the box direction at fixed d "
                 "does not"},
        {"claim": "7 seed-level ladder E_int(d)",
         "claimed": "E_single 18.8316/21.2630, anti negative all d, same "
                    "positive, slope ~7.0/6.2, approx linear, n48 ratio "
                    "2.05->1.22 monotone",
         "audit": "8 states re-computed with own functional: rel dev <= "
                  "2e-16; signs confirmed; avg slopes 7.02/6.21; ratios "
                  "2.049->1.226 monotone; R2 0.988/0.989 BUT local slope "
                  "drifts 8.9->5.2 (n32), 8.0->5.2 (n48)",
         "verdict": "CONFIRMED",
         "note": "'approximately linear' carries a systematic 44-52% "
                 "local-slope drift (concave); the average-slope quote "
                 "is exact"},
        {"claim": "8 string-confinement interpretation",
         "claimed": "same-sector linear E(d) = physical string (tension "
                    "~7.0/6.2) required topologically; smooth +z escape "
                    "tube",
         "audit": "linear growth lives 100% in the tube region (region "
                  "slope: tube 6.95 of total 6.71) and survives r_t "
                  "2->4.5; BUT tension varies 16% with r_t (7.28/6.71/"
                  "6.22), 12% across boxes, 40-50% across the d-window; "
                  "AND the seed axis is NOT smoothly escaped: exact "
                  "in-plane magnitude 0.707 (z=0) / 0.36 (z=4) as "
                  "rho->0, a genuine singular winding-2 line",
         "verdict": "PARTIAL",
         "note": "string FORM robust (structure is topologically "
                 "required; k=2 transverse winding must cross every "
                 "inter-core plane); the note's 'smooth escape blend "
                 "(integer winding escapes)' misdescribes its own seed, "
                 "the tension value is ansatz-grade and cutoff-"
                 "sensitive, and the axis argument is a lift-level "
                 "argument (apolar +/-z are one state) without saying "
                 "so"},
    ]
    R["unclaimed_hazards"] = [
        "SINGULAR LINE: the same-pair seed's escape blend does not "
        "desingularize the inter-core axis (in-plane winding-2 component "
        "magnitude 0.36-0.71 persists at rho->0 along the whole "
        "segment). Seed-level string tension therefore includes lattice-"
        "regularized singular-line energy: h-dependent, so the 7.0/6.2 "
        "number is cutoff-sensitive (no h-refinement was run on it).",
        "BLEND ARTIFACT IN THE ANTI NEAR ZONE: E_int at small d moves "
        "21-26% with the construction blend radius r_c (d10: -19.88 at "
        "r_c=4 vs -15.64 at r_c=2.5; d12: -14.66 vs -10.86): part of "
        "the 'near-zone enhancement' (the largest table ratios) is "
        "construction, not physics.",
        "SAME-SECTOR OFFSET: the charge-2 far texture vs 2-singles "
        "reference contributes a d-independent-ish +13 (far region) to "
        "E_int; harmless for slopes, box-function for magnitudes.",
        "CONVENTION FRAGILITY: per-state global lift sign is arbitrary; "
        "the gates table prints absolute per-core signs without the "
        "convention attached (my analytic-aligned lift reads the "
        "antipair as (-1, +1) top/bot).",
    ]
    R["meta"]["wall_s"] = time.time() - t00
    with open(os.path.join(DATA, "m5_21_4_audit.json"), "w") as f:
        json.dump(R, f, indent=1, default=float)
    print(f"audit done {R['meta']['wall_s']:.0f}s -> m5_21_4_audit.json",
          flush=True)


# =====================================================================
# ================= ROUND 2: the capture-run claims ===================
# =====================================================================
DEC = P.DEC
SCRATCH = os.environ.get(
    "AUDIT_SCRATCH",
    "/private/tmp/claude-501/-Users-xrodz-Documents-source-code-NEPTUNYA-"
    "SABER/2a442551-3079-4299-8564-9369478bed29/scratchpad")


def my_leap_step_tracked(M, V, cfg, free, gam, dt):
    """bit-parity replica of DEC.leap_step that ALSO returns the exact
    kinetic energy removed by the two implicit damping divisions
    (V -> V/(1+g), g = 0.5 dt gamma(r)); per division the removal is
    0.5 h^3 (|V_pre|^2 - |V_post|^2) = h^3 sum g(1+g/2)|V_post|^2."""
    h3 = cfg["h"] ** 3
    w = 0.0
    F = -INS.grad(M, cfg) / h3
    V = (V + 0.5 * dt * F) * free
    ke = 0.5 * h3 * float(np.sum(V * V))
    V = V / (1.0 + 0.5 * dt * gam)
    w += ke - 0.5 * h3 * float(np.sum(V * V))
    M = M + dt * V
    F = -INS.grad(M, cfg) / h3
    V = (V + 0.5 * dt * F) * free
    ke = 0.5 * h3 * float(np.sum(V * V))
    V = V / (1.0 + 0.5 * dt * gam)
    w += ke - 0.5 * h3 * float(np.sum(V * V))
    return M, V, w


def r2_1_closure(steps=400, n=32, d0=18.0, heal_it=30, gam_in=0.02,
                 dt=0.025):
    """independent ledger closure: by-difference absorbed vs the exact
    discrete damping work vs the naive Riemann work integral, on a short
    damped antipair evolution (same seed builder, same integrator
    scheme, own tracking)."""
    cfg = INS.base_cfg(term="T2", n=n, L=1.5 * n, bc="free",
                       stencil="sym")
    free_b = np.ones((n, n, n), dtype=bool)
    M0 = P.seed_pair(cfg, "anti", d0)
    M, _, _ = INS.fire(M0, cfg, free_b, heal_it, log_every=heal_it,
                       tag="r2heal")
    fr = free_b[..., None, None].astype(float)
    gam = DEC.sponge(cfg) + gam_in
    out = {"n": n, "d0": d0, "heal_it": heal_it, "gam_in": gam_in,
           "dt": dt, "steps": steps}
    # step parity vs DEC.leap_step (must be bit-exact)
    Ma, Va = M.copy(), np.zeros_like(M)
    Mb, Vb = M.copy(), np.zeros_like(M)
    dev = 0.0
    for _ in range(3):
        Ma, Va, _ = my_leap_step_tracked(Ma, Va, cfg, fr, gam, dt)
        Mb, Vb = DEC.leap_step(Mb, Vb, cfg, fr, gam, dt)
        dev = max(dev, float(np.abs(Ma - Mb).max()),
                  float(np.abs(Va - Vb).max()))
    out["step_parity_maxdev"] = dev
    # tracked run
    V = np.zeros_like(M)
    E0, K0 = DEC.e_tot(M, V, cfg)
    myE0 = my_energy(M, cfg["h"])
    out["E0_dec"], out["E0_mine"] = E0, myE0
    absorbed = w_exact = w_riem = 0.0
    Eb, Kb = E0, K0
    marks = {}
    t0 = time.time()
    for it in range(1, steps + 1):
        M, V, w = my_leap_step_tracked(M, V, cfg, fr, gam, dt)
        Ea, Ka = DEC.e_tot(M, V, cfg)
        absorbed += (Eb + Kb) - (Ea + Ka)
        w_exact += w
        w_riem += dt * cfg["h"] ** 3 * float(np.sum(gam * V * V))
        Eb, Kb = Ea, Ka
        if it in (200, steps):
            marks[str(it)] = {
                "absorbed_by_diff": absorbed,
                "W_exact_discrete": w_exact,
                "W_riemann_endV": w_riem,
                "gap_exact_rel": (absorbed - w_exact)
                / max(abs(absorbed), 1e-300),
                "gap_riemann_rel": (absorbed - w_riem)
                / max(abs(absorbed), 1e-300),
                "E_plus_KE": Ea + Ka}
            print(f"r2_1 it {it}: " + json.dumps(marks[str(it)]),
                  flush=True)
    out["marks"] = marks
    out["E_end_mine"] = my_energy(M, cfg["h"])
    out["E_drop_mine"] = myE0 - out["E_end_mine"]
    out["E_drop_dec"] = E0 - Eb
    # gamma = 0 drift control from the same healed start
    Mg, Vg = Ma.copy(), np.zeros_like(M)   # reuse post-parity state? no:
    Mg, _, _ = INS.fire(P.seed_pair(cfg, "anti", d0), cfg, free_b,
                        heal_it, log_every=heal_it, tag="r2heal0")
    Vg = np.zeros_like(Mg)
    Eg0, Kg0 = DEC.e_tot(Mg, Vg, cfg)
    for _ in range(200):
        Mg, Vg = DEC.leap_step(Mg, Vg, cfg, fr, 0.0, dt)
    Eg1, Kg1 = DEC.e_tot(Mg, Vg, cfg)
    out["gamma0_drift_abs_200"] = abs((Eg1 + Kg1) - (Eg0 + Kg0))
    out["gamma0_drift_rel_200"] = out["gamma0_drift_abs_200"] \
        / max(abs(Eg0 + Kg0), 1e-300)
    out["wall_s"] = time.time() - t0
    return out


def r2_1_rows_circularity():
    """structural check ON THE DATA: the cap rows' absorbed column is
    exactly (E0+KE0) - (E+KE) at every row (pure telescoping), and the
    ledger column is identically zero."""
    j = jload("m5_21_4_ev_cap_rows.json")
    tot0 = j["E0"] + j["KE0"]
    led = max(abs(r["ledger"]) for r in j["rows"])
    dev = max(abs(r["absorbed"] - (tot0 - (r["E"] + r["KE"])))
              for r in j["rows"])
    return {"n_rows": len(j["rows"]), "ledger_absmax": led,
            "absorbed_minus_telescope_absmax": dev}


# ---------- snapshot instruments (own) ----------
def cube_flux_faces(B, h, center, half):
    """my_cube_flux split into its 6 face contributions."""
    n = B.shape[0]
    c = (n - 1) / 2.0
    ic = [int(round(center[a] / h + c)) for a in range(3)]
    k = int(round(half / h))
    lo = [ic[a] - k for a in range(3)]
    hi = [ic[a] + k for a in range(3)]
    if min(lo) < 1 or max(hi) > n - 2:
        return None
    m = 2 * k + 1
    w1 = np.ones(m)
    w1[0] = w1[-1] = 0.5
    W = np.outer(w1, w1)
    out = {}
    labs = (("xm", "xp"), ("ym", "yp"), ("zm", "zp"))
    for ax in range(3):
        sl_hi = [slice(lo[a], hi[a] + 1) for a in range(3)] + [ax]
        sl_lo = [slice(lo[a], hi[a] + 1) for a in range(3)] + [ax]
        sl_hi[ax], sl_lo[ax] = hi[ax], lo[ax]
        out[labs[ax][1]] = float((B[tuple(sl_hi)] * W).sum()) \
            * h * h / (4.0 * np.pi)
        out[labs[ax][0]] = -float((B[tuple(sl_lo)] * W).sum()) \
            * h * h / (4.0 * np.pi)
    return out


def plane_profile(B, h, rho2d, rho_min=None, rho_max=None):
    """F(z) = (h^2/4pi) sum_{x,y} B_z, optionally rho-restricted."""
    m = np.ones(rho2d.shape, dtype=bool)
    if rho_min is not None:
        m &= rho2d >= rho_min
    if rho_max is not None:
        m &= rho2d < rho_max
    return (B[..., 2] * m[..., None]).sum(axis=(0, 1)) \
        * h * h / (4.0 * np.pi)


def slab_lobes(zax, F):
    """slab charge density dF between adjacent planes; per half-space
    the dominant-sign lobe's integral, |d|-weighted centroid, peak z."""
    zc = 0.5 * (zax[1:] + zax[:-1])
    dF = np.diff(F)
    out = {}
    for lab, sel in (("top", zc > 1.0), ("bot", zc < -1.0)):
        d, zz = dF[sel], zc[sel]
        q = float(d.sum())
        s = 1.0 if q >= 0 else -1.0
        dp = np.clip(s * d, 0.0, None)
        ws = float(dp.sum())
        out[lab] = {
            "q_half": q,
            "centroid": float((zz * dp).sum() / ws) if ws > 1e-12
            else None,
            "z_peak": float(zz[int(np.argmax(dp))]) if ws > 1e-12
            else None,
            "q_domsign": s * ws}
    return out


def gap_z_profile(M, rho2d, rho_max=3.2):
    lam = np.linalg.eigvalsh(M)
    gap = np.minimum(lam[..., 2] - lam[..., 1],
                     lam[..., 1] - lam[..., 0])
    m = rho2d < rho_max
    return np.where(m[..., None], gap, np.inf).min(axis=(0, 1))


def r2_snapshot_pass(npz_name="m5_21_4_ev_cap.npz", d0=24.0):
    cfg = INS.base_cfg(term="T2", n=48, L=72.0, bc="free",
                       stencil="sym")
    n, h = cfg["n"], cfg["h"]
    X, Y, Z3 = INS.coords(n, h)
    rr = np.sqrt(X * X + Y * Y + Z3 * Z3)
    rho2d = np.hypot(X[:, :, 0], Y[:, :, 0])
    rho3d = np.hypot(X, Y)
    zax = (np.arange(n) - (n - 1) / 2.0) * h
    far_m = rr > 20.0
    z = np.load(os.path.join(DATA, npz_name))
    keys = sorted(z.files, key=lambda k: int(k[len("M_it"):]))
    zcs = np.arange(-18.0, 18.1, 1.5)
    rows = []
    for key in keys:
        it = int(key[len("M_it"):])
        t0 = time.time()
        M = z[key].astype(np.float64)
        vo, confs = my_orient(M)
        # sign anchor: analytic anti far field is (0,0,-1); one
        # convention for every snapshot so signs are time-comparable
        if float(vo[far_m][:, 2].sum()) > 0.0:
            vo = -vo
        B = my_B(vo, h)
        row = {"it": it, "t": it * 0.025, "conflicts": len(confs)}
        # R2-2 fixed cubes (two sizes) + spheres + far
        row["Q_top_c8"] = my_cube_flux(B, h, (0, 0, +d0 / 2), 8.0)
        row["Q_bot_c8"] = my_cube_flux(B, h, (0, 0, -d0 / 2), 8.0)
        row["Q_top_c6"] = my_cube_flux(B, h, (0, 0, +d0 / 2), 6.0)
        row["Q_bot_c6"] = my_cube_flux(B, h, (0, 0, -d0 / 2), 6.0)
        row["Q_top_sR8"] = sphere_degree_lattice(vo, h, (0, 0, d0 / 2),
                                                 8.0, nth=90, nph=180)
        row["Q_bot_sR8"] = sphere_degree_lattice(vo, h, (0, 0, -d0 / 2),
                                                 8.0, nth=90, nph=180)
        row["Q_far_c27"] = my_cube_flux(B, h, (0, 0, 0), 27.0)
        row["Q_far_c30"] = my_cube_flux(B, h, (0, 0, 0), 30.0)
        row["far_faces_c27"] = cube_flux_faces(B, h, (0, 0, 0), 27.0)
        # R2-3/R2-5 plane-flux profiles (full + rho-resolved)
        F = plane_profile(B, h, rho2d)
        row["F"] = [float(f) for f in F]
        row["F_ax6"] = [float(f) for f in
                        plane_profile(B, h, rho2d, rho_max=6.0)]
        row["F_ann6_15"] = [float(f) for f in
                            plane_profile(B, h, rho2d, 6.0, 15.0)]
        row["F_out15"] = [float(f) for f in
                          plane_profile(B, h, rho2d, rho_min=15.0)]
        inner = slice(1, -1)
        row["F_stats"] = {
            "min": float(F[inner].min()), "max": float(F[inner].max()),
            "z_at_min": float(zax[inner][int(np.argmin(F[inner]))]),
            "z_at_max": float(zax[inner][int(np.argmax(F[inner]))])}
        row["slab"] = slab_lobes(zax[inner], F[inner])
        # moving-window flux scan
        scan = [my_cube_flux(B, h, (0, 0, zc), 4.5) for zc in zcs]
        row["scan_zc"] = [float(v) for v in zcs]
        row["scan_Q"] = [float(v) for v in scan]
        for lab, sel in (("top", zcs > 2.0), ("bot", zcs < -2.0)):
            qs = np.array(scan)[sel]
            zz = zcs[sel]
            i = int(np.argmax(np.abs(qs)))
            row[f"scan_{lab}_absmax"] = {"zc": float(zz[i]),
                                         "Q": float(qs[i])}
        # gap + energy localization
        g = gap_z_profile(M, rho2d)
        for lab, sel in (("top", zax > 2.0), ("bot", zax < -2.0)):
            gg = g[sel]
            i = int(np.argmin(gg))
            row[f"gapmin_{lab}"] = {"z": float(zax[sel][i]),
                                    "gap": float(gg[i])}
        dens = my_edens(M, h)
        e_ax = (dens * (rho3d < 4.5)).sum(axis=(0, 1)) * h ** 3
        row["e_ax"] = [float(v) for v in e_ax]
        for lab, sel in (("top", zax > 2.0), ("bot", zax < -2.0)):
            i = int(np.argmax(e_ax[sel]))
            row[f"e_ax_peak_{lab}"] = {"z": float(zax[sel][i]),
                                       "e": float(e_ax[sel][i])}
        row["E_r_split"] = {
            "inner_r12": float(dens[rr < 12.0].sum() * h ** 3),
            "mid": float(dens[(rr >= 12.0) & (rr < 23.4)].sum()
                         * h ** 3),
            "sponge_r23.4": float(dens[rr >= 23.4].sum() * h ** 3)}
        rows.append(row)
        print(f"r2 snap it{it}: Qtop {row['Q_top_c8']:+.4f} "
              f"Qbot {row['Q_bot_c8']:+.4f} far {row['Q_far_c30']:+.1e} "
              f"F[{row['F_stats']['min']:+.4f},"
              f"{row['F_stats']['max']:+.4f}] "
              f"[{time.time() - t0:.0f}s]", flush=True)
    return rows


def r2_5_floor_probe(d0=24.0):
    """f32 floor of the plane-flux chain: analytic antipair seed at
    f64 through lift+B+profile vs the same state quantized to f32;
    plus lift-choice sensitivity (my BFS vs the task's orient_v1) on
    the final snapshot."""
    cfg = INS.base_cfg(term="T2", n=48, L=72.0, bc="free",
                       stencil="sym")
    n, h = cfg["n"], cfg["h"]
    X, Y, Z3 = INS.coords(n, h)
    rr = np.sqrt(X * X + Y * Y + Z3 * Z3)
    rho2d = np.hypot(X[:, :, 0], Y[:, :, 0])
    far_m = rr > 20.0

    def chain(M):
        vo, _ = my_orient(M)
        if float(vo[far_m][:, 2].sum()) > 0.0:
            vo = -vo
        return plane_profile(my_B(vo, h), h, rho2d)

    M64 = P.seed_pair(cfg, "anti", d0)
    F64 = chain(M64)
    F32 = chain(M64.astype(np.float32).astype(np.float64))
    out = {"f32_floor_maxdF": float(np.abs(F64 - F32)[1:-1].max()),
           "F_seed64": [float(v) for v in F64]}
    # lift-choice sensitivity at it4000
    z = np.load(os.path.join(DATA, "m5_21_4_ev_cap.npz"))
    M = z["M_it4000"].astype(np.float64)
    F_bfs = chain(M)
    vo2, _ = P.orient_v1(M)
    if float(vo2[far_m][:, 2].sum()) > 0.0:
        vo2 = -vo2
    F_v1 = plane_profile(my_B(vo2, h), h, rho2d)
    out["lift_choice_maxdF_it4000"] = float(
        np.abs(F_bfs - F_v1)[1:-1].max())
    out["F_it4000_bfs_minmax"] = [float(F_bfs[1:-1].min()),
                                  float(F_bfs[1:-1].max())]
    out["F_it4000_v1_minmax"] = [float(F_v1[1:-1].min()),
                                 float(F_v1[1:-1].max())]
    return out


def main_round2():
    t00 = time.time()
    os.makedirs(SCRATCH, exist_ok=True)
    res = {"date": "2026-07-21"}
    res["r2_1_rows_circularity"] = r2_1_rows_circularity()
    print("circularity", json.dumps(res["r2_1_rows_circularity"]),
          flush=True)
    res["r2_1"] = r2_1_closure()
    res["r2_5_floor"] = r2_5_floor_probe()
    print("floor", json.dumps({k: v for k, v in
                               res["r2_5_floor"].items()
                               if not k.startswith("F_seed")}),
          flush=True)
    res["snapshots"] = r2_snapshot_pass()
    res["wall_s"] = time.time() - t00
    with open(os.path.join(SCRATCH, "m5_21_4_r2_results.json"),
              "w") as f:
        json.dump(res, f, indent=1, default=float)
    print(f"round2 compute done {res['wall_s']:.0f}s", flush=True)


def main_merge2():
    """fold the round-2 computed results + verdicts into
    ../data/m5_21_4_audit.json under the 'round2' key."""
    with open(os.path.join(SCRATCH, "m5_21_4_r2_results.json")) as f:
        res = json.load(f)
    # brief item: snapshot E re-read vs the live rows (own functional)
    rows_j = jload("m5_21_4_ev_cap_rows.json")
    z = np.load(os.path.join(DATA, "m5_21_4_ev_cap.npz"))
    h48 = 1.5
    ere = {}
    for it in (2000, 4000):
        e_mine = my_energy(z[f"M_it{it}"].astype(np.float64), h48)
        e_rows = next(r["E"] for r in rows_j["rows"] if r["it"] == it)
        ere[str(it)] = {"E_rows": e_rows, "E_mine": e_mine,
                        "rel_dev": abs(e_mine - e_rows) / abs(e_rows)}
    # KE-shape check on the rows (note section 5 ledger row)
    t = np.array([r["t"] for r in rows_j["rows"]])
    ke = np.array([r["KE"] for r in rows_j["rows"]])
    sel = t <= 15.0
    ke_exp = float(np.polyfit(np.log(t[sel]), np.log(ke[sel]), 1)[0])
    ipk = int(np.argmax(ke))
    # compact per-snapshot summary (full arrays stay in the scratch file)
    snaps = []
    for s in res["snapshots"]:
        F = np.array(s["F"])[1:-1]
        snaps.append({k: s[k] for k in
                      ("it", "t", "conflicts", "Q_top_c8", "Q_bot_c8",
                       "Q_top_c6", "Q_bot_c6", "Q_top_sR8", "Q_bot_sR8",
                       "Q_far_c27", "Q_far_c30", "far_faces_c27",
                       "F_stats", "slab", "scan_top_absmax",
                       "scan_bot_absmax", "gapmin_top", "gapmin_bot",
                       "e_ax_peak_top", "e_ax_peak_bot", "E_r_split")}
                     | {"F_mean": float(F.mean())})
    zax = (np.arange(48) - 23.5) * 1.5
    sel_out = (np.abs(zax) > 16) & (np.abs(zax) < 34)
    F0 = np.array(next(s["F"] for s in res["snapshots"]
                       if s["it"] == 0))
    F4k = np.array(next(s["F"] for s in res["snapshots"]
                        if s["it"] == 4000))
    Fs = np.array(res["r2_5_floor"]["F_seed64"])
    bg = {"seed_outside_minmaxmean":
          [float(Fs[sel_out].min()), float(Fs[sel_out].max()),
           float(Fs[sel_out].mean())],
          "it0_outside_minmaxmean":
          [float(F0[sel_out].min()), float(F0[sel_out].max()),
           float(F0[sel_out].mean())],
          "it4000_outside_minmaxmean":
          [float(F4k[sel_out].min()), float(F4k[sel_out].max()),
           float(F4k[sel_out].mean())],
          "it4000_between_minmax":
          [float(F4k[np.abs(zax) < 10].min()),
           float(F4k[np.abs(zax) < 10].max())]}
    R2 = {
        "meta": {"date": res["date"], "wall_s": res["wall_s"],
                 "scope": ["R2-1", "R2-2", "R2-3", "R2-4", "R2-5"],
                 "rerun": "n=32 d0=18 heal=30 free bc, 400 damped steps",
                 "full_arrays": "scratch m5_21_4_r2_results.json "
                                "(session-local)"},
        "r2_1_rows_circularity": res["r2_1_rows_circularity"],
        "r2_1_closure": {k: res["r2_1"][k] for k in
                         ("step_parity_maxdev", "marks",
                          "gamma0_drift_rel_200", "E0_dec", "E0_mine")},
        "snapshot_E_reread": ere,
        "rows_KE_shape": {"loglog_exp_t2.5_15": ke_exp,
                          "KE_peak": {"t": float(t[ipk]),
                                      "KE": float(ke[ipk])}},
        "snapshots": snaps,
        "residual_background": bg,
        "r2_5_floor": {k: v for k, v in res["r2_5_floor"].items()
                       if not k.startswith("F_seed")},
    }
    R2["verdicts"] = [
        {"claim": "R2-1 ledger circularity self-catch + closure fix "
                  "sufficiency",
         "claimed": "rows 'absorbed' is by-difference (ledger column "
                    "trivially zero); independent closure = gamma=0 "
                    "gate + audit work integral",
         "audit": "circularity PROVEN on the data: ledger == 0.0 and "
                  "absorbed == (E0+KE0)-(E+KE) EXACTLY at all 40 rows. "
                  "Independent closure on my own tracked rerun (bit-"
                  "exact step replica, parity dev 0.0): absorbed-by-"
                  "difference matches the exact discrete damping work "
                  "h^3 sum g(1+g/2)|V_post|^2 (two implicit divisions "
                  "per step) to 6.9e-4 rel at 400 steps (1.5e-3 at "
                  "200), and the naive Riemann integral dt h^3 sum "
                  "gamma|V|^2 to 1.3e-3 (3.6e-3); gamma=0 drift "
                  "control 5.1e-6 rel/200 steps (their gate: 4.6e-6)",
         "verdict": "CONFIRMED",
         "note": "'absorbed' is physically the damping/sponge removal "
                 "to ~0.1%, not integrator drift; the honesty label + "
                 "closure route in the note are sufficient"},
        {"claim": "R2-2 charge dissolution sigmoid +/-0.48 -> +/-0.02, "
                  "far <= 1e-4",
         "claimed": "fixed cubes at the seed core positions: sigmoid "
                    "+/-0.48 (post-heal t=0) -> +/-0.02 by t=100, main "
                    "drop t 20-70; far-sphere charge <= 1e-4 "
                    "throughout",
         "audit": "my own lift (BFS, analytic-anchored) + trapezoid "
                  "cubes on all 11 snapshots: |Q| 0.408/0.437 (t=0) -> "
                  "0.0034/0.0030 (t=100), smooth monotone, half-loss "
                  "t ~ 43-47, 90% gone by t ~ 62-65; the t=0 offset vs "
                  "their 0.48 is the known 5-16% cube-quadrature "
                  "spread; endpoint is CLEANER than claimed (<=0.005 "
                  "on cubes at two sizes AND spheres vs their 0.02). "
                  "Far: their instrument <= 1.11e-4 at all 40 rows "
                  "(reproduced arithmetic); MY far cube reads up to "
                  "7.3e-3 while the pair lives, traced ENTIRELY to the "
                  "two z-face punctures of the pre-existing axial line "
                  "(constant +0.052/-0.040 through it2000, then "
                  "<= 5.1e-5)",
         "verdict": "CONFIRMED",
         "note": "the sigmoid + interior-zero content confirmed and "
                 "strengthened; the specific '<= 1e-4' bound is "
                 "instrument-fragile (hazard list)"},
        {"claim": "R2-3 conduit annihilation, no ballistic core-walk",
         "claimed": "the +/- steps dissolve IN PLACE, the inter-core "
                    "line flux unwinds to zero by t ~ 50-60; no "
                    "localized core translation toward the center",
         "audit": "three independent localizers agree: (a) deep-gap "
                  "core locus stays |z| 15-17 through the entire "
                  "charge collapse, never migrates to 0; (b) slab-"
                  "charge centroid moves only 10.5 -> 9.7 (top) / "
                  "-9.2 -> -6.3 (bot) while amplitude falls 0.39 -> "
                  "0.18: >50% of the charge cancels with the centroids "
                  "still at ~80% of their initial separation, and ~68% "
                  "cancels before any centroid passes |z| = 6; (c) "
                  "moving-window flux peak fixed at |zc| = 9 through "
                  "t = 50 with amplitude 0.32 -> 0.15. The between-"
                  "core plateau 0.364 -> 0.094 by t = 60, gone at "
                  "t = 70 (matches 'unwinds by t ~ 50-60'). The walk-"
                  "then-merge alternative REQUIRES amplitude survival "
                  "until the lobes meet: refuted. Refinement: the LAST "
                  "~30% of the residual density does slide inward "
                  "(centroid -> |z| 3.5-4.4) during t = 50-70 while "
                  "dying",
         "verdict": "CONFIRMED",
         "note": "conduit annihilation is what the snapshots support; "
                 "the late-remnant inward slide is a refinement, not a "
                 "walk (amplitude already down 70+%)"},
        {"claim": "R2-4 interior annihilation vs boundary drain",
         "claimed": "far charge never moves, boundary planes carry "
                    "only the weak residual line flux, energy exits as "
                    "absorbed radiation; drain signature absent",
         "audit": "far-cube FACE histories flat (+0.052/-0.040 z-face "
                  "line punctures, x/y faces ~1e-3) through it2000, "
                  "collapsing to <= 5.1e-5 after annihilation: no "
                  "transit event on any face; the every-100-step rows "
                  "(2.5 t-units) bound |far| <= 1.11e-4 continuously, "
                  "closing the 400-step snapshot gap (a core exit "
                  "would swing a face by ~0.4 for many consecutive "
                  "rows at the observed dynamics speed ~0.24 units/t); "
                  "radial energy split declines in ALL shells (sponge "
                  "0.865 -> 0.634, no late pile-up). Discriminator "
                  "sound",
         "verdict": "CONFIRMED",
         "note": "interior annihilation; the absorbed channel (31%) is "
                 "consistent with the R2-1 work-integral closure"},
        {"claim": "R2-5 residual all-positive axial line flux "
                  "0.03-0.08 per plane at t=100: line debris, not "
                  "vacuum",
         "claimed": "a weak all-positive axial line flux persists at "
                    "t = 100; endpoint = line debris, not clean "
                    "vacuum",
         "audit": "reality: CONFIRMED far above floor (f32 chain floor "
                  "6.6e-9, signal 7 orders larger; two independent "
                  "lifts give bit-identical profiles at it4000); "
                  "uniform sign + range confirmed (mine 0.031-0.080 "
                  "per plane, all one sign in the anchored lift; sign "
                  "itself is lift-convention). NOT VACUUM: confirmed. "
                  "'DEBRIS' REFUTED: the flux is seed-inherited "
                  "background, not an annihilation product: outside-"
                  "core planes carry the SAME flux at t=0 (healed "
                  "mean -0.038; raw analytic seed -0.045) as at t=100 "
                  "(-0.046); the alpha-map far tail threads the whole "
                  "box from the start. The annihilation adds only a "
                  "modest deepening on the formerly-between-cores "
                  "planes (-0.076..-0.080 vs -0.046 background). And "
                  "within the window it does NOT anneal: net per-"
                  "plane flux still growing at t=100 (mean -0.0591 -> "
                  "-0.0599 over the last 10 t-units) while spreading "
                  "radially (rho<6 channel -0.094 -> -0.072, annulus "
                  "6-15 deepens -0.058 -> -0.086)",
         "verdict": "PARTIAL",
         "note": "'persistent seed-inherited background + weak line "
                 "remnant, radially spreading, non-annealing in this "
                 "window' is what the data shows; 'debris' (transient "
                 "annihilation byproduct) is not supported"},
    ]
    R2["unclaimed_hazards"] = [
        "KE SHAPE: the note's 'KE grows ~ t^2 early (constant-force "
        "phase)' is not what the rows show: local log-log exponent "
        "1.45 over t = 2.5-15 (1.56 at the earliest row pair, falling); "
        "any t^2 regime lives below the first recorded row t = 2.5. "
        "Peak: 0.446 at the t = 80 row (note says 0.45 at t ~ 77: row-"
        "cadence rounding, fine).",
        "FAR-BOUND FRAGILITY: '|Q_far| <= 1e-4' is specific to the "
        "task's cube quadrature; an independent trapezoid cube reads "
        "up to 7.3e-3 while the pair lives, entirely the two z-face "
        "punctures of the box-threading axial line. Any 'far charge' "
        "bound read through a surface the line pierces is a line-flux "
        "mixture, not a clean charge bound.",
        "SEED FAR-FIELD BACKGROUND: the antipair alpha-map's far tail "
        "leaves a uniform-sign per-plane flux -0.03..-0.06 threading "
        "the ENTIRE box already at t = 0 (raw seed AND healed): the "
        "t = 100 'residual line' is mostly this pre-existing "
        "background. Future antipair arenas wanting a clean vacuum "
        "endpoint need a compactified far field (or the background "
        "subtracted).",
        "ENDPOINT CHARGE: the claimed +/-0.02 endpoint is the task "
        "instrument's floor in the presence of the line; my cubes "
        "(two sizes) and spheres read <= 0.005 at t = 100: the "
        "annihilation endpoint is cleaner than claimed.",
        "LINE-CONTAMINATED PER-CORE READS: while the conduit is "
        "active, ALL closed surfaces around one core are pierced "
        "twice by the axial line: my R = 8 spheres read +0.32 -> "
        "+0.52 (RISING while the cube reads fall) on the same states. "
        "The per-core 'charge' during the capture is a core+line flux "
        "mixture on every instrument (the task's +/-0.48 at t = 0 "
        "included); only the far-zone additivity and the endpoint are "
        "surface-shape independent.",
    ]
    ap = os.path.join(DATA, "m5_21_4_audit.json")
    with open(ap) as f:
        R = json.load(f)
    R["meta"]["skipped"]["4"] = "resolved in round 2"
    R["meta"]["skipped"]["5"] = "resolved in round 2"
    R["round2"] = R2
    with open(ap, "w") as f:
        json.dump(R, f, indent=1, default=float)
    print("merged round2 ->", ap)
    print(json.dumps({"snapshot_E_reread": ere,
                      "rows_KE_shape": R2["rows_KE_shape"]}, indent=1))
    for v in R2["verdicts"]:
        print(f"{v['claim'][:60]:60s}  {v['verdict']}")


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "round1"
    if mode == "round2":
        main_round2()
    elif mode == "round1":
        main()
    elif mode == "merge2":
        main_merge2()
    else:
        raise SystemExit(f"unknown mode {mode}")
