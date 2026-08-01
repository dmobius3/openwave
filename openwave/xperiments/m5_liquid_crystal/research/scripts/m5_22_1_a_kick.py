"""M5.22.1 opening move: the kick-apart identity probe on the two-ring
neutral census state (the author's dineutron test, 2026-07-30 reply).

Target: the m5_22 census endpoint P-1 (Q = 0, central column + two
vortex rings at Z ~ +-12.5). The author's ask: kick the two rings
apart; expected outcome two baryons (ring count = baryon number).

Two machine-checkable branches, both run (series exhaustion rule):

  split  static separation: resample M along z with the odd shift map
         z_src = z - d*tanh(z/w0) (w0 >= 1.5 d keeps it monotone), the
         pin shell re-imposed exactly, then FIRE-relax. Outcomes:
         return (one bound object) / split (separated structures) /
         dissolve.
  kick   dynamic separation: velocity kick V0 = -v*tanh(z/w0)*dzM
         (rigid +z translation of the upper half, -z of the lower),
         damped leapfrog (M5.21.6 form; dt = 0.025 from the verified
         gates) with the absorbing sponge pushed out to r0 = 0.8*half
         so the rings at r ~ 16 sit in the undamped interior, then
         FIRE from the evolved state to land the product basin.

Reads per endpoint (and per evolve snap, light form): ring ledger =
connected components of the low-eigengap set (thr ladder), classified
column (rho < 3) vs ring (rho >= 3); per-ring fragment charge =
oriented Mermin-Ho cube flux centered on the ring z; total Q_far;
E parts. Instruments reused from the census driver (m5_22_b_census).

Convention: signed charges below are MATHEMATICAL topological charge;
the electric reading is the NEGATED value (the author's 2026-07-30
hedgehog convention, census note § 7).

Modes:
  split tag=P-1_plane_sc6_n32_pinned_d0.3 d=3 w0=6 maxit=12000
  kick  tag=... v=0.5 steps=3000 w0=6 maxit=8000 [r0frac=0.8]
Out: ../data/m5_22_1_end_<newtag>.npz + m5_22_1_row_<newtag>.json
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
from scipy import ndimage
from scipy.ndimage import map_coordinates

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")


def _load(name, fname):
    sp = importlib.util.spec_from_file_location(
        name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(mod)
    return mod


CEN = _load("cen", "m5_22_b_census.py")
INS = CEN.INS
PAIR = CEN.PAIR
W2_T2 = CEN.W2_T2
DT = 0.025                       # GL1-verified leapfrog step (M5.21.6)


# ================= geometry helpers =================
def rgrid(cfg):
    X, Y, Z = INS.coords(cfg["n"], cfg["h"])
    return np.sqrt(X * X + Y * Y + Z * Z)


def separate(M, cfg, d, w0):
    """resample along z with z_src = z - d*tanh(z/w0): structure at
    +z moves up by ~d*tanh, at -z down; equator fixed. Monotone for
    d < w0. Shell values re-imposed by the caller."""
    n, h = cfg["n"], cfg["h"]
    c = (n - 1) / 2.0
    z_phys = (np.arange(n) - c) * h
    src_idx = (z_phys - d * np.tanh(z_phys / w0)) / h + c
    ii = np.broadcast_to(np.arange(n)[:, None, None], (n, n, n))
    jj = np.broadcast_to(np.arange(n)[None, :, None], (n, n, n))
    kk = np.broadcast_to(src_idx[None, None, :], (n, n, n))
    coords = np.stack([ii, jj, kk], axis=0).astype(float)
    out = np.empty_like(M)
    for a in range(3):
        for b in range(a, 3):
            v = map_coordinates(M[..., a, b], coords, order=1,
                                mode="nearest")
            out[..., a, b] = v
            out[..., b, a] = v
    return out


def sponge(cfg, g_max=0.5, r0frac=0.8):
    r = rgrid(cfg)
    half = cfg["L"] / 2.0
    r0, r1 = r0frac * half, 0.98 * half
    s = np.clip((r - r0) / max(r1 - r0, 1e-9), 0.0, 1.0)
    return (g_max * s * s)[..., None, None]


def leap_step(M, V, cfg, free, gam, dt):
    h3 = cfg["h"] ** 3
    F = -INS.grad(M, cfg) / h3
    V = (V + 0.5 * dt * F) * free
    V = V / (1.0 + 0.5 * dt * gam)
    M = M + dt * V
    F = -INS.grad(M, cfg) / h3
    V = (V + 0.5 * dt * F) * free
    V = V / (1.0 + 0.5 * dt * gam)
    return M, V


# ================= the ring ledger =================
def ring_read(M, cfg, thr):
    """low-eigengap components inside r < 0.9*half, pin shell
    excluded; column parts (rho_mean < 3) vs ring parts (>= 3)."""
    n, h = cfg["n"], cfg["h"]
    lam = np.linalg.eigvalsh(M)
    gap = np.minimum(lam[..., 1] - lam[..., 0],
                     lam[..., 2] - lam[..., 1])
    X, Y, Z = INS.coords(n, h)
    rho = np.sqrt(X * X + Y * Y)
    interior = (~INS.pin_shell(n, h)) & (rgrid(cfg) < 0.9 * cfg["L"] / 2)
    mask = (gap < thr) & interior
    lab, ncomp = ndimage.label(mask, structure=np.ones((3, 3, 3), int))
    cols, rings = [], []
    for k in range(1, ncomp + 1):
        sel = lab == k
        size = int(sel.sum())
        if size < 2:
            continue
        entry = {"voxels": size,
                 "rho_mean": float(rho[sel].mean()),
                 "z_centroid": float(Z[sel].mean()),
                 "z_min": float(Z[sel].min()),
                 "z_max": float(Z[sel].max()),
                 "gap_min": float(gap[sel].min())}
        (cols if entry["rho_mean"] < 3.0 else rings).append(entry)
    rings.sort(key=lambda e: e["z_centroid"])
    cols.sort(key=lambda e: e["z_centroid"])
    return {"thr": thr, "n_rings": len(rings), "rings": rings,
            "n_column_parts": len(cols), "columns": cols[:6]}


def slab_flux(B, cfg, z0, z1, half_lat):
    """Q = (1/4pi) sum_faces B.dS over the box |x|,|y| <= half_lat,
    z0 <= z <= z1 (physical units): the fragment box, laterally wide
    enough to ENCLOSE a ring at rho ~ 10 (a centered cube is not)."""
    n, h = cfg["n"], cfg["h"]
    c = (n - 1) / 2.0
    k = int(round(half_lat / h))
    i0, i1 = int(round(c)) - k, int(round(c)) + k
    k0 = max(1, int(round(z0 / h + c)))
    k1 = min(n - 2, int(round(z1 / h + c)))
    if i0 < 1 or i1 > n - 2 or k1 <= k0:
        return float("nan")
    s = 0.0
    s += B[i1, i0:i1 + 1, k0:k1 + 1, 0].sum() - \
        B[i0, i0:i1 + 1, k0:k1 + 1, 0].sum()
    s += B[i0:i1 + 1, i1, k0:k1 + 1, 1].sum() - \
        B[i0:i1 + 1, i0, k0:k1 + 1, 1].sum()
    s += B[i0:i1 + 1, i0:i1 + 1, k1, 2].sum() - \
        B[i0:i1 + 1, i0:i1 + 1, k0, 2].sum()
    return float(s * h * h / (4.0 * np.pi))


def fragment_charges(M, cfg, ring_zs, dz=6.0):
    nhat, ncf = PAIR.orient_v1(M)
    B = PAIR.mermin_B(nhat, cfg["h"])
    far0 = 0.5 * cfg["L"] - 4.0 * cfg["h"]
    out = {"n_conflicts": int(ncf), "fragments": []}
    for zc in ring_zs:
        out["fragments"].append({
            "z_center": zc,
            "q_slab": slab_flux(B, cfg, zc - dz, zc + dz, far0)})
    h = cfg["h"]
    out["q_upper_half"] = slab_flux(B, cfg, 0.5 * h, far0, far0)
    out["q_lower_half"] = slab_flux(B, cfg, -far0, -0.5 * h, far0)
    out["q_far"] = PAIR.cube_flux(B, cfg, 0.0, far0)
    return out


def endpoint_reads(M, cfg, thresholds=(0.06, 0.09, 0.15)):
    reads = {f"thr{t:g}": ring_read(M, cfg, t) for t in thresholds}
    zs = [r["z_centroid"] for r in reads["thr0.09"]["rings"]]
    if not zs:
        zs = [r["z_centroid"] for r in reads["thr0.15"]["rings"]]
    e_u, e_d, e_v = INS.e_parts(M, cfg)
    return {"E_end": float(e_u + e_d + e_v), "E_u": float(e_u),
            "E_d": float(e_d), "E_v": float(e_v),
            "rings": reads,
            "charges": fragment_charges(M, cfg, zs)}


# ================= branches =================
def _load_parent(tag):
    z = np.load(os.path.join(DATA, f"m5_22_end_{tag}.npz"))
    cfg = INS.base_cfg(term="T2", stencil="sym", eps=0.0, w2=W2_T2,
                       n=int(z["n"]), delta=float(z["delta"]),
                       bc="pinned")
    return z["M"].astype(np.float64), cfg


def _finish(row, M, M_start, cfg, newtag, t0):
    row.update(endpoint_reads(M, cfg))
    row["wall_s"] = time.time() - t0
    np.savez_compressed(os.path.join(DATA, f"m5_22_1_end_{newtag}.npz"),
                        M=M.astype(np.float32),
                        M0=M_start.astype(np.float32),
                        delta=cfg["delta"], h=cfg["h"], n=cfg["n"])
    with open(os.path.join(DATA, f"m5_22_1_row_{newtag}.json"),
              "w") as f:
        json.dump(row, f, indent=1)
    r9 = row["rings"]["thr0.09"]
    print(json.dumps({
        "tag": newtag, "E_end": row["E_end"], "stop": row.get("stop"),
        "n_rings_thr0.09": r9["n_rings"],
        "ring_zs": [r["z_centroid"] for r in r9["rings"]],
        "q_far": row["charges"]["q_far"],
        "fragments": row["charges"]["fragments"]}))
    return row


def split(kw):
    tag = kw.pop("tag")
    d = float(kw.pop("d", 3.0))
    w0 = float(kw.pop("w0", max(6.0, 1.5 * d)))
    maxit = int(kw.pop("maxit", 12000))
    newtag = f"{tag}_split_d{d:g}"
    M, cfg = _load_parent(tag)
    cfg["maxit"] = maxit
    free = ~INS.pin_shell(cfg["n"], cfg["h"])
    Ms = separate(M, cfg, d, w0)
    Ms = np.where(free[..., None, None], Ms, M)
    Ms = 0.5 * (Ms + np.swapaxes(Ms, -1, -2))
    e0 = float(sum(INS.e_parts(Ms, cfg)))
    t0 = time.time()
    Mr, _, info = INS.fire(Ms, cfg, free, max_iter=maxit,
                           log_every=1000, tag=newtag)
    row = {"branch": "split", "parent": tag, "d": d, "w0": w0,
           "maxit": maxit, "tag": newtag, "E_separated": e0,
           "stop": info["stop"], "trace": info["trace"][-3:]}
    return _finish(row, Mr, Ms, cfg, newtag, t0)


def kick(kw):
    tag = kw.pop("tag")
    v = float(kw.pop("v", 0.5))
    steps = int(kw.pop("steps", 3000))
    w0 = float(kw.pop("w0", 6.0))
    maxit = int(kw.pop("maxit", 8000))
    r0frac = float(kw.pop("r0frac", 0.8))
    newtag = f"{tag}_kick_v{v:g}"
    M, cfg = _load_parent(tag)
    n, h = cfg["n"], cfg["h"]
    free_b = ~INS.pin_shell(n, h)
    free = free_b[..., None, None].astype(float)
    X, Y, Z = INS.coords(n, h)
    dzM = np.empty_like(M)
    dzM[:, :, 1:-1] = (M[:, :, 2:] - M[:, :, :-2]) / (2 * h)
    dzM[:, :, 0] = dzM[:, :, 1]
    dzM[:, :, -1] = dzM[:, :, -2]
    V = -v * np.tanh(Z / w0)[..., None, None] * dzM * free
    ke0 = 0.5 * h ** 3 * float(np.sum(V * V))
    gam = sponge(cfg, r0frac=r0frac)
    e_start = float(sum(INS.e_parts(M, cfg)))
    hist = []
    absorbed = 0.0
    t0 = time.time()
    Mev = M.copy()
    for it in range(1, steps + 1):
        absorbed += float(np.sum(gam * V * V)) * h ** 3 * DT
        Mev, V = leap_step(Mev, V, cfg, free, gam, DT)
        if it % 500 == 0 or it == steps:
            e = float(sum(INS.e_parts(Mev, cfg)))
            ke = 0.5 * h ** 3 * float(np.sum(V * V))
            rr = ring_read(Mev, cfg, 0.09)
            row_s = {"it": it, "E": e, "KE": ke, "absorbed": absorbed,
                     "n_rings": rr["n_rings"],
                     "ring_zs": [r["z_centroid"] for r in rr["rings"]]}
            hist.append(row_s)
            print(f"  {newtag} it {it:5d} E {e:9.4f} KE {ke:8.4f} "
                  f"abs {absorbed:8.4f} rings {rr['n_rings']} "
                  f"zs {[round(z, 1) for z in row_s['ring_zs']]} "
                  f"[{time.time() - t0:.0f}s]", flush=True)
    cfg["maxit"] = maxit
    Mr, _, info = INS.fire(Mev, cfg, free_b, max_iter=maxit,
                           log_every=1000, tag=newtag)
    row = {"branch": "kick", "parent": tag, "v": v, "steps": steps,
           "w0": w0, "dt": DT, "r0frac": r0frac, "maxit": maxit,
           "tag": newtag, "E_start": e_start, "KE_kick": ke0,
           "hist": hist, "stop": info["stop"],
           "trace": info["trace"][-3:]}
    return _finish(row, Mr, Mev, cfg, newtag, t0)


def reads(kw):
    tag = kw.pop("tag")
    src = kw.pop("src", "m5_22")
    z = np.load(os.path.join(DATA, f"{src}_end_{tag}.npz"))
    cfg = INS.base_cfg(term="T2", stencil="sym", eps=0.0, w2=W2_T2,
                       n=int(z["n"]), delta=float(z["delta"]),
                       bc="pinned")
    out = endpoint_reads(z["M"].astype(np.float64), cfg)
    r9 = out["rings"]["thr0.09"]
    print(json.dumps({
        "tag": tag, "E_end": out["E_end"],
        "n_rings_thr0.09": r9["n_rings"],
        "ring_zs": [r["z_centroid"] for r in r9["rings"]],
        "q_far": out["charges"]["q_far"],
        "fragments": out["charges"]["fragments"]}, indent=1))
    return out


def parse_kv(argv):
    kw = {}
    for a in argv:
        k, val = a.split("=", 1)
        kw[k] = val
    return kw


if __name__ == "__main__":
    ARGV = sys.argv[1:]
    mode = ARGV[0]
    if mode == "split":
        split(parse_kv(ARGV[1:]))
    elif mode == "kick":
        kick(parse_kv(ARGV[1:]))
    elif mode == "reads":
        reads(parse_kv(ARGV[1:]))
