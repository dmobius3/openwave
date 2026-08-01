"""M5.22.1: the deuteron construction forks (two-knot bound state).

Rung 3 of the author ladder: a |Q| = 1, two-baryon bound state from
the census constituents (proton-analog P-0.5: column + 1 ring; the
neutral P-1: column + 2 rings, the rings carrying OPPOSITE unit
topological slab charges per the m5_22_1_a_kick baseline read).

Construction forks, ALL self-run (series exhaustion rule; the author's
analytic deuteron seed folds in on arrival, construction does not
wait):

  seed2  analytic two-center composite: the author's P-family
         cross-section angles ADDED with axial offsets,
         ang(x, y) = ang_P(x - x1, y; s1) + ang_P(x - x2, y; s2),
         lifted by the census rotation (m5_22_a_seeds lift), per-term
         ring-core blends, then FIRE. Fork axes: (s1, s2) pairing and
         the separation a (author units; lattice offset = a * scale).
  graft  endpoint composition: the relaxed P-0.5 field shifted down
         meets the relaxed P-1 field shifted up across a tanh(z/w0)
         blend, shell re-pinned to the blend values, then FIRE.

Reads per endpoint: E vs the constituent sums (binding sign), ring
ledger (ring count = baryon number, the author's 2026-07-30
interpretation, tested against additivity), fragment slab charges,
Q_far, and the MOMENTS of the topological charge density
rho = div B / 4pi (B = oriented Mermin-Ho flux): dipole p_z and the
quadrupole Q2_zz = sum rho (3 Z^2 - r^2) h^3.

Convention: signed values are MATHEMATICAL topological charge; the
electric reading NEGATES them (the author's 2026-07-30 hedgehog
convention, census note § 7), so electric p_z / Q2_zz = -p_z / -Q2_zz.

Modes:
  seed2   s1=-0.5 s2=-1 a=2 [scale=6 n=32 delta=0.3 maxit=12000 rc=2]
  graft   ptag=P-0.5_plane_sc6_n32_pinned_d0.3
          ntag=P-1_plane_sc6_n32_pinned_d0.3 zoff=6 w0=4 [maxit=12000]
  moments tag=<endpoint tag> [src=m5_22|m5_22_1]
Out: ../data/m5_22_1_end_<tag>.npz + m5_22_1_row_<tag>.json
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

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")


def _load(name, fname):
    sp = importlib.util.spec_from_file_location(
        name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(mod)
    return mod


KICK = _load("kick", "m5_22_1_a_kick.py")
CEN = KICK.CEN
INS = KICK.INS
PAIR = KICK.PAIR
SEEDS = CEN.SEEDS
W2_T2 = KICK.W2_T2


# ================= moments of the charge density =================
def moments(M, cfg):
    """dipole + quadrupole of rho = div B / 4pi over the interior
    (pin shell excluded). Signed TOPOLOGICAL values; electric = -."""
    n, h = cfg["n"], cfg["h"]
    nhat, _ = PAIR.orient_v1(M)
    B = PAIR.mermin_B(nhat, h)
    div = np.zeros(B.shape[:3])
    for ax in range(3):
        d = np.zeros(B.shape[:3])
        sl = [slice(None)] * 3
        sp, sm = list(sl), list(sl)
        sl[ax] = slice(1, -1)
        sp[ax] = slice(2, None)
        sm[ax] = slice(None, -2)
        d[tuple(sl)] = (B[tuple(sp) + (ax,)]
                        - B[tuple(sm) + (ax,)]) / (2 * h)
        div += d
    rho = div / (4.0 * np.pi)
    X, Y, Z = INS.coords(n, h)
    interior = ~INS.pin_shell(n, h)
    r2 = X * X + Y * Y + Z * Z
    w = rho * interior
    h3 = h ** 3
    return {"q_vol": float(w.sum() * h3),
            "p_z": float((w * Z).sum() * h3),
            "Q2_zz": float((w * (3 * Z * Z - r2)).sum() * h3),
            "Q2_xx": float((w * (3 * X * X - r2)).sum() * h3)}


# ================= seed2: analytic two-center composite ============
def seed2_field(cfg, pairs, conv="plane", scale=6.0, r_c=2.0):
    """pairs = [(s1, x1), (s2, x2)] in author units; P family only."""
    n, h, delta = cfg["n"], cfg["h"], cfg["delta"]
    X, Y, Z = INS.coords(n, h)
    rho = np.sqrt(X * X + Y * Y)
    rhos = np.where(rho < 1e-12, 1e-12, rho)
    rhohat = np.stack([X / rhos, Y / rhos, np.zeros_like(Z)], axis=-1)
    near = rho < 1e-9
    if np.any(near):
        rhohat[near] = np.array([1.0, 0.0, 0.0])
    zhat = np.zeros_like(rhohat)
    zhat[..., 2] = 1.0
    ang = np.zeros_like(rho)
    w = np.ones_like(rho)
    for s, x0 in pairs:
        ang = ang + SEEDS.ang_family("P", s, Z / scale - x0, rho / scale)
        d_ring = np.sqrt((rho - SEEDS.R_B * scale) ** 2
                         + (Z - x0 * scale) ** 2)
        w = w * (1.0 - np.exp(-((d_ring / r_c) ** 2)))
    ca, sa = np.cos(ang), np.sin(ang)
    n1 = ca[..., None] * zhat + sa[..., None] * rhohat
    n2 = -sa[..., None] * zhat + ca[..., None] * rhohat
    S = (n1[..., :, None] * n1[..., None, :]
         + delta * n2[..., :, None] * n2[..., None, :])
    w = w * (1.0 - sa * sa * np.exp(-((rho / r_c) ** 2)))
    a_iso = (1.0 + delta) / 3.0
    return w[..., None, None] * S \
        + (1.0 - w[..., None, None]) * (a_iso * np.eye(3))


def _relax_and_finish(M0, cfg, tag, meta):
    free = ~INS.pin_shell(cfg["n"], cfg["h"])
    t0 = time.time()
    M, _, info = INS.fire(M0, cfg, free, max_iter=cfg["maxit"],
                          log_every=1000, tag=tag)
    row = dict(meta)
    row.update({"tag": tag, "E_seed": float(sum(INS.e_parts(M0, cfg))),
                "stop": info["stop"], "trace": info["trace"][-3:],
                "n": cfg["n"], "delta": cfg["delta"],
                "maxit": cfg["maxit"]})
    row.update(KICK.endpoint_reads(M, cfg))
    row["moments"] = moments(M, cfg)
    row["wall_s"] = time.time() - t0
    np.savez_compressed(os.path.join(DATA, f"m5_22_1_end_{tag}.npz"),
                        M=M.astype(np.float32),
                        M0=M0.astype(np.float32),
                        delta=cfg["delta"], h=cfg["h"], n=cfg["n"])
    with open(os.path.join(DATA, f"m5_22_1_row_{tag}.json"), "w") as f:
        json.dump(row, f, indent=1)
    r9 = row["rings"]["thr0.09"]
    r15 = row["rings"]["thr0.15"]
    print(json.dumps({
        "tag": tag, "E_end": row["E_end"], "stop": row["stop"],
        "rings_thr0.09": r9["n_rings"], "rings_thr0.15": r15["n_rings"],
        "ring_zs_thr0.15": [r["z_centroid"] for r in r15["rings"]],
        "q_far": row["charges"]["q_far"],
        "q_upper": row["charges"]["q_upper_half"],
        "q_lower": row["charges"]["q_lower_half"],
        "p_z": row["moments"]["p_z"],
        "Q2_zz": row["moments"]["Q2_zz"]}))
    return row


def seed2(kw):
    s1 = float(kw.pop("s1", -0.5))
    s2 = float(kw.pop("s2", -1.0))
    a = float(kw.pop("a", 2.0))
    scale = float(kw.pop("scale", 6.0))
    r_c = float(kw.pop("rc", 2.0))
    cfg = INS.base_cfg(term="T2", stencil="sym", eps=0.0, w2=W2_T2,
                       n=int(kw.pop("n", 32)),
                       delta=float(kw.pop("delta", 0.3)),
                       bc="pinned", maxit=int(kw.pop("maxit", 12000)))
    tag = (f"d2_s{s1:+g}_s{s2:+g}_a{a:g}_sc{scale:g}_n{cfg['n']}"
           f"_d{cfg['delta']:g}")
    pairs = [(s1, -a / 2.0), (s2, a / 2.0)]
    M0 = seed2_field(cfg, pairs, scale=scale, r_c=r_c)
    return _relax_and_finish(M0, cfg, tag,
                             {"branch": "seed2", "s1": s1, "s2": s2,
                              "a": a, "scale": scale, "r_c": r_c})


def seedn(kw):
    """arbitrary multi-center composite: pairs=s:x,s:x,... in author
    units (P family). The one-consistent-far-field fix to the measured
    graft/additivity obstructions: pick the s_i so the TOTAL 2D
    winding is odd (no 3D escape) while interior rings carry the
    target content."""
    spec = kw.pop("pairs")
    pairs = []
    for item in spec.split(","):
        s, x0 = item.split(":")
        pairs.append((float(s), float(x0)))
    scale = float(kw.pop("scale", 6.0))
    r_c = float(kw.pop("rc", 2.0))
    cfg = INS.base_cfg(term="T2", stencil="sym", eps=0.0, w2=W2_T2,
                       n=int(kw.pop("n", 32)),
                       delta=float(kw.pop("delta", 0.3)),
                       bc="pinned", maxit=int(kw.pop("maxit", 12000)))
    tag = ("dn_" + "_".join(f"{s:+g}at{x0:+g}" for s, x0 in pairs)
           + f"_n{cfg['n']}_d{cfg['delta']:g}")
    M0 = seed2_field(cfg, pairs, scale=scale, r_c=r_c)
    return _relax_and_finish(M0, cfg, tag,
                             {"branch": "seedn", "pairs": pairs,
                              "scale": scale, "r_c": r_c})


def graft(kw):
    ptag = kw.pop("ptag")
    ntag = kw.pop("ntag")
    zoff = float(kw.pop("zoff", 6.0))
    w0 = float(kw.pop("w0", 4.0))
    maxit = int(kw.pop("maxit", 12000))
    zp = np.load(os.path.join(DATA, f"m5_22_end_{ptag}.npz"))
    zn = np.load(os.path.join(DATA, f"m5_22_end_{ntag}.npz"))
    cfg = INS.base_cfg(term="T2", stencil="sym", eps=0.0, w2=W2_T2,
                       n=int(zp["n"]), delta=float(zp["delta"]),
                       bc="pinned", maxit=maxit)
    n, h = cfg["n"], cfg["h"]
    Mp = KICK.separate(zp["M"].astype(np.float64), cfg, -zoff,
                       3.0 * zoff)
    Mn = KICK.separate(zn["M"].astype(np.float64), cfg, zoff,
                       3.0 * zoff)
    _, _, Z = INS.coords(n, h)
    wz = 0.5 * (1.0 + np.tanh(Z / w0))[..., None, None]
    M0 = (1.0 - wz) * Mp + wz * Mn
    M0 = 0.5 * (M0 + np.swapaxes(M0, -1, -2))
    tag = f"graft_p_below_n_above_z{zoff:g}_n{n}_d{cfg['delta']:g}"
    return _relax_and_finish(M0, cfg, tag,
                             {"branch": "graft", "ptag": ptag,
                              "ntag": ntag, "zoff": zoff, "w0": w0})


def extend(kw):
    """continue a saved m5_22_1 endpoint: extend tag=... maxit=..."""
    tag = kw.pop("tag")
    maxit = int(kw.pop("maxit", 12000))
    z = np.load(os.path.join(DATA, f"m5_22_1_end_{tag}.npz"))
    cfg = INS.base_cfg(term="T2", stencil="sym", eps=0.0, w2=W2_T2,
                       n=int(z["n"]), delta=float(z["delta"]),
                       bc="pinned", maxit=maxit)
    M0 = z["M"].astype(np.float64)
    return _relax_and_finish(M0, cfg, tag + "_ext",
                             {"branch": "extend", "parent": tag})


def moments_mode(kw):
    tag = kw.pop("tag")
    src = kw.pop("src", "m5_22")
    z = np.load(os.path.join(DATA, f"{src}_end_{tag}.npz"))
    cfg = INS.base_cfg(term="T2", stencil="sym", eps=0.0, w2=W2_T2,
                       n=int(z["n"]), delta=float(z["delta"]),
                       bc="pinned")
    out = moments(z["M"].astype(np.float64), cfg)
    print(json.dumps({"tag": tag} | out, indent=1))
    return out


if __name__ == "__main__":
    ARGV = sys.argv[1:]
    mode = ARGV[0]
    kw = KICK.parse_kv(ARGV[1:])
    if mode == "seed2":
        seed2(kw)
    elif mode == "seedn":
        seedn(kw)
    elif mode == "extend":
        extend(kw)
    elif mode == "graft":
        graft(kw)
    elif mode == "moments":
        moments_mode(kw)
