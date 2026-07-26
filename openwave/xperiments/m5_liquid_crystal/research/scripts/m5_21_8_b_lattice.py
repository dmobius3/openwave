"""M5.21.8 P2/P3: the analytical ansatz confronted on the 4D lattice.

Builds the author's dressed family M(x) = Qb(m) Qh d Qh^T Qb(m)^T
(boost hedgehog x twisting hedgehog, eigenvalues pinned to the
vacuum spectrum) ON THE GRID of the certified M5.21.3 4D stack
(imported), and measures with that stack's audited instruments:

    E_stat(m)   e_parts on the dressed static state (V4 = 0 exactly
                on the family: eigenvalue-pinned); pure E_u curve
    kin(m)      kin_of with the ansatz's own unit-omega clock flow
                a0 = dM/dt|_{t=0} (analytic family, FD in t);
                E_kin = omega^2 * kin under the eta-reading
                (P0 bridge: = 2x his convention, same signs)
    box probe   kin at L = 48 (n = 32) vs L = 36 (n = 24, same h):
                his constant-omega kinetic is IR-extensive (~L)
    relax       FIRE survival probe from the dressed minimum
    gladder     m*_lattice(g) + kin(m*, g) vs his
                m* = (1/2) ln((1+g)/(g-1)) and the 1/g shrink claim

Modes: mcurve s=-1 | relax m=<val> s=-1 maxit=3000 | gladder | box
Out: data/m5_21_8_lat_*.json
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")

_spec = importlib.util.spec_from_file_location(
    "ins4", os.path.join(HERE, "m5_21_3_a_4d.py"))
INS4 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(INS4)

G1 = np.zeros((4, 4)); G1[2, 3] = -1.0; G1[3, 2] = 1.0
G2 = np.zeros((4, 4)); G2[1, 3] = 1.0; G2[3, 1] = -1.0
G3 = np.zeros((4, 4)); G3[1, 2] = -1.0; G3[2, 1] = 1.0


def rot_field(G, A):
    """per-cell Rodrigues: I + sin(A) G + (1-cos(A)) G^2."""
    G2m = G @ G
    R = (np.eye(4)[None, None, None]
         + np.sin(A)[..., None, None] * G[None, None, None]
         + (1 - np.cos(A))[..., None, None] * G2m[None, None, None])
    return R


def dressed(cfg, m, t=0.0, om=1.0):
    """the author's ansatz on the grid (d = vac4: his d at s=-1)."""
    n, h = cfg["n"], cfg["h"]
    X, Y, Z = INS4.coords(n, h)
    R = np.sqrt(X * X + Y * Y + Z * Z)
    rho = np.sqrt(X * X + Y * Y)
    phi = np.arctan2(Y, X)
    th = -np.arctan2(Z, rho)
    Qh = np.einsum("...ab,...bc,...cd->...ad",
                   rot_field(G3, phi), rot_field(G2, th),
                   rot_field(G1, om * t * np.ones_like(phi)))
    nx, ny, nz = X / R, Y / R, Z / R
    K = np.zeros(X.shape + (4, 4))
    K[..., 0, 1], K[..., 0, 2], K[..., 0, 3] = nx, ny, nz
    K[..., 1, 0], K[..., 2, 0], K[..., 3, 0] = nx, ny, nz
    K2 = np.zeros_like(K)
    K2[..., 0, 0] = 1.0
    for i, a in enumerate((nx, ny, nz)):
        for j, b in enumerate((nx, ny, nz)):
            K2[..., 1 + i, 1 + j] = a * b
    Qb = (np.eye(4)[None, None, None] + np.sinh(m) * K
          + (np.cosh(m) - 1.0) * K2)
    Q = np.einsum("...ab,...bc->...ac", Qb, Qh)
    d4 = INS4.vac4(cfg)
    M = np.einsum("...ab,bc,...dc->...ad", Q, d4, Q)
    return INS4.sym4(M)


def a0_unit(cfg, m, dt=1e-4):
    """unit-omega clock flow of the family at t = 0."""
    return (dressed(cfg, m, t=dt, om=1.0)
            - dressed(cfg, m, t=-dt, om=1.0)) / (2 * dt)


def mcurve(s=-1.0, g=8.0, n=32, L=48.0, tag=None,
           m0=-0.35, m1=0.35, nm=29):
    cfg = INS4.base_cfg(s=s, g=g, n=n, L=L)
    ms = np.linspace(m0, m1, nm)
    rows = []
    for m in ms:
        M = dressed(cfg, m)
        e_u, e_v = INS4.e_parts(M, cfg)
        a0 = a0_unit(cfg, m)
        k = INS4.kin_of(M, a0, cfg)
        rows.append({"m": float(m), "E_u": float(e_u),
                     "E_v": float(e_v), "kin": float(k)})
        print(f"  s={s:+.0f} g={g:g} m {m:+.3f} E_u {e_u:10.4f} "
              f"E_v {e_v:.2e} kin {k:+.4f}", flush=True)
    Es = [r["E_u"] + r["E_v"] for r in rows]
    i = int(np.argmin(Es))
    mfit = ms[i]
    if 0 < i < len(ms) - 1:
        a, b, c = Es[i - 1], Es[i], Es[i + 1]
        mfit = ms[i] - 0.5 * (ms[1] - ms[0]) * (c - a) / (c - 2 * b + a)
    out = {"s": s, "g": g, "n": n, "L": L, "rows": rows,
           "m_star_lattice": float(mfit),
           "m_star_his": float(0.5 * np.log((1 + g) / (g - 1))),
           "E_min": float(min(Es))}
    tg = tag or f"mcurve_s{int(s)}_g{g:g}_n{n}"
    with open(os.path.join(DATA, f"m5_21_8_lat_{tg}.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps({k: out[k] for k in
                      ("m_star_lattice", "m_star_his", "E_min")}))
    return out


def box():
    rows = []
    for (n, L) in ((32, 48.0), (24, 36.0)):
        cfg = INS4.base_cfg(s=-1.0, g=8.0, n=n, L=L)
        m = 0.5 * np.log(9.0 / 7.0)
        M = dressed(cfg, m)
        a0 = a0_unit(cfg, m)
        rows.append({"n": n, "L": L,
                     "kin": float(INS4.kin_of(M, a0, cfg))})
        print(rows[-1], flush=True)
    rows.append({"kin_ratio": rows[0]["kin"] / rows[1]["kin"],
                 "L_ratio": 48.0 / 36.0})
    print(rows[-1])
    with open(os.path.join(DATA, "m5_21_8_lat_box.json"), "w") as f:
        json.dump(rows, f, indent=1)
    return rows


def relax(m, s=-1.0, maxit=3000, dt0=0.02, dt_max=0.2):
    cfg = INS4.base_cfg(s=s, g=8.0)
    M0 = dressed(cfg, m)
    n = cfg["n"]
    free = ~INS4.pin_shell(n, cfg["h"])
    tag = f"relax_m{m:+.3f}_s{int(s)}_dt{dt0:g}"
    M, info = INS4.fire(M0, cfg, free, max_iter=maxit,
                        log_every=250, tag=tag, dt0=dt0,
                        dt_max=dt_max)
    e0 = INS4.e_parts(M0, cfg)
    e1 = INS4.e_parts(M, cfg)
    dv = float(np.sqrt(np.sum((M - M0) ** 2))
               / max(np.sqrt(np.sum((M0 - INS4.vac4(cfg)) ** 2)),
                     1e-300))
    a0 = a0_unit(cfg, m)
    out = {"m": m, "s": s, "E0": float(sum(e0)),
           "E_end": float(sum(e1)),
           "E_u_end": float(e1[0]), "E_v_end": float(e1[1]),
           "rel_move": dv, "kin_end": float(INS4.kin_of(M, a0, cfg)),
           "stop": info["stop"]}
    with open(os.path.join(DATA, f"m5_21_8_lat_{tag}.json"),
              "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out))
    return out


def gladder():
    outs = []
    for g in (8.0, 16.0, 32.0, 64.0):
        o = mcurve(s=-1.0, g=g, tag=f"glad_g{g:g}")
        i = int(np.argmin([abs(r["m"] - o["m_star_lattice"])
                           for r in o["rows"]]))
        outs.append({"g": g, "m_star_lattice": o["m_star_lattice"],
                     "m_star_his": o["m_star_his"],
                     "kin_at_mstar": o["rows"][i]["kin"],
                     "E_min": o["E_min"]})
        print(outs[-1], flush=True)
    with open(os.path.join(DATA, "m5_21_8_lat_gladder.json"),
              "w") as f:
        json.dump(outs, f, indent=1)
    return outs


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "mcurve"
    kw = dict(a.split("=", 1) for a in sys.argv[2:])
    if mode == "mcurve":
        mcurve(s=float(kw.get("s", -1)), g=float(kw.get("g", 8)),
               n=int(kw.get("n", 32)), L=float(kw.get("L", 48)),
               tag=kw.get("tag"), m0=float(kw.get("m0", -0.35)),
               m1=float(kw.get("m1", 0.35)), nm=int(kw.get("nm", 29)))
    elif mode == "box":
        box()
    elif mode == "relax":
        relax(float(kw["m"]), s=float(kw.get("s", -1)),
              maxit=int(kw.get("maxit", 3000)),
              dt0=float(kw.get("dt0", 0.02)),
              dt_max=float(kw.get("dt_max", 0.2)))
    elif mode == "gladder":
        gladder()
