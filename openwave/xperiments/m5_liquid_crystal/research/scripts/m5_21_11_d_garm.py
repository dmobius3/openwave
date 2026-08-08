"""M5.21.11 g-arm: rigid Qb(m) dressing READS on the delta = 0.30 endpoints.

Framework § 2 (FROZEN): evaluation only, no relaxation, on the N = 48
production endpoints at delta = 0.3: g in {8, 16, 32}, fine m-grids
(the coarse-grid anti-recipe respected), both (-g)^p signs at g = 8
(s = -1 the paper vacuum branch, s = +1 the flipped sign), per branch.
Measures: m*_lattice vs artanh(1/g), the dressing gain
DeltaE_dress(g) = E(m*) - E(0), and the fall exponent q from
log(-gain) vs log(artanh(1/g)).

F4 (frozen § 6): FAILS only if the gain falls SLOWER than artanh(1/g)^2
across the g-arm rungs; q >= 2 (any faster fall included) PASSES.

Energy code: the certified 4D instrument m5_21_3_a_4d.py imported
unchanged (embed34 + e_parts); the boost hedgehog Qb(m) is the
m5_21_8_b_lattice construction (K, K2 from the radial hedgehog), a
rigid dressing of the loaded endpoint, not new energy code.

Mode: all
Out: ../data/m5_21_11_garm.json + per-arm curves inside it
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")

_spec = importlib.util.spec_from_file_location(
    "ins4", os.path.join(HERE, "m5_21_3_a_4d.py"))
INS4 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(INS4)

BRANCHES = ("A", "C", "B")
DELTA = 0.3
N = 48


def qb_field(cfg, m):
    """rigid boost-hedgehog Qb(m) on the grid (m5_21_8 construction)."""
    n, h = cfg["n"], cfg["h"]
    X, Y, Z = INS4.coords(n, h)
    R = np.sqrt(X * X + Y * Y + Z * Z)
    Rs = np.where(R < 1e-12, 1e-12, R)
    nx, ny, nz = X / Rs, Y / Rs, Z / Rs
    K = np.zeros(X.shape + (4, 4))
    K[..., 0, 1], K[..., 0, 2], K[..., 0, 3] = nx, ny, nz
    K[..., 1, 0], K[..., 2, 0], K[..., 3, 0] = nx, ny, nz
    K2 = np.zeros_like(K)
    K2[..., 0, 0] = 1.0
    for i, a in enumerate((nx, ny, nz)):
        for j, b in enumerate((nx, ny, nz)):
            K2[..., 1 + i, 1 + j] = a * b
    return (np.eye(4)[None, None, None] + np.sinh(m) * K
            + (np.cosh(m) - 1.0) * K2)


def e_dressed(M4, cfg, m):
    Qb = qb_field(cfg, m)
    Md = np.einsum("...ab,...bc,...dc->...ad", Qb, M4, Qb)
    e_u, e_v = INS4.e_parts(INS4.sym4(Md), cfg)
    return float(e_u + e_v)


def arm(branch, g, s, nm=121, span=3.0):
    tag = f"t11lad_{branch}_n{N}_d{DELTA:g}"
    Z = np.load(os.path.join(DATA, f"m5_21_11_end_{tag}.npz"))
    M3 = Z["M"].astype(np.float64)
    cfg = INS4.base_cfg(s=s, g=g, n=N, L=48.0, delta=DELTA)
    M4 = INS4.embed34(M3, cfg)
    m_his = float(np.arctanh(1.0 / g))
    ms = np.linspace(-span * m_his, span * m_his, nm)
    Es = np.array([e_dressed(M4, cfg, m) for m in ms])
    E0 = e_dressed(M4, cfg, 0.0)
    i = int(np.argmin(Es))
    m_star, E_star = float(ms[i]), float(Es[i])
    if 0 < i < nm - 1:
        a, b, c = Es[i - 1], Es[i], Es[i + 1]
        dm = ms[1] - ms[0]
        m_star = float(ms[i] - 0.5 * dm * (c - a) / (c - 2 * b + a))
        E_star = float(b - 0.125 * (c - a) ** 2 / (c - 2 * b + a))
    out = {"branch": branch, "g": g, "s": s, "nm": nm,
           "m_his": m_his, "m_star": m_star, "E0": E0,
           "E_star": E_star, "gain": E_star - E0,
           "curve": [{"m": float(m), "E": float(E)}
                     for m, E in zip(ms[::6], Es[::6])]}
    print(json.dumps({k: out[k] for k in
                      ("branch", "g", "s", "m_his", "m_star",
                       "gain")}), flush=True)
    return out


def all_mode():
    arms = []
    for br in BRANCHES:
        for g in (8.0, 16.0, 32.0):
            arms.append(arm(br, g, -1.0))
        arms.append(arm(br, 8.0, +1.0))
    # the fall exponent per branch (s = -1 arms)
    fits = {}
    for br in BRANCHES:
        rows = [a for a in arms if a["branch"] == br and a["s"] == -1]
        x = np.log([np.arctanh(1.0 / a["g"]) for a in rows])
        y = np.log([max(-a["gain"], 1e-300) for a in rows])
        slope, icpt = np.polyfit(x, y, 1)
        pair_slopes = [float((y[i + 1] - y[i]) / (x[i + 1] - x[i]))
                       for i in range(len(x) - 1)]
        kappa = float(np.exp(icpt))
        fits[br] = {"q_lsq": float(slope), "pair_slopes": pair_slopes,
                    "kappa": kappa,
                    "gains": {f"g{a['g']:g}": a["gain"]
                              for a in rows},
                    "f4_pass": bool(slope >= 2.0 - 0.05)}
        print(json.dumps({"branch": br} | fits[br]), flush=True)
    out = {"arms": arms, "fits": fits,
           "f4_pass_all": bool(all(f["f4_pass"]
                                   for f in fits.values()))}
    with open(os.path.join(DATA, "m5_21_11_garm.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps({"f4_pass_all": out["f4_pass_all"]}))
    return out


if __name__ == "__main__":
    all_mode()
