"""M5.21.4 films: the antipair capture run on the TRUE m5_film.py
templates (basic + thermal), the series film rule.

Same meridional embed as m5_21_6_c_films.py (y = 0 half-plane, 3x3
embedded 4x4 block-diag, display-constant time eigenvalue). Row
selection: t = 0 (the healed seed) + up to 6 rows chosen to bracket the
capture (dense around the fastest d(t) change, from the rows JSON).

Run:  python3 m5_21_4_c_films.py <tag>          (tag = cap, ...)
Out:  ../plots/m5_21_4_film_{basic,thermal}_<tag>.png
"""
import json
import os
import runpy
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
PLOTS = os.path.join(HERE, "..", "plots")
sys.path.insert(0, HERE)
g = runpy.run_path(os.path.join(HERE, "m5_21_2b_a_instrument.py"),
                   run_name="not_main")
from m5_film import film_strip  # noqa: E402

G_DISP = 8.0


def to_axisym4(M3):
    n = M3.shape[0]
    j = n // 2
    half = M3[n // 2:, j]
    M4 = np.zeros(half.shape[:2] + (4, 4))
    M4[..., 1:4, 1:4] = half
    M4[..., 0, 0] = G_DISP
    return M4


def dens3_slice(M3, cfg):
    h = cfg["h"]
    ud = np.zeros(M3.shape[:3])
    for br, wt in g["branches"](cfg["stencil"]):
        A = [g["d1"](M3, ax, h, br) for ax in range(3)]
        for i in range(3):
            for k in range(i + 1, 3):
                C = A[i] @ A[k] - A[k] @ A[i]
                ud += wt * 4.0 * np.einsum("...kl,...kl->...", C, C)
    vd = g["v_density"](M3, cfg)
    n = M3.shape[0]
    return (ud + vd)[n // 2:, n // 2]


def pick_rows(snap_its, rows, n_mid=5):
    """t=0 + the last snap + n_mid rows spread over the fastest d(t)
    window (largest per-interval |d dsep|)."""
    if len(snap_its) <= n_mid + 2:
        return snap_its
    dse = {r["it"]: r["dsep"] for r in rows}
    avail = [it for it in snap_its if it in dse]
    if len(avail) < 3:
        step = max(len(snap_its) // (n_mid + 1), 1)
        return snap_its[::step][:n_mid + 1] + [snap_its[-1]]
    scored = sorted(avail[1:-1],
                    key=lambda it: -abs(dse.get(it, 0.0)
                                        - dse.get(avail[avail.index(it) - 1],
                                                  0.0)))
    mids = sorted(scored[:n_mid])
    return [snap_its[0]] + mids + [snap_its[-1]]


def descent_strip(kind, d, n=32):
    """rows = the seed + the saved descent depths (120/400/1500 it):
    the sector's descent-fate film on the true templates."""
    import importlib.util
    sp = importlib.util.spec_from_file_location(
        "pair", os.path.join(HERE, "m5_21_4_a_pair.py"))
    P = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(P)
    delta = 0.3
    cfg = g["base_cfg"](term="T2", n=n, L=1.5 * n, stencil="sym",
                        delta=delta)
    states3 = [{"it": 0, "t": 0.0,
                "M3": P.seed_pair(cfg, kind, float(d))}]
    for it, name in [(120, f"m5_21_4_lad_{kind}_d{d:g}_n{n}_it120.npz"),
                     (400, f"m5_21_4_lad_{kind}_d{d:g}_n{n}_it400.npz"),
                     (1500, f"m5_21_4_lad_{kind}_d{d:g}_n{n}.npz")]:
        p = os.path.join(DATA, name)
        if os.path.exists(p):
            states3.append({"it": it, "t": float(it),
                            "M3": np.load(p)["M"].astype(np.float64)})
    dens = [dens3_slice(st["M3"], cfg) for st in states3]
    states = [{"it": st["it"], "t": st["t"],
               "M": to_axisym4(st["M3"])} for st in states3]
    it = iter(dens)

    def density_fn(M4):
        return next(it)
    tag = f"descent_{kind}_d{d:g}"
    film_strip(states, os.path.join(
        PLOTS, f"m5_21_4_film_basic_{tag}.png"),
        template="basic", delta=delta, h=cfg["h"], g=G_DISP,
        wscale=g["W1"], density_fn=density_fn,
        suptitle=f"M5.21.4 {kind} pair d={d:g} DESCENT ladder (n={n} "
                 f"pinned; rows = FIRE depth) on the TRUE basic template")
    film_strip(states, os.path.join(
        PLOTS, f"m5_21_4_film_thermal_{tag}.png"),
        template="thermal", delta=delta, h=cfg["h"], g=G_DISP,
        suptitle=f"M5.21.4 {kind} pair d={d:g} descent ladder on the "
                 f"TRUE thermal template")
    print("descent films written for", tag)


def main(tag="cap"):
    Z = np.load(os.path.join(DATA, f"m5_21_4_ev_{tag}.npz"))
    meta = json.load(open(os.path.join(DATA,
                                       f"m5_21_4_ev_{tag}_rows.json")))
    n = int(meta["cfg_n"])
    dt = float(meta["dt"])
    delta = 0.3
    cfg = g["base_cfg"](term="T2", n=n, L=1.5 * n, stencil="sym",
                        delta=delta)
    snap_its = sorted(int(k[len("M_it"):]) for k in Z.files
                      if k.startswith("M_it"))
    keep = pick_rows(snap_its, meta["rows"])
    states3 = [{"it": it, "t": float(it),
                "M3": Z[f"M_it{it}"].astype(np.float64)} for it in keep]
    dens = [dens3_slice(st["M3"], cfg) for st in states3]
    states = [{"it": st["it"], "t": st["t"],
               "M": to_axisym4(st["M3"])} for st in states3]
    it = iter(dens)

    def density_fn(M4):
        return next(it)
    film_strip(states, os.path.join(
        PLOTS, f"m5_21_4_film_basic_{tag}.png"),
        template="basic", delta=delta, h=cfg["h"], g=G_DISP,
        wscale=g["W1"], density_fn=density_fn,
        suptitle=f"M5.21.4 antipair capture {tag} (n={n} free, dt={dt}) "
                 f"on the TRUE basic template (y = 0 half-plane; time "
                 f"row = display constant g)")
    film_strip(states, os.path.join(
        PLOTS, f"m5_21_4_film_thermal_{tag}.png"),
        template="thermal", delta=delta, h=cfg["h"], g=G_DISP,
        suptitle=f"M5.21.4 antipair capture {tag} on the TRUE thermal "
                 f"template (y = 0 half-plane; time row = display "
                 f"constant g)")
    print("films written for", tag, "rows", keep)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "descent":
        descent_strip(sys.argv[2], float(sys.argv[3]),
                      int(sys.argv[4]) if len(sys.argv) > 4 else 32)
    else:
        main(sys.argv[1] if len(sys.argv) > 1 else "cap")
