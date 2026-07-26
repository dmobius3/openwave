"""M5.21.3 films: the 4D states on the TRUE m5_film.py templates.

Unlike the 2b adapter no embed is needed: the field IS 4x4. Two film
axes, labeled in the suptitle:
  - mode=period: the rotating state over ONE internal-clock period,
    7 shots at phase 0, pi/3, ..., 2pi (the conjugation orbit
    Lambda(theta) M Lambda^T with the run's generator); time axis =
    CLOCK PHASE, not descent iteration.
  - mode=descent: stored M_it snapshots (7-shot standard).
Density panel = this task's true density (u_eta + V4 per cell).

Run: python3 m5_21_3_c_films.py <npz> <tag> [mode] [gen]
Out: ../plots/m5_21_3_film_basic_<tag>.png + _thermal_
"""
import os
import runpy
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
from scipy.linalg import expm  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
PLOTS = os.path.join(HERE, "..", "plots")
sys.path.insert(0, HERE)
g4 = runpy.run_path(os.path.join(HERE, "m5_21_3_a_4d.py"),
                    run_name="not_main")
from m5_film import film_strip  # noqa: E402


def dens_slice(M4, cfg):
    h = cfg["h"]
    ETA = g4["ETA"]
    ud = np.zeros(M4.shape[:3])
    for br, wt in g4["branches"](cfg["stencil"]):
        A = [g4["d1"](M4, ax, h, br) for ax in range(3)]
        for i in range(3):
            for j in range(i + 1, 3):
                F = A[i] @ ETA @ A[j] - A[j] @ ETA @ A[i]
                ud += wt * 4.0 * np.einsum(
                    "...ab,...cd,ac,bd->...", F, F, ETA, ETA)
    Me = M4 @ ETA
    P = Me
    t = []
    for p in range(4):
        if p:
            P = P @ Me
        t.append(np.einsum("...kk->...", P))
    cp = g4["c4_of"](cfg)
    vd = g4["W1"] * sum((t[p] - cp[p]) ** 2 for p in range(4))
    n = M4.shape[0]
    return (ud + vd)[n // 2:, n // 2]


def half(M4):
    n = M4.shape[0]
    return M4[n // 2:, n // 2]


def main(path, tag, mode="period", gen_name="clock_local"):
    Z = np.load(path)
    M = Z["M"].astype(np.float64)
    s = float(Z["s"]) if "s" in Z.files else 1.0
    cfg = g4["base_cfg"](s=s, delta=float(Z["delta"]),
                         n=M.shape[0], L=M.shape[0] * float(Z["h"]))
    states, dens = [], []
    if mode == "period":
        a0s = g4["gen_catalog"](cfg, M)
        # the orbit: global generators exponentiate exactly; local
        # generators rotate per-cell about the local axis
        lam, V = np.linalg.eigh(M[..., 1:4, 1:4])
        for kk in range(7):
            th = 2.0 * np.pi * kk / 6.0
            if gen_name in ("rot_z", "rot_x", "boost_z", "boost_x"):
                Gm = {"rot_z": (1, 2), "rot_x": (2, 3)}.get(gen_name)
                W = np.zeros((4, 4))
                if Gm:
                    W[Gm[0], Gm[1]], W[Gm[1], Gm[0]] = -1.0, 1.0
                else:
                    i = {"boost_z": 3, "boost_x": 1}[gen_name]
                    W[0, i] = W[i, 0] = 1.0
                L = expm(th * W)
                Mk = np.einsum("ab,...bc,dc->...ad", L, M, L)
            else:
                axis = V[..., :, 2] if gen_name == "clock_local" \
                    else V[..., :, 0]
                W = g4["gen_catalog"].__globals__  # noqa: F841
                # per-cell Rodrigues rotation of the spatial block
                K = np.zeros(M.shape[:3] + (3, 3))
                n1, n2, n3 = axis[..., 0], axis[..., 1], axis[..., 2]
                K[..., 0, 1], K[..., 0, 2] = -n3, n2
                K[..., 1, 0], K[..., 1, 2] = n3, -n1
                K[..., 2, 0], K[..., 2, 1] = -n2, n1
                R = (np.eye(3) + np.sin(th) * K
                     + (1 - np.cos(th)) * (K @ K))
                Mk = M.copy()
                Mk[..., 1:4, 1:4] = np.einsum(
                    "...ab,...bc,...dc->...ad", R, M[..., 1:4, 1:4], R)
            states.append({"it": kk, "t": th, "M": half(Mk)})
            dens.append(dens_slice(Mk, cfg))
        sup = (f"M5.21.3 {tag}: ONE internal-clock period "
               f"({gen_name}; frame axis = CLOCK PHASE 0..2pi, "
               f"y = 0 half-plane, true 4x4 field)")
    else:
        snaps = sorted(int(k[len("M_it"):]) for k in Z.files
                       if k.startswith("M_it"))
        seq = [(0, Z["M_it0"] if "M_it0" in Z.files else None)]
        seq = []
        for it in snaps:
            seq.append((it, Z[f"M_it{it}"].astype(np.float64)))
        seq.append((int(Z["maxit"]) if "maxit" in Z.files else -1, M))
        for it, Mk in seq:
            states.append({"it": it, "t": float(it), "M": half(Mk)})
            dens.append(dens_slice(Mk, cfg))
        sup = (f"M5.21.3 {tag}: descent (y = 0 half-plane, "
               f"true 4x4 field)")
    it = iter(dens)

    def density_fn(_):
        return next(it)
    film_strip(states, os.path.join(
        PLOTS, f"m5_21_3_film_basic_{tag}.png"),
        template="basic", delta=cfg["delta"], h=cfg["h"],
        g=cfg["g"], wscale=g4["W1"], density_fn=density_fn,
        suptitle=sup + " on the TRUE basic template")
    film_strip(states, os.path.join(
        PLOTS, f"m5_21_3_film_thermal_{tag}.png"),
        template="thermal", delta=cfg["delta"], h=cfg["h"],
        g=cfg["g"], suptitle=sup + " on the TRUE thermal template")
    print("films written for", tag)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else "run",
         sys.argv[3] if len(sys.argv) > 3 else "period",
         sys.argv[4] if len(sys.argv) > 4 else "clock_local")
