"""M5.21.6 films: decay-run states on the TRUE m5_film.py templates
(basic + thermal), the series film rule.

Same embed as m5_21_2b_c_films.py (this is its t-named twin for the
M5.21.6 npz family): y = 0 meridional half-plane slice, 3x3 embedded
4x4 block-diag with the DISPLAY-CONSTANT time eigenvalue g; the
basic template's energy panel gets the TRUE 3D density (T2 term, sym
stencil) via the density_fn hook.

Run:  python3 m5_21_6_c_films.py <npz> <tag> [seedkind] [end_it]
Out:  ../plots/m5_21_6_film_basic_<tag>.png + _thermal_
"""
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
    nr, nz = half.shape[:2]
    M4 = np.zeros((nr, nz, 4, 4))
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
        for i in range(3):
            ud += wt * cfg["eps"] * np.einsum("...kl,...kl->...",
                                              A[i], A[i])
    vd = g["v_density"](M3, cfg)
    n = M3.shape[0]
    return (ud + vd)[n // 2:, n // 2]


def main(path, tag, seedkind="A", end_it=None):
    Z = np.load(path)
    M = Z["M"].astype(np.float64)
    h = float(Z["h"]) if "h" in Z else 1.0
    delta = float(Z["delta"])
    term = str(Z["term"]) if "term" in Z.files else "T2"
    n = M.shape[0]
    cfg = g["base_cfg"](seed=seedkind, term=term, n=n, L=n * h,
                        delta=delta)
    seed = g["make_seed"](cfg)
    snap_its = sorted(int(k[len("M_it"):]) for k in Z.files
                      if k.startswith("M_it"))
    if end_it is None:
        end_it = int(Z["maxit"]) if "maxit" in Z.files else \
            (snap_its[-1] + 1 if snap_its else 1)
    states3 = [{"it": 0, "t": 0.0, "M3": seed}]
    for it in snap_its:                     # the 7-shot film standard
        states3.append({"it": it, "t": float(it),
                        "M3": Z[f"M_it{it}"].astype(np.float64)})
    states3.append({"it": end_it, "t": float(end_it), "M3": M})
    dens = [dens3_slice(st["M3"], cfg) for st in states3]
    states = [{"it": st["it"], "t": st["t"],
               "M": to_axisym4(st["M3"])} for st in states3]
    it = iter(dens)

    def density_fn(M4):
        return next(it)
    film_strip(states, os.path.join(
        PLOTS, f"m5_21_6_film_basic_{tag}.png"),
        template="basic", delta=delta, h=h, g=G_DISP,
        wscale=g["W1"], density_fn=density_fn,
        suptitle=f"M5.21.6 {tag} (3D {term} sym) on the TRUE basic "
                 f"template (y = 0 half-plane; time row = display "
                 f"constant g)")
    film_strip(states, os.path.join(
        PLOTS, f"m5_21_6_film_thermal_{tag}.png"),
        template="thermal", delta=delta, h=h, g=G_DISP,
        suptitle=f"M5.21.6 {tag} on the TRUE thermal template "
                 f"(y = 0 half-plane; time row = display constant g)")
    print("films written for", tag)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else "run",
         sys.argv[3] if len(sys.argv) > 3 else "A",
         int(sys.argv[4]) if len(sys.argv) > 4 else None)
