"""M5.21.2e: the TRUE-template film adapter for 3D states (the film
compliance fix, user catch at the M5.21.2 review).

The m5_film.py basic/thermal templates consume axisym (nr, nz, 4, 4)
states. This adapter renders a 3D cube state on them: the y = 0
meridional HALF-plane slice (rho = x >= 0, z), the 3x3 spatial block
embedded 4x4 block-diag with the display-constant time eigenvalue g,
and the energy density panel fed the TRUE 3D density (our fwd-stencil
u + V on the cube, sliced) via the density_fn hook, never the axisym
4D density.

Run:  python3 m5_21_2_e_film_adapter.py [endpoint_npz] [label]
Out:  ../plots/m5_21_2_film_basic_<label>_template.png (+ _thermal_)
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
os.environ.setdefault("M5212_STENCIL", "fwd")
g = runpy.run_path(os.path.join(HERE, "m5_21_2_a_scan3d.py"),
                   run_name="not_main")
from m5_film import film_strip  # noqa: E402

DELTA = 0.3
G_DISP = 8.0                     # display-constant time eigenvalue


def to_axisym4(M3):
    """y = 0 half-plane slice -> (nr, nz, 4, 4), rho = x >= 0."""
    n = M3.shape[0]
    j = n // 2
    half = M3[n // 2:, j]                     # (nr, nz, 3, 3)
    nr, nz = half.shape[:2]
    M4 = np.zeros((nr, nz, 4, 4))
    M4[..., 1:4, 1:4] = half
    M4[..., 0, 0] = G_DISP
    return M4


def dens3_slice(M3):
    """true 3D density (fwd stencil), y = 0 half-plane slice."""
    A = [g["d_ax"](M3, ax) for ax in range(3)]
    ud = np.zeros(M3.shape[:3])
    for i in range(3):
        for k in range(i + 1, 3):
            C = A[i] @ A[k] - A[k] @ A[i]
            ud += 4.0 * np.einsum("...kl,...kl->...", C, C)
    M2 = M3 @ M3
    t1 = np.einsum("...kk->...", M3)
    t2 = np.einsum("...kk->...", M2)
    t3 = np.einsum("...kk->...", M2 @ M3)
    cp = g["c_p"](DELTA)
    vd = g["WSCALE"] * ((t1 - cp[0]) ** 2 + (t2 - cp[1]) ** 2
                        + (t3 - cp[2]) ** 2)
    n = M3.shape[0]
    return (ud + vd)[n // 2:, n // 2]


def main(path, label):
    Mend = np.load(path)["M"].astype(np.float64)
    n = Mend.shape[0]
    seed = g["seed3"](n, DELTA, "A") if label.startswith("A") else \
        g["seed_ring"](n, DELTA) if label.startswith("R") else \
        g["seed3"](n, DELTA, label[0])
    states3 = [{"it": 0, "t": 0.0, "M3": seed},
               {"it": 8000, "t": 8000.0, "M3": Mend}]
    dens = [dens3_slice(st["M3"]) for st in states3]
    states = [{"it": st["it"], "t": st["t"], "M": to_axisym4(st["M3"])}
              for st in states3]
    it = iter(dens)

    def density_fn(M4):
        return next(it)
    film_strip(states, os.path.join(
        PLOTS, f"m5_21_2_film_basic_{label}_template.png"),
        template="basic", delta=DELTA, h=1.0, g=G_DISP,
        wscale=g["WSCALE"], density_fn=density_fn,
        suptitle=f"M5.21.2 {label} (3D, fwd instrument) on the "
                 f"TRUE basic template (y = 0 half-plane; time row "
                 f"= display constant g)")
    print("basic template film written for", label)


if __name__ == "__main__":
    p = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        DATA, "m5_21_2_end_A_pinned_fwd.npz")
    lab = sys.argv[2] if len(sys.argv) > 2 else "A_pinned_fwd"
    main(p, lab)
