"""M5.22 endpoint slices: meridional-plane portraits of census states.

For each tag: the y = 0 plane (lattice X-Z), three panels:
  (a) director portrait: leading-eigenvector in-plane projection
      (line segments, director-style) over the eigengap heatmap
  (b) min eigengap (log color): cores = dark
  (c) Mermin-Ho |B| of the oriented lift: charge density location
Usage: python3 m5_22_d_slices.py tag1 [tag2 ...]
Out: ../plots/m5_22_slice_<tag>.png
"""
from __future__ import annotations

import importlib.util
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
PLOTS = os.path.join(HERE, "..", "plots")


def _load(name, fname):
    sp = importlib.util.spec_from_file_location(
        name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(mod)
    return mod


INS = _load("ins", "m5_21_2b_a_instrument.py")
PAIR = _load("pair", "m5_21_4_a_pair.py")


def slices(tag):
    z = np.load(os.path.join(DATA, f"m5_22_end_{tag}.npz"))
    M = z["M"].astype(np.float64)
    n, h = int(z["n"]), float(z["h"])
    j = n // 2                       # y = 0 plane (closest row)
    ax_coord = (np.arange(n) - (n - 1) / 2.0) * h
    lam, vec = np.linalg.eigh(M)
    gap = np.minimum(lam[..., 1] - lam[..., 0],
                     lam[..., 2] - lam[..., 1])
    nhat, _ = PAIR.orient_v1(M)
    B = PAIR.mermin_B(nhat, h)
    Bmag = np.sqrt(np.sum(B * B, axis=-1))
    v1 = vec[..., :, 2]
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.6))
    ex = [ax_coord[0], ax_coord[-1], ax_coord[0], ax_coord[-1]]
    # (a) director portrait over gap
    a0 = axs[0]
    a0.imshow(gap[:, j, :].T, origin="lower", extent=ex,
              cmap="viridis", aspect="equal")
    st = max(1, n // 24)
    for i in range(0, n, st):
        for k in range(0, n, st):
            vx, vz = v1[i, j, k, 0], v1[i, j, k, 2]
            L = 0.42 * st * h
            a0.plot([ax_coord[i] - L * vx, ax_coord[i] + L * vx],
                    [ax_coord[k] - L * vz, ax_coord[k] + L * vz],
                    "w-", lw=0.7)
    a0.set_title(f"(a) director (X-Z plane) over eigengap: {tag}",
                 fontsize=9)
    a0.set_xlabel("X")
    a0.set_ylabel("Z (rotation axis)")
    # (b) gap log
    im1 = axs[1].imshow(np.log10(np.maximum(gap[:, j, :].T, 1e-6)),
                        origin="lower", extent=ex, cmap="magma",
                        aspect="equal")
    axs[1].set_title("(b) log10 min eigengap (cores dark)", fontsize=9)
    plt.colorbar(im1, ax=axs[1], shrink=0.8)
    # (c) |B|
    im2 = axs[2].imshow(np.log10(np.maximum(Bmag[:, j, :].T, 1e-8)),
                        origin="lower", extent=ex, cmap="cividis",
                        aspect="equal")
    axs[2].set_title("(c) log10 |Mermin-Ho B| (charge density)",
                     fontsize=9)
    plt.colorbar(im2, ax=axs[2], shrink=0.8)
    fig.tight_layout()
    p = os.path.join(PLOTS, f"m5_22_slice_{tag}.png")
    fig.savefig(p, dpi=110)
    plt.close(fig)
    print(f"saved {p}")


if __name__ == "__main__":
    for t in sys.argv[1:]:
        slices(t)
