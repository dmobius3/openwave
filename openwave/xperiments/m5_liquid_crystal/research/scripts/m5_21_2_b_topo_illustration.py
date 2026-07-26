"""M5.21.2 topology-catalog illustration: vortex loop vs hedgehog vs
charged ring, three cross-section rows (meridional xz, equatorial xy,
far-sphere read). Director drawn nematic-style (headless segments);
yellow = defect core (the slide convention). Out-of-plane directors in
the xy cut drawn as dots.

The mental picture: the charged ring bridges the two green cells of
the Toulouse-Kleman triangle: locally the s=1, d=1 vortex-line cell
(a closed cord, half winding around it), globally the s=2, d=0 erizo
cell (the far sphere sees q = 1). The chargeless vortex loop (the
neutrino-side object) is its inside-out partner: interior transverse,
exterior uniform, so the far sphere sees q = 0.

Out: ../plots/m5_21_2_topo_catalog.png
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
PLOTS = os.path.join(HERE, "..", "plots")
A = 5.0            # ring / loop radius
L = 11.0           # half box
YEL = "#f5c400"


def seg_panel(ax, X, Z, nx, nz, mask=None, seg=1.15, color="#1a4a8a"):
    x1 = X - 0.5 * seg * nx
    x2 = X + 0.5 * seg * nx
    z1 = Z - 0.5 * seg * nz
    z2 = Z + 0.5 * seg * nz
    segs = np.stack([np.stack([x1, z1], -1), np.stack([x2, z2], -1)],
                    axis=-2).reshape(-1, 2, 2)
    if mask is not None:
        segs = segs[mask.ravel()]
    ax.add_collection(LineCollection(segs, colors=color, lw=1.7))


def style(ax, title=None, ylab=None):
    ax.set_xlim(-L, L)
    ax.set_ylim(-L, L)
    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    if title:
        ax.set_title(title, fontsize=10)
    if ylab:
        ax.set_ylabel(ylab, fontsize=10)


def grid2(n=15):
    x = np.linspace(-L + 0.6, L - 0.6, n)
    return np.meshgrid(x, x, indexing="ij")


def psi_ring(rho, z, a=A):
    return 0.5 * (np.arctan2(rho - a, z) + np.arctan2(rho + a, z))


def psi_loop(x, z, a=A):
    """chargeless half-disclination loop, uniform (+z) far field."""
    return 0.5 * (np.arctan2(z, x - a) - np.arctan2(z, x + a))


fig, axs = plt.subplots(3, 3, figsize=(12.6, 12.9))

# ---------- column 0: the chargeless vortex loop ----------
X, Z = grid2()
ps = psi_loop(X, Z)
seg_panel(axs[0, 0], X, Z, np.sin(ps), np.cos(ps))
axs[0, 0].plot([A, -A], [0, 0], "o", color=YEL, ms=10, mec="k")
style(axs[0, 0], "vortex loop (chargeless)\nthe neutrino-side object",
      "meridional cut (xz)")
axs[0, 0].annotate("half winding\naround each cord", (A, 0),
                   (A - 8.6, -8.6), fontsize=8,
                   arrowprops=dict(arrowstyle="->", lw=0.8))

Xq, Yq = grid2()
rho = np.sqrt(Xq ** 2 + Yq ** 2)
inside = rho < A
rr = np.where(rho < 1e-9, 1e-9, rho)
seg_panel(axs[1, 0], Xq, Yq, Xq / rr, Yq / rr, mask=inside)
axs[1, 0].scatter(Xq[~inside], Yq[~inside], s=6, color="#1a4a8a")
th = np.linspace(0, 2 * np.pi, 200)
axs[1, 0].plot(A * np.cos(th), A * np.sin(th), color=YEL, lw=3)
style(axs[1, 0], None, "equatorial cut (xy)")
axs[1, 0].text(0, -9.9, "dots = director out of plane (+z)",
               ha="center", fontsize=8)

ax = axs[2, 0]
thc = np.linspace(0, 2 * np.pi, 25)[:-1]
Rc = 8.5
xs, zs = Rc * np.cos(thc), Rc * np.sin(thc)
seg_panel(ax, xs, zs, np.zeros_like(thc), np.ones_like(thc), seg=1.6)
ax.plot(Rc * np.cos(th), Rc * np.sin(th), color="0.6", lw=1, ls="--")
style(ax, None, "far-sphere read")
ax.text(0, 0, "q = 0\nuniform far field\nno Coulomb tail",
        ha="center", va="center", fontsize=10)

# ---------- column 1: the hedgehog / erizo ----------
X, Z = grid2()
r = np.sqrt(X ** 2 + Z ** 2)
rs = np.where(r < 1e-9, 1e-9, r)
seg_panel(axs[0, 1], X, Z, X / rs, Z / rs)
axs[0, 1].plot([0], [0], "o", color=YEL, ms=10, mec="k")
style(axs[0, 1], "hedgehog / erizo (q = 1)\nseed A, point core")

seg_panel(axs[1, 1], Xq, Yq, Xq / rr, Yq / rr)
axs[1, 1].plot([0], [0], "o", color=YEL, ms=10, mec="k")
style(axs[1, 1])

ax = axs[2, 1]
seg_panel(ax, xs, zs, xs / Rc, zs / Rc, seg=1.6)
ax.plot(Rc * np.cos(th), Rc * np.sin(th), color="0.6", lw=1, ls="--")
style(ax)
ax.text(0, 0, "q = 1\nradial winding\nCoulomb tail", ha="center",
        va="center", fontsize=10)

# ---------- column 2: the charged ring ----------
X, Z = grid2()
ps = psi_ring(np.abs(X), Z)
nx = np.sign(X) * np.sin(ps)
nz = np.cos(ps)
seg_panel(axs[0, 2], X, Z, nx, nz)
axs[0, 2].plot([A, -A], [0, 0], "o", color=YEL, ms=10, mec="k")
style(axs[0, 2], "charged ring (q = 1)\nthe hedgehog opened")
axs[0, 2].annotate("escaped interior\n(director = +z)", (0, 0),
                   (-8.8, 8.0), fontsize=8,
                   arrowprops=dict(arrowstyle="->", lw=0.8))
axs[0, 2].annotate("half winding\naround the cord", (-A, 0),
                   (-9.6, -9.2), fontsize=8,
                   arrowprops=dict(arrowstyle="->", lw=0.8))

seg_panel(axs[1, 2], Xq, Yq, Xq / rr, Yq / rr, mask=~inside)
axs[1, 2].scatter(Xq[inside], Yq[inside], s=6, color="#1a4a8a")
axs[1, 2].plot(A * np.cos(th), A * np.sin(th), color=YEL, lw=3)
style(axs[1, 2])
axs[1, 2].text(0, -9.9, "the loop panel inside-out", ha="center",
               fontsize=8)

ax = axs[2, 2]
seg_panel(ax, xs, zs, xs / Rc, zs / Rc, seg=1.6)
ax.plot(Rc * np.cos(th), Rc * np.sin(th), color="0.6", lw=1, ls="--")
style(ax)
ax.text(0, 0, "q = 1\nsame far sphere\nas the hedgehog", ha="center",
        va="center", fontsize=10)

fig.suptitle(
    "M5.21.2 topology catalog: loop vs hedgehog vs charged ring "
    "(yellow = core; segments = director; dots = out of plane)\n"
    "Toulouse-Kleman placement: the charged ring is LOCALLY the "
    "s=1, d=1 line cell (a closed cord) and GLOBALLY the s=2, d=0 "
    "erizo cell (far sphere q = 1)", fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.94))
os.makedirs(PLOTS, exist_ok=True)
out = os.path.join(PLOTS, "m5_21_2_topo_catalog.png")
fig.savefig(out, dpi=130)
print("wrote", out)
