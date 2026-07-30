"""M5.22 seed factory: the author's analytic cross-section families.

Source: the 2026-07-28 13:20 seed PDF (electron/positron family) + the
2026-07-29 charge-calc PDF (both baryon families re-issued, R = 1,
s in {-1,-1/2,0,1/2,1}, built-in far-field winding). Task doc:
tasks/m5_22_task_details.md § The seeding update. Mathematica
ArcTan[a, b] == atan2(b, a); step[t] = ArcTan[t]/pi + 1/2.

Families (half-plane cross-section angle ang(x, y), y >= 0):
  E  electron/positron, R = 2, s in {-1/2, +1/2}:
       ang = -s (atan2(-x, R - y) - atan2(-x, R))
  N  central-pi-step family (neutron candidate at s = -1/2 since the
     2026-07-29 charge flip), R = 1, charge law ~ -2s - 1:
       ang = s (atan2(y - R, x) - atan2(-R, x)) - pi step(5 x)
  P  fractional-step family (proton candidate at s = -1/2), R = 1,
     charge law ~ -2s (the +-pi/3 steps cancel far-field):
       ang = s (atan2(y - R, x) - atan2(-R, x)) - (2 pi/3) step(5 x)
             + (pi/3)(step(5 (x - 1)) + step(5 (x + 1)))

2D validation (the author's own instrument): Q2d = (1/pi) *
[ang(path end) - ang(path start)] along the semicircle
(x, y) = (-d cos phi, d sin phi), phi in [0, pi], d = 10, dense-sampled
+ unwrapped. GATE A: reproduce the author's printed charges to 2e-3.

3D lift (rotation around the cross-section's y = 0 axis -> lattice z):
  x_author = Z / scale, y_author = rho / scale, rho = sqrt(X^2 + Y^2)
  n1 = cos(ang) zhat + sin(ang) rhohat        (the long axis)
  second-axis convention (the F1 fork):
    'plane': n2 = -sin(ang) zhat + cos(ang) rhohat  (author: "perpendicular
             in this plane")
    'phi':   n2 = phihat                       (the M5 charged-ring
             convention, m5_21_4 _tensor_from_nhat)
  M = n1 n1^T + delta n2 n2^T                 (spectrum (1, delta, 0))
  isotropic blend (a I, a = (1+delta)/3) at the ring core
  (d_ring = sqrt((rho - R*scale)^2 + Z^2) < ~r_c) and on the axis defect
  segment (where sin(ang) != 0 at small rho):
     w = [1 - exp(-(d_ring/r_c)^2)] * [1 - sin(ang)^2 exp(-(rho/r_c)^2)]

Modes: validate | gallery
Out: ../data/m5_22_seed_charges.json, ../plots/m5_22_seed_gallery.png
"""
from __future__ import annotations

import importlib.util
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
PLOTS = os.path.join(HERE, "..", "plots")

_s1 = importlib.util.spec_from_file_location(
    "ins", os.path.join(HERE, "m5_21_2b_a_instrument.py"))
INS = importlib.util.module_from_spec(_s1)
_s1.loader.exec_module(INS)

SR = (-1.0, -0.5, 0.0, 0.5, 1.0)
R_E, R_B = 2.0, 1.0
# the author's printed far-field charges (d = 10), the GATE A targets
AUTHOR_Q = {
    "N": {-1.0: 1.01273, -0.5: 0.0127307, 0.0: -0.987269,
          0.5: -1.98727, 1.0: -2.98727},
    "P": {-1.0: 1.99991, -0.5: 0.999914, 0.0: -0.0000856704,
          0.5: -1.00009, 1.0: -2.00009},
}


def step(t):
    return np.arctan(t) / np.pi + 0.5


def ang_family(fam, s, x, y, qshift=0.0):
    """the cross-section angle in AUTHOR units. qshift (P family only,
    the Q38 arm): displace the two side steps symmetrically outward,
    x = +-1 -> +-(1 + qshift); the far-field charge is unchanged (the
    steps still cancel)."""
    if fam == "E":
        return -s * (np.arctan2(-x, R_E - y) - np.arctan2(-x, R_E))
    pair = s * (np.arctan2(y - R_B, x) - np.arctan2(-R_B, x))
    if fam == "N":
        return pair - np.pi * step(5.0 * x)
    if fam == "P":
        d = 1.0 + qshift
        return pair - (2.0 * np.pi / 3.0) * step(5.0 * x) \
            + (np.pi / 3.0) * (step(5.0 * (x - d))
                               + step(5.0 * (x + d)))
    raise ValueError(fam)


def q2d(fam, s, d=10.0, m=400001):
    """far-field winding along the semicircle, the author's instrument.
    Each atan2 term is unwrapped SEPARATELY (a continuous angle of a
    moving point; whole-angle unwrap eats genuine 2pi accumulations when
    a term crosses its branch cut)."""
    phi = np.linspace(0.0, np.pi, m)
    x, y = -d * np.cos(phi), d * np.sin(phi)
    if fam == "E":
        a = -s * (np.unwrap(np.arctan2(-x, R_E - y))
                  - np.unwrap(np.arctan2(-x, R_E * np.ones_like(y))))
    else:
        a = s * (np.unwrap(np.arctan2(y - R_B, x))
                 - np.unwrap(np.arctan2(-R_B * np.ones_like(x), x)))
        if fam == "N":
            a = a - np.pi * step(5.0 * x)
        else:
            a = a - (2.0 * np.pi / 3.0) * step(5.0 * x) \
                + (np.pi / 3.0) * (step(5.0 * (x - 1.0))
                                   + step(5.0 * (x + 1.0)))
    return float((a[-1] - a[0]) / np.pi)


def seeds_list():
    """the pinned 2026-07-29 grid: 12 seeds."""
    out = [("E", s) for s in (-0.5, 0.5)]
    out += [("N", s) for s in SR]
    out += [("P", s) for s in SR]
    return out


# ================= 3D lift =================
def seed_3d(fam, s, cfg, conv="plane", scale=6.0, r_c=2.0, qshift=0.0):
    n, h, delta = cfg["n"], cfg["h"], cfg["delta"]
    X, Y, Z = INS.coords(n, h)
    rho = np.sqrt(X * X + Y * Y)
    rhos = np.where(rho < 1e-12, 1e-12, rho)
    rhohat = np.stack([X / rhos, Y / rhos, np.zeros_like(Z)], axis=-1)
    phihat = np.stack([-Y / rhos, X / rhos, np.zeros_like(Z)], axis=-1)
    near = rho < 1e-9
    if np.any(near):
        rhohat[near] = np.array([1.0, 0.0, 0.0])
        phihat[near] = np.array([0.0, 1.0, 0.0])
    zhat = np.zeros_like(rhohat)
    zhat[..., 2] = 1.0
    ang = ang_family(fam, s, Z / scale, rho / scale, qshift=qshift)
    ca, sa = np.cos(ang), np.sin(ang)
    n1 = ca[..., None] * zhat + sa[..., None] * rhohat
    if conv == "plane":
        n2 = -sa[..., None] * zhat + ca[..., None] * rhohat
    elif conv == "phi":
        n2 = phihat
    else:
        raise ValueError(conv)
    S = (n1[..., :, None] * n1[..., None, :]
         + delta * n2[..., :, None] * n2[..., None, :])
    R_lat = (R_E if fam == "E" else R_B) * scale
    d_ring = np.sqrt((rho - R_lat) ** 2 + Z * Z)
    w = (1.0 - np.exp(-((d_ring / r_c) ** 2))) \
        * (1.0 - sa * sa * np.exp(-((rho / r_c) ** 2)))
    a = (1.0 + delta) / 3.0
    return w[..., None, None] * S \
        + (1.0 - w[..., None, None]) * (a * np.eye(3))


# ================= modes =================
def validate():
    out, worst = {}, 0.0
    for fam in ("N", "P"):
        for s in SR:
            q = q2d(fam, s)
            tgt = AUTHOR_Q[fam][s]
            err = abs(q - tgt)
            worst = max(worst, err)
            out[f"{fam}_s{s:g}"] = {"q2d": q, "author": tgt,
                                    "abs_err": err}
    for s in (-0.5, 0.5):
        out[f"E_s{s:g}"] = {"q2d": q2d("E", s), "author": None}
    out["gateA_worst_abs_err"] = worst
    out["gateA_pass"] = bool(worst < 2e-3)
    os.makedirs(DATA, exist_ok=True)
    with open(os.path.join(DATA, "m5_22_seed_charges.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))
    return out


def gallery():
    """quiver gallery of all 12 cross-sections (visual check vs the
    author's VectorPlots)."""
    slist = seeds_list()
    fig, axes = plt.subplots(3, 5, figsize=(16, 10))
    gx, gy = np.meshgrid(np.linspace(-2, 2, 25),
                         np.linspace(0, 4, 25), indexing="ij")
    row = {"E": 0, "N": 1, "P": 2}
    used = set()
    for fam, s in slist:
        col = {(-1.0): 0, (-0.5): 1, (0.0): 2, (0.5): 3, (1.0): 4}[s]
        ax = axes[row[fam]][col]
        a = ang_family(fam, s, gx, gy)
        ax.quiver(gx, gy, np.cos(a), np.sin(a), scale=35,
                  headwidth=2, headlength=3, width=2.4e-3)
        ax.set_title(f"{fam}  s = {s:g}", fontsize=9)
        ax.set_aspect("equal")
        used.add((row[fam], col))
    for i in range(3):
        for j in range(5):
            if (i, j) not in used:
                axes[i][j].axis("off")
    fig.suptitle("M5.22 seed cross-sections (author formulas, "
                 "2026-07-29 set): rotate around y = 0", fontsize=11)
    fig.tight_layout()
    os.makedirs(PLOTS, exist_ok=True)
    p = os.path.join(PLOTS, "m5_22_seed_gallery.png")
    fig.savefig(p, dpi=110)
    print(f"saved {p}")


if __name__ == "__main__":
    import sys
    mode = sys.argv[1] if len(sys.argv) > 1 else "validate"
    if mode == "validate":
        validate()
    elif mode == "gallery":
        gallery()
