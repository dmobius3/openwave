"""M5.21.2 stencil-artifact quantification (the M5.21.1e anti-recipe,
re-fired here by the film read: interior checkerboard from it ~2000).

The 2h central-difference curvature stencil has a null mode (odd-even
eigenvector alternation costs zero measured u). This script re-reads
every saved endpoint with a COMPACT 1h forward-difference stencil,
which the sawtooth cannot hide from:
    sawtooth ratio  rho_st = ||D_fwd M|| / ||D_2h M||   (~1 smooth,
        >> 1 checkerboarded)
    E_u^fwd vs E_u^2h : the artifact's hidden curvature load
V is pointwise (stencil-free) and is not re-read.

Run: python3 m5_21_2_c_stencil_check.py
Out: ../data/m5_21_2_stencil.json (+ printed table)
"""
import glob
import json
import os
import runpy

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
g = runpy.run_path(os.path.join(HERE, "m5_21_2_a_scan3d.py"),
                   run_name="not_main")


def d_fwd(f, ax):
    sl = [slice(None)] * f.ndim

    def at(i):
        s = list(sl); s[ax] = i; return tuple(s)
    out = np.zeros_like(f)
    out[at(slice(0, -1))] = f[at(slice(1, None))] - f[at(slice(0, -1))]
    return out


def e_u_fwd(M):
    A = [d_fwd(M, ax) for ax in range(3)]
    e = 0.0
    for i in range(3):
        for j in range(i + 1, 3):
            C = A[i] @ A[j] - A[j] @ A[i]
            e += 4.0 * np.sum(np.einsum("...kl,...kl->...", C, C))
    return float(e)


def sawtooth(M):
    num, den = 0.0, 0.0
    for ax in range(3):
        num += float(np.sum(d_fwd(M, ax) ** 2))
        den += float(np.sum(g["d_ax"](M, ax) ** 2))
    return float(np.sqrt(num / max(den, 1e-300)))


rows = {}
files = sorted(glob.glob(os.path.join(DATA, "m5_21_2_end_*.npz")))
seed = g["seed3"](48, 0.3, "A")
e2h_s, _ = g["e_split"](seed, 0.3)
rows["seed_A_baseline"] = {"E_u_2h": float(e2h_s),
                           "E_u_fwd": e_u_fwd(seed),
                           "sawtooth": sawtooth(seed)}
for p in files:
    k = os.path.basename(p)[len("m5_21_2_end_"):-len(".npz")]
    M = np.load(p)["M"].astype(np.float64)
    e2h, ev = g["e_split"](M, 0.3)
    rows[k] = {"E_u_2h": float(e2h), "E_u_fwd": e_u_fwd(M),
               "E_v_pointwise": float(ev), "sawtooth": sawtooth(M)}
with open(os.path.join(DATA, "m5_21_2_stencil.json"), "w") as f:
    json.dump(rows, f, indent=1)
print(f"{'endpoint':>18} {'E_u_2h':>10} {'E_u_fwd':>10} "
      f"{'ratio':>7} {'sawtooth':>9}")
for k, r in rows.items():
    print(f"{k:>18} {r['E_u_2h']:>10.3f} {r['E_u_fwd']:>10.3f} "
          f"{r['E_u_fwd'] / max(r['E_u_2h'], 1e-300):>7.2f} "
          f"{r['sawtooth']:>9.3f}")
