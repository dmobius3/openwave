"""M5.21.4 capture post-pass: proper d(t) + charge ledger from the
evolution snapshots (the live rows' on-axis tracker chased pin/line
artifacts and mid-run cube centers; this is the record-grade re-read).

Core positions: per half-space (z > 0 / z < 0), the centroid of the
lowest-gap connected component restricted to the axial region rho < 7.5
(cores may open into small rings; the centroid still tracks them).
Charges: FIXED cubes at the seed core positions (+/- d0/2, half 8) plus
the far cube, Mermin-Ho flux on the continuity-oriented lift.

Run:  python3 m5_21_4_e_post.py [tag]      (default cap)
Out:  ../data/m5_21_4_cap_post.json
"""
import importlib.util
import json
import os
import sys

import numpy as np
from scipy import ndimage

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
_s1 = importlib.util.spec_from_file_location(
    "pair", os.path.join(HERE, "m5_21_4_a_pair.py"))
P = importlib.util.module_from_spec(_s1)
_s1.loader.exec_module(P)
INS = P.INS


def core_locus(M, cfg):
    n, h = cfg["n"], cfg["h"]
    x = (np.arange(n) - (n - 1) / 2.0) * h
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    rho = np.sqrt(X * X + Y * Y)
    lam = np.linalg.eigvalsh(M)
    gap = np.minimum(lam[..., 2] - lam[..., 1], lam[..., 1] - lam[..., 0])
    out = {}
    for lab, sel in [("top", Z > 2 * h), ("bot", Z < -2 * h)]:
        reg = sel & (rho < 7.5) & (np.abs(Z) < 0.5 * cfg["L"] - 6 * h)
        g = np.where(reg, gap, np.inf)
        gmin = float(g.min())
        mask = g < min(gmin * 1.5 + 0.01, 0.12)
        labs, nl = ndimage.label(mask)
        if nl == 0:
            out[lab] = {"z": None, "gap": gmin}
            continue
        best, bs = None, np.inf
        for li in range(1, nl + 1):
            w = np.argwhere(labs == li)
            gm = gap[tuple(w.T)].min()
            if gm < bs:
                bs, best = gm, w
        out[lab] = {"z": float(x[best[:, 2]].mean()),
                    "rho": float(np.sqrt(x[best[:, 0]] ** 2
                                         + x[best[:, 1]] ** 2).mean()),
                    "gap": float(bs), "nvox": int(len(best))}
    return out


def main(tag="cap"):
    Z = np.load(os.path.join(DATA, f"m5_21_4_ev_{tag}.npz"))
    meta = json.load(open(os.path.join(DATA,
                                       f"m5_21_4_ev_{tag}_rows.json")))
    n, d0 = int(meta["cfg_n"]), float(meta["d0"])
    cfg = INS.base_cfg(term="T2", n=n, L=1.5 * n, bc="free",
                       stencil="sym")
    rows = []
    for key in sorted(Z.files, key=lambda k: int(k[len("M_it"):])):
        it = int(key[len("M_it"):])
        M = Z[key].astype(np.float64)
        loc = core_locus(M, cfg)
        nh, ncf = P.orient_v1(M)
        B = P.mermin_B(nh, cfg["h"])
        far0 = 0.5 * cfg["L"] - 4.0 * cfg["h"]
        row = {"it": it, "t": it * float(meta["dt"]), "locus": loc,
               "Q_top_fix": P.cube_flux(B, cfg, +d0 / 2, 8.0),
               "Q_bot_fix": P.cube_flux(B, cfg, -d0 / 2, 8.0),
               "Q_mid": P.cube_flux(B, cfg, 0.0, d0 / 2 + 8.0),
               "Q_far": P.cube_flux(B, cfg, 0.0, far0),
               "conflicts": ncf}
        zt = loc["top"].get("z")
        zb = loc["bot"].get("z")
        row["dsep"] = (zt - zb) if (zt is not None and zb is not None) \
            else None
        rows.append(row)
        print(json.dumps(row), flush=True)
    with open(os.path.join(DATA, f"m5_21_4_{tag}_post.json"), "w") as f:
        json.dump(rows, f, indent=1)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "cap")
