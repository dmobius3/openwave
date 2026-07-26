"""M5.21.2b-I5: the split-line geometry, BOTH hemispheres + angles
(Q31, Duda 2026-07-17 14:16: "four 1/4 have more sense - what are the
angles between them?" + the delta-trend ask).

Extends the M5.21.2d read (m5_21_2_d_vortex_split.py): plaquette
winding of the mid-eigenvector director angle (mod pi) on xy
cross-section planes, now at PHYSICAL z-offsets on BOTH sides of the
core (his four-line picture: two half-lines up + two down from the
point core), with per-plane core azimuths and the pairwise angle
between the two clean cores of each plane.

Per plane: eigh -> ascending; phi_m = angle of the mid eigenvector's
xy projection (mod pi); plaquette winding sum of minimal mod-pi
differences; +-pi plaquette = a +-1/2 core; gap_lo floor 0.05 flags
"clean". Physical winding-accounting radius RHO_IN.

HARDENING (2026-07-17 23:1x, the axis + band-crossing catches): the
plaquette sum is ambiguous AT the axis (the seed's co-located package
steps by exactly pi/2 across it, the mod-pi branch boundary) and
inside deep-restructured interiors the eigenframe develops
NEAR-DEGENERATE ANNULI (gap -> 1e-3) across which the band identity
reshuffles, minting fake opposite-sign cores. So the record includes
CONTOUR winding: the net mod-pi winding of the mid-band angle sampled
on circles rho_c in {6, 10, 14, 18}, with the min eigen-gap along
each contour as its reliability flag. The seed reads +2 half-units at
every (plane, radius); a pinned endpoint must read +2 at the shell.

Run:  python3 m5_21_2b_b_split.py <endpoint_npz> <tag> [seedkind]
Out:  ../data/m5_21_2b_split_<tag>.json,
      ../plots/m5_21_2b_split_<tag>.png
"""
import json
import os
import runpy
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
PLOTS = os.path.join(HERE, "..", "plots")
g = runpy.run_path(os.path.join(HERE, "m5_21_2b_a_instrument.py"),
                   run_name="not_main")

Z_PLANES = (-16.0, -12.0, -8.0, -4.0, 4.0, 8.0, 12.0, 16.0)  # physical
RHO_IN = 12.0
GAP_CLEAN = 0.05


def plane_analysis(M, h, z_phys):
    n = M.shape[0]
    c = (n - 1) / 2.0
    k = int(round(c + z_phys / h))
    k = max(0, min(n - 1, k))
    P = M[:, :, k]
    lam, V = np.linalg.eigh(P)
    gap_lo = lam[..., 1] - lam[..., 0]
    vm = V[..., :, 1]
    phi = np.arctan2(vm[..., 1], vm[..., 0])
    x = (np.arange(n) - c) * h
    X, Y = np.meshgrid(x, x, indexing="ij")
    rho = np.hypot(X, Y)

    def dmodpi(a):
        return (a + np.pi / 2) % np.pi - np.pi / 2
    w = (dmodpi(phi[1:, :-1] - phi[:-1, :-1])
         + dmodpi(phi[1:, 1:] - phi[1:, :-1])
         + dmodpi(phi[:-1, 1:] - phi[1:, 1:])
         + dmodpi(phi[:-1, :-1] - phi[:-1, 1:]))
    half = np.round(w / np.pi).astype(int)
    rho_p = 0.25 * (rho[1:, :-1] + rho[1:, 1:] + rho[:-1, 1:]
                    + rho[:-1, :-1])
    inside = rho_p < RHO_IN
    cores = []
    for (i, j) in zip(*np.nonzero((half != 0) & inside)):
        cx = 0.5 * (x[i] + x[i + 1])
        cy = 0.5 * (x[j] + x[j + 1])
        gmin = float(min(gap_lo[i, j], gap_lo[i + 1, j],
                         gap_lo[i, j + 1], gap_lo[i + 1, j + 1]))
        cores.append({"xy": [float(cx), float(cy)],
                      "rho": float(np.hypot(cx, cy)),
                      "angle_deg": float(np.degrees(
                          np.arctan2(cy, cx))),
                      "half_units": int(half[i, j]),
                      "gap_min": gmin,
                      "clean": bool(gmin > GAP_CLEAN)})
    clean = [cr for cr in cores if cr["clean"]]
    pair_angle = None
    if len(clean) == 2:
        da = abs(clean[0]["angle_deg"] - clean[1]["angle_deg"])
        pair_angle = float(min(da, 360.0 - da))
    return ({"z": float(z_phys), "kz": int(k),
             "n_cores": len(cores),
             "n_clean": len(clean),
             "net_half_units": int(half[inside].sum()),
             "pair_angle_deg": pair_angle,
             "clean_azimuths_deg": [cr["angle_deg"] for cr in clean],
             "clean_rhos": [cr["rho"] for cr in clean],
             "cores": cores},
            gap_lo, phi, half, x)


CONTOUR_RHOS = (6.0, 10.0, 14.0, 18.0)


def contour_winding(P, h, rho_c, band, nseg=720):
    """net mod-pi winding of the band's eigenvector angle on a circle
    of radius rho_c; returns (half_units, min_gap_along_contour)."""
    n = P.shape[0]
    c = (n - 1) / 2.0
    lam, V = np.linalg.eigh(P)
    th = np.linspace(0, 2 * np.pi, nseg, endpoint=False)
    ii = np.clip(np.round(c + rho_c * np.cos(th) / h).astype(int),
                 0, n - 1)
    jj = np.clip(np.round(c + rho_c * np.sin(th) / h).astype(int),
                 0, n - 1)
    vv = V[ii, jj][:, :, band]
    phi = np.arctan2(vv[:, 1], vv[:, 0])
    d = np.diff(np.concatenate([phi, phi[:1]]))
    d = (d + np.pi / 2) % np.pi - np.pi / 2
    gaps = np.minimum(lam[ii, jj, 1] - lam[ii, jj, 0],
                      lam[ii, jj, 2] - lam[ii, jj, 1])
    return int(np.round(d.sum() / np.pi)), float(gaps.min())


def contours(M, h, z_phys):
    n = M.shape[0]
    k = int(round((n - 1) / 2.0 + z_phys / h))
    k = max(0, min(n - 1, k))
    P = M[:, :, k]
    out = []
    for rho_c in CONTOUR_RHOS:
        w_half, gmin = contour_winding(P, h, rho_c, 1)
        out.append({"rho_c": rho_c, "half_units": w_half,
                    "min_gap": round(gmin, 4)})
    return out


def run(path, tag, seedkind="A"):
    Z = np.load(path)
    M = Z["M"].astype(np.float64)
    h = float(Z["h"]) if "h" in Z else 1.0
    delta = float(Z["delta"])
    n = M.shape[0]
    # float32-parity seed baseline (same shift handling as the run)
    term = str(Z["term"]) if "term" in Z else "T1"
    cfg = g["base_cfg"](seed=seedkind, term=term, n=n,
                        L=n * h, delta=delta)
    seed = g["make_seed"](cfg).astype(np.float32).astype(np.float64)
    out = {"tag": tag, "h": h, "delta": delta, "term": term,
           "endpoint": [], "seed": []}
    for z in Z_PLANES:
        r_end = plane_analysis(M, h, z)[0]
        r_seed = plane_analysis(seed, h, z)[0]
        r_end["contours"] = contours(M, h, z)
        r_seed["contours"] = contours(seed, h, z)
        out["endpoint"].append(r_end)
        out["seed"].append(r_seed)
    # hemisphere summary: mean clean azimuths above vs below
    for side, sel in (("up", lambda z: z > 0), ("down", lambda z: z < 0)):
        az = [a for r in out["endpoint"] if sel(r["z"])
              for a in r["clean_azimuths_deg"]]
        cnt = [r["n_clean"] for r in out["endpoint"] if sel(r["z"])]
        out[f"{side}_n_clean_per_plane"] = cnt
        out[f"{side}_azimuths"] = az
    # figure: two z planes (one per hemisphere), seed top row
    fig, axs = plt.subplots(2, 3, figsize=(13.5, 8.6))
    for row, zp in enumerate((8.0, -8.0)):
        res, gap_lo, phi, half, x = plane_analysis(M, h, zp)
        ax = axs[row, 0]
        im = ax.imshow(gap_lo.T, origin="lower", cmap="viridis",
                       vmin=0, vmax=0.4)
        ax.set_title(f"z = {zp:+.0f}: gap (lam_mid - lam_min)",
                     fontsize=9)
        fig.colorbar(im, ax=ax, fraction=0.045)
        ax = axs[row, 1]
        im = ax.imshow((np.mod(phi, np.pi)).T, origin="lower",
                       cmap="twilight", vmin=0, vmax=np.pi)
        ax.set_title(f"z = {zp:+.0f}: mid-eigvec angle mod pi",
                     fontsize=9)
        fig.colorbar(im, ax=ax, fraction=0.045)
        ax = axs[row, 2]
        im = ax.imshow(half.T, origin="lower", cmap="bwr", vmin=-2,
                       vmax=2)
        c = (M.shape[0] - 1) / 2.0
        for cr in res["cores"]:
            ax.plot(cr["xy"][0] / h + c, cr["xy"][1] / h + c, "o",
                    ms=10, mfc="none",
                    mec="k" if cr["clean"] else "gray")
        pa = res["pair_angle_deg"]
        ax.set_title(f"z = {zp:+.0f}: winding; clean = "
                     f"{res['n_clean']}, net = "
                     f"{res['net_half_units']}/2"
                     + (f", pair {pa:.0f} deg" if pa else ""),
                     fontsize=9)
        fig.colorbar(im, ax=ax, fraction=0.045)
    fig.suptitle(f"M5.21.2b split geometry, BOTH hemispheres: {tag} "
                 f"(delta = {delta:g})", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS, f"m5_21_2b_split_{tag}.png"),
                dpi=120)
    with open(os.path.join(DATA, f"m5_21_2b_split_{tag}.json"),
              "w") as f:
        json.dump(out, f, indent=1)
    for side in ("endpoint", "seed"):
        for r in out[side]:
            print(side, "z", f"{r['z']:+.0f}", "clean",
                  r["n_clean"], "net", f"{r['net_half_units']}/2",
                  "pair", r["pair_angle_deg"],
                  [(round(c0["rho"], 1), round(c0["angle_deg"]))
                   for c0 in r["cores"] if c0["clean"]][:6])


if __name__ == "__main__":
    p = sys.argv[1]
    tag = sys.argv[2] if len(sys.argv) > 2 else "run"
    kind = sys.argv[3] if len(sys.argv) > 3 else "A"
    run(p, tag, kind)
