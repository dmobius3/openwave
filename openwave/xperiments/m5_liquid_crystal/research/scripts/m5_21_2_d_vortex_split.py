"""M5.21.2d: the vortex-split read (Duda 2026-07-17 11:05 question).

HIS QUESTION: "does the biaxial hedgehog ansatz with degree 1
topological vortex (in cross-section) split into two vortices of 1/2?
... Could you tell more about structure of topological vortices
finally going out of minimizing energy hedgehog - how many, of what
degree, in what angles?"

METHOD. On xy cross-section planes (perpendicular to the seeded vortex
axis z) of the fwd-instrument seed-A endpoint vs the seed:
    per cell: eigh -> ascending (lam_min, lam_mid, lam_max).
    transverse-frame angle phi_m = angle of the MID eigenvector's xy
        projection, mod pi (director).
    plaquette winding: sum of minimal mod-pi differences around each
        unit plaquette; a +-pi plaquette = a +-1/2 disclination core.
    gap maps (lam_mid - lam_min), (lam_max - lam_mid): core
        localization + the near-degeneracy caveat regions.
Totals per plane: core count, positions (rho, polar angle), net
winding inside rho < 12. Seed baseline expectation: the degree-1 axis
vortex = 2 half-units co-located at the axis (net winding 2pi of the
mod-pi director).

Run:  python3 m5_21_2_d_vortex_split.py [endpoint_npz]
Out:  ../data/m5_21_2_split.json, ../plots/m5_21_2_split.png
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
g = runpy.run_path(os.path.join(HERE, "m5_21_2_a_scan3d.py"),
                   run_name="not_main")

DELTA = 0.3
N = 48
PLANES = (4, 8, 12, 16)          # z offsets from center (axis cells)
RHO_IN = 12.0                     # winding-accounting radius


def plane_analysis(M, kz):
    """analyze the xy-plane at z-offset kz (cells); M (N,N,N,3,3)."""
    c = (N - 1) / 2.0
    k = int(round(c + kz))
    P = M[:, :, k]
    lam, V = np.linalg.eigh(P)
    gap_lo = lam[..., 1] - lam[..., 0]
    gap_hi = lam[..., 2] - lam[..., 1]
    vm = V[..., :, 1]                       # mid eigenvector
    phi = np.arctan2(vm[..., 1], vm[..., 0])   # mod pi below
    x = np.arange(N) - c
    X, Y = np.meshgrid(x, x, indexing="ij")
    rho = np.hypot(X, Y)

    def dmodpi(a):
        return (a + np.pi / 2) % np.pi - np.pi / 2
    w = (dmodpi(phi[1:, :-1] - phi[:-1, :-1])
         + dmodpi(phi[1:, 1:] - phi[1:, :-1])
         + dmodpi(phi[:-1, 1:] - phi[1:, 1:])
         + dmodpi(phi[:-1, :-1] - phi[:-1, 1:]))
    half = np.round(w / np.pi).astype(int)     # winding in half-units
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
                      "clean": bool(gmin > 0.05)})
    net = int(half[inside].sum())
    n_clean = sum(1 for cr in cores if cr["clean"])
    return {"kz": int(kz), "n_cores": len(cores), "n_clean": n_clean,
            "net_half_units": net, "cores": cores}, gap_lo, phi, half


def run(path, tag):
    M = np.load(path)["M"].astype(np.float64)
    seed = g["seed3"](N, DELTA, "A").astype(np.float32).astype(
        np.float64)                     # float32 parity with the npz
    out = {"endpoint": [], "seed": []}
    for kz in PLANES:
        r_end, _, _, _ = plane_analysis(M, kz)
        r_seed, _, _, _ = plane_analysis(seed, kz)
        out["endpoint"].append(r_end)
        out["seed"].append(r_seed)
    # the figure: the kz = 8 plane, seed vs endpoint
    fig, axs = plt.subplots(2, 3, figsize=(13.5, 8.6))
    for row, (Mf, lab) in enumerate(((seed, "seed"), (M, tag))):
        res, gap_lo, phi, half = plane_analysis(Mf, 8)
        ax = axs[row, 0]
        im = ax.imshow(gap_lo.T, origin="lower", cmap="viridis",
                       vmin=0, vmax=0.4)
        ax.set_title(f"{lab}: gap (lam_mid - lam_min), z = +8",
                     fontsize=9)
        fig.colorbar(im, ax=ax, fraction=0.045)
        ax = axs[row, 1]
        im = ax.imshow((np.mod(phi, np.pi)).T, origin="lower",
                       cmap="twilight", vmin=0, vmax=np.pi)
        ax.set_title(f"{lab}: mid-eigenvector angle mod pi",
                     fontsize=9)
        fig.colorbar(im, ax=ax, fraction=0.045)
        ax = axs[row, 2]
        im = ax.imshow(half.T, origin="lower", cmap="bwr", vmin=-2,
                       vmax=2)
        for cr in res["cores"]:
            ax.plot(cr["xy"][0] + (N - 1) / 2,
                    cr["xy"][1] + (N - 1) / 2, "o", ms=10, mfc="none",
                    mec="k")
        ax.set_title(f"{lab}: plaquette winding (half-units); "
                     f"cores = {res['n_cores']}, net = "
                     f"{res['net_half_units']}/2", fontsize=9)
        fig.colorbar(im, ax=ax, fraction=0.045)
    fig.suptitle("M5.21.2d: transverse-vortex structure, seed vs "
                 "relaxed endpoint (fwd instrument, seed A pinned)",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS, "m5_21_2_split.png"), dpi=120)
    with open(os.path.join(DATA, "m5_21_2_split.json"), "w") as f:
        json.dump(out, f, indent=1)
    for side in ("seed", "endpoint"):
        for r in out[side]:
            print(side, "kz", r["kz"], "cores", r["n_cores"],
                  f"(clean {r['n_clean']})", "net",
                  f"{r['net_half_units']}/2",
                  [(round(c['rho'], 1), round(c['angle_deg']),
                    c['half_units'], "c" if c['clean'] else "x")
                   for c in r["cores"]][:8])


if __name__ == "__main__":
    p = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        DATA, "m5_21_2_end_A_pinned_fwd.npz")
    run(p, "endpoint")
