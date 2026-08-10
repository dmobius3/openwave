"""M5.22.2 panels: the calibration figure, the kick ladder, the
before/after structure read.

Out: ../plots/m5_22_2_calib_coulomb.png
     ../plots/m5_22_2_kick_ladder.png
     ../plots/m5_22_2_k3_before_after.png
"""
from __future__ import annotations

import glob
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


def _load(name, fname):
    sp = importlib.util.spec_from_file_location(
        name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(mod)
    return mod


DIVE = _load("dive", "m5_22_2_a_dive.py")
INS = DIVE.INS
PAIR = DIVE.PAIR
W2_T2 = DIVE.W2_T2


def fig_calib():
    cfg = INS.base_cfg(term="T2", stencil="sym", eps=0.0, w2=W2_T2,
                       n=32, delta=0.3, bc="pinned")
    n, h = cfg["n"], cfg["h"]
    X, Y, Z = INS.coords(n, h)
    r = np.sqrt(X * X + Y * Y + Z * Z)
    rs = np.where(r < 1e-9, 1e-9, r)
    nhat = np.stack([X / rs, Y / rs, Z / rs], axis=-1)
    Ef = DIVE.e_full(nhat, h)
    Ec = DIVE.e_curv(nhat, h)
    mag_f = np.linalg.norm(Ef, axis=-1).ravel()
    mag_c = np.linalg.norm(Ec, axis=-1).ravel()
    rr = r.ravel()
    sel = (rr > 1.5) & (rr < 22.0)
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    ax.loglog(rr[sel], mag_f[sel], ".", ms=2, alpha=0.25,
              color="tab:blue",
              label=r"$|E_{\rm full}|$ (dual curvature), hedgehog")
    ax.loglog(rr[sel], np.clip(mag_c[sel], 1e-9, None), ".", ms=2,
              alpha=0.25, color="tab:orange",
              label=r"$|E_{\rm curv}|=|(\hat n\cdot\nabla)\hat n|$")
    rline = np.linspace(1.5, 22.0, 100)
    ax.loglog(rline, 1.0 / rline ** 2, "k-", lw=1.5,
              label=r"Coulomb $1/r^2$")
    ax.set_xlabel("r (lattice units)")
    ax.set_ylabel("|E|")
    ax.set_title("div E calibration on the analytic hedgehog "
                 "(n = 32): the dual curvature IS Coulomb;\n"
                 "the field-line curvature is numerical zero "
                 "(flux charge 1.06 vs 0.09)")
    ax.legend(loc="lower left", fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS, "m5_22_2_calib_coulomb.png"),
                dpi=150)
    plt.close(fig)
    print("calib figure done")


def fig_ladder():
    rows = []
    for f in sorted(glob.glob(os.path.join(
            DATA, "m5_22_2_row_*.json"))):
        if f.endswith("_ext.json"):
            continue
        d = json.load(open(f))
        if "kick" not in d:
            continue
        ext = f.replace(".json", "_ext.json")
        if os.path.exists(ext):
            d["E_end"] = json.load(open(ext))["E_end"]
        rows.append(d)
    fig, ax = plt.subplots(figsize=(6.8, 4.6))
    mk = {"K1": "o", "K2": "s", "K3": "^"}
    done_labels = set()
    for d in rows:
        fam = d["kick"].split(":")[0]
        state = "P-1 (neutron-analog)" if d["start_tag"].startswith(
            "P-1") else "pp cousin"
        col = "tab:blue" if state.startswith("P-1") else "tab:red"
        lbl = f"{state}, {fam}"
        ax.plot(max(d["E_injected"], 1e-2),
                abs(d["E_end"] - d["E_start"]),
                mk[fam], color=col, ms=8, alpha=0.8,
                label=lbl if lbl not in done_labels else None)
        done_labels.add(lbl)
    ax.axhline(0.15, color="gray", ls="--", lw=1)
    ax.text(0.013, 0.17, "returned-band bound 0.15", fontsize=7,
            color="gray")
    ax.set_xscale("log")
    ax.set_yscale("symlog", linthresh=1e-3)
    ax.set_xlabel(r"$E_{\rm injected}$ (state energies: 12.73 / 15.05)")
    ax.set_ylabel(r"$|E_{\rm end}-E_{\rm start}|$ after evolve+FIRE")
    ax.set_title("The kick ladder, 18 runs: every endpoint returns\n"
                 "(K1 core-random, K2 core-twist, K3 ring-localized; "
                 "K1:0.4 shown after FIRE extension)")
    ax.legend(fontsize=7, loc="upper left")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS, "m5_22_2_kick_ladder.png"),
                dpi=150)
    plt.close(fig)
    print(f"ladder figure done ({len(rows)} runs)")


def fig_before_after():
    z0 = np.load(os.path.join(
        DATA, "m5_22_end_P-1_plane_sc6_n32_pinned_d0.3.npz"))
    z1 = np.load(os.path.join(
        DATA, "m5_22_2_end_P-1_plane_sc6_n32_pinned_d0.3"
        "__K30.4_s216.npz"))
    cfg = INS.base_cfg(term="T2", stencil="sym", eps=0.0, w2=W2_T2,
                       n=int(z0["n"]), delta=float(z0["delta"]),
                       bc="pinned")
    n, h = cfg["n"], cfg["h"]
    ext = [-(n - 1) / 2 * h, (n - 1) / 2 * h] * 2
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.4))
    for ax, z, ttl in [
            (axes[0], z0, "P-1 start (E = 12.730)"),
            (axes[1], z1, "after K3:0.4 (E_inj = 671) + evolve + "
             "FIRE (E = 12.749)")]:
        M = z["M"].astype(np.float64)
        lam = np.linalg.eigvalsh(M)
        gap = np.minimum(lam[..., 1] - lam[..., 0],
                         lam[..., 2] - lam[..., 1])
        im = ax.imshow(gap[:, n // 2, :].T, origin="lower",
                       extent=ext, cmap="viridis", vmin=0.0,
                       vmax=0.35)
        ax.set_title(ttl, fontsize=9)
        ax.set_xlabel("x")
        ax.set_ylabel("z")
    fig.colorbar(im, ax=axes, shrink=0.85,
                 label="min eigengap (defect cores dark)")
    fig.suptitle("The ring-antiring pair survives the largest "
                 "ring-localized kick unchanged", fontsize=10)
    fig.savefig(os.path.join(PLOTS, "m5_22_2_k3_before_after.png"),
                dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("before/after figure done")


if __name__ == "__main__":
    os.makedirs(PLOTS, exist_ok=True)
    fig_calib()
    fig_ladder()
    fig_before_after()
