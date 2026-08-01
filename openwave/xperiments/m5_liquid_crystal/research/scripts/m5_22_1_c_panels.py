"""M5.22.1 figures: opening-move panels.

  slices  4-panel meridional portrait of any saved endpoint
          (m5_22 census tags or m5_22_1 task tags): (a) director over
          eigengap, (b) log eigengap, (c) log |Mermin-Ho B|,
          (d) SIGNED topological charge density div B / 4pi (RdBu).
          Panel (d) is labeled with the author's 2026-07-30 hedgehog
          convention: the ELECTRIC reading is the NEGATED value
          (census note § 7).
  ladder  the delta ladder: E(delta) per state + the n/p and p/e
          ratios from the census + M5.22.1 rows.
  kick    the kick-apart summary: ring z trajectories (dynamic
          branches) + endpoint E vs the parent.

Usage:
  python3 m5_22_1_c_panels.py slices <tag> [src=m5_22|m5_22_1] ...
  python3 m5_22_1_c_panels.py ladder
  python3 m5_22_1_c_panels.py kick
Out: ../plots/m5_22_1_*.png
"""
from __future__ import annotations

import glob
import importlib.util
import json
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


def slices(tag, src="m5_22_1"):
    z = np.load(os.path.join(DATA, f"{src}_end_{tag}.npz"))
    M = z["M"].astype(np.float64)
    n, h = int(z["n"]), float(z["h"])
    j = n // 2
    ax_coord = (np.arange(n) - (n - 1) / 2.0) * h
    lam, vec = np.linalg.eigh(M)
    gap = np.minimum(lam[..., 1] - lam[..., 0],
                     lam[..., 2] - lam[..., 1])
    nhat, _ = PAIR.orient_v1(M)
    B = PAIR.mermin_B(nhat, h)
    Bmag = np.sqrt(np.sum(B * B, axis=-1))
    div = np.zeros(B.shape[:3])
    for ax in range(3):
        sl = [slice(None)] * 3
        sp_, sm = list(sl), list(sl)
        sl[ax] = slice(1, -1)
        sp_[ax] = slice(2, None)
        sm[ax] = slice(None, -2)
        d = np.zeros(B.shape[:3])
        d[tuple(sl)] = (B[tuple(sp_) + (ax,)]
                        - B[tuple(sm) + (ax,)]) / (2 * h)
        div += d
    rho = div / (4.0 * np.pi)
    v1 = vec[..., :, 2]
    fig, axs = plt.subplots(1, 4, figsize=(19, 4.6))
    ex = [ax_coord[0], ax_coord[-1], ax_coord[0], ax_coord[-1]]
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
    a0.set_title(f"(a) director (X-Z) over eigengap: {tag}",
                 fontsize=9)
    a0.set_xlabel("X")
    a0.set_ylabel("Z (rotation axis)")
    im1 = axs[1].imshow(np.log10(np.maximum(gap[:, j, :].T, 1e-6)),
                        origin="lower", extent=ex, cmap="magma",
                        aspect="equal")
    axs[1].set_title("(b) log10 min eigengap (cores dark)", fontsize=9)
    plt.colorbar(im1, ax=axs[1], shrink=0.8)
    im2 = axs[2].imshow(np.log10(np.maximum(Bmag[:, j, :].T, 1e-8)),
                        origin="lower", extent=ex, cmap="cividis",
                        aspect="equal")
    axs[2].set_title("(c) log10 |Mermin-Ho B|", fontsize=9)
    plt.colorbar(im2, ax=axs[2], shrink=0.8)
    r = rho[:, j, :].T
    vmax = max(1e-6, float(np.abs(r).max()))
    im3 = axs[3].imshow(r, origin="lower", extent=ex, cmap="RdBu_r",
                        vmin=-vmax, vmax=vmax, aspect="equal")
    axs[3].set_title("(d) TOPOLOGICAL charge density div B / 4pi\n"
                     "(electric reading = NEGATED, note § 7)",
                     fontsize=9)
    plt.colorbar(im3, ax=axs[3], shrink=0.8)
    fig.tight_layout()
    p = os.path.join(PLOTS, f"m5_22_1_slice_{tag}.png")
    fig.savefig(p, dpi=110)
    plt.close(fig)
    print(f"saved {p}")


def _row(src, tag):
    p = os.path.join(DATA, f"{src}_row_{tag}.json")
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return json.load(f)


def ladder():
    states = {
        "P-0.5 (proton-analog)": "P-0.5_plane_sc6_n32_pinned",
        "P-1 (neutral)": "P-1_plane_sc6_n32_pinned",
        "E-0.5 (lepton ref)": "E-0.5_plane_sc6_n32_pinned",
    }
    deltas = (0.1, 0.2, 0.3)
    fig, (a0, a1) = plt.subplots(1, 2, figsize=(11, 4.4))
    table = {}
    for label, stem in states.items():
        es = []
        for d in deltas:
            r = _row("m5_22", f"{stem}_d{d:g}")
            es.append(r["E_end"] if r else np.nan)
        table[label] = es
        a0.plot(deltas, es, "o-", label=label)
    a0.set_xlabel("delta (biaxiality)")
    a0.set_ylabel("E_end (n = 32, pinned)")
    a0.set_title("the delta ladder: headline-state energies",
                 fontsize=10)
    a0.legend(fontsize=8)
    a0.invert_xaxis()
    ks = list(table)
    ratio_np = [table[ks[1]][i] / table[ks[0]][i]
                for i in range(len(deltas))]
    ratio_pe = [table[ks[0]][i] / table[ks[2]][i]
                for i in range(len(deltas))]
    a1.plot(deltas, ratio_np, "s-", label="E(neutral) / E(proton)")
    a1.plot(deltas, ratio_pe, "d-", label="E(proton) / E(lepton)")
    a1.axhline(1.0, color="k", lw=0.6)
    a1.set_xlabel("delta")
    a1.set_ylabel("ratio")
    a1.set_title("mass-ordering ratios vs delta (physical: "
                 "n/p 1.0014, p/e 1836)", fontsize=10)
    a1.legend(fontsize=8)
    a1.invert_xaxis()
    fig.tight_layout()
    p = os.path.join(PLOTS, "m5_22_1_delta_ladder.png")
    fig.savefig(p, dpi=110)
    plt.close(fig)
    print(f"saved {p}")
    print(json.dumps({"E": table, "ratio_neutral_over_proton":
                      ratio_np, "ratio_proton_over_lepton": ratio_pe}))


def kick():
    parent_tag = "P-1_plane_sc6_n32_pinned_d0.3"
    parent = _row("m5_22", parent_tag)
    rows = []
    for p in sorted(glob.glob(os.path.join(
            DATA, f"m5_22_1_row_{parent_tag}_*.json"))):
        with open(p) as f:
            rows.append(json.load(f))
    fig, (a0, a1) = plt.subplots(1, 2, figsize=(12, 4.4))
    for r in rows:
        if r["branch"] != "kick" or "hist" not in r:
            continue
        its = [s["it"] for s in r["hist"]]
        for k in (0, 1):
            zs = [s["ring_zs"][k] if len(s["ring_zs"]) > k else np.nan
                  for s in r["hist"]]
            a0.plot(its, zs, "o-", ms=3,
                    label=f"v={r['v']:g} ring {k}" if k == 0 else None)
    a0.set_xlabel("evolve step")
    a0.set_ylabel("ring z centroid")
    a0.set_title("dynamic kick: ring trajectories", fontsize=10)
    a0.legend(fontsize=8)
    labels, es = [], []
    for r in rows:
        lab = (f"{r['branch']} "
               + (f"d={r['d']:g}" if r["branch"] == "split"
                  else f"v={r['v']:g}"))
        labels.append(lab)
        es.append(r["E_end"])
    a1.bar(range(len(es)), es, color="tab:blue")
    a1.axhline(parent["E_end"], color="k", ls="--", lw=1,
               label=f"parent E = {parent['E_end']:.3f}")
    a1.set_xticks(range(len(labels)))
    a1.set_xticklabels(labels, rotation=20, fontsize=8)
    a1.set_ylabel("E_end after re-relax")
    a1.set_title("kick-apart branches: endpoint energies", fontsize=10)
    a1.legend(fontsize=8)
    fig.tight_layout()
    p = os.path.join(PLOTS, "m5_22_1_kick_panel.png")
    fig.savefig(p, dpi=110)
    plt.close(fig)
    print(f"saved {p}")


if __name__ == "__main__":
    ARGV = sys.argv[1:]
    mode = ARGV[0]
    if mode == "slices":
        args = [a for a in ARGV[1:] if not a.startswith("src=")]
        src = next((a.split("=", 1)[1] for a in ARGV[1:]
                    if a.startswith("src=")), "m5_22_1")
        for t in args:
            slices(t, src=src)
    elif mode == "ladder":
        ladder()
    elif mode == "kick":
        kick()
