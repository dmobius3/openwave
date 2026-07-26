"""M5.21.5 panel: localization, state sensitivity, g assembly, S
conventions. Reads data/m5_21_5_mu_<tag>.json + m5_21_5_bridge.json;
writes plots/m5_21_5_panel.png."""
from __future__ import annotations

import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
PLOTS = os.path.join(HERE, "..", "plots")

TAGS = ["t24A", "t32A", "t48A", "fjom0.2", "fjom0.5", "fjom1"]
COLS = {"t24A": "tab:red", "t32A": "tab:blue", "t48A": "tab:cyan",
        "fjom0.2": "tab:green", "fjom0.5": "tab:olive",
        "fjom1": "tab:purple"}


def main():
    ds = {}
    for t in TAGS:
        p = os.path.join(DATA, f"m5_21_5_mu_{t}.json")
        if os.path.exists(p):
            ds[t] = json.load(open(p))
    br = json.load(open(os.path.join(DATA, "m5_21_5_bridge.json")))

    fig, axs = plt.subplots(2, 2, figsize=(12.5, 9.5))
    ax = axs[0, 0]
    for t, d in ds.items():
        rad = d["mu_tilt12"]["radial"]
        Rs, vals = [], []
        for k, v in rad.items():
            Rs.append(float(k[1:]))
            vals.append(v["mu_clock"])
        Rs = [min(r, 24) for r in Rs]
        ax.plot(Rs, vals, "o-", color=COLS[t], label=t)
    ax.set_xlabel("R (integration radius)")
    ax.set_ylabel(r"$\mu_{clock}(<R)$ tilt12, env flow")
    ax.set_yscale("symlog", linthresh=1e-3)
    ax.set_title("localization: envelope-clock moment plateaus")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axs[0, 1]
    xs, ys = [], []
    for i, t in enumerate([t for t in TAGS if t in ds]):
        xs.append(i)
        ys.append(ds[t]["mu_tilt12"]["mu_clock"])
    ax.bar(xs, ys, color=[COLS[t] for t in TAGS if t in ds])
    ax.axhline(0.221, color="k", ls="--", lw=1,
               label="analytic seed 0.221 (24^3, rigid)")
    ax.set_xticks(xs)
    ax.set_xticklabels([t for t in TAGS if t in ds], fontsize=8)
    ax.set_yscale("log")
    ax.set_ylabel(r"$\mu_{clock}$ (tilt12, env)")
    ax.set_title("state sensitivity of the clock moment")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis="y")

    ax = axs[1, 0]
    labels, gfp, gleg = [], [], []
    for t, s in br["states"].items():
        m = s["pairings"].get("mixed", {})
        if "g_first_principles" in m:
            labels.append(t)
            gfp.append(m["g_first_principles"])
            gleg.append(m["g_legacy_4_over_alpha"])
    x = np.arange(len(labels))
    ax.bar(x - 0.18, gfp, 0.36, label="g first-principles (k=1/4pi)")
    ax.bar(x + 0.18, gleg, 0.36, label="g legacy (K=4/alpha)")
    ax.axhline(2.0023, color="k", ls="--", lw=1, label="2.0023")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_yscale("log")
    ax.set_ylabel("g")
    ax.set_title("the g assembly: no state closes at 2")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis="y")

    ax = axs[1, 1]
    labels = [t for t in TAGS if t in ds]
    sr = [ds[t]["S_twist"]["S_rigid"] for t in labels]
    se = [ds[t]["S_twist"]["S_env"] for t in labels]
    x = np.arange(len(labels))
    ax.bar(x - 0.18, sr, 0.36, label="S rigid (IR-extensive)")
    ax.bar(x + 0.18, se, 0.36, label="S env (state-of-record)")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("S (twist, unit rate)")
    ax.set_title("the S conventions")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis="y")

    fig.suptitle("M5.21.5: mu + g-factor closure under the verified"
                 " L", fontsize=13)
    fig.tight_layout()
    os.makedirs(PLOTS, exist_ok=True)
    out = os.path.join(PLOTS, "m5_21_5_panel.png")
    fig.savefig(out, dpi=130)
    print("saved", out)


if __name__ == "__main__":
    main()
