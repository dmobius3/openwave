"""M5.22.4 panels.

  fullf   the add-on panel: per-state per-read flux bars + the
          identity diff (comp3 vs basic)
  ladder  E*(omega) per state (clock_local), against the matched
          static control; per-rung q_far annotated

Out: ../plots/m5_22_4_{fullf,ladder}.png
"""
from __future__ import annotations

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

KEYS = ["prot", "neut", "d2", "deut"]
LABEL = {"prot": "proton-analog", "neut": "neutron-analog",
         "d2": "d2 neutral", "deut": "deuteron cand."}
TAGS = {
    "prot": "P-0.5_plane_sc6_n32_pinned_d0.3",
    "neut": "P-1_plane_sc6_n32_pinned_d0.3",
    "d2": "d2_s-0.5_s-0.5_a2_sc6_n32_d0.3",
    "deut": "dn_-0.5at-2_+0.5at+0_-0.5at+2_n32_d0.3",
}


def fullf():
    with open(os.path.join(DATA, "m5_22_4_fullf_all.json")) as f:
        allr = json.load(f)
    reads = ["basic", "comp3", "comp2", "comp1", "norm3"]
    cols = ["0.3", "tab:blue", "tab:orange", "tab:green",
            "tab:red"]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    ax = axes[0]
    x = np.arange(len(KEYS))
    wdt = 0.16
    for i, rd in enumerate(reads):
        v = [allr[TAGS[k]][rd]["half18"] for k in KEYS]
        ax.bar(x + (i - 2) * wdt, v, wdt, label=rd, color=cols[i])
    ax.axhline(0, color="k", lw=0.6)
    for q in (-1, 1):
        ax.axhline(q, color="0.8", lw=0.6, ls="--")
    ax.set_xticks(x, [LABEL[k] for k in KEYS], fontsize=8)
    ax.set_ylabel("Gauss flux (half = 18)")
    ax.legend(fontsize=7)
    ax.set_title("full-F reads vs the basic instrument", fontsize=9)
    ax = axes[1]
    dif = [allr[TAGS[k]]["identity_diff"]["max_abs"] for k in KEYS]
    rms = [allr[TAGS[k]]["identity_diff"]["rms"] for k in KEYS]
    brms = [allr[TAGS[k]]["identity_diff"]["basic_rms"]
            for k in KEYS]
    ax.semilogy(x, dif, "o-", label="|comp3 - basic| max")
    ax.semilogy(x, rms, "s-", label="|comp3 - basic| rms")
    ax.semilogy(x, brms, "^--", color="0.5", label="basic rms (scale)")
    ax.set_xticks(x, [LABEL[k] for k in KEYS], fontsize=8)
    ax.legend(fontsize=7)
    ax.set_title("the identity: full-F long-axis component = basic",
                 fontsize=9)
    fig.tight_layout()
    out = os.path.join(PLOTS, "m5_22_4_fullf.png")
    fig.savefig(out, dpi=140)
    print(out)


def ladder():
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5), sharex=True)
    for k, ax in zip(KEYS, axes.ravel()):
        p = os.path.join(DATA, f"m5_22_4_row_p3_{k}_clock_local.json")
        if not os.path.exists(p):
            ax.set_title(f"{LABEL[k]}: (no ladder)")
            continue
        with open(p) as f:
            row = json.load(f)
        lad = row["ladder"]
        om = [r["omega"] for r in lad]
        E = [r["E"] for r in lad]
        ax.plot(om, E, "o-", label="E*(omega) clock_local")
        cp = os.path.join(DATA, f"m5_22_4_row_ctrl_{k}.json")
        if os.path.exists(cp):
            with open(cp) as f:
                c = json.load(f)
            ax.axhline(c["E_end"], color="tab:red", ls="--",
                       label=f"static ctrl ({c['maxit']} it)")
        ax.axhline(row["E_static"], color="0.6", ls=":",
                   label="E_static (P1)")
        for r in lad[1:]:
            if r.get("q_far") is not None:
                ax.annotate(f"q {r['q_far']:+.2f}",
                            (r["omega"], r["E"]), fontsize=6,
                            textcoords="offset points",
                            xytext=(4, 4))
        ax.set_title(LABEL[k], fontsize=9)
        ax.legend(fontsize=7)
    for ax in axes[1]:
        ax.set_xlabel("omega")
    for ax in axes[:, 0]:
        ax.set_ylabel("E")
    fig.suptitle("M5.22.4 the omega-twist ladder, four baryon-sector "
                 "states (eta-signed functional)", fontsize=10)
    fig.tight_layout()
    out = os.path.join(PLOTS, "m5_22_4_ladder.png")
    fig.savefig(out, dpi=140)
    print(out)


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "fullf"
    if mode == "fullf":
        fullf()
    elif mode == "ladder":
        ladder()
