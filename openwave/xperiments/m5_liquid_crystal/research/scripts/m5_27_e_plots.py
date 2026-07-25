"""M5.27 plots: the tongue map + the baseline/ledger panels.

Out: plots/m5_27_tongue_map.png, plots/m5_27_baseline_panel.png
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import m5_27_a_harness as H  # noqa: E402

SCAN = os.path.join(H.DATA, "m5_27_lockscan.json")
GATES = os.path.join(H.DATA, "m5_27_gates.json")
VCOL = {"SUSTAINED": "#1a9850", "DRIVEN": "#fdae61",
        "UNSTABLE": "#d73027", "NULL": "#cccccc"}


VERD = os.path.join(H.DATA, "m5_27_verdicts.json")


def tongue_map():
    d = json.load(open(SCAN))
    # verdicts come from the REFUTATION pass, not the raw scan: the raw
    # SUSTAINED threshold was frozen before the control's noise band was known
    # and turned out to sit below it (m5_27_h_refute.py)
    ref = json.load(open(VERD))
    d["rows"] = ref["rows"]
    d["refute"] = {k: v for k, v in ref.items() if k != "rows"}
    rows = d["rows"]
    eps_g = sorted({r["eps"] for r in rows})
    rat_g = sorted({round(r["ratio"], 3) for r in rows})
    gain = np.full((len(eps_g), len(rat_g)), np.nan)
    verd = np.empty((len(eps_g), len(rat_g)), dtype=object)
    for r in rows:
        i, j = eps_g.index(r["eps"]), rat_g.index(round(r["ratio"], 3))
        gain[i, j] = r["gain"]
        verd[i, j] = r["verdict"]

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 4.6))
    ax = axes[0]
    vmax = max(0.12, float(np.nanmax(np.abs(gain))))
    im = ax.imshow(gain, origin="lower", aspect="auto", cmap="RdBu_r",
                   vmin=-vmax, vmax=vmax,
                   extent=[-0.5, len(rat_g) - 0.5, -0.5, len(eps_g) - 0.5])
    ax.set_xticks(range(len(rat_g)))
    ax.set_xticklabels([f"{r:g}" for r in rat_g])
    ax.set_yticks(range(len(eps_g)))
    ax.set_yticklabels([f"{e:g}" for e in eps_g])
    ax.set_xlabel(r"$\bar\omega/\omega^*$")
    ax.set_ylabel(r"drive excursion $\epsilon=\kappa A/g$")
    ax.set_title("J-retention gain vs control (the Arnold-tongue read)")
    for i in range(len(eps_g)):
        for j in range(len(rat_g)):
            ax.text(j, i, f"{gain[i, j]:+.3f}", ha="center", va="center", fontsize=6.5)
    for x in (1.0, 2.0):
        if x in rat_g:
            ax.axvline(rat_g.index(x), color="k", ls="--", lw=1, alpha=0.55)
    plt.colorbar(im, ax=ax, label="retention gain (positive = drive sustains J)")

    ax = axes[1]
    for i, e in enumerate(eps_g):
        for j, rr in enumerate(rat_g):
            ax.scatter(rr, e, s=190, marker="s",
                       c=VCOL.get(verd[i, j], "#999999"),
                       edgecolors="k", linewidths=0.4)
    ax.set_xlabel(r"$\bar\omega/\omega^*$")
    ax.set_ylabel(r"$\epsilon$")
    ax.set_yscale("log")
    ax.set_title("verdict map after the refutation pass\n"
                 f"({ref['n_refuted']}/7 raw SUSTAINED flags refuted)")
    for x in (1.0, 2.0):
        ax.axvline(x, color="k", ls="--", lw=1, alpha=0.55)
    ax.set_ylim(eps_g[0] * 0.45, eps_g[-1] * 3.2)
    ax.text(1.0, eps_g[-1] * 1.9, "1:1", ha="center", fontsize=8)
    ax.text(2.0, eps_g[-1] * 1.9, "2:1 parametric", ha="center", fontsize=8)
    handles = [plt.Line2D([], [], marker="s", ls="", ms=9, mec="k",
                          mfc=VCOL[k], label=k) for k in VCOL]
    ax.legend(handles=handles, fontsize=7.5, loc="lower left", ncol=2)
    fig.suptitle("M5.27 phase A: prescribed background scalar $g\\to g+\\kappa\\chi$, "
                 f"$\\chi=A\\cos\\bar\\omega t$ (tag {d['tag']}, "
                 f"$\\omega^*={d['om_star']:.4f}$)", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    p = os.path.join(H.PLOTS, "m5_27_tongue_map.png")
    fig.savefig(p, dpi=140)
    print("wrote", p)
    return d


def baseline_panel(d):
    g = json.load(open(GATES))
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))

    ax = axes[0]
    jt = np.array(g["free_baseline"]["J_trace"], float)
    kt = np.array(g["free_baseline"]["kin_trace"], float)
    ax.plot(jt[:, 0], jt[:, 1] / jt[0, 1], "o-", lw=1.4, ms=3, label="carried $J/J_0$")
    ax2 = ax.twinx()
    ax2.plot(kt[:, 0], kt[:, 1], "s--", color="crimson", lw=1.2, ms=3,
             label="kinetic")
    ax2.set_ylabel("kinetic ledger", color="crimson")
    ax.set_xlabel(r"$t$ (canonical units)")
    ax.set_ylabel(r"$J/J_0$")
    ax.set_title("P0 control: free release does NOT hold the clock\n"
                 "(the M5.21.3 result on this harness)")
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(alpha=0.3)

    ax = axes[1]
    ctrl = d["control"]
    ct = np.array(ctrl["J_trace"], float)
    ax.plot(ct[:, 0], ct[:, 1] / ct[0, 1], "k-", lw=2, label="control (no drive)")
    for raw in d["raw"]:
        if raw["eps"] == max(d["eps_grid"]):
            rt = np.array(raw["J_trace"], float)
            ax.plot(rt[:, 0], rt[:, 1] / rt[0, 1], lw=1,
                    label=rf"$\bar\omega/\omega^*$={raw['om_bar']/d['om_star']:.1f}")
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$J/J_0$")
    ax.set_title(rf"driven runs at $\epsilon={max(d['eps_grid'])}$ vs control")
    ax.legend(fontsize=6.5, ncol=2)
    ax.grid(alpha=0.3)

    ax = axes[2]
    rows = d["rows"]
    for e in sorted({r["eps"] for r in rows}):
        sel = sorted([r for r in rows if r["eps"] == e], key=lambda r: r["ratio"])
        ax.plot([r["ratio"] for r in sel], [r["gain"] for r in sel],
                "o-", ms=3.5, lw=1.2, label=rf"$\epsilon$={e:g}")
    ax.axhline(0, color="k", lw=1)
    ax.axhline(0.10, color="g", ls=":", lw=1.2, label="raw threshold (0.10)")
    ax.axhspan(-d["refute"]["noise_band_retention"], d["refute"]["noise_band_retention"],
               color="grey", alpha=0.22,
               label=f"control noise band (+/-{d['refute']['noise_band_retention']:.3f})")
    for x in (1.0, 2.0):
        ax.axvline(x, color="k", ls="--", lw=1, alpha=0.5)
    ax.set_xlabel(r"$\bar\omega/\omega^*$")
    ax.set_ylabel("J-retention gain vs control")
    ax.set_title("no tongue at any registered $(\\epsilon,\\bar\\omega)$:\n"
                 "every positive gain sits inside the control's own wander")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)

    fig.tight_layout()
    p = os.path.join(H.PLOTS, "m5_27_baseline_panel.png")
    fig.savefig(p, dpi=140)
    print("wrote", p)


if __name__ == "__main__":
    d = tongue_map()
    baseline_panel(d)
