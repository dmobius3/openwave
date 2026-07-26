"""M5.21.10 panel: ledger + rotation + loop census for the three
f64 dynamics windows (C mu-candidate, B tau-candidate, A control)."""
from __future__ import annotations

import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
PLOTS = os.path.join(HERE, "..", "plots")

RUNS = ("C_free64", "B_free64", "A_free64")
LBL = {"C_free64": "C (mu-candidate)", "B_free64": "B (tau-candidate)",
       "A_free64": "A (electron control)"}
COL = {"C_free64": "tab:orange", "B_free64": "tab:red",
       "A_free64": "tab:blue"}

fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
for tag in RUNS:
    ev = json.load(open(os.path.join(DATA, f"m5_21_10_ev_{tag}.json")))
    h = ev["hist"]
    t = [r["t"] for r in h]
    ax[0].plot(t, [r["E"] + r["KE"] + r["absorbed"] for r in h],
               color=COL[tag], lw=2, label=f"{LBL[tag]} E+KE+abs")
    ax[0].plot(t, [r["absorbed"] for r in h], color=COL[tag], ls="--",
               lw=1, label=f"{LBL[tag]} absorbed")
    ax[1].plot(t, [r["rot_core_deg"] for r in h], color=COL[tag],
               lw=2, label=LBL[tag])
    lp = json.load(open(os.path.join(DATA,
                                     f"m5_21_10_lp_{tag}.json")))
    ks = sorted(lp["snaps"], key=lambda k: int(k[4:]))
    tt = [lp["snaps"][k]["t"] for k in ks]
    for i, thr, ls in ((0, 0.06, ":"), (1, 0.09, "-"), (2, 0.15, "--")):
        ax[2].plot(tt, [lp["snaps"][k]["thr"][i]["n_compact"]
                        for k in ks], color=COL[tag], ls=ls, lw=1.5,
                   label=f"{LBL[tag]} thr {thr}" if i == 1 else None)
ax[0].set_xlabel("t"); ax[0].set_ylabel("energy")
ax[0].set_title("ledger: E+KE+absorbed (solid) / absorbed (dashed)")
ax[0].legend(fontsize=7)
ax[1].set_xlabel("t"); ax[1].set_ylabel("core rotation [deg]")
ax[1].set_title("frame rotation (rotation-vs-melt read)")
ax[1].legend(fontsize=8)
ax[2].set_xlabel("t"); ax[2].set_ylabel("compact components")
ax[2].set_title("loop census (solid thr 0.09; : 0.06; -- 0.15)")
ax[2].legend(fontsize=8)
fig.suptitle("M5.21.10: n=64 free-arena dynamics, t=150 "
             "(dt 0.025, snap M_it1000 starts)")
fig.tight_layout()
out = os.path.join(PLOTS, "m5_21_10_panel.png")
fig.savefig(out, dpi=110)
print(out)
