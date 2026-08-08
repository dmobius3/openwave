"""M5.21.11 close-out panel: the measured record, gates visible.

Three panels over the 42-run ladder record (data-only, no fit claims:
the frozen fit was unreachable at 1 usable point):
  (1) E(delta) per branch at N = 48, filled markers = § 4-usable,
      open = gate-excluded (the F3 story visible at a glance)
  (2) the § 3 refinement subset E(N) at delta in {0.30, 0.12, 0.05}
  (3) the g-arm dressing gain vs artanh(1/g) with the slope-2 reference
      (the F4 story: measured flat, required falling)

Out: ../plots/m5_21_11_panel.png
"""
from __future__ import annotations

import glob
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
PLOTS = os.path.join(HERE, "..", "plots")

COLS = {"A": "tab:blue", "C": "tab:orange", "B": "tab:red"}
REFINED = (0.30, 0.12, 0.05)


def main():
    rows = {}
    for p in glob.glob(os.path.join(DATA,
                                    "m5_21_11_row_t11lad_*.json")):
        r = json.load(open(p))
        rows[(r["branch"], int(r["n"]), round(r["delta"], 3))] = r
    fit = json.load(open(os.path.join(DATA, "m5_21_11_fit.json")))
    garm = json.load(open(os.path.join(DATA, "m5_21_11_garm.json")))

    fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.4))

    for br in ("A", "C", "B"):
        ds_u, es_u, ds_x, es_x = [], [], [], []
        for (b, n, d), r in rows.items():
            if b != br or n != 48:
                continue
            st = fit["rung_state"].get(f"{br}_n48_d{d:g}", {})
            if st.get("state") == "USABLE":
                ds_u.append(d)
                es_u.append(r["E_end"])
            else:
                ds_x.append(d)
                es_x.append(r["E_end"])
        o = np.argsort(ds_u + ds_x)
        alld = np.array(ds_u + ds_x)[o]
        alle = np.array(es_u + es_x)[o]
        ax[0].plot(alld, alle, "-", color=COLS[br], lw=0.8, alpha=0.5)
        ax[0].plot(ds_u, es_u, "o", color=COLS[br],
                   label=f"{br} usable")
        ax[0].plot(ds_x, es_x, "o", mfc="none", color=COLS[br],
                   label=f"{br} excluded (§ 4)")
    ax[0].set_xlabel("delta")
    ax[0].set_ylabel("E_end (N = 48)")
    ax[0].set_yscale("log")
    ax[0].invert_xaxis()
    ax[0].set_title("production ladder: E(delta), gates visible")
    ax[0].legend(fontsize=7)

    for br in ("A", "C", "B"):
        for d in REFINED:
            es = [(n, rows[(br, n, round(d, 3))]["E_end"])
                  for n in (32, 48, 64)
                  if (br, n, round(d, 3)) in rows]
            if len(es) == 3:
                ax[1].plot([e[0] for e in es], [e[1] for e in es],
                           "o-", color=COLS[br], alpha=0.75)
                ax[1].annotate(f"{br} d={d:g}", (64, es[-1][1]),
                               fontsize=6, color=COLS[br])
    ax[1].set_xlabel("N (L = 48 fixed)")
    ax[1].set_ylabel("E_end")
    ax[1].set_yscale("log")
    ax[1].set_title("refinement subset (§ 3)")

    for br in ("A", "C", "B"):
        f = garm["fits"][br]
        gs = [8.0, 16.0, 32.0]
        x = [np.arctanh(1.0 / g) for g in gs]
        y = [-f["gains"][f"g{g:g}"] for g in gs]
        ax[2].loglog(x, y, "o-", color=COLS[br],
                     label=f"{br}  q = {f['q_lsq']:+.2f}")
    xr = np.array([np.arctanh(1 / 32.0), np.arctanh(1 / 8.0)])
    ax[2].loglog(xr, 0.65 * (xr / xr[-1]) ** 2, "k--", lw=0.8,
                 label="slope 2 (F4 bar)")
    ax[2].set_xlabel("artanh(1/g)")
    ax[2].set_ylabel("-gain (dressing)")
    ax[2].set_title("g-arm: measured FLAT, F4 requires slope >= 2")
    ax[2].legend(fontsize=7)

    fig.tight_layout()
    os.makedirs(PLOTS, exist_ok=True)
    out = os.path.join(PLOTS, "m5_21_11_panel.png")
    fig.savefig(out, dpi=140)
    print("wrote", out)


if __name__ == "__main__":
    main()
