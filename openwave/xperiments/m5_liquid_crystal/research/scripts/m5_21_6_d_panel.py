"""M5.21.6 panel: the decay-run summary figure (4 quadrants).

(a) free-arena E(t) traces (A holds, C -> A level, B under-drains)
(b) Arena-1 kick ladder: E_kicked vs E_end per kick (return lines)
(c) endpoint core spectra vs the vacuum spectrum
(d) compact biaxial-component count vs iteration (the two-loop read)

Run:  python3 m5_21_6_d_panel.py
Out:  ../plots/m5_21_6_panel.png
"""
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
PLOTS = os.path.join(HERE, "..", "plots")

fig, ax = plt.subplots(2, 2, figsize=(13, 10))

# (a) E(t)
tr = json.load(open(os.path.join(DATA, "m5_21_6_f48_traces.json")))
for t, c in (("A", "tab:green"), ("B", "tab:red"), ("C", "tab:blue")):
    its = [r[0] for r in tr[t]]
    es = [r[1] for r in tr[t]]
    ax[0, 0].plot(its, es, c=c, label=f"seed {t}")
ax[0, 0].set_yscale("log")
ax[0, 0].set_xlabel("FIRE iteration")
ax[0, 0].set_ylabel("E")
ax[0, 0].set_title("(a) free arena n=48: C joins the A level, "
                   "B drains under it")
ax[0, 0].legend()
ax[0, 0].grid(alpha=0.3)

# (b) kick ladder
for i, (tg, c) in enumerate((("p32_B", "tab:red"),
                             ("p32_C", "tab:blue"))):
    d = json.load(open(os.path.join(DATA, f"m5_21_6_p1_{tg}.json")))
    E0 = d["E_start"]
    xs = np.arange(len(d["rows"])) + i * 0.35
    ek = [r["E_kicked"] for r in d["rows"]]
    ee = [r["E_end"] for r in d["rows"]]
    ax[0, 1].bar(xs, ek, width=0.3, color=c, alpha=0.35,
                 label=f"{tg} kicked")
    ax[0, 1].bar(xs, ee, width=0.3, color=c, alpha=1.0,
                 label=f"{tg} end")
    ax[0, 1].axhline(E0, color=c, ls=":", lw=1)
labels = [f"{r['family']}:{r['amp']:g}" for r in d["rows"]]
ax[0, 1].set_xticks(np.arange(len(labels)) + 0.17, labels,
                    rotation=45)
ax[0, 1].set_yscale("log")
ax[0, 1].set_title("(b) pinned kicks: every kick relaxes back "
                   "(dotted = start level)")
ax[0, 1].legend(fontsize=7)
ax[0, 1].grid(alpha=0.3)

# (c) core spectra
specs = json.load(open(os.path.join(DATA, "m5_21_6_spec_core.json")))
w = 0.2
for i, lam in enumerate(("lam1", "lam2", "lam3")):
    vals = [specs[t][i] for t in ("A", "C", "B")]
    ax[1, 0].bar(np.arange(3) + i * w, vals, width=w,
                 label=f"eig {i + 1}")
for i, v in enumerate((1.0, 0.3, 0.0)):
    ax[1, 0].axhline(v, color="gray", ls=":", lw=1)
ax[1, 0].set_xticks(np.arange(3) + w, ["A (ref)", "C", "B"])
ax[1, 0].set_title("(c) endpoint core spectra (r<6): C == A; "
                   "B near-vacuum (dotted = vacuum 1, 0.3, 0)")
ax[1, 0].legend()
ax[1, 0].grid(alpha=0.3)

# (d) compact component counts
for t, c in (("A", "tab:green"), ("B", "tab:red"), ("C", "tab:blue")):
    p = os.path.join(DATA, f"m5_21_6_lp_end_f48_{t}_loops.json")
    d = json.load(open(p))
    pts = []
    for k, reads in d["snaps"].items():
        it = 12000 if k == "M" else int(k[len("M_it"):])
        r9 = [r for r in reads if abs(r["thr"] - 0.09) < 1e-9][0]
        pts.append((it, r9["n_compact"]))
    pts.sort()
    ax[1, 1].plot([p[0] for p in pts], [p[1] for p in pts],
                  "o-", c=c, label=f"seed {t}")
ax[1, 1].axhline(2, color="k", ls="--", lw=1,
                 label="2 (the electron pair / the author's conjecture)")
ax[1, 1].set_xlabel("FIRE iteration")
ax[1, 1].set_ylabel("compact biaxial components (thr 0.09)")
ax[1, 1].set_title("(d) topology bookkeeping: A holds 2, C -> 2, "
                   "B -> 0")
ax[1, 1].legend(fontsize=8)
ax[1, 1].grid(alpha=0.3)

fig.suptitle("M5.21.6: the 3-lepton decay runs (T2 sym instrument; "
             "descent = basin topography, not dynamics)")
fig.tight_layout()
out = os.path.join(PLOTS, "m5_21_6_panel.png")
fig.savefig(out, dpi=110)
print("wrote", out)
