"""M5.32 R29-2: the dispersion map figure from m5_32_r29_2_dispersion.json (no computation; headless matplotlib).
Left: omega^2 against k_lat^2 per M-field term on the toy vacuum (k along x1), colored by norm sign, with the
helicity-2-content modes marked. Right: the torsion scan, negative-norm count and unstable-mode count per point.
Run from this directory: python3 m5_32_r29_2_plot.py -> ../plots/m5_32_r29_2_map.png
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(HERE, "..", "data", "m5_32_r29_2_dispersion.json")))
OUT = os.path.join(HERE, "..", "plots", "m5_32_r29_2_map.png")
os.makedirs(os.path.dirname(OUT), exist_ok=True)

fig, axes = plt.subplots(1, 2, figsize=(14, 5.6))
ax = axes[0]
vn = "toy (32, 1, 0.3, 0), w W1"
terms = ["K_T (coef +1)", "K_cov (coef +1)", "K_h (coef +1)", "K_spec (coef +1)"]
colors = {"K_T (coef +1)": "tab:blue", "K_cov (coef +1)": "tab:red", "K_h (coef +1)": "tab:green", "K_spec (coef +1)": "tab:purple"}
for lab in terms:
    rows = d["M_field"][vn][lab]["per_k"]["x1"]
    for r in rows:
        for m in r["modes"]:
            if abs(m["omega2"]) > 50:          # the diag0 mode at 2.5e7 is off the chart
                continue
            mk = "o" if m["norm_sign"] > 0 else ("x" if m["norm_sign"] < 0 else ".")
            ax.plot(r["k_lat2"], m["omega2"], mk, color=colors[lab], ms=5 if m["norm_sign"] else 3, alpha=0.8)
            if m["h"]["|h|=2"] > 0.9:
                ax.plot(r["k_lat2"], m["omega2"], "s", mfc="none", mec="k", ms=10)
kk = np.linspace(0, 1.7, 50)
ax.plot(kk, kk, "k--", lw=0.8, label="omega^2 = k_lat^2 (c = 1)")
for lab in terms:
    ax.plot([], [], "o", color=colors[lab], label=lab)
ax.plot([], [], "x", color="gray", label="negative norm")
ax.plot([], [], "s", mfc="none", mec="k", label="helicity-2 content > 0.9")
ax.set_xlabel("k_lat^2 (sym stencil, h = 1.5), k along x1")
ax.set_ylabel("omega^2")
ax.set_ylim(-0.05, 1.9)
ax.set_title("R29-2: M-field terms about the toy vacuum (diag0 mode off-chart)")
ax.legend(fontsize=7, loc="upper left")

ax = axes[1]
s = d["summary"]["torsion"]
names = list(s.keys())
keys = ["mu_K 0, rank3_pseudoinverse", "mu_K 0, null_leg_dropped", "mu_K 30, rank3_pseudoinverse", "mu_K 30, null_leg_dropped"]
x = np.arange(len(names))
w = 0.2
for i, k in enumerate(keys):
    neg = [s[n][k]["max_negative_norm"] for n in names]
    unst = [s[n][k]["unstable"] for n in names]
    ax.bar(x + (i - 1.5) * w, neg, w, label=f"neg norm, {k.replace('_', ' ')}", alpha=0.85)
    ax.plot(x + (i - 1.5) * w, [min(u, 12) for u in unst], "k_", ms=6)
ax.set_xticks(x)
ax.set_xticklabels(names, rotation=70, fontsize=6.5)
ax.set_ylabel("max negative-norm modes over k (bars); unstable modes, capped at 12 (ticks)")
ax.set_title("R29-2: torsion scan (a1, a2, a3), locked variant, V4 on")
ax.legend(fontsize=6.5)
plt.tight_layout()
plt.savefig(OUT, dpi=130)
print("wrote", os.path.relpath(OUT, HERE))
