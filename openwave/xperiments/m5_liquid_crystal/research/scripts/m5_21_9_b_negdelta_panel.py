"""M5.21.9 pre-gate panel: the negative-delta suggestion, verified.

(a) sign map of the author's Hm omega^2 coefficient over (delta, g),
    now including delta < 0 (the author's 2026-07-19 suggestion):
    the whole delta < 0 half-plane is BOUNDED (positive coefficient)
(b) E(omega) at g = 8: delta = +0.3 runaway vs delta = -0.3 bounded
    (minimum at omega -> 0: bounded but STATIC, no free clock)
(c) the cone ladder E(Delta): negatively divergent at +0.3 (collapse
    channel) vs POSITIVELY divergent at -0.3 (normal vortex-core mass)

Run:  python3 m5_21_9_b_negdelta_panel.py
Out:  ../plots/m5_21_9_negdelta.png
"""
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
PLOTS = os.path.join(HERE, "..", "plots")
sys.path.insert(0, HERE)

from m5_21_8_a_verify import E_cone, mstar, om2_coeff  # noqa: E402

nd = json.load(open(os.path.join(DATA, "m5_21_9_negdelta.json")))
p0 = json.load(open(os.path.join(DATA, "m5_21_8_p0.json")))

fig, ax = plt.subplots(1, 3, figsize=(16.5, 5.0))

# (a) sign map over (delta symlog, g log)
dls = np.concatenate([-np.logspace(np.log10(3.5), -12, 140),
                      np.logspace(-12, np.log10(3.5), 140)])
gs = np.logspace(np.log10(1.05), 10, 160)
S = np.array([[np.sign(om2_coeff(g, d)) for d in dls] for g in gs])
ax[0].pcolormesh(dls, gs, S, cmap="RdBu", vmin=-1, vmax=1,
                 shading="nearest")
gb = np.logspace(np.log10(1.6), 10, 80)
ax[0].plot(2.0 * gb / (gb - 1.0), gb, "k-", lw=1.2,
           label="delta = 2g/(g-1)")
ax[0].axvline(0.0, color="k", lw=1.2, ls="--", label="delta = 0")
for d, g, mk in ((-0.3, 8.0, "s"), (0.3, 8.0, "o"),
                 (-1e-10, 1e10, "^"), (1e-10, 1e10, "v")):
    ax[0].plot([d], [g], mk, color="k", ms=7, mfc="yellow")
ax[0].set_xscale("symlog", linthresh=1e-10)
ax[0].set_yscale("log")
ax[0].set_xticks([-1.0, -1e-5, -1e-10, 1e-10, 1e-5, 1.0])
ax[0].set_xticklabels(["-1", "-1e-5", "-1e-10", "1e-10", "1e-5", "1"])
ax[0].set_xlabel("delta (symlog)")
ax[0].set_ylabel("g")
ax[0].set_title("(a) sign of the author's Hm omega^2 coefficient\n"
                "= 8 delta^3 [(g-1) delta - 2g]/(g-1): the WHOLE\n"
                "delta < 0 half is bounded (blue), incl. (1e10, -1e-10)")
ax[0].legend(fontsize=8, loc="lower left")

# (b) E(omega) at g = 8, delta = +/-0.3
row = nd["N4_linearity_omstar0"]["pipeline"]
oms = np.array(row["oms"])
Es_n = np.array(row["E"])
m = -abs(mstar(8.0))
Es_p = np.array([E_cone(8.0, 0.3, m, om, r=1.3) for om in oms])
ax[1].plot(oms[1:], Es_n[1:], "s-", color="tab:blue",
           label=f"delta = -0.3 (slope +{row['slope']:.3f}: bounded)")
ax[1].plot(oms[0], Es_n[0], "s", mfc="none", color="tab:blue",
           label="omega = 0 point (t = 0 snapshot evaluator)")
ax[1].plot(oms[1:], Es_p[1:], "o-", color="tab:red",
           label="delta = +0.3 (negative slope: the -inf branch)")
ax[1].plot(oms[0], Es_p[0], "o", mfc="none", color="tab:red")
ax[1].set_xlabel("omega")
ax[1].set_ylabel("E (pipeline, r = 1.3, m = -|m*|)")
ax[1].set_title("(b) E(omega) at g = 8: negative delta flips the\n"
                "omega channel to bounded, minimum at omega -> 0\n"
                "(STATIC: still no free de Broglie clock)")
ax[1].legend(fontsize=8)
ax[1].grid(alpha=0.3)

# (c) cone ladders at +/-0.3
c4 = next(r for r in p0["C4_cone_divergence"]
          if r.get("dl") == 0.3)
Dsp = sorted(float(k) for k in c4["E_at_Delta"])
Ep = [c4["E_at_Delta"][f"{d:g}"] for d in Dsp]
n8 = nd["N8_cone_channel"]
Dsn = sorted(float(k) for k in n8["E_at_Delta"])
En = [n8["E_at_Delta"][f"{d:g}"] for d in Dsn]
ax[2].semilogx(Dsp, Ep, "o-", color="tab:red",
               label="delta = +0.3: E FALLS as the cone shrinks\n"
                     "(negatively divergent core channel)")
ax[2].semilogx(Dsn, En, "s-", color="tab:blue",
               label="delta = -0.3: E RISES as the cone shrinks\n"
                     "(positive vortex-core mass, log slope +0.28)")
ax[2].invert_xaxis()
ax[2].set_xlabel("cone cut Delta (shrinking to the core ->)")
ax[2].set_ylabel("E(Delta)")
ax[2].set_title("(c) the dropped cone channel: negative delta\n"
                "flips the core divergence from collapse to\n"
                "normal positive core energy")
ax[2].legend(fontsize=8)
ax[2].grid(alpha=0.3)

fig.suptitle("M5.21.9 pre-gate: the author's negative-delta "
             "suggestion verified: BOTH runaway channels flip to "
             "bounded; the free minimum stays static (omega* = 0)")
fig.tight_layout()
out = os.path.join(PLOTS, "m5_21_9_negdelta.png")
fig.savefig(out, dpi=110)
print("wrote", out)
