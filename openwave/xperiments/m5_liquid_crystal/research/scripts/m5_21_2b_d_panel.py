"""M5.21.2b summary panel: the six headline reads in one figure.

Run:  python3 m5_21_2b_d_panel.py
Out:  ../plots/m5_21_2b_panel.png
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
J = json.load(open(os.path.join(DATA, "m5_21_2b_all.json")))

fig, axs = plt.subplots(2, 3, figsize=(16, 9))

# (1) the instrument: xstencil ratio per I1 arm (log)
ax = axs[0, 0]
arms = ["i1_A_fwd_e0", "i1_R_fwd_e0", "i1_A_sym_e0", "i1_R_sym_e0",
        "i1_A_sym_e0.001", "i1_R_sym_e0.001"]
vals = [J[a]["consistency"]["xstencil_ratio"] for a in arms]
cols = ["crimson"] * 2 + ["seagreen"] * 4
ax.bar(range(len(arms)), vals, color=cols)
ax.set_yscale("log")
ax.axhline(1.5, color="k", ls="--", lw=1, label="pass bar 1.5")
ax.set_xticks(range(len(arms)),
              [a[3:].replace("_", "\n") for a in arms], fontsize=7)
ax.set_title("I1: cross-stencil ratio (log): fwd funnels x294-386,\n"
             "sym is consistent (1.10-1.11), both defect types")
ax.legend(fontsize=8)

# (2) the eps-family E(eps) linearity -> eps = 0 intercept
ax = axs[0, 1]
eps = [0.0, 0.0003, 0.001, 0.003]
for s, mk in (("A", "o"), ("R", "s")):
    E = [J[f"i1_{s}_sym_e{e:g}" if e else f"i1_{s}_sym_e0"]["E_end"]
         for e in eps]
    ax.plot(eps, E, mk + "-", label=f"seed {s}")
ax.set_xlabel("eps (Dirichlet)")
ax.set_ylabel("E_end")
ax.set_title("I1: the eps-family is linear and lands on the\n"
             "bare eps = 0 endpoint (well-posedness certificate)")
ax.legend(fontsize=8)

# (3) census E (log) at N = 48
ax = axs[0, 2]
runs = ["c48_A_d0.3", "c48_R_d0.3", "c48_C_d0.3", "c48_B_d0.3"]
E = [J[r]["E_end"] for r in runs]
xr = [J[r]["consistency"]["xstencil_ratio"] for r in runs]
cols = ["seagreen" if x < 1.5 else "goldenrod" for x in xr]
ax.bar(range(4), E, color=cols)
ax.set_yscale("log")
for i, (e, x) in enumerate(zip(E, xr)):
    ax.text(i, e * 1.05, f"xr {x:.2f}", ha="center", fontsize=8)
ax.set_xticks(range(4), ["A (1,d,0)", "ring R", "C (0,1,d)",
                         "B (d,0,1)"], fontsize=8)
ax.set_title("I3 census N = 48: A < C < B stationary minima;\n"
             "A = R (basin merger); gold = lattice-flagged energy")

# (4) virial residual: term sets + the T1 eps crossing
ax = axs[1, 0]
labs = ["T1 e0", "T1 e3e-4", "T1 e1e-3", "T1 e3e-3", "T2 e0", "T3 e0"]
keys = ["i1_A_sym_e0", "i1_A_sym_e0.0003", "i1_A_sym_e0.001",
        "i1_A_sym_e0.003", "i2_A_T2", "i2_A_T3"]
v = [J[k]["virial_resid"] for k in keys]
cols = ["steelblue"] * 4 + ["seagreen", "indianred"]
ax.bar(range(len(v)), v, color=cols)
ax.axhline(0, color="k", lw=1)
ax.set_xticks(range(len(v)), labs, fontsize=8)
ax.set_title("I2: Derrick residual (0 = scale-stationary):\n"
             "T1 crosses on the eps ladder; T2 lands at eps = 0")

# (5) the refinement ladder E(h)
ax = axs[1, 1]
hs = [1.5, 1.0, 0.75]
EA = [J["i1_A_sym_e0"]["E_end"], J["c48_A_d0.3"]["E_end"],
      J["lad64_A"]["E_end"]]
ER = [J["i1_R_sym_e0"]["E_end"], J["c48_R_d0.3"]["E_end"],
      J["lad64_R"]["E_end"]]
ax.plot(hs, EA, "o-", label="seed A")
ax.plot(hs, ER, "s--", label="ring R")
for h, e, n in zip(hs, EA, (32, 48, 64)):
    ax.annotate(f"N={n}", (h, e), textcoords="offset points",
                xytext=(5, 5), fontsize=8)
ax.invert_xaxis()
ax.set_xlabel("h (fixed physical box L = 48)")
ax.set_ylabel("E_end")
ax.set_title("Ladder: consistency holds every rung (xr 1.1-1.2);\n"
             "absolute E still drifts: quote orderings, not values")
ax.legend(fontsize=8)

# (6) the selected transverse scale rho_w(delta), all seedings
ax = axs[1, 2]
ds = [0.3, 0.2, 0.1]
series = {"point": ["c48_A_d0.3", "c48_A_d0.2", "c48_A_d0.1"],
          "a=3": ["i4_a3_d0.3", "i4_a3_d0.2", "i4_a3_d0.1"],
          "a=4.5": ["i4_a4.5_d0.3", "i4_a4.5_d0.2", "i4_a4.5_d0.1"],
          "a=6": ["i4_a6_d0.3", "i4_a6_d0.2", "i4_a6_d0.1"]}
for lab, ks in series.items():
    ax.plot(ds, [J[k]["ring"]["rho_w"] for k in ks], "o-", label=lab)
dd = np.linspace(0.1, 0.3, 50)
ax.plot(dd, 2.49 * (dd / 0.3) ** -0.25, "k:", lw=1,
        label="~delta^-0.25 guide")
ax.invert_xaxis()
ax.set_xlabel("delta")
ax.set_ylabel("rho_w (excess-weighted)")
ax.set_title("I4 (Q30/Q31): every seeded radius anneals to ONE\n"
             "delta-selected scale; grows as delta falls")
ax.legend(fontsize=8)

fig.suptitle("M5.21.2b: the well-posed 3D instrument, hedgehog + "
             "charged ring (sym stencil; N = 32/48/64, L = 48, "
             "pinned)", fontsize=13)
fig.tight_layout()
fig.savefig(os.path.join(PLOTS, "m5_21_2b_panel.png"), dpi=115)
print("panel written")
