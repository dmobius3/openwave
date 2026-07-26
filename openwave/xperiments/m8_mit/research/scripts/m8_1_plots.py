#!/usr/bin/env python3
"""M8.1 designer-side comparison plots.

Reads the blind solver's JSONs (data/m8_1_spectrum.json, data/m8_1_delta_scan.json)
and overlays the author's claimed formulas (Theorems 1.1/1.2 of the first-eigenvalue
paper), which the solver never saw:
  lambda1+(W) = 2/R^2                     for 0 < W <= pi R/2
  lambda1+(W) = a0(a0+1)/R^2, a0 = pi R/2W, for W > pi R/2
  stability threshold delta0 = 2R/e; bound state lambda_b ~ -4 e^(-2 gamma)/delta0^2
Run from this scripts/ directory or the repo root; outputs land in ../plots/.
"""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
PLOTS = HERE.parent / "plots"
PLOTS.mkdir(exist_ok=True)

GAMMA = 0.57721566490153286
R = 1.0

spec = json.load(open(DATA / "m8_1_spectrum.json"))
scan = json.load(open(DATA / "m8_1_delta_scan.json"))

# ---- Plot 1: lambda1+ vs W, solver points vs claimed curves -------------------
wkeys = sorted(spec["per_W"], key=float)
Ws = [float(k) for k in wkeys]
lam1 = [spec["per_W"][k]["lambda1plus"]["lambda"] for k in wkeys]

Wf = np.linspace(0.05, 3.1, 400)
claim = np.where(Wf <= np.pi / 2, 2.0, (np.pi / (2 * Wf)) * (np.pi / (2 * Wf) + 1))

fig, (ax, axr) = plt.subplots(
    2, 1, figsize=(7.2, 5.6), height_ratios=[3, 1], sharex=True
)
ax.plot(Wf, claim, "-", color="tab:orange", lw=1.5, label="claimed: 2 (narrow) / α₀(α₀+1) (wide)")
ax.plot(Ws, lam1, "o", color="tab:blue", ms=6, label="blind solver λ₁⁺ (Friedrichs)")
ax.axvline(np.pi / 2, color="gray", ls=":", lw=1, label="W = πR/2 (claimed crossing)")
ax.set_ylabel("λ₁⁺ · R²")
ax.legend(fontsize=8)
ax.set_title("M8.1 gate: first positive eigenvalue vs claimed formula (solver was blind)")

res = [
    l - (2.0 if w <= np.pi / 2 else (np.pi / (2 * w)) * (np.pi / (2 * w) + 1))
    for w, l in zip(Ws, lam1)
]
axr.axhline(0, color="gray", lw=0.8)
axr.plot(Ws, res, "s", color="tab:red", ms=4)
axr.set_ylabel("residual")
axr.set_xlabel("W / R")
axr.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
fig.tight_layout()
fig.savefig(PLOTS / "m8_1_lambda1_vs_W.png", dpi=150)
plt.close(fig)

# ---- Plot 2: the delta0 scan --------------------------------------------------
d0 = np.array([r["delta0"] for r in scan["per_delta0"]])
neg = np.array([r["negative_eigenvalues"][0] for r in scan["per_delta0"]])
pos = np.array([r["first_positive"] for r in scan["per_delta0"]])

fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.6, 4.0))
a1.loglog(d0, -neg, "o", color="tab:blue", label="solver: −λ_b(δ₀)")
d0f = np.logspace(-3.1, 1.05, 200)
a1.loglog(d0f, 4 * np.exp(-2 * GAMMA) / d0f**2, "-", color="tab:orange", lw=1.2,
          label="claimed: 4e⁻²ᵞ/δ₀²")
a1.set_xlabel("δ₀ / R")
a1.set_ylabel("−λ_b · R²")
a1.legend(fontsize=8)
a1.set_title("bridging bound state vs claimed asymptotic")

a2.semilogx(d0, pos, "o", color="tab:blue", label="solver: first positive")
a2.axhline(2.0, color="tab:orange", lw=1.2, label="claimed stable value 2")
a2.axvline(2 / np.e, color="gray", ls=":", lw=1, label="claimed threshold 2R/e")
a2.set_xlabel("δ₀ / R")
a2.set_ylabel("λ₁⁺ · R²")
a2.legend(fontsize=8)
a2.set_title("extension-stability of the first positive level (W = 1)")
fig.tight_layout()
fig.savefig(PLOTS / "m8_1_delta_scan.png", dpi=150)
plt.close(fig)

# ---- Plot 3: convergence ladder ----------------------------------------------
fig, ax = plt.subplots(figsize=(6.4, 4.2))
bykey = {k: v for k, v in spec["per_W"].items()}
def wrec(x):
    for k, v in bykey.items():
        if abs(float(k) - x) < 1e-9:
            return k, v
    raise KeyError(x)
for wval, color in [(0.9, "tab:blue"), (1.5, "tab:green"), (2.2, "tab:red"), (3.0, "tab:purple")]:
    wkey, node = wrec(wval)
    rec = node["lambda1plus"]
    raw = np.array(rec["raw_by_resolution"])
    ns = np.array(node["resolutions_elements_total"], dtype=float)[: len(raw)]
    ax.loglog(ns, np.abs(raw - rec["lambda"]), "o-", color=color,
              label=f"W = {wkey} (q_fit = {rec['q_fit']:.2f})")
nn = np.array([2.5e2, 2.5e3])
ax.loglog(nn, 3e2 * nn**-2.0, "--", color="gray", lw=1, label="slope −2 guide")
ax.set_xlabel("total elements")
ax.set_ylabel("|λ₁⁺(N) − λ₁⁺(extrap)|")
ax.legend(fontsize=8)
ax.set_title("λ₁⁺ convergence ladder (solver resolutions)")
fig.tight_layout()
fig.savefig(PLOTS / "m8_1_convergence.png", dpi=150)
plt.close(fig)

print("wrote", *[p.name for p in sorted(PLOTS.glob("m8_1_*.png"))])
