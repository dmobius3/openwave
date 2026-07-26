"""M6.2 post-gate diagnostic: r_max sensitivity of H/Q (window-dependence).

Runs AFTER the locked gate (m6_2_hq_decision.py). Not part of the decision:
the gate number is fixed at the production r_max = 12. This diagnostic tests
the STABILITY of the production claim itself: the N1 profile shows a zero
crossing near rho ~ 6 and a rising tail at the domain edge (a slowly damped
radial oscillation, not a nodeless bound state), which predicts that H and Q
are integration-window quantities. If H/Q drifts with r_max, the benchmark
number is window-defined, the same property the author's July 8 paper
concedes for the M7 3D measurement.

Output: data/m6_2_rmax_sensitivity.json + plots/m6_2_rmax_sensitivity.png.
"""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent.parent
OUT = HERE / "data" / "m6_2_rmax_sensitivity.json"
PLOT = HERE / "plots" / "m6_2_rmax_sensitivity.png"

G, LAM, W, A0, B0 = 1.0, 1.0, 1.0, 0.1, 0.1
R_INNER = 0.02

def charged_odes(r, y, g, omega, lam):
    a, da, b, db = y
    d2a = -da / r + a / r**2 - omega**2 * a + b
    d2b = -db / r + b / r**2 - omega**2 * b + a - lam * b - 4 * g * b**3
    return [da, d2a, db, d2b]

def run(r_max, n_grid):
    r_eval = np.linspace(R_INNER, r_max, n_grid)
    y0 = [A0 * R_INNER, A0, B0 * R_INNER, B0]
    sol = solve_ivp(charged_odes, [R_INNER, r_max], y0, args=(G, W, LAM),
                    t_eval=r_eval, method="RK45", rtol=1e-9, atol=1e-11,
                    max_step=0.02)
    if not sol.success:
        return None
    r, a, b = sol.t, sol.y[0], sol.y[2]
    dr = np.diff(r)
    rm, am, bm = (0.5 * (v[:-1] + v[1:]) for v in (r, a, b))
    dam, dbm = np.diff(a) / dr, np.diff(b) / dr
    I = lambda f: float(np.sum(f * rm * dr))
    H_p = I(dam**2 + dbm**2 + am**2 / rm**2 + bm**2 / rm**2 + 4 * G * bm**4)
    Q_p = I(bm**2)
    quad = (W**2 * am**2 + dam**2 + am**2 / rm**2
            + W**2 * bm**2 + dbm**2 + bm**2 / rm**2)
    H_d = 0.25 * I(quad) - 0.5 * I(am * bm) + 0.375 * G * I(bm**4)
    Q_d = W * I(am**2 + bm**2)
    tail = float(abs(a[-1]) + abs(b[-1]))
    zeros_b = int(np.sum(np.sign(bm[:-1]) != np.sign(bm[1:])))
    return {"r_max": r_max, "tail": round(tail, 5),
            "beta_sign_changes": zeros_b,
            "H_p": round(H_p, 6), "Q_p": round(Q_p, 6),
            "HQ_production": round(H_p / Q_p, 6),
            "HQ_derived": round(H_d / Q_d, 6)}

rows = [run(rm_, int(800 * rm_ / 12)) for rm_ in (12, 18, 24, 36, 48)]
rows = [r_ for r_ in rows if r_]
drift = (max(r_["HQ_production"] for r_ in rows)
         - min(r_["HQ_production"] for r_ in rows))
result = {
    "rows": rows,
    "HQ_production_drift_absolute": round(drift, 4),
    "HQ_production_drift_pct_of_1p689": round(drift / 1.689 * 100, 2),
    "reading": ("if the drift is comparable to or larger than the claimed "
                "0.09%-0.56% accuracy, the benchmark H/Q is a window-defined "
                "quantity at the production r_max, matching the July 8 "
                "concession for the 3D case"),
}

fig, ax = plt.subplots(figsize=(7, 4.2))
ax.plot([r_["r_max"] for r_ in rows], [r_["HQ_production"] for r_ in rows],
        "o-", label="H_p/Q_p (production-coded)")
ax.axhline(1.6875, color="tab:green", ls="--", lw=1, label="1.6875 target")
ax.axhline(1.6969, color="tab:orange", ls=":", lw=1, label="1.6969 internal")
ax.set_xlabel("r_max (integration window)")
ax.set_ylabel("H/Q")
ax.set_title("M6.2 diagnostic: window dependence of the benchmark H/Q")
ax.legend(fontsize=8); ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(PLOT, dpi=150)

OUT.write_text(json.dumps(result, indent=2))
print(json.dumps(result, indent=2))
print(f"\nwritten: {OUT}\nplot: {PLOT}")
