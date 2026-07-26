"""M6.2 numeric gate (N1-N2): runs AFTER the pre-registration lock.

N1: solve the production ODE exactly as archive/sandbox_v8/ouroboros_benchmark.py
    does (same reduction, solver, tolerances, BCs, localization test) at the
    pinned parameters g = 1, lam = 1, w = 1, A0 = B0 = 0.1.
N2: evaluate on that profile every observable fixed in m6_2_preregistration.md
    section 3, and apply the section-4 decision rule. No knob exists here: the
    functionals and parameters are the pre-registered ones.

Outputs: data/m6_2_hq_decision.json, plots/m6_2_profiles.png (seeded + solved
field states per the simulation-evidence rule), plots/m6_2_term_decomposition.png.
Runtime: seconds.
"""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent.parent
OUT_JSON = HERE / "data" / "m6_2_hq_decision.json"
PLOT_PROF = HERE / "plots" / "m6_2_profiles.png"
PLOT_TERMS = HERE / "plots" / "m6_2_term_decomposition.png"

# ---- pinned configuration (pre-registration section 1; production verbatim)
G, LAM, W, A0, B0 = 1.0, 1.0, 1.0, 0.1, 0.1
R_INNER, R_MAX, N_GRID = 0.02, 12.0, 800
TAIL_THRESH = 0.15
TARGET_PHYS, TARGET_INTERNAL = 1.6875, 1.6969

def charged_odes(r, y, g, omega, lam):
    a, da, b, db = y
    d2a = -da / r + a / r**2 - omega**2 * a + b
    d2b = -db / r + b / r**2 - omega**2 * b + a - lam * b - 4 * g * b**3
    return [da, d2a, db, d2b]

# ------------------------------------------------------------- N1: solve
r_eval = np.linspace(R_INNER, R_MAX, N_GRID)
y0 = [A0 * R_INNER, A0, B0 * R_INNER, B0]
sol = solve_ivp(charged_odes, [R_INNER, R_MAX], y0, args=(G, W, LAM),
                t_eval=r_eval, method="RK45", rtol=1e-9, atol=1e-11,
                max_step=0.02)
assert sol.success, "production solve failed"
r, a, b = sol.t, sol.y[0], sol.y[2]
tail = float(abs(a[-1]) + abs(b[-1]))

# midpoint discretization identical to the production observables code
dr = np.diff(r)
rm = 0.5 * (r[:-1] + r[1:])
am = 0.5 * (a[:-1] + a[1:])
bm = 0.5 * (b[:-1] + b[1:])
dam = np.diff(a) / dr
dbm = np.diff(b) / dr

def I(f):
    """integral of f * rho drho on the midpoint grid"""
    return float(np.sum(f * rm * dr))

# ------------------------------------------------- N2: the observables
# production-coded pair (configuration certificate)
H_p = I(dam**2 + dbm**2 + am**2 / rm**2 + bm**2 / rm**2 + 4 * G * bm**4)
Q_p = I(bm**2)

# derived pair (pre-registration DF1/DF2, matched phases)
quad = (W**2 * am**2 + dam**2 + am**2 / rm**2
        + W**2 * bm**2 + dbm**2 + bm**2 / rm**2)
H_d = 0.25 * I(quad) - 0.5 * I(am * bm) + 0.375 * G * I(bm**4)
Q_d = W * I(am**2 + bm**2)
# labeled branches
Q_Jonly = W * I(bm**2)
H_d_printed_phases = 0.25 * I(quad) + 0.375 * G * I(bm**4)  # <J.A> = 0

results = {
    "N1_config": {"g": G, "lam": LAM, "omega": W, "A0": A0, "B0": B0,
                  "tail": round(tail, 5), "localized": tail < TAIL_THRESH},
    "reference_production_coded": {
        "H_p": round(H_p, 6), "Q_p": round(Q_p, 6),
        "HQ": round(H_p / Q_p, 6),
        "expected_from_v11_sec8": 1.6890,
        "reproduces_config": abs(H_p / Q_p - 1.6890) / 1.6890 < 0.02,
    },
    "PRIMARY_derived": {
        "H_d": round(H_d, 6), "Q_d": round(Q_d, 6),
        "HQ": round(H_d / Q_d, 6),
        "gap_vs_1.6875_pct": round(abs(H_d / Q_d - TARGET_PHYS) / TARGET_PHYS * 100, 3),
        "gap_vs_1.6969_pct": round(abs(H_d / Q_d - TARGET_INTERNAL) / TARGET_INTERNAL * 100, 3),
    },
    "branch_Jonly_Q": {
        "HQ": round(H_d / Q_Jonly, 6),
        "gap_vs_1.6875_pct": round(abs(H_d / Q_Jonly - TARGET_PHYS) / TARGET_PHYS * 100, 3),
        "gap_vs_1.6969_pct": round(abs(H_d / Q_Jonly - TARGET_INTERNAL) / TARGET_INTERNAL * 100, 3),
    },
    "secondary_printed_phases": {
        "HQ_with_Qd": round(H_d_printed_phases / Q_d, 6),
        "HQ_with_QJonly": round(H_d_printed_phases / Q_Jonly, 6),
    },
    "term_decomposition_H_d": {
        "quadratic_quarter": round(0.25 * I(quad), 6),
        "cross_minus_half_ab": round(-0.5 * I(am * bm), 6),
        "quartic_3_8_g": round(0.375 * G * I(bm**4), 6),
    },
}

# decision rule (pre-registration section 4)
gap = abs(H_d / Q_d - TARGET_PHYS) / TARGET_PHYS * 100
verdict = ("branch (a): survives" if gap <= 1.0
           else "partial: user decides" if gap <= 5.0
           else "branch (b): fails")
results["DECISION"] = {"primary_gap_pct_vs_1.6875": round(gap, 3),
                       "rule": "<=1% survives; 1-5% partial; >5% fails",
                       "verdict": verdict}

# ------------------------------------------------------------------ plots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))
ax1.plot(r, A0 * r * np.exp(-r), ":", color="gray", lw=1,
         label="seed slope-BC shape (illustr.)")
ax1.plot(r, a, label="alpha(rho) [A-field]", lw=1.6)
ax1.plot(r, b, label="beta(rho) [J-field]", lw=1.6)
ax1.axhline(0, color="k", lw=0.5)
ax1.set_xlabel("rho (code units)"); ax1.set_ylabel("field")
ax1.set_title(f"M6.2 N1: production solve at g=1, w=1, lam=1, A0=B0=0.1\n"
              f"tail={tail:.4f} (<{TAIL_THRESH} localized)")
ax1.legend(fontsize=8); ax1.grid(alpha=0.3)

labels = ["quad/4", "-ab/2", "3g b^4/8"]
vals = [results["term_decomposition_H_d"][k] for k in
        ("quadratic_quarter", "cross_minus_half_ab", "quartic_3_8_g")]
ax2.bar(labels, vals, color=["tab:blue", "tab:red", "tab:green"])
ax2.axhline(0, color="k", lw=0.5)
ax2.set_ylabel("contribution to H_d")
ax2.set_title("H_d term decomposition (derived functional, matched phases)")
ax2.grid(alpha=0.3, axis="y")
fig.tight_layout()
PLOT_PROF.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(PLOT_PROF, dpi=150)

fig2, ax = plt.subplots(figsize=(7, 4.2))
rows = [("H_p/Q_p (production-coded)", H_p / Q_p),
        ("H_d/Q_d (PRIMARY derived)", H_d / Q_d),
        ("H_d/Q_Jonly (branch)", H_d / Q_Jonly),
        ("printed-phase H/Q_d", H_d_printed_phases / Q_d)]
ax.barh([x[0] for x in rows], [x[1] for x in rows], color="tab:blue")
ax.axvline(TARGET_PHYS, color="tab:green", ls="--", lw=1.2,
           label="physical target 1.6875")
ax.axvline(TARGET_INTERNAL, color="tab:orange", ls=":", lw=1.2,
           label="model-internal 1.6969")
ax.set_xlabel("H/Q"); ax.legend(fontsize=8)
ax.set_title("M6.2 decision gate: pre-registered observables vs targets")
fig2.tight_layout()
fig2.savefig(PLOT_TERMS, dpi=150)

OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
OUT_JSON.write_text(json.dumps(results, indent=2))
print(json.dumps(results, indent=2))
print(f"\nwritten: {OUT_JSON}\nplots: {PLOT_PROF}, {PLOT_TERMS}")
