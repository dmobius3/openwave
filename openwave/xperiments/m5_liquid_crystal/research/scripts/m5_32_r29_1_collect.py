"""M5.32 R29-1 collect: the ladder table, the 1/L fits and the figure from the rung files of m5_32_r29_1_gap.py.

EQUATIONS FIRST
---------------
Per rung (box L at h 1.5; the h 1.0 rung is the mesh check, not a ladder point): the lowest K-positive omega^2, bare
(the converged LOBPCG pass) and with the K-null auxiliary elimination (the CG post-correction), and the K-negative
witness (the Lanczos read). The pre-registered read: the lowest K-positive omega^2 fitted linearly in 1/L over the
rungs, intercept c0 and slope c1, residual sigma_fit = sqrt(sum r^2 / (n - 2)); the error of the intercept carries
sigma_fit, the solver drift of each rung and the mesh difference (h 1.0 against h 1.5 on the same box) in quadrature;
GAPPED if c0 > 2 sigma, GAPLESS if |c0| < 2 sigma and the lowest mode is an exterior texture, otherwise reported as
not resolved. A 1/L^2 fit (the scaling of a box mode) is given beside it as the diagnostic.
Run from this directory: python3 m5_32_r29_1_collect.py -> ../data/m5_32_r29_1_gap.json, ../plots/m5_32_r29_1_ladder.png
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "..", "data", "m5_32_r29_1")
OUT = os.path.join(HERE, "..", "data", "m5_32_r29_1_gap.json")
PLOT = os.path.join(HERE, "..", "plots", "m5_32_r29_1_ladder.png")
os.makedirs(os.path.dirname(PLOT), exist_ok=True)

# the file of record per rung: the converged bare pass, its post-correction, its Lanczos ghost read
RUNGS = [
    dict(tag="rad_pin_d0.3_w25_c0_n32_L48", n=32, L=48.0, h=1.5, gap="gap_rad_pin_d0.3_w25_c0_n32_L48_pass2.json", post="schurpost_rad_pin_d0.3_w25_c0_n32_L48_pass2.json", ghost="ghost_rad_pin_d0.3_w25_c0_n32_L48.json", role="continuity"),
    dict(tag="rad_pin_d0.3_w25_c0_n36_L54", n=36, L=54.0, h=1.5, gap="gap_rad_pin_d0.3_w25_c0_n36_L54_bare.json", post="schurpost_rad_pin_d0.3_w25_c0_n36_L54.json", ghost="ghost_rad_pin_d0.3_w25_c0_n36_L54.json", role="ladder"),
    dict(tag="rad_pin_d0.3_w25_c0_n48_L72", n=48, L=72.0, h=1.5, gap="gap_rad_pin_d0.3_w25_c0_n48_L72_bare2.json", post="schurpost_rad_pin_d0.3_w25_c0_n48_L72_pass2.json", ghost="ghost_rad_pin_d0.3_w25_c0_n48_L72.json", role="ladder"),
    dict(tag="rad_pin_d0.3_w25_c0_n72_L108", n=72, L=108.0, h=1.5, gap="gap_rad_pin_d0.3_w25_c0_n72_L108_bare.json", post="schurpost_rad_pin_d0.3_w25_c0_n72_L108.json", ghost="ghost_rad_pin_d0.3_w25_c0_n72_L108.json", role="ladder"),
    dict(tag="rad_pin_d0.3_w25_c0_n48_L48", n=48, L=48.0, h=1.0, gap="gap_rad_pin_d0.3_w25_c0_n48_L48_bare.json", post="schurpost_rad_pin_d0.3_w25_c0_n48_L48.json", ghost="ghost_rad_pin_d0.3_w25_c0_n48_L48.json", role="mesh check"),
]


def load(name):
    p = os.path.join(D, name)
    return json.load(open(p)) if os.path.exists(p) else None


rows = []
for r in RUNGS:
    g, s, gh = load(r["gap"]), load(r["post"]), load(r["ghost"])
    row = dict(r)
    if g:
        kp = g["K_positive"]
        m0 = kp["modes"][0]
        row.update(omega2_bare=m0["omega2_S2_over_2K"], K=m0["K"], S2=m0["S2"], r_mean=m0["r_weighted_mean"],
                   drift=float(kp["last_step_rel_change"][0]) if kp.get("last_step_rel_change") else None,
                   iterations=kp["iterations"], richardson=m0["richardson_rel_err_HD"], profile=m0["radial_weight_profile_(bins to L/2)"],
                   E=g["E"], maxG=g["max|G| spatial free"], static_stiff=g.get("static_stiff_directions"))
        hist = kp.get("lowest_eigenvalue_history_(every_1/20th_iteration)")
        if hist and len(hist) > 2:
            row["history_drift_last_tenth"] = float(abs(hist[-1] - hist[-3]) / max(1e-300, abs(hist[-1])))
    if s:
        row["omega2_cg"] = s["modes"][0]["omega2_cg"]
        row["omega2_cg_second"] = s["modes"][1]["omega2_cg"] if len(s["modes"]) > 1 else None
    if gh:
        m = gh["K_negative_ghost_read"]["modes"][0]
        row.update(ghost_S2=m["S2"], ghost_K=m["K"], ghost_r_mean=m["r_weighted_mean"])
    rows.append(row)

ladder = [r for r in rows if r.get("omega2_bare") is not None and r["h"] == 1.5]
mesh = next((r for r in rows if r["role"] == "mesh check" and r.get("omega2_bare") is not None), None)
out = {"rungs": rows}


def fit(xs, ys):
    xs, ys = np.array(xs), np.array(ys)
    A = np.vstack([np.ones_like(xs), xs]).T
    c, res, _, _ = np.linalg.lstsq(A, ys, rcond=None)
    pred = A @ c
    sig = float(np.sqrt(np.sum((ys - pred) ** 2) / max(1, len(xs) - 2))) if len(xs) > 2 else 0.0
    # the intercept's standard error from the fit alone
    if len(xs) > 2:
        cov = sig ** 2 * np.linalg.inv(A.T @ A)
        se0 = float(np.sqrt(cov[0, 0]))
    else:
        se0 = 0.0
    return dict(intercept=float(c[0]), slope=float(c[1]), residual=sig, intercept_se_fit=se0, n=len(xs))


for key in ("omega2_bare", "omega2_cg"):
    pts = [(r["L"], r[key]) for r in ladder if r.get(key) is not None]
    if len(pts) >= 2:
        Ls, ys = zip(*pts)
        f1 = fit([1 / L for L in Ls], ys); f2 = fit([1 / L ** 2 for L in Ls], ys)
        # the error budget on the intercept: fit se, the rungs' solver drift (max), the mesh difference
        drift = max((r.get("drift") or 0.0) * r[key] for r in ladder if r.get(key) is not None)
        mesh_diff = abs(mesh[key] - next(r[key] for r in ladder if r["L"] == mesh["L"])) if (mesh and mesh.get(key) is not None) else None
        err = float(np.sqrt(f1["intercept_se_fit"] ** 2 + drift ** 2 + (mesh_diff or 0.0) ** 2))
        word = "not resolved"
        if mesh_diff is not None and mesh_diff > abs(f1["intercept"]):
            word = "NOT RESOLVED at h 1.5 (the mesh difference exceeds the intercept)"
        elif f1["intercept"] > 2 * err:
            word = "GAPPED"
        elif abs(f1["intercept"]) < 2 * err:
            word = "GAPLESS (if the lowest mode is an exterior texture)"
        out[f"fit_{key}"] = {"points_(L, omega2)": pts, "fit_in_1_over_L": f1, "fit_in_1_over_L2": f2, "solver_drift_max": drift,
                             "mesh_difference_(h 1.0 vs 1.5, same box)": mesh_diff, "intercept_error_quadrature": err, "pre_registered_read": word}

json.dump(out, open(OUT, "w"), indent=1, default=float)

# the figure
fig, axes = plt.subplots(1, 2, figsize=(12.5, 5))
ax = axes[0]
for key, c, lab in (("omega2_bare", "tab:blue", "bare pencil"), ("omega2_cg", "tab:green", "with the K-null elimination")):
    pts = [(r["L"], r[key]) for r in ladder if r.get(key) is not None]
    if pts:
        Ls, ys = zip(*pts); ax.plot([1 / L for L in Ls], ys, "o", color=c, label=lab)
        fk = out.get(f"fit_{key}")
        if fk:
            x = np.linspace(0, max(1 / L for L in Ls) * 1.05, 50)
            ax.plot(x, fk["fit_in_1_over_L"]["intercept"] + fk["fit_in_1_over_L"]["slope"] * x, "-", color=c, lw=0.8, label=f"1/L fit: c0 = {fk['fit_in_1_over_L']['intercept']:.2e}")
            ax.plot(x, fk["fit_in_1_over_L2"]["intercept"] + fk["fit_in_1_over_L2"]["slope"] * x ** 2 * 0 + fk["fit_in_1_over_L2"]["slope"] * (x ** 2), "--", color=c, lw=0.8, label=f"1/L^2 fit: c0 = {fk['fit_in_1_over_L2']['intercept']:.2e}")
if mesh and mesh.get("omega2_bare") is not None:
    ax.plot(1 / mesh["L"], mesh["omega2_bare"], "s", color="tab:red", label="h 1.0 mesh check (shell branch)")
ax.axhline(0, color="k", lw=0.5)
ax.set_xlabel("1 / L"); ax.set_ylabel("lowest K-positive omega^2"); ax.set_title("R29-1: the ladder at h 1.5"); ax.legend(fontsize=7)
ax = axes[1]
for r in rows:
    if r.get("ghost_S2") is not None:
        ax.plot(r["L"] if r["h"] == 1.5 else r["L"] + 2, r["ghost_S2"], "v", color="tab:red" if r["h"] == 1.5 else "tab:orange", ms=8)
        ax.annotate(f"h {r['h']}, r {r['ghost_r_mean']:.0f}", (r["L"] if r["h"] == 1.5 else r["L"] + 2, r["ghost_S2"]), fontsize=7, textcoords="offset points", xytext=(4, 4))
ax.axhline(0, color="k", lw=0.5)
ax.set_xlabel("L"); ax.set_ylabel("S2 of the K-negative witness (unit norm)"); ax.set_title("R29-1: the GHOST witness per rung")
plt.tight_layout(); plt.savefig(PLOT, dpi=130)
print(json.dumps({k: v for k, v in out.items() if k.startswith("fit_")}, indent=1, default=float))
print("wrote", os.path.relpath(OUT, HERE), os.path.relpath(PLOT, HERE))
