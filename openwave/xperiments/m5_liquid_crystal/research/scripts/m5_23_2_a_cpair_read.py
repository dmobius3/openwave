"""M5.23.2 arm (1) — the tracer's first physics read: the M5.21.10 C-decay
ejection pair (the two-loop identity question the census could not close).

The M5.21.10 record (m5_21_10_note.md): the C (tau-candidate) free decay on
the 64^3 arena shows a PAIRED SYMMETRIC ejection at t ~ 80 (two off-center
features, directions 108.9 deg apart), but the audit refuted the compact-
fragment picture — the census fragments were 1-cell filament doublets whose
identity was CUT-SENSITIVE (blob components inside a radius cut). The
tracer replaces the census read: connectivity-assembled lines with explicit
Euler-characteristic closure verdicts and no arena radii.

READ QUESTIONS (pre-registered, m5_23_2_task_details TASK PLANNING T-phys —
a measured read, no pass/fail): per snapshot t = 10..150 of
m5_21_10_ev_C_free64.npz: how many defect lines; are the ejecta CLOSED
loops (the two-released-loops picture needs closed) or boundary run-outs /
open filaments; A-control (electron holds) as the null.

Outputs: data/m5_23_2_cpair_read.json + plots/m5_23_2_tracer_panel.png.
Runtime ~2 min (16 x 64^3 eigvalsh + the instrument-validation traces).
"""

import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
PLOTS = os.path.join(HERE, "..", "plots")
sys.path.insert(0, HERE)
import m5_23_2_tracer as tracer  # noqa: E402

C_NPZ = os.path.join(DATA, "m5_21_10_ev_C_free64.npz")
A_NPZ = os.path.join(DATA, "m5_21_10_ev_A_free64.npz")
E_NPZ = os.path.join(DATA, "m5_21_9_fixedj_conj_om0.2_end.npz")

zc = np.load(C_NPZ)
dt = float(zc["dt"])
h = float(zc["h"])

rows = []
snaps = {}
for key in [f"M_it{it}" for it in range(400, 6001, 400)] + ["M"]:
    it = int(key[4:]) if key.startswith("M_it") else 6000
    t = it * dt
    if key == "M" and rows and rows[-1]["t"] == t:
        continue  # final M duplicates it6000
    r = tracer.trace(zc[key], h=h)
    ln_sum = [
        {k: v for k, v in ln.items() if k not in ()}
        for ln in r["lines"]
    ]
    rows.append(
        {
            "t": t,
            "n_lines": len(r["lines"]),
            "bulk_split": r["bulk_split"],
            "thr": r["thr"],
            "verdicts": [(ln["n_vox"], ln["verdict"]) for ln in r["lines"]],
            "lines": ln_sum,
        }
    )
    snaps[t] = r
    print(
        f"t = {t:6.1f}: {len(r['lines'])} line(s): "
        + ", ".join(f"{ln['n_vox']}vox {ln['verdict']}" for ln in r["lines"][:6])
        + (" ..." if len(r["lines"]) > 6 else "")
    )

# the A-control null (the electron HOLDS in free dynamics)
za = np.load(A_NPZ)
r_a = tracer.trace(za["M"], h=h)
print(
    f"A control (t = 150): {len(r_a['lines'])} line(s): "
    + ", ".join(f"{ln['n_vox']}vox {ln['verdict']}" for ln in r_a["lines"])
)

# instrument-validation traces for the panel
ze = np.load(E_NPZ)
r_e = tracer.trace(ze["M"][..., 1:, 1:], h=1.5)

out = {
    "source": os.path.basename(C_NPZ),
    "read": "tracer line census per snapshot (m5_23_2 arm 1 first physics read)",
    "rows": [{k: v for k, v in row.items() if k != "lines"} for row in rows],
    "rows_full": rows,
    "a_control_final": {
        "n_lines": len(r_a["lines"]),
        "verdicts": [(ln["n_vox"], ln["verdict"]) for ln in r_a["lines"]],
    },
    "electron_endpoint_validation": {
        "n_lines": len(r_e["lines"]),
        "verdicts": [(ln["n_vox"], ln["verdict"]) for ln in r_e["lines"]],
    },
}
os.makedirs(DATA, exist_ok=True)
with open(os.path.join(DATA, "m5_23_2_cpair_read.json"), "w") as f:
    json.dump(out, f, indent=1)

# ---------------- the evidence panel ----------------
fig = plt.figure(figsize=(16, 12))
gs = fig.add_gridspec(3, 4, hspace=0.35, wspace=0.3)

# row 1: instrument validation — electron endpoint + the C seed state
ax = fig.add_subplot(gs[0, 0])
sm = r_e["split"][:, 16, :]
im = ax.imshow(sm.T, origin="lower", cmap="viridis")
ax.contour((r_e["mask"][:, 16, :]).T, levels=[0.5], colors="r", linewidths=1)
ax.set_title("electron endpoint: split map y-mid\n(red = detected core; polar rods)")
plt.colorbar(im, ax=ax, shrink=0.8)

ax = fig.add_subplot(gs[0, 1], projection="3d")
for ln in r_e["lines"]:
    pass
cols = plt.cm.tab10.colors
idx_all = np.argwhere(r_e["mask"])
lab = np.zeros(len(idx_all), int)
for li, ln in enumerate(r_e["lines"]):
    lo, hi = np.array(ln["bbox"][0]), np.array(ln["bbox"][1])
    inbb = np.all((idx_all >= lo) & (idx_all <= hi), axis=1)
    lab[inbb] = li + 1
for li in range(1, lab.max() + 1):
    p = idx_all[lab == li]
    ax.scatter(p[:, 0], p[:, 1], p[:, 2], s=6, color=cols[(li - 1) % 10])
ax.set_title(f"electron: {len(r_e['lines'])} traced line(s)\n(the Stage D polar rods, FOUND)")
ax.set_xlim(0, 32), ax.set_ylim(0, 32), ax.set_zlim(0, 32)

# row 1 cont: C at t = 10 (pre-decay) and t = 80 (the split)
for col, t_show in ((2, 10.0), (3, 80.0)):
    ax = fig.add_subplot(gs[0, col])
    r_s = snaps[t_show]
    sm = r_s["split"][:, 32, :]
    im = ax.imshow(sm.T, origin="lower", cmap="viridis")
    ax.contour((r_s["mask"][:, 32, :]).T, levels=[0.5], colors="r", linewidths=1)
    ax.set_title(f"C state t = {t_show:.0f}: split y-mid\n({r_s['thr']:.3f} thr)")
    plt.colorbar(im, ax=ax, shrink=0.8)

# row 2: the line census vs time
ax = fig.add_subplot(gs[1, :2])
ts = [row["t"] for row in rows]
ax.plot(ts, [row["n_lines"] for row in rows], "o-", label="traced lines")
ax.plot(
    ts,
    [sum(1 for _, v in row["verdicts"] if v == "closed-loop") for row in rows],
    "s-",
    label="closed loops",
)
ax.plot(
    ts,
    [sum(1 for _, v in row["verdicts"] if v == "boundary") for row in rows],
    "^-",
    label="boundary run-outs",
)
ax.plot(
    ts,
    [sum(1 for _, v in row["verdicts"] if v == "open") for row in rows],
    "v-",
    label="open filaments",
)
ax.axvline(80, color="gray", ls="--", lw=1)
ax.text(80, ax.get_ylim()[1] * 0.95, " t=80 split (M5.21.10)", fontsize=8, va="top")
ax.set_xlabel("t (tau)"), ax.set_ylabel("count")
ax.set_title("C-decay line census (tracer, connectivity + Euler-chi closure)")
ax.legend(fontsize=8)

# row 2 cont: total core voxels vs time (defect material budget)
ax = fig.add_subplot(gs[1, 2:])
ax.plot(ts, [sum(n for n, _ in row["verdicts"]) for row in rows], "o-", color="k")
ax.axvline(80, color="gray", ls="--", lw=1)
ax.set_xlabel("t (tau)"), ax.set_ylabel("core voxels (all lines)")
ax.set_title("detected core material vs time")

# row 3: 3D scatter of the components at three times
for col, t_show in ((0, 60.0), (1, 100.0), (2, 150.0)):
    ax = fig.add_subplot(gs[2, col], projection="3d")
    r_s = snaps[t_show]
    labels_arr = np.zeros(r_s["mask"].shape, int)
    from scipy import ndimage

    labels_arr, _ = ndimage.label(r_s["mask"], structure=tracer.STRUCT26)
    for li, ln in enumerate(r_s["lines"][:8]):
        p = np.argwhere(labels_arr == ln["id"])
        ax.scatter(p[:, 0], p[:, 1], p[:, 2], s=5, color=cols[li % 10],
                   label=f"{ln['n_vox']}vox {ln['verdict']}")
    ax.set_title(f"C t = {t_show:.0f}")
    ax.set_xlim(0, 64), ax.set_ylim(0, 64), ax.set_zlim(0, 64)
    ax.legend(fontsize=6, loc="upper left")

# row 3 cont: A control final
ax = fig.add_subplot(gs[2, 3], projection="3d")
labels_a, _ = ndimage.label(r_a["mask"], structure=tracer.STRUCT26)
for li, ln in enumerate(r_a["lines"][:8]):
    p = np.argwhere(labels_a == ln["id"])
    ax.scatter(p[:, 0], p[:, 1], p[:, 2], s=5, color=cols[li % 10],
               label=f"{ln['n_vox']}vox {ln['verdict']}")
ax.set_title("A control t = 150 (electron holds)")
ax.set_xlim(0, 64), ax.set_ylim(0, 64), ax.set_zlim(0, 64)
ax.legend(fontsize=6, loc="upper left")

fig.suptitle(
    "M5.23.2 arm (1): the disclination-line tracer — instrument validation + "
    "the M5.21.10 C-pair read (m5_21_10_ev_C_free64 snapshots)",
    fontsize=12,
)
os.makedirs(PLOTS, exist_ok=True)
fig.savefig(os.path.join(PLOTS, "m5_23_2_tracer_panel.png"), dpi=110, bbox_inches="tight")
print("saved plots/m5_23_2_tracer_panel.png + data/m5_23_2_cpair_read.json")
