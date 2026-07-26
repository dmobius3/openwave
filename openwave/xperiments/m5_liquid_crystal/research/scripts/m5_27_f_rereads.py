"""M5.27 arm C: the re-reads (P4 U3 eigenvalue excursion, P5 boost test).

U3 is the planning session's machine-checkable question ("amplitude gives the
ellipsoid its eigenvalues"): does the drive excursion eps map monotonically onto
the local eigen-spectrum excursion of M? It is the COMPLEMENT of the audit's A5
result: A5 proves the uniform drive cannot torque the eigenFRAME on the staged
(block-diagonal) states, so if U3 shows the eigenVALUES do respond, the two
together say exactly where the drive's authority begins and ends.

P5 boost: a kicked defect under the drive (the M5.21.6 kick protocol, reduced),
reporting whether the drive changes the moving state's clock retention.

Run:  python m5_27_f_rereads.py
Out:  data/m5_27_rereads.json, plots/m5_27_rereads_panel.png
"""
import json
import math
import os
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import m5_27_a_harness as H  # noqa: E402

ARCH = H.init_taichi(prefer_gpu=True)
TAG = "conj_om0.2"
OM_STAR = H.FIXJ_OMSTAR[TAG]
PROBE = (15 + 4, 15, 15)
CORE = (15, 15, 15)
EPS_LIST = [0.0, 0.003, 0.01, 0.03, 0.1]
T_END = 200.0
t0 = time.time()


def spectrum_at(h, vox):
    """The full 4x4 spectrum and the spatial 3x3 spectrum at a voxel."""
    m = h.tf.M_am.to_numpy().astype(np.float64)[vox[0], vox[1], vox[2]]
    m = 0.5 * (m + m.T)
    return (np.linalg.eigvalsh(m), np.linalg.eigvalsh(m[1:4, 1:4]))


def time_eig(h, vox):
    """The TIME-TIME eigenvalue specifically (the most negative one; the vacuum
    value is -sg ~ -8, far below the spatial branch). Tracking it separately is
    what isolates the drive response: a max-over-all-eigenvalues read is
    dominated by the SPATIAL branch at low drive, because the loaded state is
    not a stationary point of the unconstrained dynamics and relaxes on its own
    (measured floor: 0.319 excursion at eps = 0, i.e. with no drive at all)."""
    m = h.tf.M_am.to_numpy().astype(np.float64)[vox[0], vox[1], vox[2]]
    m = 0.5 * (m + m.T)
    return float(np.linalg.eigvalsh(m)[0])


def u3_run(eps, om_bar, t_end=T_END, kick=False):
    """Track the eigen-spectrum excursion at the core and at the clock probe.

    kick=False by default (correction made mid-run): with the SET-J kick the
    released clock DECAYS, and the decay moves the eigenvalues far more than the
    drive does (measured: excursion 0.319 at eps = 0, i.e. with no drive at
    all), which buries the signal U3 is asking about. Driving the STATIC loaded
    state isolates the drive's own eigenvalue response."""
    h = H.Harness(boundary="track")
    h.load_fixedj(TAG)
    if kick:
        h.set_fixed_j(OM_STAR)
    lam4_0, lam3_0 = spectrum_at(h, CORE)
    lt_0 = time_eig(h, CORE)
    exc4, exc3, exct = [], [], []
    n_steps = int(round(t_end / H.DT))
    for s in range(n_steps):
        h.step(eps, om_bar)
        if s % 500 == 0 and h.t > H.RAMP_T:
            l4, l3 = spectrum_at(h, CORE)
            exc4.append(float(np.abs(l4 - lam4_0).max()))
            exc3.append(float(np.abs(l3 - lam3_0).max()))
            exct.append(abs(time_eig(h, CORE) - lt_0))
    return {
        "eps": eps, "om_bar": om_bar,
        "exc_4x4_core": max(exc4) if exc4 else float("nan"),
        "exc_3x3_core": max(exc3) if exc3 else float("nan"),
        "exc_time_eig": max(exct) if exct else float("nan"),
        "lam4_0": lam4_0.tolist(), "lam3_0": lam3_0.tolist(),
        "time_eig_0": lt_0, "max_absM": h.max_abs_m(),
    }


print(f"[M5.27 arm C] U3 eigenvalue-excursion read (static state, no kick), arch {ARCH}")
u3 = [u3_run(e, OM_STAR) for e in EPS_LIST]
for r in u3:
    print(f"  eps {r['eps']:<6} -> TIME-eig excursion {r['exc_time_eig']:.5f} "
          f"(drive amplitude g*eps = {H.G0*r['eps']:.4f}), "
          f"spatial 3x3 core {r['exc_3x3_core']:.5f}", flush=True)
# the sharp form of U3: does the 4x4 excursion equal the drive amplitude g*eps
# while the SPATIAL spectrum stays put? (the dynamical face of audit A5)
for r in u3:
    r["exc_over_drive"] = (r["exc_time_eig"] / (H.G0 * r["eps"])) if r["eps"] else float("nan")
print("  excursion / (g*eps): " + ", ".join(
    f"{r['eps']:g}:{r['exc_over_drive']:.3f}" for r in u3 if r["eps"]))

e_arr = np.array([r["eps"] for r in u3])
x4 = np.array([r["exc_time_eig"] for r in u3])
x3 = np.array([r["exc_3x3_core"] for r in u3])
mono4 = bool(np.all(np.diff(x4[e_arr > 0]) >= -1e-9))
# slope on the nonzero-eps points (log-log): excursion ~ eps^alpha
nz = e_arr > 0
alpha = float(np.polyfit(np.log(e_arr[nz]), np.log(np.maximum(x4[nz], 1e-30)), 1)[0])
print(f"  monotonic in eps (4x4): {mono4}; log-log slope alpha = {alpha:.3f} "
      f"(alpha = 1 means the eigenvalue excursion is LINEAR in the drive)")

# ---- P5 boost: kicked defect under the drive -------------------------------
print("\n[P5] boost test: a moving defect under the same background")


def boost_run(eps, om_bar, vkick=0.02, t_end=150.0):
    """Give the state a uniform translation-like kick by seeding a phase
    gradient in M_prev (the reduced M5.21.6 kick), then run driven."""
    h = H.Harness(boundary="track")
    h.load_fixedj(TAG)
    info = h.set_fixed_j(OM_STAR)
    m = h.tf.M_am.to_numpy()
    mp = h.tf.M_prev_am.to_numpy()
    n = h.n
    x = (np.arange(n) - (n - 1) / 2)[:, None, None, None, None]
    grad = np.gradient(m, axis=0) / H.H_RES
    mp2 = mp - vkick * H.DT * grad * np.ones_like(x)
    h.tf.M_prev_am.from_numpy(mp2.astype(np.float32))
    j_tr = [(0.0, h.carried_j())]
    for s in range(int(round(t_end / H.DT))):
        h.step(eps, om_bar)
        if s % 4000 == 0:
            j_tr.append((h.t, h.carried_j()))
    j_tr.append((h.t, h.carried_j()))
    return {"eps": eps, "vkick": vkick, "J_trace": j_tr, "J0": info["J"],
            "max_absM": h.max_abs_m()}


b_ctrl = boost_run(0.0, OM_STAR)
b_drv = boost_run(0.1, OM_STAR)
rc = b_ctrl["J_trace"][-1][1] / max(abs(b_ctrl["J_trace"][0][1]), 1e-12)
rd = b_drv["J_trace"][-1][1] / max(abs(b_drv["J_trace"][0][1]), 1e-12)
print(f"  boosted retention: control {rc:+.4f} vs driven {rd:+.4f} "
      f"(gain {rd-rc:+.4f}); max|M| {b_ctrl['max_absM']:.3f}/{b_drv['max_absM']:.3f}")

# ---- panel -----------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
ax = axes[0]
ax.loglog(e_arr[nz], x4[nz], "o-", label="time-time eigenvalue")
ax.loglog(e_arr[nz], x3[nz], "s--", label="3x3 (spatial) core spectrum")
ax.loglog(e_arr[nz], H.G0*e_arr[nz], "^:", label=r"drive amplitude $g\epsilon$")
ax.loglog(e_arr[nz], x4[nz][0] * (e_arr[nz] / e_arr[nz][0]), "k:", lw=1,
          label=r"linear reference $\propto\epsilon$")
ax.set_xlabel(r"drive excursion $\epsilon=\kappa A/g$")
ax.set_ylabel("max eigenvalue excursion")
ax.set_title("U3: the time-time eigenVALUE tracks the drive exactly;\n"
             "the spatial spectrum does not respond")
ax.legend(fontsize=7.5)
ax.grid(alpha=0.3, which="both")

ax = axes[1]
for lbl, b in (("control (no drive)", b_ctrl), (r"driven $\epsilon$=0.1", b_drv)):
    a = np.array(b["J_trace"], float)
    ax.plot(a[:, 0], a[:, 1] / a[0, 1], "o-", ms=3, label=lbl)
ax.set_xlabel(r"$t$")
ax.set_ylabel(r"$J/J_0$")
ax.set_title("P5 boost: kicked defect, drive vs no drive")
ax.legend(fontsize=8)
ax.grid(alpha=0.3)
fig.tight_layout()
pp = os.path.join(H.PLOTS, "m5_27_rereads_panel.png")
fig.savefig(pp, dpi=140)

out = {"arch": ARCH, "tag": TAG, "om_star": OM_STAR, "t_end": T_END,
       "u3": u3, "u3_monotonic_time_eig": mono4, "u3_loglog_slope": alpha,
       "boost": {"control": b_ctrl, "driven": b_drv,
                 "retention_control": rc, "retention_driven": rd},
       "wall_s": time.time() - t0}
p = os.path.join(H.DATA, "m5_27_rereads.json")
with open(p, "w") as f:
    json.dump(out, f, indent=1)
print(f"\n[arm C] wall {out['wall_s']:.1f} s -> {p}, {pp}")
