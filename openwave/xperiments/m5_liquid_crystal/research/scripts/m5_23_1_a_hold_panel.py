"""M5.23.1 — the hold-evidence panel: field-state cross-sections + the
carried-charge trace of the PRODUCTION fixed-J evolution (the simulation-
prints rule: seed + endpoint states shown as direct evidence the run ran).

Reruns the selftest S4 hold arm (production kernels only, no f64 reference):
the omega = 0.2 conjugation rung, SET-J kick (Mdot(0) = omega*.a0), 400
leapfrog steps at dt = 0.005 on the 32^3 research arena. Saves:
    plots/m5_23_1_hold_panel.png   (seed/end lambda_max cross-sections,
                                    J_self(t), the energy ledger)
    data/m5_23_1_selftest.json     (the summary record of the gated run)

USAGE (repo root):
    python -m openwave.xperiments.m5_liquid_crystal.research.scripts.m5_23_1_a_hold_panel
"""

import importlib.util
import json
import os

import matplotlib
import numpy as np
import taichi as ti

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ti.init(arch=ti.cpu, log_level=ti.WARN)

import openwave.xperiments.m5_liquid_crystal.engine2_pde as pde  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
PLOTS = os.path.join(HERE, "..", "plots")

spec = importlib.util.spec_from_file_location("ref", os.path.join(HERE, "m5_21_3_a_4d.py"))
ref = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ref)

cfg = ref.base_cfg(s=-1.0, tag="p1_s-1")
N, H, SG, DELTA, RENV = cfg["n"], cfg["h"], cfg["sg"], cfg["delta"], cfg["renv"]
W1 = pde.W1_SPECTRAL
DT = 0.005
N_HOLD = 400


class DuckField:
    def __init__(self, n):
        self.nx = self.ny = self.nz = n
        self.M_am = ti.Matrix.field(4, 4, dtype=ti.f32, shape=(n, n, n))
        self.M_prev_am = ti.Matrix.field(4, 4, dtype=ti.f32, shape=(n, n, n))
        self.M_new_am = ti.Matrix.field(4, 4, dtype=ti.f32, shape=(n, n, n))
        self.Md_am = ti.Matrix.field(4, 4, dtype=ti.f32, shape=(n, n, n))
        self.curv_flux_x = ti.Matrix.field(4, 4, dtype=ti.f32, shape=(n, n, n))
        self.curv_flux_y = ti.Matrix.field(4, 4, dtype=ti.f32, shape=(n, n, n))
        self.curv_flux_z = ti.Matrix.field(4, 4, dtype=ti.f32, shape=(n, n, n))
        self.fire_partials = ti.field(dtype=ti.f32, shape=(3, n))

    def load(self, m_np):
        self.M_am.from_numpy(m_np.astype(np.float32))
        self.M_prev_am.from_numpy(m_np.astype(np.float32))
        self.M_new_am.from_numpy(m_np.astype(np.float32))

    def swap_matrix_buffers(self):
        self.M_prev_am.copy_from(self.M_am)
        self.M_am.copy_from(self.M_new_am)


def lam_max_slice(m_np):
    """lambda_max of the spatial 3x3 block on the y = mid meridional plane."""
    sl = m_np[:, N // 2, :, 1:4, 1:4]
    return np.linalg.eigvalsh(sl)[..., -1]


Z = np.load(os.path.join(DATA, "m5_21_9_fixedj_conj_om0.2_end.npz"))
with open(os.path.join(DATA, "m5_21_9_fixedj_conj_om0.2.json")) as f:
    rec = json.load(f)
M0 = Z["M"].astype(np.float64)
om_star = rec["final"]["omega_star_final"]

tf = DuckField(N)
tf.load(M0)
info = pde.set_fixed_j(tf, H, RENV, om_star, DT, shell=2)
seed_slice = lam_max_slice(tf.M_am.to_numpy().astype(np.float64))

ts, js, es = [], [], []
for step in range(N_HOLD + 1):
    if step % 20 == 0:
        js.append(pde.read_carried_j(tf, H, RENV, DT))
        ts.append(step * DT)
    if step % 100 == 0:
        m_np = tf.M_am.to_numpy().astype(np.float64)
        mt_np = (m_np - tf.M_prev_am.to_numpy().astype(np.float64)) / DT
        eu, ev = ref.e_parts(m_np, cfg)
        es.append((step * DT, float(eu) + float(ev) + 0.5 * H**3 * float(np.sum(mt_np * mt_np))))
    if step < N_HOLD:
        pde.compute_eta_flux(tf, 0, H)
        pde.evolve_M_eta_start(tf, DT, H)
        pde.compute_eta_flux(tf, 1, H)
        pde.evolve_M_eta_finish(tf, DT, H, SG, DELTA, W1)
        tf.swap_matrix_buffers()

end_slice = lam_max_slice(tf.M_am.to_numpy().astype(np.float64))

fig, axes = plt.subplots(2, 2, figsize=(10.5, 8.5))
for ax, sl, ttl in (
    (axes[0, 0], seed_slice, "t = 0 (post SET-J kick)"),
    (axes[0, 1], end_slice, f"t = {N_HOLD * DT:g} (endpoint)"),
):
    im = ax.imshow(sl.T, origin="lower", cmap="viridis")
    ax.set_title(f"lambda_max spatial block, y-mid plane\n{ttl}")
    ax.set_xlabel("x [cell]")
    ax.set_ylabel("z [cell]")
    fig.colorbar(im, ax=ax, shrink=0.85)
axes[1, 0].plot(ts, js, "o-", ms=3)
axes[1, 0].axhline(om_star, ls="--", lw=1, label=f"omega* = {om_star:.5f}")
axes[1, 0].set_xlabel("t [tau]")
axes[1, 0].set_ylabel("carried J_self = <Mdot, a0(M)>")
axes[1, 0].set_title("the carried isorotation charge (hold)")
axes[1, 0].legend()
e_t, e_v = zip(*es)
axes[1, 1].plot(e_t, e_v, "s-", ms=4)
axes[1, 1].set_xlabel("t [tau]")
axes[1, 1].set_ylabel("E_u + V4 + KE")
axes[1, 1].set_title("the energy ledger (f32 floor scale)")
fig.suptitle(
    "M5.23.1 production fixed-J hold: omega = 0.2 rung, 32^3 research arena, dt = 0.005",
    fontsize=11,
)
fig.tight_layout()
os.makedirs(PLOTS, exist_ok=True)
out_png = os.path.join(PLOTS, "m5_23_1_hold_panel.png")
fig.savefig(out_png, dpi=110)
print("saved", out_png)

summary = {
    "arena": {"n": N, "h": H, "sg": SG, "delta": DELTA, "renv": RENV, "dt": DT},
    "set_j": {k: info[k] for k in ("kin", "J", "om_star")},
    "hold": {
        "steps": N_HOLD,
        "J_self_start": js[0],
        "J_self_end": js[-1],
        "retention": js[-1] / js[0],
        "E_start": es[0][1],
        "E_end": es[-1][1],
    },
    "provenance": "m5_23_1_fixedj_engine_selftest.py = the gate record (ALL 12 GREEN 2026-07-24)",
}
with open(os.path.join(DATA, "m5_23_1_selftest.json"), "w") as f:
    json.dump(summary, f, indent=1)
print(json.dumps(summary["hold"], indent=1))
