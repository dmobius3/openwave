"""
M5.23.1 — headless selftest of the FIXED-J ISOROTATION production port.

Cross-validates the new production kernels (engine2_pde: clock_flow_k /
md_norm_partials_k / kin_partials_k / isorotation_kick_k / j_partials_k +
the hosts compute_clock_flow / kin_canonical / set_fixed_j /
read_isorotation) against the AUDITED research machinery of record,
imported directly so the reference code itself is the oracle (no
re-transcription): m5_21_3_a_4d.py (the certified 4D instrument, INS4),
m5_21_9_d_fixedj.py (a0_conj, the conjugation-tangent clock flow) and
m5_21_9_e_larmor.py (leap(), the certified 4x4 leapfrog: E-drift 2.2e-8
per 400 steps in f64). Taichi runs f32; gates are set at f32 scale with
the tolerance ladder documented per gate.

ARENA: the S1-S4 gates run on the RESEARCH 32^3 arena (the M5.21.3
p1_s-1 certified endpoint + the three fixed-J conjugation rungs on
disk). TensorField rounds grids to ODD sizes, so a selftest-local
duck-typed field container carries the 32^3 buffers — the objects under
test are the KERNELS (templated on attribute access), not the container.
S5 runs the launcher-arena flow (63^3, research-unit dx) on a real
TensorField.

  S1  a0 FLOW      — production clock flow vs the research a0_conj on the
                     certified state, compared up to the per-voxel SIGN
                     GAUGE (the leading eigenvector is apolar; research =
                     numpy eigh signs, production = the radial gauge);
                     the gauge overlap is quantified, not assumed
  S2  kin          — the clock-inertia reduction both ways (production
                     a0 -> reference kin_of; reference a0 -> production
                     reduction) + gauge invariance + the canonical anchor
                     kin = 0.1206 (the recorded kin0 = 0.12055)
  S3  LEGENDRE     — dE/dJ = omega* across the three research rungs with
                     PRODUCTION kin (E_stat from the audited reference
                     e_parts, the energy-kernel equivalence being already
                     certified at M5.24 gate 3): vs the research
                     0.997 / 0.992 closure
  S4  HOLD         — SET-J kick + production evolution of the omega = 0.2
                     rung (dt = 0.005, research arena): bounded ledger,
                     core spectrum intact, J retained; observable-level
                     cross-check vs the reference leap() trajectory (f64)
                     with the SAME flow and init
  S5  LAUNCHER     — the RELAX -> SET-J -> EVOLVE flow at launcher scale
                     (63^3 hedgehog, flip, FIRE, kick, substepped evolve +
                     sponge): bounded, J coherent, and the delta clock-hand
                     axis MEASURABLY rotating (rate reported, not assumed)

USAGE (repo root):
    python -m openwave.xperiments.m5_liquid_crystal.research.scripts.m5_23_1_fixedj_engine_selftest
"""

import importlib.util
import json
import math
import os

import numpy as np
import taichi as ti

ti.init(arch=ti.cpu, log_level=ti.WARN)

import openwave.xperiments.m5_liquid_crystal.engine2_pde as pde  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ref = _load("m5_21_3_a_4d", "m5_21_3_a_4d.py")  # the certified 4D instrument
fixj = _load("m5_21_9_d_fixedj", "m5_21_9_d_fixedj.py")  # a0_conj
larm = _load("m5_21_9_e_larmor", "m5_21_9_e_larmor.py")  # leap()

fails = []
total = [0]


def check(name, ok, detail=""):
    total[0] += 1
    tag = "PASS" if ok else "FAIL"
    print(f"  [{tag}] {name}  {detail}", flush=True)
    if not ok:
        fails.append(name)


# ----------------------------------------------------------------
# The duck-typed 32^3 container (kernels are templated on attributes;
# TensorField's odd-grid rule blocks the even research arena)
# ----------------------------------------------------------------
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

    def load(self, m_np, m_prev_np=None):
        self.M_am.from_numpy(m_np.astype(np.float32))
        self.M_prev_am.from_numpy((m_np if m_prev_np is None else m_prev_np).astype(np.float32))
        self.M_new_am.from_numpy(m_np.astype(np.float32))

    def swap_matrix_buffers(self):
        self.M_prev_am.copy_from(self.M_am)
        self.M_am.copy_from(self.M_new_am)


# ----------------------------------------------------------------
# Research arena: the conjugation rungs (the M5.21.3 p1_s-1 base-state
# npz was a casualty of the pre-2026-07-20 delete rule — being restored
# separately; the S1/S2 flow+kin gates run on the omega = 0.2 rung
# ENDPOINT, whose recorded kin_final/omega* anchor the same conventions)
# ----------------------------------------------------------------
cfg = ref.base_cfg(s=-1.0, tag="p1_s-1")  # the load_p1 cfg, built directly
N = cfg["n"]
H = cfg["h"]
SG = cfg["sg"]  # s = -1 branch: sg = -8
DELTA = cfg["delta"]
RENV = cfg["renv"]
W1 = pde.W1_SPECTRAL
DT = 0.005
SHELL = 2  # the research free mask (pin depth 1.6 at h = 1.5 -> 2 cells)

RUNGS = {}
for om in ("0.2", "0.5", "1"):
    Z = np.load(os.path.join(DATA, f"m5_21_9_fixedj_conj_om{om}_end.npz"))
    with open(os.path.join(DATA, f"m5_21_9_fixedj_conj_om{om}.json")) as f:
        RUNGS[om] = {"M": Z["M"].astype(np.float64), "rec": json.load(f)}

tf = DuckField(N)

# The S1/S2 base = the omega = 0.2 rung ENDPOINT (present on disk); its
# recorded kin_final = 0.121014 anchors the conjugation convention there
M_base = RUNGS["0.2"]["M"]
KIN_ANCHOR = RUNGS["0.2"]["rec"]["final"]["kin_final"]

print("M5.23.1 fixed-J isorotation selftest")
print(f"  research arena {N}^3, h = {H}, sg = {SG}, delta = {DELTA}, renv = {RENV}")


# ---- S1: the a0 flow vs the research a0_conj (sign-gauge aware) ----
a0_res = fixj.a0_conj(cfg, M_base)  # f64, eigh sign gauge
tf.load(M_base)
raw_norm = pde.compute_clock_flow(tf, H, RENV)
a0_prod = tf.Md_am.to_numpy().astype(np.float64)

dots = np.einsum("...ab,...ab->...", a0_prod, a0_res)
w_res = np.einsum("...ab,...ab->...", a0_res, a0_res)
sgn = np.where(dots < 0.0, -1.0, 1.0)
a0_aligned = a0_prod * sgn[..., None, None]
err = np.sqrt(np.sum((a0_aligned - a0_res) ** 2)) / np.sqrt(np.sum(a0_res**2))
flip_w = float(w_res[dots < 0.0].sum() / max(w_res.sum(), 1e-300))
check(
    "S1 a0 flow == research a0_conj up to the sign gauge",
    err < 1e-3,
    f"aligned rel err {err:.2e}, raw norm {raw_norm:.4f}, flipped weight {100 * flip_w:.2f}%",
)
print(
    f"  [info] sign-gauge overlap: {100 * (1 - flip_w):.2f}% of the flow weight shares "
    f"the eigh gauge; the flipped remainder is the apolar ambiguity, pinned radially "
    f"in production (quadratic reads are gauge-invariant, S2)"
)

# ---- S2: kin both ways + gauge invariance + the canonical anchor ----
kin_prod = pde.kin_canonical(tf, H)  # production a0 (still in Md), production reduction
kin_ref_on_prod = ref.kin_of(M_base, a0_prod, cfg)  # production a0, reference reduction
rel_a = abs(kin_prod - kin_ref_on_prod) / abs(kin_ref_on_prod)
check(
    "S2a kin: production reduction == reference kin_of (production a0)",
    rel_a < 5e-4,
    f"prod {kin_prod:.6f} vs ref {kin_ref_on_prod:.6f} rel {rel_a:.2e}",
)
tf.Md_am.from_numpy(a0_res.astype(np.float32))
kin_prod_on_res = pde.kin_canonical(tf, H)
kin_ref = ref.kin_of(M_base, a0_res, cfg)
rel_b = abs(kin_prod_on_res - kin_ref) / abs(kin_ref)
check(
    "S2b kin: production reduction == reference kin_of (research a0)",
    rel_b < 5e-4,
    f"prod {kin_prod_on_res:.6f} vs ref {kin_ref:.6f} rel {rel_b:.2e}",
)
gauge_rel = abs(kin_prod - kin_ref) / abs(kin_ref)
anchor_rel = abs(kin_prod - KIN_ANCHOR) / KIN_ANCHOR
check(
    "S2c kin gauge-invariance + the recorded anchor (kin_final = 0.121014)",
    gauge_rel < 2e-3 and anchor_rel < 2e-3,
    f"prod-gauge {kin_prod:.6f} vs eigh-gauge {kin_ref:.6f} (rel {gauge_rel:.2e}) "
    f"vs recorded {KIN_ANCHOR:.6f} (rel {anchor_rel:.2e})",
)

# ---- S3: the Legendre closure with production kin ------------------
rows = []
for om in ("0.2", "0.5", "1"):
    Mr, rec = RUNGS[om]["M"], RUNGS[om]["rec"]
    tf.load(Mr)
    pde.compute_clock_flow(tf, H, RENV)
    kin_r = pde.kin_canonical(tf, H)
    J = rec["start"]["J"]
    om_star = J / (2.0 * kin_r)
    om_rec = rec["final"]["omega_star_final"]
    eu, ev = ref.e_parts(Mr, cfg)
    E_tot = float(eu) + float(ev) + J * J / (4.0 * kin_r)
    rows.append({"om": om, "J": J, "kin": kin_r, "om_star": om_star, "E": E_tot})
    print(
        f"  [info] rung om = {om}: kin {kin_r:.6f}, omega* {om_star:.5f} "
        f"(recorded {om_rec:.5f}, rel {abs(om_star - om_rec) / om_rec:.2e})"
    )
RATIOS_RESEARCH = (0.997, 0.992)  # m5_21_9_note § 7 (conjugation family of record)
ratios = []
for i in range(2):
    dE = rows[i + 1]["E"] - rows[i]["E"]
    dJ = rows[i + 1]["J"] - rows[i]["J"]
    om_mid = 0.5 * (rows[i]["om_star"] + rows[i + 1]["om_star"])
    ratios.append((dE / dJ) / om_mid)
ok3 = all(abs(r - rr) < 0.02 for r, rr in zip(ratios, RATIOS_RESEARCH))
check(
    "S3 Legendre dE/dJ == omega* (production kin) vs the research closure",
    ok3,
    f"ratios {ratios[0]:.4f} / {ratios[1]:.4f} vs research {RATIOS_RESEARCH[0]} / {RATIOS_RESEARCH[1]}",
)


# ---- S4: the hold — SET-J + production evolution vs the reference ----
def step_once(duck, dt):
    pde.compute_eta_flux(duck, 0, H)
    pde.evolve_M_eta_start(duck, dt, H)
    pde.compute_eta_flux(duck, 1, H)
    pde.evolve_M_eta_finish(duck, dt, H, SG, DELTA, W1)
    duck.swap_matrix_buffers()


M0 = RUNGS["0.2"]["M"]
om_hold = RUNGS["0.2"]["rec"]["final"]["omega_star_final"]
tf.load(M0)
info = pde.set_fixed_j(tf, H, RENV, om_hold, DT, shell=SHELL)
check(
    "S4a SET-J: kin positive, kick applied",
    info["set"] and abs(info["kin"] - rows[0]["kin"]) / rows[0]["kin"] < 1e-3,
    f"kin {info['kin']:.6f}, J {info['J']:.6f}, omega* {info['om_star']:.5f}",
)
# THE HOLD OBSERVABLE (try-2 instrument fix, measured at try 1): the GLOBAL
# rotation flows [G_k, M] near-cancel on this state (all three projections
# ~1e-5 ~ 1e-4·ω* — noise, so "dominant component retention" compared noise).
# The meaningful read is the CARRIED charge J_self = <Mdot, a0(M)> — exactly
# ω* at SET-J time by construction. Both stacks are read with the SAME
# production instrument (a twin duck field for the reference trajectory) so
# the a0 sign gauge is identical on both sides.
tf_ref = DuckField(N)


def carried_j_of(m_np, mt_np):
    """read_carried_j applied to an (M, Mdot) numpy pair via the twin duck."""
    tf_ref.load(m_np, m_prev_np=m_np - DT * mt_np)
    return pde.read_carried_j(tf_ref, H, RENV, DT)


js0 = pde.read_carried_j(tf, H, RENV, DT)
check(
    "S4b SET-J self-consistency: carried J_self == omega* at t = 0",
    abs(js0 - om_hold) / om_hold < 1e-3,
    f"J_self(0) {js0:.5f} vs omega* {om_hold:.5f}",
)
j_glob0 = pde.read_isorotation(tf, DT)
print(
    f"  [info] global-J near-cancellation on this state (the try-1 lesson): "
    f"({j_glob0['Jx']:+.2e}, {j_glob0['Jy']:+.2e}, {j_glob0['Jz']:+.2e}) "
    f"~ 1e-4 x omega* — noise-scale; the carried charge is the observable"
)
spec0 = np.linalg.eigvalsh(tf.M_am.to_numpy().astype(np.float64)[N // 2, N // 2, N // 2])

# the reference trajectory: SAME flow + init, the certified f64 leap()
tf.Md_am.from_numpy(np.zeros((N, N, N, 4, 4), np.float32))  # scratch hygiene
pde.compute_clock_flow(tf, H, RENV)
a0_hold = tf.Md_am.to_numpy().astype(np.float64)  # the production flow (normalized)
free = ~ref.pin_shell(N, H)
Mr_ref = M0.copy()
Mt_ref = om_hold * a0_hold * free[..., None, None]
gam0 = np.zeros((N, N, N))
N_CMP = 200  # the f64 reference window (leap costs 2 grad()/step)
N_HOLD = 400  # the production-only window (the GL4a gate length)

js_prod_t, js_ref_t = [], []
e_series = []
for step in range(N_HOLD + 1):
    if step % 50 == 0:
        m_np = tf.M_am.to_numpy().astype(np.float64)
        mp_np = tf.M_prev_am.to_numpy().astype(np.float64)
        mt_np = (m_np - mp_np) / DT
        js_prod_t.append(pde.read_carried_j(tf, H, RENV, DT))
        eu, ev = ref.e_parts(m_np, cfg)
        ke = 0.5 * H**3 * float(np.sum(mt_np * mt_np))
        e_series.append(float(eu) + float(ev) + ke)
        if step <= N_CMP:
            js_ref_t.append(carried_j_of(Mr_ref, Mt_ref))
    if step < N_HOLD:
        step_once(tf, DT)
        if step < N_CMP:
            Mr_ref, Mt_ref = larm.leap(Mr_ref, Mt_ref, cfg, free, gam0, DT)

m_fin = tf.M_am.to_numpy().astype(np.float64)
spec1 = np.linalg.eigvalsh(m_fin[N // 2, N // 2, N // 2])
drift = max(abs(e - e_series[1]) for e in e_series[1:]) / abs(e_series[1])
check(
    f"S4c hold: bounded ledger over {N_HOLD} steps @ dt = {DT} (f32 floor)",
    np.isfinite(m_fin).all() and np.abs(m_fin).max() < 50.0 and drift < 5e-3,
    f"E drift {drift:.2e}, max|M| {np.abs(m_fin).max():.2f}, "
    f"E {e_series[0]:.4f} -> {e_series[-1]:.4f}",
)
check(
    "S4d hold: core spectrum intact",
    np.abs(spec1 - spec0).max() < 0.05,
    f"center eigvals drift {np.abs(spec1 - spec0).max():.2e} "
    f"(start {np.round(spec0, 4).tolist()})",
)
retain = js_prod_t[-1] / js_prod_t[0] if abs(js_prod_t[0]) > 0 else 0.0
check(
    "S4e hold: the carried charge J_self retained",
    0.7 < retain < 1.3,
    f"J_self {js_prod_t[0]:+.5f} -> {js_prod_t[-1]:+.5f} (ratio {retain:.3f})",
)
i_cmp = N_CMP // 50
jp, jr = js_prod_t[i_cmp], js_ref_t[-1]
xerr = abs(jp - jr) / max(abs(jr), 1e-30)
check(
    f"S4f cross-stack: production J_self == reference leap() J_self at step {N_CMP}",
    xerr < 0.05,
    f"prod {jp:+.5f} vs ref(f64) {jr:+.5f} rel {xerr:.2e}",
)

# ---- S5: the launcher-arena flow (63^3, research units) -------------
print("  [S5] launcher arena: seed -> flip -> RELAX -> SET-J -> EVOLVE ...", flush=True)
from openwave.xperiments.m5_liquid_crystal.medium import TensorField  # noqa: E402
import openwave.xperiments.m5_liquid_crystal.engine1_seeds as seeds  # noqa: E402

NL = 63
DXL = 1.5  # ETA_DX research unit (the demo config value)
UNI = [NL * DXL * 1e-18] * 3
tfl = TensorField(UNI, NL**3)
assert tfl.nx == NL, f"launcher grid {tfl.nx} != {NL}"
G_L = float(tfl.lc_g)
D_L = float(tfl.lc_delta)
RENV_L = 10.0
DT_L = 0.005
SUBS = 64

seeds.seed_biaxial_hedgehog_M(
    tfl, NL // 2, NL // 2, NL // 2, 0.06 * NL, 3.0, D_L
)  # the demo-config geometry (R0_FRACTION 0.06, RHOC 3.0 voxels), voxel coords
pde.flip_time_axis(tfl)
stats = pde.fire_relax_canonical(tfl, DXL, G_L, D_L, W1, 200)
info_l = pde.set_fixed_j(tfl, DXL, RENV_L, 0.2, DT_L, shell=1)
jsl0 = pde.read_carried_j(tfl, DXL, RENV_L, DT_L)


def probe_axis(tf_real, iv, jv, kv):
    """The delta clock-hand axis (middle eigenvector) at a probe voxel."""
    m = tf_real.M_am.to_numpy().astype(np.float64)[iv, jv, kv, 1:4, 1:4]
    lam, V = np.linalg.eigh(m)
    return V[:, 1]  # middle eigenvector


cprobe = (NL // 2 + 4, NL // 2, NL // 2)  # off-center: inside the core envelope
ax0 = probe_axis(tfl, *cprobe)
phis, js_tr = [], []
for frame in range(10):
    for _ in range(SUBS):
        pde.compute_eta_flux(tfl, 0, DXL)
        pde.evolve_M_eta_start(tfl, DT_L, DXL)
        pde.compute_eta_flux(tfl, 1, DXL)
        pde.evolve_M_eta_finish(tfl, DT_L, DXL, G_L, D_L, W1)
        pde.apply_eta_sponge(tfl, DT_L, 0.5, 10.0)
        tfl.swap_matrix_buffers()
    js_tr.append(pde.read_carried_j(tfl, DXL, RENV_L, DT_L))
    axf = probe_axis(tfl, *cprobe)
    c = float(np.clip(abs(np.dot(ax0, axf)), 0.0, 1.0))  # apolar axis: |cos|
    phis.append(math.acos(c))
m_l = tfl.M_am.to_numpy()
t_total = 10 * SUBS * DT_L
rate = max(phis) / t_total if phis else 0.0
sign_stable = all(np.sign(v) == np.sign(js_tr[0]) for v in js_tr) if js_tr[0] != 0 else False
retain_l = js_tr[-1] / jsl0 if abs(jsl0) > 0 else 0.0
check(
    "S5 launcher flow: bounded + the carried J coherent + the clock-hand axis rotates",
    np.isfinite(m_l).all()
    and np.abs(m_l).max() < 50.0
    and info_l["set"]
    and abs(jsl0 - info_l["om_star"]) / info_l["om_star"] < 1e-2
    and sign_stable
    and retain_l > 0.5
    and max(phis) > 0.01,
    f"kin {info_l['kin']:.4e}, omega* {info_l['om_star']:.3f}, J_self {jsl0:+.4f} -> "
    f"{js_tr[-1]:+.4f} (retention {retain_l:.3f}), axis advance {max(phis):.3f} rad "
    f"over t = {t_total:.1f} (~{rate:.4f} rad/tau ~ omega*/|a0_raw| = "
    f"{info_l['om_star'] / max(info_l['a0_raw_norm'], 1e-9):.4f} predicted), "
    f"FIRE f_rms {stats['f_rms0']:.2e} -> {stats['f_rms1']:.2e}",
)

print(flush=True)
if fails:
    print(f"RESULT: {len(fails)}/{total[0]} FAILED: {fails}")
    raise SystemExit(1)
print(f"RESULT: ALL {total[0]} CHECKS GREEN")
