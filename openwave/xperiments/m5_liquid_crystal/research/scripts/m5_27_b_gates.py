"""M5.27 arm A gates: the pre-scan battery (DoD item 2).

Every gate here must be GREEN before any (kappa, om_bar) scan point is trusted.
The battery implements the tripwires and blindspots frozen in the task plan:

  G-vac    defect-free driven box tracks the analytic adiabatic vacuum with no
           spurious interior gradients — AND DECIDES the boundary handling (B7)
  G-phase  the clock-phase instrument reproduces the known omega* on the
           fixed-J live hold, with the apolar mod-pi fold (B9)
  G-power  the drive-power ledger dV4/dsg matches a finite-difference of V4
           (the analytic kernel is exact, so this is a machine check)
  G-box    the box mode spectrum + the stiff M00 mode are pre-computed, so an
           om_bar scan can mark the frequencies where structure is spurious (T4)
  G-dt     the dt cost of the high-om_bar (Kapitza) window is measured (T5)
  G-reg    kappa = 0 regression: the live hold reproduces the M5.23.2 anchor
           and the free release reproduces the M5.21.3 decay baseline
  G-static statics at kappa != 0 with the drive ON (T3)

Run:  python m5_27_b_gates.py
Out:  data/m5_27_gates.json
"""
import json
import math
import os
import time

import numpy as np

import m5_27_a_harness as H

ARCH = H.init_taichi(prefer_gpu=True)
OUT = os.path.join(H.DATA, "m5_27_gates.json")
RESULTS = []


def check(name, ok, detail=""):
    RESULTS.append({"gate": name, "pass": bool(ok), "detail": detail})
    print(f"{'PASS' if ok else 'FAIL'}  {name}\n      {detail}", flush=True)
    return ok


PROBE = (15 + 4, 15, 15)          # the M5.23.2 D1-proven probe offset (r = 4 vox)
TAG = "conj_om0.2"                # the primary endpoint
OM_STAR = H.FIXJ_OMSTAR[TAG]

t_start = time.time()
print(f"[M5.27 gates] arch = {ARCH}\n")

# ================================================================
# G-box — the spectra an om_bar scan must not confuse with physics (T4)
# ================================================================
print("[G-box] box + stiff-mode spectrum")
# The stiff M00 mode: uniform field, V4 curvature about the vacuum in the (0,0)
# direction. V4 = w Sum_p (t_p - C_p)^2 with t_p from M = diag(-x,1,delta,0)
# (eta = diag(-1,1,1,1) so (M eta)_00 = +x): t_p = x^p + 1 + delta^p, so
# V4(x) = w Sum_p (x^p - g^p)^2. The mode frequency is sqrt(V4''(x=g)) since the
# kinetic term is 1/2 ||Mdot||^2 with unit mass per component.
g, dl, w = H.G0, H.DELTA, H.W1


def v4_uniform(x):
    return w * sum((x**p - g**p) ** 2 for p in (1, 2, 3, 4))


hstep = 1e-3
v4pp = (v4_uniform(g + hstep) - 2 * v4_uniform(g) + v4_uniform(g - hstep)) / hstep**2
om_m00 = math.sqrt(max(v4pp, 0.0))
# box modes of the finite grid: k_n = n pi / L with L = (n_grid - 1) h, and the
# curvature sector propagates at the canonical speed 1 in these units
L = (H.N_GRID - 1) * H.H_RES
box_modes = [n * math.pi / L for n in range(1, 6)]
check(
    "G-box: stiff M00 mode + box modes pre-computed",
    om_m00 > 10.0 and len(box_modes) == 5,
    f"omega_M00 = {om_m00:.2f} (V4'' = {v4pp:.4g}); box modes {['%.3f' % b for b in box_modes]}; "
    f"clock band omega* = {OM_STAR:.4f} (ratio to M00 = {OM_STAR/om_m00:.2e})",
)

# ================================================================
# G-vac — the defect-free driven box (B7: boundary handling DECIDED here)
# ================================================================
print("\n[G-vac] defect-free driven box: analytic tracking + boundary handling")
EPS_V, OMB_V = 0.03, 1.0
N_V = int(round(6 * 2 * math.pi / OMB_V / H.DT))   # 6 drive cycles
vac_rows = {}
for mode in ("track", "pin"):
    hv = H.Harness(boundary=mode, project_mixed=True)
    hv.load_vacuum()
    dev, grad_max = [], 0.0
    for s in range(N_V):
        sg = hv.step(EPS_V, OMB_V, ramp_t=20.0)
        if s % 200 == 0:
            m = hv.tf.M_am.to_numpy()
            c = H.N_GRID // 2
            m00_c = float(m[c, c, c, 0, 0])
            dev.append(abs(m00_c - H.adiabatic_vacuum_m00(sg)) / abs(sg))
            inner = m[2:-2, 2:-2, 2:-2, 0, 0]
            grad_max = max(grad_max, float(np.abs(np.diff(inner, axis=0)).max()))
    vac_rows[mode] = {
        "max_rel_dev_from_adiabatic": float(max(dev)),
        "max_interior_gradient": grad_max,
    }
    print(f"      {mode}: rel dev {max(dev):.3e}, max interior |d M00| {grad_max:.3e}")

best = min(vac_rows, key=lambda k: vac_rows[k]["max_interior_gradient"])
check(
    "G-vac: driven vacuum tracks the analytic adiabatic response, no spurious interior gradients",
    vac_rows[best]["max_rel_dev_from_adiabatic"] < 0.05
    and vac_rows[best]["max_interior_gradient"] < 1e-3,
    f"BOUNDARY DECISION = '{best}'. track: dev {vac_rows['track']['max_rel_dev_from_adiabatic']:.3e} "
    f"grad {vac_rows['track']['max_interior_gradient']:.3e} | "
    f"pin: dev {vac_rows['pin']['max_rel_dev_from_adiabatic']:.3e} "
    f"grad {vac_rows['pin']['max_interior_gradient']:.3e}",
)
BOUNDARY = best

# ================================================================
# G-power — the drive-power ledger vs a finite difference of V4
# ================================================================
print("\n[G-power] dV4/dsg analytic kernel vs an independent f64 numpy reference")
hp = H.Harness(boundary=BOUNDARY)
hp.load_fixedj(TAG)


def v4_and_dv4_numpy(M, sg, delta, w1):
    """Independent f64 reference for BOTH V4 and dV4/dsg (also the seed of the
    audit arm). eta = diag(-1,1,1,1); t_p = tr((M eta)^p)."""
    eta = np.diag([-1.0, 1.0, 1.0, 1.0])
    me = M @ eta
    p2 = me @ me
    p3 = p2 @ me
    p4 = p3 @ me
    tr = [np.trace(x, axis1=-2, axis2=-1) for x in (me, p2, p3, p4)]
    cp = [sg**p + 1.0 + delta**p for p in (1, 2, 3, 4)]
    v4 = w1 * sum((tr[p - 1] - cp[p - 1]) ** 2 for p in (1, 2, 3, 4))
    dv4 = -2.0 * w1 * sum(
        p * sg ** (p - 1) * (tr[p - 1] - cp[p - 1]) for p in (1, 2, 3, 4)
    )
    return float(v4.sum()), float(dv4.sum())


M_np = hp.tf.M_am.to_numpy().astype(np.float64)
v4_ref, dv4_ref = v4_and_dv4_numpy(M_np, g, dl, w)
# Self-check of the reference by CENTRAL DIFFERENCE in f64, at two step sizes:
# a correct analytic derivative leaves only O(h^2) truncation, so halving h must
# cut the residual ~4x. (A wrong analytic form would leave an h-independent
# floor.) A direct f32 FD on the grid sums is unusable: catastrophic cancellation.
def fd_at(hfd):
    return (v4_and_dv4_numpy(M_np, g + hfd, dl, w)[0]
            - v4_and_dv4_numpy(M_np, g - hfd, dl, w)[0]) / (2 * hfd)


err_h = abs(dv4_ref - fd_at(1e-3))
err_h2 = abs(dv4_ref - fd_at(5e-4))
conv_ratio = err_h / max(err_h2, 1e-30)
d_an = hp.dv4_dsg(g)
v4_an = hp.v4_total(g)
rel_p = abs(d_an - dv4_ref) / max(abs(dv4_ref), 1e-12)
rel_v = abs(v4_an - v4_ref) / max(abs(v4_ref), 1e-12)
check(
    "G-power: the taichi ledger matches an independent f64 numpy reference (V4 and dV4/dsg)",
    rel_p < 5e-4 and rel_v < 5e-4 and 3.0 < conv_ratio < 5.0,
    f"dV4/dsg taichi {d_an:.4f} vs numpy-f64 {dv4_ref:.4f} (rel {rel_p:.2e}, f32 grid-sum level); "
    f"V4 taichi {v4_an:.5f} vs numpy {v4_ref:.5f} (rel {rel_v:.2e}); "
    f"numpy FD truncation falls {conv_ratio:.2f}x on halving h (O(h^2) expected 4x) "
    f"-> the analytic dV4/dsg is confirmed exact",
)

# ================================================================
# G-phase — the clock-phase instrument on the live hold (B9)
# ================================================================
print("\n[G-phase] phase instrument vs the known omega* on the fixed-J live hold")
hph = H.Harness(boundary=BOUNDARY)
hph.load_fixedj(TAG)
info_j = hph.set_fixed_j(OM_STAR)
rate_pred = info_j["om_star"] / max(info_j["a0_raw_norm"], 1e-9)   # the D1 visible rate
# the rotation AXIS is the local LEADING spatial eigenvector (a0_conj rotates
# about it); the clock hand is the MIDDLE eigenvector rotating in its plane.
axis_hint = hph.probe_axis(PROBE, which=2)
pt = H.PhaseTracker(hph.probe_axis(PROBE, which=1), axis_hint=axis_hint)
# Validate in the SAME regime the delivered D1 gate certified (t ~ 1.6): the
# released state decays (the M5.21.3 no-free-clock result, quantified by the
# free baseline below), so a long window measures decay, not the clock rate.
N_PH = 700
for s in range(N_PH):
    hph.step(0.0, 0.0)
    if s % 5 == 0:
        pt.update(hph.probe_axis(PROBE, which=1), hph.t)
rate_meas = pt.rate()
ratio = rate_meas / rate_pred if rate_pred else float("nan")
check(
    "G-phase: unwrapped apolar phase advances at the visible clock rate",
    abs(pt.phi_acc) > 0.01 and 0.3 < abs(ratio) < 3.0,
    f"phi advance {pt.phi_acc:+.5f} rad over t = {hph.t:.2f}; rate {rate_meas:+.5f} "
    f"vs predicted omega*/|a0_raw| = {rate_pred:.5f} (ratio {ratio:+.3f}); "
    f"rotation axis = local leading spatial eigenvector",
)

# ================================================================
# G-reg — kappa = 0 regression: live hold + free release (B4 noise floor)
# ================================================================
print("\n[G-reg] kappa = 0 regression: the live hold anchor + the free-release baseline")
hr = H.Harness(boundary=BOUNDARY)
hr.load_fixedj(TAG)
info_r = hr.set_fixed_j(OM_STAR)
j0 = hr.carried_j()
for _ in range(100):
    hr.step(0.0, 0.0)
j100 = hr.carried_j()
drift_100 = (j0 - j100) / max(abs(j0), 1e-12)
check(
    "G-reg: the fixed-J live hold reproduces the delivered anchor over 100 steps",
    abs(drift_100) < 0.02,
    f"J {j0:.5f} -> {j100:.5f} ({drift_100*100:+.3f}% over 100 steps); "
    f"M5.23.2 anchor 0.19923 -> 0.19865 (-0.29%)",
)

# the undriven long-run baseline = the P0 control AND the B4 phase-noise floor
hb = H.Harness(boundary=BOUNDARY)
hb.load_fixedj(TAG)
info_b = hb.set_fixed_j(OM_STAR)
rate_pred_b = info_b["om_star"] / max(info_b["a0_raw_norm"], 1e-9)
ptb = H.PhaseTracker(hb.probe_axis(PROBE, which=1),
                     axis_hint=hb.probe_axis(PROBE, which=2))
jb, kb = [], []
N_B = 40000
for s in range(N_B):
    hb.step(0.0, 0.0)
    if s % 50 == 0:
        ptb.update(hb.probe_axis(PROBE), hb.t)
    if s % 4000 == 0:
        jb.append((hb.t, hb.carried_j()))
        kb.append((hb.t, hb.kinetic()))
jb.append((hb.t, hb.carried_j()))
kb.append((hb.t, hb.kinetic()))
rate_free = ptb.rate()
j_decay = (jb[0][1] - jb[-1][1]) / max(abs(jb[0][1]), 1e-12)
check(
    "G-reg: the free-release baseline is measured (the P0 control + noise floor)",
    np.isfinite(rate_free) and hb.max_abs_m() < 50.0,
    f"t = {hb.t:.1f} ({N_B} steps): J {jb[0][1]:.5f} -> {jb[-1][1]:.5f} ({j_decay*100:+.2f}%), "
    f"kin {kb[0][1]:.5f} -> {kb[-1][1]:.5f}, free phase rate {rate_free:+.5f} "
    f"(pred {rate_pred_b:.5f}), max|M| {hb.max_abs_m():.3f}, leak {hb.mixed_leak:.2e}",
)

# ================================================================
# G-static — statics at kappa != 0 with the drive ON (T3)
# ================================================================
print("\n[G-static] the inherited static sector under the drive (T3)")
hs = H.Harness(boundary=BOUNDARY)
hs.load_fixedj(TAG)
m_before = hs.tf.M_am.to_numpy().astype(np.float64)
c = H.N_GRID // 2
spec_before = np.linalg.eigvalsh(m_before[c, c, c, 1:4, 1:4])
EPS_S, OMB_S = 0.03, OM_STAR
N_S = int(round(6 * 2 * math.pi / OMB_S / H.DT))
for _ in range(N_S):
    hs.step(EPS_S, OMB_S)
m_after = hs.tf.M_am.to_numpy().astype(np.float64)
spec_after = np.linalg.eigvalsh(m_after[c, c, c, 1:4, 1:4])
spec_shift = float(np.abs(spec_after - spec_before).max())
check(
    "G-static: the spatial core spectrum survives the drive (topology-protected sector intact)",
    np.isfinite(spec_shift) and spec_shift < 0.5 and hs.max_abs_m() < 50.0,
    f"core spatial eigenvalues {np.round(spec_before,4).tolist()} -> "
    f"{np.round(spec_after,4).tolist()} (max shift {spec_shift:.4f}) at eps = {EPS_S}",
)

# ================================================================
# G-dt — the cost/stability of the high-om_bar (Kapitza) window (T5)
# ================================================================
print("\n[G-dt] the high-om_bar window numerics pre-check")
OMB_K = 10.0
steps_per_cycle = (2 * math.pi / OMB_K) / H.DT
hk = H.Harness(boundary=BOUNDARY)
hk.load_fixedj(TAG)
hk.set_fixed_j(OM_STAR)
for _ in range(4000):
    hk.step(0.03, OMB_K)
check(
    "G-dt: the Kapitza window is resolvable at the certified dt",
    steps_per_cycle > 20 and hk.max_abs_m() < 50.0 and np.isfinite(hk.max_abs_m()),
    f"om_bar = {OMB_K}: {steps_per_cycle:.0f} steps/cycle at dt = {H.DT}; "
    f"4000 steps stable, max|M| {hk.max_abs_m():.3f}",
)

# ================================================================
out = {
    "arch": ARCH,
    "tag": TAG,
    "om_star": OM_STAR,
    "probe": list(PROBE),
    "boundary_decision": BOUNDARY,
    "omega_M00": om_m00,
    "box_modes": box_modes,
    "vacuum_gate": vac_rows,
    "free_baseline": {
        "steps": N_B,
        "t_end": hb.t,
        "J_trace": jb,
        "kin_trace": kb,
        "phase_rate_free": rate_free,
        "phase_rate_pred": rate_pred_b,
    },
    "phase_gate": {
        "phi_acc": pt.phi_acc,
        "rate_meas": rate_meas,
        "rate_pred": rate_pred,
        "ratio": ratio,
    },
    "results": RESULTS,
    "all_pass": all(r["pass"] for r in RESULTS),
    "wall_s": time.time() - t_start,
}
with open(OUT, "w") as f:
    json.dump(out, f, indent=1)
n_pass = sum(r["pass"] for r in RESULTS)
print(f"\n[M5.27 gates] {n_pass}/{len(RESULTS)} pass, boundary = {BOUNDARY}, "
      f"wall {out['wall_s']:.1f} s -> {OUT}")
