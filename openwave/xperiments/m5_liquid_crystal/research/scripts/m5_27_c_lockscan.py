"""M5.27 arm B: the lock scan (protocols P0-P3).

P0  control      release the fixed-J endpoint, drive OFF (the decay baseline)
P1  capture      release + drive ON near the clock frequency
P2  tongue map   the (eps, om_bar) grid, each point classified
P3  drive-off    after any candidate lock, switch the drive off

WHAT THE BASELINE CHANGED (measured at the gates, and it sets the primary
observable): the released state does NOT hold its clock under the verified L.
Free release loses 96.6% of the carried J by t = 200 while the kinetic ledger
grows ~100x (the M5.21.3 no-free-clock result, quantified on this harness).
So the operative question is the task doc SS 3 "decisive detail": does the
prescribed background SUPPLY the J-budget that the fixed-J constraint imposes
by hand? The verdicts below are therefore read against the CONTROL, not against
an assumed persistent oscillator.

Verdicts (pre-registered, task_details TASK PLANNING, adapted to the baseline):
  SUSTAINED  J retained well above control at matched t AND mean drive power
             within the control noise band (the injection-lock signature)
  DRIVEN     J retained above control WITH sustained nonzero mean drive power
  UNSTABLE   monotone ledger growth / bounded-energy guard tripped (T2 Mathieu)
  NULL       indistinguishable from control

Usage:
  python m5_27_c_lockscan.py probe            one diagnostic point + control
  python m5_27_c_lockscan.py scan             the full (eps, om_bar) map
  python m5_27_c_lockscan.py boxcheck         the box-mode discriminator
"""
import json
import math
import os
import sys
import time

import numpy as np

import m5_27_a_harness as H

ARCH = H.init_taichi(prefer_gpu=True)

TAG = "conj_om0.2"
OM_STAR = H.FIXJ_OMSTAR[TAG]
PROBE = (15 + 4, 15, 15)
BOUNDARY = "track"          # decided by the G-vac gate
RAMP_T = H.RAMP_T          # drive-amplitude ramp duration in TIME units (B8)

# frozen registers (task_details)
EPS_GRID = [0.003, 0.01, 0.03, 0.1]
# the window covers the 1:1 AND the parametric 2:1 tongue (blindspot B10)
OMB_RATIOS = [0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 1.8, 2.0, 2.2, 2.6]

J_EVERY = 2000
PH_EVERY = 50
PW_EVERY = 200


def t_end_for(om_bar, min_t=200.0, cycles=15.0, cap=400.0):
    """Run long enough for the drive to act (>= `cycles` drive cycles) and long
    enough to sit in the window where the control decays (t ~ 200)."""
    return float(min(max(min_t, cycles * 2 * math.pi / om_bar), cap))


def run_point(eps, om_bar, t_end, n_grid=H.N_GRID, drive_off_at=None, tag=TAG):
    """One driven release. Returns the trace dict."""
    h = H.Harness(n=n_grid, boundary=BOUNDARY)
    h.load_fixedj(tag)
    info = h.set_fixed_j(H.FIXJ_OMSTAR[tag])
    pt = H.PhaseTracker(h.probe_axis(PROBE, which=1),
                        axis_hint=h.probe_axis(PROBE, which=2))
    n_steps = int(round(t_end / H.DT))
    j_tr, k_tr, p_tr, v_tr = [], [], [], []
    off_step = int(round(drive_off_at / H.DT)) if drive_off_at else None
    for s in range(n_steps):
        e_now = 0.0 if (off_step is not None and s >= off_step) else eps
        sg = h.step(e_now, om_bar, ramp_t=H.RAMP_T)
        if s % PH_EVERY == 0:
            pt.update(h.probe_axis(PROBE, which=1), h.t)
        if s % PW_EVERY == 0:
            p_tr.append((h.t, h.drive_power(e_now, om_bar, H.RAMP_T)))
        if s % J_EVERY == 0:
            j_tr.append((h.t, h.carried_j()))
            k_tr.append((h.t, h.kinetic()))
            v_tr.append((h.t, h.v4_total(sg)))
            if not np.isfinite(j_tr[-1][1]) or h.max_abs_m() > 50.0:
                break
    j_tr.append((h.t, h.carried_j()))
    k_tr.append((h.t, h.kinetic()))
    v_tr.append((h.t, h.v4_total(H.sg_of(h.t, eps, om_bar, H.RAMP_T))))
    return {
        "eps": eps, "om_bar": om_bar, "t_end": h.t, "n_grid": n_grid, "tag": tag,
        "J0": info["J"], "om_star": info["om_star"],
        "J_trace": j_tr, "kin_trace": k_tr, "power_trace": p_tr, "v4_trace": v_tr,
        "phase_hist": pt.hist[::4],
        "phase_rate_late": pt.rate(t_from=0.5 * h.t),
        "phi_total": pt.phi_acc,
        "max_absM": h.max_abs_m(), "mixed_leak": h.mixed_leak,
        "drive_off_at": drive_off_at,
    }


def summarize(r, ctrl_j_at, ctrl_kin_at):
    """Classify one run against the control at matched time, with kin and drive
    power VACUUM-REFERENCED (see `vacuum_reference`). J and phase are read raw:
    the driven vacuum carries exactly zero J and leaves the spatial block fixed."""
    t_end = r["t_end"]
    j_final = r["J_trace"][-1][1]
    j0 = r["J_trace"][0][1]
    j_ctrl = ctrl_j_at(t_end)
    ret = j_final / max(abs(j0), 1e-12)
    ret_ctrl = j_ctrl / max(abs(j0), 1e-12)
    gain = ret - ret_ctrl                      # the PRIMARY (clean) observable

    if r["eps"] > 0.0:
        vac = vacuum_reference(r["eps"], r["om_bar"], t_end, r["n_grid"])
        k_vac = float(np.interp(t_end, vac["kin"][:, 0], vac["kin"][:, 1]))
        pw_v = vac["pw"]
        pw_vac_late = pw_v[pw_v[:, 0] >= 0.5 * t_end][:, 1]
    else:
        k_vac, pw_vac_late = 0.0, np.array([0.0])
    k_defect = r["kin_trace"][-1][1] - k_vac   # vacuum-referenced kinetic
    k_ctrl = ctrl_kin_at(t_end)
    pw = np.array([p for _, p in r["power_trace"][len(r["power_trace"]) // 2:]])
    n_c = min(len(pw), len(pw_vac_late))
    pw_ref = pw[-n_c:] - pw_vac_late[-n_c:] if n_c else np.array([float("nan")])
    pw_mean = float(pw_ref.mean())
    pw_rms = float(np.sqrt((pw_ref**2).mean()))

    unstable = (not np.isfinite(j_final)) or r["max_absM"] > 50.0 or \
               (k_defect > 20.0 * max(k_ctrl, 1e-12))
    if unstable:
        verdict = "UNSTABLE"
    elif gain > 0.10:
        verdict = "SUSTAINED" if abs(pw_mean) < 0.05 * max(pw_rms, 1e-30) else "DRIVEN"
    else:
        verdict = "NULL"
    return {
        "eps": r["eps"], "om_bar": r["om_bar"], "ratio": r["om_bar"] / OM_STAR,
        "t_end": t_end, "J_final": j_final, "J_ctrl": j_ctrl,
        "retention": ret, "retention_ctrl": ret_ctrl, "gain": gain,
        "kin_raw": r["kin_trace"][-1][1], "kin_vac": k_vac,
        "kin_defect": k_defect, "kin_ctrl": k_ctrl,
        "power_mean_ref": pw_mean, "power_rms_ref": pw_rms,
        "phase_rate_late": r["phase_rate_late"], "phi_total": r["phi_total"],
        "max_absM": r["max_absM"], "verdict": verdict,
    }


_VAC_CACHE = {}


def vacuum_reference(eps, om_bar, t_end, n_grid=H.N_GRID):
    """The DEFECT-FREE box under the IDENTICAL drive (blindspot B3, extended to
    the energy ledger).

    Measured 2026-07-24 (deviation log): at eps = 0.1, om_bar = om*, a defect
    run reaches kin = 225 / P_mean = -9115 while the defect-free box under the
    same drive reaches kin = 208 / P = -9449. The whole box breathes, so raw kin
    and drive power are ~90% common-mode and CANNOT be read as defect physics.
    They are referenced against this run. The carried J needs no reference:
    measured EXACTLY 0.0 on the driven vacuum at every sample (the clock flow a0
    is off-diagonal in the spatial block while the drive moves only M00), and
    the phase probe likewise reads spatial eigenvectors the uniform drive leaves
    untouched. J and phase are therefore the clean primary observables.
    """
    key = (round(eps, 6), round(om_bar, 6), n_grid)
    if key in _VAC_CACHE and _VAC_CACHE[key]["t_end"] >= t_end - 1e-9:
        c = _VAC_CACHE[key]
    else:
        hv = H.Harness(n=n_grid, boundary=BOUNDARY)
        hv.load_vacuum()
        n_steps = int(round(t_end / H.DT))
        k_tr, p_tr = [(0.0, hv.kinetic())], [(0.0, 0.0)]
        for s in range(n_steps):
            hv.step(eps, om_bar, ramp_t=H.RAMP_T)
            if s % J_EVERY == 0:
                k_tr.append((hv.t, hv.kinetic()))
            if s % PW_EVERY == 0:
                p_tr.append((hv.t, hv.drive_power(eps, om_bar, H.RAMP_T)))
        k_tr.append((hv.t, hv.kinetic()))
        c = {"t_end": hv.t, "kin": np.array(k_tr, float), "pw": np.array(p_tr, float)}
        _VAC_CACHE[key] = c
    return c


def make_control(t_max, n_grid=H.N_GRID, tag=TAG):
    c = run_point(0.0, 1.0, t_max, n_grid=n_grid, tag=tag)
    jt = np.array(c["J_trace"], float)
    kt = np.array(c["kin_trace"], float)

    def j_at(t):
        return float(np.interp(t, jt[:, 0], jt[:, 1]))

    def k_at(t):
        return float(np.interp(t, kt[:, 0], kt[:, 1]))

    return c, j_at, k_at


# ================================================================
if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "probe"
    t0 = time.time()

    if mode == "probe":
        print(f"[P1 probe] arch {ARCH}, tag {TAG}, om* = {OM_STAR:.5f}")
        T = 200.0
        ctrl, j_at, k_at = make_control(T)
        print(f"  control: J {ctrl['J_trace'][0][1]:.5f} -> {ctrl['J_trace'][-1][1]:.5f}, "
              f"kin {ctrl['kin_trace'][0][1]:.5f} -> {ctrl['kin_trace'][-1][1]:.5f}, "
              f"phi {ctrl['phi_total']:+.4f}")
        for eps, ratio in ((0.1, 1.0), (0.1, 2.0), (0.03, 1.0)):
            r = run_point(eps, ratio * OM_STAR, T)
            s = summarize(r, j_at, k_at)
            print(f"  eps {eps:<6} om_bar/om* {ratio:<4} -> J_final {s['J_final']:+.5f} "
                  f"(ctrl {s['J_ctrl']:+.5f}, gain {s['gain']:+.4f}) "
                  f"kin_def {s['kin_defect']:.3f} (ctrl {s['kin_ctrl']:.3f}) "
                  f"P_ref {s['power_mean_ref']:+.4g} rms {s['power_rms_ref']:.4g} "
                  f"max|M| {s['max_absM']:.3f} -> {s['verdict']}", flush=True)
        print(f"[probe] wall {time.time()-t0:.1f} s")

    elif mode == "scan":
        print(f"[P2 tongue map] arch {ARCH}, tag {TAG}, om* = {OM_STAR:.5f}")
        t_max = max(t_end_for(r * OM_STAR) for r in OMB_RATIOS)
        ctrl, j_at, k_at = make_control(t_max)
        print(f"  control to t = {t_max:.0f}: J {ctrl['J_trace'][0][1]:.5f} -> "
              f"{ctrl['J_trace'][-1][1]:.5f}", flush=True)
        rows, raw = [], []
        for eps in EPS_GRID:
            for ratio in OMB_RATIOS:
                om_bar = ratio * OM_STAR
                r = run_point(eps, om_bar, t_end_for(om_bar))
                s = summarize(r, j_at, k_at)
                rows.append(s)
                raw.append({k: v for k, v in r.items() if k != "phase_hist"})
                print(f"  eps {eps:<6} ratio {ratio:<4} t {s['t_end']:.0f} -> "
                      f"gain {s['gain']:+.4f} P_ref {s['power_mean_ref']:+.3g} "
                      f"phi_rate {s['phase_rate_late']:+.5f} {s['verdict']}", flush=True)
        out = {
            "arch": ARCH, "tag": TAG, "om_star": OM_STAR, "boundary": BOUNDARY,
            "eps_grid": EPS_GRID, "omb_ratios": OMB_RATIOS, "ramp_t": RAMP_T,
            "control": {k: v for k, v in ctrl.items() if k != "phase_hist"},
            "rows": rows, "raw": raw, "wall_s": time.time() - t0,
        }
        p = os.path.join(H.DATA, "m5_27_lockscan.json")
        with open(p, "w") as f:
            json.dump(out, f, indent=1)
        from collections import Counter
        print(f"\n[scan] {Counter(x['verdict'] for x in rows)} wall "
              f"{out['wall_s']:.1f} s -> {p}")

    elif mode == "boxcheck":
        # T4 discriminator: box modes are dense (spacing pi/L ~ 0.07) and one
        # sits within 5% of om* — so any candidate structure is re-run at a
        # DIFFERENT box size, which moves box modes but not the clock.
        print(f"[T4 boxcheck] arch {ARCH}")
        T = 200.0
        res = {}
        for n_grid in (31, 39):
            ctrl, j_at, k_at = make_control(T, n_grid=n_grid)
            L = (n_grid - 1) * H.H_RES
            r = run_point(0.1, OM_STAR, T, n_grid=n_grid)
            s = summarize(r, j_at, k_at)
            res[n_grid] = {"L": L, "box1": math.pi / L, "summary": s,
                           "ctrl_J_final": ctrl["J_trace"][-1][1]}
            print(f"  n {n_grid} (L {L:.1f}, box1 {math.pi/L:.4f}): gain {s['gain']:+.4f} "
                  f"J_final {s['J_final']:+.5f} ctrl {s['J_ctrl']:+.5f} {s['verdict']}",
                  flush=True)
        p = os.path.join(H.DATA, "m5_27_boxcheck.json")
        with open(p, "w") as f:
            json.dump(res, f, indent=1, default=str)
        print(f"[boxcheck] wall {time.time()-t0:.1f} s -> {p}")
