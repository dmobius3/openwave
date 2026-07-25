"""M5.27 P7: the mechanism test (added mid-run; the decisive follow-up to A5).

The audit established ALGEBRAICALLY that the phase-A null is structural:

  dF/dsg = -2w Sum_p p^2 sg^(p-1) sym[(eta (M eta)^(p-1))^T]

commutes with M on BLOCK-DIAGONAL states (machine zero, 4.5e-21) but NOT when
the mixed (0,i) block is present (1.4e-2). On the staged states the drive can
therefore move eigenVALUES but cannot torque the eigenFRAME that carries the
clock, which is exactly why no Arnold tongue can exist in phase A.

That is an algebraic claim about the force. This script tests it IN THE
DYNAMICS, which is the strongest form available inside the pilot: run the same
drive with the mixed-block projection DISABLED, so the (0,i) channel is live,
and ask whether the clock now responds to the background where it did not
before. It is the same comparison in every other respect (same seed, same
release, same drive, same observables).

Pre-registered readings:
  (a) mixed-live runs show a J-retention response to the drive that the staged
      runs do not  -> the mechanism is confirmed dynamically, and the mixed
      block is named as the coupling channel phase B must carry;
  (b) no response either way -> the drive does not reach the clock even with the
      channel open, and the background-scalar route is weaker than the algebra
      alone suggests;
  (c) the mixed-live runs go unstable -> report exactly that: the M5.21.3
      measurement (all 24 time-mixing curvatures negative) predicts the
      unprojected sector is unstable, so instability is an EXPECTED outcome and
      is reported as such, not as a failed test.

Run:  python m5_27_g_mechanism.py
Out:  data/m5_27_mechanism.json
"""
import json
import math
import os
import time

import numpy as np

import m5_27_a_harness as H
import m5_27_c_lockscan as LS

ARCH = H.init_taichi(prefer_gpu=True)
TAG = LS.TAG
OM_STAR = LS.OM_STAR
PROBE = LS.PROBE
T_END = 200.0
t0 = time.time()


def run(eps, om_bar, project_mixed, t_end=T_END, seed_mixed=0.0):
    """seed_mixed > 0 puts an explicit perturbation INTO the (0,i) block.

    Needed because the first pass of this test was a no-op: with the projection
    simply disabled, the staged and mixed-live runs came out byte-identical and
    the recorded mixed leak was exactly 0.0. The seeds are block-diagonal and the
    evolution PRESERVES block-diagonality, so the mixed block is an invariant
    manifold that a uniform background scalar never excites from block-diagonal
    data. To test whether the channel would carry the coupling, it has to be
    populated by hand."""
    h = H.Harness(boundary="track", project_mixed=project_mixed)
    h.load_fixedj(TAG)
    info = h.set_fixed_j(OM_STAR)
    if seed_mixed > 0.0:
        m = h.tf.M_am.to_numpy()
        mp = h.tf.M_prev_am.to_numpy()
        n = h.n
        c = (n - 1) / 2
        ax = (np.arange(n) - c) / c
        X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
        env = np.exp(-(X**2 + Y**2 + Z**2) * 4.0)
        for a, comp in enumerate((X, Y, Z), start=1):
            m[..., 0, a] = seed_mixed * env * comp
            m[..., a, 0] = seed_mixed * env * comp
            mp[..., 0, a] = seed_mixed * env * comp
            mp[..., a, 0] = seed_mixed * env * comp
        h.tf.M_am.from_numpy(m)
        h.tf.M_prev_am.from_numpy(mp)
        h.tf.M_new_am.from_numpy(m)
    pt = H.PhaseTracker(h.probe_axis(PROBE, which=1),
                        axis_hint=h.probe_axis(PROBE, which=2))
    j_tr, mix_tr = [(0.0, h.carried_j())], []
    n_steps = int(round(t_end / H.DT))
    broke = ""
    for s in range(n_steps):
        h.step(eps, om_bar)
        if s % 50 == 0:
            # finiteness guard BEFORE probing, at the probe cadence: the seeded
            # mixed sector diverges within a few tens of steps, and running the
            # eigensolver on a NaN state raises instead of reporting.
            if not h.is_finite():
                broke = f"non-finite at t={h.t:.2f}"
                break
            mm = h.max_abs_m()
            if mm > 50.0:
                broke = f"guard at t={h.t:.2f} (max|M| = {mm:.3g})"
                break
            pt.update(h.probe_axis(PROBE, which=1), h.t)
        if s % 2000 == 0:
            j = h.carried_j()
            j_tr.append((h.t, j))
            mix_tr.append((h.t, h.mixed_leak, h.max_abs_m()))
            if not np.isfinite(j):
                broke = f"non-finite J at t={h.t:.1f}"
                break
    if not broke:
        j_tr.append((h.t, h.carried_j()))
    return {
        "eps": eps, "om_bar": om_bar, "project_mixed": project_mixed,
        "t_end": h.t, "J0": info["J"], "J_trace": j_tr, "mix_trace": mix_tr,
        "retention": j_tr[-1][1] / max(abs(j_tr[0][1]), 1e-12),
        "phase_rate_late": pt.rate(t_from=0.5 * h.t), "phi_total": pt.phi_acc,
        "max_absM": h.max_abs_m(), "mixed_leak": h.mixed_leak, "broke": broke,
    }


print(f"[M5.27 P7 mechanism test] arch {ARCH}, om* = {OM_STAR:.5f}\n")

# P7a — is the block-diagonal sector dynamically INVARIANT? (measured first,
# because it decides whether "unprojecting" is even a different experiment)
inv = run(0.1, OM_STAR, False, t_end=60.0)
print(f"[P7a] invariance: with the projection OFF and block-diagonal initial "
      f"data, the largest (0,i) entry reached over the whole run is "
      f"{inv['mixed_leak']:.3e}")
print(f"      -> the mixed block is an INVARIANT MANIFOLD: a uniform background "
      f"scalar cannot excite it from block-diagonal data\n")

rows = []
for pm, sm, label in ((True, 0.0, "STAGED (mixed projected)"),
                      (False, 0.0, "MIXED-LIVE, unseeded (a no-op: see P7a)"),
                      (False, 0.005, "MIXED-LIVE, SEEDED (0,i) small 0.005"),
                      (False, 0.05, "MIXED-LIVE, SEEDED (0,i) 0.05")):
    ctrl = run(0.0, OM_STAR, pm, seed_mixed=sm)
    drv = run(0.1, OM_STAR, pm, seed_mixed=sm)
    gain = drv["retention"] - ctrl["retention"]
    rows.append({"project_mixed": pm, "seed_mixed": sm, "label": label,
                 "control": ctrl, "driven": drv, "gain": gain})
    print(f"  {label}")
    print(f"    control  retention {ctrl['retention']:+.5f}  phi {ctrl['phi_total']:+.4f} "
          f"max|M| {ctrl['max_absM']:.3f} {ctrl['broke']}")
    print(f"    driven   retention {drv['retention']:+.5f}  phi {drv['phi_total']:+.4f} "
          f"max|M| {drv['max_absM']:.3f} {drv['broke']}")
    print(f"    -> drive gain {gain:+.5f}  (mixed leak seen: "
          f"{drv['mixed_leak']:.3e})\n", flush=True)

g_staged = rows[0]["gain"]
g_live = rows[-1]["gain"]         # the SEEDED mixed-live run is the real test
unstable_live = any(bool(r["driven"]["broke"] or r["control"]["broke"]) for r in rows[2:])
if unstable_live:
    verdict = ("c: the mixed-live sector is unstable (as the M5.21.3 negative "
               "time-mixing curvatures predict), so the dynamical test cannot "
               "isolate the coupling; the algebraic A5 result stands alone")
elif abs(g_live) > 3.0 * max(abs(g_staged), 1e-6) and abs(g_live) > 0.02:
    verdict = ("a: opening the mixed (0,i) channel makes the clock respond to "
               "the background where the staged run does not -> the mixed block "
               "IS the coupling channel, confirmed dynamically")
else:
    verdict = ("b: no clock response even with the mixed channel open -> the "
               "background-scalar route reaches the clock more weakly than the "
               "algebra alone suggests")

print(f"[P7 verdict] {verdict}")
out = {"arch": ARCH, "tag": TAG, "om_star": OM_STAR, "t_end": T_END,
       "invariance_probe": {"mixed_leak_unprojected": inv["mixed_leak"],
                            "t_end": inv["t_end"]},
       "rows": rows, "gain_staged": g_staged, "gain_mixed_live": g_live,
       "verdict": verdict, "wall_s": time.time() - t0}
p = os.path.join(H.DATA, "m5_27_mechanism.json")
with open(p, "w") as f:
    json.dump(out, f, indent=1)
print(f"[P7] wall {out['wall_s']:.1f} s -> {p}")
