"""M5.27: the refutation pass over the tongue map (adversarial, pre-plot).

The raw scan flagged 7 of 40 points SUSTAINED against the pre-registered
threshold `gain > 0.10`. That threshold was frozen BEFORE the control's
late-time noise amplitude was known, and it turned out to sit BELOW it. This
pass tries to refute every candidate on four independent grounds; only a
candidate that survives all four is reported as a lock.

R1 NOISE BAND     the control's own late-time |J| excursion, in retention units.
                  A gain smaller than the control's own wander is not a signal.
R2 ZERO CROSSING  by t = 400 the control's J has decayed THROUGH zero (it ends
                  at -0.021). A positive "gain" over a negative control can just
                  be two numbers that are both consistent with zero, so the
                  candidate must have |J_final| itself above the noise band.
R3 TONGUE SHAPE   an Arnold tongue WIDENS monotonically with drive amplitude.
                  Candidates at small eps with none at the largest eps are the
                  signature of noise, not of entrainment.
R4 PHASE COHERENCE a locked clock advances coherently at the drive rate (or a
                  rational fraction of it). A phase rate consistent with zero
                  means no clock is running to be locked.

Run:  python m5_27_h_refute.py
Out:  data/m5_27_verdicts.json
"""
import json
import os

import numpy as np

import m5_27_a_harness as H

SCAN = os.path.join(H.DATA, "m5_27_lockscan.json")
OUT = os.path.join(H.DATA, "m5_27_verdicts.json")

d = json.load(open(SCAN))
rows = d["rows"]
J0 = d["control"]["J_trace"][0][1]
ct = np.array(d["control"]["J_trace"], float)

# R1: the control's own late-time wander sets the noise band
late = ct[ct[:, 0] >= 0.5 * ct[-1, 0]]
band_J = float(np.abs(late[:, 1]).max())
band_ret = band_J / abs(J0)
# the undriven VISIBLE clock rate at release (gates: omega*/|a0_raw|)
VIS_RATE0 = 0.01603

print(f"[refute] control noise band: |J| <= {band_J:.5f} "
      f"= {band_ret:.4f} in retention units (J0 = {J0:.5f})")
print(f"[refute] pre-registered SUSTAINED threshold was gain > 0.10 "
      f"-> {'BELOW' if 0.10 < band_ret else 'above'} the noise band\n")

n_eps = {e: 0 for e in d["eps_grid"]}
for r in rows:
    if r["verdict"] == "SUSTAINED":
        n_eps[r["eps"]] += 1
eps_max = max(d["eps_grid"])
tongue_widens = n_eps[eps_max] >= max(n_eps.values()) and n_eps[eps_max] > 0

out_rows = []
for r in rows:
    fails = []
    if r["verdict"] == "SUSTAINED":
        if r["gain"] <= band_ret:
            fails.append(f"R1 gain {r['gain']:+.4f} <= noise band {band_ret:.4f}")
        if abs(r["J_final"]) <= band_J:
            fails.append(f"R2 |J_final| {abs(r['J_final']):.5f} <= band {band_J:.5f} "
                         f"(indistinguishable from zero)")
        if not tongue_widens:
            fails.append(f"R3 no tongue widening: SUSTAINED per eps = "
                         f"{[n_eps[e] for e in d['eps_grid']]} over eps "
                         f"{d['eps_grid']} (zero at the LARGEST drive)")
        exp_rate = r["om_bar"] / 12.4377     # drive rate in visible-probe units
        if abs(r["phase_rate_late"]) < 0.25 * VIS_RATE0:
            fails.append(f"R4 phase rate {r['phase_rate_late']:+.6f} is consistent "
                         f"with zero (undriven release rate {VIS_RATE0:.5f}, "
                         f"drive would give {exp_rate:.5f})")
    final = "NULL" if (r["verdict"] == "SUSTAINED" and fails) else r["verdict"]
    out_rows.append({**r, "verdict_raw": r["verdict"], "verdict": final,
                     "refutations": fails})

n_ref = sum(1 for x in out_rows if x["verdict_raw"] == "SUSTAINED"
            and x["verdict"] == "NULL")
survivors = [x for x in out_rows if x["verdict"] == "SUSTAINED"]
for x in out_rows:
    if x["verdict_raw"] == "SUSTAINED":
        print(f"  eps {x['eps']:<6} ratio {x['ratio']:.1f}: "
              f"{'REFUTED' if x['verdict'] == 'NULL' else 'SURVIVES'}")
        for f in x["refutations"]:
            print(f"      - {f}")

from collections import Counter
print(f"\n[refute] {n_ref} of {sum(1 for x in out_rows if x['verdict_raw']=='SUSTAINED')} "
      f"candidates refuted; {len(survivors)} survive")
print(f"[refute] corrected verdict map: {Counter(x['verdict'] for x in out_rows)}")

res = {
    "noise_band_J": band_J, "noise_band_retention": band_ret, "J0": J0,
    "sustained_per_eps": {str(k): v for k, v in n_eps.items()},
    "tongue_widens_with_eps": tongue_widens,
    "n_refuted": n_ref, "n_survivors": len(survivors),
    "rows": out_rows,
    "conclusion": (
        "NULL everywhere: no Arnold tongue at any registered (eps, om_bar). "
        "All raw SUSTAINED flags fall inside the control's own late-time J "
        "wander, do not widen with drive amplitude (zero candidates at the "
        "largest drive), and carry phase rates consistent with zero."
        if not survivors else
        f"{len(survivors)} candidate lock(s) survive refutation; see rows."
    ),
}
with open(OUT, "w") as f:
    json.dump(res, f, indent=1)
print(f"[refute] -> {OUT}")
