"""M5.21.15 A2 secondary blocks driver: re-runs the s = +1 (g = 32)
and s = -1 (g = 8) clock-channel blocks with the PATCHED
symmetry-breaking starts (the original in-process run predated the
even-symmetry fix, deviations log #5), merging into the existing
m5_21_15_coupled.json without touching the main block.
"""
import importlib.util
import json
import os
import time

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")

_sb = importlib.util.spec_from_file_location(
    "b15", os.path.join(HERE, "m5_21_15_b_coupled.py"))
B15 = importlib.util.module_from_spec(_sb)
_sb.loader.exec_module(B15)
C14 = B15.C14

t0 = time.time()
with open(os.path.join(DATA, "m5_21_15_coupled.json")) as f:
    out = json.load(f)

grid_opt = C14.make_grid(48, 8, 16)
grid_full = C14.make_grid(72, 12, 24)
for s, g, tag in ((+1.0, 32.0, "s+1_g32"), (-1.0, 8.0, "s-1_g8")):
    C14.S_SIGN = s
    cc = B15.ChanCorr(grid_opt, g)
    cc_full = B15.ChanCorr(grid_full, g)
    blk = {"s": s, "g": g, "kin_base": cc_full.kin_base}
    blk["scan_clock"] = B15.scan_channel(cc, cc_full, "clock",
                                         [0.01] + [0.0] * 9,
                                         f"{tag}_clock")
    blk["minkin_clock"] = B15.min_kin_channel(cc, cc_full, "clock")
    print(json.dumps({"tag": tag,
                      "minkin_clock":
                          blk["minkin_clock"]["kin_tot_min"]}),
          flush=True)
    out[tag] = blk
    with open(os.path.join(DATA, "m5_21_15_coupled.json"), "w") as f:
        json.dump(out, f, indent=1)
out["runtime_s_secondary"] = round(time.time() - t0, 1)
with open(os.path.join(DATA, "m5_21_15_coupled.json"), "w") as f:
    json.dump(out, f, indent=1)
print(json.dumps({"done": True,
                  "runtime_s": out["runtime_s_secondary"]}))
