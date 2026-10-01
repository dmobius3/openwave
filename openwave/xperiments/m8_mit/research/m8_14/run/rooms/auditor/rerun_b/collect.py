"""Merge exact.json (T3/T5) and t4.json (T4) into results.json; print the T4 tables."""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(HERE, "exact.json")) as fh:
    exact = json.load(fh)
with open(os.path.join(HERE, "t4.json")) as fh:
    t4 = json.load(fh)

out = {"R": 1, "t4": t4["t4"], "t4_extra_all_six": t4["t4_extra_all_six"],
       "t3_t5_exact": exact, "t4_minus_exact_extra": {}}
for lab, d in t4["t4"].items():
    ex = exact[lab]["eigs"]
    out["t4_minus_exact_extra"][lab] = {
        "finest_minus_exact": [f - e["value"] for f, e in zip(d["levels"][-1], ex)],
        "extrap_minus_exact_bottom": d["extrap"] - ex[0]["value"]}
    print(f"\n### W = {lab}\n")
    print("| level | " + " | ".join(f"λ{i+1}" for i in range(6)) + " |")
    print("|---" * 7 + "|")
    for i, row in enumerate(d["levels"]):
        print(f"| {i} | " + " | ".join(f"{v:.10f}" for v in row) + " |")
    print(f"| exact | " + " | ".join(f"{e['value']:.10f}" for e in ex) + " |")
    print(f"\np = {d['p']:.6f}, extrap = {d['extrap']:.12f}, err = {d['err']:.3e}, "
          f"extrap - exact = {out['t4_minus_exact_extra'][lab]['extrap_minus_exact_bottom']:.3e}")
with open(os.path.join(HERE, "results.json"), "w") as fh:
    json.dump(out, fh, indent=1)
