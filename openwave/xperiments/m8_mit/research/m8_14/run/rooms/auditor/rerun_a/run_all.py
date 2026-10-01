"""Run every script of this room in order through ./py, then the planted-defect runs, then write results.json.

Order: t1_second_variation.py, t3_spectrum.py, t4_fem.py, planted runs, collation.
A planted run 'fires' when its output contains at least one FAIL line for the check it targets.
"""
import os, sys, json, subprocess

here = os.path.dirname(os.path.abspath(__file__))
# run_all.py itself is started as `./py run_all.py`; ./py is a sandbox wrapper that cannot be nested,
# so the children are started with the very interpreter ./py launched (sys.executable).
PY = sys.executable


def run(args):
    r = subprocess.run([PY] + args, cwd=here, capture_output=True, text=True)
    return r.returncode, r.stdout + r.stderr


log = {}
for script in ("t1_second_variation.py", "t3_spectrum.py", "t4_fem.py"):
    print(f"===== {script}")
    rc, out = run([script])
    print(out, end="")
    log[script] = rc
    if rc != 0:
        print(f"!! {script} exited with code {rc}")

PLANTS = [
    ("t1_second_variation.py", "straight", None, ["C1.1", "C1.2"]),
    ("t1_second_variation.py", "seamsign", None, ["C1.3"]),
    ("t3_spectrum.py", "lam", None, ["C3.1", "C3.2"]),
    ("t3_spectrum.py", "seam", None, ["C3.3"]),
    ("t4_fem.py", "abscos", "7/5", ["C4.1"]),
    ("t4_fem.py", "ygrad", "7/5", ["C4.2"]),
    ("t4_fem.py", "asym", "7/5", ["C4.3"]),
    ("t4_fem.py", "fibidx", "7/5", ["C4.4"]),
    ("t4_fem.py", "seam", "7/5", ["C4.5"]),
]
print("===== planted defects")
planted = []
for script, plant, only, targets in PLANTS:
    args = [script, "--plant", plant] + (["--only", only] if only else [])
    rc, out = run(args)
    fired = {t: any(t in line and "FAIL" in line for line in out.splitlines()) for t in targets}
    status = "FIRED" if all(fired.values()) else "DID NOT FIRE"
    print(f"--- {script} --plant {plant}{' --only ' + only if only else ''}: targets {targets} -> {status} (rc={rc})")
    for line in out.splitlines():
        if "FAIL" in line or "Error" in line:
            print("    " + line.replace(here + os.sep, ""))
    planted.append({"script": script, "plant": plant, "only": only, "targets": targets, "fired": fired, "rc": rc})

with open(os.path.join(here, "t3_out.json")) as fh:
    t3 = json.load(fh)
with open(os.path.join(here, "t4_out.json")) as fh:
    t4 = json.load(fh)

results = {
    "R": 1,
    "t3_closed_form": t3,
    "t4": {lab: {"levels": v["levels"], "p": v["p"], "extrap": v["extrap"], "err": v["err"],
                 "ndof": v["ndof"], "ratio_of_differences": v["ratio"]} for lab, v in t4.items()},
    "extra_t4_vs_closed_form": {lab: {"extrap_minus_exact": t4[lab]["extrap"] - t3[lab]["lowest6"][0],
                                      "finest_minus_exact": [a - b for a, b in zip(t4[lab]["levels"][-1], t3[lab]["lowest6"])]}
                                for lab in t4},
    "t5": {lab: {"jacobi_lowest": t3[lab]["jacobi_lowest"], "negative_count": t3[lab]["jacobi_negative_count"],
                 "nullity": t3[lab]["jacobi_nullity"],
                 "jacobi_lowest_from_t4_extrap": t4[lab]["extrap"] - 2.0} for lab in t3},
    "planted": planted,
    "exit_codes": log,
}
with open(os.path.join(here, "results.json"), "w") as fh:
    json.dump(results, fh, indent=1)

print("===== extra: T4 extrapolation minus closed form")
for lab in t4:
    e = results["extra_t4_vs_closed_form"][lab]
    print(f"W={lab:6s} extrap-exact={e['extrap_minus_exact']:+.3e}  finest-exact=" +
          " ".join(f"{d:+.2e}" for d in e["finest_minus_exact"]))
print("all planted defects fired:", all(all(p["fired"].values()) for p in planted))
print("wrote results.json")
print("===== make_tables.py")
rc, out = run(["make_tables.py"])
print(out, end="")
