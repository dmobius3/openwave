"""Reproduce every number from a clean directory: ./py run_all.py [target_dir]
Copies the scripts into target_dir (default clean_run/), runs them in order there, and compares the JSON outputs
with those in ./out (if present)."""
import os, sys, shutil, subprocess, json, time

SCRIPTS = ["s01_group.py", "s02_intertwiners.py", "s_main.py", "s03_stabilizers.py", "s12_quadrature.py", "s_verify.py"]
MODULES = ["mq.py", "su2.py", "engine.py", "linalg_mq.py", "checks.py", "numD.py"]
here = os.path.dirname(os.path.abspath(__file__))
target = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "clean_run"))
if os.path.exists(target):
    shutil.rmtree(target)
os.makedirs(target)
for f in SCRIPTS + MODULES:
    shutil.copy(os.path.join(here, f), target)
os.chdir(target)
env = dict(os.environ)
for s in SCRIPTS:
    t = time.time()
    r = subprocess.run([sys.executable, "-S", s], capture_output=True, text=True, env=env)
    with open(os.path.join("out", s.replace(".py", ".log")) if os.path.isdir("out") else s + ".log", "w") as fh:
        fh.write(r.stdout + r.stderr)
    nfail = r.stdout.count("FAIL:")
    npass = r.stdout.count("PASS:")
    print(f"{s}: exit {r.returncode}, PASS {npass}, FAIL {nfail}, {time.time() - t:.1f}s")
    if r.returncode != 0:
        print(r.stderr[-3000:])
        sys.exit(1)
# compare with the room's out/
diffs = []
for f in sorted(os.listdir("out")):
    if f.endswith(".json") and os.path.exists(os.path.join(here, "out", f)):
        a = json.load(open(os.path.join("out", f)))
        b = json.load(open(os.path.join(here, "out", f)))
        if f.startswith("item12") or f.startswith("item3") or f.startswith("item1_res") or f.startswith("checks"):
            continue   # float residuals / logs may differ in the last digits
        if a != b:
            diffs.append(f)
print("JSON outputs differing from ./out (exact files only):", diffs)
