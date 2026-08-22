#!/usr/bin/env python3
"""M8.5-B Phase A: the single qualification entry point.

    python3 research/m8_5b/run_qualification.py

Runs, in order, every check that backs the Phase A claim, and ends with one
verdict.  Nonzero exit on anything short of a full pass.

    0  provenance   subprocess path confinement (PYTHONPATH cleared,
                    PYTHONNOUSERSITE set, targets under the tree,
                    inherited) plus inspected first-party origins in
                    the controlling process
    1  environment  interpreter and library versions
    2  manifest     every shipped file matches its recorded hash
    3  schema       Packet I and Packet II gate batteries
    4  structural   target-scored battery, 8 of 8, plus non-vacuity
    4b evaluator    the rung-3b evaluator's own validation and mutations
    5  rehearsal    Q4 integration rehearsal, including the deletion limb
    6  integrated   the Q1/Q2/Q3/Q5 battery
    7  records      the fresh run reproduces the shipped qualification records

WHAT THIS DOES NOT SHOW.  No qualification run in this tree is adjudication
evidence; every case exercised by this command is synthetic or a frozen
tuning case.
"""

import hashlib
import json
import os
import platform
import subprocess
import sys

ROOT = os.path.dirname(os.path.abspath(__file__))
QUAL = os.path.join(ROOT, "qualification")
ok = True
results = []


def step(n, title):
    print(f"\n{'=' * 70}\n{n}. {title}\n{'=' * 70}")


def record(name, passed, detail=""):
    global ok
    ok &= bool(passed)
    results.append((name, bool(passed), detail))
    print(f"   {'PASS' if passed else 'FAIL'}  {name}" + (f"   {detail}" if detail else ""))


def run(script, *args):
    """Run a battery under subprocess path confinement.

    `PYTHONPATH` is cleared and `PYTHONNOUSERSITE` is set, and the target is
    required to live under this tree.  The environment is inherited, so the
    confinement also covers the batteries' own subprocesses, which a
    parent-side `sys.modules` scan cannot reach.

    SCOPE, stated exactly.  This establishes that no qualification subprocess
    inherits an external first-party path through `PYTHONPATH` or the user
    site directory.  It is not a claim of absolute filesystem confinement:
    system site-packages and any `.pth` configuration remain visible, as they
    must be for NumPy and SciPy.  The complementary check is step 0, which
    affirmatively inspects the origins of the modules the controlling process
    imports.
    """
    target = os.path.join(ROOT, script)
    if not os.path.abspath(target).startswith(ROOT + os.sep):
        raise RuntimeError(f"refusing to run a script outside the tree: {script}")
    env = {**os.environ, "PYTHONPATH": "", "PYTHONNOUSERSITE": "1"}
    r = subprocess.run([sys.executable, target, *args],
                       cwd=ROOT, capture_output=True, text=True, env=env)
    return r.returncode, r.stdout + r.stderr


# --- 0 provenance ------------------------------------------------------------
step(0, "provenance: imports confined to this tree")
sys.path.insert(0, ROOT)
for sub in ("gates", "production", "pilot", "eval3b"):
    sys.path.insert(0, os.path.join(ROOT, sub))
before = set(sys.modules)
import packet_schema                      # noqa: E402
sys.path_importer_cache.clear()
import adjudication_gates, adapter_3b, step3_schema, lauret_evaluator  # noqa: E402,F401
outside = []
for name, mod in sys.modules.items():
    if name in before:
        continue
    f = getattr(mod, "__file__", None) or ""
    if f.endswith(".py") and not f.startswith(ROOT) and "site-packages" not in f \
       and "lib/python" not in f:
        outside.append((name, f))
record("no first-party module imported from outside the tree "
       "(parent closure, inspected)", not outside, str(outside[:3]))
# Probed, not asserted.  This line used to pass a literal True, so the gate
# could not detect its own regression: deleting `env=env` from `run()` left it
# printing PASS.  The probe goes THROUGH `run()`, so it observes the environment
# the batteries themselves receive, and it runs under a deliberately poisoned
# parent, because a probe run under a clean environment cannot fail: with no
# hostile PYTHONPATH to inherit, dropping the confinement would be
# indistinguishable from enforcing it.
POISON = "/nonexistent/m8_5b_confinement_probe_must_not_be_inherited"
_saved = {k: os.environ.get(k) for k in ("PYTHONPATH", "PYTHONNOUSERSITE")}
os.environ["PYTHONPATH"] = POISON
os.environ["PYTHONNOUSERSITE"] = ""
try:
    _rc, _out = run("confinement_probe.py")
finally:
    for _k, _v in _saved.items():
        if _v is None:
            os.environ.pop(_k, None)
        else:
            os.environ[_k] = _v

_report = [l for l in _out.splitlines() if l.startswith("CONFINEMENT_PROBE_JSON ")]
if _rc == 0 and _report:
    _seen = json.loads(_report[0].split(" ", 1)[1])
    _bad = []
    for _who, _v in _seen.items():
        if _v["PYTHONPATH"] != "":
            _bad.append(f"{_who} inherited PYTHONPATH={_v['PYTHONPATH']!r}")
        if _v["PYTHONNOUSERSITE"] != "1":
            _bad.append(f"{_who} PYTHONNOUSERSITE={_v['PYTHONNOUSERSITE']!r}")
        if any(POISON in _p for _p in _v["sys_path"]):
            _bad.append(f"{_who} carries the poisoned path on sys.path")
    record("subprocess path confinement observed in a child AND a grandchild, "
           "probed under a poisoned parent environment", not _bad, str(_bad[:3]))
else:
    record("subprocess path confinement probe reported", False,
           f"rc={_rc}, no probe report on stdout")

# --- 1 environment -----------------------------------------------------------
step(1, "environment")
import numpy, scipy                       # noqa: E402
env = {"python": sys.version.split()[0], "numpy": numpy.__version__,
       "scipy": scipy.__version__, "platform": platform.platform()}
for k, v in env.items():
    print(f"   {k:10} {v}")
shipped = os.path.join(QUAL, "ENVIRONMENT.json")
record("ENVIRONMENT.json present", os.path.exists(shipped))
if os.path.exists(shipped):
    was = json.load(open(shipped))
    same = {k: was.get(k) for k in ("python", "numpy", "scipy")} == \
           {k: env[k] for k in ("python", "numpy", "scipy")}
    # Informational, never a verdict: qualification has reproduced off the
    # recorded stack, so an exact version match is not required for a pass.
    print(f"   INFO  recorded environment {'matches' if same else 'DIFFERS from'} "
          f"the current one; qualification comparison continues regardless")
else:
    print("   INFO  no recorded environment to compare against")

# --- 2 manifest --------------------------------------------------------------
step(2, "manifest: shipped bytes match their recorded hashes")
man_path = os.path.join(QUAL, "MANIFEST.json")
if os.path.exists(man_path):
    man = json.load(open(man_path))
    bad = []
    for rel, h in sorted(man["files"].items()):
        p = os.path.join(ROOT, rel)
        got = hashlib.sha256(open(p, "rb").read()).hexdigest() if os.path.exists(p) else None
        if got != h:
            bad.append(rel)
    record(f"{len(man['files'])} files match MANIFEST.json", not bad, str(bad[:4]))
    # the other half of the chain: the freeze record must pin THIS manifest
    frz = os.path.join(QUAL, "PHASE_A_FREEZE.md")
    if os.path.exists(frz):
        import re as _re
        pinned = _re.search(r"SHA-256 ([0-9a-f]{64})", open(frz).read())
        live = hashlib.sha256(open(man_path, "rb").read()).hexdigest()
        record("PHASE_A_FREEZE.md pins this manifest",
               bool(pinned) and pinned.group(1) == live,
               f"freeze {pinned.group(1)[:12] if pinned else 'absent'}, manifest {live[:12]}")
    else:
        record("PHASE_A_FREEZE.md present", False, "missing")
else:
    record("MANIFEST.json present", False, "missing")

# --- 3 schema ----------------------------------------------------------------
step(3, "Packet I and Packet II gate batteries")
rc, out = run("packet_schema.py")
record("packet_schema.py suite exits 0", rc == 0)
record("suite demonstrates it can fail", "suite behaves: the gate can fail" in out)

# --- 4 structural ------------------------------------------------------------
step(4, "structural coverage and non-vacuity")
rc, out = run("verify_packet_structural.py")
record("verify_packet_structural.py exits 0", rc == 0)
record("coverage 8 of 8 predicates target-exercised", "(8/8)" in out)

# --- 4b eval3b ---------------------------------------------------------------
step("4b", "rung-3b evaluator: Gamma = 1 validation and mutation battery")
rc, out = run("eval3b/gate_gamma1.py")
record("gate_gamma1.py exits 0 (recovers the rung-2 tower)", rc == 0)
rc, out = run("eval3b/mutations.py")
record("eval3b mutation battery exits 0", rc == 0)

# --- 4c route (a) group-closure battery (Addendum 12.3) ----------------------
step("4c", "route (a) group-closure battery: central equivalence, both regimes")
rc, out = run("route_a_closure_battery.py")
record("route (a) closure battery exits 0", rc == 0)
# The prefix is load-bearing: `ck` prints its label under PASS or FAIL
# alike, so matching the bare label would be satisfied by the arm failing.
record("closure battery demonstrates the repair is load-bearing",
       "PASS  mutation: removing central equivalence makes L(7;1,2) fail" in out
       and "BATTERY PASS" in out)

# --- 5 rehearsal -------------------------------------------------------------
step(5, "Q4 integration rehearsal, including the deletion limb")
rc, out = run("rehearsal_q4.py")
record("rehearsal_q4.py exits 0", rc == 0)
record("3a and 3b GREEN on both routes",
       out.count("GREEN route a") >= 2 and out.count("GREEN route b") >= 2)
record("deletion limb completed from committed artifacts alone",
       "adjudication completed from committed artifacts alone" in out)

# --- 6 integrated ------------------------------------------------------------
step(6, "integrated qualification battery")
rc, out = run("qualify_integration.py")
record("qualify_integration.py exits 0", rc == 0)
fresh = os.path.join(ROOT, "rehearsal", "QUALIFY_RECORD.json")
n_items = n_fail = None
record("the fresh battery wrote its record", os.path.exists(fresh))
if os.path.exists(fresh):
    q = json.load(open(fresh))
    n_items, n_fail = len(q["results"]), q["failed"]
    record(f"{n_items} items, {n_fail} failed", n_fail == 0)

# --- 7 records ---------------------------------------------------------------
step(7, "the fresh run reproduces the shipped qualification records")
ship = os.path.join(QUAL, "QUALIFY_RECORD.json")
record("the shipped qualification record is present",
       os.path.exists(ship) and n_items is not None)
if os.path.exists(ship) and n_items is not None:
    s = json.load(open(ship))
    fresh_q = json.load(open(fresh))
    record("fresh results array equals the shipped one element-for-element",
           s["results"] == fresh_q["results"],
           f"shipped {len(s['results'])} items, fresh {n_items}")
    record("shipped record also reports zero failures", s["failed"] == 0)

# --- verdict -----------------------------------------------------------------
print(f"\n{'=' * 70}")
if ok:
    print(f"PHASE A QUALIFICATION: PASS - structural 8/8; integrated "
          f"{n_items}/{n_items}; Q4 GREEN; deletion GREEN")
    print("No qualification run in this tree is adjudication evidence; every case")
    print("exercised by this command is synthetic or a frozen tuning case.")
else:
    print("PHASE A QUALIFICATION: FAIL")
    for name, passed, detail in results:
        if not passed:
            print(f"   FAILED: {name}   {detail}")
print("=" * 70)
sys.exit(0 if ok else 1)
