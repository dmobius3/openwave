#!/usr/bin/env python3
"""M8.5-C2 room launch gate, implementing protocol § 12's import-closure enforcement:
every supplied module proves its closure in a PER-MODULE SUBPROCESS with a SENTINEL, so a
sys.exit inside an import cannot false-green.

Two probe classes, because the supplied python falls in two kinds:

LIBRARY (imported by in-room code): the subprocess imports the module; the import must
RETURN (SystemExit of any code is a FAIL) and the sentinel must print AFTER it returns.

SCRIPT (the five design-input evidence scripts; they execute their whole check at import,
one exits via a computed verdict): the subprocess RUNS the script from the room root with
the README's documented dependency path; pass requires exit 0 AND the script's TERMINAL
SIGNATURE, a fragment taken verbatim from each script's actual final stdout line, found
within the LAST THREE nonblank lines of the captured output. Anywhere-in-stdout matching
was a presence test that an early exit 0 after a mid-report print could false-green (the
redline finding); positional matching on the true terminal lines is what makes the claim
hold. This is stronger than import closure: the supplied evidence reproduces green, to
its own last line, in THIS room.

Self-arming (--selftest): builds two throwaway defectives under build/selftest_gate/, a
library importing a nonexistent module and a script that exits 0 EARLY without its verdict
line; the gate must FAIL both or the gate itself is broken. Run the selftest FIRST.

Usage, from the room root:
    python3 room_import_gate.py --selftest   # arms first; must end SELFTEST GREEN
    python3 room_import_gate.py              # the launch gate; must end GATE GREEN, exit 0

Preflight results land in the COMMISSIONING record, never the output ledger (§ 11).
"""
import os, subprocess, sys, shutil
os.environ["PYTHONDONTWRITEBYTECODE"] = "1"   # supplied dirs stay byte-inventory clean

SENTINEL = "M85C2-IMPORT-OK-7f3a"
ROOT = os.path.dirname(os.path.abspath(__file__))
R = "openwave/xperiments/m8_mit/research"
PILOT = os.path.join(ROOT, R, "m8_5b/pilot")           # README's documented dependency dir
DESIGN = os.path.join(R, "m8_5c/design_inputs")

LIBRARIES = [(PILOT, "route_a_nonabelian")]
SCRIPTS = [  # (path from room root, TERMINAL SIGNATURE: verbatim from each script's actual final line)
    (f"{DESIGN}/mode_count.py",              "asymptotic: ~ d_rho * N^3 / 360"),
    (f"{DESIGN}/jacobian_check.py",          "VERDICT: CONFIRMED: i*phi is an exact kernel vector"),
    (f"{DESIGN}/exact_quad_check.py",        "non-orthonormal spanning set per level."),
    (f"{DESIGN}/cascade_quad_check.py",      "VERDICT: CONFIRMED: the monitor needs its own 6N-exact rule"),
    (f"{DESIGN}/right_translation_check.py", "remain to be qualified by protocol gates 3 and 5."),
]
THIRD_PARTY = ["numpy", "scipy"]

def probe_library(pkgdir, mod):
    code = (f"import sys; sys.path.insert(0, {pkgdir!r})\n"
            f"try:\n    import {mod}\nexcept SystemExit as e:\n"
            f"    raise RuntimeError(f'SystemExit({{e.code}}) during import') from None\n"
            f"print({SENTINEL!r})")
    p = subprocess.run([sys.executable, "-c", code], cwd=ROOT,
                       capture_output=True, text=True, timeout=900)
    ok = (p.returncode == 0) and (SENTINEL in p.stdout)
    err = p.stderr.strip().splitlines()[-1] if p.stderr.strip() else ""
    return ok, f"rc={p.returncode}" + (f"  {err}" if not ok and err else "")

def probe_script(relpath, fragment):
    env = dict(os.environ, PYTHONPATH=PILOT)
    p = subprocess.run([sys.executable, relpath], cwd=ROOT, env=env,
                       capture_output=True, text=True, timeout=1800)
    tail = [l for l in p.stdout.splitlines() if l.strip()][-3:]
    has = any(fragment in l for l in tail)
    ok = (p.returncode == 0) and has
    why = "" if ok else (f"rc={p.returncode}" + ("" if has else f", terminal signature ('{fragment[:40]}') not in the last three nonblank lines"))
    err = p.stderr.strip().splitlines()[-1] if p.stderr.strip() else ""
    return ok, why + (f"  {err}" if not ok and err else "")

def run_gate():
    fails = 0
    print("room import gate (§ 12 semantics: per-module subprocess + sentinel; two probe classes)")
    print(f"interpreter: {sys.version.split()[0]} at {sys.executable}")
    for pkgdir, mod in LIBRARIES:
        ok, note = probe_library(pkgdir, mod); fails += (not ok)
        print(f"  {'PASS' if ok else 'FAIL'}  LIB    {mod:28s} {note}")
    for relpath, frag in SCRIPTS:
        ok, note = probe_script(relpath, frag); fails += (not ok)
        print(f"  {'PASS' if ok else 'FAIL'}  SCRIPT {os.path.basename(relpath):28s} {note}")
    for tp in THIRD_PARTY:
        p = subprocess.run([sys.executable, "-c",
            f"import {tp}; print({tp}.__version__); print({SENTINEL!r})"],
            capture_output=True, text=True, timeout=300)
        ok = (p.returncode == 0) and (SENTINEL in p.stdout); fails += (not ok)
        ver = p.stdout.splitlines()[0] if ok else "?"
        print(f"  {'PASS' if ok else 'FAIL'}  3RDPTY {tp:28s} version {ver}")
    print(f"GATE: {'GREEN, room launches' if fails == 0 else f'RED, {fails} failure(s); the room does NOT launch'}")
    return 0 if fails == 0 else 1

def run_selftest():
    d = os.path.join(ROOT, "build", "selftest_gate")
    shutil.rmtree(d, ignore_errors=True); os.makedirs(d)
    with open(os.path.join(d, "broken_lib.py"), "w") as f:
        f.write("import module_that_does_not_exist_m85c\n")
    with open(os.path.join(d, "early_exit_script.py"), "w") as f:
        f.write("import sys\nprint('setup fine')\nsys.exit(0)   # early, BEFORE any verdict\nprint('VERDICT: never reached')\n")
    with open(os.path.join(d, "early_token_script.py"), "w") as f:
        f.write("import sys\nprint('VERDICT: preliminary chatter, token appears EARLY')\n"
                "print('line 2 of later work')\nprint('line 3')\nprint('line 4')\n"
                "sys.exit(0)   # exits 0 with the token far above the tail\n"
                "print('VERDICT: the true terminal line, never reached')\n")
    ok1, _ = probe_library(d, "broken_lib")
    ok2, _ = probe_script(os.path.relpath(os.path.join(d, "early_exit_script.py"), ROOT), "VERDICT:")
    ok3, _ = probe_script(os.path.relpath(os.path.join(d, "early_token_script.py"), ROOT), "VERDICT:")
    print(f"  {'RED as required' if not ok1 else 'FALSE-GREEN, gate broken'}: library with unresolvable import")
    print(f"  {'RED as required' if not ok2 else 'FALSE-GREEN, gate broken'}: script exiting 0 early, no verdict line")
    print(f"  {'RED as required' if not ok3 else 'FALSE-GREEN, gate broken'}: token printed EARLY then exit 0 above the true terminal line (the presence-test defect)")
    shutil.rmtree(d, ignore_errors=True)
    green = (not ok1) and (not ok2) and (not ok3)
    print(f"SELFTEST: {'GREEN, all three arms fire' if green else 'RED, the gate cannot fail and proves nothing'}")
    return 0 if green else 1

if __name__ == "__main__":
    sys.exit(run_selftest() if "--selftest" in sys.argv else run_gate())
