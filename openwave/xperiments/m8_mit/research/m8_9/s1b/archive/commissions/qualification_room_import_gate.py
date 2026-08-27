#!/usr/bin/env python3
"""Launch precondition: the room's import closure, its origins, and its manifest.

Copying a directory is not supplying its dependencies. This room shipped
incomplete three times in this program, and the gate that was supposed to catch
that was itself false-green once: modules in this tree run a self-test battery
at import and call sys.exit(0), which terminated the checker with status 0 and
no output. Exit code 0 is therefore NOT the success condition here.

    PYTHONPATH=.:m8_5b:m8_5b/pilot:m8_5b/production:m8_5b/gates \
        python3 room_import_gate.py

Three things must hold, and each can fail independently:

  1. SENTINEL.  Every REQUIRED module must print IMPORT_COMPLETE:<name> AFTER
     importlib.import_module() returns. Exit 0 without the sentinel is RED.
  2. ORIGIN.    Every REQUIRED module's resolved __file__ must live under this
     room, so a missing dependency cannot silently satisfy itself from the live
     checkout or another PYTHONPATH entry.
  3. MANIFEST.  Every supplied file is hashed. The manifest is written once and
     re-verified afterwards; the commissioned unit checks it before writing code
     and again at the end.

REQUIRED means: in the transitive import closure of the S1b entry points. A
module outside that closure that exits at import is a SCRIPT and is classified,
not failed. One inside it would be an import defect.
"""
import ast
import collections
import hashlib
import json
import os
import pathlib
import subprocess
import sys

ROOT = pathlib.Path(__file__).parent.resolve()
PKGS = ("p0", "p1a")
FLAT = ("m8_5b", "m8_5b/production", "m8_5b/pilot", "m8_5b/gates")
ENTRY = ["p0.bundle_operator", "p0.cloud", "p0.group", "p0.representations",
         "p1a.mass_matrix", "p1a.subspace", "p1a.diagnostics",
         "route_a_repn", "route_a_nonabelian", "route_a_twosided"]
# Both frozen documents govern. The addendum SUPERSEDES the parent wherever both
# carry a number for the same quantity, so a gate that verified only the parent
# would leave the governing document unchecked.
FROZEN = (
    ("contract/S1B_DECISION_RULE.md", "<!-- FREEZE-BOUNDARY -->",
     "c44c603a8058ed8529e5bb0f42ec168b443ff22b050d3748f1f1de3537c7d297"),
    ("contract/S1B_ADDENDUM_1.md", "<!-- ADDENDUM-BOUNDARY -->",
     "6da36a1c672772e8c731d4d66f20a77e81779d19d5cb883e617c028d43f46746"),
    ("contract/S1B_ADDENDUM_2.md", "<!-- ADDENDUM2-BOUNDARY -->",
     "14011c338ce331d9bb3424f0ed5aaff31a1c89a9ceafbc2bae41a91ebb41a222"),
    ("contract/S1B_ADDENDUM_3.md", "<!-- ADDENDUM3-BOUNDARY -->",
     "e3304fe993f7a2ea523cfcbc4110c0f37174caea95383b17e465c8fa85b6e28c"),
)
PRIOR_CODE = ("prior/s1b_qualification.py",
              "5a9e04845375c4d12c3a475607f20a1d5f13cc82829d64caa022c82d5e784802")
SENTINEL = "IMPORT_COMPLETE:"


def module_index():
    idx = {}
    for pkg in PKGS:
        for f in sorted((ROOT / pkg).glob("*.py")):
            if f.stem != "__init__":
                idx[f"{pkg}.{f.stem}"] = f
    for d in FLAT:
        for f in sorted((ROOT / d).glob("*.py")):
            idx.setdefault(f.stem, f)
    return idx


def imports_of(path):
    try:
        tree = ast.parse(path.read_text(encoding="utf-8"))
    except Exception:
        return set()
    out = set()
    for n in ast.walk(tree):
        if isinstance(n, ast.Import):
            out.update(a.name for a in n.names)
        elif isinstance(n, ast.ImportFrom):
            out.add(n.module or "")
    return out


def closure(idx):
    """Transitive closure from the DECLARED entry points.

    An entry point absent from `idx` is not skipped: it means the room does not
    supply a module the contract names, which is the exact case a PYTHONPATH
    entry could silently satisfy from outside. main() reports it.
    """
    seen, q = set(), collections.deque(ENTRY)
    while q:
        m = q.popleft()
        if m in seen or m not in idx:
            continue
        seen.add(m)
        for dep in imports_of(idx[m]):
            for c in (dep, dep.split(".")[-1], f"p0.{dep}", f"p1a.{dep}"):
                if c in idx and c not in seen:
                    q.append(c)
    return seen


def probe(mod):
    """Import in a subprocess and demand a sentinel plus an in-room origin."""
    code = (f"import importlib, sys\n"
            f"m = importlib.import_module({mod!r})\n"
            f"print({SENTINEL!r} + {mod!r} + '|' + str(getattr(m, '__file__', '')))\n")
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                       cwd=str(ROOT), env={**os.environ})
    line = next((l for l in r.stdout.splitlines() if l.startswith(SENTINEL)), None)
    if line is None:
        why = "exited without the sentinel" if r.returncode == 0 else \
              (r.stderr.strip().splitlines() or ["(no stderr)"])[-1]
        return False, None, why
    return True, line.split("|", 1)[1], None


# The manifest covers EVERY supplied file, this gate and TASK.md included, and
# excludes only itself. An earlier version excluded TASK.md because TASK.md
# quoted the manifest hash, which was a cycle: the commission authenticated the
# file that was supposed to authenticate the commission, and TASK.md was
# therefore supplied but unpinned. The cycle is broken at the other end instead:
# the manifest hash is given to the unit OUT OF BAND in the handoff prompt, no
# literal appears in any file inside the room, and everything in the room is
# pinned.
MANIFEST_EXCLUDE = {"ROOM_MANIFEST.json"}


def manifest():
    out = {}
    for f in sorted(ROOT.rglob("*")):
        if f.is_file() and "__pycache__" not in f.parts and f.name not in MANIFEST_EXCLUDE:
            out[str(f.relative_to(ROOT))] = hashlib.sha256(f.read_bytes()).hexdigest()
    return out


def main():
    idx = module_index()
    req = closure(idx)
    print(f"  entry points {len(ENTRY)}; REQUIRED closure {len(req)} module(s); "
          f"supplied {len(idx)} module(s)\n")

    bad, scripts = [], []

    # A declared entry point that the room does NOT supply is the silent-
    # substitution case: it drops out of the room scan, so probing only what the
    # scan found would never look at it, and a PYTHONPATH entry outside the room
    # can satisfy the import while the gate reports READY. An earlier version of
    # this gate had exactly that hole, and its own arm caught it. Every declared
    # entry point is therefore probed BY NAME whether or not the scan found it.
    for m in ENTRY:
        if m not in idx:
            ok, origin, why = probe(m)
            if ok:
                bad.append((m, f"NOT SUPPLIED BY THE ROOM, resolved from {origin}"))
            else:
                bad.append((m, f"NOT SUPPLIED BY THE ROOM and does not import: {why}"))

    # Probe ONLY the required closure. Importing the supplied scripts has side
    # effects: three of them run a rehearsal at import and REWRITE JSON files
    # under m8_5b/rehearsal/, so a gate that probed everything mutated the very
    # room it was pinning and the manifest hash moved on every run. Modules
    # outside the closure are listed, not imported.
    for m in sorted(idx):
        if m not in req:
            scripts.append(m)
            continue
        ok, origin, why = probe(m)
        if not ok:
            bad.append((m, why))
        elif not str(pathlib.Path(origin).resolve()).startswith(str(ROOT)):
            bad.append((m, f"ORIGIN OUTSIDE ROOM: {origin}"))

    for m, why in bad:
        print(f"  FAIL  {m}\n        {why}")
    print(f"  {len(scripts)} supplied module(s) outside the required closure: listed, NOT imported "
          f"(three of them rewrite files at import)")
    print(f"\n  {len(bad)} required-module failure(s)")

    chash_ok = True
    for rel, marker, want in FROZEN:
        f = ROOT / rel
        if not f.exists():
            print(f"  MISSING {rel}"); chash_ok = False; continue
        body = f.read_bytes().split(marker.encode())[0]
        got = hashlib.sha256(body).hexdigest()
        ok = got == want
        chash_ok &= ok
        print(f"  frozen {rel:34s} {'OK' if ok else 'MISMATCH: ' + got}")
    rel, want = PRIOR_CODE
    f = ROOT / rel
    if f.exists():
        got = hashlib.sha256(f.read_bytes()).hexdigest()
        ok = got == want
        chash_ok &= ok
        print(f"  prior  {rel:34s} {'OK' if ok else 'MISMATCH: ' + got}")
    else:
        print(f"  MISSING {rel}"); chash_ok = False

    man = manifest()
    (ROOT / "ROOM_MANIFEST.json").write_text(json.dumps(man, indent=1, sort_keys=True) + "\n")
    man_sha = hashlib.sha256((ROOT / "ROOM_MANIFEST.json").read_bytes()).hexdigest()
    print(f"  room manifest: {len(man)} file(s), manifest SHA-256 {man_sha}")

    ok = not bad and chash_ok
    print(f"\n  ROOM {'READY' if ok else 'INCOMPLETE, LAUNCH BLOCKED'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
