#!/usr/bin/env python3
"""Recovery wrapper for the M8.8 construction-packet provenance material, candidate 2.

Unlike candidate 1, whose sources hard-coded a repo path on the author's
machine, the frozen sources here contain no machine-specific paths, so the
derivation is location-portable and needs no relocation step. This wrapper
demonstrates that rather than asserting it: it copies the archive to a fresh
temporary directory and runs every gate there, with the whole manifest re-hashed
on exit. (Deliberate portable paths in the documentation, /tmp scratch dirs and
the container mount, are not machine-specific and are not what portability is
about.)

THIS IS AN AUTHOR-SUPPLIED SELF-CHECK. Its verdict is informational only and
cannot serve as the audit: it runs inside the artifact it verifies. Everything
it checks is recomputable independently from MANIFEST.md; see ENVIRONMENT.md
for the commands.

Usage:  python3 recover.py [--keep] [--skip-mutations]
        --skip-mutations drops the slowest step (the maintainer audit's
        mutation suite, ~11 full gate runs); the default runs it, because a
        verification pass that never demonstrates its checks CAN fail is
        exactly the shape this project keeps getting caught by.
"""
import hashlib, json, os, pathlib, shutil, subprocess, sys, tempfile
HERE = pathlib.Path(__file__).resolve().parent
HASHED = ["qphi.py", "complex.py", "kernel.py", "sat.py", "search_sym.py", "build_packet_v2.py"]
A_PUBLISHED = "e3b0c945bbbb15b4549fa641234c9461062c2337b3d1e372af621b614d4883a9"
ENV = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
fail = []
def chk(label, cond, detail=""):
    print(f"  [{'PASS' if cond else '**FAIL**'}] {label}" + (f"   {detail}" if detail else ""))
    if not cond: fail.append(label)

pk = json.loads((HERE / "m8_8_construction_packet.json").read_text())
import re
MAN = (HERE / "MANIFEST.md").read_text()
manifest_rows = re.findall(r"`([^`]+\.(?:py|json|md|txt|sh))` \| `([0-9a-f]{64})`", MAN)
before = {f: hashlib.sha256((HERE / f).read_bytes()).hexdigest() for f, _ in manifest_rows}

print("\n== 1. frozen derivation source ==")
got = hashlib.sha256(b"".join((HERE / "build" / f).read_bytes() for f in HASHED)).hexdigest()
chk("the six canonical files reproduce the packet's source_content_sha256",
    got == pk["provenance_id"]["source_content_sha256"], got[:32] + "...")

print("\n== 2. group packet ==")
a = hashlib.sha256((HERE / "m8_5a_packet.json").read_bytes()).hexdigest()
chk("local M8.5-A copy matches its PUBLISHED hash", a == A_PUBLISHED, a[:32] + "...")
chk("and matches group_packet_sha256 in the construction packet", a == pk["group_packet_sha256"])

print("\n== 3. location portability, demonstrated not asserted ==")
tmp = pathlib.Path(tempfile.mkdtemp(prefix="m88-recover-"))
shutil.copytree(HERE / "build", tmp / "build")
for f in ("m8_8_construction_packet.json", "m8_5a_packet.json"):
    shutil.copy(HERE / f, tmp / "build" / f)
# Needles are built at runtime so this file does not itself contain the strings it
# hunts, which also lets it scan itself. Scans EVERY text file in the archive, root
# included, not only *.py: an earlier version globbed *.py only and walked straight
# past a .txt carrying an absolute path in the same directory.
#
# Scope, stated exactly: this hunts MACHINE-SPECIFIC paths (user/home directories and
# macOS temp), which are what break location portability. Portable absolute paths that
# the docs use deliberately (/tmp scratch dirs, the container mount /w) are in scope
# for a human read, not for this needle set, and the check label says so.
needles = [chr(47) + w + chr(47) for w in ("Users", "home", "opt", "private")] \
        + [chr(47) + "var" + chr(47) + "folders"] + [chr(126) + chr(47)]   # tilde-home prefix: machine-specific in spirit, missed once
residual = []
for f in sorted(HERE.rglob("*")):
    if not f.is_file() or f.suffix not in (".py", ".md", ".txt", ".json", ".sh"):
        continue
    txt = f.read_text(errors="replace")
    if any(n in txt for n in needles):
        residual.append(str(f.relative_to(HERE)))
chk("no machine-specific path (user/home dirs) in any text file, self included",
    not residual, "" if not residual else f"found in: {', '.join(residual)}")

print("\n== 4. gates, from the relocated copy ==")
for label, script, args in (
        ("certify.py", "certify.py", []),
        ("maintainer audit A1-A11", "m8_8_packet_audit.py",
         ["--packet", "m8_8_construction_packet.json",
          "--group-packet", "m8_5a_packet.json", "--no-write"])):
    r = subprocess.run([sys.executable, script, *args], cwd=tmp / "build",
                       capture_output=True, text=True, env=ENV)
    tail = [l for l in r.stdout.splitlines() if l.strip()]
    chk(f"{label} from the relocated copy", r.returncode == 0,
        tail[-1].strip() if tail else "(no output)")

if "--skip-mutations" in sys.argv:
    print("\n== 4b. mutation suite: SKIPPED on request ==")
    print("  the verdict below therefore does NOT demonstrate the gates can fail")
else:
    print("\n== 4b. the maintainer audit's mutation suite, from the relocated copy ==")
    r = subprocess.run([sys.executable, "m8_8_packet_audit.py",
                        "--packet", "m8_8_construction_packet.json",
                        "--group-packet", "m8_5a_packet.json",
                        "--no-write", "--mutation-tests"],
                       cwd=tmp / "build", capture_output=True, text=True, env=ENV)
    det = [l.strip() for l in r.stdout.splitlines()
           if l.strip().startswith(("DETECTED", "MISSED"))]   # per-mutation lines only,
    # not the "MUTATIONS ALL DETECTED" summary, which once inflated this count to 11
    chk("every mutation detected, so the gates demonstrably CAN fail",
        r.returncode == 0 and det and all(not l.startswith("MISSED") for l in det),
        f"{len(det)} mutation(s) exercised")

print("\n== 5. nothing in the tree was mutated by this run ==")
after = {f: hashlib.sha256((HERE / f).read_bytes()).hexdigest() for f, _ in manifest_rows}
chk(f"all {len(manifest_rows)} manifest-listed files are byte-identical to before this run",
    after == before,
    "" if after == before else
    "changed: " + ", ".join(f for f in after if after[f] != before[f]))
bad = [f for f, h in manifest_rows if after.get(f) != h]
chk("and every one still matches the manifest itself", not bad,
    "" if not bad else f"mismatch: {', '.join(bad)}")

if "--keep" in sys.argv:
    print(f"\n  copy kept at {tmp}")
else:
    shutil.rmtree(tmp); print("\n  temporary copy removed")
print(f"\n  {'RECOVERY VERIFIED (informational; recompute independently per ENVIRONMENT.md)' if not fail else 'FAILED: ' + '; '.join(fail)}\n")
sys.exit(0 if not fail else 1)
