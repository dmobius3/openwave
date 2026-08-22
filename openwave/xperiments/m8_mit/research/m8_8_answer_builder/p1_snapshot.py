"""Build the P1 redline attachment bundle: the exact bytes P1 read.

Target-free machinery. The bundle it writes contains the three public source
artifacts, the two pinned packets, and the governing protocol sections, each
with its SHA-256 recorded, so a reviewer can establish that the cross-check
table faithfully represents the source bytes rather than taking that on trust.

The first P1 redline could not do that: it received the extractor and the
generated table but not the inputs, so "the table matches the sources" was the
one claim it had no way to test.

Section slicing is by heading, not by line number, so the bundle does not go
stale against an edited protocol; a missing heading is an error, not a silent
empty file.
"""

import hashlib
import json
import pathlib
import re
import shutil
import sys
import tempfile

# THE GOVERNING HANDLING SENTENCE for build/bundle/. Defined once, here, and
# repeated verbatim in P1_REPORT.md. Two files previously carried mutually
# exclusive instructions: the report said the bundle was safe for P5, the
# manifest said it was barred from P5. Whichever travelled with the bundle
# decided the answer, which is not a thing a handling label may leave open.
HANDLING = [
    "**HANDLING.** The three source records and both packets copied here are",
    "public and ARE part of the P5 independent reader's declared input set. This",
    "MANIFEST, and every other piece of author-generated packaging, is NOT:",
    "author-side informed review only.",
    "",
    "Public does not mean target-free. The source records STATE THE TARGET VALUES",
    "and trigger leak-scan hits. Nothing here is a § 4.3 permitted clean-room",
    "input, and neither are the originals it copies.",
]


SECTIONS = ("5.2", "5.3", "5.4", "5.5")


def sha256_file(p):
    return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()


def strip_role(name, role):
    """Drop an existing `<role>__` prefix so re-running on a bundle is a fixed
    point. Without this, rebuilding from the bundle produced
    `group__group__...` and a protocol excerpt whose metadata line and hash
    both changed, so the documented tool could not reproduce its own layout."""
    pre = f"{role}__"
    return name[len(pre):] if name.startswith(pre) else name


def slice_section(text, number):
    """The text of `### <number> ...` up to the next heading of any level."""
    m = re.search(rf"^###\s+{re.escape(number)}\s", text, re.M)
    if not m:
        raise ValueError(f"protocol section {number} not found")
    nxt = re.search(r"^##+\s", text[m.end():], re.M)
    return text[m.start():m.end() + (nxt.start() if nxt else len(text))]


def capture_gate_ids(path):
    """Run a source program in a scratch tree and record its gate ids in order."""
    import subprocess
    import tempfile
    p = pathlib.Path(path)
    sand = pathlib.Path(tempfile.mkdtemp(prefix="p1pin-"))
    try:
        (sand / "scripts").mkdir()
        shutil.copy(p, sand / "scripts" / p.name)
        r = subprocess.run([sys.executable, p.name], cwd=str(sand / "scripts"),
                           capture_output=True, text=True, timeout=900,
                           env={"PYTHONDONTWRITEBYTECODE": "1", "PATH": "/usr/bin:/bin"})
        if r.returncode != 0:
            raise RuntimeError(f"{p.name} exited {r.returncode}; refusing to pin it")
        ids = re.findall(r"^\s*(?:PASS|FAIL)\s+\[([^\]]+)\]", r.stdout + r.stderr, re.M)
        if not ids or len(set(ids)) != len(ids):
            raise RuntimeError(f"{p.name} reported {len(ids)} gate ids, not all unique")
        return ids
    finally:
        shutil.rmtree(sand, ignore_errors=True)


def main():
    if len(sys.argv) != 8:
        print("usage: p1_snapshot.py <out_dir> <protocol> <group> <construction> "
              "<theory> <note> <repro>", file=sys.stderr)
        return 2
    dest = pathlib.Path(sys.argv[1])
    protocol = sys.argv[2]
    named = dict(zip(("group", "construction", "theory", "note", "repro"), sys.argv[3:]))

    # Build into a scratch directory and replace the destination wholesale. The
    # previous version reused the destination with exist_ok and never removed
    # extras, so a pre-existing file survived into the finished bundle, was
    # absent from the manifest, and broke the bundle's "exact bytes P1 read"
    # meaning. It could equally have carried stale answer-bearing material.
    stage = pathlib.Path(tempfile.mkdtemp(prefix="p1bundle-"))
    out = stage / "bundle"
    out.mkdir(parents=True)
    lines = ["# P1 attachment bundle: the exact bytes P1 read", ""] + HANDLING + [
             "",
             "Every file below is a byte copy of the artifact P1 used, with its",
             "SHA-256. The two packets are pinned by `adjudicates` and the",
             "protocol sections are the governing text.",
             "",
             "| role | file | sha256 |",
             "| --- | --- | --- |"]

    for role, src in named.items():
        p = pathlib.Path(src)
        dst = out / f"{role}__{strip_role(p.name, role)}"
        shutil.copy(p, dst)
        lines.append(f"| {role} | `{dst.name}` | `{sha256_file(dst)}` |")

    # The FULL protocol ships too. With only the sliced sections present, the
    # manifest named and hashed a source that was absent, so the slicing step
    # could not be re-executed from the bundle alone and the reviewer had to
    # take the excerpt on trust.
    pdst = out / f"protocol__{strip_role(pathlib.Path(protocol).name, 'protocol')}"
    shutil.copy(protocol, pdst)
    lines.append(f"| protocol, full | `{pdst.name}` | `{sha256_file(pdst)}` |")

    ptext = pathlib.Path(protocol).read_text()
    body = "\n".join(slice_section(ptext, s).rstrip() + "\n" for s in SECTIONS)
    sec = out / "protocol__sections_5_2_to_5_5.md"
    sec.write_text(f"# Protocol § {', § '.join(SECTIONS)}, verbatim\n\n"
                   f"Sliced from `{pdst.name}`, sha256 `{sha256_file(protocol)}`.\n"
                   f"Reproduce with: `python3 p1_snapshot.py <out> {pdst.name} ...`\n\n{body}")
    lines.append(f"| protocol sections | `{sec.name}` | `{sha256_file(sec)}` |")

    env = {}
    try:
        import platform
        env = {"python": platform.python_version(), "platform": platform.platform()}
        for mod in ("numpy", "mpmath"):
            env[mod] = __import__(mod).__version__
    except Exception:                                          # noqa: BLE001
        pass

    # The pin file the driver reads. Deliberate state: a source revision changes
    # only when someone re-pins it on purpose.
    #
    # It pins the GATE INVENTORY as well as the hash. Pinning the hash alone left
    # `GATES: n/n pass` as the only evidence a program ran, and a revision that
    # exits early after printing `GATES: 1/1 pass` satisfies that while having
    # run almost nothing.
    gate_ids = {}
    for role in ("theory", "repro"):
        gate_ids[role] = capture_gate_ids(named[role])
        print(f"  captured {len(gate_ids[role])} gate ids for {role}")
    pins = {"_note": "Expected revisions AND gate inventories of the source "
                     "records. p1.py gates SRC_PINNED and SRC_RUN_* against this. "
                     "EXTERNAL TRUST ROOT: this file is caller-selectable, so a "
                     "synchronized source-and-pin replacement is not detectable "
                     "from inside P1. Its hash is recorded in the bundle manifest "
                     "so a reviewer can pin the pin.",
            "sources": {r: sha256_file(named[r]) for r in ("theory", "note", "repro")},
            "gate_ids": gate_ids,
            "environment_when_pinned": env}
    (dest.parent / "sources.pinned.json").write_text(
        json.dumps(pins, indent=2, sort_keys=True) + "\n")

    pinfile = dest.parent / "sources.pinned.json"
    lines += ["", "## Environment when this bundle was built", "",
              "```json", json.dumps(env, indent=2, sort_keys=True), "```", "",
              "Recorded, NOT gated: P1 runs green on other runtimes. Provenance,",
              "not a compatibility guarantee.", "",
              "## External trust root", "",
              f"Source pins live in `../sources.pinned.json`, sha256 "
              f"`{sha256_file(pinfile)}`. That file is caller-selectable, so P1",
              "cannot itself distinguish an approved repin from a synchronized",
              "source-and-pin replacement. Pin it here and review changes to it."]
    (out / "MANIFEST.md").write_text("\n".join(lines) + "\n")

    # The finished filename set must equal the manifest's, exactly.
    listed = set(re.findall(r"^\| [^|]+ \| `([^`]+)`", (out / "MANIFEST.md").read_text(), re.M))
    listed.add("MANIFEST.md")
    present = {f.name for f in out.iterdir()}
    if present != listed:
        raise RuntimeError(f"bundle contents do not match the manifest: "
                           f"unmanifested {sorted(present - listed)}, "
                           f"missing {sorted(listed - present)}")

    if dest.exists():
        shutil.rmtree(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(out), str(dest))
    shutil.rmtree(stage, ignore_errors=True)

    print(f"P1 attachment bundle written to {dest}")
    print(f"  filename set verified against MANIFEST.md: {len(listed)} files, no extras")
    for f in sorted(dest.iterdir()):
        print(f"  {sha256_file(f)[:12]}  {f.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
