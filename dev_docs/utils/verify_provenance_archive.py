#!/usr/bin/env python3
"""Verify an out-of-band provenance archive against its own manifest.

A clean-room reproduction protocol delivers its derivation to the maintainer
outside the repository, because the derivation bytes must not enter the base the
room opens from. What arrives is an archive plus an author-supplied manifest of
hashes. The maintainer's job is to recompute all of it, and specifically NOT to
run the self-check that ships inside the archive: a checker running inside the
artifact it verifies is the one gate an audit cannot delegate, since it can
report success by construction.

This script is that recomputation. It never imports, executes, or trusts
anything from the archive. It reads the manifest as data.

What it checks:

  1. the archive file's own hash, when one is given to check against;
  2. every manifest-listed file, hash and byte count;
  3. closure both ways, no unlisted file present and no listed file absent;
  4. the ordered subset's concatenation, compared against the value declared in
     the PACKET rather than the value restated in the manifest;
  5. every cross-hash the packet declares about itself, resolved by content;
  6. orphan hash references, any 64-hex value written anywhere in the archive
     that corresponds to nothing the archive contains.

Check 4 is the reason this is worth automating. The manifest supplies the recipe
(which files, in what order) and the packet supplies the target. They are
separate documents produced at separate times, so their agreement is evidence,
while a manifest agreeing with itself is not.

Check 6 catches the failure that a per-file hash sweep cannot see: an artifact
left in the archive from a superseded run still pins the object it was computed
against, and that pin resolves to nothing. Every listed file can hash correctly
while one of them describes a different construction entirely.

ORDER MATTERS. Run this BEFORE rerunning anything from the archive. Derivation
scripts commonly write their outputs next to themselves, so a rerun in the
extracted tree can overwrite a manifest-listed file and turn a clean archive
into a failing one, with the audit as the cause.

Manifest schema, kept loose on purpose so a future packet format still parses:
a row of `| path | sha256 | bytes |` in backticks is a listed file, and a
leading integer column marks membership in the ordered set, in that order.

Usage:
    python3 dev_docs/utils/verify_provenance_archive.py DIR
    python3 dev_docs/utils/verify_provenance_archive.py --archive FILE.tar.gz
    python3 dev_docs/utils/verify_provenance_archive.py --archive FILE.tar.gz \\
        --expect-sha 0123abcd... --no-quantities --strict
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import tarfile
import tempfile
from pathlib import Path

HEX64 = re.compile(r"\b([0-9a-f]{64})\b")
ROW = re.compile(
    r"^\|\s*(?:(\d+)\s*\|\s*)?`([^`]+)`\s*\|\s*`([0-9a-f]{64})`\s*\|\s*(\d+)\s*\|\s*$"
)
TEXT_SUFFIXES = {".md", ".py", ".txt", ".json", ".csv", ".toml", ".cfg", ".yml", ".yaml"}

# Decimal renderings with at least this many places count as a reproduced
# quantity under --no-quantities. Version strings are excluded separately.
DECIMAL_PLACES = 3
VERSIONISH = re.compile(r"\d+\.\d+\.\d+")


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


class Result:
    def __init__(self) -> None:
        self.failures: list[str] = []
        self.warnings: list[str] = []

    def check(self, ok: bool, label: str, detail: str = "") -> bool:
        mark = "✅" if ok else "❌"
        print(f"  {mark} {label}" + (f"   {detail}" if detail else ""))
        if not ok:
            self.failures.append(label)
        return ok

    def warn(self, label: str, detail: str = "") -> None:
        print(f"  ⚠️  {label}" + (f"   {detail}" if detail else ""))
        self.warnings.append(label)


def parse_manifest(text: str) -> tuple[list[tuple[str, str, int]], list[str]]:
    """Return (listed files, ordered subset paths in order)."""
    listed: list[tuple[str, str, int]] = []
    ordered: list[tuple[int, str]] = []
    for line in text.splitlines():
        m = ROW.match(line)
        if not m:
            continue
        index, path, digest, size = m.group(1), m.group(2), m.group(3), int(m.group(4))
        listed.append((path, digest, size))
        if index is not None:
            ordered.append((int(index), path))
    return listed, [p for _, p in sorted(ordered)]


def find_packet(root: Path, listed: list[tuple[str, str, int]]) -> tuple[Path, dict] | None:
    """The packet is whichever listed JSON declares a provenance id. Found by
    schema rather than by filename, so a renamed packet still resolves."""
    for path, _, _ in listed:
        f = root / path
        if f.suffix != ".json" or not f.is_file():
            continue
        try:
            doc = json.loads(f.read_text())
        except (ValueError, UnicodeDecodeError):
            continue
        if isinstance(doc, dict) and isinstance(doc.get("provenance_id"), dict):
            return f, doc
    return None


def verify(root: Path, args: argparse.Namespace, r: Result) -> None:
    manifest = root / "MANIFEST.md"
    if not manifest.is_file():
        r.check(False, "MANIFEST.md present at the archive root")
        return

    listed, ordered = parse_manifest(manifest.read_text())
    print(f"\n== manifest ==\n  {len(listed)} file(s) listed, {len(ordered)} in the ordered set")

    print("\n== 1. every listed file, recomputed ==")
    digests: dict[str, str] = {}
    bad = 0
    for path, want_sha, want_len in listed:
        f = root / path
        if not f.is_file():
            print(f"  ❌ missing: {path}")
            bad += 1
            continue
        blob = f.read_bytes()
        digests[path] = sha256(blob)
        if digests[path] != want_sha or len(blob) != want_len:
            print(f"  ❌ {path}")
            print(f"       sha want {want_sha}")
            print(f"       sha got  {digests[path]}")
            print(f"       len want {want_len}  got {len(blob)}")
            bad += 1
    r.check(bad == 0, f"{len(listed) - bad} of {len(listed)} listed file(s) verify")

    print("\n== 2. closure, the direction a manifest cannot check on itself ==")
    on_disk = {str(p.relative_to(root)) for p in root.rglob("*") if p.is_file()}
    names = {p for p, _, _ in listed} | {"MANIFEST.md"}
    extra, absent = sorted(on_disk - names), sorted(names - on_disk)
    r.check(not extra, "no unlisted file rode along", "" if not extra else ", ".join(extra))
    r.check(not absent, "no listed file is absent", "" if not absent else ", ".join(absent))

    print("\n== 3. the ordered set, against the packet's declared value ==")
    if not ordered:
        r.warn("the manifest declares no ordered set, so nothing pins the source content")
        concat = None
    else:
        h = hashlib.sha256()
        total = 0
        for rel in ordered:
            blob = (root / rel).read_bytes()
            h.update(blob)
            total += len(blob)
        concat = h.hexdigest()
        print(f"  concatenation of {len(ordered)} file(s), {total} bytes")
        print(f"  {concat}")

    found = find_packet(root, listed)
    if found is None:
        r.warn("no listed JSON declares a provenance id, so check 3 has no independent target")
    else:
        packet, doc = found
        rel = str(packet.relative_to(root))
        declared = doc["provenance_id"].get("source_content_sha256")
        print(f"  packet: {rel}  id {doc['provenance_id'].get('id', '(none)')}")
        if concat is not None and declared:
            r.check(
                concat == declared,
                "the ordered set reproduces the packet's source_content_sha256",
            )
        elif not declared:
            r.warn("the packet declares no source_content_sha256")

        print("\n== 4. the packet's other declared cross-hashes, resolved by content ==")
        by_digest = {v: k for k, v in digests.items()}
        cross = {k: v for k, v in doc.items() if k.endswith("_sha256") and isinstance(v, str)}
        if not cross:
            print("  (none declared)")
        for key, value in sorted(cross.items()):
            owner = by_digest.get(value)
            r.check(owner is not None, f"{key} resolves to a file in the archive", owner or value)

    print("\n== 5. orphan hash references ==")
    known = set(digests.values()) | ({concat} if concat else set())
    if args.expect_sha:
        known.add(args.expect_sha)
    orphans: list[str] = []
    for path, _, _ in listed:
        f = root / path
        if f.suffix.lower() not in TEXT_SUFFIXES:
            continue
        try:
            content = f.read_text(errors="replace")
        except OSError:
            continue
        for n, line in enumerate(content.splitlines(), 1):
            for digest in HEX64.findall(line):
                if digest not in known:
                    orphans.append(f"{path}:{n}: {digest}")
    if orphans:
        for o in orphans[:20]:
            print(f"     {o}")
        if len(orphans) > 20:
            print(f"     ... and {len(orphans) - 20} more")
        msg = "hash reference(s) resolve to nothing in this archive"
        if args.strict:
            r.check(False, f"{len(orphans)} {msg}")
        else:
            r.warn(f"{len(orphans)} {msg}", "a superseded artifact may have been left behind")
    else:
        r.check(True, "every hash written in the archive resolves to something it contains")

    if args.no_quantities:
        print("\n== 6. no reproduced quantity ==")
        hits: list[str] = []
        for path, _, _ in listed:
            f = root / path
            if f.suffix.lower() not in TEXT_SUFFIXES:
                continue
            try:
                content = f.read_text(errors="replace")
            except OSError:
                continue
            for n, line in enumerate(content.splitlines(), 1):
                stripped = VERSIONISH.sub("", line)
                if re.search(rf"\d+\.\d{{{DECIMAL_PLACES},}}", stripped):
                    hits.append(f"{path}:{n}: {line.strip()[:100]}")
        for h in hits[:20]:
            print(f"     {h}")
        r.check(not hits, f"no decimal rendering of {DECIMAL_PLACES}+ places in any listed file")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("directory", nargs="?", help="extracted archive root, holding MANIFEST.md")
    ap.add_argument(
        "--archive", help="tar.gz to hash, screen and extract to a temporary directory"
    )
    ap.add_argument("--expect-sha", help="sha256 the --archive file must have")
    ap.add_argument(
        "--no-quantities",
        action="store_true",
        help="assert the archive contains no decimal rendering, for clean-room deliveries",
    )
    ap.add_argument("--strict", action="store_true", help="treat warnings as failures")
    args = ap.parse_args()

    if bool(args.directory) == bool(args.archive):
        ap.error("give exactly one of DIRECTORY or --archive")

    r = Result()
    tmp: tempfile.TemporaryDirectory | None = None

    if args.archive:
        blob = Path(args.archive).read_bytes()
        got = sha256(blob)
        print(f"\n== 0. the archive file ==\n  {got}   {len(blob)} bytes")
        if args.expect_sha:
            r.check(got == args.expect_sha, "archive sha256 matches the expected value")
        with tarfile.open(args.archive) as tar:
            members = tar.getnames()
            unsafe = [m for m in members if m.startswith("/") or ".." in Path(m).parts]
            r.check(not unsafe, "no member escapes the extraction root", ", ".join(unsafe[:3]))
            if unsafe:
                return 1
            tmp = tempfile.TemporaryDirectory(prefix="provenance-verify-")
            tar.extractall(tmp.name, filter="data")
        roots = [p for p in Path(tmp.name).iterdir() if p.is_dir()]
        root = roots[0] if len(roots) == 1 else Path(tmp.name)
    else:
        root = Path(args.directory).resolve()

    verify(root, args, r)

    if tmp is not None:
        tmp.cleanup()

    failed = r.failures or (r.warnings and args.strict)
    print()
    if failed:
        print(
            f"❌ archive verification FAILED: {len(r.failures)} failure(s), "
            f"{len(r.warnings)} warning(s)"
        )
        for f in r.failures:
            print(f"    {f}")
        return 1
    if r.warnings:
        print(f"⚠️  archive verified with {len(r.warnings)} warning(s)")
        return 0
    print("✅ archive verified: every hash recomputed, closure both ways, nothing orphaned")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
