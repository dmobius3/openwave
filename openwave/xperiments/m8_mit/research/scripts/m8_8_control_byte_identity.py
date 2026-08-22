#!/usr/bin/env python3
"""M8.8 byte-identity control: the repaired harness changed one string and no record byte.

WHAT THIS ESTABLISHES.  The first adjudication attempt ran the committed pre-reveal harness
at `ea3452d3` and was REFUSED at the section 8 step-7 packet-domain gate (STRUCTURAL
FAILURE; record `data/m8_8_adjudication_attempt1_refusal.json`).  The repair merged as
PR #457 replaced the source-domain literal and extracted the terminal printer.  This control
shows the rerun's record is the result of the one-string correction ALONE: it takes the
`ea3452d3` harness bytes from the commit graph, substitutes ONLY the literal, runs that
scratch copy on the published answer packet, and compares the record it writes, byte for
byte, against the committed rerun record.  The scratch copy still carries the `ea3452d3`
printer defect, so it writes the record and then exits nonzero on `KeyError: 'recomputed'`;
that exit is expected and is itself evidence the printer fix is downstream of the record.

It does NOT make attempt 1 retroactively green, and it claims no result category.

    python3 m8_8_control_byte_identity.py          # exit 0 iff the records are identical

Reproducible from public objects only: the commit `ea3452d3`, the published packet, and
the committed record.  No timestamp enters the record.
"""
from __future__ import annotations

import hashlib
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parent
REPO = RESEARCH.parents[2]
HARNESS_REL = "openwave/xperiments/m8_mit/research/scripts/m8_8_adjudication.py"
REFUSED_HARNESS_COMMIT = "ea3452d3"
PLACEHOLDER = "answer_packet.rows"               # the PR #453 maintainer placeholder
SOURCE_DOMAIN = "packet_rows_canonical_position"  # the frozen pre-reveal builder literal
PACKET = RESEARCH / "data" / "m8_8_answer_packet.json"
RECORD = RESEARCH / "data" / "m8_8_adjudication.json"
EXPECTED_SUBSTITUTIONS = 4


def sha256(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def main() -> int:
    src = subprocess.run(
        ["git", "-C", str(REPO), "show", f"{REFUSED_HARNESS_COMMIT}:{HARNESS_REL}"],
        check=True, capture_output=True, text=True,
    ).stdout
    n = src.count(PLACEHOLDER)
    if n != EXPECTED_SUBSTITUTIONS:
        print(f"FAIL: expected {EXPECTED_SUBSTITUTIONS} placeholder sites at "
              f"{REFUSED_HARNESS_COMMIT}, found {n}")
        return 1
    patched = src.replace(PLACEHOLDER, SOURCE_DOMAIN)
    with tempfile.TemporaryDirectory() as td:
        # the harness resolves its inputs relative to its own location, so it runs in place
        scratch = HERE / "_m8_8_control_scratch_ea3452d3.py"
        out = Path(td) / "control_record.json"
        try:
            scratch.write_text(patched)
            proc = subprocess.run(
                [sys.executable, str(scratch), "--packet", str(PACKET), "--json", str(out)],
                capture_output=True, text=True,
            )
        finally:
            scratch.unlink(missing_ok=True)
        if not out.exists():
            print("FAIL: the scratch harness wrote no record")
            print(proc.stderr[-2000:])
            return 1
        got, want = out.read_bytes(), RECORD.read_bytes()
    printer_died = "KeyError: 'recomputed'" in proc.stderr
    print(f"scratch harness   {REFUSED_HARNESS_COMMIT} + {n} substitutions "
          f"{PLACEHOLDER!r} -> {SOURCE_DOMAIN!r}")
    print(f"scratch exit      {proc.returncode} "
          f"({'printer KeyError after the write, as at ea3452d3' if printer_died else 'no printer defect seen'})")
    print(f"control record    sha256 {sha256(got)}  ({len(got)} bytes)")
    print(f"committed record  sha256 {sha256(want)}  ({len(want)} bytes)")
    same = got == want
    print("BYTE-IDENTICAL" if same else "DIFFER")
    return 0 if same and printer_died else 1


if __name__ == "__main__":
    sys.exit(main())
