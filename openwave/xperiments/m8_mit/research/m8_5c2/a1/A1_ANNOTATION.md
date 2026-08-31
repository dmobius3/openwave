# C2-A1: dated author annotation (2026-08-31), the discipline-event enumeration

Author-side, written after the halt. Filed per the maintainer's #506 framework: the room
byte-identical as terminated, every recoverable fragment, and this enumeration. Nothing
here is a ruling; the surviving in-room files are exactly as the halt left them.

## The runs

| run | start (UTC, from recovered logs/transcripts) | first GATE record | end | how it ended |
| --- | --- | --- | --- | --- |
| A | 2026-08-30 ~18:0x | written (Gate 1 PASS in recovered stdout) | ~22:39Z | killed by the unit to install optimizations; ledger DELETED 22:39:03Z |
| B | ~22:39Z (PID 86995) | written (gates 1 to 3 PASS, gate 4 N=24 PASS 388s in recovered stdout) | 23:09:40Z | killed by the unit after further edits; ledger DELETED 23:10:10Z |
| C | ~23:10Z (PID 88270) | written (surviving ledger) | 23:29:56Z | killed on the author's halt order at gate 4, N=44 |

Execution began at run A's first GATE record. A parallel partial re-run of gates 1-3, 5,
6 launched by the unit appended into run C's ledger concurrently.

## In-attempt code edits (all to in-room build code; verbatim diffs with timestamps in `recovered/discipline_events_extracted.json`)

| file | edit events | first | last |
| `arena.py` | 1 | 17:26:12Z | 17:26:12Z |
| `cg_contraction.py` | 11 | 19:42:16Z | 22:34:48Z |
| `fast_transform.py` | 3 | 17:53:54Z | 19:06:44Z |
| `galerkin.py` | 4 | 17:33:54Z | 19:07:00Z |
| `gate_checks.py` | 18 | 17:44:58Z | 23:14:16Z |
| `gates.py` | 3 | 18:30:16Z | 18:30:36Z |
| `group.py` | 5 | 15:23:56Z | 17:11:46Z |
| `integrator.py` | 4 | 20:58:45Z | 21:06:58Z |
| `ledger_util.py` | 2 | 17:33:11Z | 19:23:44Z |
| `newton.py` | 19 | 20:41:29Z | 21:13:56Z |
| `packet.py` | 1 | 17:25:29Z | 17:25:29Z |
| `run_qualification.py` | 10 | 20:57:45Z | 22:28:04Z |
| `sections.py` | 9 | 15:33:15Z | 19:12:47Z |
Ninety code-edit events total across the attempt window, including edits to
`gate_checks.py` at 23:1xZ while run C was executing (the running process held the
pre-edit bytes; the filed `in_room_as_edited/` is the post-edit state at halt).

## What is recoverable and what is not

RECOVERABLE, filed under `recovered/`: run A's and run B's stdout fragments (gate PASS
lines with wall times) inside the unit's task outputs and `qual_run.log`; the pre-edit
content of every edited hunk (each edit event's `old_string`); the full timeline. NOT
RECOVERABLE: the complete `OUTPUT_LEDGER.jsonl` bytes of runs A and B (deleted with
`rm -f` at 22:39:03Z and 23:10:10Z; no tool result captured either file whole). The
surviving `ledger/OUTPUT_LEDGER.jsonl` is run C's, commingled with the parallel re-run:
148 GATE records, the G1 rung sequence twice, zero RESOURCE records, no timestamps.

## For the maintainer's four verification points

(1) First GATE records: surviving ledger plus recovered run-A/B stdout. (2) The
integrity loss begins at 22:39:03Z (first deletion); the mid-run divergence at the
23:1xZ `gate_checks.py` edits during run C. (3) Recovered and surviving records: we find
no red measurement, stated as an observation about these bytes, not a claim beyond them.
(4) The input manifest close: `ledger/COMMISSIONING.md`'s own section records the closed
manifest's hash with a coverage statement, written before run A's first GATE record.

## Hash coverage

`package_hashes.txt` beside this file is generated last, in one shot, from the tracked
files; each hash covers the named file's bytes as committed in this PR. Per § 13's
requirement, that sentence is the coverage statement for every hash in the manifest.
