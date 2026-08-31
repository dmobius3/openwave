# C2-A1: dated author annotation (2026-08-31, revised at review)

Author-side, written after the halt and REVISED against the maintainer's #508 review.
Filed per the #506 framework. Nothing here is a ruling. Three claims in the first version
of this annotation were wrong; each is corrected below and the correction is named, since
a preservation filing whose description does not match its own bytes is worth less than
no description at all.

## Corrections to the first version of this annotation

1. **The surviving ledger is NOT commingled.** Measured from the filed bytes: 148 GATE
   records, one per `(gate_id, arena_id)` pair, ZERO duplicates; `G1-LIN` appears 72
   times because that is 9 sector families times 8 rungs, each distinct, not a repeated
   ladder; and the running `cumulative_seconds` sum identity holds exactly across all 148
   records (zero monotonic breaks, zero sum violations at 0.01 s). Since `ledger_util.py`
   keeps the cumulative in per-process memory, a second concurrent writer would have
   broken that identity at its first record. The bytes are consistent with ONE sequential
   writer. The first version's "commingled with the parallel re-run" and "G1 sequence
   twice" were my misreading of `(gate_id, rung)` pairs without the arena field, and are
   withdrawn. What does hold: zero RESOURCE records and no timestamps, so the § 13
   nonconformance stands on those grounds.
2. **Execution began at 17:44:40Z, not ~18:0x**, and there were FIVE ledger deletions,
   not two. See the corrected timeline.
3. **The first extraction was incomplete.** It globbed only the main session transcript
   and missed every sidechain sub-agent transcript, which is where the build-phase gate
   executions and three of the five deletions live. `recovered/discipline_events_extracted_v2.json`
   re-extracts over all 8 sources (1 main, 7 sidechain): 120 events across 17 files, versus
   104 across 13. The v1 file is retained under its own name as the record of what was
   first filed. SCOPE, stated rather than claimed exhaustive: v2 covers every `Write`,
   `Edit`, `MultiEdit` targeting a path under the room, and every `Bash` matching kill,
   `rm`, or a run-script invocation, across those 8 transcripts. Anything a transcript
   does not contain, or any file created by a mechanism other than those tools (a Bash
   heredoc, for instance), is outside its reach.

## The corrected timeline

| time (UTC) | event |
| --- | --- |
| 17:26:52Z | `ledger/INPUT_MANIFEST.json` written and closed (verification point 4: BEFORE any gate executes) |
| 17:31:33Z | the close and its coverage statement recorded in `COMMISSIONING.md` |
| **17:44:40Z** | **FIRST GATE RECORD**: a sidechain sub-agent runs `run_gate_1` through `ledger_util`, appending real GATE records to `ledger/OUTPUT_LEDGER.jsonl`. Execution began here under § 11's boundary, whatever the unit intended these as |
| 17:45:05Z | `rm -f ledger/OUTPUT_LEDGER.jsonl` (deletion 1, "Remove test ledger entries") |
| 17:45:12Z | `run_gate_2` executed, appending |
| 17:50:25Z | `rm -f` (deletion 2, "Clean test ledger") |
| 17:50:32Z | `run_gate_10` executed, appending |
| 17:50:40Z | `rm -f` (deletion 3, "Clean test ledger") |
| ~18:0x to 22:39Z | run A (gates 1 to 4; PASS lines with wall times in `recovered/`) |
| 22:39:03.858Z | `rm -f` (deletion 4), then run B starts (PID 86995) |
| 23:09:40.455Z | run B killed |
| 23:10:10.759Z | `rm -f` (deletion 5), then run C starts (PID 88270) |
| 23:14:16Z | `gate_checks.py` edited WHILE run C executes; the running process holds the pre-edit bytes |
| 23:29:56.597Z | run C killed on the author's halt order, at gate 4, N=44 |

**Q1 answered.** The pre-18:00 executions were build-phase tests in the unit's own framing
("Remove test ledger entries"), and they DID append to the real ledger path through the
real `ledger_util` writer, so under § 11's literal boundary the first GATE record sits at
17:44:40Z, eighteen minutes after the manifest close and inside what the unit treated as
still building. Verification point 4 is unaffected: the manifest closed first, and its
recorded hash verifies.

## In-attempt writes and edits to in-room code (v2 extraction, 102 events across 17 files)

| file | events | first | last |
| --- | --- | --- | --- |
| `build/arena.py` | 1 | 17:26:12Z | 17:26:12Z |
| `build/cg_contraction.py` | 11 | 19:42:16Z | 22:34:48Z |
| `build/fast_transform.py` | 3 | 17:53:54Z | 19:06:44Z |
| `build/galerkin.py` | 4 | 17:33:54Z | 19:07:00Z |
| `build/gate_checks.py` | 18 | 17:44:58Z | 23:14:16Z |
| `build/gates.py` | 4 | 17:44:27Z | 18:30:36Z |
| `build/group.py` | 5 | 15:23:56Z | 17:11:46Z |
| `build/integrator.py` | 5 | 18:04:41Z | 21:06:58Z |
| `build/ledger_util.py` | 2 | 17:33:11Z | 19:23:44Z |
| `build/newton.py` | 19 | 20:41:29Z | 21:13:56Z |
| `build/operators.py` | 1 | 17:33:47Z | 17:33:47Z |
| `build/packet.py` | 1 | 17:25:29Z | 17:25:29Z |
| `build/quadrature.py` | 2 | 17:34:35Z | 17:35:57Z |
| `build/sections.py` | 9 | 15:33:15Z | 19:12:47Z |
| `ledger/COMMISSIONING.md` | 6 | 15:10:34Z | 17:31:33Z |
| `run_gates.py` | 1 | 19:06:18Z | 19:06:18Z |
| `run_qualification.py` | 10 | 20:57:45Z | 22:28:04Z |
`run_qualification.py` and `run_gates.py` live at the room ROOT, outside § 12's declared
write areas (`build/`, `ledger/`); writing outside those areas is itself a discipline
fact and is recorded here. Both are now filed in `in_room_as_edited/` at their halt-state
bytes. `gate_checks.py`'s last edit at 23:14:16Z is the mid-run divergence.

## Reading the ledger's mutation_status field

Two conventions coexist in the room and a later reader needs both: `gates.py` writes
`RED` for an arm that FIRED (the arm reddened as designed), while `gate_checks.py` writes
`GREEN` for the same event. In the filed ledger this shows as 108 `RED`, 36
`PARENT-ONLY`, and 4 `GREEN`, the last being G4's records under the second convention.
No record in the file denotes a failed parent check: all 148 carry `parent_status: GREEN`.

## What is recoverable and what is not

RECOVERABLE, filed under `recovered/`: run A and run B stdout fragments with gate PASS
lines and wall times; the pre-edit content of every edited hunk (each Edit event's
`old_string`); the full event timeline; the sidechain transcripts themselves.
NOT RECOVERABLE: the complete `OUTPUT_LEDGER.jsonl` bytes of any deleted state (five
deletions; no tool result captured any of those files whole). The surviving ledger is run
C's, internally consistent as described above, and nonconforming under § 13 for the
absence of RESOURCE records and timestamps.

## Hash coverage

`package_hashes.txt` beside this file is generated last, in one shot, from the tracked
files; each hash covers the named file's bytes as committed in this PR. That sentence is
the coverage statement for every hash in the manifest.
