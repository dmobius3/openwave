# S14-C2: the M8.5-C2 pre-freeze obligations, execution record

Author-side, begun 2026-08-29 on F1's D-7-closed GO. Draft under execution:
`../M8_5C2_PROTOCOL_DRAFT.md` (the § 14 dispositions ruled at D-5: items 1, 2, 3, 8
FRESH; 9 a-through-m and 10 NEW; 4 through 7 INHERITED, item 5 mechanically).
Environment: python 3.13.13, numpy 2.5.0, scipy 1.18.0, same machine as the A1 room and
the C1 s14 package.

## Status

| item | disposition | record | verdict |
| --- | --- | --- | --- |
| 1 | FRESH | `raw/c2_rerun_right_translation_check.out` | GREEN, exit 0, all nine report lines, script at pinned hash `752baa2d…` |
| 2 | FRESH | `raw/c2_rerun_{mode_count,jacobian_check,exact_quad_check,cascade_quad_check}.out` | GREEN, exit 0 all four, scripts at pinned hashes |
| 3 | FRESH (D-5: non-inheritable by its own text) | `c2_item3_wallclocks.py` (scratch, C1's code repointed to a C2-only output name), `raw/c2_item3_wallclocks.json`, `raw/c2_item3_run.out` | GREEN: composed agreement 0.72 h + monitor 1.74 h + Control B 1.74 h = 4.20 h against 48, margin 11x, FEASIBILITY-ONLY |
| 4 | INHERITED | C1 package `research/m8_5c/s14/raw/lattice_tables.json` at § 15 pin `efe6c2d6…` | stands |
| 5 | INHERITED MECHANICALLY | `raw/c2_item5_mechanical.out` | GREEN: block digest recomputed on the C2 draft = § 15 row = the reference file's own stamp, all three `823e9066…` |
| 6 | INHERITED | C1 `raw/item6_v2_amended.json` at § 15 pin `6615b396…` | stands |
| 7 | INHERITED | C1 `raw/item7_denominators.json` | stands |
| 8 | FRESH re-read | `raw/c2_item8_crosscheck_v2.out`, checker `c2_item8_check.py` (self-arming, three mutations all fire) | GREEN at the C1-audit surface assignment: (a)(b) at §§ 2 and 5, (c) at § 2 and gate 5's row, TRIGGER L at § 2 and the § 14 classes |
| 9 | NEW, a through m | `c2_item9_arm_demos.py` at `42ff8da8a3cf0d808ac4dad7bec8fac5a41cf5a686fbb82a576e47d3c6d2001e` (FULL value recomputed from the shipped bytes at pinning, per the anchor's own instruction and F2's transcribed-hash caution); PINNED record `raw/c2_item9_record_run5_PINNED.json` (`8ce491c2771813ea51f2708f2d1d4040585af3de485247e8f3126abbf2c9c302`), stdout `raw/c2_item9_run5_PINNED.out` (measurement record) | GREEN: all thirteen demonstrations, wall 2470 s; five-run chain preserved (1 RED by design, 2 CANDIDATE, 3 KILLED-thrash, 4 STOPPED-at-anchor, 5 PINNED); gate-6 leak per-level margins 5.2e-4 / 3.9e-3 / 1.0e-2 settle the marginal-member question for `{N+1, N+3, N+5}` |
| 10 | NEW | `c2_item10_audit.md` v2 (redline fixes: row 5's layer-1 re-runs moved to a credit-free provenance column per § 2's rule, the same inferential bridge removed once already; row 3's resolution named explicitly) | GREEN: every arm arena-named, every citation resolving to a green record, both regression controls named and exempt, row 11 exempt by construction |

## Dispositions of note

`raw/c2_item8_crosscheck.out` (v1) is RETAINED as the record of a CHECKER defect, not a
text defect: its probe truncated § 2's trigger text at 700 characters (missing clause (c)
beyond the window) and asserted § 5 carries clause (c), which the C1 audit record says
§ 5 omits BY DESIGN. The v2 checker keys on the true surface assignment and is
self-arming; v1 was never evidence about the text and is kept under its own name per the
never-redirect-a-record rule.

Item 3's fresh measurement (4.20 h) sits beside the C1 record's 4.35 h, same machine,
consistent; the C1 record is untouched and the C2 scratch writes only C2-named files.

## § 14 CLOSED (2026-08-29, all ten items disposed)

Items 1, 2, 3, 8 FRESH and green; 5 INHERITED MECHANICALLY (three-way digest equality); 4,
6, 7 INHERITED at the C1 pins; 9 and 10 NEW and green with the pinned records above. The
ANCHOR section below is now historical; its single-next-step was executed as run 5.

## ANCHOR (2026-08-29): item 9 paused mid-execution, state for pickup

The author stopped work deliberately (away from the machine); nothing is wrong. The single
next step: relaunch the item-9 script AS-IS and let it finish; on all-green it pins as the
item-9 record and item 10's audit follows.

Command, from this directory:
    python3 c2_item9_arm_demos.py > raw/c2_item9_run.out 2>&1

Script bytes at anchor: `c2_item9_arm_demos.py` SHA-256 `42ff8da8a3cf0d808ac4…` (recompute
the full value before pinning). Expected wall 60 to 90 minutes; peak memory ~2 GB with the
chunk guard.

The run chain so far, every record preserved under its own name in `raw/`:
| run | files | outcome |
| --- | --- | --- |
| 1 | `c2_item9_record_run1_RED.json`, `c2_item9_run1_RED.out` | RED by design: caught the § 4.2(b) false sentence (leak read levels), the kron-convention quotient bug, the distinct-level injection bug, and the low-read exactness fact |
| 2 | `c2_item9_record_run2_CANDIDATE.json`, `c2_item9_run2_CANDIDATE.out` | CANDIDATE, 13 of 14 green; exposed the masked coefficient-rotation transpose; measured the rung-realization leak margin ~1e-2 |
| 3 | `c2_item9_record_run3_KILLED.json`, `c2_item9_run3_KILLED.out` | KILLED mid-(f): memory thrash (9.6 GB RSS on 24 GB; the level-29 read table's chunk array); first 10 sections all green incl. (l) post-transpose-fix |
| 4 | `c2_item9_record_run4_STOPPED.json`, `c2_item9_run4_STOPPED.out` | STOPPED by the author's anchor, not a defect; bytes identical to the relaunch |

Draft state at anchor: § 4.2(b) and row 6 carry the RUNG-RELATIVE read levels
`n ∈ {N+1, N+3, N+5}` per both redlines' blocker; the orthonormalizer is the symmetric
Löwdin map with a symmetric second-pass refinement per F1's ruling (no triangular factor);
register/diff regenerated at 19 hunks, 197 additions, 51 deletions. F2's row-6 D-7 close
is superseded by the amendment; the row-6 re-read completes against the amended text with
the pinned run's record. After item 9 pins: item 10's audit, § 15 pin fill (script hash +
record hash), then the freeze sequence.
