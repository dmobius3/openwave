# M7.21, complete + govern the MODELS.md column

> Task detail doc for **M7.21** (Phase 3). Roadmap row: [`m7_roadmap.md`](../m7_roadmap.md) (DONE, § Phase 3). The publication target: repo-root [`MODELS.md`](../../../../../MODELS.md). Canonical spec: [`m7_theory_canonical.md`](../m7_theory_canonical.md).

## TASK PLANNING

**Scope**: publish the 21-cell HydroBoros (M7) column, drafted + staged at the M7.7 milestone (as a preview draft, removed at publication), into the repo-root [`MODELS.md`](../../../../../MODELS.md) cross-model coverage matrix. The M7 research program is (2026-07-07, post-M7.10), pending the author's external pass of the Phase-1 walkthrough); the maintainer's call is to publish the results-so-far now rather than wait, so the column enters honestly (0✅ / 8⚠️ / 13🚧) with every filled cell script-backed. This also re-applies the **new column-ordering rule** (below) across all four models.

**Column-ordering rule (maintainer, 2026-07-18)**: sequence the model columns by their **validated + partial count (✅ + ⚠️)**, descending. Tie-break 1: more **✅** (validated) wins. Tie-break 2: fewer **❌** (honest negatives) wins.

Applying it with the current counts:

| Model | ✅ | ⚠️ | ❌ | ✅+⚠️ | Rank |
| --- | --- | --- | --- | --- | --- |
| Liquid Crystal (M5) | 8 | 8 | 0 | 16 | 1 |
| Ouroboros (M6) | 4 | 5 | 0 | 9 | 2 |
| HydroBoros (M7) | 0 | 8 | 0 | 8 | 3 (ties M4 on 8 and on ✅ = 0; wins on ❌: 0 < 3) |
| EWT (M4) | 0 | 8 | 3 | 8 | 4 |

New column order: **M5, M6, M7, M4** (M4 moves from 3rd to 4th; M7 inserts at 3rd).

**Definition of done**:
- The M7 column is present in all three domain tables (PARTICLES 13 rows, FORCES 5 rows, WAVES 3 rows = 21) plus the Summary Count table, in the M5·M6·M7·M4 order.
- Every filled M7 cell links to a runnable `scripts/m7_N_*.py` or research note, with repo-root-relative paths (MODELS.md sits at repo root).
- The M7 icons match the staged preview exactly (0✅ / 8⚠️ / 0❌ / 0🔶 / 13🚧); no icon inflation.
- The prose that references the model set is updated: the models list, the column-order note, the "Reading the table" framing (three → four frameworks), and the "Per-model results of record" table (M7 graduates out of the "beyond the three scored columns" line into a scored row).
- `check_docs.py` exits 0 over MODELS.md and every touched `.md`.

**Gating**: none (housekeeping / publication). Upstream: M7.7 (drafted + staged the column) and M7.8-M7.10 (the Phase-1-extension work behind the caveats). No `Gated By` dependency.

**Blindspot pass**: skipped: routine (a documentation-migration task; the physics is already earned cell by cell in M7.0-M7.10, this only moves the staged column into the live matrix and re-sorts columns).

**Deviation from the roadmap's original M7.21 scope**: the roadmap framed M7.21 as "open the new-model issue; PR + DCO + light review per governance". Since M7 is now waiting model author, this run does the **content migration** (the column lands in MODELS.md, reordered); the GitHub issue + PR governance flow stays the maintainer's to run (git stays the maintainer's). The canonical-spec refresh #2 named in the roadmap row is not triggered by this migration (no new physics); noted, not done.

**Research body**: the edit lands in repo-root [`MODELS.md`](../../../../../MODELS.md); this task_details is the record; the source of the cells was the M7.7 preview draft (`preview_models.md`, removed at publication to avoid duplication).

**Sub-experiments / preconditions**: none (no scripts run; a re-verification that every M7 cell's linked script/doc path exists is the only check, plus the doc checker).

## TASK REVIEW (2026-07-18)

**Task Duration:** 00:17 (from 11:54 to 12:11 EDT; the post-approval finalization sweep ran after)
**Usage Cap Triggered:** NO

**Results:**
- ✅ M7 HydroBoros column published into repo-root [`MODELS.md`](../../../../../MODELS.md), all 21 cells (13 particles / 5 forces / 3 waves).
- ✅ Column reorder applied by the ✅+⚠️ rule (tie-break: more ✅, then fewer ❌) → **M5 · M6 · M7 · M4**. M7 lands 3rd (ties M4 on 8 and on ✅ = 0, wins on ❌: 0 < 3).
- ✅ M7 icon counts preserved from the M7.7 draft: **0✅ / 8⚠️ / 0❌ / 0🔶 / 13🚧** (script-extracted the M7 column and counted).
- ✅ All 10 M7 link targets resolve, repo-root-relative; 0 missing across all 52 `openwave/` link targets (no M5/M6/M4 link broke).
- ✅ 5-column consistency across the whole coverage region + summary count; doc checker `--changed` clean; no SABER/DHC/market leak (OpenWave public, total invisibility upheld).
- ✅ Prose synced: models list, column-order note, "Reading the table" (three → four), per-model record row, "beyond the N scored columns", and the stale "table spans three → four" line.

**Issues / blockers:** none.

**Deviations from plan:**
- Roadmap framed M7.21 as issue + PR governance; run did the **content migration only** (M7 waiting, publish-directly call). GitHub issue + PR flow stays the maintainer's.
- Fixed en route: task_details MODELS.md links were one level too shallow (`../../../../` → `../../../../../`); removed "SABER" + "user" phrasing from this OpenWave-public doc.
- **Post-approval addition (maintainer directive):** hard-deleted the M7.7-era `preview_models.md` to remove duplication risk, then repointed / reworded all **~20 inbound references** across 8 files (active: briefing, roadmap, question tracker, this doc; frozen history reworded to past tense: `m7_7_consolidation.md` ×5, `m7_phase1_report.md` ×2, `m7_7_call_prep.md` ×1). The gitignored `checkpoints/m7_7_progress.md` mention was left (local-only, plain text).
- **Post-approval cleanup (maintainer directive):** swept the M7 folder's pre-existing dangling `theory/*` links (13, across `m7_background.md`, `m7_0_bootstrap.md`, `m7_1_infra.md`, `m7_2_fleury_torus.md`, `m7_6_observables.md`, `m7_question_tracker.md`), repointing each to the sibling `_CITATIONS.md` manifest (theory docs are local-only for copyright compliance; the manifest is the public index). Result: **zero broken relative links across all 42 M7 files**.

**Action needed:** none outstanding. Commit is the maintainer's (proposed message in the terminal review). GitHub new-model issue + PR governance optional, maintainer's call.

**Findings:** The M7 HydroBoros column is live in the OpenWave cross-model matrix as the 3rd column, entered honestly at 0✅/8⚠️/13🚧 with every filled cell script-backed; the new ✅+⚠️ ordering rule ranks M7 above M4 purely because M4 carries 3 honest negatives to M7's zero. No physics changed, this is publication of the program's results-so-far; the duplicate staging preview was removed so the matrix has a single source.

**Research docs created / updated:**
- [`m7_21_models_column.md`](m7_21_models_column.md) (this record) · [`MODELS.md`](../../../../../MODELS.md) (column + reorder + summary + prose) · [`m7_roadmap.md`](../m7_roadmap.md) (M7.21 → Done, Phase 3) · [`__M7_model_briefing.md`](../../__M7_model_briefing.md) (status → published) · [`m7_question_tracker.md`](../m7_question_tracker.md) (reference repoint)
- Frozen-history reword: [`m7_7_consolidation.md`](m7_7_consolidation.md) · [`m7_phase1_report.md`](m7_phase1_report.md) · [`m7_7_call_prep.md`](m7_7_call_prep.md)
- Removed: `preview_models.md` (the M7.7 staging draft)
