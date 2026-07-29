# T3: MODELS.md cell budget, re-derive the rule

> Roadmap row: [`../platform_roadmap.md`](../platform_roadmap.md) T3. Proposed by the user
> 2026-07-28 while [PR #362](https://github.com/openwave-labs/openwave/pull/362) was open
> with one cell over budget, and run the same day once both gates cleared (§ Gating).
>
> **Outcome: the budget moves 55 → 65**, on the column-geometry argument rather than the
> parity argument the task was filed on, and the number now sits inside a measured bracket
> instead of being inherited (§ FINDINGS F1).

## TASK PLANNING

### Scope

**Re-derive the cell-budget rule for the table structure MODELS.md actually has today,
then update every place that states it.** The expected outcome is 65 words, matching the
`Description` budget in [`ROADMAP_STANDARDS.md`](../ROADMAP_STANDARDS.md) so the platform
carries ONE table-cell budget, but the task is a review with an expected answer, not a
rubber stamp: if the re-derivation says something other than 65, that is the result.

The review has three parts:

| Part | Question |
| --- | --- |
| The NUMBER | what budget does a two-column `RESULTS BY MODEL` table actually want, given today's criteria set and today's five columns? Re-derived, not inherited (§ "The number is inherited") |
| The WORDING | does every statement of the rule describe the CURRENT structure? The rule was written for a layout that no longer exists, and stale phrasing survives elsewhere: `ONBOARDING_MODELS.md` and `REPRODUCE.md` still speak of a "matrix cell" carrying evidence, which is now the icon-only matrix and carries none |
| The SCOPE | the budget governs the per-model tables only. The at-a-glance matrix (icons), the simplest-test companion table, and the score-board have no stated budget. Confirm that is deliberate and say so in the rule, rather than leaving it as an absence a reader has to infer |

Files that change at minimum:

| File | Change |
| --- | --- |
| [`utils/check_models_md.py`](../utils/check_models_md.py) | `LIMIT` default 55 → 65, plus the docstring's opening line and its check-1 description |
| [`MODELS.md`](../../MODELS.md) | the `Cell format (the 55-word rule)` paragraph: the rule name and the number inside it |
| [`PR_REVIEW_STANDARDS.md`](../PR_REVIEW_STANDARDS.md) | § 7.1 check 1, "over 55 words of prose" |

Historical references to 55 in closed task records ([T1](t1_task_details.md) and the
change-logs) are NOT rewritten: they record what the budget was at the time.

⚠️ **The table above is the PLAN, and the run went wider.** The three files did change as
listed, but so did what the rule says about itself, five stale-wording sites the plan did
not anticipate, and three stale criteria counts folded in mid-run. See § DEVIATIONS LOG.

### Why this is a rendering rule, not a content rule

The reasoning that motivates the raise, recorded because it is the part worth keeping:

| Point | Statement |
| --- | --- |
| What the number actually governs | how much prose fits in a markdown table cell before the table stops being scannable. It is a **visualization** constraint on the rendered table, not a judgement about how much a physics claim deserves to say |
| Why the same number can serve both files | the per-model `RESULTS BY MODEL` tables are two-column, and a roadmap row's `Description` sits in a four-column table. Both are read as tables, both fail the same way when a cell turns into a paragraph |
| Why the number must be re-derived rather than nudged | see the lineage below: 55 is the budget of a table layout that no longer exists, restored by hand when the layout changed |

### The number is inherited, not derived (git-verified, 2026-07-28)

The user's account of the history is correct and the commits confirm it:

| Date | Commit | What happened to the structure | Budget |
| --- | --- | --- | --- |
| before 2026-07-25 | | ONE wide table per section, `\| Criteria \| Liquid Crystal (M5) \| HydroBoros (M7) \| EWT (M4) \| Ouroboros (M6) \| MIT (M8) \|`: every model's EVIDENCE PROSE in its own column, six columns wide | 55 |
| 2026-07-25 10:19 | `bbe27e7f` | unchanged, but the wide table was unreadable at 55 words per cell across five model columns | 55 → **45** |
| 2026-07-25 | `f80c44ef` | THE SPLIT: an icon-only at-a-glance matrix plus one two-column `\| Criteria \| Status + result summary \|` table per model | 45 → **55** |

So the current 55 was not calculated for the two-column layout. It is the pre-split
number, restored when the split removed the pressure that had forced 45. The structure
then kept moving underneath it: [M8](../../openwave/xperiments/m8_mit/) joined on
2026-07-21 as a fifth column, and [T1](t1_task_details.md) took the criteria set from 22
to 31 rows on 2026-07-28.

⚠️ This corrects a claim made when T3 was proposed, that the budget-governed tables "were
always two-column". Narrowly that is true, since the per-model tables were born
two-column at the split. But it misses the point that matters here: the NUMBER predates
the split and was inherited across it. A budget calibrated on a six-column table, walked
down to 45 and back up to 55 by hand, has no derivation behind it for the layout it now
governs. That is the actual case for T3, and it is stronger than the "one number across
two files" tidiness argument that opened this task.

⚠️ **Do NOT generalize this into a platform-wide "all table descriptions are 65 words"
rule as part of this task.** The user's instruction (2026-07-28) is to note the
observation and stop there: other table kinds may well want different budgets, and the
right time to define a general rule is when a third case actually asks for one. Two
files agreeing is a coincidence worth having, not a principle worth declaring.

### Gating

| Gate | Why |
| --- | --- |
| [PR #362](https://github.com/openwave-labs/openwave/pull/362) merged FIRST, with its over-budget cell trimmed to 55 | the raise must not land while a contribution is failing the old limit. In a public repo the git history would show the standard being widened to pass a pending PR, which spends credibility the standard needs. The trim is cheap (10 words, and the detail already lives in the linked method note) |
| User "go" | a standards change, so it is the user's call, not a maintenance sweep |

### Definition of done

| # | Item |
| --- | --- |
| 1 | The budget is RE-DERIVED for the current two-column layout and the derivation is written down, so the next structural change can test the number instead of inheriting it again |
| 2 | Every file stating a cell budget carries the same value, and a repository-wide grep finds no other number |
| 3 | `python3 dev_docs/utils/check_models_md.py` exits 0 at the new limit |
| 4 | Mutation check: a cell one word over the limit fails, a cell exactly at it passes (the limit is enforced at the new value, not merely renamed) |
| 5 | The word-count measurement is re-run and recorded, so the change rests on evidence at execution time rather than at proposal time |
| 6 | The WORDING review is done: every statement of the rule describes the current structure, and the stale "matrix cell" phrasing in [`ONBOARDING_MODELS.md`](../../ONBOARDING_MODELS.md) and [`REPRODUCE.md`](../../REPRODUCE.md) either matches the icon-only matrix or is repointed at the per-model tables |
| 7 | The rule states explicitly which tables it governs and which carry no budget (matrix, simplest-test companion, score-board), so the absence reads as a decision rather than an oversight |

### The measurement at proposal time (2026-07-28)

Taken on `main` after the T1 restructure, over all 155 per-model summary cells:

| Words of prose | Cells |
| --- | --- |
| 0-20 | 94 |
| 21-35 | 29 |
| 36-45 | 20 |
| 46-50 | 10 |
| 51-55 | 2 |
| 56-65 | 0 |

Median 17, mean 20.5, max 54. **No existing cell needs the raise**, which is the honest
statement of what this task is: an alignment for future writing, not a repair. If the
re-run at execution time shows cells crowding 55, that strengthens the case; if it still
shows a median near 17, the task is optional and should be closed as such rather than
executed for tidiness.

### Blindspot pass

| Blindspot | Mitigation |
| --- | --- |
| A wider budget quietly becomes the target, and cells drift toward 65 | the `MODELS.md` paragraph keeps its "condensed summary, never a report" framing and the "sharper sentence plus a better link" instruction, which is the part doing the real work. Only the number moves |
| Model authors read the raise as permission to move detail out of their notes and into the matrix | the change-log entry states the intent (rendering alignment) explicitly, and the linked-record principle is restated in the same edit |
| A stale 55 survives somewhere and the two disagree | the DoD requires a repository-wide grep for a stated budget, not just the three known files |

## DEVIATIONS LOG

| # | Deviation | Action taken |
| --- | --- | --- |
| 1 | **The first pass reached the wrong answer by measuring half the problem.** Reading load alone says the two caps already render the same cell, so the run initially recommended keeping 55. The user's stated reason for the raise, that the split into two-column tables bought room, is a COLUMN-GEOMETRY argument that pass never measured | Measured it (F3). The summary column takes 1.36× the width of a roadmap `Description` column, which moves the equivalent cap from 49 to 74. The honest answer is a bracket, and 65 sits inside it. Raise implemented |
| 1b | The plan's parity argument ("one number across two files") is still wrong, and is now MORE misleading, because after the raise the two numbers coincide while still counting different things | The rule text no longer says the numbers match. It says they are the same number in different units, states the 1.36× width ratio as the reason, and tells the next reader not to "fix" either to match the other |
| 2 | The plan named [`ONBOARDING_MODELS.md`](../../ONBOARDING_MODELS.md) and [`REPRODUCE.md`](../../REPRODUCE.md) as the two stale-wording sites. ONBOARDING is NOT stale ("the four icons a matrix cell can carry" is exactly right for the icon-only matrix), and four other root docs are | Left ONBOARDING alone; fixed `REPRODUCE.md` (two places), `AI_HYGIENE.md`, `TUTORIAL.md` (two places) and `CLAUDE.md` |
| 3 | [`PR_REVIEW_STANDARDS.md § 7.1`](../PR_REVIEW_STANDARDS.md) listed five linter checks; the linter has carried six since [T1](t1_task_details.md) added the simplest-test check | Added the missing row, renumbered, and added the do-not-equalize warning where a reviewer reads § 7.1 and § 7.2 back to back |
| 4 | Platform tasks had no checkpoint home: `.gitignore` covered `openwave/xperiments/*/research/checkpoints/` only | Added `dev_docs/checkpoints/` to the same rule and used it for this run |
| 5 | Stale criteria counts (21 where the set has been 31 since [T1](t1_task_details.md)) found in three model-owned docs, a T1 tail | Surfaced at review as out of scope, then folded into this task on the user's instruction. Fixed in all three (F8) |

## FINDINGS

### F1. The budget moves 55 → 65, and the number is now bracketed by two measurements

The re-derivation ran over all 155 per-model summary cells and all 88 roadmap
`Description` cells, with [`../utils/models_cell_stats.py`](../utils/models_cell_stats.py).
A cell budget is a limit on how many rendered LINES a cell occupies, and lines are
`characters / column width`, so the derivation needs both factors. They pull opposite ways:

| Factor | What it says about MODELS.md's equivalent cap | Section |
| --- | --- | --- |
| Reading load: this file's rule excludes freight the roadmap rule counts, so a cell renders larger than its number claims | **49** counted words | F2 |
| Column geometry: the two-column table gives its summary column 1.36× the width, so the same characters wrap into fewer lines | **74** counted words | F3 |

**The honest answer is the bracket [49, 74], not a point**, because neither factor is
measured to better than the crudeness of a column-width model. **65 sits inside it**, near
the midpoint, and is the value the maintainer set. Recording the bracket rather than a
single derived number is the durable part: the next structural change re-runs the script
and checks whether the standing number is still inside.

⚠️ **The first pass of this task recommended keeping 55**, having measured reading load
and not geometry. The geometry half came from the user's own reason for the raise (the
split into two-column tables bought room), and measuring it moved the answer. Recorded
because the failure mode is worth remembering: a derivation that measures one of two
opposing factors will land confidently on the wrong side.

### F2. Reading load: the two linters do not count the same thing

"One budget across two files" assumed they do:

| Excluded from the count | [`check_models_md.py`](../utils/check_models_md.py) | [`check_roadmaps.py`](../utils/check_roadmaps.py) |
| --- | --- | --- |
| Link targets | yes | yes |
| Link **labels** | yes | **no, they count** |
| The leading `<icon> [status]` tag | yes | no such element |
| `<br>→` pointer tails | yes | no such element |

Measured effect: a MODELS.md cell reads as a median **1.33×** its counted words (max
2.33×), a roadmap cell exactly **1.00×** by construction. So the two 65s are the same
number in different units, and a MODELS.md cell at 65 puts about a third more on the page
than a roadmap row at 65. **That is now the most important thing the rule says about
itself**, because the coincidence of numbers is more misleading after the raise than
before it.

### F3. Column geometry: the half that decided the task

Column widths are set per column for the whole table, roughly in proportion to
max-content, so the comparison must be per table and never per row:

| Table kind | Columns | Share of table width taken by the budgeted column |
| --- | --- | --- |
| MODELS.md per-model results (n = 5) | `Criteria` (max 39 chars) + summary | **0.910** |
| Roadmap scoped tables (n = 22) | TaskID + Title + `Description` + Gated By | **0.668** |

A **1.36× width ratio**. That is exactly the room the 2026-07-25 split bought: before it,
the summary prose sat in one of five model columns; after it, each model's prose has a
column almost to itself. At equal rendered lines that ratio carries a roadmap cell's ~413
characters up to ~563, which at the measured 7.6 characters per counted word is 74.

⚠️ Stated, not modelled precisely: this uses max-content proportional shrink-to-fit as a
stand-in for the browser's auto table layout, and ignores min-content floors. Measuring
per row instead gives only 1.06×, because a roadmap `Description` carries most of its
row's characters even where its column is narrow. The per-table number is the right one,
and the gap between the two is why the bracket is wide.

### F4. Nothing needed the room, so the risk is drift, not breakage

| Corpus | Median counted | Share of its cap | Cells in the top band |
| --- | --- | --- | --- |
| Roadmap `Description` (cap 65) | 52 | 80% | 34 of 88 sit at 56-65 |
| MODELS.md summary (old cap 55) | 17 | 31% | 5 of 155 sat at 51-55 |

**0 of 155 cells** were in the 56-65 band the raise opens, so this is an allowance for
future writing, not a repair. The roadmap corpus shows what happens when a budget becomes
a target: its rows sit at 80% of cap. MODELS.md's sit at 31%. The mitigation is the part
of the rule that did NOT change: "a condensed summary, never a report", and "the fix is a
sharper sentence plus a better link, not a longer cell". Only the number moved.

### F5. What the rule now says

[`MODELS.md`](../../MODELS.md) states the rule in **one condensed paragraph** (the cap, the
exclusions, the warning that this 65 is not the roadmap's 65, and the re-measure command),
plus **which tables carry a budget** as a four-row table so the absence reads as a
decision:

| Table | Budget |
| --- | --- |
| Per-model results tables | 65 words of prose per summary cell |
| At-a-glance matrix | none: icons, and no prose at all, which the linter enforces |
| Simplest-test companion | none needed: one short clause per row, observed maximum 12 words |
| Score-board | none: counts only |

The same warning is repeated in [`check_models_md.py`](../utils/check_models_md.py)'s
docstring and in [`PR_REVIEW_STANDARDS.md § 7.1`](../PR_REVIEW_STANDARDS.md), which is
where a reviewer comparing the two linters would otherwise draw the wrong conclusion.

### F6. The boundary is enforced at 65, not merely renamed (mutation check)

Run on a scratch copy of the repo (`MODELS.md` plus the linter at its expected depth), with
the M8 quarks cell, the one trimmed to exactly 55 in
[PR #362](https://github.com/openwave-labs/openwave/pull/362), padded to sit on the new
boundary:

```text
cell padded to 66 words  ->  limit 65 | max 66 ... L352 (Quarks, M8) over budget: 66 > 65
                             exit 1
cell padded to 65 words  ->  limit 65 | max 65 ... clean
                             exit 0
```

One word over fails, exactly at passes. The live file is untouched by the test and stays
clean at 155 cells, max 55.

### F7. Wording review: one plan target was fine, four other docs were not

| File | State | Change |
| --- | --- | --- |
| [`ONBOARDING_MODELS.md`](../../ONBOARDING_MODELS.md) | correct | none. "The four icons a matrix cell can carry" describes the icon-only matrix accurately, and the `regime` reference at line 246 is current |
| [`REPRODUCE.md`](../../REPRODUCE.md) | stale, twice | the traversal table routed readers from a "coverage-matrix cell" to "the script it links"; the matrix links nothing. Repointed at the per-model results table, and the same fix applied to the "good first reproduction" line |
| [`AI_HYGIENE.md`](../../AI_HYGIENE.md) | stale | "no cell in the coverage matrix rests on prose" now reads "no icon ... the per-model results table behind it links every claim" |
| [`TUTORIAL.md`](../../TUTORIAL.md) | stale, twice | the agent-prompt example and the reproduction code comment both assumed a matrix cell carries a script link |
| [`CLAUDE.md`](../../CLAUDE.md) | stale | the doc-map row claimed "every cell links the runnable script"; now names the two-layer structure |

### F8. Stale criteria counts in three model docs, fixed

[T1](t1_task_details.md) took the criteria set from 22 to 31 and did not sweep the docs
that state the count. Surfaced at review as out of scope, then folded in on the user's
instruction:

| File | Was | Now |
| --- | --- | --- |
| [`m8_platform_pointers.md`](../../openwave/xperiments/m8_mit/research/m8_platform_pointers.md) | "the shared 21-criteria coverage matrix" | 31-criteria |
| [`m6_roadmap.md`](../../openwave/xperiments/m6_ouroboros/research/m6_roadmap.md) | "21 criteria ... 3 ✅ / 3 ⚠️ / 3 ❌ / 12 🚧" | 31 criteria, 22 🚧, matching the live score-board |
| [`__M7_model_briefing.md`](../../openwave/xperiments/m7_hydroboros/__M7_model_briefing.md) | "the 21-cell column is PUBLISHED", twice, with entry-time icons read as current | "the M7 column"; the Implementation Status line keeps the 21-cell entry figures but dates them to M7.21 and points at the live score-board |

Not touched: historical statements in closed task records ("published the 21-cell
column", [`m7_21_models_column.md`](../../openwave/xperiments/m7_hydroboros/research/tasks/m7_21_models_column.md),
the M7 roadmap DONE rows). Those were true when written and are the record of what
happened.

### F9. The derivation is re-runnable

[`../utils/models_cell_stats.py`](../utils/models_cell_stats.py) is the artifact that
makes DoD 1 real: the next structural change to MODELS.md can test the number instead of
inheriting it a second time. It reports both corpora under three measures, the histogram
of each, the simplest-test companion table for the scope claim, the column-geometry shares,
and the equivalent-cap bracket that decided this task. `--csv` dumps every cell.

```bash
python3 dev_docs/utils/models_cell_stats.py
```

### F10. Cross-check of the load-bearing measurement

The bracket rests on measured rendered sizes. Per the adversarial-audit rule
([`AI_HYGIENE.md § 1`](../../AI_HYGIENE.md)), the reading-load end was re-measured by a
second route that imports none of the stats script: different regex order, link labels
substituted before `<br>` handling, characters counted on the raw table line.

| Cell | Visible words | Rendered characters |
| --- | --- | --- |
| MODELS.md M8 quarks, counted at the old 55 cap | 60 | 388 |
| MODELS.md M5 quarks, counted at 39 | 61 | 396 |
| M8 roadmap M8.1.1 row, counted at 64 of the 65 cap | 64 | 451 |
| M8 roadmap M8.1 row, counted at 57 | 57 | 399 |

Confirms the F2 end of the bracket: at the OLD 55 cap a MODELS.md cell already rendered
the same size as a roadmap row at 65, so the raise is a genuine widening in reading load
and is justified only by the column-width finding (F3), not by parity.

⚠️ **The geometry end (F3) has no independent cross-check.** It rests on one
proportional-shrink-to-fit model of table layout, and the per-row alternative disagrees
with it by a factor of 1.3. That is the weakest link in this task, which is why the result
is recorded as a bracket with 65 inside it rather than as a derived point value.

**No method note is produced for this task**, deliberately: [`METHOD_NOTE.md`](../METHOD_NOTE.md)
governs results reported to a model owner or an internal physics audit, and this is a
documentation standard carrying no physics claim. The re-runnable script (F9) plus this
cross-check are the equivalent audit surface.

### Definition of done

| # | Item | State |
| --- | --- | --- |
| 1 | Budget re-derived for the current layout, derivation written down | ✅ F1-F3, recorded as a bracket, re-runnable via the script (F9) |
| 2 | Every file stating a cell budget carries the same value, repo-wide grep clean | ✅ 65 in `MODELS.md` (twice), the linter default and § 7.1; no other number states a cell budget |
| 3 | `check_models_md.py` exits 0 at the new limit | ✅ clean, 155 cells, limit 65, max 55 |
| 4 | Mutation check: one word over fails, exactly at passes | ✅ F6, at the new boundary on a scratch copy |
| 5 | Measurement re-run and recorded at execution time | ✅ F1-F4 |
| 6 | Wording review done | ✅ F7, with the plan's own target list corrected |
| 7 | The rule states which tables carry no budget | ✅ F5 |
| + | (added mid-run) Stale criteria counts fixed | ✅ F8 |

## TASK REVIEW (2026-07-28)

`Task Duration: 00:40 (from 22:13 to 22:53 EDT)`
`Usage Cap Triggered: NO`

### Results

| Item | Verdict |
| --- | --- |
| The budget | ✅ **55 → 65**, on column geometry, not on the parity argument the task was filed on |
| The derivation | ✅ a measured bracket [49, 74] with 65 inside it, re-runnable, replacing an inherited number |
| The rule's self-description | ✅ the exclusions, the not-the-same-unit warning and the width ratio now stated in all three places a reader meets the rule |
| Wording review | ✅ 5 stale sites in 4 root docs repointed at the per-model tables; the plan's own target list was itself wrong (F7) |
| Stale criteria counts | ✅ 3 model docs corrected, 21 → 31 (F8) |
| Page prose | ✅ the three long `MODELS.md` paragraphs cut to half length at the user's call: the criteria-versus-models tables are the page, the prose around them is not |
| Enforcement | ✅ mutation-tested at the new boundary (66 fails, 65 passes); all three checkers clean |

### Issues and blockers

None blocking. One stated weakness: the geometry end of the bracket rests on a single
proportional-shrink-to-fit model of table layout and has no independent cross-check, and
the per-row alternative disagrees by 1.3× (F10). It is the softest number in the task and
it is the one carrying the raise, which is why the result is a bracket rather than a point.

### Deviations from plan

Five, logged as they happened (§ DEVIATIONS LOG). The load-bearing one: the first pass
measured reading load only, concluded "keep 55", and was wrong because it never measured
column width. The user's own reason for the raise named the missing variable.

### Action needed

None outstanding. Two commits by design, `4b457177` carrying the derivation that kept 55
and its successor carrying the raise, so the history shows a standard measured first and
changed with a stated reason rather than bent to fit.

### Findings

The `MODELS.md` cell budget moves 55 → 65 on a column-geometry derivation: the 2026-07-25
split into two-column per-model tables gave the summary column 0.910 of its table's width
against a roadmap `Description` column's 0.668, a 1.36× ratio that puts the equivalent cap
in the bracket [49, 74]. The parity argument that motivated the task was wrong both before
and after the change, because the two linters exclude different things, so the newly
coincident numbers are now explicitly flagged as the same number in different units. No
existing cell needed the room (0 of 155 above 55), making this an allowance for future
writing whose real risk is drift, which is why the "condensed summary, never a report"
framing was left standing while everything around it was cut in half.

### Research documents created / updated

| Document | What changed |
| --- | --- |
| [`t3_task_details.md`](t3_task_details.md) | this record: F1-F10, deviations, DoD, review |
| [`../utils/models_cell_stats.py`](../utils/models_cell_stats.py) | **new**: the re-runnable derivation, both corpora, three measures, column geometry, the bracket |
| [`../../MODELS.md`](../../MODELS.md) | the 65-word rule, what it counts, which tables carry a budget; three paragraphs halved |
| [`../PR_REVIEW_STANDARDS.md`](../PR_REVIEW_STANDARDS.md) | check 1 → 65, the missing check 5, the do-not-align warning |
| [`../utils/check_models_md.py`](../utils/check_models_md.py) | `LIMIT` 55 → 65, docstring |
| [`../../REPRODUCE.md`](../../REPRODUCE.md) · [`../../AI_HYGIENE.md`](../../AI_HYGIENE.md) · [`../../TUTORIAL.md`](../../TUTORIAL.md) · [`../../CLAUDE.md`](../../CLAUDE.md) | matrix-cell routes repointed at the per-model tables |
| [`m6_roadmap.md`](../../openwave/xperiments/m6_ouroboros/research/m6_roadmap.md) · [`__M7_model_briefing.md`](../../openwave/xperiments/m7_hydroboros/__M7_model_briefing.md) · [`m8_platform_pointers.md`](../../openwave/xperiments/m8_mit/research/m8_platform_pointers.md) | criteria counts 21 → 31 |
| [`../platform_roadmap.md`](../platform_roadmap.md) | T3 row closed into DONE |
