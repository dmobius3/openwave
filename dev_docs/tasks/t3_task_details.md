# T3: MODELS.md cell budget, re-derive the rule

> Roadmap row: [`../platform_roadmap.md`](../platform_roadmap.md) T3. Proposed by the user
> 2026-07-28 while [PR #362](https://github.com/openwave-labs/openwave/pull/362) was open
> with one cell over budget, and run the same day once both gates cleared (§ Gating).
>
> **Outcome: the budget stays 55.** The expected raise to 65 rested on a premise the
> measurement refutes, and the number is now derived rather than inherited (§ FINDINGS F1).

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

⚠️ **The table above is the PLAN, and the run did not follow it.** The re-derivation kept
the budget at 55, so no `LIMIT` moved and no number changed; what changed instead is what
the rule says about itself, plus five stale-wording sites the plan did not anticipate.
See § DEVIATIONS LOG entry 1 and § FINDINGS F1.

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
| 1 | **The re-derivation refutes the task's own expected answer.** The plan expected 55 → 65 for parity with the roadmap `Description` cap. The measurement says the two caps already produce the same rendered cell, so equalizing the numbers would widen this one by about a third | Kept 55, wrote the derivation down, and rewrote the rule so the next reader cannot repeat the inference. Reported as the headline result, not buried |
| 2 | The plan named [`ONBOARDING_MODELS.md`](../../ONBOARDING_MODELS.md) and [`REPRODUCE.md`](../../REPRODUCE.md) as the two stale-wording sites. ONBOARDING is NOT stale ("the four icons a matrix cell can carry" is exactly right for the icon-only matrix), and four other root docs are | Left ONBOARDING alone; fixed `REPRODUCE.md` (two places), `AI_HYGIENE.md`, `TUTORIAL.md` (two places) and `CLAUDE.md` |
| 3 | [`PR_REVIEW_STANDARDS.md § 7.1`](../PR_REVIEW_STANDARDS.md) listed five linter checks; the linter has carried six since [T1](t1_task_details.md) added the simplest-test check | Added the missing row, renumbered, and added the do-not-equalize warning where a reviewer reads § 7.1 and § 7.2 back to back |
| 4 | Platform tasks had no checkpoint home: `.gitignore` covered `openwave/xperiments/*/research/checkpoints/` only | Added `dev_docs/checkpoints/` to the same rule and used it for this run |
| 5 | Stale criteria counts (21 where the set has been 31 since T1) found in three model-owned docs | NOT fixed. Out of scope for a standards task and they are author-facing surfaces; listed at review for a follow-up decision (§ FINDINGS F7) |

## FINDINGS

### F1. The number stays 55, and it is now derived

The re-derivation was run over all 155 per-model summary cells and all 88 roadmap
`Description` cells, with [`../utils/models_cell_stats.py`](../utils/models_cell_stats.py):

| Question | Answer |
| --- | --- |
| Does any cell need more than 55? | No. Median 17, mean 21, max 55, and **0 of 155 cells** fall in the 56-65 band the raise would open |
| Do the two caps already agree? | Yes, on the only thing a reader sees. MODELS.md cells at 51-55 counted words render at a median 404 characters; roadmap cells at 61-65 render at a median 413. A 2% difference |
| So why are the numbers different? | Because the counting rules are. See F2 |

The budget is therefore **confirmed at 55**, not replaced. Its history is still an
inheritance (§ "The number is inherited"), but the number it inherited lands within 2% of
the roadmap's on the rendered measure, so the honest verdict is that it was right by
luck and is now right on evidence.

### F2. The premise of the raise was arithmetically false

"One budget across two files" assumed the two linters count the same thing. They do not:

| Excluded from the count | [`check_models_md.py`](../utils/check_models_md.py) | [`check_roadmaps.py`](../utils/check_roadmaps.py) |
| --- | --- | --- |
| Link targets | yes | yes |
| Link **labels** | yes | **no, they count** |
| The leading `<icon> [status]` tag | yes | no such element |
| `<br>→` pointer tails | yes | no such element |

Measured effect: a MODELS.md cell reads as a median **1.33×** its counted words (max
2.33×), a roadmap cell exactly **1.00×** by construction. The worst two cells in
MODELS.md today count 45 and 51 words and read as **68**, already past the roadmap's
65-word cap while sitting comfortably inside a 55-word budget. Raising to 65 would have
permitted roughly 87 read-words, a third more than the roadmap allows, while the change
was being justified as making the two equal.

### F3. The roadmap corpus shows a budget becoming a target

The blindspot the plan listed as a risk is already visible in the other file:

| Corpus | Median counted | Share of its cap | Cells in the top band |
| --- | --- | --- | --- |
| Roadmap `Description` (cap 65) | 52 | 80% | 34 of 88 sit at 56-65 |
| MODELS.md summary (cap 55) | 17 | 31% | 5 of 155 sit at 51-55 |

MODELS.md cells are not pressing against their budget; roadmap rows are pressing against
theirs. That is an argument against moving the MODELS.md number at all, in either
direction, and it is the reason the rule keeps its "condensed summary, never a report"
framing unchanged.

### F4. What the rule now says

[`MODELS.md`](../../MODELS.md) gained two short blocks under the existing cell-format
paragraph: **what the 55 counts** (the exclusions, why they are deliberate, and the
explicit warning that the number is not comparable to the roadmap's 65), and **which
tables carry a budget**, as a four-row table so the absence reads as a decision:

| Table | Budget |
| --- | --- |
| Per-model results tables | 55 words of prose per summary cell |
| At-a-glance matrix | none: icons, and no prose at all, which the linter enforces |
| Simplest-test companion | none needed: one short clause per row, observed maximum 12 words |
| Score-board | none: counts only |

The same warning is repeated in [`check_models_md.py`](../utils/check_models_md.py)'s
docstring and in [`PR_REVIEW_STANDARDS.md § 7.1`](../PR_REVIEW_STANDARDS.md), which is
where a reviewer comparing the two linters would otherwise draw the wrong conclusion.

### F5. The boundary is enforced at the stated value (mutation check)

Not a rename: the limit bites at exactly 55.

```text
python3 dev_docs/utils/check_models_md.py       ->  clean, exit 0   (max 55, median 17)
python3 dev_docs/utils/check_models_md.py 54    ->  1 violation, exit 1
                                                   L339 (Quarks, M8) over budget: 55 > 54
```

The cell that trips at 54 is the M8 quarks cell trimmed to exactly 55 during the
[PR #362](https://github.com/openwave-labs/openwave/pull/362) maintainer edit, so the
mutation runs against a real cell sitting on the boundary rather than a synthetic one.

### F6. Wording review: one plan target was fine, four other docs were not

| File | State | Change |
| --- | --- | --- |
| [`ONBOARDING_MODELS.md`](../../ONBOARDING_MODELS.md) | correct | none. "The four icons a matrix cell can carry" describes the icon-only matrix accurately, and the `regime` reference at line 246 is current |
| [`REPRODUCE.md`](../../REPRODUCE.md) | stale, twice | the traversal table routed readers from a "coverage-matrix cell" to "the script it links"; the matrix links nothing. Repointed at the per-model results table, and the same fix applied to the "good first reproduction" line |
| [`AI_HYGIENE.md`](../../AI_HYGIENE.md) | stale | "no cell in the coverage matrix rests on prose" now reads "no icon ... the per-model results table behind it links every claim" |
| [`TUTORIAL.md`](../../TUTORIAL.md) | stale, twice | the agent-prompt example and the reproduction code comment both assumed a matrix cell carries a script link |
| [`CLAUDE.md`](../../CLAUDE.md) | stale | the doc-map row claimed "every cell links the runnable script"; now names the two-layer structure |

### F7. Found, not fixed: stale criteria counts in three model docs

[T1](t1_task_details.md) took the criteria set from 22 to 31 and did not sweep the docs
that state the count. Present tense and wrong today:

| File | Line | Says | Should say |
| --- | --- | --- | --- |
| [`m8_platform_pointers.md`](../../openwave/xperiments/m8_mit/research/m8_platform_pointers.md) | 20 | "the shared 21-criteria coverage matrix" | 31 |
| [`m6_roadmap.md`](../../openwave/xperiments/m6_ouroboros/research/m6_roadmap.md) | 7 | "21 criteria ... 3 ✅ / 3 ⚠️ / 3 ❌ / 12 🚧" | 31 criteria, 22 🚧 |
| [`__M7_model_briefing.md`](../../openwave/xperiments/m7_hydroboros/__M7_model_briefing.md) | 8 | "the 21-cell column" | 31 |

Left alone deliberately: a standards task should not quietly edit three model authors'
pages, and the fix belongs with whoever sweeps T1's tail. Historical statements ("published
the 21-cell column", in closed task records) are correct as written and are not touched.

### F8. The derivation is re-runnable

[`../utils/models_cell_stats.py`](../utils/models_cell_stats.py) is the artifact that
makes DoD 1 real: the next structural change to MODELS.md can test the number instead of
inheriting it a second time. It reports both corpora under three measures, the histogram
of each, the simplest-test companion table for the scope claim, and the head-to-head that
decided this task. `--csv` dumps every cell.

```bash
python3 dev_docs/utils/models_cell_stats.py
```

### F9. Cross-check of the load-bearing claim

Everything above rests on one measurement: that a MODELS.md cell at 55 and a roadmap row
at 65 render the same size. Per the adversarial-audit rule
([`AI_HYGIENE.md § 1`](../../AI_HYGIENE.md)), it was re-measured by a second route that
imports none of the stats script: different regex order, link labels substituted before
`<br>` handling, characters counted on the raw table line.

| Cell | Visible words | Rendered characters |
| --- | --- | --- |
| MODELS.md M8 quarks, counted at exactly the 55 cap | 60 | 388 |
| MODELS.md M5 quarks, counted at 39 | 61 | 396 |
| M8 roadmap M8.1.1 row, counted at 64 of the 65 cap | 64 | 451 |
| M8 roadmap M8.1 row, counted at 57 | 57 | 399 |

Same size class, with the MODELS.md cell at its cap landing slightly BELOW the roadmap
cell at its cap. The claim survives, and the direction of the residual argues against a
raise rather than for one.

**No method note is produced for this task**, deliberately: [`METHOD_NOTE.md`](../METHOD_NOTE.md)
governs results reported to a model owner or an internal physics audit, and this is a
documentation standard carrying no physics claim. The re-runnable script (F8) plus this
cross-check are the equivalent audit surface.

### Definition of done

| # | Item | State |
| --- | --- | --- |
| 1 | Budget re-derived for the current layout, derivation written down | ✅ F1, F2, and the script (F8) |
| 2 | Every file stating a cell budget carries the same value, repo-wide grep clean | ✅ 55 in `MODELS.md`, the linter default and § 7.1; no other number states a cell budget |
| 3 | `check_models_md.py` exits 0 | ✅ clean, 155 cells |
| 4 | Mutation check at the boundary | ✅ F5 |
| 5 | Measurement re-run and recorded at execution time | ✅ F1, F3 |
| 6 | Wording review done | ✅ F6, with the plan's own target list corrected |
| 7 | The rule states which tables carry no budget | ✅ F4 |
