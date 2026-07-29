# T3: MODELS.md cell budget, raise 55 to 65

> Roadmap row: [`../platform_roadmap.md`](../platform_roadmap.md) T3 (BACKLOG, ahead of
> [T2](t2_task_details.md)). Proposed by the user 2026-07-28 while
> [PR #362](https://github.com/openwave-labs/openwave/pull/362) was open with one cell
> over budget. Deliberately NOT executed that day; see § Gating.

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

(not started)

## FINDINGS

(not started)
