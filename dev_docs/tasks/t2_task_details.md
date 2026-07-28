# T2: MODELS.md normalized model score: a fair single number per column

**Status**: 🔶 DEFERRED 2026-07-28, the same day it was proposed (user decision). Parked at the Backlog tail with an explicit re-open trigger: **a critical mass of tested rows**. Platform-wide, docs-plus-linter, no physics runs. The intent was to add a derived score to the [`MODELS.md`](../../MODELS.md) SCORE-BOARD so a reader could compare columns without counting icons; the analysis below is kept complete so the revisit starts from here rather than from scratch. When it does run it goes AFTER [T1](t1_task_details.md), which changes the criteria set and therefore every count the score is computed from.

## DECISION: deferred (2026-07-28, user)

| Point | Content |
| --- | --- |
| What stands today | The score-board **already does the two jobs that matter**: it separates covered from not-covered, and it gives the column ordering. Nothing is broken and nothing is being replaced |
| Why not now | A scalar leaderboard at today's coverage buys more discussion and author complaint than information. Four of five columns sit under 50% attempted, so any single number is mostly reporting how little has been run |
| Re-open trigger | A critical mass of tested rows across the columns (user call). At that point the small-sample defect below stops dominating, and the two cures listed become testable rather than theoretical |
| What is preserved | The formula analysis, the three defects, the sensitivity band, the regime cut, the linter design, and the ordering pre-commitment. All of it below, none of it in `MODELS.md` |

## The proposal as received

The user's spreadsheet, verbatim in structure (the two bottom rows are the new ones; the four icon rows and `Total criteria` already exist in MODELS.md):

| MODEL SCORE-BOARD | WEIGHT | LC (M5) | HydroBoros (M7) | EWT (M4) | Ouroboros (M6) | MIT (M8) |
| --- | --- | --- | --- | --- | --- | --- |
| ✅ validated in-platform | 1 | 9 | 0 | 0 | 3 | 0 |
| ⚠️ partial / with caveats | 0.5 | 7 | 9 | 8 | 3 | 1 |
| ❌ honest negative | -1 | 1 | 0 | 3 | 3 | 0 |
| 🚧 planned / not tested | 0 | 5 | 13 | 11 | 13 | 21 |
| **Total criteria** | | **22** | **22** | **22** | **22** | **22** |
| weighted total | `SUMPRODUCT($B2:$B5,C2:C5)` | 11.5 | 4.5 | 1.0 | 1.5 | 0.5 |
| NORMALIZED SCORE | `(C7/C6+1)/0.2` | 7.6 | 6.0 | 5.2 | 5.3 | 5.1 |

(The SUMPRODUCT row is a sum of products, so `weighted total` is the accurate label; the spreadsheet called it `weighted product`.)

**What the formula is, decoded**: `credit = Σ wᵢnᵢ / N` is the mean credit per criterion, which lives in `[-1, +1]`; `(credit + 1) / 0.2` is the linear rescale of that interval onto `[0, 10]`. So all-✅ = 10.0, all-❌ = 0.0, and **all-🚧 = 5.0**. Reproduced independently and it matches the spreadsheet to the printed digit for all five columns.

## The two problems with it (the user's own question, sharpened)

| # | Problem | Evidence in the current numbers |
| --- | --- | --- |
| 1 | **Untested scores as average, not as unknown.** 🚧 has weight 0 but the denominator is the FULL criteria set, so every untested row drags a column toward 5.0, which reads on a 0-10 scale as "middling model" rather than "not yet examined" | M8 has ONE attempted criterion out of 22 and scores **5.1**, one tenth of a point below M4, which has attempted 11 and carries 3 documented negatives. The score is measuring silence |
| 2 | **Honest negatives score below silence.** ❌ = -1 is strictly worse than 🚧 = 0, so running a test and reporting failure costs a column more than never running it | Directly contradicts the file's own stated ethos ([`MODELS.md § Summary Count`](../../MODELS.md#summary-count)): "A ❌ is a result, not an embarrassment ... documented negatives (with the scripts that produced them) are part of the platform's value." A published score that punishes them makes that sentence untrue |

Both problems have one root: **coverage and performance are being collapsed into a single number**, and once collapsed, no choice of weights can separate "has not been tested" from "was tested and scored in the middle".

## The fix for those two, and the third defect it brings with it

| Quantity | Formula | Reads as |
| --- | --- | --- |
| **Attempted coverage** | `A / N` where `A = ✅ + ⚠️ + ❌`, `N` = total criteria | how much of the arena this column has actually been put through |
| **Score on attempted** | `10 × (✅ + 0.5·⚠️) / A`, undefined when `A = 0` | of what it attempted, how much it earned |

The denominator change (`A` instead of `N`) fixes problem 1. Dropping the negative weight fixes problem 2, **and it becomes the coherent choice once untried rows leave the denominator**: with `A` as denominator, ❌ already costs the column by sitting in the denominator and contributing nothing to the numerator. Keeping `w = -1` on top of that penalizes the same failure twice, which is what makes reporting a negative worse than staying quiet.

Applied to today's counts:

| Model | attempted | coverage | score on attempted | current single number |
| --- | --- | --- | --- | --- |
| M5 | 17 | 77% | **7.4** | 7.6 |
| M7 | 9 | 41% | **5.0** | 6.0 |
| M4 | 11 | 50% | **3.6** | 5.2 |
| M6 | 9 | 41% | **5.0** | 5.3 |
| M8 | 1 | 5% | **5.0** | 5.1 |

M8 now reads `5.0 on 5% coverage`, which is the honest statement (one partial result, nothing else attempted) instead of a number that looks like a mid-table finish. M4 drops to 3.6, which is what its evidence says: half the arena attempted, mostly partials, three failures.

**Publication rule that would go with it**: the score NEVER appears without its coverage companion, in the same cell or the adjacent row, in every place it is quoted. A number without its sample size is the thing that gets screenshotted.

### Defect 3 (user catch, 2026-07-28): the fix is unstable at small samples

Dividing by attempted rows cures the silence problem by making the denominator small, and a small denominator is its own failure mode. **If M8's single gravity row were upgraded from ⚠️ to ✅, M8 would print 10.0 out of 10, top of the board, on 5% coverage.** One validated criterion out of 22 would outrank a column with nine validations. The publication rule (never without coverage) mitigates the reading but does not fix the number, and the number is what travels.

So the two-number form is a strict improvement on the original and still not publishable as-is. The candidate cures, both deferred to the revisit:

| Cure | Mechanism | Cost |
| --- | --- | --- |
| **Minimum-coverage gate** | no score printed below a floor (say 25% attempted); the cell reads `n/a, X% coverage` until the column clears it | one threshold, transparent and arguable, but a hard cliff: a column crossing the floor appears from nowhere at a high score |
| **Shrinkage toward a prior** (option B below, promoted) | `10 × (✅ + 0.5⚠️ + k·p) / (A + k)` with prior `p` and pseudo-count `k`; few-sample columns are pulled toward the prior and earn their way out | two knobs instead of zero, and the score stops being hand-checkable from the visible counts. Worked example at `p = 0.5, k = 5`: M8-with-one-✅ prints 5.8 instead of 10.0, M5 prints 6.8 |

Neither is worth adopting while most columns sit under half coverage, which is the deferral above.

## Sensitivity: the middle of the table is not stable, and the file must say so

The ⚠️ half-credit is a judgment, not a measurement, so the score was recomputed at `w = 0.3` and `w = 0.7`:

| Model | w = 0.3 | w = 0.5 | w = 0.7 |
| --- | --- | --- | --- |
| M5 | 6.5 | 7.4 | 8.2 |
| M7 | 3.0 | 5.0 | 7.0 |
| M4 | 2.2 | 3.6 | 5.1 |
| M6 | 4.3 | 5.0 | 5.7 |
| M8 | 3.0 | 5.0 | 7.0 |

Only the ends are stable: M5 first and M4 last at every weight. **M6 versus M7 flips**, because M6's credit comes from 3 validated rows and M7's comes entirely from partials, so M6 wins at low partial-credit and loses at high. At `w = 0.5` they are exactly tied, which is a knife edge, not a result. Deliverable: publish the score with a one-line note that ranks inside the tied band are weight-dependent, or publish the band itself. Do not present a 0.1-point gap as an ordering.

## Alternatives considered

| Option | What it does | Verdict |
| --- | --- | --- |
| A: two numbers (above) | separate coverage from performance | **The best of these, and still not enough**: fixes defects 1 and 2, opens defect 3. Would need the coverage gate bolted on |
| B: shrinkage / Bayesian prior | one number, `10 × (✅ + 0.5⚠️ + k·p)/(A + k)`, few-sample columns pulled toward the prior (the IMDb-top-250 / Wilson-interval trick) | **PROMOTED to leading candidate at revisit** (it is the principled cure for defect 3). Cost unchanged: it hides the sample size inside the number instead of showing it, and `k` is a knob a column author can argue with |
| C: keep the single number, drop ❌ to 0 | fixes problem 2 only | Rejected: untested still scores 5.0 |
| D: difficulty weights per criterion | harder rows worth more | Rejected: it puts maintainer judgment into a tally that is currently mechanical and auditable. The mechanical substitute is the regime cut below |
| E: per-regime cut (secondary, optional) | report coverage + score separately over `static` / `dynamic` / `both` rows, using the existing `regime` column | **Worth doing** if it fits: it is fully mechanical and it explains the shape of a column. M5 today reads static 8.3 on 55%, dynamic 6.9 on 100%; M8 reads 0% static coverage, 0% dynamic, one attempted `both` row. That is a fair description of a top-down column that has no engine yet, which a single blended number erases |

## Scope (if and when it runs)

| Piece | Content |
| --- | --- |
| Small-sample handling | Defect 3 must be closed BEFORE anything is published: a minimum-coverage gate, shrinkage, or a demonstration that coverage has risen far enough that neither is needed. No score ships while a one-row column can top the board |
| MODELS.md rows | Two new SCORE-BOARD rows (`Attempted coverage` + `Score on attempted (0-10)`), a WEIGHT column on the icon rows so the formula is visible in the file, and a short paragraph under the table stating the formula, the ⚠️-weight sensitivity, and what the score is not |
| Linter extension | [`dev_docs/utils/check_models_md.py`](../utils/check_models_md.py) currently rejects any score-board row that is neither an icon nor a total (L151-L153). Teach it the two derived rows, RECOMPUTE them from the tallies it already builds, and fail with the correct values printed on mismatch. The numbers must never be hand-maintained: T1 alone will move every one of them |
| Criteria-set stamp | The score is not comparable across criteria-set revisions (N changes when rows split). The row label or the caption carries the criteria-set date, so an old screenshot cannot be read as current |
| Column ordering decision | State explicitly whether the score governs column order. Today the order is coverage-flavored (✅+⚠️ count, ties to more ✅ then fewer ❌), and the score DISAGREES with it: M4 is 3rd by the existing rule and last by score. Recommendation: **the existing order stands, the score does not reorder columns**, and the reason is written down so nobody later "fixes" the inconsistency by demoting a column for reporting negatives |
| Regime cut | Option E above, if it fits the table width; otherwise recorded as a follow-on |

## What this task is not

No physics runs, no new evidence, no icon changes, no cell rewrites, and **nothing in `MODELS.md` changes on the strength of this document**: the deferral means the file stays exactly as it is. When the task does run it adds a derived view of counts that already exist, and even then, if the score and a reader's intuition disagree, the per-model tables are the record and the score is the summary, not the other way round.

## Blindspots

| Risk | Guard |
| --- | --- |
| **Arena neutrality**: the file states the platform has no house model and no stake in which column wins. A published leaderboard number where the maintainer-developed column leads by 2+ points invites exactly the opposite reading | The formula is mechanical, lives in the linter, and is computed from cells every column author can contest; coverage is published beside it; the "what this is not" line sits under the table. Flag to the user before it lands |
| Column authors first see their model scored when it is already public | The columns belong to their authors. Route the formula to them before merge, as a proposal, not as a fait accompli (user call on channel and timing) |
| Hand-typed numbers drifting from the tallies | The linter recomputes them; mismatch is an error, not a warning |
| The score being quoted without coverage | Publication rule stated above; the two rows sit adjacent so a crop of one is visibly missing the other |
| Reading a 0.1 gap as an ordering | The sensitivity band is published with the score |
| N changing under T1 mid-flight | This task is gated on T1 closing; the criteria-set stamp handles later revisions |
| **Author friction exceeding the information gained**: a scalar that ranks other people's models is a standing invitation to dispute the weights rather than the physics, and at today's coverage the number mostly reports how little has been run | The 2026-07-28 deferral. The score-board's existing covered/not-covered split and column ordering already carry the useful content without the argument |
| A one-row column topping the board | Defect 3 above; the coverage gate or shrinkage must be settled before publication, not after |

**Gated by**: 🔶 DEFERRED. Re-open on a critical mass of tested rows (user call), then [T1](t1_task_details.md) (the criteria set must be final first) + user "go" (docs-only, can interleave).

## Cross-links

| Doc | Why |
| --- | --- |
| [`MODELS.md`](../../MODELS.md) | The target doc (SCORE-BOARD + Summary Count) |
| [T1](t1_task_details.md) | The gate: it changes the criteria set and every count feeding the score |
| [`dev_docs/utils/check_models_md.py`](../utils/check_models_md.py) | The linter that must own the computation |
| [`../m5_roadmap.md`](../../openwave/xperiments/m5_liquid_crystal/research/m5_roadmap.md) | Roadmap row + the change-log entry |
