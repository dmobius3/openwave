# PLATFORM ROADMAP

> The task roadmap for **platform-wide** work: the artifacts every model column shares, and the tooling that keeps them honest. Live work: [IN PROGRESS](#in-progress) → [BACKLOG](#backlog); the record: [DONE](#done). Reading rules (what belongs here, the T-ID scheme, task-doc anatomy): [§ CONVENTIONS](#conventions). Dated decisions and migrations: [§ CHANGE-LOG](#change-log).
>
> Model-specific physics does NOT live here. Each model keeps its own roadmap ([M5](../openwave/xperiments/m5_liquid_crystal/research/m5_roadmap.md), [M8](../openwave/xperiments/m8_mit/research/m8_roadmap.md), and so on); this file is for the work that would otherwise be filed under whichever model happened to notice it.

---

## IN PROGRESS

> The single currently-running platform task (one at a time; empty = between tasks).

| TaskID | Title | Description | Gated By |
| --- | --- | --- | --- |
| [T1](tasks/t1_task_details.md) | MODELS.md criteria housekeeping: split the coarse rows, name the simplest passing test per row | Author-proposed 2026-07-27 ([`tasks/m5_22_convo.md`](../openwave/xperiments/m5_liquid_crystal/research/tasks/m5_22_convo.md)), staged same day (user decision), migrated to this roadmap as T1 on 2026-07-28 (was M5.29). Platform-wide, docs-only: split the bundled criteria (Weak force → muon decay \| beta decay; Baryons → mass ordering + core/shell \| exact masses \| beta decay; sweep the rest), add a per-row simplest-passing-test entry (Larmor precession = the magnetic-force test; ortho/para-positronium the harder follow-on), add the missing Lorentz covariance row, re-derive the score-board counts + every model's detail table, checker green. No new runs, no icon upgrades; bundled icons may honestly split into better + worse sub-rows. Input #2 (2026-07-28, [`tasks/t1_convo.md`](tasks/t1_convo.md)): the author's priority-ranked properties-to-test slide (group post, [`MODELS.md`](../MODELS.md) embedded in the public talk deck) supplies most per-row tests + priority tiers + three candidate new rows (running coupling; deuteron binding + quadrupole moment; larger nuclei / halos); transcription + row mapping in the task doc. **Run-order PROMOTED ahead of [M5.22](../openwave/xperiments/m5_liquid_crystal/research/tasks/m5_22_task_details.md) (user decision 2026-07-28)**: no dependency on census results (missing-evidence rows enter 🚧), and freezing the criteria set first means M5.22 files into stable rows and [T2](tasks/t2_task_details.md) reads a stable N. | user "go" (docs-only; promoted ahead of M5.22, awaiting go) |

## BACKLOG

> Future platform tasks. **Row order = the run sequence.** Deferred rows carry their re-open trigger in `Gated By`.

| TaskID | Title | Description | Gated By |
| --- | --- | --- | --- |
| [T2](tasks/t2_task_details.md) | MODELS.md normalized score: a fair single number per column | 🔶 DEFERRED 2026-07-28, same day as the proposal (user decision), migrated to this roadmap the same day (was M5.30): **the score-board as it stands works**, it separates covered from not-covered and gives the column ordering, so a scalar score waits for a critical mass of tested rows; before that it buys more author friction than information. Docs-plus-linter, no runs. Full analysis kept in the task doc: the proposed formula (weights 1 / 0.5 / -1 / 0 over all criteria) reproduces exactly but scores silence as average (M8 = 5.1 on ONE attempted criterion) and puts ❌ below 🚧; the two-number fix (coverage + score-on-attempted) repairs both and adds a small-sample defect (user catch: M8 would read **10.0** if its single gravity row validated), cured at revisit by a coverage gate or shrinkage. | critical mass of tested rows (user call) + [T1](tasks/t1_task_details.md) final criteria set + user "go" |

## CONVENTIONS

**WHAT BELONGS HERE.** The test is ownership, not subject matter: if the artifact is shared by every model column, or is the tooling that guards one, it is platform work and it files here. If it is evidence about one model's physics, it files on that model's roadmap even when it touches a shared file.

| Platform (this roadmap) | Model roadmap |
| --- | --- |
| [`MODELS.md`](../MODELS.md) structure: the criteria set, the score-board, the cell rules | a single model's icons and evidence cells |
| [`ONBOARDING_MODELS.md`](../ONBOARDING_MODELS.md), [`REPRODUCE.md`](../REPRODUCE.md), [`AI_HYGIENE.md`](../AI_HYGIENE.md), [`dev_docs/`](.) standards | a model's own findings, notes, question tracker |
| [`dev_docs/utils/`](utils) checkers and procedural tooling | a model's research scripts under its `research/scripts/` |
| Cross-model conventions (status vocabulary, regime column, claim labels) | a model's task ladder and gates |

**T-ID SCHEME.** Platform tasks are numbered `T<n>`, assigned in order of creation and never reused. They are deliberately outside the `M<model>.<n>` space so a platform task cannot be mistaken for one model's work. A task migrated in from a model roadmap keeps a `(was M5.29)` note in its row and its change-log entry; the old ID is not reused.

**TASK DOCS.** One doc per task at [`tasks/t<n>_task_details.md`](tasks), same anatomy as the model roadmaps use (PLANNING → DEVIATIONS LOG → FINDINGS, blindspots table, `Gated by` line, cross-links table). Conversation records for a task, when the input came from a group thread or an author exchange, sit beside it as `t<n>_convo.md`.

**DONE ORDER.** Newly closed tasks are APPENDED at the END of the DONE list, so it reads in order of completion.

## DONE

> Closed platform tasks, in order of completion.

| TaskID | Title | Result | Closed |
| --- | --- | --- | --- |
| | | (none yet) | |

## CHANGE-LOG

**PLATFORM ROADMAP CREATED, TWO TASKS MIGRATED IN (2026-07-28, user decision).** MODELS.md work was being tracked on the M5 roadmap because M5 is where the proposals happened to arrive, which mis-filed it: the criteria set and the score-board belong to every column, not to the liquid-crystal one. This roadmap now owns platform-wide work, with its own `T<n>` ID space and its own [`tasks/`](tasks) folder. Migrated: **M5.29 → [T1](tasks/t1_task_details.md)** (criteria housekeeping, with its convo record [`tasks/t1_convo.md`](tasks/t1_convo.md)) and **M5.30 → [T2](tasks/t2_task_details.md)** (the normalized score, deferred). Same move: the MODELS.md linter went from `dev_docs/check_models_md.py` to [`dev_docs/utils/check_models_md.py`](utils/check_models_md.py), a procedural test rather than a standards document, with every reference updated ([`MODELS.md`](../MODELS.md), [`PR_REVIEW_STANDARDS.md`](PR_REVIEW_STANDARDS.md) § 7.1 and the review-command block, the two task docs) and its `MODELS.md` path resolution fixed for the extra directory level. The [M5 roadmap](../openwave/xperiments/m5_liquid_crystal/research/m5_roadmap.md) keeps its 2026-07-27 and 2026-07-28 change-log entries as the provenance record of how these two tasks arrived.

**T1 STAGED, THEN EXPANDED BY THE PROPERTIES-TO-TEST SLIDE (2026-07-27 afternoon → 2026-07-28 morning, as M5.29).** The author's 2026-07-27 reply proposed the criteria restructure (split the coarse rows, per-row simplest passing test, the missing Lorentz covariance row, Larmor precession as the magnetic-force test) and it was staged the same day. On 2026-07-28 a new group thread ([`tasks/t1_convo.md`](tasks/t1_convo.md)) added a "Properties to test for SM + gravity?" slide to the public talk deck: a 16-line test list with explicit priority shading set beside a screenshot of the [`MODELS.md`](../MODELS.md) matrix itself, plus an open question to the group. That converts the test column from a design exercise into a transcription-and-fit exercise, confirms every 2026-07-27 ask, and adds three candidate new rows (running coupling; deuteron binding + quadrupole moment; larger nuclei / halos). It also raises the stakes: the deck's opening slide names OpenWave and carries the matrix icons into a colloquium context, so the re-derivation has a public external audience. Any reply to the group question is user-gated and follows the T1 output. Same day, T1 was promoted ahead of [M5.22](../openwave/xperiments/m5_liquid_crystal/research/tasks/m5_22_task_details.md) in the run order (docs-only, no dependency on census results, and freezing the criteria set first gives the census stable rows to file into).

**T2 PROPOSED AND DEFERRED THE SAME DAY (2026-07-28, user proposal + user decision, as M5.30).** A per-column 0-10 score was proposed for the score-board (icon weights 1 / 0.5 / -1 / 0, denominator = the full criteria set, linear rescale onto 0-10), analysed in full, and parked with an explicit re-open trigger. The formula reproduces exactly; the analysis surfaced two defects in it, one defect in the proposed replacement, and one governance question. **Defect 1**: with the full criteria set as denominator, untested rows pull a column to 5.0, the midpoint, so silence reads as "average model" rather than "unknown" (M8 scores 5.1 on ONE attempted criterion out of 22). **Defect 2**: ❌ = -1 sits strictly below 🚧 = 0, so running a test and publishing the failure costs a column more than never running it, contradicting the file's own "a ❌ is a result, not an embarrassment" stance. The proposed fix (attempted coverage `A/N` beside score-on-attempted `10 × (✅ + 0.5⚠️)/A`, negative weight dropped as double-counting) repairs both, but carries **defect 3, the user's catch: it is unstable at small samples**, since M8 would print 10.0 if its single gravity row validated, so a small-denominator column can top the board on one result. The governance question: any score DISAGREES with the existing column-ordering rule (M4 sits 3rd by ✅+⚠️ count and last by score). **The deferral decision**: the score-board already does the two jobs that matter, it separates covered from not-covered and it orders the columns, so a scalar score is not worth the author friction it would buy at today's coverage; revisit when a critical mass of tested rows exists, with a minimum-coverage gate or shrinkage toward a prior as the candidate cures for defect 3. Parked with it: the ⚠️-weight sensitivity band (only the ends are stable, M6 vs M7 flips across w = 0.3 to 0.7 and ties exactly at 0.5), the mechanical `regime` cut as the fair substitute for difficulty weighting, linter-owned recomputation so no number is ever hand-typed, and the pre-commitment that a score must not reorder columns so no later edit can demote a column for reporting negatives.
