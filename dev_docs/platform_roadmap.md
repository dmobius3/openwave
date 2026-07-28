# PLATFORM ROADMAP

> The task roadmap for **platform-wide** work: the artifacts every model column shares, and the tooling that keeps them honest. Live work: [IN PROGRESS](#in-progress) → [BACKLOG](#backlog); the record: [DONE](#done). Reading rules (what belongs here, the T-ID scheme, task-doc anatomy): [§ CONVENTIONS](#conventions). Dated decisions and migrations: [§ CHANGE-LOG](#change-log).
>
> Model-specific physics does NOT live here. Each model keeps its own roadmap ([M5](../openwave/xperiments/m5_liquid_crystal/research/m5_roadmap.md), [M8](../openwave/xperiments/m8_mit/research/m8_roadmap.md), and so on); this file is for the work that would otherwise be filed under whichever model happened to notice it.

---

## IN PROGRESS

> The single currently-running platform task (one at a time; empty = between tasks).

| TaskID | Title | Description | Gated By |
| --- | --- | --- | --- |
| | | (between tasks) | |

## BACKLOG

> Future platform tasks. **Row order = the run sequence.** Deferred rows carry their re-open trigger in `Gated By`.

| TaskID | Title | Description | Gated By |
| --- | --- | --- | --- |
| [T2](tasks/t2_task_details.md) | MODELS.md normalized score: a fair single number per column | 🔶 DEFERRED 2026-07-28, same day as the proposal (was M5.30): the score-board already separates covered from not-covered and orders the columns, so a scalar score waits for a critical mass of tested rows. The proposed formula scores silence as average and ranks ❌ below 🚧; the coverage-plus-score-on-attempted repair is unstable at small samples. Analysis and candidate cures in the task doc. | critical mass of tested rows (user call) + [T1](tasks/t1_task_details.md) final criteria set + user "go" |

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

| TaskID | Title | Description | Completed |
| --- | --- | --- | --- |
| [T1](tasks/t1_task_details.md) | MODELS.md criteria housekeeping: split the coarse rows, name the simplest passing test per row | Criteria set 22 → 31 (was M5.29): bundled rows split (weak, baryons, neutrinos, gravity, strong), Lorentz covariance + running coupling + deuteron + nuclear structure added, every row names its simplest passing test (tier decode in the task doc), linter extended to guard the test column, all counts re-derived, no icon upgraded. | 2026-07-28 |

## CHANGE-LOG

**PLATFORM ROADMAP CREATED, TWO TASKS MIGRATED IN (2026-07-28, user decision).** MODELS.md work was being tracked on the M5 roadmap because M5 is where the proposals happened to arrive, which mis-filed it: the criteria set and the score-board belong to every column, not to the liquid-crystal one. This roadmap now owns platform-wide work, with its own `T<n>` ID space and its own [`tasks/`](tasks) folder. Migrated: **M5.29 → [T1](tasks/t1_task_details.md)** (criteria housekeeping, with its convo record [`tasks/t1_convo.md`](tasks/t1_convo.md)) and **M5.30 → [T2](tasks/t2_task_details.md)** (the normalized score, deferred). Same move: the MODELS.md linter went from `dev_docs/check_models_md.py` to [`dev_docs/utils/check_models_md.py`](utils/check_models_md.py), a procedural test rather than a standards document, with every reference updated ([`MODELS.md`](../MODELS.md), [`PR_REVIEW_STANDARDS.md`](PR_REVIEW_STANDARDS.md) § 7.1 and the review-command block, the two task docs) and its `MODELS.md` path resolution fixed for the extra directory level. The [M5 roadmap](../openwave/xperiments/m5_liquid_crystal/research/m5_roadmap.md) keeps its 2026-07-27 and 2026-07-28 change-log entries as the provenance record of how these two tasks arrived.

**T1 STAGED, THEN EXPANDED BY THE PROPERTIES-TO-TEST SLIDE (2026-07-27 → 2026-07-28, as M5.29).** The author's 2026-07-27 reply proposed the criteria restructure (split the coarse rows, per-row simplest passing test, the missing Lorentz covariance row, Larmor precession as the magnetic-force test) and it was staged the same day. On 2026-07-28 a new group thread ([`tasks/t1_convo.md`](tasks/t1_convo.md)) added a "Properties to test for SM + gravity?" slide to the public talk deck: a 16-line test list with explicit priority shading, set beside a screenshot of the [`MODELS.md`](../MODELS.md) matrix itself, plus an open question to the group. That converts the test column from a design exercise into a transcription-and-fit exercise, confirms every 2026-07-27 ask, and adds three candidate new rows (running coupling; deuteron binding + quadrupole moment; larger nuclei / halos). It also raises the stakes: the deck's opening slide names OpenWave and carries the matrix icons into a colloquium context, so the re-derivation has a public external audience. Any reply to the group question is user-gated and follows the T1 output. Same day, T1 was promoted ahead of [M5.22](../openwave/xperiments/m5_liquid_crystal/research/tasks/m5_22_task_details.md): docs-only, no dependency on census results, and freezing the criteria set first gives the census stable rows to file into.

**T2 PROPOSED AND DEFERRED THE SAME DAY (2026-07-28, user proposal + user decision, as M5.30).** A per-column 0-10 score was proposed for the score-board (icon weights 1 / 0.5 / -1 / 0 over the full criteria set, rescaled onto 0-10), analysed in full, and parked with a re-open trigger. It reproduces exactly and has two defects: untested rows pull a column to the 5.0 midpoint, so silence reads as "average model" rather than "unknown" (M8 scores 5.1 on ONE attempted criterion out of 22), and ❌ = -1 sits below 🚧 = 0, so publishing a failure costs more than never running the test. The repair (coverage `A/N` beside score-on-attempted `10 × (✅ + 0.5⚠️)/A`) fixes both and adds a third defect, the user's catch: at small samples it is unstable, since M8 would print 10.0 on one validated row. **Deferred** because the score-board already separates covered from not-covered and orders the columns, so a scalar is not worth the author friction at today's coverage. Revisit on a critical mass of tested rows, with a minimum-coverage gate or shrinkage as the candidate cures. Full analysis, sensitivity band and the no-reordering pre-commitment: [T2](tasks/t2_task_details.md).
