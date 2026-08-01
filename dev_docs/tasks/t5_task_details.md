# T5: Evaluate issue-based task tracking (GitHub Issues + Projects v2)

> Roadmap row: [`../platform_roadmap.md`](../platform_roadmap.md). Status: ✅ **CLOSED
> 2026-08-01, DECLINED** (user decision). Owner: maintainers. Filed 2026-07-31 as an evaluation
> in which declining to migrate was a first-class outcome, and that is how it closed: task
> tracking stays in the repository roadmaps, and GitHub issues are reserved for platform
> defects. The decision, its price, and the issue-by-issue coverage audit are in
> [§ FINDINGS](#findings). Everything under PLANNING is the pre-decision record and is not
> rewritten.

## PLANNING

### The problem

Every edit to a roadmap or a task-details file is a repository change, so it costs a full
review cycle. In the 24 hours around the M8.5-A protocol freeze, the coordination edits alone
took three reviewed merges ([#382](https://github.com/openwave-labs/openwave/pull/382),
[#383](https://github.com/openwave-labs/openwave/pull/383),
[#384](https://github.com/openwave-labs/openwave/pull/384)): status stamps, settled working
assumptions, a landing map. For a maintainer that is merge overhead; for a fork-based author
every such edit is fork → branch → commit → PR → DCO → merge. The cost is linear in edits ×
authors, and the platform's direction is more authors.

The precedent already exists in miniature: contributor-facing physics tasks are tracked as
GitHub issues (the `theory` label, [#197](https://github.com/openwave-labs/openwave/issues/197)
onward) precisely so that picking one up needs no repository write.

### The proposal

Move the LIVE COORDINATION STATE to GitHub Issues plus a Projects v2 board (kanban columns
mirroring Backlog / In Progress / Done). The research record stays in git, unchanged.

| Migrates to issues + board | Stays in git |
| --- | --- |
| roadmap task rows: status, owner, gating, run order | everything under each column's `research/`: scripts, data, plots, findings |
| task-details PLANNING content (scope, plan, blindspots) → the issue body | frozen documents: pre-registrations, locked protocols, method notes |
| deviations log → dated issue comments | roadmap CONVENTIONS sections (standing rules other documents cite) |
| status discussion that today rides PR threads | roadmap CHANGE-LOG sections (dated records, append-only) |
| | `MODELS.md`, the standards documents, the linters |

**The boundary rule that decides every case:** an issue may hold only what nothing needs to
cite. Anything a frozen document pins, a paper links, or a findings note references lives in
git. Issues are coordination residue; git is the record.

At task close, the citable outcome lands in the column's `findings/` as it already does; the
closed issue is tracking history, not the record of the result.

### What this knowingly gives up

Filed as a proposal because the trade is real. Each loss with its mitigation:

| Loss | Mitigation |
| --- | --- |
| Commit-pinning: a roadmap row's state at a SHA is provable; issue bodies edit in place | the boundary rule. Nothing pinnable lives in an issue; frozen documents keep pinning git paths |
| The [`ROADMAP_STANDARDS.md`](../ROADMAP_STANDARDS.md) budgets, enforced by [`check_roadmaps.py`](../utils/check_roadmaps.py); nothing lints issue text | not actually a loss: the budgets exist to keep row TABLES scannable, and issues have no row tables (details live in the body and comments). Standard and linter retire from this repository at completion; see § Retirement |
| The review gate: every tracking edit lands through a reviewed, DCO-signed merge | under the intake model below, every STATE change stays a maintainer action, so the mediation survives; what disappears is diff review of maintainers' own tracking edits, which are largely self-merged today anyway |
| Single-commit atomicity and clone completeness: a fresh clone holds the whole project state | a scripted snapshot (`gh` export of issues + board into an archive folder) at milestones, so the repo periodically re-absorbs the tracking state |
| Automatic agent context: files load by read, and the repository is self-describing | `gh` CLI access documented in the orientation docs; agents fetch tracking state deliberately |

### Measured constraints (2026-07-31)

| Fact | Consequence |
| --- | --- |
| The most active model author's repository permission is **read**; the other active contributor is fork-based with none | shapes the intake model below: it is designed to need NO access grant at all |
| An org project named `openwave` already exists with 20 items | audited before anything is created: reuse, archive, or replace is a block-one decision, not an assumption |
| `gh` reads org projects with the present token scopes | the tooling path is open; write-scope confirmation happens at execution |
| GitHub issue and comment edit history is viewable but deletable by the content author | issues are not tamper-evident the way git is. This single fact routes every integrity-bearing artifact to git under the boundary rule |

### Intake and permission model (decided: no author triage or project write)

**Fixed before any consultation, not negotiated in it: authors do not receive triage or
project write.** Completion validation and board state are the maintainers'
quality-assurance role, and granting them away would dissolve the QA gate the platform runs
on. The consultation asks whether migration is worth it GIVEN this model, never whether the
model bends.

No access grant is needed under it. The split follows action frequency: the high-frequency
author actions become friction-free, and every STATE change stays a maintainer action.

| Actor | Can, with no grant beyond a GitHub account | Cannot |
| --- | --- | --- |
| Author | open task issues in their column, edit their own issue bodies (planning, proposed gates), comment on any issue (deviations, discussion, asks) | add labels, add an issue to the board, order the backlog, move cards, close a task as complete |
| Maintainer | label the issue, add it to the board under the author's column, order per the author's stated intent, move cards, validate completion before Done | |

Compared against the PR flow this is strictly less author friction: today an author's
tracking edit is fork → branch → PR → DCO → merge; under this model the same edit is an issue
comment, with maintainer mediation only at state changes.

**Two costs follow from the fixed line, recorded as named consequences rather than
negotiated away:**

| Cost | Consequence |
| --- | --- |
| The author's column-driving is mediated | an author cannot move their own cards or close their own tasks. The column is author-INTENT-driven (the intent, ordering, gates, and planning are the author's text) but maintainer-OPERATED, a real narrowing of "the author drives the column" as the model roadmaps state it |
| Maintainer overhead | every author-side state change costs a maintainer transcription + validation action. The migration RELOCATES state-change friction to the maintainers rather than deleting it, and that load scales with author count, the same variable the migration is meant to relieve |

That second row is the evaluation's sharpest question. The net depends on the actual edit
mix: details, deviations, and discussion (the high-frequency traffic) become free, while
state changes (rare per task) stay mediated at roughly today's maintainer cost. If the real
traffic turns out to be state-heavy, the migration buys little.

### Gating strength

Two kinds of gate ride the `Gated By` column today, and the migration treats them differently:

| Kind | Example | Where it lives after migration |
| --- | --- | --- |
| Scheduling gate | task B runs after task A | the issue and board (dependency links, the gate named in the issue body). Public, timestamped, sufficient |
| Evidential gate | a pre-registered outcome unlocks a task; a frozen protocol pins a path | the gate's DEFINITION stays in git (it is citable, so the boundary rule already routes it there); the issue carries a pointer to the git artifact |

Enforcement strength is unchanged under the intake model: marking a gate satisfied is a card
move, which stays a maintainer action, exactly as a `Gated By` edit today passes a maintainer
merge. What is NOT acceptable is an evidential gate whose only definition is issue text, since
edit history there is deletable; that is the one configuration the boundary rule exists to
forbid.

### Pins and freezes constrain the order

Frozen documents pin live task paths today: the M8.5-A protocol's header pairs it with
[`m8_5_task_details.md`](../../openwave/xperiments/m8_mit/research/tasks/m8_5_task_details.md)
and its § 10 pins that path. A frozen document cannot be edited to follow a move.

| Rule | Consequence |
| --- | --- |
| A column migrates only when no frozen document pins one of its live task paths | the M8 column waits for M8.5-A to close |
| Retired roadmap and task files become redirect stubs, never deletions | every existing link and pin keeps resolving |

### Rollout

| Step | What |
| --- | --- |
| 1 | consult the model authors on whether migration is worth it under the fixed intake model; a decline closes T5 |
| 2 | audit the existing org project; decide reuse vs replace |
| 3 | pilot on this platform roadmap only (the T-space), the maintainers' own tracking |
| 4 | evaluate against the losses table above, in writing |
| 5 | model columns opt in per author, each behind its pins check |
| 6 | update the scaffolding procedure: [`ONBOARDING_MODELS.md`](../../ONBOARDING_MODELS.md) § 3.3's artifact table creates a per-column roadmap FILE ("maintainer drafts, the author reorders"); post-migration it creates the column's issue set and board view instead. Sweep the other standards that reference roadmap files ([`AUTHOR_INVITATION.md`](../AUTHOR_INVITATION.md), [`PR_REVIEW_STANDARDS.md`](../PR_REVIEW_STANDARDS.md), [`CROSS_MODEL_TESTING.md`](../CROSS_MODEL_TESTING.md), [`REPRODUCE.md`](../../REPRODUCE.md)) |

### Retirement and handoff

When the last in-repo roadmap migrates, [`ROADMAP_STANDARDS.md`](../ROADMAP_STANDARDS.md) and
[`check_roadmaps.py`](../utils/check_roadmaps.py) retire from this repository: their subject,
budgeted row tables, no longer exists here. They are handed off rather than deleted: other
projects that adopted the format continue with local roadmaps, and the standard and linter
migrate to them. What stays in-repo (CONVENTIONS standing rules, append-only change-logs) is
prose, not row tables, and needs no format linter.

### Suggested definition of done

| # | Item |
| --- | --- |
| 1 | Authors consulted and the intake model recorded before any card exists |
| 2 | The boundary rule holds: nothing in any issue is cited by any repo document |
| 3 | Redirect stubs at every retired path; every pre-existing link resolves |
| 4 | The snapshot export script exists and has run once |
| 5 | The issue-anatomy standard exists (successor to the retired row-table standard) |
| 6 | The scaffolding procedure and the standards sweep (rollout step 6) are done |
| 7 | `ROADMAP_STANDARDS.md` + `check_roadmaps.py` retired per § Retirement and handoff |

### Blindspots

| Risk | Guard |
| --- | --- |
| Two-system drift: the board says done, the repo says pending | the board is authoritative for STATUS only; anything else appearing in an issue that contradicts a repo document is a defect by the boundary rule |
| Record content quietly accumulating in issue threads | at task close, the closing comment names the in-repo findings note; a close without one is incomplete |
| The migration itself churns a standard adopted 2026-07-28 | the pilot is scoped to one roadmap; nothing else moves until the written evaluation |
| The maintainer transcription queue lags and the board goes stale | the board is authoritative for STATUS only, so a lag is visible on its face rather than silently wrong; state changes batch at the review touchpoints that already exist |
| GitHub outage or account loss takes the tracking layer | the snapshot export; the repo remains sufficient to reconstruct state |

### Ownership + gating

**Gated by**: model-author consultation (a discussion thread, both active authors invited),
then the pilot sequence above. Per-column migration additionally waits on that column having
no frozen document pinning a live task path.

## DEVIATIONS LOG

**2026-08-01. The rollout never reached step 1, and the author consultation did not happen.**
§ Rollout made the consultation the first gate and a decline the closing condition. The decline
arrived from the user instead, on the evaluation's own trade-off tables, which makes the
consultation moot rather than skipped: there is nothing to consult about once the migration is
off. What the authors were going to be asked (is this worth it under a no-triage intake model)
is answered in the negative from the maintainer side, and the answer costs them nothing, since
declining leaves their workflow exactly as it is.

## FINDINGS

**DECISION, 2026-08-01 (user): DECLINED. Task tracking stays in the repository roadmaps, and
GitHub issues are reserved for platform issues.**

### Why

Three reasons, and the first is the one the evaluation had underweighted:

| Reason | Content |
| --- | --- |
| Version-tracked diffs are the product | A roadmap edit is a diff: attributable, reviewable, revertable, and readable as history, which is what matters when several authors touch the same files. An issue body edits in place and its edit history is deletable by the content author (§ Measured constraints). The evaluation had already routed every integrity-bearing artifact to git for that reason; the decision extends the same reasoning to the tracking layer rather than splitting it |
| Author permissions were never the whole cost | The intake model fixed before any consultation (no author triage, no project write) means the author's column stays maintainer-OPERATED under either system. The migration was therefore never going to hand an author their own board, only a cheaper comment channel |
| The maintainer cost survives the migration | § Intake called this "the evaluation's sharpest question" and it decides the outcome: every state change stays a maintainer transcription plus validation action, so the friction the task was filed on is relocated, not deleted, and a second system arrives to keep in sync |

### What it costs, stated plainly

The friction that opened the task is real and it stays. A fork-based author's tracking edit is
fork → branch → commit → PR → DCO → merge, and that is the price of the tracked diff. Two things
soften it and neither was a reason to migrate: the fork is one-time per author, and a tracking
edit normally rides the same pull request that carries the work it describes.

### What replaces the boundary rule

The migration's boundary rule (an issue holds only what nothing cites) is void with the
migration. The rule that replaces it is the reverse split, and it lives in the platform
roadmap's [§ CONVENTIONS](../platform_roadmap.md#conventions) as the single source: **tasks are
roadmap rows with a task document behind them; a GitHub issue is a platform defect or request,
something reproducible and closable that is wrong with the platform itself.** Proposals,
questions and offers of contribution belong in Discussions, and become rows if they become work.

### What the decision cancels

| Item in PLANNING | Status |
| --- | --- |
| § Rollout steps 2-6 (project audit, pilot, per-column opt-in, scaffolding sweep) | void |
| § Retirement and handoff | void. [`ROADMAP_STANDARDS.md`](../ROADMAP_STANDARDS.md) and [`check_roadmaps.py`](../utils/check_roadmaps.py) stay in this repository, with their subject intact |
| § Suggested definition of done, items 1-7 | void, superseded by the close-out below |
| § Pins and freezes, and the M8-column gate that [T4](t4_task_details.md) released | moot: nothing migrates, so nothing has to wait for a frozen document to stop pinning a live path |

### Close-out: every open issue got a roadmap home

The condition on closing the issues was that no open one is the only record of its work. Audited
2026-08-01 across all twelve open issues; four needed a new row, and one is not a task.

| Issue | Roadmap home | Kind |
| --- | --- | --- |
| [#197](https://github.com/openwave-labs/openwave/issues/197) effective Dirac | M5.10, backlog (body already archived 2026-07-02) | existing row |
| [#198](https://github.com/openwave-labs/openwave/issues/198) 3D pair-annihilation | M5.14, retired 2026-07-23; scope delivered by M5.21.4, M5.21.6, M5.20 | done, at higher rigor |
| [#200](https://github.com/openwave-labs/openwave/issues/200) lepton mass spectrum | M5.9, absorbed into the M5.21 electron hunt 2026-07-19; the live successor is M5.21.11 | existing home |
| [#201](https://github.com/openwave-labs/openwave/issues/201) K-selectivity | **new**: [M4.1](../../openwave/xperiments/m4_ewt/research/tasks/m4_1_task_details.md) | new row + archive |
| [#202](https://github.com/openwave-labs/openwave/issues/202) emergent Coulomb | **new**: [M4.2](../../openwave/xperiments/m4_ewt/research/tasks/m4_2_task_details.md) | new row + archive |
| [#209](https://github.com/openwave-labs/openwave/issues/209) electron gravitational mass | M5.23.6, backlog | existing row |
| [#213](https://github.com/openwave-labs/openwave/issues/213) dynamic demo suite | M5.23.3, backlog | existing row |
| [#247](https://github.com/openwave-labs/openwave/issues/247) particle field configurations | not a task: adopted as a standing requirement, the per-particle field-config section every briefing carries ([`ONBOARDING_MODELS.md`](../../ONBOARDING_MODELS.md) § 3.3), and the prescription recorded in the M7 roadmap | standard, not a row |
| [#248](https://github.com/openwave-labs/openwave/issues/248) collect the M6 author's papers | satisfied by the 2026-07-20 M6 refresh harvest: 29 latest-version records in [`_CITATIONS.md`](../../openwave/xperiments/m6_ouroboros/theory/_CITATIONS.md), including the primary dynamics paper | done |
| [#290](https://github.com/openwave-labs/openwave/issues/290) regularized vortex-loop | M5.19, closed COMPLETE 2026-07-10 with three adversarial audits | done |
| [#298](https://github.com/openwave-labs/openwave/issues/298) SMT solvers for combinatorial species | **new**: [T6](t6_task_details.md), deferred on arrival. The thread was already answered and acknowledged; the row exists because that answer promised a follow-up if a consumer appeared, and a closed issue is not where a promise survives | new row, deferred |
| [#324](https://github.com/openwave-labs/openwave/issues/324) stiffness ladder | **new**: [M5.21.13](../../openwave/xperiments/m5_liquid_crystal/research/tasks/m5_21_13_task_details.md), with a premise check on the δ direction | new row + archive |

Each new row's task document carries the issue body verbatim, so closing the issue destroys no
record. The one non-task ([#247](https://github.com/openwave-labs/openwave/issues/247)) closes on
the standard it became.

### One thing the exercise proved on its own

The audit found three open issues describing work already delivered (#198, #248, #290) and one
whose scope had been absorbed under a different ID (#200). None of that drift was visible from
the issue side, and all of it was visible from the roadmaps. That is the same argument as the
first reason above, arriving as evidence instead of principle.
