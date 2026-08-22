# M9 Agent Orientation: the agent front door

> **If you are an AI agent and were told "read the m9_agent_orientation.md": this page
> is your bootstrap.** Load every document in § 1 into context, in order. Then follow
> the completion protocol in § 4: confirm you are oriented, summarize what you read,
> and declare yourself ready. From that point the author can ask you anything about
> the M9 model, plan next moves with you, or start work with a simple
> **"go task m9.2"** (or any roadmap task).
>
> Humans are welcome to read this page too; the model's human front door is
> [`__M9_model_briefing.md`](../__M9_model_briefing.md).

## 1. The orientation reading list (load ALL, in this order)

| # | Doc | What it gives you |
| --- | --- | --- |
| 1 | [`__M9_model_briefing.md`](../__M9_model_briefing.md) | the column at a glance: identity, model profile, honest status, help wanted |
| 2 | [`research/m9_theory_canonical.md`](m9_theory_canonical.md) | THE SPEC OF RECORD (canonical when docs disagree): the arena, the effective action, the locked conventions, the declared opens, and the consumption rules that keep particle rows unscored |
| 3 | [`research/m9_background.md`](m9_background.md) | the gap map (what NSM has, what it lacks), the evidence-weight grading, the onboarding evaluation of record, and the admission terms |
| 4 | [`research/m9_roadmap.md`](m9_roadmap.md) | the program: tasks, gates, current status, what is startable now. The author creates and manages these tasks directly |
| 5 | [`research/m9_platform_pointers.md`](m9_platform_pointers.md) | the cross-model reading map and the house standards; ALSO load its § 1 platform-contract docs, [`AI_HYGIENE.md`](../../../../AI_HYGIENE.md) above all |
| 6 | [`research/tasks/m9_1_task_details.md`](tasks/m9_1_task_details.md) + [`research/findings/m9_1_hehl_datta_note.md`](findings/m9_1_hehl_datta_note.md) | **the worked template**: a closed task run with pre-registered criteria, a mutation-tested gate, a second-method audit, and an honest FAIL shipped as a scored result; § 2 below tells you what to copy from it |
| 7 | [`theory/_CITATIONS.md`](../theory/_CITATIONS.md) | the source registry: which record backs which claim; papers live in the author's repo, never here; the never-fabricate identifier policy |

## 2. How tasks run here (M9.1 is the template)

| Phase | What M9.1 did, to copy |
| --- | --- |
| PLAN | scope + definition of done; **pre-registered pass/fail criteria written BEFORE any numerics** (its C1-C5 table); conventions locked in the canonical so the computation cannot drift toward the target |
| EXECUTE | scripts / data / findings under `research/` with the task's `m9_<id>_` prefix, where `<id>` is a registered roadmap row (Gate A10 in [`dev_docs/PR_REVIEW_STANDARDS.md`](../../../../dev_docs/PR_REVIEW_STANDARDS.md)); every claim backed by a runnable script; deviations logged as they happen |
| FINISH | a finding note per [`dev_docs/METHOD_NOTE.md`](../../../../dev_docs/METHOD_NOTE.md) shape: equations first, equation-to-code map, the audit record; findings written into the task doc; honest status flips (a negative is a result and syncs the same day) |
| REVIEW | a short review block in the task doc: results per pre-registered criterion, issues found and dispositioned, findings takeaway, docs touched |

## 3. The working rules this column runs under

These bind the author's agents exactly as they bind the maintainers' (they are the
admission terms of [discussion #442](https://github.com/openwave-labs/openwave/discussions/442)
plus the platform contract, in one place):

| Rule | What it means at the keyboard |
| --- | --- |
| GitHub thread FIRST, at every session start | before any other work, and before any commit, read the open PR's conversation and review state: `gh pr view 441 --repo openwave-labs/openwave --comments` and `gh pr view 441 --repo openwave-labs/openwave --json reviews,isDraft` (adjust the number to the live PR). Maintainer decisions, verified merge slices, and draft conversions land there; work done blind to that thread can void itself. If anything in the thread conflicts with the current plan, STOP and surface it to the author before proceeding |
| One task, one PR | each roadmap task lands as its own pull request; nothing else rides along |
| Draft until quiescent | while commits are still landing, the PR is marked draft; review runs on a settled head, and commits pushed mid-review void the parts of the review they touch |
| Task ids resolve | a new campaign gets its roadmap row and `tasks/m9_<id>_task_details.md` BEFORE its scripts take the id; no artifact carries an id that names a different task or no task |
| Run before write | no new process document (packet, addendum, freeze, audit) until the experiment the last one governs has RUN; pre-registrations cap near 8,000 words ([`PR_REVIEW_STANDARDS.md`](../../../../dev_docs/PR_REVIEW_STANDARDS.md) § 12.2). Rules are priced in reader-hours, not tokens |
| Maintainer pacing | review and merging follow maintainer schedule and resource availability; the queue absorbs throughput as per-task PRs, not as one moving head |
| Author gates | questions of intent, provenance, or definition are the author's alone; an agent never resolves them by inference. If a paper underdetermines a computation, that is a question to ask, not a choice to make silently |
| Adversarial audit | before any substantive claim is trusted or shared: an independent second pass tries to REFUTE it with its own script and its own method; the audit is recorded in the deliverable. M9.1's second-method sympy audit is the local example |
| Precision then efficiency | **Priority 1 is precision. Priority 2 is efficiency.** Physics, logic, and math first. Prefer exact bases and multiprecision (`mpmath`) over float64 LAPACK when the object admits it. Do not thin a pre-register, drop digits, or switch to a noisier eigensolve to save wall time. After the result is exact enough to trust, compiled helpers are allowed; the Nuitka and python-flint/ARB helpers (`compile_nuitka.sh`, `m9_flint.py`) are parked with the campaign at `a5640709` and return with the per-task PR that needs them. Do not replay Papers 70--73. Do not replace mpmath with float64 for speed |

And the one non-negotiable, from [`AI_HYGIENE.md`](../../../../AI_HYGIENE.md) (the
platform-wide contract, not a suggestion): **model output is a draft or hypothesis,
never a result, until verified by something that is not a language model**: a runnable
script, a hand-checked derivation, a measurement, or the author's own authority.
Agents show the script and the number, never a verdict alone.

## 4. Completion protocol (what to print after reading)

Once §§ 1-3 are consumed, print, in this order:

| # | Print |
| --- | --- |
| 1 | A confirmation that you are ORIENTED on the M9 column and the platform contract |
| 2 | A one-line summary of EACH document you read (so the author can verify nothing was skipped) |
| 3 | The column's current status from the roadmap: what is done, what is startable, what is gated |
| 4 | A readiness statement: you can now (a) answer the author's questions about the M9 model and its OpenWave context, (b) help plan next moves, and (c) execute a roadmap task on command: **"go task m9.2"** means open [`research/tasks/m9_2_task_details.md`](tasks/m9_2_task_details.md), confirm its pre-registered gates with the author, and run it per §§ 2-3 |
