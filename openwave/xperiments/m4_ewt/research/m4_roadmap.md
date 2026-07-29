# M4 / EWT, ROADMAP

> The M4 program: what is running, what is queued, what is closed. This file is a
> **skeleton**, deliberately empty of tasks. The people extending M4 fill it. Read
> [`## CONVENTIONS`](#conventions) first: it explains the ID scheme, how to add a task,
> and which M8 documents to copy from. Model orientation lives in
> [`README.md`](README.md) and [`../__M4_model_briefing.md`](../__M4_model_briefing.md).

## IN PROGRESS

| TaskID | Title | Description | Gated By |
| --- | --- | --- | --- |

## BACKLOG

> Queued tasks, row order is the run sequence. The open problems worth turning into
> rows are listed in the briefing's Status and Help Wanted tables.

| TaskID | Title | Description | Gated By |
| --- | --- | --- | --- |

## CONVENTIONS

**Task IDs.** `M4.<n>`, optionally `M4.<n>.<m>` for a sub-task. Assigned in creation
order, never reused, including after a renumber. A task keeps its ID for life, so
scripts, data and plots can be prefixed `m4_<n>_` and stay grouped forever.

**Row shape.** Every roadmap in this repository follows
[`dev_docs/ROADMAP_STANDARDS.md`](../../../../dev_docs/ROADMAP_STANDARDS.md), which the
checker enforces:

```bash
python3 dev_docs/utils/check_roadmaps.py openwave/xperiments/m4_ewt/research/m4_roadmap.md
```

The premise behind it: **a roadmap row is a preview, the task document is the record.**
A `Description` cell is capped at 65 words. Anything that does not fit belongs in
`tasks/m4_<n>_task_details.md`, not in a longer row.

**How to add a task.** Four steps, none of which need permission:

| # | Step |
| --- | --- |
| 1 | Add a row to `BACKLOG` with the next free ID, a title under 15 words, and a description under 65 that says what the task does and what decides it |
| 2 | Create `tasks/m4_<n>_task_details.md`: the scope, the method, and the **pass/fail criteria written before any numbers exist**. This is the part that makes the result mean something later |
| 3 | Run it headless. Scripts to `scripts/`, data to `data/`, plots to `plots/`, all prefixed `m4_<n>_` |
| 4 | Write the findings into the task document, move the row to `DONE` with the verdict and the date, and sync any status that changed in the briefing or in [`MODELS.md`](../../../../MODELS.md) |

**Worked examples to copy.** M8 is the closest reference because its tasks are recent and
complete end to end:

| To see | Read |
| --- | --- |
| A roadmap with live, gated and closed rows | [`m8_roadmap.md`](../../m8_mit/research/m8_roadmap.md) |
| A task document, from plan through execution to review | [`m8_1_task_details.md`](../../m8_mit/research/tasks/m8_1_task_details.md) |
| Pre-registered pass/fail criteria, written before the numbers | the C1-C5 table in that same task document |
| A finished result write-up: equations first, equation-to-code map, embedded figures, adversarial audit | [`m8_1_method_note.md`](../../m8_mit/research/findings/m8_1_method_note.md) |
| A task that gates itself on something else finishing | M8.7 in the roadmap, the rendering port |

M5 and M7 carry longer histories if a second example helps:
[`m5_roadmap.md`](../../m5_liquid_crystal/research/m5_roadmap.md),
[`m7_roadmap.md`](../../m7_hydroboros/research/m7_roadmap.md).

**A negative closes a task.** "This does not work, here is the script, here is why" is a
result and moves its row to `DONE` exactly like a positive one. The briefing's Status
table already carries several honest ❌ rows; each is a task waiting to be written.

## DONE

| TaskID | Title | Description | Completed |
| --- | --- | --- | --- |

## CHANGE-LOG

**2026-07-29.** Roadmap created as a skeleton alongside the
[`research/`](README.md) scaffold. No tasks planned: the program belongs to whoever
extends the model, and the structure is here so the first task has somewhere to land.
