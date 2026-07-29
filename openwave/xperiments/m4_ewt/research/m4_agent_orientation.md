# M4 Agent Orientation: the agent front door

> **If you are an AI agent and were told "read the m4_agent_orientation.md": this page is
> your bootstrap.** Load every document in § 1 into context, in order. Then follow the
> completion protocol in § 4: confirm you are oriented, summarize what you read, and
> declare yourself ready. From that point you can answer questions about M4, plan next
> moves, or run a task on a simple **"go task m4.1"** once the roadmap has rows.
>
> Humans are welcome here too; the model's human front door is
> [`__M4_model_briefing.md`](../__M4_model_briefing.md), and the research folder explains
> itself in [`README.md`](README.md).

## 1. The orientation reading list (load ALL, in this order)

| # | Doc | What it gives you |
| --- | --- | --- |
| 1 | [`__M4_model_briefing.md`](../__M4_model_briefing.md) | the column at a glance: identity, model profile, honest per-criterion status, help wanted |
| 2 | [`README.md`](README.md) | the two solvers and which one a question belongs to, the folder layout, what counts as a finished result |
| 3 | [`m4_theory_canonical.md`](m4_theory_canonical.md) | the spec of record. It is a skeleton: its § 0 maps where M4's theory record actually lives, and you should read those sources rather than treat the empty sections as absence of theory |
| 4 | [`M4_engine_upgrade.md`](M4_engine_upgrade.md) | the implementation record, phase by phase, with the validation evidence for each: the closest thing M4 has to a worked task history |
| 5 | [`M4_k_selectivity_Formalization.md`](M4_k_selectivity_Formalization.md) | the contributed nonlinear-stabilisation formalization behind the current K-selectivity work |
| 6 | [`m4_roadmap.md`](m4_roadmap.md) | the program. It is a skeleton with no tasks: its `CONVENTIONS` section is the part to read, since it defines IDs, row shape and how a task gets created |
| 7 | [`../theory/_CITATIONS.md`](../theory/_CITATIONS.md) | the source registry: which document backs which claim, and the never-fabricate identifier policy |
| 8 | [`AI_HYGIENE.md`](../../../../AI_HYGIENE.md) | the platform contract. Not optional, and it outranks anything else on this list |

Read alongside, for the lineage M4 inherits:
[`../../m3_wolff_lafreniere/research/0a_equations.md`](../../m3_wolff_lafreniere/research/0a_equations.md)
and
[`../../m3_wolff_lafreniere/research/0_STATUS.md`](../../m3_wolff_lafreniere/research/0_STATUS.md).

## 2. How tasks run here (M8.1 is the template)

M4 has no closed task in the platform's task format yet, so the reference execution lives
in another column. Read it as a worked example and reuse the shape:

| Phase | What to copy from [M8.1](../../m8_mit/research/tasks/m8_1_task_details.md) |
| --- | --- |
| PLAN | scope and definition of done; **pass/fail criteria pre-registered BEFORE any numerics** (its C1-C5 table); a blindspot pass naming what could silently go wrong; sub-experiments named with their artifact paths |
| EXECUTE | scripts, data and plots under `research/` with `m4_<n>_` prefixes; every claim backed by a runnable script; deviations from the plan logged as they happen, never reconstructed afterwards |
| FINISH | a method note per [`dev_docs/METHOD_NOTE.md`](../../../../dev_docs/METHOD_NOTE.md): equations first, an equation-to-code map, embedded plots, the audit record. Findings written into the task document. Honest status flips in the roadmap, the briefing and the [`MODELS.md`](../../../../MODELS.md) column, a negative syncing the same day as a positive would |
| REVIEW | a short review block in the task document: result per pre-registered criterion, issues found and dispositioned, the takeaway, the docs touched |

The finished-result standard is in [`README.md § 4`](README.md), and
[`m8_1_method_note.md`](../../m8_mit/research/findings/m8_1_method_note.md) is what one
looks like when it is done.

## 3. The one rule that is not a suggestion

From [`AI_HYGIENE.md`](../../../../AI_HYGIENE.md), the platform-wide contract:
**model output is a draft or a hypothesis, never a result, until verified by something
that is not a language model**: a runnable script, a hand-checked derivation, a
measurement, or the author's own authority. Show the script and the number, never a
verdict alone.

Two practices carry it:

| Practice | What it means |
| --- | --- |
| Adversarial audit | before a substantive claim is trusted or shared, an independent second agent tries to REFUTE it using its own script and its own method (different discretization, different derivation), returning a verdict per claim. The audit is recorded in the deliverable |
| Author gates | questions of intent, provenance or definition can only be answered by the model's author or the contributor who wrote the work. An agent never resolves them by inference. If a published definition underdetermines a computation, that is a question to ask, not a choice to make silently |

Route every unknown by who or what can settle it:

| Unknown type | Route |
| --- | --- |
| Machine-checkable | a script can decide it, so write the script |
| Author-gated | only the author can answer, so ask |
| Nature-gated | only an experiment or observation can decide, so register it as a falsifier with a threshold |
| Unknown-unknowns | no route exists, so build tripwires: the blindspot pass at PLAN, the deviations log during EXECUTE, and a self-quiz before anything is shared outward ("what would a hostile reader attack first?") |

## 4. M4-specific cautions

Things an agent working this column will hit, and which have burned time before:

| Caution | Detail |
| --- | --- |
| Which solver | a parameter study, a stability question, or any number destined for [`MODELS.md`](../../../../MODELS.md) is headless work under `research/`. The rendered launcher comes after, and runs the kernels the headless work validated. [`README.md § 2`](README.md) carries the reasoning and the citations |
| Imposed versus emergent | several M4 results are imposed by construction, charge sign among them. Never report an imposed property as a derived one; when unsure which a given property is, that is an author-gated question |
| Boundary contamination | Dirichlet reflections are a known contaminant of the soliton search. A stability claim should say what the boundary was doing |
| The papers are not in git | [`../theory/`](../theory/) holds third-party documents local-only and gitignored. Cite through [`_CITATIONS.md`](../theory/_CITATIONS.md), and never invent a DOI, page or quotation |
| The roadmap is empty on purpose | do not populate it speculatively. Propose tasks to the people extending the model and let them decide the sequence |

## 5. Completion protocol (what to print after reading)

Once §§ 1-4 are consumed, print, in this order:

| # | Print |
| --- | --- |
| 1 | A confirmation that you are ORIENTED on the M4 column and the platform contract |
| 2 | A one-line summary of EACH document you read, so a reader can verify nothing was skipped |
| 3 | The column's honest status from the briefing: what is built, what is imposed, what is open |
| 4 | A readiness statement: you can now (a) answer questions about M4 and its OpenWave context, (b) help plan next moves, and (c) execute a roadmap task on command once rows exist, where **"go task m4.1"** means open `tasks/m4_1_task_details.md`, complete its pre-registration, and run it per § 2 |
