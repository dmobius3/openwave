# M9 platform pointers: where to read, for the author and the author's AI agents

> **Purpose.** A curated map of the OpenWave repository written to be CONSUMED BY AI
> AGENTS working for the M9 author (and by the author directly). Each entry gives a
> linked short name and what it teaches for the M9 program (gravity certification on
> an installed-matter effective theory: [`m9_background.md`](m9_background.md)).
> Every link was verified at scaffold time (2026-08-15).
>
> **If you are an AI agent:** load § 1 before generating any claim about this repo,
> then read only the sections the current task needs. Model output is a draft, never
> a result, until backed by a runnable script (the repo-wide contract).

## 1. Read first (the platform contract)

| Doc | Why |
| --- | --- |
| [`AI_HYGIENE.md`](../../../../AI_HYGIENE.md) | the working contract for AI-assisted research here: script-backed claims only, adversarial audit of substantive derivations, author-gated questions stay with the author |
| [`ONBOARDING_MODELS.md`](../../../../ONBOARDING_MODELS.md) | the onboarding standard M9 was admitted under, including § 3.6 headless-first: research scripts answer questions, rendering waits for canonical physics |
| [`MODELS.md`](../../../../MODELS.md) | the shared coverage matrix and its status semantics. M9's native rows are the two gravity rows; cells flip only via runnable script + honest research note, and per [`m9_theory_canonical.md`](m9_theory_canonical.md) § 7 the particle rows are never scored from this column |
| [`dev_docs/METHOD_NOTE.md`](../../../../dev_docs/METHOD_NOTE.md) | the reporting standard for any substantive result: equations first, equation-to-code map, embedded figures, adversarial audit recorded |
| [`dev_docs/PR_REVIEW_STANDARDS.md`](../../../../dev_docs/PR_REVIEW_STANDARDS.md) | what a review here verifies and how findings are tiered (§ 12.1: computation may block, wording never costs a round). Intake row 6 and Gate A10 encode the two working rules this column was admitted under: a PR under review is quiescent (draft while work continues), and every artifact id resolves to a registered roadmap task |
| [`REPRODUCE.md`](../../../../REPRODUCE.md) | the task-id glue convention: one id prefixes the task doc, its scripts, its data, and its findings, so any artifact resolves to the task doc that holds its commands |
| [`CONTRIBUTING.md`](../../../../CONTRIBUTING.md) | setup, fork → branch → PR flow, DCO sign-off, Black / PEP 8 as the style target |

M9's own spec of record lives in [`m9_theory_canonical.md`](m9_theory_canonical.md);
the program in [`m9_roadmap.md`](m9_roadmap.md).

## 2. The platform's gravity sector (what exists to compare against)

M9 is the platform's first gravity-native column. The existing gravity evidence is
thin by design (the emergence columns reach gravity last), which is exactly the
comparison surface M9's certificates land on.

| Pointer | Content |
| --- | --- |
| [`MODELS.md`](../../../../MODELS.md) gravity rows | Newton limit (GEM) and metric phenomena: the two rows M9 can natively score, and the current cross-model state of each (M5 ⚠️ on the coupling mechanism, the force itself not computed; everything else 🚧) |
| [`m5_8_2q_delta_scaling.py`](../../m5_liquid_crystal/research/scripts/m5_8_2q_delta_scaling.py) | the platform's measured gravity artifact: in M5, gravity enters via the boost tilt of the time axis (GEM \(\propto (b\cdot g)^2\)). A useful contrast: M5 measures a coupling mechanism on a lattice; M9 certifies an action-level limit. The two meet, eventually, at the same Newton row |
| [`dev_docs/CROSS_MODEL_TESTING.md`](../../../../dev_docs/CROSS_MODEL_TESTING.md) | scoring rules if M9 work ever touches another column's structures, and the routing convention for author-gated questions |
| The M9.2 route, stated | the roadmap's Newton task derives \(1/r^2\) from 3-d Poisson of the inherited Einstein equations, deliberately NOT the M5 GEM route; when both rows are eventually filled, the matrix carries two independent routes to the same limit, which is the platform working as intended |

## 3. The certification-column precedent (M8 is the sibling)

M8 (Mode Identity Theory) is the other externally-authored, certification-first
column, and its records are the worked precedent for how M9 tasks should run.

| Pointer | What it teaches |
| --- | --- |
| [`m8_1_task_details.md`](../../m8_mit/research/tasks/m8_1_task_details.md) + [`m8_1_method_note.md`](../../m8_mit/research/findings/m8_1_method_note.md) | the certification-gate template: pre-registered pass/fail criteria before numerics, an independent second-method audit, a defect found and dispositioned in the audit record. M9.1 already followed this shape; keep following it |
| [`m8_roadmap.md`](../../m8_mit/research/m8_roadmap.md) | a mature externally-authored roadmap under [`dev_docs/ROADMAP_STANDARDS.md`](../../../../dev_docs/ROADMAP_STANDARDS.md): row budgets, gates, the task-doc-is-the-record premise |
| The M8 review history (PRs [#402](https://github.com/openwave-labs/openwave/pull/402), [#436](https://github.com/openwave-labs/openwave/pull/436)) | what review rounds cost and how the platform prices them now: pre-registrations cap near 8,000 words, no new process artifact until the pending experiment has RUN, computation rigor stays, prose rigor is comment-tier ([`dev_docs/PR_REVIEW_STANDARDS.md`](../../../../dev_docs/PR_REVIEW_STANDARDS.md) § 12.2). An agent drafting M9 documents weighs rules in reader-hours, not tokens |

## 4. What M9 consumes from shared code (deliberately little)

| Item | Status |
| --- | --- |
| `openwave/common/` engines and the GGUI stack | not consumed. M9's instruments are algebra (numpy/sympy) and small lattice-fermion diagonalizations; there is no field evolution to render. Headless-first resolves to headless-only until a gravity cell exists ([`m9_roadmap.md`](m9_roadmap.md) conventions) |
| `numpy`, `scipy`, `sympy` | already platform dependencies (`pyproject.toml`); research scripts stay dependency-light and self-contained, one script one claim |
| Shared docs layer | fully consumed: the § 1 contract, the roadmap linter (`dev_docs/utils/check_roadmaps.py`), and the review standards are the interfaces this column actually plugs into |

## 5. House standards for M9 tasks (so contributions merge smoothly)

| Standard | Rule |
| --- | --- |
| Naming | scripts / data / findings under `research/` with `m9_<id>_` prefixes, where `<id>` is a registered row on [`m9_roadmap.md`](m9_roadmap.md) (Gate A10). No campaign files outside the roadmap's id space; shared helper modules take the id of the task that owns them or live in `scripts/utils/` |
| Task ownership | the author creates and manages roadmap tasks directly (admission term, discussion #442); registration is a roadmap row plus a `tasks/m9_<id>_task_details.md`, in the same PR as the artifacts |
| PR shape | one task, one PR, draft until quiescent; review and merging follow maintainer schedule and resource availability |
| Status honesty | `MODELS.md` cells flip only with a runnable script + a research note documenting pass/fail; negatives are results and are wired into the column the same day |
| Pre-registration | gates, conventions, and success criteria written BEFORE numerics, in the task doc; forks reported with all numbers, never tuned toward published values |
| Audit | substantive claims get an independent adversarial pass (own script, own method) before they are trusted; record the audit in the deliverable |
| Sources | papers and PDFs stay in the author's repo; [`../theory/_CITATIONS.md`](../theory/_CITATIONS.md) is the tracked registry; never fabricate an identifier |
| Style | markdown tables for structured content; status icons only ✅ ⚠️ ❌ 🔶 🚧; escape literal pipes as `\|` inside table cells; relative links; English throughout |
