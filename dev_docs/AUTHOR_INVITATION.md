# Inviting a Model Author

> **What this page is.** One real invitation, sent 2026-07-29, kept as a **reference**: a maintainer answering a model author who asked, on the Models-of-Particles list, what a platform like this could possibly test about a theory whose Lagrangian is already published.
>
> **It is not a template.** No future invitation should be this message with the names swapped. Every author arrives with a different question, a different framework, and a different reason for being skeptical, and the answer that lands is the one that starts from what they actually asked. What this page carries is the set of **pieces an invitation should contain**, the **links to embed**, and one worked example of the register that works: direct, specific, and honest about what the platform does not supply.
>
> The path being invited into is [`ONBOARDING_MODELS.md`](../ONBOARDING_MODELS.md). This page is the outreach that precedes it.

## The pieces an invitation should carry

Ordered as they worked in the reference message: answer the question first, invite second. An author who has not been told what the platform *does* has no reason to read four steps of how to join it.

| Piece | What it has to establish |
| --- | --- |
| What OpenWave is | Open-source, volunteer-maintained, for testing candidate field-theoretic models of matter. Say what "model" means here, since the word is overloaded: a framework claiming to build particles, forces and waves out of an underlying field, written as a closed set of equations |
| How a model gets tested | The loop, concretely: implement the Lagrangian in a GPU solver, energy-minimize to find which stable configurations the model actually has, run dynamics on them, read off observables. The point to land: these are properties of *solutions*, not of the action, so they cannot be read off the paper by inspection |
| Where results land | [`MODELS.md`](../MODELS.md): one column per framework, shared criteria rows, every cell earned by a runnable script plus a research note or left marked "not yet tested". Reproducibility, not orthodoxy; a documented negative counts as much as a positive |
| That the arena is not fixed | Name the existing substrates and how unalike they are, so "my framework does not look like yours" stops being a reason not to apply. Pair it with the boundary: the equations stay the author's, and the platform cannot supply them |
| That an incomplete model still qualifies | The entry bar is not a finished theory. M8 was accepted with no equation of motion at all, and supplying one became the program: the author brings the arena and the target structure, the platform brings Lagrangian candidates, simulation engineering and grading. Say it, because an author whose framework has a known gap will otherwise assume the gap disqualifies it |
| The formal invitation | Say the slot is open and the invitation is formal. Mechanics in one line: one discussion post in the New Model category, no form, a maintainer picks it up |
| The one page to read | [`ONBOARDING_MODELS.md`](../ONBOARDING_MODELS.md), plus the highest-leverage instruction available: point an AI assistant at it. Summarize its four steps rather than listing headings |
| Bring your own compute | State it plainly and early enough to be fair. Runs are author-driven, on the author's own AI tokens and hardware; the platform supplies criteria, scaffold, review, standards and prior art. The column carries the author's name, so the author answers challenges to it |
| A live worked example | A named column, with dates and specifics. What makes it persuasive is the unflattering half: the certification gate that ran *before* any endorsement, and the cells still reading "not yet tested" |
| The AI-hygiene contract | [`AI_HYGIENE.md`](../AI_HYGIENE.md) in one sentence: an AI model's output is a draft, never a result, until something that is not a language model confirms it. This is what separates the platform from the venues that make careful people wary |
| A first step that costs nothing | End on an action with no commitment attached. Here: hand the paper and the onboarding page to an AI assistant and run the STEP 1 self-evaluation, which produces the application material as a side effect |

**On the worked example.** Use whichever column is furthest along in onboarding at the time of writing, and cite it by its record rather than by adjectives: the application discussion, the scaffold date, what was verified and how, what the author has contributed since, and what remains untested. An invitation that only lists successes reads like recruitment; one that reports the open cells reads like the platform it is describing.

**Register notes.** Answer the recipient's question before making any offer. Do not oversell: the sentence about what the platform cannot supply does more work than any claim about what it can. Name a model's gap as a program with a task ID, never as a verdict on the model or its author; "no dynamical field at start" and "no dynamical field" read very differently to the person who wrote the theory. Never offer a call or a meeting; end on the technical next step. Keep every number checkable against the repository on the day it is sent.

## Links to embed

The reference message carries these as inline links on the anchor text shown, rather than as pasted URLs.

| Anchor in the message | Target |
| --- | --- |
| MODELS.md | `https://github.com/openwave-labs/openwave/blob/main/MODELS.md` |
| ONBOARDING_MODELS.md | `https://github.com/openwave-labs/openwave/blob/main/ONBOARDING_MODELS.md` |
| "New Model" category | `https://github.com/openwave-labs/openwave/discussions/categories/new-model` |
| AI_HYGIENE.md | `https://github.com/openwave-labs/openwave/blob/main/AI_HYGIENE.md` |
| M8-MIT | `https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m8_mit/__M8_model_briefing.md` |

[`CONTRIBUTING.md`](../CONTRIBUTING.md) belongs in the set too when the recipient has not already been given it; in this thread another list member had linked it a few hours earlier, which is why the message does not repeat it.

## Facts the reference message asserts, and where they come from

Any reused claim gets re-checked against the repository on the day it is sent. These were current on 2026-07-29.

| Claim | Source of record |
| --- | --- |
| M8 applied with a single discussion post and was scaffolded from it | [Discussion #312](https://github.com/openwave-labs/openwave/discussions/312), 2026-07-21; [`m8_roadmap.md` § CHANGE-LOG](../openwave/xperiments/m8_mit/research/m8_roadmap.md) |
| Two blind agents re-derived the bedrock eigenvalue theorem to 10-digit precision, unprompted by the printed constants | [M8.1](../openwave/xperiments/m8_mit/research/tasks/m8_1_task_details.md) certification gate, 2026-07-21, and its [method note](../openwave/xperiments/m8_mit/research/findings/m8_1_method_note.md) |
| Two author pull requests merged in the first week | [PR #350](https://github.com/openwave-labs/openwave/pull/350) (M8.2) and [PR #362](https://github.com/openwave-labs/openwave/pull/362) (M8.3) |
| One cell earned, the rest "not yet tested" | The M8 column in [`MODELS.md`](../MODELS.md) and its score-board |
| M8 was onboarded without an equation of motion, and supplying one is the program | [`m8_roadmap.md`](../openwave/xperiments/m8_mit/research/m8_roadmap.md) opening line ("a field-dynamics collaboration"), [M8.2](../openwave/xperiments/m8_mit/research/tasks/m8_2_task_details.md) pre-registration locked 2026-07-27, [M8.4](../openwave/xperiments/m8_mit/research/tasks/m8_4_task_details.md) the Lagrangian-family survey |
| Substrates named for M5, M6, M8 | The three model briefings; [`MODELS.md`](../MODELS.md) column headers |

A fourth M8 fact was drafted and cut for length, worth restoring when the audience needs the strongest available evidence that the process corrects rather than confirms: the author's own reproducer script ([M8.3](../openwave/xperiments/m8_mit/research/tasks/m8_3_task_details.md)) found a dropped scalar-zeta term in the author's own published mass-formula page, which was then corrected upstream.

## The reference message (sent 2026-07-29)

Sent to the Models-of-Particles list in reply to a model author's question, "Can someone explain to me how this system simulates or tests physical theories? My Lagrangian is in my paper but I don't understand how this system would test anything interesting."

Reproduced as sent, with two changes: personal names replaced by `<<placeholders>>`, and the third-person pronouns for the worked example's author replaced by "the author" (the repository convention for referring to a model author).

One phrase is deliberate and should survive reuse: "no dynamical field **at start**". That is the state M8 was accepted in, not a verdict on the model. Supplying the equation of motion is the program itself, and the author is working through it now: the [field-dynamics pre-registration](../openwave/xperiments/m8_mit/research/findings/m8_2_preregistration.md) is locked and was the author's own first pull request, with the Lagrangian-family survey ([M8.4](../openwave/xperiments/m8_mit/research/tasks/m8_4_task_details.md)) as the decisive science behind it.

```text
Hi <<model_author>>, here is a detailed description and instructions for you, I hope it helps:

OpenWave is an open-source platform, maintained by volunteers, for testing candidate field-theoretic models of matter: frameworks that claim to build particles, forces, and waves out of an underlying field, written down as a closed set of equations. It tests a model by running it: the model author supplies the equations, and the workflow turns them into a lattice implementation whose behavior is then measured rather than asserted.

The loop <<fellow_author>> described is the standard one: implement the Lagrangian in a GPU solver, relax an ansatz by energy minimization to find which stable field configurations the model actually has, then run dynamics on those configurations (kick them, collide them, boost them) and read off observables: energies, spectra, dispersion relations, decay products. What the model paper says could exist, the lattice shows does or does not exist, and how it behaves. The page <<fellow_author>> linked is exactly one of the research notes such a run produces.

The results land in a shared comparison table, MODELS.md: each framework is a column scored against the same criteria (rows covering particles, forces, waves, and quantum emergence), and every cell is earned by a runnable script plus a research note, or it stays marked "not yet tested". The bar is reproducibility, not orthodoxy, and a documented negative counts as much as a positive.

The platform is not one engine with one substrate. M5 is a 4x4 real symmetric tensor field, M6 is two coupled vector fields, M8 is spectral geometry on a fixed quotient manifold with no dynamical field at start. Each model gets its own directory, its own solver, and its own roadmap, so your framework does not have to look like any of the existing ones. What the platform cannot supply is the equations; that part stays yours, as the model author.

On <<fellow_author>>'s M9 suggestion: yes, the slot is open and the invitation is now formal. The application is one discussion post in the "New Model" category, no form; a maintainer picks it up from there. The one page to read is ONBOARDING_MODELS.md, and the most efficient thing you can do is point your AI assistant at it. It covers four steps: how to drive the work with an AI coding agent and what a model author owns, a self-evaluation to run on your own framework before anything else, the application itself, and what gets scaffolded if it lands.

One part worth stating plainly, because it decides what you can plan: this is a bring-your-own-compute offer. The runs behind your column are author-driven, on your own AI tokens and hardware. The platform supplies the shared criteria, the scaffold, the review, the standards, and the cross-model prior art; the column carries your name, so when someone challenges a claim on it, you answer.

For what onboarding a model to OpenWave looks like in practice: M8-MIT is the live worked example, one week in. <<another_model_author>> (model author) applied with a single discussion post and the column was scaffolded for the author. The first thing that happened was a certification gate: two blind agents independently re-derived the headline eigenvalue theorem of the bedrock paper to 10-digit precision without being shown the printed constants. Since then the author contributes through the normal fork-and-PR flow (two pull requests merged in the first week). One cell of the column is earned so far and the rest still reads "not yet tested"; that honesty is exactly what makes the earned cell worth something.

One more page worth knowing, because it shapes how everything gets done there: AI_HYGIENE.md, the working contract for AI-assisted research in the repo. The short version: an AI model's output is a draft, never a result, until something that is not a language model confirms it: a runnable script, a measurement, or the person who holds the authority to say so.

If you decide to try it, the natural first step costs nothing but a prompt: give your AI assistant your paper and the onboarding page above, and ask it to walk you through the STEP 1 self-evaluation. Whatever comes out of that is already the raw material for the application.

Hope this helps,
```

The message closed with the sender's standard signature, which lists the tool stack the maintainer side runs (IDE, coding agent, terminal, intelligence models, session hooks, Python, Taichi, Git). Worth keeping in an invitation: it makes the bring-your-own-compute paragraph concrete by showing what the other side of it looks like in practice.

## See also

| Doc | For |
| --- | --- |
| [`ONBOARDING_MODELS.md`](../ONBOARDING_MODELS.md) | The path this message invites into: STEP 0 drive it with an AI agent, STEP 1 self-evaluation, STEP 2 application, STEP 3 scaffold and first PR |
| [`MODELS.md`](../MODELS.md) | The comparison table an invited column joins, and the shared criteria it is scored against |
| [`AI_HYGIENE.md`](../AI_HYGIENE.md) | The AI-collaboration contract the message summarizes in one sentence |
| [`CONTRIBUTING.md`](../CONTRIBUTING.md) | Setup, fork and pull-request flow, DCO sign-off: the mechanics an author meets right after onboarding |
| [`m8_mit/__M8_model_briefing.md`](../openwave/xperiments/m8_mit/__M8_model_briefing.md) | The worked example the message cites, and what a scaffolded column looks like from outside |
