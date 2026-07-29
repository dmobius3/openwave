# M4 / EWT, THEORY CANONICAL (spec of record)

> **Status: skeleton.** The sections below are the shape a spec of record takes on this
> platform, with a note in each on what belongs there. M4's theory record exists, but it
> is spread across the source library, the M3 equation docs and the briefing (§ 0 maps
> it). Consolidating it here is a job for the people who work on the model, not one a
> maintainer should do by paraphrase.
>
> **What this file is for once it is filled.** It becomes the document that wins when
> other docs disagree: the equations M4 actually implements, the constants and where each
> comes from, what the model claims, and what is assumed rather than derived. Everything
> else in the column cites it. The M8 equivalent,
> [`m8_theory_canonical.md`](../../m8_mit/research/m8_theory_canonical.md), is the
> reference for tone and depth.

## 0. Where the M4 record lives today

Until the sections below are written, these are the sources of record. Anything stated
here must trace to one of them.

| Source | Holds |
| --- | --- |
| [`../theory/_CITATIONS.md`](../theory/_CITATIONS.md) | the source registry: every EWT document and contributed paper behind this column, with DOIs where they exist. The papers themselves are local-only, gitignored, obtained from their venues |
| [`../__M4_model_briefing.md`](../__M4_model_briefing.md) | the column at a glance: model profile, free parameters, honest per-criterion status |
| [`../../m3_wolff_lafreniere/research/0a_equations.md`](../../m3_wolff_lafreniere/research/0a_equations.md) | the shared wave-structure equation set M4 inherits |
| [`../../m3_wolff_lafreniere/research/0_STATUS.md`](../../m3_wolff_lafreniere/research/0_STATUS.md) | the honest-blockers status record for the EWT lineage |
| [`M4_engine_upgrade.md`](M4_engine_upgrade.md) | the implementation record: what the vector-PDE engine does, phase by phase, with its validation evidence |
| [`M4_k_selectivity_Formalization.md`](M4_k_selectivity_Formalization.md) | the contributed nonlinear-stabilisation formalization |
| [`../wave_engine.py`](../wave_engine.py), [`../medium.py`](../medium.py), [`../force_motion.py`](../force_motion.py) | the code of record: what is actually solved, as opposed to what is written down |

## 1. The medium and the wave equation

What belongs here: the field, its domain, and the evolution equation M4 integrates,
written once, in the form the code implements. State the discretization separately from
the continuum equation, and say which one a given claim is about.

## 2. The potential and its variants

What belongs here: `V(ψ)` and each implemented mode, with the parameter that selects it,
the regime each is meant for, and the failure mode each was introduced to fix. A mode
that collapses or fails to arrest is documented here, not removed.

## 3. Wave centers, phase, and the particle construction

What belongs here: how a particle is built out of wave centers, what fixes their count
and arrangement, and what phase assignment means physically. Mark clearly which parts are
imposed by construction and which emerge from evolution: that boundary is the single most
important thing this document records.

## 4. Constants and their provenance

What belongs here: every constant the model uses, with its value, its source, and its
kind. A constant is one of: measured input, derived from other constants, or fitted. A
fitted constant is not a defect, but an undeclared one is.

## 5. What is assumed rather than derived

What belongs here: the honest ledger. The briefing already names several, including
charge sign imposed by a phase offset rather than emerging. Each entry names the
assumption, why it is currently needed, and what would remove it.

## 6. Consumption rules (standing)

What belongs here: how other columns and other documents may cite M4. The M8 version
states, for example, which numbers may be quoted as established and which carry a
tension. Until this section exists, cite the briefing's status icons, which are kept
honest per criterion.

## OPEN QUESTIONS

Seeded from the briefing's Status and Roadmap tables. Each is a candidate task for
[`m4_roadmap.md`](m4_roadmap.md); give it a stable `Qn` id when this section is adopted,
and keep resolved entries in place with their resolution, so citations to them stay
valid.

| # | Question | Where it stands |
| --- | --- | --- |
| Q1 | Is there a localized, stable soliton core in the vector PDE, and what arrests the collapse? | open, the column's headline problem: pure cubic collapses, saturating quintic does not arrest at CFL dt |
| Q2 | Is any wave-center count K selected, or are all K equally stable at perfect placement? | open, and the subject of contributed formalization and numerical work |
| Q3 | Can charge quantization emerge rather than be imposed through a phase offset? | open, imposed by construction today |
| Q4 | Does the L / T divergence-curl split produce the Coulomb far field? | open, not implemented |
| Q5 | Do the Dirichlet boundary reflections change any soliton-search conclusion, and does an absorbing boundary remove the effect? | open, a known contaminant of current results |

Route each by who or what can answer it: a script, the model's author, or an experiment.
The routing discipline is in
[`m8_agent_orientation.md § 3`](../../m8_mit/research/m8_agent_orientation.md).
