# METHOD NOTE: the audit page standard

> **Standing rule (2026-07-03, mandatory; renamed from "method report" 2026-07-04).** Every report, summary doc, or email addressed to a model's theory owner (M5: Dr. Jarek Duda; M7: Marc Fleury; any external physicist or advisor) carries a method note built to this standard. The bar: the reader must be able to **audit the result by reading**, without trusting us and without reverse-engineering Python. A result whose method cannot be found is not a result to that audience.
>
> Why a NOTE, not a report: the method note is a **component, not the report itself**. Sometimes it ships attached to a results report alongside findings; sometimes it stands alone as the internal audit page (the maintainer re-verifying what an instrument computes, months later, is the same reader as the theory owner). Either way it is the page where the equations and their code locations live.

## The motivating incident (2026-07-03, M5.16)

The M5.16 calibration report led with results and buried the physics: the energy functional and potential lived inside a ~1000-line script's docstring, and the report's inspection list pointed at a driver file containing no physics. The theory owner read the code and replied, verbatim: *"this median energy radius requires choice of potential (to be found) - I don't see anywhere there, and this minimized energy is integral of Hamiltonian - I also don't see ... so honestly still I have no idea what does it calculate."* None of the five questions attached to the report were answered. The failure was **legibility of the verification surface, not correctness**: his same reply restated a prescription that largely matched what the task had implemented. Full record: [`m5_4h_convo_2026.07.03.md`](../openwave/xperiments/m5_liquid_crystal/research/tasks/m5_4h_convo_2026.07.03.md).

The lesson generalizes: physicists verify by reading the Lagrangian / Hamiltonian / potential in mathematical notation and checking that the code implements it term by term. A rendered simulation, a results table, or a well-written prose summary is not a verification surface for this audience.

## The required shape (in this order)

| # | Element | Requirement |
| --- | --- | --- |
| 1 | **Equations first** | The report OPENS with the physics in math notation: the Lagrangian / Hamiltonian (or free-energy) density, the potential with its explicit form and coefficients, the ansatz / boundary conditions, the discretization, the minimization or evolution scheme, and the DEFINITION of every reported observable. All before any result. |
| 2 | **Equation-to-code map** | A table mapping every term of the functional (and every observable) to `function name` + `script.py:line` + a clickable GitHub **permalink** (commit-pinned `https://github.com/.../blob/<commit>/...#L<n>-L<m>`, so line anchors never drift). Test: any term findable in under a minute by a human clicking links. |
| 3 | **Small auditable physics module** | The energy functional / Hamiltonian lives in a single-purpose module (target ~200 lines, not a monolith), whose docstring carries the same equations as the report. Drivers, sweeps, and plotting import it; they never re-implement physics. |
| 4 | **Results after methods** | Each result states its pre-registered gate and convergence evidence next to the number. |
| 5 | **Minimal inspection set** | At most ~4 artifacts, ordered physics-first: the functional module FIRST, the driver LAST. Never send a reader into a driver to find physics. |
| 6 | **Scope honesty** | State explicitly what was NOT computed (observables, regimes, anchors), in the owner's own vocabulary for the model. |

## Pre-send checklist

| Check | Pass condition |
| --- | --- |
| Can the owner find `V(...)` and the Hamiltonian density from the report in one click? | yes, via the equation-to-code map |
| Does the report open with equations, not results? | yes |
| Is every code link clickable (absolute https) with a file:line anchor? | yes. `blob/main` links are the DEFAULT for task-scoped files (frozen after task close, anchors stable; The maintainer 2026-07-05); commit-pin ONLY when linking into live/evolving files whose lines will drift |
| Is the functional in a small module whose docstring matches the report's equations? | yes |
| Is each headline number next to its gate + convergence evidence? | yes |
| Is the not-computed list explicit? | yes |
| REQUIRED (upgraded from RECOMMENDED 2026-07-05): has an independent second agent ADVERSARIALLY audited the note (own script, different method/seed, tried to REFUTE each claim, re-walked the equation-to-code map)? | yes, before sending; adopted from the model owner's own advice ("careful small steps, maybe multiple agents verifying each other", Duda 2026-07-03). Upgraded after the audit's first real catch (M5.18: a witness overclaim refuted + a missed vacuum-branch structure found, `m5_18_verification_note.md § 10`). Record the audit outcome IN the note (the § 10 pattern): the audit record is part of the rigor surface |

## Scope

Applies to every model in `openwave/xperiments/` when reporting outward to its theory owner or any external physicist. It does not constrain internal working notes; it constrains anything that crosses the repo boundary toward a scientist.
