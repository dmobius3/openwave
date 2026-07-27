# M8 / MIT, ROADMAP

> **Scaffold baseline (2026-07-21).** The M8 column was scaffolded by the maintainers
> from the author's onboarding proposal
> ([discussion #312](https://github.com/openwave-labs/openwave/discussions/312)). The
> program is a **field-dynamics collaboration**: MIT supplies the arena (S³/2I + the
> Möbius edge) and the target spectrum (the McKay ladder); the platform supplies
> Lagrangian candidates, simulation engineering, and grading standards. The spec of
> record: [`m8_theory_canonical.md`](m8_theory_canonical.md) (FIRST READ; canonical
> when docs disagree). Rationale and gap map: [`m8_background.md`](m8_background.md);
> cross-model reading map: [`m8_platform_pointers.md`](m8_platform_pointers.md). AI
> agents bootstrap via [`m8_agent_orientation.md`](m8_agent_orientation.md) ("read the
> m8_agent_orientation.md"), then run tasks with "go task m8.<n>".
>
> **Mode of work.** Research mode FIRST: headless scripts + research notes + plots, no
> GUI (the suggested per-task layout: `tasks/m8_<id>_task_details.md`, scripts / data /
> plots under `research/` with `m8_<id>_` prefixes). The 3D rendering port (M5-style
> interactive launcher) is a LATER stage, gated on field dynamics validating
> in-platform: M8.7 below, pointers in
> [`m8_platform_pointers.md § 7`](m8_platform_pointers.md).
>
> **Standing rules (platform-standard, from day one):** pre-register gates and
> conventions BEFORE each run; no calibrated conventions (derive and pre-register, or
> record as a fit with its search space); honest negatives are results; substantive
> claims get an adversarial audit ([`AI_HYGIENE.md`](../../../../AI_HYGIENE.md) § 1);
> author-facing reports follow
> [`dev_docs/METHOD_NOTE.md`](../../../../dev_docs/METHOD_NOTE.md) (equations first +
> equation-to-code map). Ownership is marked per task: the author drives the column;
> maintainer-run tasks are labeled.
>
> **Self-checks must be able to fail (M8.2, 2026-07-27).** Any line a script prints as
> PASS should be mutation-tested before it ships: change the thing it checks to something
> wrong and confirm the check goes red. A check whose two sides evaluate the same
> expression always passes and reads, to a later maintainer, exactly like a verified
> result. M8.2 shipped one (the coexact table's trivial column compared against the very
> rule that produced it) and it was replaced at merge by two that fail under mutation:
> `dims` against the affine E8 mark condition, and `dist` against BFS on the McKay graph.
> Where a table has no independent published target to check against, label it ASSERTED
> rather than giving it a self-check that cannot discriminate.
>
> **Borrowing other columns' families (platform ruling, 2026-07-24).** The M8 program
> tests candidate families drawn from M4 / M5 / M7 inside the MIT arena, so it runs under
> [`dev_docs/CROSS_MODEL_TESTING.md`](../../../../dev_docs/CROSS_MODEL_TESTING.md): a
> borrowed family is native and untwisted unless its author declares an internal
> representation or a soldering prescription is pre-registered; "not applicable" is a
> neutral status; a soldered family is scored as its own object; author silence is a
> valid terminal state and never blocks a pre-registration lock. Author-gated questions
> go from author to author directly (Q&A discussion), not through a maintainer relay.

## IN PROGRESS

| TaskID | Title | Description | Gated By |
| --- | --- | --- | --- |
| (none) | | M8.1 ✅ (2026-07-21) and M8.2 ✅ (2026-07-27) are both closed; the column is the author's to drive: M8.3, M8.5 and M8.6 are all startable, and M8.4 now waits only on M8.5 | |

## BACKLOG

| TaskID | Title | Description | Owner | Gated By |
| --- | --- | --- | --- | --- |
| [M8.1.1](tasks/m8_1_1_task_details.md) | Second blind run: the remaining bedrock theorems | Blind independent verification (the M8.1 protocol reused) of the two bedrock papers the author shared on [#312](https://github.com/openwave-labs/openwave/discussions/312#discussioncomment-17758091) (2026-07-24, "gaps and asymetry on S3": [SSRN 6968698](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=6968698) + [SSRN 7129118](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=7129118)). Maintainer-run BY CONSTRUCTION (the author cannot blind-verify the author's own theorems). Trigger: any analytic number that becomes a pre-registered target of the dynamics program is blind-verified before it can bear a target. The M8.2 lock left this open deliberately (its § 2 timing guard: a T-b2 verification landing before the first T-b2-eligible run enters through a newly signed module version, one landing after is a separate prospective test). Does NOT gate M8.3 / M8.6; runs at maintainer pace | maintainers | (none; scheduling only) |
| [M8.3](tasks/m8_3_task_details.md) | Mass-formula reproducer script | `m = μ_Λ · C_geom · (√Ω)^(dist/30) · T²` with every constant (McKay distance, Reidemeister torsion, C_geom weight) recomputed from its own definition, never quoted; PDG comparison; residuals reported at the ledger's weight. Grades the analytic sector the way the platform scores EWT's analytic masses | author | (none) |
| [M8.4](tasks/m8_4_task_details.md) | Lagrangian-family survey on S³/2I | The central question: can a nonlinear field equation on S³/2I have topological-defect or standing-wave solutions whose energies realize the McKay slot structure? Candidates drawn from the platform's columns (M5 Landau-de Gennes + Frank, M4 nonlinear scalar/vector wave, M7 two-vector; pointer map § 2), written on the compact quotient with the background clock; includes the Derrick analysis on a compact arena (the scaling argument changes: R sets a scale) and the anti-periodic (double-cover) sector. **Why this survey is Lagrangian-framed:** because the target is *energies*. A mass ladder needs an energy functional to compare slots against, and the Derrick analysis is itself a scaling argument about such a functional. That is a constraint imposed by the target, not a platform requirement: the admission bar is a closed set of equations, and tier-2 dynamics with no action behind it is admissible ([`ONBOARDING_MODELS.md § the four tiers`](../../../../ONBOARDING_MODELS.md#step-1-self-evaluation)) | author, platform support | M8.2 ✅ (2026-07-27), M8.5 |
| [M8.5](tasks/m8_5_task_details.md) | Quotient-manifold simulation engineering | The real infrastructure work: simulate on S³/2I via (a) a 2I-equivariant grid on S³ (identification maps) or (b) a spectral method in 2I-symmetric harmonics (the basis IS the representation theory). Prototype both far enough to pick one; document trade-offs (pointer map § 6) | author, platform support | M8.2 ✅ (2026-07-27) |
| [M8.6](tasks/m8_6_task_details.md) | The McKay-distance rule vs M5's measured lepton hierarchy | Bounded cross-check, no simulation needed: does the McKay-distance rule reproduce M5's measured eigenvalue hierarchy (1 : 5.9 : 15.1, the open "hierarchy origin" of the M5 lepton row) and/or the physical mass ratios? Pre-registered mapping BEFORE computing; either outcome closes a live question in TWO columns | joint | (none) |

## LATER (gated)

| TaskID | Title | Description | Gate |
| --- | --- | --- | --- |
| [M8.7](tasks/m8_7_task_details.md) | The 3D rendering port | Port the M5-style interactive stack (per-model `_launcher.py` + `engine1-4` split + the shared GGUI rendering in `openwave/i_o/`) so the validated S³/2I dynamics runs as a live demo. Do NOT start this before the gate: rendering an unvalidated dynamics showcases nothing. Porting instructions for AI agents: [`m8_platform_pointers.md § 7`](m8_platform_pointers.md) | field dynamics validated in-platform (an M8.4-lineage result with an audited method note) |

## STATUS AT A GLANCE (2026-07-27)

| Question | Answer |
| --- | --- |
| Where is M8? | Scaffold merged AND the certification gate passed: M8.1 ✅ (2026-07-21) verified the arena's headline eigenvalue theorem at 10-digit precision (blind two-agent, audited), flipping the first MODELS.md cell (gravity → ⚠️, count 1 ⚠️ / 20 🚧). M8.2 ✅ (2026-07-27) then locked the field-dynamics pre-registration, the author's first contribution through the normal fork → branch → PR flow ([`CONTRIBUTING.md`](../../../../CONTRIBUTING.md), DCO). M8.3, M8.5 and M8.6 are startable; M8.4 waits on M8.5 |
| What kind of column is it? | The platform's first top-down structural model: strong on the origin of the numbers (representation theory on S³/2I), absent on dynamics. The M8 program exists to supply the dynamics half |
| What decides the program? | M8.1 ✅ certified the arena's headline eigenvalue (the certification-gate philosophy paid off: both blind agents re-derived the paper's constants 2/e and −4e^(−2γ) to 10 digits without seeing them); M8.2 ✅ then froze what M8.4 will be graded against, before any numerics existed to tune it toward; M8.4 is the decisive science (does ANY reasonable Lagrangian on S³/2I realize the McKay slot structure?); M8.6 pays off in both the M8 and M5 columns regardless of sign |
| Evidence discipline | The author's own claim ledger is adopted as the grading baseline (structural results = the core; the numeric mass table = low weight, capped by the author's own pre-registered nulls); platform standards (pre-registration, adversarial audit, method notes) apply from day one |

## DONE

> Newly done tasks APPEND at the end.

| TaskID | Title | Description | Completed |
| --- | --- | --- | --- |
| (scaffold) | The M8 column scaffold | Onboarding evaluation (provenance + § 1 / § 4 / § 5.1 checks) passed; [`__M8_model_briefing.md`](../__M8_model_briefing.md) + [`m8_background.md`](m8_background.md) + [`m8_platform_pointers.md`](m8_platform_pointers.md) + [`../theory/_CITATIONS.md`](../theory/_CITATIONS.md) written; MODELS.md column added at 21 🚧 | 2026-07-21 |
| [M8.1](tasks/m8_1_task_details.md) | THE CERTIFICATION GATE: independent eigensolve of the twisted Möbius Laplacian | ✅ PASSED, all five pre-registered claims CONFIRMED: λ₁⁺ = 2/R² (narrow) and α₀(α₀+1)/R² (wide) exact across two mutually independent blind implementations (≤ 1.4e-8 agreement, two extra blind W points); the extension-stability threshold bisected blind to 2R/e (1e-12) and the bridging defect-state coefficient extrapolated blind to −4e^(−2γ) (10 digits); Friedrichs zero mode exact, exactly one bridging bound state; audit 6/6 fidelity checks, one summary-wording defect (AF-1) dispositioned; the auditor independently re-derived the paper's Legendre ladder structure. Canonical § 3 + briefing + MODELS.md gravity cell (🚧 → ⚠️) synced. [Method note](findings/m8_1_method_note.md) | 2026-07-21 |
| [M8.2](tasks/m8_2_task_details.md) | Pre-registration lock for the field-dynamics program | ✅ LOCKED: [`findings/m8_2_preregistration.md`](findings/m8_2_preregistration.md), a MODULAR contract (immutable core §§ 1-5 + per-family modules § 6 + signed per-family execution appendices § 7, which carry the numerics and gate M8.4 rather than this lock). Targets = the structural ladder only (the 24-entry numeric table stays OUT); four-axis outcome language; no-search rule; provenance pinned at `c9dc3796` / `1ec2cd97` + SHA-256 source hashes, pin audit verified against upstream. The index question became the platform standing rule [`CROSS_MODEL_TESTING.md`](../../../../dev_docs/CROSS_MODEL_TESTING.md), under which all three native families record T-c "not applicable" (neutral), with @jeffsyee ([#333](https://github.com/openwave-labs/openwave/discussions/333), M4 geometric displacement) and @JarekDuda ([#334](https://github.com/openwave-labs/openwave/discussions/334), M5 level-3 spacetime tensor) both declaring in confirmation of the default; "M5 + P" and `M7_ad` stay M8-owned objects. First-occurrence table independently reproduced 9/9 by the maintainer (explicit-quaternion + Burnside, no shared method). Author-contributed via [PR #350](https://github.com/openwave-labs/openwave/pull/350); maintainer edits at merge announced in that thread. Remaining and NOT part of this lock: the native quotient operators (M4, M5, M7) and the Zenodo-deposit byte-check, both § 7 appendix work | 2026-07-27 |
