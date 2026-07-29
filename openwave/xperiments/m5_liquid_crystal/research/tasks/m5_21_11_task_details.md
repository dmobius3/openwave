# M5.21.11: the realistic-parameter bridge + unit calibration (the mass-ratio read)

**Status**: 🚧 PLANNED STUB (2026-07-19, the electron-status gap review). Lineage: this is the former **M5.9.0** ("Duda delta/g calibration + unit-scale prep"), re-scoped and renumbered into the M5.21 series the same day: it consumes the series outputs (the census levels + the fixed-J state) and delivers the lepton mass-ratio read, so it belongs to the hunt. The old record stays as archive: [`m5_9_0_task_details.md`](m5_9_0_task_details.md). Full PLAN at go.

## Scope (stub level)

| Piece | Content | Notes |
| --- | --- | --- |
| The regime walk | Ladder the toy results toward the physical regime (δ ~ 1e-10, g ~ 1e10, from the author's paper hints g⁴ ~ ke²/Gm² ≈ 1e38, δ² ~ ħc): measure at affordable rungs, verify each follows an analytic law, extrapolate the law | Direct simulation is out of lattice reach (the author-flagged numerics obstacle, [Q33](../m5_question_tracker.md#q33-detail); specialist contacts activated 2026-07-19) |
| Existing rungs | The M5.21.1 P4 asymptotic laws (the author's 2-equal-vortex + 3-equal-core structure exact in the physical limit); the [M5.21.8](m5_21_8_task_details.md) g-ladder (dressed-minimum ratio 0.82-0.84 stable across g = 8-64, tracking the 1/g law); the m\* formula verified to 0.009% | The bridge pattern is already demonstrated; this task systematizes it |
| The unit half | Anchor lattice numbers to real units: the M5.16 c₂ = αħc/64π Coulomb lock, the 511 keV Faber anchor, the ω ∝ m scale-covariance (#220) | The Coulomb-unit + LdG-to-rest-energy axis from the M5.9.0 lineage |
| The target read | The 1 : 206.8 : 3477 mass ratios (and the Koide check) on the census levels A < C < B (toy ratios 1 : 4.2 : 16.0) | The sharpest falsifiable number the lepton hunt can produce |
| Cross-column consumer | The M8 column's lepton-hierarchy cross-check (M8.6) is GATED on this task: its [readiness note](../../../m8_mit/research/findings/m8_6_readiness_note.md) § 8 states what the bridge has to deliver before that comparison is admissible | A frozen functional-level spec (parameters, spacing, couplings) derived independently of the lepton target, then a NEW census AT those parameters: no transform of the existing lattice energies, since a free `E_physical = f(E_lattice)` map would be a three-point mass fit |

**Gated by**: the M5.21 core results (the census levels + the [M5.21.9](m5_21_9_task_details.md) fixed-J state) + [Q33](../m5_question_tracker.md#q33-detail) (the specialist contacts) + user "go".

## OPEN DESIGN QUESTION (2026-07-29): does the extrapolated law count as "a census at the physical parameters"?

Raised when the M8.6 gate was cross-linked (row above), by reading the two specs against each other. **Settle this at PLAN, before any rung is measured**, because it decides what the target read is allowed to be.

| Side | What it requires |
| --- | --- |
| [M8.6 § 8](../../../m8_mit/research/findings/m8_6_readiness_note.md) condition 3 | the comparison uses stationary-state energies from a NEW census run AT the frozen physical parameters, directly; no state-by-state map or free nonlinear transformation of the three existing lattice energies (5.2611, 22.059, 84.085) |
| This task's regime walk (§ Scope) | measure at AFFORDABLE rungs, verify each follows an analytic law, extrapolate the law, because [Q33](../m5_question_tracker.md#q33-detail) says a run at δ ~ 1e-10, g ~ 1e10 is out of lattice reach |
| This task's target read (§ Scope) | the mass ratios "on the census levels A < C < B", which read literally is the transformation condition 3 rules out |

The tension is real and not a wording slip: condition 3 asks for a run at parameters Q33 says cannot be run, so taken literally it is unsatisfiable by any route M5 has. Three ways out, to be chosen and pre-registered rather than discovered mid-run:

| Option | What it commits to |
| --- | --- |
| **A. Extrapolated law as the census** (likely) | the ladder IS the measurement: each rung is a genuine census run, the law is fitted on rungs and validated out-of-sample, and the physical point is quoted with an extrapolation uncertainty. Requires the law's functional form and the rung set frozen BEFORE the ratios are looked at, or the fit becomes the three-point mass fit condition 3 exists to prevent |
| B. Narrow the claim | report the physical-regime read as an M5-internal result and leave M8.6 gated, with the cross-column comparison explicitly out of reach |
| C. Renegotiate condition 3 | take option A to the M8 column and get the amended condition written into the readiness note before it is relied on. The condition is the M8 author's, so this is author-gated either way |

Whichever is chosen, two constraints from § 8 hold unchanged and cost nothing to honor: the `A→e, C→μ, B→τ` assignment stays frozen from the pre-existing stability/decay rationale, never re-derived from a mass match, and no `1 : 5.9 : 15.1` Yukawa-derived figure enters the derivation or the calibration.
