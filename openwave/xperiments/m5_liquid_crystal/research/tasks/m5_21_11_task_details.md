# M5.21.11: the realistic-parameter bridge + unit calibration (the mass-ratio read)

**Status**: 🚧 PLANNED STUB (2026-07-19, the electron-status gap review). Lineage: this is the former **M5.9.0** ("Duda delta/g calibration + unit-scale prep"), re-scoped and renumbered into the M5.21 series the same day: it consumes the series outputs (the census levels + the fixed-J state) and delivers the lepton mass-ratio read, so it belongs to the hunt. The old record stays as archive: [`m5_9_0_task_details.md`](m5_9_0_task_details.md). Full PLAN at go.

## Scope (stub level)

| Piece | Content | Notes |
| --- | --- | --- |
| The regime walk | Ladder the toy results toward the physical regime (δ ~ 1e-10, g ~ 1e10, from the author's paper hints g⁴ ~ ke²/Gm² ≈ 1e38, δ² ~ ħc): measure at affordable rungs, verify each follows an analytic law, extrapolate the law | Direct simulation is out of lattice reach (the author-flagged numerics obstacle, [Q33](../m5_question_tracker.md#q33-detail); specialist contacts activated 2026-07-19) |
| Existing rungs | The M5.21.1 P4 asymptotic laws (the author's 2-equal-vortex + 3-equal-core structure exact in the physical limit); the [M5.21.8](m5_21_8_task_details.md) g-ladder (dressed-minimum ratio 0.82-0.84 stable across g = 8-64, tracking the 1/g law); the m\* formula verified to 0.009% | The bridge pattern is already demonstrated; this task systematizes it |
| The unit half | Anchor lattice numbers to real units: the M5.16 c₂ = αħc/64π Coulomb lock, the 511 keV Faber anchor, the ω ∝ m scale-covariance (#220) | The Coulomb-unit + LdG-to-rest-energy axis from the M5.9.0 lineage |
| The target read | The 1 : 206.8 : 3477 mass ratios (and the Koide check) on the census levels A < C < B (toy ratios 1 : 4.2 : 16.0) | The sharpest falsifiable number the lepton hunt can produce |
| Cross-column consumer | The M8 column's lepton-hierarchy cross-check (M8.6) is GATED on this task: its [readiness note](../../../m8_mit/research/findings/m8_6_readiness_note.md) § 8 states what the bridge has to deliver before that comparison is admissible | A frozen functional-level spec (parameters, spacing, couplings) derived independently of the lepton target, then a NEW census: either directly at those parameters (route a) or extrapolated from a fresh pre-registered ladder (route b, § 8 amended 2026-07-29). No transform of the existing lattice energies either way, since a free `E_physical = f(E_lattice)` map would be a three-point mass fit |

**Gated by**: the M5.21 core results (the census levels + the [M5.21.9](m5_21_9_task_details.md) fixed-J state) + [Q33](../m5_question_tracker.md#q33-detail) (the specialist contacts) + user "go".

## DESIGN QUESTION, RESOLVED (2026-07-29): does the extrapolated law count as "a census at the physical parameters"?

Raised when the M8.6 gate was cross-linked (row above), by reading the two specs against each other. **Resolved the same day by the M8 author**, who amended the condition rather than reinterpreting it: [readiness note § 8](../../../m8_mit/research/findings/m8_6_readiness_note.md) now admits an extrapolation route explicitly, in the author's own summary `"A scientifically, C procedurally"`. The question and the options are kept below, because the amended condition has to be read against them.

### Why it was raised

| Side | What it required |
| --- | --- |
| [M8.6 § 8](../../../m8_mit/research/findings/m8_6_readiness_note.md) condition 3, as originally written | the comparison uses stationary-state energies from a NEW census run AT the frozen physical parameters, directly; no state-by-state map or free nonlinear transformation of the three existing lattice energies (5.2611, 22.059, 84.085) |
| This task's regime walk (§ Scope) | measure at AFFORDABLE rungs, verify each follows an analytic law, extrapolate the law, because [Q33](../m5_question_tracker.md#q33-detail) says a run at δ ~ 1e-10, g ~ 1e10 is out of lattice reach |
| This task's target read (§ Scope) | the mass ratios "on the census levels A < C < B", which read literally is the transformation condition 3 ruled out |

The tension was real and not a wording slip: condition 3 as written asked for a run at parameters Q33 says cannot be run, so taken literally it was unsatisfiable by any route M5 has.

| Option | What it committed to | Outcome |
| --- | --- | --- |
| **A. Extrapolated law as the census** | the ladder IS the measurement: each rung is a genuine census run, the law is fitted on rungs and validated out-of-sample, and the physical point is quoted with an extrapolation uncertainty | ✅ taken, on the science |
| B. Narrow the claim | report the physical-regime read as an M5-internal result and leave M8.6 gated, with the cross-column comparison explicitly out of reach | ❌ not taken |
| **C. Renegotiate condition 3** | get the amended condition written into the readiness note before it is relied on, rather than reading extrapolation into the old wording | ✅ taken, on the procedure |

One correction to option A as first stated here: it asked for the law and rung set frozen "before the ratios are looked at", and the amendment deliberately rejects that as unachievable blindness, since the charged-lepton ratios are public and known to everyone involved. The enforceable discipline is **no-refit pre-registration**, not blindness. That is the stronger requirement, and this task adopts it.

### What the amended condition 3 requires of this task

Route (b), the extrapolation route, is admitted under guardrails. These are PLAN preconditions: each has to be settled and frozen BEFORE any rung is measured, or the route closes.

| Precondition | Requirement |
| --- | --- |
| Asymptotic form | derived from M5-side theory independently of the lepton target, never chosen for its fit quality |
| Frozen framework | rung set, fitting procedure, holdout tests, branch-tracking rules, and uncertainty model all fixed before any performance evaluation against the charged-lepton ratios is run |
| Uniform application | the same frozen framework applied to all three branches (A, C, B), no per-branch tuning |
| Barred inputs | the three existing toy energies (5.2611, 22.059, 84.085) may not enter the new fit as data points, an exponent search, or a post-hoc transformation |
| Uncertainty bar | a usable uncertainty on `E_C/E_A` and `E_B/E_A`, not merely a stable ordering (condition 4) |
| Discretization term | inside that uncertainty model, a per-rung discretization term: grid refinement of ALL THREE branches on a pre-registered subset of rungs spanning the ladder, plus a frozen rule propagating the term to unrefined rungs, carried alongside the extrapolation error (condition 4 as amended, [#378](https://github.com/openwave-labs/openwave/pull/378)) |
| Failure is terminal | if the pre-registered scaling law or the holdout gates fail, M8.6 stays gated. No second framework |
| Claim ceiling | route (b)'s output is a model-based physical-regime PREDICTION, not a directly simulated physical census, and it cannot make `1 : 5.9 : 15.1` independent evidence |

Two constraints from § 8 hold unchanged and cost nothing to honor: the `A→e, C→μ, B→τ` assignment stays frozen from the pre-existing stability/decay rationale, never re-derived from a mass match, and no `1 : 5.9 : 15.1` Yukawa-derived figure enters the derivation or the calibration.

**The uncertainty gap, raised in review of the amendment and closed the same day.** Condition 4 as first amended covered extrapolation error along the ladder but not discretization error at each rung. The readiness note's own § 5 records `E_A` moving `4.920 → 5.261 → 5.921` across three grid rungs at fixed δ (~20%), with B (2.61) and C (2.35) exceeding the consistency bar where A does not, so branch-dependent error of that size cannot be assumed to cancel in a ratio. [#378](https://github.com/openwave-labs/openwave/pull/378) closed the gap, going further than the review asked: all three branches refined, a pre-registered subset spanning the ladder, and a frozen propagation rule (the `Discretization term` row above).

**One scoping consequence for this task's PLAN.** Refining a rung costs roughly 8× that rung in 3D, and the ladder's top rung sits at the affordability ceiling by construction, so a refinement subset read as literally including the top rung meets the same affordability wall that made the original condition 3 unsatisfiable. The satisfiable reading is the propagation rule: measure the convergence order where three resolutions are affordable for all three branches, then propagate it upward to the unrefined rungs. Budget the ladder so that reading holds, and if the wording needs pinning down, raise it as a one-clause amendment before the ladder is frozen, not after.
