# M5.23.1, fixed-J isorotation dynamics → production port (the former M5.24 round 4)

> 🚧 PLANNED STUB (2026-07-19, staged at the [M5.24](m5_24_task_details.md) close, user-approved re-home of its round-4 scope). Renumbered M5.26 → M5.23.1 (2026-07-23, the rendering-series consolidation: M5.23.x = the official rendering series). Roadmap row: [`m5_roadmap.md § RENDERING TASKS`](../m5_roadmap.md). The full `## TASK PLANNING` lands at "go".

**Why**: [M5.24](m5_24_task_details.md) brought the launcher to the verified-L canonical stack, and the live demo now shows (correctly) that free evolution never spins up the ZBW clock: the two-stack consensus ([M5.21.3](m5_21_3_task_details.md) + [M5.21.8](m5_21_8_task_details.md)) is that the clock exists only as a FIXED-J isorotation state (ω* = J/2kin, measured clock inertia kin ≈ 0.119). Once [M5.21.9](m5_21_9_task_details.md) builds + verifies that state research-side (with the author's Larmor protocol as the acceptance observable), this task ports the validated fixed-J evolution into the production engines: the same catch-up pattern as M5.24 (taichi-first, per-gap selftests against the research reference, launcher wiring, live verification).

**Scope sketch (to be firmed at go)**:

| Step | Deliverable |
| --- | --- |
| 1 | Port the M5.21.9-validated fixed-J machinery (the constraint-carried J evolution) into `engine2_pde`, gated against the research scripts |
| 2 | Launcher wiring on the canonical path (a fixed-J xperiment config; the RELAX → set-J → EVOLVE flow) |
| 3 | Live verification: the stable rotating state on screen = the first honest ZBW δ-sweep (simulated dynamics only, per the standing no-display-only-kinematics directive, [`m5_visualization.md`](../m5_visualization.md)) |

**Gating**: user "go" ([M5.21.9](m5_21_9_task_details.md) ✅; the physics-first hold cleared at the [M5.21.4](m5_21_task_details.md) close 2026-07-21; the former RENDERING UNLOCK marker retired at the 2026-07-23 Backlog reorder, row order = run sequence). First row of the sequenced Backlog. Feeds [M5.23.2](m5_23_2_task_details.md) (the J/μ twist demo arm rides this port).

## Consumes from M5.21.9 (wired 2026-07-20 at the run close)

| Input | Where it lives |
| --- | --- |
| The fixed-J states (three J rungs, all holding) | `research/data/m5_21_9_fixedj_om{0.2,0.5,1}_end.npz` (local, manifest-listed; regen `python3 m5_21_9_d_fixedj.py om=<v> maxit=1200 refresh=300`, ~6 min each). ω\* = J/(2kin) per state: the port's live clock spins at THESE measured rates, nothing display-only |
| The certified 4×4 leapfrog | `m5_21_9_e_larmor.py leap()` (the M5.21.6 form lifted to 4×4: velocity masked every kick, implicit γ): E-conservation 2.2e-8 per 400 steps, dt = 0.02 certified post-audit (4× margin): the research reference the port's per-gap selftests run against |
| The clock thermodynamics | dE/dJ = ω\* at ~1% ([`../findings/m5_21_9_note.md § 7`](../findings/m5_21_9_note.md)): the port's energy ledger must reproduce the Legendre closure as a selftest |
| The Larmor round-2 protocol | The J-flip discriminator + the modeled-floor ±B pair ladder (note § 6): long windows are native here; whatever the M5.21.9 addendum leaves open of the linear read, this task inherits the measurement-grade ladder |
| ⚠️ the kin-convention flag | 0.297 (probe variant, this run) vs 0.1206 (conjugation-tangent, adopted at M5.21.3; the stub's "kin ≈ 0.119" above is that convention): the port must PIN one convention in code and document it (audit CL9; absolute J and ħ/2 numbers depend on the factor 2.46). **QUANTIFIED at [M5.21.5](m5_21_5_task_details.md) (2026-07-21)**: the m5_21_9 carried J (unit-Frobenius flow parametrization, 0.048 at the ω\* = 0.2 rung) is NOT the physical angular momentum; the PHYSICAL-RATE S_env(twist) = 37.2 on the same state (ratio 771). The port pins the physical-rate convention ([`../findings/m5_21_5_note.md § 5`](../findings/m5_21_5_note.md)) for any g or ħ/2 statement |

## Round 3 carries Q36 (wired 2026-07-21)

The author's reply to the 21.9 note poses the quadratic-effect origin question ([Q36](../m5_question_tracker.md#q36-detail); [`m5_21_convo.md § 2026-07-21 03:30`](m5_21_convo.md)): the measured instant-on systematic (ours) vs a formulation-truncation effect ("3x3 formulation (4x4 might change)") vs a genuine quadratic coupling ("they might search experimentally"). Round 3 (adiabatic ramp-on + body-frame read) is the discriminator for the first; its design should also record which functional terms the port's field coupling carries, so the truncation hypothesis is answerable by inspection rather than a new run.
