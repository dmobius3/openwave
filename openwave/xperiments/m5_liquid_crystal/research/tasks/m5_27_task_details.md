# M5.27: the background-scalar time sector (the staged 4×4), entrainment pilot

**Status**: 🚧 PROPOSED 2026-07-21, **pending author go**. Staged in the roadmap Backlog above the RENDERING UNLOCK marker for visibility; it does NOT hold the rendering block (M5.23.1 → M5.23.2 run on user "go" regardless). Roadmap: [`m5_roadmap.md § BACKLOG`](../m5_roadmap.md). Offered as the THIRD branch of the fork left open at the [M5.21.3](m5_21_3_task_details.md) close (the fixed-J reading vs the T2η potential lift vs this), proposed to the author in the M5.21-series close round (sent 2026-07-21 17:06 EDT; [convo record](m5_21_convo.md)).

## 1. WHY NOW: the M5.21-series close diagnosis

The M5.21 electron series closed with a diagnostic pattern.

- Everything static-topological closes ✅ (charge quantization, spin-½, the annihilation ledgers, the lepton census),
- and everything that needs the TIME sector fails to close ❌  (free-clock existence, μ, g, force-law dynamics).

The Lagrangian has been frozen since the verified-L era by design; M5.17 → M5.21.4 was instrument-building, audits, and honest re-grading on a fixed L. The record now points at exactly one missing organ: the time sector.

| Task | Headline outcome | Model-evolution read |
| --- | --- | --- |
| [M5.21.3](m5_21_3_task_details.md) (the 4D lift, closed 07-18) | THE SADDLE: all 24 time-mixing curvatures negative (time derivatives WANT to be nonzero), BUT free 4×4 minimization finds NO stable dynamical electron: every rotation velocity costs energy, no finite stationary ω, statics wins outright. The clock survives only as a fixed-J CONSTRAINED state (ω\* = J/2kin, J imposed by hand). [Note](../findings/m5_21_3_note.md) | The first true 4×4 dynamics attempt came back existence-negative |
| [M5.21.5](m5_21_5_task_details.md) (μ + g, closed 07-21) | μ channel exists but is a parity-cancellation residue spanning 4 orders across preparations; g NO closure (8.5e-4 to 1.45); the canonical g = 1.97 RETRO-FLAGGED. [Note](../findings/m5_21_5_note.md) | The observable that needs the clock degraded |
| [M5.21.4](m5_21_4_task_details.md) (2-particle, closed 07-21) | Conduit annihilation mechanism audited (genuinely new, the antimatter row delivered); but NO static pair regime exists, Coulomb 1/d NOT confirmed at reachable d, string tension ansatz-grade. [Note](../findings/m5_21_4_note.md) | Mixed: dynamics instrument matured, force law open |

## 2. THE PROPOSAL

| Element | Content |
| --- | --- |
| Fields | M(x,t) 3×3 traceless symmetric, T2 term set, verified-L UNCHANGED (charge quantization, census, annihilation all inherited) + ONE real scalar χ(x,t) |
| Vacuum | χ = A·cos(ω̄t), spatially uniform: a vacuum rest-frame background oscillation |
| Coupling | Promote the g-eigenvalue of the author's own spectrum `D = diag(g, 1, δ, 0)` to `g + κχ(x,t)`: the 4×4's time-time entry made dynamical FIRST, the (0,i) mixed block deferred. A STAGED 4×4, not an alternative to the author's 4×4 selection |
| Economy | +1 field component vs +4 (M₀₀ + three M₀ᵢ); a standard sign-definite kinetic term, while the all-negative time-mixing curvatures ([M5.21.3](m5_21_3_task_details.md)) and the divergent regions of the author's analytical notebook live in the deferred mixed block (hypothesis, checkable against [M5.21.8](m5_21_8_task_details.md)) |
| First observables | Lock-or-not (Arnold tongue in κ, ω̄); J and μ re-read under lock (the [M5.21.5](m5_21_5_task_details.md) protocol re-run); two-defect force under a shared background (the in-phase attraction channel); the drive-power ledger |

The economy argument in one line: promoting M to a 4×4 spacetime tensor adds ~4 awkwardly-signed components; a scalar background adds 1. The 4×4 route's persistent trouble may simply mean the time sector wants to be a SEPARATE FIELD, not extra indices on the order parameter.

## 3. DESIGN DECISION: ENTRAIN, not REPLACE

| Argument | Why it lands on ENTRAIN |
| --- | --- |
| Energy honesty | At lock, an injection-locked oscillator draws ZERO average power from the drive; the soliton stays self-sustaining and Derrick escape is not secretly pumping. REPLACE makes rest energy drive-dependent: fake stability, broken ledger |
| The record | M5.8's self-starting clock is a ✅ in-platform result; REPLACE would contradict the column's own record |
| The M5.21.3 twist | On the verified-L stack the free clock does NOT self-start (rotation costs energy). The regime question is genuinely open, which makes it the FIRST pre-registered observable of this pilot, not an assumption |
| The decisive detail | The fixed-J constraint is a hand-imposed stand-in for exactly what a background drive supplies dynamically. If lock delivers a J-budget without imposing J, the constrained-state crutch disappears |

Pre-registered outcomes (all three are results): (a) an intrinsic mode locks to the background = ENTRAIN confirmed; (b) oscillation exists only while driven, with measured drive power = replace-like, the ledger flags it; (c) no coupling = null.

## 4. THE BACKGROUND ENERGY LEDGER (declared)

The background field's energy is booked with the classical energy density of an oscillating medium: `E/V = ½ρω²A²` (equivalently `ρ(fA)²` up to convention factors), quadratic in both amplitude and frequency. This is declared, not derived: it is the textbook harmonic form (acoustics, elasticity, cavity modes), and the pilot's job is to measure what the coupled system does with it, not to justify it.

| Reference anchor | Role in the ledger |
| --- | --- |
| Classical oscillating-medium energy density (Rayleigh; ½ρω²A² per mode) | The bookkeeping form itself; standard for any harmonic background |
| de Broglie's internal clock (1924): `m₀c² = hν₀` | The target identification the pilot tests: rest energy ↔ lock frequency, mass = coupling strength to the background mode; the column's ZBW clock re-reads are the measured side of the same identification |
| Planck-Einstein `E = ħω` | The quantum limit the entrained ledger must stay consistent with when the realistic-parameter bridge is crossed ([Q33](../m5_question_tracker.md#q33-detail) lineage); not used at toy parameters |
| Adler injection-locking power relation | The drive-power instrument: zero average power at lock, nonzero off-lock; the § 3 regime discriminator is read directly off this ledger line |

Discipline: all quantities in M5 program units; no physical constants enter the pilot; the three new knobs (κ, A, ω̄) come pre-registered with search spaces (§ 6), nothing tuned per observable.

## 5. CHANCES BY SECTOR (pre-run assessment)

The structural reason the import is attractive: each side's validated strength covers the other's recorded negative. M5's open sectors (free-clock existence, hierarchy origin, the 4×4 route, gravity) are exactly where a background wave has known mechanisms.

| Sector | What the background brings | Known-physics anchor | Chance |
| --- | --- | --- | --- |
| Clock | Entrainment of the defect's oscillation to a universal reference; mass = coupling strength to the background mode | Driven-oscillator entrainment (Arnold tongues); M5.8 measured an intrinsic clock, so this is locking an existing oscillator, the easy case | HIGH |
| Particle stability | The drive supplies the time-periodicity M5's Derrick escape needs, externally sustained: parametric stabilization of otherwise-unstable configurations | Kapitza pendulum, Paul traps, ponderomotive trapping; all textbook | HIGH |
| Propulsion / de Broglie wave | A boosted standing wave Doppler-splits into carrier + modulation, which is how de Broglie derived the phase wave (1924) | Couder-Bush walking droplets: a particle propelled and guided by its own background standing wave, pilot-wave dynamics emerging | MEDIUM-HIGH |
| Gravity (weak field, clock sector) | Two coinciding mechanisms: (a) a defect near-field depresses the local background frequency, neighboring clocks slow, refractive attraction; (b) secondary Bjerknes force: oscillators phase-locked to the SAME background attract when in phase, amplitude-proportional (mass-proportional), universally attractive | Bjerknes forces between bubbles in ultrasound (measured, classical); analog-gravity refractive metrics | MEDIUM-HIGH for Newtonian + redshift |
| Gravity (full GR grade) | Scalar-background gravity historically fails light bending by the factor 2 (Nordström); recovering Schwarzschild needs the moving-medium metric (Gordon / Painlevé-Gullstrand flow), not a static index | Analog-gravity literature | LOW-MEDIUM without extra structure |
| μ + magnetic curl | Circulation of the near-field modification around the defect axis; could supply a first-principles bridge behind the currently underived K = 4/α factor | Vortex-carrying acoustic/optical near fields | MEDIUM (double-counting risk: M5 already has a μ channel) |
| Spin ½ | Little to add: the apolar mechanism already closed this ✅ machine-exact; a vector background could add photon-like helicity | | LOW need, keep M5's |
| Lepton hierarchy | A cavity background has a DISCRETE harmonic ladder; defects locking to different harmonics gives a discrete mass ladder, a mechanism M5 currently lacks entirely | LaFreniere mode numbers; the M8 column's spectral ladder | SPECULATIVE, but a mechanism where M5 has none |

## 6. HAZARDS

| Hazard | Why it bites |
| --- | --- |
| Preferred frame | A universal standing wave defines a rest frame; Lorentz invariance must emerge from Doppler-consistent wave kinematics. The deepest theoretical risk; needs a pre-registered boost test early |
| Energy bookkeeping | If defects draw stability from the background, the background is an energy reservoir; drive input vs soliton energy must be booked honestly or "stability" is just pumping |
| Replace vs entrain | Decided upfront (§ 3): ENTRAIN. An undecided hybrid would double-count the clock energy |
| Parameter growth | κ, A, ω̄ are new knobs; platform standards apply (derive or pre-register with search spaces, nothing tuned per observable) |

## 7. ANTICIPATED OBJECTIONS (named upfront, answers on record)

| Objection | Answer |
| --- | --- |
| "This abandons the selected 4×4 route" | No: it is a STAGED 4×4 (the time-time entry first, the 0-row later) and the third branch of the fork the author left open at the M5.21.3 close, motivated by the author's own results (the saddle + the notebook divergences) |
| Preferred frame | Named in § 6, boost test pre-registered; de Broglie's 1924 phase-wave construction is the anchor for covariance emerging from a boosted oscillation |
| "Another field, another knob" | The § 2 economy row + § 6 parameter discipline: κ, A, ω̄ pre-registered, no per-observable tuning |
| "This imports oscillation theory" | The mechanism anchors (de Broglie, injection locking, Kapitza stabilization) are model-independent classical oscillator physics; the energy form is textbook. The convergence, if it appears, is platform triangulation, not adoption (§ 4) |

## 8. RELATION TO THE PROGRAM + SEQUENCING

| Item | Relation |
| --- | --- |
| The M5.21.3-close fork | This is the third branch (fixed-J reading vs T2η lift vs background scalar); the author's call |
| [M5.21.8](m5_21_8_task_details.md) (notebook verification) | Natural companion: the § 2 economy hypothesis (the pathology lives in the mixed block) is checkable against its results |
| RENDERING UNLOCK / [M5.23.1](m5_23_1_task_details.md) | Rendering proceeds regardless. Note: M5.23.1 ports the fixed-J clock into production; if this pilot delivers a lock, the clock's origin story changes (a dynamically supplied J-budget vs the hand-imposed constraint), so the pilot tells us whether the port target is the final clock or the crutch |
| Gate | Pending author go; nothing runs before it |

## 9. DEFINITION OF DONE (pilot scope, draft; frozen at go)

| # | Item |
| --- | --- |
| 1 | Pre-registration lock: κ, A, ω̄ search spaces + the § 3 three-outcome table + the boost test protocol, frozen BEFORE numerics |
| 2 | Gates: 3D regression exact on the verified-L stack at κ = 0; χ-sector energy conservation; the drive-power ledger instrumented |
| 3 | The lock scan: defect response vs (κ, ω̄), Arnold-tongue map, regime verdict per § 3 |
| 4 | J and μ re-read under lock (the M5.21.5 protocol) + the two-defect force under a shared background |
| 5 | Adversarial audit (independent implementation) before anything author-facing; method note per the house standard |
