# M7.11: the variable-λ vacuum-cure check (Q14)

> **Status: PLANNED** (defined 2026-07-07 at the post-M7.10 park; first task of the reserved Maxwell band M7.11-M7.14). m7 waits until the author returns from the conference and reviews the [walkthrough § 7](m7_phase1_walkthrough.md) package; an adversarial pass may re-scope this plan before it runs. Roadmap row: [`m7_roadmap.md`](../m7_roadmap.md) § BACKLOG (top = next task).

## 0. Planning notes (2026-07-07, the post-M7.10 next-task analysis)

**Why this is the recommended next run.** It is exactly where M7.10 points (a Q14 cure must fix `k → 0` while keeping soliton-scale binding); it interpolates the two tracks (coupled inside, pure Maxwell outside); it unblocks M7.18 and every real-time program; and nonlinearity makes orbits ISOLATED, so the M7.9 `find_cycle` machinery finally hunts for real. The roadmap lists this work as the M7.18 pre-step; it runs here as the first Maxwell-band task because M7.10 bootstrapped it by design.

**The two-track scoreboard this task must thread** (the M7.10 closing finding): the harmonic/functional track has particles with a sick vacuum; the pure-Maxwell track has a healthy vacuum with no particles (no localization, no attractor, no isolated orbits: all three measured). The variable-λ construction is the natural candidate precisely because it is the interpolation.

**Sequencing at the park:** first action on the author's return = the maintainer sends the walkthrough package (the maintainer's voice; content bullets preparable on request); this task is the next "go" trigger after that, unless the author's pass re-scopes the band. Parallel options that do NOT need the Q14 cure: [M7.15](m7_15_magnetic.md) / [M7.19](m7_19_lepton_family.md) / [M7.20](m7_20_dark_matter.md) (harmonic-frame measurements); [M7.16](m7_16_gravity.md) any time (honest pass/fail); the M7.9 E6 book-walkthrough sessions are user-scheduled with the [reading digest](m7_9_reading_digest.md) and [benchmark report](m7_9_benchmark_report.md) ready.

## 1. Objective

Answer the Q14 cure question in-model, by measurement: **can the vacuum tachyon be fixed while keeping the soliton?** [M7.10](m7_10_pure_maxwell.md) sealed the constraint any cure must satisfy: the bilinear `A·J` term carries BOTH the tachyon (`det M(0) = −ε_x²`, rates on `0.786√ε` at 0.1-0.3%) and the confinement ((1,0) binds the most compact state; (0,1) binds nothing). A cure must therefore fix the `k → 0` sector while keeping the bilinear's binding at soliton scale. Two candidates are computable now:

| Candidate | The check |
| --- | --- |
| The scalar/Gauss sector | linearize the vacuum WITH the scalar components active (a₀/j₀ unfrozen, the full 16-component doublet): does the constraint sector lift `det M(0) < 0`? (the M7.4 § 1 finding says the naive averaged Hamiltonian is unstable in the timelike components: this check settles whether the full sector cures or worsens) |
| The variable-λ mass profile (the author's construction, 2026-07-05 notes) | a spatially-varying coupling `ε_x(r)` (equivalently `m_A(r) = ℏλ(r)/c`): canonical inside the soliton, `→ 0` outside (`λ(r)² = λ₀²(1 − r²/w²)` as the reference profile). Prediction to test: a localized bound interior + a massless, healthy, tachyon-free exterior vacuum |

## 2. Pre-registered experiments (draft; finalize at PLAN on go)

| # | Experiment | Gate / known answer |
| --- | --- | --- |
| E1 | scalar-sector linearization: the full-doublet analogue of the `M(k)` matrix, analytic + lattice noise-probe cross-check | the M7.5/M7.10 pattern: closed form first, measured growth rates against it (the eps-ladder machinery from [`m7_10_pure_maxwell.py`](../scripts/m7_10_pure_maxwell.py) reruns with scalars unfrozen) |
| E2 | variable-`ε_x(r)` engine (one new kernel argument: the coupling becomes a field) + vacuum probe OUTSIDE the profile | exterior noise does NOT grow (pure-Maxwell healthy vacuum, the M7.10 (0,0) anchor); interior supports the bound state |
| E3 | the relaxed state on the variable-λ background: does a localized minimizer exist at ω = 1 with a healthy exterior? | r50 at soliton scale (the M7.10 ladder's ε = 1 anchor, 3.4) AND zero exterior growth over ≥ 10 T: the needle threaded, or an honest failure mode documented |
| E4 | orbit isolation: `find_cycle` on the variable-λ flow from a perturbed period | with nonlinear/inhomogeneous structure the period equation is no longer degenerate (the M7.10 non-isolation demo is the null reference): genuine period recovery = the first isolated orbit of the program, the "island of superstability" made literal |

## 3. Unknowns pass (blindspots at definition; redo at go)

| Unknown | Route |
| --- | --- |
| What SUSTAINS `λ(r)` (back-reaction)? | nature-gated and out of scope here: this task tests the profile as a fixed background (the author's own framing); the self-consistency question is the M7.12+ physics and the author's post-August lead |
| Does the scalar sector help or hurt? | machine-checkable: the M7.4 § 1 timelike-instability finding cuts both ways; the linearization decides |
| The reference profile's edge (`λ → 0` at finite r) may be C⁰-sharp like the CK cutoff | anticipate the M7.10 sheet lesson: measure the discrete profile's own spectrum first; smooth the edge if the corner radiates |
| Author re-scope risk | the park exists for this: the plan is a draft until the author's walkthrough pass lands |

## 4. Cross-refs

[M7.10 findings](m7_10_pure_maxwell.md) (the constraint + the ladder machinery) · [tracker Q14 detail](../m7_question_tracker.md) (the cure-candidate history) · [roadmap](../m7_roadmap.md) (M7.18 lists this as its pre-step; it runs here as the first Maxwell-band task) · the author's variable-λ construction: `theory/electron_beltrami/` via [`theory/_CITATIONS.md`](../../theory/_CITATIONS.md).
