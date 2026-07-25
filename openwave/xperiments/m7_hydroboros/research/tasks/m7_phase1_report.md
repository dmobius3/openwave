# M7 / HydroBoros - PHASE 1 REPORT

  ![hydroboros_icon](../images/hydroboros_icon_small.jpg)

> The human-readable summary of Phase 1 (tasks M7.0-M7.7, finished on 2026-07-04). Everything here is a pointer: each claim links its detail doc, script, and data. Written for the two physics parents (Marc Fleury, Paul Werbos) and for anyone auditing the program.

## 1. What M7 is

**HydroBoros** fuses Fleury's toroidal-EM electron (analytic, [arXiv:2510.22384](https://arxiv.org/abs/2510.22384)) with Werbos's Ouroboros self-confinement (M6, 1D radial) into one full-3D-PDE lattice program: the electron's field configuration is both *specified* (a self-linked toroidal Beltrami vortex) and *earned* (as the constrained energy-minimizer of the written functional). One-page orientation: [`__M7_model_briefing.md`](../../__M7_model_briefing.md). Full background (the two parents, the mathematical boundaries, the dynamics frame): [`m7_background.md`](../m7_background.md).

## 2. What Phase 1 did (M7.0 → M7.7)

| Task | What landed | Detail |
| --- | --- | --- |
| M7.0 bootstrap | 64-doc theory corpus + planning scaffolding | [`m7_0_bootstrap.md`](m7_0_bootstrap.md) |
| M7.1 infra | 3D harmonic lattice, Taichi-AD gradients (2.3e-15 vs complex-step), FIRE + constraints; Woltjer-Taylor known-answer gate (`λ → 2π/L` at 5.5e-6) | [`m7_1_infra.md`](m7_1_infra.md) |
| M7.2 reproduce Fleury | the FLDB torus quadrature to closed forms at 1.4e-4; the printed solution reconstructed digit-for-digit; the corrected-convention `U = 0.958 m_ec²` (Q10) | [`m7_2_fleury_torus.md`](m7_2_fleury_torus.md) |
| M7.3 reproduce M6 in 3D | the benchmark ODE pinned as a verbatim EL reduction (sympy scan); the M6 charged ledger `H/Q` reproduced in full 3D to **4.7e-5**; three discoveries: the ledger is window-defined (Q11), the `§ 2.2` form is not an EL reduction (Q12), the M6 electron is a 3D constrained saddle (Q13) | [`m7_3_ouroboros_3d.md`](m7_3_ouroboros_3d.md) |
| M7.4 the charged soliton (new physics) | **the program's first stable finite-size 3D soliton family** (Taylor-dressed, `E = 0.802\|H_A\|` universal); helicity measurably load-bearing; the written `f = gs²` beats the benchmark signs (Q6 resolved empirically) | [`m7_4_charged_soliton.md`](m7_4_charged_soliton.md) |
| M7.5 clock + stability | translation exactness (`⟨E_real⟩ = E_ω` to 1.85e-14); **the vacuum tachyon discovered** (measured rate 0.785 vs analytic 0.786) + the existence threshold **`ω* = 0.786`** (the clock IS the stabilizer); the ω-`Q_can` Legendre conjugacy verified to 1-2% (Q14 opened) | [`m7_5_clock_stability.md`](m7_5_clock_stability.md) |
| M7.6 observables | **the rotating electron**: a clean `j_z = 1` per-quantum wave (0.6%); **Coulomb landed** (Gauss 99.1%, far field −2.14, two-charge `1/d` reference-matched, a measured **1.17 coupling dressing**); the RKKY-style neutral exchange discovered | [`m7_6_observables.md`](m7_6_observables.md) |
| M7.7 consolidation (milestone) | the canonical spec + a ~200-line physics module + a one-script reproduction (**all gates first-try at two resolutions; engine-vs-reference 1.4e-14**); the 21-cell column drafted + staged | [`m7_7_consolidation.md`](m7_7_consolidation.md) |

## 3. The canonical spec (the equations)

**[`m7_theory_canonical.md`](../m7_theory_canonical.md)** is the results-of-record spec: it OPENS with the Lagrangian and every pinned convention (with the provenance of each pin), then the equation-to-code map, then results after methods, then the explicit not-computed list. The physics itself lives in one auditable ~200-line module, [`scripts/m7_functional.py`](../scripts/m7_functional.py) (its docstring carries the same equations), and the whole electron reproduces from one script, [`scripts/m7_7_canonical.py`](../scripts/m7_7_canonical.py): ~3 min quick mode, ~10 min record mode, with the Taichi engine cross-validated against the reference module to 1.4e-14 as a first-class gate.

## 4. Results and findings (the staged column)

The 21-cell coverage column was drafted at M7.7: **0 ✅ / 8 ⚠️ / 0 ❌ / 13 🚧**, every partial carrying its named caveat. It was deliberately **staged, not yet published to [`MODELS.md`](../../../../../MODELS.md)** at that point (published at M7.21, 2026-07-18): the program still carries its own open theory question (Q14, the vacuum tachyon) and convention questions (Q15, Q12; Q10 resolved 2026-07-04 by the author, corrected figure 0.95 vs our 0.958) that move icons either way; **promotion comes after Phase 2, gated by the question tables below** (the M7.21 publication milestone, Phase 3). What the eight partials say, in three plots:

**The soliton exists and is grid-convergent.** Three helicity-carrying seeds relax to one Taylor-dressed family; the winner is localized (r50 = 3.4), approximately Beltrami in the core (`λ_eff ≈ −1` near-constant), dilation-stable (a measured constrained-Derrick interior minimum), with the J-condensate co-located; both zero-helicity parent seeds die, so helicity is measurably the guard:

![the M7.4 winner sections](../plots/m7_4_winner_sections.png)

**The clock is load-bearing, quantitatively.** The linearized vacuum of the truncation is unconditionally tachyonic at long wavelength (left: both quartic conventions, `det M(0) = −1`); the measured vacuum growth rate matches the analytic band edge to 0.15% with no amplitude threshold (middle); and harmonic solitons exist only above `ω* = 0.786`, exactly the band edge (right: runaway at ω = 0.70/0.75, clean solitons from 0.79 up):

![the M7.5 tachyon + threshold triptych](../plots/m7_5_tachyon_scan.png)

**Coulomb is real once a monopole exists.** With a fixed charge reservoir (`j₀`) and the scalar potential unfrozen, Gauss's law closes at 99.1% (left), the far field is `1/r²` (middle, slope −2.14 vs the chargeless −3.7), and the two-charge splitting matches the Poisson reference in the same grounded box at a constant **1.17 ± 0.02** dressing (right), the program's first fine-structure-flavored number:

![the M7.6 Coulomb triptych](../plots/m7_6_coulomb.png)

## 5. Open questions for Marc Fleury, most critical first

Tracked with stable IDs in [`m7_question_tracker.md`](../m7_question_tracker.md) (full statement + history at each anchor).

> **Updated 2026-07-06, post-call.** Q10 ✅ resolved by the author's reply (the Eq 122/124/127 bug confirmed; corrected figure 0.95 vs our 0.958). **Q15 ✅ resolved by author directive at the 2026-07-06 call**: no frequency mapping pinned ("the frequency is emergent"); the target is the **observable `S_z = ℏ/2`**, measured by [M7.8](m7_8_helicity_pair.md). **M7.8 measured it (2026-07-07)**: the pair asymmetry reads **one quantum, not ℏ/2**, at every stable state of the fixed-`(Q_can, H_A)` frame (the two-mode pair is expelled toward single-helicity; adversarially audited); the ℏ/2 route shifts to the frequency mapping ([walkthrough § 7.1](m7_phase1_walkthrough.md)). The remaining asks (Q7, Q3) ride the author's post-August static→dynamic Beltrami window; Q1/Q4 as below. Priority order matches the [tracker](../m7_question_tracker.md).

| ID | Question | Why we ask |
| --- | --- | --- |
| [Q7](../m7_question_tracker.md#q7-detail) | The charge-carrying construction, reframed by the notes' **no-go theorem**: a localized oscillating interior carries no net charge; the Coulomb tail lives exactly where `λ → 0`, and closure = the turning-point matching at `r = w` plus a **tail-strength quantization** (what fixes `e`?) | your "start Trkalian, take off the training wheels" is validated (approximately-Beltrami cores, `\|align\| = 0.96`; charge measured, never imposed) and the fixed-reservoir prescription landed Coulomb, which reads as an imposed tail; the matching prescription and the quantization condition are the pieces we cannot pick ourselves |
| [Q3](../m7_question_tracker.md#q3-detail) | Divergence charge vs linking charge: we measured them **independent** (linking gates existence; the RMS charge coexists, unslaved); the notes agree structurally (charge rides the λ-gradient, helicity rides the winding) and propose the **tail-quantization conjecture** as the deeper slaving | decides whether charge quantization should be sought in topology (then we design that test) or in the source sector |
| [Q1](../m7_question_tracker.md#q1-detail) | Substrate reading: the `(A_μ, J_μ)` doublet vs a single Riemann-Silberstein field; target manifold S² (Pisello) vs S³ (Faber)? (the notes reconfirm A-primary + `φ = 0`, silent on the manifold) | fixes which topological invariants Phase 2 should constrain (the knot-sector work needs the right manifold) |
| [Q4](../m7_question_tracker.md#q4-detail) | Mostly delivered: the 2026-07-05 working notes ARE the promised Beltrami material; the residual ask is the Enciso / Peralta-Salas contact status | the rigidity/existence theorems shaped the whole relax-and-measure design; more of that literature de-risks Phase 2 |

## 6. Open questions for Paul Werbos, most critical first

> **Updated 2026-07-05:** the first round-trip brought no question resolutions, so the list was re-scoped. **Q14 / Q11 / Q13 moved to SELF-DETERMINE**: they are physics/measurement questions the lattice answers in-model during Phase 2, and the results get shared rather than asked. The remaining asks are **Q12 + Q9**, both author-gated by nature (intent + definitions) and micro-scoped (one smallest unit of authorial input per contact).

The two asks:

| ID | Question | Why we ask |
| --- | --- | --- |
| [Q12](../m7_question_tracker.md#q12-detail) | ODE provenance (+ the Q6 residual): the benchmark ODE is the verbatim EL of the same-phase doublet with FOCUSING signs, while the written `f = gs²` (the stable branch, which we run) differs; and `0d_canonical § 2.2`'s `2ωα` form matches no EL reduction we scanned: which ansatz + ODE is canonical, and do you confirm the benchmark sign slip + the `(g/4)` transcription fix? | author-gated (intent); your confirmation closes the convention file and fixes our M6 doc; nothing is blocked on it (we proceed on the empirically stable branch) |
| [Q9](../m7_question_tracker.md#q9-detail) | The `(Ω, G)` dictionary of the v5 scan: the micro-ask is **one verified electron-island point** in `(ω, g, λ, m_J)` | author-gated (definitions); our measured threshold `ω* = 0.786` is a candidate anchor, and the Q11 micro-slice added a constraint: a decaying (`k > 0`) island is impossible at the canonical `(λ, m_c)`, so a consistent dictionary must move `λ` and/or `m_J` |

Self-determined during Phase 2 (results will be shared as they land):

| ID | Question | Self-route |
| --- | --- | --- |
| [Q14](../m7_question_tracker.md#q14-detail) | **The vacuum tachyon**: `det M(0) = −1` for ANY `f` (measured 0.785 vs analytic 0.786; no `β*` threshold in the vector truncation); the gift: the harmonic frame survives exactly because `ω > ω* = 0.786`, **the de Broglie clock is the vacuum-stability mechanism** | linearize WITH the scalar/Gauss sector on our lattice (does it lift `det M(0)`?) + test the variable-λ mass profile; the M7.18 pre-step (real-time runs stay blocked until cured) |
| [Q11](../m7_question_tracker.md#q11-detail) | The charged `H/Q = 1.6890` is window-defined; does a genuinely localized charged branch exist? **2026-07-05 micro-slice** ([script](../scripts/m7_13_q11_window.py)): window-divergent for every `(ω, g)` tested, and the localization boundary at canonical `(λ, m_c)` is exactly `ω < ω*`, so localization and existence are mutually exclusive there | the branch hunt scans `(λ, m_c)`, not `(ω, g)`; rides M7.19 |
| [Q13](../m7_question_tracker.md#q13-detail) | The M6 electron is a 3D constrained **saddle**; helicity stabilizes the same torus (our M7.4 fix); solitons exist only above `ω*` | the physics is measured in-model; the residual 1D-conjugate-point provenance question rides with Q12 as a share/confirm item |

## 7. What comes next (Phase 2)

Full roadmap: [`m7_roadmap.md`](../m7_roadmap.md), restructured after the 2026-07-06 Phase-1-review call: **Phase 1 extends** through the author-collaboration track before Phase 2 opens.

| Next | Task | One line |
| --- | --- | --- |
| 1 | [M7.8](m7_8_helicity_pair.md) helicity-pair 3:1 test + the Phase 1 walkthrough | **GO (2026-07-06 call)**: the independent lattice check of the closure notes (`U₊/U₋` vs `3 + α/2 + 4f_bb`) + the pair-asymmetry spin vs the ℏ/2 directive (Q15) + the under-the-hood [walkthrough](m7_phase1_walkthrough.md) + the canonical-spec refresh; the package goes to the author for an adversarial external pass |
| 2 | [M7.9](m7_9_chaosbook.md) ChaosBook benchmark | canonical dynamical-systems exercises with published solutions, implemented + unit-tested; the orbit-hunting toolkit; the self-test the author asked for |
| 3 | [M7.10](m7_10_pure_maxwell.md) pure-Maxwell sector: the no-Lagrangian test | coupling terms off, plain Maxwell evolution: Theorem 2 (Trkalian cavity mode) as a known-answer gate, free-space evaporation + electron destruction time (the honest boundary as a measurement), the coupling ladder (tachyon attribution); uses the M7.9 orbit toolkit; results join M7.8 + M7.9 in the one [walkthrough](m7_phase1_walkthrough.md) report for the author |
| 4 | M7.11-M7.14 (reserved) | the rest of the pure-Maxwell / Beltrami orbit track (variable-λ seeds, the orbit hunt, static→dynamic charge), firming up with the author's post-August re-entry |
| 5 | Phase 2: M7.15-M7.20 | magnetic, gravity, nuclear, antimatter + annihilation (Q14 pre-step), leptons + neutrinos (the Q11 `(λ, m_c)` scan), dark matter |

After Phase 2: **M7.21 publishes the column to MODELS.md** for cross-model benchmarking (Phase 3, including the second canonical-spec refresh), gated by the question tables above; then composites (M7.22, Phase 4) and production rendering (M7.23, Phase 5).

---

Cross-refs: briefing [`__M7_model_briefing.md`](../../__M7_model_briefing.md) · background [`m7_background.md`](../m7_background.md) · spec [`m7_theory_canonical.md`](../m7_theory_canonical.md) · published column [`MODELS.md`](../../../../../MODELS.md) · tracker [`m7_question_tracker.md`](../m7_question_tracker.md) · roadmap [`m7_roadmap.md`](../m7_roadmap.md) · the milestone task [`m7_7_consolidation.md`](m7_7_consolidation.md).
