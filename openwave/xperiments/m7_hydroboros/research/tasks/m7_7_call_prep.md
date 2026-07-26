# M7.7 call prep, 2026-07-06: Phase 1 review + next steps

> **ARCHIVED (2026-07-06, post-call).** The call happened; every open item here was folded forward, so nothing needs to be read from this document anymore: outcomes in the [tracker chronology](../m7_question_tracker.md) (Q15 resolved-directive; Q7/Q3 riding the author's post-August window), the M7.8 GO + scope in [`m7_8_helicity_pair.md`](m7_8_helicity_pair.md) (§ 6 there maps every section of this sheet to its new home), the ChaosBook program in [`m7_9_chaosbook.md`](m7_9_chaosbook.md), and the restructured [roadmap](../m7_roadmap.md). Kept as the pre-call record.
>
> Agenda + reference sheet for the 2026-07-06 call. Nothing to prepare or answer in advance; everything below gets walked through live. Companion reading (optional): the one-page [Phase 1 report](m7_phase1_report.md).

## 1. Call goals

| # | Goal |
| --- | --- |
| A | Review Phase 1 (M7.0-M7.7) and the [Phase 1 report](m7_phase1_report.md) |
| B | Settle the **units contract** (Q15): `ω = ω_Compton` vs `ω = 2ω_C` (Zitter), the one choice the mass + spin readings need |
| C | Decide **run / no-run** on the proposed M7.8 lattice test of the closure-notes spin structure (§ 3 below) |
| D | Walk the remaining open questions (§ 6) |
| E | Agree the Phase 2 sequence ([roadmap](../m7_roadmap.md)) |

## 2. Marc's question, the 3:1 ratio, in two lines

The ratio is not a free claim; it is forced by the notes' own two closure postulates. With `U₊ + U₋ = ℏω` (one quantum, P2) and `(U₊ − U₋)/ω = ℏ/2` (spin, P3):

```text
U₊ = 3ℏω/4,   U₋ = ℏω/4   ⟹   U₊/U₋ = 3        (exactly, at δ = 0)
```

Deducting the charge-tail + backbone fraction `δ = α/8 + f_bb` from the pair budget shifts it at first order:

```text
U₊/U₋ = (3 − 2δ)/(1 − 2δ) ≈ 3 + α/2 + 4f_bb ≈ 3.004 + 4f_bb
```

So the open physics is not the arithmetic; it is whether the postulates P2 + P3 survive in a full nonlinear 3D theory. That is a measurable question, which is § 3.

## 3. Proposed decision: the M7.8 helicity-pair test

An independent, cross-formalism check: the notes predict the number analytically in the thin-torus ansatz; the lattice measures the same quantity by full 3D relaxation, with no ansatz assumed in the final state (seeds relax to the true minimizer, then observables are read off).

| Aspect | Content |
| --- | --- |
| Build | the repaired helicity-pair state from the notes: CK modes `(m, n, s) = (1, ±1, ±1)` with `A_r ≠ 0`, LG profile under the closures `λ₀σ = 2`, `w = λ₀σ²` |
| Relax | fixed `Q_can` (+ helicity), the same constrained-relaxation frame validated in Phase 1 (M7.3-M7.6) |
| Measure | the per-helicity energy split `U₊/U₋` |
| Prediction | `3 + α/2 + 4f_bb` (§ 2) |
| If it lands | the closure postulates survive the full nonlinear theory; the spin-½ reading gains an in-model measurement |
| If it fails | either a postulate or the thin-torus truncation breaks in 3D, and the failure mode says which |

It is queued as **M7.8, the first task of Phase 2** ([roadmap](../m7_roadmap.md)), pending this call: run / no-run + scope.

## 4. Systematic ansatz search in the Beltrami register: what the pipeline already does

The machinery for trying new ansatz families systematically exists and is validated; Phase 1 used the old ansatz only as **seeds**, never as assumed answers (states are relaxed to the functional's true minimizer, then measured).

| Capability | Where |
| --- | --- |
| Seed gallery, 6 families (ABC / Trkalian, CK spheromak, Bateman / Hopf, FLDB torus, M6 embedding, Ceperley rotating mode) | [`m7_1_infra.md`](m7_1_infra.md), gallery plot [`m7_1_seeder_gallery.png`](../plots/m7_1_seeder_gallery.png) |
| Constrained relaxation (fixed `Q_can` + helicity, Taichi-AD gradients) | [`m7_functional.py`](../scripts/m7_functional.py) (~200-line auditable physics module) |
| Observables battery (energy budget, `j_z` per quantum, spin sector split, Gauss / Coulomb, KG dispersion) | [`m7_6_observables.md`](m7_6_observables.md) |
| One-script reproduction of everything of record | [`m7_7_canonical.py`](../scripts/m7_7_canonical.py) (all gates, minutes at N = 48) |
| The canonical spec (equations first, equation-to-code map, commit-pinned) | [`m7_theory_canonical.md`](../m7_theory_canonical.md) |

The new CK/LG closure family is fully parameterized (`m, n, s, σ, w, R`), so it slots into the same sweep machinery; a chaotic-flow direction has a corpus anchor as well (Dombre et al. 1986, ABC chaotic streamlines, already in the source library).

## 5. The units contract (Q15), the one choice to pin

Full decision table: [`m7_theory_canonical.md § 4`](../m7_theory_canonical.md). Short form, using the measured `j_z = 1` per quantum and `ωL_z/E = 2.07`:

| Choice | Consequence |
| --- | --- |
| `ω = ω_Compton` (the notes' Case B: P1 + P2, tail budget excludes Case A; spin ℏ/2 read from the pair asymmetry `(U₊ − U₋)/ω`) | total field `L_z ≈ 2ℏ`; the `2ω_C` Zitterbewegung appears in the bilinears (`cos(2φ − 2ωt)`) |
| `ω = 2ω_C` (Zitter; our earlier on-record recommendation) | total `L_z ≈ ℏ` within 3%; survives in the notes only on the flagged half-quantum branch |

The two readings agree at the observable level (the bilinears oscillate at `2ω` either way); the choice pins the absolute mass anchor (`E = 6.3246` program units → `m_ec²`) and the ℏ/2-vs-ℏ spin cell. The call decides which is canonical.

## 6. Open questions to discuss (priority order)

Full statements + evidence: [question tracker](../m7_question_tracker.md).

| Q | One line | What moved on 2026-07-05 |
| --- | --- | --- |
| [Q15](../m7_question_tracker.md#q15-detail) | the units contract (§ 5) | the notes pin Case B; needs live confirmation |
| [Q7](../m7_question_tracker.md#q7-detail) | the charge construction | the no-go theorem relocates net charge to the exterior tail (`λ → 0` region); the discussion item is the turning-point matching at `r = w` and what quantizes the tail to exactly `e` |
| [Q3](../m7_question_tracker.md#q3-detail) | the two charges | we measured divergence charge and linking charge independent (M7.4); the notes make them different carriers (λ-gradient vs winding) and propose the tail-quantization condition as the deeper link |
| [Q1](../m7_question_tracker.md#q1-detail) | target manifold | S² (Pisello) vs S³ (Faber); the notes reconfirm A-primary + `φ = 0` but are silent here |
| [Q4](../m7_question_tracker.md#q4-detail) | source material | Enciso / Peralta-Salas: still worth pursuing? |

## 7. Already answered since the Phase 1 report

| Item | Status |
| --- | --- |
| Q10, the FLDB energy convention | ✅ Resolved by your reply: the Eq 122/124/127 bug confirmed; your corrected figure 0.95 vs our independently-computed 0.958 (the ~1% gap sits inside the thin-torus approximations). The lattice reproduction and the corrected convention are now the numbers of record. |
| Q4, the promised Beltrami material | ✅ Mostly delivered: the two working-note documents (ansatz recommendation + closure) are archived in the local source library ([manifest](../../theory/_CITATIONS.md), author copies, kept off the public repo). Remaining: the Enciso / Peralta-Salas question, logistics only. |
| The working notes themselves | **Independently audited on receipt, and the core math checks out.** We re-derived, by hand and separately from the notes: Theorem 1 (constant-λ fields are chargeless), Theorem 2 (oscillating Trkalian at `ω = cλ₀` solves free Maxwell exactly), the LG closure relations (`λ₀² − k² = 4/σ²`, `w = λ₀σ²`, `ε = 1/4` exact), the tail-energy integral (`U_tail = e²/8πε₀w = (α/8)m_ec²` at `w = 4λ̄_C`), and the corrected 3:1 algebra (§ 2). The assumed parts are cleanly separable: the three pinning postulates P1-P3, the harmonic `λ²(r)` profile, and the secular-tail reading. |
| New measurement (2026-07-05 night): the charged-profile window scan | The windowed charge **diverges** (`Q(2W)/Q(W) ≈ 2`, the radiation signature) for every `(ω, g)` combination tested, `ω ∈ {0.8, 1.0} × g ∈ {0, 0.25, 0.5, 1.0}`, windows 8-128 ([script](../scripts/m7_13_q11_window.py) · [plot](../plots/m7_13_q11_window.png)): the quartic coupling rescales amplitude but never localizes (it does not enter the far field). **Bonus finding:** the decaying-channel condition at the canonical couplings solves to `ω < ω* = 0.786`, the same golden-ratio number as the tachyon band edge and the existence threshold, so **localization and soliton existence are mutually exclusive** in the vector truncation there. Your variable-λ profile (the twist dying at the edge, leaving a massless exterior) is a natural candidate escape, which makes the § 6 tail-matching discussion (Q7) even more central. |

The notes also moved three of our open questions in useful ways, folded into § 6 above.

## 8. After the call

Phase 2 ([roadmap](../m7_roadmap.md), M7.15-M7.20): the helicity-pair test (pending § 3), then magnetic force, gravity, nuclear, annihilation (with the in-model vacuum-cure check as its pre-step: linearizing WITH the scalar/Gauss sector, plus the variable-λ mass profile from the working notes), lepton / neutrino family (carrying the localized-branch scan over `(λ, m_c)`), dark matter; the coverage column stayed staged as a draft until the research matured (published to [`MODELS.md`](../../../../../MODELS.md) at M7.21).
