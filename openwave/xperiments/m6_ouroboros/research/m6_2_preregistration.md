# M6.2 pre-registration: the H/Q decision gate, locked before the number

> **LOCK STATEMENT.** Task [M6.2](tasks/m6_2_task_details.md), 2026-07-20, locked at 11:05pm BEFORE `scripts/m6_2_hq_decision.py` was first run (file mtimes witness the order). Everything below (functionals, conventions, protocol, decision rule) is fixed here; nothing is adjusted after the numbers appear. Any post-lock change would go to the task's deviations log with its reason; none is intended. Derivation evidence: [`scripts/m6_2_derive_functionals.py`](scripts/m6_2_derive_functionals.py) → [`data/m6_2_derivation.json`](data/m6_2_derivation.json) (symbolic only, produces no H/Q number). Conventions inherited from the certified sheet: [`m6_1_v11_convention_sheet.md`](m6_1_v11_convention_sheet.md).

## 1. The pinned configuration (production form of record)

From [`archive/sandbox_v8/ouroboros_benchmark.py`](archive/sandbox_v8/ouroboros_benchmark.py), used verbatim:

| Item | Value |
| --- | --- |
| Reduction | azimuthal fields `A = φ̂ α(ρ)·T_A(t)`, `J = φ̂ β(ρ)·T_J(t)`, temporal gauge, z-independent, 2D-cylindrical measure `ρ dρ` (2π and per-unit-z cancel in all ratios) |
| ODE solved (N1) | `α″ + α′/ρ − α/ρ² + ω²α = β` · `β″ + β′/ρ − β/ρ² + ω²β = α − λβ − 4gβ³` |
| Parameters | g = 1.0, λ = 1.0, ω = 1.0, A₀ = B₀ = 0.1 |
| Solver | RK45, rtol 1e-9, atol 1e-11, max_step 0.02, r ∈ [0.02, 12], 800 pts, slope BCs (α ≈ A₀ρ, β ≈ B₀ρ near origin), localization tail < 0.15 |
| Phases | matched `e^{iωt}` (real parts); the printed cos/sin prescription is recorded as a labeled secondary (it kills the ⟨J·A⟩ term and makes (2.1) unsatisfiable; M6.1 F2 already found the printed ansatz section unreliable) |

## 2. What the derivation established (script-verified, before any number)

| # | Fact | Consequence |
| --- | --- | --- |
| DF1 | The real two-field theory has NO exact internal U(1): `J·A` is not invariant under the joint (A,J) rotation; the exact symmetries are spacetime ones (H, L, P) | **Q is a convention, not a Noether charge.** The only U(1)-derivable charge (complexified fields, joint phase) is `Q_joint = ω∫(α² + β²)ρdρ`. A J-only charge exists only by dropping the coupling term |
| DF2 | The period-averaged energy of the certified ℒ_ref on the pinned reduction is `⟨H⟩ = ¼∫[ω²α² + α′² + α²/ρ² + ω²β² + β′² + β²/ρ²]ρdρ − ½∫αβ ρdρ + (3/8)g∫β⁴ρdρ` (matched phases; gradient cross terms reduce by the boundary identity) | the production-coded `H_p = ∫[α′² + α²/ρ² + β′² + β²/ρ² + 4gβ⁴]ρdρ` is NOT the energy of ℒ_ref under any convention examined: it omits the ω² electric terms and the cross term, and its quartic coefficient (4g vs ⅜g) matches no normalization |
| DF3 | The `⟨H⟩ − λQ` constrained variation yields the OPPOSITE ω² sign from the production ODE (Q-ball ω-elimination subtlety); the production pair is therefore not that variational system | the production ODE's variational credential fails under the certified conventions |
| DF4 | The time reduction of the pinned dynamics (2.1)/(2.2) reproduces the production pair (sans λ-term) ONLY under the mostly-MINUS reading of (J·J), the opposite of the signature M6.1 certified from the printed interaction signs; the `−λβ` term arises from NO reading of the dynamics | the production ODE embodies a signature convention inconsistent with the printed-equation reconstruction, plus an injected multiplier term |

## 3. The observables evaluated on the N1 profile (all fixed here)

| Label | Functional | Status |
| --- | --- | --- |
| **PRIMARY: H_d/Q_d** | `H_d = ⟨H⟩` of DF2; `Q_d = ω∫(α² + β²)ρdρ` (DF1, the only derivable pairing) | the gate's number |
| Branch: H_d/Q_J | same H_d; `Q_J = ω∫β²ρdρ` | labeled convention branch (J-only charge is NOT derivable; reported for completeness) |
| Reference: H_p/Q_p | the production-coded pair (`H_p` above; `Q_p = ∫β²ρdρ`) | configuration certificate: expected ≈ 1.689 per v11 § 8; if this fails to reproduce, HALT and diagnose the configuration before any verdict |
| Secondary phases | H_d recomputed with the printed cos/sin phases (⟨J·A⟩ = 0 variant) | reported, labeled |

## 4. The decision rule (fixed before the number)

Let `gap = |H_d/Q_d − 1.6875| / 1.6875` (PRIMARY vs the physical target; the model-internal 1.6969 is reported alongside as comparison-only, canonical § 4).

| Outcome | Verdict |
| --- | --- |
| gap ≤ 1% | **branch (a)**: the electron benchmark survives a no-search derivation; the cell is re-earned |
| 1% < gap ≤ 5% | partial: reported as-is, user decides the branch at REVIEW |
| gap > 5% | **branch (b)**: the benchmark does not survive; the electron sector closes honestly (MODELS.md + briefing + hunt synced at REVIEW), M6 holds on the DM sector + M7 lineage |

The no-search rule: no functional, convention, parameter, or solver setting may be changed after this lock to move any number toward any target. Convention forks are resolved by derivation argument above, never by outcome. All four numbers in § 3 are reported regardless of the verdict.

## 5. Not computed (scope pins)

| Not computed | Why |
| --- | --- |
| Localized solutions of the DERIVED EL system (DF3) | different task; tonight's gate evaluates derived observables on the production profile of record |
| L, g-factor | v11's own footnote concedes L = ωQ is an identity; out of gate scope |
| Any 3D/toroidal generalization | M7's territory |

Cross-links: task [`tasks/m6_2_task_details.md`](tasks/m6_2_task_details.md) · derivation [`scripts/m6_2_derive_functionals.py`](scripts/m6_2_derive_functionals.py) + [`data/m6_2_derivation.json`](data/m6_2_derivation.json) · numeric gate [`scripts/m6_2_hq_decision.py`](scripts/m6_2_hq_decision.py) (runs after this lock) · sheet [`m6_1_v11_convention_sheet.md`](m6_1_v11_convention_sheet.md) · roadmap [`m6_roadmap.md`](m6_roadmap.md).
