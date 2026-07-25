# M7.7, consolidate the M7 column (MILESTONE)

> Task **M7.7** (M7 / HydroBoros). taskID = M7.N iteration. Status: **Done** (2026-07-04, review approved with the preview-staging decision) · Roadmap: [`../m7_roadmap.md`](../m7_roadmap.md)

This doc is the task's full record: planning + findings + documentation. **M7.7 is the Phase 1 milestone**: fold the winning recipe (M7.1-M7.6) into a canonical spec + a small auditable physics module (the [`METHOD_NOTE.md`](../../../../../dev_docs/METHOD_NOTE.md) standard), land the **HydroBoros (M7) column in MODELS.md** with honest icons, decide the units contract, and prepare the milestone stop's two comms packages (content bullets; the maintainer phrases everything).

---

## Plan

| Step | What | Gate |
| --- | --- | --- |
| 1 | **`m7_functional.py`**: the small physics module (~200 lines): the doublet convention, `E_ω` density term by term (numpy reference), both constraints, the core observables; docstring = the equations | the Taichi engine agrees with the module's reference energy on the canonical state (cross-validation, machine-ish) |
| 2 | **`m7_7_canonical.py`**: ONE runnable canonical script: rotating blend seed → fixed-`Q_can`+`H_A` relax → gate table (E, `\|g\|`, `H_A`, `Q_can`, j_z, Q_ρ) vs recorded values with tolerances | **reproducible first-try** (fresh run passes all gates) |
| 3 | **`m7_theory_canonical.md`**: the canonical spec (equations FIRST, pinned conventions with their provenance, the recipe, the gates, the units-contract decision table, honest limits incl. Q14) | METHOD_NOTE-shaped; every claim links its task doc |
| 4 | **MODELS.md column**: all 21 cells with honest icons + backing scripts; summary count; per-model results-of-record row; the briefing updated from pre-implementation | cells match the task-doc findings; no icon inflation |
| 5 | **Comms packages**: Marc (Q7, Q3, Q10, Q1, Q4) + Werbos (Q14 first, Q11, Q13, Q12+Q6 residual, Q9) as content bullets, METHOD_NOTE-compliant skeletons (permalinks pinned after the maintainer commits) | packages ready for the maintainer's voice; the new-model governance issue text drafted (creation = the maintainer's call) |

Honest-icon pre-commitment (drafted at PLAN, finalized against the docs): ⚠️ for charge / mass / clock / stability / μ+spin / Coulomb / EM waves / KG (each with its named caveat: windowed ledger Q11, real-time tachyon Q14, fixed-reservoir monopole, charge-unit-blocked μ, ℏ/2-vs-ℏ contract); 🚧 for the 13 Phase 2 + 4 cells. **No ✅ claimed at this stage**: the column enters honestly at 8 partials, and upgrades land cell by cell as the caveats close.

Artifacts: `scripts/m7_functional.py` + `scripts/m7_7_canonical.py` + `m7_theory_canonical.md` (research root) + the MODELS.md edits + `data/m7_7_canonical.json`.

---

## FINDINGS (2026-07-03/04, execution)

Artifacts: [`../m7_theory_canonical.md`](../m7_theory_canonical.md) (the spec) · [`../scripts/m7_functional.py`](../scripts/m7_functional.py) (the physics module) · [`../scripts/m7_7_canonical.py`](../scripts/m7_7_canonical.py) (the one-script reproduction) · [`../data/m7_7_canonical_full.json`](../data/m7_7_canonical_full.json) · [`../data/m7_7_canonical_quick.json`](../data/m7_7_canonical_quick.json) · [`MODELS.md`](../../../../../MODELS.md) (the published column; the M7.7-era preview draft was removed at M7.21) · [`__M7_model_briefing.md`](../../__M7_model_briefing.md) (status updated)

### 1. The canonical reproduction ✅ first-try

| Gate (N = 64 record run) | Value | Verdict |
| --- | --- | --- |
| `E` | 6.32462 (= M7.6's 6.3246) | ✅ |
| `‖∇E‖` | 1.6e-7 | ✅ |
| `j_z` per quantum (A / J) | 0.9939 / 0.9934 | ✅ |
| `Q_ρ(<0.3L)` | 0.0259 | ✅ |
| **METHOD_NOTE cross-validation** (Taichi engine vs `m7_functional.py` reference) | **1.43e-14** | ✅ machine agreement |
| constraints held (`H_A`, `Q_can`) | 4.1e-8 / 0.0 | ✅ |

Runtime 566 s; the physics now lives in ONE ~200-line module whose docstring carries the same equations as the spec, and the driver contains no physics of its own (the M5.16 lesson, structurally enforced). The **quick mode** (N = 48, the anyone-reproduces path) also passes all gates first-try: E = 6.2866 (grid-shifted, in bracket), `j_z` = 0.9885/0.9888, cross-validation 4.9e-15, 201 s.

### 2. The 21-cell column: STAGED as a preview (the milestone, revised at review)

The **HydroBoros (M7)** column was drafted in full (**0 ✅ / 8 ⚠️ / 0 ❌ / 13 🚧**; no icon inflation: every ⚠️ carries its named caveat: windowed ledger Q11, real-time tachyon Q14, fixed-reservoir monopole, charge-unit-blocked μ, the open units contract) and, **per the review decision (2026-07-04), staged as a draft instead of entering MODELS.md now**: the research program is still in flight and not yet ready for the cross-model benchmark. MODELS.md stays untouched; the column enters via the governance flow (issue + script-backed PR) at the M7.21 readiness call, where it would currently rank fourth (8, tied with M4). The briefing graduated from pre-implementation to Phase-A-complete-with-staged-preview. The units-contract decision table (ℏ/2 vs ℏ ↔ ω_C vs ω_D) is in [`../m7_theory_canonical.md § 4`](../m7_theory_canonical.md) with the Zitter-mapping recommendation on record, decision = user's call.

### 3. The comms packages (content bullets; the maintainer phrases everything)

**Package W (Werbos): Q14 → Q11 → Q13 → Q12 (+Q6 residual) → Q9.** Content bullets:

- the vacuum tachyon: `det M(0) = −1` for any `f` (one-line determinant fact); band `k² < 0.618` (repulsive) / `1.618` (focusing); measured growth rate 0.785 vs analytic 0.786, amplitude-independent → **no `β*` threshold exists in the vector truncation** (the Werbos v5 Test 1); ask = what cures the vacuum in the full model (A-mass / condensate vacuum / scalar-Gauss sector / parameter islands)
- the constructive face: `E_ω` PSD iff `ω ≥ ω* = k* = 0.786`; soliton existence threshold MEASURED (runaway at 0.70/0.75, solitons from 0.79, `Q_J` diverging toward ω\*): **the de Broglie clock is the vacuum-stability mechanism** (particles exist because they tick fast enough): a strong Ouroboros-thesis statement if the full theory confirms it
- Q11: `H/Q = 1.6890` is window-defined (both far-field roots propagating at the canonical point); the same dispersion's upper branch is exact KG (`m_eff² = (1+√5)/2`); consequence measured: neutral solitons interact via an oscillatory RKKY-style exchange, period `π/k`
- Q13: the M6 electron is a 3D constrained saddle (focusing collapse); adding helicity to the same torus stabilizes it; `ω*` sharpens the 1D-vs-3D conjugate-point ask
- Q12 + Q6 residual: the benchmark ODE = the same-phase doublet EL with focusing pins; the WRITTEN `f = gs²` is the branch with stable 3D solitons; `0d_canonical § 2.2`'s `2ωα` form matches no EL reduction; confirm intent + the `(g/4)` doc fix
- Q9: (Ω, G) definitions still needed; candidate anchor offered: Werbos's stability islands may be the `ω > ω*(g, λ)` region; one concrete electron-island parameter point would let us check on the lattice
- deliverables behind it: the M7 MODELS.md column + [`m7_theory_canonical.md`](../m7_theory_canonical.md) (equations first + equation-to-code map + the reproduction script)

**Package M (Marc): Q7 → Q3 → Q10 → Q1 → Q4.** Content bullets:

- the headline gift: the author's torus, made 3D-dynamical: the rotating blend relaxes to a stable soliton that is a clean **`j_z = 1` per-quantum wave (0.6%)**: the circulating-photon `L = ℏ` reading, measured; ball-vs-torus: the ball is the ground state at canonical parameters, the shape question is now quantitative
- Q7: the author's "start Trkalian, take off the training wheels" implemented as relax-and-measure and VALIDATED (approximately-Beltrami cores, `\|align\| = 0.96`); the charge prescription self-resolved at the fixed-reservoir level (Gauss 99.1%, far field −2.14, two-charge `1/d` reference-matched, **1.17 ± 0.02 dressing measured**); asks = the S&Y variable-h recipe; does the author's picture single out a `j₀` profile
- Q3: divergence charge and linking charge measured INDEPENDENT at this level (linking gates existence; RMS charge coexists): the author's read vs Pisello/Gauss-Bonnet?
- Q10: the Eq 122/124/127 algebra (M7.2): corrected `U ≈ 0.958 m_ec²`: confirm intended
- Q1: substrate reading (doublet vs single RS field) + target manifold (S² vs S³)
- Q4: the promised Beltrami material + the Enciso / Peralta-Salas collaboration status
- deliverables behind it: same column + spec; the M7.2 digit-for-digit reproduction of the author's printed solution

Both packages follow [`METHOD_NOTE.md`](../../../../../dev_docs/METHOD_NOTE.md) (equations first; the equation-to-code map gets commit-pinned permalinks once the milestone commit lands; ≤4 artifacts each: the spec, the module, the canonical script, the relevant plot).

### 4. Governance (deferred with the column; issue text drafted for the M7.21 entry)

Draft title: `New model column: HydroBoros (M7), toroidal-Beltrami electron on the Ouroboros doublet`. Draft body bullets: what M7 is (two named parents, credited per cell); Phase A complete with the one-script reproduction (`m7_7_canonical.py`, all gates first-try, engine-vs-reference 1.4e-14); the column enters at its then-current icons with named caveats (honest-icons rule); deep links (canonical spec, preview, roadmap, tracker). Per [`MODELS.md § Contributing`](../../../../../MODELS.md): issue first, then the script-backed PR with DCO. **Per the review decision the whole entry is deferred to M7.21**; until then the column lived as a draft and MODELS.md was untouched (published at M7.21, 2026-07-18).

---

Cross-refs: roadmap [`../m7_roadmap.md`](../m7_roadmap.md) (M7.7, Phase D) · upstream [`m7_1_infra.md`](m7_1_infra.md) · [`m7_2_fleury_torus.md`](m7_2_fleury_torus.md) · [`m7_3_ouroboros_3d.md`](m7_3_ouroboros_3d.md) · [`m7_4_charged_soliton.md`](m7_4_charged_soliton.md) · [`m7_5_clock_stability.md`](m7_5_clock_stability.md) · [`m7_6_observables.md`](m7_6_observables.md) · tracker [`../m7_question_tracker.md`](../m7_question_tracker.md) (comms plan + Q14) · [`MODELS.md`](../../../../../MODELS.md) · standard [`dev_docs/METHOD_NOTE.md`](../../../../../dev_docs/METHOD_NOTE.md).

---

## TASK REVIEW (2026-07-04)

**Task Duration:** 00:19 (from 21:52 to 22:11 EDT; closeout adjustments 2026-07-04)
**Usage Cap Triggered:** NO (ping armed for 2:15am, parked without firing)

**Results** (full detail: [`§ FINDINGS`](#findings-2026-070304-execution)):

| Gate | Outcome |
| --- | --- |
| one runnable canonical script, first-try | ✅ ALL GATES PASS at BOTH N = 64 (E = 6.32462, 566 s) and N = 48 (201 s) |
| METHOD_NOTE enforced structurally | ✅ physics in ONE ~200-line module; driver physics-free; **engine-vs-reference cross-validation 1.4e-14 / 4.9e-15** as a first-class gate |
| canonical spec | ✅ [`m7_theory_canonical.md`](../m7_theory_canonical.md): equations first, provenance-pinned conventions, equation-to-code map, units-contract decision table, not-computed list |
| the 21-cell column | ✅ drafted at 0 ✅ / 8 ⚠️ / 13 🚧 (named caveats) and **staged as a draft** per the review decision: MODELS.md untouched at M7.7, benchmark entry deferred to M7.21 (published 2026-07-18) |
| comms packages | ✅ content bullets ready (W: Q14-led · M: Q7-led, § 3); governance issue text drafted for the M7.21 entry |

**Issues / blockers:** none; permalink pinning follows the maintainer's commit; no >1MB data created.

**Question audit:** consolidation task, no resolutions; count 5/5/0/4.

**Action needed (the milestone stop, the maintainer's items):** audit the diff, decide the units-contract row (spec § 4), commit + push (message proposed in the terminal review), then say the word for permalink pinning; the packages go out in the maintainer's voice from § 3.

**Findings**: Phase 1 closed with the consolidation the METHOD_NOTE standard demands (a spec auditable by reading, a 200-line term-by-term module, a one-script reproduction that passed every gate first-try at two resolutions with machine-level engine-vs-reference agreement), and with the benchmark discipline the program demanded of itself: the 21-cell column was drafted honestly and then deliberately held back from MODELS.md as a staged preview, because a research program still carrying its own open theory question (Q14) does not benchmark against mature columns until it is ready.

**Research docs created / updated**:

- [this task doc](m7_7_consolidation.md) (plan + FINDINGS §§ 1-4 + both comms packages)
- [`../m7_theory_canonical.md`](../m7_theory_canonical.md) (the canonical spec) · [`MODELS.md`](../../../../../MODELS.md) (the published column; M7.7 preview draft removed at M7.21)
- [`../scripts/m7_functional.py`](../scripts/m7_functional.py) · [`../scripts/m7_7_canonical.py`](../scripts/m7_7_canonical.py)
- [`../data/m7_7_canonical_full.json`](../data/m7_7_canonical_full.json) · [`../data/m7_7_canonical_quick.json`](../data/m7_7_canonical_quick.json)
- [`../../__M7_model_briefing.md`](../../__M7_model_briefing.md) (Phase-A-complete, preview-staged)
- [`../m7_question_tracker.md`](../m7_question_tracker.md) (M7.7 chronology; count 5/5/0/4); MODELS.md restored untouched
