# M5.20.4: the formulation search: a usable dynamics for the quartic L, found by us

**Status**: ✅ CLOSED COMPLETE 2026-07-14 (same-day: go 14:55, review approved 16:21+; roadmap row in [DONE](../m5_roadmap.md)). Deliberately NOT gated on a Duda answer: the strategy (user, 2026-07-14) was to search the formulation space ourselves first and gate on him only where a question is genuinely author-gated: he himself said the regularization is the hardest part and that he has tried many: we searched where he has not. Outcome: arm A (his least-action branch) WINS with a concrete periodic-orbit family; full record in [FINDINGS](#findings-2026-07-14-method-note--findingsm5_20_4_method_notemd) + the [method note](../findings/m5_20_4_method_note.md).

**Lineage**: [M5.20.3](m5_20_3_task_details.md) (the verdict this task answers: free-EL IVP ill-posed, every regularization blows up, E → −∞, his BVP branch fires; the true-L instrument `m5_20_3_a_constraint.py` with GC0a-e green) · tracker [Q24](../m5_question_tracker.md#q24-detail) (the formulation question: rides the next outbound in parallel; his answer, if it lands, folds in mid-task) · the M5.12 phase-D time-periodic BVP instrument (the arm-A machinery) · [`../m5_particle_hunt.md`](../m5_particle_hunt.md) (both hunts funnel into this task's outcome).

## TASK PLANNING (2026-07-14, at the M5.20.3 close)

### Scope

Search three self-directed formulation routes for a USABLE dynamics of the purely-quartic verified L on the vortex-loop background, each with pre-registered kill criteria; deliver a per-arm verdict and, if an arm survives, demonstrate it on the loop (the M5.21.1 hedgehog run then inherits the winner). One contingency arm held in reserve.

| Arm | Route | The question it answers | Kill criterion (pre-registered) |
| --- | --- | --- | --- |
| **A** (primary) | **Time-periodic least-action BVP under the true L**: extremize the action over periodic orbits (the period IS the clock); re-derive the M5.12 phase-D instrument (`Shat = S0 − ω²Q2`, balance root, the `H_mean = 0` identity) with Q2 built from the TRUE kinetic form K(M) instead of the canonical norm; start at the rigid-orbit level (the M5.20.2 H(ω) machinery re-read through the action), then Fourier-in-time profile relaxation if the rigid level signals | Does a finite-ω periodic orbit (a particle clock) EXIST in the formulation he himself endorses (least action, two-sided BCs)? If yes: containment + the radius-breathing observable become measurable at last | No balance root on any orbit family at budget, or every extremum is a runaway direction (S unbounded along the family): then the BVP search is honestly negative at the seed-family level and the result is the sharpest possible Q24 follow-up |
| **B** | **Dirac consistent initial data**: the from-rest EL-inconsistency (98.6% null-force) says (M, V = 0) at loop seeds is off the theory's constraint surface. Find data ON it: (i) config route: solve P_null(M)·G_static(M) = 0 near the loop (q intact); (ii) velocity route: at fixed M solve RHS_null(M, V) = 0 for V (quadratic in V, Newton). Then evolve and measure surface preservation (nff(t)) and whether the blowup persists on-surface | Is the IVP well-posed on the theory's OWN consistent sector (the mathematically required Dirac-Bergmann treatment, not an added restriction)? Almost certainly unsearched by the author | The surface is empty near loop states (decisive, reportable), or evolution leaves it at a rate that makes it useless, or the blowup persists on-surface (which strengthens the BVP-only conclusion) |
| **C** (diagnostic) | **The sanctioned-term kinetic fix**: enumerate the Lorentz-invariant gradient-term classes he has sanctioned (the Q13 redirect: Skyrme-like terms on the 4×4 tensor; note the η-metric QUADRATIC term tr(ηṀηṀ) is itself indefinite in the time-mixing sector: the derivation decides what is admissible), derive each candidate's contribution to K(M), and MEASURE whether any closes the indefiniteness (K positive-semidefinite on the measured backgrounds) while preserving the statics anchors (r_half, the Coulomb lock, the core gate: gate re-runs, not production) | Does a term HE already sanctioned fix the kinetic form at the source? Nobody has computed the K-effect of the Q13 class | No sanctioned term closes K on the measured backgrounds, or every term that does breaks a statics anchor gate |
| **D** (contingency, only if A-C all kill) | **The ε-canonical bridge**: L + ε·canonical kinetic; t*(ε) and stability-window scan; if a stable window exists, measure the oscillation observables and study the ε → 0 limit | A pragmatic numerical bridge to the observables, documented as regularization | No stable window, or observables do not converge as ε shrinks |

Run order: **C-census first** (cheapest; if a term fixes K the other arms change meaning), then **B**, then **A** (the deepest and highest-payoff). D only on triple kill. All arms on the loop at 64×128 (the M5.20.3 grids), δ = 0.3, the recipe seed as primary.

NOT in scope: production oscillation runs under a winning formulation (that is the successor, with M5.21.1 taking the hedgehog); any new outbound content beyond folding (the Q24 email is already drafted and user-gated); the AMBer fit.

### Definition of done

| ✅ when | Bar |
| --- | --- |
| Per-arm verdict delivered | Each arm reaches its pre-registered survive/kill criterion, machine-checked (try cap 3 per gate); timeboxes: C ≈ 2 h, B ≈ 2-3 h, A ≈ 4-5 h; the tail ships partials |
| Instruments gated before physics | Arm A: the quartic-Q2 derivation gated complex-step against the discrete action (the GC0 pattern); Arm B: the surface residual + its Jacobian gated; Arm C: each candidate term's K-contribution gated against a from-scratch build (the GC0d pattern) |
| A winner demonstrated (if any survives) | Arm A: a balance-root orbit with its ω, energy, and containment read on the loop; Arm B: an on-surface evolution to T ≥ 50 with nff ≈ 0 and the ledger clean; Arm C: a K-closing term with ALL statics anchor gates re-run green |
| Honest disposition if all kill | The triple negative + arm-D scan = the sharpest possible Q24 follow-up (we searched his space AND ours); the outbound updates only through the user |
| Records | Method note `../findings/m5_20_4_method_note.md` (equations first, code map, figures per the film standard where states evolve); independent adversarial audit (cardinal rule); tracker Q24 + convo + roadmap routing; film strips per [`../m5_visualization.md`](../m5_visualization.md) (basic template; blowup spacing where applicable) |

### Gating

Roadmap `Gated By`: user "go" only. Q24's outbound rides in parallel; **folding rule**: if Duda's Q24 answer lands mid-task, fold it as a branch decision at the next arm boundary (the M5.20.3 folding-table pattern), log the fold in this file, and continue: his answer refines, never blocks.

### Blindspot pass (unfamiliar territory: action-level search + constraint surfaces + term enumeration)

| # | Unknown unknown surfaced | Fold into plan |
| --- | --- | --- |
| 1 | The quartic Q2 may itself be degenerate/indefinite on orbit families (the K_eff census in action language): the balance-root machinery assumed a sign | Pre-register the reading: an indefinite Q2 is the H(ω) census re-derived, reported as structure, not failure; the orbit family census comes FIRST in arm A |
| 2 | Arm B's config-route solve may converge to the unwound remnant (the statics slide) instead of a loop state | Constrain the solve (winding read in the loop, the M5.19 corepin precedent as fallback); an empty ON-loop surface is itself the decisive answer |
| 3 | Arm C's admissible-term enumeration is a derivation risk (the η-quadratic is indefinite; Skyrme-like on a symmetric tensor has several contraction patterns) | Enumerate ALL contraction patterns at quadratic and quartic order first, machine-check each for Lorentz invariance with the M5.18 check-1 instrument, THEN derive K-contributions only for the invariant ones |
| 4 | Wall-clock: three deep arms in one task | Hard timeboxes per arm (DoD); the run order puts the cheapest decisives first; hard cap ~1 day, ship partials + the sharpest question |
| 5 | The winner may behave differently on the hedgehog (different K census: rank 8 core vs rank 5) | Out of scope here: M5.21.1 inherits the winner AND re-runs its gates on the hedgehog before physics (the series re-capture rule) |

### Research-body destinations

| Artifact | Destination |
| --- | --- |
| Scripts | `../scripts/m5_20_4_a_bvp.py` · `m5_20_4_b_dirac.py` · `m5_20_4_c_terms.py` · (`m5_20_4_d_bridge.py` contingency) · `m5_20_4_plots.py` |
| Data / plots | `../data/m5_20_4_*.json` (+ npz ≤ 1 MB; larger deleted at FINISH with regen docs) · `../plots/m5_20_4_*.png` (film standard) |
| Findings | `../findings/m5_20_4_method_note.md` |
| Records | This file (FINDINGS + TASK REVIEW) · tracker Q24 · [`m5_20_convo.md`](m5_20_convo.md) if outbound content changes (user-gated) · `../checkpoints/m5_20_4_progress.md` (opens at go, resume-complete) |

### Preconditions

| Precondition | State |
| --- | --- |
| M5.20.3 instrument + verdicts (K builds, accel, the blowup card) | ✅ closed 2026-07-14 |
| M5.12 phase-D BVP scripts on disk | ✅ frozen at the M5.12 close (regen documented there) |
| The Q13 sanctioned-term record (his 2026-07-06 slides + the redirect decode) | ✅ tracker Q13 + `m5_18_convo` record |
| Q24 outbound | drafted, user-gated, independent of this task |
| Resume ping + checkpoint | 🚧 at go (user supplies resets_at) |

### Model + effort

Fable 5 high for the arm-A action re-derivation and the arm-C term enumeration (novel derivations); default effort for the mechanical gate loops; independent agent with own instruments for the audit (cardinal rule).

### Deviations log (EXECUTE, as they happen)

| # | Deviation | Why / disposition |
| --- | --- | --- |
| 1 | Arms ran partly in PARALLEL (C foreground first, then B and A interleaved while C/B long jobs ran in background) instead of strictly C → B → A | wall-clock efficiency; the arms use independent instruments and share only read-only backgrounds: no cross-contamination |
| 2 | Arm A gained an unplanned discovery channel: the ELLIPTIC (boost-conjugated rotation) orbit families (a1b) | the a1 census showed rotations Q2 > 0 / boosts Q2 < 0 with boosts non-periodic; conjugation Λ J Λ⁻¹ preserves periodicity while importing boost content: this turned out to be THE root channel |
| 3 | Arm A gate AG3 was redefined mid-run | first cut wrongly expected Q2 = 0 on the co-rotating vacuum; the equivariant vacuum carries the frame texture (A_φ ≠ 0, K ≠ 0): the gate now uses a genuinely gradient-free (rotation-invariant transverse) uniform state, and the co-rotating vacuum's Q2 = 1706 is recorded as physics |
| 4 | Arm C added the texture-attribution measurement (U_C on the PURE vacuum texture) | needed to attribute the non-perturbative U_C(β*) ≈ 2e5 to the zero-energy Goldstone texture rather than the loop: it became the obstruction lemma's measured half |
| 5 | Arm B's b2 needed multi-start (r(c) = r₀ + Q(c,c): zero Jacobian at the origin) | discovered at implementation: both grad_M T and kdot are exactly quadratic in V; documented in the method note § 2.1 |
| 6 | The audit round produced two same-day corrective measurements not in the plan: a1c (the φ-averaged Q2 elliptic re-scan) and b2c (the structured-velocity cancellation scan) | the audit refuted C6c and qualified C7; rather than shipping refuted claims, both were re-measured under the corrected reading before FINISH (the roots survive; arm B re-scoped to zero-energy-closed / full-V-open) |

### Contingencies + comms

| Trigger | Action |
| --- | --- |
| Duda's Q24 answer lands mid-task | Fold at the next arm boundary (folding rule above), log, continue |
| Hard cap (~1 day) | Ship per-arm partials + ONE sharp follow-up question |
| An arm produces an outbound-worthy positive early (e.g. a balance-root clock) | Checkpoint it, finish the arms, present at REVIEW; outbound stays user-gated |

## FINDINGS (2026-07-14; method note = [`../findings/m5_20_4_method_note.md`](../findings/m5_20_4_method_note.md))

The one-line verdict (audit-corrected): **arm C kills as scoped (with a structural obstruction lemma AND an audit-discovered escape candidate), arm B's zero-energy route is closed by measurement (the full-velocity route re-opened by the audit), and arm A SURVIVES at the rigid level on the audit-corrected functional: free-period least-action roots exist on the loop, on exactly periodic orbits, with the H = 0 identity.** The formulation winner for M5.21.1 inheritance is the least-action BVP route, now with a concrete orbit family in hand plus the audit's γ = −1 Skyrme-subtraction completion as the hypothesis-status runner-up.

| # | Finding | Status |
| --- | --- | --- |
| 1 | **The null is class-structural (lemma L1)**: NO commutator-class ("Skyrme-like") term, dressed or not, can lift the exact Ṁ ∝ η null (X₀ = I commutes with everything); machine-verified at 1.1e-13. The Q13-sanctioned class cannot fix the primary constraint even in principle | ✅ derived + machine-checked |
| 2 | **Exactly one invariant family closes the kinetic form**: the two-sided dressing −β η^μν tr(X_μPX_νP), P = aI − ηM, window a ∈ (1, g) (measured (1.25, 7.75)); it closes ALL cells at β* ≈ 1.306 on every background (quartic min-eig −32 vs dressing +24/unit-β); one-sided, trace-split, bare, and dressed-Skyrme candidates all fail (derived + measured; s2 min-eig < 0 at every sampled cell and a) | ✅ measured |
| 3 | **The closure-vs-statics obstruction**: any single closing quadratic term must charge the (1,2) off-diagonal sector (PSD sign algebra, incl. the P ≠ Q generalization), which is exactly where the equivariant vacuum texture lives: U_C(pure texture) = 9.1e4 at β* while its true energy is EXACTLY 0. At β* the statics minimizer destroys the loop (audit-verified: \|M13\| collapses to 7e-12, 1218-center scan finds no winding anywhere): **arm C KILL per pre-registration**, with a structural explanation for why single-term regularizations fail. ⚠️ The audit REFUTED the "mutual exclusivity" headline AS STATED: the combination **L − 1·s2 + β→0⁺·qc** (from this task's own term list) is PSD-marginal at machine zero on all backgrounds with vanishing texture cost: recorded as the hypothesis-status escape (γ = −1 sign admissibility author-gated; statics gates + ghosts + band-confinement unmeasured) | ✅ measured + derived; escape 🔶 |
| 4 | **The closed dynamics is well-posed but is a different theory**: full-rank mass matrix, T = 10 reached, E conserved to 3.7e-8 over 8000 steps on the dressed minimizer (vs every true-L run blowing at t* ≤ 7.2). The LOOP seed at β*: no blowup to T = 10 (texture slosh, ring persists ~17, winding read churns). Sub-closure β = 0.01: blowup at t = 6.87 (delayed 13×, not removed): no stable window below closure: arm D (the ε-bridge) is subsumed and answered by the invariant β-ladder | ✅ measured |
| 5 | **From-rest inconsistency is total and localized**: nff = 0.9999977 at strict cutoff; the null force lives in the b0 time-diagonal sector (0.43 vs ≤ 2.8e-4 elsewhere) | ✅ measured |
| 6 | **Arm B (audit-corrected): the ZERO-ENERGY route is closed, the full-velocity route re-opens**: r(c) = r₀ + Q(c,c) exactly (no linear term); 6 JFNK multi-starts stall at relative reduction 2e-5; b2c: null-projected structured combos cancel exactly NOTHING (\|r\|/\|r₀\| = 1.0000) while full-V structured velocities (the audit's (0,2) bump, align −0.479) reach 0.8855: consistent data with kinetic energy MAY exist (hard underdetermined solve, successor); the from-rest inconsistency (nff = 0.9999977) stands | ✅ measured |
| 7 | **Arm A SURVIVES (rigid level, on the audit-corrected functional): balance roots EXIST on periodic orbits**: the audit caught the equivariance defect (slice Q2 ≠ 3D orbit kinetic outside the J12-commutant); the corrected Q2_avg (φ-averaged generator; T and U constant along the true 3D orbit) STILL crosses zero on ALL SIX exactly-periodic elliptic families (J23^K2 ≡ J13^K1 at χ = 0.75, the rest at 0.25): with U = 0.344 > 0 the free-period root ω* = sqrt(−U/Q2_avg) is real, H = 0 exact | ✅ measured (corrected) |
| 8 | **The descent machinery works and the loop survives it**: at the deepest slice-functional root the residual drops ×75 in 120 steps with q staying exactly 0.5 and U at loop scale (0.344 → 0.315). ⚠️ Audit: the descended functional is the SLICE Q2 (variationally unfounded as the 3D orbit action): the corrected-functional extremal solve (Q2_avg, 16-24× cost) is the successor's first step | ✅ measured, re-scoped |
| 9 | **The U < 0 route does not exist in the probed family**: static clock-sector dressing never drives U below 0 (V4 dominates all 15 amp × width points, q intact): the roots live ONLY on the elliptic (negative-Q2) route: the negative kinetic directions ARE the clock's enabler in the BVP formulation, quantifying the author's "negative Hamiltonian terms propel angular momentum" | ✅ measured |
| 10 | **Adversarial audit: 5 CONFIRMED / 1 QUALIFIED / 2 REFUTED** (independent agent, own instruments: `../scripts/m5_20_4_audit_check.py` + `../data/m5_20_4_audit.json`): C1-C5 confirmed at machine precision (incl. the β* factor-pairing to 1e-12 and the statics-gate loop-destruction read); C6c REFUTED (structured velocities reach the b0 sector: → the b2c re-scope, row 6); C7 QUALIFIED (the slice-equivariance defect: → the q2_avg correction + a1c re-scan, row 7); C8 REFUTED as stated (the γ = −1 escape: → row 3). Both refutations were converted into same-day corrective measurements; the full verdict table is the method note § 5 | ✅ audited |

**Route forward (feeds the roadmap)**: M5.21.1's gate stays on this task's winner = the least-action BVP formulation; the named successor steps are (i) the saddle-aware extremal solve under the CORRECTED Q2_avg from the a4 state, (ii) stationarity over χ coupled in, (iii) the Fourier (non-rigid) profile container rebuilt on the CURRENT stack (the M5.12 one is era-drifted), THEN the loop's oscillation observables (radius breathing, gap-ladder spectrum) read on the converged orbit, and only then the hedgehog inheritance. Parallel candidates from the audit round: the γ = −1 Skyrme-subtraction completion (gate program: sign admissibility ask, statics anchors, ghosts, band confinement) and the full-V Dirac solve with structured initialization.

## TASK REVIEW (2026-07-14)

**Task Duration:** 01:26 (from 14:55 to 16:21 EDT)
**Usage Cap Triggered:** NO

**Results** (full tables in [FINDINGS](#findings-2026-07-14-method-note--findingsm5_20_4_method_notemd) above and the [method note](../findings/m5_20_4_method_note.md)):

| Result | Status |
| --- | --- |
| Arm A SURVIVES: free-period balance roots (ω*² = −U/Q2, H = 0 exact) exist on exactly-periodic boost-conjugated rotation orbits; confirmed under the audit-corrected φ-averaged kinetic (all six families cross); descent holds q = 0.5 with residual ×75 down | ✅ measured |
| Arm C KILL as scoped: the η-null lemma + the closure-vs-statics obstruction (β* = 1.306; the texture charge is forced; loop destroyed at closing β, audit-verified) | ✅ measured + derived |
| The audit's escape candidate: L − 1·s2 + β→0⁺·qc (PSD-marginal machine zero, texture cost → 0; γ = −1 admissibility author-gated, gates unmeasured) | 🔶 hypothesis |
| Arm B: zero-energy Dirac route closed (null-projected cancellation exactly 0); full-V route re-opened by the audit (successor) | ✅ measured |
| Closed dynamics well-posed (E drift 3.7e-8) but a different theory; sub-closure still blows (t* = 6.87): arm D subsumed | ✅ measured |
| Adversarial audit 5 CONFIRMED / 1 QUALIFIED / 2 REFUTED; both refutations re-measured same-day (a1c, b2c) | ✅ audited |

**Issues / blockers**: the β* slosh run's 2% E-drift is dt-marginal (diagnostic only); a4 = machinery demonstration on the slice functional, not an orbit solve; the escape and the full-V Dirac solve are successor-scale.

**Deviations from plan**: six, logged in the [Deviations log](#deviations-log-execute-as-they-happen).

**Action needed**: successor task = the corrected-functional extremal-orbit solve (**[M5.20.5](m5_20_5_task_details.md)**, PLAN written same day, top of Backlog); the Q24 outbound draft now carries three sharpened asks (BVP confirmation, elliptic-orbit class, γ = −1 admissibility), user-gated (the M5.20.3-close pattern).

**Findings**: the formulation search eliminated the alternatives with structure instead of fizzles: single sanctioned terms cannot fix the kinetic form without destroying the loop (lemma + measurement), zero-energy Dirac data does not exist at the loop, and Duda's own least-action branch positively signals: free-period particle-clock roots exist on exactly periodic boost-conjugated orbits, the theory's negative kinetic directions doing the enabling, all surviving an adversarial audit that also contributed a candidate well-posed completion (γ = −1 Skyrme subtraction). The oscillation program is unblocked without waiting on Duda.

**Research docs created / updated**: [`../findings/m5_20_4_method_note.md`](../findings/m5_20_4_method_note.md) · this file · [`../m5_question_tracker.md`](../m5_question_tracker.md) (Q24) · [`../m5_roadmap.md`](../m5_roadmap.md) · scripts `../scripts/m5_20_4_{a_bvp,b_dirac,c_terms,plots,audit_check}.py` · data `../data/m5_20_4_*.json` · plots `../plots/m5_20_4_{a1b_roots,c_closure,a4_descent,film_recipe_closed}.png`.
