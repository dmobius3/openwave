# M5.21.11: the route-(b) pre-registration framework (derivation + frozen protocol)

> **Status: PRE-FREEZE DRAFT (2026-08-06).** This document is the pre-registration artifact required by the amended M8.6 reopening condition 3, route (b) ([readiness note § 8](../../../m8_mit/research/findings/m8_6_readiness_note.md)). It contains the derived asymptotic forms (§ 1), the frozen ladder protocol (§§ 2-6), and the compliance audits (§§ 7-9). **No rung has been measured, no fit run, no ratio computed**: the ladder compute is a separate follow-on run gated on this document being FROZEN. Freeze mechanics: after the author sanity-check of § 1 (or the user's call to proceed without it), this header changes to FROZEN with the pinning commit SHA; any edit after that point voids route (b) (failure is terminal, no second framework).

| Provenance | Task record: [`m5_21_11_task_details.md`](../tasks/m5_21_11_task_details.md) (planning + the resolved design question) |
| --- | --- |
| Consumes (measured inputs) | [M5.21.2b](m5_21_2b_note.md) (instrument of record, core-scale law, merge lesson) · [M5.21.1](m5_21_1_method_note.md) (P4 scaling laws) · [M5.21.8](m5_21_8_note.md) (m\* law, convention bridge) · [theory canonical](../m5_theory_canonical.md) §§ 1, 4, 5 |
| Governs | the follow-on ladder-compute run (not scheduled by this doc) |

## 1. The derivation (D1): asymptotic forms E_branch(δ, g) from M5-side theory

Everything in this section is derived from the M5 functional and from measured M5-internal scaling laws. No lepton mass, lepton ratio, or Yukawa figure enters anywhere (§ 7 audits this).

### 1.1 The functional and instrument of record

The census instrument ([M5.21.2b § 1/§ 4](m5_21_2b_note.md), the only term set whose bare minimum is simultaneously stencil-consistent and virial-balanced):

```text
E[M] = E_u[M] + E_V[M]                    M(x): real symmetric 3x3, spatial block

E_u  = h^3 SUM_x 4 SUM_{i<j} <F_ij, F_ij>          F_ij = [d_i M, d_j M]
E_V  = h^3 SUM_x w2 SUM_k (lam_k(M) - v_k)^2       v = (0, delta, 1) ascending  (T2)
```

evaluated on the stencil-symmetrized discretization E = ½(E_fwd + E_bwd), FIRE to f_tol, pinned boundary, with the consistency gates of record (cross-stencil ratio ≤ 1.5, virial residual, eigen-gap guard). The three lepton-candidate branches A, C, B are the three distinct protected minima of this functional. Equation-to-code map: both terms and every gate are implemented in [`m5_21_2b_a_instrument.py`](../scripts/m5_21_2b_a_instrument.py) (E_u in the stencil-symmetrized `e_parts`, T2 in `v_density`; the gate battery per its header), the SAME code the ladder will run; this document introduces no new energy code.

### 1.2 Derrick/virial structure

Under a spatial dilation M_s(x) = M(x/s): each factor ∂_i contributes s⁻¹ and the volume element s³, and E_u is quartic in first derivatives while E_V carries no derivatives:

```text
E_u(s) = E_u / s          E_V(s) = s^3 E_V
dE/ds |_{s=1} = 0   =>    E_u = 3 E_V        (the Faber virial balance)
=>  E_min = (4/3) E_u = 4 E_V
```

Measured at the T2 bare minimum: virial residual +0.034 at N = 48 ([M5.21.2b § 4](m5_21_2b_note.md)). Consequence used below: at every rung the minimum is scale-stationary, so the per-rung virial residual is a valid convergence-quality gate (§ 4), and the core size is dynamically selected by a curvature-vs-potential competition (§ 1.3), not imposed.

### 1.3 The δ-axis: far/core decomposition and the asymptotic form

**Far field.** Outside the core, M = O(x) D O(x)ᵀ with D = diag(1, δ, 0) and O ∈ SO(3). With Γ_i = Oᵀ∂_iO (antisymmetric), ∂_iM = O [Γ_i, D] Oᵀ and ([Γ_i, D])_ab = (Γ_i)_ab (d_b − d_a): every gradient entry is weighted by an eigenvalue gap, and the gaps are (1 − δ, 1, δ). F_ij is bilinear in two such factors and ⟨F_ij, F_ij⟩ is their square; since M is LINEAR in δ, the far-field energy density is EXACTLY a polynomial of degree ≤ 4 in δ (audit: a random smooth texture fits a quartic to relative residual 9e-16), with a finite δ → 0 limit (the δ-gap channel, the e₂-e₃ twist, decouples; the long-axis winding that carries the charge sits on the gap 1 − δ → 1). Hence

```text
E_far(delta) = E_far(0) * (1 + c1*delta + ...)     (exactly polynomial, degree <= 4)
```

**Core.** The measured cores are vortex LINES (the braided +½ pair for A; line bundles for the heavier branches, [M5.21.2b § 6](m5_21_2b_note.md)). Per unit line length, with core radius a and an O(1) eigenvalue rearrangement amplitude: the quartic curvature contributes ~ κ_u/a² (four gradient factors ~1/a over a cross-section a²) and the T2 penalty of the rearranged core contributes ~ κ_V w δˢ a², where δˢ is the penalty scale of the core's eigenvalue rearrangement (a clean 2-equal merge of the targets δ and 0 costs (λ − v)² ~ δ², i.e. s = 2; a partial rearrangement gives s < 2). Minimizing e(a) = κ_u/a² + κ_V w δˢ a²:

```text
a*(delta)      ∝ delta^(-s/4)  =: delta^(-nu)
e(a*)          ∝ delta^(s/2)   =  delta^(2*nu)        (line-core relation: theta = 2*nu)
```

The clean-merge bound s ≤ 2 gives ν ≤ ½ and θ = 2ν ≤ 1. The measured core-scale law is ν = 0.2 ± 0.1 as published ([M5.21.2b § 5](m5_21_2b_note.md), a = 2.5/2.9/3.4 at δ = 0.3/0.2/0.1); the audit's own least-squares refit of those three radii gives ν ≈ 0.27, so the working expectation is θ = 2ν ≈ 0.4-0.55 (effective s ≈ 0.8-1.1), a partial rearrangement either way. Two independent consistency anchors, cited not fitted: the 2-equal axis split closes ∝ δ^1.03 ([M5.21.1 P4](m5_21_1_method_note.md), linear-order core-structure corrections), and the archived M5.8-era statics was δ-flat to ∝ δ^−0.04 over 0.3 → 0.001 (older functional; supports a finite δ → 0 limit with weak total variation).

**The derived asymptotic form** (each branch, δ → 0):

```text
E_br(delta) = E∞_br * (1 + b_br * delta^theta + c_br * delta)
theta shared across branches (same local core physics),  0 < theta <= 1
```

The δ^θ term is the core-opening approach (θ = 2ν, expected ≈ 0.4); the linear term is the analytic far-field correction. At θ → 1 the two terms are degenerate; the frozen handling of that degeneracy is in § 3.

### 1.4 The g-axis: separability is derived, and the ladder is one-dimensional

**Exact statement.** The census functional (§ 1.1) contains no g: the T2 targets are (0, δ, 1) and the curvature term is parameter-free. Therefore on the instrument of record

```text
dE_3x3 / dg = 0        (identically, at every delta)
```

(the 4D-era g-blindness gate G8 measured this as exact; here it is structural). All g-dependence of a physical state enters through the 4D time row, and the measured channel is the boost-hedgehog dressing ([M5.21.8](m5_21_8_note.md)): the dressing energy E(m) is even in the boost parameter m with twin minima at ±m\*, where the author's law, verified to 0.009% and lattice-tracked at 0.82-0.84× across g = 8-64, is

```text
m*(g) = (1/2) ln((g+1)/(g-1)) = artanh(1/g) = 1/g + 1/(3g^3) + ...
```

**The gain law is an expectation, not a derivation** (adversarial-audit catch, § 9 claim 5, adopted). Evenness plus the m\* law do NOT fix how the dressing gain ΔE_dress = E(m\*) − E(0) scales with g: a stiffening-quartic family with E''(0) bounded gives gain ∝ m\*² ∝ g⁻², while a pitchfork family with E''(0) → 0 gives gain ∝ m\*⁴ ∝ g⁻⁴ (both audit-built with minima exactly at ±artanh(1/g)). The generic window is therefore

```text
E_br(delta, g) = E_br^{3x3}(delta) + ΔE_dress,br(delta, g),
ΔE_dress ~ -amplitude_br * artanh(1/g)^q   with q in [2, 4]   (to be MEASURED by the g-arm)
```

What IS derived: the correction is additive (the 3×3 sector is exactly g-free), and under any measured q ≥ 2 with the amplitude anchored at toy g, its relative size at the author-anchor g_phys is ≤ ~1e-20, negligible against every other term in the uncertainty budget for the energies and a fortiori for the ratios. **Consequence: the pre-registered ladder is one-dimensional in δ**, conditional on the g-arm verifying that the gain falls at least as fast as artanh(1/g)² (the § 6 F4 criterion; a faster fall, e.g. q = 4, PASSES and only strengthens the conclusion). The g-arm also carries the (−g)^p sign knob both ways (a ~1% 4D statics effect, [M5.21.1 § 4](m5_21_1_method_note.md); it never touches the 3×3 ladder).

### 1.5 The physical point (condition-1 specification)

The toy-to-physical specification, fixed from the author's own paper anchors and M5-side considerations only (no lepton input): δ_phys ~ 1e-10 (from the author's δ² ~ ħc hint) and g_phys ~ 1e10 (from the author's g⁴ ~ ke²/Gm² ≈ 1e38 hint; the audit notes both are ORDER-OF-MAGNITUDE anchors, (1e38)^(1/4) ≈ 3.2e9, and the prediction is insensitive to them: over g_phys ∈ [1e9, 1e11] and δ_phys ∈ [1e-11, 1e-9] every correction term moves by amounts far below σ). Functional: § 1.1 exactly (T2, w2 as pinned at calibration, w fixed across δ per the M5.12 pin). Lattice: N = 48, L = 48, h = 1, pinned bc, with the § 3 refinement subset.

**The extrapolation gap, stated plainly** (audit catch adopted): the rungs span δ ∈ [0.05, 0.30] (0.8 decades) while δ_phys sits 8+ decades below, where the core scale a\*(δ_phys) ≈ 200-2000 lattice units could fit no feasible box: E(δ_phys) is pure functional-form extrapolation from the measured decade, which is exactly why route (b) exists and why its output is claim-ceiling'd as a model-based prediction. Numerically the physical point is the E∞ limit (at θ ≈ 0.4, b·δ_phys^θ ~ b × 1e-4); the residual is carried, not dropped.

### 1.6 The target observables and outcome honesty

```text
R_C = E_C(delta_phys, g_phys) / E_A(delta_phys, g_phys)
R_B = E_B(delta_phys, g_phys) / E_A(delta_phys, g_phys)
```

with total uncertainties per § 6. Pre-committed reading, both directions: if the δ-dependence of the ratios is weak (as the derived form allows), the physical-regime ratios land near the toy-regime ones and the M8.6 comparison simply fails against its targets; that negative is reported under this same framework (condition 7). Nothing in this framework is steered toward any target value, and under the claim ceiling even a positive lands only as a model-based physical-regime prediction.

## 2. The frozen ladder protocol (D2): rung set and budget

| Element | Frozen specification |
| --- | --- |
| Rung set (δ) | {0.30, 0.25, 0.20, 0.15, 0.12, 0.09, 0.07, 0.05}, 8 rungs, each a genuine census run of all three branches at N = 48, L = 48, pinned, T2/sym, FIRE depth 12000 iterations (the census production depth) |
| Branch runs | 3 branches × 8 rungs = 24 production relaxations + § 3 refinement (18 additional) + the g-arm reads |
| δ-continuation | rung k + 1 seeds each branch from its rung-k endpoint; branch identity then verified by § 5 (never assumed from continuation) |
| Rung floor | δ_min = 0.05 sits above the expected instrument floor: the bulk eigen-gap ~ δ and the T2 gap guard held 0.012-0.038 at δ = 0.3; predicted core scale a(0.05) ≈ 3.6-4.1 (ν = 0.2 published / 0.27 audit refit) keeps L/a ≈ 12-13 (gate: L/a\* ≥ 10, measured per endpoint) |
| Box policy | L = 48 fixed across all rungs (fixed physical box); the L/a\* ≥ 10 gate excludes any rung the core outgrows; no mid-ladder box change |
| The g-arm | rigid Qb(m) dressing READS (evaluation only, no relaxation) on the N = 48 endpoints at δ = 0.3: g ∈ {8, 16, 32}, fine m-grids (the coarse-grid anti-recipe respected), both (−g)^p signs at g = 8, per branch; verifies gain ∝ artanh(1/g)² and measures κ_br |
| Affordability | MEASURED, not assumed (`scripts/m5_21_11_a_timing.py`, two short instrument-of-record relaxations, endpoints discarded, tags `t11timing_*`; `data/m5_21_11_timing.json`): 0.311 s/iter at N = 48 and 0.727 s/iter at N = 64 (ratio 2.34 vs the volume prediction 2.37), so one production rung ≈ 62 min, one N = 64 refinement run ≈ 145 min, and the whole frozen program (24 production + 9 N = 32 + 9 N = 64 runs) ≈ 50 CPU-hours, embarrassingly parallel across runs (an overnight batch at ~8 processes) |
| No silent extension | the rung set does not grow or shrink after freeze; a rung lost to its gates is EXCLUDED (recorded), never replaced |

## 3. The refinement subset and the discretization propagation rule (D2, per #378)

| Element | Frozen specification |
| --- | --- |
| Refinement subset | δ ∈ {0.30, 0.12, 0.05} (top, middle, bottom of the ladder) × N ∈ {32, 48, 64} at fixed L = 48 × all three branches = 27 runs (18 beyond the production rungs) |
| Per-rung continuum fit | on each refined (branch, δ): fit E_N = E_cont + k·N^(−p) through the three points (exact 3-parameter solve); σ_disc = \|E_cont − E_48\|. Stability guard (frozen): if the solved p falls outside (0.5, 6), the solve is deemed unstable for that (branch, δ) and σ_disc = \|E_64 − E_48\| is used instead (conservative direct difference). Recorded caveat: three points over a factor 2 in h is a zero-df solve; p is unvalidated by construction and the guard is the accepted handling |
| Propagation rule | for unrefined rungs: the relative discretization term σ_disc/E_48 is interpolated linearly in log δ between the two nearest refined rungs OF THE SAME BRANCH, then applied to that rung's E_48 (branch-dependent by construction: the heavier branches measured less consistency-converged than A) |
| Continuum correction | the ratio ρ_br = mean(E_cont/E_48) over that branch's three refined rungs multiplies E∞_br once at the end; the spread of the three values enters σ (§ 6). The fit itself runs on E_48 values |
| Affordability reading | the subset spans the ladder without touching an unaffordable top rung: per-iteration cost at N = 64 is timing-measured (§ 2); the propagation rule, not top-rung refinement, is what #378's amended wording requires (the scoping consequence recorded in the task doc), and no further amendment is needed under this reading. If the M8 side reads it otherwise, the one-clause amendment is raised BEFORE freeze |

**Fit degeneracy rule (frozen, deterministic).** The § 1.3 form is fit jointly: shared θ plus per-branch (E∞, b, c), 10 parameters. Determinism pins (audit catches adopted): θ is UNCONSTRAINED in the optimizer (the derived bound θ ≤ 1 is interpretive, never enforced as a fit constraint), and the θ interval is the PROFILE-LIKELIHOOD 68% interval (Δχ² = 1 on the θ profile), not a delta-method estimate. If that interval includes 1 (the δ^θ and linear terms degenerate), the frozen fallback is the single-correction form E∞_br(1 + c_br δ) refit deterministically; this rule is part of the framework, decided now, and is not a second framework.

## 4. Rung usability gates (frozen)

A rung's (branch, δ) result is USABLE only if all of: FIRE reaches f_tol at depth; cross-stencil consistency ratio ≤ 1.5; virial residual \|E_u − 3E_V\|/E_u ≤ 0.05; the eigen-gap guard silent on the endpoint read; L/a\* ≥ 10. A rung failing any gate is excluded with the failure recorded. These gates judge instrument health only; no gate consumes an energy comparison between branches.

## 5. The branch-identification rule (frozen; topology, never energy)

| Signature | Read |
| --- | --- |
| Electric charge class | surface-flux charge of the oriented long axis (the instrument of record, [canonical § 5.1](../m5_theory_canonical.md); quantized on the 1 − δ gap, which SURVIVES δ → 0) |
| Line census | the defect-line tracer (count + closure verdicts, [canonical § 5.4](../m5_theory_canonical.md)) with per-contour gap flags |
| Core spectrum class | which eigenvalue pair equalizes at the core |

Branch identity at every rung = the (charge class, line census, core class) triple matching the branch's rung-0.30 triple. **Energy ordering is never used for identification** (the M5.21.2b merge lesson). A MERGE is a defined outcome, not a failure: two branches landing within the merge metric of record (ΔE ≤ 0.04% and field distance ≤ 6%) are recorded as one state carrying both labels from that rung down; the framework verdict then goes through F3 (§ 6). Degeneracy of the reader at small δ is handled by the gap flags: a read with active flags is an EXCLUDED rung (§ 4), not a guessed identity.

## 6. The uncertainty model and the terminal failure criteria (frozen)

**Model.** Weighted least squares of the § 1.3 form on the usable non-holdout rungs, per-point weight σ_disc (§ 3). Parameter covariance → σ_extrap on E_br(δ_phys) by the delta method (the shared θ correlates branches; the ratio uncertainties use the full joint covariance). Total, per branch:

```text
sigma_br^2 = sigma_extrap^2 + sigma_rho^2 + sigma_g^2
```

with σ_ρ from the spread of the § 3 continuum-correction values and σ_g the (bounded, ~1e-20 relative) g-correction residual from the § 2 g-arm. Ratios R_C, R_B carry the propagated joint uncertainty: the condition-4 deliverable is (R_C ± σ, R_B ± σ), not an ordering.

**Holdouts (pre-registered, out of sample).** δ ∈ {0.20, 0.07}, all three branches (6 points): measured with the ladder but excluded from every fit; after the fit, each holdout is predicted with its 95% interval.

**Terminal failure criteria** (any ⇒ route (b) fails, M8.6 stays gated, no second framework):

| # | Criterion |
| --- | --- |
| F1 | joint fit quality: χ²/df > 3 under the frozen form (after the § 3 degeneracy rule, if triggered) |
| F2 | holdouts: ≥ 2 of the 6 points outside their 2σ prediction intervals, or any single point outside 3σ. Accepted risk, stated: with calibrated intervals this fires by noise alone with probability ≈ 3.3% (audit estimate); that false-failure rate is accepted as the price of a bright-line rule |
| F3 | branch integrity: any branch (or merged label-pair) with fewer than 6 usable rungs remaining after § 4/§ 5 exclusions and merges. Accepted marginal case, stated: at exactly the floor the fit runs thin (as low as df = 2 if all three branches sit at 6), making F1 noisy; the floor is kept as the bright line and the thinness is reported with the result |
| F4 | the g-arm: the measured dressing gain falls with g SLOWER than artanh(1/g)² across the g-arm rungs (audit-respecified: any fall at least as fast as artanh², including the pitchfork-like q = 4 case, PASSES; only a slower-than-quadratic fall breaks the § 1.4 negligibility bound and with it the one-dimensional-ladder reduction) |

**Reported diagnostics (pre-stated, non-terminal).** (i) θ̂ vs 2ν̂: the ladder measures a(δ) per rung; the line-core mechanism predicts θ = 2ν (§ 1.3) and the comparison is reported either way. (ii) Per-branch θ refits as a consistency read on the shared-θ choice. Neither diagnostic can alter the prediction or rescue a failure.

## 7. Barred-inputs audit (D3)

| Barred input | Status in this document and in the frozen procedure |
| --- | --- |
| The three existing toy census energies (N = 48) | absent from this document; never data points, exponent sources, or transformation anchors in any § 2-6 step (the ladder re-measures every rung fresh under the frozen spec, including δ = 0.30) |
| Charged-lepton mass ratios | absent; no target value appears anywhere in §§ 1-6 |
| Yukawa-derived 1 : 5.9 : 15.1 | absent |
| A→e, C→μ, B→τ assignment | consumed as FROZEN from the pre-existing stability/decay rationale; never re-derived here |
| Free E_physical = f(E_lattice) map | structurally impossible under §§ 2-3: every fitted quantity comes from new rung measurements |

Every numeric constant in § 1 is one of: an M5-internal measured law with its provenance link (core-scale ν, axis-split exponent, m\* law, virial residual, merge metric), a derived bound (θ ≤ 1, s ≤ 2), or an author-paper order-of-magnitude anchor for the physical point (δ_phys, g_phys: gravity/EM coupling scale and ħc, mass-independent; § 1.5 states their order-of-magnitude character and the prediction's insensitivity to them). No constant in any category descends from a lepton mass.

## 8. Mapping to the amended M8.6 § 8 preconditions

| § 8 requirement | Where satisfied |
| --- | --- |
| Asymptotic form derived from M5-side theory, independent of the lepton target | § 1 (with § 7 auditing independence) |
| Rung set, fitting, holdouts, branch rules, uncertainty model frozen before any performance evaluation | §§ 2-6, frozen at the § 0 freeze event, before any rung runs |
| Uniform application to A, C, B | §§ 2-5 are branch-uniform by construction (no per-branch tuning; branch-dependence only where derived, e.g. the § 3 propagation rule's branch separation) |
| Barred inputs | § 7 |
| Usable uncertainty on E_C/E_A, E_B/E_A (condition 4) | § 6 model |
| Per-rung discretization term + frozen propagation (condition 4 as amended, #378) | § 3 |
| Failure terminal | § 6 criteria; § 0 freeze mechanics |
| Claim ceiling | § 1.6: model-based physical-regime prediction, stated pre-commitment both ways |
| A→e/C→μ/B→τ frozen; no Yukawa figure | § 7 |

## 9. Adversarial audit record (D3)

Independent second agent, own algebra and own numeric checks throughout ([`m5_21_11_e_audit.py`](../scripts/m5_21_11_e_audit.py), results [`m5_21_11_audit.json`](../data/m5_21_11_audit.json)), run 2026-08-06 against the pre-fix draft. Verdicts:

| # | Claim | Verdict | Auditor's numbers |
| --- | --- | --- | --- |
| 1 | Derrick scaling + virial (§ 1.2) | ✅ CONFIRMED | measured exponents on a random smooth field: exactly −1.0 and +3.0; caveat carried: the dilation family is boundary-free, the pinned box breaks exact stationarity (handled by the § 4 virial gate) |
| 2 | Far-field δ-polynomial (§ 1.3) | ✅ CONFIRMED + strengthened | gap identity exact (err 0.0; path check 2.7e-11); quartic polynomial fit residual 9e-16; finite E(0); nonzero linear term. Strengthening ADOPTED: exactly degree ≤ 4, since M is linear in δ |
| 3 | Line-core competition, θ = 2ν (§ 1.3) | ✅ CONFIRMED | own minimization at s = 0.8: slopes −0.2000/+0.4000 (line), −0.2000/+0.2000 (point); the quartic-density premises verified for THIS functional (a Dirichlet functional would break the line premise) |
| 4 | artanh identity + series (§ 1.4) | ✅ CONFIRMED | exact to 3.6e-15 (g ≤ 64); series residual at g = 8 matches the next term to 1% |
| 5 | evenness ⇒ gain ∝ m\*² (§ 1.4 as drafted) | ❌ REFUTED | counterexample families with minima exactly at ±artanh(1/g): gain slope −2.00 (stiffening quartic) vs −4.00 (pitchfork); the m\*² law silently assumed E''(0) bounded. FIX ADOPTED: § 1.4 rewritten as a q ∈ [2, 4] window to be measured; F4 re-specified so a faster fall passes |
| 6 | g-freeness of the 3×3 instrument (§ 1.4) | ✅ CONFIRMED | no g in the instrument's configuration or source; structural |
| 7 | statistics (§§ 3, 6) | 🔶 PARTIAL | arithmetic confirmed (18 points, 10 params, df 8; 27 refinement runs); determinism gaps FIXED: profile-likelihood interval pinned, θ unconstrained in the optimizer |
| 8 | internal consistency numbers | 🔶 PARTIAL | a(0.05), L/a, δ_phys^0.4, artanh² all reproduced; catches ADOPTED: ν refit of the published radii ≈ 0.27 (θ expectation widened to 0.4-0.55), g_phys order-of-magnitude character stated, the 8-decade extrapolation gap stated plainly in § 1.5 |

Additional audit findings adopted: the § 3 three-point-solve stability guard; the F2 false-alarm rate stated as accepted (~3.3%); the F3 floor thinness stated as accepted; the § 7 closing sentence corrected (the physical-point anchors are author-paper anchors, not M5-internal measurements). Barred-inputs sweep: clean (no lepton ratio, no Yukawa figure outside its own barred-item label, no toy census energy anywhere).
