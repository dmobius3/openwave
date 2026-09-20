# M5 Summary Report: Results & Canonical Implementation

> ⚠️ **ARCHIVED 2026-07-14.** This is the M5.8-era results-of-record (last substantive update 2026-06-10), written on the PRE-M5.16 stack: the Eq.18 signed Hamiltonian + the saturating quartic `V_u = u + βu²` + the constrained spectral-projection integrator. The model of record has since moved to the M5.18-verified Lagrangian and its formulation program: the living registry is [`../m5_theory_canonical.md`](../m5_theory_canonical.md) (its § 4b carries this report's still-load-bearing anchors; its § 2 carries the historical-precedent note on the § 4.2 integrator). Numbers here remain the M5.8-era record; reproduction commands remain valid for that stack.
>
> **Index convention (2026-06-21).** The M5 engine now stores the order parameter at INDEX-0: `D = diag(g, 1, δ, 0)`, `eta = diag(-1, 1, 1, 1)` (time/g axis = array index 0, Duda's convention). Sections written before the flip use the legacy index-3 labels (`D = diag(1, δ, 0, g)`; the time-coupled curvature block called the `(α,3)` block) , the SAME physics under the relabel `index k -> (k+1) mod 4` (proven physics-neutral: [`_convention_refactor/CONVENTION.md`](../_convention_refactor/CONVENTION.md)). Read `(α,3)` as `(α,0)`.

**Purpose:** the results-of-record for the M5 Liquid-Crystal program (OpenWave). A narrative report of the 3+1D time-crystal field test, then every headline test as `topic | verdict + one-sentence finding + source script`, organized into the three validation bodies (**Liquid Crystal**, **Time Crystal**, **Zitterbewegung-clock existence**), then the **canonical implementation** (what works: the build that produced these results). The per-phase narrative lives in [`m5_roadmap.md`](../m5_roadmap.md) and the hardest-pieces board in [`m5_question_tracker.md`](../m5_question_tracker.md).

**Reading the verdict column:** the first line of each descriptive cell is the result: ✅ PASS / POSITIVE, ❌ NEGATIVE, ⚠️ NEGATIVE-INFORMATIVE or caveat. Scope throughout (unless noted): single defect, natural/lattice units, Duda's Eq.18-class matrix dynamics `M = ODOᵀ`, sandbox grids 24³/48³, f32 GPU with f64 numpy cross-checks.

**Orientation:** "M5" is OpenWave's implementation of Duda's liquid-crystal model; M5.7 / M5.8 are its phases (the resonance hunt; the Zitterbewegung clock). "N-1…N-6e" are the clock-measurement sub-experiments; the remaining labels (NG-, G-, Track C) index the project tracker and are not load-bearing for the results.

**Last updated:** 2026-06-08 (the full N-1…N-6e ZBW program, 2026-06-07).

---

## Report: the 3+1D time-crystal field test

Following exchange with Jarek Duda, this is the finished 3+1D field test of the time-crystal mechanism from arXiv:2501.04036. It is a results report: just what came out. Everything stated here is reproducible from the OpenWave repository (scripts plus per-experiment notes).

**Setup.** The field is the matrix `M = O D Oᵀ` with `D = diag(1, δ, 0, g)`, δ = 0.3, time as the fourth index. The dynamics is the Eq.18 signed Hamiltonian with the Minkowski (α,3) block. The only stable integrator found is a constrained spectral-projection scheme that keeps the positive-inertia directions per voxel and projects the momentum onto that subspace every step: every cheap positive-inertia kinetic tried has slow growing modes. All runs are headless on 24³ and 48³ grids (with a 63³ check on the runaway onset), f32 on GPU with f64 numpy cross-checks (they agree to about 1% everywhere checked), and a full-Euclidean twin ran as a kill-control at every stage.

**The quadratic action does not saturate at 3+1D.** The bare `−ΣF²` action runs away at a fixed phase, about two clock periods, and the onset is dt-invariant (it scales exactly with dt under halving), identical in f64, f32 and Metal, grid-independent across 24³ and 63³, and independent of the potential and the momentum clamp. That the runaway onsets at a fixed phase rather than a fixed time indicates a phase-locked (parametric) instability, not a fixed-clock numerical limit. This is the field-level version of the toy model's own requirement for a βR⁴ floor.

**A saturating quartic restores bounded dynamics.** Adding `V = u + β u²` on the signed spatial-curvature density `u = 2 Σ_(i<j) ⟨F_ij,F_ij⟩_s`, with β fixed by `2β|u_min| ≈ 3`, bounds the runaway: H breathes over the full horizon, the floor-bounce is real, and it is dt-converged and f64-confirmed.

**The signature is what drives it, not the nonlinearity.** With the same seed and kick, the Euclidean twin's clock decays while the Minkowski-plus-quartic clock grows by about 3×. The (α,3) block is doing the work.

**The saturated state self-starts.** A configuration damped to rest, restarted with momentum exactly zero, regrows kinetic energy from zero (to 5.76 by τ = 4000, H conserved to 0.5%), reproduced at half dt to four significant figures: a field starting at exactly zero momentum develops motion at conserved energy. The collective-coordinate analysis behind this finds the dressed minimum to be un-sittable: energetically minimal but carrying ghost breathing inertia (K_bb < 0): which reads as the 3+1D face of the negative-energy auto-propulsion in Fig 10 of the paper.

**Its frequency is a property of the state.** The breathing fundamental is the same: within one FFT bin: whether the run is kicked, exactly unkicked, or seeded with small random jitter. One correction belongs here: an earlier reading claimed a strictly periodic ω₀ = 0.262 with an exact harmonic comb; a window-length sweep showed that was an FFT-window artifact: the peak moved with the window length and the comb was just bin arithmetic. The state carries a resolved fundamental plus a 2ω harmonic over a broadband, aperiodic background: a coherent comb, not a strictly periodic line (the background is classified below as low-dimensional chaos).

**The defect holds, and the holding is resolution-robust.** This is the M5.7 dispersal question asked again under the quartic with full backreaction. The structural alignment decays slower at 48³ than at 24³ at every matched time, settling toward a plateau well above the random baseline, and the result does not depend on the masking. This is the reverse of the M5.7 free-defect behaviour, where the wash-out strengthened with resolution.

**The frequency belongs to the core, not the dressing.** Across a family of states spanning a 2.6× range in rest energy (varying the dressing width), the fundamental is constant to within one bin: the exponent of ω against energy is about 0.03, not 1. So the naive ω ∝ mass law does not hold for the dressing energy; the frequency is set by the core, which the V = 0 frame cannot vary because it is scale-free.

**The absolute number, and where it falls short.** Under two explicit postulates: the rest energy mapped to 0.511 MeV and the lattice action unit mapped to ℏ: the clock runs at about 5.5×10¹⁹ rad/s, roughly 28× below the electron Zitterbewegung 2mc²/ℏ. Because the frequency is rigid against energy, that gap is structural rather than something closable by energy bookkeeping; it points at physics this V = 0 sandbox does not contain: the V-on/Faber-r₀ core, the faithful kinetic, or the action-to-ℏ postulate itself. This is reported as the first measurement, not as a pass or fail.

**Two more readings.** A chaos battery classifies the saturated state as a molten clock: a persistent coherent fundamental on low-dimensional chaotic dressing whose intensity grows with excitation: and the deep-settled cold state reads near-regular, so the clock appears to tick cleanly toward the ground state and melt as it heats. Full coldness is forbidden by the same un-sittability that makes it self-start. Separately, the rotation Noether charge of the clock kick is zero: the Θ-twist cancels over the hedgehog sphere, so the breathing channel carries no net frame angular momentum at this scale: though box torque on the 24³ grid bounds that measurement, so it is not read as a spin-½ statement either way.

**What is not claimed.** No particle-stability result; the state is a molten clock, not a strict single-line oscillator. The absolute frequency is 28× off. There is no intrinsic spin above the box-torque floor. It is a single defect in natural units, and the explicit stepper stiffens near the quartic floor so the very-late-time horizons are numerically limited.

**Standing.** The repository reproduces everything above. The 2+1D moving-defect (pilot-wave) test is the natural next step.

The detailed per-test evidence: every claim above plus the Liquid-Crystal substrate validations and the canonical build: follows.

---

## 1. Liquid-Crystal validations (the substrate is faithful)

These establish that the matrix substrate reproduces the established physics Duda's framework requires: charge, forces, EM, geometric mass. All re-verified on the production 4×4 engine after the M5.8.1 promotion (the reproduction of the original 3×3 numbers is itself proof the promotion is physics-preserving).

| Topic | Result · finding · source |
| --- | --- |
| Topological charge | ✅ POSITIVE<br>Charge is an integer, additive winding of the director: Q = ±0.996 for single hedgehogs, ±1 around each core of a pair and 0.000 on the enclosing sphere: quantization and additivity both emergent, not imposed.<br>`m5_8_1_topo_charge.py` |
| Coulomb force (1/d) | ✅ POSITIVE<br>Opposite-sign defects attract with a 1/d law fit at R² = 0.9781 (b < 0), matching the M5.1 reference (0.978) digit-for-digit on the matrix substrate.<br>`m5_4_coulomb_matrix.py` |
| Coulomb: Duda page-18 cross-check | ✅ POSITIVE<br>Duda's own Mathematica page-18 configuration reproduces an attractive 1/d fit at R² = 0.9959: an independent geometry confirming the force law.<br>`m5_4_coulomb_page18.py` |
| Maxwell from tilts (hydro↔EM) | ✅ PASS<br>Defining E = ∇·n̂ and B = ∇×n̂ recovers the full Maxwell set: ∇·B = 0, the Gauss charge identity, Faraday, and the Lorentz force: to machine precision via the hydro↔EM dictionary.<br>`m5_6_4a_hydro_em.py` |
| Maxwell from tilts (Faber curvature) | ✅ PASS<br>The Faber field strength R = Γ×Γ is a closed 2-form (Maurer-Cartan ⇒ homogeneous Maxwell), giving abelian Coulomb ‖R‖ ~ 1/r² far-field with a running-coupling onset at the core radius r₀.<br>`m5_6_4b_faber_curvature_em.py` |
| KG mass is GEOMETRIC | ✅ POSITIVE<br>The Klein-Gordon mass emerges from minimal coupling to the hedgehog connection Â (Fig.9 reproduced): the explicit mass term cancels exactly, so mass is geometry, not an added V_ψ.<br>`m5_6_1_kg_operator_check.py` |
| Biaxial hedgehog sources its own twist | ✅ POSITIVE<br>The biaxial defect dynamically generates C_μν ≠ 0 with a ‖C‖ ~ 1/r² restoring-mass profile: it cannot sit static, which is the seed of the M5.8 clock.<br>`m5_6_2a_biaxial_hedgehog.py` |
| Faber mass calibration (E ∝ 1/r₀) | ✅ POSITIVE<br>Porting Faber's MTF regularization pins the defect mass to the core radius, E·r₀ = const (CV = 0.0%), giving the physical knob that maps r₀ = 2.2132 fm → 0.511 MeV: the M5.9 calibration handle.<br>`m5_6_3a_faber_*` / `m5_6_3b_*` |
| Eq.18 matrix evolution conserves energy | ✅ PASS<br>The Duda Eq.18 leapfrog conserves energy: the drift falls 2.15% → 1.13% → 0.03% as dt → 0 (O(dt²)), confirming the production EOM conserves the Hamiltonian to integrator order.<br>`m5_5_4_matrix_evolution_check.py` |
| LdG potential confinement | ✅ PASS<br>The Landau-de Gennes V(M) confines the amplitude Tr(M²) ~3.3× across k ∈ [0.5, 25] with no blow-up: the well is stable across two decades of stiffness.<br>`m5_6_5c_prod_scale.py` |
| Free 3D defect disperses (M5.7) | ❌ NEGATIVE (expected: Derrick)<br>Both a seeded l=1 perturbation and the defect's own intrinsic oscillation disperse their orientation energy in pure 3D (localization washed out 32³→48³): confirming a static/oscillatory 3D particle is forbidden and the stable object must be 4D.<br>`m5_7_1_l1_resonance_seed.py` · `m5_7_2_intrinsic_oscillation.py` |
| Driven defect sustains an (A,ω) excess (M5.7.3) | ✅ POSITIVE<br>A continuous EM-like drive holds a bounded, frequency-selective (A,ω) excess at ~3× the free baseline at the resonant mode (N=48-confirmed): a resonant drive holds the defect's (A,ω) state in a maintained excess.<br>`m5_7_3_driven_oscillation.py` |

---

## 2. Time-Crystal validations (the mechanism works)

These are the M5.8 core: a time-periodic resonance: bounded, self-starting, signature-driven: exists in 3+1D. The arc is "quadratic fails → quartic saturates → the saturated state is a self-starting attractor."

| Topic | Result · finding · source |
| --- | --- |
| 1D toy-model clock (Duda 2501.04036) | ✅ POSITIVE<br>The 1+1D time crystal reproduces to 4-5 digits (ω\* = √(70/61) = 1.0712, E* < E(0)) with machine-precision energy conservation, and the GHOST finding: the ψ-sector is negative-kinetic ⇒ a constrained integrator is mandatory: set the rule for everything 4D.<br>`m5_8_0b_*` · `m5_8_0c0d_propulsion_robustness.py` |
| 4D Minkowski clock fuel is real | ✅ POSITIVE<br>The exact 4D Hamiltonian gives E(ω,b) = A(b) + ω²·C(b) with C(b) < 0 for every clock-plane × boost-axis once the time axis is dressed: and exactly zero undressed (the M5.7 null in functional form); the Euclidean flip kills it, so the signature IS the mechanism.<br>`m5_8_2a_4d_hamiltonian.py` |
| The quadratic action does NOT saturate | ❌ NEGATIVE (decisive)<br>The bare `−ΣF²` action runs away at a fixed phase τ ≈ 2 clock periods: dt-INVARIANT (ratio exactly 2.0 under dt-halving), f64 ≡ f32 ≡ Metal, 24³ ≡ 63³, V/clamp-independent: the field-level twin of the 1D toy's own requirement for a βR⁴ floor.<br>`m5_8_2cb_taichi_constrained.py` |
| The saturating quartic restores bounded dynamics | ✅ POSITIVE<br>`V_u = u + β·u²` on the signed spatial-curvature density (β = 1.558, data-anchored) bounds the runaway: H breathes over the full horizon (≥8 fundamental periods), floor-bounce real, dt-converged, f64-confirmed.<br>`m5_8_2d_quartic_saturation.py` |
| The signature IS the engine (kill-control) | ✅ POSITIVE<br>Same seed + kick: the Euclidean twin's clock DECAYS (1.4e-2 → 6.7e-3) while the Minkowski+quartic clock GROWS 3× (7.6e-3 → 2.4e-2): the Minkowski signature, not generic nonlinearity, sustains the time crystal.<br>`m5_8_2e_invariant_matrix.py` (euclid mode) |
| Static reductions all closed (Track C) | ❌ NEGATIVE (decisive, ×3)<br>The global-clock BVP has no interior ghost minimum (C1), the static clock-phase twist is pure cost (C2-B), and the dressed minimum is UN-SITTABLE: K_bb(a*) = −67.6 < 0 with U″ > 0 (C3): so the breathing state is irreducibly time-dependent: ḃ ≠ 0 IS the bounce.<br>`m5_8_2f_breathing_bvp.py` · `m5_8_2f2_localized_clock.py` · `m5_8_2f3_breather_orbit.py` |
| Spontaneity confirmed at field level | ✅ POSITIVE<br>A damped-settled configuration restarted with P = 0 EXACT regrows kinetic energy 0 → 5.76 by τ = 4000, H conserved to 0.5%, dt/2-converged to 4 significant digits: from exactly zero momentum the field develops motion at conserved H: the clock self-starts.<br>`m5_8_2g_spontaneity.py` |
| The invariant matrix is COMPLETE | ✅ PASS (all candidates classified)<br>Signed-u is the working invariant; Skyrme `β_E·u_E²` also saturates but damps the clock 10× (sign-blind); the M6-style amplitude quartic A4 is geometrically pinned (no lever arm: the fuel is curvature-sector); the covariant `𝒮+β𝒮²` is deferred-with-reason (cubic Legendre, static sector ≡ the validated quartic).<br>`m5_8_2e_invariant_matrix.py` · `m5_8_2l_invariant_completion.py` |
| The defect HOLDS, resolution-robust (N-2) | ✅ POSITIVE (G-2c-1)<br>Under the 4D quartic stack with full backreaction the free defect's structural alignment decays SLOWER at 48³ than 24³ at every matched t: the anti-M5.7 signature (its wash-out strengthened with resolution; this reverses): settling to a plateau 0.8-0.88 ≫ the 0.5 random baseline, mask-insensitive.<br>`scripts/m5_8_2i_dispersal_gate.py` |

---

## 3. Zitterbewegung-clock existence + identification (N-1…N-6e)

The split that matters: **existence** (rows 1-6, all ✅: a self-starting clock with a frequency that is a property of the state) vs **identification with the electron** (rows 7-10: the calibration program, where the gaps are). "Does the ZBW clock exist?": yes. "Is it the electron's?": not yet, and the gap is localized.

| Topic | Result · finding · source |
| --- | --- |
| ω is an ATTRACTOR (N-1) | ✅ POSITIVE<br>Kicked, exactly-unkicked, and jittered starts all settle to the same breathing fundamental ω₁ ≈ 1.1 + a 2ω₁ harmonic (detrend-stable across filter windows): the frequency is a property of the STATE, not the preparation.<br>`scripts/m5_8_2h_omega_attractor.py` |
| The "0.262 comb" was an FFT artifact (N-1) | ⚠️ NEGATIVE-INFORMATIVE (self-caught)<br>A window-length sweep showed the earlier "strictly periodic ω₀ = 0.262 + exact comb" was an FFT-window artifact (the peak moves with window length; the comb was bin arithmetic): the real state is a resolved fundamental over a broadband, aperiodic background (classified below as low-dimensional chaos, N-6d).<br>`scripts/m5_8_2h_omega_attractor.py` (ref mode) |
| The breather is NOT a disclination-rod mode (N-1) | ✅ POSITIVE<br>The motion is near-isotropic (anisotropy A_z ≈ 1.0-1.3 vs the static rod baseline 5.58) and core-centered (3× core-shell concentration early), with no excess on the z-axis rod: the clock lives on the defect core, not the line.<br>`scripts/m5_8_2h_omega_attractor.py` (W4) |
| First ZBW ratio recorded (N-3) | ✅ MEASURED (no pass/fail)<br>ω₁/(2H_rest) = 0.0326-0.0343 across arms (5.4% spread, start-independent) in lattice-natural units: the first measurement of the ratio for the saturated breather, with the factor-2 apolar-doubling bookkeeping respected (no double-correction).<br>`scripts/m5_8_2j_zbw_ratio.py` |
| ω is RIGID across a mass family (N-6a) | ✅ POSITIVE (a property of the core)<br>Across a knob-gated dressing-width family spanning 2.6× rest / 3.9× state energy, ω₁ is constant within one FFT bin (exponent ω ∝ H^0.03 vs the naive law's 1): the frequency belongs to the defect core, which the scale-free V=0 frame cannot vary.<br>`scripts/m5_8_2m_zbw_law.py` |
| Matter + light share a radial cone (N-4) | ✅ POSITIVE<br>The twist (KG) and both tilt (EM) channels hit the same maximal radial signal speed c(r̂)/2 = 1.0000 at N=48 and 64: matter and light propagate on a common radial cone on the defect background (a necessary condition for emergent local Lorentz invariance, not a full demonstration); the tangential cone shapes differ per channel (one tilt polarization rank-1, non-propagating along one direction).<br>`scripts/m5_8_2k_tilt_cone.py` |
| Molten clock classification (N-6d) | ✅ POSITIVE (classified)<br>The saturated breather is a MOLTEN CLOCK: a persistent ω₁ comb riding low-dimensional chaotic dressing (λ_max +0.4-0.7/τ, D₂ 2.7-3.0, control-anchored) whose intensity grows with excitation: the DTC literature's middle "melting" phase: a coherent clock dressed by low-dimensional chaos, not destroyed by it.<br>`scripts/m5_8_2n_chaos_battery.py` |
| The clock regularizes toward the ground state (N-6c) | ✅ POSITIVE<br>The coldest measured state (deep-settle, T = 2.36) reads near-regular (λ_max +0.11 ≈ the estimator bias floor, D₂ 1.68 vs 2.7-3.0 hot): cold = ticking clock, hot = molten: and full coldness is FORBIDDEN by the un-sittability that drives the self-start (the clock and its excitation are one object).<br>`scripts/m5_8_2o_omega_of_E.py` |
| First absolute ω (N-6b) | ⚠️ NEGATIVE-INFORMATIVE (structural gap)<br>Under two explicit postulates (H_rest ↔ 0.511 MeV; lattice action ↔ ℏ) the clock runs at 5.5×10¹⁹ rad/s: ~28× below the electron ZBW 2m_ec²/ℏ = 1.55×10²¹: and because ω is rigid the gap is STRUCTURAL (not closable by energy bookkeeping), pointing at the V-on/Faber-r₀ core, the faithful kinetic, or the action↔ℏ postulate.<br>`scripts/m5_8_2j_zbw_ratio.py` (N-6b addendum) |
| First spin readout (N-6e) | ⚠️ NEGATIVE-INFORMATIVE (clean null + bound)<br>The rotation Noether charge (gated ⟨P,Ṁ⟩/(2T) = 1.000000 exact) reads J = 0 on the clock kick: the Θ-twist cancels over the hedgehog sphere, so the breathing channel carries no net frame angular momentum; the measurement is bounded by box torque at 24³ (an ℏ/2-class result needs torque-free boundaries), and spin-½, if carried, lives in the polarized-seed / hopfion / far-field sectors.<br>`scripts/m5_8_2p_spin_readout.py` |

---

## 4. Canonical implementation: what works

The build that produced every result above. This is the load-bearing configuration; anything not listed here is either superseded (the ψ engine, retired) or backlog (the V-on core, the faithful kinetic).

### 4.1 The substrate + action

| Element | Canonical form |
| --- | --- |
| Field | Matrix `M = O D Oᵀ`, `O ∈ SO(3)` embedded 4×4 (time = index 3), `D = diag(g, 1, δ, 0)`, δ = 0.3: biaxial, NOT Vector(3) ψ (the M5.2 closed negative). |
| Action | Duda Eq.18 `ℒ = −Σ F_μν F^μν − V(M)`, `F = [∂_μM, ∂_νM]`; the 4D clock sector uses the signed (Minkowski `(α,3)`-block) flux `F → ηFη`. |
| Saturating invariant | `V_u = u + β·u²` on the SIGNED spatial-curvature density `u = 2Σ_(i<j)⟨F_ij,F_ij⟩_s`, β = 1.558 (anchored 2β·\|u_min(seed)\| ≈ 3): the 1D βR⁴ floor, transplanted to the channel that runs. |
| Potential | LdG `V(M)` on the SPATIAL 3×3 block only (the g axis is decoupled: the M5.8.1 spatial-block restriction); Faber regularization for the r₀ mass scale. |
| Seed | Boost-dressed biaxial hedgehog `O = O_hh(x)·exp(b·w(r)·B)·R(ωt)`, w = exp(−(r/R_W)²), R_W = 3.5; the dressed minimum (NOT bare hedgehog + kick) is the constraint-manifold seed. |

### 4.2 The integrator (the only stable scheme found)

| Element | Canonical form |
| --- | --- |
| Kinetic | The faithful F-commutator metric at CC level; production ships the simple `½‖Ṁ‖²` (frequency-corrected only: the 5d diagnosis). |
| Stepper | The constrained spectral-projection integrator: per-voxel 10×10 inertia operator, keep positive-inertia directions (λ > 0.05·max\|λ\|), project P onto the kept subspace EVERY step, global-(α,3)-momentum clamp. Every cheap positive-inertia kinetic has slow growing modes: this is the one that holds. |
| Precision | f32 on GPU (Metal), f64 numpy probes as the independent cross-check: f32 ≡ f64 to ~1% everywhere checked; the deep-floor cascade (the late-time stiffness limit) is the known numerical bound, kept off the load-bearing windows by dt-convergence. |
| Guards | No free global dressing amplitude (the 2b ghost channel); bounded-energy monitor (mass can't go ≤ 0); momentum-aware det𝕂 zero-mode guard. |

### 4.3 The validated method discipline (earned, reusable)

| Rule | Why it exists |
| --- | --- |
| Surrogate GUIDES, direct quadrature DECIDES | Track C: spline/4-pt surrogates carry mesh artifacts; only machine-exact direct evaluation is load-bearing. |
| The dt-invariance discriminator | Real dynamics reproduces at matched τ under dt-halving; stepper-driven growth shifts with dt: the test that killed the quadratic-runaway-is-numerical hypothesis. |
| The exact 9-pt phase average | Densities have bounded harmonics (V,K ≤ e^{i8ψ}, Q ≤ e^{i16ψ}) ⇒ 9 uniform phase points average exactly; the 4-pt set is not shift-invariant under twist. |
| FFT-window discipline | No peak within ~2 bins of DC; combs at bin multiples are automatic; vary the detrend window as a dial; pair signed probes with rectifying envelopes (caught the 0.262 artifact). |
| The knob gate | Validate a parameter family actually changes the seed (fixed-domain invariant spread > 5%) BEFORE spending compute (caught the RC no-op). |
| Headless-first validation | Every M5.8 result was produced headless (gates + npz + plots); rendering gates nothing scientific (the 2026-06-07 decision). |

### 4.4 What is NOT in the canonical build (scoped out)

| Item | Status |
| --- | --- |
| The electron mass/ω calibration | NG-1: needs the V-on/Faber-r₀ core (the 28× absolute-ω gap lives here); the ratio infrastructure exists, the units fix does not. |
| The faithful kinetic in production | NG-1 residual: the 5d diagnosis says it corrects the clock FREQUENCY only; re-measure ω₁ under it when production-scale clock runs need it. |
| Spin-½ as a measured charge | Open: the clock channel is J-neutral; spin-½ (if carried) needs the polarized-seed / hopfion / far-field sectors + torque-free boundaries. |
| Multi-defect / composites | 15a: everything above is single-defect; ensemble behavior is the many-body question. |
| Strict single-line periodicity | Not claimed: the state is a molten clock (a coherent comb on low-dimensional chaotic dressing), regularizing toward the cold ground state. |

---

## Reproduction commands (the evidentiary index)

Every script above lives under `research/scripts/` (the M5.8.0-2g arc) or `research/scripts/` (the N-1…N-6e ladder). The full command table: prerequisites, run times, expected outputs: is in [`m5_roadmap.md`](../m5_roadmap.md) § *Reproduction commands*. Each script carries its own PREREQUISITE + RESULTS docstring with the caveats inline.

**M5.8.2c/d/e/f REPRODUCTION COMMANDS** (from `research/scripts/`; all results above re-derivable from these: scripts carry PREREQUISITE + RESULTS notes in their docstrings):

| Result | Command | Expected |
| --- | --- | --- |
| Seed npz (PREREQUISITE for 2d/2e) | `CB_STEPS=900 python m5_8_2cb_taichi_constrained.py ref` (~5 min) | writes `_m5_8_2cb_ref.npz` (uncommitted derived data; seed arrays run-length-independent) |
| Option B spike gates J1,B1,B2,B2m,B3,B4 | `python m5_8_2cb_taichi_constrained.py all` (~15 min) | PASS: machine-exact vs f64, 14.5 ms/step at 64³ |
| The quadratic runaway (any precision) | `CB_STEPS=6000 python m5_8_2cb_taichi_constrained.py ti metal f32` / `... ref` (f64, ~38 min) | onset ~1150, H → −8.6×10⁹, align → 0.55 |
| The dt-invariance (fixed-τ) proof | `CB_STEPS=12000 CB_DT=0.001 python m5_8_2cb_taichi_constrained.py ti metal f32` | onset ~2350 (ratio exactly 2.0 vs dt=0.002) |
| Production signed-GUI repro (4 variants) | `python m5_8_2cb2_signed_gui_repro.py` (63³ Metal, ~10 min) | all bounded at the derated dt; align probe shows the runaway timing |
| 2d quartic scan (D1-D5) | `python m5_8_2d_quartic_saturation.py` (~5 min) | D1/D2 PASS; β ladder −156/−39/+38.8; ⚠️ read the per-β table, not the align-0.90 auto-headline |
| 2d long-horizon + bounce | `python m5_8_2d_quartic_saturation.py long 1.558 24000` / `long 1.558 48000 0.5` (dt/2) | bounce at τ≈32-36; dt/2 bounded end-to-end (full-dt late heating = stepper stiffness) |
| 2e Skyrme / kill-control / ω / f64 | `python m5_8_2e_invariant_matrix.py skyrme\|euclid\|omega\|f64` | Skyrme saturates (clock 10× weaker); Euclid clock decays vs Mink grows; f64 no-runaway. ⚠️ its `omega` print (ω₀=0.262 + comb) is the RETIRED window artifact: the valid ω measurement is the 2h row ↓ |
| Production headless (33 gates) | `python m5_8_1_headless_check.py` (repo root, ~2 min) | 33/33 PASS |
| 2f Track C C1 (the BVP attack, ansatz verdict) | `python m5_8_2f_breathing_bvp.py all` (~70 s; tables auto-cached to `_m5_8_2f_tables.npz`) | hard gates PASS (canary 0.00%, F1a −48.8408 exact); NO interior ghost minimum: the slice rises 8.4 → 1.5×10¹⁰; 2b E* refined to 2.93 direct |
| 2f2 Track C C2-B (the twisted-frame clock) | `python m5_8_2f2_localized_clock.py all` (~50 s + one-time ~110 s 9-phase tabulation) | hard gates PASS (t-canary 0.01% post-9-pt-fix); NO twisted state below the t≡0 control; the marginal ΔH=−0.002 candidate self-refutes via the built-in twist-differential dial |
| 2f3 Track C C3 (the reduced breather dynamics) | `python m5_8_2f3_breather_orbit.py all` (~36 s + one-time ~110 s kinetic tabulation) | all gates PASS; K_bΘ≡0 exact; the un-sittable minimum (K_bb(a*)=−67.6 direct); H3 = reduction boundary at det 𝕂 = 0 (compulsory motion, containment many-mode) |
| 2g the spontaneity test (the field handoff) | `python m5_8_2g_spontaneity.py` (4×12k, ~6 min) then `settle 24000 16000` (~5 min) then `restart 0.5 8000` (~1 min) | unkicked = kicked end-state with \|s\|~1e-18 (blindness); S5 from settled, P=0 exact: T 0→5.76 at τ=4000, dt/2-converged to 4 digits: SPONTANEITY CONFIRMED |
| 2h the ω-attractor + rod readout (from `research/scripts/`) | `python m5_8_2h_omega_attractor.py ref` (free, saved-data W1 scrutiny) then `run 48000` (4 arms dt/2, ~22 min) then `analyze` (re-print from npz) | W1: the 2e slow "ω₀" moves with window length (artifact); W2: ω₁ ≈ 1.09/1.15/1.07 + 2ω₁ across kicked/unkicked/jittered (detrend-stable): ATTRACTOR; W4: A_z ≈ 1.0-1.3, frac_rod = volume fraction (not rod-riding); rig checks: dM 0.2501 @τ-equiv 4000 = 2g to 4 digits |
| 2j the first ZBW ratio (from `research/scripts/`) | `python m5_8_2j_zbw_ratio.py` (CPU-only, seconds; needs the 2h npz) | ω₁/(2H_static) = 0.0326-0.0343 across arms (5.4%); H_static = 16.74; H_dyn/H_static ≈ 2.7-2.8 |
| 2k the EM-tilt cone check (from `research/scripts/`) | `python m5_8_2k_tilt_cone.py` (CPU-only, ~3 min) | shared radial ceiling c(r̂)/2 = 1.0000 all channels; twist cone 1.743 N-stable; Gy near-rank-1; Gz rank-1 exact |
| 2i the G-2c-1 dispersal gate (from `research/scripts/`) | `python m5_8_2i_dispersal_gate.py run 48000` (24³ 3 arms ~13 min) then `M58_N=48 ... run 48000` (~3.5 h; seed: `M58_N=48 CB_STEPS=2 ... ref`) then `analyze` | align decays SLOWER at 48³ (anti-M5.7) → plateau 0.8-0.88; mask-insensitive; G-2c-1 ✅; 48³ cascade at t≈13 |
| 2l the invariant completion (from `research/scripts/`) | `python m5_8_2l_invariant_completion.py classify` (CPU sympy, instant) + `a4 6000` (~3 min GPU) | covariant: P cubic in q̇ ⇒ deferred; A4: dev p95 = 0.0000 ⇒ amplitude channel pinned, no lever |
| 2m the ZBW-law family (from `research/scripts/`) | `python m5_8_2m_zbw_law.py run` (knob gate + 3 × ~33 min) | knob gate 177.8% spread; ω₁ = 1.152/1.188/1.191 across H 11/16.7/29.1: RIGID (exponent ≈ 0) |
| 2n the chaos battery (from `research/scripts/`) | `python m5_8_2n_chaos_battery.py` (CPU, ~3 min; needs the 2h npz) | controls calibrate the battery (K unreliable both ways); arms: MOLTEN CLOCK (λ +0.4-0.7, D₂ 2.7-3.0) |
| 2o ω(E) + the ground clock (from `research/scripts/`) | `python m5_8_2o_omega_of_E.py run` (~70 min) then `settle 48000` + `restart 48000` (~11 + 35 min) | excitation floor ~1.7× H_rest (kicks can't cool); ground readout: λ +0.110, D₂ 1.68 at T = 2.36: REGULARIZES |
| 2p the first spin readout (from `research/scripts/`) | `python m5_8_2p_spin_readout.py 6000` (~2 min; needs the 2o settled npz) | G1 pairing 1.000000 exact; the clock kick is J-NEUTRAL (J < 1e-4 at t→0); G2 = box torque (secular J growth: the 24³ walls) |

## DUDA 2026-06-08 FOLLOW-UP: the δ / g calibration

**Full correspondence record:** [`m5_4c_convo_2026.06.08.md`](../tasks/m5_4c_convo_2026.06.08.md) (the email thread + both rounds of runs). This section is the findings summary.

After the report, Duda read `m5_summary_report.md` (now seeing the 3+1D **4×4** tensor sim correctly: `D = diag(1, δ, 0, g)`) and added Manfried Faber. His directions, with the QED Lagrangian he attached as the decoder ring:

| QED term (his image) | coefficient | our `D` axis | his prescription |
| --- | --- | --- | --- |
| Dirac kinetic `ψ̄γ∂ψ` (quantum phase) | `ℏc` ("tiny") | `δ` (minor), enters as `δ²` | `δ ~ 10⁻¹⁰`, not 0.3 |
| mass `mc²ψ̄ψ` (→ gravity) | `mc²` | `g` (time axis) | `g ~ 1/δ ~ 10¹⁰`, ideally `g·δ = 1` |
| gauge `F²/4` (EM) | `1/4` | the `1` (major) | order unity: we have it |

Plus: fix units by comparing to Coulomb or the clock (ours are dimensionless); tune the LdG potential to the particle rest energies. He found the 5.5×10¹⁹ rad/s clock "surprisingly close: especially for using much too large δ."

**What we ran** (`scripts/m5_8_2q_delta_scaling.py`: seed-level, exact; gate reproduces H_static = 16.74 at δ=0.3, g=8):

| Finding | Result |
| --- | --- |
| `g=1/δ` (his hierarchy) | H_static **diverges** `∝ δ^(−6.8)` (→9×10¹⁶ at δ=0.001): the literal physical scale is numerically impossible, matching his own "too large for full simulations → neglect gravity" |
| g-sensitivity at δ=0.3 (⚠️ SUPERSEDED, see 2026-06-09 addendum) | raising `g` at FIXED boost gives `∝ g⁸`, first read as "gravity-dominated." Duda corrected this 2026-06-09: gravity enters ONLY through the boost tilt `b·g`, not the eigenvalue `g`. Raising `g` at fixed boost is an unphysically large tilt, not gravity. The physical picture is in the addendum below |
| δ-sweep, gravity **decoupled** (g=8 fixed) | H_static is **δ-flat**: 16.74→20.64 across δ 0.3→0.001 (fit `δ^(−0.04)`): the quantum-phase δ-axis is only a ~23% correction to a g/major-axis core |

**The clock-ratio consequence**: combining the measured δ-flat `H_rest` with the N-6a ZBW law (`ω ∝ H_rest`, so `ω/2H ≈ 0.033`):

`R(δ) = ω·δ/(2·H_rest) ≈ 0.033·δ` (the m5_8_2j ℏ↔δ convention). At δ=0.3 → 0.0099 ✅ (matches the recorded 0.0098-0.0103); R=1 would need δ≈30; δ=10⁻¹⁰ gives R≈3×10⁻¹². **So the δ knob does NOT calibrate the clock at fixed functional: smaller δ drives the ratio AWAY from 1.** The δ=0.3 "28× close" is set by the two unit postulates (P1: H↔0.511 MeV, P2: action↔ℏ), not a δ-calibration, and does not improve at physical δ.

**Takeaway:** the calibration lives in Duda's *other two* directions: fix units via Coulomb (calibrate on one observable, predict the clock) and tune the LdG potential to rest energy (the Faber-`r₀` handle, NG-1): **not** in δ at fixed functional. Gravity must be decoupled (`g=1/δ` diverges), which independently confirms "put clock propulsion by hand."

**Reproduce:** `python m5_8_2q_delta_scaling.py` (seed-level, ~10 s). The direct ω(δ) on the settled state (the `2h` run_dense path) is deferred. `m5_8_2q_omega.py` is the scaffold, but its fresh-seed probe reads a still-settling state (ω 0.60 ≠ canonical 1.10); the conclusion does not depend on it (it follows from the measured δ-flat H + the validated ZBW law).

### 2026-06-09 addendum: Duda's boost correction + the EM/GEM split

Duda disagreed with "rest energy dominated by gravity," and he is right. Verified against the seed code: the hedgehog frame `O4 = block-diag(O3, 1)` leaves the time axis fixed, and the boost `boost_field(b·w, a=1)` is the ONLY channel mixing the g/time axis into the spatial gradients. So `g` enters the energy only as the boost tilt `b·g`. The earlier `∝g⁸` came from raising `g` at fixed boost `b=0.13`, i.e. inflating the tilt `b·g` from ~1 to ~130, an unphysically large boost, not gravity.

**What we ran** (`2q` Phase D/E): split the quadratic seed energy into Duda's two sectors via the signed η-blocks (`SP_PAIRS` spatial = EM = curvature of rotations, η-positive; `TM_PAIRS` time-mixing = GEM = curvature of boosts, η-negative). `EM + GEM = H_quad` exactly (split validated).

| Finding | Result |
| --- | --- |
| boost = 0 | GEM (boost block) is **exactly 0**: gravity contributes nothing without the time-axis tilt |
| physical knob | `GEM ∝ (b·g)²` in the small-tilt regime (the `b=0.13,g=8` and `b=0.013,g=80` pair both give GEM ≈ −9.3; large pairs diverge via the `sinh` nonlinearity) |
| EM/GEM ratio (Duda's question) | NOT a constant: **210:1** at a physical small boost (`b=0.01`), **2:1** at the clock dressing (`b=0.13`); scales as `1/(b·g)²`. EM (the 1-axis, Faber unit vector) dominates the rest energy in every physical case |
| GEM sign | **negative** (the Minkowski clock-fuel block), so the boost/clock REDUCES the rest energy by `\|GEM\|` |
| mass-reduction (Duda's "how much?") | `\|GEM\|/EM ∝ (b·g)²`: ~0.5% at a physical boost, up to ~50% at the (large) clock dressing. Stopping the Zitterbewegung removes the negative GEM and the mass rises |

**Caveat:** this split is on the STATIC seed, so it weighs EM and gravity but NOT the quantum-phase δ sector, whose energy lives in the fast (~10²¹ Hz) clock evolution (Duda's point). That dynamical weighing is the open piece (NG-12).

**Reproduce:** `python m5_8_2q_delta_scaling.py` (Phase D/E, seed-level, ~15 s). Sent to Duda 2026-06-09.

## ELECTRON-ID PROJECT

The active program. Source: Duda's round-3 suggestion (`m5_4c_convo_2026.06.08.md §8-9`): "try getting such electron with 3x3 field and fixed clock, to get proper magnetic dipole moment and angular momentum of electron." The model has two of the four electron identifiers (mass = rest energy, charge = effective Coulomb from the hedgehog); this project targets the remaining two, **μ** (magnetic dipole moment, target one Bohr magneton, g ≈ 2) and **J** (spin, target ℏ/2).

Why the fixed-clock route: the `2p` readout on the dynamical 4D clock came out J-neutral (`J < 1e-4`, swamped by 24³ box torque), so μ/J wash out dynamically. A static 3×3 hedgehog with the clock pinned at definite phase/winding carries a definite circulating current (integrate μ) and a definite field angular momentum (integrate J), with the divergent boost sector dropped.

| Phase | What | Outcome (2026-06-10) |
| --- | --- | --- |
| EID-B: 3-way sector split | refine the `2q` Phase E split to Duda's exact `F_μν` figure: EM (tilt-tilt `R¹+g²R̃¹`) vs QM (tilt-twist `δR²−δ²R̃²`) vs GEM (boosts), separating the spatial block our 2-way split lumped | ✅ COMPLETE, all gates pass (sum exact to 1e-11; the 16.7379 gate holds). EM 16.34 / QM 2.23 / GEM −9.37 at clock dressing = Duda's hierarchy confirmed. KEY CORRECTION from the first run: tilt×tilt curvature points ALONG the major generator (`R=Γ×Γ`), so EM = component pair (1,2), QM = pairs (0,a), the reverse of the naive map |
| EID-C: the 3×3 fixed-clock electron | static hedgehog (3×3 field, no 4th axis), clock PINNED at definite phase/winding; integrate μ = ½∫r×J_curr dV from the circulating current and J from the field momentum density | ⚠️ COMPLETE WITH STRUCTURE FINDINGS: μ exists ONLY via the TILT (precession) channel (0.221 at 24³, linear response, b-independent); the TWIST clock (Duda's Γ¹) is EM-silent (abelian projection blind to twist, μ=0 structurally). ORBITAL J = 0 structurally (r×p, localized, Poynting all vanish: the hedgehog is dyon-like, E∥B kills Thomson); spin lives in the Noether clock charge L_int = Σ⟨P,Mth⟩ (the L/Q=ω family): 61.6, φ-flat to 0.03% |

Connects: NG-8 (the magnetic-dipole placeholder this makes real), NG-12(a) (EID-B's tilt-twist energy is the δ-sector's static weight), NG-1/NG-3 (the unit calibration that converts μ/J ratios into absolute statements). Report back to the Duda thread on results or for advice.

### EID implementation spec (frozen 2026-06-09, run-ready)

Both phases are SEED-LEVEL numpy (no evolution, no Taichi, no npz dependencies): the clock tangent comes analytically from `seed_M`, so the fixed-clock momentum needs no time stepping. New script: `scripts/m5_8_2r_electron_id.py`.

**Assets to reuse (verified on disk 2026-06-09):**

| Asset | Where | Role |
| --- | --- | --- |
| `seed_M(g, b) → (M, Mth)` | `scripts/m5_8_2c1_full_evolution.py` | the seed AND the analytic clock tangent `Mth = conj(W, G·D4 − D4·G)`, `G = gen4(PLANE)`; fixed clock ⇒ `Ṁ = ω·Mth`, no evolution |
| `rot4(plane, ψ)`, `gen4`, `conj`, `boost_field` | `scripts/m5_8_2a_4d_hamiltonian.py` | pin the clock at phase φ: `W(φ) = W·rot4(PLANE, φ)`; φ-sweep checks φ-independence of \|μ\|, \|J\| |
| `build_grid_n(n, box)` | `scripts/m5_8_2cb_taichi_constrained.py` | grids at 24³/32³/48³ for the box-robustness gate |
| seed constants `DELTA, RC, RHOC` / `L, B_STAR, R_W, PLANE, A_BOOST` | `scripts/m5_6_2a_biaxial_hedgehog.py` / `2c1` | the validated N-3 stack; `M58_DELTA`/`M58_G` env knobs exist |
| 2-way sector split `u_sectors` + gate `H_static = 16.74` | `scripts/m5_8_2q_delta_scaling.py` | EID-B extends this; the gate must keep passing |
| μ route A (abelian) | `scripts/m5_6_4a_hydro_em.py` | hydro↔EM dictionary: `B = ∇×A` from tilts, charge/current from the Lamb-vector divergence |
| μ route B (primary, Faber) | `scripts/m5_6_4b_faber_curvature_em.py` | `Γ_i = q0∂_iq − (∂_iq0)q + q×∂_iq`, `R_ij = Γ_i×Γ_j`, `*F ∝ R`; current from static Ampère `j = ∇×B` |

**EID-B recipe (the 3-way split):** the 2-way split partitions lab-frame index pairs, which cannot separate tilt from twist. Rotate `F` into the local eigenframe (exact on the seed: `F_eig = O4ᵀ F O4`, the analytic frame from `build_grid_n`'s `O4` + the boost). In the eigenframe with axes (0=major/1, 1=minor/δ, 2=zero, 3=time): TILT (EM) = pairs (0,1),(0,2) (they move the major axis); TWIST (QM) = pair (1,2) (rotation about the major axis, Duda's `Γ¹`); BOOST (GEM) = pairs (a,3). Gate: EM+QM+GEM = H_quad to machine precision, and at `b=0` the QM sector should carry the `δ²` weight Duda's figure predicts (`δR² − δ²R̃²`).

**EID-C recipe (μ and J at fixed clock):**

| Step | Formula | Note |
| --- | --- | --- |
| 1. Build the 3×3 electron | the 2c1 seed at `b_star = 0` (boost off), read the spatial 3×3 block; clock pinned at φ via `rot4(PLANE, φ)` | "3×3 + fixed clock" = Duda's prescription; gravity sector exactly absent (the `2q` Phase D result) |
| 2. Field momentum density | `p_i = −⟨Ṁ, ∂_iM⟩` with `Ṁ = ω·Mth(φ)` | ω is a free overall factor; J direction + the μ/J RATIO are ω-independent |
| 3. Spin | `J = Σ_vox r × p · h³` over the act mask | the 2p J-neutrality came from the DYNAMICAL kick; the pinned twist has definite winding |
| 4. Magnetic moment | route B primary: Faber `R = Γ×Γ` → `B`-field → `j = ∇×B` → `μ = ½ Σ r×j · h³`; route A as the abelian cross-check | the two routes agreeing is itself a result (M5.6.4 found both give Maxwell) |
| 5. The g-factor target | `g_factor = (μ/J)·(2m/q)` with m = H_static (lattice), q = the topological charge (winding = 1, M5.1) | ⚠️ CORRECTED 2026-06-10: NOT unit-free after all. μ comes out in director-curvature (EM-sector) units while the spin L_int is in action units; their relative normalization IS the Coulomb `e_scale` calibration. **g ≈ 2 becomes testable only after the NG-1/NG-3 unit fix**, exactly Duda's "fix units by comparing with Coulomb" |

**Gates + pitfalls (pre-registered):**

| Gate | Pass condition |
| --- | --- |
| G-EID-1 box robustness | J and μ stable across 24³ → 32³ → 48³ and under radial-mask sweep (the 2p box-torque lesson: if the box dominates, the number is the box's, not the electron's) |
| G-EID-2 φ-independence | \|μ\|, \|J\| independent of the pinned clock phase φ (sweep 4-8 phases); direction co-rotates with the winding |
| G-EID-3 factor-2 bookkeeping | apolar doubling: M-field probes read `ω_M = 2ω_clock` (the G7 rule, `m5_8_2j`); decide ONCE whether μ/J use ω_clock or ω_M and state it (the ratio μ/J is immune, the absolute values are not) |
| G-EID-4 split sanity | EID-B sectors sum to H_quad exactly; the existing `2q` gate (16.7379) keeps passing |
| Pitfall: eigenvector sign | use the analytic seed frame `O4`, never per-voxel `eigh` eigenvector signs (the NG-7 gauge lesson); numpy `eigh` is fine for VALUES |
| Pitfall: δ degeneracy at small δ | at δ→0 axes 1,2 degenerate and the eigenframe is ill-defined off-seed; stay on the analytic frame |

**Run plan:** (1) EID-B at the `2q` operating points (~minutes), (2) EID-C at 24³ with both μ routes + φ-sweep, (3) the box-robustness ladder 32³/48³, (4) the g-factor assembly, (5) document in this section + `4c` §9 + report to the thread. Total estimated compute: well under an hour, all CPU numpy.

### EID results (2026-06-10, `scripts/m5_8_2r_electron_id.py`)

Run-of-record numbers live in the script's RESULTS docstring; the outcome table above carries the headlines. The structural findings, one line each:

| Finding | Statement |
| --- | --- |
| The sector map | EM = the (1,2) component (tilt×tilt curvature along the major generator); the naive moved-axes map is backwards. With correct labels, Duda's hierarchy lands: EM dominant, QM δ-weighted small, GEM small negative |
| EM floor | EM(δ→0) → the `2q` δ-flat hedgehog floor (~19): the floor IS the EM curvature |
| μ channel | The dipole requires the TILT (precession) component of the Zitterbewegung; the pure twist phase (Γ¹) is EM-silent. μ(tilt, φ=0) = 0.221/0.248/0.277 at 24³/32³/48³ |
| Orbital J | Zero structurally: the centered hedgehog is dyon-like (emergent E ∥ B, no Thomson angular momentum) and the rigid clock is equivariant. The static face of the 2p J-neutrality |
| Spin | Lives in the Noether charge of the clock rotation, `L_int = Σ⟨P, Mth⟩` (the M6 L/Q=ω identity family): 61.6 (twist, 24³), φ-flat to 0.03% |
| Tilt-at-finite-φ | A finite tilt rotation destroys the hedgehog (φ=π/2 → disclination texture): the tilt channel is meaningful only as linear response |

**Residuals (fold into NG-1/NG-3):** (a) box convergence: μ and L_int both grow ~11%/step across the ladder (tail-dominated integrals; bigger boxes or radial windowing); (b) the g-factor needs the cross-sector `e_scale` normalization (the Coulomb unit fix) before g ≈ 2 is testable; (c) whether the spin-½ magnitude itself requires the Q₈ spinor structure (NG-9) rather than any classical integral remains open.
