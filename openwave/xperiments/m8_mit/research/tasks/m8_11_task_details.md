# M8.11: LOCAL BRANCH GERMS AT THE LEVEL-6 CRITICAL RAYS, and the order-a² tilt at the two rays the symmetry does not pin

> Roadmap row: [`../m8_roadmap.md`](../m8_roadmap.md) M8.11 (**author-proposed, maintainer-run**).
> Standing: [#512](https://github.com/openwave-labs/openwave/discussions/512#discussioncomment-18415036).
> Parents: [`m8_10_task_details.md`](m8_10_task_details.md) and [`m8_1_2_task_details.md`](m8_1_2_task_details.md).
> Status: DONE 2026-09-18: every frozen value reproduces blind; T1, S1 and C2's lemma pass as audited arguments. Maintainer-run, headless rooms.

## TASK PLANNING

### The question

M8.10 carried the fourth bedrock paper's range equation one order at its four symmetry-pinned critical rays (coherent `v_3`, zonal `v_0`, octahedron, hexagon) in sectors `3′` and `4`. It stayed within the paper's ceiling: the critical rays "are then the critical directions of the leading reduced problem, and not solution branches" (§ 6.1).

Its [method note](../findings/m8_10_method_note.md) (§ 8) leaves open finite-amplitude branches and the critical rays of Lemma 5.3(b) other than the pentagonal pyramid, there only a negative control.

The pyramid and the trigonal prism each lie on a projective fixed locus of complex dimension one, so their direction can move inside the locus at the next order. At the pyramid it must: the order-`a⁵` block equation has a component orthogonal to `Φ` (M8.10's D2).

This task asks whether the formal expansions are the Taylor expansions of actual branch germs at all six rays and, at the pyramid and the prism, what the order-`a²` correction to the direction (the tilt) is.

One local existence theorem answers the first. Finite harmonic algebra answers the second, with candidate exact values identified at high precision. The ceiling is local: germs for sufficiently small amplitude, with no radius, no stability and no claim at finite amplitude. T rests on the audit; the blind solvers verify its hypotheses (T2 to T4).

### Standing: the #512 ruling

The two conditions of #512 ([first](https://github.com/openwave-labs/openwave/discussions/512#discussioncomment-18222368), [second](https://github.com/openwave-labs/openwave/discussions/512#discussioncomment-18415036)) carry forward unchanged.

| Point | Ruling |
| --- | --- |
| Run-before-write | met on its plain reading: the last pre-registration, M8.10's, has run ([#547](https://github.com/openwave-labs/openwave/pull/547)) and closed ([#550](https://github.com/openwave-labs/openwave/pull/550)); the maintainer confirms at filing |
| Identity | the maintainer's to assign; M8.11 is the next free one, filed as a BACKLOG row |
| Budget | this task doc and its worklist together, within § 12.2's budget of about 8,000 words; nothing further lands until it has run |
| Defect claims | reproduced by the maintainer with independent code before ratification |
| Solver side | the analytic side: T is analysis on the full section space, and its hypotheses and every value come from finite algebra in the sectors' harmonics, with no truncation (the cubic stops at level 18), continuation, time-stepping, quadrature, C-room code or data, or solve at finite amplitude |
| The law | derived and pinned here, from nothing in `m8_5c/design_inputs/` ([#508](https://github.com/openwave-labs/openwave/issues/508)) |

### What this task is not

- **No radius.** T gives "there exists `ε > 0`" with no lower bound on `ε`, and none is sought.
- **No stability of any kind.** Deferred to their own decision: the reduced problem's second variation, the normal form's relative equilibria and a branch's spectral stability.
- **Not the rest:** the global critical set (§ 6.2's ceiling stands), the weight states `v_1` and `v_2` (covered by T, not expanded or scored), slot selection, OQ1, the stress tensor, the metric, a MODELS.md cell, or [M8.7](m8_7_task_details.md)'s gate, which it neither satisfies nor weakens.

**Short forms.** The BACKLOG row, the PR body and any comment say "local branch germs, for sufficiently small amplitude". Never write "branches" alone, and never call the result finite-amplitude.

### Ownership and run format

As in M8.10, the author supplies candidate values and an instrument, then stops. Recommended: two blind agents with separate implementations, then an adversarial audit that commits its own method first and grades T1 and S1 step by step, then the designer's comparison. Two blind agents and an audit earn blind; one maintainer reproduction earns [independent-method reproduction](../m8_roadmap.md#conventions).

### Sources of record

| Item | Pin |
| --- | --- |
| Paper | *The Surviving Ray*, Zenodo version DOI [10.5281/zenodo.22681502](https://doi.org/10.5281/zenodo.22681502), PDF MD5 `d2405316e20c060f53d338c1516298bc`, as pinned in M8.1.2. Used: § 2.5 (`σ_3′` and `σ_4` are real), Lemmas 4.1, 5.3 and 5.3(b), Corollary 5.6 and its remark, § 5.7, Proposition 3.3 and §§ 5.8, 6.1 and 6.2 |
| Parent verification | M8.1.2 ([#540](https://github.com/openwave-labs/openwave/pull/540)): its D3, D4, D7 and A5, reproduced blind. M8.10 ([#547](https://github.com/openwave-labs/openwave/pull/547), closed at [#550](https://github.com/openwave-labs/openwave/pull/550)): the pinned rays' values, now exact, and the pyramid values below |
| Engine | `openwave/xperiments/m4_ewt/wave_engine.py`, blob SHA-256 `6520ca41762cc4078660c00fd7bc4279094927f9182e09b0b373aa5dbe72a001`, last changed in `b910b5e40e145bc53ddb4c4d7a122f2f30876f09` (2026-07-27), read at `main` `59973a6499cd7c9c784e12d7bb1c29cbc821e8a3`: the potential table (lines 142-176, `v_mode 1` at line 150), the restoring force (line 339) and the continuous PDE (line 507) |

## THE LAW, DERIVED AND PINNED

As derived in [M8.10](m8_10_task_details.md#the-law-derived-and-pinned), read again at the `main` commit above, where the engine's bytes are unchanged:
- **(a) Engine.** `v_mode 1` (`cubic_nls`): `V = (c1/4)u²` and `dV = c1·u·ψ`, with `u = ψ·ψ`, evolved by `∂²ψ/∂t² = c²∇²ψ − dV(ψ)`; the sign of `c1` is the caller's.
- **(b) Modeling extension, called that.** The Hermitian `ψ†ψ` on sections of `E_σ` over `S³/2I`, so `V = (g/4)(ψ†ψ)²`, with gradients in the real Hilbert metric; the Wirtinger derivative is half of that.
- **(c) Reduction.** `ψ = e^{iωt}φ` gives the paper's § 2.5 equation with `c = 1`, `λ = ω²` and `g = c1`.

## SETTING AND DERIVATION

**Setting**, as in M8.10:
- the equation `(−Δ − λ)ψ + g|ψ|²ψ = 0` on sections of `E_σ`, with `σ ∈ {3′, 4}`, `g ≠ 0`, `R = 1` and `∫_X 1 = 1`;
- level `n` at `−Δ = n(n+2)`, the seven-dimensional level-6 block `K`, and `A = −Δ − 48`, inverted on `K^⊥`;
- the cubic `N(ψ) = |ψ|²ψ`, with real derivative `DN_Φ[h] = 2 Re(Φ†h)Φ + |Φ|²h`.

A fibre vector `u ∈ V_3` lifts to the block section `Φ = Φ_{σ,u}`, normalized as a section: `∫_X |Φ|² = 1`. The inner product `⟨·,·⟩` is antilinear in its first slot, and gradients and Hessians are taken in `⟨·,·⟩_ℝ = Re⟨·,·⟩`.

`Q_σ([u]) = ∫_X |Φ|⁴` is `λ₂/g` at a critical ray. By M8.1.2's D7, `Q_σ = 1 + w₆(σ) r̂₆`, with `w₆(3′) = 28/39`, `w₆(4) = 21/52` and `r̂₆ = ‖ρ₆(u)‖²/‖u‖⁴`.

**Symmetry.** Rotations act by left translation, commuting with `Δ` and `N`, and on block sections through the fibre vector: in the worklist's lift, `R·Φ_{σ,u} = Φ_{σ,D³(R̄)u}`, where `R̄` is the entrywise conjugate; `U(1)` acts by phase. Since `σ` is real (§ 2.5), fibre conjugation in a real basis of the sector is antiunitary and commutes with `Δ` and `N`; on the block's fibre vector it is time reversal up to a phase: `Θ(Σ u_m v_m) = Σ (−1)^m conj(u_m) v_{−m}`.

**The fixed space.**
- For a critical ray `[u]`, `H` is its stabilizer among rotations, acting by a character: `h·u = χ(h)u`.
- The twisted fixed space `X_H = {ψ : h·ψ = χ(h)ψ for all h ∈ H}` is closed in every Sobolev space. `K_H = K ∩ X_H` has complex dimension `m`.
- The kernel variable is `k = a·u(θ)`, with `r = a²`. Here `θ` is a coordinate on a local slice of `ℙ(K_H)` through `[u]`, transverse to the continuous symmetries that preserve `ℙ(K_H)`.
- `T` is the slice's tangent space at `[u]`, realized by unit sections `E ∈ K_H` with `⟨Φ, E⟩ = 0`. `Π_T` projects onto it real-orthogonally.

**Theorem T.**
- **T.a Function spaces.** `F(ψ, λ) = (−Δ − λ)ψ + g|ψ|²ψ` maps `H^{s+2}(E_σ) × ℝ` to `H^s(E_σ)`, for `s > 3/2`, between real Banach spaces. `H^s` is an algebra for `s > 3/2` in dimension 3, so the cubic is real-analytic. `A` is Fredholm of index zero, with kernel `K` and a bounded inverse on `K^⊥`. The nearest other level is 8 in sector `4` and 10 in sector `3′`, a gap of 32 or 72.
- **T.b The fixed space.** `F` maps `X_H × ℝ` into the fixed space of `H^s`, since `N(χψ) = χN(ψ)` for `|χ| = 1`. `N(Φ)` carries `χ` and `A⁻¹` is equivariant, so `ξ ∈ X_H`, and everything below lives in one space.
- **T.c Reduction.** Lyapunov–Schmidt in `X_H` solves the range equation for an analytic `w ∈ X_H ∩ K^⊥`. The reduced equations are the gradient of the reduced functional.

  The phase `−1 ∈ U(1)` sends `k` to `−k`, so the reduced map is odd in `a`. Its radial part divided by `a`, and its tangential part divided by `a³`, are therefore even analytic functions of `a`, hence analytic in `r = a²`:
  - **Radial:** `q(θ, r, λ) = (48 − λ) + g·Q(θ)·r + O(r²)`, with `∂_λ q = −1` at `r = 0`.
  - **Tangential,** divided also by `g`: `G(θ, r) = ¼∇_T Q(θ) + O(r)`, the `¼` coming from the bridge below.

  Dividing by `a` removes the trivial branch `ψ = 0`. Before that division, the radial equation's `λ`-derivative is `−a`, which vanishes at `a = 0`.

  At `(θ₀, 0, 48)`, where `θ₀` is the ray, the Jacobian in `(λ, θ)` is block-triangular, with diagonal blocks `−1` and `¼Hess_T Q(θ₀)`. If `Hess_T Q(θ₀)` is invertible on the slice, the implicit function theorem gives analytic `λ(r)` and `θ(r) = θ₀ + r·θ₂ + O(r²)`. The reduced functional is invariant, so its derivative along the symmetry orbits vanishes, which is why the slice suffices.
- **T.d Statement.** There are a neighborhood `W` of `(0, 48)` in `H^{s+2}(E_σ) × ℝ`, a neighborhood `U` of `θ₀` in the chosen slice and an `ε > 0` such that, for each `a` with `0 < |a| < ε`, exactly one solution `(ψ, λ) ∈ W`, with `ψ ∈ X_H` in the gauge of T.e, has `[Π_K ψ] ∈ U`. `λ` and `ψ/a` are analytic in `r = a²`, and rotations carry this germ to equivalent germs at the conjugate rays. No lower bound on `ε` is given.
- **T.e Identification.** The germ's Taylor coefficients satisfy the formal expansion's order-by-order equations. Those equations have a unique solution under the gauge: `⟨Φ, ψ⟩ = a` real, the range part orthogonal to `K`, and the tilt in the slice. So the formal coefficients are the germ's Taylor coefficients.

**E1: the four pinned rays (`m = 1`).** By Lemma 5.3, the stabilizer's isotypic space in the block is one-dimensional. So `K_H = ℂΦ`, which `U(1)` reduces to one real amplitude. There is no `θ`, nondegeneracy is vacuous, and T gives the germ. By T.e, M8.10's `λ₂`, `ξ` and `λ₄` are its Taylor coefficients, and its eight formal expansions become eight local branch germs.

**E2: the pyramid and the prism (`m = 2`).** Expand

`ψ = aΦ + a³(v + ξ) + O(a⁵)`,  `λ = 48 + λ₂a² + λ₄a⁴ + O(a⁶)`,

where the tilt `v` lies in `K_H`, with `⟨Φ, v⟩ = 0`, in the slice, and `ξ ⊥ K`.
- Order `a³` is as in M8.10: `λ₂ = gQ`, and `ξ = −gA⁻¹f` with `f = Π_⊥N(Φ)`.
- The block part of order `a⁵` is `g Π_K DN_Φ[v + ξ] − λ₂v − λ₄Φ = 0`. Its tangential part fixes the tilt:

`L_T v = −Π_T DN_Φ[ξ]`,  where `L_T = Π_T DN_Φ − Q` on `T`,

so `v ∝ g`.

**The bridge.** Take a unit tangent `E` in the block, real-orthogonal to `Φ`. Then `d²Q(cos s·Φ + sin s·E)/ds² = 4(⟨E, DN_Φ E⟩_ℝ − Q)` at `s = 0`. So `L_T = ¼Hess_T Q_σ = (w₆(σ)/4) Hess_T r̂₆`, which is the Jacobian block of T.c: one set of numbers decides both nondegeneracy and the tilt.

**The `λ₄` lemma.** Pair the order-`a⁵` block equation with `Φ`. The tilt contributes `g⟨Φ, DN_Φ[v]⟩ = g(2X + X̄)`, where criticality gives `X = ⟨N(Φ), v⟩ = ⟨Π_K N(Φ), v⟩ = Q⟨Φ, v⟩ = 0`. So

`λ₄ = g⟨Φ, DN_Φ[ξ]⟩ = −3g² Σ_n ‖Π_n N(Φ)‖²/(n(n+2) − 48)`,

which is M8.10's formula, as if the direction were held fixed. By M8.10's argument `λ₄ < 0`, since `Q ≠ 1` at both rays. Criticality is essential: at a non-critical point the pairing survives, which is the negative control N1.

**The pyramid.**
- **Locus and orbit.** The locus `ℙ span{v_2, v_−3}` is fixed by the fivefold rotation about `z`, under which both weights carry the same phase. `R_z(θ)` preserves the locus for every `θ` and moves along it. So the horizontal direction `i·e_t` is a rotation orbit, null for `Q`.
- **Slice and tangent.** The slice is the real curve `u(t) = cos t·v_2 + sin t·v_−3`, of unit speed and critical at `sin²t = 12/25` (M8.1.2's D3). `e_t` is the unit section tangent along increasing `t`: the lift of `−sin t·v_2 + cos t·v_−3`, normalized as a section.
- **Second variation.** M8.1.2's D3 quartic is `r̂₆ = −(125/132)s² + (10/11)s + 3/77`, with `s = sin²t`, stationary at the ray. So `Hess_T r̂₆(e_t, e_t) = r̂₆''(s)·sin²2t = −104/55`, in both sectors.
- **Relation to M8.10's record.** Its § 5.4 records the orthogonal component `R` along the unit fibre vector `e = sin t·v_2 − cos t·v_−3`, which points along decreasing `t`. A section's norm is `√(d/7)` times its fibre vector's, so `F = −√(d/7)·R` (A3).

**The prism.**
- **Locus.** The locus is the paper's `D₃` chart `u(z) = v_3 + z·v_0 + v_−3`, with `z = x + iy`. It is critical at `z = ±√(23/10)` (M8.1.2's D4), and the representative is `z₀ = +√(23/10)`.
- **No continuous symmetry along it.** `R_z(θ)` preserves the chart only when `e^{6iθ} = 1`. `R_z(π/3)` sends `z ↦ −z`, exchanging the two prism points.
- **Slice.** The slice is the whole chart, of real dimension two, and `K_H` is the lift of `span{v_3 + v_−3, v_0}`.
- **Tangents.** `τ_x` is the unit section tangent along increasing `x` at `z₀`: the horizontal lift of `∂_x u = v_0`, normalized as a section. `τ_y = i·τ_x` points along increasing `y`.
- **Metric.** The Fubini–Study metric there is `g_xx = g_yy = 2/(2 + |z|²)² = 200/1849`, so a unit tangent is `43/(10√2)` coordinate units.
- **Second variation.** `Hess_T r̂₆` in `(τ_x, τ_y)` is the Hessian of M8.1.2's D4 chart function divided by that metric: `diag(920/473, 8/11)`, in both sectors.

**The S lemma (prism).** `Θ` acts on the chart as `z ↦ −z̄`, and `R_z(π/3)` as `z ↦ −z`. So `S = R_z(π/3)∘Θ` acts as `z ↦ z̄`, and it fixes the prism points.

`S` lifts to sections as a rotation composed with fibre conjugation. So it maps solutions to solutions, because `−Δ`, `λ` and `g` are real and the conjugate of `|ψ|²ψ` is `|ψ|²ψ̄`. Being antiunitary, it fixes the prism's `Φ` once the phase of `Φ` is chosen. Then:
- by uniqueness, `Sξ = ξ`, so the forcing is `S`-even;
- `τ_x` is `S`-even and `τ_y` is `S`-odd, so the forcing has no `τ_y` component, and `L_T`, an `S`-even form, has no cross term.

Hence `v_y = 0` exactly. `v_x` is computed, and it is nonzero in both sectors.

**Levels,** as in M8.10 and keyed by `n`, not `j = n/2`. None vanishes at the two new rays (B1): level 16, carried by the spin-8 channel, vanishes only on time-reversal-invariant rays (M8.1.2's A5), and neither ray is one.

## THE FIREWALL

M8.10's firewall ([its requirements](m8_10_task_details.md#the-firewall)) applies, and the contexts also have no access to the records of M8.1.2 and M8.10 (task docs, method notes, packets, returns and packages), since M8.10's record holds the pyramid values graded again here.

**Answer key and withheld terms.** This file is the answer key and the leak path of record; at go the designer states where it is held and how solvers are kept from it.
- **Withheld:** the author, the model, the repositories, the paper's title, "surviving ray", every frozen value, every claim ID and M8.10's ray labels.
- **Permitted:** the vocabulary the problem needs, including the six rays as explicit fibre vectors with their critical parameters, the slice tangents with their orientations, and `Θ`.

**Designer-only material: the author's package**: the derivation's code, outputs and clean-directory logs, the 100-digit repeat, the failed first pass's log, and `NOTES.md`, whose post-derivation observations are not claimed. No agent sees it; the designer opens it after adjudicating, for provenance comparison. Every file passed a privacy scan before hashing, and `m810_core.py` and `m810_exact.py` are byte-identical to #550's copies. It lands after the verdict under `research/scripts/m8_11_author/`, byte-identical to these hashes:

| file | SHA-256 |
| --- | --- |
| `m810_core.py` | `9339d87117a4fe524b913f5fe63972f9a19dfaa0d708f059323bde7553f7ecda` |
| `m810_exact.py` | `8436ef9ffc465e6a91878b4b994c1385766167cc8ebc4047d482cf65c72463a9` |
| `m811_ops.py` | `6753b9827e69624180d3eb8d3b627fc242fd41b4a81a3c1c8d5e6ec36bdae1c0` |
| `m811_pyramid.py` | `d72c58e18763e67488c531540deedccaca2d139f0993e435af9c409a29582eff` |
| `m811_prism.py` | `58e60e8bbe2f1621ae2c3f765d9bb6f5f934cbec5703f5967537e8486c77dd6c` |
| `m811_exact.py` | `5ea22c908f6501aa03fea8c42ec328d8b7ec4fbc1207e6143840a2a631b01bd5` |
| `NOTES.md` | `3687e056f8c0a10f7c0b4c11ac14b25df75eb1ad2706c251ea1e33f7737cef47` |
| `out/pyramid.json` | `cd9a8c5204af9cbbca439b590e3898248562205fea779ae56e3ca37b300299b0` |
| `out/prism.json` | `018e8274490ab6b8c69d5833ebbe0766c820759561ce716bb349936dfa124b01` |
| `out/exact.json` | `16e1ffc83b20bb8fba5c47dc6d5815357d2c1da45a158af681e31e96588008c4` |
| `out/step1_pyramid_clean.log` | `7cb5eabd9939c9cc90ad436e655ac35b57d8fe0b8b6e33a5abd49af510f39368` |
| `out/step2_prism_clean.log` | `9cb6875655c0bfd2674d53c94aab23067d95f69ee6805a590fba6117a1f86f73` |
| `out/step3_exact_clean.log` | `980b34f179f014e3b02118fb9c920f5e657048060be433a581ee8414f88a125b` |
| `out/step3_exact_dps100.log` | `63e7bf37618414e1de89808af416fa79461f394b341c3cfbb8a664c1392f73ee` |
| `out/step3_exact_firstpass_sqrt230_failed.log` | `6e47dfb50115a80f4e4026266cb71d2000f68be7367a70ef6aa557ed25af8074` |

## DISCLOSURE

**Exposure before the derivation.** For E1, the prior numbers are M8.10's, now exact. For E2, every prior number was published or verified:
- the paper's formulas, as M8.1.2's D3 and D4 record them;
- M8.10's orthogonal component, and #547's pyramid values (`Q`, the level norms and `λ₄`), which the author had read;
- the in-locus second variations, computed exactly from those formulas in an author-side feasibility check that both review units verified. Their signs were `(−, −)` at the hexagon, `(−, +)` at the octahedron and `(+, +)` at the prism;
- a review unit's numerical check that `S` fixes the prism ray.

The prism's forcing, tilt, level norms and `λ₄` had not been computed before the derivation.

**Exploratory work (2026-08-31 to 2026-09-05), for the deferred stability question:** the numerical Hessian of `Q` at the hexagon (nine negative, five null, none positive; sector ratio `16/9`), exact spectra of the block multiplication operator at six critical orbits, and the hexagon branch by Newton continuation to `a ≈ 8.7` (code as in M8.10's disclosure). None enters a claim, and the continuation is not evidence for T. Label trap: those notes call the sectors `R4 = 3′` and `R5 = 4`, M8.10 calls the rays `R1` to `R5`; this document names rays by shape.

**The derivation**, author-side (the author and the author's AI agents), used a local folder importing numpy, scipy, mpmath and two pinned M8.10 files (the group, `D³`, sector projectors, exact Racah coefficients), and nothing from the rooms, design inputs, M8.10's agents or exploratory scripts. Its three steps: the bridge and the pyramid; the prism, S lemma first; a separate high-precision route sharing only those two files with the float pipeline and reading its outputs once, at the end. Fresh, not blind: two author-side review units redlined the outline (three revisions) and each step, which is review, not verification.

**Caught in the derivation, and kept.** Three catches share one failure mode, a basis vector's normalization: step 1's first comparison with M8.10's orthogonal component failed by exactly `√(d/7)`, section against fibre (now claim A3); the prism field stated before step 3, a rational times `√230` (a review unit's suggestion, adopted), failed in both sectors because the unit tangent's `√2` was dropped, so the field is `ℚ·√115`, with the failed log kept and `√230` an arm that must fail; and the `w₆/4` and Fubini–Study conversions ride on two-route gates whose arms fail without them. Hence every value row below names its basis.

## CANDIDATE PRE-REGISTERED CLAIMS

Frozen by the designer at go. Values are stated per power of `g`: `λ₂/g`, `F/g`, the tilt over `g`, `‖v‖²/g²`, `‖Π_n ξ‖²/g²` and `λ₄/g²`.

**Status: candidate exact.** Except for `Q` and the second variations, which are derived from the pinned parent formulas, and `‖v‖²`, which follows exactly from the tilt, each value frozen as candidate exact was computed by finite algebra at 80 digits and identified as a rational, or a rational times the field element predicted in advance (`√39` at the pyramid, `√115` at the prism), accepted only with denominator at most `10²⁵` and residual below `10⁻⁶⁸`. Distinct rationals with such denominators differ by at least `10⁻⁵⁰`, so each identification is unique in that class, and the route returns byte-identical values at 100 digits.

The largest identified denominator is 19,152,322,028,017,152.

This is identification, not symbolic derivation. A blind reproduction that satisfies [#547's exactness rule](m8_10_task_details.md#the-go-time-checklist-answered) promotes a value at adjudication, labeled by its route. The bridge, T, the `λ₄` lemma, the S lemma and `v_y = 0` are arguments, and are graded as arguments.

**New or reproduced**, row by row: at the pyramid, `Q`, the level norms and `λ₄` (reported exactly by #547's two blind agents, ungraded) and the orthogonal component (graded as M8.10's D2) are reproductions, while its second variation, tilt, existence and `λ₄` lemma are new; at the prism, `Q` follows from M8.1.2's D4 and D7 and everything else is new; at the four pinned rays, only the existence statement and its hypotheses are new.

The items of #547 cited below are those of its [as-run worklist](../m8_10/worklist_as_run.md).

### Group A: parents, reproduced rather than new

| ID | Claim | Standing | Pass condition | Fail condition |
| --- | --- | --- | --- | --- |
| A1 | the critical data of M8.1.2's D3 and D4: on the pyramid's curve, stationary at `s = 12/25` with value `9/35`; on the prism's chart, critical at `z = ±√(23/10)` with value `200/903`; both points critical in the whole block | M8.1.2 D3 and D4; Lemma 5.3(b) | both values reproduced, and whole-block criticality shown by the stabilizer argument or by exact evaluation | a value differs, or criticality is shown only on the locus |
| A2 | `λ₂/g = Q = 1 + w₆(σ) r̂₆` at both rays in both sectors, as tabulated | M8.1.2 D7 with A1; at the pyramid, also #547's `R5` item 3, ungraded | the four fractions | any differs |
| A3 | the pyramid's forcing is M8.10's graded orthogonal component: `F = −√(d/7)·R`, with `R` the coefficient along the unit fibre vector `e = sin t·v_2 − cos t·v_−3` recorded in M8.10's method note (§ 5.4): `−7188839√91/81061695000` in `3′` and `+19565553√273/384292480000` in `4` | M8.10's D2 (#547), graded; exact values from both agents | the identity `F = −√(d/7)·R` holds between D1's pyramid forcing and the recorded `R`, in magnitude and sign | any mismatch, including one by `√(d/7)` or by a sign |

| ray | normalization | `λ₂/g`, sector `3′` | `λ₂/g`, sector `4` |
| --- | --- | --- | --- |
| pyramid, `sin²t = 12/25` | unit section `Φ` | `77/65` | `287/260` |
| prism, `z₀ = +√(23/10)` | unit section `Φ` | `5831/5031` | `609/559` |

### Group T: the existence theorem

| ID | Claim | Standing | Pass condition | Fail condition |
| --- | --- | --- | --- | --- |
| T1 | Theorem T, T.a to T.e, as stated in the setting: the author's proof, graded by the audit | new; a first: no filed M8 task grades an author's proof | the auditor's own existence argument, committed at stage 1 before it sees T.a to T.e, establishes the result, saying where the trivial branch is divided out and what fails without that step; at stage 2 the auditor grades each step of T.a to T.e, validating the stated route or supplying an equivalent derivation. The solvers' item-10 arguments are a diagnostic, scored on seven elements (function spaces, Fredholm structure, fixed space, reduction, nondegeneracy, trivial branch, identification; five including the trivial branch make it complete) | the auditor's argument does not establish the result, or applies the implicit function theorem before removing the trivial branch; a false statement in the frozen theorem text, or a necessary hypothesis it leaves unstated, even if the theorem can otherwise be rescued; or a step the audit cannot establish |
| T2 | `dim_ℂ K_H = 1` at the four pinned rays and `2` at the pyramid and the prism, in both sectors | new | computed from each stabilizer and its character, each group of rays serving as the other's control | any other dimension |
| T3 | `ξ ∈ X_H` at all six rays in both sectors, and `Sξ = ξ` at the prism | new | shown for the computed `ξ`, by equivariance or by exact evaluation | any failure |
| T4 | nondegeneracy on the slice, as tabulated: `L_T = (w₆/4) Hess_T r̂₆`, the Hessians of `r̂₆` the same in both sectors, the orbit direction `i·e_t` null and no cross term at the prism | new; from the formulas of M8.1.2's D3 and D4, which state no second variation | the tabulated values, the null orbit direction, and no cross term at the prism | a value differs, a cross term appears, or an eigenvalue on the slice is zero |

| ray | direction (unit section tangent) | `Hess_T r̂₆`, both sectors | `L_T`, sector `3′` | `L_T`, sector `4` |
| --- | --- | --- | --- | --- |
| pyramid | `e_t`, along increasing `t` | `-104/55` | `-56/165` | `-21/110` |
| pyramid | `i·e_t`, the rotation orbit | `0` | `0` | `0` |
| prism | `τ_x`, along increasing `x` at `z₀ = +√(23/10)` | `920/473` | `6440/18447` | `2415/12298` |
| prism | `τ_y = i·τ_x`, along increasing `y` | `8/11` | `56/429` | `21/286` |
| prism | the cross term `(τ_x, τ_y)` | `0` | `0` | `0` |

### Group S: the prism's S lemma

| ID | Claim | Standing | Pass condition | Fail condition |
| --- | --- | --- | --- | --- |
| S1 | `S = R_z(π/3)∘Θ` acts on the chart as `z ↦ z̄` and fixes the prism points. It maps solutions to solutions, and with the phase of `Φ` chosen so that `SΦ = Φ`, `ξ` and the forcing are `S`-even, `τ_x` is `S`-even and `τ_y` is `S`-odd, and `v_y = 0` exactly | new | the argument, with the section-level compatibility derived rather than assumed, and each consequence checked on the computed objects | a false statement in the claim, even if `v_y = 0` holds another way; an unsupported step; or a nonzero `S`-odd component |

### Group B: the level norms at the two new rays

| ID | Claim | Standing | Pass condition | Fail condition |
| --- | --- | --- | --- | --- |
| B1 | the 18 per-level norms `‖Π_n ξ‖²/g²` at the pyramid and the prism, as tabulated; none is zero | pyramid: #547's `R5` item 4, ungraded; prism: new | every entry matches its frozen candidate value under the go-time exactness rule, and none is zero | an entry differs, or an entry is zero |

Sector `3′`, `‖Π_n ξ‖²/g²`:

| ray | normalization | `n = 10` | `n = 14` | `n = 16` | `n = 18` |
| --- | --- | --- | --- | --- | --- |
| pyramid, `sin²t = 12/25` | unit section `Φ` | `77/5475600` | `281211/189112352000` | `553/1591200000` | `30233/56458242360` |
| prism, `z₀ = +√(23/10)` | unit section `Φ` | `319550/88158077163` | `12920075/7517877885232` | `16583/68316230736` | `2355742375/6059914391677302` |

Sector `4`, `‖Π_n ξ‖²/g²`:

| ray | normalization | `n = 8` | `n = 12` | `n = 14` | `n = 16` | `n = 18` |
| --- | --- | --- | --- | --- | --- | --- |
| pyramid, `sin²t = 12/25` | unit section `Φ` | `3591/346112000` | `63/83200000` | `361557/931014656000` | `34839/147097600000` | `2371131/8029616691200` |
| prism, `z₀ = +√(23/10)` | unit section `Φ` | `7875/859947712` | `3521/1587595776` | `16611525/37011091127296` | `116081/701717332992` | `4105722425/19152322028017152` |

### Group C: `λ₄` at the two new rays

| ID | Claim | Standing | Pass condition | Fail condition |
| --- | --- | --- | --- | --- |
| C1 | `λ₄/g²` at the pyramid and the prism in both sectors, as tabulated | pyramid: #547's `R5` item 7, ungraded; prism: new | all four match their frozen candidate values under the go-time exactness rule | any differs |
| C2 | the `λ₄` lemma: at a critical ray the tilt does not change `λ₄` at this order, so `λ₄ = −3g² Σ_n ‖Π_n N(Φ)‖²/(n(n+2) − 48)`, and `λ₄ < 0` for `g ≠ 0` | new; the sign argument is M8.10's C2 | the lemma derived from criticality, `λ₄` agreeing with and without the tilt, and the sign argued rather than read off | the lemma asserted without criticality, a disagreement, or the sign read from the totals |

| ray | normalization | `λ₄/g²`, sector `3′` | decimal | `λ₄/g²`, sector `4` | decimal |
| --- | --- | --- | --- | --- | --- |
| pyramid, `sin²t = 12/25` | unit section `Φ` | `-267786421/58544557500` | `-4.574062e-03` | `-4797453339/2497901120000` | `-1.920594e-03` |
| prism, `z₀ = +√(23/10)` | unit section `Φ` | `-336158940460/150812349114141` | `-2.228988e-03` | `-11093213145/4965015608696` | `-2.234276e-03` |

The sector ratios follow from C1 and are not a separate claim: `λ₄(3′)/λ₄(4)` is `699523712/293721633` (`2.381587`) at the pyramid, #547's `R5` item 8, ungraded, and `10976618464/11002656303` (`0.997633`) at the prism.

### Group D: the order-a² tilt

| ID | Claim | Standing | Pass condition | Fail condition |
| --- | --- | --- | --- | --- |
| D1 | the forcing `F = ⟨E, DN_Φ[ξ]⟩_ℝ` over `g`, along `e_t` at the pyramid and along `τ_x` at the prism; zero along `i·e_t` and along `τ_y` | pyramid: reproduced (A3); prism: new | matches the frozen candidate value in the stated basis, under the go-time exactness rule | a value differs, or a component along `i·e_t` or `τ_y` is nonzero |
| D2 | the tilt over `g`: `β` along `e_t` at the pyramid; `v_x` along `τ_x`, and `v_y = 0`, at the prism | new | matches the frozen candidate values, with the signs in the stated orientations, under the go-time exactness rule | any differs |
| D3 | `‖v‖²/g²`, which does not depend on the basis: the primary tilt claim | new | matches the frozen candidate value under the go-time exactness rule | any differs |
| D4 | with the tilt, the tangential part of the order-`a⁵` block equation vanishes in every slice direction; without it, that part equals D1 | new | shown exactly or to a stated precision | a nonzero residual with the tilt |

| ray | sector | basis: unit section tangents | `F/g` | tilt over `g` | `‖v‖²/g²` |
| --- | --- | --- | --- | --- | --- |
| pyramid | `3′` | `e_t`, increasing `t` | `7188839*sqrt(39)/81061695000` (`5.538286e-04`); `0` along `i·e_t` | `β`: `1026977*sqrt(39)/3930264000` (`1.631816e-03`) | `1054681758529/396076284864000000` (`2.662825e-06`) |
| pyramid | `4` | `e_t`, increasing `t` | `-19565553*sqrt(39)/192146240000` (`-6.359054e-04`); `0` along `i·e_t` | `β`: `-931693*sqrt(39)/1746784000` (`-3.330933e-03`) | `2604155538747/234711872512000000` (`1.109512e-05`) |
| prism | `3′` | `τ_x`, increasing `x` at `z₀`; `τ_y = i·τ_x` | `-120576113*sqrt(115)/1288994436873` along `τ_x` (`-1.003134e-03`); `0` along `τ_y` | `v_x`: `17225159*sqrt(115)/64285514280` (`2.873420e-03`); `v_y`: `0` | `296706102575281/35935889967339860160` (`8.256540e-06`) |
| prism | `4` | `τ_x`, increasing `x` at `z₀`; `τ_y = i·τ_x` | `214430223*sqrt(115)/3055394220736` along `τ_x` (`7.526060e-04`); `0` along `τ_y` | `v_x`: `-10210963*sqrt(115)/28571339680` (`-3.832525e-03`); `v_y`: `0` | `104263765387369/7098447400956021760` (`1.468825e-05`) |

### Negative control

| ID | Claim | Standing | Pass condition | Fail condition |
| --- | --- | --- | --- | --- |
| N1 | at `sin²t = 1/4` on the pyramid's curve, the point is not critical (`dr̂₆/dt = 115√3/528`), and the lemma's pairing `Re⟨Φ, DN_Φ[e_t]⟩` is nonzero there, so criticality is load-bearing. Its value, `3·Re⟨N(Φ), e_t⟩ = (3/4)·w₆(σ)·115√3/528`, is recorded for reference: `805*sqrt(3)/6864` in `3′` and `2415*sqrt(3)/36608` in `4` | new; a derivative of M8.1.2's D3 quartic | nonzero in both sectors | zero in either. A nonzero value that misses the recorded one is logged under X1, not as an N1 failure |

### Diagnostics (not claims)

| ID | What it detects |
| --- | --- |
| X0 (packet check, before solving) | the two quaternions are units and generate a perfect group of order 120 |
| X1 (normalization) | the three places a basis vector's length enters. **The state:** a unit fibre vector in place of a unit section scales `λ₂` by `d/7`, and everything after it. **A block section read through its fibre vector:** its coefficient along a unit fibre direction is `√(7/d)` times its coefficient along the unit section direction, which is A3's factor. **The prism's tangent:** a coordinate tangent `∂_x` in place of `τ_x` multiplies the forcing by `10√2/43`, which moves its field to `ℚ·√230`; it multiplies the tilt by `43/(10√2)` and the second variation by `200/1849`. A nonzero N1 pairing off its recorded value points to the factor 3, the `¼` bridge or `w₆` |
| X2 (gradient convention) | taking the Wirtinger derivative as the force halves `g`, which halves `λ₂` and quarters `λ₄`; dropping the conjugate term of `DN` puts `1/8` in the bridge in place of `1/4` |
| X3 (the level-18 plane) | a one-copy Schur scalar reused at level 18 in sector `4` misses part of that level |
| X4 (the cubic, pointwise) | the order-`a³` range equation checked with an explicit Casimir, and `N(Φ)` also evaluated pointwise, since a coefficient-space residual cannot see an error in `N` itself (#547's audit) |
| X5 (orientation) | at `z = −√(23/10)`, or with `τ_x` along decreasing `x`, `v_x` changes sign; with M8.10's `e = −e_t`, `β` and `F` change sign; `‖v‖²` changes in neither case |

### Feasibility

T1, S1 and C2's lemma are arguments; the rest are finite evaluations. A partial audit verdict on one step of T1 is a legitimate outcome, not a packet defect, and so is an S1 failure, since the worklist leaves the symmetry for blind agents to find.

### Instrument qualification (2026-09-13)

As M8.10's was, the worklist was dry-run once by a fresh headless Claude session (Opus 4.6), the author's AI agent, in a room outside every repository holding only the worklist and a brief, the answer locations and the network denied to every tool. A canary loaded nothing answer-bearing but did load the account email and username. The transcript audit, its planted calls firing, found no call outside the room and no refused attempt in 187 tool calls. This is not a verification: the agent and the instrument are the author's.

Against criteria frozen before the run, item 10 missed where the trivial branch is divided out, so T1 is graded as in its row; item 9 found an antiunitary symmetry unprompted but assumed that it preserves the sector. Of 104 value comparisons, 65 mismatched and 6 were missing; all but three, the agent's own slips, trace to one convention the worklist left implicit: each level's intertwiner was built for its own realization of `σ`, which scaled every level norm by a fixed factor and dropped the level of multiplicity two. Stabilizers were also taken inside `Γ`. The agent explained both symptoms away: the dropped level as a symmetry zero, item 12's misfit as truncation.

The worklist now defines rotations and fixes one `σ` at every level, both copies kept, item 1 reporting the equivariance residuals; item 12 checks the range equation as a vector identity, `N(Φ)` projected pointwise through level 18, a residual above the stated precision unresolved; item 9 keeps missing values apart from zeros and requires a symmetry to preserve the equation and the sector; item 10 asks what fails without each step. No expected value changed. Instead of a second run, each new check was run on the run's own code: it fails there and on planted errors, and passes after the repair, which recovers all 54 frozen level norms. In drafting, four false statements in T were caught, three by the author and one by a review unit. The run used the worklist with SHA-256 `1dd5842e11e66efe04a4e32782dc99be059dc23a2ec3ea8b8823dd51f78c8d7f`; its brief and return join the author's package, its transcript and checks only after deterministic privacy treatment.

### To be fixed at go (before numerics)

M8.10's [go-time items](m8_10_task_details.md#to-be-fixed-at-go-before-numerics) apply, plus:

| Item | Note |
| --- | --- |
| Exactness rule | #547's unless changed: symbolic derivation, or identification stating precision, denominator bound and field, repeated at a second precision |
| Argument grading | what counts as a gap in T1 and S1 |

### Definition of done

[M8.10's seven items](m8_10_task_details.md#definition-of-done-skeleton-finalized-at-go), finalized at go, with `m8_11_` prefixes and the method note at `findings/m8_11_method_note.md`.

### What is frozen

From the designer's freeze, the values stop moving on the author's side: no edit because a computed result disagrees with a frozen value. A disagreement is a result, and a claim that this text is defective goes to maintainer reproduction under #512. This binds the author, not the maintainers.

## GO-TIME PRE-REGISTRATION (2026-09-18, go 11:50 EDT)

Written by the designer and frozen before any room received a packet. Nothing in this section is edited after the rooms launch; anything the run forces off it goes in the deviations log below.

### The go-time checklist, answered

| Item | Decision |
| --- | --- |
| Deposit pin | unchanged from M8.10's go ([#546](https://github.com/openwave-labs/openwave/pull/546) review): MD5 `d2405316e20c060f53d338c1516298bc` |
| Engine pin | `wave_engine.py` blob SHA-256 `6520ca41762cc4078660c00fd7bc4279094927f9182e09b0b373aa5dbe72a001` at `main` `3e116e3d`, unchanged |
| Claims frozen | the claims tables above, as written, with no value edited. A designer transcription (`frozen_claims.json`, SHA-256 `272a79100acaad739011eba6eae15fc9372891951ee0095ee99f55a3b8d6117f`) was checked against this file and against the tables' own identities before launch: `λ₂/g = 1 + w₆ r̂₆`, `L_T = (w₆/4) Hess_T r̂₆`, A3's `F = −√(d/7)·R`, tilt `= −F/L_T`, `‖v‖² =` tilt², `λ₄/g² = −3 Σ (n(n+2) − 48) ‖Π_n ξ‖²/g²` at all four new expansions, both sector ratios and N1's `(3/4) w₆ · 115√3/528`, all exact |
| Instrument | the offered worklist, adopted unchanged (SHA-256 `d0895cf4...`, the file reviewed at [#554](https://github.com/openwave-labs/openwave/pull/554)). It already asks items 4 to 8 at all six fibre vectors, so no ray is marked |
| Handout audit | maintainer-side: the semantic read at the #554 review, and a go-time gate over the 53 fractions, 18 decimals and 123 integers of three or more digits in the claims section, plus 56 withheld terms (claim IDs, ray names and labels, author, model, repositories, paper). The only matches are `1/4`, `12/25` and `23/10`, which define the fibre vectors and N1's point. A planted value and term fire the gate |
| Run format | two blind rooms, then the auditor's second stage over the solver's work and the theorem text, per the roles table below. Earns blind for the values; the arguments earn the label below |
| Exactness rule | #547's, unchanged: exact means a symbolic derivation, or an identification that states its precision, denominator bound and field, repeated at a second precision. Equality with a frozen value is decided symbolically, so any equivalent radical form matches. The record labels each value by the route that produced it |
| Argument grading | T1, S1 and C2's lemma are graded as arguments, by the auditor at stage 2, one verdict per step: **ESTABLISHED** (the stated route holds, or the auditor supplies an equivalent derivation, shown), **GAP** (the step cannot be established) or **DEFECT** (a false statement, or a necessary hypothesis left unstated, even if the result can be rescued). T1 passes only if the auditor's own stage-1 argument establishes the result, locating the trivial branch's division and what fails without it, and every step of T.a to T.e is ESTABLISHED. A GAP gives a partial verdict and a DEFECT fails the claim, per its row. The designer reads the grading and may overrule a verdict only with a stated reason in the method note. The solvers' item-10 arguments are the diagnostic of T1's row, scored on its seven elements |
| Label | a passing argument is recorded as an **audited argument**: checked step by step by one AI auditor and read by the designer. It is never recorded as a verified or proven theorem, since AI agents sit on both the drafting and the checking side ([`AI_HYGIENE.md`](../../../../../AI_HYGIENE.md)) |
| Compute | each room's interpreter pins the math libraries to one thread and runs at `nice 15`: Python 3.12.14, numpy 2.5.3, scipy 1.18.1, sympy 1.14.0, mpmath 1.3.0 |
| Answer-key containment | the table below |
| Author's package | not in the repository. The provenance comparison waits for its landing PR, after the verdict |

| Packet file | Room | SHA-256 | Bytes |
| --- | --- | --- | --- |
| `worklist.md` | both | `d0895cf4327ef162df58b62a322fceac3a3f2095a78992bc3e88443cc4e35c17` | 10,638 |
| `BRIEF.md` (solver) | solver | `f5fa720e6453184b09c1b7cd750c03c4126faf6d807db9285c555b42451e2216` | 1,908 |
| `BRIEF.md` (auditor, stage 1) | auditor | `7b16836de313aad826009a6a2719c2eff3a14f8c42b0807d1c08ede8e29e8b51` | 2,432 |
| `BRIEF_stage2.md` | auditor, stage 2 | `e4822e54a180cf150d96e62511ec87554705a0569b0830f50115fc4ab4c9d1f7` | 2,347 |
| `THEOREM.md` | auditor, stage 2 | `2b25f48ace681cff392eb547f3913f7eb0b589d64da927bf7c695a7df90e9e74` | 7,283 |

`THEOREM.md` transcribes the setting's symmetry and fixed-space paragraphs, Theorem T (T.a to T.e), E1, E2's expansion and loci, the bridge, the `λ₄` lemma and the S lemma, with references to the paper and to parent tasks removed and no frozen value except the gap `32` or `72` and the prism parameter. The packet files land beside the returns at FINISH, under `m8_11/`.

### Containment: where the answer key lives, and what walls it off

The frozen claims live in this file, in a maintainer worktree outside the rooms. Both rooms run headless, per [`CLEAN_ROOM_STANDARDS.md` § 3.4](../../../../../dev_docs/CLEAN_ROOM_STANDARDS.md), launched by [`clean_room_launch.sh`](../../../../../dev_docs/utils/clean_room_launch.sh) in neutral folders under `/tmp/cr-20260918-m811/`. No in-session subagent is used at any step.

| Route to the answers | Guard |
| --- | --- |
| This file, M8.1.2's and M8.10's records, every repository | the session's file tools are restricted to the room; the shell runs only `./py`, whose OS sandbox refuses every read under `/Users` and `/private/tmp` except the room |
| Web search and fetch, connectors, skills, spawning agents | the session has four tools (read, write, edit, shell), no MCP server and no slash command |
| Instruction files and memory | `--restricted` loads none; the launcher aborts on any instruction file on the room's ancestor path |
| The Python route | `import openwave` must fail from the room, and the sandbox must refuse an outside read that succeeds unsandboxed, or the launch aborts |
| Network | refused to the interpreter by the sandbox |
| What still reaches a room | the account email and the room's own path, disclosed in § 3.4 |
| After the run | every tool call in each stream-json transcript is classified against the room; a call outside it is a protocol failure |

### Roles and ordering

| Step | Who | Receives | Before the next step |
| --- | --- | --- | --- |
| 1 | solver and auditor, in parallel, separate rooms | the worklist and a brief. The auditor writes its method, including its item-10 outline, before computing, on a route with no library coupling tables | each return snapshotted and hashed outside the room |
| 2 | auditor, a second headless launch in the same room | its own saved stage-1 files, the solver's work and `THEOREM.md` | per-item verdicts on the solver, a hunt for checks that cannot fail, and the per-step grades of `THEOREM.md` |
| 3 | designer | everything | comparison against the frozen claims, X0 as a precondition and X1 to X5 recorded apart from the verdicts; the argument grades read and recorded |
| 4 | designer, only after the verdict is recorded | the author's package, when it lands | provenance comparison, with the agreement and disagreement asymmetry stated |

### Definition of done, finalized

| # | Item |
| --- | --- |
| 1 | Solver and auditor returns, scripts and data in the repository, under `scripts/m8_11_solver/` and `scripts/m8_11_audit/` |
| 2 | Adversarial audit with its own method, per-item verdicts, a hunt for checks that cannot fail, and the per-step argument grades |
| 3 | Designer comparison against the frozen claims, every number stated; diagnostics recorded apart from the verdicts |
| 4 | Transcript audit of every room session |
| 5 | Method note `findings/m8_11_method_note.md`: equations first, equation-to-code map, audit record, manifests |
| 6 | Author package landed byte-identical to its hashes and compared for provenance only after the verdict, with the asymmetry stated (its own PR) |
| 7 | Doc sync (roadmap row and briefing), doc checker and roadmap linter exit 0, TASK REVIEW presented |

## DEVIATIONS LOG

| Date | Deviation | Disposition |
| --- | --- | --- |
| 2026-09-18 | The containment table says the shell runs only `./py`. The permission layer also allowed shell writes into the room by redirection and a `time` prefix, which run outside the interpreter sandbox; four other commands were refused automatically | Every such call is classified in the transcript audit: all paths stayed in the room ([method note § 6](../findings/m8_11_method_note.md#6-containment-record)) |
| 2026-09-18 | The solver read the output files of three of its own background commands, and the auditor at stage 2 read one of its own oversized tool outputs; the harness writes both to per-session folders outside the room and points the session at them | Each folder held only that session's own output, and no other file outside a room was touched. The route is open in principle; § 3.4 of the clean-room standards should name it |
| 2026-09-18 | A desktop file browser, not an agent, created `.DS_Store` files in the run folder and the auditor's room while the rooms ran | The auditor only listed it; no agent read it. Excluded at landing |
| 2026-09-18 | The stage-2 packet ([`../m8_11/theorem_stage2.md`](../m8_11/theorem_stage2.md)) removed references to the paper and to parent tasks. It rendered "By M8.10's argument" as "By an earlier argument, not supplied here", and it omitted the value paragraphs, including the second variations that E2's nondegeneracy rests on | The auditor graded the `λ₄` lemma a DEFECT for an unstated hypothesis (every non-block level above 6) and noted E2 does not state nondegeneracy. The designer overruled the first and cleared the second, both as packet artifacts, with the reasons in [method note § 5.3](../findings/m8_11_method_note.md#53-the-arguments-graded) |
| 2026-09-18 | The auditor re-imported one stage-1 script at stage 2, which rewrote three of its stage-1 output files | Disclosed by the auditor; byte-identical to its stage-1 clean run, and all 72 stage-1 files match their stage-1 hashes |
| 2026-09-18 | The agents' console logs carried a `.log` extension, which the repository ignores, and their returns were named `RETURN.md`, `AUDIT_STAGE1.md` and `AUDIT_STAGE2.md` | Landed as `*_log.txt` and `*_return.md`, content unchanged. Binary caches (`*.pkl`), the auditor's clean-run copy, its mutation copies and its rerun of the solver are not landed; the unedited bytes of every agent file are hashed in the maintainer's run checkpoints |


## FINDINGS

Full record with the equations, the code map and the audit: [`../findings/m8_11_method_note.md`](../findings/m8_11_method_note.md).

| ID | Finding |
| --- | --- |
| F1 | **Every frozen value reproduces blind.** A1 to D4, N1 and X0 match in both agents: 104 of 104 scripted checks each, by separate implementations that never saw a claimed value. The auditor then confirmed the solver on 490 of 490 paired values and reran its code byte-identically |
| F2 | **The candidate values are now derived.** The auditor computed every value exactly, in multi-quadratic arithmetic with its own coupling coefficients, including the tilt in `√39·ℚ` and `√115·ℚ` and the eighteen new level norms; the solver derived most exactly and identified the rest at two precisions |
| F3 | **Local branch germs exist at all six rays, as an audited argument.** Both agents gave their own complete existence argument, blow-up and division by `a³` before the implicit function theorem, on a slice transverse to the symmetry orbit. The auditor graded T.a to T.e ESTABLISHED step by step. T1, S1 and C2's lemma pass as audited arguments, not as verified theorems |
| F4 | ⚠️ **The `λ₄` lemma's sign step does not stand alone.** "`λ₄ < 0` since `Q ≠ 1`" needs every non-block level above 6, which the frozen text carries only by citing M8.10. The auditor's DEFECT on the transcribed text was overruled as a packet artifact, and the claim passes on both agents' own complete arguments |
| F5 | **The tilt is forced and correct.** Both agents found the order-`a⁵` block equation's tangential part equal to the forcing without the tilt and exactly zero with it, the orbit direction `i·e_t` null, and no `τ_y` component at the prism, with the handout asking the same questions at all six rays |
| F6 | ⚠️ **Several solver checks cannot fail.** Its item 9 reason lines attach labels without reading a verification, its spectral-gap line cannot see a level below 6, and two `results.json` fields are literals. No number is affected, since the auditor reproduced each by its own exact route |
| F7 | **Containment held on the record, with one open route.** No call reached answer-bearing material. The only reads outside a room were a session's own outputs, reached through harness pointers that `--restricted` allows |

## TASK REVIEW (2026-09-18)

Task Duration: 01:31 (from 11:50 to 13:21)
Usage Cap Triggered: NO

Approved by the maintainer on 2026-09-18.

| Result | Status |
| --- | --- |
| Every frozen value, A1 to D4 and N1, reproduced blind: 104 of 104 scripted checks in each agent | ✅ |
| Every value derived exactly by the auditor, including the tilt and the 18 new level norms; the solver confirmed on 490 of 490 paired values and its code rerun byte-identically | ✅ |
| T1, local branch germs at all six rays: the author's proof graded ESTABLISHED at every step, and the auditor's own stage-1 argument establishes the result | ✅ audited argument |
| S1, the prism's symmetry, with the section-level lift derived | ✅ audited argument |
| C2, the `λ₄` lemma: the auditor's DEFECT overruled as a packet artifact, and both agents' own sign arguments complete | ✅ audited argument, ⚠️ F4 |
| Several solver checks cannot fail | ⚠️ no value moves |
| Harness pointers let a room read its own outputs outside the room, and `--restricted` allows it | ⚠️ for the clean-room standards |

| Remaining | Where |
| --- | --- |
| Provenance comparison against the author's package | its landing PR, under `research/scripts/m8_11_author/` |

**Findings.** At the four symmetry-pinned rays and at the pentagonal pyramid and trigonal prism, the formal expansions are Taylor expansions of local branch germs, for sufficiently small amplitude, established as an audited argument rather than a verified theorem. The order-`a²` tilt at the two new rays is now an exact result, derived blind by two agents. No radius, stability or finite-amplitude claim is made.

**Research docs created/updated.** [Task doc](m8_11_task_details.md), [method note](../findings/m8_11_method_note.md), [roadmap](../m8_roadmap.md), [briefing](../../__M8_model_briefing.md), [canonical](../m8_theory_canonical.md), [solver scripts](../scripts/m8_11_solver/), [audit scripts](../scripts/m8_11_audit/), [packet as run](../m8_11/).
