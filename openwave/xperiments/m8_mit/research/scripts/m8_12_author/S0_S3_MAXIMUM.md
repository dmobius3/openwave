# At spin 3, the hexagon is the unique maximum of the top multipole

2026-09-19. A lemma note for the S0 search. It closes the one conditional item in the S0 memo. Every step is checked exactly, or armed numerically, by `check_s3_maximum.py` (30 checks, 0 failures; log `check_s3_maximum_log.txt`). It computes no Hessian and no index.

## Statement

For unit `u` in `V₃`, the top density multipole `ρ₆(u) = [u ⊗ Θu]₆` satisfies

`1/924 ≤ ‖ρ₆(u)‖² ≤ 463/924`.

The minimum is attained exactly on coherent states. The maximum is attained exactly on the orbit of `(v₃ + v₋₃)/√2` under rotations and phase: the hexagon, which quantum optics calls the NOON state. **Status: DERIVED.**

The minimum was already derived in the S0 memo. The maximum is new relative to every source read for S0. Romero et al. (2024) pose the general-`S` question and report that no random constellation beats NOON; at `S = 3` this note settles it. That is a search report, not a novelty claim: no search targeted this statement.

## Proof of the maximum

**Step 1: the invariant form.** The `U(1)×SU(2)`-invariant quartics on `V₃` form a four-dimensional space, because `Sym²V₃ = V₀ ⊕ V₂ ⊕ V₄ ⊕ V₆` is multiplicity-free. Take four functions on the unit sphere, with Kawaguchi and Ueda's definitions of the last three:

- the constant 1;
- the magnetization squared `|f|²`, where `f = ⟨u, f u⟩`;
- the singlet-pair weight `|a₀₀|²`, where `a₀₀ = ⟨Θu, u⟩/√7`;
- the nematic invariant `TrN̄²`, where `N̄ᵢⱼ = Re⟨u, fᵢfⱼu⟩`.

All four are invariant quartics. Their values at the weight states `v₃, v₂, v₁, v₀` form a matrix with determinant `180/7 ≠ 0`, so they are a basis. The weight-state values of `‖ρ₆‖²` (1, 36, 225 and 400 over 924) then force

`‖ρ₆‖² = −5/231 − |f|²/22 + (7/11)|a₀₀|² + TrN̄²/198`.

**Step 2: two easy bounds.** `|f|² ≥ 0`, with equality iff `⟨f⟩ = 0`. By Cauchy–Schwarz, `|a₀₀|² ≤ 1/7`, with equality iff `Θu ∝ u`, that is, iff `u` is time-reversal invariant.

**Step 3: the nematic bound, `TrN̄² ≤ 171/2`.** Kawaguchi and Ueda's review states this bound without proof. Write `N̄ = 4I + Q`, with `Q` traceless. Then `TrN̄² = 48 + ‖Q‖²`. For any density matrix,

`‖Q(ρ)‖ = max_E tr(ρ A_E)`, taken over unit real symmetric traceless `E`, where `A_E = Σ Eᵢⱼ{fᵢ, fⱼ}/2`.

Rotating `E` to diagonal form `e = (e₁, e₂, e₃)` does not change `λ_max(A_E)`. So the bound is equivalent to

`λ_max(e₁fₓ² + e₂f_y² + e₃f_z²) ≤ 15/√6` for every unit traceless `e`.

Set `e₃ = 2a/3` and `e₁ − e₂ = 4b`. The constraint becomes `(2/3)a² + 8b² = 1`, and the operator becomes `a(f_z² − 4) + b(f₊² + f₋²)`. The parity of `m` and the flip `m ↦ −m` split it into a zero eigenvalue and three 2×2 blocks:

`M₁ = [[5a, √60·b], [√60·b, −3a + 12b]]`, `M₂` = `M₁` at `−b`, `M₃ = [[0, √240·b], [√240·b, −4a]]`.

The three blocks are bounded one at a time:

- **`M₃`:** `λ_max = −2a + √(30 − 16a²)`. Also `(15/√6 + 2a)² − (30 − 16a²) = 20(a + √6/4)² ≥ 0`, and `15/√6 + 2a ≥ 9/√6 > 0` on the ellipse. So `λ_max(M₃) ≤ 15/√6`.
- **`M₁`:** `λ_max = a + 6b + √(16a² − 48ab + 96b²)`. The right side of `15/√6 − a − 6b ≥ √(…)` is positive, since `a + 6b ≤ √6` on the ellipse, so squaring is valid. Write `N = √((2/3)a² + 8b²)`, which equals 1 on the ellipse. After squaring, and homogenizing with `(15/√6)²N² = 25a² + 300b²`, the bound reads `(3/√6)(a + 6b)·N ≤ a² + 6ab + 24b²`. This holds trivially when `a + 6b ≤ 0`, because the right side is `(a + 3b)² + 15b² ≥ 0`. When `a + 6b > 0`, squaring again gives `(3/2)(a + 6b)²N² ≤ (a² + 6ab + 24b²)²`, and the difference of the two sides is exactly `36b²(a + 2b)² ≥ 0`.
- **`M₂`:** it has `M₁`'s spectrum at `−b`, and the constraint is even in `b`.

So `λ_max ≤ 15/√6` on the whole circle. Tracing the equality cases (`b = 0` with `a > 0`, `a = −2b` with `b > 0`, and `a = −√6/4` in `M₃`) shows equality only at the three permutations of `(2, −1, −1)/√6`. That is, `E = (3nnᵀ − I)/√6` for a unit vector `n`.

Hence `TrN̄² ≤ 48 + 225/6 = 171/2` for every state. Equality holds iff `u` lies in the top eigenspace of `(3(n·f)² − 12)/√6` for some `n`, that is, iff `u ∈ span{|3,3⟩ₙ, |3,−3⟩ₙ}`.

**Step 4: the equality case.** Equality in `‖ρ₆‖² ≤ −5/231 + 1/11 + 171/396 = 463/924` needs all three bounds at once. So `u = α|3,3⟩ₙ + β|3,−3⟩ₙ`. Then `⟨f⟩ = 3(|α|² − |β|²)n`, so `⟨f⟩ = 0` forces `|α| = |β|`, and `|a₀₀|² = 1/7` follows. The remaining relative phase is a rotation about `n`. The equality set is therefore the hexagon's orbit. □

## What the check script does

- **Step 1:** solves for the coefficients in exact rationals. The weight-state inputs come from the aid's own implementation, not only from the paper.
- **Step 3, exactly (sympy):** the operator identity as 7×7 matrices, the constraint, the factorization of the characteristic polynomial into the three blocks, the square completion, the sum of squares, and the side conditions. One arm perturbs the sum of squares (`24b² → 23b²`); another shows that `M₁` and `M₂` themselves differ, so the `b ↦ −b` step is needed.
- **Step 3, numerically:** a 20001-point sweep of the circle, with equality exactly at the three axis angles and a margin of more than `10⁻⁵` below `15/√6` everywhere more than 0.01 rad from an axis. The arm is a perturbed operator (`+0.02 fₓ²`), which the same sweep must reject.
- **Step 4:** checks the equality case on the computed objects.
- **Arithmetic:** checks the arithmetic the note quotes: `−5/231 + 1/11 + 171/396 = 463/924`, `48 + 225/6 = 171/2`, and the `M₃` equality point.

Two first drafts of the numerical gates failed. The equality check was limited by the grid, and one arm was a constant `True`. Both were replaced before this note was written.

## What it changes

**In S.** The hexagon is the global maximum of `r̂₆`, so its Hessian has no positive transverse direction. The exact transverse nullity, and therefore the exact Morse index, remains for S1; the exploratory count supports nondegeneracy modulo symmetry (S0 memo, § 3). With the coherent minimum, both global extrema of the reduced quartic are now derived.

**In the MIT paper (Blake's call).** The fourth bedrock paper's § 7.2 records [RK]'s top-multipole question as open. At `S = 3` it now has an answer, and the proof is about a page long.
