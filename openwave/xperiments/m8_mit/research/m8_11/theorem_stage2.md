# Statements to grade

These are statements about the equation and objects of `worklist.md`, written by someone else. Grade them; do not assume any of them is correct. The notation follows `worklist.md` except where defined here. `U1` to `U6` are the fibre vectors of its § 1.5, and `e_t`, `τ_x`, `τ_y` the tangents of its § 1.6.

## Setting

- The equation `(−Δ − λ)ψ + g|ψ|²ψ = 0` on sections in one sector, `g ≠ 0`.
- The level-6 block `K` (seven-dimensional), and `A = −Δ − 48`, inverted on `K^⊥`.
- `N(ψ) = |ψ|²ψ`, with real derivative `DN_Φ[h] = 2 Re(Φ†h)Φ + |Φ|²h`.
- A fibre vector `u` gives the block section `Φ = Φ_{σ,u}`, normalized as a section. Gradients and Hessians are taken in `⟨·,·⟩_ℝ`.
- `Q([u]) = ∫_X |Φ|⁴`, which is `λ₂/g` at a critical ray.

**Symmetry.** Rotations act by left translation, commuting with `Δ` and `N`, and on block sections through the fibre vector: `R·Φ_{σ,u} = Φ_{σ,D³(R̄)u}`, where `R̄` is the entrywise conjugate; `U(1)` acts by phase. Since `σ` is real, fibre conjugation in a real basis of the sector is antiunitary and commutes with `Δ` and `N`; on the block's fibre vector it is time reversal up to a phase: `Θ(Σ u_m v_m) = Σ (−1)^m conj(u_m) v_{−m}`.

**The fixed space.**
- For a critical ray `[u]`, `H` is its stabilizer among rotations, acting by a character: `h·u = χ(h)u`.
- The twisted fixed space `X_H = {ψ : h·ψ = χ(h)ψ for all h ∈ H}` is closed in every Sobolev space. `K_H = K ∩ X_H` has complex dimension `m`.
- The kernel variable is `k = a·u(θ)`, with `r = a²`. Here `θ` is a coordinate on a local slice of `ℙ(K_H)` through `[u]`, transverse to the continuous symmetries that preserve `ℙ(K_H)`.
- `T` is the slice's tangent space at `[u]`, realized by unit sections `E ∈ K_H` with `⟨Φ, E⟩ = 0`. `Π_T` projects onto it real-orthogonally.

## Theorem T

- **T.a Function spaces.** `F(ψ, λ) = (−Δ − λ)ψ + g|ψ|²ψ` maps `H^{s+2}(E_σ) × ℝ` to `H^s(E_σ)`, for `s > 3/2`, between real Banach spaces. `H^s` is an algebra for `s > 3/2` in dimension 3, so the cubic is real-analytic. `A` is Fredholm of index zero, with kernel `K` and a bounded inverse on `K^⊥`. The nearest other level is 8 in the 4-dimensional sector and 10 in the 3-dimensional sector, a gap of 32 or 72.
- **T.b The fixed space.** `F` maps `X_H × ℝ` into the fixed space of `H^s`, since `N(χψ) = χN(ψ)` for `|χ| = 1`. `N(Φ)` carries `χ` and `A⁻¹` is equivariant, so `ξ ∈ X_H`, and everything below lives in one space.
- **T.c Reduction.** Lyapunov–Schmidt in `X_H` solves the range equation for an analytic `w ∈ X_H ∩ K^⊥`. The reduced equations are the gradient of the reduced functional.

  The phase `−1 ∈ U(1)` sends `k` to `−k`, so the reduced map is odd in `a`. Its radial part divided by `a`, and its tangential part divided by `a³`, are therefore even analytic functions of `a`, hence analytic in `r = a²`:
  - **Radial:** `q(θ, r, λ) = (48 − λ) + g·Q(θ)·r + O(r²)`, with `∂_λ q = −1` at `r = 0`.
  - **Tangential,** divided also by `g`: `G(θ, r) = ¼∇_T Q(θ) + O(r)`, the `¼` coming from the bridge below.

  Dividing by `a` removes the trivial branch `ψ = 0`. Before that division, the radial equation's `λ`-derivative is `−a`, which vanishes at `a = 0`.

  At `(θ₀, 0, 48)`, where `θ₀` is the ray, the Jacobian in `(λ, θ)` is block-triangular, with diagonal blocks `−1` and `¼Hess_T Q(θ₀)`. If `Hess_T Q(θ₀)` is invertible on the slice, the implicit function theorem gives analytic `λ(r)` and `θ(r) = θ₀ + r·θ₂ + O(r²)`. The reduced functional is invariant, so its derivative along the symmetry orbits vanishes, which is why the slice suffices.
- **T.d Statement.** There are a neighborhood `W` of `(0, 48)` in `H^{s+2}(E_σ) × ℝ`, a neighborhood `U` of `θ₀` in the chosen slice and an `ε > 0` such that, for each `a` with `0 < |a| < ε`, exactly one solution `(ψ, λ) ∈ W`, with `ψ ∈ X_H` in the gauge of T.e, has `[Π_K ψ] ∈ U`. `λ` and `ψ/a` are analytic in `r = a²`, and rotations carry this germ to equivalent germs at the conjugate rays. No lower bound on `ε` is given.
- **T.e Identification.** The germ's Taylor coefficients satisfy the formal expansion's order-by-order equations. Those equations have a unique solution under the gauge: `⟨Φ, ψ⟩ = a` real, the range part orthogonal to `K`, and the tilt in the slice. So the formal coefficients are the germ's Taylor coefficients.

**E1: U1 to U4 (`m = 1`).** At each of these the stabilizer's isotypic space in the block is one-dimensional. So `K_H = ℂΦ`, which `U(1)` reduces to one real amplitude. There is no `θ`, nondegeneracy is vacuous, and T gives the germ, whose Taylor coefficients are the formal ones by T.e.

**E2: U5 and U6 (`m = 2`).** Expand `ψ = aΦ + a³(v + ξ) + O(a⁵)`, `λ = 48 + λ₂a² + λ₄a⁴ + O(a⁶)`, where the tilt `v` lies in `K_H`, with `⟨Φ, v⟩ = 0`, in the slice, and `ξ ⊥ K`.
- Order `a³`: `λ₂ = gQ`, and `ξ = −gA⁻¹f` with `f = Π_⊥N(Φ)`.
- The block part of order `a⁵` is `g Π_K DN_Φ[v + ξ] − λ₂v − λ₄Φ = 0`. Its tangential part fixes the tilt: `L_T v = −Π_T DN_Φ[ξ]`, where `L_T = Π_T DN_Φ − Q` on `T`, so `v ∝ g`.
- At U5, the locus `ℙ span{v_2, v_−3}` is fixed by the fivefold rotation about `z`, under which both weights carry the same phase; `R_z(θ)` preserves the locus for every `θ` and moves along it, so the direction `i·e_t` is a rotation orbit, null for `Q`. The slice is the real curve `u(t)`.
- At U6, the locus is the chart `u(z) = v_3 + z·v_0 + v_−3`. `R_z(θ)` preserves the chart only when `e^{6iθ} = 1`, and `R_z(π/3)` sends `z ↦ −z`. The slice is the whole chart, of real dimension two, and `K_H` is the lift of `span{v_3 + v_−3, v_0}`.

## The bridge

Take a unit tangent `E` in the block, real-orthogonal to `Φ`. Then `d²Q(cos s·Φ + sin s·E)/ds² = 4(⟨E, DN_Φ E⟩_ℝ − Q)` at `s = 0`. So `L_T = ¼Hess_T Q`, which is the Jacobian block of T.c: one set of numbers decides both nondegeneracy and the tilt.

## The λ₄ lemma

Pair the order-`a⁵` block equation with `Φ`. The tilt contributes `g⟨Φ, DN_Φ[v]⟩ = g(2X + X̄)`, where criticality gives `X = ⟨N(Φ), v⟩ = ⟨Π_K N(Φ), v⟩ = Q⟨Φ, v⟩ = 0`. So `λ₄ = g⟨Φ, DN_Φ[ξ]⟩ = −3g² Σ_n ‖Π_n N(Φ)‖²/(n(n+2) − 48)`: the tilt does not change `λ₄` at this order. By an earlier argument, not supplied here, `λ₄ < 0` since `Q ≠ 1` at both rays (U5 and U6).

## The S lemma (U6)

`Θ` acts on the chart as `z ↦ −z̄`, and `R_z(π/3)` as `z ↦ −z`. So `S = R_z(π/3)∘Θ` acts as `z ↦ z̄`, and it fixes the U6 points `z = ±√(23/10)`.

`S` lifts to sections as a rotation composed with fibre conjugation. So it maps solutions to solutions, because `−Δ`, `λ` and `g` are real and the conjugate of `|ψ|²ψ` is `|ψ|²ψ̄`. Being antiunitary, it fixes U6's `Φ` once the phase of `Φ` is chosen. Then:
- by uniqueness, `Sξ = ξ`, so the forcing is `S`-even;
- `τ_x` is `S`-even and `τ_y` is `S`-odd, so the forcing has no `τ_y` component, and `L_T`, an `S`-even form, has no cross term.

Hence `v_y = 0` exactly.
