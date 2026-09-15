# Conventions, and a worklist

This file is a conventions extract and a list of questions; nothing below states a value you are asked for. Work from these conventions only. If a question seems underdetermined, say so and state the reading you took rather than choosing one silently.

Work by finite algebra: Clebsch-Gordan coefficients, the representation theory of the finite group generated below, and exact or high-precision linear algebra. Do not solve the nonlinear equation at finite amplitude, by continuation, time-stepping or any other numerical method.

Answer item 0 first, so that a mismatch can be localized before anything else is judged.

---

## 1. Setup

### 1.1 Spin representations

`V_j` is the irreducible `SU(2)` representation of spin `j` and dimension `2j+1`, with orthonormal weight basis `v_m`, `−j ≤ m ≤ j`, and `D^j(g)` is the matrix of `g ∈ SU(2)` in that basis. Clebsch-Gordan coefficients `⟨j₁ m₁; j₂ m₂ | J M⟩` are in the Condon-Shortley convention: real, with `⟨j j; j −j | 2j 0⟩ > 0`. For `x ∈ V_{j₁}` and `y ∈ V_{j₂}`, `[x ⊗ y]_J ∈ V_J` has components `Σ ⟨j₁ m₁; j₂ m₂ | J M⟩ x_{m₁} y_{m₂}`.

Time reversal on `V_3` is the antilinear map `Θ(Σ u_m v_m) = Σ (−1)^m conj(u_m) v_{−m}`.

### 1.2 The group and the quotient

Identify `S³` with `SU(2)` carrying the round metric of radius 1. `Γ ⊂ SU(2)` is the finite subgroup generated, as unit quaternions `w + xi + yj + zk`, by

`q₁ = ½(1 + i + j + k)`,  `q₂ = ½(φ + φ⁻¹ i + j)`,  `φ = (1 + √5)/2`.

`Γ` acts on `S³` by right translation and `X = S³/Γ`. The group is not named here and none of its representation theory is supplied: derive what you need from the two generators. Any standard identification of unit quaternions with `SU(2)` will do; no question depends on which.

### 1.3 Sections, levels and the Laplacian

A unitary representation `σ` of `Γ` on `ℂ^d` determines a flat bundle over `X`, whose sections are the functions `ψ: SU(2) → ℂ^d`, written as rows, with `ψ(gh) = ψ(g) σ(h)` for `h ∈ Γ`. An intertwiner `η ∈ Hom_Γ(σ, V_j)` satisfies `D^j(h) η = η σ(h)`. For each `σ` used, fix one matrix representation and use that same `σ(h)` at every level. Normalize every intertwiner so that `η†η = I`, and where `Hom_Γ(σ, V_j)` has dimension two, keep two intertwiners `η, η′` with `η†η′ = 0`. By Peter-Weyl the sections at level `n = 2j` are `V_j ⊗ Hom_Γ(σ, V_j)`: sums, over the intertwiners kept, of

`ψ_a(g) = Σ_{m,k} u_m D^j_{mk}(g) η_{ka}`,

so the left index is free and the right index carries the sector. The round Laplacian satisfies `−Δψ = n(n+2)ψ` at level `n`. Rotations `r ∈ SU(2)` act on sections by left translation, `(r·ψ)(g) = ψ(r⁻¹g)`. Integrals use the Haar measure of total mass 1, equivalently `∫_X 1 = 1`, and `⟨φ, ψ⟩ = ∫ Σ_a conj(φ_a) ψ_a`, antilinear in the first slot; `⟨φ, ψ⟩_ℝ = Re⟨φ, ψ⟩`. Write `|ψ|² = Σ_a |ψ_a|²`, and `Π_n` for the orthogonal projection onto level `n`.

### 1.4 The two sectors and the block

Restricted to `Γ` through `D^3`, `V_3` has two irreducible constituents, one of dimension 3 and one of dimension 4, each occurring once; call them the 3-dimensional and the 4-dimensional sector. For either, `Hom_Γ(σ, V_3)` is one-dimensional, with normalized generator `η`. The block is the level-6 section space, and a fibre vector `u ∈ V_3` gives the block section `Φ_{σ,u}(g) = Σ_{m,k} u_m D³_{mk}(g) η_{k·}`. Where a question asks for `Φ` normalized, rescale `u` so that `∫_X |Φ_{σ,u}|² = 1`: the section is normalized, not the fibre vector.

### 1.5 Six fibre vectors

`U1 = v_3`, `U2 = v_0`, `U3 = v_2 + v_{−2}`, `U4 = v_3 + v_{−3}`, `U5 = cos t · v_2 + sin t · v_{−3}` with `sin² t = 12/25` and `0 < t < π/2`, and `U6 = v_3 + z₀ v_0 + v_{−3}` with `z₀ = √(23/10)`.

### 1.6 Tangents at U5 and U6

A tangent at a normalized block section `Φ` is a block section `E` with `⟨Φ, E⟩ = 0` and `∫_X |E|² = 1`.
- At U5, `e_t` is the block section of `−sin t · v_2 + cos t · v_{−3}`, normalized, oriented by increasing `t` along `u(t) = cos t · v_2 + sin t · v_{−3}`.
- At U6, `τ_x` is the block section of the component of `v_0` orthogonal to `U6`, normalized, oriented by increasing `x` along `u(z) = v_3 + z v_0 + v_{−3}`, `z = x + iy`. Set `τ_y = i·τ_x`.

For a tangent `E` at `Φ`, the second derivative of a function `f` on fibre vectors along `E` means `d²/ds² f(cos s · û + sin s · ŵ)` at `s = 0`, where `û` and `ŵ` are the fibre vectors underlying the section-normalized `Φ` and `E`, each rescaled to unit norm in `V_3`.

## 2. The equation and the expansion

The equation is `(−Δ − λ)ψ + g N(ψ) = 0` with `N(ψ) = |ψ|²ψ` pointwise, `g` real and nonzero, and `ψ` a section in one sector. Near `λ = 48`, seek formal solutions

`ψ = aΦ + a³(κ + ξ) + a⁵ζ + …`,  `λ = 48 + λ₂a² + λ₄a⁴ + …`,

with `Φ = Φ_{σ,u}` normalized, `a > 0`, `ξ` and `ζ` orthogonal to the block, and the tilt `κ` a block section with `⟨Φ, κ⟩ = 0`, so that `a = ⟨Φ, ψ⟩`. `DN_Φ[h]` is the real derivative `d/dε N(Φ + εh)` at `ε = 0`, for real `ε`.

For fibre vectors: `ρ_K(u) = [u ⊗ Θu]_K` for `0 ≤ K ≤ 6`, and `r̂₆(u) = ‖ρ₆(u)‖²/‖u‖⁴`.

## 3. Worklist

**0.** (a) Confirm `‖q₁‖ = ‖q₂‖ = 1`, compute the order of `Γ`, and report whether `Γ` equals its derived subgroup. (b) Report the value of `⟨3 3; 3 −3 | 6 0⟩`.

**1.** For each sector and each level `n ≤ 18`, `dim Hom_Γ(σ, V_{n/2})`. Say whether you DERIVED these from the generators or RECOGNIZED the group and used known tables, and how. For every intertwiner used in later items, verify at `h = q₁` and `h = q₂` that `D^j(h)η = ησ(h)` for the one `σ` of § 1.3, and that `η†η = I` and, at dimension two, `η†η′ = 0`: report the largest residual of each, with the precision you require.

**2.** For each of U1 to U6 in each sector, with `Φ` normalized: whether `Π_6 N(Φ)` is a multiple of `Φ`, and if so the multiple, exactly. Report `r̂₆` at each of the six, and whether each is a stationary point of `r̂₆` on the unit sphere of `V_3`.

**3.** For each of U1 to U6: the rotations `r ∈ SU(2)` that fix `Φ` up to a phase, the character by which they act on `Φ`, and the complex dimension of the subspace of the block on which they act by that character.

**4.** At U1 to U6 in each sector: solve the order-`a³` equation for `ξ`, and report `‖Π_n ξ‖²/g²` at every level `n` of the sector, exactly, including any level where it vanishes. Say whether `ξ` transforms under the rotations of item 3 by the same character as `Φ`.

**5.** Project the order-`a⁵` equation onto the block. At each of U1 to U6 and in each sector, report the component of `Π_6 DN_Φ[ξ]` along `Φ` and the norm of its component orthogonal to `Φ`, per `g`. At U5, report its `⟨·,·⟩_ℝ` components along `e_t` and `i·e_t`; at U6, along `τ_x` and `τ_y`.

**6.** At U5 and U6 in each sector: the quadratic form `E ↦ ⟨E, DN_Φ E⟩_ℝ − Q` on tangents, with `Q` item 2's multiple, exactly: at U5 on `e_t` and on `i·e_t`; at U6 on `τ_x`, on `τ_y`, and the cross term. Also the second derivative of `r̂₆` along the same tangents, as defined in § 1.6.

**7.** At U1 to U6 in each sector: the tilt `κ` of § 2 for which the order-`a⁵` block equation has no component orthogonal to `Φ`. Where that equation leaves part of `κ` free, take `κ` `⟨·,·⟩_ℝ`-orthogonal to the free directions and say which they are. Report `‖κ‖²/g²` and the `⟨·,·⟩_ℝ` components of `κ/g`, exactly: at U5 along `e_t` and `i·e_t`, at U6 along `τ_x` and `τ_y`, with the norm of any remainder. Then report the equation's component orthogonal to `Φ` with and without `κ`.

**8.** At U1 to U6 in each sector, `λ₄/g²` exactly, computed with `κ` included and again without it, with an exact reason if the two agree; and the sign of `λ₄` for `g ≠ 0`, with an argument rather than only the values.

**9.** For every quantity in items 4 to 8 that vanishes, say whether it was computed and found zero; a value your computation did not reach is reported as missing, never as zero. For each computed zero, give an exact argument for why it vanishes, or report that zero as unresolved. If your reason for a zero is a symmetry, prove that the symmetry preserves both the equation and the sector in which you are working.

**10.** For each of U1 to U6: does an actual solution `(ψ, λ)` of the equation exist for every sufficiently small `a > 0` whose expansion begins as in § 2? Give the argument as prose, list every hypothesis it uses, say where in items 1 to 9 each hypothesis is verified, and for each step state what would fail if it were omitted. If yes, state the local uniqueness your argument proves, including any symmetry quotient and neighborhood restriction; state the regularity of the solution in the amplitude parameter; and identify any degeneracy at `a = 0` that your argument must handle. If you cannot give an argument, say what prevents it.

**11.** On U5's curve `u(t)` at `sin² t = 1/4`, with `Φ` normalized and `e_t` defined as in § 1.6 at that point: is it a stationary point of `r̂₆` on the unit sphere of `V_3`? Report `Re⟨Φ, DN_Φ[e_t]⟩` in each sector, exactly.

**12.** At U1 to U6, check the order-`a³` equation `(−Δ − 48)ξ + g(N(Φ) − Π_6 N(Φ)) = 0` as a vector identity whose two sides come from different operations. Obtain `−Δξ` by applying the `SU(2)` Casimir, built from explicit spin-matrix generators, to the `ξ` of item 4; substituting `n(n+2)` does not count. Evaluate `N(Φ)` pointwise at enough group elements and project it onto every section level through 18, keeping every copy. Say why `N(Φ)` has no component above level 18. Report the residual with the precision you require; a residual above it is reported as unresolved, and naming its cause requires showing it.

## 4. What to return

For each item, the values, exact where the item asks, with the method. Where an item asks for an exact value, either derive it symbolically or, for an identification, state the precision, the denominator bound and the field, and repeat the identification at a second precision. Scripts that reproduce every number from a clean directory. A consulted-material manifest listing every source you read, with item 1's DERIVED or RECOGNIZED declaration. Any question you found underdetermined, with the reading you took.
