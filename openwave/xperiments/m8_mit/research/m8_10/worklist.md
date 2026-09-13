# Conventions, and a worklist

This file is a conventions extract plus a list of questions. It carries the definitions and none of the results: nothing below states a value you are asked for. Work from these conventions only. If a question seems underdetermined by what is written here, say so and state the reading you took rather than choosing one silently.

Work by finite algebra: Clebsch-Gordan coefficients, the representation theory of the finite group generated below, and exact or high-precision linear algebra. Do not solve the nonlinear equation at finite amplitude, by continuation, time-stepping or any other numerical method. Every question is about the linear equations of successive orders.

Item 0 exists so that a mismatch can be localized before anything else is judged. Answer it first.

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

A unitary representation `σ` of `Γ` on `ℂ^d` determines a flat bundle over `X`, whose sections are the functions `ψ: SU(2) → ℂ^d`, written as rows, with `ψ(gh) = ψ(g) σ(h)` for `h ∈ Γ`. An intertwiner `η ∈ Hom_Γ(σ, V_j)` satisfies `D^j(h) η = η σ(h)`. By Peter-Weyl the sections at level `n = 2j` are `V_j ⊗ Hom_Γ(σ, V_j)`, realized as

`ψ_a(g) = Σ_{m,k} u_m D^j_{mk}(g) η_{ka}`,

so the left index is free and the right index carries the sector. The round Laplacian satisfies `−Δψ = n(n+2)ψ` at level `n`. Integrals use the Haar measure of total mass 1, equivalently `∫_X 1 = 1`, and `⟨φ, ψ⟩ = ∫ Σ_a conj(φ_a) ψ_a`, antilinear in the first slot. Write `|ψ|² = Σ_a |ψ_a|²`, and `Π_n` for the orthogonal projection onto level `n`.

### 1.4 The two sectors and the block

Restricted to `Γ` through `D^3`, `V_3` has two irreducible constituents, one of dimension 3 and one of dimension 4, each occurring once; call them the 3-dimensional and the 4-dimensional sector. For either, `Hom_Γ(σ, V_3)` is one-dimensional; normalize its generator so that `η†η = I_σ`. The block is the level-6 section space, and a fibre vector `v ∈ V_3` gives the block section `Φ_{σ,v}(g) = Σ_{m,k} v_m D³_{mk}(g) η_{k·}`. Where a question asks for `Φ` normalized, rescale `v` so that `∫_X |Φ_{σ,v}|² = 1`: the section is normalized, not the fibre vector.

### 1.5 Five fibre vectors

`R1 = v_3`, `R2 = v_0`, `R3 = v_2 + v_{−2}`, `R4 = v_3 + v_{−3}`, and `R5 = cos t · v_2 + sin t · v_{−3}` with `sin² t = 12/25` and `0 < t < π/2`.

## 2. The equation and the expansion

The equation is `(−Δ − λ)ψ + g N(ψ) = 0` with `N(ψ) = |ψ|²ψ` pointwise, `g` real and nonzero, and `ψ` a section in one sector. Near `λ = 48`, seek formal solutions

`ψ = aΦ + a³ξ + a⁵ζ + …`,  `λ = 48 + λ₂a² + λ₄a⁴ + …`,

with `Φ = Φ_{σ,v}` normalized, `a > 0`, and `ξ` and `ζ` orthogonal to the block, so that `a = ⟨Φ, ψ⟩`. `DN_Φ[h]` is the real derivative `d/dε N(Φ + εh)` at `ε = 0`, for real `ε`.

For the multipole maps used in item 2: `ρ_K(v) = [v ⊗ Θv]_K` and `M_K(v) = [ρ_K(v) ⊗ v]_3`, for `0 ≤ K ≤ 6`.

## 3. Worklist

**0.** (a) Confirm `‖q₁‖ = ‖q₂‖ = 1`, compute the order of `Γ`, and report whether `Γ` equals its derived subgroup. (b) Report the value of `⟨3 3; 3 −3 | 6 0⟩`.

**1.** For each sector and each level `n ≤ 18`, `dim Hom_Γ(σ, V_{n/2})`. Say whether you DERIVED these from the generators or RECOGNIZED the group and used known tables, and how.

**2.** For a general, unnormalized `v ∈ V_3`, the fibre vector `w(v)` with `Π_6 N(Φ_{σ,v}) = Φ_{σ,w(v)}`, in closed form, in both sectors. Express it through the `M_K` above, with every coefficient exact; if you work in another basis, give its exact relation to the `M_K`.

**3.** For each of `R1` to `R5` in each sector, with `Φ` normalized: whether `Π_6 N(Φ)` is a multiple of `Φ`, and if so the multiple, exactly.

**4.** At `R1` to `R4`: solve the order-`a³` equation for `ξ`, and report `‖Π_n ξ‖²/g²` at every level `n` of the sector, exactly, including any level where it vanishes. State your method, and for rational identification the precision and the denominator bound. At `R5`, solve the same equation for `ξ` only as far as item 6 needs it.

**5.** For every level at which item 4 found `‖Π_n ξ‖ = 0`, give an exact argument for why it vanishes. If you cannot establish one, report that zero as unresolved.

**6.** Project the order-`a⁵` equation onto the block. At `R1` to `R5` and in each sector, report the component of `Π_6 DN_Φ[ξ]` along `Φ` and its component orthogonal to `Φ`.

**7.** At `R1` to `R4` and in each sector, `λ₄/g²` exactly, with the method; and the sign of `λ₄` for `g ≠ 0`, with an argument rather than only the values.

**8.** At `R1` to `R4`, the ratio of `λ₄` in the 3-dimensional sector to `λ₄` in the 4-dimensional sector, exactly. Then say whether one number `c` gives `λ₄(3-dimensional) = c · λ₄(4-dimensional)` at all four, and answer the same question for item 3's multiples minus 1.

**9.** At `R1` to `R4`, check the order-`a³` equation for `ξ` as a vector identity: build the `SU(2)` Casimir from explicit spin-matrix generators, apply it to the `ξ` you computed to obtain `−Δξ`, and report the numerical residual. Substituting the formula `n(n+2)` does not count as the check.

## 4. What to return

For each item, the values, exact where the item asks, with the method. Scripts that reproduce every number from a clean directory. A consulted-material manifest listing every source you read, with item 1's DERIVED or RECOGNIZED declaration. Any question you found underdetermined, with the reading you took.
