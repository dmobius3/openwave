# Conventions, and a worklist

Everything below happens in one seven-dimensional complex space. No quotient, no sections, no bundles are needed, and none are supplied.

## 1. Setup

### 1.1 The space

`V₃ = ℂ⁷` with basis `v₃, v₂, v₁, v₀, v₋₁, v₋₂, v₋₃`, written `u = Σ_m c_m v_m`. The Hermitian product is `⟨u, w⟩ = Σ_m conj(c_m) d_m`, and the real inner product is `Re⟨u, w⟩`.

### 1.2 The rotation action

`J_z v_m = m·v_m`, and `J_± v_m = √(12 − m(m ± 1))·v_{m±1}`, with `J_x = (J_+ + J_−)/2` and `J_y = (J_+ − J_−)/2i`. A rotation by angle `θ` about a unit axis `n` acts by `D³(n, θ) = exp(−iθ·(n_x J_x + n_y J_y + n_z J_z))`. This is a representation of `SO(3)`: the rotation by `2π` acts as the identity.

### 1.3 Time reversal and the quartic

`(Θu)_m = (−1)^{3−m}·conj(c_{−m})`.

With the Clebsch-Gordan coefficients `⟨3 m₁; 3 m₂ | 6 Q⟩` in the Condon-Shortley convention, fixed by `⟨3 3; 3 3 | 6 6⟩ = +1`, define, for `Q = −6, …, 6`,

`ρ₆(u)_Q = Σ_{m₁} ⟨3 m₁; 3 (Q − m₁) | 6 Q⟩ · c_{m₁} · (Θu)_{Q − m₁}`,

and then

`r̂₆(u) = Σ_Q |ρ₆(u)_Q|² / ‖u‖⁴`.

`r̂₆` is real, invariant under `u ↦ e^{iφ}u` and under every rotation, and homogeneous of degree 0.

### 1.4 Fixed spaces

For a closed subgroup `H` of `SO(3)` and a one-dimensional character `χ : H → U(1)`, write

`V₃^{(H, χ)} = {u ∈ V₃ : D³(h)u = χ(h)·u for every h in H}`.

This is the fixed space of the subgroup `{(χ(h)⁻¹, h) : h ∈ H}` of `U(1) × SO(3)`, acting by `u ↦ χ(h)⁻¹·D³(h)u`. Its image in `ℙ(V₃)` has projective dimension one less than its complex dimension.

### 1.5 Tangents and the transverse space

For a unit `u`:

- `T_u = {e ∈ V₃ : Re⟨u, e⟩ = 0}`, of 13 real dimensions;
- `O_u = span_ℝ{i·u, −i·J_x u, −i·J_y u, −i·J_z u}`, the phase and rotation directions;
- `N_u = T_u ⊖ O_u`, the orthogonal complement of `O_u` inside `T_u`, with respect to `Re⟨·,·⟩`.

For `e` in `T_u`, write `H_u(e, e) = d²/ds² r̂₆(u·cos s + e·sin s)` at `s = 0`, extended to a symmetric bilinear form by polarization.

### 1.6 One input from elsewhere

In each of two sectors the quantity actually wanted is `Q_σ = 1 + w₆(σ)·r̂₆`, with `w₆ = 28/39` in one sector and `21/52` in the other. Both are positive. You are not asked to derive this; it is given.

## 2. Worklist

Answer in order. Where an item asks for an argument, give it as prose with every hypothesis listed, and for each step say what would fail if that step were omitted.

**0.** Report `⟨3 3; 3 −3 | 6 0⟩`, and `Θ(Θu)` in terms of `u`. Report `r̂₆` at `u = 2v₃ + v₁ − 3v₋₂` and at `u = v₃ + i·v₀ + 2v₋₁`, exactly, and say what precision your route carries.

**1.** Classify the fixed spaces of § 1.4 of complex dimension 1 and of complex dimension 2, up to rotation. Say how you enumerated the subgroups and their characters, whether you covered the continuous subgroups, and what would fail if a conjugacy class of subgroups were omitted. For each class you keep, give a representative space and say which subgroup and character produce it. Report the complex dimensions of the fixed spaces you found that are neither 1 nor 2. Show that your classes are pairwise distinct, that is, that no rotation carries one onto another. Wherever two of them have the same stabilizer, give that separation exactly rather than numerically, and say for each separation which it is.

**2.** For each class of complex dimension 2 from item 1: restrict `r̂₆` to it, exactly, in whatever chart you choose, stating the chart and the point it omits. Then give the complete critical set of the restriction, exactly, including the omitted point. State the argument by which the set is complete, say whether that argument is an elimination or a call to a solver, and say what would fail if you relied on a solver alone. If you use an index sum or Euler characteristic check, say what it can and cannot detect.

**3.** For each class of complex dimension 1 from item 1: is it critical for `r̂₆` on the unit sphere of `V₃`? Give the value of `r̂₆` there, exactly, and the argument for criticality.

**4.** Unite the critical points of items 2 and 3 modulo rotations and phase. Report the distinct orbits with their exact values of `r̂₆`. Wherever two of your critical points share a value, decide whether they are one orbit or two, and say what decided it: a rotation you exhibit, or a rotation-invariant quantity that separates them.

**5.** At each orbit of item 4: report `dim N_u`, the signature `(n₋, n₀, n₊)` of `H_u` restricted to `N_u`, and the characteristic polynomial of that restriction, exactly. State and prove the formula you use for `H_u`, and say where the degree-0 homogeneity of `r̂₆` enters. Check that `H_u` annihilates `O_u`, report the residual, and exhibit a point where that check fails.

**5b.** For each orbit of item 4 that is an interior point of a class of complex dimension 2: take a real basis of the tangent space to that class at the point, inside `T_u`, and say how many independent directions that is. For each of them, first say whether its projection onto `N_u` is zero. Where it is zero, say which direction of `O_u` it is, and report `H_u` along the corresponding unit unprojected tangent, as an orbit-null control. Normalize the rest in `Re⟨·,·⟩` and report `H_u` on them exactly, every diagonal entry and every off-diagonal entry, with their Gram matrix in `Re⟨·,·⟩`, so that the off-diagonal entries are unambiguous.

**6.** Interpret the signatures. For a reduced energy proportional to `g·r̂₆` on the unit sphere, report the Morse index at each orbit of item 4 for `g > 0` and for `g < 0`, and say which orbits are local extrema for each sign.

**7.** For any orbit whose restricted form has a kernel: identify those kernel directions with something structural, such as a direction inside one of the loci of item 1 or a symmetry, or report the kernel as unresolved. Say what distinguishes a genuine kernel from one produced by your numerics.

**8.** Over the whole unit sphere of `V₃`, is the minimum of `r̂₆` attained at one of your orbits? Is the maximum? Answer each with an argument rather than a search, list every hypothesis, and say where each hypothesis is verified. If you cannot give an argument, say what prevents it and what a search alone does and does not establish.

**9.** For every quantity in items 2 to 7 that vanishes, say whether it was computed and found zero. A value your computation did not reach is reported as missing, never as zero. For each computed zero, give an exact reason, or report that zero as unresolved. If your reason is a symmetry, show that the symmetry preserves `r̂₆`.

**10.** Using § 1.6, state what changes in your answers to items 4, 5, 5b and 6 when `r̂₆` is replaced by `Q_σ` in each of the two sectors, and why. Report the numbers of item 5b multiplied by each `w₆`, exactly.

**11.** At `u = (v₃ + v₁ − v₋₂)/√3`: report `924·r̂₆` exactly, the norm of the tangential gradient of `r̂₆` on the unit sphere, and `‖M_u d‖` for each of the four generators of `O_u` exactly as § 1.5 writes them, without normalizing, where `M_u = Hess N(u) − 4N(u)·I` is the matrix you used in item 5. Say what this tells you about the checks in item 5.

**12.** Pick one orbit from item 4 and one tangent direction in its `N_u`. Compute `H_u` along it in two ways: as the second derivative of `s ↦ r̂₆(u cos s + e sin s)` at `s = 0`, and by the matrix formula of item 5. Report the largest disagreement and the precision of each route.

## 3. What to return

One markdown file with:

- a table per item, with exact values as rationals, or as algebraic numbers with their minimal polynomials;
- for every numerical route, the precision it carries and the residual it reached;
- every argument requested in items 1, 2, 3, 5, 7, 8 and 9, written out, with its hypotheses;
- the tangent normalization of item 5b, stated explicitly;
- every "what would fail if omitted" answer;
- the tools you used, and anything you looked up rather than derived.

Report disagreement with your own earlier steps rather than smoothing it. If an item cannot be completed, say which and why; an incomplete item reported as incomplete is worth more than a filled one that is not verified.
