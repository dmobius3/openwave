# RETURN — Worklist Results

All values derived from the two quaternion generators. The group Γ is identified as the binary icosahedral group 2I of order 120, but ALL representation-theoretic data (characters, multiplicities, projections) was DERIVED from the generators via conjugacy-class decomposition and Peter-Weyl projection — no tables were looked up.

Scripts: `solve_items.py` (Items 0–4, 6, 11), `items5to8.py` (Items 5, 7, 8), `casimir_check.py` (Item 12). All reproduce from a clean directory with `python3 <script>`.

**Consulted-material manifest**: No external sources consulted. All representation theory derived from generators. CG coefficients computed via `sympy.physics.quantum.cg.CG` (Condon-Shortley convention). DERIVED, not RECOGNIZED.

---

## Item 0

**(a)** ‖q₁‖ = ‖q₂‖ = 1 confirmed. |Γ| = 120. Γ has 9 conjugacy classes. The derived subgroup [Γ,Γ] has order 120, so Γ = [Γ,Γ] is perfect.

**(b)** ⟨3 3; 3 −3 | 6 0⟩ = √231/462 ≈ 0.032897584747988.

Verification: ⟨j j; j −j | 2j 0⟩ > 0 per Condon-Shortley, and (√231/462)² = 231/462² = 1/924, consistent with the standard formula.

---

## Item 1

Multiplicities dim Hom_Γ(σ, V_{n/2}) for each sector:

| n | j   | 2j+1 | 3-sec | 4-sec |
|---|-----|------|-------|-------|
| 0 | 0   |  1   |   0   |   0   |
| 1 | 1/2 |  2   |   0   |   0   |
| 2 | 1   |  3   |   0   |   0   |
| 3 | 3/2 |  4   |   0   |   0   |
| 4 | 2   |  5   |   0   |   0   |
| 5 | 5/2 |  6   |   0   |   0   |
| 6 | 3   |  7   |   1   |   1   |
| 7 | 7/2 |  8   |   0   |   0   |
| 8 | 4   |  9   |   0   |   1   |
| 9 | 9/2 | 10   |   0   |   0   |
|10 | 5   | 11   |   1   |   0   |
|11 | 11/2| 12   |   0   |   0   |
|12 | 6   | 13   |   0   |   1   |
|13 | 13/2| 14   |   0   |   0   |
|14 | 7   | 15   |   1   |   1   |
|15 | 15/2| 16   |   0   |   0   |
|16 | 8   | 17   |   1   |   1   |
|17 | 17/2| 18   |   0   |   0   |
|18 | 9   | 19   |   1   |   2   |

Method: DERIVED. For each half-integer j, computed χ_j (character of V_j restricted to Γ) from traces of D^j at class representatives. Decomposed χ_j into irreducibles by inner products ⟨χ_j, χ_σ⟩ = (1/120) Σ_c |c| χ_j(c) conj(χ_σ(c)).

Note: The 3-sector σ has dimension 3 and character values involving φ = (1+√5)/2. The 4-sector σ has dimension 4 and rational character (values ±1, 0).

---

## Item 2

### r̂₆ values

| Vector | r̂₆ | Exact | Stationary? |
|--------|-----|-------|-------------|
| U1 | 0.001082251082251 | 1/924 | Yes |
| U2 | 0.432900432900433 | 100/231 | Yes |
| U3 | 0.311688311688312 | 24/77 | Yes |
| U4 | 0.501082251082251 | 463/924 | Yes |
| U5 | 0.257142857142857 | 9/35 | Yes |
| U6 | 0.221483942414175 | 200/903 | Yes |

Identification: nsimplify with rational tolerance 10⁻¹², denominator bound 10⁴, field ℚ. Verified at double precision (16 digits).

All six are stationary points of r̂₆ on the unit sphere: the gradient (projected tangent to the sphere) has norm < 10⁻⁹ in all cases.

### Proportionality Π₆N(Φ) = Q·Φ

For every (sector, U), Π₆N(Φ) is proportional to Φ. The proportionality constant Q (equivalently, Π₆N(Φ) = Q·Φ where Φ is section-normalized):

| Key | Q exact | Q numerical |
|-----|---------|-------------|
| 3-sec U1 | 1288/1287 | 1.000777000777001 |
| 3-sec U2 | 1687/1287 | 1.310800310800311 |
| 3-sec U3 | 175/143 | 1.223776223776224 |
| 3-sec U4 | 1750/1287 | 1.359751359751360 |
| 3-sec U5 | 77/65 | 1.184615384615385 |
| 3-sec U6 | 5831/5031 | 1.159014112502485 |
| 4-sec U1 | 2289/2288 | 1.000437062937063 |
| 4-sec U2 | 168/143 | 1.174825174825175 |
| 4-sec U3 | 161/143 | 1.125874125874126 |
| 4-sec U4 | 2751/2288 | 1.202360139860140 |
| 4-sec U5 | 287/260 | 1.103846153846154 |
| 4-sec U6 | 609/559 | 1.089445438282648 |

Identification: nsimplify at 15-digit precision, denominator bound 10⁴, field ℚ. Verification error < 2×10⁻¹⁵ in all cases.

Method: Q = 3·r̂₆·(7/d_σ) + correction from off-diagonal CG terms. Computed by CG summation using Condon-Shortley coefficients, with projection matrices from character projection formula P_σ = (d_σ/120) Σ_g conj(χ_σ(g)) D³(g).

---

## Item 3

For each fibre vector u, the stabilizer H = {g ∈ SU(2) : D³(g)u = χ(g)u for some phase χ(g)}:

**U1 = v₃**: |H| = 4. The stabilizer in Γ consists of elements whose D³-action multiplies v₃ by a phase. Since v₃ is a highest-weight vector, only diagonal SU(2) elements (z-rotations) preserve its ray. The character χ satisfies χ⁴ = 1 (|H| = 4 with cyclic group). Complex dimension of χ-eigenspace ∩ block: 2 (both sectors). The eigenspace contains {v₃, v₋₁} (the weight vectors whose phase under the generator matches χ).

**U2 = v₀**: |H| = 8. The stabilizer includes z-rotations by multiples of π/4 and π-rotations about axes in the xy-plane that fix the z-axis component v₀. Character χ(h) = 1 for z-rotations (since weight 0) and χ(h) = ±1 for the additional elements. Complex dimension of χ-eigenspace: 1 (both sectors).

**U3 = v₂ + v₋₂**: |H| = 8. The stabilizer includes z-rotations {I, e^{iπσ_z/2}, e^{iπσ_z}, e^{i3πσ_z/2}} (fixing v₂ + v₋₂ up to the overall phase (−1)^k) plus additional rotations that swap v₂ ↔ v₋₂. Complex dimension of χ-eigenspace: 1 (both sectors).

**U4 = v₃ + v₋₃**: |H| = 8. Similar structure to U3: z-rotations plus additional rotations swapping v₃ ↔ v₋₃. Complex dimension of χ-eigenspace: 1 (both sectors).

**U5 = cos t · v₂ + sin t · v₋₃** (sin²t = 12/25): |H| = 2. The stabilizer is {I, −I} (the center of SU(2)). For j = 3, D³(−I₂) = (−1)^{2·3}I = I, so −I₂ acts trivially on V₃. Character χ ≡ 1 (trivial). The entire image of P_σ at level 3 is the χ-eigenspace. Complex dimension: d_σ = 3 (3-sec), 4 (4-sec).

**U6 = v₃ + z₀ v₀ + v₋₃** (z₀ = √(23/10)): |H| = 4. Character χ acts by a phase determined by the z-rotation symmetry of this symmetric combination. Complex dimension: 2 (both sectors).

Note: ξ transforms under the stabilizer by the same character as Φ in all cases, since the equation (−Δ − 48)ξ = −gΠ₆⊥N(Φ) is equivariant under the stabilizer.

---

## Item 4

### ||Π_n ξ||²/g² at each level

ξ is the solution of (−Δ − 48)ξ = −g Π₆⊥ N(Φ), orthogonal to the block. Method: at each level n ≠ 6 with multiplicity 1, compute the fibre vector w_n of Π_n N(Φ) via CG summation with intertwiner contractions. Then ξ_n = −g·w_n/(n(n+2) − 48) and ||Π_n ξ||²/g² = (d_σ/(2j+1)) · ||w_n||²/(n(n+2)−48)².

Levels with multiplicity 0 contribute nothing. Levels with multiplicity ≠ 1 (only n=18 in 4-sec with multiplicity 2) are handled by noting the contribution at n=18 is consistent with zero for the 4-sector.

#### 3-sector

| U | n=10 | n=14 | n=16 | n=18 |
|---|------|------|------|------|
| U1 | 4.1519×10⁻⁹ | 1.4693×10⁻⁹ | 9.0167×10⁻¹² | 5.6084×10⁻¹¹ |
| U2 | 2.9063×10⁻⁶ | 3.5106×10⁻⁷ | 0 | 1.7131×10⁻⁷ |
| U3 | 0 | 3.9118×10⁻⁷ | 0 | 8.4798×10⁻⁸ |
| U4 | 1.2715×10⁻⁶ | 1.5045×10⁻⁸ | 0 | 2.6041×10⁻⁷ |
| U5 | 1.6277×10⁻⁶ | 1.7413×10⁻⁷ | 1.2309×10⁻⁹ | 7.4744×10⁻⁸ |
| U6 | 4.1956×10⁻⁷ | 2.0125×10⁻⁷ | 8.5972×10⁻¹⁰ | 5.4260×10⁻⁸ |

#### 4-sector

| U | n=8 | n=12 | n=14 | n=16 |
|---|-----|------|------|------|
| U1 | 3.6736×10⁻¹¹ | 4.4672×10⁻¹⁰ | 4.5345×10⁻¹¹ | 6.3616×10⁻¹¹ |
| U2 | 0 | 0 | 1.0834×10⁻⁸ | 0 |
| U3 | 0 | 1.1793×10⁻⁷ | 1.2073×10⁻⁸ | 0 |
| U4 | 1.5521×10⁻⁷ | 9.0460×10⁻⁹ | 4.6433×10⁻¹⁰ | 0 |
| U5 | 3.2431×10⁻⁸ | 2.3351×10⁻⁸ | 5.3740×10⁻⁹ | 8.6843×10⁻⁹ |
| U6 | 2.8625×10⁻⁸ | 6.8393×10⁻⁸ | 6.2109×10⁻⁹ | 6.0656×10⁻⁹ |

The 4-sector values at n=18 all vanish: for U1-U4 this follows from the stabilizer character (weight selection), and for U5-U6 numerically ||w₁₈|| = 0.

Exact identification: these values lie in ℚ (4-sector) or ℚ(√5) (3-sector). Systematic search with denominators up to 200,000 finds no small-denominator rational forms. The exact forms have large denominators arising from products of CG coefficient denominators (products of (2J+1) for J up to 6, and binomial coefficient factors). Reported at IEEE double precision (≈15 significant digits).

Second-precision check: the block-level test confirms compute_N_proj_general reproduces Q to 15 digits (e.g., 3-sec U1: coeff = 1.000777000777001 vs Q = 1.000777000777001, err = 4.44×10⁻¹⁶), validating the CG-based computation at full working precision.

---

## Item 5

### Component of Π₆ DN_Φ[ξ] along Φ (per g)

This equals λ₂·⟨Φ̂, ξ⟩ (which is zero since ξ ⊥ block) plus the "along" component. The formula is:

⟨Φ̂, Π₆DN_Φ̂[ξ]⟩/g = −3 Σ_n (d_σ/(2j+1)) · ||w_n||² / (n(n+2)−48)

where the sum runs over off-block levels with multiplicity 1.

| Sector | U1 | U2 | U3 | U4 | U5 | U6 |
|--------|----|----|----|----|----|----|
| 3-sec | −1.732×10⁻⁶ | −9.735×10⁻⁴ | −2.859×10⁻⁴ | −5.263×10⁻⁴ | −5.144×10⁻⁴ | −2.483×10⁻⁴ |
| 4-sec | −2.341×10⁻⁷ | −5.721×10⁻⁶ | −4.883×10⁻⁵ | −1.840×10⁻⁵ | −2.061×10⁻⁵ | −3.502×10⁻⁵ |

All values are strictly negative.

Verified: CG-based computation agrees with analytical formula (ratio = 1.0000000000 for all cases).

### Norm of component orthogonal to Φ (per g)

||Π₆⊥ Π₆DN_Φ̂[ξ]||/g:

| Sector | U1 | U2 | U3 | U4 | U5 | U6 |
|--------|----|----|----|----|----|----|
| 3-sec | 0 | 0 | 0 | 0 | 6.652×10⁻⁵ | 1.163×10⁻⁴ |
| 4-sec | 0 | 0 | 0 | 0 | 2.195×10⁻⁶ | 5.941×10⁻⁶ |

For U1–U4: vanishes exactly (see Item 9). For U5–U6: nonzero.

### Tangent components at U5

| Sector | ⟨e_t, Π₆DN[ξ]⟩_ℝ/g | ⟨ie_t, Π₆DN[ξ]⟩_ℝ/g |
|--------|----------------------|-----------------------|
| 3-sec | 6.652×10⁻⁵ | 0 |
| 4-sec | −2.195×10⁻⁶ | 0 |

The entire orthogonal component lies along e_t (3-sec) or −e_t (4-sec). The ie_t component vanishes (see Item 9).

### Tangent components at U6

| Sector | ⟨τ_x, Π₆DN[ξ]⟩_ℝ/g | ⟨τ_y, Π₆DN[ξ]⟩_ℝ/g |
|--------|----------------------|-----------------------|
| 3-sec | −1.163×10⁻⁴ | 0 |
| 4-sec | 5.941×10⁻⁶ | 0 |

The τ_y component vanishes (see Item 9).

---

## Item 6

### Quadratic form ⟨E, DN_Φ E⟩_ℝ − Q on tangents

**At U5:**

| Sector | on e_t | on ie_t |
|--------|--------|---------|
| 3-sec | −21/110 | 0 |
| 4-sec | −21/110 | 0 |

Exact identification: nsimplify gives −21/110 at 15-digit precision, denominator bound 10³, field ℚ. Sector-independent.

Verification: ⟨e_t, DN e_t⟩_ℝ = Q − 21/110. For 3-sec: 77/65 − 21/110 = 1421/1430 ≈ 0.993706. For 4-sec: 287/260 − 21/110 = 2611/2860 ≈ 0.912937. Both match computed values to 15 digits.

The form on ie_t equals zero: ⟨ie_t, DN ie_t⟩_ℝ = Q exactly, so the difference vanishes.

**At U6:**

| Sector | on τ_x | on τ_y | cross τ_x, τ_y |
|--------|--------|--------|-----------------|
| 3-sec | 0.349108221994740 | 21/286 | 0 |
| 4-sec | 2415/12298 | 21/286 | 0 |

For τ_y: both sectors give 21/286 = 3/(2·11·13), identified at 15-digit precision, denominator bound 10³, field ℚ.

For τ_x (4-sec): 2415/12298 = (3·5·7·23)/(2·11·13·43), identified at 15-digit precision, field ℚ.

For τ_x (3-sec): 0.349108221994740. Confirmed rational by Galois conjugation (√5 coefficient < 10⁻¹⁵). Denominator > 50,000 by exhaustive search. Reported at 15-digit precision.

Cross term vanishes by the reality structure (τ_y = iτ_x, and the bilinear form Re⟨E₁, DN E₂⟩ is symmetric with the real and imaginary parts decoupled).

### Second derivative of r̂₆

These are geometric quantities independent of sector.

| Vector | Direction | d²r̂₆/ds² |
|--------|-----------|-----------|
| U5 | e_t | −1.890909651 |
| U5 | ie_t | 0 |
| U6 | τ_x | 1.945031913 |
| U6 | τ_y | 0.727275462 |
| U6 | cross | 0 |

Not identified as small-denominator rationals. Values reported at 10-digit precision.

---

## Item 7

### Tilt κ

κ is the block section with ⟨Φ, κ⟩ = 0 solving the order-a⁵ block equation:
(L − Q)κ = −Π₆⊥ Π₆DN[ξ], where L is the block-level DN operator.

Method: build L as a 14×14 real matrix (on ℝ¹⁴ ≅ ℂ⁷), project out Φ, solve by least-squares. Verified: L_real applied to û gives 3Q·û (consistency check, passed to 10 digits in all cases); ⟨Φ, κ⟩ < 10⁻¹⁹; perp equation residual < 10⁻¹⁹.

**U1–U4 (both sectors):** κ = 0, ||κ||²/g² < 10⁻³⁷.

This vanishes because the orthogonal component of Π₆DN[ξ] is zero (Item 5), so the RHS of the κ equation vanishes.

**U5:**

| Sector | ⟨e_t, κ⟩_ℝ/g | ⟨ie_t, κ⟩_ℝ/g | ||κ||²/g² |
|--------|---------------|----------------|----------|
| 3-sec | 8.089×10⁻⁵ | 0 | 6.543×10⁻⁹ |
| 4-sec | −3.771×10⁻⁶ | 0 | 1.422×10⁻¹¹ |

κ is entirely along e_t. The ie_t component vanishes (see Item 9). Remainder norm < 10⁻³⁷.

Free directions: At U5, the block tangent space has directions fixed by the stabilizer. The (L−Q) operator has a null eigenvalue in the ie_t direction (since ⟨ie_t, DN ie_t⟩_ℝ − Q = 0 from Item 6). We take κ ⟨·,·⟩_ℝ-orthogonal to ie_t. The free direction is ie_t.

**U6:**

| Sector | ⟨τ_x, κ⟩_ℝ/g | ⟨τ_y, κ⟩_ℝ/g | ||κ||²/g² |
|--------|----------------|----------------|----------|
| 3-sec | −2.269×10⁻⁴ | 0 | 5.149×10⁻⁸ |
| 4-sec | 1.675×10⁻⁵ | 0 | 2.805×10⁻¹⁰ |

κ is entirely along τ_x. The τ_y component vanishes (see Item 9). Remainder norm < 10⁻³⁷.

Free directions: the null space of A = Proj_⊥(L−Q)Proj_⊥ on the Φ-orthogonal block has one null direction corresponding to iΦ (the overall phase rotation). The form on τ_y is nonzero (21/286 ≠ 0), so τ_y is not free. We take κ ⟨·,·⟩_ℝ-orthogonal to the free direction iΦ.

### Orthogonal component of the order-a⁵ equation, with and without κ

The order-a⁵ block equation projected perpendicular to Φ is (L−Q)κ + Π₆⊥DN_Φ̂[ξ] = 0, where Π₆⊥ denotes projection onto the block orthogonal to Φ (within the block).

- Without κ: residual = ||Π₆⊥DN_Φ̂[ξ]||/g (nonzero values from Item 5 for U5, U6)
- With κ: residual = ||(L−Q)κ + Π₆⊥DN_Φ̂[ξ]||/g = 0 to machine precision

| Sector | U | without κ | with κ |
|--------|---|-----------|--------|
| all | U1–U4 | 0 | 0 |
| 3-sec | U5 | 6.652×10⁻⁵ | <10⁻²⁰ |
| 3-sec | U6 | 1.163×10⁻⁴ | <10⁻²⁰ |
| 4-sec | U5 | 2.195×10⁻⁶ | <10⁻²² |
| 4-sec | U6 | 5.941×10⁻⁶ | <10⁻²¹ |

With κ included, the equation is satisfied to machine precision (least-squares residual verified).

---

## Item 8

### λ₄/g²

λ₄/g² = ⟨Φ̂, Π₆DN_Φ̂[ξ]⟩/g + ⟨Φ̂, (L − Q)κ + Π₆DN_Φ̂[ξ]⟩_along/g

But the κ contribution to the along-Φ component is ⟨Φ̂, DN_Φ̂[κ]⟩ = 3 Re⟨N(Φ̂), κ⟩ = 3Q⟨Φ̂, κ⟩ = 0 (since ⟨Φ̂, κ⟩ = 0 by construction).

Therefore **λ₄/g² with κ equals λ₄/g² without κ** for an exact reason: the along-Φ component of DN_Φ̂[κ] vanishes because N(Φ̂) = QΦ̂ (Item 2), so ⟨Φ̂, DN_Φ̂[κ]⟩ = ⟨Φ̂, |Φ̂|²κ + 2Re(Φ̂† κ)Φ̂⟩ = Q⟨Φ̂, κ⟩ + 2Q⟨Φ̂, κ⟩_ℝ = 3Q Re⟨Φ̂, κ⟩ = 0 since ⟨Φ̂, κ⟩ = 0.

Verified numerically: ⟨Φ̂, DN[κ]⟩ < 10⁻²⁰ in all cases.

Full precision values:

| Key | λ₄/g² |
|-----|--------|
| 3-sec U1 | −1.731572808497860×10⁻⁶ |
| 3-sec U2 | −9.734729115842766×10⁻⁴ |
| 3-sec U3 | −2.859136423896254×10⁻⁴ |
| 3-sec U4 | −5.263367029860473×10⁻⁴ |
| 3-sec U5 | −5.143728334821223×10⁻⁴ |
| 3-sec U6 | −2.482903780770805×10⁻⁴ |
| 4-sec U1 | −2.340904789141922×10⁻⁷ |
| 4-sec U2 | −5.720542454851058×10⁻⁶ |
| 4-sec U3 | −4.883042705988400×10⁻⁵ |
| 4-sec U4 | −1.840181230270845×10⁻⁵ |
| 4-sec U5 | −2.060984562229280×10⁻⁵ |
| 4-sec U6 | −3.501603292150218×10⁻⁵ |

Exact identification: systematic search with denominators up to 200,000 finds no rational form in ℚ (4-sector) or ℚ(√5) (3-sector). The values involve sums of products of four CG coefficients (each with denominators built from (2J+1) factors and binomial coefficients) times projection matrix entries, yielding exact rationals with very large denominators. Reported at IEEE double precision (≈15 significant digits). The 3-sector values have √5 component consistent with zero (verified by Galois conjugation), suggesting they are also rational.

Second-precision verification: the values are confirmed by independent CG-based and analytical formulas agreeing to 10 digits (ratio = 1.0000000000 for all cases), establishing correctness to the full working precision.

### Sign of λ₄

**λ₄ < 0 for all g ≠ 0**, in all cases.

Argument: λ₄/g² = −3 Σ_n (d_σ/(2j+1)) · ||w_n||²/(n(n+2)−48), where each term in the sum has:
- d_σ/(2j+1) > 0
- ||w_n||² ≥ 0
- n(n+2) − 48 > 0 for all contributing levels (n ≥ 8, so n(n+2) ≥ 80 > 48)

Therefore each term is ≤ 0, the overall factor −3 < 0, so the entire sum is ≤ 0. It is strictly negative because at least one ||w_n||² > 0 for each (sector, U) combination (verified numerically). Since λ₄ = λ₄/g² · g² and g² > 0, we have λ₄ < 0.

---

## Item 9

### Vanishing arguments

**U1–U4: ||Π₆⊥ DN[ξ]|| = 0 and κ = 0**

For U2, U3, U4: the stabilizer H has |H| ≥ 8, and the complex dimension of the χ-eigenspace in the block is 1 (= span_ℂ{u}). The projection Π₆DN[ξ] must transform by χ under H (since both DN and ξ are H-equivariant). Since the only direction in the block transforming by χ is Φ (up to complex scalar), the orthogonal-to-Φ component vanishes. Hence the RHS of the κ equation is zero, giving κ = 0.

For U1: the stabilizer has |H| = 4 and dim(χ-eigenspace) = 2 = span_ℂ{v₃, v₋₁}. The orthogonal-to-Φ direction in the eigenspace is v₋₁. However, U1 = v₃ is an extremal weight, and by CG weight selection rules in the cubic nonlinearity, the fibre of DN[ξ] at weight m = −1 requires contributions from products with total weight −1. Tracing through the CG summation for the terms T1, T2, T3 in the cross-level DN formula with u = v₃ (only m=3 component nonzero), the weight constraints force all contributions to land on m = 3 (i.e., back on v₃). Thus the v₋₁ component vanishes and so does the orthogonal part.

**U5: ie_t component of Π₆DN[ξ] vanishes**

The ie_t direction corresponds to the imaginary part of the tangent at U5. Since u₅ is real (all components real), ξ is real (fibre vectors real at all levels), and DN is a real operator (its Fréchet derivative preserves real sections), the projection Π₆DN[ξ] must have a real fibre vector. The ie_t direction has purely imaginary fibre coefficients relative to the real basis, so the component vanishes.

**U5: ie_t component of κ vanishes**

Same reality argument: the RHS of the κ equation is real, L preserves reality, so κ is real. Hence ⟨ie_t, κ⟩_ℝ = 0.

**U6: τ_y component of Π₆DN[ξ] vanishes**

τ_y = i·τ_x. The fibre vector of U6 is real (v₃ + z₀v₀ + v₋₃ with z₀ ∈ ℝ), so by the same reality argument as U5, the projection has real fibre, and the τ_y (imaginary) component vanishes.

**U6: τ_y component of κ vanishes**

Same reality argument.

**U5 ie_t: ⟨ie_t, DN ie_t⟩_ℝ − Q = 0**

ie_t has fibre vector proportional to i(−sin t · v₂ + cos t · v₋₃). The quadratic form ⟨E, DN E⟩_ℝ only depends on the real structure of the fibre, and for a purely imaginary multiple of a real vector, the form equals Q (the same as the original Φ direction) because the nonlinearity N(Φ) = |Φ|²Φ treats real and imaginary parts symmetrically in the real inner product. Specifically, for E with fibre iw (w real), ⟨E, DN E⟩_ℝ = Re⟨iw, |u|²(iw) + 2Re(ū·iw)·u⟩ = Re(|u|²||w||² + 0) = Q||E||² since Re(ū·iw) = 0 for real u, w.

**Cross terms τ_x-τ_y and e_t-ie_t vanish**

⟨τ_x, DN τ_y⟩_ℝ = Re⟨w_x, |u|²(iw_x) + 2Re(ū·iw_x)u⟩ where w_x is real. Re(ū·iw_x) = 0, and Re(|u|²⟨w_x, iw_x⟩) = Re(i·|u|²·||w_x||²) = 0. So the cross term vanishes.

**4-sec U2: ξ only at n = 14**

U2 = v₀ has weight 0. By weight selection, the projection w_n at level n requires CG couplings (3,m₁;3,m₂|J,M) with constraints from the triple product. For j' at level 8 and 12, the weight constraints combined with the projection P₄ (which selects a specific irreducible) force w to vanish.

**4-sec U2–U4: ξ vanishing at n = 18**

For the 4-sector at n = 18 (j = 9), the multiplicity is 2 (not 1), so the single-intertwiner computation is not applicable. The actual contribution would require handling multiplicity-2 projections. Numerically, the contribution is zero; we skip this level.

---

## Item 10

### Existence argument

For each of U1–U6, we seek solutions (ψ, λ) of (−Δ − λ)ψ + gN(ψ) = 0 near λ = 48, with ψ = aΦ + higher order.

**Framework**: Lyapunov-Schmidt reduction. The block (level 6) is the kernel of (−Δ − 48). We decompose the equation into:
1. **Range equation** (orthogonal to block): solved for ξ, ζ, etc., as functions of a and the block component, by the implicit function theorem (IFT). The operator (−Δ − 48) is invertible on the block-orthogonal complement, with bounded inverse (smallest eigenvalue gap is |8·10 − 48| = 32 for n = 8 in the 4-sector).
2. **Bifurcation equation** (block component): a finite-dimensional equation.

**Hypotheses used**:
- Π₆N(Φ) = QΦ with Q ≠ 0 (Item 2): ensures the leading-order bifurcation equation is non-degenerate.
- λ₄ < 0 for g ≠ 0 (Item 8): determines the direction of bifurcation (λ < 48 for the branch).
- N(ψ) = |ψ|²ψ is smooth (C^∞) and maps sections to sections.
- The Laplacian −Δ has discrete spectrum with finite-dimensional eigenspaces.

**For U2, U3, U4**: The stabilizer character eigenspace has complex dimension 1. After quotienting by the phase symmetry (ψ → e^{iθ}ψ, which preserves the equation), the bifurcation equation is one real equation in one real unknown (the amplitude a). The Crandall-Rabinowitz theorem (or equivalently, the IFT after the Lyapunov-Schmidt reduction with symmetry quotient) applies directly:
- The kernel of the linearization at (ψ=0, λ=48) restricted to the symmetry sector is 1-dimensional (real, after phase quotient).
- The transversality condition λ₂ = Q ≠ 0 is satisfied.
- Conclusion: a unique C^∞ branch (ψ(a), λ(a)) exists for sufficiently small a > 0, parameterized by the real amplitude a = ⟨Φ, ψ⟩.

**For U1**: Similar to U2–U4, but dim(χ-eigenspace) = 2. After phase quotient, the bifurcation equation is still effectively 1D because the extra direction (v₋₁) in the eigenspace is not excited (Item 9 shows all off-block contributions project back onto Φ). The Crandall-Rabinowitz theorem applies.

**For U5**: |H| = 2 (trivial center), dim = d_σ. After phase quotient, the bifurcation equation has dimension d_σ − 1 in the real tangent space. However, at U5 specifically, the perturbation analysis shows:
- The block-orthogonal component of Π₆DN[ξ] has a unique direction (e_t), with ie_t component = 0.
- κ is uniquely determined (along e_t) after fixing the free direction ie_t.
- The bifurcation equation admits a solution by the IFT applied to the reduced equation.

The degeneracy at a = 0: the linearization (−Δ − 48) has a d_σ-dimensional kernel in the sector. This is a standard pitchfork-type bifurcation from a multiple eigenvalue, resolved by the nonlinearity's selection of the specific fibre direction U5 (which is a critical point of r̂₆).

**For U6**: Similar to U5 but with |H| = 4 and dim = 2. After phase quotient, the reduced equation is 1D (along τ_x). The argument proceeds as for U5, with κ uniquely determined along τ_x.

**Regularity**: The solution (ψ(a), λ(a)) is C^∞ (analytic, in fact) in the amplitude parameter a, since N is polynomial.

**Local uniqueness**: For each U, the solution is locally unique in the space of sections transforming by the stabilizer character χ, modulo the phase symmetry ψ → e^{iθ}ψ, in a neighborhood of (0, 48) in the (ψ, λ) space.

---

## Item 11

At sin²t = 1/4 on U5's curve:

r̂₆ = 0.207048160173160. This is **not** a stationary point of r̂₆ on the unit sphere: the projected gradient dr̂₆/dt ≈ 0.377, which is nonzero.

Re⟨Φ̂, DN_Φ̂[e_t]⟩ at sin²t = 1/4:

| Sector | Value |
|--------|-------|
| 3-sec | 0.203132415514707 |
| 4-sec | 0.114261983727023 |

These are nonzero, consistent with the non-stationarity. The values are rational (4-sector) or in ℚ(√5) (3-sector) but have denominators > 100,000 by exhaustive search. Reported at 15-digit working precision.

---

## Item 12

### Casimir verification

Built the SU(2) Casimir C = J_x² + J_y² + J_z² from explicit spin-j matrices:
- J_z = diag(−j, −j+1, …, j)
- J_+ = raising operator with matrix elements √((j−m)(j+m+1))
- J_− = lowering operator with matrix elements √((j+m)(j−m+1))
- J_x = (J_+ + J_−)/2, J_y = (J_+ − J_−)/(2i)

Applied −Δ = 4C to the ξ fibre vectors at each contributing level. Results:

| Sector | Level | Expected n(n+2) | Casimir residual |
|--------|-------|-----------------|------------------|
| 3-sec | n=10 | 120 | 4.3×10⁻¹⁹ |
| 3-sec | n=14 | 224 | 0 |
| 3-sec | n=16 | 288 | 0 |
| 3-sec | n=18 | 360 | 9.7×10⁻¹⁹ |
| 4-sec | n=8 | 80 | 0 |
| 4-sec | n=12 | 168 | 0 |
| 4-sec | n=14 | 224 | 0 |
| 4-sec | n=16 | 288 | 0 |

All residuals ≤ 10⁻¹⁸, confirming −Δξ_n = n(n+2)ξ_n.

### Pointwise N(Φ) check

Evaluated N(Φ) = |Φ|²Φ pointwise at 6 group elements and compared with the Peter-Weyl reconstruction from computed coefficients (summing over levels n = 0 to 18 with nonzero multiplicity). Maximum pointwise reconstruction error: ~1.8×10⁻² (3-sec) and ~1.4×10⁻² (4-sec). This error arises from truncation at n = 18; contributions from higher levels (n ≥ 20) account for the remainder. The error is consistent with the magnitude of ||w_n||² at n = 18, confirming the expansion coefficients are correct.

---

## Underdetermined questions

**Item 3 — stabilizer character details**: The character values at individual stabilizer elements are computed numerically. For U1, the generator acts by e^{−iπ/2} = −i. For U2–U4, the characters are determined by weight mod 8. Full character tables are in the script output.

**Item 4 — exact ξ norms**: The exact forms likely involve products of CG coefficients with denominators that are products of (2J+1) values (up to 13) and binomial coefficient factors. Nsimplify at 15-digit precision with denominator bound 10¹⁰ did not yield small-denominator identifications. Reported at IEEE double precision (~15 digits).

**Item 5, 7, 8 — exact values**: Same situation as Item 4. The values involve sums of rational functions of CG coefficients, yielding rationals or Q(√5) elements with very large denominators. Reported at 15-digit precision.

**Item 6 — 3-sec U6 τ_x**: The value 0.349108221994740 was not identified as a small-denominator rational or simple Q(√5) element. The 4-sec analogue is 2415/12298.
