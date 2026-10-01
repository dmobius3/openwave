# RETURN — solver_b

`R = 1` in every number below. All numbers come from the scripts in this room. `./py run_all.py` runs them in order: `t1_t3_symbolic.py` normally and with `--plant`, `t3_t5_exact.py` normally and with `--plant`, `t4_fem.py --checks` normally and with `--plant`, then `t4_fem.py` and `collect.py`. The full log of the last run is `run_all_log.txt`. `results.json` holds every number in the tables (block `"t4"`, plus `"t3_t5_exact"`, `"t4_extra_all_six"` and `"t4_minus_exact_extra"`, both labelled extra).

## 0. Reduction used throughout (stated once, used in T1–T5)

**Unfolding.** For `π/2 < y ≤ π`, put `y' = y − π ∈ (−π/2, 0]` and `w' = −w`. Then `−X(y, −w) = X(y', w')` (direct computation: `cos y = −cos y'`, `sin y = −sin y'`). So in `ℝP³` the second half of the rectangle is the antipodal image of the piece `X(y', w')`, `y' ∈ (−π/2, 0]`, of the same great sphere `x₄ = 0`. The whole of `M(W)` is therefore the image of the spherical **lune**
`L_W = {X(y', w) : −π/2 ≤ y' ≤ π/2, −W ≤ w ≤ W}` on `S²(1) = S³ ∩ {x₄ = 0}`,
with its two vertices `y' = ±π/2` (north and south poles, which are antipodal) identified to the cone point `p`. The seam `y = 0 ∼ y = π` becomes the interior equator arc `y' = 0` of the lune.

**Functions.** The antipodal map sends `e₄` to `−e₄`, so the normal field `φν` on the second half corresponds on the lune to `ψ e₄` with
`ψ(y', w') = φ(y, w)` for `y' = y ∈ [0, π/2]`, and `ψ(y', w') = −φ(y'+π, −w')` for `y' < 0`.
The seam condition `φ(π, −w) = −φ(0, w)` is exactly `ψ(0⁻, w) = ψ(0⁺, w)`, i.e. continuity of `ψ` across the equator. The map `U : φ ↦ ψ` is a unitary map of `L²(|cos y| dy dw)` onto `L²(L_W, cos y' dy' dw)`. It also preserves the Dirichlet integral, because the change of variables is an isometry on each half. The edges `w = ±W` go to the two sides of the lune. The "vanish near the cone point" condition becomes "vanish near both vertices".

**Remarks on the geometry, not used in the analysis.** (i) For `W = π/2` the lune is a hemisphere whose boundary great circle is identified antipodally in `ℝP³`, so the edge `w = W` and the edge `w = −W` are the same projective line in `ℝP³`. (ii) For `W = 17/10 > π/2` the lune overlaps its antipodal image in `ℝP³`. The spec says to treat both intrinsically ("as a surface in its own right"), and I do, with Dirichlet data on both sides.

---

## T1. Second variation of area

**Hypotheses.**
- (H1) `S³(1)` is round, and the projection `S³ → ℝP³` is a local isometry, so areas may be computed on the lift, piece by piece.
- (H2) The image of `X` lies in the great sphere `x₄ = 0`.
- (H3) `ν = e₄` at every point. It is a unit vector tangent to `S³` at every point of `x₄ = 0` (`X·e₄ = 0`) and normal to the surface (`∂X·e₄ = 0`).
- (H4) `φ` is smooth on the closed rectangle, satisfies the seam condition, vanishes on `w = ±W`, and vanishes on a neighbourhood of `y = π/2`.

**Step 1 (geodesics).** For `|x| = 1` and a unit vector `e ⊥ x`, the geodesic of `S³` is `γ(s) = cos s · x + sin s · e`. For general `R`: `cos(s/R) x + R sin(s/R) e`.
*Omitted:* without `e ⊥ x` (H3) this is not a curve on `S³`, and nothing below holds.

**Step 2 (the varied surface).** `F_t(y, w) = cos u · X + sin u · e₄` with `u = tφ(y, w)`, using the half-parametrizations of section 2.
Well-defined in `ℝP³`: at the seam, `F_t(π, −w) = cos(−u₀)(−X(0, w)) + sin(−u₀) e₄ = −F_t(0, w)` with `u₀ = tφ(0, w)`. So the two halves meet in `ℝP³`.
*Omitted:* without the seam condition the varied surface tears along the seam. The area integral below still makes sense, but it is no longer the area of a variation of `M(W)`.
The edge is fixed because `φ = 0` there. The cone point is fixed, and its neighbourhood untouched, because `φ = 0` near `y = π/2`.
*Omitted:* if `φ ≠ 0` on the collapsed segment, the single point `p` is sent to the curve `w ↦ F_t(π/2, w)`. The conic structure is then changed, against section 3.

**Step 3 (induced metric).** We have `∂ᵢF = cos u ∂ᵢX − sin u (∂ᵢu) X + cos u (∂ᵢu) e₄`. Using `X ⊥ ∂ᵢX` (from `|X|² = 1`), `X ⊥ e₄` and `∂ᵢX ⊥ e₄` (H2, H3), with `|X| = |e₄| = 1`:
`g_t = cos²u · g + du ⊗ du`.

**Step 4 (density).** By the matrix determinant lemma,
`det g_t = cos⁴u · det g · (1 + |∇u|²_g / cos²u) = cos²u · det g · (cos²u + t²|∇φ|²)`.
So the area density is `cos u · √(cos²u + t²|∇φ|²) · |cos y|` for `|t| < t₀`, where `t₀ max|φ| < π/2` (so that `cos u > 0`; H4 gives `φ` bounded). `t1_t3_symbolic.py` checks this determinant symbolically against a generic `φ(y, w)` (PASS line 1). This formula is exact, not an expansion.

**Step 5 (differentiate).** `A(t) = ∫_rect cos u √(cos²u + t²|∇φ|²) |cos y| dy dw`. The integrand is smooth in `t`, and its `t`-derivatives up to order 2 are bounded uniformly on the compact rectangle (H4: `φ` and `∇φ` are bounded; `|∇φ|² = φ_y² + φ_w²/cos²y` is bounded because `φ ≡ 0` near `cos y = 0`). So differentiation under the integral is allowed (dominated convergence).
*Omitted:* if `φ_w ≠ 0` near `y = π/2`, the term `φ_w²/cos²y` makes `|∇φ|²|cos y|` non-integrable in general, and the interchange fails.
Expanding: `d/dt|₀ = 0` (the surface is totally geodesic) and `d²/dt²|₀ = |∇φ|² − 2φ²`; both checked symbolically. Hence

> **`A''(0) = Q(φ) = ∫_M ( |∇φ|² − 2φ²/R² ) dA`**, with `dA = |cos(y/R)| dy dw`.

This is the general formula `∫ |∇φ|² − (Ric(ν, ν) + |A|²) φ²` with `Ric(ν, ν) = 2/R²` and `A = 0`.

**Boundary, seam and vertex terms: none.**
- Step 5 is a pointwise identity integrated over the rectangle. No integration by parts is done anywhere.
- So no edge term appears. In the general formula, edge terms come from tangential components of the variation field (none: it is purely normal), from the acceleration `∇_t∂_t F` (zero: the curves are geodesics) and from `H` (zero). Note honestly: the Dirichlet condition is *not needed for the formula itself*. It is what holds the edge fixed.
- No seam term appears, even though section 3's list allows `∂_yφ` to jump across the seam. The seam is a null set, and no derivative of `φ` across it is ever taken.
- No vertex term appears: `φ ≡ 0` near the cone point, so the integrand there equals its `t = 0` value.

**Extension by continuity.** Write `B(φ, χ) = ∫ ∇φ·∇χ − 2φχ`. Then `|B(φ, χ)| ≤ 2‖φ‖_{H¹} ‖χ‖_{H¹}`, where `‖φ‖²_{H¹} = ∫ |∇φ|² + φ²` is section 3's norm. Hence `|Q(φ) − Q(χ)| = |B(φ − χ, φ + χ)| ≤ 2‖φ − χ‖ ‖φ + χ‖`. So `Q` is uniformly continuous on bounded sets of the smooth class, and it extends uniquely to a continuous quadratic form on the completion, given by the same integral with weak gradients.
- The completion is a genuine space of functions: the norm dominates `L²` and the weak gradient is closable. If this were omitted, the "extension" would act on abstract Cauchy classes.
- What is extended is the *form*, not `A''(0)`. For a non-smooth admissible `φ` the varied surface need not be `C¹`, and `t ↦ A(t)` need not be twice differentiable.

**Assessment T1.** Established: the exact area formula, `Q` and its extension. The symbolic checks cover steps 3–5 for generic `φ`.

---

## T2. Self-adjointness, sector by sector

**Reading (underdetermined point).** Read literally, "the smooth functions of section 3's list" constrains only *values* at the seam (`φ(π, −w) = −φ(0, w)`), not `y`-derivatives. On that literal domain, `−Δ` applied pointwise on the rectangle is not even symmetric. Green's formula leaves the seam term `∫_seam (φ ∂ₙχ − χ ∂ₙφ)`, which is nonzero for functions with a kink at the seam, so asking about self-adjoint extensions makes no sense there.
I take the reading **D₀ = functions whose image `ψ = Uφ` is smooth on the lune** (all `y`-derivatives antiperiodic across the seam), vanishing on the sides and near both vertices. Under this reading `−Δ` on `D₀` is symmetric and non-negative: Green's formula has no seam term (the seam is interior and `ψ` is smooth), no side term (`ψ = 0`), and no vertex term (`ψ = 0` near the vertices). The form closure is unchanged by the reading, because kinked functions lie in the `H¹`-closure of `D₀`. So the Friedrichs operator, and hence T3–T5, are the same under either reading.

**Sectors.** Take `s_k(w) = W^{−1/2} sin(kπ(w+W)/(2W))`, `k ≥ 1`. These form an orthonormal basis of `L²(−W, W)`, with `−s_k'' = ν_k² s_k`, `ν_k = kπ/(2W)`. The area element is `cos y' dy' dw`, a product, so
`L²(L_W) = ⊕_k L²((−π/2, π/2), cos y' dy') ⊗ s_k`.
In sector `k`, `−Δ` acts as `l_k f = −(cos y')⁻¹ (cos y' f')' + ν_k² f / cos²y'`.

**Deficiency indices split over sectors.** Let `A` be `−Δ` on `D₀`.
- (a) For `f ∈ C_c^∞(−π/2, π/2)`, `f s_k ∈ D₀` and `A(f s_k) = (l_k f) s_k`. Testing `A*u` against these gives `(A*u)_k = l_k u_k` in the distributional sense, so `u_k ∈ D(l_k,max)`.
- (b) Conversely, for `v ∈ D₀` the coefficient `v_k(y') = ∫ v s_k dw` is smooth with compact support in `(−π/2, π/2)`. It has compact support because `v` vanishes near the vertices; without that hypothesis `v_k` would not be compactly supported and the step fails. Also `(Av)_k = l_k v_k`, by two integrations by parts in `w` whose boundary terms vanish since `v = s_k = 0` at `w = ±W`.
- (c) Hence `D(A*) = {u : u_k ∈ D(l_k,max), Σ‖l_k u_k‖² < ∞}`, and `ker(A* ∓ i) = ⊕_k ker(l_k,max ∓ i)`. So `n± = Σ_k n±(k)`, and the closure of `A` is self-adjoint iff every sector has `n±(k) = 0`.

**Each sector at the cone point.** Interior points are regular. Both ends `y' = ±π/2` are singular: the weight `cos y' → 0` and the potential blows up. The two ends are symmetric under `y' ↦ −y'`. Let `r` be the distance to a vertex, so `cos y' = sin r`.
- At `λ = 0` the equation has the exact solutions `f± = tan^{±ν}(r/2)`, checked symbolically (PASS). The indicial exponents are `±ν` (also checked).
- Near `r = 0`, `f±² sin r ~ 2^{∓2ν} r^{1 ± 2ν}`. So `f₊ ∈ L²` always, and `f₋ ∈ L²(sin r dr)` near 0 **iff `ν < 1`**.
- By Weyl's alternative (for one `λ`, all solutions are `L²` iff the endpoint is limit-circle), each end is **limit-circle iff `ν_k < 1`** and **limit-point iff `ν_k ≥ 1`**.
- The borderline `ν = 1` is limit-point: `∫_ε r⁻¹ dr = log(1/ε)` diverges (PASS line). Writing "`ν ≤ 1`" instead of "`ν < 1`" would misclassify `W = π/2`.
- In the limit-circle case both ends are limit-circle, so `n±(k) = 2`. In the limit-point case `n±(k) = 0`.

**Classification.** `ν_k < 1 ⇔ k < 2W/π`.
- **`0 < W ≤ π/2`:** `ν_k ≥ ν₁ = π/(2W) ≥ 1` for all `k`. **Every sector is limit-point at the cone point**; `W = π/2`, `k = 1` is the borderline case `ν = 1`, still limit-point. So `n± = 0`: the closure is self-adjoint, and it is the only self-adjoint extension.
- **`W = 17/10`:** `ν₁ = π/3.4 = 0.92399…`. This is `< 1` because `π < 3.4`, an exact inequality. So **sector `k = 1` is limit-circle at both ends** (both are the cone point). For `k ≥ 2`, `ν_k ≥ 2π/3.4 > 1`, so those sectors are limit-point. In total `n± = (2, 2)`: the closure is *not* self-adjoint. Its self-adjoint extensions form a `U(2)` family of boundary conditions on the coefficients of `(f₊, f₋)` at the two ends. These include conditions that couple the two ends, which is geometrically natural since both ends are the single point `p`.

**Which realization the admissible class selects.** On `D₀`, `⟨(A+1)φ, φ⟩ = ∫ |∇φ|² + φ²`, which is section 3's norm. The admissible class contains `D₀` and lies in the `H¹`-closure of `D₀`, so it equals the form domain of the Friedrichs extension `A_F`. In sector 1 the form `∫ (|f'|² + ν² f² / sin² r) sin r dr` is infinite for a cutoff of `f₋` (`∫ r^{−2ν−1} dr = ∞`) and finite for a cutoff of `f₊` (`∫ r^{2ν−1} dr < ∞`, `ν > 0`). So the admissible class selects **the Friedrichs extension**: in the limit-circle sector it imposes, separately at each end, that the coefficient of `tan^{−ν}(r/2)` vanishes (the principal / regular solution). This is the decoupled condition. Among all self-adjoint extensions it is the only one whose operator domain lies in the admissible class.

**Assessment T2.** Established under the stated reading. The classification rests on Weyl's alternative and standard deficiency theory, both used from memory (see manifest).

---

## T3. Spectrum of `−Δ` on the admissible class (= `A_F`)

**Claim.** The spectrum is exactly
> **`λ_{k,m} = (ν_k + m)(ν_k + m + 1) / R²`, with `ν_k = kπR/(2W)`, `k = 1, 2, …`, `m = 0, 1, 2, …`**,

with eigenfunctions (in lune coordinates, colatitude `θ = π/2 − y'`)
`F_{k,m} = sin^{ν_k}θ · C_m^{(ν_k + 1/2)}(cos θ) · s_k(w)`.
Each eigenvalue has multiplicity equal to the number of labels `(k, m)` giving it. The **lowest is `λ_{1,0} = (π/(2W)) (π/(2W) + 1)`**, eigenfunction `cos^{π/(2W)}y' · cos(πw/(2W))`. In rectangle coordinates this is `±|cos y|^{π/(2W)} cos(πw/(2W))`, with `+` on `y < π/2` and `−` on `y > π/2`. The formula holds unchanged for `W = 17/10`.

**Derivation.**
1. *Sector ODE.* In colatitude, `l_k f = −(sin θ)⁻¹ (sin θ f')' + ν² f / sin²θ`. Substituting `f = sin^ν θ · g(x)`, `x = cos θ`, gives the Gegenbauer equation `(1−x²)g'' − (2ν+2) x g' + (λ − ν(ν+1)) g = 0`. This has the polynomial solution `C_m^{(ν+1/2)}` iff `λ = ν(ν+1) + m(m+2ν+1) = (ν+m)(ν+m+1)`. Checked symbolically for `m = 0..5` with symbolic `ν` (PASS lines).
2. *Each `F_{k,m}` is in the form domain.*
   - `F` is smooth in the interior of the lune, including across the seam, and vanishes on the sides.
   - Near a vertex, `|F| ≲ r^ν` and `|∇F| ≲ r^{ν−1}`, so `∫ |∇F|² dA < ∞` because `ν > 0`.
   - The cutoffs `χ(r/ε) F ∈ D₀` converge to `F` in `H¹`: the error is `∫ |∇χ_ε|² F² ≲ ε⁻² ε^{2ν} ε² → 0`.
   - *Omitted:* for `ν ≤ 0` this would fail. Here `ν > 0` because the `w`-data are Dirichlet.
3. *Each `F_{k,m}` is in `D(A_F)` with `A_F F = λF`.* For `v ∈ D₀`, Green's formula gives `∫ ∇F·∇v = λ ∫ F v`. There is no side term (`v = 0` on the sides), no vertex term (`v = 0` near the vertices) and no seam term (`F` smooth). By density this extends to all `v` in the form domain, which is the defining property of `D(A_F)`.
4. *Completeness.* `{s_k}` is an orthonormal basis in `w`. For fixed `k`, the map `f ↦ g = f / sin^ν θ` is an isometry `L²(sin θ dθ) → L²((1−x²)^ν dx)`. In that space `C_m^{(ν+1/2)}` (degree exactly `m`, orthogonal for the weight `(1−x²)^ν`) span the polynomials, which are dense by Weierstrass on a finite measure on `[−1, 1]`. Orthogonality is checked symbolically for sample pairs (PASS lines). So `{F_{k,m}}` is a complete orthogonal set in `L²(L_W)`.
5. *Spectrum.* A self-adjoint operator with a complete orthonormal set of eigenvectors has spectrum equal to the closure of the eigenvalues. Since `λ_{k,m} ≥ max(ν_k², m²) → ∞`, the set is discrete and closed, so the spectrum is exactly `{λ_{k,m}}`. *Omitted:* without completeness (step 4) there could be further spectrum.

**Values** (`t3_t5_exact.py`, exact expressions in `exact.json` and `results.json`). The values were evaluated at 30 and at 50 digits and agree to `10⁻²⁸`. They are evaluations of closed forms, not identifications.

| W | ν₁ | lowest six λ (label k,m) |
|---|---|---|
| 1/4 | 2π | 45.7616029115 (1,0), 60.3279735259 (1,1), 76.8943441403 (1,2), 95.4607147546 (1,3), 116.0270853690 (1,4), 138.5934559833 (1,5) |
| 1/2 | π | 13.0111970547 (1,0), 21.2943823619 (1,1), 31.5775676690 (1,2), 43.8607529762 (1,3), 45.7616029115 (2,0), 58.1439382834 (1,4) |
| 1 | π/2 | 4.0381974271 (1,0), 9.1797900807 (1,1), 13.0111970547 (2,0), 16.3213827342 (1,2), 21.2943823619 (2,1), 25.4629753878 (1,3) |
| 7/5 | 5π/14 | 2.3808754887 (1,0), 6.6248702412 (1,1), 7.2795072021 (2,0), 12.8688649938 (1,2), 13.7674967072 (2,1), 14.6958951403 (3,0) |
| 3/2 | π/3 | 2.1438202624 (1,0), 6.2382153648 (1,1), 6.4808859473 (2,0), 12.3326104672 (1,2), 12.6696761521 (2,1), 13.0111970547 (3,0) |
| π/2 | 1 | 2 (1,0), 6 (1,1), 6 (2,0), 12 (1,2), 12 (2,1), 12 (3,0) |
| 17/10 | 5π/17 | 1.7777698463 (1,0), 5.2630837066 (2,0), 5.6257655249 (1,1), 10.4559415810 (3,0), 10.9590750638 (2,1), 11.4737612035 (1,2) |

**Assessment T3.** Established, by a symbolic derivation plus a completeness argument. The Gegenbauer orthogonality/weight fact is from memory, spot-checked symbolically.

---

## T4. Finite elements (protocol as specified)

**Implementation (`t4_fem.py`).**
- P1 elements on the rectangle `[0, π] × [−W, W]`.
- Stiffness `∫ (|cos y| φ_y² + φ_w²/|cos y|)`, mass `∫ |cos y| φ²`.
- The node `(π, −w_j)` is the dof of `(0, w_j)` with sign `−1`.
- Zero data on `w = ±W` and on the fiber `y = π/2`.
- `y`-nodes at `π/2 ∓ (π/2)(k/16)²`, so 32 cells in `y`; 8 uniform cells in `w`. Four uniform refinements, giving levels of 32×8, 64×16, 128×32, 256×64 and 512×128 cells (217 → 64 897 dofs).
- Lowest six eigenvalues by shift-invert Lanczos at `σ = 0`.
- `p = log₂(|λ₂ − λ₃| / |λ₃ − λ₄|)` from levels 2, 3, 4 (0-based) for the bottom eigenvalue. `λ_extrap = λ₄ + (λ₄ − λ₃)/(2^p − 1)`, and `err = |λ_extrap − λ₄|`.

**Readings for underdetermined points of the protocol.**
- (a) "Uniform refinement" means bisecting every cell in `y` and in `w` at its `(y, w)`-midpoint. This keeps the initial grading; it does not regenerate `(k/2N)²` nodes.
- (b) Each rectangle cell is split into two triangles by its `(y₀, w₀)–(y₁, w₁)` diagonal.
- (c) Quadrature is 12-point Gauss in `y` × 3-point Gauss in `w` on each triangle. The stiffness needs only the weight integrals, since gradients are constant per triangle. On triangles with an edge on the fiber, the only free node has an exactly zero `w`-gradient, so the large `1/|cos y|` integral there multiplies 0. On triangles touching the fiber only at a vertex, the integrand `width(y)/|cos y|` is bounded.

**Checks (can fail; planted defects shown firing in `run_all_log.txt`).**
- Weighted area equals `4W`. Planted `cos y` in place of `|cos y|` gives area 0 → FAIL.
- The seam dof map reproduces an antiperiodic test function at `y = π`. Planted periodic seam gives error 2.0 → FAIL.
- `K`, `M` are symmetric and `M` is positive definite. The planted `cos y` weight makes `M` indefinite → FAIL.
- Also, the first run of `t1_t3_symbolic.py` printed a real (unplanted) FAIL on the indicial check. I had declared the exponent symbol positive, so sympy dropped `−ν`. I fixed the symbol declaration, not the expectation. This is further evidence that the line can fail.

**Tables** (lowest six eigenvalues per level; the `exact` row is from T3 and is **extra**, shown for comparison only).

### W = 1/4
| level | λ1 | λ2 | λ3 | λ4 | λ5 | λ6 |
|---|---|---|---|---|---|---|
| 0 | 47.1476815893 | 64.3714547952 | 84.6915018190 | 108.6934750033 | 136.4154980252 | 168.6695606584 |
| 1 | 46.1158097873 | 61.3603140582 | 78.8877413710 | 98.8254717310 | 121.2388143889 | 146.2524103195 |
| 2 | 45.8507607068 | 60.5876239882 | 77.3958886140 | 96.3059030453 | 117.3368920021 | 140.5159219711 |
| 3 | 45.7839334422 | 60.3929885130 | 77.0199405617 | 95.6722742148 | 116.3549771052 | 139.0745629759 |
| 4 | 45.7671881775 | 60.3442337480 | 76.9257565912 | 95.5136212459 | 116.1090860824 | 138.7137636496 |
| exact (extra) | 45.7616029115 | 60.3279735259 | 76.8943441403 | 95.4607147546 | 116.0270853690 | 138.5934559833 |

p = 1.996684, extrap = 45.761589281931, err = 5.599e-03 (extra: extrap − exact = −1.363e-05)

### W = 1/2
| level | λ1 | λ2 | λ3 | λ4 | λ5 | λ6 |
|---|---|---|---|---|---|---|
| 0 | 13.3126078573 | 22.2194759002 | 33.5685384037 | 47.5942447618 | 49.0137421855 | 64.5646145956 |
| 1 | 13.0868167765 | 21.5264939056 | 32.0818941344 | 44.8064220477 | 46.5844116303 | 59.7650030564 |
| 2 | 13.0301243623 | 21.3524704626 | 31.7040809952 | 44.0979832459 | 45.9683874194 | 58.5501859128 |
| 3 | 13.0159304072 | 21.3089083119 | 31.6092235266 | 43.9201125170 | 45.8133785900 | 58.2455624934 |
| 4 | 13.0123804907 | 21.2980140974 | 31.5854833648 | 43.8755961305 | 45.7745520368 | 58.1693482552 |
| exact (extra) | 13.0111970547 | 21.2943823619 | 31.5775676690 | 43.8607529762 | 45.7616029115 | 58.1439382834 |

p = 1.999420, extrap = 13.011196550217, err = 1.184e-03 (extra: extrap − exact = −5.045e-07)

### W = 1
| level | λ1 | λ2 | λ3 | λ4 | λ5 | λ6 |
|---|---|---|---|---|---|---|
| 0 | 4.1201239976 | 9.4810950163 | 13.8719077087 | 17.1000342615 | 23.2813603638 | 27.0897906115 |
| 1 | 4.0586293941 | 9.2550724839 | 13.2272330236 | 16.5131874034 | 21.7956917617 | 25.8624441097 |
| 2 | 4.0433027231 | 9.1986088226 | 13.0652528523 | 16.3691893425 | 21.4199131927 | 25.5625075009 |
| 3 | 4.0394735914 | 9.1844946696 | 13.0247140008 | 16.3333259491 | 21.3257768571 | 25.4878394837 |
| 4 | 4.0385164584 | 9.1809662222 | 13.0145764817 | 16.3243680206 | 21.3022317097 | 25.4691902608 |
| exact (extra) | 4.0381974271 | 9.1797900807 | 13.0111970547 | 16.3213827342 | 21.2943823619 | 25.4629753878 |

p = 2.000226, extrap = 4.038197480573, err = 3.190e-04 (extra: extrap − exact = 5.351e-08)

### W = 7/5
| level | λ1 | λ2 | λ3 | λ4 | λ5 | λ6 |
|---|---|---|---|---|---|---|
| 0 | 2.4282858380 | 6.8267585220 | 7.7650797853 | 13.4222387659 | 15.0293106368 | 16.7663490529 |
| 1 | 2.3926822162 | 6.6752602421 | 7.3997509968 | 13.0064837665 | 14.0806251808 | 15.2071295106 |
| 2 | 2.3838245317 | 6.6374635061 | 7.3094955941 | 12.9032375064 | 13.8456190717 | 14.8231821257 |
| 3 | 2.3816126002 | 6.6280183693 | 7.2869998074 | 12.8774565389 | 13.7870172098 | 14.7276833601 |
| 4 | 2.3810597597 | 6.6256572743 | 7.2813800739 | 12.8710128223 | 13.7723762014 | 14.7038400940 |
| exact (extra) | 2.3808754887 | 6.6248702412 | 7.2795072021 | 12.8688649938 | 13.7674967072 | 14.6958951403 |

p = 2.000371, extrap = 2.380875542763, err = 1.842e-04 (extra: extrap − exact = 5.410e-08)

### W = 3/2
| level | λ1 | λ2 | λ3 | λ4 | λ5 | λ6 |
|---|---|---|---|---|---|---|
| 0 | 2.1865024844 | 6.4263245361 | 6.9140518166 | 12.8549195651 | 13.8285296756 | 14.8551467853 |
| 1 | 2.1544479192 | 6.2851612727 | 6.5879966135 | 12.4625954396 | 12.9567352005 | 13.4647893694 |
| 2 | 2.1464748164 | 6.2499481908 | 6.5075889694 | 12.3650833310 | 12.7412629814 | 13.1240244452 |
| 3 | 2.1444837946 | 6.2411485492 | 6.4875570675 | 12.3407279985 | 12.6875616560 | 13.0393673763 |
| 4 | 2.1439861476 | 6.2389487041 | 6.4825534386 | 12.3346399375 | 12.6741468272 | 13.0182373457 |
| exact (extra) | 2.1438202624 | 6.2382153648 | 6.4808859473 | 12.3326104672 | 12.6696761521 | 13.0111970547 |

p = 2.000314, extrap = 2.143820313410, err = 1.658e-04 (extra: extrap − exact = 5.098e-08)

### W = π/2
| level | λ1 | λ2 | λ3 | λ4 | λ5 | λ6 |
|---|---|---|---|---|---|---|
| 0 | 2.0398446446 | 6.1798114043 | 6.4016403433 | 12.5035338084 | 13.0966272424 | 13.7071188919 |
| 1 | 2.0099207147 | 6.0448742180 | 6.0992291697 | 12.1253672917 | 12.2713472267 | 12.4190106864 |
| 2 | 2.0024781130 | 6.0112158016 | 6.0247325760 | 12.0313247125 | 12.0676504257 | 12.1041670094 |
| 3 | 2.0006194788 | 6.0028041737 | 6.0061784976 | 12.0078313970 | 12.0169008856 | 12.0260043340 |
| 4 | 2.0001548856 | 6.0007011524 | 6.0015443348 | 12.0019581322 | 12.0042244883 | 12.0064987382 |
| exact (extra) | 2 | 6 | 6 | 12 | 12 | 12 |

p = 2.000203, extrap = 2.000000050300, err = 1.548e-04 (extra: extrap − exact = 5.030e-08). The discretization splits the exact degeneracies (6 twice, 12 three times).

### W = 17/10
| level | λ1 | λ2 | λ3 | λ4 | λ5 | λ6 |
|---|---|---|---|---|---|---|
| 0 | 1.8132823372 | 5.6164590398 | 5.7928633358 | 11.9463898560 | 11.9549594879 | 11.9599043885 |
| 1 | 1.7866129810 | 5.3502722352 | 5.6674720810 | 10.8221677892 | 11.2062932234 | 11.5920255262 |
| 2 | 1.7799794312 | 5.2848076529 | 5.6361926219 | 10.5469091989 | 11.0206837269 | 11.5033239002 |
| 3 | 1.7783224013 | 5.2685101338 | 5.6283734585 | 10.4786459308 | 10.9744649295 | 11.4811549184 |
| 4 | 1.7779080592 | 5.2644400290 | 5.6264178922 | 10.4616153144 | 10.9629217611 | 11.4756106715 |
| exact (extra) | 1.7777698463 | 5.2630837066 | 5.6257655249 | 10.4559415810 | 10.9590750638 | 11.4737612035 |

p = 1.999706, extrap = 1.777769907611, err = 1.382e-04 (extra: extrap − exact = 6.131e-08)

**Extra (labelled):**
- `results.json → "t4_extra_all_six"` gives `p`, extrapolation and error for all six eigenvalues by the same formulas.
- `"t4_minus_exact_extra"` gives finest-level minus exact.
- Nothing was tuned. All the levels converge from above, as a conforming method must.
- The protocol's `err` (= `|extrap − finest|`) overstates the actual error of the extrapolated value by 2–4 orders of magnitude at every width. It is an error estimate for the finest level, not for `λ_extrap`.

**Assessment T4.** All seven widths were completed under the stated readings, and a full `run_all.py` rerun reproduced every number. The observed order is `p ≈ 2` at all widths, including `W = 17/10` where `ν₁ < 1`; the quadratic grading evidently suffices there. The FEM bottoms agree with T3 to `≤ 1.4·10⁻⁵` after extrapolation. This is consistent with, but not a proof of, T3.

---

## T5. Sign of the lowest eigenvalue, index, nullity

**Argument.**
- The second-variation form `Q` (T1) has the same form domain as `A_F` (the admissible class, T2) and differs from `A_F`'s form by `−2‖φ‖²`. So its operator is the Jacobi operator `J = A_F − 2`.
- By T3, `J` has the eigenvalues `λ_{k,m} − 2 = (ν_k + m − 1)(ν_k + m + 2)`, discrete with finite multiplicities.
- Min–max: the index of `Q` on the admissible class (the maximal dimension of a subspace on which `Q < 0`) equals the number of negative eigenvalues of `J` with multiplicity. The nullity (`{u : B(u, v) = 0 ∀v}`) equals `dim ker J`.
- The second factor is `> 0`, so the sign of each eigenvalue is the sign of `ν_k + m − 1`. A negative eigenvalue needs `m = 0` and `ν_k < 1`, i.e. `k < 2W/π`. A zero needs `m = 0` and `k = 2W/π`; the other option, `m = 1` with `ν_k = 0`, is impossible.
- Lowest: `λ₁(J) = (π/(2W) − 1)(π/(2W) + 2)`.
- *Omitted:* without T2's identification of the admissible class with the Friedrichs form domain, `J` could be another extension, and in the limit-circle case `W = 17/10` its spectrum, and possibly its index, would differ.

The exact sign decisions use `3 < π < 3.4`; sympy decides them on exact expressions.

| W | λ₁(J) exact | λ₁(J) | sign | index | nullity |
|---|---|---|---|---|---|
| 1/4 | 4π² + 2π − 2 | 43.761602911537 | + | 0 | 0 |
| 1/2 | π² + π − 2 | 11.011197054679 | + | 0 | 0 |
| 1 | (π − 2)(π + 4)/4 | 2.038197427067 | + | 0 | 0 |
| 7/5 | 25π²/196 + 5π/14 − 2 | 0.380875488666 | + | 0 | 0 |
| 3/2 | (π − 3)(π + 6)/9 | 0.143820262429 | + | 0 | 0 |
| π/2 | 0 | 0 | zero | 0 | 1 |
| 17/10 | 25π²/289 + 5π/17 − 2 | −0.222230153694 | − | 1 | 0 |

- At `W = π/2` the null Jacobi field is `cos y' cos w`, the restriction of `x₁`. This is the normal component of the rotation in the `x₁x₄`-plane.
- At `W = 17/10` the single negative direction is the `(1, 0)` mode, `ν₁ = 5π/17`.
- The FEM numbers (T4) agree: the extrapolated `λ₁ − 2` has the same sign at every width. At `W = π/2` it is `+5·10⁻⁸`, which is zero within the method's accuracy; the sign there comes from T3, not from the numerics.

**Assessment T5.** Established, given T2 and T3.

---

## Overall assessment and open points
- **Established:** T1, T3 and T5, by derivations with symbolic checks; T2 under the stated reading of D₀; T4 under readings (a)–(c), carried out and reproduced.
- **Not settled / caveats:**
  - The literal domain in T2 is not symmetric, hence the reading.
  - The two geometric remarks in §0 (edges coincide in `ℝP³` at `W = π/2`; self-overlap at `W = 17/10`) were treated intrinsically, as instructed.
  - Weyl's alternative, deficiency theory, Friedrichs theory and min–max are used from memory, not re-derived.
- **Disagreement with my own steps:** the only one was the indicial-exponent check's first FAIL (a symbol-assumption bug), described in T4's checks.

## Consulted-material manifest
- **Files read:** `BRIEF.md`, `spec_sheet.md`, and my own outputs `t4_log.txt` and `run_all_log.txt`.
- **Directory listing:** one listing of the room showed the names `.room.sb` and `py`. I did not open either.
- **Loaded without asking:** a session context note (the user's account e-mail), not used. A sandbox error message from the first `run_all.py` attempt printed the absolute path of `.room.sb`.
- **Outside the room:** no reads. Harness task-output files exist outside the room; I did not read them.
- **No network or search.**
- **Results used from memory, not derived here:**
  - the geodesic formula on the round sphere;
  - the matrix determinant lemma;
  - `Ric = 2/R²` on `S³(R)`, and the general second-variation formula (used only for comparison);
  - Weyl's limit-point/limit-circle alternative, and deficiency indices for Sturm–Liouville operators with limit-circle ends;
  - the Friedrichs extension and its characterization by the form domain;
  - the Frobenius method;
  - the zero-energy solutions `tan^{±ν}(r/2)` (recalled, then verified symbolically);
  - the Gegenbauer equation and orthogonality weight `(1−x²)^{α−1/2}` (recalled, then verified symbolically for samples);
  - the Weierstrass density of polynomials;
  - the cutoff / capacity argument in 2D;
  - the min–max principle;
  - the Richardson extrapolation formula.
