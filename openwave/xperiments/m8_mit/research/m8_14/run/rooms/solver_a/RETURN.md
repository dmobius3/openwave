# RETURN — solver_a

All numbers use `R = 1`; the derivations keep `R` where it costs nothing. Everything numeric comes from `run_all.py`, run as `./py run_all.py`. It runs `t1_second_variation.py`, `t3_spectrum.py`, `t4_fem.py`, the planted-defect runs, and `make_tables.py`, and writes `results.json`. The log of the last run is `run_all.log`.

## 0. Readings taken (where the sheet is underdetermined)

- **R1 (seam smoothness, T2).** In T2, "smooth functions" means smooth *on the surface*. So the seam matching holds for all `y`-derivatives: `∂_y^m φ(πR,−w) = −∂_y^m φ(0,w)`. The seam is an interior curve of the surface, not part of the edge. Under the alternative reading, where only values match, the operator is not even symmetric, because the boundary terms `∫(φ_y ψ − φ ψ_y)` at `y = 0, πR` do not cancel. So no self-adjointness question can be asked. The form domain (the `H¹` completion) is the same under both readings, because the completion only sees the trace.
- **R2 (refinement, T4).** "Uniform refinement" means midpoint (red) refinement of every triangle. The graded nodes of level 0 are kept and midpoints are inserted. Level ℓ therefore has `16·2^ℓ` intervals per side, but these are not regraded with `(k/N)²` at the larger `N`.
- **R3 (triangulation, T4).** Each `(y,w)` cell is cut along the diagonal `(y0,w0)–(y1,w1)`. Red refinement reproduces this structured mesh (checked by hand: the 4 children of each triangle are again cells of the refined grid with the same diagonal).
- **R4 (estimates, T4).** `p = log2((λ⁽²⁾−λ⁽³⁾)/(λ⁽³⁾−λ⁽⁴⁾))` from the bottom eigenvalue on levels 2, 3, 4 (0-based), since `h` halves per level. `λ_extrap = λ⁽⁴⁾ + (λ⁽⁴⁾−λ⁽³⁾)/(2^p − 1)` uses that estimated `p`. `err = |λ_extrap − λ⁽⁴⁾|`.
- **R5 (quadrature, T4).** The weights are integrated with 16-point Gauss–Legendre in `y` and exactly (Simpson, exact for quadratics) in `w`. On every element the integrands are analytic: `|cos y| = sin s` with `s = |y − π/2|`, and the only singular weight `1/|cos y|` appears multiplied by a width that vanishes linearly at the fiber. See T4 for the one exception, where the free-node `w`-gradient is exactly zero.
- **R6 (T1 extension).** For non-smooth admissible `φ` the area of the varied surface is not defined pointwise. So "second variation" on the whole class *means* the unique continuous extension of the form.
- **R7 (W/R = 17/10).** The surface is only immersed in `ℝP³` (the unfolded lune is wider than a hemisphere). As the sheet says, I treat it as the abstract surface with metric `ds²` on the rectangle and the same data.

## 1. Geometry used throughout: the unfolding

**Hypotheses.** (G1) `X` from section 2. (G2) The seam identification `(0,w) ∼ (πR,−w)` with the sign rule `φ(πR,−w) = −φ(0,w)` (section 3). (G3) The metric `ds² = dy² + cos²(y/R) dw²`.

Define latitude/longitude `(ψ, θ)` on the rectangle minus the fiber:
- for `0 ≤ y < πR/2`, set `(ψ,θ) = (y/R, w/R)` and `u = φ`;
- for `πR/2 < y ≤ πR`, set `(ψ,θ) = (y/R − π, −w/R)` and `u = −φ`.

Then:
1. `ds² = R²(dψ² + cos²ψ dθ²)`, because `cos²(y/R − π) = cos²(y/R)`. *Without (G3)'s form this isometry fails.*
2. The line `y = 0` goes to `ψ = 0⁺`, and the line `y = πR` (at `w`) goes to `ψ = 0⁻` (at `θ = −w/R`). So the seam becomes the interior line `ψ = 0`. The condition `φ(πR,−w) = −φ(0,w)` says exactly that `u` is continuous across `ψ = 0`. Under R1, `u` is smooth across it. *Without the minus sign in the definition of `u` on the lower sheet, the seam condition would turn into a jump of `u` across `ψ = 0`. This is planted defect C3.3 `seam`, which fires.*
3. The two sides of the fiber become the two poles `ψ = ±π/2`. Both map to the cone point `p`. The edge lines `w = ±W` become `θ = ±W/R`.
4. Integrands are preserved: `φ² = u²`, `|∇φ|² = |∇u|²`, and `dA = R² cos ψ dψ dθ`.

So every problem below is a problem on the "lune" `L = (−π/2, π/2) × [−W/R, W/R]` with metric `R²(dψ² + cos²ψ dθ²)`. It has Dirichlet data on `θ = ±W/R`, and the two poles are excluded (both are the cone point). For `W ≤ πR/2` this is an honest lune in the great 2-sphere `x₄ = 0`, between longitudes `±W/R`.

## 2. T1 — the second variation

**Hypotheses.**
- (H1) `S³(R)` is round, and `M` lies in the great sphere `Σ = {x₄ = 0}`.
- (H2) `ν = e₄` is a unit vector orthogonal to `X`, `X_y` and `X_w` at every point.
- (H3) `φ` is smooth on the rectangle, vanishes on `w = ±W`, vanishes on a neighbourhood `{|y − πR/2| < δ}` of the fiber, and satisfies `φ(πR,−w) = −φ(0,w)`.
- (H4) The variation follows normal geodesics: `exp_x(v) = cos(|v|/R) x + R sin(|v|/R) v/|v|` for `v ⟂ x`.

**Step 1 (the varied map).** By (H2) and (H4), `F_t = cos(tφ/R) X + R sin(tφ/R) e₄` on each branch. Write `c = cos(tφ/R)`. *If (H4) were replaced by the straight line `X + tφe₄`, the map would leave `S³` and step 3 would give `∫|∇φ|²` with no `−2φ²/R²`. This is planted defect `straight`, and C1.1 and C1.2 fire.*

**Step 2 (well defined on ℝP³; edge and cone point fixed).** At the seam, branch 2 at `(πR, −w)` is `X(πR, w) = −X(0, w)`. Then `F_t(πR,−w) = c(−φ(0,w))·(−X(0,w)) + R sin(−tφ(0,w)/R) e₄ = −F_t(0,w)`. This is the same point of `ℝP³` (symbolic check C1.3). *Without the sign rule `φ(πR,−w) = −φ(0,w)`, the two seam points would separate and the varied surface would tear. This is planted defect `seamsign`, and C1.3 fires.* On the edge `φ = 0`, so `F_t = X` and the edge is fixed pointwise. On `{|y − πR/2| < δ}` again `F_t = X`, so the cone point and its conic structure are untouched. *Without the "vanish near the cone point" hypothesis, the cone point would move, or the varied surface would be singular differently there.*

**Step 3 (exact metric).** Differentiate:
- `F_y = c X_y + t φ_y (−sin(tφ/R) X + cos(tφ/R) R e₄)/R`, and similarly for `F_w`.
- Use `X·X = R²` (hence `X·X_y = X·X_w = 0`), `X·e₄ = X_y·e₄ = X_w·e₄ = 0` and `e₄·e₄ = 1`.
- The cross terms vanish and `(−sin X + R cos e₄)·(same)/R² = 1`.

So exactly, `g_t = c² g + t² dφ⊗dφ` with `g = dy² + cos²(y/R) dw²`. Check C1.1 is a symbolic identity on both branches. *If `ν` were not orthogonal to `X` and to `TM` (H2), cross terms `t·(X_y·ν)` would survive and a first-order term would appear. If `|X| = R` failed, `X·X_y ≠ 0` would add terms.*

**Step 4 (exact area).** By the matrix-determinant lemma, `det g_t = c⁴ det g (1 + t²|∇φ|²_g/c²)`. So

`A(t) = ∫ √(c⁴ + t² c² |∇φ|²_g) dA`, with `|∇φ|²_g = φ_y² + φ_w²/cos²(y/R)` and `dA = |cos(y/R)| dy dw`.

This holds for `|t| < πR/(2 max|φ|)`, where `c > 0`.

**Step 5 (differentiate under the integral).** By (H3), `|∇φ|²_g` is bounded on the rectangle: `φ_w = 0` where `cos(y/R)` is small. So the integrand and its `t`-derivatives are bounded uniformly, and differentiation under the integral is legitimate. *Without "vanish near the cone point", `φ_w²/cos²` can be non-integrable against `|cos| dy` (for example, `φ_w ≠ 0` on the fiber gives `∫ dy/|cos y| = ∞`).* Expand:
- `c⁴ = 1 − 2t²φ²/R² + O(t⁴)` and `c²t²|∇φ|² = t²|∇φ|² + O(t⁴)`;
- `√(1+ε) = 1 + ε/2 + O(ε²)`.

Hence `A(0) = Area(M)`, `A'(0) = 0`, and

**`A''(0) = Q(φ) = ∫_M (|∇φ|² − (2/R²) φ²) dA = ∫∫ [φ_y² + φ_w²/cos²(y/R) − 2φ²/R²] |cos(y/R)| dy dw.`**

C1.2 checks the `t⁰`, `t¹` and `t²` coefficients symbolically on both branches.

**Boundary or vertex term: none.** The derivation is pointwise: no integration by parts is ever done, so nothing is produced on the edge, at the seam, at the edge corners, or at the cone point. In the general formula `∫(|∇φ|² − (|A|² + Ric(ν,ν))φ²) + boundary terms`, the terms are `|A|² = 0` (totally geodesic, H1) and `Ric(ν,ν) = 2/R²`. The boundary term involves the acceleration `∇_t ∂_t F`, which is zero for geodesic variations (H4). It is also multiplied by `φ = 0` on the edge. *If (H1) failed, `|A|²` would appear. If (H4) failed, an acceleration term `∫ ⟨H, ∇_t∂_tF⟩` and an edge term would appear.*

**Extension to the admissible class.** Let `‖φ‖² = ∫(|∇φ|² + φ²) dA`. Then `Q(φ) = ‖φ‖² − (1 + 2/R²)‖φ‖²_{L²}` and `|Q(φ) − Q(χ)| ≤ (2 + 2/R²)‖φ − χ‖ (‖φ‖ + ‖χ‖)`. So `Q` is uniformly continuous on bounded sets in the norm and extends uniquely to the completion `H`. Realize `H` as a function space, `H ⊂ L²(dA)`, with weak gradient: the gradient is a closable operator, so a Cauchy sequence that converges to 0 in `L²` has gradients converging to 0. Then the extension is `Q(φ) = ∫(|∇φ|² − 2φ²/R²) dA` with `φ` and `∇φ` the `L²` limits. *Without the continuity estimate, the extension would not be unique. Without closability, `H` would not be a space of functions and "`Q(φ)` for `φ ∈ H`" would be ambiguous.* By R6, this extension is what "second variation" means off the smooth class.

**Established.** The formula, the absence of boundary and vertex terms, and the extension. Checks C1.1–C1.3 are exact symbolic identities.

## 3. T2 — self-adjointness per sector

**Hypotheses.** R1; section 1 (the unfolding); (S1) `e_k(θ) = sin(ν_k(θ + W/R))/√(W/R)` with `ν_k = kπR/(2W)`, `k ≥ 1`, are the orthonormal Dirichlet eigenfunctions on `[−W/R, W/R]` (in `w` they are the same functions rescaled).

**Step 1 (sector reduction).** Let `T = −Δ` on the core `D` (smooth on the unfolded strip, zero on `θ = ±W/R`, zero near both poles). `T` is symmetric:
- the edge terms vanish because both functions are 0 there;
- the seam is interior by R1;
- near the poles everything is 0.

Let `P_k u = e_k ⟨u, e_k⟩_θ`. Then `P_k D = C_c^∞(−π/2, π/2) ⊗ e_k`: the projection of a compactly supported smooth `u` is compactly supported and smooth in `ψ`, and every `f ⊗ e_k` lies in `D`. Also `P_k T = T P_k` on `D`: integrate by parts twice in `θ`, and the boundary terms vanish because `u = e_k = 0` at `θ = ±W/R`. For `u ∈ D`, the partial sums `Σ_{k≤K} P_k u → u` and `T Σ P_k u = Σ P_k T u → T u` in `L²`. So `T̄ = ⊕_k T̄_k`, where `T_k = ℓ_k` on `C_c^∞(−π/2, π/2)`, `ℓ_k f = −(cos ψ f')'/cos ψ + ν_k² f/cos²ψ` (with `R = 1`), and the weight is `cos ψ dψ`. Deficiency spaces add: `ker(T* ∓ i) = ⊕ ker(T_k* ∓ i)`. *Without the commutation, which uses the Dirichlet data on `θ = ±W/R`, the sectors would couple. Without the partial-sum argument, the closure of `T` could be larger than `⊕ T̄_k`.*

**Step 2 (endpoints).** In the unfolded picture `ψ = 0` is an interior regular point, so the only endpoints are `ψ = ±π/2`. Both are the cone point, one from each sheet. *Without the unfolding, `y = 0` and `y = πR` would look like two regular endpoints. The R1 matching joins them, which is the same thing.*

**Step 3 (Frobenius at the poles).** Put `r = π/2 − ψ`. Then `ℓ_k f = λ f` becomes `f'' + (cos r/sin r) f' + (λ − ν²/sin²r) f = 0`, a regular singular point at `r = 0`:
- `r·cos r/sin r → 1` and `r²(λ − ν²/sin²r) → −ν²`, both analytic;
- the indicial roots are `±ν`;
- there is a solution `f₊ = r^ν(1 + O(r²))` and a second one `f₋ ~ r^{−ν}`. If `2ν ∈ ℤ`, `f₋` may carry a `log r · f₊` term, which does not change its size at 0.

With weight `sin r dr ~ r dr`, `∫₀ r^{−2ν} r dr < ∞` iff `ν < 1`, and `f₊` is always square integrable. So by Weyl's alternative (which does not depend on `λ`), the endpoint is **limit-circle iff `ν_k < 1`** and **limit-point iff `ν_k ≥ 1`**. The end `ψ = −π/2` is the same, since the coefficients are even in `ψ`. *If the weight `cos ψ` were dropped (that is, the wrong area element), the threshold would move to `ν < 1/2`. If the indicial roots were not `±ν`, the threshold would move too.*

**Step 4 (the classification).** `ν_k = kπR/(2W)`.
- **`0 < W ≤ πR/2`.** `ν_k ≥ k ≥ 1` for every `k`. Every sector is limit-point at both ends; the borderline `W = πR/2`, `k = 1`, `ν = 1` gives `∫ r^{−1} dr = ∞`, which is still limit-point. Each `T̄_k` has deficiency `(0,0)`. **The closure `T̄` is self-adjoint** (`T` is essentially self-adjoint) and there are no other self-adjoint extensions.
- **`W/R = 17/10`.** `ν₁ = π/3.4 = 0.92400 < 1`, so **sector k = 1 is limit-circle at both ends**. `ν_k ≥ 2π/3.4 = 1.848` for `k ≥ 2`, so those sectors are limit-point. Deficiency indices are `(2,2)` in sector 1 and `(2,2)` in total. **`T̄` is not self-adjoint.** Its self-adjoint extensions are a `U(2)` family of boundary conditions on the coefficients of `f₊`, `f₋` at the two ends. This includes conditions coupling the two sheets, which is natural here because both ends are the same point `p`. In general, sectors with `k < 2W/(πR)` are limit-circle.

**Step 5 (the realization selected by the admissible class, W/R = 17/10, k = 1).** The admissible class `H` is the completion of `D` in the form norm. So the closed form `a(u) = ∫|∇u|²` on `H` is the closure of the form of `T`, and its operator `A` (Kato's representation theorem) is the **Friedrichs extension**. Concretely:
- **`D(A) ⊂ D(T*)`.** By limit-circle theory, near each end `f = α f₊ + β f₋ + (an element of D(T̄))`.
- **`D(T̄) ⊂ H`** (because `T̄ ⊂ A`).
- **`χ f₊ ∈ H`**, where `χ` is a cutoff. Take `η_ε` equal to 0 for `r < ε`, 1 for `r > 2ε`, with slope `1/ε`. The error in the energy is `∫_ε^{2ε} r^{2ν} ε^{−2} r dr + ν²∫_0^{2ε} r^{2ν−1} dr = O(ε^{2ν}) → 0`, which uses `ν > 0`.
- **`χ f₋ ∉ H`**: its sector energy `∫(f₋'² + ν² f₋²/sin²r) sin r dr ~ 2ν² ∫ r^{−2ν−1} dr = ∞`.

Hence `β = 0` at each end separately, and `α` is free at each end. **The admissible class selects the separated condition "no `r^{−ν}` component on either sheet", i.e. `f = O(r^ν)` (in particular `f → 0`) at `p` from both sides.** This is the Friedrichs realization. It does not couple the two sheets. *Without the cutoff estimate, `f₊` might not be admissible and the condition would be Dirichlet-like in a stronger sense. Without the divergence of the `f₋` energy, other extensions in the `U(2)` family would be compatible with the form.*

**Established.** The full classification. All of it is by argument; the computed part is only `ν_k` (in `t3_out.json`, field `limit_circle_sectors`, which is `[1]` for 17/10 and empty otherwise).

## 4. T3 — the spectrum

**Hypotheses.** Sections 1 and 3: the operator `A` is the one defined by the admissible class (the Friedrichs extension; for `W ≤ πR/2` it is the unique self-adjoint extension). Plus (S1).

**Step 1 (Legendre form).** In sector `k`, substitute `x = sin ψ ∈ (−1,1)`. Then `cos ψ dψ = dx` and `(cos ψ f_ψ)_ψ / cos ψ = ((1 − x²) f_x)_x`. So `ℓ_k f = −((1−x²) f')' + ν² f/(1−x²)` on `L²((−1,1), dx)`. *Without the weight `cos ψ`, the measure would not become `dx`.*

**Step 2 (reduction to Gegenbauer).** Write `f = (1−x²)^{ν/2} g`. Direct differentiation gives

`((1−x²)f')' = (1−x²)^{ν/2} [(1−x²)g'' − 2(ν+1)x g' − ν g + ν²x² g/(1−x²)]`.

Hence `ℓ_k f = (1−x²)^{ν/2} [−(1−x²)g'' + 2(ν+1)x g' + ν(ν+1) g]`. So `ℓ_k f = λ f` iff `(1−x²)g'' − (2α+1)x g' + (λ − ν(ν+1)) g = 0` with `α = ν + 1/2`. This is Gegenbauer's equation. It has the polynomial solution `C_n^{(α)}` (degree exactly `n`, since `α > 0`) iff `λ − ν(ν+1) = n(n + 2α) = n(n + 2ν + 1)`. That is:

**`λ_{k,n} = (ν_k + n)(ν_k + n + 1)/R²`, with `ν_k = kπR/(2W)`, `k = 1, 2, …`, `n = 0, 1, 2, …`,**

and the eigenfunction is `f_{k,n} = cos^{ν_k}ψ · C_n^{(ν_k+1/2)}(sin ψ)`. Checks:
- C3.1 is a symbolic check in `ν` for `n = 0..5`.
- C3.2 checks the full separated function `u = f e_k` against `−Δ` in `(ψ,θ)`, symbolically in `ν`, for `n = 0..3`.
- The planted wrong eigenvalue `(ν+n)(ν+n+2)` makes both fire.

In `(y,w)` with `R = 1`:
- `φ = cos^ν y · C_n(sin y) · e_k(w)` for `y < π/2`;
- `φ = (−1)^{n+k} |cos y|^ν C_n(sin y) e_k(w)` for `y > π/2`.

These use `C_n(−x) = (−1)^n C_n(x)` and `e_k(−w) = (−1)^{k+1} e_k(w)`. C3.3 checks the seam matching of the `y`-picture function, values and first three derivatives, for `k = 1..4` and `n = 0..2`.

**Step 3 (these are eigenfunctions of A, not of some other extension).**
- Near `x = ±1`, `1 ∓ x ≈ r²/2`, so `f_{k,n} ~ c r^ν`. This is of `f₊` type and lies in `H` by the cutoff estimate of T2 step 5.
- For `g ∈ C_c^∞`, `a(f_{k,n}, g) = ⟨ℓ_k f_{k,n}, g⟩ = λ ⟨f_{k,n}, g⟩`. This is integration by parts against a compactly supported `g`, so there are no boundary terms.
- Both sides are continuous in the form norm and `C_c^∞` is a form core, so the identity holds for all `g ∈ H`.
- By the representation theorem, `f_{k,n} ∈ D(A)` and `A f_{k,n} = λ_{k,n} f_{k,n}`.

*Without the membership `f_{k,n} ∈ H`, the representation theorem would not apply. For `W/R = 17/10, k = 1`, the solutions with an `r^{−ν}` part are excluded exactly by this step.*

**Step 4 (completeness, hence nothing else in the spectrum).** Suppose `h ∈ L²(dx)` is orthogonal to every `f_{k,n}` (fixed `k`). Put `H = h (1−x²)^{−ν/2}`. Then `H ∈ L²((1−x²)^ν dx)` and `H` is orthogonal there to every `C_n^{(α)}`, so to every polynomial. Polynomials are dense in `L²` of a finite measure on `[−1,1]` (Weierstrass, plus density of continuous functions), so `H = 0` and `h = 0`. So `{f_{k,n}}_n` is an orthogonal eigenbasis of `A_k`. Over all `k` (with `{e_k}` complete) we get an eigenbasis of `A`. Each value `s = ν_k + n` is attained by finitely many `(k,n)` and the eigenvalues tend to `∞`. So **`A` has compact resolvent and its spectrum is exactly `{λ_{k,n}}`**. The multiplicity of `λ` is the number of `(k,n)` with `ν_k + n = s`, where `λ = s(s+1)`. *Without completeness, eigenvalues outside this family could not be ruled out.*

**Lowest eigenvalue.** It is `λ_{1,0} = ν₁(ν₁+1)/R² = π(πR + 2W)/(4W²R)`, with label `(k,n) = (1,0)` and eigenfunction `cos^{ν₁}ψ · sin(ν₁(θ + W/R))`. With `R = 1` this is `π(π + 2W)/(4W²)`.

Labels and closed-form values for the seven widths are in the table in section 5. At `W = π/2`, `ν_k = k` and the spectrum is `l(l+1)` with multiplicity `l`.

**Established.** The closed form. It is exact: a symbolic derivation, checked by symbolic residuals. No float was identified.

## 5. T4 — the finite-element protocol

**Implementation (`t4_fem.py`).**
- P1 elements in `(y,w)`, stiffness `∫(φ_yψ_y + φ_wψ_w/cos²y)|cos y|` and mass `∫φψ|cos y|`.
- `y`-nodes are `π/2 ∓ (π/2)(k/16)²`, `k = 0..16`. There are 8 cells across, and four red refinements follow (R2, R3).
- Dirichlet data on `w = ±W` and on the fiber nodes.
- The seam is imposed by mapping the node `(π, w_j)` to the unknown at `(0, −w_j)` with sign `−1`. The `w`-grid is symmetric, so `−w_j` is a node.
- The six lowest eigenpairs come from `scipy.sparse.linalg.eigsh` with shift-invert at 0, sorted ascending.

On a triangle with an edge on the fiber, `∫1/|cos y|` diverges. There the only free node has zero `w`-gradient (the P1 function vanishes along the fiber edge), so the term is dropped exactly. C4.4 asserts this on every such element at every level.

**Checks, and how each can fail.** Each check was made to fire by a planted defect, run by `run_all.py`; see `run_all.log`, section "planted defects".

| check | what it tests | planted defect that fires it |
| --- | --- | --- |
| C4.1 | total of the unconstrained mass matrix `= ∫|cos y| dy dw = 4W` (rel. 1e-12) | `abscos` (signed `cos`): total 0.000 |
| C4.2 | `y`-stiffness of `φ = y` `= 4W` (P1 reproduces linears) | `ygrad` (drop `1/dy`): 0.137 vs 5.6 |
| C4.3 | `K`, `M` symmetric, `M` SPD | `asym` (one local entry ×1.01) |
| C4.4 | dropped divergent terms meet zero gradients | `fibidx` (check the wrong node) |
| C4.5 | finest-level Rayleigh quotient of the interpolated closed-form bottom eigenfunction within 5% of `ν₁(ν₁+1)` (gross check of seam sign and weights; tunes nothing) | `seam` (periodic instead of anti-periodic): 264.85 vs 2.38 |

**Tables.** The tables are generated by `make_tables.py` from `results.json`. The "extra" columns compare with T3 and are not part of the protocol.

#### W/R = 1/4

| level | N per side | cells across | ndof | λ1 | λ2 | λ3 | λ4 | λ5 | λ6 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 16 | 8 | 217 | 47.1476815893 | 64.3714547952 | 84.6915018190 | 108.6934750033 | 136.4154980252 | 168.6695606584 |
| 1 | 32 | 16 | 945 | 46.1158097873 | 61.3603140582 | 78.8877413710 | 98.8254717310 | 121.2388143889 | 146.2524103195 |
| 2 | 64 | 32 | 3937 | 45.8507607068 | 60.5876239882 | 77.3958886140 | 96.3059030453 | 117.3368920021 | 140.5159219711 |
| 3 | 128 | 64 | 16065 | 45.7839334422 | 60.3929885130 | 77.0199405617 | 95.6722742148 | 116.3549771052 | 139.0745629759 |
| 4 | 256 | 128 | 64897 | 45.7671881775 | 60.3442337480 | 76.9257565912 | 95.5136212459 | 116.1090860824 | 138.7137636496 |

p = 1.996684, extrapolated bottom = 45.7615892819, error estimate = 5.599e-03

#### W/R = 1/2

| level | N per side | cells across | ndof | λ1 | λ2 | λ3 | λ4 | λ5 | λ6 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 16 | 8 | 217 | 13.3126078573 | 22.2194759002 | 33.5685384037 | 47.5942447618 | 49.0137421855 | 64.5646145956 |
| 1 | 32 | 16 | 945 | 13.0868167765 | 21.5264939056 | 32.0818941344 | 44.8064220477 | 46.5844116303 | 59.7650030564 |
| 2 | 64 | 32 | 3937 | 13.0301243623 | 21.3524704626 | 31.7040809952 | 44.0979832459 | 45.9683874194 | 58.5501859128 |
| 3 | 128 | 64 | 16065 | 13.0159304072 | 21.3089083119 | 31.6092235266 | 43.9201125170 | 45.8133785900 | 58.2455624934 |
| 4 | 256 | 128 | 64897 | 13.0123804907 | 21.2980140974 | 31.5854833648 | 43.8755961305 | 45.7745520368 | 58.1693482552 |

p = 1.999420, extrapolated bottom = 13.0111965502, error estimate = 1.184e-03

#### W/R = 1

| level | N per side | cells across | ndof | λ1 | λ2 | λ3 | λ4 | λ5 | λ6 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 16 | 8 | 217 | 4.1201239976 | 9.4810950163 | 13.8719077087 | 17.1000342615 | 23.2813603638 | 27.0897906115 |
| 1 | 32 | 16 | 945 | 4.0586293941 | 9.2550724839 | 13.2272330236 | 16.5131874034 | 21.7956917617 | 25.8624441097 |
| 2 | 64 | 32 | 3937 | 4.0433027231 | 9.1986088226 | 13.0652528523 | 16.3691893425 | 21.4199131927 | 25.5625075009 |
| 3 | 128 | 64 | 16065 | 4.0394735914 | 9.1844946696 | 13.0247140008 | 16.3333259491 | 21.3257768571 | 25.4878394837 |
| 4 | 256 | 128 | 64897 | 4.0385164584 | 9.1809662222 | 13.0145764817 | 16.3243680206 | 21.3022317097 | 25.4691902608 |

p = 2.000226, extrapolated bottom = 4.0381974806, error estimate = 3.190e-04

#### W/R = 7/5

| level | N per side | cells across | ndof | λ1 | λ2 | λ3 | λ4 | λ5 | λ6 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 16 | 8 | 217 | 2.4282858380 | 6.8267585220 | 7.7650797853 | 13.4222387659 | 15.0293106368 | 16.7663490529 |
| 1 | 32 | 16 | 945 | 2.3926822162 | 6.6752602421 | 7.3997509968 | 13.0064837665 | 14.0806251808 | 15.2071295106 |
| 2 | 64 | 32 | 3937 | 2.3838245317 | 6.6374635061 | 7.3094955941 | 12.9032375064 | 13.8456190717 | 14.8231821257 |
| 3 | 128 | 64 | 16065 | 2.3816126002 | 6.6280183693 | 7.2869998074 | 12.8774565389 | 13.7870172098 | 14.7276833601 |
| 4 | 256 | 128 | 64897 | 2.3810597597 | 6.6256572743 | 7.2813800739 | 12.8710128223 | 13.7723762014 | 14.7038400940 |

p = 2.000371, extrapolated bottom = 2.3808755428, error estimate = 1.842e-04

#### W/R = 3/2

| level | N per side | cells across | ndof | λ1 | λ2 | λ3 | λ4 | λ5 | λ6 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 16 | 8 | 217 | 2.1865024844 | 6.4263245361 | 6.9140518166 | 12.8549195651 | 13.8285296756 | 14.8551467853 |
| 1 | 32 | 16 | 945 | 2.1544479192 | 6.2851612727 | 6.5879966135 | 12.4625954396 | 12.9567352005 | 13.4647893694 |
| 2 | 64 | 32 | 3937 | 2.1464748164 | 6.2499481908 | 6.5075889694 | 12.3650833310 | 12.7412629814 | 13.1240244452 |
| 3 | 128 | 64 | 16065 | 2.1444837946 | 6.2411485492 | 6.4875570675 | 12.3407279985 | 12.6875616560 | 13.0393673763 |
| 4 | 256 | 128 | 64897 | 2.1439861476 | 6.2389487041 | 6.4825534386 | 12.3346399375 | 12.6741468272 | 13.0182373457 |

p = 2.000314, extrapolated bottom = 2.1438203134, error estimate = 1.658e-04

#### W/R = pi/2

| level | N per side | cells across | ndof | λ1 | λ2 | λ3 | λ4 | λ5 | λ6 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 16 | 8 | 217 | 2.0398446446 | 6.1798114043 | 6.4016403433 | 12.5035338084 | 13.0966272424 | 13.7071188919 |
| 1 | 32 | 16 | 945 | 2.0099207147 | 6.0448742180 | 6.0992291697 | 12.1253672917 | 12.2713472267 | 12.4190106864 |
| 2 | 64 | 32 | 3937 | 2.0024781130 | 6.0112158016 | 6.0247325760 | 12.0313247125 | 12.0676504257 | 12.1041670094 |
| 3 | 128 | 64 | 16065 | 2.0006194788 | 6.0028041737 | 6.0061784976 | 12.0078313970 | 12.0169008856 | 12.0260043340 |
| 4 | 256 | 128 | 64897 | 2.0001548856 | 6.0007011524 | 6.0015443348 | 12.0019581322 | 12.0042244883 | 12.0064987382 |

p = 2.000203, extrapolated bottom = 2.0000000503, error estimate = 1.548e-04

#### W/R = 17/10

| level | N per side | cells across | ndof | λ1 | λ2 | λ3 | λ4 | λ5 | λ6 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 16 | 8 | 217 | 1.8132823372 | 5.6164590398 | 5.7928633358 | 11.9463898560 | 11.9549594879 | 11.9599043885 |
| 1 | 32 | 16 | 945 | 1.7866129810 | 5.3502722352 | 5.6674720810 | 10.8221677892 | 11.2062932234 | 11.5920255262 |
| 2 | 64 | 32 | 3937 | 1.7799794312 | 5.2848076529 | 5.6361926219 | 10.5469091989 | 11.0206837269 | 11.5033239002 |
| 3 | 128 | 64 | 16065 | 1.7783224013 | 5.2685101338 | 5.6283734585 | 10.4786459308 | 10.9744649295 | 11.4811549184 |
| 4 | 256 | 128 | 64897 | 1.7779080592 | 5.2644400290 | 5.6264178922 | 10.4616153144 | 10.9629217611 | 11.4756106715 |

p = 1.999706, extrapolated bottom = 1.7777699076, error estimate = 1.382e-04

#### Summary

| W/R | p | extrapolated bottom | error estimate | closed form ν1(ν1+1) (T3) | extra: extrap − closed form |
| --- | --- | --- | --- | --- | --- |
| 1/4 | 1.996684 | 45.7615892819 | 5.599e-03 | 45.7616029115 | -1.363e-05 |
| 1/2 | 1.999420 | 13.0111965502 | 1.184e-03 | 13.0111970547 | -5.045e-07 |
| 1 | 2.000226 | 4.0381974806 | 3.190e-04 | 4.0381974271 | +5.351e-08 |
| 7/5 | 2.000371 | 2.3808755428 | 1.842e-04 | 2.3808754887 | +5.410e-08 |
| 3/2 | 2.000314 | 2.1438203134 | 1.658e-04 | 2.1438202624 | +5.098e-08 |
| pi/2 | 2.000203 | 2.0000000503 | 1.548e-04 | 2.0000000000 | +5.030e-08 |
| 17/10 | 1.999706 | 1.7777699076 | 1.382e-04 | 1.7777698463 | +6.131e-08 |

#### Closed-form lowest six (T3), labels (k, n)

| W/R | ν1 | lowest six λ_{k,n} | labels |
| --- | --- | --- | --- |
| 1/4 | 6.283185307180 | 45.7616029115, 60.3279735259, 76.8943441403, 95.4607147546, 116.0270853690, 138.5934559833 | (1,0) (1,1) (1,2) (1,3) (1,4) (1,5) |
| 1/2 | 3.141592653590 | 13.0111970547, 21.2943823619, 31.5775676690, 43.8607529762, 45.7616029115, 58.1439382834 | (1,0) (1,1) (1,2) (1,3) (2,0) (1,4) |
| 1 | 1.570796326795 | 4.0381974271, 9.1797900807, 13.0111970547, 16.3213827342, 21.2943823619, 25.4629753878 | (1,0) (1,1) (2,0) (1,2) (2,1) (1,3) |
| 7/5 | 1.121997376282 | 2.3808754887, 6.6248702412, 7.2795072021, 12.8688649938, 13.7674967072, 14.6958951403 | (1,0) (1,1) (2,0) (1,2) (2,1) (3,0) |
| 3/2 | 1.047197551197 | 2.1438202624, 6.2382153648, 6.4808859473, 12.3326104672, 12.6696761521, 13.0111970547 | (1,0) (1,1) (2,0) (1,2) (2,1) (3,0) |
| pi/2 | 1.000000000000 | 2.0000000000, 6.0000000000, 6.0000000000, 12.0000000000, 12.0000000000, 12.0000000000 | (1,0) (1,1) (2,0) (1,2) (2,1) (3,0) |
| 17/10 | 0.923997839291 | 1.7777698463, 5.2630837066, 5.6257655249, 10.4559415810, 10.9590750638, 11.4737612035 | (1,0) (2,0) (1,1) (3,0) (2,1) (1,2) |

#### T5 table

| W/R | lowest Jacobi eigenvalue ν1(ν1+1) − 2 | sign | from T4 extrap − 2 | negative eigenvalues | nullity |
| --- | --- | --- | --- | --- | --- |
| 1/4 | +43.761602911537 | + | +43.761589282 | 0 | 0 |
| 1/2 | +11.011197054679 | + | +11.011196550 | 0 | 0 |
| 1 | +2.038197427067 | + | +2.038197481 | 0 | 0 |
| 7/5 | +0.380875488666 | + | +0.380875543 | 0 | 0 |
| 3/2 | +0.143820262429 | + | +0.143820313 | 0 | 0 |
| pi/2 | +0.000000000000 | 0 | +0.000000050 | 0 | 1 |
| 17/10 | -0.222230153694 | − | -0.222230092 | 1 | 0 |


**Observations.**
- `p ≈ 2` for every width, including `W/R = 17/10`, where the bottom eigenfunction is only `~ r^{0.924}` at the cone point. The quadratic grading appears sufficient there; I did not test other gradings.
- The extrapolated bottom agrees with the closed form to ≤ 1.4e-5 (at 1/4) and ≤ 6.2e-8 (for `W ≥ 1`). That is far inside the protocol's error estimate. That estimate measures the finest-level error, not the error of the extrapolated value.
- At the finest level, `λ₂…λ₆` lie above the closed-form values by 1e-4 to 1.2e-1 (`finest_minus_exact` in `results.json`). This is consistent in sign with a conforming method.

**Established.** The numbers in the tables, under readings R2–R5. Not settled: whether a different reading of "uniform refinement" (regrading with `N` doubled) would change `p`. I did not run it.

## 6. T5 — sign, index and nullity of the Jacobi operator

**Hypotheses.** T1 (`Q = a − (2/R²)‖·‖²` on `H`); T3 (the spectrum of `A`, compact resolvent).

**Step 1.** The operator of the closed form `Q` on `H` is `J = A − 2/R²`: a bounded perturbation of the form changes the operator by the same constant. Its spectrum is `{λ_{k,n} − 2/R²}`. *If the admissible class selected a different extension, for example a coupled one at 17/10, `J` would be a different operator.*

**Step 2.** The lowest eigenvalue is `(ν₁(ν₁+1) − 2)/R² = (ν₁ − 1)(ν₁ + 2)/R²`. Its sign is `sign(ν₁ − 1) = sign(πR/2 − W)`.

**Step 3.** `s ↦ s(s+1)` is increasing for `s > 0` and equals 2 at `s = 1`. So the negative eigenvalues are the `(k,n)` with `ν_k + n < 1`, i.e. `n = 0` and `k < 2W/(πR)`. The kernel is `ν_k + n = 1`, i.e. `n = 0` and `k = 2W/(πR)` (`n = 1` would need `ν_k = 0`). By min-max, the number of negative eigenvalues equals the Morse index of `Q` on `H`, and the multiplicity of 0 is its nullity. *Without compact resolvent, min-max could count essential spectrum.*

**Result.**
- **`0 < W < πR/2`:** lowest Jacobi eigenvalue positive; 0 negative eigenvalues; nullity 0 (strictly stable).
- **`W = πR/2`:** lowest eigenvalue `0` exactly (`ν₁ = 1` symbolically); 0 negative eigenvalues; **nullity 1**. The Jacobi field is `u = cos ψ cos θ = x₁/R`, the rotation in the `x₁x₄`-plane, which fixes the edge circle `x₁ = 0`.
- **`W/R = 17/10`:** lowest eigenvalue `ν₁(ν₁+1) − 2 = −0.22223015369` (negative); **1 negative eigenvalue** (`2W/π = 1.0823`, so only `k = 1`); nullity 0 (`2W/π` is not an integer).

The values are in the T5 table above. The T4 extrapolations give the same signs; at `π/2` the numerics alone (`+5.0e-8` with error estimate `1.5e-4`) cannot decide the sign, and the zero comes from T3.

**Established.** All of T5 rests on T1 and T3. The 17/10 index assumes the realization selected in T2 (the Friedrichs extension).

## 7. Assessment summary

- **T1.** Established: an exact formula, no boundary or vertex terms, and a unique continuous extension.
- **T2.** Established:
  - for `W ≤ πR/2`, essentially self-adjoint, with every sector limit-point;
  - for `W/R = 17/10`, sector `k = 1` is limit-circle at both sheets, deficiency `(2,2)`, and the Friedrichs (`r^{ν}`, separated) realization is selected.
  - This depends on reading R1 for the operator, but not for the form.
- **T3.** Established: `λ_{k,n} = (ν_k+n)(ν_k+n+1)/R²` with `ν_k = kπR/(2W)`, and lowest `π(πR+2W)/(4W²R)`.
- **T4.** Done per protocol under readings R2–R5. The results agree with T3; that comparison is labelled extra.
- **T5.** Established from T1 and T3:
  - for `W < π/2`: positive, index 0, nullity 0;
  - for `W = π/2`: zero, index 0, nullity 1;
  - for `W = 17/10`: negative, index 1, nullity 0.
- **Disagreements with my own earlier steps.** None found. One process note: the first `run_all.py` invoked `./py` from inside `./py`. The sandbox refused this (exit 65), so that run produced no results. I switched to `sys.executable` and reran; the log in the room is from the rerun. One shell command (a `for` loop not beginning with `./py`) was denied by the harness and was not run; its planted runs were moved into `run_all.py`.

## 8. Consulted-material manifest

**Files read:**
- `BRIEF.md`, `spec_sheet.md`;
- my own outputs `t4_tables.md`, `run_all.log`, `results.json`, `t3_out.json`, `t4_out.json`, read back through scripts.

**Loaded without my asking:**
- a directory listing that showed the name `.room.sb`, which I did not open;
- a sandbox error message naming it;
- the harness-supplied environment context.

Nothing outside the room was read.

**Results used from memory, not derived here:**
- the exponential map of the round sphere;
- the matrix-determinant lemma;
- the general second-variation formula, with `Ric_{S³(R)} = 2/R²` (used only as a cross-check of T1);
- closability of the weak gradient;
- the Weyl limit-point/limit-circle alternative and its independence of `λ`;
- the Frobenius method at a regular singular point;
- von Neumann deficiency theory and the `U(n)` parametrization of extensions;
- the limit-circle decomposition of `D(T*)`;
- Kato's first representation theorem and the Friedrichs extension;
- the Gegenbauer differential equation, and that `C_n^{(α)}` has degree exactly `n` for `α > 0`;
- the symmetry `C_n(−x) = (−1)^n C_n(x)`;
- Weierstrass approximation and density of polynomials in `L²` of a finite measure;
- the min-max principle;
- the Richardson extrapolation formula;
- standard P1 finite-element assembly;
- the library routines `sympy.gegenbauer` and `scipy.sparse.linalg.eigsh`.
