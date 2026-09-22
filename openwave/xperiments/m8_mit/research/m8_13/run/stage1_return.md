# Return, stage 1

Answers to `worklist.md`, derived in this room. Everything called *exact* below is a
symbolic sympy computation over `Q` (or over `Q` adjoined square roots), or a
hand derivation written out here; everything called *numeric* is flagged with its
precision and with what it can and cannot establish.

---

## 0. Item 0

| Quantity | Value | Route | Precision |
| --- | --- | --- | --- |
| `⟨3 3; 3 −3 \| 6 0⟩` | `√231/462 = 1/√924` | sympy `CG(3,3,3,-3,6,0).doit()`, cross-checked against the stretched closed form for all 49 pairs `(m₁,m₂)` | exact (rational multiple of a surd); 50-digit value `0.032897584747988449419196416030388842072689018315373` |
| its square | `1/924` | same | exact |
| `Θ(Θu)` | `u`, for every `u` | symbolic on generic `c_m` | exact |

`Θ² = 1` because `(−1)^{3−m}·(−1)^{3+m} = (−1)^6 = +1` for every integer `m`, and the two
conjugations compose to the identity. Since `3` is an integer spin this is the
`Θ² = (−1)^{2j} = +1` case; I derived it here rather than quoting it.

The stretched closed form I used as the second route is

```text
<j1 m1; j2 m2 | (j1+j2) M>
  = sqrt( (2j1)!(2j2)!(J+M)!(J-M)! / [ (2J)!(j1+m1)!(j1-m1)!(j2+m2)!(j2-m2)! ] )
```

This formula I RECALLED rather than derived; I then verified it against sympy at all
49 pairs `(m₁,m₂)` with `j₁=j₂=3, J=6` (zero mismatches), so nothing rests on the
recollection. With `j₁=j₂=3, m₁=3, m₂=−3, M=0` it gives `sqrt(6!6!/12!) = 1/√C(12,6) = 1/√924`.

Two further exact facts established at the same time, both used later:

| Fact | Status |
| --- | --- |
| `ρ_J(u)_Q = ⟨u\|T^J_Q\|u⟩` for an explicit HS-orthonormal operator set `T^J_Q`, `J=0..6`, `49` operators | orthonormality checked exactly, all `49×49` inner products |
| `Σ_{J=0..6} r̂_J(u) = 1` and `r̂₀(u) = 1/7` **identically**, and `Θ_J ρ_J(u) = (−1)^J ρ_J(u)` for every `J` | checked exactly on symbolic `u` (14 real variables) |

`r̂₀ ≡ 1/7` gives the first, crude bound `r̂₆ ≤ 6/7`. It is **not** attained (see § 5).

---

## 1. The maximum

> **`max r̂₆ = 463/924 = 1/2 + 1/924 = 0.501082251082251082251082251082…`**

### Hypotheses

H1. `V₃ = C⁷` with the worklist's basis, product, `J` action, `Θ` and `ρ₆`; nothing else.
H2. `u ≠ 0`; `r̂₆` is degree-0 homogeneous so the unit sphere is the whole story.
H3. Standard real analysis (compactness), Schur's lemma, the Maclaurin/Newton inequality
`e₂ ≤ e₁²/3` for non-negative reals (proved inline as a sum of squares), Cauchy–Schwarz,
and the classification of closed subgroups of `SO(3)` (recognised, see § 6).

### The argument, in five steps

**Step 1 — `r̂₆` is a Hermitian form on `Sym²(V₃)`, scalar on each isotypic block.**

`r̂₆(u) = Σ_Q |⟨u|T^6_Q|u⟩|²` is of bidegree `(2,2)` in `(u, ū)`, hence
`r̂₆(u) = ⟨u⊗u| N |u⊗u⟩` for a Hermitian `N`, explicitly

```text
|<u|T|u>|^2 = sum_{i,j,k,l} conj(u_i u_k) (T)_{ij} conj((T)_{lk}) (u_j u_l)
N_{(i,k),(j,l)} = sum_Q (T_Q)_{ij} conj((T_Q)_{lk})
```

`r̂₆(D(g)u) = r̂₆(u)` (unitarity of the rank-6 Wigner matrix), and a Hermitian form on
`Sym²` is determined by its values on the Veronese `{u⊗u}` (both sides have `28²` real
dimensions and the map is injective). Hence `P_sym N P_sym` commutes with `D⊗D`.
`Sym²(V₃) = W₆ ⊕ W₄ ⊕ W₂ ⊕ W₀` has dimensions `13+9+5+1 = 28` and is **multiplicity
free** (verified constructively: the 28 CG-built block vectors are exactly orthonormal).
By Schur, `N` is a scalar `n_J` on each block:

> **`r̂₆(u) = n₀ p₀ + n₂ p₂ + n₄ p₄ + n₆ p₆`,  `Σ_J p_J = 1`,  `p_J ≥ 0`**

where `p_J(u)` is the weight of the **holomorphic** square `u⊗u` in block `J` (no `Θ`).

| `J` | `n_J` | `924·n_J` | which is |
| --- | --- | --- | --- |
| 0 | `13/7` | 1716 | `C(13,6)` |
| 2 | `65/84` | 715 | `C(13,4)` |
| 4 | `13/154` | 78 | `C(13,2)` |
| 6 | `1/924` | 1 | `C(13,0)` |

Each `n_J` computed exactly twice, on the highest-weight vector `Q = +J` and on the
`Q = 0` vector of its block; the two agree. Block-scalarity itself is forced by Schur and
was additionally checked numerically over the full `28×28` matrix (max off-block entry
`8.3e−17`). The `C(13, 6−J)` pattern is an observation, not used.

**Step 2 — the maximum is attained on the real form `{Θu = u}`.**

`Θ` is an antiunitary involution, so `V₃ = V_R ⊕ i V_R` with `V_R = {Θu = u}` a real
7-dimensional `SO(3)`-invariant subspace. Write `u = x + iy`, `x,y ∈ V_R`. Then
`uu†` has symmetric part `S = xx^T + yy^T` and antisymmetric part `i(yx^T − xy^T)`.
`K = 6` is **even**, so `T^6_Q` lies in the complexified symmetric part, so

```text
rhat_6(u) = || Pi_6 S ||^2 ,     S >= 0,  tr S = 1,  rank S <= 2,
```

and `S` ranges over exactly `{PSD, trace 1, rank ≤ 2}`. `S ↦ ‖Π₆S‖²` is a **convex**
quadratic form. Its maximum over the convex compact set `{PSD, trace 1}` is attained at an
extreme point, i.e. at a rank-one `S = xx^T`, which already lies in the rank-≤2 set.
Hence the two maxima coincide and the maximum is attained at `Θu = u` (up to phase).
For such `u`, `u⊗Θu = u⊗u`, so `r̂₆ = p₆` there.

**Step 3 — the three-axes formula.**

`Θu ∝ u` means the Majorana constellation of `u` (the 6 roots of its spinor polynomial
`P_u`) is invariant under the antipodal map, i.e. it is 3 axes `±n₁, ±n₂, ±n₃`. Then
`P_u = q₁q₂q₃` with `q_k` the quadratic whose roots are the pair `±n_k`, and for a unit
spinor `z` with Bloch vector `m`, `|q_k(z)|² = (1 − (m·n_k)²)/4`. Using
`‖P‖²_F = n!(n+1)⟨|P|²⟩_{S³}` (derived: checked on every monomial) together with
`p₆ = ‖P²‖²_F / (924 ‖P‖⁴_F)` and the Hopf pushforward of the normalised measure of `S³`
to that of `S²`:

> **`r̂₆ = (13/49) ⟨F²⟩ / ⟨F⟩²`,  `F(m) = Π_k (1 − (m·n_k)²)`, averages over `S²`.**

With `a = n₁·n₂`, `b = n₁·n₃`, `c = n₂·n₃`, `q₁ = a²+b²+c²`, `q₂ = a²b²+b²c²+c²a²`,
`p = abc`, the exact sphere moments give

```text
<F>   = 4 (3 q1 - 2 p + 5) / 105                                       = 4 Den / 105
<F^2> = 64 ( 15 q1^2 + 7 q2 - 20 p q1 + 9 p^2 + 71 q1 - 108 p + 30 )
        * 16 / 45045 ... collected as:   rhat_6 = (20/77) * N / Den^2
  Den = 3 q1 - 2 p + 5
  N   = 15 q1^2 + 7 q2 - 20 p q1 + 9 p^2 + 71 q1 - 108 p + 30
```

The domain is the **elliptope**: `(a,b,c)` is the Gram data of three unit vectors in `R³`
iff `|a|,|b|,|c| ≤ 1` and `det G = 1 + 2p − q₁ ≥ 0`.

**Step 4 — the inequality.** Put `Φ := 463·Den² − 240·N`, so `r̂₆ ≤ 463/924 ⟺ Φ ≥ 0`.
Exactly,

```text
Phi = -1680 q2 + G(q1,p),   G = 567 q1^2 - 308 p^2 - 756 p q1 + 16660 p - 3150 q1 + 4375
```

| Step | Statement | Why |
| --- | --- | --- |
| S1 | `q₂ ≤ q₁²/3`, so `Φ ≥ Ψ(q₁,p) := 7q₁² − 308p² − 756pq₁ + 16660p − 3150q₁ + 4375` | Maclaurin on `x=a², y=b², z=c² ≥ 0`; `e₁²/3 − e₂ = ((x−y)²+(y−z)²+(z−x)²)/6`, verified identically |
| S2 | `∂Ψ/∂q₁ = 14q₁ − 756p − 3150 ≤ −6027/2 < 0` on `q₁∈[0,3], p∈[−1/8,1]`, and `q₁ ≤ 1+2p`, so `Ψ ≥ Ψ(1+2p, p) =: ψ(p)` | Gram PSD; `q₁ ≤ 3` from `\|a\|,\|b\|,\|c\| ≤ 1` |
| S3 | `ψ(p) = −1792p² + 9632p + 1232 = 1792 (p + 1/8)(11/2 − p)` | exact factorisation; roots `−1/8` and `11/2` |
| S4 | `p = abc ≥ −1/8` for any three unit vectors | AM-GM `q₁ ≥ 3\|p\|^{2/3}` with `q₁ ≤ 1+2p`; for `p = −σ³` this is `1 − 2σ³ − 3σ² ≥ 0`, i.e. `−(σ+1)²(2σ−1) ≥ 0`, i.e. `σ ≤ 1/2` |

`p ≤ 1 < 11/2`, so `ψ(p) ≥ 0`, so `Φ ≥ 0`. **`r̂₆ ≤ 463/924` for every unit `u ∈ C⁷`.**

**Step 5 — attained.** `u = (v₃ + v₋₃)/√2` satisfies `Θu = u` exactly and gives, exactly,

```text
rhat_K = ( 1/7 , 0 , 25/84 , 0 , 9/154 , 0 , 463/924 )   K = 0..6,  sum = 1
rho_6  = ( 1/2 , 0,0,0,0,0 , sqrt(231)/462 , 0,0,0,0,0 , 1/2 ),   ||rho_6||^2 = 1/2 + 1/924
```

Its constellation is the regular hexagon on a great circle = 3 coplanar axes at 60°,
Gram `(a,b,c) = (−1/2,−1/2,−1/2)`, `det G = 0`. Note `1/2 + 1/924 = 463/924`, and the
`Q = 0` entry is precisely item 0's coefficient.

---

## 2. The complete set of maximisers

> **The maximum is attained exactly on**
> **`M = { e^{iφ} · D³(g) · (v₃ + v₋₃)/√2  :  φ ∈ [0,2π), g ∈ SO(3) }`**,
> equivalently: the unit vectors whose Majorana constellation is a **regular hexagon
> inscribed in a great circle**. This is a single `SO(3)×U(1)` orbit, a compact
> 4-real-dimensional manifold `(SO(3)/D₆) × U(1)`; the projective stabiliser is the
> dihedral group `D₆` of order 12.

### Shape of the completeness argument

It splits into three branches, which together cover every unit `u ∈ C⁷`.

```text
  every unit u
      |
  (B1) is rank(S) = 1, i.e. Theta u = (phase) u ?
      |-- no  -> (B2) rank 2: excluded, no maximiser there
      |-- yes -> real form, 3 axes, point of the elliptope
                    |
                (B3) equality analysis of Steps S1-S4 pins the point uniquely
```

**B1 / B2 — only rank-one `S` can be a maximiser.** With `S = Σ_i w_i x_i x_i^T` its
spectral decomposition (`x_i` real orthonormal, `w_i > 0`, `Σ w_i = 1`, at most two terms),

```text
f(S) = || sum_i w_i Pi_6(x_i x_i^T) ||^2
     <= ( sum_i w_i ||Pi_6(x_i x_i^T)|| )^2  <= max_i f(x_i x_i^T) <= 463/924
```

Equality forces (i) every `x_i` to be a real maximiser and (ii) all the vectors
`Π₆(x_i x_i^T) = ρ₆(x_i)` to be **equal** (equal norms plus equality in the triangle
inequality). By B3 every real maximiser is `x = ±D(g)h` with `h` the hexagon state, and
`ρ₆(D(g)h) = D⁶(g) ρ₆(h)`, so `ρ₆(x₁) = ρ₆(x₂)` forces `g₂^{-1}g₁ ∈ Stab(ρ₆(h))`.
Computed exactly: `ρ₆(h)` has non-zero entries only at `Q = ±6` (value `1/2`) and `Q = 0`
(value `1/√924`). Its stabiliser is `D₆`:

- rotation by `θ` about `z` multiplies the `Q` entry by `e^{−iQθ}`; invariance needs
  `e^{−6iθ} = 1`, giving `C₆`; the `π` rotation about `x` sends `Q → −Q` and the vector is
  `Q → −Q` symmetric, giving `D₆`;
- a closed subgroup `H ⊇ D₆` contains an element of order 6, so `H ∈ {C_{6k}, D_{6k},
  SO(2)_z, O(2)_z, SO(3)}` (`T`, `O`, `I` have no order-6 element). `C₁₂` fails (rotation by
  `π/6` negates the `Q = ±6` entries; checked numerically, deviation `1.414`),
  `SO(2)_z/O(2)_z` would force all `Q ≠ 0` entries to vanish, `SO(3)` would force `ρ₆ = 0`.
  The 12 elements of `D₆` were checked to fix `ρ₆(h)` to `2.1e−15`.

`Stab(ρ₆(h)) = D₆` is also the projective stabiliser of `h`, so `x₁` and `x₂` span the same
ray; both being real, `x₁ = ±x₂` and `x₁x₁^T = x₂x₂^T`, contradicting orthogonality in a
rank-2 decomposition. Hence rank `S = 1`: every maximiser is `Θ`-invariant up to phase.

**B3 — equality in the chain S1–S4 pins the configuration.** Equality in Step 4 requires
all three of:

| Condition | Consequence |
| --- | --- |
| `q₂ = q₁²/3` | `a² = b² = c²` (Maclaurin equality) |
| `q₁ = 1 + 2p` | `det G = 0`: the three axes are **coplanar** (`∂Ψ/∂q₁ < 0` strictly, so the boundary is forced) |
| `ψ(p) = 0` | `p = −1/8` (the other root `11/2` is impossible since `\|abc\| ≤ 1`) |

Together: `3t² = 1 + 2(−1/8) = 3/4`, so `|a| = |b| = |c| = 1/2` and `abc = −1/8`. The four
sign patterns `(1/2,1/2,−1/2)`, `(1/2,−1/2,1/2)`, `(−1/2,1/2,1/2)`, `(−1/2,−1/2,−1/2)` all
give `R = 463/924` and `det G = 0`, and all describe the **same** axis configuration:
flipping the sign of one `n_k` flips exactly two of `(a,b,c)`. The Gram matrix
`[[1,−½,−½],[−½,1,−½],[−½,−½,1]]` has eigenvalues `3/2, 3/2, 0`, so rank 2: three coplanar
unit vectors at `120°`, i.e. three coplanar **axes at 60°**, unique up to `SO(3)` (planar, so
the mirror image is a rotation). The constellation determines the state up to phase.

### What would fail if each case were omitted

| Case | What is specifically lost if omitted |
| --- | --- |
| B1 (convex reduction to rank one) | The three-axes formula holds **only** for `Θ`-invariant states. Without B1 there is no reduction at all: the problem stays a quartic on a 13-real-dimensional sphere with no exact handle, and no bound better than the crude `6/7` survives. |
| B2 (rank-two exclusion) | The **value** `463/924` would still be proved, but the reported maximiser set would be wrong: one could not exclude a continuum `S = w·xx^T + (1−w)·yy^T`, `0 < w < 1`, of genuinely complex maximisers. What is lost is exactly the assertion "every maximiser satisfies `Θu = e^{iφ}u`". |
| B3-S1 (Maclaurin) | `Φ` contains `−1680 q₂` and `q₂` is not bounded above by `q₁` and `p` alone without it; **no** inequality follows, in either direction. |
| B3-S2 (`q₁ ≤ 1+2p`, Gram PSD) | `Ψ` is not sign-definite on the box: e.g. `Ψ(q₁=3, p=−1/8) = −6815.8 < 0`. Omitting the PSD constraint makes the final step false, not merely weaker. |
| B3-S3/S4 (`p ≥ −1/8`) | `ψ(p) = 1792(p+1/8)(11/2−p)` is **negative** for `p < −1/8`. Without S4 the claimed bound is unproved and the last step of the chain reverses sign. |
| Degenerate constellations (repeated axes, `\|a\|=1`) | These are interior to nothing but are points of the closed elliptope and are covered by the same chain. Omitting them would leave unchecked e.g. `(1,1,1)` (state `v₀`, `r̂₆ = 100/231`) and `(1,0,0)` (state `(v₁−v₋₁)/√2`, `r̂₆ = 145/308`); neither is a maximiser, but the bound must hold there and does. |
| Merging the four sign patterns | One would report four maximiser families where there is one orbit; the structure of the answer (a single `SO(3)×U(1)` orbit with stabiliser `D₆`) would be misstated. |

### Does anything rest on a solver?

The Gröbner-basis route to the interior critical points (`s09_maximise.py`) produced **no
output within the session budget and was killed**; nothing in the argument above uses it.
sympy is used for polynomial expansion, factorisation and CG coefficients only; every
factorisation it produced (`ψ(p) = −112(2p−11)(8p+1)`, `1−2σ³−3σ² = −(σ+1)²(2σ−1)`,
`e₁²/3 − e₂` as a sum of squares) is verifiable by hand in one line, and I re-checked each
by substitution. **If I relied on a solver alone** the exposure would be: a wrong CG
convention or a wrong selection rule silently producing a different quartic — which is
exactly the failure that did occur once (§ 7), and it was caught only by cross-checking
against an independently computed value, not by the solver.

---

## 3. The minimum and its minimisers

> **`min r̂₆ = 1/924 = 0.00108225108225108225108225108225…`**
> **attained exactly on `C = { e^{iφ} · D³(g) · v₃ : φ ∈ [0,2π), g ∈ SO(3) }`,**
> the spin-3 **coherent states**: Majorana constellation = 6 coincident points.
> A single `SO(3)×U(1)` orbit, `S² × U(1)`, 3 real dimensions; projective stabiliser
> `SO(2)` (the rotations about the constellation's own axis).

### Argument

From Step 1, `r̂₆ = Σ_J n_J p_J` with `p_J ≥ 0`, `Σ p_J = 1`, and

```text
n_6 = 1/924  <  n_4 = 13/154  <  n_2 = 65/84  <  n_0 = 13/7      (strict)
```

so `r̂₆ ≥ n₆ = 1/924`, with equality **iff** `p₆ = 1` (strictness of the minimum is what
makes the equality set a single face). So the minimisers are exactly the `u` with
`u⊗u ∈ W₆`.

**`p₆ = 1 ⟺ u coherent.`** With `b_k = c_{k−3}` and `P_u = Σ_k b_k z₁^k z₂^{6−k}/√(k!(6−k)!)`,

```text
p_6 = sum_n |r_n|^2 n!(12-n)! / [924 (sum_k |b_k|^2)^2],
      r_n = sum_k b_k b_{n-k} / sqrt(k!(6-k)!(n-k)!(6-n+k)!)
```

Cauchy–Schwarz on each `n`, with `sum_k 1/(k!(6−k)!(n−k)!(6−n+k)!) = C(12,n)/(6!6!)`
by **Vandermonde** (verified for `n = 0..12`), gives
`|r_n|² n!(12−n)! ≤ 924 Σ_k |b_k|²|b_{n−k}|²`; summing over `n` gives `p₆ ≤ 1`.

Equality for every `n` means, with `β_k := b_k/√C(6,k)`, that `β_k β_{n−k}` is independent
of `k` for each `n`. Sub-cases:

| Sub-case | Conclusion |
| --- | --- |
| exactly one `β_k ≠ 0`, `1 ≤ k ≤ 5` | contradiction: for `n = 2k` the index `j = k−1` is in range and gives `β_{k−1}β_{k+1} = 0 ≠ β_k²` |
| exactly one `β_k ≠ 0`, `k ∈ {0,6}` | allowed: `P ∝ z₂⁶` or `z₁⁶`, which **are** coherent (`v₋₃`, `v₃`) |
| two or more non-zero | support is a contiguous block `[i,j]`; `n = 2i` forces `i = 0` and `n = 2j` forces `j = 6`; then `β_kβ_{k+2} = β_{k+1}²` makes `β` geometric, `β_k = βq^k`, so `a_k ∝ q^k C(6,k)` and `P ∝ (qz₁ + z₂)⁶` — a perfect sixth power, i.e. coherent |

Checked exactly: `p₆ = 1` for `b_k = q^k√C(6,k)` (all `q > 0`), for `b = e₆` and for `b = e₀`;
and `p₆ = 100/231` for `b = e₃`, `6/11` for `b = e₅`, `463/924` for the hexagon — so the
excluded monomials really are excluded, not merely unproved.

Attainment: `r̂₆(v₃) = ⟨3 3; 3 −3 | 6 0⟩² = 1/924`, computed exactly. Its full spectrum is
`r̂_K(v₃) = ⟨3 3;3 −3|K 0⟩²`, exactly `(1/7, 9/28, 25/84, 1/6, 9/154, 1/84, 1/924)`, summing
to 1. The minimiser is **not** `Θ`-invariant: `Θv₃ = v₋₃` and `⟨v₃, Θv₃⟩ = 0` exactly.

### What would fail if each case were omitted

| Case | What is specifically lost |
| --- | --- |
| Step-1 reduction (Schur, multiplicity-freeness) | There is no lower bound at all; `r̂₆ ≥ 0` is all one has, and `0` is not attained. |
| Strictness `n₆ < n₄` | If the smallest `n_J` were tied, the equality set would be a larger face (`p` supported on the tied blocks) and the minimiser set would be strictly bigger than the coherent orbit. Strictness is the whole content of "the minimisers are exactly `p₆ = 1`". |
| `p₆ = 1 ⟹ coherent` | The answer would read "the `u` with `u⊗u ∈ W₆`" — a correct but unresolved description, not an explicit set. |
| The sub-case "one non-zero `β_k`, `1 ≤ k ≤ 5`" | These are the states `v₂, v₁, v₀, v₋₁, v₋₂` and their rotations, which are **not** coherent. Omitting the sub-case leaves a gap through which non-coherent states could enter the minimiser set. (It turns out empty, but that is a result, not an assumption: `p₆(v₀) = 100/231`, `p₆(v₂) = 6/11`.) |
| The sub-case `k ∈ {0,6}` | These *are* coherent; omitting them would wrongly exclude `v₃` and `v₋₃`, i.e. the very states that attain the minimum. |
| Attainment | Without exhibiting `v₃`, `1/924` is a lower bound and the "minimum" is unproved to exist. (Compactness gives existence, not the value.) |

---

## 4. Item 4 — what was computed, what equals a bound, what is missing

### Computed and found to vanish (each with its exact reason)

| Quantity | Reason |
| --- | --- |
| `Θ(Θu) − u = 0` | `(−1)^{3−m}(−1)^{3+m} = 1` |
| `r̂₀(u) − 1/7 = 0` identically | the `J=0` multipole of `u⊗Θu` is `‖u‖²/√7` because `Θ` supplies exactly the invariant pairing; checked symbolically |
| `Σ_J r̂_J − 1 = 0` identically | the CG matrix is a unitary change of basis of `V₃⊗V₃`; the 49 `T^J_Q` are HS-orthonormal (all `49×49` checked) |
| `Θ_J ρ_J − (−1)^J ρ_J = 0`, all `J` | `ρ_J(Θu,Θw) = Θ_J ρ_J(u,w)` plus `ρ_J(Θu,u) = (−1)^J ρ_J(u,Θu)`; checked symbolically for `J = 0..6` |
| `r̂₁ = r̂₃ = r̂₅ = 0` at every maximiser | maximisers are `Θ`-invariant, so `uu†` is a real symmetric matrix in the real basis, while the odd-`K` multipole operators span the **antisymmetric** part; exact, and confirmed exactly at the hexagon |
| `p₀ = p₂ = p₄ = 0` at every minimiser | `p₆ = 1` exactly and the `p_J` are non-negative and sum to 1 |
| `⟨v₃, Θ v₃⟩ = 0` | `Θv₃ = v₋₃`, orthogonal |
| `det Gram = 0` at the maximising configuration | eigenvalues `3/2, 3/2, 0`: the three axes are coplanar |
| `Φ = 463·Den² − 240·N = 0` at the maximising Gram | exact substitution |
| `49p₂ − 77p₆ + 24‖u‖⁴ = 0` and `49p₄ + 126p₆ − 66‖u‖⁴ = 0` on `V_R` | verified **exactly** on the 7-parameter real form (these were first read off numerically in `s06`; see § 7) |
| off-block entries of `N` on `Sym²` | forced to vanish by Schur given multiplicity-freeness; **numerically** `8.3e−17`, not computed exactly entry by entry — see "missing" below |

### Found equal to a bound

| Equality | Status |
| --- | --- |
| `r̂₆(coherent) = n₆ = 1/924` | computed both ways (directly as `⟨33;3−3\|60⟩²` and as the block constant); they agree exactly |
| `r̂₆(hexagon) = 463/924` | computed three ways: CG route on the state; Fock-norm route `‖P²‖²_F/(924‖P‖⁴)`; three-axes formula. All exactly equal |
| `p₀ = 1/7` on the real form | equals its own maximum: `p₀ = \|u^T G u\|²/7 ≤ ‖u‖⁴/7` with equality iff `Θu = e^{iφ}u` |

### Bounds derived here that are **not** attained (reported as bounds, never as values)

| Bound | Where it comes from | Gap to the truth |
| --- | --- | --- |
| `r̂₆ ≤ 1` | `Σ_J r̂_J = 1` | far |
| `r̂₆ ≤ 6/7 = 792/924` | `r̂₀ ≡ 1/7` | not attained: it would require the little group of `u` to kill `J = 1..5`, which among closed subgroups of `SO(3)` only the icosahedral group does, and spin 3 contains **no** icosahedral invariant (`I` is perfect, so a projective invariant is a genuine one; `I`-invariants occur at `j = 0,6,10,…`) |
| `r̂₆ ≤ 13/14` | LP with only `p₀ ≤ 1/7` | not attained |
| `r̂₆ ≤ 11/21 = 484/924` | LP over `{p ≥ 0, Σp = 1, r̂_K(p) ≥ 0 ∀K}`; on the real form the same number is exactly the `p₄ ≥ 0` bound | not attained; the true max is `463/924`, so the LP relaxation has a genuine gap of `21/924` |
| `r̂₆ ≤ max_J n_J = 13/7` | the convex-combination bound | vacuous (`> 1`): `p₀ = 1` is impossible for a rank-one tensor |
| `r̂₆ ≥ 24/77` on the real form | `p₂ ≥ 0` | attained there (orthogonal axes, the octahedron state), but it is **not** the global minimum — the global minimiser is not `Θ`-invariant |

### Missing — computed nowhere, reported as missing rather than as zero

| Item | Status |
| --- | --- |
| The full `28×28` matrix `N` in exact arithmetic | **not computed**. Only (i) the four diagonal block constants, each twice, exactly, and (ii) the whole matrix numerically at `1e−16`. Block-scalarity is asserted from Schur plus multiplicity-freeness, not from an exact entry-by-entry computation. |
| Interior critical points of `R` on the elliptope, by Gröbner basis | **not computed**: `s09_maximise.py` ran past the session budget with no output and was killed. The exact chain in `s12` does not need it, but I did not independently enumerate the critical points. |
| A symbolic proof of `\|q_n(z)\|² = (1−(m·n)²)/4` and of the Hopf pushforward | **not symbolically proved here**. Verified numerically to `1.4e−16` on random inputs, and the end-to-end three-axes formula was verified **exactly** at seven configurations (§ 5). |
| An SOS certificate for the quartic `(463/924)‖u‖⁴ − r̂₆` | **not found**, and shown not to exist at the natural degree: the invariant degree-4 SOS cone is `{Σ_J c_J q_J : c_J ≥ 0}` and the needed form has a negative `q₂` coefficient. A higher-degree Positivstellensatz certificate was not attempted. |
| The value of `r̂₆` on the second-largest critical plateau | **not computed**. All 200 Powell restarts and all 400 complex restarts landed on the single plateau `463/924`; no second plateau was resolved, so I report none rather than guessing. |

---

## 5. Numerical routes: precision, and what each can and cannot establish

| Script | Route | Precision | Establishes | Does **not** establish |
| --- | --- | --- | --- | --- |
| `s03_numeric.py` | self-consistent eigen-iteration `u ← normalise((±H(u)+sI)u)` on the complex unit sphere, 400 random starts | IEEE double, `~1e−15` | existence of states with `r̂₆ = 0.501082251082252`; rotation and phase invariance to `1.1e−15`; `D³(n,2π) = I` to `2.0e−15`; a single ascent plateau over 400 starts | that `463/924` is the maximum; nothing about completeness. **Descent did not reach the minimum**: best `0.001082296499157` vs `1/924 = 0.001082251082251`, an excess of `4.5e−8` (slow convergence at the degenerate coherent critical point). Reported as found, not rounded to `1/924`. |
| `s04b` block check | full `28×28` `N` numerically | `~1e−16` | block-scalarity consistent with Schur; `Σ n_J p_J = r̂₆` on random `u` to `1e−15` | exactness of the off-block zeros |
| `s06_realform.py` | rank test of `{1, p₂, p₄, p₆}` on `V_R`, 60 samples | double; singular values `8.75, 1.10, 1.6e−15, 4.9e−16` | that two affine relations exist (rank 2) | the relations themselves — those were then **verified exactly** in `s15` |
| `s07` CHECK 5 | Monte Carlo for `⟨F⟩, ⟨F²⟩`, `4e5` samples | `~1e−3` | that the exact sphere moments are right to `1e−3` | anything sharper |
| `s08` | formula vs direct state evaluation, random axes | `~1e−16` | agreement of the three-axes formula with the CG definition | a proof; superseded by the seven exact checks |
| `s10_rank2_case.py` | 60 stored near-maximisers; 300-start search for orthogonal real maximisers with equal `Π₆` image | double | all 60 satisfy `\|⟨u,Θu⟩\| = 1.000000000000000`; no rank-two maximiser pair found (best residual `4.3e−1`, i.e. the search never got near one) | non-existence — that is settled by the stabiliser argument in § 2 |
| `s12` | `min abc` over `6e5` random unit triples | double | `−0.1248821870` vs the exact `−1/8` | the bound itself, which is proved in S4 |
| `s14` | 200 Powell restarts over axis triples | double | single plateau `0.501082251082251`; `\|Gram\| = (1/2,1/2,1/2)`; `det G = −3.3e−16` — exactly the configuration the equality analysis predicts | the maximum, or that there is no other plateau |

**Second-precision repeat.** Every reported number was re-evaluated at 40 decimal digits
with mpmath and agrees with its exact rational value, e.g.
`463/924 = 0.501082251082251082251082251082`,
`1/924 = 0.00108225108225108225108225108225`,
`1/√924 = 0.0328975847479884494191964160304`.
The identifications are not float-to-rational guesses: each was obtained as an exact
symbolic value first, and the float is the check, not the source.

**Exact checks of the three-axes formula** (seven configurations, all exact):

| Configuration | Gram `(a,b,c)` | `r̂₆` direct (CG on the state) | formula | equal |
| --- | --- | --- | --- | --- |
| hexagon (3 coplanar, 60°) | `(1/2,−1/2,1/2)` | `463/924` | `463/924` | yes |
| three coincident axes (state `v₀`) | `(1,1,1)` | `100/231` | `100/231` | yes |
| orthogonal (octahedron) | `(0,0,0)` | `24/77` | `24/77` | yes |
| two coincident + one orthogonal | `(1,0,0)` | `145/308` | `145/308` | yes |
| asymmetric 1 | `(0,0,3/5)` | `44925/111188` | `44925/111188` | yes |
| asymmetric 2 | `(4/5,0,9/25)` | `45950425/100420628` | `45950425/100420628` | yes |
| coherent (not `Θ`-invariant, outside the formula's scope) | — | `1/924` | n/a | — |

---

## 6. Tools, and what was looked up rather than derived

**Tools.** `python3` (the room's `python3`; there is no `./py` in this directory and the
brief names `python3`, whose `sympy 1.14.0 / numpy 2.5.3 / mpmath 1.3.0` match the brief).
`sympy` for exact CG coefficients, symbolic algebra, factorisation and the one attempted
Gröbner basis; `numpy` and `scipy.optimize` (Nelder–Mead, Powell) for the numerical
searches; `mpmath` for the 40-digit repeat.

**Recognised from prior knowledge rather than derived here**, each flagged with how it was
handled:

| Fact | How used | Verified here? |
| --- | --- | --- |
| Closed-form for the **stretched** CG coefficient | as a second route for item 0 | yes, against sympy at all 49 pairs |
| **Majorana** representation: a spin-`j` state ↔ `2j` points on the sphere via the roots of its spinor polynomial; `Θ` sends the constellation to its antipode | the structural basis of Step 3 | consequences verified exactly at seven configurations and numerically to `1e−16` |
| **Wick/Isserlis pairing rule for sphere moments**, `⟨(m·v₁)…(m·v_{2k})⟩ = Σ_matchings ∏v_i·v_j / (2k+1)!!` | to compute `⟨F⟩, ⟨F²⟩` | yes, against Monte Carlo (`1e−3`) and against the hand-computed `⟨F⟩ = (20+12S−8D)/105` |
| **Hopf pushforward**: the normalised measure of `S³` maps to that of `S²` | Step 3 | indirectly, via the seven exact checks |
| **Schur's lemma** and the decomposition `Sym²(spin 3) = 6⊕4⊕2⊕0` | Step 1 | multiplicity-freeness verified constructively (28 orthonormal CG block vectors, dims `1+5+9+13`) |
| **Classification of closed subgroups of `SO(3)`** (`C_n, D_n, T, O, I, SO(2), O(2), SO(3)`) and their invariant degrees | the `D₆` stabiliser argument and the "`6/7` is unattainable" remark | not re-derived; the specific facts used (`T,O,I` have no order-6 element; `I`-invariants start at `j=6`) were used only in arguments whose conclusions were independently confirmed numerically |
| Maclaurin `e₂ ≤ e₁²/3`, Vandermonde `Σ_k C(6,k)C(6,n−k) = C(12,n)`, Cauchy–Schwarz | the two inequality chains | both identities verified symbolically here |

**Do I recall a published treatment of this exact maximisation?** I recognise the general
setting — `r̂_K` are the **state multipoles / statistical tensors** of the pure spin-3
density matrix, and states with vanishing low multipoles are called **anticoherent**;
there is a literature on extremising multipole content (the "Kings of Quantumness" line of
work). I do **not** recall a specific published statement of the value `463/924`, of the
hexagon as the maximiser, or of the four constants `n_J = C(13,6−J)/924`. Nothing above
was taken from that recollection; I mention it because the brief asks, and because a reader
should discount any agreement with that literature accordingly. I did **not** look anything
up during this run — the room has no network and I read no file outside it.

**One recollection that I tried and discarded.** I attempted to identify the recoupling
matrix `A[K,J]` with `(2J+1)(2K+1){3 3 K; 3 3 J}²` (a Racah crossing matrix). The ratios
came out as `7, −7, 7, −7, …` in the `J=0` column and as `28/15, 84/19, 42/5, …` elsewhere,
i.e. **not** constant, so the identification is wrong as written and I dropped it. It plays
no role. I report it because it was a real step that failed.

---

## 7. Disagreement with my own earlier steps

Reported as they happened, not smoothed.

| # | What I first had | Why it changed |
| --- | --- | --- |
| 1 | **Index-pairing bug in `N`.** In `s04_reduction.py` I paired the bra/ket indices as `N_{(i,l),(j,k)} = Σ_Q (T_Q)_{ij} conj((T_Q)_{lk})` and evaluated it as `A† · (T B T†)`. | The block constants came out **inconsistent**: `n₂` took two values, `{125/924, −25/462}`, on the same block, which is impossible if Schur applies. Re-expanding `\|⟨u\|T\|u⟩\|²` by hand showed the bra pair is `(i,k)` and the ket pair `(j,l)`, and the contraction is `A† · (T B T̄)` with `T̄` the **entrywise** conjugate (one transpose too many in the original). Corrected in `s04b_reduction_fix.py`; the constants then came out single-valued and the identity `r̂₆ = Σ n_J p_J` verified on random `u`. |
| 2 | **`s07` CHECK 4 appeared to refute the three-axes formula** (discrepancies of `0.1`–`0.5`). | The failure was in the *checker*, not the formula. `s07` defined a local CG helper whose selection guard was `m₁+m₂ = Q` while it called `CG(3,m₁,3,−m₂,…)`, which needs `m₁−m₂ = Q`; the resulting `T⁶` matrix was wrong. `s08_debug_convention.py` re-tested the formula without any axes at all (Monte Carlo over `S³`) and then with a corrected `T⁶`, getting agreement to `1e−16`. The other checks in `s07` used the exact symbolic route and were unaffected. This is the one place where a bug nearly reversed a conclusion. |
| 3 | I first conjectured `max r̂₆ = 6/7`, reasoning that a state with icosahedral little group would kill `J = 1..5`. | Spin 3 contains **no** icosahedral invariant (and `I = A₅` is perfect, so projective invariance collapses to genuine invariance). The bound `6/7` is therefore not attained. Kept in § 4 as a derived-but-unattained bound. |
| 4 | I then tried an LP over `{p ≥ 0, Σp = 1, r̂_K ≥ 0}` and got `11/21 = 484/924`. | Also not tight (`463/924` is the truth). The LP relaxation has a real gap; recorded rather than discarded. |
| 5 | I tried to prove `μ_max(M) ≤ ½‖T‖²` for the contraction `M_{ab} = T_{acd}T_{bcd}` of a harmonic 3-tensor, as a route to the bound. | **False.** The axially symmetric state (`T₁₁₃ = T₂₂₃ = 1, T₃₃₃ = −2`, i.e. `v₀`) has `μ_max = 3/5 > 1/2`. Abandoned; the three-axes route replaced it. |
| 6 | The affine relations `p₂ = (11/7)p₆ − 24/49`, `p₄ = 66/49 − (18/7)p₆` on `V_R` were first read off **numerically** from a rank test plus two sample points. | I flagged them as numerically identified and kept them out of the proof chain; then verified both **exactly** on the 7-parameter real form in `s15`. They now stand as exact, and as an independent cross-check (on `V_R`, `r̂₆ = (7/11)p₂ + 24/77`, so the maximum of `r̂₆` is the maximum of the quadrupole weight `p₂`). |
| 7 | `s09_maximise.py` (Gröbner + heavy Nelder–Mead) was my first plan for completeness. | It exceeded the session budget with no output and was killed. Replaced by the closed-form chain in `s11`/`s12` and the fast search in `s14`. The script is kept, unedited, as a record of the route that did not finish. |

---

## 8. Summary table

| Item | Answer | Attained on |
| --- | --- | --- |
| `⟨3 3;3 −3\|6 0⟩` | `√231/462 = 1/√924` | — |
| `Θ(Θu)` | `u` | — |
| **max `r̂₆`** | **`463/924`** | `{e^{iφ}D³(g)(v₃+v₋₃)/√2}`: regular-hexagon constellation on a great circle; one `SO(3)×U(1)` orbit, 4-dimensional, stabiliser `D₆` |
| **min `r̂₆`** | **`1/924`** | `{e^{iφ}D³(g)v₃}`: the spin-3 coherent states, 6 coincident Majorana points; one `SO(3)×U(1)` orbit, 3-dimensional, stabiliser `SO(2)` |

Both extremal values are `1/924` apart from a round number: `min = 1/924`,
`max = 1/2 + 1/924`. The `1/924 = ⟨3 3;3 −3\|6 0⟩²` of item 0 appears in both, as the
whole of the minimum and as the `Q = 0` component `1/√924` of `ρ₆` at the maximum.

### Scripts written (all kept, in this directory)

| Script | Role |
| --- | --- |
| `s01_setup.py` | exact CG, item 0, `Θ²`, operator orthonormality |
| `s02_structure.py` | `Σ_J r̂_J = 1`, `r̂₀ ≡ 1/7`, `Θ_J ρ_J = (−1)^J ρ_J`, all exact |
| `s03_numeric.py` | complex-sphere ascent/descent, 400 starts (corroboration) |
| `s04_reduction.py` | first, **buggy** build of `N` — kept as the record of § 7 #1 |
| `s04b_reduction_fix.py` | corrected `N`; the four constants `n_J` |
| `s05_recoupling.py` | full `p → r̂` matrix; the LP relaxation (bound `11/21`) |
| `s06_realform.py` | real-form reduction checks; numerical rank test |
| `s07_axes_formula.py` | exact `⟨F⟩`, `⟨F²⟩`; the three-axes formula (its CHECK 4 is the bug of § 7 #2) |
| `s08_debug_convention.py` | isolation and repair of that bug |
| `s09_maximise.py` | Gröbner route; **did not finish**, killed, kept unedited |
| `s10_rank2_case.py` | the rank-two case |
| `s11_exact_bound.py` | `Φ = 463Den² − 240N`; discovery of the chain |
| `s12_chain_verify.py` | exact verification of every step S1–S4 and the equality analysis |
| `s13_extremal_states.py` | exact multipole spectra; `p₆ = 1 ⟺ coherent`; the `D₆` stabiliser |
| `s14_numeric_close.py` | fast numerical closure; 40-digit repeat of every number |
| `s15_affine_relations.py` | exact proof of the two affine relations on `V_R` |
| `s16_axes_exact_check.py` | exact three-axes checks on asymmetric rational configurations |

### What did not close

Nothing in items 0–4 is left open. The weakest links, stated plainly, are the two entries
in § 4 "missing": the off-block vanishing of `N` rests on Schur plus a `1e−16` numerical
check rather than an exact `28×28` computation, and the spinor identity
`|q_n(z)|² = (1−(m·n)²)/4` together with the Hopf pushforward is verified (exactly at seven
configurations, numerically to `1e−16` elsewhere) rather than proved symbolically in this
room. Both are places where an error would change the answer, so they are named here rather
than absorbed.
