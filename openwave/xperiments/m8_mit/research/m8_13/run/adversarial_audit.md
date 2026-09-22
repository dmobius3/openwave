# Adversarial audit of eight claims about `r̂₆` on `V₃ = ℂ⁷`

Everything below was built from `conventions.md` alone. No number in `claims.md` was
used as an input to any computation; every claimed quantity was recomputed from the
definitions and only then compared. The room contains no `./py`, so the interpreter
used is the `python3` named in the task brief (3.12.14, sympy 1.14.0, numpy 2.5.3,
mpmath 1.3.0, scipy present).

## 0. Verdict table

| # | Claim | Verdict | Strength of the check |
| --- | --- | --- | --- |
| C1 | `max r̂₆ = 463/924` | SURVIVED | Exact proof (rational arithmetic, no floating point) |
| C2 | maximizers = rotation+phase orbit of `(v₃+v₋₃)/√2`, nowhere else | SURVIVED | Exact completeness certificate + 4000-start scan found no counterexample |
| C3 | `min r̂₆ = 1/924` on the coherent-state orbit | SURVIVED | Exact proof by a route the claims do not use (Sym² spectrum), plus the C8 route |
| C4 | `TrN̄² ≤ 171/2`, equivalently `λ_max ≤ 15/√6` | SURVIVED | Operator bound proved exactly; "equivalently" is a theorem with hypotheses, stated below |
| C5 | zero eigenvalue + three 2×2 blocks `M₁, M₂, M₃` | SURVIVED with one wording nuance | Exact symbolic block reduction; see § 5 on `M₂ = M₁ at −b` |
| C6 | `λ_max = 15/√6` at exactly three points, no fourth | SURVIVED | Exact closed-form factorization of `det(λ*I − M_k)` on the ellipse; the search is closed, not scanned |
| C7 | multiplicity exactly two, reached by exactly two of three blocks | SURVIVED | Exact spectra at all three points |
| C8 | the quartic identity | SURVIVED | Exact polynomial identity in 14 independent symbols |

Nothing was refuted. Three things the claims do not mention are in § 9, one of which
(the `span{v₃,v₋₃}` plane) is the concrete form of the trap the brief warns about.

## 1. Arming the instrument first

Before testing any claim I checked my own implementation against facts with known
answers (`s0_instrument.py`, output `out_s0.txt`). If any of these had failed, every
disagreement I reported afterwards would have been my bug, not someone else's.

| Check | Result |
| --- | --- |
| `[Jx,Jy] = iJz`, `[Jy,Jz] = iJx`, `J² = 12·I`, all `Jᵢ` Hermitian | ≤ 1.8e-15 |
| `D(n,2π) = +I` (integer spin, `SO(3)` rep as conventions 1.2 asserts) | 2.5e-15 |
| `D` unitary; `D_z(θ)v_m = e^{-imθ}v_m`; `D_z(a)D_z(b)=D_z(a+b)` | ≤ 8.9e-16 |
| `Θ` antilinear, `⟨Θu,Θw⟩ = conj⟨u,w⟩`, `ΘJΘ⁻¹ = −J`, `ΘD = DΘ` | ≤ 4.3e-15 |
| `‖Θu‖ = ‖u‖` exactly (sympy, symbolic `c₀..c₆`) | difference `= 0` |
| CG Condon-Shortley anchor `⟨3 3;3 3\|6 6⟩` | `= +1` |
| full 49×49 CG matrix over `J=0..6` is orthogonal | `max\|MMᵀ − I\| = 0` exactly |
| `Σ_{J,Q}\|ρ_J,Q\|²/‖u‖⁴` | `= 1.000000000000000` |
| `r̂₆` invariant under phase, scale, and 40 random rotations | ≤ 3.6e-16 |
| `N̄(Ru) = R N̄(u)Rᵀ` and `f(Ru) = R f(u)` for 20 random rotations | ≤ 9.9e-14 |
| `Tr N̄ = 12` for unit `u` | 5.3e-15 |

**Worklist item 0.** `⟨3 3; 3 −3 \| 6 0⟩ = √231/462 = 1/√924 = 0.03289758474798844941919642…`,
computed exactly by sympy and independently reproduced by the textbook stretched-coupling
formula `⟨j m; j −m\|2j 0⟩ = √(C(2j,j+m)C(2j,j−m)/C(4j,2j)) = √(1/C(12,6))`; the two agree
with difference exactly `0`. Its square is exactly `1/924`, which is where the denominator
in every claim comes from. `Θ(Θu) = +u` (verified symbolically, residual the zero vector),
as it must be for integer spin. Precision: exact rational/radical arithmetic, no rounding.

The sum rule deserves emphasis because it is a structural fact none of the claims state:
`r̂₆` is exactly the weight of the `J=6` channel in the coupling of `u` with `Θu`, so
`0 ≤ r̂₆ ≤ 1` identically, before any of this analysis.

## 2. C8: refuted? No, it is an exact identity

The brief suggests identities are cheap to refute. I did better than spot-checking. Treat
`c_a` and `d_a` (standing for `conj(c_a)`) as **14 independent symbols**. Every object in
C8 is the natural polarization of a bidegree-(2,2) form:

- `‖u‖² = Σ d_a c_a`
- `f_i = Σ d_a (f_i)_{ab} c_b`
- `N̄_ij = Σ d_a ({f_i,f_j}/2)_{ab} c_b` (`Re⟨u,f_if_ju⟩` *is* the anticommutator expectation, since each `f_i` is Hermitian)
- `a₀₀ = Σ_m (−1)^{3−m} c_{−m}c_m/√7`, which is **holomorphic** of degree 2
- `ρ_Q = Σ_{m₁} ⟨3m₁;3(Q−m₁)\|6Q⟩(−1)^{3−(Q−m₁)} c_{m₁} d_{m₁−Q}`

Expanding `lhs − rhs` in these 14 symbols gives **exactly `0`** (`s2_identity.py`; 72
monomials on each side). This is stronger than the claim: it holds for `d` unrelated to
`conj(c)`, hence in particular for every `u ∈ ℂ⁷`. The claim survives with no caveat.

Twelve awkward numerical spot checks (basis vectors, all-ones, repeated entries with
interleaved zeros, pure phases `1,i,−1,−i,…`, a `1e-13` entry next to a `1`, a
near-coherent state, a generic complex vector) agree to ≤ 2.1e-14 absolute on values up
to 112. No awkward point broke it.

## 3. C1: the maximum

The bound follows from C8 plus three separate facts, each established on its own:

| Input | Status |
| --- | --- |
| `\|f\|² ≥ 0` | trivial, it is a squared norm of a real 3-vector |
| `\|a₀₀\|² ≤ 1/7` for unit `u` | `Θ` antiunitary ⟹ `‖Θu‖=‖u‖` (verified **symbolically**, difference `0`), then Cauchy-Schwarz. Equality iff `Θu = λu`, `\|λ\|=1` |
| `TrN̄² ≤ 171/2` | C4, proved exactly in § 4 |

Exact rational arithmetic: `−5/231 − 0 + (7/11)(1/7) + (171/2)/198 = 463/924`, confirmed
by sympy (`sp.simplify(bound − 463/924) == 0`). Attainment: `r̂₆((v₃+v₋₃)/√2) = 463/924`
**exactly** by symbolic evaluation of the definition, and again at 50 digits by an
independent mpmath implementation:

```text
mpmath, 50 dps:  0.5010822510822510822510822510822510822511
463/924       :  0.5010822510822510822510822510822510822511
difference    :  1.34e-51     (mpmath round-off at 50 dps)
```

Corroborating search: 4000 L-BFGS runs in the full 14 real dimensions with analytic
gradient (gradient checked against central differences to 1e-10), from five start
families (Gaussian; sparse with ~65 % exact zeros; purely real; heavy-tailed with a
1e-3…1e1 dynamic range; 1-to-3-nonzero-coordinate starts). Best value found
`0.5010822510822516`, exceeding `463/924` by `4.4e-16`, i.e. round-off. **No start
exceeded the bound.** Verdict: SURVIVED.

## 4. C4: the bound, and what "equivalently" hides

The operator form is proved exactly in § 6. The equivalence is **not** a restatement; it is
a small theorem, and I state its hypotheses because the claim does not:

1. `TrN̄ = ⟨u,f²u⟩ = 12‖u‖²` (verified, 5.3e-15), so with `N̄₀ := N̄ − 4I` traceless,
   `TrN̄² = 48 + ‖N̄₀‖_F²`, and `48 + (15/√6)² = 48 + 75/2 = 171/2` exactly. The two forms
   of the claim are numerically consistent.
2. `N̄₀` is **real symmetric** (because `Re⟨u,f_if_ju⟩` is the anticommutator expectation),
   so it is diagonalized by some orthogonal `O`; if `det O = −1` use `−O`, which is in
   `SO(3)` and diagonalizes equally. This step is needed to reduce a general `e` to a
   *diagonal* one, which is all the claim's operator form covers.
3. `N̄` is rotation-covariant: `N̄(Du) = R N̄(u) Rᵀ`. Verified numerically to 9.9e-14 for
   20 random axes, and exactly for `z`-rotations (`D(z,t)†f_iD(z,t) = R(t)_{ij}f_j`,
   symbolic residual the zero matrix for all three `i`). This is the standard
   vector-operator identity (RECOGNISED from prior knowledge, and DERIVED here for the
   generating `z`-rotation).
4. Then `‖N̄₀‖_F = Tr(N̄₀ e)` with `e := N̄₀/‖N̄₀‖_F` unit traceless diagonal,
   `= ⟨u, A(e) u⟩ ≤ λ_max(A(e)) ≤ 15/√6`.

Without step 2 or 3 the claim's diagonal-only `e` would not control the off-diagonal part
of `N̄₀` and the equivalence would fail.

Independent global check: 1500 Powell maximizations of `TrN̄²` over the unit sphere.
Maximum found `85.50000000000013` (excess 1.3e-13, round-off); **no run exceeded
`85.5 + 1e-9`**. Verdict: SURVIVED.

## 5. C5: the block decomposition, and one wording nuance

Exact symbolic reduction (`s4_blocks.py`). `A = e₁f_x² + e₂f_y² + e₃f_z²` with
`e₃ = 2a/3`, `e₁ − e₂ = 4b` and `e₁+e₂+e₃ = 0` forces `e₁ = −a/3+2b`, `e₂ = −a/3−2b`, and
`\|e\|² − ((2/3)a² + 8b²) = 0` exactly, so the claimed ellipse is right.

In the parity-adapted basis (`P: v_m ↦ v_{−m}`) `A` is exactly block diagonal:

| Block | Basis | Matrix (exact) | Claim |
| --- | --- | --- | --- |
| 1×1 | `(v₂−v₋₂)/√2` | `0`, and its coupling to all six other basis vectors is `0` | zero eigenvalue ✅ |
| `M₃` | `((v₂+v₋₂)/√2, v₀)` | `[[0, 4√15 b],[4√15 b, −4a]]` | `[[0,√240 b],[√240 b,−4a]]`, and `4√15 = √240`; difference matrix exactly `0` ✅ |
| `M₁` | `((v₃+v₋₃)/√2, (v₁+v₋₁)/√2)` | `[[5a, 2√15 b],[2√15 b, −3a+12b]]` | `2√15 = √60`; difference matrix exactly `0` ✅ |
| `M₂` | `((v₃−v₋₃)/√2, (v₁−v₋₁)/√2)` | `[[5a, +2√15 b],[+2√15 b, −3a−12b]]` | see below |

**The nuance.** `M₁` evaluated at `−b` is `[[5a, −√60 b],[−√60 b, −3a−12b]]`. The natural
antisymmetric block carries `+√60 b`, not `−√60 b`. The two matrices are conjugate by
`diag(1,−1)` and their characteristic polynomials are equal (verified: difference of
charpolys simplifies to `0`), and flipping the sign of one antisymmetric basis vector
makes them literally equal. Since C6 and C7 are statements about **eigenvalues**, this
costs nothing. But "`M₂ = M₁` at `−b`" is true only up to that basis orientation, not as
an equality of matrices in the obvious basis. I record it as a wording nuance, not a
refutation.

**Convention probes** (the brief asks specifically). Reading `e₁−e₂ = 4b` as `e₂−e₁ = 4b`
is exactly `b → −b`. The ellipse is invariant under `b → −b`, and the claimed contact set
`{(√6/2,0), (−√6/4,±√6/8)}` is invariant **as a set** (it swaps the last two, and swaps the
labels `M₁ ↔ M₂`). Reading `e₃ = 2a/3` with the distinguished axis on `x` or `y` instead of
`z` is an axis relabelling, and the contact set is the full set of permutations of
`√6·e = (2,−1,−1)`, which is permutation-invariant. **Neither defensible re-reading changes
the answer to C6.** Verdict: SURVIVED.

## 6. C6: the completeness claim, closed exactly rather than scanned

This is where a scan would have been worthless, so I did not rely on one. Substituting
`a = √6 α`, `b = √6 β` turns the ellipse into `4α² + 48β² = 1`, and `λ* = 15/√6 = 5√6/2`.
For a real symmetric 2×2 block, `λ_top ≤ λ*` **iff** `tr(λ*I−M) ≥ 0` and `det(λ*I−M) ≥ 0`,
and `λ_top = λ*` iff additionally `det(λ*I−M) = 0`. Reducing each determinant modulo the
ellipse gives, exactly:

```text
det(λ*I − M₁)|ellipse = 30 (1 − 2α)(α + 1 − 6β)
det(λ*I − M₂)|ellipse = 30 (1 − 2α)(α + 1 + 6β)
det(λ*I − M₃)|ellipse = (15/2)(4α + 1)²
tr(λ*I − M₁) = √6 (5 − 2α − 12β),  tr(λ*I − M₂) = √6 (5 − 2α + 12β),  tr(λ*I − M₃) = √6 (5 + 4α)
```

On `4α²+48β²=1`, Cauchy-Schwarz gives `max(uα+vβ) = √(u²/4 + v²/48)`, computed exactly:

| linear form | max on the ellipse | consequence |
| --- | --- | --- |
| `2α` | `1` | `1 − 2α ≥ 0` |
| `6β − α` | `1` | `α + 1 − 6β ≥ 0` |
| `−6β − α` | `1` | `α + 1 + 6β ≥ 0` |
| `2α + 12β` | `2` | `tr(λ*I−M₁) ≥ 3√6 > 0` |
| `4α` | `2` | `tr(λ*I−M₃) ≥ 3√6 > 0` |

So all three determinants are non-negative and all three traces strictly positive on the
whole ellipse: `λ*I − M_k ⪰ 0` everywhere, i.e. `λ_max(A) ≤ 15/√6` **everywhere**, which is
C4's operator form, proved. The zero set is read off the factorizations and confirmed by
exact solve:

| block | contact points `(α,β)` |
| --- | --- |
| `M₁` | `(1/2, 0)` and `(−1/4, 1/8)` |
| `M₂` | `(1/2, 0)` and `(−1/4, −1/8)` |
| `M₃` | `(−1/4, 1/8)` and `(−1/4, −1/8)` |

Union: exactly three points, `(a,b) = (√6/2, 0), (−√6/4, √6/8), (−√6/4, −√6/8)`, whose
`√6·e` are `(−1,−1,2)`, `(2,−1,−1)`, `(−1,2,−1)`, precisely the three permutations of
`(2,−1,−1)`, computed, not assumed. **There is no fourth point**, and this is a certificate:
the two factors and the perfect square are non-negative on the ellipse and I have all their
zeros in closed form. A fourth point would require one of `(1−2α)`, `(α+1∓6β)`, `(4α+1)²` to
vanish elsewhere on the ellipse, which the Cauchy-Schwarz equality analysis excludes.

One trap avoided: a block could in principle have `λ*` as its *lower* eigenvalue. It cannot,
because that would force its top eigenvalue above `λ*`, contradicting `λ*I − M_k ⪰ 0`.
Verdict: SURVIVED.

## 7. C7: "exactly two", tested for exactly the failure modes named

Full 7×7 spectra at each contact point, exact:

| `(a,b)` | spectrum of `A` | top | mult | `M₁` | `M₂` | `M₃` |
| --- | --- | --- | --- | --- | --- | --- |
| `(√6/2, 0)` | `5√6/2`(×2), `0`(×2), `−3√6/2`(×2), `−2√6`(×1) | `5√6/2 = λ*` | **2** | reaches, mult 1 | reaches, mult 1 | does not |
| `(−√6/4, √6/8)` | `5√6/2`(×2), `0`(×2), `−3√6/2`(×2), `−2√6`(×1) | `λ*` | **2** | reaches, mult 1 | does not | reaches, mult 1 |
| `(−√6/4, −√6/8)` | same | `λ*` | **2** | does not | reaches, mult 1 | reaches, mult 1 |

No point has three blocks at the bound; no multiplicity other than two occurs; each
reaching block contributes multiplicity 1 (a 2×2 block could contribute 2 only if it were
`λ*I`, which none is). The `(v₂−v₋₂)/√2` vector is annihilated exactly at all three points.
Verdict: SURVIVED, including the "exactly two" that the brief flagged.

## 8. C2 and C3: the two completeness claims

### 8.1 What the scan covered, and what it cannot establish

4000 optimizations for the maximum and 4000 for the minimum, L-BFGS-B with analytic
gradient, `ftol=1e-18`, `gtol=1e-14`, each polished by two restarts, in the full 14 real
dimensions, from the five start families listed in § 3.

| Run | converged at the extremum | worst deviation | `\|f\|²` | `\|a₀₀\|²` | `TrN̄²` |
| --- | --- | --- | --- | --- | --- |
| MAX | 3534 / 4000 | 4.4e-16 | `[0, 6.5e-16]` | `0.142857142857 = 1/7` | `85.5` |
| MIN | 3658 / 4000 | 3.0e-18 | `9.000000` | `0.000000` | `85.5` |

The equality pattern C8 predicts is exactly what the numerics show, independently.

Orbit membership was tested **deterministically for every single one of those 7192 points**,
not by an inner search: for a maximizer, rotate into the frame where `N̄₀` is diagonal with
its distinct eigenvalue on `z` (trying all 24 frame choices) and measure the mass outside
`{v₃,v₋₃}` plus the `\|c₃\|²−\|c₋₃\|²` imbalance; for a minimizer, rotate `⟨f⟩` onto `+z` and
measure `1 − \|c₃\|²`. The test was armed first: it returns ≤ 1.3e-15 on 200 random
rotations-and-phases of each reference, and returns large values (0.97, 0.63, 1.00, 0.10)
on vectors known to be outside.

```text
MAX: orbit gap over ALL 3534 points: max 8.1e-09, median 2.9e-11, count(>1e-6) = 0
MIN: orbit gap over ALL 3658 points: max 2.4e-15, median ~0,     count(>1e-6) = 0
```

Local second-order check at the two references (central-difference Hessian on the sphere):
the cat state has 5 zero modes (1 radial + 4 orbit) and **all nine other eigenvalues
negative**; `v₃` has 4 zero modes (1 radial + 3 orbit) and **all ten others positive**. So
each extremum is strict transverse to its orbit, and the extremizer set is locally exactly
the orbit.

**What this could not have caught.** A basin of attraction that L-BFGS never enters from any
of those 20000 starts; an extremizer set component of measure zero in the sampling
distribution; anything at a scale below the 1e-9 clustering tolerance; and, most
importantly, it cannot establish completeness at all, because a scan that finds nothing is weak
evidence. That is why both completeness claims below are closed by certificates instead.

### 8.2 C2: exact completeness certificate

By C8 (exact) equality in `r̂₆ = 463/924` forces **all three** of `\|f\|² = 0`,
`\|a₀₀\|² = 1/7`, `TrN̄² = 171/2` simultaneously, since each is at its own extreme with a
fixed nonzero coefficient. Then:

1. `\|a₀₀\|² = 1/7` forces `Θu = λu`, `\|λ\| = 1` (Cauchy-Schwarz equality, `‖Θu‖=‖u‖` proved
   symbolically). Rephasing `u' = e^{iφ}u` gives `Θu' = e^{−2iφ}λu'`; choosing
   `e^{2iφ} = λ` gives `Θu' = u'`. `r̂₆` is phase-invariant and the claimed orbit contains
   all phases, so this is WLOG. *Omit this and the final linear algebra has no handle on
   the relative phase of `α` and `β`.*
2. `TrN̄² = 171/2` forces equality all along the chain in § 4: with `u` rotated so `N̄₀` is
   diagonal and `e := N̄₀/‖N̄₀‖_F`, both `⟨u,A(e)u⟩ = λ_max(A(e))` and
   `λ_max(A(e)) = 15/√6`. By C6 (exactly proved), `e` is one of the three contact points;
   by C7, `u` lies in the 2-dimensional top eigenspace. *Omit the "`e` is a contact point"
   half and you are left with the trap in § 9.3: saturating C4 alone permits every value
   of `r̂₆` in `[1/924, 463/924]`.*
3. The three contact `e` are one `SO(3)` orbit: the cyclic axis permutation
   `(x,y,z)→(z,x,y)` has determinant `+1`, checked. So WLOG `(a,b) = (√6/2, 0)`. *Omit
   this and two of the three branches would have to be redone verbatim.*
4. At `(√6/2, 0)`, sympy's exact `eigenvects` gives the top eigenspace as exactly
   `span{v₃, v₋₃}` (multiplicity 2). So `u = αv₃ + βv₋₃`.
5. On that plane, exactly: `⟨f_z⟩ = 3\|α\|² − 3\|β\|²`, and `⟨f_x⟩ = ⟨f_y⟩ = 0` identically
   (`f_±` shift `m` by 1, and `3` and `−3` differ by 6). So `⟨f⟩ = 0 ⟺ \|α\| = \|β\|`.
   And `Θu = (conj(β), 0,0,0,0,0, conj(α))`, so `Θu = u ⟺ α = conj(β)`. With
   `\|α\|²+\|β\|² = 1` this gives exactly `α = e^{iγ}/√2`, `β = e^{−iγ}/√2`.
   Sympy confirms `r̂₆` of that family `= 463/924` for every `γ`.
6. `D_z(θ)v_m = e^{−imθ}v_m`, so `D_z(−γ/3)(v₃+v₋₃)/√2 = (e^{iγ}v₃ + e^{−iγ}v₋₃)/√2`
   (checked numerically to 6.5e-16 for 200 random `γ`). Every such `u` is in the orbit.

Maximizer set `⊆` orbit by 1-6, and `⊇` orbit because `r̂₆` is rotation- and
phase-invariant and equals `463/924` at the reference. Hence **equality**. Verdict:
SURVIVED, with a certificate, not a scan.

The orbit is 4 real dimensions inside the 13-sphere (3 rotation + 1 phase, discrete
stabilizer: `C₃` about `z` and the `π`-rotation about `x`), matching the 5 zero Hessian
modes (4 + radial).

### 8.3 C3: exact completeness by a route the claims never use

I did not want the minimum to depend on the same C8/C4/C6 chain as the maximum, because a
shared error would survive both. So I built a second instrument.

Write `U = u ⊗ u ∈ Sym²(ℂ⁷)` (28 complex dimensions). Then
`r̂₆(u)‖u‖⁴ = ⟨U, H U⟩` for a Hermitian `H` on `Sym²`, built from `ρ_Q = ⟨u, B_Q u⟩`
(verified: `⟨U,HU⟩ − r̂₆‖u‖⁴` is ≤ 2.1e-14 for 30 random `u`; `‖H − H†‖ = 0` exactly on
`Sym²`). Since `r̂₆` is rotation-invariant and `Sym²(V₃) = V₆ ⊕ V₄ ⊕ V₂ ⊕ V₀`, Schur's
lemma forces `H` to be a scalar on each irrep. Computed **exactly** (sympy, on the
highest-weight vector of each `V_K`), and independently by numerical diagonalization with
the right multiplicities `13, 9, 5, 1`:

| irrep | `h_K` exact | `×924` | multiplicity |
| --- | --- | --- | --- |
| `V₆` | `1/924` | `1` | 13 |
| `V₄` | `13/154` | `78` | 9 |
| `V₂` | `65/84` | `715` | 5 |
| `V₀` | `13/7` | `1716` | 1 |

So **identically**, with `w_K = ‖P_K U‖²/‖U‖²` and `Σ w_K = 1`:

```text
924 · r̂₆(u) = w₆ + 78 w₄ + 715 w₂ + 1716 w₀
```

Therefore `r̂₆ ≥ 1/924` for **every** `u`, with equality **iff** `u⊗u` lies entirely in the
`V₆` component. That is a one-line exact proof of C3's value, independent of C8, C4 and C6.
(Note it does **not** prove C1: `λ_max(H) = 13/7 ≈ 1.857 ≫ 463/924`, because the top
eigenvector is not of rank one. The maximum genuinely needs the rank-1 constraint and hence
the other route. Reporting this asymmetry rather than smoothing it.)

The equality set: `u⊗u ∈ V₆` means the 15 quadrics spanning `(V₀⊕V₂⊕V₄)` all vanish at `u`.

1. Parametrize the perfect-6th-power cone by `c_m = √C(6,3+m)·x^{3+m}y^{3−m}`. Verified
   numerically that this reproduces the rotation orbit of `v₃` (overlap deficit 2.2e-16 for
   50 random directions).
2. **Exact symbolic check:** all 15 of those quadrics vanish identically in `(x,y)` on the
   cone; all 13 of the `V₆` quadrics do not.
3. In the coordinates `a_k = c_{k−3}/√C(6,k)`, the restriction of the quadric `a_ia_j` to
   the cone is `x^{i+j}y^{12−i−j}`, and `i+j` runs over `0..12`, so restriction onto
   degree-12 binary forms is **onto**; its kernel has dimension `28 − 13 = 15`. The 15
   quadrics are independent (`V₀⊕V₂⊕V₄` is 15-dimensional and `Sym² → quadrics` is an
   isomorphism; numerical rank confirmed `15`), lie in that kernel, hence **span** it.
4. Every `2×2` Hankel minor `a_ia_j − a_ka_l` (`i+j = k+l`) vanishes on the cone, verified
   symbolically for all 91 of them, so each lies in that kernel, hence in the span of the
   15. Numerically the combined rank of `[15 quadrics ; 91 minors]` is exactly `15`, so the
   two spans coincide.
5. Rank-≤1 Hankel is elementary: if `a₀ ≠ 0` set `t = a₁/a₀`, then `a_k a₀ = a_{k−1}a₁`
   gives `a_k = a₀t^k`; if `a₀ = 0` then `a₁² = a₀a₂ = 0` forces `a₁ = 0`, and inductively
   `a = (0,…,0,a₆)`. Both are perfect 6th powers `l⁶`.
6. `SU(2)` is transitive on `P¹`, so every `l⁶` is a scalar times a rotated `v₃`. Hence
   `{u : u⊗u ∈ V₆}` is exactly the cone over the coherent states, and the unit vectors in
   it are exactly the rotation-and-phase orbit of `v₃`.

Independent numerical confirmation: 600 least-squares solves of `P_{V₀⊕V₂⊕V₄}(u⊗u) = 0`
from random starts found 600 solutions, **all** within 1.6e-15 of the coherent orbit.
Verdict: SURVIVED, complete.

Step 5 and the classical statement that the rational normal curve is cut out by the 2×2
minors of its Hankel matrix are the only places I lean on RECOGNISED prior knowledge, and I
derived the specific instances (steps 2-4) by computation rather than citing them.

## 9. Things the claims do not mention

### 9.1 `r̂₆` is a probability weight

`Σ_J ‖ρ_J‖²/‖u‖⁴ = 1` (verified to 1e-15). So `r̂₆` is the `J = 6` share of a normalized
multipole decomposition and lies in `[0,1]` before any argument. The interesting content is
that the achievable range is only `[1/924, 463/924] ≈ [0.00108, 0.50108]`.

### 9.2 The full critical spectrum

The scan converged onto exactly seven critical levels, all with tangential gradient below
2e-8 (so genuine critical points, not solver stalls), identified as exact rationals and
confirmed by exact symbolic evaluation at simple vectors:

| value | `×924` | attained at | `\|f\|²` | `TrN̄²` | `\|a₀₀\|²` |
| --- | --- | --- | --- | --- | --- |
| `1/924` | 1 | `v₃` (the minimum) | 9 | 85.5 | 0 |
| `3/77` | 36 | `v₂` | 4 | 48 | 0 |
| `75/308` | 225 | `v₁` | 1 | 61.5 | 0 |
| `9/35` | **237.6, not an integer** | a saddle | 4/25 | 1416/25 | 0 |
| `24/77` | 288 | `(v₂+v₋₂)/√2` | 0 | 48 | 1/7 |
| `100/231` | 400 | `v₀` | 0 | 72 | 1/7 |
| `463/924` | 463 | `(v₃+v₋₃)/√2` (the maximum) | 0 | 85.5 | 1/7 |

`9/35` is the one critical value that is **not** of the form `n/924`; C8 reproduces it
exactly from `−5/231 − (4/25)/22 + (1416/25)/198 = 9/35`. Also `r̂₆((v₁+v₋₁)/√2) = 145/308
= 435/924`, which is below the maximum.

### 9.3 The trap, made concrete: the plane `span{v₃, v₋₃}`

This is the sharpest thing I found, and no claim mentions it. For any `α, β` (exact
symbolic result, with `α = p+iq`, `β = r+is` and the residual simplifying to `0`):

```text
unit u = α v₃ + β v₋₃   ⟹   r̂₆(u) = 2|α|²|β|² + 1/924        and        Tr N̄² = 171/2
```

`N̄` on that whole plane is exactly `diag(3/2, 3/2, 9)`, independent of `α, β`. So **every
state in a 2-complex-dimensional plane saturates C4's bound**, while `r̂₆` sweeps the entire
achievable interval `[1/924, 463/924]` across the same plane: `α=1, β=0` gives the global
minimum, `α=β` and `α=1, β=i` both give the global maximum, and `α=1, β=1/2` gives
`7417/23100`, strictly between. Counting the rotation orbit of the plane, the
`TrN̄² = 171/2` set is at least 5-dimensional, and a direct scan confirms it: 1496 of 1500
random maximizations of `TrN̄²` land on `85.5`, and the `r̂₆` values there spread over the
full interval (44 distinct values at 9-decimal rounding).

This is precisely the confusion the brief warns about. Saturating C4 is necessary for both
extrema and sufficient for neither. Any argument that stops at "`TrN̄²` is maximal, therefore
maximal `r̂₆`" is wrong, and the `⟨f⟩ = 0` and `Θu = λu` conditions are doing the real work
in C2.

## 10. Tools, and what was looked up rather than derived

| Item | How obtained |
| --- | --- |
| Clebsch-Gordan values | `sympy.physics.wigner.clebsch_gordan`, Condon-Shortley. Anchor `⟨33;33\|66⟩ = +1` and full 49×49 orthogonality verified exactly, so the library's convention is confirmed, not assumed |
| Stretched-coupling closed form `⟨j m;j −m\|2j 0⟩` | RECOGNISED (standard angular-momentum identity), then DERIVED-equal to the library value, difference exactly `0` |
| `Θ² = +1` for integer spin | DERIVED symbolically from the given `Θ` (not assumed from the standard `(−1)^{2j}` result) |
| Vector-operator relation `D†f_iD = R_{ij}f_j` | RECOGNISED, and DERIVED exactly for `z`-rotations plus verified numerically for 20 random axes |
| `Sym²(V₃) = V₆⊕V₄⊕V₂⊕V₀` and Schur's lemma | RECOGNISED (representation theory); the four eigenvalues and multiplicities were DERIVED by computation |
| Rational normal curve cut out by `2×2` Hankel minors; binary form is a perfect power iff its Hessian vanishes | RECOGNISED; the specific instances needed (15 quadrics vanish on the cone; 91 minors vanish; spans coincide; surjectivity of the restriction) were all DERIVED here |
| `SU(2)` transitive on `P¹` | RECOGNISED, elementary |
| Optimizers | `scipy.optimize` L-BFGS-B (analytic gradient, checked to 1e-10 against central differences), Powell, `least_squares` |
| High precision | mpmath at 50 dps, independent re-implementation of `r̂₆` |
| Exact algebra | sympy rational + radical arithmetic throughout; no floating point in any certificate |

Scripts, all in this directory, all using `pathlib.Path(__file__).parent`:
`core.py`, `grad.py`, `s0_instrument.py`, `s1_scan.py`, `s2_identity.py`, `s3_sos.py`,
`s4_blocks.py`, `s5_orbits.py`, `s6_certificate.py`, `s7_extras.py`, `s8_exactsteps.py`,
`s9_fix.py`, `s10_saddles.py`; outputs `out_s0.txt` … `out_s9.txt`.

## 11. Honest residue

- Two sympy `simplify` calls initially returned non-zero residues (the `span{v₃,v₋₃}`
  formula, and the `f_y` case of the vector-operator relation). Both were **simplifier
  failures**, not defects: re-running with `α,β` split into real and imaginary parts, and
  with a `rewrite(cos)` before expansion, gave exactly `0` in both cases. I report this
  because the intermediate state looked exactly like a disagreement with C5/C4 and was not.
- The rotation-covariance of `N̄`, needed for C4's "equivalently", is proved exactly only
  for `z`-rotations here; for general axes it rests on numerical verification to 1e-13 plus
  the recognised vector-operator identity. This is the weakest link in my chain and I flag
  it rather than claim more.
- Three rank computations in § 8.3 (`15`, `15`, `15`) are numerical SVD ranks at tolerance
  1e-9. The conclusion does not actually depend on the third of them, because the dimension
  count in step 3 is exact and forces the spans to coincide; but the numbers themselves are
  numerical.
- No claim was refuted. Given that the brief asked me to break them, I want to be explicit
  that I looked for breaks in the places named (a fourth ellipse point; three blocks at the
  bound; the identity at awkward points; the convention re-readings; off-orbit extremizers
  in 14 real dimensions) and found none, and that for C2, C3 and C6 the negative result is
  backed by a certificate rather than by an exhausted search.
