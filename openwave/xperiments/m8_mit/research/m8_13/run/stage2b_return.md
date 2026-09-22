# Return, stage 2b: grading `S0_S3_MAXIMUM.md`

Ten verdicts, numbered as `stage2b_grading.md` § "What to grade" numbers them, then the
reconcile. Every number below comes from one of the six scripts listed in § 8; none is
read off the author's files and none is arithmetic done in my head. `check_s3_maximum.py`
was read as text and never executed.

## 0. Verdict table

| # | Part | Verdict |
| --- | --- | --- |
| 1 | Step 1, the invariant form | **ESTABLISHED, SUPPLIED** |
| 2 | Step 2, the two bounds and their equality cases | **ESTABLISHED** |
| 3a | Reduction to `λ_max(e₁fₓ²+e₂f_y²+e₃f_z²)`, and the split | **ESTABLISHED, SUPPLIED** |
| 3b(`M₁`) | The `M₁` bound | **ESTABLISHED, SUPPLIED** |
| 3b(`M₂`) | The `M₂` bound | **ESTABLISHED** |
| 3b(`M₃`) | The `M₃` bound | **ESTABLISHED, SUPPLIED** |
| 3c | The equality set | **ESTABLISHED, SUPPLIED** |
| 3d | From `‖Q‖ = 15/√6` to the top eigenspace | **ESTABLISHED, SUPPLIED** |
| 3e | That eigenspace as `span{\|3,3⟩ₙ, \|3,−3⟩ₙ}` | **ESTABLISHED, SUPPLIED** |
| 4 | Step 4, the combination | **ESTABLISHED, SUPPLIED** |

No GAP and no DEFECT. Every claim in the note that I could check came out true; what is
repeatedly missing is derivation, not correctness. The argument is sound and the omissions
are all one-liners or short computations, which is why nothing is graded GAP; and I found
no false statement and no hypothesis that could fail, which is why nothing is graded
DEFECT. § 7 records the three places where I considered DEFECT and why I did not take it.

---

## 1. Step 1, the invariant form — ESTABLISHED, SUPPLIED

**Verified.** Everything the note asserts here is correct.

| Claim | My computation |
| --- | --- |
| `‖ρ₆‖²` at `v₃,v₂,v₁,v₀` is `1, 36, 225, 400` over `924` | `1/924, 3/77 = 36/924, 75/308 = 225/924, 100/231 = 400/924` |
| `\|f\|²` there | `9, 4, 1, 0` |
| `\|a₀₀\|²` there | `0, 0, 0, 1/7` |
| `TrN̄²` there | `171/2, 48, 123/2, 72` (and `N̄` diagonal at each, off-diagonals exactly `0`) |
| independence determinant | `180/7`, exactly as the note states |
| solved coefficients | `[−5/231, −1/22, 7/11, 1/198]`, exactly as the note states |
| the identity itself | residual `9.643e−18` over 12 fresh random unit states in 30-digit arithmetic |

The stretched CG closed form I used to build `ρ₆` was checked against sympy's `CG` at all
49 pairs: zero mismatches.

**Supplied.**

| # | What I supplied | Why it is needed |
| --- | --- | --- |
| S1.1 | The bridge from "invariant quartic" to "Hermitian form on `Sym²`", **with injectivity**. The note's dimension count only gives `dim = 4` once one knows that bidegree-`(2,2)` real forms correspond bijectively to Hermitian forms on `Sym²V₃`. I supplied the reason: the `784` functions `conj(m_α)m_β` built from the `28` degree-2 monomials carry pairwise distinct (antiholomorphic, holomorphic) bidegrees (computed: `28` monomials, `784` ordered pairs, `0` duplicates), so they are linearly independent, so the map is injective. Dimension counting alone does **not** give injectivity. | Without it, "the invariant quartics form a four-dimensional space" does not follow from multiplicity-freeness. |
| S1.2 | That `‖ρ₆‖²` is itself a `U(1)×SU(2)`-invariant quartic, i.e. that it lies in the span at all. The note goes straight to "the weight-state values then force". | Without it the four weight-state equations determine nothing. |
| S1.3 | That `U(1)`-invariance kills every bidegree except `(2,2)` (a monomial `u^α ū^β` picks up `e^{i(\|α\|−\|β\|)φ}`). | This is what restricts "quartic" to the `Sym²` picture in the first place. |

---

## 2. Step 2, the two bounds and their equality cases — ESTABLISHED

Both sentences are complete and correct as written.

- `|f|² ≥ 0` with equality iff `⟨f⟩ = 0`: immediate, it is a squared norm.
- `|a₀₀|² ≤ 1/7` with equality iff `Θu ∝ u`: Cauchy–Schwarz on `⟨Θu, u⟩` with
  `‖Θu‖ = ‖u‖`. I verified `‖Θu‖² − ‖u‖² ≡ 0` symbolically on generic `u`, and `Θ² = id`.

Two compressions I record without upgrading the verdict, because neither is a step the
reader has to fill in to reach the conclusion:

| Compression | Status |
| --- | --- |
| `‖Θu‖ = ‖u‖` is used and never stated. It is a property of the object `Θ` (immediate from its definition), not an extra hypothesis, so it is not a DEFECT; I verified it. | verified |
| "`Θu ∝ u`, that is, iff `u` is time-reversal invariant" identifies `Θu = λu` with `Θu = u`. These agree **as rays**, not as vectors: if `Θu = λu` with `\|λ\| = 1` then `Θ(e^{iφ}u) = e^{iφ}u` exactly when `e^{2iφ} = λ`, which is always solvable. The whole note works modulo phase ("under rotations and phase"), so the gloss is legitimate. | verified |

---

## 3a. The reduction and the split — ESTABLISHED, SUPPLIED

**Verified.**

| Claim | My computation |
| --- | --- |
| `Tr(QE) = ⟨u, A_E u⟩` for unit traceless symmetric `E` | max error `2.08e−15` over 200 random (state, `E`) draws |
| `e₃ = 2a/3`, `e₁−e₂ = 4b` with `e` traceless; the operator identity `e₁fₓ²+e₂f_y²+e₃f_z² = a(f_z²−4)+b(f₊²+f₋²)` | exact, `True` |
| the constraint `‖e‖² = (2/3)a² + 8b²` | exact, `True` |
| the split, in the parity/flip adapted basis | exact: all inter-block entries `0`; `M₁ = [[5a, 2√15 b],[2√15 b, −3a+12b]]` (and `2√15 = √60`), `M₃ = [[0, 4√15 b],[4√15 b, −4a]]` (and `4√15 = √240`), seventh diagonal entry exactly `0` |

**Supplied.**

| # | Supplied | Note |
| --- | --- | --- |
| S3a.1 | **`TrN̄ = 12`**, which the sentence "Write `N̄ = 4I + Q`, with `Q` traceless. Then `TrN̄² = 48 + ‖Q‖²`" silently requires. I verified `TrN̄ − 12‖u‖² ≡ 0` symbolically on generic `u`. The note never states it; `check_s3_maximum.py` has **no gate** for it (the log's check of it, line 31, is from a different script, see § 6). | hinge of the whole step |
| S3a.2 | **`‖Q‖ = max_E Tr(QE)`**: `E = Q/‖Q‖` is unit, symmetric and traceless and achieves it; Cauchy–Schwarz in the 5-dimensional space of traceless symmetric matrices gives `≤`. The note states the identity without either half. | |
| S3a.3 | **The equivariance `A_{RER^T} = D(R) A_E D(R)†`**, which is what licenses "rotating `E` to diagonal form does not change `λ_max(A_E)`". Verified: max error `1.91e−14` over 50 random rotations and `E`. The note asserts the consequence only. | |
| S3a.4 | **The max-swap** `max_ρ max_E tr(ρA_E) = max_E max_ρ tr(ρA_E) = max_E λ_max(A_E)`, valid because both ranges are compact and the pairing is continuous. Unstated. | |
| S3a.5 | The block entries themselves. The note exhibits `M₁, M₂, M₃` with no derivation; I rebuilt them in an explicit adapted basis. | |

One observation on wording, not graded. The note writes "`M₂` = `M₁` at `−b`". As the block
actually arises in the natural basis `(t₃, t₁) = ((\|3⟩−\|−3⟩)/√2, (\|1⟩−\|−1⟩)/√2)` it is
`[[5a, +2√15 b],[2√15 b, −3a−12b]]`, which is **not** literally `M₁(−b) = [[5a, −2√15 b],
[−2√15 b, −3a−12b]]`. Flipping the sign of the basis vector `t₁` makes it literally so, and
the two have identical characteristic polynomials (verified exactly), so the note's
statement is realisable and nothing downstream depends on the sign. Worth knowing because
`check_s3_maximum.py` line 64 defines `M2` with the **opposite** sign convention from the
note's sentence and then gates the charpoly equality (line 82), which is the careful thing
to do.

---

## 3b. Each block's bound

### `M₁` — ESTABLISHED, SUPPLIED

**Verified.** `λ_max(M₁) = a + 6b + √(16a² − 48ab + 96b²)`: `tr/2 = a+6b` exactly and
`(tr/2)² − det − (16a²−48ab+96b²) = 0` exactly, so it is the larger root by the quadratic
formula; numerically `max |λ_max − formula| = 1.421e−14` over 20000 random `(a,b)`.
The three identities the argument turns on are exact:
`(a²+6ab+24b²)² − (3/2)(a+6b)²N² = 36b²(a+2b)²`;
`a²+6ab+24b² = (a+3b)² + 15b²`;
`K²N² = 25a² + 300b²` with `K = 15/√6`.

**Supplied.** The note says "After squaring, and homogenizing with `(15/√6)²N² = 25a²+300b²`,
the bound reads `(3/√6)(a+6b)·N ≤ a²+6ab+24b²`" and shows no intermediate. I wrote that
algebra out: squaring `K − a − 6b ≥ √(16a²−48ab+96b²)` gives
`K² − 2K(a+6b) − 15a² + 60ab − 60b² ≥ 0`; replacing `K² → K²N² = 25a²+300b²` and
`K → KN` makes it homogeneous of degree 2, giving `10a² + 60ab + 240b² ≥ 2KN(a+6b)`, and
dividing by 10 with `K/5 = 3/√6` gives the stated form. I also supplied the two
reversibility conditions for the second squaring (both sides nonnegative on the branch
`a+6b > 0`), and that the homogeneous inequality holds for **all** `(a,b)`, not only on the
ellipse, which is what the sum of squares actually delivers.

I also record the garbled sentence: "The right side of `15/√6 − a − 6b ≥ √(…)` is positive,
since `a + 6b ≤ √6` on the ellipse" calls the **left** side the right side. Each clause is
individually true (a square root is nonnegative too), and the justification given is
exactly the one the left side needs, so this is a wording slip and not a false statement.

### `M₂` — ESTABLISHED

"It has `M₁`'s spectrum at `−b`, and the constraint is even in `b`." That is a complete
inference: `λ_max(M₂(a,b)) = λ_max(M₁(a,−b)) ≤ K` because `(a,−b)` is on the ellipse
whenever `(a,b)` is. I verified the spectral half exactly (identical characteristic
polynomials) and the evenness is immediate from `(2/3)a² + 8b² = 1`. Nothing had to be
written that the note does not contain. The verdict is about `M₂`'s own two clauses; the
`M₁` bound it leans on carries its own verdict above.

### `M₃` — ESTABLISHED, SUPPLIED

**Verified.** In general `λ_max(M₃) = −2a + √(4a² + 240b²)` (same quadratic-formula check;
numerically `1.421e−14` over 20000 random `(a,b)`), and the square completion
`(K + 2a)² − (30 − 16a²) = 20(a + √6/4)²` is exact.

**Supplied.**

| # | Supplied |
| --- | --- |
| S3b3.1 | That the note's `λ_max = −2a + √(30 − 16a²)` is the **on-ellipse** form: `240b² = 30 − 20a²` there, verified exactly. The note writes the reduced formula without saying the ellipse was used. |
| S3b3.2 | `max\|a\| = √(3/2) = 1.22474487139159` on the ellipse, which is what makes `K + 2a ≥ K − 2√(3/2) = 3.67423461417477 = 9/√6 > 0`. The note asserts the inequality; the bound on `\|a\|` is unstated. |
| S3b3.3 | `30 − 16a² ≥ 30 − 24 = 6 > 0` on the ellipse, so the square root is real. Unstated. |

---

## 3c. The equality set — ESTABLISHED, SUPPLIED

**This is the part where this text differs from the excerpt I graded in stage 2a.** The note
here reads "Tracing the equality cases (`b = 0` with `a > 0`, `a = −2b` with `b > 0`, **and
`a = −√6/4` in `M₃`**) shows equality only at the three permutations of `(2,−1,−1)/√6`". The
stage-2a excerpt carried only the first two loci. See § 6 of the reconcile; I flag it here
because it changes the finding.

**Whether every block's equality points are found, and how I know.**

Computed, by an exact decision procedure rather than by inspection: for a real symmetric
`2×2` block, `λ_max(M) = K` is *equivalent* to `det(M − KI) = 0` **and** `tr M ≤ 2K`. That
turns the question into intersecting one conic with the ellipse, a zero-dimensional system.
In the rationalising coordinates `a = √6 u`, `b = √6 v` (ellipse `4u² + 48v² = 1`), reducing
each conic modulo the ellipse eliminates `v²` and leaves a polynomial that factors
completely over `Q`:

| Block | `det(M − KI)` reduced mod the ellipse | Equality points `(u,v)` |
| --- | --- | --- |
| `M₁` | `−30(2u − 1)(u − 6v + 1)` | `(1/2, 0)`, `(−1/4, 1/8)` |
| `M₂` | `−30(2u − 1)(u + 6v + 1)` | `(1/2, 0)`, `(−1/4, −1/8)` |
| `M₃` | `15(4u + 1)²` (a **double** root: `M₃` is tangent to the bound) | `(−1/4, 1/8)`, `(−1/4, −1/8)` |
| seventh eigenvalue `≡ 0` | — | none, since `K = 6.123724356957945 > 0` |

Each factor is linear, so it meets the ellipse in at most two points, and `factor_list` over
`Q` is a complete factorisation: there is nowhere for a root to hide. That is how I know the
list is complete; it does not rest on a scan. The union is **3** distinct points, carrying
`e = (−1,−1,2)/√6`, `(2,−1,−1)/√6`, `(−1,2,−1)/√6` — the three permutations, as the note says.

Answer to the question: **the three loci the note names do deliver all three points** (I
computed this: `b=0,a>0 → (1/2,0)`; `a=−2b,b>0 → (−1/4,1/8)`; `a=−√6/4` meets the ellipse at
**both** `(−1/4,±1/8)`, the second of which is the third permutation, and nothing in the
union is missed). But **`M₂`'s equality locus is never found**: it is the one block whose
equality points the note does not trace. Its points happen to lie inside the union, so the
conclusion survives, but the note gives no argument for that.

**Supplied.**

| # | Supplied | What would fail without it |
| --- | --- | --- |
| S3c.1 | `eq(M₂) ⊆ {the three points}`, via the reflection already stated in the `M₂` bullet: `eq(M₂)` is the `b ↦ −b` image of `eq(M₁)`, namely `{(1/2,0), (−1/4,−1/8)}`. Independently confirmed by the factorisation `−30(2u−1)(u+6v+1)`. | The word "only" is an upper bound on a **union over blocks**; leaving one block untraced leaves open that it contributes a fourth point, which would make the claim false. |
| S3c.2 | Derivation of the sign conditions `a > 0` and `b > 0`, which the note asserts bare. They are in fact **forced**: `2u − 1 = 0` gives `v = 0` and `u = 1/2 > 0`; `u − 6v + 1 = 0` meets the ellipse only at the double root `v = 1/8 > 0`. | Otherwise one must check the discarded signs separately; I did, and they are strictly below: `λ_max ∈ {3.674234614174767, 4.898979485566356}` against `K = 6.123724356957945`. |
| S3c.3 | The equality analysis of the branch `a + 6b ≤ 0` in `M₁`, which the note handles only for the **bound**. Equality there would need `a+6b = 0` and `(a+3b)²+15b² = 0` at once, hence `a = b = 0`, which is off the ellipse. | Without it, `eq(M₁)` is not shown to be confined to the two named loci. |
| S3c.4 | That `eq(M₃)` follows from the square completion (`20(a+√6/4)² = 0 ⟺ a = −√6/4`), and that this locus meets the ellipse at two points. The note names the locus but never says it vanishes anywhere on the ellipse, i.e. never says `M₃` **attains** the bound. | The third permutation comes from exactly there. |

---

## 3d. From `‖Q‖ = 15/√6` to the top eigenspace — ESTABLISHED, SUPPLIED

The note states the conclusion ("Equality holds iff `u` lies in the top eigenspace of
`(3(n·f)²−12)/√6` for some `n`") and no reasoning at all. **Supplied, in full:**

Take `E* = Q/‖Q‖`, unit traceless symmetric. Then `tr(ρ A_{E*}) = Tr(Q E*) = ‖Q‖ = 15/√6`.
By 3a–3c, `λ_max(A_E) ≤ 15/√6` for every unit traceless `E`. Hence
`15/√6 = tr(ρA_{E*}) ≤ λ_max(A_{E*}) ≤ 15/√6`, so both are equalities: a pure state whose
expectation equals the top eigenvalue lies in the top eigenspace, and `E*` is itself an
equality point, i.e. `E* = (3nnᵀ − I)/√6` for some unit `n`.

Verified on the hexagon: `N̄ = diag(1.5, 1.5, 9)`, `‖Q‖ = 6.123724356957944` against
`15/√6 = 6.123724356957946`; `E*` has eigenvalues `(−1,−1,2)/√6`;
`λ_max(A_{E*}) = 6.123724356957944` and `⟨u|A_{E*}|u⟩ = 6.123724356957942`; the weight of
`u` in the top eigenspace of `A_{E*}` is `1.000000000000000`. I also verified
`A_E = (3(n·f)² − 12)/√6` for `E = (3nnᵀ−I)/√6`, max error `2.22e−15` over 40 random `n`,
and that this `E` is traceless and of unit norm.

`check_s3_maximum.py` has **no gate** for this step.

---

## 3e. The eigenspace as `span{|3,3⟩ₙ, |3,−3⟩ₙ}` — ESTABLISHED, SUPPLIED

Correct, and a one-line computation the note does not perform. **Supplied:** for `n = ẑ`
the operator `(3f_z² − 12)/√6` is diagonal with entries `(3m² − 12)/√6`, computed as

```text
m  =   3        2    1          0        -1         -2    -3
val = 5√6/2     0   -3√6/2   -2√6     -3√6/2        0    5√6/2
```

so the spectrum is `{5√6/2 (×2), 0 (×2), −3√6/2 (×2), −2√6 (×1)}`, the top eigenvalue is
`5√6/2 = 15/√6` (verified exactly), and `m = ±3` are the only `m` with `3m² − 12 = 15`.
The top eigenspace is therefore exactly `span{|3,3⟩, |3,−3⟩}`, two-dimensional. A general
`n` follows by rotating. No gate covers this either.

---

## 4. Step 4, the combination — ESTABLISHED, SUPPLIED

**Verified.** Every sub-claim is true.

| Claim | My computation |
| --- | --- |
| the sign pattern that makes "all three bounds at once" the right reading | coefficient of `\|f\|²` is `−1/22 < 0` (so the bound uses `\|f\|² ≥ 0`), of `\|a₀₀\|²` is `+7/11 > 0`, of `TrN̄²` is `+1/198 > 0`; all three nonzero, so equality in the sum forces equality in each |
| `−5/231 + 1/11 + 171/396 = 463/924` | exact, `True` |
| `48 + 225/6 = 171/2`, `(15/√6)² = 75/2` | exact, `True` |
| `⟨f⟩ = 3(\|α\|²−\|β\|²)n` | exact: `⟨f⟩ = (0, 0, 3r_a² − 3r_b²)` for `u = r_a e^{iφ₁}\|3,3⟩ + r_b e^{iφ₂}\|3,−3⟩` |
| "`\|a₀₀\|² = 1/7` follows" | exact: at `\|α\|=\|β\|=1/√2`, `⟨Θu,u⟩ = e^{i(φ₁+φ₂)}`, so `\|a₀₀\|² = 1/7` for **every** pair of phases |
| the relative phase is a rotation about `n` | rotation by `χ` about `n` sends `\|3,m⟩ → e^{−imχ}\|3,m⟩`, shifting the relative phase by `−6χ`; `χ ∈ [0, π/3)` covers the full circle |
| the value on the family | `‖ρ₆‖² = 463/924` on the whole `\|α\|=\|β\|` family, all phases, exactly |

**Supplied.** The note asserts the `⟨f⟩` formula, the "`|a₀₀|² = 1/7` follows", and the
phase-as-rotation claim with no computation; all three are above. I also supplied the
strictness the step needs but does not state: on the span, exactly

```text
r6(s) = 463/924 - cos^2(2s)/2      for u = cos(s)|3,3> + sin(s)|3,-3>
      = 463/924 - (|alpha|^2 - |beta|^2)^2 / 2
```

(verified exactly; e.g. `s = π/8 → 58/231`, `s = π/6 → 695/1848`, `s = π/4 → 463/924`,
`s = 0 → 1/924`), so `‖ρ₆‖² < 463/924` **strictly** unless `|α| = |β|`. On that same span
`TrN̄² = 171/2` identically and `|f|² = 9cos²(2s)`, which is why `⟨f⟩ = 0` is the only
condition doing work there. The checker tests the strictness on 8 sampled points
(`fam`, line 117) rather than deriving it.

---

## 5. What `check_s3_maximum.py` does **not** check

Read as text, never run. Computed rather than counted by eye: the file has **27** `gate()`
call sites (`L1` 3, `L2` 17, `L3` 4, `L4` 3); one of the `L1` sites sits in a four-iteration
loop, which is how 27 sites produce the **30** executed checks the note claims and the
**30** `PASS` lines the log carries from `(L1)` onward.

**Not gated at all** (I supplied each of these in §§ 1–4):

| Unchecked | Where it is needed |
| --- | --- |
| `TrN̄ = 12` | the whole of `N̄ = 4I + Q`, `TrN̄² = 48 + ‖Q‖²` |
| `‖Q‖ = max_E Tr(QE)` | the variational rewriting, 3a |
| the equivariance `A_{RER^T} = D A_E D†` | "rotating `E` to diagonal does not change `λ_max`", 3a |
| the max-swap over `ρ` and `E` | 3a |
| the **equality set** as a computation | 3c: no gate computes any block's equality locus |
| the passage of 3d | 3d |
| the identification of the top eigenspace | 3e |
| `TrN̄² = 171/2 ⟹ u ∈ span` (the `⊆` direction) | Step 4; gate `L3`/line 116 checks only the `⊇` direction |

**Gated more weakly than the label reads:**

- line 86 is labelled "`K − a − 6b > 0` on the ellipse: `max(a+6b) = sqrt((3/2) + 36/8) = √6`"
  but asserts only the arithmetic `3/2 + 36/8 == 6` and `√6 < K`. That `√6` **is** the
  maximum on the ellipse is not tested. I computed it: `max(a+6b) = √6 = 2.44948974278318`,
  and `K − √6 = 3.67423461417477 > 0`.
- line 87 is labelled "`K + 2a > 0` on the ellipse" but asserts only the constant identity
  `K − 2√(3/2) = 9/√6`. That `|a| ≤ √(3/2)` on the ellipse is not tested. I computed it.
- line 106, the strictness gate, **excludes** every sweep point within `0.01` rad of an axis.
  Computed: grid spacing `3.141593e−04` rad, so `192` of the `20001` points (about `31.8`
  grid points per axis per side) are excluded. Inside those bands the sweep asserts nothing
  beyond the axis value itself. The equality claim of 3c therefore rests entirely on the one
  parenthetical sentence, with no exact gate behind it.

**Auditability limits, both by the coordinator's design:**

- The checker imports `rho_sq, vec, ms, idx, j` from `translate_rho6` and `inv, F3, jz, jp`
  from `check_extrema`, neither supplied. All four `L3` gates and the four weight-state
  `L1` gates run entirely on those unsupplied objects. I recomputed every value they report
  from the room's own conventions and they agree, so nothing is wrong, but those gates
  cannot be audited from what was handed over.
- The log holds `55` `PASS` lines under ten section headers; only the `30` under
  `(L1)`–`(L4)` come from this script (its `print` statements are exactly those four).
  Headers `(A)`–`(F)` are from other, unsupplied runs. So the log's line-31 check
  "`TrN = f(f+1) = 12`" is **not** evidence that the supplied checker tests it; it does not.
- Line 61, `Bm = rhs.subs({e3: 2*a/3}).subs({e1 - e2: 4*b})`, is inert: `e1 - e2` is not an
  atom, so that substitution cannot fire. Line 62 overwrites `Bm` correctly, so the result
  is right, but a reader could mistake line 61 for a check that the substitution works.

The note's line 3, "Every step is checked exactly, or armed numerically, by
`check_s3_maximum.py`", overstates the coverage on the evidence above. The count
("30 checks, 0 failures") is accurate.

---

## 6. Reconcile with my stage-1 answer

### The short answer: they agree, exactly.

| | Stage 1 (my route) | This note |
| --- | --- | --- |
| maximum | `463/924` | `463/924` |
| maximiser set | `{e^{iφ} D³(g)(v₃+v₋₃)/√2}`, Majorana constellation a regular hexagon on a great circle | `u = (e^{iφ₁}\|3,3⟩ₙ + e^{iφ₂}\|3,−3⟩ₙ)/√2` over `n ∈ S²` and the phases |
| minimum | `1/924`, coherent states | `1/924`, coherent states (asserted, derived elsewhere) |

The two descriptions are the same set: the relative phase is absorbed by a rotation about
`n` (verified above), leaving `n ∈ S²` plus a global phase, which is the rotation-and-phase
orbit of the hexagon. Computed: over 200 random rotations and phases of `(v₃+v₋₃)/√2`,
`max |‖ρ₆‖² − 463/924| = 2.11e−15`, with `‖ρ₆‖²` evaluated through the author's Step-1
identity built from my own operators. **No disagreement in either direction.**

### Two independent cross-links I did not have in stage 1

1. **My `n_J` reappear as the author's pair-channel constants.** I recomputed them here by a
   route independent of stage 1's Schur argument: `p_J(v_m) = |⟨3m;3m|J,2m⟩|²` gives a
   triangular `4×4` system on `v₃,v₂,v₁,v₀` (row sums all `1`, verified), and solving it
   against the four `‖ρ₆‖²` values yields `n_J = [13/7, 65/84, 13/154, 1/924]` — exactly my
   stage-1 constants, and exactly the `c_F` the author's log reports.
2. **The author's Step 3 and my stage-1 Step 3 are the same optimisation.** On the real form
   `{Θu = u}` I verified exactly that `|f|² ≡ 0` and that
   `TrN̄² = 48‖u‖⁴ + 126·p₂`, where `p₂` is the quadrupole weight I maximised in stage 1.
   So "maximise `TrN̄²`" and "maximise `p₂`" are literally the same problem, with
   `‖Q‖² = 126 p₂`; at the hexagon `126 · 25/84 = 75/2 = (15/√6)²`. The two routes reach the
   same bottleneck by different roads and agree on its value.

### Where this text shows my stage-1 reasoning was weaker than I thought

**One substantive item, and I did not flag it in stage 1.** My stage-1 Step 1 argued that
`r̂₆(u) = ⟨u⊗u|N|u⊗u⟩` and then that "a Hermitian form on `Sym²` is determined by its values
on the Veronese `{u⊗u}`", justifying it with: *"both sides have `28²` real dimensions and the
map is injective."* That is not a justification. Equality of dimensions does **not** imply
injectivity, and I asserted the injectivity I was supposed to prove. The claim is true, and
the proof is one line, which I have now written and checked here (S1.1): the `784` functions
`conj(m_α)m_β` are pairwise distinct monomials in `(u, ū)` — computed, `0` duplicates among
`784` ordered pairs of the `28` degree-2 monomials — hence linearly independent, hence the
map is injective. That sub-step is load-bearing: without it the Schur conclusion that made
my whole stage-1 reduction work does not follow. Reading this note forced me to look at it,
because the author's Step 1 leans on exactly the same bridge and states it no better.

**Two further items, smaller, and to my stage-1 credit I did flag both.** My stage-1 return
named its two weakest links as (i) the off-block vanishing of `N`, checked numerically at
`8.3e−17` rather than exactly, and (ii) the Majorana/Hopf ingredients of my three-axes
formula, verified at seven configurations rather than proved. This note's Step 1 shows that
(i) was avoidable: four weight-state evaluations plus one `4×4` determinant pin the whole
invariant identity in exact rationals with no off-block computation anywhere. And the
author's Step 3 shows (ii) was avoidable too: the eigenvalue lemma reaches the same bound
with no Majorana representation at all. Neither is an error in stage 1, but both are places
where I took an exposed route when a cleaner one existed.

**One disagreement running the other way, about the texts rather than the mathematics.** In
stage 2a I graded part 3c of an excerpt whose equality parenthetical named **two** loci, and
I reported, with the computation, that those two deliver only two of the three permutations.
The full note here names **three** loci, adding "`a = −√6/4` in `M₃`", and I have computed
that the three do deliver all three points with nothing missed. So the specific finding that
carried my stage-2a verdict does **not** apply to this text. My stage-2a verdict stands for
the text I was given; against this text the equality set is in better shape, and what
remains missing in 3c is narrower: `M₂`'s locus untraced, the sign conditions underived, and
the `a+6b ≤ 0` branch not analysed for equality. Both texts land on ESTABLISHED, SUPPLIED,
but for different reasons, and I would not want the stage-2a reasoning read as applying here.

---

## 7. Where I considered DEFECT and did not take it

| Candidate | Why not |
| --- | --- |
| "`M₂` = `M₁` at `−b`" (line 40) | False as a matrix identity for the block in the natural basis, but true after a sign flip of one basis vector, and identical characteristic polynomials either way (verified). The author never fixes a basis, so the statement is realisable. Downstream only spectra are used. |
| "The right side of `15/√6 − a − 6b ≥ √(…)` is positive" (line 45) | Left side called the right side. Every clause is individually true, and the justification given is the one the left side needs. A wording slip, not a false statement. |
| `TrN̄ = 12` never stated (line 30) | It is a theorem about the objects (`Σf_i² = 12` at spin 3), not a hypothesis that could fail, so it is an omitted derivation rather than an unstated hypothesis. Supplied as S3a.1. Flagged prominently because the supplied checker has no gate for it. |

---

## 8. Method note, and a correction to my own working

| Item | Detail |
| --- | --- |
| Provenance | Every object rebuilt from this room's `worklist.md` conventions. Nothing imported from the author's files; `check_s3_maximum.py` read as text, never executed, and I did not look for its missing modules. |
| Exact vs numeric | Exact (sympy over `Q`, `Q[√6]`): the weight-state table, the determinant `180/7`, the coefficients, `TrN̄ = 12`, `‖Θu‖ = ‖u‖`, the operator identity and the block split, the three `M₁`/`M₃` identities, the per-block equality sets and their factorisations, the `3e` spectrum, the Step-4 formulas, `r6(s) = 463/924 − cos²(2s)/2`, the `n_J` solve, and `TrN̄² = 48‖u‖⁴ + 126p₂`. Double precision, quoted to 15 digits: the variational identity (`2.08e−15`), the equivariance (`1.91e−14`), the `λ_max` formulas (`1.421e−14` over 20000 draws), `3d` on the hexagon, `A_E` for random `n` (`2.22e−15`), the 200 rotated hexagons (`2.11e−15`). 30-digit: the Step-1 identity residual (`9.643e−18`). |
| **A `False` in my own run, reported not smoothed** | `s2b_02_steps23.py` printed `False` for "`M₁` top eigenvalue formula is a root of `char(M₁)`" and for the general `M₃` formula. Both are **my** script's failure, not the author's: `charpoly(...).subs(surd)` followed by `simplify` did not close. `s2b_03_eigformulas.py` settles both three ways — `tr/2` and the discriminant match exactly, `expand(λ² − Tλ + D) = 0` identically, and `max |λ_max − formula| = 1.421e−14` over 20000 random `(a,b)` for each block. The note's formulas are correct. |
| Containment | Read: `stage2b_grading.md`, `S0_S3_MAXIMUM.md`, `check_s3_maximum.py`, `check_s3_maximum_log.txt`, and my own files. Nothing outside `/tmp/wl-audit/s1`; no network. `return.md` and `stage2a/return2a.md` were not modified. |

### Scripts written for this stage

| Script | Role |
| --- | --- |
| `s2b_01_step1.py` | Step 1: CG closed form vs sympy, weight-state table, `180/7`, the coefficients, the identity at 30 digits, `TrN̄ = 12`, the injectivity bridge |
| `s2b_02_steps23.py` | Step 2; 3a (variational identity, equivariance, max-swap, substitution, block split); 3b identities and side conditions; 3c complete per-block equality sets |
| `s2b_03_eigformulas.py` | runs down the two `False` results in `s2b_02`; confirms both `λ_max` formulas three ways |
| `s2b_04_step4_reconcile.py` | 3d, 3e, Step 4, and the reconcile (independent `n_J` solve, `TrN̄² = 48 + 126p₂` on the real form, hexagon orbit) |
| `s2b_05_span_formula.py` | `r6` on `span{\|3,3⟩,\|3,−3⟩}` exactly, and the strictness Step 4 needs |
| `s2b_06_checker_coverage.py` | gate counts, log attribution, sweep spacing and the excluded band, the inert line 61, the unsupplied imports |
