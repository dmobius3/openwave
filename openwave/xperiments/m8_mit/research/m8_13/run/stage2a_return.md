# Return, stage 2a: grading part 3c of `step3_text.md`

## Verdict

> **3c: ESTABLISHED, SUPPLIED.**

The text's equality analysis for part 3c is the single sentence on line 19:

> "Tracing the equality cases (`b = 0` with `a > 0` and `a = −2b` with `b > 0`) shows
> equality only at the three permutations of `(2, −1, −1)/√6`."

Its conclusion is **correct** (I confirm the equality set is exactly those three
permutations). Its stated route does **not** reach that conclusion. Two of the three
blocks never have their equality points traced at all, and the two loci the sentence
names deliver two points, not three. Everything missing is a case or a derivation that
I can supply and do supply below, so this is not a GAP; nothing the text asserts is
mathematically false, so I did not grade it DEFECT (§ 6 argues that call explicitly).

---

## 1. Are every block's equality points found? No.

Answer, block by block. `L₀ := 15/√6`; I use the rationalising substitution
`a = √6·u`, `b = √6·v`, under which the constraint ellipse `(2/3)a² + 8b² = 1`
becomes `4u² + 48v² = 1`.

| Block | Its equality points on the ellipse (computed, § 2) | Traced in the text? |
| --- | --- | --- |
| `M₁` | `(u,v) = (1/2, 0)` and `(−1/4, 1/8)` | **yes** — these are exactly the two named loci |
| `M₂` | `(u,v) = (1/2, 0)` and `(−1/4, −1/8)` | **no** — the `M₂` bullet establishes only the *bound* ("it has `M₁`'s spectrum at `−b`, and the constraint is even in `b`"); its equality locus is never written down |
| `M₃` | `(u,v) = (−1/4, 1/8)` and `(−1/4, −1/8)` | **no** — the `M₃` bullet stops at `20(a + √6/4)² ≥ 0` and never remarks that this vanishes, i.e. never notices that `M₃` *attains* `L₀` |
| seventh eigenvalue (`≡ 0`) | none | **no** — never mentioned; trivially fine since `L₀ = 6.1237243569579452455 > 0` |

So: **one block of three has its equality points found.**

Two concrete consequences, both computed:

**(a) The named loci give two points, not three.** Intersecting each named locus with
the ellipse:

| Named locus | Point on the ellipse | `e` |
| --- | --- | --- |
| `b = 0` with `a > 0` | `(u,v) = (1/2, 0)`, i.e. `(a,b) = (√6/2, 0)` | `(−1, −1, 2)/√6` |
| `a = −2b` with `b > 0` | `(u,v) = (−1/4, 1/8)`, i.e. `(a,b) = (−√6/4, √6/8)` | `(2, −1, −1)/√6` |

That is **2** of the 3 permutations. The third, `e = (−1, 2, −1)/√6`, sits at
`(u,v) = (−1/4, −1/8)`, i.e. `(a,b) = (−√6/4, −√6/8)`, on the locus `a = 2b` with
`b < 0`, which is **not** among the named loci. Verified: at that point
`λ_max = 6.123724356957945 = L₀`.

**(b) The omission of `M₃` is load-bearing for the text's own closing sentence.**
Line 21 asserts the equality eigenspace is `span{|3,3⟩ₙ, |3,−3⟩ₙ}`, a **2-dimensional**
space, so `L₀` must occur with multiplicity 2 at every equality point. Computed spectra
of the full `7×7` operator `A = a(f_z²−4) + b(f₊²+f₋²)`:

| Point `(u,v)` | `e` | spectrum | multiplicity of `L₀` | blocks supplying the two copies |
| --- | --- | --- | --- | --- |
| `(1/2, 0)` | `(−1,−1,2)/√6` | `−4.898979486, −3.674234614, −3.674234614, 0, 0, 6.123724357, 6.123724357` | 2 | `M₁`, `M₂` |
| `(−1/4, 1/8)` | `(2,−1,−1)/√6` | same multiset | 2 | `M₁`, `M₃` |
| `(−1/4, −1/8)` | `(−1,2,−1)/√6` | same multiset | 2 | `M₂`, `M₃` |

At `(−1/4, −1/8)` **both** copies of `L₀` come from `M₂` and `M₃`, the two blocks whose
equality loci the text never traces. So the gap is not cosmetic: a reader following the
text's tracing literally recovers neither that point nor the 2-dimensionality that
line 21 depends on.

---

## 2. How I know my list is complete

This is the part the grading note flags, so I state the mechanism rather than the list.

**The decision procedure.** For a real symmetric `2×2` matrix `M`, the statement
`λ_max(M) = L₀` is **exactly equivalent** to

```text
   det(M - L0 I) = 0        (L0 is an eigenvalue)
   AND   tr M <= 2 L0       (it is the larger of the two)
```

This is an equivalence, not an implication, so no equality point can escape it. It
converts "find every point where this block reaches `L₀`" into "intersect one conic with
one ellipse", a **zero-dimensional** system.

**The elimination is explicit, not a search.** `det(M − L₀I) = 0` is a conic in `(u,v)`
containing `v` only through `v²` and `uv`; reducing it modulo the ellipse
(`48v² = 1 − 4u²`) removes `v²` and leaves a polynomial that **factors over `Q`**:

| Block | `det(M − L₀I)` in `(u,v)` | reduced modulo the ellipse |
| --- | --- | --- |
| `M₁` | `−90u² + 360uv − 30u − 360v² − 180v + 75/2` | `−30·(2u − 1)·(u − 6v + 1)` |
| `M₂` | `−90u² − 360uv − 30u − 360v² + 180v + 75/2` | `−30·(2u − 1)·(u + 6v + 1)` |
| `M₃` | `60u − 1440v² + 75/2` | `(15/2)·(4u + 1)²` |

Each factor is linear, so intersecting it with the ellipse is a single quadratic whose
**complete** root set sympy returns. `factor_list` over `Q` is a complete factorisation,
and a linear factor meets the ellipse in at most 2 points; there is nowhere for a root to
hide. Working out:

- `2u − 1 = 0` ⟹ `u = 1/2` ⟹ `48v² = 1 − 1 = 0` ⟹ `v = 0`. One point.
- `u − 6v + 1 = 0` ⟹ `4(6v−1)² + 48v² = 1` ⟹ `192v² − 48v + 3 = 0` ⟹ `(8v − 1)² = 0`
  ⟹ `v = 1/8`, `u = −1/4`. One point (a double root).
- `u + 6v + 1 = 0` ⟹ `(8v + 1)² = 0` ⟹ `v = −1/8`, `u = −1/4`. One point.
- `(4u + 1)² = 0` ⟹ `u = −1/4` ⟹ `48v² = 3/4` ⟹ `v = ±1/8`. Two points.

Then the trace test: at all of these `tr M = √6` and `2L₀ = 5√6`, so `tr M ≤ 2L₀` holds
and `L₀` really is the **larger** root in every case; confirmed numerically,
`λ_max = 6.123724356957945` at each.

**The parametrisation loses nothing.** The map `(a,b) ↦ e = (−a/3 + 2b, −a/3 − 2b, 2a/3)`
is linear with the explicit inverse `a = 3e₃/2`, `b = (e₁ − e₂)/4`, so it is a
**bijection** from the ellipse onto the unit traceless diagonal `e`. Enumerating ellipse
points therefore enumerates all such `e`; no `e` is unreachable. (Verified symbolically.)

**Independent corroboration, stated with its precision.** A dense sweep of
**4,000,001** equally spaced points of the ellipse (`a = √(3/2)cos t`,
`b = sin t/(2√2)`), diagonalising the full `7×7` operator at each, gives

- `max λ_max = 6.123724356957945` against `L₀ = 6.123724356957946`, overshoot
  `−8.882e−16` (double precision; the bound is never exceeded);
- exactly **3** clusters of points within `1e−9` of `L₀`, centred at
  `(u,v) = (0.500000000, 0.000000000)`, `(−0.250000453, 0.124999924)`,
  `(−0.249999773, −0.125000038)` — the three exact points, to sweep resolution.

The sweep is corroboration only: at `1e−9` tolerance and `4e6` samples it cannot
establish that no fourth equality point exists. The completeness claim rests on the
exact equivalence plus the complete factorisation above; the sweep is there to catch an
algebra slip, and it did not fire.

**The discarded signs, checked rather than assumed.** The text restricts to `a > 0` and
`b > 0` without derivation. Those restrictions are in fact **automatic** from the
factorisations (`2u − 1 = 0` forces `v = 0` and `u = 1/2 > 0`; `u − 6v + 1 = 0` forces
`v = 1/8 > 0`), and the discarded sign choices really are strictly below the bound:

| Discarded point | `λ_max(M₁)` | `λ_max(M₂)` | `λ_max(M₃)` | `L₀` |
| --- | --- | --- | --- | --- |
| `b = 0, a < 0`, i.e. `(u,v) = (−1/2, 0)` | `3.674234614174767` | `3.674234614174767` | `4.898979485566356` | `6.123724356957945` |
| `a = −2b, b < 0`, i.e. `(u,v) = (1/4, −1/8)` | `3.674234614174767` | `4.898979485566356` | `3.674234614174767` | `6.123724356957945` |

---

## 3. The equality set, as I find it

| `(u,v)` | `(a,b)` | `e` | attained by |
| --- | --- | --- | --- |
| `(1/2, 0)` | `(√6/2, 0)` | `(−1, −1, 2)/√6` | `M₁`, `M₂` |
| `(−1/4, 1/8)` | `(−√6/4, √6/8)` | `(2, −1, −1)/√6` | `M₁`, `M₃` |
| `(−1/4, −1/8)` | `(−√6/4, −√6/8)` | `(−1, 2, −1)/√6` | `M₂`, `M₃` |

Three distinct points, the three permutations of `(2, −1, −1)/√6`. The text's stated
answer is right; its stated route reaches only the first two rows.

---

## 4. Exactly what I supplied

Named precisely, and shown above rather than asserted:

| # | Supplied | Where |
| --- | --- | --- |
| S1 | **The decision procedure** `λ_max(M) = L₀ ⟺ det(M − L₀I) = 0 ∧ tr M ≤ 2L₀`, which is what turns "every block's equality points" into a finite, complete computation rather than an inspection. The text has no analogue: it reasons through two nested squarings for `M₁` and through a one-sided inequality for `M₃`. | § 2 |
| S2 | **`M₂`'s equality locus.** The text's `M₂` bullet gives the ingredient (`spec M₂(a,b) = spec M₁(a,−b)`, ellipse even in `b`) but never draws the equality consequence. Supplied: `eq(M₂) = ` reflection of `eq(M₁) = {(1/2,0), (−1/4,−1/8)}`, and independently by the factorisation `−30(2u−1)(u+6v+1)`. This is what produces the **third** permutation. | § 1(a), § 2 |
| S3 | **`M₃`'s equality locus, absent from the text entirely.** Supplied: reduced modulo the ellipse, `det(M₃ − L₀I) = (15/2)(4u+1)²`, a **double** root at `u = −1/4` (so `M₃` is tangent to the bound), giving `(u,v) = (−1/4, ±1/8)`; and the verification that both already lie in the union, which is what rescues the word "only" in the text's claim. | § 2 |
| S4 | **The zero-eigenvalue case**, never mentioned: `0 ≠ L₀` because `L₀ = 6.1237243569579452455 > 0`. | § 1 |
| S5 | **Derivation of the sign conditions** `a > 0`, `b > 0`, which the text asserts bare. Supplied: they are forced by the factorisations, and the discarded signs give `λ_max ∈ {3.674234614174767, 4.898979485566356}`, both `< L₀`. | § 2 |
| S6 | **The bijection** `(a,b) ↔ e`, so that enumerating ellipse points enumerates all unit traceless diagonal `e`. | § 2 |
| S7 | **The multiplicity count** at each equality point (`L₀` doubly degenerate everywhere), which links 3c to the text's own line 21 and shows the `M₃` omission is load-bearing. | § 1(b) |
| S8 | **Re-derivation of the whole setup** the text states without proof, so that none of the above rests on the text being right about its own objects: `f² = 12I`; `e₁f_x² + e₂f_y² + e₃f_z² = a(f_z²−4) + b(f₊²+f₋²)` under `e₃ = 2a/3`, `e₁−e₂ = 4b`; `‖e‖² = (2/3)a² + 8b²`; the parity/flip basis is orthonormal and block-diagonalises `A` into `M₁ ⊕ M₂ ⊕ M₃ ⊕ (0)`; `M₁ = [[5a, 2√15 b],[2√15 b, −3a+12b]]` with `2√15 = √60`, `M₃ = [[0, 4√15 b],[4√15 b, −4a]]` with `4√15 = √240`; seventh diagonal entry `= 0`. All returned `True`. | `s2a_01_blocks.py` |

---

## 5. What I did **not** need to supply, and what in the text is right

To be fair to the text, these check out as written and I verified each:

| Text's claim | Status |
| --- | --- |
| `M₁`'s two equality loci as named | correct: they are precisely `eq(M₁)` |
| `M₃` bound via `(15/√6 + 2a)² − (30 − 16a²) = 20(a + √6/4)²` and `15/√6 + 2a ≥ 9/√6 > 0` | correct as a bound; `|a| ≤ √6/2` on the ellipse so `15/√6 + 2a ≥ 15/√6 − √6 = 9/√6` |
| `a + 6b ≤ √6` on the ellipse, so the first squaring is valid | correct: the Cauchy–Schwarz maximum of `a + 6b` on the ellipse is `√(3/2 + 36/8) = √6`, and `L₀ ≈ 6.124 > √6 ≈ 2.449` |
| `M₁`'s final difference `36b²(a + 2b)²` | consistent with my factor `u − 6v + 1` meeting the ellipse only at `a = −2b`; both name the same single point |
| The final answer: three permutations of `(2,−1,−1)/√6` | correct |

One wording slip, noted but not graded: the text writes "The right side of
`15/√6 − a − 6b ≥ √(…)` is positive", while `15/√6 − a − 6b` is the **left** side of the
displayed inequality. The intended quantity is the one that is positive, so this is a
typo, not an error of substance.

---

## 6. Why this verdict and not another

| Verdict | Why not |
| --- | --- |
| **ESTABLISHED** | Ruled out. Two of the three blocks' equality loci are never traced, and the two loci that are named yield 2 points against an asserted 3. The grading rule is explicit that a derivation I write out is a supplied part even when the text holds its ingredients, so "as written" is unavailable here in any case. |
| **GAP** | Ruled out. Every missing piece is completable, and I completed it in § 2: the block equality loci come out of a complete factorisation of a zero-dimensional system. Nothing here is beyond reach. |
| **DEFECT** | The closest competitor, and I record the case for it rather than burying it. On a literal reading, "Tracing the equality cases (`b = 0` with `a > 0` and `a = −2b` with `b > 0`) shows equality only at the three permutations" is false as a claim about what that tracing shows: I computed that those two loci meet the ellipse in exactly **2** points, giving 2 of the 3 permutations. I did not take DEFECT because (i) every individual mathematical statement in the text is true — the named loci *are* equality loci, the bounds *are* valid, the final set *is* right; (ii) no necessary hypothesis is hidden — the `M₂` reflection principle the reader needs is already stated in the `M₂` bullet; (iii) what fails is that the case analysis is incomplete, and "omits a case or a derivation, and you supply it" is the rubric's own description of ESTABLISHED, SUPPLIED. A grader who weights the count mismatch as a false statement rather than as an elliptical one would land on DEFECT, and the computed facts in § 1(a) are what that decision should turn on. |

---

## 7. Method note

| Item | Detail |
| --- | --- |
| Source of the spin-3 operators | Built from the conventions in this room's `worklist.md` § 1.2 (`J_z v_m = m v_m`, `J_± v_m = √(12 − m(m±1)) v_{m±1}`), not taken from `step3_text.md`. `f² = 12I` verified. |
| Every number asserted | Produced by the four scripts below and quoted from their output; none is arithmetic done in my head and none is read off `step3_text.md`. |
| Exact vs numeric | The block equality sets, the factorisations, the setup identities, the `(a,b) ↔ e` bijection and the trace tests are **exact** (sympy over `Q` and `Q[√6]`). The `λ_max` values, the spectra, the multiplicities and the ellipse sweep are **double precision**, quoted to 15 digits; the sweep's tolerance is `1e−9` over `4,000,001` samples and it corroborates only. |
| Kawaguchi–Ueda | The text attributes the bound to a review that "states this bound without proof". I did not and could not consult it; I graded the text's own reasoning only. |
| Containment | I read `stage2a_grading.md` and `step3_text.md` and my own files. Nothing outside `/tmp/wl-audit/s1` was read, listed or searched; no network. My committed stage-1 return was not edited. |

### Scripts written for this stage

| Script | Role |
| --- | --- |
| `s2a_01_blocks.py` | rebuilds the spin-3 operators, the substitution, the constraint and the block decomposition from scratch, exactly; confirms `M₁`, `M₂`, `M₃` and the seventh (zero) entry |
| `s2a_02_equality_set.py` | first attempt at the equality set via `sympy.solve` on the radical system; **killed**, ran past budget, kept unedited as the record |
| `s2a_02b_equality_set.py` | second attempt with surds cleared; also **killed**, same reason, kept unedited |
| `s2a_03_equality_fast.py` | the working route: explicit elimination modulo the ellipse, complete factorisation, per-block equality points, the two named loci, the discarded signs, and the 4,000,001-point sweep |
| `s2a_04_multiplicity.py` | the `(a,b) ↔ e` bijection and the multiplicity of `L₀` at each equality point |
