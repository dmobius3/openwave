# The maintainer's own grading of the argument under audit

Written and frozen **before** the stage-1 room is given the author's argument,
so that the adjudication has a second grading that the room's cannot anchor.
Same scale as the room's: ESTABLISHED · ESTABLISHED, SUPPLIED · GAP · DEFECT.

Every number asserted below comes from `check_maint.py` (69 passed, 0 failed,
2 deliberately-dead arms) or `poly_identity.py`, both written here.

Calibration used: **ESTABLISHED** where the text's own sentences carry the
inference and a reader following them reaches the conclusion without adding an
idea. **ESTABLISHED, SUPPLIED** where I had to write a computation or a case
the text does not contain, even when the text contains its ingredients.

## Step 1, the invariant form: ESTABLISHED, SUPPLIED

The route is: the invariant quartics are four-dimensional because `Sym²V₃` is
multiplicity-free; the four listed functions are invariant quartics; their
values at `v₃, v₂, v₁, v₀` have determinant `180/7 ≠ 0`, so they are a basis;
`‖ρ₆‖²`'s weight-state values then force the coefficients.

**Supplied.**

1. The dimension count itself. Multiplicity-freeness gives a four-dimensional
   commutant by Schur, but the passage to *quartics* needs the `U(1)` factor:
   invariance under `u ↦ e^{iφ}u` restricts to type `(2,2)`, since type `(4,0)`
   carries `e^{4iφ}` and type `(3,1)` carries `e^{2iφ}`. Only then is the space
   in question the commutant of `Sym²V₃`. The note states the conclusion and
   the reason but not this restriction.
2. That each of the four functions is invariant. For `|a₀₀|²` and `‖ρ₆‖²` this
   needs `Θ` to intertwine the rotation action, which is true and unstated.

**Not needed after all.** I verified the identity

`‖ρ₆‖²·‖u‖⁴ = −(5/231)‖u‖⁴ − |f|²/22 + (7/11)|a₀₀|² + TrN̄²/198`

as an **exact polynomial identity in the 14 real coordinates** (`poly_identity.py`:
residual `0`, and the arm `1/198 → 1/200` leaves 257 terms). So the conclusion
stands on its own, with no dimension count and no fit. That is a stronger
verification than either the note or its checker carries: the checker tests the
identity on 60 random states to `1.7e-16`.

Recomputed independently: determinant `180/7`, coefficients
`(−5/231, −1/22, 7/11, 1/198)`, weight-state values `1, 36, 225, 400` over `924`.

## Step 2, the two bounds: ESTABLISHED

`|f|² ≥ 0` with equality iff `⟨f⟩ = 0` is a tautology. Cauchy-Schwarz gives
`|a₀₀|² = |⟨Θu, u⟩|²/7 ≤ ‖Θu‖²‖u‖²/7 = 1/7`, with equality iff `Θu ∝ u`.

The one unstated ingredient, `‖Θu‖ = ‖u‖`, is immediate from the definition
(`Θ` permutes `m ↦ −m` and multiplies by a unit phase), so I do not count it as
supplied. `Θ² = 1` at integer spin, checked here.

The phrase "time-reversal invariant" is loose: equality gives `Θu = λu` with
`|λ| = 1`, and a global phase then makes `Θu = u`. Not a defect; `u` is defined
up to phase throughout.

## Step 3a, reduction and block split: ESTABLISHED, SUPPLIED

**Supplied.**

1. The duality `‖Q‖ = max_E Tr(QE)` over unit traceless symmetric `E`, and
   `Tr(N̄E) = ⟨u, A_E u⟩` for traceless `E`. Checked here (7.4). The note
   asserts the first and does not state the second at all, although the second
   is what turns a statement about `N̄` into one about an operator on `V₃`.
2. Rotation equivariance `A_{RER^T} = D(R) A_E D(R)†`, which is what makes
   "rotating `E` to diagonal form does not change `λ_max`" true.
3. The block basis: which pairs of `m` sit in which block. The note names the
   two symmetries and gives the matrices, but not the pairing.

Recomputed independently and exactly: `e₁fₓ² + e₂f_y² + e₃f_z² = a(f_z²−4) +
b(f₊²+f₋²)` as 7×7 matrices under `e₃ = 2a/3`, `e₁−e₂ = 4b`; the unit-traceless
condition is exactly `(2/3)a² + 8b² = 1`; and the characteristic polynomial
factors as `λ·det(M₁−λ)·det(M₂−λ)·det(M₃−λ)`. Two arms fire on this: a
`(1/50)fₓ²` perturbation of the operator, and `√60 → √61` in `M₁`.

## Step 3b, each block's bound

| Block | Grade | Note |
| --- | --- | --- |
| `M₃` | ESTABLISHED | `λ_max = −2a + √(30−16a²)` uses `240b² = 30 − 20a²`, which is the ellipse, and the note has already put the reader on the ellipse. The completion `(15/√6+2a)² − (30−16a²) = 20(a+√6/4)²` is exact (checked), and `15/√6 + 2a ≥ 9/√6 > 0` follows from `\|a\| ≤ √6/2` |
| `M₁` | ESTABLISHED, SUPPLIED | The homogenization and the sum of squares `(a²+6ab+24b²)² − (3/2)(a+6b)²N² = 36b²(a+2b)²` are exact (checked), and the arm `24b² → 23b²` fires. **Supplied:** `max(a+6b) = √6` on the ellipse, which the note asserts. One line of Cauchy-Schwarz: `a+6b = √(3/2)cos t + (6/(2√2))sin t` has amplitude `√(3/2 + 36/8) = √6`. It is load-bearing, because it is what makes the first squaring valid |
| `M₂` | ESTABLISHED | `M₂` is `M₁` at `−b` and the ellipse is even in `b`; both checked exactly |

## Step 3c, the equality set: ESTABLISHED, SUPPLIED

This is the substantive part, and it is where the note is thinnest. It reads:
"Tracing the equality cases (`b = 0` with `a > 0`, `a = −2b` with `b > 0`, and
`a = −√6/4` in `M₃`) shows equality only at the three permutations of
`(2, −1, −1)/√6`."

**Supplied, in four parts.**

1. **Where the sign conditions come from.** On the branch `a + 6b > 0`, equality
   in the twice-squared inequality is `36b²(a+2b)² = 0`, so `b = 0` or
   `a = −2b`. That is the derivation of `M₁`'s two case lines; the note states
   them as conditions rather than deriving them.
2. **Which point of each line survives.** `b = 0` meets the ellipse at
   `a = ±√6/2`, and `a + 6b > 0` keeps `a = +√6/2`. `a = −2b` meets it at
   `b = ±√6/8`, and `a + 6b = 4b > 0` keeps `b = +√6/8`. Each line meets the
   ellipse twice and only one of each pair is an equality point, which is
   exactly what the sign conditions are doing and what the note does not say.
3. **The other branch.** `a + 6b ≤ 0` is not addressed at all. There the bound
   reads `(3/√6)(a+6b)N ≤ (a+3b)² + 15b²` with the left side `≤ 0` and the
   right `≥ 0`, so equality needs both to vanish, and `(a+3b)² + 15b² = 0` has
   **no real solution on the ellipse** (checked: 0 solutions). Without this the
   tracing is not exhaustive, whatever the sign conditions say.
4. **`M₂`'s points**, by `b ↦ −b`, and `M₃`'s, from `a = −√6/4` with the ellipse
   giving `b = ±√6/8`.

**Recomputed independently, by a different route from the author's.** Rather
than tracing cases, I solved, for each block separately, the characteristic
equation `det(M − (15/√6)I) = 0` together with the ellipse, and kept a solution
only when `15/√6` is that block's **larger** eigenvalue there. This returns

| point `(a, b)` | reached by | `√6·e` |
| --- | --- | --- |
| `(√6/2, 0)` | `M₁`, `M₂` | `(−1, −1, 2)` |
| `(−√6/4, √6/8)` | `M₁`, `M₃` | `(2, −1, −1)` |
| `(−√6/4, −√6/8)` | `M₂`, `M₃` | `(−1, 2, −1)` |

and nothing else: **exactly three points, the three permutations of
`(2, −1, −1)/√6`, each reached by exactly two blocks.** Three arms fire on this
solve (`15/√6 → 14/√6`, ellipse constant `1 → 51/50`, `√60 → √61`).

⭐ This is the one place where the note's own checker does not reach the claim.
Its `(L2)` verifies the inequalities exactly, but its equality evidence is a
20001-point sweep plus "strictly below away from them: max over points more
than 0.01 rad from an axis". A sweep rejects candidates; it cannot show there
are no others, and the `±0.01 rad` neighborhoods are left to the grid. The
exact per-block solve above is what closes it.

**A trap I fell into and record rather than smooth.** My first version of the
branch check asserted that *no* equality point has `a + 6b ≤ 0`. That is false:
the point reached by `M₂` and `M₃` has `a + 6b = −√6`. `a + 6b > 0` is **`M₁`'s**
squaring branch, not a property of the equality set. `M₂`'s mirror branch is
`a − 6b > 0`, and `M₃` has no such branch at all. A grader who conflates them
will either reject the true set or accept a false argument for it.

## Step 3d, from `‖Q‖ = 15/√6` to the top eigenspace: ESTABLISHED, SUPPLIED

**Supplied:** the chain itself. Take `E* = Q/‖Q‖`, which is unit, symmetric and
traceless. Then `⟨u, A_{E*}u⟩ = Tr(QE*) = ‖Q‖ = 15/√6`, while
`⟨u, A_{E*}u⟩ ≤ λ_max(A_{E*}) ≤ 15/√6` by step 3b. Both are therefore
equalities: `E*` is itself an extremal `e`, and `u` lies in the top eigenspace
of `A_{E*}`. The note states the conclusion and leaves this out. That
`⟨u,Au⟩ = λ_max(A)` forces `u` into the top eigenspace is standard.

`A_E` for `E = (3nnᵀ − I)/√6` is `(3(n·f)² − 12)/√6`, using `Σf_i² = 12·I`
(checked here, 7.1).

## Step 3e, the eigenspace is `span{|3,3⟩ₙ, |3,−3⟩ₙ}`: ESTABLISHED, SUPPLIED

Asserted in the note. **Supplied:** the computation. At `√6·e = (−1,−1,2)`,
that is `n = ẑ`, the top eigenvalue of `a(f_z²−4) + b(f₊²+f₋²)` is `15/√6` with
**multiplicity two**, and the eigenspace is exactly `span{v₃, v₋₃}` (checked:
every top eigenvector has zero component on `v₂ … v₋₂`). The general `n`
follows by rotating. The multiplicity is two at all three points, which is
diagnostic **D1** and is consistent with each point being reached by two
blocks: one vector from each.

The author's checker verifies the converse direction (every `u` on the span has
`TrN̄² = 171/2`), not this one.

## Step 4, the combination: ESTABLISHED, SUPPLIED

The sign structure is right: the identity carries `−|f|²/22`, so the maximum
wants `|f|²` smallest and the other two largest, and all three are attainable
together. `−5/231 + 1/11 + 171/396 = 463/924`, checked.

**Supplied:** that `|a₀₀|² = 1/7` "follows". On the span,
`|a₀₀|² = 4|α|²|β|²/(7(|α|²+|β|²)²)`, so `|α| = |β|` gives `1/7`. Checked here,
along with `|f|² = 9(|α|²−|β|²)²` and `TrN̄² = 171/2` identically on the span.

Worth recording, and absent from the note: on that span **either** condition
closes it alone. `|a₀₀|² = 1/7` is `−(4/7)(|α|²−1/2)² = 0`, which forces
`|α|² = 1/2` by itself, without using `⟨f⟩ = 0`.

"The remaining relative phase is a rotation about `n`" is right: `R_z(χ)` sends
the relative phase to itself plus `6χ`, so every relative phase is reached, and
the equality set is one orbit. Checked. `r̂₆` on the span is `−2t² + 2t + 1/924`
in `t = |α|²`, with its only critical point at `t = 1/2` and value `463/924`.

## Summary

| Step | Grade | What I supplied |
| --- | --- | --- |
| 1 | ESTABLISHED, SUPPLIED | the `U(1)` type restriction behind the dimension count; invariance of each function. Conclusion independently proved as a polynomial identity |
| 2 | ESTABLISHED | nothing |
| 3a | ESTABLISHED, SUPPLIED | the duality and `Tr(N̄E) = ⟨u, A_E u⟩`; rotation equivariance; the block basis |
| 3b `M₃` | ESTABLISHED | nothing |
| 3b `M₁` | ESTABLISHED, SUPPLIED | `max(a+6b) = √6` on the ellipse, which validates the first squaring |
| 3b `M₂` | ESTABLISHED | nothing |
| 3c | ESTABLISHED, SUPPLIED | the derivation of the sign conditions; which point of each line survives; **the branch `a+6b ≤ 0`**; `M₂`'s and `M₃`'s points |
| 3d | ESTABLISHED, SUPPLIED | the `E* = Q/‖Q‖` chain |
| 3e | ESTABLISHED, SUPPLIED | the eigenspace computation and its multiplicity two |
| 4 | ESTABLISHED, SUPPLIED | `\|a₀₀\|² = 4\|α\|²\|β\|²/7` on the span |

**No GAP and no DEFECT.** Every step is completable from what the note gives,
and the conclusion, the maximum `463/924` attained exactly on the hexagon
orbit, is independently established here.

The heaviest supplied part is 3c's branch analysis, which is what the
pre-registration predicted, and it is load-bearing: without it the tracing is
not exhaustive.

## My grade for the N2 variant's 3c

The variant differs from the note in exactly one line (diffed: one changed
line, the tracing sentence), which drops `and a = −√6/4 in M₃`.

**ESTABLISHED, SUPPLIED**, with `M₃`'s equality case supplied, plus everything
supplied for the real 3c above.

The conclusion survives the omission, and I can say why from my own solve
rather than from the note: `M₃`'s two equality points are `(−√6/4, ±√6/8)`,
which are already `M₁`'s and `M₂`'s, so the **set** loses nothing. What is lost
is the argument that `M₃` contributes nothing further: with the case dropped,
the text bounds `M₃` and then never revisits it, so a reader has no reason to
believe `M₃` does not reach `15/√6` somewhere the other two do not. The
tracing sentence claims exhaustiveness over all three blocks while having
examined two.

It is not a GAP, because the case takes two lines: `20(a+√6/4)² = 0` pins
`a = −√6/4`, and the ellipse then gives `b² = 3/32`.

It is **not** ESTABLISHED as written, for the same reason the real 3c is not.

⚠️ The asymmetry that makes this control sharp: dropping `M₃`'s case leaves the
set intact, while dropping `b = 0` with `a > 0` would lose a point outright,
because that case carries `(√6/2, 0)` for `M₁` and, through the symmetry
sentence, for `M₂`. A grader reading conclusions rather than arguments cannot
tell the two omissions apart.
