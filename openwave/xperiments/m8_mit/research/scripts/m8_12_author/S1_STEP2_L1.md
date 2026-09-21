# S1 derivation, step 2 of 3: L1, the critical set on the seven lines

2026-09-19. Author-side and local, for the units' review before step 3. It computes no transverse Hessian; the in-line second derivative appears only as the Poincaré–Hopf index and in the exposure below.

## Result

Together with L0, this is the partial classification the outline promised: **every critical orbit of `r̂₆` whose projective stabilizer has a fixed locus of projective dimension at most one is one of ten orbits.** They are exactly the ten of the S0 memo's § 3, so the lines carry no further orbit.

| `924·r̂₆` | orbit | rotation orders | where it sits |
|---|---|---|---|
| 1 | coherent `v₃` | one-parameter | endpoints of the C₄ `{v₃, v₋₁}`, C₅ and C₆ lines |
| 36 | `v₂` | one-parameter | endpoints of the C₃, C₄ `{v₂, v₋₂}` and C₅ lines |
| `8800/43` | prism | 2, 3 | D₃ line at `z = ±√230/10`, in-line minima |
| 225 | `v₁` | one-parameter | endpoints of the C₃ and C₄ `{v₃, v₋₁}` lines |
| 225 | D2 ray | 2 | D₂ line at `z = ±i√30/5`, in-line minima |
| `1188/5` | C3 ray | 3 | C₃ line, interior circle at `s* = 1/5` |
| `1188/5` | pyramid | 5 | C₅ line, interior circle at `s* = 12/25` |
| 288 | octahedron | 2, 3, 4 | C₄ `{v₂, v₋₂}` circle; D₃ at `z = ±i√10/2`, saddles; D₂ at `z = 0`, saddle |
| 400 | zonal `v₀` | one-parameter | D₃ at `z = ∞`, in-line maximum; D₂ at `z = ∞`, in-line saddle |
| 463 | hexagon | 2, 3, 6 | C₆ circle; D₃ at `z = 0`, maximum; D₂ at `z = ±√30/3`, maxima |

## The five cyclic lines

With `u = v_a + r·e^{iφ}·v_b` and `s = |c_a|²`, the restriction is a quadratic in `s`, independent of the relative phase. The endpoints are the weight states; an interior vertex is a circle of critical points, which is one orbit.

| line | `q(s)` | vertex |
|---|---|---|
| C₃ `{v₂, v₋₁}` | `−(15/44)s² + (3/22)s + 75/308` | `s* = 1/5`, the C3 ray |
| C₄ `{v₃, v₋₁}` | `−(8/33)s² + 75/308` | `s* = 0`, an endpoint, so no interior critical point |
| C₄ `{v₂, v₋₂}` | `−(12/11)s² + (12/11)s + 3/77` | `s* = 1/2`, the octahedron |
| C₅ `{v₃, v₋₂}` | `−(125/132)s² + (10/11)s + 3/77` | `s* = 12/25`, the pyramid |
| C₆ `{v₃, v₋₃}` | `−2s² + 2s + 1/924` | `s* = 1/2`, the hexagon |

The C₅ row is the paper's pyramid-line quadratic, up to the flip, and it is checked as a control. The C₄ `{v₃, v₋₁}` row is where Kawaguchi and Ueda's C₄ phase J would sit. Its linear coefficient is **exactly zero**, not small, so the vertex falls on the endpoint `s = 0` and the line contributes no interior orbit. The phrase "at the M8 couplings" is the right one: `r̂₆` is one member of the four-parameter family, and their phase J lives elsewhere in it. Every leading coefficient is negative, so each interior circle is an in-line maximum.

## The two dihedral lines

The chart is `u = b₁ + z·b₂` with `z = x + iy`, plus the point `z = ∞`, which is `v₀`.

- **D₃, `b₁ = v₃ + v₋₃`, `b₂ = v₀`:** `f = (100|z|⁴ − 20x² + 148y² + 463) / (231(|z|² + 2)²)`, which is the paper's § 5.8 formula (control). Its critical points are the hexagon at `z = 0` (maximum), the prisms at `z = ±√230/10` (minima), the octahedra at `z = ±i√10/2` (saddles), and the zonal at `z = ∞` (maximum). The indices sum to `1 + 1 + 1 − 1 − 1 + 1 = 2`.
- **D₂, `b₁ = v₂ + v₋₂`, `b₂ = v₀`:** `f = 4(25|z|⁴ + 142x² + 30y² + 72) / (231(|z|² + 2)²)`. Its critical points are the octahedron at `z = 0` (saddle), hexagons at `z = ±√30/3` (maxima), the D2 rays at `z = ±i√30/5` (minima), and the zonal at `z = ∞` (saddle). The indices sum to `−1 + 1 + 1 + 1 + 1 − 1 = 2`.

The hexagons on the D₂ line are the ones L0's exposure said were forced by `D₂ ⊂ D₆`. They also match Kawaguchi and Ueda's remark that their D₂ calculation returned a state of D₆ symmetry.

**Completeness on each dihedral line rests on the elimination route.** Each `f` is even in `x` and in `y` separately, so `f_x = x·A(X, Y)` and `f_y = y·B(X, Y)` with `X = x²` and `Y = y²`. That leaves four exhaustive cases:

1. `x = 0`, `y = 0`: the origin;
2. `y = 0` and `A(X, 0) = 0` with `X > 0`: a real-root count;
3. `x = 0` and `B(0, Y) = 0` with `Y > 0`: a real-root count;
4. `A = B = 0` with `X, Y > 0`: a resultant, which has no positive solution on either line.

Sympy's polynomial solve corroborates it and returns the same five chart points, but it is a tool with no completeness guarantee on a polynomial system, whereas the case analysis is an argument. So the elimination route carries the claim.

**The Poincaré–Hopf sum is a consistency check, not part of the completeness argument**, since a missed maximum and saddle would cancel in it.

## The union modulo rotation

Critical points are collected line by line and then united modulo rotation, as F1 required, since one orbit sits on several lines. The hexagon appears on the C₆, D₃ and D₂ lines, the octahedron on the C₄ `{v₂, v₋₂}`, D₃ and D₂ lines, the zonal on both dihedral lines, and each weight state at two or three cyclic endpoints.

Two values are shared by different orbits, and both ties are separated:

- **225:** `v₁` and the D2 ray. Their invariant triples `(|f|², |a₀₀|², TrN̄²)` differ, `(1, 0, 123/2)` against `(0, 1/112, 411/8)`, and so do their isotropy orders. This tie is a coincidence of one linear functional at the M8 couplings, not a structural identity.
- **`1188/5`:** the C3 ray and the pyramid. Their invariant triples are **equal**, which is the structural degeneracy the S0 memo recorded, so every invariant quartic agrees on them and the separation is by isotropy: orders 3 against 5.

## Controls and arms

- The paper's D₃ chart formula and the pyramid quadratic reproduce.
- With `u ⊗ u` in place of `u ⊗ Θu`, the D₃ formula differs.
- Every point found is critical on the whole sphere, which Palais' principle guarantees and the run checks: the largest tangential gradient over all 26 points is `1.6e−10`. A non-critical point on the C₃ line fails that check at `2e−1`.
- Dropping a point breaks each index sum.
- The union merges two rotated copies of the hexagon.
- Every exact value agrees with the numerical value at its point.
- 43 checks pass, 0 fail.

## Exposure for step 3

The in-line data already fixes part of step 3's answer, and the pre-registration's disclosure will say so.

- The zonal is an in-line maximum on the D₃ line and an in-line saddle on the D₂ line. The two lines meet `v₀` in different directions, `v₃ + v₋₃` and `v₂ + v₋₂`, and neither is an orbit direction, since the zonal's orbit tangents are the phase and the `m = ±1` rotation directions. So in the **transverse tangent space** it has both a negative and a positive direction: it is a saddle, whatever the exact signature turns out to be. That is an H-level fact, visible here.
- The prism is an in-line minimum, the octahedron an in-line saddle on both dihedral lines and an in-line maximum on its cyclic line, and the hexagon an in-line maximum wherever it appears.
- These agree with M8.11's published in-locus signs: `(−, −)` at the hexagon, `(−, +)` at the octahedron and `(+, +)` at the prism.
- On the cyclic lines every interior circle is an in-line maximum, with the orbit direction null. The C₅ case is M8.11's `−104/55` along `e_t`.

## Next

Step 3 (H): the exact transverse signature at each of these ten orbits, by `LDLᵀ` inertia, with the characteristic polynomial and the isotypic blocks, and the two global-extremum arguments as the ends.

## Review (2026-09-19)

Both units passed step 2 and its mathematics independently: F1 rederived both dihedral gradient factorizations and the absence of off-axis points, F2 checked the cyclic quadratics, the D₂ formula at all four of its special points and the `v₀` saddle argument. Their comment-tier fixes are applied above: the transverse tangent space (F1), the elimination route named as the carrier of completeness with its four cases written out (F2), the 225 tie labeled a coincidence (F2), and the exactly-zero linear coefficient on the C₄ line (F2).
