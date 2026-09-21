# S1 derivation, step 3 of 3: H, the transverse signatures

2026-09-19. Author-side and local, for the units' review. This is the step carrying the values S1 would freeze.

## The object, and why the Hessian is this simple

For a unit critical `u`, `H_u(e, e) = d²/ds² r̂₆(u cos s + e sin s)` at `s = 0`, on `T_u = {e : Re⟨u, e⟩ = 0}`.

`r̂₆ = N/D²` is homogeneous of degree 0, with `N(x) = ‖ρ₆‖²` and `D = ‖u‖²`. A degree-0 function has zero radial derivative, so at a critical point on the sphere its full gradient vanishes, and the great-circle second derivative is the ambient Hessian. With `‖u‖ = 1`,

`Hess(N·D^{−2})(u) = Hess N(u) − 4N(u)·I − 8N(u)·u uᵀ`,

and the last term drops on vectors orthogonal to `u`. So

**`H_u = Hess N(u) − 4N(u)·I` on `T_u`.**

`O_u = span{iu, −iJ_x u, −iJ_y u, −iJ_z u}` is null for it, and `N_u = T_u ⊖ O_u` has dimension 10 at the weight states and 9 elsewhere.

## The census

| orbit | `924·r̂₆` | `dim N_u` | `(n₋, n₀, n₊)` | index, `g > 0` | index, `g < 0` |
|---|---|---|---|---|---|
| coherent `v₃` | 1 | 10 | (0, 0, 10) | 0 | 10 |
| `v₂` | 36 | 10 | (2, 0, 8) | 2 | 8 |
| prism | `8800/43` | 9 | (3, 0, 6) | 3 | 6 |
| `v₁` | 225 | 10 | (4, **2**, 4) | 4 | 4 |
| D2 ray | 225 | 9 | (5, 0, 4) | 5 | 4 |
| C3 ray | `1188/5` | 9 | (5, 0, 4) | 5 | 4 |
| pyramid | `1188/5` | 9 | (5, 0, 4) | 5 | 4 |
| octahedron | 288 | 9 | (6, 0, 3) | 6 | 3 |
| zonal `v₀` | 400 | 10 | (8, 0, 2) | 8 | 2 |
| hexagon | 463 | 9 | (9, 0, 0) | 9 | 0 |

The signature is sector-independent by D7; in a sector the form is `w₆(σ)` times this one.

**What the table says.**

- **The two global extrema are transversely nondegenerate, modulo phase and rotations.** The coherent orbit is a minimum with no flat direction in `N_u`, `(0, 0, 10)`, and the hexagon a maximum, `(9, 0, 0)`. The global arguments G1 and G2 could not settle that; the census does.
- **No other orbit in the class is a local extremum, for either sign of `g`.** Every other one has both signs present.
- **`v₁` is the only degenerate orbit**, with a two-dimensional kernel. The kernel is exactly the direction of the C₄ line through it, `span{v₋₃, i·v₋₃}`, which is the line whose restriction has its vertex at that endpoint (step 2). So the degeneracy is the line's, not an accident.
  - **Explanatory, not a frozen claim (F1):** those two directions are quartic-flat, not numerically flat. On that line the restriction moves by `Δr̂₆ = −(8/33)s²` with `s = |c|²`, so along the kernel the first change is `−(8/33)|c|⁴`.
- **The zonal is a saddle**, as step 2 predicted from its two in-line types, and the census puts it at `(8, 0, 2)`.
- **Three different orbits share the signature `(5, 0, 4)`:** the D2 ray, the C3 ray and the pyramid. The C3 ray and the pyramid also share every invariant quartic, so for them even the eigenvalues might have been expected to agree; they do not (see below).
- **The index grows with the value**, from 0 at the minimum to 9 at the maximum, with no inversion.

## The exact spectra

Characteristic polynomials of the transverse operator, factored (the log carries them in full):

| orbit | eigenvalues (multiplicity in brackets) |
|---|---|
| coherent `v₃` | `2/11 [2], 3/11 [2], 32/33 [2], 65/33 [2], 4 [2]` |
| `v₂` | `−8/33 [2], 12/11 [2], 4/3 [2], 20/11 [2], 24/11 [2]` |
| `v₁` | `−5/3 [2], roots of 363λ² − 374λ − 600 [2], 0 [2], 3/11 [2]` |
| zonal `v₀` | `−100/33 [2], −20/11 [2], −40/33 [2], −12/11 [2], 8/11 [2]` |
| octahedron | `−24/11 [3], −8/11 [3], 40/33 [3]` |
| hexagon | `−4, −23/11, −67/33 [2], −64/33 [2], −15/11, −10/11 [2]` |
| pyramid | `−104/55, roots of 1815λ² − 1012λ − 3360 [2], −24/55 [2], 4/165 [2]` |
| prism | `roots of 671187λ² + 438944λ − 339200 [2], −696/473, 560/1419 [2], 8/11, 920/473` |
| C3 ray | `−24/55, roots of a quartic [2]` |
| D2 ray | `roots of 31944λ³ + 46948λ² − 13530λ − 7125 [2], −7/11, 10/11, 64/33` |

The eigenvalues are rational at the weight states, the octahedron and the hexagon, and algebraic of degree 2, 3 or 4 elsewhere. The multiplicities are the isotypic pattern: threefold at the octahedron, twofold nearly everywhere, and simple in the directions the stabilizer fixes.

## Method, exactness, controls

- **Primary:** exact inertia, by symmetric elimination with a 2×2 block whenever a pivot vanishes. That is Sylvester's law, so the signature is basis-independent.
- **Secondary:** the exact characteristic polynomial of the transverse operator in an exact orthonormal basis.
- **Cross-check:** numerical eigenvalues at 30 digits. Exact and numerical signatures agree at all ten orbits.
- **Parents.** M8.11's in-locus values reproduce from the full transverse form: `−104/55` at the pyramid along the C₅-line tangent, and `920/473` and `8/11` at the prism along `τ_x` and `τ_y`. The sector bridge then reproduces M8.11's `L_T` values `−56/165` and `−21/110`.
- **Gates that can fail, and did.** At every orbit, the orbit directions must be annihilated exactly; at a non-critical point they are not (arm). Dropping the `−4N·I` term changes the hexagon signature (arm). The ends must match G1 and G2. The zonal must be a saddle. 29 checks pass, 0 fail.

**Two catches by the step's own gates.**

- The prism's `τ_y` control failed first time. My tangent convention projected out only the radial direction, while M8.11's `τ_y = i·τ_x` is also orthogonal to the phase direction. The control was right and the convention wrong. Both conventions give unit tangents, so only a parent value could tell them apart.
- The kernel gate at `v₁` failed first time, because it pointed at `v₃`. The canonical C₄ line ends at `v₋₁`, so at `v₁` the corresponding direction is `v₋₃`. Fixed, and the kernel is now identified exactly.

## The parity diagnostic

The indices at the four weight-state orbits for `g > 0` are 0, 2, 4 and 8. Their orbits are `S²`, `S²`, `S²` and `ℝℙ²`, with Euler characteristics 2, 2, 2 and 1, while every orbit with a finite stabilizer is a closed 3-manifold and contributes 0. The alternating sum is `2 + 2 + 2 + 1 = 7 = χ(ℂℙ⁶)`, exactly.

The `v₁` group orbit is transversely degenerate, with two kernel directions beyond symmetry. Whether those directions integrate into a larger critical component lies outside S1's partial classification, so the global Morse–Bott hypothesis is not established. The parity observation is therefore diagnostic only, and it is not adjudicated against L or H.

It is informative precisely at `v₁` (F2). That is the one orbit whose index is not rigid: a perturbation resolving its kernel could send it to 4, 5 or 6, and an odd value would break the sum.

## Two observations on the table

- **The index sequence skips 1 and 7.** For `g > 0` it runs 0, 2, 3, 4, 5, 5, 5, 6, 8, 9, monotone in the value with no inversion. The three deferred orbits sit at values 215.78, 238.02 and 238.36, which monotonicity would place between indices 3 and 6, so they do not fill either gap. Either the critical set has more in it than the deferred three, or monotonicity fails somewhere outside the class. Both are open, and the gaps are visible in the frozen table, so the pre-registration says this rather than waiting to be asked.
- **Three geometrically different orbits share `(5, 0, 4)`:** the D2 ray, the C3 ray and the pyramid. The last two share every invariant quartic, so agreement there is less surprising than at the D2 ray, whose invariants differ.

## What this does and does not say

It is the Morse index of the leading reduced quartic on the block, transverse to the symmetry orbit. It is not stability of a branch germ, which is not particle stability. No MODELS.md cell moves and M8.7's gate is untouched. The orbits with isotropy C2, the C3 plane, or none, including the three census orbits, stay outside the claim.

## Exposure

Part of this census was visible before the inertia computation ran, and the pre-registration states the split row by row, as M8.11's did.

| orbit | what was already visible | what the census adds |
|---|---|---|
| coherent `v₃` | a global minimum (G1), so `n₋ = 0` | the nullity, and the full spectrum |
| hexagon | a global maximum (G2), so `n₊ = 0`; M8.11's in-locus `(−, −)`; an in-line maximum on three lines (step 2); the exploratory transverse count 9 negative, 5 null of 14 ambient | the exact spectrum, and that the count is exact and nondegenerate |
| octahedron | M8.11's in-locus `(−, +)`; in-line saddle on both dihedral lines, in-line maximum on its cyclic line (step 2) | the full signature `(6, 0, 3)` and a threefold spectrum |
| prism | M8.11's in-locus `(+, +)`, that is `920/473` and `8/11`; in-line minimum (step 2) | the other seven directions, and `(3, 0, 6)` |
| pyramid | M8.11's `−104/55` along `e_t`, and the null orbit direction; in-line maximum (step 2) | the other eight directions, and `(5, 0, 4)` |
| zonal `v₀` | a saddle, from the two in-line types (step 2); the "saddle" label in the exploratory rounds, without a Hessian | the exact `(8, 0, 2)` |
| C3 ray, D2 ray | in-line types only (step 2) | the full signatures |
| `v₂`, `v₁` | nothing beyond block-operator spectra | everything, including `v₁`'s degeneracy |

So a blind agent reproducing the hexagon's `(9, 0, 0)` is partly confirming something already on the public record, while `v₂` and `v₁` are unexposed.

## Next

The three steps are done. What remains before filing is the pre-registration itself: the task doc and worklist, with this census as the frozen values, the disclosure inventory attached, and the blind-run handout in the fibre model.

## Review (2026-09-19)

Both units reproduced the census independently, F1 all ten signatures from a fresh Hessian and F2 three full spectra plus the Hessian identity, and both passed it. Their fixes are applied above: the Morse–Bott wording and the quartic-flat explanation (F1), "transversely nondegenerate modulo phase and rotations" (F1), the exposure section (F2), and the two observations on the table (F2).
