# S1 derivation, step 1 of 3: L0, the fixed-locus classification

2026-09-19. Author-side and local. This is the first of the three steps that settled outline revision 3 names (L0, then L1, then H); each is reviewed by the units before the next starts. Step 1 is pure representation theory on `V₃`. It computes no Hessian and uses no critical data of `r̂₆`, except that the invariants separating the lines include `‖ρ₆‖²` (see *Exposure*).

## Statement

Let `H` be a closed subgroup of `SO(3)` and `χ` a one-dimensional character of `H`. Write

`V₃^{(H,χ)} = {u : D³(h)u = χ(h)u for all h in H}`.

This is the fixed space of the graph subgroup `{(χ(h)⁻¹, h)}` of `U(1) × SO(3)` acting by `u ↦ χ(h)⁻¹·D³(h)u`, so Palais' principle applies to it in step 2.

**L0.** Up to rotation, the fixed spaces of complex dimension 1 and 2 (projective points and lines) are exactly:

- **six points:** `v₃`, `v₂`, `v₁`, `v₀`, the octahedron `(v₂ + v₋₂)/√2` and the hexagon `(v₃ + v₋₃)/√2`;
- **five cyclic lines:**
  - C₃ on `span{v₂, v₋₁}`;
  - C₄ on `span{v₃, v₋₁}`;
  - C₄ on `span{v₂, v₋₂}`;
  - C₅ on `span{v₃, v₋₂}`;
  - C₆ on `span{v₃, v₋₃}`;
- **two dihedral lines:**
  - D₃ on `span{v₃ + v₋₃, v₀}`;
  - D₂ on `span{v₂ + v₋₂, v₀}`.

Every other nonzero character-fixed space of a finite subgroup has complex dimension at least 3, that is, projective dimension at least 2. There are four, from three groups: `C₁`'s (7), `C₂`'s two (4 and 3), and `C₃`'s `span{v₃, v₀, v₋₃}` (3).

## Proof

**Setup.** `R_z(θ)` acts by `v_m ↦ e^{−imθ}v_m`, and the half-turn `s = R_x(π)` by `v_m ↦ −v₋ₘ`. Closed subgroups of `SO(3)` are finite (cyclic, dihedral, T, O, I) or conjugate to `SO(2)`, `O(2)` or `SO(3)`.

**Cyclic groups.** For `C_n = ⟨R_z(2π/n)⟩`, the character spaces are `span{v_m : m ≡ k (mod n)}`, with these sizes:

| n | sizes |
|---|---|
| 1 | 7 |
| 2 | 4, 3 |
| 3 | 3, 2, 2 |
| 4 | 2, 2, 2, 1 |
| 5 | 2, 2, 1, 1, 1 |
| 6 | 2, 1, 1, 1, 1, 1 |
| 7 or more | all 1, since `m` in `[−3, 3]` are distinct mod `n` |

The lines are:

- `{2, −1}` and `{1, −2}`, for `n = 3`;
- `{3, −1}`, `{1, −3}` and `{2, −2}`, for `n = 4`;
- `{3, −2}` and `{2, −3}`, for `n = 5`;
- `{3, −3}`, for `n = 6`.

The flip `s` exchanges the members of each of the pairs `{2, −1}` and `{1, −2}`, `{3, −1}` and `{1, −3}`, and `{3, −2}` and `{2, −3}`, which leaves five cyclic lines.

**Dihedral groups.** Take `D_n = ⟨r, s⟩` with `r = R_z(2π/n)`. A one-dimensional character has `χ(r) = χ(srs⁻¹) = χ(r)⁻¹`, so `χ(r) = ±1`, and `χ(r) = 1` when `n` is odd. Its space is therefore the `χ(s)`-eigenspace of `s` inside `span{v_m : m ≡ 0}` or `span{v_m : m ≡ n/2}` (mod `n`).

| n | how the spaces split |
|---|---|
| 2 | even `m`: `v₂ − v₋₂` and `span{v₂ + v₋₂, v₀}`; odd `m`: `span{v₃ − v₋₃, v₁ − v₋₁}` and `span{v₃ + v₋₃, v₁ + v₋₁}`. That is 1 + 2 + 2 + 2. |
| 3 | `{3, 0, −3}` splits into `v₃ − v₋₃` and `span{v₃ + v₋₃, v₀}` |
| 4 | `v₀` and `v₂ ± v₋₂`, all points |
| 5 | `v₀` only |
| 6 | `v₀` and `v₃ ± v₋₃`, all points |
| 7 or more | `v₀` only, since `n/2 > 3` |

The 3-fold rotation about `(1, 1, 1)` normalizes `D₂` and permutes its three nontrivial characters, so the three D₂ lines form one class. That leaves two dihedral lines.

**T, O, I.**

- **T ≅ A₄.** Its one-dimensional characters factor through `T/D₂ ≅ C₃`. The character formula gives multiplicities 1, 0 and 0, so there is one point. It is `D₂`-fixed with trivial `D₂` character, so it is `D₂`'s point `v₂ − v₋₂`, the octahedron.
- **O ≅ S₄.** Its characters are trivial (multiplicity 0) and the sign `A₂` (multiplicity 1): the octahedron again.
- **I ≅ A₅.** It is perfect, so it has only the trivial character, with multiplicity 0, since `V₃|_I = 3′ ⊕ 4`. There are no fixed spaces.

**Continuous groups.** `SO(2)` gives the weight points `C·v_m`; `O(2)` gives `v₀`; `SO(3)` gives nothing, since `V₃` is irreducible and nontrivial.

**Identifications and distinctness.**

- Points: `v_m` and `v₋ₘ` are related by `s`. `v₃ − v₋₃` is `R_z(π/6)` of the hexagon, up to phase. `v₂ − v₋₂` is `R_z(π/4)` of the octahedron. T's and O's points are the octahedron.
- The six points are distinct: their invariants `(|f|², TrN̄², |a₀₀|²)` differ. They are `(9, 171/2, 0)`, `(4, 48, 0)`, `(1, 123/2, 0)`, `(0, 72, 1/7)`, `(0, 48, 1/7)` and `(0, 171/2, 1/7)`.
- **The seven lines are distinct, exactly.** For a line `L`, let `Fix(L)` be the rotations whose `D³` acts on `L` as a scalar. A rotation carrying `L₁` to `L₂` conjugates `Fix(L₁)` onto `Fix(L₂)` and carries one character to the other.
  - Every cyclic line contains a weight ray `[v_a]` with `a ≠ 0`. Its stabilizer is `SO(2)_z`, because it must fix `⟨f⟩ = a·ẑ`. So `Fix(L)` lies in `SO(2)_z`, and the scalar condition gives C₃, C₄, C₄, C₅ and C₆.
  - Every dihedral line contains the zonal ray `[v₀]`, whose stabilizer is `O(2)_z`. The scalar condition gives D₃ and D₂.
  - These groups are pairwise non-conjugate, except the two C₄'s. On those, the generator `R_z(π/2)` acts by `i` on `span{v₃, v₋₁}` and by `−1` on `span{v₂, v₋₂}`. The normalizer of C₄ conjugates the generator to its inverse, which sends `i` to `−i` and fixes `−1`. So the character classes `{i, −i}` and `{−1}` differ.
  - Independently, `|f|² = (4s − 1)²` on the first C₄ line (range `[0, 9]`) and `4(2s − 1)²` on the second (range `[0, 4]`).

  The grid table of multipole ranges in the log corroborates this. It is also the source of the exposure below. □

## The computational check (`s1_L0.py`)

The check covers the groups `C₁…C₁₂`, `D₂…D₁₂`, T, O and I, each generated by closure with its order verified. The general argument above covers `n > 12` and the continuous groups.

- **Route 1, characters.** Each group's one-dimensional characters are found by closure-consistency on its generators, trying only roots of unity of each generator's order. Multiplicities come from `χ_V₃(θ) = 1 + 2·Σ_{m=1}^{3} cos(mθ)`.
- **Route 2, joint eigenspaces.** Common eigenspaces of the generators' spin-3 matrices, computed explicitly.
- **The routes agree on all 26 groups.** The space dimensions are the tables above, for example `[1, 2, 2, 2]` for `D₂`, `[1, 2]` for `D₃`, `[1]` for T and O, and none for I.
- **Classification.** For each point and line from route 2, a numerical rotation search finds a rotation carrying it to its canonical representative. The log records the rotation vector and the residual on recheck, below `10⁻¹²`. That yields exactly six point classes and seven line classes. All three `D₂` lines go to the D₂ class, and T's and O's points to the octahedron. The only spaces of dimension 3 or more are `C₁`'s 7, `C₂`'s 3 and 4, and `C₃`'s 3.
- **Distinctness.**
  - The six points have pairwise distinct multipole-norm vectors.
  - The exact argument for the lines is checked directly. Within `O(2)_z`, the log counts the rotations about `z` (in 6° steps) and the horizontal half-turn axes (in 3° steps) that act as scalars on each line. It also checks the C₄ generator's action, `i` against `−1`, and the two `|f|²` formulas.
  - The grid table corroborates this. Its closest pair is D₃ against D₂, at 0.028, and they stay separated with `K = 6` left out (0.018, mainly on `‖ρ₄‖²`).
- **Arms:**
  - a spin-2 character table in route 1 must break agreement at `D₂` (`[1, 1, 1, 2]` against `[1, 2, 2, 2]`). Spin 2 is used because it is a genuine SO(3) representation, whereas a half-integer character is not a class function on SO(3);
  - a planted wrong line, `span{v₃ + v₋₃, v₁}`, must fail the D₃ counts;
  - with `D₂` removed, its line class must be lost (6 classes);
  - T alone must give no line;
  - the two C₄ lines must not be carried to each other;
  - a control asserts that a half-turn has eigenvalues `±1` only.
- **Two failures, both fixed.**
  - The first full run failed its own route-agreement gate at `C₂`, all dihedral groups and T. The axis-and-angle routine took `arccos` of the trace, which loses half its digits near a half-turn, so every half-turn came out as nearly the identity. It now uses scipy's rotation vector, and the half-turn control guards the fix.
  - The first run of the fixing-group check, added in review, counted 60 horizontal half-turn axes on three cyclic lines. The expected count was 0. Its predicate accepted a zero matrix as a scalar, which is what a half-turn gives when it maps a line to an orthogonal one. The predicate now requires a scalar of modulus 1, and an arm plants exactly that case. The gate failed on the wrong count, which is the job it was written for.

## Exposure

The line invariants include the range of `‖ρ₆‖² = r̂₆` on each line. So the log's table already shows, to grid accuracy, each line's minimum and maximum of `r̂₆`, which is step 2's territory:

- on the C₃ line, a maximum of `9/35`, the C3 ray;
- on the D₃ line, a minimum of `200/903`, the prism, which is already published;
- on the D₂ line, a minimum of `75/308`, the D2 ray, and a maximum equal to the hexagon's `463/924`.

**The last item is forced, not an observation.** `D₂ ⊂ D₆`, since `R_z(π/3)³ = R_z(π)` and `R_x(π)` lies in `D₆`. The hexagon's character restricted to `D₂` is `(−1, −1)` on `(R_z(π), R_x(π))`. So the hexagon lies on the D₂ line `span{v₃ + v₋₃, v₁ + v₋₁}`, which the 3-fold rotation carries onto the canonical one. Step 2 confirms it rather than determines it.

The octahedron and the zonal state also lie on the canonical D₂ line, at `z = 0` and `z = ∞`. The octahedron lies on the D₃ line and the C₄ line `span{v₂, v₋₂}` as well, and the hexagon on the C₆ and D₃ lines.

Step 2 derives the complete critical sets exactly. The pre-registration's disclosure lists this log.

## Files

`s1_L0.py` and `s1_L0_log.txt` are in `M8_S/`, with their hashes in the message that carries this note.

## Next

Step 2 (L1): the exact critical set of `r̂₆` on each of the seven lines, with the Poincaré–Hopf gate on the two dihedral lines and the quadratic on the five cyclic ones. As F1 requires, the critical points are reported line by line and then united modulo rotation, with a deduplication gate, before H counts orbits. The overlaps listed under *Exposure* make that gate necessary.

## Review (2026-09-19)

Both units passed the L0 result. Their fixes are applied above:

- "nonzero" in the statement (F1);
- the count, now four spaces from three groups (F2);
- the exact distinctness argument, with fixing groups and the C₄ character from F2, and the `|f|²` ranges from F1;
- the spin-2 arm (F1). F2 had praised the spin-7/2 arm's role, which spin 2 keeps without the conceptual defect;
- the rotation recorded per class (F1);
- the finer-grid claim dropped, since it was not in the script (F1);
- the hexagon on the D₂ line stated as forced (F2);
- the deduplication requirement carried into step 2 (F1).
