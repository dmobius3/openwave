# M8.15: THE SOFT-SLOT BRANCHES AT FINITE AMPLITUDE, a pre-registered continuation and persistence run under the frozen slot action

Pre-registration, for review. It governs nothing until the maintainer's go. **No M8.15 target has been run.** The solver exists, and only its controls have been run: the bench dry runs, the closing control runs that fix the values below, and a rehearsal of the full pipeline on the controls. The values marked "fixed at the freeze" are set from the closed control runs before the freeze and never after any target has run.

## TASK PLANNING

### The question

The author's slot action ([`slot-action.md`](https://github.com/dmobius3/mode-identity-theory/blob/341d78925dbba8682537a0b33057c531981e415a/files/framework/files/working/files/slot-action.md), §§ I to VIII frozen at commit `341d789`) gives the eight slot fields one law on `ℝ × S³/2I`. In its three soft slots, `R2`, `R4` and `R5`, the lowest level is not rigid, and the soft blocks' reduced quartic has critical orbits beyond its minimum. The author's [The Soft-Slot Branches](https://github.com/dmobius3/mode-identity-theory/blob/d0de8ca61f8e166586731b6167ab895c09368630/files/framework/files/working/files/soft-slot-branches.md) (commit `d0de8ca`) proves two things:
- every critical orbit that is nondegenerate modulo the stationary symmetry group continues to a branch of standing waves at small amplitude (Lemma G);
- the branch's small-amplitude linear spectrum is decided by the reduced quartic (Theorem R).

It also computes exactly, or certifies, the orbits named below.

This task asks one question: **do the standing-wave branches born at the named non-minimizing orbits continue to finite amplitude, and along them, do they stay linearly stable and pass a pre-registered nonlinear persistence test?** The small-amplitude facts are theorems and are not re-tested. The run's content is what mathematics leaves open:
- the finite-amplitude continuation;
- where, if anywhere, ellipticity is lost;
- the full linearized spectrum along each continued branch;
- nonlinear persistence.

### Standing

- **The ruling it files under** is [#512](https://github.com/openwave-labs/openwave/discussions/512#discussioncomment-18415036)'s.
  - A program that goes beyond the closed chassis says what separates its solver from it (§ THE SOLVER).
  - It pins its own pointwise law, and takes nothing from `m8_5c/design_inputs/`.
  - A new program gets a new identity, filed as a BACKLOG row with its task doc in the pull request that carries the pre-registration.
- **Run-before-write** (`PR_REVIEW_STANDARDS.md` § 12.2) is met. M8.14 ran and was adjudicated on 2026-10-01 ([#610](https://github.com/openwave-labs/openwave/pull/610)). Since then the column has landed only that record's dated note ([#611](https://github.com/openwave-labs/openwave/pull/611)) and living-docs syncs ([#613](https://github.com/openwave-labs/openwave/pull/613), [#617](https://github.com/openwave-labs/openwave/pull/617)).
- **Not the closed chassis.** The solver is a physical-space finite-element instrument. It uses no truncated harmonic basis, no Galerkin projection onto harmonics, and no chassis code or input.
- **Word budget.** About 6,800 words, inside § 12.2's cap of about 8,000.
- **ID.** M8.15 is the next free ID in creation order, allocated by the row this pull request creates (`ROADMAP_STANDARDS.md` § 6). If a row taking M8.15 lands first, whoever merges second renumbers this one.

### What this task is not

- **No slot relation.** M8.2's row 7 needs a filed relation, and none is filed, so the ceiling is row 6. Small-amplitude continuation carries M8.4's label "nonlinear persistence of the installed free structure". Whether finite-amplitude stability of 2I's non-minimizing branches earns more under M8.4's second ceiling is the adjudication's to decide, not this document's.
- **No mass, no particle state, no physical matter action.** The slot action is MOTIVATED, an effective model.
- **No re-test of theorems.** The small-amplitude existence and types are the author's theorems. A run disagreeing with them at the smallest amplitudes indicts the instrument, and the controls exist to catch exactly that.
- **No census claim.** The targets are named orbits. Whether the critical sets are complete stays open, and nothing here depends on it.
- **Records.** No `MODELS.md` cell moves, M8.7's gate is unchanged, and no earlier M8 verdict is reread.

### Ownership and run format

- **The author's agents build and run** the solver under the frozen terms, and publish the scripts, data and records with a reproduction route (`REPRODUCE.md`, row C6).
- **The frozen terms' SHA-256 lands on `main` before the run starts**, so the run can be checked against terms fixed in advance.
- **No obligation on the maintainers.** Reproduction is optional and declinable per run. A claim of a defect in the frozen text is reproduced by the maintainer with independent code before it is ratified, as in M8.14.
- **Adjudication** is the maintainer's. It scores the run against the claims and rules below.

### Sources of record

| Source | What it supplies |
| --- | --- |
| MIT [`slot-action.md`](https://github.com/dmobius3/mode-identity-theory/blob/341d78925dbba8682537a0b33057c531981e415a/files/framework/files/working/files/slot-action.md) §§ I to VIII at `341d789` | the frozen law: the field equation, the energy, the charge, the slots |
| the same page, § IX's dated notes, at [`d0de8ca`](https://github.com/dmobius3/mode-identity-theory/blob/d0de8ca61f8e166586731b6167ab895c09368630/files/framework/files/working/files/slot-action.md) | the fibre symmetry, O(2) in the real slots and Sp(1) in the quaternionic ones; stability claims stated relative to the whole stationary family, assessed point by point along each Sp(1) family |
| MIT [`soft-slot-branches.md`](https://github.com/dmobius3/mode-identity-theory/blob/d0de8ca61f8e166586731b6167ab895c09368630/files/framework/files/working/files/soft-slot-branches.md) at `d0de8ca` | the targets and controls: representatives, values, signatures, germs (Lemma G), small-amplitude types (Theorem R), the exact slow spectra, and the interval certificate of the third outside orbit |
| MIT [`slot-ground-states.md`](https://github.com/dmobius3/mode-identity-theory/blob/d0de8ca61f8e166586731b6167ab895c09368630/files/framework/files/working/files/slot-ground-states.md) at `d0de8ca` | Theorems A, B and D; the reduced quartic in closed form |
| [M8.12](https://github.com/openwave-labs/openwave/blob/344f4162eb2bc256728c240c31382db277f399b1/openwave/xperiments/m8_mit/research/tasks/m8_12_task_details.md) at `344f416` | the census orbits in `R4` |
| M8.2, M8.4 | the ladder (row 6 ceiling), and the persistence label |
| This task | the claims, controls, rules and frozen values below |

## SETTING

- **Units.** `R = c = g = 1` in every reported number.
- **The field equation.** `ψ` is a section of the slot bundle `E_ρ` over `X = S³/2I`, equivalently a `ρ`-equivariant field on `S³`, `ψ(γx) = ρ(γ)ψ(x)` for `γ ∈ 2I` acting on the left. It obeys the frozen equation `ψ_TT − Δψ + |ψ|²ψ = 0`.
- **Standing waves** are `e^{iωT}φ`, with `σ = ω²`.
- **The soft slots.** Their lowest levels are `K = 6` for `R4` and `R5` (`λ₀ = 48`) and `K = 7` for `R2` (`λ₀ = 63`). The next levels where each slot occurs are `K = 10` (`λ₁ = 120`), `K = 8` (`80`) and `K = 13` (`195`). The amplitude parameter is `ε = σ − λ̄₀`, with `λ̄₀` the mean of the slot's lowest discrete level at the level run, the solver's convention. In `R4` it exceeds `λ₀` by `2.5×10⁻²` at level 2 and `1.5×10⁻³` at level 3, so the levels compare at equal amplitude.
- **Exact symmetries of the stationary problem:**
  - right `SU(2)`;
  - the phase, and charge conjugation in `R4` and `R5`;
  - `Sp(1)` in `R2`, generated by `i` and the quaternionic structure `J`, with `J v = ε_{ab} conj(v)` on the fibre.

  Stability is stated relative to the stationary family, as the slot action's § IX requires.

## THE SOLVER, AND THE SEPARATION STATEMENT

**Geometry.**
- The 600-cell, whose 120 vertices are 2I itself.
- **Refinement.** It is refined canonically: each tetrahedron splits into its 4 red-refinement corners plus the inner octahedron cut from its projected centre into 8. No choice is made, so the mesh is invariant under left 2I (the bundle, exactly) and under right 2I. Red refinement cannot keep right 2I, because each regular tetrahedron's symmetry permutes the three octahedron diagonals.
- **Elements.** Exact-geometry P2 finite elements: flat tetrahedra radially projected onto `S³`, with the metric evaluated exactly at the quadrature points. The rule is a collapsed Gauss rule averaged over the vertex permutations, so it commutes with each element's symmetries.
- **Unknowns.** Values at one canonical point per left-orbit of nodes. Elements are integrated once per left-orbit, with blocks `K_ij ρ(g_i)†ρ(g_j)`. Node orbits number 84, 1008 and 12096 at levels 1, 2 and 3.

**The law.**
- The stiffness `K`, the consistent mass `M`, and the quartic `¼∫|ψ|⁴` by the same quadrature, with its exact gradient and real-linear Jacobian.
- The semi-discrete system is `MΨ'' + KΨ + F(Ψ) = 0`, Hamiltonian, with the phase, `J` (in `R2`) and right 2I exact.

**Standing waves, and the continuation algorithm.**
- Newton at fixed `Φ†MΦ = m`, with `σ` unknown, bordered by gauge constraints for the exact continuous fibre symmetries. The gauge multipliers vanish at a solution. Newton stops at a relative residual of 10⁻¹¹.
- **Continuation is in `m`, and `ε = σ − λ̄₀` is an output.** From the accepted point at `ε_{k−1}`, the `m` for `ε_k` is predicted by a secant through the last two accepted points (by `m ∝ ε` at the first step). Newton is run at that `m`, and the secant is iterated on `m` until `|ε − ε_k| ≤ 10⁻⁶ ε_k`, at most 8 iterations.
- **Accepted solutions.** Every converged Newton solution is an accepted solution, secant iterates included, and its `(m, ε)` is recorded. If 8 secant iterations do not reach `ε_k`, the fold search below runs only if the accepted points meet the fold criterion. Otherwise the branch ends as a secant target failure, at its last accepted point.
- **A fold** is declared at the first accepted solution whose `ε` is smaller than the previous accepted one while `m` increased. The turning point is then located by a golden-section search in `m` that maximizes `ε(m)` within the bracket of the last three accepted points, to 10⁻³ relative in `m`. The branch ends there, and the turning point is recorded.
- **Newton failure.** The step in `m` from the last accepted point is retried at one half, then one quarter. A failure after both ends the branch.
- **The seed rule.**
  - At `ε_min`: the oriented nodal interpolant of the germ direction, scaled to the leading-order norm, solved inside the pinned fixed space.
  - At every later ladder point and continuation step: the previous accepted solution, scaled to the new `m`.
  - There are no other seeds, and the method is never switched for a branch.

**Pinned seeds.**
- The mesh keeps right 2I but not right `SU(2)`, so a continuum critical orbit breaks into isolated discrete critical points, and Newton is nearly singular along the lifted rotation directions.
- Each representative is oriented so that the subgroup in the table below lies in the mesh's right 2I, with the matching phase or fibre map. Newton is solved inside that subgroup's fixed space.
- Where a one-parameter near-degeneracy survives (the third outside orbit, rotations about its 2-fold axis), the step along it is constrained by a slice, and the residual force along it is recorded with the floors.
- Seeds are the nodal interpolants of the oriented germ directions.

**The linearization.**
- The rotating frame, `M(ζ'' + 2iωζ') + Lζ = 0`, with `L` the stationary Jacobian, in real first-order form.
- Slow eigenvalues by shift-invert at shifts off the exact zeros, each kept only with a residual below the threshold.
- The negative index of `L`, which is that of the phase-space form `B`, for the Krein count below.

**Time integration.**
- **The scheme.** Störmer–Verlet with the consistent mass, factorized once.
- **The start.** The exact discrete standing wave: the two-step scheme carries `e^{inθ}Φ` exactly for `θ = 2 arcsin(dt √σ / 2)`, given the initial velocity `Π₀ = (e^{iθ} − 1 + dt²σ/2) Φ / dt`.
- **The charge** `Im(Ψ†MΨ')` is a bilinear invariant, conserved to round-off.

**Separation statement.**
- No truncated harmonic basis, no Galerkin projection onto harmonics, and nothing from the closed chassis or `m8_5c/design_inputs/`.
- The finite-element weak formulation is a physical-space Galerkin discretization. What the closed-route ruling excludes is truncation in harmonic space, or Galerkin projection onto a harmonic basis.
- Harmonics appear only in the validation reference (the exact levels `K(K+2)`) and in the seed interpolants.

**Stated deviation from the first design.** The first design named a lumped mass. P2 row-sum lumping on tetrahedra gives nonpositive vertex weights, so the consistent mass is used. This is an instrument correction made before any target ran.

**The instrument, fixed.** The solver runs single-threaded, with fixed ARPACK starting vectors, so it is deterministic under the frozen environment that the manifest records: the library versions, the BLAS and ARPACK builds, and the thread counts.

- The linear record (statics, spectra, Krein counts) runs at level 3.
- Dynamics run at level 2, since level 3's time steps are impractical.
- The cross-level checks below are stated in advance.

## FROZEN TARGETS

| ID | slot | orbit (representative, Condon–Shortley `v_m`) | `(n₋, n₀, n₊)` | small-amplitude type (Theorem R) | pinned subgroup |
| --- | --- | --- | --- | --- | --- |
| T1 | R4, R5 | `v₂` | (2, 0, 8) | spectrally stable, (e) | `C₅` about a 5-fold axis |
| T2 | R4, R5 | prism, `v₃ + √(23/10) v₀ + v₋₃` | (3, 0, 6) | spectrally stable, (e) | `D₃`: the 3-fold about z; the half-turn about x, with phase −1 |
| T3 | R4, R5 | hexagon, `v₃ + v₋₃` | (9, 0, 0) | spectrally stable, (d) | `D₃ ⊂ D₆`, as for T2 |
| T4 | R4, R5 | the third outside orbit, interval-certified | (7, 0, 2) | spectrally stable, (e) | `C₂` about its `J` axis, with a slice |
| T5 | R2 | the member with complex structure `i` of the orbit of `√(3/10) v_{7/2} + √(7/10) v_{−3/2}`, `Q = 22/13` | (9, 0, 0) | spectrally stable at small amplitude, (d) | `D₅`: the 5-fold about z with a phase; the 2-folds at 18° + 36°k in the xy-plane, y among them, with `J` on the fibre |

Each target counts once per slot, so there are nine branches. In the `h` normalization, `R5`'s small-amplitude slow spectra are `R4`'s scaled by 9/16 (OpenWave's D7). Its `Q`, and so its `κ`, are its own (§ FROZEN VALUES). At finite amplitude the two slots run separately.

**T5 is one member.** Its Sp(1) family holds distinct standing waves, whose finite-amplitude spectra may differ point by point (the slot action, § IX). Theorem R(d) makes the whole family elliptic at small amplitude, but the run tests only the member named, and its verdict is not extended to the other members.

## CONTROLS (armed: a failed control stops the run before any target is scored)

| ID | case | known answer | what it shows |
| --- | --- | --- | --- |
| C1 | `R3`'s exact branch at `ε = 6` (`v` a weight vector; pin `C₅`) | `ω² = λ̄₀ + \|v\|²`, with the discrete lowest level (Theorem B); spectrum `±2iω` twice and `±i√(4ω² + 2\|v\|²)` (Theorem D); Persists under the frozen persistence protocol | the instrument reproduces exact branches and their spectra, and measures paired lifted modes for the floor |
| C2a | `R4`'s coherent branch, `ε ∈ {0.5, 1, 2}` (pin `C₅`) | slow frequencies `εκτ`, `κ = 1/(8ω₀Q)`, `τ = 28/143, 56/429, 140/99, 112/39, 896/1287`; elliptic; persists (Theorem R(c)) | soft branches, and the small-amplitude reduction |
| C2b | `R2`'s weight-7/2 branch, the member with complex structure `i`, the same `ε` (pin `C₅`) | `τ = 112/143 … 560/143`; elliptic; persists | the Sp(1) handling: the exact `J` and `iJ` zero modes, and the gauge |
| C3 | `R4`'s octahedron, the same `ε` (pin `T`, the tetrahedral group, inside both the octahedral and the icosahedral groups; `v₂ + v₋₂` is `z(x² − y²)`, so `T`'s 2-fold axes are z and `(1, ±1, 0)/√2`) | a real pair at `εκ√ν`, `ν = 250880/552123`, three times; the perturbation grows at that rate (Theorem R(b)) | the instrument sees an instability when there is one; it also carries radical lifted modes (`J = 0`) |

**Gates.**
- **Linear, C2a, C2b and C3.** At `ε = 0.5`, every tracked nonzero slow frequency, and C3's rate, is within 5% of the leading order. With `r_ε` the ratio computed/leading-order, `|r_{0.5} − 1| < |r_1 − 1| < |r_2 − 1|` for each of them.
- **Linear, C1:** the frequency law to 10⁻⁶ relative, and Theorem D's eigenvalues to 10⁻⁶.
- **The cluster, C1 to C3.** At every control point, `G_tol` passes the cluster checks (§ THE LINEAR VERDICT), with C1's rigid-family pair included, and every lifted mode meets its bound.
- **Dynamic, under the persistence protocol below, at `ε = 1`, and C1 at its `ε = 6`:**
  - C1, C2a and C2b give Persists;
  - C3 gives Fails, and its `d_rot` at the Fails observation is reported;
  - **C3's rate.** A separate C3 run with `η = 10⁻⁷`, along the same seeded direction, is fitted on `d_⊥ ∈ [10η, 10⁴η]`. The fit must lie within 5% of the computed level-2 linear rate at the same point. The linear gate above separately compares that computed rate with Theorem R's leading order.

None of the known answers leans on the targets.

## NUMERICAL FLOORS AND TOLERANCES (fixed at the freeze, from the closed control runs)

- **Residuals.** Newton: relative residual ≤ 10⁻¹¹. Eigenpairs: backward error `‖Ax − λBx‖ / ((‖A‖₁ + |λ| ‖B‖₁) ‖x‖) ≤ 10⁻¹⁰`.
- **The floor, `r_f⁽ʰ⁾`,** at level `h = 3` for the linear record and `h = 2` for the dynamics: the largest modulus among the paired lifted modes (below) of C1, C2a and C2b at that level. Closing runs, pinned: `r_f⁽³⁾ = 1.040×10⁻⁶` (C2a's pair at `ε = 0.5`; C2b 6.8×10⁻⁷; C1 below 10⁻¹²) and `r_f⁽²⁾ = 3.970×10⁻⁵` (C2b's pair at `ε = 0.5`; C2a 9.2×10⁻⁶). All are imaginary pairs.
- **The zero tolerance,** `τ₀⁽ʰ⁾ = 10 r_f⁽ʰ⁾`: the displacement from 0 that the identified paired modes may show. It never defines the cluster by itself (§ THE LINEAR VERDICT).
- **The ellipticity tolerance,** `τ_Re⁽ʰ⁾ = 10 r_f⁽ʰ⁾`. A candidate instability needs a real part above `τ_Re⁽³⁾` and level stability (below).
- **Which level's values apply.** The linear record uses `τ₀⁽³⁾` and `τ_Re⁽³⁾`. The level-2 cluster, which the persistence perturbation is projected off, uses `τ₀⁽²⁾`. The dynamics margin is `r_f⁽²⁾ T`.
- **Two kinds of lifted mode.** The mesh breaks right `SU(2)` to 2I. A rotation direction lifts in one of two ways, read from the Result's `σ`-structure.
  - **Paired.** The direction is `σ`-paired with another, through `J(u) ≠ 0`. It lifts linearly in the mesh's breaking, and `τ₀` bounds it. This covers `v₂`, two of the third outside orbit's three directions, C1, C2a and C2b.
  - **C1's rigid family.** It adds two lifted directions: those of the family's tangent that the symmetry orbit does not span. They are paired with each other, and `τ₀` bounds them.
  - **Radical.** The direction lies in the `σ`-radical and has a Jordan partner. The Jordan pair splits at about the square root of the breaking, and grows with `ε`. This covers all three directions at the prism, the hexagon, T5 and C3, and the third outside orbit's direction along `J`.
  - **At C3's pinned points, level 2,** the radical modes are imaginary pairs along the whole ladder, from `7.3×10⁻⁴` at `ε = 0.5` to `4.5×10⁻²` at `ε = 72`. That is two to three orders above C2a's paired pair, and at most 0.15 of C3's slow rate. At level 3 they run from 2.3×10⁻⁴ (`ε = 0.5`) to 5.5×10⁻⁴ (`ε = 2`), still imaginary. They shrink about 3.2 times per level, near the square root of the mesh splitting's 9.1. Unpinned branches show real pairs of the same size.
  - **The radical bound.** A radical mode is not held to `τ₀`. At either level, it must stay below half the largest modulus among the branch's frozen slow eigenvalues at that `ε`, `½ ε κ τ_max`. Otherwise the point is Unresolved. Its modulus and real part are reported.
- **Resolved modes, decided at the freeze.** A slow mode is resolved at `ε` if `τ₀⁽³⁾ ≤ ε κ τ / 10`, from the frozen spectra with its branch's own `κ`: ten times the zero tolerance.
  - Every branch starts at `ε_min = 0.25`. A mode not resolved there is left out of the validity check at `ε_min`.
  - **No early promotion.** Until a mode's pre-registered first-resolved point, its branch's local type is Unresolved at every ladder point, whatever the mode's measured position. The measured spectrum is reported. From that point on, the measured mode takes over. Target output cannot promote a mode early.
  - The freeze block lists each unresolved mode and the ladder point where it is first resolved, before any target runs.
  - **The binding mode** is the third outside orbit's smallest, `τ ≈ 5.74×10⁻³`. At `ε = 0.25` it is `2.2×10⁻⁵` in `R4` and `1.3×10⁻⁵` in `R5`, resolved only if `r_f⁽³⁾ ≤ 2.2×10⁻⁷` and `1.3×10⁻⁷`.
  - Next is `v₂`'s Krein-negative frequency, `1.9×10⁻⁴` and `1.1×10⁻⁴` at `ε = 0.25`, resolved if `r_f⁽³⁾ ≤ 1.9×10⁻⁶` and `1.1×10⁻⁶`.
  - Every other target mode is resolved at 0.25 if `r_f⁽³⁾ ≤ 3.2×10⁻⁶`, the bound set by the third outside orbit's second-smallest mode in `R5`.
  - **At the closing floor** `r_f⁽³⁾ = 1.040×10⁻⁶`, every target mode is resolved at `ε = 0.25` except the third outside orbit's smallest, which is first resolved at `ε ≈ 1.41` in `R4` and 2 in `R5`. `v₂`'s Krein-negative mode in `R5` clears its bound, 1.09×10⁻⁶, by 5%.
- **The reference floor.** The unperturbed reference's distance from its stationary family must stay at most `η_ref = η` up to the run's first event, or `T` if there is none. Otherwise that test is Invalid, never Fails.
- **Conservation floors.**
  - The relative charge drift stays ≤ 10⁻¹⁰ over every run.
  - **The energy** is read against the perturbation's own energy `E_η = E(0) − E_ref`, which the fast half keeps positive. `E_η` must be strictly positive at initialization; otherwise the test is Invalid (energy normalization), with no redraw and no reweighting.
    - **The oscillation:** `max |E(t) − E(0)| ≤ 0.25 E_η`. Verlet's bounded fluctuation of the fast half is about 0.09 of it at the registered `dt`.
    - **The secular drift:** the mean of `E` over the last tenth of the observations, minus its mean over the first tenth, is at most `10⁻² E_η` in modulus. A correct Verlet run has none.
    - `10⁻⁶` of the total would fail a correct run at `η = 10⁻³`.
  - These are read up to the run's first event, or `T` if there is none. A run that breaks any of them there is invalid; what it does after its verdict is decided cannot change it. The floors are armed: a control whose run breaks one fails its gate.

## THE AMPLITUDE LADDER

- **The grid.** `ε_k = 0.25 × 2^{k/2}`, every grid value from `ε_min` up to `ε_max = λ₁ − λ₀`, with `ε_max` itself appended as the last requested point.
  - `ε_min = 0.25` for every branch. Modes the level-3 floor does not resolve there are handled per mode (§ NUMERICAL FLOORS).
  - `ε_max` is a fixed cutoff equal to the slot's exact continuum gap `λ₁ − λ₀`: 72 in `R4`, 32 in `R5` and 132 in `R2`. At a finite level it is a common cutoff for both levels, not the exact place of the mesh's next discrete level. Slow–fast interactions start in that range.
- **Branch end.** A branch ends at the first of:
  - a fold, or a secant target failure (§ THE SOLVER);
  - a Newton failure after the frozen retries (half the step, then a quarter);
  - `ε_max`.

  A change of type does not end the branch: it is recorded and located by bisection to 2% in `ε`, and the record continues.

## THE LINEAR VERDICT (at every ladder point, at level 3)

- **The tracked symmetry cluster `G_tol`.** It is defined by structure, not by size. It holds:
  - the exact fibre-symmetry vectors (`iΦ`, and in `R2` also `JΦ` and `iJΦ`) and their generalized partners;
  - the lifted right-rotation modes: the eigenvectors whose static part has `M`-overlap at least 0.9 with `span{X_aΦ}`, with `X_a` the discrete right-rotation derivatives: the `L²`-projection of the derivative along the left-invariant fields `x ↦ x e_a/2`, which generate the right rotations, assembled with a 5-point rule. A radical direction's Jordan pair splits on the mesh into two such eigenvectors, and each identified eigenvector brings its Hamiltonian partners;
  - in C1, also the rigid family's pair: the eigenvectors whose static part has `M`-overlap at least 0.9 with the interpolants of `π(x)δv`, for `δv` in the two directions of the family's tangent at `v` that the symmetry orbit does not span.

  Membership is read against the union of these structural spans, with the phase partner's static part included. A degenerate cluster, such as C1's six zeros, comes back from the eigensolver in an arbitrary basis.

  Lifted modes are tracked from one ladder point to the next by overlap. `τ₀` bounds how far the identified paired modes may sit from 0, and the radical bound (§ NUMERICAL FLOORS) bounds the radical ones. Any other eigenvalue in `|λ| ≤ τ₀` is a zero-crossing event: the local type is Unresolved until it is resolved, and it is never absorbed into `G_tol`. The lifted modes' real parts are reported in every case.
- **The cluster checks.** If any of these fails, the point is Unresolved. `G_tol` must:
  - have the frozen algebraic dimension `k + r_u`, the slow matrix's zero multiplicity in the Result's tables: 4 for `v₂`, 8 for the prism and the hexagon, 6 for the third outside orbit, 10 for T5; for the controls, 6 for C1 (Theorem D), 4 for C2a, 6 for C2b and 8 for C3;
  - be closed under `λ ↦ −λ, λ̄, −λ̄` and under its Jordan chains;
  - have an Ω-Gram matrix that is nondegenerate, with its smallest singular value at least 10⁻⁶ of its largest, in an `M`-orthonormal basis.
- **The Krein count, on the mesh's kernel.**
  - `B` is the full phase-space second variation of `E − ωN` at the standing wave. Its velocity block is positive, so `n(B) = n(L)`, the stationary Jacobian's negative index, which is how the total is computed. `n(B|G_tol)` is the negative index of `B` restricted to the phase-space subspace `G_tol`, from its Gram matrix there. `G_tol` and its Ω-complement are spectral subspaces, so the difference is the count on the complement.
  - The identity `k_r + 2k_c + 2k_i⁻ = n(B) − n(B|G_tol)` must balance with the computed eigenvalues outside `G_tol`. Here `k_i⁻` counts imaginary pairs of negative Krein signature, and `k_r`, `k_c` count real pairs and quartets beyond `τ_Re`.
  - **Representatives.** Each threshold-classified Hamiltonian pair or quartet counts once, through its member with `Re λ ≥ 0` and `Im λ ≥ 0`. Both signs are read at `τ_Re`, as the classes are: a computed eigenvalue is a representative when `Re λ ≥ −τ_Re` and `Im λ ≥ −τ_Re`. An imaginary pair's computed real part is round-off of either sign, and so is the imaginary part of a degenerate real pair that the eigensolver returns as a conjugate pair.
  - If the slow eigenvalues alone do not balance it, the fast spectrum is searched by shift-invert on a frozen grid: shifts `i jω/20` for `j = 0, …, 40`, plus the same shifted off the axis by `ω/40` to catch quartets. The 24 eigenvalues nearest each shift are kept if their residuals pass, they are representatives, and they lie beyond the slow set in modulus, de-duplicated at 10⁻⁶ relative. If the count still does not balance, the point is Unresolved.
- **Level stability.** A candidate instability, a real part above `τ_Re` outside `G_tol`, is declared only if its real part at levels 2 and 3 agrees within the level-stability tolerance `δ_lev = max(10⁻², 2 d₃)`. Here `d₃` is the largest relative difference between C3's level-2 and level-3 rates over its three `ε`, fixed from C3 before any target runs. Closing: `d₃ = 1.25×10⁻⁴`, so `δ_lev = 10⁻²`. The calibration can fail: if `δ_lev > 5×10⁻²`, this version does not freeze. A discretization artifact shrinks about 16 times per level; a genuine instability does not. A disagreement makes the point Unresolved.
- **Elliptic** if `k_r = k_c = 0` with the count balanced. **Hyperbolic** if a level-stable `k_r + k_c > 0`. **Unresolved** if the count cannot be balanced, at a zero-crossing event, on a level disagreement, or before a mode's pre-registered first-resolved point.
- **Reported:** the slow eigenvalues with Krein signs; Krein collisions and the `ε` where they occur; `n(B)` and `n(B|G_tol)`; the floors, and every lifted mode's modulus and real part.
- **Cross-level check.** It runs at every third ladder point, and at every point where an instability or a change of type is declared.
  - **Matching.** At each level, every Hamiltonian pair or quartet of slow eigenvalues outside `G_tol` is represented as in the Krein count, with both levels' signs and classes read at `τ_Re⁽³⁾`. The representatives are matched one-to-one between the levels within each class, by the minimum total distance `Σ |λ⁽³⁾ − λ⁽²⁾|` (the assignment problem). The classes are imaginary pairs of each Krein sign, real pairs, and quartets. Unequal counts in a class are a disagreement.
  - **Agreement.** Each matched pair satisfies `|λ⁽³⁾ − λ⁽²⁾| ≤ 0.01 max(|λ⁽³⁾|, |λ⁽²⁾|) + s₂/ω₀`. Here `s₂` is the slot's measured level-2 spread of its lowest level (frozen: `6.739×10⁻⁴` in `R4`, `9.547×10⁻⁴` in `R5`, `2.136×10⁻³` in `R2`). So `s₂/ω₀` is twice level 2's splitting in frequency, and level 2's resolution of the smallest slow frequencies is not read as a disagreement.
  - A disagreement makes that point Unresolved.

## THE NONLINEAR PERSISTENCE VERDICT (at level 2)

- **Where.** At two scheduled points, the ladder points nearest `ε_max/8` and `ε_max/2`.
  - A scheduled point not reached as an elliptic point, failing the drift margin (below), or outside the slow subspace's lineage, is replaced by the last accepted elliptic ladder point below it that meets the margin and lies inside the lineage. That covers a branch that turns hyperbolic, folds, fails first, or loses its slow subspace.
  - Two scheduled tests that land on one point run once.
  - If no elliptic ladder point above `ε_min` exists, the test is recorded Not reached. If elliptic points exist but none lies inside the lineage, it is recorded Not run (slow subspace unresolved); if some do but none meets the drift margin, Not run (floor). None of these is Persists or Fails.
- **The perturbation.**
  - **Size.** Each half has energy norm `η/√2`, with `η = 10⁻³` relative to `‖Φ‖_M`, in the energy norm of the rotating-frame phase space, `‖(ζ, ζ')‖²_E = ζᵀMζ + ζ'ᵀMζ'/c²`. The halves are not orthogonal in that norm, so their sum is only about `η`.
  - **The draw.** For each persistence test, a fresh generator `numpy.random.default_rng(20261006)` draws `standard_normal` for the real parts, then the imaginary parts, of the reduced unknowns, in the solver's node-orbit order. Execution order therefore cannot change a perturbation. The draw is taken as a static phase-space vector `v`.
  - **The slow subspace `S`,** tracked at every accepted level-2 ladder point.
    - **The candidate** is unique and has the frozen dimension, 14 in `R4` and `R5` and 16 in `R2`: `G_tol` plus the lowest-`|λ|` non-cluster Hamiltonian modes.
    - **The gap:** `|λ_next| ≥ 2 max |λ|` over `S` minus `G_tol`, where `λ_next` is the next eigenvalue in modulus.
    - **Continuity:** `σ_min(U_prevᵀ M U_S) ≥ 0.9` for the phase-space `M`-orthonormal bases of the previous point's `S` and this one's: every principal cosine. At `ε_min`, the branch's validity check and the cross-level check take its place.
    - **The lineage** ends at the first point that fails a check, with no re-entry. Persistence tests are placed inside the lineage only (below).
    - C1, whose slow spectrum is its symmetry cluster alone, has `S = G_tol`.
  - **The declared split.** `P_G` and `P_S` are the Ω-orthogonal projections onto the level-2 cluster `G_tol` and onto `S ⊃ G_tol`. Both are built at level 2 by the same structural rules, over the full space rather than the pinned fixed space.
    - **The slow half** is `(P_S − P_G)v`, normalized to `η/√2`: the genuine slow modes, where the question lives.
    - **The fast half** is `(1 − P_S)v`, normalized to `η/√2`. It exercises fast–slow interaction and keeps the perturbation's energy positive.
    - **With no genuine slow mode,** as in C1, whose slow spectrum is its symmetry cluster alone, the whole `η` goes to the fast half.
    - The lifted rotation modes are not excited at linear order. An unsplit white draw puts only about 5% of `η` into the slow modes at level 2 (3% at C3, 8% at C2a), so the test would probe mostly the fast modes, which are coercive and stable.
  - Position and velocity are perturbed together: `(Φ + δ, Π₀(Φ + δ) + p)`, with `(δ, p)` the projected direction's position and rotating-frame velocity parts, and `Π₀` the standing formula. This is one restricted phase-space direction, and the verdict is stated for it, not as generic phase-space stability.
- **Duration** `T = 1000`, with `dt` 0.8 times the level-2 stability limit `2/ω_max`: 0.008197 in `R2`, 0.008158 in `R3`, 0.008187 in `R4` and `R5`. The margin `r_f⁽²⁾ T` is 0.040.
- **The measure.** The perturbed and reference trajectories are integrated in lockstep, and compared every unit of time, at the nearest step. Their difference is aligned modulo the exact continuous symmetries (the phase in `R4` and `R5`, Sp(1) in `R2`), in closed form, and `d(t)` is its `M`-norm. Two parts of it are read with the level-2 `G_tol`. All three are relative to `‖Φ‖_M`, as `η` is:
  - `d_⊥`, the `M`-norm of its component off `G_tol`, from the Ω-orthogonal spectral projection;
  - `d_rot`, the `M`-norm of its `M`-orthogonal projection onto `span{X_aΦ}`: how far the run has turned along the lifted rotation directions.
- **Drift along the lifted modes.** At level 2 a radical lifted mode can be a real pair. On unpinned C3 branches the largest rate runs from `8×10⁻⁴` at `ε = 0.5` to `1.0×10⁻²` at `ε = 16`; at C3's pinned points the radical modes are imaginary. The perturbation is projected off it, but nonlinear terms seed it at about `η²`. So the drift reaches about `η² e^{rT}`, with `r` the largest real part in the branch's level-2 cluster at that point.
  - **The margin.** A test runs only where `r T ≤ ln(1/(10η)) ≈ 4.6`, which keeps that drift below `η/10`. A point that fails it falls under the substitution rule above.
  - This rule can move or withhold a test, and the verdict below can invalidate one. Neither ever scores one.
  - **Expected coverage.** The drift margin applies to each target's own pinned lifted modes.
    - At C3's pinned points they are imaginary across the ladder, so the margin is not expected to bind there.
    - T4, pinned only by `C₂`, is where real pairs may still appear.
    - The linear verdict at level 3 covers the full ladder.
- **The verdict reads `d_⊥`.** The lifted modes are left out of the Persists and Fails thresholds, so their drift can never manufacture a Fails. Excessive rotation drift instead invalidates the test, under the first-event rule:
  - **Invalid (drift)** if `d_rot` exceeds `√η` before `d_⊥` exceeds `100η`. Up to that turn, the leak into `d_⊥`, second order in the turn, stays within about `η/2`.
  - **Fails** if `d_⊥` exceeds `100η` first, whatever the later drift. A tie at one observation is Invalid (drift).
  - Otherwise, over `[0, T]`: **Persists** if `max d_⊥ ≤ 10η`, and **Unresolved** if not.

  `d`, `d_⊥`, `d_rot`, and the reference trajectory's own distance from the family are all reported. The last is the floor of the measure.
- **Cross-level check.** The same test at level 1 runs only where level 1's floor rate times `T` is at most 0.1, and there it must give the same verdict. A disagreement makes the point Unresolved. Where the margin fails, the check is not run, and that is recorded.

## THE VERDICT RULE, AND STOP CONDITIONS

- **Controls first.** A failed control gate ends that frozen version of the run.
  - A repair is made only while no target has run. The solver and the terms then get new hashes on `main`, with a dated note of what changed, and the controls are rerun under the new frozen version.
  - No target result from a failed version is retained.
  - Target output produced before all controls have passed is quarantined, and cannot be scored under this pre-registration.
- **Validity.** A point that breaks a conservation floor, a residual bound or the seed rule is invalid and reported as such. It is never scored as stable or unstable.
- **Validity at `ε_min`.** Each branch's slow spectrum at `ε_min`, as `ε κ τ` with its slot's own `κ`, must match its frozen small-amplitude spectrum within 5%, with the types and Krein signs exactly, in every mode resolved there (§ NUMERICAL FLOORS). A mismatch makes the branch invalid, an instrument defect, never a counterexample to the theorems. A mode unresolved at `ε_min` is reported at every ladder point. Its type and Krein sign are recorded from its first resolved point, with no tolerance on its value.
- **Per branch.** The record holds:
  - the last `ε` reached and why the branch ended;
  - the linear verdict at each ladder point, with every change of type located;
  - the persistence verdicts.
- **The answer to the question** is the per-branch record. No branch's outcome is aggregated into a verdict about another.
- **No search.** Every computed number is reported, and nothing is tuned toward a target. Any change to these terms after the go is a new pre-registration.

## FROZEN VALUES

- **Small-amplitude slow spectra,** for unit `u` in `R4`'s `h = Q|w|⁴` normalization (the Result, §IV).
  - `R5`'s are these times 9/16.
  - The physical slow eigenvalues are `ε κ` times these, with `κ = 1/(8 ω₀ Q(u))` and each slot's own `Q`.

  | orbit | `Q` in `R4` | `Q` in `R5` | slow spectrum (`R4`) |
  | --- | --- | --- | --- |
  | `v₂` | 147/143 | 581/572 | `τ = 56/1287` (Krein −); `112/143, 112/117, 224/143, 560/429` (Krein +) |
  | prism | 5831/5031 | 609/559 | `τ² = 5770240/7913763`; `19066880/278420571` twice; all Krein + |
  | hexagon | 1750/1287 | 2751/2288 | `τ² = 90160/61347`; `501760/552123` twice; all Krein − |
  | third outside orbit | 1.185203461338846… | 1.104176947003101… | `τ ≈ 1.01209160489` (Krein +); `0.310692742362, 0.139251832964, 0.00573685922309` (Krein −); all simple. The roots are certified (the Result); values and Krein signs are computed from the certified point at 60 digits by the frozen check (§ TO BE FIXED AT THE FREEZE) |
  | T5, `R2`'s line point | 22/13 (`R2`) | n/a | `τ² = 153664/20449`, and the roots of `418161601ν² + 6989958976ν + 28677390336`; all Krein − |

- **Controls:** `Q = 1288/1287` (C2a), `144/143` (C2b), `175/143` (C3); their spectra are in the controls table.

## FIXED AT THE FREEZE (from the closing control runs)

- **Floors.**
  - `r_f⁽³⁾ = 1.040×10⁻⁶`, so `τ₀⁽³⁾ = τ_Re⁽³⁾ = 1.040×10⁻⁵`.
  - `r_f⁽²⁾ = 3.970×10⁻⁵`, so `τ₀⁽²⁾ = τ_Re⁽²⁾ = 3.970×10⁻⁴`.
- **The level-stability tolerance:** `δ_lev = 10⁻²` (`d₃ = 1.25×10⁻⁴`).
- **`s₂`:** 6.739×10⁻⁴ (`R4`), 9.547×10⁻⁴ (`R5`), 2.136×10⁻³ (`R2`).
- **Persistence:** `η = η_ref = 10⁻³`, `T = 1000`, `dt` as above.
- **The unresolved modes:** only the third outside orbit's smallest, first resolved at `ε = 1.41` in `R4` and 2 in `R5`.
- **Levels:** 3 and 2 stand. The controls and floors pass.
- **The freeze block, in full:**
  - the solver commit and the SHA-256 of every script (the manifest);
  - the terms' SHA-256;
  - the random seed `20261006`;
  - the orientation matrices, pinned subgroups and fibre maps of every case at both levels, with each seed's `Q` against the Result's;
  - the numbers above;
  - the check that prints the third outside orbit's frozen frequencies and Krein signs from the certified point in `branches_check.py` at `d0de8ca`, with its output;
  - the continuation algorithm as written above.

  All of it lands before the run. The solver commit and a terms-hash line on the MIT parent page are provenance; the pre-registration itself is OpenWave's merge of this task and its BACKLOG row, and no target is solved before it.

## DEFINITION OF DONE

- **The record:**
  - every branch's ladder, with its linear verdicts, Krein data and change-of-type locations;
  - the persistence tests;
  - the controls' gates;
  - the floors.

  It is published with scripts, data and a reproduction route.
- **Adjudication** by the maintainer, under the rule above.
- **Records:**
  - a method note, equations first, with an equation-to-code map;
  - the roadmap Done row;
  - after the verdict, the author records a dated note on the Result page.
- **What does not move.** No `MODELS.md` cell, and M8.7's gate is unchanged.
