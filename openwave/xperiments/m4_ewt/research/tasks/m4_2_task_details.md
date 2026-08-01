# M4.2: Emergent Coulomb force (far-field direction, sign, 1/r)

> Task **M4.2** (M4 / EWT model). Status: **BACKLOG**, no owner · Roadmap:
> [`../m4_roadmap.md`](../m4_roadmap.md) · Migrated from GitHub issue
> [#202](https://github.com/openwave-labs/openwave/issues/202) on 2026-08-01, when
> [T5](../../../../../dev_docs/tasks/t5_task_details.md) settled that tasks live in roadmaps and
> issues are reserved for platform defects. The issue body is archived below in full, so closing
> it loses nothing.

This doc is the task's full record: planning, then findings when someone runs it.

## PLANNING

### The target

Numerical evidence that the electric force **emerges from wave interference**: two wave centers
give a charge-dependent force direction (same repels, opposite attracts) with `~1/r` scaling,
through wave physics alone, with the charge sign **not** imposed as a `±1` label. This is the ❌
honest negative on the Coulomb row of [`MODELS.md`](../../../../../MODELS.md) for this lineage:
sinc envelope barriers block far-field attraction and repulsion, and the signed envelope is a
modeling choice.

### Why it is the pivotal problem

| Reason | Detail |
| --- | --- |
| First force in the hierarchy | Coulomb is the first far-field force to emerge (longitudinal `∇E_L`); magnetic, strong and orbital forces build on it. No emergent Coulomb means no validated far-field force at all |
| The core emergence claim | The program set out to show charge and Coulomb come from waves rather than a hand-placed sign. If the sign must be imposed, the forces-from-waves thesis stays unproven for EWT |
| Coupled to K-selectivity | Same suspected missing physics (variable `λ(r)`, spin `L → T`, non-linearity) as [M4.1](m4_1_task_details.md) |

### The precise blocker

The interaction energy of two coherent monochromatic spherical wave centers goes as
`cos(k·Δr + Δφ)`, so the sinc oscillation flips the absolute force direction every `λ/2` of
separation. It is proven math rather than a numerical artifact: any model with a single frequency,
spherical `1/r` propagation and coherent superposition produces it. Two sub-problems follow: the
far-field direction never settles, and the isolated cross-term gives same → attract /
opposite → repel, inverted from Coulomb (and possibly the strong-force capture mechanism instead).

Ten approaches were already eliminated in the M3 Phase 1 sweep (all five wave equations, Gaussian
smoothing, signed disturbance, the wave-center disturbance models, spin `L → T`, variable `λ(r)`,
2D/3D gradient integration, statistical averaging); the eliminations and the surviving leads are
tabulated in the archive below and in
[`0_STATUS.md`](../../../m3_wolff_lafreniere/research/0_STATUS.md) and
[`1z_results.md`](../../../m3_wolff_lafreniere/research/1z_results.md).

### The two best first contributions

| Route | Why |
| --- | --- |
| 2D flux → 3D spherical flux | Radiation pressure `S = -c²·ψ·∇ψ` already gives 100% charge discrimination and `~1/r`, Coulomb-correct at half-integer `λ`; full solid-angle averaging on a sphere may smooth the residual `λ/2` oscillation. Untested |
| The M4 vector route | A scalar substrate collapses directional information, so an emergent charge SIGN plausibly needs the vector field: divergence, curl, Poynting-like flux, elliptical-rotation handedness. The substrate exists ([`../M4_engine_upgrade.md`](../M4_engine_upgrade.md) P0-P4) |

### What would resolve it

A wave-physics computation that, at **all** far-field separations, yields a consistent Coulomb
direction (same = repel, opposite = attract) with `~1/r` scaling and an emergent sign, while
leaving the validated near-field results (lock-in, annihilation) intact. A demonstration that no
such computation exists within the model is equally a result and closes the task.

**Gated by**: an owner. The M4 column has no active author, and the program belongs to whoever
extends the model.

## GitHub issue archive (#202)

> Migrated from OpenWave GitHub issue
> [#202](https://github.com/openwave-labs/openwave/issues/202) on 2026-08-01. Title: "Solve
> emergent Coulomb force in EWT: far-field direction/sign + emergent charge". Opened 2026-06-15 by
> `xrodz`. State at migration: OPEN. Labels: `help wanted`. Body verbatim, with its headings
> demoted one level so they nest under this section.

### Goal

Produce numerical evidence that the **electric (Coulomb) force emerges from wave interference** in
the EWT / Wolff-LaFreniere model: two wave centers (WCs) yield a charge-dependent force direction
(same repels, opposite attracts) with `~1/r` scaling, through wave physics alone, **without
imposing charge as a `±1` label**.

This is currently the **❌ honest-negative** on the Coulomb row of
[`MODELS.md`](https://github.com/openwave-labs/openwave/blob/main/MODELS.md) for M3: "sinc envelope
barriers block far-field attraction/repulsion; signed envelope is a modeling choice."

### Why this is the pivotal M3 problem

| Reason | Detail |
| --- | --- |
| It is the **first** force in the EWT hierarchy | Coulomb / electric is the 1st measurable force to emerge (longitudinal `∇E_L`), per [`0_OVERVIEW.md`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m3_wolff_lafreniere/research/0_OVERVIEW.md). Magnetic (`∇E_T`), strong (sub-λ electric+magnetic), and orbital forces all build on it. No emergent Coulomb → no validated far-field force at all |
| It is the **core emergence claim** | The whole M3 program "set out to" show charge + Coulomb come from waves, not from a hand-placed sign ([`1z_results.md §1`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m3_wolff_lafreniere/research/1z_results.md)). If charge sign must be imposed (`cos(source_offset)`), the forces-from-waves thesis is unproven for EWT |
| It **gates the model's credibility vs M5/M6** | On `MODELS.md` M5 already has Coulomb ✅ (1/r from pure topology); M3 is ❌. Closing this is what moves EWT from "near-field mechanism only" to "produces real forces" |
| It is **coupled to the K-selectivity problem** | Shares the same suspected missing physics (variable `λ(r)`, spin `L→T`, non-linearity). Companion issue: #201 |

> Honest framing (from `0_STATUS.md`): M3 scalar with monochromatic waves demonstrates the
> near-field MECHANISM (lock-in) but not the far-field FORCE (emergent Coulomb) or the SELECTIVITY
> (#201). Both likely need the missing physics.

### The precise blocker: the sinc oscillation

The interaction energy of two coherent monochromatic spherical WCs is `∝ cos(k·Δr + Δφ)`. This
**sinc oscillation flips the absolute force direction every `λ/2`** of separation. It is proven
math, not a numerical artifact: any model with (1) a single frequency, (2) spherical `1/r`
propagation, and (3) coherent superposition produces it.

| The sinc DOES produce ✅ | The sinc does NOT produce ❌ |
| --- | --- |
| Near-field lock-in (same-phase wells at `λ/2`) = strong-force binding | Consistent far-field Coulomb **direction** at all separations |
| Annihilation (opposite-phase well at `Δr = 0`) | Correct Coulomb **sign** in the energy-gradient approach |
| Relative charge discrimination (same ≠ opposite) | Emergent charge sign (currently imposed via `cos(source_offset)`) |

Two specific unsolved sub-problems: (a) the **far-field direction oscillation** (`cos(k·Δr)` never
settles to "same always REP, opposite always ATT"); (b) the **sign problem**, the isolated
interaction cross-term gives same→ATT / opposite→REP, which is 180° inverted from Coulomb (and may
instead be the strong-force capture mechanism).

### Already systematically eliminated (Phase 1)

| Approach | Why it failed |
| --- | --- |
| All 5 wave equations (Wolff, LaFreniere-Marcotte, phase-warped, Combined W-L, weighted PSW) | Oscillation intrinsic to coherent interference |
| Gaussian smoothing | Removes the oscillation AND the charge info |
| Signed disturbance (`1a`) | Charge imposed, not emergent |
| 10 WC disturbance models (`1b`) | Only `L→T` distinguishes charges; passive models fail in isotropic fields |
| Spin `L→T` conversion (`1c`) | Creates **magnetic** (transverse), not electric (radial) |
| Variable `λ(r)` energy (`1d`) | Charge-blind: `λ` depends on K (structure), not phase (charge) |
| 2D/3D gradient integration | Base wave dominates WC-WC term |
| Statistical averaging, broadband shells, in/out-wave decomposition | Sinc symmetric / beating / present in all components |

### Most promising open leads (carry-over to Phase 2)

| Lead | Status |
| --- | --- |
| **Flux / radiation pressure** `S = -c²·ψ·∇ψ` (not energy gradient) | Closest result so far: 100% charge discrimination, `~1/r`, Coulomb-correct at half-integer `λ`; remaining issue = the `λ/2` absolute oscillation. LaFreniere uses radiation pressure, not `F = -∇E` |
| **3D spherical flux integration** | The 2D test used a circle; full solid-angle averaging on a sphere may smooth the residual oscillation (untested) |
| **LaFreniere `λ/2` core phase shift** | The compressed core (1λ, 7× smaller volume → shorter inner `λ`) may supply the missing `π` offset that flips the sign |
| **Variable `λ(r)` profiles** | Yee & Hauger shells (`λ` longest near core) vs LaFreniere (shortest inside core) is an unresolved contradiction; test both |
| **Non-linear `ψ³`** (Smoliński NLS soliton stabilizer) | Changes the spatial structure away from pure sinc; ties to #201 |
| **K-scale + vector field** | Coulomb requires a K≥10 standalone particle whose spin drives `L→T` → far-field traveling wave; a single WC (K=1 neutrino) is correctly neutral. Scalar M3 collapses directional info, so emergent charge **direction** likely needs the M4 vector field (divergence / curl / Poynting-like flux, elliptical-rotation handedness as the charge-sign indicator) |

### What would resolve this

A wave-physics computation that, at **all far-field separations**, yields a consistent Coulomb
direction (same = REP, opposite = ATT) with `~1/r` scaling, where the **charge sign is emergent**
(from geometry / spin / phase, not an imposed `±1` or `source_offset` convention). Bonus: the same
model should leave the validated near-field results (lock-in, annihilation) intact.

### References

- [`0_STATUS.md`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m3_wolff_lafreniere/research/0_STATUS.md)
  (open blockers + 7 K-selectivity avenues)
- [`1z_results.md`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m3_wolff_lafreniere/research/1z_results.md)
  (Phase 1 Coulomb results, eliminated mechanisms, the sign problem, the 2D-flux lead)
- [`0_OVERVIEW.md`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m3_wolff_lafreniere/research/0_OVERVIEW.md)
  (far-field force hierarchy)
- [`0b_additional.md`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m3_wolff_lafreniere/research/0b_additional.md)
  (two-domain Energy/EMC force separation, Smoliński non-linear terms)
- Companion open problem: #201 (K-selectivity & perturbation robustness)

Help welcome. The 2D-flux → 3D-spherical-flux path and the M4 vector route are the two best first
contributions.

## DEVIATIONS LOG

(none)

## FINDINGS

(pending: the task has not been run)
