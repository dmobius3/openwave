# M8.4: The Lagrangian-family survey on S³/2I

> Roadmap row: [`../m8_roadmap.md`](../m8_roadmap.md). Status: 🚧 PLANNED, gated by
> M8.2 + M8.5. This is a scaffold-stage planning aid written by the maintainers
> (2026-07-21); the author owns the column and may amend everything here.

## PLANNING

### Scope (the decisive science of the M8 program)

The central question ([`../m8_theory_canonical.md`](../m8_theory_canonical.md) OQ1):
can a nonlinear field equation on S³/2I have topological-defect or standing-wave
solutions whose energies realize the McKay SLOT STRUCTURE, without per-slot tuning?
Candidates come from the platform's working families
([`../m8_platform_pointers.md § 2`](../m8_platform_pointers.md)): M5 Landau-de Gennes
matrix field + Frank energy, M4 nonlinear vector wave, M7 two-vector; each written on
the compact quotient with the background clock. MIT supplies the target structure;
this task supplies (or refutes) the missing dynamics. Either outcome is a result.

### Suggested per-family pipeline (order matters)

| # | Stage | Why this order |
| --- | --- | --- |
| 1 | Write the family on the arena (covariant derivatives on S³/2I; the anti-periodic double-cover sector included) | the arena changes the operator, not just the domain |
| 2 | VACUUM SPECTRUM FIRST: linearize about the vacuum and compute the band structure on the arena | the M7 lesson: a truncation's vacuum can be tachyonic; hunting solitons on an unstable vacuum wastes everything downstream |
| 3 | Derrick / scaling analysis ON THE COMPACT ARENA, explicitly | flat-space Derrick conclusions do not transfer: R provides a scale; the background clock adds the oscillation escape (the M6-validated third route) |
| 4 | Defect-sector census (OQ5): which homotopy sectors exist for this family's target space on the quotient; does anti-periodicity create new ones | knowing what CAN exist bounds the search before compute is spent |
| 5 | Relax / evolve candidate states per sector; measure energies with window-robust observables | the M6 lesson: window-defined observables silently manufacture spectra |
| 6 | Compare against the M8.2 lock; verdict per family | proportionality without per-slot tuning, or refuted |

### Suggested definition of done

| # | Item |
| --- | --- |
| 1 | Per-family verdict table (vacuum stable? sectors exist? localized states found? energies vs slots?) with every number reported |
| 2 | Adversarial audit of any positive claim (independent script, different method) BEFORE it is trusted |
| 3 | Method note per [`dev_docs/METHOD_NOTE.md`](../../../../../dev_docs/METHOD_NOTE.md); MODELS.md cells flipped honestly (a clean negative across all families is a publishable close, not a failure) |

### Blindspots

| Risk | Guard |
| --- | --- |
| Per-slot tuning creeping in via "reasonable" coupling adjustments | the M8.2 lock enumerates couplings + bounds in advance; every point run is reported |
| Discretization-selected minimizers (the M5 stencil lesson) | cross-check stationary states on a second discretization before claiming them |
| Energy windows (the M6 lesson) | report energies vs window size; a window-drifting energy is not a result |
| Compute sprawl | the sector census (stage 4) prunes before the expensive stage 5 |

### Ownership + gating

Author-driven with platform support (this is the collaboration's core). Gated by M8.2
(the lock) and M8.5 (the engine).

## DEVIATIONS LOG

**2026-08-18, stage 6 retargeted; stages 1-5 decoupled.** As scaffolded, stage 6 presupposed
that quotient states carry deck-isotypic slot content. The task's first result (FINDINGS 1)
proves they cannot, for every native family and every transport convention, so the slot
comparison retargets onto the M8-owned twisted object `M4_int`
([`m8_2_preregistration.md § 6.1`](../findings/m8_2_preregistration.md)) across the three
frozen flat connections, with `σ_0` the mandatory null control and a per-connection contract:
each connection gets its own complete record and verdict, no post-hoc connection selection,
no pooling that lets one connection rescue another, and any cross-connection claim
pre-declared before the first target-bearing run. Stages 1-5 decouple from OQ1 and survive as
an optional descriptive native-family survey (OQ5 census, existence, stability, Derrick on
the compact arena). The frozen free per-connection tables become calibration gates earning
zero evidentiary credit; target-bearing observables must be genuinely nonlinear. `M5 + P`
and `M7_ad` stay outside the coming preregistration, with no contingency clause admitting
them on an `M4_int` outcome.

## FINDINGS

1. **THE KINEMATIC CLOSE (2026-08-18).** Native single-valued fields on S³/2I carry none of
   the 8 nontrivial McKay slots at any harmonic level: every quotient field lifts
   2I-invariant, and the isotypic projectors annihilate invariant content. The spinorial
   half, the electron slot included, is absent by center parity alone; the control realizes
   the full `n = d` ladder in the twisted sectors. OQ1's native branch closes negative by
   theorem, with no dynamics run. Record:
   [`../findings/m8_4_kinematic_close.md`](../findings/m8_4_kinematic_close.md) +
   [`../scripts/m8_4_kinematic_check.py`](../scripts/m8_4_kinematic_check.py) (two
   computationally separate routes, mutation-armed, exit 0).
