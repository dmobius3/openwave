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

**2026-08-23, stage 6 retargeted again, and the 2026-08-18 nomination corrected.** The entry
above named `M4_int` across the three connections as the slot survey. That was too quick, and
this task's own calculation is what corrects it: a field in `Ω⁰(X; E_{τ_σ})` has a section at
level `n` exactly when the TRIVIAL isotype occurs in `V_n ⊗ τ_σ`, which is the certification
table's `R_0` ROW read across the three columns. The eight nontrivial rows are the ambient
decomposition of `V_n ⊗ τ_σ`, a fact about the cover, and not eight physical sectors of an
`E_{τ_σ}`-valued field. Checked against the frozen table: the computed first section levels are
0, 2 and 6 for the trivial, standard and Galois bundles, matching that row exactly. M8.2 § 3
already says the gate "certifies the coefficient bundle only" and does not furnish a fluctuation
spectrum.

So **`M4_int` receives a STRUCTURAL N/A for the eight-slot physical-section question**, and the
target-bearing object is `M4L_Erho`, one field per flat bundle `E_ρ` with `E_R0` the mandatory
null control, filed as [`m8_4_preregistration.md`](../findings/m8_4_preregistration.md). Its
ceiling is explicit and narrower than OQ1: the eight sectors are installed by construction, so it
tests whether ONE nonlinear law acts coherently across them without per-sector tuning, never
whether the dynamics selects them. The per-connection contract above does not carry over, because
the connection axis does not survive the substitution: `σ` entered only through `κ_σ`, and in
`E_ρ` the coefficient system is `ρ` itself. A genuine selection experiment needs a nonlinearity
with a representation-changing term and is a separate object and preregistration.

**2026-08-24, P1A closed: estimator qualified, substrate not.** The engineering pilot's first
phase ran and is recorded in [`m8_4_p1a_closeout.md`](../findings/m8_4_p1a_closeout.md). The
invariant-subspace estimator the frozen § 7 observable needs is QUALIFIED: projector identities at
1e-16, Schur against an independently preregistered Riesz contour agreeing to 4e-08, the three-way
leakage identity closing to 8e-16, and spectral projector norms of 2.2 to 3.1 on the target
clusters, so the operator's strong non-normality does not reach the subspaces M8.4 scores.

What FAILED is the substrate's spectral-contamination qualification, globally, under a power gate
frozen before the fit. The continuum operator is self-adjoint, so imaginary parts in the computed
spectrum are error; they run from 5.95e-14 on the `R_0` control to 1.8e-01 on `R_2`, and the envelope needed to
cover that spread is about 58,688 times its own trend, which is a rule with no discriminating
power. No sector verdicts issued, P1A.5 was never opened, and **no nonlinear target configuration
was ever run: no target sector has been spent.**

A precision-ladder localization test, mutation-armed and reproducible from
[`m8_4_p1a/localization/`](../m8_4_p1a/localization/), establishes that the failure is
substantive rather than an instrument defect: with the assembled `L_h` held byte-for-byte fixed,
the target imaginary parts are unchanged from float64 through 50-digit arithmetic. They belong to
the matrix, not the eigensolver, so no numerical floor can reclassify them. The test separates
solver from matrix; it does NOT identify whether the equivariant assembly or the underlying
RBF-FD discretization is responsible, and nothing here claims to.

**The nonlinear pilot is therefore blocked on the present discrete operator.** Continuing would
require separately commissioned, independently justified structure-preserving assembly or
discretization work, which is not a repair to P1A and does not reopen its outcome.

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
