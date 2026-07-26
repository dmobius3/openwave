# M8 background: what MIT has, what it lacks, and why the combination is the program

> **Purpose.** This is the gap map behind the M8 column: an honest statement of what
> Mode Identity Theory brings, what it is missing, and why a field-dynamics program on
> its arena is worth running on this platform. It was written by the maintainers at
> scaffold time (2026-07-21) from the author's onboarding submission
> ([discussion #312](https://github.com/openwave-labs/openwave/discussions/312)), the
> author's repo artifacts (claim ledger, null tests, calculator), and the platform's
> own onboarding evaluation. The author owns the science; corrections via PR or
> discussion are welcome.

## 1. What MIT is

MIT is a **top-down structural model**: a fixed arena (the three-sphere S³, its binary
icosahedral quotient S³/2I, and a Möbius edge), a single background standing wave
`Ψ = cos(t/2)` sampled on it, and a phase-to-time map (the Waltz clock
`dt/dτ = S^(-1/2)`). Particles are representation-theory slots: an irrep of 2I at a
McKay-lattice position, read in one of the three flat connections (the three
generations), with color from the Z₃ face stabilizers. The couplings, the fermion mass
ratios, and Λ are read off this structure (McKay distance, Reidemeister torsion, the
120-cell), not evolved from a dynamics.

The existing OpenWave columns run the same program from the opposite end: they start
from a nonlinear field Lagrangian and evolve it until particles emerge, and their
hardest open item is exactly where MIT is strongest, the **origin of the numbers**
(M5's lepton row names its open piece as "the eigenvalue-hierarchy origin,
1 : 5.9 : 15.1"). MIT is strong where the emergent models stay honestly open, and
absent where they are strong.

## 2. The gap map (what the M8 program must supply)

| Gap | The honest statement | What supplies it |
| --- | --- | --- |
| No field equation | MIT has no Lagrangian and no equation of motion; Einstein's equations are kept unchanged as the geometry's dynamics. Masses are assigned by structure, not evolved | M8.4: the Lagrangian-family survey, drawing candidates from M4/M5/M7 ([`m8_platform_pointers.md`](m8_platform_pointers.md) § 2) |
| No defect concept | MIT's topology is the ARENA (the manifold itself), not topology IN a field (M5's sense: defects, windings, vortex strings). Its particles are explicitly not vortices | M8.4: put a dynamical field ON S³/2I and ask whether its defect or standing-wave sectors realize the McKay slots |
| Stability asserted, not derived | "The first positive level is stable across the cone's self-adjoint extensions; matter modes return under the double cover" is a spectral statement, not a dynamical one. No Derrick analysis exists because there is no soliton | M8.4: Derrick analysis on a compact arena (see § 4 below) |
| Clock assumed, not derived | The Waltz exponent −1/2 is empirically forced (Δχ² > 60 vs integer alternatives) but not derived from the embedding. The opposite of M5, where the clock is the measured energy-minimizing state | The M8 dynamics question: does the background mode + field coupling SELECT the clock |
| Radiation unmodeled | photon/gluon are "massless at the edge-only level", a layer label with no wave mechanics behind it | downstream of a validated field dynamics; not an early target |

## 3. Evidence weight (the author's own ledger, adopted as the column's grading)

The author's `claim-ledger.md` runs the platform's
[`ONBOARDING_MODELS.md`](../../../../ONBOARDING_MODELS.md) § 4 parameter count on the
framework itself, and the M8 column adopts its verdicts as the grading baseline:

| Layer | Weight | Why |
| --- | --- | --- |
| Structural results | the durable core | 3 generations from the 3 flat connections of S³/2I, the coexact mass gap 4/R², color from Z₃ stabilizers, chirality: retrodictions with zero adjustable freedom (the ledger's "bones") |
| The 2/R² eigenvalue | decidable, unverified | analytic result awaiting independent recomputation; adopted as the M8.1 certification gate |
| The 24-slot mass table | LOW evidential weight | the author's own pre-registered null: random torsion reassignment reproduces or beats the within-3× coverage at p = 0.174; the well set {13, 21, 34, 55} is never extremal under its own variational test (12.2 percentile); 8 slots unassigned; down 3.2×, top 3.9× residuals listed |
| The calibration web | consistency, not prediction | α is input-and-output (Cycle 2); m_e ↔ Λ closes to ~11%; the two routes to R disagree ~4× (Cycle 7, the framework's flagged master-scale tension) |
| Forward bets | the falsifiable layer | Euclid DR1 (mid 2027) pre-registered thresholds; ν₂ = 8.6 meV (JUNO/DUNE) |

The practical rule for M8 tasks: **target the structural ladder (slot structure, gap
ratios, generation count), not the numeric mass table.** A dynamics whose spectrum
reproduces the McKay SLOT STRUCTURE without per-slot tuning would be a result; chasing
the 24 numbers would inherit freedom the author's own audit already flagged.

## 4. Why the combination is technically attractive

| Feature | Why it matters |
| --- | --- |
| Compact arena vs Derrick | On a compact manifold the Derrick scaling argument loses its free dilation (the radius R sets a scale), so the collapse pressure that constrains flat-space solitons is structurally weakened. And the platform's measured lesson (M6: oscillation is a genuine third Derrick escape; M5: stability requires time-periodicity) rides naturally on a background clock: the two mechanisms stack |
| The Möbius anti-periodicity | The twisted boundary condition (return only under the double cover) is a geometric home for 720-degree return, the spin-statistics row that most columns leave 🚧. A dynamical field on this arena inherits it by construction |
| Spectral basis = target ladder | Simulating on S³/2I naturally expands fields in 2I-symmetric harmonics, the same representation theory the McKay ladder lives in. The simulation basis and the target spectrum speak one language, so "do defect energies land on McKay slots" is a well-posed, pre-registrable question |
| Clean division of labor | MIT supplies the arena and the target spectrum (top-down); the platform supplies Lagrangian candidates, simulation engineering, and grading standards (bottom-up). Either outcome is a result: energies on the slots validates both halves; energies off the slots is a clean documented negative, which this author demonstrably accepts |
| Cross-model payoff regardless | M8.6 (the McKay rule vs M5's measured 1 : 5.9 : 15.1 hierarchy) is a bounded test of the ONE candidate mechanism currently on the table for M5's open lepton-hierarchy origin. It needs no simulation and pays off whichever way it lands |

## 5. Onboarding evaluation of record (2026-07-21)

| Check | Result |
| --- | --- |
| Provenance | author account + repo real and active; all artifacts named in the submission exist (`calculator.html`, `claim-ledger.md`, `mass-null-test.py` + frozen inputs/results, per-paper test scripts); all 10 Zenodo DOIs machine-verified resolving (SSRN IDs author-attested, see [`../theory/_CITATIONS.md`](../theory/_CITATIONS.md)) |
| [`ONBOARDING_MODELS.md`](../../../../ONBOARDING_MODELS.md) § 1 fit | partial criteria coverage (particle rows yes, dynamics rows structurally absent); genuine forward predictions exist; reproducible (calculator + scripts + public data); strongly falsifiable (dated pre-registrations with thresholds, 3 documented negatives on record) |
| § 4 parameter count | run by the author in advance; freedom found non-trivial and DISCLOSED (the nulls and cycles above), which is the honest shape the test exists to find |
| § 5.1 red flags | one real hit: unrefereed sourcing (all venues preprint-tier), recorded in the citations registry; the other flags largely pass BECAUSE residuals and negatives are preserved rather than smoothed |
| Rigor culture | pre-registration with frozen tags, one-shot runs, negatives reported unprompted: compatible with [`AI_HYGIENE.md`](../../../../AI_HYGIENE.md) and the platform's method-note standard from day one |

Scaffold decision: admit as a research-mode column at 21 🚧, with the certification gate
(M8.1) maintainer-run before any deeper commitment, and the author driving the column
from there ([`m8_roadmap.md`](m8_roadmap.md)).
