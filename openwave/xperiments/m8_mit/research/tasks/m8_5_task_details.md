# M8.5: Quotient-manifold simulation engineering

> Roadmap row: [`../m8_roadmap.md`](../m8_roadmap.md). Status: 🚧 PLANNED, gated by
> M8.2. This is a scaffold-stage planning aid written by the maintainers (2026-07-21);
> the author owns the column and may amend everything here.

## PLANNING

### Scope

Build and certify the simulation infrastructure M8.4 needs: fields evolving on the
compact quotient S³/2I. No existing OpenWave column runs a curved compact arena, so
this is genuinely new platform ground. Two candidate routes
([`../m8_platform_pointers.md § 6`](../m8_platform_pointers.md)); prototype BOTH far
enough to choose one on evidence.

| Route | Sketch | Known risks |
| --- | --- | --- |
| (a) 2I-equivariant grid | an S³ grid (embedding or intrinsic charts) with the 120-element identification imposed as an equivariance/ghost-cell map | the identification map bookkeeping; chart seams; where the Möbius edge / cone structure of the MIT arena lives on the grid |
| (b) Spectral in 2I-symmetric harmonics | expand fields in S³ harmonics restricted to 2I-invariant (or covariant) subspaces; evolve coefficients | nonlinear terms need convolution handling (cost grows fast with band limit); but the basis IS the McKay representation theory, so slot structure is manifest |

### The certification benchmark (fix before building)

Certify each prototype on a problem with a KNOWN answer, not on the target problem:
the free Laplacian on S³ has eigenvalues `l(l+2)/R²` with known multiplicities, and on
S³/2I the multiplicities restrict by 2I-invariance (computable independently by
character theory). A prototype that reproduces that spectrum + multiplicity pattern is
certified; one that cannot is refuted before any physics rides on it. This mirrors the
M8.1 gate philosophy one level up.

### Suggested definition of done

| # | Item |
| --- | --- |
| 1 | Both prototypes pass the certification benchmark (spectrum + multiplicities), scripts + JSON in the repo |
| 2 | Trade-off table measured, not argued: accuracy vs cost vs implementation complexity at matched resolution |
| 3 | Route decision recorded with its rationale; the losing prototype kept as the cross-check tool for M8.4 |
| 4 | Prototypes are research scripts (NumPy/SciPy fine); Taichi-first applies only if/when this graduates to production per-frame kernels |

### Blindspots

| Risk | Guard |
| --- | --- |
| Certifying on the target problem (circular) | the benchmark is fixed above, with an independent character-theory multiplicity check |
| Silent symmetry breaking by the grid (route a) | measure the certified spectrum's degeneracy splitting as the resolution ladder climbs; report it |
| Band-limit truncation masquerading as physics (route b) | convergence in the band limit reported for every observable |

### Ownership + gating

Author-driven with platform support. Gated by M8.2 (so the engine is built against
locked requirements, not drifting ones).

## DEVIATIONS LOG

(none)

## FINDINGS

(pending)
