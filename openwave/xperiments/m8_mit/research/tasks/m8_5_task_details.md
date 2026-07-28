# M8.5: Quotient-manifold simulation engineering

> Roadmap row: [`../m8_roadmap.md`](../m8_roadmap.md). Status: 🚧 PLANNED; its gate
> M8.2 closed ✅ 2026-07-27. This is a scaffold-stage planning aid written by the
> maintainers (2026-07-21, § "Independent reproduction" added 2026-07-28); the author
> owns the column and may amend everything here.

## PLANNING

### Scope

Two deliverables that the M8.2 lock put in the same task, tracked separately below
because they have different owners and different independence rules:

| # | Deliverable | Owner |
| --- | --- | --- |
| A | Independent reproduction of the M8.2 certification tables (lock § 3) | protocol: author; implementation: maintainers |
| B | The simulation engine for the quotient arena (routes a/b below) | author + platform support |

Build and certify the simulation infrastructure M8.4 needs: fields evolving on the
compact quotient S³/2I. No existing OpenWave column runs a curved compact arena, so
this is genuinely new platform ground. Two candidate routes
([`../m8_platform_pointers.md § 6`](../m8_platform_pointers.md)); prototype BOTH far
enough to choose one on evidence.

| Route | Sketch | Known risks |
| --- | --- | --- |
| (a) 2I-equivariant grid | an S³ grid (embedding or intrinsic charts) with the 120-element identification imposed as an equivariance/ghost-cell map | the identification map bookkeeping; chart seams; where the Möbius edge / cone structure of the MIT arena lives on the grid |
| (b) Spectral in 2I-symmetric harmonics | expand fields in S³ harmonics restricted to 2I-invariant (or covariant) subspaces; evolve coefficients | nonlinear terms need convolution handling (cost grows fast with band limit); but the basis IS the McKay representation theory, so slot structure is manifest |

### Independent reproduction, deliverable A (added 2026-07-28)

The M8.2 lock § 3 puts an obligation on this task:
[`../findings/m8_2_preregistration.md`](../findings/m8_2_preregistration.md) requires every
certification table to be reproduced through an INDEPENDENTLY implemented decomposition,
which may compare against
[`../scripts/m8_2_first_occurrence.py`](../scripts/m8_2_first_occurrence.py) but may not call
it, import its tables, or share its derived fixtures. The author added the procedure for it in
the [PR #350](https://github.com/openwave-labs/openwave/pull/350) close-out thread:

| Step | What | Who |
| --- | --- | --- |
| 1 | A frozen reproduction protocol: context firewall, operator conventions, result categories, mutation-tested gates, provenance requirements | the author, submitted as a document before any implementation exists |
| 2 | Implementation of that protocol in a FRESH context with no M8.2 internals loaded | maintainers |

Why the roles are this way round: the context that produced M8.2 holds the target tables and
the derived fixtures, so it cannot serve as its own reproducer no matter how separately the
second implementation is written. The author raised this rather than being asked.

**The claim ceiling is independent-method reproduction, not blind.** Blind means the verifying
agent has not seen the claimed values
([`ONBOARDING_MODELS.md § 3.2`](../../../../../ONBOARDING_MODELS.md#32-the-maintainer-sequence)
step 5), and that is already spent here twice over: the author's context built the tables, and
the maintainer reproduced all nine rows by explicit quaternions plus Burnside class-sums during
the PR review. An independent method still carries real weight, and the write-up says which one
it is. Do NOT let the word "blind" into the M8.5 deliverable.

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
| Deliverable A written up as blind, or A's implementer reusing the M8.2 context | the claim ceiling and the firewall are stated above; the protocol (step 1) fixes both before step 2 opens |
| A PASS line that cannot go red | mutation-test every gate before it ships (roadmap standing rules; the M8.2 defect) |

### Ownership + gating

Deliverable B is author-driven with platform support. Deliverable A splits: the author
writes the protocol, the maintainers implement it. Gated by M8.2 ✅ (so the engine is
built against locked requirements, not drifting ones); A's step 2 additionally waits on
A's step 1.

## DEVIATIONS LOG

(none)

## FINDINGS

(pending)
