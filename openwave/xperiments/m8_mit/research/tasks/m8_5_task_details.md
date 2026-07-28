# M8.5: Quotient-manifold simulation engineering

> Roadmap row: [`../m8_roadmap.md`](../m8_roadmap.md). Status: 🚧 PLANNED; its gate
> M8.2 closed ✅ 2026-07-27. This is a scaffold-stage planning aid written by the
> maintainers (2026-07-21, § "Independent reproduction" added 2026-07-28); the author
> owns the column and may amend everything here.

## PLANNING

### Scope

Two sub-deliverables that the M8.2 lock put in the same task, named and tracked separately
below because they have different owners and different independence rules. They stay inside
M8.5 by agreement (author, [PR #350](https://github.com/openwave-labs/openwave/pull/350)
close-out 2026-07-28), and split into their own task IDs only if ownership or scheduling
diverge enough to make the roadmap unclear:

| # | Sub-deliverable | Owner |
| --- | --- | --- |
| M8.5-A | Independent table and decomposition reproduction (lock § 3) | protocol: author; implementation: maintainers |
| M8.5-B | The quotient backend: the simulation engine for the arena (routes a/b below) | author + platform support |

**A gates any claim that B is certified.** B can be built and can run before A closes; what
it cannot do before then is carry a certified-against-the-McKay-structure claim.

Build and certify the simulation infrastructure M8.4 needs: fields evolving on the
compact quotient S³/2I. No existing OpenWave column runs a curved compact arena, so
this is genuinely new platform ground. Two candidate routes
([`../m8_platform_pointers.md § 6`](../m8_platform_pointers.md)); prototype BOTH far
enough to choose one on evidence.

| Route | Sketch | Known risks |
| --- | --- | --- |
| (a) 2I-equivariant grid | an S³ grid (embedding or intrinsic charts) with the 120-element identification imposed as an equivariance/ghost-cell map | the identification map bookkeeping; chart seams; where the Möbius edge / cone structure of the MIT arena lives on the grid |
| (b) Spectral in 2I-symmetric harmonics | expand fields in S³ harmonics restricted to 2I-invariant (or covariant) subspaces; evolve coefficients | nonlinear terms need convolution handling (cost grows fast with band limit); but the basis IS the McKay representation theory, so slot structure is manifest |

### Independent reproduction, M8.5-A (added 2026-07-28)

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

#### The three objects, and what M8.5-A's implementer may read

The protocol distinguishes three objects rather than two, because the review-time verification
is now a repository artifact rather than a throwaway (landed 2026-07-28 at the author's request):

| # | Object | Method | Quarantine status for M8.5-A |
| --- | --- | --- | --- |
| 1 | [`../scripts/m8_2_first_occurrence.py`](../scripts/m8_2_first_occurrence.py) | McKay / affine-E8 recursion; irrep labels, dims, distances as literals | FORBIDDEN until step 2's source + raw output are committed |
| 2 | [`../scripts/m8_2_indep_reconstruction.py`](../scripts/m8_2_indep_reconstruction.py) + [note](../findings/m8_2_indep_reconstruction_note.md) | explicit quaternions → conjugacy classes → Burnside class-sums; dims, adjacency, distances derived | FORBIDDEN until step 2's source + raw output are committed |
| 3 | the M8.5-A implementation | the author's frozen protocol, fresh context | the deliverable |

After 3 is committed, 1 and 2 become **adjudication references only**. Agreement of all three
is then reported as **three-way agreement**: it strengthens provenance, and it does NOT move the
claim label off independent-method reproduction.

**What object 2 does not cover** (the author's scope point, 2026-07-28): it reconstructs 2I, its
characters, the McKay distances and the scalar (0-form) first-occurrence table. It does not derive
the coexact one-form entry rule, which stays **ASSERTED**. If M8.5-A's object needs the coexact
table, its entry rule has no independent target in the repository and the protocol has to say
what standing it can be given.

Object 2 was landed to the requirements the author set for it: source, environment, raw output,
the commit verified against, a short method note, a mutation test for every PASS line, and an
explicit statement of what it does and does not verify.

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
| M8.5-A written up as blind, or A's implementer reusing the M8.2 context | the claim ceiling and the firewall are stated above; the protocol (step 1) fixes both before step 2 opens |
| A's implementer reading either existing artifact early | both are named and quarantined above; they open only after A's own source + raw output are committed |
| M8.5-B claiming certification while A is open | A gates that claim, stated in § Scope; B may run, it may not claim |
| A PASS line that cannot go red | mutation-test every gate before it ships ([roadmap § CONVENTIONS](../m8_roadmap.md#conventions); the M8.2 defect) |

### Ownership + gating

M8.5-B is author-driven with platform support. M8.5-A splits: the author writes the
protocol, the maintainers implement it. Gated by M8.2 ✅ (so the engine is built against
locked requirements, not drifting ones); A's step 2 additionally waits on A's step 1, and
any certification claim for B waits on A.

## DEVIATIONS LOG

(none)

## FINDINGS

(pending)
