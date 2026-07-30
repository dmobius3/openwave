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

**Step 1 ✅ LOCKED 2026-07-30**: [`../findings/m8_5a_reproduction_protocol.md`](../findings/m8_5a_reproduction_protocol.md),
author-written and filed before any implementation existed, landed through
[PR #380](https://github.com/openwave-labs/openwave/pull/380). It is the binding spec for step 2;
where this planning doc and the protocol differ, the protocol governs. Step 2 is startable.

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
the coexact one-form entry rule, which stays **ASSERTED**.

**That standing question is answered** (protocol § 6, 2026-07-30). The coexact rule is not a
certification target: it appears only as a separately labeled ASSERTED adjudication module,
optional and pre-declared in the § 9 commitment, with four pre-declared verdicts. The hard rule
is that numerical agreement in any amount never upgrades the rule's standing and is never
reported as independent verification; only a general argument from the operator and
representation structure does, and a pattern interpolated from the computed range does not
qualify however clean the fit.

Object 2 was landed to the requirements the author set for it: source, environment, raw output,
the commit verified against, a short method note, a mutation test for every PASS line, and an
explicit statement of what it does and does not verify.

#### What step 2 adds to the frozen minimum gate set (2026-07-30)

The protocol's § 8 minimum set (G1-G10) covers the 2I group theory completely, and it maps onto
six of object 2's eight checks. It does not gate the two constructions that sit between the group
and the answer, which are exactly the two the earlier mutation suite found worth testing. § 8
permits the implementer to add gates, so step 2 adds these and discloses them:

| Gate | Content | The defect it catches |
| --- | --- | --- |
| G11 | `χ_{Sym²σ}(g) = (χ_σ(g)² + χ_σ(g²))/2` on every class | object 2's `sym2_as_square`: `τ` built as `χ²` instead of `Sym²(χ)`. A dimension-only check does not discriminate; this identity does |
| G12 | `χ_{V_n}(e) = n+1` over the range actually searched | object 2's `chiv_offbyone`: the SU(2) character summed over `n` weights instead of `n+1` |

Why it matters that these are gated rather than left to chance: both defects pass G1-G10 and
land as a systematically shifted table, so they report as **partial disagreement** rather than
structural failure. That is a false negative on reproduction, the expensive direction, since the
fresh implementation reads as having failed to reproduce a correct table.

Raised in the [PR #380](https://github.com/openwave-labs/openwave/pull/380) review. Whether they
also enter the frozen floor as a § 11 dated addendum is the author's call; step 2 runs them
either way.

#### Operational items to settle before the clean room opens

| Item | What needs deciding | Working assumption |
| --- | --- | --- |
| Clean-room location, and who commits | § 4 forbids access to the M8 tree beyond the audited packet, while § 9 requires source, output and manifest committed before anything is unsealed, and committing needs the repository | the implementation runs OUTSIDE the working tree and a maintainer performs the commit. For an AI implementer this is load-bearing: a session started inside the repo auto-loads `CLAUDE.md` and its M8 pointers |
| Generic references | § 4 permits generic SU(2) representation-theory and harmonic-analysis references as construction inputs, then defines the audited packet as the protocol plus the group input, and § 6's module needs the implementer's own harmonic analysis | either they enter the packet and the audit covers them, or the maintainer-approved allowlist is the stated route |
| The § 7 comparison harness | it transcribes the quarantined § 6.1, so it is necessarily written after the commitment, and § 9 covers post-commitment changes only as a dated rerun | the harness is a second, separately dated commitment, so the ordering record stays unambiguous |

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
