# M8.2: Pre-registration lock for the field-dynamics program

> Roadmap row: [`../m8_roadmap.md`](../m8_roadmap.md). Status: ✅ DONE (2026-07-27).
> The deliverable is [`../findings/m8_2_preregistration.md`](../findings/m8_2_preregistration.md),
> LOCKED at that date via [PR #350](https://github.com/openwave-labs/openwave/pull/350);
> later changes go in a dated addendum there, never in-place. Everything below the
> PLANNING heading is the original scaffold-stage planning aid written by the
> maintainers (2026-07-21), kept for the record; the outcome is in
> [`## TASK REVIEW`](#task-review-2026-07-27).

## PLANNING

### Scope

Write and freeze the pre-registration document that the whole field-dynamics program
(M8.4) will be graded against. This task produces no numerics; its deliverable is a
locked set of targets, criteria, and conventions, so that whatever M8.4 finds is a
result rather than a fit ([`AI_HYGIENE.md`](../../../../../AI_HYGIENE.md);
[`ONBOARDING_MODELS.md parameter-count test`](../../../../../ONBOARDING_MODELS.md)).

### Why it exists as its own task

The platform's recorded failure modes (the M6 provenance ledger) all trace to
conventions chosen after seeing numbers: calibrated conventions, window-defined
observables, targets that drift with the result. One locked page prevents the whole
class.

### What must be fixed in the lock (the go-time checklist)

| Item | To fix BEFORE any run |
| --- | --- |
| Target observables | the STRUCTURAL ladder only: McKay slot structure, gap ratios, generation count (3 flat connections). The 24-entry numeric mass table is explicitly OUT of scope ([`../m8_background.md § 3`](../m8_background.md)) |
| Success criterion | defect / standing-wave energies proportional to the McKay slot values WITHOUT per-slot tuning: one global scale + the Lagrangian's own couplings, nothing per-particle |
| Comparison level | eigenvalue-level vs energy-level comparison, chosen and justified in advance |
| Lagrangian families in scope | the candidate list (from [`../m8_platform_pointers.md § 2`](../m8_platform_pointers.md)) with each family's free couplings enumerated and bounded |
| The no-search rule | every (family, coupling) point run is reported; forks reported with all numbers; a scan is declared as a scan with its grid stated up front |
| Vacuum gate | for each family, the vacuum spectrum on the arena is computed and published BEFORE soliton hunting (the M7 tachyonic-band lesson) |
| Failure criteria | what counts as the family REFUTED on this arena (so negatives close cleanly) |

### Suggested definition of done

| # | Item |
| --- | --- |
| 1 | [`../findings/m8_2_preregistration.md`](../findings/m8_2_preregistration.md) written, reviewed by both the author and a maintainer, and frozen (dated; later edits go in a dated addendum, never in-place) |
| 2 | Each M8.4-scope family has its couplings + bounds + failure criteria enumerated |
| 3 | The lock is linked from the roadmap and from `m8_theory_canonical.md` |

### Blindspots

| Risk | Guard |
| --- | --- |
| Registering the mass table anyway ("just as a secondary check") | banned outright; the author's own null tests already cap its evidential weight |
| Vague success wording ("roughly matches the ladder") | the criterion must be numeric: a stated statistic with a stated threshold |
| Silent post-hoc edits to the lock | frozen file + dated addenda only |

### Ownership + gating

Author-driven, maintainer-reviewed. Gated by M8.1 (no point locking a program against
an arena whose headline eigenvalue failed certification).

## DEVIATIONS LOG

| Deviation from the scaffold plan | Why |
| --- | --- |
| The lock is MODULAR (immutable core §§ 1-5 + per-family modules § 6 + signed per-family execution appendices § 7) rather than one flat page | A single page cannot freeze three families whose operator questions resolve at different times. The core locks now; each family's numerics freeze in its own signed appendix before its first target-bearing run |
| "Comparison level" and "vacuum gate" are not single pre-declared choices | Both became per-family: § 5 records background admissibility and vacuum stability as separate outcome axes per sector, and § 4's ladder makes the vacuum a computed output rather than a gate to declare in advance |
| Three-connection compatibility is recorded as "not applicable" for all three native families | Downstream of the platform standing rule this task produced ([`CROSS_MODEL_TESTING.md`](../../../../../dev_docs/CROSS_MODEL_TESTING.md)): a borrowed family is native and untwisted unless its author declares otherwise. Two authors declared, and both declarations confirmed the default |

## FINDINGS

| # | Finding |
| --- | --- |
| 1 | **The lock exists and is modular.** Core contract (arena, bundles, the two preconditions, the certification gate, the success ladder, four-axis outcome language, the no-search rule) frozen 2026-07-27; per-family numerics deferred to signed execution appendices (§ 7) that gate M8.4, not this lock |
| 2 | **The index question became a platform rule.** The audit behind this task produced [`dev_docs/CROSS_MODEL_TESTING.md`](../../../../../dev_docs/CROSS_MODEL_TESTING.md) (2026-07-24): uniform default, what counts as a declaration, the soldering clause, "not applicable" as neutral, author silence terminal. No family blocks on an inbox |
| 3 | **Two author declarations landed, both confirming the default.** @jeffsyee: M4's `ψ` is a geometric spatial displacement, not an internal triplet ([#333](https://github.com/openwave-labs/openwave/discussions/333)). @JarekDuda: M5's `M` is a level-3 spacetime/frame tensor, with latitude to use M5 unmodified ([#334](https://github.com/openwave-labs/openwave/discussions/334)). The latitude is read as level-choice, NOT a soldering endorsement, so "M5 + P" stays M8-owned |
| 4 | **The certification table is independently reproduced.** All 9 rows of the § 6.1 first-occurrence table verified by a maintainer implementation sharing no method with the shipped script (2I built as 120 explicit unit quaternions, conjugacy classes by brute-force conjugation, irreducible characters by Burnside class-sum diagonalization, versus the script's McKay recursion). Class sizes, irrep dimensions, `Sym²Q` / `Sym²Q'` and the McKay distances all recovered independently |
| 5 | **Provenance verified.** Both pinned SHAs resolve as described, and the § 1 pin audit is accurate: the M4 engine, both M5 engine files, the M7 functional and the M8.1 note are byte-identical at `c9dc3796`, with the M5 canonical's drift correctly characterized as documentation-only |

## TASK REVIEW (2026-07-27)

Contributed by the column's author (Blake Shatto) as
[PR #350](https://github.com/openwave-labs/openwave/pull/350), reviewed and closed by the
maintainer under [`dev_docs/PR_REVIEW_STANDARDS.md`](../../../../../dev_docs/PR_REVIEW_STANDARDS.md).

| Result | Status |
| --- | --- |
| Core contract + family modules locked | ✅ |
| Certification table independently reproduced, 9/9 rows | ✅ |
| Author declarations obtained and correctly scoped (M4, M5) | ✅ |
| Provenance pins verified against upstream | ✅ |
| Native quotient operators for M4, M5, M7 | 🚧 execution-appendix work under § 7; gates M8.4, not this lock |
| Zenodo-deposit byte-check of the mass-spectrum source | 🚧 local SHA-256 recorded; deposit cross-check open |

**Maintainer edits applied at merge.** Announced in the PR thread rather than sent back as a
review round: the pre-registration moved to `findings/` (frozen artifacts live there; the
living research docs stay at `research/` root) with its relative links repointed; two date
slips corrected (a pin audit dated a day forward, and the @JarekDuda follow-up dated 07-25
when [#334](https://github.com/openwave-labs/openwave/discussions/334) shows 07-26); the
status header changed from WORK IN PROGRESS to LOCKED with an explicit scope of what the lock
does not cover; `black` applied to the script; the script docstring rewritten to match §§ 6.1
and 6.3, which disown both tables as native certification targets; and one self-check removed.

**On the removed self-check (the one substantive finding).** `m7_trivial_ok` compared the
coexact table's trivial column against `coexact_level(dist[...])`, but for the trivial
connection `m7_first_level` reduces to exactly that call, so both sides were the same
expression and the line printed PASS for any entry rule. Confirmed by replacing the rule with
deliberate nonsense: the check still passed while the table filled with wrong values. It is
replaced by two checks that can fail, and do: `dims` against the affine E8 mark condition
`A·dims = 2·dims`, and `dist` against BFS distance on the McKay graph. Both were mutation-tested.
The coexact entry rule itself is now labeled ASSERTED rather than certified, since unlike the
0-form table it has no published target to check against.

**Author close-out (2026-07-28, in the PR thread).** The author confirmed `m7_trivial_ok` as the
one real defect, accepted the two replacement checks and the ASSERTED label for the coexact entry
rule, and adopted "mutation-test every PASS line before it ships" as a standing rule of the
column. The thread also produced one forward constraint that is NOT an amendment to the lock: the
context that built M8.2 holds the target tables and fixtures, so it cannot serve as the
independent reproducer that lock § 3 requires. M8.5 therefore runs protocol-first (author writes
the frozen protocol including a context firewall; maintainers implement it in a fresh context)
and its result is reported as independent-method reproduction rather than blind. Recorded in
[`m8_5_task_details.md § Independent reproduction`](m8_5_task_details.md#independent-reproduction-deliverable-a-added-2026-07-28)
and the roadmap's [§ CONVENTIONS](../m8_roadmap.md#conventions). The lock's § 3 wording already required implementation
independence and stands unchanged, so no addendum was opened.

### Findings

M8.2 delivers a modular pre-registration whose core contract is frozen and whose per-family
numerics are deliberately not: the operator selections that M8.4 will actually test move to
signed execution appendices under § 7. The task's most reusable output is not the lock itself
but the platform rule it forced into existence, `CROSS_MODEL_TESTING.md`, which converted an
open cross-model index question into a default that resolves without blocking on author
replies, and then drew two author declarations that independently confirmed that default.

### Research docs created/updated

| Doc | What |
| --- | --- |
| [`../findings/m8_2_preregistration.md`](../findings/m8_2_preregistration.md) | The lock (core contract §§ 1-5, family modules § 6, execution appendices § 7, definition of done § 8) |
| [`../scripts/m8_2_first_occurrence.py`](../scripts/m8_2_first_occurrence.py) | The first-occurrence tables and their self-checks |
| [`../m8_roadmap.md`](../m8_roadmap.md) | M8.2 row moved to DONE; M8.4 / M8.5 remain gated |
| [`dev_docs/CROSS_MODEL_TESTING.md`](../../../../../dev_docs/CROSS_MODEL_TESTING.md) | The platform standing rule this task produced (landed separately, 2026-07-24) |
