# M8.2: Pre-registration lock for the field-dynamics program

> Roadmap row: [`../m8_roadmap.md`](../m8_roadmap.md). Status: 🚧 PLANNED. This is a
> scaffold-stage planning aid written by the maintainers (2026-07-21); the author owns
> the column and may amend everything here. The go-time pre-registration below MUST be
> completed and frozen before any M8.4 numerics run.

## PLANNING

### Scope

Write and freeze the pre-registration document that the whole field-dynamics program
(M8.4) will be graded against. This task produces no numerics; its deliverable is a
locked set of targets, criteria, and conventions, so that whatever M8.4 finds is a
result rather than a fit ([`AI_HYGIENE.md`](../../../../../AI_HYGIENE.md);
[`ONBOARDING_MODELS.md § 4`](../../../../../ONBOARDING_MODELS.md)).

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
| 1 | `m8_2_preregistration.md` written, reviewed by both the author and a maintainer, and frozen (dated; later edits go in a dated addendum, never in-place) |
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

(none)

## FINDINGS

(pending)
