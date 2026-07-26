# M8.1.1: SECOND BLIND RUN, the remaining bedrock theorems (gaps + asymmetry on S³)

> Roadmap row: [`../m8_roadmap.md`](../m8_roadmap.md) M8.1.1 (**maintainer-run**).
> Parent template: [`m8_1_task_details.md`](m8_1_task_details.md) (the worked blind-run
> protocol, reused as-is). Sources: the two bedrock papers the author shared on
> [#312](https://github.com/openwave-labs/openwave/discussions/312#discussioncomment-17758091)
> (2026-07-24): [SSRN 6968698](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=6968698)
> and [SSRN 7129118](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=7129118).

## TASK PLANNING (2026-07-24, registered; go pending)

### Scope

Blind independent verification, in the M8.1 sense (script and number by agents who have
not seen the derivations or the claimed values), of the two remaining MIT bedrock papers,
described by the author as establishing "gaps and asymetry on S3". Expected identification
from the briefing's key-inputs list: the **coexact spectral gap from McKay distance** (the
Yang-Mills mass-gap result on S³/2I) and the **E8-filling Galois pair** (the asymmetry
structure behind the three flat connections). This identification is a planning guess and
is CONFIRMED only at claim extraction: at planning time the papers were deliberately left
unread (see the quarantine note below).

### Why this task exists (the principled trigger, not courtesy)

| Point | Statement |
| --- | --- |
| Platform relevance | These theorems are target structure for the field-dynamics program: [M8.2](m8_2_task_details.md) pre-registers gap ratios among its observables and [M8.4](m8_4_task_details.md) aims the Lagrangian survey at the resulting slot structure. A pre-registered target should stand on verified structure. |
| The rule | Any analytic number that becomes a pre-registered target of the dynamics program gets blind-verified BEFORE M8.2 locks it. |
| Why maintainer-run BY CONSTRUCTION | The author's agents have read the papers; a self-run "verification" is a reproducer, not an independent recompute. Blindness is the one thing the platform can supply that the author structurally cannot. |
| What this task is NOT | A gate on the author's roadmap: M8.2, M8.3 and M8.6 are startable now and none of them wait for this run. It is also not a standing free-validation service; M8.1 was the certification gate, this run exists because of the pre-registered-target rule above. |

### Quarantine status (in effect since registration)

At registration (2026-07-24) the maintainers' planning session read ONLY the discussion
thread and the briefing's one-line paper descriptions, NOT the papers themselves; the
SSRN pages were not fetched. The M8.1 independence protocol therefore reuses cleanly at
go time:

| Role | Sees | Does not see |
| --- | --- | --- |
| Designer (go-time orchestrator) | the full papers (fetched to the session scratchpad, OUTSIDE the repo) | n/a |
| Solver agent | ONLY a self-contained spec sheet per paper (objects, operators, group data, boundary conditions) + task list | the claimed values, theorem statements' numeric content, the author's name/repo/papers, ALL repository docs |
| Audit agent | the same spec sheet + the solver's outputs + script | same blindness to the claims |
| Comparison to the claims | by the designer AFTER both agents return numbers | |

No-search rule: every computed number is reported; nothing is tuned toward the claims;
if the numbers land somewhere else, that IS the result.

### To be fixed at go time (BEFORE numerics)

| Item | Note |
| --- | --- |
| Paper identification confirmed | match the two SSRN IDs to the briefing's named bedrock inputs; update this doc |
| Pre-registered claims tables (C1..Cn per paper) | the M8.1 § "Pre-registered claims" format, one table per paper; sub-runs S-A (gap paper) and S-B (asymmetry paper) |
| Feasibility triage | the Möbius operator was 1D-reducible; a coexact gap on S³/2I carries representation theory and may need a spectral method in 2I-symmetric harmonics (see [`../m8_platform_pointers.md § 6`](../m8_platform_pointers.md)); if a claim is not boundable as a numerical check, say so and scope it out honestly |
| Blindspot pass | redo the M8.1 blindspot table for the new operators (sector completeness, quotient identification maps, seam/gluing analogues) |
| Citations sync | add both papers to [`../../theory/_CITATIONS.md`](../../theory/_CITATIONS.md) (or adopt the author's PR if it lands first) |

### Definition of done (skeleton, finalized at go)

| # | Item |
| --- | --- |
| 1 | Per-paper solver runs with converged numbers, scripts + JSON in the repo (`m8_1_1_` prefixes) |
| 2 | Adversarial audit with its own method, per-claim verdicts |
| 3 | Designer comparison against the pre-registered claims, all numbers stated |
| 4 | Method note `findings/m8_1_1_method_note.md` (equations first, eq-to-code map, embedded plots, audit record) |
| 5 | Doc sync: canonical + briefing + MODELS.md cells + roadmap row |
| 6 | Doc checker exit 0; TASK REVIEW presented |

### Scheduling

Maintainer-run at maintainer pace; registered 2026-07-24, run when maintainer resources
free up (communicated on #312: not in the next days). Not on the author's critical path.

## DEVIATIONS LOG

(none yet)

## FINDINGS

(pending run)
