# M5.21.13: the stiffness ladder (does the census survive δ = 1, 3, 10?)

> Task **M5.21.13** (M5 / Liquid-Crystal model). Status: 🚧 **PLANNED STUB** · Roadmap:
> [`../m5_roadmap.md`](../m5_roadmap.md) · Migrated from GitHub issue
> [#324](https://github.com/openwave-labs/openwave/issues/324) (opened 2026-07-23 by
> `tekemperor`) on 2026-08-01, when
> [T5](../../../../../dev_docs/tasks/t5_task_details.md) settled that tasks live in roadmaps and
> issues are reserved for platform defects. The proposal is archived below in full, so closing the
> issue loses nothing.

This doc is the task's full record: planning, then findings at the run.

## PLANNING

### The proposal

A cheap robustness probe of the lepton census in the stiffness direction, run on the existing
instrument before any 4×4 upgrade is paid for. Everything is held fixed from the
[M5.21.2b](m5_21_2b_task_details.md) converged census (48³ grid, both stencils, pinned and free
boundaries, FIRE settings, seeds A/B/C); only δ moves, with the potential trace targets updated to
the spectrum `{1, δ, 0}` at each rung. First pass at fixed `w`; a second pass rescales `w` if the
Derrick balance `u/3V` drifts far from 1.

| Rung | Purpose |
| --- | --- |
| δ = 0.3 | control, should bit-reproduce the census |
| δ = 1.0, 3.0, 10 | the ladder, as far as descent still converges |

Recorded at each rung: the A/C/B energy ordering, the energy ratios, the compact-component
topology signatures, whether any seed drains to vacuum, and the cross-stencil energy ratio.

### Verdicts, as the proposal pre-registers them

| Outcome | Reading |
| --- | --- |
| Ordering stable across all rungs and both stencils | continuation evidence toward the stiff vacuum |
| Ordering flips, or a seed drains at some δ | the toy corner is not predictive, and the rung where it breaks is the diagnostic |
| Cross-stencil disagreement growing with δ | the discretization inconsistency worsens with stiffness, which is worth knowing before building the 4×4 |

Cost as estimated by the proposal: three to six census re-runs on the existing grid, no new code
beyond parameter changes.

### PREMISE CHECK, to resolve before running

The proposal opens with "the physical regime needs g, delta ~ 10^10". The regime this column works
in is `0 < δ ≪ 1 ≪ g`, and [M5.21.11](m5_21_11_task_details.md) states the physical corner as
δ ~ 1e-10 with g ~ 1e10, so **δ = 1, 3, 10 walks away from that corner rather than toward it**: at
δ → 1 the δ-axis becomes degenerate with the middle axis and the biaxial hierarchy that produces
three lepton levels is gone by construction. Two readings, and the run is scoped differently under
each:

| Reading | What the ladder then is |
| --- | --- |
| The proposal means the stiffness SCALE (g, or the ratio), and δ is a slip | then the rungs belong on the g axis, where [M5.21.8](m5_21_8_task_details.md) already walked g = 8-64 and found the dressed-minimum ratio stable at 0.82-0.84 tracking the 1/g law |
| The proposal means δ literally | then it is a degeneracy probe, not a physical-regime continuation: a measurement of how the census ordering dies as biaxiality is removed. Still informative, and the honest label for it is different |

Related evidence already on disk, in the opposite direction: the [M5.22.1](m5_22_1_task_details.md)
ladder walked δ = 0.3 → 0.2 → 0.1 on the nuclear observables and found the ratios do NOT converge
by δ-steps alone ([note](../findings/m5_22_1_note.md) § 2). Whoever runs this task reads that first,
since it is the nearest measurement of what a δ walk does to a census result.

Route: ask the proposer which reading is meant before the first rung. It costs one comment and
decides the axis.

### POST-M5.21.11 RE-SCOPE (2026-08-08, user-approved): the premise check is largely answered

The [M5.21.11](m5_21_11_task_details.md) ladder run (2026-08-07, terminal on F3 + F4) resolves most
of this task's premise questions from data rather than from the proposer:

| Element above | State after M5.21.11 |
| --- | --- |
| Reading (a): the g axis | Fully covered without running anything here: [M5.21.8](m5_21_8_task_details.md) walked the m\* position law and the M5.21.11 g-arm measured the dressing GAIN, flat in g across 8-32. Nothing left for this task on that axis |
| Reading (b): literal δ upward | The surviving scope, and the honest label is a DEGENERACY PROBE: at δ → 1 the target spectrum (0, δ, 1) degenerates and past δ = 1 the hierarchy inverts; the run measures how the census dies as biaxiality is removed, nothing about the physical corner |
| Verdict "ordering stable = continuation evidence toward the stiff vacuum" | NOT earnable at any outcome: M5.21.11 measured the 4D dressing term O(1), branch-dependent, and g-flat, so no 3×3-only census walk in any δ direction is continuation evidence toward the physical vacuum. That verdict row is dead and must not be claimed from a run of this task |
| Verdict "cross-stencil growth with δ" | Partially pre-answered downward: the ladder measured xstencil at 8 rungs (δ 0.30 → 0.05): B/C sit at 2.1-2.5 (beyond the 1.5 bar) at EVERY rung, roughly δ-flat, A at ~1.03-1.1. Only the upward direction stays unmeasured |
| The "cheap test before the 4×4 investment" motivation | Moot: M5.21.11 settled that the 4×4 (dressing carried inside the ladder) is REQUIRED for any physical-regime claim regardless of stiffness-direction behavior |
| Nearest related measurement | Now the M5.21.11 E(δ) record itself ([§ LADDER RUN](m5_21_11_task_details.md)): 8 gated rungs, E rising as δ falls on all branches, A < C < B holding at every rung, and the C-B margin collapsing downward (C/B 0.28 → 0.85 over δ 0.3 → 0.05, uncertified energies); the [M5.22.1](m5_22_1_task_details.md) pointer stays as secondary |

Net: the task stays a valid, cheap, non-gating stub, but its scope shrinks to the reading-(b)
degeneracy probe and two of the three pre-registered verdict rows are relabeled as above before any
run. The proposer courtesy question stands, with the answer no longer changing the scope much.

**Gated by**: the premise check above + user "go". Non-gating for the electron program: the census
verdicts stand on their own arena, and this measures how far they travel.

## GitHub issue archive (#324)

> Migrated from OpenWave GitHub issue
> [#324](https://github.com/openwave-labs/openwave/issues/324) on 2026-08-01. Title: "Stiffness
> ladder: does the M5.21.2 lepton census survive delta = 1, 3, 10?". Opened 2026-07-23 by
> `tekemperor`. State at migration: OPEN. Labels: `help wanted`. Body verbatim.

The physical regime needs g, delta ~ 10^10 while current runs use delta = 0.3. Before investing in
the 4x4 upgrade, it's cheap to test whether the census results are stable in the stiffness
direction at all.

Protocol:

1. Hold everything fixed from the M5.21.2 census: 48^3 grid, both stencils, pinned and free
   boundaries, FIRE settings, seeds A/B/C.
2. Re-run the census at delta = 0.3 (control, should bit-reproduce), then 1.0, 3.0, and 10 if
   descent still converges. Update the potential trace targets to the new spectrum {1, delta, 0}.
   First pass with w fixed; if the Derrick balance u/3V drifts far from 1, a second pass rescaling
   w.
3. At each rung record: the A/C/B energy ordering, the energy ratios, the compact-component
   topology signatures, whether any seed drains to vacuum, and the cross-stencil energy ratio.
4. Verdicts, written before running: ordering stable across all rungs and both stencils = first
   real continuation evidence toward the stiff vacuum. Ordering flips or a seed drains at some
   delta = the toy corner is not predictive, and the rung where it breaks is the diagnostic.
   Cross-stencil disagreement growing with delta = the discretization inconsistency worsens with
   stiffness, important to know before building the 4x4.

Cost: three to six census re-runs on the existing grid, no new code beyond parameter changes.

## DEVIATIONS LOG

(none)

## FINDINGS

(pending: the task has not been run)
