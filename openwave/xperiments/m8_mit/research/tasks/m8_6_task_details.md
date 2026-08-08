# M8.6: MIT-M5 lepton-hierarchy comparison (gated readiness audit)

> Roadmap row: [`../m8_roadmap.md`](../m8_roadmap.md). Status: ❌ CLOSED WITHOUT RUNNING
> (2026-08-07; row retired to [Done](../m8_roadmap.md#done) 2026-08-08, maintainer call).
> Lineage: moved from Backlog to LATER (gated) 2026-07-29 after a readiness audit found
> the named target circular ([`findings/m8_6_readiness_note.md`](../findings/m8_6_readiness_note.md));
> closed when the last admissible route failed terminally (finding 9).
> This is a scaffold-stage planning aid written by the maintainers (2026-07-21); the
> author owns the MIT side; the M5 side is graded by the platform's M5 record. Joint
> task.

## ORIGINAL PLAN (NOW GATED)

> The plan below is the maintainers' original scaffold (2026-07-21). A readiness audit
> (2026-07-29, see FINDINGS below) found it cannot run as written: the named target is
> circular and the fallback is non-independent. Kept verbatim as the historical record;
> do not read it as the live plan: [`findings/m8_6_readiness_note.md`](../findings/m8_6_readiness_note.md) is.

### Original scope

A bounded cross-check needing NO simulation: does MIT's McKay-distance rule reproduce
the lepton hierarchy that M5 measures but cannot yet derive? M5's record
([`../m8_platform_pointers.md § 5`](../m8_platform_pointers.md)): three rotation
minima with the eigenvalue hierarchy `1 : 5.9 : 15.1` (the open "hierarchy origin" of
the M5 lepton row), the mass law `E ∝ Λ³` already fixed, physical ratios
`1 : 206.8 : 3477.2`. MIT's candidate mechanism: mass ratios from
`(√Ω)^(dist/30)`-type structure at McKay distances. This is the ONE candidate
mechanism currently on the table for exactly that open item; either outcome closes a
live question in TWO columns.

### Original preregistration design (not executed)

All three choices below are made and frozen BEFORE any number is computed:

| Choice | To fix in advance |
| --- | --- |
| The mapping | which McKay slots/distances correspond to (e, μ, τ); justified structurally (generation = flat connection, per the MIT spec), not selected by fit |
| The comparison level | eigenvalue-level (`1 : 5.9 : 15.1`, then cubed by M5's `E ∝ Λ³`) vs mass-level (`1 : 206.8 : 3477.2`); one is primary, stated in advance |
| The tolerance | what counts as "reproduces" (a stated relative-error threshold) and what counts as refuted |

### Original definition of done (superseded)

| # | Item |
| --- | --- |
| 1 | The frozen pre-registration block written into this doc BEFORE numerics |
| 2 | A short script (`scripts/m8_6_mckay_hierarchy.py`) computing the McKay-side ratios from group theory (no quoted constants) |
| 3 | Verdict either way, adversarially audited, wired into BOTH columns (the M5 question tracker's hierarchy item AND the M8 lepton cell) |

### Blindspots

| Risk | Guard |
| --- | --- |
| Post-hoc mapping (trying assignments until one lands) | the mapping is frozen first; if the frozen mapping fails, that IS the result; alternative mappings may be reported afterwards but only labeled as exploratory |
| Double freedom (choosing mapping AND comparison level after seeing numbers) | both frozen; the secondary level is reported but cannot rescue a failed primary |
| Carrying MIT's torsion-map weight into this test | the McKay DISTANCE structure is the ledger's stronger layer; keep T out of the primary comparison or justify its role in advance |

### Ownership + gating

Joint (author supplies the MIT-side mapping rationale; the platform supplies the M5
numbers and the audit). Originally ungated; GATED 2026-07-29 pending M5.21.11 (see
FINDINGS below and [`../m8_roadmap.md § LATER (gated)`](../m8_roadmap.md#later-gated)).

## DEVIATIONS LOG

| # | Deviation | Reason |
| --- | --- | --- |
| 1 | Before writing the pre-registration, the task's first work was a provenance/readiness audit of the named M5 target, not the numerical run itself. | The planning doc's own integrity requirement ("justified structurally, not selected by fit") made tracing the target's provenance a precondition, not an optional check; the audit found the plan could not proceed as scoped. |

## FINDINGS

1. **The named target is circular.** `1:5.9:15.1` (cited via [`../m8_platform_pointers.md § 5`](../m8_platform_pointers.md)) is not independent M5 output: M5's own findings note defines it as `Λ := m^(1/3)`, states reproducing the masses this way is "near-tautological... a consistency check, not a parameter-free prediction," and that the eigenvalue values "remain Yukawa-like input." `5.9` and `15.1` are the cube roots of the known muon/tau-to-electron mass ratios, rounded to two figures.

2. **A genuine measured alternative exists, but isn't usable yet.** M5.21.2/2b independently measured three stationary-state energies (`A<C<B`, `C/A≈4.2`, `B/A≈16.0`) with a physically-motivated ordering (lowest-energy candidate state; a separately-measured decay mechanism identifying the μ- and τ-candidates). But the source itself states these are "consistency-converged, not value-converged" (E drifts 7-13% between grid rungs) and carry no frozen physical-parameter and units bridge to physical mass ("the voxel → fm anchor is Q17, unset"); the toy-to-physical calibration is explicitly deferred to M5.21.11.

3. **An ordering-only fallback was considered and rejected.** MIT's own charged-lepton identities (`e=(R7,triv), μ=(R8,std), τ=(R4,gal)`) were assigned in M8.3 by matching to measured PDG masses (`mass-spectrum.md`: "the gates fix the kind; the mass fixes the generation"). `m_e<m_μ<m_τ` is therefore already built into which slot carries which label; checking whether the inherited triple comes out light-middle-heavy is true by construction, not a finding.

4. **Verdict: not yet well-posed, on either side.** Full provenance table, exact source citations, and the reopening conditions: [`findings/m8_6_readiness_note.md`](../findings/m8_6_readiness_note.md).

5. **Condition 3 amended, 2026-07-29.** Cross-linking the gate surfaced a real tension: the original condition 3 admitted only a direct run at physical parameters, but M5.21.11's own scope states that regime (δ ~ 1e-10, g ~ 1e10) is out of lattice reach by any direct method: unsatisfiable as written. The readiness note's § 8 now also admits a preregistered extrapolation route, under guardrails (theory-derived asymptotic form, frozen rung set/fitting/holdout/uncertainty model, the three existing toy energies barred from the new fit) that keep the same anti-circularity requirement the original condition enforced.

6. **Condition 4 strengthened, 2026-07-29.** Review of the condition-3 amendment found its uncertainty requirement, as first drafted, covered only extrapolation error along the g/δ ladder, leaving out the per-rung discretization error already visible at fixed δ (the readiness note's own finding, § 5, that `E_A` drifts ~20% across three grid resolutions at fixed δ, with B and C less consistency-converged than A). The readiness note's § 8 condition 4 now requires, under route (b), a per-rung discretization term established by grid refinement of all three branches on a preregistered subset of rungs spanning the ladder, together with a frozen rule for propagating that term to unrefined rungs, alongside the extrapolation-uncertainty term, since branch-dependent discretization error cannot be assumed to cancel in the ratios.

7. **The route-(b) pre-registration artifact delivered, 2026-08-06.** M5.21.11's P0 run produced [`m5_21_11_framework.md`](../../../m5_liquid_crystal/research/findings/m5_21_11_framework.md): the derived asymptotic forms, the frozen ladder protocol (rung set, holdouts, the condition-4 refinement subset + propagation rule, branch rules, uncertainty model, terminal failure criteria), a barred-inputs audit, and a § 8 table mapping each reopening condition to the section satisfying it. Status at delivery: PRE-FREEZE (the artifact's own header is the live status; it freezes with a pinning commit SHA after an M5-side author sanity-check of the derived forms). Two notes for the M8 side: the refinement subset {δ = 0.30, 0.12, 0.05} × N ∈ {32, 48, 64} × all three branches satisfies condition 4's "subset spanning the ladder" literally, including both endpoint rungs, because rung cost is δ-independent (the affordability tension anticipated in the M5.21.11 task doc did not materialize; no further amendment is requested); and this gate does NOT reopen at the freeze: reopening still requires the ladder run to complete under the frozen framework with its holdout and failure gates passed.

8. **The framework FROZEN, 2026-08-07.** The M5-side author sanity-check was answered 2026-08-07 (a public reply endorsing the practical-approximation premise of route (b), with no objection to any derived form; the three specific form questions went unanswered, and the frozen fit covers every branch of them), and the framework froze on the M5-side user call the same day: the artifact's § 0 header is FROZEN with the pinning commit SHA recorded there. Gate state after the freeze: the pre-registration precondition of the amended condition 3 is now satisfied and immutable (any further edit voids route (b)); what remains before this gate reopens is unchanged from item 7: the ladder run itself, completed under the frozen framework with its holdout and terminal failure gates passed. On that outcome the roadmap row is released from [LATER (gated)](../m8_roadmap.md#later-gated) back to the Backlog (recorded here when it happens); a terminal ladder failure leaves it gated, route (b) being the last admissible route.

9. **The ladder ran and route (b) failed TERMINALLY, 2026-08-07: this gate is closed permanently.** The frozen program executed the same day the framework froze (M5-side record: [`m5_21_11_task_details.md § LADDER RUN`](../../../m5_liquid_crystal/research/tasks/m5_21_11_task_details.md)). Two pre-committed terminal criteria fired: the g-arm measured the boost-dressing gain FLAT in g (the framework's one-dimensional-in-δ reduction fails, so the extrapolation object itself was ill-founded), and branch integrity collapsed under the frozen instrument-health gates (1/0/0 usable rungs for A/C/B against the 6-per-branch floor; the two heavier branches, the μ- and τ-candidates this comparison needed most, were uncertifiable at every rung). No fit ran, no ratio was produced, and the failure is terminal under the framework's own § 6: route (a) was already out of lattice reach (finding 5), route (b) was the last admissible route, and the amended condition 3 admits no third. The M8.6 comparison is therefore not merely gated but CLOSED on the M5 side as instrumented; it could only reopen if a future M5 platform delivered a fundamentally different certified route (a new instrument generation, not a new analysis of this data), which no existing task promises. This finding is the closing record; the roadmap row was retired to [Done](../m8_roadmap.md#done) on 2026-08-08 (maintainer call, CLOSED WITHOUT RUNNING).
