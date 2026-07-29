# M8.6: MIT-M5 lepton-hierarchy comparison — gated readiness audit

> Roadmap row: [`../m8_roadmap.md`](../m8_roadmap.md). Status: 🚧 GATED (2026-07-29),
> pending M5.21.11. Moved from Backlog to [LATER (gated)](../m8_roadmap.md#later-gated)
> after a readiness audit found the named target circular; see
> [`findings/m8_6_readiness_note.md`](../findings/m8_6_readiness_note.md).
> This is a scaffold-stage planning aid written by the maintainers (2026-07-21); the
> author owns the MIT side; the M5 side is graded by the platform's M5 record. Joint
> task.

## ORIGINAL PLAN — NOW GATED

> The plan below is the maintainers' original scaffold (2026-07-21). A readiness audit
> (2026-07-29, see FINDINGS below) found it cannot run as written: the named target is
> circular and the fallback is non-independent. Kept verbatim as the historical record;
> do not read it as the live plan — [`findings/m8_6_readiness_note.md`](../findings/m8_6_readiness_note.md) is.

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

### Original preregistration design — not executed

All three choices below are made and frozen BEFORE any number is computed:

| Choice | To fix in advance |
| --- | --- |
| The mapping | which McKay slots/distances correspond to (e, μ, τ); justified structurally (generation = flat connection, per the MIT spec), not selected by fit |
| The comparison level | eigenvalue-level (`1 : 5.9 : 15.1`, then cubed by M5's `E ∝ Λ³`) vs mass-level (`1 : 206.8 : 3477.2`); one is primary, stated in advance |
| The tolerance | what counts as "reproduces" (a stated relative-error threshold) and what counts as refuted |

### Original definition of done — superseded

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
