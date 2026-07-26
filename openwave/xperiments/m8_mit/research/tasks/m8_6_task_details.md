# M8.6: The McKay-distance rule vs M5's measured lepton hierarchy

> Roadmap row: [`../m8_roadmap.md`](../m8_roadmap.md). Status: 🚧 PLANNED, ungated.
> This is a scaffold-stage planning aid written by the maintainers (2026-07-21); the
> author owns the MIT side; the M5 side is graded by the platform's M5 record. Joint
> task.

## PLANNING

### Scope

A bounded cross-check needing NO simulation: does MIT's McKay-distance rule reproduce
the lepton hierarchy that M5 measures but cannot yet derive? M5's record
([`../m8_platform_pointers.md § 5`](../m8_platform_pointers.md)): three rotation
minima with the eigenvalue hierarchy `1 : 5.9 : 15.1` (the open "hierarchy origin" of
the M5 lepton row), the mass law `E ∝ Λ³` already fixed, physical ratios
`1 : 206.8 : 3477.2`. MIT's candidate mechanism: mass ratios from
`(√Ω)^(dist/30)`-type structure at McKay distances. This is the ONE candidate
mechanism currently on the table for exactly that open item; either outcome closes a
live question in TWO columns.

### The pre-registration (the task's entire integrity)

All three choices below are made and frozen BEFORE any number is computed:

| Choice | To fix in advance |
| --- | --- |
| The mapping | which McKay slots/distances correspond to (e, μ, τ); justified structurally (generation = flat connection, per the MIT spec), not selected by fit |
| The comparison level | eigenvalue-level (`1 : 5.9 : 15.1`, then cubed by M5's `E ∝ Λ³`) vs mass-level (`1 : 206.8 : 3477.2`); one is primary, stated in advance |
| The tolerance | what counts as "reproduces" (a stated relative-error threshold) and what counts as refuted |

### Suggested definition of done

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
numbers and the audit). Ungated; can run any time, including before M8.1 closes.

## DEVIATIONS LOG

(none)

## FINDINGS

(pending)
