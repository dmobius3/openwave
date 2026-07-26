# M8.3: The mass-formula reproducer script

> Roadmap row: [`../m8_roadmap.md`](../m8_roadmap.md). Status: 🚧 PLANNED. This is a
> scaffold-stage planning aid written by the maintainers (2026-07-21); the author owns
> the column and may amend everything here.

## PLANNING

### Scope

A standalone script that reproduces MIT's 24-entry fermion mass spectrum from first
definitions:

```text
m = μ_Λ · C_geom · (√Ω)^(dist/30) · T²
```

with EVERY constant recomputed from its own definition, never quoted from the papers:
McKay distances from the 2I McKay graph (buildable from the character table of the
binary icosahedral group), Reidemeister torsions from the three flat connections of
S³/2I, Kostant weights C_geom from their definition, and the anchors (μ_Λ from
measured Λ, m_e as normalization) declared as INPUTS. PDG values are the comparison
set. This is the [`ONBOARDING_MODELS.md § 6`](../../../../../ONBOARDING_MODELS.md)
"independent recomputer" pass applied to the analytic sector, the same category the
platform scores EWT's analytic masses under.

### Suggested sub-steps

| # | Step | Note |
| --- | --- | --- |
| 1 | Recompute the 2I character table + McKay graph from group theory; extract the McKay distance of each irrep | pure derivation, no author input needed |
| 2 | Recompute the Reidemeister torsion of each (irrep, flat connection) pair from its definition | if the published definition underdetermines a choice (basis, normalization), that is an AUTHOR-GATED question: log it, ask, do not guess |
| 3 | Recompute C_geom (Kostant) weights | as above |
| 4 | Assemble the formula; compare against PDG; produce the residual table | anchors (μ_Λ, m_e) declared as calibration inputs, per the ledger |
| 5 | Adversarial audit: a second agent recomputes 1-3 independently | disagreement between recomputers is itself a finding |

### Suggested definition of done

| # | Item |
| --- | --- |
| 1 | `scripts/m8_3_mass_reproducer.py` runs standalone and regenerates every number; `data/m8_3_masses.json` |
| 2 | Residual table published at the ledger's weight: the within-3× hit rate is REPORTED but explicitly NOT counted as evidence for the torsion map (the author's own pre-registered null, p = 0.174) |
| 3 | Any constant that could NOT be recomputed from a published definition is listed by name (that list is a deliverable, not a failure) |
| 4 | Method note per [`dev_docs/METHOD_NOTE.md`](../../../../../dev_docs/METHOD_NOTE.md); MODELS.md mass-cell updated honestly |

### Blindspots

| Risk | Guard |
| --- | --- |
| Quoting a printed constant "temporarily" | defeats the task's entire point; the script must derive or fail loudly |
| Over-reading the hit rate | the null-test context is restated next to every residual table |
| Definition drift between papers | pin each definition to its record id ([`../../theory/_CITATIONS.md`](../../theory/_CITATIONS.md)); mismatches between records are findings |

### Ownership + gating

Author-driven (the platform's pointer map and standards support it). Ungated: can run
before or in parallel with M8.1/M8.2.

## DEVIATIONS LOG

(none)

## FINDINGS

(pending)
