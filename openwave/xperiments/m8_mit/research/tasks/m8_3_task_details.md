# M8.3: The mass-formula reproducer script

> Roadmap row: [`../m8_roadmap.md`](../m8_roadmap.md). Status: Implementation complete; PR
> pending (2026-07-28, author-run, normal fork → PR → review workflow, no independence
> firewall). Method note: [`../findings/m8_3_method_note.md`](../findings/m8_3_method_note.md).
> This is a scaffold-stage planning aid written by the maintainers (2026-07-21); the author
> owns the column and has amended it below with the executed findings.

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
set. This is the [`ONBOARDING_MODELS.md STEP 0`](../../../../../ONBOARDING_MODELS.md)
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
| 2 | Residual table published at the ledger's weight: the within-3× hit rate is REPORTED but explicitly NOT counted as evidence for the torsion map (the author's own pre-registered null, p = 0.174 [original planning value; superseded by the corrected table's `mass-null-v1.1`, `p_A = 0.690`, see FINDINGS below]) |
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
before or in parallel with M8.1/M8.2. No independence firewall applies here (the author's
2026-07-28 note): unlike M8.5's deliverable A, this task has no target the author's own
context must be held away from, so it runs the normal fork → PR → review workflow.

## DEVIATIONS LOG

| # | Deviation | Reason |
| --- | --- | --- |
| 1 | `data/m8_3_masses.json` also carries the ranked table, the PDG scorecard, and null-test provenance, not only "every number" the formula outputs. | One artifact as the source of truth for anything downstream (MODELS.md, the method note) cites. |
| 2 | The mutation-test registry (23 gates, coverage-enforced: every gate id attacked, every mutation red) goes beyond the docstring's original placeholder. | Already a standing M8 roadmap rule ("every PASS gate must be demonstrated to fail," precedent `m7_trivial_ok`); closing it finished the task to its own written standard. |

## FINDINGS

1. **Half-integer torsion defect discovered.** Building against the then-published mass-spectrum page, this script's own spectral-zeta computation of the four half-integer torsion singles (R1, R2, R6, R8) disagreed with the page. Root cause: those values were the coexact one-form contribution only, omitting the scalar tower's `-2*zeta'_scalar(0)` term (supported at half-integer j). Full diagnosis: [method note](../findings/m8_3_method_note.md) § 3.1.

2. **Corrected result and null context.** Restoring the scalar term gives an exact closed form for every irrep, not just the four integer-spin ones the page previously had: `T^2(R1)=phi^-4/4`, `T^2(R2)=phi^4/4`, `T^2(R6)=1`, `T^2(R8)=4`. Reported upstream; mode-identity-theory corrected the page and re-ran its pre-registered null test as `mass-null-v1.1` (`p_A = 0.690`, corrected table, superseding `mass-null-v1.0`'s `p_A = 0.174` on the pre-correction table). Re-deposited at Zenodo [10.5281/zenodo.21652153](https://doi.org/10.5281/zenodo.21652153) (concept DOI 10.5281/zenodo.18603975).

3. **Current scorecard.** With m_e as the calibration benchmark (not counted), of the remaining 8 charged fermions: 5 have a compatible entry within ×3 (μ, s, t, τ, b); 4 survive sector-first adjudication (μ, s, t, τ). Down is assigned but outside ×3 (3.22); up and charm are unassigned (up's former ~6% hit was the coexact-only artifact); bottom is compatible (1.17) but outside its own structural sector (R2). Full figures: `data/m8_3_masses.json`.

## Links

- Script: [`scripts/m8_3_mass_reproducer.py`](../scripts/m8_3_mass_reproducer.py)
- Data: [`data/m8_3_masses.json`](../data/m8_3_masses.json)
- Method note: [`findings/m8_3_method_note.md`](../findings/m8_3_method_note.md)
- PR: pending from branch `m8.3-mass-reproducer`
