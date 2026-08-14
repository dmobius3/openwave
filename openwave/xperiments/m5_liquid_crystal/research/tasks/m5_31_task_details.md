# M5.31 — The M5/Faber curvature scale curve: dense C(rho) instrument (external contribution)

**Contributor:** [@vantasnerdan](https://github.com/vantasnerdan) (first external delivery on the running-coupling debt path). **Delivered via:** [PR #437](https://github.com/openwave-labs/openwave/pull/437) (rewritten from its first submission after a request-changes review), preregistered in [Discussion #438](https://github.com/openwave-labs/openwave/discussions/438).

## What the task delivered

A dense, audited measurement of the dimensionless shell observable `C(rho) = <r^2 sqrt(sum |R_ij|^2)>` on the M5.6.4b regularized hedgehog, 25 logarithmic radii over `0.6 <= rho <= 5.5`, replacing the five-shell onset statement of [`m5_6_4b_faber_curvature_em.py`](../scripts/m5_6_4b_faber_curvature_em.py) with a resolved scale curve. Full record with equations, gates, and results: [the method note](../findings/m5_31_coupling_curvature_note.md).

| Artifact | File |
| --- | --- |
| Physics module | [`m5_31_coupling_curvature_field.py`](../scripts/m5_31_coupling_curvature_field.py) |
| Scan driver (11 gates) | [`m5_31_coupling_curvature_scan.py`](../scripts/m5_31_coupling_curvature_scan.py) |
| Independent audit (16 gates, no contribution imports) | [`m5_31_coupling_curvature_audit.py`](../scripts/m5_31_coupling_curvature_audit.py) |
| Tracked evidence | [`m5_31_coupling_curvature_scan.json`](../data/m5_31_coupling_curvature_scan.json), [`m5_31_coupling_curvature_audit.json`](../data/m5_31_coupling_curvature_audit.json), [plot](../plots/m5_31_coupling_curvature_scan.png) |

## Scope boundary (the honest part)

The result is a classical single-ansatz form factor, not a renormalized coupling: the contribution's own independent audit REFUTES the stronger reading. Both `C -> g_R` conventions (energy/action and field-amplitude) are reported without selection; choosing one is author-gated on the M5 source/action dictionary (asked of the author in Discussion #438 and the PR body). `MODELS.md` is untouched; the Running-coupling debt row stays ⚠️ until the convention is licensed or a two-core scheme measures the coupling directly.

## Review record (2026-08-14)

Maintainer verification, run independently before merge:

| Check | Result |
| --- | --- |
| Closed forms (curvature norm, log-slope) | ✅ re-derived in an independent sympy script, exact zero residual both |
| Gate reproduction | ✅ all 11 scan + 16 audit gates reproduce; shipped JSON matches regeneration to machine epsilon |
| Mutation discrimination | ✅ corrupting one stored value fails 2 audit gates; a 5% stencil error fails 3 scan gates |
| Audit independence | ✅ imports neither contribution module (numpy/sympy/scipy only) |
| Code intent (A9) | ✅ no network, no environment reads, writes only under `research/data` and `research/plots` |

Maintainer edits at merge (announced in the review): task-ID rename `m5_coupling_curvature_*` to `m5_31_*`, fork-branch permalinks repointed to repo-relative paths, this task record and the roadmap row added. Review under [`PR_REVIEW_STANDARDS.md`](../../../../../dev_docs/PR_REVIEW_STANDARDS.md) § 12.1: no blocking findings; approve-with-comments.

## Successors (contributor-proposed, not yet allocated)

Discussion #438 preregisters two further measurement tasks: a confinement/string-tension instrument and a boosted-clock two-body force test. Each gets its own task ID when its PR materializes.
