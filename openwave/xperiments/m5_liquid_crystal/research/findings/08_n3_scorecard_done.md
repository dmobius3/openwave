# N3 checkpoint 08 , S4 scorecard + N3 COMPLETE

`m5_11_n3_scorecard.py` (summary `m5_11_n3_scorecard.json`, figure `m5_11_n3_summary.png`). N3 is complete.

## Scorecard vs NuFIT 6.0 (NO)

| Parameter | TBM gate (predicted) | theta13 fit (delta=0.1, eps*=0.639) | NuFIT 6.0 |
| --- | --- | --- | --- |
| theta12 | 35.26 | 34.93 | 33.68 |
| theta23 | 45.00 | 40.59 | 43.30 |
| theta13 | 0.00 | 8.56 | 8.56 |
| delta_CP | 180 | 180 | 212 |

mass spectrum ratios at the gate: 1.00 : 1.15 : 1.68.

## N3 status , all DoD met

| DoD | Status |
| --- | --- |
| ★ TBM gate (3 angles from dynamics) | ✅ PASS robust 81/81 |
| theta13 derived / honest report | ✅ resolved: O(0.1) mu-tau asymmetry, NOT bare delta~1e-10 |
| search + document (geometry, delta, asymmetry) | ✅ alpha* locus + theta13(delta,eps) map |
| compare to NuFIT 6.0 | ✅ scorecard |
| reproducible | ✅ findings + scripts |

## FINISH bookkeeping
- findings doc: `n3_findings.md` (cross-linked <-> foundation findings <-> 10a <-> #199) ✅
- raw data > 1MB: NONE produced (all artifacts < 100 KB) ✅ nothing to delete
- GitHub #236: HELD (stays In progress; nothing posted; The maintainer: don't touch issues) ✅
- git: NOT committed (the maintainer's, later) ✅
- resume ping: never armed (the maintainer: 11pm reset = cap release, not deadline; one completion ping) ✅

## The 7 N3 scripts (sandbox_v10)
m5_11_n3_derisk.py (scaffold) -> m5_11_n3_mass_matrix.py (loop->matrix) -> m5_11_n3_search.py (TBM gate) ->
m5_11_n3_theta13.py (crux) -> m5_11_n3_scorecard.py (NuFIT + figure). All headless numpy f64, 16-core where parallel.
