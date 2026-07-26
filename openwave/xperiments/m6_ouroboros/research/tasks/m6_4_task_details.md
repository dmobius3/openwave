# M6.4: The stability census rerun + the ω-selection hunt

Roadmap row: [`m6_roadmap.md`](../m6_roadmap.md) (In Progress). Consumes: [M6.1](m6_1_task_details.md) (the certified v11 convention sheet + ℒ_ref), [M6.2](m6_2_task_details.md) (the far-field localization machinery + the branch-(b) close). Feeds: canonical OQ3 + OQ7 ([`m6_theory_canonical.md`](../m6_theory_canonical.md)), the hunt μ/τ + stability rows ([`m6_particle_hunt.md`](../m6_particle_hunt.md)), the MODELS.md stability cell.

## TASK PLANNING

**Scope**: two arms, both branch-independent of the M6.2 outcome. **(a) The Gelfand-Fomin census rerun** (canonical OQ7): Zenodo 20866581 v2 claims 319 of 360 parameter families pass a Wronskian conjugate-point stability test; the corpus elsewhere says 62 (z20044392, May; internally also 42) and 62 again (z21268405, July 8). Reproduce the published scan from the record's own inline material, reconcile the three counts, state the Lean-vs-integrated ODE mismatch term-by-term, and extend the test beyond radial perturbations (the record's admitted gap). **(b) The ω-selection hunt** (canonical OQ3, open since the archive era): scan the frozen term set for ANY mechanism that makes the lepton ω values (1, 11-12.8, 40.7-50, 15) discrete; a clean negative bounds the lepton-spectrum claims and decides whether the μ/τ rows are physics or curve-labeling.

**Definition of done**:

| # | Criterion |
| --- | --- |
| 1 | The published census reproduced from the record's own ODE + test + grid (or the failure to reproduce characterized with the exact divergence point); the count reported against 319/360 |
| 2 | The Lean-formalized ODE vs the integrated ODE vs the M6.2-certified v11 reduction: a term-by-term comparison table, script-verified |
| 3 | The 62 vs 319 vs 62 reconciliation: each count attributed to its object (what was counted, under which criteria), script-backed where the sources permit |
| 4 | The GF test extended to non-radial (l ≥ 1) perturbation sectors on the radially-passing families; the census restated with the extension |
| 5 | OQ3: the localization window in ω mapped for the frozen term set; a discreteness verdict (mechanism found, or a clean negative) with the explicit μ/τ-ladder implication |
| 6 | Adversarial audit of every substantive claim; method note in `findings/`; doc checker clean; TASK REVIEW presented in the terminal |

**Pre-registered criteria (fixed BEFORE numerics, per the standing rigor bar)**:

| Item | Criterion fixed in advance |
| --- | --- |
| "Reproduced" (census) | the same 319/360 count under the record's own ODE, test, and grid as printed in the docx; if the docx under-specifies the grid, the reconstruction is documented and the count reported for the reconstruction (labeled as such), not tuned toward 319 |
| "Stable" (record's sense) | no Wronskian zero crossing before r_max, exactly the record's stated test; our extensions (localization class, l ≥ 1) are reported as SEPARATE columns, never silently merged into the record's criterion |
| "Localized" | the far-field linear system around the vacuum admits a decaying channel at the family's (ω, λ) AND the solved profile's tail follows it; states without a decaying channel are non-localized regardless of any stability count |
| "Discreteness" (OQ3) | a mechanism counts ONLY if it selects isolated ω values from a continuum via the frozen term set itself (spectral condition, topological constraint, variational selection); grid artifacts, solver tolerances, and curve-labeling do not count |
| No-search rule | no parameter is tuned toward any published number (319, 62, L/Q = 0.5, m_μ/m_e); scans run on pre-declared grids; forks resolved by derivation and reported with all numbers |

**Gating**: M6.1 ✅ (certified spec) + user "go" (given 2026-07-20 23:32 EDT); independent of the M6.2 branch by design.

**Blindspot pass** (unknown unknowns surfaced):

| Risk | Routing |
| --- | --- |
| The record's "attached script" does not exist as a file on Zenodo (verified: v1 = `chaoiton_theorem.lean.txt` only, v2 = the 13 KB docx only), so the rerun rests on the docx's inline listing; if the listing is partial, the reconstruction gap must be documented, not papered over | machine-checkable (D0 extracts the docx verbatim; canonical § 4 provenance row gets corrected either way) |
| Wronskian conjugate-point testing on a COUPLED 2-component Jacobi system needs a 2×2 fundamental-matrix determinant (4 independent solutions), not 2; if the record used the scalar recipe on the coupled system, the 319 count may be methodologically void before any physics | machine-checkable (both constructions run side by side) |
| Stability-vs-localization conflation: the far-field analysis of the Lean ODE gives decaying channels only for ω² < (λ + √(λ²+4))/2, which EXCLUDES the record's own best-fit point (ω = 2.20, λ = 0.5); the census may be counting non-localized states (the M6.2 pathology again) | machine-checkable (localization class per family) |
| r_max window dependence (the M6.2 lesson: 16% drift): census verdicts could flip with the window | machine-checkable (window-sensitivity column) |
| The 62-family counts may count DIFFERENT objects than the 360-grid scan (solution families vs parameter combos vs "distinct" post-dedup families) | provenance reading first (D1); anything not decidable from the record text is logged author-gated, never assumed |
| ω-scan existence: shooting/BVP solvers can fabricate or miss solutions near thresholds; discreteness could be a solver artifact | machine-checkable (two independent methods on any candidate discreteness signal, per the audit rule) |

**Research body destinations**: task record here; scripts `../scripts/m6_4_*.py`; data `../data/m6_4_*.json` (+ local `.npz` per the dataset policy, manifest regenerated); plots `../plots/m6_4_*.png`; method note `../findings/m6_4_method_note.md`; checkpoints `../checkpoints/m6_4_progress.md`; canonical OQ3/OQ7 + § 4 ledger rows, hunt + MODELS.md stability cells at REVIEW.

**Sub-experiments (execution order)**:

| ID | What | Output |
| --- | --- | --- |
| D0 | Extract the 20866581 v2 docx verbatim (pandoc); pin the integrated ODE, the test recipe, the 360-point grid, the claimed numbers; verify the Zenodo file inventory of both versions | reference text + `m6_4_record_extract.json` |
| D1 | Provenance of the three counts: pull the exact claim wording + criteria from z20044392 (May), z20866581 v2 (June), z21268405 (July 8) | reconciliation table (inputs) |
| N1 | Reimplement the record's census exactly as printed (ODE + K-Bessel start + Wronskian test + grid); reproduce or diverge | `m6_4_gf_census.py` + count |
| N2 | Per-family classification columns: far-field channel (localized / oscillatory), profile tail behavior, window sensitivity, coupled-system (2×2) Wronskian vs the record's recipe | census table + plots |
| N3 | Non-radial extension: the l ≥ 1 Jacobi sectors on radially-passing families | extended census |
| N4 | OQ3 ω-scan: localization window + solution-family structure vs ω for the frozen term set (the certified v11 reduction AND the record's ODE); discreteness verdict | `m6_4_omega_selection.py` + verdict |
| A | Adversarial audit (independent agent, own scripts) of every claim above | audit section in the method note |

**Model/effort**: Fable 5 / high (research default; novel numerics + audit).

## FINDINGS

Full record with equations, code map, and embedded figures: [`m6_4_method_note.md`](../findings/m6_4_method_note.md). The headline table (all ✅ script-backed unless noted):

| # | Finding | Status |
| --- | --- | --- |
| F1 | The Lean-vs-integrated mismatch is RESOLVED and fatal to the record's self-consistency: the Lean ODE of z20866581 § 1 is the May system verbatim (correct Klein-Gordon +ω²), while the same record's § 2-4 "numerical verification" integrates a structurally different system (flipped ω² inside m² = 1+ω², no centrifugal, extra α-cubic, halved coupling); no sign redefinition maps one onto the other (sympy) | ✅ |
| F2 | The record's stability test is doubly broken: the printed Jacobi Q₁₂ = −c/4 is half the true linearization (−c/2), and the value-IC Wronskian recipe is not the Gelfand-Fomin conjugate-point construction (the two disagree on ~20/360 families) | ✅ |
| F3 | The 62 vs 319 vs 62 counts are RECONCILED as: an irreproducible May count (0/1280 pass the printed criteria under faithful reconstruction, audit-confirmed on a 136-combo independent subsample; the same paper prints both 62 and 42), a June count reproduced to within 1-2 across integrators (318 primary / 317 audit vs published 319, integrator-unstable, AF-b) ONLY under the tail-multiplier convention on near-vacuum backgrounds whose localization is imposed by construction, and a July citation of the dead May count as "Hard Quantization" while the June record's own § 3.2 says "not fine-tuned ... broad regions pass" (89% of the grid = the opposite of quantization) | ✅ |
| F4 | The June census is window-defined: 94/360 stability verdicts flip between r_max = 10 and 15 (the M6.2 window pathology at census scale); best-point L/Q reproduces as 0.471 vs the published 0.567, and the published L/Q range and "340 localized" do not reproduce under any amplitude convention | ✅ |
| F5 | The non-radial gap closes WITHIN the reduction: +(l²/r²)·I is positive semidefinite, and by the matrix/Morse-index comparison for self-adjoint systems (the coupled-system form of the Sturm argument; audit A7) conjugate points are monotone non-increasing in l (5 spot checks, audit-reproduced); radial verdicts are binding; ansatz-breaking modes remain M7's domain | ✅ |
| F6 | OQ3 CLOSES NEGATIVE, strongest form: the frozen-spec charged-sector reduction admits NO localized state at any in-window ω (48-ω envelope scan, both panels, amplitudes to 1.2; the 3 near-window-top candidate zeros were chased and refuted as aliasing artifacts of nonlinear oscillatory tails: fields O(1) at r = 25 vs 3e-4 required); every ladder ω (11-50) sits above every window (ω_max ≤ 1 for λ ≥ 0); the electron's ω = 1 sits exactly ON the λ = 0 boundary; and the flipped-sign system that produced the ladder localizes at EVERY ω by construction, so no selection mechanism is possible in it either. The term set selects amplitudes, never frequencies | ✅ |
| F7 | Provenance: `chaoiton_gf_verification.py` does not exist on either Zenodo version of the record lineage (v1 = the Lean txt, v2 = the docx), while the v2 docx claims it is attached; every rerun is a documented reconstruction | ✅ |
| F8 | Adversarial audit: 9/10 claims CONFIRMED by independent methods (different integrators, exact-Bessel projections, own algebra, full-360 independent June rerun), A7 PARTIAL (conclusion holds; justification upgraded to the matrix/Morse-index form), none refuted; auditor findings AF-a..AF-g all dispositioned in the method note § 6 | ✅ |

## DEVIATIONS LOG (running, appended during EXECUTE)

| # | Deviation | Action taken |
| --- | --- | --- |
| 1 | The first June-census run ground for 8+ min inside saturated inward integrations (fields to ~1e5 with cubic curvature; DOP853 stepping microscopically) | added the divergence-classification guard at \|field\| > 100 (documented in-code: it only classifies rows, never alters surviving rows' numbers); rerun clean |
| 2 | The Nelder-Mead ω-scan (N4) reported non-existence at every in-window ω, but the optimizer-free contour probe found 31 candidate zero-crossing cells at ω = 0.9 (g=0.5, λ=0), min \|F\| = 6.5e-3: an optimizer funnel-miss, exactly at the May record's representative point | conservative option taken: added N4c (`m6_4_root_refine.py`): 2D root-finding from every candidate cell + ω-continuation of any converged root + GF test on it; the § 3.5 verdict is stated from the refined result, not the scan |

## TASK REVIEW (2026-07-21)

**Task Duration:** 01:10 (from 2026-07-20 23:32 to 2026-07-21 00:42)
**Usage Cap Triggered:** NO

**Results**

| # | Result | Status |
| --- | --- | --- |
| 1 | Lean-vs-integrated mismatch RESOLVED: the Lean statement formalizes the May system; the same record's numerics integrate a structurally different (flipped-ω²) system; no redefinition maps them | ✅ measured (sympy + audit A1) |
| 2 | The record's stability test doubly broken: Q₁₂ = −c/4 is half the true linearization; the value-IC Wronskian is not the GF construction (disagrees on ~21/360) | ✅ measured (audit A2/A6) |
| 3 | 62 vs 319 vs 62 reconciled: May 0/1280 reproduce (paper prints 62 AND 42); June 318/360 only under the tail-multiplier convention with localization imposed by construction (V and S conventions: 0/360); July cites the dead May count as "Hard Quantization" against June's own "not fine-tuned" | ✅ measured |
| 4 | The June census is window-defined (94/360 flips between r_max 10/15) and integrator-unstable (319/318/317); best-point L/Q 0.469-0.471 vs published 0.567; 285/360 L/Q negative | ✅ measured |
| 5 | Non-radial extension closed within the reduction: conjugate points monotone non-increasing in l (matrix/Morse-index + 5 audited spot checks); ansatz-breaking modes = M7's domain | ✅ measured |
| 6 | OQ3 CLOSED, clean negative: no localized charged-sector state at ANY in-window ω (48-ω scan + dense grids + BVP; 3 candidate zeros refuted as aliasing artifacts, fields O(1) at r=25 vs 3e-4 required); every ladder ω above every window (ω_max ≤ 1); electron ω = 1 exactly ON the boundary; the flipped-sign system localizes at every ω by construction and cannot select either | ✅ measured |
| 7 | Provenance: `chaoiton_gf_verification.py` on NO Zenodo version despite the docx claiming attachment | ✅ measured (audit A10) |
| 8 | Adversarial audit: 9/10 CONFIRMED, A7 PARTIAL (justification upgraded), 0 refuted; AF-a..g dispositioned in the method note § 6 | ✅ |

**Issues / blockers**: none open. Conv S ran 50 min (stiff Jacobi solves on oscillatory backgrounds); the audit began before the census JSON existed (AF-a), resolved by finalizing § 3.3 after completion.

**Deviations from plan**: two, logged in § DEVIATIONS: (1) divergence-classification guard added after the first census run ground at machine-small steps; (2) the optimizer's non-existence verdict was not trusted when the contour probe found candidate zeros: the chase (root-finding → continuation → profile diagnostic) settled them as artifacts.

**Action needed**: none pending; the approved sync is applied (canonical OQ3/OQ7 resolved + 3 § 4 ledger rows + the stability-count row corrected; hunt μ/τ + pion + stability rows + closing; MODELS.md 3 cell rewords, statuses unchanged; briefing sector rows + roadmap row + help-wanted; roadmap M6.4 → Done appended, M6 back on HOLD). M6 has no live tasks.

**Findings**: The corpus's "most reproducible artifact" is not reproducible: its three stability counts trace to two incompatible ODE systems (the record's own numerics do not implement its own Lean-formalized theorem), a broken test, and one amplitude convention whose localization is imposed rather than found. The ω-selection question (OQ3, open since the archive era) closes with a clean negative of the strongest form: the frozen charged-sector term set admits no localized state at any frequency, so the lepton ladder (μ, τ, pion) is curve-labeling on non-localized configurations of a sign-flipped system. Both closures are adversarially audited (9/10 confirmed + 1 partial, none refuted) with all audit scripts preserved.

**Research docs created / updated**:

- [`m6_4_task_details.md`](m6_4_task_details.md) (this record) · [`m6_4_method_note.md`](../findings/m6_4_method_note.md) (the full record: equations, code map, results, audit)
- scripts: [`m6_4_gf_census.py`](../scripts/m6_4_gf_census.py) · [`m6_4_may_census.py`](../scripts/m6_4_may_census.py) · [`m6_4_ode_comparison.py`](../scripts/m6_4_ode_comparison.py) · [`m6_4_nonradial_extension.py`](../scripts/m6_4_nonradial_extension.py) · [`m6_4_omega_selection.py`](../scripts/m6_4_omega_selection.py) · [`m6_4_existence_probe.py`](../scripts/m6_4_existence_probe.py) · [`m6_4_root_refine.py`](../scripts/m6_4_root_refine.py) · [`m6_4_root_diagnostic.py`](../scripts/m6_4_root_diagnostic.py) + 6 preserved audit scripts (`m6_4_audit_*.py`)
- data: `m6_4_*.json` (14 tracked) + record fulltext extracts + [`_DATASETS.md`](../data/_DATASETS.md) (census/may row npz local-only)
- plots: [`m6_4_omega_window.png`](../plots/m6_4_omega_window.png) · [`m6_4_existence_contours.png`](../plots/m6_4_existence_contours.png) · [`m6_4_root_profiles.png`](../plots/m6_4_root_profiles.png) · [`m6_4_census_best_point.png`](../plots/m6_4_census_best_point.png) · [`m6_4_may_representative.png`](../plots/m6_4_may_representative.png) · [`m6_4_branch_continuation.png`](../plots/m6_4_branch_continuation.png)
- synced: [`m6_theory_canonical.md`](../m6_theory_canonical.md) (OQ3/OQ7 + § 4) · [`m6_particle_hunt.md`](../m6_particle_hunt.md) · [`../../__M6_model_briefing.md`](../../__M6_model_briefing.md) · repo-root `MODELS.md` · [`m6_roadmap.md`](../m6_roadmap.md)
