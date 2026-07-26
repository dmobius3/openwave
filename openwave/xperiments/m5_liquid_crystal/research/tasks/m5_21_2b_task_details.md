# M5.21.2b: the well-posed 3D instrument (term-set discrimination + converged census)

**Status**: 🔶 IN PROGRESS (go 2026-07-17 20:59 EDT, reset 1:30am, resume ping armed; roadmap row in [In Progress](../m5_roadmap.md)).

**2026-07-17 14:16 update**: Duda's reply to the checkpoint package ([`m5_21_convo.md § 2026-07-17 14:16`](m5_21_convo.md)) DOUBLY confirms this task as the next run and adds two headline observables, his asks back to us: [Q30](../m5_question_tracker.md#q30-detail) (the ring radius a(δ) + extrapolation toward δ ~ 1e-10 + the ~2fm scattering-bound comparison) and [Q31](../m5_question_tracker.md#q31-detail) (the split-line count in BOTH hemispheres, four-vs-two, the inter-line angles, the δ-trend). No regularization or term-set preference was given, so the I1 order below stands self-determined. The 1+1D toy lineage is formally retired on his side; 3D minima → 4D seeds reaffirmed.

## TASK PLANNING

**Scope**: build the discretization-consistent 3D instrument that the M5.21.2 instrument finding demands (no stencil-consistent minimizer at toy parameters on the bare quartic functional), discriminate the Q25 term sets on it, then produce the CONVERGED census and the answers to his two 14:16 asks (Q30 ring radius, Q31 split geometry). **Standing user directive (at go): every instrument-level verdict is read on BOTH defect types, the point hedgehog AND the charged ring**: I1 consistency gates, I2 term-set runs, I3 census, and films all run the pair.

**Definition of done** (testable):

| # | Criterion |
| --- | --- |
| 1 | I1: an instrument whose deep endpoints pass the consistency gate: E_u re-read under fwd / bwd / 2h stencils + the factor-2 subsample probe agree within ×1.5, AND the fixed-physical-box refinement ladder (N ∈ {32, 48, 64}, L = 48) shows a shrinking inter-rung ΔE; measured for seed A AND ring R. If NO candidate passes, the honest negative with the per-candidate failure mode (itself a Q25 answer) |
| 2 | I2: T1 (trace-target) vs T2 (eigenvalue penalty) vs T3 (shifted-spectrum det penalty) each relaxed on the winning instrument at the working point, A + R, with the virial accounting (E_u vs 3E_V + the ε-term) and a discrimination verdict |
| 3 | I3: converged census A/B/C + R at N = 48 on the winner: energies, the ordering vs the qualitative A < C < B, the ring-vs-point verdict |
| 4 | I4 (Q30): relaxed ring cord radius across seed a ∈ {3, 4, 6} × δ ∈ {0.3, 0.2, 0.1}: the a(δ) read + trend statement (extrapolation qualitative; the voxel → fm scale caveat named, Q17) |
| 5 | I5 (Q31): half-line count in BOTH hemispheres + inter-line angles on relaxed A endpoints at δ ∈ {0.3, 0.2, 0.1}: the "how many, of what degree, in what angles" answer + the δ-trend |
| 6 | Films: final A + R endpoints on the TRUE `m5_film.py` basic + thermal templates via the [adapter](../scripts/m5_21_2_e_film_adapter.py) |
| 7 | Independent adversarial audit of the substantive claims BEFORE the review; doc checker clean over every touched `.md` |

**Gating**: [M5.21.2](m5_21_2_task_details.md) ✅ (met). Runs ahead of [M5.21.3](m5_21_3_task_details.md), which consumes THIS task's converged minima.

**Blindspot pass** (unknown unknowns surfaced at PLAN):

| Risk | Fold into scope |
| --- | --- |
| The ε-Dirichlet term changes the Derrick balance (quadratic gradient scales +λ: quartic + quadratic = Skyrme-type stabilization) | Track the FULL virial identity per ε rung (E_D − E_u + 3E_V, sign convention checked at gate time); never assume the bare u = 3V form under ε > 0 |
| T2 eigenvalue-penalty gradient needs eigenvector derivatives; near-degenerate eigenvalues at cores make it noisy | Gap-floor guard + monitor the min inter-eigenvalue gap along descent; report if the guard ever activates |
| Fixed-physical-box refinement changes the pin-shell thickness in voxels | Define the pin width in PHYSICAL units (else the ladder confounds pin depth with resolution) |
| T3 det gradient via M⁻ᵀ is singular when an eigenvalue crosses 0 mid-descent | Use the adjugate (cofactor) form, finite always |
| Cross-stencil re-reads at mixed h conventions corrupt the ×-ratios | h explicit in every difference and every h³ cell weight; gates re-run on the h-aware stack before physics |
| 64³ deep runs are the slow wall-clock tail | Ladder rungs at reduced depth (4k, plateau-checked); full depth (8k+) only at N = 48 |

**Research body**: findings → `research/findings/m5_21_2b_note.md` (method-note grade, the next-package carrier); scripts `research/scripts/m5_21_2b_*.py`; data `research/data/m5_21_2b_*`; plots `research/plots/m5_21_2b_*`; checkpoints → `research/checkpoints/m5_21_2b_progress.md` (gitignored).

**Preconditions**: the M5.21.2 gate battery re-run on the new h-aware stack (complex-step gradient check, SO(3) conjugation invariance, vacuum zero) BEFORE any physics run; per-run row files (race-free) + eager checkpoints per the series standard.

**Model/effort**: Fable 5 / high (research default; set at go).

## FINDINGS (2026-07-17 night run)

Full record with equations, code map, gates, and every table: [`findings/m5_21_2b_note.md`](../findings/m5_21_2b_note.md). The headline results:

| # | Finding | Status |
| --- | --- | --- |
| 1 | **THE INSTRUMENT FIX**: the stencil-symmetrized functional E = ½(E_fwd + E_bwd) removes the M5.21.2 pathology outright: true stationary points (f_tol) with cross-stencil ratio 1.10-1.14 (the fwd control reproduces the funnel on demand: ×294-386); the ε-Dirichlet family is linear in ε and lands on the bare ε = 0 endpoint (the well-posedness certificate); Q25's well-posedness reframe is RESOLVED constructively | ✅ measured, audited |
| 2 | **THE Q25 TERM-SET DISCRIMINATION**: T2 (his Eq-12 eigenvalue-penalty form) is the only term set whose bare (ε = 0) minimum is simultaneously consistent AND virial-balanced (resid +0.034 at N = 48, retention 0.89, compact, eigen-gap held 0.012-0.038); T1 needs ε* ≈ 2e-3 to balance; T3 (shifted det) never lands (resid −0.58, gap-collapsed cores) | ✅ measured |
| 3 | **THE CENSUS on the clean instrument**: A < C < B survives with B and C now CERTIFIED stationary (f_tol) = three distinct protected minima; ratios C/A ≈ 4.2, B/A ≈ 16.0; A/R pass consistency (xr 1.14), B/C lattice-flagged (xr 2.4-2.6); refinement ladder N = 32/48/64: consistency holds every rung, absolute E drifts +7-13% (quote orderings, not values) | ✅ measured |
| 4 | **THE BASIN MERGER (Q29/Q30 sharpened)**: point-hedgehog and charged-ring seeds relax to ONE state at every rung tested (ΔE ≤ 0.04%, field distance ≤ 6%; N = 32, 48; T1 and T2): the M5.21.2 "tie" was two seedings of the same object | ✅ measured |
| 5 | **THE Q30 ANSWER**: the core's transverse scale is DYNAMICALLY SELECTED: seeded cord radii {3, 4.5, 6} and the point seed all anneal to one δ-dependent scale (2.5 / 2.9 / 3.4 at δ = 0.3 / 0.2 / 0.1; ±1.5% at δ = 0.1), growing as δ falls (~δ^−0.2±0.1); no compact torus cord survives; fm-scale comparison awaits the Q17 anchor | ✅ measured (δ = 0.3 large-cord runs still annealing: trend-supported there) |
| 6 | **THE Q31 ANSWER (split geometry)**: TWO +½ vortex lines, both hemispheres, net +2 half-units on every cross-section (contour-verified), braiding along z (near-core pair azimuth rotates 72° → 169° → 135°), transverse scale = finding 5's; the four-line picture is NOT observed at δ = 0.3; measured band-clean on T2 (T1 interiors are gap-collapsed tangles: stated) | ✅ measured |
| 7 | **The split-reader hardening** (anti-recipe): plaquette winding fails AT the axis (exact π/2 steps = the mod-π branch boundary; the M5.21.2d float32-parity read was masking this) and across band-crossing annuli (fake opposite-sign cores); the contour winding with per-contour gap flags is the valid instrument (seed +2 everywhere, gap 0.28) | ✅ measured |
| 8 | Films: every N = 48 final on BOTH TRUE templates (basic + thermal); the film-rule catch from the M5.21.2 review is closed | ✅ |

## DEVIATIONS LOG

| When | Deviation |
| --- | --- |
| 22:0x | The planned "winner selection" collapsed early: sym passed at ε = 0 (the plan expected ε to be load-bearing); the ε ladder was kept as the certificate rather than the fix |
| 22:1x | I2 T2's virial landing promoted T2 from a discrimination arm to a co-instrument: a T2 N = 48 pair was added to the census fleet (2 unplanned runs) |
| 23:0x | The I5 plaquette reader (planned as a straight extension of M5.21.2d) was found BROKEN on deep endpoints (axis ambiguity + band crossings): replaced mid-task by the contour-winding instrument; the M5.21.2d two-line result survives but its method is retro-flagged (see finding 7) |
| 23:2x | I4's δ = 0.3 large-cord runs (a = 4.5, 6) did not finish annealing at 8k iters: the collapse statement there rests on the trend + δ = 0.2/0.1 completions (honest caveat in § 7 of the note) |

## ARTIFACTS + CLEANUP

Scripts: [`m5_21_2b_a_instrument.py`](../scripts/m5_21_2b_a_instrument.py) (stack + gates + relax), [`m5_21_2b_b_split.py`](../scripts/m5_21_2b_b_split.py) (split census, contour-hardened), [`m5_21_2b_c_films.py`](../scripts/m5_21_2b_c_films.py) (TRUE-template films), [`m5_21_2b_d_panel.py`](../scripts/m5_21_2b_d_panel.py) (summary panel), [`m5_21_2b_audit_check.py`](../scripts/m5_21_2b_audit_check.py) (the independent audit's own checks). Data kept: `m5_21_2b_all.json` (every run row), `m5_21_2b_gates.json`, `m5_21_2b_calib.json`, `m5_21_2b_split_*.json`, `m5_21_2b_audit.json`, per-run row JSONs, N = 32 endpoint npz (< 1 MB each). Plots: panel + 10 films + 2 split figures.

Deleted raw data (> 1 MB rule; ALL regenerate deterministically via `cd research/scripts && python3 m5_21_2b_a_instrument.py relax <keys>` with the keys recorded in each row of `m5_21_2b_all.json`; exact commands below):

| File(s) | Size | Regen |
| --- | --- | --- |
| `m5_21_2b_end_lad64_{A,R}.npz` | 5.7 MB each | `relax seed=<A\|R> term=T1 stencil=sym eps=0 n=64 maxit=6000 tag=lad64_<S>` |
| `m5_21_2b_end_c48_{A,B,C,R}_d0.3.npz` | 2.4 MB each | `relax seed=<S> term=T1 stencil=sym eps=0 n=48 maxit=12000 tag=c48_<S>_d0.3` |
| `m5_21_2b_end_c48_A_d{0.2,0.1}.npz` | 2.4 MB each | `relax seed=A term=T1 stencil=sym eps=0 n=48 delta=<d> maxit=12000 tag=c48_A_d<d>` |
| `m5_21_2b_end_c48_{A,R}_T2.npz` | 2.4 MB each | `relax seed=<S> term=T2 stencil=sym eps=0 n=48 maxit=12000 w2=0.002758100 tag=c48_<S>_T2` |
| `m5_21_2b_end_i4_a{3,4.5,6}_d{0.3,0.2,0.1}.npz` (9) | 2.4 MB each | `relax seed=R term=T1 stencil=sym eps=0 n=48 delta=<d> aring=<a> maxit=8000 tag=i4_a<a>_d<d>` |

## Scope (stub level)

M5.21.2 ended instrument-blocked: at toy parameters the quartic commutator functional has NO stencil-consistent minimizer (each stencil's descent hides curvature in its blind directions, ×7-128 cross-stencil disagreement, audited). This task builds the instrument that CAN converge a 3D electron, then re-runs the census on it for converged minima.

| Arm | What | Kill / survive read |
| --- | --- | --- |
| I1 the instrument | Candidate fixes, cheapest first: (a) stencil-symmetrized functional (average fwd/bwd stencil energies: kills both known blind families); (b) a tiny Dirichlet regularizer ε‖∇M‖² with ε → 0 extrapolation (kills the aligned-gradient soft directions); (c) grid-refinement ladder (h-halving at fixed physical box) as the convergence CHECK on whichever wins | an instrument whose deep endpoints read consistently across stencils AND refinement rungs; else the honest negative that the term set itself is deficient |
| I2 term-set discrimination (the Q25 arms, self-determined) | On the winning instrument: the trace-target vs the Eq 12 eigenvalue-penalty potential; the det-constraint variant with the author's own caveat folded (spectrum must shift off zero first, his 2026-07-17 reply) | which term set pins a compact virial-balanced minimum (u = 3V landed, not just bracketed) |
| I3 the converged census | Seeds A/B/C + the charged ring R on the winning (instrument, term set): converged minima, the real lepton-hierarchy read (three minima? energy ratios?) and the ring-vs-point verdict at conviction. He endorsed the triplet reading ("clear candidates for 3 leptons"), so the converged re-verdict carries the quantitative burden | the M5.21.2 qualitative ordering (A < C < B, electron lowest) either survives to convergence or the honest re-verdict |
| I4 the ring radius ([Q30](../m5_question_tracker.md#q30-detail)) | On the converged instrument: relax the ring across the δ-ladder with FREED radius (multiple seeded a, watch where the cord settles), fit a(δ), extrapolate toward δ ~ 1e-10, state the ~2fm scattering-bound comparison (voxel → fm needs the Q17 anchors; say so) | a(δ) trend with extrapolation, or the honest null (radius pinned by the box / instrument) |
| I5 the split geometry ([Q31](../m5_question_tracker.md#q31-detail)) | On the converged instrument: half-line count in BOTH hemispheres (four-vs-two, his picture is FOUR from the point core), inter-line angles (azimuthal + polar), (count, angles) vs the δ-ladder with the extrapolation statement; extend `m5_21_2_d_vortex_split.py` to z < 0 planes + angle extraction | the "how many, of what degree, in what angles" answer at conviction (the M5.21.2 read was two ½-lines on z > 0 only, cell-scale separation) |
| Films | THE FILM-TEMPLATE FIX (user catch at the M5.21.2 review): the adapter EXISTS and rendered its first compliant strip ([`m5_21_2_e_film_adapter.py`](../scripts/m5_21_2_e_film_adapter.py), y = 0 half-plane slice, 3×3 embedded 4×4 block-diag with display g, true 3D density via the density_fn hook); regenerate ALL census films on it here | series film rule compliant |

Consolidation goal (updated 2026-07-17 14:16): the M5.21.2 package went out and LANDED ("looks great!"); the NEXT package = this task's converged answers to his two asks back (Q30 ring radius, Q31 split geometry) plus the converged census and whichever term set proves well-posed (the Q25 constructive answer).

**Gated by**: M5.21.2 ✅ + user "go" (runs ahead of [M5.21.3](m5_21_3_task_details.md), which re-gates on THIS task's converged minima). His 14:16 reply confirms the order; "I will read/think and respond further" may still add guidance, but nothing in this task waits on it.

## TASK REVIEW (2026-07-17)

**Task Duration:** 2:31 (from 20:59 to 23:30)
**Usage Cap Triggered:** NO

**Results** (all independently audited, 8/8 CONFIRMED, [`§ 11 of the note`](../findings/m5_21_2b_note.md)):

| # | Result | Status |
| --- | --- | --- |
| 1 | The instrument fix: sym stencil → genuine stationary points, xratio 1.10-1.14 (fwd control reproduces the funnel ×294/×386); ε-family linear, lands on ε = 0 (well-posedness certificate); Q25's well-posedness reframe RESOLVED constructively | ✅ measured, audited |
| 2 | Q25 term sets: T2 (Eq-12 eigenvalue penalty) = the only bare virial-balanced consistent minimum (resid +0.034, retention 0.89); T1 needs ε* ≈ 2e-3; T3 never lands | ✅ measured |
| 3 | Census: A < C < B with B/C CERTIFIED stationary (C/A ≈ 4.2, B/A ≈ 16.0); A/R pass consistency, B/C lattice-flagged; ladder consistency holds every rung, absolute E drifts +7-13% (orderings quotable, values not) | ✅ measured |
| 4 | Basin merger: point and ring → ONE state everywhere tested (ΔE ≤ 0.04%, distance 6.0-7.7%) | ✅ measured |
| 5 | Q30: transverse core scale DYNAMICALLY SELECTED (2.5/2.9/3.4 at δ 0.3/0.2/0.1, ±1.5% at δ = 0.1, ~δ^−0.2±0.1); no compact cord survives; fm comparison waits on Q17 | ✅ measured (δ = 0.3 large cords trend-supported) |
| 6 | Q31: TWO +½ lines (not four), braiding (pair azimuth 72° → 169° → 135° along z), net +2 every plane both hemispheres; band-clean on T2 | ✅ measured |
| 7 | Split-reader hardening: plaquette method retro-flagged (axis + band-crossing artifacts); contour winding + gap flags = the valid instrument | ✅ anti-recipe recorded |
| 8 | Films on both TRUE templates; independent audit 8/8 with 3 wording tightenings adopted | ✅ |

**Issues / blockers**: absolute energies not continuum-converged (ladder drift +7-13%); B/C core energies lattice-contaminated; δ → 1e-10 extrapolations qualitative pending the Q17 anchor.

**Deviations from plan**: four, logged above (ε became certificate not fix; T2 promoted to co-instrument; split reader replaced mid-task; δ = 0.3 large-cord anneals incomplete).

**Review issue (user catch at approval): the 7-shot film standard.** The deep runs persisted only endpoint fields, so the first film set had 2 rows (seed + end). Fix applied at close: `relax` now stores snapshot fields (`M_it500/1000/2000/4000/8000`) in the endpoint npz when `snaps=1`, the film script renders every stored frame (7 shots: 0, 500, 1000, 2000, 4000, 8000, 12000), and the four census finals were deterministically re-run with snapshots to regenerate all films compliant (same tags, byte-identical descent). The regenerated npz (with snapshots, ~8-17 MB each) are deleted after film rendering per the > 1 MB rule; regen = the same commands + `snaps=1`.

**Action taken at close**: roadmap row → Done (appended at end); model-doc sweep (theory canonical + briefing + particle hunt); tracker Q30/Q31 → answered-by-measurement; the Duda package drafted (terminal-only) for the user's send.

**Findings**: The 3D electron statics is now a well-posed measurement: the symmetrized stencil turns the M5.21.2 no-minimizer negative into audited stationary minima, and the author's own Eq-12 eigenvalue potential is the term set that makes the electron scale-stationary without any regularizer. On that instrument his two questions get concrete answers: the charged ring and the point hedgehog are one object (the basins merge), whose core is a braided pair of +½ vortex lines with a dynamically selected transverse scale that grows as δ falls; the four-line picture is not observed at toy δ.

**Research docs created / updated**: [`findings/m5_21_2b_note.md`](../findings/m5_21_2b_note.md) (the audited record) · this task_details · [`m5_question_tracker.md`](../m5_question_tracker.md) (Q30/Q31) · scripts `m5_21_2b_{a_instrument, b_split, c_films, d_panel, audit_check}.py` · data `m5_21_2b_{all, gates, calib, split_*, audit}.json` · plots [`m5_21_2b_panel.png`](../plots/m5_21_2b_panel.png), `m5_21_2b_split_c48_{A,R}_T2.png`, the film set.
