# M5.22.4: dynamical baryons, the ω-twist probe

**Status**: 🚧 GO-READY. Staged 2026-08-05 from the author's 2026-08-02 checkpoint reply ([`m5_22_convo.md § 2026-08-02`](m5_22_convo.md)): the dynamical-minima directive. REFRESHED 2026-08-06 from the author's M5.22.2 checkpoint reply ([`m5_22_convo.md § 2026-08-06`](m5_22_convo.md)): the dynamics route re-endorsed ("mainly twists of long axis"), and the full-F electric instrument staged as the opening add-on. The series goals, borders, and checkpoint policy live in [M5.22 § TASK PLANNING](m5_22_task_details.md).

## Scope: the author's directive, run on the baryon states

The directive (2026-08-02, verbatim core): the 3×3 case "misses looking crucial angular momentums", such local energy minima "should not be static, but very dynamical", suspected "mainly by twists of the long axis"; the full route is the complete 4×4 case, "or maybe some its approximations - e.g. radius dependent hedgehog of boosts and twists of constant frequency omega".

| Piece | Content |
| --- | --- |
| OPENING ADD-ON: the full-F electric instrument | The author's 2026-08-06 definition ([convo § 2026-08-06](m5_22_convo.md)): the full electric read is the spatial F-tensor components E = (F_23, F_31, F_12), containing contributions from ALL eigenvalues and eigenvectors, vs the calibrated longest-axis dual-curvature read (the [M5.22.2](m5_22_2_task_details.md) instrument, the author's "basic" form). Build the literal full-F object and diff the two on the M5.22.2 calibration set + the four target states. Expected near-identity here (the M5.22.2 axis fork measured middle/short flux ≤ 0.34/0.01); the diff converts that expectation into a measured statement and closes the author's definition exactly. Minutes of compute |
| Machinery | The [M5.21.3](m5_21_3_task_details.md) omega-ladder (twists of constant frequency ω on a re-relaxed state, E\* = min_M [E_stat + ω² kin], the exact object the author names) + the [M5.21.9](m5_21_9_task_details.md) fixed-J clock as the dynamical-state instrument |
| Targets | The census proton-analog and neutron-analog ([M5.22](m5_22_task_details.md)), the ring-antiring neutral state and the deuteron candidate ([M5.22.1](m5_22_1_task_details.md)) |
| Reads | E\*(ω) ladder per state; whether any state's minimum sits at ω > 0 (a measured DYNAMICAL minimum, the directive's claim); the effect on the mass ordering, the n/p ratio, and the identity reads; the twist-axis choice reported (the long axis per the directive, others as controls) |
| Caution discharged | Until this runs, every static M5.22.x energy/identity read carries the author's 2026-08-02 caution that neglected angular momenta "could qualitatively change behavior" |

## TASK PLANNING (2026-08-06, at go 10:51 EDT)

**Model/effort**: Fable / high (per the series plan). **The four target states** (all existing endpoints, no re-relax of the 3D record): census proton-analog `P-0.5_plane_sc6_n32_pinned_d0.3` and neutron-analog (the ring-antiring pair) `P-1_plane_sc6_n32_pinned_d0.3` ([M5.22](m5_22_task_details.md)); the second neutral basin (pp-control cousin) `d2_s-0.5_s-0.5_a2_sc6_n32_d0.3` and the deuteron candidate `dn_-0.5at-2_+0.5at+0_-0.5at+2_n32_d0.3` ([M5.22.1](m5_22_1_task_details.md)). The scope table's "ring-antiring neutral state" resolves to the d2 basin: the census neutron-analog IS the primary ring-antiring pair, and the [M5.22.2](m5_22_2_task_details.md) decay probe treated exactly these two as "both neutral states".

| Phase | Content | Verdict artifact |
| --- | --- | --- |
| A0 the full-F add-on | Build the literal paper object (arXiv 2108.07896 eq 3-8): the oriented full eigenframe O, the connection Γ_μ = O^T ∂_μ O, Γ⃗_μ = ((Γ)_32, (Γ)_13, (Γ)_21), the curvature R⃗_ij = Γ⃗_i × Γ⃗_j, and E = (F_23, F_31, F_12) from the space-space pairs. Key derivation made at build (to verify numerically, not assume): the internal long-axis component of R⃗_ij is IDENTICALLY the basic instrument (e₃·(∂ᵢe₃ × ∂ⱼe₃), the M5.22.2 dual-curvature read), and the internal short/middle components are the Mermin-Ho densities of the other two frame axes, i.e. exactly what the M5.22.2 axis fork measured (flux ≤ 0.34/0.01). Diff protocol: pointwise + flux-ladder diff of the internal-3 read vs the calibrated `mermin_B` on the analytic hedgehog, the two calibration states, and the four targets; per-component flux table; a signed-norm combined read as the labeled extra arm | `m5_22_4_fullf.json` |
| B1 static 4D lift | Per state: `embed34` → static FIRE relax (the audited M5.21.3 V4 stack, s = +1 primary, maxit 6000; s = −1 spot-check on the neutron-analog); RING-SURVIVAL read at the endpoint (census charge profile + core ledger) since the 4D potential is the trace-target, not the census T2 (the known M5.21.3 caveat, now with structured states at risk: measured, not assumed) | per-state P1 rows |
| B2 the generator catalog | Per state: kin sign for all 6 named generators (`clock_local` = the directive's long-axis twist, `plane_1d`, global rots, boosts) + the axial-twist channel E(k) linear term b (a nonzero b = SPONTANEOUS twist selection, the read the electron nulled in [M5.21.3](m5_21_3_task_details.md) finding 7) | per-state P2 rows |
| B3 the omega-ladder | Per state: E\*(ω) = min_M [E_stat + ω² kin] on `clock_local` (ω = 0.05 → 0.8, the M5.21.3 rungs) + the matched-depth static control (15000 it); verdict per state = a measured minimum at ω > 0, or its absence, read against the control | per-state P3 rows + ctrl |
| B4 close | Panels, findings note, independent adversarial audit, method-note-grade record, tracker/roadmap/canonical/briefing sweep | note + review |

**Try caps**: 3 per machine-checkable gate. **Parameters**: n = 32, L = 48, δ = 0.3, g = 8, sym stencil, both-reading labeling as in M5.21.3 (the η-signed functional is the ladder object; the Hamiltonian reading bounds E\*(ω) ≥ E\*(0) for rotations analytically, so the ladder's information is in the re-relaxation and the twist linear term). **Blindspots logged at plan**: (i) the potential swap (census T2 → 4D trace-target) can deform ring states: the B1 ring-survival read is the guard, a lost ring is a first-class reported outcome; (ii) toy ω is not a physical rate (Q33 bridge only); (iii) the d2/“ring-antiring” naming ambiguity resolved above, flagged for the review.

## Definition of done

The full-F electric instrument built and diffed against the longest-axis read (the opening add-on); the ω-ladder measured on the four baryon-sector states with a static-vs-dynamical verdict per state (a minimum at ω > 0 measured, or its absence characterized honestly); full closeout per the series policy (method note with inline images, tracker update, roadmap move).

## Gated by

[M5.22.2](m5_22_2_task_details.md) close + user "go". User pick may pull it ahead of the conditional [M5.22.3](m5_22_3_task_details.md) stretch: the directive says the angular momenta could change the answers the stretch would build on.

## FINDINGS (2026-08-06; the record is [`../findings/m5_22_4_note.md`](../findings/m5_22_4_note.md))

| Deliverable | Result |
| --- | --- |
| The full-F instrument (opening add-on) | ✅ built as the literal paper object (eq 5-8: full oriented eigenframe, Γ_μ = O^T ∂_μ O, R⃗ = Γ⃗ × Γ⃗, E from the space-space pairs) with the DERIVED identity verified: the internal long-axis component IS the basic instrument, identically; the flux diff on the four targets is ~1e-3 (box-convention-pinned), the middle/short contributions non-quantized (≤ 0.30/0.014), the norm read not quantized. The author's full form changes no charge read ([note § 1](../findings/m5_22_4_note.md)); his independent same-day "rather negligible" assessment now carries a measured backing |
| The static 4D lift (5 runs) | ✅ all four states + the s = −1 spot-check survive the trace-target lift: charges exact (−1.029/−0.002/+0.002/−1.002), rings preserved and deepened, block-diagonality exact ([note § 2](../findings/m5_22_4_note.md)) |
| The generator catalog | ✅ the M5.21.3 electron sign table REPRODUCES on the baryons: kin > 0 for every rotation incl. the directive's `clock_local` twist (0.154 to 0.773), boosts the only negative channel; the spontaneous-twist linear term is numerically zero on all five runs ([note § 3](../findings/m5_22_4_note.md)) |
| The omega-ladder, 4 states | ✅ the honest negative: NO minimum at ω > 0 on any state; all 20 rungs sit above the matched-depth static control (+0.095 to +0.498); the decoupling exact (static parts = control to 3-4 decimals at ω = 0.8, offset = ω²·kin); charges preserved along every rung ([note § 4](../findings/m5_22_4_note.md)) |
| Adversarial audit | ✅ 5 PASS, 0 refuted (own script, own implementations); all four catches adopted into the note ([note § 5](../findings/m5_22_4_note.md)) |

![the full-F panel](../plots/m5_22_4_fullf.png)

![the omega ladder](../plots/m5_22_4_ladder.png)

Datasets (heavy arrays LOCAL-ONLY, gitignored, kept; tracked = JSONs/plots/scripts + the `_DATASETS.md` manifest): 9 endpoint npz `data/m5_22_4_p1_{prot,neut,d2,deut}[_s-1].npz` + `m5_22_4_p3_{key}_clock_local.npz` (~870 KB each; regen: `python3 m5_22_4_b_omega.py p1 <key> [s]` ~13 min each, `p3 <key>` ~45 min each, Apple Silicon parallel-4); rows `data/m5_22_4_row_*.json` + `m5_22_4_all.json` (regen: the same modes + `collect`, seconds); instrument JSONs `data/m5_22_4_fullf_{calib,all}.json` (regen: `python3 m5_22_4_a_fullf.py calib|all`, ~1 min); audit `data/m5_22_4_audit.json` (regen: `python3 m5_22_4_e_audit.py`, ~1 s); panels `plots/m5_22_4_{fullf,ladder}.png` (regen: `python3 m5_22_4_c_panels.py fullf|ladder`).

## Deviations log (live, per the flow doc)

| When | Deviation | Why |
| --- | --- | --- |
| 2026-08-06 11:05 | The generator envelope widened from the M5.21.3 electron's renv = 10 to renv = 18 for all baryon runs | The ring cores sit at r ≈ 16-17; the origin-centered renv = 10 envelope would give them ~1% weight (the M5.22.2 K1 lesson). Set at PLAN, logged here for visibility |
| 2026-08-06 11:10 | The 60-iteration p1 smoke test left a stale `row_p1_prot` + npz until the full run overwrote them ~16 min later | No effect on results (overwritten); noted because a mid-window checkpoint read in that window would have seen the stale row |

## TASK REVIEW (2026-08-06)

Task Duration: 01:52 (from 10:51 to 12:43 EDT)
Usage Cap Triggered: NO

Approved by the user 2026-08-06 (terminal review). Results: the full-F instrument built as the literal paper object with the DERIVED identity (its long-axis component IS the calibrated basic read, O(h²)-verified; flux diff ~1e-3 on all targets; middle/short non-quantized ≤ 0.30/0.014; norm read not quantized) ✅ · the static 4D lift survives on all four states + the s = −1 spot-check (charges exact, rings preserved and deepened, block-diagonality exact) ✅ · the generator catalog reproduces the electron sign table (every rotation positive incl. `clock_local` +0.154 to +0.773, boosts the only negative channel, spontaneous-twist term numerically zero) ✅ · the ω-ladder honest negative (no minimum at ω > 0 on any state; 20/20 rungs above matched controls; decoupling exact; charges preserved) ✅ · audit 5 PASS / 0 refuted with all 4 catches adopted ✅ · the author's 11:45 same-day reply captured (basic-sufficient scoping, the parameters caution, the time-crystal positioning; no pushback on the decline) ✅.

Issues: two logged deviations (renv = 18 envelope, set at plan; the 16-minute stale smoke-test row window, no effect). The proton is the boundary case of the decoupling claim (1.3e-3 vs 4-6e-4 elsewhere).

Action taken at approval: roadmap row → Done (appended at the end); "What happens next" + change-log updated (next = user pick among the conditional [M5.22.3](m5_22_3_task_details.md), [M5.22.5](m5_22_5_task_details.md), and the [M5.21.11](m5_21_11_task_details.md) parameter bridge both negatives point at); model-doc sweep already run at FINISH (canonical instrument row + 2 anti-recipe rows, briefing weak-force cell, MODELS.md beta-decay cell, tracker Q40). Checkers: docs ✅ 9 files, roadmaps clean, MODELS.md clean. Git stays the user's.

**Findings**: The author's full electric definition was built literally and measured to CONTAIN the calibrated instrument as its exact long-axis component, closing the instrument question with no numeric change; and the constant-ω long-axis-twist route to dynamical baryons is a clean measured negative on all four states (the electron's positive-kin decoupling generalizes unchanged), routing the dynamical-baryon question to free full-4×4 dynamics and/or physical parameters, exactly the two ingredients the author's same-day reply named.

**Research docs created/updated**: this task_details · [`../findings/m5_22_4_note.md`](../findings/m5_22_4_note.md) (the record) · [`m5_22_convo.md`](m5_22_convo.md) (the 11:45 capture) · [`../m5_question_tracker.md`](../m5_question_tracker.md) (Q40) · [`../m5_theory_canonical.md`](../m5_theory_canonical.md) · [`../../__M5_model_briefing.md`](../../__M5_model_briefing.md) · [`../../../../../MODELS.md`](../../../../../MODELS.md) · [`../m5_roadmap.md`](../m5_roadmap.md) · scripts [`m5_22_4_a_fullf.py`](../scripts/m5_22_4_a_fullf.py) · [`m5_22_4_b_omega.py`](../scripts/m5_22_4_b_omega.py) · [`m5_22_4_c_panels.py`](../scripts/m5_22_4_c_panels.py) · [`m5_22_4_e_audit.py`](../scripts/m5_22_4_e_audit.py) (auditor's) · data [`../data/m5_22_4_fullf_calib.json`](../data/m5_22_4_fullf_calib.json) + `_all` + 18 row JSONs + [`../data/m5_22_4_all.json`](../data/m5_22_4_all.json) + [`../data/m5_22_4_audit.json`](../data/m5_22_4_audit.json) + 9 local npz + the `_DATASETS.md` manifest · plots [`../plots/m5_22_4_fullf.png`](../plots/m5_22_4_fullf.png) · [`../plots/m5_22_4_ladder.png`](../plots/m5_22_4_ladder.png)

## Cross-links

| Doc | Why |
| --- | --- |
| [M5.22](m5_22_task_details.md) | The series plan home + the census states |
| [M5.21.3](m5_21_3_task_details.md) | The omega-ladder machinery being transplanted (also the object of the author's 2026-08-04 article-1 proposal on the lepton side) |
| [`m5_22_convo.md`](m5_22_convo.md) | The 2026-08-02 directive + the caution's exact wording |
| [`../m5_roadmap.md`](../m5_roadmap.md) | Backlog row + the 2026-08-05 re-plan change-log entry |
