# M5.22.2: beta decay of a kicked neutron-analog

**Status**: ✅ CLOSED 2026-08-05 (approved same day; review below). Staged 2026-07-30 at the M5.22 series re-plan (rung 4 of the author ladder as its own checkpoint-bounded subtask); the series goals, borders, and checkpoint policy live in [M5.22 § TASK PLANNING](m5_22_task_details.md). REFRESHED 2026-08-05 from the author's 2026-08-02 checkpoint reply ([`m5_22_convo.md § 2026-08-02`](m5_22_convo.md)): the expected decay channels named, the div E electric instrument staged as the opening move, the target set widened to both neutral states.

## Scope: rung 4, the second quantitative anchor

| Piece | Content |
| --- | --- |
| OPENING MOVE: the div E electric instrument | The author's 2026-08-02 correction: the M5.22.1 moments read (ρ = div B/4π) is MAGNETIC charge; the electric read needs div E, with E = the curvature of the long axis, or the full form of arXiv 2108.07896. Build BOTH variants, cross-check on the proton-analog (known electric +1), then RECOMPUTE the census + deuteron-candidate moments (the negative quadrupole sign is unsettled until this runs; [note § 8](../findings/m5_22_1_note.md)) |
| Protocol | The [M5.21.6](m5_21_6_task_details.md) kick protocol transplanted to the census neutron-analog; kick parameters self-scanned per the series exhaustion rule |
| Mechanism anchor | The author's 2026-07-27 panel: neutron → shift → split with energy release → reconnection → proton + electron + neutrino ([M5.22 § The charge-quantization frame](m5_22_task_details.md); sketch local-only in `theory/`) |
| The expected channels (author, 2026-08-02) | Neutron target: beta decay to proton + electron + neutrino (the neutrino as a topological vortex knot). Ring-antiring target: bineutron → proton + neutron + electron, or two neutrons. Caveat attached: toy parameters + the missing 3×3 angular momenta "could qualitatively change behavior", so a non-decay is a reportable outcome, not a failure |
| Reads | Decay products (proton-analog + fast escaping charged vortex + neutral ejecta); the β energy probability distribution SHAPE vs the experimentally known spectrum (shape only at toy parameters); the charge classes of all products by the winding/degree instrument + the new div E read |
| The Sulich composite read | Does the neutron-analog decompose under the kick as proton + electron (+ antineutrino), per the 2026-07-28 09:05 group question? Reported qualitatively from the decay products |
| The target set | BOTH neutral states: the census neutron-analog (primary; the M5.22.1 kick-apart did NOT split the two-ring state, so it survives as its own object) AND the ring-antiring state itself, probed against the author's named bineutron channels; each outcome routes the [Q40](../m5_question_tracker.md#q40-detail) identity by its decay products |

## Definition of done

A kicked-neutron decay run with products identified and the spectrum-shape read (or the non-decay characterized honestly); full closeout per the series policy (method note with inline images, tracker + [`MODELS.md`](../../../../../MODELS.md) "Weak force: beta decay (n → p)" cell update, roadmap move).

## Gated by

[M5.22.1](m5_22_1_task_details.md) close per the author ladder order + user "go". The physics needs only the census neutron-analog, so the user can pull this ahead of the deuteron on pick.

## FINDINGS (2026-08-05; the record is [`../findings/m5_22_2_note.md`](../findings/m5_22_2_note.md))

| Deliverable | Result |
| --- | --- |
| The div E instrument (both variants) | ✅ built + calibrated ([note § 2](../findings/m5_22_2_note.md)): only the arXiv 2108.07896 dual-curvature form passes Gauss-law calibration (analytic hedgehog → Coulomb r̂/r², 0.7% median error, flux +1.06; proton-analog/lepton quantized −0.98/−0.97); the literal field-line curvature (n̂·∇)n̂ FAILS quantization on every control (+7.0 / −4.4, flux-ladder inconsistent) |
| The instrument identity | ✅ E_full is numerically the SAME array the M5.22.1 moments called B: under the paper's dual mapping the M5.22.1 moment values stand as electric-sector reads. 🔶 The dual mapping itself is author-gated at the checkpoint |
| The frame-axis fork | ✅ only the LONG axis carries quantized flux (−1.02/−1.01/−1.01 charged, ~0 neutral); middle ≤ 0.34, short ≤ 0.01 |
| The moment recompute | ✅ the deuteron candidate's electric quadrupole sign stays NEGATIVE (−21.8 n32, −61.5 n48): the sign tension vs the physical deuteron is REAL under the calibrated instrument ([M5.22.1 note § 8](../findings/m5_22_1_note.md) updated); magnitudes not citable (resolution-drifting) |
| The decay probe | ✅ measured, the honest negative: NO decay in 20 runs on BOTH neutral states (K1 to 23×, K2 twists ≈ free by near-axisymmetry, K3 ring-localized to 53× the state energy); every endpoint returns with structure + slab charges intact. Neither named channel realized; [Q40](../m5_question_tracker.md#q40-detail) unresolved by products; the missing angular momenta ([M5.22.4](m5_22_4_task_details.md)) now the live suspect |
| The β-spectrum SHAPE anchor | ❌ not reachable this rung: no decaying channel to ensemble ([note § 5](../findings/m5_22_2_note.md)) |
| The Sulich composite read | ✅ qualitative NO: the neutron-analog does not decompose as p + e under any probed kick |

![div E calibration: the dual curvature is Coulomb on the hedgehog](../plots/m5_22_2_calib_coulomb.png)

![the kick ladder: every endpoint returns](../plots/m5_22_2_kick_ladder.png)

![the ring-antiring pair survives the largest ring-localized kick](../plots/m5_22_2_k3_before_after.png)

Datasets (heavy arrays LOCAL-ONLY, gitignored, kept; tracked = rows/JSON/plots): 20 kick rows + 2 extend rows `data/m5_22_2_row_*.json`, endpoints `data/m5_22_2_end_*.npz` (regen: `python3 m5_22_2_b_decay.py stage1` ~2.7 h + `stage2` ~2.6 h, sequential n = 32 runs, Apple Silicon); instrument JSONs `data/m5_22_2_dive_{calib,all,axes}.json` (regen: `python3 m5_22_2_a_dive.py calib|all|axes`, seconds-minutes each); panels `plots/m5_22_2_{calib_coulomb,kick_ladder,k3_before_after}.png` (regen: `python3 m5_22_2_c_panels.py`).

## Deviations log (live, per the flow doc)

| When | Deviation | Why |
| --- | --- | --- |
| 2026-08-05 11:39 | The stage1 kick ladder WIDENED from the staged K1:{0.4, 0.8} + K2:90 to K1:{0.05, 0.15, 0.4} + K2:{30, 90} | The smoke test measured the injection scale: K1:0.4 injects E ≈ 295 on E_start ≈ 12.7 (23×), so the staged ladder only blasted the state; a decay-barrier probe must bracket from below. Conservative widening, no rung dropped |
| 2026-08-05 14:35 | NEW kick family K3 (ring-localized envelope on one ring torus) added for stage2, ladder K3:{0.02, 0.05, 0.15, 0.4} | Stage1 measured the gap: the K1 envelope exp(−(r/8)²) is origin-centered, so the rings at r ≈ 17 received weight ~0.01: the decay targets were barely kicked while the column took the blast. K3 is the convert-one-neutron probe the channel prediction actually needs. Stage1 verdict so far: ALL kicks returned (K2 twists inject ≤ 0.19: the state is near-axisymmetric); the two K1:0.4 blasts converge back toward the start basins (extends running) |

## TASK REVIEW (2026-08-05)

Task Duration: 05:48 (from 11:28 to 17:16 EDT)
Usage Cap Triggered: NO

| Result | Status |
| --- | --- |
| The div E instrument, both variants + calibration | ✅ only the dual-curvature form quantizes (analytic hedgehog = exact Coulomb, audit derivative-exact; controls −0.98/−0.97); the literal field-line curvature fails every control |
| The instrument identity | ✅ E_full = the M5.22.1 "B" array (diff 0.0): the correction lands as a REREAD; 🔶 the dual mapping is author-gated |
| The moment recompute | ✅ the deuteron candidate's electric quadrupole stays NEGATIVE (−21.8 n32, −61.5 n48); magnitudes non-citable |
| Frame-axis fork | ✅ only the long axis carries the Gauss charge (audit-reproduced independently) |
| The decay probe | ✅ the honest negative: 20 runs, three kick families (to 53×), both neutral states, all returned with charges intact; no named channel opens |
| β-spectrum anchor / Sulich read | ❌ unreachable (no channel) / ✅ qualitative NO |
| Adversarial audit | ✅ 7 PASS, 1 PARTIAL, 0 refuted; all catches adopted (note § 6) |

Issues: two measured deviations (ladder widening; the K3 family), logged live. One 3-voxel residual on the P-1 K1:0.4 extension = the single literal ring-count exception (filtered, stated).

Action taken at approval: roadmap row → Done; Backlog resequenced ([M5.22.4](m5_22_4_task_details.md) ahead of [M5.22.3](m5_22_3_task_details.md)); model-doc sweep run (canonical recipes/anti-recipes + briefing weak-force line); MODELS.md beta-decay cell ❌ with scoreboard resync. The checkpoint comms: user-handled (the article-grade-negative item explicitly reserved to the user; the technical reply drafted terminal-only without it).

**Findings**: The div E correction resolves by measurement (only the dual-curvature form quantizes, it is the array already computed, the negative deuteron quadrupole is real under it), and the beta-decay probe returns a clean negative: both neutral ring-antiring states are deep minima to 53× their energy, so the burden moves to the missing angular momenta ([M5.22.4](m5_22_4_task_details.md)).

**Research docs created/updated**: this task_details · [`../findings/m5_22_2_note.md`](../findings/m5_22_2_note.md) (the record) · [`../findings/m5_22_1_note.md § 8`](../findings/m5_22_1_note.md) · [`../m5_question_tracker.md`](../m5_question_tracker.md) (Q40) · [`../../../../../MODELS.md`](../../../../../MODELS.md) · [`../m5_theory_canonical.md`](../m5_theory_canonical.md) · [`../../__M5_model_briefing.md`](../../__M5_model_briefing.md) · scripts [`m5_22_2_a_dive.py`](../scripts/m5_22_2_a_dive.py) · [`m5_22_2_b_decay.py`](../scripts/m5_22_2_b_decay.py) · [`m5_22_2_c_panels.py`](../scripts/m5_22_2_c_panels.py) · [`m5_22_2_e_audit.py`](../scripts/m5_22_2_e_audit.py) (auditor's) · data [`../data/m5_22_2_dive_calib.json`](../data/m5_22_2_dive_calib.json) + `_all` + `_axes` + 22 row JSONs + [`../data/m5_22_2_audit.json`](../data/m5_22_2_audit.json) + the `_DATASETS.md` manifest · plots [`../plots/m5_22_2_calib_coulomb.png`](../plots/m5_22_2_calib_coulomb.png) · [`../plots/m5_22_2_kick_ladder.png`](../plots/m5_22_2_kick_ladder.png) · [`../plots/m5_22_2_k3_before_after.png`](../plots/m5_22_2_k3_before_after.png)

## Cross-links

| Doc | Why |
| --- | --- |
| [M5.22](m5_22_task_details.md) | The census supplying the neutron-analog; the series plan home |
| [M5.21.6](m5_21_6_task_details.md) | The kick protocol + decay instruments being transplanted |
| [`m5_22_convo.md`](m5_22_convo.md) | The author-channel record, including the beta-decay mechanism rounds |
| [`../m5_roadmap.md`](../m5_roadmap.md) | Backlog row + the 2026-07-30 series re-plan change-log entry |
