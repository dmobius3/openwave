# M5.22.1: deuteron, the first nucleus (binding + quadrupole)

**Status**: ✅ DONE (run + closed 2026-08-01; the record is [`../findings/m5_22_1_note.md`](../findings/m5_22_1_note.md); review approved same day). Staged 2026-07-30 at the M5.22 series re-plan as the first checkpoint-bounded successor of the baryon census; the series goals, borders, and checkpoint policy live in [M5.22 § TASK PLANNING](m5_22_task_details.md). **The census checkpoint was answered 2026-07-30** ([`m5_22_convo.md § 2026-07-30`](m5_22_convo.md)): the reply re-scoped this task's constituents and staged the opening moves below. **This ID names a different task than pre-2026-07-30 records**: until the renumber, "M5.22.1" meant the far arc, now [M5.22.9](m5_22_9_task_details.md).

## Scope: rung 3 of the author ladder

| Piece | Content |
| --- | --- |
| Constituents | The proton-analog is author-CONFIRMED (the 2026-07-30 reply). The neutral constituent is OPEN: the author rereads the census two-ring state as candidate DINEUTRON (ring count = baryon number), so the opening moves below decide what the deuteron pairs with: the two-ring state as-is, or the single-ring neutral product of its kick-apart split |
| Construction | Two-knot bound-state attempt; ALL construction forks self-run per the series exhaustion rule (superposition scheme, initial separation, relative orientation, relaxation schedule): every branch executed, the pick made from data. The author's analytic deuteron seed is ANNOUNCED (2026-07-30) but not yet received: it folds in on arrival, and construction does not wait for it |
| Reads | Binding sign (m_d < m_p + m_n) and the electric quadrupole moment sign, the first quantitative anchor; both reported with magnitudes, both CLAIMED as signs only at toy parameters (declared default 6 of the census) |
| Escalation | The author's seed-formalization channel is used only if every self-run construction fork fails |

## Opening moves (staged by the 2026-07-30 reply, all self-runnable)

| Move | Content |
| --- | --- |
| The δ = 0.1 ladder | The author's direct ask: headline-state energies at δ = 0.1 (extends the census δ-probe, which read ordering 1.55 at δ = 0.2); reads whether the energy ladder tightens in the physical-δ direction |
| The kick-apart identity probe | The author's dineutron test: kick the two rings of the neutral state apart ([M5.21.6](m5_21_6_task_details.md) protocol), expected outcome two baryons; decides neutron vs dineutron, and a split hands this task its single-ring neutral constituent (and [M5.22.2](m5_22_2_task_details.md) its decay target) |
| Ring count as baryon number | The core-ledger ring count reported on EVERY endpoint as the baryon-number instrument (the author's proposed interpretation, testable against additivity in the binding runs) |
| The reversed electric convention | All signed charge plots/labels here and in later M5.22.x tasks carry electric = negated topological degree (the author's hedgehog-electron convention; [census note § 7](../findings/m5_22_note.md)) |

## Definition of done

A bound two-knot state exists (with binding + quadrupole read) or the failure is characterized honestly; either outcome is the checkpoint content (first nucleus, article-grade either way); full closeout per the series policy (method note with inline images, tracker + [`MODELS.md`](../../../../../MODELS.md) "Deuteron (binding + quadrupole)" cell update, roadmap move).

## Gated by

[M5.22](m5_22_task_details.md) closed ✅ and the census checkpoint answered 2026-07-30 (proton confirmed; the neutral identity routed to the opening probes). Remaining gate: user "go".

## RUN RECORD (2026-08-01)

The research body is [`../findings/m5_22_1_note.md`](../findings/m5_22_1_note.md) (method-note form: equations + code map + measured results + figures); this section is the task-side summary and the deviations log.

| Opening move | Outcome |
| --- | --- |
| The δ = 0.1 ladder | ✅ MEASURED (all f_tol, worst xr 1.45): ratios n/p 1.498 → 1.548 → 1.577 and p/e 1.409 → 1.446 → 1.464 across δ 0.3/0.2/0.1: existence δ-robust, ratios NOT converging by δ-steps alone ([note § 2](../findings/m5_22_1_note.md)); [Q33](../m5_question_tracker.md#q33-detail) receipt added |
| The kick-apart identity probe | ✅ MEASURED, honest negative vs the author's expectation: all 4 branches (split d = 3/6, kicks v = 0.5/1) return to the same topological sector + ring geometry (audit-qualified: the endpoints are nearby distinct stationary states within 0.35% in E); no decay into two baryons at toy parameters, pinned arena ([note § 4](../findings/m5_22_1_note.md)) |
| Ring count + internal structure | ✅ NEW STRUCTURAL FINDING: the neutral state's rings carry EXACTLY ±1 topological charges (audit exact degrees; FD slab ±1.04), dipole p_z = +24.0: a bound ring-antiring pair ([note § 3](../findings/m5_22_1_note.md)); identity question routed as [Q40](../m5_question_tracker.md#q40-detail) (author-gated, next checkpoint) |
| Reversed electric convention | ✅ applied: every signed panel labeled topological with electric = negated (note § 1 + panel (d) captions) |

| Construction fork | Outcome |
| --- | --- |
| seed2 p+n composites (a = 2, 3) | ❌ charge-ambiguous (q_far −0.75/−0.68); a = 3 stationary but ABOVE the constituent sum |
| seed2 pp control | 🔶 ESCAPED to a second ring-antiring neutral basin (E = 15.047, f_tol): ring count ADDS, charge CANCELS (the 3D escape law) |
| grafts (zoff = 6, 9) | ❌ stationary but unusable: the blended pinned far field is LAYERED (interior exactly Q = 0; the unit winding expelled to a frustrated pin-boundary seam; the fractional q_far reads are FD smears of the seam, audit-sharpened) |
| seedn three-center (one consistent \|Q\| = 1 far field) | 🔶 THE CANDIDATE: f_tol at E = 15.245, \|Q\| exactly 1 (electric +1), 2 rings with an internal +1/−1/+1 charge stack, 5.98 BELOW the constituent sum; electric quadrupole NEGATIVE (sign tension vs the physical deuteron); **xr 1.70 fails the n = 32 citation bar**; the first n = 48 pass is INCONCLUSIVE (E → 14.19 at max_iter, geometry unsettled, electric-quadrupole sign NEGATIVE at both resolutions; extension launched at close) ([note §§ 5b-5c](../findings/m5_22_1_note.md)) |

**Deviations from plan**: (1) two instruments ADDED beyond the staged list when the centered-cube fragment read proved geometrically wrong for rings at ρ ≈ 10: the lateral slab flux and the charge-density moments (both now in the note's equation map); (2) one fork ADDED mid-run (the three-center seedn composite) after the first five forks measured the two obstructions (escape-law charge cancellation, graft far-field inconsistency) that motivate it; (3) nothing from the staged plan was dropped.

**Artifacts**: scripts [`m5_22_1_a_kick.py`](../scripts/m5_22_1_a_kick.py) · [`m5_22_1_b_deuteron.py`](../scripts/m5_22_1_b_deuteron.py) · [`m5_22_1_c_panels.py`](../scripts/m5_22_1_c_panels.py) (+ the auditor's [`m5_22_1_e_audit.py`](../scripts/m5_22_1_e_audit.py)); data `m5_22_1_end_*.npz` (local, gitignored) + `m5_22_1_row_*.json` + [`m5_22_1_audit.json`](../data/m5_22_1_audit.json); plots `m5_22_1_slice_*.png`, [`m5_22_1_delta_ladder.png`](../plots/m5_22_1_delta_ladder.png), [`m5_22_1_kick_panel.png`](../plots/m5_22_1_kick_panel.png). Regen: each row JSON carries its CLI line; ladder runs ~160 s each, kick evolutions ~380 s, forks ~300-600 s (n = 32, M4).

## TASK REVIEW (2026-08-01)

**Task Duration:** 02:40 (from go 10:00 EDT to the terminal review 12:40 EDT; the n = 48 confirmation arm ran detached to 14:13 and was folded at 16:10, before approval)
**Usage Cap Triggered:** NO (the resume ping was armed, pushed forward once at the 13:30 watchdog, and parked unfired)

**Results**

| Result | Status |
| --- | --- |
| δ = 0.1 ladder: all three headline states persist as f_tol minima at δ = 0.3/0.2/0.1; n/p ratio 1.498 → 1.548 → 1.577, p/e 1.409 → 1.446 → 1.464: δ-steps alone do NOT tighten the ladder ([note § 2](../findings/m5_22_1_note.md)) | ✅ measured |
| The neutral state's rings carry EXACTLY opposite unit topological charges (audit degrees ±1.0000, p_z = +24.0, column charge 0): a bound ring-antiring pair | ✅ measured |
| Kick-apart: all 4 branches return to the same topological sector + two-ring geometry; NO decay into two baryons (endpoints nearby stationary states within 0.35% in E) | ✅ measured, honest negative |
| Construction obstructions: seed-level charge additivity dies by the 3D escape law (pp control → a second neutral basin at 15.047); grafted far fields expel winding to a frustrated pin-boundary seam | ✅ measured |
| Deuteron candidate (three-center seed): f_tol at E = 15.245, exactly \|Q\| = 1 (electric +1), 2 rings, internal +1/−1/+1 stack, 5.98 below the constituent sum; xr 1.70 fails the n = 32 bar; n = 48 (24000 it) flattens at 14.148 without f_tol or the two-ring geometry | 🔶 candidate, n = 32 evidence only |
| Electric quadrupole sign NEGATIVE in all three reads (−21.8 / −85.9 / −61.5) vs the physical deuteron's positive | ⚠️ sign tension, the run's most resolution-robust quantitative read |
| Adversarial audit: 4 CONFIRMED / 2 QUALIFIED / 0 refuted; every catch adopted ([note § 7](../findings/m5_22_1_note.md)) | ✅ |

**Issues / blockers**: the author's announced analytic deuteron seed never arrived during the run (folds in at the next checkpoint); the n = 48 convergence of the candidate stays open (slow anneal).

**Deviations from plan**: see [§ RUN RECORD](#run-record-2026-08-01) (two instruments added, one fork added, nothing dropped).

**Action needed**: the checkpoint outbound (Q40 + the candidate + the quadrupole-sign read, one batched note); [M5.22.2](m5_22_2_task_details.md) next on user "go"; MODELS.md deuteron cell unchanged until the candidate is citable.

**Findings**: The census neutral state is one deeply bound ring-antiring object (exact ±1 charges) that does not decay into two baryons under any probed kick, an honest negative for the dineutron reading, routed to the author as [Q40](../m5_question_tracker.md#q40-detail). A deuteron candidate with the right charge and ring count exists at n = 32, bound by 28% relative to its constituents, found only after the run measured why naive compositions fail (escape-law additivity, blended-far-field seams); its electric quadrupole sign is negative at every resolution probed, against the physical deuteron's positive, and its citability waits on n = 48 convergence. The δ ladder showed the mass ratios do not converge by δ-steps alone, reinforcing [Q33](../m5_question_tracker.md#q33-detail) as the quantitative bottleneck.

**Research docs created / updated**: [`../findings/m5_22_1_note.md`](../findings/m5_22_1_note.md) (the record) · this task_details · [`../m5_question_tracker.md`](../m5_question_tracker.md) (Q40 new, Q33 receipt) · scripts [`m5_22_1_a_kick.py`](../scripts/m5_22_1_a_kick.py) · [`m5_22_1_b_deuteron.py`](../scripts/m5_22_1_b_deuteron.py) · [`m5_22_1_c_panels.py`](../scripts/m5_22_1_c_panels.py) · [`m5_22_1_e_audit.py`](../scripts/m5_22_1_e_audit.py) (auditor's) · [`../data/m5_22_1_audit.json`](../data/m5_22_1_audit.json) + row JSONs + the `_DATASETS.md` manifest · plots [`m5_22_1_slice_P-1_plane_sc6_n32_pinned_d0.3.png`](../plots/m5_22_1_slice_P-1_plane_sc6_n32_pinned_d0.3.png) (the ring-antiring figure) · [`m5_22_1_slice_dn_-0.5at-2_+0.5at+0_-0.5at+2_n32_d0.3.png`](../plots/m5_22_1_slice_dn_-0.5at-2_+0.5at+0_-0.5at+2_n32_d0.3.png) (the candidate) · [`m5_22_1_delta_ladder.png`](../plots/m5_22_1_delta_ladder.png) · [`m5_22_1_kick_panel.png`](../plots/m5_22_1_kick_panel.png)

## Cross-links

| Doc | Why |
| --- | --- |
| [M5.22](m5_22_task_details.md) | The census this task consumes; the series plan + checkpoint policy home |
| [`m5_22_convo.md`](m5_22_convo.md) | The author-channel record for the whole nuclear program |
| [M5.22.2](m5_22_2_task_details.md) | The next rung (beta decay), gated on this close per the author ladder order |
| [`../m5_roadmap.md`](../m5_roadmap.md) | Backlog row + the 2026-07-30 series re-plan change-log entry |
