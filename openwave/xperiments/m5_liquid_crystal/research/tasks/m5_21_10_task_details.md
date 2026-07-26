# M5.21.10: the decay-grade extension (N = 64 free arena + longer window)

**Status**: ✅ CLOSED COMPLETE 2026-07-20, review approved same day (~23:30 EDT), AUDITED (4 CONFIRMED / 3 NUANCE / 1 interpretation refuted + adopted) (go 15:24 EDT, user: "go 21.10 carrying 21.12"; staged 2026-07-19 from the electron-status gap review: the [M5.21.6](m5_21_6_task_details.md) open items get a task home). Rider [M5.21.12](m5_21_12_task_details.md) ran in the quiet slots, ✅ complete. Roadmap: [`m5_roadmap.md § IN PROGRESS`](../m5_roadmap.md).

## TASK PLANNING

**Arena decision (logged at PLAN)**: n = 64, L = 64, h = 1.0: the SAME lattice spacing as the M5.21.6 free arena (f48: n = 48, L = 48, h = 1.0), box grown 48 → 64 (+33% linear). This isolates the BOX variable (the named caveat on both M5.21.6 verdicts); a spacing change would confound it. The pinned census arena (n = 32, L = 48, h = 1.5) stays untouched as the statics reference.

| Phase | Content | Instrument | Machine-checkable gate (try cap 3) |
| --- | --- | --- | --- |
| P0 arena certification | The certified stack at n = 64 free: gate reruns (leapfrog γ = 0 E_tot conservation; sponge monotone drain; seed symmetry) + the descent regressions: A/B/C relaxed at f64 (snaps on), levels compared to the f48 record | `m5_21_6_a_decay.py` machinery via a thin n = 64 wrapper (`m5_21_10_a_decay64.py`), physics config IDENTICAL to f48 except L | gates green at n = 64; the A (electron) level consistent with f48 within the stated h-effect |
| P1 Read 1: the two-loop count + release kinematics | The μ-candidate C at dynamics grade in the bigger arena: evolve (damped wave + sponge) over a LONGER window; per-snap loop census (the M5.21.6 26-connected biaxial-core bookkeeping, thr {0.06, 0.09, 0.15} pre-registered) + the NEW release-velocity + direction read (per-component centroid track across snaps → velocity vector, ejection direction; the author's 2026-07-19 refinement: merger = parameter artifact, the two loops should eject in DIFFERENT directions) | evolve + loops modes extended with the kinematics read | loop count reported per thr with the connectivity convention disclosed; release velocity/direction measured for every departing component, or the honest "no separated component in-window" verdict |
| P2 Read 2: the B discrimination | The τ-candidate B at the bigger box: descent + dynamics window; discriminate DRAIN (energy to sponge, no core transition) vs DECAY (the rotation mechanism, core transition + release) vs HOLD; compare against the f48 drain record | same stack; the M5.21.6 rotation-vs-melt reads (spec_dev, frame rotation) | a stated B verdict at n = 64 with the f48 → f64 box trend; if B still drains, the box ladder is extended honestly, not closed |
| P3 Read 3: the ring tie-breaker | The point-vs-charged-ring core discrimination at next rigor (the M5.21.2 tie was instrument-limited: ring −3.7% fwd vs +23% 2h): both seeds relaxed on the certified symmetrized instrument at the census arena with a maxit/f_tol upgrade + the 2h re-read on the CONVERGED endpoints | `m5_21_2b_a_instrument.py` (T2, sym stencil), both seeds, cross-stencil re-read | the E ordering stated with cross-stencil consistency, or the tie honestly re-affirmed with the instrument limit quantified |
| P4 rider: the Q35 full read | [M5.21.12](m5_21_12_task_details.md) in the quiet slots (while f64 compute runs): sub-agent sweeps over the three arms (GR positive-energy; Ostrogradsky; constraint-carried/isorotation stabilization), citations verified before use, synthesis mapped onto the measured structure (ω\* = 0 / −∞ dichotomy, IR-extensive kinetic, the δ < 0 bounded family) | sub-agent fan-out + web verification; NO simulation | every citation chased (no unchased references); the synthesis lands in the M5.21.12 details + upgrades the 21.9/21.10 notes |
| P5 close | Findings note (method-note grade), panel plot, adversarial audit (independent agent, own reads), doc sweep | | audit verdicts recorded; doc checker exit 0 |

**Definition of done**: Reads 1-3 each carry a stated verdict (or the honest characterization of why the read fails) at the n = 64 arena + the Q35 synthesis landed; audited before anything author-facing.

**Blindspot pass**: (a) runtime/memory at 64³ (~2.4× f48): relax ~2.5 h/seed, evolve hours: run DETACHED (cap-surviving), checkpoint each JSON on arrival, floats32 snapshots; (b) the loop count is convention-sensitive: the M5.21.6 conventions (26-connectivity, thr ladder, edge-disjoint = closed candidate) are pre-registered here unchanged, sensitivity disclosed per thr; (c) the release-kinematics read is NEW instrumentation (flagged ours): centroid of each 26-connected component tracked across snaps, velocity by finite difference, direction as unit vector; degeneracy (components merging/splitting between snaps) handled by overlap matching and disclosed; (d) a bigger box EXTENDS the ladder, it cannot close it: if B drains at f64 too, the verdict is "drain persists at 64" not "B drains, period"; (e) config identity: f64 runs differ from f48 ONLY in n/L (asserted in-script).

**Unknowns routing**: machine-checkable → the phase gates; author-gated → none new this run (the two-loop conjecture and ejection-direction refinement are the author's own 2026-07-19 statements, tested as-is); nature-gated → none at toy parameters.

**Research-body destinations**: scripts `scripts/m5_21_10_*.py`, data `data/m5_21_10_*`, plots `plots/m5_21_10_*`, findings `findings/m5_21_10_note.md`, checkpoints `checkpoints/m5_21_10_progress.md`; the Q35 rider documents into `m5_21_12_task_details.md` (+ a findings home if the synthesis warrants one).

## RESULTS (2026-07-20 run; canonical record = [`../findings/m5_21_10_note.md`](../findings/m5_21_10_note.md))

| Phase | Verdict |
| --- | --- |
| P0 arena certification | ✅ gates GREEN at n = 64 (GK exact; GL1 dt² ladder 1.4e-3/3.6e-4/8.9e-5; GL2 monotone); dynamics ledger closes ≤ 3.1e-5 over t = 150 on all three windows; descents at the f48 budget: A retention 0.90 (better than f48), C the same two-eig decay signature but milder, B spreads with winding retained |
| P1 Read 1 (two-loop + kinematics) | 🔶 RE-BASED AT AUDIT: the paired SYMMETRIC ejection is real (split 1 → 3 at t ≈ 80, two off-center features departing at 0.052 in directions 108.9° apart, the same physical pair across thr 0.06/0.09; consistent with the author's different-directions refinement and unlike the control's mask erosion), but the audit's raw-snapshot recompute REFUTED the compact-fragment picture (late-time 1-cell z-column filaments, edge-connected beyond the census cut): the TWO-LOOP COUNT STAYS OPEN; the [M5.23.2](m5_23_2_task_details.md) arm-(1) tracer is the closure instrument |
| P2 Read 2 (B drain vs decay) | 🔶 the f48 "drain" verdict does NOT survive the bigger box: B DISINTEGRATES through the rotation channel (rot 20.1°, absorbed 3.79 ≈ half its budget, core 4298 → 316 cells migrating out, no in-place remnant > 11 cells); endpoint beyond t = 150 open |
| P2b A control | ✅ dual: energetic HOLD (absorbed 4.4% of budget; E −2.8% row-to-row / −5.6% vs the t = 0 budget, audit-corrected) = the stability clock-half box extension; + the INSTRUMENT NULL (both "departing" tracks are eroding masks near center) that calibrates Reads 1/2; the soft 12° late frame-rotation read confirmed no-melt (specdev monotone down) |
| P3 Read 3 (ring tie-breaker) | ✅ RESOLVED as a measured near-degeneracy, audit-CONFIRMED by independent recompute: all three seeds converge f_tol to ONE level (spread 3.78e-4; fields 1-2% apart; sym = fwd ordering, 2h flips it at that spread → below the cross-stencil floor; 2h absolute offset −7.64%); the M5.21.2 +23% 2h swing = convergence artifact |
| P4 rider M5.21.12 | ✅ COMPLETE: three verified arms + synthesis ([`../findings/m5_21_12_q35_read.md`](../findings/m5_21_12_q35_read.md)); headline: ghost-not-Ostrogradsky classification, fixed-J = energy-Casimir (HMRW 1985), radiation-window falsifier imported |
| P5 close | ✅ adversarial audit run (7 claims: 4 CONFIRMED / 3 NUANCE / 1 interpretation refuted + adopted; note § 8); doc sweep + manifest done |

**Deviations from plan** (logged live in [`../checkpoints/m5_21_10_progress.md`](../checkpoints/m5_21_10_progress.md)): (1) evolve dt = 0.025 instead of the gate auto-pick 0.05 (f48 comparability; 7× margin); (2) raw f64/t32 relax endpoints keep instrument-owned `m5_21_6_end_*` naming (the certified script owns its file naming; all analysis outputs are `m5_21_10_*`); (3) kin thr-sensitivity runs added beyond plan (thr 0.06/0.15 replicas) after the census proved threshold-noisy; (4) process note: two background completion notifications were swallowed (session freeze, cap suspected), costing ~100 idle minutes at the analysis handoff; detached compute lost nothing; (5) AUDIT CATCHES ADOPTED: the ev JSONs were missing the `snap=M_it1000` start-state field (added post-hoc with a disclosure note; a reproducer using the stored final `M` would run a different experiment) and the Read-1 compact-fragment interpretation was refuted by raw-snapshot recompute (note § 3 re-based to the filament-honest verdict; the census radii cut-sensitivity routed to the M5.23.2 tracer requirements).

**Datasets (local, gitignored, manifest regenerated at close)**: `m5_21_6_end_f64_{A,B,C}.npz` (~58 MB each, endpoints + descent snapshots; regen `python3 m5_21_6_a_decay.py relax seed=<S> n=64 L=64.0 bc=free maxit=12000 snaps=1 tag=f64_<S>`, ~150 min each) · `m5_21_10_ev_{A,B,C}_free64.npz` (~100 MB each, dynamics snapshots; regen `python3 m5_21_10_a_decay64.py evolve f64_<S> steps=6000 dt=0.025 snap_every=200 snap=M_it1000 out_tag=<S>_free64`, ~2.5 h each) · `m5_21_6_end_t32_{A,R4,R6}.npz` (small; regen `python3 m5_21_6_a_decay.py relax seed=<A|R> [aring=4|6] n=32 L=48.0 bc=pinned maxit=16000 tag=t32_*`, ~30 min). Tracked: all `m5_21_10_*.json` + logs + the panel.

## Scope (stub level)

| Piece | Content | Notes |
| --- | --- | --- |
| The arena | The honest free arena at N = 64 + a longer dynamics window (M5.21.6 ran N = 48) | Box size + horizon are the two caveats on the M5.21.6 verdicts |
| Read 1: the two-loop count | The released structure at dynamics grade: M5.21.6 measured ONE equatorial ring at descent grade vs the author's conjectured TWO loops ("four 1/2-vortices... two loops - hence two neutrinos") | The sharp discriminator on the author's decay-product conjecture; feeds the neutrino-side construction ([M5.20.8](m5_20_8_task_details.md)) |
| Read 2: the B discrimination | The τ-candidate B: drain vs true decay at the bigger box (B drained at N = 48; a box artifact is not excluded) | Decides whether the hierarchy's third level decays by the same rotation mechanism |
| Read 3: the ring tie-breaker | The point-vs-charged-ring core discrimination re-run at next rigor (the M5.21.2 tie was instrument-limited: ring −3.7% under fwd vs +23% under the 2h re-read) | The one same-charge alternative in the sanctioned term set ([`m5_particle_hunt.md`](../m5_particle_hunt.md) core-TYPE row) |

**Rider (user decision 2026-07-20)**: [M5.21.12](m5_21_12_task_details.md) (the Q35 full literature read) runs in THIS task's quiet slots (desk-only sub-agent sweeps while the N = 64 compute runs): its verified bibliography upgrades the M5.21.9/M5.21.10 notes before the consolidated send, and it produces no separate author-facing report.

**Reporting (user decision 2026-07-20)**: this task's note is one of the FOUR linked notes in the consolidated author send at the RENDERING UNLOCK milestone (21.9 + addendum, 21.10, 21.5, 21.4; packaging proposed at send time, grouped by concern: stability + spin / decay + neutrino count / films).

**Gated by**: the M5.21 series run order (after the fixed-J round, [M5.21.9](m5_21_9_task_details.md) incl. its round-2 addendum) + user "go".

## TASK REVIEW (2026-07-20)

Task Duration: 07:56 (from 15:24 go to 23:20 review; ~100 min lost to a suspected cap freeze that swallowed two completion notifications, detached compute unaffected)
Usage Cap Triggered: NO (resume ping pushed forward at both watchdog wakes, parked at finish, never fired)

| Phase | Verdict |
| --- | --- |
| P0 arena certification | ✅ gates GREEN at n = 64; ledger conserves at the 1e-4 grade over t = 150 on all three windows; descents at the f48 budget with clean cross-box signatures |
| P1 Read 1 two-loop | 🔶 RE-BASED BY THE AUDIT: paired SYMMETRIC ejection real (t ≈ 80 split, two off-center features at speed 0.052, directions 108.9° apart, same pair across thr 0.06/0.09, unlike the control null; consistent with the author's different-directions refinement); compact-fragment picture REFUTED (late-time 1-cell filament doublets, edge-connected beyond the census cut): the two-loop COUNT stays OPEN; [M5.23.2](m5_23_2_task_details.md) tracer = the closure instrument |
| P2 Read 2 B | 🔶 the f48 "drain" does NOT survive the bigger box: B DISINTEGRATES through the rotation channel (rot 20.1°, half the budget radiated, core volume ×0.15, no in-place remnant); endpoint beyond t = 150 open |
| P2b A control | ✅ electron HOLDS in free dynamics at n = 64 (4.4% absorbed) + the instrument null that calibrated Reads 1/2 |
| P3 Read 3 ring | ✅ RESOLVED, audit-confirmed: one level (spread 3.8e-4, below the cross-stencil floor); the M5.21.2 +23% 2h swing = convergence artifact |
| P4 rider M5.21.12 | ✅ COMPLETE (24 verified citations + synthesis): ghost-not-Ostrogradsky; fixed-J = energy-Casimir (HMRW 1985); dE/dJ = ω\* the textbook Q-ball identity; the RADIATION WINDOW imported as the sharpest falsifier |
| Audit | ✅ 7 claims: 4 CONFIRMED / 3 NUANCE / 1 interpretation refuted + adopted; `snap` provenance catch fixed with disclosure |

**Issues**: loop identity instrument-bound (census cut-sensitivity routed to the M5.23.2 tracer requirements); B endpoint needs a longer window; swallowed-notification failure mode cost idle time twice (process).

**Deviations**: five, logged live (§ RESULTS).

**Decision at review (user)**: checkpoint with the author BEFORE [M5.21.5](m5_21_5_task_details.md), sending TWO method notes only ([`m5_21_9_note.md`](../findings/m5_21_9_note.md) + [`m5_21_10_note.md`](../findings/m5_21_10_note.md); four would be too much to analyse); the M5.21.5/M5.21.4 notes follow at the rendering-unlock milestone.

**Findings**: the box extension flipped the τ-candidate's story from quiet drain to a violent rotation-channel disintegration, and promoted the two-loop question from unmeasurable to sharply instrumented (a real symmetric direction-separated ejection whose loop identity waits on the M5.23.2 tracer). The ring tie-breaker closed a two-week instrument ambiguity; the electron banked its free-dynamics hold; the Q35 read gave fixed-J its literature spine and sharpest falsifier; the audit caught a would-be overclaim before any author-facing text.

**Research docs created/updated**: [this task_details](m5_21_10_task_details.md) · [`../findings/m5_21_10_note.md`](../findings/m5_21_10_note.md) (audited § 8) · [`../findings/m5_21_12_q35_read.md`](../findings/m5_21_12_q35_read.md) · [M5.21.12 details](m5_21_12_task_details.md) · [`../m5_roadmap.md`](../m5_roadmap.md) · [`../m5_particle_hunt.md`](../m5_particle_hunt.md) · [`../m5_theory_canonical.md`](../m5_theory_canonical.md) · [`../../__M5_model_briefing.md`](../../__M5_model_briefing.md) · [`../m5_question_tracker.md`](../m5_question_tracker.md) (Q35) · [M5.23.2](m5_23_2_task_details.md) tracer requirements · scripts `m5_21_10_{a..d}` + data + `../plots/m5_21_10_panel.png` + manifest.
