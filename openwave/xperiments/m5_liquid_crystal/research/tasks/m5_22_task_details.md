# M5.22: nuclei as vortex knots (the toy baryon census, on current hardware)

**Status**: 🚧 PLANNED, **THE MAIN PRIORITY since 2026-07-22** (the author's group reply to the M5.21 close-out: "I believe should be the main priority now", the convince-others vehicle for recruiting soliton-specialist help; user committed 2026-07-23 to a next-week start; thread decode: [`m5_22_convo.md`](m5_22_convo.md)). Staged 2026-07-18 from the author's autonomy goal list, goal (e); **TOY PHASE PULLED INTO 3×3 REACH 2026-07-19** by his 05:09 group message: "returning to 3x3 and studying simple vortex knots - which should correspond to proton and neutron ... with even toy simulations definitely deserving article"; [`m5_21_convo.md § 2026-07-19 05:09`](m5_21_convo.md)). This is the baryon/nuclei program (proton/neutron hunt); the toy phase is un-gated by the 4D fork. **Sequencing (2026-07-24 user reorder)**: ran AFTER [M5.27](m5_27_task_details.md) in the queue; M5.27 ✅ closed 2026-07-24 (structural null), so this task is now the Backlog head; the main-priority standing and the week-of-2026-07-27 start commitment are unchanged. **Carries the M5.27 hand-off** (§ The M5.27 hand-off below): the sharpened 4×4 time-sector problem statement rides this task's recruitment brief. **Series re-plan 2026-07-30 (user decision, pre-run review)**: re-scoped to the BARYON CENSUS subtask of the checkpoint-bounded M5.22.x series; goals, subtask borders, and the checkpoint policy in § TASK PLANNING below.

**SPLIT 2026-07-27 (user decision)**: this task is now the toy census ONLY, everything that runs on the hardware the project actually has. The upper rungs (larger nuclei, halo nuclei, alpha clusters, the table of nuclides, fusion/fission) moved to **[M5.22.9](m5_22_9_task_details.md)** (named M5.22.1 until the 2026-07-30 renumber), gated on compute the project does not own. Reason for the split: what is promised to the author should be exactly what can be delivered, and the rest should be named as a funded cost center rather than sitting inside a half-runnable task. Records written before the split call this scope "M5.22 Phase 1"; that is this task.

## The 3×3 vortex-knot baryon census, toy

His picture ([`../../theory/duda_2026-07-19_baryons_vortex_knots.png`](../../theory/duda_2026-07-19_baryons_vortex_knots.png)): proton = +e hedgehog with a vortex loop around it (lighter: encloses the charge); neutron = charge-compensated arrangement (−e/3 \| +2e/3 \| −e/3), heavier, positive core / negative shell; deuteron next.

| Piece | Content | Notes |
| --- | --- | --- |
| Arena | The certified 3D instrument ([M5.21.2b](m5_21_2b_task_details.md) T2 sym stack), n ≥ 48; seed construction = the author's three analytic cross-section families rotated around y = 0 (2026-07-28 13:20, charge-flip revision 2026-07-29; § The seeding update below, superseding both the hedgehog+loop phrasing and the transcription route) | 3×3, consumes delivered instruments; pin-mask minimizer extension now OPTIONAL |
| Pre-registered reads | (i) do the knot seeds relax to protected minima (the baryon analogs exist)?; (ii) **mass ordering: neutron-analog > proton-analog** (his compensation argument); (iii) **charge-distribution profile: positive core / negative shell for the neutron-analog** vs Wilson, PhysRevLett 7, 144 (his experimental anchor); (iv) deuteron-analog: two-knot bound state + electric quadrupole moment sign; (v) **same-charge mass hierarchy: proton-analog ≫ charged-lepton vortex at equal topological charge** (the "nearly 2000×" question, qualitative direction only; added 2026-07-28, § The seeding update); (vi) topology comparison: same winding/degree for the proton-analog and the charged-lepton vortex, DIFFERENT knot content by the component ledger (the Sulich challenge, added 2026-07-28) | Each read is a falsifiable claim from his picture; report whatever comes |
| Charge instrument | The winding/charge bookkeeping from the census era + the M5.21.6 topology instruments (component ledger, loop geometry) | Delivered |
| Article flag | His: "even toy simulations definitely deserving article" (the LC-topology venue family) | Human-owned prose over script-backed results, per AI_HYGIENE |

**The skeptic gauntlet (pre-registered honest caveats)**, from the 2026-07-18 X round ([`m5_21_convo.md`](m5_21_convo.md)): any vortex-knot nuclei claim must eventually answer (i) what plays SU(3) confinement / asymptotic freedom, (ii) sub-barrier fusion needs tunneling (the Gamow factor): classical knots deform, they do not tunnel, so either the quantized theory or an honest scope statement, (iii) nucleon statistics ([Q34](../m5_question_tracker.md#q34-detail)). This is the falsifier checklist for the claim language of this task AND of every successor subtask through [M5.22.9](m5_22_9_task_details.md), carried openly rather than buried. Item (i) is now COMMUNITY-POSED (2026-07-28, A. Sulich on the group list, [`m5_22_convo.md § 2026-07-28`](m5_22_convo.md)): if topology alone protects the proton, what stability role is left for the strong force? In this picture it must EMERGE as the knot/string tension (the Cornell arm, [Q38](../m5_question_tracker.md#q38-detail)), not be silently removed; the census note answers it head-on.

**Where this task stops (re-scoped 2026-07-30, § TASK PLANNING below).** This task runs the CENSUS ONLY: rungs 1-2 plus the two charge-sector arms ([Q37](../m5_question_tracker.md#q37-detail)/[Q38](../m5_question_tracker.md#q38-detail)), closing at the census checkpoint (the ranked lightest-baryon answer to the author's direct ask). The remaining toy rungs are checkpoint-bounded successor subtasks: [M5.22.1](m5_22_1_task_details.md) (deuteron), [M5.22.2](m5_22_2_task_details.md) (beta decay), [M5.22.3](m5_22_3_task_details.md) (the conditional A = 3-4 stretch). Everything above (rungs 7-8, the table of nuclides, fusion/fission) is [M5.22.9](m5_22_9_task_details.md) (the far arc) and needs compute the project does not have.

**Gated by**: user "go" (the M5.21 series is closed; the author-named main priority; start committed for the week of 2026-07-27).

**Pre-run step (user, 2026-07-28)**: at "go", FIRST review the then-current state of the [`MODELS.md`](../../../../../MODELS.md) criteria table before locking the run list. [T1](../../../../../dev_docs/tasks/t1_task_details.md) (platform roadmap) ✅ RAN 2026-07-28 (criteria set 22 → 31): the Baryons split landed (mass ordering + core/shell \| exact masses \| beta decay), plus the deuteron-binding + quadrupole row, the nuclear-structure row, and per-row author-named tests with priority tiers ([`t1_convo.md`](../../../../../dev_docs/tasks/t1_convo.md)). The pre-registered reads and the claim language of this task must target the criteria rows as they stand on run day, not as they stood at staging; the run-day review confirms no further restructure landed in between.

## TASK PLANNING (2026-07-30 pre-run review)

**Scope**: the baryon census, rungs 1-2 of the author ladder plus the two charge-sector arms. Implement the 2026-07-29 seed set (12 seeds: electron/positron at R = 2 in both signs; the central-π-step and fractional-step families at R = 1 for s ∈ {-1, -1/2, 0, 1/2, 1}), rotate each half-plane cross-section around y = 0, relax by gradient descent on the certified 3D stack with FULL 3D freedom (cylindrical symmetry watched, not imposed); dedup the relaxed states by invariants (energy, 3D degree, component ledger); rank per measured charge class under the selection principle; read the mass ordering + the core/shell profile; run the lepton-reference control ([Q37](../m5_question_tracker.md#q37-detail) rides it) and the quark-shift scan ([Q38](../m5_question_tracker.md#q38-detail)). **The seed-set version is PINNED at "go"** (the 2026-07-29 charge-calc PDF); author revisions arriving mid-run fold in at the next checkpoint, not mid-run.

**Definition of done**: every seed relaxed to convergence or its failure characterized, with restart + resolution robustness; a deduped census table (state, energy with error bar, charge class, topology ledger); the proton/neutron-analog assignment per the selection principle, or an honest tie report; pre-registered reads (i)-(vi) answered; Q37/Q38 measured; the method note (equations first, inline images, equation-to-code map) presented at the checkpoint; [`MODELS.md`](../../../../../MODELS.md) cells + trackers updated at close.

**Gating**: user "go". Instruments delivered ([M5.21.2b](m5_21_2b_task_details.md) T2 sym stack, n ≥ 48; winding/degree + component ledger; the [M5.21.4](m5_21_4_task_details.md) tension instrument for Q38); the pin-mask extension OPTIONAL.

**Blindspot pass**: the unknown-unknowns row of § The unknowns map below.

**Research body**: scripts/data/plots in `research/` under `m5_22_` naming; findings note `findings/m5_22_note.md`; census table + key plots embedded inline in this doc at close.

### The series goals

Claim language series-wide: **qualitative at toy parameters** (the 2026-07-23 posture); quantitative reads wait on the [Q33](../m5_question_tracker.md#q33-detail)/[M5.21.11](m5_21_11_task_details.md) bridge. The skeptic gauntlet rides every note.

| Goal | Detail | Scored by |
| --- | --- | --- |
| Baryon existence + identity | The 12 seeds relax to protected minima; deduped, ranked per measured charge class; proton-analog = lightest charged, neutron-analog = lightest neutral; heavier states reported as candidate excited baryons | [`MODELS.md`](../../../../../MODELS.md) "Baryons: bound state (p, n)"; the author's direct ask answered |
| Mass ordering + charge profile | Neutron-analog HEAVIER (the compensation prediction); positive core / negative shell vs Wilson PRL 7, 144 | "Baryons: mass ordering + charge profile" |
| Charge-sector demonstrations | [Q37](../m5_question_tracker.md#q37-detail): far-field degree(proton-analog) = degree(positron-analog) exactly; [Q38](../m5_question_tracker.md#q38-detail): quark-shift energy scan vs the [M5.21.4](m5_21_4_task_details.md) string tension, aimed at the Cornell ~1 GeV/fm frame | "Charge quantization" + "Quarks" + "Strong force: confinement"; both group-visible asks |
| Deuteron | Two-knot bound state: binding sign (m_d < m_p + m_n) + electric quadrupole moment sign = the first quantitative anchor | "Deuteron (binding + quadrupole)"; [M5.22.1](m5_22_1_task_details.md) |
| Beta decay | Kicked neutron-analog → proton-analog + fast escaping charged vortex (+ neutral ejecta); β-spectrum SHAPE vs the known distribution = the second anchor; the Sulich composite question read from the decay products | "Weak force: beta decay (n → p)"; [M5.22.2](m5_22_2_task_details.md) |
| Recruitment portfolio | Article-grade per-rung notes ("even toy simulations definitely deserving article"); the sharpened 4×4 problem statement (§ The M5.27 hand-off) rides the brief | The author's convince-others strategy, the reason this series is the main priority |

### The nuclear-hunt series plan (borders = author checkpoints)

Each subtask CLOSES with the full workflow (finish, review, approval, task-details documentation, method note with inline images, question-tracker + criteria-cell updates, roadmap move). **That closeout IS the checkpoint**: the note is the guidance request.

| Subtask | Scope | Why the border sits here |
| --- | --- | --- |
| **M5.22** (this task) | The baryon census: rungs 1-2 + the Q37/Q38 arms | Triple trigger: rungs 1-2 complete; the author's direct ask ("search especially for the lightest baryons") answered; and the review of the ranked census (which states are p/n vs excited baryons vs artifacts, whether the seed set extends) is SUBSTANTIAL guidance gating every downstream construction |
| [M5.22.1](m5_22_1_task_details.md) | Deuteron (rung 3): two-knot construction from census-confirmed constituents; binding + quadrupole signs | First NUCLEUS + first quantitative anchor, article-grade either way; constituents must be checkpoint-confirmed first |
| [M5.22.2](m5_22_2_task_details.md) | Beta decay (rung 4): the [M5.21.6](m5_21_6_task_details.md) kick transplanted; products + β-spectrum shape; the Sulich composite read | Dynamic decay result, article-grade; the mechanism panel exists, so no guidance needed to start, only to close |
| [M5.22.3](m5_22_3_task_details.md) | The A = 3, 4 stretch (rungs 5-6): tritium → He3, He4 | CONDITIONAL on feasibility measured in the earlier subtasks; otherwise these rungs fold into the far arc and the toy series truncates |
| [M5.22.9](m5_22_9_task_details.md) | The far arc (rungs 7-8 + its A-F ladder), compute-gated | Unchanged: the named cost center |

### The checkpoint policy (series-wide)

| Rule | Content |
| --- | --- |
| Exhaustion rule | A fork decidable by running scripts (even slow ones) is RUN, not asked: every branch executed, the pick made from data, documented in the task details. The checkpoint note carries fork RESULTS, never fork questions |
| Trigger 1: substantial guidance | The next construction consumes an author-gated answer (state identity, seed extension, structure interpretation) |
| Trigger 2: article-grade finding | A positive OR negative finding that would change the author's picture goes out at the next checkpoint, not sat on |
| Trigger 3: rung completion | Each ladder rung landing is a checkpoint by default |
| Escalation channel | The author's standing "let me know if there are problems" is used only AFTER our own forks fail, never instead of them |
| Batching | One outbound per checkpoint, method-note-grade context per question (the tracker's outbound policy, [`m5_question_tracker.md § Email cadence`](../m5_question_tracker.md)) |
| Substrate constraint | Series-wide 3×3 ONLY: no 4×4 runs inside M5.22.x; the 4×4 material rides only as the recruitment-brief hand-off (§ The M5.27 hand-off) |

### The unknowns map (quadrant-routed; successors inherit the pattern)

| Quadrant | Content → route |
| --- | --- |
| Known knowns | The certified instrument stack; the 12 analytic seeds with charge class known pre-relaxation; the selection principle; pre-registered reads (i)-(vi); the six declared defaults (§ The seed prescription); the skeptic gauntlet |
| Known unknowns, machine-checkable (exhaust pre-checkpoint) | Cylindrical-symmetry survival under full-3D relaxation (author-flagged "or not?"); lattice resolution for R = 1 cores + box-size scan; boundary conditions for net-charged states in a finite box; relaxation hyperparameters + restart robustness; dedup of relaxed states by invariants; what the s = 0 and \|s\| = 1 seeds relax to; the lepton control (family 1 → the M5 charged ring, the author's own "(?)"); 3D degree vs 2D pre-rotation winding per state; ranking error bars + a pre-registered tie protocol; a small (g, δ) scan so the mass ordering is not a parameter accident; [Q37](../m5_question_tracker.md#q37-detail); [Q38](../m5_question_tracker.md#q38-detail) |
| Known unknowns, author-gated (the checkpoint batch) | Identity of the ranked states against the baryon list; whether the seed grid extends; the cores-vs-shells reading of the relaxed neutron candidate; the Cornell-vs-string-tension identity (the author's own open 14:06 question); deuteron seed guidance ONLY if all self-run construction forks fail |
| Known unknowns, nature-gated | A = 3-4 feasibility on owned hardware (measured during the census, decides the truncation point); quantitative masses/moments ([Q33](../m5_question_tracker.md#q33-detail) bridge); the far arc's compute |
| Unknown unknowns (blindspot pass) | Lattice-pinning artifacts masquerading as census states (mitigation: resolution + restart ensembles); the π₁ ½-line split contaminating the π₂ charge read ([Q31](../m5_question_tracker.md#q31-detail) tension, rides Q37); mid-run seed revisions (mitigation: the pinned seed-set version); two lightest states tied within error (mitigation: report both, no forced pick) |

## The author's public framing question (2026-07-20, captured)

On the 2026-07-20 group thread ([`m5_21_convo.md § 2026-07-20 13:44`](m5_21_convo.md)) the author posed Phase 1's program to the group verbatim: "Then there are various size knots, and various size nuclei - is there correspondence between them? How does it look like?". The knots-nuclei correspondence is now a publicly seeded community question, which strengthens the article case for the Phase 1 toy census ("even toy simulations definitely deserving article") and gives the eventual note its opening frame.

## The direct request (2026-07-21, captured)

In the reply to the two-note checkpoint on our own thread ([`m5_21_convo.md § 2026-07-21 03:30`](m5_21_convo.md)) the author named this program the priority interest for the 3×3 stack: "as there is still problem with 4x4, the most interesting for 3x3 would be studying nuclei as vortex knots - starting with getting proton lighter then neutron, then deuteron with electric quadrupole momentu ... then larger vortices hopefully getting close 'periodic table' with some qualitative agreements e.g. for electromagnetic moments ... halo nuclei". Two additions to the Phase 1 frame: the ladder now runs explicitly proton/neutron mass ordering → deuteron quadrupole → toy periodic table (electromagnetic moments as the qualitative scorecard) → HALO NUCLEI (new far rung), and the mass-ordering read is restated in the same direction as the 2026-07-19 picture (proton lighter, neutron heavier by charge compensation). Sequencing unchanged: user "go" after the M5.21 series (the staged M5.21.5 → M5.21.4 → RENDERING UNLOCK route stands; this arc is the author-endorsed next physics program after it).

## The nuclei-first directive (2026-07-22 group reply + 2026-07-23 follow-ups, captured)

The author's group reply to the M5.21 close-out (full decode: [`m5_22_convo.md`](m5_22_convo.md)) upgrades this task from "next physics program" to **the main priority**, with an explicit strategy: the 4×4 blockers (physical g/δ beyond numerics; the gravity+clock case) are author-open and may need soliton specialists, and qualitative 3×3 nuclear results are the vehicle "to have a chance to convince others" to bring that help. The 4×4 case stays green-lit for self-directed work in parallel ("if Fable has idea how to move forward with 4x4 case, please work on it"), with Skyrme-style extra Lagrangian terms floated as a candidate cure (screened before any spend by the M5.20.4 arm-C lemmas + the M5.20.5 γ = −1 statics kill).

**The run ladder** (his 2026-07-22 sequence, superseding the 2026-07-21 ordering; sketch filed local-only at [`../../theory/duda_2026-07-22_nuclei_beta_decay_sketch.png`](../../theory/duda_2026-07-22_nuclei_beta_decay_sketch.png)):

| Rung | Target | Read / anchor |
| --- | --- | --- |
| 1 | Proton analog (seed = the author's fractional-step family at s = -1/2, charge +1 per the 2026-07-29 flip; § The seeding update below, superseding the two-vortex-prescription seeding of § The seed prescription) | protected minimum exists |
| 2 | Neutron analog (charge-compensated; seed = the author's central-π-step family at s = -1/2, charge 0 per the 2026-07-29 flip; § The seeding update below, superseding the field-rotation route) | HEAVIER than the proton analog (compensation argument); positive core / negative shell (Wilson PRL 7, 144) |
| 3 | Deuteron → **[M5.22.1](m5_22_1_task_details.md)** (2026-07-30 series re-plan) | simplest nuclear binding; **electric quadrupole moment**, experimentally known = the first quantitative anchor |
| 4 | **Beta decay of a kicked neutron** (NEW arm) → **[M5.22.2](m5_22_2_task_details.md)** | the [M5.21.6](m5_21_6_task_details.md) kick protocol transplanted; decay must also release a fast electron; the experimentally known β energy probability distribution = the second quantitative anchor |
| 5 | Tritium → He3 → **[M5.22.3](m5_22_3_task_details.md)** (conditional) | the analogous beta decay at A = 3 |
| 6 | He4 → **[M5.22.3](m5_22_3_task_details.md)** (conditional) | "the crucial alpha" |
| 7 | Larger nuclei (far rung) → **[M5.22.9](m5_22_9_task_details.md)** | stable halo nuclei; alpha clusters (C12, O16), compared against the Skyrmion models of nuclei (the named baseline literature) |
| 8 | Strange-baryon decay (far rung, added 2026-07-26) → **[M5.22.9](m5_22_9_task_details.md)** | decay releasing a pion or a kaon; the meson side of the same reconnection mechanism (the standing paper-level reading: pion = twist/reconnection of a vortex loop, kaon = a Möbius-like twisted loop with strangeness as the twist) |

The 2026-07-23 1:1 follow-ups (summarized in the convo doc): if the 4×4 case stays stuck, nuclei-first is the sequencing; the user committed to starting the week of 2026-07-27; and the expected grade is restated honestly: even with inappropriate (toy) parameters, MANY QUALITATIVE agreements should be in reach, and that is the claim language this task carries (quantitative reads wait on the [Q33](../m5_question_tracker.md#q33-detail)/[M5.21.11](m5_21_11_task_details.md) bridge).

## The charge-quantization frame (2026-07-26/27 group thread, captured)

The 2026-07-27 exchange with Adrian Sulich on the group list ([`m5_22_convo.md § 2026-07-27`](m5_22_convo.md)) supplies the WHY behind rungs 1-3, and it is the frame the Phase-1 note should be written in. The author's argument, in one line: electric charge is the Gauss-law integral read as the degree of a deeper field (`∮_S E·dA = (e₀/4π) ∮ du dv (∂_u n̂ × ∂_v n̂)·n̂`, so charge is topological and cannot halve), the vortex interaction inside a baryon deforms the field toward a hedgehog, the proton closes that into one full elementary unit while the neutron has to compensate it, and the compensation is what makes the neutron heavier, gives it a positive core with a negative shell, and gives the deuteron its quadrupole moment. Three consequences for this task:

| Consequence | What it changes |
| --- | --- |
| The mass ordering is not a separate posit | Rung 2's pre-registered read (neutron-analog heavier) is a PREDICTION of the charge-quantization structure, not an independent assumption. If the census finds the ordering, the note can say why; if it finds the opposite, the failure is structural, not cosmetic |
| Two external questions now ride the ladder | [Q37](../m5_question_tracker.md#q37-detail) (why is a composite's charge magnitude exactly the electron's, and can the winding instrument show it) and [Q38](../m5_question_tracker.md#q38-detail) (the ~1 GeV/fm Cornell cost of fractional charge, against the linear string term already measured at [M5.21.4](m5_21_4_task_details.md)). Both are answerable on the certified 3D stack, and both are group-visible asks |
| The beta-decay mechanism is specified, not just named | The author's 2026-07-27 panel draws rung 4 as a sequence: neutron → shift → split with energy release → reconnection → proton + electron + neutrino, with the electron as a vortex INSIDE and the muon as a vortex loop AROUND. That is the protocol shape for the kicked-neutron arm (local-only sketch: [`../../theory/duda_2026-07-27_baryons_knots_beta_decay.png`](../../theory/duda_2026-07-27_baryons_knots_beta_decay.png)) |

## The seed prescription (2026-07-27 author reply, captured): the go-pack for rungs 1-2

The 2026-07-27 afternoon reply ([`m5_22_convo.md § 2026-07-27 afternoon`](m5_22_convo.md)) answered the pre-run questions with an operational recipe rather than point-by-point. This section is what the run consumes.

| Piece | Prescription |
| --- | --- |
| Proton seed (rung 1) | Two vortices along DIFFERENT axes: the straight vortex (the lightest excitation) along the LONGEST axis, the second-lightest vortex looping around it. Their interaction suggests the different-axes arrangement. Prepare both, run energy minimization; the expected outcome is a stable configuration corresponding to the proton. This supersedes the earlier hedgehog+loop seed reading of the 2026-07-19 picture: the hedgehog is what the vortex interaction DEFORMS the field toward (the charge-quantization frame above), not what gets seeded. Cross-section panel local-only: [`../../theory/duda_2026-07-27_proton_two_vortex_seed.png`](../../theory/duda_2026-07-27_proton_two_vortex_seed.png) |
| Neutron route (rung 2) | Built AFTER the proton lands, by adding field rotations to bring the total electric charge to 0. Detail deferred by the author to that point |
| The Cornell arm (NEW, feeds [Q38](../m5_question_tracker.md#q38-detail)) | On the neutron analog, try shifting its quarks and read the energy cost of the displacement, aiming at the ~1 GeV/fm Cornell scale. This is the author's own in-model version of the fractional-charge cost question, and it composes with the [M5.21.4](m5_21_4_task_details.md) string-tension instrument. The 14:06 question (is our measured linear string term the same object as the Cornell tension?) is still open; the quark-shift scan is how it gets an in-model answer |
| Escalation path | "If there are issues, please let me know and I will try to formalize it tomorrow": a same-week formalization channel is explicitly open for the seed construction. Use it instead of guessing past a blocker |

**Status of the six pre-run questions (asked 2026-07-27 14:06).** None was answered point by point; the ledger below states what replaces each, so the run starts on declared defaults instead of silent ones.

| # | Question | Status → what the run does |
| --- | --- | --- |
| 1 | What separates the two vortex types in the field | OPEN as theory; OPERATIONALLY REPLACED: the recipe says "the lightest" and "the second lightest" vortex excitations, so the run identifies them by relaxed energy on the certified stack, not by a priori classification |
| 2 | Expected loop radius vs core scale | ANSWERED implicitly: energy minimization selects it; no target imposed |
| 3 | Neutron: three cores vs one core + shells | DEFERRED 2026-07-27, SUPERSEDED 2026-07-28 13:20 (concrete neutron seed supplied), REVISED 2026-07-29: the charge flip makes the neutron candidate the central-π-step seed (ONE regularized core), while the three-fractional-core seed carries charge +1 (the proton candidate); the relaxed states answer cores-vs-shells empirically either way (§ The seeding update below) |
| 4 | Sector parameters for nuclear runs | OPEN → DEFAULT: same (1, δ, 0) spectrum and toy δ as the lepton work, stated as an assumption in the note; revisit only if the proton seed fails to stabilize |
| 5 | Charge instrument = far-field degree / Gauss law | OPEN → DEFAULT: yes, the winding/degree instrument as in the author's sketch; neutron reads 0 in total, core/shell structure read from the radial profile |
| 6 | Deuteron binding readable at toy parameters, or sign only | OPEN → DEFAULT: report both the sign and the magnitude, claim only the sign (qualitative grade per the 2026-07-23 posture) |

**Unblocked verdict**: rung 1 is fully unblocked (seed recipe + escalation path). Rung 2 is directionally unblocked (route named, detail deferred until rung 1 lands). Rungs 3-4 are unchanged (their anchors were already pre-registered). The declared defaults above are the only assumptions in play.

**Close-out routing (user, 2026-07-27)**: the method note for this task goes to the author cc'ing the colleague from the 2026-07-27 reply (the nuclei-interested third party; address kept off the public repo, on local file). Resource-contributor outreach stays HELD until the user calls its moment.

## The seeding update (2026-07-28/29, captured): the automatic protocol, the ansatz retraction, the three seed families, then the charge flip

The 2026-07-28/29 rounds ([`m5_22_convo.md § 2026-07-28`](m5_22_convo.md) and following) OPERATIONALIZE the rung-1 construction above in three steps: the morning replaced the closed-formula hope with a constrained-relaxation route, the 13:20 mail (now cc'd to the group) DELIVERED analytic cross-sections after all, three seed families ready to rotate and relax, and the 2026-07-29 follow-up added a seed-level total-charge calculation that FLIPS the two baryon candidate labels. This section is what the run now consumes for seeding; the two-vortex prescription above remains the physics picture it implements.

| Piece | Content |
| --- | --- |
| The automatic protocol (PRIMARY route, 01:17) | (1) Transcribe the cross-section panel using ONLY its black segments = long-axis angles (the eigenvector of the largest eigenvalue) at their positions, and extend the central vortex parallel to the end of the box; (2) impose cylindrical symmetry by rotating around the bottom edge; (3) fill the rest of the lattice by energy minimization with ONLY those long axes CONSTRAINED; (4) release the constraints, minimize further: gradient descent should reach at least the proton (the neutron is harder). Panel to transcribe, local-only: [`../../theory/duda_2026-07-28_seed_cross_section_black_segments.png`](../../theory/duda_2026-07-28_seed_cross_section_black_segments.png) |
| New instrument requirement | Constrained relaxation: a pin mask holding selected director long axes fixed while the rest of the lattice minimizes. Not in the delivered toolkit; small extension of the certified minimizer. Demoted to OPTIONAL by the 13:20 analytic seed families below (they fill the whole cross-section, no constrained fill needed) |
| The analytic ansatz, RETRACTED same morning | 06:15: `ang = -(ArcTan[-y+R, -x] - ArcTan[R, -x])/2` (R = outer-loop radius), rotated around y = 0, second axis perpendicular in plane, claimed to give an external hedgehog + hidden vortex loop ([`../../theory/duda_2026-07-28_proton_ansatz_mathematica.pdf`](../../theory/duda_2026-07-28_proton_ansatz_mathematica.pdf)). 08:02, the author's own correction: rotating that configuration gives the charged vortex of the ELECTRON/POSITRON, not the proton; for the proton the hedgehog must sit BETWEEN the vortices, per the arXiv diagram ([`../../theory/duda_2026-07-28_arxiv_baryon_diagram.png`](../../theory/duda_2026-07-28_arxiv_baryon_diagram.png)) |
| What the retraction buys | The retracted ansatz is RETAINED as the **charged-lepton reference seed**: same topological charge as the proton target with far less structure. It is the control object for pre-registered reads (v) (same-charge mass hierarchy, the author's "nearly 2000 heavier" frame) and (vi) (same degree, different knot content, the Sulich challenge) |
| Hedgehog placement is load-bearing | The one structural invariant the retraction establishes: hedgehog BETWEEN the vortex cores = baryon; hedgehog closed by a single vortex = charged lepton. The seed transcription must respect it, and the relaxed configurations should be checked against it |
| Escalation channel | "Let me know if there are problems" (01:17): the same-week formalization channel stays open; the author is also thinking about the neutron modification |

**The three seed families (13:20, group-public; verbatim record + formulas in [`m5_22_convo.md § 2026-07-28 13:20`](m5_22_convo.md), PDF local-only at [`../../theory/duda_2026-07-28_electron_proton_neutron_seeds.pdf`](../../theory/duda_2026-07-28_electron_proton_neutron_seeds.pdf))**: half-plane cross-sections of the director angle, each in both signs s = ±1/2, to be rotated around the bottom y = 0 axis (cylindrical symmetry, author-flagged as an open assumption: "or not?"), second axis perpendicular in the plane, then energy-minimized by gradient descent on the 3×3 stack.

| Family | Construction | What the run does with it |
| --- | --- | --- |
| Electron/positron (charged vortex) | `ang = -s (ArcTan[-y+R, -x] - ArcTan[R, -x])`, R = 2: the retracted morning ansatz relabeled to what the retraction said it is, now in both signs | The charged-lepton reference seed is now AUTHOR-AUTHORED: relax it and check it reproduces the M5 charged ring (the author asks this with a question mark); it anchors reads (v)-(vi) |
| Central-π-step family (13:20 "proton?"; **NEUTRON candidate since the 2026-07-29 charge flip**) | The same hedgehog pair plus a smoothed central π step (`step[x] = ArcTan[x]/π + 1/2`, term `-π·step[5x]`): the central charge ACTIVATED, the hedgehog-between-cores structure in closed form. Far-field charge ≈ -2s - 1, so **s = -1/2 is NEUTRAL** | Rung-2 seed (supersedes "field rotations after the proton lands"): evaluate the formula on the half-plane grid, rotate, relax; no transcription needed |
| Fractional-step family (13:20 "neutron?"; **PROTON candidate since the 2026-07-29 charge flip**) | Fractional steps instead of the central π: `-2π/3` at center, `+π/3` at x = ±1, the -e/3 \| +2e/3 \| -e/3 bar of the arXiv diagram, regularized in the centers. The fractional steps CANCEL in the far field: charge ≈ -2s exactly, so **s = -1/2 carries +1** | Rung-1 seed: relax alongside the neutron candidate; the d-u-d bar seed turns out to be the CHARGED one |

**The charge-calculation revision (2026-07-29, [`m5_22_convo.md § 2026-07-29`](m5_22_convo.md); PDF local-only at [`../../theory/duda_2026-07-29_proton_neutron_charge_calc.pdf`](../../theory/duda_2026-07-29_proton_neutron_charge_calc.pdf))**: the author added a seed-level total-charge instrument to the notebook, the far-field winding `(1/π) ∫₀^π ∂ang/∂φ dφ` along a semicircle of radius d = 10, and reported "it seems the opposite: previous upper for s = -1/2 seems neutron, lower for s = -1/2 seems proton". The same revision re-issues both baryon families with **R = 1** (was 2) and the sign range extended to **s ∈ {-1, -1/2, 0, 1/2, 1}**, so the seed grid is now **12 seeds** (2 lepton signs + 5 signs × 2 baryon families), superseding "the six seeds". The charge class of every seed is now known analytically BEFORE relaxation; the run's winding/degree instrument still verifies each relaxed 3D state, since the 2D cross-section winding is a pre-rotation diagnostic, not the 3D degree read.

Two framing changes ride the same mail. **The selection principle**: these are baryon CANDIDACIES (the author links the full baryon list), not necessarily proton/neutron; the identity is assigned by ranking, proton-analog = the LIGHTEST CHARGED relaxed state, neutron-analog = the LIGHTEST NEUTRAL relaxed state, and the winding/degree instrument measures each state's charge class rather than assuming it. Heavier relaxed states are candidate excited baryons, not failures. The 2026-07-29 charge flip is this principle already doing work at seed level: the labels swapped and the protocol did not change. **The direct ask**: "Please take a look at them, search especially for the lightest baryons, and I will think further tomorrow", so the run has an explicit author request attached (relax the 12 seeds, rank per charge class) and further input keeps arriving: the seed set is live, not final.

**Net effect on the unblocked verdict**: rungs 1 AND 2 are now both seed-unblocked with author-authored analytic seeds (the neutron's field-rotation deferral is superseded). Operational decode (ours, not the author's words): the analytic cross-sections fill the whole half-plane, so the automatic protocol's constrained-fill step is unnecessary for rungs 1-2 and the pin-mask minimizer extension demotes to OPTIONAL (still useful for transcribing author sketches without formulas, e.g. in [M5.22.9](m5_22_9_task_details.md)). Cylindrical symmetry is an explicit, author-flagged assumption: relax with full 3D freedom and report whether the rotated symmetry survives.

## The M5.27 hand-off: the sharpened 4×4 problem statement (2026-07-24)

[M5.27](m5_27_task_details.md) (the background-scalar entrainment pilot, the third fork branch of the 4×4 case) closed with a structural null that SHARPENS the 4×4 blocker this task exists to recruit help for. This section is the carrier: the recruitment brief and any 4×4 side-work run under this task inherit it, so the Lagrangian-level scope has a live trigger and does not sit only in a change-log.

| Hand-off item | Content |
| --- | --- |
| The measured constraint | A background scalar coupled through the spectral targets (`g → g + κχ`) cannot entrain or power the clock: the drive force commutes with `M` on block-diagonal states (commutator 4.5e-21, machine zero), so it owns the eigenVALUES exactly (time-time eigenvalue tracks the drive with slope 1.000) and exerts zero torque on the eigenFRAME that carries the clock. Tongue map NULL 40/40; audit 5/5 ([`../findings/m5_27_note.md`](../findings/m5_27_note.md)) |
| The load-bearing channel | The mixed (0,i) block is the ONLY coupling channel to the clock (commutator 1.38e-2 once populated), is dynamically invariant from block-diagonal data (the drive cannot excite it), and is UNSTABLE when hand-seeded (non-finite by t ≈ 15, the [M5.21.3](m5_21_3_task_details.md) all-negative time-mixing curvatures made dynamical) |
| The consequence | Any working time-sector coupling must be NON-COMMUTING with `M`, which means changing the Lagrangian (the mixed block stabilized, or derivative/frame couplings), NOT designing a better background drive. The [M5.27.x series is PARKED](../m5_roadmap.md) on exactly this ground; the anti-recipe row is in [`../m5_theory_canonical.md § 6`](../m5_theory_canonical.md) |
| The candidate cure | Skyrme-style extra Lagrangian terms (author-floated 2026-07-22, this task's § above), screened BEFORE any spend by the [M5.20.4](m5_20_4_task_details.md) arm-C lemmas + the [M5.20.5](m5_20_5_task_details.md) γ = −1 statics kill |
| How this task triggers it | (a) The convince-others brief now states the concrete blocker in one sentence: the time sector needs a non-commuting coupling and the only known channel is unstable, soliton-specialist territory. (b) The M5.27 negative rides the next author-channel outbound batch (user-gated, queued in the [M5.27 review](m5_27_task_details.md)). (c) Any Skyrme-terms screening or mixed-block stabilization attempt during this task's arc runs as an M5.22-adjacent sub-task consuming this table as its problem statement |

## FINDINGS (2026-07-30, the census run)

Full record: [`findings/m5_22_note.md`](../findings/m5_22_note.md) (equations, gates, forks, census table, reads, Q38, audit). Data: [`data/m5_22_census.json`](../data/m5_22_census.json) (36 rows) + [`data/m5_22_audit.json`](../data/m5_22_audit.json). Adversarial audit: **9/9 claims confirmed, 0 refuted** (independent script, own methods).

| Headline | Result |
| --- | --- |
| The baryon analogs EXIST | Proton-analog = the fractional-family s = −1/2 state: E = 8.250 (n = 48, residual force 5e-7), \|Q\| = 1, a central vortex column + one equatorial ring. Neutron-analog = the fractional-family s = −1 state: E = 12.719 at exact f_tol, Q = 0, scale-stable to 0.1%, a column + TWO rings. Both return exactly after a 2% perturbation kick |
| Mass ordering ✅ (read ii) | Neutron-analog heavier: ratio 1.54 at δ = 0.3, 1.55 at δ = 0.2 (δ-robust, both f_tol); qualitative grade only |
| Same-charge hierarchy ✅ (read v) | Proton-analog (8.250) heavier than the protocol-matched lepton state (6.253), ratio 1.32 at toy parameters |
| Q37 answered ✅ (read vi) | Far-field degrees of the proton-analog and the lepton EXACTLY equal integers (solid-angle instrument, audit-confirmed), with DIFFERENT knot content: the Sulich challenge simulated |
| HONEST NEGATIVE: the author's neutron candidate | The central-π-step s = −1/2 seed DISSOLVES (E slides 2.34 → 2.06 over 24000 iterations, cross-stencil ratio ~5; free BC → vacuum): not protected at toy parameters. The selection principle reassigned the neutron title to the census's own protected neutral state |
| The \|s\| = 1 escape | Integer-wound cross-sections (2D charge ±2, ±3) do NOT carry their charge into 3D (exact analytic degree 0, audit-confirmed): no charge-2 candidates from this seed set; the escaped P −1 seed is HOW the census found its neutron-analog |
| Q38 first data point | In-baryon quark-shift cost 0.123/lattice unit vs the M5.21.4 string tension 6.2-7.0: a 50-57× separation: the two linear terms are NOT the same object at this rung |
| Side-finding → [Q39](../m5_question_tracker.md#q39-detail) | The lepton seed relaxes 8.6% BELOW the certified M5.21.2b electron state at identical config: the lightest known charged state moved |
| Pre-run question #3 (cores vs shells) answered | Empirically: the neutron-analog is a compound of one axis column + two vortex rings with quadrupolar charge lobes netting zero, not a simple core/shell onion; the Wilson profile signature is not cleanly resolved at these sizes |

![proton-analog, n = 48](../plots/m5_22_slice_P-0.5_plane_sc6_n48_pinned_d0.3.png)

![neutron-analog, n = 48](../plots/m5_22_slice_P-1_plane_sc6_n48_pinned_d0.3.png)

![the seed gallery](../plots/m5_22_seed_gallery.png)

![the census panel, n = 32](../plots/m5_22_census_panel.png)

## TASK REVIEW (2026-07-30)

**Task Duration:** 02:52 (from 11:28 to 14:20 EDT)
**Usage Cap Triggered:** NO

**Results**: ✅ measured: protected proton-analog (fractional s = -1/2: E = 8.250 n = 48, \|Q\| = 1, residual force 5e-7, column + one ring) and neutron-analog (fractional s = -1: E = 12.719 exact f_tol, Q = 0, column + two rings), both perturbation-stable; mass ordering 1.54 (δ-robust: 1.55 at δ = 0.2); same-charge hierarchy 1.32; [Q37](../m5_question_tracker.md#q37-detail) RESOLVED (exactly equal integer degrees, different knots); the \|s\| = 1 escape (2D winding does not lift, exact analytic 0); [Q38](../m5_question_tracker.md#q38-detail) first data point (quark-shift 0.123/lattice unit vs tension 6.2-7.0: 50-57× apart); ✅ HONEST NEGATIVE: the central-π-step neutron candidate dissolves (24000 it, free BC → vacuum), the selection principle reassigned the title; 🔶 partial: core/shell profile (compound, Wilson signature unresolved); ✅ adversarial audit 9/9 confirmed, 0 refuted, two claims sharpened. New [Q39](../m5_question_tracker.md#q39-detail) (the lepton seed lands 8.6% below the certified 2b electron state).

**Issues / blockers**: P+1/2 keeps xr ≈ 2 (reported as the conjugate partner, splitting 1.3% = discretization); N+0 / N+1/2 under-resolved, not-citable. None blocking.

**Deviations from plan**: n = 32 census + n = 48 confirmations replaced "all 12 at n = 48" (seeds converged ~5× faster than calibrated; headliners got full n = 48 + extensions).

**Action needed**: the census checkpoint note to the author (the direct ask answered; [M5.22.1](m5_22_1_task_details.md) gated on the state-identity review); the two MODELS.md baryon cells updated at this close; [Q39](../m5_question_tracker.md#q39-detail) is a self-contained machine-checkable follow-up.

**Findings**: The toy baryon census delivered the author's two headline predictions at qualitative grade: protected proton and neutron analog states exist on the certified 3×3 stack, the neutral one is heavier (ratio 1.54, δ-robust), and the proton carries exactly the lepton's quantized charge with different knot content. The one honest negative is structural: the central-π-step neutron seed is not protected at toy parameters; the census's own selection principle recovered a protected neutral state from the escaped s = -1 fractional seed instead.

**Research docs created / updated**: [`findings/m5_22_note.md`](../findings/m5_22_note.md) (the census note) · this doc § FINDINGS · [`m5_question_tracker.md`](../m5_question_tracker.md) (Q37 ✅, Q38 data point, Q39 new) · scripts [`m5_22_a_seeds.py`](../scripts/m5_22_a_seeds.py) · [`m5_22_b_census.py`](../scripts/m5_22_b_census.py) · [`m5_22_c_rank.py`](../scripts/m5_22_c_rank.py) · [`m5_22_d_slices.py`](../scripts/m5_22_d_slices.py) · [`m5_22_e_audit.py`](../scripts/m5_22_e_audit.py) · data [`m5_22_census.json`](../data/m5_22_census.json) + [`m5_22_audit.json`](../data/m5_22_audit.json) + [`_DATASETS.md`](../data/_DATASETS.md) · plots: [`m5_22_seed_gallery.png`](../plots/m5_22_seed_gallery.png) · [`m5_22_census_panel.png`](../plots/m5_22_census_panel.png) · the state slices (`m5_22_slice_*.png`)

## CHECKPOINT OUTCOME (the author's 2026-07-30 reply)

The census checkpoint closed with the author's same-evening reply (full capture + decode: [`m5_22_convo.md § 2026-07-30`](m5_22_convo.md)): the census publicly endorsed with the audience widened to the vortex-knot lab authors and the note cited as the "Fable audit"; the proton-analog identity CONFIRMED, with the fractional structure read as quarks-by-vortex-interaction; the neutral state REREAD as candidate DINEUTRON (ring count = baryon number, matching the measured column + two rings); the electric-charge sign convention fixed as reversed (note § 7 addendum); and four work items routed: the δ = 0.1 ladder and the kick-apart identity probe into [M5.22.1](m5_22_1_task_details.md) opening moves, the analytic deuteron seed announced (pending arrival), and the fusion/collision route into [M5.22.9](m5_22_9_task_details.md) rung F.
