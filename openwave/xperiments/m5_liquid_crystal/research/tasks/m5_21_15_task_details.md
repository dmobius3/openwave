# M5.21.15: the coupled ω-minimum (electron analog with angular momentum)

> Task **M5.21.15** (M5 / Liquid-Crystal model). Status: 🚧 **STAGED, AWAITING USER "GO"**
> (2026-08-10) · Roadmap: [`../m5_roadmap.md`](../m5_roadmap.md) · Staged from the author's
> 1:1 reply to the [M5.21.14](m5_21_14_task_details.md) close-out (decode:
> [`m5_21_convo.md § 2026-08-10`](m5_21_convo.md)). Series: M5.21.x, the electron hunt.
> The reply to the author is HELD until this task runs (user call, 2026-08-10).

This doc is the task's full record: planning, then findings at the run.

## PLANNING

### Why this task exists

The author's 2026-08-10 reply accepted the [M5.21.14](m5_21_14_task_details.md) close ("looks good", the dressed minimum read as "clear energy minimum suggesting gravitational mass") and named the next quest: "(coupled) nonzero omega from energy minimization" for the electron, pointing at our own [M5.21.3 record](../findings/m5_21_3_note.md) with the read "it stopped minimization process, but visually there should be minimum for positive energy and nonzero omega (?) - worth finding this energy minimum: for electron analog with angular momentum".

The two halves already on disk make this a well-posed measurement rather than an open hunt:

| Baseline fact | Record |
| --- | --- |
| Free 4×4 descent lands NO stationary ω at toy parameters: E\*(ω) monotone decreasing (boost channel, shallow, profile-decoupled), every rotation channel kin > 0; the P3 rungs are depth-bounded at max_iter (the "stopped" the author sees), and the equal-depth control subtraction certifies the whole ω-advantage as the quadratic kinetic margin | [M5.21.3 § 6](../findings/m5_21_3_note.md) |
| The constraint-carried route works: fixed-J states exist and hold with exact clock thermodynamics dE/dJ = ω\* | [M5.21.9](m5_21_9_task_details.md) |
| The (1/g) dressing drives the ω² coefficient DOWN in the bulk (T1_kin = −8Σ\|Ṁ3v_i\|²·ω², negative-semidefinite; per-unit-radius slopes base +13.5 vs dressing −19.6), with the certified 48-box exactly at the flip threshold; favorability-not-onset, and T1 alone is unbounded below so the dressing must enter GUARDED | [M5.21.14](../findings/m5_21_14_note.md) |

What is genuinely NEW in the ask: whether the landscape holds a free interior MINIMUM at positive energy and nonzero ω once ω is COUPLED to the responding profile (texture + dressing b(r)), rather than an ω-slope read at fixed profile (M5.21.3) or a constraint-carried construction (M5.21.9). At fixed profile E(ω) is exactly quadratic (the M5.21.3 decoupling), so any interior minimum can only come from the coupling. This is also exactly the onset question M5.21.14's caveat left open: if a stabilized minimum exists, "favorability" upgrades to "onset".

### The author's ask (2026-08-10, pinned before any computation)

| Item | Content |
| --- | --- |
| The object | the energy minimum at POSITIVE energy and NONZERO ω, "for electron analog with angular momentum", out of energy minimization with ω coupled |
| The reference record | the [M5.21.3](../findings/m5_21_3_note.md) ω-ladder ("it stopped minimization process"); the author's "(?)" marks it a genuine open question, not a claim |
| The context signals | the reply reads the M5.21.14 dressed minimum as gravitational mass; the standing negative-Hamiltonian frame ([Q35](../m5_question_tracker.md#q35-detail)): gravity's negative terms are what "allow nonzero time derivatives by energy minimization (angular momenta of the electron...)" |

### The plan (staged; firmed at go)

**The instrument.** The certified 4D stack ([`../scripts/m5_21_3_a_4d.py`](../scripts/m5_21_3_a_4d.py): e_parts + kin) plus the M5.21.14 exact-dressing machinery ([`../scripts/m5_21_14_c_minimize.py`](../scripts/m5_21_14_c_minimize.py): `ExactCorr`, the smooth guarded b(r) families, the corner hygiene). No new physics enters: the run measures the landscape the two records jointly define.

| Arm | Content |
| --- | --- |
| A1 the undressed retro-gate | reproduce the M5.21.3 baseline on the current stack (E\*(ω) monotone, no stationary ω, rotation kin > 0): the hard gate before anything new, and the honest answer to the author's "(?)" on the UNDRESSED functional |
| A2 the coupled dressed scan | E\*(ω) = min over (texture deformation, guarded b(r)) of the EXACT dressed functional at finite g, ω swept through and past the probed range; the profile re-optimized per rung (the coupling); box-radius as a controlled axis (the flip is BULK: R below/at/above the threshold the 48-box sits on); both s signs |
| A3 the angular-momentum read | the rotation (clock) channels, not only the boost channel: does the dressed functional turn any rotation-sector ω² coefficient negative anywhere in (R, g, profile) space? J measured on every candidate minimum (the ask is an electron analog WITH angular momentum) |
| A4 the fixed-J bridge | the constraint-carried map E(J) with the dressing: ω\* = J/(2·kin_eff) as kin_eff crosses zero; where the free minimum (if any) sits relative to the fixed-J family; dE/dJ = ω\* re-verified on the dressed states |
| A5 the verdict | pre-registered both ways: a positive-energy nonzero-ω minimum EXISTS (report location, depth, J, stability reads) or DOES NOT in the probed ranges (report the exhausted ranges + which structural fact blocks it); either answer is the deliverable and the checkpoint payload |

**Verification**: the A1 retro-gate; the M5.21.14 threshold numbers reproduced before the scan (kin_corr −426.3 vs base +426.5 on the certified box); every minimum candidate re-checked on the lattice instrument at the affordable n; independent adversarial audit (own route) before anything is trusted (cardinal rule).

**Guard discipline** (the M5.21.14 mandatory guard, unanswered fork): the dressing enters ONLY through constrained smooth b(r) families (the exact functional at finite g, never bare T1); the guard-choice question stays open with the author, so the family is declared a PROVISIONAL guard in every output, revisable when the regularization potential's details land ([Q25](../m5_question_tracker.md#q25-detail)).

**Blindspot pass** (run at go): (1) runaway vs minimum: a descending E\*(ω) that never turns is the M5.21.3 negative repeated, not a discovery; the stop rule is a stationary-point bracket or an exhausted declared range, never "deep = good"; (2) depth-bounded relaxation masquerading as a minimum (the author's own "stopped" trap): convergence certified per rung (force norm, not iteration cap) or the rung is labeled contained-not-converged; (3) the R-axis confound: the bulk flip means the verdict can depend on box radius; R is a declared scan axis, not a fixed choice; (4) positive energy: E > 0 checked against the correct vacuum reference on the dressed functional (E_V invariant, but E_u references matter); (5) the two functional readings (η vs Hamiltonian) both carried, as in M5.21.3.

**Definition of done**: A1 gate green; the coupled scan run with per-rung convergence certification; the angular-momentum and fixed-J reads delivered; the pre-registered verdict stated both-ways-honest; audit recorded; method-note-grade findings (`../findings/m5_21_15_note.md`) + the checkpoint outbound drafted (the HELD reply rides on this).

**Artifacts** (all `m5_21_15_` named): `scripts/m5_21_15_a_baseline.py` (A1), `scripts/m5_21_15_b_coupled.py` (A2-A3), `scripts/m5_21_15_c_fixedj.py` (A4), `scripts/m5_21_15_e_audit.py`, `data/m5_21_15_*.json`, `plots/m5_21_15_panel.png`, `findings/m5_21_15_note.md`, checkpoint `checkpoints/m5_21_15_progress.md`.

**Model/effort**: Fable / high (landscape-measurement task on certified machinery; the compute is lattice re-relaxations at moderate n plus continuum quadrature).

### What this task does NOT do

| Non-goal | Why |
| --- | --- |
| Stage the corrected-3×3 ladder | the guard-choice fork is the author's open question; this task runs 4×4 where the dressing is native |
| Claim onset from favorability | the upgrade happens only if a certified stationary minimum is measured; otherwise the M5.21.14 caveat stands verbatim |
| Quote mass ratios or physical ω | B/C remain uncertifiable at N = 48; toy parameters only, the realistic-parameter bridge stays [Q33](../m5_question_tracker.md#q33-detail) |
| Answer the author before the run | the reply is HELD (user call, 2026-08-10); the checkpoint outbound is drafted at close and sent by the user |

**Gated by**: user "go".
