# M5.21.8: verify + confront the analytical twisting massive hedgehog

**Status**: ✅ CLOSED 2026-07-19 (review approved same day; § TASK REVIEW below). From the author's 05:09 group message ([`m5_21_convo.md § 2026-07-19 05:09`](m5_21_convo.md)). The author's direct ask: "Maybe please ask Fable to verify it, and if agreeing analyze why its simulations lead to infinites for electron ... I suspect due to use of inadequate g, delta parameters (?)". Source: [`../../theory/duda_2026-07-19_3p1D_analytical_twisting_massive_hedgehog.pdf`](../../theory/duda_2026-07-19_3p1D_analytical_twisting_massive_hedgehog.pdf).

## TASK PLANNING (2026-07-19)

**Verification strategy (P0)**: independent pipeline, not transcription: build the ansatz with closed-form group elements (Rodrigues rotations, cosh/sinh boosts, no symbolic matrix-exp), exact derivatives, high-precision quadrature for the time average and cone-regularized spherical integral; verify the author's DOWNSTREAM claims (the m* formula symbolically/numerically for several g; the transcribable Hm and ω-min formulas; both divergences; the −∞ regions of the author's own Minimize) rather than the author's long intermediate expressions. The author's approximations (single vortex, dropped cone divergence, rigid eigenvalue-pinned profile) are read back as stated scope, not defects.

**The key structural facts to measure, pre-registered**:

| # | Fact | Arm |
| --- | --- | --- |
| 1 | The author's H includes TIME-derivative commutators with component-split ± signs; our stack separates E_stat (spatial) + E_kin (generator flow). The bridge must express one in the other per (i, j) block on the same ansatz | P1 |
| 2 | The author's kinetic sector at constant ω is IR-divergent (∫ ω² dr); our M5.21.3 kinetic was envelope-localized (IR-finite by construction): the two setups ask DIFFERENT ω questions; the lattice box makes the author's finite but box-dependent (measure kin(L)) | P1/P2 |
| 3 | The author's ansatz adds the boost-hedgehog dressing Qb(m), absent from every M5.21.3 run: E_stat(m) may have the finite m* minimum the author computes even on our functional | P2 |
| 4 | The author's vacuum d = diag(g, 1, δ, 0) = our s = −1 branch (−sg = +g); run s = +1 as robustness | P2 |
| 5 | The author's g-critique: m*(g) ~ 1/g and the boost channel shrinking at large g: lattice-measurable NOW (g-ladder) without any 1e10 run | P3 |

| Phase | DoD item |
| --- | --- |
| P0 | verdict per claim (m*, Hm, ω-min value + regions, divergence 1 (cone), divergence 2 (r-integral)): CONFIRMED / REFUTED / PARTIAL with our own numbers |
| P1 | the per-block relation table (the author's H blocks vs our η-trace blocks) + the IR statement measured |
| P2 | E(m, ω) surface on the lattice (rigid) at (g = 8, δ = 0.3) both branches + stationary-point verdict + profile-relaxation survival of any minimum found |
| P3 | m*_lattice(g) and kin_boost(g) across the g-ladder + the extrapolation statement (feeds [Q33](../m5_question_tracker.md#q33-detail)) |
| P4 | the premise-correction paragraph (no infinities measured: bounded slope, no stationary point) + the consolidated outbound draft (terminal-only) carrying the held M5.21.6 note |
| All | eager checkpoints; independent adversarial audit before anything author-facing; doc checker clean; review in terminal + mobile ping |

**Blindspot pass**: machine-checkable = every P0-P3 number. Author-gated = the author's Hamiltonian sign convention as THE intended energy (we measure the bridge, the author owns the choice); the regularization the author deferred ("not certain about boost radius dependence"). Nature-gated = physical parameter values. Risk pre-registered: SymPy blowup on full symbolics → the numeric-sampling fallback is the plan of record, symbolic only where cheap (the m* solve).

## Scope (stub level)

| Phase | Content | Notes |
| --- | --- | --- |
| P0 symbolic reproduction | Reproduce the notebook end-to-end with our own CAS route (SymPy): the ansatz M = Qb(m)·Qh(ωt)·d·(...)ᵀ, F_ij = [∂_i M, ∂_j M]_ξ, the author's component-split H, the time average, **m* = ½ ln((1+g)/(g−1))**, the Hm density, the two divergences (vortex cone; constant-ω r-integral), the piecewise ω-minimum incl. its −∞ regions | The verification the author asked for, done adversarially (the author's own −∞ branches and dropped divergences are part of the honest read-back) |
| P1 convention bridge | The author's H (component-wise ± split of F entries, ξ-sandwich) vs our η-trace u_η on the SAME ansatz: identical / proportional / different, term-by-term; locate the author's twist × spatial cross terms (the claimed ω-stabilizers) in our functional | Decides whether the M5.21.3 verdict and the author's benchmark even talk about the same energy |
| P2 lattice confrontation | Implement the Qb(m) boost-hedgehog dressing in the 4D stack; scan E(m, ω) at (g = 8, δ = 0.3): rigid (the author's approximation) AND profile-relaxed; finite (m*, ω*) or not; explain any difference (the author's rigid eigenvalue-pinned family + dropped divergences vs our relaxing lattice; IR: the author's constant-ω kinetic diverges, our envelope-localized kinetic is IR-finite by construction) | The M5.21.3 ω-ladder had NO boost-hedgehog dressing; the author's m-parameter is a genuinely new direction in the ansatz space |
| P3 the g-ladder scaling test | The author's critique measured: kin_boost(g) and m*(g) along g = 8/16/32/64...: does the negative boost channel shrink ~ 1/g as the author's m* formula implies? Extrapolation statement toward g ~ 1e10 | The [Q33](../m5_question_tracker.md#q33-detail) bridge in action; the cheapest honest answer to "resolved at real g" |
| P4 the premise correction + outbound | The reply paragraph: M5.21.3 measured NO infinities (bounded shallow slope, no stationary point; the divergence era was the pre-symmetrization signature dive, cured); what is missing is a turning point, and where the author's cross-terms/IR structure could supply one; package with the held [M5.21.6](m5_21_6_task_details.md) note | The combined outbound closes both threads |

Series rules: independent adversarial audit before anything author-facing; method-note grade; [Q35](../m5_question_tracker.md#q35-detail) (negative-Hamiltonian literature) can ride as a P0 side-read.

**Staged follow-on (the author's Larmor protocol, spelled out 2026-07-19)**: "introduce constant external magnetic field by constant field derivative in its direction, and temporal field derivative - should lead to electron precession with frequency proportional to magnetic field strength"; runs on whichever stable-oscillation state survives this task (a true (m*, ω*) minimum if P2 finds one, else the fixed-J/enforced state, which the author now explicitly allows: "maybe enforced if numerical problems remain"). The linear-in-B precession slope is the acceptance observable.

**Reference map for P1**: the author's F_μν structural dictionary ([`../../theory/duda_2026-07-18_Fmunu_decomposition.png`](../../theory/duda_2026-07-18_Fmunu_decomposition.png), decoded in [`m5_21_convo.md § 2026-07-18 19:44`](m5_21_convo.md)): EM = tilt-tilt (R¹ + g²Γ̃¹), QM = tilt-twist (δR² − g²Γ̃²), gravity = boosts, "clock propulsion?" on the twist-boost cross blocks; the bridge should express our η-trace u in exactly these blocks.

**Gated by**: user "go" (all inputs delivered; no 4D re-run needed for P0/P1).

## FINDINGS (2026-07-19)

Full method-note record: [`../findings/m5_21_8_note.md`](../findings/m5_21_8_note.md). The headline rows:

| # | Finding | Where |
| --- | --- | --- |
| 1 | **The convention bridge: the author's component-split Hamiltonian = EXACTLY ½ × our η-trace, block by block, at every sample: the SAME functional.** No sign-convention gap between the author's benchmark and the M5.21.3 verdict | note § 3 |
| 2 | **The author's m* = ½ ln((1+g)/(g−1)) CONFIRMED** analytically (g = 3: 0.3456 vs 0.3466; g = 8: 0.1231 vs 0.1257, cone-robust) AND on the lattice: twin even minima at \|m*\| = 0.1027, E_u = 0.968 (65× below the undressed hedgehog), h-robust; **the g-ladder tracks the author's 1/g law at a stable 0.82-0.84× across g = 8-64**: the gravitational-mass mechanism is REAL on our audited functional | note § 4, § 5, § 5b |
| 3 | **The author's "ω minimization → finite" REFUTED AS STATED, by the author's own formulas + our pipeline + the lattice**: the ω² coefficient of the author's Hm is negative for ALL g > 1, δ ≤ 1.9 (incl. (1e10, 1e-10): the author's own Minimize routes toy AND physical parameters to −∞); where positive (g < 1 or δ ≥ 2.5, outside the construction) the minimum sits at ω* = 0. **No finite NONZERO ω exists anywhere in the notebook.** On the lattice at the energy minimum kin = +75.5 > 0 → ω* = 0: finite, static, no spontaneous clock | note § 4, § 5 |
| 4 | The author's g-critique SPLIT verdict: RIGHT for the m-sector (toy boosts are rapidity-gigantic; m* ~ 1/g → 1e-10 at real g), does NOT heal the free-ω sector (the sign map is g-independent for g > 1) | note § 4 |
| 5 | The dropped cone divergence is NEGATIVE at the author's chosen branch (an ω-independent unbounded-below channel at the vortex core): "needs regularization" is load-bearing | note § 4 |
| 6 | The constant-ω kinetic is IR-extensive (analytic: per-r density → negative const; lattice: kin ∝ L, 1.285 vs 4/3): the author's margin note ("ω needs to reduce with distance?") is structurally right; an ω(r)-profiled clock is a different variational family, designed-not-run | note § 4, § 7 |
| 7 | The rigid family is NOT directly relaxable at n = 32 (FIRE non-finite at dt0 down to 1e-4, 3 attempts): the author's own "singularities should be regularized" caveat, measured | note § 5 |
| 8 | The premise correction stands: M5.21.3 measured NO infinities (bounded slope, no stationary point); the author's own analytics now agree (ω* = 0 or −∞); the constructive convergence = constraint-carried ω (fixed-J), which the author's "maybe enforced" already sanctions | note § 6 |

![panel](../plots/m5_21_8_panel.png)

## DEVIATIONS LOG

| # | Deviation | Why |
| --- | --- | --- |
| 1 | Checkpoint interim claims corrected mid-run: "the lattice breaks ±m degeneracy" (compared different \|m\|) and "kin < 0 across the family" (read from the \|m\| > 0.125 tail + a box probe at m = +0.1256, past the sign flip) were both WRONG; the full curve is exactly even with kin > 0 in the minimum band | caught on the full-curve read; the note carries only the corrected picture |
| 2 | relax unpack bug (4D fire returns 2 values, not 3, unlike the 2b fire) cost one crashed run | fixed; relax then measured non-finite at 3 dt rungs = the § 5 finding |
| 3 | g = 16/32/64 coarse m-grids under-resolved their minima (grid 0.025 vs features ≤ 0.031) | fine grids added for all rungs; the coarse artifact documented |

## TASK REVIEW (2026-07-19)

Task Duration: 1:18 (from 10:08 to 11:26)
Usage Cap Triggered: NO (finished before the 2:30pm reset; the armed ping parked unfired)

| Result | Status |
| --- | --- |
| The convention bridge: the author's Hamiltonian = ½ × our η-trace exactly, generic theorem (audit-sharpened) | ✅ |
| The author's m* law: essentially exact (0.009%); lattice twin minima track 1/g at 0.82-0.84× across g = 8-64; gravitational-mass mechanism real | ✅ |
| The author's "ω → finite" claim: refuted as stated by the author's own formulas (ω* = 0 or −∞ everywhere, incl. (1e10, 1e-10)); lattice agrees (kin(m*) = +75.5 → ω* = 0) | ⚠️ honest headline |
| g-critique split verdict; cone term negatively divergent; kinetic IR-extensive (kin ∝ L); rigid family not lattice-relaxable (3 dt rungs) | ✅ scope findings |
| Independent adversarial audit 7/7 CONFIRMED, 0 refuted, 2 nuances adopted (one crediting the author beyond our first pass) | ✅ |
| Consolidated outbound draft presented (carries the held M5.21.6 results); author-formality standard applied repo-wide to the touched docs | ✅ |

Issues: the ω(r)-profiled clock family and the regularized-core variant are designed-not-run (both author-gated on the core profile); the g = 50 analytic scan stays under-resolved (flagged, non-load-bearing).

Action at close: roadmap row → Done (appended at END); the author-formality rule ("the author", never "his/he"; "your" only in direct messages) adopted as the standing standard for any model author and saved to memory.

**Findings**: The author's analytical notebook verifies essentially exactly on the gravitational-mass sector (the m* law is real on the audited lattice, tracking 1/g at 0.83×) and refutes-as-stated on the clock sector by the author's own formulas (no finite nonzero ω exists at any physical parameters), putting the analytics and the lattice in full agreement that the clock must come from regularization, an ω(r) profile, or a carried constraint rather than free minimization.

**Research docs created/updated**: [`m5_21_8_task_details.md`](m5_21_8_task_details.md) · [`../findings/m5_21_8_note.md`](../findings/m5_21_8_note.md) · [`../m5_roadmap.md`](../m5_roadmap.md) · scripts `m5_21_8_a_verify.py` / `m5_21_8_b_lattice.py` / `m5_21_8_d_panel.py` / `m5_21_8_audit_check.py` · data `m5_21_8_p0.json` + lattice/audit JSONs · plots `m5_21_8_panel.png`
