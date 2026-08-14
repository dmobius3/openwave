# M5.21.16: The sign-reversed boost Hamiltonian: (Γ̃)² flip, Coulomb vs Newton

## TASK PLANNING

**Scope**: run the author's 2026-08-13 quest (email replying to the [M5.21.15](m5_21_15_task_details.md) send; notebook [`theory/duda_fmunu_4d_hamiltonian_imaginary_2026_08_13.pdf`](../../theory/)): "modify Lagrangian and Hamiltonian for negative terms of (Γ̃)², or imaginary g and conjugation - checking Coulomb vs Newton with signs, and if electron could get finite gravitational mass and frequency by energy minimization?" The author's own symbolic finding (the notebook): the ω² energy coefficient at leading order g→1/δ is exactly `−2Σ(Γ̃)²` (negative-definite, confirming our measured boost-channel runaway), and his variant with imaginary g plus a conjugated second F makes it positive-definite.

**Definition of done** (five arms):

| Arm | Deliverable | Testable criterion |
| --- | --- | --- |
| A symbolic | `m5_21_16_a_symbolic.py`: reproduce the notebook in sympy (his exact conventions: Γ_μ 4×4, ξ = diag(−1,1,1,1), coms, d = diag(g,1,δ,0), H with spatial-internal² minus time-internal², series in ω = Γ_0^1) for baseline, variant A (imaginary Γ̃), variant B (imaginary g + conjugation); plus the BRIDGE LEMMA: our `inner_eta → Frobenius` substitution reproduces the sign-reversed (Γ̃)² Hamiltonian | zero symbolic residuals vs the notebook's three ω² coefficients; bridge lemma exact |
| B field kin | `m5_21_16_b_field.py`: the FLIP functional (same comm_η, same V4, Frobenius contraction on curvature) on the certified M5.21.3 stack; per-channel kin (clock_local, plane_1d, rot_z, rot_x, boost_z, boost_x) baseline vs flip at g = 32, δ = 0.3, s = −1; invariance audit (SO(3) kept, SO(1,3) expected broken, SO(4) expected gained); vacuum + complex-step gradient gates under flip; the 3×3-embed identity (static block fields: flip == baseline, so the Coulomb sector is untouched by construction) | boost-channel kin sign under flip measured; all gates green; the invariance table measured not assumed |
| C boundedness + fixed-J | `m5_21_16_c_fixedj.py`: the M5.21.14 dressing family under flip: T1_kin sign, the unbounded-below E_corr well (the [0.02, 0.05] guard bracket) re-probed with the family bound OPENED; fixed-J E(ω) minimum under flip: ω*, E sign, guard-width sensitivity | does the guard become moot (bounded without it); ω* finite and E sign stable vs family width |
| D two-body signs | `m5_21_16_d_twobody.py`: seed-level two-center energy vs separation: charge channel (3×3 hedgehog pair) and mass channel (two boost-dressed clocks), sign of dE_int/dd, baseline vs flip | the 2×2 sign table (charge/mass × baseline/flip) measured at seed level; Coulomb repulsion preserved, Newton sign read |
| E audit | `m5_21_16_e_audit.py`: independent recompute (no imports from A-D): finite-difference ω² coefficients vs the symbolic forms on random single points; mutation tests (a wrong sign in the flip must redden a gate) | audit green with discriminating gates |

**Gating**: [M5.21.15](m5_21_15_task_details.md) (instruments + the envelope-concavity theorem). The theorem SURVIVES the flip (E affine in ω² per configuration is metric-independent), so the honest success criterion for "preferred frequency" is: vacuum prefers ω = 0 free (all kin ≥ 0 under flip) and the electron's finite ω comes from the fixed-J branch; free minimization can still never select interior nonzero ω on this functional class.

**Blindspot pass** (unfamiliar territory: signature surgery):

| Blindspot | Route |
| --- | --- |
| The flip breaks SO(1,3) invariance of the energy (Frobenius is compact-group invariant): Duda's "time imaginary" is a Euclideanization, the internal symmetry becomes SO(4)-like | machine-checkable: measure invariance under SO(1,3) and SO(4) conjugations in arm B; report as a structural consequence, do not hide it |
| His notebook H is a single-point symbolic density (field-independent Γ symbols); our functional adds stencils, boundaries, V4 | the bridge lemma (arm A) is at the density level; field gates (arm B) certify the lattice side |
| Variant B (imaginary g + conjugation) complexifies M; not directly implementable on the real-symmetric stack | symbolic-only in this task; arm A checks whether its leading-order sign pattern equals the plain flip; field implementation = follow-up if they differ |
| Seed-level two-body signs are not force curves (no relaxation) | label honestly; the full 1/r² curve is out of scope here |
| "Hamiltonian" vs "Lagrangian" naming in the notebook | the parabola-in-ω analysis fixes the observable: the ω² coefficient maps to our kin up to a positive factor; verified numerically in the bridge lemma |

**Research body**: findings note [`findings/m5_21_16_note.md`](../findings/m5_21_16_note.md) (method-note shape); scripts/data/plots under `research/` with the `m5_21_16_` prefix; checkpoints in `research/checkpoints/m5_21_16_progress.md`; convo record for the author's email: [`m5_21_convo.md`](m5_21_convo.md) (the email is the M5.21.15 thread reply; this task links it).

**Preconditions**: the author's PDF filed in `theory/` ✅; the certified stack (`m5_21_3_a_4d.py`) and the analytic dressed family (`m5_21_14_c_minimize.py`) present ✅.

## TASK REVIEW (2026-08-14)

**Task Duration:** 00:21 (from 15:56 EDT go to 16:17 EDT review)
**Usage Cap Triggered:** NO

**Results**

| Arm | Result |
| --- | --- |
| A symbolic | ✅ measured: the author's notebook reproduces EXACTLY (independent sympy): baseline ω² lead `−2Σ(Γ̃)²`; variant A `+2Σ`; variant B positive-definite with extra Γ̃₀ terms. Bridge lemma EXACT: notebook H == `Σ_{μ<ν}⟨F,F⟩_η`; the field-level flip (Frobenius contraction) == variant A exactly at leading order |
| B field | ✅ measured: boost channels kin flips sign exactly (pure time-row curvature, spatial part 0); all channels kin ≥ 0 under flip; charge sector flip == η identically (rel 0.0); COST: SO(1,3) energy invariance breaks 26% at b = 0.25 (SO(3) exact) |
| C boundedness + fixed-J | ✅ measured: wide-family wells η −2354 (all 6 starts < −1600) vs flip −98.8 clustered; sawtooth runaway closes (−2339 → +21950): the M5.21.14 guard MOOT under flip. Fixed-J electron: interior optimum every J, E positive (440/457/735), ω* finite (0.046/0.185/0.742) |
| D two-body | 🔶 partial: Coulomb closed by construction; Newton sign NOT decidable at seed level, and why is measured: vacuum boost response QUARTIC in amplitude (ratio 16.0 = 2⁴, no quadratic kernel; matches the M5.20.2 purely-quartic L), two-bump overlap metric-insensitive (0.09%). Routed to the relaxed two-defect successor, gated on the canonical two-center construction |
| E audit | ✅ 11/11 gates, fully independent (own matrices/stencils/quadrature); mutations redden; η well recomputes −2308 vs claim −2354 (3.6%, independent quadrature) |

**Issues / blockers**: none blocking. Two audit-design defects found and fixed inside the audit itself (float cancellation at δ = 1e-4; over-tight tolerance across intentionally different stencils).

**Deviations from plan**: arm D's planned two boost-dressed clocks replaced by the two-boost-bump vacuum probe + the quartic-kernel measurement, after the seed-level construction proved incapable of carrying the sign question (measured, not assumed).

**Action needed**: [Q48](../m5_question_tracker.md#q48-detail) (author-gated: canonical variant A vs B; is losing Lorentz invariance of the energy accepted) rides the next author note; successor candidate = the relaxed two-defect Newton-sign task; variant B at field level if the author names it canonical.

**Model-doc sweep**: canonical registry updated (§ 5 candidate-variant row + § 6 anti-recipe row); briefing unchanged (the flip is not adopted, nothing the card states moved).

**Findings**: The author's sign diagnosis is right and his fix works where he predicted: the (Γ̃)² reversal (== his variant A, implemented as the Frobenius contraction on the certified stack) makes every kinetic channel non-negative, closes the M5.21.14 unbounded dressing channel with no guard, and gives the fixed-J electron a finite frequency at positive energy, while leaving the Coulomb sector bit-identical. The measured price is Lorentz invariance of the energy (26% break at b = 0.25), and the Newton-vs-Coulomb sign cannot be read from any vacuum-level probe (the vacuum boost response is quartic): it needs relaxed two-defect states on the canonical two-center construction only the author can specify.

**Research docs created / updated**: [findings/m5_21_16_note.md](../findings/m5_21_16_note.md) (method note) · this task record · [m5_21_convo.md](m5_21_convo.md) (2026-08-13 entry) · [m5_question_tracker.md Q48](../m5_question_tracker.md#q48-detail) · [m5_roadmap.md](../m5_roadmap.md) · scripts `m5_21_16_{a..e}_*.py` · data 5 tracked JSONs · plots [m5_21_16_panel.png](../plots/m5_21_16_panel.png), [m5_21_16_twobody.png](../plots/m5_21_16_twobody.png)
