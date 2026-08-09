# M5.21.14: the (1/g) dressing term (symbolic 4×4 boost hedgehog)

> Task **M5.21.14** (M5 / Liquid-Crystal model). Status: ✅ **CLOSED COMPLETE + APPROVED**
> (2026-08-09, review below) · Roadmap: [`../m5_roadmap.md`](../m5_roadmap.md) · Staged
> 2026-08-09 from the author's 1:1 reply to the [M5.21.11](m5_21_11_task_details.md) ladder
> close (decode: [`m5_21_convo.md § 2026-08-09`](m5_21_convo.md)). Series: M5.21.x, the
> electron hunt (lepton-sector instrument work).

This doc is the task's full record: planning, then findings at the run.

## PLANNING

### Why this task exists

The [M5.21.11](m5_21_11_task_details.md) ladder closed route (b) terminally: the 4D boost-dressing
correction is O(1), branch-dependent, and FLAT in g (measured on the ladder endpoints and
retro-read in the [M5.21.8](m5_21_8_task_details.md) family record), so no 3×3-only ladder can
reach physical-regime energies or ratios. The close-out named the required shape of any successor:
the dressing carried INSIDE the ladder. The author's 2026-08-09 reply supplies the concrete first
rung, symbolic rather than numeric, so it needs no lattice, no certified branches, and no new
relaxations.

### The author's recipe (2026-08-09, pinned before any computation)

| Step | Content |
| --- | --- |
| 1 | Work symbolically in the 4×4 case with at least radial dependence, assuming the spherical boost hedgehog `MatrixExp[b(r) {x, y, z} . boostgenerators]` |
| 2 | Find the FIRST NONTRIVIAL TERM in the (1/g) expansion |
| 3 | Include that term in the further-considered 3×3 case |
| The criterion (the author's, pre-stated) | the term "needs to have negative Hamiltonian contribution to get oscillations" |
| The instrument correction it carries | the rigid uniform-m read is too crude (the measured m\* "seems much too large, and the radius dependence seems more complicated"): start from minimization of a GENERAL at-least-radius-dependent b(r), not a single rapidity knob |

### The plan (firmed at go, 2026-08-09)

**The object.** The 4D Hamiltonian of the instrument of record ([`../scripts/m5_21_3_a_4d.py`](../scripts/m5_21_3_a_4d.py): E_u = 4 Σ_{i<j} ⟨F_ij, F_ij⟩_η with F_ij = ∂_iM η ∂_jM − ∂_jM η ∂_iM, η = diag(−1,1,1,1); E_V = the trace-power potential, targets c_p = sg^p + 1 + δ^p), evaluated on the dressed family M_d = Qb M4 Qbᵀ with Qb = exp(b(r) r̂·K) (the general-b(r) form of the rigid `qb_field` in [`../scripts/m5_21_11_d_garm.py`](../scripts/m5_21_11_d_garm.py); closed form I + sinh(b)K + (cosh(b)−1)K², K³ = K). One structural fact fixed at PLAN: Qb is a pointwise η-congruence, so tr((Mη)^p) is invariant and **E_V never fires under any b(r)**; all b-dependence lives in E_u (to be verified symbolically as gate V0).

| Arm | Content |
| --- | --- |
| S1a self-gate | vacuum base D = diag(−sg, 1, δ, 0), CONSTANT b = m: the derived E(m) landscape must reproduce the verified m\* = artanh(1/g) law ([M5.21.8 § 4](../findings/m5_21_8_note.md), 0.009% agreement of record); hard gate before anything new |
| S1b general b(r), vacuum base | sympy: derive the exact radial density h[b, b′; r, g, δ, s], expand in 1/g; the double-expansion ordering is carried BOTH ways (b fixed; b = β/g scaled) and the "first nontrivial term" is read in the scaling the author's m\* ~ 1/g regime selects |
| S1c the texture coupling | replace the vacuum spatial block with the hedgehog texture (the author's static frame, [M5.21.8 § 1](../findings/m5_21_8_note.md)): derive the term COUPLING b(r) to the 3×3 sector: this is the object step 3 hands to the 3×3 functional |
| S1d the oscillation read | add the clock (the ωt rotation of the author's Qh) at quadratic order in ω (the Hm-coefficient structure verified in [M5.21.8 § 4](../findings/m5_21_8_note.md)): derive the dressing correction to the ω² coefficient: the author's criterion reads HERE (negative contribution → oscillations favored) |
| S2 the sign verdict | minimize the derived functional over b(r) (Euler-Lagrange symbolically; scipy on the radial functional if closed form fails), both s = ±1; the pre-registered read: does the first nontrivial term contribute NEGATIVE Hamiltonian density at the minimizer, and does it turn the ω² coefficient negative? Either answer is the deliverable |
| S3 the 3×3 handoff | write the term as an additive correction to the 3×3 functional and state what it retrodicts for (i) the FLAT-in-g gain (the M5.21.11 g-arm + the M5.21.8 family record) and (ii) the author's m\*-too-large read (variational b(r) vs the rigid uniform-m read) |
| Numerical sequel (OUT of scope) | lattice implementation of the corrected 3×3 functional + an instrument certifying B/C: a separate task if S1-S3 land |

**Verification (AI_HYGIENE: nothing trusted on symbolics alone).**

| Gate | Check |
| --- | --- |
| V0 | E_V invariance under Qb(b(r)): symbolic, exact |
| V1 | the symbolic density vs independent NUMERIC differentiation of the dressed field at random off-axis points (numpy, no sympy in the check path) |
| V2 | the constant-b limit vs the lattice instrument (`e_dressed` on the analytic family, the [M5.21.8](../findings/m5_21_8_note.md) pipeline pattern) |
| V3 | the m\* self-gate (S1a): structural agreement with artanh(1/g) across g |
| Audit | independent adversarial agent, own derivation route, tries to REFUTE the leading-term form and the sign verdict before anything is trusted (cardinal rule) |

**Blindspot pass** (unfamiliar-territory task, run at PLAN): (1) the double-expansion ordering ambiguity (b vs 1/g): routed machine-checkable, both regimes carried; (2) r → 0 regularity: r̂ is singular at the origin, so b(0) = 0 is a required boundary condition and the functional's core behavior is stated explicitly (the author's own "boost radius dependence" caveat); (3) the base-texture question (vacuum / hedgehog / census endpoint): vacuum + hedgehog run symbolically, the census endpoint stays numeric-only (S3); (4) the clock must be IN scope (the criterion reads on oscillations): carried at O(ω²) exactly as the author's own Hm; (5) sympy tractability: staged reduction (K-algebra, radial/angular split, series-in-b before simplify), try cap 3 strategies per stage, on cap checkpoint + descope with the deviation logged; (6) the s = ±1 sign knob: both carried, never averaged.

**Definition of done**: the first nontrivial (1/g) term derived + numerically verified (V0-V3 green); the sign verdict at the variational minimizer read against the author's pre-stated criterion (both s signs); the 3×3 correction functional written with both retrodictions; adversarial audit recorded; method-note-grade findings ([`../findings/m5_21_14_note.md`](../findings/m5_21_14_note.md), equations first + equation-to-code map + embedded figures + not-computed list).

**Artifacts** (all `m5_21_14_` named): `scripts/m5_21_14_a_symbolic.py` (S1 + V0/V3), `scripts/m5_21_14_b_verify.py` (V1/V2), `scripts/m5_21_14_c_minimize.py` (S2), `scripts/m5_21_14_e_audit.py` (the audit), `data/m5_21_14_*.json`, `plots/m5_21_14_panel.png`, `findings/m5_21_14_note.md`, checkpoint `checkpoints/m5_21_14_progress.md`.

**Model/effort**: Fable / high (symbolic-derivation-heavy, compute-light; no lattice relaxations).

### What this task does NOT do

| Non-goal | Why |
| --- | --- |
| Reopen route (b) | the [framework](../findings/m5_21_11_framework.md) is frozen and its terminal verdict stands; this task builds a NEW instrument generation, it does not re-analyze the ladder data |
| Quote any mass ratio | B/C remain uncertifiable at N = 48; a symbolic term changes nothing about the certification problem |
| Wait on the potential details | the author states the current eigenvalue potential is "just a first guess" ([Q25](../m5_question_tracker.md#q25-detail)); S1 runs on the verified quartic L + the T2 base of record, and the V-dependence of the term's sign is REPORTED, not assumed away |

**Gated by**: none (picked by the user 2026-08-09 over [M5.22.3](m5_22_3_task_details.md) /
[M5.22.5](m5_22_5_task_details.md)).

## DEVIATIONS LOG

| When | Deviation |
| --- | --- |
| 2026-08-09 S2 first pass | The planned direct minimization of T1_static over β(r) is ILL-POSED: the term alone is unbounded below (a pure-radial UV channel: grid-scale β oscillation kills the quartic \|W\|² identically while the negative −2\|G_iv_j−G_jv_i\|² channel grows as β′²; measured runaway to −1.04e7 with a ±11 sawtooth). Conservative option taken: the runaway is recorded as a FINDING (it constrains how the term can enter a 3×3 ladder), and the sign-verdict minimization moved to the EXACT dressing functional at finite g (the author's own "minimization of general b(r)" object) within a smooth radial family, with an explicit sawtooth boundedness probe of the exact functional |

## FINDINGS

Run 2026-08-09 (go 12:29 EDT). Full method-note record: [`../findings/m5_21_14_note.md`](../findings/m5_21_14_note.md) (equations, code map, gates, verdicts). The author's three recipe steps all delivered; every derivation claim machine-gated; adversarial audit 7 CONFIRMED / 2 PARTIAL / 0 REFUTED with every fix adopted (note § 7).

**The term (the author's step 2, derived + twice-verified + numerically confirmed at O(1/g))**:

```text
T1_static[beta; M3] = 4 SUM_{i<j} [ 2<[G_i,G_j], v_j v_i^T - v_i v_j^T>
                      + |v_j v_i^T - v_i v_j^T|^2 - 2|G_i v_j - G_j v_i|^2 ]
T1_kin[beta; M3]    = -8 SUM_i |Mdot3 v_i|^2        (omega^2 multiplier)
G_i = d_iM3,  v_i = d_i(beta(r) nhat),  beta = g*b  (the only nontrivial ordering:
at fixed b the density grows as g^4),  E_V exactly invariant under any b(r)
```

| Finding | Content |
| --- | --- |
| 1. The author's criterion is MET at the term level | T1_kin is NEGATIVE-SEMIDEFINITE: any dressing with Ṁ3v_i ≠ 0 for some i lowers the ω² coefficient, for any texture, both (−g)^p signs; measured −401 at the g = 32 minimizer, g-flat, instrument-robust (1.6-6% lattice agreement). Favored, not proven: the flipped-sign regime needs a stabilizer (fixed-J) before ω\* is finite |
| 2. The kin flip is BULK and sits at threshold | per unit radius (r ∈ 10-20): base +13.5 vs dressing −19.6, net −6.1: the dressed constant-ω ledger descends with box radius; the certified L = 48 lattice measures kin_corr −426.3 vs base +426.5 (threshold). Favorability, not proven oscillation: the ω-stabilizer is the standing [Q35](../m5_question_tracker.md#q35-detail) item |
| 3. The flat-gain record is DERIVED | the M5.21.11 F4 flat gain (q ≈ 0) = the O(1) depth of the collapsed β = g·m functional; V2 reproduces the M5.21.8 m\* ratios to 3 digits (0.8168/0.8275/0.8329) and measures the gain flat to 1.3% |
| 4. The m\*-critique is MEASURED | variational b(r) beats the rigid constant ×2.6 ON the certified lattice (−160 vs −61 at n = 32) and ×10 on the continuum quadrature (−4333 vs −424); the optimal profile is core-weighted and oscillatory |
| 5. T1 alone is UNBOUNDED BELOW | the pure-radial UV channel (β′ large, β small) kills the quartic \|W\|² identically while −2\|G_iv_j−G_jv_i\|² grows as β′² (runaway to −1.04e7 measured); the exact functional restores a per-wavelength floor that deepens ∝ λ^−1.8: regulator-set depth ([Q25](../m5_question_tracker.md#q25-detail) connection) |
| 6. The 3×3 handoff (the author's step 3) | E[M3, β] = E_3×3[M3] + T1_static[β; M3], kin correspondingly, with β(r) a NEW g-free O(1) radial field: the corrected ladder stays g-free at leading order. MANDATORY guard from finding 5: T1 must never be minimized freely (constrained profiles, finite-g stabilizers, or the future regularization potential) |

![panel](../plots/m5_21_14_panel.png)

**Instrument agreement**: masked (r > 3) static deviation 37% at h = 1.5 → 9.9% at h = 0.75 (the h² factor); kin 1.6-6% everywhere: signs and orderings are instrument-robust; continuum magnitudes are quadrature statements about a regulator-set channel. **Not run**: the corrected-3×3 ladder itself (needs the finding-6 guard + B/C certification, [M5.21.11](m5_21_11_task_details.md) finding 2 stands); the census-texture minimization; time-averaged clocks; the s-odd next order.

Artifacts: `scripts/m5_21_14_{a_symbolic, b_verify, c_minimize, d_resolution, e_audit}.py` · `data/m5_21_14_{symbolic, verify, minimize, resolution, audit}.json` · [`../plots/m5_21_14_panel.png`](../plots/m5_21_14_panel.png) · [`../findings/m5_21_14_note.md`](../findings/m5_21_14_note.md). All data JSON-small and tracked (no heavy arrays; every number regenerable by the named script, runtimes 6 s to 4 min).

## TASK REVIEW (2026-08-09)

Task Duration: 01:28 (from 12:29 to 13:57 EDT)
Usage Cap Triggered: NO

| # | Result | Detail |
| --- | --- | --- |
| 1 | ✅ The term DERIVED | T1_static + T1_kin as exact identities, β = g·b the only nontrivial ordering (fixed-b density ~ g⁴), E_V exactly invariant; gates exact / 1e-48 / O(1/g)-clean |
| 2 | ✅ The criterion MET at the term level | T1_kin negative-semidefinite (strict when Ṁ3v_i ≠ 0); −401 at the g = 32 minimizer, g-flat; favored-not-proven, fixed-J stabilizer still required |
| 3 | ✅ Retrodictions | the M5.21.11 flat gain DERIVED (O(1) depth of the collapsed functional, flat to 1.3%); the M5.21.8 m\* ratios reproduced to ~0.002; the m\*-critique measured ×2.6 on-lattice |
| 4 | 🔶 The kin flip is BULK, at threshold | net slope −6.1 per unit R (base +13.5, dressing −19.6); the certified L = 48 box at threshold (−426.3 vs +426.5) |
| 5 | ⚠️ T1 alone UNBOUNDED below | pure-radial UV channel; exact-functional floor ∝ λ^−1.88: regulator-set; the 3×3 handoff carries a MANDATORY guard |
| 6 | ✅ Audit | 7 CONFIRMED / 2 PARTIAL / 0 REFUTED, independent routes; every fix adopted + the evenness sharpening (exact at finite g) |

Issues: none blocking; honest caveats carried in the note (favorability ≠ onset; regulator-set continuum magnitudes; two-point h-trend, masked-region-only attribution).

Deviations from plan: one, logged as it happened (the T1-alone minimization ill-posed → recorded as finding 5, sign-verdict moved to the exact functional).

Action taken at review (user-approved): roadmap row → Done; change-log entry; Q33 + Q35 + Q25 receipts; author close-out message drafted (terminal-only, user sends).

Findings
The author's recipe delivered end-to-end: the first nontrivial (1/g) term of the boost-hedgehog dressing is derived, triple-verified, and MEETS the author's pre-stated criterion (negative Hamiltonian contribution in the oscillation sector at the term level), with the constant-ω kin flip measured as a bulk effect sitting exactly at threshold on the certified box. The term also retrodicts the M5.21.11 flat-gain terminal finding and the author's m\*-too-large critique, and its standalone unboundedness (a scale-free negative channel whose floor is regulator-set) is the concrete reason any corrected 3×3 ladder must carry the dressing guarded, tying directly to the author's own open regularization item.

Research docs created/updated
- [m5_21_14_task_details.md](m5_21_14_task_details.md) (this record)
- [../findings/m5_21_14_note.md](../findings/m5_21_14_note.md) (method note: equations, code map, gates, verdicts, audit)
- [../plots/m5_21_14_panel.png](../plots/m5_21_14_panel.png) + 5 scripts + 5 data JSONs (see FINDINGS artifacts)
- [../m5_roadmap.md](../m5_roadmap.md) (Done row, STATUS, change-log) · [../m5_question_tracker.md](../m5_question_tracker.md) (Q33/Q35/Q25 receipts)
