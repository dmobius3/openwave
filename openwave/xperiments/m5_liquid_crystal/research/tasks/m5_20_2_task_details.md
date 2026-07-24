# M5.20.2 - Full 4×4 (g, 1, δ, 0) dynamics: the clock sector on the loop

**Status**: ✅ DONE 2026-07-12 (review approved; roadmap DONE row; the combined email awaits the maintainer's send). **Spec**: Duda 2026-07-11/12 ([`m5_20_convo.md`](m5_20_convo.md)): "finally to get oscillations requires full 4x4 tensor with potential preferring (1,g,delta,0) eigenvalues" (his 07-12 ordering read as a misspelling of the standing `(g, 1, δ, 0)`, user call 2026-07-12) + "clock propulsion with negative Hamiltonian terms require full 4x4 tensor field with (g, 1, delta, 0) spectrum" (07-11). Predecessor: [M5.20.1](m5_20_1_task_details.md) (the spatial (1, δ, 0) verdict: UNWOUND across the grid; the spectrum objection fails in-framework; the time sector is protection's last in-framework door). Delivery intent (user call 2026-07-12): M5.20.1 + M5.20.2 results go to Duda in ONE combined email, no intermediate ask round, PROVIDED no load-bearing fork requires his answer first; the § Q23-resolution plan below is how we avoid the turnaround.

## TASK PLANNING

### Scope

Take the M5.18-VERIFIED 4D Lagrangian (his own, machine-checked: Lorentz invariance to 1.3e-11, Legendre transform exact, `m5_18_lorentz_check.py`) and RUN it: derive the M-variable equations of motion in the axisym reduction, gate them, and evolve the 4×4 field with the potential pinning `(g, 1, δ, 0)`, dynamical time block included, on (a) the uniform vacuum (clock baseline) and (b) the M5.20.1 loop/remnant configurations (the clock on the loop). The physics question: does the time sector carry the oscillation his neutrino picture requires, and does it change the protection verdict?

### The Q23 resolution path (why no author turnaround is needed)

Q23 asked for his time term in M variables. We DERIVE it instead of asking: the M5.18 file already holds L in field variables,

```text
L = − Σ_{μ<ν} F_{μν} : F^{μν} − V,   F_{μν} = [∂_μ M, ∂_ν M]_η,   η = diag(−1,1,1,1)
V = w Σ_{p=1..4} (Tr_η(M^p) − C_p)²,  C_p = (−g)^p + 1 + δ^p + 0^p  (mixed-index trace; sign per branch)
```

so the kinetic sector is the F_{0i} block: velocity enters ONLY through `[∂_t M, ∂_i M]_η` (no standalone Mdot² term). The EOM derivation is mechanical and MACHINE-CHECKABLE (FD-gate the derived force against the action, the standing S1/GF pattern). His confirmation is then requested in the results email (evidence, not blocking).

### Known landmines (M5.18 findings 3/4b/4c, the honest hazard list)

| Landmine | M5.18 evidence | Mitigation in this plan |
| --- | --- | --- |
| Degenerate Legendre map ([Q18](../m5_question_tracker.md#q18-detail)) | the Mdot ~ η direction never enters L (primary constraint, check 3) | constrained (M, P) evolution, the M5.8.2c Option B precedent; the constraint set derived at phase A, gated |
| H unbounded below on the vacuum manifold | boost×rotation product textures give strictly negative density (check 4b, audit-corrected) | phase-B linearized mode scan around the branch vacuum BEFORE any production run; if the loop/remnant configurations couple to the runaway sector, STOP and report (that itself answers "does the clock sector protect": it destabilizes instead); energy-monitor abort guards on every run |
| Branch choice ([Q19](../m5_question_tracker.md#q19-detail)) | V = 0 is a UNION of 4 Lorentz orbits: which eigenvalue rides timelike (check 4c) | primary runs on the g-timelike branch (his stated convention); ONE alternate-branch control; possible domain walls noted as out of scope |
| Static-sector consistency | check 6: uniform-time-row fields reduce EXACTLY to the M5.17 static stack | gate GB0: the M5.20.1 spatial results must be reproduced bit-level when the time block is frozen |

### Definition of done

1. **EOM derivation + gates**: the axisym 4×4 force (curvature-η + spectral-4 targets) FD-gated against the action (< 1e-8 directional); GB0 static-sector reduction (frozen time block == the M5.20.1 stack, ≤ 1e-12); the constrained integrator re-gated on the spectral potential (energy ledger closure, dt² scaling).
2. **4×4 gap map**: the 10×10 V-Hessian at diag(g, 1, δ, 0) per branch: zero modes, ω ladder, the negative-direction census (which modes are runaway vs oscillatory at the linear level): the clock-frequency prediction ω_clock(g, δ).
3. **Clock baseline**: small dynamical time-block excitation on the uniform branch vacuum: measure the oscillation spectrum vs the gap map (the Zitterbewegung analog of the M5.8.2 line, now on his verified Lagrangian).
4. **The clock on the loop**: time-block excitation on (a) the M5.20.1 biax loop seed and (b) the pair_1d unwound remnant endpoint: does a coherent clock line appear, does it shift the unwinding, any sign of the "clock propulsion" his 07-11 reply names.
5. **Protection re-probe**: one (δ = 0.3, pair_1d) loop run with the FULL 4×4 dynamics: does the time sector change the UNWOUND verdict (the M5.20.1 classifier + instruments re-used verbatim).
6. Method-note-grade close page (FIELD CONTENT box: now the time block is DYNAMICAL) + adversarial audit + tracker updates (Q18/Q19/Q23) + the combined-email deliverable with M5.20.1.
7. Seed + endpoint cross-section maps rendered for the 4×4 runs (the M5.20.1 seed-map/endpoint-map pattern, per the user's 2026-07-12 request).

### Blindspot pass

| # | Blindspot | Mitigation |
| --- | --- | --- |
| 1 | The η-commutator curvature makes the FAST gradient wrong if naively reused (η insertions change the adjoint scatter) | derive the adjoints fresh; gate against FD; keep the slow einsum path as the audited reference (the GF pattern) |
| 2 | Runaway may be SLOW (negative density is texture-specific): a run can look stable for hundreds of t then dive | energy-floor abort guard + the phase-B mode census sets the expected runaway rate; budget runs accordingly |
| 3 | The uniform time row of every M5.20.x seed means t-block excitation must be INJECTED; the injection profile could itself be a boost texture (the 4b sector) | rotation-sector injection only (time-mixing generator with spatial rotation, not boost); the injection classified against the 4b census before use |
| 4 | wscale for the 4-target potential is a new calibration question | reuse w = 7.24e-4 for comparability + one autochi control (the M5.20.1 recal pattern) |
| 5 | Axisym reduction may not accommodate the time-mixing channels cleanly (extra equivariant background terms) | the phase-A reduction derivation enumerates ALL channels; the A_φ lesson from M5.20.1 applied to the (0, i) block |
| 6 | Combined-email scope creep: M5.20.2 could balloon and delay the M5.20.1 delivery | HARD CAP: if phase B declares the sector ill-posed without author input, or production cannot land within ~1 day, ship M5.20.1 + the derivation + the ONE sharp question instead (the user's stated fallback) |

### Phased plan

| Phase | Content | Budget |
| --- | --- | --- |
| A: derivation + gates | axisym 4×4 EOM (η-curvature adjoints + 4-target spectral force), FD gates, GB0 static reduction, constrained-integrator re-gate; the 10×10 gap map + negative-mode census per branch | afternoon block (derivation-heavy, CPU-light) |
| B: well-posedness triage | linearized runs on the branch vacuum: mode census confirmed dynamically; the injection protocol classified; GO/NO-GO on production | ~1h |
| C: clock baseline + clock-on-loop | DoD 3-4 runs (short horizons first, T ~ 200-500) | evening |
| D: protection re-probe + production | DoD 5 + the oscillation measurement matrix (g fixed 8.0; δ ∈ {0.3} primary, {0.1, 0.5} if budget allows; one alternate-branch control) | overnight |
| E: synthesis | verdicts, maps, method note, audit, the combined email draft (FABLE VOICE, 4×4 usage reinforced) | tomorrow morning |

### Preconditions

| Precondition | State |
| --- | --- |
| His 4D Lagrangian verified (L, H, constraint, witnesses) | ✅ [`m5_18_lorentz_check.py`](../scripts/m5_18_lorentz_check.py) + `../findings/m5_18_verification_note.md` |
| 4×4 substrate + spectral pinning | ✅ engine since M5.8.1; `Σ_p (Tr M^p − c_p)²` pins exactly (M5.18 S3) |
| Constrained 4D integrator precedent | ✅ M5.8.2c Option B (needs re-gating here) |
| The M5.20.1 instrument stack (seeds, winding, core diagnostic, classifier v2) | ✅ this task's predecessors, all gates green |
| The M5.20.1 verdict (the object the clock dresses) | ✅ UNWOUND grid: the remnant + seed configurations are the test articles |
| Blockers | none author-gated: Q23 resolved by derivation (confirmation requested with results); Q18/Q19 handled per the landmine table |

### Research-body destinations

Scripts `m5_20_2_*` in [`../scripts/`](../scripts/), data/plots `m5_20_2_*`, close page `../findings/m5_20_2_method_note.md`, checkpoints `../checkpoints/m5_20_2_progress.md`; roadmap row [`../m5_roadmap.md`](../m5_roadmap.md) (Backlog → In Progress at go); tracker [Q18](../m5_question_tracker.md#q18-detail) / [Q19](../m5_question_tracker.md#q19-detail) / [Q23](../m5_question_tracker.md#q23-detail); convo [`m5_20_convo.md`](m5_20_convo.md).

### Model + effort

Same as M5.20.1 (main-loop model; high effort on the derivation, the mode census, and the audit; standard on run babysitting). Headless NumPy; the 4×4 axisym stack stays at 128×256.

## FINDINGS (2026-07-12, all phases run; adversarial audit DONE: all six claims CONFIRMED, zero refutations; sharpenings folded into the combined note)

**The clock sector is real, quantified, and ill-posed without one structural input from the author.** Everything below was derived and run from HIS OWN machine-verified Lagrangian (M5.18), no new physics assumptions; the combined method note carries the full detail and is the single Duda-facing document ([`../findings/m5_20_2_method_note.md`](../findings/m5_20_2_method_note.md)).

| Result | Status |
| --- | --- |
| Q23 resolved by derivation: his L is PURELY QUARTIC in derivatives (kinetic sector = the F₀ᵢ block); the kinetic form K(M) vanishes on gradient-free states, has the η null direction everywhere (his verified primary constraint), and is rank 8/10 with 3 negative directions at the loop core | ✅ measured (phase A gates GB0-GB4 all PASS) |
| 4×4 gap map at diag(−g, 1, δ, 0): 6 exact boost/rotation flats + ω = {0.0093, 0.0466, 0.1349, 78.28} at δ = 0.3; all four timelike-assignment branches exact V = 0 (Q19 structure) | ✅ measured (analytic, FD cross-check 1.7e-11) |
| Boost-sector RUNAWAY under the canonical completion: time-mixing injections on the biaxial vacuum diverge at t ≈ 21 (audit: dt-robust across 4×, physical; requires spatial gradients: inert on gradient-free backgrounds); rotation textures exactly stable: the M5.18 4b indefiniteness is ACTIVE on physical states: field-level time integration is obstructed | ✅ measured (pre-registered GO/NO-GO → NO-GO; the pivot executed) |
| The ansatz-level clock (exact, no integration): H(ω) = H(0) + ω²K_eff on rigid Lorentz orbits; ALL rotation orbits cost energy (+8e2..+4e4), ALL boost orbits are negative (−3e5..−2.5e6, 100% negative density): "clock propulsion with negative Hamiltonian terms" quantified; H unbounded along boosts ⟹ no stable finite-ω clock without a constraint | ✅ measured (vacuum J12/B03 zeros exact; alternate branch ~100× stiffer) |
| Protection re-probe: UNWOUND under his FULL 4-target potential (p = 1..4, +1.70 seed energy), time row frozen, T = 2000, drift 2.3e-6: the p = 4 target does not rescue the spatial verdict | ✅ measured |
| The one author-gated question: what constraint closes the boost sector (Dirac treatment? rotation-orbit restriction? another branch/sign?) so the 4×4 oscillation run is well-posed; given it, the run is ready same-day | 🚧 outbound with the combined email |

![the ansatz-level clock](../plots/m5_20_2_clock_homega.png)

![triage](../plots/m5_20_2_clock_baseline.png)

![4-target protection re-probe](../plots/m5_20_2_protection_maps.png)

**Deviations from plan**: (1) the pre-registered NO-GO branch fired at phase B (boost runaway), so the original production oscillation matrix was replaced by the ansatz-level H(ω) measurement (blindspot 6's hard cap, logged as it happened); (2) DoD 3-4's "clock baseline + clock on the loop" were delivered through the exact K_eff census instead of time integration (stronger: exact, not integrated); (3) dt production = 0.02 for the frozen-time run (the stiff mode lives in the frozen sector), 0.01 was the triage pick for full-4×4 dynamics.

Checkpoint trail: [`../checkpoints/m5_20_2_progress.md`](../checkpoints/m5_20_2_progress.md). Data: `m5_20_2_a_eom.json`, `m5_20_2_b_triage.json`, `m5_20_2_c_clock.json`, `m5_20_2_run_protection.json`.

## TASK REVIEW (2026-07-12)

**Task Duration:** 04:10 (from 11:03 to 15:13 EDT)
**Usage Cap Triggered:** NO

Approved by the maintainer 2026-07-12 (terminal review). Results: Q23 resolved by derivation (purely-quartic L, degenerate K) ✅; 4×4 gap map analytic ✅; boost runaway physical (dt-robust, gradient-requiring) ✅; the exact H(ω) clock census (rotations positive, boosts negative, unbounded ⟹ the one constraint question) ✅; protection re-probe UNWOUND under the full 4-target V ✅; adversarial audit 6/6 CONFIRMED ✅. Deviations: the pre-registered NO-GO fired at phase B, production replaced by the exact H(ω) measurement (logged). Delivery: ONE combined method note ([`../findings/m5_20_2_method_note.md`](../findings/m5_20_2_method_note.md)) + the FABLE VOICE block handed to the maintainer (he sends after commit/merge); next = await Duda's constraint answer, which unblocks the oscillation run same-day.
