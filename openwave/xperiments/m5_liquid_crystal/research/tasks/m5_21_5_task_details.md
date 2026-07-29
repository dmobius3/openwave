# M5.21.5: the μ + g-factor closure under the verified L (the 4-observable electron)

**Status**: ✅ CLOSED COMPLETE 2026-07-21 (run 10:55-12:00 EDT; terminal review approved same day; was: 🚧 PLANNED STUB 2026-07-16, created as M5.21.3 at the electron-hunt scorecard sync: [`../m5_particle_hunt.md § ELECTRON HUNT`](../m5_particle_hunt.md) coverage map; the g ≈ 2 row was the one non-validated row no task covered. **Renumbered M5.21.5 on 2026-07-17**, run-order renumber, user decision: it runs LAST in the physics series, after the 3D scan [M5.21.2](m5_21_2_task_details.md), the 4D sequel [M5.21.3](m5_21_3_task_details.md), and [M5.21.9](m5_21_9_task_details.md)/[M5.21.10](m5_21_10_task_details.md)/[M5.21.12](m5_21_12_task_details.md); the 2-particle films M5.21.4 follow it, then the RENDERING UNLOCK). Best-state selection resolved at PLAN: the [M5.21.9](m5_21_9_task_details.md) fixed-J conjugation-tangent states (preferred per the roadmap gate) + the T2 sym census family for the ladder.

**Lineage**: [M5.21.1](m5_21_1_task_details.md) (supplies the converged verified-L state this task measures on) · the canonical-era baseline: g = 1.97 at 24³ via the K = 4/α bridge, box ladder [1.97, 2.22] bracketing 2.0023, #219 ([`../findings/m5_8_2za_findings.md`](../findings/m5_8_2za_findings.md), [`../scripts/m5_8_2za_g_factor.py`](../scripts/m5_8_2za_g_factor.py)) · the μ channel mechanism: clock tilt/precession, pure twist EM-silent, EID-C ([`../scripts/m5_8_2r_electron_id.py`](../scripts/m5_8_2r_electron_id.py)) · the e_scale anchor: c₂ = αħc/64π locked at M5.16 ([`../findings/m5_16_report.md`](../findings/m5_16_report.md)).

## Scope (stub level)

| Step | Content | Kill / survive (pre-registered at PLAN) |
| --- | --- | --- |
| G1 | Re-establish the μ channel (clock tilt/precession) on the M5.21.1 converged 4D state: verified L, not the canonical stack | channel present/absent measured either way |
| G2 | The μ box-convergence ladder (the #219 caveat: μ was not box-converged) | monotone ladder or the honest divergence record |
| G3 | Replace the structurally-motivated K = 4/α bridge with the locked Coulomb e_scale (c₂ = αħc/64π, M5.16) for a first-principles g read | the bridge either derives or stays flagged as the residual assumption |
| G4 | Verdict: g vs 2.0023 + the 4-observable electron (mass, charge, μ, J) closed under the verified L | closure, or the measured characterization of the gap |

**Pre-registered contingency**: if M5.21.1 ends WITHOUT a converged verified-L clock state, this task re-plans on the best measured state (partials labeled 🔶) or holds, user's call at that review.

## Records + gating

| Item | Content |
| --- | --- |
| Destinations | scripts `../scripts/m5_21_5_*.py` · data/plots `m5_21_5_*` (film standard, both templates where states evolve) · findings `../findings/m5_21_5_method_note.md` · checkpoint at go |
| Standards | method note (equations first, code map) + independent adversarial audit (cardinal rule) |
| Scorecard | closes the electron-hunt g ≈ 2 row ([`../m5_particle_hunt.md`](../m5_particle_hunt.md)); MODELS.md's μ + spin row updates in the post-program sweep (user 2026-07-16: MODELS.md stays as-is until the particle-hunt program completes) |
| Gated By | the best stable state from the series (the [M5.21.2](m5_21_2_task_details.md) 3D minima / the [M5.21.3](m5_21_3_task_details.md) 4D rotating state; the M5.21.1-negative contingency clause above stands) + user "go" |

## TASK PLANNING (2026-07-21, at go)

**Go-time**: 2026-07-21 10:55 EDT. Resume ping armed (`Resume: Task M5.21.5`, fires 14:55 EDT = reset 14:50 + 5); reset-time watchdog running.

### The measurement, mapped onto the verified-L era

The canonical-era result ([`../findings/m5_8_2za_findings.md`](../findings/m5_8_2za_findings.md)) paired the EM-active TILT μ (0.2209 at 24³) with the twist-clock Noether spin (61.61) and crossed sectors with the structurally-motivated K = 4/α, landing g = 1.97 with the ladder [1.97, 2.22] bracketing 2.0023. Three things have changed since: (a) the state is now the RELAXED T2 eigenvalue-penalty electron, not the analytic conjugated seed; (b) the clock is now the [M5.21.9](m5_21_9_task_details.md) fixed-J construction on the conjugation tangent (rotation about the LOCAL leading eigenvector = the twist channel per-voxel), with J carried exactly and kin measured; (c) the M5.16 Coulomb lock (c₂ = αħc/64π) is available as a first-principles anchor. The key algebraic observation made at PLAN (to verify in G3): with c₂ = αħc/64π the far-field hedgehog energy density 8c₂/r⁴ equals EXACTLY αħc/8πr⁴, the Gaussian-units energy density of the point charge e: the emergent-EM field normalization is therefore DERIVABLE, and in the g assembly the length unit cancels, leaving g = 2·k·μ_lat·E_lat/S_lat with E_lat the anchored rest energy and k a pure convention factor traceable by hand through the definitions (B = ½εF, μ = ½∫r×j, the F² coefficient). Whether k derives cleanly decides G3's kill/survive.

### Steps (pre-registered)

| Step | Content | Kill / survive |
| --- | --- | --- |
| P0 | Instrument gates: (i) my Mermin-Ho pipeline reproduces the closed-form hedgehog curvature density 8c₂/r⁴ shape on the analytic seed; (ii) EID-REPRODUCTION gate: on the analytic 24³ pinned-clock seed, my reimplementation reproduces μ_tilt ≈ 0.221 and L_int(twist) ≈ 61.6 (the m5_8_2r record) before touching any relaxed state | gates fail → fix before any physics claim |
| G1 | The μ channel on the verified-L states: native fixed-J clock read (expect μ ≈ 0, the twist EM-silence, MEASURED not assumed) + the tilt linear-response μ (global internal rotation probe, 3 axes) on the fixed-J rungs (om 0.2/0.5/1) and the t32_A census electron; defect-core masking by the λ₂−λ₃ criterion (relaxed states carry the braided ½-lines; Mermin-Ho is ill-defined on the cores; mask disclosed) | channel present/absent measured either way |
| G2 | The μ box-convergence ladder on the MATCHED T2 sym pinned family: n = 24 (fresh descent), 32 (t32_A on disk), 48 (fresh descent), h = 1.5; plus the radial-window profile μ(<R) at each n (the 2za caveat named tail-domination: measure it, not just the endpoint); f48/f64 free-endpoint reads as a preparation-robustness check if time allows | monotone ladder or the honest divergence record |
| G3 | The first-principles bridge: derive the μ and S physical conversions from the M5.16 lock (c₂ = αħc/64π, length unit, m_e anchor) with every convention factor traced by hand; assemble g with NO free factor; compare vs the legacy K = 4/α on the same raw numbers | the bridge either derives or stays flagged as the residual assumption (pre-registered honest outcome) |
| G4 | Verdict: g vs 2.0023 + the 4-observable table (mass, charge, μ, J) under the verified L | closure, or the measured characterization of the gap |
| P5 | Independent adversarial audit (cardinal rule): own implementations of the director extraction, Mermin-Ho curvature, μ integral, and the bridge algebra; verdicts per claim recorded in the note | audit catches adopted before FINISH |

### Blindspot pass + unknowns routing

| Unknown | Class | Routing |
| --- | --- | --- |
| Director sign field on relaxed states (headless n; random eigenvector sign flips corrupt derivatives) | machine-checkable | align to the radial hemisphere (n·r̂ ≥ 0), measure the flip-boundary artifact by comparing against a continuity-propagated alignment on one state |
| Mermin-Ho on the ½-line cores (Q₈ topology; the abelian projection has 2π ambiguities there) | machine-checkable | mask by the biaxial split criterion (M5.23 Stage D: λ₂−λ₃), report masked volume; μ(<R) profiles show core-vs-tail weight |
| Tilt probe at finite angle destroys the hedgehog (EID-C finding 4) | known structure | linear response only (analytic tangent, no finite rotation) |
| The time-unit / action-sector conversion in G3 (c_lat = 1 assumption enters S_phys) | derivable-with-assumption | state the assumption explicitly; it is THE candidate residual if k does not close |
| The pairing choice (mixed μ_tilt/S_twist vs matched) | ours, EID-established | keep both, mixed = physical per the 2za record; flagged in the note |
| Whether the fixed-J J and the ladder L_int agree as the S in g | machine-checkable | compute both on the same state; discrepancy = the kin-convention flag (M5.23.1 stub ⚠️) made quantitative |
| f48/f64 endpoint key formats | trivial | inspect at load; fall back to the matched pinned family (the ladder never depends on them) |

**Contingency (pre-registered)**: if the fresh 24³/48³ descents do not converge in the window, the ladder reports {32} + radial windows with the descents queued; if the tilt μ is ZERO on relaxed states too (channel absent under the verified L), that IS the G1 result and G3/G4 run on the canonical-era contrast honestly.

**Destinations**: scripts `../scripts/m5_21_5_{a_mu,b_ladder,c_bridge,d_audit}.py` · data/plots `m5_21_5_*` · findings note `../findings/m5_21_5_note.md` (recent-series naming; the stub's `m5_21_5_method_note.md` name superseded, logged) · checkpoint `../checkpoints/m5_21_5_progress.md`.

## EXECUTE log: deviations from plan (logged as they happened)

| When | Deviation | Why + status |
| --- | --- | --- |
| 11:20 | Findings file named `m5_21_5_note.md`, not the stub's `m5_21_5_method_note.md` | recent-series naming consistency (m5_21_9/10 pattern); links updated here |
| 11:40 | The v2 transverse gauge became load-bearing (not in plan) | eigh's per-voxel sign ambiguity decoheres every frame-based tangent: ungauged = 2.5× low, azimuthal reference = 22× high, POLAR reference = the EID record to 8 digits; the polar gauge is now part of the instrument definition |
| 12:05 | Instrument v3: envelope-localized clock flow replaced the rigid flow as PRIMARY | measured: the rigid moment on relaxed pinned states is pin-layer/corner dominated (fjom: half the moment outside r = 16; t24A cancels 60 → 3); the envelope flow (the m5_21_9 state-of-record clock shape) converges radially by R12-14; rigid kept as diagnostic |
| 12:20 | E of record per state + S_env convention pinned mid-run | the 4D-embedded energy of a 3×3 state (13.45) is not its T2 energy of record (4.77); the rigid S is IR-extensive (the M5.21.8 constant-ω pathology): both conventions were pre-registered as machine-checkable unknowns and resolved by measurement |
| 12:40 | seed32 comparator added (not in plan) | the suppression claim needed the grid-matched unrelaxed reference; it landed the texture-tracking headline (seed-A tilt μ = 6.6e-4 vs EID-seed 0.221 vs relaxed 0.018-23) |

## Datasets (heavy arrays LOCAL + gitignored; JSON + plots tracked)

| Artifact | Regen | Runtime |
| --- | --- | --- |
| `data/m5_21_5_end_t24_A.npz` + row | `python3 scripts/m5_21_5_b_ladder.py n=24` | 45 s |
| `data/m5_21_5_end_t48_A.npz` + row | `python3 scripts/m5_21_5_b_ladder.py n=48` | ~80 min |
| `data/m5_21_5_seed32.npz` | `INS.make_seed(base_cfg(term="T2", n=32, L=48, bc="pinned", seed="A"))` (see the seed32 block in the note § 3) | seconds |
| `data/m5_21_5_gates.json` | `python3 scripts/m5_21_5_a_mu.py gates` | ~2 min |
| `data/m5_21_5_mu_<tag>.json` (seed32, t24A, t32A, t48A, fjom0.2/0.5/1) | `python3 scripts/m5_21_5_a_mu.py read tag=<t> path=<npz>` | ~2 min each |
| `data/m5_21_5_bridge.json` | `python3 scripts/m5_21_5_c_bridge.py` | seconds |
| `data/m5_21_5_audit.json` | `python3 scripts/m5_21_5_d_audit.py` | audit-owned |

## RESULTS (2026-07-21; canonical record = [`../findings/m5_21_5_note.md`](../findings/m5_21_5_note.md))

| Step | Result |
| --- | --- |
| P0 gates | ✅ closed-form hedgehog 2% at h = 1.5; EID reproduction EXACT (0.22090514 vs 0.22090513, 8 digits); twist EM-silence 2.4e-8; the polar transverse gauge is load-bearing (part of the instrument definition) |
| G1 μ channel | 🔶 EXISTS and radially converged per state (envelope-localized clock, plateau R12-14), but PREPARATION-FRAGILE: μ_clock(tilt12) spans seed-A 6.6e-4 / census 0.018-23.1 / fixed-J 3.38; the audit's mechanism: a PARITY-CANCELLATION residue (net/gross 5.3e-3 vs 5.1e-5; gross responses only 1.8× apart); gauge-free globalz uniformly small (2e-5 to 0.7); twist stays EM-silent on states; a STATIC ∇×B moment exists (0.118 fixed-J, stable across rungs) |
| G2 ladder | ❌ divergence record, sharper than #219: NON-MONOTONE 23.1 → 0.018 → 16.4 across n = 24/32/48 (same recipe → different basins; family E non-monotone 4.77 → 4.91; r_half tracks the box 9.3 → 17.7); μ does not box-converge on descent endpoints |
| G3 bridge | ✅ DERIVED with no free factor: g = μ_lat·E_lat/(2π·S_lat) (the Coulomb anchor makes the hedgehog charge exactly e; the length unit cancels; audit C1 exact to 4e-16). The legacy K = 4/α remains underived; k_needed/k_derived = 20.5 on the state of record |
| G4 verdict | ❌ NO closure: g_fp spans 8.5e-4 (t32A) to 1.45 (t24A), 0.097 on the fixed-J state of record; the canonical 1.97 was a box-truncated μ read × the underived bridge (audit C4 retro-refuted the record's convergence: μ ∝ R^3.22, 27.8% from corners). The 4-observable electron stands at mass ✅ charge ✅ J ✅ μ 🔶 g ❌ with the closure routes named (dynamically-selected texture via M5.23.1; measured c_lat; h-refinement) |
| P5 audit | ✅ 6 claims: C1/C2/C5/C6 CONFIRMED, C3 confirmed with the mechanism upgrade, C4 values-confirmed/convergence-REFUTED; 2 catches adopted (k_needed 2× bug; IR growth linear-in-box); new structure: the fjom moment is equatorial, 4.9° from the mean nematic axis, 61% from shell r = 6-8, rung-stable to 0.5° |

**Scorecard consequence**: the electron-hunt g ≈ 2 row does NOT close; it re-bases from "right-order ✅ (1.97-2.22 brackets 2.0023)" to the honest characterization above. MODELS.md stays untouched until the program-complete sweep (user 2026-07-16).

## TASK REVIEW (2026-07-21)

Task Duration: 01:05 (from 10:55 to 12:00 EDT)
Usage Cap Triggered: NO (resume ping parked unfired; watchdog stopped)

| Step | Result |
| --- | --- |
| P0 instrument gates | ✅ EID record reproduced EXACTLY (0.22090514 vs 0.22090513, 8 digits, same seed); closed-form hedgehog 2% at h = 1.5; twist EM-silence 2.4e-8; the POLAR transverse gauge is load-bearing (ungauged 2.5× low, azimuthal 22× high) |
| G1 the μ channel | 🔶 exists + radially converged within any fixed state (envelope clock, plateau R12-14) but NOT an electron invariant: 6.6e-4 (seed) → 0.018 (census 32³) → 3.38 (fixed-J) → 23.1 (squeezed 24³); mechanism (audit): a PARITY-CANCELLATION residue (gross responses only 1.8× apart; net/gross 5.3e-3 vs 5.1e-5) |
| G2 the ladder | ❌ divergence record: non-monotone 23.1 → 0.018 → 16.4 (n = 24/32/48, same recipe → different basins; E 4.77 → 4.91 non-monotone; r_half 9.3 → 17.7 tracks the box) |
| G3 the bridge | ✅ DERIVED, no free factor: g = μ_lat·E_lat/(2π·S_lat); q = e exact by the Coulomb anchor; length unit cancels; audit-exact 4e-16; K = 4/α stays underived (k_needed/k_derived = 20.5 on the state of record) |
| G4 verdict | ❌ NO closure: g spans 8.5e-4 to 1.45 (0.097 on the state of record); the canonical 1.97 retro-flagged (the 0.221 grows as R^3.22 with cutoff: box-truncated observable through an underived bridge) |
| P5 audit | ✅ 6 claims independent: 4 confirmed, C3 mechanism-upgraded, C4 convergence-REFUTED; 2 catches adopted (k_needed 2× bug; linear-in-box IR growth); new structure: equatorial moment, 4.9° from the mean nematic axis, rung-stable ≤ 0.5° |

Issues: t24A box-squeezed (labeled); census h = 1.5 resolution unseparated at the restructured core; c_lat = 1 assumption in the S conversion (→ M5.21.11).

Action taken at close (approved): review written here; roadmap row → Done (appended at the end); high-level sweep (particle hunt g row re-based, canonical recipe + anti-recipe, briefing live-arc, M5.23.1 kin-convention quantification).

Findings: Under the verified L with a first-principles bridge (derived from the Coulomb anchor, no free factor, audit-exact), the electron's g-factor does NOT close at 2: the μ observable is a parity-cancellation residue tracking preparation basin and transverse texture across 4 orders of magnitude, and the canonical-era g = 1.97 landing is retro-flagged as a box-truncated read through an underived bridge. The 4-observable electron stands at mass ✅ charge ✅ (exactly e by anchoring) J ✅ μ 🔶 g ❌; the strongest closure route is reading μ on the dynamically-selected texture of a long fixed-J evolution ([M5.23.1](m5_23_1_task_details.md)).

Research docs created/updated: [`m5_21_5_task_details.md`](m5_21_5_task_details.md) (this record) · [`../findings/m5_21_5_note.md`](../findings/m5_21_5_note.md) (method note: equations, code map, gates, divergence record, bridge derivation, audit, panel) · scripts `m5_21_5_{a_mu,b_ladder,c_bridge,d_audit,e_panel}.py` · data `m5_21_5_{gates,bridge,audit}.json` + 7 μ-read JSONs + 2 descents (npz local, manifest regenerated) · [`../plots/m5_21_5_panel.png`](../plots/m5_21_5_panel.png) · sweep: [`../m5_particle_hunt.md`](../m5_particle_hunt.md), [`../m5_theory_canonical.md`](../m5_theory_canonical.md), [`../../__M5_model_briefing.md`](../../__M5_model_briefing.md), [`m5_23_1_task_details.md`](m5_23_1_task_details.md), [`../m5_roadmap.md`](../m5_roadmap.md).
