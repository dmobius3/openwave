# M5.21.2: the 3D 3-lepton axis-permutation scan (his 2026-07-17 prescription)

**Status**: ✅ CLOSED COMPLETE 2026-07-17 (go 09:37, duration 1:55, no cap; review approved same day; audited 8/8). Opened 2026-07-17 at the Q25-Q28 reply decode ([`m5_21_convo.md § 2026-07-17`](m5_21_convo.md); staged as M5.21.4 at creation, renumbered M5.21.2 same day in the run-order renumber). Roadmap row: [`m5_roadmap.md § DONE`](../m5_roadmap.md). Successor: [M5.21.2b](m5_21_2b_task_details.md) (the well-posed instrument).

## TASK PLANNING

**Scope**: run the author's 2026-07-17 3D-FIRST prescription end-to-end. Build the genuine 3×3 spatial stack (no time row, no g: vacuum spectrum (1, δ, 0), symmetry group SO(3), full-3D cubic grid), gate it, control the ONE numerical risk he flagged (small box + δ = 0.3 destroying the topological constraints through the boundary), then the scan: energy minimization from the three axis-permutation hedgehog seeds diag(1, δ, 0) / (δ, 0, 1) / (0, 1, δ), censusing the local minima ("global should be electron"; three distinct minima = the lepton-hierarchy signature, FRAME p. 15 = the [M5.9](m5_9_task_details.md) lepton row's first concrete run). The headline physics question: does the 3D hedgehog CONVERGE (topological protection real at statics) where the 4D one diluted ([M5.21.1](m5_21_1_task_details.md) P1 negative, mechanism bracketed at [M5.21.1e](m5_21_1e_task_details.md))?

**Definition of done** (testable):

| # | Criterion |
| --- | --- |
| 1 | Gates green: complex-step gradient check ≤ 1e-10 rel; global SO(3) energy invariance ≤ 1e-10 rel + a negative control; vacuum energy ≈ 0; seed spectrum verified (far-field = the permuted vacuum, two-equal on the vortex axis, three-equal center) |
| 2 | Boundary-integrity ladder measured: hedgehog retention vs (N, δ, BC pinned/free) with a stated working point, or the honest negative that no rung protects |
| 3 | The 3-seed census at the working point: per-seed endpoint characterized (E, virial ratio u/3V, r_half, radial alignment/charge, spectrum profiles) + the three-minima verdict |
| 4 | Independent adversarial audit of every substantive claim before the review |
| 5 | Doc checker clean; TASK REVIEW in the terminal with the two log lines; >1 MB raw data deleted + documented |

**Gating**: user "go" ✅ 2026-07-17. Consumes: the [M5.21.1](m5_21_1_task_details.md) full-3D 4×4 machinery (`m5_21_1_c_4d.py`, ported to 3×3 with η → identity, time row dropped) + the M5.21.1/1e results as comparators. Feeds: [M5.21.3](m5_21_3_task_details.md) (the 3 minima = its 4D descent seeds), [M5.21.5](m5_21_5_task_details.md), the Q25/Q26 discrimination, the Duda thread ([`m5_21_convo.md`](m5_21_convo.md)).

**Blindspot pass** (unknown unknowns surfaced):

| Blindspot | Route |
| --- | --- |
| The 3D potential p-range (p = 1..3 pins all three invariants; p = 1..4 is the 4D habit) | machine-checkable: p = 1..3 is the primary (3 eigenvalues ↔ 3 invariants); a p = 1..4 control at the working point, cheap |
| Boundary-condition semantics: pinned ENFORCES the winding (protection by construction), free TESTS it | run BOTH; the interpretation caution goes in the findings, never blended |
| The M5.21.1 amplitude-dilution escape may reappear in 3D | watch the same diagnostics (spectral force fraction, coreball, r90); pre-registered fallback = a wscale ×100 stiffness arm at the working point |
| Director sign ambiguity (±v identified) breaks naive degree sums | primary observable = the radial-alignment profile a(r) = mean((v·r̂)²) (1 = hedgehog, 1/3 = unwound); the lattice degree is secondary |
| Near-degenerate (δ, 0) pair churns eigenvector ORDER at small δ (the M5.21.1e audit lesson: gap-dip radii unreliable) | census reads use whole profiles, never single radii |
| FIRE dt stability on the 3D h = 1 grid | dt0 = 0.02, force refreshed EVERY iteration, non-finite backoff (the M5.21.1e toy lesson) |
| N = 64 runtime blowup | N = 48 is the workhorse; N = 64 only as the boundary-ladder spot-check |

**Research body**: script `../scripts/m5_21_2_a_scan3d.py` (one script, phase modes: gates / ladder / scan / plots, the M5.21.1e pattern) + the audit agent's own script `../scripts/m5_21_2_audit_check.py`. Data `../data/m5_21_2_*.json` (+ transient `.npz` endpoints, deleted if > 1 MB and documented). Plots `../plots/m5_21_2_*.png` (seed maps + endpoint maps REQUIRED per the field-state rule; ladder heatmap; census panel; film strips adapted to 3×3 meridional cross-sections, both templates where portable). Findings `../findings/m5_21_2_census.md` (method-note-grade: this goes back to the author). Checkpoint `../checkpoints/m5_21_2_progress.md`.

**Sub-experiments** (pre-registered):

| Phase | What | Kill / survive read |
| --- | --- | --- |
| P0 stack + gates | Port the 4×4 3D statics to 3×3 (curvature u = commutator norm over the 3 derivative channels; V = wscale·Σ_p (Tr(M^p) − C_p)², C_p = 1 + δ^p, wscale = 7.2402e-4 carried for continuity); complex-step + SO(3) + vacuum + seed gates | gates green before ANY physics; goal-loop cap 3 tries then stop and surface |
| P1 boundary ladder | Seed A (electron); N ∈ {32, 48, 64} × δ ∈ {0.3, 0.1, 0.03} × BC {pinned, free}; ~800-iter probes; read a(r)/q retention, boundary-layer energy, far-field spectrum drift | any rung where FREE BC keeps the winding → protection measured; only-pinned-holds → protection is boundary-enforced (honest label); pick the working point |
| P2 the 3-seed census | Seeds A/B/C at the working point; long FIRE (≤ 30k iters, convergence = force max + E plateau); snapshots → films; endpoint census (E, u/3V, r_half, a(r), spectra) | three DISTINCT minima vs collapse-to-one vs dilution; the virial balance read per endpoint |
| P3 the charged-ring arm (user amendment 2026-07-17 mid-run, replacing the optional barrier read) | The charged disclination ring: seed A's exact far field (same unit winding sector), core opened into a half-disclination ring (radius a = 4, escaped interior), pinned, working point; `seed_ring` + `core_locus` in the script | E_ring vs E_A endpoints = does the electron ground state prefer a ring core (the m5_particle_hunt synthesis nuance + the M5.16 Q8 off-origin melt, made testable); the barrier read drops unless time allows |
| Audit | independent agent tries to REFUTE: minima genuineness (stationarity), gate reproduction, census numbers | verdicts per claim recorded in the findings |

**Preconditions**: numpy (vectorized einsum, float64) suffices at N ≤ 64; no GPU dependency. Model/effort: Fable 5 / high (research default). `resets_at` = 13:50 EDT is the RELEASE time, not a deadline: checkpoint eagerly, work normally.

## FINDINGS (2026-07-17)

Full record + equations + audit: [`../findings/m5_21_2_census.md`](../findings/m5_21_2_census.md) (the method note for this task). Headline table:

| Question | Answer | Where |
| --- | --- | --- |
| Is the 3D hedgehog topologically protected? | ✅ YES and box-fed: N = 32 free-boundary drains at every δ (his "small box + huge delta" warning measured), N = 48 free survives at probe depth but DRAINS by 24k iters, N = 64 free ≈ pinned at probe depth (retention 0.93-0.97); protection becomes intrinsic only as the box grows | § 4, § 5a |
| Does δ = 0.3 need reducing? | Not at N ≥ 48 for probe-depth work; deep free-boundary work needs bigger boxes than N = 48 | § 4 |
| The 3-lepton signature? | ✅ three DISTINCT levels, electron family lowest, A < C < B, ordering CROSS-INSTRUMENT ROBUST (fwd ×13/×4; 2h re-read 1.64×/2.65×); honest caveat: orders lattice-contaminated energies, not converged continuum minima | § 5c, § 8 C3 |
| The charged ring (third type, P3)? | 🔶 instrument-limited TIE with the point hedgehog (−3.7% fwd vs +23% 2h re-read); ring stays a live co-candidate; source + construction + tests documented | § 9, § 8 C4 |
| Does a compact electron converge at toy parameters? | ❌ NOT with this term set: 2h deep runs are checkerboard artifact (self-caught, demoted); fwd runs are clean but not grid-converged (audit subsample probe); the virial balance is BRACKETED (0.60 < 1 < 2.47 between wscale ×100 and ×1) but never landed | § 5b, § 5c, § 8 C6/C7 |
| THE INSTRUMENT FINDING | at toy parameters the quartic commutator functional has NO stencil-consistent minimizer in reach: each stencil's descent hides curvature in its blind directions (×7-128 cross-stencil disagreement); sharpens Q25 (the term set may need the missing ingredient before the 3D statics is well-posed) | § 5c, § 8 C6 |

## ARTIFACTS + CLEANUP (>1 MB rule)

| Artifact | Status |
| --- | --- |
| `../scripts/m5_21_2_a_scan3d.py` (stack, env knobs M5212_STENCIL/WMULT/TAG/MAXIT) · `m5_21_2_b_topo_illustration.py` · `m5_21_2_c_stencil_check.py` · `m5_21_2_audit_check.py` (audit agent's) | kept |
| `../data/m5_21_2_scan3d.json` (gates + ladder + census after collect) · per-run `m5_21_2_lad_*.json` / `m5_21_2_scanrow_*.json` · `m5_21_2_stencil.json` · `m5_21_2_audit.json` | kept (all < 1 MB) |
| Plots: `m5_21_2_panel.png`, `m5_21_2_topo_catalog.png`, films ×18 (`m5_21_2_film_{basic,thermal}_*`) | kept (plots exempt) |
| Nine endpoint fields `../data/m5_21_2_end_*.npz` (1.7-2.5 MB each, ~19 MB) | DELETED at FINISH (audit complete). Regen (deterministic descents, from `research/scripts/`): 2h era: `python3 m5_21_2_a_scan3d.py scan <A\|B\|C> pinned 48 0.3` and `scan A free 48 0.3` (24k iters default); fwd era: `M5212_STENCIL=fwd M5212_TAG=_fwd M5212_MAXIT=8000 python3 m5_21_2_a_scan3d.py scan <A\|B\|C\|R> pinned 48 0.3 4.0`; stiffness arm: `M5212_STENCIL=fwd M5212_TAG=_w100fwd M5212_WMULT=100 M5212_MAXIT=8000 python3 m5_21_2_a_scan3d.py scan A pinned 48 0.3` |

## DEVIATIONS LOG

(running; entries added as they happen)

| When | Entry |
| --- | --- |
| 10:55 | USER AMENDMENT (not a deviation): P3 redefined from the optional barrier read to the charged-ring arm (the topology-catalog discussion while P2 ran: the ring is the one same-winding-sector alternative that could beat the hedgehog on stability). Ring seed verified (far lams (0, δ, 1) exact; core locus ρ ≈ 3.8; boundary values within 0.9% of A's, noted for the audit) and launched pinned at the working point (bg bqtpl97dw). The A/B/C census rows predate `core_locus` (added with the ring); it will be post-computed from their endpoint npz at collect |
| 11:20 | USER AMENDMENT 2: a static-viewing launcher xperiment for the charged ring (production, physics-only, user-requested): `seed_charged_ring_M` kernel in `engine1_seeds.py` (production frame convention: δ on e_Θ; the sandbox census family, δ on the azimuth, noted in the docstring), `charged_ring` branch in `_launcher.py` (no relax, no 4D activation), `xparameters/_topo_charged_ring.py` ("Charged Ring (static)", boots PAUSED, viewing-only warning in the docstring). Headless-validated at 23³: time row exact, far spectrum (0, 0.29, 1), far radial alignment 0.9995, director +z at origin/interior (escaped), radial far. Cord-shell spectrum not fully isotropic at the smoke grid (a = 2.3 < melt 3.0, a small-grid overlap; clean at the launcher's 64³ where a = 6.4) |
| 11:50 | USER-REPORTED GGUI ARTIFACT on the viewing xperiment (below-plane EM-div/deviation weirdness) → diagnosed as the DIRECTOR-representation seam of the half-integer ring (non-orientable around the cord; EM div/curl + div-colored glyphs are director-derived, engine3 lines 205-237): the first branch-cut choice hung the seam on a half-cylinder at ρ = a, z < 0. FIXED same day: cuts rotated so the seam sits on the minimal spanning disk (ψ = ½[atan2(−z, ρ−a) + atan2(−z, ρ+a)] + π/2); verified 410/410 high-div voxels on the disk, 0 below-plane, far alignment 0.9993, tensor float32-exactly mirror-symmetric under S M(−z) S (3e-7). Side lesson: my first mirror check compared raw components (false-alarmed 0.70 = the sign-flipping xz off-diagonal); conjugate by the mirror operator |
| 12:00 | AUDIT LANDED: 8/8 CONFIRMED; corrections adopted (A_free ratio 1.65e6 outside the quoted pinned range + its milder sawtooth 2.07; cross-stencil factors ×7-128 not ×10-100; the C3 margins caveat: 2h re-read gaps only 1.64×/2.65×, the ordering orders lattice-contaminated energies). The audit's subsample probe STRENGTHENED C6 (coarse-sampled fwd endpoints reproduce the 2h readings) |
| 11:05 | **THE CHECKERBOARD CATCH (major, self-caught from our own film before the audit)**: the A-pinned basic film shows interior single-cell checkerboard from it ~2000; the compact-stencil re-read ([`../scripts/m5_21_2_c_stencil_check.py`](../scripts/m5_21_2_c_stencil_check.py)) measures E_u^fwd/E_u^2h = 10478 (A pinned), 4014 (B), 3299 (C), 1.65e6 (A free) vs 1.26 at the seed; sawtooth ratios 2.4-3.5 vs 1.02. ALL FOUR deep 2h-endpoints are null-mode artifact states; the 24k census E values + ordering are DEMOTED to instrument-artifact status. This GENERALIZES the M5.21.1e anti-recipe: on this functional any deep descent over the 2h central stencil with a weak pointwise potential funnels into the odd-even null mode (the M5.21.1 4D free-FIRE endpoint measured clean at 1.07 there, so the 4D-era results are not automatically contaminated; flagged for a retro-check anyway). ACTION: both running 2h arms stopped (ring, w100); `STENCIL=fwd` compact-difference mode added (exact adjoint; complex-step 2.1e-15, SO(3) 2.2e-16, vacuum 0 re-gated); the census RELAUNCHED on fwd: A/B/C/R pinned + A wscale×100, 8k iters each (bg bjxxxw3me) |

## TASK REVIEW (2026-07-17)

**Task Duration:** 1:55 (from 09:37 to 11:32 EDT)
**Usage Cap Triggered:** NO

**Results** (approved by the user same day):

| # | Result | Status |
| --- | --- | --- |
| 1 | 3×3 3D stack built + gated, both stencils (complex-step ≤ 5.2e-14; SO(3) ≤ 7.6e-16 + negative controls; vacuum 0) | ✅ |
| 2 | Topological protection real and box-fed (N32 free drains at every δ; N48 free survives probe depth, drains by 24k; N64 free ≈ pinned 0.93-0.97); δ = 0.3 kept at N ≥ 48 | ✅ |
| 3 | The 3-lepton signature: three distinct levels, electron lowest, A < C < B, cross-instrument robust (fwd ×13/×4; 2h re-read 1.64×/2.65×); orders lattice-contaminated energies, not converged minima | ✅ with the audit caveat |
| 4 | The checkerboard catch (self-caught pre-audit): 2h deep census demoted (E_u ratios to 1e4); fwd instrument built + re-gated; census re-run clean | ✅ |
| 5 | THE INSTRUMENT FINDING: no stencil-consistent minimizer at toy parameters (×7-128 cross-stencil disagreement; audit subsample probe confirmed the period-2 structure is real); sharpens Q25 | ✅ (interpretation ⚠️) |
| 6 | The charged ring (third defect type, user mid-run addition): instrument-limited TIE with the hedgehog (−3.7% fwd vs +23% 2h re-read); live co-candidate; Alexander RMP 2012 source filed | 🔶 |
| 7 | Virial balance bracketed: 0.60 (w×100) < 1 < 2.47 (base); never landed | ✅ |
| 8 | Independent adversarial audit 8/8 CONFIRMED, corrections adopted | ✅ |
| 9 | Extras: topology-catalog illustration; "Charged Ring (static)" launcher xperiment + the same-day director-seam fix (Dirac sheet → spanning disk) | ✅ |

**Issues / blockers**: the deep-statics program is instrument-blocked (no grid-converged endpoints on either stencil); ISSUE (user catch at review): the census films used 3×3-ADAPTED panel strips, NOT the standard basic/thermal templates of `m5_film.py` (which are 4×4-stack-specific): non-compliant with the series film rule; fix assigned to M5.21.2b (a 4×4-embed adapter so the true templates render 3D states; the endpoint states regenerate there).

**Deviations from plan**: four, logged live in § DEVIATIONS LOG (P3 ring redefinition; the checkerboard catch + fwd re-run; the launcher xperiment + seam fix; audit corrections).

**Action needed**: [M5.21.2b](m5_21_2b_task_details.md) staged (the well-posed instrument + term-set discrimination + converged census, films on the true templates); [M5.21.3](m5_21_3_task_details.md) re-gated on 2b's converged minima. User decisions at review: NO Duda checkpoint yet (consolidate more content across runs first); the Q25 sharpening waits with it.

**Findings**: the 3D-first prescription delivered the qualitative 3-lepton signature (three distinct rotation-family levels, electron lowest, ordering robust across instruments) and confirmed topological protection as real but box-fed; no compact electron converges at toy parameters because the quartic commutator functional currently has no stencil-consistent minimizer (each discretization hides curvature in its own blind directions; self-caught, audited 8/8): the strongest concrete sharpening of Q25 to date. The charged ring, added mid-run as the third defect type, ties the point hedgehog within instrument limits and stays a live co-candidate for the electron core.

**Research docs created / updated**: this file · [`../findings/m5_21_2_census.md`](../findings/m5_21_2_census.md) (the method note) · [`../../theory/_CITATIONS.md`](../../theory/_CITATIONS.md) (+ the filed Alexander RMP 2012 PDF) · scripts [`m5_21_2_a_scan3d.py`](../scripts/m5_21_2_a_scan3d.py) · [`m5_21_2_b_topo_illustration.py`](../scripts/m5_21_2_b_topo_illustration.py) · [`m5_21_2_c_stencil_check.py`](../scripts/m5_21_2_c_stencil_check.py) · [`m5_21_2_audit_check.py`](../scripts/m5_21_2_audit_check.py) · data JSONs (gates / ladder / census / stencil / audit) · plots [`m5_21_2_panel.png`](../plots/m5_21_2_panel.png) + [`m5_21_2_topo_catalog.png`](../plots/m5_21_2_topo_catalog.png) + films ×18 · production [`engine1_seeds.py`](../../engine1_seeds.py) + [`_launcher.py`](../../_launcher.py) + [`_topo_charged_ring.py`](../../xparameters/_topo_charged_ring.py)
