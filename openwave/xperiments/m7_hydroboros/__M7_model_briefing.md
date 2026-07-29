# M7 HydroBoros: Model Briefing

> **What M7 brings.** The rigorous full-3D-PDE toroidal Beltrami electron that neither parent
> built: it fuses Fleury's toroidal-EM electron (analytic) with Werbos's Ouroboros
> self-confinement (1D radial) on OpenWave's M5-proven Taichi lattice, so the electron's
> field configuration is both *specified* and *earned* as the energy-minimizer. **Phase 1
> is complete (M7.0-M7.7, 2026-07-04)**: the canonical spec is
> [`research/m7_theory_canonical.md`](research/m7_theory_canonical.md) and the M7
> column is PUBLISHED in
> [`MODELS.md`](../../../MODELS.md) (M7.21, 2026-07-18); its position and its icon
> tally are the score-board's, read there.

![M7_hero](research/images/hydroboros_icon.jpg)

## Identity

| Field | Value |
| --- | --- |
| Model ID | M7 |
| Name | HydroBoros (Hydrodynamics + Ouroboros; evokes the Hydra water-snake) |
| Author | Marc Fleury (toroidal-EM electron) + Paul Werbos (Ouroboros / M6): the two physics parents |
| Author contact | GitHub [@marcf999](https://github.com/marcf999) (toroidal-EM electron) + [@pwerbos](https://github.com/pwerbos) (Ouroboros), for author-gated questions (definitions, intent, what the model does and does not claim); routing convention in [`dev_docs/CROSS_MODEL_TESTING.md`](../../../dev_docs/CROSS_MODEL_TESTING.md) § 6. Note: this column is PARKED by its own roadmap pending the model author's return from conference, so author-gated rows here record "not applicable" under the § 1 default in the meantime |
| Blend | Fleury's toroidal EM electron (arXiv:2510.22384) fused with Werbos's Ouroboros self-confinement (M6) |
| Lineage | force-free / Beltrami (Trkalian → variable-λ) + knotted-EM / Clebsch + Faber geometric soliton |
| Primary sources | Fleury / dos Santos arXiv:2510.22384; Werbos M6; Faber arXiv:2201.13262 + Faber & Golubich [Phys. Rev. D 114, 014510](https://doi.org/10.1103/8zn4-rwth) (2026; preprint arXiv:2604.12021); Sato-Yamada arXiv:1809.03136; Ceperley; Pisello 1977 |
| In-repo | [`research/m7_roadmap.md`](research/m7_roadmap.md) (roadmap) + [`m7_background.md`](research/m7_background.md) (background) + [`m7_question_tracker.md`](research/m7_question_tracker.md) (Q1-Q7); `theory/` (66-doc electron-Beltrami corpus + notes); `images/` icon |

## Model Profile (what it brings, short form)

| Attribute | M7 |
| --- | --- |
| Substrate | A-primary ontology; leading candidate = the Ouroboros doublet `(A_μ, J_μ)` read via Riemann-Silberstein (diagnostic); target manifold S² / S³ open (Q1); worked in the **time-harmonic (fixed-ω) frame**, both parents are time-periodic |
| Particle | electron = self-linked toroidal Beltrami vortex (approximately Beltrami: the relaxed minimizer of the full functional) |
| Charge | the **measured** divergence `∇·F` the Ouroboros coupling drives, unified with Chern-Simons linking / helicity `∫A·B` (the Q3 experiment) |
| Stabilization | **helicity anti-collapse** (Arnold bound) + **Ouroboros confinement anti-expansion**; Nadirashvili's theorem forces the blend (no finite-energy pure-Beltrami field exists); a 4th-order term is optional (Q2, off by default) |
| Energy functional | period-averaged Maxwell kinetic (both fields) + Ouroboros coupling `m_J²A·J − f(J·J)` (Werbos), minimized at fixed ω + helicity |
| Solve method | fixed-ω minimization (FIRE / L-BFGS) + fixed-helicity relaxation + Taichi reverse-mode AD gradient (vs numpy FD to ~1e-12) + constrained Minkowski leapfrog; the M5.11 machinery in the time-harmonic frame |
| Free parameters | `m_J`, `g` (M6 canonical g = 1.0; Werbos-v5 point g = 1.0625, dictionary = Q9); ω swept for the spectrum |
| Reproduction targets | Fleury `U ≈ 0.795 m_e c²`, `ω = 2c/R₀`; M6 `H/Q = 1.6890` at the same `(g=1.0, ω=1)` calibration, in full 3D; Faber & Golubich α⁻¹ ≈ 137 (M7.6 stretch) |

## Field Configuration of Particles

The electron's field configuration *is* the self-linked toroidal Beltrami vortex (a
topological vortex).

| Test | M7's answer | Phase |
| --- | --- | --- |
| Coulomb | two charge configs → read the `1/r` interaction energy `E(d)` | M7.4 (single-charge `1/r`) + M7.6 (two-charge `E(d) ~ 1/d`) |
| Clock | the de Broglie frequency = the energy-minimizing one | M7.5 (= the M5.8 mechanism) |

Other particles ride the same substrate later: the neutral knot = dark matter (M7.20), the
lepton / neutrino family (M7.19), quarks / baryons / mesons (M7.22).

## Implementation Status

**Phase 1 complete (M7.0-M7.7, 2026-07-04); the M7 column is PUBLISHED in [`MODELS.md`](../../../MODELS.md)** (entered at M7.21 (2026-07-18) as the 3rd column at 21 cells, 0 ✅ / 8 ⚠️ / 13 🚧; the shared criteria set has since grown to 31, so read the live score-board for today's icons)<!-- count:historical -->. Both parents
reproduced from one lattice code (Fleury's torus at 1.4e-4, M7.2; M6's charged `H/Q` in
full 3D at 4.7e-5, M7.3), then the new physics earned: the first stable finite-size 3D
soliton family (Taylor-dressed, `E = 0.802\|H_A\|`, M7.4), the vacuum-tachyon discovery +
the `ω* = 0.786` existence threshold (the clock IS the stabilizer, M7.5), and the rotating
`j_z = 1` electron with the Coulomb sector (Gauss 99.1%, two-charge `1/d` reference-matched,
1.17 dressing measured, M7.6). Canonical spec + one-script reproduction:
[`research/m7_theory_canonical.md`](research/m7_theory_canonical.md) + `research/scripts/m7_7_canonical.py`
(physics module: `research/scripts/m7_functional.py`). Honest open items: the vacuum
tachyon (Q14, the top theory question), the self-consistent charge, absolute μ, the
units contract (ℏ/2 vs ℏ).

## Roadmap

Five phases, M7.0 → M7.23 (full detail in [`research/m7_roadmap.md`](research/m7_roadmap.md)):

| Phase | Tasks | What lands |
| --- | --- | --- |
| 1 , electron + the column ✅ DONE | M7.0-M7.7 | bootstrap → infra → reproduce Fleury → reproduce M6 in 3D → charged soliton (the new physics) + Coulomb → clock + stability → observables → consolidate (milestone; column drafted, published to [`MODELS.md`](../../../MODELS.md) at M7.21) |
| 1 extension , the author-collaboration track ✅ | M7.8 ✅, M7.9 ✅, M7.10 ✅ + reserved M7.11-M7.14 | the helicity-pair 3:1 test + the Phase 1 walkthrough (M7.8 ✅ 2026-07-07: the closure pair not stationary, spin reads one quantum); the ChaosBook study + canonical-exercise benchmark (M7.9 ✅ 2026-07-07: 5/5 benchmarks incl. the full Rössler cycle table 14/14, audit CONFIRMED, orbit toolkit ready); the pure-Maxwell no-Lagrangian test (M7.10 ✅ 2026-07-07: Theorem 2 at the integrator floor, the honest boundary as measurements (walls 0.371 vs free space 0.042; one state, two deaths), the Q14 tachyon attributed to the bilinear via the 0.786√ε ladder at 0.1-0.3%, audit CONFIRMED); the Beltrami orbit track: [M7.11 defined](research/tasks/m7_11_variable_lambda.md) (the variable-λ vacuum-cure check, top of backlog), M7.12-M7.14 reserved (post-August); **the walkthrough § 7 report is COMPLETE** and **M7 is PARKED** until the author returns from the conference and takes an adversarial external pass (sequencing: [M7.11 § 0](research/tasks/m7_11_variable_lambda.md)) |
| 2 , forces + sectors | M7.15-M7.20 | magnetic, gravity (hard), nuclear (strong / weak), antimatter / annihilation, lepton + neutrino family, dark matter |
| 3 , MODELS.md publication ✅ | M7.21 ✅ | published the 21-cell HydroBoros column to MODELS.md (2026-07-18; benchmark entry, 3rd column) |
| 4 , composites | M7.22 | quarks, baryons, mesons, orbital quantization |
| 5 , production | M7.23 | graduate the canonical recipe to `medium.py` + engines + `_launcher.py` rendering |

Open questions Q1-Q15 (priority-sorted in the tracker): the charge-carrying construction (Q7, now a
tail-matching + tail-quantization program per the 2026-07-05 Fleury notes; rides the author's
post-August window), the two charges (Q3, measured independent at M7.4; same window), substrate +
target manifold (Q1), source-material logistics (Q4), the vacuum tachyon (Q14, the top theory
question, self-determined), the windowed charged calibration (Q11, self-determine), chaoiton 3D
stability (Q13, self-determine), ODE provenance (Q12), and the v5 calibration dictionary (Q9).
Resolved so far: Q8, Q2, Q5, Q6, Q10 (the FLDB energy convention, author-confirmed 2026-07-04), and
Q15 (the units contract, resolved by author directive at the 2026-07-06 call: target the observable
`S_z = ℏ/2`; M7.8 measures it).

## Help Wanted

M7 is the newest and most open program, actively recruiting. It is the natural home for the
hydrodynamics / Beltrami and toroidal-electron communities.

| You can contribute | How |
| --- | --- |
| Build a phase | pick a sandbox phase (M7.1 infra, M7.2 Fleury reproduction, ...) and land its gate |
| Bring Beltrami expertise | the variable-λ (charge-carrying) construction is the research core (the Spanish Beltrami school is being invited) |
| Contribute a source or seeder | knotted-EM / Bateman / Trkalian seed generators |
| Answer an open question | Q1-Q7, especially the substrate manifold (S² vs S³) and the `λ(x)` profile |

Flow: open an issue or discussion → fork → branch → PR with a DCO sign-off
(`git commit -s`), under Apache 2.0. New-model governance: open an issue first so a
maintainer adds the column at the M7.21 publication milestone (Phase 3). See
[`../../../MODELS.md`](../../../MODELS.md) § Contributing,
[`../../../ONBOARDING_MODELS.md`](../../../ONBOARDING_MODELS.md),
[`../../../CONTRIBUTING.md`](../../../CONTRIBUTING.md). Model discussion runs in the
Models-of-Particles group.

## Rich Context for Deep Reader

This is a top level documentation and orientation content. For additional context on this model, a detailed read in the /theory and /research folders is recommended, as well as the production files in this model root folder, that may contain the canonical full PDE implementation of the theory.
