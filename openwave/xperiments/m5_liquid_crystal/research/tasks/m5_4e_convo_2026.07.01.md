# Convo record 2026-07-01: Duda on serious simulations, the g/δ parameter-lock, and the lepton-vs-neutrino clock slide

Participants: Jarek Duda (Jagiellonian), the maintainer (OpenWave). Thread: preparing for heavier simulations (Fable 5 available again), what "serious" verification requires, and the two basic parameters (`g`, `δ`) that still need pinning. Continues the [`m5_4d_convo_2026.06.11.md`](m5_4d_convo_2026.06.11.md) themes (ZBW energy-minimum, regularization, neutrino SO(3)/δ_CP, gravity-to-propel-oscillations) and the [`m5_4c_convo_2026.06.08.md`](m5_4c_convo_2026.06.08.md) δ/g calibration thread. Slide from Duda in [`../../theory/duda_2026-07-01_lepton_neutrino_clock.png`](../../theory/duda_2026-07-01_lepton_neutrino_clock.png).

> **Why this record matters.** It reframes the calibration program: the honest bar for convincing mainstream is a **static energy-minimization** simulation (lattice/FEM, cylindrical symmetry, a regularizing potential), not a scale-free dynamical clock. That is the design fork behind the new task [`m5_16_task_details.md`](m5_16_task_details.md). Clean for the public repo.

---

## 1. The methodology bar — "serious simulations"

Duda's honest read on what verification actually requires:

| Point (Duda, 2026-07-01) | Verbatim / close paraphrase |
| --- | --- |
| Serious sims are heavy | not seconds but weeks; "Faber said about weeks for running coupling" |
| The method | "needs lattice or FEM, assuming initial field configuration like hedgehog for electron or vortex loop for neutrino, and numerically performing **energy minimization**" |
| The hardest part | "the most difficult is **regularization in the center of hedgehog or vortex**, also requiring potential, which details are still to be established" |
| Why it is required | "without serious simulations I don't think we could convince mainstream e.g. in article" |
| The honest self-assessment | "personally I still don't know what to think about these results" |
| What the current results ARE good for | "sufficient to filter out the bad models, also stimulating the authors for serious verification, discussion, going out of own bubble, comparing with nature and others" |

**Dimension reduction (the enabling trick):** "for both electron and neutrino we can assume **cylindrical symmetry to reduce dimension**." This is what makes the physical-regime run (below) tractable.

**The OpenWave read.** This is a paradigm fork, not a footnote. M5.8 established the clock by **dynamical leapfrog** at V=0, which is scale-free (the `2m` rigidity + the 28× absolute-ω gap, `m5_summary_report.md` §3, N-6b). The "serious" route Duda describes is **relaxation to a minimum-energy configuration** of an axisymmetric ansatz **with** a regularizing potential, which is exactly where the missing length anchor lives (`#218`: the ZBW scale is recoverable only by anchoring energy AND length jointly, the `E·r₀` line). OpenWave already built this minimizer inside M5.11 P0 (`m5_11_p0_minimizer.py`) and it reproduced the electron (511.00 keV at `r₀=2.2132 fm`) + `α⁻¹→137.03`, but only at placeholder `g=8, δ=0.3`. Graduating it to the physical regime is [`m5_16_task_details.md`](m5_16_task_details.md).

---

## 2. The parameter-lock — g, δ, the potential, regularization

Duda's ordered prescription (2026-07-01, second email):

| Step | Duda's statement | Anchor |
| --- | --- | --- |
| **First: fix g, δ** | "the first step should be **establishing two basic parameters: g, delta** — seems we still don't know" | — |
| δ source | "QED Lagrangian suggests delta" | the `ℏc` Dirac-kinetic coefficient (`m5_4c` decoder-ring: `δ` enters as `δ²`) |
| second anchor | "**Coulomb can give another anchor**" | the fixed-units-via-Coulomb axis (NG-1/NG-3) |
| **Then: the potential** | "further simulations also require potential, for which we only know that **minima should be eigenspectrum (g,1,delta,0)**, can be obtained in various ways, like **Landau-de Gennes with traces of powers**, still requiring to find its parameters" | Q7 (the `a·Tr(M²) − b·Tr(M³) + c·(Tr M²)²` coefficients) |
| Time scale | "**electron clock seems a good anchor for time scale**, but requires regularization with potential" | `#220` absolute-ω |
| g source | "**maybe g parameter can be obtained from electron clock, neutrino oscillations**. Otherwise might need **gravitational mass — we are certain of only for baryons**" | `#220` / M5.11 neutrino / the gravity sector |

**Key structural point.** The potential + regularization are not free add-ons to be guessed: their parameters are **fixed by requiring the minimizer to hit the known anchors** (electron 511 keV via Faber `r₀`, Coulomb, the electron clock time-scale). That closes the loop Duda is describing: minimize energy under `V`, tune `(g, δ, a, b, c)` until the anchors land, read out everything else. This is M5.16's core method.

**The numerical obstacle (from M5.11 planning, same regime):** the physical `g~1e10`, `δ~1e-10` is a ~`1e20` dynamic range that exceeds f32/f64, so the build needs a **non-dimensionalized formulation or a perturbative δ-expansion**. Cylindrical symmetry (§1) is the other half of making it run.

---

## 3. The slide — charged leptons vs 3 neutrinos, and the clock topology

[`../../theory/duda_2026-07-01_lepton_neutrino_clock.png`](../../theory/duda_2026-07-01_lepton_neutrino_clock.png). Duda's comparison table + beta-decay cartoon + the clock-topology punchline.

| Property | electron / muon / taon | 3 neutrinos |
| --- | --- | --- |
| size | `e⁻e⁺` scattering: at most ~2 fm | Nature: wavepacket ~6200 fm |
| mass | constant 511 keV | <1 eV?, could vary |
| self | decay into the lightest: electron | oscillates between 3 neutrinos |
| interaction | Coulomb: electric charge quant(ized) | non-interaction: cross Earth |
| clock | fixed angular momentum **SO(2)** | **SU(2)_L ~ SO(3)** oscillations |

- **Beta decay** drawn as a quark-string sequence: neutron → shift → split (energy release) → reconnect quark string → proton + electron + (anti)neutrino. Consistent with the `m5_4d` §2 reading (beta decay = topology-reconnection / defect-class transition, not an EM-style force).
- **The punchline (bottom of the slide):** the **electron de Broglie clock = SO(2)~U(1) 2D rotations** (the twist / `clock_twist`), while **neutrino oscillations = SO(3)~SU(2)_L 3D rotations**. This is the clean group-theoretic split between the charged-lepton single-clock (one phase, one frequency) and the neutrino 3-axis flavour swing (`#199` SO(3) route, the δ_CP = 180° falsifier).
- e/μ/τ are rendered as textured hedgehog tori (family = axis-choice of the biaxial Λ); the 3 neutrinos as vortex-loop / topological-string tori.

**Consistency check with the platform.** The SO(2) vs SO(3) split is the group-theory face of what the model already carries: the charged lepton is the single-`clock_twist` de Broglie clock (M5.8 existence ✅), the neutrino is the δ-0 axis-2↔3 swing (`#199`/M5.11). Nothing here contradicts the record; it sharpens the neutrino target as a genuinely 3D (SO(3)) object vs the electron's 2D (SO(2)) clock.

---

## 4. Where this lands (placements)

| Insight | Home | Action |
| --- | --- | --- |
| Static energy-minimization + cylindrical symmetry = the serious-sim bar; physical-regime g/δ | **new [`m5_16_task_details.md`](m5_16_task_details.md)** | created; graduates the M5.11-P0 minimizer into the calibration instrument; gates the calibration cluster |
| "First step = fix g, δ" via clock / neutrino-osc / baryon gravitational mass | **[`m5_9_0_task_details.md`](m5_9_0_task_details.md)** | cross-ref added; the real-calibration-axis residual now routes through M5.16 |
| Potential minima = (g,1,δ,0); LdG traces-of-powers; parameters TBD | **Q7** in [`../m5_question_tracker.md`](../m5_question_tracker.md) | 2026-07-01 note added |
| Regularization at hedgehog/vortex center = the hardest part | **Q8** in [`../m5_question_tracker.md`](../m5_question_tracker.md) | 2026-07-01 note added |
| Cylindrical symmetry + the parameter-lock for the neutrino vortex loop; SO(3) slide | **[`m5_11_task_details.md`](m5_11_task_details.md)** | cross-ref added; the parameter-lock blocker = M5.16 |

---

## 5. Open follow-ups

- [ ] **Build M5.16** (the physical-regime axisymmetric minimizer) as the recommended next task, ahead of resuming M5.11 (whose own blocker is the same parameter-lock).
- [ ] Decide the non-dimensionalization vs perturbative-δ route for the `1e20` dynamic range (M5.16 §Approach).
- [ ] Report back to Duda once g/δ are pinned by the minimizer (his explicit "please write if I can help").

## Out of scope for M5 (noted, not filed here)

The maintainer's reply mentioned tooling prep for Fable 5 and helping Marc Fleury implement Paul's Ouroboros model in a new **M7** environment. The Ouroboros/M7 thread is an M7 seed, not M5 physics; not captured in this M5 record.
