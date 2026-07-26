# M7.2, reproduce Fleury's torus on the lattice (quadrature, the trust-builder)

> Task **M7.2** (M7 / HydroBoros). taskID = M7.N iteration. Status: **Done** (2026-07-02, all gates pass) · Roadmap: [`../m7_roadmap.md`](../m7_roadmap.md)

This doc is the task's full record: planning + findings + future planning + documentation. **M7.2 reproduces the dos Santos & Fleury electron ([arXiv:2510.22384](https://arxiv.org/abs/2510.22384), corpus #1) on the M7 lattice.** It is **pure quadrature**: the ansatz is analytic, the lattice evaluates it and integrates observables. **No minimizer is needed, so M7.2 can run in parallel with the M7.1 minimizer build** (it needs only the lattice, the field containers, and the integration kernels). This is M7's equivalent of M5.11's "reproduce Faber first" trust-builder.

---

## 1. What the paper actually does (read 2026-07-02, full text)

The ansatz satisfies **Maxwell-with-sources where the sources are read off the fields** (the "geometrical charge" ontology): Gauss-B and Faraday are verified (Faraday forces `ω = 2c/R₀`, monochromatic), and Ampère-Maxwell *defines* `J` from the fields; the consistency check is `∇·J` (from Ampère) == `−∂ρ/∂t` (continuity). So the lattice job is:

```text
1. evaluate  E, B  from the ansatz on the grid          (complex fields; physical = Re[·])
2. compute   ρ = ε₀ ∇·E ,   J = (∇×B − c⁻² ∂ₜE)/μ₀     (finite differences)
3. integrate Q_rms, μ, L (spin), U over the torus       (the paper's four constraint observables)
4. compare to the paper's solved values                  (all dimensionless ratios)
```

## 2. Conventions contract (the reproduction traps, pinned)

| # | Trap | Contract |
| --- | --- | --- |
| 1 | **The Heaviside is inverted**: the paper defines `H(r−r₀) = 1 for r < r₀`, 0 outside (their Eq 1 text), the opposite of the standard `H` | implement their stated definition (fields live INSIDE the tube), never the formula's face value |
| 2 | **Instantaneous charge integrates to zero** (`ρ ∝ sin(φ−ωt)` over one wavelength around the ring) | the electron charge is `Q_rms` (RMS over volume), their Eq 19: `√2 π² ε₀ E₀ r₀² = e`; compute RMS, not the plain integral |
| 3 | **Complex-field energy convention**: the paper's `u(t) = ε₀E₀²(1 + R/4R₀)` (Eq 31) comes from the appendix's specific convention | pin the convention against the **FLDB Appendix** (corpus #1, Appendix PDF) before comparing `U`; if ambiguous, compute both real-part-instantaneous and complex-modulus averages and report which matches Eq 32; tracked as **Q10** in [`../m7_question_tracker.md`](../m7_question_tracker.md) (a direct ask to Marc Fleury, the author) |
| 4 | **The mask's delta shell**: `∇·E` at the tube boundary is a distribution; the lattice smears it | grid-convergence study is mandatory (§ 4); report observables vs `h` with Richardson extrapolation |
| 5 | `ω` is not free | enforce `ω = 2c/R₀` (Faraday); phase velocity `v_p = 2c` is a check, not a violation |

## 3. Reproduction targets (paper Eqs 37-41, 43; all dimensionless)

| Observable | Target |
| --- | --- |
| amplitude | `E₀/E_S ≈ 0.286` (Schwinger scale) |
| major radius | `R₀/r_c ≈ 1.573` (`= π/2`) |
| minor radius | `r₀/r_c ≈ 0.152` (`= √(πα)`... per Eq 39) |
| rest energy | `U/m_ec² ≈ 0.795` (`= 5/2π`) |
| frequency | `ω/ω_D ≈ 0.636`, `ℏω/m_ec² ≈ 1.27` |
| charge | `Q_rms = e` (constraint Eq 19 reproduced by the lattice integral) |
| magnetic moment | `μ = μ_B(1 + α/2π)` (constraint Eq 21) |
| spin | `L_z = ℏ/2` (constraint Eq 29) |

Procedure: set `(E₀, R₀, r₀, ω)` to the paper's solved values in `r_c = 1` natural units (M7.1 units contract) and verify the lattice integrals return the four constraints. This validates our integration kernels against a published result, exactly what "reproduce Faber" did for M5.11.

## 4. Grid plan

The tube is thin (`r₀ ≈ 0.152 r_c` vs `R₀ ≈ 1.573 r_c`): resolving it needs `h ≲ r₀/8 ≈ 0.02 r_c`. Box `~4 r_c` across → `~200³` at the finest level; run `h ∈ {r₀/4, r₀/8, r₀/16}` and Richardson-extrapolate. Report every observable as (value, convergence order, extrapolated value).

## 5. Stretch: the Bessel-envelope variant (first beyond-the-paper result)

Fleury § 5.5 names the fix for the unphysical Heaviside mask: Bessel envelopes. Ceperley's rotating mode ([`../../theory/ceperley_rotating_waves.md`](../../theory/_CITATIONS.md) § 4, Eq 15) is the same phase structure with a `J_m(k_c r)` radial envelope. The stretch: replace the mask with the `J_1` envelope, re-integrate the four observables, report how the constraint solution shifts (does `U/m_ec²` move toward 1?). Cheap (same quadrature machinery), and it is M7's **first new result** beyond the published paper.

## 6. Gates

| Gate | Criterion |
| --- | --- |
| primary | all § 3 targets recovered at demonstrated grid convergence |
| consistency | `∇·J` (Ampère) == `−∂ρ/∂t` (continuity) on the lattice to grid accuracy |
| honest limits restated | the mask is unphysical at the boundary, far field suppressed, `U = 0.795 m_ec²` short of 1: reproduced AS the paper states them, not repaired silently |
| stretch | the Bessel-envelope observables table |

Artifacts: `scripts/m7_2_fleury_torus.py` + `data/m7_2_*.npz` + `plots/m7_2_*.png` (field sections, convergence curves, observables table).

---

## FINDINGS (2026-07-02, the quadrature run)

Code: [`../scripts/m7_2_fleury_torus.py`](../scripts/m7_2_fleury_torus.py) · data [`../data/m7_2_observables.json`](../data/m7_2_observables.json) · plots [`../plots/m7_2_field_sections.png`](../plots/m7_2_field_sections.png), [`../plots/m7_2_convergence.png`](../plots/m7_2_convergence.png). Every convention was pinned against the FLDB Main + Appendix full texts before coding; the run takes 5.8 s.

### 1. The paper's printed solution, reconstructed digit-for-digit

Three constraint solutions coexist; the paper's § 3.7 printed values (Eqs 37-40) are **exactly thin-torus + the Schwinger factor** `s = α/2π` (no `r₀²/R₀²` corrections): `R₀ = (π/2)(1+s) = 1.5726`, `r₀ = 0.15141(1+s) = 0.1516`, `E₀/E_S = 0.28658/(1+s)² = 0.2859`, `U = 0.79577/(1+s) = 0.7949`. Every printed digit lands.

| Solution | `E₀/E_S` | `R₀/r_c` | `r₀/r_c` | `U_Eq32/m_ec²` |
| --- | --- | --- | --- | --- |
| thin-torus exact | 0.28658 | 1.57080 | 0.15141 | 0.79577 |
| paper printed = thin + Schwinger | 0.2859 | 1.5726 | 0.1516 | 0.7949 |
| **full corrections (this run's basis)** | 0.28789 | 1.56529 | 0.15106 | 0.79707 |

### 2. Lattice reproduction (the trust-builder gate) ✅ measured

Interface-weighted voxel quadrature (linear boundary reconstruction), grids `h = r₀/4, r₀/8, r₀/16`, fitted-order extrapolation:

| Observable | Extrapolated | Closed form | Rel err | Order |
| --- | --- | --- | --- | --- |
| `Q_rms` | 0.302864 | 0.302822 (= e) | 1.4e-4 | 2.46 |
| `μ` | 0.151609 | 0.151587 (= μ_B(1+α/2π)) | 1.5e-4 | 2.49 |
| `\|L_z\|` | 0.500072 | 0.500000 (= ℏ/2) | 1.4e-4 | 2.48 |
| `U` (Eq 31 density, as published) | 0.797183 | 0.797069 | 1.4e-4 | 2.47 |
| `U` (corrected convention, § 3 below) | 0.958028 | 0.957892 | 1.4e-4 | 2.47 |

Machinery gates (stencil column): bulk `ρ = ∇·E` vs analytic `−iE₀e^{iφ}/R₀` at 1.2e-5 (O(h²)); continuity `∇·J = iωρ` at 2e-14 (a discrete identity, machine); Faraday residual 7.4e-6 (O(h²)); instantaneous total energy == period average at every sampled phase; `L_z < 0` as the left-handed triad demands (Eq 23's `S ∝ −φ̂`), magnitude gated. The closed forms themselves were re-verified by independent 2D quadrature to 2.3e-16.

### 3. The energy-convention finding (the Q10 evidence) ✅ measured, interpretation pending Marc

The appendix's own convention is the standard phasor average `E²(t) = ½E·E*` (its Eq 113-115), but **Eq 122/124 drop the square on `(1+R/R₀)` in `E_φE_φ*`, and Eq 127 drops the ½ on the B term**; the two slips produce Eq 31's `u = ε₀E₀²(1+R/4R₀)`. Applying the appendix's own convention exactly:

```text
u      = (ε₀E₀²/4) [2 + (1+R/R₀)²]
U_phys = 3π²ε₀R₀r₀²E₀² (1 + (5/24) r₀²/R₀²)  =  (6/5) × Eq 32   (thin torus)
```

The three constraints (Q, μ, L) are slip-free, so the solved parameters stand and only the energy PREDICTION moves: **`U/m_ec² = 0.795 → 0.958`**, closing most of the paper's acknowledged gap to `m_ec²` (its § 3.11). Both densities were integrated on the lattice against their closed forms (table above). Status: the algebra is ✅ measured (three independent routes agree: closed form, 2D quadrature, 3D lattice); calling it a slip in the paper is 🔶 pending Marc's confirmation (Q10, [`../m7_question_tracker.md`](../m7_question_tracker.md)).

### 4. Bessel-envelope stretch: the mask's hidden surface charge, exposed ✅ measured

The naive scalar `J₀(j₀₁s/r₀)` envelope (the literal reading of the paper's § 5.5 fix) is **not viable**, for a physical reason worth recording:

| Quantity | Value | Reading |
| --- | --- | --- |
| `Q_rms` at paper params, Bessel | 2.78 = **9.2×** the Heaviside `Q` | dominated by the envelope-gradient charge |
| gradient term vs core term | 2.77 vs 0.16 | `ρ ∝ g′(s)E₀ ~ (j₀₁/r₀)E₀` swamps the model's `E₀/R₀` |
| the Heaviside **shell** charge (column B) | 3.32 | the SAME object: the sharp mask hides it as a boundary delta the paper's bulk integrals drop; smoothing makes it explicit and un-droppable |
| Faraday residual, Bessel | 19.8 | the naive envelope is nowhere near monochromatic |
| constraint re-solve | RUNAWAY (`R₀ → 33` at iter 0) | no nearby smooth solution of the naive scalar-envelope form |

Conclusion: Fleury § 5.5's Bessel fix requires the full per-component mode structure (Ceperley's construction) or a relaxation of the full functional, which is exactly M7.4's job. This sharpens the paper's own § 5.2 honesty note into numbers: the "unphysical mask" is quantitatively a hidden surface charge ~18× the bulk RMS charge.

### 5. Honest limits restated (gate requirement)

As the paper states them: the mask is unphysical at the boundary (now quantified, § 4); the far field is suppressed (fields identically zero outside the tube, no Coulomb tail); `U` falls short of `m_ec²` (0.795 as published; 0.958 under the corrected convention, § 3, still 4.2% short). Nothing was repaired silently; the published numbers are reproduced AS published alongside the corrected column.

---

## TASK REVIEW (2026-07-02)

**Task Duration:** 00:31 (from 19:21 to 19:52 EDT)
**Usage Cap Triggered:** NO

**Results**: all 6 gates ✅ measured (full tables in § Findings above): the Fleury electron reproduced on the lattice (`Q/μ/L/U` to closed forms at 1.4e-4, order ~2.5; machinery gates machine-exact to O(h²)); the paper's printed solution reconstructed digit-for-digit (thin-torus + Schwinger); the Q10 energy finding (Eq 122/124/127 slips → corrected `U = 0.958 m_ec²`, algebra measured three ways, attribution 🔶 pending Marc); the Bessel stretch exposed the mask's hidden surface charge (an honest negative with a physical reading).

**Issues / blockers**: none; two implementation bugs (a 2π in μ, erosion starving thin-tube integrals) caught by the gates and fixed via the interface-weighted quadrature.

**Action needed**: Q10 ask sharpened for Marc (tracker updated); the § 3 finding corrects the authors' paper, outbound phrasing is the maintainer's. Next: M7.3 (verbatim-ODE pre-gate + the 3D chaoiton); after M7.3, the consolidated report + ask round to Marc.

**Findings**: The lattice reproduces every published observable of the Fleury electron at demonstrated O(h²·⁵) convergence (the "reproduce Faber" trust-builder passed), and finds two results beyond the paper: its energy formula rests on two identifiable algebra slips (corrected with the paper's own convention, `U` rises 0.795 → 0.958 m_ec² at unchanged parameters), and the "unphysical mask" is quantitatively a hidden surface charge ~18× the bulk RMS charge, pointing at M7.4's relaxation as the principled fix.

**Research docs created / updated**: this doc (§ Findings) · script [`../scripts/m7_2_fleury_torus.py`](../scripts/m7_2_fleury_torus.py) · data [`../data/m7_2_observables.json`](../data/m7_2_observables.json) · tracker [`../m7_question_tracker.md`](../m7_question_tracker.md) (Q10) · roadmap row → Done.

### Plots

![`../plots/m7_2_convergence.png`](../plots/m7_2_convergence.png) (key plot)

![`../plots/m7_2_field_sections.png`](../plots/m7_2_field_sections.png)

---

Cross-refs: roadmap [`../m7_roadmap.md`](../m7_roadmap.md) (M7.2) · background [`../m7_background.md § 3`](../m7_background.md) (the Fleury model + targets) · conventions sources: FLDB Main + Appendix (corpus #1), [`../../theory/ceperley_rotating_waves.md`](../../theory/_CITATIONS.md) (Bessel stretch), [`../../theory/feynman_maxwell_equations.md`](../../theory/_CITATIONS.md) (the Maxwell baseline) · upstream [`m7_1_infra.md`](m7_1_infra.md) (lattice + units; minimizer NOT required) · downstream [`m7_4_charged_soliton.md`](m7_4_charged_soliton.md) (the Fleury seed).
