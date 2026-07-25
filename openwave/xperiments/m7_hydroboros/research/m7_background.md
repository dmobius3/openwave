# M7 / HydroBoros - Background

> **Purpose.** Build a rigorous **full 3D PDE simulation** of the electron as a
> self-linked toroidal vortex, blending the **hydrodynamics school** (Marc Fleury's toroidal
> electromagnetic electron + the Beltrami / force-free field tradition) with the **Ouroboros
> model** (Paul Werbos's coupled two-vector-field chaoiton, M6). HydroBoros = Hydrodynamics +
> Ouroboros, also evoking the mythological **Hydra**, the water snake. Both are built on **M5 (Liquid
> Crystal, Duda) as the third source**: M7 inherits OpenWave's M5-proven lattice + energy-minimizer,
> M5's physics primitives, and M5's development workflow (the [`m7_roadmap.md`](m7_roadmap.md) + [`m7_question_tracker.md`](m7_question_tracker.md)).
>
> The work runs headless in
> `research/` (scripts/, data/, plots/, tasks/; Taichi GPU lattice, energy minimizer, matplotlib diagnostics),
> converges to a canonical form, then graduates into the production `medium.py` + `engine1-4` +
> `_launcher.py` for visual rendering, exactly the path M5 followed.
> A primary deliverable is the **HydroBoros (M7) column** of the repo-root
> [`MODELS.md`](../../../../MODELS.md) coverage matrix, every cell backed by a runnable in-platform
> script, the same bar M5 meets (see the [roadmap](m7_roadmap.md) Phase D).
> Rigor standard: [`../../m5_liquid_crystal/research/11a_vortex_loop.md`](../../m5_liquid_crystal/research/tasks/m5_11a_vortex_loop.md) (M5.11).
> Theory sources: [`../theory/`](../theory/). Findings as we execute: `research/scripts/`, `research/data/`, `research/plots/`, tracked in [`m7_roadmap.md`](m7_roadmap.md).

---

## 0. The sources , two physics parents + the M5 method

| Source | Substrate | How it is solved today | Charge origin | Role in M7 |
| --- | --- | --- | --- | --- |
| **Fleury** (hydrodynamics school), [arXiv:2510.22384](https://arxiv.org/abs/2510.22384) | EM field `E, B` confined to a torus | closed-form **analytic** ansatz + Heaviside mask | `ρ = ∇·E` (geometric / divergence charge) | physics parent; limit: mask unphysical at boundary, energy `0.795 m_e c²`, **no dynamics/PDE** |
| **Werbos Ouroboros** (M6), [`0d_canonical.md`](../../m6_ouroboros/research/archive/0d_canonical.md) | **double vector field** `(A_μ, J_μ)`, **time-periodic** (`cos ωt / sin ωt` pair) | **reduced 1D ODE / BVP** (cylindrical `(α,β)`; spherical `l=1`) | `Q_CS` = Chern-Simons **linking** number | physics parent; limit: never a full **3D** field, electron only a radial profile (`H/Q` ledger in § 4) |

The honest gap the two physics parents share: **neither has ever been evolved as a full 3D nonlinear
PDE on a lattice.** Fleury is analytic; M6 is a 1D radial reduction. That is exactly the gap M5.11
closed for the tensor model (Taichi-AD lattice + FIRE / L-BFGS minimizer, reproduce Faber's electron
from the relaxed 3D field). **M7's contribution is to build the rigorous 3D PDE simulation neither
parent did**, with the toroidal vortex as the soliton sector.

## 1. M5 as a reference structure

**M5 as a reference** (Liquid Crystal, Duda) , the **method + rigor** source | 3×3 / 4×4 tensor field `M = ODO^T` | **full 3D Taichi-AD lattice + FIRE / L-BFGS minimizer** (M5.11) | topological winding (Gauss-Bonnet integer) | M7 borrows the **method + engine + workflow + physics primitives**, not the substrate

**M5 is the third source , for its proven method, not its substrate.** M7 is built on OpenWave's
M5-validated lattice + energy-minimizer and inherits M5's physics primitives (the **4th-order term
that evades Derrick**, the **energy-minimizing de Broglie clock** of M5.8, **Coulomb-from-topology**,
**particle-as-field-configuration**), and adopts M5's development **workflow** (roadmap + GitHub
issues + question tracker + research→canonical→production + AI-reviewer passes + the 11a rigor standard,
detailed in the [`m7_roadmap.md`](m7_roadmap.md) + [`m7_question_tracker.md`](m7_question_tracker.md)). M5 is OpenWave's most validated model (16 MODELS.md cells), so it is the natural template for
getting M7 to the same bar.

| Rigor | Details |
| --- | --- |
| **AI-reviewer pass** | every canonical recipe + MODELS.md cell goes through an AI code-review (`/code-review`) before merge , the quality gate M5 uses |
| **11a rigor standard** ([`m5_11a_vortex_loop.md`](../../m5_liquid_crystal/research/tasks/m5_11a_vortex_loop.md)) | energy-minimizer to `‖∇E‖ → 0`, AD-validated gradient (`1e-12`), each task gated against a known result, honest pass/fail (§ 5, the [roadmap](m7_roadmap.md)) |

**The wider OpenWave wave-physics library (M1-M4) , mine for insights, not a source.** HydroBoros is
fundamentally a *wave* model (rotating waves, EM waves, Beltrami flows, knotted light), and the other
OpenWave models carry a deep body of wave physics a coding agent should consult: **M1**
([`m1_granule_motion`](../../m1_granule_motion/)) , granule-motion wave propagation (standing /
traveling / spherical / energy waves); **M2** ([`m2_free_wave`](../../m2_free_wave/)) , free-wave
dynamics; **M3** ([`m3_wolff_lafreniere`](../../m3_wolff_lafreniere/)) , Wolff-LaFreniere
wave-structure-of-matter (standing-wave centers, interference → particles); **M4**
([`m4_ewt`](../../m4_ewt/)) , Jeff Yee's Energy Wave Theory (wave constants + equations). These are
**not parents** (they do not define HydroBoros's substrate), but their wave-decomposition,
interference, and wave-center primitives bear directly on the M7 **seeders** (task M7.1, [roadmap](m7_roadmap.md)) and the
wave-propagation layer of the PDE , consult them when building those.

## 2. The HydroBoros thesis

The Ouroboros "snake eating its tail" is a **self-linked toroidal vortex**. Fleury's toroidal EM
wave is what that vortex looks like in the hydrodynamic / Maxwell reading. The **Beltrami**
(force-free) condition `∇×F = λ F` is the steady, self-sustaining circulation that ties the two
together: it is simultaneously Fleury's monochromatic toroidal eigenmode (`ω = 2c/R₀`) and the
Ouroboros self-confinement. The bridge is already documented inside M5 (Fleury's Navier-Stokes ≡
generalized-Maxwell equivalence, Duda's superfluid mapping): see
[`../../m5_liquid_crystal/research/1b_topological_defect.md`](../../m5_liquid_crystal/research/tasks/m5_1b_topological_defect.md)
§ "EM-hydrodynamics formal equivalence".

### Trkalian vs Beltrami: where the charge lives (Marc, 2026-06-30)

A Beltrami field `∇×w = λw` splits on whether the proportionality factor `λ` is constant. Taking the
divergence of both sides (`∇·(∇×w) = 0`) gives `∇λ·w + λ ∇·w = 0`, so:

| Case | `λ` | Divergence | Charge |
| --- | --- | --- | --- |
| **Trkalian** | constant | `∇·w = 0` (solenoidal, forced) | none |
| **general Beltrami** | variable `λ(x)` | `∇·w = −(∇λ·w)/λ ≠ 0` | **yes , the charge IS the variable-λ divergence** |

This is the bridge to Fleury's `charge = ∇·E`: the **electric charge is the non-constant-λ Beltrami
divergence**. The **endorsed build approach** (Marc, 2026-06-30, "the way to build it"): use the
**modern decomposition of Beltrami** (Sato-Yamada, [M7.0 corpus](tasks/m7_0_bootstrap.md) #11) and **start the code on Trkalian
(constant-λ) exact solutions, then take off the training wheels on variable-λ** , generalize S&Y to
non-constant λ, find the charge-carrying ansatz, and evolve it. The non-constant case is exactly the
hard part ("the math gets very hairy"), so the Trkalian start is the foothold. This maps directly
onto the tasks: M7.1-M7.3 stay in the charge-free Trkalian/neutral regime, and **M7.4 takes off the
training wheels , the constant-λ → variable-λ transition where the charge `∇·F` first appears** (the [roadmap](m7_roadmap.md), task M7.4).
Marc may bring in the **Spanish Beltrami school** (Enciso & Peralta-Salas, [M7.0 corpus](tasks/m7_0_bootstrap.md) #6) as collaborators.

### What the Beltrami mathematics allows (corpus review, 2026-07-02)

The corpus pins four hard boundaries on the variable-λ program; M7.4 is designed around them
(see the reframe in [`tasks/m7_4_charged_soliton.md`](tasks/m7_4_charged_soliton.md)):

| Fact | Statement | Consequence for M7 |
| --- | --- | --- |
| **Divergence identity** | `∇·w = −(w·∇λ)/λ` (take `∇·` of `∇×w = λw`) | charge lives exactly where λ varies **along** field lines; solenoidal fields force `w·∇λ = 0` (λ constant on field lines) |
| **Rigidity** (Clelland & Klotz 2020, ARMA; Enciso-Peralta-Salas; cited in the Enciso 2023 survey, corpus) | nonconstant-λ Beltrami fields face many obstructions, **even locally**; arbitrary `λ(x)` profiles admit no solution | do NOT hunt an exact charged variable-λ Beltrami ansatz; relax the full functional and **measure** the deviation |
| **Perturbative existence** (Kaiser-Neudert-von Wahl 2000, corpus) | variable-α force-free fields exist for **small** α, near vacuum fields, in exterior domains | small-charge regime is mathematically safe ground; large-charge knots are terra incognita, exactly M7.4's honest experiment |
| **No finite-energy Beltrami** (Nadirashvili; Enciso 2023 survey Rmk 3.2, corpus) | no Beltrami field in ℝ³ has finite `L²` norm, even variable-λ | the pure-Maxwell Beltrami electron is **impossible**; the Ouroboros confinement term is **forced** (§ 5b), the sharp form of the HydroBoros thesis |

The flexible sector is **constant-λ**: Enciso & Peralta-Salas (Ann. Math. 2012, corpus #6) realize
**arbitrary knots and links** as vortex tubes of Trkalian fields. So the Trkalian sector is rich
enough to seed every topology M7 needs, and the charge enters as the **measured** `∇·F` deviation
the Ouroboros coupling drives, not as an imposed exact ansatz. The diagnostic field
`λ_eff(x) = F·(∇×F)/\|F\|²` (local alignment eigenvalue) replaces the "find λ(x)" hunt.

**Pisello's precedent (via Álvaro García López, Models-of-Particles thread, 2026-06-30).** Daniel
Pisello's 1979 book *Gravitation, Electromagnetism and Quantised Charge* ([M7.0 corpus](tasks/m7_0_bootstrap.md) #16) builds a **toroidal
Hopf-knot** EM solution whose **charge is quantized ∝ the homotopy class of the knot** , not Rañada's
helicity-only zero-charge hopfion, but a genuinely **charged** knot from an electrovacuum Lagrangian
yielding **non-homogeneous equations á la Fleury-Dos Santos** (charge density ∝ the potential and its
derivatives, i.e. gauge-fixed, **A-primary**). This is the **closest published precedent to the M7
target** (a charged toroidal Beltrami knot) and sharpens three open questions: it **unifies** the
divergence charge with the topological / homotopy charge (Q3), shows a divergence-ful field **does**
hold a quantized-charge knot (Q5), and gives a concrete gauge-fixed charge-carrying construction (Q7).
Álvaro García López's own *Massive wave solutions to the Einstein-Maxwell equations* ([M7.0 corpus](tasks/m7_0_bootstrap.md) #17, Zitter
Institute) is the partner mechanism: the **electrovacuum as a charged nonlinear optical medium**, with
**Klein-Gordon from gauge-fixing** (M7.6) and an EM origin of mass.

The primary math is in **Pisello's 1977 paper** ([M7.0 corpus](tasks/m7_0_bootstrap.md) #18, the Lagrangian + **target-space S²**), now in
the corpus, alongside **Faber's geometric model** ([M7.0 corpus](tasks/m7_0_bootstrap.md) #19, arXiv:2201.13262), which builds the same
object on **target-space S³** (a soliton = electron with quantized charge, mass = field energy, spin),
plus **Faber & Golubich's precision lattice electron** ([M7.0 corpus](tasks/m7_0_bootstrap.md) #20, arXiv:2604.12021), the `α⁻¹ ≈ 137`
SU(2)-dipole result that is the M7.2 / M5.11 validation target. The **S² (Pisello) vs S³ (Faber)**
target-space choice is itself a live input to the substrate question **Q1** (what manifold the M7
field maps into).

---

## 3. Fleury's toroidal model , the analytic target to reproduce first

dos Santos & Fleury, *An electromagnetic model of the electron*, [arXiv:2510.22384](https://arxiv.org/abs/2510.22384).
The electron is a rotating EM wave confined to a torus:

```text
E = i E₀ H(r−r₀) e^{i(φ−ωt)} [ â_R + i (1 + R/R₀) â_φ ]
B = i B₀ H(r−r₀) e^{i(φ−ωt)} â_z ,        B₀ = E₀/c
phase ψ = (φ − ωt)  ,  one wavelength around the ring  ,  left-handed triad
```

Free parameters: amplitude `E₀`, major radius `R₀`, minor radius `r₀`, frequency `ω`. The model
identifies **charge with the field divergence** (`ρ = ∇·E`, non-zero in free space, the novelty vs
divergence-free Rañada/Irvine knots), and pins the four parameters with three QED constraints:

| Constraint | Equation (thin-torus limit) |
| --- | --- |
| RMS charge | `√2 π² ε₀ E₀ r₀² = e` |
| Magnetic moment | `√2 ε₀ π c E₀ R₀ r₀² (1 + r₀²/2R₀²) = μ_B (1 + α/2π)` |
| Spin | `(1/c) ε₀ E₀² π² R₀² r₀² (1 + r₀²/4R₀²) = ℏ/2` |
| Faraday (fixes ω) | `ω = 2c/R₀` (monochromatic), phase velocity `v_p = 2c` |

Solved numerical values (the **M7.2 reproduction targets**):

```text
E₀ ≈ 0.286 E_S  (Schwinger limit)        R₀ ≈ 1.573 r_c  ≈ 6.073e−13 m  (Compton scale)
r₀ ≈ 0.152 r_c  ≈ 5.854e−14 m            U  ≈ 0.795 m_e c² ≈ 0.406 MeV
ω  ≈ 0.636 ω_D  ≈ 9.86e20 rad/s          (ω_D = Dirac / zitterbewegung frequency)
```

Acknowledged limits (Fleury, § 5.2): the Heaviside mask is unphysical at the torus boundary
(suggests Bessel-function envelopes), and the energy falls short of `m_e c²`. **M7 replaces the
analytic mask with a relaxed lattice field**, which is precisely where the PDE treatment earns its
keep.

Fleury's ansatz is a **rotating wave** in Ceperley's sense (Fleury's ref [13], [M7.0 corpus](tasks/m7_0_bootstrap.md) #14): Ceperley's
circularly-polarized-EM rotating mode `E_r = E₀ e^{i(κz+φ−ωt)}` at `m=1` (Ceperley Eq 15) is literally this
torus phase, its **`J_m(κr)` Bessel envelope** is exactly the § 5.2 fix for the Heaviside mask, and
its angular-momentum law **`L_z = m(U/ω)`** (quantized when `U/ω = ℏ`) is the structural origin of the
spin `= ℏ/2` constraint. See [`../theory/ceperley_rotating_waves.md`](../theory/_CITATIONS.md) § 4b.

---

## 4. Ouroboros (M6) , the self-confinement the blend inherits

Werbos's two coupled vector fields `A_μ, J_μ` on Minkowski `(−,+,+,+)`
([`0d_canonical.md`](../../m6_ouroboros/research/archive/0d_canonical.md), Werbos's
*Evaluating Universe Model Alternatives v5* shared doc):

```text
L = −¼ F_μν F^μν − ¼ G_μν G^μν + m_J² A_μ J^μ − f(J_μ J^μ)
    F = dA  ,  G = dJ  ,  f(s) = (g/4) s²  (canonical electron choice, g = 1.0)
```

In the linear limit (`J→0`) the A-field reproduces Maxwell exactly; the nonlinear `m_J² A·J −
f(J·J)` coupling produces localized time-periodic solitons ("chaoitons"), spin lock-in `2L/Q = 2ω`,
charge from the Chern-Simons linking number `Q_CS = 1`. **All of this exists today only as a 1D
radial reduction** (the `(α,β)` cylindrical ODE for the charged sector; an `l=1` spherical BVP for
the neutral sector). M7 carries the same Lagrangian, `m_J`, `g`, and the `H/Q` benchmark into a full
3D lattice.

**The `H/Q` calibration ledger (pinned 2026-07-02, drives the M7.3 gate).** Three numbers circulate
and must not be conflated:

| Calibration | `g` | `H/Q` | Source |
| --- | --- | --- | --- |
| **M6 canonical (M7.3 primary target)** | 1.0 | **1.6890** | [`0d_canonical.md § 2.5`](../../m6_ouroboros/research/archive/0d_canonical.md) (repo-validated 1D BVP) |
| Werbos v5 canonical point | 1.0625 | 1.6969 | *Evaluating Universe Model Alternatives v5* (corpus #10, June 2026) |
| physical target | , | 1.6875 | v5 (`≈ 27/16`) |

M7.3 compares **like with like**: the 3D lattice vs the M6 1D BVP at the **same** `(g, ω, f)`; the
Werbos-v5 calibration is tracked secondary until the dictionary question (Q9) resolves. Werbos v5
also adds new structure M7 adopts as targets: the **`(Ω, G)` bifurcation islands** (electron
`Ω = 1.050` stable `k > 0`; muon `Ω = 0.914` **resonant/metastable** `k < 0`, a falsifiable
stable-vs-resonant distinction for M7.19), the **`β*` vacuum-stability threshold** (an M7.5 probe),
and the 319-family parameter scan (Zenodo 20866581). The `(Ω, G) ↔ (ω, g, m_J)` map is **Q9**.

**Structural note (an explicit M7.4 design item).** The two parents' electrons are *different field
configurations*: Fleury's torus lives in the `A`-sector with `B = ∇×A ≠ 0`; M6's chaoiton has
`A⃗ = 0` (so `B = 0`), the oscillation living in `A₀` and the azimuthal `J`. On M6's own ansatz the
coupling `m_J² A_μJ^μ` even vanishes pointwise (it acts through the EOM/derivative structure). The
blended HydroBoros electron presumably activates both sectors; which components carry the torus and
which the confinement is a design decision M7.4 must make explicitly, protected by the M7.3
verbatim-ODE gate (§ 5a).

---

## 5. The dynamics (what makes it rigorous, not analytic)

> **Reworked 2026-07-02** after the full theory review (Fleury FLDB, Werbos v5, M6 `0d_canonical.md`,
> Woltjer 1958, Enciso 2023, Kaiser 2000, Sutcliffe, M5.11). Two corrections vs the earlier draft:
> the frame is **time-harmonic** (there is no static electron in this model), and the stabilizer is
> **helicity + confinement**, not the drafted Faddeev-Niemi term (which is inert on Beltrami fields).

### 5a. The frame is time-harmonic (no static electron exists in this model)

Both parents are time-periodic: Fleury's torus is a rotating wave `e^{i(φ−ωt)}` (phase velocity `2c`)
and M6's chaoiton oscillates as the pair `A₀ = α(r)cos(ωt)`, `J = β(r)sin(ωt) φ̂`. A static relaxation
(the M5.11 frame, correct for Faber's static hedgehog) targets an object that does not exist here. M7
therefore builds the **time-harmonic reduced energy functional** from the start: every field component
carries one global frequency ω through a `(cos ωt, sin ωt)` component pair; the reduced functional
`E_ω[fields]` is the period-averaged energy; the soliton is the **minimizer of `E_ω` at fixed ω**
(then sweep ω for the spectrum, M6's lepton scan generalized to 3D). The de Broglie clock is thereby
IN the soliton from M7.1; M7.5 validates the reduction in real time rather than adding a clock later.

The exact reduction (component bookkeeping, constraint structure) is an M7.1 deliverable, gated at
M7.3 by the **verbatim-ODE check**: restricted to M6's ansatz, the 3D functional must reproduce the
M6 `(α,β)` ODE exactly, term by term (including the `2ωα` chiral cross-term). M6's
[`0d_canonical.md § 6`](../../m6_ouroboros/research/archive/0d_canonical.md) (ten sandbox versions of subtle
reduction failures: signs, regularity classes, Laplacians, measures) is the cautionary record. The
Lagrangian itself is M6's, unchanged (§ 4).

### 5b. Stabilization: helicity + confinement (the Faddeev-Niemi term is out)

The earlier draft added `κ ∫ \|F×(∇×F)\|²/\|F\|²` as the Derrick-evading stabilizer. The theory review
kills it: the term **vanishes identically on every Beltrami configuration** (`∇×F ∥ F` makes the cross
product zero, for constant AND variable λ, and rescaled Beltrami fields remain Beltrami), so along
exactly the family M7 targets it provides zero outward pressure; it is also singular wherever
`\|F\| → 0` (knot cores, far field). The stabilization structure that survives the corpus:

| Direction | Guard | Source |
| --- | --- | --- |
| collapse (`μ → 0`) | **helicity**: at fixed `H = ∫A·B`, magnetic energy obeys `E_B ~ H/μ` (Arnold's bound `E ≥ λ₁\|H\|`, `λ₁ ~ 1/size`), diverging under shrinkage | Arnold / Moffatt; Woltjer frame |
| expansion (`μ → ∞`) | **Ouroboros confinement** `m_J² A·J − f(J·J)`, scaling `~ μ³`; free Maxwell knots expand (Rañada's ball-lightning problem; M5.11's P2 measured smooth knots expanding) | Werbos (M6); M5.11 P2 |
| existence | **Nadirashvili's theorem**: no finite-`L²` Beltrami field exists in ℝ³, even variable-λ (Enciso 2023 survey, Remark 3.2, corpus) , a finite-energy pure-Maxwell Beltrami electron is **impossible**, so the confinement term is mathematically **forced**, not decorative | Enciso 2023 (corpus) |

The last row is the HydroBoros thesis with a theorem behind it: the Fleury-Werbos blend is not a
design choice, it is what existence requires. Caveat kept honest: the confinement term's net sign on
the soliton (it is not positive-definite) is established empirically by M6's stable 1D solutions;
the 3D balance is exactly what M7.3/M7.4 test. A 4th-order term remains available as an optional
M7.4 experiment (Q2), **off by default**.

### 5c. Solve methods

| Method | Use |
| --- | --- |
| **Fixed-ω minimization** of the harmonic functional (FIRE / L-BFGS) with **fixed-helicity relaxation** (project or penalize `dH` during descent) | the soliton (M7.2-M7.4, M7.6, Phase B-C) |
| **Reverse-mode Taichi AD** for `δE_ω/δ(fields)`, validated against a numpy finite-difference gradient to `~1e-12` **before trusting any run** | the gradient for relaxation |
| **Woltjer-Taylor known-answer**: fixed-helicity relaxation of `∫\|B\|²` on the periodic box must converge to the constant-λ curl eigenfield (ABC flow), `λ → 2π/L` and `E → λH` to grid accuracy | the M7.1 machinery gate (theorem-anchored) |
| **Minkowski leapfrog** (constrained integrator) | M7.5 real-time validation of the harmonic reduction + stability; M7.18 annihilation |

### 5d. Boundary conditions, gauge, units (M7.1 design decisions)

| Decision | Resolution |
| --- | --- |
| **BCs for charged sectors** | a net charge on a periodic lattice is Gauss-inconsistent (total boundary flux 0 ≠ `Q/ε₀`): charged configs use **vacuum-fixed boundaries** (Faber & Golubich's and Sutcliffe's practice); periodic boxes remain for neutral / net-zero configs (the M7.6 two-charge dipole is net-neutral) |
| **Gauge handling** | gauge orbits of `A` are flat directions of the minimizer, and `m_J² A·J` is gauge-sensitive off-shell; the concrete scheme (Coulomb gauge on `a⃗` + kept `a₀`, vs projection, vs penalty) is **Q8**, decided at M7.1 |
| **Units contract** | M6's natural units adopted wholesale (`c = 1`, electron `ω = 1`, `m_e` anchor; conversion table in [`0d_canonical.md § 5`](../../m6_ouroboros/research/archive/0d_canonical.md)); Fleury's targets are all dimensionless ratios (`E₀/E_S`, `R₀/r_c`, `U/m_ec²`, `ω/ω_D`), so M7.2 works in `r_c = 1` units and compares ratios directly |

---

## Cross-references

- Theory (our notes are in-repo; the source PDFs are local-only / gitignored, cited by public DOI/arXiv):
  [arXiv:2510.22384](https://arxiv.org/abs/2510.22384) (Fleury torus) ·
  Werbos *Evaluating Universe Model Alternatives v5* (shared doc, local only) ·
  [`../theory/sato_yamada_beltrami.md`](../theory/_CITATIONS.md) ([arXiv:1809.03136](https://arxiv.org/abs/1809.03136)) ·
  [`../theory/ceperley_rotating_waves.md`](../theory/_CITATIONS.md) (Ceperley rotating-wave equations; [AJP DOI 10.1119/1.17020](https://doi.org/10.1119/1.17020) / [IEEE DOI 10.1109/22.216476](https://doi.org/10.1109/22.216476)) ·
  [`../theory/feynman_maxwell_equations.md`](../theory/_CITATIONS.md) (Feynman Lectures II Ch 18, the Maxwell baseline for M7.2)
- Source corpus (the evolving M7 library, 64 docs): [`../theory/`](../theory/) , Marc's original 53 in [`electron_beltrami/`](../theory/electron_beltrami/) plus 11 later additions in the theory root, manifest [`../theory/_CITATIONS.md`](../theory/_CITATIONS.md) (Beltrami / force-free / knotted-EM / ball-lightning / LENR / Pisello-Faber canon, the M7.0 corpus)
- Rigor standard: [`../../m5_liquid_crystal/research/11a_vortex_loop.md`](../../m5_liquid_crystal/research/tasks/m5_11a_vortex_loop.md) (M5.11 vortex-loop)
- Wave-physics library to mine (§ 0, not parents): [`m1_granule_motion`](../../m1_granule_motion/) · [`m2_free_wave`](../../m2_free_wave/) · [`m3_wolff_lafreniere`](../../m3_wolff_lafreniere/) · [`m4_ewt`](../../m4_ewt/)
- Ouroboros canonical spec: [`../../m6_ouroboros/research/archive/0d_canonical.md`](../../m6_ouroboros/research/archive/0d_canonical.md) ·
  background [`../../m6_ouroboros/research/archive/0a_background.md`](../../m6_ouroboros/research/archive/0a_background.md)
- EM ≡ hydrodynamics bridge (already in M5): [`../../m5_liquid_crystal/research/1b_topological_defect.md`](../../m5_liquid_crystal/research/tasks/m5_1b_topological_defect.md) § "EM-hydrodynamics formal equivalence"
- Production-engine template: [`../../m5_liquid_crystal/`](../../m5_liquid_crystal/) (`medium.py`, `engine1_seeds.py` … `_launcher.py`)
- Comparison table (the M7 column goal): [`MODELS.md`](../../../../MODELS.md) · new-model flow [`ONBOARDING_MODELS.md`](../../../../ONBOARDING_MODELS.md), [`CONTRIBUTING.md`](../../../../CONTRIBUTING.md)
- Model briefing stub: [`../__M7_model_briefing.md`](../__M7_model_briefing.md)
- Question Tracker: [m7_question_tracker.md](m7_question_tracker.md)
- Roadmap: [m7_roadmap.md](m7_roadmap.md)
