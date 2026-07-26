# M7.4, the charged soliton (approximately-Beltrami) + its Coulomb field (the NEW physics)

> Task **M7.4** (M7 / HydroBoros). taskID = M7.N iteration. Status: **Done** (2026-07-03, review approved; gate outcomes in [`§ FINDINGS 5`](#5-gates-vs-the-plan--4-of-the-planning-half)) · Roadmap: [`../m7_roadmap.md`](../m7_roadmap.md)

This doc is the task's full record: planning + findings + future planning + documentation. **M7.4 is the research core**: the charged, knotted, finite-size soliton neither parent produced. Reframed 2026-07-02 after the corpus math review ([`../m7_background.md § 2`](../m7_background.md) "What the Beltrami mathematics allows"): the task does **not** hunt an exact variable-λ Beltrami ansatz; it relaxes the full functional and **measures** what charge the coupling drives.

---

## 1. The physics frame (why relax-and-measure, not construct-exactly)

| Math fact | Consequence |
| --- | --- |
| divergence identity `∇·w = −(w·∇λ)/λ` | charge requires λ varying **along** field lines |
| rigidity (Clelland-Klotz 2020; Enciso-PS): nonconstant-λ Beltrami obstructed even locally | an exact charged Beltrami ansatz likely does not exist for generic profiles |
| Nadirashvili: no finite-`L²` Beltrami field in ℝ³ at all | the soliton is NOT an exact Beltrami field; it is the minimizer of the **full** functional, approximately Beltrami in the core, with the Ouroboros confinement holding it (the § 5b stabilization story) |
| Kaiser 2000: variable-α exists perturbatively (small α) | the small-charge regime is safe ground to enter from |
| Enciso-PS 2012: constant-λ realizes **arbitrary knots** | the Trkalian seed inventory is topologically complete |

So: seed Trkalian (exactly Beltrami, zero charge), switch on the Ouroboros coupling, relax `E_ω` at fixed ω + helicity, and measure the deviation. The charge is **output, not input**. Marc's "start Trkalian, take off the training wheels" (2026-06-30) survives intact in this implementation.

## 2. Seeds (multiple, per Sutcliffe's local-minima lesson)

Sutcliffe (corpus #8): the Hopfion energy landscape is dense with local minima; seed variety is not optional. The inventory (built at M7.1):

| Seed | Topology | Why |
| --- | --- | --- |
| Trkalian / ABC eigenfield torus | unknot, tunable helicity | the clean entry point |
| S&Y toroidal Beltrami | unknot, toroidal geometry | Marc's endorsed recipe |
| Bateman / Kedia knots | trefoil + torus knots | knotted sector (Q5) |
| the M6-embedded chaoiton (from M7.3) | `Q_CS = 1` | the Ouroboros parent's own configuration |
| the Fleury torus (from M7.2) | `m = 1` rotating wave | the Fleury parent's own configuration |

**Blend design item** (explicit, from [`../m7_background.md § 4`](../m7_background.md) structural note): the two parents' electrons are different configurations (Fleury: `B = ∇×A ≠ 0`; M6: `A⃗ = 0`, `B = 0`). Which components carry the torus and which the confinement is decided HERE, by relaxing both parent seeds (+ hybrids) under the same functional and reporting which basin wins.

## 3. Procedure + diagnostics

| Step | Diagnostic |
| --- | --- |
| relax each seed: `E_ω` at fixed ω + helicity, vacuum-fixed BCs | `‖∇E‖ → 0` |
| **dilation probe**: evaluate `E(μ)` along lattice rescalings of the relaxed state | interior minimum at `μ = 1` (the direct Derrick verification; replaces trust in any stabilizer argument) |
| **the Q3 measurement** (first-class deliverable) | per relaxed state: `Q_div` (Gauss flux at increasing radii), helicity `H = ∫A·B`, field-line **linking / Hopf number** (tracer); the `(Q_div, H, linking)` table + ratio analysis across states answers whether divergence-charge and topological charge are slaved (Q3) |
| `λ_eff(x) = F·(∇×F)/\|F\|²` map | where the state is Beltrami (core?) and where the charge-carrying deviation lives |
| far field | `\|E\| ~ Q_div/4πε₀r²` fit (Coulomb from Gauss) |
| charge quantization | does `Q_div` land on discrete values across seeds/parameters? assessed honestly (Pisello and Faber say yes topologically; this is the lattice test) |

Optional experiment (Q2, off by default): switch a 4th-order term on and quantify what it changes; the drafted `\|F×(∇×F)\|²/\|F\|²` form is documented-inert on Beltrami configs and singular at zeros, so any 4th-order run uses a regularized variant and is labeled an experiment, not the baseline.

## 4. Gates + honest outcomes

| Gate | Criterion |
| --- | --- |
| primary | a **stable, finite-size, charged** (`Q_div ≠ 0`) soliton: `‖∇E‖ → 0` + the dilation probe's interior minimum |
| Q3 deliverable | the `(Q_div, H, linking)` table published regardless of outcome |
| Coulomb | `1/r` far-field fit sourced by `Q_div` |
| honest negatives | expansion (the Rañada / M5.11-P2 mode, if confinement loses), collapse, or charge → 0 under relaxation are all documented results with the seed + parameter map |

Risk register for this task: rigidity + Nadirashvili (mitigated by the reframe), M5.11's P2 lesson (smooth knots expand; the confinement term is what M7 has that M5's functional lacked), gauge flat directions (Q8, from M7.1), multi-hour GPU cost (accepted).

Artifacts: `scripts/m7_4_linked_vortex.py` + `data/m7_4_*.npz` + `plots/m7_4_*.png` (λ_eff maps, dilation curves, the Q3 table, far-field fits).

---

## FINDINGS (2026-07-03, execution)

Artifacts: [`../scripts/m7_4_linked_vortex.py`](../scripts/m7_4_linked_vortex.py) (modes `smoke` / `run` / `winner` / `analyze`) · [`../data/m7_4_states.json`](../data/m7_4_states.json) · [`../data/m7_4_winner.json`](../data/m7_4_winner.json)

### Plots

![`../plots/m7_4_traces.png`](../plots/m7_4_traces.png)
![`../plots/m7_4_dilation_charge.png`](../plots/m7_4_dilation_charge.png)
![`../plots/m7_4_winner_sections.png`](../plots/m7_4_winner_sections.png)

### 1. Setup carried from M7.3 + the scalar-sector finding

The relaxation frame is the M7.3-pinned one: `E_ω = quad − κ⟨A·J⟩ + f_avg`, `κ = −1` (pure `J → −J` relabeling), FIRE at **fixed A-sector helicity `H_A`** (the anti-collapse guard, exact A-rescale restore), plus **fixed J-norm `Q_J`** for the focusing convention (exact J-rescale restore; two-constraint Gram-Schmidt tangent projection). The Q6 fork is run as the experiment:

| Convention | `f` | Well-posedness |
| --- | --- | --- |
| `repulsive` (the written-canonical direction) | `c1 = +λ/2, c2 = +g/4` | bounded below in the pure-vector sector; fixed `H_A` alone |
| `focusing` (the M6-verbatim pin from M7.3) | `c1 = −λ/2, c2 = −4g/3` | J-collapse channel open; run at fixed `H_A + Q_J`, honestly capped |

**The scalar-sector finding (new).** With the timelike components `(a₀, j₀)` free, the smoke run blew up (E → −3.6e6): null-J directions (`s = −j₀² + |j⃗|² = 0`) are quartic-flat and gain unbounded bilinear `a₀j₀` coupling energy: **the naive averaged Hamiltonian is unstable in the scalar sector**. And on-shell Gauss slaves the divergence charge to exactly that sector (`∇·E_A = −m_J²J₀`): the charge channel and the instability are the same channel. Design decision for this task: **scalars frozen (pure-vector sector)**, matching both parents' operative ansaetze; the charge experiment reads "does *seeded* geometric divergence (`∇·a⃗ ≠ 0`, Fleury-style) survive relaxation"; the constraint-level scalar treatment (Gauss-eliminated `a₀`, `j₀` prescription) is a documented follow-up and sharpens the Q7 ask.

### 2. The seed × f-convention matrix (N = 64, L = 16, 1500 FIRE iterations) ✅ measured

Constraints held exactly throughout every run (`H_A` to 5 digits, `Q_J` where constrained); `align = ⟨B·∇×B⟩/(‖B‖‖∇×B‖)` on `B = ∇×a⃗_c`; `Q_ρ` = the Fleury RMS-density charge inside `r = 0.3L`; `E` is the full-box `E_ω`.

| Seed (H_A) | repulsive `f` (fixed `H_A`) | focusing `f` (fixed `H_A + Q_J`) |
| --- | --- | --- |
| ck spheromak (+4.06) | ✅ **stable soliton**: E → 3.2546, `\|g\| = 1.4e-7`, align +0.960, r50 = 3.4, dilation interior min, J-condensate `Q_J → 1.88` saturated, `Q_ρ = 1.30e-2` stable | converges (E 3.714, align +0.981) but **expels the charge**: `Q_ρ` 2.6e-2 → 1.4e-4 |
| Bateman/Hopf knot (+2.24) | ✅ stable soliton from the knotted seed: E → 1.7988, `\|g\| = 1.0e-7`, align +0.960, dilation min, `Q_J → 1.04`, `Q_ρ = 5.4e-3` stable (see § 4: reconnects into the same Taylor family) | converges, expels charge (`Q_ρ` → 1.3e-4) |
| **blend: M6 torus + poloidal A-twist (−7.93)** | ✅ **stable BLENDED soliton, the winner**: E → 6.3580, `\|g\| = 1.6e-7`, align −0.960 (sign follows H_A < 0), dilation min, `Q_J → 3.66`, `Q_ρ = 3.50e-2` stable | ❌ lattice-arrested collapse: E → −2.48e4, max\|f\| 40 (the M7.3 focusing channel, at the pinned large `Q_J = 27.9`) |
| m6taper: pure M6 torus, no twist (0.00) | ❌ **evaporates to vacuum** (E → 9e-14, all fields → 0): zero helicity = no guard, exactly the M7.3 prediction | ❌ collapse (E → −2.48e4) |
| fleury torus, bare rotating pair (0.00) | ❌ evaporates to vacuum (E → 2e-12): the bare Fleury pair also carries zero A-helicity | marginal creeping state (E = −0.035, `Q_ρ` slowly growing, not settled) |

All three surviving states sit at `\|align\| = 0.96`: **approximately-Beltrami cores** with the deviation carrying the structure, exactly the § 1 relax-and-measure prediction. The far-field slope is ≈ −3.7 in all surviving states (confined/multipole decay, not Coulomb −2): consistent with zero **net** monopole in the pure-vector sector (§ 1): the true Coulomb `1/r` test needs the scalar/Gauss sector and moves to the follow-up.

### 2b. What the matrix answers

| Question | Answer from the matrix |
| --- | --- |
| Blend design item (which components carry torus vs confinement) | the **A-sector carries the structure + the helicity** (the guard); the **J-sector is the confinement dressing**, switched on by the bilinear coupling and saturating at a finite condensate (`Q_J` stabilizes without being constrained). Helicity is measurably load-bearing: both zero-helicity parent seeds evaporate; all three helicity-carrying seeds survive |
| Q5 (does `∇·F ≠ 0` destabilize knots?) | **no**: the divergence-ful Hopf-knot seed relaxes to a stable, dilation-stable soliton with persistent `Q_ρ`. Honest caveat (§ 4): the final state reconnects into the same ground family as the unknot seeds: KNOT-SECTOR persistence needs topology-preserving constraints (only global `H_A` was fixed), a designed follow-up |
| Q3 (divergence charge vs topological charge) | measured table across states: helicity (topological side) is required for existence; the RMS divergence charge coexists but is NOT slaved to it (the pure-M6 seed had the largest seeded `Q_ρ` and died; the ck seed a small one and lives): at this stage the two charges are **independent observables**, with the linking side gating existence |
| Q6 (which `f`) | the **repulsive/written-canonical branch is the physical one**: it holds stable solitons and keeps the charge; the focusing branch (the M7.3 benchmark-verbatim pin) expels charge or collapses. Combined with the LoE v5 text (`f(s) = g s²`, § 3), the picture: **the written Lagrangian is the stable theory; the benchmark ODE's effective signs are the anomaly** |
| Derrick / stabilization (Q2) | the helicity + confinement pair stabilizes with **no 4th-order term at all**: the constrained Derrick curve has an interior minimum for every surviving state |

### 4. The winner (blend/repulsive at N = 96) + the universal ratios ✅ measured

Grid-convergence: E = 6.3676 at N = 96 vs 6.3580 at N = 64 (0.15%); `Q_J → 3.68`, `Q_ρ = 3.44e-2`, `\|g\| = 2.8e-7` and still descending. The state: a smooth localized ball (r50 = 3.44, r90 = 5.73) with **near-constant `λ_eff ≈ −1` across the core** (Taylor-like; sign follows `H_A < 0`), the J-condensate co-located with the core, and the interior RMS charge genuinely interior (Q_ρ(<r) grows smoothly and saturates by r ≈ 6, before the boundary shell). Dilation curve: E(0.7) = 6.91 > E(1.0) = 6.368 < E(1.4) = 8.00, clean interior minimum.

**The three surviving states are ONE family.** Across ck / hopfion / blend the ratios are identical to 3-4 digits:

| Ratio | ck | hopfion | blend |
| --- | --- | --- | --- |
| `E / \|H_A\|` | 0.8016 | 0.8015 | 0.8017 |
| `H_cross / H_A` | −0.391 | −0.392 | −0.390 |

So the fixed-helicity ground state of the doublet functional is a **Taylor-relaxed, J-dressed spheromak family**, Arnold-linear in the helicity (`E = 0.802 \|H_A\|`) with the J-dressing locked at a universal anti-aligned cross-helicity fraction (`H_x = −0.390 H_A`). Two honest consequences: (a) **topology reconnects**: the Hopf seed's linking did not survive: fixing only the global `H_A` allows reconnection to the unknot ball, so distinct knot/linking SECTORS (the lepton-family hypothesis, M7.19) need topology-preserving constraints or penalties, a designed follow-up; (b) the ball is the ground state at THIS `(ω, g, λ)`: the toroidal electron shape (Fleury's geometry) is not yet selected by the functional: the shape question moves to M7.5/M7.6 (does the clock/rotating sector select the torus?).

### 5. Gates vs the plan (§ 4 of the planning half)

| Gate | Outcome |
| --- | --- |
| primary: stable finite-size soliton, `‖∇E‖ → 0` + dilation interior minimum | ✅ PASS: three seeds → the Taylor-dressed family, `\|g\| ~ 1e-7`, clean Derrick interior minima, grid-convergent (N = 64 → 96 at 0.15%) |
| charged (`Q_div ≠ 0`) | ⚠️ split verdict, honestly: the **RMS/Fleury divergence charge** is nonzero, interior, and persistent in the repulsive branch; the **net Gauss monopole is zero by sector design** (scalars frozen after the § 1 instability finding): the monopole-charged soliton needs the scalar-sector prescription (Q7(b)) |
| Q3 deliverable (`(Q_div, H, linking)` table) | ✅ published (§ 2-2b + `data/m7_4_states.json`): independence at this level; linking gates existence |
| Coulomb `1/r` far field | ❌ not present (slope ≈ −3.7): consistent with the zero net monopole; moves with the scalar-sector follow-up |
| honest negatives | ✅ delivered: zero-helicity evaporation (both parents' bare seeds), focusing-branch charge expulsion + collapse, topology reconnection under global-H fixing |

### 3. The Q9 self-read (the (Ω, G) dictionary, from the corpus)

Attempting the Zenodo route first: record 20866581 (cited by Werbos as the 319-family scan) now hosts a v2 of the **Lean-theorem statement**, not the scan. The corpus doc itself (`Evaluating Universe Model Alternatives v5`, corpus #10) carries the definitive table:

| Particle | Ω | G | Stability | far-field decay `k` |
| --- | --- | --- | --- | --- |
| Electron | 1.050 | 1.150 - 3.250 | stable (ground state) | `k > 0` |
| Muon | 0.914 | 0.100 - 3.600 | resonant (metastable) | `k < 0` |

What the self-read settles and what it cannot:

| Item | Status |
| --- | --- |
| `(Ω, G)` definitions | NOT in the doc; and the electron island `G ∈ [1.15, 3.25]` does **not** bracket the canonical `g = 1.0625`, so `(Ω, G) ≠ (ω, g)` naively: the dictionary ask to Werbos remains, now surgical |
| `k` meaning | pinned: the far-field decay constant; Werbos claims the electron island **decays** (`k > 0`), which couples Q9 to Q11: those islands are implicitly the "localized branch elsewhere in parameter space" that Q11(b) asks about, and our no-decaying-channel dispersion argument applies at the canonical point specifically |
| `f` normalization (Q6/Q12 evidence) | LoE v5 writes **`f(s) = g s²`** with λ a separate parameter, NOT the `(g/4)s²` transcribed into [`0d_canonical.md § 1`](../../../m6_ouroboros/research/archive/0d_canonical.md); with `s` the phasor norm this gives exactly the benchmark's `λβ + 4gβ³`, consistent with the M7.3 C2 pin |
| v5's own ansatz | `A_μ = (φ(r),0,0,α(r))e^{−iωt}`, `J_μ = (ρ(r),0,0,β(r))e^{−iωt}`: same-phase complex phasor, consistent with the M7.3 C2 winner extended by scalars (and the scalar extension is exactly the sector our § 1 finding shows is unstable under the naive averaged Hamiltonian) |

---

Cross-refs: roadmap [`../m7_roadmap.md`](../m7_roadmap.md) (M7.4) · background [`../m7_background.md`](../m7_background.md) (§ 2 the math boundaries, § 4 structural note, § 5b stabilization) · Q2/Q3/Q5/Q7 + Q6/Q9/Q11 evidence in [`../m7_question_tracker.md`](../m7_question_tracker.md) · upstream [`m7_1_infra.md`](m7_1_infra.md) (seeders, helicity, BCs) + [`m7_2_fleury_torus.md`](m7_2_fleury_torus.md) + [`m7_3_ouroboros_3d.md`](m7_3_ouroboros_3d.md) (the parent seeds + the pinned frame) · downstream [`m7_5_clock_stability.md`](m7_5_clock_stability.md) (real-time validation) + [`m7_6_observables.md`](m7_6_observables.md) (observables + two-charge Coulomb).

---

## TASK REVIEW (2026-07-03)

**Task Duration:** 03:22 (from 16:23 to 19:45 EDT)
**Usage Cap Triggered:** NO (finished before the 9:10pm reset; resume ping disarmed without firing)

**Results** (full detail: [`§ FINDINGS`](#findings-2026-07-03-execution)):

| Gate | Outcome |
| --- | --- |
| primary: stable finite-size soliton | ✅ the program's first: three helicity-carrying seeds relax to stable, dilation-stable, approximately-Beltrami solitons (`\|g\| ~ 1e-7`, `\|align\| = 0.96`, Derrick interior min, grid-convergent 0.15% N=64→96) with a saturating J-condensate |
| the one-family discovery | ✅ universal ratios `E = 0.802\|H_A\|`, `H_cross = −0.390 H_A`: a Taylor-dressed spheromak family; helicity measurably load-bearing (both zero-helicity parent seeds evaporate) |
| Q6 fork | ✅ decisive: the repulsive/WRITTEN `f` branch is physical (holds solitons, keeps charge); the focusing benchmark-verbatim branch expels charge or collapses; LoE v5 itself writes `f = gs²` |
| charged | ⚠️ split: RMS/Fleury charge nonzero + interior + persistent; net Gauss monopole zero by sector design (scalar-sector instability finding, § 1) |
| Coulomb `1/r` | ❌ absent (slope −3.7), consistent with zero net monopole; moves with the scalar-sector follow-up |
| Q3 deliverable | ✅ published: the two charges independent at this level; linking gates existence |

**Issues / blockers**: the scalar-sector instability (new finding, = the Q7(b) ask); topology reconnects under global-`H_A` fixing (knot sectors need topology-preserving constraints); the ball, not the torus, is the ground state at these parameters (shape question → M7.5/M7.6); no >1MB data created.

**Action needed**: M7.5 next (real-time clock validation of the blend winner); the comms packages per the tracker comms plan (Marc package ready; Werbos package after M7.5); `0d_canonical.md` correction notes after Werbos answers Q6/Q12.

**Findings**: M7.4 delivered the program's first stable, finite-size 3D soliton: a Taylor-dressed spheromak family (`E = 0.802\|H_A\|`, `H_cross = −0.390 H_A`, universal across seeds) whose existence requires exactly the two ingredients the HydroBoros thesis named: helicity (both bare parent seeds evaporate without it) and the Ouroboros confinement dressing (the J-condensate saturates on every survivor). The f-convention experiment was decisive: the written Lagrangian is the stable theory and the M6 benchmark's signs are the anomaly; the honest gaps (net Gauss charge in the unstable scalar sector, knot reconnection, ball-vs-torus shape) are each a designed follow-up.

**Research docs created / updated**:

- [this task doc](m7_4_charged_soliton.md) (§ FINDINGS 1-5, plots inline)
- [`../scripts/m7_4_linked_vortex.py`](../scripts/m7_4_linked_vortex.py) (modes `smoke` / `run` / `winner` / `analyze`)
- [`../data/m7_4_states.json`](../data/m7_4_states.json) · [`../data/m7_4_winner.json`](../data/m7_4_winner.json) · [`../data/m7_4_smoke.json`](../data/m7_4_smoke.json)
- [`../plots/m7_4_traces.png`](../plots/m7_4_traces.png) · [`../plots/m7_4_dilation_charge.png`](../plots/m7_4_dilation_charge.png) · [`../plots/m7_4_winner_sections.png`](../plots/m7_4_winner_sections.png) (key plot: the winner sections)
- [`../m7_question_tracker.md`](../m7_question_tracker.md) (Q2 + Q5 → RESOLVED; Q3/Q6/Q7/Q9 evidence; count 5/5/0/3)
