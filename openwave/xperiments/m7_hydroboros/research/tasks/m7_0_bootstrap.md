# M7.0, bootstrap (collect the theory corpus + planning scaffolding)

> Task **M7.0** (M7 / HydroBoros). taskID = M7.N iteration. Status: **Done** (2026-07-02) · Roadmap: [`../m7_roadmap.md`](../m7_roadmap.md)

This doc is the task's full record: planning + findings + future planning + documentation. **M7.0 (bootstrap) is complete** (2026-07-02): the theory-source corpus is assembled (64 docs in `electron_beltrami/`, evolving from Marc's library) and the planning scaffolding stands up (background, question-tracker, roadmap, `research/` structure). The build proceeds at **M7.1**.

---

## Scope

The bootstrap collects two things:

1. **The theory corpus** (Task 0 sources below): the primary Fleury + Ouroboros + Beltrami papers and Marc's consolidated `electron_beltrami/` library, all landed in [`../../theory/`](../../theory/).
2. **The planning scaffolding**: the implementation background ([`../m7_background.md`](../m7_background.md): the two physics parents + M5 method + dynamics), the question tracker ([`../m7_question_tracker.md`](../m7_question_tracker.md)), the roadmap ([`../m7_roadmap.md`](../m7_roadmap.md)), and the `research/{scripts,data,plots,tasks}/` folder structure.

---

## Task 0, collect the theory sources (do this first)

| # | Source | ID / venue | Role in M7 | Status |
| --- | --- | --- | --- | --- |
| 1 | dos Santos & Fleury, *EM model of the electron* | arXiv:**2510.22384** | toroidal ansatz + the 3 QED constraints | ✅ in `theory/` |
| 2 | Fleury & Rousselle, *Critical Review of Zitterbewegung Electron Models* | Symmetry **17**, 360 (2025) | why point-particle zitter models fail; the design brief | ✅ in `electron_beltrami/` ([DOI 10.3390/sym17030360](https://doi.org/10.3390/sym17030360)) |
| 3 | Rañada, *Topological theory of the EM field* | Lett. Math. Phys. **18**, 97-106 (1989) | helicity = linking; the Hopfion construction | ✅ in `electron_beltrami/` (purchased PDF) |
| 4 | Rañada, *Knotted solutions of Maxwell in vacuum* | J. Phys. A **23**, L815 (1990) | knotted **divergence-free** EM (the foil Fleury breaks); **helicity = linking = Hopf index** `n = (1/𝒶)∫A·B d³r` (Eq 3), the velocity↔A / vorticity↔B hydrodynamic analogy | ✅ in `electron_beltrami/` |
| 5 | Kedia, Bialynicki-Birula, Peralta-Salas, Irvine, *Tying knots in light fields* | PRL **111**, 150404 (2013) | Bateman construction = the **seeder** for knotted initial data | ✅ in `electron_beltrami/` ([arXiv:1302.0342](https://arxiv.org/abs/1302.0342)) |
| 6 | Enciso & Peralta-Salas, *Knots and links in steady Euler flows* | Ann. Math. **175** (2012) | **Beltrami** fields realize any knot as a steady fluid vortex; the **Spanish Beltrami school** (ICMAT Madrid) Marc may enroll as collaborators | ✅ in `electron_beltrami/` ([arXiv:1003.3122](https://arxiv.org/abs/1003.3122)) |
| 7 | Faddeev & Niemi, *Knots and particles* | Nature **387** (1997); hep-th/9610193 | the 4th-order **Hopfion stabilizer** (Derrick-evading, parallels M5's Skyrme term) | ✅ in `electron_beltrami/` ([hep-th/9610193](https://arxiv.org/abs/hep-th/9610193)) |
| 8 | Sutcliffe, *Knots in the Skyrme-Faddeev model* | arXiv:0705.1468 | the canonical **numerical** Hopfion relaxation recipe | ✅ in `electron_beltrami/` ([arXiv:0705.1468](https://arxiv.org/abs/0705.1468)) |
| 9 | Werbos, Claude & DeepSeek, *Formal Statement of the Chaoiton Existence Theorem in Lean 4* | [Zenodo 20030162](https://zenodo.org/records/20030162) (2026) | the Chaoiton existence theorem, Lean-4 formalization; the Ouroboros Lagrangian + soliton baseline **narrative lives in #10** (the record itself is a `.lean.txt`, not fetched) | ⚠️ citation corrected; substance in #10 |
| 10 | Werbos, *Evaluating Universe Model Alternatives v5* | shared `.docx` (June 2026) | Ouroboros / TUFT params, `H/Q = 1.6969` | ✅ in `electron_beltrami/` (updated to newest v5, 2026-07) |
| 11 | **Sato & Yamada, *Local Representation and Construction of Beltrami Fields*** (Marc's Beltrami source) | arXiv:**1809.03136** (Physica D, 2019) | the construction recipe for the toroidal-Beltrami **seeder** (eikonal + equal-scale-factor rule); inhomogeneous + non-solenoidal `h` | ✅ [arXiv:1809.03136](https://arxiv.org/abs/1809.03136) (local PDF gitignored) + note [`sato_yamada_beltrami.md`](../../theory/_CITATIONS.md) |
| 12 | Duda, *Framework for liquid crystal based particle models* (superfluid / EM ≡ hydrodynamics mapping) | arXiv:2108.07896 | the EM ≡ hydrodynamics equivalence already cited in M5 | ✅ in `electron_beltrami/` ([arXiv:2108.07896](https://arxiv.org/abs/2108.07896)) |
| 13 | Marc's **2 non-constant-λ (variable-λ) Beltrami papers** = **Kaiser** (*Force-free fields, nonconstant α*, CMP 211, 2000) + **Kravchenko** (*Beltrami fields, nonconstant proportionality factor*, J.Phys.A 36, 2003) | in the corpus (§ 3b) | the **variable-λ** case where the charge lives (`∇·w ≠ 0`); the target to generalize S&Y onto | ✅ both in `electron_beltrami/` (+ the Kravchenko & Oviedo *...on the plane* companion); the core of Q7 / M7.4 |
| 14 | Ceperley, *Rotating Waves* (Fleury's ref [13]) | Am. J. Phys. **60**, 938 (1992); DOI 10.1119/1.17020 | the **phase-vortex / rotating-wave** formalism; **circularly-polarized EM at `m=1` = Fleury's torus** (Eq 15), spin `L_z=mU/ω` + QM bridge (`U/ω=ℏ`), spherical + radiating forms, **Bessel envelope** (= Fleury's § 5.2 fix) | ✅ [DOI 10.1119/1.17020](https://doi.org/10.1119/1.17020) (paywalled; local PDF gitignored) + note [`ceperley_rotating_waves.md`](../../theory/_CITATIONS.md) |
| 15 | Velazco & Ceperley, *Rotating Wave Fields for Microwave Applications* | IEEE Trans. MTT **41**(2), 330 (1993) | applications companion to #14: the **full cylindrical-cavity `E,H` field set** + the `ω_rot=ω/m` derivation | ✅ [DOI 10.1109/22.216476](https://doi.org/10.1109/22.216476) (paywalled; local PDF gitignored) |
| 16 | Pisello, *Gravitation, Electromagnetism and Quantised Charge: The Einstein Insight* | book, Ann Arbor Science (1979), ISBN 9780250402861 | a **toroidal Hopf-knot** EM solution with **charge quantized ∝ the knot's homotopy class** (vs Rañada's helicity / zero charge); electrovacuum Lagrangian → **non-homogeneous (charged) equations á la Fleury-Dos Santos**, the **closest published precedent to the M7 charged toroidal knot** | ✅ in `electron_beltrami/` (scanned PDF acquired; image-only, no OCR layer) |
| 17 | García López, *Massive wave solutions to the Einstein-Maxwell equations* | Álvaro García López (URJC); [Preprints.org 202504.0927 v2](https://www.preprints.org/manuscript/202504.0927/v2) (2025) | the **electrovacuum as a charged nonlinear optical medium**; mass + current from self-coupling + gravitational back-reaction; **Klein-Gordon from gauge-fixing**, the partner charge + KG mechanism (M7.6) | ✅ in `electron_beltrami/` |
| 18 | Pisello, *Nonlinear Classical Theory of Electromagnetism* | Int. J. Theor. Phys. **16**(11), 863-866 (1977) | the **original "donut electron"**, the primary paper behind #16, with the actual Lagrangian + **target space S²** | ✅ in `electron_beltrami/` |
| 19 | Faber, *A geometric model in 3+1D space-time for electrodynamic phenomena* | arXiv:**2201.13262**; Universe **8**(2):73 (2022), [DOI 10.3390/universe8020073](https://doi.org/10.3390/universe8020073) | a spherically-symmetric **topological soliton = electron**: quantized charge, **mass = field energy**, spin; **target space S³**, the Faber model we already reproduce in M5.11 | ✅ in `electron_beltrami/` (open access) |
| 20 | Faber & Golubich, *High-precision lattice determination of the SU(2) solitonic-dipole potential vs QED* | arXiv:**2604.12021** (2026) | the **precision Faber electron**: `E(d)=2m_ec²−α_sol ℏc/d`, `α_sol⁻¹=137.1`, QED running of `α(d)` via CG lattice minimization, the **M7.2 / M5.11 validation target** | ✅ in `electron_beltrami/` (also in M5 theory) |
| 21 | Faber, *A Model for Topological Fermions* | [hep-th/9910221](https://arxiv.org/abs/hep-th/9910221) (v4, 2000) | the **foundational S³ topological-soliton model** (sine-Gordon 1+1D → 3+1D); the **parent** of the Faber electron in #19 / #20, the lineage we reproduce in M5.11 | ✅ in `electron_beltrami/` |
| 22 | Duda, *Hydrodynamical Analogues of Quantum Phenomena* (walking-droplet / Couder-Fort notes) | Duda slide-notes (2026) | **hydrodynamic QM analogs** (Couder-Fort walkers, de Broglie-Bohm pilot wave, tunneling / orbit quantization); the **fluid-intuition source for HydroBoros** (the EM ≡ hydrodynamics bridge, ties to #12) | ✅ in `electron_beltrami/` |

Note on #11: Marc shared the Sato-Yamada writeup directly (the Gemini share link is Google-auth gated and cannot be fetched headless). The paper + Marc's summary are now in `theory/`; the summary note ties each result to a specific M7 piece.

### 3b. Marc's electron-Beltrami corpus (consolidated + complete, 2026-06-30)

The author sent a **complete consolidated library**: **52 PDFs + 1 docx** in [`../../theory/electron_beltrami/`](../../theory/electron_beltrami/) (year-prefixed `YEAR - Author - Title.pdf`). This is the **source-of-record**, a superset of the original 20-source NotebookLM list ("HydroBoros: the Beltrami-OpenWave implementation") plus ~28 further papers. **All originals are now present**, including the ones we could not fetch openly (Kaiser, Dombre, Kravchenko-2003, Rañada-Maxwell). The earlier lossy Gemini-**ingested text** folder was **removed** (superseded by the full PDFs); the two originals we had that Marc's set lacked (Fleury *Zilch-Zitter*, Kravchenko & Oviedo *...on the plane* 2008) were folded into `electron_beltrami/`. **The corpus is an evolving M7 library: it started from Marc's consolidated 53 (the source-of-record) and has since grown to 64** (2026-07) with **11 separately-sourced additions**: the Task-0 supplements #2 / #5 / #6 / #7 / #8 / #12 (fetched from arXiv + MDPI), plus Pisello-1979 #16 (scanned), Rañada-1989 #3 (purchased), García López #17, Faber *A Model for Topological Fermions* #21, and Duda's *Hydrodynamical Analogues* #22. Full manifest + per-file provenance: [`../../theory/_CITATIONS.md`](../../theory/_CITATIONS.md).

Highlights by M7 relevance (the full 48 are in the folder):

| Theme | Standout papers in `electron_beltrami/` |
| --- | --- |
| **Electron model (toroidal / clock)** | dos Santos *EM Model of the Electron* (FLDB Main + Appendix) + *Poloidal-Toroidal System v8e*; Fleury *Zilch-Zitter Electron*; Ceperley *Rotating Waves*; Hestenes *Zitterbewegung Structure*; **Catillon** *de Broglie particle internal clock* (the measured ZBW); **Wan** *Toroidal Vortices of Light*; Kovacs *What is Inside an Electron* |
| **Force-free / Beltrami / Trkalian** | Sato *Local Rep. of Beltrami Fields*; Woltjer; Marsh *Force-Free Magnetic Fields*; **Kaiser** + **Kravchenko** (the variable-λ pair, both now present) + Kravchenko & Oviedo *...on the plane*; Reed *Beltrami-Trkalian*; Dombre *ABC Flows*; **Enciso** *Beltrami Fields & Knotted Vortex Structures*; Mochizuki *Zero-Poynting EH Beltrami*; Takahashi *Beltrami vortex solutions* |
| **Knotted EM / topological / Clebsch** | Rañada *Knotted Maxwell*; Trueba *Topological EM*; Irvine *Linked & Knotted Beams*; **Arrayás** *Fibration by field lines*; Nastase *Fluid-EM Helicities*; Migdal *Clebsch Confinement*; **Wiegmann** *Chiral Anomaly in Euler Fluid & Beltrami Flow*; Yoshida *Nambu / Clebsch* |
| **Hopfions / chiral solitons** | **Kent** *Hopfions in Magnetic Multilayers* (experimental); Hall *Fusion & Fission of Particle-Like Chiral Structures* |
| **Ball lightning / LENR (applications)** | Rañada *Ball Lightning as Force-Free Knot*; Donoso *Riddle of Ball Lightning*; Biberian *Pd transmutation*; Lewis; Fleury *SuperColdFusion LENR*; Fomitchev *Neutron emission* |
| **Foundations / context** | Cheng *Field & Wave EM*; Barrett *Advanced EM*; Puthoff *Vacuum ZPE*; Celani *Maxwell & Occam*; Hillion; Rapoport *Cartan-Weyl Dirac*; Barbarosie *Divergence-free fields*; Brady *Emergent QM*; Morguer *QED test* |

The corpus confirms and sharpens Marc's framing: **toroidal-EM electron** + **force-free / variable-λ Beltrami** (where the charge lives) + **knotted-EM / Clebsch** topology, with a **ball-lightning-as-force-free-knot** thread and an **LENR** applications thread. New high-value finds for the build: **Catillon** (the measured de Broglie clock, ties to Fleury's `ω`), **Wan** *Toroidal Vortices of Light* (a toroidal EM structure directly adjacent to the M7 soliton), **Enciso** (the Spanish-school Beltrami-knot existence results, Q4 / Q6 collaborators), and **Wiegmann** (the chiral-anomaly bridge between Euler fluid and Beltrami flow). Provenance caveat: the **Reed** piece is a viXra non-reviewed essay, read critically (seeds-only discipline, [`../m7_question_tracker.md`](../m7_question_tracker.md)).

---

Cross-refs: roadmap [`../m7_roadmap.md`](../m7_roadmap.md) (task M7.0) · full background [`../m7_background.md`](../m7_background.md) (§ 3 Fleury, § 4 Ouroboros, § 5 dynamics) · corpus manifest [`../../theory/_CITATIONS.md`](../../theory/_CITATIONS.md) · corpus folder [`../../theory/electron_beltrami/`](../../theory/electron_beltrami/).
