# M9 Emergent Gravity / NSM: Model Briefing (draft)

> **What M9 brings.** Gravity as entanglement bookkeeping, assembled as the
> Standard Model minimally coupled to Einstein-Cartan gravity. Matter is
> installed, not emerged. The one interaction beyond SM + Einstein gravity is
> the Hehl-Datta axial-axial contact, coefficient fixed by \(G\). This is a
> gravity-certification column, closer to M8 than to M5. The model ID is
> proposed; a maintainer assigns the official number.

## Identity

| Field | Value |
| --- | --- |
| Model ID | M9 (proposed) |
| Name | Emergent Gravity / New Standard Model (NSM) |
| Author | Dr. Robert W. McGwier, PhD, CTO, Cohere Technology Group (sole author) |
| Author contact | GitHub [@n4hy](https://github.com/n4hy), for author-gated questions (definitions, intent, what the model does and does not claim); routing convention in [`dev_docs/CROSS_MODEL_TESTING.md`](../../../dev_docs/CROSS_MODEL_TESTING.md) § 6 |
| Lineage | Faulkner-Guica-Hartman-Myers-Van Raamsdonk 2014; Casini-Huerta-Myers; Einstein-Cartan-Sciama-Kibble; Hehl-Datta 1971; Jacobson entanglement equilibrium |
| Primary sources | Author repo [github.com/n4hy/New_Model_Emergent_Gravity](https://github.com/n4hy/New_Model_Emergent_Gravity) (CC-BY-4.0 PDFs under `research/`); registry in [`theory/_CITATIONS.md`](theory/_CITATIONS.md). Specification of the first task: Paper III action (2) |
| In-repo | Headless only. No launcher. First task: [`research/scripts/hehl_datta.py`](research/scripts/hehl_datta.py) |
| Application | [discussion #442](https://github.com/openwave-labs/openwave/discussions/442); [`research/APPLICATION.md`](research/APPLICATION.md); first PR [openwave-labs/openwave#441](https://github.com/openwave-labs/openwave/pull/441) |

## Model Profile (what it brings, short form)

| Attribute | M9 |
| --- | --- |
| Substrate | Coframe \(e^a\) and independent Lorentz connection \(\omega^{ab}\), plus the installed SM fields. Not a lattice defect medium |
| Vacuum / dynamics | Einstein-Cartan: torsion algebraic, locked to fermionic spin, zero in vacuum. Linearized metric dynamics claimed equivalent to the entanglement first law in holographic AdS only |
| Particle | SM field excitation. Not a soliton and not a topological defect |
| Charge | SM hypercharge / electric charge, installed, not derived |
| Derrick escape | Not applicable: no soliton |
| Clock | No native de Broglie clock. Dirac / QFT time evolution is installed |
| EM | Maxwell from the installed \(U(1)_{\mathrm{EM}}\), not an emergent tilt mode |
| Quantum | SM quantum field theory is input. Geometric superposition of the metric is a holographic theorem in AdS, open beyond |
| Gravity | Einstein-Cartan with cosmological constant. The program's claim is that this gravitational sector is selected by entanglement bookkeeping |
| Free parameters | The full SM (gauge group, three generations, Yukawas, Higgs, \(G\), \(\Lambda\)). Neutrino Dirac vs Majorana open. No extra knob in the Hehl-Datta coefficient |
| Lab anchor | \(G\); torsion / Lorentz-violation bounds (Kostelecky-Russell-Tasson; Shapiro reviews) |
| Formal artifacts | Paper series + this repo's M9.1 algebra scripts. \(I_B\) campaign scripts exist and are documented negatives, not a multi-digit coefficient |
| Next falsifier | A null Bose / Marletto-Vedral gravitationally-mediated-entanglement result at sufficient sensitivity (off-lattice). Near-term in-platform: none sharper than the HD algebra itself |

## Decision-Relevant Attributes

| Attribute | M9 |
| --- | --- |
| Free parameters | SM content is input. The HD coefficient is not fit. The load-bearing choice is "Einstein-Cartan is the unique modular-selected coupling," which M9.1 does not test |
| Honest residuals | Multi-digit \(I_B\) open and now a documented negative on three routes; de Sitter open; nonlinear Einstein-Cartan open as a positive theorem and closed as an axial obstruction; UV completion open by construction; no dark-matter particle; masses and generations unexplained. Author Final Status table, 14 August 2026 |
| Formal artifacts | M9.1 extracts \(\mathcal{L}_{\mathrm{int}}/(-\kappa J_5\cdot J_5)\) by stationarity, compares to \(3/16\) only after extraction |
| Falsifiable near-term tests | HD contact is Planck-suppressed (author: undetectable at colliders). GME experiments are the live gravitational-sector falsifier and are not an OpenWave lattice test |

## Field Configuration of Particles

Standing demand: state the field configuration of each particle, and whether it
uses topological vortices. This column's honest answer is that it does **not**
supply defect configurations. Particles are the installed SM fields.

| Particle | Configuration in NSM | Topological vortex? |
| --- | --- | --- |
| Electron / leptons / quarks | Dirac (or Weyl) SM spinors, minimally coupled to \(e,\omega,A,H\) | No |
| Gauge bosons | SM connections; field strengths \(F=dA+A\wedge A\), no torsion coupling | No |
| Higgs | Complex scalar doublet, no spin, does not source torsion | No |
| Graviton / torsion | Metric perturbation of Einstein-Cartan; torsion non-propagating | No |

## Implementation Status

Proposed column. No `MODELS.md` icons have been moved. Almost every particle
row will stay 🚧 because matter is installed. Closed attempts are recorded as
negatives; they are not cells.

| Sector | Status |
| --- | --- |
| Hehl-Datta contact coefficient \(3\kappa/16\) | M9.1 gate PASS: on-shell ratio \(3/16\) in both signatures (scatter \(<10^{-15}\)), second-method audit CONFIRMED. Paper's printed \(s=-\frac14\varepsilon J_5\) FAIL: measured \(-\frac12\). Note: [`research/findings/m9_1_hehl_datta_note.md`](research/findings/m9_1_hehl_datta_note.md) |
| Gravity: Newton limit (GEM) | 🚧 pre-registered 2026-08-15, **no run**. Attractive \(1/r^2\) from 3-d Poisson / geodesic of inherited Einstein, not the M5 GEM route. Gates locked: [`research/tasks/m9_2_task_details.md`](research/tasks/m9_2_task_details.md) |
| Gravity: metric phenomena | 🚧 domain note written, **no cell**. Einstein+\(\Lambda\); FGHMV cited not re-proved; dS `[O]`; torsion out of scope. [`research/findings/m9_metric_phenomena_note.md`](research/findings/m9_metric_phenomena_note.md) |
| \(I_B\) multi-digit coefficient | ❌ documented negative (campaign, not a `MODELS.md` row). Hard-cutoff residue moves (`FAILED_MULTI_DIGIT`); hole-scheme \(r\) is source-dependent (`NOT_UNIVERSAL`); Mittag-Leffler / polygamma \(H(\tau)\) is complete and not proportional to a local kernel (`HADAMARD_COMPLETE_NOT_UNIVERSAL`). Scripts: `m9_2_ib_hadamard.py`, `m9_3_ib_analytic.py`, `m9_4_ib_hadamard_complete.py` |
| Second-order Einstein-Cartan from entanglement | ❌ documented negative as a positive EC theorem. Metric Einstein through second order is cited (FHHPRV 2017). Axial matching is obstructed: CFT \(\langle J_5 J_5\rangle\) is nonlocal, algebraic torsion has no kinetic term (`STRUCTURE_ONLY`). Script: `m9_5_ec_symplectic.py`. Paper 14 |
| Torsion (status, not a cell) | EC theorems: algebraic, vacuum-vanishing, non-propagating. Spacetime HD \(\sim G\) is not lab-visible. Spintronic Berry / SOC is not Palatini \(\omega\). A late-FLRW spin average is an estimate, not a cosmological no-go |
| Charge, masses, clock, \(\mu\), spectrum, confinement, weak decays, DM | 🚧 not derived. Do not score as emergence |
| FGHMV / Condition NL | Not in-platform. Author `[P]`/`[O]` only |

## Roadmap

| Task | What lands |
| --- | --- |
| M9.1 | Certification gate: independent Hehl-Datta elimination (closed 2026-08-15) |
| M9.2 | Linearized Einstein / Newton \(1/r^2\). Pre-registered; not run |
| M9.3 | Gravity-metric note: Einstein+\(\Lambda\), domain labels, no holographic overclaim (written) |
| Rendering | Not before a gravity cell exists. Headless first |

Full row preview: [`research/m9_roadmap.md`](research/m9_roadmap.md).

## Help Wanted

| Contribution | What it would settle |
| --- | --- |
| Application discussion | Posted: [discussion #442](https://github.com/openwave-labs/openwave/discussions/442). Maintainer admission and official ID still wanted |
| Hostile recompute of Papers IV-VII | Whether Condition NL and the pure-information HD magnitude survive a second implementation |
| Grid Newton-limit script | Execute the locked M9.2 gates. That, not this briefing, is what could earn the first gravity cell |

Flow: [discussion #442](https://github.com/openwave-labs/openwave/discussions/442)
→ [PR #441](https://github.com/openwave-labs/openwave/pull/441) with DCO
(`git commit -s`). Start here:
[`ONBOARDING_MODELS.md`](../../../ONBOARDING_MODELS.md),
[`CONTRIBUTING.md`](../../../CONTRIBUTING.md).

## Rich Context for Deep Reader

Spec of record: [`research/m9_theory_canonical.md`](research/m9_theory_canonical.md).
Application body: [`research/APPLICATION.md`](research/APPLICATION.md).
Citations: [`theory/_CITATIONS.md`](theory/_CITATIONS.md).
M9.1 method note: [`research/findings/m9_1_hehl_datta_note.md`](research/findings/m9_1_hehl_datta_note.md).
Metric domain note: [`research/findings/m9_metric_phenomena_note.md`](research/findings/m9_metric_phenomena_note.md).
Author-gated questions stay with the author.
