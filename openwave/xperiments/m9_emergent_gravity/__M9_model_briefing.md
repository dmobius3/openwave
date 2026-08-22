# M9 NSM (Emergent Gravity): Model Briefing

> **What M9 brings.** Gravity as entanglement bookkeeping, assembled as the
> Standard Model minimally coupled to Einstein-Cartan gravity. Matter is
> installed, not emerged. The one interaction beyond SM + Einstein gravity is
> the Hehl-Datta axial-axial contact, coefficient fixed by \(G\). This is a
> gravity-certification column, closer to M8 than to M5. Model ID M9 was
> assigned at admission (discussion #442, 2026-08-15).
>
> **Archive pointer (2026-08-22).** The campaign records this document cites (M9.3 to M9.73: notes, scripts, data, latex, APPLICATION) are not in the tree; they are parked in PR #441's branch history at `a5640709` and return as per-task PRs against registered roadmap rows. Their links below were reduced to plain text at the merge trim.

## Identity

| Field | Value |
| --- | --- |
| Model ID | M9 |
| Name | NSM (New Standard Model): gravity as entanglement bookkeeping on SM + Einstein-Cartan |
| Author | Dr. Robert W. McGwier, PhD, CTO, Cohere Technology Group (sole author) |
| Author contact | GitHub [@n4hy](https://github.com/n4hy), for author-gated questions (definitions, intent, what the model does and does not claim); routing convention in [`dev_docs/CROSS_MODEL_TESTING.md`](../../../dev_docs/CROSS_MODEL_TESTING.md) § 6 |
| Lineage | Faulkner-Guica-Hartman-Myers-Van Raamsdonk 2014; Casini-Huerta-Myers; Einstein-Cartan-Sciama-Kibble; Hehl-Datta 1971; Jacobson entanglement equilibrium |
| Primary sources | Author repo [github.com/n4hy/New_Model_Emergent_Gravity](https://github.com/n4hy/New_Model_Emergent_Gravity) (CC-BY-4.0 PDFs under `research/`); registry in [`theory/_CITATIONS.md`](theory/_CITATIONS.md). Specification of the first task: Paper III action (2) |
| In-repo | Headless only. No launcher. First task: [`research/scripts/hehl_datta.py`](research/scripts/hehl_datta.py) |
| Application | [discussion #442](https://github.com/openwave-labs/openwave/discussions/442); `research/APPLICATION.md`; first PR [openwave-labs/openwave#441](https://github.com/openwave-labs/openwave/pull/441) |

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
| Honest residuals | Parked campaign records, not in the tree (`a5640709`), author-side statements: multi-digit \(I_B\) documented negative; FGHMV copy onto the cosmological horizon obstructed; the spinless vacuum of Einstein+\(\Lambda>0\) \emph{is} de Sitter (metric, not a holographic selection); nonlinear EC closed as an axial obstruction; UV pair not selected; an \(S^3\)/Bloch picture is a guess; as a description of \(K\) it fails (Paper 32); point-source \(\delta S\) vs site-energy is not a CHM first law even in 1d (Paper 36); a valid fixed-\(H\) probe still selects enclosed energy in 1d and CHM on small 3d balls (Paper 37); a non-sphere native weight beats exported CHM on cubes (Paper 38, guess not a theorem); hop-area and cut-correlator area are energy proxies, not a Clausius \(\eta\) (Papers 39--43); finite \(\delta S\) is enclosure of \(2h(\alpha)\), not Clausius (Papers 63--64); the sea leftover is not an area term (Paper 65); hop stress of the transfer is dust, not \(\Lambda\) (Paper 66); a thermal scan is not de Sitter (Paper 67); Gauss as used is Stokes \(g=-M/A\), not a hop force (Paper 68); pair direction is \(\nabla_c M_{AB}\) (Paper 69); no DM particle; masses and generations unexplained |
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
| Gravity: Newton limit (GEM) | ⚠️ inherited, C1 PASS / C2 FAIL. Attractive \(GM/r^2\) to \(\le 1.3\%\) on a Dirichlet cube; isolated \(\Phi=-GM/r\) fails (\(26\)--\(35\%\)) because the locked box is not isolated space. Parked campaign records, not in the tree (`a5640709`): M9.38 sources the same C1 from an entanglement \(M_{\mathrm{hat}}\) (\(n=65\) residuals \(\le 3.0\%\)); Poisson is still inherited, not derived. M9.39: the actual \(\delta e\) map sources far-field Newton with mass \(\sum\delta e\), not \(M_{\mathrm{hat}}\) (C_hat FAIL \(41\)--\(43\%\); auditor CONFIRMED). M9.40--42: two packets, midpoint cancel, unequal-mass \(M/r^2\) null (auditor precision REFUTED). Wide-source first law grows; compact plateaus; not Newtonian \(\Lambda\). Not GEM, not FGHMV, no `MODELS.md` column. `research/findings/m9_42_bary_note.md` |
| Gravity: metric phenomena | 🚧 parked campaign record, not in the tree (author-side status was 🚧): domain note written, **no cell**. Einstein+\(\Lambda\); FGHMV cited not re-proved; torsion out of scope. `research/findings/m9_metric_phenomena_note.md` |
| de Sitter vacuum of the metric sector | 🚧 parked campaign record, not in the tree (author-side status was ✅): `[P]` as Einstein+\(\Lambda\): spinless vacuum with \(\Lambda>0\) is de Sitter (\(a=e^{Ht}\), \(H^2=\Lambda/3\)). Torsion is not in the problem. Not a holographic selection. `research/latex/19_deSitter_Is_the_Vacuum.tex` |
| de Sitter / FGHMV-standard cosmology | 🚧 parked campaign record, not in the tree (author-side status was ❌): FGHMV copy obstructed (M9.6). Metric bar only: sign opposite, isometries too few. Torsion is not the missing piece. `research/findings/m9_6_ds_closure_note.md` |
| Cosmological first-law sign | 🚧 parked campaign record, not in the tree (author-side status was ✅): `[P]` derived from SdS (M9.9): \(T\mathrm{d}S+\mathrm{d}M=0\), \(\mathrm{d}r_c/\mathrm{d}M\|_0=-1\). Minus sign is Einstein+\(\Lambda\), not an AdS import. `research/findings/m9_9_sds_sign_note.md` |
| A2 local-\(X\) modular ansatz | 🚧 parked campaign record, not in the tree (author-side status was ⚠️): local hop is CHM-type (Paper 25). Horizon first law: a radial weight beats flat (Papers 27, 29); the parabola is not selected over \(R-r\) (Paper 29, tie). Tracking floors fail. A 3d point *Hamiltonian* source selects enclosed energy (Paper 35) but fails the vacuum first law (Paper 36). At fixed \(H\), the first law holds; 1d prefers enclosed energy, 3d balls prefer CHM (Paper 37). Off the ball, a shape-native weight beats an exported CHM (Paper 38); that is a measurement of this probe, still a guess about gravity. A hop conformal bump makes \(\delta S\) track hop-length area (Paper 39) but \(\eta=\delta S/\delta A\) is not constant. Dual-face and weak-field area are linearly the same test (Paper 41). Hop-area is collinear with CHM energy (\(\rho=0.996\), Paper 42). A cut-correlator area moves at fixed \(H\) but is still an energy proxy (Paper 43). Two masses well inside \(R=3\) balls: \(\delta S\) tracks enclosed energy (Paper 45) with a universal \(\kappa\approx 0.97\) (Paper 46). \(\kappa\) weighs enclosed mass and locates the source (Paper 47). Feeding \(M_{\mathrm{hat}}\) into the M9.2 DST Poisson solver passes inherited C1 (Paper 48). The actual \(\delta e\) map sources far-field Newton with \(\sum\delta e\), not \(M_{\mathrm{hat}}\) (Paper 49). Two real packets: one enclosing ball reads the pair mass; the midpoint cancels; exterior matches two-point Coulomb (Paper 50). Wide vs compact: the first law grows with an extended source and plateaus on a star (Paper 51); that is not Newtonian \(\Lambda\). Unequal masses: the interior null is \(M/r^2\), not the centre of mass (Paper 52). First-law \(\kappa\) on the \(3+1\)D diamond waist moves \(0.6\%\) when the staggered mass is turned on (Paper 53). Periodic band-edge transfer: uniform \(\delta e\) to \(10^{-12}\), volume first law, inherited \(a\propto r\) (Paper 54). First-law mass plus Gauss, no Poisson solver: star slope \(-1.997\), sea slope \(+1.266\) (Paper 55). \(\sigma\)-scan: slope interpolates \(-1.87\to +1.00\); every \(a\) is inward. Paper 54's ``Newtonian \(\Lambda\)'' is withdrawn --- this is dust, not de Sitter (Paper 56). Fermi-sea vacuum: \(E_{\mathrm{vac}}<0\) but \(S\) is area-law, not a volume \(\Lambda\) (Paper 57). The first law tracks \(\delta e\), not raw \(e\): the sea is subtracted and does not gravitate (Paper 58). The complement of an enclosing ball is not a cosmological horizon: \(\delta S(B^c)/\delta S(B)\sim 10^{-4}\) (Paper 59). \(\kappa(\alpha)\) runs \(1.34\to 0.75\); it is \(2h(\alpha)/(\alpha\Delta E)\), not a coupling (Paper 60). Gauss rebuilt from \(\sum\delta e\) only: star \(-1.998\), sea \(+0.967\), independent of \(\alpha\) (Paper 61). Finite \(\delta S\) on 216 balls tracks \(P_{\mathrm{flat}}\) at every \(\alpha\); CHM and \(\mathrm{Tr}(K\Delta C)\) lose (Paper 62). \(\delta S/S_{\mathrm{global}}\) tracks \(P_{\mathrm{flat}}/M_{\mathrm{global}}\) (\(\rho=0.991\)--\(0.999\)); \(f=1\) only on the centered ball; source-inside balls leak (Paper 63). Sea-transfer \(\delta S\) is extensive: slab grow \(2.55\) (area fixed); not a packet artifact (Paper 64). The leftover \(2.5\neq 3\) is not a universal area term: two-term \(\delta S=aV+bA\) does not beat volume, and held-out rods miss by \(5\times\) (Paper 65). Hop kinetic stress of the transfer: sea \(p/\rho=0\), star \(0.16\); not \(\Lambda\) (Paper 66). Thermal excess \(r(T)\) runs and is lattice-dependent; a one-point crossing of \(-1\) is not de Sitter (Paper 67). Unique isotropic flux is \(g=-M/A\); hop \(E_{\mathrm{int}}\) is overlap, not \(1/r\) (Paper 68). Two-packet direction is \(\hat n=\nabla_c M_{AB}\); it tracks \(M/r^2\) superposition, not the CM (Paper 69). Exact open-hop basis, no LAPACK: the same \(10.746^\circ\) (Paper 70). Pair mass maps add; derived \(\mathbf g\)s do not, except far away (Paper 71). Pair \(S_{\mathrm{global}}=2h(\alpha_A)+2h(\alpha_B)\); \(\delta S\) tracks \(M_{AB}\); \(\nabla\cdot\mathbf g\) is not Poisson (Paper 72). Midpoint \(R=5\) pair ball: \(f_S=f_E=0.99\); a site-inside axial offset still leaks (Paper 73). `research/findings/m9_63_enclose_pair_note.md` |
| A1 UV coefficient | 🚧 parked campaign record, not in the tree (author-side status was ✅): `[P]` for this 1d fermion (\(\alpha=0.323\)) and for the \(3+1\)D diamond area law (\(\alpha=0.245\), UV drift \(\le 4.3\%\), auditor CONFIRMED). Not \(1/4G\), not foam. `research/findings/m9_14_A1_diamond_note.md` |
| Jacobson as Q2 substitute | 🚧 parked campaign record, not in the tree (author-side status was ❌): not `[P]` (M9.7). 1995 \(\Rightarrow\) Einstein, \(\Lambda\) free, no HD. 2016 conformal half does not apply to the SM (\(b_3=-7\)). `research/findings/m9_7_jacobson_note.md` |
| \(I_B\) multi-digit coefficient | 🚧 parked campaign record, not in the tree (author-side status was ❌): documented negative (campaign, not a `MODELS.md` row). Hard-cutoff residue moves (`FAILED_MULTI_DIGIT`); hole-scheme \(r\) is source-dependent (`NOT_UNIVERSAL`); Mittag-Leffler / polygamma \(H(\tau)\) is complete and not proportional to a local kernel (`HADAMARD_COMPLETE_NOT_UNIVERSAL`). Scripts: `m9_2_ib_hadamard.py`, `m9_3_ib_analytic.py`, `m9_4_ib_hadamard_complete.py` |
| Second-order Einstein-Cartan from entanglement | 🚧 parked campaign record, not in the tree (author-side status was ❌): documented negative as a positive EC theorem. Metric Einstein through second order is cited (FHHPRV 2017). Axial matching is obstructed: CFT \(\langle J_5 J_5\rangle\) is nonlocal, algebraic torsion has no kinetic term (`STRUCTURE_ONLY`). Script: `m9_5_ec_symplectic.py`. Paper 14 |
| Torsion (status, not a cell) | EC theorems: algebraic, vacuum-vanishing, non-propagating. Spacetime HD \(\sim G\) is not lab-visible. Spintronic Berry / SOC is not Palatini \(\omega\). A late-FLRW spin average is an estimate, not a cosmological no-go |
| Charge, masses, clock, \(\mu\), spectrum, confinement, weak decays, DM | 🚧 not derived. Do not score as emergence |
| UV completion (selected pair) | 🚧 parked campaign record, not in the tree (author-side status was ❌): selection-uniqueness answered negative (M9.5): SM is not a holographic CFT (\(b_i\neq 0\), 118 dof); certified first law is blind to \(G_{\mathrm{SM}}\). Existence of some other pair still `[O]`. `research/findings/m9_5_q4a_pair_note.md` |
| UV axial deformation (Q4b) | 🚧 parked campaign record, not in the tree (author-side status was ⚠️): candidate, not a cell. Unique quadratic massive axial action recovers \(r=3/16\) as \(M\to\infty\) and is Yukawa (not a contact) at finite \(M\). Audited. Not EC at finite \(M\). `research/findings/m9_4_uv_deformation_note.md` |
| FGHMV / Condition NL | Not in-platform. Author `[P]`/`[O]` only |

## Roadmap

| Task | What lands |
| --- | --- |
| M9.1 | Certification gate: independent Hehl-Datta elimination (closed 2026-08-15) |
| M9.2 | Linearized Einstein / Newton \(1/r^2\). Closed: C1 PASS / C2 FAIL, inherited Einstein, not GEM |
| M9.3 | Gravity-metric domain note: parked campaign record (`a5640709`), returns as a per-task PR |
| M9.4 | Axial UV deformation Q4b: parked campaign record (`a5640709`), returns as a per-task PR |
| M9.5 | Q4a pair selection: parked campaign record (`a5640709`), returns as a per-task PR |
| M9.6 | de Sitter at FGHMV standard: parked campaign record (`a5640709`), returns as a per-task PR |
| M9.7 | Jacobson as Q2 substitute: parked campaign record (`a5640709`), returns as a per-task PR |
| Rendering | Not before a gravity cell exists. Headless first |

Full row preview: [`research/m9_roadmap.md`](research/m9_roadmap.md).

## Help Wanted

| Contribution | What it would settle |
| --- | --- |
| Application discussion | Posted: [discussion #442](https://github.com/openwave-labs/openwave/discussions/442). Maintainer admission and official ID still wanted |
| Hostile recompute of Papers IV-VII | Whether Condition NL and the pure-information HD magnitude survive a second implementation |
| Official ID / MODELS.md | M9.2 C1 passed as inherited Newton; C2 failed. A cell needs a maintainer ID |

Flow: [discussion #442](https://github.com/openwave-labs/openwave/discussions/442)
→ [PR #441](https://github.com/openwave-labs/openwave/pull/441) with DCO
(`git commit -s`). Start here:
[`ONBOARDING_MODELS.md`](../../../ONBOARDING_MODELS.md),
[`CONTRIBUTING.md`](../../../CONTRIBUTING.md).

## Rich Context for Deep Reader

Spec of record: [`research/m9_theory_canonical.md`](research/m9_theory_canonical.md).
Application body: `research/APPLICATION.md`.
Citations: [`theory/_CITATIONS.md`](theory/_CITATIONS.md).
M9.1 method note: [`research/findings/m9_1_hehl_datta_note.md`](research/findings/m9_1_hehl_datta_note.md).
Metric domain note: `research/findings/m9_metric_phenomena_note.md`.
UV deformation (Q4b): `research/findings/m9_4_uv_deformation_note.md`.
Q4a pair selection: `research/findings/m9_5_q4a_pair_note.md`.
de Sitter obstruction: `research/findings/m9_6_ds_closure_note.md`.
Jacobson is not a [P] substitute: `research/findings/m9_7_jacobson_note.md`.
Author-gated questions stay with the author.
