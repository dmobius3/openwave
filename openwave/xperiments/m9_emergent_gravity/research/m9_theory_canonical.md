# M9 Emergent Gravity / NSM — theory canonical (draft)

> Proposed column. Canonical when later documents disagree. This draft records
> only what the first task needs: the effective action, the conventions, and
> what is *not* claimed. Holographic certificates stay outside this page until
> they have in-platform scripts. Author papers:
> [github.com/n4hy/New_Model_Emergent_Gravity](https://github.com/n4hy/New_Model_Emergent_Gravity).
>
> **Archive pointer (2026-08-22).** The campaign records this document cites (M9.3 to M9.73: notes, scripts, data, latex, APPLICATION) are not in the tree; they are parked in PR #441's branch history at `a5640709` and return as per-task PRs against registered roadmap rows. Their links below were reduced to plain text at the merge trim.

## 1. Arena

Effective field theory on a 4-dimensional orientable spacetime. Gravitational
variables are the coframe \(e^a\) and an independent Lorentz connection
\(\omega^{ab}\). Matter is the empirically installed Standard Model (gauge
group, representations, Yukawas, Higgs, three generations). Particles are not
emergent defects of a lattice field.

## 2. Effective action (specification of the first task)

\[
S_{\mathrm{NSM}}
=
\frac{1}{4\kappa}
\int \varepsilon_{abcd}\, e^a\wedge e^b\wedge R^{cd}
-
\frac{\Lambda}{6}
\int \varepsilon_{abcd}\, e^a\wedge e^b\wedge e^c\wedge e^d
+
\int e\,\mathcal{L}_{\mathrm{SM}}[e,\omega,\psi,A,H],
\]

with \(\kappa = 8\pi G\), \(R^{ab} = d\omega^{ab} + \omega^a{}_c\wedge\omega^{cb}\),
and \(\varepsilon^{0123} = +1\). The 4-form Palatini term equals
\(\int e\, R/(2\kappa)\) when torsion vanishes.

Fermions couple through the full connection

\[
D_\mu
=
\partial_\mu
+
\frac{1}{4}\omega_\mu^{ab}\gamma_{ab}
+
ig\cdot A_\mu,
\qquad
\gamma_{ab} := \tfrac{1}{2}[\gamma_a,\gamma_b].
\]

Gauge field strengths use the exterior derivative alone (torsion does not enter
\(F\)). Scalars do not source torsion.

The Hermitian Dirac kinetic term used by the certification scripts is

\[
\mathcal{L}_D
=
\frac{i}{2}\,e\,
\bar\psi\gamma^\mu\overleftrightarrow{D}_\mu\psi
-
e\,m\,\bar\psi\psi.
\]

## 3. Algebraic Cartan equation

Independent variation of \(\omega\) is algebraic. For a totally antisymmetric
Dirac spin tensor the traces vanish and

\[
T^{\lambda\mu\nu} = \kappa\, s^{\lambda\mu\nu}.
\]

Contorsion of a totally antisymmetric torsion is \(K = T/2\).

## 4. Claim under test (M9.1)

Eliminating \(T\) produces Riemannian Einstein-Hilbert plus the Hehl-Datta
contact term

\[
\mathcal{L}_{\mathrm{HD}}
=
-\frac{3\kappa}{16}\,J_5^\mu J_{5\mu},
\qquad
J_5^\mu = \sum_f \bar\psi_f\gamma^5\gamma^\mu\psi_f.
\]

M9.1 tests this algebra. It does not test why Einstein-Cartan was selected.

## 5. Conventions locked for M9.1

| Item | Lock |
| --- | --- |
| \(\kappa\) | \(8\pi G\) |
| \(\varepsilon^{0123}\) | \(+1\) |
| \(\gamma^5\) (mostly minus) | \(i\gamma^0\gamma^1\gamma^2\gamma^3\) |
| Dirac kinetic term | Hermitian, as above |
| Reported observable | dimensionless ratio \(\mathcal{L}_{\mathrm{int}}/(-\kappa\, J_5\cdot J_5)\) |

## 6. Known tensions / declared opens (author)

- Matter content is empirical input, not derived.
- Linearized Einstein-Cartan from entanglement is claimed only in holographic AdS.
- A2 (local modular kernel) is not refuted on free lattice Dirac in
  \(d=1,2,3\) or on a \(3+1\)D diamond waist at two spacings
  (`findings/m9_13_A2_diamond_note.md`).
  That is not \(a\to 0\), not the SM, and not a selection of de Sitter.
- A1 (UV entanglement coefficient) is IR-stable for the 1d fermion
  log and for the \(3+1\)D diamond area law
  (`findings/m9_14_A1_diamond_note.md`).
  That is not \(\eta=1/4G\) and not mean-zero foam.
- The local modular *hop* on the diamond waist tracks the CHM
  envelope (\(\rho=-0.987\);
  `findings/m9_15_chm_shape_note.md`).
  That is not a derivation of \(G\).
- The entanglement first law on hop perturbations tracks a local
  kernel; a flat NN kernel beats the CHM envelope
  (`findings/m9_16_first_law_note.md`).
  Clausius here does not select the geometric weight on the
  mixed-hop set. On the cut alone, C3 prefers CHM to flat
  (`findings/m9_17_18_horizon_first_law_note.md`);
  tracking floors still fail. On a second radius the parabola is
  not selected over a linear \(R-r\) weight
  (`findings/m9_20_horizon_shape_note.md`).
  The occupancy-stable hop probe dies at \(R=5\)
  (`findings/m9_21_larger_horizon_note.md`).
  That is not a Planck-scale result. Half-filling scores C4 at
  \(R=5\) as a tie
  (`findings/m9_22_halffill_horizon_note.md`).
- A Bloch-like dimer coordinate does not predict the modular hop;
  CHM does, on a complete covering of the ball
  (`findings/m9_23_bloch_note.md`).
- At fixed \(H\), region shape changes \(S\) beyond cut area
  (`findings/m9_24_region_deform_note.md`).
- \(\delta S\) of a family of balls is a linear functional of
  local energy; the kernel is not uniquely CHM
  (`findings/m9_25_linear_functional_note.md`).
  That is not linearized Einstein. A point source makes the
  kernel \emph{enclosed energy}, not CHM
  (`findings/m9_26_point_source_note.md`).
  The same operators fail the 1d CHM theorem case
  (`findings/m9_27_1d_point_chm_note.md`):
  \(\delta S=\mathrm{Tr}(K_{\mathrm{mid}}\Delta C)\), not a
  local functional of vacuum modular energy. Paper 35 is not
  a 3d no-go. At fixed \(H\), a localized occupation transfer
  *does* obey \(\delta S=\mathrm{Tr}(K_{\mathrm{vac}}\Delta C)\);
  1d still selects enclosed energy, 3d balls select CHM
  (`findings/m9_28_fixedh_state_note.md`).
  Off the ball, a shape-native weight beats an exported
  CHM kernel
  (`findings/m9_29_shape_note.md`).
  That measurement does not make ``need a shape'' a gravity
  theorem. A hop conformal bump makes \(\delta S\) track a
  hop-length area
  (`findings/m9_30_area_note.md`)
  without a constant Clausius \(\eta\). Dual-face area is
  linearly the same test
  (`findings/m9_31_proper_area_note.md`).
  Hop-area is collinear with CHM energy
  (`findings/m9_32_two_term_note.md`).
  A cut-correlator area at fixed \(H\) is still an energy
  proxy
  (`findings/m9_33_state_area_note.md`).
  Two masses well inside a ball: \(\delta S=\kappa M_{\mathrm{enc}}\)
  with \(\kappa\approx 0.97\) universal to \(2\%\)
  (`findings/m9_36_kappa_note.md`).
  That \(\kappa\) weighs enclosed mass and locates the source
  (`findings/m9_37_weigh_note.md`).
  Feeding \(M_{\mathrm{hat}}\) into the locked DST Poisson
  solver passes inherited C1
  (`findings/m9_38_from_kappa_note.md`).
  The actual \(\delta e\) map sources far-field Newton with
  mass \(\sum\delta e\), not \(M_{\mathrm{hat}}\)
  (`findings/m9_39_tail_note.md`).
  Two packets: one enclosing ball and the one-mass \(\kappa\)
  read the pair; the midpoint cancels
  (`findings/m9_40_pair_note.md`).
  An extended source keeps growing; a star plateaus
  (`findings/m9_41_uniform_note.md`).
  Unequal masses cancel at the inverse-square point, not
  the centre of mass
  (`findings/m9_42_bary_note.md`).
  The same \(\kappa\) on the \(3+1\)D diamond waist moves
  \(0.6\%\) under staggered mass
  (`findings/m9_43_diamond_kappa_note.md`).
  A periodic band-edge transfer is a valid fermion state
  with uniform \(\delta e\); inherited Newton of that
  density is \(a\propto r\)
  (`findings/m9_44_uniform_pbc_note.md`).
  First-law \(M(R)\) plus spherical Gauss, with no Poisson
  solver, already gives \(1/R^2\) for a star and \(a\propto R\)
  for the sea
  (`findings/m9_45_gauss_force_note.md`).
  The Gauss slope interpolates with packet width; every
  \(a\) is inward. That is dust, not \(\Lambda\)
  (`findings/m9_46_sigma_note.md`).
  The unperturbed Fermi-sea vacuum has \(E_{\mathrm{vac}}<0\)
  but area-law \(S\); the first law does not promote
  \(E_{\mathrm{vac}}\) to \(\Lambda\)
  (`findings/m9_47_vacuum_note.md`).
  \(\delta S\) tracks \(\delta e\), not raw \(e\): the sea
  is subtracted and is not repulsive Newton
  (`findings/m9_48_subtract_note.md`).
  The complement of an enclosing ball is not a cosmological
  horizon (\(\delta S(B^c)\) is four orders smaller)
  (`findings/m9_49_complement_note.md`).
  \(\kappa=\delta S/P_{\mathrm{flat}}\) runs with the
  transfer size \(\alpha\); it is \(2h(\alpha)/(\alpha\Delta E)\),
  not an \(\alpha\)-independent coupling
  (`findings/m9_50_alpha_note.md`).
  The reusable mass is \(\sum\delta e\).
  Gauss rebuilt from that energy alone keeps star
  \(1/R^2\) and sea \(a\propto R\), independent of \(\alpha\)
  (`findings/m9_51_energy_gauss_note.md`).
  Finite \(\delta S\) on the ball cloud tracks
  \(P_{\mathrm{flat}}\) at every \(\alpha\); CHM and
  \(\mathrm{Tr}(K\Delta C)\) lose
  (`findings/m9_52_predictor_note.md`).
  The two are the same fraction:
  \(\delta S/S_{\mathrm{global}}\) tracks
  \(P_{\mathrm{flat}}/M_{\mathrm{global}}\)
  (\(\rho=0.991\)--\(0.999\)); \(f=1\) only on the
  ball that actually encloses the packet
  (`findings/m9_53_enclose_note.md`).
  That fraction is not a packet artifact: on the
  periodic band-edge transfer, slabs hold area fixed
  and \(\delta S\) still grows with thickness
  (grow \(2.55\) vs volume \(3\), area \(1\))
  (`findings/m9_54_sea_enclose_note.md`).
  The leftover \(2.5\neq 3\) is not a universal
  \(aV+bA\) area term: volume-only is as good, and
  held-out rods miss by five
  (`findings/m9_55_twoterm_note.md`).
  The hop-dispersion virial of the same transfer
  has \(p/\rho=0\) on the sea and \(0.16\) on a
  star; it is dust, not \(\Lambda\)
  (`findings/m9_56_stress_note.md`).
  Heating that sea, \(\mu=0\), the excess \(r(T)\)
  runs and is lattice-dependent; a crossing of
  \(-1\) is not a cosmological constant
  (`findings/m9_57_thermal_note.md`).
  The radial field used from here is the unique
  isotropic flux completion \(g=-M/A\); the hop
  two-packet energy is overlap, not \(1/r\)
  (`findings/m9_58_flux_note.md`).
  For two packets the direction is
  \(\hat n=\nabla_c M_{AB}\); it tracks \(M/r^2\)
  superposition, not the centre of mass
  (`findings/m9_59_direction_note.md`).
  That direction is unchanged when the hop basis
  is the exact 1d product, with no eigensolve
  (`findings/m9_60_exact_note.md`).
  Pair mass maps add; the derived \(\mathbf g\)s
  do not, except far from both packets
  (`findings/m9_61_super_note.md`).
  Pair \(S_{\mathrm{global}}=2h(\alpha_A)+2h(\alpha_B)\);
  finite \(\delta S\) tracks \(M_{AB}\);
  \(\nabla\cdot\mathbf g\) is not a local Poisson law
  (`findings/m9_62_pair_ent_note.md`).
  On a midpoint \(R=5\) ball that holds the pair
  mass, \(f_S=f_E=0.99\); a ball that only
  contains both source sites can still leak
  (`findings/m9_63_enclose_pair_note.md`).
  That is Stokes plus a measured mass, not a
  derivation of Einstein.
- An \(S^3\) of radius \(\rho\) has \(\mathcal{R}=6/\rho^2\);
  \(\rho\to i\ell\) flips the sign
  (`findings/m9_19_s3_curvature_note.md`).
  Putting virtual modes on that sphere, and reading \(X_4\) as
  curvature, is a guess. It does not close Q2.

- FGHMV-standard de Sitter closure is obstructed (sign and isometries:
  `findings/m9_6_ds_closure_note.md`).
  Einstein+\(\Lambda\) from a cosmological CFT is not a theorem of this
  program. Jacobson is not a [P] substitute
  (`findings/m9_7_jacobson_note.md`). Nonlinear Einstein-Cartan as a positive theorem, and the
  multi-digit coefficient of \(I_B\), remain open.
  Q4a selection-uniqueness of the holographic pair is answered in the
  negative (`findings/m9_5_q4a_pair_note.md`);
  existence of some pair is still open. A quadratic axial deformation
  (Q4b) is recorded in
  `findings/m9_4_uv_deformation_note.md`.
- The Hehl-Datta term is a theorem of Einstein-Cartan plus Dirac (Hehl-Datta
  1971). The program's distinctive claim is the *selection* of Einstein-Cartan,
  which this page does not certify.

## 7. Consumption rules

- Cosmological questions are metric (Einstein+\(\Lambda\)). Torsion and
  Hehl-Datta are not cosmological limitations: they vanish in vacuum.
- Do not score a particle-row cell from this action: those fields are installed.
- Do not cite FGHMV, Condition NL, or \(I_B\) as in-platform results.
- A `MODELS.md` gravity cell requires a later lattice or analytic note of its
  own. M9.1 does not move one. M9.3 is a domain note and does not move one.
  M9.2 has been run: C1 PASS, C2 FAIL. No `MODELS.md` column is
  added without an official ID.
- Official M9.2 is the Newton-limit task. Campaign files
  `m9_2_ib_hadamard.py`, `m9_3_ib_analytic.py`,
  `m9_4_ib_hadamard_complete.py`, and `m9_5_ec_symplectic.py` are
  documented negatives, not gravity cells.
