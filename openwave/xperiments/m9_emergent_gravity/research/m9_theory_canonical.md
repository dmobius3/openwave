# M9 Emergent Gravity / NSM — theory canonical (draft)

> Proposed column. Canonical when later documents disagree. This draft records
> only what the first task needs: the effective action, the conventions, and
> what is *not* claimed. Holographic certificates stay outside this page until
> they have in-platform scripts. Author papers:
> [github.com/n4hy/New_Model_Emergent_Gravity](https://github.com/n4hy/New_Model_Emergent_Gravity).

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
  ([`findings/m9_13_A2_diamond_note.md`](findings/m9_13_A2_diamond_note.md)).
  That is not \(a\to 0\), not the SM, and not a selection of de Sitter.
- A1 (UV entanglement coefficient) is IR-stable for the 1d fermion
  log and for the \(3+1\)D diamond area law
  ([`findings/m9_14_A1_diamond_note.md`](findings/m9_14_A1_diamond_note.md)).
  That is not \(\eta=1/4G\) and not mean-zero foam.
- The local modular *hop* on the diamond waist tracks the CHM
  envelope (\(\rho=-0.987\);
  [`findings/m9_15_chm_shape_note.md`](findings/m9_15_chm_shape_note.md)).
  That is not a derivation of \(G\).
- FGHMV-standard de Sitter closure is obstructed (sign and isometries:
  [`findings/m9_6_ds_closure_note.md`](findings/m9_6_ds_closure_note.md)).
  Einstein+\(\Lambda\) from a cosmological CFT is not a theorem of this
  program. Jacobson is not a [P] substitute
  ([`findings/m9_7_jacobson_note.md`](findings/m9_7_jacobson_note.md)). Nonlinear Einstein-Cartan as a positive theorem, and the
  multi-digit coefficient of \(I_B\), remain open.
  Q4a selection-uniqueness of the holographic pair is answered in the
  negative ([`findings/m9_5_q4a_pair_note.md`](findings/m9_5_q4a_pair_note.md));
  existence of some pair is still open. A quadratic axial deformation
  (Q4b) is recorded in
  [`findings/m9_4_uv_deformation_note.md`](findings/m9_4_uv_deformation_note.md).
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
  M9.2 may propose a Newton-limit cell only after its locked gates are run.
- Official M9.2 is the Newton-limit task. Campaign files
  `m9_2_ib_hadamard.py`, `m9_3_ib_analytic.py`,
  `m9_4_ib_hadamard_complete.py`, and `m9_5_ec_symplectic.py` are
  documented negatives, not gravity cells.
