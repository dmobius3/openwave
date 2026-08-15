# M9 metric-phenomena domain note

> Written domain statement for the proposed Emergent Gravity / NSM column.
> Equations first. This note is **not** a `MODELS.md` cell, not a lattice
> measurement, and not a new holographic theorem. Task record:
> [`../tasks/m9_3_task_details.md`](../tasks/m9_3_task_details.md).
> Newton \(1/r^2\) is a different task and is only pre-registered
> ([`../tasks/m9_2_task_details.md`](../tasks/m9_2_task_details.md)).

## 1. Metric sector of the effective action

The NSM gravitational variables are the coframe \(e^a\) and an independent
Lorentz connection \(\omega^{ab}\). The Palatini plus cosmological term is
the specification of record
([`../m9_theory_canonical.md`](../m9_theory_canonical.md) § 2):

\[
S_{\mathrm{met}}
=
\frac{1}{4\kappa}
\int \varepsilon_{abcd}\, e^a\wedge e^b\wedge R^{cd}
-
\frac{\Lambda}{6}
\int \varepsilon_{abcd}\, e^a\wedge e^b\wedge e^c\wedge e^d,
\qquad
\kappa=8\pi G,
\qquad
\varepsilon^{0123}=+1.
\]

When the spin density vanishes, the Cartan equation forces torsion to
zero and \(\omega=\mathring{\omega}(e)\). The 4-form then equals the
Einstein-Hilbert plus \(\Lambda\) action,

\[
S_{\mathrm{met}}\big|_{T=0}
=
\frac{1}{2\kappa}
\int e\bigl(\mathring{R}-2\Lambda\bigr),
\]

and the metric field equation is Einstein's equation with cosmological
constant,

\[
\mathring{G}_{\mu\nu}+\Lambda g_{\mu\nu}=\kappa\, T_{\mu\nu}.
\]

That is the content of this cell's *domain*: light, clocks, and
\(\Lambda\) are metric phenomena of Einstein-Hilbert plus \(\Lambda\).
They are not torsion phenomena.

## 2. Weak-field metric observables (definitions only)

Mostly-minus signature, static weak field, \(\lvert\Phi\rvert\ll 1\):

\[
\mathrm{d}s^2
=
-(1+2\Phi)\,\mathrm{d}t^2
+(1-2\Phi)\,\delta_{ij}\,\mathrm{d}x^i\mathrm{d}x^j,
\qquad
\nabla^2\Phi=4\pi G\rho.
\]

The three quantities named by the platform's simplest test for
**Gravity: metric phenomena** (light bending, time dilation, \(\Lambda\))
are then the standard Einstein-plus-\(\Lambda\) statements:

| Observable | Einstein+\(\Lambda\) content | Status in this note |
| --- | --- | --- |
| Light deflection by a mass \(M\) | \(\Delta\varphi=4GM/b\) | inherited definition; not a lattice run |
| Gravitational time dilation | \(\Delta\nu/\nu=\Delta\Phi\) | inherited definition; not a lattice run |
| Cosmological constant | \(G_{\mu\nu}+\Lambda g_{\mu\nu}=\kappa T_{\mu\nu}\); de Sitter radius \(\ell=\sqrt{3/\lvert\Lambda\rvert}\) when \(T=0\) and \(\Lambda>0\) | written as the action term; holographic *selection* of \(\Lambda\) in de Sitter is `[O]` |

No numerical value of \(\Delta\varphi\), \(\Delta\nu\), or \(\Lambda\) is
claimed here. A later script that measured those on a background would
be a different artifact.

## 3. What is cited, not re-derived: FGHMV (AdS, linearized)

Faulkner, Guica, Hartman, Myers, and Van Raamsdonk, JHEP 03 (2014) 051,
arXiv:1312.7856, prove that the first law of entanglement for ball-shaped
regions in a holographic CFT is equivalent to the linearized Einstein
equations about AdS. The author's Paper I invokes that theorem for the
metric sector. Epistemic status in the author's label table: `[P]`
**inside holographic AdS, at linear order**.

This note **cites** FGHMV. It does not re-prove it, does not implement
it on a lattice, and does not treat it as an in-platform OpenWave result.
The consumption rule in the canonical stands: do not cite FGHMV as an
in-platform result.

Through second order, Faulkner, Haehl, Hijano, Parrikar, Rabideau, and
Van Raamsdonk, JHEP 08 (2017) 057, arXiv:1705.03026, have Einstein's
equations from ball relative entropy. That is again a citation, and it
is Einstein, not Einstein-Cartan. The axial / Cartan matching at that
standard is an obstruction, recorded separately
(`scripts/m9_5_ec_symplectic.py`, Paper 14, `STRUCTURE_ONLY`).

## 4. de Sitter is `[O]`

Paper IX isolates the domain question. Every `[P]`-grade holographic
step uses ingredients that are special to AdS/CFT:

- the Casini-Huerta-Myers modular Hamiltonian for balls
- a global conformal structure on the boundary
- a clean asymptotic boundary and a holographic dictionary
- a Ryu-Takayanagi surface whose modular flow is a bulk Killing horizon

None of those is known in the same form in de Sitter. In particular
there is no accepted microscopic dual that would supply an analogue of
FGHMV, and static-patch modular Hamiltonians are not a CHM local
integral of \(T_{\mu\nu}\) in the same way. A future closure would have
to supply those objects and a first-law identity whose sign is derived,
not chosen.

This note does **not** invent a de Sitter theorem, a first-law sign, or
a dual. The cosmological extrapolation remains `[O]`.

## 5. Torsion is out of scope for this cell

Einstein-Cartan theorems, used here only to *exclude* torsion from the
metric-phenomena cell:

- the Cartan equation is algebraic: torsion is locked pointwise to spin
  density
- torsion vanishes wherever the spin density vanishes (vacuum, Maxwell,
  a spinless clock, a photon)
- torsion does not propagate: there is no torsional wave and no torsional
  substitute for \(\Lambda\)

Light bending, time dilation, and \(\Lambda\) are therefore read from
\(\mathring{g}_{\mu\nu}\) and \(\Lambda\), not from \(T^{\lambda\mu\nu}\).

What is **not** a theorem, and is not used as one:

- spacetime Hehl-Datta is \(\sim G\). It is not a laboratory field, in
  a collider or in a spintronic device
- Berry curvature and spin-orbit fields in spintronics are not the
  Palatini \(\omega\)
- a late-universe spatial average of spin density is an estimate, not a
  cosmological no-go and not a replacement for \(\Lambda\)

Those distinctions are recorded in Paper 14 § "What is proved about
torsion, and what is not". They do not move this cell.

## 6. Equation-to-document map

This note has no new solver. The action and the Cartan statement live in
the canonical; the holographic citations live in the author's papers.

| Object | Where |
| --- | --- |
| Palatini + \(\Lambda\) action | [`m9_theory_canonical.md`](../m9_theory_canonical.md) § 2 |
| Algebraic Cartan / vacuum \(T=0\) | same, § 3 |
| FGHMV linearized Einstein from balls | JHEP 03 (2014) 051; author's Paper I, Paper IX |
| FHHPRV second-order Einstein | JHEP 08 (2017) 057, arXiv:1705.03026; author's Paper 14 (cited, not re-proved) |
| de Sitter `[O]` ingredients | author's Paper IX |
| Axial obstruction (not this cell) | `scripts/m9_5_ec_symplectic.py`, Paper 14 |
| Newton \(1/r^2\) (not this cell) | [`../tasks/m9_2_task_details.md`](../tasks/m9_2_task_details.md), locked, no run |

## 7. What this note is not

- A `MODELS.md` icon. The metric-phenomena row stays 🚧 until a script
  measures light bending, time dilation, or \(\Lambda\) against a
  pre-registered gate
- A holographic lattice
- A proof that entanglement selects Einstein-Cartan outside AdS
- A UV completion
- A Newton-limit result (M9.2)

## 8. Adversarial self-check (this page)

This is a domain note, not a numerical claim, so there is no second-method
digit audit. The refutation targets that were walked before filing:

| Temptation | Disposition |
| --- | --- |
| Treat FGHMV as in-platform | refused; cited only |
| Treat FHHPRV as second-order Einstein-Cartan | refused; FHHPRV is Einstein; axial matching is obstructed |
| Write a de Sitter theorem | refused; `[O]` |
| Fold torsion into light / clocks / \(\Lambda\) | refused; vacuum \(T=0\) |
| Call spintronics a proof that torsion is laboratory-visible | refused; not a theorem |
| Move `MODELS.md` because the action contains Einstein+\(\Lambda\) | refused; the platform cell needs a script |

No `[O]` is moved to `[P]`.
