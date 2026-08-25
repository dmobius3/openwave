# M4.7 Enhanced EWT Geometric Consistency

## Purpose

This artifact is a complete Python port of the Enhanced EWT
geometric consistency suite originally written in Scilab.

It is not tied to a single matrix criterion. Its role is to:

- introduce the Enhanced EWT framework to OpenWave as a runnable
  calculation,
- provide a single entry point to the model’s core logic,
- serve as a shared geometric foundation for future validation
  tasks.

## What was computed

A complete Python port of the Enhanced EWT numerical suite was
executed. The suite is a single deterministic chain in which one
geometric parameter,

\[
\epsilon_M = \frac{1}{N_{\text{final}}\pi^3}
\]

with

\[
N_{\text{final}} \to 8\pi^4
\]

drives all physical outputs.

The script reproduces the full Scilab suite and demonstrates that:

- \(G\),
- \(\alpha^{-1}\),
- \(a_e\),
- \(a_\mu\),
- \(a_\tau\),
- particle mass estimates,
- Weinberg and Cabibbo mixing,
- atomic scales,
- the neutrino radius \(r_\nu\),
- and the geometric fixed point \(g_v\)

are all consequences of the same BCC lattice geometry.

## Relation to criteria

This artifact supports several future or already existing tasks.
It can be reused to extract isolated calculations for:

- static geometric origin of \(G\),
- fine-structure constant derivation,
- lepton AMM predictions,
- local metric phenomena,
- mass hierarchy scans.

It is not itself a pass/fail validation against a single row.

## Structure of the artifact

### Part I — Gravitational constant from geometry

The script derives \(G_{\text{Base}}\) from the electron soliton
parameters and then applies the dimensionless push-out attenuation:

\[
G_{\text{EWT}} =
\frac{c^2 r_e}{m_e}
\frac{1}{A_\pi}
\left(\frac{1}{N_{\text{final}}A_\pi}\right)^3
\frac{1}{K_{WC}\sqrt{N_{\nu,\text{eff}}}}
\]

It also evaluates the pure geometric lattice variant

\[
L_p^{\text{geom}} = \frac{2}{\sqrt{3}}
\]

which alone yields \(G\) to about \(4.78\) ppm.

This is the geometric origin of the GEM criterion.

### Part II — Neutrino radius and the \(1:10^{10}\) hierarchy

The ratio \(r_e/r_\nu \approx 100\) is verified. The implied energy
density ratio is

\[
\left(\frac{r_e}{r_\nu}\right)^5 \approx 10^{10}
\]

which anchors the electron–neutrino hierarchy.

### Part III — Electron anomalous magnetic moment

The geometric base anomaly is computed from

\[
a_{\text{Base}} =
\frac{\alpha}{2\pi}
\left(1 - \epsilon_M \pi^3\right)
\]

and compared with CODATA.

### Part IV — Fine-structure constant

The inverse fine-structure constant is derived as

\[
\alpha^{-1} =
\left(4\pi^3 + \pi^2 + \pi\right) - \epsilon_M
\]

The same modulator \(\epsilon_M\) that shapes \(G\) also shapes
\(\alpha\).

### Part V — Recursive lepton hierarchy

The nodal growth law

\[
K_n = K_{n-1} + \text{round}(10^{n-1} 2\pi^2)
\]

generates the three generations:

\[
K_e = 10,\quad K_\mu = 207,\quad K_\tau = 2181
\]

The same \(\epsilon_M\) drives the anomalous moments of all three
leptons. The full \(a_\mu\) and \(a_\tau\) predictions are compared
with CODATA and PDG values.

Crucially, the AMM predictions do not use the lepton masses as input.
They are purely geometric.

### Part VI — Mass prediction modes

The script computes masses in spherical, orbital, and meson modes.
It reproduces the electron, muon, tau, quark, and heavy boson mass
scales.

The orbital mode for \(\mu\) and \(\tau\) uses calibrated amplitude
factors derived from the electron rest mass. This is a mass-sector
calibration, not an AMM input, and it is kept separate from the
geometric AMM derivation.

### Part VII — Mixing angles and dimensional hierarchy

The geometric ladder is made explicit:

- \(\pi^6\) for volumetric bosonic coupling,
- \(\pi^5\) for fermionic surface mixing,
- \(C_{\text{local}} = \epsilon_M/(2\sqrt{2})\) as the common seed.

This yields the Weinberg and Cabibbo angles.

### Part VIII — Statutory radius and decadic resonance

The neutrino radius is derived from Planck charge and Euler’s number.
The resulting \(r_e/r_\nu\) ratio is verified independently.

### Part IX — Heavy boson radii

Masses obtained in spherical mode are used to predict the radii of
the Z and Higgs solitons.

### Part X — Zero-parameter proof

The stiffness constant is reduced to pure topology:

\[
N_{\text{geometric}} = 8\pi^4
\]

and the deficit becomes

\[
\epsilon_M = \frac{1}{8\pi^7}
\]

This gives a parameter-free derivation of \(\alpha^{-1}\).

### Part XI — Unified geometric AMM identity

The electron anomaly is written in the closed form:

\[
a_e = \frac{8\pi^4 - 1}{2\pi\left(8\pi^4 A_\pi - \pi^{-3}\right)}
\]

showing that \(a_e\) is a static geometric property.

### Part XII — Atomic scales

The Rydberg constant, Bohr radius, and Compton wavelength are derived
by combining the geometrically derived \(\alpha\) with the geometric
electron radius \(r_e\), using the standard atomic relations.

These are therefore not independent first-principles derivations, but
direct consequences of the geometric \(\alpha\) and \(r_\nu\)
inherited from the earlier parts.

### Part XIII — Full particle scan

A large scan over baryons and mesons is performed. It identifies
near-integer \(K\) resonances without parameter adjustment.

### Part XIV — Geometric derivation of \(r_\nu\)

The scaling factor \(K = r_\nu/q_P\) is decomposed into three parts:

- static lattice projection,
- dynamic wave expansion,
- discrete lattice impedance.

The resulting quadratic equation for \(g_v\) has one physical root:

\[
g_v = 0.983594444559461
\]

confirming that \(g_v\) is not fitted but fixed by geometry.

## Result summary

Key exact matches between Scilab and Python:

| Quantity | Value |
|---|---|
| \(G_{\text{EWT, unified}}\) | \(6.674305000000013 \times 10^{-11}\) |
| \(G_{\text{EWT, geo}}\) | \(6.674336927110799 \times 10^{-11}\) |
| \(\alpha^{-1}_{\text{EWT}}\) | \(137.036262365010\) |
| \(a_e^{\text{EWT}}\) | \(0.001159918486472\) |
| \(a_\mu^{\text{EWT}}\) | \(0.00116620603122654\) |
| \(a_\tau^{\text{EWT}}\) | \(0.00117684332510945\) |
| \(r_\nu\) | \(2.817932844758866 \times 10^{-17}\) m |
| \(g_v\) predicted | \(0.983594444559461\) |

## Model assumptions

The model follows the Enhanced EWT manuscript, version 4.5.8 or later.
The gravitational, fine-structure, and lepton AMM sectors are
parameter-free.

The muon and tau masses in orbital mode use calibrated amplitude
factors; these are mass-sector inputs and do not enter the geometric
AMM calculations.

## Attribution

This numerical suite incorporates foundational Energy Wave Theory
calculations originally developed by Jeff Yee. The Enhanced EWT
extension, geometric formalization, and the BCC lattice interpretation
are the contribution of £ukasz Smoliñski. Details are available in the
manuscript.

## Source Scilab script

The Python port is based on the original Enhanced EWT Scilab script
archived at:

[DOI: 10.5281/zenodo.21503571](https://doi.org/10.5281/zenodo.21503571)

## Reference

Full derivation in the Enhanced EWT manuscript, version 4.5.8 or later:
[DOI: 10.5281/zenodo.22100322](https://doi.org/10.5281/zenodo.22100322)

Relevant sections:

- The Geometric Identity of \(G\)
- Geometric Equation of the Fine-Structure Constant
- Recursive Lepton Hierarchy
- Numerical Verification