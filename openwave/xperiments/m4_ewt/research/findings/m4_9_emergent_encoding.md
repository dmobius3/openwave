# M4.9 Emergent Encoding from Lattice Dynamics

## Criterion
Gravity: local metric phenomena — foundational encoding derivation

## What was computed

A one-dimensional spring-mass lattice was simulated. The lattice
spacing was treated as a physical length scale \(a(\eta)=1/\eta\), and
the per-site mass was kept fixed.

The local stiffness \(k(\eta)\) was derived from a microscopic
one-dimensional pair potential

\[
V(r)=\frac{V_0}{r}
\]

evaluated at the equilibrium spacing.

The resulting wave speed was measured in physical length units, not in
site-index units.

## Emergence chain

### 1. Microscopic potential to stiffness

The stiffness is the second derivative of the pair potential at the
equilibrium spacing:

\[
k(\eta)
=
\left.
\frac{d^2V}{dr^2}
\right|_{r=1/\eta}
\propto \eta^3
\]

### 2. Stiffness and mass to physical wave speed

With fixed mass per lattice site \(m_0\), the long-wavelength speed is

\[
v_{\text{phys}}
=
a(\eta)\sqrt{\frac{k(\eta)}{m_0}}
\]

which scales as

\[
v_{\text{phys}} \propto \eta^{+1/2}
\]

### 3. Refractive index

The corresponding effective refractive index is

\[
n_\gamma(\eta)
=
\frac{c_0}{v_{\text{phys}}}
\propto
\eta^{-1/2}
\]

This is the same encoding used in M4.3 and M4.5.

### 4. Clock frequency

For a standing-wave clock of fixed physical length \(L\),

\[
f \propto \frac{v_{\text{phys}}}{L}
\propto \eta^{+1/2}
\]

This is the same encoding used in M4.4.

## Measured results

### Stiffness exponent \(\alpha\)

| \(\eta\) | \(k\) |
|---|---|
| 0.4 | 0.128000 |
| 0.5 | 0.250000 |
| 0.6 | 0.432000 |
| 0.7 | 0.686001 |
| 0.8 | 1.024000 |
| 0.9 | 1.458000 |

Fitted exponent:

\[
\alpha = 3.000002
\]

Expected: \(3.0\). This is the analytic derivative, not an independent
measurement.

### Wave-speed exponent \(\beta\)

| \(\eta\) | \(v_{\text{phys}}\) |
|---|---|
| 0.4 | 0.889680 |
| 0.5 | 0.998259 |
| 0.6 | 1.091577 |
| 0.7 | 1.180412 |
| 0.8 | 1.260294 |
| 0.9 | 1.337968 |

Fitted exponent:

\[
\beta = 0.501714
\]

Expected: \(0.5\).

### Oscillator-frequency exponent \(\gamma\)

The oscillator frequency follows from \(f = v_{\text{phys}}/L\).
Therefore

\[
\gamma = \beta = 0.501714
\]

Expected: \(0.5\).

Phase 3 is a consistency check on phase 2, not an independent
simulation.

## Model assumptions

- One-dimensional reduction along a high-symmetry axis of the BCC lattice.
- Fixed mass per lattice site \(m_0\).
- Pair potential \(V(r)=V_0/r\).
- Lattice spacing \(a(\eta)=1/\eta\).

These are the corrected assumptions from the reviewer’s B1/B2 findings.

## Interpretation

The two Schwarzschild encodings used in M4.3–M4.5 now follow from
microscopic lattice dynamics:

- \(v_{\text{clock}} \propto \sqrt{\eta}\)
- \(n_\gamma \propto \eta^{-1/2}\)

The exponents were not assumed. They emerge from the corrected
microscopic model.

## Reference

Enhanced EWT manuscript, version 4.5.12:
[DOI: 10.5281/zenodo.22140646](https://doi.org/10.5281/zenodo.22140646)

Relevant section:

- „Emergent Metric Encodings from Lattice Dynamics”