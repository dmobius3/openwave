# M4.9 Ab-Initio Emergent Encoding from Lattice Dynamics

## Criterion
Gravity: local metric phenomena — foundational encoding derivation

## What was computed

A one-dimensional spring-mass lattice was simulated. The local spring
stiffness was not assumed. It was derived from a microscopic pair
potential and then used to measure:

1. the wave-speed exponent \(\beta\) in \(v \propto \eta^{\beta}\),
2. the oscillator-frequency exponent \(\gamma\) in \(f \propto \eta^{\gamma}\).

The two resulting exponents correspond to the encodings used in
M4.3–M4.5:

- \(n_\gamma \propto \eta^{-1/2}\),
- \(v_{\text{clock}} \propto \sqrt{\eta}\).

## Emergence chain

### 1. Microscopic potential to stiffness

The lattice stiffness \(k(\eta)\) is derived from a logarithmic
repulsive pair potential,

\[
V(r) = -V_0 \ln r,
\]

evaluated at the equilibrium spacing \(a = 1/\eta\):

\[
k(\eta)
=
\left.
\frac{d^2V}{dr^2}
\right|_{r=1/\eta}
=
V_0\,\eta^2 .
\]

Thus \(k \propto \eta^2\) is not an input. It is a consequence of the
microscopic potential.

### 2. Stiffness and mass to wave speed

The mass per lattice site is proportional to the local EMC density,

\[
m \propto \eta .
\]

The speed of sound in the lattice is then

\[
v(\eta)
=
\sqrt{\frac{k(\eta)}{m(\eta)}}
=
\sqrt{\frac{\eta^2}{\eta}}
=
\sqrt{\eta} .
\]

This is the same clock-speed encoding used in M4.4.

The corresponding effective refractive index is

\[
n_\gamma(\eta)
=
\frac{c_0}{v(\eta)}
\propto
\eta^{-1/2} .
\]

This is the same index used in M4.3 and M4.5.

### 3. Oscillator frequency

A local oscillator in a density well has the same harmonic frequency
scale,

\[
f \propto \sqrt{\frac{k}{m}} \propto \sqrt{\eta} .
\]

## Measured results

### Stiffness exponent

| \(\eta\) | \(k\) |
|---|---|
| 0.4 | 0.160001 |
| 0.5 | 0.249999 |
| 0.6 | 0.360000 |
| 0.7 | 0.490000 |
| 0.8 | 0.640000 |
| 0.9 | 0.810000 |

Fitted exponent:

\[
\alpha = 1.999997
\]

Expected from the pair potential: \(2.0\).

### Wave-speed exponent

| \(\eta\) | \(v\) |
|---|---|
| 0.4 | 0.630149 |
| 0.5 | 0.704725 |
| 0.6 | 0.773320 |
| 0.7 | 0.833313 |
| 0.8 | 0.893784 |
| 0.9 | 0.947254 |

Fitted exponent:

\[
\beta = 0.502855
\]

Expected: \(0.5\).

### Oscillator-frequency exponent

| \(\eta_{\text{core}}\) | \(f\) |
|---|---|
| 0.4 | 0.040218 |
| 0.5 | 0.044942 |
| 0.6 | 0.049221 |
| 0.7 | 0.053151 |
| 0.8 | 0.056811 |
| 0.9 | 0.060249 |

Fitted exponent:

\[
\gamma = 0.498444
\]

Expected: \(0.5\).

## Model assumptions

The only microscopic inputs are:

- logarithmic pair potential \(V(r) = -V_0 \ln r\),
- mass density \(m \propto \eta\).

Both follow from the BCC lattice description used in the manuscript.
The exponents were not assumed and were not fitted to the target
values.

## Interpretation

The two Schwarzschild encodings used in M4.3–M4.5 now have an
in-platform lattice-dynamical derivation.

The exponents \(1/2\) emerge from:

1. the microscopic pair potential,
2. the mass-density scaling,
3. the wave equation on the resulting spring-mass lattice.

No independent optical or clock hypothesis is required.

## Reference

Enhanced EWT manuscript, version 4.5.11:
[DOI: 10.5281/zenodo.22133680](https://doi.org/10.5281/zenodo.22133680)

Relevant section:

- „Emergent Metric Encodings from Lattice Dynamics”