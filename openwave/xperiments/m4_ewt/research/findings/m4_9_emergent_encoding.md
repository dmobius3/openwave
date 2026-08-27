# M4.9 Emergent Encoding from Lattice Dynamics

## Criterion
Gravity: local metric phenomena — foundational encoding derivation

## What was computed

A one-dimensional spring-mass lattice was simulated. The microscopic
parameters were taken from the BCC continuum picture:

- mass per lattice site: \(m \propto \eta\)
- spring stiffness: \(k \propto \eta^2\)

No assumption was made about the final form of the refractive index
or the clock rate.

Two properties were measured:

1. The speed of a propagating pulse as a function of uniform density
   \(\eta\).
2. The frequency of an oscillator in a wide density well as a function
   of the core density \(\eta_{\text{core}}\).

## Method

For the wave-speed test:

- A Gaussian pulse was initialised on a uniform lattice.
- The pulse position was tracked with sub-grid precision using a
  squared centroid.
- A unidirectional initialisation prevented pulse splitting.
- The speed was measured as the centroid displacement over time.

For the oscillator test:

- A static Gaussian density well was used.
- A compact carrier wave packet was initialised in the well centre.
- Exact zero-crossing times were obtained by linear interpolation.
- The frequency was extracted from the first clean half-periods.

The exponents were obtained by linear regression in log-log space.

## Result

Wave-speed test:

| \(\eta\) | measured speed |
|---|---|
| 0.4 | 0.630147 |
| 0.5 | 0.704726 |
| 0.6 | 0.773320 |
| 0.7 | 0.833313 |
| 0.8 | 0.893784 |
| 0.9 | 0.947254 |

Fitted exponent:

\[
\beta = 0.502856
\]

Oscillator test:

| \(\eta_{\text{core}}\) | measured frequency |
|---|---|
| 0.4 | 0.040220 |
| 0.5 | 0.044945 |
| 0.6 | 0.049220 |
| 0.7 | 0.053151 |
| 0.8 | 0.056811 |
| 0.9 | 0.060249 |

Fitted exponent:

\[
\gamma = 0.498371
\]

Both measured exponents are close to the expected value \(0.5\).

## Model assumptions

The microscopic scalings are:

\[
m \propto \eta,
\qquad
k \propto \eta^2
\]

They follow from the continuum BCC lattice picture and are the same
scalings used in the manuscript. They are not fitted to the target
exponent.

## Interpretation

The two Schwarzschild encodings used in M4.3–M4.5,

\[
n_\gamma \propto \eta^{-1/2},
\qquad
v_{\text{clock}} \propto \sqrt{\eta},
\]

are recovered as emergent properties of the lattice dynamics.

The exponents were not assumed. They were measured from the
microscopic equations.

## Reference

Enhanced EWT manuscript, version 4.5.9 or later:
[DOI: 10.5281/zenodo.22110605](https://doi.org/10.5281/zenodo.22110605)

Relevant sections:

- The Two Faces of EMC Displacement
- Mechanical Origin of Gravitational Redshift in the EMC Lattice