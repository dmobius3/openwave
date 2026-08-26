# M4.8 Monopole Flux from Lattice Geometry with Derived Zeta

## Purpose

This artifact closes the remaining logical gap in the M4 local
metric chain.

Previous artifacts M4.3–M4.6 established that the far-field profile

\[
\eta(r) = 1 - \frac{A}{r}
\]

reproduces light bending, gravitational redshift, and Shapiro delay
when the monopole amplitude \(A\) is known.

This artifact derives \(A\) from the BCC lattice geometry itself,
without using the measured gravitational constant as an input.

## What was computed

Two tests were performed:

### Test 1 — Pure BCC geometry

The ideal geometric stiffness was used:

\[
N_{\text{ideal}} = 8\pi^4 .
\]

This gave:

- \(G_{\text{EWT, pure}} = 6.662662892091293 \times 10^{-11}\)
- Amplitude \(A_{\text{pure}} = 2.948975829211922 \times 10^{3}\) m
- Difference from CODATA: \(1743.57\) ppm

This is the raw scale of the problem. The pure geometric limit is
close, but not sufficient.

### Test 2 — BCC geometry corrected by the packing impedance

The estimated packing impedance was used:

\[
\zeta_{\text{est}}
=
\frac{1-\eta_{\text{BCC}}}{\eta_{\text{BCC}}\,N_{\text{ideal}}}
\]

where

\[
\eta_{\text{BCC}} = \frac{\sqrt{3}\,\pi}{8}
\]

is the BCC sphere packing fraction.

This gave:

- \(N_{\zeta} = N_{\text{ideal}}(1-\zeta_{\text{est}})\)
- \(G_{\text{EWT, zeta}} = 6.674738142638409 \times 10^{-11}\)
- Amplitude \(A_{\zeta} = 2.954320482329130 \times 10^{3}\) m
- Difference from CODATA: \(65.65\) ppm

This is the key result. The packing impedance correction, derived
solely from the BCC sphere-packing fraction, brings the amplitude
to within about \(66\) ppm of the measured value.

No empirical value of \(G\) and no calibration to the measured
fine-structure constant were used.

## Role in the platform

This artifact demonstrates that the flux condition used in M4.6 is
not an independent boundary input.

The monopole amplitude \(A\) is first derived from BCC lattice
parameters and then used as the boundary condition for the Laplace
equation.

This completes the inversion required for the local metric
validations.

## Dependency on M4.7

The full geometric derivation of \(G_{\text{EWT}}\) is implemented
in:

- `m4_7_enhanced_ewt_geometric_consistency.py`

M4.8 does not re-derive the gravitational identity from scratch.
It uses the same BCC lattice parameters and the same geometric
chain, then extends it to the far-field monopole flux condition.

This keeps the artifact focused and traceable.

## Relation to criteria

This artifact is a foundational contribution to:

- `Gravity: local metric phenomena`
- `Fundamental constants: gravitational constant G`

It is not itself a pass/fail validation against a single row, but it
supplies the missing link between the geometric derivation of \(G\)
and the local metric observables.

## Result summary

| Test | \(G_{\text{EWT}}\) | \(A\) | Difference from CODATA |
|---|---|---|---|
| Pure BCC geometry | \(6.662662892091293 \times 10^{-11}\) | \(2.948975829211922 \times 10^{3}\) m | \(1743.57\) ppm |
| BCC + packing \(\zeta\) | \(6.674738142638409 \times 10^{-11}\) | \(2.954320482329130 \times 10^{3}\) m | \(65.65\) ppm |

The radial Laplace equation was solved using the corrected amplitude.
The profile matched the analytic monopole form to a maximum relative
error of \(5.7 \times 10^{-14}\), and the asymptotic Robin condition
was satisfied to \(5.6 \times 10^{-14}\).

## Model assumptions

The model follows the Enhanced EWT manuscript, version 4.5.9 or later.

The only lattice parameters used in the derivation chain are:

- the BCC coordination number,
- the ideal sphere-packing fraction,
- the geometric lattice projection factor \(L_p^{\text{geom}}\).

No free numerical parameters were introduced.

## Reference

Full derivation in the Enhanced EWT manuscript, version 4.5.9 or later:
[DOI: 10.5281/zenodo.22110605](https://doi.org/10.5281/zenodo.22110605)

Relevant section:

- „From Microscopic EMC Displacement to the Gravitational Radius”