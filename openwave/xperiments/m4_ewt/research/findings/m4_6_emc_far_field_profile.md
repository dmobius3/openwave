# M4.6 In-Platform Derivation of the Far-Field EMC Density Profile

## Criterion
Gravity: local metric phenomena — foundational density profile

## What was computed

The weak-field EMC density profile around a spherical mass was generated
by numerically solving the radial Laplace equation with the monopole
strain flux condition at the solar surface and the asymptotic Robin
boundary condition at the outer boundary.

The resulting normalised density was compared with the analytic form

\[
\eta(r) = 1 - \frac{r_s}{r}
\]

used in M4.3–M4.5.

## Result

- Monopole amplitude:
  \(A_{\text{numerical}} = 2.954126555055404 \times 10^{3}\ \text{m}\)
- Expected gravitational radius:
  \(r_s = 2.954126555055405 \times 10^{3}\ \text{m}\)
- Max relative error of \(\eta(r)\):
  \(5.7 \times 10^{-14}\)

## Model assumptions (derived in the manuscript, not fitted here)

- Continuum limit \(r \gg \lambda_l\).
- Linearised weak-field response.
- Spherical symmetry.

No free numerical parameters are used.

## Interpretation

The numerical boundary-value solution reproduces the analytic monopole
profile to machine precision. This closes the gap between the assumed
encoding used in M4.3–M4.5 and the far-field lattice equilibrium: the
profile is generated from the equation, not inserted by hand.

## Reference

Full derivation in the Enhanced EWT manuscript, version 4.5.8:
[DOI: 10.5281/zenodo.22100322](https://doi.org/10.5281/zenodo.22100322)

Relevant sections:

- "Macroscopic EMC Density Profile of a Star"
- "Derivation of the Far-Field EMC Density Profile"