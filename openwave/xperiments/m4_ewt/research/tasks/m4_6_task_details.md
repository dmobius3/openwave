# M4.6 - In-Platform Derivation of the Far-Field EMC Density Profile

## Status
DONE (post-hoc)

## Criterion
`Gravity: local metric phenomena` - foundational density profile

## Objective
Numerically solve the far-field EMC lattice equilibrium equations and
verify that the resulting normalised density profile matches the
analytic monopole form used in M4.3-M4.5.

## Method
1. Model the BCC lattice as an isotropic elastic continuum for
   \(r \gg \lambda_l\).
2. Solve the radial Laplace equation
   \[
   \frac{d}{dr}\left(r^2 \frac{d\,\delta\eta}{dr}\right) = 0
   \]
   outside the solar surface.
3. Apply the Gauss strain flux condition at \(r = R_\odot\).
4. Apply the asymptotic Robin boundary condition at
   \(r = R_{\text{outer}}\).
5. Determine the monopole amplitude \(A\) numerically.
6. Compare the resulting \(\eta(r)=1+\delta\eta(r)\) with
   \(\eta(r)=1-r_s/r\).

## Result
- \(A_{\text{numerical}} = 2.954126555055404 \times 10^{3}\ \text{m}\)
- \(r_s = 2.954126555055405 \times 10^{3}\ \text{m}\)
- Max relative error of \(\eta(r)\): \(5.7 \times 10^{-14}\)

## Interpretation
The numerical solution reproduces the analytic monopole profile to
machine precision. The profile is therefore not an arbitrary fitting
ansatz but the far-field equilibrium solution of the EMC lattice.

## Artifacts
- `research/scripts/m4_6_emc_far_field_profile.py`
- `research/findings/m4_6_emc_far_field_profile.md`

## Reference
Enhanced EWT manuscript, version 4.5.8:
[DOI: 10.5281/zenodo.22100322](https://doi.org/10.5281/zenodo.22100322)

Relevant sections:

- "Macroscopic EMC Density Profile of a Star"
- "Derivation of the Far-Field EMC Density Profile"