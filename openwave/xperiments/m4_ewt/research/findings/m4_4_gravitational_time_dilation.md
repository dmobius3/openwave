# M4 Gravitational Time Dilation from the Internal EMC Soliton Clock

## Criterion
Gravity: metric phenomena — time dilation component

## Status
✅ validated numerically

## Mechanism
A clock is a standing-wave soliton. Its period is the time for an
internal longitudinal EMC lattice wave to travel from one side of
the soliton to the opposite side and back.

In natural lattice units (\(c=1\)), the internal wave speed is

\[
v_{\text{clock}}(r) = \sqrt{\eta(r)},
\]

where

\[
\eta(r) = \frac{N_\nu(r)}{N_{\text{stat}}} = 1 - \frac{r_s}{r}.
\]

In the weak-field limit this gives the standard gravitational
redshift:

\[
\frac{\Delta f}{f} = -\Phi_N(r) = -\frac{GM}{c^2 r}.
\]

## Method
- Compute the EMC density ratio at the solar limb:
  \(\eta = 1 - r_s/R_\odot\).
- Compute the internal clock speed:
  \(v_{\text{clock}}/c = \sqrt{\eta}\).
- Predict the gravitational redshift:
  \(\Delta f/f = -\Phi_N\).
- Compare with the standard GR value.

## Result
- \(\eta = 0.999995753735\)
- \(v_{\text{clock}}/c = 0.999997876865\)
- Predicted \(\Delta f/f = -2.123132 \times 10^{-6}\)
- Target \(\Delta f/f = -2.123132 \times 10^{-6}\)
- Relative difference: \(0.0000\%\)

## Free choices
- Linear internal wave-speed response to the local density ratio is
  assumed in the weak-field limit.

## Reference

Full derivation in the Enhanced EWT manuscript, version 4.5.6:
[DOI: 10.5281/zenodo.17654657](https://doi.org/10.5281/zenodo.17654657)

Relevant section:

- “Mechanical Origin of Gravitational Redshift in the EMC Lattice”