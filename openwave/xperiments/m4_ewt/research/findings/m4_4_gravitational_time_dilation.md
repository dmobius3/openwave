# M4.4 Gravitational Time Dilation from the Internal EMC Soliton Clock

## Criterion
Gravity: metric phenomena — time dilation component

## Status
✅ validated numerically

This validates the gravitational time dilation component only; the
light-bending component is covered separately in M4.3.

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

The predicted fractional frequency shift is computed directly as

\[
\frac{\Delta f}{f} = \frac{v_{\text{clock}}(r)}{c} - 1,
\]

not as an externally inserted metric formula.

The script evaluates the formula in SI units for numerical
convenience. The natural-unit argument in the manuscript establishes
the dimensionless structure \(v_{\text{clock}}/c = \sqrt{\eta}\),
so the conversion factors \(c\) and \(G\) cancel consistently.

## Method
- Compute the EMC density ratio at the solar limb:
  \(\eta = 1 - r_s/R_\odot\).
- Compute the internal clock speed:
  \(v_{\text{clock}}/c = \sqrt{\eta}\).
- Predict the gravitational redshift:
  \(\Delta f/f = v_{\text{clock}}/c - 1\).
- Compare with the standard GR value: \(-\Phi_N\).

## Result
- \(\eta = 0.999995753735\)
- \(v_{\text{clock}}/c = 0.999997876865\)
- Predicted \(\Delta f/f = -2.123135 \times 10^{-6}\)
- Target \(\Delta f/f = -2.123132 \times 10^{-6}\)
- Relative difference: \(0.000106\%\)

## Model assumptions (derived, not fitted)

The following are not free parameters but structural consequences of
the Enhanced EWT lattice model:

- The local EMC density ratio is
  \(\eta(r) = N_\nu(r)/N_{\text{stat}} = 1 - r_s/r\),
  the same weak-field profile used in the light-bending artifact.
- The internal clock speed is
  \(v_{\text{clock}}(r) = \sqrt{\eta(r)}\),
  following from the same lattice-strain mechanism that defines the
  external optical index.
- The redshift is computed directly from the clock-speed cycle as
  \(\Delta f/f = v_{\text{clock}}/c - 1\),
  without inserting the GR formula by hand.

No free numerical parameters are introduced.

## Reference

Full derivation in the Enhanced EWT manuscript, version 4.5.6:
[DOI: 10.5281/zenodo.17654657](https://doi.org/10.5281/zenodo.17654657)

Relevant section:

- “Mechanical Origin of Gravitational Redshift in the EMC Lattice”