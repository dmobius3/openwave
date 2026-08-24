# M4.4 — Gravitational Time Dilation from the Internal EMC Soliton Clock

## Status
DONE (post-hoc)

## Criterion
`Gravity: metric phenomena` — time dilation component

## Objective
Test whether the EMC lattice model reproduces the standard
gravitational redshift at the solar limb using the internal
soliton clock mechanism.

## Method
1. Define the dimensionless EMC density ratio:
   \(\eta(r) = N_\nu(r)/N_{\text{stat}} = 1 - r_s/r\).
2. Model a clock as a standing-wave soliton whose internal round-trip
   time is governed by the longitudinal EMC wave speed:
   \(v_{\text{clock}}(r) = \sqrt{\eta(r)}\).
3. Compute the fractional frequency shift at the solar limb:
   \(\Delta f/f = -\Phi_N = -GM/(c^2 R_\odot)\).
4. Compare with the standard GR value.

## Result
- Predicted \(\Delta f/f = -2.123132 \times 10^{-6}\)
- Target \(\Delta f/f = -2.123132 \times 10^{-6}\)
- Relative difference: \(0.0000\%\)

## Interpretation
This is a consistency test of the EMC clock-speed encoding, not yet
a derivation from the full BCC lattice elasticity. The model-specific
step is the relation \(v_{\text{clock}} = \sqrt{\eta}\). Its
derivation from the underlying EMC lattice remains open.

## Artifacts
- `research/scripts/m4_4_gravitational_time_dilation.py`
- `research/findings/m4_4_gravitational_time_dilation.md`

## Reference

Manuscript: Enhanced EWT, version 4.5.6,
DOI: 10.5281/zenodo.17654657.

Relevant section:

- “Mechanical Origin of Gravitational Redshift in the EMC Lattice”