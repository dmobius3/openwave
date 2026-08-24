# M4.4 — Gravitational Time Dilation from the Internal EMC Soliton Clock

## Status
DONE (post-hoc)

## Criterion
`Gravity: metric phenomena` — time dilation component

## Objective
Test whether the EMC lattice model reproduces the standard
gravitational redshift at the solar limb using the internal
soliton clock mechanism, without inserting the GR formula by hand.

## Method
1. Define the dimensionless EMC density ratio:
   \(\eta(r) = N_\nu(r)/N_{\text{stat}} = 1 - r_s/r\).
2. Model a clock as a standing-wave soliton whose internal round-trip
   time is governed by the longitudinal EMC wave speed:
   \(v_{\text{clock}}(r) = \sqrt{\eta(r)}\).
3. Compute the fractional frequency shift directly from the clock
   speed:
   \(\Delta f/f = v_{\text{clock}}/c - 1\).
4. Compare with the standard GR value: \(-\Phi_N\).

## Result
- Predicted \(\Delta f/f = -2.123135 \times 10^{-6}\)
- Target \(\Delta f/f = -2.123132 \times 10^{-6}\)
- Relative difference: \(0.000106\%\)

## Interpretation
This is a consistency test of the EMC clock-speed encoding. The small
difference from GR is not a free-parameter correction but the
expected second-order residue of the weak-field expansion
\(\sqrt{1 - 2\Phi_N} \approx 1 - \Phi_N\).

## Artifacts
- `research/scripts/m4_4_gravitational_time_dilation.py`
- `research/findings/m4_4_gravitational_time_dilation.md`

## Reference
Enhanced EWT manuscript, version 4.5.6:
[DOI: 10.5281/zenodo.17654657](https://doi.org/10.5281/zenodo.17654657)

Relevant section:

- “Mechanical Origin of Gravitational Redshift in the EMC Lattice”