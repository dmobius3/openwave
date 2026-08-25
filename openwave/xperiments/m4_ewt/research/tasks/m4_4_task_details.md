# M4.4 - Gravitational Time Dilation from the Internal EMC Soliton Clock

## Status
DONE (post-hoc)

## Criterion
`Gravity: local metric phenomena` - time dilation component

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
4. Compare with the exact GR value for a static clock,
   \(\sqrt{1 - 2\Phi_N} - 1\) (identity check), and print the
   first-order value \(-\Phi_N\) as a labelled sanity line.

## Result
- Predicted \(\Delta f/f = -2.123135 \times 10^{-6}\)
- Exact GR \(\Delta f/f = -2.123135 \times 10^{-6}\); difference \(0\)
- First-order reference \(-\Phi_N = -2.123132 \times 10^{-6}\);
  relative difference \(0.000106\% = x/4\), \(x = r_s/R_\odot\)

## Interpretation
This is a consistency test of the EMC clock-speed encoding. Once
\(\eta = 1 - r_s/r\) is granted, \(\sqrt{\eta} - 1\) is identically
the exact Schwarzschild redshift factor, so the model-vs-GR difference
is zero. The \(0.000106\%\) figure is the truncation error of the
first-order reference \(-\Phi_N\) (exactly \(x/4\)), not a model
residue. The derivation of \(\eta(r)\) and \(v_{\text{clock}} =
\sqrt{\eta}\) from lattice elasticity is given in the manuscript
(Reference); it is not recomputed in-platform.

## Artifacts
- `research/scripts/m4_4_gravitational_time_dilation.py`
- `research/findings/m4_4_gravitational_time_dilation.md`

## Reference
Enhanced EWT manuscript, version 4.5.6:
[DOI: 10.5281/zenodo.22086668](https://doi.org/10.5281/zenodo.22086668) (version DOI; concept DOI 10.5281/zenodo.17654657)

Relevant section:

- "Mechanical Origin of Gravitational Redshift in the EMC Lattice"
