# M4.10 - Newtonian Force from EMC Push-Out Pressure

## Status
DONE (post-hoc)

## Criterion
`Gravity: Newton limit (GEM)`

## Objective
Show that the attractive \(1/r^2\) force law follows from the EMC
push-out mechanism using exact-domain numerical field integration.

## Method

1. Define monopole density deficits \(\delta\eta = -A/r\).

2. Use amplitude \(A = 2G_{\text{EWT, geo}}M/c^2\).

3. Evaluate the gradient-overlap energy
   \(\int \nabla\eta_1\cdot\nabla\eta_2\,dV\)
   over the entire infinite domain.

4. Use the coordinate mapping \(r=R/(1-t)\) to avoid truncation
   errors at large radius.

5. Compute the force by central numerical differentiation.

6. Convert the geometric result to physical force with \(K_{\text{emc}}\).

7. Compare with Newton’s law.

## Result

- \(F_{\text{EMC}} = 3.542516096914 \times 10^{22}\ \text{N}\)
- \(F_{\text{Newton}} = 3.542516523099 \times 10^{22}\ \text{N}\)
- Relative difference: \(1.203 \times 10^{-5}\%\)

## Interpretation

The Newtonian force law emerges from the EMC pressure mechanism.
The \(4\pi\) factor comes from the angular integral, and the radial
integral is performed numerically over the full domain.

## Artifacts

- `research/scripts/m4_10_newtonian_force_emc.py`
- `research/findings/m4_10_newtonian_force_emc.md`

## Reference

Enhanced EWT manuscript, version 4.5.12 or later:
[DOI: 10.5281/zenodo.22140646](https://doi.org/10.5281/zenodo.22140646)