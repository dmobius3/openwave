# M4.8 - Monopole Flux from Lattice Geometry with Derived Zeta

## Status
DONE (post-hoc)

## Purpose
Close the logical gap between the geometric derivation of \(G\) and
the local metric boundary condition.

Previous artifacts used the monopole amplitude \(A\) as an input.
This task derives \(A\) from the BCC lattice geometry and then uses
it as the boundary flux.

## Method

1. Compute the ideal geometric stiffness \(N_{\text{ideal}} = 8\pi^4\).
2. Compute the BCC sphere-packing fraction
   \(\eta_{\text{BCC}} = \sqrt{3}\,\pi/8\).
3. Estimate the packing impedance
   \(\zeta_{\text{est}} = (1-\eta_{\text{BCC}})/(\eta_{\text{BCC}}\,N_{\text{ideal}})\).
4. Compute the corrected stiffness
   \(N_{\zeta} = N_{\text{ideal}}(1-\zeta_{\text{est}})\).
5. Compute \(G_{\text{EWT, zeta}}\) using \(N_{\zeta}\) and the
   ideal lattice projection \(L_p^{\text{geom}} = 2/\sqrt{3}\).
6. Derive the monopole amplitude
   \(A = 2\,G_{\text{EWT, zeta}}\,M/c^2\).
7. Use \(A\) as the boundary condition for the radial Laplace
   equation.
8. Verify the profile and the asymptotic Robin condition.

## Result

| Test | \(G_{\text{EWT}}\) | Difference from CODATA |
|---|---|---|
| Pure BCC geometry | \(6.662662892091293 \times 10^{-11}\) | \(1743.57\) ppm |
| BCC + packing \(\zeta\) | \(6.674738142638409 \times 10^{-11}\) | \(65.65\) ppm |

The profile matched the analytic monopole solution to \(5.7 \times
10^{-14}\), and the Robin condition was satisfied to \(5.6 \times
10^{-14}\).

## Interpretation

The packing impedance correction, derived only from the BCC
sphere-packing fraction, reduces the raw geometric discrepancy from
about 1744 ppm to about 66 ppm.

This shows that the flux condition is not an arbitrary input but a
consequence of the lattice geometry.

## Dependency

The geometric derivation of \(G\) used in this task is available in:

- `m4_7_enhanced_ewt_geometric_consistency.py`
- `m4_7_enhanced_ewt_geometric_consistency.md`

That artifact is the reference implementation of the EWT geometric
core.

## Artifacts

- `research/scripts/m4_8_monopole_flux_from_lattice.py`
- `research/findings/m4_8_monopole_flux_from_lattice.md`

## Reference

Enhanced EWT manuscript, version 4.5.9 or later:
[DOI: 10.5281/zenodo.22110605](https://doi.org/10.5281/zenodo.22110605)

Relevant section:

- „From Microscopic EMC Displacement to the Gravitational Radius”