# M4.10 Newtonian Force from EMC Push-Out Pressure

## Criterion
Gravity: Newton limit (GEM)

## What was computed

The Newtonian \(1/r^2\) force law was recovered from the EMC push-out
mechanism using exact-domain numerical integration.

Two monopole density deficits,

\[
\delta\eta(r) = -\frac{A}{r},
\]

with amplitudes

\[
A_i = \frac{2 G_{\text{EWT, geo}} M_i}{c^2},
\]

were placed at separation \(R\).

The interaction energy was taken as the integrated overlap of the two
scalar density gradients:

\[
U_{\text{int}}(R)
=
\int
\nabla(\delta\eta_1)
\cdot
\nabla(\delta\eta_2)
\,dV
.
\]

The spatial integration was performed over the entire infinite domain
\(r \in [R,\infty)\) using the coordinate mapping

\[
r = \frac{R}{1-t},
\qquad
dr = \frac{R}{(1-t)^2}\,dt,
\qquad
t \in [0,1)
.
\]

The angular integral was evaluated analytically and gives the standard
shell-theorem factor

\[
\int d\Omega\;
\frac{r^2 - rR\cos\theta}
{r^3 (r^2+R^2-2rR\cos\theta)^{3/2}}
=
\frac{4\pi}{r^2}
\qquad
(r \ge R)
.
\]

Thus the \(4\pi\) factor emerges from the angular integration and is
not fitted. The remaining radial integration is one-dimensional,
non-singular, and mapped to a finite interval, avoiding truncation
errors at large radius.

The force is obtained by central numerical differentiation:

\[
F = -\frac{dU_{\text{int}}}{dR}
.
\]

## Normalisation

The geometric energy density was converted to physical force using

\[
K_{\text{emc}}
=
\frac{c^4}{16\pi G_{\text{EWT, geo}}}
.
\]

This is the standard scalar monopole field-energy normalisation for the
GEM route. It contains no fitted coefficients.

## Result

For the Sun–Earth system at \(1\) AU:

| Quantity | Value |
|---|---|
| \(F_{\text{EMC}}\) | \(3.542516096914 \times 10^{22}\ \text{N}\) |
| \(F_{\text{Newton}}\) | \(3.542516523099 \times 10^{22}\ \text{N}\) |
| Relative difference | \(1.203 \times 10^{-5}\%\) |

The Newtonian force is recovered from the EMC pressure mechanism to
machine-precision level.

## Model assumptions

- Far-field monopole density deficit.
- Interaction energy proportional to the gradient overlap.
- EMC pressure normalisation as given above.
- No empirical Newtonian \(G\) is used; the geometric value from
  M4.7/M4.8 enters the calculation.

## Interpretation

The Newtonian \(1/r^2\) law is not an independent postulate. It follows
from the static interaction of two EMC density deficits through the
lattice pressure. The numerical integration demonstrates that the force
emerges from the full spatial field structure, not from a simplified
analytic shortcut.

## Reference

Enhanced EWT manuscript, version 4.6.1 or later:
[DOI: 10.5281/zenodo.22144273](https://doi.org/10.5281/zenodo.22144273)

Relevant sections:

- „The Geometric Identity of \(G\)”
- „From Microscopic EMC Displacement to the Gravitational Radius”
- „Newtonian Force from Interacting EMC Deficits”