# M4.6 In-Platform Derivation of the Far-Field EMC Density Profile

> This is an Enhanced EWT extension, authored by Łukasz Smoliński, as registered in
> [`_CITATIONS.md`](../../theory/_CITATIONS.md) (The Geometric Identity of Gravity and
> Dimensional Unification, v4.5.8, DOI
> [10.5281/zenodo.22100322](https://doi.org/10.5281/zenodo.22100322)).

## Criterion
Gravity: local metric phenomena (foundational density profile)

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
- Max relative error of \(\delta\eta(r)\) itself:
  \(2.5 \times 10^{-8}\) (RK4 truncation at \(\Delta r = 0.025\,R_\odot\); the
  \(\eta\)-relative figure above is this number diluted by \(r_s/R_\odot = 4.2 \times 10^{-6}\))

## Model assumptions (derived in the manuscript, not fitted here)

- Continuum limit \(r \gg \lambda_l\).
- Linearised weak-field response.
- Spherical symmetry.

No free numerical parameters are used.

## Interpretation

The boundary-value solution reproduces the analytic monopole profile: the
\(1/r\) form is the harmonic solution of the radial Laplace equation with
the constant killed by the Robin condition, and the RK4 solve confirms it
to \(2.5 \times 10^{-8}\) of \(\delta\eta\).

What the artifact does and does not derive: the \(1/r\) shape comes out of
the equation; the amplitude does not. \(r_s = 2GM/c^2\) enters through the
Gauss flux condition at the solar surface,
\(d\delta\eta/dr\,(R_\odot) = r_s/R_\odot^2\), and comes back out as
\(A\). The derivation of that flux (and of the encoding of \(\eta\) into
the refractive index and the clock speed used in M4.3 to M4.5) from the
lattice elasticity lives in the manuscript, not in-platform. This is a
consistency check of the weak-field profile's form, the same standing as
M4.3 to M4.5.

## Reference

Full derivation in the Enhanced EWT manuscript, version 4.5.8:
[DOI: 10.5281/zenodo.22100322](https://doi.org/10.5281/zenodo.22100322)

Relevant sections:

- "Macroscopic EMC Density Profile of a Star"
- "Derivation of the Far-Field EMC Density Profile"
