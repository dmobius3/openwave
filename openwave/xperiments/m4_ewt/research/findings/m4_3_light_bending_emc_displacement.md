# M4 Light Bending from EMC Density Gradient

> This is an Enhanced EWT extension, authored by Łukasz Smoliński, as registered in
> [`_CITATIONS.md`](../../theory/_CITATIONS.md) (The Geometric Identity of Gravity and
> Dimensional Unification, v4.5.5, DOI
> [10.5281/zenodo.22042784](https://doi.org/10.5281/zenodo.22042784)).

## Criterion
Gravity: metric phenomena (light bending, with gravitational time
dilation treated as the same EMC-deformation effect; Lambda omitted)

## Status
⚠️ partial validation candidate

## Mechanism

In the Enhanced EWT framework, the speed of light is the structural
conversion factor between the spatial and temporal steps of the BCC
lattice:

\[
c \equiv \frac{\lambda_l}{t_p}.
\]

In the natural units of the lattice, \(c = 1\) and therefore
\([m] = [s]\). Consequently, a single EMC-density deformation
manifests simultaneously as:

- **light bending** — the ray follows the deformed lattice geometry
  produced by the EMC displacement field
  \(\vec{u}(r) = -\chi \nabla N_\nu(r)\),
- **gravitational time dilation** — a clock ticks more slowly
  because the geometric path required for each internal signal
  changes in the same density gradient.

Thus, within this model, a test of light bending is also a test of
the geometric mechanism underlying gravitational time dilation.
They are not independent phenomena.

## Method

- Assumed EMC density profile:
  \[
  N_\nu(r) = N_{\text{stat}} \left(1 - \frac{2r_s}{r}\right),
  \]
  where \(r_s = 2GM_{\odot}/c^2\).

- The displacement field is encoded by the scalar
  \[
  n(r) = 1 / \sqrt{1 - 2r_s/r},
  \]
  which is not an independent optical assumption but a convenient
  representation of the deformed EMC geometry.

- The bending angle is obtained from the standard ray integral in
  the variable \(u = R_{\odot}/r\):

  \[
  \Delta\theta
  = \frac{2r_s}{R_{\odot}}
    \int_0^1
    \frac{
      u\left(1 - \frac{2r_s u}{R_{\odot}}\right)^{-3/2}
    }{
      \sqrt{1-u^2}
    }
    \,du .
  \]

## Result

- Solar-limb bending angle:
  \(\Delta\theta = 1.751728\) arcsec
- Reference value (general-relativistic prediction): \(1.7517\) arcsec
- Relative difference: \(0.0016\%\)

## Free choices

- The profile \(N_\nu(r)\) is chosen as the simplest weak-field model
  consistent with the EMC push-out mechanism.
- \(n(r)\) is treated as the scalar encoding of the EMC displacement
  field; its direct derivation from BCC lattice elasticity remains
  open.

## Interpretation

The numerical result shows that once the EMC density deficit is
encoded as \(n(r) = (N_\nu/N_{\text{stat}})^{-1/2} = (1 - 2r_s/r)^{-1/2}\),
the standard ray integral reproduces the observed solar-limb bending.
The encoding, not the bending, is the model-specific step: deriving
\(n(r)\) from BCC lattice elasticity remains open (Free choices), so
this is a consistency test of the encoding, not yet a derivation of
light bending from the EMC mechanism. The same encoding is argued,
via \(c \equiv \lambda_l/t_p\), to control clock rates; gravitational
time dilation is not separately computed here.

Formally, the ray bends because of the gradient of the phase
velocity. In EWT, the physical carrier of this gradient is the
lattice deformation field \(\vec{u}(r)\), not an abstract optical
property.

## Reference

Full derivation in the Enhanced EWT manuscript, version 4.5.6:
[DOI: 10.5281/zenodo.17654657](https://doi.org/10.5281/zenodo.17654657)

Relevant sections:

- “The Two Faces of EMC Displacement: Speed and Trajectory”
- “Bridging the Vector Displacement to the Scalar Refractive Index”
- “Asymptotic Continuous Limit and Schwarzschild Equivalence”
