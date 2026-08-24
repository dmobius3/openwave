# M4.3 — Static solar light bending from the EMC density encoding

> Contributed by Łukasz Smoliński in
> [PR #462](https://github.com/openwave-labs/openwave/pull/462); this document's method
> text was provided by the contributor in the PR thread and transcribed at merge, with
> the relative-difference figure corrected to match the 1.7517 reference (announced in
> the same thread).

## Status

DONE (post-hoc: the criteria below were recorded after the numbers existed, which the
review noted and accepted as honest filing for contributed work)

## Criterion

`Gravity: metric phenomena` — light bending component

## Objective

Test whether the EMC density-deficit encoding

\[
n(r) = \left(\frac{N_\nu(r)}{N_{\nu,\text{stat}}}\right)^{-1/2}
\]

reproduces the observed solar-limb light bending when inserted into
the standard ray integral.

## Method

1. Assume the simplest weak-field EMC density profile outside the Sun:

   \[
   N_\nu(r) = N_{\text{stat}}\left(1 - \frac{2r_s}{r}\right),
   \qquad
   r_s = \frac{2GM_\odot}{c^2}.
   \]

2. Encode the lattice deformation as the scalar index:

   \[
   n(r) = \frac{1}{\sqrt{1 - 2r_s/r}}.
   \]

3. Compute the bending angle from the standard ray integral using the
   substitution \(u = R_\odot/r\):

   \[
   \Delta\theta
   =
   \frac{2r_s}{R_\odot}
   \int_0^1
   \frac{
     u\left(1 - \frac{2r_s u}{R_\odot}\right)^{-3/2}
   }{
     \sqrt{1-u^2}
   }
   \,du.
   \]

4. Compare with the general-relativistic solar-limb value
   \(1.7517\) arcsec.

## Result

- Computed bending angle: \(1.751728\) arcsec
- Reference value: \(1.7517\) arcsec
- Relative difference: \(0.0016\%\)

## Review note (maintainer)

With the script's own constants the first-order GR value \(2r_s/R_\odot\) is
\(1.750885\) arcsec; the integral's higher-order factor contributes the remaining
\(+0.048\%\). The conventionally quoted \(1.7517\) arcsec is the first-order value
with standard constants, so the \(0.0016\%\) agreement partly reflects that choice of
reference. Recorded here so the comparison's meaning is on file with the number.

## Interpretation

This is a consistency test of the EMC density encoding, not yet a
derivation of light bending from BCC lattice elasticity. The
model-specific step is the encoding

\[
n(r) = \left(\frac{N_\nu(r)}{N_{\text{stat}}}\right)^{-1/2}.
\]

The derivation of \(n(r)\) from the underlying EMC displacement field
remains open.

## Artifacts

- [`research/scripts/m4_3_light_bending_emc_displacement.py`](../scripts/m4_3_light_bending_emc_displacement.py)
- [`research/findings/m4_3_light_bending_emc_displacement.md`](../findings/m4_3_light_bending_emc_displacement.md)

- 
## Reference

Manuscript: Enhanced EWT, version 4.5.6,
DOI: 10.5281/zenodo.17654657.

Relevant sections:

- “The Two Faces of EMC Displacement: Speed and Trajectory”
- “Bridging the Vector Displacement to the Scalar Refractive Index”
- “Asymptotic Continuous Limit and Schwarzschild Equivalence”
