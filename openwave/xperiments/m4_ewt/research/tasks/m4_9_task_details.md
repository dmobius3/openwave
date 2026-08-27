# M4.9 - Emergent Encoding from Lattice Dynamics

## Status
DONE (post-hoc)

## Criterion
`Gravity: local metric phenomena` — foundational encoding derivation

## Objective
Show in-platform that the two encodings used in M4.3–M4.5 emerge from
lattice dynamics without assuming them.

The required encodings are:

- \(n_\gamma \propto \eta^{-1/2}\)
- \(v_{\text{clock}} \propto \sqrt{\eta}\)

## Method

1. Build a one-dimensional spring-mass lattice.
2. Set microscopic parameters:
   \(m = \eta\), \(k = \eta^2\).
3. Measure pulse speed as a function of uniform \(\eta\).
4. Fit the exponent \(\beta\) in \(v \propto \eta^{\beta}\).
5. Build a density well and excite an oscillator.
6. Measure frequency as a function of core \(\eta\).
7. Fit the exponent \(\gamma\) in \(f \propto \eta^{\gamma}\).

## Result

\[
\beta = 0.502856
\]

\[
\gamma = 0.498371
\]

Both are close to \(0.5\).

## Interpretation

The encodings are emergent. They are not inserted by hand.

## Artifacts

- `research/scripts/m4_9_emergent_encoding.py`
- `research/findings/m4_9_emergent_encoding.md`

## Reference

Enhanced EWT manuscript, version 4.5.9 or later:
[DOI: 10.5281/zenodo.22110605](https://doi.org/10.5281/zenodo.22110605)