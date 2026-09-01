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

3. Evaluate the geometric overlap integral
   \(I(R) = \int \nabla\eta_1\cdot\nabla\eta_2\,dV\)
   over the entire infinite domain.

4. Use the coordinate mapping \(r=R/(1-t)\) to avoid truncation
   errors at large radius.

5. Compute the force by central numerical differentiation.

6. Convert the geometric result to physical force with \(K_{\text{emc}}\),
   using the negative-definite energy functional
   \(U_{\text{int}} = -K_{\text{emc}} \int |\nabla\delta\eta|^2 dV\).
   The script reports magnitudes only; the attractive sign is carried
   by this functional.

7. Compare with Newton's law.

## Result

- \(F_{\text{EMC}} = 3.542516096914 \times 10^{22}\ \text{N}\)
- \(F_{\text{Newton}} = 3.542516523099 \times 10^{22}\ \text{N}\)
- Relative difference: \(1.203 \times 10^{-5}\%\)

## Interpretation

The \(1/r^2\) force law follows from the EMC field-overlap geometry.
The \(4\pi\) factor comes from the angular integral, and the radial
integral is performed numerically over the full domain.

The minus sign in the physical force is derived from the negative-definite
elastic energy functional, not imposed by hand. It is the field-theoretic
analogue of the negative Newtonian field energy
\(-\frac{|\nabla\Phi|^2}{8\pi G}\).

Because \(G\), \(c\), \(M\) and \(R\) cancel identically in
\(F_{\text{EMC}} \equiv F_{\text{Newton}}\), the artifact is a
normalization consistency gate, not a derivation of the strength.

## Artifacts

- `research/scripts/m4_10_newtonian_force_emc.py`
- `research/findings/m4_10_newtonian_force_emc.md`

## Reference

Enhanced EWT manuscript, version 4.6.1 or later:
[DOI: 10.5281/zenodo.22144273](https://doi.org/10.5281/zenodo.22144273)
