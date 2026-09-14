# M4.10 - Newtonian Force from EMC Push-Out Pressure

## Status
DONE (post-hoc)

## Criterion
`Gravity: Newton limit (GEM)`

## Objective
Show that the attractive $1/r^2$ force law follows from the EMC
push-out mechanism using exact-domain numerical field integration, with
both the monopole amplitude and the EMC coupling built from the
self-consistent geometric trinity $G_{\text{geom}}, \hbar_{\text{geom}},
\lambda_l$ derived in the M4.7 emergence engine.

## Method

1. Derive $G_{\text{geom}}$ from the M4.7 self-consistent trinity.

2. State the accuracy half of the strength clause: the $G_{\text{geom}}$
   residual against CODATA is $0.048169\%$, or $21.9\times$ the CODATA
   2022 relative uncertainty on $G$ (22 ppm).

3. Define monopole density deficits $\delta\eta = -A/r$ with
   $A = 2 G_{\text{geom}} M / c^2$.

4. Evaluate the geometric overlap integral
   $I(R) = \int \nabla\eta_1\cdot\nabla\eta_2\,dV$
   over the entire infinite domain.

5. Use the coordinate mapping $r=R/(1-t)$ to avoid truncation
   errors at large radius.

6. Compute the force by central numerical differentiation.

7. Convert the geometric result to physical force with
   $K_{\text{emc}} = c^4 / (16\pi G_{\text{geom}})$, using the
   negative-definite energy functional
   $E[\delta\eta] = -\frac{1}{2} K_{\text{emc}} \int |\nabla\delta\eta|^2 dV$.
   Its $R$-dependent cross term is $U_{\text{int}} = -K_{\text{emc}} I(R)$.
   The script reports magnitudes only; the attractive sign is carried by
   this functional.

8. Compare with Newton's law using the same $G_{\text{geom}}$.

## Result

- $G_{\text{geom}} = 6.677519975508460 \times 10^{-11}\ \text{m}^3\,\text{kg}^{-1}\,\text{s}^{-2}$
- $\lambda_l$ (derived) $= 1.616646406608278 \times 10^{-35}\ \text{m}$
- $G_{\text{geom}}$ residual vs CODATA: $0.048169\%$
- residual / CODATA uncertainty: $21.9\times$
- $|F_{\text{EMC}}| = 3.544205553218 \times 10^{22}\ \text{N}$
- $|F_{\text{Newton}}| = 3.544205979545 \times 10^{22}\ \text{N}$
- Relative difference: $1.203 \times 10^{-5}\%$

## Interpretation

The $1/r^2$ force law follows from the EMC field-overlap geometry.
The $4\pi$ factor comes from the angular integral, and the radial
integral is performed numerically over the full domain.

The minus sign in the physical force is derived from the negative-definite
field energy functional $E[\delta\eta] = -\frac{1}{2} K_{\text{emc}} \int |\nabla\delta\eta|^2 dV$,
not imposed by hand. It is the field-theoretic analogue of the negative Newtonian field energy
$-\frac{1}{8\pi G}\int |\nabla\Phi|^2 dV$.

The circularity half of the strength clause is discharged, and the accuracy half is stated and not met: the residual is 21.9x the CODATA uncertainty on $G$. Circularity: $G_{\text{geom}}$ is derived from BCC geometry without re-entering the chain through an input, since the Planck length has moved from input to output in the M4.7 trinity.

Because $G_{\text{geom}}$ enters both $A$ and $K_{\text{emc}}$, it cancels identically in $F_{\text{EMC}} \equiv F_{\text{Newton}}$, so the force comparison remains a normalization consistency gate, not a derivation of the strength.

The dimensional anchors $r_e, m_e, c$ are the measured quantities used to build the dimensionless geometric ratio. The derivation produces $G_{\text{geom}}$ as a dimensionless ratio against $c^2 r_e / m_e$, not $G$ from nothing.

## Artifacts

- `research/scripts/m4_10_newtonian_force_emc.py`
- `research/findings/m4_10_newtonian_force_emc.md`

## Reference

Enhanced EWT manuscript, version 5.0.0 or later:
[DOI: 10.5281/zenodo.22540635](https://doi.org/10.5281/zenodo.22540635)
