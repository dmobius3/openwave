# M4.12: Newtonian Force with Geometric Amplitude from the EMC Deficit

## Status
DONE (post-hoc)

## Criterion
`Gravity: Newton limit (GEM)`

## Abstract
This artifact derives the monopole amplitude $A$ from the geometric EMC deficit — using the ratio $N_{\nu,\text{stat}} / N_{\nu,\text{eff}}$ explicitly — rather than writing it as $2 G_{\text{geom}} M / c^2$. In this form the amplitude does not carry $G_{\text{geom}}$, so the Newtonian force gate actually tests $G_{\text{geom}}$ instead of cancelling it. The baseline gate passes at the same precision as M4.10, and a mutation of $G_{\text{geom}}$ in the coupling now breaks the gate at 100%.

## Geometric Amplitude

The amplitude is written as

$$A = \frac{2 M r_e}{m_e} \cdot \frac{\sqrt{X_{\text{eff}}}}{A_\pi^4 \, N_{\text{geom}}^3 \, K_{WC} \sqrt{N_{\nu,\text{stat}}}}$$

with $X_{\text{eff}} = N_{\nu,\text{stat}} / N_{\nu,\text{eff}}$. Every factor on the right is geometric: $A_\pi = 4\pi^3 + \pi^2 + \pi$, $N_{\text{geom}} = 8\pi^4(1-\zeta)$, $K_{WC} = 10$, and $X_{\text{eff}}, N_{\nu,\text{stat}}$ from the M4.7 chain. The anchors $r_e, m_e$ are the measured dimensional quantities; $G_{\text{geom}}$ does not appear.

At $N_{\text{geom}} = 778.8025$:

| Quantity | Value |
| :--- | :--- |
| $G_{\text{geom}}$ | $6.677519975508460 \times 10^{-11}\ \text{m}^3\,\text{kg}^{-1}\,\text{s}^{-2}$ |
| $N_{\text{geom}}$ | 778.8025 |
| $N_{\nu,\text{stat}}$ | $3.295904 \times 10^{52}$ |
| $N_{\nu,\text{eff}}$ | $6.247249 \times 10^{48}$ |
| $X_{\text{eff}}$ | 5275.768324 |
| $A_1$ (geometric) | $2.955552 \times 10^{3}\ \text{m}$ |
| $A_2$ (geometric) | $8.874085 \times 10^{-3}\ \text{m}$ |

The geometric form agrees with $2 G_{\text{geom}} M / c^2$ to $1.5 \times 10^{-16}$ for $A_1$ and $2.0 \times 10^{-16}$ for $A_2$ — machine precision, as expected from algebraic equivalence.

## Force Gate

Baseline: geometric amplitude with the derived $G_{\text{geom}}$ in $K_{\text{emc}}$ and $F_{\text{Newton}}$:

| Quantity | Value |
| :--- | :--- |
| $F_{\text{EMC}}$ | $3.544205553218 \times 10^{22}\ \text{N}$ |
| $F_{\text{Newton}}$ | $3.544205979545 \times 10^{22}\ \text{N}$ |
| Relative difference | $1.203 \times 10^{-5}\%$ |

Mutation test: replace $G_{\text{geom}}$ with $1.0$ in $K_{\text{emc}}$ and $F_{\text{Newton}}$, keeping the geometric amplitude:

| Quantity | Value |
| :--- | :--- |
| $F_{\text{EMC}}$ (mutated) | $2.366650337892 \times 10^{12}\ \text{N}$ |
| $F_{\text{Newton}}$ (mutated) | $5.307668105141 \times 10^{32}\ \text{N}$ |
| Relative difference | $100.000000\%$ |

## Interpretation

The gate now depends on the value of $G_{\text{geom}}$. With the geometric amplitude in place, the factor $G_{\text{geom}}$ appears only in $K_{\text{emc}}$ and in $F_{\text{Newton}}$. Mutating it by orders of magnitude breaks the agreement, which is what makes the gate a test of the geometric coupling rather than a normalisation-consistency check.

**Strength clause.** Circularity is discharged: $G_{\text{geom}}$ enters the force test directly, derived from BCC geometry, without re-entering through an input. Accuracy is stated and not met: the $G_{\text{geom}}$ residual against CODATA is $0.048169\%$, or $21.9\times$ the CODATA 2022 relative uncertainty on $G$. The gate tests $G_{\text{geom}}$ now, but the value it tests does not agree with the measured one within its uncertainty.

**Dimensional anchors.** $r_e, m_e, c$ are the measured anchors used to build the dimensionless geometric ratio. The amplitude $A$ is a length built from these anchors and the BCC geometry.

## Artifacts

- `research/scripts/m4_12_newtonian_force_geometric_amplitude.py`
- `research/findings/m4_12_newtonian_force_geometric_amplitude.md`

## Reference

Enhanced EWT manuscript, version 5.0.0:
[DOI: 10.5281/zenodo.22540635](https://doi.org/10.5281/zenodo.22540635)

Relevant section:

- "Newtonian Force from Interacting EMC Deficits"