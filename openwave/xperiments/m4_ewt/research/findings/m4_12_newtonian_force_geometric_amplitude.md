# M4.12: Newtonian Force with the Amplitude in M4.7 Chain Factors

## Status
BACKLOG. This record documents the amplitude written in M4.7 chain factors; the task's amplitude-derivation clause stays open (see Interpretation).

## Criterion
`Gravity: Newton limit (GEM)`

## Abstract
This artifact writes the monopole amplitude $A$ in the M4.7 chain factors, using the ratio $N_{\nu,\text{stat}} / N_{\nu,\text{eff}}$ explicitly, without the symbol $G_{\text{geom}}$. The product is the engine's $G_{\text{EWT}}$ formula multiplied by $2M/c^2$, term for term, so $A = 2 G_{\text{geom}} M / c^2$ by value and $G_{\text{geom}}$ still cancels identically in $F_{\text{EMC}} \equiv F_{\text{Newton}}$. The force comparison therefore remains a normalization-consistency gate, as in M4.10, not a test of $G_{\text{geom}}$. The baseline gate passes at the same precision as M4.10. Mutating $G_{\text{geom}}$ in the coupling only breaks the gate at 100%, which shows that the two copies of the formula agree, not that the value of $G_{\text{geom}}$ is right.

## Amplitude in Chain Factors

The amplitude is written as

$$A = \frac{2 M r_e}{m_e} \cdot \frac{\sqrt{X_{\text{eff}}}}{A_\pi^4 \, N_{\text{geom}}^3 \, K_{WC} \sqrt{N_{\nu,\text{stat}}}}$$

with $X_{\text{eff}} = N_{\nu,\text{stat}} / N_{\nu,\text{eff}}$. Every factor on the right is geometric: $A_\pi = 4\pi^3 + \pi^2 + \pi$, $N_{\text{geom}} = 8\pi^4(1-\zeta)$, $K_{WC} = 10$, and $X_{\text{eff}}, N_{\nu,\text{stat}}$ from the M4.7 chain. The anchors $r_e, m_e$ are the measured dimensional quantities. The symbol $G_{\text{geom}}$ does not appear, but the product is the engine's $G_{\text{EWT}}$ formula (`gravity_sector` in the M4.7 engine) multiplied by $2M/c^2$, term for term, so $A = 2 G_{\text{geom}} M / c^2$ by value.

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

The geometric form agrees with $2 G_{\text{geom}} M / c^2$ to $1.5 \times 10^{-16}$ for $A_1$ and $2.0 \times 10^{-16}$ for $A_2$, machine precision: the two are the same expression, so this agreement is an identity, not a confirmation.

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

The gate does not depend on the value of $G_{\text{geom}}$. The amplitude carries $G_{\text{geom}}$ by value through the chain factors and the coupling carries it by symbol, so it sits on both sides and cancels identically in $F_{\text{EMC}} \equiv F_{\text{Newton}}$, as in M4.10: changing any input to the chain moves $G_{\text{geom}}$ and the gate still passes. The mutation above changes $G_{\text{geom}}$ in the coupling only, one copy of the formula and not the other, a state the chain cannot produce, so it confirms the identity and does not test the value. The force comparison remains a normalization-consistency gate. Making it test $G_{\text{geom}}$ needs an amplitude whose value is fixed without passing through the $G_{\text{EWT}}$ formula, which is the task's open clause.

**Strength clause.** Circularity is discharged: $G_{\text{geom}}$ enters the force test directly, derived from BCC geometry, without re-entering through an input. Accuracy is stated and not met: the $G_{\text{geom}}$ residual against CODATA is $0.048169\%$, or $21.9\times$ the CODATA 2022 relative uncertainty on $G$. At that distance the derived value does not agree with the measured $G$ within its uncertainty.

**Dimensional anchors.** $r_e, m_e, c$ are the measured anchors used to build the dimensionless geometric ratio. The amplitude $A$ is a length built from these anchors and the BCC geometry.

## Artifacts

- `research/scripts/m4_12_newtonian_force_geometric_amplitude.py`
- `research/findings/m4_12_newtonian_force_geometric_amplitude.md`

## Reference

Enhanced EWT manuscript, version 5.0.0:
[DOI: 10.5281/zenodo.22540635](https://doi.org/10.5281/zenodo.22540635)

Relevant section:

- "Newtonian Force from Interacting EMC Deficits"
