# M4.12 - Newtonian Force with Geometric Amplitude

## Status
DONE (post-hoc)

## Criterion
`Gravity: Newton limit (GEM)`

## Objective
Make the Newtonian force gate actually test $G_{\text{geom}}$ by deriving
the monopole amplitude $A$ from the geometric EMC deficit instead of
writing it as $2 G_{\text{geom}} M / c^2$. M4.10 showed that with the old
form the factor $G_{\text{geom}}$ cancels identically in
$F_{\text{EMC}} \equiv F_{\text{Newton}}$, so the gate could not test the
geometric coupling.

## Method

1. Derive $G_{\text{geom}}$ and the geometric quantities $N_{\text{geom}}$, $N_{\nu,\text{stat}}$, $N_{\nu,\text{eff}}$, $X_{\text{eff}}$, $A_\pi$ from the M4.7 emergence engine.

2. Build the amplitude from the deficit ratio:
   $$A = \frac{2 M r_e}{m_e} \cdot \frac{\sqrt{X_{\text{eff}}}}{A_\pi^4 \, N_{\text{geom}}^3 \, K_{WC} \sqrt{N_{\nu,\text{stat}}}}$$
   with no $G$ in the formula.

3. Verify algebraically that the geometric form equals $2 G_{\text{geom}} M / c^2$ to machine precision.

4. State the accuracy half of the strength clause: the $G_{\text{geom}}$ residual against CODATA is $0.048169\%$, or $21.9\times$ the CODATA 2022 relative uncertainty on $G$.

5. Run the baseline force gate with the geometric amplitude and the derived $G_{\text{geom}}$ in $K_{\text{emc}}$ and $F_{\text{Newton}}$.

6. Run a mutation test: replace $G_{\text{geom}}$ with $1.0$ in the coupling only, and confirm the gate fails.

## Result

Baseline:

- $F_{\text{EMC}} = 3.544205553218 \times 10^{22}\ \text{N}$
- $F_{\text{Newton}} = 3.544205979545 \times 10^{22}\ \text{N}$
- Relative difference: $1.203 \times 10^{-5}\%$ (PASS)

Mutation ($G_{\text{geom}} \to 1.0$):

- Relative difference: $100.000000\%$ (FAIL, as expected)

Geometric amplitude:

- $A_1 = 2.955552 \times 10^{3}\ \text{m}$
- $A_2 = 8.874085 \times 10^{-3}\ \text{m}$
- Agreement with $2 G_{\text{geom}} M / c^2$ at machine precision.

## Interpretation

The gate now depends on $G_{\text{geom}}$: the amplitude is geometric, the
coupling carries $G_{\text{geom}}$, and mutating $G_{\text{geom}}$ breaks
the agreement. M4.10 could not do this because $G_{\text{geom}}$ sat on
both sides and cancelled.

Strength clause: circularity is discharged (the trinity enters directly
without re-entering through an input); accuracy is stated and not met
($21.9\times$ the CODATA uncertainty on $G$).

The dimensional anchors $r_e, m_e, c$ are named at the point of use.

## Artifacts

- `research/scripts/m4_12_newtonian_force_geometric_amplitude.py`
- `research/findings/m4_12_newtonian_force_geometric_amplitude.md`

## Reference

Enhanced EWT manuscript, version 5.0.0:
[DOI: 10.5281/zenodo.22540635](https://doi.org/10.5281/zenodo.22540635)