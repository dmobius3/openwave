# M4.10: Mutual Consistency of Newtonian Force from EMC Field Overlap

## Abstract
This artifact verifies the mutual mathematical and physical consistency of the EMC push-out mechanism with Newton's gravitational force law. By integrating the field overlap of two density deficits over space, the framework demonstrates that the interaction energy yields an inverse-square force law consistent with Newtonian dynamics.

## Angular Integral Correction & Field Overlap Formulation
The 3D volumetric interaction energy between two displaced monopole deficits $\delta\eta_1 = -A_1/r_1$ and $\delta\eta_2 = -A_2/r_2$ separated by distance $R$ reduces via angular integration to:

$$\int d\Omega \frac{r - R \cos\theta}{\left(r^2 + R^2 - 2rR \cos\theta\right)^{3/2}} = \begin{cases} 0 & \text{for } r < R \\ \frac{4\pi}{r^2} & \text{for } r \ge R \end{cases}$$

The total geometric overlap integral $I(R)$ is evaluated over the spatial domain $r \in [R, \infty)$:

$$I(R) = \int \nabla \delta\eta_1 \cdot \nabla \delta\eta_2 \, dV = A_1 A_2 \int_{R}^{\infty} \frac{4\pi}{r^2} \, dr = \frac{4\pi A_1 A_2}{R}$$

Differentiating with respect to $R$ gives the geometric gradient magnitude:

$$F_{\text{geom}} = -\frac{dI}{dR} = \frac{4\pi A_1 A_2}{R^2}$$

*(Note: $F_{\text{geom}}$ is a purely geometric intermediate gradient magnitude, superseded by the physical field-energy force $F_{\text{EMC}}$ in the next section.)*

## Physical Interaction Energy and Sign Convention

The physical interaction energy $U_{\text{int}}(R)$ follows from the negative-definite field energy functional $E[\delta\eta]$ of the overlapping EMC deficits, serving as the exact field-theoretic analogue of the negative Newtonian field energy $-\frac{1}{8\pi G}\int |\nabla\Phi|^2 dV$:

$$
E[\delta\eta] = -\frac{1}{2} K_{\text{emc}} \int |\nabla \delta\eta|^2 \, dV
$$

For two superposed monopole deficits $\delta\eta = \delta\eta_1 + \delta\eta_2$ the two self-energy terms are independent of $R$, so the interaction energy is the cross term of $E$:

$$
U_{\text{int}}(R) = -\frac{1}{2} K_{\text{emc}} \cdot 2 \int \nabla \delta\eta_1 \cdot \nabla \delta\eta_2 \, dV = -K_{\text{emc}} I(R) = -\frac{4\pi K_{\text{emc}} A_1 A_2}{R}
$$

The minus sign is therefore derived from the energy functional, not imposed by hand.

Differentiating with respect to $R$ defines the attractive physical force $F_{\text{EMC}}$:

$$
F_{\text{EMC}} = -\frac{dU_{\text{int}}}{dR} = -\frac{4\pi K_{\text{emc}} A_1 A_2}{R^2}
$$

which is strictly attractive ($F_{\text{EMC}} < 0$ pointing inward towards decreasing $R$). The shipped script reports magnitudes only; the sign is carried by the energy functional.

Coupling this to the EMC pressure constant $K_{\text{emc}} = \frac{c^4}{16\pi G}$ and substituting the monopole amplitudes $A_i = \frac{2G M_i}{c^2}$ yields:

$$
F_{\text{EMC}} = -\frac{G M_1 M_2}{R^2}
$$

which matches the attractive Newtonian force $F_{\text{Newton}}$.

## Structural Discrimination Analysis
Because $G$, $c$, $M$, and $R$ cancel identically in the formal equality $F_{\text{EMC}} \equiv F_{\text{Newton}}$, the artifact operates as a consistency gate verifying that the combined normalization of $A$, $K_{\text{emc}}$, and the spherical overlap factor $4\pi$ is correct.

Mutating any single structural input breaks the agreement:

| Parameter Mutation | Observed Rel. Diff. | Gate Result |
| :--- | :--- | :--- |
| Baseline (Shipped) | $1.203 \times 10^{-5}\%$ | **PASS** |
| Angular Factor ($4\pi \to 3\pi$) | $25.0\%$ | **FAIL** |
| Monopole Amplitude ($2GM/c^2 \to GM/c^2$) | $75.0\%$ | **FAIL** |
| Coupling Denominator ($16\pi \to 8\pi$) | $100.0\%$ | **FAIL** |

## Precision Note
The residual of $\sim 1.203 \times 10^{-5}\%$ reported under the coordinate mapping $r(t) = \frac{R}{1-t}$ is due to accumulated floating-point roundoff from midpoint summation over a constant transformed integrand $\frac{4\pi}{R}$, rather than a physical error.

## Reference

Enhanced EWT manuscript, version 4.6.1 or later:
[DOI: 10.5281/zenodo.22144273](https://doi.org/10.5281/zenodo.22144273)

Relevant section:

- "Newtonian Force from Interacting EMC Deficits"
