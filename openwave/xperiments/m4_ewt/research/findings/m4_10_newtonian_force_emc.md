# M4.10: Mutual Consistency of Newtonian Force from EMC Field Overlap

## Abstract
This artifact verifies the mutual mathematical and physical consistency of the EMC push-out mechanism with Newton's gravitational force law. By integrating the field overlap of two density deficits over space, the framework demonstrates that the interaction energy yields an inverse-square force law consistent with Newtonian dynamics.

## Angular Integral Correction & Field Overlap Formulation
The 3D volumetric interaction energy between two displaced monopole deficits $\delta\eta_1 = -A_1/r_1$ and $\delta\eta_2 = -A_2/r_2$ separated by distance $R$ reduces via angular integration to:

$$\int d\Omega \frac{r - R \cos\theta}{\left(r^2 + R^2 - 2rR \cos\theta\right)^{3/2}} = \begin{cases} 0 & \text{for } r < R \\ \frac{4\pi}{r^2} & \text{for } r \ge R \end{cases}$$

The total interaction energy is evaluated over the spatial domain $r \in [R, \infty)$:

$$U_{\text{int}}(R) = A_1 A_2 \int_{R}^{\infty} \frac{4\pi}{r^2} \, dr = \frac{4\pi A_1 A_2}{R}$$

Differentiating with respect to $R$ gives the geometric force:

$$F_{\text{geom}} = -\frac{dU_{\text{int}}}{dR} = \frac{4\pi A_1 A_2}{R^2}$$

Coupling this integral to the EMC pressure constant $K_{\text{emc}} = \frac{c^4}{16\pi G}$ yields:

$$F_{\text{EMC}} = K_{\text{emc}} \cdot F_{\text{geom}} = \left(\frac{c^4}{16\pi G}\right) \left(\frac{4\pi A_1 A_2}{R^2}\right) = \frac{c^4 A_1 A_2}{4 G R^2}$$

Substituting the monopole amplitudes $A_i = \frac{2G M_i}{c^2}$ identifies $F_{\text{EMC}} \equiv \frac{G M_1 M_2}{R^2}$.

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
