# Formalization of Nonlinear Stabilisation Mechanisms for the Electron Soliton in the Energy Wave Theory Framework

### Recommendation for OpenWave M4 Implementation

**Author:** Łukasz Smoliński  
**Affiliation:** Independent Researcher, 61-160 Czapury, Poland  
**Version:** 1.0  
**Date:** 2026-07-26  
**DOI:** [10.5281/zenodo.21557369](https://doi.org/10.5281/zenodo.21557369)

---

## Download PDF

The full PDF is available on Zenodo:  
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21557369.svg)](https://doi.org/10.5281/zenodo.21557369)

---

## Abstract

The Energy Wave Theory (EWT) provides a geometric foundation for particle physics, wherein the electron is modelled as a stable soliton composed of ten wave centres arranged in a 1-3-6 tetrahedral geometry. While the geometric derivation of coupling constants ($G$, $\alpha$, $\epsilon_M$) has been established, the dynamical stabilisation of this configuration requires a nonlinear term in the wave equation that counteracts dispersion. This paper formalises the nonlinear stabilisation mechanisms, extending the existing framework with Gaussian density profiles, quintic saturation, modular density profiles, and shell structures for higher-generation leptons. The formalisms are presented as extensions to the core EWT theory, providing a foundation for numerical validation.

---

## Table of Contents

1. Introduction
2. The BCC Vacuum Lattice and Density Profiles
   - 2.1 Statutory Background
   - 2.2 General Density Profile
   - 2.3 Profile Functions
3. Nonlinear Wave Equation
   - 3.1 General Form
   - 3.2 The Magnetic Deficit and Coupling Coefficient
4. Variants of the Nonlinear Term
   - 4.1 Variant A: Pure Cubic Nonlinearity
   - 4.2 Variant B: Quintic Saturation
   - 4.3 Variant C: Density-Modulated Cubic
   - 4.4 Variant D: Gaussian Density Profile
   - 4.5 Variant E: Gaussian Profile + Quintic Saturation (Recommended)
   - 4.6 Variant F: 1-3-6 Source-Driven Dynamics
5. The 1-3-6 Geometry and Topological Selection
   - 5.1 Structural Emergence
   - 5.2 Topological Selection Rule
6. Pressure Force from Vacuum Density Gradient
7. Shell Structures for Higher Generations
   - 7.1 Core + Shell (Onion Model)
   - 7.2 Spherical Shell (Golden Angle)
   - 7.3 Toroidal Shell
8. Modular Density Profiles
9. The Degraded EMC Wall
10. Conclusion
11. References

---

## 1. Introduction

The Energy Wave Theory (EWT) [1, 2, 3] posits that fundamental particles are stable solitons in an elastic medium—the Body-Centered Cubic (BCC) lattice of Elastic Medium Constituents (EMCs). The electron is modelled as a composite soliton of ten wave centres arranged in a 1-3-6 tetrahedral geometry [1, 4].

The geometric derivation of coupling constants is well-established [1]. However, the dynamical stabilisation of the 1-3-6 configuration requires a nonlinear term $\mathcal{F}$ in the wave equation that counteracts dispersion. This paper formalises the nonlinear stabilisation mechanisms, extending the existing framework with new profiles and structures.

The present document constitutes a proposal for the implementation of these formalisms within the M4 vector-field model of the OpenWave simulation platform [5]. OpenWave is an open-source computational framework for exploring fundamental physics through classical field theory enriched with topology and nonlinearity, providing a unified environment for testing candidate field-theoretic models against particle-scale phenomena [5]. The M4 model, as a vector-field PDE solver, serves as the natural substrate for the nonlinear stabilisation mechanisms described herein, enabling direct numerical validation of the electron soliton stability and K-selectivity predictions.

---

## 2. The BCC Vacuum Lattice and Density Profiles

### 2.1 Statutory Background

The vacuum is modelled as a BCC lattice of spherical EMCs. The statutory background density is [1]:

$$N_{\nu,\text{stat}} = \left( \frac{r_{\nu}}{2\lambda_l e} \right)^3 \approx 3.30 \times 10^{52},$$

where $r_{\nu}$ is the neutrino soliton radius and $\lambda_l$ is the Planck length.

### 2.2 General Density Profile

The soliton creates a local density deficit:

$$\rho(r) = \rho_0 - \Delta\rho \cdot f(r),$$

where $f(r)$ is a radial profile function. The deficit fraction (modulation factor) is:

$$\text{mod}(r) = 1 - \frac{\rho(r)}{\rho_0} = \frac{\Delta\rho}{\rho_0} f(r).$$

### 2.3 Profile Functions

Several radial profile functions $f(r)$ are considered:

**Gaussian profile:**
$$f_{\text{gauss}}(r) = \exp\left(-\frac{r^2}{R^2}\right),$$

where $R$ is the characteristic soliton radius. This profile is smooth everywhere, has no wall, and asymptotically approaches the statutory background.

**Exponential deficit:**
$$f_{\text{exp}}(r) = \exp\left(-\frac{r}{R}\right).$$

**Flat-bottom (sigmoid) profile:**
$$f_{\text{flat}}(r) = \frac{1}{2}\left(1 - \tanh\left(\frac{r - R}{\sigma}\right)\right),$$

where $\sigma$ controls the transition width. This profile has a constant deficit inside the soliton.

**Flat-bottom profile with EMC wall:**
$$f_{\text{wall}}(r) = f_{\text{flat}}(r) + \frac{\rho_{\text{wall}} - \rho_0}{\Delta\rho} \exp\left(-\frac{(r - R_{\text{wall}})^2}{2\sigma_{\text{wall}}^2}\right).$$

---

## 3. Nonlinear Wave Equation

### 3.1 General Form

The vector displacement field $\boldsymbol{\Psi}(\mathbf{r},t)$ satisfies [1]:

$$\left( \frac{\partial^2}{\partial t^2} - c^2 \nabla^2 \right) \boldsymbol{\Psi}(\mathbf{r},t) + \mathcal{F}(\boldsymbol{\Psi}) = 0,$$

where $\mathcal{F}(\boldsymbol{\Psi})$ is the nonlinear restoring force.

### 3.2 The Magnetic Deficit and Coupling Coefficient

The magnetic deficit $\epsilon_M$ is the fundamental lattice response parameter [1]:

$$\epsilon_M = \frac{1}{N_{\text{final}}\pi^3} \approx \frac{1}{8\pi^7}.$$

The nonlinear coupling coefficient $\gamma$ is the inverse:

$$\gamma = \frac{1}{\epsilon_M} = N_{\text{final}}\pi^3 \approx 2.414 \times 10^4.$$

Two values are available:
- $\gamma_{\text{geo}} = 8\pi^7$ — the pure geometric limit.
- $\gamma_{\text{final}} = N_{\text{final}}\pi^3$ — calibrated to the measured fine-structure constant.

---

## 4. Variants of the Nonlinear Term

### 4.1 Variant A: Pure Cubic Nonlinearity

The simplest form is the cubic (self-focusing) nonlinearity:

$$\mathcal{F}(\boldsymbol{\Psi}) = \gamma \|\boldsymbol{\Psi}\|^2 \boldsymbol{\Psi}.$$

### 4.2 Variant B: Quintic Saturation

To prevent amplitude blow-up, a quintic (fifth-order) saturation term is added:

$$\mathcal{F}(\boldsymbol{\Psi}) = \gamma \|\boldsymbol{\Psi}\|^2 \boldsymbol{\Psi} - \kappa \|\boldsymbol{\Psi}\|^4 \boldsymbol{\Psi},$$

where $\kappa > 0$ is the saturation coefficient.

### 4.3 Variant C: Density-Modulated Cubic

The nonlinearity is modulated by the EMC density deficit:

$$\mathcal{F}(\boldsymbol{\Psi}) = \gamma \cdot \text{mod}(r) \cdot \|\boldsymbol{\Psi}\|^2 \boldsymbol{\Psi}.$$

### 4.4 Variant D: Gaussian Density Profile

Using the Gaussian profile:

$$\mathcal{F}(\boldsymbol{\Psi}) = \gamma \cdot \text{mod}_{\text{gauss}}(r) \cdot \|\boldsymbol{\Psi}\|^2 \boldsymbol{\Psi},$$

where:

$$\text{mod}_{\text{gauss}}(r) = \frac{\Delta\rho}{\rho_0} \exp\left(-\frac{r^2}{R^2}\right).$$

### 4.5 Variant E: Gaussian Profile + Quintic Saturation (Recommended)

This variant combines the smooth Gaussian profile with quintic saturation:

$$\mathcal{F}(\boldsymbol{\Psi}) = \gamma \cdot \text{mod}_{\text{gauss}}(r) \cdot \|\boldsymbol{\Psi}\|^2 \boldsymbol{\Psi} - \kappa \cdot \text{mod}_{\text{gauss}}(r) \cdot \|\boldsymbol{\Psi}\|^4 \boldsymbol{\Psi}.$$

### 4.6 Variant F: 1-3-6 Source-Driven Dynamics

In the native EWT model, the wave centres are discrete sources:

$$\left( \frac{\partial^2}{\partial t^2} - c^2 \nabla^2 \right) \boldsymbol{\Psi} + \gamma \|\boldsymbol{\Psi}\|^2 \boldsymbol{\Psi} = \mathbf{J}_{1-3-6}(\mathbf{r},t),$$

where:

$$\mathbf{J}_{1-3-6}(\mathbf{r},t) = \mathbf{j}_{\text{core}}(\mathbf{r} - \mathbf{r}_0) + \sum_{i=1}^{3} \mathbf{j}_{\text{inner}}(\mathbf{r} - \mathbf{r}_i(t)) + \sum_{j=1}^{6} \mathbf{j}_{\text{outer}}(\mathbf{r} - \mathbf{r}_j(t)).$$

---

## 5. The 1-3-6 Geometry and Topological Selection

### 5.1 Structural Emergence

The 1-3-6 arrangement is the unique self-consistent candidate that satisfies:

1. The topological winding number $Q = 10$ on $S^2$.
2. The cubic nonlinearity with coefficient $\gamma$.
3. The octahedral symmetry of the $Im\bar{3}m$ BCC lattice.

The winding number is defined as [1]:

$$Q = \frac{1}{4\pi} \int_{S^2} \boldsymbol{\Psi}^* \cdot \omega \in \mathbb{Z},$$

where $\omega$ is the normalised area form on $S^2$.

### 5.2 Topological Selection Rule

The selection of $K_{\text{WC}} = 10$ is enforced by the topological class of the soliton. Because $Q$ is a topological invariant, it cannot change under continuous deformations of the wave field.

---

## 6. Pressure Force from Vacuum Density Gradient

For density-modulated profiles, a pressure force arises from the gradient of the EMC density [1]:

$$\mathbf{F}_{\text{pressure}} = -\alpha \nabla \rho(r),$$

where $\alpha$ is the pressure coupling strength. This force acts on each wave centre, pulling it toward the centre of the density deficit.

---

## 7. Shell Structures for Higher Generations

### 7.1 Core + Shell (Onion Model)

For higher-generation leptons (muon, tau), the native electron core is preserved and additional wave centres are placed in shells [1]:

$$K_{\text{total}} = K_{\text{core}} + \sum_{n=1}^{N_{\text{shells}}} K_{\text{shell},n},$$

where $K_{\text{core}} = 10$.

### 7.2 Spherical Shell (Golden Angle)

Points are distributed on a sphere of radius $R_{\text{shell}}$ using golden angle phyllotaxis:

$$\theta_i = \arccos\left(1 - \frac{2i}{K_{\text{shell}}}\right), \quad \phi_i = \pi(3 - \sqrt{5})\,i.$$

### 7.3 Toroidal Shell

Points are distributed on a torus surface:

$$x_i = (R + r\cos\psi_i)\cos\theta_i,$$
$$y_i = (R + r\cos\psi_i)\sin\theta_i,$$
$$z_i = r\sin\psi_i,$$

where $\theta_i = 2\pi i/K_{\text{shell}}$ and $\psi_i$ follows the golden angle.

---

## 8. Modular Density Profiles

A modular system is proposed for systematic comparison of density profiles:

| Profile | Name | Formula for $f(r)$ |
|---------|------|--------------------|
| 1 | Gaussian | $\exp(-r^2/R^2)$ |
| 2 | Exponential | $\exp(-r/R)$ |
| 3 | Flat-bottom (sigmoid) | $\frac{1}{2}(1 - \tanh((r-R)/\sigma))$ |
| 4 | Flat-bottom + wall | $f_3(r) + \frac{\rho_{\text{wall}} - \rho_0}{\Delta\rho} \exp(-(r-R_{\text{wall}})^2/(2\sigma_{\text{wall}}^2))$ |

---

## 9. The Degraded EMC Wall

The Degraded EMC Wall is a spherical shell at radius $r_{\text{wall}}$ where the EMC density transitions from the soliton interior to the statutory background [1]. The wall radius is estimated as:

$$r_{\text{wall}} \sim \frac{\lambda_C}{2\pi} = \frac{r_e}{\alpha},$$

where $\lambda_C$ is the electron Compton wavelength.

---

## 10. Conclusion

The formalisms presented in this paper extend the Energy Wave Theory framework with:

- Gaussian density profiles for smooth vacuum response.
- Quintic saturation for amplitude limiting.
- Modular density profiles for systematic comparison.
- Shell structures for higher-generation leptons.
- The Degraded EMC Wall as a natural boundary condition.

These extensions provide a foundation for numerical validation of the electron soliton stability and K-selectivity.

---

## References

[1] Smoliński, Ł. (2026). *The Geometric Identity of Gravity and Dimensional Unification Resolving $\alpha$, Lepton $(g-2)_l$, Weinberg, and Cabibbo Mixing*. Version 4.5.2. DOI: [10.5281/zenodo.17654657](https://doi.org/10.5281/zenodo.17654657).

[2] Yee, J. (2019). *The Geometry of Particles and the Explanation of their Creation and Decay*. ResearchGate Preprint. DOI: [10.13140/RG.2.2.14966.14401](https://doi.org/10.13140/RG.2.2.14966.14401).

[3] Yee, J. and Gardi, L. (2019). *The Geometry of Spacetime and the Unification of the Electromagnetic, Gravitational and Strong Forces*. ResearchGate Preprint. DOI: [10.13140/RG.2.2.23094.24642](https://doi.org/10.13140/RG.2.2.23094.24642).

[4] Yee, J. (2020). *The Geometry of Particle Standing Waves*. ResearchGate Preprint. DOI: [10.13140/RG.2.2.27401.88169](https://doi.org/10.13140/RG.2.2.27401.88169).

[5] OpenWave Collaboration (2026). *OpenWave Simulation Platform*. GitHub Repository: [https://github.com/openwave-labs/openwave](https://github.com/openwave-labs/openwave).

[6] Yee, J. and Smoliński, Ł. (2025). *The Geometric Black Hole: The Role of $\epsilon_{G}$ in Extreme Wave Geometries*. DOI: [10.5281/zenodo.17397981](https://doi.org/10.5281/zenodo.17397981).

---

## Companion Documents

- **Numerical Test Report (K-Selectivity Validation):** [10.5281/zenodo.21591354](https://doi.org/10.5281/zenodo.21591354)
- **Main EWT Manuscript (The Geometric Identity of Gravity...):** [10.5281/zenodo.17654657](https://doi.org/10.5281/zenodo.17654657)