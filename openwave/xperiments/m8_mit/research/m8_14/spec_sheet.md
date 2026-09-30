# Spec sheet and worklist

You are one of two independent rooms. Work offline. Do not consult any document, repository or web page: everything you need is below, and every hypothesis is stated here. Report every number you compute, and tune nothing toward any expected value. Set `R = 1` in every number you report.

## 1. The ambient space

`S³(R)` is the round sphere of radius `R` in `ℝ⁴`, with coordinates `(x₁, x₂, x₃, x₄)`. `ℝP³(R) = S³(R)/{±I}` is its quotient by the antipodal map `x ↦ −x`, with the quotient metric.

## 2. The surface

Fix a width `W` with `0 < W ≤ πR/2`, and write

`X(y, w) = R (cos(y/R) cos(w/R), cos(y/R) sin(w/R), sin(y/R), 0)`.

The surface `M(W)` is parametrized on the rectangle `0 ≤ y ≤ πR`, `−W ≤ w ≤ W`, by:
- `X(y, w)` for `0 ≤ y ≤ πR/2`;
- `X(y, −w)` for `πR/2 < y ≤ πR`.

It is then passed to `ℝP³(R)`. In these coordinates:
- The segment `y = πR/2` maps to the single point `p = (0, 0, R, 0)`. Call it the cone point.
- The quotient by `±I` identifies the point `(0, w)` with `(πR, −w)`, since the parametrization sends the second to the antipode of the first. Call this the seam.
- The induced metric is `ds² = dy² + cos²(y/R) dw²`, and the area element is `|cos(y/R)| dy dw`.
- The edge of `M(W)` is the image of the two lines `w = W` and `w = −W`, joined through the seam and meeting at the cone point.

Sections 3 to 5 also use the widths `W/R = π/2` and `W/R = 17/10`. The second lies outside the range above, and the formulas of this section apply to it unchanged, as a surface in its own right.

## 3. The variational problem

- **The normal.** At every point of the rectangle take `ν` to be the constant vector `e₄ = (0, 0, 0, 1)`, the unit normal of the great sphere `x₄ = 0` in `S³(R)`. The seam is the antipodal map, which carries `e₄` to `−e₄`. So `φν` is well defined on `ℝP³(R)` exactly when `φ(πR, −w) = −φ(0, w)`.
- **The variation.** For such a `φ`, the normal-geodesic variation with speed `φ` moves each point `x` to `exp_x(tφ(x)ν(x))`.
- **The admissible class.** The completion, in the norm `∫(|∇φ|² + φ²) dA`, of the smooth functions `φ` on the rectangle that:
  - satisfy `φ(πR, −w) = −φ(0, w)`;
  - vanish on the lines `w = ±W`;
  - vanish near the cone point.
- **What is held fixed.** The edge is held fixed pointwise, and the surface's conic structure at the cone point is not resolved or changed.

## 4. Tasks

- **T1.** Compute the second derivative in `t`, at `t = 0`, of the area of the varied surface, first for smooth `φ` of the kind that the admissible class completes. Write it as a quadratic form in `φ`, state whether any boundary or vertex term appears, and show the steps. Then justify the form's extension to the whole admissible class by continuity.
- **T2.** Consider `−Δ` on the smooth functions of section 3's list, the ones that vanish near the cone point. Is its closure self-adjoint, or does it have other self-adjoint extensions?
  - Answer per transverse sector, meaning the sectors obtained by expanding in the Dirichlet eigenfunctions of the variable `w` on `[−W, W]`.
  - Classify each sector as limit-point or limit-circle at the cone point, for each `W` with `0 < W ≤ πR/2` and for `W/R = 17/10`.
  - Where a sector is limit-circle, state which realization the admissible class of section 3 selects.
- **T3.** Derive, analytically and by any route, the spectrum of `−Δ` on the admissible class, where `Δ` is the Laplace-Beltrami operator of `ds²` in section 2, with the data of section 3. Give the eigenvalues in closed form as functions of `W`, with their labels, and the lowest one.
- **T4.** Compute numerically the lowest eigenvalues of `−Δ` on the admissible class, by this protocol:
  - **Discretization.** P1 finite elements in the coordinates `(y, w)`. The domain is the rectangle `[0, πR] × [−W, W]` with the anti-periodic seam `(0, w) ∼ (πR, −w)`, the area weight `|cos(y/R)|`, and zero data on the arcs and on the collapsed fiber.
  - **Mesh.** Nodes are graded toward the fiber, at distance `(πR/2)(k/N)²` on each side. The initial mesh has `N = 16` per side and 8 cells across, followed by exactly four uniform refinements.
  - **Estimates.** The order `p` is estimated from the last three levels, and the bottom is Richardson-extrapolated from the last two. The error estimate is `|λ_extrap − λ_finest|`.
  - **Widths.** `W/R ∈ {1/4, 1/2, 1, 7/5, 3/2}`, and also `W/R = π/2` and `W/R = 17/10`.
  - **Report, per width.** The lowest six eigenvalues at every level, `p`, the extrapolated bottom, and the error estimate.
- **T5.** From T1 to T4, state for each width the sign of the lowest eigenvalue of the operator whose quadratic form is the second variation, and give its number of negative eigenvalues and its nullity on the admissible class.

## 5. What to return

- **A report.** A markdown file with the derivations for T1, T2, T3 and T5, each step shown, and the tables for T4.
- **Scripts and data.** The scripts that produced every number, and one JSON file holding every number in the tables.
- **Your own assessment.** For each task, what you consider established, and anything you could not settle.
