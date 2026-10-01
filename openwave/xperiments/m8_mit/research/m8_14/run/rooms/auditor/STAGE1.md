# Stage 1 (V2): independent recomputation of the T4 eigenproblem

R = 1 throughout. Every number below comes from the scripts in this room. `./py run_all.py` runs all of them in order: `fem2d_q2.py`, `fem1d_sectors.py`, `closed_form.py`, `checks.py`, `make_report.py`. Every number is in `stage1.json`. The log of the last full run is `run_all.log`.

## Problem as solved

On the rectangle `0 ≤ y ≤ π`, `−W ≤ w ≤ W`, the problem is `−Δφ = λφ`, with `Δ` the Laplace–Beltrami operator of `dy² + cos²y dw²`. Equivalently, it is the pencil

`a(φ,φ) = ∫ |cos y| φ_y² + φ_w²/|cos y| dy dw`, `m(φ,φ) = ∫ |cos y| φ² dy dw`,

on the admissible class of the spec sheet (§3), with these conditions:
- the seam `φ(π, −w) = −φ(0, w)`;
- `φ = 0` on `w = ±W`;
- `φ = 0` on the collapsed fiber `y = π/2`.

### Readings taken where the spec is underdetermined

1. **Zero data on the fiber, and the realization.** I impose zero data on the fiber, as T4 states. A finite-energy discrete function is in any case forced to be constant on the fiber, because `φ_w²/|cos y|` is not integrable there otherwise, and that constant must be 0 by the data at `w = ±W`. In the `k = 1` sector at `W = 17/10` the transverse index is `a₁ = π/3.4 ≈ 0.924 < 1`. The discrete spaces are conforming subspaces of the form domain, so what they approximate is the form (Friedrichs) realization. This is my reading of what the admissible class selects.
2. **The seam.** It is imposed as a nodal identification. Matching the derivative across the seam is the natural condition of the form, so it is not imposed separately.
3. **The lowest six eigenvalues** are counted with multiplicity. At `W = π/2` they are 2, 6, 6, 12, 12, 12.
4. **`W = 17/10`** lies outside `(0, π/2]` and is treated by the same formulas, as §2 of the spec says.
5. **The bottom** is the lowest eigenvalue. My best values are method B's, because it is the more accurate method. Method A is a fully independent 2D check that does not use separation of variables. The two agree within their stated errors at every width and every index (check C5).
6. **`run_all.py`** is run as `./py run_all.py` and executes the scripts in-process with `runpy`. `./py` cannot be nested: its sandbox refuses a second `sandbox-exec`. A first attempt to call `./py` from inside `run_all.py` failed with "Operation not permitted", and that message included the absolute path of the room's sandbox profile. That path was printed by the system, not by my scripts. I did not open the file.

## Method A: Q2 finite elements directly in (y, w) (`fem2d_q2.py`)

- **Elements.** Biquadratic tensor-product Lagrange elements, so `K = A_y⊗M_w + B_y⊗K_w` and `M = M_y⊗M_w`.
- **y-mesh.** Nodes at distance `s_k = (π/2)(k/N)³` from the fiber on each side, a cubic grading that differs from T4's quadratic one. Q2 midpoints sit at element centres. The 1D integrals use 8-point Gauss quadrature, with weights evaluated in `s` (`|cos y| = sin s`) to avoid cancellation near the fiber.
- **w-mesh.** `M = N/2` uniform Q2 elements on `[−W, W]`. The mesh is symmetric, so `w ↦ −w` maps nodes to nodes.
- **Seam.** Imposed through a prolongation `P`: the row for `(π, w_j)` is `−1` times the dof at `(0, w_{2M−j})`.
- **Eigensolver.** `scipy.sparse.linalg.eigsh` in shift-invert mode at σ = 0, with tolerance 1e-13. The relative residual of every returned pair is stored per level (max 1.1e-11).
- **Levels.** `N = 8, 16, 32, 64, 128, 256` per side. Each refinement doubles `N`, and the meshes are nested. The finest level has 260,865 dofs.
- **Estimates.** Per eigenvalue index, `p` comes from the last three levels and the Richardson extrapolation from the last two. The error estimate is `max(|λ_extrap − λ_finest|, |λ_extrap − λ_extrap(previous triple)|)`.

**Convergence evidence.**
- `p` from the last triple lies between 3.973 and 4.000 for all 42 (width, index) pairs, and `p` from the preceding triple between 3.88 and 3.998 (check C4).
- The successive differences shrink monotonically with a fixed sign.
- This is the O(h⁴) rate expected for Q2 eigenvalues, and it holds even at `W = 17/10`, where the eigenfunction behaves like `s^0.924` at the fiber. The cubic grading is enough there.
- Errors on the bottom range from 6.1e-8 (`W = 1/4`) to 6.7e-10 (`W = 17/10`).

## Method B: exact sine modes in w, high-order FEM in y (`fem1d_sectors.py`)

- **Sectors.** For `φ = u(y) S_k(w)`, with `S_k = sin(kπ(w+W)/(2W))`, `μ_k = (kπ/(2W))²` and `S_k(−w) = (−1)^{k+1} S_k(w)`, the seam becomes `u(π) = (−1)^k u(0)`. The sectors are invariant and the sines are complete, so the spectrum is the union of the sector spectra.
- **Sector cutoff.** Every eigenvalue in sector `k` is at least `μ_k`, because `1/|cos y| ≥ |cos y|`. So sectors are added until `μ_k` exceeds the current sixth-lowest value. Check C6 confirms that adding one more sector changes nothing.
- **Elements.** Degree-8 Lagrange elements on Gauss–Lobatto nodes, with nodes graded as `(π/2)(k/N)⁴` and `N = 2, 4, 8, 16, 32` per side. The degree-6 sequence is kept as a cross-check in the p-direction.
- **Eigensolver.** Dense, solved as `M x = θ K x` for `θ = 1/λ`. The first version solved `K x = λ M x` directly; the graded mesh makes `M` nearly singular, and that version lost accuracy as the mesh was refined. It was discarded.
- **Convergence evidence.**
  - The bottom is stable to 1e-12 relative or better from `N = 4` (`W ≥ 1/2`) or `N = 8` (`W = 1/4`) onward.
  - The degree-6 and degree-8 finest values agree to at most 7.5e-13 on the bottom.
  - Because the finest levels sit at roundoff, `p` is undefined (differences not monotone) at five widths. It is 2.24 and 3.39 at 7/5 and 3/2, where it measures roundoff noise, not an asymptotic order.
- **Error estimate.** It is formed as `max(|λ_extrap − λ_finest|, |λ_extrap − λ_extrap(previous triple)|, |λ_finest(q=8) − λ_finest(q=6)|, 1e-12·λ)`. For every bottom the `1e-12·λ` roundoff floor dominates: the bottom errors run from 1.8e-12 to 4.6e-11. The largest error over all indices is 1.5e-7, for λ₆ at `W = 1/4`. It comes from the previous-triple term, whose levels `N = 2, 4, 8` are preasymptotic there.

## Closed-form cross-check (`closed_form.py`; not used in any estimate)

1. **Unfolding.** In sector `k`, set `t = y` on `(0, π/2)` and continue across the seam by `v(t) = σ u(t+π)` on `(−π/2, 0)`, with `σ = (−1)^k`. This turns the sector problem into one smooth problem on `(−π/2, π/2)` with weight `cos t`.
2. **Reduction to Gegenbauer.** With `x = sin t` it becomes `−((1−x²)u')' + a²u/(1−x²) = λu`, where `a = kπ/(2W)`. Substituting `u = (1−x²)^{a/2} v` gives `(1−x²)v'' − 2(a+1)x v' + (λ − a(a+1))v = 0`.
3. **Eigenvalues.** The polynomial solutions of degree `m` (Gegenbauer `C_m^{(a+1/2)}`) give `λ_{k,m} = (m+a)(m+a+1)`.

Numerically, methods A and B agree with these values within their own error estimates at every width and index (C5′). Two facts here are taken from memory, not derived: the Gegenbauer solutions themselves, and their completeness in `L²((1−x²)^a)`.

## Seam variants: an observation relevant to stage 2

Check C3 shows that these three seam defects leave the eigenvalues unchanged:
- flipping the seam sign (`φ(π,−w) = +φ(0,w)`);
- omitting the reflection (`φ(π,w) = −φ(0,w)`);
- cutting the seam, with zero data on `y = π` only and `y = 0` free.

At `W = 7/5`, `N = 16`, the sign flip and the omitted reflection both reproduce the correct values to 12 digits: 2.380937161975, 6.625355708626, …. The cut seam, extrapolated from `N = 32, 64, 128`, is indistinguishable from method B (ratio 0.011).

The reasons:
- the sign of the unfolding is free;
- `w ↦ −w` is a symmetry of the problem;
- the surface is symmetric about the seam line, so the even and odd modes split between the two halves.

So eigenvalues alone cannot detect seam errors. They have to be read from the scripts. C3 does detect them, from the eigenfunctions.

## Checks and planted defects (`checks.py`)

On the real computation, 73 checks print PASS and none print FAIL. Each check type was also run against a planted defect, and all 8 planted defects fired. That is how I know each check can fail.

| Check | What it tests | Planted defect that made it FAIL |
| --- | --- | --- |
| C1 | Assembly identities: `ΣM_y = 2`, `A_y·1 = 0`, `K_w·1 = 0`, `ΣM_w = 2W`, to 1e-12 | 1-point quadrature (defect 4.1e-3) |
| C2 | Eigen-residual below 1e-9 at the solver's own level | Eigenvalue shifted by 1e-6 relative (residual 1.0e-6) |
| C3 | Reconstructed eigenfunctions satisfy `φ(π,−w) = −φ(0,w)` and vanish on `w = ±W` and on the fiber | Seam sign flipped (violation 2.0); reflection omitted (violation 2.0) |
| C4 | Asymptotic regime for each A sequence: monotone differences and stable `p` (`|Δp| < 0.1`) | Alternating ±1e-6 noise on the `W = 17/10` bottom sequence |
| C5 / C5′ | A vs B, and each vs the closed form, within `errA + errB` | A solved at `W/2` (width misread), through the same pipeline (ratio 7.5e7); B with `μ_k = (kπ/W)²` (ratio 2.0e10) |
| C6 | Adding the next sector does not change the lowest six | Only `k = 1` kept at `W = π/2` (change 22) |
| — | `make_report.py` checks that this file contains `table_stage1.md` verbatim | It printed FAIL on the run made before this file existed |

## Results

"err(bottom)" for B is the bottom error estimate described above. For A it is the combined Richardson estimate. The per-index errors and every level of both methods are in `stage1.json` (`method_A`, `method_B`).

| W/R | method | λ1 (bottom) | λ2 | λ3 | λ4 | λ5 | λ6 | err(bottom) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1/4 | **B (best)** | 45.761602911538 | 60.327973525896 | 76.894344140255 | 95.460714754614 | 116.027085368971 | 138.593455983332 | 4.6e-11 |
| | A extrap. | 45.7616029114 | 60.3279735234 | 76.8943441280 | 95.4607146827 | 116.0270851734 | 138.5934553564 | 6.1e-08 |
| | A p (last 3) | 3.997 | 3.991 | 3.990 | 3.983 | 3.981 | 3.973 | |
| | closed form (check only) | 45.761602911537 | 60.327973525896 | 76.894344140255 | 95.460714754615 | 116.027085368974 | 138.593455983333 | |
| | labels (k,m) | (1,0) | (1,1) | (1,2) | (1,3) | (1,4) | (1,5) | |
| 1/2 | **B (best)** | 13.011197054679 | 21.294382361858 | 31.577567669038 | 43.860752976218 | 45.761602911538 | 58.143938283397 | 1.3e-11 |
| | A extrap. | 13.0111970547 | 21.2943823616 | 31.5775676678 | 43.8607529663 | 45.7616029110 | 58.1439382536 | 9.3e-09 |
| | A p (last 3) | 3.999 | 3.995 | 3.994 | 3.989 | 3.998 | 3.988 | |
| | closed form (check only) | 13.011197054679 | 21.294382361859 | 31.577567669038 | 43.860752976218 | 45.761602911537 | 58.143938283397 | |
| | labels (k,m) | (1,0) | (1,1) | (1,2) | (1,3) | (2,0) | (1,4) | |
| 1 | **B (best)** | 4.038197427067 | 9.179790080656 | 13.011197054679 | 16.321382734245 | 21.294382361858 | 25.462975387836 | 4.0e-12 |
| | A extrap. | 4.0381974271 | 9.1797900806 | 13.0111970546 | 16.3213827341 | 21.2943823615 | 25.4629753857 | 1.9e-09 |
| | A p (last 3) | 3.999 | 3.997 | 3.998 | 3.996 | 3.996 | 3.992 | |
| | closed form (check only) | 4.038197427067 | 9.179790080657 | 13.011197054679 | 16.321382734247 | 21.294382361859 | 25.462975387837 | |
| | labels (k,m) | (1,0) | (1,1) | (2,0) | (1,2) | (2,1) | (1,3) | |
| 7/5 | **B (best)** | 2.380875488666 | 6.624870241229 | 7.279507202100 | 12.868864993794 | 13.767496707228 | 14.695895140301 | 2.4e-12 |
| | A extrap. | 2.3808754887 | 6.6248702412 | 7.2795072020 | 12.8688649937 | 13.7674967071 | 14.6958951388 | 9.5e-10 |
| | A p (last 3) | 3.999 | 3.997 | 3.998 | 3.997 | 3.997 | 3.995 | |
| | closed form (check only) | 2.380875488666 | 6.624870241230 | 7.279507202100 | 12.868864993794 | 13.767496707228 | 14.695895140301 | |
| | labels (k,m) | (1,0) | (1,1) | (2,0) | (1,2) | (2,1) | (3,0) | |
| 3/2 | **B (best)** | 2.143820262429 | 6.238215364819 | 6.480885947323 | 12.332610467214 | 12.669676152094 | 13.011197054679 | 2.1e-12 |
| | A extrap. | 2.1438202624 | 6.2382153648 | 6.4808859473 | 12.3326104671 | 12.6696761520 | 13.0111970533 | 8.4e-10 |
| | A p (last 3) | 3.999 | 3.997 | 3.998 | 3.997 | 3.997 | 3.995 | |
| | closed form (check only) | 2.143820262429 | 6.238215364822 | 6.480885947322 | 12.332610467215 | 12.669676152108 | 13.011197054679 | |
| | labels (k,m) | (1,0) | (1,1) | (2,0) | (1,2) | (2,1) | (3,0) | |
| pi/2 | **B (best)** | 2.000000000000 | 5.999999999999 | 6.000000000000 | 11.999999999999 | 11.999999999999 | 12.000000000000 | 2.0e-12 |
| | A extrap. | 2.0000000000 | 6.0000000000 | 5.9999999999 | 11.9999999999 | 11.9999999999 | 11.9999999988 | 7.7e-10 |
| | A p (last 3) | 4.000 | 3.997 | 3.998 | 3.997 | 3.997 | 3.995 | |
| | closed form (check only) | 2.000000000000 | 6.000000000000 | 6.000000000000 | 12.000000000000 | 12.000000000000 | 12.000000000000 | |
| | labels (k,m) | (1,0) | (1,1) | (2,0) | (1,2) | (2,1) | (3,0) | |
| 17/10 | **B (best)** | 1.777769846306 | 5.263083706641 | 5.625765524887 | 10.455941581004 | 10.959075063805 | 11.473761203470 | 1.8e-12 |
| | A extrap. | 1.7777698463 | 5.2630837066 | 5.6257655249 | 10.4559415799 | 10.9590750637 | 11.4737612034 | 6.7e-10 |
| | A p (last 3) | 3.999 | 3.998 | 3.998 | 3.995 | 3.997 | 3.997 | |
| | closed form (check only) | 1.777769846306 | 5.263083706641 | 5.625765524888 | 10.455941581005 | 10.959075063805 | 11.473761203470 | |
| | labels (k,m) | (1,0) | (2,0) | (1,1) | (3,0) | (2,1) | (1,2) | |

The bottom is positive at every width. It equals `a₁(a₁+1)` with `a₁ = π/(2W)`, which is 2 at `W = π/2` and 1.7777698463 at `W = 17/10`. The labels `(k,m)` come from the closed form. At degenerate eigenvalues the ordering within a cluster is arbitrary.

## Consulted-material manifest

**Files read in this room:** `BRIEF.md`, `auditor_brief.md`, `spec_sheet.md`.

**Files that only appeared in a directory listing (contents not opened):** `py`, `.room.sb`. The sandbox error message showed the absolute path of `.room.sb`.

**Files I created and read back:** `common.py`, `fem2d_q2.py`, `fem1d_sectors.py`, `closed_form.py`, `checks.py`, `make_report.py`, `run_all.py`, `results_A.json`, `results_B.json`, `results_C.json`, `results_checks.json`, `table_stage1.md`, `run_all.log`, `stage1.json`.

**Loaded without my asking:** the session context supplied the user's account e-mail address. It was not used. The Python standard library, numpy and scipy were imported by the scripts. No other file was read, and nothing outside the room was read.

**Results used from memory rather than derived here:**
- Gauss–Legendre and Gauss–Lobatto quadrature, and Lagrange interpolation;
- conforming FEM theory: O(h^{2q}) eigenvalue convergence, algebraic grading toward an `s^a` singularity, and the fact that conforming discretizations approximate the form (Friedrichs) realization;
- Richardson extrapolation;
- the Frobenius exponents `±a` at the pole, and limit-circle for `a < 1`;
- the Gegenbauer equation, its polynomial solutions and their completeness, and the associated-Legendre form of the operator;
- completeness of the Dirichlet sine basis on `[−W, W]`;
- the Rayleigh-quotient lower bound used for the sector cutoff;
- the generalized symmetric eigenproblem (Cholesky of the better-conditioned matrix) and shift-invert Lanczos (scipy).
