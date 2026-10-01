# Stage 2: audit of the returns of Room A and Room B

R = 1 throughout. My comparison values are my stage-1 values in `stage1.json`, which are unchanged. **No stage-1 correction was needed**, and none was made. Every verdict below is in `stage2.json`, which `make_stage2.py` writes from the outputs of the audit scripts.

## Scripts written at stage 2 (all new; no stage-1 file was edited)

| Script | What it does | Output |
| --- | --- | --- |
| `s2_rerun.py` | Copies each room's scripts to `rerun_a/` and `rerun_b/` and runs each room's `run_all.py` there, with this interpreter (both exit 0) | `rerun_a.log`, `rerun_b.log` |
| `s2_compare.py` | C1: rerun JSON vs return JSON. C2: report T4 tables vs JSON. C3: report T3/T5 numbers vs JSON | `s2_compare.json` |
| `s2_audit.py` | K1–K7 on the rooms' own assembly code (copied to `lib_a/` and `lib_b/`), at all 7 widths and all 5 levels | `s2_audit.json`, `s2_audit.log` |
| `s2_fiber_b.py` | F1: Room B's "exactly zero" fiber gradient | `s2_fiber_b.json` |
| `s2_symbolic.py` | G1–G6: independent checks behind the grades | `s2_symbolic.json`, `s2_symbolic.log` |
| `s2_regrade.py` | Sensitivity to the other reading of "uniform refinement". Informational; it prints no PASS | `s2_regrade.json`, `s2_regrade.log` |
| `make_stage2.py` | Writes `stage2.json` | `stage2.json` |

## Checks, and how I know each can fail

On the real returns every check line printed PASS: 238 in `s2_audit`, 6 in `s2_compare`, 2 in `s2_fiber_b` and 7 in `s2_symbolic`. None printed FAIL. Every check type has a planted twin, and every twin printed FAIL:

| Check | What it tests | Planted defect, and what it printed |
| --- | --- | --- |
| C1 | Every numeric leaf of the rerun JSON equals the return JSON (rel 1e-9, absolute floor 1e-12). A: 560 leaves, max 8.4e-11. B: 558 leaves, max 4.2e-11 | One level value altered by 1e-6: FAIL |
| C2 | Every printed T4 number equals the JSON at printed precision (A: 231 numbers, B: 231) | Last digit of one eigenvalue changed: FAIL |
| C3 | Every printed T3 lowest-six value and T5 value equals the JSON (49 per room) | One digit changed: FAIL |
| K1 | The mesh matches the protocol: level-0 nodes `π/2 ∓ (π/2)(k/16)²`, 8 w-cells, levels nested by midpoints, ndof 217…64897 | Expecting a cubic grading: FAIL |
| K2 | Seam read from the rooms' own dof maps. The rebuilt bottom eigenvector satisfies `φ(π,−w)=−φ(0,w)` (residual 0), is zero on `w=±W` and on the fiber, and correlates with the closed-form bottom eigenfunction (≥ 0.999997) | Each room's periodic seam: residual 2.0, correlation 0.000: FAIL (both rooms) |
| K3 | Eigen-residual below 1e-8 for all six pairs at every level and width. Max 4.1e-12 (A) and 4.1e-12 (B) | Eigenvalue shifted by 1e-6: residual 1.0e-6: FAIL |
| K4 | Sylvester inertia (dense LDLᵀ, levels 0–2): 0 eigenvalues below `0.999999·λ₁`, exactly 6 below `(λ₆+λ₇)/2` | A solver that missed λ₁: count 7: FAIL |
| K5 | My eigsh on the rooms' matrices reproduces their JSON (max rel 5.0e-15) | JSON altered by 1e-7: FAIL |
| K6 | `p`, extrapolation and error recomputed from the JSON levels with the spec formulas (exact match) | Using levels 1,2,3: FAIL |
| K7 | V3 relative difference below 1e-3 | A bottom offset by 2e-3: FAIL |
| F1 | Room B's free-node w-gradient on fiber-edge triangles is zero to roundoff (`|g_w|·h_w < 1e-12`, leak below 1e-12 of the K diagonal) | Gradient + 1e-6: FAIL. A first version tested only the leak ratio, and its planted twin did **not** fire (5.5e-13 < 1e-12). I restated the criterion and reran |
| G1 | My own finite-difference check, built from the embedding in ℝ⁴: `A''(0) = Q(φ)` for an admissible test φ. Rel 3.0e-8 at W=7/5 and 2.5e-8 at W=17/10 | Q without `−2φ²`: rel 0.59: FAIL |
| G2 | Room A's T3 step-2 identity, symbolically | `−νg` term dropped: FAIL |
| G3 | Gegenbauer sector eigenfunctions, `λ=(ν+m)(ν+m+1)`, at exact and numerical ν | `(ν+m)(ν+m+2)`: FAIL |
| G4 | Indicial roots ±ν. Without the weight the roots are `(1±√(1+4ν²))/2`, so the limit-circle threshold is `ν<√3/2` | 3D-radial operator: roots not ±ν: FAIL |
| G5 | The lowest six equal my stage-1 closed form (`results_C.json`) and both rooms' values (rel < 1e-12). T5 sign, index and nullity decided exactly, and equal to both rooms' | `ν_k = kπ/W`: FAIL |
| G6 | At W=π/2 the edges `w=±W` are the same set in ℝP³ (Room B §0 remark): distance 5.6e-16 | W=1.4: distance 0.34: FAIL |

My own scripts had bugs during development, all of them in the harness and none in the rooms. C1 at first had no absolute floor: run-to-run noise of 1e-14 on the derived `extrap − exact` values failed it. Two parsers misread Room B's `π/2` headings and its `(k,m)` commas. `s2_regrade.py` lacked `__file__`. All of these were fixed, and the scripts were rerun.

## V1: validity

### Room A: **VALID**

- **Explicit requirements: all met.** I read them from `t4_fem.py` and confirmed them with K1, K2 and K6:
  - P1 elements in `(y,w)` on `[0,π]×[−W,W]`;
  - the seam maps node `(π,w_j)` to dof `(0,−w_j)` with sign −1;
  - weight `|cos y|`, with stiffness `|cos y|φ_y² + φ_w²/|cos y|`;
  - Dirichlet data on `w=±W` and on the fiber;
  - grading `(π/2)(k/16)²`, 8 cells across, and exactly four midpoint refinements;
  - `p` from levels 2,3,4, Richardson from 3,4 with that `p`, `err = |λ_extrap − λ_finest|`;
  - all seven widths, six eigenvalues per level.

  Section 3's geometry, normal `e₄`, seam sign, admissible class and data are used as stated.
- **Readings.** R2 reads "uniform refinement" as midpoint insertion without regrading. That is the plain meaning of refining a given mesh uniformly, and Room B made the same choice. The regraded alternative moves no bottom by more than 3.1e-7 relative (`s2_regrade.json`). R1 (seam smoothness in T2) does not touch any explicit requirement, and the form domain is the same under either reading.
- **Sources.** No sign of any source beyond the spec sheet. I scanned the scripts, logs and manifest: no outside paths, no network, and nothing opened beyond the room's own outputs.
- **Instrument defects.** None that are material:
  - the rerun reproduces the JSON (C1);
  - the report reproduces the JSON (C2, C3);
  - the solver is converged at every level (K3, K4).

  Noted as immaterial: C4.5 is only a 5% gross check. It can fail; its periodic-seam plant gave 264.85 vs 2.38.

### Room B: **VALID**

- **Explicit requirements: all met**, as for Room A (K1, K2, K6). The seam is `(π,−w_j) → −(0,w_j)` (`dof_map`).
- **Sources.** No sign of any.
- **Instrument defects.** None material. Two are noted as immaterial:
  1. The report and a script comment say the free node's w-gradient on fiber-edge triangles is "exactly zero". It is computed through `np.linalg.inv`, so it is zero only to roundoff: `|g_w|·h_w ≤ 1.1e-16`, and the leaked stiffness is at most 1.6e-31 of the smallest K diagonal (F1). Its effect: Room B's eigenvalues equal Room A's to at most 8.2e-14 relative at every width and level, and Room A drops the term exactly.
  2. `--plant` switches on two defects at once (periodic seam and signed `cos`), and the T4 checks run only on the level-0 mesh at three widths. Each check was still shown to fail. I also planted the seam sign alone through Room B's own dof map (K2: residual 2.0).

Rooms A and B, written independently, agree at every width and level to at most 8.2e-14 relative, and their extrapolations to at most 2.1e-13.

## V3: comparison with my stage-1 bottoms

The rooms' values are identical to the digits shown; Room A's W=3/2 extrapolation ends …3411 rather than …3410. "|Δ|/err" is the absolute difference over the room's own error estimate.

| W/R | room bottom (A and B) | `p` | room err | my stage-1 bottom (err) | rel. diff | \|Δ\|/err | > 1e-3? |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1/4 | 45.761589281931 | 1.996684 | 5.599e-03 | 45.761602911538 (4.6e-11) | −2.98e-07 | 2.4e-03 | no |
| 1/2 | 13.011196550217 | 1.999420 | 1.184e-03 | 13.011197054679 (1.3e-11) | −3.88e-08 | 4.3e-04 | no |
| 1 | 4.038197480573 | 2.000226 | 3.190e-04 | 4.038197427067 (4.0e-12) | +1.32e-08 | 1.7e-04 | no |
| 7/5 | 2.380875542763 | 2.000371 | 1.842e-04 | 2.380875488666 (2.4e-12) | +2.27e-08 | 2.9e-04 | no |
| 3/2 | 2.143820313411 (A), …410 (B) | 2.000314 | 1.658e-04 | 2.143820262429 (2.1e-12) | +2.38e-08 | 3.1e-04 | no |
| π/2 | 2.000000050300 | 2.000203 | 1.548e-04 | 2.000000000000 (2.0e-12) | +2.52e-08 | 3.2e-04 | no |
| 17/10 | 1.777769907611 | 1.999706 | 1.382e-04 | 1.777769846306 (1.8e-12) | +3.45e-08 | 4.4e-04 | no |

All 14 comparisons (2 rooms × 7 widths) fall below 1e-3; the largest is 3.0e-7, at W=1/4. **There is no defect to locate, either in the returns or in my stage-1 computation.** The rooms' error estimates bound their actual error by 2–4 orders of magnitude, as both reports themselves say. Their finest-level λ₁…λ₆ lie above my values by 7.7e-5 to 8.7e-4 relative (`s2_audit.json`), consistent with a conforming P1 method.

## G: grading

These are graded as **audited arguments**, not as verified or proven theorems.

### Room A

**T1.**

| Step | Grade |
| --- | --- |
| 1 Varied map `cos(tφ/R)X + R sin(tφ/R)e₄` | ESTABLISHED |
| 2 Well defined on ℝP³; edge and cone point fixed | ESTABLISHED |
| 3 Exact metric `g_t = c²g + t²dφ⊗dφ` | ESTABLISHED |
| 4 Exact area via the determinant lemma, for `|t| < πR/(2 max|φ|)` | ESTABLISHED |
| 5 Differentiation under the integral, giving `A''(0) = Q` | ESTABLISHED |
| No boundary or vertex term | ESTABLISHED |
| Extension by continuity (estimate plus closability) | ESTABLISHED |

It establishes `Q(φ) = ∫(|∇φ|² − 2φ²/R²) dA`, with `dA = |cos(y/R)| dy dw`. No edge, seam, corner or cone-point term appears, and the form extends uniquely to the admissible class. G1 confirms the formula independently (rel 3e-8).

**T2.**

| Step | Grade |
| --- | --- |
| Reading R1: the literal value-only domain is not symmetric; true, as checked by Green's formula | ESTABLISHED |
| 1 Sector reduction, closure = ⊕ of the sector closures | ESTABLISHED |
| 2 Endpoints = the two poles | ESTABLISHED |
| 3 Frobenius exponents ±ν; LC iff ν<1 | ESTABLISHED |
| 4 Classification | ESTABLISHED |
| 5 Friedrichs realization selected | ESTABLISHED |

Remark on step 3, outside the argument: the italic aside says that if the weight `cos ψ` were dropped, "the threshold would move to ν < 1/2". That holds only if the measure alone changes. If the area element is dropped from both the operator and the measure, the threshold is `ν < √3/2` (G4). It carries no weight in the derivation, so the step stands.

It establishes:
- **0 < W ≤ πR/2:** every sector is limit-point (ν_k ≥ 1; the borderline ν=1 at W=π/2, k=1 is limit-point). The closure is self-adjoint, with no other extensions.
- **W/R = 17/10:** sector k=1 (ν₁ = 5π/17 ≈ 0.924) is limit-circle at both ends, and k ≥ 2 are limit-point. The deficiency indices are (2,2), with a U(2) family of extensions. The admissible class selects the Friedrichs extension: separated conditions with no r^{−ν} part at either end.

**T3.**

| Step | Grade |
| --- | --- |
| 1 Legendre form | ESTABLISHED |
| 2 Gegenbauer, `λ = (ν+n)(ν+n+1)` (identity confirmed by G2) | ESTABLISHED |
| 3 `f_{k,n}e_k ∈ D(A)` | **ESTABLISHED, SUPPLIED** |
| 4 Completeness, spectrum exactly `{λ_{k,n}}` | ESTABLISHED |
| Lowest `π(πR+2W)/(4W²R)` | ESTABLISHED |

Supplied part for step 3: the step tests against sector functions with the sector form, but the conclusion concerns the 2D Friedrichs operator, which needs `A = ⊕_k A_k`. For `u` in the core, integrate twice by parts in θ (`u = e_k = 0` at θ = ±W) and apply Parseval. This gives `a(u) = Σ_k a_k(u_k)`, so the closure is `H = {u : u_k ∈ H_k, Σ a_k(u_k) < ∞}`, and its operator is `⊕A_k`. Hence `a(f e_k, g) = a_k(f, g_k) = λ⟨f, g_k⟩ = λ⟨f e_k, g⟩` for all `g ∈ H`.

It establishes the spectrum `{(ν_k+n)(ν_k+n+1)/R²}` with `ν_k = kπR/(2W)`, labels (k,n), and multiplicity equal to the number of labels. The list is complete (eigenbasis, compact resolvent), and the lowest eigenvalue is (1,0). G5 checks agreement with my stage-1 closed form.

**T5.**

| Step | Grade |
| --- | --- |
| 1 `J = A − 2/R²` | ESTABLISHED |
| 2 Lowest eigenvalue `(ν₁−1)(ν₁+2)/R²`, sign = sign(πR/2 − W) | ESTABLISHED |
| 3 Index and nullity by counting, via min-max | ESTABLISHED |

What it establishes is in the T5 table below.

### Room B

**T1.**

| Step | Grade |
| --- | --- |
| 1 Geodesics | ESTABLISHED |
| 2 Varied surface well defined; edge and cone point fixed | ESTABLISHED |
| 3 `g_t = cos²u·g + du⊗du` | ESTABLISHED |
| 4 Exact density | ESTABLISHED |
| 5 Differentiation, `A''(0) = Q` | ESTABLISHED |
| No boundary, seam or vertex terms | ESTABLISHED |
| Extension, with `|B(φ,χ)| ≤ 2‖φ‖‖χ‖` | ESTABLISHED |

It establishes the same as Room A's T1.

**T2.**

| Step | Grade |
| --- | --- |
| Reading D₀; the form closure is unchanged | ESTABLISHED |
| Sectors | ESTABLISHED |
| (a)–(c) Deficiency indices split over sectors | ESTABLISHED |
| Endpoints: `tan^{±ν}(r/2)`, exponents ±ν, LC iff ν<1 | ESTABLISHED |
| Classification | ESTABLISHED |
| Friedrichs selected | ESTABLISHED |

It establishes the same as Room A's T2. Room B's §0 geometric remarks are correct but unused: G6 confirms that at W=π/2 the two edges are the same set in ℝP³.

**T3.**

| Step | Grade |
| --- | --- |
| 1 Gegenbauer | ESTABLISHED |
| 2 Form-domain membership via cutoffs | ESTABLISHED |
| 3 `F ∈ D(A_F)` via the 2D Green formula on D₀ | ESTABLISHED |
| 4 Completeness in `L²(L_W)` | ESTABLISHED |
| 5 Spectrum discrete and exact | ESTABLISHED |

It establishes the same spectrum as Room A, with labels (k,m), complete.

**T5.**

| Step | Grade |
| --- | --- |
| `J = A_F − 2`, eigenvalues `(ν_k+m−1)(ν_k+m+2)` | ESTABLISHED |
| Min-max | ESTABLISHED |
| Sign analysis | ESTABLISHED |
| Exact table | ESTABLISHED |

### T5 as established by both arguments (and recomputed exactly in G5)

| W/R | lowest Jacobi eigenvalue | sign | negative eigenvalues | nullity |
| --- | --- | --- | --- | --- |
| 1/4 | +43.761602911537 | + | 0 | 0 |
| 1/2 | +11.011197054679 | + | 0 | 0 |
| 1 | +2.038197427067 | + | 0 | 0 |
| 7/5 | +0.380875488666 | + | 0 | 0 |
| 3/2 | +0.143820262429 | + | 0 | 0 |
| π/2 | 0 | 0 | 0 | 1 |
| 17/10 | −0.222230153694 | − | 1 | 0 |

The index of 1 at 17/10 depends on the Friedrichs realization that the admissible class selects. Both rooms state this dependence.

**Replacement rooms:** none are needed, since both returns are VALID.

## Consulted-material manifest

**Files read in this room at stage 2:**
- the briefs: `BRIEF.md`, `auditor_brief.md`, `spec_sheet.md`;
- my stage-1 files: `STAGE1.md`, plus `stage1.json` and `results_C.json` read through scripts;
- Room A: `RETURN.md`, `t4_fem.py`, `run_all.py`, `t1_second_variation.py`, `t3_spectrum.py`, `make_tables.py`, `run_all.log`, and `results.json` through scripts;
- Room B: `RETURN.md`, `t4_fem.py`, `run_all.py`, `collect.py`, `t1_t3_symbolic.py`, `t3_t5_exact.py`, `out_symbolic.txt`, `out_symbolic_planted.txt`, `run_all_log.txt`, `t4_log.txt`, and `results.json` through scripts;
- the rest of `room_a/` and `room_b/` (`t4_tables.md`, `t3_out.json`, `t4_out.json`, `exact.json`, `t4.json`) were read only by a pattern-scan script that looked for outside paths and network use. Only the matching lines were displayed.

**Files I created and read back:** everything under `rerun_a/`, `rerun_b/`, `lib_a/` and `lib_b/`; the `s2_*` scripts, logs and JSON; `make_stage2.py`; `stage2.json`.

**Loaded without my asking:**
- the session context, which supplied the user's account e-mail address (not used);
- directory listings showing `py` and `.room.sb` (neither opened);
- background-task notifications that named an output file outside the room (not read; I used my own logs inside the room);
- one Python traceback that printed the interpreter's library path outside the room (nothing opened there);
- harness notes that my own scripts had changed on disk.

Nothing outside the room was read. There was no network use.

**Results used from memory rather than derived here:**
- the exponential map of the round sphere;
- the matrix-determinant lemma;
- the general second-variation formula with `Ric = 2/R²`;
- Weyl's limit-point/limit-circle alternative;
- von Neumann deficiency theory and the limit-circle decomposition of `D(T*)`;
- Kato's representation theorem and the Friedrichs extension;
- the Frobenius method;
- the Gegenbauer equation and its weight `(1−x²)^{α−1/2}`;
- Weierstrass density;
- the min-max principle;
- Sylvester's law of inertia and the block structure of the Bunch–Kaufman LDLᵀ;
- shift-invert Lanczos (ARPACK);
- Richardson extrapolation;
- Gauss–Legendre quadrature and central differences;
- conforming-FEM convergence from above.
