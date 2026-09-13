# AUDIT, stage 1 (independent implementation)

All values below come from scripts in this room: exact values from `algebra.py` (sympy rationals and radicals, no numerical identification needed anywhere), cross-checked by `functional.py` (a CG-free route with explicit functions and exact cubature), in two quaternion identifications. Every value is also in `audit_results.json`. To reproduce everything: `./run_all.sh`.

## Route in one paragraph

Γ is built by exact closure of `{q1, q2}` over `Q(sqrt5)`. `D^j` is realized on binary forms (`v_m = x^(j+m) y^(j-m)/sqrt((j+m)!(j-m)!)`), and CG coefficients come from my own Racah formula, checked against a highest-weight/lowering construction and against the product formula `D^j1 D^j2 = sum CG CG D^J`. The sector projectors are exact: the class sum of the 12-element class `2w = -(1+sqrt5)/2` lies in the 2-dimensional commutant, its eigenvalues are `-3` and `2 - 2 sqrt5`, and its eigenspaces have dimensions 4 and 3. Sections are handled through `psi -> psi eta^dagger`, which is isometric and commutes with `N`, so only `P = eta eta^dagger` is needed. Products of `D` functions factor into (left coupling) times (right coupling), so every quantity is a finite sum over coupling paths of (fibre-only left vector) times (sector-only right matrix).

## Item 0

| Quantity | Value |
| --- | --- |
| `\|q1\|^2`, `\|q2\|^2` | exactly 1 and 1 (exact `Q(sqrt5)` arithmetic) |
| order of Γ | 120 |
| derived subgroup | order 120 (closure of all commutators; 120 distinct commutators), so Γ equals its derived subgroup |
| element orders (derived) | 1:1, 2:1, 3:20, 4:30, 5:24, 6:20, 10:24; 9 classes of sizes 1, 12, 20, 12, 30, 12, 20, 12, 1 |
| `<3 3; 3 -3 \| 6 0>` | `sqrt(231)/462 = 1/(2 sqrt231) = 0.032897584747988449419...` (Racah, equal to the lowering construction) |

## Item 1: `dim Hom_Γ(σ, V_(n/2))`, n = 0..18

| n | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 17 | 18 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 3-dim | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 1 | 0 | 1 | 0 | 1 |
| 4-dim | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 1 | 0 | 1 | 0 | 1 | 0 | 2 |

**DERIVED.** Method: closure of the generators, classes by conjugation, `chi_j(h) = U_(2j)(Re h)` (Chebyshev, the trace over weights, checked against traces of the explicit `D^j(h)` on all 120 elements to 1.6e-14), exact sector characters `tr(P_sigma D^3(h))`. The 3-dim sector's character on the classes is `3, (1-sqrt5)/2, 0, (1+sqrt5)/2, -1, (1+sqrt5)/2, 0, (1-sqrt5)/2, 3`; the 4-dim sector's is `4, -1, 1, -1, 0, -1, 1, -1, 4`; each satisfies `<chi,chi> = 1` exactly. The dims are exact character inner products in `Q(sqrt5)`, asserted to be non-negative integers, and they agree with numeric averages of explicit matrices in both identifications. Before computing, I recognized Γ from prior knowledge as the binary icosahedral group (so `X` is the Poincaré homology sphere). No table was used and no value above depends on that recognition.

## Item 2: `Pi_6 N(Phi_(sigma,v)) = Phi_(sigma, w(v))`

| Sector | `w(v)` exactly |
| --- | --- |
| 3-dim | `w(v) = -(3/sqrt7) M_0(v) - (4/sqrt91) M_6(v)` |
| 4-dim | `w(v) = -(4/sqrt7) M_0(v) - (3/sqrt91) M_6(v)` |

Both are `w = -(d/sqrt7) M_0 - ((7-d)/sqrt91) M_6`. The coefficients of `M_1 .. M_5` are exactly 0. In general `w = sum_K (-1)^K c_K M_K`, where `c_K` is the scalar by which `sum_b [[Theta eta_b (x) eta_b]_K (x) eta_a]_3` acts on the one-dimensional Hom. `sum_b Theta eta_b (x) eta_b` is Γ-invariant, and `V_K^Γ = 0` for `K = 1..5` (computed: invariants of Γ in `V_K` for K ≤ 12 occur only at K = 0, 6, 10, 12). Also `M_0(v) = -|v|^2 v / sqrt7`, so the `M_0` term is `(d/7)|v|^2 v`, i.e. `(integral |Phi|^2) Phi`. Checks: the identity `[[Theta v (x) v]_K (x) v]_3 = (-1)^K M_K(v)` holds for all K = 0..6, `w` from the `M_K` formula equals `w` from the coupling paths, and the CG-free route gives the same fibre vector.

## Item 3: `Pi_6 N(Phi) = mu Phi`, Φ normalized (`|v|^2 = 7/d`)

| Ray | 3-dim: multiple? | 3-dim multiple | 4-dim: multiple? | 4-dim multiple |
| --- | --- | --- | --- | --- |
| R1 | yes | 1288/1287 | yes | 2289/2288 |
| R2 | yes | 1687/1287 | yes | 168/143 |
| R3 | yes | 175/143 | yes | 161/143 |
| R4 | yes | 1750/1287 | yes | 2751/2288 |
| R5 | yes | 77/65 | yes | 287/260 |

The residual `w - mu v` is exactly 0 in all 10 cases. The multiple equals `lambda_2/g = integral |Phi|^4`.

## Item 4: `||Pi_n xi||^2 / g^2`

`xi = -g sum_(n != 6) Pi_n N(Phi) / (n(n+2) - 48)`, exact. All even n ≤ 18 other than 6 are listed. Sector levels are marked `*`; the other entries are levels with `Hom = 0`, where the value is 0 automatically. All odd n and all n > 18 are 0 (item 5).

3-dimensional sector (sector levels 10, 14, 16, 18):

| Ray | n=0,2,4 | n=8 | n=10* | n=12 | n=14* | n=16* | n=18* |
| --- | --- | --- | --- | --- | --- | --- | --- |
| R1 | 0 | 0 | 7/195150384 | 0 | 18375/1464486053888 | 7/2749593600 | 245/609749017488 |
| R2 | 0 | 0 | 1225/48787596 | 0 | 8575/2860324324 | **0** | 171500/139734149841 |
| R3 | 0 | 0 | **0** | 0 | 735/220024948 | **0** | 1715/2822912118 |
| R4 | 0 | 0 | 8575/780601536 | 0 | 735/5720648648 | **0** | 1516795/812998689984 |
| R5 | 0 | 0 | 77/5475600 | 0 | 281211/189112352000 | 553/1591200000 | 30233/56458242360 |

4-dimensional sector (sector levels 8, 12, 14, 16, 18):

| Ray | n=0,2,4 | n=8* | n=10 | n=12* | n=14* | n=16* | n=18* |
| --- | --- | --- | --- | --- | --- | --- | --- |
| R1 | 0 | 63/5360582656 | 0 | 7/483225600 | 23625/7209777496064 | 49/28242739200 | 427/1927108005888 |
| R2 | 0 | **0** | 0 | **0** | 11025/14081596672 | **0** | 74725/110407229504 |
| R3 | 0 | **0** | 0 | 7/1830400 | 945/1083199744 | **0** | 26901/80296166912 |
| R4 | 0 | 1575/31719424 | 0 | 189/644300800 | 945/28163193344 | **0** | 2643557/2569477341184 |
| R5 | 0 | 3591/346112000 | 0 | 63/83200000 | 361557/931014656000 | 34839/147097600000 | 2371131/8029616691200 |

Method: exact, no numerical identification. `||Pi_J f||^2 = |coefficient tensor|^2 / (2J+1)` (Schur orthogonality). Off the block only the K = 6 path survives, so `||Pi_(2J) N||^2 = l_J(v) r_J(sigma)/(2J+1)`, where `l_J(v) = |[[Theta v (x) v]_6 (x) v]_J|^2` and `r_J = |[I_6 (x) eta]_J|_F^2`. Exact sector constants: `r_J(3-dim) = 48/91, 0, 132/91, 0, 270/221, 6/7, 240/221` and `r_J(4-dim) = 36/91, 108/91, 0, 12/7, 90/119, 18/13, 2196/1547` for J = 3..9. One entry was hand-checked from the factorization: 3-dim sector, R1, J = 5 gives `(7/3)^3 (1/9009)(132/91)/11 = 28/150579`, matching the route. The factorization was then machine-checked exactly against the route for every J = 4..9, every ray and both sectors, with `l_J` from `item5b.py` and `r_J` from `algebra.py`: 0 mismatches. The same check confirmed the ratio list in item 8 and the value 16/9. The 4-dim sector's level 18 has `dim Hom = 2`, but ξ's level-18 part is rank one (a single right factor). The CG-free route reproduces every entry to a relative 1.4e-14 (primary) and 2.2e-14 (secondary identification).

## Item 5: exact arguments for every zero

| Zero | Argument (all checked exactly in `item5.py` / `item5b.py`) |
| --- | --- |
| n > 18, and odd n | `N(Phi)` is cubic in level-6 functions, so it lives on levels n ≤ 18; half-integer spins have `Hom = 0` since `-1 ∈ Γ` acts as -1 there and trivially on σ ⊂ V_3 |
| n = 0, 2, 4 | `\|Phi\|^2` has right factor `sum_b [Theta eta_b (x) eta_b]_K ∈ V_K^Γ`, which is nonzero only for K = 0, 6. The K = 0 path reaches only J = 3; the K = 6 path reaches J = 3..9, i.e. n = 6..18 |
| non-sector levels (3-dim n = 8, 12; 4-dim n = 10) | the right factor lies in `V_J (x) Hom_Γ(σ, V_J) = 0` |
| R2 = v_0 at J even (n = 8, 12, 16) | `[Theta v_0 (x) v_0]_6` has only an M = 0 component, and `<6 0; 3 0 \| J 0> = 0` exactly for J even (6+3+J odd); equivalently the O(2) stabilizer with `R_y(π)` eigenvalue -1 forces J odd |
| J = 8 (n = 16) for R2, R3, R4 | these rays satisfy `Theta v = ±v` (R2: +, R3: +, R4: -). Coupling to top spin is the product of forms (`[a x b]_6 = (sqrt231/462) f_a f_b`), and coupling one below top is the Jacobian (`[F x f]_8 = (sqrt546/6552)(F_x f_y - F_y f_x)`). Both identities were verified exactly on generic vectors, and the constants are vector-independent. So `l_8 ∝ (f^2, f)_1 = 2f(f_x f_y - f_y f_x) = 0`. R1 and R5 are not Θ-eigenvectors, and `l_8 ≠ 0` there |
| R3 at J = 4, 5 (3-dim n = 10; also 4-dim n = 8) and J = 8 | R3 ∝ `xy(x^4 + y^4)` is fixed up to a character by a binary octahedral group of order 48 (generated explicitly, 48 elements checked). `V_J` has no vector with that character for J = 4, 5, 8 (character inner products over the 48 elements: 0). By equivariance, `l_J(R3) ∈` that space, so it is 0 |
| r_J at sector levels | never zero: `r_J ≠ 0` at every sector level (list above), so each remaining zero is a zero of `l_J(v)` |

No zero is unresolved. For R4 the D6 stabilizer alone does not force the J = 8 zero (multiplicity 1), which is why the Jacobian argument is needed there.

## Item 6: block projection of the order-a^5 equation

`Pi_6 DN_Phi[xi]` is computed exactly from `(conj Phi . xi)Phi + (conj xi . Phi)Phi + |Phi|^2 xi` by coupling paths. Every right part is verified to equal `c P` exactly. Values per g (ξ is proportional to g). "Along" is `<Phi, Pi_6 DN_Phi[xi]>`; "orthogonal" is the fibre vector `u_perp` of `Phi_(sigma,u_perp)`, with `||Phi_(sigma,u_perp)||^2 = (d/7)|u_perp|^2`.

| Sector, ray | along Φ (/g) | orthogonal component (/g) |
| --- | --- | --- |
| 3-dim R1 | -608931967/36722893315680 | 0 |
| 3-dim R2 | -1871763250/229518083223 | 0 |
| 3-dim R3 | -2203040/944518861 | 0 |
| 3-dim R4 | -1921938515/459036166446 | 0 |
| 3-dim R5 | -267786421/58544557500 | **nonzero**: `u_perp = (7188839 sqrt7/31177575000) v_(-3) - (7188839 sqrt273/202654237500) v_2`, norm^2 `51679406167921/168487138365975000000 = 3.0672611968556856e-07` (per g^2) |
| 4-dim R1 | -368688201/38687492546560 | 0 |
| 4-dim R2 | -15820875/15112301776 | 0 |
| 4-dim R3 | -162530109/75561508880 | 0 |
| 4-dim R4 | -226441775133/38687492546560 | 0 |
| 4-dim R5 | -4797453339/2497901120000 | **nonzero**: `u_perp = -(19565553 sqrt21/147804800000) v_(-3) + (58696659 sqrt91/960731200000) v_2`, norm^2 `1148432592587427/2840013657395200000000 = 4.043757288268624e-07` (per g^2) |

The CG-free route gives the same orthogonal norms (3.067e-07, 4.044e-07 at R5; below 1e-18 elsewhere). **Reading:** at R5 the block equation `lambda_4 Phi = g Pi_6 DN_Phi[xi]` has no solution with Φ fixed and `xi, zeta ⊥ block`. The expansion as posed breaks down at order a^5 on R5, in both sectors, even though the order-a^3 condition (item 3) holds there. For R1 to R4 the orthogonal part vanishes exactly, consistent with their stabilizers: the relevant semi-invariant subspace of `V_3` is one-dimensional there. For R5 it is two-dimensional (C5 only).

## Item 7: `lambda_4 / g^2`

| Ray | 3-dim | 4-dim |
| --- | --- | --- |
| R1 | -608931967/36722893315680 = -1.65818080227e-05 | -368688201/38687492546560 = -9.52990686993e-06 |
| R2 | -1871763250/229518083223 = -8.15518857476e-03 | -15820875/15112301776 = -1.04688718069e-03 |
| R3 | -2203040/944518861 = -2.33244680542e-03 | -162530109/75561508880 = -2.15096431251e-03 |
| R4 | -1921938515/459036166446 = -4.18689997758e-03 | -226441775133/38687492546560 = -5.85310032333e-03 |
| R5 | -267786421/58544557500 = -4.57406174776e-03 | -4797453339/2497901120000 = -1.92059377394e-03 |

Method: `lambda_4 = g <Phi, DN_Phi[xi]>` from the Φ-component of the block projection (item 6), exact. Sign argument: `<Phi, DN_Phi[xi]> = 2<N(Phi), xi> + conj<N(Phi), xi>`. With `xi = -g sum Pi_n N / mu_n`, `<N, xi>` is real, so `lambda_4 = 3g<N(Phi), xi> = -3 g^2 sum_(n=8..18) ||Pi_n N(Phi)||^2 / (n(n+2) - 48)`. Every level present has `n ≥ 8`, so `mu_n ≥ 32 > 0`, and `||Pi_18 N|| > 0` at every ray (`l_9 > 0`, `r_9 > 0`). Hence `lambda_4 < 0` for every real `g ≠ 0`, independent of the sign of g. The identity was verified exactly for all 10 cases, and numerically by the CG-free route. At R5 (item 6) the value is the Φ-component only, since the full block equation is inconsistent there.

## Item 8: ratios

| Ray | `lambda_4(3-dim)/lambda_4(4-dim)` | `(mu_3 - 1)/(mu_4 - 1)` |
| --- | --- | --- |
| R1 | 3181358848/1828392507 = 1.73997587269701 | 16/9 |
| R2 | 4889504/627669 = 7.78994023920251 | 16/9 |
| R3 | 3596800/3316941 = 1.08437261922959 | 16/9 |
| R4 | 803291852800/1122966354231 = 0.715330294423727 | 16/9 |
| R5 | 699523712/293721633 = 2.38158730378569 | 16/9 |

For λ_4 there is **no** single c. For item 3's multiples minus 1 there **is** one: c = 16/9 at every ray. Reason: `mu - 1 = (7/d) kappa_6(sigma) <v^, M_6(v^)>` with `v^` the unit fibre vector, so the ratio is `(kappa_6(3)/3)/(kappa_6(4)/4) = ((4/sqrt91)/3)/((3/sqrt91)/4) = 16/9`, independent of the ray. λ_4 instead weights the levels by `r_J(3)/r_J(4)`, which varies with J (J=3: 4/3; J=7: 21/13; J=8: 13/21; J=9: 140/183; and it is 0 or infinite at J = 4, 5, 6).

## Item 9: order-a^3 equation as a vector identity

| Check | Result (max over rays, both sectors) |
| --- | --- |
| exact: Casimir `J^2 = Jz^2 + (J+J- + J-J+)/2` built from explicit spin-J matrices, `-Delta = 4 J^2` applied to the exact coefficient tensors of ξ (left and right index), residual of `(-Delta - 48) xi - lambda_2 Phi + g N(Phi)` | exactly 0 in all 10 cases |
| pointwise at random points, `-Delta xi` from spin-matrix Casimir, `N(Phi)` evaluated directly from the explicit section | ≤ 1.67e-14 (primary), ≤ 2.12e-14 (secondary); scale `max\|N\|` from 0.78 to 2.64 |
| pointwise, `-Delta xi` by 9-point finite differences along the geodesics `g(cos t + e sin t)`, e = i, j, k (the spin-1/2 generators through the identification), h = 1e-3 | ≤ 1.79e-10 (primary), ≤ 1.89e-10 (secondary) |

The first line alone is close to tautological (on one level `J^2` is scalar), so the finite-difference line is the independent one. It uses no level decomposition and no formula `n(n+2)`, and its mutation (48 changed to 47) fails as it should.

## Checks and mutations

168 checks, 168 PASS. Each has a recorded mutation that was run and failed (records in the `*.json` files). Some mutations only move the threshold ("add 1 to the error"; used in `check_cg_D.py` for homomorphism, unitarity and raising-operator checks, and in some item 1 numeric lines). These are weak. The substantive checks have structural mutations: sign conventions dropped, a non-group element, wrong λ_2, 48 changed to 47, a different sector's value, dims shifted by two levels, the factor-d bug reintroduced.

Defects found and fixed during the run (all before the final numbers):

| Defect | Effect | Fix |
| --- | --- | --- |
| `(-1)**(i-j)` with a negative exponent gives a Python float inside the exact Θ | floats leaked into the exact route; first run aborted | integer sign |
| sector projector had unrationalized `1/(a + b sqrt5)` entries | exact zero tests failed | `radsimp` |
| functional route extracted the fibre vector without dividing by `sum \|eta\|^2 = d` | multiples d times too large, item 9 residual O(1) | divide by d; this bug is now a mutation |
| invalid mutations (identity element, only even K, Schur at J=10 beyond the cubature degree, sympy failing to simplify a constant ratio) | false FAIL or uncaught mutant | replaced (see scripts) |

## Readings taken (underdetermined points)

| Point | Reading |
| --- | --- |
| quaternion to SU(2) identification | two identifications run: `w - i(x sx + y sy + z sz)` and `i -> i sz, j -> i sy, k -> i sx`. All exact values are identical, and the CG-free route agrees in both |
| "every level of the sector" (item 4) | all even n ≤ 18 listed, sector levels marked; ξ has no content above 18 |
| "component orthogonal to Φ" (item 6) | reported as the fibre vector and the squared norm of the section, per g |
| λ_4 at R5 | the block equation at order a^5 is inconsistent there (item 6); λ_4 is reported as the Φ-component |

## Reproduction from a clean directory

I copied the scripts and `py` alone into a fresh subfolder `clean_repro/` and ran `./run_all.sh` there. Result: 168/168 checks PASS. Items 0, 1, 2, 3, 4, 5, 6, 7, 8 are identical to the room run (JSON equality), and the item 9 residuals are identical (max difference 0; the random points are seeded). The first clean attempt failed because `run_all.sh` ran the report before `item5b.py`; the order is fixed.

## CONSULTED-MATERIAL MANIFEST

| Kind | Item |
| --- | --- |
| Files read | `BRIEF.md`, `worklist.md`, `py` (room interpreter), and my own `METHOD.md`, scripts, logs and JSON outputs in this room |
| Context content disregarded | the harness context contained user/project instruction and memory text from outside the room; per the rules I did not use it |
| From memory: classical results | Racah formula for CG; Condon-Shortley phase; realization of `V_j` on binary forms; `conj D^j_mk = (-1)^(m-k) D^j_(-m,-k)`; Peter-Weyl and Schur orthogonality; Schur's lemma / commutant; Weyl character of SU(2) as Chebyshev `U_2j`; Haar measure on S^3 in Hopf coordinates (`\|a\|^2` uniform); exactness of Gauss-Legendre and trapezoid rules; `-Delta = 4 J^2` on the unit S^3 with eigenvalue n(n+2); transvectants of binary forms (top coupling = product, next = Jacobian) and uniqueness of equivariant maps; 9-point central second-difference stencil; `xy(x^4 ± y^4)` as the octahedral sextic |
| Group-theoretic facts | order 120, perfectness, classes, sector characters, Hom dims, `V_K^Γ` dims, stabilizer multiplicities: all **DERIVED** by computation from the two generators. Recognized from prior knowledge but not used as input: Γ is the binary icosahedral group (X the Poincaré homology sphere), and the invariant ring of the icosahedral group has generators in degrees 12, 20, 30 (Klein). The computed invariant spins for K ≤ 12 (0, 6, 10, 12, i.e. degrees 0, 12, 20, 24) are consistent with that |
| Item 1 declaration | **DERIVED** (exact characters from the generators, class-sum projectors, exact inner products; numerically cross-checked with explicit matrices in two identifications) |
