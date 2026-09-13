# METHOD (written before any computation; not edited afterwards, only appended with dated notes)

## Conventions fixed by me

| Choice | Value |
| --- | --- |
| Quaternion to SU(2), primary | `w + x i + y j + z k -> w I - i(x sx + y sy + z sz) = [[w - i z, -y - i x], [y - i x, w + i z]]` |
| Quaternion to SU(2), secondary (check only) | `i -> i sz, j -> i sy, k -> i sx`, i.e. `[[w + i x, y + i z], [-y + i z, w - i x]]` |
| `V_j` realisation | binary forms of degree `2j`, `v_m = x^(j+m) y^(j-m) / sqrt((j+m)!(j-m)!)`, `U = [[a,b],[c,d]]` acts by `P(x,y) -> P(ax+cy, bx+dy)`; `D^(1/2)(U) = U`; `J+ v_m = sqrt((j-m)(j+m+1)) v_(m+1)` (Condon-Shortley phase), checked in code |
| CG | my own Racah-formula implementation in exact arithmetic (sympy `Rational`, `sqrt` only); second construction by highest-weight vector plus lowering operator; third check: the product formula `D^j1 D^j2 = sum CG CG D^J` at random group elements |

## Route per item

| Item | Route |
| --- | --- |
| 0 | Exact quaternion arithmetic over `Q(sqrt5)` (own class, `Fraction` pairs). Norms exactly. Closure of `{q1, q2}` under multiplication gives `Gamma`; order = size. Derived subgroup = closure of all commutators; compare sizes. `<3 3; 3 -3 \| 6 0>` from my Racah routine, cross-checked by the lowering construction |
| 1 | Characters from the generators only: `chi_j(h) = U_(2j)(Re h)` (Chebyshev, i.e. the trace over weights; also checked against the trace of the explicit `D^j(h)`). Sector projectors obtained without any table: a class sum of `Gamma` in `D^3` lies in the commutant, which is 2-dimensional here; its two eigenspaces (dims 3 and 4) are the sectors, the projector is then exact. `chi_sigma(h) = tr(P_sigma D^3(h))` exactly. `dim Hom = (1/\|Gamma\|) sum conj(chi_sigma) chi_j`, exact in `Q(sqrt5)`, checked to be a non-negative integer. Declaration will be DERIVED |
| 2 | Factorisation: products of `D` functions factor through CG as left coupling times right coupling. `N(Phi) = (conj Phi . Phi) Phi` level-3 part is `sum_K (-1)^K c_K(sigma) M_K(v)` where `c_K` is the scalar by which the right-hand contraction `sum_b [[Theta eta_b (x) eta_b]_K (x) eta_a]_3` acts on the one-dimensional `Hom`. `c_K` computed exactly from the exact projector (I replace `eta` by `P = eta eta^dagger`, which is an isometric embedding commuting with `N`). Cross-check: the explicit-function route below |
| 3 | `w(v) = mu v` test exactly per ray and sector, with `\|v\|^2 = 7/d` |
| 4 | `xi = -g sum_(n != 6) Pi_n N(Phi) / (n(n+2) - 48)`; `\|\|Pi_n xi\|\|^2/g^2 = \|\|Pi_n N\|\|^2/(n(n+2)-48)^2`, `\|\|Pi_J f\|\|^2 = \|coef tensor\|^2/(2J+1)`, exact. All sector levels `n <= 18` listed. The block equation's consistency is item 3; if it fails I still solve the off-block part and say so |
| 5 | Exact arguments: degree bound, Gamma-invariance of `\|Phi\|^2` (only `K` with `V_K^Gamma != 0` survive), stabiliser / selection-rule arguments on the left factor, and exact vanishing of the right factor where it occurs |
| 6 | `Pi_6 DN_Phi[xi]` level-3 part computed from the three terms `(conj Phi . xi)Phi + (conj xi . Phi)Phi + \|Phi\|^2 xi`, each by left/right coupling paths, exact; fibre vector `u`; component along `Phi` and the orthogonal fibre vector with its norm |
| 7 | `lambda_4 = g <Phi, DN_Phi[xi]>` (from the block projection); expected identity `= 3 g <N(Phi), xi> = -3 g^2 sum \|\|Pi_n N\|\|^2 / mu_n`, both computed; sign argument from `mu_n > 0` for all levels present |
| 8 | Ratios from item 7 and item 3 exact |
| 9 | Spin matrices `Jz, J+, J-` built explicitly, `J^2`, `-Delta = 4 J^2` applied to the coefficient tensors of `xi` (on the left index and separately on the right index), residual of `(-Delta - 48) xi - lambda_2 Phi + g N(Phi)` exactly (symbolic) and numerically on explicit functions; plus an independent geodesic finite-difference Laplacian at random points in high precision |

## Independent second route (explicit functions)

Sections are built as explicit functions of a unit quaternion: `Phi(g) = v^T D^3(g) eta` with `eta` an orthonormal basis of the sector obtained by group averaging (commutant of `{D^3(h)}`), equivariance `Phi(gh) = Phi(g) sigma(h)` checked for all 120 `h`. Level projections by exact cubature on `SU(2)` (Hopf coordinates: `|a|^2` uniform, Gauss-Legendre in `|a|^2`, trapezoid in the two phases, exact for the polynomial degrees involved: degree <= 36). This route uses no CG at all, so it cross-checks items 2-8 numerically (float64). Run in both quaternion identifications.

## Exactness policy

Primary values exact (sympy radicals / rationals; no numerical identification needed if the exact route closes). If any value is only obtained numerically, I state the precision and denominator bound, and repeat at a second precision.

## Checks

Every PASS line has a recorded mutation (a deliberately wrong target or perturbed object) that is run and must FAIL.
