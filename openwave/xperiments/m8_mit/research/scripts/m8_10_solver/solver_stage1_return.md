# RETURN

Every number below comes from a script in this room. `sh run_all.sh` reproduces all of them from a clean directory (about 90 s). The same numbers are in `results.json`. Scripts:

| Script | Role |
| --- | --- |
| `common.py` | exact Q(sqrt5) quaternions, group closure, Racah CG (exact as sign*sqrt(rational)), D^j by symmetric powers, spin matrices, identification helpers |
| `item0_1.py` | items 0, 1 (exact) |
| `engine.py DPS` | items 2 to 9 by Peter-Weyl coefficients and CG products, at 60, 110 and 160 digits |
| `quadcheck.py` | independent float64 check: exact-degree Euler-angle quadrature on SU(2), no CG used |
| `exact_parts.py` | exact symbolic route: exact projector P, exact T_L, r_L, F_J, exact ray parts, exact assembly of items 3, 4, 7, 8 |
| `relations.py DPS` | item 2: relations among the M_K |
| `item5_zeros.py` | item 5: exact vanishing arguments, stabilizer multiplicities (also used for items 3 and 6) |
| `identify.py` | identifies every engine number at two precisions, cross-checks against the exact route and the quadrature, writes `results.json` |

## Checks that can fail, and how I know they can

| Check | Result | Why it can fail |
| --- | --- | --- |
| Racah CG vs sympy CG (426 coefficients) | 0 mismatches | two independent implementations; a sign or factor bug gives mismatches |
| trace D^j(q) = U_2j(Re q), all classes, n <= 18 | max dev 3.8e-40 | wrong D^j or wrong character formula gives O(1) deviations |
| D-product rule and conj rule, random g | 4e-112, 4e-111 (110 digits) | a CG-convention mismatch between my D and my CG gives O(1) errors |
| exact P: P^2 = P, [P, D'(q1)] = [P, D'(q2)] = 0, tr P = d, Hermitian, P3 + P4 = I | all True, exact rationals in Q(sqrt5)(i) | exact equality tests over Q(sqrt5)(i) |
| generator of D^J(exp(t i sigma_k)) equals 2 i J_k | dev 1.1e-69 (110 digits; finite-difference error ~t^2) | a wrong spin-matrix normalization gives O(1) |
| exact route vs 110-digit engine, all ray numbers | max dev 6.0e-110 | two different computations (exact structure formula vs full CG products) |
| float64 quadrature vs engine | max dev 6.7e-16 | quadrature uses no CG at all |
| identification at 60 and 110 digits agrees, and agrees with exact | 0 failures | a wrong identification or non-rational value fails the tolerance or disagrees between precisions |
| item 9 control: Casimir with factor 3 instead of 4 | residual 0.007 to 0.2 | shows the residual test is not vacuous |
| exact zero tests at item 5 levels, with neighbouring-level controls | zeros exact; controls nonzero | sympy returns nonzero when the sum is nonzero (controls) |

## Item 0

| Quantity | Value | Method |
| --- | --- | --- |
| norm(q1)^2, norm(q2)^2 | 1, 1 (exact) | exact arithmetic in Q(sqrt5) |
| order of Gamma | 120 | exact closure of the generators under quaternion multiplication; closure and inverses rechecked |
| derived subgroup | order 120: Gamma equals its derived subgroup (perfect) | 120 distinct commutators; closure of the commutator set |
| conjugacy classes | 9 classes, sizes 1, 1, 12, 12, 12, 12, 20, 20, 30; element orders 1, 2, 10, 5, 10, 5, 6, 3, 4; distinct Re q per class | exact conjugation |
| <3 3; 3 -3 \| 6 0> | 1/sqrt(924) = 0.03289758474798845 (positive) | Racah formula, exact square 1/924; sympy agrees |

## Item 1: dim Hom_Gamma(sigma, V_{n/2}), n = 0..18

| n | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 17 | 18 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 3-dim sector | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 1 | 0 | 1 | 0 | 1 |
| 4-dim sector | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 1 | 0 | 1 | 0 | 1 | 0 | 2 |
| trivial (used in item 2) | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 |

**Declaration: DERIVED.** Method:

1. I built Gamma exactly from the two generators.
2. Characters of V_j come from chi_j = U_2j(Re q) in exact Q(sqrt5), and I checked them numerically against traces of D^j.
3. I split V_3 into the 3- and 4-dimensional constituents numerically, as eigenspaces of a Gamma-averaged Hermitian matrix at 40 digits.
4. I read off their characters per class and matched each value to p + q*phi with integers \|p\|, \|q\| <= 10.
5. I verified the matched characters exactly: <chi,chi> = 1, chi(1) = d, chi_3 + chi_4 = chi_V3 on every class.
6. The multiplicities are exact character inner products in Q(sqrt5), each required to be a non-negative integer.

No character table was used. I recognized Gamma as the binary icosahedral group from prior knowledge (order 120, perfect), but that recognition was not an input to any number.

Sector characters per class (class order as in `results.json`): 3-dim: 3, 3, (1-sqrt5)/2, (1+sqrt5)/2, (1+sqrt5)/2, (1-sqrt5)/2, 0, 0, -1. 4-dim: 4, 4, -1, -1, -1, -1, 1, 1, 0.

## Item 2: Pi_6 N(Phi_{sigma,v}) = Phi_{sigma,w(v)}

**Closed form, both sectors, v unnormalized.**

```text
w(v) = sum_K r_K M_K(v),  r_1 = ... = r_5 = 0
3-dim:  r_0 = -3/sqrt7,  r_6 = -4/sqrt91
4-dim:  r_0 = -4/sqrt7,  r_6 = -3/sqrt91
M_0(v) = -|v|^2 v / sqrt7   (exact; checked numerically)
hence  w(v) = (d/7) |v|^2 v + r_6 M_6(v)
```

**Derivation.**

1. Write conj(Phi_b) = sum (Theta u)_m D_mk (Theta eta_b)_k.
2. Couple the three D^3 factors, (1,2) to L and then to J, on both indices. This gives Pi_J N(Phi_{sigma,u}) the coefficients M^{L,J}(u)_M (R^{L,J})_{Ka}, summed over L.
3. Here R^{L,J}_{.a} = [T_L (x) eta_a]_J and T_L = sum_b [eta_b (x) Theta eta_b]_L.
4. T_L is the spin-L part of the Gamma-invariant tensor corresponding to the projector P = eta eta^dagger. So T_L is a Gamma-invariant vector of V_L.
5. By item 1 the only invariants with L <= 6 are at L = 0 and L = 6. Hence T_1 = ... = T_5 = 0, confirmed exactly: \|T_L\|^2 = 0.
6. At J = 3, R^{L,3} lies in the 1-dimensional Hom, so R^{L,3} = r_L eta with r_L = tr(Op_3(T_L) P)/d.
7. I evaluated r_0 and r_6 exactly from the exact P, with \|T_0\|^2 = d^2/7 and \|T_6\|^2 = 12/7 in both sectors.

**Checks.**

- The full CG engine on random unnormalized v reproduces w(v) = r_0 M_0 + r_6 M_6 to 1e-110.
- The high-precision values of r_0^2 and r_6^2 identify to 9/7 and 16/91 (3-dim) and 16/7 and 9/91 (4-dim) at 60 and 110 digits.
- r_6(4-dim) = (3/4) r_6(3-dim) follows from P_4 = I - P_3.

**Non-uniqueness (reading).** The seven maps M_K span only a 4-dimensional space (rank 4 from the singular values 54.9, 22.0, 15.9, 14.1, then about 1e-110). So "through the M_K" is not unique; the expression above is the one the derivation produces, and any other differs by the three linear relations among the M_K.

Exact relation to the independent basis B_K'(v) = [[v (x) v]_K' (x) Theta v]_3, K' = 0, 2, 4, 6 ([v (x) v]_K' = 0 for odd K'). M_K = sum W_KK' B_K', with W_KK'^2 identified at 60 and 110 digits (bound 1e12) and the exact W verified on fresh vectors to 1e-110:

| M_K | B_0 | B_2 | B_4 | B_6 |
| --- | --- | --- | --- | --- |
| M_0 | 1/7 | sqrt5/7 | 3/7 | sqrt13/7 |
| M_1 | sqrt3/7 | sqrt(135/784) | sqrt(3/196) | -sqrt(351/784) |
| M_2 | sqrt5/7 | 19/84 | -sqrt(45/196) | sqrt(1625/7056) |
| M_3 | 1/sqrt7 | -sqrt(5/252) | -sqrt(1/28) | -sqrt(13/252) |
| M_4 | 3/7 | -sqrt(45/196) | 97/154 | sqrt(117/23716) |
| M_5 | sqrt11/7 | -sqrt(1375/7056) | -sqrt(289/2156) | -sqrt(13/77616) |
| M_6 | sqrt13/7 | sqrt(1625/7056) | sqrt(117/23716) | 1/924 |

## Item 3: is Pi_6 N(Phi) a multiple of Phi (Phi normalized)?

It is a multiple at all five rays in both sectors. Exact multiples kappa = <Phi, N(Phi)> = integral of \|Phi\|^4:

| Ray | 3-dim | 4-dim |
| --- | --- | --- |
| R1 | 1288/1287 | 2289/2288 |
| R2 | 1687/1287 | 168/143 |
| R3 | 175/143 | 161/143 |
| R4 | 1750/1287 | 2751/2288 |
| R5 | 77/65 | 287/260 |

**Method.**

- With the section normalized, \|v\|^2 = 7/d. Then w(v) = v + r_6 (7/d)^{3/2} M_6(v_hat), so Pi_6 N is parallel to Phi iff M_6(v_hat) is parallel to v_hat.
- For all five rays, \|M_6(v) - <v,M_6 v> v\|^2 = 0 exactly (sympy).
- kappa - 1 = r_6 (7/d) <v_hat, M_6(v_hat)>. The inner products are exactly -sqrt91/12012, -100 sqrt91/3003, -24 sqrt91/1001, -463 sqrt91/12012 and -9 sqrt91/455 for R1 to R5.
- Engine parallel residuals: 3e-111 to 1.3e-110 at 110 digits and 1e-161 to 8e-161 at 160 digits.
- Why R1 to R4 must be parallel: the construction commutes with left translations. The ray is an eigenline of its stabilizer, and the matching character occurs once in V_3:
  - R1, R2 by weight;
  - R3: 48-element binary octahedral stabilizer, multiplicity 1;
  - R4: 24-element binary dihedral stabilizer, multiplicity 1.
- For R5 the stabilizer is cyclic of order 10 and the character occurs twice in V_3. So parallelism at R5 is not forced by symmetry; it holds because sin^2 t = 12/25, which the exact evaluation confirms.
- Also kappa - 1 = \|Pi_12 \|Phi\|^2\|^2: \|Phi\|^2 has components only at n = 0 (value 1) and n = 12, as expected from its right-Gamma-invariance.

## Item 4: norm(Pi_n xi)^2 / g^2 at the sector levels

Method: xi = -g sum_{n != 6} Pi_n N(Phi) / (n(n+2) - 48). This is the unique solution of the order-a^3 equation orthogonal to the block; it exists because item 3 holds.

**Exact route.** For J != 3, Pi_J N has coefficients M^{6,J}(u) (x) R^{6,J}_sigma (the L = 0 term only reaches J = 3). So

```text
||Pi_J N||^2 = (7/d)^3 |M^{6,J}(v_hat)|^2 F_J / (2J+1),   F_J = ||R^{6,J}_sigma||^2 = tr(Op_J(T_6)^dag Op_J(T_6) P)
```

Exact F_J, n = 2J = 6..18:

| F_J | n=6 | 8 | 10 | 12 | 14 | 16 | 18 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 3-dim | 48/91 | 0 | 132/91 | 0 | 270/221 | 6/7 | 240/221 |
| 4-dim | 36/91 | 108/91 | 0 | 12/7 | 90/119 | 18/13 | 2196/1547 |

The zeros of F_J fall exactly at levels absent from the sector.

**Numeric confirmation.**

- Every value was also identified independently from the engine: best rational with denominator <= 1e20, accepted within 1e-50 at 60 digits and within 1e-100 at 110 digits.
- Both precisions gave the same rational, equal to the exact value.

3-dim sector (levels n = 6, 10, 14, 16, 18; n = 6 is 0 by construction):

| Ray | n=10 | n=14 | n=16 | n=18 |
| --- | --- | --- | --- | --- |
| R1 | 7/195150384 | 18375/1464486053888 | 7/2749593600 | 245/609749017488 |
| R2 | 1225/48787596 | 8575/2860324324 | 0 | 171500/139734149841 |
| R3 | 0 | 735/220024948 | 0 | 1715/2822912118 |
| R4 | 8575/780601536 | 735/5720648648 | 0 | 1516795/812998689984 |
| R5 | 77/5475600 | 281211/189112352000 | 553/1591200000 | 30233/56458242360 |

4-dim sector (levels n = 6, 8, 12, 14, 16, 18; n = 6 is 0 by construction):

| Ray | n=8 | n=12 | n=14 | n=16 | n=18 |
| --- | --- | --- | --- | --- | --- |
| R1 | 63/5360582656 | 7/483225600 | 23625/7209777496064 | 49/28242739200 | 427/1927108005888 |
| R2 | 0 | 0 | 11025/14081596672 | 0 | 74725/110407229504 |
| R3 | 0 | 7/1830400 | 945/1083199744 | 0 | 26901/80296166912 |
| R4 | 1575/31719424 | 189/644300800 | 945/28163193344 | 0 | 2643557/2569477341184 |
| R5 | 3591/346112000 | 63/83200000 | 361557/931014656000 | 34839/147097600000 | 2371131/8029616691200 |

Every level n <= 18 outside the sector is 0 (no sections there), and every level n > 18 is 0 because N(Phi) is cubic in D^3 entries, so it has no component above J = 9.

## Item 5: why the zeros vanish

Since F_J != 0 at every sector level, a sector level vanishes iff M^{6,J}(v) = 0. The zeros are 3-dim: R2 n=16; R3 n=10, 16; R4 n=16. 4-dim: R2 n=8, 12, 16; R3 n=8, 16; R4 n=16.

| Zero | Exact argument | Status |
| --- | --- | --- |
| R2 = v_0, J = 4, 6, 8 (n = 8, 12, 16) | rho_6(v_0) has only an M = 0 component (10 sqrt231/231). So M^{6,J}(v_0) is proportional to <6 0; 3 0 \| J 0>, which vanishes when 6 + 3 + J is odd, i.e. J even. The printed values are 0 at J = 4, 6, 8 and nonzero at J = 3, 5, 7, 9 | resolved: selection rule |
| R3 = v_2 + v_-2, J = 4, 5, 8 (n = 8, 10, 16) | R3 is an eigenline of a 48-element subgroup H of SU(2) (verified: closed, eigenline dev 5e-16). Equivariance forces D^J(h) M^{6,J}(R3) = chi(h) M^{6,J}(R3), and the chi-multiplicity in V_J is 0 at J = 4, 5, 8 (1 at J = 3, 6, 7, 9). Multiplicities are integer character averages over 48 elements, float error below 1e-15. Independently, all components vanish in exact sympy arithmetic | resolved: symmetry, plus exact evaluation |
| R4 = v_3 + v_-3, J = 8 (n = 16) | explicit cancellation: M^{6,8}(R4)_3 = <6 6;3 -3\|8 3> rho_6,6 + <6 0;3 3\|8 3> rho_6,0 = (sqrt273/546)(-1) + (-sqrt143/26)(-sqrt231/231) = 0, because sqrt(143*231) = 11 sqrt273. The M = -3 component cancels the same way, and only M = +-3 can appear | resolved: exact CG cancellation. It is not forced by symmetry: R4's 24-element stabilizer allows a 1-dimensional chi-subspace in V_8 (multiplicity 1), so this is a coefficient identity, not a selection rule |

No zero is unresolved.

## Item 6: Pi_6 DN_Phi[xi] projected on the block

DN_Phi[h] = \|Phi\|^2 h + 2 Re(Phi^dag h) Phi. The component along Phi equals g times the value below; the orthogonal component has norm \|g\| times the value below.

| Ray | 3-dim along / g | 3-dim orthogonal norm / \|g\| | 4-dim along / g | 4-dim orthogonal norm / \|g\| |
| --- | --- | --- | --- | --- |
| R1 | -608931967/36722893315680 | 0 | -368688201/38687492546560 | 0 |
| R2 | -1871763250/229518083223 | 0 | -15820875/15112301776 | 0 |
| R3 | -2203040/944518861 | 0 | -162530109/75561508880 | 0 |
| R4 | -1921938515/459036166446 | 0 | -226441775133/38687492546560 | 0 |
| R5 | -267786421/58544557500 | (7188839/2078505000)/sqrt39 = 5.538e-4 | -4797453339/2497901120000 | (58696659/14780480000)/sqrt39 = 6.359e-4 |

- **Along:** <Phi, DN_Phi[xi]> = 2<N,xi> + <xi,N> = 3<N,xi>, which is real. So the along value is -3 g sum \|\|Pi_n N\|\|^2/(n(n+2)-48), exact from item 4. The engine's direct CG computation agrees to 1e-110.
- **Orthogonal, R1 to R4:** zero exactly. The fibre vector y(v) of Pi_6 DN_Phi[xi] is left-equivariant and homogeneous of degree (3,2), so it lies in the stabilizer-character subspace of V_3, which is 1-dimensional for R1 to R4 (item 3). Numerically below 3e-162 at 160 digits.
- **Orthogonal, R5 (anomaly):** the orthogonal component is nonzero in both sectors. It points along sin t v_2 - cos t v_-3 (overlap 1 to 30 digits), which is allowed because R5's character occurs twice in V_3.
  - Norm^2: 3-dim 51679406167921/168487138365975000000, 4-dim 1148432592587427/2840013657395200000000.
  - These were identified numerically, not derived exactly: at 110 and 160 digits, denominator bound 1e40, tolerance 1e-100 and 1e-150, same rational at both. At 60 digits the bound 1e20 is too small for these denominators, so no identification was attempted there.
  - Each norm^2 times 39 is a perfect rational square (primes <= 19 in the denominator, like every exact value above), consistent with the q/sqrt39 form.
  - The engine and the independent quadrature agree (5.5e-4, 6.4e-4).
  - The sign of the coefficient along sin t v_2 - cos t v_-3 was not recorded.

## Item 7: lambda_4 / g^2

The values are the "along / g" column of item 6: lambda_4 = g <Phi, DN_Phi[xi]> = -3 g^2 sum_{n != 6} \|\|Pi_n N(Phi)\|\|^2/(n(n+2)-48).

| Ray | 3-dim | 4-dim |
| --- | --- | --- |
| R1 | -608931967/36722893315680 (-1.658e-5) | -368688201/38687492546560 (-9.530e-6) |
| R2 | -1871763250/229518083223 (-8.155e-3) | -15820875/15112301776 (-1.047e-3) |
| R3 | -2203040/944518861 (-2.332e-3) | -162530109/75561508880 (-2.151e-3) |
| R4 | -1921938515/459036166446 (-4.187e-3) | -226441775133/38687492546560 (-5.853e-3) |
| R5 | -267786421/58544557500 (-4.574e-3), see caveat | -4797453339/2497901120000 (-1.921e-3), see caveat |

**Sign: lambda_4 < 0 for every g != 0, in both sectors, at every ray.** Argument:

1. Neither sector has sections at n < 6 (item 1), and N(Phi) has no component above n = 18. So every term has n(n+2) - 48 > 0, and lambda_4/g^2 = -3 sum (non-negative)/(positive) <= 0.
2. The sum is strictly positive iff N(Phi) has an off-block part, i.e. iff integral \|Phi\|^6 > kappa^2.
3. By Cauchy-Schwarz, kappa = integral \|Phi\|^3 \|Phi\| <= (integral \|Phi\|^6)^{1/2}, with equality only if \|Phi\| is constant. \|Phi\|^2 is a real-analytic function, and it is non-constant because kappa - 1 = \|\|Pi_12 \|Phi\|^2\|\|^2 > 0 (item 3).
4. Hence lambda_4/g^2 < 0. Since g^2 > 0, the sign does not depend on the sign of g.

**Caveat (R5).** At R5 the orthogonal part of the projected order-a^5 equation is nonzero (item 6). So under the stated ansatz (xi and zeta orthogonal to the block, Phi fixed) the order-a^5 equation has no solution at R5. The R5 entry is only the Phi-component of that equation.

## Item 8: ratios, 3-dim over 4-dim

| Ray | lambda_4(3)/lambda_4(4) | (kappa-1)(3)/(kappa-1)(4) |
| --- | --- | --- |
| R1 | 3181358848/1828392507 (1.7400) | 16/9 |
| R2 | 4889504/627669 (7.7899) | 16/9 |
| R3 | 3596800/3316941 (1.0844) | 16/9 |
| R4 | 803291852800/1122966354231 (0.7153) | 16/9 |
| R5 | 699523712/293721633 (2.3816) | 16/9 |

**For lambda_4: no single c.** The ratios differ from ray to ray. The reason is structural: the sector constants F_J are not proportional between sectors (for example F_8-level: 0 vs 108/91, F_10-level: 132/91 vs 0), so each ray weights them differently.

**For the item-3 multiples minus 1: yes, c = 16/9 at every ray.** Derivation:

1. P_sigma = (d/7) I + P_sigma^(6), where P^(6) is the spin-6 part, and P_4^(6) = -P_3^(6) because P_3 + P_4 = I.
2. With \|u\|^2 = 7/d, \|Phi\|^2 = 1 + (7/d) v_hat^T D P_sigma^(6) D^dag conj(v_hat).
3. So kappa - 1 = \|\|Pi_12 \|Phi\|^2\|\|^2 scales as (7/d)^2, and the ratio is (7/3)^2/(7/4)^2 = 16/9.
4. Equivalently, r_6(3)(7/3) / (r_6(4)(7/4)) = 16/9, using r_6(4) = (3/4) r_6(3).

## Item 9: order-a^3 equation checked with an explicit Casimir

For each level J, I built Jz, J+ and J- from the Condon-Shortley formulas, set Jx = (J+ + J-)/2 and Jy = (J+ - J-)/(2i), formed C = Jx^2 + Jy^2 + Jz^2, and applied -Delta = 4C to the free (left) index of the computed xi. The factor 4 comes from the check that d/dt D^J(exp(t i sigma_k)) = 2 i J_k, where exp(t i sigma_k) are unit-speed geodesics on the radius-1 S^3 (dev 1e-69). I then evaluated the full residual \|\|(-Delta - 48) xi - lambda_2 Phi + g N(Phi)\|\| with g = 1 and lambda_2 = kappa, including the level-6 part.

| Ray | 3-dim residual (60 / 110 / 160 digits) | 4-dim residual (60 / 110 / 160) | control, factor 3 (110 digits), 3-dim / 4-dim |
| --- | --- | --- | --- |
| R1 | 6.0e-61 / 1.1e-110 / 4.8e-161 | 4.2e-61 / 5.6e-111 / 2.0e-161 | 0.0094 / 0.0072 |
| R2 | 4.7e-61 / 8.1e-111 / 4.6e-161 | 4.1e-61 / 4.4e-111 / 3.2e-161 | 0.205 / 0.089 |
| R3 | 4.8e-61 / 1.0e-110 / 7.3e-161 | 3.8e-61 / 3.3e-111 / 6.6e-161 | 0.124 / 0.110 |
| R4 | 5.7e-61 / 1.3e-110 / 7.6e-161 | 5.6e-61 / 4.4e-111 / 4.0e-161 | 0.159 / 0.170 |
| R5 | 5.6e-61 / 6.9e-111 / 3.9e-161 | 2.7e-61 / 3.5e-111 / 1.4e-161 | 0.153 / 0.102 |

The residual falls with working precision. The wrong-normalization control does not vanish, so the check is not circular.

## Readings taken (underdetermined points)

| Point | Reading |
| --- | --- |
| "every level n of the sector" (item 4) | all sector levels n <= 18 are listed; levels above 18 are 0 by degree; n = 6 is 0 by construction |
| "through the M_K" (item 2) | the M_K are linearly dependent (rank 4 of 7), so the representation is not unique. I give the derived one (only M_0, M_6) plus the exact table relating the M_K to an independent basis |
| item 6 "component orthogonal to Phi" | reported as the norm of that block section, its exact-form identification, and its fibre direction |
| items 6 and 7 scale with g | reported as along/g, orthogonal/\|g\|, lambda_4/g^2 |
| lambda_4 at R5 | reported as the Phi-component of the projected order-a^5 equation, which has no solution at R5 under the stated ansatz |
| item 8 "multiples minus 1" | ratio (kappa_3 - 1)/(kappa_4 - 1) |
| identification of quaternions with SU(2) | q = w + xi + yj + zk maps to [[w + ix, y + iz], [-y + iz, w - ix]]. Results are invariant: all reported quantities are SU(2)-invariant contractions of the sector projector, and left translations act only on the free index |

## Things that looked wrong or surprised me

- **R5 fails at order a^5.** R5 passes item 3 exactly (Pi_6 N is parallel to Phi), but Pi_6 DN_Phi[xi] has a nonzero component orthogonal to Phi in both sectors. Two independent computations (the CG engine at three precisions, and the float64 quadrature with no CG) agree on this.
- **My first exact script stopped on a wrong assumption.** It asserted that <v, M_6 v> is rational. It is actually a rational multiple of sqrt91, which cancels against r_6 = -4/sqrt91 or -3/sqrt91. I fixed the assertion; no result depended on it.
- **The R5 orthogonal norm^2 has denominators near 1e20.** This is beyond a safe bound at 60 digits, so its identification uses 110 and 160 digits only. It is the one reported value that is identified numerically but not derived exactly.

## Consulted-material manifest

| Item | Detail |
| --- | --- |
| Files read | `BRIEF.md`, `worklist.md`, and the scripts, logs and JSON outputs I wrote in this room. Nothing outside the room. My context also contained instruction and memory text I did not request (a user CLAUDE.md, a project CLAUDE.md, a memory index); I disregarded it as the rules require |
| Textbook facts used (RECOGNIZED, then checked numerically where marked) | Peter-Weyl and Schur orthogonality (integral of D conj D = delta/(2j+1)); the D-matrix product rule via CG (checked); conj D^j_mk = (-1)^(m-k) D^j_-m,-k (checked); the Racah formula for CG (checked against sympy); the symmetric-power realization of V_j in the Condon-Shortley basis; characters chi_j = U_2j(cos theta) (checked against traces); Laplacian on a group with bi-invariant metric = sum of squares of orthonormal invariant fields, giving 4 x Casimir on the unit S^3 (normalization checked); CG symmetry making <j1 0 j2 0 \| J 0> = 0 for odd j1+j2+J (confirmed by the printed exact values); Cauchy-Schwarz; real-analytic functions taking finitely many values are constant; Gauss-Legendre exactness and the Euler-angle Haar measure (quadrature check only); linear independence of square roots of square-free integers (implicit in sympy's canonical radicals) |
| Group-theoretic facts | order 120, perfectness, classes, sector characters, all item-1 multiplicities, invariants of V_L (L <= 9 only at 0 and 6), exact sector projectors: DERIVED by computation from the two generators. That Gamma is the binary icosahedral group: RECOGNIZED from prior knowledge, not used as input. The 48-element stabilizer of R3 (a rotated binary octahedral group built from Hurwitz units and (+-a+-b)/sqrt2): RECOGNIZED as a candidate from prior knowledge (the xyz harmonic), then VERIFIED by computation (closure, eigenline test, multiplicities). The stabilizers of R4 (binary dihedral, 24 elements) and R5 (cyclic, 10 elements): generators chosen by hand from weight considerations, the groups and multiplicities DERIVED by computation |
| Item 1 declaration | DERIVED |
| Software | the room interpreter `./py` with numpy, mpmath, sympy |
