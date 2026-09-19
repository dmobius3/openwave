# M8.11 derivation log

Local. F1 and F2 below are the author's two review units; the letters name units, not findings. The fresh derivation was authorized on 2026-09-13 by both units, after revision 3 of the author-side M8.11 outline, which is not filed. Nothing is pushed. Steps 1 to 3 are done, and both units closed the derivation after step 3. The next surface is the pre-registration; there is no fourth calculation round (F1).

## Files

- `m810_core.py` and `m810_exact.py`: byte-identical copies of the files pinned in #546 (sha256 `9339d871…` and `8436ef9f…`), read-only.
- `m811_ops.py`: section algebra for the float route, reusing M8.10's product, fibre and norm routines, plus F1's support gate.
- `m811_pyramid.py`: step 1, the normalization bridge and the pyramid; 42 checks.
- `m811_prism.py`: step 2, the prism, with the S lemma first; 58 checks.
- `m811_exact.py`: step 3, the independent high-precision route, 80 digits by default, sharing no code with `m811_ops.py`; 80 checks. The file name is historical: the route identifies algebraic values, it does not evaluate them symbolically.
- `out/`:
  - `pyramid.json`, `prism.json` and `exact.json`;
  - the clean-directory logs `step1_pyramid_clean.log`, `step2_prism_clean.log` and `step3_exact_clean.log`;
  - `step3_exact_firstpass_sqrt230_failed.log`, the record of the failed field prediction.

All three scripts run in order from a clean directory, a fresh folder holding only the six sources with no `out/`. Every check passes, every arm fires, and each exits 0; F1 reproduced this. F1 also reran step 3 at 100 digits and `out/exact.json` came out byte-identical. F1 found that step 1's first version crashed on a clean checkout, because my folder already had `out/`; it now creates the directory.

## What "exact" means here (F1's redline)

Step 3 identifies each value with `Fraction(...).limit_denominator(10^25)` and accepts it only if the residual is below 10^-68 at 80 digits. Two distinct rationals with denominators up to 10^25 differ by at least about 10^-50, so each identification is unique within that class, and the 100-digit rerun reproduced every one. These are very strong candidate exact values. They are not symbolic proofs, and the pre-registration freezes them as candidate exact, for the blind agents' derivations to promote.

Claims stay exact only where an analytic argument exists:
- the 1/4 bridge;
- the existence theorem;
- v_y = 0 from S;
- λ₄'s independence of the tilt;
- the support selection;
- the factorized form, once its Clebsch–Gordan identity is written out as an identity rather than inferred from sample points.

The step-3 checks at 1e-162 on the forbidden levels corroborate the representation-theoretic zero; they do not prove it.

## Step 1: the normalization bridge and the pyramid

The tilt operator is Π_T DN_Φ − Q = (1/4)·Hess_T Q_σ = (w₆(σ)/4)·Hess_T r̂₆, in the real Hilbert structure. That identity is analytic. It was checked by three routes through the shared float algebra: first order at a non-critical point, second order at the pyramid, and at the hexagon along random tangents. They share code and are not independent (F1). Dropping DN's conjugate term gives 1/8. F2's 1/2 is the Wirtinger convention.

On the first pass, a parent gate caught a section-versus-fibre mismatch in my own comparison with M8.10's audited residual: the ratio was exactly √(d/7) in both sectors. The comparison now runs in the audit's fibre basis. #547's solver and auditor used the two conventions (5.54e-4 and 8.46e-4, the same number), so every frozen number states its basis.

## Standing of #547's R5 values, from the filed as-run worklist

- Items 4, 7, 8 and 9 were asked at R1 to R5; that was the designer's go-time change.
- Only item 6, the residual, was graded, as D2.

So the residual is a parent, while Q, the ξ level norms and λ₄ at R5 are blind, asked, ungraded diagnostics. Reproducing them is a consistency check that adds standing.

## Step 2: the prism

**The S lemma is green in all three steps:** the fibre side; the bundle level by two routes (the projector identity, and a real model of σ checked on all 120 elements); and Sξ = ξ with a vanishing S-odd forcing. v_y = 0 rests on two consequences of S, both S-forced and both checked (F2):
- the forcing's S-odd component vanishes;
- the cross term H_xy vanishes, because the second-variation form is S-even while τ_x and τ_y have opposite S-parity.

**Provenance of the second variation (F1).** Q and the chart restriction ‖ρ₆‖²(x, y) are paper parents. The in-locus second-variation values are exact derivatives of that restriction, converted by w₆/4 and the Fubini–Study metric 200/1849. They were exposed in the feasibility check and are rederived here; the paper states none of them.

## Step 3: the independent high-precision route

**Algorithm.** It does not use M8.10's cubing through three-fold coupling channels.
- The mixed density of two block sections is c₀ + d₁₂. Here c₀ = (d/7)⟨a, b⟩, and d₁₂ is one level-12 matrix coefficient X_d R⁶(g) y_d(a, b). Lemma 4.1, that the density has no multipoles K = 1 to 5, is checked rather than assumed.
- A level-12 density times a block section has, at level 2J, the factorized form X_J R^J(g) Y_J. X_J = (X_d ⊗ P) C_J^T is fixed by the sector, and Y_J = C_J(y_d ⊗ c) lives on the left index.
- So everything reduces to left-index Clebsch–Gordan algebra times the sector constants C_J = tr(X_J†X_J)/(2J+1).

The factorized cubic is checked pointwise at two group elements to 1e-81, and dropping the W twist in y_d makes that fail.

**Shared, disclosed:** `m810_core`'s exact Racah coefficients, and `m810_exact`'s mpmath group, D³ and isotypic projectors. These are the pinned M8.10 primitives; nothing else is shared. The float JSON is read only at the end, as a posterior cross-check.

**Candidate exact sector constants, identified at 80 digits and reproduced at 100.** They are the post-derivation note's C_{σ,n}.

| | sector 3′ | sector 4 |
| --- | --- | --- |
| C by level | 10: 12/91; 14: 18/221; 16: 6/119; 18: 240/4199 | 8: 12/91; 12: 12/91; 14: 6/119; 16: 18/221; 18: 2196/29393 |
| forbidden levels | 8, 12: C below 1e-160 | 10: C below 1e-160 |

This identifies F2's full candidate family, including 2196/29393 = 2196/(7·13·17·19) at the multiplicity-two slot. As a by-product the route derives the factorization the M8.10 note left as an observation; that stays author-side and unscored.

**The pyramid, independently reproduced at high precision with the same candidate exact forms.** Q, the second variation along e_t, the forcing, the tilt, the ξ level norms and λ₄ all come out equal to the step-1 table. That table's entries are the audited residual, the inherited tilt forms and #547's blind diagnostics. So the author now has an independent route to the pyramid tilt. Its exact form still rests on #547's symbolic forcing and the paper's exact second variation; this route confirms it to 80 digits, which is identification, not symbolic proof.

**The prism: candidate exact values from the independent high-precision route.** These values are new. Conventions: g = 1, unit sections, τ_x the S-even chart tangent oriented by increasing x, τ_y = iτ_x, at the representative z = +√(23/10).

| | sector 3′ | sector 4 |
| --- | --- | --- |
| forcing along τ_x | −120576113√115/1288994436873 | 214430223√115/3055394220736 |
| tilt v_x, per g (v_y = 0 by S) | 17225159√115/64285514280 ≈ +2.87342e-3 | −10210963√115/28571339680 ≈ −3.83253e-3 |
| ‖v‖²/g² | 296706102575281/35935889967339860160 | 104263765387369/7098447400956021760 |
| λ₄/g² | −336158940460/150812349114141 | −11093213145/4965015608696 |
| ‖Πₙξ‖²/g² | 10: 319550/88158077163; 14: 12920075/7517877885232; 16: 16583/68316230736; 18: 2355742375/6059914391677302 | 8: 7875/859947712; 12: 3521/1587595776; 14: 16611525/37011091127296; 16: 116081/701717332992; 18: 4105722425/19152322028017152 |

The float pipeline agrees with every entry to 1e-12, and F2 checked the table's internal consistency in exact arithmetic: forcing = −H_xx·v_x, and ‖v‖² = 115 × (the rational in v_x)².

**Recorded: the field prediction failed on the first run and was corrected.** The stated prediction, made before the run, was a rational times √230 at the prism; it was F2's suggested field, which I adopted. It failed in both sectors: the best rational fits had 25-digit denominators. The parity argument itself held, F(−x₀) = −F(x₀). The error was the unit tangent's norm √(4 + 2x₀²) = √2·√(2 + x₀²), which was dropped. With it, F = (x₀/√2) times a rational function of x₀², and x₀/√2 = √115/10. The corrected field, a rational times √115, identifies with 9- to 13-digit denominators, and √230 is kept as an arm that must fail. The pyramid's prediction, a rational times √39, held as stated. The failed log is kept in the package (F1).

## One failure mode, three instances (F2)

Three of this program's catches share one failure mode, the normalization of a basis vector:
- the section-versus-fibre factor √(d/7) in step 1;
- the unit tangent's √2 in step 3;
- the w₆ and Fubini–Study conversions, which the two-route gates had to carry.

The rule for the pre-registration: every frozen number names its basis in its own row (unit section or unit fibre vector, and which tangent), and a field is predicted with every normalization included.

## Observations, author-side and unscored

- **The prism's λ₄ sector ratio is close to 1 but distinct (F2).** λ₄(3′)/λ₄(4) = 10976618464/11002656303 ≈ 0.99763. It differs from 1 by 26037839/11002656303 ≈ 2.37e-3, and the two fractions differ throughout, not in a late digit. It is not an identity. Across the six rays the ratios run 0.72 to 7.79, so this near-miss strengthens the finding that no single multiplier exists.
- **Why the sectors tilt in opposite directions.** On the levels both sectors share (14, 16, 18), the forcing's per-level contributions have the same sign in both sectors, at both rays. The sign of each total is set by the nearest level, which the sectors do not share: 10 in 3′ and 8 in 4.
- **The heuristic's domain (F2).** The nearest-level heuristic governs signed first-order quantities, such as the forcing. It fails for sums of squares, such as λ₄, where every term has one sign and the numerators decide.

## Next: the pre-registration (F2's two requirements)

1. Every frozen table states its normalization in each row.
2. Each row says whether it is new or reproduced, ray by ray:
   - at the pyramid, the new content is the tilt, the existence statement and the λ₄ lemma, while Q, the level norms, the residual and λ₄ are reproductions (one graded, three ungraded diagnostics);
   - at the prism, everything is new.
