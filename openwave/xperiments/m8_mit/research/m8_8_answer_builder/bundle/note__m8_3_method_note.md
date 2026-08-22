# M8.3 method note: the mass-formula reproducer

> Task: [`../tasks/m8_3_task_details.md`](../tasks/m8_3_task_details.md). Spec of
> record: [mass-spectrum.md](https://github.com/dmobius3/mode-identity-theory/blob/13a5c2b/files/spectrum/files/mass-spectrum.md)
> at repo commit `13a5c2b`, archived as Zenodo version
> [10.5281/zenodo.21652153](https://doi.org/10.5281/zenodo.21652153) (concept DOI
> 10.5281/zenodo.18603975, all versions). Status: ✅ COMPLETE (2026-07-28). No
> independence firewall applies to this task (author, normal fork → PR → review
> workflow; roadmap row); the finding below was discovered in the course of ordinary
> reproduction, not a blind protocol.

## 1. Equations first (the formula under test, transcribed from the source page)

### 1.1 The mass formula

```text
m(rho, sigma) = mu_Lambda * C_geom(rho) * (sqrt(Omega_Lambda))^(dist(rho)/30) * T^2(rho x sigma)
```

Four factors: `mu_Lambda` (vacuum energy floor, recomputed from measured Lambda; the calibration
input is Lambda's own upstream measurement, not mu_Lambda itself), `C_geom(rho)` (a geometric
phase factor), `(sqrt(Omega_Lambda))^(dist/30)` (a
hierarchy exponent from McKay graph distance), `T^2(rho x sigma)` (Reidemeister torsion of a
vacuum-twisted flat bundle on `S^3/2I`).

### 1.2 C_geom (Kostant exponents)

```text
C_geom(rho) = ( prod_e 2 sin^2(pi e / D) )^(1 / (2 dim rho))
D = 60 for integer-spin rho, 120 for half-integer rho (the Z4-stabilizer spin rule)
```

The exponents `e` are defined by the generating function

```text
sum_n mult_rho(V_n) z^n = ( sum_i z^{e_i} ) / ( (1-z^12)(1-z^20) )
```

### 1.3 The anchors

```text
mu_Lambda   = rho_Lambda^(1/4),        rho_Lambda = Lambda c^4 / (8 pi G)
sqrt(Omega) = R / l_P,                  R = sqrt(3 / Lambda)
Lambda      = 3 Omega_L H0^2 / c^2      (H0, Omega_L: Planck 2018)
```

`H0`, `Omega_L`, `G`, `c`, and `hbar` are declared measurement inputs (Planck 2018 + CODATA);
`Lambda`, `mu_Lambda`, and `sqrt(Omega_Lambda)` are recomputed from them, not quoted.

### 1.4 The torsion (the corrected definition; this is the finding, § 3)

```text
log T^2 = zeta'_coexact(0) - 2 * zeta'_scalar(0)      (Ray-Singer combination, S^3/2I)
```

For a half-integer irrep the scalar tower is supported at half-integer j (odd n; `V_1|_2I = R1`
is the first occupant). Every irrep has an exact closed form in `phi = (1+sqrt(5))/2`:

```text
T^2(R0)=1        T^2(R1)=phi^-4/4      T^2(R2)=phi^4/4     T^2(R3)=(4/5)phi^-2
T^2(R4)=(4/5)phi^2   T^2(R5)=25/9      T^2(R6)=1           T^2(R7)=9/4        T^2(R8)=4
```

with the Galois identities `T^2(R3)/T^2(R4) = phi^-4`, `T^2(R1)/T^2(R2) = phi^-8`, and the sector
products `T^2(R3) T^2(R7) T^2(R5) T^2(R4) = 4`, `T^2(R1) T^2(R2) T^2(R6) T^2(R8) = 1/4` (exact
inverses). The 24 vacuum-twisted products follow from
`log T^2(rho x sigma) = sum_tau N_{rho sigma tau} log T^2(tau)`.

### 1.5 The comparison (declared data inputs, not re-derived here)

The `(rho, sigma) -> fermion` assignment, each fermion's structural sector, and the PDG comparison
masses are declared data inputs (page's own ranked-table reading; § 4 lists what this script does
NOT recompute). `m_e` is the calibration benchmark, excluded from the scorecard count. A
pre-registered null test (`mass-null-v1.1`, corrected table) reports `p_A = 0.690`: the ×3
proximity count is reproduced by random torsion reassignment and is therefore not evidence for the
specific torsion values; that context is restated next to every residual below (task DoD #2).

## 2. Equation-to-code map

Permalinks are commit-pinned to `3786428c` (the script is new in this PR, so `blob/main` would not
resolve until merge). Every row was mechanically re-verified against this exact commit, after two
rows were found wrong during the adversarial audit (§ 5, item 3).

| Equation / quantity (§ 1) | Code |
| --- | --- |
| 2I character table (McKay recursion `V_{n+1} = V1 x Vn - V_{n-1}`), McKay graph, distances | [`m8_3_mass_reproducer.py:62-147`](https://github.com/openwave-labs/openwave/blob/3786428c/openwave/xperiments/m8_mit/research/scripts/m8_3_mass_reproducer.py#L62-L147) `build_chartable`, `mult_in_Vn`, `tensor_mult`, `bfs_distances` |
| Kostant exponents (§ 1.2 generating function) | [`m8_3_mass_reproducer.py:150-167`](https://github.com/openwave-labs/openwave/blob/3786428c/openwave/xperiments/m8_mit/research/scripts/m8_3_mass_reproducer.py#L150-L167) `kostant_exponents` |
| C_geom (§ 1.2) | [`m8_3_mass_reproducer.py:171-177`](https://github.com/openwave-labs/openwave/blob/3786428c/openwave/xperiments/m8_mit/research/scripts/m8_3_mass_reproducer.py#L171-L177) `c_geom` |
| Anchors: H0/Omega_L/G/c/hbar declared, Lambda/mu_Lambda/sqrt(Omega) recomputed (§ 1.3) | [`m8_3_mass_reproducer.py:180-200`](https://github.com/openwave-labs/openwave/blob/3786428c/openwave/xperiments/m8_mit/research/scripts/m8_3_mass_reproducer.py#L180-L200) `anchors` |
| Twisted scalar + coexact spectra, quasi-quadratic fit, Hurwitz-zeta reduction, `zeta'(0)` (§ 1.4) | [`m8_3_mass_reproducer.py:203-292`](https://github.com/openwave-labs/openwave/blob/3786428c/openwave/xperiments/m8_mit/research/scripts/m8_3_mass_reproducer.py#L203-L292) `quasi_quadratic_coeffs`, `HurwitzCache`, `zeta_from_coeffs`, `zeta_prime_at_0`, `torsion_T2` |
| Sign convention fixed on R7 | [`m8_3_mass_reproducer.py:393-399`](https://github.com/openwave-labs/openwave/blob/3786428c/openwave/xperiments/m8_mit/research/scripts/m8_3_mass_reproducer.py#L393-L399) `sign_fix` |
| 8 gated closed forms (R1-R8); R0=1 asserted; both Galois ratios; both sector products (§ 1.4) | [`m8_3_mass_reproducer.py:507-534`](https://github.com/openwave-labs/openwave/blob/3786428c/openwave/xperiments/m8_mit/research/scripts/m8_3_mass_reproducer.py#L507-L534) Part E gate block in `main` |
| 24-product propagation via `N_{rho sigma tau}` (§ 1.4) | [`m8_3_mass_reproducer.py:402-419`](https://github.com/openwave-labs/openwave/blob/3786428c/openwave/xperiments/m8_mit/research/scripts/m8_3_mass_reproducer.py#L402-L419) `products_24` |
| Declared assignment, sectors, PDG masses, null-test provenance, rank order (§ 1.5) | [`m8_3_mass_reproducer.py:313-370`](https://github.com/openwave-labs/openwave/blob/3786428c/openwave/xperiments/m8_mit/research/scripts/m8_3_mass_reproducer.py#L313-L370) `PDG_GEV`, `SM_ASSIGN`, `DECLARED_RATIO`, `NULL_V11`/`NULL_V10`, `CANONICAL_RANK_ORDER`, `NOT_RECOMPUTED` |
| Scorecard counting (compatible / adjudicated / out-of-sector), ranking, JSON output | [`m8_3_mass_reproducer.py:555-638`](https://github.com/openwave-labs/openwave/blob/3786428c/openwave/xperiments/m8_mit/research/scripts/m8_3_mass_reproducer.py#L555-L638) Part G block in `main` |
| Mutation-test registry + coverage/duplicate-id meta-guard | [`m8_3_mass_reproducer.py:644-754`](https://github.com/openwave-labs/openwave/blob/3786428c/openwave/xperiments/m8_mit/research/scripts/m8_3_mass_reproducer.py#L644-L754) `if args.mutation_tests:` block in `main` |

Data: [`m8_3_masses.json`](https://github.com/openwave-labs/openwave/blob/3786428c/openwave/xperiments/m8_mit/research/data/m8_3_masses.json)
(24 ranked entries, the scorecard, both null-test runs, and the `NOT_RECOMPUTED` list, all in one
artifact).

## 3. Results

### 3.1 The finding: a real defect in the reproduction target

Building Parts A-F against the (then-current) published mass-spectrum page, this script's own
from-scratch spectral-zeta computation of the four half-integer torsion singles (R1, R2, R6, R8)
*disagreed* with the published page values. The disagreement was exact and diagnosable: the page
values equalled `exp(+zeta'_coexact(0))`, the coexact one-form contribution alone, to 4/4 figures.
Root cause: the scalar tower of a half-integer bundle is supported at half-integer j, and the
original page computation, finding zero scalar multiplicity at *integer* j, had recorded "scalar
sector identically zero" and omitted the `-2*zeta'_scalar(0)` term. Restoring it collapses every
half-integer torsion to an exact golden-ratio closed form (§ 1.4), on equal footing with the four
integer-spin irreps that already had closed forms; the earlier reading that "the phi mechanism is
structurally absent for half-integer irreps" was an artifact of the truncation.

This was reported upstream. mode-identity-theory corrected the page, landed a reproducible
verification artifact (`torsion-correction.test.py`: 12 gates, a coverage-enforced mutation
registry, and a diagnosis gate that reproduces the pre-correction values as exactly the coexact-only
quantity), re-ran its pre-registered null test as `mass-null-v1.1`, and re-deposited the archived
manuscript. This script's own validation targets are updated to the corrected values below.

### 3.2 Gate results (23/23, standard mode)

| Part | Gates | Result |
| --- | --- | --- |
| A: character table, McKay graph, distances | 3 | 3/3 |
| B: Kostant exponents (9/9 irreps) | 1 | 1/1 |
| C: C_geom (8/8 irreps) | 1 | 1/1 |
| D: anchors (mu_Lambda, sqrt(Omega)) | 2 | 2/2 |
| E: torsion (8 gated closed forms R1-R8, 2 Galois ratios, 1 sector-products check) | 11 | 11/11 |
| F: 24-product table (24/24 entries) | 1 | 1/1 |
| G: scorecard, declared ratios, null-test provenance, rank order | 4 | 4/4 |
| **Total** | **23** | **23/23** |

(`T^2(R0) = 1` is a declared convention, asserted rather than independently gated; it is not one
of E's 11 gates. 8 gated closed forms + 2 Galois ratios + 1 sector-products check = 11.)

`--precise` mode (dps 60, jmax 90) reproduces the same 23/23. Full numeric output:
`m8_3_masses.json`.

### 3.3 The corrected scorecard

With `m_e` as the calibration benchmark (loop-closure cross-check, not counted), of the remaining 8
charged fermions:

| Fermion | Address | Ratio | Status | Null-test context |
| --- | --- | --- | --- | --- |
| μ | (R8, std) | 1.03 | assigned, within ×3 | |
| s | (R8, std) | 1.10 | assigned, within ×3 | |
| τ | (R4, gal) | 2.74 | assigned, within ×3 (weakest counted) | |
| t | (R2, triv) | 0.93 (raw pred/obs; 7%) | assigned, within ×3 | |
| d | (R8, gal) | 3.22 | assigned, **outside** ×3 | |
| b | (R4, gal), colored channel | 1.17 | compatible, but out-of-sector (its own sector is R2) | |
| u | none | — | unassigned (former ~6% hit was the coexact-only artifact) | |
| c | none | — | unassigned | |

**5 of 8 compatible within ×3; 4 of 8 survive sector-first adjudication.** Restated per DoD #2:
this ×3 hit rate is REPORTED ONLY. `mass-null-v1.1` (corrected table) gives `p_A = 0.690`: random
reassignment of the torsion factors across the fixed quantum-number slots reproduces or exceeds the
observed compatible-coverage score in 69% of draws (null S1 mean 5.02 vs observed 5), so the ×3
count is not evidence for the specific torsion values. The pre-correction run, `mass-null-v1.0`
(`p_A = 0.174`), is superseded and kept only as history.

### 3.4 Rank ordering (cross-check against the canonical page)

All 24 formula-computed masses, sorted, reproduce the canonical page's rank order exactly (gate
`G4`), including the five featured addresses at ranks 10 (e), 13 (d), 15 (μ/s), 17 (τ/b), and 22
(t). Full 24-row table: `m8_3_masses.json`.

## 4. Not computed

| Item | Why |
| --- | --- |
| The `(rho, sigma) -> fermion` assignment itself, and each fermion's structural sector (e.g. the bottom quark's sector is R2, making its numerically compatible (R4, gal) address out-of-sector) | Declared data input (the page's own ranked-table reading). The representation-theoretic gate machinery that derives quantum numbers and sectors (the Coxeter-Galois gate, the Z3/Z4/Z5 stabilizer decompositions) is a separate concern from reproducing the mass formula's arithmetic, and is not re-derived here |
| `T^2(R0) = 1`, the non-acyclic ledger convention | ASSERTED per the page's own topological argument (the canonical integral-cohomology value on the two diagonal products) |
| The overall torsion sign convention | Fixed once on the R7 closed form (declared); every other closed form (the remaining 8 irreps, both Galois ratios, both sector products) is then a genuine, independent check, not circular |
| `H0`, `Omega_L`, `G`, `c`, `hbar` | Declared measurement inputs (Planck 2018 + CODATA). `Lambda`, `mu_Lambda`, and `sqrt(Omega_Lambda)` are NOT in this row: they are recomputed from these five (§ 1.3, Part D), not quoted, and not derived from the postulate |
| PDG comparison masses | Measurement inputs (the page's own declared figures) |
| The dead-zone entries' physical status (sterile neutrino / warm dark matter candidates vs. structural residuals) | Open question on the source page itself; outside this script's scope, which is the formula's arithmetic and its comparison, not particle-physics interpretation |
| Independent-method reproduction of the corrected torsion closed forms | The corrected closed forms are established by a single implementation (this script and the mode-identity-theory artifact use closely related but not disjoint methods); a genuinely independent-method reproduction is tracked as [M8.8](../m8_roadmap.md#backlog), placed outside both M8.5 sub-deliverables by the [M8.5-A protocol § 0](m8_5a_reproduction_protocol.md) |

## 5. Adversarial audit record

Independent agent, given the script and an earlier draft of this note but not the derivation
history, tasked to refute every claim and re-walk the equation-to-code map, before this PR's
history was cleaned up into its final commits. Every finding below was addressed and independently
re-verified before this note was finalized (this section's permalinks and § 2/§ 3.2 already reflect
the corrected state, at `3786428c`, the commit containing the fixes).

| Audit item | Verdict |
| --- | --- |
| Recompute: does the script run standalone and reproduce 23/23 in both standard and `--precise` modes? | CONFIRMED |
| Mutation coverage: does every gate id have a matching mutation, and does every mutation actually go red? | CONFIRMED for the literal run, but ISSUE FOUND in the mechanism (item 1 below) |
| Equation-to-code map: does every § 1 equation resolve to the claimed function at the claimed lines? | REFUTED: 2 of 10 rows cited wrong line ranges in the audited draft (item 3 below) |
| Scorecard arithmetic: recompute the 5-compatible/4-adjudicated count independently from `m8_3_masses.json`'s raw entries, without trusting the script's own count | CONFIRMED, independent hand-count matches exactly (5: s, mu, tau, b, t; 4: s, mu, tau, t) |
| Scope honesty: does § 4's not-computed list match what the code actually declares vs. computes (no silent assertion dressed as a derivation)? | ISSUE FOUND: two declared inputs were not itemized (item 2 below); no hidden fudge factor found |
| Null-test framing: is `p_A = 0.690` presented as removing evidential weight from the ×3 count, not as validating the torsion map? | CONFIRMED |

**Findings and disposition:**

| # | Finding | Disposition |
| --- | --- | --- |
| 1 | The coverage meta-guard (`gate_ids == target_ids`, both sets) cannot see two distinct `gate()` call sites accidentally sharing one id: the duplicate collapses into one set entry and reads as "covered" by whichever mutation targets that id, even though one call site's own logic was never attacked. Empirically demonstrated by the auditor via injection (a vacuous always-true gate reusing an existing id: `GATES` grew to 24, `gate_ids` stayed at 23, coverage still reported `True`). No such collision exists in the audited script. | Accepted. Added `no_dupe_ids = len(GATES) == len(gate_ids)` to the meta-guard (now required alongside coverage and all-red). Re-ran the auditor's own injection against the fixed script: `NO DUPLICATE GATE IDS: False (24 gate() calls, only 23 distinct ids)`, `META-GUARD FAILED`, exit 2. |
| 2 | Two declared inputs were not itemized in `NOT_RECOMPUTED` / § 4, though they are the same category as the items that were: the `mass-null-v1.1`/`v1.0` `p_A` and S1 figures (quoted from the external null-test analysis, not computed by this script), and each fermion's `SM_ASSIGN` role bucket (a declared label alongside the address and sector, not derived from the ratio at runtime). No hidden fudge factor was found; the closed-form torsion values were confirmed to be genuinely computed (~26s of real spectral-zeta computation, not copied through). | Accepted. Both items added to `NOT_RECOMPUTED` and § 4. |
| 3 | Two of ten equation-to-code rows cited wrong line ranges in the audited draft: the "sign convention / closed forms / Galois / sector products" row pointed at Part B/C validation code (Kostant exponents), not the sign-fix function or the Part E gate block; the "mutation registry + meta-guard" row's cited range excluded the meta-guard itself (the exact code the row is about) and the G3/G4 mutations. | Accepted. § 2 rebuilt against commit `3786428c` with every boundary re-verified by `grep`/`sed` against the actual file (not eyeballed), including splitting the sign-fix and Part-E-gate-block mappings into two rows since they are non-adjacent in the file. |
| 4 | § 3.2's Part E row said "9 closed forms," which doesn't arithmetically sum to its own stated row total of 11 (9+2+... would be 12); the real count is 8 gated closed forms (R1-R8) plus 2 Galois ratios plus 1 sector-products check = 11. `T^2(R0) = 1` is asserted by convention, not gated. Same off-by-one appeared in the script's own docstring and console header. | Accepted. Fixed in § 3.2 (with the arithmetic spelled out), the script's docstring, and the Part E console header (commit `3786428c`). |

No finding required a change to the underlying mathematics, the corrected torsion values, the
scorecard result, or the null-test interpretation; all four are about the audit surface's own
legibility and completeness, which is exactly what this method note standard exists to catch.

---

*Correction history and scope statements above are the author's own record, consistent with the
platform's AI-hygiene standard: the defect was found in the course of ordinary reproduction, not
concealed, and is reported exactly as diagnosed.*
