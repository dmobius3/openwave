# M8 / MIT, THEORY CANONICAL (specs of record)

> **Purpose.** The single place that pins the MIT spec the M8 column is graded
> against: the arena, the background mode, the operator layer, the particle map, the
> mass formula, and the known tensions, each version-pinned to a published record
> ([`../theory/_CITATIONS.md`](../theory/_CITATIONS.md)). When other docs disagree,
> this file is canonical; when THIS file changes, the change is dated and the record
> id updated. The author owns the science; edits to the spec itself come from the
> author (PR or discussion), edits to verification status come from platform runs.
>
> **Scaffold-stage caveat (2026-07-21).** Everything below is TRANSCRIBED from the
> author's records and onboarding submission
> ([discussion #312](https://github.com/openwave-labs/openwave/discussions/312)),
> not yet verified in-platform. Verification status is tracked per item; the first
> flip is the M8.1 certification gate ([`m8_roadmap.md`](m8_roadmap.md)).

## 1. The arena (topology is the input)

| Object | Definition | Record |
| --- | --- | --- |
| S³ | the three-sphere of curvature radius R (the master length scale; see the R-problem, § 6) | framework deposit [10.5281/zenodo.18064856](https://doi.org/10.5281/zenodo.18064856) |
| S³/2I | the quotient by the binary icosahedral group 2I (\|2I\| = 120; the Poincaré homology sphere) | framework deposit; coexact-gap paper |
| The Möbius edge | a conic Möbius band whose boundary circle `S¹ = ∂(Möbius)` carries the edge layer; the cone point admits a family of self-adjoint extensions | first-eigenvalue paper (SSRN 6968741, author registry) |
| Boundary condition | anti-periodic (twisted): matter modes return only under the double cover. The structural home of 720-degree return | first-eigenvalue paper |

## 2. The background mode + the clock

| Object | Definition | Status |
| --- | --- | --- |
| The standing wave | `Ψ = cos(t/2)` sampled on the fixed geometry; not a dynamical field, the sampling structure particles are read from | assumed (the postulate) |
| The Waltz clock | `dt/dτ = S^(−1/2)` with `S = sin(t/2)`, mapping phase to time | assumed; the exponent −1/2 is empirically pinned (integer alternatives excluded at Δχ² > 60) but NOT derived from the embedding (OQ3) |
| Dynamics | none native: no Lagrangian, no equation of motion; Einstein's field equations kept unchanged as the geometry's dynamics | the defining gap; supplying it is the M8 program ([`m8_background.md`](m8_background.md)) |

## 3. The operator layer (the decidable mathematics)

| Result | Statement | Record | Platform status |
| --- | --- | --- | --- |
| First positive eigenvalue | the twisted Laplacian on the conic Möbius band has first positive eigenvalue `λ₁ = 2/R²`, stable across the cone's self-adjoint extensions | SSRN 6968741 (author registry); operator statement: [first-eigenvalue.md](https://github.com/dmobius3/mode-identity-theory/blob/main/files/framework/files/bedrock/files/first-eigenvalue.md) (author repo, shared on #312, 2026-07-21) | ✅ VERIFIED IN-PLATFORM (M8.1, 2026-07-21): blind two-agent verification, adversarially audited; λ₁⁺ = 2/R² (narrow), α₀(α₀+1)/R² (wide), the 2R/e stability threshold and the −4e^(−2γ)/δ₀² defect state confirmed at 10-digit precision ([`findings/m8_1_method_note.md`](findings/m8_1_method_note.md)) |
| Λ from the edge mode | `Λ = 3/R²`, the eigenvalue `2/R²` carried by the Gauss-Codazzi factor 3/2 (3 = Ricci trace derived, 1/2 imported from GR) | Λ ground-mode paper [10.5281/zenodo.18307314](https://doi.org/10.5281/zenodo.18307314) | 🚧 the 2/R² input is now verified (M8.1 ✅); the Gauss-Codazzi 3/2 step and the R-problem (§ 6) remain untested |
| Coexact spectral gap | Yang-Mills-sector mass gap `4/R²` on S³/2I from McKay distance for flat bundles | SSRN 6968698 (author registry); YM paper [10.5281/zenodo.18463584](https://doi.org/10.5281/zenodo.18463584) | 🚧 eigensolve planned after M8.1 |
| Galois pair | the affine rho-index conversion and the E8-filling Galois pair on the Poincaré homology sphere | SSRN 7129118 (author registry) | 🚧 has an author-side test script (`galois-pair.test.py`); not yet reproduced here |

## 4. The particle map (representation theory, not field configurations)

| Object | Definition | Status |
| --- | --- | --- |
| Particle | an irreducible representation of 2I at a McKay-lattice position (a mode slot, NOT a soliton or defect) | structural |
| Generations | the same irrep read in the three flat connections of S³/2I (trivial, standard, Galois): exactly three, zero freedom | structural, the ledger's strongest class |
| Color | the Z₃ face stabilizers of 2I: triplet vs singlet per irrep | structural |
| Torsion-to-slot map | Reidemeister-torsion assignment of irreps to mass slots | ❌ author's pre-registered null (`mass-null-v1.0`: p = 0.174); the map is NOT evidence-bearing; slots are graded structurally only |
| Unassigned | charm unplaced; 8 of 24 slots empty; rank-16 / dead-zone states named but uncomputed | open, listed honestly |

## 5. The mass formula (analytic sector)

```text
m = μ_Λ · C_geom · (√Ω)^(dist/30) · T²

μ_Λ  = ρ_Λ^(1/4) ≈ 2.25 meV     (the mass-sector anchor, from measured Λ: an INPUT)
C_geom = Kostant geometric weight of the irrep
Ω     = (R/ℓ_P)² ≈ 10¹²²          (the calibration-web hub)
dist  = McKay distance of the slot
T     = Reidemeister torsion of the assignment (evidence-null, § 4)
```

Calibration anchors, one per sector: H₀ (edge), Λ or α (surface), m_e (mass
normalization). Record: mass-spectrum paper
[10.5281/zenodo.18603975](https://doi.org/10.5281/zenodo.18603975). Platform status:
🚧 the M8.3 reproducer (every constant recomputed from its definition) is the entry
bar; graded at the ledger weight regardless of hit rate (§ 6).

## 6. Known tensions + evidence grading (adopted from the author's claim ledger)

| Item | The honest statement |
| --- | --- |
| Cycle 2 (α) | α is input AND output: it calibrates R, then the scaling law "predicts" α. The 0.5% match is a consistency check; the non-circular content is Λ at ~24% downstream |
| Cycle 7 (the R-problem) | the two independent routes to R disagree ~4× (coupling route ≈ 5.3 Gpc, mass route ≈ 20 Gpc; ~14× in Λ). Every R-dependent number inherits this until one route is shown correct (OQ2) |
| m_e ↔ Λ | the closure loop holds only to ~11% |
| Residuals | down quark 3.2×, top 3.9×, tau a soft miss; not smoothed |
| Documented negatives | torsion null (p = 0.174), SPARC a₀(z) coherence-scale mechanism (falsified, pre-registered pipeline [10.5281/zenodo.20271702](https://doi.org/10.5281/zenodo.20271702)), H₀ bimodality (unimodal, p = 0.217) |
| Grading rule | platform tasks target the STRUCTURAL ladder (slot structure, gap ratios, generation count), never the 24-entry numeric table ([`m8_background.md § 3`](m8_background.md)) |

## 7. Consumption rules (standing)

| Rule | Meaning |
| --- | --- |
| Version-pin | every result cites the record id it was computed against; the author's repo is context, the deposits are the citable layer |
| Pre-register | gates, conventions, and success criteria BEFORE numerics; forks reported with all numbers, never tuned toward published values |
| No calibrated conventions | derive and pre-register, or record as a fit with its search space |
| Audit | substantive claims get an independent adversarial pass before they are trusted ([`AI_HYGIENE.md`](../../../../AI_HYGIENE.md) § 1) |
| Author-gated | intent, provenance, and definition questions go to the author; the platform does not resolve them by inference |

## OPEN QUESTIONS

| ID | Question | Route | Status |
| --- | --- | --- | --- |
| OQ1 | Does ANY reasonable nonlinear field equation on S³/2I have defect or standing-wave solutions whose energies realize the McKay slot structure? | M8.4 (the decisive science) | 🚧 open |
| OQ2 | The R-problem: which of the two routes (coupling ≈ 5.3 Gpc vs mass ≈ 20 Gpc) is the correct determination of the master scale, and what breaks in the other? | author-side; platform can bound it once M8.1/M8.3 land | 🚧 open |
| OQ3 | Can the Waltz exponent −1/2 be DERIVED from the embedding (or from a field dynamics) rather than empirically pinned? | M8.4 lineage | 🚧 open |
| OQ4 | Does the McKay-distance rule map onto M5's measured lepton hierarchy (1 : 5.9 : 15.1 at eigenvalue level, `E ∝ Λ³` fixed)? | M8.6 (bounded, no simulation) | 🚧 open |
| OQ5 | Which topological-defect sectors exist at all on S³/2I (homotopy of the target space restricted to the quotient), and does the anti-periodic double-cover condition create sectors flat space lacks? | M8.4 prerequisite analysis | 🚧 open |
| OQ6 | Does the double-cover return supply a genuine spin-½ mechanism once a dynamical field lives on the arena (the 720-degree row most columns leave 🚧)? | downstream of OQ1 | 🚧 open |
