# M8.10: FIRST CORRECTION AT THE LEVEL-6 CRITICAL RAYS, the range equation one order past the fourth bedrock paper

> Roadmap row: [`../m8_roadmap.md`](../m8_roadmap.md) M8.10 (**author-proposed, maintainer-run**).
> Standing: [#512](https://github.com/openwave-labs/openwave/discussions/512#discussioncomment-18415036).
> Parent verification and run-format template: [`m8_1_2_task_details.md`](m8_1_2_task_details.md).
> Status: 📋 PROPOSED 2026-09-12, author-drafted; registration and go pending.

## TASK PLANNING (author-proposed; registration and go pending)

### The question

The fourth bedrock paper, verified blind under M8.1.2, stops at the leading reduced problem. Its critical rays "are then the critical directions of the leading reduced problem, and not solution branches: an actual branch requires the range equation described in Section 6.2, which is not carried out." This task carries the range equation one order, at the four critical rays of the paper's § 5.8 (coherent `v_3`, zonal `v_0`, octahedron, hexagon) in the sectors `3′` and `4`: eight formal branch expansions. It computes the first correction `ξ` to the state and the next coefficient `λ₄` of the eigenvalue, and checks that each ray persists at that order. It establishes no branch.

### Standing: the #512 ruling

| Point | Ruling |
| --- | --- |
| Run-before-write | M8.1.2 meets the first condition; nothing else needs to land first |
| Identity | a new identity, M8.10, filed as a BACKLOG row with its task doc in the PR that carries the pre-registration; not an M8.1.x sub-ID, since the range equation is work the paper names and does not do |
| Budget | one pre-registration within the § 12.2 budget of about 8,000 words; nothing further lands on it until it has run |
| Defect claims | a claim that this frozen text is defective is reproduced by the maintainer with independent code before it is ratified |
| Solver side | the pre-registration states which side of the closed spectral route its solver is on |
| The law | derived and pinned here, starting from nothing in `m8_5c/design_inputs/`, whose `jacobian_check.py` carries the global law ([#508](https://github.com/openwave-labs/openwave/issues/508)) |

### Which side of the closed spectral route

The analytic side, to a stated order: the state to order `a³` and the eigenvalue to order `a⁴`. Everything is finite algebra in the sectors' harmonics: Clebsch-Gordan coefficients and averages over the 120 elements of the group. There is no truncation, because the cubic of a level-6 state has no content above level 18 and every sum is finite; no continuation, no time-stepping and no quadrature; and no code or data from the C rooms. What separates this from the chassis gate 7 was built for is the equation being solved: never the nonlinear equation at finite amplitude, only the linear equations of successive orders.

### What this task is not

Not branches at finite amplitude, existence or stability. Not slot selection: the sectors are chosen, not selected. Nothing about OQ1, the stress tensor or the metric. It is not an in-platform field-dynamics validation for [M8.7](m8_7_task_details.md) and neither satisfies nor weakens that gate; the M8.7 gate-line correction filed with this task is bookkeeping and leaves the condition unchanged.

### Ownership and run format

The author supplies candidate values from a fresh derivation and an offered instrument, and then stops; the designer may edit the claims and freezes them at go. The recommended run follows M8.1.2's format: two blind agents with separate implementations, an adversarial audit that commits its own method first, then the designer's comparison against the frozen claims. The author has seen the September estimates and all eight fresh values, and the author's package shares one convention stack, so only separate blind implementations can tell a convention reproducing itself from independent recovery. The format is the designer's to set, and the label follows it: two blind agents and an audit earn blind, and a single maintainer reproduction earns [independent-method reproduction](../m8_roadmap.md#conventions).

### Sources of record

| Item | Pin |
| --- | --- |
| Paper | *The Surviving Ray*, Zenodo version DOI [10.5281/zenodo.22681502](https://doi.org/10.5281/zenodo.22681502), PDF MD5 `d2405316e20c060f53d338c1516298bc`, as pinned in M8.1.2. Statements used: § 2.5 (the equation), Lemma 4.1 (the level-6 density carries only ranks 0 and 6, restated there as `R_K(P) = 0` for `K = 1..5`), Proposition 3.3 (the spin-8 channel vanishes exactly on the time-reversal-invariant rays), Lemma 5.3 and Lemma 5.3(b) (criticality), Corollary 5.6 and the remark after it (section normalization and `Q_d`), § 5.8 (the four rays), § 6.2 (the range equation) |
| Parent verification | M8.1.2 ([#540](https://github.com/openwave-labs/openwave/pull/540)), whose 21 claims include C3, D1 and D7, all reproduced blind |
| Engine | `openwave/xperiments/m4_ewt/wave_engine.py`, blob SHA-256 `6520ca41762cc4078660c00fd7bc4279094927f9182e09b0b373aa5dbe72a001`, last changed in `b910b5e40e145bc53ddb4c4d7a122f2f30876f09` (2026-07-27) and read at `main` `5dd2cc048ed91bb1d8bff8b675903ae03b6f2803`: the potential table (lines 142-176, `v_mode 1` at line 150), the restoring force (line 339) and the continuous PDE (line 507) |
| Author's package | the fresh derivation, pinned file by file in [the firewall](#the-firewall) and quarantined until the verdict |

## THE LAW, DERIVED AND PINNED

| Step | Statement |
| --- | --- |
| (a) Engine | `v_mode 1` (`cubic_nls`): `V = (c1/4)u²` and `dV = c1·u·ψ` with `u = ψ·ψ`, on one real three-component field, evolved by `∂²ψ/∂t² = c²∇²ψ − dV(ψ)`. The engine leaves the sign of `c1` to the caller ("For a focusing (attractive) potential, use c1 < 0") |
| (b) Modeling extension, called that | the engine's real inner product `ψ·ψ` becomes the Hermitian `ψ†ψ` on sections of the flat unitary bundle `E_σ` over `X = S³/2I`, so `V = (g/4)(ψ†ψ)²`, with the Laplacian of `X` in place of the engine's `∇²`. The complexification lives here, not in step (c). Gradients are taken with respect to the underlying real Hilbert metric, for which `V` has gradient `g(ψ†ψ)ψ`; the Wirtinger derivative `∂V/∂ψ̄` is half of that. The engine's bytes supply the local form, not `3′` or `4` sections |
| (c) Reduction | with `ψ = e^{iωt}φ`, the PDE of (a) under (b) gives `(−c²Δ − ω²)φ + g(φ†φ)φ = 0`, which is the paper's § 2.5 equation with `c = 1`, `λ = ω²` and `g = c1` |

Every value is stated per power of `g`, which covers both signs.

## SETTING AND DERIVATION (equations first)

**Setting.** `(−Δ − λ)ψ + g|ψ|²ψ = 0` on sections of `E_σ`, `σ ∈ {3′, 4}`, with `R = 1` and the paper's normalized Haar measure (`∫_X 1 = 1`), the conventions under which the paper's `Q_d` and M8.1.2's values hold. Level `n` has `−Δ` eigenvalue `n(n+2)`. The kernel `K` is the seven-dimensional level-6 block. A fibre ray `v ∈ V_3` lifts to a block section `Φ = Φ_{σ,v}`, normalized by `∫_X |Φ|² = 1`; the paper's remark after Corollary 5.6 warns that normalizing the fibre instead gives a different coefficient, so every value is stated for `Φ`.

**Notation, fixed against the paper's.** `N(ψ) = |ψ|²ψ` is the pointwise cubic (the paper's `𝓝` is its block projection). `Π_K`, `Π_⊥ = I − Π_K` and `Π_n` project onto the block, its complement and level `n`. `A = −Δ − 48`, inverted on the range of `Π_⊥`. The first correction is `ξ`. The paper's `χ` (a stabilizer character), its `P = ηη†` and its `Π_{u⊥}` (a projector inside the fibre) keep their own meanings and are not used here.

**Expansion.** `ψ = aΦ + a³ξ + O(a⁵)` and `λ = 48 + λ₂a² + λ₄a⁴ + O(a⁶)`, with `a > 0`, the phase fixed so that `⟨Φ, ψ⟩` is real, and `Π_K ξ = 0`. At order `a³`, `Aξ − λ₂Φ + gN(Φ) = 0`. Its block part gives `λ₂ = g·Q_d([v])`, since at a critical ray `Π_K N(Φ) = Q_d Φ` (Corollary 5.6); its range part gives `ξ = −g A⁻¹ f` with `f = Π_⊥ N(Φ)`. At order `a⁵`, pairing with `Φ` gives `λ₄ = g⟨Φ, DN_Φ[ξ]⟩`. With `DN_Φ[h] = 2 Re(Φ†h) Φ + |Φ|² h`, the pairing is `2X + X̄` with `X = ⟨N(Φ), ξ⟩ = −g⟨f, A⁻¹f⟩`, which is real, so

`λ₄ = −3g² ⟨f, A⁻¹f⟩ = −3g² Σ_n ‖Π_n N(Φ)‖² / (n(n+2) − 48)`,

summed over the sector's levels above 6. Both `‖Π_n ξ‖² = g² ‖Π_n N(Φ)‖² / (n(n+2) − 48)²` and `λ₄` are therefore fixed by the level norms of one cubic.

**Strict sign.** Above level 6 every denominator is positive and every summand nonnegative, so `λ₄ ≤ 0`. At least one summand is nonzero. If `f = 0`, then `N(Φ) = Q_d Φ` pointwise, so `|Φ|²` takes only the values `0` and `Q_d`; being continuous on the connected `X`, it is constant, and the normalization makes it 1, forcing `Q_d = 1`. But `Q_d = 1 + w₆(d) r̂₆` with `w₆(d) > 0`, and `r̂₆ = ‖ρ₆‖²/‖v‖⁴`, the normalized top multipole, is nonzero at all four rays. So `λ₄ < 0` whenever `g ≠ 0`.

**Persistence at the four rays.** By Lemma 5.3, the stabilizer `H` of each ray acts on `Φ` by a character whose isotypic space in the block is one-dimensional. `N`, `Π_K`, `Π_⊥` and `A⁻¹` all commute with `H`, so `f` and `ξ` carry the same character as `Φ`, and `Π_K DN_Φ[ξ]` lies in `ℂΦ`. The order-`a⁵` block equation therefore has no component orthogonal to `Φ`, and no kernel correction `a³v` with `v ⊥ Φ` is needed. The two further critical rays of Lemma 5.3(b), the pentagonal pyramid and the trigonal prism on the `D₃` locus, each lie on a fixed locus of dimension one. There the argument fails and the direction can move at the next order, which is why they are out of scope; claim D2 checks it at the pentagonal pyramid.

**Levels.** From the 2I characters, the sectors occupy different levels: `3′` at 6, 10, 14, 16 and 18, and `4` at 6, 8, 12, 14, 16 and 18, each once except level 18 in sector 4, which occurs twice (`Hom_2I(4, V_9)` is two-dimensional). A reproduction must project onto the whole two-copy space there, because the one-copy Schur reduction used at level 6 is unavailable.

**Why some levels vanish.** `Π_n N(Φ)` carries the character of `Φ`'s stabilizer, so its fibre factor lies in that character's isotypic part of `V_{n/2}`, and a level vanishes identically where that part is zero. Three rules follow. At the zonal ray (stabilizer `O(2)`, character `−1` on the flip), level `n` vanishes exactly when `n ≡ 0 (mod 4)`. At the octahedron (the `A₂` character of the octahedral rotation group), level `n` vanishes exactly when `A₂` is absent from `V_{n/2}`, which in range is `n = 8, 10, 16`. At the coherent ray (`U(1)`), no level vanishes. The hexagon's stabilizer, dihedral of order 12, forces nothing in range. Level 16 has a second cause, common to the three time-reversal-invariant rays: with `Θv = λv`, the spin-8 channel `[ρ₆(v) ⊗ v]_8` that carries level 16 (the paper's § 5.7 and Proposition 3.3) becomes `λ[[v ⊗ v]_6 ⊗ v]_8`, a holomorphic cubic, which factors through `Sym³V_3` and vanishes because `Sym³V_3` contains no `V_8` (M8.1.2's A1). At the zonal and octahedron rays the stabilizer forces the same zero independently.

## THE FIREWALL

The requirements follow M8.1.2's; containment and roles are the designer's to fix at go.

| Requirement | Statement |
| --- | --- |
| Repository and paper | solver and audit contexts have no access to the OpenWave M8 tree, the `mode-identity-theory` repository, the deposited paper, or M8.1.2's records |
| Network | no unrestricted web search; any external reference is a generic representation-theory or computer-algebra source approved in advance by the designer |
| Manifest | each agent returns a consulted-material manifest; an empty or unreturned manifest is a protocol failure |
| Audit ordering | the auditor commits its own method and results before it sees the solver's work |
| Group presentation | the group is given by two generating unit quaternions, not by name, and the packet check X0 runs before any claim |

**Where the frozen values live.** In this file, which is an explicit answer key, so the leak path of record is this file and not the handout. At go the designer states where it is held and what prevents a solver session from reading that path.

**Withheld terms.** The author, the model, the repositories, the paper's title, "surviving ray", every frozen value and every claim ID. Permitted: the vocabulary the problem needs, including the four rays as explicit fibre vectors.

**Designer-only material: the author's package.** The fresh derivation's code, `RESULTS.md` (the candidate values with their provenance, regenerated by one command), a post-derivation note recording an observation made after the values were known, which is not claimed here, and the record of the dry run described under instrument qualification. None of it reaches an agent. The designer adjudicates against the frozen claims first and opens the package only afterwards, for provenance comparison. Agreement between the package and the agents is weak evidence, since a shared convention error survives it; disagreement is strong evidence that something is wrong. The package lands in the repository after the verdict, under `research/scripts/m8_10_author/`, byte-identical to these hashes:

| file | SHA-256 |
| --- | --- |
| `m810_core.py` | `9339d87117a4fe524b913f5fe63972f9a19dfaa0d708f059323bde7553f7ecda` |
| `m810_selftest.py` | `d0570144b051b1435c1194773ed10c3158e3b294e38107db541f2e8421ff9251` |
| `m810_parent.py` | `bf4aeb556be3d73ad8adbadf8109c251ec32f8cce93fd90d1f249a6c5a38d965` |
| `m810_levels.py` | `9d00d7ab200247ab806b45b93a80847b011fccb13467addfdaaf1d4623c41de4` |
| `m810_exact.py` | `8436ef9ffc465e6a91878b4b994c1385766167cc8ebc4047d482cf65c72463a9` |
| `m810_gates.py` | `70c6490296272f9be60aec1ce8e45772bd9e93afe746eb13f5a51b57bbf58d1a` |
| `m810_zeros.py` | `531a5f3e87c063e11e6d770d84d259f691e4697a741ebea02201a345acadc21b` |
| `m810_results.py` | `8fe498ff0cf59917bb346ba35a028fcef1557f6cdb6f1edf8d55c4e19d3f4d62` |
| `RESULTS.md` | `9986195a4199d76db970c052bf5ab93b9ad373f3b190bc6dcb31226821a742ae` |
| `POST_DERIVATION_NOTE.md` | `4f6b96a0330ea982bfcfd29e48356897a5a3e8fb114d540b8eb169ade4007620` |
| `audit_dryrun.py` | `e0cbae2dea1750bf0692ea8f6698e7451b983cf0261ff7489b3f2a439ec4d371` |
| `compare_dryrun.py` | `dc45e596664a6fdb3d7bf8c955a97af54ea2c2b37c01bbc10b1e8034acaff4cd` |
| `dryrun/BRIEF.md` | `4918be4276d2a5e3849bfdbe4e997df5ab4c05f081bb00aea1702c8299032541` |
| `dryrun/worklist.md` | `b08d23c6ceba1ba18cc9ba53b3801d3654e11f495c4b9945e97e353692348248` |
| `dryrun/RETURN.md` | `c048a1268a76ea9efcdafa8f1d6614351d8d8b92160312639ad24490b897a2c7` |
| `dryrun/transcript.jsonl` | `a96cb211671baaaf72e8237e52e6244d02ae9c8a53c29b6e734eaa2296d6ddd7` |

**Offered instrument (adopt, modify, or discard).** [`../m8_10/worklist.md`](../m8_10/worklist.md): the conventions and a worklist, carrying the definitions and none of the results. It is the only file a solver receives. Before the room opens, the designer audits it against the frozen claims independently; the author's own review of it is not evidence of non-leakage. The worklist labels the rays `R1` to `R5`: `v_3`, `v_0`, `v_2 + v_−2` and `v_3 + v_−3` (the coherent, zonal, octahedron and hexagon rays, unnormalized), and the pentagonal pyramid. Coverage, for the designer to confirm rather than take from this sentence: items 2 and 3 reach A1 and A2; item 1 reaches B1, whose level-6 row the handout's § 1.4 supplies, so B1 is graded from level 8 up; items 4 and 5 reach B2, B3a and B3b; items 7 and 8 reach Group C; item 6 reaches Group D, with `R5` as D2's negative control; and items 0 and 9 reach X0 and X4.

## DISCLOSURE

**Exploratory runs before this task (2026-08-31 and 2026-09-01).** Local runs followed the hexagon branch in both sectors by Newton continuation in a truncated harmonic basis (levels up to 18, and in sector 4 also up to 30). They ran on the preserved C2 room's code (fifteen modules byte-identical to `m8_5c2/a1/in_room_as_edited/` on `main`, plus an empty package marker) with a script-local pointwise Jacobian: the room's `real_jacobian` was never reached and nothing in `m8_5c/design_inputs/` was imported. With `g = 1` they returned `λ₂ = 1.202360` (sector `4`) and `1.359751` (sector `3′`), equal to the paper's `Q_d` at the hexagon; `‖ξ‖ = 7.142118×10⁻³` and `3.602677×10⁻³`; and `(λ − 48 − λ₂a²)/a⁴ = −5.852×10⁻³` and `−4.186×10⁻³` at `a = 0.1`, the smallest amplitude, with `a` the solution's norm. The fresh derivation reproduces both `‖ξ‖` values to `8×10⁻⁸` relative and both `λ₄` readings to within `2×10⁻⁴` relative, the size of the `O(a²)` offset at `a = 0.1`.

**The derivation.** Author-side, carried out by the author and the author's AI agents, in a folder that imports numpy, scipy, mpmath and its own modules only: nothing from the rooms, the design inputs or the September scripts. It is fresh but not blind, since it was written with the September numbers known. Two author-side review units redlined the outline and then the derivation, three rounds each; that review is not verification. The worklist was also dry-run once before proposal by a fresh author-side context; see [instrument qualification](#instrument-qualification-2026-09-12).

## CANDIDATE PRE-REGISTERED CLAIMS

Frozen by the designer at go, not by the author. Conventions as in the setting above. Tabulated values are `λ₂/g`, `‖Π_n ξ‖²/g²` and `λ₄/g²`, so no sign of `g` is assumed.

**Status of the values: candidate exact rationals.** Each level norm `‖Π_n N(Φ)‖²` was computed by finite algebra in 60- and 100-digit arithmetic and identified as a rational with denominator at most `10¹⁶`. The two precisions give identical fractions; the worst residual at 100 digits is 2.57e-100, and the largest identified denominator is 237,562,624. The level norms need not be rational a priori, since the sector projectors involve the golden ratio, so their rationality is itself a finding. Conditional on those identified fractions, the `λ₄` and `‖ξ‖²` fractions follow by exact rational arithmetic. This is rational identification, not a symbolic derivation; a blind agent that derives the fractions symbolically promotes them to exact results at adjudication.

### Group A: the parent, reproduced rather than new

| ID | Claim | Pass condition | Fail condition |
| --- | --- | --- | --- |
| A1 | C3 in both sectors: the block part of the cubic, written on the fibre, is `(d/7)‖v‖²v − ((7−d)/√91) M₆(v)` | blind derivation returns both coefficients in both sectors | either coefficient differs |
| A2 | `λ₂/g = Q_d` at all eight expansions, as tabulated | all eight fractions match | any differs |

| ray | `λ₂/g`, sector `3′` | `λ₂/g`, sector `4` |
| --- | --- | --- |
| coherent v3 | `1288/1287` | `2289/2288` |
| zonal v0 | `1687/1287` | `168/143` |
| octahedron | `175/143` | `161/143` |
| hexagon | `1750/1287` | `2751/2288` |

### Group B: level structure

| ID | Claim | Pass condition | Fail condition |
| --- | --- | --- | --- |
| B1 | the level sets above 6: `3′` at 10, 14, 16, 18 and `4` at 8, 12, 14, 16, 18, each with multiplicity 1 except level 18 in sector `4`, with multiplicity 2 | derived from the group rather than recognized, with the level-18 multiplicity stated; level 6 is supplied by the handout and not graded | a level missing or extra, or multiplicity 1 at level 18 in sector `4` |
| B2 | the 36 per-level squared norms `‖Π_n ξ‖²/g²`, zeros included | every nonzero entry reproduced exactly as a fraction and every zero as an exact zero | any entry differs, or a zero reported as a small number without an exact argument |
| B3a | the zero census: exactly the tabulated zeros vanish, and no other entry | the census matches, every zero reproduced as an exact zero | a zero missing or extra |
| B3b | every zero has an exact reason: each zero outside level 16 is forced by the ray's stabilizer, and level 16 vanishes at the three time-reversal-invariant rays (at the zonal and octahedron rays also by the stabilizer) | every zero carries an exact argument, by these routes or any other | a zero without an exact argument, including one reported as unresolved |

Sector `3′`, `‖Π_n ξ‖²/g²`:

| ray | `n = 10` | `n = 14` | `n = 16` | `n = 18` |
| --- | --- | --- | --- | --- |
| coherent v3 | `7/195150384` | `18375/1464486053888` | `7/2749593600` | `245/609749017488` |
| zonal v0 | `1225/48787596` | `8575/2860324324` | `0` | `171500/139734149841` |
| octahedron | `0` | `735/220024948` | `0` | `1715/2822912118` |
| hexagon | `8575/780601536` | `735/5720648648` | `0` | `1516795/812998689984` |

Sector `4`, `‖Π_n ξ‖²/g²`:

| ray | `n = 8` | `n = 12` | `n = 14` | `n = 16` | `n = 18` |
| --- | --- | --- | --- | --- | --- |
| coherent v3 | `63/5360582656` | `7/483225600` | `23625/7209777496064` | `49/28242739200` | `427/1927108005888` |
| zonal v0 | `0` | `0` | `11025/14081596672` | `0` | `74725/110407229504` |
| octahedron | `0` | `7/1830400` | `945/1083199744` | `0` | `26901/80296166912` |
| hexagon | `1575/31719424` | `189/644300800` | `945/28163193344` | `0` | `2643557/2569477341184` |

| ray | zero levels, `3′` | zero levels, `4` | reason |
| --- | --- | --- | --- |
| coherent v3 | none | none | none: `U(1)` forces no zero, and the ray is not time-reversal invariant |
| zonal v0 | 16 | 8, 12, 16 | `n ≡ 0 (mod 4)`, forced by `O(2)`; level 16 also by time reversal |
| octahedron | 10, 16 | 8, 16 | `A₂` absent from `V_{n/2}`, forced by the octahedral group; level 16 also by time reversal |
| hexagon | 16 | 16 | level 16 only, by time reversal; the stabilizer forces nothing in range |

### Group C: the new coefficient

| ID | Claim | Pass condition | Fail condition |
| --- | --- | --- | --- |
| C1 | `λ₄/g²` at the eight expansions, as tabulated | all eight fractions reproduced exactly | any differs |
| C2 | `λ₄ < 0` at all eight for `g ≠ 0`, as a derived claim | nonnegative summands shown from the structure, and at least one nonzero summand at each ray, by a general argument or by an exact nonzero entry; reading the sign off the eight totals does not pass | a nonnegative value, or the sign asserted from the totals alone |
| C3 | no common sector multiplier: the four ratios `λ₄(3′)/λ₄(4)` are pairwise distinct, whereas at leading order `(λ₂(3′) − g)/(λ₂(4) − g) = w₆(3)/w₆(4) = 16/9` at every ray | the four ratios reproduced | any ratio differs, or two coincide |

| ray | `λ₄/g²`, sector `3′` | decimal | `λ₄/g²`, sector `4` | decimal |
| --- | --- | --- | --- | --- |
| coherent v3 | `-608931967/36722893315680` | `-1.6581808023e-05` | `-368688201/38687492546560` | `-9.5299068699e-06` |
| zonal v0 | `-1871763250/229518083223` | `-8.1551885748e-03` | `-15820875/15112301776` | `-1.0468871807e-03` |
| octahedron | `-2203040/944518861` | `-2.3324468054e-03` | `-162530109/75561508880` | `-2.1509643125e-03` |
| hexagon | `-1921938515/459036166446` | `-4.1868999776e-03` | `-226441775133/38687492546560` | `-5.8531003233e-03` |

| ray | `λ₄(3′)/λ₄(4)` | decimal |
| --- | --- | --- |
| coherent v3 | `3181358848/1828392507` | `1.739976` |
| zonal v0 | `4889504/627669` | `7.789940` |
| octahedron | `3596800/3316941` | `1.084373` |
| hexagon | `803291852800/1122966354231` | `0.715330` |

### Group D: persistence, with one negative control

| ID | Claim | Pass condition | Fail condition |
| --- | --- | --- | --- |
| D1 | at the four rays, `Π_K DN_Φ[ξ]` lies in `ℂΦ`, and `⟨Φ, DN_Φ[ξ]⟩/g`, computed directly rather than through the level sum, is real and equals C1's `λ₄/g²` (equivalently, `g⟨Φ, DN_Φ[ξ]⟩ = λ₄`) | the component orthogonal to `Φ` vanishes at all eight, by the stabilizer argument or to a stated precision, and the direct pairing returns C1 | a nonzero orthogonal component, or a mismatch |
| D2 (negative control) | at `R5`, the pentagonal pyramid of Lemma 5.3(b) (`v = cos t·v_2 + sin t·v_−3`, `sin²t = 12/25`), the component of `Π_K DN_Φ[ξ]` orthogonal to `Φ` is nonzero in both sectors, so the fixed-ray expansion fails there and no `λ₄` is defined for it; the worklist asks at `R5` only for that component | nonzero in both sectors | zero in either |

### Diagnostics (not claims)

| ID | What it detects |
| --- | --- |
| X0 (packet check, before solving) | the two quaternions are units and generate a perfect group of order 120, as in M8.1.2 |
| X1 (normalization) | normalizing the fibre, `‖v‖ = 1`, instead of the section, `∫_X \|Φ\|² = 1`, scales `λ₂` by `d/7`, the artifact the paper's remark after Corollary 5.6 describes, and every later value moves with it |
| X2 (gradient convention) | taking the Wirtinger derivative as the force halves `g`, which halves `λ₂` and quarters `λ₄` |
| X3 (the level-18 plane) | a one-copy Schur scalar reused at level 18 in sector `4` misses part of that level's norm |
| X4 (range equation) | `Π_⊥(Aξ + gN(Φ)) = 0` checked coefficient by coefficient, with `A` applied through the SU(2) Casimir rather than the eigenvalue table, so that a mislabeled level fails it |

### Feasibility, stated honestly

Most claims are finite evaluations. Two are not spot checks: C2 is a sign argument, not eight signs, and B3b needs each ray's stabilizer and its character, or time reversal at level 16. The dry run reached only the zonal arguments, so a B3b failure in the real run would be a legitimate outcome, not a defect of the packet. Exact reproduction may come from symbolic algebra or from rational identification at a precision and denominator bound the agent states.

### Instrument qualification (2026-09-12)

Before proposal, the worklist was run once by a fresh author-side context: a headless Claude session (Opus 4.6), the author's AI agent, in a room outside every repository holding only the worklist and a short brief. Web tools were disabled and reads of the answer locations denied; a canary run showed nothing answer-bearing loaded at startup, and the transcript audit found no call outside the room in 55 tool calls. The run took 2 h 7 min. The mechanical comparison against the frozen values was 103 comparisons with no discrepancy: the level multiplicities (ten rows, the level-18 multiplicity included), the eight `λ₂/g` (and `R5`'s two against the paper's `r̂₆`), all 36 per-level norms of `ξ` and of the cubic with their zeros, the eight `λ₄/g²` and the four ratios. It also reported the expected qualitative behavior at `R5`, a nonzero orthogonal component and no `λ₄` defined there, without being told that `R5` was a control. This is not a verification: the agent is the author's, and the author wrote the instrument. Four repairs followed, and no expected value changed.

| change | why |
| --- | --- |
| item 2 asks for exact coefficients through the `M_K` | the agent gave C3's coefficients to four digits, exact only in a basis of its own |
| item 9 requires the Casimir to be built and applied | the agent substituted `n(n+2)` and reported the resulting tautology as the check |
| B3 split into the census (B3a) and the arguments (B3b) | the census reproduced exactly, but only the zonal zeros at levels 8 and 12 got exact arguments; level 16 was reduced to one channel and then checked numerically, and the octahedron zeros were checked numerically and put down to an accident. The split separates a census from a derivation that proved to differ in difficulty; B3b's bar is unchanged, since an unresolved report still fails |
| C2's pass condition states what counts for strictness | the agent derived nonnegativity from the structure and strictness from exact nonzero entries, which the old wording left unclear |

The run used the worklist with SHA-256 `b08d23c6ceba1ba18cc9ba53b3801d3654e11f495c4b9945e97e353692348248`; items 2 and 9 were reworded after it. Its brief, return and transcript are in the author's package.

### To be fixed at go (before numerics)

| Item | Note |
| --- | --- |
| Claims frozen | the designer freezes the tables above, edited freely, before any solver launches |
| Answer-key containment | where this file is held, and what prevents a solver reading it |
| Instrument | adopt the offered worklist after confirming it reaches every claim, modify it, or build a fresh one |
| Handout audit | maintainer-side and independent of the author's review |
| Run format | two blind agents and an audit (recommended; earns blind), or a single maintainer reproduction (earns independent-method reproduction) |
| Exactness rule | what counts as exact reproduction: symbolic derivation, or rational identification at a stated precision and bound |

### Definition of done (skeleton, finalized at go)

| # | Item |
| --- | --- |
| 1 | Solver and auditor returns, scripts and data in the repository with `m8_10_` prefixes |
| 2 | Adversarial audit with its own method, per-claim verdicts, and a hunt for checks that cannot fail |
| 3 | Designer comparison against the frozen claims, every number stated; diagnostics recorded apart from the verdicts |
| 4 | Transcript audit of every agent |
| 5 | Method note `findings/m8_10_method_note.md`: equations first, equation-to-code map, audit record, manifests |
| 6 | Author package landed byte-identical to its hashes and compared for provenance only after the verdict, with the asymmetry stated |
| 7 | Doc sync (roadmap row and briefing), doc checker and roadmap linter exit 0, TASK REVIEW presented |

### What is frozen

From the designer's freeze the values stop moving on the author's side: no edit because a computed result disagrees with a frozen value. A disagreement is a result, and a claim that this text is defective goes to maintainer reproduction under #512. This binds the author, not the maintainers.

### Scheduling

Maintainer-run at maintainer pace. Nothing else in the column waits on it.

## DEVIATIONS LOG

(none)

## FINDINGS

(pending)
