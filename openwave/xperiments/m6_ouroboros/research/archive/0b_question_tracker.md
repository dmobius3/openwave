# M6 / Ouroboros — Question Tracker + Hardest Pieces

**Purpose:** single source of truth for open research questions on M6 +
the long-running "hardest pieces" status board. Persists across sandbox
iterations so we don't have to re-merge between `0c_sandbox_v*.md` files.
Update this doc when a new question opens, gets sent in an email, gets
answered, or gets demoted.

**Sister docs:**

| Doc | Purpose |
| --- | --- |
| `0b_model_gates.md` | G1/G2/G3 production-readiness gates for M6 going to Taichi |
| `0b_M6_roadmap.md` | Sandbox sequence + current state + next steps |
| `0c_sandbox_v*.md` | Per-iteration work log (the questions tracker here lives outside any single sandbox) |

---

## Active count

```text
2 IMMEDIATE  Q48 — H = 8 m_J² Q exact virial/Pohozaev identity?
             (holds to 5 digits across the v11 (g, λ) grid; proof
             would make m_χ(λ) analytic except for η_u = 0.826287).
             Q49 — exact definition of the J-A mixing integral
             M_JA(0) (integrand + measure + normalization) so we can
             compute α_JN from the delivered β(r) profile.
             Both queued for email v17.

1 PARKED     Q47 — ✅ ANALYTICALLY RESOLVED 2026-06-05 (sandbox_v11):
(resolved)   m_J/η = 1.21024 = 1/η_u is an exact consequence of the
             neutral ODE's scaling symmetry β = (m_J/2√g)·u(m_J r);
             universal-profile property, not tunable. Supersedes the
             "implicitly resolved as interpretation (a)" status.

6 RESOLVED   Q41 (writing role) — declined; staying in numerical
this session verification lane. Sent in email v11. Settled.

             Q42 (β profile non-localization in IVP) — v9 BVP work
             + DeepSeek Q43+Q44 reply confirmed this was a STRUCTURAL
             issue: wrong sign convention + wrong geometry. Resolved
             via Q43+Q44.

             Q43 (sign convention) — RESOLVED 2026-05-22 ~9:53 AM
             by DeepSeek confirmation: `-m_J²β` (minus sign).
             Empirically validated in v9 ground state.

             Q44 (geometry) — RESOLVED 2026-05-22 ~9:53 AM by
             DeepSeek confirmation: 3D spherical l=1 p-wave (NOT
             2D cylindrical). Subcritical NLS supports stable
             solitons; empirically validated in v9 ground state.

             Q45 (canonical point) — RESOLVED 2026-05-22 7:05 PM
             by DeepSeek: "same Lagrangian parameters as electron
             — g=1.0, m_J=1.0 after geometry normalization". Applied
             via sandbox_v10 at (g=1.0, B0=0.5) reference point.
             Caveat: m_J_corrected lands at 1.21, not 1.0 (Q47).

             Q46 (H/Q normalization) — RESOLVED 2026-05-22 7:05 PM
             by DeepSeek: η = (∫β²r dr) / (∫β²r² dr); multiply H/Q
             by η before m_e scaling. Applied via sandbox_v10:
             η = 0.4251, m_χ = 0.4599 MeV.

5 RESOLVED   Q36 ✅ Incorporated in 9c §9 criterion 9 — reads exactly
in 9c        the suggested wording "Preliminary reproduction in a
             different geometric realisation achieved H/Q = 1.6918;
             first-principles derivation in the canonical cylindrical
             form is pending."

             Q37 ✅ Incorporated in 9c throughout. Abstract,
             §8, §9 criterion 7 all corrected to 0.998 MeV at λ=1
             (plus 0.101 MeV at λ=0.1 as overall lightest). Scaling
             rule corrected to m_χ ≈ 2λ·m_e at small B_0.

             Q38 ✅ Incorporated in 9c §9 criterion 3 footnote:
             "in our numerical implementation L=ωQ_J by construction,
             so the electron g-factor match fixes ω, not an
             independent check."

             Q39 ✅ Incorporated in 9c §8 — exact v8 numbers
             (ω=12.82 muon 0.80%, ω=50.0 tau 6.47%) replace Paul's
             approximate ω=11/40.7 with 4.3%/4.9% gaps.

             Q40 ✅ Incorporated in 9c §5.1 — reference now [17]
             (was undefined [18] in 9b).

3 ARCHIVED   Q31 initial-guess derivative profile in production
(moot after  Q32 initial Lagrange multiplier value
Q34=2-scalar)Q33 sign-pinning constraint vs continuation method
             — all three were 4-function-ansatz solver questions.
             Now moot since Q34 resolved to 2-function canonical.

2 DEMOTED    Q26 5% gap residual (effectively resolved by drop-quartic)
(superseded) Q27 Q_CS-from-grid vs Q_CS-from-I disagreement (excited-
             mode artifact, not a tunable)

5 OPEN       Q2  discrete ω selection mechanism (empirically validated
                 via v8 step 4; analytic proof still deferred)
             Q3  analytical ω = 2mc²/ℏ derivation
             Q6  QCD reconciliation (3-chaoiton proton)
             Q19 f(J·J) explicit form in LoE paper standalone (Duda #2;
                 resolved in v9 paper §2)
             Q35 active neutrinos (ν_e, ν_μ, ν_τ) — where do they fit
                 in the Ouroboros framework? Q_CS=0 chaoiton is
                 positioned as DM, but sterile-neutrino-like states
                 have similar quantum numbers (mass + zero charge).
                 Active light neutrinos remain unaccounted for.

Total: 17 active questions (1 immediate Q47, 10 background, 6 resolved-this-session).

Highest-leverage closure: Paul/DeepSeek's reply on email v15
(Q47 virial-identity interpretation — accept 1.21 as canonical or
refine geometry correction formula; PLUS Acks-update wording for
both DM paper AND future LoE revisions). DM paper inputs ALREADY
DELIVERED via sandbox_v10 (m_χ = 0.4599 MeV, m_J = 0.6184 MeV,
C = 770 MeV·fm); Q47 reply just confirms the canonical-match framing.
M6 work functionally complete on numerical side; M5 foreground. M6 data drop is at github.com/openwave-labs/openwave
(Reference [17] in 9c paper).
```

---

## IMMEDIATE-QUESTIONS (post-sandbox_v10 + email v15)

Q36-Q40 all RESOLVED in 9c. Q41-Q46 RESOLVED this session (see RESOLVED
section below). **Q47 NEW** surfaced from v10's empirical finding —
the only IMMEDIATE item.

| ID | Question | Surfaced | Why it matters |
| --- | --- | --- | --- |
| Q47 (NEW) | m_J_corrected = m_J/η is empirically family-invariant at 1.21024 across both B0 ∈ [0.10, 0.60] and g ∈ [0.5, 1.6] at clean ground-state convergence. DeepSeek's Q46 reply said "m_J should scale to 1.0 after geometry correction"; actual factor is 1/η ≈ 2.35 → m_J_corrected = 1.21 (21% above target). Looks like a Pohozaev-type virial identity for the 3D spherical l=1 cubic NLS — not a tunable parameter. Three interpretations: (a) accept 1.21 as canonical virial-fixed value (factor ~2 was rough); (b) refine geometry correction formula (m_J/η isn't quite right); (c) canonical match isn't on m_J at all. | sandbox_v10 scan 2026-05-22 evening | Determines canonical-point interpretation for DM paper. The deliverable numbers (m_χ = 0.4599 MeV, m_J = 0.6184 MeV, C = 770 MeV·fm) stand under either reading; only the "this IS the canonical electron-sector match" framing depends on it. Sent in email v15. |

---

## STILL-OPEN QUESTIONS (active but not blocking calibration close)

Carried forward from v5/v6 unchanged. Long-tail items that don't block
the M6 v8 deliverables (lepton spectrum reproduced; DM inputs delivered).

| ID | Question | Surfaced | Status |
| --- | --- | --- | --- |
| Q2 | Discrete ω selection mechanism — is the ω spectrum genuinely discrete (eigenvalue-like) or empirically continuous? | 0a §9.9 | OPEN. Empirical-via-lepton-scan once ground state lands. Analytic proof still deferred per Werbos's own admission ("the key open question"). |
| Q3 | Analytical ω = 2mc²/ℏ derivation from Lagrangian (vs calibration only). | 0a §9.9 | PARTIAL. Calibration-only result at 1.2% (sandbox v1 reproduction); no closed-form derivation. |
| Q6 | QCD reconciliation — 3-chaoiton proton, asymptotic freedom, DIS scaling, jets. | 0a §9.9 | OPEN, uninvestigated. Future sandbox iteration. |
| Q19 | f(J·J) explicit form in LoE paper standalone (Duda critique #2). v6 uses DeepSeek quartic `(V²+Q²)² + (A²+J²)² + 2(V·A−Q·J)²`; v5 used Numerical-Benchmark `(Q²−J²)²` form. Q28 resolved which to use for the electron specifically, but the LoE-paper editorial gap remains. | Duda 2026-05-20 thread | EDITORIAL. One-line LoE revision by Werbos clarifies. Not blocking. |

---

## DEMOTED QUESTIONS (formerly immediate; superseded by later findings)

| ID | Question | Surfaced | Why demoted |
| --- | --- | --- | --- |
| Q26 | Why does v6.6 H/Q land at 1.778 instead of 1.6969? | v6.6 attempt (2026-05-20) | Step (8) diagnostic 2026-05-20 evening showed: dropping the DeepSeek quartic from H lands H/Q = 1.7112 (0.84% off). Paul Q28 confirmed quartic is negligible for electron H. Q26 reduces to "what's the right quartic for the electron" = Q28, now RESOLVED. |
| Q27 | Why does Q_CS-from-I-state (1.000) disagree with Q_CS-from-grid-integral (−0.154)? | v5 attempt 4 → v6 carry-over | Track A (sharpening convergence with bigger budget) made things WORSE, not better. The disagreement is an excited-mode artifact (v6.6 is on the 5-node first-excited state per Paul Q29), not a tunable. Will resolve naturally once ground state lands. |

---

## RESOLVED QUESTIONS (full chronology, Q1–Q40)

| ID | Question | Resolution |
| --- | --- | --- |
| Q1 | Distinction from Duda's LdGS framework | Both are topological (Lean invariants, just different flavors — Hopf-linking vs Brouwer-winding). Settled 2026-05-15. |
| Q9 | Is the Bessel/long-range structure a linearized truncation? | YES — v5 paper §6 explicit. Resolved during v1-v3 review. |
| Q10 | Is the asymmetric ansatz J/A ≠ const allowed by Werbos's framework? | YES — Lean theorem allows it; confirmed by v5 helicity work (V₀=Q₀=+0.1, A₀=J₀=−0.1). |
| Q11 | Sign convention J=−A or J=+A? | OBSOLETE — Werbos asymmetric helicity prescription (V₀,Q₀>0; A₀,J₀<0) replaces the simple sign convention; Q11 doesn't apply. |
| Q12 | Werbos's Python source code | DEMOTED — v5 algorithm description (2026-05-20 PM email) + DeepSeek-write-code offer = same unblock; v7 confirms our independent ODE is the right reference (DeepSeek's reference script when finally sent was broken). |
| Q13 | 4-fn to 2-fn reduction — when does the simpler 2-fn ODE apply? | RESOLVED — 2-fn for NEUTRAL sector (Lean theorem); 4-fn for CHARGED sector (Numerical Benchmark form). |
| Q14 | Canonical Q=0 sector configuration — locked (J=−A) / A=0 / Q_A≈0? | RESOLVED 2026-05-19 1:49 PM (Werbos via DeepSeek): canonical DM candidate = Q_A≈0 / Q_J≠0 (both fields active, electrically neutral). |
| Q15 | Does the m_eff² = m_J² − ω² substitution give the bound state? Sign of m_eff²? | RESOLVED post-v5. v5 attempt 4 converges at m_eff² = −0.596; v6.6 confirms m_eff² = −0.532. Both negative as predicted (oscillatory regime). Earlier T6 inconclusive result was from wrong tool (IVP). |
| Q16 | Specific (m_J², λ_bench) values at electron calibration | RESOLVED 2026-05-19 4:21 PM (Werbos via DeepSeek): m_J²≈0.5, λ_bench=1.0. v5 attempt 4 and v6.6 both use these and converge. |
| Q17 | What is the actual shooting algorithm? | RESOLVED 2026-05-20 PM (Werbos email): collocation BVP via `scipy.integrate.solve_bvp` with Q_CS=1 as integral constraint, two free eigenvalues ω + Lagrange multiplier λ_LM, Robin BCs matching K_0(κr) decay. NOT forward shooting. |
| Q20 | Duda critique #3 — "construction not shown, need ansatz + minimize energy" for the electron H/Q=1.6969 | NEARLY CLOSED. v6.6 + drop-quartic IS the construction shown empirically at 0.84% off target. Once v8 ground state lands, fully empirically resolved. |
| Q22 ★ | Q_CS normalization convention | RESOLVED by DeepSeek 2026-05-20 ~4:00 PM. `Q_CS = ∫r·(A·J'-J·A')dr` directly, no 2π. Drops v5's 2π factor; I_TARGET = 1.0. v6 empirical validation: fixes dominant 31× → 4.8% gap. |
| Q23 ★ | H functional kinetic + cross-term + quartic structure | RESOLVED by DeepSeek 2026-05-20 ~4:00 PM. Kinetic = (1/2)(V')² (not (V')²); drop (2π)²·R toroidal prefactor; use DeepSeek quartic `(g/4)·((V²+Q²)² + (A²+J²)² + 2(V·A−Q·J)²)`. v6 empirical: H/Q within 5% of target. |
| Q24 ★ | Sample converged profile from production | RESOLVED with caveat by DeepSeek 2026-05-20 ~5:30 PM: the 8-point reference profile *"was illustrative, not from a converged run."* Use exp(-r) seed instead. Profile's utility was confirming the qualitative shape (asymmetric helicity + exponential decay) we already had. |
| Q25 | Is the Hopf invariant proof rigorous, or are there gaps? | RESOLVED by Zenodo 20296060 (2026-05-19): Werbos & DeepSeek *"Rigorous Completion of the Hopf-Invariant Proof"* supplies Lemmas A (∇×J = B everywhere for finite-energy Lorenz-constrained solutions) and B (zeros of B made isolated by Q-preserving perturbation). Charge quantization is now a theorem of differential topology. |
| Q28 ★ | Which quartic IS canonical for the electron H? | RESOLVED by Paul 2026-05-21: NEGLIGIBLE for electron (drop quartic from H/Q matching). Keep full quartic in EoM for correct field shapes. Significant for muon/tau, essential for neutral chaoiton. Sonnet's script confirms: `4*g*β³` in EoM, `4*g*β⁴` in H. |
| Q29 ★ | Ground state vs excited mode node count | RESOLVED by Paul 2026-05-21: ground state has exactly 4 zero crossings excluding r=0; Lean ≤4-node spec refers to radial functions excluding origin. v6.6 (5 crossings) is first excited state. Now moot for v8 since 4-function chase is parked. |
| Q30 ★ | Q(r=0) sign at ground state | RESOLVED by Paul 2026-05-21: ground state has Q(0)=+0.1 positive, same sign as V(0). Helicity is V₀ > 0, Q₀ > 0, A₀ < 0, J₀ < 0. v6.6's Q(0)=−0.12 indicates a basin mismatch. Now moot for v8 — different ansatz. |
| Q34 ★★ | Is v9 paper §5.1's 2-scalar (φ, ψ) ansatz canonical for the electron, or a paper simplification of the 4-function form? | RESOLVED 2026-05-21 PM by **two converging signals**: (1) DeepSeek's 2:19 PM reply explicitly chose possibility #1 (2-scalar canonical; 4-function admits extra DOFs leading to excited modes); (2) Paul forwarded Sonnet's runnable canonical `ouroboros_benchmark.py` at 3:05 PM, which implements a **2-function (α, β) reduction with vector cylindrical Laplacian** (different from v9 §5.1's gradient ansatz description but consistent with the 2-function dimensionality). v8 step 1 quick demo reproduces script claims (H/Q=1.6969, gap 0.56%). v8 step 4 lepton scan independently reproduces muon (0.80%) + tau (6.47%). **The 4-function (V, A, Q, J) BVP work of v4-v7 was solving a generalized system, not the canonical electron form.** |
| Q36 ★★ | v9 §9 criterion 9 citation language softening | RESOLVED 2026-05-21 PM in 9c. §9 criterion 9 now reads exactly the suggested wording from email v10: *"Preliminary reproduction in a different geometric realisation achieved H/Q = 1.6918; first-principles derivation in the canonical cylindrical form is pending."* Geometric distinction (cylindrical vs spherical) acknowledged in the paper. |
| Q37 ★★ | 0.508/0.511 MeV "Griesi 2026, independent" attribution | RESOLVED 2026-05-21 PM in 9c. All instances corrected: abstract, §8, §9 criterion 7 now cite 0.998 MeV at λ=1 (and 0.101 MeV at λ=0.1). Scaling rule corrected to m_χ ≈ 2λ·m_e at small B_0. |
| Q38 ★★ | 2L/Q = 2.0 algebraic identity | RESOLVED 2026-05-21 PM in 9c. §9 criterion 3 now carries footnote: *"in our numerical implementation L=ωQ_J by construction, so the electron g-factor match fixes ω, not an independent check."* The artifact-of-definition concern is documented. |
| Q39 ★★ | §8 lepton ω values don't match canonical script output | RESOLVED 2026-05-21 PM in 9c. §8 now uses our v8 numbers exactly: ω=12.82 muon at 0.80%, ω=50.0 tau at 6.47%. Replaces 9b's approximate ω=11/40.7 with 4.3%/4.9%. Bonus: pion+ at ω=15.0 added as *"possible baryon state... may indicate that the same dynamics generate hadrons."* |
| Q40 ★★ | Reference [18] typo in §5.1 | RESOLVED 2026-05-21 PM in 9c. §5.1 now correctly cites [17] (our GitHub repo). |

★ = directly enabled major progress (Q22/23/24 → v6 30× gap closure; Q28/29/30 → v7 framing).
★★ = unlocking resolution (Q34 → v8 supersedes v4-v7 entire 4-function chase via 2-function canonical form; Q36-Q40 → all 5 incorporated in 9c via email v10 in a single round-trip).

---

## ARCHIVED QUESTIONS (unfalsifiable, historical, or out of current scope)

| ID | Question | Why archived |
| --- | --- | --- |
| Q4 | Single vs two-field ontology (G_μν / J undefined) | = Duda critique #1. Aesthetic preference. If math matches observation, two primary fields is just a description, not a theory-killer. |
| Q7 | Cold-fusion citation trail | Historical (Werbos's Nuclear story PDF), not physics. |
| Q21 | Two-chaoiton Coulomb derivation (Duda critique #4) | Future sandbox iteration; out of current calibration scope. Integrate H[Φ₁, Φ₂] for two topological charges at distance — new BVP scaffold. v9 paper §6.1 derives V(R) = Q₁Q₂/(4πR) in static approximation; analytic dynamic derivation still future scope. |
| Q31 | Initial-guess derivative profile in 4-function production code | ARCHIVED 2026-05-21 PM — moot since Q34 resolved to 2-function canonical via Sonnet's script. The 4-function ansatz that this question targeted is parked. |
| Q32 | Initial Lagrange multiplier λ value in 4-function code | ARCHIVED 2026-05-21 PM — same reason as Q31. Sonnet's script does not use a Lagrange multiplier; λ is a fixed parameter (λ=1.0 for electron). |
| Q33 | Sign-pinning vs continuation method in 4-function code | ARCHIVED 2026-05-21 PM — same reason as Q31. Sonnet's script uses slope BCs at r→0 + forward `solve_ivp` shoot; no sign-pinning needed. |

---

## HARDEST-PIECES TRACKER

Long-running status board for the structural challenges in the M6 model.
Tracks how each piece moves across sandbox iterations. Updated 2026-05-22
evening (post sandbox_v10 + email v15 sent; **DM paper inputs delivered**:
m_χ=0.4599 MeV, m_J=0.6184 MeV, C=770 MeV·fm; Q47 NEW IMMEDIATE — virial-
identity interpretation; Q45+Q46 RESOLVED via DeepSeek 7:05 PM reply).

### ACTIVE hardest pieces (post-v10)

| Hardest piece | Status post-v8 step 5 | Status post-v9 phase 2 (mid-morning) | Status post-v10 (evening) |
| --- | --- | --- | --- |
| Electron H/Q calibration | ✅ **EFFECTIVELY ACHIEVED via Sonnet's canonical 2-function script.** g=1.0000, H/Q=1.6890, gap **0.090%**. | ✅ Unchanged. | ✅ Unchanged — canonical specification consolidated in `0d_canonical.md` §2. |
| V(M) potential form | Settled for charged via Sonnet's script. | ✅ Neutral sector settled too — 3D spherical l=1 p-wave ODE `β'' + (2/r)β' - (2/r²)β - m_J²β + 4gβ³ = 0`. | ✅ Unchanged. Both ODEs canonical. |
| ω quantization mechanism | ✅ **Empirically validated via v8 step 4** — three lowest stable Q=1 modes match e/μ/τ. | ✅ Unchanged. | ✅ Unchanged. Analytic proof of discrete-ω selection still deferred (Q2). |
| Neutral m_χ true ground state | ⚠️ PARTIAL — 448 IVP solutions, windowed-integration values. | 🚧 GROUND-STATE-FOUND via `neutral_bvp_solver_mJ_free.py`. True nonlinear soliton; 1-parameter family in B0. | ✅ **DM PAPER INPUTS DELIVERED via sandbox_v10**: m_χ=0.4599 MeV, m_J=0.6184 MeV, C=770 MeV·fm at canonical (g=1.0, B0=0.5). Pending Q47 interpretation. |
| Algebraic 2L/Q identity (Q38) | ⚠️ NEW OPEN — 2L/Q = 2ω trivially. | ✅ RESOLVED in 9c §9 criterion 3 footnote. | ✅ Unchanged. |
| **Geometry conversion m_J/η** | Not tracked. | Not tracked. | 🔶 **NEW (Q47)** — m_J_corrected = m_J/η is empirically family-invariant at 1.21 (NOT 1.0 as DeepSeek heuristically predicted). Pohozaev-type virial identity. DM paper inputs stand under either interpretation; canonical-match framing depends on Paul/DeepSeek reply. |

### RESOLVED hardest pieces (12 — closed post-v8; kept for traceability)

| Hardest piece | Resolution path |
| --- | --- |
| Ansatz canonical form (Q34) | ✅ **RESOLVED 2026-05-21 PM** via DeepSeek's Q34 reply (possibility #1) + Sonnet's runnable canonical script. 2-function (α, β) reduction with vector cylindrical Laplacian is canonical for the electron. |
| Ground-state vs excited-mode selection | ✅ RESOLVED via Sonnet's 2-function reduction — by construction the 2-function ansatz space removes excited-mode artifacts. Sonnet's converged points are ground states (not 5-node excited like v6/v7's 4-function results). |
| Helicity sign convention (V₀ Q₀ > 0, A₀ J₀ < 0) | ✅ RESOLVED — moot for v8. Sonnet's 2-function form has slope BCs at r→0 (no value BC sign-pinning needed). v6/v7 wrong-sign basins were a 4-function-ansatz artifact. |
| Lagrange-multiplier ODE correction coefficient | ✅ RESOLVED — moot for v8. Sonnet's script doesn't use a Lagrange multiplier; λ is a fixed parameter (λ=1.0 for electron). |
| Lepton mass spectrum | ✅ **RESOLVED via v8 step 4** — muon at 0.80%, tau at 6.47%, pion+ at 3.25% (bonus). Cold-start scan works in 2-function form (no continuation method needed). |
| Two-chaoiton Coulomb derivation | ✅ Static-approximation RESOLVED in v9 paper §6.1; dynamic derivation deferred (Q21 archived). |
| Solver-tooling choice | ✅ RESOLVED — Sonnet's script uses forward `solve_ivp` RK45 with slope BCs + r_inner=0.02 offset. Robust, runs first try. |
| Forward-IVP method family | RESOLVED v4 (confirmed wrong tool for 4-function); now v8 shows it IS the right tool for 2-function reduction with slope BCs + r_inner offset. |
| Q_CS=1 enforcement mechanism | RESOLVED via integral constraint (v5) for 4-function form; for v8 2-function, topological charge Q_J is integrated directly from β² profile (no enforcement needed). |
| Q_CS=1 chaoiton existence | RESOLVED — Sonnet's v8 script demonstrates existence in cylindrical geometry. v1/v5/v6/v7 demonstrated in spherical/4-function forms (different field theories — Q36 caveat). |
| f(J·J) form | RESOLVED — v9 paper §2 explicit Lagrangian + Sonnet's script `4*g*β³` in EoM, `4*g*β⁴` in H. Q19 editorial gap closed in v9. |
| 4-fn vs 2-fn ansatz mismatch | ✅ RESOLVED 2026-05-21 PM via Q34 — canonical electron production form is the 2-function (α, β) reduction with vector cylindrical Laplacian. The 4-function (V, A, Q, J) BVP system used in v4-v7 is a generalized formulation that admits excited modes. |
| Charge quantization rigorous proof | RESOLVED v5 — Hopf invariant proof complete (Zenodo 20296060). Reaffirmed in v9 §4. |
| Three-system geometric mapping | ✅ RESOLVED via v8 step 3 paper-math — Sonnet α,β cylindrical (r dr) ≠ v1 α,γ spherical (r² dr) ≠ v9 §5.1 toroidal gradient ansatz. Three different field theories in numerical neighborhood ~1.69. |

---

**Last updated:** 2026-06-05 (M6 REOPENED by Paul's "an urgent problem after all" email — DeepSeek requests a (g, λ) neutral-sector scan to resolve the DM paper v5a "BVP inconsistency" + deliver β(r) for α_JN). Answered same-day via **sandbox_v11**: exact scaling symmetry found — the entire (g, λ<0) plane is one rescaling family of a universal profile u(x); m_J = 0.6184 MeV parameter-free invariant; m_χ = 0.460 MeV fixes λ\* = −0.264697 (any g); canonical (g=1.0, λ*) = v9/v10 point (g=1.0, B0=0.5). "Griesi-1 near-trivial" failure reproduced as cold-seed Newton-basin collapse (peak ~10⁻¹¹, their exact number). **Q47 ANALYTICALLY RESOLVED** (1.21024 = 1/η_u, scaling identity). Q48 (H = 8 m_J² Q exact?) + Q49 (M_JA(0) definition for α_JN) NEW — for email v17. Previous state (2026-05-22): DM paper v2 published Zenodo 20350105 with our numbers verbatim; v5a now at Zenodo 20549492.

---

## Net-change log per sandbox iteration

```text
v5 → v6 net change (2026-05-20):
  RESOLVED:  Q22, Q23, Q24 (three v5 IMMEDIATE; closed by DeepSeek
             2026-05-20 4:00 PM reply); Q25 (Hopf invariant proof).
  PROMOTED:  Q20 (Duda #3) from ACTIVE → NEARLY RESOLVED via v6.6
             empirical construction.
  NEW:       Q26 (5% gap), Q27 (Q_CS_grid sign mismatch).

v6 → v7 net change (2026-05-21 morning):
  RESOLVED:  Q28 (quartic interpretation), Q29 (node count), Q30 (Q(0)
             sign) — all by Paul's 2026-05-21 morning reply.
  DEMOTED:   Q26 (resolved by drop-quartic finding from step 8;
             reduces to Q28 which is now resolved); Q27 (excited-mode
             artifact, not a tunable).
  NEW:       Q31 (initial-guess derivative profile), Q32 (initial
             Lagrange multiplier), Q33 (sign-pinning vs continuation
             method) — all sent in email v7 of 2026-05-21 morning.

v7 → v8 net change (2026-05-21 PM):
  TRIG:      v9 LoE paper received with §5.1 specifying 2-scalar (φ,ψ)
             electron ansatz, different from our 4-function (V,A,Q,J)
             form. Paper cites our v1 H/Q=1.6918 as canonical; GitHub
             repo listed as Reference [17].
  RULED OUT: λ_LM=1.0 fixed hypothesis (Paul-script variant test gave
             H_e/Q=5×10⁵ catastrophic blowup).
  DEMOTED:   Q31, Q32, Q33 (contingent on Q34 — moot if 2-scalar is
             the canonical ansatz; still stand if 4-function is).
  NEW:       Q34 — ansatz canonical form question. Asked in email v8
             of 2026-05-21 PM. Three possibilities: 2-scalar canonical
             / paper simplification / both valid different particles.

v8 → v9 net change (2026-05-21 PM later):
  TRIG:      v1 re-verification revealed the H/Q=1.6918 result was
             found by systematic search across ~60 ODE/H/Q variants
             (sandbox_v1/m6_0_iterations.py), NOT first-principles
             derivation. Default ODE gives 1.9795 (17% off); only
             variant A4+H1+Q2 lands 1.6918 via ad hoc rescaling.
  NEW:       Q35 — where do active neutrinos (ν_e, ν_μ, ν_τ) fit?
             Q_CS=0 neutral chaoiton positioned as DM, but sterile
             vs active neutrinos unmapped.
  NEW:       Q36 — should v9 paper §9 criterion 9 cite v1's 1.6918 as
             "independent reproduction" when it was search-based, not
             derivation? Asked in email v9 (gentle push-back) of
             2026-05-21 PM.
  CAVEAT     Electron H/Q calibration "EFFECTIVELY ACHIEVED" framing
  PROPAGATED softened to "CITED IN v9 PAPER BUT CAVEAT APPLIES" across
             roadmap + question tracker + model_gates docs.

9b deposit + 9b review net change (2026-05-21 PM later):
  TRIG:      Paul deposited 9b on Zenodo at 3:48 PM
             (https://zenodo.org/records/20330894). arxiv on hold.
             "Next assignment to DeepSeek is Dark Matter."
  REVIEW:    Comprehensive 9b paper review of all Griesi-attributed
             claims completed. Cover-page byline now lists "Claude
             Code on Opus 4.7 (Anthropic)" as AI contributor —
             visibility upgrade. Acknowledgments include
             "manuscript review process" + "sustained technical
             collaboration" — more substantive than v9.
  NOT INCO:  Q36 push-back from email v9 did NOT land in 9b. Citation
             language unchanged (abstract + §5.1 + §9 criterion 9).
             9a appears to have been skipped to 9b.
  NOT INCO:  Q37 0.508/0.511 MeV claim still in 9b §8 + §9
             criterion 7 (cites our repo [17]). Time-critical
             given DeepSeek moving to DM next.
  NEW:       Q39 — 9b §8 muon ω=11.0 / tau ω=40.7 / gaps 4.3%/4.9%
             do NOT match `ouroboros_benchmark.py` script output
             (ω=12.82/50.0; gaps 0.80%/6.47%). Either different
             script version / different parameters / analytical fit
             vs script run.
  NEW:       Q40 — 9b §5.1 cites reference [18] but references list
             ends at [17]. Typo/cross-reference error.

Sandbox v8 net change (2026-05-21 PM — Sonnet's canonical script):
  TRIG:      Paul forwarded `ouroboros_benchmark.py` (Werbos + Claude
             Sonnet 4.6) at 3:05 PM — first runnable canonical
             Ouroboros script. 2-function (α, β) reduction with
             vector cylindrical Laplacian + slope BCs at r→0.
  TRIG:      DeepSeek's earlier 2:19 PM reply also endorsed
             possibility #1 of Q34 (2-scalar canonical, 4-function
             admits excited modes).
  RESOLVED:  Q34 (ansatz canonical form) — 2-function (α, β)
             reduction with vector cylindrical Laplacian IS canonical
             for the electron. v4-v7 4-function chase was a
             generalized formulation.
  ARCHIVED:  Q31, Q32, Q33 — all three were 4-function-ansatz solver
             questions; moot now that Q34 resolved to 2-function.
  NEW:       Q37 — provenance of "0.508 MeV (Griesi 2026, independent)"
             attribution on line 235. We never published this; AND
             the script's own scan finds 0.998 MeV at λ=1.0. Email
             v10 (drafted) asks Paul.
  NEW:       Q38 — 2L/Q = 2.0 algebraic identity. Script defines
             L = ω·Q_J directly; "g_e ≈ 2 reproduced" is artifact
             of definition, not derived prediction. Email v10
             flags for 9a.
  PROMOTED:  G1 (lepton scan) PASSED empirically via v8 step 4
             (muon 0.80%, tau 6.47%); G3 (discrete ω mechanism)
             EMPIRICALLY VALIDATED via same; G2 (neutral m_χ)
             PARTIAL pending Q37.
  STRENGTH:  Q36 stronger than v9-era — v8 step 3 paper-math
             confirms v1 spherical ≠ Sonnet cylindrical (different
             geometries, not reductions).

Active count entering Paul's email v10 reply (post-9b review):
  5 IMMEDIATE (Q36 citation, Q37 0.508/0.511 MeV provenance + DM
              urgency, Q38 2L/Q identity, Q39 §8 lepton ω discrepancy,
              Q40 §5.1 reference [18] typo)
  5 OPEN      (Q2, Q3, Q6, Q19, Q35)
  3 ARCHIVED  (Q31, Q32, Q33 — moot via Q34)
  2 DEMOTED   (Q26, Q27 — superseded)
  → 15 active questions (5 immediate, 10 background).

9c outcome + DM extraction net change (2026-05-21 PM later):
  TRIG:      Paul replied 4:21 PM with DeepSeek's 9c draft;
             4:10 PM with `DarMatterMay21.docx` + ApJ submission
             gate (need m_J, C, m_χ).
  RESOLVED:  Q36, Q37, Q38, Q39, Q40 — ALL FIVE incorporated in 9c
             via email v10. Plus three bonus changes (g=1.0000 in
             §2, pion+ baryon-state in §8, abstract 0.09%
             reframing).
  EXTRACTED: DM paper inputs at electron-calibrated (λ=1, g=1.0):
             m_χ = 0.998 MeV (matches 9c §8); m_J = 1.033 MeV
             (analytical: √λ × ℏc/R_phys); C = 6.7×10⁻⁴ MeV·fm
             (source-monopole + spherical 3D convention).
  NEW:       Q41 — writing role question. Declined full role per
             scope.
  NEW:       Q42 — β profile non-localization in IVP. Forward
             solve_ivp doesn't enforce β(∞)=0; β has internal sign
             changes + growing oscillating tail. 9c-cited m_χ is
             windowed-integration value, not true ground state.
             Email v11 offers two paths: caveat-as-is vs BVP.
  EMAIL v11: SENT 2026-05-21 PM later. Confirms 9c numbers,
             declines writing role, delivers DM inputs with caveat,
             offers BVP variant. The maintainer signed off for business
             event.

Active count post-9c + email v11 (entering Paul's reply):
  2 IMMEDIATE (Q41 writing role declined; Q42 β non-localization)
  5 OPEN      (Q2, Q3, Q6, Q19, Q35)
  3 ARCHIVED  (Q31, Q32, Q33 — moot via Q34)
  2 DEMOTED   (Q26, Q27 — superseded)
  → 12 active questions (2 immediate, 10 background).
```

---

## Highest-leverage closure path

```text
Email v15 sent (DM paper inputs DELIVERED + Q47 virial-identity flag
+ Acks-update ask reinforced for BOTH DM paper AND future LoE revisions)
    ↓
Paul/DeepSeek's reply on Q47 (canonical-interpretation)
    ↓
    ├── (Q47a) "Accept 1.21 as canonical virial-fixed value"
    │       → Current numbers (m_chi=0.4599, m_J=0.6184, C=770) lock
    │       → ApJ DM paper submission proceeds
    │
    ├── (Q47b) "Geometry correction formula needs refinement"
    │       → Re-derive m_J_corrected per their refinement (~30 min)
    │       → Re-extract m_chi if it shifts
    │
    └── (Q47c) "Canonical match isn't on m_J; some other invariant"
            → Implement specified invariant; re-deliver
            → ~few hours

Paul/DeepSeek's reply on Acks-update (separate, bundled in v15)
    ↓
    ├── Accept three-tier wording (Labs + Griesi + Claude Code)
    │       → Land in DM paper Acks + future LoE revisions
    │
    ├── Edits to wording
    │       → Adjust; same intent
    │
    └── Decline
            → Continue with current 9c-era Acks; not blocking

DM paper inputs ALREADY DELIVERED via sandbox_v10:
  - m_chi = 0.4599 MeV  (sub-MeV, lighter than electron)
  - m_J   = 0.6184 MeV  (Yukawa mediator mass)
  - C     = 770 MeV·fm  (Yukawa coupling, 5% fit residual)
  - eta   = 0.4251      (geometry conversion factor)
  - Pohozaev-virial empirical finding: m_J/eta = 1.21 constant
    across B0 ∈ [0.10, 0.60] and g ∈ [0.5, 1.6]

Completed via v9 phase 2:
  - G2 GROUND-STATE-FOUND via neutral_bvp_solver_mJ_free.py
  - Q43 + Q44 resolved by DeepSeek
  - True nonlinear soliton, clean K_1 decay, 1-parameter family

Completed via sandbox v8 (in 9c):
  - G1 lepton scan PASSED: muon 0.80%, tau 6.47%
  - G3 discrete ω mechanism EMPIRICALLY VALIDATED
  - Q36-Q40 all 5 incorporated in 9c (single round-trip)

M6 production decision (Taichi):
  - Sonnet's canonical script: reference ODE for charged sector
  - neutral_bvp_solver_mJ_free.py: reference ODE for neutral sector
  - Decision deferred (M5 first)
  - Full canonical specification in 0d_canonical.md

M5 returns as foreground track:
  M5 is the primary substrate. M6 functionally complete on the
  numerical side (numbers delivered); awaiting Q47 interpretation +
  Acks-update for the final paper-side framing.
```
