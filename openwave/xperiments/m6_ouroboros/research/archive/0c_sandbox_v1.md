# 2026-05-17 update — parallel M5/M6 build plan

This section supersedes §9.10. New plan triggered by:

- 2017 foundation paper reading (§11) — confirms framework maturity
- Duda's own acknowledgment that V(M) is "the most difficult" piece
- The matched-substrate observation — both M5 and M6 have V/f as the same bottleneck
- Cross-validation as a risk-mitigation strategy for the overall research thrust

## Decision context — both frameworks have the same hardest piece

Reframing the parallel-build question in light of §10.3:

| Framework | "Hardest piece" admission |
| --- | --- |
| Duda LdGS | "Choosing the details especially of potential is very difficult, will rather require PDE simulations" — Duda paper §III |
| Duda LdGS | "There is potential with minimum in this `diag(g, 1, δ, 0)` — getting EM+QM+GEM vacuum dynamics, and activating this potential especially to regularize infinite energy of e.g. charge. But finding its details seems the most difficult — could be like in Landau-de Gennes or slightly different." — Duda 2026-05-15 reply |
| Duda LdGS | "Usually potentials are effective — there might be even deeper ~anisotropic fluid, effectively described by such liquid crystal-like potential." — Duda 2026-05-15 reply |
| Werbos Ouroboros | f(J · J) "can be any nonnegative function" — 2017 paper. The specific form is a calibration choice. |

**This is the SAME bottleneck.** Both frameworks rest on the right choice of nonlinear potential. Both creators acknowledge they don't know the exact form. The numerical work on either side could give intuition for the other.

This is the strongest reason to build both: **diversification of failure modes on the shared bottleneck.**

### The plan — 4 steps

```text
   1. NumPy sandbox to reproduce Werbos's mass-frequency scan       (1-2 weeks)
                  │
                  ▼ pass → proceed; fail → drop M6 priority
                  │
   2. Decision on building M6                                       (1 day)
                  │
                  ▼ go
                  │
   3. Inform Werbos of decision + cc Models-of-Particles            (1 day)
                  │
                  ▼
                  │
   4. Scaffold M6 from current M5 vector substrate + parallel build (ongoing)
                  │
                  ▼
   Both M5 and M6 build proceed in parallel,
   cross-validating on shared observables
```

### Step 1 — NumPy reproduction of Werbos's mass-frequency scan

**What to reproduce** (from Particle Spectrum paper §2):

| Element | Werbos's value |
| --- | --- |
| Coupling | g = 1.0625 |
| Lagrange multiplier | λ = 1.0 |
| ω scan range | 1.0 to 80 in 31 steps |
| Integrator | RK45, rtol = 10⁻⁹, atol = 10⁻¹¹ |
| Domain | r ∈ [0.02, 15.0] code units |
| Grid | 500 grid points |
| Localization criterion | `\|A(r_max)\| + \|J(r_max)\| < 0.25` |
| Expected results | Electron at ω=1.0, muon at ω≈11.0, pion at ω≈13.0, tau at ω≈40.7 |
| Expected scaling | m ∝ ω^2.22 (near-quadratic) |

**Why NumPy first, not Taichi:**

| Reason | Detail |
| --- | --- |
| Direct comparison | Werbos's code is scipy solve_ivp (RK45); using same tooling enables direct cross-check |
| Speed of implementation | 1-2 weeks vs months for Taichi production |
| Decision-grade output | the scan either reproduces or doesn't — clear pass/fail Gate 0 |
| No M5 disruption | runs alongside M5.4 substrate refactor without splitting Taichi engineering time |

**Outputs**:

- Plot of m_predicted vs ω across the scan
- Localization-acceptance map across (g, ω) parameter space
- Verification of the m ∝ ω^2.22 scaling law
- Decision-grade verdict: framework reproduces / partially reproduces / doesn't reproduce

**Suggested location**: `openwave/xperiments/m6_ouroboros/research/scripts/m6_0_werbos_reproduction.py` (creates the M6 directory in a minimal-scaffold form even if M6 build doesn't go ahead).

### Step 2 — Decision criteria

| Reproduction outcome | M6 decision |
| --- | --- |
| Reproduces within 5% across all 3 leptons | **GO** — Werbos's numerics are real; M6 worth building |
| Reproduces for electron+muon but tau gap > 10% | **Conditional GO** — start scaffolding while investigating tau discrepancy |
| Fails at electron calibration | **NO-GO** — framework's foundational numerics don't replicate; drop M6 priority |
| Reproduction shows ω is continuous (not eigenvalue-spectrum) | **Conditional** — the 3 "fitted" ω's mean Ouroboros has the same parameter-count problem as Standard Model; still worth M6 sandbox but lower priority |
| Numerical instability / non-convergent | **Investigate before deciding** — could be implementation detail or fundamental issue |

### Step 3 — Inform Werbos

Short message to Werbos (cc Models-of-Particles list per group norms). Brief, specific, technical:

```text
- Acknowledge the Particle Spectrum paper directly addresses Q2 from our prior exchange
- State the reproduction result (qualitative — was the spectrum scan reproducible?)
- Announce decision on M6 scaffolding
- Reaffirm cross-validation pitch (same engine, shared observables)
- Re-invite engagement on architecture / design discussion at the openwave-labs repo
- Keep the 3 unaddressed questions visible (Q1 mechanism, Q3 analytical ω, Q4 — even though 2017 paper sharpens the answer)
```

Same voice / tone as the prior Werbos email. Goes to the same recipient list. Public commitment to M6 (if GO) demonstrates the platform-for-multiple-models pitch is real.

### Step 4 — Scaffold M6 from M5 vector substrate

If GO from step 2:

| Step | Source | Target |
| --- | --- | --- |
| 4a | Create `openwave/xperiments/m6_ouroboros/` directory tree | mirror M5's structure |
| 4b | Copy `lagrangian_engine.py` → `ouroboros_engine.py` | adapt for Vector(4) × 2 (A and J) instead of Vector(3) ψ |
| 4c | Adapt seeders for toroidal-poloidal ansatz | reuse M5's hedgehog seeders as template |
| 4d | Implement Lorenz constraint enforcer | new — not in M5 |
| 4e | Adapt `update_director_glyphs` for vector-field rendering | reuse with minor changes |
| 4f | Run Werbos's reproduction in Taichi at 65³ as M6 Gate 1 | matches what NumPy did, but in production engine |
| 4g | Implement Chern-Simons linking number kernel | M6.1 gating test |

**Important**: M6's Vector(4) × 2 substrate is LIGHTER than M5's matrix substrate (8 DoF/voxel vs 6 for M5 with the Lorenz constraints applied — but the constraint-enforcement makes it net comparable). The Taichi scaffolding can largely mirror M5 with substrate-type changes only.

## Parallel-build engineering reality

Once both M5 and M6 are scaffolded:

| Risk | Mitigation |
| --- | --- |
| Context switching cost | Time-box weeks — alternate weeks between M5 and M6 production work, not days |
| Cognitive load | Keep cross-validation observables identical across both (same Coulomb test, same g-factor test, same resonance protocol) — this reduces mental model count |
| Risk of incomplete both | Hard prioritization: M5 stays primary; M6 is secondary; if either is at risk of stalling, finish M5.7 first |
| Shared bottleneck on V/f | Treat as a feature: insights from one might inform the other; track findings in both research folders |

**Communicating the parallel build**:

- README / WELCOME: model-agnostic positioning is already in place (per the 2026-05-16 README rewrite)
- Once M6 is scaffolded, add to the Major Theoretical Contributions table (M5: Duda; M6: Werbos)
- Update memory entries to reflect both methods are in active build

### NSF X-Labs strategic note

If Werbos secures NSF X-Labs funding for his Ouroboros sensor architecture, having OpenWave's M6 scaffolded becomes immediately valuable:

- Demonstrates independent implementation of the framework
- Positions OpenWave as the open-source numerical platform underneath his sensor work
- Could open collaboration paths (his sensor + our numerical engine)
- Or competitive paths (our M5+M6 cross-validation as the better scientific platform)

Either way, having M6 scaffolded ahead of any X-Labs award outcome is strategically sound.

## Summary — the new plan

| What changed from §9.10 | Why |
| --- | --- |
| "Don't commit before M5.7" → "NumPy sandbox now" | 2017 paper + matched-substrate observation make M6 evaluation cheaper to start earlier |
| "Sequential build" → "parallel build after reproduction" | Duda's V(M) hedging + cross-validation risk-mitigation |
| "M6 as alternative" → "M6 as backup + cross-validator" | Failure-mode diversification on the shared V/f bottleneck |
| "Wait for Werbos reply" → "Reproduce his numerics ourselves" | Faster decision; reproduction is the strongest possible technical validation |

**Action gate:** start step 1 (NumPy sandbox) as soon as you have the engineering hours. Likely scheduled to interleave with the M5.4 substrate-refactor study work currently in progress, since they're complementary (one is sandbox NumPy, one is Taichi production planning).

## SANDBOX v1 RESULTS, ROADBLOCKS, NEXT STEPS

Hi Paul,

Thanks again for the Particle Spectrum paper — it directly addresses the lepton-hierarchy question (Q2) from our earlier exchange. I spent the last few days running independent
reproductions in a NumPy sandbox environment.

My background is in mechanical engineering, currently working on an ocean thermal conversion startup (Former NSF SBIR PI and Florida NextEra Incubator Alumni), but not physics. I use AI-assisted scripting (Claude Code on Opus 4.7) to iterate faster on the numerical interpretations. The math gets checked; the framing is mine.

### OUROBOROS: WHAT REPRODUCES

Electron calibration at (g = 1.0625, omega = 1.0, A_0 = B_0 = 0.1): I get H/Q = 1.6918 vs your published 1.6969 — 0.30% gap. The ODE form that works uses quarter-coupling on sources (consistent with your L = -F·F normalization, not the standard L = -F·F/4) and both A and J in the U(1) Noether charge. The L/Q = omega identity holds exactly.

### WHAT DOESN'T REPRODUCES

I scanned omega from 1 to 80 and get H/Q ~ omega^2.04, not your published omega^2.22. Mass predictions at your stated lepton omega values:

| Lepton | ω | Predicted (MeV) | Published (MeV) | Gap |
| --- | --- | --- | --- | --- |
| electron | 1.0 | 0.51 | 0.511 | 0.3% |
| muon | 11.0 | 73 | 105.7 | 31% |
| pion+ | 13.0 | 102 | 139.6 | 27% |
| tau | 40.7 | 1003 | 1777 | 44% |

Tried to close the gap: amplitude scans (H/Q is scale-free in amplitude), 96 combinations of (ODE form, H density, Q density) (caps at omega^2.04), the full 4-function ansatz from
your 2017 paper (worsens to omega^2.00, plus 1/r^2 instability), and integration-domain sensitivity (slope stable at 2.04).

The pattern I see: omega^2 is the framework-intrinsic scaling (matches classical oscillator E = (1/2)m·omega^2·A^2), and omega^2.22 corresponds to observed lepton mass ratios —
log(207)/log(11) = 2.21 and log(3478)/log(40.7) = 2.21.

Reproduction sandbox + plots + JSON summaries:
<https://github.com/openwave-labs/openwave/tree/main/openwave/xperiments/m6_ouroboros/research/archive/sandbox_v1>

### THE QUESTION

- Is there a stability or quantization criterion in the Ouroboros framework — beyond the localization check |A| + |J| < 0.15 — that selects omega = 1, 11, 13, 40.7 as discrete eigenvalues rather than as fitted points? In my scan, all omega values 1 to 60 give localized solutions equally well.

STILL OPEN (carrying over from previous message)

- Q1. Physical distinction from Duda's LdGS beyond the choice of topological invariant (Hopf-linking vs Brouwer-winding)?
- Q3. Is omega = 2mc²/h-bar derived analytically from the Lagrangian, or confirmed only via numerical chaoiton existence?
- Q4. The 2017 paper sharpens the "why two fields" answer (toroidal/poloidal mutual confinement). Jarek's deeper objection — single deeper field with A as derived connection, plus
Aharonov-Bohm — is still unaddressed.

### M6 DECISION PENDING SANDBOX RESULTS

OpenWave M6 production is pending clarification on the discrete-spectrum question. If I missed a selection mechanism, I want to fix the sandbox and proceed — the M6-Ouroboros-model directory is already structured to be built in parallel with the M5-Liquid-Crystal-model. The cross-validation value (running both your Ouroboros and Jarek's LdGS in the same numerical engine, same observables) is real and that's still the OpenWave pitch.

Open to design / architecture engagement at the repo from anyone in the group.

### IMMEDIATE-QUESTIONS (gates v2 sandbox)

```text
ID  | Question                                | First | Status at end of v1
    |                                         | surf  |
----|-----------------------------------------|-------|-------------------
Q2  | Discrete-spectrum mechanism — is there  | 0a    | ANSWERED EMPIRICALLY
    | a stability/quantization criterion      | §9.9  | (NEGATIVE) — v1 scan
    | beyond |A|+|J|<0.15 that selects        |       | sees all ω in [1,60]
    | ω={1,11,13,40.7} as DISCRETE            |       | localized equally;
    | eigenvalues, not fitted points? AND     |       | reproduction shows
    | do the predicted lepton masses match    |       | 31-44% lepton gaps
    | Werbos's claimed 4-6% gaps, or          |       | (not 4-6%); ω^2.22
    | larger?  [merges old Q5]                |       | = log(207)/log(11)
    |                                         |       | suggests fitted.
    |                                         |       | Awaiting Werbos
    |                                         |       | response to v1 email.
```

One gating question, sent to Werbos in the email above on 2026-05-17.

### STILL OPEN-QUESTIONS (load-bearing only)

```text
ID  | Question                                | First | Status at end of v1
    |                                         | surf  |
----|-----------------------------------------|-------|-------------------
Q3  | Analytical derivation of ω = 2mc²/ℏ?    | 0a    | PARTIAL — 1.2% via
    |                                         | §9.9  | R^phys calibration,
    |                                         |       | not analytical
Q6  | QCD reconciliation (3-chaoiton proton)  | 0a    | OPEN — Priority C
    |                                         | §9.9  | in v2 plan
```

### RESOLVED-QUESTIONS

```text
ID  | Question                                | Resolution
----|-----------------------------------------|------------------------
Q1  | Physical distinction from Duda's LdGS   | BOTH ARE TOPOLOGICAL.
    | beyond topological-invariant choice?    | 0a §2.1, §2.4, §9.4 and
    |                                         | Werbos's v2 paper §3
    |                                         | confirm Ouroboros uses
    |                                         | Chern-Simons linking =
    |                                         | Hopf invariant. Duda
    |                                         | uses Brouwer degree. The
    |                                         | distinction is invariant-
    |                                         | type + field structure
    |                                         | (matrix M vs two-vector
    |                                         | A,J), NOT topology-vs-
    |                                         | no-topology.
```

### ARCHIVED-QUESTIONS (framework-semantics / historical)

```text
ID  | Question                                | Why archived
----|-----------------------------------------|------------------------
Q4  | Single-deeper-field vs two-field        | Werbos 2017 explains
    | ontology                                | two-field via toroidal-
    |                                         | poloidal mutual
    |                                         | confinement. Duda's
    |                                         | "single deeper field"
    |                                         | preference is aesthetic,
    |                                         | not falsifiable physics.
Q7  | Cold-fusion citation trail              | Historical/credibility
    |                                         | question; not framework
    |                                         | physics. Tracking
    |                                         | separately if needed.
```

Active question count at end of v1: 1 IMMEDIATE + 2 STILL OPEN = 3 load-bearing.

### HARDEST-PIECES TRACKER (state at end of v1)

```text
Declared            | Source             | Resolution status at v1
--------------------|--------------------|--------------------------
V(M) potential form | Duda §III          | UNRESOLVED — Duda himself
                    | Duda 2026-05-15    |   admits "the most difficult"
                    |   email            |
f(J·J) form         | Werbos 2017 paper  | UNRESOLVED — Werbos says
                    |   ("any nonneg     |   f(s) = gs² in Calibration,
                    |   function")       |   but 2017 paper says any
                    |                    |   nonneg function works
ω quantization      | Werbos Spectrum    | UNRESOLVED — Werbos's own
  mechanism         |   §6.1 ("the key   |   "key open question";
                    |   open question")  |   sandbox v1 confirms ω is
                    |                    |   NOT discrete (continuous in
                    |                    |   the localization scan)
```

Pattern at end of v1: three hardest-pieces, all unresolved. V/f form bottleneck is shared between M5 and M6 — same physics problem, different math objects. ω-quantization is what sandbox v1 sharpened from "open question" → "empirically continuous in our scan".
