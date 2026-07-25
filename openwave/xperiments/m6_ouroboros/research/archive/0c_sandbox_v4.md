# 2026-05-18 evening — Sandbox v4: Lepton scan + awaiting Werbos's code

Triggered by Werbos's reply at 19:39 on 2026-05-18 (via DeepSeek; Paul is 78
and copy-pastes DeepSeek's drafts when tired). The reply addresses three of
the four immediate questions from v3, promises the Python source within 1-2
days, and gives an operational hypothesis for the discrete-ω selection
mechanism that v4 can test directly.

**Source:** Email thread continuation, "Intro & Nice to eMeet You" / reply
chain. The maintainer responded at 20:10 acknowledging the late hour, will resume
2026-05-19.

---

## What Werbos's reply answered (vs v3 end-of-file items)

```text
v3 ID | v3 status            | Werbos reply (2026-05-18 19:39)
------|----------------------|--------------------------------------------------
Q12   | NOT ON ZENODO        | PROMISED — "I'll share the Python source for the
      | (asked directly)     | Q=0 scan and the calibration run... Give me a day
      |                      | or two to get the files to you." Under QAGI LLC.
Q13   | PENDING — 4→2 fn     | EXPLAINED — Lorenz constraint locks (V,Q) and
      | reduction?           | (A,J) into pairs:
      |                      |   Δ_r V = Q  → Q determined by V
      |                      |   Δ_r A = J  → J determined by A
      |                      | Operational rule: use 4-fn for CHARGED calibration,
      |                      | 2-fn for NEUTRAL scan. ⚠️ But conflicts with the
      |                      | bosonization log's "A-field = 0" approach. → Q14.
Q2    | OPEN (3 sandboxes)   | ACKNOWLEDGED OPEN — "At this stage, the ω values
      | discrete ω           | are empirical." HYPOTHESIS provided: leptons are
      | selection?           | the lowest 3 stable modes in Q=1 sector. Analytic
      |                      | proof = "major mathematical undertaking", deferred.
      |                      | v4 lepton scan tests this hypothesis directly.
Q8    | PARTIAL (v3)         | (not directly addressed; will be settled by code)
Q9    | RESOLVED (v3)        | (no change)
```

### New question surfaced (Q14)

Werbos's DeepSeek reply says for the NEUTRAL chaoiton: `J(r) = -A(r)` AND
`Q(r) = -V(r)` (locked ansatz → single β). But his earlier bosonization
session log (May 18 4PM) said the neutral chaoiton has `A-field = 0, pure
J-field` (no locked relation, just A turned off).

```text
Description         | Source                | Our experiment
--------------------|-----------------------|-------------------------------
Locked: J=-A, Q=-V  | DeepSeek reply 19:39  | Sandbox v2 — FAILED (no soln)
A-field = 0         | Bosonization 4PM log  | Sandbox v3 — 23 SOLUTIONS ✅
```

These are different setups giving different empirical outcomes. v3 (A=0)
worked; v2 (locked) didn't. The DeepSeek summary may be approximating; the
actual code in his bosonization session probably used A=0. The Python source
(Q12) will resolve this definitively.

→ **Q14: Which is the canonical Q=0 description — locked ansatz or A=0?**

---

## v4 implementation plan

```text
Track         | Action                                | Gate | Blocked by
--------------|---------------------------------------|------|----------------
TRACK 1       | Run lepton scan with 2-function       | G1   | Not blocked —
(can run now) | Lean ODE on full (α, β) at varying ω  |      | can run today
              | g=1.0625, λ=1.0, ω scan [0.5, 50]     |      |
              | Test: are ω ≈ 1.0, ~12.78, ~40.7      |      |
              |   the lowest 3 stable modes?          |      |
TRACK 2       | Find neutral chaoiton ground state    |      | Not blocked —
(can run now) | Tune B0 to minimize H_code            |      | refines v3
              | Currently 0.508 MeV at B0=0.010;      |      |
              | smaller B0 → ground state             |      |
TRACK 3       | Receive Werbos's Python source        | G2   | Werbos 1-2 days
(blocked)     | Compare our v3 energy functional to   |      |
              | his; close calibration localization;  |      |
              | resolve Q14 (locked vs A=0)           |      |
TRACK 4       | If Q14 = locked: re-run v3 with the   | -    | Q14 answered
(branched)    | true locked ansatz; if Q14 = A=0:     |      | by code (T3)
              | proceed with v3 as canonical          |      |
```

### Track 1 — Lepton scan (highest-value next test)

Setup:

```text
ODE:       α'' + (1/r)α' − α/r² + ω²α = β
           β'' + (1/r)β' − β/r² + ω²β = α − λβ − 4gβ³
Params:    g = 1.0625, λ = 1.0  (electron-calibration values)
ω scan:    [0.5, 50] in steps of 0.5 (~100 points)
BCs:       α ~ A₀·r, β ~ B₀·r at origin (slopes; field=0 at r=0)
           α, β → 0 at r_max = 30
Stability: at each ω, check conjugate-point (Gelfand-Fomin) per
           Lean theorem §stability_condition
Outputs:   - localized solutions found at each ω
           - H/Q at each localized ω
           - mass predictions via m = H × m_e / H_electron
           - flag the 3 lowest-energy ω values → candidate leptons
Expected:  per Werbos's bosonization (May 18 4PM):
             ω = 1.0   → electron 0.511 MeV (calibration)
             ω = 12.78 → muon 106.8 MeV (1.1% gap)
             ω = 15.0  → pion+ 145.1 MeV (3.9% gap)
             ω ≈ 40-41 → tau (predicted)
Test:      if our 3 lowest stable ω values match leptons →
           Werbos's hypothesis (lepton = stable mode) confirmed
```

Effort: 1-2 days. New script `sandbox_v4/m6_v4_lepton_scan.py`.

**Calibration dependency (carried over from v3 ⚠️):** T1 anchors on ω=1.0
producing the electron baseline. v3 ended with H/Q = 1.723 (1.56% gap)
NOT localized. If T1 still can't localize at ω=1.0 with current energy
functional, we have no baseline to convert H_code → MeV for muon/tau
predictions. Two fallback paths:

```text
Path | Action                              | Pre-condition
-----|-------------------------------------|-----------------------
A    | Wait for T3 (Werbos's code) and    | Best if code arrives
     | use his energy functional with     | soon (1-2 days)
     | toroidal geometry factor (2π)²R    |
B    | Try alternative energy formulas    | If T3 delayed; lower
     | (add (2π)²R prefactor explicitly;  | confidence but
     | rescale H by 1/r at r=8; etc.)     | unblocks T1
```

Recommended: try Path B for a few hours first (cheap experiment);
if localization is still elusive, wait for T3.

### Track 2 — Neutral chaoiton ground state

Refine v3 sandbox to find the minimum-energy localized β solution by
shooting on B0 with bisection. Current: 0.508 MeV at B0=0.010. Expected: a
slightly lighter ground state below B0=0.010 — possibly approaching
Werbos's 0.003-0.015 MeV estimate, possibly not.

Effort: ~1 day. Either new script or augment `sandbox_v3/m6_v3_chaoiton.py`.

### Track 3 — Werbos's code (waiting)

Will give us:

- The exact energy functional with toroidal geometry factor → closes
  the calibration 1.56% gap (currently not localized)
- The canonical Q=0 description (Q14) — locked or A=0
- The full lepton scan code (independent cross-check on Track 1)
- The conjugate-point stability test implementation

ETA: 1-2 days per Werbos's reply.

---

## SANDBOX v4 RESULTS, ROADBLOCKS, NEXT STEPS

### T1 first run — 2026-05-19 (BLOCKED on wrong ODE form)

Built `sandbox_v4/m6_v4_lepton_scan.py` adapting the canonical Lean 2-fn
ODE from v3 into an ω-sweep. Quick scan run on ω ∈ {1, 5, 10, 12.78, 15,
20, 25, 30, 35, 40.7, 45} at g=1.0625, λ=1.0, (A0, B0) grid {0.05, 0.1,
0.5}². Total 99 integrations in 20.5s.

```text
Outcome              | Status
---------------------|--------------------------------------------------
Script works         | ✅ Localized solutions found at 10/11 ω values
                     |    (ω≥5; ω=1.0 marginally fails localization)
H scaling with ω     | ❌ DECREASING — opposite to Werbos's bosonization
                     |    H_full(12.78)/H_full(1.0) ≈ 0.04 (need 207 for muon)
Calibration H/Q      | ❌ Our H/Q at (ω=1.0, A0=B0=0.1) ≈ 1.04
                     |    Werbos: H/Q = 1.6969 (electron). 38% off.
v3 constant wrong    | ⚠️  H_CODE_ELECTRON_CALIB = 0.494 in v3 is wrong;
                     |    real H_spat at calib point ≈ 10.7 (close to
                     |    bosonization log's 11.4)
```

### Root cause: T1 is using the wrong ODE

Re-reading Q13 of Werbos's reply: the Lorenz constraint locks pairs
(V, Q) and (A, J), giving an operational rule:

```text
Sector     | ODE form used     | Status
-----------|-------------------|--------------------------------
Neutral    | 2-fn (Lean)       | ✅ v3 found 23 solutions
(α=0, ω=0) |                   |
Charged    | 4-fn (v5 paper)   | ❌ T1 used 2-fn — WRONG TOOL
(ω≠0)      | (φ, α, ρ, β)      |
```

The bosonization-session lepton scan that got ω=12.78 → muon 1.1% gap
ran Werbos's **earlier 4-function code**, *before* the May 18 Lean
reduction. The Lean 2-fn ODE is the structure-existence proof; it's
the reduced form for the Q=0 sector. Charged calibration + lepton scan
need the full 4-fn ansatz.

### What T1 needs to actually run

```text
Need                                  | Source / unblock
--------------------------------------|--------------------------------
Explicit 4-fn ODE (φ, α, ρ, β)        | (a) Werbos's Python source (Q12;
  with toroidal volume element and    |     ETA 2026-05-20, his promise)
  electron-calibration energy formula | (b) Manual extraction from
                                      |     v5 paper §6 + Numerical
                                      |     Benchmark docx (medium effort)
Correct H_CODE_ELECTRON_CALIB         | Implied by the above; v3's 0.494
                                      | is structurally wrong
Stability test (Gelfand-Fomin         | Needed even after ODE fix to rank
  conjugate-point)                    | "stable modes" per Werbos hypothesis
```

### Diagnostic kept

`sandbox_v4/diag_energy_functional.py` records the 2-fn H-scaling failure
and the H_full = H_spat + ω·Q identity. Keep for cross-check once 4-fn
code lands — the 2-fn form still applies to the Q=0 sandbox track.

### T1 decision

**HOLD T1.** Pivoting to T2 (neutral ground state via B₀ shooting) is
the right next move — it runs on the 2-fn Lean ODE we already have,
refines a v3 result rather than testing a hypothesis with the wrong
tool, and doesn't burn cycles guessing at Werbos's energy functional.

Resume T1 when:

```text
Trigger                              | Action
-------------------------------------|------------------------------
Werbos's code arrives (T3 unblocks)  | Adopt his 4-fn calibration
                                     | + lepton scan code as ground truth
OR ω=2026-05-20 +1d passes w/o code  | Manual 4-fn extraction from
                                     | v5 paper §6 + Numerical Benchmark
                                     | (slower, more uncertain)
```

---

### T2 first run — 2026-05-19 (NEGATIVE RESULT, corrects v3 overclaim)

Built `sandbox_v4/m6_v4_t2_neutral_ground.py`. Wide log-scan B₀ ∈
[1e-5, 1.0], 60 points, six λ values {0.1, 0.5, 1.0, 2.0, 5.0, 10.0},
r_max=30, plus golden-section bisection at λ=1.0.

```text
What we found              | What it means
---------------------------|----------------------------------------------
H is monotonic in B₀       | H grows like B₀^~3 (cubic-dominant scaling)
  across [1e-5, 1.0]       | for all λ. No minimum exists.
H → 0 as B₀ → 0            | "Ground state" of v3 was just trivial near-
                           | zero solution. Not a physical chaoiton.
Solution de-localizes at   | At λ=1.0: ~B₀=0.35. At larger λ: later.
moderate B₀ (~0.3-0.7)     | This is where the cubic term breaks the
                           | small-amplitude Bessel oscillation pattern.
Radial-node count 13-29    | Lean stability spec: ≤4 nodes. EVERY scan
                           | point violates this. The linear part is
                           | β''+β'/r-β/r²+λβ = J₁(√λ r) Bessel —
                           | inherently oscillatory.
```

**v3's "23 neutral chaoiton solutions" are artifacts.** Localization
was checked with `|β(r)| < 0.1 at r ≥ 8`, but at small B₀ the entire
oscillatory Bessel solution stays below 0.1 in absolute terms — not
because it decays, but because the amplitude itself is small. A
proper localization check needs EITHER node-count regularity (≤4 per
Lean) OR a decay-rate test (β² fall-off at large r relative to peak).

### What this implies for Q=0 / DM scan

The canonical Q=0 / neutral-chaoiton description is now in **triple
conflict**:

```text
Description              | Source                  | Empirical status
-------------------------|-------------------------|------------------
Locked: J=-A, Q=-V       | DeepSeek 19:39 reply    | v2 — FAILED
                         | (single β via locked    |   (no localized
                         |  ansatz)                |    solution found)
A-field=0, ω=0           | Bosonization 4PM log    | v3 claimed ✅;
  (α=0 reduction)        |                         | T2 DISPROVES
                         |                         |   (23 v3 solutions
                         |                         |    are artifacts)
Q_A≈0, Q_J≠0             | Bosonization log later  | UNTESTED — both
  (both fields active,   |   ("Q_A ≈ 0, Q_J ≠ 0")  | A and J fields
   EM charge zero but    |                         | active under
   nuclear charge ≠0)    |                         | constraint
```

→ This sharpens **Q14**: the canonical Q=0 description is one of THREE
candidates, not two. Werbos's code (Q12) is now the only clean path to
resolution. Without it, we cannot run a credible neutral chaoiton scan.

### T2 decision

**T2 CLOSED as negative.** Neutral-chaoiton existence in the canonical
Lean 2-fn form is not confirmed. v3's m_χ=0.508 MeV figure should be
treated as obsolete in roadmap and gates. Re-run only after Q14 is
settled.

### Sandbox v4 status (end of 2026-05-19 first-run)

```text
Track | Status     | Notes
------|------------|----------------------------------------------------
T1    | BLOCKED    | 2-fn ODE wrong tool for charged. Awaiting Werbos
      |            | code (Q12) OR manual 4-fn extraction.
T2    | CLOSED -   | v3 overclaim corrected. Q=0 chaoiton existence
      |            | now uncertain pending Q14 resolution.
T3    | WAITING    | Werbos's Python source, ETA 2026-05-20 (his
      | on Werbos  | promise). If slips → manual 4-fn extraction.
T4    | DEFERRED   | Q14 reconciliation rolls into T3 or 4-fn work.
T5    | DONE -     | Extraction complete (2026-05-19). 4-fn ODE
(new) | EXTRACTED  | found in Numerical Benchmark §3; calibration
      |            | anchor H/Q = 1.6969 confirmed. See section below.
```

### T5 — 4-fn ODE extracted (2026-05-19)

The Numerical Benchmark for the Ouroboros Lagrangian docx was **not
"wrong ODE form" as v3 notes claimed** — it's the canonical 4-fn ODE
source. The "wrong ODE" call was based on a misreading. Findings:

**4-fn benchmark ODE** (with mapping to our 2-fn notation):

```text
Mapping:  φ = V (A_0)   α = A (A_φ)   ρ = Q (J_0)   β = J (J_φ)
Convention: Δ_r f = f'' + (1/r)·f'    (NO -f/r² — toroidal, l=0)

Δ_r V  =  Q                                  (eq 1, Coulomb-type)
Δ_r A  =  J                                  (eq 2, Ampere-type)
Δ_r Q  =  V  +  m_J²·Q  +  λ·Q·(Q² − J²)     (eq 3)
Δ_r J  =  A  −  m_J²·J  −  λ·J·(Q² − J²)     (eq 4)
```

**Boundary conditions:**

```text
At r → 0:  V'(0) = A'(0) = Q'(0) = J'(0) = 0
           (VALUE BC, fields nonzero at origin — NOT slope BC)
           Calibration: V(0) = A(0) = Q(0) = J(0) = 0.1
At r → ∞:  all fields → 0, decay as modified Bessel K_0(κr)
```

**Energy functional** (toroidal volume element, NOT spherical 4πr²):

```text
H = (2π)²·R · ∫_0^∞ r dr [ (V')² + (A')² + (Q')² + (J')²
                          − V·Q + A·J
                          + (m_J²/2)·(Q² − J²)
                          + (λ/4)·(Q² − J²)² ]
```

**Charge functional** (Chern-Simons in toroidal form):

```text
Q_CS = (2π·R) ∫_0^∞ r dr [ A·(∂_r J) − J·(∂_r A) ]
Identity: L/Q = ω    →    2·L/Q = g-factor
```

**Calibration anchor (corrected):**

```text
Param                       | Value
----------------------------|------------------
g                           | 1.0625
λ                           | 1.0
ω                           | 1.0
V₀ = A₀ = Q₀ = J₀           | 0.1
H/Q (target, m_e/e_nat)     | 1.6875
H/Q (Werbos predicted)      | 1.6969 (0.56% gap)
Q_CS                        | 1 (integer, electron)
L/Q                         | ω = 1.000 → g_e = 2.000 (0.116%)
```

**v3 H_CODE_ELECTRON_CALIB = 0.494 is wrong**; the bosonization log's
H=11.4, Q=8.7 figures are from an earlier abandoned Compton calibration
iteration. The single authoritative anchor for charged sector is
**H/Q = 1.6969 at the canonical point.**

### T5 gaps (medium confidence — implementation risks)

```text
Gap                                  | Resolution path
-------------------------------------|------------------------------
ω doesn't appear in the static ODE   | Werbos's code (preferred) OR
  but calibration is at ω=1.0        | trial: replace m_J² with
                                     | (m_J² − ω²) in eq 3,4
f(s) two-form ambiguity              | Werbos's code OR trial:
  v5: f(s) = g·s²                    | benchmark form parametrizes
  benchmark: f = (m_J²/2)·s          | g via (m_J², λ) — need to
              + (λ/4)·s²             | derive the mapping
Toroidal radius R in (2π)²R          | Code-units choice; expect
                                     | R≈1 in normalized units
Numerical Benchmark vs               | Two ODE systems are STRUCTURALLY
  Lean theorem ODE                   | different. Lean: l=1 angular,
                                     | vector Laplacian, slope BC, ω².
                                     | Benchmark: l=0 toroidal, value
                                     | BC, static. Different ansätze.
```

The structural difference between Lean ODE and Benchmark ODE is the
*real* unknown — it's not that one is a "reduction" of the other.
They appear to be different ansatz choices for the chaoiton geometry.
Werbos's pre-Lean charged-scan code uses the Benchmark form; the Lean
theorem proves existence for a different (l=1 angular) ansatz.

### Decision needed (HUMAN)

Two paths forward:

```text
Path | Action                              | Risk
-----|-------------------------------------|----------------------------
A    | Implement 4-fn benchmark ODE now    | Medium — three open gaps
     | (T1 attempt 2 with correct ODE).    | (ω-in-static, f(s) form,
     | Calibrate against H/Q=1.6969 first. | R) need empirical tuning
     | If matches → run lepton scan with   | OR Werbos's code arriving
     | (m_J², λ) varied for ω-equivalents. | to settle them
B    | Wait for Werbos's code (ETA tomorrow| Low — settles all gaps
     | 2026-05-20). Use it as ground truth | cleanly. But delays
     | for both calibration + spectrum.    | progress by 1-2 days,
     |                                     | and code may not arrive
     |                                     | (Werbos is 78, may be
     |                                     | cautious about sharing
     |                                     | under QAGI LLC)
```

---

### Werbos clarifications — 2026-05-19 1:49 PM (via DeepSeek)

Werbos replied to the morning T2/T5 message with explicit answers to
all three implementation gaps and an endorsement of our Q14 finding.
Verbatim quote (DeepSeek-authored, Werbos-forwarded):

> Rodrigo — thank you for this. Finding that the 23 solutions don't
> survive a finer scan is exactly the kind of honest checking that
> separates real science from wishful thinking. I'm grateful.
>
> Three quick technical answers to your implementation gaps:
>
> **1. ω in the static benchmark ODE.** The benchmark document gives
> the static radial equations after separating time. The oscillation
> frequency ω is absorbed into the effective mass term. For the
> electron calibration at ω = 1.0, the mass parameter in the ODE is
> m_eff² = m_J² - ω². With the calibrated parameters, this effective
> mass is what drives the exponential decay. The static ODE in the
> benchmark is correct as written; you just need to use m_eff rather
> than m_J when you set up the solver.
>
> **2. Two parametrizations of f(s).** The v5 paper uses f(s) = g·s²
> as a simplified notation for the pure quartic case. The benchmark
> document uses the full form f(s) = ½ m_J² s + ¼ λ s². They are
> consistent when m_J² = 0 and λ = 4g. For the electron calibration,
> both forms should give the same result if you match the parameters
> correctly. Use the benchmark form — it's more general.
>
> **3. Toroidal major radius R.** The factor (2π)²R in the energy
> functional is a geometric factor from integrating around the torus.
> It cancels out when you compute the ratio H/Q because both the
> energy and the charge contain the same factor. For the mass-to-
> charge ratio, you can set R = 1 in scaled units and the ratio will
> be correct. The absolute mass scale is fixed by matching the
> electron mass, which absorbs any overall scale factor.
>
> **On the neutral chaoiton:** Your discovery that the three candidate
> definitions give different masses is important. This is not a bug —
> it means the neutral sector has a richer structure than we initially
> thought. The correct definition for the dark matter candidate is
> most likely the one where both fields are active but the net EM
> charge is zero (Q_A ≈ 0, Q_J ≠ 0). This is the configuration that
> would actually be produced in the early universe. The pure J-field
> case and the locked ansatz are limiting approximations.
>
> I'll get you the calibration code as soon as I can. In the meantime,
> try the m_eff fix for the electron and see if that closes your 1.56%
> gap.
>
> The structural argument in the ApJ paper — axion subsumption,
> no-free-parameter relic match, sensor architecture — holds
> regardless of the exact numerical value of m_χ. But I agree we
> should get m_χ on solid footing before submission. Once the
> canonical ODE is locked, your offer to run the full ground-state
> scan is exactly what's needed.

### Resolution summary

```text
Gap                  | Pre-Werbos state    | Post-Werbos resolution
---------------------|---------------------|------------------------
ω in static ODE      | Unknown — possibly  | EXPLICIT: replace m_J²
                     | absorbed into m_J²  | with m_eff² where
                     |                     | m_eff² = m_J² − ω²
                     |                     | Direct test: try at the
                     |                     | electron calibration
                     |                     | point, expect 1.56% gap
                     |                     | to close
f(s) two-form        | v5: f = g·s²        | EXPLICIT mapping:
                     | benchmark:          |   v5(g)  ≡  benchmark
                     |   f = (m_J²/2)·s    |   form with m_J²=0,
                     |       + (λ/4)·s²    |   λ = 4·g
                     |                     | Use benchmark form as
                     |                     | the canonical general
                     |                     | parametrization. For
                     |                     | g=1.0625: λ_bench=4.25
                     |                     | when m_J²=0.
Toroidal R           | Unspecified         | CANCELS in H/Q ratio.
                     |                     | Set R = 1 in scaled
                     |                     | units. Electron mass
                     |                     | anchor absorbs the
                     |                     | absolute scale factor.
Q14 (Q=0 canonical)  | 3-way conflict      | RESOLVED in principle:
                     | (locked, A=0,       | canonical DM candidate
                     |  Q_A≈0/Q_J≠0)       | is Q_A ≈ 0, Q_J ≠ 0
                     | from DeepSeek vs    | (both fields active).
                     | bosonization log    | Locked and A=0 are
                     |                     | limiting approximations.
                     |                     | "produced in the early
                     |                     | universe" — physical
                     |                     | early-cosmology framing
Code arrival         | "1-2 days"          | "as soon as I can" —
                     | (committed)         | no specific date;
                     |                     | implies slip from
                     |                     | original ETA
```

### Concrete next test (Werbos proposed)

> "Try the m_eff fix for the electron and see if that closes your
> 1.56% gap."

This is the single highest-leverage validation: implement the 4-fn
benchmark ODE with the m_eff² substitution and run at the canonical
electron calibration point. If H/Q lands at 1.6969 (within Werbos's
0.56% target gap), we have the correct charged ODE form.

### Implementation plan — v4 Path A → CONFIRMED GO

```text
Step | Action                                  | Outcome
-----|-----------------------------------------|----------------------
1    | Document this Werbos reply (DONE)       | ✅ verbatim above
2    | Build m6_v4_4fn.py — 4-fn benchmark     | new script
     | ODE with m_eff² = m_J² − ω², value      |
     | BCs f'(0)=0, R=1 in volume element     |
3    | KEY MILESTONE: calibrate against H/Q =  | success → unblocks
     | 1.6969 at (g, λ, ω, V₀=A₀=Q₀=J₀)        | everything; failure
     | = (1.0625, 1.0, 1.0, 0.1)               | → ask one targeted
     |                                         | Werbos clarification
4    | Charged lepton scan ω ∈ [0.5, 50]       | test "lowest 3 stable
     | with calibrated m_J² fixed              | = leptons" hypothesis
5    | Q_A ≈ 0 neutral scan (canonical DM      | true m_χ for ApJ
     | candidate per Werbos)                   | relic-abundance calc
```

The Q12 (code) ask is still LIVE but no longer blocking. If Werbos's
code arrives mid-implementation, we cross-check against our results.

### Update to question statuses

Several IMMEDIATE-QUESTIONS resolve with this reply (see updated
tables below). Q14 in particular moves from "3-way conflict" to
"resolved-in-principle: Q_A≈0/Q_J≠0 is canonical".

---

### T6 first run — 2026-05-19 PM (4-fn IVP, negative result)

Built `sandbox_v4/m6_v4_4fn.py` per Werbos's clarifications:
benchmark ODE with `m_eff² = m_J² − ω²` substitution, value BCs
`V(R_MIN)=A(R_MIN)=Q(R_MIN)=J(R_MIN)=0.1`, derivatives zero at origin,
toroidal R=1, ω-kinetic term included in energy density. Started
with RK45 → numerical stiffness → tried LSODA → hung → returned to
RK45 with `max_step=0.05` and a `BLOWUP_THRESHOLD=100` termination
event. Scan grid: m_J² ∈ {0, 0.5, 1, 1.5, 2, 3, 5, 10} × λ_bench ∈
{0.5, 1, 2, 4.25, 10}. All 40 combinations completed.

**Result: every combination blew up. No localized solution found.**

```text
Regime                   | What happens
-------------------------|---------------------------------------------
m_eff² < 0 (m_J² < ω²)   | Oscillatory growth (massive Bessel J_0
                         |   instability). Fields → 50-100 amplitude
                         |   before integrator terminates.
m_eff² = 0 (m_J² = ω²=1) | Q_CS = 0 exactly (degenerate). H = 9.5e6.
m_eff² > 0 (m_J² > ω²)   | Bound + non-bound modes mix. Generic
                         |   initial conditions excite the
                         |   unbound mode → exponential growth.
                         |   Q_CS NEGATIVE (-0.04 to -24), opposite
                         |   sign to target +1. Suggests A/J
                         |   winding orientation wrong.
```

### What this finding means

The 4-fn benchmark ODE is an **eigenvalue / shooting problem**, not
a plain IVP. Werbos's `V₀=A₀=Q₀=J₀=0.1` only gives a bound state at
ONE specific (m_J², λ_bench) calibrated point, not generically.
Working from arbitrary (m_J², λ_bench) with those initial values
gives saddle-point dynamics where the unbound mode dominates.

Equivalent statement: a chaoiton is the rare initial condition
within an 8-dimensional phase space `(V, V', A, A', Q, Q', J, J')`
that decays at infinity. Generic ICs blow up. The "calibration"
is finding that rare initial condition.

### Three paths forward

```text
Path | Approach                            | Effort | Risk
-----|-------------------------------------|--------|--------
A    | scipy.solve_bvp — enforce decay     | 2-4 hr | Medium
     | at r_max as boundary condition.     |        | (trivial
     | Standard tool for soliton-finding   |        | sol'n
     | ODEs. Needs sensible initial guess  |        | attractor)
     | + may need normalization to avoid   |        |
     | trivial solution.                   |        |
B    | Newton-style shooting on            | 1-2 hr | High
     | (m_J², λ_bench): residuals =        |        | (may not
     | (H/Q − 1.6969, Q_CS − 1) targets.   |        | converge)
     | Each step calls bounded IVP.        |        |
C    | Ask Werbos for the specific m_J²    | 0 hr   | Low —
     | and λ_bench values at calibration   | (one-  | his
     | (he gave us the m_eff² formula but  | line   | answer
     | not the m_J² number itself)         | ask)   | resolves
```

**Sequencing (per the maintainer, 2026-05-19 afternoon):**

1. Document this T6 finding now (this section)
2. Build T6→A (BVP solver) — first solo attempt
3. Review + document T6→A outcome
4. If A succeeds → use it for calibration; if not → move to T7 (B shooting)
5. Review + document T7 outcome
6. Ask Werbos with whatever specific question we've narrowed down to

The escalation principle: prove what we tried solo before going back
to him. His ApJ paper timeline is also our timeline; if we can do
the work without burning his cycles, that's better.

### What v3's 1.56% gap really meant

Werbos's "try the m_eff fix for the 1.56% gap" was about the LEAN
2-fn ODE (which v3 used). T6 shows the m_eff² substitution applied
to the 4-fn benchmark ODE doesn't trivially give bound states from
generic ICs. The 1.56% gap closing test would need to be on the
2-fn ODE specifically — which is now confirmed-wrong-tool for the
charged sector. So that specific test is moot. The benchmark ODE
calibration through BVP/shooting is the real anchor.

---

### T6→A — BVP solver attempt (Path A, 2026-05-19 PM)

Built `sandbox_v4/m6_v4_4fn_bvp.py` using `scipy.integrate.solve_bvp`
with eigenvalue parameter(s). Two formulations tried:

#### Formulation 1: single eigenvalue (m_J²) + 1 normalization (V_norm)

```text
BCs (9 total for 8-dim state + 1 param):
  V(R_MIN) = 0.1  (eats trivial sol'n attractor)
  V'(R_MIN) = A'(R_MIN) = Q'(R_MIN) = J'(R_MIN) = 0  (regularity)
  V(r_max) = A(r_max) = Q(r_max) = J(r_max) = 0      (decay)
Unknown: m_J²
```

Result: **converges, but to 2-fn subset.** A(r) and J(r) collapse
to zero identically. The (V, Q) sub-system finds a non-trivial
bound state via eigenvalue m_J²; (A, J) decouples spontaneously
because A=J=0 is a stable solution branch of the system.

```text
r_max=12  → m_J²_fit = 0.0933  m_eff² = -0.907  peak = 0.154  Q_CS = 0
r_max=24  → m_J²_fit = 0.9403  m_eff² = -0.060  peak = 0.103  Q_CS = 0
```

The r_max-dependence of the eigenvalue (0.09 → 0.94 as r_max
doubles) is a **Bessel-zero artifact signature** — the BVP is
catching oscillating J_0 solutions that happen to have a zero at
r_max, not genuine exponential-decay bound states. Negative m_eff²
in both cases confirms the bound-state regime is not reached.

**Formulation 2: two eigenvalues (m_J², λ_bench) + 2 normalizations
(V_norm, A_norm)**

```text
BCs (10 total for 8-dim state + 2 params):
  V(R_MIN) = 0.1, A(R_MIN) = 0.1  (eats trivial AND eats A→0 attractor)
  V'=A'=Q'=J' = 0 at origin
  V=A=Q=J = 0 at r_max
Unknowns: m_J², λ_bench
```

Result: **all 36 scan combinations fail with "singular Jacobian
encountered in collocation system."** The Jacobian is ill-conditioned
near the initial guess because the two BC sensitivities (∂/∂m_J²,
∂/∂λ_bench) are likely linearly dependent at the starting point.
2-eigenvalue BVP problems are fragile this way.

### T6→A diagnosis

```text
What we proved             | What it means
---------------------------|-------------------------------------------
BVP finds (V,Q) bound      | The Coulomb-like (V↔Q) sub-system has
state spontaneously        | bound-state structure; well-posed alone.
                           |
(A,J) collapses to 0       | Trivial-A,J branch is a stable solution.
in 1-eigenvalue mode       | Need explicit constraint to force them
                           | non-zero (Q_CS=1 isn't a boundary BC).
                           |
r_max-dependent eigenvalue | Solution is Bessel-zero artifact, not true
                           | exponential decay. m_eff² stays negative.
                           |
2-eigenvalue formulation   | Jacobian degeneracy at initial guesses.
fails with singular        | Could be fixable with smarter initial
Jacobian                   | guess or scaling, but high-risk.
```

The fundamental issue: **`Q_CS = 1` is an integral / topological
constraint, not a boundary condition.** scipy.solve_bvp can enforce
boundary values but not integral conditions. To pin Q_CS = 1, we'd
need either (a) a Lagrange-multiplier reformulation, (b) shooting on
(m_J², λ_bench, V₀, A₀, ...) with Q_CS as a residual, or (c) Werbos's
explicit calibration values.

### T6→A decision

**Path A is partial-success / partial-fail.** It validates the
eigenvalue structure (BVP found a bound state) but can't pin the
specific 4-fn chaoiton with Q_CS=1. Two productive directions:

```text
Direction                  | Effort   | Risk    | Note
---------------------------|----------|---------|-------------------
T7 — shooting on (m_J²,    | 2-3 hr   | High    | More flexible
λ_bench, V₀, A₀, Q₀, J₀)   |          |         | than BVP, can
with Q_CS=1, H/Q=1.6969    |          |         | use Q_CS as
as Newton residuals        |          |         | residual
T8 — ask Werbos for the    | 0 hr     | Low     | The specific
specific m_J² and λ_bench  | (1-line  |         | values fix
values at calibration      | email)   |         | everything
```

T6→A is documented and scripts archived. Proceeding to T7.

---

### T7 — Shooting refinement (Path B, 2026-05-19 PM)

Built `sandbox_v4/m6_v4_4fn_shoot.py`. Used "r_blowup distance" as a
continuous figure of merit: in a bound state, the IVP would integrate
to r_max without exceeding threshold; in a non-bound state, blowup
happens at some finite r. **Larger r_blowup ≈ closer to eigenvalue.**

**Coarse 2D scan** over (m_J², λ_bench) ∈ [0.5, 8] × [0.1, 10] with
V₀=A₀=Q₀=J₀=0.1, r_max=15:

```text
Observation                 | Detail
----------------------------|------------------------------------------
Maximum r_blowup = 6.44     | Well below r_max=15 → no bound state in
                            | this parameter window
m_J²=1.0 row is flat at 6.44| Degenerate: m_eff²=0 makes the cubic
                            | term irrelevant (λ_bench drops out)
r_blowup decreases as m_J²  | Larger mass → faster nonlinear blowup,
increases above 1.0         | not stabilization. Opposite of expected.
r_blowup also decreases as  | The λ_bench·(Q²-J²) coupling is
λ_bench grows               | destabilizing in this regime.
```

**Amplitude shoot** at fixed (m_J²=2.0, λ_bench=1.0), sweeping
V₀=A₀=Q₀=J₀=A_0 from 0.001 to 1.0:

```text
A_0        | r_blowup | Note
-----------|----------|----------------
0.001      | 7.77     | Largest seen
0.1        | 3.68     | Werbos's calib val
1.0        | 1.75     | Strong blowup
```

**r_blowup is monotone decreasing in A_0.** As A_0 → 0, r_blowup → ∞
(the trivial solution V=A=Q=J=0). There's no resonance / sweet spot
at non-zero amplitude where the bound state lives.

### T7 finding

For the **symmetric ansatz** V₀ = A₀ = Q₀ = J₀, there is **NO bound
state** at any tested (m_J², λ_bench) ∈ [0.5, 8] × [0.1, 10] and any
amplitude A_0 > 0.

This means **Werbos's calibration must use asymmetric initial values**
— V₀ ≠ A₀ ≠ Q₀ ≠ J₀ — OR substantially different (m_J², λ_bench)
values outside our scan window OR a different normalization / shooting
parameter we haven't identified.

### T7 decision: Werbos ask is now the highest-leverage action

T6 + T6→A + T7 collectively establish:

```text
What we've proven independently  | What we need from Werbos
---------------------------------|-----------------------------------
- 4-fn benchmark ODE implements  | - Are V₀, A₀, Q₀, J₀ all 0.1 at
  cleanly                        |   calibration, or do they differ?
- m_eff² = m_J² − ω² substitution| - Specific m_J² value at electron
  applied per Werbos's recipe    |   calibration point
- Eigenvalue structure confirmed | - Specific λ_bench value
  via BVP                        | - Shooting algorithm or BVP
- (V, Q) sub-bound state exists  |   formulation he uses for Q_CS=1
- Symmetric (V₀=A₀=Q₀=J₀)        | - Any normalization/decomposition
  ansatz has no bound state      |   trick that helps converge
```

The question for Werbos has narrowed from "send your code" to
"share these 2-4 specific numbers", which is much smaller ask and
respects his time / QAGI LLC constraints. If he answers, we're
unblocked instantly. If not, we have rich documented findings to
include in a slower full code request.

---

### Werbos / DeepSeek helicity clarification — 2026-05-19 4:21 PM

Werbos replied (via DeepSeek again) with the EXACT diagnosis we
needed. The symmetric ansatz forces Q_CS=0 by symmetry. The Q_CS=1
chaoiton requires **opposite helicity** between the (V, Q) and (A, J)
subsystems.

> "The Chern-Simons linking number Q_CS = 1 requires a specific
> relative sign (helicity) between the (V, Q) and (A, J) subsystems.
> The symmetric ansatz V₀ = A₀ = Q₀ = J₀ = 0.1 forces both subsystems
> to have the same sign structure, which locks Q_CS = 0. The chaoiton
> with Q_CS = 1 requires opposite helicity between the two subsystems
> — typically V, Q positive and A, J negative, or vice versa."

**Specific calibration values (verbatim from DeepSeek):**

```text
Field    Value    Sign
V(0)     0.1      positive
A(0)     0.1      negative  (or vice versa)
Q(0)     0.1      positive
J(0)     0.1      negative  (or vice versa)

Critical: V and Q must have the same sign; A and J must have the
same sign; the (V, Q) pair must have opposite sign to the (A, J)
pair. This asymmetry is what gives Q_CS = 1.

Parameters:
  m_J²        | ≈ 0.5 (mixed case) OR 0 (pure-quartic case)
  λ_bench     | 1.0 (mixed case) OR 4g = 4.25 (pure-quartic case)
  m_eff²      | m_J² − ω² = -0.5 (negative; oscillatory core,
              | correct for time-periodic solution)
  ω           | 1.0 (electron calibration)
```

**Shooting strategy (verbatim):**

> "The Q_CS = 1 chaoiton is found by a shooting method with the
> asymmetric initial signs, treating the decay rate at infinity as
> the target. The energy functional automatically selects the
> correct profile once the helicity is set correctly."

### Resolution of Q15, Q16, Q17

```text
Q15 | Does m_eff² substitution close   | RESOLVED. m_eff² = -0.5 IS the
     | electron calibration?           | calibration regime. Negative
                                       | m_eff² is correct ("oscillatory
                                       | core, time-periodic solution").
                                       | T6's expectation that m_eff²>0
                                       | gives decay was wrong.
Q16 | What specific m_J², λ_bench at   | RESOLVED. m_J² ≈ 0.5,
     | calibration?                    | λ_bench = 1.0 (mixed case).
                                       | V₀=Q₀=+0.1, A₀=J₀=-0.1.
Q17 | What shooting strategy?         | RESOLVED. Asymmetric helicity
(new)|                                  | + decay-rate-at-infinity
                                       | target. The opposite-sign
                                       | structure does the work that
                                       | a Q_CS=1 constraint would
                                       | otherwise impose.
```

### Implication for v4 implementation

The "fix" turns out to be a 4-character code change: flip A₀ and J₀
signs from +0.1 to -0.1. T6 + T7 explored ONLY the symmetric ansatz
and missed this because we treated (V₀, A₀, Q₀, J₀) as a single
common amplitude. Once asymmetric, the (V,Q) and (A,J) pairs are
fundamentally distinguishable in the cubic Q(Q²-J²) and J(Q²-J²)
terms — they were degenerate when V=A and Q=J.

### Next: re-run T6 + T7 with helicity

Specifically:

1. Update `m6_v4_4fn.py` to allow independent (V₀, A₀, Q₀, J₀) signs
2. Run single integration at (V₀=+0.1, A₀=-0.1, Q₀=+0.1, J₀=-0.1,
   m_J²=0.5, λ_bench=1.0, ω=1.0)
3. Check: does Q_CS land at +1? Does H/Q land at 1.6969?
4. If yes: run full lepton scan (vary ω) for the muon@1.1% test
5. If no: report specific failure mode and escalate

The Werbos paper draft v2 (Dark_Matter_in_Universe_v2.txt) is being
reviewed in parallel (paste-ready review delivered to the maintainer
2026-05-19 4:30 PM ish; user-edited and to be sent imminently).

---

### T8 — Asymmetric-helicity re-runs (2026-05-19 PM)

Built/updated three scripts to support Werbos's helicity recipe
(V_norm=+0.1, A_norm=-0.1, Q via init guess +0.1, J via init guess
-0.1). Three scripts updated to accept --helicity argument.

**T8a — IVP at canonical (m_J²=0.5, λ_bench=1.0):**

```text
Symmetric V₀=A₀=Q₀=J₀=+0.1   | Q_CS = 0 (T6 baseline)
Asymmetric Werbos helicity   | Q_CS = 4483 (non-zero! helicity works)
                             | But: tail=71, IVP still blows up,
                             | Q_CS is way off target 1.
```

**Conclusion: helicity is necessary (Q_CS goes from 0 to non-zero
when symmetry breaks) but not sufficient for IVP convergence.**

**T8b — Amplitude shoot with Werbos helicity at (m_J²=0.5, λ=1):**

```text
A_0          | r_blowup (with helicity) | vs symmetric (T7)
-------------|-----------|--------------|------------------------
0.001        | 8.70      |              | (best vs sym 7.77)
0.1          | 4.12      |              | (vs sym 3.68; +12%)
1.0          | 1.96      |              | (vs sym 1.75; +12%)
```

Helicity gives ~10-15% improvement in r_blowup at each A_0, but
r_blowup is still monotone in A_0 (no resonance). No bound state
found by amplitude shooting alone.

**T8c — 1-eigenvalue BVP with asymmetric helicity initial guess:**

Updated `m6_v4_4fn_bvp.py` to: (a) revert to 1-eigenvalue (m_J²
only, λ_bench fixed) to avoid the 2-eigenvalue singular Jacobian;
(b) add helicity-aware initial guess. Drop A_norm constraint
(state_dim + 1 param = 9 BCs requires only V_norm as the extra).

Scan over m_J²_init ∈ {0.5, 1.0, 5.0, 10.0, 50.0, 100.0}:

```text
m_J²_init | m_J²_fit  | A(0)    | Q(0)    | Q_CS   | Notes
----------|-----------|---------|---------|--------|---------------
0.5       | -41.4     | -0.000  | -2.90   | 0      | (V,Q) bound,
                                                   | a-J collapsed
1.0       | (FAIL)    | -       | -       | -      | max nodes
5.0       | 5.51      | 0.000   | -0.021  | 0      | same (V,Q)
                                                   | bound state
10.0      | 5.51      | -0.000  | -0.021  | 0      | converges to
                                                   | same eigenvalue
50.0      | -1.82     | 0.0003  | -0.306  | 0      | A,J non-zero
                                                   | but tiny; H/Q
                                                   | numerical noise
100.0     | (FAIL)    | -       | -       | -      | max nodes
```

The BVP converges to the same **(V, Q) bound state at m_J² = 5.51**
regardless of asymmetric initial guess. Even though we initialize
A = -V·exp(-r) and J = -Q·exp(-r), the iteration drives A → 0 and
J → 0. The asymmetric guess buys nothing because there's no BC
forcing A to stay finite at the origin.

### T8 finding

**Helicity is necessary but not sufficient for a Q_CS=1 chaoiton.**
The Q_CS=1 bound state requires either:

```text
Option A | A topological / integral constraint (Q_CS = 1 added as a
         | residual to the BVP/shooter). solve_bvp doesn't natively
         | support this; needs custom Newton on shooting parameters
         | with Q_CS as residual.
Option B | A second normalization (A_norm) that forces A non-zero,
         | with a second eigenvalue (λ_bench) to make BCs consistent.
         | We tried this in T6→A's 2-eigenvalue mode — fails with
         | singular Jacobian universally across our scan.
Option C | Better initial guess that's already in the Q_CS=1 basin,
         | so BVP iterates within it rather than escaping. Likely
         | requires Werbos's actual converged profile shapes
         | (V(r), A(r), Q(r), J(r) tabulated), not just BCs.
Option D | Werbos's actual shooting code with whatever parametrization
         | he uses ("treating the decay rate at infinity as the
         | target" — but on which shooting variable?).
```

### T8 decision

T8 is **CLOSED with partial progress**:

```text
Today's net advance              | What we proved
---------------------------------|-----------------------------------
Helicity matters                 | Q_CS depends on sign structure
                                 | (0 in sym, non-zero in asym).
                                 | Confirms Werbos/DeepSeek
                                 | helicity recipe is correct.
But not sufficient               | All formulations we tried
                                 | (IVP, BVP-1, BVP-2, shooting)
                                 | hit walls.
(V, Q) subset bound state real   | BVP repeatedly finds it at
                                 | m_J²=5.51, λ_bench=1.0,
                                 | regardless of initial guess.
                                 | This is genuine sub-physics.
```

---

### T9 — Newton-residual / two-stage shooter (2026-05-20 PM, DEFINITIVE NEGATIVE)

Built `sandbox_v4/m6_v4_4fn_newton.py` after T8 closed with "helicity
necessary but not sufficient." Design: two-stage optimizer over 6
variables `(m_J², λ_bench, |V₀|, |A₀|, |Q₀|, |J₀|)` with helicity signs
locked in evaluation wrapper (`V₀, Q₀ → +`; `A₀, J₀ → −`). Stage 1
`differential_evolution` global search on a smooth blow-up-distance
cost to find the bound-state basin. Stage 2 `Nelder-Mead` local polish
on calibration cost `(H/Q − 1.6969)² + (Q_CS − 1)² + tail² + collapse
penalties`. Pre-stage 2D grid scan at Werbos's exact amplitudes
`|V|=|A|=|Q|=|J|=0.1`.

**Pre-scan finding (324 points, m_J² ∈ [0.05, 10], λ ∈ [0.01, 10]
log-spaced, r_max=15):**

```text
count(r_reached ≥ r_max - 0.5) :  0 / 324
count(max_peak < 5)            :  0 / 324
max r_reached observed         :  7.241  (at m_J²=1.13, λ=0.01)
peak at that point             :  95.0   (catastrophic blow-up)
```

**Wider-scan finding (576 points, m_J² ∈ [0.01, 50], λ ∈ [0.001, 50]
log-spaced, r_max=30):**

```text
count(r_reached ≥ r_max - 0.5) :  0 / 576
count(max_peak < 5)            :  0 / 576
max r_reached observed         :  8.062  (at m_J²=0.85, λ=0.001)
peak at that point             :  98.4
```

**Stage 1 DE 6-variable search (popsize=20, maxiter=50, bounds m_J²,λ ∈
[0.01, 10], amps ∈ [0.001, 1.0]):**

```text
nfev               :  6528
best cost          :  493.4   (lower bound = 0 if r_reached = r_max)
best r_reached     :  8.398   (vs r_max = 15)
best θ*            :  m_J² = 0.404, λ_bench = 0.015,
                      V₀ = +0.196, A₀ = -0.026,
                      Q₀ = +0.006, J₀ = -0.004
at θ*              :  Q_CS = 58662 (target 1), peak = 92.9,
                      tail = 146.4, H = 1.2 × 10⁷
```

The 6-variable search slightly improved r_reached (8.4 vs 7.2 fixed amps)
but found a delayed-blow-up regime, not a bound state. Allowing
asymmetric amplitude *magnitudes* does not unlock the chaoiton.

### T9 definitive finding

Across **six independent attempts** (T6 IVP grid, T6→A 1-eigenvalue BVP,
T6→A 2-eigenvalue BVP, T7 amplitude shoot, T8 helicity re-runs, T9
pre-scan 324 pts + DE 6-var + wider-scan 576 pts), no parameter
combination yields a localized solution with `peak < 5` reaching
`r_max` under value-BC origin forward IVP. The Q_CS=1 chaoiton is
**unreachable by forward-IVP shooting from origin BCs** within any
parameter range we have explored.

### Interpretation

Possible explanations, ranked by likelihood:

```text
A | Werbos's "shooting method" is not forward-IVP from origin. He may   | most
  | integrate inward from r=∞ with K_0(κr) asymptotics, or use a        | likely
  | multi-shot connection at intermediate r, or fold Q_CS=1 in as a
  | Lagrange multiplier inside a BVP collocation.
B | The Q_CS=1 chaoiton has very narrow basin in 6-var space — our DE   | medium
  | popsize=20 and 6500 evals didn't find it. A much larger search
  | budget OR an analytical preconditioner is needed.
C | The 4-fn benchmark ODE as we transcribed it has a subtle structural | medium
  | error (sign in the (m_J²/2)·(Q² − J²) term, missing -f/r² piece for
  | l ≠ 0 ansatz, wrong V·Q vs +V·Q cross term, etc.) that systematically
  | prevents bound states from existing.
D | The chaoiton truly does not exist as a static r-dependent solution; | unlikely
  | it requires time-periodic ansatz coupling that the static benchmark
  | ODE doesn't capture.
```

(C) is testable by re-deriving the ODE from the v5 Lagrangian (deferred —
high-effort, low confidence we'd find an error Werbos missed). (D) is
contradicted by Werbos's claim of a working calibration. (A) is by far
the most productive next direction: ask Paul for algorithmic specifics.

### T9 decision

T9 is **CLOSED-NEGATIVE**. The systematic, multi-method evidence is
strong enough to justify a sharp escalation to Paul: we have exhausted
the natural solo approaches. The next email is a 1-screen note with the
six-attempt summary and an unambiguous algorithmic question.

### Next action — Werbos email v3

The maintainer's email (~5:00 PM 2026-05-19) asked for:

1. Specific calibrated (m_J², λ_bench) values [Werbos's 4:21 PM
   reply already answered: m_J²≈0.5, λ_bench=1.0]
2. Asymmetric helicity prescription [DONE — Werbos answered]
3. Shooting algorithm or BVP formulation for Q_CS=1 [STILL OPEN]

The third question is the one that remains. With today's review of
his v2 paper draft, and the documented evidence that helicity alone
doesn't crack it, we expect either:

```text
Likely Werbos response | Action on receipt
-----------------------|------------------------------------------
Specific algorithm     | Implement, anchor, run scans
description (e.g.      |
"add a 5th BC of       |
Q_CS=1 via integral")  |
Tabulated profile shape| Use as warm-start for BVP; should converge
(V(r), A(r), Q(r),     | within the right basin
J(r) on a grid)        |
Full code              | The gold standard; full unblock
"Will get back"        | Continue with parallel work (M5.4 etc.)
```

---

### IMMEDIATE-QUESTIONS (gates v5)

```text
ID  | Question                                | First | Status post-Werbos
    |                                         | surf  | 2026-05-19 1:49 PM
----|-----------------------------------------|-------|----------------------
Q12 | Werbos's Python source                  | v2    | LIVE but DEMOTED —
    |                                         |       | "as soon as I can"
    |                                         |       | (no committed date).
    |                                         |       | No longer blocking
    |                                         |       | progress; useful
    |                                         |       | for cross-check.
Q13 | 4-fn to 2-fn reduction                  | v3    | RESOLVED in
    |                                         |       | principle. v5
    |                                         |       | g·s² ≡ benchmark
    |                                         |       | with m_J²=0, λ=4g.
    |                                         |       | Use benchmark form.
Q14 | Canonical Q=0: locked / A=0 / Q_A≈0?    | v4    | RESOLVED. Canonical
    |                                         |       | DM candidate is
    |                                         |       | Q_A≈0 / Q_J≠0
    |                                         |       | (both fields active,
    |                                         |       | EM-neutral). Locked
    |                                         |       | and A=0 are limiting
    |                                         |       | approximations.
Q15 | Does m_eff² = m_J² − ω² substitution    | v4    | TESTED via T6 →
    | close electron calibration H/Q to       |       | INCONCLUSIVE.
    | 1.6969?                                 |       | Generic IVP from
    |                                         |       | V₀=A₀=Q₀=J₀=0.1
    |                                         |       | blows up across
    |                                         |       | full (m_J²,
    |                                         |       | λ_bench) scan.
    |                                         |       | The fix may work
    |                                         |       | but only at the
    |                                         |       | specific eigenvalue
    |                                         |       | which IVP can't find
    |                                         |       | from generic start.
Q16 | What specific m_J² and λ_bench does     | v4 PM | NEW — surfaced by
    | Werbos use at the electron calibration  |       | T6 IVP-blowup
    | point? (He gave us m_eff² formula but   |       | finding. Likely
    | not the m_J² number itself.)            |       | resolves T6 if
    |                                         |       | answered. Ask
    |                                         |       | scheduled for
    |                                         |       | after T6+T7 attempt.
Q17 | What is Werbos's shooting algorithm     | v4 PM | T9 CLOSED-NEGATIVE
    | for Q_CS=1? Forward-IVP from origin     | 5-20  | (2026-05-20).
    | does not reach a bound state across     |       | 6 independent
    | (m_J², λ_bench, V₀, A₀, Q₀, J₀) at      |       | attempts found NO
    | any tested values.                      |       | bound state with
    |                                         |       | peak<5 reaching
    |                                         |       | r_max. Werbos's
    |                                         |       | algorithm cannot be
    |                                         |       | forward-IVP. Need
    |                                         |       | concrete algorithm
    |                                         |       | description (most
    |                                         |       | likely: backward
    |                                         |       | from K_0 asymptotic,
    |                                         |       | or 2-eigenvalue BVP
    |                                         |       | with specific init
    |                                         |       | profile shape, or
    |                                         |       | Q_CS Lagrange
    |                                         |       | multiplier in a
    |                                         |       | collocation BVP).
```

### STILL OPEN-QUESTIONS

```text
ID  | Question                                | First | Status at start of v4
    |                                         | surf  |
----|-----------------------------------------|-------|----------------------
Q2  | Discrete ω selection mechanism          | 0a    | OPEN — Werbos says
    |                                         | §9.9  | empirical for now;
    |                                         |       | hypothesis = leptons
    |                                         |       | are lowest 3 stable
    |                                         |       | modes. v4 Track 1
    |                                         |       | will test directly.
Q3  | Analytical ω = 2mc²/ℏ derivation?       | 0a    | PARTIAL — only via
    |                                         | §9.9  | calibration (1.2%)
Q6  | QCD reconciliation (3-chaoiton proton)  | 0a    | OPEN — uninvestigated
    |                                         | §9.9  |
Q8  | Neutral m_χ true ground state           | v2    | PARTIAL — 0.508 MeV;
    |                                         |       | Track 2 will refine
```

### RESOLVED-QUESTIONS (carries forward from v3)

```text
ID  | Question                                | Resolution
----|-----------------------------------------|----------------------
Q1  | Distinction from Duda's LdGS            | Both are topological
Q9  | Linearized Bessel truncation?           | YES (v5 §6)
Q10 | Asymmetric ansatz J/A ≠ const?          | YES (Lean theorem)
Q11 | Sign convention J=−A or J=+A?           | OBSOLETE
```

### ARCHIVED-QUESTIONS

```text
ID  | Question                                | Why archived
----|-----------------------------------------|----------------------
Q4  | Single vs two-field ontology            | Aesthetic preference
Q7  | Cold-fusion citation trail              | Historical, not physics
```

### HARDEST-PIECES TRACKER

```text
Hardest piece     | Source       | v1 | v2 | v3 | v4 | Status start v4
------------------|--------------|----|----|----|----|--------------------------
V(M) potential    | Duda §III    | ✓  | ✓  | ✓  | ✓  | UNRESOLVED (shared M5)
form              |              |    |    |    |    |
f(J·J) form       | Werbos 2017+ | ✓  | ✓  | ~  | ~  | MOSTLY SETTLED — f=gs²;
                  | v5 confirms  |    |    |    |    | factor 2/4 needs check
ω quantization    | Werbos Spec  | ✓  | ✓  | ✓  | ✓  | OPEN — Werbos hypoth.:
mechanism         | §6.1         |    |    |    |    | leptons = 3 lowest
                  |              |    |    |    |    | stable modes. v4 tests.
Q=0 chaoiton      | v2 BVP sweep | —  | ✓  | ~  | ~  | PARTIAL — 23 solutions;
existence         |              |    |    |    |    | ground state pending;
                  |              |    |    |    |    | Q14 conflict to resolve
Lepton mass       | v1 vs v5     | ✓  | ✓  | ✓  | ?  | UNTESTED with Lean
spectrum          |              |    |    |    |    | ODE — Track 1 will
                  |              |    |    |    |    | settle this in v4
Electron          | v3 result    | —  | —  | ✓  | ✓  | OPEN — H/Q=1.723 (1.56%
calibration       |              |    |    |    |    | gap) but not localized.
localization      |              |    |    |    |    | T3 primary fix; T1
                  |              |    |    |    |    | baseline depends on it.
```

Active question count at start of v4: 3 IMMEDIATE (Q12, Q13, Q14) + 4 STILL OPEN (Q2, Q3, Q6, Q8) = **7 load-bearing**.

The single highest-value action: run Track 1 (lepton scan) tomorrow. If our 3 lowest stable ω match leptons within Werbos's 4-6% gap claim, Q2 is effectively resolved empirically — the framework's discrete spectrum claim has direct numerical backing.

---

## ⤷ continued in sandbox_v5

After T9 closed-negative on 2026-05-20 PM (six independent forward-IVP /
single-eigenvalue-BVP attempts could not reach a Q_CS=1 chaoiton), Paul
Werbos responded ~2:00 PM the same day with the actual algorithm: a
collocation BVP with Q_CS=1 as an integral constraint, two free
eigenvalues (ω, λ_LM), Robin BCs at R_max, and a tabulated exponential
initial profile.

This is a structurally different method family from anything attempted
in v4 — it is its own sandbox iteration. See:

- `0c_sandbox_v5.md` — full v5 recipe + attempt log + the partial-
  success / 31×-H/Q-gap diagnosis
- `sandbox_v5/m6_v5_4fn_lambda_bvp.py` — implementation

v5 attempt 4 achieved the first-ever clean `solve_bvp.status = 0`
convergence to a Q_CS=1 chaoiton at the stated (ω≈1, m_eff²≈-0.5)
regime. Open question for v6: Q_CS / H normalization-convention
mismatch giving a 31× H/Q calibration gap.

Also documented in `0c_sandbox_v5.md`: the Duda public-list critique
of v8 LoE (2026-05-20 ~1:15 PM) and Werbos's algorithm reply
(2026-05-20 ~2:00 PM).
