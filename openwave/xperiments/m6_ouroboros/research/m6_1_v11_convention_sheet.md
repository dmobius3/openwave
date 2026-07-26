# The v11 convention sheet: the frozen spec M6 validates against

> **Produced by task [M6.1](tasks/m6_1_task_details.md) (2026-07-20).** Pins the M6 validation spec to **LoE v11, Zenodo [20357670](https://zenodo.org/records/20357670)** (2026-05-23), the record the archive era validated against and the M7 3D program carries, per the roadmap's pre-analysis decision ([`m6_roadmap.md`](m6_roadmap.md)). Every convention a later task needs is classified **FIXED-BY-PAPER** (with the section cite) or **GAP** (the paper does not fix it; the consuming task must DERIVE and pre-register it before running, selection-to-match-targets prohibited per the canonical § 6 rules). Machine checks behind every claim here: [`scripts/m6_1_v11_conventions_check.py`](scripts/m6_1_v11_conventions_check.py) → [`data/m6_1_v11_conventions.json`](data/m6_1_v11_conventions.json); method note: [`m6_1_method_note.md`](findings/m6_1_method_note.md).

## 1. The spec as printed (v11, verbatim anchors)

| Object | As printed | Where |
| --- | --- | --- |
| Lagrangian | `ℒ_JA = −F^{μν}F_{μν} − G^{μν}G_{μν} + J^μ A_μ − g (J^μ J_μ)²` | § 2 |
| Field strengths | `F_μν = ∂_μ A_ν − ∂_ν A_μ`, `G_μν = ∂_μ J_ν − ∂_ν J_μ` | § 2 Table 1 |
| Constraints | `∂^μ A_μ = 0`, `∂^μ J_μ = 0` (Lorenz, both fields) | § 2 Table 1 |
| EL equations | `(2.1) □A_μ = J_μ` · `(2.2) □J_μ = A_μ − 4g (J^ν J_ν) J_μ` | § 2 |
| Wave operator | `□ = ∂²/∂t² − ∇²` | § 2 Table 1 |
| Topological charge | `Q_CS = (1/4π²) ∫ ε^{μνρσ} F_μν G_ρσ d⁴x` | § 4 |
| Electron ansatz | `A_0 = 0, A = r̂ × ∇φ(r) cos(ωt)`; `J_0 = 0, J = r̂ × ∇ψ(r) sin(ωt)` | § 5.1 |
| Parameters | `g = 1.0000`, `λ = 1.0`, `ω = 1.0` (electron), `A₀ = B₀ = 0.1`; neutral sector `B₀ = 0.5`, `η = 0.4251` | § 8, Table 2 |
| Physical scale | `H_code × (ℏc/R_phys) = 0.511 MeV → R_phys = 191 fm` | § 8 Step 3 |
| Constrained minimum | chaoiton = minimum of `H′ = H − λQ` (λ the Lagrange multiplier) | § 5 |

## 2. Consistency findings on the printed spec (script-verified)

| # | Finding | Consequence for consumers |
| --- | --- | --- |
| F1 | **The printed EL pair is not the EL system of the printed Lagrangian.** Symbolic derivation (both metric signatures, residual-verified): the printed `ℒ` yields `□A = ±J/4`, `□J = ±A/4 ∓ g(J·J)J` in Lorenz gauge. The unique coefficient reconstruction reproducing (2.1)/(2.2) with the printed interaction signs is `ℒ_ref = −¼F² − ¼G² + J^μA_μ − g(J^μJ_μ)²` under **mostly-plus signature (−,+,+,+)**; under (+,−,−,−) the coupling and quartic signs must flip too. An overall rescaling of ℒ cannot repair the print: the quartic coefficient −4g in (2.2) pins the relative scale | The DYNAMICS OF RECORD is the EL pair (2.1)/(2.2): it is what every numerical artifact in the record integrates. `ℒ_ref` above is the Lagrangian of record consistent with it |
| F2 | **The printed § 5.1 ansatz is identically zero**: `r̂ × ∇φ(r) ≡ 0` for radial `φ(r)` (the gradient of a radial function is parallel to `r̂`). The printed production ansatz is vacuous as written | The actual field configuration exists only in benchmark code; any consumer must state its ansatz explicitly (§ 4 GAP row) |
| F3 | **v11 never states its metric signature**, and F1 shows the choice is load-bearing (it decides which interaction signs are consistent) | Signature (−,+,+,+) is the inferred convention of record (the only one consistent with the printed interaction signs); pre-register it |

## 3. The target table (the H/Q number, all four values)

`e_nat = √(4πα) = 0.302822` (CODATA α; matches the paper's 0.30282). **Physical target: `H/Q = m_e/e_nat = 0.511/0.30282 = 1.6875`** (v11 § 8 Step 2 states this derivation itself).

| Value | What it is | Gap vs 1.6875 | Gap vs 1.6969 |
| --- | --- | --- | --- |
| 1.6875 | the physical target, m_e/e_nat | 0 | 0.554% |
| 1.6890 | v11's own scan output at g = 1.0 (§ 8: "gap 0.090%") | 0.089% ✅ label consistent | 0.466% |
| 1.6918 | the OpenWave reproduction (§ 5.1: "0.30% match") | 0.254% | 0.301% ✅ label consistent only vs 1.6969 |
| 1.6969 | the May-15 model output at g = 1.0625, relabeled "target" by v11 § 5.1 (provenance: canonical § 4) | 0.557% | 0 |

Script-verified label inconsistencies inside v11 itself:

| Printed claim | Arithmetic | Verdict |
| --- | --- | --- |
| Table 2: "Griesi independent 1.6969 \| target 1.6875 \| gap 0.30%" | \|1.6969 − 1.6875\|/1.6875 = **0.557%** | the printed gap is inconsistent with the printed pair; 0.30% only holds vs the model-internal 1.6969 |
| § 8 Step 3: "Q_predicted = 0.3011 ... predicted (not fitted) ... 0.56%" | 0.511/1.6969 = **0.30114** | `Q_predicted = m_e/(H/Q)`: the 0.56% charge "prediction" is the H/Q-vs-physical-target gap restated, not an independent number |
| § 8 Step 3: "R_phys = 191 fm" | requires `H_code = 0.511 × 191/197.327 = 0.4946`, printed nowhere in v11 (every H value v11 does print gives 652-818 fm) | the scale-setting is not reproducible from printed numbers alone; `ℏc/191 fm = 1.0331 MeV` numerically equals the corpus's `m_J = 1.033 MeV` (appears as early as the 2026-05-15 record as the "chaoiton carrier frequency"; identity ✅, interpretation 🔶 open) |

**Convention of record for M6.2: the target is the PHYSICAL 1.6875.** The model-internal 1.6969 is reported alongside as comparison-only, never as "the target" (canonical § 4, the target-drift row).

## 4. The M6.2 pre-registration checklist: FIXED-BY-PAPER vs GAP

| Convention | Status | Content |
| --- | --- | --- |
| Dynamics | ✅ FIXED | the EL pair (2.1)/(2.2) verbatim (§ 2); `ℒ_ref` of § 2-F1 is its consistent Lagrangian |
| Metric signature | ⚠️ GAP (inferred) | not stated; (−,+,+,+) is the unique choice consistent with the printed signs; pre-register it |
| Quartic form | ✅ FIXED | `f(s) = g s²` with s = J^μJ_μ (v11); NOTE the drift: z20274505 uses `f(s) = ½m_J²s + ¼λs²`; v11 has NO explicit J mass term |
| Parameters | ✅ FIXED | g = 1.0, λ = 1.0, ω = 1.0, A₀ = B₀ = 0.1 (§ 8, Table 2) |
| H functional | ❌ GAP | **never printed in v11, in any form**; must be DERIVED canonically (Legendre) from the pinned dynamics, term-by-term, and pre-registered BEFORE any run |
| Q in "H/Q" | ❌ GAP | Table 1 says "energy per unit J-charge" but no operational definition is printed (Noether Q_J of which symmetry? Q_CS? evaluated on what surface/window?); must be derived from the pinned symmetry and pre-registered. The 4D `Q_CS` integral additionally needs a time-window convention for time-periodic fields (none stated; the author's July 8 paper concedes the quantity is "window-defined") |
| Ansatz | ❌ GAP | the printed form is identically zero (§ 2-F2); the consuming task states its exact functional form (fields, angular structure, radial profiles, phases) before running, with its relation to the archive benchmark ([`archive/sandbox_v8/ouroboros_benchmark.py`](archive/sandbox_v8/ouroboros_benchmark.py)) documented |
| Q_CS prefactor | ✅ FIXED (as printed) | `1/(4π²)` (§ 4); consumers do NOT re-derive it silently; if a derivation disagrees, that is a reportable result |
| Targets | ✅ FIXED | physical 1.6875 (derivation in § 3 above); 1.6969 comparison-only |
| Physical scale | 🔶 PARTIAL | procedure fixed (§ 8 Step 3) but not reproducible from printed numbers (needs H_code = 0.4946, unprinted); re-derive alongside |
| Coupling mass scale μ² | ⚠️ GAP (implicit) | dimensionally, `[J] = 1` (forced by G²) makes the J·A term dim 2 and (2.1) `□A = J` dim-inconsistent; both hold only with an implicit `μ² = 1` in code units (audit finding AF2, [method note § 6](findings/m6_1_method_note.md)); any physical-units mapping must state μ, and v11 § 7.1's "[g] = −4" (the superrenormalizability power-counting input) follows from no consistent assignment of the printed ℒ |
| Boundary conditions / solver | ❌ GAP | nothing printed; pre-register |

## 5. Notation-drift resolution (across the record, resolved to v11)

| Drift item | Records | Resolution of record |
| --- | --- | --- |
| `c` coefficient + dagger on the quartic | z20044392 (`− c G² ... f(J_μJ_μ†)`) → dropped by z20146920 | v11 form: no c, no dagger |
| Kinetic normalization | −F² −G² printed everywhere May 12+; DM v12 drops ¼'s "without comment" (canonical § 2) | the ¼'s were never printed but are REQUIRED by the printed EL pair (§ 2-F1): `ℒ_ref` carries them |
| Quartic form | `g s²` (v11) vs `½m_J²s + ¼λs²` (z20274505) | v11: `g s²`; any m_J-term variant is a different spec, version-pin it |
| λ meaning | Lagrange multiplier of `H′ = H − λQ` = the calibration input "λ = 1.0" (v11) vs quartic coupling λ (z20274505) vs the v4 double-use | v11: λ is the constraint multiplier, value 1.0 |
| Topological charge | 4D `Q = (1/4π²)∫εFG d⁴x` from LoE v8 on | v11 § 4 form as printed (window convention GAP, § 4 above) |
| μ/τ harmonics | ω = 11.0/40.7 (g = 1.0625 era) → 12.82/50.0 (v11, g = 1.0) | v11 values; the ω-selection mechanism itself remains open (canonical OQ3) |

## 6. Cross-links

Task: [`tasks/m6_1_task_details.md`](tasks/m6_1_task_details.md) · script: [`scripts/m6_1_v11_conventions_check.py`](scripts/m6_1_v11_conventions_check.py) · data: [`data/m6_1_v11_conventions.json`](data/m6_1_v11_conventions.json) · method note: [`m6_1_method_note.md`](findings/m6_1_method_note.md) · canonical: [`m6_theory_canonical.md`](m6_theory_canonical.md) (§ 2 version lineage, § 4 provenance ledger, § 6 consumption rules) · roadmap: [`m6_roadmap.md`](m6_roadmap.md) (M6.2 pre-registers against § 4 of this sheet).
