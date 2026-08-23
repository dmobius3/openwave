# M8.4 pre-registration: the McKay-slot question on the flat bundles `E_ρ`

> **Status: FILED.** Companion working notes are author-side and out of git:
> `m8_4_connection_semantics_note.md` (the structural finding that forced the scope election) and
> `m8_4_gate_findings.md` (the connection-axis gate, the reality types, a text/script divergence,
> and a withdrawn direct-sum proposal recorded with its defects). Every computed statement below
> was checked against the frozen [`scripts/m8_2_first_occurrence.py`](../scripts/m8_2_first_occurrence.py)
> by importing its character machinery, never by reimplementing it.

## 0. Scope, the two objects, and the field type

M8.4's first result, the landed [kinematic close](../findings/m8_4_kinematic_close.md), proved
that single-valued fields on `X = S³/2I` carry no nontrivial McKay slot at any level. The slot
survey therefore moves to a twisted object. This pre-registration fixes **two** objects and gives
them different jobs, rather than replacing one with the other.

**Object A, `M4_int`, preserved as a scope and control object.** `ψ ∈ Ω⁰(X; E_{τ_σ})` with
`τ_σ = Sym²(σ)` over the three frozen flat SU(2) connections, exactly as M8.2 § 6.1 defines it.
Its eight-slot comparison is closed N/A for the reason in § 2. **No `M4_int` dynamics are
executed under this pre-registration.** A later three-sector dynamical study, if it is still
worth doing, gets its own compact appendix.

**Object B, the target-bearing family, named `M4L_Erho`.** For each irrep `ρ` of `2I`, a field
`ψ_ρ ∈ Ω⁰(X; E_ρ)` on the flat bundle `E_ρ`, carrying M4-lineage differential dynamics. **`M4L`
means M4-lineage, not native M4**; the suffix names the defining extension.

**The naming sentence, frozen.** *The target-bearing family is an M8-owned twisted extension
using M4-lineage differential dynamics on the flat bundles `E_ρ`. It is not `M4_int` as defined
in M8.2 § 6.1, and no result obtained on it is attributed to native M4.*

**The field type, frozen here and not left to the pilot.** `E_ρ` is the **complex** flat bundle
associated to the complex irrep `ρ`, carrying its `2I`-invariant Hermitian metric, and
`u = ⟨ψ, ψ⟩` is the Hermitian norm-squared entering the potential. The complexification is part
of the M8-owned extension and is stated as such: native M4's field content was real.

The choice is forced rather than merely convenient. Frobenius-Schur indicators, computed:
`R_1, R_2, R_6, R_8` are **quaternionic** (`ν = −1`) and `R_0, R_3, R_4, R_5, R_7` are real
(`ν = +1`). Four of the eight target sectors therefore admit no real form at all, so an
"underlying real bundle" reading is not a uniform alternative: it would double the fibre in those
four and change the endomorphism algebra from `ℝ` to `ℍ`. One complex bundle per complex irrep
keeps the nine-sector McKay indexing literal.

Native M4 is the author-declared geometric displacement field and fills no cell here. A soldered
`M4 + P` remains a third object, out of scope.

## 1. Claim ceiling

Whatever `M4L_Erho` returns is a statement about an M8-owned extension on `S³/2I`. It is not
evidence for or against native M4, not a reproduction of any published table, and not a mass-
spectrum result. The torsion and mass work is a separate lineage and no result here is quoted
into it.

**The limitation that matters most, stated before any run.** OQ1 asks whether a nonlinear field
equation has solutions "whose energies realize the McKay SLOT STRUCTURE, without per-slot
tuning." `M4L_Erho` puts one field in each slot by hand: the eight target sectors are **installed
by construction**, so this family **cannot** answer whether a dynamics SELECTS the McKay slots.
It answers a strictly narrower question: with the slots installed and one action common to all
eight (§ 9), does the nonlinear evolution preserve the `n = d_rho` level structure, or deform it
sector-dependently? A preserved ladder under a common action is evidence about per-slot tuning,
which is half of OQ1's clause; it is not evidence of selection. No deliverable may state or imply
otherwise. The object that would restore the selection question is held out of this
pre-registration and needs its own; see the companion note.

**A second ceiling, on what counts as interesting.** A small-amplitude continuation that stays
near the free `n = d_ρ` cluster is NOT a striking realization of McKay. Perturbative continuation
normally preserves ancestry over some interval, so that outcome is expected and is reported as
**"nonlinear persistence of the installed free structure"**, never as "dynamical realization."
The outcomes that earn more are a common finite-amplitude branch or stability structure across
the eight under one action, a systematic and nontrivial deformation law, or the separately
pre-registered energy relation M8.2's success ladder requires. Which of those, if any, is
claimed must follow from the § 7 quantities and the § 9 common-action constraint, not from
prose.

## 2. Why `M4_int` cannot carry the slot test, stated as a result

An earlier draft argued this as "three coefficient bundles cannot index eight slots." That
phrasing conflates the coefficient bundle with the slot index and is withdrawn. The frozen 9x3
table does have `rho` as its rows. The precise claim, which is checkable, is about **which part
of that table describes sections on `X`**.

A field `psi in Omega^0(X; E_tau)` lifts to a `2I`-equivariant map into `W_tau`. Its level-`n`
content is `(V_n|_2I (x) W_tau)^{2I}`, so a section exists at level `n` exactly when the TRIVIAL
isotype occurs in `V_n (x) tau`. That condition is the table's **`R_0` row**, read across the
three columns, not any column read down.

Computed, and this is the falsifiable form of the claim:

| `tau_sigma` | first section level, computed from scratch | frozen table's `R_0`-row entry |
| --- | --- | --- |
| trivial | 0 | 0 |
| standard | 2 | 2 |
| Galois | 6 | 6 |

The full section supports are trivial `{0, 12, 20, 24, 30, ...}`, standard
`{2, 10, 12, 14, 18, ...}`, Galois `{6, 10, 14, 16, 18, ...}`.

So the `R_0` row is the section content of `M4_int`, and **the eight nontrivial rows are the
ambient decomposition of `V_n (x) tau_sigma`, a fact about the cover, not about what a field on
`X` carries.** M8.2 § 3 says the same in its own words: the gate "certifies the coefficient
bundle only" and "does not imply an admissible vacuum or a defined fluctuation spectrum exists."
The 9x3 table is a kinematic certification fixture, and the frozen contract already declines to
read it as a fluctuation spectrum.

**Therefore `M4_int`'s eight-slot McKay comparison is not applicable by structure**, and the
reason is that the eight rows are ambient, not that three is smaller than eight. If this reading
is wrong, the table above is where it fails: exhibit a level at which `E_{tau_sigma}` carries a
section outside the `R_0`-row support, and the argument collapses.

## 3. `M4L_Erho`: the slot index is the bundle

A flat bundle on `X` is a representation of `π₁(X) = 2I`, so there are exactly nine. Sections are
`(C^∞(S³) ⊗ W_ρ)^{2I}`, and level `n` carries a section exactly when `ρ` occurs in `V_n|_{2I}`.
Computed:

| bundle | `dim_ℂ W_ρ` | type | `d_ρ` | first level | `dim_ℂ H_{ρ,d_ρ}` |
| --- | --- | --- | --- | --- | --- |
| `E_R0` | 1 | real | 0 | 0 | 1 |
| `E_R1` | 2 | quaternionic | 1 | 1 | 2 |
| `E_R3` | 3 | real | 2 | 2 | 3 |
| `E_R6` | 4 | quaternionic | 3 | 3 | 4 |
| `E_R7` | 5 | real | 4 | 4 | 5 |
| `E_R8` | 6 | quaternionic | 5 | 5 | 6 |
| `E_R4` | 3 | real | 6 | 6 | 7 |
| `E_R5` | 4 | real | 6 | 6 | 7 |
| `E_R2` | 2 | quaternionic | 7 | 7 | 8 |

`n = d_ρ` for all nine, and this column is identical to the frozen table's trivial column, a
consistency check on the whole reading. First occurrence is **simple**: `⟨χ_{V_n}, χ_ρ⟩ = 1` at
`n = d_ρ` in every sector, so the first eigenspace `H_{ρ,d_ρ}` has complex dimension exactly
`d_ρ + 1`. It is degenerate for every nontrivial `ρ`, which § 7 depends on.

**`E_R0` is the mandatory null control**, rank one. The other eight are target-bearing.

**No projector, no dominance fraction, no post hoc assignment.** A field valued in `E_ρ` occupies
slot `ρ` by construction. Slot identity is bundle identity.

## 4. The connection-axis gate: resolved, and there is no 9×3

The ruling required, before any 9×3 execution, that a term in the operator be identified that
depends on the frozen connection `σ` independently of the coefficient bundle. **No such term
exists.** M8.2 § 3's formula admits `σ` only through `κ_σ`, and here the coefficient system is
`ρ` itself. Checked against M8.2 § 6.1's operator: the bundle is fixed by `ρ`; the flat
connection on `E_ρ` *is* `ρ`; the base metric is the round metric; the fibre Hermitian form is
the invariant one on `W_ρ`, unique up to scale because `ρ` is irreducible; and `c1`, `c2`,
`v_mode` are scalars and a discrete label.

**The three-connection axis is not part of `M4L_Erho`.** Nine runs per configuration, not
twenty-seven. The three frozen connections keep their role on `M4_int`, where `τ_σ` genuinely
differs across them.

## 5. Vacuum admissibility: `double_well` is excluded in all eight target sectors

M8.2 § 6.1 requires that a nonzero `v_mode` minimum be holonomy-invariant, established as a
global parallel section or an equivariant texture, not merely picked. Computed:
`dim (W_ρ)^{2I} = ⟨χ_ρ, 1⟩ = 0` for every nontrivial `ρ`, and `1` only for `R_0`.

So in each of the eight target bundles the only covariantly constant section is zero, and a fixed
nonzero minimum is not holonomy-invariant. **`v_mode = double_well` is excluded from the eight
nontrivial sectors** unless a specific equivariant texture is derived and filed as an amendment
before any run. `v_mode = density_mod` stays excluded, as M8.2 § 6.1 already excludes it. `E_R0`
may carry `double_well`; it is the one bundle with an invariant vector, and it is a control.

## 6. Calibration carries zero evidentiary credit

`n = d_ρ` at zero amplitude follows from the character table alone, with no dynamics. **G0,
calibration.** A run that recovers the ladder before nonlinear evolution has learned nothing and
is scored zero. Any claim that the free spectrum "reproduces McKay" is barred from the
deliverables.

## 7. The scored observable: the fate of the free first eigenspace

The question, per sector: does the M4-lineage nonlinear evolution **preserve, deform, split, mix,
or destroy** the `n = d_ρ` structure?

Because `H_{ρ,d_ρ}` has complex dimension `d_ρ + 1 > 1` in every target sector, there is no
distinguished single mode to follow, and after linearization about a nonconstant background the
fluctuation operator need not commute with the free Laplacian, so `n` is no longer an exact
quantum number. **The scored object is therefore the whole free first eigenspace `H_{ρ,d_ρ}`
under continuation into the fluctuation operator about the nonlinear state**, never a single
eigenvector and never a bare eigenvalue compared against a harmonic index.

Three quantities are recorded per sector, and the verdict labels are derived deterministically
from them:

| Quantity | What it is |
| --- | --- |
| **cluster position** | the continued eigenvalue cluster descending from `H_{ρ,d_ρ}`, against the free value `d_ρ(d_ρ+2)/R²` |
| **cluster splitting** | the spread within that cluster, zero in the free limit by the degeneracy above |
| **subspace overlap and leakage** | the principal angles between the continued invariant subspace and the original `H_{ρ,d_ρ}`, and the norm fraction leaving it |

These are basis-independent, which is what makes "split" and "mix" mean anything. **The numerical
estimator, the amplitude ladder, the relaxation and convergence criteria, and the thresholds that
turn the three quantities into labels are fixed by the pilot of § 8 and frozen in the execution
appendix before any target sector runs.** The mathematical observable above is frozen now, so the
pilot cannot decide after seeing a spectrum what counts as survival.

**No pooling, no rescue.** Each nontrivial `E_ρ` gets its own frozen record and its own verdict.
Any cross-slot statement must be declared in the execution appendix before the runs and may not
be composed afterwards from whichever sectors cooperated.

## 8. The pilot spends no target information

The pilot fixes estimators, thresholds, and solver engineering. It may use:

- **`E_R0`**, including nonlinear runs; it is the control and carries no target verdict.
- **Free and manufactured configurations on a nontrivial `E_ρ`**: the free `n = d_ρ` spectrum is
  already zero-credit known structure, so rank, degeneracy, eigencluster tracking, and overlap
  machinery may be exercised there. A **manufactured potential with a known synthetic answer** may
  be used to calibrate the label thresholds.

It may **not** run a nonlinear target configuration on a nontrivial `E_ρ`. Every nontrivial
bundle is target-bearing by § 3, so such a run would read a target outcome while claiming not to.
**If one proves indispensable, that sector is declared SPENT and drops from the eight**, with the
remaining seven reported as such. The preference is not to spend one.

## 9. The experiment matrix, frozen before the pilot

**Both `cubic_nls` and `saturating` are target-bearing configurations**, run and reported
separately across all eight sectors. `linear` is calibration only and answers nothing here, per
§ 6. **Best-of-two reporting is prohibited**: both configurations' verdicts are published per
sector whatever they say.

**Anti-tuning, frozen.** Within a configuration, the action, `c1`, `c2`, the amplitude ladder, the
seed-construction rule, the normalization, and the stopping rule are **common across all eight
target bundles**. No sector-specific retuning is permitted. The pilot chooses those common values
before any target run. This is M8.2's original success-ladder logic: whether one frozen action and
coupling set works across the applicable sectors.

## 10. Gates and mutation requirements

Every gate is mutation-tested under the standing rule: for each line the run prints as a pass, the
mutation that makes it wrong is constructed, the check is confirmed red, and both are recorded. A
coverage check proves exact set equality between the gate registry and the executed records, and
incomplete coverage exits nonzero.

| Gate | What it establishes |
| --- | --- |
| G-CHAR | the character table, the nine `d_ρ`, and the Frobenius-Schur types reproduce the frozen values |
| G-SECT | the section-level and eigenspace-dimension table of § 3 is recomputed, not transcribed |
| G-NULL-a | **free-limit arm, on a nontrivial `E_ρ`**: the scoring pipeline of § 7, run at zero amplitude where the answer is fixed by geometry, must report the cluster at `d_ρ(d_ρ+2)/R²`, zero splitting, and unit overlap. Any other output is an instrument defect and fails the run |
| G-NULL-b | **degeneracy arm, on `E_R0`**: `H_{R_0,0}` is 1-dimensional, so intra-cluster splitting is not a possible outcome there. A pipeline that reports splitting on a one-dimensional eigenspace is broken and fails the run |
| G-NULL-c | **manufactured arm**: a synthetic operator with a known continued cluster, **nonzero** intra-cluster splitting, and **nonzero** leakage out of the original subspace; the pipeline must recover all three within the frozen tolerance. It may not be another exact or no-change case, or it would share G-NULL-a's blind spot |
| G-VAC | the vacuum used in each sector is admissible under § 5 |
| G-CALIB | the free spectrum is recovered and recorded as zero-credit calibration |
| G-CONV | the numerical solution is converged under the appendix's frozen criterion |

**Why the null changed shape.** The kinematic close's `σ₀` control was falsifiable because a
slot-detection instrument could report false nontrivial content in the trivial sector. This
architecture has no cross-slot detector, so that exact test has no object, and a control that
merely checked "the scoring layer declines to score `R_0`" would be an architectural tautology
rather than a null. The three arms above restore falsifiability by aiming at the instrument that
actually exists here, the § 7 continuation pipeline: each arm has a fixed correct answer that the
pipeline can get wrong, and G-NULL-a runs on a nontrivial bundle at zero amplitude, spending no
target information because the free limit is zero-credit under § 6.

## 10b. The task transition, stated for the roadmap

`M4_int` receives a **structural N/A** result for the eight-slot physical-section question, on the
§ 2 calculation. `M4L_Erho` is the newly authorized M8-owned target-bearing extension, with the
explicit ceiling that it tests common nonlinear dynamics across already-defined McKay sectors and
**not** dynamical sector selection.

The roadmap's current sentence, "`M4_int` across the three frozen flat connections" as the slot
survey, predates the § 2 calculation and is amended by this filing. The kinematic close's
nomination of `M4_int` was correct that the survey must move to a twisted object and too quick in
reading the 9x3 certification table as eight physical sectors for an `E_{τ_σ}`-valued field. This
is an M8.4 finding that forces a roadmap amendment, not a preregistration departing from the
roadmap.

## 11. Deliverables

A method note with the per-sector records, the execution appendix as frozen, the gate and mutation
records, and the scripts. The `M4_int` inapplicability of § 2 is a deliverable in its own right
and is reported whether or not `M4L_Erho` returns anything interesting.

## 12. A known discrepancy in the frozen record, declared not inherited

M8.2 § 3's text gives the adjoint/vector block's coefficient system as `κ_σ = τ_σ = Sym²σ` with
`(3R_0 trivial, R3 standard, R4 Galois)`. The frozen `m8_2_first_occurrence.py` implements the
trivial entry as the 1-dimensional `R_0`. The text and the mathematics agree with each other:
`Sym²` of the rank-2 trivial representation is 3-dimensional and wholly trivial.

Verified here: the first-occurrence tables under `R_0` and under `3R_0` are identical across all
nine rows, because first occurrence is multiplicity-blind. **This does not reopen M8.5-A**, whose
certified observable was first occurrence and is unaffected.

M8.4 is the first point at which the multiplicity is physical field content, one scalar against
three. **Declared: for `M4_int` the trivial sector is `3R_0`**, following the frozen definition
literally. `M4L_Erho` is unaffected, its trivial control being `E_R0` of rank one by construction.
The script is multiplicity-blind rather than wrong about anything it certified; flagged so its
comment is not later read as a competing declaration.

## 13. Pins

`m8_2_preregistration.md` §§ 2, 3, 6.1; `m8_2_first_occurrence.py`; `m8_5a_method_note.md`
§§ A.5, A.6, B, D; `m8_4_kinematic_close.md` and `m8_4_kinematic_check.py`. Exact commit and hash
pins are filled at filing.
