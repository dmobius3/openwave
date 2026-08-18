# M8.4 finding: the kinematic close

> **The task's first result, and it is analytic.** Native, single-valued fields on the closed
> quotient `S³/2I` carry none of the eight nontrivial McKay slots, at any harmonic level, for
> every family and every transport convention. The obstruction sits in the configuration
> space, before any field equation, so OQ1's native branch
> ([`../m8_theory_canonical.md`](../m8_theory_canonical.md)) closes negative by theorem, with
> no dynamics run. The control confirms the full `n = d` first-occurrence ladder in the
> twisted sectors over the quotient, exactly where the M8.2 lock's kinematic tables put it.
> Reproduce: `python3 ../scripts/m8_4_kinematic_check.py`, exit 0, final line
> `SUMMARY: ALL PASS`. Claim ceiling: an analytic result verified by two computationally
> separate in-artifact routes with mutation-armed checks; a statement about kinematic
> accessibility, not about the physics of any borrowed family (§ 6).

## 1. Claim

Let `X = S³/2I` with `2I` the binary icosahedral group (`|2I| = 120`), acting freely by deck
transformations on the cover `S³`, and let T-a denote the frozen target structure of
[`m8_2_preregistration.md`](m8_2_preregistration.md) § 2: the 8 nontrivial `2I` irreps at
McKay distances `d(ρ)`, first occurring in the cover harmonics `V_n` at level `n = d(ρ)`
(`j_first = d/2`).

**Theorem.** Every single-valued field on `X`, valued in any native geometric bundle
(functions, `Γ(TX)`, `Ω¹(X)`, symmetric tensors, or any bundle functorially built from the
tangent bundle), has deck-isotypic content only in the trivial representation: for every
nontrivial irrep `ρ`, the isotypic projector `P_ρ` annihilates its lift, at every harmonic
level, including levels whose ambient cover space contains a `ρ`-summand alongside the
invariant one. Consequently no state of any native field theory posed on `X` occupies any
nontrivial McKay slot.

**Corollary (parity).** The four spinorial slots, `ρ(−1) = −1` at distances `d = 1, 3, 5, 7`,
among them the electron slot `R_1` at `d = 1`, are absent for the simplest reason: everything
on `X` is `(−1)`-even, while those slots live at odd levels only.

**Control.** Sections of the flat twisted bundle `E_ρ` over `X` draw exactly on the
`ρ`-isotype of the cover tower and first occur at `n = d(ρ)`: the T-a ladder is realized,
with eigenvalue ladder `d(ρ)(d(ρ)+2)/R²`, precisely in the twisted sectors. The target
structure is real; it does not live in the native function spaces.

## 2. Proof

**Step 1 (pullback invariance).** A field on `X` is, by definition of the quotient, a
2I-invariant field on `S³` (equivariance under the deck action combined with the bundle
action). Its harmonic expansion therefore contains only invariant vectors, and for any
nontrivial `ρ` the projector `P_ρ = (dim ρ / 120) Σ_g χ_ρ(g) g` annihilates it, because
distinct isotypic summands are orthogonal. A `ρ`-summand coexisting in the same ambient `V_n`
is a neighbor the invariant vector is orthogonal to, not content it acquires.

**Step 2 (parity, making half the absence transparent).** `−1 ∈ 2I` acts on cover level `n`
by `(−1)ⁿ`, and trivially on every native geometric bundle, since frames transform through
`Ad(−1) = 1`. So all `X`-content sits at even `n`. The spinorial irreps occur only at odd
`n`, first at `n = d`. Disjoint supports; no convention can move a parity.

**Step 3 (the invariant support, quantifying the rest).** The invariant scalar levels on `X`
are `n ∈ {0, 12, 20, 24, 30, 32, 36, ...}`, the classical icosahedral invariant structure.
For completeness: the integer irreps first coexist with invariant content in the same ambient
`V_n` at levels `12` (`d = 2`), `12` (`d = 4`), and, for the two `d = 6` slots, `12` (the
`4`) and `20` (the `3'`), but by Step 1 this coexistence confers no slot content on quotient
states; it is recorded as information only.

## 3. Equation-to-code map

All in [`../scripts/m8_4_kinematic_check.py`](../scripts/m8_4_kinematic_check.py):

| Claim | Check | Red condition |
| --- | --- | --- |
| T-a control: first occurrence of every irrep in `V_n\|_2I` equals `d(ρ)` | CHECK A, by two computationally separate routes: icosian character sums, and the McKay recursion `m(n+1) = A·m(n) − m(n−1)` on the affine E8 graph, asserted equal on every multiplicity, all 9 irreps, `n ≤ 60` | any first occurrence off `d`; any route disagreement |
| Parity: invariant support even-only; spinorial content odd-only | CHECK B | any odd-`n` invariant; any even-`n` spinorial content |
| The theorem: `P_ρ` annihilates invariant vectors even where the ambient level contains `ρ` | CHECK C: an explicit invariant vector at `n = 12` (where the trivial and the `3` coexist), projected through `P_3` built from an explicit 13-dimensional representation; measured `3.0e-08` relative, against `2.7e+00` for the non-invariant mutation arm | a nonvanishing invariant projection; a vanishing non-invariant one |
| The checks can fail | CHECK D: a shifted d-map breaks A; admitting odd levels exposes spinorial content; a dropped McKay edge breaks the recursion route | any mutation that fails to redden its check |

Construction sanity inside the script, each asserted: closure and order of the 120 icosians;
irreducibility and mutual orthogonality of all nine characters with `Σ dim² = 120`; the
explicit representation's traces reproducing `χ_{V_n}` on all 120 elements; projector
idempotence; integrality of every multiplicity.

## 4. Verification record

Two computationally separate in-artifact routes certify the multiplicity table: the icosian
character route and the McKay-recursion route (pure integer linear algebra on the affine E8
graph, no icosians, no characters). The routes share the McKay graph as an input, consumed
as data by the recursion and re-derived internally by the character route (each derived
character is asserted irreducible and orthonormal, so a wrong edge fails the construction
rather than propagating): separate computations, not independent derivations. Four mutation
arms confirm the checks discriminate. The
script was independently rerun, and the recursion route independently reproduced, in
author-side review before this note was filed; both agree with the embedded output.

Literature corroboration, registered in [`../../theory/_CITATIONS.md`](../../theory/_CITATIONS.md):
Lachièze-Rey 2004 (arXiv:[math/0401153](https://arxiv.org/abs/math/0401153)) establishes the
selection rule this theorem's Step 1 rests on, that eigenmodes of `S³/Γ` are the Γ-invariant
eigenmodes of `S³`. The slot-inaccessibility reading of that rule is this note's, not the
paper's.

## 5. Consequences for the task

1. **OQ1's native branch is closed.** Dynamics cannot overcome absence from the
   configuration space, so the scaffold pipeline's stages 1 to 5 are not prerequisites for
   this verdict and were not run. They survive as an optional descriptive native-family
   survey (existence, stability, Derrick on the compact arena, the OQ5 sector census),
   decoupled from OQ1.
2. **The slot-bearing continuation is a different object class, owned by this column.**
   Extending a family's differential expression to twisted coefficient bundles changes the
   configuration space, and under
   [`CROSS_MODEL_TESTING.md`](../../../../../dev_docs/CROSS_MODEL_TESTING.md) a prescription
   belongs to whoever defines it. The next object is `M4_int`, the internally-valued
   replacement model this column defined in [`m8_2_preregistration.md`](m8_2_preregistration.md)
   § 6.1 (never attributed to native M4 or its author), across the three frozen flat
   connections, with the trivial connection `σ_0` as a mandatory null control: by this same
   theorem its trivial coefficient bundle must show no nontrivial slots, so an instrument
   that finds them there is broken, and no target-bearing result is read until it is fixed.
3. **The free ladder is calibration, never evidence.** Recovering the frozen per-connection
   first-occurrence tables at zero amplitude only certifies the wiring: any Laplace-type
   kinetic term inherits `d(d+2)/R²` from the arena. Target-bearing observables for the
   continuation must be genuinely nonlinear. The preregistration for that survey is a
   separate, later document; this note pre-freezes none of it.

## 6. Not computed, not established

- No dynamics were run, and nothing here evaluates the physics of M4, M5, or M7 on any
  arena. The theorem locates the slots outside native quotient configuration spaces; it is a
  negative about slot-carrying, not about the families, and it changes no MODELS.md cell.
- Nothing here establishes how `M4_int` behaves nonlinearly, in any sector.
- The optional descriptive native survey (stages 1 to 5) is untouched and unscored.
- The theorem is stated for single-valued sections of bundles functorially built from `TX`
  on the closed quotient. Objects outside that scope, twisted bundles, cover-level fields,
  edge or cone structures, are exactly where it does not apply, and the control shows the
  ladder alive there.
