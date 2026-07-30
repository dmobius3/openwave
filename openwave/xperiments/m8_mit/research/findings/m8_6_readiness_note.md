# M8.6: MIT-M5 lepton-hierarchy comparison (gated readiness audit)

> This is a readiness audit, not a pre-registration. M8.6's intended test (does MIT's
> McKay-distance rule reproduce M5's measured lepton-hierarchy ratios) cannot presently
> run with integrity. This note records why, and the exact conditions under which it can
> be reopened. Task: [`../tasks/m8_6_task_details.md`](../tasks/m8_6_task_details.md).

## At a glance

| Object | What it is | Independent of masses? | Value-converged? | Physical bridge? | Admissible target? |
| --- | --- | --- | --- | --- | --- |
| `1:5.9:15.1` | `(m/m_e)^(1/3)`, cube roots of lepton mass ratios | No, defined from masses | N/A | N/A (circular) | **No** |
| `A<C<B` energies | Stationary-state energies, 3 seed-permutations | Yes, as numerical outputs; census was lepton-motivated | Ordering only, not the ratios | No physical-parameter bridge exists | **Not yet** (§8) |

## 1. What M8.6 set out to test

[`m8_platform_pointers.md § 5`](../m8_platform_pointers.md) names the target as M5's
"eigenvalue hierarchy `1 : 5.9 : 15.1`, with the mass law `E ∝ Λ³` already fixed (so the
physical ratios `1 : 206.8 : 3477.2` come from Λ ratios)". Before writing a
pre-registration against that target, this audit traced its provenance.

## 2. The named target is not independent M5 output

`1 : 5.9 : 15.1` is not measured; it is defined from the charged-lepton masses. From
M5's own findings note ([`m5_9_lepton_mass_clock_findings.md:122-124`](../../../m5_liquid_crystal/research/findings/m5_9_lepton_mass_clock_findings.md), commit `09da7c50`, 2026-07-02):

> "The mass law is `E ~ Λ^3`... To reproduce `E ~ m` the eigenvalue hierarchy must then
> be `Λ ~ m^(1/3)` (`Λ_mu/e = 5.9, Λ_tau/e = 15.1`)... **Reproducing the masses is
> near-tautological.** ... **This is a consistency check, not a parameter-free
> prediction.** ... the eigenvalue **values** (`1 : 5.9 : 15.1`) **remain Yukawa-like
> input**."

`5.9` is `(m_μ/m_e)^(1/3)` (≈5.913) rounded to two figures; `15.1` is
`(m_τ/m_e)^(1/3)` (≈15.150) rounded the same way. Comparing either against the real
mass ratios is comparing a number to its own definition. This target is **inadmissible**
as an independent M5 input to any cross-model test.

## 3. The genuine measured object is different, and it is real

A separate line of work, [M5.21.2](../../../m5_liquid_crystal/research/findings/m5_21_2_census.md) →
[M5.21.2b](../../../m5_liquid_crystal/research/findings/m5_21_2b_note.md) (commit
`9b8af739`, 2026-07-18), independently measures three stationary states of the same
field functional M5 uses throughout its program, seeded as three axis-permutations of
the fixed spectrum `{1, δ=0.3, 0}` (not three amplitudes of one parameter):

```text
A_i = d_i M / h            (finite-difference channel)
u   = 4 Σ_{i<j} tr(C_ij^T C_ij)      curvature density
V   = per term set (T1 trace-target, this census)
E   = h³ Σ_cells (u + ε·D + V)       total stationary-state energy
```

At N=48: `E_A = 5.2611`, `E_C = 22.059`, `E_B = 84.085`, giving `C/A ≈ 4.19`,
`B/A ≈ 15.98` (recomputed here from the raw energies, not taken from rounded prose).
Nothing here is back-solved from lepton masses: no adjustable parameter was fit to hit
206.8 or 3477.2, and the actual mass comparison is explicitly marked unattempted:
*"the 1:206:3477 / Koide mass ambition stays recorded (needs the realistic-parameter
bridge, M5.21.11)"* ([`m5_particle_hunt.md:65`](../../../m5_liquid_crystal/research/m5_particle_hunt.md)).
The numerical energies were not fitted to the known masses, but the census was
lepton-motivated from the outset: Duda's own prescription that kicked it off was
*"searching for local minima..., hopefully getting candidates for 3 leptons"*
([`m5_21_2_census.md:3`](../../../m5_liquid_crystal/research/findings/m5_21_2_census.md));
its independence is therefore numerical rather than hypothesis-blind.

**This is the only candidate M5 target for M8.6 that is not circular.**

## 4. Why cubing does not apply to A, B, C

The `E ∝ Λ^3` law relates energy to an amplitude parameter `Λ` that is swept directly in
the unrelated M5.9 toy setup. A, B, C differ by *seed permutation* at *fixed* δ=0.3;
there is no Λ being varied here. If A, B, C are ever compared to physical masses, the
comparison is direct (energy to energy), not through a cube. Note also that a uniform
units conversion cannot fix this on its own (§8): `C/A≈4.19` and `B/A≈16.0` are dimensionless
ratios, unchanged by any `E_physical = k·E_lattice` rescaling. Reaching `206.8:3477.2`
requires the underlying ratios themselves to change, which means re-running the census
at genuinely different (more physical) functional parameters, not relabeling units on
the existing toy output.

## 5. Why the numbers cannot be used yet even though they are real

Two independent blockers, both stated by the source itself:

1. **Not value-converged.** *"B (2.61) and C (2.35) exceed [the consistency bar]:
   their cores hold structure at cell scale, so their ENERGIES are lattice-contaminated"*
   and *"the ladder is consistency-converged, **not value-converged**... E still drifts
   +7%/+13% between rungs... quote ORDERINGS and geometry, not absolute E"*
   ([`m5_21_2b_note.md:163,165`](../../../m5_liquid_crystal/research/findings/m5_21_2b_note.md)).
   The refinement ladder shows `E_A` moving `4.920 → 5.261 → 5.921` across three grid
   rungs (~20% drift) at fixed δ. Only the ordering `A<C<B` is asserted as robust; the
   ratios 4.2/16.0 are one snapshot in an actively moving sequence.
2. **No physical parameters, not just no units.** These are bare lattice-energy ratios
   at toy parameters (δ=0.3, not a physical biaxiality value). No fm or MeV anchor
   exists: *"the voxel → fm anchor is Q17, unset"* ([`m5_21_2b_note.md:182`]). This is
   not merely an unlabeled axis: the toy-to-physical calibration explicitly deferred to
   M5.21.11 has to change the underlying *ratios* (4.19, 16.0), which a units conversion
   alone cannot do (§4).

A direct magnitude comparison right now would be dimensionally sound in principle but
physically unjustified in practice: there is nothing on the M5 side yet to compare
MIT's predictions against.

## 6. Why an ordering-only fallback is not a substantive MIT-side test either

M5's `A<C<B` ordering has genuine non-mass rationale: A is the lowest-energy state among
the three measured candidates (lowest = lightest, standard physics logic; the source
does not claim global minimality over the full configuration space), and separately,
M5.21.6 measured that the
μ-candidate C decays toward A while the τ-candidate B drains; an independent
stability/decay argument, not a permutation search against known masses.

But MIT's side of any ordering test is not independent. The identities
`e=(R7,triv), μ=(R8,std), τ=(R4,gal)` were assigned in M8.3 by matching to measured PDG
masses, per mass-spectrum.md's own stated rule: fermions sharing every gate-visible
quantum number are told apart only by which entry the measured mass fits:
*"the gates fix the kind; the mass fixes the generation"*
([`mass-spectrum.md:314`](https://github.com/dmobius3/mode-identity-theory/blob/main/files/spectrum/files/mass-spectrum.md)).
`m_e < m_μ < m_τ` is therefore already built into which slot carries which label. Asking
whether the inherited MIT triple comes out light-middle-heavy in its own predicted
masses is true by construction of the labeling, not a finding. At most, a match would
show that M5's independent stability-based ordering agrees with the ordinary observed
mass ordering; a possible M5-side note, not a cross-model validation of the McKay rule.

## 7. Verdict

M8.6 cannot presently perform its intended lepton-hierarchy comparison, on magnitude or
on ordering. This is not "M8.6 failed"; it is a finding that the planned test is not
yet well-posed, for reasons on the M5 side (no value-converged, physically parameterized
target)
and the MIT side (no ordering test independent of the mass-conditioned generation
labels).

## 8. Reopening conditions

> **Amended 2026-07-29.** Condition 3 originally admitted only a direct run at physical
> parameters. M5.21.11's own scope makes clear the physical regime (δ ~ 1e-10, g ~ 1e10)
> is out of lattice reach by any direct method, so condition 3 now also admits a
> preregistered extrapolation route, under guardrails that preserve the anti-circularity
> requirement it was written to enforce. Full reasoning: [PR #374](https://github.com/openwave-labs/openwave/pull/374).
> The tension was first raised from the M5 side, in the design-question section of
> [`m5_21_11_task_details.md`](../../m5_liquid_crystal/research/tasks/m5_21_11_task_details.md#design-question-resolved-2026-07-29-does-the-extrapolated-law-count-as-a-census-at-the-physical-parameters),
> which also records this amendment's adoption.

M8.6 reopens only once ALL of the following hold. Conditions 1-3 are deliberately
stronger than "attach units to the existing numbers": a uniform rescaling cannot change
a ratio, so closing the roughly 49× (`C/A`) and 218× (`B/A`) gaps between the toy ratios
and the physical ones requires new physics inputs, not new labels (§4).

1. M5.21.11 derives a common toy-to-physical specification at the level of the field
   functional, lattice spacing, couplings, and boundary conditions, using M5-side
   considerations independent of the lepton masses.
2. That specification is frozen before the new A/C/B census ladder is run and before
   any resulting physical-regime ratios are compared with the charged-lepton targets.
3. **Physical-regime target.** The comparison uses either (a) stationary-state energies
   from a new census run at independently fixed physical parameters, directly, or (b),
   where such a run is computationally inaccessible, a physical-regime extrapolation
   from a fresh, preregistered census ladder. Under route (b): the
   asymptotic functional form must be derived from M5-side theory independently of the
   lepton target; the rung set, fitting procedure, holdout tests, branch-tracking
   rules, and uncertainty model must be frozen BEFORE any M8.6 performance evaluation
   against the charged-lepton target ratios is run;
   the same frozen framework applies to all three branches (A, C, B); and the three
   EXISTING toy lattice energies (5.2611, 22.059, 84.085) may not enter the new fit as
   data points, an exponent search, or a post-hoc transformation: that would let an
   `E_physical = f(E_lattice)` curve silently become a three-point mass fit under a
   different name. Failure of the preregistered scaling law or the holdout gates
   leaves M8.6 gated. The output of route (b) is classified as a model-based
   physical-regime PREDICTION, not a directly simulated physical census: the affordable
   rungs are measured censuses, the physical point is an extrapolation from them.
4. Grid refinement (route a) or the frozen uncertainty model (route b) establishes a
   usable UNCERTAINTY on `E_C/E_A` and `E_B/E_A`, not merely a stable ordering. Under
   route (b), the uncertainty model must carry a per-rung DISCRETIZATION term
   (established by grid refinement on at least a subset of rungs) alongside the
   extrapolation error: §5's own finding that `E_A` drifts ~20% across three grid
   resolutions at fixed δ, with B and C less consistency-converged than A, means
   per-branch discretization error does not cancel in the ratio and cannot be absorbed
   into the extrapolation-uncertainty term alone.
5. The assignment `A→e, C→μ, B→τ` stays frozen from the pre-existing stability/decay
   rationale (§6), never re-derived from a mass match.
6. No `1:5.9:15.1` Yukawa-derived figure enters the derivation, calibration, or
   validation at any point.
7. The eventual comparison reports the result under the FIRST frozen specification,
   including a negative result, with no exponent search, permutation, or post-hoc
   rescaling.

**Claim ceiling under route (b).** Even a fully successful extrapolation would
establish that M5's independently derived and simulation-validated asymptotic bridge
yields stationary-state energy ratios consistent with the charged-lepton hierarchy,
under the pre-existing `A→e, C→μ, B→τ` assignment. It would not establish a direct
physical simulation, derive the MIT generation mapping, or make `1:5.9:15.1`
independent evidence.

## Not covered by this note

| Item | Why |
| --- | --- |
| Whether M5.21.11, once it lands, will pass or fail against MIT's predictions | Unknowable before the bridge exists |
| Any numeric verdict on the McKay-distance rule | No admissible comparison currently exists |
| M5's own internal physics (confiner mechanism, decay dynamics, virial balance) | M5-side findings, out of MIT's scope to adjudicate |
